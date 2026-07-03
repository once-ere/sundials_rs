/* -----------------------------------------------------------------
 * Translated from src/cvodes/cvodes_ls.c (CVODES 7.7.0).
 * CVODES' linear solver interface (CVLS).
 *
 * Part I contains routines for using CVSLS on forward problems:
 * attach a SUNLinearSolver, Jacobian and preconditioner plumbing,
 * difference-quotient Jacobian and Jacobian-times-vector
 * approximations, and the generic lsetup/lsolve routines called
 * from the nonlinear solvers.
 *
 * Part II contains wrappers for using CVSLS on adjoint (backward)
 * problems.
 *
 * The Rust phrasing of Part I is inherited from the verified CVODE
 * donor (crates/cvode_rs/src/cvode_ls.rs); the CVODES-specific
 * differences are:
 *  - cvLsSolve receives the C `weight` argument as an ewtS selector
 *    (the four CVODES nonlinear solvers pass cv_ewt for state solves
 *    and cv_ewtS[is] for sensitivity solves);
 *  - the current nonlinear iteration is read from the
 *    cv_nls_curiter mirror (in C it is queried from NLSsim / NLSstg /
 *    NLSstg1 / NLS depending on the sensitivity mode; the active
 *    corrector maintains the mirror);
 *  - CVodeGetLinReturnFlagName knows the adjoint codes.
 * -----------------------------------------------------------------*/
use std::cell::RefCell;

use crate::cvodes_impl::*;
use crate::cvodes_ls_impl::*;
use crate::nvector_serial::*;
use crate::sundials_errors::*;
use crate::sundials_linearsolver::*;
use crate::sundials_math::*;
use crate::sundials_matrix::*;
use crate::sundials_types::*;

/* Private constants */
const MIN_INC_MULT: f64 = 1000.0;
const MAX_DQITERS: i32 = 3; /* max. number of attempts to recover in DQ J*v */
const ZERO: f64 = 0.0;
const PT25: f64 = 0.25;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;

/*================================================================
  PART I - forward problems
  ================================================================*/

/*===============================================================
  Access helper (cvLs_AccessLMem)
  ===============================================================*/
pub(crate) fn cvLs_AccessLMem<'a>(
    cv_mem: &'a mut CVodeMem,
    fname: &str,
) -> Result<&'a mut CVLsMem, i32> {
    match &mut cv_mem.cv_lmem {
        LsModule::Ls(ls) => Ok(ls),
        _ => {
            cvProcessError(None, CVLS_LMEM_NULL, line!(), fname, file!(), MSG_LS_LMEM_NULL);
            Err(CVLS_LMEM_NULL)
        }
    }
}

/*===============================================================
  CVSLS Exported functions -- Required
  ===============================================================*/

/*---------------------------------------------------------------
  CVodeSetLinearSolver specifies the linear solver.
  ---------------------------------------------------------------*/
pub fn CVodeSetLinearSolver(
    cv_mem: &mut CVodeMem,
    LS: LinearSolver,
    A: Option<SUNMatrix>,
) -> i32 {
    /* Retrieve the LS type */
    let ls_type = LS.ls_type();

    /* Set flags based on LS type */
    let iterative = ls_type != SUNLINEARSOLVER_DIRECT;
    let matrixbased =
        ls_type != SUNLINEARSOLVER_ITERATIVE && ls_type != SUNLINEARSOLVER_MATRIX_EMBEDDED;

    /* Ensure that A is None when LS is matrix-embedded */
    if ls_type == SUNLINEARSOLVER_MATRIX_EMBEDDED && A.is_some() {
        cvProcessError(Some(cv_mem), CVLS_ILL_INPUT, line!(), "CVodeSetLinearSolver", file!(),
                       "Incompatible inputs: matrix-embedded LS requires NULL matrix");
        return CVLS_ILL_INPUT;
    }

    /* Check for compatible LS type and matrix */
    if iterative {
        if matrixbased && A.is_none() {
            cvProcessError(Some(cv_mem), CVLS_ILL_INPUT, line!(), "CVodeSetLinearSolver", file!(),
                           "Incompatible inputs: matrix-iterative LS requires non-NULL matrix");
            return CVLS_ILL_INPUT;
        }
    } else if A.is_none() {
        cvProcessError(Some(cv_mem), CVLS_ILL_INPUT, line!(), "CVodeSetLinearSolver", file!(),
                       "Incompatible inputs: direct LS requires non-NULL matrix");
        return CVLS_ILL_INPUT;
    }

    /* free any existing system solver attached to CVode (RAII on overwrite) */

    let has_matrix = A.is_some();

    /* Allocate memory for CVLsMemRec and set defaults */
    let mut cvls_mem = Box::new(CVLsMem {
        iterative,
        matrixbased,

        /* Set defaults for Jacobian-related fields */
        jacDQ: has_matrix,
        jac: None, /* None + jacDQ => internal cvLsDQJac */
        jbad: SUNTRUE,
        dgmax_jbad: CVLS_DGMAX,

        scalesol: SUNFALSE,

        eplifac: CVLS_EPLIN,
        nrmfac: ZERO,

        LS,
        A,
        savedJ: None, /* allocated in cvLsInitialize */
        ytemp: N_VClone(&cv_mem.cv_tempv),
        x: N_VClone(&cv_mem.cv_tempv),

        msbj: CVLS_MSBJ,
        nje: 0,
        nfeDQ: 0,
        nstlj: 0,
        npe: 0,
        nli: 0,
        nps: 0,
        ncfl: 0,
        njtsetup: 0,
        njtimes: 0,
        tnlj: ZERO,

        /* Set defaults for preconditioner-related fields */
        pset: None,
        psolve: None,
        prec_module: PrecModule::None,

        /* Jacobian-times-vector: internal DQ by default */
        jtimesDQ: SUNTRUE,
        jtsetup: None,
        jtimes: None, /* None + jtimesDQ => internal cvLsDQJtimes */
        jt_f: cv_mem.cv_f,

        user_linsys: SUNFALSE,
        linsys: None, /* None => internal cvLsLinSys */

        setup_disabled: SUNFALSE,

        last_flag: CVLS_SUCCESS,
    });

    /* For iterative LS, compute default norm conversion factor */
    if iterative {
        cvls_mem.nrmfac = SUNRsqrt(N_VGetLength(&cvls_mem.ytemp) as f64);
    }

    /* Check if solution scaling should be enabled */
    cvls_mem.scalesol = matrixbased && cv_mem.cv_lmm == CV_BDF;

    /* Attach linear solver memory to integrator memory */
    cv_mem.cv_lmem = LsModule::Ls(cvls_mem);

    CVLS_SUCCESS
}

/*===============================================================
  Optional Set routines
  ===============================================================*/

/* CVodeSetJacFn specifies the Jacobian function. */
pub fn CVodeSetJacFn(cv_mem: &mut CVodeMem, jac: Option<CVLsJacFn>) -> i32 {
    let cvls_mem = match cvLs_AccessLMem(cv_mem, "CVodeSetJacFn") {
        Ok(m) => m,
        Err(e) => return e,
    };

    /* return with failure if jac cannot be used */
    if jac.is_some() && cvls_mem.A.is_none() {
        cvProcessError(None, CVLS_ILL_INPUT, line!(), "CVodeSetJacFn", file!(),
                       "Jacobian routine cannot be supplied for NULL SUNMatrix");
        return CVLS_ILL_INPUT;
    }

    /* set the Jacobian routine pointer, and update relevant flags */
    if let Some(j) = jac {
        cvls_mem.jacDQ = SUNFALSE;
        cvls_mem.jac = Some(j);
    } else {
        cvls_mem.jacDQ = SUNTRUE;
        cvls_mem.jac = None;
    }

    /* ensure the internal linear system function is used */
    cvls_mem.user_linsys = SUNFALSE;
    cvls_mem.linsys = None;

    CVLS_SUCCESS
}

/* CVodeSetDeltaGammaMaxBadJac specifies the maximum gamma ratio change
 * after a NLS convergence failure with a potentially bad Jacobian. If
 * |gamma/gammap-1| < dgmax_jbad then the Jacobian is marked as bad */
pub fn CVodeSetDeltaGammaMaxBadJac(cv_mem: &mut CVodeMem, dgmax_jbad: f64) -> i32 {
    let cvls_mem = match cvLs_AccessLMem(cv_mem, "CVodeSetDeltaGammaMaxBadJac") {
        Ok(m) => m,
        Err(e) => return e,
    };

    /* Set value or use default */
    cvls_mem.dgmax_jbad = if dgmax_jbad <= ZERO { CVLS_DGMAX } else { dgmax_jbad };
    CVLS_SUCCESS
}

/* CVodeSetEpsLin specifies the nonlinear -> linear tolerance scale factor */
pub fn CVodeSetEpsLin(cv_mem: &mut CVodeMem, eplifac: f64) -> i32 {
    let cvls_mem = match cvLs_AccessLMem(cv_mem, "CVodeSetEpsLin") {
        Ok(m) => m,
        Err(e) => return e,
    };

    /* Check for legal eplifac */
    if eplifac < ZERO {
        cvProcessError(None, CVLS_ILL_INPUT, line!(), "CVodeSetEpsLin", file!(),
                       MSG_LS_BAD_EPLIN);
        return CVLS_ILL_INPUT;
    }

    cvls_mem.eplifac = if eplifac == ZERO { CVLS_EPLIN } else { eplifac };
    CVLS_SUCCESS
}

/* CVodeSetLSNormFactor sets or computes the factor to use when converting
   from the integrator tolerance to the linear solver tolerance (WRMS to
   L2 norm). */
pub fn CVodeSetLSNormFactor(cv_mem: &mut CVodeMem, nrmfac: f64) -> i32 {
    let cvls_mem = match cvLs_AccessLMem(cv_mem, "CVodeSetLSNormFactor") {
        Ok(m) => m,
        Err(e) => return e,
    };

    if nrmfac > ZERO {
        /* user-provided factor */
        cvls_mem.nrmfac = nrmfac;
    } else if nrmfac < ZERO {
        /* compute factor for WRMS norm with dot product */
        N_VConst(ONE, &mut cvls_mem.ytemp);
        cvls_mem.nrmfac = SUNRsqrt(N_VDotProd(&cvls_mem.ytemp, &cvls_mem.ytemp));
    } else {
        /* compute default factor for WRMS norm from vector length */
        cvls_mem.nrmfac = SUNRsqrt(N_VGetLength(&cvls_mem.ytemp) as f64);
    }

    CVLS_SUCCESS
}

/* CVodeSetJacEvalFrequency specifies the frequency for recomputing the
   Jacobian matrix and/or preconditioner */
pub fn CVodeSetJacEvalFrequency(cv_mem: &mut CVodeMem, msbj: i64) -> i32 {
    let cvls_mem = match cvLs_AccessLMem(cv_mem, "CVodeSetJacEvalFrequency") {
        Ok(m) => m,
        Err(e) => return e,
    };

    /* Check for legal msbj */
    if msbj < 0 {
        cvProcessError(None, CVLS_ILL_INPUT, line!(), "CVodeSetJacEvalFrequency", file!(),
                       "A negative evaluation frequency was provided.");
        return CVLS_ILL_INPUT;
    }

    cvls_mem.msbj = if msbj == 0 { CVLS_MSBJ } else { msbj };
    CVLS_SUCCESS
}

/* CVodeSetLinearSolutionScaling enables or disables scaling the linear
   solver solution to account for changes in gamma. */
pub fn CVodeSetLinearSolutionScaling(cv_mem: &mut CVodeMem, onoff: bool) -> i32 {
    let lmm = cv_mem.cv_lmm;
    let cvls_mem = match cvLs_AccessLMem(cv_mem, "CVodeSetLinearSolutionScaling") {
        Ok(m) => m,
        Err(e) => return e,
    };

    /* check for valid solver and method type */
    if !cvls_mem.matrixbased || lmm != CV_BDF {
        return CVLS_ILL_INPUT;
    }

    /* set solution scaling flag */
    cvls_mem.scalesol = onoff;
    CVLS_SUCCESS
}

/* CVodeSetPreconditioner specifies the user-supplied preconditioner
   setup and solve routines */
pub fn CVodeSetPreconditioner(
    cv_mem: &mut CVodeMem,
    psetup: Option<CVLsPrecSetupFn>,
    psolve: Option<CVLsPrecSolveFn>,
) -> i32 {
    let cvls_mem = match cvLs_AccessLMem(cv_mem, "CVodeSetPreconditioner") {
        Ok(m) => m,
        Err(e) => return e,
    };

    /* store function pointers for user-supplied routines in CVLs interface */
    cvls_mem.pset = psetup;
    cvls_mem.psolve = psolve;
    cvls_mem.prec_module = PrecModule::User;

    CVLS_SUCCESS
}

/* CVodeSetJacTimes specifies the user-supplied Jacobian-vector product
   setup and multiply routines */
pub fn CVodeSetJacTimes(
    cv_mem: &mut CVodeMem,
    jtsetup: Option<CVLsJacTimesSetupFn>,
    jtimes: Option<CVLsJacTimesVecFn>,
) -> i32 {
    let cv_f = cv_mem.cv_f;
    let cvls_mem = match cvLs_AccessLMem(cv_mem, "CVodeSetJacTimes") {
        Ok(m) => m,
        Err(e) => return e,
    };

    /* store function pointers for user-supplied routines in CVLs
       interface (NULL jtimes implies use of DQ default) */
    if jtimes.is_some() {
        cvls_mem.jtimesDQ = SUNFALSE;
        cvls_mem.jtsetup = jtsetup;
        cvls_mem.jtimes = jtimes;
    } else {
        cvls_mem.jtimesDQ = SUNTRUE;
        cvls_mem.jtsetup = None;
        cvls_mem.jtimes = None;
        cvls_mem.jt_f = cv_f;
    }

    CVLS_SUCCESS
}

/* CVodeSetJacTimesRhsFn specifies an alternative user-supplied ODE
   right-hand side function to use in the internal finite difference
   Jacobian-vector product */
pub fn CVodeSetJacTimesRhsFn(cv_mem: &mut CVodeMem, jtimesRhsFn: Option<CVRhsFn>) -> i32 {
    let cv_f = cv_mem.cv_f;
    let cvls_mem = match cvLs_AccessLMem(cv_mem, "CVodeSetJacTimesRhsFn") {
        Ok(m) => m,
        Err(e) => return e,
    };

    /* check if using internal finite difference approximation */
    if !cvls_mem.jtimesDQ {
        cvProcessError(None, CVLS_ILL_INPUT, line!(), "CVodeSetJacTimesRhsFn", file!(),
                       "Internal finite-difference Jacobian-vector product is disabled.");
        return CVLS_ILL_INPUT;
    }

    /* store function pointers for RHS function (NULL implies use ODE RHS) */
    match jtimesRhsFn {
        Some(f) => cvls_mem.jt_f = Some(f),
        None => cvls_mem.jt_f = cv_f,
    }

    CVLS_SUCCESS
}

/* CVodeSetLinSysFn specifies the linear system setup function. */
pub fn CVodeSetLinSysFn(cv_mem: &mut CVodeMem, linsys: Option<CVLsLinSysFn>) -> i32 {
    let cvls_mem = match cvLs_AccessLMem(cv_mem, "CVodeSetLinSysFn") {
        Ok(m) => m,
        Err(e) => return e,
    };

    /* return with failure if linsys cannot be used */
    if linsys.is_some() && cvls_mem.A.is_none() {
        cvProcessError(None, CVLS_ILL_INPUT, line!(), "CVodeSetLinSysFn", file!(),
                       "Linear system setup routine cannot be supplied for NULL SUNMatrix");
        return CVLS_ILL_INPUT;
    }

    /* set the linear system routine pointer, and update relevant flags */
    if linsys.is_some() {
        cvls_mem.user_linsys = SUNTRUE;
        cvls_mem.linsys = linsys;
    } else {
        cvls_mem.user_linsys = SUNFALSE;
        cvls_mem.linsys = None;
    }

    CVLS_SUCCESS
}

/*===============================================================
  Optional Get routines
  ===============================================================*/

pub fn CVodeGetJac<'a>(cv_mem: &'a mut CVodeMem, j: &mut Option<&'a SUNMatrix>) -> i32 {
    let cvls_mem = match cvLs_AccessLMem(cv_mem, "CVodeGetJac") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *j = cvls_mem.savedJ.as_ref();
    CVLS_SUCCESS
}

pub fn CVodeGetJacTime(cv_mem: &mut CVodeMem, t_j: &mut f64) -> i32 {
    let cvls_mem = match cvLs_AccessLMem(cv_mem, "CVodeGetJacTime") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *t_j = cvls_mem.tnlj;
    CVLS_SUCCESS
}

pub fn CVodeGetJacNumSteps(cv_mem: &mut CVodeMem, nst_j: &mut i64) -> i32 {
    let cvls_mem = match cvLs_AccessLMem(cv_mem, "CVodeGetJacNumSteps") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *nst_j = cvls_mem.nstlj;
    CVLS_SUCCESS
}

/* CVodeGetLinWorkSpace returns the length of workspace allocated for
   the CVLS linear solver interface */
pub fn CVodeGetLinWorkSpace(cv_mem: &mut CVodeMem, lenrw_ls: &mut i64, leniw_ls: &mut i64) -> i32 {
    let (lrw1, liw1) = N_VSpace(&cv_mem.cv_tempv);
    let cvls_mem = match cvLs_AccessLMem(cv_mem, "CVodeGetLinWorkSpace") {
        Ok(m) => m,
        Err(e) => return e,
    };

    /* start with fixed sizes plus vector/matrix pointers */
    *lenrw_ls = 2;
    *leniw_ls = 30;

    /* add NVector sizes */
    *lenrw_ls += 2 * lrw1;
    *leniw_ls += 2 * liw1;

    /* add SUNMatrix size (only account for the one owned by Ls interface) */
    if let Some(saved) = &cvls_mem.savedJ {
        let (lrw, liw) = saved.space();
        *lenrw_ls += lrw;
        *leniw_ls += liw;
    }

    /* add LS sizes */
    {
        let (lrw, liw) = cvls_mem.LS.space();
        *lenrw_ls += lrw;
        *leniw_ls += liw;
    }

    CVLS_SUCCESS
}

/* CVodeGetNumJacEvals returns the number of Jacobian evaluations */
pub fn CVodeGetNumJacEvals(cv_mem: &mut CVodeMem, njevals: &mut i64) -> i32 {
    let cvls_mem = match cvLs_AccessLMem(cv_mem, "CVodeGetNumJacEvals") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *njevals = cvls_mem.nje;
    CVLS_SUCCESS
}

/* CVodeGetNumLinRhsEvals returns the number of calls to the ODE function
   needed for the DQ Jacobian approximation or J*v product approximation */
pub fn CVodeGetNumLinRhsEvals(cv_mem: &mut CVodeMem, nfevals_ls: &mut i64) -> i32 {
    let cvls_mem = match cvLs_AccessLMem(cv_mem, "CVodeGetNumLinRhsEvals") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *nfevals_ls = cvls_mem.nfeDQ;
    CVLS_SUCCESS
}

/* CVodeGetNumPrecEvals returns the number of calls to the user- or
   CVode-supplied preconditioner setup routine */
pub fn CVodeGetNumPrecEvals(cv_mem: &mut CVodeMem, npevals: &mut i64) -> i32 {
    let cvls_mem = match cvLs_AccessLMem(cv_mem, "CVodeGetNumPrecEvals") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *npevals = cvls_mem.npe;
    CVLS_SUCCESS
}

/* CVodeGetNumPrecSolves returns the number of calls to the user- or
   CVode-supplied preconditioner solve routine */
pub fn CVodeGetNumPrecSolves(cv_mem: &mut CVodeMem, npsolves: &mut i64) -> i32 {
    let cvls_mem = match cvLs_AccessLMem(cv_mem, "CVodeGetNumPrecSolves") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *npsolves = cvls_mem.nps;
    CVLS_SUCCESS
}

/* CVodeGetNumLinIters returns the number of linear iterations */
pub fn CVodeGetNumLinIters(cv_mem: &mut CVodeMem, nliters: &mut i64) -> i32 {
    let cvls_mem = match cvLs_AccessLMem(cv_mem, "CVodeGetNumLinIters") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *nliters = cvls_mem.nli;
    CVLS_SUCCESS
}

/* CVodeGetNumLinConvFails returns the number of linear solver
   convergence failures (as reported by the LS object) */
pub fn CVodeGetNumLinConvFails(cv_mem: &mut CVodeMem, nlcfails: &mut i64) -> i32 {
    let cvls_mem = match cvLs_AccessLMem(cv_mem, "CVodeGetNumLinConvFails") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *nlcfails = cvls_mem.ncfl;
    CVLS_SUCCESS
}

/* CVodeGetNumJTSetupEvals returns the number of calls to the
   user-supplied Jacobian-vector product setup routine */
pub fn CVodeGetNumJTSetupEvals(cv_mem: &mut CVodeMem, njtsetups: &mut i64) -> i32 {
    let cvls_mem = match cvLs_AccessLMem(cv_mem, "CVodeGetNumJTSetupEvals") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *njtsetups = cvls_mem.njtsetup;
    CVLS_SUCCESS
}

/* CVodeGetNumJtimesEvals returns the number of calls to the
   Jacobian-vector product multiply routine */
pub fn CVodeGetNumJtimesEvals(cv_mem: &mut CVodeMem, njvevals: &mut i64) -> i32 {
    let cvls_mem = match cvLs_AccessLMem(cv_mem, "CVodeGetNumJtimesEvals") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *njvevals = cvls_mem.njtimes;
    CVLS_SUCCESS
}

/* CVodeGetLinSolveStats returns statistics related to the linear solve. */
pub fn CVodeGetLinSolveStats(
    cv_mem: &mut CVodeMem,
    njevals: &mut i64,
    nfevals_ls: &mut i64,
    nliters: &mut i64,
    nlcfails: &mut i64,
    npevals: &mut i64,
    npsolves: &mut i64,
    njtsetups: &mut i64,
    njtimes: &mut i64,
) -> i32 {
    let cvls_mem = match cvLs_AccessLMem(cv_mem, "CVodeGetLinSolveStats") {
        Ok(m) => m,
        Err(e) => return e,
    };

    *njevals = cvls_mem.nje;
    *nfevals_ls = cvls_mem.nfeDQ;
    *nliters = cvls_mem.nli;
    *nlcfails = cvls_mem.ncfl;
    *npevals = cvls_mem.npe;
    *npsolves = cvls_mem.nps;
    *njtsetups = cvls_mem.njtsetup;
    *njtimes = cvls_mem.njtimes;

    CVLS_SUCCESS
}

/* CVodeGetLastLinFlag returns the last flag set in a CVLS function */
pub fn CVodeGetLastLinFlag(cv_mem: &mut CVodeMem, flag: &mut i64) -> i32 {
    let cvls_mem = match cvLs_AccessLMem(cv_mem, "CVodeGetLastLinFlag") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *flag = cvls_mem.last_flag as i64;
    CVLS_SUCCESS
}

/* CVodeGetLinReturnFlagName translates from the integer error code
   returned by an CVLs routine to the corresponding string equivalent
   for that flag */
pub fn CVodeGetLinReturnFlagName(flag: i64) -> String {
    match flag as i32 {
        CVLS_SUCCESS => "CVLS_SUCCESS",
        CVLS_MEM_NULL => "CVLS_MEM_NULL",
        CVLS_LMEM_NULL => "CVLS_LMEM_NULL",
        CVLS_ILL_INPUT => "CVLS_ILL_INPUT",
        CVLS_MEM_FAIL => "CVLS_MEM_FAIL",
        CVLS_PMEM_NULL => "CVLS_PMEM_NULL",
        CVLS_JACFUNC_UNRECVR => "CVLS_JACFUNC_UNRECVR",
        CVLS_JACFUNC_RECVR => "CVLS_JACFUNC_RECVR",
        CVLS_SUNMAT_FAIL => "CVLS_SUNMAT_FAIL",
        CVLS_SUNLS_FAIL => "CVLS_SUNLS_FAIL",
        CVLS_NO_ADJ => "CVLS_NO_ADJ",
        CVLS_LMEMB_NULL => "CVLS_LMEMB_NULL",
        _ => "NONE",
    }
    .to_string()
}

/*=================================================================
  CVSLS private functions
  =================================================================*/

/*-----------------------------------------------------------------
  cvLsDenseDQJac

  Generates a dense difference quotient approximation to the
  Jacobian of f(t,y).  y is cv_mem.cv_y (perturbed in place and
  restored, exactly as the C code perturbs the caller's vector);
  fy is cv_mem.cv_ftemp; ftemp workspace is cv_mem.cv_vtemp1.
  (CVODES tests constraint presence through the cv_constraints
  vector itself: empty = C NULL.)
  -----------------------------------------------------------------*/
fn cvLsDenseDQJac(
    cv_mem: &mut CVodeMem,
    nfeDQ: &mut i64,
    t: f64,
    jac: &mut crate::sunmatrix_dense::DenseMatrix,
) -> i32 {
    let n = jac.n;

    /* Set minimum increment based on uround and norm of f */
    let srur = SUNRsqrt(cv_mem.cv_uround);
    let fnorm = N_VWrmsNorm(&cv_mem.cv_ftemp, &cv_mem.cv_ewt);
    let min_inc = if fnorm != ZERO {
        MIN_INC_MULT * SUNRabs(cv_mem.cv_h) * cv_mem.cv_uround * n as f64 * fnorm
    } else {
        ONE
    };

    let constraints_set = !cv_mem.cv_constraints.data.is_empty();
    let f = cv_mem.cv_f.unwrap();
    let mut retval = 0;

    for j in 0..(n as usize) {
        /* Generate the jth col of J(tn,y) */
        let yjsaved = cv_mem.cv_y.data[j];
        let mut inc = SUNMAX(srur * SUNRabs(yjsaved), min_inc / cv_mem.cv_ewt.data[j]);

        /* Adjust sign(inc) if y_j has an inequality constraint. */
        if constraints_set {
            let conj = cv_mem.cv_constraints.data[j];
            if SUNRabs(conj) == ONE {
                if (yjsaved + inc) * conj < ZERO {
                    inc = -inc;
                }
            } else if SUNRabs(conj) == TWO && (yjsaved + inc) * conj <= ZERO {
                inc = -inc;
            }
        }

        cv_mem.cv_y.data[j] += inc;

        retval = f(t, &cv_mem.cv_y, &mut cv_mem.cv_vtemp1, &mut cv_mem.cv_user_data);
        *nfeDQ += 1;
        if retval != 0 {
            break;
        }

        cv_mem.cv_y.data[j] = yjsaved;

        let inc_inv = ONE / inc;
        let col_j = jac.col_mut(j as i64);
        for (i, cj) in col_j.iter_mut().enumerate() {
            *cj = inc_inv * (cv_mem.cv_vtemp1.data[i] - cv_mem.cv_ftemp.data[i]);
        }
    }

    retval
}

/*-----------------------------------------------------------------
  cvLsBandDQJac

  Generates a banded difference quotient approximation to the
  Jacobian of f(t,y). ftemp = cv_vtemp1, ytemp = cv_vtemp2.
  -----------------------------------------------------------------*/
fn cvLsBandDQJac(
    cv_mem: &mut CVodeMem,
    nfeDQ: &mut i64,
    t: f64,
    jac: &mut crate::sunmatrix_band::BandMatrix,
) -> i32 {
    let n = jac.n;
    let mupper = jac.mu;
    let mlower = jac.ml;

    /* Load ytemp with y = predicted y vector */
    {
        let CVodeMem { cv_y, cv_vtemp2, .. } = cv_mem;
        cv_vtemp2.data.copy_from_slice(&cv_y.data);
    }

    /* Set minimum increment based on uround and norm of f */
    let srur = SUNRsqrt(cv_mem.cv_uround);
    let fnorm = N_VWrmsNorm(&cv_mem.cv_ftemp, &cv_mem.cv_ewt);
    let min_inc = if fnorm != ZERO {
        MIN_INC_MULT * SUNRabs(cv_mem.cv_h) * cv_mem.cv_uround * n as f64 * fnorm
    } else {
        ONE
    };

    /* Set bandwidth and number of column groups for band differencing */
    let width = mlower + mupper + 1;
    let ngroups = width.min(n);

    let constraints_set = !cv_mem.cv_constraints.data.is_empty();
    let f = cv_mem.cv_f.unwrap();
    let mut retval = 0;

    /* Loop over column groups. */
    for group in 1..=ngroups {
        /* Increment all y_j in group */
        let mut j = group - 1;
        while j < n {
            let ju = j as usize;
            let mut inc = SUNMAX(
                srur * SUNRabs(cv_mem.cv_y.data[ju]),
                min_inc / cv_mem.cv_ewt.data[ju],
            );

            /* Adjust sign(inc) if yj has an inequality constraint. */
            if constraints_set {
                let conj = cv_mem.cv_constraints.data[ju];
                if SUNRabs(conj) == ONE {
                    if (cv_mem.cv_vtemp2.data[ju] + inc) * conj < ZERO {
                        inc = -inc;
                    }
                } else if SUNRabs(conj) == TWO
                    && (cv_mem.cv_vtemp2.data[ju] + inc) * conj <= ZERO
                {
                    inc = -inc;
                }
            }

            cv_mem.cv_vtemp2.data[ju] += inc;
            j += width;
        }

        /* Evaluate f with incremented y */
        retval = {
            let CVodeMem { cv_vtemp2, cv_vtemp1, cv_user_data, .. } = cv_mem;
            f(t, cv_vtemp2, cv_vtemp1, cv_user_data)
        };
        *nfeDQ += 1;
        if retval != 0 {
            break;
        }

        /* Restore ytemp, then form and load difference quotients */
        let mut j = group - 1;
        while j < n {
            let ju = j as usize;
            cv_mem.cv_vtemp2.data[ju] = cv_mem.cv_y.data[ju];
            let mut inc = SUNMAX(
                srur * SUNRabs(cv_mem.cv_y.data[ju]),
                min_inc / cv_mem.cv_ewt.data[ju],
            );

            /* Adjust sign(inc) as before. */
            if constraints_set {
                let conj = cv_mem.cv_constraints.data[ju];
                if SUNRabs(conj) == ONE {
                    if (cv_mem.cv_vtemp2.data[ju] + inc) * conj < ZERO {
                        inc = -inc;
                    }
                } else if SUNRabs(conj) == TWO
                    && (cv_mem.cv_vtemp2.data[ju] + inc) * conj <= ZERO
                {
                    inc = -inc;
                }
            }

            let inc_inv = ONE / inc;
            let i1 = 0.max(j - mupper);
            let i2 = (j + mlower).min(n - 1);
            for i in i1..=i2 {
                let val = inc_inv
                    * (cv_mem.cv_vtemp1.data[i as usize] - cv_mem.cv_ftemp.data[i as usize]);
                jac.set(i, j, val);
            }
            j += width;
        }
    }

    retval
}

/* cvLsDQJac: wrapper choosing the dense or band DQ routine */
fn cvLsDQJac(cv_mem: &mut CVodeMem, nfeDQ: &mut i64, t: f64, jac: &mut SUNMatrix) -> i32 {
    match jac {
        SUNMatrix::Dense(dm) => cvLsDenseDQJac(cv_mem, nfeDQ, t, dm),
        SUNMatrix::Band(bm) => cvLsBandDQJac(cv_mem, nfeDQ, t, bm),
        _ => {
            cvProcessError(Some(cv_mem), CVLS_ILL_INPUT, line!(), "cvLsDQJac", file!(),
                           "unrecognized matrix type for cvLsDQJac");
            CVLS_ILL_INPUT
        }
    }
}

/*-----------------------------------------------------------------
  cvLsLinSys

  Setup the linear system A = I - gamma J
  -----------------------------------------------------------------*/
fn cvLsLinSys(
    cv_mem: &mut CVodeMem,
    cvls_mem: &mut CVLsMem,
    t: f64,
    jok: bool,
    jcur: &mut bool,
    gamma: f64,
) -> i32 {
    /* Check if Jacobian needs to be updated */
    if jok {
        /* Use saved copy of J */
        *jcur = SUNFALSE;

        /* Overwrite linear system matrix with saved J */
        let (saved, a) = (cvls_mem.savedJ.as_ref().unwrap(), cvls_mem.A.as_mut().unwrap());
        let retval = saved.copy_to(a);
        if retval != 0 {
            cvProcessError(Some(cv_mem), CVLS_SUNMAT_FAIL, line!(), "cvLsLinSys", file!(),
                           MSG_LS_SUNMAT_FAILED);
            cvls_mem.last_flag = CVLS_SUNMAT_FAIL;
            return cvls_mem.last_flag;
        }
    } else {
        /* Call jac() routine to update J */
        *jcur = SUNTRUE;

        /* Clear the linear system matrix if necessary */
        if cvls_mem.LS.ls_type() == SUNLINEARSOLVER_DIRECT {
            let retval = cvls_mem.A.as_mut().unwrap().zero();
            if retval != 0 {
                cvProcessError(Some(cv_mem), CVLS_SUNMAT_FAIL, line!(), "cvLsLinSys", file!(),
                               MSG_LS_SUNMAT_FAILED);
                cvls_mem.last_flag = CVLS_SUNMAT_FAIL;
                return cvls_mem.last_flag;
            }
        }

        /* Compute new Jacobian matrix */
        let retval = if cvls_mem.jacDQ {
            let a = cvls_mem.A.as_mut().unwrap();
            cvLsDQJac(cv_mem, &mut cvls_mem.nfeDQ, t, a)
        } else {
            let jac = cvls_mem.jac.unwrap();
            let a = cvls_mem.A.as_mut().unwrap();
            let CVodeMem {
                cv_y,
                cv_ftemp,
                cv_user_data,
                cv_vtemp1,
                cv_vtemp2,
                cv_vtemp3,
                ..
            } = cv_mem;
            jac(t, cv_y, cv_ftemp, a, cv_user_data, cv_vtemp1, cv_vtemp2, cv_vtemp3)
        };
        if retval < 0 {
            cvProcessError(Some(cv_mem), CVLS_JACFUNC_UNRECVR, line!(), "cvLsLinSys", file!(),
                           MSG_LS_JACFUNC_FAILED);
            cvls_mem.last_flag = CVLS_JACFUNC_UNRECVR;
            return -1;
        }
        if retval > 0 {
            cvls_mem.last_flag = CVLS_JACFUNC_RECVR;
            return 1;
        }

        /* Update saved copy of the Jacobian matrix */
        let (a, saved) = (cvls_mem.A.as_ref().unwrap(), cvls_mem.savedJ.as_mut().unwrap());
        let retval = a.copy_to(saved);
        if retval != 0 {
            cvProcessError(Some(cv_mem), CVLS_SUNMAT_FAIL, line!(), "cvLsLinSys", file!(),
                           MSG_LS_SUNMAT_FAILED);
            cvls_mem.last_flag = CVLS_SUNMAT_FAIL;
            return cvls_mem.last_flag;
        }
    }

    /* Perform linear combination A = I - gamma*J */
    let retval = cvls_mem.A.as_mut().unwrap().scale_addi(-gamma);
    if retval != 0 {
        cvProcessError(Some(cv_mem), CVLS_SUNMAT_FAIL, line!(), "cvLsLinSys", file!(),
                       MSG_LS_SUNMAT_FAILED);
        cvls_mem.last_flag = CVLS_SUNMAT_FAIL;
        return cvls_mem.last_flag;
    }

    CVLS_SUCCESS
}

/*-----------------------------------------------------------------
  cvLsInitialize

  This routine performs remaining initializations specific to the
  linear solver interface (and solver itself)
  -----------------------------------------------------------------*/
pub fn cvLsInitialize(cv_mem: &mut CVodeMem, cvls_mem: &mut CVLsMem) -> i32 {
    /* Test for valid combinations of matrix & Jacobian routines */
    if let Some(a) = &cvls_mem.A {
        /* Matrix-based case */
        if !cvls_mem.user_linsys {
            /* Internal linear system function */
            cvls_mem.linsys = None;

            if cvls_mem.jacDQ {
                /* Internal difference quotient Jacobian. Check that A is
                   dense or band, otherwise return an error */
                let ok = matches!(a, SUNMatrix::Dense(_) | SUNMatrix::Band(_));
                if !ok {
                    cvProcessError(Some(cv_mem), CVLS_ILL_INPUT, line!(), "cvLsInitialize",
                                   file!(),
                                   "No Jacobian constructor available for SUNMatrix type");
                    cvls_mem.last_flag = CVLS_ILL_INPUT;
                    return CVLS_ILL_INPUT;
                }
                cvls_mem.jac = None; /* internal DQ */
            }

            /* Allocate internally saved Jacobian if not already done */
            if cvls_mem.savedJ.is_none() {
                cvls_mem.savedJ = Some(a.clone_empty());
            }
        } /* end matrix-based case */
    } else {
        /* Matrix-free case: ensure jac and linsys are disabled */
        cvls_mem.jacDQ = SUNFALSE;
        cvls_mem.jac = None;

        cvls_mem.user_linsys = SUNFALSE;
        cvls_mem.linsys = None;
    }

    /* reset counters */
    cvLsInitializeCounters(cvls_mem);

    /* Set Jacobian-vector product related fields, based on jtimesDQ */
    if cvls_mem.jtimesDQ {
        cvls_mem.jtsetup = None;
        cvls_mem.jtimes = None; /* internal DQ */
    }

    /* if A is NULL and psetup is not present, then cvLsSetup does not
       need to be called, so disable the lsetup dispatch.
       (When the BandPre/BBDPre PrecModule variants land, this condition
       gains the donor's prec_module clause.) */
    cvls_mem.setup_disabled = cvls_mem.A.is_none() && cvls_mem.pset.is_none();

    /* When using a matrix-embedded linear solver, disable lsetup call
       and solution scaling */
    if cvls_mem.LS.ls_type() == SUNLINEARSOLVER_MATRIX_EMBEDDED {
        cvls_mem.setup_disabled = SUNTRUE;
        cvls_mem.scalesol = SUNFALSE;
    }

    /* Call LS initialize routine, and return result */
    cvls_mem.last_flag = cvls_mem.LS.initialize();
    cvls_mem.last_flag
}

pub fn cvLsReInitialize(_cv_mem: &mut CVodeMem, cvls_mem: &mut CVLsMem) -> i32 {
    /* Initialize counters */
    cvLsInitializeCounters(cvls_mem);
    CVLS_SUCCESS
}

/*-----------------------------------------------------------------
  cvLsSetup

  Conditionally calls the LS 'setup' routine, creating the system
  matrix A = I - gamma*J when a SUNMatrix is used, or calling the
  preconditioner setup routine for iterative solvers.
  ypred = cv_mem.cv_y, fpred = cv_mem.cv_ftemp (as passed by the
  cvNlsLSetup / cvNlsLSetupSens* drivers in the C code);
  vtemp1/2/3 = cv_vtemp1/2/3.
  -----------------------------------------------------------------*/
pub fn cvLsSetup(
    cv_mem: &mut CVodeMem,
    cvls_mem: &mut CVLsMem,
    convfail: i32,
    jcur_ptr: &mut bool,
) -> i32 {
    /* Immediately return when using matrix-embedded linear solver */
    if cvls_mem.LS.ls_type() == SUNLINEARSOLVER_MATRIX_EMBEDDED {
        cvls_mem.last_flag = CVLS_SUCCESS;
        return cvls_mem.last_flag;
    }

    /* Use nst, gamma/gammap, and convfail to set J/P eval. flag jok */
    let dgamma = SUNRabs((cv_mem.cv_gamma / cv_mem.cv_gammap) - ONE);
    cvls_mem.jbad = cv_mem.cv_nst == 0
        || cv_mem.first_step_after_resize
        || cv_mem.cv_nst >= cvls_mem.nstlj + cvls_mem.msbj
        || (convfail == CV_FAIL_BAD_J && dgamma < cvls_mem.dgmax_jbad)
        || convfail == CV_FAIL_OTHER;

    /* Setup the linear system if necessary */
    if cvls_mem.A.is_some() {
        /* Update J if appropriate and evaluate A = I - gamma J */
        let jok = !cvls_mem.jbad;
        let retval = if cvls_mem.user_linsys {
            let linsys = cvls_mem.linsys.unwrap();
            let a = cvls_mem.A.as_mut().unwrap();
            let CVodeMem {
                cv_y,
                cv_ftemp,
                cv_user_data,
                cv_vtemp1,
                cv_vtemp2,
                cv_vtemp3,
                cv_tn,
                cv_gamma,
                ..
            } = cv_mem;
            linsys(
                *cv_tn, cv_y, cv_ftemp, a, jok, jcur_ptr, *cv_gamma, cv_user_data, cv_vtemp1,
                cv_vtemp2, cv_vtemp3,
            )
        } else {
            let t = cv_mem.cv_tn;
            let gamma = cv_mem.cv_gamma;
            cvLsLinSys(cv_mem, cvls_mem, t, jok, jcur_ptr, gamma)
        };

        /* Update J eval count and step when J was last updated */
        if *jcur_ptr {
            cvls_mem.nje += 1;
            cvls_mem.nstlj = cv_mem.cv_nst;
            cvls_mem.tnlj = cv_mem.cv_tn;
        }

        /* Check linsys() return value and return if necessary */
        if retval != CVLS_SUCCESS {
            if cvls_mem.user_linsys {
                if retval < 0 {
                    cvProcessError(Some(cv_mem), CVLS_JACFUNC_UNRECVR, line!(), "cvLsSetup",
                                   file!(), MSG_LS_JACFUNC_FAILED);
                    cvls_mem.last_flag = CVLS_JACFUNC_UNRECVR;
                    return -1;
                } else {
                    cvls_mem.last_flag = CVLS_JACFUNC_RECVR;
                    return 1;
                }
            } else {
                return retval;
            }
        }
    } else {
        /* Matrix-free case, set jcur to jbad */
        *jcur_ptr = cvls_mem.jbad;
    }

    /* Call LS setup routine -- for direct solvers this factors A; for
       iterative solvers the generic C SUNLinSolSetup calls cvLsPSetup,
       which passes the heuristic suggestions above to the user code */
    cvls_mem.last_flag = match &mut cvls_mem.LS {
        LinearSolver::Dense(_) | LinearSolver::Band(_) => {
            let CVLsMem { LS, A, .. } = cvls_mem;
            LS.setup(A.as_mut())
        }
        _ => {
            /* iterative solver: preconditioner setup (cvLsPSetup). In C the
               jcurPtr seen here and the &cv_mem->cv_jcur handed to the
               user/module psetup are the same storage; the psetup routines
               here write cv_mem.cv_jcur, so sync it with the caller's
               jcur_ptr across the call. */
            cv_mem.cv_jcur = *jcur_ptr;
            let retval = cvLsPSetup(cv_mem, cvls_mem);
            *jcur_ptr = cv_mem.cv_jcur;
            retval
        }
    };

    /* If Matrix-free, update heuristics flags */
    if cvls_mem.A.is_none() {
        /* If user set jcur to SUNTRUE, increment npe and save nst value */
        if *jcur_ptr {
            cvls_mem.npe += 1;
            cvls_mem.nstlj = cv_mem.cv_nst;
            cvls_mem.tnlj = cv_mem.cv_tn;
        }

        /* Update jcur flag if we suggested an update */
        if cvls_mem.jbad {
            *jcur_ptr = SUNTRUE;
        }
    }

    cvls_mem.last_flag
}

/* cvLsPSetup: interface to the (user or module) preconditioner setup.
   Returns 0 when no preconditioner setup exists (the C SUNLinSolSetup
   is a no-op in that case). */
fn cvLsPSetup(cv_mem: &mut CVodeMem, cvls_mem: &mut CVLsMem) -> i32 {
    let jok = !cvls_mem.jbad;
    match &mut cvls_mem.prec_module {
        PrecModule::None => 0,
        PrecModule::User => {
            if let Some(pset) = cvls_mem.pset {
                let CVodeMem {
                    cv_y,
                    cv_ftemp,
                    cv_user_data,
                    cv_tn,
                    cv_gamma,
                    cv_jcur,
                    ..
                } = cv_mem;
                pset(*cv_tn, cv_y, cv_ftemp, jok, cv_jcur, *cv_gamma, cv_user_data)
            } else {
                0
            }
        }
    }
}

/*-----------------------------------------------------------------
  cvLsSolve

  Interfaces between CVode and the generic linear solver object,
  setting the tolerance and scaling vectors, calling the solver and
  accumulating statistics.  b holds the RHS on entry and the
  solution on return; ynow = cv_y, fnow = cv_ftemp.

  The C `weight` argument becomes the ewtS_is selector: the four
  CVODES nonlinear solvers pass cv_ewt for state solves (None here)
  and cv_ewtS[is] for staggered/staggered1/simultaneous sensitivity
  solves (Some(is) here).
  -----------------------------------------------------------------*/
pub fn cvLsSolve(
    cv_mem: &mut CVodeMem,
    cvls_mem: &mut CVLsMem,
    b: &mut NVector,
    ewtS_is: Option<usize>,
) -> i32 {
    /* get current nonlinear solver iteration. In C this queries
       NLSsim / NLSstg / NLSstg1 / NLS depending on the sensitivity
       mode; the active corrector maintains the cv_nls_curiter mirror
       (donor adaptation), so the four branches collapse. */
    let curiter = cv_mem.cv_nls_curiter;

    /* If the linear solver is iterative: test norm(b), if small, return
       x = 0 or x = b; set linear solver tolerance (in left/right scaled
       2-norm) */
    let delta;
    if cvls_mem.iterative {
        let deltar = cvls_mem.eplifac * cv_mem.cv_tq[4];
        let bnorm = match ewtS_is {
            Some(is) => N_VWrmsNorm(b, &cv_mem.cv_ewtS[is]),
            None => N_VWrmsNorm(b, &cv_mem.cv_ewt),
        };

        if bnorm <= deltar {
            if curiter > 0 {
                N_VConst(ZERO, b);
            }
            cvls_mem.last_flag = CVLS_SUCCESS;
            return cvls_mem.last_flag;
        }
        /* Adjust tolerance for 2-norm */
        delta = deltar * cvls_mem.nrmfac;
    } else {
        delta = ZERO;
    }

    /* Our iterative solvers all accept scaling vectors (s1 = s2 =
       weight), matching the C solvers that implement
       SetScalingVectors; the w_mean tolerance adjustment branch is
       therefore never taken. */

    /* Set initial guess x = 0 to LS */
    N_VConst(ZERO, &mut cvls_mem.x);
    cvls_mem.LS.set_zero_guess(SUNTRUE);

    /* If a user-provided jtsetup routine is supplied, call that here */
    if let Some(jts) = cvls_mem.jtsetup {
        let flag = {
            let CVodeMem { cv_y, cv_ftemp, cv_user_data, cv_tn, .. } = cv_mem;
            jts(*cv_tn, cv_y, cv_ftemp, cv_user_data)
        };
        cvls_mem.last_flag = flag;
        cvls_mem.njtsetup += 1;
        if flag != 0 {
            cvProcessError(Some(cv_mem), flag, line!(), "cvLsSolve", file!(),
                           MSG_LS_JTSETUP_FAILED);
            return cvls_mem.last_flag;
        }
    }

    /* Call solver, and copy x to b (in finish_solve) */
    let retval = if matches!(cvls_mem.LS, LinearSolver::Custom(_)) {
        /* matrix-embedded solver: gets (t, gamma, user_data) */
        let CVLsMem { LS, x, .. } = cvls_mem;
        if let LinearSolver::Custom(cls) = LS {
            let CVodeMem { cv_user_data, cv_tn, cv_gamma, .. } = cv_mem;
            cls.solve(x, b, delta, *cv_tn, *cv_gamma, cv_user_data)
        } else {
            unreachable!()
        }
    } else if !cvls_mem.iterative {
        let CVLsMem { LS, A, x, .. } = cvls_mem;
        match (LS, A.as_mut()) {
            (LinearSolver::Dense(dls), Some(SUNMatrix::Dense(am))) => dls.solve(am, x, b),
            (LinearSolver::Band(bls), Some(SUNMatrix::Band(am))) => bls.solve(am, x, b),
            _ => SUN_ERR_ARG_INCOMPATIBLE,
        }
    } else {
        cvLsSolveIterative(cv_mem, cvls_mem, b, delta, ewtS_is)
    };

    finish_solve(cv_mem, cvls_mem, b, retval, curiter)
}

/* Iterative Krylov solve: builds the ATimes (cvLsATimes) and PSolve
   (cvLsPSolve) callbacks over the integrator memory.  The two closures
   share cv_mem through a RefCell — the solvers never call them
   re-entrantly (checked at runtime).  The error weight vectors are
   detached from CVodeMem for the duration of the solve: cv_ewt always
   (the DQ Jacobian-vector product uses it, exactly as the C
   cvLsDQJtimes reads cv_mem->cv_ewt), plus cv_ewtS[is] when the solve
   weight is a sensitivity error weight; the solve weight is passed as
   the (read-only) scaling vectors s1 = s2. */
fn cvLsSolveIterative(
    cv_mem: &mut CVodeMem,
    cvls_mem: &mut CVLsMem,
    b: &NVector,
    delta: f64,
    ewtS_is: Option<usize>,
) -> i32 {
    let CVLsMem {
        LS,
        x,
        ytemp,
        njtimes,
        nfeDQ,
        nps,
        jtimes,
        jtimesDQ,
        jt_f,
        psolve,
        prec_module,
        ..
    } = cvls_mem;
    let jtimes = *jtimes;
    let jtimes_dq = *jtimesDQ;
    let jt_f = *jt_f;
    let psolve_fn = *psolve;

    let has_psolve = psolve_fn.is_some();

    let ewt = std::mem::take(&mut cv_mem.cv_ewt);
    let weightS = ewtS_is.map(|is| std::mem::take(&mut cv_mem.cv_ewtS[is]));
    let cvm = RefCell::new(&mut *cv_mem);
    let ewt_ref = &ewt;
    let w_ref: &NVector = match &weightS {
        Some(w) => w,
        None => &ewt,
    };

    let mut atimes = |v: &NVector, z: &mut NVector| -> i32 {
        let mut guard = cvm.borrow_mut();
        let cvr: &mut CVodeMem = &mut *guard;
        /* cvLsATimes: z = v - gamma * J v */
        let jret = if jtimes_dq {
            cvLsDQJtimes(cvr, ewt_ref, nfeDQ, jt_f, v, z, ytemp)
        } else {
            let jt = jtimes.unwrap();
            let CVodeMem { cv_y, cv_ftemp, cv_user_data, cv_tn, .. } = cvr;
            jt(v, z, *cv_tn, cv_y, cv_ftemp, cv_user_data, ytemp)
        };
        *njtimes += 1;
        if jret != 0 {
            return jret;
        }
        /* add contribution from identity matrix: z = v - gamma*z */
        z.linear_sum_with(-cvr.cv_gamma, ONE, v);
        0
    };

    let mut psolve_cb = |r: &NVector, z: &mut NVector, tol: f64, lr: i32| -> i32 {
        let mut guard = cvm.borrow_mut();
        let cvr: &mut CVodeMem = &mut *guard;
        let ret = match prec_module {
            PrecModule::User => {
                if let Some(ps) = psolve_fn {
                    let CVodeMem { cv_y, cv_ftemp, cv_user_data, cv_tn, cv_gamma, .. } = cvr;
                    ps(*cv_tn, cv_y, cv_ftemp, r, z, *cv_gamma, tol, lr, cv_user_data)
                } else {
                    0
                }
            }
            PrecModule::None => 0,
        };
        *nps += 1;
        ret
    };

    let retval = if has_psolve {
        LS.solve(
            None,
            x,
            b,
            delta,
            &mut atimes,
            Some(&mut psolve_cb),
            Some(w_ref),
            Some(w_ref),
        )
    } else {
        LS.solve(None, x, b, delta, &mut atimes, None, Some(w_ref), Some(w_ref))
    };

    drop(atimes);
    drop(psolve_cb);
    let cv_mem_back = cvm.into_inner();
    cv_mem_back.cv_ewt = ewt;
    if let Some(is) = ewtS_is {
        cv_mem_back.cv_ewtS[is] = weightS.unwrap();
    }

    retval
}

/* NOTE: the iterative branch above cannot borrow the weight vectors
   while cv_mem is inside the RefCell; cvLsSolve delegates to this
   function instead. */
fn finish_solve(
    cv_mem: &mut CVodeMem,
    cvls_mem: &mut CVLsMem,
    b: &mut NVector,
    retval: i32,
    curiter: i32,
) -> i32 {
    b.data.copy_from_slice(&cvls_mem.x.data);

    /* If using a direct or matrix-iterative solver, BDF method, and gamma
       has changed, scale the correction to account for change in gamma */
    if cvls_mem.scalesol && cv_mem.cv_gamrat != ONE {
        b.scale_inplace(TWO / (ONE + cv_mem.cv_gamrat));
    }

    /* Retrieve statistics from iterative linear solvers */
    let mut nli_inc = 0;
    if cvls_mem.iterative {
        nli_inc = cvls_mem.LS.num_iters();
    }

    /* Increment counters nli and ncfl */
    cvls_mem.nli += nli_inc as i64;
    if retval != SUN_SUCCESS {
        cvls_mem.ncfl += 1;
    }

    /* Interpret solver return value */
    cvls_mem.last_flag = retval;

    match retval {
        SUN_SUCCESS => 0,
        SUNLS_RES_REDUCED => {
            /* allow reduction but not solution on first Newton iteration,
               otherwise return with a recoverable failure */
            if curiter == 0 { 0 } else { 1 }
        }
        SUNLS_CONV_FAIL | SUNLS_ATIMES_FAIL_REC | SUNLS_PSOLVE_FAIL_REC
        | SUNLS_PACKAGE_FAIL_REC | SUNLS_QRFACT_FAIL | SUNLS_LUFACT_FAIL => 1,
        SUN_ERR_ARG_CORRUPT | SUN_ERR_ARG_INCOMPATIBLE | SUN_ERR_MEM_FAIL | SUNLS_GS_FAIL
        | SUNLS_QRSOL_FAIL => -1,
        SUN_ERR_EXT_FAIL => {
            cvProcessError(Some(cv_mem), SUN_ERR_EXT_FAIL, line!(), "cvLsSolve", file!(),
                           "Failure in SUNLinSol external package");
            -1
        }
        SUNLS_ATIMES_FAIL_UNREC => {
            cvProcessError(Some(cv_mem), SUNLS_ATIMES_FAIL_UNREC, line!(), "cvLsSolve", file!(),
                           MSG_LS_JTIMES_FAILED);
            -1
        }
        SUNLS_PSOLVE_FAIL_UNREC => {
            cvProcessError(Some(cv_mem), SUNLS_PSOLVE_FAIL_UNREC, line!(), "cvLsSolve", file!(),
                           MSG_LS_PSOLVE_FAILED);
            -1
        }
        _ => 0,
    }
}

/*-----------------------------------------------------------------
  cvLsDQJtimes

  Difference quotient approximation to the Jacobian-vector product:
  Jv = [f(y + v*sig) - f(y)]/sig with sig = 1/||v||_WRMS.
  y = cv_mem.cv_y, fy = cv_mem.cv_ftemp, work = cvls ytemp.
  -----------------------------------------------------------------*/
fn cvLsDQJtimes(
    cv_mem: &mut CVodeMem,
    ewt: &NVector,
    nfeDQ: &mut i64,
    jt_f: Option<CVRhsFn>,
    v: &NVector,
    jv: &mut NVector,
    work: &mut NVector,
) -> i32 {
    /* Initialize perturbation to 1/||v|| (ewt is cv_ewt, detached from
       cv_mem during the iterative solve) */
    let mut sig = ONE / N_VWrmsNorm(v, ewt);

    let f = jt_f.unwrap();
    let mut retval = 0;
    for _iter in 0..MAX_DQITERS {
        /* Set work = y + sig*v */
        N_VLinearSum(sig, v, ONE, &cv_mem.cv_y, work);

        /* Set Jv = f(tn, y+sig*v) */
        retval = f(cv_mem.cv_tn, work, jv, &mut cv_mem.cv_user_data);
        *nfeDQ += 1;
        if retval == 0 {
            break;
        }
        if retval < 0 {
            return -1;
        }

        /* If f failed recoverably, shrink sig and retry */
        sig *= PT25;
    }

    /* If retval still isn't 0, return with a recoverable failure */
    if retval > 0 {
        return 1;
    }

    /* Replace Jv by (Jv - fy)/sig */
    let siginv = ONE / sig;
    jv.linear_sum_with(siginv, -siginv, &cv_mem.cv_ftemp);

    0
}

/*-----------------------------------------------------------------
  cvLsInitializeCounters

  This routine resets all counters from an CVLsMem structure.
  -----------------------------------------------------------------*/
pub fn cvLsInitializeCounters(cvls_mem: &mut CVLsMem) -> i32 {
    cvls_mem.nje = 0;
    cvls_mem.nfeDQ = 0;
    cvls_mem.nstlj = 0;
    cvls_mem.npe = 0;
    cvls_mem.nli = 0;
    cvls_mem.nps = 0;
    cvls_mem.ncfl = 0;
    cvls_mem.njtsetup = 0;
    cvls_mem.njtimes = 0;
    0
}

/*================================================================
  PART II - backward problems
  ================================================================*/

/*---------------------------------------------------------------
  CVSLS Exported functions -- Required
  ---------------------------------------------------------------*/

/* CVodeSetLinearSolverB specifies the linear solver for backward
   integration */
pub fn CVodeSetLinearSolverB(
    cv_mem: &mut CVodeMem,
    which: i32,
    LS: LinearSolver,
    A: Option<SUNMatrix>,
) -> i32 {
    /* Was ASA initialized? */
    if !cv_mem.cv_adjMallocDone {
        cvProcessError(Some(cv_mem), CVLS_NO_ADJ, line!(), "CVodeSetLinearSolverB", file!(),
                       MSG_LS_NO_ADJ);
        return CVLS_NO_ADJ;
    }

    /* Check which */
    if which >= cv_mem.cv_adj_mem.as_ref().unwrap().ca_nbckpbs {
        cvProcessError(Some(cv_mem), CVLS_ILL_INPUT, line!(), "CVodeSetLinearSolverB", file!(),
                       MSG_LS_BAD_WHICH);
        return CVLS_ILL_INPUT;
    }
    let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();

    /* Find the CVodeBMem entry corresponding to which (C: linked-list
       search by cv_index; backward problems are created with
       cv_index = 0..nbckpbs-1, so the entry always exists) */
    let idx = ca_mem.cvB_mem.iter().position(|b| b.cv_index == which).unwrap();
    let cvB_mem = &mut ca_mem.cvB_mem[idx];

    /* Get memory for CVLsMemRecB and attach it (lmemB); this also frees
       any existing system solver attached to cvB (C: cv_lfree; RAII on
       overwrite here). The C lfreeB attachment is cvLsFreeB, subsumed
       by Drop. */
    cvB_mem.cv_lmem = Some(Box::new(CVLsMemB::default()));

    /* set the linear solver for this backward problem */
    let retval = CVodeSetLinearSolver(&mut cvB_mem.cv_mem, LS, A);
    if retval != CVLS_SUCCESS {
        cvB_mem.cv_lmem = None;
    }

    retval
}

/*---------------------------------------------------------------
  CVSLS Exported functions -- Optional input/output
  ---------------------------------------------------------------*/

pub fn CVodeSetJacFnB(cv_mem: &mut CVodeMem, which: i32, jacB: Option<CVLsJacFnB>) -> i32 {
    /* access relevant memory structures */
    let idx = match cvLs_AccessLMemB(cv_mem, which, "CVodeSetJacFnB") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let cvB_mem = &mut cv_mem.cv_adj_mem.as_mut().unwrap().cvB_mem[idx];

    /* set jacB function pointer */
    cvLsB_downcast(cvB_mem).jacB = jacB;

    /* call corresponding routine for cvodeB_mem structure */
    if jacB.is_some() {
        CVodeSetJacFn(&mut cvB_mem.cv_mem, Some(cvLsJacBWrapper))
    } else {
        CVodeSetJacFn(&mut cvB_mem.cv_mem, None)
    }
}

pub fn CVodeSetJacFnBS(cv_mem: &mut CVodeMem, which: i32, jacBS: Option<CVLsJacFnBS>) -> i32 {
    /* access relevant memory structures */
    let idx = match cvLs_AccessLMemB(cv_mem, which, "CVodeSetJacFnBS") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let cvB_mem = &mut cv_mem.cv_adj_mem.as_mut().unwrap().cvB_mem[idx];

    /* set jacBS function pointer */
    cvLsB_downcast(cvB_mem).jacBS = jacBS;

    /* call corresponding routine for cvodeB_mem structure */
    if jacBS.is_some() {
        CVodeSetJacFn(&mut cvB_mem.cv_mem, Some(cvLsJacBSWrapper))
    } else {
        CVodeSetJacFn(&mut cvB_mem.cv_mem, None)
    }
}

pub fn CVodeSetEpsLinB(cv_mem: &mut CVodeMem, which: i32, eplifacB: f64) -> i32 {
    /* access relevant memory structures */
    let idx = match cvLs_AccessLMemB(cv_mem, which, "CVodeSetEpsLinB") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let cvB_mem = &mut cv_mem.cv_adj_mem.as_mut().unwrap().cvB_mem[idx];

    /* call corresponding routine for cvodeB_mem structure */
    CVodeSetEpsLin(&mut cvB_mem.cv_mem, eplifacB)
}

pub fn CVodeSetLSNormFactorB(cv_mem: &mut CVodeMem, which: i32, nrmfacB: f64) -> i32 {
    /* access relevant memory structures */
    let idx = match cvLs_AccessLMemB(cv_mem, which, "CVodeSetLSNormFactorB") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let cvB_mem = &mut cv_mem.cv_adj_mem.as_mut().unwrap().cvB_mem[idx];

    /* call corresponding routine for cvodeB_mem structure */
    CVodeSetLSNormFactor(&mut cvB_mem.cv_mem, nrmfacB)
}

pub fn CVodeSetLinearSolutionScalingB(cv_mem: &mut CVodeMem, which: i32, onoffB: bool) -> i32 {
    /* access relevant memory structures */
    let idx = match cvLs_AccessLMemB(cv_mem, which, "CVodeSetLinearSolutionScalingB") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let cvB_mem = &mut cv_mem.cv_adj_mem.as_mut().unwrap().cvB_mem[idx];

    /* call corresponding routine for cvodeB_mem structure */
    CVodeSetLinearSolutionScaling(&mut cvB_mem.cv_mem, onoffB)
}

pub fn CVodeSetPreconditionerB(
    cv_mem: &mut CVodeMem,
    which: i32,
    psetupB: Option<CVLsPrecSetupFnB>,
    psolveB: Option<CVLsPrecSolveFnB>,
) -> i32 {
    /* access relevant memory structures */
    let idx = match cvLs_AccessLMemB(cv_mem, which, "CVodeSetPreconditionerB") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let cvB_mem = &mut cv_mem.cv_adj_mem.as_mut().unwrap().cvB_mem[idx];

    /* Set preconditioners for the backward problem. */
    {
        let cvlsB_mem = cvLsB_downcast(cvB_mem);
        cvlsB_mem.psetB = psetupB;
        cvlsB_mem.psolveB = psolveB;
    }

    /* Call the corresponding "set" routine for the backward problem */
    let cvls_psetup: Option<CVLsPrecSetupFn> =
        if psetupB.is_none() { None } else { Some(cvLsPrecSetupBWrapper) };
    let cvls_psolve: Option<CVLsPrecSolveFn> =
        if psolveB.is_none() { None } else { Some(cvLsPrecSolveBWrapper) };
    CVodeSetPreconditioner(&mut cvB_mem.cv_mem, cvls_psetup, cvls_psolve)
}

pub fn CVodeSetPreconditionerBS(
    cv_mem: &mut CVodeMem,
    which: i32,
    psetupBS: Option<CVLsPrecSetupFnBS>,
    psolveBS: Option<CVLsPrecSolveFnBS>,
) -> i32 {
    /* access relevant memory structures */
    let idx = match cvLs_AccessLMemB(cv_mem, which, "CVodeSetPreconditionerBS") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let cvB_mem = &mut cv_mem.cv_adj_mem.as_mut().unwrap().cvB_mem[idx];

    /* Set preconditioners for the backward problem. */
    {
        let cvlsB_mem = cvLsB_downcast(cvB_mem);
        cvlsB_mem.psetBS = psetupBS;
        cvlsB_mem.psolveBS = psolveBS;
    }

    /* Call the corresponding "set" routine for the backward problem */
    let cvls_psetup: Option<CVLsPrecSetupFn> =
        if psetupBS.is_none() { None } else { Some(cvLsPrecSetupBSWrapper) };
    let cvls_psolve: Option<CVLsPrecSolveFn> =
        if psolveBS.is_none() { None } else { Some(cvLsPrecSolveBSWrapper) };
    CVodeSetPreconditioner(&mut cvB_mem.cv_mem, cvls_psetup, cvls_psolve)
}

pub fn CVodeSetJacTimesB(
    cv_mem: &mut CVodeMem,
    which: i32,
    jtsetupB: Option<CVLsJacTimesSetupFnB>,
    jtimesB: Option<CVLsJacTimesVecFnB>,
) -> i32 {
    /* access relevant memory structures */
    let idx = match cvLs_AccessLMemB(cv_mem, which, "CVodeSetJacTimesB") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let cvB_mem = &mut cv_mem.cv_adj_mem.as_mut().unwrap().cvB_mem[idx];

    /* Set jacobian routines for the backward problem. */
    {
        let cvlsB_mem = cvLsB_downcast(cvB_mem);
        cvlsB_mem.jtsetupB = jtsetupB;
        cvlsB_mem.jtimesB = jtimesB;
    }

    /* Call the corresponding "set" routine for the backward problem */
    let cvls_jtsetup: Option<CVLsJacTimesSetupFn> =
        if jtsetupB.is_none() { None } else { Some(cvLsJacTimesSetupBWrapper) };
    let cvls_jtimes: Option<CVLsJacTimesVecFn> =
        if jtimesB.is_none() { None } else { Some(cvLsJacTimesVecBWrapper) };
    CVodeSetJacTimes(&mut cvB_mem.cv_mem, cvls_jtsetup, cvls_jtimes)
}

pub fn CVodeSetJacTimesBS(
    cv_mem: &mut CVodeMem,
    which: i32,
    jtsetupBS: Option<CVLsJacTimesSetupFnBS>,
    jtimesBS: Option<CVLsJacTimesVecFnBS>,
) -> i32 {
    /* access relevant memory structures */
    let idx = match cvLs_AccessLMemB(cv_mem, which, "CVodeSetJacTimesBS") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let cvB_mem = &mut cv_mem.cv_adj_mem.as_mut().unwrap().cvB_mem[idx];

    /* Set jacobian routines for the backward problem. */
    {
        let cvlsB_mem = cvLsB_downcast(cvB_mem);
        cvlsB_mem.jtsetupBS = jtsetupBS;
        cvlsB_mem.jtimesBS = jtimesBS;
    }

    /* Call the corresponding "set" routine for the backward problem */
    let cvls_jtsetup: Option<CVLsJacTimesSetupFn> =
        if jtsetupBS.is_none() { None } else { Some(cvLsJacTimesSetupBSWrapper) };
    let cvls_jtimes: Option<CVLsJacTimesVecFn> =
        if jtimesBS.is_none() { None } else { Some(cvLsJacTimesVecBSWrapper) };
    CVodeSetJacTimes(&mut cvB_mem.cv_mem, cvls_jtsetup, cvls_jtimes)
}

pub fn CVodeSetJacTimesRhsFnB(
    cv_mem: &mut CVodeMem,
    which: i32,
    jtimesRhsFn: Option<CVRhsFn>,
) -> i32 {
    /* access relevant memory structures */
    let idx = match cvLs_AccessLMemB(cv_mem, which, "CVodeSetJacTimesRhsFnB") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let cvB_mem = &mut cv_mem.cv_adj_mem.as_mut().unwrap().cvB_mem[idx];

    /* Call the corresponding "set" routine for the backward problem */
    CVodeSetJacTimesRhsFn(&mut cvB_mem.cv_mem, jtimesRhsFn)
}

pub fn CVodeSetLinSysFnB(cv_mem: &mut CVodeMem, which: i32, linsysB: Option<CVLsLinSysFnB>) -> i32 {
    /* access relevant memory structures */
    let idx = match cvLs_AccessLMemB(cv_mem, which, "CVodeSetLinSysFnB") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let cvB_mem = &mut cv_mem.cv_adj_mem.as_mut().unwrap().cvB_mem[idx];

    /* set linsysB function pointer */
    cvLsB_downcast(cvB_mem).linsysB = linsysB;

    /* call corresponding routine for cvodeB_mem structure */
    if linsysB.is_some() {
        CVodeSetLinSysFn(&mut cvB_mem.cv_mem, Some(cvLsLinSysBWrapper))
    } else {
        CVodeSetLinSysFn(&mut cvB_mem.cv_mem, None)
    }
}

pub fn CVodeSetLinSysFnBS(
    cv_mem: &mut CVodeMem,
    which: i32,
    linsysBS: Option<CVLsLinSysFnBS>,
) -> i32 {
    /* access relevant memory structures */
    let idx = match cvLs_AccessLMemB(cv_mem, which, "CVodeSetLinSysFnBS") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let cvB_mem = &mut cv_mem.cv_adj_mem.as_mut().unwrap().cvB_mem[idx];

    /* set linsysBS function pointer */
    cvLsB_downcast(cvB_mem).linsysBS = linsysBS;

    /* call corresponding routine for cvodeB_mem structure */
    if linsysBS.is_some() {
        CVodeSetLinSysFn(&mut cvB_mem.cv_mem, Some(cvLsLinSysBSWrapper))
    } else {
        CVodeSetLinSysFn(&mut cvB_mem.cv_mem, None)
    }
}

/*-----------------------------------------------------------------
  CVSLS private functions for backwards problems
  -----------------------------------------------------------------*/

/* Downcast helper for the CVLsMemB stored in CVodeBMem.cv_lmem
   (the C `(CVLsMemB)(cvB_mem->cv_lmem)` cast). Callers reach this
   only after cvLs_AccessLMemB(BCur) verified the attachment. */
fn cvLsB_downcast(cvB_mem: &mut CVodeBMem) -> &mut CVLsMemB {
    cvB_mem.cv_lmem.as_mut().unwrap().downcast_mut::<CVLsMemB>().unwrap()
}

/* The C wrappers receive the *forward* integrator memory through the
   void* user_data of the backward problem (cvodea.c: CVodeCreateB
   does CVodeSetUserData(cvodeB_mem, cvode_mem)). This is the Rust
   analogue of the `(CVodeMem)cvode_mem` cast: the cvodea.rs backward
   driver must place the outer CVodeMem in the nested problem's
   cv_user_data for the duration of the backward run (forward
   reference / cross-module contract with cvodea.rs). */
fn cvLs_AccessCVodeMem<'a>(
    cvode_mem: &'a mut UserData,
    fname: &str,
) -> Result<&'a mut CVodeMem, i32> {
    match cvode_mem.as_mut().and_then(|d| d.downcast_mut::<CVodeMem>()) {
        Some(m) => Ok(m),
        None => {
            cvProcessError(None, CVLS_MEM_NULL, line!(), fname, file!(), MSG_LS_CVMEM_NULL);
            Err(CVLS_MEM_NULL)
        }
    }
}

/* ca_mem->ca_IMget(cv_mem, t, y, yS): dispatch on ca_IMtype per the
   pinned interpolation-module design (cvodes_impl.rs decision 4). An
   empty yS Vec plays the role of the C NULL argument.
   FORWARD REFERENCE: cvaHermiteGetY / cvaPolynomialGetY are
   implemented by the cvodea.c port (cvodea.rs). */
fn cvLsIMget(cv_mem: &mut CVodeMem, t: f64, y: &mut NVector, yS: &mut Vec<NVector>) -> i32 {
    if cv_mem.cv_adj_mem.as_ref().unwrap().ca_IMtype == CV_HERMITE {
        crate::cvodea::cvaHermiteGetY(cv_mem, t, y, yS)
    } else {
        crate::cvodea::cvaPolynomialGetY(cv_mem, t, y, yS)
    }
}

/* cvLsJacBWrapper interfaces to the CVLsJacFnB routine provided
   by the user. cvLsJacBWrapper is of type CVLsJacFn. */
fn cvLsJacBWrapper(
    t: f64,
    yB: &NVector,
    fyB: &NVector,
    JB: &mut SUNMatrix,
    cvode_mem: &mut UserData,
    tmp1B: &mut NVector,
    tmp2B: &mut NVector,
    tmp3B: &mut NVector,
) -> i32 {
    /* access relevant memory structures */
    let cv_mem = match cvLs_AccessCVodeMem(cvode_mem, "cvLsJacBWrapper") {
        Ok(m) => m,
        Err(e) => return e,
    };
    let which = match cvLs_AccessLMemBCur(cv_mem, "cvLsJacBWrapper") {
        Ok(w) => w,
        Err(e) => return e,
    };

    /* Forward solution from interpolation */
    let mut ytmp = std::mem::take(&mut cv_mem.cv_adj_mem.as_mut().unwrap().ca_ytmp);
    let mut noS: Vec<NVector> = Vec::new();
    let retval = cvLsIMget(cv_mem, t, &mut ytmp, &mut noS);
    if retval != CV_SUCCESS {
        cv_mem.cv_adj_mem.as_mut().unwrap().ca_ytmp = ytmp;
        cvProcessError(Some(cv_mem), -1, line!(), "cvLsJacBWrapper", file!(), MSG_LS_BAD_TINTERP);
        return -1;
    }

    /* Call user's adjoint jacB routine (of type CVLsJacFnB) */
    let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
    let cvB_mem = &mut ca_mem.cvB_mem[which];
    let jacB = cvLsB_downcast(cvB_mem).jacB.unwrap();
    let retval = jacB(t, &ytmp, yB, fyB, JB, &mut cvB_mem.cv_user_data, tmp1B, tmp2B, tmp3B);
    ca_mem.ca_ytmp = ytmp;
    retval
}

/* cvLsJacBSWrapper interfaces to the CVLsJacFnBS routine provided
   by the user. cvLsJacBSWrapper is of type CVLsJacFn. */
fn cvLsJacBSWrapper(
    t: f64,
    yB: &NVector,
    fyB: &NVector,
    JB: &mut SUNMatrix,
    cvode_mem: &mut UserData,
    tmp1B: &mut NVector,
    tmp2B: &mut NVector,
    tmp3B: &mut NVector,
) -> i32 {
    /* access relevant memory structures */
    let cv_mem = match cvLs_AccessCVodeMem(cvode_mem, "cvLsJacBSWrapper") {
        Ok(m) => m,
        Err(e) => return e,
    };
    let which = match cvLs_AccessLMemBCur(cv_mem, "cvLsJacBSWrapper") {
        Ok(w) => w,
        Err(e) => return e,
    };

    /* Forward solution from interpolation */
    let (mut ytmp, mut yStmp, interp_sensi) = {
        let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
        (
            std::mem::take(&mut ca_mem.ca_ytmp),
            std::mem::take(&mut ca_mem.ca_yStmp),
            ca_mem.ca_IMinterpSensi,
        )
    };
    let retval = if interp_sensi {
        cvLsIMget(cv_mem, t, &mut ytmp, &mut yStmp)
    } else {
        let mut noS: Vec<NVector> = Vec::new();
        cvLsIMget(cv_mem, t, &mut ytmp, &mut noS)
    };
    if retval != CV_SUCCESS {
        let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
        ca_mem.ca_ytmp = ytmp;
        ca_mem.ca_yStmp = yStmp;
        cvProcessError(Some(cv_mem), -1, line!(), "cvLsJacBSWrapper", file!(), MSG_LS_BAD_TINTERP);
        return -1;
    }

    /* Call user's adjoint dense djacBS routine (of type CVLsDenseJacFnBS) */
    let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
    let cvB_mem = &mut ca_mem.cvB_mem[which];
    let jacBS = cvLsB_downcast(cvB_mem).jacBS.unwrap();
    let retval = jacBS(
        t, &ytmp, &yStmp, yB, fyB, JB, &mut cvB_mem.cv_user_data, tmp1B, tmp2B, tmp3B,
    );
    ca_mem.ca_ytmp = ytmp;
    ca_mem.ca_yStmp = yStmp;
    retval
}

/* cvLsPrecSetupBWrapper interfaces to the CVLsPrecSetupFnB routine
   provided by the user */
fn cvLsPrecSetupBWrapper(
    t: f64,
    yB: &NVector,
    fyB: &NVector,
    jokB: bool,
    jcurPtrB: &mut bool,
    gammaB: f64,
    cvode_mem: &mut UserData,
) -> i32 {
    /* access relevant memory structures */
    let cv_mem = match cvLs_AccessCVodeMem(cvode_mem, "cvLsPrecSetupBWrapper") {
        Ok(m) => m,
        Err(e) => return e,
    };
    let which = match cvLs_AccessLMemBCur(cv_mem, "cvLsPrecSetupBWrapper") {
        Ok(w) => w,
        Err(e) => return e,
    };

    /* Get forward solution from interpolation */
    let mut ytmp = std::mem::take(&mut cv_mem.cv_adj_mem.as_mut().unwrap().ca_ytmp);
    let mut noS: Vec<NVector> = Vec::new();
    let retval = cvLsIMget(cv_mem, t, &mut ytmp, &mut noS);
    if retval != CV_SUCCESS {
        cv_mem.cv_adj_mem.as_mut().unwrap().ca_ytmp = ytmp;
        cvProcessError(Some(cv_mem), -1, line!(), "cvLsPrecSetupBWrapper", file!(),
                       MSG_LS_BAD_TINTERP);
        return -1;
    }

    /* Call user's adjoint precondB routine */
    let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
    let cvB_mem = &mut ca_mem.cvB_mem[which];
    let psetB = cvLsB_downcast(cvB_mem).psetB.unwrap();
    let retval = psetB(t, &ytmp, yB, fyB, jokB, jcurPtrB, gammaB, &mut cvB_mem.cv_user_data);
    ca_mem.ca_ytmp = ytmp;
    retval
}

/* cvLsPrecSetupBSWrapper interfaces to the CVLsPrecSetupFnBS routine
   provided by the user */
fn cvLsPrecSetupBSWrapper(
    t: f64,
    yB: &NVector,
    fyB: &NVector,
    jokB: bool,
    jcurPtrB: &mut bool,
    gammaB: f64,
    cvode_mem: &mut UserData,
) -> i32 {
    /* access relevant memory structures */
    let cv_mem = match cvLs_AccessCVodeMem(cvode_mem, "cvLsPrecSetupBSWrapper") {
        Ok(m) => m,
        Err(e) => return e,
    };
    let which = match cvLs_AccessLMemBCur(cv_mem, "cvLsPrecSetupBSWrapper") {
        Ok(w) => w,
        Err(e) => return e,
    };

    /* Forward solution from interpolation */
    let (mut ytmp, mut yStmp, interp_sensi) = {
        let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
        (
            std::mem::take(&mut ca_mem.ca_ytmp),
            std::mem::take(&mut ca_mem.ca_yStmp),
            ca_mem.ca_IMinterpSensi,
        )
    };
    let retval = if interp_sensi {
        cvLsIMget(cv_mem, t, &mut ytmp, &mut yStmp)
    } else {
        let mut noS: Vec<NVector> = Vec::new();
        cvLsIMget(cv_mem, t, &mut ytmp, &mut noS)
    };
    if retval != CV_SUCCESS {
        let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
        ca_mem.ca_ytmp = ytmp;
        ca_mem.ca_yStmp = yStmp;
        cvProcessError(Some(cv_mem), -1, line!(), "cvLsPrecSetupBSWrapper", file!(),
                       MSG_LS_BAD_TINTERP);
        return -1;
    }

    /* Call user's adjoint precondB routine */
    let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
    let cvB_mem = &mut ca_mem.cvB_mem[which];
    let psetBS = cvLsB_downcast(cvB_mem).psetBS.unwrap();
    let retval = psetBS(
        t, &ytmp, &yStmp, yB, fyB, jokB, jcurPtrB, gammaB, &mut cvB_mem.cv_user_data,
    );
    ca_mem.ca_ytmp = ytmp;
    ca_mem.ca_yStmp = yStmp;
    retval
}

/* cvLsPrecSolveBWrapper interfaces to the CVLsPrecSolveFnB routine
   provided by the user */
fn cvLsPrecSolveBWrapper(
    t: f64,
    yB: &NVector,
    fyB: &NVector,
    rB: &NVector,
    zB: &mut NVector,
    gammaB: f64,
    deltaB: f64,
    lrB: i32,
    cvode_mem: &mut UserData,
) -> i32 {
    /* access relevant memory structures */
    let cv_mem = match cvLs_AccessCVodeMem(cvode_mem, "cvLsPrecSolveBWrapper") {
        Ok(m) => m,
        Err(e) => return e,
    };
    let which = match cvLs_AccessLMemBCur(cv_mem, "cvLsPrecSolveBWrapper") {
        Ok(w) => w,
        Err(e) => return e,
    };

    /* Forward solution from interpolation */
    let mut ytmp = std::mem::take(&mut cv_mem.cv_adj_mem.as_mut().unwrap().ca_ytmp);
    let mut noS: Vec<NVector> = Vec::new();
    let retval = cvLsIMget(cv_mem, t, &mut ytmp, &mut noS);
    if retval != CV_SUCCESS {
        cv_mem.cv_adj_mem.as_mut().unwrap().ca_ytmp = ytmp;
        cvProcessError(Some(cv_mem), -1, line!(), "cvLsPrecSolveBWrapper", file!(),
                       MSG_LS_BAD_TINTERP);
        return -1;
    }

    /* Call user's adjoint psolveB routine */
    let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
    let cvB_mem = &mut ca_mem.cvB_mem[which];
    let psolveB = cvLsB_downcast(cvB_mem).psolveB.unwrap();
    let retval = psolveB(
        t, &ytmp, yB, fyB, rB, zB, gammaB, deltaB, lrB, &mut cvB_mem.cv_user_data,
    );
    ca_mem.ca_ytmp = ytmp;
    retval
}

/* cvLsPrecSolveBSWrapper interfaces to the CVLsPrecSolveFnBS routine
   provided by the user */
fn cvLsPrecSolveBSWrapper(
    t: f64,
    yB: &NVector,
    fyB: &NVector,
    rB: &NVector,
    zB: &mut NVector,
    gammaB: f64,
    deltaB: f64,
    lrB: i32,
    cvode_mem: &mut UserData,
) -> i32 {
    /* access relevant memory structures */
    let cv_mem = match cvLs_AccessCVodeMem(cvode_mem, "cvLsPrecSolveBSWrapper") {
        Ok(m) => m,
        Err(e) => return e,
    };
    let which = match cvLs_AccessLMemBCur(cv_mem, "cvLsPrecSolveBSWrapper") {
        Ok(w) => w,
        Err(e) => return e,
    };

    /* Forward solution from interpolation */
    let (mut ytmp, mut yStmp, interp_sensi) = {
        let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
        (
            std::mem::take(&mut ca_mem.ca_ytmp),
            std::mem::take(&mut ca_mem.ca_yStmp),
            ca_mem.ca_IMinterpSensi,
        )
    };
    let retval = if interp_sensi {
        cvLsIMget(cv_mem, t, &mut ytmp, &mut yStmp)
    } else {
        let mut noS: Vec<NVector> = Vec::new();
        cvLsIMget(cv_mem, t, &mut ytmp, &mut noS)
    };
    if retval != CV_SUCCESS {
        let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
        ca_mem.ca_ytmp = ytmp;
        ca_mem.ca_yStmp = yStmp;
        cvProcessError(Some(cv_mem), -1, line!(), "cvLsPrecSolveBSWrapper", file!(),
                       MSG_LS_BAD_TINTERP);
        return -1;
    }

    /* Call user's adjoint psolveBS routine */
    let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
    let cvB_mem = &mut ca_mem.cvB_mem[which];
    let psolveBS = cvLsB_downcast(cvB_mem).psolveBS.unwrap();
    let retval = psolveBS(
        t, &ytmp, &yStmp, yB, fyB, rB, zB, gammaB, deltaB, lrB, &mut cvB_mem.cv_user_data,
    );
    ca_mem.ca_ytmp = ytmp;
    ca_mem.ca_yStmp = yStmp;
    retval
}

/* cvLsJacTimesSetupBWrapper interfaces to the CVLsJacTimesSetupFnB
   routine provided by the user */
fn cvLsJacTimesSetupBWrapper(t: f64, yB: &NVector, fyB: &NVector, cvode_mem: &mut UserData) -> i32 {
    /* access relevant memory structures */
    let cv_mem = match cvLs_AccessCVodeMem(cvode_mem, "cvLsJacTimesSetupBWrapper") {
        Ok(m) => m,
        Err(e) => return e,
    };
    let which = match cvLs_AccessLMemBCur(cv_mem, "cvLsJacTimesSetupBWrapper") {
        Ok(w) => w,
        Err(e) => return e,
    };

    /* Forward solution from interpolation */
    let mut ytmp = std::mem::take(&mut cv_mem.cv_adj_mem.as_mut().unwrap().ca_ytmp);
    let mut noS: Vec<NVector> = Vec::new();
    let retval = cvLsIMget(cv_mem, t, &mut ytmp, &mut noS);
    if retval != CV_SUCCESS {
        cv_mem.cv_adj_mem.as_mut().unwrap().ca_ytmp = ytmp;
        cvProcessError(Some(cv_mem), -1, line!(), "cvLsJacTimesSetupBWrapper", file!(),
                       MSG_LS_BAD_TINTERP);
        return -1;
    }

    /* Call user's adjoint jtsetupB routine */
    let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
    let cvB_mem = &mut ca_mem.cvB_mem[which];
    let jtsetupB = cvLsB_downcast(cvB_mem).jtsetupB.unwrap();
    let retval = jtsetupB(t, &ytmp, yB, fyB, &mut cvB_mem.cv_user_data);
    ca_mem.ca_ytmp = ytmp;
    retval
}

/* cvLsJacTimesSetupBSWrapper interfaces to the CVLsJacTimesSetupFnBS
   routine provided by the user */
fn cvLsJacTimesSetupBSWrapper(
    t: f64,
    yB: &NVector,
    fyB: &NVector,
    cvode_mem: &mut UserData,
) -> i32 {
    /* access relevant memory structures */
    let cv_mem = match cvLs_AccessCVodeMem(cvode_mem, "cvLsJacTimesSetupBSWrapper") {
        Ok(m) => m,
        Err(e) => return e,
    };
    let which = match cvLs_AccessLMemBCur(cv_mem, "cvLsJacTimesSetupBSWrapper") {
        Ok(w) => w,
        Err(e) => return e,
    };

    /* Forward solution from interpolation */
    let (mut ytmp, mut yStmp, interp_sensi) = {
        let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
        (
            std::mem::take(&mut ca_mem.ca_ytmp),
            std::mem::take(&mut ca_mem.ca_yStmp),
            ca_mem.ca_IMinterpSensi,
        )
    };
    let retval = if interp_sensi {
        cvLsIMget(cv_mem, t, &mut ytmp, &mut yStmp)
    } else {
        let mut noS: Vec<NVector> = Vec::new();
        cvLsIMget(cv_mem, t, &mut ytmp, &mut noS)
    };
    if retval != CV_SUCCESS {
        let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
        ca_mem.ca_ytmp = ytmp;
        ca_mem.ca_yStmp = yStmp;
        cvProcessError(Some(cv_mem), -1, line!(), "cvLsJacTimesSetupBSWrapper", file!(),
                       MSG_LS_BAD_TINTERP);
        return -1;
    }

    /* Call user's adjoint jtsetupBS routine */
    let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
    let cvB_mem = &mut ca_mem.cvB_mem[which];
    let jtsetupBS = cvLsB_downcast(cvB_mem).jtsetupBS.unwrap();
    let retval = jtsetupBS(t, &ytmp, &yStmp, yB, fyB, &mut cvB_mem.cv_user_data);
    ca_mem.ca_ytmp = ytmp;
    ca_mem.ca_yStmp = yStmp;
    retval
}

/* cvLsJacTimesVecBWrapper interfaces to the CVLsJacTimesVecFnB routine
   provided by the user */
fn cvLsJacTimesVecBWrapper(
    vB: &NVector,
    JvB: &mut NVector,
    t: f64,
    yB: &NVector,
    fyB: &NVector,
    cvode_mem: &mut UserData,
    tmpB: &mut NVector,
) -> i32 {
    /* access relevant memory structures */
    let cv_mem = match cvLs_AccessCVodeMem(cvode_mem, "cvLsJacTimesVecBWrapper") {
        Ok(m) => m,
        Err(e) => return e,
    };
    let which = match cvLs_AccessLMemBCur(cv_mem, "cvLsJacTimesVecBWrapper") {
        Ok(w) => w,
        Err(e) => return e,
    };

    /* Forward solution from interpolation */
    let mut ytmp = std::mem::take(&mut cv_mem.cv_adj_mem.as_mut().unwrap().ca_ytmp);
    let mut noS: Vec<NVector> = Vec::new();
    let retval = cvLsIMget(cv_mem, t, &mut ytmp, &mut noS);
    if retval != CV_SUCCESS {
        cv_mem.cv_adj_mem.as_mut().unwrap().ca_ytmp = ytmp;
        cvProcessError(Some(cv_mem), -1, line!(), "cvLsJacTimesVecBWrapper", file!(),
                       MSG_LS_BAD_TINTERP);
        return -1;
    }

    /* Call user's adjoint jtimesB routine */
    let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
    let cvB_mem = &mut ca_mem.cvB_mem[which];
    let jtimesB = cvLsB_downcast(cvB_mem).jtimesB.unwrap();
    let retval = jtimesB(vB, JvB, t, &ytmp, yB, fyB, &mut cvB_mem.cv_user_data, tmpB);
    ca_mem.ca_ytmp = ytmp;
    retval
}

/* cvLsJacTimesVecBSWrapper interfaces to the CVLsJacTimesVecFnBS
   routine provided by the user */
fn cvLsJacTimesVecBSWrapper(
    vB: &NVector,
    JvB: &mut NVector,
    t: f64,
    yB: &NVector,
    fyB: &NVector,
    cvode_mem: &mut UserData,
    tmpB: &mut NVector,
) -> i32 {
    /* access relevant memory structures */
    let cv_mem = match cvLs_AccessCVodeMem(cvode_mem, "cvLsJacTimesVecBSWrapper") {
        Ok(m) => m,
        Err(e) => return e,
    };
    let which = match cvLs_AccessLMemBCur(cv_mem, "cvLsJacTimesVecBSWrapper") {
        Ok(w) => w,
        Err(e) => return e,
    };

    /* Forward solution from interpolation */
    let (mut ytmp, mut yStmp, interp_sensi) = {
        let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
        (
            std::mem::take(&mut ca_mem.ca_ytmp),
            std::mem::take(&mut ca_mem.ca_yStmp),
            ca_mem.ca_IMinterpSensi,
        )
    };
    let retval = if interp_sensi {
        cvLsIMget(cv_mem, t, &mut ytmp, &mut yStmp)
    } else {
        let mut noS: Vec<NVector> = Vec::new();
        cvLsIMget(cv_mem, t, &mut ytmp, &mut noS)
    };
    if retval != CV_SUCCESS {
        let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
        ca_mem.ca_ytmp = ytmp;
        ca_mem.ca_yStmp = yStmp;
        cvProcessError(Some(cv_mem), -1, line!(), "cvLsJacTimesVecBSWrapper", file!(),
                       MSG_LS_BAD_TINTERP);
        return -1;
    }

    /* Call user's adjoint jtimesBS routine */
    let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
    let cvB_mem = &mut ca_mem.cvB_mem[which];
    let jtimesBS = cvLsB_downcast(cvB_mem).jtimesBS.unwrap();
    let retval = jtimesBS(vB, JvB, t, &ytmp, &yStmp, yB, fyB, &mut cvB_mem.cv_user_data, tmpB);
    ca_mem.ca_ytmp = ytmp;
    ca_mem.ca_yStmp = yStmp;
    retval
}

/* cvLsLinSysBWrapper interfaces to the CVLsLinSysFnB routine provided
   by the user. cvLsLinSysBWrapper is of type CVLsLinSysFn. */
fn cvLsLinSysBWrapper(
    t: f64,
    yB: &NVector,
    fyB: &NVector,
    AB: &mut SUNMatrix,
    jokB: bool,
    jcurB: &mut bool,
    gammaB: f64,
    cvode_mem: &mut UserData,
    tmp1B: &mut NVector,
    tmp2B: &mut NVector,
    tmp3B: &mut NVector,
) -> i32 {
    /* access relevant memory structures */
    let cv_mem = match cvLs_AccessCVodeMem(cvode_mem, "cvLsLinSysBWrapper") {
        Ok(m) => m,
        Err(e) => return e,
    };
    let which = match cvLs_AccessLMemBCur(cv_mem, "cvLsLinSysBWrapper") {
        Ok(w) => w,
        Err(e) => return e,
    };

    /* Forward solution from interpolation */
    let mut ytmp = std::mem::take(&mut cv_mem.cv_adj_mem.as_mut().unwrap().ca_ytmp);
    let mut noS: Vec<NVector> = Vec::new();
    let retval = cvLsIMget(cv_mem, t, &mut ytmp, &mut noS);
    if retval != CV_SUCCESS {
        cv_mem.cv_adj_mem.as_mut().unwrap().ca_ytmp = ytmp;
        cvProcessError(Some(cv_mem), -1, line!(), "cvLsLinSysBWrapper", file!(),
                       MSG_LS_BAD_TINTERP);
        return -1;
    }

    /* Call user's adjoint linsysB routine (of type CVLsLinSysFnB) */
    let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
    let cvB_mem = &mut ca_mem.cvB_mem[which];
    let linsysB = cvLsB_downcast(cvB_mem).linsysB.unwrap();
    let retval = linsysB(
        t, &ytmp, yB, fyB, AB, jokB, jcurB, gammaB, &mut cvB_mem.cv_user_data, tmp1B, tmp2B,
        tmp3B,
    );
    ca_mem.ca_ytmp = ytmp;
    retval
}

/* cvLsLinSysBSWrapper interfaces to the CVLsLinSysFnBS routine provided
   by the user. cvLsLinSysBSWrapper is of type CVLsLinSysFn. */
fn cvLsLinSysBSWrapper(
    t: f64,
    yB: &NVector,
    fyB: &NVector,
    AB: &mut SUNMatrix,
    jokB: bool,
    jcurB: &mut bool,
    gammaB: f64,
    cvode_mem: &mut UserData,
    tmp1B: &mut NVector,
    tmp2B: &mut NVector,
    tmp3B: &mut NVector,
) -> i32 {
    /* access relevant memory structures */
    let cv_mem = match cvLs_AccessCVodeMem(cvode_mem, "cvLsLinSysBSWrapper") {
        Ok(m) => m,
        Err(e) => return e,
    };
    let which = match cvLs_AccessLMemBCur(cv_mem, "cvLsLinSysBSWrapper") {
        Ok(w) => w,
        Err(e) => return e,
    };

    /* Forward solution from interpolation */
    let (mut ytmp, mut yStmp, interp_sensi) = {
        let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
        (
            std::mem::take(&mut ca_mem.ca_ytmp),
            std::mem::take(&mut ca_mem.ca_yStmp),
            ca_mem.ca_IMinterpSensi,
        )
    };
    let retval = if interp_sensi {
        cvLsIMget(cv_mem, t, &mut ytmp, &mut yStmp)
    } else {
        let mut noS: Vec<NVector> = Vec::new();
        cvLsIMget(cv_mem, t, &mut ytmp, &mut noS)
    };
    if retval != CV_SUCCESS {
        let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
        ca_mem.ca_ytmp = ytmp;
        ca_mem.ca_yStmp = yStmp;
        cvProcessError(Some(cv_mem), -1, line!(), "cvLsLinSysBSWrapper", file!(),
                       MSG_LS_BAD_TINTERP);
        return -1;
    }

    /* Call user's adjoint dense djacBS routine (of type CVLsDenseJacFnBS) */
    let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
    let cvB_mem = &mut ca_mem.cvB_mem[which];
    let linsysBS = cvLsB_downcast(cvB_mem).linsysBS.unwrap();
    let retval = linsysBS(
        t, &ytmp, &yStmp, yB, fyB, AB, jokB, jcurB, gammaB, &mut cvB_mem.cv_user_data, tmp1B,
        tmp2B, tmp3B,
    );
    ca_mem.ca_ytmp = ytmp;
    ca_mem.ca_yStmp = yStmp;
    retval
}

/* cvLsFreeB frees memory associated with the CVSLS wrapper */
pub fn cvLsFreeB(cvB_mem: &mut CVodeBMem) -> i32 {
    /* free CVLsMemB interface structure (Rust: Drop) */
    cvB_mem.cv_lmem = None;
    CVLS_SUCCESS
}

/* cvLs_AccessLMemB checks the adjoint memory, `which`, and the CVLsMemB
   attachment, and returns the index of the backward problem in
   ca_mem.cvB_mem (the C version unpacks cv_mem/ca_mem/cvB_mem/cvlsB_mem
   pointers instead). If any is missing it returns CVLS_MEM_NULL,
   CVLS_NO_ADJ, CVLS_ILL_INPUT, or CVLS_LMEMB_NULL. */
pub(crate) fn cvLs_AccessLMemB(
    cv_mem: &mut CVodeMem,
    which: i32,
    fname: &str,
) -> Result<usize, i32> {
    /* access CVadjMem structure */
    if !cv_mem.cv_adjMallocDone {
        cvProcessError(Some(cv_mem), CVLS_NO_ADJ, line!(), fname, file!(), MSG_LS_NO_ADJ);
        return Err(CVLS_NO_ADJ);
    }
    let ca_mem = cv_mem.cv_adj_mem.as_ref().unwrap();

    /* Check which */
    if which >= ca_mem.ca_nbckpbs {
        cvProcessError(None, CVLS_ILL_INPUT, line!(), fname, file!(), MSG_LS_BAD_WHICH);
        return Err(CVLS_ILL_INPUT);
    }

    /* Find the CVodeBMem entry corresponding to which */
    let idx = ca_mem.cvB_mem.iter().position(|b| b.cv_index == which).unwrap();

    /* access CVLsMemB structure */
    let has_lmemB = ca_mem.cvB_mem[idx]
        .cv_lmem
        .as_ref()
        .map(|l| l.is::<CVLsMemB>())
        .unwrap_or(false);
    if !has_lmemB {
        cvProcessError(None, CVLS_LMEMB_NULL, line!(), fname, file!(), MSG_LS_LMEMB_NULL);
        return Err(CVLS_LMEMB_NULL);
    }

    Ok(idx)
}

/* cvLs_AccessLMemBCur returns the index of the currently integrated
   backward problem (ca_bckpbCrt) after the same checks. If any piece
   is missing it returns CVLS_MEM_NULL, CVLS_NO_ADJ, or
   CVLS_LMEMB_NULL. */
pub(crate) fn cvLs_AccessLMemBCur(cv_mem: &mut CVodeMem, fname: &str) -> Result<usize, i32> {
    /* access CVadjMem structure */
    if !cv_mem.cv_adjMallocDone {
        cvProcessError(Some(cv_mem), CVLS_NO_ADJ, line!(), fname, file!(), MSG_LS_NO_ADJ);
        return Err(CVLS_NO_ADJ);
    }
    let ca_mem = cv_mem.cv_adj_mem.as_ref().unwrap();

    /* get current backward problem */
    let which = match ca_mem.ca_bckpbCrt {
        Some(w) => w,
        None => {
            cvProcessError(None, CVLS_LMEMB_NULL, line!(), fname, file!(), MSG_LS_LMEMB_NULL);
            return Err(CVLS_LMEMB_NULL);
        }
    };

    /* access CVLsMemB structure */
    let has_lmemB = ca_mem.cvB_mem[which]
        .cv_lmem
        .as_ref()
        .map(|l| l.is::<CVLsMemB>())
        .unwrap_or(false);
    if !has_lmemB {
        cvProcessError(None, CVLS_LMEMB_NULL, line!(), fname, file!(), MSG_LS_LMEMB_NULL);
        return Err(CVLS_LMEMB_NULL);
    }

    Ok(which)
}

/*---------------------------------------------------------------
  EOF
  ---------------------------------------------------------------*/
