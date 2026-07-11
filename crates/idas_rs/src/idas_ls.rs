/* -----------------------------------------------------------------
 * Translated from src/idas/idas_ls.c (IDAS 7.7.0).
 * IDAS's linear solver interface (IDALS): attach a SUNLinearSolver,
 * Jacobian and preconditioner plumbing, difference-quotient Jacobian
 * (dense + band, with IDA's two-vector (yy,yp) cj-weighted
 * perturbations) and Jacobian-times-vector approximations, the
 * generic linit/lsetup/lsolve/lperf routines dispatched from
 * IDASolve, and the PART II backward-problem (adjoint) wrappers.
 * Structural donor: ida_rs/src/ida_ls.rs (verified Phase 4); the
 * PART I half is donor-verbatim except where noted (leniwLS = 34
 * per the IDAS C source, and the ida_bbdpre dispatch points are
 * deferred until the idas_bbdpre units land — see the PrecModule
 * comments below).  The interface hooks carry the IDALsMem
 * explicitly per the workspace take() convention:
 *
 *   idaLsInit (ida_mem, idals_mem)                          -> i32
 *   idaLsSetup(ida_mem, idals_mem, y, yp, r)                -> i32
 *   idaLsSolve(ida_mem, idals_mem, b, weight, ycur, ypcur,
 *              rescur)                                      -> i32
 *   idaLsPerf (ida_mem, idals_mem, perftask)
 *
 * (C: idaLsSetup additionally receives three tmp vectors — always
 * ida_tempv1/2/3, taken from the IDAMem fields here; idaLsSolve's
 * ycur/ypcur/rescur inputs replace the borrowed ycur/ypcur/rcur
 * pointers the C code parks inside IDALsMem, see idas_ls_impl.rs.)
 * -----------------------------------------------------------------*/
use std::cell::RefCell;

use crate::idas_bbdpre::{IDABBDPrecSetup, IDABBDPrecSolve};
use crate::idas_impl::*;
use crate::idas_ls_impl::*;
use crate::nvector_serial::*;
use crate::sundials_errors::*;
use crate::sundials_linearsolver::*;
use crate::sundials_math::*;
use crate::sundials_matrix::*;
use crate::sundials_types::*;
use crate::sundials_utils::fmt_g;

/* constants */
const MAX_ITERS: i32 = 3; /* max. number of attempts to recover in DQ J*v */
const ZERO: f64 = 0.0;
const PT25: f64 = 0.25;
const PT05: f64 = 0.05;
const PT9: f64 = 0.9;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;

/*---------------------------------------------------------------
  idaLs_AccessLMem

  This routine unpacks the idals_mem structure (IDALsMem) from the
  integrator memory.  If it is missing it returns IDALS_LMEM_NULL
  (the ida_mem NULL check of the C original cannot arise: ida_mem
  is a &mut IDAMem here).
  ---------------------------------------------------------------*/
pub fn idaLs_AccessLMem<'a>(
    ida_mem: &'a mut IDAMem,
    fname: &str,
) -> Result<&'a mut IDALsMem, i32> {
    match &mut ida_mem.ida_lmem {
        LsModule::Ls(ls) => Ok(ls),
        _ => {
            IDAProcessError(None, IDALS_LMEM_NULL, line!(), fname, file!(), MSG_LS_LMEM_NULL);
            Err(IDALS_LMEM_NULL)
        }
    }
}

/*================================================================
  PART I - forward problems
  ================================================================*/

/*---------------------------------------------------------------
  IDASetLinearSolver specifies the linear solver
  ---------------------------------------------------------------*/
pub fn IDASetLinearSolver(ida_mem: &mut IDAMem, LS: LinearSolver, A: Option<SUNMatrix>) -> i32 {
    /* (The C NULL-input and missing-ops checks — gettype/solve on the
       LS, nvconst/nvwrmsnorm/nvgetlength on the vector, resid/numiters/
       setatimes on iterative solvers — cannot fail here: every
       workspace LinearSolver and the serial NVector implement them.) */

    /* Retrieve the LS type */
    let ls_type = LS.ls_type();

    /* Set flags based on LS type */
    let iterative = ls_type != SUNLINEARSOLVER_DIRECT;
    let matrixbased =
        ls_type != SUNLINEARSOLVER_ITERATIVE && ls_type != SUNLINEARSOLVER_MATRIX_EMBEDDED;

    /* Ensure that A is NULL when LS is matrix-embedded */
    if ls_type == SUNLINEARSOLVER_MATRIX_EMBEDDED && A.is_some() {
        IDAProcessError(Some(ida_mem), IDALS_ILL_INPUT, line!(), "IDASetLinearSolver", file!(),
                        "Incompatible inputs: matrix-embedded LS requires NULL matrix");
        return IDALS_ILL_INPUT;
    }

    /* Check for compatible LS type, matrix and "atimes" support */
    if iterative {
        if matrixbased && A.is_none() {
            IDAProcessError(Some(ida_mem), IDALS_ILL_INPUT, line!(), "IDASetLinearSolver", file!(),
                            "Incompatible inputs: matrix-iterative LS requires non-NULL matrix");
            return IDALS_ILL_INPUT;
        }
    } else if A.is_none() {
        IDAProcessError(Some(ida_mem), IDALS_ILL_INPUT, line!(), "IDASetLinearSolver", file!(),
                        "Incompatible inputs: direct LS requires non-NULL matrix");
        return IDALS_ILL_INPUT;
    }

    /* free any existing system solver attached to IDA
       (RAII: dropped when ida_lmem is overwritten below) */

    /* In C the four system linear solver function fields (ida_linit /
       ida_lsetup/ida_lsolve/ida_lfree) and the iterative-only ida_lperf
       hook are installed here; the Rust port dispatches them through
       the LsModule::Ls variant, with the lperf dispatch guarded by the
       `iterative` flag stored below (idas_impl.rs contract). */

    /* Set defaults for Jacobian-related fields (A non-NULL => internal
       DQ Jacobian; the C jac = idaLsDQJac / J_data = IDA_mem pairing is
       jac = None + jacDQ here) */
    let jacDQ = A.is_some();

    /* Allocate memory for ytemp, yptemp and x (C: N_VClone(ida_tempv1);
       infallible here) */
    let n = ida_mem.ida_tempv1.len();

    /* Get memory for IDALsMemRec and set defaults (C memset(0) +
       explicit assignments) */
    let idals_mem = Box::new(IDALsMem {
        /* set SUNLinearSolver pointer */
        LS,

        /* Linear solver type information */
        iterative,
        matrixbased,

        /* Set defaults for Jacobian-related fields */
        J: A,
        jacDQ,
        jac: None, /* None + jacDQ => internal idaLsDQJac */

        /* Jacobian-times-vector: internal DQ by default */
        jtimesDQ: SUNTRUE,
        jtsetup: None,
        jtimes: None, /* None + jtimesDQ => internal idaLsDQJtimes */
        jt_res: ida_mem.ida_res,

        /* Set defaults for preconditioner-related fields
           (C pdata = ida_user_data, pfree = NULL) */
        pset: None,
        psolve: None,
        prec_module: PrecModule::None,

        /* Initialize counters (idaLsInitializeCounters) */
        nje: 0,
        nreDQ: 0,
        npe: 0,
        nli: 0,
        nps: 0,
        ncfl: 0,
        njtsetup: 0,
        njtimes: 0,

        /* memset(0) remainder */
        nst0: 0,
        nni0: 0,
        ncfn0: 0,
        ncfl0: 0,
        nwarn: 0,
        nstlj: 0,
        tnlj: 0.0,

        /* Set default values for the rest of the Ls parameters */
        eplifac: PT05,
        dqincfac: ONE,
        last_flag: IDALS_SUCCESS,

        /* If LS supports ATimes/preconditioning, the IDALs routines are
           supplied as closures at solve time (idaLsSolveIterative); the C
           SUNLinSolSetATimes / SUNLinSolSetPreconditioner registration
           calls have no Rust counterpart. */

        /* Allocate memory for ytemp, yptemp and x */
        ytemp: NVector::new(n),
        yptemp: NVector::new(n),
        x: NVector::new(n),

        /* For iterative LS, compute sqrtN (else memset(0)) */
        nrmfac: if iterative { SUNRsqrt(n as f64) } else { 0.0 },

        /* For matrix-based LS, enable solution scaling */
        scalesol: matrixbased,

        setup_disabled: SUNFALSE,
    });

    /* Attach linear solver memory to integrator memory */
    ida_mem.ida_lmem = LsModule::Ls(idals_mem);

    IDALS_SUCCESS
}

/*===============================================================
  Optional Set routines
  ===============================================================*/

/* IDASetJacFn specifies the Jacobian function */
pub fn IDASetJacFn(ida_mem: &mut IDAMem, jac: Option<IDALsJacFn>) -> i32 {
    /* access IDALsMem structure */
    let idals_mem = match idaLs_AccessLMem(ida_mem, "IDASetJacFn") {
        Ok(m) => m,
        Err(e) => return e,
    };

    /* return with failure if jac cannot be used */
    if jac.is_some() && idals_mem.J.is_none() {
        IDAProcessError(None, IDALS_ILL_INPUT, line!(), "IDASetJacFn", file!(),
                        "Jacobian routine cannot be supplied for NULL SUNMatrix");
        return IDALS_ILL_INPUT;
    }

    /* set Jacobian routine pointer, and update relevant flags */
    if let Some(j) = jac {
        idals_mem.jacDQ = SUNFALSE;
        idals_mem.jac = Some(j);
    } else {
        idals_mem.jacDQ = SUNTRUE;
        idals_mem.jac = None; /* internal idaLsDQJac */
    }

    IDALS_SUCCESS
}

/* IDASetEpsLin specifies the nonlinear -> linear tolerance scale factor */
pub fn IDASetEpsLin(ida_mem: &mut IDAMem, eplifac: f64) -> i32 {
    /* access IDALsMem structure */
    let idals_mem = match idaLs_AccessLMem(ida_mem, "IDASetEpsLin") {
        Ok(m) => m,
        Err(e) => return e,
    };

    /* Check for legal eplifac */
    if eplifac < ZERO {
        IDAProcessError(None, IDALS_ILL_INPUT, line!(), "IDASetEpsLin", file!(),
                        MSG_LS_NEG_EPLIFAC);
        return IDALS_ILL_INPUT;
    }

    idals_mem.eplifac = if eplifac == ZERO { PT05 } else { eplifac };

    IDALS_SUCCESS
}

/* IDASetWRMSNormFactor sets or computes the factor to use when converting from
   the integrator tolerance to the linear solver tolerance (WRMS to L2 norm). */
pub fn IDASetLSNormFactor(ida_mem: &mut IDAMem, nrmfac: f64) -> i32 {
    /* access IDALsMem structure */
    let idals_mem = match idaLs_AccessLMem(ida_mem, "IDASetLSNormFactor") {
        Ok(m) => m,
        Err(e) => return e,
    };

    if nrmfac > ZERO {
        /* user-provided factor */
        idals_mem.nrmfac = nrmfac;
    } else if nrmfac < ZERO {
        /* compute factor for WRMS norm with dot product */
        N_VConst(ONE, &mut idals_mem.ytemp);
        idals_mem.nrmfac = SUNRsqrt(N_VDotProd(&idals_mem.ytemp, &idals_mem.ytemp));
    } else {
        /* compute default factor for WRMS norm from vector length */
        idals_mem.nrmfac = SUNRsqrt(idals_mem.ytemp.len() as f64);
    }

    IDALS_SUCCESS
}

/* IDASetLinearSolutionScaling enables or disables scaling the linear solver
   solution to account for changes in cj. */
pub fn IDASetLinearSolutionScaling(ida_mem: &mut IDAMem, onoff: bool) -> i32 {
    /* access IDALsMem structure */
    let idals_mem = match idaLs_AccessLMem(ida_mem, "IDASetLinearSolutionScaling") {
        Ok(m) => m,
        Err(e) => return e,
    };

    /* check for valid solver type (no error message in C) */
    if !idals_mem.matrixbased {
        return IDALS_ILL_INPUT;
    }

    /* set solution scaling flag */
    idals_mem.scalesol = onoff;

    IDALS_SUCCESS
}

/* IDASetIncrementFactor specifies increment factor for DQ approximations to Jv */
pub fn IDASetIncrementFactor(ida_mem: &mut IDAMem, dqincfac: f64) -> i32 {
    /* access IDALsMem structure */
    let idals_mem = match idaLs_AccessLMem(ida_mem, "IDASetIncrementFactor") {
        Ok(m) => m,
        Err(e) => return e,
    };

    /* Check for legal dqincfac */
    if dqincfac <= ZERO {
        IDAProcessError(None, IDALS_ILL_INPUT, line!(), "IDASetIncrementFactor", file!(),
                        MSG_LS_NEG_DQINCFAC);
        return IDALS_ILL_INPUT;
    }

    idals_mem.dqincfac = dqincfac;

    IDALS_SUCCESS
}

/* IDASetPreconditioner specifies the user-supplied psetup and psolve routines */
pub fn IDASetPreconditioner(
    ida_mem: &mut IDAMem,
    psetup: Option<IDALsPrecSetupFn>,
    psolve: Option<IDALsPrecSolveFn>,
) -> i32 {
    /* access IDALsMem structure */
    let idals_mem = match idaLs_AccessLMem(ida_mem, "IDASetPreconditioner") {
        Ok(m) => m,
        Err(e) => return e,
    };

    /* store function pointers for user-supplied routines in IDALs interface */
    idals_mem.pset = psetup;
    idals_mem.psolve = psolve;
    idals_mem.prec_module = PrecModule::User;

    /* issue error if LS object does not allow user-supplied
       preconditioning (C: LS->ops->setpreconditioner == NULL, true of
       the direct and matrix-embedded solvers) */
    if matches!(
        idals_mem.LS,
        LinearSolver::Dense(_) | LinearSolver::Band(_) | LinearSolver::Custom(_)
    ) {
        IDAProcessError(None, IDALS_ILL_INPUT, line!(), "IDASetPreconditioner", file!(),
                        "SUNLinearSolver object does not support user-supplied preconditioning");
        return IDALS_ILL_INPUT;
    }

    /* notify iterative linear solver to call IDALs interface routines:
       done at solve time via the idaLsSolveIterative closures */

    IDALS_SUCCESS
}

/* IDASetJacTimes specifies the user-supplied Jacobian-vector product
   setup and multiply routines */
pub fn IDASetJacTimes(
    ida_mem: &mut IDAMem,
    jtsetup: Option<IDALsJacTimesSetupFn>,
    jtimes: Option<IDALsJacTimesVecFn>,
) -> i32 {
    let ida_res = ida_mem.ida_res;

    /* access IDALsMem structure */
    let idals_mem = match idaLs_AccessLMem(ida_mem, "IDASetJacTimes") {
        Ok(m) => m,
        Err(e) => return e,
    };

    /* issue error if LS object does not allow user-supplied ATimes
       (C: LS->ops->setatimes == NULL, true of the direct and
       matrix-embedded solvers) */
    if matches!(
        idals_mem.LS,
        LinearSolver::Dense(_) | LinearSolver::Band(_) | LinearSolver::Custom(_)
    ) {
        IDAProcessError(None, IDALS_ILL_INPUT, line!(), "IDASetJacTimes", file!(),
                        "SUNLinearSolver object does not support user-supplied ATimes routine");
        return IDALS_ILL_INPUT;
    }

    /* store function pointers for user-supplied routines in IDALs
       interface (NULL jtimes implies use of DQ default) */
    if let Some(jt) = jtimes {
        idals_mem.jtimesDQ = SUNFALSE;
        idals_mem.jtsetup = jtsetup;
        idals_mem.jtimes = Some(jt);
    } else {
        idals_mem.jtimesDQ = SUNTRUE;
        idals_mem.jtsetup = None;
        idals_mem.jtimes = None; /* internal idaLsDQJtimes */
        idals_mem.jt_res = ida_res;
    }

    IDALS_SUCCESS
}

/* IDASetJacTimesResFn specifies an alternative user-supplied DAE residual
   function to use in the internal finite difference Jacobian-vector
   product */
pub fn IDASetJacTimesResFn(ida_mem: &mut IDAMem, jtimesResFn: Option<IDAResFn>) -> i32 {
    let ida_res = ida_mem.ida_res;

    /* access IDALsMem structure */
    let idals_mem = match idaLs_AccessLMem(ida_mem, "IDASetJacTimesResFn") {
        Ok(m) => m,
        Err(e) => return e,
    };

    /* check if using internal finite difference approximation */
    if !idals_mem.jtimesDQ {
        IDAProcessError(None, IDALS_ILL_INPUT, line!(), "IDASetJacTimesResFn", file!(),
                        "Internal finite-difference Jacobian-vector product is disabled.");
        return IDALS_ILL_INPUT;
    }

    /* store function pointers for Res function (NULL implies use DAE Res) */
    match jtimesResFn {
        Some(f) => idals_mem.jt_res = Some(f),
        None => idals_mem.jt_res = ida_res,
    }

    IDALS_SUCCESS
}

/*===============================================================
  Optional Get routines
  ===============================================================*/

pub fn IDAGetJac<'a>(ida_mem: &'a mut IDAMem, j: &mut Option<&'a SUNMatrix>) -> i32 {
    /* access IDALsMem structure; set output and return */
    let idals_mem = match idaLs_AccessLMem(ida_mem, "IDAGetJac") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *j = idals_mem.J.as_ref();
    IDALS_SUCCESS
}

pub fn IDAGetJacCj(ida_mem: &mut IDAMem, cj_J: &mut f64) -> i32 {
    let cjold = ida_mem.ida_cjold;

    /* access IDALsMem structure; set output and return */
    if let Err(e) = idaLs_AccessLMem(ida_mem, "IDAGetJacCj") {
        return e;
    }
    *cj_J = cjold;
    IDALS_SUCCESS
}

pub fn IDAGetJacTime(ida_mem: &mut IDAMem, t_J: &mut f64) -> i32 {
    /* access IDALsMem structure; set output and return */
    let idals_mem = match idaLs_AccessLMem(ida_mem, "IDAGetJacTime") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *t_J = idals_mem.tnlj;
    IDALS_SUCCESS
}

pub fn IDAGetJacNumSteps(ida_mem: &mut IDAMem, nst_J: &mut i64) -> i32 {
    /* access IDALsMem structure; set output and return */
    let idals_mem = match idaLs_AccessLMem(ida_mem, "IDAGetJacNumSteps") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *nst_J = idals_mem.nstlj;
    IDALS_SUCCESS
}

/* IDAGetLinWorkSpace returns the length of workspace allocated
   for the IDALS linear solver interface */
pub fn IDAGetLinWorkSpace(ida_mem: &mut IDAMem, lenrwLS: &mut i64, leniwLS: &mut i64) -> i32 {
    let (lrw1, liw1) = N_VSpace(&ida_mem.ida_tempv1);

    /* access IDALsMem structure */
    let idals_mem = match idaLs_AccessLMem(ida_mem, "IDAGetLinWorkSpace") {
        Ok(m) => m,
        Err(e) => return e,
    };

    /* start with fixed sizes plus vector/matrix pointers
       (IDAS: 34 integer words; the IDA donor counts 33) */
    *lenrwLS = 3;
    *leniwLS = 34;

    /* add N_Vector sizes */
    *lenrwLS += 3 * lrw1;
    *leniwLS += 3 * liw1;

    /* add LS sizes */
    {
        let (lrw, liw) = idals_mem.LS.space();
        *lenrwLS += lrw;
        *leniwLS += liw;
    }

    IDALS_SUCCESS
}

/* IDAGetNumJacEvals returns the number of Jacobian evaluations */
pub fn IDAGetNumJacEvals(ida_mem: &mut IDAMem, njevals: &mut i64) -> i32 {
    /* access IDALsMem structure; store output and return */
    let idals_mem = match idaLs_AccessLMem(ida_mem, "IDAGetNumJacEvals") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *njevals = idals_mem.nje;
    IDALS_SUCCESS
}

/* IDAGetNumPrecEvals returns the number of preconditioner evaluations */
pub fn IDAGetNumPrecEvals(ida_mem: &mut IDAMem, npevals: &mut i64) -> i32 {
    /* access IDALsMem structure; store output and return */
    let idals_mem = match idaLs_AccessLMem(ida_mem, "IDAGetNumPrecEvals") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *npevals = idals_mem.npe;
    IDALS_SUCCESS
}

/* IDAGetNumPrecSolves returns the number of preconditioner solves */
pub fn IDAGetNumPrecSolves(ida_mem: &mut IDAMem, npsolves: &mut i64) -> i32 {
    /* access IDALsMem structure; store output and return */
    let idals_mem = match idaLs_AccessLMem(ida_mem, "IDAGetNumPrecSolves") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *npsolves = idals_mem.nps;
    IDALS_SUCCESS
}

/* IDAGetNumLinIters returns the number of linear iterations */
pub fn IDAGetNumLinIters(ida_mem: &mut IDAMem, nliters: &mut i64) -> i32 {
    /* access IDALsMem structure; store output and return */
    let idals_mem = match idaLs_AccessLMem(ida_mem, "IDAGetNumLinIters") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *nliters = idals_mem.nli;
    IDALS_SUCCESS
}

/* IDAGetNumLinConvFails returns the number of linear convergence failures */
pub fn IDAGetNumLinConvFails(ida_mem: &mut IDAMem, nlcfails: &mut i64) -> i32 {
    /* access IDALsMem structure; store output and return */
    let idals_mem = match idaLs_AccessLMem(ida_mem, "IDAGetNumLinConvFails") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *nlcfails = idals_mem.ncfl;
    IDALS_SUCCESS
}

/* IDAGetNumJTSetupEvals returns the number of calls to the
   user-supplied Jacobian-vector product setup routine */
pub fn IDAGetNumJTSetupEvals(ida_mem: &mut IDAMem, njtsetups: &mut i64) -> i32 {
    /* access IDALsMem structure; store output and return */
    let idals_mem = match idaLs_AccessLMem(ida_mem, "IDAGetNumJTSetupEvals") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *njtsetups = idals_mem.njtsetup;
    IDALS_SUCCESS
}

/* IDAGetNumJtimesEvals returns the number of calls to the
   Jacobian-vector product multiply routine */
pub fn IDAGetNumJtimesEvals(ida_mem: &mut IDAMem, njvevals: &mut i64) -> i32 {
    /* access IDALsMem structure; store output and return */
    let idals_mem = match idaLs_AccessLMem(ida_mem, "IDAGetNumJtimesEvals") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *njvevals = idals_mem.njtimes;
    IDALS_SUCCESS
}

/* IDAGetNumLinResEvals returns the number of calls to the DAE
   residual needed for the DQ Jacobian approximation or J*v
   product approximation */
pub fn IDAGetNumLinResEvals(ida_mem: &mut IDAMem, nrevalsLS: &mut i64) -> i32 {
    /* access IDALsMem structure; store output and return */
    let idals_mem = match idaLs_AccessLMem(ida_mem, "IDAGetNumLinResEvals") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *nrevalsLS = idals_mem.nreDQ;
    IDALS_SUCCESS
}

/* IDAGetLastLinFlag returns the last flag set in a IDALS function */
pub fn IDAGetLastLinFlag(ida_mem: &mut IDAMem, flag: &mut i64) -> i32 {
    /* access IDALsMem structure; store output and return */
    let idals_mem = match idaLs_AccessLMem(ida_mem, "IDAGetLastLinFlag") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *flag = idals_mem.last_flag as i64;
    IDALS_SUCCESS
}

/* IDAGetLinReturnFlagName translates from the integer error code
   returned by an IDALs routine to the corresponding string
   equivalent for that flag */
pub fn IDAGetLinReturnFlagName(flag: i64) -> String {
    match flag as i32 {
        IDALS_SUCCESS => "IDALS_SUCCESS",
        IDALS_MEM_NULL => "IDALS_MEM_NULL",
        IDALS_LMEM_NULL => "IDALS_LMEM_NULL",
        IDALS_ILL_INPUT => "IDALS_ILL_INPUT",
        IDALS_MEM_FAIL => "IDALS_MEM_FAIL",
        IDALS_PMEM_NULL => "IDALS_PMEM_NULL",
        IDALS_JACFUNC_UNRECVR => "IDALS_JACFUNC_UNRECVR",
        IDALS_JACFUNC_RECVR => "IDALS_JACFUNC_RECVR",
        IDALS_SUNMAT_FAIL => "IDALS_SUNMAT_FAIL",
        IDALS_SUNLS_FAIL => "IDALS_SUNLS_FAIL",
        _ => "NONE",
    }
    .to_string()
}

/*===============================================================
  IDASLS Private functions
  ===============================================================*/

/*---------------------------------------------------------------
  idaLsATimes:

  This routine generates the matrix-vector product z = Jv, where
  J is the system Jacobian, by calling either the user provided
  routine or the internal DQ routine.  The return value is
  the same as the value returned by jtimes -- 0 if successful,
  nonzero otherwise.  (The C reads ycur/ypcur/rcur out of IDALsMem;
  they are passed as arguments here, see idas_ls_impl.rs.
  idaLsSolveIterative replicates this body inside the ATimes
  closure handed to the iterative solvers.)
  ---------------------------------------------------------------*/
pub fn idaLsATimes(
    ida_mem: &mut IDAMem,
    idals_mem: &mut IDALsMem,
    ycur: &NVector,
    ypcur: &NVector,
    rcur: &NVector,
    v: &NVector,
    z: &mut NVector,
) -> i32 {
    /* call Jacobian-times-vector product routine
       (either user-supplied or internal DQ) */
    let retval = if idals_mem.jtimesDQ {
        let gmres = matches!(idals_mem.LS, LinearSolver::Spgmr(_) | LinearSolver::Spfgmr(_));
        let ewt = std::mem::take(&mut ida_mem.ida_ewt);
        let tn = ida_mem.ida_tn;
        let cj = ida_mem.ida_cj;
        let IDALsMem { ytemp, yptemp, nreDQ, nrmfac, dqincfac, jt_res, .. } = idals_mem;
        let r = idaLsDQJtimes(ida_mem, &ewt, nreDQ, gmres, *nrmfac, *dqincfac, *jt_res, tn,
                              ycur, ypcur, rcur, v, z, cj, ytemp, yptemp);
        ida_mem.ida_ewt = ewt;
        r
    } else {
        let jt = idals_mem.jtimes.unwrap();
        let IDALsMem { ytemp, yptemp, .. } = idals_mem;
        let IDAMem { ida_tn, ida_cj, ida_user_data, .. } = ida_mem;
        jt(*ida_tn, ycur, ypcur, rcur, v, z, *ida_cj, ida_user_data, ytemp, yptemp)
    };
    idals_mem.njtimes += 1;
    retval
}

/*---------------------------------------------------------------
  idaLsPSetup:

  This routine interfaces between the generic iterative linear
  solvers and the user's psetup routine.  It passes to psetup all
  required state information from ida_mem.  Its return value
  is the same as that returned by psetup.  In C the generic
  iterative linear solvers guarantee that idaLsPSetup will only
  be called in the case that the user's psetup routine is
  non-NULL; here the None case returns 0 (setup no-op).
  (C pdata = ida_user_data for the user PrecModule.  IDA's psetup
  carries no jok/jcurPtr — unlike CVODE's — so the donor's cv_jcur
  aliasing fix has no counterpart here.  The donor's
  PrecModule::BBDPre arm — C pset = IDABBDPrecSetup — lands with
  the idas_bbdpre units.)
  ---------------------------------------------------------------*/
pub fn idaLsPSetup(
    ida_mem: &mut IDAMem,
    idals_mem: &mut IDALsMem,
    ycur: &NVector,
    ypcur: &NVector,
    rcur: &NVector,
) -> i32 {
    let IDALsMem { pset, prec_module, npe, .. } = idals_mem;
    match prec_module {
        /* internal idas_bbdpre module (C pset = IDABBDPrecSetup) */
        PrecModule::BBDPre(pdata) => {
            let (tn, cj) = (ida_mem.ida_tn, ida_mem.ida_cj);
            let retval = IDABBDPrecSetup(ida_mem, pdata, tn, cj, ycur, ypcur);
            *npe += 1;
            retval
        }
        _ => {
            if let Some(pset) = *pset {
                /* Call user pset routine to update preconditioner. ida_ewt /
                   ida_hh are handed in directly: they are what a C user pset
                   fetches via IDAGetErrWeights / IDAGetCurrentStep (see
                   IDALsPrecSetupFn note). */
                let IDAMem { ida_tn, ida_cj, ida_ewt, ida_hh, ida_user_data, .. } = ida_mem;
                let retval =
                    pset(*ida_tn, ycur, ypcur, rcur, *ida_cj, ida_ewt, *ida_hh, ida_user_data);
                *npe += 1;
                retval
            } else {
                0
            }
        }
    }
}

/*---------------------------------------------------------------
  idaLsPSolve:

  This routine interfaces between the generic SUNLinSolSolve
  routine and the user's psolve routine.  It passes to psolve all
  required state information from ida_mem.  Its return value is
  the same as that returned by psolve.  The generic SUNLinSol
  solver guarantees that idaLsPSolve will not be called in the
  case in which preconditioning is not done; here the None case
  returns 0.  (idaLsSolveIterative replicates this body inside
  the PSolve closure handed to the iterative solvers; the lr
  input is unused, as in C.)
  ---------------------------------------------------------------*/
#[allow(clippy::too_many_arguments)]
pub fn idaLsPSolve(
    ida_mem: &mut IDAMem,
    idals_mem: &mut IDALsMem,
    ycur: &NVector,
    ypcur: &NVector,
    rcur: &NVector,
    r: &NVector,
    z: &mut NVector,
    tol: f64,
    _lr: i32,
) -> i32 {
    let IDALsMem { psolve, nps, .. } = idals_mem;
    if let Some(psolve) = *psolve {
        /* call the user-supplied psolve routine, and accumulate count */
        let IDAMem { ida_tn, ida_cj, ida_user_data, .. } = ida_mem;
        let retval = psolve(*ida_tn, ycur, ypcur, rcur, r, z, *ida_cj, tol, ida_user_data);
        *nps += 1;
        retval
    } else {
        0
    }
}

/*---------------------------------------------------------------
  idaLsDQJac:

  This routine is a wrapper for the Dense and Band
  implementations of the difference quotient Jacobian
  approximation routines.
  (The C Jac NULL check and N_Vector capability checks always
  hold for the workspace matrix enum and the serial vector.)
---------------------------------------------------------------*/
#[allow(clippy::too_many_arguments)]
pub fn idaLsDQJac(
    ida_mem: &mut IDAMem,
    nreDQ: &mut i64,
    tt: f64,
    c_j: f64,
    yy: &NVector,
    yp: &NVector,
    rr: &NVector,
    jac: &mut SUNMatrix,
) -> i32 {
    /* Call the matrix-structure-specific DQ approximation routine */
    match jac {
        SUNMatrix::Dense(dm) => idaLsDenseDQJac(ida_mem, nreDQ, tt, c_j, yy, yp, rr, dm),
        SUNMatrix::Band(bm) => idaLsBandDQJac(ida_mem, nreDQ, tt, c_j, yy, yp, rr, bm),
        _ => {
            IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "idaLsDQJac", file!(),
                            "unrecognized matrix type for idaLsDQJac");
            IDA_ILL_INPUT
        }
    }
}

/*---------------------------------------------------------------
  idaLsDenseDQJac

  This routine generates a dense difference quotient approximation
  to the Jacobian F_y + c_j*F_y'.  It assumes a dense SUNMatrix
  input (stored column-wise, and that elements within each column
  are contiguous).  The C jthCol N_VSetArrayPointer aliasing of
  the matrix column becomes a direct write into the column slice.

  The C code perturbs the caller's yy/yp data in place and restores
  it after each column; the pinned Rust hook signature passes yy/yp
  behind &, so the perturbations act on owned copies (the values
  seen by res are identical, and the C restore leaves the originals
  equal to their entry values anyway).  rtemp (C tmp1) = ida_tempv1.
---------------------------------------------------------------*/
#[allow(clippy::too_many_arguments)]
pub fn idaLsDenseDQJac(
    ida_mem: &mut IDAMem,
    nreDQ: &mut i64,
    tt: f64,
    c_j: f64,
    yy: &NVector,
    yp: &NVector,
    rr: &NVector,
    jac: &mut crate::sunmatrix_dense::DenseMatrix,
) -> i32 {
    /* access matrix dimension */
    let n = jac.n;

    /* perturbable copies of yy / yp (see the header comment) */
    let mut ytemp = yy.clone();
    let mut yptemp = yp.clone();

    let res = ida_mem.ida_res.unwrap();
    let srur = SUNRsqrt(ida_mem.ida_uround);
    let mut retval = 0;

    for j in 0..(n as usize) {
        /* Generate the jth col of J(tt,yy,yp) as delta(F)/delta(y_j). */

        /* Save y_j and yp_j values. */
        let yj = ytemp.data[j];
        let ypj = yptemp.data[j];

        /* Set increment inc to y_j based on sqrt(uround)*abs(y_j), with
        adjustments using yp_j and ewt_j if this is small, and a further
        adjustment to give it the same sign as hh*yp_j. */

        let mut inc = SUNMAX(
            srur * SUNMAX(SUNRabs(yj), SUNRabs(ida_mem.ida_hh * ypj)),
            ONE / ida_mem.ida_ewt.data[j],
        );

        if ida_mem.ida_hh * ypj < ZERO {
            inc = -inc;
        }
        inc = (yj + inc) - yj;

        /* Adjust sign(inc) again if y_j has an inequality constraint. */
        if ida_mem.ida_constraintsSet {
            let conj = ida_mem.ida_constraints.data[j];
            if SUNRabs(conj) == ONE {
                if (yj + inc) * conj < ZERO {
                    inc = -inc;
                }
            } else if SUNRabs(conj) == TWO && (yj + inc) * conj <= ZERO {
                inc = -inc;
            }
        }

        /* Increment y_j and yp_j, call res, and break on error return. */
        ytemp.data[j] += inc;
        yptemp.data[j] += c_j * inc;

        retval = {
            let IDAMem { ida_tempv1, ida_user_data, .. } = ida_mem;
            res(tt, &ytemp, &yptemp, ida_tempv1, ida_user_data)
        };
        *nreDQ += 1;
        if retval != 0 {
            break;
        }

        /* Construct difference quotient in jthCol:
           N_VLinearSum(inc_inv, rtemp, -inc_inv, rr, jthCol) */
        let inc_inv = ONE / inc;
        let col_j = jac.col_mut(j as i64);
        for (i, cij) in col_j.iter_mut().enumerate() {
            *cij = inc_inv * (ida_mem.ida_tempv1.data[i] - rr.data[i]);
        }

        /*  reset y_j, yp_j */
        ytemp.data[j] = yj;
        yptemp.data[j] = ypj;
    }

    retval
}

/*---------------------------------------------------------------
  idaLsBandDQJac

  This routine generates a banded difference quotient approximation
  JJ to the DAE system Jacobian J.  It assumes a band SUNMatrix
  input (stored column-wise, and that elements within each column
  are contiguous).  The columns of the Jacobian are constructed
  using mupper + mlower + 1 calls to the res routine, and
  appropriate differencing.  rtemp (C tmp1) = ida_tempv1,
  ytemp (tmp2) = ida_tempv2, yptemp (tmp3) = ida_tempv3.
  The return value is either 0, or the nonzero value returned by
  the res routine, if any.
  ---------------------------------------------------------------*/
#[allow(clippy::too_many_arguments)]
pub fn idaLsBandDQJac(
    ida_mem: &mut IDAMem,
    nreDQ: &mut i64,
    tt: f64,
    c_j: f64,
    yy: &NVector,
    yp: &NVector,
    rr: &NVector,
    jac: &mut crate::sunmatrix_band::BandMatrix,
) -> i32 {
    /* access matrix dimensions */
    let n = jac.n;
    let mupper = jac.mu;
    let mlower = jac.ml;

    /* Initialize ytemp and yptemp. */
    {
        let IDAMem { ida_tempv2, ida_tempv3, .. } = ida_mem;
        N_VScale(ONE, yy, ida_tempv2);
        N_VScale(ONE, yp, ida_tempv3);
    }

    /* Compute miscellaneous values for the Jacobian computation. */
    let srur = SUNRsqrt(ida_mem.ida_uround);
    let width = mlower + mupper + 1;
    let ngroups = width.min(n);

    let res = ida_mem.ida_res.unwrap();
    let mut retval = 0;

    /* Loop over column groups. */
    for group in 1..=ngroups {
        /* Increment all yy[j] and yp[j] for j in this group. */
        let mut j = group - 1;
        while j < n {
            let ju = j as usize;
            let yj = yy.data[ju];
            let ypj = yp.data[ju];
            let ewtj = ida_mem.ida_ewt.data[ju];

            /* Set increment inc to yj based on sqrt(uround)*abs(yj), with
              adjustments using ypj and ewtj if this is small, and a further
              adjustment to give it the same sign as hh*ypj. */
            let mut inc = SUNMAX(
                srur * SUNMAX(SUNRabs(yj), SUNRabs(ida_mem.ida_hh * ypj)),
                ONE / ewtj,
            );
            if ida_mem.ida_hh * ypj < ZERO {
                inc = -inc;
            }
            inc = (yj + inc) - yj;

            /* Adjust sign(inc) again if yj has an inequality constraint. */
            if ida_mem.ida_constraintsSet {
                let conj = ida_mem.ida_constraints.data[ju];
                if SUNRabs(conj) == ONE {
                    if (yj + inc) * conj < ZERO {
                        inc = -inc;
                    }
                } else if SUNRabs(conj) == TWO && (yj + inc) * conj <= ZERO {
                    inc = -inc;
                }
            }

            /* Increment yj and ypj. */
            ida_mem.ida_tempv2.data[ju] += inc;
            ida_mem.ida_tempv3.data[ju] += c_j * inc;

            j += width;
        }

        /* Call res routine with incremented arguments. */
        retval = {
            let IDAMem { ida_tempv1, ida_tempv2, ida_tempv3, ida_user_data, .. } = ida_mem;
            res(tt, ida_tempv2, ida_tempv3, ida_tempv1, ida_user_data)
        };
        *nreDQ += 1;
        if retval != 0 {
            break;
        }

        /* Loop over the indices j in this group again. */
        let mut j = group - 1;
        while j < n {
            let ju = j as usize;

            /* Reset ytemp and yptemp components that were perturbed. */
            let yj = yy.data[ju];
            ida_mem.ida_tempv2.data[ju] = yj;
            let ypj = yp.data[ju];
            ida_mem.ida_tempv3.data[ju] = ypj;
            let ewtj = ida_mem.ida_ewt.data[ju];

            /* Set increment inc exactly as above. */
            let mut inc = SUNMAX(
                srur * SUNMAX(SUNRabs(yj), SUNRabs(ida_mem.ida_hh * ypj)),
                ONE / ewtj,
            );
            if ida_mem.ida_hh * ypj < ZERO {
                inc = -inc;
            }
            inc = (yj + inc) - yj;
            if ida_mem.ida_constraintsSet {
                let conj = ida_mem.ida_constraints.data[ju];
                if SUNRabs(conj) == ONE {
                    if (yj + inc) * conj < ZERO {
                        inc = -inc;
                    }
                } else if SUNRabs(conj) == TWO && (yj + inc) * conj <= ZERO {
                    inc = -inc;
                }
            }

            /* Load the difference quotient Jacobian elements for column j */
            let inc_inv = ONE / inc;
            let i1 = 0.max(j - mupper);
            let i2 = (j + mlower).min(n - 1);
            for i in i1..=i2 {
                let val =
                    inc_inv * (ida_mem.ida_tempv1.data[i as usize] - rr.data[i as usize]);
                jac.set(i, j, val);
            }

            j += width;
        }
    }

    retval
}

/*---------------------------------------------------------------
  idaLsDQJtimes

  This routine generates a difference quotient approximation to
  the matrix-vector product z = Jv, where J is the system
  Jacobian. The approximation is
       Jv = [F(t,y1,yp1) - F(t,y,yp)]/sigma,
  where
       y1 = y + sigma*v,  yp1 = yp + cj*sigma*v,
       sigma = sqrt(Neq)*dqincfac.
  The return value from the call to res is saved in order to set
  the return flag from idaLsSolve.

  (C signature: idaLsDQJtimes(tt, yy, yp, rr, v, Jv, c_j, ida_mem,
  work1, work2).  The bits the C body pulls out of IDALsMem —
  gmres = SUNLinSolGetID(LS) in {SPGMR, SPFGMR}, nrmfac, dqincfac,
  jt_res, the nreDQ counter — arrive as explicit parameters because
  IDALsMem is destructured at the call sites; `ewt` is
  IDA_mem->ida_ewt, passed separately so the iterative-solve closure
  can supply the caller's weight vector, which is that same vector
  at every C call site.)
  ---------------------------------------------------------------*/
#[allow(clippy::too_many_arguments)]
pub fn idaLsDQJtimes(
    ida_mem: &mut IDAMem,
    ewt: &NVector,
    nreDQ: &mut i64,
    gmres: bool,
    nrmfac: f64,
    dqincfac: f64,
    jt_res: Option<IDAResFn>,
    tt: f64,
    yy: &NVector,
    yp: &NVector,
    rr: &NVector,
    v: &NVector,
    Jv: &mut NVector,
    c_j: f64,
    work1: &mut NVector,
    work2: &mut NVector,
) -> i32 {
    let mut sig = if gmres {
        nrmfac * dqincfac
    } else {
        dqincfac / N_VWrmsNorm(v, ewt)
    };

    /* Rename work1 and work2 for readability */
    let y_tmp = work1;
    let yp_tmp = work2;

    let res = jt_res.unwrap();
    let mut retval = 0;

    for _iter in 0..MAX_ITERS {
        /* Set y_tmp = yy + sig*v, yp_tmp = yp + cj*sig*v. */
        N_VLinearSum(sig, v, ONE, yy, y_tmp);
        N_VLinearSum(c_j * sig, v, ONE, yp, yp_tmp);

        /* Call res for Jv = F(t, y_tmp, yp_tmp), and return if it failed. */
        retval = res(tt, y_tmp, yp_tmp, Jv, &mut ida_mem.ida_user_data);
        *nreDQ += 1;
        if retval == 0 {
            break;
        }
        if retval < 0 {
            return -1;
        }

        sig *= PT25;
    }

    if retval > 0 {
        return 1;
    }

    /* Set Jv to [Jv - rr]/sig and return. */
    let siginv = ONE / sig;
    Jv.linear_sum_with(siginv, -siginv, rr);

    0
}

/*---------------------------------------------------------------
 idaLsInitialize

 This routine performs remaining initializations specific
 to the iterative linear solver interface (and solver itself)

 C: int idaLsInitialize(IDAMem IDA_mem); the IDALsMem is threaded
 explicitly per the workspace take() convention.
---------------------------------------------------------------*/
pub fn idaLsInitialize(ida_mem: &mut IDAMem, idals_mem: &mut IDALsMem) -> i32 {
    /* Test for valid combinations of matrix & Jacobian routines: */
    if idals_mem.J.is_none() {
        /* If SUNMatrix A is NULL: ensure 'jac' function pointer is NULL */
        idals_mem.jacDQ = SUNFALSE;
        idals_mem.jac = None;
    } else if idals_mem.jacDQ {
        /* If J is non-NULL, and 'jac' is not user-supplied:
           - if J is dense or band, ensure that our DQ approx. is used
           - otherwise => error */
        let ok = matches!(
            idals_mem.J,
            Some(SUNMatrix::Dense(_)) | Some(SUNMatrix::Band(_))
        );
        if !ok {
            IDAProcessError(Some(ida_mem), IDALS_ILL_INPUT, line!(), "idaLsInitialize", file!(),
                            "No Jacobian constructor available for SUNMatrix type");
            idals_mem.last_flag = IDALS_ILL_INPUT;
            return IDALS_ILL_INPUT;
        }
        idals_mem.jac = None; /* internal idaLsDQJac */
    } else {
        /* If J is non-NULL, and 'jac' is user-supplied: the C J_data
           reset is a no-op here (user_data always comes from IDAMem) */
    }

    /* reset counters */
    idaLsInitializeCounters(idals_mem);

    /* Set Jacobian-related fields, based on jtimesDQ (the C jt_data
       assignments are no-ops here) */
    if idals_mem.jtimesDQ {
        idals_mem.jtsetup = None;
        idals_mem.jtimes = None; /* internal idaLsDQJtimes */
    }

    /* if J is NULL and psetup is not present, then idaLsSetup does
       not need to be called, so disable the lsetup dispatch
       (C idas_ls.c sets IDA_mem->ida_lsetup = NULL; setup_disabled
       carries that state, see idas_ls_impl.rs.  In C an internal
       preconditioner module installs a real pset fn pointer, so the
       None|User pattern below — "no internal prec module attached" —
       keeps this line correct when the idas_bbdpre PrecModule::BBDPre
       variant lands, per the donor's !BBDPre clause). */
    idals_mem.setup_disabled = idals_mem.J.is_none()
        && idals_mem.pset.is_none()
        && matches!(idals_mem.prec_module, PrecModule::None | PrecModule::User);

    /* When using a matrix-embedded linear solver disable lsetup call */
    if idals_mem.LS.ls_type() == SUNLINEARSOLVER_MATRIX_EMBEDDED {
        idals_mem.setup_disabled = SUNTRUE;
        idals_mem.scalesol = SUNFALSE;
    }

    /* Call LS initialize routine */
    idals_mem.last_flag = idals_mem.LS.initialize();
    idals_mem.last_flag
}

/* idaLsInit — the interface hook dispatched by the main solver
   (idas.rs LsModule dispatch name for idaLsInitialize; C ida_linit). */
pub fn idaLsInit(ida_mem: &mut IDAMem, idals_mem: &mut IDALsMem) -> i32 {
    idaLsInitialize(ida_mem, idals_mem)
}

/*---------------------------------------------------------------
 idaLsSetup

 This calls the Jacobian evaluation routine (if using a SUNMatrix
 object), updates counters, and calls the LS 'setup' routine to
 prepare for subsequent calls to the LS 'solve' routine.

 C: int idaLsSetup(IDAMem IDA_mem, N_Vector y, N_Vector yp,
                   N_Vector r, N_Vector vt1, N_Vector vt2,
                   N_Vector vt3); the three tmp vectors are always
 ida_tempv1/2/3 (both idaNlsLSetup and IDACalcIC pass them), taken
 from the IDAMem fields here.  Callers pass y/yp/r detached from
 IDAMem (std::mem::take) since in C they alias IDAMem fields
 (yy/yp/savres resp. yy0/yp0/delta in IDACalcIC).
---------------------------------------------------------------*/
pub fn idaLsSetup(
    ida_mem: &mut IDAMem,
    idals_mem: &mut IDALsMem,
    y: &NVector,
    yp: &NVector,
    r: &NVector,
) -> i32 {
    /* Immediately return when using matrix-embedded linear solver */
    if idals_mem.LS.ls_type() == SUNLINEARSOLVER_MATRIX_EMBEDDED {
        idals_mem.last_flag = IDALS_SUCCESS;
        return idals_mem.last_flag;
    }

    /* (C parks the IDALs ycur/ypcur/rcur pointers here for the ATimes /
       PSetup/PSolve bridges; dropped — passed along as arguments.) */

    /* Update values for last jac/pset call */
    idals_mem.nstlj = ida_mem.ida_nst;
    idals_mem.tnlj = ida_mem.ida_tn;

    /* recompute if J if it is non-NULL */
    if idals_mem.J.is_some() {
        /* Increment nje counter. */
        idals_mem.nje += 1;

        /* Clear the linear system matrix if necessary */
        if idals_mem.LS.ls_type() == SUNLINEARSOLVER_DIRECT {
            let retval = idals_mem.J.as_mut().unwrap().zero();
            if retval != 0 {
                IDAProcessError(Some(ida_mem), IDALS_SUNMAT_FAIL, line!(), "idaLsSetup", file!(),
                                MSG_LS_MATZERO_FAILED);
                idals_mem.last_flag = IDALS_SUNMAT_FAIL;
                return idals_mem.last_flag;
            }
        }

        /* Call Jacobian routine */
        let retval = if idals_mem.jacDQ {
            let tt = ida_mem.ida_tn;
            let cj = ida_mem.ida_cj;
            let j = idals_mem.J.as_mut().unwrap();
            idaLsDQJac(ida_mem, &mut idals_mem.nreDQ, tt, cj, y, yp, r, j)
        } else {
            let jac = idals_mem.jac.unwrap();
            let j = idals_mem.J.as_mut().unwrap();
            let IDAMem {
                ida_tn,
                ida_cj,
                ida_user_data,
                ida_tempv1,
                ida_tempv2,
                ida_tempv3,
                ..
            } = ida_mem;
            jac(*ida_tn, *ida_cj, y, yp, r, j, ida_user_data, ida_tempv1, ida_tempv2, ida_tempv3)
        };
        if retval < 0 {
            IDAProcessError(Some(ida_mem), IDALS_JACFUNC_UNRECVR, line!(), "idaLsSetup", file!(),
                            MSG_LS_JACFUNC_FAILED);
            idals_mem.last_flag = IDALS_JACFUNC_UNRECVR;
            return -1;
        }
        if retval > 0 {
            idals_mem.last_flag = IDALS_JACFUNC_RECVR;
            return 1;
        }
    }

    /* Call LS setup routine -- the LS will call idaLsPSetup if applicable:
       for direct solvers this factors J; for the iterative solvers the
       generic C SUNLinSolSetup calls the registered idaLsPSetup (when the
       user pset is non-NULL, or an internal prec module — idas_bbdpre —
       installed one; the !(None|User) pattern covers the future BBDPre
       variant) and maps a failure to SUNLS_PSET_FAIL_UNREC /
       SUNLS_PSET_FAIL_REC (sunlinsol_spgmr.c et al.), reproduced here */
    let direct_setup = matches!(idals_mem.LS, LinearSolver::Dense(_) | LinearSolver::Band(_));
    idals_mem.last_flag = if direct_setup {
        let IDALsMem { LS, J, .. } = idals_mem;
        LS.setup(J.as_mut())
    } else if idals_mem.pset.is_some()
        || !matches!(idals_mem.prec_module, PrecModule::None | PrecModule::User)
    {
        let retval = idaLsPSetup(ida_mem, idals_mem, y, yp, r);
        if retval != 0 {
            if retval < 0 { SUNLS_PSET_FAIL_UNREC } else { SUNLS_PSET_FAIL_REC }
        } else {
            SUN_SUCCESS
        }
    } else {
        SUN_SUCCESS
    };
    idals_mem.last_flag
}

/*---------------------------------------------------------------
 idaLsSolve

 This routine interfaces between IDA and the generic
 SUNLinearSolver object LS, by setting the appropriate tolerance
 and scaling vectors, calling the solver, accumulating
 statistics from the solve for use/reporting by IDA, and scaling
 the result if using a non-NULL SUNMatrix and cjratio does not
 equal one.

 C: int idaLsSolve(IDAMem IDA_mem, N_Vector b, N_Vector weight,
                   N_Vector ycur, N_Vector ypcur, N_Vector rescur);
 callers pass weight/ycur/ypcur/rescur detached from IDAMem
 (std::mem::take) since in C they alias IDAMem fields
 (ewt/yy/yp/savres resp. the IDACalcIC vectors).
---------------------------------------------------------------*/
pub fn idaLsSolve(
    ida_mem: &mut IDAMem,
    idals_mem: &mut IDALsMem,
    b: &mut NVector,
    weight: &NVector,
    ycur: &NVector,
    ypcur: &NVector,
    rescur: &NVector,
) -> i32 {
    /* If the linear solver is iterative: set convergence test constant tol,
       in terms of the Newton convergence test constant epsNewt and safety
       factors. The factor nrmfac assures that the convergence test is
       applied to the WRMS norm of the residual vector, rather than the
       weighted L2 norm. */
    let tol = if idals_mem.iterative {
        idals_mem.nrmfac * idals_mem.eplifac * ida_mem.ida_epsNewt
    } else {
        ZERO
    };

    /* (C parks the ycur/ypcur/rcur pointers here for use by the Atimes
       and Psolve interface routines; dropped — the iterative-solve
       closures below capture the arguments directly.) */

    /* Set scaling vectors for LS to use (if applicable): every workspace
       iterative solver accepts the weight/weight scaling vectors
       (s1 = s2 below), matching the C solvers that implement
       SUNLinSolSetScalingVectors — so the w_mean tolerance-adjustment
       branch for iterative solvers without scaling support is never
       taken. */

    /* Set initial guess x = 0 to LS */
    N_VConst(ZERO, &mut idals_mem.x);

    /* Set zero initial guess flag (infallible here) */
    idals_mem.LS.set_zero_guess(SUNTRUE);

    /* (nps_inc snapshot dropped: it only feeds the SUNLogInfo output,
       which the default build compiles out.) */

    /* If a user-provided jtsetup routine is supplied, call that here */
    if let Some(jtsetup) = idals_mem.jtsetup {
        idals_mem.last_flag = {
            let IDAMem { ida_tn, ida_cj, ida_user_data, .. } = ida_mem;
            jtsetup(*ida_tn, ycur, ypcur, rescur, *ida_cj, ida_user_data)
        };
        idals_mem.njtsetup += 1;
        if idals_mem.last_flag != 0 {
            IDAProcessError(Some(ida_mem), idals_mem.last_flag, line!(), "idaLsSolve", file!(),
                            MSG_LS_JTSETUP_FAILED);
            return idals_mem.last_flag;
        }
    }

    /* Call solver */
    let retval = if matches!(idals_mem.LS, LinearSolver::Custom(_)) {
        /* matrix-embedded solver: gets (t, cj, user_data); IDA's cj rides
           in the CustomLinSol `gamma` slot (the trait is shared with the
           CVODE-family integrators) */
        let IDALsMem { LS, x, .. } = idals_mem;
        if let LinearSolver::Custom(cls) = LS {
            let IDAMem { ida_user_data, ida_tn, ida_cj, .. } = ida_mem;
            cls.solve(x, b, tol, *ida_tn, *ida_cj, ida_user_data)
        } else {
            unreachable!()
        }
    } else if !idals_mem.iterative {
        let IDALsMem { LS, J, x, .. } = idals_mem;
        match (LS, J.as_mut()) {
            (LinearSolver::Dense(dls), Some(SUNMatrix::Dense(am))) => dls.solve(am, x, b),
            (LinearSolver::Band(bls), Some(SUNMatrix::Band(am))) => bls.solve(am, x, b),
            _ => SUN_ERR_ARG_INCOMPATIBLE,
        }
    } else {
        idaLsSolveIterative(ida_mem, idals_mem, b, weight, ycur, ypcur, rescur, tol)
    };

    /* Copy appropriate result to b (depending on solver type) */
    if idals_mem.iterative {
        /* Retrieve solver statistics (resnorm dropped: logging only) */
        let nli_inc = idals_mem.LS.num_iters();

        /* Copy appropriate result to b (C idas_ls.c idaLsSolve): when the
           solve converged in 0 iterations (and the LS is not matrix-
           embedded) the correction is the LS residual vector — the
           preconditioned residual of the zero initial guess,
           N_VScale(ONE, SUNLinSolResid(LS), b) — NOT x (which is still 0).
           Returning x=0 gives a spurious zero Newton correction and derails
           the step trajectory whenever the preconditioner alone meets the
           linear tolerance (e.g. idaHeat2D_kry). Otherwise copy x. */
        if nli_inc == 0 && !matches!(idals_mem.LS, LinearSolver::Custom(_)) {
            b.data.copy_from_slice(&idals_mem.LS.resid().data);
        } else {
            b.data.copy_from_slice(&idals_mem.x.data);
        }

        /* Increment nli counter */
        idals_mem.nli += nli_inc as i64;
    } else {
        /* Copy x to b */
        b.data.copy_from_slice(&idals_mem.x.data);
    }

    /* If using a direct or matrix-iterative solver, scale the correction to
       account for change in cj */
    if idals_mem.scalesol && ida_mem.ida_cjratio != ONE {
        b.scale_inplace(TWO / (ONE + ida_mem.ida_cjratio));
    }

    /* Increment ncfl counter */
    if retval != SUN_SUCCESS {
        idals_mem.ncfl += 1;
    }

    /* Interpret solver return value */
    idals_mem.last_flag = retval;

    match retval {
        SUN_SUCCESS => 0,
        SUNLS_RES_REDUCED | SUNLS_CONV_FAIL | SUNLS_PSOLVE_FAIL_REC
        | SUNLS_PACKAGE_FAIL_REC | SUNLS_QRFACT_FAIL | SUNLS_LUFACT_FAIL => 1,
        SUN_ERR_ARG_CORRUPT | SUN_ERR_ARG_INCOMPATIBLE | SUN_ERR_MEM_FAIL | SUNLS_GS_FAIL
        | SUNLS_QRSOL_FAIL => -1,
        SUN_ERR_EXT_FAIL => {
            IDAProcessError(Some(ida_mem), SUN_ERR_EXT_FAIL, line!(), "idaLsSolve", file!(),
                            "Failure in SUNLinSol external package");
            -1
        }
        SUNLS_PSOLVE_FAIL_UNREC => {
            IDAProcessError(Some(ida_mem), SUNLS_PSOLVE_FAIL_UNREC, line!(), "idaLsSolve",
                            file!(), MSG_LS_PSOLVE_FAILED);
            -1
        }
        /* the C switch has no default and falls through to return(0) —
           notably the SUNLS_ATIMES_FAIL_* flags land here */
        _ => 0,
    }
}

/* Iterative Krylov solve: builds the ATimes (idaLsATimes) and PSolve
   (idaLsPSolve) callbacks over the integrator memory.  The two closures
   share ida_mem through a RefCell — the solvers never call them
   re-entrantly (checked at runtime).  ycur/ypcur/rescur are the
   arguments the C code parks in IDALsMem for these bridges; `weight`
   doubles as the (read-only) scaling vectors s1 = s2 and as the ewt
   used by the DQ Jv sigma (in C the DQ routine reads
   IDA_mem->ida_ewt, which is that same vector at every call site —
   idaNlsLSolve and IDACalcIC both pass ida_ewt as the weight). */
#[allow(clippy::too_many_arguments)]
fn idaLsSolveIterative(
    ida_mem: &mut IDAMem,
    idals_mem: &mut IDALsMem,
    b: &NVector,
    weight: &NVector,
    ycur: &NVector,
    ypcur: &NVector,
    rescur: &NVector,
    tol: f64,
) -> i32 {
    let gmres = matches!(idals_mem.LS, LinearSolver::Spgmr(_) | LinearSolver::Spfgmr(_));
    let IDALsMem {
        LS,
        x,
        ytemp,
        yptemp,
        njtimes,
        nreDQ,
        nps,
        jtimes,
        jtimesDQ,
        jt_res,
        psolve,
        nrmfac,
        dqincfac,
        prec_module,
        ..
    } = idals_mem;
    let jtimes = *jtimes;
    let jtimes_dq = *jtimesDQ;
    let jt_res = *jt_res;
    let psolve_fn = *psolve;
    let nrmfac = *nrmfac;
    let dqincfac = *dqincfac;

    /* The internal idas_bbdpre module dispatches its psolve through the
       PrecModule payload rather than a psolve fn pointer. */
    let mut bbd_pdata = match prec_module {
        PrecModule::BBDPre(p) => Some(&mut **p),
        _ => None,
    };

    /* In C, PSolve is registered with the LS only when a non-NULL
       psolve was supplied (IDASetPreconditioner), or the idas_bbdpre
       module installed IDABBDPrecSolve. */
    let has_psolve = psolve_fn.is_some() || bbd_pdata.is_some();

    let idam = RefCell::new(&mut *ida_mem);

    /* idaLsATimes body (see the standalone bridge above) */
    let mut atimes = |v: &NVector, z: &mut NVector| -> i32 {
        let mut guard = idam.borrow_mut();
        let imr: &mut IDAMem = &mut *guard;
        let jret = if jtimes_dq {
            let tn = imr.ida_tn;
            let cj = imr.ida_cj;
            idaLsDQJtimes(imr, weight, nreDQ, gmres, nrmfac, dqincfac, jt_res, tn, ycur,
                          ypcur, rescur, v, z, cj, ytemp, yptemp)
        } else {
            let jt = jtimes.unwrap();
            let IDAMem { ida_tn, ida_cj, ida_user_data, .. } = imr;
            jt(*ida_tn, ycur, ypcur, rescur, v, z, *ida_cj, ida_user_data, ytemp, yptemp)
        };
        *njtimes += 1;
        jret
    };

    /* idaLsPSolve body */
    let mut psolve_cb = |r: &NVector, z: &mut NVector, ptol: f64, _lr: i32| -> i32 {
        let mut guard = idam.borrow_mut();
        let imr: &mut IDAMem = &mut *guard;
        let ret = if let Some(ps) = psolve_fn {
            let IDAMem { ida_tn, ida_cj, ida_user_data, .. } = imr;
            ps(*ida_tn, ycur, ypcur, rescur, r, z, *ida_cj, ptol, ida_user_data)
        } else if let Some(pdata) = bbd_pdata.as_deref_mut() {
            /* internal idas_bbdpre module (C psolve = IDABBDPrecSolve) */
            IDABBDPrecSolve(pdata, r, z)
        } else {
            0
        };
        *nps += 1;
        ret
    };

    /* weight/weight scaling, as in the C
       SUNLinSolSetScalingVectors(LS, weight, weight) */
    if has_psolve {
        LS.solve(None, x, b, tol, &mut atimes, Some(&mut psolve_cb), Some(weight), Some(weight))
    } else {
        LS.solve(None, x, b, tol, &mut atimes, None, Some(weight), Some(weight))
    }
}

/*---------------------------------------------------------------
 idaLsPerf: accumulates performance statistics information
 for IDA

 C: int idaLsPerf(IDAMem IDA_mem, int perftask); the int return
 (0/1) is discarded by every C call site, so the pinned hook
 signature drops it.  Dispatched only for iterative
 SUNLinearSolver objects: the C installs ida_lperf only when
 `iterative`, and idas.rs guards the dispatch on
 idals_mem.iterative (idas_impl.rs contract).
---------------------------------------------------------------*/
pub fn idaLsPerf(ida_mem: &mut IDAMem, idals_mem: &mut IDALsMem, perftask: i32) {
    /* when perftask == 0, store current performance statistics */
    if perftask == 0 {
        idals_mem.nst0 = ida_mem.ida_nst;
        idals_mem.nni0 = ida_mem.ida_nni;
        idals_mem.ncfn0 = ida_mem.ida_ncfn;
        idals_mem.ncfl0 = idals_mem.ncfl;
        idals_mem.nwarn = 0;
        return;
    }

    /* Compute statistics since last call

       Note: the performance monitor that checked whether the average
         number of linear iterations was too close to maxl has been
         removed, since the 'maxl' value is no longer owned by the
         IDALs interface.
     */
    let nstd = ida_mem.ida_nst - idals_mem.nst0;
    let nnid = ida_mem.ida_nni - idals_mem.nni0;
    if nstd == 0 || nnid == 0 {
        return;
    }

    let rcfn = ((ida_mem.ida_ncfn - idals_mem.ncfn0) as f64) / (nstd as f64);
    let rcfl = ((idals_mem.ncfl - idals_mem.ncfl0) as f64) / (nnid as f64);
    let lcfn = rcfn > PT9;
    let lcfl = rcfl > PT9;
    if !(lcfn || lcfl) {
        return;
    }
    idals_mem.nwarn += 1;
    if idals_mem.nwarn > 10 {
        return;
    }
    if lcfn {
        IDAProcessError(Some(ida_mem), IDA_WARNING, line!(), "idaLsPerf", file!(),
            &format!("Warning: at t = {}, poor iterative algorithm performance. \
Nonlinear convergence failure rate is {}.",
                     fmt_g(ida_mem.ida_tn, 0, 15), fmt_g(rcfn, 0, 15)));
    }
    if lcfl {
        IDAProcessError(Some(ida_mem), IDA_WARNING, line!(), "idaLsPerf", file!(),
            &format!("Warning: at t = {}, poor iterative algorithm performance. \
Linear convergence failure rate is {}.",
                     fmt_g(ida_mem.ida_tn, 0, 15), fmt_g(rcfl, 0, 15)));
    }
}

/*---------------------------------------------------------------
 idaLsFree frees memory associated with the IDALs system
 solver interface (RAII: dropping the LsModule drops the IDALsMem,
 its LinearSolver, its SUNMatrix J, the ytemp/yptemp/x vectors and
 any preconditioner module — the C pfree hook)
---------------------------------------------------------------*/
pub fn idaLsFree(ida_mem: &mut IDAMem) -> i32 {
    /* Return immediately if IDA_mem->ida_lmem is NULL */
    if ida_mem.ida_lmem.is_none() {
        return IDALS_SUCCESS;
    }

    /* Free N_Vector memory, nullify the ycur/ypcur/rcur and SUNMatrix
       pointers, free preconditioner memory (pfree) and the IDALs
       interface structure — all via drop */
    ida_mem.ida_lmem = LsModule::None;

    IDALS_SUCCESS
}

/*---------------------------------------------------------------
 idaLsInitializeCounters resets all counters from an
 IDALsMem structure.
---------------------------------------------------------------*/
pub fn idaLsInitializeCounters(idals_mem: &mut IDALsMem) -> i32 {
    idals_mem.nje = 0;
    idals_mem.nreDQ = 0;
    idals_mem.npe = 0;
    idals_mem.nli = 0;
    idals_mem.nps = 0;
    idals_mem.ncfl = 0;
    idals_mem.njtsetup = 0;
    idals_mem.njtimes = 0;
    0
}

/*================================================================
  PART II - backward problems

  Modeling follows the pinned cvodes_ls.rs PART II design:
  - The C static wrappers are installed as the INNER (backward)
    problem's forward-interface callbacks; their void* data is the
    OUTER ida_mem, which arrives here through the inner problem's
    UserData (idaa.rs installs it, as the C idaa.c does) and is
    recovered by downcast (idaLs_AccessIDAMem).
  - IDAB_mem.ida_lmem is Option<Box<dyn Any>> holding an IDALsMemB
    (idaLsB_downcast); idaLs_AccessLMemB / idaLs_AccessLMemBCur
    return the Vec index of the backward problem.
  - The ia_yyTmp/ia_ypTmp (and yyS/ypS) workspaces are taken out of
    the IDAadjMem as owned locals around the user-callback call and
    restored afterwards (borrow discipline, as in cvodes_ls.rs).
  ================================================================*/

/*---------------------------------------------------------------
  IDASLS Exported functions -- Required
  ---------------------------------------------------------------*/

/* IDASetLinearSolverB specifies the iterative linear solver
   for backward integration */
pub fn IDASetLinearSolverB(
    ida_mem: &mut IDAMem,
    which: i32,
    LS: LinearSolver,
    A: Option<SUNMatrix>,
) -> i32 {
    /* Was ASA initialized? */
    if !ida_mem.ida_adjMallocDone {
        IDAProcessError(Some(ida_mem), IDALS_NO_ADJ, line!(), "IDASetLinearSolverB", file!(),
                        MSG_LS_NO_ADJ);
        return IDALS_NO_ADJ;
    }
    let idaadj_mem = ida_mem.ida_adj_mem.as_mut().unwrap();

    /* Check the value of which */
    if which >= idaadj_mem.ia_nbckpbs {
        IDAProcessError(None, IDALS_ILL_INPUT, line!(), "IDASetLinearSolverB", file!(),
                        MSG_LS_BAD_WHICH);
        return IDALS_ILL_INPUT;
    }

    /* Find the IDABMem entry in the linked list corresponding to 'which'. */
    let idx = idaadj_mem.IDAB_mem.iter().position(|b| b.ida_index == which).unwrap();
    let idaB_mem = &mut idaadj_mem.IDAB_mem[idx];

    /* Get memory for IDALsMemRecB, initialize the Jacobian and
       preconditioner function pointers to NULL (Default), free any
       existing system solver attached to IDAB (drop on overwrite) and
       attach lmemB data (the C lfreeB hook is Rust Drop / idaLsFreeB) */
    idaB_mem.ida_lmem = Some(Box::new(IDALsMemB::default()));

    /* set the linear solver for this backward problem */
    let retval = IDASetLinearSolver(&mut idaB_mem.IDA_mem, LS, A);
    if retval != IDALS_SUCCESS {
        idaB_mem.ida_lmem = None;
    }

    retval
}

/*---------------------------------------------------------------
  IDASLS Exported functions -- Optional input/output
  ---------------------------------------------------------------*/

pub fn IDASetJacFnB(ida_mem: &mut IDAMem, which: i32, jacB: Option<IDALsJacFnB>) -> i32 {
    /* access relevant memory structures */
    let idx = match idaLs_AccessLMemB(ida_mem, which, "IDASetJacFnB") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let idaB_mem = &mut ida_mem.ida_adj_mem.as_mut().unwrap().IDAB_mem[idx];

    /* set jacB function pointer */
    idaLsB_downcast(idaB_mem).jacB = jacB;

    /* call corresponding routine for IDAB_mem structure */
    let wrapper: Option<IDALsJacFn> =
        if jacB.is_some() { Some(idaLsJacBWrapper) } else { None };
    IDASetJacFn(&mut idaB_mem.IDA_mem, wrapper)
}

pub fn IDASetJacFnBS(ida_mem: &mut IDAMem, which: i32, jacBS: Option<IDALsJacFnBS>) -> i32 {
    /* access relevant memory structures */
    let idx = match idaLs_AccessLMemB(ida_mem, which, "IDASetJacFnBS") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let idaB_mem = &mut ida_mem.ida_adj_mem.as_mut().unwrap().IDAB_mem[idx];

    /* set jacBS function pointer */
    idaLsB_downcast(idaB_mem).jacBS = jacBS;

    /* call corresponding routine for IDAB_mem structure */
    let wrapper: Option<IDALsJacFn> =
        if jacBS.is_some() { Some(idaLsJacBSWrapper) } else { None };
    IDASetJacFn(&mut idaB_mem.IDA_mem, wrapper)
}

pub fn IDASetEpsLinB(ida_mem: &mut IDAMem, which: i32, eplifacB: f64) -> i32 {
    /* access relevant memory structures */
    let idx = match idaLs_AccessLMemB(ida_mem, which, "IDASetEpsLinB") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let idaB_mem = &mut ida_mem.ida_adj_mem.as_mut().unwrap().IDAB_mem[idx];

    /* call corresponding routine for IDAB_mem structure */
    IDASetEpsLin(&mut idaB_mem.IDA_mem, eplifacB)
}

pub fn IDASetLSNormFactorB(ida_mem: &mut IDAMem, which: i32, nrmfacB: f64) -> i32 {
    /* access relevant memory structures */
    let idx = match idaLs_AccessLMemB(ida_mem, which, "IDASetLSNormFactorB") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let idaB_mem = &mut ida_mem.ida_adj_mem.as_mut().unwrap().IDAB_mem[idx];

    /* call corresponding routine for IDAB_mem structure */
    IDASetLSNormFactor(&mut idaB_mem.IDA_mem, nrmfacB)
}

pub fn IDASetLinearSolutionScalingB(ida_mem: &mut IDAMem, which: i32, onoffB: bool) -> i32 {
    /* access relevant memory structures */
    let idx = match idaLs_AccessLMemB(ida_mem, which, "IDASetLinearSolutionScalingB") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let idaB_mem = &mut ida_mem.ida_adj_mem.as_mut().unwrap().IDAB_mem[idx];

    /* call corresponding routine for IDAB_mem structure */
    IDASetLinearSolutionScaling(&mut idaB_mem.IDA_mem, onoffB)
}

pub fn IDASetIncrementFactorB(ida_mem: &mut IDAMem, which: i32, dqincfacB: f64) -> i32 {
    /* access relevant memory structures */
    let idx = match idaLs_AccessLMemB(ida_mem, which, "IDASetIncrementFactorB") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let idaB_mem = &mut ida_mem.ida_adj_mem.as_mut().unwrap().IDAB_mem[idx];

    /* call corresponding routine for IDAB_mem structure */
    IDASetIncrementFactor(&mut idaB_mem.IDA_mem, dqincfacB)
}

pub fn IDASetPreconditionerB(
    ida_mem: &mut IDAMem,
    which: i32,
    psetupB: Option<IDALsPrecSetupFnB>,
    psolveB: Option<IDALsPrecSolveFnB>,
) -> i32 {
    /* access relevant memory structures */
    let idx = match idaLs_AccessLMemB(ida_mem, which, "IDASetPreconditionerB") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let idaB_mem = &mut ida_mem.ida_adj_mem.as_mut().unwrap().IDAB_mem[idx];

    /* Set preconditioners for the backward problem. */
    {
        let lmemB = idaLsB_downcast(idaB_mem);
        lmemB.psetB = psetupB;
        lmemB.psolveB = psolveB;
    }

    /* Call the corresponding "set" routine for the backward problem */
    let idals_psetup: Option<IDALsPrecSetupFn> =
        if psetupB.is_some() { Some(idaLsPrecSetupB) } else { None };
    let idals_psolve: Option<IDALsPrecSolveFn> =
        if psolveB.is_some() { Some(idaLsPrecSolveB) } else { None };
    IDASetPreconditioner(&mut idaB_mem.IDA_mem, idals_psetup, idals_psolve)
}

pub fn IDASetPreconditionerBS(
    ida_mem: &mut IDAMem,
    which: i32,
    psetupBS: Option<IDALsPrecSetupFnBS>,
    psolveBS: Option<IDALsPrecSolveFnBS>,
) -> i32 {
    /* access relevant memory structures */
    let idx = match idaLs_AccessLMemB(ida_mem, which, "IDASetPreconditionerBS") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let idaB_mem = &mut ida_mem.ida_adj_mem.as_mut().unwrap().IDAB_mem[idx];

    /* Set preconditioners for the backward problem. */
    {
        let lmemB = idaLsB_downcast(idaB_mem);
        lmemB.psetBS = psetupBS;
        lmemB.psolveBS = psolveBS;
    }

    /* Call the corresponding "set" routine for the backward problem */
    let idals_psetup: Option<IDALsPrecSetupFn> =
        if psetupBS.is_some() { Some(idaLsPrecSetupBS) } else { None };
    let idals_psolve: Option<IDALsPrecSolveFn> =
        if psolveBS.is_some() { Some(idaLsPrecSolveBS) } else { None };
    IDASetPreconditioner(&mut idaB_mem.IDA_mem, idals_psetup, idals_psolve)
}

pub fn IDASetJacTimesB(
    ida_mem: &mut IDAMem,
    which: i32,
    jtsetupB: Option<IDALsJacTimesSetupFnB>,
    jtimesB: Option<IDALsJacTimesVecFnB>,
) -> i32 {
    /* access relevant memory structures */
    let idx = match idaLs_AccessLMemB(ida_mem, which, "IDASetJacTimesB") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let idaB_mem = &mut ida_mem.ida_adj_mem.as_mut().unwrap().IDAB_mem[idx];

    /* Set jacobian routines for the backward problem. */
    {
        let lmemB = idaLsB_downcast(idaB_mem);
        lmemB.jtsetupB = jtsetupB;
        lmemB.jtimesB = jtimesB;
    }

    /* Call the corresponding "set" routine for the backward problem */
    let idals_jtsetup: Option<IDALsJacTimesSetupFn> =
        if jtsetupB.is_some() { Some(idaLsJacTimesSetupB) } else { None };
    let idals_jtimes: Option<IDALsJacTimesVecFn> =
        if jtimesB.is_some() { Some(idaLsJacTimesVecB) } else { None };
    IDASetJacTimes(&mut idaB_mem.IDA_mem, idals_jtsetup, idals_jtimes)
}

pub fn IDASetJacTimesBS(
    ida_mem: &mut IDAMem,
    which: i32,
    jtsetupBS: Option<IDALsJacTimesSetupFnBS>,
    jtimesBS: Option<IDALsJacTimesVecFnBS>,
) -> i32 {
    /* access relevant memory structures */
    let idx = match idaLs_AccessLMemB(ida_mem, which, "IDASetJacTimesBS") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let idaB_mem = &mut ida_mem.ida_adj_mem.as_mut().unwrap().IDAB_mem[idx];

    /* Set jacobian routines for the backward problem. */
    {
        let lmemB = idaLsB_downcast(idaB_mem);
        lmemB.jtsetupBS = jtsetupBS;
        lmemB.jtimesBS = jtimesBS;
    }

    /* Call the corresponding "set" routine for the backward problem */
    let idals_jtsetup: Option<IDALsJacTimesSetupFn> =
        if jtsetupBS.is_some() { Some(idaLsJacTimesSetupBS) } else { None };
    let idals_jtimes: Option<IDALsJacTimesVecFn> =
        if jtimesBS.is_some() { Some(idaLsJacTimesVecBS) } else { None };
    IDASetJacTimes(&mut idaB_mem.IDA_mem, idals_jtsetup, idals_jtimes)
}

pub fn IDASetJacTimesResFnB(
    ida_mem: &mut IDAMem,
    which: i32,
    jtimesResFn: Option<IDAResFn>,
) -> i32 {
    /* access relevant memory structures */
    let idx = match idaLs_AccessLMemB(ida_mem, which, "IDASetJacTimesResFnB") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let idaB_mem = &mut ida_mem.ida_adj_mem.as_mut().unwrap().IDAB_mem[idx];

    /* call corresponding routine for IDAB_mem structure */
    IDASetJacTimesResFn(&mut idaB_mem.IDA_mem, jtimesResFn)
}

/*-----------------------------------------------------------------
  IDASLS Private functions for backwards problems
  -----------------------------------------------------------------*/

/* idaLsJacBWrapper interfaces to the IDAJacFnB routine provided
   by the user. idaLsJacBWrapper is of type IDALsJacFn. */
#[allow(clippy::too_many_arguments)]
fn idaLsJacBWrapper(
    tt: f64,
    c_jB: f64,
    yyB: &NVector,
    ypB: &NVector,
    rrB: &NVector,
    JacB: &mut SUNMatrix,
    ida_mem: &mut UserData,
    tmp1B: &mut NVector,
    tmp2B: &mut NVector,
    tmp3B: &mut NVector,
) -> i32 {
    /* access relevant memory structures */
    let ida_mem = match idaLs_AccessIDAMem(ida_mem, "idaLsJacBWrapper") {
        Ok(m) => m,
        Err(e) => return e,
    };
    let which = match idaLs_AccessLMemBCur(ida_mem, "idaLsJacBWrapper") {
        Ok(w) => w,
        Err(e) => return e,
    };

    /* Forward solution from interpolation */
    let (mut yyTmp, mut ypTmp, no_interp) = {
        let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
        (std::mem::take(&mut adj.ia_yyTmp), std::mem::take(&mut adj.ia_ypTmp), adj.ia_noInterp)
    };
    if !no_interp {
        let mut noS: Vec<NVector> = Vec::new();
        let mut noSp: Vec<NVector> = Vec::new();
        let retval = idaLsGetY(ida_mem, tt, &mut yyTmp, &mut ypTmp, &mut noS, &mut noSp);
        if retval != IDA_SUCCESS {
            let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
            adj.ia_yyTmp = yyTmp;
            adj.ia_ypTmp = ypTmp;
            /* (C error target: IDAB_mem->IDA_mem) */
            IDAProcessError(Some(ida_mem), -1, line!(), "idaLsJacBWrapper", file!(),
                            MSG_LS_BAD_T);
            return -1;
        }
    }

    /* Call user's adjoint jacB routine */
    let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
    let idaB_mem = &mut adj.IDAB_mem[which];
    let jacB = idaLsB_downcast(idaB_mem).jacB.unwrap();
    let retval = jacB(tt, c_jB, &yyTmp, &ypTmp, yyB, ypB, rrB, JacB,
                      &mut idaB_mem.ida_user_data, tmp1B, tmp2B, tmp3B);
    adj.ia_yyTmp = yyTmp;
    adj.ia_ypTmp = ypTmp;
    retval
}

/* idaLsJacBSWrapper interfaces to the IDAJacFnBS routine provided
   by the user. idaLsJacBSWrapper is of type IDALsJacFn. */
#[allow(clippy::too_many_arguments)]
fn idaLsJacBSWrapper(
    tt: f64,
    c_jB: f64,
    yyB: &NVector,
    ypB: &NVector,
    rrB: &NVector,
    JacB: &mut SUNMatrix,
    ida_mem: &mut UserData,
    tmp1B: &mut NVector,
    tmp2B: &mut NVector,
    tmp3B: &mut NVector,
) -> i32 {
    /* access relevant memory structures */
    let ida_mem = match idaLs_AccessIDAMem(ida_mem, "idaLsJacBSWrapper") {
        Ok(m) => m,
        Err(e) => return e,
    };
    let which = match idaLs_AccessLMemBCur(ida_mem, "idaLsJacBSWrapper") {
        Ok(w) => w,
        Err(e) => return e,
    };

    /* Get forward solution from interpolation. */
    let (mut yyTmp, mut ypTmp, mut yySTmp, mut ypSTmp, no_interp, interp_sensi) = {
        let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
        (
            std::mem::take(&mut adj.ia_yyTmp),
            std::mem::take(&mut adj.ia_ypTmp),
            std::mem::take(&mut adj.ia_yySTmp),
            std::mem::take(&mut adj.ia_ypSTmp),
            adj.ia_noInterp,
            adj.ia_interpSensi,
        )
    };
    if !no_interp {
        let retval = if interp_sensi {
            idaLsGetY(ida_mem, tt, &mut yyTmp, &mut ypTmp, &mut yySTmp, &mut ypSTmp)
        } else {
            let mut noS: Vec<NVector> = Vec::new();
            let mut noSp: Vec<NVector> = Vec::new();
            idaLsGetY(ida_mem, tt, &mut yyTmp, &mut ypTmp, &mut noS, &mut noSp)
        };
        if retval != IDA_SUCCESS {
            let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
            adj.ia_yyTmp = yyTmp;
            adj.ia_ypTmp = ypTmp;
            adj.ia_yySTmp = yySTmp;
            adj.ia_ypSTmp = ypSTmp;
            IDAProcessError(Some(ida_mem), -1, line!(), "idaLsJacBSWrapper", file!(),
                            MSG_LS_BAD_T);
            return -1;
        }
    }

    /* Call user's adjoint jacBS routine */
    let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
    let idaB_mem = &mut adj.IDAB_mem[which];
    let jacBS = idaLsB_downcast(idaB_mem).jacBS.unwrap();
    let retval = jacBS(tt, c_jB, &yyTmp, &ypTmp, &yySTmp, &ypSTmp, yyB, ypB, rrB, JacB,
                       &mut idaB_mem.ida_user_data, tmp1B, tmp2B, tmp3B);
    adj.ia_yyTmp = yyTmp;
    adj.ia_ypTmp = ypTmp;
    adj.ia_yySTmp = yySTmp;
    adj.ia_ypSTmp = ypSTmp;
    retval
}

/* idaLsPrecSetupB interfaces to the IDALsPrecSetupFnB
   routine provided by the user (installed as a forward
   IDALsPrecSetupFn: the ewt/hh extension args are unused) */
#[allow(clippy::too_many_arguments)]
fn idaLsPrecSetupB(
    tt: f64,
    yyB: &NVector,
    ypB: &NVector,
    rrB: &NVector,
    c_jB: f64,
    _ewtB: &NVector,
    _hhB: f64,
    ida_mem: &mut UserData,
) -> i32 {
    /* access relevant memory structures */
    let ida_mem = match idaLs_AccessIDAMem(ida_mem, "idaLsPrecSetupB") {
        Ok(m) => m,
        Err(e) => return e,
    };
    let which = match idaLs_AccessLMemBCur(ida_mem, "idaLsPrecSetupB") {
        Ok(w) => w,
        Err(e) => return e,
    };

    /* Get forward solution from interpolation. */
    let (mut yyTmp, mut ypTmp, no_interp) = {
        let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
        (std::mem::take(&mut adj.ia_yyTmp), std::mem::take(&mut adj.ia_ypTmp), adj.ia_noInterp)
    };
    if !no_interp {
        let mut noS: Vec<NVector> = Vec::new();
        let mut noSp: Vec<NVector> = Vec::new();
        let retval = idaLsGetY(ida_mem, tt, &mut yyTmp, &mut ypTmp, &mut noS, &mut noSp);
        if retval != IDA_SUCCESS {
            let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
            adj.ia_yyTmp = yyTmp;
            adj.ia_ypTmp = ypTmp;
            IDAProcessError(Some(ida_mem), -1, line!(), "idaLsPrecSetupB", file!(),
                            MSG_LS_BAD_T);
            return -1;
        }
    }

    /* Call user's adjoint precondB routine */
    let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
    let idaB_mem = &mut adj.IDAB_mem[which];
    let psetB = idaLsB_downcast(idaB_mem).psetB.unwrap();
    let retval = psetB(tt, &yyTmp, &ypTmp, yyB, ypB, rrB, c_jB, &mut idaB_mem.ida_user_data);
    adj.ia_yyTmp = yyTmp;
    adj.ia_ypTmp = ypTmp;
    retval
}

/* idaLsPrecSetupBS interfaces to the IDALsPrecSetupFnBS routine
   provided by the user */
#[allow(clippy::too_many_arguments)]
fn idaLsPrecSetupBS(
    tt: f64,
    yyB: &NVector,
    ypB: &NVector,
    rrB: &NVector,
    c_jB: f64,
    _ewtB: &NVector,
    _hhB: f64,
    ida_mem: &mut UserData,
) -> i32 {
    /* access relevant memory structures */
    let ida_mem = match idaLs_AccessIDAMem(ida_mem, "idaLsPrecSetupBS") {
        Ok(m) => m,
        Err(e) => return e,
    };
    let which = match idaLs_AccessLMemBCur(ida_mem, "idaLsPrecSetupBS") {
        Ok(w) => w,
        Err(e) => return e,
    };

    /* Get forward solution from interpolation. */
    let (mut yyTmp, mut ypTmp, mut yySTmp, mut ypSTmp, no_interp, interp_sensi) = {
        let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
        (
            std::mem::take(&mut adj.ia_yyTmp),
            std::mem::take(&mut adj.ia_ypTmp),
            std::mem::take(&mut adj.ia_yySTmp),
            std::mem::take(&mut adj.ia_ypSTmp),
            adj.ia_noInterp,
            adj.ia_interpSensi,
        )
    };
    if !no_interp {
        let retval = if interp_sensi {
            idaLsGetY(ida_mem, tt, &mut yyTmp, &mut ypTmp, &mut yySTmp, &mut ypSTmp)
        } else {
            let mut noS: Vec<NVector> = Vec::new();
            let mut noSp: Vec<NVector> = Vec::new();
            idaLsGetY(ida_mem, tt, &mut yyTmp, &mut ypTmp, &mut noS, &mut noSp)
        };
        if retval != IDA_SUCCESS {
            let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
            adj.ia_yyTmp = yyTmp;
            adj.ia_ypTmp = ypTmp;
            adj.ia_yySTmp = yySTmp;
            adj.ia_ypSTmp = ypSTmp;
            IDAProcessError(Some(ida_mem), -1, line!(), "idaLsPrecSetupBS", file!(),
                            MSG_LS_BAD_T);
            return -1;
        }
    }

    /* Call user's adjoint precondBS routine */
    let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
    let idaB_mem = &mut adj.IDAB_mem[which];
    let psetBS = idaLsB_downcast(idaB_mem).psetBS.unwrap();
    let retval = psetBS(tt, &yyTmp, &ypTmp, &yySTmp, &ypSTmp, yyB, ypB, rrB, c_jB,
                        &mut idaB_mem.ida_user_data);
    adj.ia_yyTmp = yyTmp;
    adj.ia_ypTmp = ypTmp;
    adj.ia_yySTmp = yySTmp;
    adj.ia_ypSTmp = ypSTmp;
    retval
}

/* idaLsPrecSolveB interfaces to the IDALsPrecSolveFnB routine
   provided by the user */
#[allow(clippy::too_many_arguments)]
fn idaLsPrecSolveB(
    tt: f64,
    yyB: &NVector,
    ypB: &NVector,
    rrB: &NVector,
    rvecB: &NVector,
    zvecB: &mut NVector,
    c_jB: f64,
    deltaB: f64,
    ida_mem: &mut UserData,
) -> i32 {
    /* access relevant memory structures */
    let ida_mem = match idaLs_AccessIDAMem(ida_mem, "idaLsPrecSolveB") {
        Ok(m) => m,
        Err(e) => return e,
    };
    let which = match idaLs_AccessLMemBCur(ida_mem, "idaLsPrecSolveB") {
        Ok(w) => w,
        Err(e) => return e,
    };

    /* Get forward solution from interpolation. */
    let (mut yyTmp, mut ypTmp, no_interp) = {
        let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
        (std::mem::take(&mut adj.ia_yyTmp), std::mem::take(&mut adj.ia_ypTmp), adj.ia_noInterp)
    };
    if !no_interp {
        let mut noS: Vec<NVector> = Vec::new();
        let mut noSp: Vec<NVector> = Vec::new();
        let retval = idaLsGetY(ida_mem, tt, &mut yyTmp, &mut ypTmp, &mut noS, &mut noSp);
        if retval != IDA_SUCCESS {
            let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
            adj.ia_yyTmp = yyTmp;
            adj.ia_ypTmp = ypTmp;
            IDAProcessError(Some(ida_mem), -1, line!(), "idaLsPrecSolveB", file!(),
                            MSG_LS_BAD_T);
            return -1;
        }
    }

    /* Call user's adjoint psolveB routine */
    let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
    let idaB_mem = &mut adj.IDAB_mem[which];
    let psolveB = idaLsB_downcast(idaB_mem).psolveB.unwrap();
    let retval = psolveB(tt, &yyTmp, &ypTmp, yyB, ypB, rrB, rvecB, zvecB, c_jB, deltaB,
                         &mut idaB_mem.ida_user_data);
    adj.ia_yyTmp = yyTmp;
    adj.ia_ypTmp = ypTmp;
    retval
}

/* idaLsPrecSolveBS interfaces to the IDALsPrecSolveFnBS routine
   provided by the user */
#[allow(clippy::too_many_arguments)]
fn idaLsPrecSolveBS(
    tt: f64,
    yyB: &NVector,
    ypB: &NVector,
    rrB: &NVector,
    rvecB: &NVector,
    zvecB: &mut NVector,
    c_jB: f64,
    deltaB: f64,
    ida_mem: &mut UserData,
) -> i32 {
    /* access relevant memory structures */
    let ida_mem = match idaLs_AccessIDAMem(ida_mem, "idaLsPrecSolveBS") {
        Ok(m) => m,
        Err(e) => return e,
    };
    let which = match idaLs_AccessLMemBCur(ida_mem, "idaLsPrecSolveBS") {
        Ok(w) => w,
        Err(e) => return e,
    };

    /* Get forward solution from interpolation. */
    let (mut yyTmp, mut ypTmp, mut yySTmp, mut ypSTmp, no_interp, interp_sensi) = {
        let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
        (
            std::mem::take(&mut adj.ia_yyTmp),
            std::mem::take(&mut adj.ia_ypTmp),
            std::mem::take(&mut adj.ia_yySTmp),
            std::mem::take(&mut adj.ia_ypSTmp),
            adj.ia_noInterp,
            adj.ia_interpSensi,
        )
    };
    if !no_interp {
        let retval = if interp_sensi {
            idaLsGetY(ida_mem, tt, &mut yyTmp, &mut ypTmp, &mut yySTmp, &mut ypSTmp)
        } else {
            let mut noS: Vec<NVector> = Vec::new();
            let mut noSp: Vec<NVector> = Vec::new();
            idaLsGetY(ida_mem, tt, &mut yyTmp, &mut ypTmp, &mut noS, &mut noSp)
        };
        if retval != IDA_SUCCESS {
            let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
            adj.ia_yyTmp = yyTmp;
            adj.ia_ypTmp = ypTmp;
            adj.ia_yySTmp = yySTmp;
            adj.ia_ypSTmp = ypSTmp;
            IDAProcessError(Some(ida_mem), -1, line!(), "idaLsPrecSolveBS", file!(),
                            MSG_LS_BAD_T);
            return -1;
        }
    }

    /* Call user's adjoint psolveBS routine */
    let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
    let idaB_mem = &mut adj.IDAB_mem[which];
    let psolveBS = idaLsB_downcast(idaB_mem).psolveBS.unwrap();
    let retval = psolveBS(tt, &yyTmp, &ypTmp, &yySTmp, &ypSTmp, yyB, ypB, rrB, rvecB, zvecB,
                          c_jB, deltaB, &mut idaB_mem.ida_user_data);
    adj.ia_yyTmp = yyTmp;
    adj.ia_ypTmp = ypTmp;
    adj.ia_yySTmp = yySTmp;
    adj.ia_ypSTmp = ypSTmp;
    retval
}

/* idaLsJacTimesSetupB interfaces to the IDALsJacTimesSetupFnB
   routine provided by the user */
fn idaLsJacTimesSetupB(
    tt: f64,
    yyB: &NVector,
    ypB: &NVector,
    rrB: &NVector,
    c_jB: f64,
    ida_mem: &mut UserData,
) -> i32 {
    /* access relevant memory structures */
    let ida_mem = match idaLs_AccessIDAMem(ida_mem, "idaLsJacTimesSetupB") {
        Ok(m) => m,
        Err(e) => return e,
    };
    let which = match idaLs_AccessLMemBCur(ida_mem, "idaLsJacTimesSetupB") {
        Ok(w) => w,
        Err(e) => return e,
    };

    /* Get forward solution from interpolation. */
    let (mut yyTmp, mut ypTmp, no_interp) = {
        let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
        (std::mem::take(&mut adj.ia_yyTmp), std::mem::take(&mut adj.ia_ypTmp), adj.ia_noInterp)
    };
    if !no_interp {
        let mut noS: Vec<NVector> = Vec::new();
        let mut noSp: Vec<NVector> = Vec::new();
        let retval = idaLsGetY(ida_mem, tt, &mut yyTmp, &mut ypTmp, &mut noS, &mut noSp);
        if retval != IDA_SUCCESS {
            let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
            adj.ia_yyTmp = yyTmp;
            adj.ia_ypTmp = ypTmp;
            IDAProcessError(Some(ida_mem), -1, line!(), "idaLsJacTimesSetupB", file!(),
                            MSG_LS_BAD_T);
            return -1;
        }
    }

    /* Call user's adjoint jtsetupB routine */
    let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
    let idaB_mem = &mut adj.IDAB_mem[which];
    let jtsetupB = idaLsB_downcast(idaB_mem).jtsetupB.unwrap();
    let retval = jtsetupB(tt, &yyTmp, &ypTmp, yyB, ypB, rrB, c_jB, &mut idaB_mem.ida_user_data);
    adj.ia_yyTmp = yyTmp;
    adj.ia_ypTmp = ypTmp;
    retval
}

/* idaLsJacTimesSetupBS interfaces to the IDALsJacTimesSetupFnBS
   routine provided by the user */
fn idaLsJacTimesSetupBS(
    tt: f64,
    yyB: &NVector,
    ypB: &NVector,
    rrB: &NVector,
    c_jB: f64,
    ida_mem: &mut UserData,
) -> i32 {
    /* access relevant memory structures */
    let ida_mem = match idaLs_AccessIDAMem(ida_mem, "idaLsJacTimesSetupBS") {
        Ok(m) => m,
        Err(e) => return e,
    };
    let which = match idaLs_AccessLMemBCur(ida_mem, "idaLsJacTimesSetupBS") {
        Ok(w) => w,
        Err(e) => return e,
    };

    /* Get forward solution from interpolation. */
    let (mut yyTmp, mut ypTmp, mut yySTmp, mut ypSTmp, no_interp, interp_sensi) = {
        let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
        (
            std::mem::take(&mut adj.ia_yyTmp),
            std::mem::take(&mut adj.ia_ypTmp),
            std::mem::take(&mut adj.ia_yySTmp),
            std::mem::take(&mut adj.ia_ypSTmp),
            adj.ia_noInterp,
            adj.ia_interpSensi,
        )
    };
    if !no_interp {
        let retval = if interp_sensi {
            idaLsGetY(ida_mem, tt, &mut yyTmp, &mut ypTmp, &mut yySTmp, &mut ypSTmp)
        } else {
            let mut noS: Vec<NVector> = Vec::new();
            let mut noSp: Vec<NVector> = Vec::new();
            idaLsGetY(ida_mem, tt, &mut yyTmp, &mut ypTmp, &mut noS, &mut noSp)
        };
        if retval != IDA_SUCCESS {
            let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
            adj.ia_yyTmp = yyTmp;
            adj.ia_ypTmp = ypTmp;
            adj.ia_yySTmp = yySTmp;
            adj.ia_ypSTmp = ypSTmp;
            IDAProcessError(Some(ida_mem), -1, line!(), "idaLsJacTimesSetupBS", file!(),
                            MSG_LS_BAD_T);
            return -1;
        }
    }

    /* Call user's adjoint jtimesBS routine (C comment verbatim; this
       calls the jtsetupBS member, as the C body does) */
    let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
    let idaB_mem = &mut adj.IDAB_mem[which];
    let jtsetupBS = idaLsB_downcast(idaB_mem).jtsetupBS.unwrap();
    let retval = jtsetupBS(tt, &yyTmp, &ypTmp, &yySTmp, &ypSTmp, yyB, ypB, rrB, c_jB,
                           &mut idaB_mem.ida_user_data);
    adj.ia_yyTmp = yyTmp;
    adj.ia_ypTmp = ypTmp;
    adj.ia_yySTmp = yySTmp;
    adj.ia_ypSTmp = ypSTmp;
    retval
}

/* idaLsJacTimesVecB interfaces to the IDALsJacTimesVecFnB routine
   provided by the user */
#[allow(clippy::too_many_arguments)]
fn idaLsJacTimesVecB(
    tt: f64,
    yyB: &NVector,
    ypB: &NVector,
    rrB: &NVector,
    vB: &NVector,
    JvB: &mut NVector,
    c_jB: f64,
    ida_mem: &mut UserData,
    tmp1B: &mut NVector,
    tmp2B: &mut NVector,
) -> i32 {
    /* access relevant memory structures */
    let ida_mem = match idaLs_AccessIDAMem(ida_mem, "idaLsJacTimesVecB") {
        Ok(m) => m,
        Err(e) => return e,
    };
    let which = match idaLs_AccessLMemBCur(ida_mem, "idaLsJacTimesVecB") {
        Ok(w) => w,
        Err(e) => return e,
    };

    /* Get forward solution from interpolation. */
    let (mut yyTmp, mut ypTmp, no_interp) = {
        let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
        (std::mem::take(&mut adj.ia_yyTmp), std::mem::take(&mut adj.ia_ypTmp), adj.ia_noInterp)
    };
    if !no_interp {
        let mut noS: Vec<NVector> = Vec::new();
        let mut noSp: Vec<NVector> = Vec::new();
        let retval = idaLsGetY(ida_mem, tt, &mut yyTmp, &mut ypTmp, &mut noS, &mut noSp);
        if retval != IDA_SUCCESS {
            let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
            adj.ia_yyTmp = yyTmp;
            adj.ia_ypTmp = ypTmp;
            IDAProcessError(Some(ida_mem), -1, line!(), "idaLsJacTimesVecB", file!(),
                            MSG_LS_BAD_T);
            return -1;
        }
    }

    /* Call user's adjoint jtimesB routine */
    let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
    let idaB_mem = &mut adj.IDAB_mem[which];
    let jtimesB = idaLsB_downcast(idaB_mem).jtimesB.unwrap();
    let retval = jtimesB(tt, &yyTmp, &ypTmp, yyB, ypB, rrB, vB, JvB, c_jB,
                         &mut idaB_mem.ida_user_data, tmp1B, tmp2B);
    adj.ia_yyTmp = yyTmp;
    adj.ia_ypTmp = ypTmp;
    retval
}

/* idaLsJacTimesVecBS interfaces to the IDALsJacTimesVecFnBS routine
   provided by the user */
#[allow(clippy::too_many_arguments)]
fn idaLsJacTimesVecBS(
    tt: f64,
    yyB: &NVector,
    ypB: &NVector,
    rrB: &NVector,
    vB: &NVector,
    JvB: &mut NVector,
    c_jB: f64,
    ida_mem: &mut UserData,
    tmp1B: &mut NVector,
    tmp2B: &mut NVector,
) -> i32 {
    /* access relevant memory structures */
    let ida_mem = match idaLs_AccessIDAMem(ida_mem, "idaLsJacTimesVecBS") {
        Ok(m) => m,
        Err(e) => return e,
    };
    let which = match idaLs_AccessLMemBCur(ida_mem, "idaLsJacTimesVecBS") {
        Ok(w) => w,
        Err(e) => return e,
    };

    /* Get forward solution from interpolation. */
    let (mut yyTmp, mut ypTmp, mut yySTmp, mut ypSTmp, no_interp, interp_sensi) = {
        let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
        (
            std::mem::take(&mut adj.ia_yyTmp),
            std::mem::take(&mut adj.ia_ypTmp),
            std::mem::take(&mut adj.ia_yySTmp),
            std::mem::take(&mut adj.ia_ypSTmp),
            adj.ia_noInterp,
            adj.ia_interpSensi,
        )
    };
    if !no_interp {
        let retval = if interp_sensi {
            idaLsGetY(ida_mem, tt, &mut yyTmp, &mut ypTmp, &mut yySTmp, &mut ypSTmp)
        } else {
            let mut noS: Vec<NVector> = Vec::new();
            let mut noSp: Vec<NVector> = Vec::new();
            idaLsGetY(ida_mem, tt, &mut yyTmp, &mut ypTmp, &mut noS, &mut noSp)
        };
        if retval != IDA_SUCCESS {
            let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
            adj.ia_yyTmp = yyTmp;
            adj.ia_ypTmp = ypTmp;
            adj.ia_yySTmp = yySTmp;
            adj.ia_ypSTmp = ypSTmp;
            IDAProcessError(Some(ida_mem), -1, line!(), "idaLsJacTimesVecBS", file!(),
                            MSG_LS_BAD_T);
            return -1;
        }
    }

    /* Call user's adjoint jtimesBS routine */
    let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
    let idaB_mem = &mut adj.IDAB_mem[which];
    let jtimesBS = idaLsB_downcast(idaB_mem).jtimesBS.unwrap();
    let retval = jtimesBS(tt, &yyTmp, &ypTmp, &yySTmp, &ypSTmp, yyB, ypB, rrB, vB, JvB, c_jB,
                          &mut idaB_mem.ida_user_data, tmp1B, tmp2B);
    adj.ia_yyTmp = yyTmp;
    adj.ia_ypTmp = ypTmp;
    adj.ia_yySTmp = yySTmp;
    adj.ia_ypSTmp = ypSTmp;
    retval
}

/* idaLsFreeB frees memory associated with the IDASLS wrapper */
pub fn idaLsFreeB(idaB_mem: &mut IDABMem) -> i32 {
    /* free IDALsMemB interface structure (Rust: Drop) */
    idaB_mem.ida_lmem = None;
    IDALS_SUCCESS
}

/* IDAADJ_mem->ia_getY(IDA_mem, t, yy, yp, yyS, ypS): dispatch on
   ia_interpType per the pinned interpolation-module design
   (idas_impl.rs: the C ia_storePnt/ia_getY/ia_malloc/ia_free
   function pointers become dispatch on ia_interpType).  Empty
   yyS/ypS Vecs play the role of the C NULL arguments.

   FORWARD REFERENCE (pinned): the IDAAhermiteGetY /
   IDAApolynomialGetY families are implemented by the idaa.c port
   (idaa.rs).  Until idaa.rs lands, nothing in the workspace can
   construct an IDAadjMem (idaa.c owns IDAAdjInit), so this bridge
   is statically unreachable; the idaa.c unit MUST replace this body
   with the real dispatch, mirroring cvodes_ls.rs cvLsIMget
   (PROGRESS.md pin on the idas_ls.c entry). */
#[allow(unused_variables)]
/* (C: IDAADJ_mem->ia_getY — the ia_interpType dispatch to
   IDAAhermiteGetY / IDAApolynomialGetY lives in idaa.rs::IDAAgetY,
   per the recorded PIN; empty yyS/ypS Vecs play the C NULL inputs.) */
pub(crate) fn idaLsGetY(
    ida_mem: &mut IDAMem,
    t: f64,
    yy: &mut NVector,
    yp: &mut NVector,
    yyS: &mut Vec<NVector>,
    ypS: &mut Vec<NVector>,
) -> i32 {
    crate::idaa::IDAAgetY(ida_mem, t, yy, yp, yyS, ypS)
}

/* Downcast the IDAB_mem.ida_lmem Box<dyn Any> attachment to the
   IDALsMemB installed by IDASetLinearSolverB (guarded by the
   is::<IDALsMemB> checks in idaLs_AccessLMemB / idaLs_AccessLMemBCur). */
fn idaLsB_downcast(idaB_mem: &mut IDABMem) -> &mut IDALsMemB {
    idaB_mem.ida_lmem.as_mut().unwrap().downcast_mut::<IDALsMemB>().unwrap()
}

/* Recover the OUTER integrator memory from the inner (backward)
   problem's UserData (installed by idaa.rs, as the C idaa.c passes
   the outer ida_mem as the inner problem's user data). */
pub(crate) fn idaLs_AccessIDAMem<'a>(
    ida_mem: &'a mut UserData,
    fname: &str,
) -> Result<&'a mut IDAMem, i32> {
    match ida_mem.as_mut().and_then(|d| d.downcast_mut::<IDAMem>()) {
        Some(m) => Ok(m),
        None => {
            IDAProcessError(None, IDALS_MEM_NULL, line!(), fname, file!(), MSG_LS_IDAMEM_NULL);
            Err(IDALS_MEM_NULL)
        }
    }
}

/* idaLs_AccessLMemB checks the adjoint memory, `which`, and the
   IDALsMemB attachment, and returns the index of the backward problem
   in IDAADJ_mem.IDAB_mem (the C version unpacks IDA_mem/IDAADJ_mem/
   IDAB_mem/idalsB_mem pointers instead).  If any is missing it returns
   IDALS_MEM_NULL, IDALS_NO_ADJ, IDALS_ILL_INPUT, or IDALS_LMEMB_NULL. */
pub(crate) fn idaLs_AccessLMemB(
    ida_mem: &mut IDAMem,
    which: i32,
    fname: &str,
) -> Result<usize, i32> {
    /* access IDAadjMem structure */
    if !ida_mem.ida_adjMallocDone {
        IDAProcessError(Some(ida_mem), IDALS_NO_ADJ, line!(), fname, file!(), MSG_LS_NO_ADJ);
        return Err(IDALS_NO_ADJ);
    }
    let idaadj_mem = ida_mem.ida_adj_mem.as_ref().unwrap();

    /* Check the value of which */
    if which >= idaadj_mem.ia_nbckpbs {
        IDAProcessError(None, IDALS_ILL_INPUT, line!(), fname, file!(), MSG_LS_BAD_WHICH);
        return Err(IDALS_ILL_INPUT);
    }

    /* Find the IDABMem entry in the linked list corresponding to which */
    let idx = idaadj_mem.IDAB_mem.iter().position(|b| b.ida_index == which).unwrap();

    /* access IDALsMemB structure */
    let has_lmemB = idaadj_mem.IDAB_mem[idx]
        .ida_lmem
        .as_ref()
        .map(|l| l.is::<IDALsMemB>())
        .unwrap_or(false);
    if !has_lmemB {
        IDAProcessError(None, IDALS_LMEMB_NULL, line!(), fname, file!(), MSG_LS_LMEMB_NULL);
        return Err(IDALS_LMEMB_NULL);
    }

    Ok(idx)
}

/* idaLs_AccessLMemBCur returns the index of the currently integrated
   backward problem (ia_bckpbCrt) after the same checks.  If any piece
   is missing it returns IDALS_MEM_NULL, IDALS_NO_ADJ, or
   IDALS_LMEMB_NULL. */
pub(crate) fn idaLs_AccessLMemBCur(ida_mem: &mut IDAMem, fname: &str) -> Result<usize, i32> {
    /* access IDAadjMem structure */
    if !ida_mem.ida_adjMallocDone {
        IDAProcessError(Some(ida_mem), IDALS_NO_ADJ, line!(), fname, file!(), MSG_LS_NO_ADJ);
        return Err(IDALS_NO_ADJ);
    }
    let idaadj_mem = ida_mem.ida_adj_mem.as_ref().unwrap();

    /* get current backward problem */
    let which = match idaadj_mem.ia_bckpbCrt {
        Some(w) => w,
        None => {
            IDAProcessError(None, IDALS_LMEMB_NULL, line!(), fname, file!(), MSG_LS_LMEMB_NULL);
            return Err(IDALS_LMEMB_NULL);
        }
    };

    /* access IDALsMemB structure */
    let has_lmemB = idaadj_mem.IDAB_mem[which]
        .ida_lmem
        .as_ref()
        .map(|l| l.is::<IDALsMemB>())
        .unwrap_or(false);
    if !has_lmemB {
        IDAProcessError(None, IDALS_LMEMB_NULL, line!(), fname, file!(), MSG_LS_LMEMB_NULL);
        return Err(IDALS_LMEMB_NULL);
    }

    Ok(which)
}

/*===============================================================
  Tests
  ===============================================================*/
#[cfg(test)]
mod tests {
    use super::*;
    use crate::sundials_context::SUNContext;
    use crate::sunlinsol_dense::SUNLinSol_Dense;
    use crate::sunmatrix_dense::SUNDenseMatrix;

    /* Linear DAE residual with a known Jacobian:
         F0 = yp0 + 2*y0 + 1*y1 - 1
         F1 = yp1 + 1*y0 + 3*y1 - 2
       => dF/dy = [[2, 1], [1, 3]], dF/dyp = I, so
          J = dF/dy + cj*dF/dyp = [[2 + cj, 1], [1, 3 + cj]]. */
    fn resfn(_t: f64, yy: &NVector, yp: &NVector, rr: &mut NVector, _ud: &mut UserData) -> i32 {
        rr.data[0] = yp.data[0] + 2.0 * yy.data[0] + yy.data[1] - 1.0;
        rr.data[1] = yp.data[1] + yy.data[0] + 3.0 * yy.data[1] - 2.0;
        0
    }

    fn make_ida_mem(n: usize) -> IDAMem {
        let mut ida_mem = IDAMem::default();
        ida_mem.ida_res = Some(resfn);
        ida_mem.ida_ewt = NVector::from_slice(&vec![1.0; n]);
        ida_mem.ida_tempv1 = NVector::new(n);
        ida_mem.ida_tempv2 = NVector::new(n);
        ida_mem.ida_tempv3 = NVector::new(n);
        ida_mem.ida_hh = 1.0e-2;
        ida_mem.ida_cj = 2.0;
        ida_mem.ida_cjratio = 1.0;
        ida_mem
    }

    /* IDASetLinearSolver sets the C defaults (idas_ls.c):
       iterative/matrixbased flags, jacDQ/jtimesDQ, eplifac = 0.05,
       dqincfac = 1, scalesol for matrix-based LS, counters zeroed. */
    #[test]
    fn idasetlinearsolver_defaults() {
        let sunctx = SUNContext::default();
        let mut ida_mem = make_ida_mem(2);

        let a = SUNDenseMatrix(2, 2, &sunctx);
        let ls = SUNLinSol_Dense(&ida_mem.ida_tempv1, &a, &sunctx);

        let retval = IDASetLinearSolver(&mut ida_mem, ls, Some(a));
        assert_eq!(retval, IDALS_SUCCESS);

        let idals_mem = idaLs_AccessLMem(&mut ida_mem, "test").unwrap();
        assert!(!idals_mem.iterative); /* direct solver */
        assert!(idals_mem.matrixbased);
        assert!(idals_mem.jacDQ);
        assert!(idals_mem.jtimesDQ);
        assert!(idals_mem.jtsetup.is_none());
        assert!(idals_mem.pset.is_none() && idals_mem.psolve.is_none());
        assert_eq!(idals_mem.eplifac, 0.05);
        assert_eq!(idals_mem.dqincfac, 1.0);
        assert_eq!(idals_mem.nrmfac, 0.0); /* direct: not computed */
        assert!(idals_mem.scalesol); /* matrix-based */
        assert_eq!(idals_mem.last_flag, IDALS_SUCCESS);
        assert_eq!(idals_mem.nje, 0);
        assert_eq!(idals_mem.nreDQ, 0);
        assert!(!idals_mem.setup_disabled);
        assert_eq!(idals_mem.ytemp.len(), 2);
        assert_eq!(idals_mem.yptemp.len(), 2);
        assert_eq!(idals_mem.x.len(), 2);
    }

    /* direct LS without a matrix is rejected (IDALS_ILL_INPUT) */
    #[test]
    fn idasetlinearsolver_direct_needs_matrix() {
        let sunctx = SUNContext::default();
        let mut ida_mem = make_ida_mem(2);
        let a = SUNDenseMatrix(2, 2, &sunctx);
        let ls = SUNLinSol_Dense(&ida_mem.ida_tempv1, &a, &sunctx);
        assert_eq!(IDASetLinearSolver(&mut ida_mem, ls, None), IDALS_ILL_INPUT);
        assert!(ida_mem.ida_lmem.is_none());
    }

    /* idaLsSetup drives the dense DQ Jacobian: for the linear test DAE
       the exact Jacobian is dF/dy + cj*dF/dyp = [[2+cj, 1], [1, 3+cj]];
       the difference quotient must agree to ~sqrt(uround), and the
       nje/nreDQ/nstlj/tnlj updates follow idas_ls.c idaLsSetup. */
    #[test]
    fn idalssetup_dense_dq_jacobian() {
        let sunctx = SUNContext::default();
        let mut ida_mem = make_ida_mem(2);
        ida_mem.ida_nst = 5;
        ida_mem.ida_tn = 0.25;
        let cj = ida_mem.ida_cj; /* 2.0 */

        let yy = NVector::from_slice(&[0.5, -0.25]);
        let yp = NVector::from_slice(&[0.1, 0.2]);
        let mut rr = NVector::new(2);
        {
            let IDAMem { ida_user_data, .. } = &mut ida_mem;
            resfn(0.25, &yy, &yp, &mut rr, ida_user_data);
        }

        let a = SUNDenseMatrix(2, 2, &sunctx);
        let ls = SUNLinSol_Dense(&ida_mem.ida_tempv1, &a, &sunctx);
        assert_eq!(IDASetLinearSolver(&mut ida_mem, ls, Some(a)), IDALS_SUCCESS);

        /* take the module out for the call, donor dispatch pattern */
        let mut lmem = std::mem::take(&mut ida_mem.ida_lmem);
        let idals_mem = match &mut lmem {
            LsModule::Ls(m) => m,
            _ => unreachable!(),
        };
        assert_eq!(idaLsInit(&mut ida_mem, idals_mem), 0);
        assert!(!idals_mem.setup_disabled);

        /* raw difference-quotient Jacobian (before factorization) */
        {
            let j = idals_mem.J.as_mut().unwrap();
            assert_eq!(
                idaLsDQJac(&mut ida_mem, &mut idals_mem.nreDQ, 0.25, cj, &yy, &yp, &rr, j),
                0
            );
            let jac = match &*j {
                SUNMatrix::Dense(dm) => dm,
                _ => unreachable!(),
            };
            let tol = 1.0e-6;
            assert!((jac.get(0, 0) - (2.0 + cj)).abs() < tol);
            assert!((jac.get(0, 1) - 1.0).abs() < tol);
            assert!((jac.get(1, 0) - 1.0).abs() < tol);
            assert!((jac.get(1, 1) - (3.0 + cj)).abs() < tol);
        }
        assert_eq!(idals_mem.nreDQ, 2); /* one res eval per column */
        idals_mem.nreDQ = 0;

        assert_eq!(idaLsSetup(&mut ida_mem, idals_mem, &yy, &yp, &rr), 0);
        assert_eq!(idals_mem.nje, 1);
        assert_eq!(idals_mem.nreDQ, 2);
        assert_eq!(idals_mem.nstlj, 5);
        assert_eq!(idals_mem.tnlj, 0.25);
        assert_eq!(idals_mem.last_flag, SUN_SUCCESS);

        ida_mem.ida_lmem = lmem;

        /* stats getters see the accumulated counters */
        let mut nje = -1;
        assert_eq!(IDAGetNumJacEvals(&mut ida_mem, &mut nje), IDALS_SUCCESS);
        assert_eq!(nje, 1);
        let mut nre = -1;
        assert_eq!(IDAGetNumLinResEvals(&mut ida_mem, &mut nre), IDALS_SUCCESS);
        assert_eq!(nre, 2);
        let mut nstj = -1;
        assert_eq!(IDAGetJacNumSteps(&mut ida_mem, &mut nstj), IDALS_SUCCESS);
        assert_eq!(nstj, 5);
        let mut tj = -1.0;
        assert_eq!(IDAGetJacTime(&mut ida_mem, &mut tj), IDALS_SUCCESS);
        assert_eq!(tj, 0.25);
    }

    /* idaLsSolve (direct): solves J x = b after setup; with
       cjratio == 1 no rescaling is applied, then with cjratio != 1
       the correction is scaled by 2/(1 + cjratio) exactly as
       idas_ls.c idaLsSolve. */
    #[test]
    fn idalssolve_direct_round_trip() {
        let sunctx = SUNContext::default();
        let mut ida_mem = make_ida_mem(2);
        let cj = ida_mem.ida_cj; /* J = [[4, 1], [1, 5]] */

        let yy = NVector::from_slice(&[0.5, -0.25]);
        let yp = NVector::from_slice(&[0.1, 0.2]);
        let mut rr = NVector::new(2);
        {
            let IDAMem { ida_user_data, .. } = &mut ida_mem;
            resfn(0.0, &yy, &yp, &mut rr, ida_user_data);
        }

        let a = SUNDenseMatrix(2, 2, &sunctx);
        let ls = SUNLinSol_Dense(&ida_mem.ida_tempv1, &a, &sunctx);
        assert_eq!(IDASetLinearSolver(&mut ida_mem, ls, Some(a)), IDALS_SUCCESS);

        let mut lmem = std::mem::take(&mut ida_mem.ida_lmem);
        let idals_mem = match &mut lmem {
            LsModule::Ls(m) => m,
            _ => unreachable!(),
        };
        assert_eq!(idaLsInit(&mut ida_mem, idals_mem), 0);
        assert_eq!(idaLsSetup(&mut ida_mem, idals_mem, &yy, &yp, &rr), 0);

        /* J = [[2+cj, 1], [1, 3+cj]] with cj = 2: solve J x = [9, 11]
           -> x = [(9*5 - 11)/19, (11*4 - 9)/19] = [34/19, 35/19] */
        assert_eq!(cj, 2.0);
        let mut b = NVector::from_slice(&[9.0, 11.0]);
        let weight = NVector::from_slice(&[1.0, 1.0]);
        let retval = idaLsSolve(&mut ida_mem, idals_mem, &mut b, &weight, &yy, &yp, &rr);
        assert_eq!(retval, 0);
        assert_eq!(idals_mem.last_flag, SUN_SUCCESS);
        assert_eq!(idals_mem.nli, 0);
        assert_eq!(idals_mem.ncfl, 0);
        assert!((b.data[0] - 34.0 / 19.0).abs() < 1.0e-9);
        assert!((b.data[1] - 35.0 / 19.0).abs() < 1.0e-9);

        /* cjratio != 1 with scalesol: b is scaled by 2/(1 + cjratio) */
        ida_mem.ida_cjratio = 0.5;
        let mut b2 = NVector::from_slice(&[9.0, 11.0]);
        let retval = idaLsSolve(&mut ida_mem, idals_mem, &mut b2, &weight, &yy, &yp, &rr);
        assert_eq!(retval, 0);
        let scale = 2.0 / 1.5;
        assert!((b2.data[0] - scale * 34.0 / 19.0).abs() < 1.0e-9);
        assert!((b2.data[1] - scale * 35.0 / 19.0).abs() < 1.0e-9);

        ida_mem.ida_lmem = lmem;

        let mut flag = -1;
        assert_eq!(IDAGetLastLinFlag(&mut ida_mem, &mut flag), IDALS_SUCCESS);
        assert_eq!(flag, SUN_SUCCESS as i64);
    }

    /* the internal DQ Jv product matches J*v for the linear test DAE
       (idaLsDQJtimes, idas_ls.c) */
    #[test]
    fn idalsdqjtimes_matches_analytic() {
        let mut ida_mem = make_ida_mem(2);
        let cj = ida_mem.ida_cj; /* J = [[4, 1], [1, 5]] */

        let yy = NVector::from_slice(&[0.5, -0.25]);
        let yp = NVector::from_slice(&[0.1, 0.2]);
        let mut rr = NVector::new(2);
        {
            let IDAMem { ida_user_data, .. } = &mut ida_mem;
            resfn(0.0, &yy, &yp, &mut rr, ida_user_data);
        }

        let ewt = NVector::from_slice(&[1.0, 1.0]);
        let v = NVector::from_slice(&[1.0, -2.0]);
        let mut jv = NVector::new(2);
        let mut work1 = NVector::new(2);
        let mut work2 = NVector::new(2);
        let mut nreDQ = 0i64;

        /* non-GMRES sigma: dqincfac / ||v||_WRMS */
        let retval = idaLsDQJtimes(&mut ida_mem, &ewt, &mut nreDQ, false, 0.0, 1.0,
                                   Some(resfn), 0.0, &yy, &yp, &rr, &v, &mut jv, cj,
                                   &mut work1, &mut work2);
        assert_eq!(retval, 0);
        assert_eq!(nreDQ, 1);
        /* J*v = [4*1 + 1*(-2), 1*1 + 5*(-2)] = [2, -9] (exact: linear F) */
        assert!((jv.data[0] - 2.0).abs() < 1.0e-8);
        assert!((jv.data[1] + 9.0).abs() < 1.0e-8);
    }

    #[test]
    fn return_flag_names() {
        assert_eq!(IDAGetLinReturnFlagName(IDALS_SUCCESS as i64), "IDALS_SUCCESS");
        assert_eq!(
            IDAGetLinReturnFlagName(IDALS_JACFUNC_RECVR as i64),
            "IDALS_JACFUNC_RECVR"
        );
        assert_eq!(IDAGetLinReturnFlagName(1234), "NONE");
    }

    /* PART II guards: the B set routines reject a forward problem with
       no adjoint module initialized (IDALS_NO_ADJ), per idas_ls.c
       idaLs_AccessLMemB / IDASetLinearSolverB. */
    #[test]
    fn backward_set_routines_require_adjoint() {
        let sunctx = SUNContext::default();
        let mut ida_mem = make_ida_mem(2);

        let a = SUNDenseMatrix(2, 2, &sunctx);
        let ls = SUNLinSol_Dense(&ida_mem.ida_tempv1, &a, &sunctx);
        assert_eq!(IDASetLinearSolverB(&mut ida_mem, 0, ls, Some(a)), IDALS_NO_ADJ);

        assert_eq!(IDASetJacFnB(&mut ida_mem, 0, None), IDALS_NO_ADJ);
        assert_eq!(IDASetEpsLinB(&mut ida_mem, 0, 0.1), IDALS_NO_ADJ);
        assert_eq!(IDASetLSNormFactorB(&mut ida_mem, 0, 0.0), IDALS_NO_ADJ);
        assert_eq!(IDASetLinearSolutionScalingB(&mut ida_mem, 0, true), IDALS_NO_ADJ);
        assert_eq!(IDASetIncrementFactorB(&mut ida_mem, 0, 1.0), IDALS_NO_ADJ);
        assert_eq!(IDASetPreconditionerB(&mut ida_mem, 0, None, None), IDALS_NO_ADJ);
        assert_eq!(IDASetJacTimesB(&mut ida_mem, 0, None, None), IDALS_NO_ADJ);
        assert_eq!(IDASetJacTimesResFnB(&mut ida_mem, 0, None), IDALS_NO_ADJ);
    }
}
