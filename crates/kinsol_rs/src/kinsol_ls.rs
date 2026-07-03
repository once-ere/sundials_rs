/* -----------------------------------------------------------------
 * Translated from src/kinsol/kinsol_ls.c (KINSOL 7.7.0).
 * KINSOL's linear solver interface (KINLS): attach a SUNLinearSolver,
 * Jacobian and preconditioner plumbing, difference-quotient Jacobian
 * and Jacobian-times-vector approximations, and the generic
 * linit/lsetup/lsolve routines dispatched from KINSol.
 * Translation conventions follow the donor cvode_ls.rs.
 * -----------------------------------------------------------------*/
use std::cell::RefCell;

use crate::kinsol_bbdpre::{KINBBDPrecSetup, KINBBDPrecSolve};
use crate::kinsol_impl::*;
use crate::kinsol_ls_impl::*;
use crate::nvector_serial::*;
use crate::sundials_errors::*;
use crate::sundials_linearsolver::*;
use crate::sundials_math::*;
use crate::sundials_matrix::*;
use crate::sundials_types::*;

/* constants */
const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/*---------------------------------------------------------------
  kinLs_AccessLMem

  This routine unpacks the ls_mem structure (KINLsMem) from the
  solver memory.  If it is missing it returns KINLS_LMEM_NULL
  (the kin_mem NULL check of the C original cannot arise: kinmem
  is a &mut KINMem here).
  ---------------------------------------------------------------*/
pub fn kinLs_AccessLMem<'a>(
    kin_mem: &'a mut KINMem,
    fname: &str,
) -> Result<&'a mut KINLsMem, i32> {
    match &mut kin_mem.kin_lmem {
        LsModule::Ls(ls) => Ok(ls),
        _ => {
            KINProcessError(None, KINLS_LMEM_NULL, line!(), fname, file!(), MSG_LS_LMEM_NULL);
            Err(KINLS_LMEM_NULL)
        }
    }
}

/*==================================================================
  KINLS Exported functions -- Required
  ==================================================================*/

/*---------------------------------------------------------------
  KINSetLinearSolver specifies the linear solver
  ---------------------------------------------------------------*/
pub fn KINSetLinearSolver(kin_mem: &mut KINMem, LS: LinearSolver, A: Option<SUNMatrix>) -> i32 {
    /* Retrieve the LS type */
    let ls_type = LS.ls_type();

    /* Return with error if LS has 'matrix-embedded' type */
    if ls_type == SUNLINEARSOLVER_MATRIX_EMBEDDED {
        KINProcessError(Some(kin_mem), KINLS_ILL_INPUT, line!(), "KINSetLinearSolver", file!(),
                        "KINSOL is incompatible with MATRIX_EMBEDDED LS objects");
        return KINLS_ILL_INPUT;
    }

    /* Set flags based on LS type */
    let iterative = ls_type != SUNLINEARSOLVER_DIRECT;
    let matrixbased = ls_type != SUNLINEARSOLVER_ITERATIVE;

    /* (Required serial N_Vector operations — nvconst, nvdotprod,
       nvgetlength — are always implemented; the C capability checks
       cannot fail.  Likewise every workspace iterative solver takes
       the ATimes closure, so the setatimes support check holds.) */

    /* Check for compatible LS type, matrix and "atimes" support */
    if iterative {
        if matrixbased && A.is_none() {
            KINProcessError(Some(kin_mem), KINLS_ILL_INPUT, line!(), "KINSetLinearSolver", file!(),
                            "Incompatible inputs: matrix-iterative LS requires non-NULL matrix");
            return KINLS_ILL_INPUT;
        }
    } else if A.is_none() {
        KINProcessError(Some(kin_mem), KINLS_ILL_INPUT, line!(), "KINSetLinearSolver", file!(),
                        "Incompatible inputs: direct LS requires non-NULL matrix");
        return KINLS_ILL_INPUT;
    }

    /* free any existing system solver attached to KIN
       (RAII: dropped when kin_lmem is overwritten below) */

    /* Determine if this is an iterative linear solver */
    kin_mem.kin_inexact_ls = iterative;

    /* In C the four system linear solver function fields
       (kin_linit/kin_lsetup/kin_lsolve/kin_lfree) are set here; the
       Rust port dispatches them through the LsModule::Ls variant. */

    /* Get memory for KINLsMemRec and set defaults (C memset(0) +
       explicit assignments) */
    let kinls_mem = Box::new(KINLsMem {
        iterative,
        matrixbased,

        /* Set defaults for Jacobian-related fields */
        jacDQ: A.is_some(),
        jac: None, /* None + jacDQ => internal kinLsDQJac */

        /* set SUNLinearSolver pointer */
        LS,

        /* set SUNMatrix pointer (can be None) */
        J: A,

        /* initialize tolerance scaling factor */
        tol_fac: -ONE,

        /* Initialize counters (kinLsInitializeCounters) */
        nje: 0,
        nfeDQ: 0,
        npe: 0,
        nli: 0,
        nps: 0,
        ncfl: 0,
        njtimes: 0,

        new_uu: SUNFALSE,

        /* Set default values for the rest of the LS parameters */
        last_flag: KINLS_SUCCESS,

        /* Set defaults for preconditioner-related fields */
        pset: None,
        psolve: None,
        prec_module: PrecModule::None,

        /* Jacobian-times-vector: internal DQ by default */
        jtimesDQ: SUNTRUE,
        jtimes: None, /* None + jtimesDQ => internal kinLsDQJtimes */
        jt_func: kin_mem.kin_func,

        setup_disabled: SUNFALSE,
    });

    /* If LS supports ATimes/preconditioning, the KINLs routines are
       supplied as closures at solve time (kinLsSolveIterative); the C
       SUNLinSolSetATimes / SUNLinSolSetPreconditioner registration
       calls have no Rust counterpart. */

    /* Attach linear solver memory to integrator memory */
    kin_mem.kin_lmem = LsModule::Ls(kinls_mem);

    KINLS_SUCCESS
}

/*==================================================================
  Optional Set routines
  ==================================================================*/

/*------------------------------------------------------------------
  KINSetJacFn specifies the Jacobian function
  ------------------------------------------------------------------*/
pub fn KINSetJacFn(kin_mem: &mut KINMem, jac: Option<KINLsJacFn>) -> i32 {
    /* access KINLsMem structure */
    let kinls_mem = match kinLs_AccessLMem(kin_mem, "KINSetJacFn") {
        Ok(m) => m,
        Err(e) => return e,
    };

    /* return with failure if jac cannot be used */
    if jac.is_some() && kinls_mem.J.is_none() {
        KINProcessError(None, KINLS_ILL_INPUT, line!(), "KINSetJacFn", file!(),
                        "Jacobian routine cannot be supplied for NULL SUNMatrix");
        return KINLS_ILL_INPUT;
    }

    if let Some(j) = jac {
        kinls_mem.jacDQ = SUNFALSE;
        kinls_mem.jac = Some(j);
    } else {
        kinls_mem.jacDQ = SUNTRUE;
        kinls_mem.jac = None; /* internal kinLsDQJac */
    }

    KINLS_SUCCESS
}

/*------------------------------------------------------------------
  KINSetPreconditioner sets the preconditioner setup and solve
  functions
  ------------------------------------------------------------------*/
pub fn KINSetPreconditioner(
    kin_mem: &mut KINMem,
    psetup: Option<KINLsPrecSetupFn>,
    psolve: Option<KINLsPrecSolveFn>,
) -> i32 {
    /* access KINLsMem structure */
    let kinls_mem = match kinLs_AccessLMem(kin_mem, "KINSetPreconditioner") {
        Ok(m) => m,
        Err(e) => return e,
    };

    /* store function pointers for user-supplied routines in KINLS interface */
    kinls_mem.pset = psetup;
    kinls_mem.psolve = psolve;
    kinls_mem.prec_module = PrecModule::User;

    /* issue error if LS object does not support user-supplied
       preconditioning (C: LS->ops->setpreconditioner == NULL, true of
       the direct solvers) */
    if matches!(
        kinls_mem.LS,
        LinearSolver::Dense(_) | LinearSolver::Band(_) | LinearSolver::Custom(_)
    ) {
        KINProcessError(None, KINLS_ILL_INPUT, line!(), "KINSetPreconditioner", file!(),
                        "SUNLinearSolver object does not support user-supplied preconditioning");
        return KINLS_ILL_INPUT;
    }

    /* notify iterative linear solver to call KINLs interface routines:
       done at solve time via the kinLsSolveIterative closures */

    KINLS_SUCCESS
}

/*------------------------------------------------------------------
  KINSetJacTimesVecFn sets the matrix-vector product function
  ------------------------------------------------------------------*/
pub fn KINSetJacTimesVecFn(kin_mem: &mut KINMem, jtv: Option<KINLsJacTimesVecFn>) -> i32 {
    let kin_func = kin_mem.kin_func;

    /* access KINLsMem structure */
    let kinls_mem = match kinLs_AccessLMem(kin_mem, "KINSetJacTimesVecFn") {
        Ok(m) => m,
        Err(e) => return e,
    };

    /* issue error if LS object does not support user-supplied ATimes
       (C: LS->ops->setatimes == NULL, true of the direct solvers) */
    if matches!(
        kinls_mem.LS,
        LinearSolver::Dense(_) | LinearSolver::Band(_) | LinearSolver::Custom(_)
    ) {
        KINProcessError(None, KINLS_ILL_INPUT, line!(), "KINSetJacTimesVecFn", file!(),
                        "SUNLinearSolver object does not support user-supplied ATimes routine");
        return KINLS_ILL_INPUT;
    }

    /* store function pointers for user-supplied routine in KINLs
       interface (NULL jtimes implies use of DQ default) */
    if let Some(f) = jtv {
        kinls_mem.jtimesDQ = SUNFALSE;
        kinls_mem.jtimes = Some(f);
    } else {
        kinls_mem.jtimesDQ = SUNTRUE;
        kinls_mem.jtimes = None; /* internal kinLsDQJtimes */
        kinls_mem.jt_func = kin_func;
    }

    KINLS_SUCCESS
}

/* KINSetJacTimesVecSysFn specifies an alternative user-supplied system function
   to use in the internal finite difference Jacobian-vector product */
pub fn KINSetJacTimesVecSysFn(kin_mem: &mut KINMem, jtimesSysFn: Option<KINSysFn>) -> i32 {
    let kin_func = kin_mem.kin_func;

    /* access KINLsMem structure */
    let kinls_mem = match kinLs_AccessLMem(kin_mem, "KINSetJacTimesVecSysFn") {
        Ok(m) => m,
        Err(e) => return e,
    };

    /* check if using internal finite difference approximation */
    if !kinls_mem.jtimesDQ {
        KINProcessError(None, KINLS_ILL_INPUT, line!(), "KINSetJacTimesVecSysFn", file!(),
                        "Internal finite-difference Jacobian-vector product is disabled.");
        return KINLS_ILL_INPUT;
    }

    /* store function pointers for system function (NULL implies use kin_func) */
    match jtimesSysFn {
        Some(f) => kinls_mem.jt_func = Some(f),
        None => kinls_mem.jt_func = kin_func,
    }

    KINLS_SUCCESS
}

/*==================================================================
  Optional Get routines
  ==================================================================*/

pub fn KINGetJac<'a>(kin_mem: &'a mut KINMem, j: &mut Option<&'a SUNMatrix>) -> i32 {
    /* access KINLsMem structure; set output and return */
    let kinls_mem = match kinLs_AccessLMem(kin_mem, "KINGetJac") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *j = kinls_mem.J.as_ref();
    KINLS_SUCCESS
}

pub fn KINGetJacNumIters(kin_mem: &mut KINMem, nni_J: &mut i64) -> i32 {
    let nnilset = kin_mem.kin_nnilset;

    /* access KINLsMem structure; set output and return */
    if let Err(e) = kinLs_AccessLMem(kin_mem, "KINGetJacNumIters") {
        return e;
    }
    *nni_J = nnilset;
    KINLS_SUCCESS
}

/*------------------------------------------------------------------
  KINGetLinWorkSpace returns the integer and real workspace size
  ------------------------------------------------------------------*/
pub fn KINGetLinWorkSpace(kin_mem: &mut KINMem, lenrwLS: &mut i64, leniwLS: &mut i64) -> i32 {
    let (lrw1, liw1) = N_VSpace(&kin_mem.kin_vtemp1);

    /* access KINLsMem structure */
    let kinls_mem = match kinLs_AccessLMem(kin_mem, "KINGetLinWorkSpace") {
        Ok(m) => m,
        Err(e) => return e,
    };

    /* start with fixed sizes plus vector/matrix pointers */
    *lenrwLS = 1;
    *leniwLS = 21;

    /* add N_Vector sizes */
    *lenrwLS += lrw1;
    *leniwLS += liw1;

    /* add LS sizes */
    {
        let (lrw, liw) = kinls_mem.LS.space();
        *lenrwLS += lrw;
        *leniwLS += liw;
    }

    KINLS_SUCCESS
}

/*------------------------------------------------------------------
  KINGetNumJacEvals returns the number of Jacobian evaluations
  ------------------------------------------------------------------*/
pub fn KINGetNumJacEvals(kin_mem: &mut KINMem, njevals: &mut i64) -> i32 {
    let kinls_mem = match kinLs_AccessLMem(kin_mem, "KINGetNumJacEvals") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *njevals = kinls_mem.nje;
    KINLS_SUCCESS
}

/*------------------------------------------------------------------
  KINGetNumPrecEvals returns the total number of preconditioner
  evaluations
  ------------------------------------------------------------------*/
pub fn KINGetNumPrecEvals(kin_mem: &mut KINMem, npevals: &mut i64) -> i32 {
    let kinls_mem = match kinLs_AccessLMem(kin_mem, "KINGetNumPrecEvals") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *npevals = kinls_mem.npe;
    KINLS_SUCCESS
}

/*------------------------------------------------------------------
  KINGetNumPrecSolves returns the total number of times the
  preconditioner was applied
  ------------------------------------------------------------------*/
pub fn KINGetNumPrecSolves(kin_mem: &mut KINMem, npsolves: &mut i64) -> i32 {
    let kinls_mem = match kinLs_AccessLMem(kin_mem, "KINGetNumPrecSolves") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *npsolves = kinls_mem.nps;
    KINLS_SUCCESS
}

/*------------------------------------------------------------------
  KINGetNumLinIters returns the total number of linear
  iterations
  ------------------------------------------------------------------*/
pub fn KINGetNumLinIters(kin_mem: &mut KINMem, nliters: &mut i64) -> i32 {
    let kinls_mem = match kinLs_AccessLMem(kin_mem, "KINGetNumLinIters") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *nliters = kinls_mem.nli;
    KINLS_SUCCESS
}

/*------------------------------------------------------------------
  KINGetNumLinConvFails returns the total number of convergence
  failures
  ------------------------------------------------------------------*/
pub fn KINGetNumLinConvFails(kin_mem: &mut KINMem, nlcfails: &mut i64) -> i32 {
    let kinls_mem = match kinLs_AccessLMem(kin_mem, "KINGetNumLinConvFails") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *nlcfails = kinls_mem.ncfl;
    KINLS_SUCCESS
}

/*------------------------------------------------------------------
  KINGetNumJtimesEvals returns the number of times the matrix
  vector product was computed
  ------------------------------------------------------------------*/
pub fn KINGetNumJtimesEvals(kin_mem: &mut KINMem, njvevals: &mut i64) -> i32 {
    let kinls_mem = match kinLs_AccessLMem(kin_mem, "KINGetNumJtimesEvals") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *njvevals = kinls_mem.njtimes;
    KINLS_SUCCESS
}

/*------------------------------------------------------------------
  KINGetNumLinFuncEvals returns the number of calls to the user's
  F routine by the linear solver module
  ------------------------------------------------------------------*/
pub fn KINGetNumLinFuncEvals(kin_mem: &mut KINMem, nfevals: &mut i64) -> i32 {
    let kinls_mem = match kinLs_AccessLMem(kin_mem, "KINGetNumLinFuncEvals") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *nfevals = kinls_mem.nfeDQ;
    KINLS_SUCCESS
}

/*------------------------------------------------------------------
  KINGetLastLinFlag returns the last flag set in the KINLS
  function
  ------------------------------------------------------------------*/
pub fn KINGetLastLinFlag(kin_mem: &mut KINMem, flag: &mut i64) -> i32 {
    let kinls_mem = match kinLs_AccessLMem(kin_mem, "KINGetLastLinFlag") {
        Ok(m) => m,
        Err(e) => return e,
    };
    *flag = kinls_mem.last_flag as i64;
    KINLS_SUCCESS
}

/*------------------------------------------------------------------
  KINGetLinReturnFlagName
  ------------------------------------------------------------------*/
pub fn KINGetLinReturnFlagName(flag: i64) -> String {
    match flag as i32 {
        KINLS_SUCCESS => "KINLS_SUCCESS",
        KINLS_MEM_NULL => "KINLS_MEM_NULL",
        KINLS_LMEM_NULL => "KINLS_LMEM_NULL",
        KINLS_ILL_INPUT => "KINLS_ILL_INPUT",
        KINLS_MEM_FAIL => "KINLS_MEM_FAIL",
        KINLS_PMEM_NULL => "KINLS_PMEM_NULL",
        KINLS_JACFUNC_ERR => "KINLS_JACFUNC_ERR",
        KINLS_SUNMAT_FAIL => "KINLS_SUNMAT_FAIL",
        KINLS_SUNLS_FAIL => "KINLS_SUNLS_FAIL",
        _ => "NONE",
    }
    .to_string()
}

/*==================================================================
  KINLS Private functions
  ==================================================================*/

/*------------------------------------------------------------------
  kinLsATimes

  This routine coordinates the generation of the matrix-vector
  product z = J*v by calling either kinLsDQJtimes, which uses
  a difference quotient approximation for J*v, or by calling the
  user-supplied routine KINLsJacTimesVecFn if it is non-null.
  (kinls_mem is threaded explicitly per the workspace take()
  convention; kinLsSolveIterative replicates this body inside the
  ATimes closure handed to the iterative solvers.)
  ------------------------------------------------------------------*/
pub fn kinLsATimes(
    kin_mem: &mut KINMem,
    kinls_mem: &mut KINLsMem,
    v: &NVector,
    z: &mut NVector,
) -> i32 {
    /* call Jacobian-times-vector product routine
       (either user-supplied or internal DQ) */
    let retval = if kinls_mem.jtimesDQ {
        kinLsDQJtimes(kin_mem, &mut kinls_mem.nfeDQ, kinls_mem.jt_func, v, z)
    } else {
        let jt = kinls_mem.jtimes.unwrap();
        let KINMem { kin_uu, kin_user_data, .. } = kin_mem;
        jt(v, z, kin_uu, &mut kinls_mem.new_uu, kin_user_data)
    };
    kinls_mem.njtimes += 1;
    retval
}

/*---------------------------------------------------------------
  kinLsPSetup:

  This routine interfaces between the generic iterative linear
  solvers and the user's psetup routine. It passes to psetup all
  required state information from kin_mem. Its return value
  is the same as that returned by psetup.  In C the generic
  iterative linear solvers guarantee that kinLsPSetup will only
  be called in the case that the user's psetup routine is
  non-NULL; here the None case returns 0 (setup no-op).
  ---------------------------------------------------------------*/
pub fn kinLsPSetup(kin_mem: &mut KINMem, kinls_mem: &mut KINLsMem) -> i32 {
    let KINLsMem { pset, prec_module, npe, .. } = kinls_mem;
    match prec_module {
        /* internal kinsol_bbdpre module (C pset = KINBBDPrecSetup) */
        PrecModule::BBDPre(pdata) => {
            let retval = KINBBDPrecSetup(kin_mem, pdata);
            *npe += 1;
            retval
        }
        _ => {
            if let Some(pset) = *pset {
                /* Call user pset routine to update preconditioner */
                let KINMem { kin_uu, kin_uscale, kin_fval, kin_fscale, kin_user_data, .. } =
                    kin_mem;
                let retval = pset(kin_uu, kin_uscale, kin_fval, kin_fscale, kin_user_data);
                *npe += 1;
                retval
            } else {
                0
            }
        }
    }
}

/* kinLsPSolve — the interface between the generic iterative linear
   solvers and the user's psolve routine — is realized as the psolve
   closure inside kinLsSolveIterative (donor cvode_ls.rs pattern):
   it copies the rhs into z (N_VScale(ONE, r, z)), calls
   psolve(uu, uscale, fval, fscale, z, user_data) — the 'tol' and
   'lr' inputs are not supported by KINSOL — and increments nps. */

/*------------------------------------------------------------------
  kinLsDQJac

  This routine is a wrapper for the Dense and Band implementations
  of the difference quotient Jacobian approximation routines.
  ------------------------------------------------------------------*/
fn kinLsDQJac(kin_mem: &mut KINMem, nfeDQ: &mut i64, jac: &mut SUNMatrix) -> i32 {
    /* Call the matrix-structure-specific DQ approximation routine */
    match jac {
        SUNMatrix::Dense(dm) => kinLsDenseDQJac(kin_mem, nfeDQ, dm),
        SUNMatrix::Band(bm) => kinLsBandDQJac(kin_mem, nfeDQ, bm),
        _ => {
            KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "kinLsDQJac", file!(),
                            "unrecognized matrix type for kinLsDQJac");
            KIN_ILL_INPUT
        }
    }
}

/*------------------------------------------------------------------
  kinLsDenseDQJac

  This routine generates a dense difference quotient approximation
  to the Jacobian of F(u).  u = kin_mem.kin_uu (perturbed in place
  and restored, exactly as the C code perturbs the caller's
  vector); fu = kin_mem.kin_fval; ftemp workspace (tmp1) is
  kin_mem.kin_vtemp1.  The C jthCol N_VSetArrayPointer aliasing of
  tmp2 onto the matrix column becomes a direct write into the
  column slice.

  The increment used in the finite-difference approximation
    J_ij = ( F_i(u+sigma_j * e_j) - F_i(u)  ) / sigma_j
  is
   sigma_j = max{|u_j|, |1/uscale_j|} * sqrt(uround)

  Note: uscale_j = 1/typ(u_j)

  NOTE: Any type of failure of the system function here leads to an
        unrecoverable failure of the Jacobian function and thus of
        the linear solver setup function, stopping KINSOL.
  ------------------------------------------------------------------*/
fn kinLsDenseDQJac(
    kin_mem: &mut KINMem,
    nfeDQ: &mut i64,
    jac: &mut crate::sunmatrix_dense::DenseMatrix,
) -> i32 {
    /* access matrix dimension */
    let n = jac.n;

    let func = kin_mem.kin_func.unwrap();
    let mut retval = 0;

    /* This is the only for loop for 0..N-1 in KINSOL */

    for j in 0..(n as usize) {
        /* Generate the jth col of J(u) */

        /* Save u_j values and scaling */
        let ujsaved = kin_mem.kin_uu.data[j];
        let ujscale = ONE / kin_mem.kin_uscale.data[j];

        /* Compute increment */
        let sign = if ujsaved >= ZERO { ONE } else { -ONE };
        let inc = kin_mem.kin_sqrt_relfunc * SUNMAX(SUNRabs(ujsaved), ujscale) * sign;

        /* Increment u_j, call F(u), and return if error occurs */
        kin_mem.kin_uu.data[j] += inc;

        retval = {
            let KINMem { kin_uu, kin_vtemp1, kin_user_data, .. } = kin_mem;
            func(kin_uu, kin_vtemp1, kin_user_data)
        };
        *nfeDQ += 1;
        if retval != 0 {
            break;
        }

        /* reset u_j */
        kin_mem.kin_uu.data[j] = ujsaved;

        /* Construct difference quotient in jthCol:
           N_VLinearSum(inc_inv, ftemp, -inc_inv, fu, jthCol)
           = inc_inv * (ftemp - fu)  (serial a == -b case) */
        let inc_inv = ONE / inc;
        let col_j = jac.col_mut(j as i64);
        for (i, cj) in col_j.iter_mut().enumerate() {
            *cj = inc_inv * (kin_mem.kin_vtemp1.data[i] - kin_mem.kin_fval.data[i]);
        }
    }

    retval
}

/*------------------------------------------------------------------
  kinLsBandDQJac

  This routine generates a banded difference quotient approximation
  to the Jacobian of F(u).  futemp (tmp1) = kin_vtemp1, utemp
  (tmp2) = kin_vtemp2.

  NOTE: Any type of failure of the system function here leads to an
        unrecoverable failure of the Jacobian function and thus of
        the linear solver setup function, stopping KINSOL.
  ------------------------------------------------------------------*/
fn kinLsBandDQJac(
    kin_mem: &mut KINMem,
    nfeDQ: &mut i64,
    jac: &mut crate::sunmatrix_band::BandMatrix,
) -> i32 {
    /* access matrix dimensions */
    let n = jac.n;
    let mupper = jac.mu;
    let mlower = jac.ml;

    /* Load utemp with u */
    {
        let KINMem { kin_uu, kin_vtemp2, .. } = kin_mem;
        N_VScale(ONE, kin_uu, kin_vtemp2);
    }

    /* Set bandwidth and number of column groups for band differencing */
    let width = mlower + mupper + 1;
    let ngroups = width.min(n);

    let func = kin_mem.kin_func.unwrap();

    for group in 1..=ngroups {
        /* Increment all utemp components in group */
        let mut j = group - 1;
        while j < n {
            let ju = j as usize;
            let inc = kin_mem.kin_sqrt_relfunc
                * SUNMAX(
                    SUNRabs(kin_mem.kin_uu.data[ju]),
                    ONE / SUNRabs(kin_mem.kin_uscale.data[ju]),
                );
            kin_mem.kin_vtemp2.data[ju] += inc;
            j += width;
        }

        /* Evaluate f with incremented u */
        let retval = {
            let KINMem { kin_vtemp2, kin_vtemp1, kin_user_data, .. } = kin_mem;
            func(kin_vtemp2, kin_vtemp1, kin_user_data)
        };
        if retval != 0 {
            return retval;
        }

        /* Restore utemp components, then form and load difference quotients */
        let mut j = group - 1;
        while j < n {
            let ju = j as usize;
            kin_mem.kin_vtemp2.data[ju] = kin_mem.kin_uu.data[ju];
            let inc = kin_mem.kin_sqrt_relfunc
                * SUNMAX(
                    SUNRabs(kin_mem.kin_uu.data[ju]),
                    ONE / SUNRabs(kin_mem.kin_uscale.data[ju]),
                );
            let inc_inv = ONE / inc;
            let i1 = 0.max(j - mupper);
            let i2 = (j + mlower).min(n - 1);
            for i in i1..=i2 {
                let val = inc_inv
                    * (kin_mem.kin_vtemp1.data[i as usize] - kin_mem.kin_fval.data[i as usize]);
                jac.set(i, j, val);
            }
            j += width;
        }
    }

    /* Increment counter nfeDQ */
    *nfeDQ += ngroups;

    0
}

/*------------------------------------------------------------------
  kinLsDQJtimes

  This routine generates the matrix-vector product z = J*v using a
  difference quotient approximation. The approximation is
  J*v = [func(uu + sigma*v) - func(uu)]/sigma. Here sigma is based
  on the dot products (uscale*uu, uscale*v) and
  (uscale*v, uscale*v), the L1Norm(uscale*v), and on sqrt_relfunc
  (the square root of the relative error in the function). Note
  that v in the argument list has already been both preconditioned
  and unscaled.

  u = kin_mem.kin_uu (the C new_u flag argument is unused by this
  DQ routine; the serial N_Vector supplies every required op).

  NOTE: Unlike the DQ Jacobian functions for direct linear solvers
        (which are called from within the lsetup function), this
        function is called from within the lsolve function and thus
        a recovery may still be possible even if the system function
        fails (recoverably).
  ------------------------------------------------------------------*/
fn kinLsDQJtimes(
    kin_mem: &mut KINMem,
    nfeDQ: &mut i64,
    jt_func: Option<KINSysFn>,
    v: &NVector,
    Jv: &mut NVector,
) -> i32 {
    let KINMem {
        kin_uu,
        kin_uscale,
        kin_fval,
        kin_vtemp1,
        kin_vtemp2,
        kin_user_data,
        kin_sqrt_relfunc,
        ..
    } = kin_mem;

    /* scale the vector v and put Du*v into vtemp1 */
    N_VProd(v, kin_uscale, kin_vtemp1);

    /* scale u and put into Jv (used as a temporary storage) */
    N_VProd(kin_uu, kin_uscale, Jv);

    /* compute dot product (Du*u).(Du*v) */
    let sutsv = N_VDotProd(Jv, kin_vtemp1);

    /* compute dot product (Du*v).(Du*v) */
    let vtv = N_VDotProd(kin_vtemp1, kin_vtemp1);

    /* compute differencing factor -- this is from p. 469, Brown and Saad paper */
    let sq1norm = N_VL1Norm(kin_vtemp1);
    let sign = if sutsv >= ZERO { ONE } else { -ONE };
    let sigma = sign * (*kin_sqrt_relfunc) * SUNMAX(SUNRabs(sutsv), sq1norm) / vtv;
    let sigma_inv = ONE / sigma;

    /* compute the u-prime at which to evaluate the function func */
    N_VLinearSum(ONE, kin_uu, sigma, v, kin_vtemp1);

    /* call the system function to calculate func(u+sigma*v) */
    let func = jt_func.unwrap();
    let retval = func(kin_vtemp1, kin_vtemp2, kin_user_data);
    *nfeDQ += 1;
    if retval != 0 {
        return retval;
    }

    /* finish the computation of the difference quotient:
       N_VLinearSum(sigma_inv, vtemp2, -sigma_inv, fval, Jv) */
    N_VLinearSum(sigma_inv, kin_vtemp2, -sigma_inv, kin_fval, Jv);

    0
}

/*------------------------------------------------------------------
  kinLsInitialize performs remaining initializations specific
  to the iterative linear solver interface (and solver itself)

  C: int kinLsInitialize(KINMem kin_mem); the KINLsMem is threaded
  explicitly per the workspace take() convention.
  ------------------------------------------------------------------*/
pub fn kinLsInitialize(kin_mem: &mut KINMem, kinls_mem: &mut KINLsMem) -> i32 {
    /* Test for valid combinations of matrix & Jacobian routines: */
    if kinls_mem.J.is_none() {
        /* If SUNMatrix A is NULL: ensure 'jac' function pointer is NULL */
        kinls_mem.jacDQ = SUNFALSE;
        kinls_mem.jac = None;
    } else if kinls_mem.jacDQ {
        /* If J is non-NULL, and 'jac' is not user-supplied:
           - if A is dense or band, ensure that our DQ approx. is used
           - otherwise => error */
        let ok = matches!(
            kinls_mem.J,
            Some(SUNMatrix::Dense(_)) | Some(SUNMatrix::Band(_))
        );
        if !ok {
            KINProcessError(Some(kin_mem), KINLS_ILL_INPUT, line!(), "kinLsInitialize", file!(),
                            "No Jacobian constructor available for SUNMatrix type");
            kinls_mem.last_flag = KINLS_ILL_INPUT;
            return KINLS_ILL_INPUT;
        }
        kinls_mem.jac = None; /* internal kinLsDQJac */

        /* (required serial vector operations for the kinLsDQJac routine
           — nvlinearsum, nvscale, nvget/setarraypointer — always exist) */
    } else {
        /* If J is non-NULL, and 'jac' is user-supplied: the C J_data
           reset is a no-op here (user_data always comes from KINMem) */
    }

    /* Prohibit Picard iteration with DQ Jacobian approximation or difference-quotient J*v */
    if kin_mem.kin_globalstrategy == KIN_PICARD && kinls_mem.jacDQ && kinls_mem.jtimesDQ {
        KINProcessError(Some(kin_mem), KINLS_ILL_INPUT, line!(), "kinLsInitialize", file!(),
                        MSG_NOL_FAIL);
        return KINLS_ILL_INPUT;
    }

    /* error-checking is complete, begin initializations */

    /* Initialize counters */
    kinLsInitializeCounters(kinls_mem);

    /* Set Jacobian-related fields, based on jtimesDQ (the C jt_data
       assignments are no-ops here) */
    if kinls_mem.jtimesDQ {
        kinls_mem.jtimes = None; /* internal kinLsDQJtimes */
    }

    /* if J is NULL and: NOT preconditioning or do NOT need to setup the
       preconditioner, then set the lsetup function to NULL
       (setup_disabled carries the C `kin_lsetup = NULL` state; the
       kinsol_bbdpre module supplies both pset and psolve in C) */
    if kinls_mem.J.is_none()
        && !matches!(kinls_mem.prec_module, PrecModule::BBDPre(_))
        && (kinls_mem.psolve.is_none() || kinls_mem.pset.is_none())
    {
        kinls_mem.setup_disabled = SUNTRUE;
    }

    /* Set scaling vectors assuming RIGHT preconditioning: every
       workspace iterative solver accepts the fscale/fscale scaling
       vectors (passed as s1 = s2 in kinLsSolveIterative), matching the
       C solvers that implement SUNLinSolSetScalingVectors — so the
       tol_fac = sqrt(N)/||fscale||_L2 tolerance-adjustment branch for
       solvers without scaling support is never taken. */
    kinls_mem.tol_fac = ONE;

    /* Call LS initialize routine, and return result */
    kinls_mem.last_flag = kinls_mem.LS.initialize();
    kinls_mem.last_flag
}

/* kinLsInit — the interface hook dispatched by the main solver
   (kinsol.rs LsModule dispatch name for kinLsInitialize). */
pub fn kinLsInit(kin_mem: &mut KINMem, kinls_mem: &mut KINLsMem) -> i32 {
    kinLsInitialize(kin_mem, kinls_mem)
}

/*------------------------------------------------------------------
  kinLsSetup call the LS setup routine
  ------------------------------------------------------------------*/
pub fn kinLsSetup(kin_mem: &mut KINMem, kinls_mem: &mut KINLsMem) -> i32 {
    /* recompute J if it is non-NULL */
    if kinls_mem.J.is_some() {
        /* Increment nje counter. */
        kinls_mem.nje += 1;

        /* Clear the linear system matrix if necessary */
        if kinls_mem.LS.ls_type() == SUNLINEARSOLVER_DIRECT {
            let retval = kinls_mem.J.as_mut().unwrap().zero();
            if retval != 0 {
                KINProcessError(Some(kin_mem), KINLS_SUNMAT_FAIL, line!(), "kinLsSetup", file!(),
                                MSG_LS_MATZERO_FAILED);
                kinls_mem.last_flag = KINLS_SUNMAT_FAIL;
                return kinls_mem.last_flag;
            }
        }

        /* Call Jacobian routine */
        let retval = if kinls_mem.jacDQ {
            let j = kinls_mem.J.as_mut().unwrap();
            kinLsDQJac(kin_mem, &mut kinls_mem.nfeDQ, j)
        } else {
            let jac = kinls_mem.jac.unwrap();
            let j = kinls_mem.J.as_mut().unwrap();
            let KINMem { kin_uu, kin_fval, kin_user_data, kin_vtemp1, kin_vtemp2, .. } = kin_mem;
            jac(kin_uu, kin_fval, j, kin_user_data, kin_vtemp1, kin_vtemp2)
        };
        if retval != 0 {
            KINProcessError(Some(kin_mem), KINLS_JACFUNC_ERR, line!(), "kinLsSetup", file!(),
                            MSG_LS_JACFUNC_FAILED);
            kinls_mem.last_flag = KINLS_JACFUNC_ERR;
            return kinls_mem.last_flag;
        }
    }

    /* Call LS setup routine -- for direct solvers this factors J; for
       iterative solvers the generic C SUNLinSolSetup calls kinLsPSetup
       (if applicable), which is invoked directly here */
    kinls_mem.last_flag = match &mut kinls_mem.LS {
        LinearSolver::Dense(_) | LinearSolver::Band(_) => {
            let KINLsMem { LS, J, .. } = kinls_mem;
            LS.setup(J.as_mut())
        }
        _ => kinLsPSetup(kin_mem, kinls_mem),
    };

    /* save nni value from most recent lsetup call */
    kin_mem.kin_nnilset = kin_mem.kin_nni;

    kinls_mem.last_flag
}

/*------------------------------------------------------------------
  kinLsSolve interfaces between KINSOL and the generic
  SUNLinearSolver object

  C: int kinLsSolve(KINMem kin_mem, N_Vector xx, N_Vector bb,
                    sunrealtype* sJpnorm, sunrealtype* sFdotJp)
  ------------------------------------------------------------------*/
pub fn kinLsSolve(
    kin_mem: &mut KINMem,
    kinls_mem: &mut KINLsMem,
    xx: &mut NVector,
    bb: &mut NVector,
    sJpnorm: &mut f64,
    sFdotJp: &mut f64,
) -> i32 {
    /* Set linear solver tolerance as input value times scaling factor
       (to account for possible lack of support for left/right scaling
       vectors in SUNLinSol object) */
    let tol = kin_mem.kin_eps * kinls_mem.tol_fac;

    /* Set initial guess x = 0 to LS */
    N_VConst(ZERO, xx);

    /* Set zero initial guess flag */
    kinls_mem.LS.set_zero_guess(SUNTRUE);

    /* set flag required for user-supplied J*v routine */
    kinls_mem.new_uu = SUNTRUE;

    /* Call solver */
    let retval = if !kinls_mem.iterative {
        let KINLsMem { LS, J, .. } = kinls_mem;
        match (LS, J.as_mut()) {
            (LinearSolver::Dense(dls), Some(SUNMatrix::Dense(am))) => dls.solve(am, xx, bb),
            (LinearSolver::Band(bls), Some(SUNMatrix::Band(am))) => bls.solve(am, xx, bb),
            _ => SUN_ERR_ARG_INCOMPATIBLE,
        }
    } else {
        kinLsSolveIterative(kin_mem, kinls_mem, xx, bb, tol)
    };

    /* Retrieve solver statistics (0 for the direct solvers, exactly as
       the C `if (LS->ops->numiters)` guard yields) */
    let nli_inc = kinls_mem.LS.num_iters();

    /* (The C KINPrintInfo PRNT_NLI / PRNT_EPS block is compiled only
       when SUNDIALS_LOGGING_LEVEL >= INFO, which the default build —
       and hence this port — excludes.) */

    /* Increment counters nli and ncfl */
    kinls_mem.nli += nli_inc as i64;
    if retval != SUN_SUCCESS {
        kinls_mem.ncfl += 1;
    }

    /* Interpret solver return value */
    kinls_mem.last_flag = retval;

    if retval != SUN_SUCCESS && retval != SUNLS_RES_REDUCED {
        match retval {
            SUNLS_ATIMES_FAIL_REC | SUNLS_PSOLVE_FAIL_REC => return 1,
            SUN_ERR_ARG_CORRUPT | SUN_ERR_ARG_INCOMPATIBLE | SUN_ERR_MEM_FAIL | SUNLS_GS_FAIL
            | SUNLS_CONV_FAIL | SUNLS_QRFACT_FAIL | SUNLS_LUFACT_FAIL | SUNLS_QRSOL_FAIL => {}
            SUNLS_PACKAGE_FAIL_REC => {
                KINProcessError(Some(kin_mem), SUNLS_PACKAGE_FAIL_REC, line!(), "kinLsSolve",
                                file!(), "Failure in SUNLinSol external package");
            }
            SUN_ERR_EXT_FAIL => {
                KINProcessError(Some(kin_mem), SUN_ERR_EXT_FAIL, line!(), "kinLsSolve", file!(),
                                "Failure in SUNLinSol external package");
            }
            SUNLS_ATIMES_FAIL_UNREC => {
                KINProcessError(Some(kin_mem), SUNLS_ATIMES_FAIL_UNREC, line!(), "kinLsSolve",
                                file!(), MSG_LS_JTIMES_FAILED);
            }
            SUNLS_PSOLVE_FAIL_UNREC => {
                KINProcessError(Some(kin_mem), SUNLS_PSOLVE_FAIL_UNREC, line!(), "kinLsSolve",
                                file!(), MSG_LS_PSOLVE_FAILED);
            }
            _ => {}
        }
        return retval;
    }

    /* SUNLinSolSolve returned SUN_SUCCESS or SUNLS_RES_REDUCED */

    /* Compute auxiliary values for use in the linesearch and in KINForcingTerm.
       These will be subsequently corrected if the step is reduced by constraints
       or the linesearch. */
    if kin_mem.kin_globalstrategy != KIN_FP {
        /* sJpnorm is the norm of the scaled product (scaled by fscale) of the
           current Jacobian matrix J and the step vector p (= solution vector xx) */
        if kin_mem.kin_inexact_ls && kin_mem.kin_etaflag == KIN_ETACHOICE1 {
            let retval = kinLsATimes(kin_mem, kinls_mem, xx, bb);
            if retval > 0 {
                kinls_mem.last_flag = SUNLS_ATIMES_FAIL_REC;
                return 1;
            } else if retval < 0 {
                kinls_mem.last_flag = SUNLS_ATIMES_FAIL_UNREC;
                return -1;
            }
            *sJpnorm = N_VWL2Norm(bb, &kin_mem.kin_fscale);
        }

        /* sFdotJp is the dot product of the scaled f vector and the scaled
           vector J*p, where the scaling uses fscale */
        if (kin_mem.kin_inexact_ls && kin_mem.kin_etaflag == KIN_ETACHOICE1)
            || kin_mem.kin_globalstrategy == KIN_LINESEARCH
        {
            /* N_VProd(bb, fscale, bb) — output aliases the input */
            bb.prod_with(&kin_mem.kin_fscale);
            bb.prod_with(&kin_mem.kin_fscale);
            *sFdotJp = N_VDotProd(&kin_mem.kin_fval, bb);
        }
    }

    0
}

/* Iterative Krylov solve: builds the ATimes (kinLsATimes) and PSolve
   (kinLsPSolve) callbacks over the solver memory.  The two closures
   share kin_mem through a RefCell — the solvers never call them
   re-entrantly (checked at runtime).  The fscale vector is detached
   from KINMem for the duration of the solve so it can be passed as
   the (read-only) scaling vectors s1 = s2 (and to psolve). */
fn kinLsSolveIterative(
    kin_mem: &mut KINMem,
    kinls_mem: &mut KINLsMem,
    xx: &mut NVector,
    bb: &NVector,
    tol: f64,
) -> i32 {
    let KINLsMem {
        LS,
        njtimes,
        nfeDQ,
        nps,
        jtimes,
        jtimesDQ,
        jt_func,
        psolve,
        prec_module,
        new_uu,
        ..
    } = kinls_mem;
    let jtimes = *jtimes;
    let jtimes_dq = *jtimesDQ;
    let jt_func = *jt_func;
    let psolve_fn = *psolve;

    /* In C, PSolve is registered with the LS only when a non-NULL
       psolve was supplied (KINSetPreconditioner) — either by the user
       or by the kinsol_bbdpre module. */
    let has_psolve = psolve_fn.is_some() || matches!(prec_module, PrecModule::BBDPre(_));

    let fscale = std::mem::take(&mut kin_mem.kin_fscale);
    let kinm = RefCell::new(&mut *kin_mem);
    let fscale_ref = &fscale;

    /* kinLsATimes body (see the standalone function above) */
    let mut atimes = |v: &NVector, z: &mut NVector| -> i32 {
        let mut guard = kinm.borrow_mut();
        let kmr: &mut KINMem = &mut *guard;
        let jret = if jtimes_dq {
            kinLsDQJtimes(kmr, nfeDQ, jt_func, v, z)
        } else {
            let jt = jtimes.unwrap();
            let KINMem { kin_uu, kin_user_data, .. } = kmr;
            jt(v, z, kin_uu, new_uu, kin_user_data)
        };
        *njtimes += 1;
        jret
    };

    /* kinLsPSolve body */
    let mut psolve_cb = |r: &NVector, z: &mut NVector, _tol: f64, _lr: i32| -> i32 {
        let mut guard = kinm.borrow_mut();
        let kmr: &mut KINMem = &mut *guard;

        /* copy the rhs into z before the psolve call */
        /* Note: z returns with the solution */
        N_VScale(ONE, r, z);

        /* note: user-supplied preconditioning with KINSOL does not
           support either the 'tol' or 'lr' inputs */
        let ret = match prec_module {
            PrecModule::BBDPre(pdata) => KINBBDPrecSolve(kmr, pdata, z),
            _ => {
                if let Some(ps) = psolve_fn {
                    let KINMem { kin_uu, kin_uscale, kin_fval, kin_user_data, .. } = kmr;
                    ps(kin_uu, kin_uscale, kin_fval, fscale_ref, z, kin_user_data)
                } else {
                    0
                }
            }
        };
        *nps += 1;
        ret
    };

    /* Set scaling vectors assuming RIGHT preconditioning (fscale twice,
       as in the C SUNLinSolSetScalingVectors(LS, fscale, fscale)) */
    let retval = if has_psolve {
        LS.solve(None, xx, bb, tol, &mut atimes, Some(&mut psolve_cb),
                 Some(fscale_ref), Some(fscale_ref))
    } else {
        LS.solve(None, xx, bb, tol, &mut atimes, None, Some(fscale_ref), Some(fscale_ref))
    };

    drop(atimes);
    drop(psolve_cb);
    let kin_mem_back = kinm.into_inner();
    kin_mem_back.kin_fscale = fscale;

    retval
}

/*------------------------------------------------------------------
  kinLsFree frees memory associated with the KINLs system
  solver interface (RAII: dropping the LsModule drops the KINLsMem,
  its LinearSolver, its SUNMatrix J and any preconditioner module)
  ------------------------------------------------------------------*/
pub fn kinLsFree(kin_mem: &mut KINMem) -> i32 {
    /* Return immediately if kin_mem->kin_lmem is NULL */
    if kin_mem.kin_lmem.is_none() {
        return KINLS_SUCCESS;
    }

    /* Nullify SUNMatrix pointer, free preconditioner memory (pfree)
       and the KINLs interface structure — all via drop */
    kin_mem.kin_lmem = LsModule::None;

    KINLS_SUCCESS
}

/*------------------------------------------------------------------
  kinLsInitializeCounters resets counters for the LS interface
  ------------------------------------------------------------------*/
pub fn kinLsInitializeCounters(kinls_mem: &mut KINLsMem) -> i32 {
    kinls_mem.nje = 0;
    kinls_mem.nfeDQ = 0;
    kinls_mem.npe = 0;
    kinls_mem.nli = 0;
    kinls_mem.nps = 0;
    kinls_mem.ncfl = 0;
    kinls_mem.njtimes = 0;
    0
}

/*==================================================================
  Tests
  ==================================================================*/
#[cfg(test)]
mod tests {
    use super::*;
    use crate::sundials_context::SUNContext;
    use crate::sunlinsol_dense::SUNLinSol_Dense;
    use crate::sunmatrix_dense::SUNDenseMatrix;

    /* F(u) = [ u0^2 + u1 - 3 ; u0 - u1^2 ] */
    fn sysfn(uu: &NVector, fval: &mut NVector, _user_data: &mut UserData) -> i32 {
        let u0 = uu.data[0];
        let u1 = uu.data[1];
        fval.data[0] = u0 * u0 + u1 - 3.0;
        fval.data[1] = u0 - u1 * u1;
        0
    }

    fn make_kin_mem(n: usize) -> KINMem {
        let mut kin_mem = KINMem::default();
        kin_mem.kin_func = Some(sysfn);
        kin_mem.kin_uu = NVector::new(n);
        kin_mem.kin_fval = NVector::new(n);
        kin_mem.kin_uscale = NVector::from_slice(&vec![1.0; n]);
        kin_mem.kin_fscale = NVector::from_slice(&vec![1.0; n]);
        kin_mem.kin_vtemp1 = NVector::new(n);
        kin_mem.kin_vtemp2 = NVector::new(n);
        kin_mem
    }

    /* KINSetLinearSolver sets the C defaults (kinsol_ls.c lines
       134-224): inexact_ls flag, jacDQ/jtimesDQ, tol_fac = -1,
       counters zeroed. */
    #[test]
    fn kinsetlinearsolver_defaults() {
        let sunctx = SUNContext::default();
        let mut kin_mem = make_kin_mem(2);

        let a = SUNDenseMatrix(2, 2, &sunctx);
        let ls = SUNLinSol_Dense(&kin_mem.kin_vtemp1, &a, &sunctx);

        let retval = KINSetLinearSolver(&mut kin_mem, ls, Some(a));
        assert_eq!(retval, KINLS_SUCCESS);
        assert!(!kin_mem.kin_inexact_ls); /* direct solver */

        let kinls_mem = kinLs_AccessLMem(&mut kin_mem, "test").unwrap();
        assert!(!kinls_mem.iterative);
        assert!(kinls_mem.matrixbased);
        assert!(kinls_mem.jacDQ);
        assert!(kinls_mem.jtimesDQ);
        assert_eq!(kinls_mem.tol_fac, -1.0);
        assert_eq!(kinls_mem.last_flag, KINLS_SUCCESS);
        assert_eq!(kinls_mem.nje, 0);
        assert_eq!(kinls_mem.nfeDQ, 0);
        assert!(!kinls_mem.setup_disabled);
    }

    /* direct LS without a matrix is rejected (KINLS_ILL_INPUT) */
    #[test]
    fn kinsetlinearsolver_direct_needs_matrix() {
        let sunctx = SUNContext::default();
        let mut kin_mem = make_kin_mem(2);
        let a = SUNDenseMatrix(2, 2, &sunctx);
        let ls = SUNLinSol_Dense(&kin_mem.kin_vtemp1, &a, &sunctx);
        assert_eq!(KINSetLinearSolver(&mut kin_mem, ls, None), KINLS_ILL_INPUT);
        assert!(kin_mem.kin_lmem.is_none());
    }

    /* kinLsSetup drives the dense DQ Jacobian: for the test system the
       exact Jacobian is [[2 u0, 1], [1, -2 u1]]; the difference
       quotient must agree to ~sqrt(uround), and the nje/nfeDQ counters
       and nnilset snapshot follow kinsol_ls.c lines 1160-1199. */
    #[test]
    fn kinlssetup_dense_dq_jacobian() {
        let sunctx = SUNContext::default();
        let mut kin_mem = make_kin_mem(2);
        kin_mem.kin_uu = NVector::from_slice(&[1.5, 0.5]);
        kin_mem.kin_nni = 7;
        {
            let KINMem { kin_uu, kin_fval, kin_user_data, .. } = &mut kin_mem;
            sysfn(kin_uu, kin_fval, kin_user_data);
        }

        let a = SUNDenseMatrix(2, 2, &sunctx);
        let ls = SUNLinSol_Dense(&kin_mem.kin_vtemp1, &a, &sunctx);
        assert_eq!(KINSetLinearSolver(&mut kin_mem, ls, Some(a)), KINLS_SUCCESS);

        /* take the module out for the call, donor dispatch pattern */
        let mut lmem = std::mem::take(&mut kin_mem.kin_lmem);
        let kinls_mem = match &mut lmem {
            LsModule::Ls(m) => m,
            _ => unreachable!(),
        };
        assert_eq!(kinLsInit(&mut kin_mem, kinls_mem), 0);
        assert_eq!(kinls_mem.tol_fac, 1.0);

        /* raw difference-quotient Jacobian (before factorization) */
        {
            let j = kinls_mem.J.as_mut().unwrap();
            assert_eq!(kinLsDQJac(&mut kin_mem, &mut kinls_mem.nfeDQ, j), 0);
            let jac = match &*j {
                SUNMatrix::Dense(dm) => dm,
                _ => unreachable!(),
            };
            let tol = 1.0e-6;
            assert!((jac.get(0, 0) - 3.0).abs() < tol); /* 2*u0   */
            assert!((jac.get(0, 1) - 1.0).abs() < tol);
            assert!((jac.get(1, 0) - 1.0).abs() < tol);
            assert!((jac.get(1, 1) + 1.0).abs() < tol); /* -2*u1  */
        }
        assert_eq!(kinls_mem.nfeDQ, 2); /* one F eval per column */
        kinls_mem.nfeDQ = 0;

        assert_eq!(kinLsSetup(&mut kin_mem, kinls_mem), 0);
        assert_eq!(kinls_mem.nje, 1);
        assert_eq!(kinls_mem.nfeDQ, 2);
        assert_eq!(kin_mem.kin_nnilset, 7);

        /* after kinLsSetup the direct LS has LU-factored J in place
           (C SUNLinSolSetup_Dense does the same): LU([[3,1],[1,-1]])
           = [[3, 1], [1/3, -4/3]] */
        let jac = match kinls_mem.J.as_ref().unwrap() {
            SUNMatrix::Dense(dm) => dm,
            _ => unreachable!(),
        };
        let tol = 1.0e-6;
        assert!((jac.get(1, 0) - 1.0 / 3.0).abs() < tol);
        assert!((jac.get(1, 1) + 4.0 / 3.0).abs() < tol);

        kin_mem.kin_lmem = lmem;
    }

    /* kinLsSolve (direct): solves J x = b after setup; last_flag,
       counters, and the untouched sJpnorm/sFdotJp for KIN_NONE with a
       direct (non-inexact) solver follow kinsol_ls.c lines 1207-1340. */
    #[test]
    fn kinlssolve_direct() {
        let sunctx = SUNContext::default();
        let mut kin_mem = make_kin_mem(2);
        kin_mem.kin_uu = NVector::from_slice(&[1.5, 0.5]);
        {
            let KINMem { kin_uu, kin_fval, kin_user_data, .. } = &mut kin_mem;
            sysfn(kin_uu, kin_fval, kin_user_data);
        }

        let a = SUNDenseMatrix(2, 2, &sunctx);
        let ls = SUNLinSol_Dense(&kin_mem.kin_vtemp1, &a, &sunctx);
        assert_eq!(KINSetLinearSolver(&mut kin_mem, ls, Some(a)), KINLS_SUCCESS);

        let mut lmem = std::mem::take(&mut kin_mem.kin_lmem);
        let kinls_mem = match &mut lmem {
            LsModule::Ls(m) => m,
            _ => unreachable!(),
        };
        assert_eq!(kinLsInit(&mut kin_mem, kinls_mem), 0);
        assert_eq!(kinLsSetup(&mut kin_mem, kinls_mem), 0);

        /* J ~= [[3, 1], [1, -1]]; b = [4, 0] -> x = [1, 1] */
        let mut xx = NVector::new(2);
        let mut bb = NVector::from_slice(&[4.0, 0.0]);
        let mut sJpnorm = 0.0;
        let mut sFdotJp = 0.0;
        let retval = kinLsSolve(&mut kin_mem, kinls_mem, &mut xx, &mut bb,
                                &mut sJpnorm, &mut sFdotJp);
        assert_eq!(retval, 0);
        assert_eq!(kinls_mem.last_flag, SUN_SUCCESS);
        assert_eq!(kinls_mem.nli, 0);
        assert_eq!(kinls_mem.ncfl, 0);
        assert!((xx.data[0] - 1.0).abs() < 1.0e-6);
        assert!((xx.data[1] - 1.0).abs() < 1.0e-6);
        /* KIN_NONE + direct: neither auxiliary product is computed */
        assert_eq!(sJpnorm, 0.0);
        assert_eq!(sFdotJp, 0.0);

        kin_mem.kin_lmem = lmem;

        /* stats getters see the accumulated counters */
        let mut nje = -1;
        assert_eq!(KINGetNumJacEvals(&mut kin_mem, &mut nje), KINLS_SUCCESS);
        assert_eq!(nje, 1);
        let mut flag = -1;
        assert_eq!(KINGetLastLinFlag(&mut kin_mem, &mut flag), KINLS_SUCCESS);
        assert_eq!(flag, SUN_SUCCESS as i64);
    }

    #[test]
    fn return_flag_names() {
        assert_eq!(KINGetLinReturnFlagName(KINLS_SUCCESS as i64), "KINLS_SUCCESS");
        assert_eq!(KINGetLinReturnFlagName(KINLS_JACFUNC_ERR as i64), "KINLS_JACFUNC_ERR");
        assert_eq!(KINGetLinReturnFlagName(1234), "NONE");
    }
}
