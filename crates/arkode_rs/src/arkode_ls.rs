/* -----------------------------------------------------------------
 * Translated from src/arkode/arkode_ls.c (ARKODE 7.7.0).
 * ARKLS: generic linear-solver interface between ARKODE steppers
 * and SUNLinearSolver objects.
 *
 * PART I: the system-matrix half — ARKodeSetLinearSolver, the
 * option/stat routines for the system LS, the internal DQ Jacobian
 * and Jacobian-times-vector approximations, arkLsLinSys, and the
 * linit/lsetup/lsolve/lfree interface installed into the stepper.
 *
 * Not yet ported (lands with mass-matrix support): the ARKLsMassMem
 * half (ARKodeSetMassLinearSolver, arkLsMass*, arkLsMTimes/MPSetup/
 * MPSolve and the mass stat/option routines).  All mass hooks below
 * follow the C code's massmem==NULL paths.
 *
 * Access convention (ARCHITECTURE.md Addendum C.1): the ARKLsMem box
 * is taken out of ark_mem (via the step_getlinmem op) at every entry
 * point, passed to an _inner worker together with ark_mem, and put
 * back before returning.  C's arkLs_AccessARKODELMem (void* entry)
 * and arkLs_AccessLMem collapse to the single helper below.
 * -----------------------------------------------------------------*/
use crate::arkode::{arkAllocVec, arkFreeVec};
use crate::arkode_impl::*;
use crate::arkode_ls_impl::*;
use crate::nvector_serial::*;
use crate::sundials_linearsolver::*;
use crate::sundials_math::*;
use crate::sundials_errors::*;
use crate::sundials_matrix::*;
use crate::sundials_types::*;
use std::cell::RefCell;

const MIN_INC_MULT: f64 = 1000.0;
const PT25: f64 = 0.25;

/*===============================================================
  Access helper
  ===============================================================*/

/* C: arkLs_AccessLMem / arkLs_AccessARKODELMem.  Takes the ARKLsMem
   box out through the step_getlinmem op; callers put it back by
   writing ark_mem.lmem. */
pub(crate) fn arkLs_AccessLMem(
    ark_mem: &mut ARKodeMem,
    fname: &str,
) -> Result<Box<ARKLsMem>, i32> {
    let taken = match ark_mem.step_getlinmem {
        Some(get) => get(ark_mem),
        None => None,
    };
    match taken {
        Some(l) => Ok(l),
        None => {
            arkProcessError(
                Some(ark_mem),
                ARKLS_LMEM_NULL,
                line!(),
                fname,
                file!(),
                MSG_LS_LMEM_NULL,
            );
            Err(ARKLS_LMEM_NULL)
        }
    }
}

/* rwt aliasing helper: rwt_is_ewt == SUNTRUE means the C rwt pointer
   aliases ewt (the Rust rwt is left unallocated; Addendum C.1). */
fn ark_rwt(ark_mem: &ARKodeMem) -> &NVector {
    if ark_mem.rwt_is_ewt {
        &ark_mem.ewt
    } else {
        &ark_mem.rwt
    }
}

/*===============================================================
  ARKLS Exported functions -- Required
  ===============================================================*/

/*---------------------------------------------------------------
  ARKodeSetLinearSolver specifies the linear solver.
  ---------------------------------------------------------------*/
pub fn ARKodeSetLinearSolver(
    ark_mem: &mut ARKodeMem,
    LS: LinearSolver,
    A: Option<SUNMatrix>,
) -> i32 {
    /* Guard against use for time steppers that do not need an algebraic solver */
    if !ark_mem.step_supports_implicit {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!(),
            "ARKodeSetLinearSolver",
            file!(),
            "time-stepping module does not require an algebraic solver",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* (C checks LS/vector ops tables for required operations; the Rust
       LinearSolver enum and serial NVector always provide them.) */

    /* Retrieve the LS type */
    let LSType = LS.ls_type();

    /* Set flags based on LS type */
    let iterative = LSType != SUNLINEARSOLVER_DIRECT;
    let matrixbased =
        LSType != SUNLINEARSOLVER_ITERATIVE && LSType != SUNLINEARSOLVER_MATRIX_EMBEDDED;

    /* Ensure that A is NULL when LS is matrix-embedded */
    if LSType == SUNLINEARSOLVER_MATRIX_EMBEDDED && A.is_some() {
        arkProcessError(
            Some(ark_mem),
            ARKLS_ILL_INPUT,
            line!(),
            "ARKodeSetLinearSolver",
            file!(),
            "Incompatible inputs: matrix-embedded LS requires NULL matrix",
        );
        return ARKLS_ILL_INPUT;
    }

    /* Check for compatible LS type, matrix and "atimes" support */
    if iterative {
        /* (iterative Rust solvers always support the ATimes routine) */
        if matrixbased && A.is_none() {
            arkProcessError(
                Some(ark_mem),
                ARKLS_ILL_INPUT,
                line!(),
                "ARKodeSetLinearSolver",
                file!(),
                "Incompatible inputs: matrix-iterative LS requires non-NULL matrix",
            );
            return ARKLS_ILL_INPUT;
        }
    } else if A.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARKLS_ILL_INPUT,
            line!(),
            "ARKodeSetLinearSolver",
            file!(),
            "Incompatible inputs: direct LS requires non-NULL matrix",
        );
        return ARKLS_ILL_INPUT;
    }

    /* Test whether time stepper module is supplied, with required routines */
    if ark_mem.step_attachlinsol.is_none()
        || ark_mem.step_getlinmem.is_none()
        || ark_mem.step_getimplicitrhs.is_none()
        || ark_mem.step_getgammas.is_none()
    {
        arkProcessError(
            Some(ark_mem),
            ARKLS_ILL_INPUT,
            line!(),
            "ARKodeSetLinearSolver",
            file!(),
            "Missing time step module or associated routines",
        );
        return ARKLS_ILL_INPUT;
    }

    /* Set default Jt_f from the stepper's implicit RHS */
    let getrhs = ark_mem.step_getimplicitrhs.unwrap();
    let Jt_f = getrhs(ark_mem);
    if Jt_f.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARKLS_ILL_INPUT,
            line!(),
            "ARKodeSetLinearSolver",
            file!(),
            "Time step module is missing implicit RHS fcn",
        );
        return ARKLS_ILL_INPUT;
    }

    let has_matrix = A.is_some();

    /* Allocate memory for ARKLsMemRec, set defaults (memset-0 + the
       explicit C initializations) */
    let mut arkls_mem = Box::new(ARKLsMem {
        /* Linear solver type information */
        iterative,
        matrixbased,

        /* Set defaults for Jacobian-related fields */
        jacDQ: has_matrix,
        jac: None, /* None + jacDQ => internal arkLsDQJac */
        jbad: SUNTRUE,

        scalesol: SUNFALSE,

        eplifac: ARKLS_EPLIN,
        nrmfac: ZERO,

        LS,
        A,
        savedJ: None, /* allocated in arkLsInitialize */
        ytemp: NVector::new(0),
        x: NVector::new(0),

        msbj: ARKLS_MSBJ,
        tcur: ZERO,
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

        /* Jacobian-times-vector: internal DQ by default */
        jtimesDQ: SUNTRUE,
        jtsetup: None,
        jtimes: None, /* None + jtimesDQ => internal arkLsDQJtimes */
        Jt_f,

        user_linsys: SUNFALSE,
        linsys: None, /* None + !user_linsys => internal arkLsLinSys */

        last_flag: ARKLS_SUCCESS,
    });

    /* Initialize counters */
    arkLsInitializeCounters(&mut arkls_mem);

    /* (C attaches arkLsATimes / NULL preconditioner hooks to the LS
       object here; the Rust iterative solvers receive those callbacks
       at solve time.) */

    /* Allocate memory for ytemp and x (arkAllocVec cannot fail here) */
    let tmpl_len = ark_mem.tempv1.data.len();
    arkAllocVec(ark_mem, tmpl_len, &mut arkls_mem.ytemp);
    arkAllocVec(ark_mem, tmpl_len, &mut arkls_mem.x);

    /* For iterative LS, compute default norm conversion factor */
    if iterative {
        arkls_mem.nrmfac = SUNRsqrt(N_VGetLength(&arkls_mem.ytemp) as f64);
    }

    /* For matrix-based LS, enable solution scaling */
    arkls_mem.scalesol = matrixbased;

    /* Attach ARKLs interface to time stepper module */
    let attach = ark_mem.step_attachlinsol.unwrap();
    let retval = attach(
        ark_mem,
        Some(arkLsInitialize),
        Some(arkLsSetup),
        Some(arkLsSolve),
        Some(arkLsFree),
        LSType,
        arkls_mem,
    );
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            retval,
            line!(),
            "ARKodeSetLinearSolver",
            file!(),
            "Failed to attach to time stepper module",
        );
        return retval;
    }

    ARKLS_SUCCESS
}

/*===============================================================
  Optional Set routines
  ===============================================================*/

/* Shared entry guard for the ARKLS option routines:
   supports_implicit check + lmem take. */
macro_rules! arkls_option_fn {
    ($ark_mem:ident, $fname:literal) => {{
        if !$ark_mem.step_supports_implicit {
            arkProcessError(
                Some($ark_mem),
                ARK_STEPPER_UNSUPPORTED,
                line!(),
                $fname,
                file!(),
                "time-stepping module does not require an algebraic solver",
            );
            return ARK_STEPPER_UNSUPPORTED;
        }
        match arkLs_AccessLMem($ark_mem, $fname) {
            Ok(l) => l,
            Err(e) => return e,
        }
    }};
}

/*---------------------------------------------------------------
  ARKodeSetJacFn specifies the Jacobian function.
  ---------------------------------------------------------------*/
pub fn ARKodeSetJacFn(ark_mem: &mut ARKodeMem, jac: Option<ARKLsJacFn>) -> i32 {
    let mut arkls_mem = arkls_option_fn!(ark_mem, "ARKodeSetJacFn");

    /* return with failure if jac cannot be used */
    if jac.is_some() && arkls_mem.A.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARKLS_ILL_INPUT,
            line!(),
            "ARKodeSetJacFn",
            file!(),
            "Jacobian routine cannot be supplied for NULL SUNMatrix",
        );
        ark_mem.lmem = Some(arkls_mem);
        return ARKLS_ILL_INPUT;
    }

    /* set the Jacobian routine pointer, and update relevant flags */
    if jac.is_some() {
        arkls_mem.jacDQ = SUNFALSE;
        arkls_mem.jac = jac;
    } else {
        arkls_mem.jacDQ = SUNTRUE;
        arkls_mem.jac = None; /* internal arkLsDQJac */
    }

    /* ensure the internal linear system function is used */
    arkls_mem.user_linsys = SUNFALSE;
    arkls_mem.linsys = None; /* internal arkLsLinSys */

    ark_mem.lmem = Some(arkls_mem);
    ARKLS_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetEpsLin specifies the nonlinear -> linear tolerance
  scale factor.
  ---------------------------------------------------------------*/
pub fn ARKodeSetEpsLin(ark_mem: &mut ARKodeMem, eplifac: f64) -> i32 {
    let mut arkls_mem = arkls_option_fn!(ark_mem, "ARKodeSetEpsLin");

    /* store input and return */
    arkls_mem.eplifac = if eplifac <= ZERO { ARKLS_EPLIN } else { eplifac };

    ark_mem.lmem = Some(arkls_mem);
    ARKLS_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetLSNormFactor sets or computes the factor to use when
  converting from the integrator tolerance (WRMS norm) to the
  linear solver tolerance (L2 norm).
  ---------------------------------------------------------------*/
pub fn ARKodeSetLSNormFactor(ark_mem: &mut ARKodeMem, nrmfac: f64) -> i32 {
    let mut arkls_mem = arkls_option_fn!(ark_mem, "ARKodeSetLSNormFactor");

    /* store input and return */
    if nrmfac > ZERO {
        /* set user-provided factor */
        arkls_mem.nrmfac = nrmfac;
    } else if nrmfac < ZERO {
        /* compute factor for WRMS norm with dot product */
        N_VConst(ONE, &mut ark_mem.tempv1);
        arkls_mem.nrmfac = SUNRsqrt(N_VDotProd(&ark_mem.tempv1, &ark_mem.tempv1));
    } else {
        /* compute default factor for WRMS norm from vector length */
        arkls_mem.nrmfac = SUNRsqrt(N_VGetLength(&ark_mem.tempv1) as f64);
    }

    ark_mem.lmem = Some(arkls_mem);
    ARKLS_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetJacEvalFrequency specifies the frequency for
  recomputing the Jacobian matrix and/or preconditioner.
  ---------------------------------------------------------------*/
pub fn ARKodeSetJacEvalFrequency(ark_mem: &mut ARKodeMem, msbj: i64) -> i32 {
    let mut arkls_mem = arkls_option_fn!(ark_mem, "ARKodeSetJacEvalFrequency");

    /* store input and return */
    arkls_mem.msbj = if msbj <= 0 { ARKLS_MSBJ } else { msbj };

    ark_mem.lmem = Some(arkls_mem);
    ARKLS_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetLinearSolutionScaling enables or disables scaling the
  linear solver solution to account for changes in gamma.
  ---------------------------------------------------------------*/
pub fn ARKodeSetLinearSolutionScaling(ark_mem: &mut ARKodeMem, onoff: bool) -> i32 {
    let mut arkls_mem = arkls_option_fn!(ark_mem, "ARKodeSetLinearSolutionScaling");

    /* check for valid solver type */
    if !arkls_mem.matrixbased {
        ark_mem.lmem = Some(arkls_mem);
        return ARKLS_ILL_INPUT;
    }

    /* set solution scaling flag */
    arkls_mem.scalesol = onoff;

    ark_mem.lmem = Some(arkls_mem);
    ARKLS_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetPreconditioner specifies the user-supplied
  preconditioner setup and solve routines.
  ---------------------------------------------------------------*/
pub fn ARKodeSetPreconditioner(
    ark_mem: &mut ARKodeMem,
    psetup: Option<ARKLsPrecSetupFn>,
    psolve: Option<ARKLsPrecSolveFn>,
) -> i32 {
    let mut arkls_mem = arkls_option_fn!(ark_mem, "ARKodeSetPreconditioner");

    /* issue error if LS object does not allow user-supplied preconditioning */
    if !arkls_mem.iterative {
        arkProcessError(
            Some(ark_mem),
            ARKLS_ILL_INPUT,
            line!(),
            "ARKodeSetPreconditioner",
            file!(),
            "SUNLinearSolver object does not support user-supplied preconditioning",
        );
        ark_mem.lmem = Some(arkls_mem);
        return ARKLS_ILL_INPUT;
    }

    /* store function pointers for user-supplied routines */
    arkls_mem.pset = psetup;
    arkls_mem.psolve = psolve;

    /* (C notifies the LS object here via SUNLinSolSetPreconditioner;
       the Rust iterative solvers receive the psolve closure at solve
       time.) */

    ark_mem.lmem = Some(arkls_mem);
    ARKLS_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetJacTimes specifies the user-supplied Jacobian-vector
  product setup and multiply routines.
  ---------------------------------------------------------------*/
pub fn ARKodeSetJacTimes(
    ark_mem: &mut ARKodeMem,
    jtsetup: Option<ARKLsJacTimesSetupFn>,
    jtimes: Option<ARKLsJacTimesVecFn>,
) -> i32 {
    let mut arkls_mem = arkls_option_fn!(ark_mem, "ARKodeSetJacTimes");

    /* issue error if LS object does not allow user-supplied ATimes */
    if !arkls_mem.iterative {
        arkProcessError(
            Some(ark_mem),
            ARKLS_ILL_INPUT,
            line!(),
            "ARKodeSetJacTimes",
            file!(),
            "SUNLinearSolver object does not support user-supplied ATimes routine",
        );
        ark_mem.lmem = Some(arkls_mem);
        return ARKLS_ILL_INPUT;
    }

    /* store function pointers for user-supplied routines in ARKLs
       interface (NULL jtimes implies use of DQ default) */
    if jtimes.is_some() {
        arkls_mem.jtimesDQ = SUNFALSE;
        arkls_mem.jtsetup = jtsetup;
        arkls_mem.jtimes = jtimes;
    } else {
        arkls_mem.jtimesDQ = SUNTRUE;
        arkls_mem.jtsetup = None;
        arkls_mem.jtimes = None; /* internal arkLsDQJtimes */
        let getrhs = ark_mem.step_getimplicitrhs.unwrap();
        arkls_mem.Jt_f = getrhs(ark_mem);

        if arkls_mem.Jt_f.is_none() {
            arkProcessError(
                Some(ark_mem),
                ARKLS_ILL_INPUT,
                line!(),
                "ARKodeSetJacTimes",
                file!(),
                "Time step module is missing implicit RHS fcn",
            );
            ark_mem.lmem = Some(arkls_mem);
            return ARKLS_ILL_INPUT;
        }
    }

    ark_mem.lmem = Some(arkls_mem);
    ARKLS_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetJacTimesRhsFn specifies an alternative user-supplied
  ODE right-hand side function to use in the internal finite
  difference Jacobian-vector product.
  ---------------------------------------------------------------*/
pub fn ARKodeSetJacTimesRhsFn(ark_mem: &mut ARKodeMem, jtimesRhsFn: Option<ARKRhsFn>) -> i32 {
    let mut arkls_mem = arkls_option_fn!(ark_mem, "ARKodeSetJacTimesRhsFn");

    /* check if using internal finite difference approximation */
    if !arkls_mem.jtimesDQ {
        arkProcessError(
            Some(ark_mem),
            ARKLS_ILL_INPUT,
            line!(),
            "ARKodeSetJacTimesRhsFn",
            file!(),
            "Internal finite-difference Jacobian-vector product is disabled.",
        );
        ark_mem.lmem = Some(arkls_mem);
        return ARKLS_ILL_INPUT;
    }

    /* store function pointers for RHS function (NULL implies use ODE RHS) */
    if jtimesRhsFn.is_some() {
        arkls_mem.Jt_f = jtimesRhsFn;
    } else {
        let getrhs = ark_mem.step_getimplicitrhs.unwrap();
        arkls_mem.Jt_f = getrhs(ark_mem);

        if arkls_mem.Jt_f.is_none() {
            arkProcessError(
                Some(ark_mem),
                ARKLS_ILL_INPUT,
                line!(),
                "ARKodeSetJacTimesRhsFn",
                file!(),
                "Time step module is missing implicit RHS fcn",
            );
            ark_mem.lmem = Some(arkls_mem);
            return ARKLS_ILL_INPUT;
        }
    }

    ark_mem.lmem = Some(arkls_mem);
    ARKLS_SUCCESS
}

/* ARKodeSetLinSysFn specifies the linear system setup function. */
pub fn ARKodeSetLinSysFn(ark_mem: &mut ARKodeMem, linsys: Option<ARKLsLinSysFn>) -> i32 {
    let mut arkls_mem = arkls_option_fn!(ark_mem, "ARKodeSetLinSysFn");

    /* return with failure if linsys cannot be used */
    if linsys.is_some() && arkls_mem.A.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARKLS_ILL_INPUT,
            line!(),
            "ARKodeSetLinSysFn",
            file!(),
            "Linear system setup routine cannot be supplied for NULL SUNMatrix",
        );
        ark_mem.lmem = Some(arkls_mem);
        return ARKLS_ILL_INPUT;
    }

    /* set the linear system routine pointer, and update relevant flags */
    if linsys.is_some() {
        arkls_mem.user_linsys = SUNTRUE;
        arkls_mem.linsys = linsys;
    } else {
        arkls_mem.user_linsys = SUNFALSE;
        arkls_mem.linsys = None; /* internal arkLsLinSys */
    }

    ark_mem.lmem = Some(arkls_mem);
    ARKLS_SUCCESS
}

/*===============================================================
  Optional Get routines
  ===============================================================*/

/* ARKodeGetJac returns a copy of the internally stored Jacobian
   (C returns the raw savedJ pointer). */
pub fn ARKodeGetJac(ark_mem: &mut ARKodeMem, J: &mut Option<SUNMatrix>) -> i32 {
    /* Return NULL for incompatible steppers */
    if !ark_mem.step_supports_implicit {
        *J = None;
        return ARK_SUCCESS;
    }

    /* access ARKLsMem structure */
    let arkls_mem = match arkLs_AccessLMem(ark_mem, "ARKodeGetJac") {
        Ok(l) => l,
        Err(e) => return e,
    };

    /* set output and return */
    *J = arkls_mem.savedJ.as_ref().map(SUNMatClone_Copy);
    ark_mem.lmem = Some(arkls_mem);
    ARKLS_SUCCESS
}

/* helper for ARKodeGetJac: clone-with-contents */
fn SUNMatClone_Copy(m: &SUNMatrix) -> SUNMatrix {
    let mut out = SUNMatClone(m);
    SUNMatCopy(m, &mut out);
    out
}

pub fn ARKodeGetJacTime(ark_mem: &mut ARKodeMem, t_J: &mut f64) -> i32 {
    /* Return an error for incompatible steppers */
    if !ark_mem.step_supports_implicit {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!(),
            "ARKodeGetJacTime",
            file!(),
            "time-stepping module does not require an algebraic solver",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    let arkls_mem = match arkLs_AccessLMem(ark_mem, "ARKodeGetJacTime") {
        Ok(l) => l,
        Err(e) => return e,
    };

    *t_J = arkls_mem.tnlj;
    ark_mem.lmem = Some(arkls_mem);
    ARKLS_SUCCESS
}

pub fn ARKodeGetJacNumSteps(ark_mem: &mut ARKodeMem, nst_J: &mut i64) -> i32 {
    /* Return 0 for incompatible steppers */
    if !ark_mem.step_supports_implicit {
        *nst_J = 0;
        return ARK_SUCCESS;
    }

    let arkls_mem = match arkLs_AccessLMem(ark_mem, "ARKodeGetJacNumSteps") {
        Ok(l) => l,
        Err(e) => return e,
    };

    *nst_J = arkls_mem.nstlj;
    ark_mem.lmem = Some(arkls_mem);
    ARKLS_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeGetLinWorkSpace returns the length of workspace allocated
  for the ARKLS linear solver interface.
  ---------------------------------------------------------------*/
pub fn ARKodeGetLinWorkSpace(ark_mem: &mut ARKodeMem, lenrw: &mut i64, leniw: &mut i64) -> i32 {
    /* Return 0 for incompatible steppers */
    if !ark_mem.step_supports_implicit {
        *lenrw = 0;
        *leniw = 0;
        return ARK_SUCCESS;
    }

    let arkls_mem = match arkLs_AccessLMem(ark_mem, "ARKodeGetLinWorkSpace") {
        Ok(l) => l,
        Err(e) => return e,
    };

    /* start with fixed sizes plus vector/matrix pointers */
    *lenrw = 3;
    *leniw = 30;

    /* add NVector sizes */
    let lrw1 = arkls_mem.x.data.len() as i64;
    let liw1 = 1i64; /* N_VSpace_Serial: (n, 1) */
    *lenrw += 2 * lrw1;
    *leniw += 2 * liw1;

    /* add SUNMatrix size (only account for the one owned by Ls interface) */
    if let Some(saved) = &arkls_mem.savedJ {
        let mut lrw = 0i64;
        let mut liw = 0i64;
        if SUNMatSpace(saved, &mut lrw, &mut liw) == 0 {
            *lenrw += lrw;
            *leniw += liw;
        }
    }

    /* add LS sizes */
    let (lrw, liw) = arkls_mem.LS.space();
    *lenrw += lrw;
    *leniw += liw;

    ark_mem.lmem = Some(arkls_mem);
    ARKLS_SUCCESS
}

/* Shared boilerplate for the scalar statistics getters: guard +
   access + copy-out. */
macro_rules! arkls_stat_get {
    ($name:ident, $out_ty:ty, $field:ident) => {
        pub fn $name(ark_mem: &mut ARKodeMem, out: &mut $out_ty) -> i32 {
            /* Return 0 for incompatible steppers */
            if !ark_mem.step_supports_implicit {
                *out = 0 as $out_ty;
                return ARK_SUCCESS;
            }

            let arkls_mem = match arkLs_AccessLMem(ark_mem, stringify!($name)) {
                Ok(l) => l,
                Err(e) => return e,
            };

            *out = arkls_mem.$field;
            ark_mem.lmem = Some(arkls_mem);
            ARKLS_SUCCESS
        }
    };
}

arkls_stat_get!(ARKodeGetNumJacEvals, i64, nje);
arkls_stat_get!(ARKodeGetNumLinRhsEvals, i64, nfeDQ);
arkls_stat_get!(ARKodeGetNumPrecEvals, i64, npe);
arkls_stat_get!(ARKodeGetNumPrecSolves, i64, nps);
arkls_stat_get!(ARKodeGetNumLinIters, i64, nli);
arkls_stat_get!(ARKodeGetNumLinConvFails, i64, ncfl);
arkls_stat_get!(ARKodeGetNumJTSetupEvals, i64, njtsetup);
arkls_stat_get!(ARKodeGetNumJtimesEvals, i64, njtimes);

pub fn ARKodeGetLastLinFlag(ark_mem: &mut ARKodeMem, flag: &mut i64) -> i32 {
    /* Return 0 for incompatible steppers */
    if !ark_mem.step_supports_implicit {
        *flag = 0;
        return ARK_SUCCESS;
    }

    let arkls_mem = match arkLs_AccessLMem(ark_mem, "ARKodeGetLastLinFlag") {
        Ok(l) => l,
        Err(e) => return e,
    };

    *flag = arkls_mem.last_flag as i64;
    ark_mem.lmem = Some(arkls_mem);
    ARKLS_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeGetLinReturnFlagName translates from the integer error
  code returned by an ARKLs routine to the corresponding string
  equivalent for that flag
  ---------------------------------------------------------------*/
pub fn ARKodeGetLinReturnFlagName(flag: i64) -> String {
    match flag as i32 {
        ARKLS_SUCCESS => "ARKLS_SUCCESS",
        ARKLS_MEM_NULL => "ARKLS_MEM_NULL",
        ARKLS_LMEM_NULL => "ARKLS_LMEM_NULL",
        ARKLS_ILL_INPUT => "ARKLS_ILL_INPUT",
        ARKLS_MEM_FAIL => "ARKLS_MEM_FAIL",
        ARKLS_MASSMEM_NULL => "ARKLS_MASSMEM_NULL",
        ARKLS_JACFUNC_UNRECVR => "ARKLS_JACFUNC_UNRECVR",
        ARKLS_JACFUNC_RECVR => "ARKLS_JACFUNC_RECVR",
        ARKLS_MASSFUNC_UNRECVR => "ARKLS_MASSFUNC_UNRECVR",
        ARKLS_MASSFUNC_RECVR => "ARKLS_MASSFUNC_RECVR",
        ARKLS_SUNMAT_FAIL => "ARKLS_SUNMAT_FAIL",
        ARKLS_SUNLS_FAIL => "ARKLS_SUNLS_FAIL",
        _ => "NONE",
    }
    .to_string()
}

/*===============================================================
  ARKLS Private functions
  ===============================================================*/

/*---------------------------------------------------------------
  arkLSSetUserData sets user_data pointers in arkLS interface.

  In C this re-points J_data/Jt_data/P_data at the new user_data
  when the corresponding defaults are in use; in Rust those data
  pointers collapsed onto ark_mem.user_data itself, so there is
  nothing to update.
  ---------------------------------------------------------------*/
pub fn arkLSSetUserData(_ark_mem: &mut ARKodeMem) -> i32 {
    ARKLS_SUCCESS
}

/*---------------------------------------------------------------
  arkLsDQJac: internal wrapper choosing the dense or band DQ
  Jacobian approximation (the jac == None + jacDQ dispatch target).
  ---------------------------------------------------------------*/
fn arkLsDQJac(
    ark_mem: &mut ARKodeMem,
    arkls_mem: &mut ARKLsMem,
    fi: ARKRhsFn,
    t: f64,
    y: &mut NVector,
    fy: &NVector,
    jac: &mut SUNMatrix,
    tmp1: &mut NVector,
    tmp2: &mut NVector,
) -> i32 {
    /* Call the matrix-structure-specific DQ approximation routine */
    match jac {
        SUNMatrix::Dense(dj) => arkLsDenseDQJac(t, y, fy, dj, ark_mem, arkls_mem, fi, tmp1),
        SUNMatrix::Band(bj) => arkLsBandDQJac(t, y, fy, bj, ark_mem, arkls_mem, fi, tmp1, tmp2),
        _ => {
            arkProcessError(
                Some(ark_mem),
                ARKLS_ILL_INPUT,
                line!(),
                "arkLsDQJac",
                file!(),
                "arkLsDQJac not implemented for this SUNMatrix type!",
            );
            ARKLS_ILL_INPUT
        }
    }
}

/*---------------------------------------------------------------
  arkLsDenseDQJac:

  This routine generates a dense difference quotient approximation
  to the Jacobian of f(t,y).  y is perturbed in place and restored,
  exactly as the C code perturbs the caller's vector; ftemp = tmp1.
  ---------------------------------------------------------------*/
fn arkLsDenseDQJac(
    t: f64,
    y: &mut NVector,
    fy: &NVector,
    jac: &mut crate::sunmatrix_dense::DenseMatrix,
    ark_mem: &mut ARKodeMem,
    arkls_mem: &mut ARKLsMem,
    fi: ARKRhsFn,
    tmp1: &mut NVector,
) -> i32 {
    /* access matrix dimension */
    let n = jac.n;

    /* Rename work vector for readability */
    let ftemp = tmp1;

    /* Set minimum increment based on uround and norm of f */
    let srur = SUNRsqrt(ark_mem.uround);
    let fnorm = N_VWrmsNorm(fy, ark_rwt(ark_mem));
    let minInc = if fnorm != ZERO {
        MIN_INC_MULT * SUNRabs(ark_mem.h) * ark_mem.uround * n as f64 * fnorm
    } else {
        ONE
    };

    let mut retval = 0;

    for j in 0..(n as usize) {
        /* Generate the jth col of J(tn,y) */
        let yjsaved = y.data[j];
        let mut inc = SUNMAX(srur * SUNRabs(yjsaved), minInc / ark_mem.ewt.data[j]);

        /* Adjust sign(inc) if y_j has an inequality constraint. */
        if let Some(constraints) = &ark_mem.constraints {
            let conj = constraints.data[j];
            if SUNRabs(conj) == ONE {
                if (yjsaved + inc) * conj < ZERO {
                    inc = -inc;
                }
            } else if SUNRabs(conj) == TWO && (yjsaved + inc) * conj <= ZERO {
                inc = -inc;
            }
        }

        y.data[j] += inc;

        /* call the user-supplied pre-RHS function (if supplied), then call RHS */
        if let Some(pre_rhs) = ark_mem.PreRhsFn {
            retval = pre_rhs(t, y, &mut ark_mem.user_data);
            if retval != 0 {
                return ARK_PRERHSFN_FAIL;
            }
        }
        retval = fi(t, y, ftemp, &mut ark_mem.user_data);
        arkls_mem.nfeDQ += 1;
        if retval != 0 {
            break;
        }

        y.data[j] = yjsaved;

        let inc_inv = ONE / inc;
        let col_j = jac.col_mut(j as i64);
        for (i, cj) in col_j.iter_mut().enumerate() {
            /* N_VLinearSum(inc_inv, ftemp, -inc_inv, fy, jthCol) */
            *cj = inc_inv * ftemp.data[i] + (-inc_inv) * fy.data[i];
        }
    }

    retval
}

/*---------------------------------------------------------------
  arkLsBandDQJac:

  This routine generates a banded difference quotient approximation
  to the Jacobian of f(t,y).  ftemp = tmp1, ytemp = tmp2.
  ---------------------------------------------------------------*/
fn arkLsBandDQJac(
    t: f64,
    y: &NVector,
    fy: &NVector,
    jac: &mut crate::sunmatrix_band::BandMatrix,
    ark_mem: &mut ARKodeMem,
    arkls_mem: &mut ARKLsMem,
    fi: ARKRhsFn,
    tmp1: &mut NVector,
    tmp2: &mut NVector,
) -> i32 {
    /* access matrix dimensions */
    let n = jac.n;
    let mupper = jac.mu;
    let mlower = jac.ml;

    /* Rename work vectors for use as temporary values of y and f */
    let ftemp = tmp1;
    let ytemp = tmp2;

    /* Load ytemp with y = predicted y vector */
    ytemp.data.copy_from_slice(&y.data);

    /* Set minimum increment based on uround and norm of f */
    let srur = SUNRsqrt(ark_mem.uround);
    let fnorm = N_VWrmsNorm(fy, ark_rwt(ark_mem));
    let minInc = if fnorm != ZERO {
        MIN_INC_MULT * SUNRabs(ark_mem.h) * ark_mem.uround * n as f64 * fnorm
    } else {
        ONE
    };

    /* Set bandwidth and number of column groups for band differencing */
    let width = mlower + mupper + 1;
    let ngroups = width.min(n);

    let mut retval = 0;

    /* Loop over column groups. */
    for group in 1..=ngroups {
        /* Increment all y_j in group */
        let mut j = group - 1;
        while j < n {
            let ju = j as usize;
            let mut inc = SUNMAX(
                srur * SUNRabs(y.data[ju]),
                minInc / ark_mem.ewt.data[ju],
            );

            /* Adjust sign(inc) if yj has an inequality constraint. */
            if let Some(constraints) = &ark_mem.constraints {
                let conj = constraints.data[ju];
                if SUNRabs(conj) == ONE {
                    if (ytemp.data[ju] + inc) * conj < ZERO {
                        inc = -inc;
                    }
                } else if SUNRabs(conj) == TWO && (ytemp.data[ju] + inc) * conj <= ZERO {
                    inc = -inc;
                }
            }

            ytemp.data[ju] += inc;
            j += width;
        }

        /* call the user-supplied pre-RHS function (if supplied), then call RHS */
        if let Some(pre_rhs) = ark_mem.PreRhsFn {
            retval = pre_rhs(t, ytemp, &mut ark_mem.user_data);
            if retval != 0 {
                return ARK_PRERHSFN_FAIL;
            }
        }
        retval = fi(t, ytemp, ftemp, &mut ark_mem.user_data);
        arkls_mem.nfeDQ += 1;
        if retval != 0 {
            break;
        }

        /* Restore ytemp, then form and load difference quotients */
        let mut j = group - 1;
        while j < n {
            let ju = j as usize;
            ytemp.data[ju] = y.data[ju];
            let mut inc = SUNMAX(
                srur * SUNRabs(y.data[ju]),
                minInc / ark_mem.ewt.data[ju],
            );

            /* Adjust sign(inc) as before. */
            if let Some(constraints) = &ark_mem.constraints {
                let conj = constraints.data[ju];
                if SUNRabs(conj) == ONE {
                    if (ytemp.data[ju] + inc) * conj < ZERO {
                        inc = -inc;
                    }
                } else if SUNRabs(conj) == TWO && (ytemp.data[ju] + inc) * conj <= ZERO {
                    inc = -inc;
                }
            }

            let inc_inv = ONE / inc;
            let i1 = 0.max(j - mupper);
            let i2 = (j + mlower).min(n - 1);
            for i in i1..=i2 {
                let val = inc_inv * (ftemp.data[i as usize] - fy.data[i as usize]);
                jac.set(i, j, val);
            }
            j += width;
        }
    }

    retval
}

/*---------------------------------------------------------------
  arkLsDQJtimes:

  This routine generates a difference quotient approximation to
  the Jacobian-vector product fi_y(t,y) * v. The approximation is
  Jv = [fi(y + v*sig) - fi(y)]/sig, where sig = 1 / ||v||_WRMS,
  i.e. the WRMS norm of v*sig is 1.  (The jtimes == None +
  jtimesDQ dispatch target; ark_mem is handed in by the caller.)
  ---------------------------------------------------------------*/
fn arkLsDQJtimes(
    ark_mem: &mut ARKodeMem,
    ewt: &NVector,
    nfeDQ: &mut i64,
    Jt_f: ARKRhsFn,
    v: &NVector,
    Jv: &mut NVector,
    t: f64,
    y: &NVector,
    fy: &NVector,
    work: &mut NVector,
) -> i32 {
    /* Initialize perturbation to 1/||v|| (the iterative solve detaches
       ark_mem.ewt for the scaling vectors, so it arrives as an
       argument — donor cvode_ls pattern) */
    let mut sig = ONE / N_VWrmsNorm(v, ewt);
    let mut retval = 0;

    for _iter in 0..MAX_DQITERS {
        /* Set work = y + sig*v */
        N_VLinearSum(sig, v, ONE, y, work);

        /* Set Jv = f(tn, y+sig*v), after calling pre-RHS function (if supplied) */
        if let Some(pre_rhs) = ark_mem.PreRhsFn {
            retval = pre_rhs(t, work, &mut ark_mem.user_data);
            if retval != 0 {
                return ARK_PRERHSFN_FAIL;
            }
        }
        retval = Jt_f(t, work, Jv, &mut ark_mem.user_data);
        *nfeDQ += 1;
        if retval == 0 {
            break;
        }
        if retval < 0 {
            return -1;
        }

        /* If fi failed recoverably, shrink sig and retry */
        sig *= PT25;
    }

    /* If retval still isn't 0, return with a recoverable failure */
    if retval > 0 {
        return 1;
    }

    /* Replace Jv by (Jv - fy)/sig */
    let siginv = ONE / sig;
    Jv.linear_sum_with(siginv, -siginv, fy);

    0
}

/*-----------------------------------------------------------------
  arkLsLinSys

  Setup the linear system A = I - gamma J or A = M - gamma J
  (the linsys == None + !user_linsys dispatch target).
  -----------------------------------------------------------------*/
fn arkLsLinSys(
    ark_mem: &mut ARKodeMem,
    arkls_mem: &mut ARKLsMem,
    t: f64,
    y: &mut NVector,
    fy: &NVector,
    M: Option<&SUNMatrix>,
    jok: bool,
    jcur: &mut bool,
    gamma: f64,
    vtemp1: &mut NVector,
    vtemp2: &mut NVector,
    _vtemp3: &mut NVector,
) -> i32 {
    let mut retval;

    /* Check if Jacobian needs to be updated */
    if jok {
        /* Use saved copy of J */
        *jcur = SUNFALSE;

        /* Overwrite linear system matrix with saved J */
        {
            let ARKLsMem { A, savedJ, .. } = arkls_mem;
            retval = SUNMatCopy(savedJ.as_ref().unwrap(), A.as_mut().unwrap());
        }
        if retval != 0 {
            arkProcessError(
                Some(ark_mem),
                ARKLS_SUNMAT_FAIL,
                line!(),
                "arkLsLinSys",
                file!(),
                MSG_LS_SUNMAT_FAILED,
            );
            arkls_mem.last_flag = ARKLS_SUNMAT_FAIL;
            return arkls_mem.last_flag;
        }
    } else {
        /* Call jac() routine to update J */
        *jcur = SUNTRUE;

        /* Clear the linear system matrix if necessary (direct linear solvers) */
        if !arkls_mem.iterative {
            retval = SUNMatZero(arkls_mem.A.as_mut().unwrap());
            if retval != 0 {
                arkProcessError(
                    Some(ark_mem),
                    ARKLS_SUNMAT_FAIL,
                    line!(),
                    "arkLsLinSys",
                    file!(),
                    MSG_LS_SUNMAT_FAILED,
                );
                arkls_mem.last_flag = ARKLS_SUNMAT_FAIL;
                return arkls_mem.last_flag;
            }
        }

        /* Compute new Jacobian matrix */
        retval = if arkls_mem.jacDQ {
            /* internal difference-quotient Jacobian */
            let fi = match ark_mem.step_getimplicitrhs {
                Some(get) => match get(ark_mem) {
                    Some(f) => f,
                    None => {
                        arkProcessError(
                            Some(ark_mem),
                            ARKLS_ILL_INPUT,
                            line!(),
                            "arkLsDQJac",
                            file!(),
                            "Time step module is missing implicit RHS fcn",
                        );
                        return ARKLS_ILL_INPUT;
                    }
                },
                None => return ARKLS_ILL_INPUT,
            };
            let mut a = arkls_mem.A.take().unwrap();
            let ret = arkLsDQJac(ark_mem, arkls_mem, fi, t, y, fy, &mut a, vtemp1, vtemp2);
            arkls_mem.A = Some(a);
            ret
        } else {
            let jac_fn = arkls_mem.jac.unwrap();
            jac_fn(
                t,
                y,
                fy,
                arkls_mem.A.as_mut().unwrap(),
                &mut ark_mem.user_data,
                vtemp1,
                vtemp2,
                _vtemp3,
            )
        };
        if retval < 0 {
            arkProcessError(
                Some(ark_mem),
                ARKLS_JACFUNC_UNRECVR,
                line!(),
                "arkLsLinSys",
                file!(),
                MSG_LS_JACFUNC_FAILED,
            );
            arkls_mem.last_flag = ARKLS_JACFUNC_UNRECVR;
            return -1;
        }
        if retval > 0 {
            arkls_mem.last_flag = ARKLS_JACFUNC_RECVR;
            return 1;
        }

        /* Update saved copy of the Jacobian matrix */
        {
            let ARKLsMem { A, savedJ, .. } = arkls_mem;
            retval = SUNMatCopy(A.as_ref().unwrap(), savedJ.as_mut().unwrap());
        }
        if retval != 0 {
            arkProcessError(
                Some(ark_mem),
                ARKLS_SUNMAT_FAIL,
                line!(),
                "arkLsLinSys",
                file!(),
                MSG_LS_SUNMAT_FAILED,
            );
            arkls_mem.last_flag = ARKLS_SUNMAT_FAIL;
            return arkls_mem.last_flag;
        }
    }

    /* Perform linear combination A = I - gamma*J or A = M - gamma*J */
    retval = match M {
        None => SUNMatScaleAddI(-gamma, arkls_mem.A.as_mut().unwrap()),
        Some(m) => SUNMatScaleAdd(-gamma, arkls_mem.A.as_mut().unwrap(), m),
    };

    /* Check matrix operation return value */
    if retval != 0 {
        arkProcessError(
            Some(ark_mem),
            ARKLS_SUNMAT_FAIL,
            line!(),
            "arkLsLinSys",
            file!(),
            MSG_LS_SUNMAT_FAILED,
        );
        arkls_mem.last_flag = ARKLS_SUNMAT_FAIL;
        return arkls_mem.last_flag;
    }

    ARKLS_SUCCESS
}

/*---------------------------------------------------------------
  arkLsInitialize performs remaining initializations specific
  to the linear solver interface (and solver itself)
  ---------------------------------------------------------------*/
pub fn arkLsInitialize(ark_mem: &mut ARKodeMem) -> i32 {
    let mut arkls_mem = match arkLs_AccessLMem(ark_mem, "arkLsInitialize") {
        Ok(l) => l,
        Err(e) => return e,
    };
    let ret = arkLsInitialize_inner(ark_mem, &mut arkls_mem);
    ark_mem.lmem = Some(arkls_mem);
    ret
}

fn arkLsInitialize_inner(ark_mem: &mut ARKodeMem, arkls_mem: &mut ARKLsMem) -> i32 {
    /* access ARKLsMassMem (if applicable) */
    let arkls_massmem = if ark_mem.step_getmassmem.is_some() {
        ark_mem.mass_mem.take()
    } else {
        None
    };

    /* Test for valid combinations of matrix & Jacobian routines: */
    if arkls_mem.A.is_some() {
        /* Matrix-based case */

        if !arkls_mem.user_linsys {
            /* Internal linear system function, reset pointers (just in case) */
            arkls_mem.linsys = None; /* internal arkLsLinSys */

            /* Check if an internal or user-supplied Jacobian function is used */
            if arkls_mem.jacDQ {
                /* Internal difference quotient Jacobian. Check that A is dense
                   or band, otherwise return an error */
                let ok = matches!(
                    arkls_mem.A,
                    Some(SUNMatrix::Dense(_)) | Some(SUNMatrix::Band(_))
                );
                if ok {
                    arkls_mem.jac = None; /* internal arkLsDQJac */
                } else {
                    arkProcessError(
                        Some(ark_mem),
                        ARKLS_ILL_INPUT,
                        line!(),
                        "arkLsInitialize",
                        file!(),
                        "No Jacobian constructor available for SUNMatrix type",
                    );
                    arkls_mem.last_flag = ARKLS_ILL_INPUT;
                    return ARKLS_ILL_INPUT;
                }
            }

            /* Allocate internally saved Jacobian if not already done */
            if arkls_mem.savedJ.is_none() {
                if let Some(a) = &arkls_mem.A {
                    arkls_mem.savedJ = Some(SUNMatClone(a));
                }
            }
        } /* end matrix-based case */
    } else {
        /* Matrix-free case: ensure 'jac' and 'linsys' function pointers are NULL */
        arkls_mem.jacDQ = SUNFALSE;
        arkls_mem.jac = None;

        arkls_mem.user_linsys = SUNFALSE;
        arkls_mem.linsys = None;
    }

    /* Test for valid combination of system matrix and mass matrix
       (if applicable) */
    if let Some(mm) = &arkls_massmem {
        /* A and M must both be NULL or non-NULL */
        if arkls_mem.A.is_some() != mm.M.is_some() {
            ark_mem.mass_mem = arkls_massmem;
            arkProcessError(
                Some(ark_mem),
                ARKLS_ILL_INPUT,
                line!(),
                "arkLsInitialize",
                file!(),
                "Cannot combine NULL and non-NULL System and mass matrices",
            );
            arkls_mem.last_flag = ARKLS_ILL_INPUT;
            return ARKLS_ILL_INPUT;
        }

        /* If A is non-NULL, A and M must have matching types */
        if let (Some(a), Some(m)) = (&arkls_mem.A, &mm.M) {
            if SUNMatGetID(a) != SUNMatGetID(m) {
                ark_mem.mass_mem = arkls_massmem;
                arkProcessError(
                    Some(ark_mem),
                    ARKLS_ILL_INPUT,
                    line!(),
                    "arkLsInitialize",
                    file!(),
                    "System and mass matrices have incompatible types",
                );
                arkls_mem.last_flag = ARKLS_ILL_INPUT;
                return ARKLS_ILL_INPUT;
            }
        }

        /* If either system or mass matrix solver is matrix-embedded,
           then both must be */
        if (arkls_mem.LS.ls_type() == SUNLINEARSOLVER_MATRIX_EMBEDDED)
            != (mm.LS.ls_type() == SUNLINEARSOLVER_MATRIX_EMBEDDED)
        {
            ark_mem.mass_mem = arkls_massmem;
            arkProcessError(
                Some(ark_mem),
                ARKLS_ILL_INPUT,
                line!(),
                "arkLsInitialize",
                file!(),
                "mismatched matrix-embedded LS types (system and mass must match)",
            );
            arkls_mem.last_flag = ARKLS_ILL_INPUT;
            return ARKLS_ILL_INPUT;
        }
    }
    ark_mem.mass_mem = arkls_massmem;

    /* reset counters */
    arkLsInitializeCounters(arkls_mem);

    /* Set Jacobian-vector product related fields, based on jtimesDQ */
    if arkls_mem.jtimesDQ {
        arkls_mem.jtsetup = None;
        arkls_mem.jtimes = None; /* internal arkLsDQJtimes */
    }

    /* If A is NULL and psetup is not present, then arkLsSetup does
       not need to be called, so set the lsetup function to NULL (if possible) */
    if arkls_mem.A.is_none() && arkls_mem.pset.is_none() {
        if let Some(disable) = ark_mem.step_disablelsetup {
            disable(ark_mem);
        }
    }

    /* When using a matrix-embedded linear solver, disable lsetup call and
       solution scaling */
    if arkls_mem.LS.ls_type() == SUNLINEARSOLVER_MATRIX_EMBEDDED {
        if let Some(disable) = ark_mem.step_disablelsetup {
            disable(ark_mem);
        }
        arkls_mem.scalesol = SUNFALSE;
    }

    /* Call LS initialize routine, and return result */
    arkls_mem.last_flag = arkls_mem.LS.initialize();
    arkls_mem.last_flag
}

/*---------------------------------------------------------------
  arkLsSetup conditionally calls the LS 'setup' routine.

  When using a SUNMatrix object, this determines whether
  to update a Jacobian matrix (or use a stored version), based
  on heuristics regarding previous convergence issues, the number
  of time steps since it was last updated, etc.; it then creates
  the system matrix from this, the 'gamma' factor and the
  mass/identity matrix, A = M-gamma*J.

  This routine then calls the LS 'setup' routine with A.
  ---------------------------------------------------------------*/
pub fn arkLsSetup(
    ark_mem: &mut ARKodeMem,
    convfail: i32,
    tpred: f64,
    ypred: &mut NVector,
    fpred: &NVector,
    jcurPtr: &mut bool,
    vtemp1: &mut NVector,
    vtemp2: &mut NVector,
    vtemp3: &mut NVector,
) -> i32 {
    let mut arkls_mem = match arkLs_AccessLMem(ark_mem, "arkLsSetup") {
        Ok(l) => l,
        Err(e) => return e,
    };
    let ret = arkLsSetup_inner(
        ark_mem,
        &mut arkls_mem,
        convfail,
        tpred,
        ypred,
        fpred,
        jcurPtr,
        vtemp1,
        vtemp2,
        vtemp3,
    );
    ark_mem.lmem = Some(arkls_mem);
    ret
}

#[allow(clippy::too_many_arguments)]
fn arkLsSetup_inner(
    ark_mem: &mut ARKodeMem,
    arkls_mem: &mut ARKLsMem,
    convfail: i32,
    tpred: f64,
    ypred: &mut NVector,
    fpred: &NVector,
    jcurPtr: &mut bool,
    vtemp1: &mut NVector,
    vtemp2: &mut NVector,
    vtemp3: &mut NVector,
) -> i32 {
    /* Immediately return when using matrix-embedded linear solver */
    if arkls_mem.LS.ls_type() == SUNLINEARSOLVER_MATRIX_EMBEDDED {
        arkls_mem.last_flag = ARKLS_SUCCESS;
        return arkls_mem.last_flag;
    }

    /* Set ARKLs time to current time (ycur/fcur are the ypred/fpred
       arguments here) */
    arkls_mem.tcur = tpred;

    /* get gamma values from time step module */
    let getgammas = ark_mem.step_getgammas.unwrap();
    let mut gamma = ZERO;
    let mut gamrat = ZERO;
    let mut jcur = SUNFALSE;
    let mut dgamma_fail = SUNFALSE;
    arkls_mem.last_flag = getgammas(ark_mem, &mut gamma, &mut gamrat, &mut jcur, &mut dgamma_fail);
    if arkls_mem.last_flag != 0 {
        arkProcessError(
            Some(ark_mem),
            arkls_mem.last_flag,
            line!(),
            "arkLsSetup",
            file!(),
            "An error occurred in ark_step_getgammas",
        );
        return arkls_mem.last_flag;
    }

    /* Use initsetup, gamma/gammap, and convfail to set J/P eval. flag jok;
       Note: the "ARK_FAIL_BAD_J" test is asking whether the nonlinear
       solver converged due to a bad system Jacobian AND our gamma was
       fine, indicating that the J and/or P were invalid */
    arkls_mem.jbad = ark_mem.initsetup
        || (ark_mem.nst >= arkls_mem.nstlj + arkls_mem.msbj)
        || (convfail == ARK_FAIL_BAD_J && !dgamma_fail)
        || (convfail == ARK_FAIL_OTHER);

    /* Check for mass matrix module and setup mass matrix */
    let have_mass = ark_mem.step_getmassmem.is_some() && ark_mem.mass_mem.is_some();
    if have_mass {
        /* Setup mass matrix linear solver (including recomputation of
           mass matrix) */
        arkls_mem.last_flag = arkLsMassSetup(ark_mem, tpred, vtemp1, vtemp2, vtemp3);
        if arkls_mem.last_flag != 0 {
            arkProcessError(
                Some(ark_mem),
                ARKLS_SUNMAT_FAIL,
                line!(),
                "arkLsSetup",
                file!(),
                "Error setting up mass-matrix linear solver",
            );
            return arkls_mem.last_flag;
        }
    }

    /* Setup the linear system if necessary */
    if arkls_mem.A.is_some() {
        /* Set shortcut to the mass matrix (NULL if matrix-free); take it
           out of ark_mem for the duration of the linsys call */
        let mass_taken = if have_mass { ark_mem.mass_mem.take() } else { None };
        let m_ref: Option<&SUNMatrix> = mass_taken.as_ref().and_then(|mm| mm.M.as_ref());

        /* Update J if appropriate and evaluate A = I-gamma*J or A = M-gamma*J */
        let jok = !arkls_mem.jbad;
        let retval = if arkls_mem.user_linsys {
            let linsys = arkls_mem.linsys.unwrap();
            linsys(
                tpred,
                ypred,
                fpred,
                arkls_mem.A.as_mut().unwrap(),
                m_ref,
                jok,
                jcurPtr,
                gamma,
                &mut ark_mem.user_data,
                vtemp1,
                vtemp2,
                vtemp3,
            )
        } else {
            arkLsLinSys(
                ark_mem, arkls_mem, tpred, ypred, fpred, m_ref, jok, jcurPtr, gamma, vtemp1,
                vtemp2, vtemp3,
            )
        };
        if mass_taken.is_some() {
            ark_mem.mass_mem = mass_taken;
        }

        /* Update J eval count and step when J was last updated */
        if *jcurPtr {
            arkls_mem.nje += 1;
            arkls_mem.nstlj = ark_mem.nst;
            arkls_mem.tnlj = tpred;
        }

        /* Check linsys() return value and return if necessary */
        if retval != ARKLS_SUCCESS {
            if arkls_mem.user_linsys {
                if retval < 0 {
                    arkProcessError(
                        Some(ark_mem),
                        ARKLS_JACFUNC_UNRECVR,
                        line!(),
                        "arkLsSetup",
                        file!(),
                        MSG_LS_JACFUNC_FAILED,
                    );
                    arkls_mem.last_flag = ARKLS_JACFUNC_UNRECVR;
                    return -1;
                } else {
                    arkls_mem.last_flag = ARKLS_JACFUNC_RECVR;
                    return 1;
                }
            } else {
                return retval;
            }
        }
    } else {
        /* Matrix-free case, set jcur to jbad */
        *jcurPtr = arkls_mem.jbad;
    }

    /* Call LS setup routine -- for direct solvers this factors A; for
       iterative solvers the generic C SUNLinSolSetup calls arkLsPSetup,
       which passes the heuristic suggestions above to the user code */
    arkls_mem.last_flag = match &arkls_mem.LS {
        LinearSolver::Dense(_) | LinearSolver::Band(_) => {
            let ARKLsMem { LS, A, .. } = arkls_mem;
            LS.setup(A.as_mut())
        }
        _ => {
            /* iterative solver: preconditioner setup (arkLsPSetup); it is
               only invoked when a user psetup routine exists.  C's pset
               writes through the stepper's &jcur from step_getgammas; the
               Rust write-back goes through the step_setjcur op. */
            if let Some(pset) = arkls_mem.pset {
                let mut jcur_ls = jcur;
                let retval = pset(
                    arkls_mem.tcur,
                    ypred,
                    fpred,
                    !arkls_mem.jbad,
                    &mut jcur_ls,
                    gamma,
                    &mut ark_mem.user_data,
                );
                if let Some(setjcur) = ark_mem.step_setjcur {
                    setjcur(ark_mem, jcur_ls);
                }
                *jcurPtr = jcur_ls;
                retval
            } else {
                SUN_SUCCESS
            }
        }
    };

    /* If the SUNMatrix was NULL, update heuristics flags */
    if arkls_mem.A.is_none() {
        /* If user set jcur to SUNTRUE, increment npe and save nst value */
        if *jcurPtr {
            arkls_mem.npe += 1;
            arkls_mem.nstlj = ark_mem.nst;
            arkls_mem.tnlj = tpred;
        }

        /* Update jcurPtr flag if we suggested an update */
        if arkls_mem.jbad {
            *jcurPtr = SUNTRUE;
        }
    }

    arkls_mem.last_flag
}

/*---------------------------------------------------------------
  arkLsSolve: interfaces between ARKODE and the generic
  SUNLinearSolver object LS, by setting the appropriate tolerance
  and scaling vectors, calling the solver, and accumulating
  statistics from the solve for use/reporting by ARKODE.

  When using a non-NULL SUNMatrix, this will additionally scale
  the solution appropriately when gamrat != 1.
  ---------------------------------------------------------------*/
pub fn arkLsSolve(
    ark_mem: &mut ARKodeMem,
    b: &mut NVector,
    tnow: f64,
    ynow: &NVector,
    fnow: &NVector,
    eRNrm: f64,
    mnewt: i32,
) -> i32 {
    let mut arkls_mem = match arkLs_AccessLMem(ark_mem, "arkLsSolve") {
        Ok(l) => l,
        Err(e) => return e,
    };
    let ret = arkLsSolve_inner(ark_mem, &mut arkls_mem, b, tnow, ynow, fnow, eRNrm, mnewt);
    ark_mem.lmem = Some(arkls_mem);
    ret
}

#[allow(clippy::too_many_arguments)]
fn arkLsSolve_inner(
    ark_mem: &mut ARKodeMem,
    arkls_mem: &mut ARKLsMem,
    b: &mut NVector,
    tnow: f64,
    ynow: &NVector,
    fnow: &NVector,
    eRNrm: f64,
    mnewt: i32,
) -> i32 {
    /* Set scalar tcur for use by the Atimes and Psolve interface
       routines (ycur/fcur are captured by the closures below) */
    arkls_mem.tcur = tnow;

    /* If the linear solver is iterative:
       test norm(b), if small, return x = 0 or x = b;
       set linear solver tolerance (in left/right scaled 2-norm) */
    let delta;
    if arkls_mem.iterative {
        let deltar = arkls_mem.eplifac * eRNrm;
        let bnorm = N_VWrmsNorm(b, ark_rwt(ark_mem));

        if bnorm <= deltar {
            if mnewt > 0 {
                N_VConst(ZERO, b);
            }
            arkls_mem.last_flag = ARKLS_SUCCESS;
            return arkls_mem.last_flag;
        }
        /* Adjust tolerance for 2-norm */
        delta = deltar * arkls_mem.nrmfac;
    } else {
        delta = ZERO;
    }

    /* (Scaling vectors: the Rust iterative solvers all accept the
       s1/s2 scaling vectors at solve time — the C
       "solver does not support scaling vectors" rwt_mean fallback
       branch is unreachable.) */

    /* Set initial guess x = 0 to LS, and zero initial guess flag */
    N_VConst(ZERO, &mut arkls_mem.x);
    arkls_mem.LS.set_zero_guess(SUNTRUE);

    /* If a user-provided jtsetup routine is supplied, call that here */
    if let Some(jtsetup) = arkls_mem.jtsetup {
        arkls_mem.last_flag = jtsetup(tnow, ynow, fnow, &mut ark_mem.user_data);
        arkls_mem.njtsetup += 1;
        if arkls_mem.last_flag != 0 {
            arkProcessError(
                Some(ark_mem),
                arkls_mem.last_flag,
                line!(),
                "arkLsSolve",
                file!(),
                MSG_LS_JTSETUP_FAILED,
            );
            return arkls_mem.last_flag;
        }
    }

    /* Call solver, and copy x to b */
    let retval = if arkls_mem.LS.ls_type() == SUNLINEARSOLVER_DIRECT {
        let ARKLsMem { LS, A, x, .. } = arkls_mem;
        match (LS, A.as_mut()) {
            (LinearSolver::Dense(dls), Some(SUNMatrix::Dense(am))) => dls.solve(am, x, b),
            (LinearSolver::Band(bls), Some(SUNMatrix::Band(am))) => bls.solve(am, x, b),
            _ => SUN_ERR_ARG_INCOMPATIBLE,
        }
    } else {
        arkLsSolveIterative(ark_mem, arkls_mem, b, delta, tnow, ynow, fnow)
    };
    {
        let ARKLsMem { x, .. } = arkls_mem;
        N_VScale(ONE, x, b);
    }

    /* If using a direct or matrix-iterative solver, scale the correction to
       account for change in gamma (this is only beneficial if M==I) */
    if arkls_mem.scalesol {
        let getgammas = ark_mem.step_getgammas.unwrap();
        let mut gamma = ZERO;
        let mut gamrat = ZERO;
        let mut jcur = SUNFALSE;
        let mut dgamma_fail = SUNFALSE;
        arkls_mem.last_flag =
            getgammas(ark_mem, &mut gamma, &mut gamrat, &mut jcur, &mut dgamma_fail);
        if arkls_mem.last_flag != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                arkls_mem.last_flag,
                line!(),
                "arkLsSolve",
                file!(),
                "An error occurred in ark_step_getgammas",
            );
            return arkls_mem.last_flag;
        }
        if gamrat != ONE {
            b.scale_inplace(TWO / (ONE + gamrat));
        }
    }

    /* Retrieve statistics from iterative linear solvers */
    let mut nli_inc = 0;
    if arkls_mem.iterative {
        nli_inc = arkls_mem.LS.num_iters();
    }

    /* Increment counters nli and ncfl */
    arkls_mem.nli += nli_inc as i64;
    if retval != SUN_SUCCESS {
        arkls_mem.ncfl += 1;
    }

    /* Interpret solver return value  */
    arkls_mem.last_flag = retval;

    match retval {
        SUN_SUCCESS => 0,
        SUNLS_RES_REDUCED => {
            /* allow reduction but not solution on first nonlinear iteration,
               otherwise return with a recoverable failure */
            if mnewt == 0 {
                0
            } else {
                1
            }
        }
        SUNLS_CONV_FAIL | SUNLS_ATIMES_FAIL_REC | SUNLS_PSOLVE_FAIL_REC
        | SUNLS_PACKAGE_FAIL_REC | SUNLS_QRFACT_FAIL | SUNLS_LUFACT_FAIL => 1,
        SUN_ERR_ARG_CORRUPT | SUN_ERR_ARG_INCOMPATIBLE | SUN_ERR_MEM_FAIL | SUNLS_GS_FAIL
        | SUNLS_QRSOL_FAIL => -1,
        SUN_ERR_EXT_FAIL => {
            arkProcessError(
                Some(ark_mem),
                SUN_ERR_EXT_FAIL,
                line!(),
                "arkLsSolve",
                file!(),
                "Failure in SUNLinSol external package",
            );
            -1
        }
        SUNLS_ATIMES_FAIL_UNREC => {
            arkProcessError(
                Some(ark_mem),
                SUNLS_ATIMES_FAIL_UNREC,
                line!(),
                "arkLsSolve",
                file!(),
                MSG_LS_JTIMES_FAILED,
            );
            -1
        }
        SUNLS_PSOLVE_FAIL_UNREC => {
            arkProcessError(
                Some(ark_mem),
                SUNLS_PSOLVE_FAIL_UNREC,
                line!(),
                "arkLsSolve",
                file!(),
                MSG_LS_PSOLVE_FAILED,
            );
            -1
        }
        _ => 0,
    }
}

/* Iterative Krylov solve: builds the ATimes (arkLsATimes) and PSolve
   (arkLsPSolve) callbacks over the integrator memory.  The two
   closures share ark_mem through a RefCell — the solvers never call
   them re-entrantly.  The error weight vectors are detached from
   ARKodeMem for the duration of the solve so they can be passed as
   the (read-only) scaling vectors s1 = rwt, s2 = ewt. */
#[allow(clippy::too_many_arguments)]
fn arkLsSolveIterative(
    ark_mem: &mut ARKodeMem,
    arkls_mem: &mut ARKLsMem,
    b: &NVector,
    delta: f64,
    tnow: f64,
    ynow: &NVector,
    fnow: &NVector,
) -> i32 {
    let ARKLsMem {
        LS,
        x,
        ytemp,
        njtimes,
        nfeDQ,
        nps,
        jtimes,
        jtimesDQ,
        Jt_f,
        psolve,
        ..
    } = arkls_mem;
    let jtimes = *jtimes;
    let jtimes_dq = *jtimesDQ;
    let jt_f = *Jt_f;
    let psolve_fn = *psolve;

    let ewt = std::mem::take(&mut ark_mem.ewt);
    /* s1 = rwt: alias of ewt when rwt_is_ewt (Addendum C.1) */
    let rwt_detached = if ark_mem.rwt_is_ewt {
        None
    } else {
        Some(std::mem::take(&mut ark_mem.rwt))
    };
    let arm = RefCell::new(&mut *ark_mem);
    let ewt_ref = &ewt;
    let rwt_ref: &NVector = match &rwt_detached {
        Some(r) => r,
        None => ewt_ref,
    };

    let mut atimes = |v: &NVector, z: &mut NVector| -> i32 {
        let mut guard = arm.borrow_mut();
        let ar: &mut ARKodeMem = &mut guard;

        /* get gamma values from time step module */
        let getgammas = ar.step_getgammas.unwrap();
        let mut gamma = ZERO;
        let mut gamrat = ZERO;
        let mut jcur = SUNFALSE;
        let mut dgamma_fail = SUNFALSE;
        let gret = getgammas(ar, &mut gamma, &mut gamrat, &mut jcur, &mut dgamma_fail);
        if gret != ARK_SUCCESS {
            arkProcessError(
                Some(ar),
                gret,
                line!(),
                "arkLsATimes",
                file!(),
                "An error occurred in ark_step_getgammas",
            );
            return gret;
        }

        /* call Jacobian-times-vector product routine
           (either user-supplied or internal DQ) */
        let jret = if jtimes_dq {
            arkLsDQJtimes(ar, ewt_ref, nfeDQ, jt_f.unwrap(), v, z, tnow, ynow, fnow, ytemp)
        } else {
            let jt = jtimes.unwrap();
            jt(v, z, tnow, ynow, fnow, &mut ar.user_data, ytemp)
        };
        *njtimes += 1;
        if jret != 0 {
            return jret;
        }

        /* Compute mass matrix vector product and add to result */
        if ar.step_getmassmem.is_some() && ar.mass_mem.is_some() {
            let mut mm = ar.mass_mem.take().unwrap();
            let mret = arkLsMTimes_inner(ar, &mut mm, v, ytemp);
            ar.mass_mem = Some(mm);
            if mret != 0 {
                return mret;
            }
            /* z = ytemp - gamma*z */
            z.linear_sum_with(-gamma, ONE, ytemp);
        } else {
            /* z = v - gamma*z */
            z.linear_sum_with(-gamma, ONE, v);
        }
        0
    };

    let mut psolve_cb = |r: &NVector, z: &mut NVector, tol: f64, lr: i32| -> i32 {
        let mut guard = arm.borrow_mut();
        let ar: &mut ARKodeMem = &mut guard;

        /* get gamma values from time step module */
        let getgammas = ar.step_getgammas.unwrap();
        let mut gamma = ZERO;
        let mut gamrat = ZERO;
        let mut jcur = SUNFALSE;
        let mut dgamma_fail = SUNFALSE;
        let gret = getgammas(ar, &mut gamma, &mut gamrat, &mut jcur, &mut dgamma_fail);
        if gret != ARK_SUCCESS {
            arkProcessError(
                Some(ar),
                gret,
                line!(),
                "arkLsPSolve",
                file!(),
                "An error occurred in ark_step_getgammas",
            );
            return gret;
        }

        /* call the user-supplied psolve routine, and accumulate count */
        let ret = if let Some(ps) = psolve_fn {
            ps(tnow, ynow, fnow, r, z, gamma, tol, lr, &mut ar.user_data)
        } else {
            0
        };
        *nps += 1;
        ret
    };

    let retval = if psolve_fn.is_some() {
        LS.solve(
            None,
            x,
            b,
            delta,
            &mut atimes,
            Some(&mut psolve_cb),
            Some(rwt_ref),
            Some(ewt_ref),
        )
    } else {
        LS.solve(None, x, b, delta, &mut atimes, None, Some(rwt_ref), Some(ewt_ref))
    };

    /* end the closures' shared borrow of the RefCell before unwrapping */
    let _ = (atimes, psolve_cb);
    let ark_mem_back = arm.into_inner();
    ark_mem_back.ewt = ewt;
    if let Some(r) = rwt_detached {
        ark_mem_back.rwt = r;
    }

    retval
}

/*---------------------------------------------------------------
  arkLsFree frees memory associates with the ARKLs system
  solver interface.
  ---------------------------------------------------------------*/
pub fn arkLsFree(ark_mem: &mut ARKodeMem) -> i32 {
    /* Return immediately if ARKodeMem, ARKLsMem are NULL */
    let taken = match ark_mem.step_getlinmem {
        Some(get) => get(ark_mem),
        None => None,
    };
    let mut arkls_mem = match taken {
        Some(l) => l,
        None => return ARKLS_SUCCESS,
    };

    /* Free N_Vector memory (with lrw/liw accounting) */
    arkFreeVec(ark_mem, &mut arkls_mem.ytemp);
    arkFreeVec(ark_mem, &mut arkls_mem.x);

    /* savedJ / A / LS memory is dropped with the box (C frees savedJ,
       nullifies the borrowed pointers, and calls any pfree here) */
    ARKLS_SUCCESS
}

/*---------------------------------------------------------------
  arkLsInitializeCounters resets all counters from an ARKLsMem
  structure.
  ---------------------------------------------------------------*/
pub fn arkLsInitializeCounters(arkls_mem: &mut ARKLsMem) -> i32 {
    arkls_mem.nje = 0;
    arkls_mem.nfeDQ = 0;
    arkls_mem.nstlj = 0;
    arkls_mem.npe = 0;
    arkls_mem.nli = 0;
    arkls_mem.nps = 0;
    arkls_mem.ncfl = 0;
    arkls_mem.njtsetup = 0;
    arkls_mem.njtimes = 0;
    0
}

/*===============================================================
  Mass-matrix linear solver interface (the ARKLsMassMem half)
  ===============================================================*/

/* C: arkLs_AccessMassMem / arkLs_AccessARKODEMassMem.  Takes the
   ARKLsMassMem box out through the step_getmassmem op; callers put
   it back by writing ark_mem.mass_mem. */
pub(crate) fn arkLs_AccessMassMem(
    ark_mem: &mut ARKodeMem,
    fname: &str,
) -> Result<Box<crate::arkode_ls_impl::ARKLsMassMem>, i32> {
    let taken = match ark_mem.step_getmassmem {
        Some(get) => get(ark_mem),
        None => None,
    };
    match taken {
        Some(m) => Ok(m),
        None => {
            arkProcessError(
                Some(ark_mem),
                ARKLS_MASSMEM_NULL,
                line!(),
                fname,
                file!(),
                "Mass matrix solver memory is NULL.",
            );
            Err(ARKLS_MASSMEM_NULL)
        }
    }
}

/*---------------------------------------------------------------
  ARKodeSetMassLinearSolver specifies the iterative mass-matrix
  linear solver and user-supplied routine to perform the
  mass-matrix-vector product.
  ---------------------------------------------------------------*/
pub fn ARKodeSetMassLinearSolver(
    ark_mem: &mut ARKodeMem,
    LS: LinearSolver,
    M: Option<SUNMatrix>,
    time_dep: bool,
) -> i32 {
    use crate::arkode_ls_impl::ARKLsMassMem;

    /* Guard against use for time steppers that do not support mass matrices */
    if !ark_mem.step_supports_massmatrix {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!(),
            "ARKodeSetMassLinearSolver",
            file!(),
            "time-stepping module does not support non-identity mass matrices",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* Retrieve the LS type; set flags based on LS type */
    let LSType = LS.ls_type();
    let iterative = LSType != SUNLINEARSOLVER_DIRECT;
    let matrixbased =
        LSType != SUNLINEARSOLVER_ITERATIVE && LSType != SUNLINEARSOLVER_MATRIX_EMBEDDED;

    /* Ensure that M is NULL when LS is matrix-embedded */
    if LSType == SUNLINEARSOLVER_MATRIX_EMBEDDED && M.is_some() {
        arkProcessError(
            Some(ark_mem),
            ARKLS_ILL_INPUT,
            line!(),
            "ARKodeSetMassLinearSolver",
            file!(),
            "Incompatible inputs: matrix-embedded LS requires NULL matrix",
        );
        return ARKLS_ILL_INPUT;
    }

    /* Check for compatible LS type, matrix and "atimes" support */
    if iterative {
        if matrixbased && M.is_none() {
            arkProcessError(
                Some(ark_mem),
                ARKLS_ILL_INPUT,
                line!(),
                "ARKodeSetMassLinearSolver",
                file!(),
                "Incompatible inputs: matrix-iterative LS requires non-NULL matrix",
            );
            return ARKLS_ILL_INPUT;
        }
    } else if M.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARKLS_ILL_INPUT,
            line!(),
            "ARKodeSetMassLinearSolver",
            file!(),
            "Incompatible inputs: direct LS requires non-NULL matrix",
        );
        return ARKLS_ILL_INPUT;
    }

    /* Test whether time stepper module is supplied, with required routines */
    if ark_mem.step_attachmasssol.is_none() || ark_mem.step_getmassmem.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARKLS_ILL_INPUT,
            line!(),
            "ARKodeSetMassLinearSolver",
            file!(),
            "Missing time step module or associated routines",
        );
        return ARKLS_ILL_INPUT;
    }

    /* When using a non-NULL SUNMatrix object, for direct linear solvers
       create M_lu to store the factorization of M (C aliases M_lu = M
       for iterative solvers; the Rust iterative solve never touches the
       matrix operand) */
    let M_lu = match (&M, iterative) {
        (Some(m), false) => Some(SUNMatClone(m)),
        _ => None,
    };

    /* Allocate memory for ARKLsMassMemRec, set defaults (memset-0 +
       the explicit C initializations) */
    let mut arkls_mem = Box::new(ARKLsMassMem {
        iterative,
        matrixbased,
        mass: None,
        M,
        M_lu,
        eplifac: ARKLS_EPLIN,
        nrmfac: ZERO,
        time_dependent: time_dep,
        msetuptime: ZERO,
        nmsetups: 0,
        nmsolves: 0,
        nmtsetup: 0,
        nmtimes: 0,
        nmvsetup: 0,
        npe: 0,
        nli: 0,
        nps: 0,
        ncfl: 0,
        LS,
        x: NVector::new(0),
        pset: None,
        psolve: None,
        mtsetup: None,
        mtimes: None,
        last_flag: ARKLS_SUCCESS,
    });

    /* Initialize counters */
    arkLsInitializeMassCounters(&mut arkls_mem);

    /* (C attaches NULL ATimes / preconditioner hooks to the LS object
       here; the Rust iterative solvers receive those callbacks at
       solve time.) */

    /* Allocate memory for x (arkAllocVec cannot fail here) */
    let tmpl_len = ark_mem.tempv1.data.len();
    arkAllocVec(ark_mem, tmpl_len, &mut arkls_mem.x);

    /* For iterative LS, compute default norm conversion factor */
    if iterative {
        arkls_mem.nrmfac = SUNRsqrt(N_VGetLength(&arkls_mem.x) as f64);
    }

    /* Attach ARKLs interface to time stepper module */
    let attach = ark_mem.step_attachmasssol.unwrap();
    let retval = attach(
        ark_mem,
        Some(arkLsMassInitialize),
        Some(arkLsMassSetup),
        Some(arkLsMTimes),
        Some(arkLsMassSolve),
        Some(arkLsMassFree),
        time_dep,
        LSType,
        arkls_mem,
    );
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            retval,
            line!(),
            "ARKodeSetMassLinearSolver",
            file!(),
            "Failed to attach to time stepper module",
        );
        return retval;
    }

    ARKLS_SUCCESS
}

/* Shared entry guard for the mass option routines:
   supports_massmatrix check + mass_mem take. */
macro_rules! arkls_mass_option_fn {
    ($ark_mem:ident, $fname:literal) => {{
        if !$ark_mem.step_supports_massmatrix {
            arkProcessError(
                Some($ark_mem),
                ARK_STEPPER_UNSUPPORTED,
                line!(),
                $fname,
                file!(),
                "time-stepping module does not support non-identity mass matrices",
            );
            return ARK_STEPPER_UNSUPPORTED;
        }
        match arkLs_AccessMassMem($ark_mem, $fname) {
            Ok(m) => m,
            Err(e) => return e,
        }
    }};
}

/*---------------------------------------------------------------
  ARKodeSetMassFn specifies the mass matrix function.
  ---------------------------------------------------------------*/
pub fn ARKodeSetMassFn(
    ark_mem: &mut ARKodeMem,
    mass: Option<crate::arkode_ls_impl::ARKLsMassFn>,
) -> i32 {
    let mut arkls_mem = arkls_mass_option_fn!(ark_mem, "ARKodeSetMassFn");

    /* return with failure if mass cannot be used */
    if mass.is_none() {
        ark_mem.mass_mem = Some(arkls_mem);
        arkProcessError(
            Some(ark_mem),
            ARKLS_ILL_INPUT,
            line!(),
            "ARKodeSetMassFn",
            file!(),
            "Mass-matrix routine must be non-NULL",
        );
        return ARKLS_ILL_INPUT;
    }
    if arkls_mem.M.is_none() {
        ark_mem.mass_mem = Some(arkls_mem);
        arkProcessError(
            Some(ark_mem),
            ARKLS_ILL_INPUT,
            line!(),
            "ARKodeSetMassFn",
            file!(),
            "Mass-matrix routine cannot be supplied for NULL SUNMatrix",
        );
        return ARKLS_ILL_INPUT;
    }

    /* set mass matrix routine pointer and return */
    arkls_mem.mass = mass;

    ark_mem.mass_mem = Some(arkls_mem);
    ARKLS_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetMassEpsLin specifies the nonlinear -> linear tolerance
  scale factor for mass matrix linear systems.
  ---------------------------------------------------------------*/
pub fn ARKodeSetMassEpsLin(ark_mem: &mut ARKodeMem, eplifac: f64) -> i32 {
    let mut arkls_mem = arkls_mass_option_fn!(ark_mem, "ARKodeSetMassEpsLin");

    /* store input and return */
    arkls_mem.eplifac = if eplifac <= ZERO { ARKLS_EPLIN } else { eplifac };

    ark_mem.mass_mem = Some(arkls_mem);
    ARKLS_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetMassLSNormFactor sets or computes the factor to use
  when converting from the integrator tolerance (WRMS norm) to the
  linear solver tolerance (L2 norm).
  ---------------------------------------------------------------*/
pub fn ARKodeSetMassLSNormFactor(ark_mem: &mut ARKodeMem, nrmfac: f64) -> i32 {
    let mut arkls_mem = arkls_mass_option_fn!(ark_mem, "ARKodeSetMassLSNormFactor");

    /* store input and return */
    if nrmfac > ZERO {
        /* set user-provided factor */
        arkls_mem.nrmfac = nrmfac;
    } else if nrmfac < ZERO {
        /* compute factor for WRMS norm with dot product */
        N_VConst(ONE, &mut ark_mem.tempv1);
        arkls_mem.nrmfac = SUNRsqrt(N_VDotProd(&ark_mem.tempv1, &ark_mem.tempv1));
    } else {
        /* compute default factor for WRMS norm from vector length */
        arkls_mem.nrmfac = SUNRsqrt(N_VGetLength(&ark_mem.tempv1) as f64);
    }

    ark_mem.mass_mem = Some(arkls_mem);
    ARKLS_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetMassPreconditioner specifies the user-supplied
  preconditioner setup and solve routines.
  ---------------------------------------------------------------*/
pub fn ARKodeSetMassPreconditioner(
    ark_mem: &mut ARKodeMem,
    psetup: Option<crate::arkode_ls_impl::ARKLsMassPrecSetupFn>,
    psolve: Option<crate::arkode_ls_impl::ARKLsMassPrecSolveFn>,
) -> i32 {
    let mut arkls_mem = arkls_mass_option_fn!(ark_mem, "ARKodeSetMassPreconditioner");

    /* issue error if LS object does not allow user-supplied preconditioning */
    if !arkls_mem.iterative {
        ark_mem.mass_mem = Some(arkls_mem);
        arkProcessError(
            Some(ark_mem),
            ARKLS_ILL_INPUT,
            line!(),
            "ARKodeSetMassPreconditioner",
            file!(),
            "SUNLinearSolver object does not support user-supplied preconditioning",
        );
        return ARKLS_ILL_INPUT;
    }

    /* store function pointers for user-supplied routines */
    arkls_mem.pset = psetup;
    arkls_mem.psolve = psolve;

    /* (the Rust iterative solvers receive the psolve closure at solve time) */

    ark_mem.mass_mem = Some(arkls_mem);
    ARKLS_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetMassTimes specifies the user-supplied mass
  matrix-vector product setup and multiply routines.
  ---------------------------------------------------------------*/
pub fn ARKodeSetMassTimes(
    ark_mem: &mut ARKodeMem,
    mtsetup: Option<crate::arkode_ls_impl::ARKLsMassTimesSetupFn>,
    mtimes: Option<crate::arkode_ls_impl::ARKLsMassTimesVecFn>,
) -> i32 {
    let mut arkls_mem = arkls_mass_option_fn!(ark_mem, "ARKodeSetMassTimes");

    /* issue error if mtimes function is unusable */
    if mtimes.is_none() {
        ark_mem.mass_mem = Some(arkls_mem);
        arkProcessError(
            Some(ark_mem),
            ARKLS_ILL_INPUT,
            line!(),
            "ARKodeSetMassTimes",
            file!(),
            "non-NULL mtimes function must be supplied",
        );
        return ARKLS_ILL_INPUT;
    }

    /* issue error if LS object does not allow user-supplied ATimes */
    if !arkls_mem.iterative {
        ark_mem.mass_mem = Some(arkls_mem);
        arkProcessError(
            Some(ark_mem),
            ARKLS_ILL_INPUT,
            line!(),
            "ARKodeSetMassTimes",
            file!(),
            "SUNLinearSolver object does not support user-supplied ATimes routine",
        );
        return ARKLS_ILL_INPUT;
    }

    /* store pointers for user-supplied routines in ARKLs interface
       (C also stores mtimes_data; collapsed onto ark_mem.user_data) */
    arkls_mem.mtsetup = mtsetup;
    arkls_mem.mtimes = mtimes;

    /* (the Rust iterative solvers receive the ATimes closure at solve time) */

    ark_mem.mass_mem = Some(arkls_mem);
    ARKLS_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeGetMassWorkSpace
  ---------------------------------------------------------------*/
pub fn ARKodeGetMassWorkSpace(ark_mem: &mut ARKodeMem, lenrw: &mut i64, leniw: &mut i64) -> i32 {
    /* Return 0 for incompatible steppers */
    if !ark_mem.step_supports_massmatrix {
        *lenrw = 0;
        *leniw = 0;
        return ARK_SUCCESS;
    }

    let arkls_mem = match arkLs_AccessMassMem(ark_mem, "ARKodeGetMassWorkSpace") {
        Ok(m) => m,
        Err(e) => return e,
    };

    /* start with fixed sizes plus vector/matrix pointers */
    *lenrw = 2;
    *leniw = 23;

    /* add NVector sizes */
    *lenrw += ark_mem.tempv1.data.len() as i64;
    *leniw += 1; /* N_VSpace_Serial: (n, 1) */

    /* add SUNMatrix size (only account for the one owned by Ls interface) */
    if !arkls_mem.iterative {
        if let Some(m_lu) = &arkls_mem.M_lu {
            let (mut lrw, mut liw) = (0i64, 0i64);
            if SUNMatSpace(m_lu, &mut lrw, &mut liw) == 0 {
                *lenrw += lrw;
                *leniw += liw;
            }
        }
    }

    /* add LS sizes */
    let (lrw, liw) = arkls_mem.LS.space();
    *lenrw += lrw;
    *leniw += liw;

    ark_mem.mass_mem = Some(arkls_mem);
    ARKLS_SUCCESS
}

/* Shared boilerplate for the scalar mass statistics getters. */
macro_rules! arkls_mass_stat_get {
    ($name:ident, $field:ident) => {
        pub fn $name(ark_mem: &mut ARKodeMem, out: &mut i64) -> i32 {
            /* Return 0 for incompatible steppers */
            if !ark_mem.step_supports_massmatrix {
                *out = 0;
                return ARK_SUCCESS;
            }

            let arkls_mem = match arkLs_AccessMassMem(ark_mem, stringify!($name)) {
                Ok(m) => m,
                Err(e) => return e,
            };

            *out = arkls_mem.$field;
            ark_mem.mass_mem = Some(arkls_mem);
            ARKLS_SUCCESS
        }
    };
}

arkls_mass_stat_get!(ARKodeGetNumMassSetups, nmsetups);
arkls_mass_stat_get!(ARKodeGetNumMassMult, nmtimes);
arkls_mass_stat_get!(ARKodeGetNumMassSolves, nmsolves);
arkls_mass_stat_get!(ARKodeGetNumMassPrecEvals, npe);
arkls_mass_stat_get!(ARKodeGetNumMassPrecSolves, nps);
arkls_mass_stat_get!(ARKodeGetNumMassIters, nli);
arkls_mass_stat_get!(ARKodeGetNumMassConvFails, ncfl);
arkls_mass_stat_get!(ARKodeGetNumMTSetups, nmtsetup);
arkls_mass_stat_get!(ARKodeGetNumMassMultSetups, nmvsetup);

/*---------------------------------------------------------------
  ARKodeGetCurrentMassMatrix returns a copy of the current mass
  matrix (C hands out the internal pointer).
  ---------------------------------------------------------------*/
pub fn ARKodeGetCurrentMassMatrix(ark_mem: &mut ARKodeMem, M: &mut Option<SUNMatrix>) -> i32 {
    /* Return NULL for incompatible steppers */
    if !ark_mem.step_supports_massmatrix {
        *M = None;
        return ARK_SUCCESS;
    }

    let arkls_mem = match arkLs_AccessMassMem(ark_mem, "ARKodeGetCurrentMassMatrix") {
        Ok(m) => m,
        Err(e) => return e,
    };

    /* set output and return */
    *M = arkls_mem.M.as_ref().map(SUNMatClone_Copy);
    ark_mem.mass_mem = Some(arkls_mem);
    ARKLS_SUCCESS
}

pub fn ARKodeGetLastMassFlag(ark_mem: &mut ARKodeMem, flag: &mut i64) -> i32 {
    /* Return 0 for incompatible steppers */
    if !ark_mem.step_supports_massmatrix {
        *flag = 0;
        return ARK_SUCCESS;
    }

    let arkls_mem = match arkLs_AccessMassMem(ark_mem, "ARKodeGetLastMassFlag") {
        Ok(m) => m,
        Err(e) => return e,
    };

    *flag = arkls_mem.last_flag as i64;
    ark_mem.mass_mem = Some(arkls_mem);
    ARKLS_SUCCESS
}

/*---------------------------------------------------------------
  arkLSSetMassUserData sets user_data pointers in the mass ARKLS
  interface (Rust: the data pointers collapsed onto
  ark_mem.user_data, so there is nothing to update).
  ---------------------------------------------------------------*/
pub fn arkLSSetMassUserData(_ark_mem: &mut ARKodeMem) -> i32 {
    ARKLS_SUCCESS
}

/*---------------------------------------------------------------
  arkLsMTimes:

  This routine generates the matrix-vector product z = Mv, where
  M is the system mass matrix, by calling the user-supplied mtimes
  routine (default) or asking the SUNMatrix to do the multiply.
  ---------------------------------------------------------------*/
pub fn arkLsMTimes(ark_mem: &mut ARKodeMem, v: &NVector, z: &mut NVector) -> i32 {
    let mut arkls_mem = match arkLs_AccessMassMem(ark_mem, "arkLsMTimes") {
        Ok(m) => m,
        Err(e) => return e,
    };
    let ret = arkLsMTimes_inner(ark_mem, &mut arkls_mem, v, z);
    ark_mem.mass_mem = Some(arkls_mem);
    ret
}

fn arkLsMTimes_inner(
    ark_mem: &mut ARKodeMem,
    arkls_mem: &mut crate::arkode_ls_impl::ARKLsMassMem,
    v: &NVector,
    z: &mut NVector,
) -> i32 {
    /* perform multiply by either calling the user-supplied routine
       (default), or asking the SUNMatrix to do the multiply */
    if let Some(mtimes) = arkls_mem.mtimes {
        /* call user-supplied mtimes routine, increment counter and return */
        let retval = mtimes(v, z, ark_mem.tcur, &mut ark_mem.user_data);
        if retval == 0 {
            arkls_mem.nmtimes += 1;
        } else {
            arkProcessError(
                Some(ark_mem),
                retval,
                line!(),
                "arkLsMTimes",
                file!(),
                "Error in user mass matrix-vector product routine",
            );
        }
        return retval;
    } else if let Some(m) = &arkls_mem.M {
        /* ask SUNMatrix to do the multiply; increment counter and return */
        let retval = SUNMatMatvec(m, v, z);
        if retval == 0 {
            arkls_mem.nmtimes += 1;
        } else {
            arkProcessError(
                Some(ark_mem),
                retval,
                line!(),
                "arkLsMTimes",
                file!(),
                "Error in SUNMatrix mass matrix-vector product routine",
            );
        }
        return retval;
    }

    /* if we made it here, then no matrix-vector product is available */
    arkProcessError(
        Some(ark_mem),
        -1,
        line!(),
        "arkLsMTimes",
        file!(),
        "Missing mass matrix-vector product routine",
    );
    -1
}

/*---------------------------------------------------------------
  arkLsMassInitialize performs remaining initializations specific
  to the mass matrix solver interface (and solver itself)
  ---------------------------------------------------------------*/
pub fn arkLsMassInitialize(ark_mem: &mut ARKodeMem) -> i32 {
    let mut arkls_mem = match arkLs_AccessMassMem(ark_mem, "arkLsMassInitialize") {
        Ok(m) => m,
        Err(e) => return e,
    };
    let ret = arkLsMassInitialize_inner(ark_mem, &mut arkls_mem);
    ark_mem.mass_mem = Some(arkls_mem);
    ret
}

fn arkLsMassInitialize_inner(
    ark_mem: &mut ARKodeMem,
    arkls_mem: &mut crate::arkode_ls_impl::ARKLsMassMem,
) -> i32 {
    /* reset counters */
    arkLsInitializeMassCounters(arkls_mem);

    /* perform checks for matrix-based mass system */
    if arkls_mem.M.is_some() {
        /* check for user-provided mass matrix constructor */
        if arkls_mem.mass.is_none() {
            arkProcessError(
                Some(ark_mem),
                ARKLS_ILL_INPUT,
                line!(),
                "arkLsMassInitialize",
                file!(),
                "Missing user-provided mass-matrix routine",
            );
            arkls_mem.last_flag = ARKLS_ILL_INPUT;
            return arkls_mem.last_flag;
        }
        /* (matrix-vector products are always available: the serial
           SUNMatrix implementations all provide matvec) */
    }

    /* perform checks for matrix-free mass system */
    if arkls_mem.M.is_none()
        && arkls_mem.mtimes.is_none()
        && arkls_mem.LS.ls_type() != SUNLINEARSOLVER_MATRIX_EMBEDDED
    {
        arkProcessError(
            Some(ark_mem),
            ARKLS_ILL_INPUT,
            line!(),
            "arkLsMassInitialize",
            file!(),
            "Missing user-provided mass matrix-vector product routine",
        );
        arkls_mem.last_flag = ARKLS_ILL_INPUT;
        return arkls_mem.last_flag;
    }

    /* if M is NULL and neither pset or mtsetup are present, then
       arkLsMassSetup does not need to be called, so set the
       msetup function to NULL */
    if arkls_mem.M.is_none() && arkls_mem.pset.is_none() && arkls_mem.mtsetup.is_none() {
        if let Some(disable) = ark_mem.step_disablemsetup {
            disable(ark_mem);
        }
    }

    /* When using a matrix-embedded linear solver, disable lsetup call */
    if arkls_mem.LS.ls_type() == SUNLINEARSOLVER_MATRIX_EMBEDDED {
        if let Some(disable) = ark_mem.step_disablemsetup {
            disable(ark_mem);
        }
    }

    /* Call LS initialize routine */
    arkls_mem.last_flag = arkls_mem.LS.initialize();
    arkls_mem.last_flag
}

/*---------------------------------------------------------------
  arkLsMassSetup calls the LS 'setup' routine.
  ---------------------------------------------------------------*/
pub fn arkLsMassSetup(
    ark_mem: &mut ARKodeMem,
    t: f64,
    vtemp1: &mut NVector,
    vtemp2: &mut NVector,
    vtemp3: &mut NVector,
) -> i32 {
    let mut arkls_mem = match arkLs_AccessMassMem(ark_mem, "arkLsMassSetup") {
        Ok(m) => m,
        Err(e) => return e,
    };
    let ret = arkLsMassSetup_inner(ark_mem, &mut arkls_mem, t, vtemp1, vtemp2, vtemp3);
    ark_mem.mass_mem = Some(arkls_mem);
    ret
}

fn arkLsMassSetup_inner(
    ark_mem: &mut ARKodeMem,
    arkls_mem: &mut crate::arkode_ls_impl::ARKLsMassMem,
    t: f64,
    vtemp1: &mut NVector,
    vtemp2: &mut NVector,
    vtemp3: &mut NVector,
) -> i32 {
    /* Immediately return when using matrix-embedded linear solver */
    if arkls_mem.LS.ls_type() == SUNLINEARSOLVER_MATRIX_EMBEDDED {
        arkls_mem.last_flag = ARKLS_SUCCESS;
        return arkls_mem.last_flag;
    }

    /* if the most recent setup essentially matches the current time,
       just return with success */
    if SUNRabs(arkls_mem.msetuptime - t) < FUZZ_FACTOR * ark_mem.uround {
        arkls_mem.last_flag = ARKLS_SUCCESS;
        return arkls_mem.last_flag;
    }

    /* Determine whether to call user-provided mtsetup routine */
    let call_mtsetup = arkls_mem.mtsetup.is_some()
        && (arkls_mem.time_dependent || arkls_mem.nmtsetup == 0);

    /* call user-provided mtsetup routine if applicable */
    if call_mtsetup {
        let mtsetup = arkls_mem.mtsetup.unwrap();
        arkls_mem.last_flag = mtsetup(t, &mut ark_mem.user_data);
        arkls_mem.nmtsetup += 1;
        arkls_mem.msetuptime = t;
        if arkls_mem.last_flag != 0 {
            arkProcessError(
                Some(ark_mem),
                arkls_mem.last_flag,
                line!(),
                "arkLsMassSetup",
                file!(),
                "The mass matrix x vector setup routine failed in an unrecoverable manner.",
            );
            return arkls_mem.last_flag;
        }
    }

    /* Perform user-facing setup based on whether this is matrix-free */
    let call_lssetup;
    if arkls_mem.M.is_none() {
        /*** matrix-free -- only call LS setup if preconditioner setup
             exists (for the Rust iterative solvers that setup is the
             MPSetup call below) ***/
        call_lssetup = arkls_mem.pset.is_some();
    } else {
        /*** matrix-based ***/

        /* If mass matrix is not time dependent, and if it has been set up
           previously, then just reuse existing matrix and factorization */
        if !arkls_mem.time_dependent && arkls_mem.nmsetups > 0 {
            arkls_mem.last_flag = ARKLS_SUCCESS;
            return arkls_mem.last_flag;
        }

        /* Clear the mass matrix if necessary (direct linear solvers) */
        if !arkls_mem.iterative {
            let retval = SUNMatZero(arkls_mem.M.as_mut().unwrap());
            if retval != 0 {
                arkProcessError(
                    Some(ark_mem),
                    ARKLS_SUNMAT_FAIL,
                    line!(),
                    "arkLsMassSetup",
                    file!(),
                    MSG_LS_SUNMAT_FAILED,
                );
                arkls_mem.last_flag = ARKLS_SUNMAT_FAIL;
                return arkls_mem.last_flag;
            }
        }

        /* Call user-supplied routine to fill the mass matrix */
        let mass = arkls_mem.mass.unwrap();
        let retval = mass(
            t,
            arkls_mem.M.as_mut().unwrap(),
            &mut ark_mem.user_data,
            vtemp1,
            vtemp2,
            vtemp3,
        );
        arkls_mem.msetuptime = t;
        if retval < 0 {
            arkProcessError(
                Some(ark_mem),
                ARKLS_MASSFUNC_UNRECVR,
                line!(),
                "arkLsMassSetup",
                file!(),
                "The mass matrix routine failed in an unrecoverable manner.",
            );
            arkls_mem.last_flag = ARKLS_MASSFUNC_UNRECVR;
            return -1;
        }
        if retval > 0 {
            arkls_mem.last_flag = ARKLS_MASSFUNC_RECVR;
            return 1;
        }

        /* Copy M into M_lu for factorization (direct linear solvers) */
        if !arkls_mem.iterative {
            let crate::arkode_ls_impl::ARKLsMassMem { M, M_lu, .. } = arkls_mem;
            let retval = SUNMatCopy(M.as_ref().unwrap(), M_lu.as_mut().unwrap());
            if retval != 0 {
                arkProcessError(
                    Some(ark_mem),
                    ARKLS_SUNMAT_FAIL,
                    line!(),
                    "arkLsMassSetup",
                    file!(),
                    MSG_LS_SUNMAT_FAILED,
                );
                arkls_mem.last_flag = ARKLS_SUNMAT_FAIL;
                return arkls_mem.last_flag;
            }
        }

        /* (matvec setup: the serial SUNMatrix implementations have no
           matvecsetup routine, matching C where call_mvsetup stays
           SUNFALSE for them) */

        /* signal call to LS setup routine */
        call_lssetup = true;
    }

    /* Call LS setup routine if applicable, and return */
    if call_lssetup {
        arkls_mem.last_flag = match &arkls_mem.LS {
            LinearSolver::Dense(_) | LinearSolver::Band(_) => {
                let crate::arkode_ls_impl::ARKLsMassMem { LS, M_lu, .. } = arkls_mem;
                LS.setup(M_lu.as_mut())
            }
            _ => {
                /* iterative solver: mass preconditioner setup (arkLsMPSetup);
                   only proceed if the mass matrix is time-dependent or if
                   pset has not been called previously */
                if let Some(pset) = arkls_mem.pset {
                    if arkls_mem.time_dependent || arkls_mem.npe == 0 {
                        let r = pset(ark_mem.tcur, &mut ark_mem.user_data);
                        arkls_mem.npe += 1;
                        r
                    } else {
                        0
                    }
                } else {
                    0
                }
            }
        };
        arkls_mem.nmsetups += 1;
    }

    arkls_mem.last_flag
}

/*---------------------------------------------------------------
  arkLsMassSolve: interfaces between ARKODE and the generic
  SUNLinearSolver object LS, by setting the appropriate tolerance
  and scaling vectors, calling the solver, and accumulating
  statistics from the solve for use/reporting by ARKODE.
  ---------------------------------------------------------------*/
pub fn arkLsMassSolve(ark_mem: &mut ARKodeMem, b: &mut NVector, nlscoef: f64) -> i32 {
    let mut arkls_mem = match arkLs_AccessMassMem(ark_mem, "arkLsMassSolve") {
        Ok(m) => m,
        Err(e) => return e,
    };
    let ret = arkLsMassSolve_inner(ark_mem, &mut arkls_mem, b, nlscoef);
    ark_mem.mass_mem = Some(arkls_mem);
    ret
}

fn arkLsMassSolve_inner(
    ark_mem: &mut ARKodeMem,
    arkls_mem: &mut crate::arkode_ls_impl::ARKLsMassMem,
    b: &mut NVector,
    nlscoef: f64,
) -> i32 {
    /* Set input tolerance for iterative solvers (in 2-norm) */
    let delta = if arkls_mem.iterative {
        arkls_mem.eplifac * nlscoef * arkls_mem.nrmfac
    } else {
        ZERO
    };

    /* (Scaling vectors: the Rust iterative solvers all accept the s1/s2
       scaling vectors at solve time — the C "solver does not support
       scaling vectors" rwt_mean fallback branch is unreachable.) */

    /* Set initial guess x = 0 for LS, and zero initial guess flag */
    N_VConst(ZERO, &mut arkls_mem.x);
    arkls_mem.LS.set_zero_guess(SUNTRUE);

    /* Call solver, copy x to b, and increment mass solver counter */
    let retval = if arkls_mem.LS.ls_type() == SUNLINEARSOLVER_DIRECT {
        let crate::arkode_ls_impl::ARKLsMassMem { LS, M_lu, x, .. } = arkls_mem;
        match (LS, M_lu.as_mut()) {
            (LinearSolver::Dense(dls), Some(SUNMatrix::Dense(am))) => dls.solve(am, x, b),
            (LinearSolver::Band(bls), Some(SUNMatrix::Band(am))) => bls.solve(am, x, b),
            _ => SUN_ERR_ARG_INCOMPATIBLE,
        }
    } else {
        arkLsMassSolveIterative(ark_mem, arkls_mem, b, delta)
    };
    {
        let crate::arkode_ls_impl::ARKLsMassMem { x, .. } = arkls_mem;
        N_VScale(ONE, x, b);
    }
    arkls_mem.nmsolves += 1;

    /* Retrieve statistics from iterative linear solvers */
    let mut nli_inc = 0;
    if arkls_mem.iterative {
        nli_inc = arkls_mem.LS.num_iters();
    }

    /* Increment counters nli and ncfl */
    arkls_mem.nli += nli_inc as i64;
    if retval != SUN_SUCCESS {
        arkls_mem.ncfl += 1;
    }

    /* Interpret solver return value  */
    arkls_mem.last_flag = retval;

    match retval {
        SUN_SUCCESS => 0,
        SUNLS_RES_REDUCED | SUNLS_CONV_FAIL | SUNLS_ATIMES_FAIL_REC | SUNLS_PSOLVE_FAIL_REC
        | SUNLS_PACKAGE_FAIL_REC | SUNLS_QRFACT_FAIL | SUNLS_LUFACT_FAIL => 1,
        SUN_ERR_ARG_CORRUPT | SUN_ERR_ARG_INCOMPATIBLE | SUN_ERR_MEM_FAIL | SUNLS_GS_FAIL
        | SUNLS_QRSOL_FAIL => -1,
        SUN_ERR_EXT_FAIL => {
            arkProcessError(
                Some(ark_mem),
                SUN_ERR_EXT_FAIL,
                line!(),
                "arkLsMassSolve",
                file!(),
                "Failure in SUNLinSol external package",
            );
            -1
        }
        SUNLS_ATIMES_FAIL_UNREC => {
            arkProcessError(
                Some(ark_mem),
                SUNLS_ATIMES_FAIL_UNREC,
                line!(),
                "arkLsMassSolve",
                file!(),
                "The mass matrix x vector routine failed in an unrecoverable manner.",
            );
            -1
        }
        SUNLS_PSOLVE_FAIL_UNREC => {
            arkProcessError(
                Some(ark_mem),
                SUNLS_PSOLVE_FAIL_UNREC,
                line!(),
                "arkLsMassSolve",
                file!(),
                MSG_LS_PSOLVE_FAILED,
            );
            -1
        }
        _ => 0,
    }
}

/* Iterative mass solve: builds the MTimes (arkLsMTimes) and MPSolve
   (arkLsMPSolve) callbacks over the integrator memory (same RefCell
   pattern as arkLsSolveIterative); scaling vectors s1 = rwt,
   s2 = ewt. */
fn arkLsMassSolveIterative(
    ark_mem: &mut ARKodeMem,
    arkls_mem: &mut crate::arkode_ls_impl::ARKLsMassMem,
    b: &NVector,
    delta: f64,
) -> i32 {
    let crate::arkode_ls_impl::ARKLsMassMem {
        LS,
        x,
        nmtimes,
        nps,
        mtimes,
        psolve,
        M,
        ..
    } = arkls_mem;
    let mtimes_fn = *mtimes;
    let psolve_fn = *psolve;

    let ewt = std::mem::take(&mut ark_mem.ewt);
    let rwt_detached = if ark_mem.rwt_is_ewt {
        None
    } else {
        Some(std::mem::take(&mut ark_mem.rwt))
    };
    let arm = RefCell::new(&mut *ark_mem);
    let ewt_ref = &ewt;
    let rwt_ref: &NVector = match &rwt_detached {
        Some(r) => r,
        None => ewt_ref,
    };

    let mut atimes = |v: &NVector, z: &mut NVector| -> i32 {
        let mut guard = arm.borrow_mut();
        let ar: &mut ARKodeMem = &mut guard;

        /* arkLsMTimes: user mtimes (default) or SUNMatrix multiply */
        if let Some(mt) = mtimes_fn {
            let ret = mt(v, z, ar.tcur, &mut ar.user_data);
            if ret == 0 {
                *nmtimes += 1;
            }
            ret
        } else if let Some(m) = M.as_ref() {
            let ret = SUNMatMatvec(m, v, z);
            if ret == 0 {
                *nmtimes += 1;
            }
            ret
        } else {
            -1
        }
    };

    let mut psolve_cb = |r: &NVector, z: &mut NVector, tol: f64, lr: i32| -> i32 {
        let mut guard = arm.borrow_mut();
        let ar: &mut ARKodeMem = &mut guard;

        /* arkLsMPSolve: call the user-supplied psolve routine */
        let ret = if let Some(ps) = psolve_fn {
            ps(ar.tcur, r, z, tol, lr, &mut ar.user_data)
        } else {
            0
        };
        *nps += 1;
        ret
    };

    let retval = if psolve_fn.is_some() {
        LS.solve(
            None,
            x,
            b,
            delta,
            &mut atimes,
            Some(&mut psolve_cb),
            Some(rwt_ref),
            Some(ewt_ref),
        )
    } else {
        LS.solve(None, x, b, delta, &mut atimes, None, Some(rwt_ref), Some(ewt_ref))
    };

    /* end the closures' shared borrow of the RefCell before unwrapping */
    let _ = (atimes, psolve_cb);
    let ark_mem_back = arm.into_inner();
    ark_mem_back.ewt = ewt;
    if let Some(r) = rwt_detached {
        ark_mem_back.rwt = r;
    }

    retval
}

/*---------------------------------------------------------------
  arkLsMassFree frees memory associated with the ARKLs mass
  matrix solver interface.
  ---------------------------------------------------------------*/
pub fn arkLsMassFree(ark_mem: &mut ARKodeMem) -> i32 {
    /* Return immediately if ARKodeMem, ARKLsMassMem are NULL */
    let taken = match ark_mem.step_getmassmem {
        Some(get) => get(ark_mem),
        None => None,
    };
    let mut arkls_mem = match taken {
        Some(m) => m,
        None => return ARKLS_SUCCESS,
    };

    /* Free N_Vector memory (with lrw/liw accounting) */
    arkFreeVec(ark_mem, &mut arkls_mem.x);

    /* M_lu / M / LS memory is dropped with the box (C destroys M_lu for
       direct solvers, nullifies the borrowed pointers, and calls any
       pfree here) */
    ARKLS_SUCCESS
}

/*---------------------------------------------------------------
  arkLsInitializeMassCounters resets all counters from an
  ARKLsMassMem structure.
  ---------------------------------------------------------------*/
pub fn arkLsInitializeMassCounters(arkls_mem: &mut crate::arkode_ls_impl::ARKLsMassMem) -> i32 {
    arkls_mem.nmsetups = 0;
    arkls_mem.nmsolves = 0;
    arkls_mem.nmtsetup = 0;
    arkls_mem.nmtimes = 0;
    arkls_mem.nmvsetup = 0;
    arkls_mem.npe = 0;
    arkls_mem.nli = 0;
    arkls_mem.nps = 0;
    arkls_mem.ncfl = 0;
    arkls_mem.msetuptime = -crate::sundials_types::SUN_BIG_REAL;
    0
}
