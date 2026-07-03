/* -----------------------------------------------------------------
 * Translated from src/cvodes/cvodes_io.c (CVODES 7.7.0).
 * Optional input and output functions for the CVODES solver.
 *
 * The C functions take `void* cvode_mem` and start with a NULL
 * check; here the memory is `&mut CVodeMem`, which cannot be null,
 * so those checks vanish. All other checks, defaults and messages
 * are translated line-for-line. Function order follows the C file.
 *
 * Documented adaptations (donor precedent: cvode_rs/src/cvode_io.rs):
 *  - CVodeSetMonitorFn / CVodeSetMonitorFrequency: monitoring is not
 *    enabled in this build (C build without SUNDIALS_ENABLE_MONITORING,
 *    and cvodes_impl.rs omits cv_monitorfun/cv_monitor_interval), so
 *    both setters take the C #else branch and return CV_ILL_INPUT,
 *    exactly as the donor's cvode_io.rs does. CVMonitorFn is defined
 *    here only to keep the CVodeSetMonitorFn signature.
 *  - CVodeSetConstraints: the C NULL-pointer tests on the
 *    `constraints` argument and on cv_mem->cv_constraints become
 *    Option<&NVector> and the empty-NVector state (empty = C NULL,
 *    cvodes_impl.rs pinned decision 1). cvodes_impl.h has no
 *    cv_constraintsSet flag (unlike the CVODE donor), so the port
 *    tests cv_constraints.is_empty() directly, mirroring the C
 *    pointer test. The C vector-ops check is dropped: the serial
 *    NVector implements every required op.
 *  - CVodeSetSensParams: C stores the user's `p` POINTER
 *    (cv_mem->cv_p = p); cv_p is an owned Vec here, so the port
 *    copies the slice (None -> empty Vec = C NULL). pbar/plist are
 *    likewise Option<&[..]> and are copied element-wise as in C.
 *  - CVodeSetMaxNonlinIters / CVodeSetSensMaxNonlinIters: the C
 *    SUNNonlinSolSetMaxIters(NLS*, ...) calls become the
 *    NonlinearSolver::set_max_iters method on the Option fields.
 *  - CVodeGetCurrentState / CVodeGetCurrentStateSens return
 *    references (&NVector / &[NVector]) instead of writing through
 *    N_Vector* / N_Vector** out-parameters (donor pattern).
 *  - CVodeGetUserData returns &mut UserData (donor pattern).
 *  - CVodePrintAllStats: FILE* -> &mut dyn std::io::Write; the
 *    sunfprintf_real/sunfprintf_long/sunfprintf_long_array helpers
 *    from src/sundials/sundials_utils.h are ported below using
 *    sundials_utils::fmt_g / fmt_e (SUN_FORMAT_G = "%.15g",
 *    SUN_FORMAT_E = "% .15e").
 *  - CVodeGetReturnFlagName returns String (C mallocs a char*).
 * -----------------------------------------------------------------*/
use crate::cvodes_impl::*;
use crate::nvector_serial::{NVector, N_VMaxNorm, N_VScale};
use crate::sundials_math::SUNRabs;
use crate::sundials_types::*;
use crate::sundials_utils::{fmt_e, fmt_g};

const ZERO: f64 = 0.0;
const HALF: f64 = 0.5;
const ONE: f64 = 1.0;
const TWOPT5: f64 = 2.5;

/* CVMonitorFn (cvodes.h): user function for monitoring integrator
   progress. Monitoring is not enabled in this build (the C build
   without SUNDIALS_ENABLE_MONITORING), so the type is only used to
   keep the CVodeSetMonitorFn signature. */
pub type CVMonitorFn = fn(cv_mem: &CVodeMem, user_data: &mut UserData);

/*
 * =================================================================
 * CVODES optional input functions
 * =================================================================
 */

/*
 * CVodeSetDeltaGammaMaxLSetup
 *
 * Specifies the gamma ratio threshold to signal for a linear solver setup
 */

pub fn CVodeSetDeltaGammaMaxLSetup(cv_mem: &mut CVodeMem, dgmax_lsetup: f64) -> i32 {
    /* Set value or use default */
    if dgmax_lsetup < ZERO {
        cv_mem.cv_dgmax_lsetup = DGMAX_LSETUP_DEFAULT;
    } else {
        cv_mem.cv_dgmax_lsetup = dgmax_lsetup;
    }

    CV_SUCCESS
}

/*
 * CVodeSetUserData
 *
 * Specifies the user data pointer for f
 */

pub fn CVodeSetUserData(cv_mem: &mut CVodeMem, user_data: UserData) -> i32 {
    cv_mem.cv_user_data = user_data;

    CV_SUCCESS
}

/*
 * CVodeSetMonitorFn
 *
 * Specifies the user function to call for monitoring
 * the solution and/or integrator statistics.
 *
 * This build corresponds to a C build without
 * SUNDIALS_ENABLE_MONITORING.
 */

pub fn CVodeSetMonitorFn(cv_mem: &mut CVodeMem, _fn_: CVMonitorFn) -> i32 {
    cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSetMonitorFn", file!(),
                   "SUNDIALS was not built with monitoring enabled.");
    CV_ILL_INPUT
}

/*
 * CVodeSetMonitorFrequency
 *
 * Specifies the frequency with which to call the user function.
 */

pub fn CVodeSetMonitorFrequency(cv_mem: &mut CVodeMem, nst: i64) -> i32 {
    if nst < 0 {
        cvProcessError(None, CV_ILL_INPUT, line!(), "CVodeSetMonitorFrequency", file!(),
                       "step interval must be >= 0\n");
        return CV_ILL_INPUT;
    }

    cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSetMonitorFrequency", file!(),
                   "SUNDIALS was not built with monitoring enabled.");
    CV_ILL_INPUT
}

/*
 * CVodeSetMaxOrd
 *
 * Specifies the maximum method order
 */

pub fn CVodeSetMaxOrd(cv_mem: &mut CVodeMem, maxord: i32) -> i32 {
    if maxord <= 0 {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSetMaxOrd", file!(),
                       MSGCV_NEG_MAXORD);
        return CV_ILL_INPUT;
    }

    /* Cannot increase maximum order beyond the value that
       was used when allocating memory */
    let mut qmax_alloc = cv_mem.cv_qmax_alloc;
    qmax_alloc = qmax_alloc.min(cv_mem.cv_qmax_allocQ);
    qmax_alloc = qmax_alloc.min(cv_mem.cv_qmax_allocS);

    if maxord > qmax_alloc {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSetMaxOrd", file!(),
                       MSGCV_BAD_MAXORD);
        return CV_ILL_INPUT;
    }

    cv_mem.cv_qmax = maxord;

    CV_SUCCESS
}

/*
 * CVodeSetMaxNumSteps
 *
 * Specifies the maximum number of integration steps
 */

pub fn CVodeSetMaxNumSteps(cv_mem: &mut CVodeMem, mxsteps: i64) -> i32 {
    /* Passing mxsteps=0 sets the default. Passing mxsteps<0 disables the test. */
    if mxsteps == 0 {
        cv_mem.cv_mxstep = MXSTEP_DEFAULT;
    } else {
        cv_mem.cv_mxstep = mxsteps;
    }

    CV_SUCCESS
}

/*
 * CVodeSetMaxHnilWarns
 *
 * Specifies the maximum number of warnings for small h
 */

pub fn CVodeSetMaxHnilWarns(cv_mem: &mut CVodeMem, mxhnil: i32) -> i32 {
    cv_mem.cv_mxhnil = mxhnil;

    CV_SUCCESS
}

/*
 * CVodeSetStabLimDet
 *
 * Turns on/off the stability limit detection algorithm
 */

pub fn CVodeSetStabLimDet(cv_mem: &mut CVodeMem, sldet: bool) -> i32 {
    if sldet && (cv_mem.cv_lmm != CV_BDF) {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSetStabLimDet", file!(),
                       MSGCV_SET_SLDET);
        return CV_ILL_INPUT;
    }

    cv_mem.cv_sldeton = sldet;

    CV_SUCCESS
}

/*
 * CVodeSetInitStep
 *
 * Specifies the initial step size
 */

pub fn CVodeSetInitStep(cv_mem: &mut CVodeMem, hin: f64) -> i32 {
    cv_mem.cv_hin = hin;

    CV_SUCCESS
}

/*
 * CVodeSetMinStep
 *
 * Specifies the minimum step size
 */

pub fn CVodeSetMinStep(cv_mem: &mut CVodeMem, hmin: f64) -> i32 {
    if hmin < ZERO {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSetMinStep", file!(),
                       MSGCV_NEG_HMIN);
        return CV_ILL_INPUT;
    }

    /* Passing 0 sets hmin = zero */
    if hmin == ZERO {
        cv_mem.cv_hmin = HMIN_DEFAULT;
        return CV_SUCCESS;
    }

    if hmin * cv_mem.cv_hmax_inv > ONE {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSetMinStep", file!(),
                       MSGCV_BAD_HMIN_HMAX);
        return CV_ILL_INPUT;
    }

    cv_mem.cv_hmin = hmin;

    CV_SUCCESS
}

/*
 * CVodeSetMaxStep
 *
 * Specifies the maximum step size
 */

pub fn CVodeSetMaxStep(cv_mem: &mut CVodeMem, hmax: f64) -> i32 {
    if hmax < ZERO {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSetMaxStep", file!(),
                       MSGCV_NEG_HMAX);
        return CV_ILL_INPUT;
    }

    /* Passing 0 sets hmax = infinity */
    if hmax == ZERO {
        cv_mem.cv_hmax_inv = HMAX_INV_DEFAULT;
        return CV_SUCCESS;
    }

    let hmax_inv = ONE / hmax;
    if hmax_inv * cv_mem.cv_hmin > ONE {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSetMaxStep", file!(),
                       MSGCV_BAD_HMIN_HMAX);
        return CV_ILL_INPUT;
    }

    cv_mem.cv_hmax_inv = hmax_inv;

    CV_SUCCESS
}

/*
 * CVodeSetEtaFixedStepBounds
 *
 * Specifies the bounds for retaining the current step size
 */

pub fn CVodeSetEtaFixedStepBounds(cv_mem: &mut CVodeMem, eta_min_fx: f64, eta_max_fx: f64) -> i32 {
    /* set allowed value or use default */
    if eta_min_fx >= ZERO && eta_min_fx <= ONE {
        cv_mem.cv_eta_min_fx = eta_min_fx;
    } else {
        cv_mem.cv_eta_min_fx = ETA_MIN_FX_DEFAULT;
    }

    if eta_max_fx >= ONE {
        cv_mem.cv_eta_max_fx = eta_max_fx;
    } else {
        cv_mem.cv_eta_max_fx = ETA_MAX_FX_DEFAULT;
    }

    CV_SUCCESS
}

/*
 * CVodeSetEtaMaxFirstStep
 *
 * Specifies the maximum step size change on the first step
 */

pub fn CVodeSetEtaMaxFirstStep(cv_mem: &mut CVodeMem, eta_max_fs: f64) -> i32 {
    /* set allowed value or use default */
    if eta_max_fs <= ONE {
        cv_mem.cv_eta_max_fs = ETA_MAX_FS_DEFAULT;
    } else {
        cv_mem.cv_eta_max_fs = eta_max_fs;
    }

    CV_SUCCESS
}

/*
 * CVodeSetEtaMaxEarlyStep
 *
 * Specifies the maximum step size change on steps early in the integration
 * when nst <= small_nst
 */

pub fn CVodeSetEtaMaxEarlyStep(cv_mem: &mut CVodeMem, eta_max_es: f64) -> i32 {
    /* set allowed value or use default */
    if eta_max_es <= ONE {
        cv_mem.cv_eta_max_es = ETA_MAX_ES_DEFAULT;
    } else {
        cv_mem.cv_eta_max_es = eta_max_es;
    }

    CV_SUCCESS
}

/*
 * CVodeSetNumStepsEtaMaxEarlyStep
 *
 * Specifies the maximum number of steps for using the early integration change
 * factor
 */

pub fn CVodeSetNumStepsEtaMaxEarlyStep(cv_mem: &mut CVodeMem, small_nst: i64) -> i32 {
    /* set allowed value or use default */
    if small_nst < 0 {
        cv_mem.cv_small_nst = SMALL_NST_DEFAULT;
    } else {
        cv_mem.cv_small_nst = small_nst;
    }

    CV_SUCCESS
}

/*
 * CVodeSetEtaMax
 *
 * Specifies the maximum step size change on a general steps (nst > small_nst)
 */

pub fn CVodeSetEtaMax(cv_mem: &mut CVodeMem, eta_max_gs: f64) -> i32 {
    /* set allowed value or use default */
    if eta_max_gs <= ONE {
        cv_mem.cv_eta_max_gs = ETA_MAX_GS_DEFAULT;
    } else {
        cv_mem.cv_eta_max_gs = eta_max_gs;
    }

    CV_SUCCESS
}

/*
 * CVodeSetEtaMin
 *
 * Specifies the minimum change on a general steps
 */

pub fn CVodeSetEtaMin(cv_mem: &mut CVodeMem, eta_min: f64) -> i32 {
    /* set allowed value or use default */
    if eta_min <= ZERO || eta_min >= ONE {
        cv_mem.cv_eta_min = ETA_MIN_DEFAULT;
    } else {
        cv_mem.cv_eta_min = eta_min;
    }

    CV_SUCCESS
}

/*
 * CVodeSetEtaMinErrFail
 *
 * Specifies the minimum step size change after an error test failure
 */

pub fn CVodeSetEtaMinErrFail(cv_mem: &mut CVodeMem, eta_min_ef: f64) -> i32 {
    /* set allowed value or use default */
    if eta_min_ef <= ZERO || eta_min_ef >= ONE {
        cv_mem.cv_eta_min_ef = ETA_MIN_EF_DEFAULT;
    } else {
        cv_mem.cv_eta_min_ef = eta_min_ef;
    }

    CV_SUCCESS
}

/*
 * CVodeSetEtaMaxErrFail
 *
 * Specifies the maximum step size change after multiple (>= small_nef) error
 * test failures
 */

pub fn CVodeSetEtaMaxErrFail(cv_mem: &mut CVodeMem, eta_max_ef: f64) -> i32 {
    /* set allowed value or use default */
    if eta_max_ef <= ZERO || eta_max_ef >= ONE {
        cv_mem.cv_eta_max_ef = ETA_MAX_EF_DEFAULT;
    } else {
        cv_mem.cv_eta_max_ef = eta_max_ef;
    }

    CV_SUCCESS
}

/*
 * CVodeSetNumFailsEtaMaxErrFail
 *
 * Specifies the maximum number of error test failures necessary to enforce
 * eta_max_ef
 */

pub fn CVodeSetNumFailsEtaMaxErrFail(cv_mem: &mut CVodeMem, small_nef: i32) -> i32 {
    /* set allowed value or use default */
    if small_nef < 0 {
        cv_mem.cv_small_nef = SMALL_NEF_DEFAULT;
    } else {
        cv_mem.cv_small_nef = small_nef;
    }

    CV_SUCCESS
}

/*
 * CVodeSetEtaConvFail
 *
 * Specifies the step size change after a nonlinear solver failure
 */

pub fn CVodeSetEtaConvFail(cv_mem: &mut CVodeMem, eta_cf: f64) -> i32 {
    /* set allowed value or use default */
    if eta_cf <= ZERO || eta_cf >= ONE {
        cv_mem.cv_eta_cf = ETA_CF_DEFAULT;
    } else {
        cv_mem.cv_eta_cf = eta_cf;
    }

    CV_SUCCESS
}

/*
 * CVodeSetStopTime
 *
 * Specifies the time beyond which the integration is not to proceed.
 */

pub fn CVodeSetStopTime(cv_mem: &mut CVodeMem, tstop: f64) -> i32 {
    /* If CVode was called at least once, test if tstop is legal
     * (i.e. if it was not already passed).
     * If CVodeSetStopTime is called before the first call to CVode,
     * tstop will be checked in CVode. */
    if cv_mem.cv_nst > 0 && (tstop - cv_mem.cv_tn) * cv_mem.cv_h < ZERO {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSetStopTime", file!(),
            &format!("The value tstop = {} is behind current t = {} in the direction of integration.",
                     tstop, cv_mem.cv_tn));
        return CV_ILL_INPUT;
    }

    cv_mem.cv_tstop = tstop;
    cv_mem.cv_tstopset = SUNTRUE;

    CV_SUCCESS
}

/*
 * CVodeSetInterpolateStopTime
 *
 * Specifies to use interpolation to fill the output solution at
 * the stop time (instead of a copy).
 */

pub fn CVodeSetInterpolateStopTime(cv_mem: &mut CVodeMem, interp: bool) -> i32 {
    cv_mem.cv_tstopinterp = interp;

    CV_SUCCESS
}

/*
 * CVodeClearStopTime
 *
 * Disable the stop time.
 */

pub fn CVodeClearStopTime(cv_mem: &mut CVodeMem) -> i32 {
    cv_mem.cv_tstopset = SUNFALSE;

    CV_SUCCESS
}

/*
 * CVodeSetMaxErrTestFails
 *
 * Specifies the maximum number of error test failures during one
 * step try.
 */

pub fn CVodeSetMaxErrTestFails(cv_mem: &mut CVodeMem, maxnef: i32) -> i32 {
    cv_mem.cv_maxnef = maxnef;

    CV_SUCCESS
}

/*
 * CVodeSetMaxConvFails
 *
 * Specifies the maximum number of nonlinear convergence failures
 * during one step try.
 */

pub fn CVodeSetMaxConvFails(cv_mem: &mut CVodeMem, maxncf: i32) -> i32 {
    cv_mem.cv_maxncf = maxncf;

    CV_SUCCESS
}

/*
 * CVodeSetMaxNonlinIters
 *
 * Specifies the maximum number of nonlinear iterations during
 * one solve.
 */

pub fn CVodeSetMaxNonlinIters(cv_mem: &mut CVodeMem, maxcor: i32) -> i32 {
    /* Are we computing sensitivities with the simultaneous approach? */
    let sensi_sim = cv_mem.cv_sensi && (cv_mem.cv_ism == CV_SIMULTANEOUS);

    if sensi_sim {
        /* check that the NLS is non-NULL */
        if cv_mem.NLSsim.is_none() {
            cvProcessError(None, CV_MEM_FAIL, line!(), "CVodeSetMaxNonlinIters", file!(),
                           MSGCV_MEM_FAIL);
            return CV_MEM_FAIL;
        }

        cv_mem.NLSsim.as_mut().unwrap().set_max_iters(maxcor)
    } else {
        /* check that the NLS is non-NULL */
        if cv_mem.NLS.is_none() {
            cvProcessError(None, CV_MEM_FAIL, line!(), "CVodeSetMaxNonlinIters", file!(),
                           MSGCV_MEM_FAIL);
            return CV_MEM_FAIL;
        }

        cv_mem.NLS.as_mut().unwrap().set_max_iters(maxcor)
    }
}

/*
 * CVodeSetNonlinConvCoef
 *
 * Specifies the coefficient in the nonlinear solver convergence
 * test
 */

pub fn CVodeSetNonlinConvCoef(cv_mem: &mut CVodeMem, nlscoef: f64) -> i32 {
    cv_mem.cv_nlscoef = nlscoef;

    CV_SUCCESS
}

/*
 * CVodeSetLSetupFrequency
 *
 * Specifies the frequency for calling the linear solver setup function to
 * recompute the Jacobian matrix and/or preconditioner
 */

pub fn CVodeSetLSetupFrequency(cv_mem: &mut CVodeMem, msbp: i64) -> i32 {
    /* check for a valid input */
    if msbp < 0 {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSetLSetupFrequency", file!(),
                       "A negative setup frequency was provided");
        return CV_ILL_INPUT;
    }

    /* use default or user provided value */
    cv_mem.cv_msbp = if msbp == 0 { MSBP_DEFAULT } else { msbp };

    CV_SUCCESS
}

/*
 * CVodeSetRootDirection
 *
 * Specifies the direction of zero-crossings to be monitored.
 * The default is to monitor both crossings.
 */

pub fn CVodeSetRootDirection(cv_mem: &mut CVodeMem, rootdir: &[i32]) -> i32 {
    let nrt = cv_mem.cv_nrtfn;
    if nrt == 0 {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSetRootDirection", file!(),
                       MSGCV_NO_ROOT);
        return CV_ILL_INPUT;
    }

    for i in 0..nrt as usize {
        cv_mem.cv_rootdir[i] = rootdir[i];
    }

    CV_SUCCESS
}

/*
 * CVodeSetNoInactiveRootWarn
 *
 * Disables issuing a warning if some root function appears
 * to be identically zero at the beginning of the integration
 */

pub fn CVodeSetNoInactiveRootWarn(cv_mem: &mut CVodeMem) -> i32 {
    cv_mem.cv_mxgnull = 0;

    CV_SUCCESS
}

/*
 * CVodeSetConstraints
 *
 * Setup for constraint handling feature
 */

pub fn CVodeSetConstraints(cv_mem: &mut CVodeMem, constraints: Option<&NVector>) -> i32 {
    /* Disable constraints */
    let constraints = match constraints {
        None => {
            if !cv_mem.cv_constraints.is_empty() {
                cv_mem.cv_constraints = NVector::default();
                cv_mem.cv_lrw -= cv_mem.cv_lrw1;
                cv_mem.cv_liw -= cv_mem.cv_liw1;
            }
            return CV_SUCCESS;
        }
        Some(c) => c,
    };

    /* (The C code tests here that the required vector ops are defined;
       the serial NVector implements them all.) */

    /* Check the constraints vector */
    let temptest = N_VMaxNorm(constraints);
    if (temptest > TWOPT5) || (temptest < HALF) {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSetConstraints", file!(),
                       MSGCV_BAD_CONSTR);
        return CV_ILL_INPUT;
    }

    /* Enable constraints */
    if cv_mem.cv_constraints.is_empty() {
        cv_mem.cv_constraints = NVector::new(constraints.len());
        cv_mem.cv_lrw += cv_mem.cv_lrw1;
        cv_mem.cv_liw += cv_mem.cv_liw1;
    }

    /* Load the constraints vector */
    N_VScale(ONE, constraints, &mut cv_mem.cv_constraints);

    CV_SUCCESS
}

/*
 * CVodeSetMaxNumConstraintFails
 *
 * Set the maximum number of constraint failure allowed in a step
 */

pub fn CVodeSetMaxNumConstraintFails(cv_mem: &mut CVodeMem, max_fails: i32) -> i32 {
    if max_fails <= 0 {
        cv_mem.max_constraint_fails = MAX_CONSTRAINT_FAILS;
    } else {
        cv_mem.max_constraint_fails = max_fails;
    }

    CV_SUCCESS
}

/*
 * CVodeGetNumConstraintFails
 *
 * Get the number of failed steps due to constraint violation
 */

pub fn CVodeGetNumConstraintFails(cv_mem: &mut CVodeMem, num_fails_out: &mut i64) -> i32 {
    *num_fails_out = cv_mem.constraint_fails;

    CV_SUCCESS
}

/*
 * CVodeGetNumConstraintCorrections
 *
 * Get the number of constraint corrections
 */

pub fn CVodeGetNumConstraintCorrections(
    cv_mem: &mut CVodeMem,
    num_corrections_out: &mut i64,
) -> i32 {
    *num_corrections_out = cv_mem.constraint_corrections;

    CV_SUCCESS
}

/*
 * =================================================================
 * Quadrature optional input functions
 * =================================================================
 */

pub fn CVodeSetQuadErrCon(cv_mem: &mut CVodeMem, errconQ: bool) -> i32 {
    cv_mem.cv_errconQ = errconQ;

    CV_SUCCESS
}

/*
 * =================================================================
 * FSA optional input functions
 * =================================================================
 */

pub fn CVodeSetSensDQMethod(cv_mem: &mut CVodeMem, DQtype: i32, DQrhomax: f64) -> i32 {
    if (DQtype != CV_CENTERED) && (DQtype != CV_FORWARD) {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSetSensDQMethod", file!(),
                       MSGCV_BAD_DQTYPE);
        return CV_ILL_INPUT;
    }

    if DQrhomax < ZERO {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSetSensDQMethod", file!(),
                       MSGCV_BAD_DQRHO);
        return CV_ILL_INPUT;
    }

    cv_mem.cv_DQtype = DQtype;
    cv_mem.cv_DQrhomax = DQrhomax;

    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn CVodeSetSensErrCon(cv_mem: &mut CVodeMem, errconS: bool) -> i32 {
    cv_mem.cv_errconS = errconS;

    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn CVodeSetSensMaxNonlinIters(cv_mem: &mut CVodeMem, maxcorS: i32) -> i32 {
    /* Are we computing sensitivities with a staggered approach? */
    let sensi_stg = cv_mem.cv_sensi && (cv_mem.cv_ism == CV_STAGGERED);

    if sensi_stg {
        /* check that the NLS is non-NULL */
        if cv_mem.NLSstg.is_none() {
            cvProcessError(None, CV_MEM_FAIL, line!(), "CVodeSetSensMaxNonlinIters", file!(),
                           MSGCV_MEM_FAIL);
            return CV_MEM_FAIL;
        }

        cv_mem.NLSstg.as_mut().unwrap().set_max_iters(maxcorS)
    } else {
        /* check that the NLS is non-NULL */
        if cv_mem.NLSstg1.is_none() {
            cvProcessError(None, CV_MEM_FAIL, line!(), "CVodeSetSensMaxNonlinIters", file!(),
                           MSGCV_MEM_FAIL);
            return CV_MEM_FAIL;
        }

        cv_mem.NLSstg1.as_mut().unwrap().set_max_iters(maxcorS)
    }
}

/*-----------------------------------------------------------------*/

pub fn CVodeSetSensParams(
    cv_mem: &mut CVodeMem,
    p: Option<&[f64]>,
    pbar: Option<&[f64]>,
    plist: Option<&[i32]>,
) -> i32 {
    /* Was sensitivity initialized? */

    if cv_mem.cv_SensMallocDone == SUNFALSE {
        cvProcessError(Some(cv_mem), CV_NO_SENS, line!(), "CVodeSetSensParams", file!(),
                       MSGCV_NO_SENSI);
        return CV_NO_SENS;
    }

    let Ns = cv_mem.cv_Ns;

    /* Parameters (C stores the user's pointer, cv_mem->cv_p = p;
       cv_p is an owned Vec here, so the slice is copied -- None
       leaves the empty Vec, the C NULL state) */

    cv_mem.cv_p = match p {
        Some(p) => p.to_vec(),
        None => Vec::new(),
    };

    /* pbar */

    if let Some(pbar) = pbar {
        for is in 0..Ns as usize {
            if pbar[is] == ZERO {
                cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSetSensParams",
                               file!(), MSGCV_BAD_PBAR);
                return CV_ILL_INPUT;
            }
            cv_mem.cv_pbar[is] = SUNRabs(pbar[is]);
        }
    } else {
        for is in 0..Ns as usize {
            cv_mem.cv_pbar[is] = ONE;
        }
    }

    /* plist */

    if let Some(plist) = plist {
        for is in 0..Ns as usize {
            if plist[is] < 0 {
                cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSetSensParams",
                               file!(), MSGCV_BAD_PLIST);
                return CV_ILL_INPUT;
            }
            cv_mem.cv_plist[is] = plist[is];
        }
    } else {
        for is in 0..Ns as usize {
            cv_mem.cv_plist[is] = is as i32;
        }
    }

    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn CVodeSetQuadSensErrCon(cv_mem: &mut CVodeMem, errconQS: bool) -> i32 {
    /* Was sensitivity initialized? */

    if cv_mem.cv_SensMallocDone == SUNFALSE {
        cvProcessError(Some(cv_mem), CV_NO_SENS, line!(), "CVodeSetQuadSensErrCon", file!(),
                       MSGCV_NO_SENSI);
        return CV_NO_SENS;
    }

    /* Check if quadrature sensitivity was initialized? */

    if cv_mem.cv_QuadSensMallocDone == SUNFALSE {
        cvProcessError(Some(cv_mem), CV_NO_QUADSENS, line!(), "CVodeSetQuadSensErrCon", file!(),
                       MSGCV_NO_QUADSENSI);
        /* (sic -- the C source returns CV_NO_QUAD here, not CV_NO_QUADSENS) */
        return CV_NO_QUAD;
    }

    cv_mem.cv_errconQS = errconQS;

    CV_SUCCESS
}

/*
 * =================================================================
 * CVODES optional output functions
 * =================================================================
 */

/*
 * CVodeGetNumSteps
 *
 * Returns the current number of integration steps
 */

pub fn CVodeGetNumSteps(cv_mem: &mut CVodeMem, nsteps: &mut i64) -> i32 {
    *nsteps = cv_mem.cv_nst;

    CV_SUCCESS
}

/*
 * CVodeGetNumRhsEvals
 *
 * Returns the current number of calls to f
 */

pub fn CVodeGetNumRhsEvals(cv_mem: &mut CVodeMem, nfevals: &mut i64) -> i32 {
    *nfevals = cv_mem.cv_nfe;

    CV_SUCCESS
}

/*
 * CVodeGetNumLinSolvSetups
 *
 * Returns the current number of calls to the linear solver setup routine
 */

pub fn CVodeGetNumLinSolvSetups(cv_mem: &mut CVodeMem, nlinsetups: &mut i64) -> i32 {
    *nlinsetups = cv_mem.cv_nsetups;

    CV_SUCCESS
}

/*
 * CVodeGetNumErrTestFails
 *
 * Returns the current number of error test failures
 */

pub fn CVodeGetNumErrTestFails(cv_mem: &mut CVodeMem, netfails: &mut i64) -> i32 {
    *netfails = cv_mem.cv_netf;

    CV_SUCCESS
}

/*
 * CVodeGetLastOrder
 *
 * Returns the order on the last successful step
 */

pub fn CVodeGetLastOrder(cv_mem: &mut CVodeMem, qlast: &mut i32) -> i32 {
    *qlast = cv_mem.cv_qu;

    CV_SUCCESS
}

/*
 * CVodeGetCurrentOrder
 *
 * Returns the order to be attempted on the next step
 */

pub fn CVodeGetCurrentOrder(cv_mem: &mut CVodeMem, qcur: &mut i32) -> i32 {
    *qcur = cv_mem.cv_next_q;

    CV_SUCCESS
}

/*
 * CVodeGetCurrentGamma
 *
 * Returns the value of gamma for the current step.
 */

pub fn CVodeGetCurrentGamma(cv_mem: &mut CVodeMem, gamma: &mut f64) -> i32 {
    *gamma = cv_mem.cv_gamma;

    CV_SUCCESS
}

/*
 * CVodeGetNumStabLimOrderReds
 *
 * Returns the number of order reductions triggered by the stability
 * limit detection algorithm
 */

pub fn CVodeGetNumStabLimOrderReds(cv_mem: &mut CVodeMem, nslred: &mut i64) -> i32 {
    if cv_mem.cv_sldeton == SUNFALSE {
        *nslred = 0;
    } else {
        *nslred = cv_mem.cv_nor;
    }

    CV_SUCCESS
}

/*
 * CVodeGetActualInitStep
 *
 * Returns the step size used on the first step
 */

pub fn CVodeGetActualInitStep(cv_mem: &mut CVodeMem, hinused: &mut f64) -> i32 {
    *hinused = cv_mem.cv_h0u;

    CV_SUCCESS
}

/*
 * CVodeGetLastStep
 *
 * Returns the step size used on the last successful step
 */

pub fn CVodeGetLastStep(cv_mem: &mut CVodeMem, hlast: &mut f64) -> i32 {
    *hlast = cv_mem.cv_hu;

    CV_SUCCESS
}

/*
 * CVodeGetCurrentStep
 *
 * Returns the step size to be attempted on the next step
 */

pub fn CVodeGetCurrentStep(cv_mem: &mut CVodeMem, hcur: &mut f64) -> i32 {
    *hcur = cv_mem.cv_next_h;

    CV_SUCCESS
}

/*
 * CVodeGetCurrentState
 *
 * Returns the current state vector (in C: *y = cv_mem->cv_y).
 */

pub fn CVodeGetCurrentState(cv_mem: &CVodeMem) -> &NVector {
    &cv_mem.cv_y
}

/*
 * CVodeGetCurrentStateSens
 *
 * Returns the current sensitivity state vector array
 * (in C: *yS = cv_mem->cv_yS).
 */

pub fn CVodeGetCurrentStateSens(cv_mem: &CVodeMem) -> &[NVector] {
    &cv_mem.cv_yS
}

/*
 * CVodeGetCurrentSensSolveIndex
 *
 * Returns the current index of the sensitivity solve when using
 * the staggered1 nonlinear solver.
 */

pub fn CVodeGetCurrentSensSolveIndex(cv_mem: &mut CVodeMem, index: &mut i32) -> i32 {
    *index = cv_mem.sens_solve_idx;

    CV_SUCCESS
}

/*
 * CVodeGetCurrentTime
 *
 * Returns the current value of the independent variable
 */

pub fn CVodeGetCurrentTime(cv_mem: &mut CVodeMem, tcur: &mut f64) -> i32 {
    *tcur = cv_mem.cv_tn;

    CV_SUCCESS
}

/*
 * CVodeGetTolScaleFactor
 *
 * Returns a suggested factor for scaling tolerances
 */

pub fn CVodeGetTolScaleFactor(cv_mem: &mut CVodeMem, tolsfact: &mut f64) -> i32 {
    *tolsfact = cv_mem.cv_tolsf;

    CV_SUCCESS
}

/*
 * CVodeGetErrWeights
 *
 * This routine returns the current weight vector.
 */

pub fn CVodeGetErrWeights(cv_mem: &mut CVodeMem, eweight: &mut NVector) -> i32 {
    N_VScale(ONE, &cv_mem.cv_ewt, eweight);

    CV_SUCCESS
}

/*
 * CVodeGetEstLocalErrors
 *
 * Returns an estimate of the local error
 */

pub fn CVodeGetEstLocalErrors(cv_mem: &mut CVodeMem, ele: &mut NVector) -> i32 {
    N_VScale(ONE, &cv_mem.cv_acor, ele);

    CV_SUCCESS
}

/*
 * CVodeGetWorkSpace
 *
 * Returns integrator work space requirements
 */

pub fn CVodeGetWorkSpace(cv_mem: &mut CVodeMem, lenrw: &mut i64, leniw: &mut i64) -> i32 {
    *leniw = cv_mem.cv_liw;
    *lenrw = cv_mem.cv_lrw;

    CV_SUCCESS
}

/*
 * CVodeGetIntegratorStats
 *
 * Returns integrator statistics
 */

pub fn CVodeGetIntegratorStats(
    cv_mem: &mut CVodeMem,
    nsteps: &mut i64,
    nfevals: &mut i64,
    nlinsetups: &mut i64,
    netfails: &mut i64,
    qlast: &mut i32,
    qcur: &mut i32,
    hinused: &mut f64,
    hlast: &mut f64,
    hcur: &mut f64,
    tcur: &mut f64,
) -> i32 {
    *nsteps = cv_mem.cv_nst;
    *nfevals = cv_mem.cv_nfe;
    *nlinsetups = cv_mem.cv_nsetups;
    *netfails = cv_mem.cv_netf;
    *qlast = cv_mem.cv_qu;
    *qcur = cv_mem.cv_next_q;
    *hinused = cv_mem.cv_h0u;
    *hlast = cv_mem.cv_hu;
    *hcur = cv_mem.cv_next_h;
    *tcur = cv_mem.cv_tn;

    CV_SUCCESS
}

/*
 * CVodeGetNumGEvals
 *
 * Returns the current number of calls to g (for rootfinding)
 */

pub fn CVodeGetNumGEvals(cv_mem: &mut CVodeMem, ngevals: &mut i64) -> i32 {
    *ngevals = cv_mem.cv_nge;

    CV_SUCCESS
}

/*
 * CVodeGetRootInfo
 *
 * Returns pointer to array rootsfound showing roots found
 */

pub fn CVodeGetRootInfo(cv_mem: &mut CVodeMem, rootsfound: &mut [i32]) -> i32 {
    let nrt = cv_mem.cv_nrtfn;

    for i in 0..nrt as usize {
        rootsfound[i] = cv_mem.cv_iroots[i];
    }

    CV_SUCCESS
}

/*
 * CVodeGetNumNonlinSolvIters
 *
 * Returns the current number of iterations in the nonlinear solver
 */

pub fn CVodeGetNumNonlinSolvIters(cv_mem: &mut CVodeMem, nniters: &mut i64) -> i32 {
    *nniters = cv_mem.cv_nni;

    CV_SUCCESS
}

/*
 * CVodeGetNumNonlinSolvConvFails
 *
 * Returns the current number of convergence failures in the
 * nonlinear solver
 */

pub fn CVodeGetNumNonlinSolvConvFails(cv_mem: &mut CVodeMem, nnfails: &mut i64) -> i32 {
    *nnfails = cv_mem.cv_nnf;

    CV_SUCCESS
}

/*
 * CVodeGetNonlinSolvStats
 *
 * Returns nonlinear solver statistics
 */

pub fn CVodeGetNonlinSolvStats(cv_mem: &mut CVodeMem, nniters: &mut i64, nnfails: &mut i64) -> i32 {
    *nniters = cv_mem.cv_nni;
    *nnfails = cv_mem.cv_nnf;

    CV_SUCCESS
}

/*
 * CVodeGetNumStepSolveFails
 *
 * Returns the current number of failed steps due to a nonlinear solver
 * convergence failure
 */

pub fn CVodeGetNumStepSolveFails(cv_mem: &mut CVodeMem, nncfails: &mut i64) -> i32 {
    *nncfails = cv_mem.cv_ncfn;

    CV_SUCCESS
}

/*
 * =================================================================
 * Quadrature optional output functions
 * =================================================================
 */

/*-----------------------------------------------------------------*/

pub fn CVodeGetQuadNumRhsEvals(cv_mem: &mut CVodeMem, nfQevals: &mut i64) -> i32 {
    if cv_mem.cv_quadr == SUNFALSE {
        cvProcessError(Some(cv_mem), CV_NO_QUAD, line!(), "CVodeGetQuadNumRhsEvals", file!(),
                       MSGCV_NO_QUAD);
        return CV_NO_QUAD;
    }

    *nfQevals = cv_mem.cv_nfQe;

    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn CVodeGetQuadNumErrTestFails(cv_mem: &mut CVodeMem, nQetfails: &mut i64) -> i32 {
    if cv_mem.cv_quadr == SUNFALSE {
        cvProcessError(Some(cv_mem), CV_NO_QUAD, line!(), "CVodeGetQuadNumErrTestFails", file!(),
                       MSGCV_NO_QUAD);
        return CV_NO_QUAD;
    }

    *nQetfails = cv_mem.cv_netfQ;

    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn CVodeGetQuadErrWeights(cv_mem: &mut CVodeMem, eQweight: &mut NVector) -> i32 {
    if cv_mem.cv_quadr == SUNFALSE {
        cvProcessError(Some(cv_mem), CV_NO_QUAD, line!(), "CVodeGetQuadErrWeights", file!(),
                       MSGCV_NO_QUAD);
        return CV_NO_QUAD;
    }

    if cv_mem.cv_errconQ {
        N_VScale(ONE, &cv_mem.cv_ewtQ, eQweight);
    }

    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn CVodeGetQuadStats(cv_mem: &mut CVodeMem, nfQevals: &mut i64, nQetfails: &mut i64) -> i32 {
    if cv_mem.cv_quadr == SUNFALSE {
        cvProcessError(Some(cv_mem), CV_NO_QUAD, line!(), "CVodeGetQuadStats", file!(),
                       MSGCV_NO_QUAD);
        return CV_NO_QUAD;
    }

    *nfQevals = cv_mem.cv_nfQe;
    *nQetfails = cv_mem.cv_netfQ;

    CV_SUCCESS
}

/*
 * =================================================================
 * Quadrature FSA optional output functions
 * =================================================================
 */

/*-----------------------------------------------------------------*/

pub fn CVodeGetQuadSensNumRhsEvals(cv_mem: &mut CVodeMem, nfQSevals: &mut i64) -> i32 {
    if cv_mem.cv_quadr_sensi == SUNFALSE {
        cvProcessError(Some(cv_mem), CV_NO_QUADSENS, line!(), "CVodeGetQuadSensNumRhsEvals",
                       file!(), MSGCV_NO_QUADSENSI);
        return CV_NO_QUADSENS;
    }

    *nfQSevals = cv_mem.cv_nfQSe;

    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn CVodeGetQuadSensNumErrTestFails(cv_mem: &mut CVodeMem, nQSetfails: &mut i64) -> i32 {
    if cv_mem.cv_quadr_sensi == SUNFALSE {
        cvProcessError(Some(cv_mem), CV_NO_QUADSENS, line!(), "CVodeGetQuadSensNumErrTestFails",
                       file!(), MSGCV_NO_QUADSENSI);
        return CV_NO_QUADSENS;
    }

    *nQSetfails = cv_mem.cv_netfQS;

    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn CVodeGetQuadSensErrWeights(cv_mem: &mut CVodeMem, eQSweight: &mut [NVector]) -> i32 {
    if cv_mem.cv_quadr_sensi == SUNFALSE {
        cvProcessError(Some(cv_mem), CV_NO_QUADSENS, line!(), "CVodeGetQuadSensErrWeights",
                       file!(), MSGCV_NO_QUADSENSI);
        return CV_NO_QUADSENS;
    }
    let Ns = cv_mem.cv_Ns;

    if cv_mem.cv_errconQS {
        for is in 0..Ns as usize {
            N_VScale(ONE, &cv_mem.cv_ewtQS[is], &mut eQSweight[is]);
        }
    }

    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn CVodeGetQuadSensStats(
    cv_mem: &mut CVodeMem,
    nfQSevals: &mut i64,
    nQSetfails: &mut i64,
) -> i32 {
    if cv_mem.cv_quadr_sensi == SUNFALSE {
        cvProcessError(Some(cv_mem), CV_NO_QUADSENS, line!(), "CVodeGetQuadSensStats", file!(),
                       MSGCV_NO_QUADSENSI);
        return CV_NO_QUADSENS;
    }

    *nfQSevals = cv_mem.cv_nfQSe;
    *nQSetfails = cv_mem.cv_netfQS;

    CV_SUCCESS
}

/*
 * =================================================================
 * FSA optional output functions
 * =================================================================
 */

/*-----------------------------------------------------------------*/

pub fn CVodeGetSensNumRhsEvals(cv_mem: &mut CVodeMem, nfSevals: &mut i64) -> i32 {
    if cv_mem.cv_sensi == SUNFALSE {
        cvProcessError(Some(cv_mem), CV_NO_SENS, line!(), "CVodeGetSensNumRhsEvals", file!(),
                       MSGCV_NO_SENSI);
        return CV_NO_SENS;
    }

    *nfSevals = cv_mem.cv_nfSe;

    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn CVodeGetNumRhsEvalsSens(cv_mem: &mut CVodeMem, nfevalsS: &mut i64) -> i32 {
    if cv_mem.cv_sensi == SUNFALSE {
        cvProcessError(Some(cv_mem), CV_NO_SENS, line!(), "CVodeGetNumRhsEvalsSens", file!(),
                       MSGCV_NO_SENSI);
        return CV_NO_SENS;
    }

    *nfevalsS = cv_mem.cv_nfeS;

    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn CVodeGetSensNumErrTestFails(cv_mem: &mut CVodeMem, nSetfails: &mut i64) -> i32 {
    if cv_mem.cv_sensi == SUNFALSE {
        cvProcessError(Some(cv_mem), CV_NO_SENS, line!(), "CVodeGetSensNumErrTestFails", file!(),
                       MSGCV_NO_SENSI);
        return CV_NO_SENS;
    }

    *nSetfails = cv_mem.cv_netfS;

    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn CVodeGetSensNumLinSolvSetups(cv_mem: &mut CVodeMem, nlinsetupsS: &mut i64) -> i32 {
    if cv_mem.cv_sensi == SUNFALSE {
        cvProcessError(Some(cv_mem), CV_NO_SENS, line!(), "CVodeGetSensNumLinSolvSetups",
                       file!(), MSGCV_NO_SENSI);
        return CV_NO_SENS;
    }

    *nlinsetupsS = cv_mem.cv_nsetupsS;

    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn CVodeGetSensErrWeights(cv_mem: &mut CVodeMem, eSweight: &mut [NVector]) -> i32 {
    if cv_mem.cv_sensi == SUNFALSE {
        cvProcessError(Some(cv_mem), CV_NO_SENS, line!(), "CVodeGetSensErrWeights", file!(),
                       MSGCV_NO_SENSI);
        return CV_NO_SENS;
    }

    let Ns = cv_mem.cv_Ns;

    for is in 0..Ns as usize {
        N_VScale(ONE, &cv_mem.cv_ewtS[is], &mut eSweight[is]);
    }

    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn CVodeGetSensStats(
    cv_mem: &mut CVodeMem,
    nfSevals: &mut i64,
    nfevalsS: &mut i64,
    nSetfails: &mut i64,
    nlinsetupsS: &mut i64,
) -> i32 {
    if cv_mem.cv_sensi == SUNFALSE {
        cvProcessError(Some(cv_mem), CV_NO_SENS, line!(), "CVodeGetSensStats", file!(),
                       MSGCV_NO_SENSI);
        return CV_NO_SENS;
    }

    *nfSevals = cv_mem.cv_nfSe;
    *nfevalsS = cv_mem.cv_nfeS;
    *nSetfails = cv_mem.cv_netfS;
    *nlinsetupsS = cv_mem.cv_nsetupsS;

    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn CVodeGetSensNumNonlinSolvIters(cv_mem: &mut CVodeMem, nSniters: &mut i64) -> i32 {
    if cv_mem.cv_sensi == SUNFALSE {
        cvProcessError(Some(cv_mem), CV_NO_SENS, line!(), "CVodeGetSensNumNonlinSolvIters",
                       file!(), MSGCV_NO_SENSI);
        return CV_NO_SENS;
    }

    *nSniters = cv_mem.cv_nniS;

    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn CVodeGetSensNumNonlinSolvConvFails(cv_mem: &mut CVodeMem, nSnfails: &mut i64) -> i32 {
    if cv_mem.cv_sensi == SUNFALSE {
        cvProcessError(Some(cv_mem), CV_NO_SENS, line!(), "CVodeGetSensNumNonlinSolvConvFails",
                       file!(), MSGCV_NO_SENSI);
        return CV_NO_SENS;
    }

    *nSnfails = cv_mem.cv_nnfS;

    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn CVodeGetSensNonlinSolvStats(
    cv_mem: &mut CVodeMem,
    nSniters: &mut i64,
    nSnfails: &mut i64,
) -> i32 {
    if cv_mem.cv_sensi == SUNFALSE {
        cvProcessError(Some(cv_mem), CV_NO_SENS, line!(), "CVodeGetSensNonlinSolvStats",
                       file!(), MSGCV_NO_SENSI);
        return CV_NO_SENS;
    }

    *nSniters = cv_mem.cv_nniS;
    *nSnfails = cv_mem.cv_nnfS;

    CV_SUCCESS
}

pub fn CVodeGetNumStepSensSolveFails(cv_mem: &mut CVodeMem, nSncfails: &mut i64) -> i32 {
    if cv_mem.cv_sensi == SUNFALSE {
        cvProcessError(Some(cv_mem), CV_NO_SENS, line!(), "CVodeGetNumStepSensSolveFails",
                       file!(), MSGCV_NO_SENSI);
        return CV_NO_SENS;
    }

    *nSncfails = cv_mem.cv_ncfnS;

    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn CVodeGetStgrSensNumNonlinSolvIters(cv_mem: &mut CVodeMem, nSTGR1niters: &mut [i64]) -> i32 {
    let Ns = cv_mem.cv_Ns;

    if cv_mem.cv_sensi == SUNFALSE {
        cvProcessError(Some(cv_mem), CV_NO_SENS, line!(), "CVodeGetStgrSensNumNonlinSolvIters",
                       file!(), MSGCV_NO_SENSI);
        return CV_NO_SENS;
    }

    if cv_mem.cv_ism == CV_STAGGERED1 {
        for is in 0..Ns as usize {
            nSTGR1niters[is] = cv_mem.cv_nniS1[is];
        }
    }

    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn CVodeGetStgrSensNumNonlinSolvConvFails(
    cv_mem: &mut CVodeMem,
    nSTGR1nfails: &mut [i64],
) -> i32 {
    let Ns = cv_mem.cv_Ns;

    if cv_mem.cv_sensi == SUNFALSE {
        cvProcessError(Some(cv_mem), CV_NO_SENS, line!(),
                       "CVodeGetStgrSensNumNonlinSolvConvFails", file!(), MSGCV_NO_SENSI);
        return CV_NO_SENS;
    }

    if cv_mem.cv_ism == CV_STAGGERED1 {
        for is in 0..Ns as usize {
            nSTGR1nfails[is] = cv_mem.cv_nnfS1[is];
        }
    }

    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn CVodeGetStgrSensNonlinSolvStats(
    cv_mem: &mut CVodeMem,
    nSTGR1niters: &mut [i64],
    nSTGR1nfails: &mut [i64],
) -> i32 {
    let Ns = cv_mem.cv_Ns;

    if cv_mem.cv_sensi == SUNFALSE {
        cvProcessError(Some(cv_mem), CV_NO_SENS, line!(), "CVodeGetStgrSensNonlinSolvStats",
                       file!(), MSGCV_NO_SENSI);
        return CV_NO_SENS;
    }

    if cv_mem.cv_ism == CV_STAGGERED1 {
        for is in 0..Ns as usize {
            nSTGR1niters[is] = cv_mem.cv_nniS1[is];
        }
        for is in 0..Ns as usize {
            nSTGR1nfails[is] = cv_mem.cv_nnfS1[is];
        }
    }

    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn CVodeGetNumStepStgrSensSolveFails(
    cv_mem: &mut CVodeMem,
    nSTGR1ncfails: &mut [i64],
) -> i32 {
    let Ns = cv_mem.cv_Ns;

    if cv_mem.cv_sensi == SUNFALSE {
        cvProcessError(Some(cv_mem), CV_NO_SENS, line!(), "CVodeGetNumStepStgrSensSolveFails",
                       file!(), MSGCV_NO_SENSI);
        return CV_NO_SENS;
    }

    if cv_mem.cv_ism == CV_STAGGERED1 {
        for is in 0..Ns as usize {
            nSTGR1ncfails[is] = cv_mem.cv_ncfnS1[is];
        }
    }

    CV_SUCCESS
}

/* -----------------------------------------------------------------
 * Counterparts of sunfprintf_real / sunfprintf_long /
 * sunfprintf_long_array (src/sundials/sundials_utils.h).
 * SUN_FORMAT_G is "%.15g" and SUN_FORMAT_E is "% .15e" for double
 * precision.
 * -----------------------------------------------------------------*/

const SUN_TABLE_WIDTH: usize = 29;

fn sunfprintf_real(
    outfile: &mut dyn std::io::Write,
    fmt: SUNOutputFormat,
    start: bool,
    name: &str,
    value: f64,
) {
    if fmt == SUN_OUTPUTFORMAT_TABLE {
        let _ = writeln!(outfile, "{:<width$} = {}", name, fmt_g(value, 0, 15),
                         width = SUN_TABLE_WIDTH);
    } else {
        if !start {
            let _ = write!(outfile, ",");
        }
        /* C "% .15e": a space is printed in place of a plus sign */
        let e = fmt_e(value, 0, 15);
        let e = if e.starts_with('-') { e } else { format!(" {}", e) };
        let _ = write!(outfile, "{},{}", name, e);
    }
}

fn sunfprintf_long(
    outfile: &mut dyn std::io::Write,
    fmt: SUNOutputFormat,
    start: bool,
    name: &str,
    value: i64,
) {
    if fmt == SUN_OUTPUTFORMAT_TABLE {
        let _ = writeln!(outfile, "{:<width$} = {}", name, value, width = SUN_TABLE_WIDTH);
    } else {
        if !start {
            let _ = write!(outfile, ",");
        }
        let _ = write!(outfile, "{},{}", name, value);
    }
}

/* (C signature: (FILE*, fmt, start, name, long* value, size_t count);
   here the count is the slice length.) */
fn sunfprintf_long_array(
    outfile: &mut dyn std::io::Write,
    fmt: SUNOutputFormat,
    start: bool,
    name: &str,
    value: &[i64],
) {
    if value.is_empty() {
        return;
    }

    if fmt == SUN_OUTPUTFORMAT_TABLE {
        let _ = write!(outfile, "{:<width$} = {}", name, value[0], width = SUN_TABLE_WIDTH);
        for v in &value[1..] {
            let _ = write!(outfile, ", {}", v);
        }
        let _ = writeln!(outfile);
    } else {
        if !start {
            let _ = write!(outfile, ",");
        }
        /* (sic -- as in C, successive "name i,value" groups of one
           array are printed with no separating comma between them) */
        for (i, v) in value.iter().enumerate() {
            let _ = write!(outfile, "{} {},{}", name, i, v);
        }
    }
}

/*
 * CVodePrintAllStats
 *
 * Print all integrator statistics
 */

pub fn CVodePrintAllStats(
    cv_mem: &mut CVodeMem,
    outfile: &mut dyn std::io::Write,
    fmt: SUNOutputFormat,
) -> i32 {
    if fmt != SUN_OUTPUTFORMAT_TABLE && fmt != SUN_OUTPUTFORMAT_CSV {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodePrintAllStats", file!(),
                       "Invalid formatting option.");
        return CV_ILL_INPUT;
    }

    /* step and method stats */
    sunfprintf_real(outfile, fmt, SUNTRUE, "Current time", cv_mem.cv_tn);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Steps", cv_mem.cv_nst);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Error test fails", cv_mem.cv_netf);
    sunfprintf_long(outfile, fmt, SUNFALSE, "NLS step fails", cv_mem.cv_ncfn);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Constraint fails", cv_mem.constraint_fails);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Constraint corrections",
                    cv_mem.constraint_corrections);
    sunfprintf_real(outfile, fmt, SUNFALSE, "Initial step size", cv_mem.cv_h0u);
    sunfprintf_real(outfile, fmt, SUNFALSE, "Last step size", cv_mem.cv_hu);
    sunfprintf_real(outfile, fmt, SUNFALSE, "Current step size", cv_mem.cv_next_h);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Last method order", cv_mem.cv_qu as i64);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Current method order", cv_mem.cv_next_q as i64);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Stab. lim. order reductions", cv_mem.cv_nor);
    /* function evaluations */
    sunfprintf_long(outfile, fmt, SUNFALSE, "RHS fn evals", cv_mem.cv_nfe);
    /* nonlinear solver stats */
    sunfprintf_long(outfile, fmt, SUNFALSE, "NLS iters", cv_mem.cv_nni);
    sunfprintf_long(outfile, fmt, SUNFALSE, "NLS fails", cv_mem.cv_nnf);
    if cv_mem.cv_nst > 0 {
        sunfprintf_real(outfile, fmt, SUNFALSE, "NLS iters per step",
                        cv_mem.cv_nni as f64 / cv_mem.cv_nst as f64);
    }

    /* linear solver stats */
    sunfprintf_long(outfile, fmt, SUNFALSE, "LS setups", cv_mem.cv_nsetups);
    /* (In C this block casts cv_lmem to CVLsMem whenever a linear solver
       module is attached; the counters below only exist for the CVLS
       interface, so the CVDiag module prints no additional stats here.) */
    if let LsModule::Ls(cvls_mem) = &cv_mem.cv_lmem {
        sunfprintf_long(outfile, fmt, SUNFALSE, "Jac fn evals", cvls_mem.nje);
        sunfprintf_long(outfile, fmt, SUNFALSE, "LS RHS fn evals", cvls_mem.nfeDQ);
        sunfprintf_long(outfile, fmt, SUNFALSE, "Prec setup evals", cvls_mem.npe);
        sunfprintf_long(outfile, fmt, SUNFALSE, "Prec solves", cvls_mem.nps);
        sunfprintf_long(outfile, fmt, SUNFALSE, "LS iters", cvls_mem.nli);
        sunfprintf_long(outfile, fmt, SUNFALSE, "LS fails", cvls_mem.ncfl);
        sunfprintf_long(outfile, fmt, SUNFALSE, "Jac-times setups", cvls_mem.njtsetup);
        sunfprintf_long(outfile, fmt, SUNFALSE, "Jac-times evals", cvls_mem.njtimes);
        if cv_mem.cv_nni > 0 {
            sunfprintf_real(outfile, fmt, SUNFALSE, "LS iters per NLS iter",
                            cvls_mem.nli as f64 / cv_mem.cv_nni as f64);
            sunfprintf_real(outfile, fmt, SUNFALSE, "Jac evals per NLS iter",
                            cvls_mem.nje as f64 / cv_mem.cv_nni as f64);
            sunfprintf_real(outfile, fmt, SUNFALSE, "Prec evals per NLS iter",
                            cvls_mem.npe as f64 / cv_mem.cv_nni as f64);
        }
    }

    /* rootfinding stats */
    sunfprintf_long(outfile, fmt, SUNFALSE, "Root fn evals", cv_mem.cv_nge);

    /* projection stats */
    if let Some(cvproj_mem) = cv_mem.proj_mem.as_deref() {
        sunfprintf_long(outfile, fmt, SUNFALSE, "Projection fn evals", cvproj_mem.nproj);
        sunfprintf_long(outfile, fmt, SUNFALSE, "Projection fails", cvproj_mem.npfails);
    }

    /* quadrature stats */
    if cv_mem.cv_quadr {
        sunfprintf_long(outfile, fmt, SUNFALSE, "Quad fn evals", cv_mem.cv_nfQe);
        sunfprintf_long(outfile, fmt, SUNFALSE, "Quad error test fails", cv_mem.cv_netfQ);
    }

    /* sensitivity stats */
    if cv_mem.cv_sensi {
        sunfprintf_long(outfile, fmt, SUNFALSE, "Sens fn evals", cv_mem.cv_nfSe);
        sunfprintf_long(outfile, fmt, SUNFALSE, "Sens RHS fn evals", cv_mem.cv_nfeS);
        sunfprintf_long(outfile, fmt, SUNFALSE, "Sens error test fails", cv_mem.cv_netfS);
        if cv_mem.cv_ism != CV_SIMULTANEOUS {
            sunfprintf_long(outfile, fmt, SUNFALSE, "Sens NLS iters", cv_mem.cv_nniS);
            sunfprintf_long(outfile, fmt, SUNFALSE, "Sens NLS fails", cv_mem.cv_nnfS);
            sunfprintf_long(outfile, fmt, SUNFALSE, "Sens NLS step fails", cv_mem.cv_ncfnS);
        }
        if cv_mem.cv_ism == CV_STAGGERED1 {
            let ns = cv_mem.cv_Ns as usize;
            sunfprintf_long_array(outfile, fmt, SUNFALSE, "Sens stgr1 NLS iters",
                                  &cv_mem.cv_nniS1[..ns]);
            sunfprintf_long_array(outfile, fmt, SUNFALSE, "Sens stgr1 NLS fails",
                                  &cv_mem.cv_nnfS1[..ns]);
            sunfprintf_long_array(outfile, fmt, SUNFALSE, "Sens stgr1 NLS step fails",
                                  &cv_mem.cv_ncfnS1[..ns]);
        }
        sunfprintf_long(outfile, fmt, SUNFALSE, "Sens LS setups", cv_mem.cv_nsetupsS);
    }

    /* quadrature-sensitivity stats */
    if cv_mem.cv_quadr_sensi {
        sunfprintf_long(outfile, fmt, SUNFALSE, "QuadSens fn evals", cv_mem.cv_nfQSe);
        sunfprintf_long(outfile, fmt, SUNFALSE, "QuadSens error test fails", cv_mem.cv_netfQS);
    }

    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn CVodeGetUserData<'a>(cv_mem: &'a mut CVodeMem) -> &'a mut UserData {
    &mut cv_mem.cv_user_data
}

/*-----------------------------------------------------------------*/

pub fn CVodeGetReturnFlagName(flag: i64) -> String {
    let flag_i32 = i32::try_from(flag).unwrap_or(i32::MIN); /* out-of-range -> "NONE" */
    let name = match flag_i32 {
        CV_SUCCESS => "CV_SUCCESS",
        CV_TSTOP_RETURN => "CV_TSTOP_RETURN",
        CV_ROOT_RETURN => "CV_ROOT_RETURN",
        CV_TOO_MUCH_WORK => "CV_TOO_MUCH_WORK",
        CV_TOO_MUCH_ACC => "CV_TOO_MUCH_ACC",
        CV_ERR_FAILURE => "CV_ERR_FAILURE",
        CV_CONV_FAILURE => "CV_CONV_FAILURE",
        CV_LINIT_FAIL => "CV_LINIT_FAIL",
        CV_LSETUP_FAIL => "CV_LSETUP_FAIL",
        CV_LSOLVE_FAIL => "CV_LSOLVE_FAIL",
        CV_RHSFUNC_FAIL => "CV_RHSFUNC_FAIL",
        CV_FIRST_RHSFUNC_ERR => "CV_FIRST_RHSFUNC_ERR",
        CV_REPTD_RHSFUNC_ERR => "CV_REPTD_RHSFUNC_ERR",
        CV_UNREC_RHSFUNC_ERR => "CV_UNREC_RHSFUNC_ERR",
        CV_RTFUNC_FAIL => "CV_RTFUNC_FAIL",
        CV_MEM_FAIL => "CV_MEM_FAIL",
        CV_MEM_NULL => "CV_MEM_NULL",
        CV_ILL_INPUT => "CV_ILL_INPUT",
        CV_NO_MALLOC => "CV_NO_MALLOC",
        CV_BAD_K => "CV_BAD_K",
        CV_BAD_T => "CV_BAD_T",
        CV_BAD_DKY => "CV_BAD_DKY",
        CV_NO_QUAD => "CV_NO_QUAD",
        CV_QRHSFUNC_FAIL => "CV_QRHSFUNC_FAIL",
        CV_FIRST_QRHSFUNC_ERR => "CV_FIRST_QRHSFUNC_ERR",
        CV_REPTD_QRHSFUNC_ERR => "CV_REPTD_QRHSFUNC_ERR",
        CV_UNREC_QRHSFUNC_ERR => "CV_UNREC_QRHSFUNC_ERR",
        CV_BAD_IS => "CV_BAD_IS",
        CV_NO_SENS => "CV_NO_SENS",
        CV_SRHSFUNC_FAIL => "CV_SRHSFUNC_FAIL",
        CV_FIRST_SRHSFUNC_ERR => "CV_FIRST_SRHSFUNC_ERR",
        CV_REPTD_SRHSFUNC_ERR => "CV_REPTD_SRHSFUNC_ERR",
        CV_UNREC_SRHSFUNC_ERR => "CV_UNREC_SRHSFUNC_ERR",
        CV_TOO_CLOSE => "CV_TOO_CLOSE",
        CV_NLS_INIT_FAIL => "CV_NLS_INIT_FAIL",
        CV_NLS_SETUP_FAIL => "CV_NLS_SETUPT_FAIL", /* (sic -- typo kept from the C source) */
        CV_NO_ADJ => "CV_NO_ADJ",
        CV_NO_FWD => "CV_NO_FWD",
        CV_NO_BCK => "CV_NO_BCK",
        CV_BAD_TB0 => "CV_BAD_TB0",
        CV_REIFWD_FAIL => "CV_REIFWD_FAIL",
        CV_FWD_FAIL => "CV_FWD_FAIL",
        CV_GETY_BADT => "CV_GETY_BADT",
        CV_NLS_FAIL => "CV_NLS_FAIL",
        CV_PROJ_MEM_NULL => "CV_PROJ_MEM_NULL",
        CV_PROJFUNC_FAIL => "CV_PROJFUNC_FAIL",
        CV_REPTD_PROJFUNC_ERR => "CV_REPTD_PROJFUNC_ERR",
        _ => "NONE",
    };

    name.to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cvodes::CVodeCreate;
    use crate::sundials_context::SUNContext_Create;

    fn fresh_mem() -> Box<CVodeMem> {
        let sunctx = SUNContext_Create();
        CVodeCreate(CV_BDF, &sunctx)
    }

    #[test]
    fn setter_validations() {
        let mut cv_mem = fresh_mem();

        /* CVodeSetMaxOrd: maxord <= 0 illegal */
        assert_eq!(CVodeSetMaxOrd(&mut cv_mem, 0), CV_ILL_INPUT);
        /* BDF: qmax_alloc = 5; cannot increase beyond it */
        assert_eq!(CVodeSetMaxOrd(&mut cv_mem, 6), CV_ILL_INPUT);
        /* CVODES-specific: the bound is min(qmax_alloc, qmax_allocQ,
           qmax_allocS) */
        cv_mem.cv_qmax_allocQ = 3;
        assert_eq!(CVodeSetMaxOrd(&mut cv_mem, 4), CV_ILL_INPUT);
        assert_eq!(CVodeSetMaxOrd(&mut cv_mem, 3), CV_SUCCESS);
        assert_eq!(cv_mem.cv_qmax, 3);

        /* CVodeSetEtaMin: out-of-range values select the default */
        assert_eq!(CVodeSetEtaMin(&mut cv_mem, 5.0), CV_SUCCESS);
        assert_eq!(cv_mem.cv_eta_min, ETA_MIN_DEFAULT);
        assert_eq!(CVodeSetEtaMin(&mut cv_mem, 0.5), CV_SUCCESS);
        assert_eq!(cv_mem.cv_eta_min, 0.5);

        /* CVodeSetMinStep: hmin < 0 illegal */
        assert_eq!(CVodeSetMinStep(&mut cv_mem, -1.0), CV_ILL_INPUT);

        /* CVodeSetStopTime: tstop behind tn in the direction of
           integration is illegal once a step was taken */
        cv_mem.cv_nst = 1;
        cv_mem.cv_tn = 1.0;
        cv_mem.cv_h = 1.0;
        assert_eq!(CVodeSetStopTime(&mut cv_mem, 0.0), CV_ILL_INPUT);
        assert_eq!(CVodeSetStopTime(&mut cv_mem, 2.0), CV_SUCCESS);
        assert!(cv_mem.cv_tstopset);

        /* CVodeSetSensDQMethod: bad DQtype / negative DQrhomax illegal */
        assert_eq!(CVodeSetSensDQMethod(&mut cv_mem, 0, 0.0), CV_ILL_INPUT);
        assert_eq!(CVodeSetSensDQMethod(&mut cv_mem, CV_CENTERED, -1.0), CV_ILL_INPUT);
        assert_eq!(CVodeSetSensDQMethod(&mut cv_mem, CV_FORWARD, 0.5), CV_SUCCESS);
        assert_eq!(cv_mem.cv_DQtype, CV_FORWARD);
        assert_eq!(cv_mem.cv_DQrhomax, 0.5);

        /* CVodeSetSensParams / CVodeSetQuadSensErrCon before
           sensitivity initialization */
        assert_eq!(CVodeSetSensParams(&mut cv_mem, None, None, None), CV_NO_SENS);
        assert_eq!(CVodeSetQuadSensErrCon(&mut cv_mem, true), CV_NO_SENS);

        /* quadrature / sensitivity getters before initialization */
        let mut n = 0i64;
        assert_eq!(CVodeGetQuadNumRhsEvals(&mut cv_mem, &mut n), CV_NO_QUAD);
        assert_eq!(CVodeGetSensNumRhsEvals(&mut cv_mem, &mut n), CV_NO_SENS);
    }

    #[test]
    fn return_flag_names() {
        assert_eq!(CVodeGetReturnFlagName(0), "CV_SUCCESS");
        assert_eq!(CVodeGetReturnFlagName(-14), "CV_NLS_SETUPT_FAIL"); /* C typo kept */
        assert_eq!(CVodeGetReturnFlagName(-30), "CV_NO_QUAD");
        assert_eq!(CVodeGetReturnFlagName(-40), "CV_NO_SENS");
        assert_eq!(CVodeGetReturnFlagName(-101), "CV_NO_ADJ");
        assert_eq!(CVodeGetReturnFlagName(-107), "CV_GETY_BADT");
        /* the CVODES switch has no case for the quad-sens codes */
        assert_eq!(CVodeGetReturnFlagName(-50), "NONE");
        assert_eq!(CVodeGetReturnFlagName(12345), "NONE");
    }

    /* Fill a Default-free CVodeMem (built via CVodeCreate) with known
       statistics that exercise the CVODES quad / sens (staggered1) /
       quad-sens sections of CVodePrintAllStats. */
    fn stats_mem() -> Box<CVodeMem> {
        let mut cv_mem = fresh_mem();
        cv_mem.cv_tn = 0.5;
        cv_mem.cv_nst = 4;
        cv_mem.cv_netf = 1;
        cv_mem.cv_ncfn = 2;
        cv_mem.constraint_fails = 0;
        cv_mem.constraint_corrections = 0;
        cv_mem.cv_h0u = 0.001;
        cv_mem.cv_hu = 0.25;
        cv_mem.cv_next_h = 0.5;
        cv_mem.cv_qu = 2;
        cv_mem.cv_next_q = 3;
        cv_mem.cv_nor = 0;
        cv_mem.cv_nfe = 10;
        cv_mem.cv_nni = 8;
        cv_mem.cv_nnf = 1;
        cv_mem.cv_nsetups = 3;
        cv_mem.cv_nge = 0;
        /* quadrature */
        cv_mem.cv_quadr = true;
        cv_mem.cv_nfQe = 6;
        cv_mem.cv_netfQ = 1;
        /* sensitivities, staggered1 */
        cv_mem.cv_sensi = true;
        cv_mem.cv_ism = CV_STAGGERED1;
        cv_mem.cv_Ns = 2;
        cv_mem.cv_nfSe = 12;
        cv_mem.cv_nfeS = 24;
        cv_mem.cv_netfS = 2;
        cv_mem.cv_nniS = 16;
        cv_mem.cv_nnfS = 1;
        cv_mem.cv_ncfnS = 1;
        cv_mem.cv_nniS1 = vec![7, 9];
        cv_mem.cv_nnfS1 = vec![1, 0];
        cv_mem.cv_ncfnS1 = vec![0, 1];
        cv_mem.cv_nsetupsS = 5;
        /* quadrature sensitivities */
        cv_mem.cv_quadr_sensi = true;
        cv_mem.cv_nfQSe = 4;
        cv_mem.cv_netfQS = 0;
        cv_mem
    }

    #[test]
    fn print_all_stats_table() {
        let mut cv_mem = stats_mem();
        let mut buf: Vec<u8> = Vec::new();
        assert_eq!(
            CVodePrintAllStats(&mut cv_mem, &mut buf, SUN_OUTPUTFORMAT_TABLE),
            CV_SUCCESS
        );
        let out = String::from_utf8(buf).unwrap();

        /* Reference: sunfprintf_* with SUN_TABLE_WIDTH = 29
           ("%-29s = ..."); values are independent literals. */
        fn line(name: &str, val: &str) -> String {
            format!("{:<29} = {}\n", name, val)
        }
        let expected: String = [
            line("Current time", "0.5"),
            line("Steps", "4"),
            line("Error test fails", "1"),
            line("NLS step fails", "2"),
            line("Constraint fails", "0"),
            line("Constraint corrections", "0"),
            line("Initial step size", "0.001"),
            line("Last step size", "0.25"),
            line("Current step size", "0.5"),
            line("Last method order", "2"),
            line("Current method order", "3"),
            line("Stab. lim. order reductions", "0"),
            line("RHS fn evals", "10"),
            line("NLS iters", "8"),
            line("NLS fails", "1"),
            line("NLS iters per step", "2"),
            line("LS setups", "3"),
            line("Root fn evals", "0"),
            line("Quad fn evals", "6"),
            line("Quad error test fails", "1"),
            line("Sens fn evals", "12"),
            line("Sens RHS fn evals", "24"),
            line("Sens error test fails", "2"),
            line("Sens NLS iters", "16"),
            line("Sens NLS fails", "1"),
            line("Sens NLS step fails", "1"),
            line("Sens stgr1 NLS iters", "7, 9"),
            line("Sens stgr1 NLS fails", "1, 0"),
            line("Sens stgr1 NLS step fails", "0, 1"),
            line("Sens LS setups", "5"),
            line("QuadSens fn evals", "4"),
            line("QuadSens error test fails", "0"),
        ]
        .concat();

        assert_eq!(out, expected);
    }

    #[test]
    fn print_all_stats_csv() {
        let mut cv_mem = stats_mem();
        let mut buf: Vec<u8> = Vec::new();
        assert_eq!(
            CVodePrintAllStats(&mut cv_mem, &mut buf, SUN_OUTPUTFORMAT_CSV),
            CV_SUCCESS
        );
        let out = String::from_utf8(buf).unwrap();

        /* C reference: "%s," SUN_FORMAT_E for reals ("% .15e") and
           "%s,%ld" for longs, comma-joined; the long-array CSV form
           prints "name i,value" groups with no separator between the
           entries of one array (as in sundials_utils.h). */
        let expected = concat!(
            "Current time, 5.000000000000000e-01",
            ",Steps,4",
            ",Error test fails,1",
            ",NLS step fails,2",
            ",Constraint fails,0",
            ",Constraint corrections,0",
            ",Initial step size, 1.000000000000000e-03",
            ",Last step size, 2.500000000000000e-01",
            ",Current step size, 5.000000000000000e-01",
            ",Last method order,2",
            ",Current method order,3",
            ",Stab. lim. order reductions,0",
            ",RHS fn evals,10",
            ",NLS iters,8",
            ",NLS fails,1",
            ",NLS iters per step, 2.000000000000000e+00",
            ",LS setups,3",
            ",Root fn evals,0",
            ",Quad fn evals,6",
            ",Quad error test fails,1",
            ",Sens fn evals,12",
            ",Sens RHS fn evals,24",
            ",Sens error test fails,2",
            ",Sens NLS iters,16",
            ",Sens NLS fails,1",
            ",Sens NLS step fails,1",
            ",Sens stgr1 NLS iters 0,7Sens stgr1 NLS iters 1,9",
            ",Sens stgr1 NLS fails 0,1Sens stgr1 NLS fails 1,0",
            ",Sens stgr1 NLS step fails 0,0Sens stgr1 NLS step fails 1,1",
            ",Sens LS setups,5",
            ",QuadSens fn evals,4",
            ",QuadSens error test fails,0",
        );

        assert_eq!(out, expected);
    }

    #[test]
    fn print_all_stats_both_formats_on_fresh_mem() {
        /* The C checks fmt against both enumerators; with the
           two-variant Rust enum that guard cannot fail, but it is
           kept line-for-line. Exercise both valid paths on a fresh
           memory (no quad/sens sections). */
        let mut cv_mem = fresh_mem();
        let mut sink: Vec<u8> = Vec::new();
        assert_eq!(
            CVodePrintAllStats(&mut cv_mem, &mut sink, SUN_OUTPUTFORMAT_TABLE),
            CV_SUCCESS
        );
        assert_eq!(
            CVodePrintAllStats(&mut cv_mem, &mut sink, SUN_OUTPUTFORMAT_CSV),
            CV_SUCCESS
        );
    }
}
