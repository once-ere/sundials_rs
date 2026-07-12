/* -----------------------------------------------------------------
 * Translation of src/arkode/arkode_sunstepper.c (SUNDIALS 7.7.0).
 *
 * This is the implementation file for ARKODE's interfacing with
 * SUNStepper.
 *
 * Ownership adaptation (documented deviation): C stores the raw
 * `void* arkode_mem` pointer as the SUNStepper content while the
 * caller retains ownership.  In safe Rust the wrapped integrator is
 * MOVED into the stepper: ARKodeCreateSUNStepper takes the
 * `Box<ARKodeMem>` by value and stores it in `stepper.content`
 * (`Option<Box<dyn Any>>`).  Callers that need the inner integrator
 * afterwards (e.g. to print statistics) downcast the content back
 * out.  Consequently arkSUNStepperSelfDestruct — which in C frees an
 * internally-created integrator — simply drops the content here.
 *
 * The resetcheckpointindex op (C arkSUNStepperResetCheckpointIndex ->
 * ARKodeSetAdjointCheckpointIndex) is not installed: the adjoint
 * module is outside this port.  Ops that C would reach through a NULL
 * function pointer (step_fullrhs/step_setforcing missing) return
 * SUN_ERR_OP_FAIL instead of invoking undefined behavior.
 * -----------------------------------------------------------------*/

use crate::arkode::{ARKodeEvolve, ARKodeReset};
use crate::arkode_impl::{
    ARKodeMem, ARK_FULLRHS_END, ARK_FULLRHS_OTHER, ARK_FULLRHS_START, ARK_NORMAL, ARK_ONE_STEP,
    ARK_SUCCESS,
};
use crate::arkode_io::{ARKodeSetStepDirection, ARKodeSetStopTime};
use crate::nvector_serial::NVector;
use crate::sundials_errors::{SUNErrCode, SUN_ERR_OP_FAIL, SUN_SUCCESS};
use crate::sundials_stepper::{
    SUNFullRhsMode, SUNStepper, SUNStepper_SetContent, SUNStepper_SetEvolveFn,
    SUNStepper_SetFullRhsFn, SUNStepper_SetGetNumStepsFn, SUNStepper_SetOneStepFn,
    SUNStepper_SetResetFn, SUNStepper_SetStepDirectionFn, SUNStepper_SetStopTimeFn,
    SUNStepper_SetForcingFn, SUN_FULLRHS_END, SUN_FULLRHS_OTHER, SUN_FULLRHS_START,
};
use crate::sundials_types::{suncountertype, sunrealtype};

/// Shortcut to downcast the stepper content back to the wrapped
/// ARKodeMem (C SUNStepper_GetContent + cast).
fn arkstepper_mem(stepper: &mut SUNStepper) -> Option<&mut ARKodeMem> {
    stepper
        .content
        .as_mut()
        .and_then(|c| c.downcast_mut::<ARKodeMem>())
}

fn arkSUNStepperEvolveHelper(
    stepper: &mut SUNStepper,
    tout: sunrealtype,
    y: &mut NVector,
    tret: &mut sunrealtype,
    mode: i32,
) -> SUNErrCode {
    /* extract the ARKODE memory struct */
    let retval = match arkstepper_mem(stepper) {
        Some(ark_mem) => ARKodeEvolve(ark_mem, tout, y, tret, mode),
        None => return SUN_ERR_OP_FAIL,
    };

    /* evolve inner ODE */
    stepper.last_flag = retval;
    if retval < 0 {
        return SUN_ERR_OP_FAIL;
    }

    SUN_SUCCESS
}

fn arkSUNStepperEvolve(
    stepper: &mut SUNStepper,
    tout: sunrealtype,
    y: &mut NVector,
    tret: &mut sunrealtype,
) -> SUNErrCode {
    arkSUNStepperEvolveHelper(stepper, tout, y, tret, ARK_NORMAL)
}

fn arkSUNStepperOneStep(
    stepper: &mut SUNStepper,
    tout: sunrealtype,
    y: &mut NVector,
    tret: &mut sunrealtype,
) -> SUNErrCode {
    arkSUNStepperEvolveHelper(stepper, tout, y, tret, ARK_ONE_STEP)
}

/*------------------------------------------------------------------------------
  Implementation of SUNStepperFullRhsFn to compute the full inner
  (fast) ODE IVP RHS.
  ----------------------------------------------------------------------------*/
fn arkSUNStepperFullRhs(
    stepper: &mut SUNStepper,
    t: sunrealtype,
    y: &NVector,
    f: &mut NVector,
    mode: SUNFullRhsMode,
) -> SUNErrCode {
    let ark_mode = match mode {
        SUN_FULLRHS_START => ARK_FULLRHS_START,
        SUN_FULLRHS_END => ARK_FULLRHS_END,
        SUN_FULLRHS_OTHER => ARK_FULLRHS_OTHER,
    };

    /* extract the ARKODE memory struct */
    let retval = match arkstepper_mem(stepper) {
        Some(ark_mem) => match ark_mem.step_fullrhs {
            Some(step_fullrhs) => step_fullrhs(ark_mem, t, y, f, ark_mode),
            None => return SUN_ERR_OP_FAIL,
        },
        None => return SUN_ERR_OP_FAIL,
    };

    stepper.last_flag = retval;
    if retval != ARK_SUCCESS {
        return SUN_ERR_OP_FAIL;
    }

    SUN_SUCCESS
}

/*------------------------------------------------------------------------------
  Implementation of SUNStepperResetFn to reset the stepper state.
  ----------------------------------------------------------------------------*/
fn arkSUNStepperReset(stepper: &mut SUNStepper, tR: sunrealtype, yR: &NVector) -> SUNErrCode {
    /* extract the ARKODE memory struct */
    let retval = match arkstepper_mem(stepper) {
        Some(ark_mem) => ARKodeReset(ark_mem, tR, yR),
        None => return SUN_ERR_OP_FAIL,
    };

    stepper.last_flag = retval;
    if retval != ARK_SUCCESS {
        return SUN_ERR_OP_FAIL;
    }

    SUN_SUCCESS
}

/*------------------------------------------------------------------------------
  Implementation of SUNStepperStopTimeFn to set the tstop time
  ----------------------------------------------------------------------------*/
fn arkSUNStepperSetStopTime(stepper: &mut SUNStepper, tstop: sunrealtype) -> SUNErrCode {
    /* extract the ARKODE memory struct */
    let retval = match arkstepper_mem(stepper) {
        Some(ark_mem) => ARKodeSetStopTime(ark_mem, tstop),
        None => return SUN_ERR_OP_FAIL,
    };

    stepper.last_flag = retval;
    if retval != ARK_SUCCESS {
        return SUN_ERR_OP_FAIL;
    }

    SUN_SUCCESS
}

fn arkSUNStepperSetStepDirection(stepper: &mut SUNStepper, stepdir: sunrealtype) -> SUNErrCode {
    /* extract the ARKODE memory struct */
    let retval = match arkstepper_mem(stepper) {
        Some(ark_mem) => ARKodeSetStepDirection(ark_mem, stepdir),
        None => return SUN_ERR_OP_FAIL,
    };

    stepper.last_flag = retval;
    if retval != ARK_SUCCESS {
        return SUN_ERR_OP_FAIL;
    }

    SUN_SUCCESS
}

fn arkSUNStepperSetForcing(
    stepper: &mut SUNStepper,
    tshift: sunrealtype,
    tscale: sunrealtype,
    forcing: &[NVector],
    nforcing: i32,
) -> SUNErrCode {
    /* extract the ARKODE memory struct */
    let retval = match arkstepper_mem(stepper) {
        Some(ark_mem) => match ark_mem.step_setforcing {
            Some(step_setforcing) => step_setforcing(ark_mem, tshift, tscale, forcing, nforcing),
            None => return SUN_ERR_OP_FAIL,
        },
        None => return SUN_ERR_OP_FAIL,
    };

    stepper.last_flag = retval;
    if retval != ARK_SUCCESS {
        return SUN_ERR_OP_FAIL;
    }

    SUN_SUCCESS
}

/// C arkSUNStepperSelfDestruct: useful when we create a
/// ARKodeMem/SUNStepper internally and want it destroyed with the
/// SUNStepper.  Here the content owns the integrator, so dropping the
/// content is C's ARKodeFree.
pub fn arkSUNStepperSelfDestruct(stepper: &mut SUNStepper) -> SUNErrCode {
    stepper.content = None;
    SUN_SUCCESS
}

fn arkSUNStepperGetNumSteps(stepper: &mut SUNStepper, nst: &mut suncountertype) -> SUNErrCode {
    match arkstepper_mem(stepper) {
        Some(ark_mem) => {
            *nst = ark_mem.nst;
            SUN_SUCCESS
        }
        None => SUN_ERR_OP_FAIL,
    }
}

/// C ARKodeCreateSUNStepper(arkode_mem, &stepper): wraps an ARKODE
/// integrator as a SUNStepper.  The Rust port takes ownership of the
/// integrator (see module header); none of the C failure paths
/// (allocation, op registration) can fail here, so the stepper is
/// returned directly.
pub fn ARKodeCreateSUNStepper(arkode_mem: Box<ARKodeMem>) -> SUNStepper {
    /* the C ops-table registration sequence, kept call-for-call */
    let has_setforcing = arkode_mem.step_setforcing.is_some();

    let mut stepper = SUNStepper {
        content: None,
        ops: Default::default(),
        last_flag: SUN_SUCCESS,
    };

    let _ = SUNStepper_SetContent(&mut stepper, Some(arkode_mem));
    let _ = SUNStepper_SetEvolveFn(&mut stepper, arkSUNStepperEvolve);
    let _ = SUNStepper_SetOneStepFn(&mut stepper, arkSUNStepperOneStep);
    let _ = SUNStepper_SetFullRhsFn(&mut stepper, arkSUNStepperFullRhs);
    let _ = SUNStepper_SetResetFn(&mut stepper, arkSUNStepperReset);
    /* C: SUNStepper_SetResetCheckpointIndexFn(...) — adjoint module
    excluded from this port, op left unset */
    let _ = SUNStepper_SetStopTimeFn(&mut stepper, arkSUNStepperSetStopTime);
    let _ = SUNStepper_SetStepDirectionFn(&mut stepper, arkSUNStepperSetStepDirection);

    if has_setforcing {
        let _ = SUNStepper_SetForcingFn(&mut stepper, arkSUNStepperSetForcing);
    }

    let _ = SUNStepper_SetGetNumStepsFn(&mut stepper, arkSUNStepperGetNumSteps);

    stepper
}
