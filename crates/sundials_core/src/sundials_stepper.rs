/* -----------------------------------------------------------------
 * Translation of sundials-7.7.0/src/sundials/sundials_stepper.c
 * (+ include/sundials/sundials_stepper.h and
 *  src/sundials/sundials_stepper_impl.h).
 *
 * SUNStepper is the generic "evolve an IVP" vtable object used by
 * ARKODE's SplittingStep/ForcingStep/adjoint machinery. Unlike the
 * other core base classes (which became enum dispatch), SUNStepper's
 * ops are registered one at a time by the wrapping integrator via the
 * SUNStepper_Set*Fn functions, so the ops table stays a table: each
 * op is an `Option<fn>` field taking `&mut SUNStepper`, and each
 * dispatcher returns SUN_ERR_NOT_IMPLEMENTED when the op is None,
 * exactly as C does for a NULL pointer.
 *
 * `content` (C `void*`, the wrapped integrator memory) is the usual
 * `Option<Box<dyn Any>>`; op implementations downcast it. C's
 * SUNStepper_GetContent(stepper, &ptr) copies the raw pointer out —
 * inexpressible here, so it returns `&mut UserData` instead.
 * The `python` field and SUNDIALS_ENABLE_PYTHON block are excluded
 * (python bindings are outside this port), and `sunctx` is carried
 * only by the constructor signature like the other core objects.
 *
 * Note: C's SUNStepper_Create leaves ops->reinit,
 * ops->resetcheckpointindex and ops->getnumsteps UNINITIALIZED
 * (malloc'd, never NULLed — a latent C bug). Rust initializes every
 * op to None.
 * -----------------------------------------------------------------*/

use crate::nvector_serial::NVector;
use crate::sundials_context::SUNContext;
use crate::sundials_errors::{SUNErrCode, SUN_ERR_NOT_IMPLEMENTED, SUN_SUCCESS};
use crate::sundials_matrix::SUNMatrix;
use crate::sundials_types::{suncountertype, sunrealtype, UserData};

/* -----------------------------------------------------------------
 * Types from include/sundials/sundials_stepper.h
 * ----------------------------------------------------------------- */

/// enum SUNFullRhsMode
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SUNFullRhsMode {
    SUN_FULLRHS_START,
    SUN_FULLRHS_END,
    SUN_FULLRHS_OTHER,
}
pub use SUNFullRhsMode::{SUN_FULLRHS_END, SUN_FULLRHS_OTHER, SUN_FULLRHS_START};

/// SUNRhsJacFn (declared in sundials_stepper.h; no in-tree C caller)
pub type SUNRhsJacFn = fn(
    t: sunrealtype,
    y: &NVector,
    fy: &NVector,
    jac: &mut SUNMatrix,
    user_data: &mut UserData,
    tmp1: &mut NVector,
    tmp2: &mut NVector,
    tmp3: &mut NVector,
) -> i32;

/// SUNRhsJacTimesFn (declared in sundials_stepper.h; no in-tree C caller)
pub type SUNRhsJacTimesFn = fn(
    v: &NVector,
    jv: &mut NVector,
    t: sunrealtype,
    y: &NVector,
    fy: &NVector,
    user_data: &mut UserData,
    tmp: &mut NVector,
) -> i32;

pub type SUNStepperEvolveFn =
    fn(stepper: &mut SUNStepper, tout: sunrealtype, vret: &mut NVector, tret: &mut sunrealtype) -> SUNErrCode;

pub type SUNStepperOneStepFn =
    fn(stepper: &mut SUNStepper, tout: sunrealtype, vret: &mut NVector, tret: &mut sunrealtype) -> SUNErrCode;

pub type SUNStepperFullRhsFn = fn(
    stepper: &mut SUNStepper,
    t: sunrealtype,
    v: &NVector,
    f: &mut NVector,
    mode: SUNFullRhsMode,
) -> SUNErrCode;

pub type SUNStepperReInitFn =
    fn(stepper: &mut SUNStepper, t0: sunrealtype, v0: &NVector) -> SUNErrCode;

pub type SUNStepperResetFn =
    fn(stepper: &mut SUNStepper, tR: sunrealtype, vR: &NVector) -> SUNErrCode;

pub type SUNStepperResetCheckpointIndexFn =
    fn(stepper: &mut SUNStepper, ckptIdxR: suncountertype) -> SUNErrCode;

pub type SUNStepperSetStopTimeFn =
    fn(stepper: &mut SUNStepper, tstop: sunrealtype) -> SUNErrCode;

pub type SUNStepperSetStepDirectionFn =
    fn(stepper: &mut SUNStepper, stepdir: sunrealtype) -> SUNErrCode;

pub type SUNStepperSetForcingFn = fn(
    stepper: &mut SUNStepper,
    tshift: sunrealtype,
    tscale: sunrealtype,
    forcing: &[NVector],
    nforcing: i32,
) -> SUNErrCode;

pub type SUNStepperGetNumStepsFn =
    fn(stepper: &mut SUNStepper, nst: &mut suncountertype) -> SUNErrCode;

pub type SUNStepperDestroyFn = fn(stepper: &mut SUNStepper) -> SUNErrCode;

/* -----------------------------------------------------------------
 * struct SUNStepper_Ops_ / struct SUNStepper_ (sundials_stepper_impl.h)
 * ----------------------------------------------------------------- */

#[derive(Default)]
pub struct SUNStepper_Ops {
    pub evolve: Option<SUNStepperEvolveFn>,
    pub onestep: Option<SUNStepperOneStepFn>,
    pub fullrhs: Option<SUNStepperFullRhsFn>,
    pub reinit: Option<SUNStepperReInitFn>,
    pub reset: Option<SUNStepperResetFn>,
    pub resetcheckpointindex: Option<SUNStepperResetCheckpointIndexFn>,
    pub setstoptime: Option<SUNStepperSetStopTimeFn>,
    pub setstepdirection: Option<SUNStepperSetStepDirectionFn>,
    pub setforcing: Option<SUNStepperSetForcingFn>,
    pub getnumsteps: Option<SUNStepperGetNumStepsFn>,
    pub destroy: Option<SUNStepperDestroyFn>,
}

pub struct SUNStepper {
    /// stepper specific content (C `void* content`)
    pub content: UserData,
    /// stepper operations
    pub ops: SUNStepper_Ops,
    /// last stepper return flag
    pub last_flag: i32,
}

/* -----------------------------------------------------------------
 * sundials_stepper.c
 * ----------------------------------------------------------------- */

/// C SUNStepper_Create(sunctx, &stepper): allocation cannot fail, so
/// the object is returned directly.
pub fn SUNStepper_Create(_sunctx: &SUNContext) -> SUNStepper {
    SUNStepper {
        content: None,
        ops: SUNStepper_Ops::default(),
        last_flag: SUN_SUCCESS,
    }
}

/// C SUNStepper_Destroy: invokes ops->destroy (content cleanup), then
/// frees; the free is ownership drop in Rust.
pub fn SUNStepper_Destroy(stepper: &mut Option<SUNStepper>) -> SUNErrCode {
    if let Some(s) = stepper.as_mut() {
        if let Some(destroy) = s.ops.destroy {
            destroy(s);
        }
        *stepper = None;
    }
    SUN_SUCCESS
}

pub fn SUNStepper_Evolve(
    stepper: &mut SUNStepper,
    tout: sunrealtype,
    y: &mut NVector,
    tret: &mut sunrealtype,
) -> SUNErrCode {
    if let Some(evolve) = stepper.ops.evolve {
        return evolve(stepper, tout, y, tret);
    }
    SUN_ERR_NOT_IMPLEMENTED
}

pub fn SUNStepper_OneStep(
    stepper: &mut SUNStepper,
    tout: sunrealtype,
    y: &mut NVector,
    tret: &mut sunrealtype,
) -> SUNErrCode {
    if let Some(onestep) = stepper.ops.onestep {
        return onestep(stepper, tout, y, tret);
    }
    SUN_ERR_NOT_IMPLEMENTED
}

pub fn SUNStepper_FullRhs(
    stepper: &mut SUNStepper,
    t: sunrealtype,
    v: &NVector,
    f: &mut NVector,
    mode: SUNFullRhsMode,
) -> SUNErrCode {
    if let Some(fullrhs) = stepper.ops.fullrhs {
        return fullrhs(stepper, t, v, f, mode);
    }
    SUN_ERR_NOT_IMPLEMENTED
}

pub fn SUNStepper_ReInit(
    stepper: &mut SUNStepper,
    t0: sunrealtype,
    y0: &NVector,
) -> SUNErrCode {
    if let Some(reinit) = stepper.ops.reinit {
        return reinit(stepper, t0, y0);
    }
    SUN_ERR_NOT_IMPLEMENTED
}

pub fn SUNStepper_Reset(stepper: &mut SUNStepper, tR: sunrealtype, yR: &NVector) -> SUNErrCode {
    if let Some(reset) = stepper.ops.reset {
        return reset(stepper, tR, yR);
    }
    SUN_ERR_NOT_IMPLEMENTED
}

pub fn SUNStepper_ResetCheckpointIndex(
    stepper: &mut SUNStepper,
    ckptIdxR: suncountertype,
) -> SUNErrCode {
    if let Some(resetcheckpointindex) = stepper.ops.resetcheckpointindex {
        return resetcheckpointindex(stepper, ckptIdxR);
    }
    SUN_ERR_NOT_IMPLEMENTED
}

pub fn SUNStepper_SetStopTime(stepper: &mut SUNStepper, tstop: sunrealtype) -> SUNErrCode {
    if let Some(setstoptime) = stepper.ops.setstoptime {
        return setstoptime(stepper, tstop);
    }
    SUN_ERR_NOT_IMPLEMENTED
}

pub fn SUNStepper_SetStepDirection(stepper: &mut SUNStepper, stepdir: sunrealtype) -> SUNErrCode {
    if let Some(setstepdirection) = stepper.ops.setstepdirection {
        return setstepdirection(stepper, stepdir);
    }
    SUN_ERR_NOT_IMPLEMENTED
}

pub fn SUNStepper_SetForcing(
    stepper: &mut SUNStepper,
    tshift: sunrealtype,
    tscale: sunrealtype,
    forcing: &[NVector],
    nforcing: i32,
) -> SUNErrCode {
    if let Some(setforcing) = stepper.ops.setforcing {
        return setforcing(stepper, tshift, tscale, forcing, nforcing);
    }
    SUN_ERR_NOT_IMPLEMENTED
}

pub fn SUNStepper_SetContent(stepper: &mut SUNStepper, content: UserData) -> SUNErrCode {
    stepper.content = content;
    SUN_SUCCESS
}

/// C SUNStepper_GetContent(stepper, &content) copies the `void*` out;
/// the Rust port hands back a mutable borrow of the content slot.
pub fn SUNStepper_GetContent(stepper: &mut SUNStepper) -> &mut UserData {
    &mut stepper.content
}

pub fn SUNStepper_GetNumSteps(stepper: &mut SUNStepper, nst: &mut suncountertype) -> SUNErrCode {
    if let Some(getnumsteps) = stepper.ops.getnumsteps {
        return getnumsteps(stepper, nst);
    }
    SUN_ERR_NOT_IMPLEMENTED
}

pub fn SUNStepper_SetLastFlag(stepper: &mut SUNStepper, last_flag: i32) -> SUNErrCode {
    stepper.last_flag = last_flag;
    SUN_SUCCESS
}

pub fn SUNStepper_GetLastFlag(stepper: &SUNStepper, last_flag: &mut i32) -> SUNErrCode {
    *last_flag = stepper.last_flag;
    SUN_SUCCESS
}

pub fn SUNStepper_SetEvolveFn(stepper: &mut SUNStepper, f: SUNStepperEvolveFn) -> SUNErrCode {
    stepper.ops.evolve = Some(f);
    SUN_SUCCESS
}

pub fn SUNStepper_SetOneStepFn(stepper: &mut SUNStepper, f: SUNStepperOneStepFn) -> SUNErrCode {
    stepper.ops.onestep = Some(f);
    SUN_SUCCESS
}

pub fn SUNStepper_SetFullRhsFn(stepper: &mut SUNStepper, f: SUNStepperFullRhsFn) -> SUNErrCode {
    stepper.ops.fullrhs = Some(f);
    SUN_SUCCESS
}

pub fn SUNStepper_SetReInitFn(stepper: &mut SUNStepper, f: SUNStepperReInitFn) -> SUNErrCode {
    stepper.ops.reinit = Some(f);
    SUN_SUCCESS
}

pub fn SUNStepper_SetResetFn(stepper: &mut SUNStepper, f: SUNStepperResetFn) -> SUNErrCode {
    stepper.ops.reset = Some(f);
    SUN_SUCCESS
}

pub fn SUNStepper_SetResetCheckpointIndexFn(
    stepper: &mut SUNStepper,
    f: SUNStepperResetCheckpointIndexFn,
) -> SUNErrCode {
    stepper.ops.resetcheckpointindex = Some(f);
    SUN_SUCCESS
}

pub fn SUNStepper_SetStopTimeFn(
    stepper: &mut SUNStepper,
    f: SUNStepperSetStopTimeFn,
) -> SUNErrCode {
    stepper.ops.setstoptime = Some(f);
    SUN_SUCCESS
}

pub fn SUNStepper_SetStepDirectionFn(
    stepper: &mut SUNStepper,
    f: SUNStepperSetStepDirectionFn,
) -> SUNErrCode {
    stepper.ops.setstepdirection = Some(f);
    SUN_SUCCESS
}

pub fn SUNStepper_SetForcingFn(stepper: &mut SUNStepper, f: SUNStepperSetForcingFn) -> SUNErrCode {
    stepper.ops.setforcing = Some(f);
    SUN_SUCCESS
}

pub fn SUNStepper_SetGetNumStepsFn(
    stepper: &mut SUNStepper,
    f: SUNStepperGetNumStepsFn,
) -> SUNErrCode {
    stepper.ops.getnumsteps = Some(f);
    SUN_SUCCESS
}

pub fn SUNStepper_SetDestroyFn(stepper: &mut SUNStepper, f: SUNStepperDestroyFn) -> SUNErrCode {
    stepper.ops.destroy = Some(f);
    SUN_SUCCESS
}

/* ----------------------------------------------------------------- */

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nvector_serial::{N_VConst, N_VNew_Serial};
    use crate::sundials_context::SUNContext_Create;

    fn dummy_evolve(
        stepper: &mut SUNStepper,
        tout: sunrealtype,
        vret: &mut NVector,
        tret: &mut sunrealtype,
    ) -> SUNErrCode {
        N_VConst(tout, vret);
        *tret = tout;
        stepper.last_flag = 42;
        SUN_SUCCESS
    }

    fn dummy_getnumsteps(_stepper: &mut SUNStepper, nst: &mut suncountertype) -> SUNErrCode {
        *nst = 7;
        SUN_SUCCESS
    }

    fn dummy_destroy(stepper: &mut SUNStepper) -> SUNErrCode {
        stepper.content = None;
        SUN_SUCCESS
    }

    #[test]
    fn test_create_ops_empty() {
        let ctx = SUNContext_Create();
        let mut stepper = SUNStepper_Create(&ctx);
        assert_eq!(stepper.last_flag, SUN_SUCCESS);

        let mut y = N_VNew_Serial(3, &ctx);
        let mut tret = 0.0;
        assert_eq!(
            SUNStepper_Evolve(&mut stepper, 1.0, &mut y, &mut tret),
            SUN_ERR_NOT_IMPLEMENTED
        );
        assert_eq!(
            SUNStepper_OneStep(&mut stepper, 1.0, &mut y, &mut tret),
            SUN_ERR_NOT_IMPLEMENTED
        );
        let v = N_VNew_Serial(3, &ctx);
        let mut f = N_VNew_Serial(3, &ctx);
        assert_eq!(
            SUNStepper_FullRhs(&mut stepper, 0.0, &v, &mut f, SUN_FULLRHS_START),
            SUN_ERR_NOT_IMPLEMENTED
        );
        assert_eq!(SUNStepper_ReInit(&mut stepper, 0.0, &v), SUN_ERR_NOT_IMPLEMENTED);
        assert_eq!(SUNStepper_Reset(&mut stepper, 0.0, &v), SUN_ERR_NOT_IMPLEMENTED);
        assert_eq!(
            SUNStepper_ResetCheckpointIndex(&mut stepper, 0),
            SUN_ERR_NOT_IMPLEMENTED
        );
        assert_eq!(SUNStepper_SetStopTime(&mut stepper, 1.0), SUN_ERR_NOT_IMPLEMENTED);
        assert_eq!(
            SUNStepper_SetStepDirection(&mut stepper, 1.0),
            SUN_ERR_NOT_IMPLEMENTED
        );
        assert_eq!(
            SUNStepper_SetForcing(&mut stepper, 0.0, 1.0, &[], 0),
            SUN_ERR_NOT_IMPLEMENTED
        );
        let mut nst: suncountertype = 0;
        assert_eq!(
            SUNStepper_GetNumSteps(&mut stepper, &mut nst),
            SUN_ERR_NOT_IMPLEMENTED
        );
    }

    #[test]
    fn test_registered_ops_dispatch() {
        let ctx = SUNContext_Create();
        let mut stepper = SUNStepper_Create(&ctx);
        assert_eq!(SUNStepper_SetEvolveFn(&mut stepper, dummy_evolve), SUN_SUCCESS);
        assert_eq!(
            SUNStepper_SetGetNumStepsFn(&mut stepper, dummy_getnumsteps),
            SUN_SUCCESS
        );

        let mut y = N_VNew_Serial(2, &ctx);
        let mut tret = 0.0;
        assert_eq!(SUNStepper_Evolve(&mut stepper, 2.5, &mut y, &mut tret), SUN_SUCCESS);
        assert_eq!(tret, 2.5);
        assert_eq!(y.data[0], 2.5);
        assert_eq!(stepper.last_flag, 42);

        let mut flag = 0;
        assert_eq!(SUNStepper_GetLastFlag(&stepper, &mut flag), SUN_SUCCESS);
        assert_eq!(flag, 42);
        assert_eq!(SUNStepper_SetLastFlag(&mut stepper, 0), SUN_SUCCESS);
        assert_eq!(stepper.last_flag, 0);

        let mut nst: suncountertype = 0;
        assert_eq!(SUNStepper_GetNumSteps(&mut stepper, &mut nst), SUN_SUCCESS);
        assert_eq!(nst, 7);
    }

    #[test]
    fn test_content_and_destroy() {
        let ctx = SUNContext_Create();
        let mut stepper = SUNStepper_Create(&ctx);
        assert_eq!(
            SUNStepper_SetContent(&mut stepper, Some(Box::new(3.5f64))),
            SUN_SUCCESS
        );
        {
            let content = SUNStepper_GetContent(&mut stepper);
            let val = content
                .as_mut()
                .and_then(|c| c.downcast_mut::<f64>())
                .expect("content downcast");
            assert_eq!(*val, 3.5);
            *val = 4.5;
        }
        assert_eq!(SUNStepper_SetDestroyFn(&mut stepper, dummy_destroy), SUN_SUCCESS);

        let mut slot = Some(stepper);
        assert_eq!(SUNStepper_Destroy(&mut slot), SUN_SUCCESS);
        assert!(slot.is_none());
        /* destroying an already-empty slot is a no-op, like C's NULL check */
        assert_eq!(SUNStepper_Destroy(&mut slot), SUN_SUCCESS);
    }
}
