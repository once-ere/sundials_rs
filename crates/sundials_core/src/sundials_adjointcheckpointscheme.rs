/* -----------------------------------------------------------------
 * Translation of sundials-7.7.0/src/sundials/
 * sundials_adjointcheckpointscheme.c
 * (+ include/sundials/sundials_adjointcheckpointscheme.h and
 *  src/sundials/sundials_adjointcheckpointscheme_impl.h).
 *
 * Like SUNStepper, the ops are registered one at a time via the
 * Set*Fn functions, so the ops table stays an Option<fn> table and
 * each dispatcher returns SUN_ERR_NOT_IMPLEMENTED when the op is
 * None. `content` is the usual Option<Box<dyn Any>> (implementations
 * downcast); GetContent returns &mut UserData since C copies the raw
 * pointer out. SUNAdjointCheckpointScheme_NewEmpty is the
 * constructor (allocation cannot fail, so it returns the object).
 * -----------------------------------------------------------------*/

use crate::nvector_serial::NVector;
use crate::sundials_context::SUNContext;
use crate::sundials_errors::{SUNErrCode, SUN_ERR_NOT_IMPLEMENTED, SUN_SUCCESS};
use crate::sundials_types::{suncountertype, sunbooleantype, sunrealtype, UserData};

pub type SUNAdjointCheckpointSchemeNeedsSavingFn = fn(
    check_scheme: &mut SUNAdjointCheckpointScheme,
    step_num: suncountertype,
    stage_num: suncountertype,
    t: sunrealtype,
    yes_or_no: &mut sunbooleantype,
) -> SUNErrCode;

pub type SUNAdjointCheckpointSchemeInsertVectorFn = fn(
    check_scheme: &mut SUNAdjointCheckpointScheme,
    step_num: suncountertype,
    stage_num: suncountertype,
    t: sunrealtype,
    y: &NVector,
) -> SUNErrCode;

/// C: `N_Vector* yout` (some schemes may hand back a reference); the
/// in-tree scheme fills the caller's vector, so `&mut NVector` here.
pub type SUNAdjointCheckpointSchemeLoadVectorFn = fn(
    check_scheme: &mut SUNAdjointCheckpointScheme,
    step_num: suncountertype,
    stage_num: suncountertype,
    peek: sunbooleantype,
    yout: &mut NVector,
    tout: &mut sunrealtype,
) -> SUNErrCode;

pub type SUNAdjointCheckpointSchemeDestroyFn =
    fn(check_scheme: &mut Option<SUNAdjointCheckpointScheme>) -> SUNErrCode;

pub type SUNAdjointCheckpointSchemeEnableDenseFn = fn(
    check_scheme: &mut SUNAdjointCheckpointScheme,
    on_or_off: sunbooleantype,
) -> SUNErrCode;

/// struct SUNAdjointCheckpointScheme_Ops_ (impl header)
#[derive(Default)]
pub struct SUNAdjointCheckpointScheme_Ops {
    pub needssaving: Option<SUNAdjointCheckpointSchemeNeedsSavingFn>,
    pub insertvector: Option<SUNAdjointCheckpointSchemeInsertVectorFn>,
    pub loadvector: Option<SUNAdjointCheckpointSchemeLoadVectorFn>,
    pub destroy: Option<SUNAdjointCheckpointSchemeDestroyFn>,
    pub enableDense: Option<SUNAdjointCheckpointSchemeEnableDenseFn>,
}

/// struct SUNAdjointCheckpointScheme_ (impl header)
pub struct SUNAdjointCheckpointScheme {
    pub ops: SUNAdjointCheckpointScheme_Ops,
    pub content: UserData,
}

/// C SUNAdjointCheckpointScheme_NewEmpty(sunctx, &scheme).
pub fn SUNAdjointCheckpointScheme_NewEmpty(_sunctx: &SUNContext) -> SUNAdjointCheckpointScheme {
    SUNAdjointCheckpointScheme {
        ops: SUNAdjointCheckpointScheme_Ops::default(),
        content: None,
    }
}

pub fn SUNAdjointCheckpointScheme_NeedsSaving(
    check_scheme: &mut SUNAdjointCheckpointScheme,
    step_num: suncountertype,
    stage_num: suncountertype,
    t: sunrealtype,
    yes_or_no: &mut sunbooleantype,
) -> SUNErrCode {
    if let Some(needssaving) = check_scheme.ops.needssaving {
        return needssaving(check_scheme, step_num, stage_num, t, yes_or_no);
    }
    SUN_ERR_NOT_IMPLEMENTED
}

pub fn SUNAdjointCheckpointScheme_InsertVector(
    check_scheme: &mut SUNAdjointCheckpointScheme,
    step_num: suncountertype,
    stage_num: suncountertype,
    t: sunrealtype,
    state: &NVector,
) -> SUNErrCode {
    if let Some(insertvector) = check_scheme.ops.insertvector {
        return insertvector(check_scheme, step_num, stage_num, t, state);
    }
    SUN_ERR_NOT_IMPLEMENTED
}

pub fn SUNAdjointCheckpointScheme_LoadVector(
    check_scheme: &mut SUNAdjointCheckpointScheme,
    step_num: suncountertype,
    stage_num: suncountertype,
    peek: sunbooleantype,
    out: &mut NVector,
    tout: &mut sunrealtype,
) -> SUNErrCode {
    if let Some(loadvector) = check_scheme.ops.loadvector {
        return loadvector(check_scheme, step_num, stage_num, peek, out, tout);
    }
    SUN_ERR_NOT_IMPLEMENTED
}

pub fn SUNAdjointCheckpointScheme_Destroy(
    check_scheme_ptr: &mut Option<SUNAdjointCheckpointScheme>,
) -> SUNErrCode {
    if let Some(scheme) = check_scheme_ptr.as_ref() {
        if let Some(destroy) = scheme.ops.destroy {
            return destroy(check_scheme_ptr);
        }
        /* no destroy op: free ops + object (drop) */
        *check_scheme_ptr = None;
    }
    SUN_SUCCESS
}

pub fn SUNAdjointCheckpointScheme_EnableDense(
    check_scheme: &mut SUNAdjointCheckpointScheme,
    on_or_off: sunbooleantype,
) -> SUNErrCode {
    if let Some(enableDense) = check_scheme.ops.enableDense {
        return enableDense(check_scheme, on_or_off);
    }
    SUN_ERR_NOT_IMPLEMENTED
}

pub fn SUNAdjointCheckpointScheme_SetContent(
    check_scheme: &mut SUNAdjointCheckpointScheme,
    content: UserData,
) -> SUNErrCode {
    check_scheme.content = content;
    SUN_SUCCESS
}

/// C SUNAdjointCheckpointScheme_GetContent(scheme, &content) copies
/// the `void*` out; the Rust port lends the content slot.
pub fn SUNAdjointCheckpointScheme_GetContent(
    check_scheme: &mut SUNAdjointCheckpointScheme,
) -> &mut UserData {
    &mut check_scheme.content
}

pub fn SUNAdjointCheckpointScheme_SetNeedsSavingFn(
    check_scheme: &mut SUNAdjointCheckpointScheme,
    f: SUNAdjointCheckpointSchemeNeedsSavingFn,
) -> SUNErrCode {
    check_scheme.ops.needssaving = Some(f);
    SUN_SUCCESS
}

pub fn SUNAdjointCheckpointScheme_SetInsertVectorFn(
    check_scheme: &mut SUNAdjointCheckpointScheme,
    f: SUNAdjointCheckpointSchemeInsertVectorFn,
) -> SUNErrCode {
    check_scheme.ops.insertvector = Some(f);
    SUN_SUCCESS
}

pub fn SUNAdjointCheckpointScheme_SetLoadVectorFn(
    check_scheme: &mut SUNAdjointCheckpointScheme,
    f: SUNAdjointCheckpointSchemeLoadVectorFn,
) -> SUNErrCode {
    check_scheme.ops.loadvector = Some(f);
    SUN_SUCCESS
}

pub fn SUNAdjointCheckpointScheme_SetDestroyFn(
    check_scheme: &mut SUNAdjointCheckpointScheme,
    f: SUNAdjointCheckpointSchemeDestroyFn,
) -> SUNErrCode {
    check_scheme.ops.destroy = Some(f);
    SUN_SUCCESS
}

pub fn SUNAdjointCheckpointScheme_SetEnableDenseFn(
    check_scheme: &mut SUNAdjointCheckpointScheme,
    f: SUNAdjointCheckpointSchemeEnableDenseFn,
) -> SUNErrCode {
    check_scheme.ops.enableDense = Some(f);
    SUN_SUCCESS
}
