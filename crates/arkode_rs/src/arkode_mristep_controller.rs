/* -----------------------------------------------------------------------------
 * Translated from src/arkode/arkode_mristep_controller.c (ARKODE 7.7.0).
 * MRIStep's multirate adaptivity controller layer.
 *
 * Storage adaptation (pinned): C wraps the user's MRI-H-TOL controller
 * in a SUNAdaptController_MRIStep object whose content carries raw
 * ark_mem/step_mem back-pointers, and hands that wrapper to ARKODE.
 * Safe Rust cannot store those back-pointers, so the "wrapper" here
 * collapses onto storing the MRI-H-TOL controller itself in
 * hadapt_mem.hcontroller (with ownership, matching C's SUNTRUE): the
 * type-forwarding ops (gettype/reset/write/space) then act on the
 * stored controller directly through the generic dispatch, while the
 * two ops that need the MRIStep memory — estimatestep and updateh —
 * become the functions below, invoked from arkAdapt/arkCompleteStep
 * where ark_mem is in scope (the controller is dispatched there by
 * its SUN_ADAPTCONTROLLER_MRI_H_TOL type).
 * ---------------------------------------------------------------------------*/
use crate::arkode_impl::ARKodeMem;
use crate::arkode_mristep::mriStep_AccessStepMem;
use crate::sundials_adaptcontroller::{
    SUNAdaptController, SUNAdaptController_EstimateStepTol, SUNAdaptController_UpdateMRIHTol,
};
use crate::sundials_errors::{SUNErrCode, SUN_ERR_ARG_OUTOFRANGE, SUN_ERR_MEM_FAIL, SUN_SUCCESS};

/*--------------------------------------------
  MRIStep SUNAdaptController wrapper functions
  --------------------------------------------*/

/* C SUNAdaptController_EstimateStep_MRIStep (including the generic
   SUNAdaptController_EstimateStep argument checks that C applies to
   the wrapper object before invoking its op) */
pub fn SUNAdaptController_EstimateStep_MRIStep(
    ark_mem: &mut ARKodeMem,
    cmri: &mut SUNAdaptController,
    h: f64,
    p: i32,
    dsm: f64,
    hnew: &mut f64,
) -> SUNErrCode {
    if !h.is_finite() || p < 0 || dsm < 0.0 {
        return SUN_ERR_ARG_OUTOFRANGE;
    }
    *hnew = h; /* initialize output with identity */

    /* Shortcuts to ARKODE and MRIStep memory */
    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "SUNAdaptController_EstimateStep_MRIStep")
    {
        None => return SUN_ERR_MEM_FAIL,
        Some(sm) => sm,
    };

    /* Estimate slow stepsize from MRI controller */
    let retval = SUNAdaptController_EstimateStepTol(
        cmri,
        h,
        step_mem.inner_rtol_factor,
        p,
        dsm,
        step_mem.inner_dsm,
        hnew,
        &mut step_mem.inner_rtol_factor_new,
    );
    ark_mem.step_mem = Some(step_mem);
    retval
}

/* C SUNAdaptController_UpdateH_MRIStep (including the generic
   SUNAdaptController_UpdateH argument checks) */
pub fn SUNAdaptController_UpdateH_MRIStep(
    ark_mem: &mut ARKodeMem,
    cmri: &mut SUNAdaptController,
    h: f64,
    dsm: f64,
) -> SUNErrCode {
    if !h.is_finite() || dsm < 0.0 {
        return SUN_ERR_ARG_OUTOFRANGE;
    }

    /* Shortcuts to ARKODE and MRIStep memory */
    let mut step_mem =
        match mriStep_AccessStepMem(ark_mem, "SUNAdaptController_UpdateH_MRIStep") {
            None => return SUN_ERR_MEM_FAIL,
            Some(sm) => sm,
        };

    /* Update MRI controller */
    let retval = SUNAdaptController_UpdateMRIHTol(
        cmri,
        h,
        step_mem.inner_rtol_factor,
        dsm,
        step_mem.inner_dsm,
    );
    if retval != SUN_SUCCESS {
        ark_mem.step_mem = Some(step_mem);
        return retval;
    }

    /* Update inner controller parameter to most-recent prediction */
    step_mem.inner_rtol_factor = step_mem.inner_rtol_factor_new;
    ark_mem.step_mem = Some(step_mem);

    /* return with success */
    SUN_SUCCESS
}

/* (C's GetType/Reset/Write/Space forwarders collapse: the stored
   controller IS the wrapped MRI-H-TOL controller, so the generic
   SUNAdaptController dispatch already reaches it.) */
