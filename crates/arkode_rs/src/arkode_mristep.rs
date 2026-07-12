/* -----------------------------------------------------------------
 * Translation of src/arkode/arkode_mristep.c (SUNDIALS 7.7.0).
 *
 * This is the implementation file for ARKODE's MRI time stepper
 * module, including the MRIStepInnerStepper base-class functions.
 *
 * step_mem follows the crate's take/downcast/put-back Access
 * pattern; TakeStep* additionally mem::take the ark_mem vectors
 * they pass to the inner stepper and restore them on every return
 * path.  Xvecs operand lists are assembled at call sites.  The
 * MRIStepInnerStepper owns its content (the wrapped inner
 * integrator); C keeps a shared pointer.
 *
 * Deferred with their subsystems (documented deviations):
 *  - mriStep_Resize (ARKodeResize itself is not yet ported)
 *  - the MRI-H-TOL SUNAdaptController wrapper
 *    (arkode_mristep_controller.c stores raw ark_mem/step_mem
 *    back-pointers that safe Rust cannot express;
 *    mriStep_SetAdaptController reports MRI-H-TOL controllers as
 *    unsupported until that layer is redesigned)
 * -----------------------------------------------------------------*/

use crate::arkode::{
    arkAllocVec, arkAllocVecArray, arkCreate, arkFreeVec, arkFreeVecArray, arkHandleFailure,
    arkInit, arkPredict_Bootstrap, arkPredict_CutoffOrder, arkPredict_MaximumOrder,
    arkPredict_VariableOrder,
};
use crate::arkode_impl::{
    arkProcessError, ARKLinsolFreeFn, ARKLinsolInitFn, ARKLinsolSetupFn, ARKLinsolSolveFn,
    ARKRhsFn, ARKodeMem, ARK_FULLRHS_END, ARK_FULLRHS_OTHER,
    ARK_FULLRHS_START, ARK_ILL_INPUT, ARK_INNERSTEP_FAIL, ARK_INNERTOOUTER_FAIL, ARK_INTERP_NONE,
    ARK_INVALID_TABLE, ARK_LINIT_FAIL, ARK_MEM_FAIL, ARK_MEM_NULL, ARK_NLS_INIT_FAIL, ARK_OUTERTOINNER_FAIL,
    ARK_POSTPROCESS_STAGE_FAIL, ARK_POSTPROCESS_STEP_FAIL, ARK_PRERHSFN_FAIL, ARK_RHSFUNC_FAIL,
    ARK_SUCCESS, ARK_TOO_CLOSE, ARK_UNREC_RHSFUNC_ERR, ARK_USER_PREDICT_FAIL, CONV_FAIL, FIRST_INIT, H0_BIAS, H0_UBFACTOR, MSG_ARK_MISSING_FULLRHS, MSG_ARK_NO_MALLOC,
    MSG_ARK_NULL_F, RESET_INIT, TRY_AGAIN, ZERO,
};
use crate::arkode_ls_impl::ARKLsMem;
use crate::arkode_mri_tables::{
    mriStepCoupling_GetStageMap, mriStepCoupling_GetStageType, MRIStepCoupling_Free,
    MRIStepCoupling_LoadTable, MRIStepCoupling_Space, MRIStepCoupling_Write, MRISTEP_EXPLICIT,
    MRISTEP_IMEX, MRISTEP_IMPLICIT, MRISTEP_MERK, MRISTEP_SR, ARKODE_MRI_NONE,
    MRISTEP_DEFAULT_EXPL_1, MRISTEP_DEFAULT_EXPL_2, MRISTEP_DEFAULT_EXPL_2_AD,
    MRISTEP_DEFAULT_EXPL_3, MRISTEP_DEFAULT_EXPL_3_AD, MRISTEP_DEFAULT_EXPL_4,
    MRISTEP_DEFAULT_EXPL_4_AD, MRISTEP_DEFAULT_EXPL_5_AD, MRISTEP_DEFAULT_IMEX_SD_1,
    MRISTEP_DEFAULT_IMEX_SD_2, MRISTEP_DEFAULT_IMEX_SD_2_AD, MRISTEP_DEFAULT_IMEX_SD_3,
    MRISTEP_DEFAULT_IMEX_SD_3_AD, MRISTEP_DEFAULT_IMEX_SD_4, MRISTEP_DEFAULT_IMEX_SD_4_AD,
    MRISTEP_DEFAULT_IMPL_SD_1, MRISTEP_DEFAULT_IMPL_SD_2, MRISTEP_DEFAULT_IMPL_SD_3,
    MRISTEP_DEFAULT_IMPL_SD_4,
};
use crate::arkode_mristep_impl::{
    ARKodeMRIStepMem, MRIStepInnerEvolveFn, MRIStepInnerFullRhsFn,
    MRIStepInnerGetAccumulatedError, MRIStepInnerResetAccumulatedError, MRIStepInnerResetFn,
    MRIStepInnerSetRTol, MRIStepInnerStepper, MRISTAGE_DIRK_FAST, MRISTAGE_DIRK_NOFAST,
    MRISTAGE_ERK_FAST, MRISTAGE_ERK_NOFAST, MRISTAGE_STIFF_ACC, MSG_MRISTEP_NO_MEM,
};
use crate::nvector_serial::{NVector, N_VConst, N_VLinearCombination, N_VLinearSum, N_VScale,
    N_VWrmsNorm};
use crate::sundials_adaptcontroller::{
    SUNAdaptController_GetType, SUN_ADAPTCONTROLLER_H, SUN_ADAPTCONTROLLER_MRI_H_TOL,
    SUN_ADAPTCONTROLLER_NONE,
};
use crate::sundials_errors::SUN_SUCCESS;
use crate::sundials_math::{SUNMAX, SUNRabs};
use crate::sundials_stepper::{
    SUNStepper, SUNStepper_Evolve, SUNStepper_FullRhs, SUNStepper_Reset, SUNStepper_SetForcing,
    SUNStepper_SetStopTime, SUNFullRhsMode, SUN_FULLRHS_END, SUN_FULLRHS_OTHER, SUN_FULLRHS_START,
};
use crate::sundials_types::{UserData, SUN_UNIT_ROUNDOFF};
use crate::sunnonlinsol_newton::SUNNonlinSol_Newton;

const ONE: f64 = 1.0;
const TWO: f64 = 2.0;

/*===============================================================
  Shortcut routines to unpack step_mem structure from ark_mem
  (take semantics: caller must put the box back)
  ===============================================================*/
pub(crate) fn mriStep_AccessStepMem(
    ark_mem: &mut ARKodeMem,
    fname: &str,
) -> Option<Box<ARKodeMRIStepMem>> {
    let taken = ark_mem.step_mem.take();
    match taken {
        Some(b) => match b.downcast::<ARKodeMRIStepMem>() {
            Ok(sm) => Some(sm),
            Err(other) => {
                ark_mem.step_mem = Some(other);
                arkProcessError(
                    Some(ark_mem),
                    ARK_MEM_NULL,
                    line!(),
                    fname,
                    file!(),
                    MSG_MRISTEP_NO_MEM,
                );
                None
            }
        },
        None => {
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_NULL,
                line!(),
                fname,
                file!(),
                MSG_MRISTEP_NO_MEM,
            );
            None
        }
    }
}

/*===============================================================
  MRIStep Exported functions -- Required
  ===============================================================*/

pub fn MRIStepCreate(
    fse: Option<ARKRhsFn>,
    fsi: Option<ARKRhsFn>,
    t0: f64,
    y0: &NVector,
    stepper: MRIStepInnerStepper,
    sunctx: &crate::sundials_context::SUNContext,
) -> Option<Box<ARKodeMem>> {
    /* Check that at least one of fse, fsi is supplied and is to be used */
    if fse.is_none() && fsi.is_none() {
        arkProcessError(
            None,
            ARK_ILL_INPUT,
            line!(),
            "MRIStepCreate",
            file!(),
            MSG_ARK_NULL_F,
        );
        return None;
    }

    /* C also checks y0 == NULL, stepper == NULL and sunctx == NULL
    (inexpressible here) */

    /* Create ark_mem structure and set default values */
    let mut ark_mem = arkCreate(sunctx);

    /* Allocate ARKodeMRIStepMem structure, and initialize to zero */
    let mut step_mem = Box::new(ARKodeMRIStepMem::default());

    /* Attach step_mem structure and function pointers to ark_mem */
    ark_mem.step_attachlinsol = Some(mriStep_AttachLinsol);
    ark_mem.step_disablelsetup = Some(mriStep_DisableLSetup);
    ark_mem.step_getlinmem = Some(mriStep_GetLmem);
    ark_mem.step_getimplicitrhs = Some(mriStep_GetImplicitRHS);
    ark_mem.step_getgammas = Some(mriStep_GetGammas);
    ark_mem.step_setjcur = Some(mriStep_SetJcur);
    ark_mem.step_init = Some(mriStep_Init);
    ark_mem.step_fullrhs = Some(mriStep_FullRHS);
    ark_mem.step = Some(mriStep_TakeStepMRIGARK);
    ark_mem.step_setuserdata = Some(crate::arkode_mristep_io::mriStep_SetUserData);
    ark_mem.step_printallstats = Some(crate::arkode_mristep_io::mriStep_PrintAllStats);
    ark_mem.step_writeparameters = Some(crate::arkode_mristep_io::mriStep_WriteParameters);
    ark_mem.step_setusecompensatedsums = None;
    /* C: ark_mem->step_resize = mriStep_Resize (deferred with ARKodeResize) */
    ark_mem.step_reset = Some(mriStep_Reset);
    ark_mem.step_free = Some(mriStep_Free);
    ark_mem.step_printmem = Some(mriStep_PrintMem);
    ark_mem.step_setdefaults = Some(crate::arkode_mristep_io::mriStep_SetDefaults);
    ark_mem.step_computestate = Some(mriStep_ComputeState);
    ark_mem.step_setoptions = Some(crate::arkode_mristep_io::mriStep_SetOptions);
    ark_mem.step_setorder = Some(crate::arkode_mristep_io::mriStep_SetOrder);
    ark_mem.step_setnonlinearsolver =
        Some(crate::arkode_mristep_nls::mriStep_SetNonlinearSolver);
    ark_mem.step_setlinear = Some(crate::arkode_mristep_io::mriStep_SetLinear);
    ark_mem.step_setnonlinear = Some(crate::arkode_mristep_io::mriStep_SetNonlinear);
    ark_mem.step_setnlsrhsfn = Some(crate::arkode_mristep_nls::mriStep_SetNlsRhsFn);
    ark_mem.step_setdeduceimplicitrhs =
        Some(crate::arkode_mristep_io::mriStep_SetDeduceImplicitRhs);
    ark_mem.step_setnonlincrdown = Some(crate::arkode_mristep_io::mriStep_SetNonlinCRDown);
    ark_mem.step_setnonlinrdiv = Some(crate::arkode_mristep_io::mriStep_SetNonlinRDiv);
    ark_mem.step_setdeltagammamax = Some(crate::arkode_mristep_io::mriStep_SetDeltaGammaMax);
    ark_mem.step_setlsetupfrequency = Some(crate::arkode_mristep_io::mriStep_SetLSetupFrequency);
    ark_mem.step_setpredictormethod = Some(crate::arkode_mristep_io::mriStep_SetPredictorMethod);
    ark_mem.step_setmaxnonliniters = Some(crate::arkode_mristep_io::mriStep_SetMaxNonlinIters);
    ark_mem.step_setnonlinconvcoef = Some(crate::arkode_mristep_io::mriStep_SetNonlinConvCoef);
    ark_mem.step_setstagepredictfn = Some(crate::arkode_mristep_io::mriStep_SetStagePredictFn);
    ark_mem.step_getnumrhsevals = Some(crate::arkode_mristep_io::mriStep_GetNumRhsEvals);
    ark_mem.step_getnumlinsolvsetups =
        Some(crate::arkode_mristep_io::mriStep_GetNumLinSolvSetups);
    ark_mem.step_getcurrentgamma = Some(crate::arkode_mristep_io::mriStep_GetCurrentGamma);
    ark_mem.step_setadaptcontroller =
        Some(crate::arkode_mristep_io::mriStep_SetAdaptController);
    ark_mem.step_getestlocalerrors = Some(crate::arkode_mristep_io::mriStep_GetEstLocalErrors);
    ark_mem.step_getnonlinearsystemdata =
        Some(crate::arkode_mristep_nls::mriStep_GetNonlinearSystemData);
    ark_mem.step_getnumnonlinsolviters =
        Some(crate::arkode_mristep_io::mriStep_GetNumNonlinSolvIters);
    ark_mem.step_getnumnonlinsolvconvfails =
        Some(crate::arkode_mristep_io::mriStep_GetNumNonlinSolvConvFails);
    ark_mem.step_getnonlinsolvstats = Some(crate::arkode_mristep_io::mriStep_GetNonlinSolvStats);
    ark_mem.step_setforcing = Some(mriStep_SetInnerForcing);
    ark_mem.step_getstageindex = Some(crate::arkode_mristep_io::mriStep_GetStageIndex);
    /* C declares mriStep_ComputeH0 but 7.7.0 never installs step_H0 */
    ark_mem.step_supports_adaptive = true;
    ark_mem.step_supports_implicit = true;
    ark_mem.step_mem = Some(step_mem);

    /* Set default values for optional inputs */
    let retval = crate::arkode_mristep_io::mriStep_SetDefaults(&mut ark_mem);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(&ark_mem),
            retval,
            line!(),
            "MRIStepCreate",
            file!(),
            "Error setting default solver options",
        );
        return None;
    }

    /* re-take step_mem for direct field initialization */
    step_mem = mriStep_AccessStepMem(&mut ark_mem, "MRIStepCreate")?;

    /* Allocate the general MRI stepper vectors using y0 as a template */
    /* NOTE: Fse, Fsi, inner_forcing, sdata, zpred and zcor will be allocated
       later on (based on the MRI method) */

    /* Copy the slow RHS functions into stepper memory */
    step_mem.fse = fse;
    step_mem.fsi = fsi;
    step_mem.fse_is_current = false;
    step_mem.fsi_is_current = false;

    /* Set implicit/explicit problem based on function pointers */
    step_mem.explicit_rhs = fse.is_some();
    step_mem.implicit_rhs = fsi.is_some();

    /* Update the ARKODE workspace requirements */
    ark_mem.liw += 49; /* fcn/data ptr, int, long int, sunindextype, sunbooleantype */
    ark_mem.lrw += 14;

    /* Create a default Newton NLS object (just in case; will be deleted if
       the user attaches a nonlinear solver) */
    step_mem.NLS = None;
    step_mem.ownNLS = false;

    if step_mem.implicit_rhs {
        let NLS = SUNNonlinSol_Newton(y0, sunctx);
        ark_mem.step_mem = Some(step_mem);
        let retval = crate::arkode_mristep_nls::mriStep_SetNonlinearSolver(&mut ark_mem, NLS);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(&ark_mem),
                ARK_MEM_FAIL,
                line!(),
                "MRIStepCreate",
                file!(),
                "Error attaching default Newton solver",
            );
            return None;
        }
        step_mem = mriStep_AccessStepMem(&mut ark_mem, "MRIStepCreate")?;
        step_mem.ownNLS = true;
    }

    /* Set the linear solver addresses to NULL (we check != NULL later) */
    step_mem.linit = None;
    step_mem.lsetup = None;
    step_mem.lsolve = None;
    step_mem.lfree = None;

    /* Initialize error norm  */
    step_mem.eRNrm = ONE;

    /* Initialize all the counters */
    step_mem.nfse = 0;
    step_mem.nfsi = 0;
    step_mem.nsetups = 0;
    step_mem.nstlp = 0;
    step_mem.nls_iters = 0;
    step_mem.nls_fails = 0;
    step_mem.inner_fails = 0;

    /* Initialize fused op work space with sufficient storage for at least
       filling the full RHS on an ImEx problem -- must be allocated here as the
       full RHS is called before mriStep_Init when nesting MRI methods */
    step_mem.nfusedopvecs = 3;
    step_mem.cvals = vec![ZERO; step_mem.nfusedopvecs as usize];
    ark_mem.lrw += step_mem.nfusedopvecs as i64;
    /* (Xvecs assembled at call sites; keep the C liw accounting) */
    ark_mem.liw += step_mem.nfusedopvecs as i64;

    /* Initialize adaptivity parameters */
    step_mem.inner_rtol_factor = ONE;
    step_mem.inner_dsm = ONE;
    step_mem.inner_rtol_factor_new = ONE;

    /* Initialize pre and post inner evolve functions */
    step_mem.pre_inner_evolve = None;
    step_mem.post_inner_evolve = None;

    /* Initialize external polynomial forcing data */
    step_mem.expforcing = false;
    step_mem.impforcing = false;
    step_mem.forcing = Vec::new();
    step_mem.nforcing = 0;

    /* Attach the inner stepper memory (moved before arkInit so the
    step ops that arkInit may re-enter see a complete step_mem) */
    step_mem.stepper = stepper;

    /* Check for required stepper functions */
    let retval = mriStepInnerStepper_HasRequiredOps(&step_mem.stepper);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(&ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "MRIStepCreate",
            file!(),
            "A required inner stepper function is NULL",
        );
        return None;
    }

    ark_mem.step_mem = Some(step_mem);

    /* Initialize main ARKODE infrastructure (allocates vectors) */
    let retval = arkInit(&mut ark_mem, t0, y0, FIRST_INIT);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(&ark_mem),
            retval,
            line!(),
            "MRIStepCreate",
            file!(),
            "Unable to initialize main ARKODE infrastructure",
        );
        return None;
    }

    /* return ARKODE memory */
    Some(ark_mem)
}

/*---------------------------------------------------------------
  MRIStepReInit:

  This routine re-initializes the MRIStep module to solve a new
  problem of the same size as was previously solved (all counter
  values are set to 0).

  NOTE: the inner stepper needs to be reinitialized before
  calling this function.
  ---------------------------------------------------------------*/
pub fn MRIStepReInit(
    ark_mem: &mut ARKodeMem,
    fse: Option<ARKRhsFn>,
    fsi: Option<ARKRhsFn>,
    t0: f64,
    y0: &NVector,
) -> i32 {
    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "MRIStepReInit") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* Check if ark_mem was allocated */
    if !ark_mem.MallocDone {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            crate::arkode_impl::ARK_NO_MALLOC,
            line!(),
            "MRIStepReInit",
            file!(),
            MSG_ARK_NO_MALLOC,
        );
        return crate::arkode_impl::ARK_NO_MALLOC;
    }

    /* Check that at least one of fse, fsi is supplied and is to be used */
    if fse.is_none() && fsi.is_none() {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "MRIStepReInit",
            file!(),
            MSG_ARK_NULL_F,
        );
        return ARK_ILL_INPUT;
    }

    /* C also checks y0 == NULL (inexpressible) */

    /* Set implicit/explicit problem based on function pointers */
    step_mem.explicit_rhs = fse.is_some();
    step_mem.implicit_rhs = fsi.is_some();

    /* Create a default Newton NLS object (just in case; will be deleted if
       the user attaches a nonlinear solver) */
    if step_mem.implicit_rhs && step_mem.NLS.is_none() {
        let sunctx = crate::sundials_context::SUNContext_Create();
        let NLS = SUNNonlinSol_Newton(y0, &sunctx);
        ark_mem.step_mem = Some(step_mem);
        let retval = crate::arkode_mristep_nls::mriStep_SetNonlinearSolver(ark_mem, NLS);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_FAIL,
                line!(),
                "MRIStepReInit",
                file!(),
                "Error attaching default Newton solver",
            );
            return ARK_MEM_FAIL;
        }
        step_mem = match mriStep_AccessStepMem(ark_mem, "MRIStepReInit") {
            Some(sm) => sm,
            None => return ARK_MEM_NULL,
        };
        step_mem.ownNLS = true;
    }

    ark_mem.step_mem = Some(step_mem);

    /* ReInitialize main ARKODE infrastructure */
    let retval = arkInit(ark_mem, t0, y0, FIRST_INIT);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            retval,
            line!(),
            "MRIStepReInit",
            file!(),
            "Unable to reinitialize main ARKODE infrastructure",
        );
        return retval;
    }

    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "MRIStepReInit") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* Copy the input parameters into ARKODE state */
    step_mem.fse = fse;
    step_mem.fsi = fsi;
    step_mem.fse_is_current = false;
    step_mem.fsi_is_current = false;

    /* Initialize all the counters */
    step_mem.nfse = 0;
    step_mem.nfsi = 0;
    step_mem.nsetups = 0;
    step_mem.nstlp = 0;
    step_mem.nls_iters = 0;
    step_mem.nls_fails = 0;
    step_mem.inner_fails = 0;

    ark_mem.step_mem = Some(step_mem);

    if let Some(lmem) = ark_mem.lmem.as_mut() {
        crate::arkode_ls::arkLsInitializeCounters(lmem);
    }

    ARK_SUCCESS
}

/*===============================================================
  Interface routines supplied to ARKODE
  ===============================================================*/

/* (mriStep_Resize deferred with ARKodeResize) */

/*---------------------------------------------------------------
  mriStep_Reset:

  This routine resets the MRIStep module state to solve the same
  problem from the given time with the input state (all counter
  values are retained).  It is called after the main ARKODE
  infrastructure is reset.
  ---------------------------------------------------------------*/
fn mriStep_Reset(ark_mem: &mut ARKodeMem, tR: f64, yR: &NVector) -> i32 {
    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_Reset") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* Reset the inner integrator with this same state */
    let retval = mriStepInnerStepper_Reset(&mut step_mem.stepper, tR, yR);
    ark_mem.step_mem = Some(step_mem);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            ARK_INNERSTEP_FAIL,
            line!(),
            "mriStep_Reset",
            file!(),
            "Unable to reset the inner stepper",
        );
        return ARK_INNERSTEP_FAIL;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_ComputeState:

  Computes y based on the current prediction and given correction.
  ---------------------------------------------------------------*/
fn mriStep_ComputeState(ark_mem: &mut ARKodeMem, zcor: &NVector, z: &mut NVector) -> i32 {
    let step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_ComputeState") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    N_VLinearSum(ONE, &step_mem.zpred, ONE, zcor, z);

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_Free frees all MRIStep memory.
  ---------------------------------------------------------------*/
fn mriStep_Free(ark_mem: &mut ARKodeMem) {
    let taken = ark_mem.step_mem.take();
    if let Some(b) = taken {
        if let Ok(mut step_mem) = b.downcast::<ARKodeMRIStepMem>() {
            /* free the coupling structure and derived quantities */
            if step_mem.MRIC.is_some() {
                let mut Cliw: i64 = 0;
                let mut Clrw: i64 = 0;
                MRIStepCoupling_Space(step_mem.MRIC.as_deref(), &mut Cliw, &mut Clrw);
                MRIStepCoupling_Free(&mut step_mem.MRIC);
                ark_mem.liw -= Cliw;
                ark_mem.lrw -= Clrw;
                if !step_mem.stagetypes.is_empty() {
                    step_mem.stagetypes = Vec::new();
                    ark_mem.liw -= (step_mem.stages + 1) as i64;
                }
                if !step_mem.stage_map.is_empty() {
                    step_mem.stage_map = Vec::new();
                    ark_mem.liw -= step_mem.stages as i64;
                }
                if !step_mem.Ae_row.is_empty() {
                    step_mem.Ae_row = Vec::new();
                    ark_mem.lrw -= step_mem.stages as i64;
                }
                if !step_mem.Ai_row.is_empty() {
                    step_mem.Ai_row = Vec::new();
                    ark_mem.lrw -= step_mem.stages as i64;
                }
            }

            /* free the nonlinear solver memory (if applicable) */
            step_mem.NLS = None;
            step_mem.ownNLS = false;

            /* free the linear solver memory */
            if let Some(lfree) = step_mem.lfree {
                lfree(ark_mem);
            }

            /* free the sdata, zpred and zcor vectors */
            arkFreeVec(ark_mem, &mut step_mem.sdata);
            arkFreeVec(ark_mem, &mut step_mem.zpred);
            arkFreeVec(ark_mem, &mut step_mem.zcor);

            /* free the RHS vectors */
            if !step_mem.Fse.is_empty() {
                let (lrw1, liw1) = (ark_mem.lrw1, ark_mem.liw1);
                let ARKodeMem { lrw, liw, .. } = ark_mem;
                arkFreeVecArray(step_mem.nstages_allocated, &mut step_mem.Fse, lrw1, lrw, liw1, liw);
            }

            if !step_mem.Fsi.is_empty() {
                let (lrw1, liw1) = (ark_mem.lrw1, ark_mem.liw1);
                let ARKodeMem { lrw, liw, .. } = ark_mem;
                arkFreeVecArray(step_mem.nstages_allocated, &mut step_mem.Fsi, lrw1, lrw, liw1, liw);
            }

            /* free the reusable arrays for fused vector interface */
            if !step_mem.cvals.is_empty() {
                step_mem.cvals = Vec::new();
                ark_mem.lrw -= step_mem.nfusedopvecs as i64;
                ark_mem.liw -= step_mem.nfusedopvecs as i64;
            }
            step_mem.nfusedopvecs = 0;

            /* the time stepper module itself is dropped here */
        }
    }
}

/*---------------------------------------------------------------
  mriStep_PrintMem:

  This routine outputs the memory from the MRIStep structure to
  a specified file pointer (useful when debugging).
  ---------------------------------------------------------------*/
fn mriStep_PrintMem(ark_mem: &mut ARKodeMem, outfile: &mut dyn std::io::Write) {
    use crate::sundials_utils::fmt_g;

    let step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_PrintMem") {
        Some(sm) => sm,
        None => return,
    };

    /* output integer quantities */
    let _ = writeln!(outfile, "MRIStep: q = {}", step_mem.q);
    let _ = writeln!(outfile, "MRIStep: p = {}", step_mem.p);
    let _ = writeln!(outfile, "MRIStep: istage = {}", step_mem.istage);
    let _ = writeln!(outfile, "MRIStep: cur_stage = {}", step_mem.cur_stage);
    let _ = writeln!(outfile, "MRIStep: stages = {}", step_mem.stages);
    let _ = writeln!(outfile, "MRIStep: maxcor = {}", step_mem.maxcor);
    let _ = writeln!(outfile, "MRIStep: msbp = {}", step_mem.msbp);
    let _ = writeln!(outfile, "MRIStep: predictor = {}", step_mem.predictor);
    let _ = writeln!(outfile, "MRIStep: convfail = {}", step_mem.convfail);
    let _ = write!(outfile, "MRIStep: stagetypes =");
    for i in 0..=step_mem.stages as usize {
        let _ = write!(outfile, " {}", step_mem.stagetypes[i]);
    }
    let _ = writeln!(outfile);

    /* output long integer quantities */
    let _ = writeln!(outfile, "MRIStep: nfse = {}", step_mem.nfse);
    let _ = writeln!(outfile, "MRIStep: nfsi = {}", step_mem.nfsi);
    let _ = writeln!(outfile, "MRIStep: nsetups = {}", step_mem.nsetups);
    let _ = writeln!(outfile, "MRIStep: nstlp = {}", step_mem.nstlp);
    let _ = writeln!(outfile, "MRIStep: nls_iters = {}", step_mem.nls_iters);
    let _ = writeln!(outfile, "MRIStep: nls_fails = {}", step_mem.nls_fails);
    let _ = writeln!(outfile, "MRIStep: inner_fails = {}", step_mem.inner_fails);

    /* output boolean quantities */
    let _ = writeln!(outfile, "MRIStep: user_linear = {}", step_mem.linear as i32);
    let _ = writeln!(
        outfile,
        "MRIStep: user_linear_timedep = {}",
        step_mem.linear_timedep as i32
    );
    let _ = writeln!(
        outfile,
        "MRIStep: user_explicit = {}",
        step_mem.explicit_rhs as i32
    );
    let _ = writeln!(
        outfile,
        "MRIStep: user_implicit = {}",
        step_mem.implicit_rhs as i32
    );
    let _ = writeln!(outfile, "MRIStep: jcur = {}", step_mem.jcur as i32);
    let _ = writeln!(outfile, "MRIStep: ownNLS = {}", step_mem.ownNLS as i32);

    /* output sunrealtype quantities */
    let _ = writeln!(outfile, "MRIStep: Coupling structure:");
    if let Some(MRIC) = step_mem.MRIC.as_deref() {
        MRIStepCoupling_Write(MRIC, outfile);
    }

    let _ = writeln!(outfile, "MRIStep: gamma = {}", fmt_g(step_mem.gamma, 0, 6));
    let _ = writeln!(outfile, "MRIStep: gammap = {}", fmt_g(step_mem.gammap, 0, 6));
    let _ = writeln!(outfile, "MRIStep: gamrat = {}", fmt_g(step_mem.gamrat, 0, 6));
    let _ = writeln!(outfile, "MRIStep: crate = {}", fmt_g(step_mem.crate_, 0, 6));
    let _ = writeln!(outfile, "MRIStep: delp = {}", fmt_g(step_mem.delp, 0, 6));
    let _ = writeln!(outfile, "MRIStep: eRNrm = {}", fmt_g(step_mem.eRNrm, 0, 6));
    let _ = writeln!(outfile, "MRIStep: nlscoef = {}", fmt_g(step_mem.nlscoef, 0, 6));
    let _ = writeln!(outfile, "MRIStep: crdown = {}", fmt_g(step_mem.crdown, 0, 6));
    let _ = writeln!(outfile, "MRIStep: rdiv = {}", fmt_g(step_mem.rdiv, 0, 6));
    let _ = writeln!(outfile, "MRIStep: dgmax = {}", fmt_g(step_mem.dgmax, 0, 6));
    let _ = write!(outfile, "MRIStep: Ae_row =");
    for i in 0..step_mem.nstages_active as usize {
        let _ = write!(outfile, " {}", fmt_g(step_mem.Ae_row[i], 0, 6));
    }
    let _ = writeln!(outfile);
    let _ = write!(outfile, "MRIStep: Ai_row =");
    for i in 0..step_mem.nstages_active as usize {
        let _ = write!(outfile, " {}", fmt_g(step_mem.Ai_row[i], 0, 6));
    }
    let _ = writeln!(outfile);

    /* print the inner stepper memory */
    mriStepInnerStepper_PrintMem(&step_mem.stepper, outfile);

    ark_mem.step_mem = Some(step_mem);
}

/*---------------------------------------------------------------
  mriStep_AttachLinsol:

  This routine attaches the various set of system linear solver
  interface routines, data structure, and solver type to the
  MRIStep module.
  ---------------------------------------------------------------*/
fn mriStep_AttachLinsol(
    ark_mem: &mut ARKodeMem,
    linit: Option<ARKLinsolInitFn>,
    lsetup: Option<ARKLinsolSetupFn>,
    lsolve: Option<ARKLinsolSolveFn>,
    lfree: Option<ARKLinsolFreeFn>,
    _lsolve_type: crate::sundials_linearsolver::SUNLinearSolver_Type,
    lmem: Box<ARKLsMem>,
) -> i32 {
    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_AttachLinsol") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* free any existing system solver */
    if let Some(old_lfree) = step_mem.lfree {
        old_lfree(ark_mem);
    }

    /* Attach the provided routines, data structure and solve type */
    step_mem.linit = linit;
    step_mem.lsetup = lsetup;
    step_mem.lsolve = lsolve;
    step_mem.lfree = lfree;
    ark_mem.lmem = Some(lmem);

    /* Reset all linear solver counters */
    step_mem.nsetups = 0;
    step_mem.nstlp = 0;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_DisableLSetup:

  This routine NULLifies the lsetup function pointer in the
  MRIStep module.
  ---------------------------------------------------------------*/
fn mriStep_DisableLSetup(ark_mem: &mut ARKodeMem) {
    if let Some(b) = ark_mem.step_mem.as_mut() {
        if let Some(step_mem) = b.downcast_mut::<ARKodeMRIStepMem>() {
            step_mem.lsetup = None;
        }
    }
}

/*---------------------------------------------------------------
  mriStep_GetLmem:

  This routine returns the system linear solver interface memory
  (take semantics; put-back writes ark_mem.lmem).
  ---------------------------------------------------------------*/
fn mriStep_GetLmem(ark_mem: &mut ARKodeMem) -> Option<Box<ARKLsMem>> {
    ark_mem.lmem.take()
}

/*---------------------------------------------------------------
  mriStep_GetImplicitRHS:

  This routine returns the implicit RHS function pointer, fi.
  ---------------------------------------------------------------*/
fn mriStep_GetImplicitRHS(ark_mem: &mut ARKodeMem) -> Option<ARKRhsFn> {
    if let Some(b) = ark_mem.step_mem.as_mut() {
        if let Some(step_mem) = b.downcast_mut::<ARKodeMRIStepMem>() {
            if step_mem.implicit_rhs {
                return step_mem.fsi;
            }
        }
    }
    None
}

/*---------------------------------------------------------------
  mriStep_GetGammas:

  This routine fills the current value of gamma, and states
  whether the gamma ratio fails the dgmax criteria.
  ---------------------------------------------------------------*/
fn mriStep_GetGammas(
    ark_mem: &mut ARKodeMem,
    gamma: &mut f64,
    gamrat: &mut f64,
    jcur: &mut bool,
    dgamma_fail: &mut bool,
) -> i32 {
    let step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_GetGammas") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* set outputs */
    *gamma = step_mem.gamma;
    *gamrat = step_mem.gamrat;
    *jcur = step_mem.jcur;
    *dgamma_fail = SUNRabs(*gamrat - ONE) >= step_mem.dgmax;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/* Rust-only companion op: write-back for the jcur pointer C hands
out via step_getgammas (Addendum C.1) */
fn mriStep_SetJcur(ark_mem: &mut ARKodeMem, jcur: bool) {
    if let Some(b) = ark_mem.step_mem.as_mut() {
        if let Some(step_mem) = b.downcast_mut::<ARKodeMRIStepMem>() {
            step_mem.jcur = jcur;
        }
    }
}

/*---------------------------------------------------------------
  mriStep_Init:

  This routine is called just prior to performing internal time
  steps (after all user "set" routines have been called) from
  within arkInitialSetup.

  With initialization type RESET_INIT, this routine does nothing.
  ---------------------------------------------------------------*/
fn mriStep_Init(ark_mem: &mut ARKodeMem, init_type: i32) -> i32 {
    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_Init") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };
    let retval = mriStep_Init_inner(ark_mem, &mut step_mem, init_type);
    ark_mem.step_mem = Some(step_mem);
    retval
}

fn mriStep_Init_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeMRIStepMem>,
    init_type: i32,
) -> i32 {
    use crate::arkode_impl::{ARK_ACCUMERROR_NONE, ARK_INTERP_HERMITE};

    /* immediately return if reset */
    if init_type == RESET_INIT {
        return ARK_SUCCESS;
    }

    /* initializations/checks for (re-)initialization call */
    if init_type == FIRST_INIT {
        /* enforce use of arkEwtSmallReal if using a fixed step size for
           an explicit method, an internal error weight function, and not
           performing accumulated temporal error estimation */
        let mut reset_efun = true;
        if step_mem.implicit_rhs {
            reset_efun = false;
        }
        if !ark_mem.fixedstep {
            reset_efun = false;
        }
        if ark_mem.user_efun {
            reset_efun = false;
        }
        if ark_mem.AccumErrorType != ARK_ACCUMERROR_NONE {
            reset_efun = false;
        }
        if reset_efun {
            ark_mem.user_efun = false;
            ark_mem.efun = Some(mri_ewt_small_real);
        }

        /* Create coupling structure (if not already set) */
        let retval = mriStep_SetCoupling(ark_mem, step_mem);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "mriStep_Init",
                file!(),
                "Could not create coupling table",
            );
            return ARK_ILL_INPUT;
        }

        /* Check that coupling structure is OK */
        let retval = mriStep_CheckCoupling(ark_mem, step_mem);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "mriStep_Init",
                file!(),
                "Error in coupling table",
            );
            return ARK_ILL_INPUT;
        }

        /* Attach correct TakeStep routine for this coupling table */
        match step_mem.MRIC.as_ref().unwrap().type_ {
            MRISTEP_EXPLICIT | MRISTEP_IMPLICIT | MRISTEP_IMEX => {
                ark_mem.step = Some(mriStep_TakeStepMRIGARK);
            }
            MRISTEP_MERK => {
                ark_mem.step = Some(mriStep_TakeStepMERK);
            }
            MRISTEP_SR => {
                ark_mem.step = Some(mriStep_TakeStepMRISR);
            }
            _ => {
                arkProcessError(
                    Some(ark_mem),
                    ARK_ILL_INPUT,
                    line!(),
                    "mriStep_Init",
                    file!(),
                    "Unknown method type",
                );
                return ARK_ILL_INPUT;
            }
        }

        /* Request arkode ensure that ycur==yn upon entry to TakeStep function */
        ark_mem.ensure_ycur = true;

        /* Retrieve/store method and embedding orders now that tables are
        finalized */
        {
            let MRIC = step_mem.MRIC.as_ref().unwrap();
            step_mem.stages = MRIC.stages;
            step_mem.q = MRIC.q;
            step_mem.p = MRIC.p;
            let hadapt_mem = ark_mem.hadapt_mem.as_mut().unwrap();
            hadapt_mem.q = step_mem.q;
            hadapt_mem.p = step_mem.p;
        }

        /* Ensure that if adaptivity or error accumulation is enabled, then
           method includes embedding coefficients */
        if (!ark_mem.fixedstep || ark_mem.AccumErrorType != ARK_ACCUMERROR_NONE)
            && step_mem.p <= 0
        {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "mriStep_Init",
                file!(),
                "Temporal error estimation cannot be performed without embedding coefficients",
            );
            return ARK_ILL_INPUT;
        }

        /* allocate/fill derived quantities from MRIC structure */

        /* stage map */
        if !step_mem.stage_map.is_empty() {
            step_mem.stage_map = Vec::new();
            ark_mem.liw -= step_mem.stages as i64;
        }
        step_mem.stage_map = vec![0; step_mem.MRIC.as_ref().unwrap().stages as usize];
        ark_mem.liw += step_mem.MRIC.as_ref().unwrap().stages as i64;
        {
            let ARKodeMRIStepMem { MRIC, stage_map, nstages_active, .. } = &mut **step_mem;
            let retval =
                mriStepCoupling_GetStageMap(MRIC.as_deref().unwrap(), stage_map, nstages_active);
            if retval != ARK_SUCCESS {
                arkProcessError(
                    Some(ark_mem),
                    ARK_ILL_INPUT,
                    line!(),
                    "mriStep_Init",
                    file!(),
                    "Error in coupling table",
                );
                return ARK_ILL_INPUT;
            }
        }

        /* stage types */
        if !step_mem.stagetypes.is_empty() {
            step_mem.stagetypes = Vec::new();
            ark_mem.liw -= step_mem.stages as i64;
        }
        step_mem.stagetypes = vec![0; (step_mem.MRIC.as_ref().unwrap().stages + 1) as usize];
        ark_mem.liw += (step_mem.MRIC.as_ref().unwrap().stages + 1) as i64;
        for j in 0..=step_mem.MRIC.as_ref().unwrap().stages {
            step_mem.stagetypes[j as usize] =
                mriStepCoupling_GetStageType(step_mem.MRIC.as_deref().unwrap(), j);
        }

        /* explicit RK coefficient row */
        if !step_mem.Ae_row.is_empty() {
            step_mem.Ae_row = Vec::new();
            ark_mem.lrw -= step_mem.stages as i64;
        }
        step_mem.Ae_row = vec![ZERO; step_mem.MRIC.as_ref().unwrap().stages as usize];
        ark_mem.lrw += step_mem.MRIC.as_ref().unwrap().stages as i64;

        /* implicit RK coefficient row */
        if !step_mem.Ai_row.is_empty() {
            step_mem.Ai_row = Vec::new();
            ark_mem.lrw -= step_mem.stages as i64;
        }
        step_mem.Ai_row = vec![ZERO; step_mem.MRIC.as_ref().unwrap().stages as usize];
        ark_mem.lrw += step_mem.MRIC.as_ref().unwrap().stages as i64;

        /* Allocate reusable arrays for fused vector operations */
        let fused_workspace_size = std::cmp::max(
            3,
            2 * step_mem.MRIC.as_ref().unwrap().stages + 2 + step_mem.nforcing,
        );

        if step_mem.nfusedopvecs < fused_workspace_size {
            if !step_mem.cvals.is_empty() {
                step_mem.cvals = Vec::new();
                ark_mem.lrw -= step_mem.nfusedopvecs as i64;
                ark_mem.liw -= step_mem.nfusedopvecs as i64;
            }
            step_mem.nfusedopvecs = 0;

            step_mem.cvals = vec![ZERO; fused_workspace_size as usize];
            step_mem.nfusedopvecs = fused_workspace_size;
            ark_mem.lrw += fused_workspace_size as i64;
            ark_mem.liw += fused_workspace_size as i64;
        }

        /* Retrieve/store method and embedding orders now that tables are
        finalized */
        {
            let MRIC = step_mem.MRIC.as_ref().unwrap();
            step_mem.stages = MRIC.stages;
            step_mem.q = MRIC.q;
            step_mem.p = MRIC.p;
        }

        /* If an MRISR method is applied to a non-ImEx problem, we "unify"
           the Fse and Fsi vectors to point at the same memory */
        step_mem.unify_Fs = false;
        if step_mem.MRIC.as_ref().unwrap().type_ == MRISTEP_SR
            && ((step_mem.explicit_rhs && !step_mem.implicit_rhs)
                || (!step_mem.explicit_rhs && step_mem.implicit_rhs))
        {
            step_mem.unify_Fs = true;
        }

        /* Allocate MRI RHS vector memory, update storage requirements */
        /*   Allocate Fse[0] ... Fse[nstages_active - 1] and           */
        /*   Fsi[0] ... Fsi[nstages_active - 1] if needed              */
        if step_mem.nstages_allocated < step_mem.nstages_active {
            let tmpl_len = ark_mem.ewt.data.len();
            let (lrw1, liw1) = (ark_mem.lrw1, ark_mem.liw1);
            if step_mem.nstages_allocated > 0 {
                if step_mem.explicit_rhs {
                    let ARKodeMem { lrw, liw, .. } = ark_mem;
                    arkFreeVecArray(
                        step_mem.nstages_allocated,
                        &mut step_mem.Fse,
                        lrw1,
                        lrw,
                        liw1,
                        liw,
                    );
                }
                if step_mem.implicit_rhs && !step_mem.unify_Fs {
                    let ARKodeMem { lrw, liw, .. } = ark_mem;
                    arkFreeVecArray(
                        step_mem.nstages_allocated,
                        &mut step_mem.Fsi,
                        lrw1,
                        lrw,
                        liw1,
                        liw,
                    );
                }
            }
            if step_mem.explicit_rhs && !step_mem.unify_Fs {
                let ARKodeMem { lrw, liw, .. } = ark_mem;
                if !arkAllocVecArray(
                    step_mem.nstages_active,
                    tmpl_len,
                    &mut step_mem.Fse,
                    lrw1,
                    lrw,
                    liw1,
                    liw,
                ) {
                    return ARK_MEM_FAIL;
                }
            }
            if step_mem.implicit_rhs && !step_mem.unify_Fs {
                let ARKodeMem { lrw, liw, .. } = ark_mem;
                if !arkAllocVecArray(
                    step_mem.nstages_active,
                    tmpl_len,
                    &mut step_mem.Fsi,
                    lrw1,
                    lrw,
                    liw1,
                    liw,
                ) {
                    return ARK_MEM_FAIL;
                }
            }
            if step_mem.unify_Fs {
                /* unified storage lives in Fse (see impl header note) */
                let ARKodeMem { lrw, liw, .. } = ark_mem;
                if !arkAllocVecArray(
                    step_mem.nstages_active,
                    tmpl_len,
                    &mut step_mem.Fse,
                    lrw1,
                    lrw,
                    liw1,
                    liw,
                ) {
                    return ARK_MEM_FAIL;
                }
            }

            step_mem.nstages_allocated = step_mem.nstages_active;
        }

        /* if any slow stage is implicit, allocate sdata, zpred, zcor vectors;
           if all stages explicit, free default NLS object, and detach all
           linear solver routines. */
        if step_mem.implicit_rhs {
            let tmpl_len = ark_mem.ewt.data.len();
            if !arkAllocVec(ark_mem, tmpl_len, &mut step_mem.sdata) {
                return ARK_MEM_FAIL;
            }
            if !arkAllocVec(ark_mem, tmpl_len, &mut step_mem.zpred) {
                return ARK_MEM_FAIL;
            }
            if !arkAllocVec(ark_mem, tmpl_len, &mut step_mem.zcor) {
                return ARK_MEM_FAIL;
            }
        } else {
            if step_mem.NLS.is_some() && step_mem.ownNLS {
                step_mem.NLS = None;
                step_mem.ownNLS = false;
            }
            step_mem.linit = None;
            step_mem.lsetup = None;
            step_mem.lsolve = None;
            step_mem.lfree = None;
            ark_mem.lmem = None;
        }

        /* Allocate inner stepper data */
        let tmpl_len = ark_mem.ewt.data.len();
        let retval = mriStepInnerStepper_AllocVecs(
            &mut step_mem.stepper,
            step_mem.MRIC.as_ref().unwrap().nmat,
            tmpl_len,
        );
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "mriStep_Init",
                file!(),
                "Error allocating inner stepper memory",
            );
            return ARK_MEM_FAIL;
        }

        /* Override the interpolant degree (if needed), used in arkInitialSetup */
        if step_mem.q > 1 && ark_mem.interp_degree > (step_mem.q - 1) {
            /* Limit max degree to at most one less than the method global order */
            ark_mem.interp_degree = step_mem.q - 1;
        } else if step_mem.q == 1 && ark_mem.interp_degree > 1 {
            /* Allow for linear interpolant with first order methods to ensure
               solution values are returned at the time interval end points */
            ark_mem.interp_degree = 1;
        }

        /* Higher-order predictors require interpolation */
        if ark_mem.interp_type == ARK_INTERP_NONE && step_mem.predictor != 0 {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "mriStep_Init",
                file!(),
                "Non-trival predictors require an interpolation module",
            );
            return ARK_ILL_INPUT;
        }
        let _ = ARK_INTERP_HERMITE; /* (referenced by C via UpdateF0 docs) */
    }

    /* Call linit (if it exists) */
    if let Some(linit) = step_mem.linit {
        let retval = with_step_mem_installed(ark_mem, step_mem, linit);
        if retval != 0 {
            arkProcessError(
                Some(ark_mem),
                ARK_LINIT_FAIL,
                line!(),
                "mriStep_Init",
                file!(),
                "The linear solver's init routine failed.",
            );
            return ARK_LINIT_FAIL;
        }
    }

    /* Initialize the nonlinear solver object (if it exists) */
    if step_mem.NLS.is_some() {
        let retval = crate::arkode_mristep_nls::mriStep_NlsInit(ark_mem, step_mem);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_NLS_INIT_FAIL,
                line!(),
                "mriStep_Init",
                file!(),
                "Unable to initialize SUNNonlinearSolver object",
            );
            return ARK_NLS_INIT_FAIL;
        }
    }

    /* get timestep adaptivity type */
    let adapt_type = match ark_mem
        .hadapt_mem
        .as_ref()
        .and_then(|h| h.hcontroller.as_ref())
    {
        Some(hc) => SUNAdaptController_GetType(hc),
        None => SUN_ADAPTCONTROLLER_NONE,
    };

    if ark_mem.fixedstep {
        /* Fixed step sizes: user must supply the initial step size */
        if ark_mem.hin == ZERO {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "mriStep_Init",
                file!(),
                "Timestep adaptivity disabled, but missing user-defined fixed stepsize",
            );
            return ARK_ILL_INPUT;
        }
    } else {
        /* ensure that a compatible adaptivity controller is provided */
        if adapt_type != SUN_ADAPTCONTROLLER_MRI_H_TOL && adapt_type != SUN_ADAPTCONTROLLER_H {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "mriStep_Init",
                file!(),
                "SUNAdaptController type is unsupported by MRIStep",
            );
            return ARK_ILL_INPUT;
        }

        /* Controller provides adaptivity (at least at the slow time scale):
           - verify that the MRI method includes an embedding */
        if step_mem.MRIC.as_ref().unwrap().p <= 0 {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "mriStep_Init",
                file!(),
                "Timestep adaptivity enabled, but non-embedded MRI table specified",
            );
            return ARK_ILL_INPUT;
        }
    }

    /* Perform additional setup for (H,tol) controller */
    if adapt_type == SUN_ADAPTCONTROLLER_MRI_H_TOL {
        /* Verify that adaptivity type is supported by inner stepper */
        if !mriStepInnerStepper_SupportsRTolAdaptivity(&step_mem.stepper) {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "mriStep_Init",
                file!(),
                "MRI H-TOL SUNAdaptController provided, but unsupported by inner stepper",
            );
            return ARK_ILL_INPUT;
        }

        /* initialize fast stepper to use the same relative tolerance as MRIStep */
        step_mem.inner_rtol_factor = ONE;
    }

    ARK_SUCCESS
}

/* local 3-arg shim installing arkEwtSetSmallReal as the internal efun
(the ARKEwtFn type takes (y, w, e_data); same pattern as the erk/lsrk
steppers) */
fn mri_ewt_small_real(_y: &NVector, w: &mut NVector, _e_data: &mut UserData) -> i32 {
    N_VConst(crate::sundials_types::SUN_SMALL_REAL, w);
    ARK_SUCCESS
}

/* helper: temporarily re-install step_mem into ark_mem around an
   ARKLS linit/lsetup/lsolve call (those re-enter the step_getgammas /
   step_getimplicitrhs / step_setjcur ops, which need step_mem in
   place — same pattern as arkode_arkstep_nls.rs). */
pub(crate) fn with_step_mem_installed<R>(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeMRIStepMem>,
    f: impl FnOnce(&mut ARKodeMem) -> R,
) -> R {
    let owned = std::mem::replace(step_mem, Box::new(ARKodeMRIStepMem::default()));
    ark_mem.step_mem = Some(owned);
    let r = f(ark_mem);
    let owned = ark_mem
        .step_mem
        .take()
        .unwrap()
        .downcast::<ARKodeMRIStepMem>()
        .unwrap();
    *step_mem = owned;
    r
}

/*------------------------------------------------------------------------------
  mriStep_ComputeH0:

  This utility routine computes the initial slow step size for MRI methods.
  ----------------------------------------------------------------------------*/
pub fn mriStep_ComputeH0(ark_mem: &mut ARKodeMem, tout: f64, hin: &mut f64) -> i32 {
    /*   tempv1 = fs(t0, y0) */
    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_ComputeH0") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };
    let mut tempv1 = std::mem::take(&mut ark_mem.tempv1);
    let yn = std::mem::take(&mut ark_mem.yn);
    let retval = mriStep_SlowRHS_inner(
        ark_mem,
        &mut step_mem,
        ark_mem.tn,
        &yn,
        &mut tempv1,
        ARK_FULLRHS_START,
    );
    ark_mem.yn = yn;
    if retval != ARK_SUCCESS {
        ark_mem.tempv1 = tempv1;
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_RHSFUNC_FAIL,
            line!(),
            "mriStep_ComputeH0",
            file!(),
            "error calling slow RHS function(s)",
        );
        return ARK_RHSFUNC_FAIL;
    }
    let retval = mriStep_Hin(ark_mem, ark_mem.tn, tout, &tempv1, hin);
    ark_mem.tempv1 = tempv1;
    ark_mem.step_mem = Some(step_mem);
    if retval != ARK_SUCCESS {
        return arkHandleFailure(ark_mem, retval);
    }

    ARK_SUCCESS
}

/*------------------------------------------------------------------------------
  mriStep_FullRHS:

  This is just a wrapper to call the user-supplied RHS functions,
  f(t,y) = fse(t,y) + fsi(t,y)  + ff(t,y).
  ----------------------------------------------------------------------------*/
fn mriStep_FullRHS(ark_mem: &mut ARKodeMem, t: f64, y: &NVector, f: &mut NVector, mode: i32) -> i32 {
    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_FullRHS") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };
    let retval = mriStep_FullRHS_inner(ark_mem, &mut step_mem, t, y, f, mode);
    ark_mem.step_mem = Some(step_mem);
    retval
}

fn mriStep_FullRHS_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeMRIStepMem>,
    t: f64,
    y: &NVector,
    f: &mut NVector,
    mode: i32,
) -> i32 {
    /* ensure that inner stepper provides fullrhs function */
    if step_mem.stepper.ops.fullrhs.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_RHSFUNC_FAIL,
            line!(),
            "mriStep_FullRHS",
            file!(),
            MSG_ARK_MISSING_FULLRHS,
        );
        return ARK_RHSFUNC_FAIL;
    }

    /* perform RHS functions contingent on 'mode' argument */
    match mode {
        ARK_FULLRHS_START | ARK_FULLRHS_END => {
            /* update the internal storage for Fse[0] and Fsi[0] */
            let retval = mriStep_UpdateF0(ark_mem, step_mem, t, y, mode);
            if retval != 0 {
                arkProcessError(
                    Some(ark_mem),
                    ARK_RHSFUNC_FAIL,
                    line!(),
                    "mriStep_FullRHS",
                    file!(),
                    &rhsfunc_failed_msg(t),
                );
                return ARK_RHSFUNC_FAIL;
            }

            /* evaluate fast component */
            let retval =
                mriStepInnerStepper_FullRhs(&mut step_mem.stepper, t, y, f, ARK_FULLRHS_OTHER);
            if retval != ARK_SUCCESS {
                arkProcessError(
                    Some(ark_mem),
                    ARK_RHSFUNC_FAIL,
                    line!(),
                    "mriStep_FullRHS",
                    file!(),
                    &rhsfunc_failed_msg(t),
                );
                return ARK_RHSFUNC_FAIL;
            }

            /* combine RHS vectors into output */
            if step_mem.explicit_rhs && step_mem.implicit_rhs {
                /* ImEx: z == X[0] with c0 == 1 -> in-place accumulate form */
                mri_lincomb_accumulate(&[ONE, ONE], &[&step_mem.Fse[0], &step_mem.Fsi[0]], f);
            } else if step_mem.implicit_rhs {
                /* implicit (C N_VLinearSum(1,Fsi[0],1,f,f) aliases z with y) */
                let fsi0 = if step_mem.unify_Fs {
                    &step_mem.Fse[0]
                } else {
                    &step_mem.Fsi[0]
                };
                f.linear_sum_with(ONE, ONE, fsi0);
            } else {
                /* explicit */
                f.linear_sum_with(ONE, ONE, &step_mem.Fse[0]);
            }
        }

        ARK_FULLRHS_OTHER => {
            /* compute the fast component (force new RHS computation) */
            let retval =
                mriStepInnerStepper_FullRhs(&mut step_mem.stepper, t, y, f, ARK_FULLRHS_OTHER);
            if retval != ARK_SUCCESS {
                arkProcessError(
                    Some(ark_mem),
                    ARK_RHSFUNC_FAIL,
                    line!(),
                    "mriStep_FullRHS",
                    file!(),
                    &rhsfunc_failed_msg(t),
                );
                return ARK_RHSFUNC_FAIL;
            }

            /* call the user-supplied pre-RHS function (if supplied) */
            if let Some(pre_rhs) = ark_mem.PreRhsFn {
                let retval = pre_rhs(t, y, &mut ark_mem.user_data);
                if retval != 0 {
                    return ARK_PRERHSFN_FAIL;
                }
            }

            /* compute the implicit component and store in sdata */
            if step_mem.implicit_rhs {
                let fsi = step_mem.fsi.unwrap();
                let retval = fsi(t, y, &mut step_mem.sdata, &mut ark_mem.user_data);
                step_mem.nfsi += 1;
                if retval != 0 {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_RHSFUNC_FAIL,
                        line!(),
                        "mriStep_FullRHS",
                        file!(),
                        &rhsfunc_failed_msg(t),
                    );
                    return ARK_RHSFUNC_FAIL;
                }
                /* accumulate f += sdata (z == X[0] c0 == 1 form of the C
                fused op) */
                f.linear_sum_with(ONE, ONE, &step_mem.sdata);
            }

            /* compute the explicit component and store in ark_tempv2 */
            if step_mem.explicit_rhs {
                let fse = step_mem.fse.unwrap();
                let retval = fse(t, y, &mut ark_mem.tempv2, &mut ark_mem.user_data);
                step_mem.nfse += 1;
                if retval != 0 {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_RHSFUNC_FAIL,
                        line!(),
                        "mriStep_FullRHS",
                        file!(),
                        &rhsfunc_failed_msg(t),
                    );
                    return ARK_RHSFUNC_FAIL;
                }
                f.linear_sum_with(ONE, ONE, &ark_mem.tempv2);
            }

            /* Add external forcing components to linear combination */
            if step_mem.expforcing || step_mem.impforcing {
                let vals = mriStep_ApplyForcing_coeffs(step_mem, t, ONE);
                mri_accumulate_forcing(step_mem, &vals, f);
            }
        }

        _ => {
            /* return with RHS failure if unknown mode is passed */
            arkProcessError(
                Some(ark_mem),
                ARK_RHSFUNC_FAIL,
                line!(),
                "mriStep_FullRHS",
                file!(),
                "Unknown full RHS mode",
            );
            return ARK_RHSFUNC_FAIL;
        }
    }

    ARK_SUCCESS
}

/* MSG_ARK_RHSFUNC_FAILED with the time formatted like the other steppers */
fn rhsfunc_failed_msg(t: f64) -> String {
    format!(
        "At t = {}, the right-hand side routine failed in an unrecoverable manner.",
        crate::sundials_utils::fmt_g(t, 0, 15)
    )
}

/// z += sum_k c[k] * x[k] — the z == X[0], c[0] == 1 branch of the
/// C N_VLinearCombination kernel.
fn mri_lincomb_accumulate(cvals: &[f64], xvecs: &[&NVector], z: &mut NVector) {
    for (k, val) in cvals.iter().enumerate() {
        for e in 0..z.data.len() {
            z.data[e] += val * xvecs[k].data[e];
        }
    }
}

/// z += sum_k vals[k] * forcing[k] — forcing tail of a C fused op with
/// z == X[0], c0 == 1.
fn mri_accumulate_forcing(step_mem: &ARKodeMRIStepMem, vals: &[f64], f: &mut NVector) {
    for (k, val) in vals.iter().enumerate() {
        for e in 0..f.data.len() {
            f.data[e] += val * step_mem.forcing[k].data[e];
        }
    }
}

/*------------------------------------------------------------------------------
  mriStep_UpdateF0:

  This routine is called by mriStep_FullRHS to update the internal storage for
  Fse[0] and Fsi[0], incorporating forcing from a slower time scale as
  necessary.
  ----------------------------------------------------------------------------*/
#[allow(clippy::nonminimal_bool)] /* C's condition shape kept verbatim */
pub(crate) fn mriStep_UpdateF0(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeMRIStepMem>,
    t: f64,
    y: &NVector,
    mode: i32,
) -> i32 {
    /* perform RHS functions contingent on 'mode' argument */
    match mode {
        ARK_FULLRHS_START => {
            /* update the RHS components */

            /* call the user-supplied pre-RHS function (if supplied) */
            if let Some(pre_rhs) = ark_mem.PreRhsFn {
                if (!step_mem.fse_is_current || !ark_mem.fn_is_current)
                    || (!step_mem.fsi_is_current || !ark_mem.fn_is_current)
                {
                    let retval = pre_rhs(t, y, &mut ark_mem.user_data);
                    if retval != 0 {
                        return ARK_PRERHSFN_FAIL;
                    }
                }
            }

            /*   implicit component */
            if step_mem.implicit_rhs {
                /* if either ARKODE or MRIStep consider Fsi[0] stale, then recompute */
                if !step_mem.fsi_is_current || !ark_mem.fn_is_current {
                    let fsi = step_mem.fsi.unwrap();
                    {
                        let ARKodeMRIStepMem { Fse, Fsi, unify_Fs, .. } = &mut **step_mem;
                        let fsi_arr = if *unify_Fs { Fse } else { Fsi };
                        let retval = fsi(t, y, &mut fsi_arr[0], &mut ark_mem.user_data);
                        if retval != 0 {
                            step_mem.nfsi += 1;
                            arkProcessError(
                                Some(ark_mem),
                                ARK_RHSFUNC_FAIL,
                                line!(),
                                "mriStep_UpdateF0",
                                file!(),
                                &rhsfunc_failed_msg(t),
                            );
                            return ARK_RHSFUNC_FAIL;
                        }
                    }
                    step_mem.nfsi += 1;
                    step_mem.fsi_is_current = true;

                    /* Add external forcing, if applicable */
                    if step_mem.impforcing {
                        let vals = mriStep_ApplyForcing_coeffs(step_mem, t, ONE);
                        let ARKodeMRIStepMem { Fse, Fsi, unify_Fs, forcing, .. } =
                            &mut **step_mem;
                        let fsi_arr = if *unify_Fs { Fse } else { Fsi };
                        for (k, val) in vals.iter().enumerate() {
                            for e in 0..fsi_arr[0].data.len() {
                                fsi_arr[0].data[e] += val * forcing[k].data[e];
                            }
                        }
                    }
                }
            }

            /*   explicit component */
            if step_mem.explicit_rhs {
                /* if either ARKODE or MRIStep consider Fse[0] stale, then recompute */
                if !step_mem.fse_is_current || !ark_mem.fn_is_current {
                    let fse = step_mem.fse.unwrap();
                    let retval = fse(t, y, &mut step_mem.Fse[0], &mut ark_mem.user_data);
                    step_mem.nfse += 1;
                    if retval != 0 {
                        arkProcessError(
                            Some(ark_mem),
                            ARK_RHSFUNC_FAIL,
                            line!(),
                            "mriStep_UpdateF0",
                            file!(),
                            &rhsfunc_failed_msg(t),
                        );
                        return ARK_RHSFUNC_FAIL;
                    }
                    step_mem.fse_is_current = true;

                    /* Add external forcing, if applicable */
                    if step_mem.expforcing {
                        let vals = mriStep_ApplyForcing_coeffs(step_mem, t, ONE);
                        let ARKodeMRIStepMem { Fse, forcing, .. } = &mut **step_mem;
                        for (k, val) in vals.iter().enumerate() {
                            for e in 0..Fse[0].data.len() {
                                Fse[0].data[e] += val * forcing[k].data[e];
                            }
                        }
                    }
                }
            }
        }

        ARK_FULLRHS_END => {
            /* compute the full RHS */
            if !ark_mem.fn_is_current {
                /* call the user-supplied pre-RHS function (if supplied) */
                if let Some(pre_rhs) = ark_mem.PreRhsFn {
                    let retval = pre_rhs(t, y, &mut ark_mem.user_data);
                    if retval != 0 {
                        return ARK_PRERHSFN_FAIL;
                    }
                }

                /* compute the implicit component */
                if step_mem.implicit_rhs {
                    let fsi = step_mem.fsi.unwrap();
                    {
                        let ARKodeMRIStepMem { Fse, Fsi, unify_Fs, .. } = &mut **step_mem;
                        let fsi_arr = if *unify_Fs { Fse } else { Fsi };
                        let retval = fsi(t, y, &mut fsi_arr[0], &mut ark_mem.user_data);
                        if retval != 0 {
                            step_mem.nfsi += 1;
                            arkProcessError(
                                Some(ark_mem),
                                ARK_RHSFUNC_FAIL,
                                line!(),
                                "mriStep_UpdateF0",
                                file!(),
                                &rhsfunc_failed_msg(t),
                            );
                            return ARK_RHSFUNC_FAIL;
                        }
                    }
                    step_mem.nfsi += 1;
                    step_mem.fsi_is_current = true;

                    /* Add external forcing, as appropriate */
                    if step_mem.impforcing {
                        let vals = mriStep_ApplyForcing_coeffs(step_mem, t, ONE);
                        let ARKodeMRIStepMem { Fse, Fsi, unify_Fs, forcing, .. } =
                            &mut **step_mem;
                        let fsi_arr = if *unify_Fs { Fse } else { Fsi };
                        for (k, val) in vals.iter().enumerate() {
                            for e in 0..fsi_arr[0].data.len() {
                                fsi_arr[0].data[e] += val * forcing[k].data[e];
                            }
                        }
                    }
                }

                /* compute the explicit component */
                if step_mem.explicit_rhs {
                    let fse = step_mem.fse.unwrap();
                    let retval = fse(t, y, &mut step_mem.Fse[0], &mut ark_mem.user_data);
                    step_mem.nfse += 1;
                    if retval != 0 {
                        arkProcessError(
                            Some(ark_mem),
                            ARK_RHSFUNC_FAIL,
                            line!(),
                            "mriStep_UpdateF0",
                            file!(),
                            &rhsfunc_failed_msg(t),
                        );
                        return ARK_RHSFUNC_FAIL;
                    }
                    step_mem.fse_is_current = true;

                    /* Add external forcing, as appropriate */
                    if step_mem.expforcing {
                        let vals = mriStep_ApplyForcing_coeffs(step_mem, t, ONE);
                        let ARKodeMRIStepMem { Fse, forcing, .. } = &mut **step_mem;
                        for (k, val) in vals.iter().enumerate() {
                            for e in 0..Fse[0].data.len() {
                                Fse[0].data[e] += val * forcing[k].data[e];
                            }
                        }
                    }
                }
            }
        }

        _ => {
            /* return with RHS failure if unknown mode is requested */
            arkProcessError(
                Some(ark_mem),
                ARK_RHSFUNC_FAIL,
                line!(),
                "mriStep_UpdateF0",
                file!(),
                "Unknown full RHS mode",
            );
            return ARK_RHSFUNC_FAIL;
        }
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_TakeStepMRIGARK:

  This routine serves the primary purpose of the MRIStep module:
  it performs a single MRI step (with embedding, if possible).
  See the C source for the dsmPtr/nflagPtr conventions.
  ---------------------------------------------------------------*/
fn mriStep_TakeStepMRIGARK(ark_mem: &mut ARKodeMem, dsmPtr: &mut f64, nflagPtr: &mut i32) -> i32 {
    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_TakeStepMRIGARK") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };
    let retval = mriStep_TakeStepMRIGARK_inner(ark_mem, &mut step_mem, dsmPtr, nflagPtr);
    ark_mem.step_mem = Some(step_mem);
    retval
}

fn mriStep_TakeStepMRIGARK_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeMRIStepMem>,
    dsmPtr: &mut f64,
    nflagPtr: &mut i32,
) -> i32 {
    use crate::arkode_impl::ARK_ACCUMERROR_NONE;

    /* initialize algebraic solver convergence flag to success;
       error estimate to zero */
    *nflagPtr = ARK_SUCCESS;
    *dsmPtr = ZERO;

    /* determine whether embedding stage is needed */
    let do_embedding =
        !ark_mem.fixedstep || (ark_mem.AccumErrorType != ARK_ACCUMERROR_NONE);

    /* initialize the current stage index */
    step_mem.istage = 0;
    step_mem.cur_stage = 0;

    /* if MRI adaptivity is enabled: reset fast accumulated error,
       and send appropriate control parameter to the fast integrator */
    let adapt_type = match ark_mem
        .hadapt_mem
        .as_ref()
        .and_then(|h| h.hcontroller.as_ref())
    {
        Some(hc) => SUNAdaptController_GetType(hc),
        None => SUN_ADAPTCONTROLLER_NONE,
    };
    let mut need_inner_dsm = false;
    if adapt_type == SUN_ADAPTCONTROLLER_MRI_H_TOL {
        need_inner_dsm = true;
        step_mem.inner_dsm = ZERO;
        let retval = mriStepInnerStepper_ResetAccumulatedError(&mut step_mem.stepper);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_INNERSTEP_FAIL,
                line!(),
                "mriStep_TakeStepMRIGARK",
                file!(),
                "Unable to reset the inner stepper error estimate",
            );
            return ARK_INNERSTEP_FAIL;
        }
        let retval = mriStepInnerStepper_SetRTol(
            &mut step_mem.stepper,
            step_mem.inner_rtol_factor * ark_mem.reltol,
        );
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_INNERSTEP_FAIL,
                line!(),
                "mriStep_TakeStepMRIGARK",
                file!(),
                "Unable to set the inner stepper tolerance",
            );
            return ARK_INNERSTEP_FAIL;
        }
    }

    /* for adaptive computations, reset the inner integrator to the beginning
    of this step */
    if !ark_mem.fixedstep {
        let retval =
            mriStepInnerStepper_Reset(&mut step_mem.stepper, ark_mem.tcur, &ark_mem.ycur);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_INNERSTEP_FAIL,
                line!(),
                "mriStep_TakeStepMRIGARK",
                file!(),
                "Unable to reset the inner stepper",
            );
            return ARK_INNERSTEP_FAIL;
        }
    }

    /* call nonlinear solver setup if it exists (the crate's inlined
    Newton/fixed-point solvers have no setup phase; C calls
    SUNNonlinSolSetup here) */

    /* Evaluate the slow RHS functions if needed. */
    let nested_mri = step_mem.expforcing || step_mem.impforcing;
    if ark_mem.fn_.data.is_empty() || nested_mri {
        let ycur = std::mem::take(&mut ark_mem.ycur);
        let retval = mriStep_UpdateF0(ark_mem, step_mem, ark_mem.tcur, &ycur, ARK_FULLRHS_START);
        ark_mem.ycur = ycur;
        if retval != 0 {
            return ARK_RHSFUNC_FAIL;
        }

        /* For a nested MRI configuration we might still need fn to create a
        predictor but it should be fn only for the current nesting level */
        if !ark_mem.fn_.data.is_empty() && nested_mri && step_mem.implicit_rhs {
            if step_mem.implicit_rhs && step_mem.explicit_rhs {
                N_VLinearSum(ONE, &step_mem.Fsi[0], ONE, &step_mem.Fse[0], &mut ark_mem.fn_);
            } else {
                let fsi0 = if step_mem.unify_Fs {
                    &step_mem.Fse[0]
                } else {
                    &step_mem.Fsi[0]
                };
                N_VScale(ONE, fsi0, &mut ark_mem.fn_);
            }
        }
    } else if !ark_mem.fn_.data.is_empty() && !ark_mem.fn_is_current {
        let ycur = std::mem::take(&mut ark_mem.ycur);
        let mut fn_ = std::mem::take(&mut ark_mem.fn_);
        let retval = mriStep_FullRHS_inner(
            ark_mem,
            step_mem,
            ark_mem.tcur,
            &ycur,
            &mut fn_,
            ARK_FULLRHS_START,
        );
        ark_mem.ycur = ycur;
        ark_mem.fn_ = fn_;
        if retval != 0 {
            return ARK_RHSFUNC_FAIL;
        }
    }
    ark_mem.fn_is_current = true;

    /* The first stage is the previous time-step solution, so its RHS
       is the [already-computed] slow RHS from the start of the step */

    /* Loop over remaining internal stages */
    let mut is: i32 = 1;
    while is < step_mem.stages - 1 {
        /* Set relevant stage times (including desired stage time for implicit
        solves) and stage index */
        let t0 = ark_mem.tn + step_mem.MRIC.as_ref().unwrap().c[(is - 1) as usize] * ark_mem.h;
        let tf = ark_mem.tn + step_mem.MRIC.as_ref().unwrap().c[is as usize] * ark_mem.h;
        ark_mem.tcur = tf;
        step_mem.istage = is;
        step_mem.cur_stage = is;

        /* Determine current stage type, and call corresponding routine */
        let retval = match step_mem.stagetypes[is as usize] {
            MRISTAGE_ERK_FAST => {
                let retval = mriStep_ComputeInnerForcing(ark_mem, step_mem, is, t0, tf);
                if retval != ARK_SUCCESS {
                    return retval;
                }
                let retval = mriStep_StageERKFast(ark_mem, step_mem, t0, tf, need_inner_dsm);
                if retval != ARK_SUCCESS {
                    *nflagPtr = CONV_FAIL;
                }
                retval
            }
            MRISTAGE_ERK_NOFAST => mriStep_StageERKNoFast(ark_mem, step_mem, is),
            MRISTAGE_DIRK_NOFAST => mriStep_StageDIRKNoFast(ark_mem, step_mem, is, nflagPtr),
            MRISTAGE_DIRK_FAST => mriStep_StageDIRKFast(ark_mem, step_mem, is, nflagPtr),
            _ => ARK_SUCCESS, /* MRISTAGE_STIFF_ACC */
        };
        if retval != ARK_SUCCESS {
            return retval;
        }

        /* apply user-supplied stage postprocessing function (if supplied) */
        if let (Some(post), true) = (
            ark_mem.PostProcessStageFn,
            step_mem.stagetypes[is as usize] != MRISTAGE_STIFF_ACC,
        ) {
            let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
            let retval = post(*tcur, ycur, user_data);
            if retval != 0 {
                return ARK_POSTPROCESS_STAGE_FAIL;
            }
        }

        /* conditionally reset the inner integrator with the modified stage
        solution */
        if step_mem.stagetypes[is as usize] != MRISTAGE_STIFF_ACC
            && (step_mem.stagetypes[is as usize] != MRISTAGE_ERK_FAST
                || ark_mem.PostProcessStageFn.is_some())
        {
            let retval = mriStepInnerStepper_Reset(&mut step_mem.stepper, tf, &ark_mem.ycur);
            if retval != ARK_SUCCESS {
                arkProcessError(
                    Some(ark_mem),
                    ARK_INNERSTEP_FAIL,
                    line!(),
                    "mriStep_TakeStepMRIGARK",
                    file!(),
                    "Unable to reset the inner stepper",
                );
                return ARK_INNERSTEP_FAIL;
            }
        }

        /* Compute updated slow RHS, except:
           1. if the stage is excluded from stage_map
           2. if the next stage has "STIFF_ACC" type, and temporal estimation
              is disabled */
        let mut calc_fslow = true;
        if step_mem.stage_map[is as usize] == -1 {
            calc_fslow = false;
        }
        if !do_embedding && step_mem.stagetypes[(is + 1) as usize] == MRISTAGE_STIFF_ACC {
            calc_fslow = false;
        }
        if calc_fslow {
            /* call the user-supplied pre-RHS function (if supplied) */
            if let Some(pre_rhs) = ark_mem.PreRhsFn {
                if step_mem.explicit_rhs
                    || (step_mem.implicit_rhs
                        && (!step_mem.deduce_rhs
                            || step_mem.stagetypes[is as usize] != MRISTAGE_DIRK_NOFAST))
                {
                    let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
                    let retval = pre_rhs(*tcur, ycur, user_data);
                    if retval != 0 {
                        return ARK_PRERHSFN_FAIL;
                    }
                }
            }

            /* store implicit slow rhs  */
            if step_mem.implicit_rhs {
                if !step_mem.deduce_rhs
                    || step_mem.stagetypes[is as usize] != MRISTAGE_DIRK_NOFAST
                {
                    let fsi = step_mem.fsi.unwrap();
                    let smap = step_mem.stage_map[is as usize] as usize;
                    let retval = {
                        let ARKodeMRIStepMem { Fse, Fsi, unify_Fs, .. } = &mut **step_mem;
                        let fsi_arr = if *unify_Fs { Fse } else { Fsi };
                        let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
                        fsi(*tcur, ycur, &mut fsi_arr[smap], user_data)
                    };
                    step_mem.nfsi += 1;

                    if retval < 0 {
                        return ARK_RHSFUNC_FAIL;
                    }
                    if retval > 0 {
                        return ARK_UNREC_RHSFUNC_ERR;
                    }

                    /* Add external forcing to Fsi, if applicable */
                    if step_mem.impforcing {
                        let vals = mriStep_ApplyForcing_coeffs(step_mem, tf, ONE);
                        let ARKodeMRIStepMem { Fse, Fsi, unify_Fs, forcing, .. } =
                            &mut **step_mem;
                        let fsi_arr = if *unify_Fs { Fse } else { Fsi };
                        for (k, val) in vals.iter().enumerate() {
                            for e in 0..fsi_arr[smap].data.len() {
                                fsi_arr[smap].data[e] += val * forcing[k].data[e];
                            }
                        }
                    }
                } else {
                    let smap = step_mem.stage_map[is as usize] as usize;
                    let gamma = step_mem.gamma;
                    let ARKodeMRIStepMem { Fse, Fsi, unify_Fs, zcor, sdata, .. } =
                        &mut **step_mem;
                    let fsi_arr = if *unify_Fs { Fse } else { Fsi };
                    N_VLinearSum(ONE / gamma, zcor, -ONE / gamma, sdata, &mut fsi_arr[smap]);
                }
            }

            /* store explicit slow rhs */
            if step_mem.explicit_rhs {
                let fse = step_mem.fse.unwrap();
                let smap = step_mem.stage_map[is as usize] as usize;
                let retval = {
                    let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
                    fse(*tcur, ycur, &mut step_mem.Fse[smap], user_data)
                };
                step_mem.nfse += 1;

                if retval < 0 {
                    return ARK_RHSFUNC_FAIL;
                }
                if retval > 0 {
                    return ARK_UNREC_RHSFUNC_ERR;
                }

                /* Add external forcing to Fse, if applicable */
                if step_mem.expforcing {
                    let vals = mriStep_ApplyForcing_coeffs(step_mem, tf, ONE);
                    let ARKodeMRIStepMem { Fse, forcing, .. } = &mut **step_mem;
                    for (k, val) in vals.iter().enumerate() {
                        for e in 0..Fse[smap].data.len() {
                            Fse[smap].data[e] += val * forcing[k].data[e];
                        }
                    }
                }
            }
        } /* compute slow RHS */

        is += 1;
    } /* loop over stages */

    /* perform embedded stage (if needed) */
    if do_embedding {
        let is = step_mem.stages;
        step_mem.istage = is;
        step_mem.cur_stage = is;

        /* Temporarily swap ark_mem->ycur and ark_mem->tempv4 pointers, copying
           data so that both hold the current ark_mem->ycur value.  This
           ensures that during this embedding "stage":
             - ark_mem->ycur will be the correct initial condition for the
               final stage.
             - ark_mem->tempv4 will hold the embedded solution vector. */
        {
            let ARKodeMem { ycur, tempv4, .. } = ark_mem;
            N_VScale(ONE, ycur, tempv4);
        }
        std::mem::swap(&mut ark_mem.ycur, &mut ark_mem.tempv4);

        /* Set relevant stage times (including desired stage time for implicit
        solves) */
        let t0 = ark_mem.tn + step_mem.MRIC.as_ref().unwrap().c[(is - 2) as usize] * ark_mem.h;
        let tf = ark_mem.tn + ark_mem.h;
        ark_mem.tcur = tf;

        /* Determine embedding stage type, and call corresponding routine */
        let retval = match step_mem.stagetypes[is as usize] {
            MRISTAGE_ERK_FAST => {
                let retval = mriStep_ComputeInnerForcing(ark_mem, step_mem, is, t0, tf);
                if retval != ARK_SUCCESS {
                    return retval;
                }
                let retval = mriStep_StageERKFast(ark_mem, step_mem, t0, tf, false);
                if retval != ARK_SUCCESS {
                    *nflagPtr = CONV_FAIL;
                }
                retval
            }
            MRISTAGE_ERK_NOFAST => mriStep_StageERKNoFast(ark_mem, step_mem, is),
            MRISTAGE_DIRK_NOFAST => mriStep_StageDIRKNoFast(ark_mem, step_mem, is, nflagPtr),
            MRISTAGE_DIRK_FAST => mriStep_StageDIRKFast(ark_mem, step_mem, is, nflagPtr),
            _ => ARK_SUCCESS,
        };
        if retval != ARK_SUCCESS {
            return retval;
        }

        /* Swap back ark_mem->ycur with ark_mem->tempv4, and reset the inner
        integrator */
        std::mem::swap(&mut ark_mem.ycur, &mut ark_mem.tempv4);
        let retval = mriStepInnerStepper_Reset(&mut step_mem.stepper, t0, &ark_mem.ycur);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_INNERSTEP_FAIL,
                line!(),
                "mriStep_TakeStepMRIGARK",
                file!(),
                "Unable to reset the inner stepper",
            );
            return ARK_INNERSTEP_FAIL;
        }
    }

    /* Compute final stage (for evolved solution), along with error estimate */
    {
        let is = step_mem.stages - 1;
        step_mem.istage = is;
        step_mem.cur_stage = is;

        /* Set relevant stage times (including desired stage time for implicit
        solves) */
        let t0 = ark_mem.tn + step_mem.MRIC.as_ref().unwrap().c[(is - 1) as usize] * ark_mem.h;
        let tf = ark_mem.tn + ark_mem.h;
        ark_mem.tcur = tf;

        /* Determine final stage type, and call corresponding routine */
        let retval = match step_mem.stagetypes[is as usize] {
            MRISTAGE_ERK_FAST => {
                let retval = mriStep_ComputeInnerForcing(ark_mem, step_mem, is, t0, tf);
                if retval != ARK_SUCCESS {
                    return retval;
                }
                let retval = mriStep_StageERKFast(ark_mem, step_mem, t0, tf, need_inner_dsm);
                if retval != ARK_SUCCESS {
                    *nflagPtr = CONV_FAIL;
                }
                retval
            }
            MRISTAGE_ERK_NOFAST => mriStep_StageERKNoFast(ark_mem, step_mem, is),
            MRISTAGE_DIRK_NOFAST => mriStep_StageDIRKNoFast(ark_mem, step_mem, is, nflagPtr),
            MRISTAGE_DIRK_FAST => mriStep_StageDIRKFast(ark_mem, step_mem, is, nflagPtr),
            _ => ARK_SUCCESS,
        };
        if retval != ARK_SUCCESS {
            return retval;
        }

        /* apply user-supplied step postprocessing function (if supplied) */
        if let (Some(post), true) = (
            ark_mem.PostProcessStepFn,
            step_mem.stagetypes[is as usize] != MRISTAGE_STIFF_ACC,
        ) {
            let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
            let retval = post(*tcur, ycur, user_data);
            if retval != 0 {
                return ARK_POSTPROCESS_STEP_FAIL;
            }
        }

        /* conditionally reset the inner integrator with the modified stage
        solution */
        if step_mem.stagetypes[is as usize] != MRISTAGE_STIFF_ACC
            && (step_mem.stagetypes[is as usize] != MRISTAGE_ERK_FAST
                || ark_mem.PostProcessStepFn.is_some())
        {
            let retval = mriStepInnerStepper_Reset(&mut step_mem.stepper, tf, &ark_mem.ycur);
            if retval != ARK_SUCCESS {
                arkProcessError(
                    Some(ark_mem),
                    ARK_INNERSTEP_FAIL,
                    line!(),
                    "mriStep_TakeStepMRIGARK",
                    file!(),
                    "Unable to reset the inner stepper",
                );
                return ARK_INNERSTEP_FAIL;
            }
        }

        /* Compute temporal error estimate via difference between step
           solution and embedding, store in ark_mem->tempv1, and take norm. */
        if do_embedding {
            let ARKodeMem { ycur, tempv4, tempv1, .. } = ark_mem;
            N_VLinearSum(ONE, tempv4, -ONE, ycur, tempv1);
            *dsmPtr = N_VWrmsNorm(&ark_mem.tempv1, &ark_mem.ewt);
        }
    } /* loop over stages */

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_TakeStepMRISR:

  This routine performs a single MRISR step.
  ---------------------------------------------------------------*/
fn mriStep_TakeStepMRISR(ark_mem: &mut ARKodeMem, dsmPtr: &mut f64, nflagPtr: &mut i32) -> i32 {
    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_TakeStepMRISR") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };
    let retval = mriStep_TakeStepMRISR_inner(ark_mem, &mut step_mem, dsmPtr, nflagPtr);
    ark_mem.step_mem = Some(step_mem);
    retval
}

fn mriStep_TakeStepMRISR_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeMRIStepMem>,
    dsmPtr: &mut f64,
    nflagPtr: &mut i32,
) -> i32 {
    use crate::arkode_impl::ARK_ACCUMERROR_NONE;
    let tol: f64 = 100.0 * SUN_UNIT_ROUNDOFF;

    /* initialize algebraic solver convergence flag to success;
       error estimate to zero */
    *nflagPtr = ARK_SUCCESS;
    *dsmPtr = ZERO;

    /* initialize the current stage index */
    step_mem.istage = 0;
    step_mem.cur_stage = 0;

    /* if MRI adaptivity is enabled: reset fast accumulated error,
       and send appropriate control parameter to the fast integrator */
    let adapt_type = match ark_mem
        .hadapt_mem
        .as_ref()
        .and_then(|h| h.hcontroller.as_ref())
    {
        Some(hc) => SUNAdaptController_GetType(hc),
        None => SUN_ADAPTCONTROLLER_NONE,
    };
    let mut need_inner_dsm = false;
    if adapt_type == SUN_ADAPTCONTROLLER_MRI_H_TOL {
        need_inner_dsm = true;
        step_mem.inner_dsm = ZERO;
        let retval = mriStepInnerStepper_ResetAccumulatedError(&mut step_mem.stepper);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_INNERSTEP_FAIL,
                line!(),
                "mriStep_TakeStepMRISR",
                file!(),
                "Unable to reset the inner stepper error estimate",
            );
            return ARK_INNERSTEP_FAIL;
        }
        let retval = mriStepInnerStepper_SetRTol(
            &mut step_mem.stepper,
            step_mem.inner_rtol_factor * ark_mem.reltol,
        );
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_INNERSTEP_FAIL,
                line!(),
                "mriStep_TakeStepMRISR",
                file!(),
                "Unable to set the inner stepper tolerance",
            );
            return ARK_INNERSTEP_FAIL;
        }
    }

    /* for adaptive computations, reset the inner integrator to the beginning
    of this step */
    if !ark_mem.fixedstep {
        let retval =
            mriStepInnerStepper_Reset(&mut step_mem.stepper, ark_mem.tcur, &ark_mem.ycur);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_INNERSTEP_FAIL,
                line!(),
                "mriStep_TakeStepMRISR",
                file!(),
                "Unable to reset the inner stepper",
            );
            return ARK_INNERSTEP_FAIL;
        }
    }

    /* Evaluate the slow RHS functions if needed. */
    let nested_mri = step_mem.expforcing || step_mem.impforcing;
    if ark_mem.fn_.data.is_empty() || nested_mri {
        let ycur = std::mem::take(&mut ark_mem.ycur);
        let retval = mriStep_UpdateF0(ark_mem, step_mem, ark_mem.tcur, &ycur, ARK_FULLRHS_START);
        ark_mem.ycur = ycur;
        if retval != 0 {
            return ARK_RHSFUNC_FAIL;
        }

        /* For a nested MRI configuration we might still need fn to create a
        predictor but it should be fn only for the current nesting level */
        if !ark_mem.fn_.data.is_empty() && nested_mri && step_mem.implicit_rhs {
            if step_mem.implicit_rhs && step_mem.explicit_rhs {
                N_VLinearSum(ONE, &step_mem.Fsi[0], ONE, &step_mem.Fse[0], &mut ark_mem.fn_);
            } else {
                let fsi0 = if step_mem.unify_Fs {
                    &step_mem.Fse[0]
                } else {
                    &step_mem.Fsi[0]
                };
                N_VScale(ONE, fsi0, &mut ark_mem.fn_);
            }
        }
    }
    if !ark_mem.fn_.data.is_empty() && !ark_mem.fn_is_current {
        let ycur = std::mem::take(&mut ark_mem.ycur);
        let mut fn_ = std::mem::take(&mut ark_mem.fn_);
        let retval = mriStep_FullRHS_inner(
            ark_mem,
            step_mem,
            ark_mem.tcur,
            &ycur,
            &mut fn_,
            ARK_FULLRHS_START,
        );
        ark_mem.ycur = ycur;
        ark_mem.fn_ = fn_;
        if retval != 0 {
            return ARK_RHSFUNC_FAIL;
        }
    }
    ark_mem.fn_is_current = true;

    /* combine both RHS into FSE for ImEx problems, since MRISR fast forcing
       function only depends on Omega coefficients  */
    if step_mem.implicit_rhs && step_mem.explicit_rhs {
        let ARKodeMRIStepMem { Fse, Fsi, .. } = &mut **step_mem;
        Fse[0].linear_sum_with(ONE, ONE, &Fsi[0]);
    }

    /* Determine how many stages will be needed */
    let max_stages = if ark_mem.fixedstep && ark_mem.AccumErrorType == ARK_ACCUMERROR_NONE {
        step_mem.stages
    } else {
        step_mem.stages + 1
    };

    /* Loop over stages */
    for stage in 1..max_stages {
        /* Determine if this is an "embedding" or "solution" stage */
        let solution = stage == step_mem.stages - 1;
        let embedding = stage == step_mem.stages;

        /* Set initial condition for this stage (all but first stage) */
        if stage > 1 {
            let ARKodeMem { yn, ycur, .. } = ark_mem;
            N_VScale(ONE, yn, ycur);
        }

        /* Set current stage abscissa */
        let cstage = if embedding {
            ONE
        } else {
            step_mem.MRIC.as_ref().unwrap().c[stage as usize]
        };

        /* Set current stage time and index */
        ark_mem.tcur = ark_mem.tn + cstage * ark_mem.h;
        step_mem.istage = stage;
        step_mem.cur_stage = stage;

        /* Compute forcing function for inner solver */
        let retval =
            mriStep_ComputeInnerForcing(ark_mem, step_mem, stage, ark_mem.tn, ark_mem.tcur);
        if retval != ARK_SUCCESS {
            return retval;
        }

        /* Reset the inner stepper on all but the first stage due to
           "stage-restart" structure */
        if stage > 1 {
            let retval =
                mriStepInnerStepper_Reset(&mut step_mem.stepper, ark_mem.tn, &ark_mem.ycur);
            if retval != ARK_SUCCESS {
                arkProcessError(
                    Some(ark_mem),
                    ARK_INNERSTEP_FAIL,
                    line!(),
                    "mriStep_TakeStepMRISR",
                    file!(),
                    "Unable to reset the inner stepper",
                );
                return ARK_INNERSTEP_FAIL;
            }
        }

        /* Evolve fast IVP for this stage, potentially get inner dsm on
           all non-embedding stages */
        let tf = ark_mem.tcur;
        let retval = mriStep_StageERKFast(
            ark_mem,
            step_mem,
            ark_mem.tn,
            tf,
            need_inner_dsm && !embedding,
        );
        if retval != ARK_SUCCESS {
            *nflagPtr = CONV_FAIL;
            return retval;
        }

        /* perform MRISR slow/implicit correction */
        let mut impl_corr = false;
        if step_mem.implicit_rhs {
            /* determine whether implicit RHS correction will require an
            implicit solve */
            impl_corr = SUNRabs(
                step_mem.MRIC.as_ref().unwrap().G[0][stage as usize][stage as usize],
            ) > tol;

            /* perform implicit solve for correction */
            if impl_corr {
                /* update stage index for prediction and nonlinear solver if
                this is an "embedded" stage */
                if embedding {
                    step_mem.istage = stage - 1;
                }

                /* Call predictor for current stage solution (result placed in
                zpred) */
                let istage = step_mem.istage;
                let mut zpred = std::mem::take(&mut step_mem.zpred);
                let retval = mriStep_Predict(ark_mem, step_mem, istage, &mut zpred);
                step_mem.zpred = zpred;
                if retval != ARK_SUCCESS {
                    return retval;
                }

                /* If a user-supplied predictor routine is provided, call that
                here */
                if let Some(stage_predict) = step_mem.stage_predict {
                    let retval =
                        stage_predict(ark_mem.tcur, &mut step_mem.zpred, &mut ark_mem.user_data);
                    if retval < 0 {
                        return ARK_USER_PREDICT_FAIL;
                    }
                    if retval > 0 {
                        return TRY_AGAIN;
                    }
                }

                /* fill sdata with explicit contributions to correction:
                   sdata = ycur - zpred + h*sum_j G[0][stage][j]*Fsi[j]
                   (dest sdata distinct from all operands -> free fused op) */
                {
                    let mut cv: Vec<f64> = Vec::with_capacity(stage as usize + 2);
                    cv.push(ONE);
                    cv.push(-ONE);
                    for j in 0..stage as usize {
                        cv.push(ark_mem.h * step_mem.MRIC.as_ref().unwrap().G[0][stage as usize][j]);
                    }
                    let ARKodeMRIStepMem { Fse, Fsi, unify_Fs, zpred, sdata, .. } =
                        &mut **step_mem;
                    let fsi_arr = if *unify_Fs { &*Fse } else { &*Fsi };
                    let mut xr: Vec<&NVector> = Vec::with_capacity(stage as usize + 2);
                    xr.push(&ark_mem.ycur);
                    xr.push(zpred);
                    for f in fsi_arr.iter().take(stage as usize) {
                        xr.push(f);
                    }
                    N_VLinearCombination(cv.len() as i32, &cv, &xr, sdata);
                }

                /* Update gamma for implicit solver */
                step_mem.gamma =
                    ark_mem.h * step_mem.MRIC.as_ref().unwrap().G[0][stage as usize][stage as usize];
                if ark_mem.firststage {
                    step_mem.gammap = step_mem.gamma;
                }
                step_mem.gamrat = if ark_mem.firststage {
                    ONE
                } else {
                    step_mem.gamma / step_mem.gammap
                };

                /* perform implicit solve (result is stored in ark_mem->ycur);
                return with positive value on anything but success */
                *nflagPtr = crate::arkode_mristep_nls::mriStep_Nls(ark_mem, step_mem, *nflagPtr);
                if *nflagPtr != ARK_SUCCESS {
                    return TRY_AGAIN;
                }
            }
            /* perform explicit update for correction */
            else {
                let mut cv: Vec<f64> = Vec::with_capacity(stage as usize + 1);
                cv.push(ONE);
                for j in 0..stage as usize {
                    cv.push(ark_mem.h * step_mem.MRIC.as_ref().unwrap().G[0][stage as usize][j]);
                }
                /* z == X[0] with c0 == 1: in-place accumulate form */
                let ARKodeMRIStepMem { Fse, Fsi, unify_Fs, .. } = &mut **step_mem;
                let fsi_arr = if *unify_Fs { &*Fse } else { &*Fsi };
                let mut xr: Vec<&NVector> = Vec::with_capacity(stage as usize);
                for f in fsi_arr.iter().take(stage as usize) {
                    xr.push(f);
                }
                mri_lincomb_accumulate(&cv[1..], &xr, &mut ark_mem.ycur);
            }
        }

        /* apply user-supplied stage or step postprocessing function (if
        supplied), and reset the inner integrator with the modified stage
        solution */
        if let (Some(post), true) = (ark_mem.PostProcessStageFn, !solution && !embedding) {
            {
                let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
                let retval = post(*tcur, ycur, user_data);
                if retval != 0 {
                    return ARK_POSTPROCESS_STAGE_FAIL;
                }
            }
            let retval =
                mriStepInnerStepper_Reset(&mut step_mem.stepper, ark_mem.tcur, &ark_mem.ycur);
            if retval != ARK_SUCCESS {
                arkProcessError(
                    Some(ark_mem),
                    ARK_INNERSTEP_FAIL,
                    line!(),
                    "mriStep_TakeStepMRISR",
                    file!(),
                    "Unable to reset the inner stepper",
                );
                return ARK_INNERSTEP_FAIL;
            }
        } else if let (Some(post), true) = (ark_mem.PostProcessStepFn, solution) {
            {
                let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
                let retval = post(*tcur, ycur, user_data);
                if retval != 0 {
                    return ARK_POSTPROCESS_STEP_FAIL;
                }
            }
            let retval =
                mriStepInnerStepper_Reset(&mut step_mem.stepper, ark_mem.tcur, &ark_mem.ycur);
            if retval != ARK_SUCCESS {
                arkProcessError(
                    Some(ark_mem),
                    ARK_INNERSTEP_FAIL,
                    line!(),
                    "mriStep_TakeStepMRISR",
                    file!(),
                    "Unable to reset the inner stepper",
                );
                return ARK_INNERSTEP_FAIL;
            }
        }

        /* Compute updated slow RHS (except for final solution or embedding) */
        if !solution && !embedding {
            /* call the user-supplied pre-RHS function (if supplied) */
            if let Some(pre_rhs) = ark_mem.PreRhsFn {
                if step_mem.explicit_rhs
                    || (step_mem.implicit_rhs && (!step_mem.deduce_rhs || !impl_corr))
                {
                    let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
                    let retval = pre_rhs(*tcur, ycur, user_data);
                    if retval != 0 {
                        return ARK_PRERHSFN_FAIL;
                    }
                }
            }

            /* store implicit slow rhs */
            if step_mem.implicit_rhs {
                if !step_mem.deduce_rhs || !impl_corr {
                    let fsi = step_mem.fsi.unwrap();
                    let retval = {
                        let ARKodeMRIStepMem { Fse, Fsi, unify_Fs, .. } = &mut **step_mem;
                        let fsi_arr = if *unify_Fs { Fse } else { Fsi };
                        let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
                        fsi(*tcur, ycur, &mut fsi_arr[stage as usize], user_data)
                    };
                    step_mem.nfsi += 1;

                    if retval < 0 {
                        return ARK_RHSFUNC_FAIL;
                    }
                    if retval > 0 {
                        return ARK_UNREC_RHSFUNC_ERR;
                    }

                    /* Add external forcing to Fsi[stage], if applicable */
                    if step_mem.impforcing {
                        let vals = mriStep_ApplyForcing_coeffs(step_mem, ark_mem.tcur, ONE);
                        let ARKodeMRIStepMem { Fse, Fsi, unify_Fs, forcing, .. } =
                            &mut **step_mem;
                        let fsi_arr = if *unify_Fs { Fse } else { Fsi };
                        for (k, val) in vals.iter().enumerate() {
                            for e in 0..fsi_arr[stage as usize].data.len() {
                                fsi_arr[stage as usize].data[e] += val * forcing[k].data[e];
                            }
                        }
                    }
                } else {
                    let gamma = step_mem.gamma;
                    let ARKodeMRIStepMem { Fse, Fsi, unify_Fs, zcor, sdata, .. } =
                        &mut **step_mem;
                    let fsi_arr = if *unify_Fs { Fse } else { Fsi };
                    N_VLinearSum(
                        ONE / gamma,
                        zcor,
                        -ONE / gamma,
                        sdata,
                        &mut fsi_arr[stage as usize],
                    );
                }
            }

            /* store explicit slow rhs */
            if step_mem.explicit_rhs {
                let fse = step_mem.fse.unwrap();
                let retval = {
                    let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
                    fse(*tcur, ycur, &mut step_mem.Fse[stage as usize], user_data)
                };
                step_mem.nfse += 1;

                if retval < 0 {
                    return ARK_RHSFUNC_FAIL;
                }
                if retval > 0 {
                    return ARK_UNREC_RHSFUNC_ERR;
                }

                /* Add external forcing to Fse[stage], if applicable */
                if step_mem.expforcing {
                    let vals = mriStep_ApplyForcing_coeffs(step_mem, ark_mem.tcur, ONE);
                    let ARKodeMRIStepMem { Fse, forcing, .. } = &mut **step_mem;
                    for (k, val) in vals.iter().enumerate() {
                        for e in 0..Fse[stage as usize].data.len() {
                            Fse[stage as usize].data[e] += val * forcing[k].data[e];
                        }
                    }
                }
            }

            /* combine both RHS into Fse for ImEx problems since
               fast forcing function only depends on Omega coefficients */
            if step_mem.implicit_rhs && step_mem.explicit_rhs {
                let ARKodeMRIStepMem { Fse, Fsi, .. } = &mut **step_mem;
                Fse[stage as usize].linear_sum_with(ONE, ONE, &Fsi[stage as usize]);
            }
        }

        /* If this is the solution stage, archive for error estimation */
        if solution {
            let ARKodeMem { ycur, tempv4, .. } = ark_mem;
            N_VScale(ONE, ycur, tempv4);
        }
    } /* loop over stages */

    /* if temporal error estimation is enabled: compute estimate via difference
       between step solution and embedding, store in ark_mem->tempv1, store
       norm in dsmPtr, and copy solution back to ycur */
    if !ark_mem.fixedstep || ark_mem.AccumErrorType != ARK_ACCUMERROR_NONE {
        {
            let ARKodeMem { ycur, tempv4, tempv1, .. } = ark_mem;
            N_VLinearSum(ONE, tempv4, -ONE, ycur, tempv1);
        }
        *dsmPtr = N_VWrmsNorm(&ark_mem.tempv1, &ark_mem.ewt);
        let ARKodeMem { ycur, tempv4, .. } = ark_mem;
        N_VScale(ONE, tempv4, ycur);
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_TakeStepMERK:

  This routine performs a single MERK step.
  ---------------------------------------------------------------*/
fn mriStep_TakeStepMERK(ark_mem: &mut ARKodeMem, dsmPtr: &mut f64, nflagPtr: &mut i32) -> i32 {
    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_TakeStepMERK") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };
    let retval = mriStep_TakeStepMERK_inner(ark_mem, &mut step_mem, dsmPtr, nflagPtr);
    ark_mem.step_mem = Some(step_mem);
    retval
}

fn mriStep_TakeStepMERK_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeMRIStepMem>,
    dsmPtr: &mut f64,
    nflagPtr: &mut i32,
) -> i32 {
    use crate::arkode_impl::ARK_ACCUMERROR_NONE;

    /* initialize algebraic solver convergence flag to success;
       error estimate to zero */
    *nflagPtr = ARK_SUCCESS;
    *dsmPtr = ZERO;

    /* initial time for step (set at the top of each stage group) */
    let mut t0: f64;

    /* initialize the current stage index */
    step_mem.istage = 0;
    step_mem.cur_stage = 0;

    /* if MRI adaptivity is enabled: reset fast accumulated error,
       and send appropriate control parameter to the fast integrator */
    let adapt_type = match ark_mem
        .hadapt_mem
        .as_ref()
        .and_then(|h| h.hcontroller.as_ref())
    {
        Some(hc) => SUNAdaptController_GetType(hc),
        None => SUN_ADAPTCONTROLLER_NONE,
    };
    let mut need_inner_dsm = false;
    if adapt_type == SUN_ADAPTCONTROLLER_MRI_H_TOL {
        need_inner_dsm = true;
        step_mem.inner_dsm = ZERO;
        let retval = mriStepInnerStepper_ResetAccumulatedError(&mut step_mem.stepper);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_INNERSTEP_FAIL,
                line!(),
                "mriStep_TakeStepMERK",
                file!(),
                "Unable to reset the inner stepper error estimate",
            );
            return ARK_INNERSTEP_FAIL;
        }
        let retval = mriStepInnerStepper_SetRTol(
            &mut step_mem.stepper,
            step_mem.inner_rtol_factor * ark_mem.reltol,
        );
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_INNERSTEP_FAIL,
                line!(),
                "mriStep_TakeStepMERK",
                file!(),
                "Unable to set the inner stepper tolerance",
            );
            return ARK_INNERSTEP_FAIL;
        }
    }

    /* for adaptive computations, reset the inner integrator to the beginning
    of this step */
    if !ark_mem.fixedstep {
        let retval =
            mriStepInnerStepper_Reset(&mut step_mem.stepper, ark_mem.tcur, &ark_mem.ycur);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_INNERSTEP_FAIL,
                line!(),
                "mriStep_TakeStepMERK",
                file!(),
                "Unable to reset the inner stepper",
            );
            return ARK_INNERSTEP_FAIL;
        }
    }

    /* Evaluate the slow RHS function if needed. */
    let nested_mri = step_mem.expforcing || step_mem.impforcing;
    if ark_mem.fn_.data.is_empty() || nested_mri {
        let ycur = std::mem::take(&mut ark_mem.ycur);
        let retval = mriStep_UpdateF0(ark_mem, step_mem, ark_mem.tcur, &ycur, ARK_FULLRHS_START);
        ark_mem.ycur = ycur;
        if retval != 0 {
            return ARK_RHSFUNC_FAIL;
        }
    } else if !ark_mem.fn_.data.is_empty() && !ark_mem.fn_is_current {
        let ycur = std::mem::take(&mut ark_mem.ycur);
        let mut fn_ = std::mem::take(&mut ark_mem.fn_);
        let retval = mriStep_FullRHS_inner(
            ark_mem,
            step_mem,
            ark_mem.tcur,
            &ycur,
            &mut fn_,
            ARK_FULLRHS_START,
        );
        ark_mem.ycur = ycur;
        ark_mem.fn_ = fn_;
        if retval != 0 {
            return ARK_RHSFUNC_FAIL;
        }
    }
    ark_mem.fn_is_current = true;

    /* The first stage is the previous time-step solution, so its RHS
       is the [already-computed] slow RHS from the start of the step */

    /* Loop over stage groups */
    for ig in 0..step_mem.MRIC.as_ref().unwrap().ngroup {
        /* Find the lowest stage number in this group. */
        let mut lowest_stage = step_mem.MRIC.as_ref().unwrap().group[ig as usize][0];
        for il in 1..step_mem.MRIC.as_ref().unwrap().stages {
            if step_mem.MRIC.as_ref().unwrap().group[ig as usize][il as usize] < 0 {
                break;
            }
            lowest_stage = std::cmp::min(
                lowest_stage,
                step_mem.MRIC.as_ref().unwrap().group[ig as usize][il as usize],
            );
        }

        /* Set up fast RHS for this stage group */
        let retval = mriStep_ComputeInnerForcing(
            ark_mem,
            step_mem,
            lowest_stage,
            ark_mem.tn,
            ark_mem.tn + ark_mem.h,
        );
        if retval != ARK_SUCCESS {
            return retval;
        }

        /* Set initial condition for this stage group (all but first group) */
        if ig > 0 {
            let ARKodeMem { yn, ycur, .. } = ark_mem;
            N_VScale(ONE, yn, ycur);
        }
        t0 = ark_mem.tn;

        /* Evolve fast IVP over each subinterval in stage group */
        for is in 0..step_mem.stages {
            /* Get stage index from group; skip to the next group if
               we've reached the end of this one */
            let stage = step_mem.MRIC.as_ref().unwrap().group[ig as usize][is as usize];
            step_mem.istage = stage;
            step_mem.cur_stage = stage;
            if stage < 0 {
                break;
            }
            let mut nextstage = -1;
            if stage < step_mem.stages {
                nextstage =
                    step_mem.MRIC.as_ref().unwrap().group[ig as usize][(is + 1) as usize];
            }

            /* Determine if this is an "embedding" or "solution" stage */
            let mut embedding = false;
            let mut solution = false;
            if ig == step_mem.MRIC.as_ref().unwrap().ngroup - 2 && stage >= 0 && nextstage < 0 {
                embedding = true;
            }
            if ig == step_mem.MRIC.as_ref().unwrap().ngroup - 1 && stage >= 0 && nextstage < 0 {
                solution = true;
            }

            /* Skip the embedding if we're using fixed time-stepping and
               temporal error estimation is disabled */
            if ark_mem.fixedstep && embedding && ark_mem.AccumErrorType == ARK_ACCUMERROR_NONE {
                break;
            }

            /* Set current stage abscissa */
            let cstage = if stage >= step_mem.stages {
                ONE
            } else {
                step_mem.MRIC.as_ref().unwrap().c[stage as usize]
            };

            /* Set desired output time for subinterval */
            let tf = ark_mem.tn + cstage * ark_mem.h;

            /* Reset the inner stepper on the first stage within all but the
               first stage group due to "stage-restart" structure */
            if stage > 1 && is == 0 {
                let retval =
                    mriStepInnerStepper_Reset(&mut step_mem.stepper, t0, &ark_mem.ycur);
                if retval != ARK_SUCCESS {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_INNERSTEP_FAIL,
                        line!(),
                        "mriStep_TakeStepMERK",
                        file!(),
                        "Unable to reset the inner stepper",
                    );
                    return ARK_INNERSTEP_FAIL;
                }
            }

            /* Evolve fast IVP for this stage, potentially get inner dsm on all
               non-embedding stages */
            let retval =
                mriStep_StageERKFast(ark_mem, step_mem, t0, tf, need_inner_dsm && !embedding);
            if retval != ARK_SUCCESS {
                *nflagPtr = CONV_FAIL;
                return retval;
            }

            /* Update "initial time" for next stage in group */
            t0 = tf;

            /* set current stage time for postprocessing and RHS calls */
            ark_mem.tcur = tf;

            /* apply user-supplied stage postprocessing function (if supplied),
               and reset the inner integrator with the modified stage solution */
            if let (Some(post), true) =
                (ark_mem.PostProcessStageFn, !solution && !embedding)
            {
                {
                    let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
                    let retval = post(*tcur, ycur, user_data);
                    if retval != 0 {
                        return ARK_POSTPROCESS_STAGE_FAIL;
                    }
                }
                let retval = mriStepInnerStepper_Reset(
                    &mut step_mem.stepper,
                    ark_mem.tcur,
                    &ark_mem.ycur,
                );
                if retval != ARK_SUCCESS {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_INNERSTEP_FAIL,
                        line!(),
                        "mriStep_TakeStepMERK",
                        file!(),
                        "Unable to reset the inner stepper",
                    );
                    return ARK_INNERSTEP_FAIL;
                }
            } else if let (Some(post), true) = (ark_mem.PostProcessStepFn, solution) {
                {
                    let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
                    let retval = post(*tcur, ycur, user_data);
                    if retval != 0 {
                        return ARK_POSTPROCESS_STEP_FAIL;
                    }
                }
                let retval = mriStepInnerStepper_Reset(
                    &mut step_mem.stepper,
                    ark_mem.tcur,
                    &ark_mem.ycur,
                );
                if retval != ARK_SUCCESS {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_INNERSTEP_FAIL,
                        line!(),
                        "mriStep_TakeStepMERK",
                        file!(),
                        "Unable to reset the inner stepper",
                    );
                    return ARK_INNERSTEP_FAIL;
                }
            }

            /* Compute updated slow RHS (except for final solution or
            embedding) */
            if !solution && !embedding {
                /* call the user-supplied pre-RHS function (if supplied) */
                if let Some(pre_rhs) = ark_mem.PreRhsFn {
                    let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
                    let retval = pre_rhs(*tcur, ycur, user_data);
                    if retval != 0 {
                        return ARK_PRERHSFN_FAIL;
                    }
                }

                /* store explicit slow rhs */
                let fse = step_mem.fse.unwrap();
                let retval = {
                    let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
                    fse(*tcur, ycur, &mut step_mem.Fse[stage as usize], user_data)
                };
                step_mem.nfse += 1;

                if retval < 0 {
                    return ARK_RHSFUNC_FAIL;
                }
                if retval > 0 {
                    return ARK_UNREC_RHSFUNC_ERR;
                }

                /* Add external forcing to Fse[stage], if applicable */
                if step_mem.expforcing {
                    let vals = mriStep_ApplyForcing_coeffs(step_mem, ark_mem.tcur, ONE);
                    let ARKodeMRIStepMem { Fse, forcing, .. } = &mut **step_mem;
                    for (k, val) in vals.iter().enumerate() {
                        for e in 0..Fse[stage as usize].data.len() {
                            Fse[stage as usize].data[e] += val * forcing[k].data[e];
                        }
                    }
                }
            }

            /* If this is the embedding stage, archive solution for error
            estimation */
            if embedding {
                let ARKodeMem { ycur, tempv4, .. } = ark_mem;
                N_VScale(ONE, ycur, tempv4);
            }
        } /* loop over stages */
    } /* loop over stage groups */

    /* if temporal error estimation is enabled: compute estimate via difference
       between step solution and embedding, store in ark_mem->tempv1, and store
       norm in dsmPtr */
    if !ark_mem.fixedstep || ark_mem.AccumErrorType != ARK_ACCUMERROR_NONE {
        {
            let ARKodeMem { ycur, tempv4, tempv1, .. } = ark_mem;
            N_VLinearSum(ONE, tempv4, -ONE, ycur, tempv1);
        }
        *dsmPtr = N_VWrmsNorm(&ark_mem.tempv1, &ark_mem.ewt);
    }

    ARK_SUCCESS
}

/*===============================================================
  Internal utility routines
  ===============================================================*/

/*---------------------------------------------------------------
  mriStep_SetCoupling

  This routine determines the MRI method to use, based on the
  desired accuracy and fixed/adaptive time stepping choice.
  ---------------------------------------------------------------*/
fn mriStep_SetCoupling(ark_mem: &mut ARKodeMem, step_mem: &mut Box<ARKodeMRIStepMem>) -> i32 {
    /* if coupling has already been specified, just return */
    if step_mem.MRIC.is_some() {
        return ARK_SUCCESS;
    }

    let mut table_id = ARKODE_MRI_NONE;

    /* select method based on order and type */
    if ark_mem.fixedstep {
        /**** fixed-step methods ****/
        if step_mem.implicit_rhs && step_mem.explicit_rhs {
            /**** ImEx methods ****/
            match step_mem.q {
                1 => table_id = MRISTEP_DEFAULT_IMEX_SD_1,
                2 => table_id = MRISTEP_DEFAULT_IMEX_SD_2,
                3 => table_id = MRISTEP_DEFAULT_IMEX_SD_3,
                4 => table_id = MRISTEP_DEFAULT_IMEX_SD_4,
                _ => {}
            }
        } else if step_mem.implicit_rhs {
            /**** implicit methods ****/
            match step_mem.q {
                1 => table_id = MRISTEP_DEFAULT_IMPL_SD_1,
                2 => table_id = MRISTEP_DEFAULT_IMPL_SD_2,
                3 => table_id = MRISTEP_DEFAULT_IMPL_SD_3,
                4 => table_id = MRISTEP_DEFAULT_IMPL_SD_4,
                _ => {}
            }
        } else {
            /**** explicit methods ****/
            match step_mem.q {
                1 => table_id = MRISTEP_DEFAULT_EXPL_1,
                2 => table_id = MRISTEP_DEFAULT_EXPL_2,
                3 => table_id = MRISTEP_DEFAULT_EXPL_3,
                4 => table_id = MRISTEP_DEFAULT_EXPL_4,
                5 => table_id = MRISTEP_DEFAULT_EXPL_5_AD,
                _ => {}
            }
        }
    } else {
        /**** adaptive methods ****/
        if step_mem.implicit_rhs && step_mem.explicit_rhs {
            /**** ImEx methods ****/
            match step_mem.q {
                2 => table_id = MRISTEP_DEFAULT_IMEX_SD_2_AD,
                3 => table_id = MRISTEP_DEFAULT_IMEX_SD_3_AD,
                4 => table_id = MRISTEP_DEFAULT_IMEX_SD_4_AD,
                _ => {}
            }
        } else if step_mem.implicit_rhs {
            /**** implicit methods ****/
            match step_mem.q {
                2 => table_id = MRISTEP_DEFAULT_IMPL_SD_2,
                3 => table_id = MRISTEP_DEFAULT_IMPL_SD_3,
                4 => table_id = MRISTEP_DEFAULT_IMPL_SD_4,
                _ => {}
            }
        } else {
            /**** explicit methods ****/
            match step_mem.q {
                2 => table_id = MRISTEP_DEFAULT_EXPL_2_AD,
                3 => table_id = MRISTEP_DEFAULT_EXPL_3_AD,
                4 => table_id = MRISTEP_DEFAULT_EXPL_4_AD,
                5 => table_id = MRISTEP_DEFAULT_EXPL_5_AD,
                _ => {}
            }
        }
    }
    if table_id == ARKODE_MRI_NONE {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "mriStep_SetCoupling",
            file!(),
            "No MRI method is available for the requested configuration.",
        );
        return ARK_ILL_INPUT;
    }

    step_mem.MRIC = MRIStepCoupling_LoadTable(table_id);
    if step_mem.MRIC.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_INVALID_TABLE,
            line!(),
            "mriStep_SetCoupling",
            file!(),
            "An error occurred in constructing coupling table.",
        );
        return ARK_INVALID_TABLE;
    }

    /* note coupling structure space requirements */
    let mut Cliw: i64 = 0;
    let mut Clrw: i64 = 0;
    MRIStepCoupling_Space(step_mem.MRIC.as_deref(), &mut Cliw, &mut Clrw);
    ark_mem.liw += Cliw;
    ark_mem.lrw += Clrw;

    /* set [redundant] stored values for stage numbers and
       method/embedding orders */
    let MRIC = step_mem.MRIC.as_ref().unwrap();
    step_mem.stages = MRIC.stages;
    step_mem.q = MRIC.q;
    step_mem.p = MRIC.p;

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_CheckCoupling

  This routine runs through the MRI coupling structure to ensure
  that it meets all necessary requirements.
  ---------------------------------------------------------------*/
fn mriStep_CheckCoupling(ark_mem: &mut ARKodeMem, step_mem: &mut Box<ARKodeMRIStepMem>) -> i32 {
    let tol: f64 = 100.0 * SUN_UNIT_ROUNDOFF;
    let MRIC = step_mem.MRIC.as_deref().unwrap();

    /* check that stages > 0 */
    if MRIC.stages < 1 {
        arkProcessError(
            Some(ark_mem),
            ARK_INVALID_TABLE,
            line!(),
            "mriStep_CheckCoupling",
            file!(),
            "stages < 1!",
        );
        return ARK_INVALID_TABLE;
    }

    /* check that method order q > 0 */
    if MRIC.q < 1 {
        arkProcessError(
            Some(ark_mem),
            ARK_INVALID_TABLE,
            line!(),
            "mriStep_CheckCoupling",
            file!(),
            "method order < 1",
        );
        return ARK_INVALID_TABLE;
    }

    /* check that embedding order p > 0 (if adaptive) */
    if MRIC.p < 1 && !ark_mem.fixedstep {
        arkProcessError(
            Some(ark_mem),
            ARK_INVALID_TABLE,
            line!(),
            "mriStep_CheckCoupling",
            file!(),
            "embedding order < 1",
        );
        return ARK_INVALID_TABLE;
    }

    /* Check that coupling table has compatible type */
    if step_mem.implicit_rhs
        && step_mem.explicit_rhs
        && MRIC.type_ != MRISTEP_IMEX
        && MRIC.type_ != MRISTEP_SR
    {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "mriStep_CheckCoupling",
            file!(),
            "Invalid coupling table for an IMEX problem!",
        );
        return ARK_ILL_INPUT;
    }
    if step_mem.explicit_rhs
        && MRIC.type_ != MRISTEP_EXPLICIT
        && MRIC.type_ != MRISTEP_IMEX
        && MRIC.type_ != MRISTEP_MERK
        && MRIC.type_ != MRISTEP_SR
    {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "mriStep_CheckCoupling",
            file!(),
            "Invalid coupling table for an explicit problem!",
        );
        return ARK_ILL_INPUT;
    }
    if step_mem.implicit_rhs
        && MRIC.type_ != MRISTEP_IMPLICIT
        && MRIC.type_ != MRISTEP_IMEX
        && MRIC.type_ != MRISTEP_SR
    {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "mriStep_CheckCoupling",
            file!(),
            "Invalid coupling table for an implicit problem!",
        );
        return ARK_ILL_INPUT;
    }

    /* Check that the matrices are defined appropriately */
    if MRIC.type_ == MRISTEP_IMEX || MRIC.type_ == MRISTEP_SR {
        /* ImEx */
        if MRIC.W.is_empty() || MRIC.G.is_empty() {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "mriStep_CheckCoupling",
                file!(),
                "Invalid coupling table for an IMEX problem!",
            );
            return ARK_ILL_INPUT;
        }
    } else if MRIC.type_ == MRISTEP_EXPLICIT || MRIC.type_ == MRISTEP_MERK {
        /* Explicit */
        if MRIC.W.is_empty() || !MRIC.G.is_empty() {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "mriStep_CheckCoupling",
                file!(),
                "Invalid coupling table for an explicit problem!",
            );
            return ARK_ILL_INPUT;
        }
    } else if MRIC.type_ == MRISTEP_IMPLICIT {
        /* Implicit */
        if !MRIC.W.is_empty() || MRIC.G.is_empty() {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "mriStep_CheckCoupling",
                file!(),
                "Invalid coupling table for an implicit problem!",
            );
            return ARK_ILL_INPUT;
        }
    }

    /* Check that W tables are strictly lower triangular */
    if !MRIC.W.is_empty() {
        let mut Wabs = ZERO;
        for k in 0..MRIC.nmat as usize {
            for i in 0..MRIC.stages as usize {
                for j in i..MRIC.stages as usize {
                    Wabs += SUNRabs(MRIC.W[k][i][j]);
                }
            }
        }
        if Wabs > tol {
            arkProcessError(
                Some(ark_mem),
                ARK_INVALID_TABLE,
                line!(),
                "mriStep_CheckCoupling",
                file!(),
                "Coupling can be up to ERK (at most)!",
            );
            return ARK_INVALID_TABLE;
        }
    }

    /* Check that G tables are lower triangular */
    if !MRIC.G.is_empty() {
        let mut Gabs = ZERO;
        for k in 0..MRIC.nmat as usize {
            for i in 0..MRIC.stages as usize {
                for j in (i + 1)..MRIC.stages as usize {
                    Gabs += SUNRabs(MRIC.G[k][i][j]);
                }
            }
        }
        if Gabs > tol {
            arkProcessError(
                Some(ark_mem),
                ARK_INVALID_TABLE,
                line!(),
                "mriStep_CheckCoupling",
                file!(),
                "Coupling can be up to DIRK (at most)!",
            );
            return ARK_INVALID_TABLE;
        }
    }

    /* Check that MERK "groups" are structured appropriately */
    if MRIC.type_ == MRISTEP_MERK {
        let mut group_counter = vec![0i32; (MRIC.stages + 1) as usize];
        for i in 0..MRIC.ngroup as usize {
            for j in 0..MRIC.stages as usize {
                let k = MRIC.group[i][j];
                if k == -1 {
                    break;
                }
                if k < 0 || k > MRIC.stages {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_INVALID_TABLE,
                        line!(),
                        "mriStep_CheckCoupling",
                        file!(),
                        "Invalid MERK group index!",
                    );
                    return ARK_INVALID_TABLE;
                }
                group_counter[k as usize] += 1;
            }
        }
        for i in 1..=MRIC.stages as usize {
            if group_counter[i] == 0 || group_counter[i] > 1 {
                arkProcessError(
                    Some(ark_mem),
                    ARK_INVALID_TABLE,
                    line!(),
                    "mriStep_CheckCoupling",
                    file!(),
                    "Duplicated/missing stages from MERK groups!",
                );
                return ARK_INVALID_TABLE;
            }
        }
    }

    /* Check that no stage has MRISTAGE_DIRK_FAST type (for now) */
    let mut okay = true;
    for i in 0..MRIC.stages {
        if mriStepCoupling_GetStageType(MRIC, i) == MRISTAGE_DIRK_FAST {
            okay = false;
        }
    }
    if !okay {
        arkProcessError(
            Some(ark_mem),
            ARK_INVALID_TABLE,
            line!(),
            "mriStep_CheckCoupling",
            file!(),
            "solve-coupled DIRK stages not currently supported",
        );
        return ARK_INVALID_TABLE;
    }

    /* check that MRI-GARK stage times are sorted */
    if MRIC.type_ == MRISTEP_IMPLICIT || MRIC.type_ == MRISTEP_EXPLICIT || MRIC.type_ == MRISTEP_IMEX
    {
        let mut okay = true;
        for i in 1..MRIC.stages as usize {
            if (MRIC.c[i] - MRIC.c[i - 1]) < -tol {
                okay = false;
            }
        }
        if !okay {
            arkProcessError(
                Some(ark_mem),
                ARK_INVALID_TABLE,
                line!(),
                "mriStep_CheckCoupling",
                file!(),
                "Stage times must be sorted.",
            );
            return ARK_INVALID_TABLE;
        }
    }

    /* check that the first stage is just the old step solution */
    let mut Gabs = SUNRabs(MRIC.c[0]);
    for k in 0..MRIC.nmat as usize {
        for j in 0..MRIC.stages as usize {
            if !MRIC.W.is_empty() {
                Gabs += SUNRabs(MRIC.W[k][0][j]);
            }
            if !MRIC.G.is_empty() {
                Gabs += SUNRabs(MRIC.G[k][0][j]);
            }
        }
    }
    if Gabs > tol {
        arkProcessError(
            Some(ark_mem),
            ARK_INVALID_TABLE,
            line!(),
            "mriStep_CheckCoupling",
            file!(),
            "First stage must equal old solution.",
        );
        return ARK_INVALID_TABLE;
    }

    /* check that the last stage is at the final time */
    if SUNRabs(ONE - MRIC.c[(MRIC.stages - 1) as usize]) > tol {
        arkProcessError(
            Some(ark_mem),
            ARK_INVALID_TABLE,
            line!(),
            "mriStep_CheckCoupling",
            file!(),
            "Final stage time must be equal 1.",
        );
        return ARK_INVALID_TABLE;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_StageERKFast

  This routine performs a single MRI stage, is, with explicit
  slow time scale and fast time scale that requires evolution.

  Operates on ark_mem->ycur (C passes ycur as an argument in every
  call; the unused ytemp argument is dropped).
  ---------------------------------------------------------------*/
fn mriStep_StageERKFast(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeMRIStepMem>,
    t0: f64,
    tf: f64,
    get_inner_dsm: bool,
) -> i32 {
    /* pre inner evolve function (if supplied) */
    if let Some(pre) = step_mem.pre_inner_evolve {
        let stepper = &step_mem.stepper;
        let retval = pre(t0, &stepper.forcing, stepper.nforcing, &mut ark_mem.user_data);
        if retval != 0 {
            return ARK_OUTERTOINNER_FAIL;
        }
    }

    /* Get the adaptivity type (if applicable) */
    let adapt_type = if get_inner_dsm {
        match ark_mem
            .hadapt_mem
            .as_ref()
            .and_then(|h| h.hcontroller.as_ref())
        {
            Some(hc) => SUNAdaptController_GetType(hc),
            None => SUN_ADAPTCONTROLLER_NONE,
        }
    } else {
        SUN_ADAPTCONTROLLER_NONE
    };

    /* advance inner method in time */
    let mut ycur = std::mem::take(&mut ark_mem.ycur);
    let retval = mriStepInnerStepper_Evolve(&mut step_mem.stepper, t0, tf, &mut ycur);
    ark_mem.ycur = ycur;

    if retval < 0 {
        arkProcessError(
            Some(ark_mem),
            ARK_INNERSTEP_FAIL,
            line!(),
            "mriStep_StageERKFast",
            file!(),
            "Failure when evolving the inner stepper",
        );
        return ARK_INNERSTEP_FAIL;
    }
    if retval > 0 {
        /* increment stepper-specific counter, and decrement ARKODE-level
           nonlinear solver counter (since that will be incremented
           automatically by ARKODE).  Return with "TRY_AGAIN" which should
           cause ARKODE to cut the step size and retry the step. */
        step_mem.inner_fails += 1;
        ark_mem.ncfn -= 1;
        return TRY_AGAIN;
    }

    /* for normal stages (i.e., not the embedding) with MRI adaptivity enabled,
       get an estimate for the fast time scale error */
    if get_inner_dsm {
        /* if the fast integrator uses adaptive steps, retrieve the error
        estimate */
        if adapt_type == SUN_ADAPTCONTROLLER_MRI_H_TOL {
            let mut inner_dsm = step_mem.inner_dsm;
            let retval =
                mriStepInnerStepper_GetAccumulatedError(&mut step_mem.stepper, &mut inner_dsm);
            step_mem.inner_dsm = inner_dsm;
            if retval != ARK_SUCCESS {
                arkProcessError(
                    Some(ark_mem),
                    ARK_INNERSTEP_FAIL,
                    line!(),
                    "mriStep_StageERKFast",
                    file!(),
                    "Unable to get accumulated error from the inner stepper",
                );
                return ARK_INNERSTEP_FAIL;
            }

            /* scale the error estimate by 1/rtol to account for different
            inner/outer tolerances */
            step_mem.inner_dsm /= ark_mem.reltol;
        }
    }

    /* post inner evolve function (if supplied) */
    if let Some(post) = step_mem.post_inner_evolve {
        let ARKodeMem { ycur, user_data, .. } = ark_mem;
        let retval = post(tf, ycur, user_data);
        if retval != 0 {
            return ARK_INNERTOOUTER_FAIL;
        }
    }

    /* return with success */
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_StageERKNoFast

  This routine performs a single MRI stage with explicit slow
  time scale only (no fast time scale evolution).
  ---------------------------------------------------------------*/
fn mriStep_StageERKNoFast(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeMRIStepMem>,
    is: i32,
) -> i32 {
    /* determine effective ERK coefficients (store in Ae_row and Ai_row) */
    {
        let ARKodeMRIStepMem { MRIC, stage_map, Ae_row, Ai_row, .. } = &mut **step_mem;
        let retval = mriStep_RKCoeffs(MRIC.as_deref().unwrap(), is, stage_map, Ae_row, Ai_row);
        if retval != ARK_SUCCESS {
            return retval;
        }
    }

    /* call fused vector operation to perform ERK update -- bound on
       j needs "SUNMIN" to handle the case of an "embedding" stage.
       z == X[0] with c0 == 1: in-place accumulate form */
    let mut cv: Vec<f64> = Vec::new();
    {
        let jmax = std::cmp::min(is, step_mem.stages);
        for j in 0..jmax as usize {
            if step_mem.explicit_rhs && step_mem.stage_map[j] > -1 {
                cv.push(ark_mem.h * step_mem.Ae_row[step_mem.stage_map[j] as usize]);
            }
            if step_mem.implicit_rhs && step_mem.stage_map[j] > -1 {
                cv.push(ark_mem.h * step_mem.Ai_row[step_mem.stage_map[j] as usize]);
            }
        }
    }
    {
        let jmax = std::cmp::min(is, step_mem.stages);
        let ARKodeMRIStepMem {
            Fse,
            Fsi,
            unify_Fs,
            stage_map,
            explicit_rhs,
            implicit_rhs,
            ..
        } = &mut **step_mem;
        let mut xr: Vec<&NVector> = Vec::new();
        for j in 0..jmax as usize {
            if *explicit_rhs && stage_map[j] > -1 {
                xr.push(&Fse[stage_map[j] as usize]);
            }
            if *implicit_rhs && stage_map[j] > -1 {
                if *unify_Fs {
                    xr.push(&Fse[stage_map[j] as usize]);
                } else {
                    xr.push(&Fsi[stage_map[j] as usize]);
                }
            }
        }
        mri_lincomb_accumulate(&cv, &xr, &mut ark_mem.ycur);
    }
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_StageDIRKFast

  This routine performs a single stage of a "solve coupled"
  MRI method -- not currently implemented.
  ---------------------------------------------------------------*/
fn mriStep_StageDIRKFast(
    ark_mem: &mut ARKodeMem,
    _step_mem: &mut Box<ARKodeMRIStepMem>,
    _is: i32,
    _nflagPtr: &mut i32,
) -> i32 {
    /* this is not currently implemented */
    arkProcessError(
        Some(ark_mem),
        ARK_INVALID_TABLE,
        line!(),
        "mriStep_StageDIRKFast",
        file!(),
        "This routine is not yet implemented.",
    );
    ARK_INVALID_TABLE
}

/*---------------------------------------------------------------
  mriStep_StageDIRKNoFast

  This routine performs a single MRI stage with implicit slow
  time scale only (no fast time scale evolution).
  ---------------------------------------------------------------*/
fn mriStep_StageDIRKNoFast(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeMRIStepMem>,
    is: i32,
    nflagPtr: &mut i32,
) -> i32 {
    /* store current stage index (for an "embedded" stage, subtract 1) */
    step_mem.istage = if is == step_mem.stages { is - 1 } else { is };

    /* Call predictor for current stage solution (result placed in zpred) */
    let istage = step_mem.istage;
    let mut zpred = std::mem::take(&mut step_mem.zpred);
    let retval = mriStep_Predict(ark_mem, step_mem, istage, &mut zpred);
    step_mem.zpred = zpred;
    if retval != ARK_SUCCESS {
        return retval;
    }

    /* If a user-supplied predictor routine is provided, call that here */
    if let Some(stage_predict) = step_mem.stage_predict {
        let retval = stage_predict(ark_mem.tcur, &mut step_mem.zpred, &mut ark_mem.user_data);
        if retval < 0 {
            return ARK_USER_PREDICT_FAIL;
        }
        if retval > 0 {
            return TRY_AGAIN;
        }
    }

    /* determine effective DIRK coefficients (store in cvals) */
    {
        let ARKodeMRIStepMem { MRIC, stage_map, Ae_row, Ai_row, .. } = &mut **step_mem;
        let retval = mriStep_RKCoeffs(MRIC.as_deref().unwrap(), is, stage_map, Ae_row, Ai_row);
        if retval != ARK_SUCCESS {
            return retval;
        }
    }

    /* Set up data for evaluation of DIRK stage residual (data stored in sdata) */
    let retval = mriStep_StageSetup(ark_mem, step_mem);
    if retval != ARK_SUCCESS {
        return retval;
    }

    /* perform implicit solve (result is stored in ark_mem->ycur); return
       with positive value on anything but success */
    *nflagPtr = crate::arkode_mristep_nls::mriStep_Nls(ark_mem, step_mem, *nflagPtr);
    if *nflagPtr != ARK_SUCCESS {
        return TRY_AGAIN;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_ComputeInnerForcing

  Constructs the 'coefficient' vectors for the forcing polynomial
  for a 'fast' outer MRI-GARK stage i (see the C source for the
  full derivation).
  ---------------------------------------------------------------*/
pub(crate) fn mriStep_ComputeInnerForcing(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeMRIStepMem>,
    stage: i32,
    t0: f64,
    tf: f64,
) -> i32 {
    let mut implicit_rhs = step_mem.implicit_rhs;
    let mut explicit_rhs = step_mem.explicit_rhs;

    /* Set inner forcing time normalization constants */
    step_mem.stepper.tshift = t0;
    step_mem.stepper.tscale = tf - t0;

    /* Adjust implicit/explicit RHS flags for MRISR methods, since these
       ignore the G coefficients in the forcing function */
    if step_mem.MRIC.as_ref().unwrap().type_ == MRISTEP_SR {
        implicit_rhs = false;
        explicit_rhs = true;
    }

    let nmat = step_mem.MRIC.as_ref().unwrap().nmat;
    let rcdiff = ark_mem.h / (tf - t0);
    let jmax = std::cmp::min(stage, step_mem.stages) as usize;

    let ARKodeMRIStepMem { Fse, Fsi, unify_Fs, MRIC, stage_map, stepper, .. } = &mut **step_mem;
    let MRIC = MRIC.as_deref().unwrap();

    /* compute inner forcing vectors (assumes cdiff != 0) */
    let mut Xvecs: Vec<&NVector> = Vec::new();
    for j in 0..jmax {
        if explicit_rhs && stage_map[j] > -1 {
            Xvecs.push(&Fse[stage_map[j] as usize]);
        }
        if implicit_rhs && stage_map[j] > -1 {
            if *unify_Fs {
                Xvecs.push(&Fse[stage_map[j] as usize]);
            } else {
                Xvecs.push(&Fsi[stage_map[j] as usize]);
            }
        }
    }

    for k in 0..nmat as usize {
        let mut cvals: Vec<f64> = Vec::new();
        for j in 0..jmax {
            if stage_map[j] > -1 {
                if explicit_rhs && implicit_rhs {
                    /* ImEx */
                    cvals.push(rcdiff * MRIC.W[k][stage as usize][j]);
                    cvals.push(rcdiff * MRIC.G[k][stage as usize][j]);
                } else if explicit_rhs {
                    /* explicit only */
                    cvals.push(rcdiff * MRIC.W[k][stage as usize][j]);
                } else {
                    /* implicit only */
                    cvals.push(rcdiff * MRIC.G[k][stage as usize][j]);
                }
            }
        }

        N_VLinearCombination(cvals.len() as i32, &cvals, &Xvecs, &mut stepper.forcing[k]);
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  Compute/return the effective RK coefficients for a "nofast"
  stage.  We may assume that "A" has already been allocated.
  ---------------------------------------------------------------*/
fn mriStep_RKCoeffs(
    MRIC: &crate::arkode_mri_tables::MRIStepCouplingMem,
    is: i32,
    stage_map: &[i32],
    Ae_row: &mut [f64],
    Ai_row: &mut [f64],
) -> i32 {
    if is < 1 || is > MRIC.stages {
        return ARK_INVALID_TABLE;
    }

    /* initialize RK coefficient array */
    for j in 0..MRIC.stages as usize {
        Ae_row[j] = ZERO;
        Ai_row[j] = ZERO;
    }

    /* compute RK coefficients -- note that bounds on j need
       "SUNMIN" to handle the case of an "embedding" stage */
    for k in 0..MRIC.nmat as usize {
        let kconst = ONE / (k as f64 + ONE);
        if !MRIC.W.is_empty() {
            let jmax = std::cmp::min(is, MRIC.stages - 1);
            for j in 0..jmax as usize {
                if stage_map[j] > -1 {
                    Ae_row[stage_map[j] as usize] += MRIC.W[k][is as usize][j] * kconst;
                }
            }
        }
        if !MRIC.G.is_empty() {
            let jmax = std::cmp::min(is, MRIC.stages - 1);
            for j in 0..=jmax as usize {
                if stage_map[j] > -1 {
                    Ai_row[stage_map[j] as usize] += MRIC.G[k][is as usize][j] * kconst;
                }
            }
        }
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_Predict

  This routine computes the prediction for a specific internal
  stage solution, storing the result in yguess.
  ---------------------------------------------------------------*/
pub(crate) fn mriStep_Predict(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeMRIStepMem>,
    istage: i32,
    yguess: &mut NVector,
) -> i32 {
    /* verify that interpolation structure is provided */
    if ark_mem.interp.is_none() && step_mem.predictor > 0 {
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_NULL,
            line!(),
            "mriStep_Predict",
            file!(),
            "Interpolation structure is NULL",
        );
        return ARK_MEM_NULL;
    }

    /* if the first step (or if resized), use initial condition as guess */
    if ark_mem.initsetup {
        N_VScale(ONE, &ark_mem.yn, yguess);
        return ARK_SUCCESS;
    }

    /* set evaluation time tau as relative shift from previous successful time */
    let tau =
        step_mem.MRIC.as_ref().unwrap().c[istage as usize] * ark_mem.h / ark_mem.hold;

    /* use requested predictor formula */
    match step_mem.predictor {
        1 => {
            /***** Interpolatory Predictor 1 -- all to max order *****/
            let retval = arkPredict_MaximumOrder(ark_mem, tau, yguess);
            if retval != ARK_ILL_INPUT {
                return retval;
            }
        }
        2 => {
            /***** Interpolatory Predictor 2 -- decrease order w/ increasing
            level of extrapolation *****/
            let retval = arkPredict_VariableOrder(ark_mem, tau, yguess);
            if retval != ARK_ILL_INPUT {
                return retval;
            }
        }
        3 => {
            /***** Cutoff predictor: max order interpolatory output for stages
            "close" to previous step, first-order predictor for subsequent
            stages *****/
            let retval = arkPredict_CutoffOrder(ark_mem, tau, yguess);
            if retval != ARK_ILL_INPUT {
                return retval;
            }
        }
        4 => {
            /***** Bootstrap predictor: if any previous stage in step has
            nonzero c_i, construct a quadratic Hermite interpolant for
            prediction; otherwise use the trivial predictor. *****/

            /* determine if any previous stages in step meet criteria */
            let MRIC = step_mem.MRIC.as_deref().unwrap();
            let mut jstage: i32 = -1;
            for i in 0..istage {
                jstage = if MRIC.c[i as usize] != ZERO { i } else { jstage };
            }

            /* if using the trivial predictor, break */
            if jstage != -1 {
                /* find the "optimal" previous stage to use */
                for i in 0..istage {
                    if MRIC.c[i as usize] > MRIC.c[jstage as usize]
                        && MRIC.c[i as usize] != ZERO
                        && step_mem.stage_map[i as usize] > -1
                    {
                        jstage = i;
                    }
                }

                /* set stage time, stage RHS and interpolation values */
                let h = ark_mem.h * MRIC.c[jstage as usize];
                let tau = ark_mem.h * MRIC.c[istage as usize];
                let smap = step_mem.stage_map[jstage as usize] as usize;

                let mut cv: Vec<f64> = Vec::new();
                let ARKodeMRIStepMem {
                    Fse,
                    Fsi,
                    unify_Fs,
                    explicit_rhs,
                    implicit_rhs,
                    ..
                } = &mut **step_mem;
                let mut xr: Vec<&NVector> = Vec::new();
                if *implicit_rhs {
                    /* Implicit piece */
                    cv.push(ONE);
                    if *unify_Fs {
                        xr.push(&Fse[smap]);
                    } else {
                        xr.push(&Fsi[smap]);
                    }
                }
                if *explicit_rhs {
                    /* Explicit piece */
                    cv.push(ONE);
                    xr.push(&Fse[smap]);
                }

                /* call predictor routine */
                let retval =
                    arkPredict_Bootstrap(ark_mem, h, tau, cv.len(), &cv, &xr, yguess);
                if retval != ARK_ILL_INPUT {
                    return retval;
                }
            }
        }
        _ => {}
    }

    /* if we made it here, use the trivial predictor (previous step solution) */
    N_VScale(ONE, &ark_mem.yn, yguess);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_StageSetup

  This routine sets up the stage data for computing the
  solve-decoupled MRI stage residual, along with the step- and
  method-related factors gamma, gammap and gamrat.
  ---------------------------------------------------------------*/
pub(crate) fn mriStep_StageSetup(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeMRIStepMem>,
) -> i32 {
    /* Set shortcut to current stage index */
    let i = step_mem.istage;

    /* Update gamma (if the method contains an implicit component) */
    step_mem.gamma = ark_mem.h * step_mem.Ai_row[step_mem.stage_map[i as usize] as usize];

    if ark_mem.firststage {
        step_mem.gammap = step_mem.gamma;
    }
    step_mem.gamrat = if ark_mem.firststage {
        ONE
    } else {
        step_mem.gamma / step_mem.gammap
    };

    /* set cvals and Xvecs for setting stage data:
       sdata = ycur - zpred + h*sum_j (Ae_row[j]*Fse[j] + Ai_row[j]*Fsi[j])
       (dest sdata distinct from all operands -> free fused op) */
    let mut cv: Vec<f64> = vec![ONE, -ONE];
    for j in 0..i as usize {
        if step_mem.explicit_rhs && step_mem.stage_map[j] > -1 {
            cv.push(ark_mem.h * step_mem.Ae_row[step_mem.stage_map[j] as usize]);
        }
        if step_mem.implicit_rhs && step_mem.stage_map[j] > -1 {
            cv.push(ark_mem.h * step_mem.Ai_row[step_mem.stage_map[j] as usize]);
        }
    }

    {
        let ARKodeMRIStepMem {
            Fse,
            Fsi,
            unify_Fs,
            stage_map,
            explicit_rhs,
            implicit_rhs,
            zpred,
            sdata,
            ..
        } = &mut **step_mem;
        let mut xr: Vec<&NVector> = vec![&ark_mem.ycur, zpred];
        for j in 0..i as usize {
            if *explicit_rhs && stage_map[j] > -1 {
                xr.push(&Fse[stage_map[j] as usize]);
            }
            if *implicit_rhs && stage_map[j] > -1 {
                if *unify_Fs {
                    xr.push(&Fse[stage_map[j] as usize]);
                } else {
                    xr.push(&Fsi[stage_map[j] as usize]);
                }
            }
        }

        /* call fused vector operation to do the work */
        N_VLinearCombination(cv.len() as i32, &cv, &xr, sdata);
    }

    /* return with success */
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_SlowRHS:

  Wrapper routine to call the user-supplied slow RHS functions,
  f(t,y) = fse(t,y) + fsi(t,y), with API matching
  ARKTimestepFullRHSFn.  This is only used to determine an
  initial slow time-step size.
  ---------------------------------------------------------------*/
pub(crate) fn mriStep_SlowRHS_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeMRIStepMem>,
    t: f64,
    y: &NVector,
    f: &mut NVector,
    _mode: i32,
) -> i32 {
    /* call the user-supplied pre-RHS function (if supplied) */
    if let Some(pre_rhs) = ark_mem.PreRhsFn {
        let retval = pre_rhs(t, y, &mut ark_mem.user_data);
        if retval != 0 {
            return ARK_PRERHSFN_FAIL;
        }
    }

    /* call fsi if the problem has an implicit component */
    if step_mem.implicit_rhs {
        let fsi = step_mem.fsi.unwrap();
        let retval = {
            let ARKodeMRIStepMem { Fse, Fsi, unify_Fs, .. } = &mut **step_mem;
            let fsi_arr = if *unify_Fs { Fse } else { Fsi };
            fsi(t, y, &mut fsi_arr[0], &mut ark_mem.user_data)
        };
        step_mem.nfsi += 1;
        step_mem.fsi_is_current = true;
        if retval != 0 {
            arkProcessError(
                Some(ark_mem),
                ARK_RHSFUNC_FAIL,
                line!(),
                "mriStep_SlowRHS",
                file!(),
                &rhsfunc_failed_msg(t),
            );
            return ARK_RHSFUNC_FAIL;
        }

        /* Add external forcing, if applicable */
        if step_mem.impforcing {
            let vals = mriStep_ApplyForcing_coeffs(step_mem, t, ONE);
            let ARKodeMRIStepMem { Fse, Fsi, unify_Fs, forcing, .. } = &mut **step_mem;
            let fsi_arr = if *unify_Fs { Fse } else { Fsi };
            for (k, val) in vals.iter().enumerate() {
                for e in 0..fsi_arr[0].data.len() {
                    fsi_arr[0].data[e] += val * forcing[k].data[e];
                }
            }
        }
    }

    /* call fse if the problem has an explicit component */
    if step_mem.explicit_rhs {
        let fse = step_mem.fse.unwrap();
        let retval = fse(t, y, &mut step_mem.Fse[0], &mut ark_mem.user_data);
        step_mem.nfse += 1;
        step_mem.fse_is_current = true;
        if retval != 0 {
            arkProcessError(
                Some(ark_mem),
                ARK_RHSFUNC_FAIL,
                line!(),
                "mriStep_SlowRHS",
                file!(),
                &rhsfunc_failed_msg(t),
            );
            return ARK_RHSFUNC_FAIL;
        }

        /* Add external forcing, if applicable */
        if step_mem.expforcing {
            let vals = mriStep_ApplyForcing_coeffs(step_mem, t, ONE);
            let ARKodeMRIStepMem { Fse, forcing, .. } = &mut **step_mem;
            for (k, val) in vals.iter().enumerate() {
                for e in 0..Fse[0].data.len() {
                    Fse[0].data[e] += val * forcing[k].data[e];
                }
            }
        }
    }

    /* combine RHS vectors into output */
    if step_mem.explicit_rhs && step_mem.implicit_rhs {
        /* ImEx */
        let ARKodeMRIStepMem { Fse, Fsi, .. } = &mut **step_mem;
        N_VLinearSum(ONE, &Fse[0], ONE, &Fsi[0], f);
    } else if step_mem.implicit_rhs {
        let fsi0 = if step_mem.unify_Fs {
            &step_mem.Fse[0]
        } else {
            &step_mem.Fsi[0]
        };
        N_VScale(ONE, fsi0, f);
    } else {
        N_VScale(ONE, &step_mem.Fse[0], f);
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_Hin

  This routine computes a tentative initial step size h0.
  ---------------------------------------------------------------*/
fn mriStep_Hin(ark_mem: &mut ARKodeMem, tcur: f64, tout: f64, fcur: &NVector, h: &mut f64) -> i32 {
    /* If tout is too close to tn, give up */
    let tdiff = tout - tcur;
    if tdiff == ZERO {
        return ARK_TOO_CLOSE;
    }
    let sign = if tdiff > ZERO { 1.0 } else { -1.0 };
    let tdist = SUNRabs(tdiff);
    let tround = ark_mem.uround * SUNMAX(SUNRabs(tcur), SUNRabs(tout));
    if tdist < TWO * tround {
        return ARK_TOO_CLOSE;
    }

    /* h0 should bound the change due to a forward Euler step, and
       include safeguard against "too-small" ||f(t0,y0)||: */
    let fnorm = N_VWrmsNorm(fcur, &ark_mem.ewt) / H0_BIAS;
    let h0_inv = SUNMAX(ONE / H0_UBFACTOR / tdist, fnorm);
    *h = sign / h0_inv;
    ARK_SUCCESS
}

/*===============================================================
  User-callable functions for a custom inner integrator
  ===============================================================*/

/// C MRIStepInnerStepper_Create(sunctx, &stepper): allocation cannot
/// fail, so the object is returned directly.
pub fn MRIStepInnerStepper_Create(
    _sunctx: &crate::sundials_context::SUNContext,
) -> MRIStepInnerStepper {
    MRIStepInnerStepper {
        last_flag: ARK_SUCCESS,
        ..Default::default()
    }
}

/// C MRIStepInnerStepper_CreateFromSUNStepper: wraps a SUNStepper
/// (taking ownership of it) as an MRIStepInnerStepper.
pub fn MRIStepInnerStepper_CreateFromSUNStepper(sunstepper: SUNStepper) -> MRIStepInnerStepper {
    let mut stepper = MRIStepInnerStepper {
        last_flag: ARK_SUCCESS,
        ..Default::default()
    };
    let _ = MRIStepInnerStepper_SetContent(&mut stepper, Some(Box::new(sunstepper)));
    let _ = MRIStepInnerStepper_SetEvolveFn(&mut stepper, mriStepInnerStepper_EvolveSUNStepper);
    let _ = MRIStepInnerStepper_SetFullRhsFn(&mut stepper, mriStepInnerStepper_FullRhsSUNStepper);
    let _ = MRIStepInnerStepper_SetResetFn(&mut stepper, mriStepInnerStepper_ResetSUNStepper);
    stepper
}

/// C MRIStepInnerStepper_Free: the drop releases the forcing vectors
/// and the owned content.
pub fn MRIStepInnerStepper_Free(stepper: &mut Option<MRIStepInnerStepper>) -> i32 {
    *stepper = None;
    ARK_SUCCESS
}

pub fn MRIStepInnerStepper_SetContent(stepper: &mut MRIStepInnerStepper, content: UserData) -> i32 {
    stepper.content = content;
    ARK_SUCCESS
}

/// C copies the `void*` out; the Rust port hands back a mutable
/// borrow of the content slot.
pub fn MRIStepInnerStepper_GetContent(stepper: &mut MRIStepInnerStepper) -> &mut UserData {
    &mut stepper.content
}

pub fn MRIStepInnerStepper_SetEvolveFn(
    stepper: &mut MRIStepInnerStepper,
    f: MRIStepInnerEvolveFn,
) -> i32 {
    stepper.ops.evolve = Some(f);
    ARK_SUCCESS
}

pub fn MRIStepInnerStepper_SetFullRhsFn(
    stepper: &mut MRIStepInnerStepper,
    f: MRIStepInnerFullRhsFn,
) -> i32 {
    stepper.ops.fullrhs = Some(f);
    ARK_SUCCESS
}

pub fn MRIStepInnerStepper_SetResetFn(
    stepper: &mut MRIStepInnerStepper,
    f: MRIStepInnerResetFn,
) -> i32 {
    stepper.ops.reset = Some(f);
    ARK_SUCCESS
}

pub fn MRIStepInnerStepper_SetAccumulatedErrorGetFn(
    stepper: &mut MRIStepInnerStepper,
    f: MRIStepInnerGetAccumulatedError,
) -> i32 {
    stepper.ops.geterror = Some(f);
    ARK_SUCCESS
}

pub fn MRIStepInnerStepper_SetAccumulatedErrorResetFn(
    stepper: &mut MRIStepInnerStepper,
    f: MRIStepInnerResetAccumulatedError,
) -> i32 {
    stepper.ops.reseterror = Some(f);
    ARK_SUCCESS
}

pub fn MRIStepInnerStepper_SetRTolFn(
    stepper: &mut MRIStepInnerStepper,
    f: MRIStepInnerSetRTol,
) -> i32 {
    stepper.ops.setrtol = Some(f);
    ARK_SUCCESS
}

pub fn MRIStepInnerStepper_AddForcing(
    stepper: &mut MRIStepInnerStepper,
    t: f64,
    f: &mut NVector,
) -> i32 {
    /* compute normalized time tau and initialize tau^i */
    let tau = (t - stepper.tshift) / stepper.tscale;
    let mut taui = ONE;

    /* always append the constant forcing term; z == X[0] with c0 == 1:
       in-place accumulate form of the C fused op */
    for i in 0..stepper.nforcing as usize {
        let val = taui;
        for e in 0..f.data.len() {
            f.data[e] += val * stepper.forcing[i].data[e];
        }
        taui *= tau;
    }

    ARK_SUCCESS
}

pub fn MRIStepInnerStepper_GetForcingData<'a>(
    stepper: &'a MRIStepInnerStepper,
    tshift: &mut f64,
    tscale: &mut f64,
    forcing: &mut &'a [NVector],
    nforcing: &mut i32,
) -> i32 {
    *tshift = stepper.tshift;
    *tscale = stepper.tscale;
    *forcing = &stepper.forcing;
    *nforcing = stepper.nforcing;
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  Internal inner integrator functions
  ---------------------------------------------------------------*/

/* Check for required operations */
pub(crate) fn mriStepInnerStepper_HasRequiredOps(stepper: &MRIStepInnerStepper) -> i32 {
    if stepper.ops.evolve.is_some() {
        ARK_SUCCESS
    } else {
        ARK_ILL_INPUT
    }
}

/* Check whether stepper supports fast/slow tolerance adaptivity */
pub(crate) fn mriStepInnerStepper_SupportsRTolAdaptivity(stepper: &MRIStepInnerStepper) -> bool {
    stepper.ops.geterror.is_some()
        && stepper.ops.reseterror.is_some()
        && stepper.ops.setrtol.is_some()
}

/* Evolve the inner (fast) ODE */
pub(crate) fn mriStepInnerStepper_Evolve(
    stepper: &mut MRIStepInnerStepper,
    t0: f64,
    tout: f64,
    y: &mut NVector,
) -> i32 {
    let evolve = match stepper.ops.evolve {
        Some(f) => f,
        None => return ARK_ILL_INPUT,
    };
    stepper.last_flag = evolve(stepper, t0, tout, y);
    stepper.last_flag
}

pub(crate) fn mriStepInnerStepper_EvolveSUNStepper(
    stepper: &mut MRIStepInnerStepper,
    _t0: f64,
    tout: f64,
    y: &mut NVector,
) -> i32 {
    let MRIStepInnerStepper { content, forcing, nforcing, tshift, tscale, .. } = stepper;
    let sunstepper = match content.as_mut().and_then(|c| c.downcast_mut::<SUNStepper>()) {
        Some(s) => s,
        None => return crate::arkode_impl::ARK_SUNSTEPPER_ERR,
    };
    let mut tret = ZERO;

    let err = SUNStepper_SetForcing(sunstepper, *tshift, *tscale, forcing, *nforcing);
    stepper.last_flag = {
        let sunstepper = stepper
            .content
            .as_mut()
            .and_then(|c| c.downcast_mut::<SUNStepper>())
            .unwrap();
        sunstepper.last_flag
    };
    if err != SUN_SUCCESS {
        return crate::arkode_impl::ARK_SUNSTEPPER_ERR;
    }

    let sunstepper = stepper
        .content
        .as_mut()
        .and_then(|c| c.downcast_mut::<SUNStepper>())
        .unwrap();
    let err = SUNStepper_SetStopTime(sunstepper, tout);
    let lf = sunstepper.last_flag;
    stepper.last_flag = lf;
    if err != SUN_SUCCESS {
        return crate::arkode_impl::ARK_SUNSTEPPER_ERR;
    }

    let sunstepper = stepper
        .content
        .as_mut()
        .and_then(|c| c.downcast_mut::<SUNStepper>())
        .unwrap();
    let err = SUNStepper_Evolve(sunstepper, tout, y, &mut tret);
    let lf = sunstepper.last_flag;
    stepper.last_flag = lf;
    if err != SUN_SUCCESS {
        return crate::arkode_impl::ARK_SUNSTEPPER_ERR;
    }

    let sunstepper = stepper
        .content
        .as_mut()
        .and_then(|c| c.downcast_mut::<SUNStepper>())
        .unwrap();
    let err = SUNStepper_SetForcing(sunstepper, ZERO, ONE, &[], 0);
    let lf = sunstepper.last_flag;
    stepper.last_flag = lf;
    if err != SUN_SUCCESS {
        return crate::arkode_impl::ARK_SUNSTEPPER_ERR;
    }

    ARK_SUCCESS
}

/* Compute the full RHS for inner (fast) time scale */
pub(crate) fn mriStepInnerStepper_FullRhs(
    stepper: &mut MRIStepInnerStepper,
    t: f64,
    y: &NVector,
    f: &mut NVector,
    mode: i32,
) -> i32 {
    let fullrhs = match stepper.ops.fullrhs {
        Some(fr) => fr,
        None => return ARK_ILL_INPUT,
    };
    stepper.last_flag = fullrhs(stepper, t, y, f, mode);
    stepper.last_flag
}

pub(crate) fn mriStepInnerStepper_FullRhsSUNStepper(
    stepper: &mut MRIStepInnerStepper,
    t: f64,
    y: &NVector,
    f: &mut NVector,
    ark_mode: i32,
) -> i32 {
    let mode: SUNFullRhsMode = match ark_mode {
        ARK_FULLRHS_START => SUN_FULLRHS_START,
        ARK_FULLRHS_END => SUN_FULLRHS_END,
        _ => SUN_FULLRHS_OTHER,
    };

    let sunstepper = match stepper
        .content
        .as_mut()
        .and_then(|c| c.downcast_mut::<SUNStepper>())
    {
        Some(s) => s,
        None => return crate::arkode_impl::ARK_SUNSTEPPER_ERR,
    };
    let err = SUNStepper_FullRhs(sunstepper, t, y, f, mode);
    let lf = sunstepper.last_flag;
    stepper.last_flag = lf;
    if err != SUN_SUCCESS {
        return crate::arkode_impl::ARK_SUNSTEPPER_ERR;
    }
    ARK_SUCCESS
}

/* Reset the inner (fast) stepper state */
pub(crate) fn mriStepInnerStepper_Reset(
    stepper: &mut MRIStepInnerStepper,
    tR: f64,
    yR: &NVector,
) -> i32 {
    if let Some(reset) = stepper.ops.reset {
        stepper.last_flag = reset(stepper, tR, yR);
        stepper.last_flag
    } else {
        /* assume stepper uses input state and does not need to be reset */
        ARK_SUCCESS
    }
}

/* Gets the inner (fast) stepper accumulated error */
pub(crate) fn mriStepInnerStepper_GetAccumulatedError(
    stepper: &mut MRIStepInnerStepper,
    accum_error: &mut f64,
) -> i32 {
    if let Some(geterror) = stepper.ops.geterror {
        stepper.last_flag = geterror(stepper, accum_error);
        stepper.last_flag
    } else {
        ARK_INNERSTEP_FAIL
    }
}

/* Resets the inner (fast) stepper accumulated error */
pub(crate) fn mriStepInnerStepper_ResetAccumulatedError(stepper: &mut MRIStepInnerStepper) -> i32 {
    if stepper.ops.geterror.is_some() {
        let reseterror = stepper.ops.reseterror.unwrap();
        stepper.last_flag = reseterror(stepper);
        stepper.last_flag
    } else {
        /* assume stepper provides exact solution and needs no reset */
        ARK_SUCCESS
    }
}

/* Sets the inner (fast) stepper relative tolerance scaling factor */
pub(crate) fn mriStepInnerStepper_SetRTol(stepper: &mut MRIStepInnerStepper, rtol: f64) -> i32 {
    if let Some(setrtol) = stepper.ops.setrtol {
        stepper.last_flag = setrtol(stepper, rtol);
        stepper.last_flag
    } else {
        /* assume stepper provides exact solution */
        ARK_SUCCESS
    }
}

pub(crate) fn mriStepInnerStepper_ResetSUNStepper(
    stepper: &mut MRIStepInnerStepper,
    tR: f64,
    yR: &NVector,
) -> i32 {
    let sunstepper = match stepper
        .content
        .as_mut()
        .and_then(|c| c.downcast_mut::<SUNStepper>())
    {
        Some(s) => s,
        None => return crate::arkode_impl::ARK_SUNSTEPPER_ERR,
    };
    let err = SUNStepper_Reset(sunstepper, tR, yR);
    let lf = sunstepper.last_flag;
    stepper.last_flag = lf;
    if err != SUN_SUCCESS {
        return crate::arkode_impl::ARK_SUNSTEPPER_ERR;
    }
    ARK_SUCCESS
}

/* Allocate MRI forcing vectors if necessary (the fused-op vals/vecs
workspace is assembled at call sites) */
pub(crate) fn mriStepInnerStepper_AllocVecs(
    stepper: &mut MRIStepInnerStepper,
    count: i32,
    tmpl_len: usize,
) -> i32 {
    /* Set space requirements for one N_Vector */
    stepper.lrw1 = tmpl_len as i64;
    stepper.liw1 = 1;

    /* Set the number of forcing vectors and allocate vectors */
    stepper.nforcing = count;

    if stepper.nforcing_allocated < stepper.nforcing {
        if stepper.nforcing_allocated > 0 {
            let (lrw1, liw1) = (stepper.lrw1, stepper.liw1);
            let MRIStepInnerStepper { forcing, lrw, liw, nforcing_allocated, .. } = stepper;
            arkFreeVecArray(*nforcing_allocated, forcing, lrw1, lrw, liw1, liw);
        }
        {
            let (lrw1, liw1) = (stepper.lrw1, stepper.liw1);
            let MRIStepInnerStepper { forcing, lrw, liw, nforcing, .. } = stepper;
            if !arkAllocVecArray(*nforcing, tmpl_len, forcing, lrw1, lrw, liw1, liw) {
                mriStepInnerStepper_FreeVecs(stepper);
                return ARK_MEM_FAIL;
            }
        }
        stepper.nforcing_allocated = stepper.nforcing;
    }

    ARK_SUCCESS
}

/* Free MRI forcing vectors if necessary */
pub(crate) fn mriStepInnerStepper_FreeVecs(stepper: &mut MRIStepInnerStepper) -> i32 {
    let (lrw1, liw1) = (stepper.lrw1, stepper.liw1);
    let MRIStepInnerStepper { forcing, lrw, liw, nforcing_allocated, .. } = stepper;
    arkFreeVecArray(*nforcing_allocated, forcing, lrw1, lrw, liw1, liw);
    ARK_SUCCESS
}

/* Print forcing vectors to output file */
pub(crate) fn mriStepInnerStepper_PrintMem(
    stepper: &MRIStepInnerStepper,
    outfile: &mut dyn std::io::Write,
) {
    /* output data from the inner stepper */
    let _ = writeln!(outfile, "MRIStepInnerStepper Mem:");
    let _ = writeln!(
        outfile,
        "MRIStepInnerStepper: inner_nforcing = {}",
        stepper.nforcing
    );
}

/*---------------------------------------------------------------
  Utility routines for MRIStep to serve as an MRIStepInnerStepper
  ---------------------------------------------------------------*/

/*------------------------------------------------------------------------------
  mriStep_ApplyForcing

  Determines the linear combination coefficients to apply forcing at a given
  value of the independent variable (t).  C appends the coefficients and
  N_Vector pointers to the cvals/Xvecs arrays; the Rust port returns the
  scaling values [s, s*tau, s*tau^2, ...] (the forcing vectors are appended
  to the operand list at the call site).
  ----------------------------------------------------------------------------*/
pub(crate) fn mriStep_ApplyForcing_coeffs(
    step_mem: &ARKodeMRIStepMem,
    t: f64,
    s: f64,
) -> Vec<f64> {
    let mut vals: Vec<f64> = Vec::with_capacity(step_mem.nforcing as usize);

    /* always append the constant forcing term */
    vals.push(s);

    /* compute normalized time tau and initialize tau^i */
    let tau = (t - step_mem.tshift) / step_mem.tscale;
    let mut taui = tau;
    for _i in 1..step_mem.nforcing {
        vals.push(s * taui);
        taui *= tau;
    }

    vals
}

/*------------------------------------------------------------------------------
  mriStep_SetInnerForcing

  Sets an array of coefficient vectors for a time-dependent external polynomial
  forcing term in the ODE RHS i.e., y' = f(t,y) + p(t).  Primarily intended
  for using MRIStep as an inner integrator within another [outer] instance of
  MRIStep.  The C code stores the caller's vector-array pointer, the Rust
  port stores owned copies.
  ----------------------------------------------------------------------------*/
fn mriStep_SetInnerForcing(
    ark_mem: &mut ARKodeMem,
    tshift: f64,
    tscale: f64,
    forcing: &[NVector],
    nvecs: i32,
) -> i32 {
    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_SetInnerForcing") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    if nvecs > 0 {
        /* enable forcing, and signal that the corresponding pre-existing RHS
           vector is no longer current, since it has a stale forcing function */
        if step_mem.explicit_rhs {
            step_mem.expforcing = true;
            step_mem.impforcing = false;
            step_mem.fse_is_current = false;
        } else {
            step_mem.expforcing = false;
            step_mem.impforcing = true;
            step_mem.fsi_is_current = false;
        }
        step_mem.tshift = tshift;
        step_mem.tscale = tscale;
        step_mem.forcing = forcing.to_vec();
        step_mem.nforcing = nvecs;

        /* Signal that any pre-existing RHS vector is no longer current, since
           it has a stale forcing function */
        ark_mem.fn_is_current = false;

        /* If the coupling table is NULL, then mriStep_Init has not been called
           and the number of stages has not been set yet.  On subsequent calls
           we check if enough space has already been allocated in case nforcing
           has increased since the original allocation. */
        if let Some(mric_stages) = step_mem.MRIC.as_ref().map(|m| m.stages) {
            /* check if there are enough reusable arrays for fused operations */
            if (step_mem.nfusedopvecs - nvecs) < (2 * mric_stages + 2) {
                /* free current work space */
                if !step_mem.cvals.is_empty() {
                    step_mem.cvals = Vec::new();
                    ark_mem.lrw -= step_mem.nfusedopvecs as i64;
                    ark_mem.liw -= step_mem.nfusedopvecs as i64;
                }

                /* allocate reusable arrays for fused vector operations */
                step_mem.nfusedopvecs = 2 * mric_stages + 2 + nvecs;
                step_mem.cvals = vec![ZERO; step_mem.nfusedopvecs as usize];
                ark_mem.lrw += step_mem.nfusedopvecs as i64;
                ark_mem.liw += step_mem.nfusedopvecs as i64;
            }
        }
    } else {
        /* disable forcing */
        step_mem.expforcing = false;
        step_mem.impforcing = false;
        step_mem.tshift = ZERO;
        step_mem.tscale = ONE;
        step_mem.forcing = Vec::new();
        step_mem.nforcing = 0;
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}
