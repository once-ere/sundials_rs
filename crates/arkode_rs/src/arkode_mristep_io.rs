/* -----------------------------------------------------------------
 * Translation of src/arkode/arkode_mristep_io.c (SUNDIALS 7.7.0).
 *
 * This is the implementation file for the optional input and
 * output functions for the ARKODE MRIStep time stepper module.
 *
 * The deprecated MRIStep* wrapper aliases at the end of the C file
 * are omitted (same convention as the other steppers).
 * -----------------------------------------------------------------*/

use crate::arkode_impl::{
    arkProcessError, ARKStagePredictFn, ARKodeMem, ARK_ILL_INPUT, ARK_MEM_NULL, ARK_NLS_OP_ERR,
    ARK_NO_FAILURES, ARK_STEPPER_UNSUPPORTED, ARK_SUCCESS, ZERO,
};
use crate::arkode_io::{sunfprintf_long, sunfprintf_real};
use crate::arkode_mri_tables::{
    MRIStepCouplingMem, MRIStepCoupling_Copy, MRIStepCoupling_Free,
    MRIStepCoupling_LoadTableByName, MRIStepCoupling_Space,
};
use crate::arkode_mristep::mriStep_AccessStepMem;
use crate::arkode_mristep_impl::{
    MRIStepPostInnerFn, MRIStepPreInnerFn, CRDOWN, DGMAX, MAXCOR, MSBP, MSG_MRISTEP_NO_COUPLING,
    NLSCOEF, RDIV,
};
use crate::nvector_serial::{NVector, N_VScale};
use crate::sundials_adaptcontroller::{
    SUNAdaptController, SUNAdaptController_GetType, SUN_ADAPTCONTROLLER_MRI_H_TOL,
};
use crate::sundials_types::{SUNOutputFormat, SUN_UNIT_ROUNDOFF};

const ONE: f64 = 1.0;

/*===============================================================
  Exported optional input functions.
  ===============================================================*/

/*---------------------------------------------------------------
  MRIStepSetCoupling:

  Specifies to use a customized coupling structure for the slow
  portion of the system.
  ---------------------------------------------------------------*/
pub fn MRIStepSetCoupling(ark_mem: &mut ARKodeMem, MRIC: &MRIStepCouplingMem) -> i32 {
    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "MRIStepSetCoupling") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* C checks MRIC == NULL (inexpressible here) */

    /* clear any existing parameters and coupling structure */
    step_mem.stages = 0;
    step_mem.q = 0;
    step_mem.p = 0;
    let mut Tliw: i64 = 0;
    let mut Tlrw: i64 = 0;
    MRIStepCoupling_Space(step_mem.MRIC.as_deref(), &mut Tliw, &mut Tlrw);
    MRIStepCoupling_Free(&mut step_mem.MRIC);
    ark_mem.liw -= Tliw;
    ark_mem.lrw -= Tlrw;

    /* set the relevant parameters */
    step_mem.stages = MRIC.stages;
    step_mem.q = MRIC.q;
    step_mem.p = MRIC.p;

    /* copy the coupling structure in step memory */
    step_mem.MRIC = MRIStepCoupling_Copy(MRIC);
    if step_mem.MRIC.is_none() {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_NULL,
            line!(),
            "MRIStepSetCoupling",
            file!(),
            MSG_MRISTEP_NO_COUPLING,
        );
        return ARK_MEM_NULL;
    }
    MRIStepCoupling_Space(step_mem.MRIC.as_deref(), &mut Tliw, &mut Tlrw);
    ark_mem.liw += Tliw;
    ark_mem.lrw += Tlrw;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  MRIStepSetPreInnerFn:

  Sets the user-supplied function called BEFORE the inner evolve
  ---------------------------------------------------------------*/
pub fn MRIStepSetPreInnerFn(ark_mem: &mut ARKodeMem, prefn: Option<MRIStepPreInnerFn>) -> i32 {
    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "MRIStepSetPreInnerFn") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* Set pre inner evolve function */
    step_mem.pre_inner_evolve = prefn;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  MRIStepSetPostInnerFn:

  Sets the user-supplied function called AFTER the inner evolve
  ---------------------------------------------------------------*/
pub fn MRIStepSetPostInnerFn(ark_mem: &mut ARKodeMem, postfn: Option<MRIStepPostInnerFn>) -> i32 {
    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "MRIStepSetPostInnerFn") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* Set pre inner evolve function */
    step_mem.post_inner_evolve = postfn;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*===============================================================
  Exported optional output functions.
  ===============================================================*/

/*---------------------------------------------------------------
  mriStep_GetNumRhsEvals:

  Returns the current number of RHS calls
  ---------------------------------------------------------------*/
pub fn mriStep_GetNumRhsEvals(
    ark_mem: &mut ARKodeMem,
    partition_index: i32,
    rhs_evals: &mut i64,
) -> i32 {
    let step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_GetNumRhsEvals") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    if partition_index > 1 {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "mriStep_GetNumRhsEvals",
            file!(),
            "Invalid partition index",
        );
        return ARK_ILL_INPUT;
    }

    match partition_index {
        0 => *rhs_evals = step_mem.nfse,
        1 => *rhs_evals = step_mem.nfsi,
        _ => *rhs_evals = step_mem.nfse + step_mem.nfsi,
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

pub fn MRIStepGetNumRhsEvals(
    ark_mem: &mut ARKodeMem,
    nfse_evals: &mut i64,
    nfsi_evals: &mut i64,
) -> i32 {
    let retval = crate::arkode_io::ARKodeGetNumRhsEvals(ark_mem, 0, nfse_evals);
    if retval != ARK_SUCCESS {
        return retval;
    }

    let retval = crate::arkode_io::ARKodeGetNumRhsEvals(ark_mem, 1, nfsi_evals);
    if retval != ARK_SUCCESS {
        return retval;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  MRIStepGetCurrentCoupling:

  Sets pointer to the slow coupling structure currently in use.
  (C hands out the internal pointer; the Rust port returns a copy.)
  ---------------------------------------------------------------*/
pub fn MRIStepGetCurrentCoupling(
    ark_mem: &mut ARKodeMem,
    MRIC: &mut Option<crate::arkode_mri_tables::MRIStepCoupling>,
) -> i32 {
    let step_mem = match mriStep_AccessStepMem(ark_mem, "MRIStepGetCurrentCoupling") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* get coupling structure from step_mem */
    *MRIC = step_mem.MRIC.as_deref().and_then(MRIStepCoupling_Copy);

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  MRIStepGetLastInnerStepFlag:

  Returns the last return value from the inner stepper.
  ---------------------------------------------------------------*/
pub fn MRIStepGetLastInnerStepFlag(ark_mem: &mut ARKodeMem, flag: &mut i32) -> i32 {
    let step_mem = match mriStep_AccessStepMem(ark_mem, "MRIStepGetLastInnerStepFlag") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* get the last return value from the inner stepper */
    *flag = step_mem.stepper.last_flag;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  MRIStepGetNumInnerStepperFails:

  Returns the number of recoverable failures encountered by the
  inner stepper.
  ---------------------------------------------------------------*/
pub fn MRIStepGetNumInnerStepperFails(ark_mem: &mut ARKodeMem, inner_fails: &mut i64) -> i32 {
    let step_mem = match mriStep_AccessStepMem(ark_mem, "MRIStepGetNumInnerStepperFails") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* set output from step_mem */
    *inner_fails = step_mem.inner_fails;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*===============================================================
  Private functions attached to ARKODE
  ===============================================================*/

/*---------------------------------------------------------------
  mriStep_SetOptions:

  Provides command-line control over MRIStep-specific "set"
  routines.
  ---------------------------------------------------------------*/
pub fn mriStep_SetOptions(
    ark_mem: &mut ARKodeMem,
    argidx: &mut usize,
    argv: &[String],
    offset: usize,
    arg_used: &mut bool,
) -> i32 {
    /* The only MRIStep-specific "Set" routine takes a custom MRIStepCoupling
       table; however, these may be specified by name. */
    if &argv[*argidx][offset..] == "coupling_table_name" {
        *argidx += 1;
        let mut Coupling = MRIStepCoupling_LoadTableByName(&argv[*argidx]);
        if Coupling.is_none() {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "mriStep_SetOptions",
                file!(),
                &format!(
                    "error setting key {} {} (invalid table name)",
                    argv[*argidx - 1],
                    argv[*argidx]
                ),
            );
            return ARK_ILL_INPUT;
        }
        let retval = MRIStepSetCoupling(ark_mem, Coupling.as_deref().unwrap());
        MRIStepCoupling_Free(&mut Coupling);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                retval,
                line!(),
                "mriStep_SetOptions",
                file!(),
                &format!(
                    "error setting key {} {} (SetCoupling failed)",
                    argv[*argidx - 1],
                    argv[*argidx]
                ),
            );
            return retval;
        }
        *arg_used = true;
        return ARK_SUCCESS;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_SetAdaptController:

  Specifies a temporal adaptivity controller for MRIStep to use.
  If a non-MRI controller is provided, this just passes that
  through to arkReplaceAdaptController.  The MRI-H-TOL wrapper
  layer (arkode_mristep_controller.c) is deferred -- see the
  module header note in arkode_mristep.rs.
  ---------------------------------------------------------------*/
pub fn mriStep_SetAdaptController(ark_mem: &mut ARKodeMem, C: Option<SUNAdaptController>) -> i32 {
    /* Retrieve the controller type */
    let ctype = match C.as_ref() {
        Some(c) => SUNAdaptController_GetType(c),
        None => crate::sundials_adaptcontroller::SUN_ADAPTCONTROLLER_NONE,
    };

    /* If this does not have MRI type, then just pass to ARKODE */
    if ctype != SUN_ADAPTCONTROLLER_MRI_H_TOL {
        return crate::arkode_io::arkReplaceAdaptController(ark_mem, C, false);
    }

    /* (deferred) SUNAdaptController_MRIStep wrapper */
    arkProcessError(
        Some(ark_mem),
        ARK_ILL_INPUT,
        line!(),
        "mriStep_SetAdaptController",
        file!(),
        "MRI-H-TOL controller wrapper is not yet supported in this port",
    );
    ARK_ILL_INPUT
}

/*---------------------------------------------------------------
  mriStep_SetUserData:

  Passes user-data pointer to attached linear solver module.
  ---------------------------------------------------------------*/
pub fn mriStep_SetUserData(ark_mem: &mut ARKodeMem) -> i32 {
    /* set user data in ARKODELS mem */
    if ark_mem.lmem.is_some() {
        let retval = crate::arkode_ls::arkLSSetUserData(ark_mem);
        if retval != 0 {
            return retval;
        }
    }
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_SetDefaults:

  Resets all MRIStep optional inputs to their default values.
  Does not change problem-defining function pointers or
  user_data pointer.
  ---------------------------------------------------------------*/
pub fn mriStep_SetDefaults(ark_mem: &mut ARKodeMem) -> i32 {
    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_SetDefaults") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* Set default values for integrator optional inputs */
    step_mem.q = 3; /* method order */
    step_mem.p = 0; /* embedding order */
    step_mem.predictor = 0; /* trivial predictor */
    step_mem.linear = false; /* nonlinear problem */
    step_mem.linear_timedep = true; /* dfs/dy depends on t */
    step_mem.deduce_rhs = false; /* deduce fi on result of NLS */
    step_mem.maxcor = MAXCOR; /* max nonlinear iters/stage */
    step_mem.nlscoef = NLSCOEF; /* nonlinear tolerance coefficient */
    step_mem.crdown = CRDOWN; /* nonlinear convergence estimate coeff. */
    step_mem.rdiv = RDIV; /* nonlinear divergence tolerance */
    step_mem.dgmax = DGMAX; /* max gamma change to recompute J or P */
    step_mem.msbp = MSBP; /* max steps between updating J or P */
    step_mem.stages = 0; /* no stages */
    step_mem.istage = 0; /* implicit solver stage index */
    step_mem.cur_stage = 0; /* current stage index */
    step_mem.jcur = false;
    step_mem.convfail = ARK_NO_FAILURES;
    step_mem.stage_predict = None; /* no user-supplied stage predictor */

    /* Remove pre-existing nonlinear solver object */
    step_mem.NLS = None;

    /* Remove pre-existing coupling table */
    if step_mem.MRIC.is_some() {
        let mut Cleniw: i64 = 0;
        let mut Clenrw: i64 = 0;
        MRIStepCoupling_Space(step_mem.MRIC.as_deref(), &mut Cleniw, &mut Clenrw);
        ark_mem.lrw -= Clenrw;
        ark_mem.liw -= Cleniw;
        MRIStepCoupling_Free(&mut step_mem.MRIC);
    }
    step_mem.MRIC = None;

    ark_mem.step_mem = Some(step_mem);

    /* Load the default SUNAdaptController */
    let retval = crate::arkode_io::arkReplaceAdaptController(ark_mem, None, true);
    if retval != 0 {
        return retval;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_SetLinear:

  Specifies that the implicit slow function, fsi(t,y), is linear
  in y, and to tighten the linear solver tolerances while taking
  only one Newton iteration.
  ---------------------------------------------------------------*/
pub fn mriStep_SetLinear(ark_mem: &mut ARKodeMem, timedepend: i32) -> i32 {
    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_SetLinear") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* set parameters */
    step_mem.linear = true;
    step_mem.linear_timedep = timedepend == 1;
    step_mem.dgmax = 100.0 * SUN_UNIT_ROUNDOFF;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_SetNonlinear:

  Specifies that the implicit slow function, fsi(t,y), is
  nonlinear in y.  Used to undo a previous call to
  mriStep_SetLinear.
  ---------------------------------------------------------------*/
pub fn mriStep_SetNonlinear(ark_mem: &mut ARKodeMem) -> i32 {
    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_SetNonlinear") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* set parameters */
    step_mem.linear = false;
    step_mem.linear_timedep = true;
    step_mem.dgmax = DGMAX;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_SetOrder:

  Specifies the method order
  ---------------------------------------------------------------*/
pub fn mriStep_SetOrder(ark_mem: &mut ARKodeMem, ord: i32) -> i32 {
    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_SetOrder") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* check for illegal inputs */
    if ord <= 0 {
        step_mem.q = 3;
    } else {
        step_mem.q = ord;
    }

    /* Clear tables, the user is requesting a change in method or a reset to
       defaults. Tables will be set in InitialSetup. */
    step_mem.stages = 0;
    step_mem.p = 0;
    let mut Tliw: i64 = 0;
    let mut Tlrw: i64 = 0;
    MRIStepCoupling_Space(step_mem.MRIC.as_deref(), &mut Tliw, &mut Tlrw);
    MRIStepCoupling_Free(&mut step_mem.MRIC);
    ark_mem.liw -= Tliw;
    ark_mem.lrw -= Tlrw;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_SetNonlinCRDown:

  Specifies the user-provided nonlinear convergence constant
  crdown.  Legal values are strictly positive; illegal values
  imply a reset to the default.
  ---------------------------------------------------------------*/
pub fn mriStep_SetNonlinCRDown(ark_mem: &mut ARKodeMem, crdown: f64) -> i32 {
    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_SetNonlinCRDown") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* if argument legal set it, otherwise set default */
    if crdown <= ZERO {
        step_mem.crdown = CRDOWN;
    } else {
        step_mem.crdown = crdown;
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_SetNonlinRDiv:

  Specifies the user-provided nonlinear convergence constant
  rdiv.  Legal values are strictly positive; illegal values
  imply a reset to the default.
  ---------------------------------------------------------------*/
pub fn mriStep_SetNonlinRDiv(ark_mem: &mut ARKodeMem, rdiv: f64) -> i32 {
    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_SetNonlinRDiv") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* if argument legal set it, otherwise set default */
    if rdiv <= ZERO {
        step_mem.rdiv = RDIV;
    } else {
        step_mem.rdiv = rdiv;
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_SetDeltaGammaMax:

  Specifies the user-provided linear setup decision constant
  dgmax.  Legal values are strictly positive; illegal values imply
  a reset to the default.
  ---------------------------------------------------------------*/
pub fn mriStep_SetDeltaGammaMax(ark_mem: &mut ARKodeMem, dgmax: f64) -> i32 {
    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_SetDeltaGammaMax") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* if argument legal set it, otherwise set default */
    if dgmax <= ZERO {
        step_mem.dgmax = DGMAX;
    } else {
        step_mem.dgmax = dgmax;
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_SetLSetupFrequency:

  Specifies the user-provided linear setup decision constant
  msbp.
  ---------------------------------------------------------------*/
pub fn mriStep_SetLSetupFrequency(ark_mem: &mut ARKodeMem, msbp: i32) -> i32 {
    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_SetLSetupFrequency") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* if argument legal set it, otherwise set default */
    if msbp == 0 {
        step_mem.msbp = MSBP;
    } else {
        step_mem.msbp = msbp;
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_SetPredictorMethod:

  Specifies the method to use for predicting implicit solutions.
  Non-default choices are {1,2,3,4}, all others will use default
  (trivial) predictor.
  ---------------------------------------------------------------*/
pub fn mriStep_SetPredictorMethod(ark_mem: &mut ARKodeMem, pred_method: i32) -> i32 {
    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_SetPredictorMethod") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* set parameter */
    step_mem.predictor = pred_method;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_SetMaxNonlinIters:

  Specifies the maximum number of nonlinear iterations during
  one solve.  A non-positive input implies a reset to the
  default value.
  ---------------------------------------------------------------*/
pub fn mriStep_SetMaxNonlinIters(ark_mem: &mut ARKodeMem, maxcor: i32) -> i32 {
    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_SetMaxNonlinIters") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* Return error message if no NLS module is present */
    if step_mem.NLS.is_none() {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_NLS_OP_ERR,
            line!(),
            "mriStep_SetMaxNonlinIters",
            file!(),
            "No SUNNonlinearSolver object is present",
        );
        return ARK_ILL_INPUT;
    }

    /* argument <= 0 sets default, otherwise set input */
    if maxcor <= 0 {
        step_mem.maxcor = MAXCOR;
    } else {
        step_mem.maxcor = maxcor;
    }

    /* send argument to NLS structure */
    let maxcor_val = step_mem.maxcor;
    step_mem.NLS.as_mut().unwrap().set_max_iters(maxcor_val);

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_SetNonlinConvCoef:

  Specifies the coefficient in the nonlinear solver convergence
  test.  A non-positive input implies a reset to the default value.
  ---------------------------------------------------------------*/
pub fn mriStep_SetNonlinConvCoef(ark_mem: &mut ARKodeMem, nlscoef: f64) -> i32 {
    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_SetNonlinConvCoef") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* argument <= 0 sets default, otherwise set input */
    if nlscoef <= ZERO {
        step_mem.nlscoef = NLSCOEF;
    } else {
        step_mem.nlscoef = nlscoef;
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_SetStagePredictFn:  Specifies a user-provided step
  predictor function having type ARKStagePredictFn.  A
  NULL input function disables calls to this routine.
  ---------------------------------------------------------------*/
pub fn mriStep_SetStagePredictFn(
    ark_mem: &mut ARKodeMem,
    PredictStage: Option<ARKStagePredictFn>,
) -> i32 {
    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_SetStagePredictFn") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    step_mem.stage_predict = PredictStage;
    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_SetDeduceImplicitRhs:

  Specifies if an optimization is used to avoid an evaluation of
  fi after a nonlinear solve for an implicit stage.
  ---------------------------------------------------------------*/
pub fn mriStep_SetDeduceImplicitRhs(ark_mem: &mut ARKodeMem, deduce: bool) -> i32 {
    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_SetDeduceImplicitRhs") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    step_mem.deduce_rhs = deduce;
    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_GetCurrentGamma: Returns the current value of gamma
  ---------------------------------------------------------------*/
pub fn mriStep_GetCurrentGamma(ark_mem: &mut ARKodeMem, gamma: &mut f64) -> i32 {
    let step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_GetCurrentGamma") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };
    *gamma = step_mem.gamma;
    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_GetEstLocalErrors: Returns the current local truncation
  error estimate vector
  ---------------------------------------------------------------*/
pub fn mriStep_GetEstLocalErrors(ark_mem: &mut ARKodeMem, ele: &mut NVector) -> i32 {
    use crate::arkode_impl::ARK_ACCUMERROR_NONE;
    let step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_GetEstLocalErrors") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* return an error if local truncation error is not computed */
    if (ark_mem.fixedstep && ark_mem.AccumErrorType == ARK_ACCUMERROR_NONE) || step_mem.p <= 0 {
        ark_mem.step_mem = Some(step_mem);
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* otherwise, copy local truncation error vector to output */
    N_VScale(ONE, &ark_mem.tempv1, ele);
    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_GetNumLinSolvSetups:

  Returns the current number of calls to the lsetup routine
  ---------------------------------------------------------------*/
pub fn mriStep_GetNumLinSolvSetups(ark_mem: &mut ARKodeMem, nlinsetups: &mut i64) -> i32 {
    let step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_GetNumLinSolvSetups") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* get value from step_mem */
    *nlinsetups = step_mem.nsetups;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_GetNumNonlinSolvIters:

  Returns the current number of nonlinear solver iterations
  ---------------------------------------------------------------*/
pub fn mriStep_GetNumNonlinSolvIters(ark_mem: &mut ARKodeMem, nniters: &mut i64) -> i32 {
    let step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_GetNumNonlinSolvIters") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    *nniters = step_mem.nls_iters;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_GetNumNonlinSolvConvFails:

  Returns the current number of nonlinear solver convergence fails
  ---------------------------------------------------------------*/
pub fn mriStep_GetNumNonlinSolvConvFails(ark_mem: &mut ARKodeMem, nnfails: &mut i64) -> i32 {
    let step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_GetNumNonlinSolvConvFails") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* set output from step_mem */
    *nnfails = step_mem.nls_fails;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_GetNonlinSolvStats:

  Returns nonlinear solver statistics
  ---------------------------------------------------------------*/
pub fn mriStep_GetNonlinSolvStats(
    ark_mem: &mut ARKodeMem,
    nniters: &mut i64,
    nnfails: &mut i64,
) -> i32 {
    let step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_GetNonlinSolvStats") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    *nniters = step_mem.nls_iters;
    *nnfails = step_mem.nls_fails;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_GetStageIndex:

  Returns the current stage index and number of stages
  ---------------------------------------------------------------*/
pub fn mriStep_GetStageIndex(ark_mem: &mut ARKodeMem, stage: &mut i32, max_stages: &mut i32) -> i32 {
    let step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_GetStageIndex") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    *stage = step_mem.cur_stage;
    *max_stages = step_mem.stages;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_PrintAllStats:

  Prints integrator statistics
  ---------------------------------------------------------------*/
pub fn mriStep_PrintAllStats(
    ark_mem: &mut ARKodeMem,
    outfile: &mut dyn std::io::Write,
    fmt: SUNOutputFormat,
) -> i32 {
    let step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_PrintAllStats") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* function evaluations */
    sunfprintf_long(outfile, fmt, true, "Explicit slow RHS fn evals", step_mem.nfse);
    sunfprintf_long(outfile, fmt, false, "Implicit slow RHS fn evals", step_mem.nfsi);

    /* inner stepper and nonlinear solver stats */
    sunfprintf_long(outfile, fmt, false, "Inner stepper failures", step_mem.inner_fails);
    sunfprintf_long(outfile, fmt, false, "NLS iters", step_mem.nls_iters);
    sunfprintf_long(outfile, fmt, false, "NLS fails", step_mem.nls_fails);
    if ark_mem.nst > 0 {
        sunfprintf_real(
            outfile,
            fmt,
            false,
            "NLS iters per step",
            step_mem.nls_iters as f64 / ark_mem.nst as f64,
        );
    }

    /* linear solver stats */
    sunfprintf_long(outfile, fmt, false, "LS setups", step_mem.nsetups);
    if let Some(arkls_mem) = ark_mem.lmem.take() {
        sunfprintf_long(outfile, fmt, false, "Jac fn evals", arkls_mem.nje);
        sunfprintf_long(outfile, fmt, false, "LS RHS fn evals", arkls_mem.nfeDQ);
        sunfprintf_long(outfile, fmt, false, "Prec setup evals", arkls_mem.npe);
        sunfprintf_long(outfile, fmt, false, "Prec solves", arkls_mem.nps);
        sunfprintf_long(outfile, fmt, false, "LS iters", arkls_mem.nli);
        sunfprintf_long(outfile, fmt, false, "LS fails", arkls_mem.ncfl);
        sunfprintf_long(outfile, fmt, false, "Jac-times setups", arkls_mem.njtsetup);
        sunfprintf_long(outfile, fmt, false, "Jac-times evals", arkls_mem.njtimes);
        if step_mem.nls_iters > 0 {
            sunfprintf_real(
                outfile,
                fmt,
                false,
                "LS iters per NLS iter",
                arkls_mem.nli as f64 / step_mem.nls_iters as f64,
            );
            sunfprintf_real(
                outfile,
                fmt,
                false,
                "Jac evals per NLS iter",
                arkls_mem.nje as f64 / step_mem.nls_iters as f64,
            );
            sunfprintf_real(
                outfile,
                fmt,
                false,
                "Prec evals per NLS iter",
                arkls_mem.npe as f64 / step_mem.nls_iters as f64,
            );
        }
        ark_mem.lmem = Some(arkls_mem);
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_WriteParameters:

  Outputs all solver parameters to the provided file pointer.
  ---------------------------------------------------------------*/
pub fn mriStep_WriteParameters(ark_mem: &mut ARKodeMem, fp: &mut dyn std::io::Write) -> i32 {
    use crate::sundials_utils::fmt_g;
    let step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_WriteParameters") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* print integrator parameters to file */
    let _ = writeln!(fp, "MRIStep time step module parameters:");
    let _ = writeln!(fp, "  Method order {}", step_mem.q);
    if step_mem.linear {
        let _ = write!(fp, "  Linear implicit problem");
        if step_mem.linear_timedep {
            let _ = writeln!(fp, " (time-dependent Jacobian)");
        } else {
            let _ = writeln!(fp, " (time-independent Jacobian)");
        }
    }
    if step_mem.explicit_rhs && step_mem.implicit_rhs {
        let _ = writeln!(fp, "  ImEx slow time scale");
    } else if step_mem.implicit_rhs {
        let _ = writeln!(fp, "  Implicit slow time scale");
    } else {
        let _ = writeln!(fp, "  Explicit slow time scale");
    }

    if step_mem.implicit_rhs {
        let _ = writeln!(fp, "  Implicit predictor method = {}", step_mem.predictor);
        let _ = writeln!(
            fp,
            "  Implicit solver tolerance coefficient = {}",
            fmt_g(step_mem.nlscoef, 0, 6)
        );
        let _ = writeln!(
            fp,
            "  Maximum number of nonlinear corrections = {}",
            step_mem.maxcor
        );
        let _ = writeln!(
            fp,
            "  Nonlinear convergence rate constant = {}",
            fmt_g(step_mem.crdown, 0, 6)
        );
        let _ = writeln!(
            fp,
            "  Nonlinear divergence tolerance = {}",
            fmt_g(step_mem.rdiv, 0, 6)
        );
        let _ = writeln!(
            fp,
            "  Gamma factor LSetup tolerance = {}",
            fmt_g(step_mem.dgmax, 0, 6)
        );
        let _ = writeln!(
            fp,
            "  Number of steps between LSetup calls = {}",
            step_mem.msbp
        );
    }
    let _ = writeln!(fp);

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}
