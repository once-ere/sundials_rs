/* -----------------------------------------------------------------
 * Translated from src/arkode/arkode_arkstep_io.c (ARKODE 7.7.0).
 * ARKStep optional input/output functions.
 *
 * PART I: the stepper-op implementations (arkStep_SetDefaults,
 * the Set and Get families, PrintAllStats, WriteParameters) plus the
 * ARKStep-specific getters used by examples.  Remaining:
 * ARKStepSetTables/SetTableNum/SetTableName, arkStep_SetOptions
 * (CLI), arkStep_SetRelaxFn (relaxation) and the deprecated
 * ARKStep* alias wrappers.
 * -----------------------------------------------------------------*/
use crate::arkode_arkstep::arkStep_AccessStepMem;
use crate::arkode_arkstep_impl::*;
use crate::arkode_butcher::{ARKodeButcherTable, ARKodeButcherTable_Space};
use crate::arkode_impl::*;
use crate::arkode_io::{sunfprintf_long, sunfprintf_real};
use crate::nvector_serial::*;
use crate::sundials_types::*;
use crate::sundials_utils::fmt_g;

/*===============================================================
  Exported ARKStep optional output functions
  ===============================================================*/

pub fn ARKStepGetNumRhsEvals(
    ark_mem: &mut ARKodeMem,
    fe_evals: &mut i64,
    fi_evals: &mut i64,
) -> i32 {
    let retval = crate::arkode_io::ARKodeGetNumRhsEvals(ark_mem, 0, fe_evals);
    if retval != ARK_SUCCESS {
        return retval;
    }

    let retval = crate::arkode_io::ARKodeGetNumRhsEvals(ark_mem, 1, fi_evals);
    if retval != ARK_SUCCESS {
        return retval;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKStepGetCurrentButcherTables:

  Returns copies of the explicit and implicit Butcher tables
  currently in use (C sets pointers).
  ---------------------------------------------------------------*/
pub fn ARKStepGetCurrentButcherTables(
    ark_mem: &mut ARKodeMem,
    Bi: &mut Option<ARKodeButcherTable>,
    Be: &mut Option<ARKodeButcherTable>,
) -> i32 {
    let step_mem = match arkStep_AccessStepMem(ark_mem, "ARKStepGetCurrentButcherTables") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* get tables from step_mem */
    *Bi = step_mem.Bi.as_ref().and_then(crate::arkode_butcher::ARKodeButcherTable_Copy);
    *Be = step_mem.Be.as_ref().and_then(crate::arkode_butcher::ARKodeButcherTable_Copy);
    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKStepGetTimestepperStats:

  Returns integrator statistics
  ---------------------------------------------------------------*/
#[allow(clippy::too_many_arguments)]
pub fn ARKStepGetTimestepperStats(
    ark_mem: &mut ARKodeMem,
    expsteps: &mut i64,
    accsteps: &mut i64,
    step_attempts: &mut i64,
    fe_evals: &mut i64,
    fi_evals: &mut i64,
    nlinsetups: &mut i64,
    netfails: &mut i64,
) -> i32 {
    let step_mem = match arkStep_AccessStepMem(ark_mem, "ARKStepGetTimestepperStats") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* set expsteps and accsteps from adaptivity structure */
    let ha = ark_mem.hadapt_mem.as_ref().unwrap();
    *expsteps = ha.nst_exp;
    *accsteps = ha.nst_acc;

    /* set remaining outputs */
    *step_attempts = ark_mem.nst_attempts;
    *fe_evals = step_mem.nfe;
    *fi_evals = step_mem.nfi;
    *nlinsetups = step_mem.nsetups;
    *netfails = ark_mem.netf;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*===============================================================
  Private functions attached to ARKODE
  ===============================================================*/

/*---------------------------------------------------------------
  arkStep_GetNumRhsEvals:

  Returns the current number of RHS calls
  ---------------------------------------------------------------*/
pub fn arkStep_GetNumRhsEvals(
    ark_mem: &mut ARKodeMem,
    partition_index: i32,
    rhs_evals: &mut i64,
) -> i32 {
    let step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_GetNumRhsEvals") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    if partition_index > 1 {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "arkStep_GetNumRhsEvals",
            file!(),
            "Invalid partition index",
        );
        return ARK_ILL_INPUT;
    }

    match partition_index {
        0 => *rhs_evals = step_mem.nfe,
        1 => *rhs_evals = step_mem.nfi,
        _ => *rhs_evals = step_mem.nfe + step_mem.nfi,
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKStepSetAdaptivityMethod: user-callable deprecated wrapper
  around arkSetAdaptivityMethod.
  ---------------------------------------------------------------*/
pub fn ARKStepSetAdaptivityMethod(
    ark_mem: &mut ARKodeMem,
    imethod: i32,
    idefault: i32,
    pq: i32,
    adapt_params: Option<&[f64; 3]>,
) -> i32 {
    crate::arkode_io::arkSetAdaptivityMethod(ark_mem, imethod, idefault, pq, adapt_params)
}

/*---------------------------------------------------------------
  ARKStepSetAdaptivityFn: user-callable deprecated wrapper around
  arkSetAdaptivityFn.
  ---------------------------------------------------------------*/
pub fn ARKStepSetAdaptivityFn(
    ark_mem: &mut ARKodeMem,
    hfun: Option<crate::arkode_impl::ARKAdaptFn>,
    h_data: crate::sundials_types::UserData,
) -> i32 {
    crate::arkode_io::arkSetAdaptivityFn(ark_mem, hfun, h_data)
}

/*---------------------------------------------------------------
  arkStep_SetRelaxFn:

  Sets up the relaxation module using ARKStep's utility routines.
  ---------------------------------------------------------------*/
pub fn arkStep_SetRelaxFn(
    ark_mem: &mut ARKodeMem,
    rfn: Option<crate::arkode_impl::ARKRelaxFn>,
    rjac: Option<crate::arkode_impl::ARKRelaxJacFn>,
) -> i32 {
    crate::arkode_relaxation::arkRelaxCreate(
        ark_mem,
        rfn,
        rjac,
        Some(crate::arkode_arkstep::arkStep_RelaxDeltaE),
        Some(crate::arkode_arkstep::arkStep_GetOrder),
    )
}

/*---------------------------------------------------------------
  arkStep_SetUserData: passed through; the Rust steppers and the
  ARKLS interface read ark_mem.user_data directly (the C version
  re-points lmem/mass user data pointers).
  ---------------------------------------------------------------*/
pub fn arkStep_SetUserData(ark_mem: &mut ARKodeMem) -> i32 {
    /* set user data in ARKODE LS mem */
    if ark_mem.lmem.is_some() {
        let retval = crate::arkode_ls::arkLSSetUserData(ark_mem);
        if retval != 0 {
            return retval;
        }
    }
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_SetDefaults:

  Resets all ARKStep optional inputs to their default values.
  Does not change problem-defining function pointers or
  user_data pointer.  Also leaves alone any data
  structures/options related to the ARKODE infrastructure itself
  (e.g., root-finding and post-process step).
  ---------------------------------------------------------------*/
pub fn arkStep_SetDefaults(ark_mem: &mut ARKodeMem) -> i32 {
    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_SetDefaults") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* Set default values for integrator optional inputs */
    step_mem.q = Q_DEFAULT; /* method order */
    step_mem.p = 0; /* embedding order */
    step_mem.predictor = 0; /* trivial predictor */
    step_mem.linear = SUNFALSE; /* nonlinear problem */
    step_mem.linear_timedep = SUNTRUE; /* dfi/dy depends on t */
    step_mem.autonomous = SUNFALSE; /* non-autonomous problem */
    step_mem.explicit = SUNTRUE; /* fe(t,y) will be used */
    step_mem.implicit = SUNTRUE; /* fi(t,y) will be used */
    step_mem.deduce_rhs = SUNFALSE; /* deduce fi on result of NLS */
    step_mem.maxcor = MAXCOR; /* max nonlinear iters/stage */
    step_mem.nlscoef = NLSCOEF; /* nonlinear tolerance coefficient */
    step_mem.crdown = CRDOWN; /* nonlinear convergence estimate coeff. */
    step_mem.rdiv = RDIV; /* nonlinear divergence tolerance */
    step_mem.dgmax = DGMAX; /* max step change before recomputing J or P */
    step_mem.msbp = MSBP; /* max steps between updates to J or P */
    step_mem.stages = 0; /* no stages */
    step_mem.istage = 0; /* current stage */
    step_mem.jcur = SUNFALSE;
    step_mem.convfail = ARK_NO_FAILURES;
    step_mem.stage_predict = None; /* no user-supplied stage predictor */

    /* Remove pre-existing Butcher tables */
    if let Some(bt) = step_mem.Be.take() {
        let (mut bliw, mut blrw) = (0i64, 0i64);
        ARKodeButcherTable_Space(&bt, &mut bliw, &mut blrw);
        ark_mem.liw -= bliw;
        ark_mem.lrw -= blrw;
    }
    if let Some(bt) = step_mem.Bi.take() {
        let (mut bliw, mut blrw) = (0i64, 0i64);
        ARKodeButcherTable_Space(&bt, &mut bliw, &mut blrw);
        ark_mem.liw -= bliw;
        ark_mem.lrw -= blrw;
    }

    /* Remove pre-existing nonlinear solver object */
    step_mem.NLS = None;

    ark_mem.step_mem = Some(step_mem);

    /* Load the default SUNAdaptController */
    let retval = crate::arkode_io::arkReplaceAdaptController(ark_mem, None, SUNTRUE);
    if retval != 0 {
        return retval;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_SetOrder:

  Specifies the method order
  ---------------------------------------------------------------*/
pub fn arkStep_SetOrder(ark_mem: &mut ARKodeMem, ord: i32) -> i32 {
    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_SetOrder") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* set user-provided value, or default, depending on argument */
    if ord <= 0 {
        step_mem.q = Q_DEFAULT;
    } else {
        step_mem.q = ord;
    }

    /* clear Butcher tables, since user is requesting a change in method
       or a reset to defaults.  Tables will be set in ARKInitialSetup. */
    step_mem.stages = 0;
    step_mem.istage = 0;
    step_mem.p = 0;

    if let Some(bt) = step_mem.Be.take() {
        let (mut bliw, mut blrw) = (0i64, 0i64);
        ARKodeButcherTable_Space(&bt, &mut bliw, &mut blrw);
        ark_mem.liw -= bliw;
        ark_mem.lrw -= blrw;
    }
    if let Some(bt) = step_mem.Bi.take() {
        let (mut bliw, mut blrw) = (0i64, 0i64);
        ARKodeButcherTable_Space(&bt, &mut bliw, &mut blrw);
        ark_mem.liw -= bliw;
        ark_mem.lrw -= blrw;
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_SetLinear:

  Specifies that the implicit portion of the problem is linear,
  and to tighten the linear solver tolerances while taking only
  one Newton iteration.  DO NOT USE IN COMBINATION WITH THE
  FIXED-POINT SOLVER.  Automatically tightens DeltaGammaMax
  to ensure that step size changes cause Jacobian recomputation.

  The argument should be 1 or 0, where 1 indicates that the
  Jacobian of fi with respect to y depends on time, and
  0 indicates that it is not time dependent.  Alternately, when
  using an iterative linear solver this flag denotes time
  dependence of the preconditioner.
  ---------------------------------------------------------------*/
pub fn arkStep_SetLinear(ark_mem: &mut ARKodeMem, timedepend: i32) -> i32 {
    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_SetLinear") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    if timedepend != 0 && step_mem.autonomous {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "arkStep_SetLinear",
            file!(),
            "Incompatible settings, the problem is autonomous but the Jacobian is time dependent",
        );
        return ARK_ILL_INPUT;
    }

    /* set parameters */
    step_mem.linear = SUNTRUE;
    step_mem.linear_timedep = timedepend == 1;
    step_mem.dgmax = 100.0 * SUN_UNIT_ROUNDOFF;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_SetNonlinear:

  Specifies that the implicit portion of the problem is nonlinear.
  Used to undo a previous call to arkStep_SetLinear.  Automatically
  loosens DeltaGammaMax back to default value.
  ---------------------------------------------------------------*/
pub fn arkStep_SetNonlinear(ark_mem: &mut ARKodeMem) -> i32 {
    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_SetNonlinear") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* set parameters */
    step_mem.linear = SUNFALSE;
    step_mem.linear_timedep = SUNTRUE;
    step_mem.dgmax = DGMAX;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_SetAutonomous:

  Indicates if the problem is autonomous (True) or non-autonomous
  (False).
  ---------------------------------------------------------------*/
pub fn arkStep_SetAutonomous(ark_mem: &mut ARKodeMem, autonomous: bool) -> i32 {
    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_SetAutonomous") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    step_mem.autonomous = autonomous;

    if autonomous && step_mem.linear {
        step_mem.linear_timedep = SUNFALSE;
    }

    ark_mem.step_mem = Some(step_mem);

    /* Reattach the nonlinear system function e.g., switching to/from an
       autonomous problem with the trivial predictor requires swapping the
       nonlinear system function provided to the nonlinear solver */
    let retval = crate::arkode_arkstep_nls::arkStep_SetNlsSysFn(ark_mem);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "arkStep_SetAutonomous",
            file!(),
            "Setting nonlinear system function failed",
        );
        return ARK_ILL_INPUT;
    }

    /* This will be better handled when the temp vector stack is added */
    if autonomous {
        /* Allocate tempv5 if needed */
        let tmpl_len = ark_mem.yn.data.len();
        let mut tempv5 = std::mem::take(&mut ark_mem.tempv5);
        crate::arkode::arkAllocVec(ark_mem, tmpl_len, &mut tempv5);
        ark_mem.tempv5 = tempv5;
    } else {
        /* Free tempv5 if necessary */
        let mut tempv5 = std::mem::take(&mut ark_mem.tempv5);
        crate::arkode::arkFreeVec(ark_mem, &mut tempv5);
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_SetNonlinCRDown:

  Specifies the user-provided nonlinear convergence constant
  crdown.  Legal values are strictly positive; illegal values
  imply a reset to the default.
  ---------------------------------------------------------------*/
pub fn arkStep_SetNonlinCRDown(ark_mem: &mut ARKodeMem, crdown: f64) -> i32 {
    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_SetNonlinCRDown") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
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
  arkStep_SetNonlinRDiv:

  Specifies the user-provided nonlinear convergence constant
  rdiv.  Legal values are strictly positive; illegal values
  imply a reset to the default.
  ---------------------------------------------------------------*/
pub fn arkStep_SetNonlinRDiv(ark_mem: &mut ARKodeMem, rdiv: f64) -> i32 {
    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_SetNonlinRDiv") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
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
  arkStep_SetDeltaGammaMax:

  Specifies the user-provided linear setup decision constant
  dgmax.  Legal values are strictly positive; illegal values imply
  a reset to the default.
  ---------------------------------------------------------------*/
pub fn arkStep_SetDeltaGammaMax(ark_mem: &mut ARKodeMem, dgmax: f64) -> i32 {
    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_SetDeltaGammaMax") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
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
  arkStep_SetLSetupFrequency:

  Specifies the user-provided linear setup decision constant
  msbp.  Positive values give the frequency for calling lsetup;
  negative values imply recomputation of lsetup at each nonlinear
  solve; a zero value implies a reset to the default.
  ---------------------------------------------------------------*/
pub fn arkStep_SetLSetupFrequency(ark_mem: &mut ARKodeMem, msbp: i32) -> i32 {
    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_SetLSetupFrequency") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
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
  arkStep_SetPredictorMethod:

  Specifies the method to use for predicting implicit solutions.
  Non-default choices are {1,2,3,4}, all others will use default
  (trivial) predictor.
  ---------------------------------------------------------------*/
pub fn arkStep_SetPredictorMethod(ark_mem: &mut ARKodeMem, pred_method: i32) -> i32 {
    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_SetPredictorMethod") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* set parameter */
    step_mem.predictor = pred_method;

    ark_mem.step_mem = Some(step_mem);

    /* Reattach the nonlinear system function e.g., switching to/from the
       trivial predictor with an autonomous problem requires swapping the
       nonlinear system function provided to the nonlinear solver */
    let retval = crate::arkode_arkstep_nls::arkStep_SetNlsSysFn(ark_mem);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "arkStep_SetPredictorMethod",
            file!(),
            "Setting nonlinear system function failed",
        );
        return ARK_ILL_INPUT;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_SetMaxNonlinIters:

  Specifies the maximum number of nonlinear iterations during
  one solve.  A non-positive input implies a reset to the
  default value.
  ---------------------------------------------------------------*/
pub fn arkStep_SetMaxNonlinIters(ark_mem: &mut ARKodeMem, maxcor: i32) -> i32 {
    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_SetMaxNonlinIters") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* Return error message if no NLS module is present */
    if step_mem.NLS.is_none() {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_NLS_OP_ERR,
            line!(),
            "arkStep_SetMaxNonlinIters",
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
    let maxcor = step_mem.maxcor;
    let retval = step_mem.NLS.as_mut().unwrap().set_max_iters(maxcor);
    if retval != 0 {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_NLS_OP_ERR,
            line!(),
            "arkStep_SetMaxNonlinIters",
            file!(),
            "Error setting maxcor in SUNNonlinearSolver object",
        );
        return ARK_NLS_OP_ERR;
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_SetNonlinConvCoef:

  Specifies the coefficient in the nonlinear solver convergence
  test.  A non-positive input implies a reset to the default value.
  ---------------------------------------------------------------*/
pub fn arkStep_SetNonlinConvCoef(ark_mem: &mut ARKodeMem, nlscoef: f64) -> i32 {
    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_SetNonlinConvCoef") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
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
  arkStep_SetStagePredictFn:  Specifies a user-provided step
  predictor function having type ARKStagePredictFn.  A NULL input
  function disables calls to this routine.
  ---------------------------------------------------------------*/
pub fn arkStep_SetStagePredictFn(
    ark_mem: &mut ARKodeMem,
    PredictStage: Option<ARKStagePredictFn>,
) -> i32 {
    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_SetStagePredictFn") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* override predictor method 5 if non-NULL PredictStage is supplied */
    if step_mem.predictor == 5 && PredictStage.is_some() {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "arkStep_SetStagePredictFn",
            file!(),
            "User-supplied predictor is incompatible with predictor method 5",
        );
        return ARK_ILL_INPUT;
    }

    step_mem.stage_predict = PredictStage;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_SetDeduceImplicitRhs:

  Specifies if an optimization is used to avoid an evaluation of
  fi after a nonlinear solve for an implicit stage.  If stage
  postprocessecing in enabled, this option is ignored, and fi is
  never deduced.

  An argument of SUNTRUE indicates that fi is deduced to compute
  fi(z_i), and SUNFALSE indicates that fi(z_i) is computed with
  an additional evaluation of fi.
  ---------------------------------------------------------------*/
pub fn arkStep_SetDeduceImplicitRhs(ark_mem: &mut ARKodeMem, deduce: bool) -> i32 {
    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_SetDeduceImplicitRhs") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* set parameter */
    step_mem.deduce_rhs = deduce;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_GetCurrentGamma: Returns the current value of gamma
  ---------------------------------------------------------------*/
pub fn arkStep_GetCurrentGamma(ark_mem: &mut ARKodeMem, gamma: &mut f64) -> i32 {
    let step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_GetCurrentGamma") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };
    *gamma = step_mem.gamma;
    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_GetEstLocalErrors: Returns the current local truncation
  error estimate vector
  ---------------------------------------------------------------*/
pub fn arkStep_GetEstLocalErrors(ark_mem: &mut ARKodeMem, ele: &mut NVector) -> i32 {
    let step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_GetEstLocalErrors") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* return an error if local truncation error is not computed */
    if (ark_mem.fixedstep && ark_mem.AccumErrorType == ARK_ACCUMERROR_NONE)
        || step_mem.p <= 0
    {
        ark_mem.step_mem = Some(step_mem);
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* otherwise, copy local truncation error vector to output */
    N_VScale(ONE, &ark_mem.tempv1, ele);
    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_GetNumLinSolvSetups:

  Returns the current number of calls to the lsetup routine
  ---------------------------------------------------------------*/
pub fn arkStep_GetNumLinSolvSetups(ark_mem: &mut ARKodeMem, nlinsetups: &mut i64) -> i32 {
    let step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_GetNumLinSolvSetups") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* get value from step_mem */
    *nlinsetups = step_mem.nsetups;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_GetNumNonlinSolvIters:

  Returns the current number of nonlinear solver iterations
  ---------------------------------------------------------------*/
pub fn arkStep_GetNumNonlinSolvIters(ark_mem: &mut ARKodeMem, nniters: &mut i64) -> i32 {
    let step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_GetNumNonlinSolvIters") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    *nniters = step_mem.nls_iters;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_GetNumNonlinSolvConvFails:

  Returns the current number of nonlinear solver convergence fails
  ---------------------------------------------------------------*/
pub fn arkStep_GetNumNonlinSolvConvFails(ark_mem: &mut ARKodeMem, nnfails: &mut i64) -> i32 {
    let step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_GetNumNonlinSolvConvFails") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* set output from step_mem */
    *nnfails = step_mem.nls_fails;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_GetNonlinSolvStats:

  Returns nonlinear solver statistics
  ---------------------------------------------------------------*/
pub fn arkStep_GetNonlinSolvStats(
    ark_mem: &mut ARKodeMem,
    nniters: &mut i64,
    nnfails: &mut i64,
) -> i32 {
    let step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_GetNonlinSolvStats") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    *nniters = step_mem.nls_iters;
    *nnfails = step_mem.nls_fails;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_GetStageIndex:

  Returns the current stage index and number of stages
  ---------------------------------------------------------------*/
pub fn arkStep_GetStageIndex(ark_mem: &mut ARKodeMem, stage: &mut i32, max_stages: &mut i32) -> i32 {
    let step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_GetStageIndex") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    *stage = step_mem.istage;
    *max_stages = step_mem.stages;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_PrintAllStats:

  Prints integrator statistics
  ---------------------------------------------------------------*/
pub fn arkStep_PrintAllStats(
    ark_mem: &mut ARKodeMem,
    outfile: &mut dyn std::io::Write,
    fmt: SUNOutputFormat,
) -> i32 {
    let step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_PrintAllStats") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* function evaluations */
    sunfprintf_long(outfile, fmt, SUNFALSE, "Explicit RHS fn evals", step_mem.nfe);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Implicit RHS fn evals", step_mem.nfi);

    /* nonlinear solver stats */
    sunfprintf_long(outfile, fmt, SUNFALSE, "NLS iters", step_mem.nls_iters);
    sunfprintf_long(outfile, fmt, SUNFALSE, "NLS fails", step_mem.nls_fails);
    if ark_mem.nst > 0 {
        sunfprintf_real(
            outfile,
            fmt,
            SUNFALSE,
            "NLS iters per step",
            step_mem.nls_iters as f64 / ark_mem.nst as f64,
        );
    }

    /* linear solver stats */
    sunfprintf_long(outfile, fmt, SUNFALSE, "LS setups", step_mem.nsetups);
    if let Some(arkls_mem) = ark_mem.lmem.take() {
        sunfprintf_long(outfile, fmt, SUNFALSE, "Jac fn evals", arkls_mem.nje);
        sunfprintf_long(outfile, fmt, SUNFALSE, "LS RHS fn evals", arkls_mem.nfeDQ);
        sunfprintf_long(outfile, fmt, SUNFALSE, "Prec setup evals", arkls_mem.npe);
        sunfprintf_long(outfile, fmt, SUNFALSE, "Prec solves", arkls_mem.nps);
        sunfprintf_long(outfile, fmt, SUNFALSE, "LS iters", arkls_mem.nli);
        sunfprintf_long(outfile, fmt, SUNFALSE, "LS fails", arkls_mem.ncfl);
        sunfprintf_long(outfile, fmt, SUNFALSE, "Jac-times setups", arkls_mem.njtsetup);
        sunfprintf_long(outfile, fmt, SUNFALSE, "Jac-times evals", arkls_mem.njtimes);
        if step_mem.nls_iters > 0 {
            sunfprintf_real(
                outfile,
                fmt,
                SUNFALSE,
                "LS iters per NLS iter",
                arkls_mem.nli as f64 / step_mem.nls_iters as f64,
            );
            sunfprintf_real(
                outfile,
                fmt,
                SUNFALSE,
                "Jac evals per NLS iter",
                arkls_mem.nje as f64 / step_mem.nls_iters as f64,
            );
            sunfprintf_real(
                outfile,
                fmt,
                SUNFALSE,
                "Prec evals per NLS iter",
                arkls_mem.npe as f64 / step_mem.nls_iters as f64,
            );
        }
        ark_mem.lmem = Some(arkls_mem);
    }

    /* (mass solve stats: mass half not ported) */

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_WriteParameters:

  Outputs all solver parameters to the provided file pointer.
  ---------------------------------------------------------------*/
pub fn arkStep_WriteParameters(ark_mem: &mut ARKodeMem, fp: &mut dyn std::io::Write) -> i32 {
    let step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_WriteParameters") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* print integrator parameters to file */
    let _ = write!(fp, "ARKStep time step module parameters:\n");
    let _ = write!(fp, "  Method order {}\n", step_mem.q);
    if step_mem.linear {
        let _ = write!(fp, "  Linear implicit problem");
        if step_mem.linear_timedep {
            let _ = write!(fp, " (time-dependent Jacobian)\n");
        } else {
            let _ = write!(fp, " (time-independent Jacobian)\n");
        }
    }
    if step_mem.explicit && step_mem.implicit {
        let _ = write!(fp, "  ImEx integrator\n");
    } else if step_mem.implicit {
        let _ = write!(fp, "  Implicit integrator\n");
    } else {
        let _ = write!(fp, "  Explicit integrator\n");
    }

    if step_mem.implicit {
        let _ = write!(fp, "  Implicit predictor method = {}\n", step_mem.predictor);
        let _ = write!(
            fp,
            "  Implicit solver tolerance coefficient = {}\n",
            fmt_g(step_mem.nlscoef, 0, 15)
        );
        let _ = write!(
            fp,
            "  Maximum number of nonlinear corrections = {}\n",
            step_mem.maxcor
        );
        let _ = write!(
            fp,
            "  Nonlinear convergence rate constant = {}\n",
            fmt_g(step_mem.crdown, 0, 15)
        );
        let _ = write!(
            fp,
            "  Nonlinear divergence tolerance = {}\n",
            fmt_g(step_mem.rdiv, 0, 15)
        );
        let _ = write!(
            fp,
            "  Gamma factor LSetup tolerance = {}\n",
            fmt_g(step_mem.dgmax, 0, 15)
        );
        let _ = write!(fp, "  Number of steps between LSetup calls = {}\n", step_mem.msbp);
    }
    let _ = write!(fp, "\n");

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKStepSetExplicit:

  Specifies that the implicit portion of the problem is disabled,
  and to use an explicit RK method.
  ---------------------------------------------------------------*/
pub fn ARKStepSetExplicit(ark_mem: &mut ARKodeMem) -> i32 {
    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "ARKStepSetExplicit") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* ensure that fe is defined */
    if step_mem.fe.is_none() {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "ARKStepSetExplicit",
            file!(),
            MSG_ARK_MISSING_FE,
        );
        return ARK_ILL_INPUT;
    }

    /* set the relevant parameters */
    step_mem.explicit = SUNTRUE;
    step_mem.implicit = SUNFALSE;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKStepSetImplicit:

  Specifies that the explicit portion of the problem is disabled,
  and to use an implicit RK method.
  ---------------------------------------------------------------*/
pub fn ARKStepSetImplicit(ark_mem: &mut ARKodeMem) -> i32 {
    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "ARKStepSetImplicit") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* ensure that fi is defined */
    if step_mem.fi.is_none() {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "ARKStepSetImplicit",
            file!(),
            MSG_ARK_MISSING_FI,
        );
        return ARK_ILL_INPUT;
    }

    /* set the relevant parameters */
    step_mem.implicit = SUNTRUE;
    step_mem.explicit = SUNFALSE;

    ark_mem.step_mem = Some(step_mem);

    /* re-attach internal error weight functions if necessary */
    ark_reattach_tolerances(ark_mem)
}

/*---------------------------------------------------------------
  ARKStepSetImEx:

  Specifies that the specifies that problem has both implicit and
  explicit parts, and to use an ARK method (this is the default).
  ---------------------------------------------------------------*/
pub fn ARKStepSetImEx(ark_mem: &mut ARKodeMem) -> i32 {
    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "ARKStepSetImEx") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* ensure that fe and fi are defined */
    if step_mem.fe.is_none() {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "ARKStepSetImEx",
            file!(),
            MSG_ARK_MISSING_FE,
        );
        return ARK_ILL_INPUT;
    }
    if step_mem.fi.is_none() {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "ARKStepSetImEx",
            file!(),
            MSG_ARK_MISSING_FI,
        );
        return ARK_ILL_INPUT;
    }

    /* set the relevant parameters */
    step_mem.explicit = SUNTRUE;
    step_mem.implicit = SUNTRUE;

    ark_mem.step_mem = Some(step_mem);

    /* re-attach internal error weight functions if necessary */
    ark_reattach_tolerances(ark_mem)
}

/* shared tail of ARKStepSetImplicit/SetImEx: re-attach internal
   error weight functions if necessary */
fn ark_reattach_tolerances(ark_mem: &mut ARKodeMem) -> i32 {
    use crate::arkode::{ARKodeSStolerances, ARKodeSVtolerances};
    use crate::arkode_impl::ARK_SV;

    if !ark_mem.user_efun {
        let retval;
        if ark_mem.itol == ARK_SV && ark_mem.Vabstol.is_some() {
            let vabstol = ark_mem.Vabstol.take().unwrap();
            retval = ARKodeSVtolerances(ark_mem, ark_mem.reltol, &vabstol);
            /* (ARKodeSVtolerances stores its own copy; restore the
               original slot only if the call failed to do so) */
            if ark_mem.Vabstol.is_none() {
                ark_mem.Vabstol = Some(vabstol);
            }
        } else {
            retval = ARKodeSStolerances(ark_mem, ark_mem.reltol, ark_mem.Sabstol);
        }
        if retval != ARK_SUCCESS {
            return retval;
        }
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKStepSetTables:

  Specifies to use customized Butcher tables for the system.

  If Bi is NULL, then this sets the integrator in 'explicit' mode.

  If Be is NULL, then this sets the integrator in 'implicit' mode.

  Returns ARK_ILL_INPUT if both Butcher tables are not supplied.
  ---------------------------------------------------------------*/
pub fn ARKStepSetTables(
    ark_mem: &mut ARKodeMem,
    q: i32,
    p: i32,
    Bi: Option<&ARKodeButcherTable>,
    Be: Option<&ARKodeButcherTable>,
) -> i32 {
    use crate::arkode_butcher::ARKodeButcherTable_Copy;

    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "ARKStepSetTables") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* check for illegal inputs */
    if Bi.is_none() && Be.is_none() {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_NULL,
            line!(),
            "ARKStepSetTables",
            file!(),
            "At least one complete table must be supplied",
        );
        return ARK_ILL_INPUT;
    }

    /* if both tables are set, check that they have the same number of stages */
    if let (Some(bi), Some(be)) = (Bi, Be) {
        if bi.stages != be.stages {
            ark_mem.step_mem = Some(step_mem);
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_NULL,
                line!(),
                "ARKStepSetTables",
                file!(),
                "Both tables must have the same number of stages",
            );
            return ARK_ILL_INPUT;
        }
    }

    /* clear any existing parameters and Butcher tables */
    step_mem.stages = 0;
    step_mem.q = 0;
    step_mem.p = 0;

    if let Some(bt) = step_mem.Be.take() {
        let (mut bliw, mut blrw) = (0i64, 0i64);
        ARKodeButcherTable_Space(&bt, &mut bliw, &mut blrw);
        ark_mem.liw -= bliw;
        ark_mem.lrw -= blrw;
    }
    if let Some(bt) = step_mem.Bi.take() {
        let (mut bliw, mut blrw) = (0i64, 0i64);
        ARKodeButcherTable_Space(&bt, &mut bliw, &mut blrw);
        ark_mem.liw -= bliw;
        ark_mem.lrw -= blrw;
    }

    /*
     * determine mode (implicit/explicit/ImEx), and perform appropriate actions
     */
    if let (None, Some(be)) = (Bi, Be) {
        /* explicit: set the relevant parameters (use table q and p) */
        step_mem.stages = be.stages;
        step_mem.q = be.q;
        step_mem.p = be.p;

        /* copy the table in step memory */
        step_mem.Be = ARKodeButcherTable_Copy(be);

        ark_mem.step_mem = Some(step_mem);

        /* set method as purely explicit */
        let retval = ARKStepSetExplicit(ark_mem);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "ARKStepSetTables",
                file!(),
                "Error in ARKStepSetExplicit",
            );
            return retval;
        }
    } else if let (Some(bi), None) = (Bi, Be) {
        /* implicit: set the relevant parameters (use table q and p) */
        step_mem.stages = bi.stages;
        step_mem.q = bi.q;
        step_mem.p = bi.p;

        /* copy the table in step memory */
        step_mem.Bi = ARKodeButcherTable_Copy(bi);

        ark_mem.step_mem = Some(step_mem);

        /* set method as purely implicit */
        let retval = ARKStepSetImplicit(ark_mem);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "ARKStepSetTables",
                file!(),
                "Error in ARKStepSetImplicit",
            );
            return ARK_ILL_INPUT;
        }
    } else {
        /* ImEx: set the relevant parameters (use input q and p) */
        let (bi, be) = (Bi.unwrap(), Be.unwrap());
        step_mem.stages = bi.stages;
        step_mem.q = q;
        step_mem.p = p;

        /* copy the explicit and implicit tables into step memory */
        step_mem.Be = ARKodeButcherTable_Copy(be);
        step_mem.Bi = ARKodeButcherTable_Copy(bi);

        ark_mem.step_mem = Some(step_mem);

        /* set method as ImEx */
        let retval = ARKStepSetImEx(ark_mem);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "ARKStepSetTables",
                file!(),
                "Error in ARKStepSetImEx",
            );
            return ARK_ILL_INPUT;
        }
    }

    /* note Butcher table space requirements */
    let step_mem = arkStep_AccessStepMem(ark_mem, "ARKStepSetTables").unwrap();
    let (mut bliw, mut blrw) = (0i64, 0i64);
    if let Some(be) = &step_mem.Be {
        ARKodeButcherTable_Space(be, &mut bliw, &mut blrw);
        ark_mem.liw += bliw;
        ark_mem.lrw += blrw;
    }
    if let Some(bi) = &step_mem.Bi {
        ARKodeButcherTable_Space(bi, &mut bliw, &mut blrw);
        ark_mem.liw += bliw;
        ark_mem.lrw += blrw;
    }
    ark_mem.step_mem = Some(step_mem);

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKStepSetTableNum:

  Specifies to use pre-existing Butcher tables for the system,
  based on the integer flags passed to
  ARKodeButcherTable_LoadERK() and ARKodeButcherTable_LoadDIRK()
  within the files arkode_butcher_erk.rs and arkode_butcher_dirk.rs
  (automatically calls ARKStepSetImEx).

  If either argument is negative (illegal), then this disables the
  corresponding table (e.g. itable = -1  ->  explicit)

  Note: this routine should NOT be used in conjunction with
  ARKodeSetOrder.
  ---------------------------------------------------------------*/
pub fn ARKStepSetTableNum(
    ark_mem: &mut ARKodeMem,
    itable: crate::arkode_butcher_dirk::ARKODE_DIRKTableID,
    etable: crate::arkode_butcher_erk::ARKODE_ERKTableID,
) -> i32 {
    use crate::arkode_butcher_dirk::{
        ARKodeButcherTable_LoadDIRK, ARKODE_ARK2_DIRK_3_1_2, ARKODE_ARK324L2SA_DIRK_4_2_3,
        ARKODE_ARK436L2SA_DIRK_6_3_4, ARKODE_ARK437L2SA_DIRK_7_3_4, ARKODE_ARK548L2SA_DIRK_8_4_5,
        ARKODE_ARK548L2SAb_DIRK_8_4_5, ARKODE_MAX_DIRK_NUM, ARKODE_MIN_DIRK_NUM,
    };
    use crate::arkode_butcher_erk::{
        ARKodeButcherTable_LoadERK, ARKODE_ARK2_ERK_3_1_2, ARKODE_ARK324L2SA_ERK_4_2_3,
        ARKODE_ARK436L2SA_ERK_6_3_4, ARKODE_ARK437L2SA_ERK_7_3_4, ARKODE_ARK548L2SA_ERK_8_4_5,
        ARKODE_ARK548L2SAb_ERK_8_4_5, ARKODE_MAX_ERK_NUM, ARKODE_MIN_ERK_NUM,
    };

    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "ARKStepSetTableNum") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* clear any existing parameters and Butcher tables */
    step_mem.stages = 0;
    step_mem.q = 0;
    step_mem.p = 0;

    if let Some(bt) = step_mem.Be.take() {
        let (mut bliw, mut blrw) = (0i64, 0i64);
        ARKodeButcherTable_Space(&bt, &mut bliw, &mut blrw);
        ark_mem.liw -= bliw;
        ark_mem.lrw -= blrw;
    }
    if let Some(bt) = step_mem.Bi.take() {
        let (mut bliw, mut blrw) = (0i64, 0i64);
        ARKodeButcherTable_Space(&bt, &mut bliw, &mut blrw);
        ark_mem.liw -= bliw;
        ark_mem.lrw -= blrw;
    }

    /* determine mode (implicit/explicit/ImEx), and perform
       appropriate actions  */

    if itable < 0 && etable < 0 {
        /*     illegal inputs */
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_NULL,
            line!(),
            "ARKStepSetTableNum",
            file!(),
            "At least one valid table number must be supplied",
        );
        return ARK_ILL_INPUT;
    } else if itable < 0 {
        /* explicit: check that argument specifies an explicit table */
        if etable < ARKODE_MIN_ERK_NUM || etable > ARKODE_MAX_ERK_NUM {
            ark_mem.step_mem = Some(step_mem);
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_NULL,
                line!(),
                "ARKStepSetTableNum",
                file!(),
                "Illegal ERK table number",
            );
            return ARK_ILL_INPUT;
        }

        /* fill in table based on argument */
        step_mem.Be = ARKodeButcherTable_LoadERK(etable);
        if step_mem.Be.is_none() {
            ark_mem.step_mem = Some(step_mem);
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_NULL,
                line!(),
                "ARKStepSetTableNum",
                file!(),
                "Error setting explicit table with that index",
            );
            return ARK_ILL_INPUT;
        }
        {
            let be = step_mem.Be.as_ref().unwrap();
            step_mem.stages = be.stages;
            step_mem.q = be.q;
            step_mem.p = be.p;
        }

        ark_mem.step_mem = Some(step_mem);

        /* set method as purely explicit */
        let flag = ARKStepSetExplicit(ark_mem);
        if flag != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "ARKStepSetTableNum",
                file!(),
                "Error in ARKStepSetExplicit",
            );
            return flag;
        }
    } else if etable < 0 {
        /* implicit: check that argument specifies an implicit table */
        if itable < ARKODE_MIN_DIRK_NUM || itable > ARKODE_MAX_DIRK_NUM {
            ark_mem.step_mem = Some(step_mem);
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_NULL,
                line!(),
                "ARKStepSetTableNum",
                file!(),
                "Illegal IRK table number",
            );
            return ARK_ILL_INPUT;
        }

        /* fill in table based on argument */
        step_mem.Bi = ARKodeButcherTable_LoadDIRK(itable);
        if step_mem.Bi.is_none() {
            ark_mem.step_mem = Some(step_mem);
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_NULL,
                line!(),
                "ARKStepSetTableNum",
                file!(),
                "Error setting table with that index",
            );
            return ARK_ILL_INPUT;
        }
        {
            let bi = step_mem.Bi.as_ref().unwrap();
            step_mem.stages = bi.stages;
            step_mem.q = bi.q;
            step_mem.p = bi.p;
        }

        ark_mem.step_mem = Some(step_mem);

        /* set method as purely implicit */
        let flag = ARKStepSetImplicit(ark_mem);
        if flag != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "ARKStepSetTableNum",
                file!(),
                "Error in ARKStepSetImplicit",
            );
            return flag;
        }
    } else {
        /* ImEx: ensure that tables match (C's chained !(..)&&!(..) shape) */
        #[allow(clippy::nonminimal_bool)]
        let incompatible = !(etable == ARKODE_ARK324L2SA_ERK_4_2_3 && itable == ARKODE_ARK324L2SA_DIRK_4_2_3)
            && !(etable == ARKODE_ARK436L2SA_ERK_6_3_4 && itable == ARKODE_ARK436L2SA_DIRK_6_3_4)
            && !(etable == ARKODE_ARK437L2SA_ERK_7_3_4 && itable == ARKODE_ARK437L2SA_DIRK_7_3_4)
            && !(etable == ARKODE_ARK548L2SA_ERK_8_4_5 && itable == ARKODE_ARK548L2SA_DIRK_8_4_5)
            && !(etable == ARKODE_ARK548L2SAb_ERK_8_4_5
                && itable == ARKODE_ARK548L2SAb_DIRK_8_4_5)
            && !(etable == ARKODE_ARK2_ERK_3_1_2 && itable == ARKODE_ARK2_DIRK_3_1_2);
        if incompatible {
            ark_mem.step_mem = Some(step_mem);
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "ARKStepSetTableNum",
                file!(),
                "Incompatible Butcher tables for ARK method",
            );
            return ARK_ILL_INPUT;
        }

        /* fill in tables based on arguments */
        step_mem.Bi = ARKodeButcherTable_LoadDIRK(itable);
        step_mem.Be = ARKodeButcherTable_LoadERK(etable);
        if step_mem.Bi.is_none() {
            ark_mem.step_mem = Some(step_mem);
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_NULL,
                line!(),
                "ARKStepSetTableNum",
                file!(),
                "Illegal IRK table number",
            );
            return ARK_ILL_INPUT;
        }
        if step_mem.Be.is_none() {
            ark_mem.step_mem = Some(step_mem);
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_NULL,
                line!(),
                "ARKStepSetTableNum",
                file!(),
                "Illegal ERK table number",
            );
            return ARK_ILL_INPUT;
        }
        {
            let bi = step_mem.Bi.as_ref().unwrap();
            step_mem.stages = bi.stages;
            step_mem.q = bi.q;
            step_mem.p = bi.p;
        }

        ark_mem.step_mem = Some(step_mem);

        /* set method as ImEx */
        if ARKStepSetImEx(ark_mem) != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "ARKStepSetTableNum",
                file!(),
                MSG_ARK_MISSING_F,
            );
            return ARK_ILL_INPUT;
        }
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKStepSetTableName:

  Specifies to use pre-existing Butcher tables for the system,
  based on the string passed to
  ARKodeButcherTable_LoadERKByName() and
  ARKodeButcherTable_LoadDIRKByName() (automatically calls
  ARKStepSetImEx).

  If itable is "ARKODE_DIRK_NONE" or etable is "ARKODE_ERK_NONE",
  then this disables the corresponding table.
  ---------------------------------------------------------------*/
pub fn ARKStepSetTableName(ark_mem: &mut ARKodeMem, itable: &str, etable: &str) -> i32 {
    ARKStepSetTableNum(
        ark_mem,
        crate::arkode_butcher_dirk::arkButcherTableDIRKNameToID(itable),
        crate::arkode_butcher_erk::arkButcherTableERKNameToID(etable),
    )
}

/*---------------------------------------------------------------
  arkStep_SetOptions:

  Provides command-line control over ARKStep-specific "set"
  routines (arkode_arkstep_io.c).
  ---------------------------------------------------------------*/
pub fn arkStep_SetOptions(
    ark_mem: &mut ARKodeMem,
    argidx: &mut usize,
    argv: &[String],
    offset: usize,
    arg_used: &mut bool,
) -> i32 {
    use crate::sundials_cli::{
        sunCheckAndSetActionArgs, sunCheckAndSetTwoCharArgs, sunKeyActionPair, sunKeyTwoCharPair,
    };

    /* Set lists of keys, and the corresponding set routines */
    let twochar_pairs: [sunKeyTwoCharPair<ARKodeMem>; 1] = [sunKeyTwoCharPair {
        key: "table_names",
        set: ARKStepSetTableName,
    }];

    let action_pairs: [sunKeyActionPair<ARKodeMem>; 3] = [
        sunKeyActionPair { key: "explicit", set: ARKStepSetExplicit },
        sunKeyActionPair { key: "implicit", set: ARKStepSetImplicit },
        sunKeyActionPair { key: "imex", set: ARKStepSetImEx },
    ];

    /* check all "twochar" keys */
    let mut j: usize = 0;
    let retval = sunCheckAndSetTwoCharArgs(
        ark_mem,
        argidx,
        argv,
        offset,
        &twochar_pairs,
        arg_used,
        &mut j,
    );
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            retval,
            line!(),
            "arkStep_SetOptions",
            file!(),
            &format!(
                "error setting command-line argument: {}",
                twochar_pairs[j].key
            ),
        );
        return retval;
    }
    if *arg_used {
        return ARK_SUCCESS;
    }

    /* check all action keys */
    let retval = sunCheckAndSetActionArgs(
        ark_mem,
        argidx,
        argv,
        offset,
        &action_pairs,
        arg_used,
        &mut j,
    );
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            retval,
            line!(),
            "arkStep_SetOptions",
            file!(),
            &format!(
                "error setting command-line argument: {}",
                action_pairs[j].key
            ),
        );
        return retval;
    }

    ARK_SUCCESS
}
