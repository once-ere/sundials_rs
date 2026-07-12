/* -----------------------------------------------------------------
 * Translation of src/arkode/arkode_sprkstep.c (ARKODE 7.7.0).
 * SPRKStep time-stepper module: symplectic partitioned Runge-Kutta
 * methods for separable Hamiltonian systems, with an optional
 * compensated-summation step form.
 *
 * step_mem access follows the erkstep take/put-back convention
 * (ARCHITECTURE.md Addendum C.1): wrappers take the Box<dyn Any>
 * out of ark_mem.step_mem, downcast, run an _inner worker, and put
 * it back.
 * -----------------------------------------------------------------*/

use crate::arkode::{arkAllocVec, arkCreate, arkFreeVec, arkInit};
use crate::arkode_impl::*;
use crate::arkode_io::ARKodeSetInterpolantType;
use crate::arkode_sprk::{ARKodeSPRKTable_Free, ARKodeSPRKTable_Load};
use crate::arkode_sprkstep_impl::*;
use crate::arkode_sprkstep_io::{
    sprkStep_GetNumRhsEvals, sprkStep_GetStageIndex, sprkStep_PrintAllStats, sprkStep_SetDefaults,
    sprkStep_SetOrder, sprkStep_SetUseCompensatedSums, sprkStep_WriteParameters,
};
use crate::nvector_serial::{N_VConst, N_VLinearSum, NVector};
use crate::sundials_context::SUNContext;
use crate::sundials_math::SUNRabs;

/*===============================================================
  Exported functions
  ===============================================================*/

pub fn SPRKStepCreate(
    f1: Option<ARKRhsFn>,
    f2: Option<ARKRhsFn>,
    t0: f64,
    y0: &NVector,
    sunctx: &SUNContext,
) -> Option<Box<ARKodeMem>> {
    /* Check that f1 and f2 are supplied */
    if f1.is_none() {
        arkProcessError(None, ARK_ILL_INPUT, line!(), "SPRKStepCreate", file!(), MSG_ARK_NULL_F);
        return None;
    }

    if f2.is_none() {
        arkProcessError(None, ARK_ILL_INPUT, line!(), "SPRKStepCreate", file!(), MSG_ARK_NULL_F);
        return None;
    }

    /* Check for legal input parameters */
    if y0.data.is_empty() {
        arkProcessError(None, ARK_ILL_INPUT, line!(), "SPRKStepCreate", file!(), MSG_ARK_NULL_Y0);
        return None;
    }

    /* Create ark_mem structure and set default values */
    let mut ark_mem = arkCreate(sunctx);

    /* Allocate ARKodeSPRKStepMem structure, and initialize to zero */
    let mut step_mem = Box::new(ARKodeSPRKStepMem::default());

    /* Allocate vectors in stepper mem */
    let tmpl_len = y0.data.len();
    arkAllocVec(&mut ark_mem, tmpl_len, &mut step_mem.sdata);

    if ark_mem.use_compensated_sums {
        arkAllocVec(&mut ark_mem, tmpl_len, &mut step_mem.yerr);
        /* Zero yerr for compensated summation */
        N_VConst(ZERO, &mut step_mem.yerr);
    }
    /* (else yerr stays empty = C NULL) */

    /* Attach step_mem structure and function pointers to ark_mem */
    ark_mem.step_init = Some(sprkStep_Init);
    ark_mem.step_fullrhs = Some(sprkStep_FullRHS);
    ark_mem.step = Some(sprkStep_TakeStep);
    ark_mem.step_printallstats = Some(sprkStep_PrintAllStats);
    ark_mem.step_writeparameters = Some(sprkStep_WriteParameters);
    ark_mem.step_setusecompensatedsums = Some(sprkStep_SetUseCompensatedSums);
    ark_mem.step_free = Some(sprkStep_Free);
    ark_mem.step_setdefaults = Some(sprkStep_SetDefaults);
    ark_mem.step_setorder = Some(sprkStep_SetOrder);
    ark_mem.step_getnumrhsevals = Some(sprkStep_GetNumRhsEvals);
    ark_mem.step_getstageindex = Some(sprkStep_GetStageIndex);
    ark_mem.step_mem = Some(step_mem);

    /* Set default values for optional inputs */
    let retval = sprkStep_SetDefaults(&mut ark_mem);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(&ark_mem),
            retval,
            line!(),
            "SPRKStepCreate",
            file!(),
            "Error setting default solver options",
        );
        return None;
    }

    /* Copy the input parameters into ARKODE state; initialize counters */
    {
        let mut sm = sprkStep_AccessStepMem(&mut ark_mem, "SPRKStepCreate").unwrap();
        sm.f1 = f1;
        sm.f2 = f2;
        sm.nf1 = 0;
        sm.nf2 = 0;
        sm.istage = 0;
        ark_mem.step_mem = Some(sm);
    }

    /* SPRKStep uses Lagrange interpolation by default, since Hermite is
       less compatible with these methods. */
    ARKodeSetInterpolantType(&mut ark_mem, ARK_INTERP_LAGRANGE);

    /* Initialize main ARKODE infrastructure */
    let retval = arkInit(&mut ark_mem, t0, y0, FIRST_INIT);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(&ark_mem),
            retval,
            line!(),
            "SPRKStepCreate",
            file!(),
            "Unable to initialize main ARKODE infrastructure",
        );
        return None;
    }

    Some(ark_mem)
}

/*---------------------------------------------------------------
  SPRKStepReInit:

  This routine re-initializes the SPRKStep module to solve a new
  problem of the same size as was previously solved. Note all
  internal counters are set to 0 on re-initialization.
  ---------------------------------------------------------------*/
pub fn SPRKStepReInit(
    ark_mem: &mut ARKodeMem,
    f1: Option<ARKRhsFn>,
    f2: Option<ARKRhsFn>,
    t0: f64,
    y0: &NVector,
) -> i32 {
    let mut step_mem = match sprkStep_AccessStepMem(ark_mem, "SPRKStepReInit") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* Check if ark_mem was allocated */
    if !ark_mem.MallocDone {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_NO_MALLOC,
            line!(),
            "SPRKStepReInit",
            file!(),
            MSG_ARK_NO_MALLOC,
        );
        return ARK_NO_MALLOC;
    }

    /* Check that f1 and f2 are supplied */
    if f1.is_none() || f2.is_none() {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "SPRKStepReInit",
            file!(),
            MSG_ARK_NULL_F,
        );
        return ARK_ILL_INPUT;
    }

    /* Check that y0 is supplied */
    if y0.data.is_empty() {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "SPRKStepReInit",
            file!(),
            MSG_ARK_NULL_Y0,
        );
        return ARK_ILL_INPUT;
    }

    /* Copy the input parameters into ARKODE state */
    step_mem.f1 = f1;
    step_mem.f2 = f2;
    ark_mem.step_mem = Some(step_mem);

    /* Initialize main ARKODE infrastructure */
    let retval = arkInit(ark_mem, t0, y0, FIRST_INIT);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            retval,
            line!(),
            "SPRKStepReInit",
            file!(),
            "Unable to reinitialize main ARKODE infrastructure",
        );
        return retval;
    }

    /* Initialize the counters; zero yerr for compensated summation */
    let mut step_mem = sprkStep_AccessStepMem(ark_mem, "SPRKStepReInit").unwrap();
    step_mem.nf1 = 0;
    step_mem.nf2 = 0;
    step_mem.istage = 0;
    if ark_mem.use_compensated_sums {
        N_VConst(ZERO, &mut step_mem.yerr);
    }
    ark_mem.step_mem = Some(step_mem);

    ARK_SUCCESS
}

/*===============================================================
  Interface routines supplied to ARKODE
  ===============================================================*/

/*---------------------------------------------------------------
  sprkStep_Free frees all SPRKStep memory.
  ---------------------------------------------------------------*/
pub fn sprkStep_Free(ark_mem: &mut ARKodeMem) {
    /* conditional frees on non-NULL SPRKStep module */
    if let Some(mut step_mem) = sprkStep_AccessStepMem(ark_mem, "sprkStep_Free") {
        let mut sdata = std::mem::take(&mut step_mem.sdata);
        arkFreeVec(ark_mem, &mut sdata);
        if !step_mem.yerr.data.is_empty() {
            let mut yerr = std::mem::take(&mut step_mem.yerr);
            arkFreeVec(ark_mem, &mut yerr);
        }
        if let Some(method) = step_mem.method.take() {
            ARKodeSPRKTable_Free(method);
        }
        /* the box is dropped here (C: free(ark_mem->step_mem)) */
    }
    ark_mem.step_mem = None;
}

/*---------------------------------------------------------------
  sprkStep_Init:

  This routine is called just prior to performing internal time
  steps (after all user "set" routines have been called) from
  within arkInitialSetup.

  With initialization type FIRST_INIT or RESIZE_INIT, this routine
  loads the default method of the selected order if necessary.

  With initialization type RESET_INIT, this routine does nothing.
  ---------------------------------------------------------------*/
pub fn sprkStep_Init(ark_mem: &mut ARKodeMem, init_type: i32) -> i32 {
    let mut step_mem = match sprkStep_AccessStepMem(ark_mem, "sprkStep_Init") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };
    let ret = sprkStep_Init_inner(ark_mem, &mut step_mem, init_type);
    ark_mem.step_mem = Some(step_mem);
    ret
}

fn sprkStep_Init_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeSPRKStepMem>,
    init_type: i32,
) -> i32 {
    /* immediately return if reset */
    if init_type == RESET_INIT {
        return ARK_SUCCESS;
    }

    /* initializations/checks for (re-)initialization call */
    if init_type == FIRST_INIT && step_mem.method.is_none() {
        step_mem.method = match step_mem.q {
            1 => ARKodeSPRKTable_Load(SPRKSTEP_DEFAULT_1),
            2 => ARKodeSPRKTable_Load(SPRKSTEP_DEFAULT_2),
            3 => ARKodeSPRKTable_Load(SPRKSTEP_DEFAULT_3),
            4 => ARKodeSPRKTable_Load(SPRKSTEP_DEFAULT_4),
            5 => ARKodeSPRKTable_Load(SPRKSTEP_DEFAULT_5),
            6 => ARKodeSPRKTable_Load(SPRKSTEP_DEFAULT_6),
            7 | 8 => ARKodeSPRKTable_Load(SPRKSTEP_DEFAULT_8),
            9 | 10 => ARKodeSPRKTable_Load(SPRKSTEP_DEFAULT_10),
            _ => {
                arkProcessError(
                    Some(ark_mem),
                    ARK_WARNING,
                    line!(),
                    "sprkStep_Init",
                    file!(),
                    "No SPRK method at requested order, using q=4.",
                );
                ARKodeSPRKTable_Load(SPRKSTEP_DEFAULT_4)
            }
        };
    }

    /* Override the interpolant degree (if needed), used in arkInitialSetup */
    let method_q = step_mem.method.as_ref().unwrap().q;
    if method_q > 1 && ark_mem.interp_degree > (method_q - 1) {
        /* Limit max degree to at most one less than the method global order */
        ark_mem.interp_degree = method_q - 1;
    } else if method_q == 1 && ark_mem.interp_degree > 1 {
        /* Allow for linear interpolant with first order methods to ensure
           solution values are returned at the time interval end points */
        ark_mem.interp_degree = 1;
    }

    /* Zero yerr for compensated summation */
    if ark_mem.use_compensated_sums {
        N_VConst(ZERO, &mut step_mem.yerr);
    }

    ARK_SUCCESS
}

/* Utility to call f1 and increment the counter */
pub(crate) fn sprkStep_f1(
    step_mem: &mut ARKodeSPRKStepMem,
    tcur: f64,
    ycur: &NVector,
    f1: &mut NVector,
    user_data: &mut crate::sundials_types::UserData,
) -> i32 {
    let retval = (step_mem.f1.unwrap())(tcur, ycur, f1, user_data);
    step_mem.nf1 += 1;
    retval
}

/* Utility to call f2 and increment the counter */
pub(crate) fn sprkStep_f2(
    step_mem: &mut ARKodeSPRKStepMem,
    tcur: f64,
    ycur: &NVector,
    f2: &mut NVector,
    user_data: &mut crate::sundials_types::UserData,
) -> i32 {
    let retval = (step_mem.f2.unwrap())(tcur, ycur, f2, user_data);
    step_mem.nf2 += 1;
    retval
}

/*------------------------------------------------------------------------------
  sprkStep_FullRHS:

  This is just a wrapper to call the user-supplied RHS,
  f1(t,y) + f2(t,y).  Since RHS values are not stored in SPRKStep
  we evaluate the RHS functions for all modes.
  ----------------------------------------------------------------------------*/
pub fn sprkStep_FullRHS(
    ark_mem: &mut ARKodeMem,
    t: f64,
    y: &NVector,
    f: &mut NVector,
    mode: i32,
) -> i32 {
    let mut step_mem = match sprkStep_AccessStepMem(ark_mem, "sprkStep_FullRHS") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };
    let ret = sprkStep_FullRHS_inner(ark_mem, &mut step_mem, t, y, f, mode);
    ark_mem.step_mem = Some(step_mem);
    ret
}

fn sprkStep_FullRHS_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeSPRKStepMem>,
    t: f64,
    y: &NVector,
    f: &mut NVector,
    mode: i32,
) -> i32 {
    /* perform RHS functions contingent on 'mode' argument */
    match mode {
        ARK_FULLRHS_START | ARK_FULLRHS_END | ARK_FULLRHS_OTHER => {
            /* call the user-supplied pre-RHS function (if supplied) */
            if let Some(pre_rhs) = ark_mem.PreRhsFn {
                let retval = pre_rhs(t, y, &mut ark_mem.user_data);
                if retval != 0 {
                    return ARK_PRERHSFN_FAIL;
                }
            }

            /* Since f1 and f2 do not have overlapping outputs the f vector
               is passed to both RHS functions. */

            let retval = sprkStep_f1(step_mem, t, y, f, &mut ark_mem.user_data);
            if retval != 0 {
                arkProcessError(
                    Some(ark_mem),
                    ARK_RHSFUNC_FAIL,
                    line!(),
                    "sprkStep_FullRHS",
                    file!(),
                    &format!(
                        "At t = {}, the right-hand side routine failed in an unrecoverable manner.",
                        t
                    ),
                );
                return ARK_RHSFUNC_FAIL;
            }

            let retval = sprkStep_f2(step_mem, t, y, f, &mut ark_mem.user_data);
            if retval != 0 {
                arkProcessError(
                    Some(ark_mem),
                    ARK_RHSFUNC_FAIL,
                    line!(),
                    "sprkStep_FullRHS",
                    file!(),
                    &format!(
                        "At t = {}, the right-hand side routine failed in an unrecoverable manner.",
                        t
                    ),
                );
                return ARK_RHSFUNC_FAIL;
            }
        }

        _ => {
            /* return with RHS failure if unknown mode is passed */
            arkProcessError(
                Some(ark_mem),
                ARK_RHSFUNC_FAIL,
                line!(),
                "sprkStep_FullRHS",
                file!(),
                "Unknown full RHS mode",
            );
            return ARK_RHSFUNC_FAIL;
        }
    }

    ARK_SUCCESS
}

/* Standard formulation of SPRK.
   This requires only 2 vectors in principle, but we use three
   since we persist the stage data. Only the stage data vector
   belongs to SPRKStep, the other two are reused from the ARKODE core.

   (C walks prev_stage/curr_stage pointers: prev = yn for the first
   stage and ycur afterwards, so the position/velocity updates
   alias ycur from stage 1 on — realized here with the in-place
   linear-sum family.) */
pub fn sprkStep_TakeStep(ark_mem: &mut ARKodeMem, dsmPtr: &mut f64, nflagPtr: &mut i32) -> i32 {
    let mut step_mem = match sprkStep_AccessStepMem(ark_mem, "sprkStep_TakeStep") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };
    let ret = sprkStep_TakeStep_inner(ark_mem, &mut step_mem, dsmPtr, nflagPtr);
    ark_mem.step_mem = Some(step_mem);
    ret
}

fn sprkStep_TakeStep_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeSPRKStepMem>,
    dsmPtr: &mut f64,
    nflagPtr: &mut i32,
) -> i32 {
    let mut ci = ZERO;
    let mut chati = ZERO;
    let stages = step_mem.method.as_ref().unwrap().stages;

    for is in 0..stages {
        /* load/compute coefficients */
        let ai = step_mem.method.as_ref().unwrap().a[is as usize];
        let ahati = step_mem.method.as_ref().unwrap().ahat[is as usize];

        ci += ai;
        chati += ahati;

        /* store current stage index */
        step_mem.istage = is;

        /* evaluate p' with the previous velocity */
        if SUNRabs(ahati) > TINY {
            N_VConst(ZERO, &mut step_mem.sdata); /* either have to do this or
                                                 ask user to set other outputs
                                                 to zero */

            /* call the user-supplied pre-RHS function (if supplied) */
            if let Some(pre_rhs) = ark_mem.PreRhsFn {
                let ARKodeMem { yn, ycur, user_data, tn, h, .. } = ark_mem;
                let prev: &NVector = if is == 0 { yn } else { ycur };
                let retval = pre_rhs(*tn + chati * *h, prev, user_data);
                if retval != 0 {
                    return ARK_PRERHSFN_FAIL;
                }
            }

            /* evaluate p' */
            let retval = {
                let f1fn = step_mem.f1.unwrap();
                let ARKodeMem { yn, ycur, user_data, tn, h, .. } = ark_mem;
                let prev: &NVector = if is == 0 { yn } else { ycur };
                let r = f1fn(*tn + chati * *h, prev, &mut step_mem.sdata, user_data);
                step_mem.nf1 += 1;
                r
            };
            if retval != 0 {
                return ARK_RHSFUNC_FAIL;
            }
        }

        /* position update */
        {
            let ARKodeMem { yn, ycur, h, .. } = ark_mem;
            if is == 0 {
                N_VLinearSum(ONE, yn, *h * ahati, &step_mem.sdata, ycur);
            } else {
                /* prev_stage == curr_stage: C's aliased
                   N_VLinearSum(1, z, b, sdata, z) => z += b*sdata */
                ycur.linear_sum_with(ONE, *h * ahati, &step_mem.sdata);
            }
        }

        /* set current stage time(s) */
        ark_mem.tcur = ark_mem.tn + chati * ark_mem.h;

        /* evaluate q' with the current positions and update velocity */
        if SUNRabs(ai) > TINY {
            N_VConst(ZERO, &mut step_mem.sdata); /* either have to do this or
                                                 ask user to set other outputs
                                                 to zero */

            /* call the user-supplied pre-RHS function (if supplied) */
            if let Some(pre_rhs) = ark_mem.PreRhsFn {
                let ARKodeMem { ycur, user_data, tn, h, .. } = ark_mem;
                let retval = pre_rhs(*tn + ci * *h, ycur, user_data);
                if retval != 0 {
                    return ARK_PRERHSFN_FAIL;
                }
            }

            /* evaluate q' */
            let retval = {
                let f2fn = step_mem.f2.unwrap();
                let ARKodeMem { ycur, user_data, tn, h, .. } = ark_mem;
                let r = f2fn(*tn + ci * *h, ycur, &mut step_mem.sdata, user_data);
                step_mem.nf2 += 1;
                r
            };
            if retval != 0 {
                return ARK_RHSFUNC_FAIL;
            }

            /* velocity update (aliased: z = z + h*ai*sdata) */
            ark_mem.ycur.linear_sum_with(ONE, ark_mem.h * ai, &step_mem.sdata);
        }

        /* apply user-supplied stage postprocessing function (if supplied) */
        if let (true, Some(post)) = (is < stages - 1, ark_mem.PostProcessStageFn) {
            let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
            let retval = post(*tcur, ycur, user_data);
            if retval != 0 {
                return ARK_POSTPROCESS_STAGE_FAIL;
            }
        } else if let (true, Some(post)) = (is == stages - 1, ark_mem.PostProcessStepFn) {
            let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
            let retval = post(*tcur, ycur, user_data);
            if retval != 0 {
                return ARK_POSTPROCESS_STEP_FAIL;
            }
        }
    }

    *nflagPtr = 0;
    *dsmPtr = ZERO;

    ARK_SUCCESS
}

/* Increment SPRK algorithm with compensated summation.
   This algorithm requires 6 vectors, but 5 of them are reused
   from the ARKODE core. */
pub fn sprkStep_TakeStep_Compensated(
    ark_mem: &mut ARKodeMem,
    dsmPtr: &mut f64,
    nflagPtr: &mut i32,
) -> i32 {
    let mut step_mem = match sprkStep_AccessStepMem(ark_mem, "sprkStep_TakeStep_Compensated") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };
    let ret = sprkStep_TakeStep_Compensated_inner(ark_mem, &mut step_mem, dsmPtr, nflagPtr);
    ark_mem.step_mem = Some(step_mem);
    ret
}

fn sprkStep_TakeStep_Compensated_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeSPRKStepMem>,
    dsmPtr: &mut f64,
    nflagPtr: &mut i32,
) -> i32 {
    let mut ci = ZERO;
    let mut chati = ZERO;
    let stages = step_mem.method.as_ref().unwrap().stages;

    /* Vector shortcuts (detached from the core for the step) */
    let mut delta_Yi = std::mem::take(&mut ark_mem.tempv1);
    let mut yn_plus_delta_Yi = std::mem::take(&mut ark_mem.tempv2);

    /* [ \Delta P_0 ] = [ 0 ]
       [ \Delta Q_0 ] = [ 0 ] */
    N_VConst(ZERO, &mut delta_Yi);

    /* if user-supplied stage preprocessing or postprocessing functions,
       we error out since those won't work with the increment form */
    if ark_mem.PreRhsFn.is_some()
        || ark_mem.PostProcessStageFn.is_some()
        || ark_mem.PostProcessStepFn.is_some()
    {
        ark_mem.tempv1 = delta_Yi;
        ark_mem.tempv2 = yn_plus_delta_Yi;
        arkProcessError(
            Some(ark_mem),
            ARK_POSTPROCESS_STAGE_FAIL,
            line!(),
            "sprkStep_TakeStep_Compensated",
            file!(),
            "Compensated summation is not compatible with stage Pre- or PostProcessing!\n",
        );
        return ARK_POSTPROCESS_STAGE_FAIL;
    }

    /* loop over internal stages to the step */
    for is in 0..stages {
        /* load/compute coefficients */
        let ai = step_mem.method.as_ref().unwrap().a[is as usize];
        let ahati = step_mem.method.as_ref().unwrap().ahat[is as usize];

        ci += ai;
        chati += ahati;

        /* store current stage index */
        step_mem.istage = is;

        /* [     ] + [            ]
           [ q_n ] + [ \Delta Q_i ] */
        N_VLinearSum(ONE, &ark_mem.yn, ONE, &delta_Yi, &mut yn_plus_delta_Yi);

        if SUNRabs(ahati) > TINY {
            /* Evaluate p' with the previous velocity */
            N_VConst(ZERO, &mut step_mem.sdata); /* either have to do this or
                                                 ask user to set other outputs
                                                 to zero */
            let retval = {
                let f1fn = step_mem.f1.unwrap();
                let ARKodeMem { user_data, tn, h, .. } = ark_mem;
                let r = f1fn(*tn + chati * *h, &yn_plus_delta_Yi, &mut step_mem.sdata, user_data);
                step_mem.nf1 += 1;
                r
            };
            if retval != 0 {
                ark_mem.tempv1 = delta_Yi;
                ark_mem.tempv2 = yn_plus_delta_Yi;
                return ARK_RHSFUNC_FAIL;
            }

            /* Incremental position update:
               [ \Delta P_i ] = [ \Delta P_{i-1} ] + [ sdata ] */
            delta_Yi.linear_sum_with(ONE, ark_mem.h * ahati, &step_mem.sdata);
        }

        /* [ p_n ] + [ \Delta P_i ]
           [     ] + [            ] */
        N_VLinearSum(ONE, &ark_mem.yn, ONE, &delta_Yi, &mut yn_plus_delta_Yi);

        /* set current stage time(s) */
        ark_mem.tcur = ark_mem.tn + chati * ark_mem.h;

        if SUNRabs(ai) > TINY {
            /* Evaluate q' with the current positions */
            N_VConst(ZERO, &mut step_mem.sdata); /* either have to do this or
                                                 ask user to set other outputs
                                                 to zero */
            let retval = {
                let f2fn = step_mem.f2.unwrap();
                let ARKodeMem { user_data, tn, h, .. } = ark_mem;
                let r = f2fn(*tn + ci * *h, &yn_plus_delta_Yi, &mut step_mem.sdata, user_data);
                step_mem.nf2 += 1;
                r
            };
            if retval != 0 {
                ark_mem.tempv1 = delta_Yi;
                ark_mem.tempv2 = yn_plus_delta_Yi;
                return ARK_RHSFUNC_FAIL;
            }

            /* Incremental velocity update:
               [ \Delta Q_i ] = [ \Delta Q_{i-1} ] + [ sdata ] */
            delta_Yi.linear_sum_with(ONE, ark_mem.h * ai, &step_mem.sdata);
        }
    }

    /*
      Now we compute the step solution via compensated summation.
       [ p_{n+1} ] = [ p_n ] + [ \Delta P_i ]
       [ q_{n+1} ] = [ q_n ] + [ \Delta Q_i ] */
    delta_Yi.linear_sum_with(ONE, -ONE, &step_mem.yerr);
    {
        let ARKodeMem { yn, ycur, .. } = ark_mem;
        N_VLinearSum(ONE, yn, ONE, &delta_Yi, ycur);
    }
    {
        /* diff = ycur - yn (in tempv3), then yerr = diff - delta_Yi */
        let mut diff = std::mem::take(&mut ark_mem.tempv3);
        {
            let ARKodeMem { yn, ycur, .. } = ark_mem;
            N_VLinearSum(ONE, ycur, -ONE, yn, &mut diff);
        }
        N_VLinearSum(ONE, &diff, -ONE, &delta_Yi, &mut step_mem.yerr);
        ark_mem.tempv3 = diff;
    }

    ark_mem.tempv1 = delta_Yi;
    ark_mem.tempv2 = yn_plus_delta_Yi;

    *nflagPtr = 0;
    *dsmPtr = ZERO;

    0
}

/*===============================================================
  Internal utility routines
  ===============================================================*/

/*---------------------------------------------------------------
  sprkStep_AccessStepMem:

  Shortcut routine to unpack the step_mem structure from ark_mem
  (take semantics; callers put it back).  If missing it reports
  ARK_MEM_NULL.
  ---------------------------------------------------------------*/
pub(crate) fn sprkStep_AccessStepMem(
    ark_mem: &mut ARKodeMem,
    fname: &str,
) -> Option<Box<ARKodeSPRKStepMem>> {
    let taken = ark_mem.step_mem.take();
    match taken {
        Some(b) => match b.downcast::<ARKodeSPRKStepMem>() {
            Ok(sm) => Some(sm),
            Err(other) => {
                ark_mem.step_mem = Some(other);
                arkProcessError(
                    Some(ark_mem),
                    ARK_MEM_NULL,
                    line!(),
                    fname,
                    file!(),
                    MSG_SPRKSTEP_NO_MEM,
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
                MSG_SPRKSTEP_NO_MEM,
            );
            None
        }
    }
}
