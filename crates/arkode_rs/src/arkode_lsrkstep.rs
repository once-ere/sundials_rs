/* -----------------------------------------------------------------
 * Translation of src/arkode/arkode_lsrkstep.c (ARKODE 7.7.0).
 * LSRKStep time-stepper module: low-storage Runge-Kutta methods —
 * super-time-stepping (RKC/RKL) with dominant-eigenvalue control,
 * and strong-stability-preserving (SSP) families.
 *
 * step_mem access follows the erkstep take/put-back convention.
 * RKC/RKL adaptation: C walks local copies of the tempv1/tempv2
 * pointers and swaps them each stage while the final embedding
 * section uses the ark_mem fields directly; here both vectors are
 * taken out, swapped as locals, and restored to their ORIGINAL
 * slots (parity-tracked) before the final-stage processing.
 * -----------------------------------------------------------------*/

use crate::arkode::{arkCreate, arkInit};
use crate::arkode_impl::*;
use crate::arkode_interp_impl::SIX;
use crate::arkode_io::ARKodeSetInterpolantType;
use crate::arkode_lsrkstep_impl::*;
use crate::arkode_lsrkstep_io::{
    lsrkStep_GetEstLocalErrors, lsrkStep_GetNumRhsEvals, lsrkStep_GetStageIndex,
    lsrkStep_PrintAllStats, lsrkStep_SetDefaults, lsrkStep_WriteParameters, LSRKStepSetSSPMethod,
    LSRKStepSetSTSMethod,
};
use crate::nvector_serial::{
    N_VConst, N_VLinearCombination, N_VLinearSum, N_VScale, N_VWrmsNorm, NVector,
};
use crate::sundials_context::SUNContext;
use crate::sundials_domeigestimator::{
    SUNDomEigEstimator_Estimate, SUNDomEigEstimator_GetNumIters, SUNDomEigEstimator_Initialize,
    SUNDomEigEstimator_SetNumPreprocessIters,
};
use crate::sundials_errors::SUN_SUCCESS;
use crate::sundials_math::{SUNMAX, SUNRabs, SUNRceil, SUNRsqrt};
use std::cell::RefCell;

/// C arkEwtSetSmallReal installed as an efun (see lsrkStep_Init).
fn lsrk_ewt_small_real(
    _ycur: &NVector,
    weight: &mut NVector,
    _e_data: &mut crate::sundials_types::UserData,
) -> i32 {
    N_VConst(crate::sundials_types::SUN_SMALL_REAL, weight);
    ARK_SUCCESS
}

/* private math helpers (arkode_lsrkstep_impl.h macros) */
fn SUNRlog(x: f64) -> f64 {
    x.ln()
}
fn SUNRsinh(x: f64) -> f64 {
    x.sinh()
}
fn SUNRcosh(x: f64) -> f64 {
    x.cosh()
}
fn SUNSQR(x: f64) -> f64 {
    x * x
}

/* N_VLinearCombination_Serial's z == X[0] branch:
   z *= c0; then z += c[i]*X[i] one operand vector at a time. */
fn lsrk_lincomb_inplace(z: &mut NVector, c0: f64, c_rest: &[f64], x_rest: &[&NVector]) {
    for e in z.data.iter_mut() {
        *e *= c0;
    }
    for (i, x) in x_rest.iter().enumerate() {
        for e in 0..z.data.len() {
            z.data[e] += c_rest[i] * x.data[e];
        }
    }
}

/*===============================================================
  Exported functions
  ===============================================================*/

pub fn LSRKStepCreateSTS(
    rhs: Option<ARKRhsFn>,
    t0: f64,
    y0: &NVector,
    sunctx: &SUNContext,
) -> Option<Box<ARKodeMem>> {
    /* Create shared LSRKStep memory structure */
    let mut ark_mem = lsrkStep_Create_Commons(rhs, t0, y0, sunctx)?;

    /* set default ARKODE_LSRK_RKC_2 method */
    let retval = LSRKStepSetSTSMethod(&mut ark_mem, ARKODE_LSRK_RKC_2);
    if retval != ARK_SUCCESS {
        lsrkStep_Free(&mut ark_mem);
        return None;
    }

    Some(ark_mem)
}

pub fn LSRKStepCreateSSP(
    rhs: Option<ARKRhsFn>,
    t0: f64,
    y0: &NVector,
    sunctx: &SUNContext,
) -> Option<Box<ARKodeMem>> {
    /* Create shared LSRKStep memory structure */
    let mut ark_mem = lsrkStep_Create_Commons(rhs, t0, y0, sunctx)?;

    /* set default ARKODE_LSRK_SSP_S_2 method */
    let retval = LSRKStepSetSSPMethod(&mut ark_mem, ARKODE_LSRK_SSP_S_2);
    if retval != ARK_SUCCESS {
        lsrkStep_Free(&mut ark_mem);
        return None;
    }

    Some(ark_mem)
}

/*---------------------------------------------------------------
  LSRKStepReInitSTS / LSRKStepReInitSSP:

  These routines re-initialize the LSRK module to solve a new
  problem of the same size as was previously solved.  All internal
  counters are set to 0 on re-initialization.
  ---------------------------------------------------------------*/
pub fn LSRKStepReInitSTS(
    ark_mem: &mut ARKodeMem,
    rhs: Option<ARKRhsFn>,
    t0: f64,
    y0: &NVector,
) -> i32 {
    lsrkStep_ReInit_Commons(ark_mem, rhs, t0, y0)
}

pub fn LSRKStepReInitSSP(
    ark_mem: &mut ARKodeMem,
    rhs: Option<ARKRhsFn>,
    t0: f64,
    y0: &NVector,
) -> i32 {
    lsrkStep_ReInit_Commons(ark_mem, rhs, t0, y0)
}

/*===============================================================
  Interface routines supplied to ARKODE
  ===============================================================*/

/*---------------------------------------------------------------
  lsrkStep_Create_Commons:

  A submodule for creating the common features of
  LSRKStepCreateSTS and LSRKStepCreateSSP.
  ---------------------------------------------------------------*/
pub(crate) fn lsrkStep_Create_Commons(
    rhs: Option<ARKRhsFn>,
    t0: f64,
    y0: &NVector,
    sunctx: &SUNContext,
) -> Option<Box<ARKodeMem>> {
    /* Check that rhs is supplied */
    if rhs.is_none() {
        arkProcessError(
            None,
            ARK_ILL_INPUT,
            line!(),
            "lsrkStep_Create_Commons",
            file!(),
            MSG_ARK_NULL_F,
        );
        return None;
    }

    /* Check for legal input parameters */
    if y0.data.is_empty() {
        arkProcessError(
            None,
            ARK_ILL_INPUT,
            line!(),
            "lsrkStep_Create_Commons",
            file!(),
            MSG_ARK_NULL_Y0,
        );
        return None;
    }

    /* Create ark_mem structure and set default values */
    let mut ark_mem = arkCreate(sunctx);

    /* Allocate ARKodeLSRKStepMem structure, and initialize to zero */
    let step_mem = Box::new(ARKodeLSRKStepMem::default());

    /* Attach step_mem structure and function pointers to ark_mem */
    ark_mem.step_init = Some(lsrkStep_Init);
    ark_mem.step_fullrhs = Some(lsrkStep_FullRHS);
    ark_mem.step = Some(lsrkStep_TakeStepRKC);
    ark_mem.step_printallstats = Some(lsrkStep_PrintAllStats);
    ark_mem.step_writeparameters = Some(lsrkStep_WriteParameters);
    ark_mem.step_free = Some(lsrkStep_Free);
    ark_mem.step_setdefaults = Some(lsrkStep_SetDefaults);
    ark_mem.step_getnumrhsevals = Some(lsrkStep_GetNumRhsEvals);
    ark_mem.step_getestlocalerrors = Some(lsrkStep_GetEstLocalErrors);
    ark_mem.step_getstageindex = Some(lsrkStep_GetStageIndex);
    ark_mem.step_mem = Some(step_mem);
    ark_mem.step_supports_adaptive = true;

    /* Set default values for optional inputs */
    let retval = lsrkStep_SetDefaults(&mut ark_mem);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(&ark_mem),
            retval,
            line!(),
            "lsrkStep_Create_Commons",
            file!(),
            "Error setting default solver options",
        );
        return None;
    }

    /* Copy the input parameters into ARKODE state; initialize the
       spectral info, flags and counters */
    {
        let mut sm = lsrkStep_AccessStepMem(&mut ark_mem, "lsrkStep_Create_Commons").unwrap();
        sm.fe = rhs;

        /* Initialize spectral radius info */
        sm.lambdaR = ZERO;
        sm.lambdaI = ZERO;
        sm.spectral_radius = ZERO;
        sm.spectral_radius_max = ZERO;
        sm.spectral_radius_min = ZERO;

        /* Initialize flags */
        sm.dom_eig_update = true;
        sm.dom_eig_is_current = false;
        sm.is_SSP = false;
        sm.init_warmup = true;

        /* Set NULL for dom_eig_fn; set NULL for DEE */
        sm.dom_eig_fn = None;
        sm.DEE = None;

        /* Initialize all the counters */
        sm.nfe = 0;
        sm.nfeDQ = 0;
        sm.stage_max = 0;
        sm.dom_eig_num_evals = 0;
        sm.stage_max_limit = STAGE_MAX_LIMIT_DEFAULT;
        sm.dom_eig_nst = 0;
        sm.num_dee_iters = 0;

        ark_mem.step_mem = Some(sm);
    }

    /* Initialize main ARKODE infrastructure */
    let retval = arkInit(&mut ark_mem, t0, y0, FIRST_INIT);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(&ark_mem),
            retval,
            line!(),
            "lsrkStep_Create_Commons",
            file!(),
            "Unable to initialize main ARKODE infrastructure",
        );
        return None;
    }

    /* Specify preferred interpolation type */
    ARKodeSetInterpolantType(&mut ark_mem, ARK_INTERP_LAGRANGE);

    Some(ark_mem)
}

/*---------------------------------------------------------------
  lsrkStep_ReInit_Commons:

  A submodule designed to reinitialize the common features of
  LSRKStepCreateSTS and LSRKStepCreateSSP.
  ---------------------------------------------------------------*/
pub(crate) fn lsrkStep_ReInit_Commons(
    ark_mem: &mut ARKodeMem,
    rhs: Option<ARKRhsFn>,
    t0: f64,
    y0: &NVector,
) -> i32 {
    let mut step_mem = match lsrkStep_AccessStepMem(ark_mem, "lsrkStep_ReInit_Commons") {
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
            "lsrkStep_ReInit_Commons",
            file!(),
            MSG_ARK_NO_MALLOC,
        );
        return ARK_NO_MALLOC;
    }

    /* Check that rhs is supplied */
    if rhs.is_none() {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "lsrkStep_ReInit_Commons",
            file!(),
            MSG_ARK_NULL_F,
        );
        return ARK_ILL_INPUT;
    }

    /* Check for legal input parameters */
    if y0.data.is_empty() {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "lsrkStep_ReInit_Commons",
            file!(),
            MSG_ARK_NULL_Y0,
        );
        return ARK_ILL_INPUT;
    }

    /* Copy the input parameters into ARKODE state */
    step_mem.fe = rhs;
    ark_mem.step_mem = Some(step_mem);

    /* Initialize main ARKODE infrastructure */
    let retval = arkInit(ark_mem, t0, y0, FIRST_INIT);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            retval,
            line!(),
            "lsrkStep_ReInit_Commons",
            file!(),
            "Unable to reinitialize main ARKODE infrastructure",
        );
        return retval;
    }

    /* Initialize all the counters, flags and stats */
    let mut sm = lsrkStep_AccessStepMem(ark_mem, "lsrkStep_ReInit_Commons").unwrap();
    sm.nfe = 0;
    sm.nfeDQ = 0;
    sm.dom_eig_num_evals = 0;
    sm.stage_max = 0;
    sm.lambdaR = ZERO;
    sm.lambdaI = ZERO;
    sm.spectral_radius = ZERO;
    sm.spectral_radius_max = 0.0;
    sm.spectral_radius_min = 0.0;
    sm.dom_eig_nst = 0;
    sm.num_dee_iters = 0;
    sm.dom_eig_update = true;
    sm.dom_eig_is_current = false;
    sm.init_warmup = true;
    ark_mem.step_mem = Some(sm);

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  lsrkStep_Init:

  This routine is called just prior to performing internal time
  steps (after all user "set" routines have been called) from
  within arkInitialSetup.

  With initialization types FIRST_INIT this routine performs
  setup and allocations needed for the method and sets
  the call_fullrhs flag.

  With other initialization types, this routine does nothing.
  ---------------------------------------------------------------*/
pub fn lsrkStep_Init(ark_mem: &mut ARKodeMem, init_type: i32) -> i32 {
    let mut step_mem = match lsrkStep_AccessStepMem(ark_mem, "lsrkStep_Init") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };
    let ret = lsrkStep_Init_inner(ark_mem, &mut step_mem, init_type);
    ark_mem.step_mem = Some(step_mem);
    ret
}

fn lsrkStep_Init_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeLSRKStepMem>,
    init_type: i32,
) -> i32 {
    /* immediately return if resize or reset */
    if init_type == RESIZE_INIT || init_type == RESET_INIT {
        return ARK_SUCCESS;
    }

    /* enforce use of arkEwtSmallReal if using a fixed step size
       and an internal error weight function */
    if ark_mem.fixedstep && !ark_mem.user_efun {
        ark_mem.user_efun = false;
        ark_mem.efun = Some(lsrk_ewt_small_real);
    }

    /* Check if user has provided dom_eig_fn or DEE */
    if !step_mem.is_SSP && step_mem.dom_eig_fn.is_none() && step_mem.DEE.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_DOMEIG_FAIL,
            line!(),
            "lsrkStep_Init",
            file!(),
            "STS methods require either a user provided dominant eigenvalue function or a SUNDomEigEstimator",
        );
        return ARK_DOMEIG_FAIL;
    }

    /* Initialize the DEE */
    if let Some(dee) = step_mem.DEE.as_mut() {
        let retval = SUNDomEigEstimator_Initialize(dee);
        if retval != SUN_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_DEE_FAIL,
                line!(),
                "lsrkStep_Init",
                file!(),
                "SUNDomEigEstimator_Initialize failed",
            );
            return ARK_DEE_FAIL;
        }

        /* Set number of DEE preprocessing iterations for the initial estimate */
        let retval = SUNDomEigEstimator_SetNumPreprocessIters(dee, step_mem.num_init_warmups);
        if retval != SUN_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_DEE_FAIL,
                line!(),
                "lsrkStep_Init",
                file!(),
                "SUNDomEigEstimator_SetNumPreprocessIters failed",
            );
            return ARK_DEE_FAIL;
        }
    }

    /* Allocate reusable arrays for fused vector interface */
    if step_mem.cvals.is_empty() {
        step_mem.cvals = vec![0.0; step_mem.nfusedopvecs as usize];
        ark_mem.lrw += step_mem.nfusedopvecs as i64;
        /* (Xvecs pointer array: operand lists are assembled at the call
           sites; keep the C liw accounting — C allocates Xvecs alongside
           cvals and only then adds liw) */
        ark_mem.liw += step_mem.nfusedopvecs as i64;
    }

    /* While LSRKStep does not currently call the full RHS function
       directly (later optimizations might) we do need the fn vector to
       always be allocated.  Signaling to the shared arkode module that
       full RHS evaluations are required will ensure fn is always
       allocated. */
    ark_mem.call_fullrhs = true;

    ARK_SUCCESS
}

/*------------------------------------------------------------------------------
  lsrkStep_FullRHS:

  This is just a wrapper to call the user-supplied RHS function,
  f(t,y).  See the C source for the mode discussion.
  ----------------------------------------------------------------------------*/
pub fn lsrkStep_FullRHS(
    ark_mem: &mut ARKodeMem,
    t: f64,
    y: &NVector,
    f: &mut NVector,
    mode: i32,
) -> i32 {
    let mut step_mem = match lsrkStep_AccessStepMem(ark_mem, "lsrkStep_FullRHS") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };
    let ret = lsrkStep_FullRHS_inner(ark_mem, &mut step_mem, t, y, f, mode);
    ark_mem.step_mem = Some(step_mem);
    ret
}

fn lsrkStep_FullRHS_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeLSRKStepMem>,
    t: f64,
    y: &NVector,
    f: &mut NVector,
    mode: i32,
) -> i32 {
    /* perform RHS functions contingent on 'mode' argument */
    match mode {
        ARK_FULLRHS_START => {
            /* compute the RHS */
            if !ark_mem.fn_is_current {
                /* call the user-supplied pre-rhs function (if supplied) */
                if let Some(pre_rhs) = ark_mem.PreRhsFn {
                    let retval = pre_rhs(t, y, &mut ark_mem.user_data);
                    if retval != 0 {
                        return ARK_PRERHSFN_FAIL;
                    }
                }
                let fe = step_mem.fe.unwrap();
                let retval = fe(t, y, f, &mut ark_mem.user_data);
                step_mem.nfe += 1;
                if retval != 0 {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_RHSFUNC_FAIL,
                        line!(),
                        "lsrkStep_FullRHS",
                        file!(),
                        &format!(
                            "At t = {}, the right-hand side routine failed in an unrecoverable manner.",
                            t
                        ),
                    );
                    return ARK_RHSFUNC_FAIL;
                }
            }
        }

        ARK_FULLRHS_END => {
            /* No further action is needed if STS since the currently
               available STS methods evaluate the RHS at the end of each
               time step.  If the stepper is an SSP, fn is updated and
               reused at the beginning of the step unless
               ark_mem->fn_is_current is changed by ARKODE. */
            if step_mem.is_SSP {
                /* call the user-supplied pre-rhs function (if supplied) */
                if let Some(pre_rhs) = ark_mem.PreRhsFn {
                    let retval = pre_rhs(t, y, &mut ark_mem.user_data);
                    if retval != 0 {
                        return ARK_PRERHSFN_FAIL;
                    }
                }
                let fe = step_mem.fe.unwrap();
                let retval = {
                    let ARKodeMem { fn_, user_data, .. } = ark_mem;
                    fe(t, y, fn_, user_data)
                };
                step_mem.nfe += 1;
                if retval != 0 {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_RHSFUNC_FAIL,
                        line!(),
                        "lsrkStep_FullRHS",
                        file!(),
                        &format!(
                            "At t = {}, the right-hand side routine failed in an unrecoverable manner.",
                            t
                        ),
                    );
                    return ARK_RHSFUNC_FAIL;
                }
                ark_mem.fn_is_current = true;
            }
            N_VScale(ONE, &ark_mem.fn_, f);
        }

        ARK_FULLRHS_OTHER => {
            /* call the user-supplied pre-rhs function (if supplied) */
            if let Some(pre_rhs) = ark_mem.PreRhsFn {
                let retval = pre_rhs(t, y, &mut ark_mem.user_data);
                if retval != 0 {
                    return ARK_PRERHSFN_FAIL;
                }
            }

            /* call f */
            let fe = step_mem.fe.unwrap();
            let retval = fe(t, y, f, &mut ark_mem.user_data);
            step_mem.nfe += 1;
            if retval != 0 {
                arkProcessError(
                    Some(ark_mem),
                    ARK_RHSFUNC_FAIL,
                    line!(),
                    "lsrkStep_FullRHS",
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
                "lsrkStep_FullRHS",
                file!(),
                "Unknown full RHS mode",
            );
            return ARK_RHSFUNC_FAIL;
        }
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  lsrkStep_TakeStepRKC:

  This routine performs a single RKC step.  ARK_RETRY_STEP
  indicates that the required stage number has reached the
  stage_max_limit with the current value of h; the step is then
  returned to adjust the step size.
  ---------------------------------------------------------------*/
pub fn lsrkStep_TakeStepRKC(ark_mem: &mut ARKodeMem, dsmPtr: &mut f64, nflagPtr: &mut i32) -> i32 {
    let mut step_mem = match lsrkStep_AccessStepMem(ark_mem, "lsrkStep_TakeStepRKC") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };
    let ret = lsrkStep_TakeStepRKC_inner(ark_mem, &mut step_mem, dsmPtr, nflagPtr);
    ark_mem.step_mem = Some(step_mem);
    ret
}

#[allow(clippy::too_many_lines)]
fn lsrkStep_TakeStepRKC_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeLSRKStepMem>,
    dsmPtr: &mut f64,
    nflagPtr: &mut i32,
) -> i32 {
    let onep54: f64 = 1.54;
    let c13: f64 = 13.0;
    let p8: f64 = 0.8;
    let p4: f64 = 0.4;

    /* initialize algebraic solver convergence flag to success,
       temporal error estimate to zero */
    *nflagPtr = ARK_SUCCESS;
    *dsmPtr = ZERO;

    /* Initialize the current stage index */
    step_mem.istage = 0;
    step_mem.req_stages = step_mem.stage_max_limit;

    /* Compute dominant eigenvalue and update stats */
    if step_mem.dom_eig_update {
        let retval = lsrkStep_ComputeNewDomEig(ark_mem, step_mem);
        if retval != ARK_SUCCESS {
            return retval;
        }
    }

    let ss = SUNRceil(SUNRsqrt(
        onep54 * SUNRabs(ark_mem.h) * step_mem.spectral_radius,
    ));
    let ss = SUNMAX(ss, 2.0);

    if ss >= step_mem.stage_max_limit as f64 {
        if !ark_mem.fixedstep {
            let hmax = ark_mem.hadapt_mem.as_ref().unwrap().safety
                * SUNSQR(step_mem.stage_max_limit as f64)
                / (onep54 * step_mem.spectral_radius);
            ark_mem.eta = hmax / ark_mem.h;
            *nflagPtr = ARK_RETRY_STEP;
            ark_mem.hadapt_mem.as_mut().unwrap().nst_exp += 1;
            return ARK_RETRY_STEP;
        } else {
            arkProcessError(
                Some(ark_mem),
                ARK_MAX_STAGE_LIMIT_FAIL,
                line!(),
                "lsrkStep_TakeStepRKC",
                file!(),
                "Unable to achieve stable results: Either reduce the step size or increase the stage_max_limit",
            );
            return ARK_MAX_STAGE_LIMIT_FAIL;
        }
    }

    step_mem.req_stages = ss as i32;
    step_mem.stage_max = std::cmp::max(step_mem.req_stages, step_mem.stage_max);

    /* Compute RHS function for the start of the step, if necessary. */
    if (!ark_mem.fn_is_current && ark_mem.initsetup) || (step_mem.step_nst != ark_mem.nst) {
        /* call the user-supplied pre-rhs function (if supplied) */
        if let Some(pre_rhs) = ark_mem.PreRhsFn {
            let ARKodeMem { yn, user_data, tn, .. } = ark_mem;
            let retval = pre_rhs(*tn, yn, user_data);
            if retval != 0 {
                return ARK_PRERHSFN_FAIL;
            }
        }

        /* call fe */
        let fe = step_mem.fe.unwrap();
        let retval = {
            let ARKodeMem { yn, fn_, user_data, tn, .. } = ark_mem;
            fe(*tn, yn, fn_, user_data)
        };
        step_mem.nfe += 1;
        if retval != ARK_SUCCESS {
            return ARK_RHSFUNC_FAIL;
        }
        ark_mem.fn_is_current = true;
    }

    /* Track the number of successful steps to determine if the previous
       step failed. */
    step_mem.step_nst = ark_mem.nst + 1;

    /* Initialize constants */
    let w0 = ONE + TWO / (c13 * SUNSQR(step_mem.req_stages as f64));
    let temp1 = SUNSQR(w0) - ONE;
    let temp2 = SUNRsqrt(temp1);
    let arg = step_mem.req_stages as f64 * SUNRlog(w0 + temp2);
    let w1 = SUNRsinh(arg) * temp1
        / (SUNRcosh(arg) * step_mem.req_stages as f64 * temp2 - w0 * SUNRsinh(arg));
    let mut bjm1 = ONE / SUNSQR(TWO * w0);
    let mut bjm2 = bjm1;
    let mut mus = w1 * bjm1;

    /* Take the tempv1/tempv2 buffers out; C swaps local pointer copies
       (parity-tracked so the buffers return to their original slots) */
    let mut tmp1 = std::mem::take(&mut ark_mem.tempv1);
    let mut tmp2 = std::mem::take(&mut ark_mem.tempv2);
    let mut swapped = false;

    macro_rules! restore_and_return {
        ($ret:expr) => {{
            if swapped {
                ark_mem.tempv1 = tmp2;
                ark_mem.tempv2 = tmp1;
            } else {
                ark_mem.tempv1 = tmp1;
                ark_mem.tempv2 = tmp2;
            }
            return $ret;
        }};
    }

    /* Begin stage 1 (store in tmp2) and initialize embedding */
    ark_mem.tcur = ark_mem.tn + ark_mem.h * mus;
    step_mem.istage = 1;
    {
        let ARKodeMem { yn, fn_, h, .. } = ark_mem;
        N_VLinearSum(ONE, yn, *h * mus, fn_, &mut tmp2);
        N_VScale(ONE, yn, &mut tmp1);
    }

    /* apply user-supplied stage postprocessing function (if supplied) */
    if let Some(post) = ark_mem.PostProcessStageFn {
        let ARKodeMem { user_data, tcur, .. } = ark_mem;
        let retval = post(*tcur, &tmp2, user_data);
        if retval != 0 {
            restore_and_return!(ARK_POSTPROCESS_STAGE_FAIL);
        }
    }

    /* Initialize constants for stage loop */
    let mut thjm2 = ZERO;
    let mut thjm1 = mus;
    let mut zjm1 = w0;
    let mut zjm2 = ONE;
    let mut dzjm1 = ONE;
    let mut dzjm2 = ZERO;
    let mut d2zjm1 = ZERO;
    let mut d2zjm2 = ZERO;

    /* Evaluate stages j = 2,...,step_mem->req_stages */
    for j in 2..=step_mem.req_stages {
        /* Complete the previous stage (evaluate the RHS and store it in
           ycur) */

        /* call the user-supplied pre-RHS function (if supplied) */
        if let Some(pre_rhs) = ark_mem.PreRhsFn {
            let ARKodeMem { user_data, tcur, .. } = ark_mem;
            let retval = pre_rhs(*tcur, &tmp2, user_data);
            if retval != 0 {
                restore_and_return!(ARK_PRERHSFN_FAIL);
            }
        }

        let fe = step_mem.fe.unwrap();
        let retval = {
            let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
            fe(*tcur, &tmp2, ycur, user_data)
        };
        step_mem.nfe += 1;

        if retval < 0 {
            restore_and_return!(ARK_RHSFUNC_FAIL);
        }
        if retval > 0 {
            restore_and_return!(RHSFUNC_RECVR);
        }

        /* Begin stage j (store in ycur) */
        let zj = TWO * w0 * zjm1 - zjm2;
        let dzj = TWO * w0 * dzjm1 - dzjm2 + TWO * zjm1;
        let d2zj = TWO * w0 * d2zjm1 - d2zjm2 + FOUR * dzjm1;
        let bj = d2zj / SUNSQR(dzj);
        let ajm1 = ONE - zjm1 * bjm1;
        let mu = TWO * w0 * bj / bjm1;
        let nu = -bj / bjm2;
        mus = mu * w1 / w0;
        let thj = mu * thjm1 + nu * thjm2 + mus * (ONE - ajm1);
        ark_mem.tcur = ark_mem.tn + ark_mem.h * thj;
        step_mem.istage = j;
        {
            /* N_VLinearCombination(5, ...) with z == X[0] == ycur */
            let ARKodeMem { ycur, yn, fn_, h, .. } = ark_mem;
            lsrk_lincomb_inplace(
                ycur,
                mus * *h,
                &[nu, ONE - mu - nu, mu, -mus * ajm1 * *h],
                &[&tmp1, yn, &tmp2, fn_],
            );
        }

        /* apply user-supplied stage or step postprocessing function (if
           supplied) */
        if let (true, Some(post)) = (j < step_mem.req_stages, ark_mem.PostProcessStageFn) {
            let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
            let retval = post(*tcur, ycur, user_data);
            if retval != 0 {
                restore_and_return!(ARK_POSTPROCESS_STAGE_FAIL);
            }
        } else if let (true, Some(post)) = (j == step_mem.req_stages, ark_mem.PostProcessStepFn) {
            let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
            let retval = post(*tcur, ycur, user_data);
            if retval != 0 {
                restore_and_return!(ARK_POSTPROCESS_STEP_FAIL);
            }
        }

        /* Shift the data for the next stage */
        if j < step_mem.req_stages {
            /* Swap tempv1 and tempv2 pointers to handle
               two-previous-stage logic */
            std::mem::swap(&mut tmp1, &mut tmp2);
            swapped = !swapped;

            N_VScale(ONE, &ark_mem.ycur, &mut tmp2);

            /* Update coefficients to handle the two-previous stage logic */
            thjm2 = thjm1;
            thjm1 = thj;
            bjm2 = bjm1;
            bjm1 = bj;
            zjm2 = zjm1;
            zjm1 = zj;
            dzjm2 = dzjm1;
            dzjm1 = dzj;
            d2zjm2 = d2zjm1;
            d2zjm1 = d2zj;
        }
    }

    /* restore the tempv buffers to their original slots before the
       final stage processing (which uses the ark_mem fields) */
    if swapped {
        ark_mem.tempv1 = tmp2;
        ark_mem.tempv2 = tmp1;
    } else {
        ark_mem.tempv1 = tmp1;
        ark_mem.tempv2 = tmp2;
    }

    /* final stage processing */
    ark_mem.tcur = ark_mem.tn + ark_mem.h;

    /* call the user-supplied pre-RHS function (if supplied) */
    if let Some(pre_rhs) = ark_mem.PreRhsFn {
        let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
        let retval = pre_rhs(*tcur, ycur, user_data);
        if retval != 0 {
            return ARK_PRERHSFN_FAIL;
        }
    }

    let fe = step_mem.fe.unwrap();
    let retval = {
        let ARKodeMem { ycur, tempv2, user_data, tcur, .. } = ark_mem;
        fe(*tcur, ycur, tempv2, user_data)
    };
    step_mem.nfe += 1;

    if retval < 0 {
        return ARK_RHSFUNC_FAIL;
    }
    if retval > 0 {
        return RHSFUNC_RECVR;
    }

    /* Compute yerr (if step adaptivity enabled) */
    if !ark_mem.fixedstep {
        /* Estimate the local error and compute its weighted RMS norm */
        {
            let ARKodeMem { yn, ycur, fn_, tempv1, tempv2, h, .. } = ark_mem;
            let c = [p8, -p8, p4 * *h, p4 * *h];
            let x: [&NVector; 4] = [yn, ycur, fn_, tempv2];
            N_VLinearCombination(4, &c, &x, tempv1);
        }
        *dsmPtr = N_VWrmsNorm(&ark_mem.tempv1, &ark_mem.ewt);
    }
    lsrkStep_DomEigUpdateLogic(ark_mem, step_mem, *dsmPtr);

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  lsrkStep_TakeStepRKL:

  This routine performs a single RKL step.
  ---------------------------------------------------------------*/
pub fn lsrkStep_TakeStepRKL(ark_mem: &mut ARKodeMem, dsmPtr: &mut f64, nflagPtr: &mut i32) -> i32 {
    let mut step_mem = match lsrkStep_AccessStepMem(ark_mem, "lsrkStep_TakeStepRKL") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };
    let ret = lsrkStep_TakeStepRKL_inner(ark_mem, &mut step_mem, dsmPtr, nflagPtr);
    ark_mem.step_mem = Some(step_mem);
    ret
}

#[allow(clippy::too_many_lines)]
fn lsrkStep_TakeStepRKL_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeLSRKStepMem>,
    dsmPtr: &mut f64,
    nflagPtr: &mut i32,
) -> i32 {
    let p8: f64 = 0.8;
    let p4: f64 = 0.4;

    /* initialize algebraic solver convergence flag to success,
       temporal error estimate to zero */
    *nflagPtr = ARK_SUCCESS;
    *dsmPtr = ZERO;

    /* Initialize the current stage index */
    step_mem.istage = 0;
    step_mem.req_stages = step_mem.stage_max_limit;

    /* Compute dominant eigenvalue and update stats */
    if step_mem.dom_eig_update {
        let retval = lsrkStep_ComputeNewDomEig(ark_mem, step_mem);
        if retval != ARK_SUCCESS {
            return retval;
        }
    }

    let ss = SUNRceil(
        (SUNRsqrt(9.0 + 8.0 * SUNRabs(ark_mem.h) * step_mem.spectral_radius) - ONE) / TWO,
    );
    let ss = SUNMAX(ss, 2.0);

    if ss >= step_mem.stage_max_limit as f64 {
        if !ark_mem.fixedstep {
            let hmax = ark_mem.hadapt_mem.as_ref().unwrap().safety
                * (SUNSQR(step_mem.stage_max_limit as f64) + step_mem.stage_max_limit as f64
                    - TWO)
                / (TWO * step_mem.spectral_radius);
            ark_mem.eta = hmax / ark_mem.h;
            *nflagPtr = ARK_RETRY_STEP;
            ark_mem.hadapt_mem.as_mut().unwrap().nst_exp += 1;
            return ARK_RETRY_STEP;
        } else {
            arkProcessError(
                Some(ark_mem),
                ARK_MAX_STAGE_LIMIT_FAIL,
                line!(),
                "lsrkStep_TakeStepRKL",
                file!(),
                "Unable to achieve stable results: Either reduce the step size or increase the stage_max_limit",
            );
            return ARK_MAX_STAGE_LIMIT_FAIL;
        }
    }

    step_mem.req_stages = ss as i32;
    step_mem.stage_max = std::cmp::max(step_mem.req_stages, step_mem.stage_max);

    /* Compute RHS function for the start of the step, if necessary. */
    if (!ark_mem.fn_is_current && ark_mem.initsetup) || (step_mem.step_nst != ark_mem.nst) {
        /* call the user-supplied pre-RHS function (if supplied) */
        if let Some(pre_rhs) = ark_mem.PreRhsFn {
            let ARKodeMem { yn, user_data, tn, .. } = ark_mem;
            let retval = pre_rhs(*tn, yn, user_data);
            if retval != 0 {
                return ARK_PRERHSFN_FAIL;
            }
        }
        let fe = step_mem.fe.unwrap();
        let retval = {
            let ARKodeMem { yn, fn_, user_data, tn, .. } = ark_mem;
            fe(*tn, yn, fn_, user_data)
        };
        step_mem.nfe += 1;
        if retval != ARK_SUCCESS {
            return ARK_RHSFUNC_FAIL;
        }
        ark_mem.fn_is_current = true;
    }

    /* Track the number of successful steps to determine if the previous
       step failed. */
    step_mem.step_nst = ark_mem.nst + 1;

    /* Initialize constants */
    let rs = step_mem.req_stages as f64;
    let w1 = FOUR / ((rs + TWO) * (rs - ONE));
    let mut bjm2 = ONE / THREE;
    let mut bjm1 = bjm2;
    let mut mus = w1 * bjm1;

    /* Take the tempv1/tempv2 buffers out (see RKC note) */
    let mut tmp1 = std::mem::take(&mut ark_mem.tempv1);
    let mut tmp2 = std::mem::take(&mut ark_mem.tempv2);
    let mut swapped = false;

    macro_rules! restore_and_return {
        ($ret:expr) => {{
            if swapped {
                ark_mem.tempv1 = tmp2;
                ark_mem.tempv2 = tmp1;
            } else {
                ark_mem.tempv1 = tmp1;
                ark_mem.tempv2 = tmp2;
            }
            return $ret;
        }};
    }

    /* Begin stage 1 (store in tmp2) and initialize embedding */
    ark_mem.tcur = ark_mem.tn + ark_mem.h * mus;
    step_mem.istage = 1;
    {
        let ARKodeMem { yn, fn_, h, .. } = ark_mem;
        N_VLinearSum(ONE, yn, *h * mus, fn_, &mut tmp2);
        N_VScale(ONE, yn, &mut tmp1);
    }

    /* apply user-supplied stage postprocessing function (if supplied) */
    if let Some(post) = ark_mem.PostProcessStageFn {
        let ARKodeMem { user_data, tcur, .. } = ark_mem;
        let retval = post(*tcur, &tmp2, user_data);
        if retval != 0 {
            restore_and_return!(ARK_POSTPROCESS_STAGE_FAIL);
        }
    }

    /* Evaluate stages j = 2,...,step_mem->req_stages */
    for j in 2..=step_mem.req_stages {
        /* Complete the previous stage (evaluate the RHS and store it in
           ycur) */

        /* call the user-supplied pre-RHS function (if supplied) */
        if let Some(pre_rhs) = ark_mem.PreRhsFn {
            let ARKodeMem { user_data, tcur, .. } = ark_mem;
            let retval = pre_rhs(*tcur, &tmp2, user_data);
            if retval != 0 {
                restore_and_return!(ARK_PRERHSFN_FAIL);
            }
        }

        let fe = step_mem.fe.unwrap();
        let retval = {
            let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
            fe(*tcur, &tmp2, ycur, user_data)
        };
        step_mem.nfe += 1;

        if retval < 0 {
            restore_and_return!(ARK_RHSFUNC_FAIL);
        }
        if retval > 0 {
            restore_and_return!(RHSFUNC_RECVR);
        }

        /* Begin stage j (store in ycur) */
        let jf = j as f64;
        let temj = (jf + TWO) * (jf - ONE);
        let bj = temj / (TWO * jf * (jf + ONE));
        let ajm1 = ONE - bjm1;
        let mu = (TWO * jf - ONE) / jf * (bj / bjm1);
        let nu = -(jf - ONE) / jf * (bj / bjm2);
        mus = w1 * mu;
        let cj = temj * w1 / FOUR;
        ark_mem.tcur = ark_mem.tn + ark_mem.h * cj;
        step_mem.istage = j;
        {
            /* N_VLinearCombination(5, ...) with z == X[0] == ycur */
            let ARKodeMem { ycur, yn, fn_, h, .. } = ark_mem;
            lsrk_lincomb_inplace(
                ycur,
                mus * *h,
                &[nu, ONE - mu - nu, mu, -mus * ajm1 * *h],
                &[&tmp1, yn, &tmp2, fn_],
            );
        }

        /* apply user-supplied stage or step postprocessing function (if
           supplied) */
        if let (true, Some(post)) = (j < step_mem.req_stages, ark_mem.PostProcessStageFn) {
            let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
            let retval = post(*tcur, ycur, user_data);
            if retval != 0 {
                restore_and_return!(ARK_POSTPROCESS_STAGE_FAIL);
            }
        } else if let (true, Some(post)) = (j == step_mem.req_stages, ark_mem.PostProcessStepFn) {
            let tpost = ark_mem.tcur + ark_mem.h * cj;
            let ARKodeMem { ycur, user_data, .. } = ark_mem;
            let retval = post(tpost, ycur, user_data);
            if retval != 0 {
                restore_and_return!(ARK_POSTPROCESS_STEP_FAIL);
            }
        }

        /* Shift the data for the next stage */
        if j < step_mem.req_stages {
            /* To avoid two data copies we swap ARKODE's tempv1 and
               tempv2 pointers */
            std::mem::swap(&mut tmp1, &mut tmp2);
            swapped = !swapped;

            N_VScale(ONE, &ark_mem.ycur, &mut tmp2);

            bjm2 = bjm1;
            bjm1 = bj;
        }
    }

    /* restore the tempv buffers (see RKC note) */
    if swapped {
        ark_mem.tempv1 = tmp2;
        ark_mem.tempv2 = tmp1;
    } else {
        ark_mem.tempv1 = tmp1;
        ark_mem.tempv2 = tmp2;
    }

    /* final stage processing */
    ark_mem.tcur = ark_mem.tn + ark_mem.h;

    /* call the user-supplied pre-RHS function (if supplied) */
    if let Some(pre_rhs) = ark_mem.PreRhsFn {
        let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
        let retval = pre_rhs(*tcur, ycur, user_data);
        if retval != 0 {
            return ARK_PRERHSFN_FAIL;
        }
    }
    let fe = step_mem.fe.unwrap();
    let retval = {
        let ARKodeMem { ycur, tempv2, user_data, tcur, .. } = ark_mem;
        fe(*tcur, ycur, tempv2, user_data)
    };
    step_mem.nfe += 1;

    if retval < 0 {
        return ARK_RHSFUNC_FAIL;
    }
    if retval > 0 {
        return RHSFUNC_RECVR;
    }

    /* Compute yerr (if step adaptivity enabled) */
    if !ark_mem.fixedstep {
        /* Estimate the local error and compute its weighted RMS norm */
        {
            let ARKodeMem { yn, ycur, fn_, tempv1, tempv2, h, .. } = ark_mem;
            let c = [p8, -p8, p4 * *h, p4 * *h];
            let x: [&NVector; 4] = [yn, ycur, fn_, tempv2];
            N_VLinearCombination(4, &c, &x, tempv1);
        }
        *dsmPtr = N_VWrmsNorm(&ark_mem.tempv1, &ark_mem.ewt);
        lsrkStep_DomEigUpdateLogic(ark_mem, step_mem, *dsmPtr);
    } else {
        lsrkStep_DomEigUpdateLogic(ark_mem, step_mem, *dsmPtr);
    }

    ARK_SUCCESS
}

/*===============================================================
  Internal utility routines
  ===============================================================*/

/*---------------------------------------------------------------
  lsrkStep_AccessStepMem:

  Shortcut routine to unpack the step_mem structure (take
  semantics; callers put it back).
  ---------------------------------------------------------------*/
pub(crate) fn lsrkStep_AccessStepMem(
    ark_mem: &mut ARKodeMem,
    fname: &str,
) -> Option<Box<ARKodeLSRKStepMem>> {
    match ark_mem.step_mem.take() {
        Some(b) => match b.downcast::<ARKodeLSRKStepMem>() {
            Ok(sm) => Some(sm),
            Err(other) => {
                ark_mem.step_mem = Some(other);
                arkProcessError(
                    Some(ark_mem),
                    ARK_MEM_NULL,
                    line!(),
                    fname,
                    file!(),
                    MSG_LSRKSTEP_NO_MEM,
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
                MSG_LSRKSTEP_NO_MEM,
            );
            None
        }
    }
}

/*---------------------------------------------------------------
  lsrkStep_Free frees all LSRKStep memory.
  ---------------------------------------------------------------*/
pub fn lsrkStep_Free(ark_mem: &mut ARKodeMem) {
    if let Some(step_mem) = lsrkStep_AccessStepMem(ark_mem, "lsrkStep_Free") {
        /* free the reusable arrays for fused vector interface */
        if !step_mem.cvals.is_empty() {
            ark_mem.lrw -= step_mem.nfusedopvecs as i64;
            ark_mem.liw -= step_mem.nfusedopvecs as i64;
        }
        /* the box (and any DEE inside) is dropped here */
    }
    ark_mem.step_mem = None;
}

/*---------------------------------------------------------------
  lsrkStep_DomEigUpdateLogic:

  This routine checks if the step is accepted or not and reassigns
  the dom_eig update flags accordingly.  (C receives fnew, which is
  always ark_mem->tempv2 at the call sites; the copy into fn is
  done through the fields here.)
  ---------------------------------------------------------------*/
pub(crate) fn lsrkStep_DomEigUpdateLogic(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeLSRKStepMem>,
    dsm: f64,
) {
    if dsm <= ONE {
        {
            let ARKodeMem { fn_, tempv2, .. } = ark_mem;
            N_VScale(ONE, tempv2, fn_);
        }
        ark_mem.fn_is_current = true;

        step_mem.dom_eig_is_current = step_mem.const_Jac;

        step_mem.dom_eig_update = false;
        if ark_mem.nst + 1 >= step_mem.dom_eig_nst + step_mem.dom_eig_freq {
            step_mem.dom_eig_update = !step_mem.dom_eig_is_current;
        }
    } else {
        step_mem.dom_eig_update = !step_mem.dom_eig_is_current;
    }
}

/*---------------------------------------------------------------
  lsrkStep_ComputeNewDomEig:

  This routine computes new dom_eig and returns SUN_SUCCESS.
  ---------------------------------------------------------------*/
pub(crate) fn lsrkStep_ComputeNewDomEig(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeLSRKStepMem>,
) -> i32 {
    let mut retval;

    if step_mem.DEE.is_some() {
        let mut dee = step_mem.DEE.take().unwrap();

        /* Build the ATimes closure over the integrator memory (C:
           lsrkStep_DQJtimes with A_data = arkode_mem, installed once
           by LSRKStepSetDomEigEstimator; here supplied per call) */
        let fe = step_mem.fe.unwrap();
        let step_nst = step_mem.step_nst;
        let ARKodeLSRKStepMem { nfeDQ, lambdaR, lambdaI, .. } = &mut **step_mem;
        {
            let arm = RefCell::new(&mut *ark_mem);
            let mut atimes = |v: &NVector, jv: &mut NVector| -> i32 {
                let mut guard = arm.borrow_mut();
                let ar: &mut ARKodeMem = &mut guard;
                lsrkStep_DQJtimes(ar, fe, step_nst, nfeDQ, v, jv)
            };
            retval = SUNDomEigEstimator_Estimate(&mut dee, &mut atimes, lambdaR, lambdaI);
        }
        step_mem.dom_eig_num_evals += 1;
        if retval != SUN_SUCCESS {
            step_mem.DEE = Some(dee);
            arkProcessError(
                Some(ark_mem),
                ARK_DEE_FAIL,
                line!(),
                "lsrkStep_ComputeNewDomEig",
                file!(),
                "SUNDomEigEstimator_Estimate failed",
            );
            return ARK_DEE_FAIL;
        }

        let mut num_iters: i64 = 0;
        retval = SUNDomEigEstimator_GetNumIters(&dee, &mut num_iters);
        if retval != SUN_SUCCESS {
            step_mem.DEE = Some(dee);
            arkProcessError(
                Some(ark_mem),
                ARK_DEE_FAIL,
                line!(),
                "lsrkStep_ComputeNewDomEig",
                file!(),
                "SUNDomEigEstimator_GetNumIters failed",
            );
            return ARK_DEE_FAIL;
        }
        step_mem.num_dee_iters += num_iters;

        /* After the first call to SUNDomEigEstimator_Estimate, the
           number of warmups is set to num_warmups; this allows the
           successive calls to use a different number of warmups. */
        if step_mem.init_warmup {
            retval = SUNDomEigEstimator_SetNumPreprocessIters(&mut dee, step_mem.num_warmups);
            if retval != SUN_SUCCESS {
                step_mem.DEE = Some(dee);
                arkProcessError(
                    Some(ark_mem),
                    ARK_DEE_FAIL,
                    line!(),
                    "lsrkStep_ComputeNewDomEig",
                    file!(),
                    "SUNDomEigEstimator_SetNumPreprocessIters failed",
                );
                return ARK_DEE_FAIL;
            }
            step_mem.init_warmup = false;
        }
        step_mem.DEE = Some(dee);
    } else if step_mem.dom_eig_fn.is_some() {
        let dom_eig_fn = step_mem.dom_eig_fn.unwrap();
        retval = {
            let ARKodeLSRKStepMem { lambdaR, lambdaI, .. } = &mut **step_mem;
            let ARKodeMem {
                yn,
                fn_,
                user_data,
                tempv1,
                tempv2,
                tempv3,
                tn,
                ..
            } = ark_mem;
            dom_eig_fn(*tn, yn, fn_, lambdaR, lambdaI, user_data, tempv1, tempv2, tempv3)
        };
        step_mem.dom_eig_num_evals += 1;
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_DOMEIG_FAIL,
                line!(),
                "lsrkStep_ComputeNewDomEig",
                file!(),
                "Unable to estimate the dominant eigenvalue",
            );
            return ARK_DOMEIG_FAIL;
        }
    } else {
        arkProcessError(
            Some(ark_mem),
            ARK_DOMEIG_FAIL,
            line!(),
            "lsrkStep_ComputeNewDomEig",
            file!(),
            "Unable to estimate the dominant eigenvalue: Either a user provided function or a SUNDomEigEstimator is required",
        );
        return ARK_DOMEIG_FAIL;
    }

    if step_mem.lambdaR * ark_mem.h > ZERO {
        arkProcessError(
            None,
            ARK_DOMEIG_FAIL,
            line!(),
            "lsrkStep_ComputeNewDomEig",
            file!(),
            "lambdaR*h must be nonpositive",
        );
        return ARK_DOMEIG_FAIL;
    } else if step_mem.lambdaR == 0.0 && SUNRabs(step_mem.lambdaI) > 0.0 {
        arkProcessError(
            None,
            ARK_DOMEIG_FAIL,
            line!(),
            "lsrkStep_ComputeNewDomEig",
            file!(),
            "DomEig cannot be purely imaginary",
        );
        return ARK_DOMEIG_FAIL;
    }

    step_mem.lambdaR *= step_mem.dom_eig_safety;
    step_mem.lambdaI *= step_mem.dom_eig_safety;
    step_mem.spectral_radius =
        SUNRsqrt(SUNSQR(step_mem.lambdaR) + SUNSQR(step_mem.lambdaI));

    step_mem.dom_eig_is_current = true;
    step_mem.dom_eig_nst = ark_mem.nst;

    step_mem.spectral_radius_max = if step_mem.spectral_radius > step_mem.spectral_radius_max {
        step_mem.spectral_radius
    } else {
        step_mem.spectral_radius_max
    };

    if step_mem.spectral_radius < step_mem.spectral_radius_min || ark_mem.nst == 0 {
        step_mem.spectral_radius_min = step_mem.spectral_radius;
    }

    step_mem.dom_eig_update = false;

    retval
}

/*---------------------------------------------------------------
  lsrkStep_DQJtimes:

  This routine generates a difference quotient approximation to
  the Jacobian-vector product f_y(t,y) * v. The approximation is
  Jv = [f(y + v*sig) - f(y)]/sig, where sig = 1 / ||v||_WRMS,
  i.e. the WRMS norm of v*sig is 1.  (The step_mem pieces arrive
  as arguments because this runs inside the DEE ATimes closure.)
  ---------------------------------------------------------------*/
fn lsrkStep_DQJtimes(
    ark_mem: &mut ARKodeMem,
    fe: ARKRhsFn,
    step_nst: i64,
    nfeDQ: &mut i64,
    v: &NVector,
    Jv: &mut NVector,
) -> i32 {
    let t = ark_mem.tn;

    /* Compute RHS function, if necessary. */
    if (!ark_mem.fn_is_current && ark_mem.initsetup) || (step_nst != ark_mem.nst) {
        /* call the user-supplied pre-RHS function (if supplied) */
        if let Some(pre_rhs) = ark_mem.PreRhsFn {
            let ARKodeMem { yn, user_data, .. } = ark_mem;
            let retval = pre_rhs(t, yn, user_data);
            if retval != 0 {
                return ARK_PRERHSFN_FAIL;
            }
        }

        let retval = {
            let ARKodeMem { yn, fn_, user_data, .. } = ark_mem;
            fe(t, yn, fn_, user_data)
        };
        *nfeDQ += 1;
        if retval != ARK_SUCCESS {
            return ARK_RHSFUNC_FAIL;
        }
        ark_mem.fn_is_current = true;
    }

    /* Initialize perturbation to 1/||v|| */
    let mut sig = ONE / N_VWrmsNorm(v, &ark_mem.ewt);
    let mut retval = 0;

    for _iter in 0..MAX_DQITERS {
        /* Set work = y + sig*v */
        {
            let ARKodeMem { yn, tempv3, .. } = ark_mem;
            N_VLinearSum(sig, v, ONE, yn, tempv3);
        }

        /* call the user-supplied pre-RHS function (if supplied) */
        if let Some(pre_rhs) = ark_mem.PreRhsFn {
            let ARKodeMem { tempv3, user_data, .. } = ark_mem;
            let retval = pre_rhs(t, tempv3, user_data);
            if retval != 0 {
                return ARK_PRERHSFN_FAIL;
            }
        }
        /* Set Jv = f(tn, y+sig*v) */
        retval = {
            let ARKodeMem { tempv3, user_data, .. } = ark_mem;
            fe(t, tempv3, Jv, user_data)
        };
        *nfeDQ += 1;
        if retval == 0 {
            break;
        }
        if retval < 0 {
            return -1;
        }

        /* If f failed recoverably, shrink sig and retry */
        sig *= 0.25;
    }

    /* If retval still isn't 0, return with a recoverable failure */
    if retval > 0 {
        return 1;
    }

    /* Replace Jv by (Jv - fn)/sig */
    let siginv = ONE / sig;
    Jv.linear_sum_with(siginv, -siginv, &ark_mem.fn_);

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  lsrkStep_TakeStepSSPs2:

  This routine performs a single SSPs2 step (with embedding).
  ---------------------------------------------------------------*/
pub fn lsrkStep_TakeStepSSPs2(ark_mem: &mut ARKodeMem, dsmPtr: &mut f64, nflagPtr: &mut i32) -> i32 {
    let mut step_mem = match lsrkStep_AccessStepMem(ark_mem, "lsrkStep_TakeStepSSPs2") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };
    let ret = lsrkStep_TakeStepSSPs2_inner(ark_mem, &mut step_mem, dsmPtr, nflagPtr);
    ark_mem.step_mem = Some(step_mem);
    ret
}

fn lsrkStep_TakeStepSSPs2_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeLSRKStepMem>,
    dsmPtr: &mut f64,
    nflagPtr: &mut i32,
) -> i32 {
    /* initialize algebraic solver convergence flag to success,
       temporal error estimate to zero */
    *nflagPtr = ARK_SUCCESS;
    *dsmPtr = ZERO;

    /* Initialize the current stage index */
    step_mem.istage = 0;

    /* Initialize method coefficients */
    let rs = step_mem.req_stages as f64;
    let sm1inv = ONE / (rs - ONE);
    let hsm1inv = ark_mem.h * sm1inv;
    let rsinv = ONE / rs;
    let hrsinv = ark_mem.h * rsinv;
    let (hbt1, hbt2, hbt3);

    /* Embedding coefficients differ when req_stages == 2 */
    if step_mem.req_stages == 2 {
        /* from https://doi.org/10.1016/j.cam.2022.114325 pg 5 */
        hbt1 = ark_mem.h * 0.694021459207626;
        hbt2 = ZERO;
        hbt3 = ark_mem.h - hbt1;
    } else {
        hbt1 = hrsinv * (ONE + rsinv);
        hbt2 = hrsinv;
        hbt3 = hrsinv * (ONE - rsinv);
    }

    /* The method is not FSAL.  Therefore, fn is computed at the
       beginning of the step unless the previous step failed or ARKODE
       updated fn. */
    if !ark_mem.fn_is_current {
        /* call the user-supplied pre-RHS function (if supplied) */
        if let Some(pre_rhs) = ark_mem.PreRhsFn {
            let ARKodeMem { yn, user_data, tn, .. } = ark_mem;
            let retval = pre_rhs(*tn, yn, user_data);
            if retval != 0 {
                return ARK_PRERHSFN_FAIL;
            }
        }

        let fe = step_mem.fe.unwrap();
        let retval = {
            let ARKodeMem { yn, fn_, user_data, tn, .. } = ark_mem;
            fe(*tn, yn, fn_, user_data)
        };
        step_mem.nfe += 1;
        if retval != ARK_SUCCESS {
            return ARK_RHSFUNC_FAIL;
        }
        ark_mem.fn_is_current = true;
    }

    /* Begin stage 1 and accumulate embedding into tempv1 */
    ark_mem.tcur = ark_mem.tn + hsm1inv;
    step_mem.istage = 1;
    {
        let ARKodeMem { yn, fn_, ycur, .. } = ark_mem;
        N_VLinearSum(ONE, yn, hsm1inv, fn_, ycur);
    }
    if !ark_mem.fixedstep {
        let ARKodeMem { yn, fn_, tempv1, .. } = ark_mem;
        N_VLinearSum(ONE, yn, hbt1, fn_, tempv1);
    }

    /* apply user-supplied stage postprocessing function (if supplied) */
    if let Some(post) = ark_mem.PostProcessStageFn {
        let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
        let retval = post(*tcur, ycur, user_data);
        if retval != 0 {
            return ARK_POSTPROCESS_STAGE_FAIL;
        }
    }

    /* Evaluate stages j = 2,...,step_mem->req_stages - 1 */
    for j in 2..step_mem.req_stages {
        /* Complete the previous stage (evaluate the RHS and store it in
           tempv2) */

        /* apply user-supplied stage preprocessing function (if supplied) */
        if let Some(pre_rhs) = ark_mem.PreRhsFn {
            let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
            let retval = pre_rhs(*tcur, ycur, user_data);
            if retval != 0 {
                return ARK_PRERHSFN_FAIL;
            }
        }

        let fe = step_mem.fe.unwrap();
        let retval = {
            let ARKodeMem { ycur, tempv2, user_data, tcur, .. } = ark_mem;
            fe(*tcur, ycur, tempv2, user_data)
        };
        step_mem.nfe += 1;

        if retval < 0 {
            return ARK_RHSFUNC_FAIL;
        }
        if retval > 0 {
            return RHSFUNC_RECVR;
        }

        /* Begin stage j (update the state and embedding) */
        step_mem.istage = j;
        ark_mem.tcur = ark_mem.tn + j as f64 * hsm1inv;
        {
            let ARKodeMem { ycur, tempv2, .. } = ark_mem;
            ycur.linear_sum_with(ONE, hsm1inv, tempv2);
        }
        if !ark_mem.fixedstep {
            let ARKodeMem { tempv1, tempv2, .. } = ark_mem;
            tempv1.linear_sum_with(ONE, hbt2, tempv2);
        }

        /* apply user-supplied stage postprocessing function (if supplied) */
        if let Some(post) = ark_mem.PostProcessStageFn {
            let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
            let retval = post(*tcur, ycur, user_data);
            if retval != 0 {
                return ARK_POSTPROCESS_STAGE_FAIL;
            }
        }
    }

    /* Complete the next-to-last stage by evaluating the RHS and storing
       it in tempv2 */
    if let Some(pre_rhs) = ark_mem.PreRhsFn {
        let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
        let retval = pre_rhs(*tcur, ycur, user_data);
        if retval != 0 {
            return ARK_PRERHSFN_FAIL;
        }
    }
    let fe = step_mem.fe.unwrap();
    let retval = {
        let ARKodeMem { ycur, tempv2, user_data, tcur, .. } = ark_mem;
        fe(*tcur, ycur, tempv2, user_data)
    };
    step_mem.nfe += 1;

    if retval < 0 {
        return ARK_RHSFUNC_FAIL;
    }
    if retval > 0 {
        return RHSFUNC_RECVR;
    }

    /* Compute the step solution */
    ark_mem.tcur = ark_mem.tn + ark_mem.h;
    step_mem.istage = step_mem.req_stages;
    {
        /* N_VLinearCombination(3, ...) with z == X[0] == ycur */
        let ARKodeMem { ycur, yn, tempv2, .. } = ark_mem;
        lsrk_lincomb_inplace(
            ycur,
            ONE / (sm1inv * rs),
            &[rsinv, hrsinv],
            &[yn, tempv2],
        );
    }

    /* apply user-supplied step postprocessing function (if supplied) */
    if let Some(post) = ark_mem.PostProcessStepFn {
        let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
        let retval = post(*tcur, ycur, user_data);
        if retval != 0 {
            return ARK_POSTPROCESS_STEP_FAIL;
        }
    }

    /* Compute yerr (if step adaptivity enabled) */
    if !ark_mem.fixedstep {
        {
            let ARKodeMem { tempv1, tempv2, .. } = ark_mem;
            tempv1.linear_sum_with(ONE, hbt3, tempv2);
        }
        {
            /* C VScaleDiff kernel: tempv1 = 1*(ycur - tempv1) */
            let ARKodeMem { ycur, tempv1, .. } = ark_mem;
            for e in 0..tempv1.data.len() {
                tempv1.data[e] = ycur.data[e] - tempv1.data[e];
            }
        }
        *dsmPtr = N_VWrmsNorm(&ark_mem.tempv1, &ark_mem.ewt);
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  lsrkStep_TakeStepSSPs3:

  This routine performs a single SSPs3 step (with embedding).
  The s = 4 case is handled separately by lsrkStep_TakeStepSSP43.
  ---------------------------------------------------------------*/
pub fn lsrkStep_TakeStepSSPs3(ark_mem: &mut ARKodeMem, dsmPtr: &mut f64, nflagPtr: &mut i32) -> i32 {
    let mut step_mem = match lsrkStep_AccessStepMem(ark_mem, "lsrkStep_TakeStepSSPs3") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };
    let ret = lsrkStep_TakeStepSSPs3_inner(ark_mem, &mut step_mem, dsmPtr, nflagPtr);
    ark_mem.step_mem = Some(step_mem);
    ret
}

#[allow(clippy::too_many_lines)]
fn lsrkStep_TakeStepSSPs3_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeLSRKStepMem>,
    dsmPtr: &mut f64,
    nflagPtr: &mut i32,
) -> i32 {
    use crate::sundials_math::SUNRround;

    /* initialize algebraic solver convergence flag to success,
       temporal error estimate to zero */
    *nflagPtr = ARK_SUCCESS;
    *dsmPtr = ZERO;

    /* Initialize the current stage index */
    step_mem.istage = 0;

    /* Initialize method coefficients */
    let rs = step_mem.req_stages as f64;
    let rn = SUNRsqrt(rs);
    let hrat = ark_mem.h / (rs - rn);
    let hrsinv = ark_mem.h / rs;
    let in_ = SUNRround(rn) as i32;

    /* The method is not FSAL.  Therefore, fn is computed at the
       beginning of the step unless ARKODE updated fn. */
    if !ark_mem.fn_is_current {
        /* call the user-supplied pre-RHS function (if supplied) */
        if let Some(pre_rhs) = ark_mem.PreRhsFn {
            let ARKodeMem { yn, user_data, tn, .. } = ark_mem;
            let retval = pre_rhs(*tn, yn, user_data);
            if retval != 0 {
                return ARK_PRERHSFN_FAIL;
            }
        }

        let fe = step_mem.fe.unwrap();
        let retval = {
            let ARKodeMem { yn, fn_, user_data, tn, .. } = ark_mem;
            fe(*tn, yn, fn_, user_data)
        };
        step_mem.nfe += 1;
        if retval != ARK_SUCCESS {
            return ARK_RHSFUNC_FAIL;
        }
        ark_mem.fn_is_current = true;
    }

    /* Begin stage 1 and accumulate embedding into tempv1 */
    ark_mem.tcur = ark_mem.tn + hrat;
    step_mem.istage = 1;
    {
        let ARKodeMem { yn, fn_, ycur, .. } = ark_mem;
        N_VLinearSum(ONE, yn, hrat, fn_, ycur);
    }
    if !ark_mem.fixedstep {
        let ARKodeMem { yn, fn_, tempv1, .. } = ark_mem;
        N_VLinearSum(ONE, yn, hrsinv, fn_, tempv1);
    }

    /* apply user-supplied stage postprocessing function (if supplied) */
    if let Some(post) = ark_mem.PostProcessStageFn {
        let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
        let retval = post(*tcur, ycur, user_data);
        if retval != 0 {
            return ARK_POSTPROCESS_STAGE_FAIL;
        }
    }

    /* Evaluate first stage group */
    for j in 2..=((in_ - 1) * (in_ - 2) / 2) {
        /* Complete the previous stage (evaluate the RHS and store it in
           tempv3) */
        if let Some(pre_rhs) = ark_mem.PreRhsFn {
            let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
            let retval = pre_rhs(*tcur, ycur, user_data);
            if retval != 0 {
                return ARK_PRERHSFN_FAIL;
            }
        }

        let fe = step_mem.fe.unwrap();
        let retval = {
            let ARKodeMem { ycur, tempv3, user_data, tcur, .. } = ark_mem;
            fe(*tcur, ycur, tempv3, user_data)
        };
        step_mem.nfe += 1;

        if retval < 0 {
            return ARK_RHSFUNC_FAIL;
        }
        if retval > 0 {
            return RHSFUNC_RECVR;
        }

        /* Begin stage j (update the state and embedding) */
        ark_mem.tcur = ark_mem.tn + j as f64 * hrat;
        step_mem.istage = j;
        {
            let ARKodeMem { ycur, tempv3, .. } = ark_mem;
            ycur.linear_sum_with(ONE, hrat, tempv3);
        }
        if !ark_mem.fixedstep {
            let ARKodeMem { tempv1, tempv3, .. } = ark_mem;
            tempv1.linear_sum_with(ONE, hrsinv, tempv3);
        }

        /* apply user-supplied stage postprocessing function (if supplied) */
        if let Some(post) = ark_mem.PostProcessStageFn {
            let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
            let retval = post(*tcur, ycur, user_data);
            if retval != 0 {
                return ARK_POSTPROCESS_STAGE_FAIL;
            }
        }
    }

    /* Copy ycur into tempv2 before looping over second stage group */
    {
        let ARKodeMem { ycur, tempv2, .. } = ark_mem;
        N_VScale(ONE, ycur, tempv2);
    }

    /* Evaluate second stage group */
    for j in ((in_ - 1) * (in_ - 2) / 2 + 1)..=(in_ * (in_ + 1) / 2 - 1) {
        /* Complete the previous stage (evaluate the RHS and store it in
           tempv3) */
        if let Some(pre_rhs) = ark_mem.PreRhsFn {
            let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
            let retval = pre_rhs(*tcur, ycur, user_data);
            if retval != 0 {
                return ARK_PRERHSFN_FAIL;
            }
        }

        let fe = step_mem.fe.unwrap();
        let retval = {
            let ARKodeMem { ycur, tempv3, user_data, tcur, .. } = ark_mem;
            fe(*tcur, ycur, tempv3, user_data)
        };
        step_mem.nfe += 1;

        if retval < 0 {
            return ARK_RHSFUNC_FAIL;
        }
        if retval > 0 {
            return RHSFUNC_RECVR;
        }

        /* Begin stage j (update the state and embedding) */
        ark_mem.tcur = ark_mem.tn + j as f64 * hrat;
        step_mem.istage = j;
        {
            let ARKodeMem { ycur, tempv3, .. } = ark_mem;
            ycur.linear_sum_with(ONE, hrat, tempv3);
        }
        if !ark_mem.fixedstep {
            let ARKodeMem { tempv1, tempv3, .. } = ark_mem;
            tempv1.linear_sum_with(ONE, hrsinv, tempv3);
        }

        /* apply user-supplied stage postprocessing function (if supplied) */
        if let Some(post) = ark_mem.PostProcessStageFn {
            let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
            let retval = post(*tcur, ycur, user_data);
            if retval != 0 {
                return ARK_POSTPROCESS_STAGE_FAIL;
            }
        }
    }

    /* apply user-supplied stage preprocessing function (if supplied) */
    if let Some(pre_rhs) = ark_mem.PreRhsFn {
        let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
        let retval = pre_rhs(*tcur, ycur, user_data);
        if retval != 0 {
            return ARK_PRERHSFN_FAIL;
        }
    }

    let fe = step_mem.fe.unwrap();
    let retval = {
        let ARKodeMem { ycur, tempv3, user_data, tcur, .. } = ark_mem;
        fe(*tcur, ycur, tempv3, user_data)
    };
    step_mem.nfe += 1;

    if retval < 0 {
        return ARK_RHSFUNC_FAIL;
    }
    if retval > 0 {
        return RHSFUNC_RECVR;
    }

    /* Begin the next stage before final stage group */
    ark_mem.tcur = ark_mem.tn + (in_ * (in_ - 1) / 2) as f64 * hrat;
    step_mem.istage = in_ * (in_ + 1) / 2;
    {
        /* N_VLinearCombination(3, ...) with z == X[0] == ycur */
        let ARKodeMem { ycur, tempv2, tempv3, .. } = ark_mem;
        lsrk_lincomb_inplace(
            ycur,
            (rn - ONE) / (TWO * rn - ONE),
            &[rn / (TWO * rn - ONE), (rn - ONE) * hrat / (TWO * rn - ONE)],
            &[tempv2, tempv3],
        );
    }
    if !ark_mem.fixedstep {
        let ARKodeMem { tempv1, tempv3, .. } = ark_mem;
        tempv1.linear_sum_with(ONE, hrsinv, tempv3);
    }

    /* apply user-supplied stage postprocessing function (if supplied) */
    if let Some(post) = ark_mem.PostProcessStageFn {
        let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
        let retval = post(*tcur, ycur, user_data);
        if retval != 0 {
            return ARK_POSTPROCESS_STAGE_FAIL;
        }
    }

    /* Evaluate final stage group */
    for j in (in_ * (in_ + 1) / 2 + 1)..=step_mem.req_stages {
        /* Complete the previous stage (evaluate the RHS and store it in
           tempv3) */
        if let Some(pre_rhs) = ark_mem.PreRhsFn {
            let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
            let retval = pre_rhs(*tcur, ycur, user_data);
            if retval != 0 {
                return ARK_PRERHSFN_FAIL;
            }
        }

        let fe = step_mem.fe.unwrap();
        let retval = {
            let ARKodeMem { ycur, tempv3, user_data, tcur, .. } = ark_mem;
            fe(*tcur, ycur, tempv3, user_data)
        };
        step_mem.nfe += 1;

        if retval < 0 {
            return ARK_RHSFUNC_FAIL;
        }
        if retval > 0 {
            return RHSFUNC_RECVR;
        }

        /* Begin stage j (update the state and embedding) */
        ark_mem.tcur = ark_mem.tn + (j - in_) as f64 * hrat;
        step_mem.istage = j;
        {
            let ARKodeMem { ycur, tempv3, .. } = ark_mem;
            ycur.linear_sum_with(ONE, hrat, tempv3);
        }
        if !ark_mem.fixedstep {
            let ARKodeMem { tempv1, tempv3, .. } = ark_mem;
            tempv1.linear_sum_with(ONE, hrsinv, tempv3);
        }

        /* apply user-supplied stage or step postprocessing function (if
           supplied) */
        if let (true, Some(post)) = (j < step_mem.req_stages, ark_mem.PostProcessStageFn) {
            let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
            let retval = post(*tcur, ycur, user_data);
            if retval != 0 {
                return ARK_POSTPROCESS_STAGE_FAIL;
            }
        } else if let (true, Some(post)) = (j == step_mem.req_stages, ark_mem.PostProcessStepFn) {
            let tpost = ark_mem.tn + ark_mem.h;
            let ARKodeMem { ycur, user_data, .. } = ark_mem;
            let retval = post(tpost, ycur, user_data);
            if retval != 0 {
                return ARK_POSTPROCESS_STEP_FAIL;
            }
        }
    }

    /* Compute yerr (if step adaptivity enabled) */
    if !ark_mem.fixedstep {
        {
            /* C VScaleDiff kernel: tempv1 = 1*(ycur - tempv1) */
            let ARKodeMem { ycur, tempv1, .. } = ark_mem;
            for e in 0..tempv1.data.len() {
                tempv1.data[e] = ycur.data[e] - tempv1.data[e];
            }
        }
        *dsmPtr = N_VWrmsNorm(&ark_mem.tempv1, &ark_mem.ewt);
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  lsrkStep_TakeStepSSP43:

  This routine performs a single SSP43 step (with embedding).
  ---------------------------------------------------------------*/
pub fn lsrkStep_TakeStepSSP43(ark_mem: &mut ARKodeMem, dsmPtr: &mut f64, nflagPtr: &mut i32) -> i32 {
    let mut step_mem = match lsrkStep_AccessStepMem(ark_mem, "lsrkStep_TakeStepSSP43") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };
    let ret = lsrkStep_TakeStepSSP43_inner(ark_mem, &mut step_mem, dsmPtr, nflagPtr);
    ark_mem.step_mem = Some(step_mem);
    ret
}

fn lsrkStep_TakeStepSSP43_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeLSRKStepMem>,
    dsmPtr: &mut f64,
    nflagPtr: &mut i32,
) -> i32 {
    /* initialize algebraic solver convergence flag to success,
       temporal error estimate to zero */
    *nflagPtr = ARK_SUCCESS;
    *dsmPtr = ZERO;

    /* Initialize the current stage index */
    step_mem.istage = 0;

    /* Initialize method coefficients */
    let rs: f64 = 4.0;
    let hp5 = ark_mem.h * 0.5;
    let hrsinv = ark_mem.h / rs;

    /* The method is not FSAL.  Therefore, fn is computed at the
       beginning of the step unless ARKODE updated fn. */
    if !ark_mem.fn_is_current {
        if let Some(pre_rhs) = ark_mem.PreRhsFn {
            let ARKodeMem { yn, user_data, tn, .. } = ark_mem;
            let retval = pre_rhs(*tn, yn, user_data);
            if retval != 0 {
                return ARK_PRERHSFN_FAIL;
            }
        }

        let fe = step_mem.fe.unwrap();
        let retval = {
            let ARKodeMem { yn, fn_, user_data, tn, .. } = ark_mem;
            fe(*tn, yn, fn_, user_data)
        };
        step_mem.nfe += 1;
        if retval != ARK_SUCCESS {
            return ARK_RHSFUNC_FAIL;
        }
        ark_mem.fn_is_current = true;
    }

    /* Begin stage 1 and accumulate embedding into tempv1 */
    ark_mem.tcur = ark_mem.tn + hp5;
    step_mem.istage = 1;
    {
        let ARKodeMem { yn, fn_, ycur, .. } = ark_mem;
        N_VLinearSum(ONE, yn, hp5, fn_, ycur);
    }
    if !ark_mem.fixedstep {
        let ARKodeMem { yn, fn_, tempv1, .. } = ark_mem;
        N_VLinearSum(ONE, yn, hrsinv, fn_, tempv1);
    }

    /* apply user-supplied stage postprocessing function (if supplied) */
    if let Some(post) = ark_mem.PostProcessStageFn {
        let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
        let retval = post(*tcur, ycur, user_data);
        if retval != 0 {
            return ARK_POSTPROCESS_STAGE_FAIL;
        }
    }

    /* call the user-supplied pre-RHS function (if supplied) */
    if let Some(pre_rhs) = ark_mem.PreRhsFn {
        let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
        let retval = pre_rhs(*tcur, ycur, user_data);
        if retval != 0 {
            return ARK_PRERHSFN_FAIL;
        }
    }

    /* Evaluate stage RHS */
    let fe = step_mem.fe.unwrap();
    let retval = {
        let ARKodeMem { ycur, tempv3, user_data, tcur, .. } = ark_mem;
        fe(*tcur, ycur, tempv3, user_data)
    };
    step_mem.nfe += 1;

    if retval < 0 {
        return ARK_RHSFUNC_FAIL;
    }
    if retval > 0 {
        return RHSFUNC_RECVR;
    }

    /* Begin stage 2 */
    ark_mem.tcur = ark_mem.tn + ark_mem.h;
    step_mem.istage = 2;
    {
        let ARKodeMem { ycur, tempv3, .. } = ark_mem;
        ycur.linear_sum_with(ONE, hp5, tempv3);
    }
    if !ark_mem.fixedstep {
        let ARKodeMem { tempv1, tempv3, .. } = ark_mem;
        tempv1.linear_sum_with(ONE, hrsinv, tempv3);
    }

    /* apply user-supplied stage postprocessing function (if supplied) */
    if let Some(post) = ark_mem.PostProcessStageFn {
        let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
        let retval = post(*tcur, ycur, user_data);
        if retval != 0 {
            return ARK_POSTPROCESS_STAGE_FAIL;
        }
    }

    /* Evaluate stage RHS */
    if let Some(pre_rhs) = ark_mem.PreRhsFn {
        let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
        let retval = pre_rhs(*tcur, ycur, user_data);
        if retval != 0 {
            return ARK_PRERHSFN_FAIL;
        }
    }

    let retval = {
        let ARKodeMem { ycur, tempv3, user_data, tcur, .. } = ark_mem;
        fe(*tcur, ycur, tempv3, user_data)
    };
    step_mem.nfe += 1;

    if retval < 0 {
        return ARK_RHSFUNC_FAIL;
    }
    if retval > 0 {
        return RHSFUNC_RECVR;
    }

    /* Begin stage 3 */
    ark_mem.tcur = ark_mem.tn + hp5;
    step_mem.istage = 3;
    {
        /* N_VLinearCombination(3, ...) with z == X[0] == ycur */
        let ARKodeMem { ycur, yn, tempv3, h, .. } = ark_mem;
        lsrk_lincomb_inplace(
            ycur,
            ONE / THREE,
            &[TWO / THREE, ONE / SIX * *h],
            &[yn, tempv3],
        );
    }
    if !ark_mem.fixedstep {
        let ARKodeMem { tempv1, tempv3, .. } = ark_mem;
        tempv1.linear_sum_with(ONE, hrsinv, tempv3);
    }

    /* apply user-supplied stage postprocessing function (if supplied) */
    if let Some(post) = ark_mem.PostProcessStageFn {
        let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
        let retval = post(*tcur, ycur, user_data);
        if retval != 0 {
            return ARK_POSTPROCESS_STAGE_FAIL;
        }
    }

    /* Evaluate stage RHS */
    if let Some(pre_rhs) = ark_mem.PreRhsFn {
        let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
        let retval = pre_rhs(*tcur, ycur, user_data);
        if retval != 0 {
            return ARK_PRERHSFN_FAIL;
        }
    }

    let retval = {
        let ARKodeMem { ycur, tempv3, user_data, tcur, .. } = ark_mem;
        fe(*tcur, ycur, tempv3, user_data)
    };
    step_mem.nfe += 1;

    if retval < 0 {
        return ARK_RHSFUNC_FAIL;
    }
    if retval > 0 {
        return RHSFUNC_RECVR;
    }

    /* Compute the time step solution and embedding */
    ark_mem.tcur = ark_mem.tn + ark_mem.h;
    step_mem.istage = 4;
    {
        let ARKodeMem { ycur, tempv3, .. } = ark_mem;
        ycur.linear_sum_with(ONE, hp5, tempv3);
    }

    /* apply user-supplied step postprocessing function (if supplied) */
    if let Some(post) = ark_mem.PostProcessStepFn {
        let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
        let retval = post(*tcur, ycur, user_data);
        if retval != 0 {
            return ARK_POSTPROCESS_STEP_FAIL;
        }
    }

    /* Compute yerr (if step adaptivity enabled) */
    if !ark_mem.fixedstep {
        {
            let ARKodeMem { tempv1, tempv3, .. } = ark_mem;
            tempv1.linear_sum_with(ONE, hrsinv, tempv3);
        }
        {
            /* C VScaleDiff kernel: tempv1 = 1*(ycur - tempv1) */
            let ARKodeMem { ycur, tempv1, .. } = ark_mem;
            for e in 0..tempv1.data.len() {
                tempv1.data[e] = ycur.data[e] - tempv1.data[e];
            }
        }
        *dsmPtr = N_VWrmsNorm(&ark_mem.tempv1, &ark_mem.ewt);
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  lsrkStep_TakeStepSSP104:

  This routine performs a single SSP104 step (with embedding).
  ---------------------------------------------------------------*/
pub fn lsrkStep_TakeStepSSP104(
    ark_mem: &mut ARKodeMem,
    dsmPtr: &mut f64,
    nflagPtr: &mut i32,
) -> i32 {
    let mut step_mem = match lsrkStep_AccessStepMem(ark_mem, "lsrkStep_TakeStepSSP104") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };
    let ret = lsrkStep_TakeStepSSP104_inner(ark_mem, &mut step_mem, dsmPtr, nflagPtr);
    ark_mem.step_mem = Some(step_mem);
    ret
}

#[allow(clippy::too_many_lines)]
fn lsrkStep_TakeStepSSP104_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeLSRKStepMem>,
    dsmPtr: &mut f64,
    nflagPtr: &mut i32,
) -> i32 {
    /* initialize algebraic solver convergence flag to success,
       temporal error estimate to zero */
    *nflagPtr = ARK_SUCCESS;
    *dsmPtr = ZERO;

    /* Initialize the current stage index */
    step_mem.istage = 0;

    /* Initialize method coefficients */
    let hsixth = ark_mem.h / SIX;
    let hfifth = ark_mem.h / FIVE;

    /* The method is not FSAL.  Therefore, fn is computed at the
       beginning of the step unless ARKODE updated fn. */
    if !ark_mem.fn_is_current {
        if let Some(pre_rhs) = ark_mem.PreRhsFn {
            let ARKodeMem { yn, user_data, tn, .. } = ark_mem;
            let retval = pre_rhs(*tn, yn, user_data);
            if retval != 0 {
                return ARK_PRERHSFN_FAIL;
            }
        }

        let fe = step_mem.fe.unwrap();
        let retval = {
            let ARKodeMem { yn, fn_, user_data, tn, .. } = ark_mem;
            fe(*tn, yn, fn_, user_data)
        };
        step_mem.nfe += 1;
        if retval != ARK_SUCCESS {
            return ARK_RHSFUNC_FAIL;
        }
        ark_mem.fn_is_current = true;
    }

    /* Copy yn into tempv2 for use in later stages */
    {
        let ARKodeMem { yn, tempv2, .. } = ark_mem;
        N_VScale(ONE, yn, tempv2);
    }

    /* Begin stage 1 and accumulate embedding into tempv1 */
    ark_mem.tcur = ark_mem.tn + hsixth;
    step_mem.istage = 1;
    {
        let ARKodeMem { yn, fn_, ycur, .. } = ark_mem;
        N_VLinearSum(ONE, yn, hsixth, fn_, ycur);
    }
    if !ark_mem.fixedstep {
        let ARKodeMem { yn, fn_, tempv1, .. } = ark_mem;
        N_VLinearSum(ONE, yn, hfifth, fn_, tempv1);
    }

    /* Evaluate stages j = 2,...,5 */
    for j in 2..=5 {
        /* Complete the previous stage (postprocesses the stage, evaluate
           the RHS, and store it in tempv3) */

        /* apply user-supplied stage postprocessing function (if supplied) */
        if let Some(post) = ark_mem.PostProcessStageFn {
            let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
            let retval = post(*tcur, ycur, user_data);
            if retval != 0 {
                return ARK_POSTPROCESS_STAGE_FAIL;
            }
        }

        /* apply user-supplied stage preprocessing function (if supplied) */
        if let Some(pre_rhs) = ark_mem.PreRhsFn {
            let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
            let retval = pre_rhs(*tcur, ycur, user_data);
            if retval != 0 {
                return ARK_PRERHSFN_FAIL;
            }
        }

        let fe = step_mem.fe.unwrap();
        let retval = {
            let ARKodeMem { ycur, tempv3, user_data, tcur, .. } = ark_mem;
            fe(*tcur, ycur, tempv3, user_data)
        };
        step_mem.nfe += 1;

        if retval < 0 {
            return ARK_RHSFUNC_FAIL;
        }
        if retval > 0 {
            return RHSFUNC_RECVR;
        }

        /* Begin stage j (update the state and embedding) */
        if j == 5 {
            ark_mem.tcur = ark_mem.tn + 2.0 * hsixth;
        } else {
            ark_mem.tcur = ark_mem.tn + j as f64 * hsixth;
        }
        step_mem.istage = j;
        {
            let ARKodeMem { ycur, tempv3, .. } = ark_mem;
            ycur.linear_sum_with(ONE, hsixth, tempv3);
        }
        if j == 4 && !ark_mem.fixedstep {
            let ARKodeMem { tempv1, tempv3, h, .. } = ark_mem;
            tempv1.linear_sum_with(ONE, 0.3 * *h, tempv3);
        }
    }

    /* no need to call RHS preprocessing here, since the stage does not
       require a RHS function evaluation */

    /* Finish stage 5 by preparing for the final stage group */
    {
        let ARKodeMem { ycur, tempv2, .. } = ark_mem;
        /* tempv2 = (1/25)*tempv2 + (9/25)*ycur (z == x aliased) */
        tempv2.linear_sum_with(1.0 / 25.0, 9.0 / 25.0, ycur);
    }
    {
        let ARKodeMem { ycur, tempv2, .. } = ark_mem;
        /* ycur = 15*tempv2 - 5*ycur (z == y aliased) */
        ycur.linear_sum_with(-5.0, 15.0, tempv2);
    }

    /* apply user-supplied stage postprocessing function (if supplied) */
    if let Some(post) = ark_mem.PostProcessStageFn {
        let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
        let retval = post(*tcur, ycur, user_data);
        if retval != 0 {
            return ARK_POSTPROCESS_STAGE_FAIL;
        }
    }

    /* Evaluate stages j = 6,...,9 */
    for j in 6..=9 {
        /* Complete the previous stage (evaluate the RHS and store in
           tempv3) */
        if let Some(pre_rhs) = ark_mem.PreRhsFn {
            let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
            let retval = pre_rhs(*tcur, ycur, user_data);
            if retval != 0 {
                return ARK_PRERHSFN_FAIL;
            }
        }

        let fe = step_mem.fe.unwrap();
        let retval = {
            let ARKodeMem { ycur, tempv3, user_data, tcur, .. } = ark_mem;
            fe(*tcur, ycur, tempv3, user_data)
        };
        step_mem.nfe += 1;

        if retval < 0 {
            return ARK_RHSFUNC_FAIL;
        }
        if retval > 0 {
            return RHSFUNC_RECVR;
        }

        /* Begin stage j (update the state and embedding) */
        ark_mem.tcur = ark_mem.tn + (j - 3) as f64 * hsixth;
        step_mem.istage = j;
        {
            let ARKodeMem { ycur, tempv3, .. } = ark_mem;
            ycur.linear_sum_with(ONE, hsixth, tempv3);
        }

        if j == 7 && !ark_mem.fixedstep {
            let ARKodeMem { tempv1, tempv3, .. } = ark_mem;
            tempv1.linear_sum_with(ONE, hfifth, tempv3);
        }
        if j == 9 && !ark_mem.fixedstep {
            let ARKodeMem { tempv1, tempv3, h, .. } = ark_mem;
            tempv1.linear_sum_with(ONE, 0.3 * *h, tempv3);
        }

        /* apply user-supplied stage postprocessing function (if supplied) */
        if let Some(post) = ark_mem.PostProcessStageFn {
            let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
            let retval = post(*tcur, ycur, user_data);
            if retval != 0 {
                return ARK_POSTPROCESS_STAGE_FAIL;
            }
        }
    }

    /* Complete the previous stage (evaluate the RHS and store it in
       tempv3) */
    if let Some(pre_rhs) = ark_mem.PreRhsFn {
        let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
        let retval = pre_rhs(*tcur, ycur, user_data);
        if retval != 0 {
            return ARK_PRERHSFN_FAIL;
        }
    }

    let fe = step_mem.fe.unwrap();
    let retval = {
        let ARKodeMem { ycur, tempv3, user_data, tcur, .. } = ark_mem;
        fe(*tcur, ycur, tempv3, user_data)
    };
    step_mem.nfe += 1;

    if retval < 0 {
        return ARK_RHSFUNC_FAIL;
    }
    if retval > 0 {
        return RHSFUNC_RECVR;
    }

    /* Compute the final time step solution */
    step_mem.istage = 10;
    {
        /* N_VLinearCombination(3, ...) with z == X[0] == ycur */
        let ARKodeMem { ycur, tempv2, tempv3, h, .. } = ark_mem;
        lsrk_lincomb_inplace(ycur, 0.6, &[ONE, 0.1 * *h], &[tempv2, tempv3]);
    }

    /* apply user-supplied step postprocessing function (if supplied) */
    if let Some(post) = ark_mem.PostProcessStepFn {
        let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
        let retval = post(*tcur, ycur, user_data);
        if retval != 0 {
            return ARK_POSTPROCESS_STEP_FAIL;
        }
    }

    /* Compute yerr (if step adaptivity enabled) */
    if !ark_mem.fixedstep {
        {
            /* C VScaleDiff kernel: tempv1 = 1*(ycur - tempv1) */
            let ARKodeMem { ycur, tempv1, .. } = ark_mem;
            for e in 0..tempv1.data.len() {
                tempv1.data[e] = ycur.data[e] - tempv1.data[e];
            }
        }
        *dsmPtr = N_VWrmsNorm(&ark_mem.tempv1, &ark_mem.ewt);
    }

    ARK_SUCCESS
}
