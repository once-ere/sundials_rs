/* -----------------------------------------------------------------
 * Translation of src/arkode/arkode_lsrkstep_io.c (ARKODE 7.7.0).
 * LSRKStep optional input/output functions.
 *
 * Not ported: lsrkStep_SetOptions (the ARKODE CLI dispatcher is a
 * separate pending module; the LSRK-specific keys land with it).
 * ---------------------------------------------------------------*/

use crate::arkode_impl::*;
use crate::arkode_io::{sunfprintf_long, sunfprintf_real};
use crate::arkode_lsrkstep::{
    lsrkStep_AccessStepMem, lsrkStep_TakeStepRKC, lsrkStep_TakeStepRKL, lsrkStep_TakeStepSSP104,
    lsrkStep_TakeStepSSP43, lsrkStep_TakeStepSSPs2, lsrkStep_TakeStepSSPs3,
};
use crate::arkode_lsrkstep_impl::*;
use crate::nvector_serial::{N_VScale, NVector};
use crate::sundials_domeigestimator::{
    SUNDomEigEstimator, SUNDomEigEstimator_SetNumPreprocessIters,
};
use crate::sundials_errors::SUN_SUCCESS;
use crate::sundials_types::{SUNOutputFormat, SUNFALSE};

/*===============================================================
  Exported optional input functions.
  ===============================================================*/

/*---------------------------------------------------------------
  LSRKStepSetSTSMethod sets method
    ARKODE_LSRK_RKC_2
    ARKODE_LSRK_RKL_2
  ---------------------------------------------------------------*/
pub fn LSRKStepSetSTSMethod(ark_mem: &mut ARKodeMem, method: ARKODE_LSRKMethodType) -> i32 {
    let mut step_mem = match lsrkStep_AccessStepMem(ark_mem, "LSRKStepSetSTSMethod") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    match method {
        ARKODE_LSRK_RKC_2 => {
            ark_mem.step = Some(lsrkStep_TakeStepRKC);
            step_mem.is_SSP = false;
            step_mem.nfusedopvecs = 5;
            step_mem.q = 2;
            ark_mem.hadapt_mem.as_mut().unwrap().q = 2;
            step_mem.p = 2;
            ark_mem.hadapt_mem.as_mut().unwrap().p = 2;
            step_mem.step_nst = 0;
        }
        ARKODE_LSRK_RKL_2 => {
            ark_mem.step = Some(lsrkStep_TakeStepRKL);
            step_mem.is_SSP = false;
            step_mem.nfusedopvecs = 5;
            step_mem.q = 2;
            ark_mem.hadapt_mem.as_mut().unwrap().q = 2;
            step_mem.p = 2;
            ark_mem.hadapt_mem.as_mut().unwrap().p = 2;
            step_mem.step_nst = 0;
        }
        ARKODE_LSRK_SSP_S_2 | ARKODE_LSRK_SSP_S_3 | ARKODE_LSRK_SSP_10_4 => {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "LSRKStepSetSTSMethod",
                file!(),
                "Invalid method option: Call LSRKStepCreateSSP to create an SSP method first.",
            );
        }
        _ => {
            ark_mem.step_mem = Some(step_mem);
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "LSRKStepSetSTSMethod",
                file!(),
                "Invalid method option.",
            );
            return ARK_ILL_INPUT;
        }
    }

    step_mem.LSRKmethod = method;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  LSRKStepSetSSPMethod sets method
    ARKODE_LSRK_SSP_S_2
    ARKODE_LSRK_SSP_S_3
    ARKODE_LSRK_SSP_10_4
  ---------------------------------------------------------------*/
pub fn LSRKStepSetSSPMethod(ark_mem: &mut ARKodeMem, method: ARKODE_LSRKMethodType) -> i32 {
    let mut step_mem = match lsrkStep_AccessStepMem(ark_mem, "LSRKStepSetSSPMethod") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    match method {
        ARKODE_LSRK_RKC_2 | ARKODE_LSRK_RKL_2 => {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "LSRKStepSetSSPMethod",
                file!(),
                "Invalid method option: Call LSRKStepCreateSTS to create an STS method first.",
            );
        }
        ARKODE_LSRK_SSP_S_2 => {
            ark_mem.step = Some(lsrkStep_TakeStepSSPs2);
            step_mem.is_SSP = true;
            step_mem.req_stages = 2;
            step_mem.nfusedopvecs = 3;
            step_mem.q = 2;
            ark_mem.hadapt_mem.as_mut().unwrap().q = 2;
            step_mem.p = 1;
            ark_mem.hadapt_mem.as_mut().unwrap().p = 1;
        }
        ARKODE_LSRK_SSP_S_3 => {
            ark_mem.step = Some(lsrkStep_TakeStepSSP43);
            step_mem.is_SSP = true;
            step_mem.req_stages = 4;
            step_mem.nfusedopvecs = 3;
            step_mem.q = 3;
            ark_mem.hadapt_mem.as_mut().unwrap().q = 3;
            step_mem.p = 2;
            ark_mem.hadapt_mem.as_mut().unwrap().p = 2;
        }
        ARKODE_LSRK_SSP_10_4 => {
            ark_mem.step = Some(lsrkStep_TakeStepSSP104);
            step_mem.is_SSP = true;
            step_mem.req_stages = 10;
            step_mem.nfusedopvecs = 3;
            step_mem.q = 4;
            ark_mem.hadapt_mem.as_mut().unwrap().q = 4;
            step_mem.p = 3;
            ark_mem.hadapt_mem.as_mut().unwrap().p = 3;
        }
        _ => {
            ark_mem.step_mem = Some(step_mem);
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "LSRKStepSetSSPMethod",
                file!(),
                "Invalid method option.",
            );
            return ARK_ILL_INPUT;
        }
    }

    step_mem.LSRKmethod = method;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

pub fn LSRKStepSetSTSMethodByName(ark_mem: &mut ARKodeMem, emethod: &str) -> i32 {
    if emethod == "ARKODE_LSRK_RKC_2" {
        return LSRKStepSetSTSMethod(ark_mem, ARKODE_LSRK_RKC_2);
    }
    if emethod == "ARKODE_LSRK_RKL_2" {
        return LSRKStepSetSTSMethod(ark_mem, ARKODE_LSRK_RKL_2);
    }
    if emethod == "ARKODE_LSRK_SSP_S_2"
        || emethod == "ARKODE_LSRK_SSP_S_3"
        || emethod == "ARKODE_LSRK_SSP_10_4"
    {
        arkProcessError(
            None,
            ARK_ILL_INPUT,
            line!(),
            "LSRKStepSetSTSMethodByName",
            file!(),
            "Invalid method option: Call LSRKStepCreateSTS to create an STS method first.",
        );
    }

    arkProcessError(
        None,
        ARK_ILL_INPUT,
        line!(),
        "LSRKStepSetSTSMethodByName",
        file!(),
        "Unknown method type",
    );

    ARK_ILL_INPUT
}

pub fn LSRKStepSetSSPMethodByName(ark_mem: &mut ARKodeMem, emethod: &str) -> i32 {
    if emethod == "ARKODE_LSRK_RKC_2" || emethod == "ARKODE_LSRK_RKL_2" {
        arkProcessError(
            None,
            ARK_ILL_INPUT,
            line!(),
            "LSRKStepSetSSPMethodByName",
            file!(),
            "Invalid method option: Call LSRKStepCreateSSP to create an SSP method first.",
        );
    }
    if emethod == "ARKODE_LSRK_SSP_S_2" {
        return LSRKStepSetSSPMethod(ark_mem, ARKODE_LSRK_SSP_S_2);
    }
    if emethod == "ARKODE_LSRK_SSP_S_3" {
        return LSRKStepSetSSPMethod(ark_mem, ARKODE_LSRK_SSP_S_3);
    }
    if emethod == "ARKODE_LSRK_SSP_10_4" {
        return LSRKStepSetSSPMethod(ark_mem, ARKODE_LSRK_SSP_10_4);
    }

    arkProcessError(
        None,
        ARK_ILL_INPUT,
        line!(),
        "LSRKStepSetSSPMethodByName",
        file!(),
        "Unknown method type",
    );

    ARK_ILL_INPUT
}

/*---------------------------------------------------------------
  LSRKStepSetDomEigFn specifies the dom_eig function.
  ---------------------------------------------------------------*/
pub fn LSRKStepSetDomEigFn(ark_mem: &mut ARKodeMem, dom_eig: Option<ARKDomEigFn>) -> i32 {
    let mut step_mem = match lsrkStep_AccessStepMem(ark_mem, "LSRKStepSetDomEigFn") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* set the dom_eig routine pointer, and update relevant flags */
    step_mem.dom_eig_fn = dom_eig;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  LSRKStepSetDomEigFrequency sets dom_eig computation frequency -
  the dominant eigenvalue is recomputed after "nsteps" successful
  steps.

  nsteps = 0 refers to constant dominant eigenvalue
  nsteps < 0 resets the default value 25 and sets nonconstant
             dominant eigenvalue
  ---------------------------------------------------------------*/
pub fn LSRKStepSetDomEigFrequency(ark_mem: &mut ARKodeMem, nsteps: i64) -> i32 {
    let mut step_mem = match lsrkStep_AccessStepMem(ark_mem, "LSRKStepSetDomEigFrequency") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    if nsteps < 0 {
        step_mem.dom_eig_freq = DOM_EIG_FREQ_DEFAULT;
        step_mem.const_Jac = false;
    }

    if nsteps == 0 {
        step_mem.const_Jac = true;
        step_mem.dom_eig_freq = 1;
    } else {
        step_mem.dom_eig_freq = nsteps;
        step_mem.const_Jac = false;
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  LSRKStepSetMaxNumStages sets the maximum number of stages
  allowed.
  ---------------------------------------------------------------*/
pub fn LSRKStepSetMaxNumStages(ark_mem: &mut ARKodeMem, stage_max_limit: i32) -> i32 {
    let mut step_mem = match lsrkStep_AccessStepMem(ark_mem, "LSRKStepSetMaxNumStages") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    if stage_max_limit < 2 {
        step_mem.stage_max_limit = STAGE_MAX_LIMIT_DEFAULT;
    } else {
        step_mem.stage_max_limit = stage_max_limit;
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  LSRKStepSetDomEigSafetyFactor sets the safety factor for the
  DomEigs.  Calling this function with dom_eig_safety < 1 resets
  the default value.
  ---------------------------------------------------------------*/
pub fn LSRKStepSetDomEigSafetyFactor(ark_mem: &mut ARKodeMem, dom_eig_safety: f64) -> i32 {
    let mut step_mem = match lsrkStep_AccessStepMem(ark_mem, "LSRKStepSetDomEigSafetyFactor") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    if dom_eig_safety < ONE {
        step_mem.dom_eig_safety = DOM_EIG_SAFETY_DEFAULT;
    } else {
        step_mem.dom_eig_safety = dom_eig_safety;
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  LSRKStepSetNumDomEigEstInitPreprocessIters sets the number of
  the preprocessing iterations before the very first estimate
  call.
  ---------------------------------------------------------------*/
pub fn LSRKStepSetNumDomEigEstInitPreprocessIters(ark_mem: &mut ARKodeMem, num_iters: i32) -> i32 {
    let mut step_mem =
        match lsrkStep_AccessStepMem(ark_mem, "LSRKStepSetNumDomEigEstInitPreprocessIters") {
            None => return ARK_MEM_NULL,
            Some(sm) => sm,
        };

    /* This value will be used in lsrkStep_Init to set the number of
       preprocessing iterations for the first dominant eigenvalue
       estimate.  If num_iters < 0, then the DEE's default will be
       used. */
    step_mem.num_init_warmups = num_iters;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  LSRKStepSetNumDomEigEstPreprocessIters sets the number of the
  preprocessing iterations before each estimate call after the
  initial estimate call.
  ---------------------------------------------------------------*/
pub fn LSRKStepSetNumDomEigEstPreprocessIters(ark_mem: &mut ARKodeMem, num_iters: i32) -> i32 {
    let mut step_mem =
        match lsrkStep_AccessStepMem(ark_mem, "LSRKStepSetNumDomEigEstPreprocessIters") {
            None => return ARK_MEM_NULL,
            Some(sm) => sm,
        };

    if num_iters < 0 {
        step_mem.num_warmups = DOM_EIG_NUM_WARMUPS_DEFAULT;
    } else {
        step_mem.num_warmups = num_iters;
    }

    /* Set the number of iterations immediately (if possible); see the
       C comment about the detach/reattach corner case. */
    if let Some(dee) = step_mem.DEE.as_mut() {
        let retval = SUNDomEigEstimator_SetNumPreprocessIters(dee, step_mem.num_init_warmups);
        if retval != SUN_SUCCESS {
            ark_mem.step_mem = Some(step_mem);
            arkProcessError(
                Some(ark_mem),
                ARK_DEE_FAIL,
                line!(),
                "LSRKStepSetNumDomEigEstPreprocessIters",
                file!(),
                "SUNDomeEigEstimator_SetNumPreprocessIters failed",
            );
            return ARK_DEE_FAIL;
        }
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  LSRKStepSetNumSSPStages sets the number of stages in the SSP
  methods.  Calling this function with num_of_stages <= 0 resets
  the default value.
  ---------------------------------------------------------------*/
pub fn LSRKStepSetNumSSPStages(ark_mem: &mut ARKodeMem, num_of_stages: i32) -> i32 {
    let mut step_mem = match lsrkStep_AccessStepMem(ark_mem, "LSRKStepSetNumSSPStages") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    if !step_mem.is_SSP {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "LSRKStepSetNumSSPStages",
            file!(),
            "Call this function only for SSP methods: Use LSRKStepSetSSPMethod to declare SSP method type first!",
        );
        return ARK_ILL_INPUT;
    }

    if num_of_stages <= 0 {
        match step_mem.LSRKmethod {
            ARKODE_LSRK_SSP_S_2 => step_mem.req_stages = 2,
            ARKODE_LSRK_SSP_S_3 => step_mem.req_stages = 4,
            ARKODE_LSRK_SSP_10_4 => step_mem.req_stages = 10,
            _ => {
                ark_mem.step_mem = Some(step_mem);
                arkProcessError(
                    Some(ark_mem),
                    ARK_ILL_INPUT,
                    line!(),
                    "LSRKStepSetNumSSPStages",
                    file!(),
                    "Call LSRKStepSetSSPMethod to declare SSP method type first!",
                );
                return ARK_ILL_INPUT;
            }
        }
        ark_mem.step_mem = Some(step_mem);
        return ARK_SUCCESS;
    } else {
        match step_mem.LSRKmethod {
            ARKODE_LSRK_SSP_S_2 => {
                if num_of_stages < 2 {
                    ark_mem.step_mem = Some(step_mem);
                    arkProcessError(
                        Some(ark_mem),
                        ARK_ILL_INPUT,
                        line!(),
                        "LSRKStepSetNumSSPStages",
                        file!(),
                        "num_of_stages must be greater than or equal to 2, or set it less than or equal to 0 to reset the default value",
                    );
                    return ARK_ILL_INPUT;
                }
            }
            ARKODE_LSRK_SSP_S_3 => {
                /* We check that num_of_stages is a perfect square (see the
                   C precision note on sqrt vs SUNRsqrt). */
                let root = (num_of_stages as f64).sqrt() as i32;
                if num_of_stages < 4 || root * root != num_of_stages {
                    ark_mem.step_mem = Some(step_mem);
                    arkProcessError(
                        Some(ark_mem),
                        ARK_ILL_INPUT,
                        line!(),
                        "LSRKStepSetNumSSPStages",
                        file!(),
                        "num_of_stages must be a perfect square greater than or equal to 4, or set it less than or equal to 0 to reset the default value",
                    );
                    return ARK_ILL_INPUT;
                }
                if num_of_stages == 4 {
                    ark_mem.step = Some(lsrkStep_TakeStepSSP43);
                } else {
                    ark_mem.step = Some(lsrkStep_TakeStepSSPs3);
                }
            }
            ARKODE_LSRK_SSP_10_4 => {
                if num_of_stages != 10 {
                    ark_mem.step_mem = Some(step_mem);
                    arkProcessError(
                        Some(ark_mem),
                        ARK_ILL_INPUT,
                        line!(),
                        "LSRKStepSetNumSSPStages",
                        file!(),
                        "SSP10_4 method has a prefixed num_of_stages = 10",
                    );
                    return ARK_ILL_INPUT;
                }
            }
            _ => {
                ark_mem.step_mem = Some(step_mem);
                arkProcessError(
                    Some(ark_mem),
                    ARK_ILL_INPUT,
                    line!(),
                    "LSRKStepSetNumSSPStages",
                    file!(),
                    "Call LSRKStepSetSSPMethod to declare SSP method type first!",
                );
                return ARK_ILL_INPUT;
            }
        }
        step_mem.req_stages = num_of_stages;
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  LSRKStepSetDomEigEstimator:

  This routine sets the dominant eigenvalue estimator DEE.  (The
  C SetATimes(DEE, arkode_mem, lsrkStep_DQJtimes) installation is
  replaced by the Estimate-time ATimes closure — see
  lsrkStep_ComputeNewDomEig.)
  ---------------------------------------------------------------*/
pub fn LSRKStepSetDomEigEstimator(
    ark_mem: &mut ARKodeMem,
    DEE: Option<SUNDomEigEstimator>,
) -> i32 {
    let mut step_mem = match lsrkStep_AccessStepMem(ark_mem, "LSRKStepSetDomEigEstimator") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* Attach the DEE to the step memory */
    step_mem.DEE = DEE;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*===============================================================
  Exported optional output functions.
  ===============================================================*/

/*---------------------------------------------------------------
  LSRKStepGetNumDomEigUpdates:

  Returns the number of dominant eigenvalue updates
  ---------------------------------------------------------------*/
pub fn LSRKStepGetNumDomEigUpdates(ark_mem: &mut ARKodeMem, dom_eig_num_evals: &mut i64) -> i32 {
    let step_mem = match lsrkStep_AccessStepMem(ark_mem, "LSRKStepGetNumDomEigUpdates") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* get values from step_mem */
    *dom_eig_num_evals = step_mem.dom_eig_num_evals;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  LSRKStepGetMaxNumStages:

  Returns the max number of stages used
  ---------------------------------------------------------------*/
pub fn LSRKStepGetMaxNumStages(ark_mem: &mut ARKodeMem, stage_max: &mut i32) -> i32 {
    let step_mem = match lsrkStep_AccessStepMem(ark_mem, "LSRKStepGetMaxNumStages") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* get values from step_mem */
    *stage_max = step_mem.stage_max;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  LSRKStepGetNumDomEigEstRhsEvals:

  Returns the number of RHS evals in DQ Jacobian computations
  ---------------------------------------------------------------*/
pub fn LSRKStepGetNumDomEigEstRhsEvals(ark_mem: &mut ARKodeMem, nfeDQ: &mut i64) -> i32 {
    let step_mem = match lsrkStep_AccessStepMem(ark_mem, "LSRKStepGetNumDomEigEstRhsEvals") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* get values from step_mem */
    *nfeDQ = step_mem.nfeDQ;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

pub fn LSRKStepGetNumDomEigEstIters(ark_mem: &mut ARKodeMem, num_iters: &mut i64) -> i32 {
    let step_mem = match lsrkStep_AccessStepMem(ark_mem, "LSRKStepGetNumDomEigEstIters") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* get values from step_mem */
    *num_iters = step_mem.num_dee_iters;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*===============================================================
  Private functions attached to ARKODE
  ===============================================================*/

/*---------------------------------------------------------------
  lsrkStep_SetDefaults:

  Resets all LSRKStep optional inputs to their default values.
  Does not change problem-defining function pointers or
  user_data pointer.
  ---------------------------------------------------------------*/
pub fn lsrkStep_SetDefaults(ark_mem: &mut ARKodeMem) -> i32 {
    let mut step_mem = match lsrkStep_AccessStepMem(ark_mem, "lsrkStep_SetDefaults") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* Set default values for integrator optional inputs
       (overwrite some adaptivity params for LSRKStep use) */
    step_mem.req_stages = 0; /* no stages */

    /* Spectral info */
    step_mem.dom_eig_safety = DOM_EIG_SAFETY_DEFAULT;
    step_mem.dom_eig_freq = DOM_EIG_FREQ_DEFAULT;
    step_mem.const_Jac = false;
    step_mem.num_init_warmups = DOM_EIG_NUM_INIT_WARMUPS_DEFAULT;
    step_mem.num_warmups = DOM_EIG_NUM_WARMUPS_DEFAULT;

    ark_mem.step_mem = Some(step_mem);

    /* Load the default SUNAdaptController */
    let retval = crate::arkode_io::arkReplaceAdaptController(ark_mem, None, true);
    if retval != 0 {
        return retval;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  lsrkStep_GetStageIndex:

  Returns the current stage index and number of stages
  ---------------------------------------------------------------*/
pub fn lsrkStep_GetStageIndex(ark_mem: &mut ARKodeMem, stage: &mut i32, max_stages: &mut i32) -> i32 {
    let step_mem = match lsrkStep_AccessStepMem(ark_mem, "lsrkStep_GetStageIndex") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    *stage = step_mem.istage;
    *max_stages = step_mem.req_stages + 1;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  lsrkStep_PrintAllStats:

  Prints integrator statistics for STS methods
  ---------------------------------------------------------------*/
pub fn lsrkStep_PrintAllStats(
    ark_mem: &mut ARKodeMem,
    outfile: &mut dyn std::io::Write,
    fmt: SUNOutputFormat,
) -> i32 {
    let step_mem = match lsrkStep_AccessStepMem(ark_mem, "lsrkStep_PrintAllStats") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    sunfprintf_long(outfile, fmt, SUNFALSE, "RHS fn evals", step_mem.nfe);
    if step_mem.is_SSP {
        sunfprintf_long(
            outfile,
            fmt,
            SUNFALSE,
            "Number of stages used",
            step_mem.req_stages as i64,
        );
    } else {
        sunfprintf_long(
            outfile,
            fmt,
            SUNFALSE,
            "Number of dom_eig updates",
            step_mem.dom_eig_num_evals,
        );
        if step_mem.DEE.is_some() {
            sunfprintf_long(
                outfile,
                fmt,
                SUNFALSE,
                "Number of fe calls for DEE",
                step_mem.nfeDQ,
            );
            sunfprintf_long(
                outfile,
                fmt,
                SUNFALSE,
                "Number of iterations for DEE",
                step_mem.num_dee_iters,
            );
        }
        sunfprintf_long(
            outfile,
            fmt,
            SUNFALSE,
            "Max. num. of stages used",
            step_mem.stage_max as i64,
        );
        sunfprintf_long(
            outfile,
            fmt,
            SUNFALSE,
            "Max. num. of stages allowed",
            step_mem.stage_max_limit as i64,
        );
        sunfprintf_real(
            outfile,
            fmt,
            SUNFALSE,
            "Max. spectral radius",
            step_mem.spectral_radius_max,
        );
        sunfprintf_real(
            outfile,
            fmt,
            SUNFALSE,
            "Min. spectral radius",
            step_mem.spectral_radius_min,
        );
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  lsrkStep_WriteParameters:

  Outputs all solver parameters to the provided file pointer.
  ---------------------------------------------------------------*/
pub fn lsrkStep_WriteParameters(ark_mem: &mut ARKodeMem, fp: &mut dyn std::io::Write) -> i32 {
    use crate::sundials_utils::fmt_g;

    let step_mem = match lsrkStep_AccessStepMem(ark_mem, "lsrkStep_WriteParameters") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* print integrator parameters to file */
    match step_mem.LSRKmethod {
        ARKODE_LSRK_RKC_2 => {
            let _ = write!(fp, "LSRKStep RKC time step module parameters:\n");
        }
        ARKODE_LSRK_RKL_2 => {
            let _ = write!(fp, "LSRKStep RKL time step module parameters:\n");
        }
        ARKODE_LSRK_SSP_S_2 => {
            let _ = write!(fp, "LSRKStep SSP(s,2) time step module parameters:\n");
        }
        ARKODE_LSRK_SSP_S_3 => {
            let _ = write!(fp, "LSRKStep SSP(s,3) time step module parameters:\n");
        }
        ARKODE_LSRK_SSP_10_4 => {
            let _ = write!(fp, "LSRKStep SSP(10,4) time step module parameters:\n");
        }
        _ => {
            ark_mem.step_mem = Some(step_mem);
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "lsrkStep_WriteParameters",
                file!(),
                "Invalid method option.",
            );
            return ARK_ILL_INPUT;
        }
    }

    let _ = write!(fp, "  Method order {}\n", step_mem.q);
    let _ = write!(fp, "  Embedding order {}\n", step_mem.p);

    if step_mem.is_SSP {
        let _ = write!(fp, "  Number of stages used = {}\n", step_mem.req_stages);
    } else {
        let _ = write!(
            fp,
            "  Maximum number of stages allowed = {}\n",
            step_mem.stage_max_limit
        );
        if step_mem.DEE.is_some() {
            let _ = write!(fp, "  Number of fe calls for DEE = {}\n", step_mem.nfeDQ);
        }
        let _ = write!(
            fp,
            "  Current spectral radius = {}\n",
            fmt_g(step_mem.spectral_radius, 0, 15)
        );
        let _ = write!(
            fp,
            "  Safety factor for the dom eig = {}\n",
            fmt_g(step_mem.dom_eig_safety, 0, 15)
        );
        let _ = write!(
            fp,
            "  Max num of successful steps before new dom eig update = {}\n",
            step_mem.dom_eig_freq
        );
        let _ = write!(
            fp,
            "  Number of first preprocessing warmups = {}\n",
            step_mem.num_init_warmups
        );
        let _ = write!(
            fp,
            "  Number of subsequent preprocessing warmups = {}\n",
            step_mem.num_warmups
        );
        let _ = write!(
            fp,
            "  Flag to indicate Jacobian is constant = {}\n",
            step_mem.const_Jac as i32
        );
    }

    let _ = write!(fp, "\n");

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  lsrkStep_GetNumRhsEvals:

  Returns the current number of RHS calls
  ---------------------------------------------------------------*/
pub fn lsrkStep_GetNumRhsEvals(
    ark_mem: &mut ARKodeMem,
    partition_index: i32,
    rhs_evals: &mut i64,
) -> i32 {
    let step_mem = match lsrkStep_AccessStepMem(ark_mem, "lsrkStep_GetNumRhsEvals") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    if partition_index > 0 {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "lsrkStep_GetNumRhsEvals",
            file!(),
            "Invalid partition index",
        );
        return ARK_ILL_INPUT;
    }

    *rhs_evals = step_mem.nfe;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  lsrkStep_GetEstLocalErrors: Returns the current local truncation
  error estimate vector
  ---------------------------------------------------------------*/
pub fn lsrkStep_GetEstLocalErrors(ark_mem: &mut ARKodeMem, ele: &mut NVector) -> i32 {
    let step_mem = match lsrkStep_AccessStepMem(ark_mem, "lsrkStep_GetEstLocalErrors") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };
    ark_mem.step_mem = Some(step_mem);

    /* return an error if local truncation error is not computed */
    if ark_mem.fixedstep {
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* otherwise, copy local truncation error vector to output */
    N_VScale(ONE, &ark_mem.tempv1, ele);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  lsrkStep_SetOptions:

  Provides command-line control over LSRKStep-specific "set"
  routines (arkode_lsrkstep_io.c).
  ---------------------------------------------------------------*/
pub fn lsrkStep_SetOptions(
    ark_mem: &mut ARKodeMem,
    argidx: &mut usize,
    argv: &[String],
    offset: usize,
    arg_used: &mut bool,
) -> i32 {
    use crate::sundials_cli::{
        sunCheckAndSetCharArgs, sunCheckAndSetIntArgs, sunCheckAndSetLongArgs,
        sunCheckAndSetRealArgs, sunKeyCharPair, sunKeyIntPair, sunKeyLongPair, sunKeyRealPair,
    };

    /* Set lists of keys, and the corresponding set routines */
    let char_pairs: [sunKeyCharPair<ARKodeMem>; 2] = [
        sunKeyCharPair { key: "sts_method_name", set: LSRKStepSetSTSMethodByName },
        sunKeyCharPair { key: "ssp_method_name", set: LSRKStepSetSSPMethodByName },
    ];

    let long_pairs: [sunKeyLongPair<ARKodeMem>; 1] = [sunKeyLongPair {
        key: "dom_eig_frequency",
        set: LSRKStepSetDomEigFrequency,
    }];

    let int_pairs: [sunKeyIntPair<ARKodeMem>; 4] = [
        sunKeyIntPair { key: "max_num_stages", set: LSRKStepSetMaxNumStages },
        sunKeyIntPair { key: "num_ssp_stages", set: LSRKStepSetNumSSPStages },
        sunKeyIntPair {
            key: "num_dom_eig_est_init_preprocess_iters",
            set: LSRKStepSetNumDomEigEstInitPreprocessIters,
        },
        sunKeyIntPair {
            key: "num_dom_eig_est_preprocess_iters",
            set: LSRKStepSetNumDomEigEstPreprocessIters,
        },
    ];

    let real_pairs: [sunKeyRealPair<ARKodeMem>; 1] = [sunKeyRealPair {
        key: "dom_eig_safety_factor",
        set: LSRKStepSetDomEigSafetyFactor,
    }];

    /* check all "char" keys */
    let mut j: usize = 0;
    let retval =
        sunCheckAndSetCharArgs(ark_mem, argidx, argv, offset, &char_pairs, arg_used, &mut j);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            retval,
            line!(),
            "lsrkStep_SetOptions",
            file!(),
            &format!("error setting key: {}", char_pairs[j].key),
        );
        return retval;
    }
    if *arg_used {
        return ARK_SUCCESS;
    }

    /* check all "long int" keys */
    let retval =
        sunCheckAndSetLongArgs(ark_mem, argidx, argv, offset, &long_pairs, arg_used, &mut j);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            retval,
            line!(),
            "lsrkStep_SetOptions",
            file!(),
            &format!("error setting key: {}", long_pairs[j].key),
        );
        return retval;
    }
    if *arg_used {
        return ARK_SUCCESS;
    }

    /* check all "int" keys */
    let retval =
        sunCheckAndSetIntArgs(ark_mem, argidx, argv, offset, &int_pairs, arg_used, &mut j);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            retval,
            line!(),
            "lsrkStep_SetOptions",
            file!(),
            &format!("error setting key: {}", int_pairs[j].key),
        );
        return retval;
    }
    if *arg_used {
        return ARK_SUCCESS;
    }

    /* check all "real" keys */
    let retval =
        sunCheckAndSetRealArgs(ark_mem, argidx, argv, offset, &real_pairs, arg_used, &mut j);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            retval,
            line!(),
            "lsrkStep_SetOptions",
            file!(),
            &format!("error setting key: {}", real_pairs[j].key),
        );
        return retval;
    }

    ARK_SUCCESS
}
