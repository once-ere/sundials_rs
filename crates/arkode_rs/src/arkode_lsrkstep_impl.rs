/* -----------------------------------------------------------------
 * Translation of src/arkode/arkode_lsrkstep_impl.h (+ the
 * ARKODE_LSRKMethodType enum and ARKDomEigFn typedef of
 * include/arkode/arkode_lsrkstep.h).  LSRKStep time-stepper module
 * memory structure (super-time-stepping RKC/RKL and SSP methods).
 *
 * The SUNDomEigEstimator is stored by value (None = C NULL); its
 * ATimes callback is supplied at Estimate time (see the pinned
 * adaptation note in sundials_domeigestimator.rs).
 * -----------------------------------------------------------------*/

use crate::arkode_impl::ARKRhsFn;
use crate::nvector_serial::NVector;
use crate::sundials_domeigestimator::SUNDomEigEstimator;
use crate::sundials_types::UserData;

/* LSRK time step module constants */
pub const STAGE_MAX_LIMIT_DEFAULT: i32 = 200;
pub const DOM_EIG_SAFETY_DEFAULT: f64 = 1.01;
pub const DOM_EIG_FREQ_DEFAULT: i64 = 25;
pub const DOM_EIG_NUM_WARMUPS_DEFAULT: i32 = 0;
pub const DOM_EIG_NUM_INIT_WARMUPS_DEFAULT: i32 = -1; /* use DEE's default value */

/* ARKODE_LSRKMethodType (arkode_lsrkstep.h) */
pub type ARKODE_LSRKMethodType = i32;
pub const ARKODE_LSRK_RKC_2: ARKODE_LSRKMethodType = 0;
pub const ARKODE_LSRK_RKL_2: ARKODE_LSRKMethodType = 1;
pub const ARKODE_LSRK_SSP_S_2: ARKODE_LSRKMethodType = 2;
pub const ARKODE_LSRK_SSP_S_3: ARKODE_LSRKMethodType = 3;
pub const ARKODE_LSRK_SSP_10_4: ARKODE_LSRKMethodType = 4;

/* ARKDomEigFn (arkode_lsrkstep.h) */
pub type ARKDomEigFn = fn(
    t: f64,
    y: &NVector,
    fn_: &NVector,
    lambdaR: &mut f64,
    lambdaI: &mut f64,
    user_data: &mut UserData,
    temp1: &mut NVector,
    temp2: &mut NVector,
    temp3: &mut NVector,
) -> i32;

/// struct ARKodeLSRKStepMemRec (arkode_lsrkstep_impl.h)
#[derive(Default)]
pub struct ARKodeLSRKStepMem {
    /* LSRK problem specification */
    pub fe: Option<ARKRhsFn>,
    pub dom_eig_fn: Option<ARKDomEigFn>,

    pub q: i32, /* method order    */
    pub p: i32, /* embedding order */

    pub istage: i32,     /* current stage            */
    pub req_stages: i32, /* number of stages in step */

    pub LSRKmethod: ARKODE_LSRKMethodType,

    /* Counters and stats */
    pub nfe: i64,   /* num fe calls */
    pub nfeDQ: i64, /* num fe calls for difference quotient approximation */
    pub dom_eig_num_evals: i64, /* num of dom_eig computations */
    pub stage_max: i32,         /* num of max stages used      */
    pub stage_max_limit: i32,   /* max allowed num of stages   */
    pub dom_eig_nst: i64, /* num of step at which the last dominant eigenvalue was computed */
    pub step_nst: i64,    /* The number of successful steps. */
    pub num_dee_iters: i64, /* number of iterations in the DEE estimates */

    /* Spectral info */
    pub lambdaR: f64,             /* Real part of the dominated eigenvalue */
    pub lambdaI: f64,             /* Imaginary part of the dominated eigenvalue */
    pub spectral_radius: f64,     /* spectral radius */
    pub spectral_radius_max: f64, /* max spectral radius */
    pub spectral_radius_min: f64, /* min spectral radius */
    pub dom_eig_safety: f64, /* some safety factor for the user provided dom_eig */
    pub dom_eig_freq: i64, /* indicates dom_eig update after dom_eig_freq successful steps */
    pub num_init_warmups: i32, /* number of warm-ups in the first DEE estimates */
    pub num_warmups: i32,      /* number of warm-ups in succeeding DEE estimates */

    pub DEE: Option<SUNDomEigEstimator>, /* DomEig estimator */

    /* Flags */
    pub dom_eig_update: bool, /* flag indicating new dom_eig is needed */
    pub const_Jac: bool,      /* flag indicating Jacobian is constant */
    pub dom_eig_is_current: bool, /* SUNTRUE if dom_eig has been evaluated at tn */
    pub is_SSP: bool,             /* flag indicating SSP method */
    pub init_warmup: bool,        /* flag indicating initial warm-up */

    /* Reusable fused vector operation arrays (Xvecs: operand lists are
       assembled at the call sites; the liw accounting is kept) */
    pub cvals: Vec<f64>,
    pub nfusedopvecs: i32, /* length of cvals and Xvecs arrays */
}

/* Initialization and I/O error messages */
pub const MSG_LSRKSTEP_NO_MEM: &str = "Time step module memory is NULL.";
