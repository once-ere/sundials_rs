/* -----------------------------------------------------------------
 * Translation of src/arkode/arkode_mristep_impl.h (+ the
 * MRIStepInnerStepper types from include/arkode/arkode_mristep.h),
 * SUNDIALS 7.7.0.
 *
 * Implementation header for ARKODE's MRI time stepper module.
 *
 * Storage adaptations (same conventions as the other steppers):
 *  - `N_Vector* Fse/Fsi` become owning Vecs.  C's `unify_Fs` mode
 *    aliases Fsi = Fse (one allocation); here the storage lives in
 *    Fse and the Fsi-indexed accesses at the MRISR call sites go to
 *    Fse when unify_Fs is set.
 *  - Xvecs operand lists are assembled at call sites; cvals is kept
 *    (with its lrw/liw accounting) like the other steppers.
 *  - the ARKLS interface memory is hoisted onto ARKodeMem.lmem
 *    (ARCHITECTURE.md Addendum C.1); step_mem keeps only the
 *    linit/lsetup/lsolve/lfree fn pointers.
 *  - `MRIStepInnerStepper` is an owned struct; its `void* content`
 *    is the usual Option<Box<dyn Any>> and owns the wrapped inner
 *    integrator (Box<ARKodeMem>) or SUNStepper.  The C `sunctx`,
 *    `python` fields and the fused-op vals/vecs workspace (assembled
 *    at call sites) are dropped.
 * -----------------------------------------------------------------*/

use crate::arkode_impl::{
    ARKLinsolFreeFn, ARKLinsolInitFn, ARKLinsolSetupFn, ARKLinsolSolveFn, ARKRhsFn,
    ARKStagePredictFn, ARK_NO_FAILURES,
};
use crate::arkode_mri_tables::MRIStepCoupling;
use crate::nvector_serial::NVector;
use crate::sundials_nonlinearsolver::NonlinearSolver;
use crate::sundials_types::UserData;

/* Stage type identifiers */
pub const MRISTAGE_FIRST: i32 = -2;
pub const MRISTAGE_STIFF_ACC: i32 = -1;
pub const MRISTAGE_ERK_FAST: i32 = 0;
pub const MRISTAGE_ERK_NOFAST: i32 = 1;
pub const MRISTAGE_DIRK_NOFAST: i32 = 2;
pub const MRISTAGE_DIRK_FAST: i32 = 3;

/* Implicit solver constants (duplicate from arkode_arkstep_impl.h) */
/*   max number of nonlinear iterations */
pub const MAXCOR: i32 = 3;
/*   constant to estimate the convergence rate for the nonlinear equation */
pub const CRDOWN: f64 = 0.3;
/*   if |gamma/gammap-1| > DGMAX then call lsetup */
pub const DGMAX: f64 = 0.2;
/*   declare divergence if ratio del/delp > RDIV */
pub const RDIV: f64 = 2.3;
/*   max no. of steps between lsetup calls */
pub const MSBP: i32 = 20;
/*   default solver tolerance factor */
pub const NLSCOEF: f64 = 0.1;

/* ------------------------------------
 * MRIStep Inner Stepper Function Types
 * ------------------------------------ */

pub type MRIStepInnerEvolveFn =
    fn(stepper: &mut MRIStepInnerStepper, t0: f64, tout: f64, y: &mut NVector) -> i32;

pub type MRIStepInnerFullRhsFn =
    fn(stepper: &mut MRIStepInnerStepper, t: f64, y: &NVector, f: &mut NVector, mode: i32) -> i32;

pub type MRIStepInnerResetFn = fn(stepper: &mut MRIStepInnerStepper, tR: f64, yR: &NVector) -> i32;

pub type MRIStepInnerGetAccumulatedError =
    fn(stepper: &mut MRIStepInnerStepper, accum_error: &mut f64) -> i32;

pub type MRIStepInnerResetAccumulatedError = fn(stepper: &mut MRIStepInnerStepper) -> i32;

pub type MRIStepInnerSetRTol = fn(stepper: &mut MRIStepInnerStepper, rtol: f64) -> i32;

/* ------------------------------
 * User-Supplied Function Types
 * ------------------------------ */

pub type MRIStepPreInnerFn =
    fn(t: f64, f: &[NVector], nvecs: i32, user_data: &mut UserData) -> i32;

pub type MRIStepPostInnerFn = fn(t: f64, y: &mut NVector, user_data: &mut UserData) -> i32;

/*===============================================================
  MRI inner time stepper data structure
  ===============================================================*/

/// struct _MRIStepInnerStepper_Ops
#[derive(Default)]
pub struct MRIStepInnerStepper_Ops {
    pub evolve: Option<MRIStepInnerEvolveFn>,
    pub fullrhs: Option<MRIStepInnerFullRhsFn>,
    pub reset: Option<MRIStepInnerResetFn>,
    pub geterror: Option<MRIStepInnerGetAccumulatedError>,
    pub reseterror: Option<MRIStepInnerResetAccumulatedError>,
    pub setrtol: Option<MRIStepInnerSetRTol>,
}

/// struct _MRIStepInnerStepper
#[derive(Default)]
pub struct MRIStepInnerStepper {
    /* stepper specific content and operations */
    pub content: UserData,
    pub ops: MRIStepInnerStepper_Ops,

    /* base class data */
    pub forcing: Vec<NVector>,   /* array of forcing vectors            */
    pub nforcing: i32,           /* number of forcing vectors active    */
    pub nforcing_allocated: i32, /* number of forcing vectors allocated */
    pub last_flag: i32,          /* last stepper return flag            */
    pub tshift: f64,             /* time normalization shift            */
    pub tscale: f64,             /* time normalization scaling          */

    /* Space requirements */
    pub lrw1: i64, /* no. of sunrealtype words in 1 N_Vector          */
    pub liw1: i64, /* no. of integer words in 1 N_Vector           */
    pub lrw: i64,  /* no. of sunrealtype words in ARKODE work vectors */
    pub liw: i64,  /* no. of integer words in ARKODE work vectors  */
}

/*===============================================================
  MRI time step module data structure
  ===============================================================*/

/// struct ARKodeMRIStepMemRec
pub struct ARKodeMRIStepMem {
    /* MRI problem specification */
    pub fse: Option<ARKRhsFn>, /* y' = fse(t,y) + fsi(t,y) + ff(t,y) */
    pub fsi: Option<ARKRhsFn>,
    pub linear: bool,         /* SUNTRUE if fi is linear        */
    pub linear_timedep: bool, /* SUNTRUE if dfi/dy depends on t */
    pub explicit_rhs: bool,   /* SUNTRUE if fse is provided     */
    pub implicit_rhs: bool,   /* SUNTRUE if fsi is provided     */
    pub deduce_rhs: bool,     /* SUNTRUE if fi is deduced after
                              a nonlinear solve              */

    /* Outer RK method storage and parameters */
    pub Fse: Vec<NVector>, /* explicit RHS at each stage (also holds
                           the unified storage when unify_Fs)      */
    pub Fsi: Vec<NVector>, /* implicit RHS at each stage               */
    pub unify_Fs: bool,    /* Fse and Fsi point at the same memory     */
    pub fse_is_current: bool,
    pub fsi_is_current: bool,
    pub MRIC: Option<MRIStepCoupling>, /* slow->fast coupling table   */
    pub q: i32,                        /* method order                             */
    pub p: i32,                        /* embedding order                          */
    pub stages: i32,                   /* total number of stages                   */
    pub nstages_active: i32,           /* number of active stage RHS vectors       */
    pub nstages_allocated: i32,        /* number of stage RHS vectors allocated    */
    pub stage_map: Vec<i32>,           /* index map for storing stage RHS vectors  */
    pub stagetypes: Vec<i32>,          /* type flags for stages                    */
    pub Ae_row: Vec<f64>,              /* equivalent explicit RK coeffs            */
    pub Ai_row: Vec<f64>,              /* equivalent implicit RK coeffs            */

    /* Algebraic solver data and parameters */
    pub sdata: NVector,                 /* old stage data in residual               */
    pub zpred: NVector,                 /* predicted stage solution                 */
    pub zcor: NVector,                  /* stage correction                         */
    pub NLS: Option<NonlinearSolver>,   /* generic SUNNonlinearSolver object        */
    pub ownNLS: bool,                   /* flag indicating ownership of NLS         */
    pub nls_fsi: Option<ARKRhsFn>,      /* fsi(t,y) used in the nonlinear solver    */
    pub gamma: f64,                     /* gamma = h * A(i,i)                       */
    pub gammap: f64,                    /* gamma at the last setup call             */
    pub gamrat: f64,                    /* gamma / gammap                           */
    pub dgmax: f64,                     /* call lsetup if |gamma/gammap-1| >= dgmax */
    pub predictor: i32,                 /* implicit prediction method to use        */
    pub crdown: f64,                    /* nonlinear conv rate estimation constant  */
    pub rdiv: f64,                      /* nonlin divergence if del/delp > rdiv     */
    pub crate_: f64,                    /* estimated nonlin convergence rate        */
    pub delp: f64,                      /* norm of previous nonlinear solver update */
    pub eRNrm: f64,                     /* estimated residual norm                  */
    pub nlscoef: f64,                   /* coefficient in nonlin. convergence test  */
    pub msbp: i32,                      /* positive => max # steps between lsetup
                                        negative => call at each Newton iter     */
    pub nstlp: i64,                     /* step number of last setup call           */
    pub maxcor: i32,                    /* max num iterations for solving the
                                        nonlinear equation                       */
    pub convfail: i32,                  /* NLS fail flag (for interface routines)   */
    pub jcur: bool,                     /* is Jacobian info for lin solver current? */
    pub stage_predict: Option<ARKStagePredictFn>, /* User-supplied stage predictor  */
    pub istage: i32,                    /* stage index used in nonlinear solve      */

    /* Informational output for mriStep_GetStageIndex */
    pub cur_stage: i32,

    /* Linear Solver Data (lmem hoisted onto ARKodeMem.lmem) */
    pub linit: Option<ARKLinsolInitFn>,
    pub lsetup: Option<ARKLinsolSetupFn>,
    pub lsolve: Option<ARKLinsolSolveFn>,
    pub lfree: Option<ARKLinsolFreeFn>,

    /* Inner stepper */
    pub stepper: MRIStepInnerStepper,

    /* User-supplied pre and post inner evolve functions */
    pub pre_inner_evolve: Option<MRIStepPreInnerFn>,
    pub post_inner_evolve: Option<MRIStepPostInnerFn>,

    /* MRI adaptivity parameters */
    pub inner_rtol_factor: f64,     /* prev control parameter */
    pub inner_dsm: f64,             /* prev inner stepper accumulated error */
    pub inner_rtol_factor_new: f64, /* upcoming control parameter */

    /* Counters */
    pub nfse: i64,         /* num fse calls                    */
    pub nfsi: i64,         /* num fsi calls                    */
    pub nsetups: i64,      /* num linear solver setup calls    */
    pub nls_iters: i64,    /* num nonlinear solver iters       */
    pub nls_fails: i64,    /* num nonlinear solver fails       */
    pub inner_fails: i64,  /* num recov. inner solver fails  */
    pub nfusedopvecs: i32, /* length of cvals array            */

    /* Data for using MRIStep with external polynomial forcing */
    pub expforcing: bool,      /* add forcing to explicit RHS */
    pub impforcing: bool,      /* add forcing to implicit RHS */
    pub tshift: f64,           /* time normalization shift    */
    pub tscale: f64,           /* time normalization scaling  */
    pub forcing: Vec<NVector>, /* array of forcing vectors    */
    pub nforcing: i32,         /* number of forcing vectors   */

    /* Reusable array for fused vector operations (Xvecs assembled
    at call sites) */
    pub cvals: Vec<f64>,
}

impl Default for ARKodeMRIStepMem {
    /* C calloc(1, ...) zero-initialization */
    fn default() -> Self {
        ARKodeMRIStepMem {
            fse: None,
            fsi: None,
            linear: false,
            linear_timedep: false,
            explicit_rhs: false,
            implicit_rhs: false,
            deduce_rhs: false,
            Fse: Vec::new(),
            Fsi: Vec::new(),
            unify_Fs: false,
            fse_is_current: false,
            fsi_is_current: false,
            MRIC: None,
            q: 0,
            p: 0,
            stages: 0,
            nstages_active: 0,
            nstages_allocated: 0,
            stage_map: Vec::new(),
            stagetypes: Vec::new(),
            Ae_row: Vec::new(),
            Ai_row: Vec::new(),
            sdata: NVector::new(0),
            zpred: NVector::new(0),
            zcor: NVector::new(0),
            NLS: None,
            ownNLS: false,
            nls_fsi: None,
            gamma: 0.0,
            gammap: 0.0,
            gamrat: 0.0,
            dgmax: 0.0,
            predictor: 0,
            crdown: 0.0,
            rdiv: 0.0,
            crate_: 0.0,
            delp: 0.0,
            eRNrm: 0.0,
            nlscoef: 0.0,
            msbp: 0,
            nstlp: 0,
            maxcor: 0,
            convfail: ARK_NO_FAILURES,
            jcur: false,
            stage_predict: None,
            istage: 0,
            cur_stage: 0,
            linit: None,
            lsetup: None,
            lsolve: None,
            lfree: None,
            stepper: MRIStepInnerStepper::default(),
            pre_inner_evolve: None,
            post_inner_evolve: None,
            inner_rtol_factor: 0.0,
            inner_dsm: 0.0,
            inner_rtol_factor_new: 0.0,
            nfse: 0,
            nfsi: 0,
            nsetups: 0,
            nls_iters: 0,
            nls_fails: 0,
            inner_fails: 0,
            nfusedopvecs: 0,
            expforcing: false,
            impforcing: false,
            tshift: 0.0,
            tscale: 0.0,
            forcing: Vec::new(),
            nforcing: 0,
            cvals: Vec::new(),
        }
    }
}

/*===============================================================
  Reusable MRIStep Error Messages
  ===============================================================*/

/* Initialization and I/O error messages */
pub const MSG_MRISTEP_NO_MEM: &str = "Time step module memory is NULL.";
pub const MSG_NLS_INIT_FAIL: &str = "The nonlinear solver's init routine failed.";
pub const MSG_MRISTEP_NO_COUPLING: &str = "The MRIStepCoupling is NULL.";
