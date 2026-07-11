/* -----------------------------------------------------------------
 * Translation of sundials-7.7.0/src/arkode/arkode_arkstep_impl.h
 * (+ the ARKSTEP_DEFAULT_* constants of include/arkode/
 * arkode_arkstep.h).
 *
 * Modeling notes (following arkode_erkstep_impl.rs):
 *  - `N_Vector* Fe/Fi/z` -> Vec<NVector> (empty = unallocated).
 *  - `Xvecs` (scratch array of N_Vector POINTERS for fused ops)
 *    cannot be stored in safe Rust; operand lists are assembled at
 *    each call site.  The liw accounting for it is kept
 *    (nfusedopvecs).
 *  - `void* lmem` lives on ARKodeMem.lmem (Addendum C.1: hoisted so
 *    the ARKLS interface avoids double-nested take/put-back); this
 *    struct keeps the linit/lsetup/lsolve/lfree pointers and
 *    lsolve_type.
 *  - `fn_implicit` (C: alias to a saved implicit RHS evaluation for
 *    the deduce_rhs/trivial-predictor-autonomous residual forms) is
 *    re-derived at the use sites in arkode_arkstep_nls.rs instead of
 *    being stored as an alias.
 *  - Mass-matrix solver data (minit/msetup/mmult/msolve/mfree,
 *    mass_mem, msolve_type) is deferred with the ARKLS mass half;
 *    mass_type is kept since the stepper logic branches on it.
 *  - `adj_fe` (SUNAdjRhsFn) is deferred with the adjoint machinery
 *    (needs the ManyVector module).
 *  - `forcing` in C aliases the caller's vectors (MRIStep); here it
 *    is an owned copy refreshed on every arkStep_SetInnerForcing
 *    call.
 * -----------------------------------------------------------------*/

use crate::arkode_butcher::ARKodeButcherTable;
use crate::arkode_butcher_dirk::{
    ARKODE_ARK2_DIRK_3_1_2, ARKODE_ARK324L2SA_DIRK_4_2_3, ARKODE_ARK437L2SA_DIRK_7_3_4,
    ARKODE_ARK548L2SAb_DIRK_8_4_5, ARKODE_BACKWARD_EULER_1_1, ARKODE_ESDIRK325L2SA_5_2_3,
    ARKODE_ESDIRK436L2SA_6_3_4, ARKODE_ESDIRK547L2SA2_7_4_5, ARKODE_DIRKTableID,
};
use crate::arkode_butcher_erk::{
    ARKODE_ARK2_ERK_3_1_2, ARKODE_ARK324L2SA_ERK_4_2_3, ARKODE_ARK437L2SA_ERK_7_3_4,
    ARKODE_ARK548L2SAb_ERK_8_4_5, ARKODE_BOGACKI_SHAMPINE_4_2_3, ARKODE_FORWARD_EULER_1_1,
    ARKODE_RALSTON_3_1_2, ARKODE_SOFRONIOU_SPALETTA_5_3_4, ARKODE_TSITOURAS_7_4_5,
    ARKODE_VERNER_10_6_7, ARKODE_VERNER_13_7_8, ARKODE_VERNER_16_8_9, ARKODE_VERNER_9_5_6,
    ARKODE_ERKTableID,
};
use crate::arkode_impl::{
    ARKLinsolFreeFn, ARKLinsolInitFn, ARKLinsolSetupFn, ARKLinsolSolveFn, ARKRhsFn,
    ARKStagePredictFn,
};
use crate::nvector_serial::NVector;
use crate::sundials_linearsolver::SUNLinearSolver_Type;
use crate::sundials_nonlinearsolver::NonlinearSolver;

/*===============================================================
  ARK time step module constants
  ===============================================================*/

/* max number of nonlinear iterations */
pub const MAXCOR: i32 = 3;
/* constant to estimate the convergence rate for the nonlinear equation */
pub const CRDOWN: f64 = 0.3;
/* if |gamma/gammap-1| > DGMAX then call lsetup */
pub const DGMAX: f64 = 0.2;
/* declare divergence if ratio del/delp > RDIV */
pub const RDIV: f64 = 2.3;
/* max no. of steps between lsetup calls */
pub const MSBP: i32 = 20;

/* Default solver tolerance factor */
pub const NLSCOEF: f64 = 0.1;

/* Mass matrix types */
pub const MASS_IDENTITY: i32 = 0;
pub const MASS_FIXED: i32 = 1;
pub const MASS_TIMEDEP: i32 = 2;

/* Default Butcher tables per order (arkode_arkstep.h) */
pub const ARKSTEP_DEFAULT_ERK_1: ARKODE_ERKTableID = ARKODE_FORWARD_EULER_1_1;
pub const ARKSTEP_DEFAULT_ERK_2: ARKODE_ERKTableID = ARKODE_RALSTON_3_1_2;
pub const ARKSTEP_DEFAULT_ERK_3: ARKODE_ERKTableID = ARKODE_BOGACKI_SHAMPINE_4_2_3;
pub const ARKSTEP_DEFAULT_ERK_4: ARKODE_ERKTableID = ARKODE_SOFRONIOU_SPALETTA_5_3_4;
pub const ARKSTEP_DEFAULT_ERK_5: ARKODE_ERKTableID = ARKODE_TSITOURAS_7_4_5;
pub const ARKSTEP_DEFAULT_ERK_6: ARKODE_ERKTableID = ARKODE_VERNER_9_5_6;
pub const ARKSTEP_DEFAULT_ERK_7: ARKODE_ERKTableID = ARKODE_VERNER_10_6_7;
pub const ARKSTEP_DEFAULT_ERK_8: ARKODE_ERKTableID = ARKODE_VERNER_13_7_8;
pub const ARKSTEP_DEFAULT_ERK_9: ARKODE_ERKTableID = ARKODE_VERNER_16_8_9;

pub const ARKSTEP_DEFAULT_DIRK_1: ARKODE_DIRKTableID = ARKODE_BACKWARD_EULER_1_1;
pub const ARKSTEP_DEFAULT_DIRK_2: ARKODE_DIRKTableID = ARKODE_ARK2_DIRK_3_1_2;
pub const ARKSTEP_DEFAULT_DIRK_3: ARKODE_DIRKTableID = ARKODE_ESDIRK325L2SA_5_2_3;
pub const ARKSTEP_DEFAULT_DIRK_4: ARKODE_DIRKTableID = ARKODE_ESDIRK436L2SA_6_3_4;
pub const ARKSTEP_DEFAULT_DIRK_5: ARKODE_DIRKTableID = ARKODE_ESDIRK547L2SA2_7_4_5;

pub const ARKSTEP_DEFAULT_ARK_ETABLE_2: ARKODE_ERKTableID = ARKODE_ARK2_ERK_3_1_2;
pub const ARKSTEP_DEFAULT_ARK_ETABLE_3: ARKODE_ERKTableID = ARKODE_ARK324L2SA_ERK_4_2_3;
pub const ARKSTEP_DEFAULT_ARK_ETABLE_4: ARKODE_ERKTableID = ARKODE_ARK437L2SA_ERK_7_3_4;
pub const ARKSTEP_DEFAULT_ARK_ETABLE_5: ARKODE_ERKTableID = ARKODE_ARK548L2SAb_ERK_8_4_5;
pub const ARKSTEP_DEFAULT_ARK_ITABLE_2: ARKODE_DIRKTableID = ARKODE_ARK2_DIRK_3_1_2;
pub const ARKSTEP_DEFAULT_ARK_ITABLE_3: ARKODE_DIRKTableID = ARKODE_ARK324L2SA_DIRK_4_2_3;
pub const ARKSTEP_DEFAULT_ARK_ITABLE_4: ARKODE_DIRKTableID = ARKODE_ARK437L2SA_DIRK_7_3_4;
pub const ARKSTEP_DEFAULT_ARK_ITABLE_5: ARKODE_DIRKTableID = ARKODE_ARK548L2SAb_DIRK_8_4_5;

/*===============================================================
  Reusable ARKStep Error Messages
  ===============================================================*/

/* Initialization and I/O error messages */
pub const MSG_ARKSTEP_NO_MEM: &str = "Time step module memory is NULL.";
pub const MSG_NLS_INIT_FAIL: &str = "The nonlinear solver's init routine failed.";

/* Other error messages */
pub const MSG_ARK_MISSING_FE: &str =
    "Cannot specify that method is explicit without providing a function pointer to fe(t,y).";
pub const MSG_ARK_MISSING_FI: &str =
    "Cannot specify that method is implicit without providing a function pointer to fi(t,y).";
pub const MSG_ARK_MISSING_F: &str =
    "Cannot specify that method is ImEx without providing function pointers to fi(t,y) and fe(t,y).";

/// struct ARKodeARKStepMemRec (arkode_arkstep_impl.h)
pub struct ARKodeARKStepMem {
    /* ARK problem specification */
    pub fe: Option<ARKRhsFn>, /* My' = fe(t,y) + fi(t,y) */
    pub fi: Option<ARKRhsFn>,
    pub autonomous: bool,     /* SUNTRUE if fi depends on t     */
    pub linear: bool,         /* SUNTRUE if fi is linear        */
    pub linear_timedep: bool, /* SUNTRUE if dfi/dy depends on t */
    pub explicit: bool,       /* SUNTRUE if fe is enabled       */
    pub implicit: bool,       /* SUNTRUE if fi is enabled       */
    pub deduce_rhs: bool,     /* SUNTRUE if fi is deduced after
                              a nonlinear solve               */

    /* (adj_fe: SUNAdjRhsFn — deferred with the adjoint machinery) */

    /* ARK method storage and parameters */
    pub Fe: Vec<NVector>,  /* explicit RHS at each stage */
    pub Fi: Vec<NVector>,  /* implicit RHS at each stage */
    pub z: Vec<NVector>,   /* stages (for relaxation)    */
    pub sdata: NVector,    /* old stage data in residual */
    pub zpred: NVector,    /* predicted stage solution   */
    pub zcor: NVector,     /* stage correction           */
    pub q: i32,            /* method order               */
    pub p: i32,            /* embedding order            */
    pub istage: i32,       /* current stage              */
    pub stages: i32,       /* number of stages           */
    pub Be: Option<ARKodeButcherTable>, /* ERK Butcher table */
    pub Bi: Option<ARKodeButcherTable>, /* IRK Butcher table */

    /* User-supplied stage predictor routine */
    pub stage_predict: Option<ARKStagePredictFn>,

    /* (Non)Linear solver parameters & data */
    pub NLS: Option<NonlinearSolver>, /* generic SUNNonlinearSolver object */
    pub ownNLS: bool,                 /* flag indicating ownership of NLS  */
    pub nls_fi: Option<ARKRhsFn>,     /* fi(t,y) used in the nonlinear solver */
    pub gamma: f64,  /* gamma = h * A(i,i)                       */
    pub gammap: f64, /* gamma at the last setup call             */
    pub gamrat: f64, /* gamma / gammap                           */
    pub dgmax: f64,  /* call lsetup if |gamma/gammap-1| >= dgmax */

    pub predictor: i32, /* implicit prediction method to use        */
    pub crdown: f64,    /* nonlinear conv rate estimation constant  */
    pub rdiv: f64,      /* nonlin divergence if del/delp > rdiv     */
    pub crate_: f64,    /* estimated nonlin convergence rate
                        (C: `crate`; renamed — reserved word)    */
    pub delp: f64,      /* norm of previous nonlinear solver update */
    pub eRNrm: f64,     /* estimated residual norm, used in nonlin
                        and linear solver convergence tests      */
    pub nlscoef: f64,   /* coefficient in nonlin. convergence test  */

    pub msbp: i32,   /* positive => max # steps between lsetup
                     negative => call at each Newton iter     */
    pub nstlp: i64,  /* step number of last setup call           */

    pub maxcor: i32, /* max num iterations for solving the
                     nonlinear equation                       */

    pub convfail: i32, /* NLS fail flag (for interface routines)   */
    pub jcur: bool,    /* is Jacobian info for lin solver current? */
    /* (fn_implicit: alias to a saved implicit RHS evaluation —
       re-derived at the arkode_arkstep_nls.rs use sites) */

    /* Linear Solver Data (the lmem box itself lives on
       ARKodeMem.lmem — Addendum C.1) */
    pub linit: Option<ARKLinsolInitFn>,
    pub lsetup: Option<ARKLinsolSetupFn>,
    pub lsolve: Option<ARKLinsolSolveFn>,
    pub lfree: Option<ARKLinsolFreeFn>,
    pub lsolve_type: SUNLinearSolver_Type,

    /* Mass matrix solver data: deferred with the ARKLS mass half;
       mass_type kept (stepper logic branches on it) */
    pub mass_type: i32, /* 0=identity, 1=fixed, 2=time-dep */

    /* Counters */
    pub nfe: i64,       /* num fe calls               */
    pub nfi: i64,       /* num fi calls               */
    pub nsetups: i64,   /* num setup calls            */
    pub nls_iters: i64, /* num nonlinear solver iters */
    pub nls_fails: i64, /* num nonlinear solver fails */

    /* Reusable arrays for fused vector operations */
    pub cvals: Vec<f64>,   /* scalar array for fused ops       */
    pub nfusedopvecs: i32, /* length of cvals and Xvecs arrays */

    /* Data for using ARKStep with external polynomial forcing */
    pub expforcing: bool,       /* add forcing to explicit RHS */
    pub impforcing: bool,       /* add forcing to implicit RHS */
    pub tshift: f64,            /* time normalization shift    */
    pub tscale: f64,            /* time normalization scaling  */
    pub forcing: Vec<NVector>,  /* array of forcing vectors    */
    pub nforcing: i32,          /* number of forcing vectors   */
    pub stage_times: Vec<f64>,  /* workspace for applying forcing */
    pub stage_coefs: Vec<f64>,  /* workspace for applying forcing */
}
