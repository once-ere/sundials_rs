/* -----------------------------------------------------------------
 * Translated from src/kinsol/kinsol_impl.h and the constants /
 * typedefs of include/kinsol/kinsol.h (KINSOL 7.7.0).
 * Main solver memory block (struct KINMemRec) and internal
 * constants. Field names keep the C names (kin_ prefix) so the
 * translation of kinsol.c reads line-for-line.
 * -----------------------------------------------------------------*/
use crate::kinsol_ls_impl::KINLsMem;
use crate::nvector_serial::NVector;
use crate::sundials_context::SUNContext;
use crate::sundials_iterative::{SUNQRAddFn, SUNQRData};
use crate::sundials_math::{SUNRpowerR, SUNRsqrt};
use crate::sundials_types::{UserData, SUNFALSE, SUNTRUE, SUN_UNIT_ROUNDOFF};

/* ===============================================================
   User-supplied function types (kinsol.h)
   =============================================================== */

pub type KINSysFn = fn(uu: &NVector, fval: &mut NVector, user_data: &mut UserData) -> i32;

pub type KINDampingFn = fn(
    iter: i64,
    u_val: &NVector,
    g_val: &NVector,
    qt_fn: &[f64],
    depth: i64,
    user_data: &mut UserData,
    damping_factor: &mut f64,
) -> i32;

pub type KINDepthFn = fn(
    iter: i64,
    u_val: &NVector,
    g_val: &NVector,
    f_val: &NVector,
    df: &[NVector],
    r_mat: &[f64],
    depth: i64,
    user_data: &mut UserData,
    new_depth: &mut i64,
    remove_indices: &mut [bool],
) -> i32;

/* ===============================================================
   KINSOL constants (kinsol.h)
   =============================================================== */

/* return values */
pub const KIN_SUCCESS: i32 = 0;
pub const KIN_INITIAL_GUESS_OK: i32 = 1;
pub const KIN_STEP_LT_STPTOL: i32 = 2;

pub const KIN_WARNING: i32 = 99;

pub const KIN_MEM_NULL: i32 = -1;
pub const KIN_ILL_INPUT: i32 = -2;
pub const KIN_NO_MALLOC: i32 = -3;
pub const KIN_MEM_FAIL: i32 = -4;
pub const KIN_LINESEARCH_NONCONV: i32 = -5;
pub const KIN_MAXITER_REACHED: i32 = -6;
pub const KIN_MXNEWT_5X_EXCEEDED: i32 = -7;
pub const KIN_LINESEARCH_BCFAIL: i32 = -8;
pub const KIN_LINSOLV_NO_RECOVERY: i32 = -9;
pub const KIN_LINIT_FAIL: i32 = -10;
pub const KIN_LSETUP_FAIL: i32 = -11;
pub const KIN_LSOLVE_FAIL: i32 = -12;
pub const KIN_SYSFUNC_FAIL: i32 = -13;
pub const KIN_FIRST_SYSFUNC_ERR: i32 = -14;
pub const KIN_REPTD_SYSFUNC_ERR: i32 = -15;
pub const KIN_VECTOROP_ERR: i32 = -16;
pub const KIN_CONTEXT_ERR: i32 = -17;
pub const KIN_DAMPING_FN_ERR: i32 = -18;
pub const KIN_DEPTH_FN_ERR: i32 = -19;

/* Anderson Acceleration Orthogonalization Choice */
pub const KIN_ORTH_MGS: i32 = 0;
pub const KIN_ORTH_ICWY: i32 = 1;
pub const KIN_ORTH_CGS2: i32 = 2;
pub const KIN_ORTH_DCGS2: i32 = 3;

/* Enumeration for eta choice */
pub const KIN_ETACHOICE1: i32 = 1;
pub const KIN_ETACHOICE2: i32 = 2;
pub const KIN_ETACONSTANT: i32 = 3;

/* Enumeration for global strategy */
pub const KIN_NONE: i32 = 0;
pub const KIN_LINESEARCH: i32 = 1;
pub const KIN_PICARD: i32 = 2;
pub const KIN_FP: i32 = 3;

/* ===============================================================
   Internal constants (kinsol_impl.h) — KINSOL default constants
   =============================================================== */

pub const MXITER_DEFAULT: i64 = 200;
pub const MXNBCF_DEFAULT: i64 = 10;
pub const MSBSET_DEFAULT: i64 = 10;
pub const MSBSET_SUB_DEFAULT: i64 = 5;

pub const OMEGA_MIN: f64 = 0.00001;
pub const OMEGA_MAX: f64 = 0.9;

/* Shortcuts KIN_PROFILER / KIN_LOGGER are not ported: the workspace
   SUNContext is a unit-like struct because the default SUNDIALS
   build compiles solver profiling/logging macros out. */

/* ===============================================================
   Linear solver module attached to KINSOL.
   In C this is the (kin_linit, kin_lsetup, kin_lsolve, kin_lfree)
   function-pointer four-tuple plus the void* kin_lmem; here it is
   an enum dispatched in kinsol.rs (donor cvode_impl.rs pattern).
   KINSOL has a single interface module (KINLS, kinsol_ls.rs). The
   module is take()n out of KINMem for the duration of a call so
   its routines may borrow the solver memory mutably.
   =============================================================== */

#[derive(Default)]
pub enum LsModule {
    #[default]
    None,
    /// kinsol_ls.rs interface (KINSetLinearSolver)
    Ls(Box<KINLsMem>),
}

impl LsModule {
    pub fn is_none(&self) -> bool {
        matches!(self, LsModule::None)
    }
}

/* ===============================================================
   Main solver memory block (struct KINMemRec)
   =============================================================== */

pub struct KINMem {
    pub kin_sunctx: SUNContext,

    /* (C void* python — Python-binding hook — is excluded with the
       foreign-runtime backends.) */
    pub kin_uround: f64, /* machine epsilon (or unit roundoff error) */

    /* problem specification data */
    pub kin_func: Option<KINSysFn>, /* nonlinear system function implementation     */
    pub kin_user_data: UserData,    /* work space available to func routine         */
    pub kin_fnormtol: f64,          /* stopping tolerance on L2-norm of function
                                       value                                        */
    pub kin_scsteptol: f64,         /* scaled step length tolerance                 */
    pub kin_globalstrategy: i32,    /* choices are KIN_NONE, KIN_LINESEARCH
                                       KIN_PICARD and KIN_FP                        */
    pub kin_mxiter: i64,            /* maximum number of nonlinear iterations       */
    pub kin_msbset: i64,            /* maximum number of nonlinear iterations that
                                       may be performed between calls to the
                                       linear solver setup routine (lsetup)         */
    pub kin_msbset_sub: i64,        /* subinterval length for residual monitoring   */
    pub kin_mxnbcf: i64,            /* maximum number of beta condition failures    */
    pub kin_etaflag: i32,           /* choices are KIN_ETACONSTANT, KIN_ETACHOICE1
                                       and KIN_ETACHOICE2                           */
    pub kin_noMinEps: bool,         /* flag controlling whether or not the value
                                       of eps is bounded below                      */
    pub kin_constraintsSet: bool,   /* flag indicating if constraints are being
                                       used                                         */
    pub kin_jacCurrent: bool,       /* flag indicating if the Jacobian info.
                                       used by the linear solver is current         */
    pub kin_callForcingTerm: bool,  /* flag set if using either KIN_ETACHOICE1
                                       or KIN_ETACHOICE2                            */
    pub kin_noResMon: bool,         /* flag indicating if the nonlinear
                                       residual monitoring scheme should be
                                       used                                         */
    pub kin_retry_nni: bool,        /* flag indicating if nonlinear iteration
                                       should be retried (set by residual
                                       monitoring algorithm)                        */
    pub kin_update_fnorm_sub: bool, /* flag indicating if the fnorm associated
                                       with the subinterval needs to be
                                       updated (set by residual monitoring
                                       algorithm)                                   */

    pub kin_mxnewtstep: f64,   /* maximum allowable scaled step length         */
    pub kin_mxnstepin: f64,    /* input (or preset) value for mxnewtstep       */
    pub kin_sqrt_relfunc: f64, /* relative error bound for func(u)             */
    pub kin_stepl: f64,        /* scaled length of current step                */
    pub kin_stepmul: f64,      /* step scaling factor                          */
    pub kin_eps: f64,          /* current value of eps                         */
    pub kin_eta: f64,          /* current value of eta                         */
    pub kin_eta_gamma: f64,    /* gamma value used in eta calculation
                                  (choice #2)                                   */
    pub kin_eta_alpha: f64,    /* alpha value used in eta calculation
                                  (choice #2)                                   */
    pub kin_noInitSetup: bool, /* flag controlling whether or not the KINSol
                                  routine makes an initial call to the
                                  linear solver setup routine (lsetup)          */
    pub kin_sthrsh: f64,       /* threshold value for calling the linear
                                  solver setup routine                          */

    /* counters */
    pub kin_nni: i64,         /* number of nonlinear iterations               */
    pub kin_nfe: i64,         /* number of calls made to func routine         */
    pub kin_nnilset: i64,     /* value of nni counter when the linear solver
                                 setup was last called                         */
    pub kin_nnilset_sub: i64, /* value of nni counter when the linear solver
                                 setup was last called (subinterval)           */
    pub kin_nbcf: i64,        /* number of times the beta-condition could not
                                 be met in KINLineSearch                       */
    pub kin_nbktrk: i64,      /* number of backtracks performed by
                                 KINLineSearch                                 */
    pub kin_ncscmx: i64,      /* number of consecutive steps of size
                                 mxnewtstep taken                              */

    /* vectors */
    pub kin_uu: NVector,          /* solution vector/current iterate (initially
                                     contains initial guess, but holds approximate
                                     solution upon completion if no errors occurred) */
    pub kin_unew: NVector,        /* next iterate (unew = uu+pp)                  */
    pub kin_fval: NVector,        /* vector containing result of nonlinear system
                                     function evaluated at a given iterate
                                     (fval = func(uu))                             */
    pub kin_gval: NVector,        /* vector containing result of the fixed point
                                     function evaluated at a given iterate;
                                     used in KIN_PICARD strategy only.
                                     (gval = uu - L^{-1}fval(uu))                  */
    pub kin_uscale: NVector,      /* iterate scaling vector                       */
    pub kin_fscale: NVector,      /* fval scaling vector                          */
    pub kin_pp: NVector,          /* incremental change vector (pp = unew-uu)     */
    pub kin_constraints: NVector, /* constraints vector                           */
    pub kin_vtemp1: NVector,      /* scratch vector #1                            */
    pub kin_vtemp2: NVector,      /* scratch vector #2                            */
    pub kin_vtemp3: NVector,      /* scratch vector #3                            */

    /* fixed point and Picard options */
    pub kin_ret_newest: bool, /* return the newest FP iteration     */
    pub kin_damping: bool,    /* flag to apply damping in FP/Picard */
    pub kin_beta: f64,        /* damping parameter for FP/Picard    */

    /* space requirements for AA, Broyden and NLEN */
    pub kin_fold_aa: NVector,     /* vector needed for AA, Broyden, and NLEN       */
    pub kin_gold_aa: NVector,     /* vector needed for AA, Broyden, and NLEN       */
    pub kin_df_aa: Vec<NVector>,  /* vector array needed for AA, Broyden, and NLEN */
    pub kin_dg_aa: Vec<NVector>,  /* vector array needed for AA, Broyden and NLEN  */
    pub kin_q_aa: Vec<NVector>,   /* vector array needed for AA                    */
    pub kin_beta_aa: f64,         /* beta damping parameter for AA                 */
    pub kin_gamma_aa: Vec<f64>,   /* array of size maa used in AA                  */
    pub kin_R_aa: Vec<f64>,       /* array of size maa*maa used in AA              */
    pub kin_T_aa: Vec<f64>,       /* array of size maa*maa used in AA with ICWY MGS */
    pub kin_m_aa: i64,            /* parameter for AA, Broyden or NLEN             */
    pub kin_m_aa_alloc: i64,      /* depth (m) used for AA memory allocations      */
    pub kin_delay_aa: i64,        /* number of iterations to delay AA              */
    pub kin_current_depth: i64,   /* current Anderson acceleration space size      */
    pub kin_damping_fn: Option<KINDampingFn>, /* function to determine the damping factor */
    pub kin_depth_fn: Option<KINDepthFn>,     /* function to determine the depth with AA  */
    pub kin_orth_aa: i32,         /* parameter for AA determining orthogonalization
                                     routine
                                     0 - Modified Gram Schmidt (standard)
                                     1 - ICWY Modified Gram Schmidt (Bjorck)
                                     2 - CGS2 (Hernandez)
                                     3 - Delayed CGS2 (Hernandez)                  */
    pub kin_orth_aa_alloc: i64,   /* depth (m) used for orthogonalization memory
                                     allocations                                   */
    pub kin_qr_func: Option<SUNQRAddFn>, /* QRAdd function for AA orthogonalization */
    pub kin_qr_data: Option<SUNQRData>,  /* Additional parameters required for QRAdd
                                            routine set for AA                      */
    pub kin_damping_aa: bool,     /* flag to apply damping in AA                   */
    pub kin_dot_prod_sb: bool,    /* use single buffer dot product                 */
    /* (C kin_cv / kin_Xv — scalar/vector alias arrays feeding the fused
       N_VLinearCombination kernels — are not stored: the serial fused
       kernels are reproduced inline, as in the donor's drop of
       cv_cvals/cv_Xvecs.) */

    /* space requirements for vector storage */
    pub kin_lrw1: i64, /* number of sunrealtype-sized memory blocks needed
                          for a single N_Vector                           */
    pub kin_liw1: i64, /* number of int-sized memory blocks needed for
                          a single N_Vecotr                               */
    pub kin_lrw: i64,  /* total number of sunrealtype-sized memory blocks
                          needed for all KINSOL work vectors              */
    pub kin_liw: i64,  /* total number of int-sized memory blocks needed
                          for all KINSOL work vectors                     */

    /* linear solver data */
    pub kin_inexact_ls: bool, /* flag set by the linear solver module
                                 (in linit) indicating whether this is an
                                 iterative linear solver (SUNTRUE), or a direct
                                 linear solver (SUNFALSE)                       */

    /* In C: kin_linit/kin_lsetup/kin_lsolve/kin_lfree function
       pointers + void* kin_lmem. Here: dispatching enum (see
       LsModule). */
    pub kin_lmem: LsModule, /* linear solver memory block */

    pub kin_fnorm: f64,   /* value of L2-norm of fscale*fval                   */
    pub kin_f1norm: f64,  /* f1norm = 0.5*(fnorm)^2                            */
    pub kin_sFdotJp: f64, /* value of scaled F(u) vector (fscale*fval)
                             dotted with scaled J(u)*pp vector (set by lsolve) */
    pub kin_sJpnorm: f64, /* value of L2-norm of fscale*(J(u)*pp)
                             (set by lsolve)                                   */

    pub kin_fnorm_sub: f64,   /* value of L2-norm of fscale*fval (subinterval) */
    pub kin_eval_omega: bool, /* flag indicating that omega must be evaluated. */
    pub kin_omega: f64,       /* constant value for real scalar used in test to
                                 determine if reduction of norm of nonlinear
                                 residual is sufficient. Unless a valid constant
                                 value is specified by the user, omega is estimated
                                 from omega_min and omega_max at each iteration.    */
    pub kin_omega_min: f64,   /* lower bound on omega                          */
    pub kin_omega_max: f64,   /* upper bound on omega                          */

    /*
     * -----------------------------------------------------------------
     * Note: The KINLineSearch subroutine scales the values of the
     * variables sFdotJp and sJpnorm by a factor rl (lambda) that is
     * chosen by the line search algorithm such that the sclaed Newton
     * step satisfies the following conditions:
     *
     *  F(u_k+1) <= F(u_k) + alpha*(F(u_k)^T * J(u_k))*p*rl
     *
     *  F(u_k+1) >= F(u_k) + beta*(F(u_k)^T * J(u_k))*p*rl
     *
     * where alpha = 1.0e-4, beta = 0.9, u_k+1 = u_k + rl*p,
     * 0 < rl <= 1, J denotes the system Jacobian, and F represents
     * the nonliner system function.
     * -----------------------------------------------------------------
     */
    pub kin_MallocDone: bool, /* flag indicating if KINMalloc has been
                                 called yet                                    */
}

/* ===============================================================
   Default construction: the C KINCreate memset(0)s KINMemRec and
   then sets the default optional inputs (kinsol.c, KINCreate).
   KINCreate in kinsol.rs builds on this Default exactly like the
   donor's CVodeCreate builds its struct literal.
   =============================================================== */
impl Default for KINMem {
    fn default() -> Self {
        let uround = SUN_UNIT_ROUNDOFF;
        KINMem {
            kin_sunctx: SUNContext::default(),
            kin_uround: uround,

            /* default values for solver optional inputs (KINCreate) */
            kin_func: None,
            kin_user_data: None,
            kin_fnormtol: SUNRpowerR(uround, 0.3333333333333333), /* ONETHIRD  */
            kin_scsteptol: SUNRpowerR(uround, 0.6666666666666667), /* TWOTHIRDS */
            kin_globalstrategy: KIN_NONE,
            kin_mxiter: MXITER_DEFAULT,
            kin_msbset: MSBSET_DEFAULT,
            kin_msbset_sub: MSBSET_SUB_DEFAULT,
            kin_mxnbcf: MXNBCF_DEFAULT,
            kin_etaflag: KIN_ETACHOICE1,
            kin_noMinEps: SUNFALSE,
            kin_constraintsSet: SUNFALSE,
            kin_jacCurrent: SUNFALSE,
            kin_callForcingTerm: SUNFALSE,
            kin_noResMon: SUNFALSE,
            kin_retry_nni: SUNFALSE,
            kin_update_fnorm_sub: SUNFALSE,

            kin_mxnewtstep: 0.0,
            kin_mxnstepin: 0.0,                  /* ZERO   */
            kin_sqrt_relfunc: SUNRsqrt(uround),
            kin_stepl: 0.0,
            kin_stepmul: 0.0,
            kin_eps: 0.0,
            kin_eta: 0.1,       /* POINT1: default for KIN_ETACONSTANT */
            kin_eta_gamma: 0.9, /* POINT9: default for KIN_ETACHOICE2  */
            kin_eta_alpha: 2.0, /* TWO:    default for KIN_ETACHOICE2  */
            kin_noInitSetup: SUNFALSE,
            kin_sthrsh: 2.0, /* TWO */

            kin_nni: 0,
            kin_nfe: 0,
            kin_nnilset: 0,
            kin_nnilset_sub: 0,
            kin_nbcf: 0,
            kin_nbktrk: 0,
            kin_ncscmx: 0,

            kin_uu: NVector::default(),
            kin_unew: NVector::default(),
            kin_fval: NVector::default(),
            kin_gval: NVector::default(),
            kin_uscale: NVector::default(),
            kin_fscale: NVector::default(),
            kin_pp: NVector::default(),
            kin_constraints: NVector::default(),
            kin_vtemp1: NVector::default(),
            kin_vtemp2: NVector::default(),
            kin_vtemp3: NVector::default(),

            kin_ret_newest: SUNFALSE,
            kin_damping: SUNFALSE,
            kin_beta: 1.0, /* ONE */

            kin_fold_aa: NVector::default(),
            kin_gold_aa: NVector::default(),
            kin_df_aa: Vec::new(),
            kin_dg_aa: Vec::new(),
            kin_q_aa: Vec::new(),
            kin_beta_aa: 1.0, /* ONE */
            kin_gamma_aa: Vec::new(),
            kin_R_aa: Vec::new(),
            kin_T_aa: Vec::new(),
            kin_m_aa: 0,
            kin_m_aa_alloc: 0,
            kin_delay_aa: 0,
            kin_current_depth: 0,
            kin_damping_fn: None,
            kin_depth_fn: None,
            kin_orth_aa: KIN_ORTH_MGS,
            kin_orth_aa_alloc: 0,
            kin_qr_func: None,
            kin_qr_data: None,
            kin_damping_aa: SUNFALSE,
            kin_dot_prod_sb: SUNFALSE,

            /* NOTE: needed since KINInit could be called after
               KINSetConstraints */
            kin_lrw1: 0,
            kin_liw1: 0,
            kin_lrw: 17,
            kin_liw: 22,

            kin_inexact_ls: SUNFALSE,
            kin_lmem: LsModule::None,

            kin_fnorm: 0.0,
            kin_f1norm: 0.0,
            kin_sFdotJp: 0.0,
            kin_sJpnorm: 0.0,

            kin_fnorm_sub: 0.0,
            kin_eval_omega: SUNTRUE,
            kin_omega: 0.0, /* ZERO: default to using min/max */
            kin_omega_min: OMEGA_MIN,
            kin_omega_max: OMEGA_MAX,

            kin_MallocDone: SUNFALSE,
        }
    }
}

/* ===============================================================
   KINSOL internal functions (prototypes in kinsol_impl.h; the
   implementations live in kinsol.rs):

     KINPrintInfo(kin_mem, info_code, module, fname, ...)
         — high level info handler (varargs formatting is done at
           each call site; output goes to &mut dyn std::io::Write
           per the workspace convention)
     KINInfoHandler(module, function, msg, user_data)
         — internal infoHandler function
     KINInitAA(kin_mem) / KINFreeAA(kin_mem)
         — Anderson acceleration utilities
     KINInitOrth(kin_mem) / KINFreeOrth(kin_mem)
         — orthogonalization utilities

   KINProcessError (high level error handler) is defined below,
   following the donor's cvProcessError.
   kinsol_user_supplied_fn_table_destroy (SUNDIALS_ENABLE_PYTHON) is
   excluded with the Python binding backend.
   =============================================================== */

/* ===============================================================
   Error handler: the C KINProcessError routes printf-style messages
   to the SUNContext error handler stack; here messages go to stderr
   (equivalent to the default SUNLogErrHandlerFn behavior).
   =============================================================== */
pub fn KINProcessError(
    _kin_mem: Option<&KINMem>,
    error_code: i32,
    line: u32,
    func: &str,
    file: &str,
    msg: &str,
) {
    if error_code == KIN_WARNING {
        eprintln!("\n[KINSOL WARNING] {file}:{line} in {func}\n  {msg}\n");
    } else {
        eprintln!("\n[KINSOL ERROR] {file}:{line} in {func}\n  {msg}\n");
    }
}

/* ===============================================================
   KINSOL error messages (kinsol_impl.h)
   =============================================================== */

pub const MSG_MEM_FAIL: &str = "A memory request failed.";
pub const MSG_NO_MEM: &str = "kinsol_mem = NULL illegal.";
pub const MSG_NULL_SUNCTX: &str = "sunctx = NULL illegal.";
pub const MSG_BAD_NVECTOR: &str = "A required vector operation is not implemented.";
pub const MSG_FUNC_NULL: &str = "func = NULL illegal.";
pub const MSG_NO_MALLOC: &str = "Attempt to call before KINMalloc illegal.";

pub const MSG_BAD_MXITER: &str = "Illegal value for mxiter.";
pub const MSG_BAD_MSBSET: &str = "Illegal msbset < 0.";
pub const MSG_BAD_MSBSETSUB: &str = "Illegal msbsetsub < 0.";
pub const MSG_BAD_ETACHOICE: &str = "Illegal value for etachoice.";
pub const MSG_BAD_ETACONST: &str = "eta out of range.";
pub const MSG_BAD_GAMMA: &str = "gamma out of range.";
pub const MSG_BAD_ALPHA: &str = "alpha out of range.";
pub const MSG_BAD_MXNEWTSTEP: &str = "Illegal mxnewtstep < 0.";
pub const MSG_BAD_RELFUNC: &str = "relfunc < 0 illegal.";
pub const MSG_BAD_FNORMTOL: &str = "fnormtol < 0 illegal.";
pub const MSG_BAD_SCSTEPTOL: &str = "scsteptol < 0 illegal.";
pub const MSG_BAD_MXNBCF: &str = "mxbcf < 0 illegal.";
pub const MSG_BAD_CONSTRAINTS: &str = "Illegal values in constraints vector.";
pub const MSG_BAD_OMEGA: &str = "scalars < 0 illegal.";
pub const MSG_BAD_MAA: &str = "maa < 0 illegal.";
pub const MSG_BAD_ORTHAA: &str = "Illegal value for orthaa.";
pub const MSG_ZERO_MAA: &str = "maa = 0 illegal.";

pub const MSG_LSOLV_NO_MEM: &str = "The linear solver memory pointer is NULL.";
pub const MSG_UU_NULL: &str = "uu = NULL illegal.";
pub const MSG_BAD_GLSTRAT: &str = "Illegal value for global strategy.";
pub const MSG_BAD_USCALE: &str = "uscale = NULL illegal.";
pub const MSG_USCALE_NONPOSITIVE: &str = "uscale has nonpositive elements.";
pub const MSG_BAD_FSCALE: &str = "fscale = NULL illegal.";
pub const MSG_FSCALE_NONPOSITIVE: &str = "fscale has nonpositive elements.";
pub const MSG_CONSTRAINTS_NOTOK: &str =
    "Constraints not allowed with fixed point or Picard iterations";
pub const MSG_INITIAL_CNSTRNT: &str = "Initial guess does NOT meet constraints.";
pub const MSG_LINIT_FAIL: &str = "The linear solver's init routine failed.";

pub const MSG_SYSFUNC_FAILED: &str =
    "The system function failed in an unrecoverable manner.";
pub const MSG_SYSFUNC_FIRST: &str = "The system function failed at the first call.";
pub const MSG_LSETUP_FAILED: &str =
    "The linear solver's setup function failed in an unrecoverable manner.";
pub const MSG_LSOLVE_FAILED: &str =
    "The linear solver's solve function failed in an unrecoverable manner.";
pub const MSG_LINSOLV_NO_RECOVERY: &str =
    "The linear solver's solve function failed recoverably, but the Jacobian \
     data is already current.";
pub const MSG_LINESEARCH_NONCONV: &str =
    "The line search algorithm was unable to find an iterate sufficiently \
     distinct from the current iterate.";
pub const MSG_LINESEARCH_BCFAIL: &str =
    "The line search algorithm was unable to satisfy the beta-condition for \
     nbcfails iterations.";
pub const MSG_MAXITER_REACHED: &str =
    "The maximum number of iterations was reached before convergence.";
pub const MSG_MXNEWT_5X_EXCEEDED: &str =
    "Five consecutive steps have been taken that satisfy a scaled step length \
     test.";
pub const MSG_SYSFUNC_REPTD: &str =
    "Unable to correct repeated recoverable system function errors.";
pub const MSG_NOL_FAIL: &str =
    "Unable to find user's Linear Jacobian, which is required for the \
     KIN_PICARD Strategy";

/* ===============================================================
   KINSOL info messages (kinsol_impl.h); C printf formats kept
   verbatim with SUN_FORMAT_G / SUN_FORMAT_E expanded to their
   double-precision forms "%.15g" / "% .15e" (KINPrintInfo renders
   them via sundials_utils::fmt_g / fmt_e).
   =============================================================== */

pub const INFO_IVAR: &str = "%s = %d";
pub const INFO_LIVAR: &str = "%s = %ld";
pub const INFO_RETVAL: &str = "Return value: %d";
pub const INFO_ADJ: &str = "no. of lambda adjustments = %ld";

pub const INFO_RVAR: &str = "%s = %.15g";
pub const INFO_NNI: &str = "nni = %4ld, nfe = %6ld, fnorm = %.15g";
pub const INFO_TOL: &str = "scsteptol = %.15g, fnormtol = %.15g";
pub const INFO_FMAX: &str = "scaled f norm (for stopping) = %.15g";
pub const INFO_PNORM: &str = "pnorm = % .15e";
pub const INFO_PNORM1: &str = "(ivio=1) pnorm = % .15e";
pub const INFO_FNORM: &str = "fnorm(L2) = % .15e";
pub const INFO_LAM: &str = "min_lam = % .15e, f1norm = % .15e, pnorm = % .15e";
pub const INFO_ALPHA: &str =
    "fnorm = % .15e, f1norm = % .15e, alpha_cond = % .15e, lam = % .15e";
pub const INFO_BETA: &str = "f1norm = % .15e, beta_cond = % .15e, lam = % .15e";
pub const INFO_ALPHABETA: &str =
    "f1norm = % .15e, alpha_cond = % .15e, beta_cond = % .15e, lam = % .15e";

#[cfg(test)]
mod tests {
    use super::*;

    /* A Default KINMem carries the C KINCreate default values
       (kinsol.c, KINCreate; constants from kinsol_impl.h). */
    #[test]
    fn kinmem_default_matches_kincreate() {
        let kin_mem = KINMem::default();

        assert_eq!(kin_mem.kin_mxiter, MXITER_DEFAULT); /* 200 */
        assert_eq!(kin_mem.kin_mxiter, 200);
        assert_eq!(kin_mem.kin_msbset, MSBSET_DEFAULT); /* 10 */
        assert_eq!(kin_mem.kin_msbset, 10);
        assert_eq!(kin_mem.kin_msbset_sub, MSBSET_SUB_DEFAULT); /* 5 */
        assert_eq!(kin_mem.kin_mxnbcf, MXNBCF_DEFAULT); /* 10 */

        assert_eq!(kin_mem.kin_omega, 0.0); /* default to using min/max */
        assert_eq!(kin_mem.kin_omega_min, OMEGA_MIN);
        assert_eq!(kin_mem.kin_omega_min, 0.00001);
        assert_eq!(kin_mem.kin_omega_max, OMEGA_MAX);
        assert_eq!(kin_mem.kin_omega_max, 0.9);
        assert!(kin_mem.kin_eval_omega);

        assert_eq!(kin_mem.kin_globalstrategy, KIN_NONE);
        assert_eq!(kin_mem.kin_etaflag, KIN_ETACHOICE1);
        assert_eq!(kin_mem.kin_eta, 0.1);
        assert_eq!(kin_mem.kin_eta_alpha, 2.0);
        assert_eq!(kin_mem.kin_eta_gamma, 0.9);
        assert_eq!(kin_mem.kin_sthrsh, 2.0);
        assert_eq!(kin_mem.kin_beta, 1.0);
        assert_eq!(kin_mem.kin_beta_aa, 1.0);
        assert_eq!(kin_mem.kin_orth_aa, KIN_ORTH_MGS);

        assert_eq!(kin_mem.kin_uround, SUN_UNIT_ROUNDOFF);
        assert_eq!(kin_mem.kin_sqrt_relfunc, SUNRsqrt(SUN_UNIT_ROUNDOFF));
        assert_eq!(
            kin_mem.kin_scsteptol,
            SUNRpowerR(SUN_UNIT_ROUNDOFF, 0.6666666666666667)
        );
        assert_eq!(
            kin_mem.kin_fnormtol,
            SUNRpowerR(SUN_UNIT_ROUNDOFF, 0.3333333333333333)
        );

        assert_eq!(kin_mem.kin_lrw, 17);
        assert_eq!(kin_mem.kin_liw, 22);
        assert_eq!(kin_mem.kin_lrw1, 0);
        assert_eq!(kin_mem.kin_liw1, 0);

        assert!(kin_mem.kin_lmem.is_none());
        assert!(!kin_mem.kin_MallocDone);
        assert!(!kin_mem.kin_constraintsSet);
        assert!(!kin_mem.kin_ret_newest);
        assert!(!kin_mem.kin_damping);
        assert!(!kin_mem.kin_damping_aa);
        assert!(!kin_mem.kin_dot_prod_sb);
        assert_eq!(kin_mem.kin_m_aa, 0);
        assert_eq!(kin_mem.kin_delay_aa, 0);
        assert_eq!(kin_mem.kin_current_depth, 0);
    }
}
