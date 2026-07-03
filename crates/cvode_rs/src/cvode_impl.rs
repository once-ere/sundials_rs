/* -----------------------------------------------------------------
 * Translated from src/cvode/cvode_impl.h (CVODE 7.7.0).
 * Main integrator memory block (struct CVodeMemRec) and internal
 * constants. Field names keep the C names (cv_ prefix) so the
 * translation of cvode.c reads line-for-line.
 * -----------------------------------------------------------------*/
use crate::cvode_ls_impl::CVLsMem;
use crate::cvode_diag_impl::CVDiagMem;
use crate::cvode_proj_impl::CVodeProjMem;
use crate::nvector_serial::NVector;
use crate::sundials_context::SUNContext;
use crate::sundials_nonlinearsolver::NonlinearSolver;
use crate::sundials_types::UserData;

/* ===============================================================
   User-supplied function types (cvode.h)
   =============================================================== */

pub type CVRhsFn = fn(t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32;
pub type CVRootFn = fn(t: f64, y: &NVector, gout: &mut [f64], user_data: &mut UserData) -> i32;
pub type CVEwtFn = fn(y: &NVector, ewt: &mut NVector, user_data: &mut UserData) -> i32;

/* ===============================================================
   CVODE constants (cvode.h)
   =============================================================== */

/* lmm */
pub const CV_ADAMS: i32 = 1;
pub const CV_BDF: i32 = 2;

/* itask */
pub const CV_NORMAL: i32 = 1;
pub const CV_ONE_STEP: i32 = 2;

/* itol (internal) */
pub const CV_NN: i32 = 0;
pub const CV_SS: i32 = 1;
pub const CV_SV: i32 = 2;
pub const CV_WF: i32 = 3;

/* return values */
pub const CV_SUCCESS: i32 = 0;
pub const CV_TSTOP_RETURN: i32 = 1;
pub const CV_ROOT_RETURN: i32 = 2;

pub const CV_WARNING: i32 = 99;

pub const CV_TOO_MUCH_WORK: i32 = -1;
pub const CV_TOO_MUCH_ACC: i32 = -2;
pub const CV_ERR_FAILURE: i32 = -3;
pub const CV_CONV_FAILURE: i32 = -4;

pub const CV_LINIT_FAIL: i32 = -5;
pub const CV_LSETUP_FAIL: i32 = -6;
pub const CV_LSOLVE_FAIL: i32 = -7;
pub const CV_RHSFUNC_FAIL: i32 = -8;
pub const CV_FIRST_RHSFUNC_ERR: i32 = -9;
pub const CV_REPTD_RHSFUNC_ERR: i32 = -10;
pub const CV_UNREC_RHSFUNC_ERR: i32 = -11;
pub const CV_RTFUNC_FAIL: i32 = -12;
pub const CV_NLS_INIT_FAIL: i32 = -13;
pub const CV_NLS_SETUP_FAIL: i32 = -14;
pub const CV_CONSTR_FAIL: i32 = -15;
pub const CV_NLS_FAIL: i32 = -16;

pub const CV_MEM_FAIL: i32 = -20;
pub const CV_MEM_NULL: i32 = -21;
pub const CV_ILL_INPUT: i32 = -22;
pub const CV_NO_MALLOC: i32 = -23;
pub const CV_BAD_K: i32 = -24;
pub const CV_BAD_T: i32 = -25;
pub const CV_BAD_DKY: i32 = -26;
pub const CV_TOO_CLOSE: i32 = -27;
pub const CV_VECTOROP_ERR: i32 = -28;

pub const CV_PROJ_MEM_NULL: i32 = -29;
pub const CV_PROJFUNC_FAIL: i32 = -30;
pub const CV_REPTD_PROJFUNC_ERR: i32 = -31;

pub const CV_CONTEXT_ERR: i32 = -32;

pub const CV_UNRECOGNIZED_ERR: i32 = -99;

/* ===============================================================
   Internal constants (cvode_impl.h)
   =============================================================== */

pub const ADAMS_Q_MAX: usize = 12; /* max value of q for lmm == ADAMS */
pub const BDF_Q_MAX: usize = 5; /* max value of q for lmm == BDF   */
pub const Q_MAX: usize = ADAMS_Q_MAX; /* max value of q for either lmm   */
pub const L_MAX: usize = Q_MAX + 1; /* max value of L for either lmm   */
pub const NUM_TESTS: usize = 5; /* number of error test quantities */

pub const HMIN_DEFAULT: f64 = 0.0;
pub const HMAX_INV_DEFAULT: f64 = 0.0;
pub const MXHNIL_DEFAULT: i32 = 10;
pub const MXSTEP_DEFAULT: i64 = 500;

pub const MSBP_DEFAULT: i64 = 20; /* max steps between lsetup calls */
pub const DGMAX_LSETUP_DEFAULT: f64 = 0.3; /* gamma threshold to call lsetup */

/* Step size change constants */
pub const ETA_MIN_FX_DEFAULT: f64 = 0.0;
pub const ETA_MAX_FX_DEFAULT: f64 = 1.5;
pub const ETA_MAX_FS_DEFAULT: f64 = 10000.0;
pub const ETA_MAX_ES_DEFAULT: f64 = 10.0;
pub const ETA_MAX_GS_DEFAULT: f64 = 10.0;
pub const ETA_MIN_DEFAULT: f64 = 0.1;
pub const ETA_MAX_EF_DEFAULT: f64 = 0.2;
pub const ETA_MIN_EF_DEFAULT: f64 = 0.1;
pub const ETA_CF_DEFAULT: f64 = 0.25;
pub const SMALL_NST_DEFAULT: i64 = 10;
pub const SMALL_NEF_DEFAULT: i32 = 2;
pub const ONEPSM: f64 = 1.000001;

/* Step size controller constants */
pub const ADDON: f64 = 0.000001;
pub const BIAS1: f64 = 6.0;
pub const BIAS2: f64 = 6.0;
pub const BIAS3: f64 = 10.0;

/* Order selection constants */
pub const LONG_WAIT: i32 = 10;

/* Failure limits */
pub const MXNCF: i32 = 10;
pub const MXNEF: i32 = 7;
pub const MXNEF1: i32 = 3;
pub const MAX_CONSTRAINT_FAILS: i32 = 10;

/* Control constants for lower-level functions used by cvStep */
pub const DO_ERROR_TEST: i32 = 2;
pub const PREDICT_AGAIN: i32 = 3;

pub const TRY_AGAIN: i32 = 5;
pub const FIRST_CALL: i32 = 6;
pub const PREV_CONV_FAIL: i32 = 7;
pub const PREV_PROJ_FAIL: i32 = 8;
pub const PREV_ERR_FAIL: i32 = 9;

pub const RHSFUNC_RECVR: i32 = 10;
pub const CONSTRFUNC_RECVR: i32 = 11;
pub const PROJFUNC_RECVR: i32 = 12;

/* Constants for convfail (input to cv_lsetup) */
pub const CV_NO_FAILURES: i32 = 0;
pub const CV_FAIL_BAD_J: i32 = 1;
pub const CV_FAIL_OTHER: i32 = 2;

/* ===============================================================
   Linear solver module attached to CVODE.
   In C this is the (cv_linit, cv_lsetup, cv_lsolve, cv_lfree)
   function-pointer four-tuple plus the void* cv_lmem; here it is an
   enum dispatched in cvode.rs. The module is take()n out of
   CVodeMem for the duration of a call so its routines may borrow
   the integrator memory mutably.
   =============================================================== */

#[derive(Default)]
pub enum LsModule {
    #[default]
    None,
    /// cvode_ls.rs interface (CVodeSetLinearSolver)
    Ls(Box<CVLsMem>),
    /// cvode_diag.rs diagonal approximation (CVDiag)
    Diag(Box<CVDiagMem>),
}

impl LsModule {
    pub fn is_none(&self) -> bool {
        matches!(self, LsModule::None)
    }
}

/* ===============================================================
   Main integrator memory block (struct CVodeMemRec)
   =============================================================== */

pub struct CVodeMem {
    pub cv_sunctx: SUNContext,

    pub cv_uround: f64, /* machine unit roundoff */

    /*--------------------------
      Problem Specification Data
      --------------------------*/
    pub cv_f: Option<CVRhsFn>, /* y' = f(t,y(t))                    */
    pub cv_user_data: UserData, /* user pointer passed to f          */
    pub cv_lmm: i32,           /* lmm = CV_ADAMS or CV_BDF          */
    pub cv_itol: i32,          /* itol = CV_SS, CV_SV, CV_WF, CV_NN */

    pub cv_reltol: f64,          /* relative tolerance                */
    pub cv_Sabstol: f64,         /* scalar absolute tolerance         */
    pub cv_Vabstol: NVector,     /* vector absolute tolerance         */
    pub cv_atolmin0: bool,       /* flag indicating that min(abstol) = 0 */
    pub cv_user_efun: bool,      /* SUNTRUE if user sets efun         */
    pub cv_efun: Option<CVEwtFn>, /* function to set ewt              */

    /*-----------------------
      Nordsieck History Array
      -----------------------*/
    pub cv_zn: Vec<NVector>, /* Nordsieck array, of size N x (q+1) */

    /*-------------------
      Vectors of length N
      -------------------*/
    pub cv_ewt: NVector,    /* error weight vector                */
    pub cv_y: NVector,      /* temporary solver storage           */
    pub cv_acor: NVector,   /* accumulated correction / local err */
    pub cv_tempv: NVector,  /* temporary storage vector           */
    pub cv_ftemp: NVector,  /* temporary storage vector           */
    pub cv_vtemp1: NVector, /* temporary storage vector           */
    pub cv_vtemp2: NVector, /* temporary storage vector           */
    pub cv_vtemp3: NVector, /* temporary storage vector           */

    /*-----------------
      Tstop information
      -----------------*/
    pub cv_tstopset: bool,
    pub cv_tstopinterp: bool,
    pub cv_tstop: f64,

    /*---------
      Step Data
      ---------*/
    pub cv_q: i32,      /* current order                         */
    pub cv_qprime: i32, /* order to be used on the next step     */
    pub cv_next_q: i32, /* order to be used on the next step     */
    pub cv_qwait: i32,  /* steps to wait before order change     */
    pub cv_L: i32,      /* L = q + 1                             */

    pub cv_hin: f64,      /* initial step size                     */
    pub cv_h: f64,        /* current step size                     */
    pub cv_hprime: f64,   /* step size to be used on the next step */
    pub cv_next_h: f64,   /* step size to be used on the next step */
    pub cv_eta: f64,      /* eta = hprime / h                      */
    pub cv_hscale: f64,   /* value of h used in zn                 */
    pub cv_tn: f64,       /* current internal value of t           */
    pub cv_tretlast: f64, /* last value of t returned by CVode     */

    pub cv_tau: [f64; L_MAX + 1], /* previous q+1 successful step sizes, 1-indexed */
    pub cv_tq: [f64; NUM_TESTS + 1], /* test quantities, 1-indexed  */
    pub cv_l: [f64; L_MAX],       /* coefficients of l(x) (degree q poly)  */

    pub cv_rl1: f64,    /* the scalar 1/l[1]            */
    pub cv_gamma: f64,  /* gamma = h * rl1              */
    pub cv_gammap: f64, /* gamma at the last setup call */
    pub cv_gamrat: f64, /* gamma / gammap               */

    pub cv_crate: f64,    /* estimated corrector convergence rate    */
    pub cv_delp: f64,     /* norm of previous nonlinear solver update */
    pub cv_acnrm: f64,    /* | acor |                                */
    pub cv_acnrmcur: bool, /* is | acor | current?                   */
    pub cv_nlscoef: f64,  /* coefficient in nonlinear convergence test */

    /*------
      Limits
      ------*/
    pub cv_qmax: i32,   /* q <= qmax                                    */
    pub cv_mxstep: i64, /* max internal steps for one user call         */
    pub cv_mxhnil: i32, /* max warning messages for t + h == t          */
    pub cv_maxnef: i32, /* max error test failures                      */
    pub cv_maxncf: i32, /* max nonlinear convergence failures           */

    pub cv_hmin: f64,       /* |h| >= hmin        */
    pub cv_hmax_inv: f64,   /* |h| <= 1/hmax_inv  */
    pub cv_etamax: f64,     /* eta <= etamax      */
    pub cv_eta_min_fx: f64, /* eta_min_fx < eta < eta_max_fx keep h */
    pub cv_eta_max_fx: f64,
    pub cv_eta_max_fs: f64, /* eta <= eta_max_fs on the first step  */
    pub cv_eta_max_es: f64, /* eta <= eta_max_es on early steps     */
    pub cv_eta_max_gs: f64, /* eta <= eta_max_gs on a general step  */
    pub cv_eta_min: f64,    /* eta >= eta_min on a general step     */
    pub cv_eta_min_ef: f64, /* eta >= eta_min_ef after an error test failure */
    pub cv_eta_max_ef: f64, /* eta on multiple error test failures  */
    pub cv_eta_cf: f64,     /* eta on a nonlinear solver convergence failure */

    pub cv_small_nst: i64, /* nst <= small_nst use eta_max_es */
    pub cv_small_nef: i32, /* nef >= small_nef use eta_max_ef */

    /*--------
      Counters
      --------*/
    pub cv_nst: i64,     /* number of internal steps taken           */
    pub cv_nfe: i64,     /* number of f calls                        */
    pub cv_ncfn: i64,    /* number of corrector convergence failures */
    pub cv_nni: i64,     /* number of nonlinear iterations performed */
    pub cv_nnf: i64,     /* number of nonlinear convergence failures */
    pub cv_netf: i64,    /* number of error test failures            */
    pub cv_nsetups: i64, /* number of setup calls                    */
    pub cv_nhnil: i32,   /* number of t + h == t messages issued     */

    /*----------------
      Step size ratios
      ----------------*/
    pub cv_etaqm1: f64, /* ratio of new to old h for order q-1 */
    pub cv_etaq: f64,   /* ratio of new to old h for order q   */
    pub cv_etaqp1: f64, /* ratio of new to old h for order q+1 */

    /*------------------
      Space requirements
      ------------------*/
    pub cv_lrw1: i64, /* no. of sunrealtype words in 1 N_Vector y  */
    pub cv_liw1: i64, /* no. of integer words in 1 N_Vector y      */
    pub cv_lrw: i64,  /* no. of sunrealtype words in CVODE vectors */
    pub cv_liw: i64,  /* no. of integer words in CVODE vectors     */

    /*---------------------
      Nonlinear Solver Data
      ---------------------*/
    pub NLS: Option<NonlinearSolver>, /* nonlinear solver object          */
    /* mirror of SUNNonlinSolGetCurIter for the linear-solver interface
       (the NLS object is detached from CVodeMem during its solve) */
    pub cv_nls_curiter: i32,
    pub ownNLS: bool,                 /* flag indicating NLS ownership    */
    pub nls_f: Option<CVRhsFn>, /* f(t,y(t)) used in the nonlinear solver */
    pub convfail: i32,          /* flag: Jacobian update may be needed    */

    /*------------------
      Linear Solver Data
      ------------------*/
    /* In C: cv_linit/cv_lsetup/cv_lsolve/cv_lfree function pointers +
       void* cv_lmem. Here: dispatching enum (see LsModule). */
    pub cv_lmem: LsModule,
    pub cv_msbp: i64,          /* max number of steps between lsetup calls */
    pub cv_dgmax_lsetup: f64,  /* gamma ratio threshold for lsetup         */

    /*------------
      Saved Values
      ------------*/
    pub cv_qu: i32,         /* last successful q value used                */
    pub cv_nstlp: i64,      /* step number of last setup call              */
    pub cv_h0u: f64,        /* actual initial stepsize                     */
    pub cv_hu: f64,         /* last successful h value used                */
    pub cv_saved_tq5: f64,  /* saved value of tq[5]                        */
    pub cv_jcur: bool,      /* is Jacobian info for linear solver current? */
    pub cv_tolsf: f64,      /* tolerance scale factor                      */
    pub cv_qmax_alloc: i32, /* value of qmax used when allocating memory   */
    pub cv_indx_acor: i32,  /* index of the zn vector with saved acor      */

    /*--------------------------------------------------------------------
      Flags turned ON by CVodeInit and read by CVodeReInit
      --------------------------------------------------------------------*/
    pub cv_VabstolMallocDone: bool,
    pub cv_MallocDone: bool,

    /*-------------------------
      Stability Limit Detection
      -------------------------*/
    pub cv_sldeton: bool,        /* is Stability Limit Detection on?   */
    pub cv_ssdat: [[f64; 4]; 6], /* scaled data array for STALD        */
    pub cv_nscon: i32,           /* counter for STALD method           */
    pub cv_nor: i64,             /* counter for number of order reductions */

    /*----------------
      Rootfinding Data
      ----------------*/
    pub cv_gfun: Option<CVRootFn>, /* function g for roots sought          */
    pub cv_nrtfn: i32,             /* number of components of g            */
    pub cv_iroots: Vec<i32>,       /* array for root information           */
    pub cv_rootdir: Vec<i32>,      /* array specifying direction of zero-crossing */
    pub cv_tlo: f64,               /* nearest endpoint of interval in root search */
    pub cv_thi: f64,               /* farthest endpoint of interval in root search */
    pub cv_trout: f64,             /* t value returned by rootfinding routine */
    pub cv_glo: Vec<f64>,          /* saved array of g values at t = tlo   */
    pub cv_ghi: Vec<f64>,          /* saved array of g values at t = thi   */
    pub cv_grout: Vec<f64>,        /* array of g values at t = trout       */
    pub cv_ttol: f64,              /* tolerance on root location trout     */
    pub cv_taskc: i32,             /* copy of parameter itask              */
    pub cv_irfnd: i32,             /* flag showing whether last step had a root */
    pub cv_nge: i64,               /* counter for g evaluations            */
    pub cv_gactive: Vec<bool>,     /* array with active/inactive event functions */
    pub cv_mxgnull: i32,           /* num. warning messages about possible g==0 */

    /*---------------------------
      Inequality Constraints Data
      ---------------------------*/
    pub cv_constraints: NVector,   /* vector of constraint flags        */
    pub cv_constraintsSet: bool,   /* constraints vector present        */
    pub constraint_corrections: i64, /* total constraint corrections    */
    pub constraint_fails: i64,     /* total constraint failures         */
    pub max_constraint_fails: i32, /* max failures allowed in a step    */

    /*---------------
      Projection Data
      ---------------*/
    pub proj_mem: Option<Box<CVodeProjMem>>, /* projection memory structure  */
    pub proj_enabled: bool,  /* flag indicating if projection is enabled  */
    pub proj_applied: bool,  /* flag indicating if projection was applied */
    pub proj_p: [f64; L_MAX], /* coefficients of p(x) (degree q poly)     */

    /*-----------------------
      Fused Vector Operations
      -----------------------*/
    pub cv_usefused: bool, /* fused CVODE kernels (never in pure-Rust build) */

    /*----------------
      Resizing History
      ----------------*/
    pub first_step_after_resize: bool, /* flag to signal a resize happened */
}

/* ===============================================================
   Error handler: the C cvProcessError routes printf-style messages
   to the SUNContext error handler stack; here messages go to stderr
   (equivalent to the default SUNLogErrHandlerFn behavior).
   =============================================================== */
pub fn cvProcessError(
    _cv_mem: Option<&CVodeMem>,
    error_code: i32,
    line: u32,
    func: &str,
    file: &str,
    msg: &str,
) {
    if error_code == CV_WARNING {
        eprintln!("\n[CVODE WARNING] {file}:{line} in {func}\n  {msg}\n");
    } else {
        eprintln!("\n[CVODE ERROR] {file}:{line} in {func}\n  {msg}\n");
    }
}

/* Error messages (cvode_impl.h) */
pub const MSGCV_NO_MEM: &str = "cvode_mem = NULL illegal.";
pub const MSGCV_CVMEM_FAIL: &str = "Allocation of cvode_mem failed.";
pub const MSGCV_MEM_FAIL: &str = "A memory request failed.";
pub const MSGCV_BAD_LMM: &str =
    "Illegal value for lmm. The legal values are CV_ADAMS and CV_BDF.";
pub const MSGCV_NULL_SUNCTX: &str = "sunctx = NULL illegal.";
pub const MSGCV_NO_MALLOC: &str = "Attempt to call before CVodeInit.";
pub const MSGCV_NEG_MAXORD: &str = "maxord <= 0 illegal.";
pub const MSGCV_BAD_MAXORD: &str = "Illegal attempt to increase maximum method order.";
pub const MSGCV_SET_SLDET: &str =
    "Attempt to use stability limit detection with the CV_ADAMS method illegal.";
pub const MSGCV_NEG_HMIN: &str = "hmin < 0 illegal.";
pub const MSGCV_NEG_HMAX: &str = "hmax < 0 illegal.";
pub const MSGCV_BAD_HMIN_HMAX: &str = "Inconsistent step size limits: hmin > hmax.";
pub const MSGCV_BAD_RELTOL: &str = "reltol < 0 illegal.";
pub const MSGCV_BAD_ABSTOL: &str = "abstol has negative component(s) (illegal).";
pub const MSGCV_NULL_ABSTOL: &str = "abstol = NULL illegal.";
pub const MSGCV_NULL_Y0: &str = "y0 = NULL illegal.";
pub const MSGCV_Y0_FAIL_CONSTR: &str = "y0 fails to satisfy constraints.";
pub const MSGCV_NULL_F: &str = "f = NULL illegal.";
pub const MSGCV_NULL_G: &str = "g = NULL illegal.";
pub const MSGCV_BAD_NVECTOR: &str = "A required vector operation is not implemented.";
pub const MSGCV_BAD_CONSTR: &str = "Illegal values in constraints vector.";
pub const MSGCV_BAD_K: &str = "Illegal value for k.";
pub const MSGCV_NULL_DKY: &str = "dky = NULL illegal.";
pub const MSGCV_NO_ROOT: &str = "Rootfinding was not initialized.";
pub const MSGCV_NLS_INIT_FAIL: &str = "The nonlinear solver's init routine failed.";
pub const MSGCV_NO_TOL: &str = "No integration tolerances have been specified.";
pub const MSGCV_LSOLVE_NULL: &str = "The linear solver's solve routine is NULL.";
pub const MSGCV_YOUT_NULL: &str = "yout = NULL illegal.";
pub const MSGCV_TRET_NULL: &str = "tret = NULL illegal.";
pub const MSGCV_BAD_EWT: &str = "Initial ewt has component(s) equal to zero (illegal).";
pub const MSGCV_BAD_ITASK: &str = "Illegal value for itask.";
pub const MSGCV_BAD_H0: &str = "h0 and tout - t0 inconsistent.";
pub const MSGCV_EWT_FAIL: &str = "The user-provide EwtSet function failed.";
pub const MSGCV_LINIT_FAIL: &str = "The linear solver's init routine failed.";
pub const MSGCV_HNIL_DONE: &str =
    "The above warning has been issued mxhnil times and will not be issued again for this problem.";
pub const MSGCV_TOO_CLOSE: &str = "tout too close to t0 to start integration.";
pub const MSGCV_RHSFUNC_FIRST: &str = "The right-hand side routine failed at the first call.";
pub const MSGCV_INACTIVE_ROOTS: &str =
    "At the end of the first step, there are still some root functions identically 0. \
     This warning will not be issued again.";
