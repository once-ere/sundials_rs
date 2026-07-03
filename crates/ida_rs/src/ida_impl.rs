/* -----------------------------------------------------------------
 * Translated from src/ida/ida_impl.h and the constants / typedefs
 * of include/ida/ida.h (IDA 7.7.0).
 * Main integrator memory block (struct IDAMemRec) and internal
 * constants for the DAE solver F(t, y, y') = 0 (BDF + Newton).
 * Field names keep the C names (ida_ prefix) so the translation of
 * ida.c reads line-for-line.
 * -----------------------------------------------------------------*/
use crate::ida_ls_impl::IDALsMem;
use crate::nvector_serial::NVector;
use crate::sundials_context::SUNContext;
use crate::sundials_math::SUNRpowerR;
use crate::sundials_nonlinearsolver::NonlinearSolver;
use crate::sundials_types::{UserData, SUNFALSE, SUNTRUE, SUN_UNIT_ROUNDOFF};

/* ===============================================================
   User-supplied function types (ida.h)
   =============================================================== */

pub type IDAResFn =
    fn(t: f64, yy: &NVector, yp: &NVector, rr: &mut NVector, user_data: &mut UserData) -> i32;

pub type IDARootFn =
    fn(t: f64, y: &NVector, yp: &NVector, gout: &mut [f64], user_data: &mut UserData) -> i32;

pub type IDAEwtFn = fn(y: &NVector, ewt: &mut NVector, user_data: &mut UserData) -> i32;

/* ===============================================================
   IDA constants (ida.h)
   =============================================================== */

/* itask */
pub const IDA_NORMAL: i32 = 1;
pub const IDA_ONE_STEP: i32 = 2;

/* icopt */
pub const IDA_YA_YDP_INIT: i32 = 1;
pub const IDA_Y_INIT: i32 = 2;

/* return values */
pub const IDA_SUCCESS: i32 = 0;
pub const IDA_TSTOP_RETURN: i32 = 1;
pub const IDA_ROOT_RETURN: i32 = 2;

pub const IDA_WARNING: i32 = 99;

pub const IDA_TOO_MUCH_WORK: i32 = -1;
pub const IDA_TOO_MUCH_ACC: i32 = -2;
pub const IDA_ERR_FAIL: i32 = -3;
pub const IDA_CONV_FAIL: i32 = -4;

pub const IDA_LINIT_FAIL: i32 = -5;
pub const IDA_LSETUP_FAIL: i32 = -6;
pub const IDA_LSOLVE_FAIL: i32 = -7;
pub const IDA_RES_FAIL: i32 = -8;
pub const IDA_REP_RES_ERR: i32 = -9;
pub const IDA_RTFUNC_FAIL: i32 = -10;
pub const IDA_CONSTR_FAIL: i32 = -11;

pub const IDA_FIRST_RES_FAIL: i32 = -12;
pub const IDA_LINESEARCH_FAIL: i32 = -13;
pub const IDA_NO_RECOVERY: i32 = -14;
pub const IDA_NLS_INIT_FAIL: i32 = -15;
pub const IDA_NLS_SETUP_FAIL: i32 = -16;
pub const IDA_NLS_FAIL: i32 = -17;

pub const IDA_MEM_NULL: i32 = -20;
pub const IDA_MEM_FAIL: i32 = -21;
pub const IDA_ILL_INPUT: i32 = -22;
pub const IDA_NO_MALLOC: i32 = -23;
pub const IDA_BAD_EWT: i32 = -24;
pub const IDA_BAD_K: i32 = -25;
pub const IDA_BAD_T: i32 = -26;
pub const IDA_BAD_DKY: i32 = -27;
pub const IDA_VECTOROP_ERR: i32 = -28;

pub const IDA_CONTEXT_ERR: i32 = -29;

pub const IDA_TOO_CLOSE: i32 = -60;

pub const IDA_UNRECOGNIZED_ERROR: i32 = -99;

/* itol (internal; ida.c) */
pub const IDA_NN: i32 = 0;
pub const IDA_SS: i32 = 1;
pub const IDA_SV: i32 = 2;
pub const IDA_WF: i32 = 3;

/* ===============================================================
   Basic IDA constants (ida_impl.h)
   =============================================================== */

pub const HMAX_INV_DEFAULT: f64 = 0.0; /* hmax_inv default value          */
pub const HMIN_DEFAULT: f64 = 0.0; /* hmin default value              */
pub const MAXORD_DEFAULT: i32 = 5; /* maxord default value            */
pub const MXORDP1: usize = 6; /* max. number of N_Vectors in phi */
pub const MXSTEP_DEFAULT: i64 = 500; /* mxstep default value            */

pub const ETA_MAX_FX_DEFAULT: f64 = 2.0; /* threshold to increase step size   */
pub const ETA_MIN_FX_DEFAULT: f64 = 1.0; /* threshold to decrease step size   */
pub const ETA_MAX_DEFAULT: f64 = 2.0; /* max step size increase factor     */
pub const ETA_MIN_DEFAULT: f64 = 0.5; /* min step size decrease factor     */
pub const ETA_LOW_DEFAULT: f64 = 0.9; /* upper bound on decrease factor    */
pub const ETA_MIN_EF_DEFAULT: f64 = 0.25; /* err test fail min decrease factor */
pub const ETA_CF_DEFAULT: f64 = 0.25; /* NLS failure decrease factor       */

pub const DCJ_DEFAULT: f64 = 0.25; /* constant for updating Jacobian/preconditioner */

pub const MAX_CONSTRAINT_FAILS: i32 = 10;

/* Return values for lower level routines used by IDASolve and functions
   provided to the nonlinear solver */

pub const IDA_RES_RECVR: i32 = 1;
pub const IDA_LSETUP_RECVR: i32 = 2;
pub const IDA_LSOLVE_RECVR: i32 = 3;
pub const IDA_NLS_SETUP_RECVR: i32 = 4;

/* ===============================================================
   Default-value constants referenced by IDACreate (defined at the
   top of src/ida/ida.c; hoisted here because the Default impl below
   is the translation of IDACreate's default-setting block).
   =============================================================== */

pub const MXNCF: i32 = 10; /* max number of convergence failures allowed */
pub const MXNEF: i32 = 10; /* max number of error test failures allowed  */
pub const MAXNH: i32 = 5; /* max. number of h tries in IC calc. */
pub const MAXNJ: i32 = 4; /* max. number of J tries in IC calc. */
pub const MAXNI: i32 = 10; /* max. Newton iterations in IC calc. */
pub const EPCON: f64 = 0.33; /* Newton convergence test constant */
pub const MAXBACKS: i32 = 100; /* max backtracks per Newton step in IDACalcIC */

/* Shortcuts IDA_PROFILER / IDA_LOGGER are not ported: the workspace
   SUNContext is a unit-like struct because the default SUNDIALS
   build compiles solver profiling/logging macros out. */

/* ===============================================================
   Linear solver module attached to IDA.
   In C this is the (ida_linit, ida_lsetup, ida_lsolve, ida_lperf,
   ida_lfree) function-pointer five-tuple plus the void* ida_lmem;
   here it is an enum dispatched in ida.rs (donor cvode_impl.rs
   pattern).  IDA has a single interface module (IDALS, ida_ls.rs).
   The module is take()n out of IDAMem for the duration of a call so
   its routines may borrow the solver memory mutably.
   The C `ida_lperf != NULL` guards (ida_lperf is only installed for
   iterative SUNLinearSolver objects) dispatch on IDALsMem.iterative;
   the `ida_lsetup != NULL` guards dispatch on IDALsMem.setup_disabled
   (see ida_ls_impl.rs).
   =============================================================== */

#[derive(Default)]
pub enum LsModule {
    #[default]
    None,
    /// ida_ls.rs interface (IDASetLinearSolver)
    Ls(Box<IDALsMem>),
}

impl LsModule {
    pub fn is_none(&self) -> bool {
        matches!(self, LsModule::None)
    }
}

/* ===============================================================
   Main integrator memory block (struct IDAMemRec)
   =============================================================== */

pub struct IDAMem {
    pub ida_sunctx: SUNContext,

    pub ida_uround: f64, /* machine unit roundoff */

    /*--------------------------
      Problem Specification Data
      --------------------------*/
    pub ida_res: Option<IDAResFn>, /* F(t,y(t),y'(t))=0; the function F     */
    pub ida_user_data: UserData,   /* user pointer passed to res            */

    pub ida_itol: i32,       /* itol = IDA_SS, IDA_SV, IDA_WF, IDA_NN */
    pub ida_rtol: f64,       /* relative tolerance                    */
    pub ida_Satol: f64,      /* scalar absolute tolerance             */
    pub ida_Vatol: NVector,  /* vector absolute tolerance             */
    pub ida_atolmin0: bool,  /* flag indicating that min(atol) = 0    */
    pub ida_user_efun: bool, /* SUNTRUE if user provides efun         */
    pub ida_efun: Option<IDAEwtFn>, /* function to set ewt            */
    /* (C ida_edata — the user pointer passed to efun — is not stored:
       a user efun receives ida_user_data, and the internal IDAEwtSet
       operates on &mut IDAMem directly, as in the donor's drop of
       cv_e_data.) */
    pub ida_suppressalg: bool, /* SUNTRUE means suppress algebraic vars
                               in local error tests                  */

    /*-----------------------------------------------
      Divided differences array and associated arrays
      -----------------------------------------------*/
    pub ida_phi: Vec<NVector>, /* phi = (maxord+1) arrays of divided differences
                               (C: N_Vector ida_phi[MXORDP1]; Vec per the
                               donor cv_zn convention)                        */

    pub ida_psi: [f64; MXORDP1], /* differences in t (sums of recent step sizes)   */
    pub ida_alpha: [f64; MXORDP1], /* ratios of current stepsize to psi values       */
    pub ida_beta: [f64; MXORDP1], /* ratios of current to previous product of psi's */
    pub ida_sigma: [f64; MXORDP1], /* product successive alpha values and factorial  */
    pub ida_gamma: [f64; MXORDP1], /* sum of reciprocals of psi values               */

    /*-------------------------
      N_Vectors for integration
      -------------------------*/
    pub ida_ewt: NVector,       /* error weight vector                            */
    pub ida_yy: NVector,        /* work space for y vector (= user's yret)        */
    pub ida_yp: NVector,        /* work space for y' vector (= user's ypret)      */
    pub ida_yypredict: NVector, /* predicted y vector                             */
    pub ida_yppredict: NVector, /* predicted y' vector                            */
    pub ida_delta: NVector,     /* residual vector                                */
    pub ida_id: NVector,        /* bit vector for diff./algebraic components      */
    pub ida_savres: NVector,    /* saved residual vector                          */
    pub ida_ee: NVector,        /* accumulated corrections to y vector, but
                                set equal to estimated local errors upon
                                successful return                              */
    pub ida_tempv1: NVector,    /* work space vector                              */
    pub ida_tempv2: NVector,    /* work space vector                              */
    pub ida_tempv3: NVector,    /* work space vector                              */
    /* (C ida_ynew / ida_ypnew / ida_delnew / ida_dtemp — the IDACalcIC
       work "vectors" — are not stored: in C they are aliases installed
       by IDACalcIC (ynew = tempv2, ypnew = ee, delnew = phi[2],
       dtemp = phi[3]); ida_ic.rs uses the underlying fields directly,
       as in the donor's drop of the cvodes senswrapper aliases.)       */

    /*------------------------------
      Variables for use by IDACalcIC
      ------------------------------*/
    pub ida_t0: f64,      /* initial t                                      */
    pub ida_yy0: NVector, /* initial y vector (user-supplied).              */
    pub ida_yp0: NVector, /* initial y' vector (user-supplied).             */

    pub ida_icopt: i32,    /* IC calculation user option                     */
    pub ida_lsoff: bool,   /* IC calculation linesearch turnoff option       */
    pub ida_maxnh: i32,    /* max. number of h tries in IC calculation       */
    pub ida_maxnj: i32,    /* max. number of J tries in IC calculation       */
    pub ida_maxnit: i32,   /* max. number of Netwon iterations in IC calc.   */
    pub ida_nbacktr: i32,  /* number of IC linesearch backtrack operations   */
    pub ida_sysindex: i32, /* computed system index (0 or 1)                 */
    pub ida_maxbacks: i32, /* max backtracks per Newton step                 */
    pub ida_epiccon: f64,  /* IC nonlinear convergence test constant         */
    pub ida_steptol: f64,  /* minimum Newton step size in IC calculation     */
    pub ida_tscale: f64,   /* time scale factor = abs(tout1 - t0)            */

    /* Tstop information */
    pub ida_tstopset: bool,
    pub ida_tstop: f64,

    /* Step Data */
    pub ida_kk: i32,    /* current BDF method order                              */
    pub ida_kused: i32, /* method order used on last successful step             */
    pub ida_knew: i32,  /* order for next step from order decrease decision      */
    pub ida_phase: i32, /* flag to trigger step doubling in first few steps      */
    pub ida_ns: i32,    /* counts steps at fixed stepsize and order              */

    pub ida_hin: f64,      /* initial step                                      */
    pub ida_h0u: f64,      /* actual initial stepsize                           */
    pub ida_hh: f64,       /* current step size h                               */
    pub ida_hused: f64,    /* step size used on last successful step            */
    pub ida_eta: f64,      /* eta = hnext / hused                               */
    pub ida_tn: f64,       /* current internal value of t                       */
    pub ida_tretlast: f64, /* value of tret previously returned by IDASolve     */
    pub ida_cj: f64,       /* current value of scalar (-alphas/hh) in Jacobian  */
    pub ida_cjlast: f64,   /* cj value saved from last successful step          */
    pub ida_cjold: f64,    /* cj value saved from last call to lsetup           */
    pub ida_cjratio: f64,  /* ratio of cj values: cj/cjold                      */
    pub ida_ss: f64,       /* scalar used in Newton iteration convergence test  */
    pub ida_oldnrm: f64,   /* norm of previous nonlinear solver update          */
    pub ida_epsNewt: f64,  /* test constant in Newton convergence test          */
    pub ida_epcon: f64,    /* coefficient of the Newton convergence test        */
    pub ida_toldel: f64,   /* tolerance in direct test on Newton corrections    */

    /*------
      Limits
      ------*/
    pub ida_maxncf: i32, /* max number of convergence failures                */
    pub ida_maxnef: i32, /* max number of error test failures                 */

    pub ida_maxord: i32,       /* max value of method order k:                      */
    pub ida_maxord_alloc: i32, /* value of maxord used when allocating memory       */
    pub ida_mxstep: i64,       /* max number of internal steps for one user call    */
    pub ida_hmax_inv: f64,     /* inverse of max. step size hmax (default = 0.0)    */
    pub ida_hmin: f64,         /* min step size hmin (default = 0.0)                */

    pub ida_eta_max_fx: f64, /* threshold to increase step size */
    pub ida_eta_min_fx: f64, /* threshold to decrease step size */
    pub ida_eta_max: f64,    /* max step size increase factor   */
    pub ida_eta_min: f64,    /* min step size decrease factor   */
    pub ida_eta_low: f64,    /* upper bound on decrease factor  */
    pub ida_eta_min_ef: f64, /* eta >= eta_min_ef after an error test failure */
    pub ida_eta_cf: f64,     /* eta on a nonlinear solver convergence failure */

    /*--------
      Counters
      --------*/
    pub ida_nst: i64,     /* number of internal steps taken                    */
    pub ida_nre: i64,     /* number of function (res) calls                    */
    pub ida_ncfn: i64,    /* number of corrector convergence failures          */
    pub ida_netf: i64,    /* number of error test failures                     */
    pub ida_nni: i64,     /* number of Newton iterations performed             */
    pub ida_nnf: i64,     /* number of Newton convergence failures             */
    pub ida_nsetups: i64, /* number of lsetup calls                            */

    /*------------------
      Space requirements
      ------------------*/
    pub ida_lrw1: i64, /* no. of sunrealtype words in 1 N_Vector            */
    pub ida_liw1: i64, /* no. of integer words in 1 N_Vector                */
    pub ida_lrw: i64,  /* number of sunrealtype words in IDA work vectors   */
    pub ida_liw: i64,  /* no. of integer words in IDA work vectors          */

    pub ida_tolsf: f64, /* tolerance scale factor (saved value)              */

    /* Flags to verify correct calling sequence */
    pub ida_SetupDone: bool, /* set to SUNFALSE by IDAMalloc and IDAReInit
                             set to SUNTRUE by IDACalcIC or IDASolve      */

    pub ida_VatolMallocDone: bool,
    pub ida_idMallocDone: bool,

    pub ida_MallocDone: bool, /* set to SUNFALSE by IDACreate
                              set to SUNTRUE by IDAMAlloc
                              tested by IDAReInit and IDASolve             */

    /*---------------------
      Nonlinear Solver Data
      ---------------------*/
    pub NLS: Option<NonlinearSolver>, /* nonlinear solver object       */
    pub ownNLS: bool,                 /* flag indicating NLS ownership */
    pub nls_res: Option<IDAResFn>,    /* F(t,y(t),y'(t))=0; used in the nonlinear
                                      solver */

    /*------------------
      Linear Solver Data
      ------------------*/
    /* In C: ida_linit/ida_lsetup/ida_lsolve/ida_lperf/ida_lfree
       function pointers + void* ida_lmem. Here: dispatching enum
       (see LsModule). */
    pub ida_lmem: LsModule, /* linear solver interface structure */
    pub ida_dcj: f64,       /* parameter that determines cj ratio thresholds for calling
                             * the linear solver setup function */

    /* Flag to indicate successful ida_linit call */
    pub ida_linitOK: bool,

    /*----------------
      Rootfinding Data
      ----------------*/
    pub ida_gfun: Option<IDARootFn>, /* Function g for roots sought                     */
    pub ida_nrtfn: i32,              /* number of components of g                       */
    pub ida_iroots: Vec<i32>,        /* array for root information                      */
    pub ida_rootdir: Vec<i32>,       /* array specifying direction of zero-crossing     */
    pub ida_tlo: f64,                /* nearest endpoint of interval in root search     */
    pub ida_thi: f64,                /* farthest endpoint of interval in root search    */
    pub ida_trout: f64,              /* t return value from rootfinder routine          */
    pub ida_glo: Vec<f64>,           /* saved array of g values at t = tlo              */
    pub ida_ghi: Vec<f64>,           /* saved array of g values at t = thi              */
    pub ida_grout: Vec<f64>,         /* array of g values at t = trout                  */
    pub ida_ttol: f64,               /* tolerance on root location                      */
    pub ida_irfnd: i32,              /* flag showing whether last step had a root       */
    pub ida_nge: i64,                /* counter for g evaluations                       */
    pub ida_gactive: Vec<bool>,      /* array with active/inactive event functions      */
    pub ida_mxgnull: i32,            /* number of warning messages about possible g==0  */

    /*---------------------------
      Inequality Constraints Data
      ---------------------------*/
    pub ida_constraints: NVector, /* vector of inequality constraint flags */
    pub ida_constraintsSet: bool, /* (adaptation) in C the "no constraints"
                                  state is ida_constraints == NULL; an owned
                                  NVector cannot be NULL, so this flag
                                  carries that state (donor CVodeMem naming) */
    pub constraint_corrections: i64, /* total constraint corrections   */
    pub constraint_fails: i64,    /* total constraint failures             */
    pub max_constraint_fails: i32, /* max failures allowed in a step        */

    /* (C ida_cvals / ida_dvals / ida_Xvecs / ida_Zvecs — the
       scalar/vector alias arrays feeding the fused
       N_VLinearCombination / N_VScaleAddMulti kernels — are not
       stored: the serial fused kernels are reproduced inline, as in
       the donor's drop of cv_cvals/cv_Xvecs.) */
}

/* ===============================================================
   Default construction: the C IDACreate memset(0)s IDAMemRec and
   then sets the default optional inputs (ida.c, IDACreate).
   IDACreate in ida.rs builds on this Default exactly like the
   donor's CVodeCreate builds its struct literal.
   =============================================================== */
impl Default for IDAMem {
    fn default() -> Self {
        let uround = SUN_UNIT_ROUNDOFF;
        IDAMem {
            ida_sunctx: SUNContext::default(),

            /* Set unit roundoff in IDA_mem */
            ida_uround: uround,

            /* Set default values for integrator optional inputs */
            ida_res: None,
            ida_user_data: None,
            ida_itol: IDA_NN,
            ida_rtol: 0.0,
            ida_Satol: 0.0,
            ida_Vatol: NVector::default(),
            ida_atolmin0: SUNTRUE,
            ida_user_efun: SUNFALSE,
            ida_efun: None,
            ida_maxord: MAXORD_DEFAULT,
            ida_mxstep: MXSTEP_DEFAULT,
            ida_hmax_inv: HMAX_INV_DEFAULT,
            ida_hmin: HMIN_DEFAULT,
            ida_eta_max_fx: ETA_MAX_FX_DEFAULT,
            ida_eta_min_fx: ETA_MIN_FX_DEFAULT,
            ida_eta_max: ETA_MAX_DEFAULT,
            ida_eta_low: ETA_LOW_DEFAULT,
            ida_eta_min: ETA_MIN_DEFAULT,
            ida_eta_min_ef: ETA_MIN_EF_DEFAULT,
            ida_eta_cf: ETA_CF_DEFAULT,
            ida_hin: 0.0,
            ida_epcon: EPCON,
            ida_maxnef: MXNEF,
            ida_maxncf: MXNCF,
            ida_suppressalg: SUNFALSE,
            ida_id: NVector::default(),
            ida_tstopset: SUNFALSE,
            ida_dcj: DCJ_DEFAULT,

            /* Initialize inequality constraint variables */
            ida_constraints: NVector::default(),
            ida_constraintsSet: SUNFALSE,
            constraint_corrections: 0,
            constraint_fails: 0,
            max_constraint_fails: MAX_CONSTRAINT_FAILS,

            /* set the saved value maxord_alloc */
            ida_maxord_alloc: MAXORD_DEFAULT,

            /* Set default values for IC optional inputs
               (ida.c: TWOTHIRDS = SUN_RCONST(0.667)) */
            ida_epiccon: 0.01 * EPCON,
            ida_maxnh: MAXNH,
            ida_maxnj: MAXNJ,
            ida_maxnit: MAXNI,
            ida_maxbacks: MAXBACKS,
            ida_lsoff: SUNFALSE,
            ida_steptol: SUNRpowerR(uround, 0.667),

            /* Initialize lrw and liw */
            ida_lrw: 25 + 5 * MXORDP1 as i64,
            ida_liw: 38,

            /* No mallocs have been done yet */
            ida_VatolMallocDone: SUNFALSE,
            ida_idMallocDone: SUNFALSE,
            ida_MallocDone: SUNFALSE,

            /* Initialize nonlinear solver variables */
            NLS: None,
            ownNLS: SUNFALSE,

            /* All remaining fields: memset(0) in IDACreate */
            nls_res: None,

            ida_phi: Vec::new(),
            ida_psi: [0.0; MXORDP1],
            ida_alpha: [0.0; MXORDP1],
            ida_beta: [0.0; MXORDP1],
            ida_sigma: [0.0; MXORDP1],
            ida_gamma: [0.0; MXORDP1],

            ida_ewt: NVector::default(),
            ida_yy: NVector::default(),
            ida_yp: NVector::default(),
            ida_yypredict: NVector::default(),
            ida_yppredict: NVector::default(),
            ida_delta: NVector::default(),
            ida_savres: NVector::default(),
            ida_ee: NVector::default(),
            ida_tempv1: NVector::default(),
            ida_tempv2: NVector::default(),
            ida_tempv3: NVector::default(),

            ida_t0: 0.0,
            ida_yy0: NVector::default(),
            ida_yp0: NVector::default(),
            ida_icopt: 0,
            ida_nbacktr: 0,
            ida_sysindex: 0,
            ida_tscale: 0.0,

            ida_tstop: 0.0,

            ida_kk: 0,
            ida_kused: 0,
            ida_knew: 0,
            ida_phase: 0,
            ida_ns: 0,

            ida_h0u: 0.0,
            ida_hh: 0.0,
            ida_hused: 0.0,
            ida_eta: 0.0,
            ida_tn: 0.0,
            ida_tretlast: 0.0,
            ida_cj: 0.0,
            ida_cjlast: 0.0,
            ida_cjold: 0.0,
            ida_cjratio: 0.0,
            ida_ss: 0.0,
            ida_oldnrm: 0.0,
            ida_epsNewt: 0.0,
            ida_toldel: 0.0,

            ida_nst: 0,
            ida_nre: 0,
            ida_ncfn: 0,
            ida_netf: 0,
            ida_nni: 0,
            ida_nnf: 0,
            ida_nsetups: 0,

            ida_lrw1: 0,
            ida_liw1: 0,

            ida_tolsf: 0.0,

            ida_SetupDone: SUNFALSE,

            ida_lmem: LsModule::None,
            ida_linitOK: SUNFALSE,

            ida_gfun: None,
            ida_nrtfn: 0,
            ida_iroots: Vec::new(),
            ida_rootdir: Vec::new(),
            ida_tlo: 0.0,
            ida_thi: 0.0,
            ida_trout: 0.0,
            ida_glo: Vec::new(),
            ida_ghi: Vec::new(),
            ida_grout: Vec::new(),
            ida_ttol: 0.0,
            ida_irfnd: 0,
            ida_nge: 0,
            ida_gactive: Vec::new(),
            ida_mxgnull: 0,
        }
    }
}

/* ===============================================================
   IDA internal functions (prototypes in ida_impl.h; the
   implementations live in the ida module files):

     IDAEwtSet(ycur, weight, &mut IDAMem)   — internal ewtSet function
                                              (ida.rs)
     IDAErrHandler                          — internal errHandler
                                              function (subsumed by
                                              IDAProcessError below)
     IDAWrmsNorm(IDA_mem, x, w, mask)       — norm function, also used
                                              for IC, so it is global
                                              (ida.rs)
     idaNlsInit(IDA_mem)                    — nonlinear solver
                                              initialization (ida_nls.rs)
   =============================================================== */

/* ===============================================================
   Error handler: the C IDAProcessError routes printf-style messages
   to the SUNContext error handler stack; here messages go to stderr
   (equivalent to the default SUNLogErrHandlerFn behavior).
   =============================================================== */
pub fn IDAProcessError(
    _ida_mem: Option<&IDAMem>,
    error_code: i32,
    line: u32,
    func: &str,
    file: &str,
    msg: &str,
) {
    if error_code == IDA_WARNING {
        eprintln!("\n[IDA WARNING] {file}:{line} in {func}\n  {msg}\n");
    } else {
        eprintln!("\n[IDA ERROR] {file}:{line} in {func}\n  {msg}\n");
    }
}

/* ===============================================================
   IDA error messages (ida_impl.h); C printf formats kept verbatim
   with SUN_FORMAT_G expanded to its double-precision form "%.15g"
   (call sites render the values via sundials_utils::fmt_g).
   =============================================================== */

pub const MSG_TIME: &str = "t = %.15g";
pub const MSG_TIME_H: &str = "t = %.15g and h = %.15g";
pub const MSG_TIME_INT: &str =
    "t = %.15g is not between tcur - hold = %.15g and tcur = %.15g";
pub const MSG_TIME_TOUT: &str = "tout = %.15g";
pub const MSG_TIME_TSTOP: &str = "tstop = %.15g";

/* General errors */

pub const MSG_MEM_FAIL: &str = "A memory request failed.";
pub const MSG_NULL_SUNCTX: &str = "sunctx = NULL illegal.";
pub const MSG_NO_MEM: &str = "ida_mem = NULL illegal.";
pub const MSG_NO_MALLOC: &str = "Attempt to call before IDAMalloc.";
pub const MSG_BAD_NVECTOR: &str = "A required vector operation is not implemented.";

/* Initialization errors */

pub const MSG_Y0_NULL: &str = "y0 = NULL illegal.";
pub const MSG_YP0_NULL: &str = "yp0 = NULL illegal.";
pub const MSG_BAD_ITOL: &str =
    "Illegal value for itol. The legal values are IDA_SS, IDA_SV, and IDA_WF.";
pub const MSG_RES_NULL: &str = "res = NULL illegal.";
pub const MSG_BAD_RTOL: &str = "rtol < 0 illegal.";
pub const MSG_ATOL_NULL: &str = "atol = NULL illegal.";
pub const MSG_BAD_ATOL: &str = "Some atol component < 0.0 illegal.";
pub const MSG_ROOT_FUNC_NULL: &str = "g = NULL illegal.";

pub const MSG_MISSING_ID: &str = "id = NULL but suppressalg option on.";
pub const MSG_NO_TOLS: &str = "No integration tolerances have been specified.";
pub const MSG_FAIL_EWT: &str = "The user-provide EwtSet function failed.";
pub const MSG_BAD_EWT: &str = "Some initial ewt component = 0.0 illegal.";
pub const MSG_Y0_FAIL_CONSTR: &str = "y0 fails to satisfy constraints.";
pub const MSG_LSOLVE_NULL: &str = "The linear solver's solve routine is NULL.";
pub const MSG_LINIT_FAIL: &str = "The linear solver's init routine failed.";
pub const MSG_NLS_INIT_FAIL: &str = "The nonlinear solver's init routine failed.";

/* IDACalcIC error messages */

pub const MSG_IC_BAD_ICOPT: &str = "icopt has an illegal value.";
pub const MSG_IC_BAD_MAXBACKS: &str = "maxbacks <= 0 illegal.";
pub const MSG_IC_MISSING_ID: &str = "id = NULL conflicts with icopt.";
pub const MSG_IC_TOO_CLOSE: &str =
    "tout1 too close to t0 to attempt initial condition calculation.";
pub const MSG_IC_BAD_ID: &str = "id has illegal values.";
pub const MSG_IC_BAD_EWT: &str = "Some initial ewt component = 0.0 illegal.";
pub const MSG_IC_RES_NONREC: &str = "The residual function failed unrecoverably. ";
pub const MSG_IC_RES_FAIL: &str = "The residual function failed at the first call. ";
pub const MSG_IC_SETUP_FAIL: &str = "The linear solver setup failed unrecoverably.";
pub const MSG_IC_SOLVE_FAIL: &str = "The linear solver solve failed unrecoverably.";
pub const MSG_IC_NO_RECOVERY: &str =
    "The residual routine or the linear setup or solve routine had a \
     recoverable error, but IDACalcIC was unable to recover.";
pub const MSG_IC_FAIL_CONSTR: &str = "Unable to satisfy the inequality constraints.";
pub const MSG_IC_FAILED_LINS: &str =
    "The linesearch algorithm failed: step too small or too many backtracks.";
pub const MSG_IC_CONV_FAILED: &str = "Newton/Linesearch algorithm failed to converge.";

/* IDASolve error messages */

pub const MSG_YRET_NULL: &str = "yret = NULL illegal.";
pub const MSG_YPRET_NULL: &str = "ypret = NULL illegal.";
pub const MSG_TRET_NULL: &str = "tret = NULL illegal.";
pub const MSG_BAD_ITASK: &str = "itask has an illegal value.";
pub const MSG_TOO_CLOSE: &str = "tout too close to t0 to start integration.";
pub const MSG_BAD_HINIT: &str = "Initial step is not towards tout.";
pub const MSG_BAD_TSTOP: &str = "The value tstop = %.15g is behind current t = %.15g\
in the direction of integration.";
pub const MSG_CLOSE_ROOTS: &str = "Root found at and very near t = %.15g.";
pub const MSG_MAX_STEPS: &str = "At t = %.15g, mxstep steps taken before reaching tout.";
pub const MSG_EWT_NOW_FAIL: &str = "At t = %.15gthe user-provide EwtSet function failed.";
pub const MSG_EWT_NOW_BAD: &str = "At t = %.15gsome ewt component has become <= 0.0.";
pub const MSG_TOO_MUCH_ACC: &str = "At t = %.15gtoo much accuracy requested.";

pub const MSG_BAD_K: &str = "Illegal value for k.";
pub const MSG_NULL_DKY: &str = "dky = NULL illegal.";
pub const MSG_NULL_DKYP: &str = "dkyp = NULL illegal.";
pub const MSG_BAD_T: &str =
    "Illegal value for t.t = %.15g is not between tcur - hold = %.15g and tcur = %.15g";
pub const MSG_BAD_TOUT: &str = "Trouble interpolating at tout = %.15g. \
tout too far back in direction of integration.";

pub const MSG_ERR_FAILS: &str =
    "At t = %.15g and h = %.15g, the error test failed repeatedly or with |h| = hmin.";
pub const MSG_CONV_FAILS: &str = "At t = %.15g and h = %.15g, \
the corrector convergence failed repeatedly or with |h| = hmin.";
pub const MSG_SETUP_FAILED: &str =
    "At t = %.15g, the linear solver setup failed unrecoverably.";
pub const MSG_SOLVE_FAILED: &str =
    "At t = %.15g, the linear solver solve failed unrecoverably.";
pub const MSG_REP_RES_ERR: &str = "At t = %.15g repeated recoverable residual errors.";
pub const MSG_RES_NONRECOV: &str =
    "At t = %.15g, the residual function failed unrecoverably.";
pub const MSG_FAILED_CONSTR: &str =
    "At t = %.15g, unable to satisfy inequality constraints.";
pub const MSG_RTFUNC_FAILED: &str =
    "At t = %.15g, the rootfinding routine failed in an unrecoverable manner.";
pub const MSG_NO_ROOT: &str = "Rootfinding was not initialized.";
pub const MSG_INACTIVE_ROOTS: &str =
    "At the end of the first step, there are still some root functions \
     identically 0. This warning will not be issued again.";
pub const MSG_NLS_INPUT_NULL: &str =
    "At t = %.15g, the nonlinear solver was passed a NULL input.";
pub const MSG_NLS_SETUP_FAILED: &str =
    "At t = %.15g, the nonlinear solver setup failed unrecoverably.";
pub const MSG_NLS_FAIL: &str =
    "At t = %.15g, the nonlinear solver failed in an unrecoverable manner.";

/* IDASet* / IDAGet* error messages */

pub const MSG_NEG_MAXORD: &str = "maxord <= 0 illegal.";
pub const MSG_BAD_MAXORD: &str = "Illegal attempt to increase maximum order.";
pub const MSG_NEG_HMAX: &str = "hmax < 0 illegal.";
pub const MSG_NEG_HMIN: &str = "hmin < 0 illegal.";
pub const MSG_NEG_EPCON: &str = "epcon <= 0.0 illegal.";
pub const MSG_BAD_CONSTR: &str = "Illegal values in constraints vector.";
pub const MSG_BAD_EPICCON: &str = "epiccon <= 0.0 illegal.";
pub const MSG_BAD_MAXNH: &str = "maxnh <= 0 illegal.";
pub const MSG_BAD_MAXNJ: &str = "maxnj <= 0 illegal.";
pub const MSG_BAD_MAXNIT: &str = "maxnit <= 0 illegal.";
pub const MSG_BAD_STEPTOL: &str = "steptol <= 0.0 illegal.";

pub const MSG_TOO_LATE: &str = "IDAGetConsistentIC can only be called before IDASolve.";

#[cfg(test)]
mod tests {
    use super::*;

    /* A Default IDAMem carries the C IDACreate default values
       (ida.c, IDACreate; constants from ida_impl.h and ida.c). */
    #[test]
    fn idamem_default_matches_idacreate() {
        let ida_mem = IDAMem::default();

        assert_eq!(ida_mem.ida_uround, SUN_UNIT_ROUNDOFF);

        assert_eq!(ida_mem.ida_itol, IDA_NN);
        assert!(ida_mem.ida_atolmin0);
        assert!(!ida_mem.ida_user_efun);

        assert_eq!(ida_mem.ida_maxord, MAXORD_DEFAULT); /* 5 */
        assert_eq!(ida_mem.ida_maxord, 5);
        assert_eq!(ida_mem.ida_maxord_alloc, MAXORD_DEFAULT);
        assert_eq!(ida_mem.ida_mxstep, MXSTEP_DEFAULT); /* 500 */
        assert_eq!(ida_mem.ida_mxstep, 500);
        assert_eq!(ida_mem.ida_hmax_inv, HMAX_INV_DEFAULT); /* 0.0 */
        assert_eq!(ida_mem.ida_hmin, HMIN_DEFAULT); /* 0.0 */

        assert_eq!(ida_mem.ida_eta_max_fx, 2.0);
        assert_eq!(ida_mem.ida_eta_min_fx, 1.0);
        assert_eq!(ida_mem.ida_eta_max, 2.0);
        assert_eq!(ida_mem.ida_eta_low, 0.9);
        assert_eq!(ida_mem.ida_eta_min, 0.5);
        assert_eq!(ida_mem.ida_eta_min_ef, 0.25);
        assert_eq!(ida_mem.ida_eta_cf, 0.25);

        assert_eq!(ida_mem.ida_hin, 0.0);
        assert_eq!(ida_mem.ida_epcon, EPCON); /* 0.33 */
        assert_eq!(ida_mem.ida_epcon, 0.33);
        assert_eq!(ida_mem.ida_maxnef, MXNEF); /* 10 */
        assert_eq!(ida_mem.ida_maxncf, MXNCF); /* 10 */
        assert!(!ida_mem.ida_suppressalg);
        assert!(!ida_mem.ida_tstopset);
        assert_eq!(ida_mem.ida_dcj, DCJ_DEFAULT); /* 0.25 */

        assert!(!ida_mem.ida_constraintsSet);
        assert_eq!(ida_mem.constraint_corrections, 0);
        assert_eq!(ida_mem.constraint_fails, 0);
        assert_eq!(ida_mem.max_constraint_fails, MAX_CONSTRAINT_FAILS); /* 10 */

        /* IC optional inputs */
        assert_eq!(ida_mem.ida_epiccon, 0.01 * EPCON);
        assert_eq!(ida_mem.ida_maxnh, MAXNH); /* 5   */
        assert_eq!(ida_mem.ida_maxnj, MAXNJ); /* 4   */
        assert_eq!(ida_mem.ida_maxnit, MAXNI); /* 10  */
        assert_eq!(ida_mem.ida_maxbacks, MAXBACKS); /* 100 */
        assert!(!ida_mem.ida_lsoff);
        assert_eq!(ida_mem.ida_steptol, SUNRpowerR(SUN_UNIT_ROUNDOFF, 0.667));

        /* lrw / liw */
        assert_eq!(ida_mem.ida_lrw, 25 + 5 * MXORDP1 as i64); /* 55 */
        assert_eq!(ida_mem.ida_lrw, 55);
        assert_eq!(ida_mem.ida_liw, 38);

        assert!(!ida_mem.ida_VatolMallocDone);
        assert!(!ida_mem.ida_idMallocDone);
        assert!(!ida_mem.ida_MallocDone);
        assert!(!ida_mem.ida_SetupDone);

        assert!(ida_mem.NLS.is_none());
        assert!(!ida_mem.ownNLS);
        assert!(ida_mem.ida_lmem.is_none());
        assert!(!ida_mem.ida_linitOK);

        /* memset(0) remainder (representative) */
        assert_eq!(ida_mem.ida_nst, 0);
        assert_eq!(ida_mem.ida_nre, 0);
        assert_eq!(ida_mem.ida_nrtfn, 0);
        assert_eq!(ida_mem.ida_psi, [0.0; MXORDP1]);
        assert!(ida_mem.ida_phi.is_empty());
    }
}
