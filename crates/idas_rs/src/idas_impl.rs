/* -----------------------------------------------------------------
 * Translated from src/idas/idas_impl.h and the constants / typedefs
 * of include/idas/idas.h (IDAS 7.7.0).
 * Main integrator memory block (struct IDAMemRec) with quadrature,
 * forward-sensitivity, quadrature-sensitivity and adjoint (IDAA)
 * extensions, plus the adjoint module memory blocks (IDAckpntMem,
 * IDAdtpntMem, IDABMem, IDAadjMem).
 * Field names keep the C names (ida_/ck_/ia_ prefixes) so the
 * translation of idas.c / idaa.c reads line-for-line.
 * Structural donor: ida_rs/src/ida_impl.rs (verified Phase 4);
 * adjoint-modeling decisions mirror cvodes_rs/src/cvodes_impl.rs
 * (pinned): linked lists -> owned Vecs, next-pointers dropped,
 * C void* dt content -> DtpntContent enum, interpolation-module
 * function pointers -> dispatch on ia_interpType, list-position
 * pointers -> Option<usize> indices.
 * -----------------------------------------------------------------*/
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

pub type IDAQuadRhsFn =
    fn(tres: f64, yy: &NVector, yp: &NVector, rrQ: &mut NVector, user_data: &mut UserData) -> i32;

pub type IDASensResFn = fn(
    Ns: i32,
    t: f64,
    yy: &NVector,
    yp: &NVector,
    resval: &NVector,
    yyS: &[NVector],
    ypS: &[NVector],
    resvalS: &mut [NVector],
    user_data: &mut UserData,
    tmp1: &mut NVector,
    tmp2: &mut NVector,
    tmp3: &mut NVector,
) -> i32;

pub type IDAQuadSensRhsFn = fn(
    Ns: i32,
    t: f64,
    yy: &NVector,
    yp: &NVector,
    yyS: &[NVector],
    ypS: &[NVector],
    rrQ: &NVector,
    rhsvalQS: &mut [NVector],
    user_data: &mut UserData,
    yytmp: &mut NVector,
    yptmp: &mut NVector,
    tmpQS: &mut NVector,
) -> i32;

/* Backward (adjoint) problem function types (idas.h) */

pub type IDAResFnB = fn(
    tt: f64,
    yy: &NVector,
    yp: &NVector,
    yyB: &NVector,
    ypB: &NVector,
    rrB: &mut NVector,
    user_dataB: &mut UserData,
) -> i32;

pub type IDAResFnBS = fn(
    t: f64,
    yy: &NVector,
    yp: &NVector,
    yyS: &[NVector],
    ypS: &[NVector],
    yyB: &NVector,
    ypB: &NVector,
    rrBS: &mut NVector,
    user_dataB: &mut UserData,
) -> i32;

pub type IDAQuadRhsFnB = fn(
    tt: f64,
    yy: &NVector,
    yp: &NVector,
    yyB: &NVector,
    ypB: &NVector,
    rhsvalBQ: &mut NVector,
    user_dataB: &mut UserData,
) -> i32;

pub type IDAQuadRhsFnBS = fn(
    t: f64,
    yy: &NVector,
    yp: &NVector,
    yyS: &[NVector],
    ypS: &[NVector],
    yyB: &NVector,
    ypB: &NVector,
    rhsvalBQS: &mut NVector,
    user_dataB: &mut UserData,
) -> i32;

/* ===============================================================
   IDA constants (ida.h)
   =============================================================== */

/* itask */
pub const IDA_NORMAL: i32 = 1;
pub const IDA_ONE_STEP: i32 = 2;

/* icopt */
pub const IDA_YA_YDP_INIT: i32 = 1;
pub const IDA_Y_INIT: i32 = 2;

/* ism (sensitivity corrector method) */
pub const IDA_SIMULTANEOUS: i32 = 1;
pub const IDA_STAGGERED: i32 = 2;

/* DQtype (sensitivity DQ approximation) */
pub const IDA_CENTERED: i32 = 1;
pub const IDA_FORWARD: i32 = 2;

/* interp (adjoint interpolation type) */
pub const IDA_HERMITE: i32 = 1;
pub const IDA_POLYNOMIAL: i32 = 2;

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

pub const IDA_NO_QUAD: i32 = -30;
pub const IDA_QRHS_FAIL: i32 = -31;
pub const IDA_FIRST_QRHS_ERR: i32 = -32;
pub const IDA_REP_QRHS_ERR: i32 = -33;

pub const IDA_NO_SENS: i32 = -40;
pub const IDA_SRES_FAIL: i32 = -41;
pub const IDA_REP_SRES_ERR: i32 = -42;
pub const IDA_BAD_IS: i32 = -43;

pub const IDA_NO_QUADSENS: i32 = -50;
pub const IDA_QSRHS_FAIL: i32 = -51;
pub const IDA_FIRST_QSRHS_ERR: i32 = -52;
pub const IDA_REP_QSRHS_ERR: i32 = -53;

pub const IDA_TOO_CLOSE: i32 = -60;

pub const IDA_UNRECOGNIZED_ERROR: i32 = -99;

/* Adjoint (IDAA) return values */
pub const IDA_NO_ADJ: i32 = -101;
pub const IDA_NO_FWD: i32 = -102;
pub const IDA_NO_BCK: i32 = -103;
pub const IDA_BAD_TB0: i32 = -104;
pub const IDA_REIFWD_FAIL: i32 = -105;
pub const IDA_FWD_FAIL: i32 = -106;
pub const IDA_GETY_BADT: i32 = -107;

/* itol (internal; ida.c) */
pub const IDA_NN: i32 = 0;
pub const IDA_SS: i32 = 1;
pub const IDA_SV: i32 = 2;
pub const IDA_WF: i32 = 3;
pub const IDA_EE: i32 = 4;

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

pub const IDA_QRHS_RECVR: i32 = 10;
pub const IDA_SRES_RECVR: i32 = 11;
pub const IDA_QSRHS_RECVR: i32 = 12;

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
    /// IDALS linear solver interface (idas_ls_impl.rs / idas_ls.rs)
    Ls(Box<crate::idas_ls_impl::IDALsMem>),
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

    /*-----------------------
      Quadrature Related Data
      -----------------------*/
    pub ida_quadr: bool,

    pub ida_rhsQ: Option<IDAQuadRhsFn>,
    /* (C ida_user_dataQ is not stored: rhsQ receives ida_user_data,
       per the donor's drop of ida_edata.) */

    pub ida_errconQ: bool,

    pub ida_itolQ: i32,
    pub ida_rtolQ: f64,
    pub ida_SatolQ: f64,     /* scalar absolute tolerance for quadratures  */
    pub ida_VatolQ: NVector, /* vector absolute tolerance for quadratures  */
    pub ida_atolQmin0: bool, /* flag indicating that min(atolQ) = 0        */

    /*------------------------
      Sensitivity Related Data
      ------------------------*/
    pub ida_sensi: bool,
    pub ida_Ns: i32,
    pub ida_ism: i32,

    pub ida_resS: Option<IDASensResFn>,
    /* (C ida_user_dataS defaults to IDA_mem itself for the internal DQ
       residual; not stored — the DQ path in idas.rs operates on
       &mut IDAMem directly and a user resS receives ida_user_data.) */
    pub ida_resSDQ: bool,

    pub ida_p: Vec<f64>,
    pub ida_pbar: Vec<f64>,
    pub ida_plist: Vec<i32>,
    pub ida_DQtype: i32,
    pub ida_DQrhomax: f64,

    pub ida_errconS: bool, /* SUNTRUE if sensitivities in err. control  */

    pub ida_itolS: i32,
    pub ida_rtolS: f64,            /* relative tolerance for sensitivities    */
    pub ida_SatolS: Vec<f64>,      /* scalar absolute tolerances for sensi.   */
    pub ida_VatolS: Vec<NVector>,  /* vector absolute tolerances for sensi.   */
    pub ida_atolSmin0: Vec<bool>,  /* flag indicating that min(atolS[is]) = 0 */

    /*-----------------------------------
      Quadrature Sensitivity Related Data
      -----------------------------------*/
    pub ida_quadr_sensi: bool, /* SUNTRUE if computing sensitivities of quadrs. */

    pub ida_rhsQS: Option<IDAQuadSensRhsFn>, /* fQS = (dfQ/dy)*yS + (dfQ/dp) */
    /* (C ida_user_dataQS defaults to IDA_mem for the internal DQ rhs;
       not stored, same convention as ida_user_dataS above.) */
    pub ida_rhsQSDQ: bool, /* SUNTRUE if using internal DQ functions        */

    pub ida_errconQS: bool, /* SUNTRUE if yQS are considered in err. con.    */

    pub ida_itolQS: i32,
    pub ida_rtolQS: f64,            /* relative tolerance for yQS                */
    pub ida_SatolQS: Vec<f64>,      /* scalar absolute tolerances for yQS        */
    pub ida_VatolQS: Vec<NVector>,  /* vector absolute tolerances for yQS        */
    pub ida_atolQSmin0: Vec<bool>,  /* flag indicating that min(atolQS[is]) = 0  */

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

    /*----------------------------
      Quadrature Related N_Vectors
      ----------------------------*/
    pub ida_phiQ: Vec<NVector>, /* (C: N_Vector ida_phiQ[MXORDP1])          */
    pub ida_yyQ: NVector,
    pub ida_ypQ: NVector,
    pub ida_ewtQ: NVector,
    pub ida_eeQ: NVector,

    /*---------------------------
      Sensitivity Related Vectors
      ---------------------------*/
    pub ida_phiS: Vec<Vec<NVector>>, /* (C: N_Vector* ida_phiS[MXORDP1]; Vec of
                                        maxcol+1 live rows per the ida_phi
                                        modeling — IDASensAllocVectors) */
    pub ida_ewtS: Vec<NVector>,

    pub ida_eeS: Vec<NVector>, /* cumulative sensitivity corrections        */

    pub ida_yyS: Vec<NVector>,        /* allocated and used for:            */
    pub ida_ypS: Vec<NVector>,        /*                 ism = SIMULTANEOUS */
    pub ida_yySpredict: Vec<NVector>, /*                 ism = STAGGERED    */
    pub ida_ypSpredict: Vec<NVector>,
    pub ida_deltaS: Vec<NVector>,

    /* (C ida_tmpS1 / ida_tmpS2 — resS work vectors — are aliases of
       ida_tempv1 / ida_tempv2 and are not stored, per the alias-drop
       convention above; only the allocated third one is a field.) */
    pub ida_tmpS3: NVector, /* allocated resS work vector (C: tmpS3)       */

    /* (C ida_savresS / ida_delnewS — IDACalcIC staggered work vectors —
       are aliases of phiS[2] / phiS[3]; C ida_yyS0new / ida_ypS0new —
       IDASensLineSrch work vectors — are aliases of phiS[4] / eeS.
       All four are dropped per the alias-drop convention above.) */

    pub ida_yyS0: Vec<NVector>, /* initial yS, ypS vectors allocated and    */
    pub ida_ypS0: Vec<NVector>, /* deallocated in IDACalcIC function        */

    /*--------------------------------------
      Quadrature Sensitivity Related Vectors
      --------------------------------------*/
    pub ida_phiQS: Vec<Vec<NVector>>, /* Mod. div. diffs. for quadr. sensi.
                                         (C: N_Vector* ida_phiQS[MXORDP1]; Vec of
                                         maxord+1 live rows — IDAQuadSensAllocVectors) */
    pub ida_ewtQS: Vec<NVector>, /* error weight vectors for sensitivities  */

    pub ida_eeQS: Vec<NVector>, /* cumulative quadr.sensi.corrections       */

    pub ida_yyQS: Vec<NVector>, /* Unlike yS, yQS is not allocated by the user */
    pub ida_tempvQS: Vec<NVector>, /* temporary storage vector (~ tempv)    */
    pub ida_savrhsQ: NVector, /* saved quadr. rhs (needed for rhsQS calls)  */

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

    pub ida_ssS: f64, /* scalar ss for staggered sensitivities             */

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

    /* (sensitivity / quadrature counters) */
    pub ida_nrQe: i64,
    pub ida_nrSe: i64,
    pub ida_nrQSe: i64, /* number of fQS calls                               */
    pub ida_nreS: i64,
    pub ida_nrQeS: i64, /* number of fQ calls from sensi DQ                  */

    pub ida_ncfnQ: i64,
    pub ida_ncfnS: i64,

    pub ida_netfQ: i64,
    pub ida_netfS: i64,
    pub ida_netfQS: i64, /* number of quadr. sensi. error test failures  */

    pub ida_nniS: i64,
    pub ida_nnfS: i64,

    pub ida_nsetupsS: i64,

    /*------------------
      Space requirements
      ------------------*/
    pub ida_lrw1: i64, /* no. of sunrealtype words in 1 N_Vector            */
    pub ida_liw1: i64, /* no. of integer words in 1 N_Vector                */
    pub ida_lrw1Q: i64,
    pub ida_liw1Q: i64,
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

    pub ida_VatolQMallocDone: bool,
    pub ida_quadMallocDone: bool,

    pub ida_VatolSMallocDone: bool,
    pub ida_SatolSMallocDone: bool,
    pub ida_sensMallocDone: bool,

    pub ida_VatolQSMallocDone: bool,
    pub ida_SatolQSMallocDone: bool,
    pub ida_quadSensMallocDone: bool,

    /*---------------------
      Nonlinear Solver Data
      ---------------------*/
    pub NLS: Option<NonlinearSolver>, /* nonlinear solver object       */
    pub ownNLS: bool,                 /* flag indicating NLS ownership */

    pub NLSsim: Option<NonlinearSolver>, /* NLS object for DAE+Sens solves with
                                         the simultaneous corrector option */
    pub ownNLSsim: bool,                 /* flag indicating NLS ownership */

    pub NLSstg: Option<NonlinearSolver>, /* NLS object for DAE+Sens solves with
                                         the staggered corrector option */
    pub ownNLSstg: bool,                 /* flag indicating NLS ownership */

    /* (The C senswrapper aliases ypredictSim/ycorSim/ewtSim and
       ypredictStg/ycorStg/ewtStg are NOT stored, per the cvodes_rs
       pinned decision: in C they merely wrap [delta,deltaS] /
       [ee,eeS] / [ewt,ewtS]; the idas_nls_sim / idas_nls_stg ports
       build SensWrapper views at the call sites.) */
    pub simMallocDone: bool,
    pub stgMallocDone: bool,
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

    /* Flag to request a call to the setup routine */
    pub ida_forceSetup: bool,

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
    /*------------------------
      Adjoint sensitivity data
      ------------------------*/
    pub ida_adj: bool, /* SUNTRUE if performing ASA              */

    pub ida_adj_mem: Option<Box<IDAadjMem>>, /* adjoint memory structure */

    pub ida_adjMallocDone: bool,

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

            /* Set default values for quad. optional inputs */
            ida_quadr: SUNFALSE,
            ida_rhsQ: None,
            ida_errconQ: SUNFALSE,
            ida_itolQ: IDA_NN,
            ida_rtolQ: 0.0,
            ida_SatolQ: 0.0,
            ida_VatolQ: NVector::default(),
            ida_atolQmin0: SUNTRUE,

            /* Set default values for sensi. optional inputs
               (C: ida_user_dataS = IDA_mem and ida_resS = IDASensResDQ;
               here the internal-DQ state is resS = None + resSDQ = true) */
            ida_sensi: SUNFALSE,
            ida_Ns: 0,
            ida_resS: None,
            ida_resSDQ: SUNTRUE,
            ida_DQtype: IDA_CENTERED,
            ida_DQrhomax: 0.0,
            ida_p: Vec::new(),
            ida_pbar: Vec::new(),
            ida_plist: Vec::new(),
            ida_errconS: SUNFALSE,
            ida_itolS: IDA_EE,
            ida_rtolS: 0.0,
            ida_SatolS: Vec::new(),
            ida_VatolS: Vec::new(),
            ida_atolSmin0: Vec::new(),
            ida_ism: -1, /* initialize to invalid option */

            /* Defaults for sensi. of quadratures */
            ida_quadr_sensi: SUNFALSE,
            ida_rhsQS: None,
            ida_rhsQSDQ: SUNTRUE,
            ida_errconQS: SUNFALSE,
            ida_itolQS: IDA_EE,
            ida_rtolQS: 0.0,
            ida_SatolQS: Vec::new(),
            ida_VatolQS: Vec::new(),
            ida_atolQSmin0: Vec::new(),
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
            ida_VatolQMallocDone: SUNFALSE,
            ida_quadMallocDone: SUNFALSE,
            ida_VatolSMallocDone: SUNFALSE,
            ida_SatolSMallocDone: SUNFALSE,
            ida_sensMallocDone: SUNFALSE,
            ida_VatolQSMallocDone: SUNFALSE,
            ida_SatolQSMallocDone: SUNFALSE,
            ida_quadSensMallocDone: SUNFALSE,

            /* Initialize nonlinear solver variables */
            NLS: None,
            ownNLS: SUNFALSE,
            NLSsim: None,
            ownNLSsim: SUNFALSE,
            NLSstg: None,
            ownNLSstg: SUNFALSE,
            simMallocDone: SUNFALSE,
            stgMallocDone: SUNFALSE,

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

            ida_phiQ: Vec::new(),
            ida_yyQ: NVector::default(),
            ida_ypQ: NVector::default(),
            ida_ewtQ: NVector::default(),
            ida_eeQ: NVector::default(),

            ida_phiS: Vec::new(),
            ida_ewtS: Vec::new(),
            ida_eeS: Vec::new(),
            ida_yyS: Vec::new(),
            ida_ypS: Vec::new(),
            ida_yySpredict: Vec::new(),
            ida_ypSpredict: Vec::new(),
            ida_deltaS: Vec::new(),
            ida_tmpS3: NVector::default(),
            ida_yyS0: Vec::new(),
            ida_ypS0: Vec::new(),

            ida_phiQS: Vec::new(),
            ida_ewtQS: Vec::new(),
            ida_eeQS: Vec::new(),
            ida_yyQS: Vec::new(),
            ida_tempvQS: Vec::new(),
            ida_savrhsQ: NVector::default(),

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
            ida_ssS: 0.0,

            ida_nst: 0,
            ida_nre: 0,
            ida_ncfn: 0,
            ida_netf: 0,
            ida_nni: 0,
            ida_nnf: 0,
            ida_nsetups: 0,
            ida_nrQe: 0,
            ida_nrSe: 0,
            ida_nrQSe: 0,
            ida_nreS: 0,
            ida_nrQeS: 0,
            ida_ncfnQ: 0,
            ida_ncfnS: 0,
            ida_netfQ: 0,
            ida_netfS: 0,
            ida_netfQS: 0,
            ida_nniS: 0,
            ida_nnfS: 0,
            ida_nsetupsS: 0,

            ida_lrw1: 0,
            ida_liw1: 0,
            ida_lrw1Q: 0,
            ida_liw1Q: 0,

            ida_tolsf: 0.0,

            ida_SetupDone: SUNFALSE,

            ida_lmem: LsModule::None,
            ida_forceSetup: SUNFALSE,
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

            ida_adj: SUNFALSE,
            ida_adj_mem: None,
            ida_adjMallocDone: SUNFALSE,
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
pub const MSG_BAD_ISM_CONSTR: &str = "Constraints can not be enforced while forward sensitivity is used with simultaneous method.";
pub const MSG_LSOLVE_NULL: &str = "The linear solver's solve routine is NULL.";
pub const MSG_LINIT_FAIL: &str = "The linear solver's init routine failed.";
pub const MSG_NLS_INIT_FAIL: &str = "The nonlinear solver's init routine failed.";

pub const MSG_NO_QUAD: &str = "Illegal attempt to call before calling IDAQuadInit.";
pub const MSG_BAD_EWTQ: &str = "Initial ewtQ has component(s) equal to zero (illegal).";
pub const MSG_BAD_ITOLQ: &str =
    "Illegal value for itolQ. The legal values are IDA_SS and IDA_SV.";
pub const MSG_NO_TOLQ: &str =
    "No integration tolerances for quadrature variables have been specified.";
pub const MSG_NULL_ATOLQ: &str = "atolQ = NULL illegal.";
pub const MSG_BAD_RTOLQ: &str = "rtolQ < 0 illegal.";
pub const MSG_BAD_ATOLQ: &str = "atolQ has negative component(s) (illegal).";

pub const MSG_NO_SENSI: &str = "Illegal attempt to call before calling IDASensInit.";
pub const MSG_BAD_EWTS: &str = "Initial ewtS has component(s) equal to zero (illegal).";
pub const MSG_BAD_ITOLS: &str =
    "Illegal value for itolS. The legal values are IDA_SS, IDA_SV, and IDA_EE.";
pub const MSG_NULL_ATOLS: &str = "atolS = NULL illegal.";
pub const MSG_BAD_RTOLS: &str = "rtolS < 0 illegal.";
pub const MSG_BAD_ATOLS: &str = "atolS has negative component(s) (illegal).";
pub const MSG_BAD_PBAR: &str = "pbar has zero component(s) (illegal).";
pub const MSG_BAD_PLIST: &str = "plist has negative component(s) (illegal).";
pub const MSG_BAD_NS: &str = "NS <= 0 illegal.";
pub const MSG_NULL_YYS0: &str = "yyS0 = NULL illegal.";
pub const MSG_NULL_YPS0: &str = "ypS0 = NULL illegal.";
pub const MSG_BAD_ISM: &str =
    "Illegal value for ism. Legal values are: IDA_SIMULTANEOUS and IDA_STAGGERED.";
pub const MSG_BAD_IS: &str = "Illegal value for is.";
pub const MSG_NULL_DKYA: &str = "dkyA = NULL illegal.";
pub const MSG_BAD_DQTYPE: &str =
    "Illegal value for DQtype. Legal values are: IDA_CENTERED and IDA_FORWARD.";
pub const MSG_BAD_DQRHO: &str = "DQrhomax < 0 illegal.";

pub const MSG_NULL_ABSTOLQS: &str = "abstolQS = NULL illegal parameter.";
pub const MSG_BAD_RELTOLQS: &str = "reltolQS < 0 illegal parameter.";
pub const MSG_BAD_ABSTOLQS: &str = "abstolQS has negative component(s) (illegal).";
pub const MSG_NO_QUADSENSI: &str =
    "Forward sensitivity analysis for quadrature variables was not activated.";
pub const MSG_NULL_YQS0: &str = "yQS0 = NULL illegal parameter.";

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

/* ===============================================================
   IDAS quadrature / sensitivity solve-time error messages
   (idas_impl.h); C printf formats kept verbatim with SUN_FORMAT_G
   expanded to "%.15g".
   =============================================================== */

pub const MSG_EWTQ_NOW_BAD: &str = "At t = %.15g, a component of ewtQ has become <= 0.";
pub const MSG_QRHSFUNC_FAILED: &str =
    "At t = %.15g, the quadrature right-hand side routine failed in an unrecoverable manner.";
pub const MSG_QRHSFUNC_UNREC: &str = "At t = %.15g, the quadrature right-hand side failed in a recoverable manner, but no recovery is possible.";
pub const MSG_QRHSFUNC_REPTD: &str =
    "At t = %.15grepeated recoverable quadrature right-hand side function errors.";
pub const MSG_QRHSFUNC_FIRST: &str =
    "The quadrature right-hand side routine failed at the first call.";

pub const MSG_NULL_P: &str =
    "p = NULL when using internal DQ for sensitivity residual is illegal.";
pub const MSG_EWTS_NOW_BAD: &str = "At t = %.15g, a component of ewtS has become <= 0.";
pub const MSG_SRES_FAILED: &str =
    "At t = %.15g, the sensitivity residual routine failed in an unrecoverable manner.";
pub const MSG_SRES_UNREC: &str = "At t = %.15g, the sensitivity residual failed in a recoverable manner, but no recovery is possible.";
pub const MSG_SRES_REPTD: &str =
    "At t = %.15grepeated recoverable sensitivity residual function errors.";

pub const MSG_NO_TOLQS: &str =
    "No integration tolerances for quadrature sensitivity variables have been specified.";
pub const MSG_NULL_RHSQ: &str = "IDAS is expected to use DQ to evaluate the RHS of quad. sensi., but quadratures were not initialized.";
pub const MSG_BAD_EWTQS: &str = "Initial ewtQS has component(s) equal to zero (illegal).";
pub const MSG_EWTQS_NOW_BAD: &str =
    "At t = %.15g, a component of ewtQS has become <= 0.";
pub const MSG_QSRHSFUNC_FAILED: &str = "At t = %.15g, the sensitivity quadrature right-hand side routine failed in an unrecoverable manner.";
pub const MSG_QSRHSFUNC_REPTD: &str = "At t = %.15grepeated recoverable sensitivity quadrature right-hand side function errors.";
pub const MSG_QSRHSFUNC_FIRST: &str =
    "The quadrature right-hand side routine failed at the first call.";

/* ===============================================================
   IDAA (adjoint) error messages (idas_impl.h)
   =============================================================== */

pub const MSGAM_NULL_IDAMEM: &str = "ida_mem = NULL illegal.";
pub const MSGAM_NO_ADJ: &str = "Illegal attempt to call before calling IDAadjInit.";
pub const MSGAM_BAD_INTERP: &str = "Illegal value for interp.";
pub const MSGAM_BAD_STEPS: &str = "Steps nonpositive illegal.";
pub const MSGAM_BAD_WHICH: &str = "Illegal value for which.";
pub const MSGAM_NO_BCK: &str = "No backward problems have been defined yet.";
pub const MSGAM_NO_FWD: &str = "Illegal attempt to call before calling IDASolveF.";
pub const MSGAM_BAD_TB0: &str = "The initial time tB0 is outside the interval over which the forward problem was solved.";
pub const MSGAM_BAD_SENSI: &str = "At least one backward problem requires sensitivities, but they were not stored for interpolation.";
pub const MSGAM_BAD_ITASKB: &str =
    "Illegal value for itaskB. Legal values are IDA_NORMAL and IDA_ONE_STEP.";
pub const MSGAM_BAD_TBOUT: &str = "The final time tBout is outside the interval over which the forward problem was solved.";
pub const MSGAM_BACK_ERROR: &str = "Error occurred while integrating backward problem # %d";
pub const MSGAM_BAD_TINTERP: &str = "Bad t = %.15g for interpolation.";
pub const MSGAM_BAD_T: &str = "Bad t for interpolation.";
pub const MSGAM_WRONG_INTERP: &str =
    "This function cannot be called for the specified interp type.";
pub const MSGAM_MEM_FAIL: &str = "A memory request failed.";
pub const MSGAM_NO_INITBS: &str = "Illegal attempt to call before calling IDAInitBS.";

/* ===============================================================
   A D J O I N T   M O D U L E   M E M O R Y   B L O C K S
   (idas_impl.h; modeling mirrors cvodes_impl.rs, pinned)
   =============================================================== */

/* IDAckpntMem (struct IDAckpntMemRec): all information at a check
   point needed to 'hot' start IDAS. The C ck_next pointer is
   dropped: IDAadjMem owns the check points as Vec<IDAckpntMem>. */
pub struct IDAckpntMem {
    /* Integration limits */
    pub ck_t0: f64,
    pub ck_t1: f64,

    /* Modified divided difference array */
    pub ck_phi: Vec<NVector>, /* (C: N_Vector ck_phi[MXORDP1]) */

    /* Do we need to carry quadratures? */
    pub ck_quadr: bool,

    /* Modified divided difference array for quadratures */
    pub ck_phiQ: Vec<NVector>,

    /* Do we need to carry sensitivities? */
    pub ck_sensi: bool,

    /* number of sensitivities */
    pub ck_Ns: i32,

    /* Modified divided difference array for sensitivities */
    pub ck_phiS: [Vec<NVector>; MXORDP1],

    /* Do we need to carry quadrature sensitivities? */
    pub ck_quadr_sensi: bool,

    /* Modified divided difference array for quadrature sensitivities */
    pub ck_phiQS: [Vec<NVector>; MXORDP1],

    /* Step data */
    pub ck_nst: i64,
    pub ck_tretlast: f64,
    pub ck_ns: i32,
    pub ck_kk: i32,
    pub ck_kused: i32,
    pub ck_knew: i32,
    pub ck_phase: i32,

    pub ck_hh: f64,
    pub ck_hused: f64,
    pub ck_eta: f64,
    pub ck_cj: f64,
    pub ck_cjlast: f64,
    pub ck_cjold: f64,
    pub ck_cjratio: f64,
    pub ck_ss: f64,
    pub ck_ssS: f64,

    pub ck_psi: [f64; MXORDP1],
    pub ck_alpha: [f64; MXORDP1],
    pub ck_beta: [f64; MXORDP1],
    pub ck_sigma: [f64; MXORDP1],
    pub ck_gamma: [f64; MXORDP1],

    /* How many phi, phiS, phiQ and phiQS were allocated? */
    pub ck_phi_alloc: i32,
    /* ck_next dropped (IDAadjMem holds Vec<IDAckpntMem>) */
}

/* IDAdtpntMem (struct IDAdtpntMemRec): all information at a data
   point needed to interpolate the forward solution. The C
   `void* content` (IDAhermiteDataMemRec / IDApolynomialDataMemRec)
   becomes the DtpntContent enum; the C interpolation-module function
   pointers (ia_storePnt / ia_getY / ia_malloc / ia_free) are
   replaced by dispatch on ia_interpType in idaa.rs (pinned). */
pub enum DtpntContent {
    None,
    /* Data for cubic Hermite interpolation (IDAhermiteDataMemRec) */
    Hermite {
        y: NVector,
        yd: NVector,
        yS: Vec<NVector>,
        ySd: Vec<NVector>,
    },
    /* Data for polynomial interpolation (IDApolynomialDataMemRec);
       yd / ySd store the derivative(s) only for the first dt point
       (C: NULL otherwise; here: None / empty Vec). */
    Polynomial {
        y: NVector,
        yS: Vec<NVector>,
        yd: Option<NVector>,
        ySd: Vec<NVector>,
        order: i32,
    },
}

pub struct IDAdtpntMem {
    pub t: f64,                /* time */
    pub content: DtpntContent, /* interpType-dependent content */
}

/* IDABMem (struct IDABMemRec): all information for ONE backward
   problem. The C ida_next pointer is dropped (IDAadjMem owns the
   backward problems as Vec<IDABMem>); the C ida_lfree/ida_pfree
   function pointers are dropped (Rust Drop frees attached memory);
   ida_lmem/ida_pmem are opaque placeholders refined by the
   idas_ls / idaa ports (pinned, as in cvodes_rs). */
pub struct IDABMem {
    /* Index of this backward problem */
    pub ida_index: i32,

    /* Time at which the backward problem is initialized. */
    pub ida_t0: f64,

    /* Memory for this backward problem */
    pub IDA_mem: Box<IDAMem>,

    /* Flags to indicate that this backward problem's RHS or quad RHS
     * require forward sensitivities */
    pub ida_res_withSensi: bool,
    pub ida_rhsQ_withSensi: bool,

    /* Residual function for backward run */
    pub ida_res: Option<IDAResFnB>,
    pub ida_resS: Option<IDAResFnBS>,

    /* Right hand side quadrature function (fQB) for backward run */
    pub ida_rhsQ: Option<IDAQuadRhsFnB>,
    pub ida_rhsQS: Option<IDAQuadRhsFnBS>,

    /* User user_data */
    pub ida_user_data: UserData,

    /* Memory block for a linear solver's interface to IDAA */
    pub ida_lmem: Option<Box<dyn std::any::Any>>,

    /* Memory block for a preconditioner's module interface to IDAA */
    pub ida_pmem: Option<Box<dyn std::any::Any>>,

    /* Time at which to extract solution / quadratures */
    pub ida_tout: f64,

    /* Workspace Nvectors */
    pub ida_yy: NVector,
    pub ida_yp: NVector,
    /* ida_next dropped (IDAadjMem holds Vec<IDABMem>) */
}

/* IDAadjMem (struct IDAadjMemRec): all information necessary for
   adjoint sensitivity analysis. Linked lists become owned Vecs;
   the C pointers ia_bckpbCrt and ia_ckpntData become Option<usize>
   indices into IDAB_mem / ck_mem; the interpolation-module function
   pointers become dispatch on ia_interpType (pinned). */
pub struct IDAadjMem {
    /* --------------------
     * Forward problem data
     * -------------------- */

    /* Integration interval */
    pub ia_tinitial: f64,
    pub ia_tfinal: f64,

    /* Flag for first call to IDASolveF */
    pub ia_firstIDAFcall: bool,

    /* Flag if IDASolveF was called with TSTOP */
    pub ia_tstopIDAFcall: bool,
    pub ia_tstopIDAF: f64,

    /* Flag if IDASolveF was called in IDA_NORMAL_MODE and encountered
       a root after tout */
    pub ia_rootret: bool,
    pub ia_troot: f64,

    /* ----------------------
     * Backward problems data
     * ---------------------- */

    /* Storage for backward problems (C: linked list headed at IDAB_mem) */
    pub IDAB_mem: Vec<IDABMem>,

    /* Number of backward problems. */
    pub ia_nbckpbs: i32,

    /* Current backward problem (C: pointer; here index into IDAB_mem) */
    pub ia_bckpbCrt: Option<usize>,

    /* Flag for first call to IDASolveB */
    pub ia_firstIDABcall: bool,

    /* ----------------
     * Check point data
     * ---------------- */

    /* Storage for check point information (C: linked list at ck_mem) */
    pub ck_mem: Vec<IDAckpntMem>,

    /* Check point for which data is available (C: pointer; here index
       into ck_mem) */
    pub ia_ckpntData: Option<usize>,

    /* Number of checkpoints. */
    pub ia_nckpnts: i32,

    /* ------------------
     * Interpolation data
     * ------------------ */

    /* Number of steps between 2 check points */
    pub ia_nsteps: i64,

    /* Last index used in IDAAfindIndex */
    pub ia_ilast: i64,

    /* Storage for data from forward runs */
    pub dt_mem: Vec<IDAdtpntMem>,

    /* Actual number of data points saved in current dt_mem */
    /* Commonly, np = nsteps+1                              */
    pub ia_np: i64,

    /* Interpolation type (IDA_HERMITE or IDA_POLYNOMIAL) */
    pub ia_interpType: i32,

    /* (The C interpolation-module function pointers ia_storePnt,
       ia_getY, ia_malloc and ia_free are replaced by dispatch on
       ia_interpType: idaa.rs implements the IDAAhermite* and
       IDAApolynomial* families and selects by ia_interpType at each
       call site, pinned as in cvodes_rs.) */

    /* Flags controlling the interpolation module */
    pub ia_mallocDone: bool,   /* IM initialized?                */
    pub ia_newData: bool,      /* new data available in dt_mem?  */
    pub ia_storeSensi: bool,   /* store sensitivities?           */
    pub ia_interpSensi: bool,  /* interpolate sensitivities?     */

    pub ia_noInterp: bool, /* interpolations are temporarily */
                           /* disabled ( IDACalcICB )        */

    /* Workspace for polynomial interpolation.
       In C, ia_Y[i]/ia_YS[i] are POINTERS into phi[i]/phiS[i]; here
       they are owned scratch storage (empty = C NULL) and the idaa.rs
       port copies data instead of aliasing (pinned, as in cvodes_rs). */
    pub ia_Y: Vec<NVector>,             /* owned scratch (C: N_Vector ia_Y[MXORDP1])    */
    pub ia_YS: [Vec<NVector>; MXORDP1], /* owned scratch (C: N_Vector* ia_YS[MXORDP1])  */
    pub ia_T: [f64; MXORDP1],

    /* Workspace for wrapper functions */
    pub ia_yyTmp: NVector,
    pub ia_ypTmp: NVector,
    pub ia_yySTmp: Vec<NVector>,
    pub ia_ypSTmp: Vec<NVector>,
}

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
