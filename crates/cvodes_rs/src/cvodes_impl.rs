/* -----------------------------------------------------------------
 * Translated from src/cvodes/cvodes_impl.h (CVODES 7.7.0).
 * Main integrator memory block (struct CVodeMemRec), the adjoint
 * module memory blocks (CVadjMemRec, CVckpntMemRec, CVdtpntMemRec,
 * CVodeBMemRec) and internal constants. Field names keep the C
 * names (cv_/ck_/ca_ prefixes) so the translations of cvodes.c and
 * cvodea.c read line-for-line. Baseline representation decisions
 * are copied verbatim from the verified CVODE donor port
 * (crates/cvode_rs/src/cvode_impl.rs).
 *
 * PINNED design decisions for the CVODES extensions:
 *
 * 1. `N_Vector*` (array-of-vectors) fields -> `Vec<NVector>`;
 *    `N_Vector cv_znS[L_MAX]`-style arrays -> `[Vec<NVector>; L_MAX]`
 *    (an empty Vec plays the role of a C NULL pointer). Optional
 *    single vectors (cv_Vabstol, cv_yQ, ...) follow the donor
 *    convention: a plain NVector whose default (empty) state is the
 *    C NULL.
 * 2. Callback typedefs are plain `fn` pointers following the
 *    workspace ABI (ARCHITECTURE.md 3.6): `user_data: &mut UserData`
 *    replaces `void* user_data`, output vectors are `&mut`, and the
 *    C parameter order is preserved otherwise.
 * 3. The four nonlinear solvers NLS / NLSsim / NLSstg / NLSstg1 are
 *    each `Option<NonlinearSolver>` plus their own ownership flag.
 *    The C N_Vector senswrapper aliases (zn0Sim, ycorSim, ewtSim,
 *    zn0Stg, ycorStg, ewtStg) are NOT stored: in C they merely wrap
 *    pointers to cv_zn[0]/cv_znS[0], cv_acor/cv_acorS and
 *    cv_ewt/cv_ewtS; per the workspace pattern the cvodes_nls_sim /
 *    cvodes_nls_stg / cvodes_nls_stg1 modules operate directly on
 *    the CVodeMem fields, so only the simMallocDone/stgMallocDone
 *    flags are kept.
 * 4. Adjoint structures: the C linked lists become owned Vecs.
 *    CVckpntMemRec -> CVckpntMem (ck_next dropped; CVadjMem holds
 *    Vec<CVckpntMem>). CVdtpntMemRec's `void* content` becomes the
 *    DtpntContent enum (Hermite/Polynomial variants replace the
 *    CVhermiteDataMemRec/CVpolynomialDataMemRec payload structs).
 *    CVodeBMemRec -> CVodeBMem (cv_next dropped; cv_lfree/cv_pfree
 *    fn pointers dropped -- Rust Drop frees the memory). CVadjMemRec
 *    -> CVadjMem with cvB_mem/ck_mem/dt_mem as Vecs and the C
 *    pointers ca_ckpntData/ca_bckpbCrt as Option<usize> indices into
 *    those Vecs. The interpolation-module function pointers
 *    (ca_IMmalloc/ca_IMfree/ca_IMstore/ca_IMget) are replaced by
 *    dispatch on ca_IMtype (CV_HERMITE/CV_POLYNOMIAL): cvodea.rs
 *    implements the cvaHermite... and cvaPolynomial... families and
 *    dispatches directly.
 *    The C workspace pointers ca_Y[i]/ca_YS[i] alias zn[i]/znS[i];
 *    here they are owned scratch vectors and the cvodea.rs port
 *    copies data instead of aliasing.
 * 5. The `void* python` field is Python-binding plumbing -- omitted.
 *    The `void* cv_e_data`, `void* cv_fS_data` and `void* cv_fQS_data`
 *    fields point at cv_mem itself when the internal DQ routines are
 *    used; as in the donor (which omits cv_e_data) they are omitted:
 *    the cv_user_efun / cv_fSDQ / cv_fQSDQ flags steer dispatch to
 *    the internal routines, which receive &mut CVodeMem directly.
 * 6. The cv_linit/cv_lreinit/cv_lsetup/cv_lsolve/cv_lfree function-
 *    pointer table plus `void* cv_lmem` become the donor's LsModule
 *    dispatching enum (see below).
 * 7. Construction: as in the donor, CVodeMem has no Default/new()
 *    here -- cvodes.rs (the cvodes.c translation) builds the struct
 *    literal inside CVodeCreate exactly like the donor's cvode.rs.
 * 8. Error message MSG* constants are ported verbatim; messages that
 *    embed printf format directives (MSG_TIME et al.) are built
 *    inline with format! at the call sites, as the donor does.
 *
 * Further donor-precedent omissions (mirroring cvode_impl.rs):
 *  - cv_monitorfun / cv_monitor_interval: monitoring is not enabled
 *    in this build (C build without SUNDIALS_BUILD_WITH_MONITORING);
 *    cvodes_io.rs defines CVMonitorFn and the erroring setters, as
 *    the donor's cvode_io.rs does.
 *  - cv_cvals / cv_Xvecs / cv_Zvecs: fused-vector-operation scratch
 *    pointer arrays; fused ops are expanded to plain loops in the
 *    pure-Rust build, so no scratch aliasing arrays exist.
 * -----------------------------------------------------------------*/
use crate::nvector_serial::NVector;
use crate::sundials_context::SUNContext;
use crate::sundials_nonlinearsolver::NonlinearSolver;
use crate::sundials_types::UserData;

/* ===============================================================
   User-supplied function types (cvodes.h)
   =============================================================== */

pub type CVRhsFn = fn(t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32;
pub type CVRootFn = fn(t: f64, y: &NVector, gout: &mut [f64], user_data: &mut UserData) -> i32;
pub type CVEwtFn = fn(y: &NVector, ewt: &mut NVector, user_data: &mut UserData) -> i32;

/* Quadrature RHS: yQdot = fQ(t, y) */
pub type CVQuadRhsFn =
    fn(t: f64, y: &NVector, yQdot: &mut NVector, user_data: &mut UserData) -> i32;

/* Sensitivity RHS (all sensitivities at once): fS = (df/dy)*yS + (df/dp) */
pub type CVSensRhsFn = fn(
    Ns: i32,
    t: f64,
    y: &NVector,
    ydot: &NVector,
    yS: &[NVector],
    ySdot: &mut [NVector],
    user_data: &mut UserData,
    tmp1: &mut NVector,
    tmp2: &mut NVector,
) -> i32;

/* Sensitivity RHS (one sensitivity at a time): fS1 = (df/dy)*yS_i + (df/dp_i) */
pub type CVSensRhs1Fn = fn(
    Ns: i32,
    t: f64,
    y: &NVector,
    ydot: &NVector,
    is: i32,
    yS: &NVector,
    ySdot: &mut NVector,
    user_data: &mut UserData,
    tmp1: &mut NVector,
    tmp2: &mut NVector,
) -> i32;

/* Quadrature sensitivity RHS: fQS = (dfQ/dy)*yS + (dfQ/dp) */
pub type CVQuadSensRhsFn = fn(
    Ns: i32,
    t: f64,
    y: &NVector,
    yS: &[NVector],
    yQdot: &NVector,
    yQSdot: &mut [NVector],
    user_data: &mut UserData,
    tmp: &mut NVector,
    tmpQ: &mut NVector,
) -> i32;

/* Backward problem RHS */
pub type CVRhsFnB =
    fn(t: f64, y: &NVector, yB: &NVector, yBdot: &mut NVector, user_dataB: &mut UserData) -> i32;

/* Backward problem RHS depending on forward sensitivities */
pub type CVRhsFnBS = fn(
    t: f64,
    y: &NVector,
    yS: &[NVector],
    yB: &NVector,
    yBdot: &mut NVector,
    user_dataB: &mut UserData,
) -> i32;

/* Backward problem quadrature RHS */
pub type CVQuadRhsFnB =
    fn(t: f64, y: &NVector, yB: &NVector, qBdot: &mut NVector, user_dataB: &mut UserData) -> i32;

/* Backward problem quadrature RHS depending on forward sensitivities */
pub type CVQuadRhsFnBS = fn(
    t: f64,
    y: &NVector,
    yS: &[NVector],
    yB: &NVector,
    qBdot: &mut NVector,
    user_dataB: &mut UserData,
) -> i32;

/* ===============================================================
   CVODES constants (cvodes.h)
   =============================================================== */

/* lmm */
pub const CV_ADAMS: i32 = 1;
pub const CV_BDF: i32 = 2;

/* itask */
pub const CV_NORMAL: i32 = 1;
pub const CV_ONE_STEP: i32 = 2;

/* ism */
pub const CV_SIMULTANEOUS: i32 = 1;
pub const CV_STAGGERED: i32 = 2;
pub const CV_STAGGERED1: i32 = 3;

/* DQtype */
pub const CV_CENTERED: i32 = 1;
pub const CV_FORWARD: i32 = 2;

/* interp */
pub const CV_HERMITE: i32 = 1;
pub const CV_POLYNOMIAL: i32 = 2;

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

pub const CV_NO_QUAD: i32 = -30;
pub const CV_QRHSFUNC_FAIL: i32 = -31;
pub const CV_FIRST_QRHSFUNC_ERR: i32 = -32;
pub const CV_REPTD_QRHSFUNC_ERR: i32 = -33;
pub const CV_UNREC_QRHSFUNC_ERR: i32 = -34;

pub const CV_NO_SENS: i32 = -40;
pub const CV_SRHSFUNC_FAIL: i32 = -41;
pub const CV_FIRST_SRHSFUNC_ERR: i32 = -42;
pub const CV_REPTD_SRHSFUNC_ERR: i32 = -43;
pub const CV_UNREC_SRHSFUNC_ERR: i32 = -44;

pub const CV_BAD_IS: i32 = -45;

pub const CV_NO_QUADSENS: i32 = -50;
pub const CV_QSRHSFUNC_FAIL: i32 = -51;
pub const CV_FIRST_QSRHSFUNC_ERR: i32 = -52;
pub const CV_REPTD_QSRHSFUNC_ERR: i32 = -53;
pub const CV_UNREC_QSRHSFUNC_ERR: i32 = -54;

pub const CV_CONTEXT_ERR: i32 = -55;

pub const CV_PROJ_MEM_NULL: i32 = -56;
pub const CV_PROJFUNC_FAIL: i32 = -57;
pub const CV_REPTD_PROJFUNC_ERR: i32 = -58;

pub const CV_BAD_TINTERP: i32 = -59;

pub const CV_UNRECOGNIZED_ERR: i32 = -99;

/* adjoint return values */
pub const CV_NO_ADJ: i32 = -101;
pub const CV_NO_FWD: i32 = -102;
pub const CV_NO_BCK: i32 = -103;
pub const CV_BAD_TB0: i32 = -104;
pub const CV_REIFWD_FAIL: i32 = -105;
pub const CV_FWD_FAIL: i32 = -106;
pub const CV_GETY_BADT: i32 = -107;

/* ===============================================================
   Internal constants (cvodes_impl.h)
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

pub const QRHSFUNC_RECVR: i32 = 13;
pub const SRHSFUNC_RECVR: i32 = 14;
pub const QSRHSFUNC_RECVR: i32 = 15;

/* nonlinear solver constants
   NLS_MAXCOR  maximum no. of corrector iterations for the nonlinear solver
   CRDOWN      constant used in the estimation of the convergence rate (crate)
               of the iterates for the nonlinear equation
   RDIV        declare divergence if ratio del/delp > RDIV
   (In the CVODE donor these live privately in cvode_nls.rs because
   cvode_impl.h does not define them; cvodes_impl.h defines them
   centrally for the four cvodes_nls* modules, so they are pub here.) */
pub const NLS_MAXCOR: i32 = 3;
pub const CRDOWN: f64 = 0.3;
pub const RDIV: f64 = 2.0;

/* Constants for convfail (input to cv_lsetup) */
pub const CV_NO_FAILURES: i32 = 0;
pub const CV_FAIL_BAD_J: i32 = 1;
pub const CV_FAIL_OTHER: i32 = 2;

/* ===============================================================
   Linear solver module attached to CVODES.
   In C this is the (cv_linit, cv_lreinit, cv_lsetup, cv_lsolve,
   cv_lfree) function-pointer five-tuple plus the void* cv_lmem;
   here it is an enum dispatched in cvodes.rs (donor: LsModule in
   cvode_impl.rs). The module is take()n out of CVodeMem for the
   duration of a call so its routines may borrow the integrator
   memory mutably.

   NOTE (phased port): the Ls(Box<CVLsMem>) and Diag(Box<CVDiagMem>)
   variants are added when cvodes_ls_impl.rs / cvodes_diag_impl.rs
   land (they mirror the donor's cvode_ls_impl.rs /
   cvode_diag_impl.rs); until then only the None state exists.
   =============================================================== */

#[derive(Default)]
pub enum LsModule {
    #[default]
    None,
    /// cvodes_ls.rs interface (CVodeSetLinearSolver)
    Ls(Box<crate::cvodes_ls_impl::CVLsMem>),
    /// cvodes_diag.rs diagonal approximation (CVDiag)
    Diag(Box<crate::cvodes_diag_impl::CVDiagMem>),
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

    /* void* python omitted: Python-binding plumbing (decision 5) */

    pub cv_uround: f64, /* machine unit roundoff */

    /*--------------------------
      Problem Specification Data
      --------------------------*/
    pub cv_f: Option<CVRhsFn>,  /* y' = f(t,y(t))                    */
    pub cv_user_data: UserData, /* user pointer passed to f          */
    pub cv_lmm: i32,            /* lmm = CV_ADAMS or CV_BDF          */
    pub cv_itol: i32,           /* itol = CV_SS, CV_SV, CV_WF, CV_NN */

    pub cv_reltol: f64,           /* relative tolerance                */
    pub cv_Sabstol: f64,          /* scalar absolute tolerance         */
    pub cv_Vabstol: NVector,      /* vector absolute tolerance         */
    pub cv_atolmin0: bool,        /* flag indicating that min(abstol) = 0 */
    pub cv_user_efun: bool,       /* SUNTRUE if user sets efun         */
    pub cv_efun: Option<CVEwtFn>, /* function to set ewt               */
    /* void* cv_e_data omitted (decision 5): internal cvEwtSet gets &mut CVodeMem */

    /*-----------------------
      Quadrature Related Data
      -----------------------*/
    pub cv_quadr: bool, /* SUNTRUE if integrating quadratures            */

    pub cv_fQ: Option<CVQuadRhsFn>, /* q' = fQ(t, y(t))                  */

    pub cv_errconQ: bool, /* SUNTRUE if quadrs. are included in error test */

    pub cv_itolQ: i32,        /* itolQ = CV_SS or CV_SV                        */
    pub cv_reltolQ: f64,      /* relative tolerance for quadratures            */
    pub cv_SabstolQ: f64,     /* scalar absolute tolerance for quadratures     */
    pub cv_VabstolQ: NVector, /* vector absolute tolerance for quadratures     */
    pub cv_atolQmin0: bool,   /* flag indicating that min(abstolQ) = 0         */

    /*------------------------
      Sensitivity Related Data
      ------------------------*/
    pub cv_sensi: bool, /* SUNTRUE if computing sensitivities           */

    pub cv_Ns: i32, /* Number of sensitivities                      */

    pub cv_ism: i32, /* ism = SIMULTANEOUS or STAGGERED              */

    pub cv_fS: Option<CVSensRhsFn>,   /* fS = (df/dy)*yS + (df/dp)      */
    pub cv_fS1: Option<CVSensRhs1Fn>, /* fS1 = (df/dy)*yS_i + (df/dp)   */
    /* void* cv_fS_data omitted (decision 5): internal DQ gets &mut CVodeMem */
    pub cv_fSDQ: bool, /* SUNTRUE if using internal DQ functions       */
    pub cv_ifS: i32,   /* ifS = ALLSENS or ONESENS                     */

    pub cv_p: Vec<f64>,     /* parameters in f(t,y,p)                       */
    pub cv_pbar: Vec<f64>,  /* scale factors for parameters                 */
    pub cv_plist: Vec<i32>, /* list of sensitivities                        */
    pub cv_DQtype: i32,     /* central/forward finite differences           */
    pub cv_DQrhomax: f64,   /* cut-off value for separate/simultaneous FD   */

    pub cv_errconS: bool, /* SUNTRUE if yS are considered in err. control */

    pub cv_itolS: i32,
    pub cv_reltolS: f64,           /* relative tolerance for sensitivities      */
    pub cv_SabstolS: Vec<f64>,     /* scalar absolute tolerances for sensi.     */
    pub cv_VabstolS: Vec<NVector>, /* vector absolute tolerances for sensi.     */
    pub cv_atolSmin0: Vec<bool>,   /* flags indicating that min(abstolS[i]) = 0 */

    /*-----------------------------------
      Quadrature Sensitivity Related Data
      -----------------------------------*/
    pub cv_quadr_sensi: bool, /* SUNTRUE if computing sensitivities of quadrs. */

    pub cv_fQS: Option<CVQuadSensRhsFn>, /* fQS = (dfQ/dy)*yS + (dfQ/dp)  */
    /* void* cv_fQS_data omitted (decision 5): internal DQ gets &mut CVodeMem */
    pub cv_fQSDQ: bool, /* SUNTRUE if using internal DQ functions       */

    pub cv_errconQS: bool, /* SUNTRUE if yQS are considered in err. con.   */

    pub cv_itolQS: i32,
    pub cv_reltolQS: f64,           /* relative tolerance for yQS                 */
    pub cv_SabstolQS: Vec<f64>,     /* scalar absolute tolerances for yQS         */
    pub cv_VabstolQS: Vec<NVector>, /* vector absolute tolerances for yQS         */
    pub cv_atolQSmin0: Vec<bool>,   /* flags indicating that min(abstolQS[i]) = 0 */

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

    /*--------------------------
      Quadrature Related Vectors
      --------------------------*/
    pub cv_znQ: Vec<NVector>, /* Nordsieck arrays for quadratures             */
    pub cv_ewtQ: NVector,     /* error weight vector for quadratures          */
    pub cv_yQ: NVector,       /* Unlike y, yQ is not allocated by the user    */
    pub cv_acorQ: NVector,    /* acorQ = yQ_n(m) - yQ_n(0)                    */
    pub cv_tempvQ: NVector,   /* temporary storage vector (~ tempv)           */

    /*---------------------------
      Sensitivity Related Vectors
      ---------------------------*/
    pub cv_znS: [Vec<NVector>; L_MAX], /* Nordsieck arrays for sensitivities  */
    pub cv_ewtS: Vec<NVector>,         /* error weight vectors for sensitivities */
    pub cv_yS: Vec<NVector>,           /* yS=yS0 (allocated by the user)      */
    pub cv_acorS: Vec<NVector>,        /* acorS = yS_n(m) - yS_n(0)           */
    pub cv_tempvS: Vec<NVector>,       /* temporary storage vector (~ tempv)  */
    pub cv_ftempS: Vec<NVector>,       /* temporary storage vector (~ ftemp)  */

    pub cv_stgr1alloc: bool, /* Did we allocate ncfS1, ncfnS1, and nniS1?    */

    /*--------------------------------------
      Quadrature Sensitivity Related Vectors
      --------------------------------------*/
    pub cv_znQS: [Vec<NVector>; L_MAX], /* Nordsieck arrays for quadr. sensitivities */
    pub cv_ewtQS: Vec<NVector>,         /* error weight vectors for sensitivities */
    pub cv_yQS: Vec<NVector>,           /* Unlike yS, yQS is not allocated by the user */
    pub cv_acorQS: Vec<NVector>,        /* acorQS = yQS_n(m) - yQS_n(0)       */
    pub cv_tempvQS: Vec<NVector>,       /* temporary storage vector (~ tempv) */
    pub cv_ftempQ: NVector,             /* temporary storage vector (~ ftemp) */

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

    pub cv_crate: f64,      /* estimated corrector convergence rate        */
    pub cv_crateS: f64,     /* estimated corrector convergence rate (Stgr) */
    pub cv_delp: f64,       /* norm of previous nonlinear solver update    */
    pub cv_acnrm: f64,      /* | acor |                                    */
    pub cv_acnrmcur: bool,  /* is | acor | current?                        */
    pub cv_acnrmQ: f64,     /* | acorQ |                                   */
    pub cv_acnrmS: f64,     /* | acorS |                                   */
    pub cv_acnrmScur: bool, /* is | acorS | current?                       */
    pub cv_acnrmQS: f64,    /* | acorQS |                                  */
    pub cv_nlscoef: f64,    /* coefficient in nonlinear convergence test   */
    pub cv_ncfS1: Vec<i64>, /* Array of Ns local counters for conv.
                             * failures (used in CVStep for STAGGERED1)
                             * (C: int*, widened to i64 per pinned decision) */

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
    pub cv_nst: i64, /* number of internal steps taken                  */

    pub cv_nfe: i64,   /* number of f calls                               */
    pub cv_nfQe: i64,  /* number of fQ calls                              */
    pub cv_nfSe: i64,  /* number of fS calls                              */
    pub cv_nfeS: i64,  /* number of f calls from sensi DQ                 */
    pub cv_nfQSe: i64, /* number of fQS calls                             */
    pub cv_nfQeS: i64, /* number of fQ calls from sensi DQ                */

    pub cv_ncfn: i64,        /* number of corrector convergence failures    */
    pub cv_ncfnS: i64,       /* number of total sensi. corr. conv. failures */
    pub cv_ncfnS1: Vec<i64>, /* number of sensi. corrector conv. failures   */

    pub cv_nni: i64,        /* number of nonlinear iterations performed     */
    pub cv_nniS: i64,       /* number of total sensi. nonlinear iterations  */
    pub cv_nniS1: Vec<i64>, /* number of sensi. nonlinear iterations        */

    pub cv_nnf: i64,        /* number of nonlinear convergence failures     */
    pub cv_nnfS: i64,       /* number of total sensi. nonlinear conv. fails */
    pub cv_nnfS1: Vec<i64>, /* number of sensi. nonlinear conv. fails       */

    pub cv_netf: i64,   /* number of error test failures                   */
    pub cv_netfQ: i64,  /* number of quadr. error test failures            */
    pub cv_netfS: i64,  /* number of sensi. error test failures            */
    pub cv_netfQS: i64, /* number of quadr. sensi. error test failures     */

    pub cv_nsetups: i64,  /* number of setup calls                          */
    pub cv_nsetupsS: i64, /* number of setup calls due to sensitivities     */

    pub cv_nhnil: i32, /* number of t + h == t messages issued            */

    /*----------------
      Step size ratios
      ----------------*/
    pub cv_etaqm1: f64, /* ratio of new to old h for order q-1 */
    pub cv_etaq: f64,   /* ratio of new to old h for order q   */
    pub cv_etaqp1: f64, /* ratio of new to old h for order q+1 */

    /*------------------
      Space requirements
      ------------------*/
    pub cv_lrw1: i64,  /* no. of sunrealtype words in 1 N_Vector y   */
    pub cv_liw1: i64,  /* no. of integer words in 1 N_Vector y       */
    pub cv_lrw1Q: i64, /* no. of sunrealtype words in 1 N_Vector yQ  */
    pub cv_liw1Q: i64, /* no. of integer words in 1 N_Vector yQ      */
    pub cv_lrw: i64,   /* no. of sunrealtype words in CVODES vectors */
    pub cv_liw: i64,   /* no. of integer words in CVODES vectors     */

    /*---------------------
      Nonlinear Solver Data
      ---------------------*/
    pub NLS: Option<NonlinearSolver>, /* nonlinear solver object          */
    /* mirror of SUNNonlinSolGetCurIter for the linear-solver interface
       (the NLS object is detached from CVodeMem during its solve);
       donor adaptation carried over */
    pub cv_nls_curiter: i32,
    pub ownNLS: bool, /* flag indicating NLS ownership    */

    pub NLSsim: Option<NonlinearSolver>, /* NLS object for the simultaneous corrector */
    pub ownNLSsim: bool,                 /* flag indicating NLS ownership             */

    pub NLSstg: Option<NonlinearSolver>, /* NLS object for the staggered corrector */
    pub ownNLSstg: bool,                 /* flag indicating NLS ownership          */

    pub NLSstg1: Option<NonlinearSolver>, /* NLS object for the staggered1 corrector */
    pub ownNLSstg1: bool,                 /* flag indicating NLS ownership           */
    pub sens_solve_idx: i32,              /* index of the current staggered1 solve   */
    pub nnip: i64,                        /* previous total number of iterations     */

    pub sens_solve: bool, /* flag indicating if the current solve is a
                             staggered or staggered1 sensitivity solve */
    pub nls_f: Option<CVRhsFn>, /* f(t,y(t)) used in the nonlinear solver */
    pub convfail: i32,    /* flag: Jacobian update may be needed         */

    /* The C N_Vector senswrapper aliases zn0Sim/ycorSim/ewtSim and
       zn0Stg/ycorStg/ewtStg are NOT stored (pinned decision 3): in C
       they wrap pointers to cv_zn[0]/cv_znS[0], cv_acor/cv_acorS and
       cv_ewt/cv_ewtS; the cvodes_nls_sim/stg/stg1 modules operate
       directly on those CVodeMem fields. Only the allocation flags
       survive: */
    pub simMallocDone: bool,
    pub stgMallocDone: bool,

    /*------------------
      Linear Solver Data
      ------------------*/
    /* In C: cv_linit/cv_lreinit/cv_lsetup/cv_lsolve/cv_lfree function
       pointers + void* cv_lmem. Here: dispatching enum (see LsModule). */
    pub cv_lmem: LsModule,
    pub cv_msbp: i64,         /* max number of steps between lsetup calls */
    pub cv_dgmax_lsetup: f64, /* gamma ratio threshold for lsetup         */
    pub cv_forceSetup: bool,  /* flag to request a call to the setup routine */

    /*------------
      Saved Values
      ------------*/
    pub cv_qu: i32,           /* last successful q value used                */
    pub cv_nstlp: i64,        /* step number of last setup call              */
    pub cv_h0u: f64,          /* actual initial stepsize                     */
    pub cv_hu: f64,           /* last successful h value used                */
    pub cv_saved_tq5: f64,    /* saved value of tq[5]                        */
    pub cv_jcur: bool,        /* is Jacobian info for linear solver current? */
    pub cv_convfail: i32,     /* flag storing previous solver failure mode   */
    pub cv_tolsf: f64,        /* tolerance scale factor                      */
    pub cv_qmax_alloc: i32,   /* value of qmax used when allocating memory   */
    pub cv_qmax_allocQ: i32,  /* qmax used when allocating quad. mem         */
    pub cv_qmax_allocS: i32,  /* qmax used when allocating sensi. mem        */
    pub cv_qmax_allocQS: i32, /* qmax used when allocating quad. sensi. mem  */
    pub cv_indx_acor: i32,    /* index of the zn vector with saved acor      */

    /*--------------------------------------------------------------------
      Flags turned ON by CVodeInit, CVodeSensMalloc, and CVodeQuadMalloc
      and read by CVodeReInit, CVodeSensReInit, and CVodeQuadReInit
      --------------------------------------------------------------------*/
    pub cv_VabstolMallocDone: bool,
    pub cv_MallocDone: bool,
    pub cv_constraintsMallocDone: bool,

    pub cv_VabstolQMallocDone: bool,
    pub cv_QuadMallocDone: bool,

    pub cv_VabstolSMallocDone: bool,
    pub cv_SabstolSMallocDone: bool,
    pub cv_SensMallocDone: bool,

    pub cv_VabstolQSMallocDone: bool,
    pub cv_SabstolQSMallocDone: bool,
    pub cv_QuadSensMallocDone: bool,

    /* cv_monitorfun / cv_monitor_interval omitted per donor precedent:
       monitoring is not enabled in this build (see module header). */

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
    pub cv_irfnd: i32,             /* flag showing whether last step had a root */
    pub cv_nge: i64,               /* counter for g evaluations            */
    pub cv_gactive: Vec<bool>,     /* array with active/inactive event functions */
    pub cv_mxgnull: i32,           /* num. warning messages about possible g==0 */

    /*---------------------------
      Inequality Constraints Data
      ---------------------------*/
    pub cv_constraints: NVector, /* vector of constraint flags        */
    pub constraint_corrections: i64, /* total constraint corrections  */
    pub constraint_fails: i64,   /* total constraint failures         */
    pub max_constraint_fails: i32, /* max failures allowed in a step  */

    /*---------------
      Projection Data
      ---------------*/
    /* PLACEHOLDER (phased port): becomes Option<Box<CVodeProjMem>>
       when cvodes_proj_impl.rs lands (mirroring the donor's
       cvode_proj_impl.rs); until then an opaque box keeps the field
       and the C name. */
    pub proj_mem: Option<Box<crate::cvodes_proj_impl::CVodeProjMem>>, /* projection memory structure  */
    pub proj_enabled: bool,   /* flag indicating if projection is enabled  */
    pub proj_applied: bool,   /* flag indicating if projection was applied */
    pub proj_p: [f64; L_MAX], /* coefficients of p(x) (degree q poly)      */

    /* Fused Vector Operations: cv_cvals / cv_Xvecs / cv_Zvecs omitted
       per donor precedent -- fused ops are expanded to plain loops in
       the pure-Rust build, so no scalar/vector scratch arrays exist. */

    /*----------------
      Resizing History
      ----------------*/
    pub first_step_after_resize: bool, /* flag to signal a resize happened */

    /*------------------------
      Adjoint sensitivity data
      ------------------------*/
    pub cv_adj: bool, /* SUNTRUE if performing ASA                */

    pub cv_adj_mem: Option<Box<CVadjMem>>, /* adjoint memory structure  */

    pub cv_adjMallocDone: bool,
}

/*
 * =================================================================
 *   A D J O I N T   M O D U L E    M E M O R Y    B L O C K
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * CVckpntMem (struct CVckpntMemRec)
 * -----------------------------------------------------------------
 * All information at a check point needed to 'hot' start cvodes.
 * The C ck_next pointer is dropped: CVadjMem owns the check points
 * as Vec<CVckpntMem> (pinned decision 4).
 * -----------------------------------------------------------------
 */

pub struct CVckpntMem {
    /* Integration limits */
    pub ck_t0: f64,
    pub ck_t1: f64,

    /* Nordsieck History Array */
    pub ck_zn: Vec<NVector>,

    /* Do we need to carry quadratures? */
    pub ck_quadr: bool,

    /* Nordsieck History Array for quadratures */
    pub ck_znQ: Vec<NVector>,

    /* Do we need to carry sensitivities? */
    pub ck_sensi: bool,

    /* number of sensitivities */
    pub ck_Ns: i32,

    /* Nordsieck History Array for sensitivities */
    pub ck_znS: [Vec<NVector>; L_MAX],

    /* Do we need to carry quadrature sensitivities? */
    pub ck_quadr_sensi: bool,

    /* Nordsieck History Array for quadrature sensitivities */
    pub ck_znQS: [Vec<NVector>; L_MAX],

    /* Was ck_zn[qmax] allocated?
       ck_zqm = 0    - no
       ck_zqm = qmax - yes      */
    pub ck_zqm: i32,

    /* Step data */
    pub ck_nst: i64,
    pub ck_tretlast: f64,
    pub ck_q: i32,
    pub ck_qprime: i32,
    pub ck_qwait: i32,
    pub ck_L: i32,
    pub ck_gammap: f64,
    pub ck_h: f64,
    pub ck_hprime: f64,
    pub ck_hscale: f64,
    pub ck_eta: f64,
    pub ck_etamax: f64,
    pub ck_tau: [f64; L_MAX + 1],
    pub ck_tq: [f64; NUM_TESTS + 1],
    pub ck_l: [f64; L_MAX],

    /* Saved values */
    pub ck_saved_tq5: f64,
    /* ck_next dropped (CVadjMem holds Vec<CVckpntMem>) */
}

/*
 * -----------------------------------------------------------------
 * CVdtpntMem (struct CVdtpntMemRec)
 * -----------------------------------------------------------------
 * All information at a data point needed to interpolate the
 * solution of forward simulations. The C `void* content` (whose
 * layout depends on IMtype: CVhermiteDataMemRec or
 * CVpolynomialDataMemRec) becomes the DtpntContent enum, and the C
 * interpolation-module function pointers cvaIMMallocFn/cvaIMFreeFn/
 * cvaIMStorePntFn/cvaIMGetYFn are replaced by dispatch on
 * ca_IMtype in cvodea.rs (pinned decision 4).
 * -----------------------------------------------------------------
 */

pub enum DtpntContent {
    None,
    /* Data for cubic Hermite interpolation (CVhermiteDataMemRec) */
    Hermite {
        y: NVector,
        yd: NVector,
        yS: Vec<NVector>,
        ySd: Vec<NVector>,
    },
    /* Data for polynomial interpolation (CVpolynomialDataMemRec) */
    Polynomial {
        y: NVector,
        yS: Vec<NVector>,
        order: i32,
    },
}

pub struct CVdtpntMem {
    pub t: f64,                /* time */
    pub content: DtpntContent, /* IMtype-dependent content */
}

/*
 * -----------------------------------------------------------------
 * CVodeBMem (struct CVodeBMemRec)
 * -----------------------------------------------------------------
 * All information for ONE backward problem. The C cv_next pointer
 * is dropped: CVadjMem owns the backward problems as
 * Vec<CVodeBMem>. The C cv_lfree/cv_pfree function pointers are
 * dropped (Rust Drop frees the attached memory). cv_lmem/cv_pmem
 * are opaque placeholders: backward problems attach CVLS through
 * the nested cv_mem, and the cvodea.c port will refine these
 * (pinned decision 4).
 * -----------------------------------------------------------------
 */

pub struct CVodeBMem {
    /* Index of this backward problem */
    pub cv_index: i32,

    /* Time at which the backward problem is initialized */
    pub cv_t0: f64,

    /* CVODES memory for this backward problem */
    pub cv_mem: Box<CVodeMem>,

    /* Flags to indicate that this backward problem's RHS or quad RHS
     * require forward sensitivities */
    pub cv_f_withSensi: bool,
    pub cv_fQ_withSensi: bool,

    /* Right hand side function for backward run */
    pub cv_f: Option<CVRhsFnB>,
    pub cv_fs: Option<CVRhsFnBS>,

    /* Right hand side quadrature function for backward run */
    pub cv_fQ: Option<CVQuadRhsFnB>,
    pub cv_fQs: Option<CVQuadRhsFnBS>,

    /* User user_data */
    pub cv_user_data: UserData,

    /* Memory block for a linear solver's interface to CVODEA */
    pub cv_lmem: Option<Box<dyn std::any::Any>>,

    /* Memory block for a preconditioner's module interface to CVODEA */
    pub cv_pmem: Option<Box<dyn std::any::Any>>,

    /* Time at which to extract solution / quadratures */
    pub cv_tout: f64,

    /* Workspace Nvector */
    pub cv_y: NVector,
    /* cv_next dropped (CVadjMem holds Vec<CVodeBMem>) */
}

/*
 * -----------------------------------------------------------------
 * CVadjMem (struct CVadjMemRec)
 * -----------------------------------------------------------------
 * All information necessary for adjoint sensitivity analysis.
 * Linked lists become owned Vecs; the C pointers ca_bckpbCrt and
 * ca_ckpntData become Option<usize> indices into cvB_mem / ck_mem;
 * the IM function pointers are replaced by dispatch on ca_IMtype
 * (pinned decision 4).
 * -----------------------------------------------------------------
 */

pub struct CVadjMem {
    /* --------------------
     * Forward problem data
     * -------------------- */

    /* Integration interval */
    pub ca_tinitial: f64,
    pub ca_tfinal: f64,

    /* Flag for first call to CVodeF */
    pub ca_firstCVodeFcall: bool,

    /* Flag if CVodeF was called with TSTOP */
    pub ca_tstopCVodeFcall: bool,
    pub ca_tstopCVodeF: f64,

    /* Flag if CVodeF was called in CV_NORMAL_MODE and encountered a
       root after tout */
    pub ca_rootret: bool,
    pub ca_troot: f64,

    /* ----------------------
     * Backward problems data
     * ---------------------- */

    /* Storage for backward problems (C: linked list headed at cvB_mem) */
    pub cvB_mem: Vec<CVodeBMem>,

    /* Number of backward problems */
    pub ca_nbckpbs: i32,

    /* Current backward problem (C: pointer; here index into cvB_mem) */
    pub ca_bckpbCrt: Option<usize>,

    /* Flag for first call to CVodeB */
    pub ca_firstCVodeBcall: bool,

    /* ----------------
     * Check point data
     * ---------------- */

    /* Storage for check point information (C: linked list headed at ck_mem) */
    pub ck_mem: Vec<CVckpntMem>,

    /* Number of check points */
    pub ca_nckpnts: i32,

    /* Check point for which data is available (C: pointer; here index
       into ck_mem) */
    pub ca_ckpntData: Option<usize>,

    /* ------------------
     * Interpolation data
     * ------------------ */

    /* Number of steps between 2 check points */
    pub ca_nsteps: i64,

    /* Last index used in CVAfindIndex */
    pub ca_ilast: i64,

    /* Storage for data from forward runs */
    pub dt_mem: Vec<CVdtpntMem>,

    /* Actual number of data points in dt_mem (typically np=nsteps+1) */
    pub ca_np: i64,

    /* Interpolation type (CV_HERMITE or CV_POLYNOMIAL) */
    pub ca_IMtype: i32,

    /* The C interpolation-module function pointers ca_IMmalloc,
       ca_IMfree, ca_IMstore and ca_IMget are replaced by dispatch on
       ca_IMtype: cvodea.rs implements cvaHermiteMalloc/Free/StorePnt/
       GetY and cvaPolynomialMalloc/Free/StorePnt/GetY and selects the
       family by ca_IMtype at each call site (pinned decision 4). */

    /* Flags controlling the interpolation module */
    pub ca_IMmallocDone: bool,  /* IM initialized? */
    pub ca_IMnewData: bool,     /* new data available in dt_mem?*/
    pub ca_IMstoreSensi: bool,  /* store sensitivities? */
    pub ca_IMinterpSensi: bool, /* interpolate sensitivities? */

    /* Workspace for the interpolation module.
       In C, ca_Y[i]/ca_YS[i] are POINTERS into zn[i]/znS[i]; here
       they are owned scratch storage (empty = C NULL) and the
       cvodea.rs port copies data instead of aliasing (pinned
       decision 4). */
    pub ca_Y: Vec<NVector>, /* owned scratch for zn[i]  (C: N_Vector ca_Y[L_MAX])   */
    pub ca_YS: [Vec<NVector>; L_MAX], /* owned scratch for znS[i] (C: N_Vector* ca_YS[L_MAX]) */
    pub ca_T: [f64; L_MAX],

    /* -------------------------------
     * Workspace for wrapper functions
     * ------------------------------- */
    pub ca_ytmp: NVector,

    pub ca_yStmp: Vec<NVector>,
}

/* ===============================================================
   Error handler: the C cvProcessError routes printf-style messages
   to the SUNContext error handler stack; here messages go to stderr
   (equivalent to the default SUNLogErrHandlerFn behavior). Same
   adaptation as the donor's cvode_impl.rs.
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
        eprintln!("\n[CVODES WARNING] {file}:{line} in {func}\n  {msg}\n");
    } else {
        eprintln!("\n[CVODES ERROR] {file}:{line} in {func}\n  {msg}\n");
    }
}

/* Error messages (cvodes_impl.h). Messages that embed printf format
   directives (MSG_TIME, MSG_TIME_H, MSG_TIME_INT, MSG_TIME_TOUT,
   MSG_TIME_TSTOP, %d) are built inline with format! at the call
   sites, as in the donor's cvode.rs. */

/* Initialization and I/O error messages */
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
pub const MSGCV_BAD_ISM_CONSTR: &str =
    "Constraints can not be enforced while forward sensitivity is used with \
     simultaneous method";
pub const MSGCV_NULL_F: &str = "f = NULL illegal.";
pub const MSGCV_NULL_G: &str = "g = NULL illegal.";
pub const MSGCV_BAD_NVECTOR: &str = "A required vector operation is not implemented.";
pub const MSGCV_BAD_CONSTR: &str = "Illegal values in constraints vector.";
pub const MSGCV_BAD_K: &str = "Illegal value for k.";
pub const MSGCV_NULL_DKY: &str = "dky = NULL illegal.";
pub const MSGCV_NO_ROOT: &str = "Rootfinding was not initialized.";
pub const MSGCV_NLS_INIT_FAIL: &str = "The nonlinear solver's init routine failed.";

pub const MSGCV_NO_QUAD: &str = "Quadrature integration not activated.";
pub const MSGCV_BAD_ITOLQ: &str =
    "Illegal value for itolQ. The legal values are CV_SS and CV_SV.";
pub const MSGCV_NULL_ABSTOLQ: &str = "abstolQ = NULL illegal.";
pub const MSGCV_BAD_RELTOLQ: &str = "reltolQ < 0 illegal.";
pub const MSGCV_BAD_ABSTOLQ: &str = "abstolQ has negative component(s) (illegal).";

pub const MSGCV_SENSINIT_2: &str = "Sensitivity analysis already initialized.";
pub const MSGCV_NO_SENSI: &str = "Forward sensitivity analysis not activated.";
pub const MSGCV_BAD_ITOLS: &str =
    "Illegal value for itolS. The legal values are CV_SS, CV_SV, and CV_EE.";
pub const MSGCV_NULL_ABSTOLS: &str = "abstolS = NULL illegal.";
pub const MSGCV_BAD_RELTOLS: &str = "reltolS < 0 illegal.";
pub const MSGCV_BAD_ABSTOLS: &str = "abstolS has negative component(s) (illegal).";
pub const MSGCV_BAD_PBAR: &str = "pbar has zero component(s) (illegal).";
pub const MSGCV_BAD_PLIST: &str = "plist has negative component(s) (illegal).";
pub const MSGCV_BAD_NS: &str = "NS <= 0 illegal.";
pub const MSGCV_NULL_YS0: &str = "yS0 = NULL illegal.";
pub const MSGCV_BAD_ISM: &str =
    "Illegal value for ism. Legal values are: CV_SIMULTANEOUS, CV_STAGGERED \
     and CV_STAGGERED1.";
pub const MSGCV_BAD_IFS: &str =
    "Illegal value for ifS. Legal values are: CV_ALLSENS and CV_ONESENS.";
pub const MSGCV_BAD_ISM_IFS: &str = "Illegal ism = CV_STAGGERED1 for CVodeSensInit.";
pub const MSGCV_BAD_IS: &str = "Illegal value for is.";
pub const MSGCV_NULL_DKYA: &str = "dkyA = NULL illegal.";
pub const MSGCV_BAD_DQTYPE: &str =
    "Illegal value for DQtype. Legal values are: CV_CENTERED and CV_FORWARD.";
pub const MSGCV_BAD_DQRHO: &str = "DQrhomax < 0 illegal.";

pub const MSGCV_BAD_ITOLQS: &str =
    "Illegal value for itolQS. The legal values are CV_SS, CV_SV, and CV_EE.";
pub const MSGCV_NULL_ABSTOLQS: &str = "abstolQS = NULL illegal.";
pub const MSGCV_BAD_RELTOLQS: &str = "reltolQS < 0 illegal.";
pub const MSGCV_BAD_ABSTOLQS: &str = "abstolQS has negative component(s) (illegal).";
pub const MSGCV_NO_QUADSENSI: &str =
    "Forward sensitivity analysis for quadrature variables not activated.";
pub const MSGCV_NULL_YQS0: &str = "yQS0 = NULL illegal.";

/* CVode Error Messages */
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

/* CVode Projection Error Messages */
pub const MSG_CV_MEM_NULL: &str = "cvode_mem = NULL illegal.";
pub const MSG_CV_MEM_FAIL: &str = "A memory request failed.";
pub const MSG_CV_PROJ_MEM_NULL: &str = "proj_mem = NULL illegal.";

pub const MSGCV_NO_TOLQ: &str =
    "No integration tolerances for quadrature variables have been specified.";
pub const MSGCV_BAD_EWTQ: &str = "Initial ewtQ has component(s) equal to zero (illegal).";
pub const MSGCV_QRHSFUNC_FIRST: &str =
    "The quadrature right-hand side routine failed at the first call.";

pub const MSGCV_NO_TOLS: &str =
    "No integration tolerances for sensitivity variables have been specified.";
pub const MSGCV_NULL_P: &str = "p = NULL when using internal DQ for sensitivity RHS illegal.";
pub const MSGCV_BAD_EWTS: &str = "Initial ewtS has component(s) equal to zero (illegal).";
pub const MSGCV_SRHSFUNC_FIRST: &str =
    "The sensitivity right-hand side routine failed at the first call.";

pub const MSGCV_NULL_FQ: &str =
    "CVODES is expected to use DQ to evaluate the RHS of quad. sensi., but \
     quadratures were not initialized.";
pub const MSGCV_NO_TOLQS: &str =
    "No integration tolerances for quadrature sensitivity variables have been \
     specified.";
pub const MSGCV_BAD_EWTQS: &str = "Initial ewtQS has component(s) equal to zero (illegal).";
pub const MSGCV_QSRHSFUNC_FIRST: &str =
    "The quadrature sensitivity right-hand side routine failed at the first \
     call.";

/* Adjoint Error Messages */
pub const MSGCV_NO_ADJ: &str = "Illegal attempt to call before calling CVodeAdjMalloc.";
pub const MSGCV_BAD_STEPS: &str = "Steps nonpositive illegal.";
pub const MSGCV_BAD_INTERP: &str = "Illegal value for interp.";
pub const MSGCV_BAD_WHICH: &str = "Illegal value for which.";
pub const MSGCV_NO_BCK: &str = "No backward problems have been defined yet.";
pub const MSGCV_NO_FWD: &str = "Illegal attempt to call before calling CVodeF.";
pub const MSGCV_BAD_SENSI: &str =
    "At least one backward problem requires sensitivities, but they were not \
     stored for interpolation.";
pub const MSGCV_BAD_ITASKB: &str =
    "Illegal value for itaskB. Legal values are CV_NORMAL and CV_ONE_STEP.";
pub const MSGCV_BAD_TBOUT: &str =
    "The final time tBout is outside the interval over which the forward \
     problem was solved.";
pub const MSGCV_WRONG_INTERP: &str =
    "This function cannot be called for the specified interp type.";

#[cfg(test)]
mod tests {
    use super::*;

    /* The donor cvode_impl.rs carries no tests; this minimal check
       pins the C constant values the CVODES modules rely on
       (including the ones whose numeric value differs from CVODE:
       CV_CONTEXT_ERR and the projection return codes). */
    #[test]
    fn c_constant_values() {
        assert_eq!(L_MAX, 13);
        assert_eq!(NUM_TESTS, 5);
        assert_eq!(CV_SIMULTANEOUS, 1);
        assert_eq!(CV_STAGGERED, 2);
        assert_eq!(CV_STAGGERED1, 3);
        assert_eq!(CV_HERMITE, 1);
        assert_eq!(CV_POLYNOMIAL, 2);
        assert_eq!(CV_NO_QUAD, -30);
        assert_eq!(CV_UNREC_QSRHSFUNC_ERR, -54);
        assert_eq!(CV_CONTEXT_ERR, -55);
        assert_eq!(CV_PROJ_MEM_NULL, -56);
        assert_eq!(CV_BAD_TINTERP, -59);
        assert_eq!(CV_NO_ADJ, -101);
        assert_eq!(CV_GETY_BADT, -107);
        assert_eq!(QSRHSFUNC_RECVR, 15);
        assert_eq!(NLS_MAXCOR, 3);
        assert_eq!(CRDOWN, 0.3);
        assert_eq!(RDIV, 2.0);
    }
}
