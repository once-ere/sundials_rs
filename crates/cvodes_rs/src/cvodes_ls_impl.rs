/* -----------------------------------------------------------------
 * Translated from src/cvodes/cvodes_ls_impl.h and the type
 * definitions of include/cvodes/cvodes_ls.h (CVODES 7.7.0).
 * CVLS linear-solver-interface memory structures (forward problem
 * CVLsMemRec and backward problem CVLsMemRecB), constants and error
 * messages. Follows the verified CVODE donor conventions
 * (crates/cvode_rs/src/cvode_ls_impl.rs) field-for-field:
 *  - the void* data pointers (J_data, P_data, jt_data, A_data,
 *    P_dataB) are dropped: internal routines receive &mut CVodeMem
 *    directly and user routines receive cv_user_data at the call
 *    sites in cvodes_ls.rs;
 *  - ycur/fcur are not stored: in C they are aliases of the
 *    integrator's cv_y/cv_ftemp set by cvLsSetup/cvLsSolve, and the
 *    cvodes_ls.rs routines read those CVodeMem fields directly;
 *  - pfree is dropped (RAII): the preconditioner attachment is the
 *    PrecModule enum below.
 * -----------------------------------------------------------------*/
use crate::cvodes_impl::CVRhsFn;
use crate::nvector_serial::NVector;
use crate::sundials_linearsolver::LinearSolver;
use crate::sundials_matrix::SUNMatrix;
use crate::sundials_types::UserData;

/* CVLS return codes (cvodes_ls.h) */
pub const CVLS_SUCCESS: i32 = 0;
pub const CVLS_MEM_NULL: i32 = -1;
pub const CVLS_LMEM_NULL: i32 = -2;
pub const CVLS_ILL_INPUT: i32 = -3;
pub const CVLS_MEM_FAIL: i32 = -4;
pub const CVLS_PMEM_NULL: i32 = -5;
pub const CVLS_JACFUNC_UNRECVR: i32 = -6;
pub const CVLS_JACFUNC_RECVR: i32 = -7;
pub const CVLS_SUNMAT_FAIL: i32 = -8;
pub const CVLS_SUNLS_FAIL: i32 = -9;

/* Return values for the adjoint module (cvodes_ls.h) */
pub const CVLS_NO_ADJ: i32 = -101;
pub const CVLS_LMEMB_NULL: i32 = -102;

/* CVSLS solver constants (cvodes_ls_impl.h) */
pub const CVLS_MSBJ: i64 = 51;
pub const CVLS_DGMAX: f64 = 0.2;
pub const CVLS_EPLIN: f64 = 0.05;

/* =================================================================
   PART I:  Forward Problems
   =================================================================*/

/* User-supplied function types (cvodes_ls.h) */
pub type CVLsJacFn = fn(
    t: f64,
    y: &NVector,
    fy: &NVector,
    jac: &mut SUNMatrix,
    user_data: &mut UserData,
    tmp1: &mut NVector,
    tmp2: &mut NVector,
    tmp3: &mut NVector,
) -> i32;

pub type CVLsPrecSetupFn = fn(
    t: f64,
    y: &NVector,
    fy: &NVector,
    jok: bool,
    jcur_ptr: &mut bool,
    gamma: f64,
    user_data: &mut UserData,
) -> i32;

pub type CVLsPrecSolveFn = fn(
    t: f64,
    y: &NVector,
    fy: &NVector,
    r: &NVector,
    z: &mut NVector,
    gamma: f64,
    delta: f64,
    lr: i32,
    user_data: &mut UserData,
) -> i32;

pub type CVLsJacTimesSetupFn =
    fn(t: f64, y: &NVector, fy: &NVector, user_data: &mut UserData) -> i32;

pub type CVLsJacTimesVecFn = fn(
    v: &NVector,
    jv: &mut NVector,
    t: f64,
    y: &NVector,
    fy: &NVector,
    user_data: &mut UserData,
    tmp: &mut NVector,
) -> i32;

pub type CVLsLinSysFn = fn(
    t: f64,
    y: &NVector,
    fy: &NVector,
    a: &mut SUNMatrix,
    jok: bool,
    jcur: &mut bool,
    gamma: f64,
    user_data: &mut UserData,
    tmp1: &mut NVector,
    tmp2: &mut NVector,
    tmp3: &mut NVector,
) -> i32;

/* Preconditioner module attached to CVLS.
   In C this is the P_data/pfree convention: P_data points at
   user_data (user-supplied pset/psolve) or at an internal
   preconditioner module (cvodes_bandpre / cvodes_bbdpre).

   NOTE (phased port, mirroring LsModule in cvodes_impl.rs): the
   donor enum carries BandPre(Box<CVBandPrecData>) and
   BBDPre(Box<CVBBDPrecData>) variants; those are added here when
   cvodes_bandpre_impl.rs / cvodes_bbdpre_impl.rs land. Until then
   only the user-supplied/none states exist. */
#[derive(Default)]
pub enum PrecModule {
    #[default]
    None,
    /// user-supplied pset/psolve get user_data
    User,
}

/* -----------------------------------------------------------------
   Types : CVLsMemRec, CVLsMem
   -----------------------------------------------------------------*/
pub struct CVLsMem {
    /* Linear solver type information */
    pub iterative: bool,   /* is the solver iterative?    */
    pub matrixbased: bool, /* is a matrix structure used? */

    /* Jacobian construction & storage */
    pub jacDQ: bool,            /* SUNTRUE if using internal DQ Jacobian approx. */
    pub jac: Option<CVLsJacFn>, /* Jacobian routine to be called                 */
    pub jbad: bool,             /* heuristic suggestion for pset                 */
    pub dgmax_jbad: f64,        /* if convfail = FAIL_BAD_J and the gamma ratio *
                                 * |gamma/gammap-1| < dgmax_jbad then J is bad  */

    /* Matrix-based solver, scale solution to account for change in gamma */
    pub scalesol: bool,

    /* Iterative solver tolerance */
    pub eplifac: f64, /* nonlinear -> linear tol scaling factor  */
    pub nrmfac: f64,  /* integrator -> LS norm conversion factor */

    /* Linear solver, matrix and vector objects/pointers */
    pub LS: LinearSolver,          /* generic linear solver object     */
    pub A: Option<SUNMatrix>,      /* A = I - gamma * df/dy            */
    pub savedJ: Option<SUNMatrix>, /* savedJ = old Jacobian            */
    pub ytemp: NVector,            /* temp vector for jtimes & psolve  */
    pub x: NVector,                /* temp vector used by CVLsSolve    */
    /* (C also aliases ycur/fcur into cv_mem; passed as args here)   */

    /* Statistics and associated parameters */
    pub msbj: i64,     /* max num steps between jac/pset calls */
    pub nje: i64,      /* nje = no. of calls to jac            */
    pub nfeDQ: i64,    /* no. of calls to f due to DQ Jacobian
                       or J*v approximations                   */
    pub nstlj: i64,    /* nstlj = nst at last jac/pset call    */
    pub npe: i64,      /* npe = total number of pset calls     */
    pub nli: i64,      /* nli = total number of linear iters   */
    pub nps: i64,      /* nps = total number of psolve calls   */
    pub ncfl: i64,     /* ncfl = total number of conv. failures*/
    pub njtsetup: i64, /* total number of calls to jtsetup     */
    pub njtimes: i64,  /* total number of calls to jtimes      */
    pub tnlj: f64,     /* tnlj = t_n at last jac/pset call     */

    /* Preconditioner computation */
    pub pset: Option<CVLsPrecSetupFn>,
    pub psolve: Option<CVLsPrecSolveFn>,
    pub prec_module: PrecModule,

    /* Jacobian times vector computation */
    pub jtimesDQ: bool,
    pub jtsetup: Option<CVLsJacTimesSetupFn>,
    pub jtimes: Option<CVLsJacTimesVecFn>,
    pub jt_f: Option<CVRhsFn>,

    /* Linear system setup function */
    pub user_linsys: bool,
    pub linsys: Option<CVLsLinSysFn>,

    /* In C, cvLsInitialize NULLs cv_mem->cv_lsetup when no setup work
       exists (matrix-free without preconditioner, or matrix-embedded);
       this flag carries that decision into the dispatch. */
    pub setup_disabled: bool,

    pub last_flag: i32, /* last error flag returned by any function */
}

/* =================================================================
   PART II:  Backward Problems
   =================================================================*/

/* User-supplied function types for backward problems (cvodes_ls.h).
   The C `N_Vector* yS` argument becomes `&[NVector]`, matching the
   CVRhsFnBS convention pinned in cvodes_impl.rs. */

pub type CVLsJacFnB = fn(
    t: f64,
    y: &NVector,
    yB: &NVector,
    fyB: &NVector,
    jacB: &mut SUNMatrix,
    user_dataB: &mut UserData,
    tmp1B: &mut NVector,
    tmp2B: &mut NVector,
    tmp3B: &mut NVector,
) -> i32;

pub type CVLsJacFnBS = fn(
    t: f64,
    y: &NVector,
    yS: &[NVector],
    yB: &NVector,
    fyB: &NVector,
    jacB: &mut SUNMatrix,
    user_dataB: &mut UserData,
    tmp1B: &mut NVector,
    tmp2B: &mut NVector,
    tmp3B: &mut NVector,
) -> i32;

pub type CVLsPrecSetupFnB = fn(
    t: f64,
    y: &NVector,
    yB: &NVector,
    fyB: &NVector,
    jokB: bool,
    jcurPtrB: &mut bool,
    gammaB: f64,
    user_dataB: &mut UserData,
) -> i32;

pub type CVLsPrecSetupFnBS = fn(
    t: f64,
    y: &NVector,
    yS: &[NVector],
    yB: &NVector,
    fyB: &NVector,
    jokB: bool,
    jcurPtrB: &mut bool,
    gammaB: f64,
    user_dataB: &mut UserData,
) -> i32;

pub type CVLsPrecSolveFnB = fn(
    t: f64,
    y: &NVector,
    yB: &NVector,
    fyB: &NVector,
    rB: &NVector,
    zB: &mut NVector,
    gammaB: f64,
    deltaB: f64,
    lrB: i32,
    user_dataB: &mut UserData,
) -> i32;

pub type CVLsPrecSolveFnBS = fn(
    t: f64,
    y: &NVector,
    yS: &[NVector],
    yB: &NVector,
    fyB: &NVector,
    rB: &NVector,
    zB: &mut NVector,
    gammaB: f64,
    deltaB: f64,
    lrB: i32,
    user_dataB: &mut UserData,
) -> i32;

pub type CVLsJacTimesSetupFnB =
    fn(t: f64, y: &NVector, yB: &NVector, fyB: &NVector, user_dataB: &mut UserData) -> i32;

pub type CVLsJacTimesSetupFnBS = fn(
    t: f64,
    y: &NVector,
    yS: &[NVector],
    yB: &NVector,
    fyB: &NVector,
    user_dataB: &mut UserData,
) -> i32;

pub type CVLsJacTimesVecFnB = fn(
    vB: &NVector,
    JvB: &mut NVector,
    t: f64,
    y: &NVector,
    yB: &NVector,
    fyB: &NVector,
    user_dataB: &mut UserData,
    tmpB: &mut NVector,
) -> i32;

pub type CVLsJacTimesVecFnBS = fn(
    vB: &NVector,
    JvB: &mut NVector,
    t: f64,
    y: &NVector,
    yS: &[NVector],
    yB: &NVector,
    fyB: &NVector,
    user_dataB: &mut UserData,
    tmpB: &mut NVector,
) -> i32;

pub type CVLsLinSysFnB = fn(
    t: f64,
    y: &NVector,
    yB: &NVector,
    fyB: &NVector,
    AB: &mut SUNMatrix,
    jokB: bool,
    jcurB: &mut bool,
    gammaB: f64,
    user_dataB: &mut UserData,
    tmp1B: &mut NVector,
    tmp2B: &mut NVector,
    tmp3B: &mut NVector,
) -> i32;

pub type CVLsLinSysFnBS = fn(
    t: f64,
    y: &NVector,
    yS: &[NVector],
    yB: &NVector,
    fyB: &NVector,
    AB: &mut SUNMatrix,
    jokB: bool,
    jcurB: &mut bool,
    gammaB: f64,
    user_dataB: &mut UserData,
    tmp1B: &mut NVector,
    tmp2B: &mut NVector,
    tmp3B: &mut NVector,
) -> i32;

/* -----------------------------------------------------------------
   Types : CVLsMemRecB, CVLsMemB

   CVodeSetLinearSolverB attaches such a structure to the lmemB
   field of CVodeBMem (here: boxed into CVodeBMem.cv_lmem, the
   Option<Box<dyn Any>> placeholder pinned in cvodes_impl.rs;
   cvodes_ls.rs downcasts it back to CVLsMemB).
   The C void* P_dataB is dropped: it is set to NULL on creation and
   never read anywhere in the CVODES sources.
   -----------------------------------------------------------------*/
#[derive(Default)]
pub struct CVLsMemB {
    pub jacB: Option<CVLsJacFnB>,
    pub jacBS: Option<CVLsJacFnBS>,
    pub jtsetupB: Option<CVLsJacTimesSetupFnB>,
    pub jtsetupBS: Option<CVLsJacTimesSetupFnBS>,
    pub jtimesB: Option<CVLsJacTimesVecFnB>,
    pub jtimesBS: Option<CVLsJacTimesVecFnBS>,
    pub linsysB: Option<CVLsLinSysFnB>,
    pub linsysBS: Option<CVLsLinSysFnBS>,
    pub psetB: Option<CVLsPrecSetupFnB>,
    pub psetBS: Option<CVLsPrecSetupFnBS>,
    pub psolveB: Option<CVLsPrecSolveFnB>,
    pub psolveBS: Option<CVLsPrecSolveFnBS>,
}

/* =================================================================
   Error Messages (cvodes_ls_impl.h)
   =================================================================*/

pub const MSG_LS_CVMEM_NULL: &str = "Integrator memory is NULL.";
pub const MSG_LS_MEM_FAIL: &str = "A memory request failed.";
pub const MSG_LS_BAD_NVECTOR: &str = "A required vector operation is not implemented.";
pub const MSG_LS_BAD_LSTYPE: &str = "Incompatible linear solver type.";
pub const MSG_LS_LMEM_NULL: &str = "Linear solver memory is NULL.";
pub const MSG_LS_BAD_SIZES: &str =
    "Illegal bandwidth parameter(s). Must have 0 <=  ml, mu <= N-1.";
pub const MSG_LS_BAD_EPLIN: &str = "eplifac < 0 illegal.";
pub const MSG_LS_BAD_PRETYPE: &str = "Illegal value for pretype. Legal values are PREC_NONE, \
                                      PREC_LEFT, PREC_RIGHT, and PREC_BOTH.";
pub const MSG_LS_PSOLVE_REQ: &str = "pretype != PREC_NONE, but PSOLVE = NULL is illegal.";
pub const MSG_LS_BAD_GSTYPE: &str =
    "Illegal value for gstype. Legal values are MODIFIED_GS and CLASSICAL_GS.";

pub const MSG_LS_PSET_FAILED: &str =
    "The preconditioner setup routine failed in an unrecoverable manner.";
pub const MSG_LS_PSOLVE_FAILED: &str =
    "The preconditioner solve routine failed in an unrecoverable manner.";
pub const MSG_LS_JTSETUP_FAILED: &str =
    "The Jacobian x vector setup routine failed in an unrecoverable manner.";
pub const MSG_LS_JTIMES_FAILED: &str =
    "The Jacobian x vector routine failed in an unrecoverable manner.";
pub const MSG_LS_JACFUNC_FAILED: &str =
    "The Jacobian routine failed in an unrecoverable manner.";
pub const MSG_LS_SUNMAT_FAILED: &str =
    "A SUNMatrix routine failed in an unrecoverable manner.";

pub const MSG_LS_NO_ADJ: &str = "Illegal attempt to call before calling CVodeAdjMalloc.";
pub const MSG_LS_BAD_WHICH: &str = "Illegal value for which.";
pub const MSG_LS_LMEMB_NULL: &str = "Linear solver memory is NULL for the backward integration.";
pub const MSG_LS_BAD_TINTERP: &str = "Bad t for interpolation.";
