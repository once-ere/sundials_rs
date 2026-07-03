/* -----------------------------------------------------------------
 * Translated from src/cvode/cvode_ls_impl.h and the type definitions
 * of include/cvode/cvode_ls.h (CVODE 7.7.0).
 * CVLS linear-solver-interface memory structure and constants.
 * -----------------------------------------------------------------*/
use crate::cvode_bandpre_impl::CVBandPrecData;
use crate::cvode_bbdpre_impl::CVBBDPrecData;
use crate::nvector_serial::NVector;
use crate::sundials_linearsolver::LinearSolver;
use crate::sundials_matrix::SUNMatrix;
use crate::sundials_types::UserData;

/* CVLS return codes (cvode_ls.h) */
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

/* CVLS solver constants (cvode_ls_impl.h) */
pub const CVLS_MSBJ: i64 = 51;
pub const CVLS_DGMAX: f64 = 0.2;
pub const CVLS_EPLIN: f64 = 0.05;

/* User-supplied function types (cvode_ls.h) */
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
   preconditioner module (cvode_bandpre / cvode_bbdpre). */
#[derive(Default)]
pub enum PrecModule {
    #[default]
    None,
    /// user-supplied pset/psolve get user_data
    User,
    /// CVBandPrecInit module data
    BandPre(Box<CVBandPrecData>),
    /// CVBBDPrecInit module data
    BBDPre(Box<CVBBDPrecData>),
}

/* -----------------------------------------------------------------
   Types : CVLsMemRec, CVLsMem
   -----------------------------------------------------------------*/
pub struct CVLsMem {
    /* Linear solver type information */
    pub iterative: bool,   /* is the solver iterative?    */
    pub matrixbased: bool, /* is a matrix structure used? */

    /* Jacobian construction & storage */
    pub jacDQ: bool,           /* SUNTRUE if using internal DQ Jacobian approx. */
    pub jac: Option<CVLsJacFn>, /* Jacobian routine to be called                */
    pub jbad: bool,            /* heuristic suggestion for pset                 */
    pub dgmax_jbad: f64,       /* |gamma/gammap-1| < dgmax_jbad => J not bad    */

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
    pub nje: i64,      /* number of calls to jac               */
    pub nfeDQ: i64,    /* number of f calls for DQ Jacobian/Jv */
    pub nstlj: i64,    /* nst at last jac/pset call            */
    pub npe: i64,      /* total number of pset calls           */
    pub nli: i64,      /* total number of linear iterations    */
    pub nps: i64,      /* total number of psolve calls         */
    pub ncfl: i64,     /* total number of convergence failures */
    pub njtsetup: i64, /* total number of calls to jtsetup     */
    pub njtimes: i64,  /* total number of calls to jtimes      */
    pub tnlj: f64,     /* t_n at last jac/pset call            */

    /* Preconditioner computation */
    pub pset: Option<CVLsPrecSetupFn>,
    pub psolve: Option<CVLsPrecSolveFn>,
    pub prec_module: PrecModule,

    /* Jacobian times vector computation */
    pub jtimesDQ: bool,
    pub jtsetup: Option<CVLsJacTimesSetupFn>,
    pub jtimes: Option<CVLsJacTimesVecFn>,
    pub jt_f: Option<crate::cvode_impl::CVRhsFn>,

    /* Linear system setup function */
    pub user_linsys: bool,
    pub linsys: Option<CVLsLinSysFn>,

    /* In C, cvLsInitialize NULLs cv_mem->cv_lsetup when no setup work
       exists (matrix-free without preconditioner, or matrix-embedded);
       this flag carries that decision into the dispatch. */
    pub setup_disabled: bool,

    pub last_flag: i32, /* last error flag returned by any function */
}
