/* -----------------------------------------------------------------
 * Translated from src/arkode/arkode_ls_impl.h and the type
 * definitions of include/arkode/arkode_ls.h (ARKODE 7.7.0).
 * ARKLS linear-solver-interface memory structure and constants.
 *
 * Storage adaptation (ARCHITECTURE.md Addendum C.1): C keeps the
 * ARKLsMem behind the stepper's `void* lmem` and reaches it through
 * the step_getlinmem op.  Here the box lives directly on
 * `ARKodeMem.lmem`; the step_getlinmem op remains (take semantics)
 * and put-back writes the field.  This avoids a double-nested
 * take/put-back (step_mem -> lmem) in every ARKLS call.
 *
 * Mass-matrix interface (ARKLsMassMemRec and the arkLsMass* family)
 * is not yet ported; no serial example verified so far requires a
 * non-identity mass matrix.  It will land with arkode_arkstep mass
 * support.
 * -----------------------------------------------------------------*/
use crate::nvector_serial::NVector;
use crate::sundials_linearsolver::LinearSolver;
use crate::sundials_matrix::SUNMatrix;
use crate::sundials_types::UserData;

/* ARKLS return codes (arkode_ls.h) */
pub const ARKLS_SUCCESS: i32 = 0;
pub const ARKLS_MEM_NULL: i32 = -1;
pub const ARKLS_LMEM_NULL: i32 = -2;
pub const ARKLS_ILL_INPUT: i32 = -3;
pub const ARKLS_MEM_FAIL: i32 = -4;
pub const ARKLS_PMEM_NULL: i32 = -5;
pub const ARKLS_MASSMEM_NULL: i32 = -6;
pub const ARKLS_JACFUNC_UNRECVR: i32 = -7;
pub const ARKLS_JACFUNC_RECVR: i32 = -8;
pub const ARKLS_MASSFUNC_UNRECVR: i32 = -9;
pub const ARKLS_MASSFUNC_RECVR: i32 = -10;
pub const ARKLS_SUNMAT_FAIL: i32 = -11;
pub const ARKLS_SUNLS_FAIL: i32 = -12;

/* ARKLS solver constants (arkode_ls_impl.h):

   ARKLS_MSBJ   default maximum number of steps between Jacobian /
                preconditioner evaluations

   ARKLS_EPLIN  default value for factor by which the tolerance
                on the nonlinear iteration is multiplied to get
                a tolerance on the linear iteration */
pub const ARKLS_MSBJ: i64 = 51;
pub const ARKLS_EPLIN: f64 = 0.05;

/* User-supplied function types (arkode_ls.h) */
pub type ARKLsJacFn = fn(
    t: f64,
    y: &NVector,
    fy: &NVector,
    jac: &mut SUNMatrix,
    user_data: &mut UserData,
    tmp1: &mut NVector,
    tmp2: &mut NVector,
    tmp3: &mut NVector,
) -> i32;

pub type ARKLsMassFn = fn(
    t: f64,
    m: &mut SUNMatrix,
    user_data: &mut UserData,
    tmp1: &mut NVector,
    tmp2: &mut NVector,
    tmp3: &mut NVector,
) -> i32;

pub type ARKLsPrecSetupFn = fn(
    t: f64,
    y: &NVector,
    fy: &NVector,
    jok: bool,
    jcur_ptr: &mut bool,
    gamma: f64,
    user_data: &mut UserData,
) -> i32;

pub type ARKLsPrecSolveFn = fn(
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

pub type ARKLsMassPrecSetupFn = fn(t: f64, user_data: &mut UserData) -> i32;

pub type ARKLsMassPrecSolveFn =
    fn(t: f64, r: &NVector, z: &mut NVector, delta: f64, lr: i32, user_data: &mut UserData) -> i32;

pub type ARKLsJacTimesSetupFn =
    fn(t: f64, y: &NVector, fy: &NVector, user_data: &mut UserData) -> i32;

pub type ARKLsJacTimesVecFn = fn(
    v: &NVector,
    jv: &mut NVector,
    t: f64,
    y: &NVector,
    fy: &NVector,
    user_data: &mut UserData,
    tmp: &mut NVector,
) -> i32;

pub type ARKLsMassTimesSetupFn = fn(t: f64, mtimes_data: &mut UserData) -> i32;

pub type ARKLsMassTimesVecFn =
    fn(v: &NVector, mv: &mut NVector, t: f64, mtimes_data: &mut UserData) -> i32;

pub type ARKLsLinSysFn = fn(
    t: f64,
    y: &NVector,
    fy: &NVector,
    a: &mut SUNMatrix,
    m: Option<&SUNMatrix>,
    jok: bool,
    jcur: &mut bool,
    gamma: f64,
    user_data: &mut UserData,
    tmp1: &mut NVector,
    tmp2: &mut NVector,
    tmp3: &mut NVector,
) -> i32;

/* -----------------------------------------------------------------
   Types : ARKLsMemRec, ARKLsMem
   -----------------------------------------------------------------*/
pub struct ARKLsMem {
    /* Linear solver type information */
    pub iterative: bool,   /* is the solver iterative?    */
    pub matrixbased: bool, /* is a matrix structure used? */

    /* Jacobian construction & storage.
       C's J_data pointer (user_data for a user jac, ark_mem for the
       internal DQ approximation) collapses to the flag convention:
       jac == None + jacDQ  =>  internal arkLsDQJac(ark_mem). */
    pub jacDQ: bool,             /* SUNTRUE if using internal DQ Jacobian approx. */
    pub jac: Option<ARKLsJacFn>, /* Jacobian routine to be called                 */
    pub jbad: bool,              /* heuristic suggestion for pset                 */

    /* Matrix-based solver, scale solution to account for change in gamma */
    pub scalesol: bool,

    /* Iterative solver tolerance */
    pub eplifac: f64, /* nonlinear -> linear tol scaling factor  */
    pub nrmfac: f64,  /* integrator -> LS norm conversion factor */

    /* Linear solver, matrix and vector objects/pointers */
    pub LS: LinearSolver,          /* generic linear solver object       */
    pub A: Option<SUNMatrix>,      /* A = M - gamma * df/dy              */
    pub savedJ: Option<SUNMatrix>, /* savedJ = old Jacobian              */
    pub ytemp: NVector,            /* temp vector for jtimes & psolve    */
    pub x: NVector,                /* solution vector for SUNLinSolSolve */
    /* (C also stashes ycur/fcur pointers here for the duration of an
       arkLsSolve; the Rust solve passes them down as arguments.)     */

    /* Statistics and associated parameters */
    pub msbj: i64,     /* max num steps between jac/pset calls */
    pub tcur: f64,     /* 'time' for current ARKLs solve       */
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

    /* Preconditioner computation.
       C's P_data == user_data (user pset/psolve) or an internal
       preconditioner module (bandpre/bbdpre, with pfree); the module
       variants land with arkode_bandpre/arkode_bbdpre. */
    pub pset: Option<ARKLsPrecSetupFn>,
    pub psolve: Option<ARKLsPrecSolveFn>,

    /* Jacobian times vector computation.
       jtimes == None + jtimesDQ  =>  internal arkLsDQJtimes. */
    pub jtimesDQ: bool,
    pub jtsetup: Option<ARKLsJacTimesSetupFn>,
    pub jtimes: Option<ARKLsJacTimesVecFn>,
    pub Jt_f: Option<crate::arkode_impl::ARKRhsFn>,

    /* Linear system setup function.
       linsys == None + !user_linsys  =>  internal arkLsLinSys. */
    pub user_linsys: bool,
    pub linsys: Option<ARKLsLinSysFn>,

    pub last_flag: i32, /* last error flag returned by any function */
}

/* Error messages (arkode_ls_impl.h); only those used by the ported
   (non-mass) half are carried over. */
pub const MSG_LS_ARKMEM_NULL: &str = "Integrator memory is NULL.";
pub const MSG_LS_MEM_FAIL: &str = "A memory request failed.";
pub const MSG_LS_BAD_NVECTOR: &str = "A required vector operation is not implemented.";
pub const MSG_LS_LMEM_NULL: &str = "Linear solver memory is NULL.";
pub const MSG_LS_BAD_SIZES: &str =
    "Illegal bandwidth parameter(s). Must have 0 <=  ml, mu <= N-1.";
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
