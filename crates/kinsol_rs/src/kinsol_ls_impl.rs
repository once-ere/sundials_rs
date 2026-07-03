/* -----------------------------------------------------------------
 * Translated from src/kinsol/kinsol_ls_impl.h and the type
 * definitions of include/kinsol/kinsol_ls.h (KINSOL 7.7.0).
 * KINLS linear-solver-interface memory structure and constants.
 * Conventions follow the donor cvode_ls_impl.rs.
 * -----------------------------------------------------------------*/
use crate::kinsol_impl::KINSysFn;
use crate::nvector_serial::NVector;
use crate::sundials_linearsolver::LinearSolver;
use crate::sundials_matrix::SUNMatrix;
use crate::sundials_types::UserData;

/* KINLS return codes (kinsol_ls.h) */
pub const KINLS_SUCCESS: i32 = 0;

pub const KINLS_MEM_NULL: i32 = -1;
pub const KINLS_LMEM_NULL: i32 = -2;
pub const KINLS_ILL_INPUT: i32 = -3;
pub const KINLS_MEM_FAIL: i32 = -4;
pub const KINLS_PMEM_NULL: i32 = -5;
pub const KINLS_JACFUNC_ERR: i32 = -6;
pub const KINLS_SUNMAT_FAIL: i32 = -7;
pub const KINLS_SUNLS_FAIL: i32 = -8;

/* keys for KINPrintInfo (do not use 1 -> conflict with PRNT_RETVAL) */
pub const PRNT_NLI: i32 = 101;
pub const PRNT_EPS: i32 = 102;

/* ===============================================================
   KINLS user-supplied function types (kinsol_ls.h)
   =============================================================== */

pub type KINLsJacFn = fn(
    u: &NVector,
    fu: &NVector,
    jac: &mut SUNMatrix,
    user_data: &mut UserData,
    tmp1: &mut NVector,
    tmp2: &mut NVector,
) -> i32;

pub type KINLsPrecSetupFn = fn(
    uu: &NVector,
    uscale: &NVector,
    fval: &NVector,
    fscale: &NVector,
    user_data: &mut UserData,
) -> i32;

pub type KINLsPrecSolveFn = fn(
    uu: &NVector,
    uscale: &NVector,
    fval: &NVector,
    fscale: &NVector,
    vv: &mut NVector,
    user_data: &mut UserData,
) -> i32;

pub type KINLsJacTimesVecFn = fn(
    v: &NVector,
    jv: &mut NVector,
    uu: &NVector,
    new_uu: &mut bool,
    user_data: &mut UserData,
) -> i32;

/* Preconditioner module attached to KINLS.
   In C this is the pdata/pfree convention: pdata points at
   user_data (user-supplied pset/psolve) or at an internal
   preconditioner module (kinsol_bbdpre; its variant is added here
   when that module lands). */
#[derive(Default)]
pub enum PrecModule {
    #[default]
    None,
    /// user-supplied pset/psolve get user_data
    User,
}

/* -----------------------------------------------------------------
   Types : KINLsMemRec, KINLsMem
   -----------------------------------------------------------------*/
pub struct KINLsMem {
    /* Linear solver type information */
    pub iterative: bool,   /* is the solver iterative?    */
    pub matrixbased: bool, /* is a matrix structure used? */

    /* Jacobian construction & storage */
    pub jacDQ: bool,             /* SUNTRUE if using internal DQ Jacobian approx. */
    pub jac: Option<KINLsJacFn>, /* Jacobian routine to be called                 */
    /* (C J_data is passed to jac; here user_data comes from KINMem)  */

    /* Linear solver, matrix and vector objects/pointers */
    pub LS: LinearSolver,     /* generic iterative linear solver object */
    pub J: Option<SUNMatrix>, /* problem Jacobian                       */

    /* Solver tolerance adjustment factor (if needed, see kinLsSolve) */
    pub tol_fac: f64,

    /* Statistics and associated parameters */
    pub nje: i64,     /* no. of calls to jac                           */
    pub nfeDQ: i64,   /* no. of calls to F due to DQ Jacobian or J*v
                         approximations                                */
    pub npe: i64,     /* npe = total number of precond calls           */
    pub nli: i64,     /* nli = total number of linear iterations       */
    pub nps: i64,     /* nps = total number of psolve calls            */
    pub ncfl: i64,    /* ncfl = total number of convergence failures   */
    pub njtimes: i64, /* njtimes = total number of calls to jtimes     */

    pub new_uu: bool, /* flag indicating if the iterate has been
                         updated - the Jacobian must be updated or
                         reevaluated (meant to be used by a
                         user-supplied jtimes function                 */

    pub last_flag: i32, /* last error return flag                      */

    /* Preconditioner computation
       (C pdata/pfree convention -> PrecModule enum, RAII) */
    pub pset: Option<KINLsPrecSetupFn>,
    pub psolve: Option<KINLsPrecSolveFn>,
    pub prec_module: PrecModule,

    /* Jacobian times vector computation
       (a) jtimes function provided by the user: jtimesDQ == SUNFALSE
       (b) internal jtimes: jtimesDQ == SUNTRUE
       (C jt_data is passed to jtimes; here user_data comes from KINMem) */
    pub jtimesDQ: bool,
    pub jtimes: Option<KINLsJacTimesVecFn>,
    pub jt_func: Option<KINSysFn>,
}

/* ===============================================================
   KINLS internal functions (prototypes in kinsol_ls_impl.h; the
   implementations land in kinsol_ls.rs):

     kinLsATimes(kinmem, v, z)               — interface to SUNLinearSolver ATimes
     kinLsPSetup(kinmem)                     — interface to SUNLinearSolver PSetup
     kinLsPSolve(kinmem, r, z, tol, lr)      — interface to SUNLinearSolver PSolve
     kinLsDQJtimes(v, Jv, u, new_u, data)    — DQ approximation to J*v
     kinLsDQJac / kinLsDenseDQJac / kinLsBandDQJac
                                             — difference-quotient Jacobians
     kinLsInitialize / kinLsSetup / kinLsSolve / kinLsFree
                                             — generic linit/lsetup/lsolve/lfree
                                               (dispatched via LsModule::Ls)
     kinLsInitializeCounters(kinls_mem)      — reset the counters above
     kinLs_AccessLMem(kinmem, fname, ...)    — guarded access to KINLsMem
   =============================================================== */

/* Error messages (kinsol_ls_impl.h) */
pub const MSG_LS_KINMEM_NULL: &str = "KINSOL memory is NULL.";
pub const MSG_LS_MEM_FAIL: &str = "A memory request failed.";
pub const MSG_LS_BAD_NVECTOR: &str = "A required vector operation is not implemented.";
pub const MSG_LS_LMEM_NULL: &str = "Linear solver memory is NULL.";
pub const MSG_LS_NEG_MAXRS: &str = "maxrs < 0 illegal.";
pub const MSG_LS_BAD_SIZES: &str =
    "Illegal bandwidth parameter(s). Must have 0 <=  ml, mu <= N-1.";

pub const MSG_LS_JACFUNC_FAILED: &str =
    "The Jacobian routine failed in an unrecoverable manner.";
pub const MSG_LS_PSET_FAILED: &str =
    "The preconditioner setup routine failed in an unrecoverable manner.";
pub const MSG_LS_PSOLVE_FAILED: &str =
    "The preconditioner solve routine failed in an unrecoverable manner.";
pub const MSG_LS_JTIMES_FAILED: &str =
    "The Jacobian x vector routine failed in an unrecoverable manner.";
pub const MSG_LS_MATZERO_FAILED: &str =
    "The SUNMatZero routine failed in an unrecoverable manner.";

/* Info messages (kinsol_ls_impl.h); C printf formats kept verbatim
   with SUN_FORMAT_G expanded to its double-precision form "%.15g"
   (KINPrintInfo renders them via sundials_utils::fmt_g). */
pub const INFO_NLI: &str = "nli_inc = %d";
pub const INFO_EPS: &str = "residual norm = %.15g  eps = %.15g";
