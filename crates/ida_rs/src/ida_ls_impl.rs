/* -----------------------------------------------------------------
 * Translated from src/ida/ida_ls_impl.h and the constants / type
 * definitions of include/ida/ida_ls.h (IDA 7.7.0).
 * IDALS linear-solver-interface memory structure and constants.
 * Conventions follow the donor cvode_ls_impl.rs and the landed
 * kinsol_ls_impl.rs.
 * -----------------------------------------------------------------*/
use crate::ida_impl::IDAResFn;
use crate::nvector_serial::NVector;
use crate::sundials_linearsolver::LinearSolver;
use crate::sundials_matrix::SUNMatrix;
use crate::sundials_types::UserData;

/* IDALS return codes (ida_ls.h) */
pub const IDALS_SUCCESS: i32 = 0;
pub const IDALS_MEM_NULL: i32 = -1;
pub const IDALS_LMEM_NULL: i32 = -2;
pub const IDALS_ILL_INPUT: i32 = -3;
pub const IDALS_MEM_FAIL: i32 = -4;
pub const IDALS_PMEM_NULL: i32 = -5;
pub const IDALS_JACFUNC_UNRECVR: i32 = -6;
pub const IDALS_JACFUNC_RECVR: i32 = -7;
pub const IDALS_SUNMAT_FAIL: i32 = -8;
pub const IDALS_SUNLS_FAIL: i32 = -9;

/* ===============================================================
   IDALS user-supplied function types (ida_ls.h).
   All carry the scalar cj (the -alphas/hh Jacobian coefficient,
   J = dF/dy + cj * dF/dy') threaded from the integrator.
   =============================================================== */

pub type IDALsJacFn = fn(
    t: f64,
    c_j: f64,
    yy: &NVector,
    yp: &NVector,
    rr: &NVector,
    jac: &mut SUNMatrix,
    user_data: &mut UserData,
    tmp1: &mut NVector,
    tmp2: &mut NVector,
    tmp3: &mut NVector,
) -> i32;

pub type IDALsPrecSetupFn = fn(
    tt: f64,
    yy: &NVector,
    yp: &NVector,
    rr: &NVector,
    c_j: f64,
    user_data: &mut UserData,
) -> i32;

pub type IDALsPrecSolveFn = fn(
    tt: f64,
    yy: &NVector,
    yp: &NVector,
    rr: &NVector,
    rvec: &NVector,
    zvec: &mut NVector,
    c_j: f64,
    delta: f64,
    user_data: &mut UserData,
) -> i32;

pub type IDALsJacTimesSetupFn = fn(
    tt: f64,
    yy: &NVector,
    yp: &NVector,
    rr: &NVector,
    c_j: f64,
    user_data: &mut UserData,
) -> i32;

pub type IDALsJacTimesVecFn = fn(
    tt: f64,
    yy: &NVector,
    yp: &NVector,
    rr: &NVector,
    v: &NVector,
    Jv: &mut NVector,
    c_j: f64,
    user_data: &mut UserData,
    tmp1: &mut NVector,
    tmp2: &mut NVector,
) -> i32;

/* Preconditioner module attached to IDALS.
   In C this is the pdata/pfree convention: pdata points at
   user_data (user-supplied pset/psolve) or at an internal
   preconditioner module (ida_bbdpre; its variant is added here when
   that module lands). */
#[derive(Default)]
pub enum PrecModule {
    #[default]
    None,
    /// user-supplied pset/psolve get user_data
    User,
}

/* -----------------------------------------------------------------
   Types : IDALsMemRec, IDALsMem
   -----------------------------------------------------------------*/
pub struct IDALsMem {
    /* Linear solver type information */
    pub iterative: bool,   /* is the solver iterative?    */
    pub matrixbased: bool, /* is a matrix structure used? */

    /* Jacobian construction & storage */
    pub jacDQ: bool,             /* SUNTRUE if using internal DQ Jacobian approx. */
    pub jac: Option<IDALsJacFn>, /* Jacobian routine to be called                 */
    /* (C J_data is passed to jac; here user_data comes from IDAMem)   */

    /* Linear solver, matrix and vector objects/pointers */
    pub LS: LinearSolver,     /* generic linear solver object                  */
    pub J: Option<SUNMatrix>, /* J = dF/dy + cj*dF/dy'                         */
    pub ytemp: NVector,       /* temp vector used by IDAAtimesDQ               */
    pub yptemp: NVector,      /* temp vector used by IDAAtimesDQ               */
    pub x: NVector,           /* temp vector used by the solve function        */
    /* (C also stores ycur/ypcur/rcur — borrowed pointers to the
       current Newton-iteration vectors, installed by idaLsSetup /
       idaLsSolve for use by the ATimes/PSetup/PSolve callbacks; here
       they are passed as arguments instead, as in the donor's drop
       of the CVLS ycur/fcur aliases.)                                */

    /* Matrix-based solver, scale solution to account for change in cj */
    pub scalesol: bool,

    /* Iterative solver tolerance */
    pub eplifac: f64, /* nonlinear -> linear tol scaling factor       */
    pub nrmfac: f64,  /* integrator -> LS norm conversion factor      */

    /* Statistics and associated parameters */
    pub dqincfac: f64, /* dqincfac = optional increment factor in Jv   */
    pub nje: i64,      /* nje = no. of calls to jac                    */
    pub npe: i64,      /* npe = total number of precond calls          */
    pub nli: i64,      /* nli = total number of linear iterations      */
    pub nps: i64,      /* nps = total number of psolve calls           */
    pub ncfl: i64,     /* ncfl = total number of convergence failures  */
    pub nreDQ: i64,    /* nreDQ = total number of calls to res         */
    pub njtsetup: i64, /* njtsetup = total number of calls to jtsetup  */
    pub njtimes: i64,  /* njtimes = total number of calls to jtimes    */
    pub nst0: i64,     /* nst0 = saved nst (for performance monitor)   */
    pub nni0: i64,     /* nni0 = saved nni (for performance monitor)   */
    pub ncfn0: i64,    /* ncfn0 = saved ncfn (for performance monitor) */
    pub ncfl0: i64,    /* ncfl0 = saved ncfl (for performance monitor) */
    pub nwarn: i64,    /* nwarn = no. of warnings (for perf. monitor)  */
    pub nstlj: i64,    /* nstlj = nst at last jac/pset call            */
    pub tnlj: f64,     /* tnlj = t_n at last jac/pset call             */

    pub last_flag: i32, /* last error return flag                       */

    /* Preconditioner computation
       (C pdata/pfree convention -> PrecModule enum, RAII) */
    pub pset: Option<IDALsPrecSetupFn>,
    pub psolve: Option<IDALsPrecSolveFn>,
    pub prec_module: PrecModule,

    /* Jacobian times vector computation
       (a) jtimes function provided by the user: jtimesDQ == SUNFALSE
       (b) internal jtimes: jtimesDQ == SUNTRUE
       (C jt_data is passed to jtimes; here user_data comes from IDAMem) */
    pub jtimesDQ: bool,
    pub jtsetup: Option<IDALsJacTimesSetupFn>,
    pub jtimes: Option<IDALsJacTimesVecFn>,
    pub jt_res: Option<IDAResFn>,

    pub setup_disabled: bool, /* In C, idaLsInitialize NULLs the
                              ida_lsetup hook when J == NULL and no
                              preconditioner setup (pset) is supplied,
                              or when the LS is matrix-embedded, and
                              ida.c/ida_ic.c guard every lsetup
                              dispatch with `ida_lsetup != NULL`.
                              With enum dispatch the hook cannot be
                              NULLed, so this flag carries that state
                              (donor: CVLsMem.setup_disabled).        */
}

/* ===============================================================
   IDALS internal functions (prototypes in ida_ls_impl.h; the
   implementations land in ida_ls.rs):

     idaLsATimes(ida_mem, v, z)            — interface to SUNLinearSolver ATimes
     idaLsPSetup(ida_mem)                  — interface to SUNLinearSolver PSetup
     idaLsPSolve(ida_mem, r, z, tol, lr)   — interface to SUNLinearSolver PSolve
     idaLsDQJtimes(tt, yy, yp, rr, v, Jv, c_j, data, work1, work2)
                                           — DQ approximation to Jv
     idaLsDQJac / idaLsDenseDQJac / idaLsBandDQJac
                                           — difference-quotient Jacobians
     idaLsInitialize / idaLsSetup / idaLsSolve / idaLsPerf / idaLsFree
                                           — generic linit/lsetup/lsolve/
                                             lperf/lfree (dispatched via
                                             LsModule::Ls)
     idaLsInitializeCounters(idals_mem)    — reset the counters above
     idaLs_AccessLMem(ida_mem, fname, ...) — guarded access to IDALsMem
   =============================================================== */

/* ===============================================================
   Error and Warning Messages (ida_ls_impl.h); C printf formats kept
   verbatim with SUN_FORMAT_G expanded to its double-precision form
   "%.15g" (call sites render the values via sundials_utils::fmt_g).
   =============================================================== */

pub const MSG_LS_TIME: &str = "at t = %.15g, ";
pub const MSG_LS_FRMT: &str = "%.15g.";

/* Error Messages */
pub const MSG_LS_IDAMEM_NULL: &str = "Integrator memory is NULL.";
pub const MSG_LS_MEM_FAIL: &str = "A memory request failed.";
pub const MSG_LS_BAD_NVECTOR: &str = "A required vector operation is not implemented.";
pub const MSG_LS_BAD_SIZES: &str =
    "Illegal bandwidth parameter(s). Must have 0 <=  ml, mu <= N-1.";
pub const MSG_LS_BAD_LSTYPE: &str = "Incompatible linear solver type.";
pub const MSG_LS_LMEM_NULL: &str = "Linear solver memory is NULL.";
pub const MSG_LS_BAD_GSTYPE: &str = "gstype has an illegal value.";
pub const MSG_LS_NEG_MAXRS: &str = "maxrs < 0 illegal.";
pub const MSG_LS_NEG_EPLIFAC: &str = "eplifac < 0.0 illegal.";
pub const MSG_LS_NEG_DQINCFAC: &str = "dqincfac < 0.0 illegal.";
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
pub const MSG_LS_MATZERO_FAILED: &str =
    "The SUNMatZero routine failed in an unrecoverable manner.";

/* Warning Messages */
pub const MSG_LS_WARN: &str =
    "Warning: at t = %.15g, poor iterative algorithm performance. ";
pub const MSG_LS_CFN_WARN: &str = "Warning: at t = %.15g, \
poor iterative algorithm performance. Nonlinear convergence failure rate is %.15g.";
pub const MSG_LS_CFL_WARN: &str = "Warning: at t = %.15g, \
poor iterative algorithm performance. Linear convergence failure rate is %.15g.";
