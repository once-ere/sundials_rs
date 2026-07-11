/* -----------------------------------------------------------------
 * Translated from src/idas/idas_ls_impl.h and the constants / type
 * definitions of include/idas/idas_ls.h (IDAS 7.7.0).
 * IDALS linear-solver-interface memory structure and constants,
 * plus the PART II backward-problem (adjoint) memory structure.
 * Structural donor: ida_rs/src/ida_ls_impl.rs (verified Phase 4);
 * the forward-problem half is donor-verbatim.
 * -----------------------------------------------------------------*/
use crate::idas_impl::IDAResFn;
use crate::nvector_serial::NVector;
use crate::sundials_linearsolver::LinearSolver;
use crate::sundials_matrix::SUNMatrix;
use crate::sundials_types::UserData;

/* IDALS return codes (idas_ls.h) */
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

/* Return values for the adjoint module */
pub const IDALS_NO_ADJ: i32 = -101;
pub const IDALS_LMEMB_NULL: i32 = -102;

/* ===============================================================
   IDALS user-supplied function types (idas_ls.h).
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

/* C's IDALsPrecSetupFn is (tt, yy, yp, rr, c_j, user_data); a user pset
   that needs the current error-weight vector or step size fetches them
   inside the callback via IDAGetErrWeights / IDAGetCurrentStep on a stored
   ida_mem handle (e.g. idaFoodWeb_kry's Precond). A pure-Rust callback
   cannot re-borrow the integrator it is running inside, so — exactly as
   the port already hands `c_j` by value — the two integrator-internal
   quantities those getters return (ida_ewt, ida_hh) are passed in
   directly. They are byte-identical to what the C getters would copy out
   at pset time. Callbacks that don't need them ignore the extra args. */
pub type IDALsPrecSetupFn = fn(
    tt: f64,
    yy: &NVector,
    yp: &NVector,
    rr: &NVector,
    c_j: f64,
    ewt: &NVector,
    hh: f64,
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
   preconditioner module. */
#[derive(Default)]
pub enum PrecModule {
    #[default]
    None,
    /// user-supplied pset/psolve get user_data
    User,
    /* The BBDPre(Box<crate::idas_bbdpre_impl::IBBDPrecData>) variant
       lands with the idas_bbdpre_impl.h / idas_bbdpre.c units (see
       PROGRESS.md), mirroring the donor's ida_bbdpre variant and the
       LsModule::Ls placeholder precedent in idas_impl.rs. */
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
                              idas.c/idas_ic.c guard every lsetup
                              dispatch with `ida_lsetup != NULL`.
                              With enum dispatch the hook cannot be
                              NULLed, so this flag carries that state
                              (donor: IDALsMem.setup_disabled).       */
}

/* ===============================================================
   IDALS internal functions (prototypes in idas_ls_impl.h; the
   implementations land in idas_ls.rs):

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

   PART II internal functions (implementations land with the idaa /
   idas_ls PART II units):

     idaLsFreeB(IDAB_mem)                  — free the IDALsMemB block
     idaLs_AccessLMemB / idaLs_AccessLMemBCur
                                           — guarded access to IDALsMemB
   =============================================================== */

/* ===============================================================
   Error and Warning Messages (idas_ls_impl.h); C printf formats kept
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

/*-----------------------------------------------------------------
  PART II - backward problems
  -----------------------------------------------------------------*/

/* ===============================================================
   IDALS backward user-supplied function types (idas_ls.h).
   C's N_Vector* yS_1d / ypS_1d one-dimensional sensitivity arrays
   map to &[NVector]; user_dataB comes from the backward problem's
   IDABMem (C void* user_dataB), per the idas_impl.rs convention.
   (The forward IDALsPrecSetupFn carries extra ewt/hh arguments for
   the pure-Rust re-borrow problem; the backward pset types keep the
   C shapes — extend the same way only if a backward-pset example
   turns out to need integrator internals.)
   =============================================================== */

pub type IDALsJacFnB = fn(
    tt: f64,
    c_jB: f64,
    yy: &NVector,
    yp: &NVector,
    yyB: &NVector,
    ypB: &NVector,
    rrB: &NVector,
    jacB: &mut SUNMatrix,
    user_dataB: &mut UserData,
    tmp1B: &mut NVector,
    tmp2B: &mut NVector,
    tmp3B: &mut NVector,
) -> i32;

pub type IDALsJacFnBS = fn(
    tt: f64,
    c_jB: f64,
    yy: &NVector,
    yp: &NVector,
    yyS: &[NVector],
    ypS: &[NVector],
    yyB: &NVector,
    ypB: &NVector,
    rrB: &NVector,
    jacB: &mut SUNMatrix,
    user_dataB: &mut UserData,
    tmp1B: &mut NVector,
    tmp2B: &mut NVector,
    tmp3B: &mut NVector,
) -> i32;

pub type IDALsPrecSetupFnB = fn(
    tt: f64,
    yy: &NVector,
    yp: &NVector,
    yyB: &NVector,
    ypB: &NVector,
    rrB: &NVector,
    c_jB: f64,
    user_dataB: &mut UserData,
) -> i32;

pub type IDALsPrecSetupFnBS = fn(
    tt: f64,
    yy: &NVector,
    yp: &NVector,
    yyS: &[NVector],
    ypS: &[NVector],
    yyB: &NVector,
    ypB: &NVector,
    rrB: &NVector,
    c_jB: f64,
    user_dataB: &mut UserData,
) -> i32;

pub type IDALsPrecSolveFnB = fn(
    tt: f64,
    yy: &NVector,
    yp: &NVector,
    yyB: &NVector,
    ypB: &NVector,
    rrB: &NVector,
    rvecB: &NVector,
    zvecB: &mut NVector,
    c_jB: f64,
    deltaB: f64,
    user_dataB: &mut UserData,
) -> i32;

pub type IDALsPrecSolveFnBS = fn(
    tt: f64,
    yy: &NVector,
    yp: &NVector,
    yyS: &[NVector],
    ypS: &[NVector],
    yyB: &NVector,
    ypB: &NVector,
    rrB: &NVector,
    rvecB: &NVector,
    zvecB: &mut NVector,
    c_jB: f64,
    deltaB: f64,
    user_dataB: &mut UserData,
) -> i32;

pub type IDALsJacTimesSetupFnB = fn(
    t: f64,
    yy: &NVector,
    yp: &NVector,
    yyB: &NVector,
    ypB: &NVector,
    rrB: &NVector,
    c_jB: f64,
    user_dataB: &mut UserData,
) -> i32;

pub type IDALsJacTimesSetupFnBS = fn(
    t: f64,
    yy: &NVector,
    yp: &NVector,
    yyS: &[NVector],
    ypS: &[NVector],
    yyB: &NVector,
    ypB: &NVector,
    rrB: &NVector,
    c_jB: f64,
    user_dataB: &mut UserData,
) -> i32;

pub type IDALsJacTimesVecFnB = fn(
    t: f64,
    yy: &NVector,
    yp: &NVector,
    yyB: &NVector,
    ypB: &NVector,
    rrB: &NVector,
    vB: &NVector,
    JvB: &mut NVector,
    c_jB: f64,
    user_dataB: &mut UserData,
    tmp1B: &mut NVector,
    tmp2B: &mut NVector,
) -> i32;

pub type IDALsJacTimesVecFnBS = fn(
    t: f64,
    yy: &NVector,
    yp: &NVector,
    yyS: &[NVector],
    ypS: &[NVector],
    yyB: &NVector,
    ypB: &NVector,
    rrB: &NVector,
    vB: &NVector,
    JvB: &mut NVector,
    c_jB: f64,
    user_dataB: &mut UserData,
    tmp1B: &mut NVector,
    tmp2B: &mut NVector,
) -> i32;

/*-----------------------------------------------------------------
  Types : IDALsMemRecB, IDALsMemB

  IDASetLinearSolverB attaches such a structure to the lmemB
  field of IDAadjMem (via IDABMem).
  -----------------------------------------------------------------*/
#[derive(Default)]
pub struct IDALsMemB {
    pub jacB: Option<IDALsJacFnB>,
    pub jacBS: Option<IDALsJacFnBS>,
    pub jtsetupB: Option<IDALsJacTimesSetupFnB>,
    pub jtsetupBS: Option<IDALsJacTimesSetupFnBS>,
    pub jtimesB: Option<IDALsJacTimesVecFnB>,
    pub jtimesBS: Option<IDALsJacTimesVecFnBS>,
    pub psetB: Option<IDALsPrecSetupFnB>,
    pub psetBS: Option<IDALsPrecSetupFnBS>,
    pub psolveB: Option<IDALsPrecSolveFnB>,
    pub psolveBS: Option<IDALsPrecSolveFnBS>,
    /* (C P_dataB is passed to psetB/psolveB; here user_dataB comes
       from the backward problem's IDABMem, per the idas_impl.rs
       drop of the user-data self-pointers.) */
}

/*-----------------------------------------------------------------
  PART II Error Messages
  -----------------------------------------------------------------*/
pub const MSG_LS_CAMEM_NULL: &str = "idaadj_mem = NULL illegal.";
pub const MSG_LS_LMEMB_NULL: &str =
    "Linear solver memory is NULL for the backward integration.";
pub const MSG_LS_BAD_T: &str = "Bad t for interpolation.";
pub const MSG_LS_BAD_WHICH: &str = "Illegal value for which.";
pub const MSG_LS_NO_ADJ: &str = "Illegal attempt to call before calling IDAAdjInit.";

/* END of idas_ls_impl.h port. */
