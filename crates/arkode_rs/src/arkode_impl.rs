/* -----------------------------------------------------------------
 * Translation of sundials-7.7.0/src/arkode/arkode_impl.h
 * (+ arkode_types_impl.h and the constants / user-supplied function
 * types of include/arkode/arkode.h).
 *
 * Modeling (mirrors the C architecture; pinned for the whole crate):
 *  - ARKodeMem's time-stepper interface is a table of Option<fn>
 *    fields taking &mut ARKodeMem, installed by each stepper's
 *    create routine exactly as in C; `step_mem` is
 *    Option<Box<dyn Any>> downcast by the owning stepper module.
 *  - The generic ARKInterp object (void* content + ops table in C)
 *    is enum dispatch over the Hermite/Lagrange contents
 *    (arkode_interp_impl.rs); dispatchers live in arkode_interp.rs.
 *  - `fn` is reserved in Rust: the ARKodeMem field is `fn_`
 *    (single documented rename; everything else keeps C names).
 *  - C's `ycur` aliases the user-provided output vector during
 *    ARKodeEvolve; the Rust field owns storage and ARKodeEvolve
 *    copies back at every return path (CLAUDE.md rule 5, same as
 *    the cvode port's cv_y).
 *  - `rwt_is_ewt` in C makes rwt literally alias ewt; here rwt is
 *    a separate vector that arkode.c keeps in sync / reads through
 *    a helper when the flag is set.
 *  - FILE* parameters are &mut dyn std::io::Write; the `python`
 *    field and XBraid hooks are excluded backends.
 *  - sunctx/profiler/logger plumbing is compiled out as in the
 *    other crates.
 * -----------------------------------------------------------------*/

use crate::arkode_adapt_impl::ARKodeHAdaptMem;
use crate::arkode_interp_impl::{ARKInterpContent_Hermite, ARKInterpContent_Lagrange};
use crate::arkode_relaxation_impl::ARKodeRelaxMem;
use crate::arkode_root_impl::ARKodeRootMem;
use crate::nvector_serial::NVector;
use crate::sundials_adaptcontroller::SUNAdaptController;
use crate::sundials_adjointcheckpointscheme::SUNAdjointCheckpointScheme;
use crate::sundials_linearsolver::SUNLinearSolver_Type;
use crate::sundials_nonlinearsolver::NonlinearSolver;
use crate::sundials_types::{suncountertype, SUNOutputFormat, UserData};

/*===============================================================
  arkode.h — ARKODE Constants
  ===============================================================*/

/* usage modes (itask) */
pub const ARK_NORMAL: i32 = 1;
pub const ARK_ONE_STEP: i32 = 2;

/* adaptivity module flags */
pub const ARK_ADAPT_CUSTOM: i32 = -1;
pub const ARK_ADAPT_PID: i32 = 0;
pub const ARK_ADAPT_PI: i32 = 1;
pub const ARK_ADAPT_I: i32 = 2;
pub const ARK_ADAPT_EXP_GUS: i32 = 3;
pub const ARK_ADAPT_IMP_GUS: i32 = 4;
pub const ARK_ADAPT_IMEX_GUS: i32 = 5;

/* Constants for evaluating the full RHS */
pub const ARK_FULLRHS_START: i32 = 0;
pub const ARK_FULLRHS_END: i32 = 1;
pub const ARK_FULLRHS_OTHER: i32 = 2;

/* interpolation module flags */
pub const ARK_INTERP_MAX_DEGREE: i32 = 5;
pub const ARK_INTERP_NONE: i32 = -1;
pub const ARK_INTERP_HERMITE: i32 = 0;
pub const ARK_INTERP_LAGRANGE: i32 = 1;

/* return values */
pub const ARK_SUCCESS: i32 = 0;
pub const ARK_TSTOP_RETURN: i32 = 1;
pub const ARK_ROOT_RETURN: i32 = 2;

pub const ARK_WARNING: i32 = 99;

pub const ARK_TOO_MUCH_WORK: i32 = -1;
pub const ARK_TOO_MUCH_ACC: i32 = -2;
pub const ARK_ERR_FAILURE: i32 = -3;
pub const ARK_CONV_FAILURE: i32 = -4;

pub const ARK_LINIT_FAIL: i32 = -5;
pub const ARK_LSETUP_FAIL: i32 = -6;
pub const ARK_LSOLVE_FAIL: i32 = -7;
pub const ARK_RHSFUNC_FAIL: i32 = -8;
pub const ARK_FIRST_RHSFUNC_ERR: i32 = -9;
pub const ARK_REPTD_RHSFUNC_ERR: i32 = -10;
pub const ARK_UNREC_RHSFUNC_ERR: i32 = -11;
pub const ARK_RTFUNC_FAIL: i32 = -12;
pub const ARK_LFREE_FAIL: i32 = -13;
pub const ARK_MASSINIT_FAIL: i32 = -14;
pub const ARK_MASSSETUP_FAIL: i32 = -15;
pub const ARK_MASSSOLVE_FAIL: i32 = -16;
pub const ARK_MASSFREE_FAIL: i32 = -17;
pub const ARK_MASSMULT_FAIL: i32 = -18;

pub const ARK_CONSTR_FAIL: i32 = -19;
pub const ARK_MEM_FAIL: i32 = -20;
pub const ARK_MEM_NULL: i32 = -21;
pub const ARK_ILL_INPUT: i32 = -22;
pub const ARK_NO_MALLOC: i32 = -23;
pub const ARK_BAD_K: i32 = -24;
pub const ARK_BAD_T: i32 = -25;
pub const ARK_BAD_DKY: i32 = -26;
pub const ARK_TOO_CLOSE: i32 = -27;

pub const ARK_VECTOROP_ERR: i32 = -28;

pub const ARK_NLS_INIT_FAIL: i32 = -29;
pub const ARK_NLS_SETUP_FAIL: i32 = -30;
pub const ARK_NLS_SETUP_RECVR: i32 = -31;
pub const ARK_NLS_OP_ERR: i32 = -32;

pub const ARK_INNERSTEP_ATTACH_ERR: i32 = -33;
pub const ARK_INNERSTEP_FAIL: i32 = -34;
pub const ARK_OUTERTOINNER_FAIL: i32 = -35;
pub const ARK_INNERTOOUTER_FAIL: i32 = -36;

/* ARK_POSTPROCESS_FAIL equals ARK_POSTPROCESS_STEP_FAIL
   for backwards compatibility. */
pub const ARK_POSTPROCESS_FAIL: i32 = -37;
pub const ARK_POSTPROCESS_STEP_FAIL: i32 = -37;
pub const ARK_POSTPROCESS_STAGE_FAIL: i32 = -38;
pub const ARK_PRESTEPFN_FAIL: i32 = -39;
pub const ARK_POSTSTEPFN_FAIL: i32 = -40;
pub const ARK_PRERHSFN_FAIL: i32 = -41;

pub const ARK_USER_PREDICT_FAIL: i32 = -42;
pub const ARK_INTERP_FAIL: i32 = -43;

pub const ARK_INVALID_TABLE: i32 = -44;

pub const ARK_CONTEXT_ERR: i32 = -45;

pub const ARK_RELAX_FAIL: i32 = -46;
pub const ARK_RELAX_MEM_NULL: i32 = -47;
pub const ARK_RELAX_FUNC_FAIL: i32 = -48;
pub const ARK_RELAX_JAC_FAIL: i32 = -49;

pub const ARK_CONTROLLER_ERR: i32 = -50;

pub const ARK_STEPPER_UNSUPPORTED: i32 = -51;

pub const ARK_DOMEIG_FAIL: i32 = -52;
pub const ARK_MAX_STAGE_LIMIT_FAIL: i32 = -53;

pub const ARK_SUNSTEPPER_ERR: i32 = -54;
pub const ARK_STEP_DIRECTION_ERR: i32 = -55;

pub const ARK_ADJ_CHECKPOINT_FAIL: i32 = -56;
pub const ARK_ADJ_RECOMPUTE_FAIL: i32 = -57;
pub const ARK_SUNADJSTEPPER_ERR: i32 = -58;

pub const ARK_DEE_FAIL: i32 = -59;

pub const ARK_STEP_H0_FAIL: i32 = -60;

pub const ARK_UNRECOGNIZED_ERROR: i32 = -99;

/*===============================================================
  arkode.h — User-Supplied Function Types
  ===============================================================*/

pub type ARKRhsFn =
    fn(t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32;

pub type ARKRootFn =
    fn(t: f64, y: &NVector, gout: &mut [f64], user_data: &mut UserData) -> i32;

pub type ARKEwtFn = fn(y: &NVector, ewt: &mut NVector, user_data: &mut UserData) -> i32;

pub type ARKRwtFn = fn(y: &NVector, rwt: &mut NVector, user_data: &mut UserData) -> i32;

pub type ARKAdaptFn = fn(
    y: &NVector,
    t: f64,
    h1: f64,
    h2: f64,
    h3: f64,
    e1: f64,
    e2: f64,
    e3: f64,
    q: i32,
    p: i32,
    hnew: &mut f64,
    user_data: &mut UserData,
) -> i32;

pub type ARKExpStabFn =
    fn(y: &NVector, t: f64, hstab: &mut f64, user_data: &mut UserData) -> i32;

pub type ARKVecResizeFn =
    fn(y: &mut NVector, ytemplate: &NVector, user_data: &mut UserData) -> i32;

pub type ARKPreStepFn =
    fn(t: f64, y: &NVector, step: i64, attempt: i32, user_data: &mut UserData) -> i32;

pub type ARKPostStepFn = fn(t: f64, y: &NVector, step: i64, user_data: &mut UserData) -> i32;

pub type ARKPostProcessFn = fn(t: f64, y: &NVector, user_data: &mut UserData) -> i32;

pub type ARKPreRhsFn = fn(t: f64, y: &NVector, user_data: &mut UserData) -> i32;

pub type ARKStagePredictFn = fn(t: f64, zpred: &mut NVector, user_data: &mut UserData) -> i32;

pub type ARKRelaxFn = fn(y: &NVector, r: &mut f64, user_data: &mut UserData) -> i32;

pub type ARKRelaxJacFn = fn(y: &NVector, j: &mut NVector, user_data: &mut UserData) -> i32;

/// enum ARKRelaxSolver (arkode.h)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ARKRelaxSolver {
    ARK_RELAX_BRENT,
    ARK_RELAX_NEWTON,
}
pub use ARKRelaxSolver::{ARK_RELAX_BRENT, ARK_RELAX_NEWTON};

/// enum ARKAccumError (arkode.h)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ARKAccumError {
    ARK_ACCUMERROR_NONE,
    ARK_ACCUMERROR_MAX,
    ARK_ACCUMERROR_SUM,
    ARK_ACCUMERROR_AVG,
}
pub use ARKAccumError::{
    ARK_ACCUMERROR_AVG, ARK_ACCUMERROR_MAX, ARK_ACCUMERROR_NONE, ARK_ACCUMERROR_SUM,
};

/*===============================================================
  ARKODE Private Constants (arkode_impl.h)
  ===============================================================*/

/* Basic ARKODE defaults */
pub const Q_DEFAULT: i32 = 4; /* method order */
pub const MXSTEP_DEFAULT: i64 = 500; /* max steps between returns */
pub const MAXNEF: i32 = 7; /* max number of error failures */
pub const MAXNCF: i32 = 10; /* max number of convergence failures */
pub const MAXCONSTRFAILS: i32 = 10; /* max number of constraint failures */
pub const MXHNIL: i32 = 10; /* max number of t+h==h warnings */
pub const MAX_DQITERS: i32 = 3; /* max number of attempts to recover in DQ J*v */

/* Numeric constants */
pub const ZERO: f64 = 0.0;
pub const TINY: f64 = 1.0e-10;
pub const TENTH: f64 = 0.1;
pub const HALF: f64 = 0.5;
pub const ONE: f64 = 1.0;
pub const TWO: f64 = 2.0;
pub const THREE: f64 = 3.0;
pub const FOUR: f64 = 4.0;
pub const FIVE: f64 = 5.0;

/* Control constants for tolerances */
pub const ARK_SS: i32 = 0;
pub const ARK_SV: i32 = 1;
pub const ARK_WF: i32 = 2;

/* Initialization types */
pub const FIRST_INIT: i32 = 0; /* first step (re-)initialization */
pub const RESET_INIT: i32 = 1; /* reset initialization           */
pub const RESIZE_INIT: i32 = 2; /* resize initialization          */

/* Control constants for lower-level time-stepping functions */
pub const PREDICT_AGAIN: i32 = 3;
pub const CONV_FAIL: i32 = 4;
pub const TRY_AGAIN: i32 = 5;
pub const FIRST_CALL: i32 = 6;
pub const PREV_CONV_FAIL: i32 = 7;
pub const PREV_ERR_FAIL: i32 = 8;
pub const RHSFUNC_RECVR: i32 = 9;
pub const CONSTR_RECVR: i32 = 10;
pub const ARK_RETRY_STEP: i32 = 11;

/* Return values for lower-level rootfinding functions */
pub const RTFOUND: i32 = 1;
pub const CLOSERT: i32 = 3;

/* Algorithmic constants */
pub const FUZZ_FACTOR: f64 = 100.0;

pub const H0_LBFACTOR: f64 = 100.0;
pub const H0_UBFACTOR: f64 = 0.1;
pub const H0_BIAS: f64 = HALF;
pub const H0_ITERS: i32 = 4;

pub const ONEPSM: f64 = 1.000001;
pub const ONEMSM: f64 = 0.999999;

/* Input flag to linear solver setup routine: CONVFAIL */
pub const ARK_NO_FAILURES: i32 = 0;
pub const ARK_FAIL_BAD_J: i32 = 1;
pub const ARK_FAIL_OTHER: i32 = 2;

/*===============================================================
  ARKODE Interface function definitions
  ===============================================================*/

/* linear solver interface functions */
pub type ARKLinsolInitFn = fn(ark_mem: &mut ARKodeMem) -> i32;
pub type ARKLinsolSetupFn = fn(
    ark_mem: &mut ARKodeMem,
    convfail: i32,
    tpred: f64,
    /* &mut: the internal DQ Jacobian perturbs ypred in place and
       restores it, exactly as the C code perturbs the caller's vector */
    ypred: &mut NVector,
    fpred: &NVector,
    jcurPtr: &mut bool,
    vtemp1: &mut NVector,
    vtemp2: &mut NVector,
    vtemp3: &mut NVector,
) -> i32;
pub type ARKLinsolSolveFn = fn(
    ark_mem: &mut ARKodeMem,
    b: &mut NVector,
    tcur: f64,
    ycur: &NVector,
    fcur: &NVector,
    client_tol: f64,
    mnewt: i32,
) -> i32;
pub type ARKLinsolFreeFn = fn(ark_mem: &mut ARKodeMem) -> i32;

/* mass matrix solver interface functions */
pub type ARKMassInitFn = fn(ark_mem: &mut ARKodeMem) -> i32;
pub type ARKMassSetupFn = fn(
    ark_mem: &mut ARKodeMem,
    t: f64,
    vtemp1: &mut NVector,
    vtemp2: &mut NVector,
    vtemp3: &mut NVector,
) -> i32;
pub type ARKMassMultFn = fn(arkode_mem: &mut ARKodeMem, v: &NVector, mv: &mut NVector) -> i32;
pub type ARKMassSolveFn = fn(ark_mem: &mut ARKodeMem, b: &mut NVector, client_tol: f64) -> i32;
pub type ARKMassFreeFn = fn(ark_mem: &mut ARKodeMem) -> i32;

/* time stepper interface functions -- general */
pub type ARKTimestepInitFn = fn(ark_mem: &mut ARKodeMem, init_type: i32) -> i32;
pub type ARKTimestepFullRHSFn =
    fn(ark_mem: &mut ARKodeMem, t: f64, y: &NVector, f: &mut NVector, mode: i32) -> i32;
pub type ARKTimestepStepFn = fn(ark_mem: &mut ARKodeMem, dsm: &mut f64, nflag: &mut i32) -> i32;
/// C passes user_data; the Rust stepper reads ark_mem.user_data.
pub type ARKTimetepSetUserDataFn = fn(ark_mem: &mut ARKodeMem) -> i32;
pub type ARKTimestepPrintAllStats =
    fn(ark_mem: &mut ARKodeMem, outfile: &mut dyn std::io::Write, fmt: SUNOutputFormat) -> i32;
pub type ARKTimestepWriteParameters =
    fn(ark_mem: &mut ARKodeMem, fp: &mut dyn std::io::Write) -> i32;
pub type ARKTimestepResize = fn(
    ark_mem: &mut ARKodeMem,
    ynew: &NVector,
    hscale: f64,
    t0: f64,
    resize: Option<ARKVecResizeFn>,
    resize_data: &mut UserData,
) -> i32;
pub type ARKTimestepReset = fn(ark_mem: &mut ARKodeMem, tR: f64, yR: &NVector) -> i32;
pub type ARKTimestepFree = fn(ark_mem: &mut ARKodeMem);
pub type ARKTimestepPrintMem = fn(ark_mem: &mut ARKodeMem, outfile: &mut dyn std::io::Write);
pub type ARKTimestepSetDefaults = fn(ark_mem: &mut ARKodeMem) -> i32;
pub type ARKTimestepSetOrder = fn(ark_mem: &mut ARKodeMem, maxord: i32) -> i32;
pub type ARKTimestepGetNumRhsEvals =
    fn(ark_mem: &mut ARKodeMem, partition_index: i32, num_rhs_evals: &mut i64) -> i32;
pub type ARKTimestepSetStepDirection = fn(ark_mem: &mut ARKodeMem, stepdir: f64) -> i32;
pub type ARKTimestepSetUseCompensatedSums = fn(ark_mem: &mut ARKodeMem, onoff: bool) -> i32;
pub type ARKTimestepSetOptions = fn(
    ark_mem: &mut ARKodeMem,
    argidx: &mut usize,
    argv: &[String],
    offset: usize,
    arg_used: &mut bool,
) -> i32;
pub type ARKTimestepGetStageIndex =
    fn(ark_mem: &mut ARKodeMem, stage: &mut i32, max_stages: &mut i32) -> i32;

/* time stepper interface functions -- temporal adaptivity */
pub type ARKTimestepComputeH0 = fn(ark_mem: &mut ARKodeMem, tout: f64, hin: &mut f64) -> i32;
pub type ARKTimestepGetEstLocalErrors = fn(ark_mem: &mut ARKodeMem, ele: &mut NVector) -> i32;
pub type ARKSetAdaptControllerFn =
    fn(ark_mem: &mut ARKodeMem, c: Option<SUNAdaptController>) -> i32;

/* time stepper interface functions -- relaxation */
pub type ARKTimestepSetRelaxFn =
    fn(ark_mem: &mut ARKodeMem, rfn: Option<ARKRelaxFn>, rjac: Option<ARKRelaxJacFn>) -> i32;

/* time stepper interface functions -- implicit solvers */
pub type ARKTimestepAttachLinsolFn = fn(
    ark_mem: &mut ARKodeMem,
    linit: Option<ARKLinsolInitFn>,
    lsetup: Option<ARKLinsolSetupFn>,
    lsolve: Option<ARKLinsolSolveFn>,
    lfree: Option<ARKLinsolFreeFn>,
    lsolve_type: SUNLinearSolver_Type,
    lmem: Box<crate::arkode_ls_impl::ARKLsMem>,
) -> i32;
pub type ARKTimestepDisableLSetup = fn(ark_mem: &mut ARKodeMem);
/// C returns `void*` (the stepper's lin-solver memory).  Storage
/// adaptation (ARCHITECTURE.md Addendum C.1): the box itself lives on
/// `ARKodeMem.lmem`; this op TAKES it out (take/put-back convention)
/// and the caller restores it by writing the field.
pub type ARKTimestepGetLinMemFn =
    fn(ark_mem: &mut ARKodeMem) -> Option<Box<crate::arkode_ls_impl::ARKLsMem>>;
/// Rust-only companion op: C's step_getgammas hands out the address of
/// the stepper's `jcur` flag, through which the ARKLS preconditioner
/// setup writes; the Rust op copies it out, so this setter carries the
/// write-back into the stepper.
pub type ARKTimestepSetJcurFn = fn(ark_mem: &mut ARKodeMem, jcur: bool);
pub type ARKTimestepGetImplicitRHSFn = fn(ark_mem: &mut ARKodeMem) -> Option<ARKRhsFn>;
pub type ARKTimestepGetGammasFn = fn(
    ark_mem: &mut ARKodeMem,
    gamma: &mut f64,
    gamrat: &mut f64,
    jcur: &mut bool,
    dgamma_fail: &mut bool,
) -> i32;
pub type ARKTimestepComputeState =
    fn(ark_mem: &mut ARKodeMem, zcor: &NVector, z: &mut NVector) -> i32;
pub type ARKTimestepSetNonlinearSolver =
    fn(ark_mem: &mut ARKodeMem, nls: NonlinearSolver) -> i32;
pub type ARKTimestepSetLinear = fn(ark_mem: &mut ARKodeMem, timedepend: i32) -> i32;
pub type ARKTimestepSetNonlinear = fn(ark_mem: &mut ARKodeMem) -> i32;
pub type ARKTimestepSetAutonomous = fn(ark_mem: &mut ARKodeMem, autonomous: bool) -> i32;
/// C allows a NULL nls_fi (reset to the stepper's fi).
pub type ARKTimestepSetNlsRhsFn =
    fn(ark_mem: &mut ARKodeMem, nls_fi: Option<ARKRhsFn>) -> i32;
pub type ARKTimestepSetDeduceImplicitRhs = fn(ark_mem: &mut ARKodeMem, deduce: bool) -> i32;
pub type ARKTimestepSetNonlinCRDown = fn(ark_mem: &mut ARKodeMem, crdown: f64) -> i32;
pub type ARKTimestepSetNonlinRDiv = fn(ark_mem: &mut ARKodeMem, rdiv: f64) -> i32;
pub type ARKTimestepSetDeltaGammaMax = fn(ark_mem: &mut ARKodeMem, dgmax: f64) -> i32;
pub type ARKTimestepSetLSetupFrequency = fn(ark_mem: &mut ARKodeMem, msbp: i32) -> i32;
pub type ARKTimestepSetPredictorMethod = fn(ark_mem: &mut ARKodeMem, method: i32) -> i32;
pub type ARKTimestepSetMaxNonlinIters = fn(ark_mem: &mut ARKodeMem, maxcor: i32) -> i32;
pub type ARKTimestepSetNonlinConvCoef = fn(ark_mem: &mut ARKodeMem, nlscoef: f64) -> i32;
pub type ARKTimestepSetStagePredictFn =
    fn(ark_mem: &mut ARKodeMem, predict_stage: Option<ARKStagePredictFn>) -> i32;
pub type ARKTimestepGetNumLinSolvSetups =
    fn(ark_mem: &mut ARKodeMem, nlinsetups: &mut i64) -> i32;
pub type ARKTimestepGetCurrentGamma = fn(ark_mem: &mut ARKodeMem, gamma: &mut f64) -> i32;
/// C fills pointers to the stepper's internal vectors (zpred, z, Fi,
/// sdata) plus tcur/gamma/user_data; safe Rust cannot lend five
/// aliased &muts — the vector out-parameters receive CLONES instead
/// (user_data omitted: it stays on ark_mem).
pub type ARKTimestepGetNonlinearSystemData = fn(
    ark_mem: &mut ARKodeMem,
    tcur: &mut f64,
    zpred: &mut NVector,
    z: &mut NVector,
    Fi: &mut NVector,
    gamma: &mut f64,
    sdata: &mut NVector,
) -> i32;
pub type ARKTimestepGetNumNonlinSolvIters = fn(ark_mem: &mut ARKodeMem, nniters: &mut i64) -> i32;
pub type ARKTimestepGetNumNonlinSolvConvFails =
    fn(ark_mem: &mut ARKodeMem, nnfails: &mut i64) -> i32;
pub type ARKTimestepGetNonlinSolvStats =
    fn(ark_mem: &mut ARKodeMem, nniters: &mut i64, nnfails: &mut i64) -> i32;

/* time stepper interface functions -- non-identity mass matrices */
pub type ARKTimestepAttachMasssolFn = fn(
    ark_mem: &mut ARKodeMem,
    minit: Option<ARKMassInitFn>,
    msetup: Option<ARKMassSetupFn>,
    mmult: Option<ARKMassMultFn>,
    msolve: Option<ARKMassSolveFn>,
    mfree: Option<ARKMassFreeFn>,
    time_dep: bool,
    msolve_type: SUNLinearSolver_Type,
    mass_mem: Box<crate::arkode_ls_impl::ARKLsMassMem>,
) -> i32;
pub type ARKTimestepDisableMSetup = fn(ark_mem: &mut ARKodeMem);
/// C returns `void*` (the stepper's mass-solver memory).
/// C returns `void*` (the stepper's mass-solver memory); like lmem
/// (Addendum C.2) the box lives on ARKodeMem.mass_mem, this op TAKES
/// it out, and put-back writes the field.
pub type ARKTimestepGetMassMemFn =
    fn(ark_mem: &mut ARKodeMem) -> Option<Box<crate::arkode_ls_impl::ARKLsMassMem>>;

/* time stepper interface functions -- forcing */
pub type ARKTimestepSetForcingFn =
    fn(ark_mem: &mut ARKodeMem, tshift: f64, tscale: f64, f: &[NVector], nvecs: i32) -> i32;

/*===============================================================
  High level error handler, used throughout ARKODE
  ===============================================================*/

/// C arkProcessError routes printf-style messages to the SUNContext
/// error handler stack; here messages go to stderr (equivalent to
/// the default SUNLogErrHandlerFn behavior), matching the
/// cvProcessError / IDAProcessError ports.
pub fn arkProcessError(
    _ark_mem: Option<&ARKodeMem>,
    error_code: i32,
    line: u32,
    func: &str,
    file: &str,
    msg: &str,
) {
    if error_code == ARK_WARNING {
        eprintln!("\n[ARKODE WARNING] {file}:{line} in {func}\n  {msg}\n");
    } else {
        eprintln!("\n[ARKODE ERROR] {file}:{line} in {func}\n  {msg}\n");
    }
}

/*===============================================================
  ARKODE interpolation module definition
  ===============================================================*/

/// The generic ARKInterp (C: void* content + ops table) as enum
/// dispatch over the two implementations; the operations
/// (arkInterpResize/Free/SetDegree/Init/Update/Evaluate dispatchers
/// and the _Hermite/_Lagrange implementations) live in
/// arkode_interp.rs.
pub enum ARKInterp {
    Hermite(ARKInterpContent_Hermite),
    Lagrange(ARKInterpContent_Lagrange),
}

/*===============================================================
  ARKODE data structures
  ===============================================================*/

/// struct ARKodeMassMemRec (arkode_impl.h): data pertaining to the
/// use of a non-identity mass matrix.
pub struct ARKodeMassMem {
    /* mass matrix linear solver interface function pointers */
    pub minit: Option<ARKMassInitFn>,
    pub msetup: Option<ARKMassSetupFn>,
    pub mmult: Option<ARKMassMultFn>,
    pub msolve: Option<ARKMassSolveFn>,
    pub mfree: Option<ARKMassFreeFn>,
    pub sol_mem: UserData, /* mass matrix solver interface data */
    pub msolve_type: i32,  /* mass matrix interface type:
                              0=iterative; 1=direct; 2=custom */
}

/// struct ARKodeMemRec (arkode_impl.h): the main ARKODE memory.
pub struct ARKodeMem {
    pub uround: f64, /* machine unit roundoff */

    /* Problem specification data */
    pub user_data: UserData, /* user ptr passed to supplied functions */
    pub itol: i32,           /* ARK_SS (scalar, default), ARK_SV (vector),
                                ARK_WF (user weight function) */
    pub ritol: i32,          /* same options for residual tolerances */
    pub reltol: f64,         /* relative tolerance                    */
    pub Sabstol: f64,        /* scalar absolute solution tolerance    */
    pub Vabstol: Option<NVector>, /* vector absolute solution tolerance */
    pub atolmin0: bool,      /* flag indicating that min(abstol) = 0  */
    pub SRabstol: f64,       /* scalar absolute residual tolerance    */
    pub VRabstol: Option<NVector>, /* vector absolute residual tolerance */
    pub Ratolmin0: bool,     /* flag indicating that min(Rabstol) = 0 */
    pub user_efun: bool,     /* SUNTRUE if user sets efun             */
    pub efun: Option<ARKEwtFn>, /* function to set ewt                */
    pub e_data: UserData,    /* user pointer passed to efun           */
    pub user_rfun: bool,     /* SUNTRUE if user sets rfun             */
    pub rfun: Option<ARKRwtFn>, /* function to set rwt                */
    pub r_data: UserData,    /* user pointer passed to rfun           */

    /* Time stepper module -- general */
    pub step_mem: UserData,
    pub step_init: Option<ARKTimestepInitFn>,
    pub step_fullrhs: Option<ARKTimestepFullRHSFn>,
    pub step: Option<ARKTimestepStepFn>,
    pub step_setuserdata: Option<ARKTimetepSetUserDataFn>,
    pub step_printallstats: Option<ARKTimestepPrintAllStats>,
    pub step_writeparameters: Option<ARKTimestepWriteParameters>,
    pub step_resize: Option<ARKTimestepResize>,
    pub step_reset: Option<ARKTimestepReset>,
    pub step_free: Option<ARKTimestepFree>,
    pub step_printmem: Option<ARKTimestepPrintMem>,
    pub step_setdefaults: Option<ARKTimestepSetDefaults>,
    pub step_setorder: Option<ARKTimestepSetOrder>,
    pub step_getnumrhsevals: Option<ARKTimestepGetNumRhsEvals>,
    pub step_setstepdirection: Option<ARKTimestepSetStepDirection>,
    pub step_setusecompensatedsums: Option<ARKTimestepSetUseCompensatedSums>,
    pub step_setoptions: Option<ARKTimestepSetOptions>,
    pub step_getstageindex: Option<ARKTimestepGetStageIndex>,

    /* Time stepper module -- temporal adaptivity */
    pub step_supports_adaptive: bool,
    pub step_H0: Option<ARKTimestepComputeH0>,
    pub step_setadaptcontroller: Option<ARKSetAdaptControllerFn>,
    pub step_getestlocalerrors: Option<ARKTimestepGetEstLocalErrors>,

    /* Time stepper module -- relaxation */
    pub step_supports_relaxation: bool,
    pub step_setrelaxfn: Option<ARKTimestepSetRelaxFn>,

    /* Time stepper module -- implicit solvers */
    pub step_supports_implicit: bool,
    pub step_attachlinsol: Option<ARKTimestepAttachLinsolFn>,
    pub step_disablelsetup: Option<ARKTimestepDisableLSetup>,
    pub step_getlinmem: Option<ARKTimestepGetLinMemFn>,
    pub step_setjcur: Option<ARKTimestepSetJcurFn>,
    pub step_getimplicitrhs: Option<ARKTimestepGetImplicitRHSFn>,
    pub step_getgammas: Option<ARKTimestepGetGammasFn>,
    /* Storage adaptation (Addendum C.1): C keeps the ARKLS interface
       memory behind the stepper's `void* lmem`; the box lives here so
       ARKLS routines avoid a double-nested take/put-back.  The
       step_getlinmem op above takes it out; put-back writes this
       field.  None = C NULL. */
    pub lmem: Option<Box<crate::arkode_ls_impl::ARKLsMem>>,
    /* Same hoisting for the ARKLS mass-solver memory (C: the
       stepper's `void* mass_mem`). */
    pub mass_mem: Option<Box<crate::arkode_ls_impl::ARKLsMassMem>>,
    pub step_computestate: Option<ARKTimestepComputeState>,
    pub step_setnonlinearsolver: Option<ARKTimestepSetNonlinearSolver>,
    pub step_setlinear: Option<ARKTimestepSetLinear>,
    pub step_setautonomous: Option<ARKTimestepSetAutonomous>,
    pub step_setnonlinear: Option<ARKTimestepSetNonlinear>,
    pub step_setnlsrhsfn: Option<ARKTimestepSetNlsRhsFn>,
    pub step_setdeduceimplicitrhs: Option<ARKTimestepSetDeduceImplicitRhs>,
    pub step_setnonlincrdown: Option<ARKTimestepSetNonlinCRDown>,
    pub step_setnonlinrdiv: Option<ARKTimestepSetNonlinRDiv>,
    pub step_setdeltagammamax: Option<ARKTimestepSetDeltaGammaMax>,
    pub step_setlsetupfrequency: Option<ARKTimestepSetLSetupFrequency>,
    pub step_setpredictormethod: Option<ARKTimestepSetPredictorMethod>,
    pub step_setmaxnonliniters: Option<ARKTimestepSetMaxNonlinIters>,
    pub step_setnonlinconvcoef: Option<ARKTimestepSetNonlinConvCoef>,
    pub step_setstagepredictfn: Option<ARKTimestepSetStagePredictFn>,
    pub step_getnumlinsolvsetups: Option<ARKTimestepGetNumLinSolvSetups>,
    pub step_getcurrentgamma: Option<ARKTimestepGetCurrentGamma>,
    pub step_getnonlinearsystemdata: Option<ARKTimestepGetNonlinearSystemData>,
    pub step_getnumnonlinsolviters: Option<ARKTimestepGetNumNonlinSolvIters>,
    pub step_getnumnonlinsolvconvfails: Option<ARKTimestepGetNumNonlinSolvConvFails>,
    pub step_getnonlinsolvstats: Option<ARKTimestepGetNonlinSolvStats>,

    /* Time stepper module -- non-identity mass matrices */
    pub step_supports_massmatrix: bool,
    pub step_attachmasssol: Option<ARKTimestepAttachMasssolFn>,
    pub step_disablemsetup: Option<ARKTimestepDisableMSetup>,
    pub step_getmassmem: Option<ARKTimestepGetMassMemFn>,
    pub step_mmult: Option<ARKMassMultFn>,

    /* Time stepper module -- forcing */
    pub step_setforcing: Option<ARKTimestepSetForcingFn>,

    /* N_Vector storage */
    pub ewt: NVector,        /* error weight vector                      */
    pub rwt: NVector,        /* residual weight vector                   */
    pub rwt_is_ewt: bool,    /* SUNTRUE if rwt is a pointer to ewt       */
    pub ycur: NVector,       /* C: pointer to user-provided solution
                                memory; here owned, with copy-back at
                                every ARKodeEvolve return path            */
    pub ensure_ycur: bool,   /* SUNTRUE if stepper expects ycur=yn on
                                entry to its takestep routine             */
    pub yn: NVector,         /* solution from the last successful step   */
    pub fn_: NVector,        /* C `fn`: full IVP right-hand side from
                                last step (renamed: fn is reserved)       */
    pub fn_is_current: bool, /* SUNTRUE if fn has been evaluated at yn   */
    pub tempv1: NVector,     /* temporary storage vectors (for local use */
    pub tempv2: NVector,     /* and by time-stepping modules)            */
    pub tempv3: NVector,
    pub tempv4: NVector,
    pub tempv5: NVector,

    /* Temporal interpolation module */
    pub interp: Option<ARKInterp>,
    pub interp_type: i32,
    pub interp_degree: i32,

    /* Tstop information */
    pub tstopset: bool,
    pub tstopinterp: bool,
    pub tstop: f64,

    /* Time step data */
    pub hin: f64,        /* initial step size                        */
    pub h: f64,          /* current step size                        */
    pub hmin: f64,       /* |h| >= hmin                              */
    pub hmax_inv: f64,   /* |h| <= 1/hmax_inv                        */
    pub hprime: f64,     /* next actual step size to be used         */
    pub next_h: f64,     /* next dynamical step size (only used in
                            getCurrentStep); note that this could
                            overtake tstop                           */
    pub eta: f64,        /* eta = hprime / h                         */
    pub tcur: f64,       /* current internal value of t
                            (changes with each stage)                */
    pub tretlast: f64,   /* value of tret last returned by ARKODE    */
    pub fixedstep: bool, /* flag to disable temporal adaptivity      */
    pub hadapt_mem: Option<ARKodeHAdaptMem>, /* time step adaptivity structure */

    /* Limits and various solver parameters */
    pub mxstep: i64, /* max number of internal steps for one user call */
    pub mxhnil: i32, /* max number of warning messages issued to the
                        user that t+h == t for the next internal step  */
    pub maxnef: i32, /* max error test fails in one step               */
    pub maxncf: i32, /* max num alg. solver conv. fails in one step    */

    /* Counters */
    pub nst_attempts: i64, /* number of attempted steps                  */
    pub nst: i64,          /* number of internal steps taken             */
    pub nhnil: i32,        /* number of messages issued to the user that
                              t+h == t for the next internal step        */
    pub ncfn: i64,         /* num corrector convergence failures         */
    pub netf: i64,         /* num error test failures                    */

    /* Space requirements for ARKODE */
    pub lrw1: i64, /* no. of sunrealtype words in 1 N_Vector          */
    pub liw1: i64, /* no. of integer words in 1 N_Vector              */
    pub lrw: i64,  /* no. of sunrealtype words in ARKODE work vectors */
    pub liw: i64,  /* no. of integer words in ARKODE work vectors     */

    /* Saved Values */
    pub h0u: f64,   /* actual initial stepsize                     */
    pub tn: f64,    /* time of last successful step                */
    pub terr: f64,  /* error in tn for compensated sums            */
    pub hold: f64,  /* last successful h value used                */
    pub tolsf: f64, /* tolerance scale factor (suggestion to user) */
    pub AccumErrorType: ARKAccumError, /* accumulated error estimation type */
    pub AccumErrorStart: f64, /* time of last accumulated error reset */
    pub AccumError: f64, /* accumulated error estimate               */
    pub VabstolMallocDone: bool,
    pub VRabstolMallocDone: bool,
    pub MallocDone: bool,
    pub initsetup: bool,    /* denotes a call to InitialSetup is needed  */
    pub init_type: i32,     /* initialization type (see constants above) */
    pub firststage: bool,   /* denotes first stage in simulation         */
    pub initialized: bool,  /* denotes arkInitialSetup has been done     */
    pub call_fullrhs: bool, /* denotes the full RHS fn will be called    */
    pub preallocated: bool, /* SUNTRUE if ARKodeInit has been called to
                               preallocate data prior to ARKodeEvolve    */

    /* Rootfinding Data */
    pub root_mem: Option<ARKodeRootMem>, /* root-finding structure */

    /* Inequality Constraints Data */
    pub constraints: Option<NVector>, /* vector of constraint flags     */
    pub nconstrfails: i64,            /* total constraint failures      */
    pub maxconstrfails: i32,          /* max failures allowed in a step */

    /* Relaxation Data */
    pub relax_enabled: bool,               /* is relaxation enabled?    */
    pub relax_mem: Option<ARKodeRelaxMem>, /* relaxation data structure */

    /* User-supplied step solution pre/post-processing functions */
    pub PreStepFn: Option<ARKPreStepFn>,
    pub PostStepFn: Option<ARKPostStepFn>,

    /* User-supplied RHS function pre-processing function */
    pub PreRhsFn: Option<ARKPreRhsFn>,

    /* User-supplied stage and step solution post-processing function */
    pub PostProcessStepFn: Option<ARKPostProcessFn>,
    pub PostProcessStageFn: Option<ARKPostProcessFn>,

    pub use_compensated_sums: bool,

    /* Adjoint solver data */
    pub load_checkpoint_fail: bool,
    pub do_adjoint: bool,
    pub adj_stage_idx: suncountertype, /* current stage index (only valid in adjoint context) */
    pub adj_step_idx: suncountertype,  /* current step index (only valid in adjoint context)  */

    /* Checkpointing data */
    pub checkpoint_scheme: Option<SUNAdjointCheckpointScheme>,
    pub checkpoint_step_idx: suncountertype, /* the step number for checkpointing */

    /* XBraid interface variables (XBraid itself is an excluded
    backend; these flags are read by the step attempt loop) */
    pub force_pass: bool, /* when true the step attempt loop will ignore the
                             return value (kflag) from arkCheckTemporalError
                             and set kflag = ARK_SUCCESS to force the step
                             attempt to always pass (if a solver failure did
                             not occur before the error test). */
    pub last_kflag: i32,  /* last value of the return flag (kflag) from a call
                             to arkCheckTemporalError. This is only set when
                             force_pass is true. */
}

/// C `memset(ark_mem, 0, sizeof(struct ARKodeMemRec))` at the top of
/// arkCreate: every field zero/NULL/false (enums take their 0 value).
/// arkCreate (arkode.rs) builds on this.
impl Default for ARKodeMem {
    fn default() -> Self {
        ARKodeMem {
            uround: 0.0,
            user_data: None,
            itol: 0,
            ritol: 0,
            reltol: 0.0,
            Sabstol: 0.0,
            Vabstol: None,
            atolmin0: false,
            SRabstol: 0.0,
            VRabstol: None,
            Ratolmin0: false,
            user_efun: false,
            efun: None,
            e_data: None,
            user_rfun: false,
            rfun: None,
            r_data: None,
            step_mem: None,
            step_init: None,
            step_fullrhs: None,
            step: None,
            step_setuserdata: None,
            step_printallstats: None,
            step_writeparameters: None,
            step_resize: None,
            step_reset: None,
            step_free: None,
            step_printmem: None,
            step_setdefaults: None,
            step_setorder: None,
            step_getnumrhsevals: None,
            step_setstepdirection: None,
            step_setusecompensatedsums: None,
            step_setoptions: None,
            step_getstageindex: None,
            step_supports_adaptive: false,
            step_H0: None,
            step_setadaptcontroller: None,
            step_getestlocalerrors: None,
            step_supports_relaxation: false,
            step_setrelaxfn: None,
            step_supports_implicit: false,
            step_attachlinsol: None,
            step_disablelsetup: None,
            step_getlinmem: None,
            step_setjcur: None,
            step_getimplicitrhs: None,
            step_getgammas: None,
            lmem: None,
            mass_mem: None,
            step_computestate: None,
            step_setnonlinearsolver: None,
            step_setlinear: None,
            step_setautonomous: None,
            step_setnonlinear: None,
            step_setnlsrhsfn: None,
            step_setdeduceimplicitrhs: None,
            step_setnonlincrdown: None,
            step_setnonlinrdiv: None,
            step_setdeltagammamax: None,
            step_setlsetupfrequency: None,
            step_setpredictormethod: None,
            step_setmaxnonliniters: None,
            step_setnonlinconvcoef: None,
            step_setstagepredictfn: None,
            step_getnumlinsolvsetups: None,
            step_getcurrentgamma: None,
            step_getnonlinearsystemdata: None,
            step_getnumnonlinsolviters: None,
            step_getnumnonlinsolvconvfails: None,
            step_getnonlinsolvstats: None,
            step_supports_massmatrix: false,
            step_attachmasssol: None,
            step_disablemsetup: None,
            step_getmassmem: None,
            step_mmult: None,
            step_setforcing: None,
            ewt: NVector::new(0),
            rwt: NVector::new(0),
            rwt_is_ewt: false,
            ycur: NVector::new(0),
            ensure_ycur: false,
            yn: NVector::new(0),
            fn_: NVector::new(0),
            fn_is_current: false,
            tempv1: NVector::new(0),
            tempv2: NVector::new(0),
            tempv3: NVector::new(0),
            tempv4: NVector::new(0),
            tempv5: NVector::new(0),
            interp: None,
            interp_type: 0,
            interp_degree: 0,
            tstopset: false,
            tstopinterp: false,
            tstop: 0.0,
            hin: 0.0,
            h: 0.0,
            hmin: 0.0,
            hmax_inv: 0.0,
            hprime: 0.0,
            next_h: 0.0,
            eta: 0.0,
            tcur: 0.0,
            tretlast: 0.0,
            fixedstep: false,
            hadapt_mem: None,
            mxstep: 0,
            mxhnil: 0,
            maxnef: 0,
            maxncf: 0,
            nst_attempts: 0,
            nst: 0,
            nhnil: 0,
            ncfn: 0,
            netf: 0,
            lrw1: 0,
            liw1: 0,
            lrw: 0,
            liw: 0,
            h0u: 0.0,
            tn: 0.0,
            terr: 0.0,
            hold: 0.0,
            tolsf: 0.0,
            AccumErrorType: ARK_ACCUMERROR_NONE,
            AccumErrorStart: 0.0,
            AccumError: 0.0,
            VabstolMallocDone: false,
            VRabstolMallocDone: false,
            MallocDone: false,
            initsetup: false,
            init_type: 0,
            firststage: false,
            initialized: false,
            call_fullrhs: false,
            preallocated: false,
            root_mem: None,
            constraints: None,
            nconstrfails: 0,
            maxconstrfails: 0,
            relax_enabled: false,
            relax_mem: None,
            PreStepFn: None,
            PostStepFn: None,
            PreRhsFn: None,
            PostProcessStepFn: None,
            PostProcessStageFn: None,
            use_compensated_sums: false,
            load_checkpoint_fail: false,
            do_adjoint: false,
            adj_stage_idx: 0,
            adj_step_idx: 0,
            checkpoint_scheme: None,
            checkpoint_step_idx: 0,
            force_pass: false,
            last_kflag: 0,
        }
    }
}

/*===============================================================
  Reusable ARKODE Error Messages (arkode_impl.h; MSG_TIME* formats
  are produced with fmt_g at the call sites)
  ===============================================================*/

pub const MSG_ARK_NO_MEM: &str = "arkode_mem = NULL illegal.";
pub const MSG_ARK_ARKMEM_FAIL: &str = "Allocation of arkode_mem failed.";
pub const MSG_ARK_MEM_FAIL: &str = "A memory request failed.";
pub const MSG_ARK_NO_MALLOC: &str = "Attempt to call before ARKODE initialized.";
pub const MSG_ARK_NULL_Y0: &str = "y0 = NULL illegal.";
pub const MSG_ARK_NULL_F: &str = "Must specify at least one of fe, fi (both NULL).";
pub const MSG_ARK_NULL_G: &str = "g = NULL illegal.";
pub const MSG_ARK_NULL_DKY: &str = "dky = NULL illegal.";
pub const MSG_ARKADAPT_NO_MEM: &str = "Adaptivity memory structure not allocated.";
pub const MSG_ARK_MISSING_FULLRHS: &str = "Time-stepping module missing fullrhs routine \
                                           (required by requested solver configuration).";

/// Pinned crate convention for C calls of the form
/// `ark_mem->step_fullrhs(ark_mem, t, ark_mem->yn, ark_mem->fn, mode)`:
/// the ark_mem-owned argument vectors are taken out for the call and
/// put back after. step_fullrhs implementations never read
/// ark_mem.yn / ark_mem.fn_ directly (guaranteed by the C code: all
/// state arrives through the y/f arguments).
pub fn ark_step_fullrhs_yn_fn(ark_mem: &mut ARKodeMem, t: f64, mode: i32) -> i32 {
    let fullrhs = ark_mem.step_fullrhs.unwrap();
    let yn = std::mem::replace(&mut ark_mem.yn, NVector::new(0));
    let mut f = std::mem::replace(&mut ark_mem.fn_, NVector::new(0));
    let retval = fullrhs(ark_mem, t, &yn, &mut f, mode);
    ark_mem.yn = yn;
    ark_mem.fn_ = f;
    retval
}
