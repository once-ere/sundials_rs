/* -----------------------------------------------------------------
 * Translated from src/idas/idas.c (IDAS 7.7.0) — PART 1 IN PROGRESS.
 * Main IDAS integrator for DAE systems F(t, y, y') = 0 with
 * quadratures, forward sensitivities and quadrature sensitivities:
 * creation/initialization/re-initialization and tolerances.
 * Structural donor: ida_rs/src/ida.rs (verified Phase 4); functions
 * whose C text is identical between ida.c and idas.c are carried
 * verbatim from the donor.
 *
 * STATUS: NOT yet registered in lib.rs (no `pub mod idas;`) and NOT
 * compile-verified — this slice references symbols that land with
 * later sections of this file (IDACheckNvector, IDAAllocVectors,
 * IDAQuadAllocVectors) and with later Phase 5 units
 * (crate::idas_nls::IDASetNonlinearSolver / idaNlsInit,
 * crate::idas_ls::idaLsInitializeCounters).  Register the module and
 * run `cargo build -p idas_rs` only after those symbols exist.
 * Conventions follow the donor cvode_rs/src/cvode.rs.
 * -----------------------------------------------------------------*/
use crate::idas_impl::*;
use crate::idas_ls::{idaLsInitialize, idaLsInitializeCounters, idaLsPerf};
use crate::idas_nls::{idaNlsInit, IDASetNonlinearSolver};
use crate::idas_nls_sim::{idaNlsInitSensSim, IDASetNonlinearSolverSensSim};
use crate::idas_nls_stg::{idaNlsInitSensStg, IDASetNonlinearSolverSensStg};
use crate::nvector_serial::*;
use crate::sundials_context::SUNContext;
use crate::sundials_math::SUNRabs;
use crate::sundials_types::*;
use crate::sundials_utils::fmt_g;
use crate::sunnonlinsol_newton::{SUNNonlinSol_Newton, SUNNonlinSol_NewtonSens};

/*
 * =================================================================
 * IDAS PRIVATE CONSTANTS
 * =================================================================
 */

const ZERO: f64 = 0.0; /* real 0.0    */
const HALF: f64 = 0.5; /* real 0.5    */
/* TWOTHIRDS (0.667) is only used by IDACreate's steptol default,
   which lives in the idas_impl::IDAMem Default impl. */
const ONE: f64 = 1.0; /* real 1.0    */
pub const ONEPT5: f64 = 1.5; /* real 1.5    */
const TWO: f64 = 2.0; /* real 2.0    */
const FOUR: f64 = 4.0; /* real 4.0    */
pub const FIVE: f64 = 5.0; /* real 5.0    */
const TEN: f64 = 10.0; /* real 10.0   */
pub const TWENTY: f64 = 20.0; /* real 20.0   */
const HUNDRED: f64 = 100.0; /* real 100.0  */
pub const PT9: f64 = 0.9; /* real 0.9    */
pub const PT1: f64 = 0.1; /* real 0.1    */
pub const PT01: f64 = 0.01; /* real 0.01   */
const PT001: f64 = 0.001; /* real 0.001  */
const PT0001: f64 = 0.0001; /* real 0.0001 */

/* real 1 + epsilon used in testing if the step size is below its bound */
pub const ONEPSM: f64 = 1.000001;

/*
 * =================================================================
 * IDAS ROUTINE-SPECIFIC CONSTANTS
 * =================================================================
 */

/* IDAStep control constants */

pub const PREDICT_AGAIN: i32 = 20;

/* Return values for lower level routines used by IDASolve */

pub const CONTINUE_STEPS: i32 = 99;

/* IDACompleteStep constants */

pub const UNSET: i32 = -1;
pub const LOWER: i32 = 1;
pub const RAISE: i32 = 2;
pub const MAINTAIN: i32 = 3;

/* IDATestError constants */

pub const ERROR_TEST_FAIL: i32 = 7;

/* Control constants for lower-level rootfinding functions */

pub const RTFOUND: i32 = 1;
pub const CLOSERT: i32 = 3;

/* Control constants for sensitivity DQ */

pub const CENTERED1: i32 = 1;
pub const CENTERED2: i32 = 2;
pub const FORWARD1: i32 = 3;
pub const FORWARD2: i32 = 4;

/* (The itol control constants IDA_NN/IDA_SS/IDA_SV/IDA_WF/IDA_EE and
   the algorithmic constants MXNCF/MXNEF/MAXNH/MAXNJ/MAXNI/EPCON/
   MAXBACKS are defined in idas_impl.rs, where the IDAMem Default
   impl uses them.) */

/*
 * -----------------------------------------------------------------
 * Message rendering helper
 *
 * The C IDAProcessError is printf-style; the idas_impl.rs message
 * constants keep the C formats verbatim with SUN_FORMAT_G expanded
 * to "%.15g".  This helper substitutes the placeholders in order,
 * rendering each value with sundials_utils::fmt_g (C's %.15g).
 * -----------------------------------------------------------------
 */
pub(crate) fn ida_msg_g(fmt: &str, vals: &[f64]) -> String {
    let mut out = String::from(fmt);
    for v in vals {
        if let Some(pos) = out.find("%.15g") {
            out.replace_range(pos..pos + 5, &fmt_g(*v, 0, 15));
        }
    }
    out
}

/*
 * =================================================================
 * EXPORTED FUNCTIONS IMPLEMENTATION
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Creation, allocation and re-initialization functions
 * -----------------------------------------------------------------
 */

/*
 * IDACreate
 *
 * IDACreate creates an internal memory block for a problem to
 * be solved by IDA.
 * If successful, IDACreate returns a pointer to the problem memory.
 * This pointer should be passed to IDAInit.
 *
 * (In C a NULL sunctx or a failed malloc return NULL; neither state
 * exists in safe Rust — the context is a unit-like struct and
 * allocation is infallible — so the error branches vanish.  The
 * default-setting body of IDACreate — including the quadrature,
 * sensitivity, quadrature-sensitivity and adjoint defaults — is the
 * idas_impl::IDAMem Default impl, mirroring how the donor's
 * CVodeCreate builds its literal.)
 */
pub fn IDACreate(sunctx: &SUNContext) -> Box<IDAMem> {
    Box::new(IDAMem {
        ida_sunctx: sunctx.clone(),
        ..IDAMem::default()
    })
}

/*-----------------------------------------------------------------*/

/*
 * IDAInit
 *
 * IDAInit allocates and initializes memory for a problem. All
 * problem specification inputs are checked for errors. If any
 * error occurs during initialization, it is reported to the
 * error handler function.
 */
pub fn IDAInit(
    ida_mem: &mut IDAMem,
    res: IDAResFn,
    t0: f64,
    yy0: &NVector,
    yp0: &NVector,
) -> i32 {
    /* Check for legal input parameters */

    if yy0.is_empty() {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAInit", file!(), MSG_Y0_NULL);
        return IDA_ILL_INPUT;
    }

    if yp0.is_empty() {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAInit", file!(), MSG_YP0_NULL);
        return IDA_ILL_INPUT;
    }

    /* (res == NULL cannot occur: IDAResFn is a plain fn pointer.) */

    /* Test if all required vector operations are implemented */

    let nvectorOK = IDACheckNvector(yy0);
    if !nvectorOK {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAInit", file!(), MSG_BAD_NVECTOR);
        return IDA_ILL_INPUT;
    }

    /* Set space requirements for one N_Vector */

    let (lrw1, liw1) = N_VSpace(yy0);
    ida_mem.ida_lrw1 = lrw1;
    ida_mem.ida_liw1 = liw1;

    /* Allocate the vectors (using yy0 as a template); infallible in
       Rust, so the MSG_MEM_FAIL branch of the C code vanishes.
       (The C fused-op work arrays ida_cvals/ida_Xvecs/ida_Zvecs are
       not stored, per the pinned alias-drop convention.) */

    IDAAllocVectors(ida_mem, yy0);

    /* Input checks complete at this point and history array allocated */

    /* Copy the input parameters into IDA memory block */
    ida_mem.ida_res = Some(res);
    ida_mem.ida_tn = t0;

    /* Initialize the phi array */
    ida_mem.ida_phi[0].data.copy_from_slice(&yy0.data);
    ida_mem.ida_phi[1].data.copy_from_slice(&yp0.data);

    /* create a Newton nonlinear solver object by default */
    let NLS = SUNNonlinSol_Newton(yy0, &ida_mem.ida_sunctx);

    /* attach the nonlinear solver to the IDA memory */
    let retval = IDASetNonlinearSolver(ida_mem, NLS);

    /* check that the nonlinear solver was successfully attached */
    if retval != IDA_SUCCESS {
        IDAProcessError(Some(ida_mem), retval, line!(), "IDAInit", file!(),
                        "Setting the nonlinear solver failed");
        return IDA_MEM_FAIL;
    }

    /* set ownership flag */
    ida_mem.ownNLS = SUNTRUE;

    /* All error checking is complete at this point */

    /* Set the linear solver addresses to NULL */

    ida_mem.ida_lmem = LsModule::None;

    /* Set forceSetup to SUNFALSE */

    ida_mem.ida_forceSetup = SUNFALSE;

    /* Initialize all the counters and other optional output values */

    ida_mem.ida_nst = 0;
    ida_mem.ida_nre = 0;
    ida_mem.ida_ncfn = 0;
    ida_mem.ida_netf = 0;
    ida_mem.ida_nni = 0;
    ida_mem.ida_nnf = 0;
    ida_mem.ida_nsetups = 0;

    ida_mem.ida_kused = 0;
    ida_mem.ida_hused = ZERO;
    ida_mem.ida_tolsf = ONE;

    ida_mem.ida_nge = 0;

    ida_mem.ida_irfnd = 0;

    /* Initialize counters specific to IC calculation. */
    ida_mem.ida_nbacktr = 0;

    /* Initialize root-finding variables */

    ida_mem.ida_glo = Vec::new();
    ida_mem.ida_ghi = Vec::new();
    ida_mem.ida_grout = Vec::new();
    ida_mem.ida_iroots = Vec::new();
    ida_mem.ida_rootdir = Vec::new();
    ida_mem.ida_gfun = None;
    ida_mem.ida_nrtfn = 0;
    ida_mem.ida_gactive = Vec::new();
    ida_mem.ida_mxgnull = 1;

    /* Initial setup not done yet */

    ida_mem.ida_SetupDone = SUNFALSE;

    /* Problem memory has been successfully allocated */

    ida_mem.ida_MallocDone = SUNTRUE;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

/*
 * IDAReInit
 *
 * IDAReInit re-initializes IDA's memory for a problem, assuming
 * it has already been allocated in a prior IDAInit call.
 * All problem specification inputs are checked for errors.
 * The problem size Neq is assumed to be unchanged since the call
 * to IDAInit, and the maximum order maxord must not be larger.
 * If any error occurs during reinitialization, it is reported to
 * the error handler function.
 * The return value is IDA_SUCCESS = 0 if no errors occurred, or
 * a negative value otherwise.
 */
pub fn IDAReInit(ida_mem: &mut IDAMem, t0: f64, yy0: &NVector, yp0: &NVector) -> i32 {
    /* Check if problem was malloc'ed */

    if !ida_mem.ida_MallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_MALLOC, line!(), "IDAReInit", file!(), MSG_NO_MALLOC);
        return IDA_NO_MALLOC;
    }

    /* Check for legal input parameters */

    if yy0.is_empty() {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAReInit", file!(), MSG_Y0_NULL);
        return IDA_ILL_INPUT;
    }

    if yp0.is_empty() {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAReInit", file!(), MSG_YP0_NULL);
        return IDA_ILL_INPUT;
    }

    /* Copy the input parameters into IDA memory block */

    ida_mem.ida_tn = t0;

    /* Set forceSetup to SUNFALSE */

    ida_mem.ida_forceSetup = SUNFALSE;

    /* Initialize the phi array */

    ida_mem.ida_phi[0].data.copy_from_slice(&yy0.data);
    ida_mem.ida_phi[1].data.copy_from_slice(&yp0.data);

    /* Initialize all the counters and other optional output values */

    ida_mem.ida_nst = 0;
    ida_mem.ida_nre = 0;
    ida_mem.ida_ncfn = 0;
    ida_mem.ida_netf = 0;
    ida_mem.ida_nni = 0;
    ida_mem.ida_nnf = 0;
    ida_mem.ida_nsetups = 0;

    ida_mem.ida_kused = 0;
    ida_mem.ida_hused = ZERO;
    ida_mem.ida_tolsf = ONE;

    ida_mem.ida_nge = 0;

    ida_mem.ida_irfnd = 0;

    ida_mem.constraint_corrections = 0;
    ida_mem.constraint_fails = 0;

    if let LsModule::Ls(ls) = &mut ida_mem.ida_lmem {
        idaLsInitializeCounters(ls);
    }

    /* Initial setup not done yet */

    ida_mem.ida_SetupDone = SUNFALSE;

    /* Problem has been successfully re-initialized */

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

/*
 * IDASStolerances
 * IDASVtolerances
 * IDAWFtolerances
 *
 * These functions specify the integration tolerances. One of them
 * MUST be called before the first call to IDA.
 *
 * IDASStolerances specifies scalar relative and absolute tolerances.
 * IDASVtolerances specifies scalar relative tolerance and a vector
 *   absolute tolerance (a potentially different absolute tolerance
 *   for each vector component).
 * IDAWFtolerances specifies a user-provides function (of type IDAEwtFn)
 *   which will be called to set the error weight vector.
 */

pub fn IDASStolerances(ida_mem: &mut IDAMem, reltol: f64, abstol: f64) -> i32 {
    if !ida_mem.ida_MallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_MALLOC, line!(), "IDASStolerances", file!(), MSG_NO_MALLOC);
        return IDA_NO_MALLOC;
    }

    /* Check inputs */

    if reltol < ZERO {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASStolerances", file!(), MSG_BAD_RTOL);
        return IDA_ILL_INPUT;
    }

    if abstol < ZERO {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASStolerances", file!(), MSG_BAD_ATOL);
        return IDA_ILL_INPUT;
    }

    /* Copy tolerances into memory */

    ida_mem.ida_rtol = reltol;
    ida_mem.ida_Satol = abstol;
    ida_mem.ida_atolmin0 = abstol == ZERO;

    ida_mem.ida_itol = IDA_SS;

    /* (C sets ida_efun = IDAEwtSet with edata = ida_mem; the internal
       efun is selected by ida_user_efun == SUNFALSE here, see
       ida_efun_dispatch.) */
    ida_mem.ida_user_efun = SUNFALSE;
    ida_mem.ida_efun = None;

    IDA_SUCCESS
}

pub fn IDASVtolerances(ida_mem: &mut IDAMem, reltol: f64, abstol: &NVector) -> i32 {
    if !ida_mem.ida_MallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_MALLOC, line!(), "IDASVtolerances", file!(), MSG_NO_MALLOC);
        return IDA_NO_MALLOC;
    }

    /* Check inputs */

    if reltol < ZERO {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASVtolerances", file!(), MSG_BAD_RTOL);
        return IDA_ILL_INPUT;
    }

    let atolmin = N_VMin(abstol);
    if atolmin < ZERO {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASVtolerances", file!(), MSG_BAD_ATOL);
        return IDA_ILL_INPUT;
    }

    /* Copy tolerances into memory */

    if !ida_mem.ida_VatolMallocDone {
        ida_mem.ida_Vatol = N_VClone(&ida_mem.ida_ewt);
        ida_mem.ida_lrw += ida_mem.ida_lrw1;
        ida_mem.ida_liw += ida_mem.ida_liw1;
        ida_mem.ida_VatolMallocDone = SUNTRUE;
    }

    ida_mem.ida_rtol = reltol;
    ida_mem.ida_Vatol.data.copy_from_slice(&abstol.data);
    ida_mem.ida_atolmin0 = atolmin == ZERO;

    ida_mem.ida_itol = IDA_SV;

    ida_mem.ida_user_efun = SUNFALSE;
    ida_mem.ida_efun = None;

    IDA_SUCCESS
}

pub fn IDAWFtolerances(ida_mem: &mut IDAMem, efun: IDAEwtFn) -> i32 {
    if !ida_mem.ida_MallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_MALLOC, line!(), "IDAWFtolerances", file!(), MSG_NO_MALLOC);
        return IDA_NO_MALLOC;
    }

    ida_mem.ida_itol = IDA_WF;

    ida_mem.ida_user_efun = SUNTRUE;
    ida_mem.ida_efun = Some(efun);
    /* (C ida_edata = NULL, set to user_data in InitialSetup — not
       stored, per the idas_impl.rs drop of ida_edata.) */

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

/*
 * IDAQuadMalloc
 *
 * IDAQuadMalloc allocates and initializes quadrature related
 * memory for a problem. All problem specification inputs are
 * checked for errors. If any error occurs during initialization,
 * it is reported to the file whose file pointer is errfp.
 * The return value is IDA_SUCCESS = 0 if no errors occurred, or
 * a negative value otherwise.
 */
pub fn IDAQuadInit(ida_mem: &mut IDAMem, rhsQ: IDAQuadRhsFn, yQ0: &NVector) -> i32 {
    /* Set space requirements for one N_Vector */
    let (lrw1Q, liw1Q) = N_VSpace(yQ0);
    ida_mem.ida_lrw1Q = lrw1Q;
    ida_mem.ida_liw1Q = liw1Q;

    /* Allocate the vectors (using yQ0 as a template); infallible in
       Rust, so the MSG_MEM_FAIL branch of the C code vanishes */
    IDAQuadAllocVectors(ida_mem, yQ0);

    /* Initialize phiQ in the history array */
    ida_mem.ida_phiQ[0].data.copy_from_slice(&yQ0.data);

    /* (C: N_VConstVectorArray(maxord, ZERO, phiQ+1); the serial fused
       kernel is reproduced inline per the workspace convention, and
       its cannot-fail retval branch vanishes.) */
    for j in 1..=(ida_mem.ida_maxord as usize) {
        N_VConst(ZERO, &mut ida_mem.ida_phiQ[j]);
    }

    /* Copy the input parameters into IDAS state */
    ida_mem.ida_rhsQ = Some(rhsQ);

    /* Initialize counters */
    ida_mem.ida_nrQe = 0;
    ida_mem.ida_netfQ = 0;

    /* Quadrature integration turned ON */
    ida_mem.ida_quadr = SUNTRUE;
    ida_mem.ida_quadMallocDone = SUNTRUE;

    /* Quadrature initialization was successful */
    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

/*
 * IDAQuadReInit
 *
 * IDAQuadReInit re-initializes IDAS's quadrature related memory
 * for a problem, assuming it has already been allocated in prior
 * calls to IDAInit and IDAQuadMalloc.
 * All problem specification inputs are checked for errors.
 * If any error occurs during initialization, it is reported to the
 * file whose file pointer is errfp.
 * The return value is IDA_SUCCESS = 0 if no errors occurred, or
 * a negative value otherwise.
 */
pub fn IDAQuadReInit(ida_mem: &mut IDAMem, yQ0: &NVector) -> i32 {
    /* Check if quadrature was initialized */
    if !ida_mem.ida_quadMallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_QUAD, line!(), "IDAQuadReInit", file!(), MSG_NO_QUAD);
        return IDA_NO_QUAD;
    }

    /* Initialize phiQ in the history array */
    ida_mem.ida_phiQ[0].data.copy_from_slice(&yQ0.data);

    /* (C: N_VConstVectorArray(maxord, ZERO, phiQ+1); serial fused
       kernel reproduced inline; cannot-fail retval branch vanishes.) */
    for j in 1..=(ida_mem.ida_maxord as usize) {
        N_VConst(ZERO, &mut ida_mem.ida_phiQ[j]);
    }

    /* Initialize counters */
    ida_mem.ida_nrQe = 0;
    ida_mem.ida_netfQ = 0;

    /* Quadrature integration turned ON */
    ida_mem.ida_quadr = SUNTRUE;

    /* Quadrature re-initialization was successful */
    IDA_SUCCESS
}

/*
 * IDAQuadSStolerances
 * IDAQuadSVtolerances
 *
 *
 * These functions specify the integration tolerances for quadrature
 * variables. One of them MUST be called before the first call to
 * IDA IF error control on the quadrature variables is enabled
 * (see IDASetQuadErrCon).
 *
 * IDASStolerances specifies scalar relative and absolute tolerances.
 * IDASVtolerances specifies scalar relative tolerance and a vector
 *   absolute tolerance (a potentially different absolute tolerance
 *   for each vector component).
 */
pub fn IDAQuadSStolerances(ida_mem: &mut IDAMem, reltolQ: f64, abstolQ: f64) -> i32 {
    /* Check if quadrature was initialized */
    if !ida_mem.ida_quadMallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_QUAD, line!(), "IDAQuadSStolerances", file!(), MSG_NO_QUAD);
        return IDA_NO_QUAD;
    }

    /* Test user-supplied tolerances */
    if reltolQ < ZERO {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAQuadSStolerances", file!(), MSG_BAD_RTOLQ);
        return IDA_ILL_INPUT;
    }

    if abstolQ < ZERO {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAQuadSStolerances", file!(), MSG_BAD_ATOLQ);
        return IDA_ILL_INPUT;
    }

    /* Copy tolerances into memory */
    ida_mem.ida_itolQ = IDA_SS;

    ida_mem.ida_rtolQ = reltolQ;
    ida_mem.ida_SatolQ = abstolQ;
    ida_mem.ida_atolQmin0 = abstolQ == ZERO;

    IDA_SUCCESS
}

pub fn IDAQuadSVtolerances(ida_mem: &mut IDAMem, reltolQ: f64, abstolQ: &NVector) -> i32 {
    /* Check if quadrature was initialized */
    if !ida_mem.ida_quadMallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_QUAD, line!(), "IDAQuadSVtolerances", file!(), MSG_NO_QUAD);
        return IDA_NO_QUAD;
    }

    /* Test user-supplied tolerances */
    if reltolQ < ZERO {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAQuadSVtolerances", file!(), MSG_BAD_RTOLQ);
        return IDA_ILL_INPUT;
    }

    if abstolQ.is_empty() {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAQuadSVtolerances", file!(), MSG_NULL_ATOLQ);
        return IDA_ILL_INPUT;
    }

    let atolmin = N_VMin(abstolQ);
    if atolmin < ZERO {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAQuadSVtolerances", file!(), MSG_BAD_ATOLQ);
        return IDA_ILL_INPUT;
    }

    /* Copy tolerances into memory */
    ida_mem.ida_itolQ = IDA_SV;
    ida_mem.ida_rtolQ = reltolQ;

    /* clone the absolute tolerances vector (if necessary) */
    if SUNFALSE == ida_mem.ida_VatolQMallocDone {
        ida_mem.ida_VatolQ = N_VClone(abstolQ);
        ida_mem.ida_lrw += ida_mem.ida_lrw1Q;
        ida_mem.ida_liw += ida_mem.ida_liw1Q;
        ida_mem.ida_VatolQMallocDone = SUNTRUE;
    }

    ida_mem.ida_VatolQ.data.copy_from_slice(&abstolQ.data);
    ida_mem.ida_atolQmin0 = atolmin == ZERO;

    IDA_SUCCESS
}

/*
 * IDASenMalloc
 *
 * IDASensInit allocates and initializes sensitivity related
 * memory for a problem. All problem specification inputs are
 * checked for errors. If any error occurs during initialization,
 * it is reported to the file whose file pointer is errfp.
 * The return value is IDA_SUCCESS = 0 if no errors occurred, or
 * a negative value otherwise.
 *
 * (C fS may be NULL to select the internal DQ residual, so fS is
 * Option<IDASensResFn> here; the C NULL/DQ branch sets
 * ida_user_dataS = IDA_mem, which is not stored per the pinned
 * convention — the DQ path in this file operates on &mut IDAMem.)
 */
pub fn IDASensInit(
    ida_mem: &mut IDAMem,
    Ns: i32,
    ism: i32,
    fS: Option<IDASensResFn>,
    yS0: &[NVector],
    ypS0: &[NVector],
) -> i32 {
    /* Check if Ns is legal */
    if Ns <= 0 {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASensInit", file!(), MSG_BAD_NS);
        return IDA_ILL_INPUT;
    }
    ida_mem.ida_Ns = Ns;

    /* Check if ism is legal */
    if ism != IDA_SIMULTANEOUS && ism != IDA_STAGGERED {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASensInit", file!(), MSG_BAD_ISM);
        return IDA_ILL_INPUT;
    }
    ida_mem.ida_ism = ism;

    /* Check if yS0 and ypS0 are non-null */
    if yS0.is_empty() {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASensInit", file!(), MSG_NULL_YYS0);
        return IDA_ILL_INPUT;
    }
    if ypS0.is_empty() {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASensInit", file!(), MSG_NULL_YPS0);
        return IDA_ILL_INPUT;
    }

    /* Store sensitivity RHS-related data */

    if let Some(fS) = fS {
        ida_mem.ida_resS = Some(fS);
        ida_mem.ida_resSDQ = SUNFALSE;
    } else {
        ida_mem.ida_resS = None;
        ida_mem.ida_resSDQ = SUNTRUE;
    }

    /* Allocate the vectors (using yS0[0] as a template); infallible
       in Rust.  (The C fused-op work-array reallocation to Ns*MXORDP1
       — ida_cvals/ida_Xvecs/ida_Zvecs — is not ported, per the pinned
       alias-drop convention.) */

    IDASensAllocVectors(ida_mem, &yS0[0]);

    /*----------------------------------------------
      All error checking is complete at this point
      -----------------------------------------------*/

    /* Initialize the phiS array
       (C: cvals[is]=ONE + N_VScaleVectorArray; serial fused kernel
       reproduced inline; cannot-fail retval branches vanish.) */
    for is in 0..(Ns as usize) {
        ida_mem.ida_phiS[0][is].data.copy_from_slice(&yS0[is].data);
        ida_mem.ida_phiS[1][is].data.copy_from_slice(&ypS0[is].data);
    }

    /* Initialize all sensitivity related counters */
    ida_mem.ida_nrSe = 0;
    ida_mem.ida_nreS = 0;
    ida_mem.ida_ncfnS = 0;
    ida_mem.ida_netfS = 0;
    ida_mem.ida_nniS = 0;
    ida_mem.ida_nnfS = 0;
    ida_mem.ida_nsetupsS = 0;

    /* Set default values for plist and pbar */
    for is in 0..(Ns as usize) {
        ida_mem.ida_plist[is] = is as i32;
        ida_mem.ida_pbar[is] = ONE;
    }

    /* Sensitivities will be computed */
    ida_mem.ida_sensi = SUNTRUE;
    ida_mem.ida_sensMallocDone = SUNTRUE;

    /* create a Newton nonlinear solver object by default */
    let NLS = if ism == IDA_SIMULTANEOUS {
        SUNNonlinSol_NewtonSens(Ns + 1, &ida_mem.ida_delta, &ida_mem.ida_sunctx)
    } else {
        SUNNonlinSol_NewtonSens(Ns, &ida_mem.ida_delta, &ida_mem.ida_sunctx)
    };

    /* attach the nonlinear solver to the IDA memory */
    let retval = if ism == IDA_SIMULTANEOUS {
        IDASetNonlinearSolverSensSim(ida_mem, NLS)
    } else {
        IDASetNonlinearSolverSensStg(ida_mem, NLS)
    };

    /* check that the nonlinear solver was successfully attached */
    if retval != IDA_SUCCESS {
        IDAProcessError(Some(ida_mem), retval, line!(), "IDASensInit", file!(),
                        "Setting the nonlinear solver failed");
        return IDA_MEM_FAIL;
    }

    /* set ownership flag */
    if ism == IDA_SIMULTANEOUS {
        ida_mem.ownNLSsim = SUNTRUE;
    } else {
        ida_mem.ownNLSstg = SUNTRUE;
    }

    /* Sensitivity initialization was successful */
    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

/*
 * IDASensReInit
 *
 * IDASensReInit re-initializes IDAS's sensitivity related memory
 * for a problem, assuming it has already been allocated in prior
 * calls to IDAInit and IDASensInit.
 * All problem specification inputs are checked for errors.
 * The number of sensitivities Ns is assumed to be unchanged since
 * the previous call to IDASensInit.
 * If any error occurs during initialization, it is reported to the
 * file whose file pointer is errfp.
 * The return value is IDA_SUCCESS = 0 if no errors occurred, or
 * a negative value otherwise.
 */
pub fn IDASensReInit(ida_mem: &mut IDAMem, ism: i32, yS0: &[NVector], ypS0: &[NVector]) -> i32 {
    /* Was sensitivity initialized? */
    if !ida_mem.ida_sensMallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_SENS, line!(), "IDASensReInit", file!(), MSG_NO_SENSI);
        return IDA_NO_SENS;
    }

    /* Check if ism is legal */
    if ism != IDA_SIMULTANEOUS && ism != IDA_STAGGERED {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASensReInit", file!(), MSG_BAD_ISM);
        return IDA_ILL_INPUT;
    }
    ida_mem.ida_ism = ism;

    /* Check if yS0 and ypS0 are non-null */
    if yS0.is_empty() {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASensReInit", file!(), MSG_NULL_YYS0);
        return IDA_ILL_INPUT;
    }
    if ypS0.is_empty() {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASensReInit", file!(), MSG_NULL_YPS0);
        return IDA_ILL_INPUT;
    }

    /*-----------------------------------------------
      All error checking is complete at this point
      -----------------------------------------------*/

    /* Initialize the phiS array
       (C: cvals[is]=ONE + N_VScaleVectorArray; serial fused kernel
       reproduced inline; cannot-fail retval branches vanish.) */
    for is in 0..(ida_mem.ida_Ns as usize) {
        ida_mem.ida_phiS[0][is].data.copy_from_slice(&yS0[is].data);
        ida_mem.ida_phiS[1][is].data.copy_from_slice(&ypS0[is].data);
    }

    /* Initialize all sensitivity related counters */
    ida_mem.ida_nrSe = 0;
    ida_mem.ida_nreS = 0;
    ida_mem.ida_ncfnS = 0;
    ida_mem.ida_netfS = 0;
    ida_mem.ida_nniS = 0;
    ida_mem.ida_nnfS = 0;
    ida_mem.ida_nsetupsS = 0;

    /* Set default values for plist and pbar */
    for is in 0..(ida_mem.ida_Ns as usize) {
        ida_mem.ida_plist[is] = is as i32;
        ida_mem.ida_pbar[is] = ONE;
    }

    /* Sensitivities will be computed */
    ida_mem.ida_sensi = SUNTRUE;

    /* Check if the NLS exists, create the default NLS if needed */
    if (ism == IDA_SIMULTANEOUS && ida_mem.NLSsim.is_none())
        || (ism == IDA_STAGGERED && ida_mem.NLSstg.is_none())
    {
        /* create a Newton nonlinear solver object by default */
        let NLS = if ism == IDA_SIMULTANEOUS {
            SUNNonlinSol_NewtonSens(ida_mem.ida_Ns + 1, &ida_mem.ida_delta, &ida_mem.ida_sunctx)
        } else {
            SUNNonlinSol_NewtonSens(ida_mem.ida_Ns, &ida_mem.ida_delta, &ida_mem.ida_sunctx)
        };

        /* attach the nonlinear solver to the IDA memory */
        let retval = if ism == IDA_SIMULTANEOUS {
            IDASetNonlinearSolverSensSim(ida_mem, NLS)
        } else {
            IDASetNonlinearSolverSensStg(ida_mem, NLS)
        };

        /* check that the nonlinear solver was successfully attached */
        if retval != IDA_SUCCESS {
            IDAProcessError(Some(ida_mem), retval, line!(), "IDASensReInit", file!(),
                            "Setting the nonlinear solver failed");
            return IDA_MEM_FAIL;
        }

        /* set ownership flag */
        if ism == IDA_SIMULTANEOUS {
            ida_mem.ownNLSsim = SUNTRUE;
        } else {
            ida_mem.ownNLSstg = SUNTRUE;
        }

        /* initialize the NLS object, this assumes that the linear solver has
           already been initialized in IDAInit */
        let retval = if ism == IDA_SIMULTANEOUS {
            idaNlsInitSensSim(ida_mem)
        } else {
            idaNlsInitSensStg(ida_mem)
        };

        if retval != IDA_SUCCESS {
            IDAProcessError(Some(ida_mem), IDA_NLS_INIT_FAIL, line!(), "IDASensReInit", file!(), MSG_NLS_INIT_FAIL);
            return IDA_NLS_INIT_FAIL;
        }
    }

    /* Sensitivity re-initialization was successful */
    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

/*
 * IDASensSStolerances
 * IDASensSVtolerances
 * IDASensEEtolerances
 *
 * These functions specify the integration tolerances for sensitivity
 * variables. One of them MUST be called before the first call to IDASolve.
 *
 * IDASensSStolerances specifies scalar relative and absolute tolerances.
 * IDASensSVtolerances specifies scalar relative tolerance and a vector
 *   absolute tolerance for each sensitivity vector (a potentially different
 *   absolute tolerance for each vector component).
 * IDASensEEtolerances specifies that tolerances for sensitivity variables
 *   should be estimated from those provided for the state variables.
 */

pub fn IDASensSStolerances(ida_mem: &mut IDAMem, reltolS: f64, abstolS: &[f64]) -> i32 {
    /* Was sensitivity initialized? */

    if !ida_mem.ida_sensMallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_SENS, line!(), "IDASensSStolerances", file!(), MSG_NO_SENSI);
        return IDA_NO_SENS;
    }

    /* Test user-supplied tolerances */

    if reltolS < ZERO {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASensSStolerances", file!(), MSG_BAD_RTOLS);
        return IDA_ILL_INPUT;
    }

    if abstolS.is_empty() {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASensSStolerances", file!(), MSG_NULL_ATOLS);
        return IDA_ILL_INPUT;
    }

    for is in 0..(ida_mem.ida_Ns as usize) {
        if abstolS[is] < ZERO {
            IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASensSStolerances", file!(), MSG_BAD_ATOLS);
            return IDA_ILL_INPUT;
        }
    }

    /* Copy tolerances into memory */

    ida_mem.ida_itolS = IDA_SS;

    ida_mem.ida_rtolS = reltolS;

    if !ida_mem.ida_SatolSMallocDone {
        ida_mem.ida_SatolS = vec![0.0; ida_mem.ida_Ns as usize];
        ida_mem.ida_atolSmin0 = vec![SUNFALSE; ida_mem.ida_Ns as usize];
        ida_mem.ida_lrw += ida_mem.ida_Ns as i64;
        ida_mem.ida_SatolSMallocDone = SUNTRUE;
    }

    for is in 0..(ida_mem.ida_Ns as usize) {
        ida_mem.ida_SatolS[is] = abstolS[is];
        ida_mem.ida_atolSmin0[is] = abstolS[is] == ZERO;
    }

    IDA_SUCCESS
}

pub fn IDASensSVtolerances(ida_mem: &mut IDAMem, reltolS: f64, abstolS: &[NVector]) -> i32 {
    /* Was sensitivity initialized? */

    if !ida_mem.ida_sensMallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_SENS, line!(), "IDASensSVtolerances", file!(), MSG_NO_SENSI);
        return IDA_NO_SENS;
    }

    /* Test user-supplied tolerances */

    if reltolS < ZERO {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASensSVtolerances", file!(), MSG_BAD_RTOLS);
        return IDA_ILL_INPUT;
    }

    if abstolS.is_empty() {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASensSVtolerances", file!(), MSG_NULL_ATOLS);
        return IDA_ILL_INPUT;
    }

    let mut atolmin = vec![0.0; ida_mem.ida_Ns as usize];
    for is in 0..(ida_mem.ida_Ns as usize) {
        atolmin[is] = N_VMin(&abstolS[is]);
        if atolmin[is] < ZERO {
            IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASensSVtolerances", file!(), MSG_BAD_ATOLS);
            return IDA_ILL_INPUT;
        }
    }

    ida_mem.ida_itolS = IDA_SV;
    ida_mem.ida_rtolS = reltolS;

    if SUNFALSE == ida_mem.ida_VatolSMallocDone {
        ida_mem.ida_VatolS = N_VCloneVectorArray(ida_mem.ida_Ns as usize, &ida_mem.ida_tempv1);
        ida_mem.ida_atolSmin0 = vec![SUNFALSE; ida_mem.ida_Ns as usize];
        ida_mem.ida_lrw += ida_mem.ida_Ns as i64 * ida_mem.ida_lrw1;
        ida_mem.ida_liw += ida_mem.ida_Ns as i64 * ida_mem.ida_liw1;
        ida_mem.ida_VatolSMallocDone = SUNTRUE;
    }

    /* (C: cvals[is]=ONE + N_VScaleVectorArray(Ns, cvals, abstolS,
       VatolS); serial fused kernel reproduced inline; cannot-fail
       retval branch vanishes.) */
    for is in 0..(ida_mem.ida_Ns as usize) {
        ida_mem.ida_atolSmin0[is] = atolmin[is] == ZERO;
        ida_mem.ida_VatolS[is].data.copy_from_slice(&abstolS[is].data);
    }

    IDA_SUCCESS
}

pub fn IDASensEEtolerances(ida_mem: &mut IDAMem) -> i32 {
    /* Was sensitivity initialized? */

    if !ida_mem.ida_sensMallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_SENS, line!(), "IDASensEEtolerances", file!(), MSG_NO_SENSI);
        return IDA_NO_SENS;
    }

    ida_mem.ida_itolS = IDA_EE;

    IDA_SUCCESS
}

/*
 * (C rhsQS may be NULL to select the internal DQ rhs, so rhsQS is
 * Option<IDAQuadSensRhsFn> here; the C branches setting
 * ida_user_dataQS to IDA_mem (DQ) or ida_user_data (user fn) are
 * not stored per the pinned convention — the DQ path in this file
 * operates on &mut IDAMem and a user rhsQS receives ida_user_data.)
 */
pub fn IDAQuadSensInit(ida_mem: &mut IDAMem, rhsQS: Option<IDAQuadSensRhsFn>, yQS0: &[NVector]) -> i32 {
    /* Check if sensitivity analysis is active */
    if !ida_mem.ida_sensi {
        IDAProcessError(Some(ida_mem), IDA_NO_SENS, line!(), "IDAQuadSensInit", file!(), MSG_NO_SENSI);
        return IDA_NO_SENS;
    }

    /* Verify yQS0 parameter. */
    if yQS0.is_empty() {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAQuadSensInit", file!(), MSG_NULL_YQS0);
        return IDA_ILL_INPUT;
    }

    /* Allocate vector needed for quadratures' sensitivities;
       infallible in Rust. */
    IDAQuadSensAllocVectors(ida_mem, &yQS0[0]);

    /* Error checking complete. */
    if let Some(rhsQS) = rhsQS {
        ida_mem.ida_rhsQSDQ = SUNFALSE;
        ida_mem.ida_rhsQS = Some(rhsQS);
    } else {
        ida_mem.ida_rhsQSDQ = SUNTRUE;
        ida_mem.ida_rhsQS = None;
    }

    /* Initialize phiQS[0] in the history array
       (C: cvals[is]=ONE + N_VScaleVectorArray; serial fused kernel
       reproduced inline; cannot-fail retval branch vanishes.) */
    for is in 0..(ida_mem.ida_Ns as usize) {
        ida_mem.ida_phiQS[0][is].data.copy_from_slice(&yQS0[is].data);
    }

    /* Initialize all sensitivities related counters. */
    ida_mem.ida_nrQSe = 0;
    ida_mem.ida_nrQeS = 0;
    ida_mem.ida_netfQS = 0;

    /* Everything all right, set the flags and return with success. */
    ida_mem.ida_quadr_sensi = SUNTRUE;
    ida_mem.ida_quadSensMallocDone = SUNTRUE;

    IDA_SUCCESS
}

pub fn IDAQuadSensReInit(ida_mem: &mut IDAMem, yQS0: &[NVector]) -> i32 {
    /* Check if sensitivity analysis is active */
    if !ida_mem.ida_sensi {
        IDAProcessError(Some(ida_mem), IDA_NO_SENS, line!(), "IDAQuadSensReInit", file!(), MSG_NO_SENSI);
        return IDA_NO_SENS;
    }

    /* Was sensitivity for quadrature already initialized? */
    if !ida_mem.ida_quadSensMallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_QUADSENS, line!(), "IDAQuadSensReInit", file!(), MSG_NO_QUADSENSI);
        return IDA_NO_QUADSENS;
    }

    /* Verify yQS0 parameter. */
    if yQS0.is_empty() {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAQuadSensReInit", file!(), MSG_NULL_YQS0);
        return IDA_ILL_INPUT;
    }

    /* Error checking complete at this point. */

    /* Initialize phiQS[0] in the history array
       (C: cvals[is]=ONE + N_VScaleVectorArray; serial fused kernel
       reproduced inline; cannot-fail retval branch vanishes.) */
    for is in 0..(ida_mem.ida_Ns as usize) {
        ida_mem.ida_phiQS[0][is].data.copy_from_slice(&yQS0[is].data);
    }

    /* Initialize all sensitivities related counters. */
    ida_mem.ida_nrQSe = 0;
    ida_mem.ida_nrQeS = 0;
    ida_mem.ida_netfQS = 0;

    /* Everything all right, set the flags and return with success. */
    ida_mem.ida_quadr_sensi = SUNTRUE;

    IDA_SUCCESS
}

/*
 * IDAQuadSensSStolerances
 * IDAQuadSensSVtolerances
 * IDAQuadSensEEtolerances
 *
 * These functions specify the integration tolerances for quadrature
 * sensitivity variables. One of them MUST be called before the first
 * call to IDAS IF these variables are included in the error test.
 *
 * IDAQuadSensSStolerances specifies scalar relative and absolute tolerances.
 * IDAQuadSensSVtolerances specifies scalar relative tolerance and a vector
 *   absolute tolerance for each quadrature sensitivity vector (a potentially
 *   different absolute tolerance for each vector component).
 * IDAQuadSensEEtolerances specifies that tolerances for sensitivity variables
 *   should be estimated from those provided for the quadrature variables.
 *   In this case, tolerances for the quadrature variables must be
 *   specified through a call to one of IDAQuad**tolerances.
 */

pub fn IDAQuadSensSStolerances(ida_mem: &mut IDAMem, reltolQS: f64, abstolQS: &[f64]) -> i32 {
    /* Check if sensitivity analysis is active */
    if !ida_mem.ida_sensi {
        IDAProcessError(Some(ida_mem), IDA_NO_SENS, line!(), "IDAQuadSensSStolerances", file!(), MSG_NO_SENSI);
        return IDA_NO_SENS;
    }

    /* Was sensitivity for quadrature already initialized? */
    if !ida_mem.ida_quadSensMallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_QUADSENS, line!(), "IDAQuadSensSStolerances", file!(), MSG_NO_QUADSENSI);
        return IDA_NO_QUADSENS;
    }

    /* Test user-supplied tolerances */

    if reltolQS < ZERO {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAQuadSensSStolerances", file!(), MSG_BAD_RELTOLQS);
        return IDA_ILL_INPUT;
    }

    if abstolQS.is_empty() {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAQuadSensSStolerances", file!(), MSG_NULL_ABSTOLQS);
        return IDA_ILL_INPUT;
    }

    for is in 0..(ida_mem.ida_Ns as usize) {
        if abstolQS[is] < ZERO {
            IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAQuadSensSStolerances", file!(), MSG_BAD_ABSTOLQS);
            return IDA_ILL_INPUT;
        }
    }

    /* Save data. */
    ida_mem.ida_itolQS = IDA_SS;
    ida_mem.ida_rtolQS = reltolQS;

    if !ida_mem.ida_SatolQSMallocDone {
        ida_mem.ida_SatolQS = vec![0.0; ida_mem.ida_Ns as usize];
        ida_mem.ida_atolQSmin0 = vec![SUNFALSE; ida_mem.ida_Ns as usize];
        ida_mem.ida_lrw += ida_mem.ida_Ns as i64;
        ida_mem.ida_SatolQSMallocDone = SUNTRUE;
    }

    for is in 0..(ida_mem.ida_Ns as usize) {
        ida_mem.ida_SatolQS[is] = abstolQS[is];
        ida_mem.ida_atolQSmin0[is] = abstolQS[is] == ZERO;
    }

    IDA_SUCCESS
}

pub fn IDAQuadSensSVtolerances(ida_mem: &mut IDAMem, reltolQS: f64, abstolQS: &[NVector]) -> i32 {
    /* Check if sensitivity analysis is active */
    if !ida_mem.ida_sensi {
        IDAProcessError(Some(ida_mem), IDA_NO_SENS, line!(), "IDAQuadSensSVtolerances", file!(), MSG_NO_SENSI);
        return IDA_NO_SENS;
    }

    /* Was sensitivity for quadrature already initialized? */
    if !ida_mem.ida_quadSensMallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_QUADSENS, line!(), "IDAQuadSensSVtolerances", file!(), MSG_NO_QUADSENSI);
        return IDA_NO_QUADSENS;
    }

    /* Test user-supplied tolerances */

    if reltolQS < ZERO {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAQuadSensSVtolerances", file!(), MSG_BAD_RELTOLQS);
        return IDA_ILL_INPUT;
    }

    if abstolQS.is_empty() {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAQuadSensSVtolerances", file!(), MSG_NULL_ABSTOLQS);
        return IDA_ILL_INPUT;
    }

    let mut atolmin = vec![0.0; ida_mem.ida_Ns as usize];
    for is in 0..(ida_mem.ida_Ns as usize) {
        atolmin[is] = N_VMin(&abstolQS[is]);
        if atolmin[is] < ZERO {
            IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAQuadSensSVtolerances", file!(), MSG_BAD_ABSTOLQS);
            return IDA_ILL_INPUT;
        }
    }

    /* Save data. */
    ida_mem.ida_itolQS = IDA_SV;
    ida_mem.ida_rtolQS = reltolQS;

    if !ida_mem.ida_VatolQSMallocDone {
        /* (C clones from abstolQS[0] here — not tempv1 as in
           IDASensSVtolerances — and counts quadrature vector sizes.) */
        ida_mem.ida_VatolQS = N_VCloneVectorArray(ida_mem.ida_Ns as usize, &abstolQS[0]);
        ida_mem.ida_atolQSmin0 = vec![SUNFALSE; ida_mem.ida_Ns as usize];
        ida_mem.ida_lrw += ida_mem.ida_Ns as i64 * ida_mem.ida_lrw1Q;
        ida_mem.ida_liw += ida_mem.ida_Ns as i64 * ida_mem.ida_liw1Q;
        ida_mem.ida_VatolQSMallocDone = SUNTRUE;
    }

    /* (C: cvals[is]=ONE + N_VScaleVectorArray(Ns, cvals, abstolQS,
       VatolQS); serial fused kernel reproduced inline; cannot-fail
       retval branch vanishes.) */
    for is in 0..(ida_mem.ida_Ns as usize) {
        ida_mem.ida_atolQSmin0[is] = atolmin[is] == ZERO;
        ida_mem.ida_VatolQS[is].data.copy_from_slice(&abstolQS[is].data);
    }

    IDA_SUCCESS
}

pub fn IDAQuadSensEEtolerances(ida_mem: &mut IDAMem) -> i32 {
    /* Check if sensitivity analysis is active */
    if !ida_mem.ida_sensi {
        IDAProcessError(Some(ida_mem), IDA_NO_SENS, line!(), "IDAQuadSensEEtolerances", file!(), MSG_NO_SENSI);
        return IDA_NO_SENS;
    }

    /* Was sensitivity for quadrature already initialized? */
    if !ida_mem.ida_quadSensMallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_QUADSENS, line!(), "IDAQuadSensEEtolerances", file!(), MSG_NO_QUADSENSI);
        return IDA_NO_QUADSENS;
    }

    ida_mem.ida_itolQS = IDA_EE;

    IDA_SUCCESS
}

/*
 * IDASensToggleOff
 *
 * IDASensToggleOff deactivates sensitivity calculations.
 * It does NOT deallocate sensitivity-related memory.
 */
pub fn IDASensToggleOff(ida_mem: &mut IDAMem) -> i32 {
    /* Disable sensitivities */
    ida_mem.ida_sensi = SUNFALSE;
    ida_mem.ida_quadr_sensi = SUNFALSE;

    IDA_SUCCESS
}

/*
 * IDARootInit
 *
 * IDARootInit initializes a rootfinding problem to be solved
 * during the integration of the DAE system.  It loads the root
 * function pointer and the number of root functions, and allocates
 * workspace memory.  The return value is IDA_SUCCESS = 0 if no
 * errors occurred, or a negative value otherwise.
 */
pub fn IDARootInit(ida_mem: &mut IDAMem, nrtfn: i32, g: Option<IDARootFn>) -> i32 {
    let nrt = if nrtfn < 0 { 0 } else { nrtfn };

    /* If rerunning IDARootInit() with a different number of root
       functions (changing number of gfun components), then free
       currently held memory resources */
    if nrt != ida_mem.ida_nrtfn && ida_mem.ida_nrtfn > 0 {
        ida_mem.ida_glo = Vec::new();
        ida_mem.ida_ghi = Vec::new();
        ida_mem.ida_grout = Vec::new();
        ida_mem.ida_iroots = Vec::new();
        ida_mem.ida_rootdir = Vec::new();
        ida_mem.ida_gactive = Vec::new();

        ida_mem.ida_lrw -= 3 * ida_mem.ida_nrtfn as i64;
        ida_mem.ida_liw -= 3 * ida_mem.ida_nrtfn as i64;
    }

    /* If IDARootInit() was called with nrtfn == 0, then set ida_nrtfn to
       zero and ida_gfun to NULL before returning */
    if nrt == 0 {
        ida_mem.ida_nrtfn = nrt;
        ida_mem.ida_gfun = None;
        return IDA_SUCCESS;
    }

    /* If rerunning IDARootInit() with the same number of root functions
       (not changing number of gfun components), then check if the root
       function argument has changed */
    /* If g != NULL then return as currently reserved memory resources
       will suffice */
    if nrt == ida_mem.ida_nrtfn {
        match g {
            None => {
                ida_mem.ida_glo = Vec::new();
                ida_mem.ida_ghi = Vec::new();
                ida_mem.ida_grout = Vec::new();
                ida_mem.ida_iroots = Vec::new();
                ida_mem.ida_rootdir = Vec::new();
                ida_mem.ida_gactive = Vec::new();

                ida_mem.ida_lrw -= 3 * nrt as i64;
                ida_mem.ida_liw -= 3 * nrt as i64;

                IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDARootInit", file!(),
                                MSG_ROOT_FUNC_NULL);
                return IDA_ILL_INPUT;
            }
            Some(gf) => {
                /* (covers both the g == ida_gfun and the changed
                   non-NULL g cases of the C code — same observable
                   behavior) */
                ida_mem.ida_gfun = Some(gf);
                return IDA_SUCCESS;
            }
        }
    }

    /* Set variable values in IDA memory block */
    ida_mem.ida_nrtfn = nrt;
    match g {
        None => {
            IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDARootInit", file!(),
                            MSG_ROOT_FUNC_NULL);
            return IDA_ILL_INPUT;
        }
        Some(gf) => ida_mem.ida_gfun = Some(gf),
    }

    /* Allocate necessary memory and return (Vec allocation is
       infallible, so the cascading free/MSG_MEM_FAIL branches of the
       C code vanish) */
    let n = nrt as usize;
    ida_mem.ida_glo = vec![ZERO; n];
    ida_mem.ida_ghi = vec![ZERO; n];
    ida_mem.ida_grout = vec![ZERO; n];
    ida_mem.ida_iroots = vec![0; n];

    /* Set default values for rootdir (both directions) */
    ida_mem.ida_rootdir = vec![0; n];

    /* Set default values for gactive (all active) */
    ida_mem.ida_gactive = vec![SUNTRUE; n];

    ida_mem.ida_lrw += 3 * nrt as i64;
    ida_mem.ida_liw += 3 * nrt as i64;

    IDA_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Linear solver module dispatch
 *
 * In C these are the ida_linit/ida_lperf function pointers (plus
 * ida_lsetup/ida_lsolve, whose only call sites are in idas_nls*.c and
 * idas_ic.c and whose dispatch therefore lives there); the module is
 * taken out of IDAMem for the duration of a call so its routine can
 * borrow the integrator memory mutably (donor cv_*_dispatch
 * pattern).
 * -----------------------------------------------------------------
 */

/* mirrors C's `ida_mem->ida_linit != NULL` guard + call */
pub(crate) fn ida_linit_dispatch(ida_mem: &mut IDAMem) -> i32 {
    let mut lmem = std::mem::take(&mut ida_mem.ida_lmem);
    let ier = match &mut lmem {
        LsModule::None => 0,
        LsModule::Ls(ls) => idaLsInitialize(ida_mem, ls),
    };
    ida_mem.ida_lmem = lmem;
    ier
}

/* mirrors C's `if (ida_mem->ida_lperf != NULL) ida_lperf(IDA_mem, perftask)`:
   idaLsSetLinearSolver installs ida_lperf only for iterative
   SUNLinearSolver objects, so the guard dispatches on
   IDALsMem.iterative (see idas_impl.rs) */
pub(crate) fn ida_lperf_dispatch(ida_mem: &mut IDAMem, perftask: i32) {
    let mut lmem = std::mem::take(&mut ida_mem.ida_lmem);
    if let LsModule::Ls(ls) = &mut lmem {
        if ls.iterative {
            idaLsPerf(ida_mem, ls, perftask);
        }
    }
    ida_mem.ida_lmem = lmem;
}

/*
 * -----------------------------------------------------------------
 * Error weight function dispatch
 *
 * In C the ewt call sites are `ida_efun(ida_phi[0], ida_ewt,
 * ida_edata)` with edata = user_data (user efun) or IDA_mem
 * (internal IDAEwtSet).  ycur is always phi[0] at every call site.
 * The weight vector is detached from IDAMem by the caller so ycur
 * may be borrowed from the same struct (donor cv_efun_* pattern).
 * -----------------------------------------------------------------
 */

pub(crate) fn ida_efun_dispatch(ida_mem: &mut IDAMem, w: &mut NVector) -> i32 {
    if ida_mem.ida_user_efun {
        let efun = ida_mem.ida_efun.unwrap();
        /* edata = user_data for a user efun */
        let IDAMem { ida_phi, ida_user_data, .. } = ida_mem;
        efun(&ida_phi[0], w, ida_user_data)
    } else {
        /* edata = IDA_mem for the internal efun */
        IDAEwtSet(ida_mem, &ida_mem.ida_phi[0], w)
    }
}

/* Apply the efun to the integrator's own ewt vector. */
pub(crate) fn ida_efun_apply_to_ewt(ida_mem: &mut IDAMem) -> i32 {
    let mut w = std::mem::take(&mut ida_mem.ida_ewt);
    let flag = ida_efun_dispatch(ida_mem, &mut w);
    ida_mem.ida_ewt = w;
    flag
}

/*
 * -----------------------------------------------------------------
 * Main solver function
 * -----------------------------------------------------------------
 */

/*
 * IDASolve
 *
 * This routine is the main driver of the IDA package.
 *
 * It integrates over an independent variable interval defined by the user,
 * by calling IDAStep to take internal independent variable steps.
 *
 * The first time that IDASolve is called for a successfully initialized
 * problem, it computes a tentative initial step size.
 *
 * IDASolve supports two modes, specified by itask:
 * In the IDA_NORMAL mode, the solver steps until it passes tout and then
 * interpolates to obtain y(tout) and yp(tout).
 * In the IDA_ONE_STEP mode, it takes one internal step and returns.
 *
 * IDASolve returns integer values corresponding to success and failure as below:
 *
 * successful returns:
 *
 * IDA_SUCCESS
 * IDA_TSTOP_RETURN
 *
 * failed returns:
 *
 * IDA_ILL_INPUT
 * IDA_TOO_MUCH_WORK
 * IDA_MEM_NULL
 * IDA_TOO_MUCH_ACC
 * IDA_CONV_FAIL
 * IDA_LSETUP_FAIL
 * IDA_LSOLVE_FAIL
 * IDA_CONSTR_FAIL
 * IDA_ERR_FAIL
 * IDA_REP_RES_ERR
 * IDA_RES_FAIL
 *
 * Aliasing note (workspace rule 5): in C, `IDA_mem->ida_yy = yret`
 * and `ida_yp = ypret` alias the user's output vectors for the whole
 * call; the internal routines (IDARcheck*, IDAStep/IDANls) then read
 * and write the user buffers through ida_yy/ida_yp.  Here ida_yy /
 * ida_yp stay owned by IDAMem: their contents are synchronized from
 * yret/ypret at the alias-establishment point, and copied back to
 * yret/ypret at every return path where the C code relied on the
 * alias to deliver the solution.  Paths where C explicitly calls
 * IDAGetSolution(..., yret, ypret) pass the user vectors directly,
 * exactly as C does.
 */
pub fn IDASolve(
    ida_mem: &mut IDAMem,
    tout: f64,
    tret: &mut f64,
    yret: &mut NVector,
    ypret: &mut NVector,
    itask: i32,
) -> i32 {
    /* Check for legal inputs in all cases. */

    /* Check if problem was malloc'ed */

    if !ida_mem.ida_MallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_MALLOC, line!(), "IDASolve", file!(), MSG_NO_MALLOC);
        return IDA_NO_MALLOC;
    }

    /* Check for legal arguments */

    if yret.is_empty() {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASolve", file!(), MSG_YRET_NULL);
        return IDA_ILL_INPUT;
    }
    /* C: IDA_mem->ida_yy = yret (alias); here: synchronize contents */
    ida_mem.ida_yy.data.copy_from_slice(&yret.data);

    if ypret.is_empty() {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASolve", file!(), MSG_YPRET_NULL);
        return IDA_ILL_INPUT;
    }
    /* C: IDA_mem->ida_yp = ypret (alias); here: synchronize contents */
    ida_mem.ida_yp.data.copy_from_slice(&ypret.data);

    /* (tret == NULL cannot occur: tret is &mut f64.) */

    if itask != IDA_NORMAL && itask != IDA_ONE_STEP {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASolve", file!(), MSG_BAD_ITASK);
        return IDA_ILL_INPUT;
    }

    /* Sensitivity-specific tests (if using internal DQ functions) */
    if ida_mem.ida_sensi && ida_mem.ida_resSDQ {
        /* (C: ida_user_dataS = ida_mem — the DQ self-pointer is not
           stored, per the pinned convention; the DQ residual operates
           on &mut IDAMem directly.) */
        /* Test if we have the problem parameters */
        if ida_mem.ida_p.is_empty() {
            IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASolve", file!(), MSG_NULL_P);
            return IDA_ILL_INPUT;
        }
    }

    if ida_mem.ida_quadr_sensi && ida_mem.ida_rhsQSDQ {
        /* (C: ida_user_dataQS = ida_mem — not stored, as above.) */
        /* Test if we have the problem parameters */
        if ida_mem.ida_p.is_empty() {
            IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASolve", file!(), MSG_NULL_P);
            return IDA_ILL_INPUT;
        }
    }

    if ida_mem.ida_nst == 0 {
        /* This is the first call */

        /* Check inputs to IDA for correctness and consistency */

        if !ida_mem.ida_SetupDone {
            let ier = IDAInitialSetup(ida_mem);
            if ier != IDA_SUCCESS {
                return ier;
            }
            ida_mem.ida_SetupDone = SUNTRUE;
        }

        /* On first call, check for tout - tn too small, set initial hh,
           check for approach to tstop, and scale phi[1], phiQ[1], and phiS[1] by hh.
           Also check for zeros of root function g at and near t0.    */

        let tdist = SUNRabs(tout - ida_mem.ida_tn);
        let troundoff = TWO * ida_mem.ida_uround * (SUNRabs(ida_mem.ida_tn) + SUNRabs(tout));
        if tdist == ZERO || tdist < troundoff {
            IDAProcessError(Some(ida_mem), IDA_TOO_CLOSE, line!(), "IDASolve", file!(), MSG_TOO_CLOSE);
            return IDA_TOO_CLOSE;
        }

        /* Set initial h */

        ida_mem.ida_hh = ida_mem.ida_hin;
        if ida_mem.ida_hh != ZERO && (tout - ida_mem.ida_tn) * ida_mem.ida_hh < ZERO {
            IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASolve", file!(), MSG_BAD_HINIT);
            return IDA_ILL_INPUT;
        }

        if ida_mem.ida_hh == ZERO {
            ida_mem.ida_hh = PT001 * tdist;
            let mut ypnorm = IDAWrmsNorm(ida_mem, &ida_mem.ida_phi[1], &ida_mem.ida_ewt,
                                         ida_mem.ida_suppressalg);
            if ida_mem.ida_errconQ {
                ypnorm = IDAQuadWrmsNormUpdate(ida_mem, ypnorm, &ida_mem.ida_phiQ[1],
                                               &ida_mem.ida_ewtQ);
            }
            if ida_mem.ida_errconS {
                ypnorm = IDASensWrmsNormUpdate(ida_mem, ypnorm, &ida_mem.ida_phiS[1],
                                               &ida_mem.ida_ewtS,
                                               ida_mem.ida_suppressalg);
            }
            if ida_mem.ida_errconQS {
                ypnorm = IDAQuadSensWrmsNormUpdate(ida_mem, ypnorm, &ida_mem.ida_phiQS[1],
                                                   &ida_mem.ida_ewtQS);
            }

            if ypnorm > HALF / ida_mem.ida_hh {
                ida_mem.ida_hh = HALF / ypnorm;
            }
            if tout < ida_mem.ida_tn {
                ida_mem.ida_hh = -ida_mem.ida_hh;
            }
        }

        /* Enforce hmax and hmin */

        let rh = SUNRabs(ida_mem.ida_hh) * ida_mem.ida_hmax_inv;
        if rh > ONE {
            ida_mem.ida_hh /= rh;
        }
        if SUNRabs(ida_mem.ida_hh) < ida_mem.ida_hmin {
            ida_mem.ida_hh *= ida_mem.ida_hmin / SUNRabs(ida_mem.ida_hh);
        }

        /* Check for approach to tstop */

        if ida_mem.ida_tstopset {
            if (ida_mem.ida_tstop - ida_mem.ida_tn) * ida_mem.ida_hh <= ZERO {
                IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASolve", file!(),
                                &ida_msg_g(MSG_BAD_TSTOP, &[ida_mem.ida_tstop, ida_mem.ida_tn]));
                return IDA_ILL_INPUT;
            }
            if (ida_mem.ida_tn + ida_mem.ida_hh - ida_mem.ida_tstop) * ida_mem.ida_hh > ZERO {
                ida_mem.ida_hh =
                    (ida_mem.ida_tstop - ida_mem.ida_tn) * (ONE - FOUR * ida_mem.ida_uround);
            }
        }

        ida_mem.ida_h0u = ida_mem.ida_hh;
        ida_mem.ida_kk = 0;
        ida_mem.ida_kused = 0; /* set in case of an error return before a step */

        /* Check for exact zeros of the root functions at or near t0. */
        if ida_mem.ida_nrtfn > 0 {
            let ier = IDARcheck1(ida_mem);
            if ier == IDA_RTFUNC_FAIL {
                IDAProcessError(Some(ida_mem), IDA_RTFUNC_FAIL, line!(), "IDASolve", file!(),
                                &ida_msg_g(MSG_RTFUNC_FAILED, &[ida_mem.ida_tn]));
                return IDA_RTFUNC_FAIL;
            }
        }

        /* set phi[1] = hh*y' */
        let hh = ida_mem.ida_hh;
        ida_mem.ida_phi[1].scale_inplace(hh);

        /* set phiQ[1] = hh*yQ' */
        if ida_mem.ida_quadr {
            ida_mem.ida_phiQ[1].scale_inplace(hh);
        }

        /* (C: cvals[is]=hh + in-place N_VScaleVectorArray for phiS[1]
           and phiQS[1]; serial fused kernels reproduced inline;
           cannot-fail retval branches vanish.) */

        if ida_mem.ida_sensi {
            /* set phiS[1][i] = hh*yS_i' */
            for is in 0..(ida_mem.ida_Ns as usize) {
                ida_mem.ida_phiS[1][is].scale_inplace(hh);
            }
        }

        if ida_mem.ida_quadr_sensi {
            for is in 0..(ida_mem.ida_Ns as usize) {
                ida_mem.ida_phiQS[1][is].scale_inplace(hh);
            }
        }

        /* Set the convergence test constants epsNewt and toldel */
        ida_mem.ida_epsNewt = ida_mem.ida_epcon;
        ida_mem.ida_toldel = PT0001 * ida_mem.ida_epsNewt;
    } /* end of first-call block. */

    /* Call lperf function and set nstloc for later performance testing. */

    ida_lperf_dispatch(ida_mem, 0);
    let mut nstloc: i64 = 0;

    /* If not the first call, perform all stopping tests. */

    if ida_mem.ida_nst > 0 {
        /* First, check for a root in the last step taken, other than the
           last root found, if any.  If itask = IDA_ONE_STEP and y(tn) was not
           returned because of an intervening root, return y(tn) now.     */

        if ida_mem.ida_nrtfn > 0 {
            let irfndp = ida_mem.ida_irfnd;

            let ier = IDARcheck2(ida_mem);

            if ier == CLOSERT {
                IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASolve", file!(),
                                &ida_msg_g(MSG_CLOSE_ROOTS, &[ida_mem.ida_tlo]));
                /* (C's yret holds IDARcheck2's interpolation through the
                   ida_yy alias; replicate) */
                yret.data.copy_from_slice(&ida_mem.ida_yy.data);
                ypret.data.copy_from_slice(&ida_mem.ida_yp.data);
                return IDA_ILL_INPUT;
            } else if ier == IDA_RTFUNC_FAIL {
                IDAProcessError(Some(ida_mem), IDA_RTFUNC_FAIL, line!(), "IDASolve", file!(),
                                &ida_msg_g(MSG_RTFUNC_FAILED, &[ida_mem.ida_tlo]));
                yret.data.copy_from_slice(&ida_mem.ida_yy.data);
                ypret.data.copy_from_slice(&ida_mem.ida_yp.data);
                return IDA_RTFUNC_FAIL;
            } else if ier == RTFOUND {
                ida_mem.ida_tretlast = ida_mem.ida_tlo;
                *tret = ida_mem.ida_tlo;
                /* ida_yy/ida_yp (root state) alias yret/ypret in C */
                yret.data.copy_from_slice(&ida_mem.ida_yy.data);
                ypret.data.copy_from_slice(&ida_mem.ida_yp.data);
                return IDA_ROOT_RETURN;
            }

            /* If tn is distinct from tretlast (within roundoff),
               check remaining interval for roots */
            let troundoff = HUNDRED * ida_mem.ida_uround *
                (SUNRabs(ida_mem.ida_tn) + SUNRabs(ida_mem.ida_hh));
            if SUNRabs(ida_mem.ida_tn - ida_mem.ida_tretlast) > troundoff {
                let ier = IDARcheck3(ida_mem, tout, itask);
                if ier == IDA_SUCCESS {
                    /* no root found */
                    ida_mem.ida_irfnd = 0;
                    if irfndp == 1 && itask == IDA_ONE_STEP {
                        ida_mem.ida_tretlast = ida_mem.ida_tn;
                        *tret = ida_mem.ida_tn;
                        let _ = IDAGetSolution(ida_mem, ida_mem.ida_tn, yret, ypret);
                        return IDA_SUCCESS;
                    }
                } else if ier == RTFOUND {
                    /* a new root was found */
                    ida_mem.ida_irfnd = 1;
                    ida_mem.ida_tretlast = ida_mem.ida_tlo;
                    *tret = ida_mem.ida_tlo;
                    /* ida_yy/ida_yp (interpolated root solution) alias
                       yret/ypret in C */
                    yret.data.copy_from_slice(&ida_mem.ida_yy.data);
                    ypret.data.copy_from_slice(&ida_mem.ida_yp.data);
                    return IDA_ROOT_RETURN;
                } else if ier == IDA_RTFUNC_FAIL {
                    /* g failed */
                    IDAProcessError(Some(ida_mem), IDA_RTFUNC_FAIL, line!(), "IDASolve", file!(),
                                    &ida_msg_g(MSG_RTFUNC_FAILED, &[ida_mem.ida_tlo]));
                    yret.data.copy_from_slice(&ida_mem.ida_yy.data);
                    ypret.data.copy_from_slice(&ida_mem.ida_yp.data);
                    return IDA_RTFUNC_FAIL;
                }
            }
        } /* end of root stop check */

        /* Now test for all other stop conditions. */

        let istate = IDAStopTest1(ida_mem, tout, tret, yret, ypret, itask);
        if istate != CONTINUE_STEPS {
            return istate;
        }
    }

    /* Looping point for internal steps. */

    let mut istate;
    loop {
        /* Check for too many steps taken. */

        if ida_mem.ida_mxstep > 0 && nstloc >= ida_mem.ida_mxstep {
            IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASolve", file!(),
                            &ida_msg_g(MSG_MAX_STEPS, &[ida_mem.ida_tn]));
            istate = IDA_TOO_MUCH_WORK;
            ida_mem.ida_tretlast = ida_mem.ida_tn;
            *tret = ida_mem.ida_tn;
            /* Here yy=yret and yp=ypret already have the current solution
               (through the alias in C; copy the owned vectors back). */
            yret.data.copy_from_slice(&ida_mem.ida_yy.data);
            ypret.data.copy_from_slice(&ida_mem.ida_yp.data);
            break;
        }

        /* Call lperf to generate warnings of poor performance. */

        ida_lperf_dispatch(ida_mem, 1);

        /* Reset and check ewt, ewtQ, ewtS and ewtQS (if not first call). */

        if ida_mem.ida_nst > 0 {
            let ier = ida_efun_apply_to_ewt(ida_mem);

            if ier != 0 {
                if ida_mem.ida_itol == IDA_WF {
                    IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASolve", file!(),
                                    &ida_msg_g(MSG_EWT_NOW_FAIL, &[ida_mem.ida_tn]));
                } else {
                    IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASolve", file!(),
                                    &ida_msg_g(MSG_EWT_NOW_BAD, &[ida_mem.ida_tn]));
                }

                istate = IDA_ILL_INPUT;
                let _ = IDAGetSolution(ida_mem, ida_mem.ida_tn, yret, ypret);
                ida_mem.ida_tretlast = ida_mem.ida_tn;
                *tret = ida_mem.ida_tn;
                break;
            }

            if ida_mem.ida_quadr && ida_mem.ida_errconQ {
                /* (weight vector detached so qcur may be borrowed from
                   the same struct — donor efun pattern) */
                let mut wQ = std::mem::take(&mut ida_mem.ida_ewtQ);
                let ier = IDAQuadEwtSet(ida_mem, &ida_mem.ida_phiQ[0], &mut wQ);
                ida_mem.ida_ewtQ = wQ;
                if ier != 0 {
                    IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASolve", file!(),
                                    &ida_msg_g(MSG_EWTQ_NOW_BAD, &[ida_mem.ida_tn]));
                    istate = IDA_ILL_INPUT;
                    let _ = IDAGetSolution(ida_mem, ida_mem.ida_tn, yret, ypret);
                    ida_mem.ida_tretlast = ida_mem.ida_tn;
                    *tret = ida_mem.ida_tn;
                    break;
                }
            }

            if ida_mem.ida_sensi {
                let mut wS = std::mem::take(&mut ida_mem.ida_ewtS);
                /* (phiS[0] also detached: IDASensEwtSet takes &mut
                   IDAMem so its EE branch can call a user efun) */
                let yS0 = std::mem::take(&mut ida_mem.ida_phiS[0]);
                let ier = IDASensEwtSet(ida_mem, &yS0, &mut wS);
                ida_mem.ida_phiS[0] = yS0;
                ida_mem.ida_ewtS = wS;
                if ier != 0 {
                    IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASolve", file!(),
                                    &ida_msg_g(MSG_EWTS_NOW_BAD, &[ida_mem.ida_tn]));
                    istate = IDA_ILL_INPUT;
                    let _ = IDAGetSolution(ida_mem, ida_mem.ida_tn, yret, ypret);
                    ida_mem.ida_tretlast = ida_mem.ida_tn;
                    *tret = ida_mem.ida_tn;
                    break;
                }
            }

            if ida_mem.ida_quadr_sensi && ida_mem.ida_errconQS {
                let mut wQS = std::mem::take(&mut ida_mem.ida_ewtQS);
                let ier = IDAQuadSensEwtSet(ida_mem, &ida_mem.ida_phiQS[0], &mut wQS);
                ida_mem.ida_ewtQS = wQS;
                if ier != 0 {
                    IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASolve", file!(),
                                    &ida_msg_g(MSG_EWTQS_NOW_BAD, &[ida_mem.ida_tn]));
                    istate = IDA_ILL_INPUT;
                    let _ = IDAGetSolution(ida_mem, ida_mem.ida_tn, yret, ypret);
                    ida_mem.ida_tretlast = ida_mem.ida_tn;
                    *tret = ida_mem.ida_tn;
                    break;
                }
            }
        }

        /* Check for too much accuracy requested. */

        let mut nrm = IDAWrmsNorm(ida_mem, &ida_mem.ida_phi[0], &ida_mem.ida_ewt,
                                  ida_mem.ida_suppressalg);
        if ida_mem.ida_errconQ {
            nrm = IDAQuadWrmsNormUpdate(ida_mem, nrm, &ida_mem.ida_phiQ[0],
                                        &ida_mem.ida_ewtQ);
        }
        if ida_mem.ida_errconS {
            nrm = IDASensWrmsNormUpdate(ida_mem, nrm, &ida_mem.ida_phiS[0],
                                        &ida_mem.ida_ewtS, ida_mem.ida_suppressalg);
        }
        if ida_mem.ida_errconQS {
            nrm = IDAQuadSensWrmsNormUpdate(ida_mem, nrm, &ida_mem.ida_phiQS[0],
                                            &ida_mem.ida_ewtQS);
        }

        ida_mem.ida_tolsf = ida_mem.ida_uround * nrm;
        if ida_mem.ida_tolsf > ONE {
            ida_mem.ida_tolsf *= TEN;
            IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASolve", file!(),
                            &ida_msg_g(MSG_TOO_MUCH_ACC, &[ida_mem.ida_tn]));
            istate = IDA_TOO_MUCH_ACC;
            ida_mem.ida_tretlast = ida_mem.ida_tn;
            *tret = ida_mem.ida_tn;
            if ida_mem.ida_nst > 0 {
                let _ = IDAGetSolution(ida_mem, ida_mem.ida_tn, yret, ypret);
            }
            break;
        }

        /* Call IDAStep to take a step. */

        let sflag = IDAStep(ida_mem);

        /* Process all failed-step cases, and exit loop. */

        if sflag != IDA_SUCCESS {
            istate = IDAHandleFailure(ida_mem, sflag);
            ida_mem.ida_tretlast = ida_mem.ida_tn;
            *tret = ida_mem.ida_tn;
            let _ = IDAGetSolution(ida_mem, ida_mem.ida_tn, yret, ypret);
            break;
        }

        nstloc += 1;

        /* If tstop is set and was reached, reset IDA_mem->ida_tn = tstop */
        if ida_mem.ida_tstopset {
            let troundoff = HUNDRED * ida_mem.ida_uround *
                (SUNRabs(ida_mem.ida_tn) + SUNRabs(ida_mem.ida_hh));
            if SUNRabs(ida_mem.ida_tn - ida_mem.ida_tstop) <= troundoff {
                ida_mem.ida_tn = ida_mem.ida_tstop;
            }
        }

        /* After successful step, check for stop conditions; continue or break. */

        /* First check for root in the last step taken. */

        if ida_mem.ida_nrtfn > 0 {
            let ier = IDARcheck3(ida_mem, tout, itask);

            if ier == RTFOUND {
                /* A new root was found */
                ida_mem.ida_irfnd = 1;
                istate = IDA_ROOT_RETURN;
                ida_mem.ida_tretlast = ida_mem.ida_tlo;
                *tret = ida_mem.ida_tlo;
                /* ida_yy/ida_yp (interpolated root solution) alias
                   yret/ypret in C */
                yret.data.copy_from_slice(&ida_mem.ida_yy.data);
                ypret.data.copy_from_slice(&ida_mem.ida_yp.data);
                break;
            } else if ier == IDA_RTFUNC_FAIL {
                /* g failed */
                IDAProcessError(Some(ida_mem), IDA_RTFUNC_FAIL, line!(), "IDASolve", file!(),
                                &ida_msg_g(MSG_RTFUNC_FAILED, &[ida_mem.ida_tlo]));
                istate = IDA_RTFUNC_FAIL;
                yret.data.copy_from_slice(&ida_mem.ida_yy.data);
                ypret.data.copy_from_slice(&ida_mem.ida_yp.data);
                break;
            }

            /* If we are at the end of the first step and we still have
             * some event functions that are inactive, issue a warning
             * as this may indicate a user error in the implementation
             * of the root function. */

            if ida_mem.ida_nst == 1 {
                let inactive_roots = ida_mem
                    .ida_gactive
                    .iter()
                    .take(ida_mem.ida_nrtfn as usize)
                    .any(|&a| !a);
                if ida_mem.ida_mxgnull > 0 && inactive_roots {
                    IDAProcessError(Some(ida_mem), IDA_WARNING, line!(), "IDASolve", file!(),
                                    MSG_INACTIVE_ROOTS);
                }
            }
        }

        /* Now check all other stop conditions. */

        istate = IDAStopTest2(ida_mem, tout, tret, yret, ypret, itask);
        if istate != CONTINUE_STEPS {
            break;
        }
    } /* End of step loop */

    istate
}

/*
 * -----------------------------------------------------------------
 * Interpolated output and extraction functions
 * -----------------------------------------------------------------
 */

/*
 * IDAGetDky
 *
 * This routine evaluates the k-th derivative of y(t) as the value of
 * the k-th derivative of the interpolating polynomial at the independent
 * variable t, and stores the results in the vector dky.  It uses the current
 * independent variable value, tn, and the method order last used, kused.
 *
 * The return values are:
 *   IDA_SUCCESS       if t is legal
 *   IDA_BAD_T         if t is not within the interval of the last step taken
 *   IDA_BAD_DKY       if the dky vector is NULL
 *   IDA_BAD_K         if the requested k is not in the range [0,order used]
 *   IDA_VECTOROP_ERR  if the fused vector operation fails
 *   (the serial linear combination below cannot fail, so the
 *    IDA_VECTOROP_ERR branch of the C code vanishes)
 */
pub fn IDAGetDky(ida_mem: &IDAMem, t: f64, k: i32, dky: &mut NVector) -> i32 {
    if dky.is_empty() {
        IDAProcessError(Some(ida_mem), IDA_BAD_DKY, line!(), "IDAGetDky", file!(), MSG_NULL_DKY);
        return IDA_BAD_DKY;
    }

    if k < 0 || k > ida_mem.ida_kused {
        IDAProcessError(Some(ida_mem), IDA_BAD_K, line!(), "IDAGetDky", file!(), MSG_BAD_K);
        return IDA_BAD_K;
    }

    /* Check t for legality.  Here tn - hused is t_{n-1}. */

    let mut tfuzz = HUNDRED * ida_mem.ida_uround *
        (SUNRabs(ida_mem.ida_tn) + SUNRabs(ida_mem.ida_hh));
    if ida_mem.ida_hh < ZERO {
        tfuzz = -tfuzz;
    }
    let tp = ida_mem.ida_tn - ida_mem.ida_hused - tfuzz;
    if (t - tp) * ida_mem.ida_hh < ZERO {
        IDAProcessError(Some(ida_mem), IDA_BAD_T, line!(), "IDAGetDky", file!(),
                        &ida_msg_g(MSG_BAD_T,
                                   &[t, ida_mem.ida_tn - ida_mem.ida_hused, ida_mem.ida_tn]));
        return IDA_BAD_T;
    }

    /* Initialize the c_j^(k) and c_k^(k-1) */
    let mut cjk = [ZERO; MXORDP1];
    let mut cjk_1 = [ZERO; MXORDP1];

    let delt = t - ida_mem.ida_tn;

    let mut i: i32 = 0;
    while i <= k {
        /* The below recurrence is used to compute the k-th derivative of the solution:
           c_j^(k) = ( k * c_{j-1}^(k-1) + c_{j-1}^{k} (Delta+psi_{j-1}) ) / psi_j

           Translated in indexes notation:
           cjk[j] = ( k*cjk_1[j-1] + cjk[j-1]*(delt+psi[j-2]) ) / psi[j-1]

           For k=0, j=1: c_1 = c_0^(-1) + (delt+psi[-1]) / psi[0]

           In order to be able to deal with k=0 in the same way as for k>0, the
           following conventions were adopted:
             - c_0(t) = 1 , c_0^(-1)(t)=0
             - psij_1 stands for psi[-1]=0 when j=1
                             for psi[j-2]  when j>1
        */
        let mut psij_1;
        if i == 0 {
            cjk[0] = 1.0;
            psij_1 = 0.0;
        } else {
            /*                                                i       i-1          1
              c_i^(i) can be always updated since c_i^(i) = -----  --------  ... -----
                                                            psi_j  psi_{j-1}     psi_1
            */
            cjk[i as usize] = cjk[i as usize - 1] * i as f64 / ida_mem.ida_psi[i as usize - 1];
            psij_1 = ida_mem.ida_psi[i as usize - 1];
        }

        /* update c_j^(i) */

        /* j does not need to go till kused */
        let mut j: i32 = i + 1;
        while j <= ida_mem.ida_kused - k + i {
            cjk[j as usize] = (i as f64 * cjk_1[j as usize - 1]
                + cjk[j as usize - 1] * (delt + psij_1))
                / ida_mem.ida_psi[j as usize - 1];
            psij_1 = ida_mem.ida_psi[j as usize - 1];
            j += 1;
        }

        /* save existing c_j^(i)'s */
        let mut j: i32 = i + 1;
        while j <= ida_mem.ida_kused - k + i {
            cjk_1[j as usize] = cjk[j as usize];
            j += 1;
        }
        i += 1;
    }

    /* Compute sum (c_j(t) * phi(t)) */

    /* Sum j=k to j<=IDA_mem->ida_kused
       (the C code uses the fused N_VLinearCombination(kused-k+1,
       cjk+k, phi+k, dky); its serial kernel computes dky = c0*X0 and
       then accumulates dky += cj*Xj, replicated here) */
    let kk = k as usize;
    let kused = ida_mem.ida_kused as usize;
    for (d, p) in dky.data.iter_mut().zip(&ida_mem.ida_phi[kk].data) {
        *d = cjk[kk] * *p;
    }
    for j in (kk + 1)..=kused {
        for (d, p) in dky.data.iter_mut().zip(&ida_mem.ida_phi[j].data) {
            *d += cjk[j] * *p;
        }
    }

    IDA_SUCCESS
}

/*
 * IDAGetQuad
 *
 * The following function can be called to obtain the quadrature
 * variables after a successful integration step.
 *
 * This is just a wrapper that calls IDAGetQuadDky with k=0.
 */
pub fn IDAGetQuad(ida_mem: &IDAMem, ptret: &mut f64, yQout: &mut NVector) -> i32 {
    *ptret = ida_mem.ida_tretlast;

    IDAGetQuadDky(ida_mem, ida_mem.ida_tretlast, 0, yQout)
}

/*
 * IDAGetQuadDky
 *
 * Returns the quadrature variables (or their
 * derivatives up to the current method order) at any time within
 * the last integration step (dense output).
 */
pub fn IDAGetQuadDky(ida_mem: &IDAMem, t: f64, k: i32, dkyQ: &mut NVector) -> i32 {
    /* Check if quadrature was initialized */
    if ida_mem.ida_quadr != SUNTRUE {
        IDAProcessError(Some(ida_mem), IDA_NO_QUAD, line!(), "IDAGetQuadDky", file!(), MSG_NO_QUAD);
        return IDA_NO_QUAD;
    }

    if dkyQ.is_empty() {
        IDAProcessError(Some(ida_mem), IDA_BAD_DKY, line!(), "IDAGetQuadDky", file!(), MSG_NULL_DKY);
        return IDA_BAD_DKY;
    }

    /* (C-exact: this bound is ida_kk, not ida_kused as in IDAGetDky) */
    if k < 0 || k > ida_mem.ida_kk {
        IDAProcessError(Some(ida_mem), IDA_BAD_K, line!(), "IDAGetQuadDky", file!(), MSG_BAD_K);
        return IDA_BAD_K;
    }

    /* Check t for legality.  Here tn - hused is t_{n-1}. */

    /* (C quirk, preserved exactly: no SUNRabs on tn/hh and no sign
       flip for hh < 0, unlike IDAGetDky) */
    let tfuzz = HUNDRED * ida_mem.ida_uround * (ida_mem.ida_tn + ida_mem.ida_hh);
    let tp = ida_mem.ida_tn - ida_mem.ida_hused - tfuzz;
    if (t - tp) * ida_mem.ida_hh < ZERO {
        IDAProcessError(Some(ida_mem), IDA_BAD_T, line!(), "IDAGetQuadDky", file!(),
                        &ida_msg_g(MSG_BAD_T,
                                   &[t, ida_mem.ida_tn - ida_mem.ida_hused, ida_mem.ida_tn]));
        return IDA_BAD_T;
    }

    /* Initialize the c_j^(k) and c_k^(k-1) */
    let mut cjk = [ZERO; MXORDP1];
    let mut cjk_1 = [ZERO; MXORDP1];

    let delt = t - ida_mem.ida_tn;

    let mut i: i32 = 0;
    while i <= k {
        let mut psij_1;
        if i == 0 {
            cjk[0] = 1.0;
            psij_1 = 0.0;
        } else {
            cjk[i as usize] = cjk[i as usize - 1] * i as f64 / ida_mem.ida_psi[i as usize - 1];
            psij_1 = ida_mem.ida_psi[i as usize - 1];
        }

        /* update c_j^(i) */
        let mut j: i32 = i + 1;
        while j <= ida_mem.ida_kused - k + i {
            cjk[j as usize] = (i as f64 * cjk_1[j as usize - 1]
                + cjk[j as usize - 1] * (delt + psij_1))
                / ida_mem.ida_psi[j as usize - 1];
            psij_1 = ida_mem.ida_psi[j as usize - 1];
            j += 1;
        }

        /* save existing c_j^(i)'s */
        let mut j: i32 = i + 1;
        while j <= ida_mem.ida_kused - k + i {
            cjk_1[j as usize] = cjk[j as usize];
            j += 1;
        }
        i += 1;
    }

    /* Compute sum (c_j(t) * phi(t)) */

    let kk = k as usize;
    let kused = ida_mem.ida_kused as usize;
    for (d, p) in dkyQ.data.iter_mut().zip(&ida_mem.ida_phiQ[kk].data) {
        *d = cjk[kk] * *p;
    }
    for j in (kk + 1)..=kused {
        for (d, p) in dkyQ.data.iter_mut().zip(&ida_mem.ida_phiQ[j].data) {
            *d += cjk[j] * *p;
        }
    }

    IDA_SUCCESS
}

/*
 * IDAGetSens
 *
 * This routine extracts sensitivity solution into yySout at the
 * time at which IDASolve returned the solution.
 * This is just a wrapper that calls IDAGetSensDky1 with k=0 and
 * is=0, 1, ... ,NS-1.
 */
pub fn IDAGetSens(ida_mem: &IDAMem, ptret: &mut f64, yySout: &mut [NVector]) -> i32 {
    /* Check the parameters */
    if yySout.is_empty() {
        IDAProcessError(Some(ida_mem), IDA_BAD_DKY, line!(), "IDAGetSens", file!(), MSG_NULL_DKY);
        return IDA_BAD_DKY;
    }

    /* are sensitivities enabled? */
    if ida_mem.ida_sensi == SUNFALSE {
        IDAProcessError(Some(ida_mem), IDA_NO_SENS, line!(), "IDAGetSens", file!(), MSG_NO_SENSI);
        return IDA_NO_SENS;
    }

    *ptret = ida_mem.ida_tretlast;

    let mut ierr = 0;
    for is in 0..ida_mem.ida_Ns {
        ierr = IDAGetSensDky1(ida_mem, *ptret, 0, is, &mut yySout[is as usize]);
        if ierr != IDA_SUCCESS {
            break;
        }
    }

    ierr
}

/*
 * IDAGetSensDky
 *
 * Computes the k-th derivative of all sensitivities of the y function at
 * time t. It repeatedly calls IDAGetSensDky1. The argument dkyS must be
 * a pointer to N_Vector and must be allocated by the user to hold at
 * least Ns vectors.
 */
pub fn IDAGetSensDky(ida_mem: &IDAMem, t: f64, k: i32, dkySout: &mut [NVector]) -> i32 {
    if ida_mem.ida_sensi == SUNFALSE {
        IDAProcessError(Some(ida_mem), IDA_NO_SENS, line!(), "IDAGetSensDky", file!(), MSG_NO_SENSI);
        return IDA_NO_SENS;
    }

    if dkySout.is_empty() {
        IDAProcessError(Some(ida_mem), IDA_BAD_DKY, line!(), "IDAGetSensDky", file!(), MSG_NULL_DKY);
        return IDA_BAD_DKY;
    }

    /* (C-exact: bound is ida_kk here) */
    if k < 0 || k > ida_mem.ida_kk {
        IDAProcessError(Some(ida_mem), IDA_BAD_K, line!(), "IDAGetSensDky", file!(), MSG_BAD_K);
        return IDA_BAD_K;
    }

    let mut ier = 0;
    for is in 0..ida_mem.ida_Ns {
        ier = IDAGetSensDky1(ida_mem, t, k, is, &mut dkySout[is as usize]);
        if ier != IDA_SUCCESS {
            break;
        }
    }

    ier
}

/*
 * IDAGetSens1
 *
 * This routine extracts the is-th sensitivity solution into ySout
 * at the time at which IDASolve returned the solution.
 * This is just a wrapper that calls IDASensDky1 with k=0.
 *
 * (C-exact: no sensi / NULL checks here — IDAGetSensDky1 performs them.)
 */
pub fn IDAGetSens1(ida_mem: &IDAMem, ptret: &mut f64, is: i32, yySret: &mut NVector) -> i32 {
    *ptret = ida_mem.ida_tretlast;

    IDAGetSensDky1(ida_mem, *ptret, 0, is, yySret)
}

/*
 * IDAGetSensDky1
 *
 * IDASensDky1 computes the kth derivative of the yS[is] function
 * at time t, where tn-hu <= t <= tn, tn denotes the current
 * internal time reached, and hu is the last internal step size
 * successfully used by the solver. The user may request
 * is=0, 1, ..., Ns-1 and k=0, 1, ..., kk, where kk is the current
 * order. The derivative vector is returned in dky. This vector
 * must be allocated by the caller. It is only legal to call this
 * function after a successful return from IDASolve with sensitivity
 * computation enabled.
 */
pub fn IDAGetSensDky1(ida_mem: &IDAMem, t: f64, k: i32, is: i32, dkyS: &mut NVector) -> i32 {
    if ida_mem.ida_sensi == SUNFALSE {
        IDAProcessError(Some(ida_mem), IDA_NO_SENS, line!(), "IDAGetSensDky1", file!(), MSG_NO_SENSI);
        return IDA_NO_SENS;
    }

    if dkyS.is_empty() {
        IDAProcessError(Some(ida_mem), IDA_BAD_DKY, line!(), "IDAGetSensDky1", file!(), MSG_NULL_DKY);
        return IDA_BAD_DKY;
    }

    /* Is the requested sensitivity index valid? */
    if is < 0 || is >= ida_mem.ida_Ns {
        IDAProcessError(Some(ida_mem), IDA_BAD_IS, line!(), "IDAGetSensDky1", file!(), MSG_BAD_IS);
        return IDA_BAD_IS;
    }

    /* Is the requested order valid? */
    if k < 0 || k > ida_mem.ida_kused {
        IDAProcessError(Some(ida_mem), IDA_BAD_K, line!(), "IDAGetSensDky1", file!(), MSG_BAD_K);
        return IDA_BAD_K;
    }

    /* Check t for legality.  Here tn - hused is t_{n-1}. */

    let mut tfuzz = HUNDRED * ida_mem.ida_uround *
        (SUNRabs(ida_mem.ida_tn) + SUNRabs(ida_mem.ida_hh));
    if ida_mem.ida_hh < ZERO {
        tfuzz = -tfuzz;
    }
    let tp = ida_mem.ida_tn - ida_mem.ida_hused - tfuzz;
    if (t - tp) * ida_mem.ida_hh < ZERO {
        IDAProcessError(Some(ida_mem), IDA_BAD_T, line!(), "IDAGetSensDky1", file!(),
                        &ida_msg_g(MSG_BAD_T,
                                   &[t, ida_mem.ida_tn - ida_mem.ida_hused, ida_mem.ida_tn]));
        return IDA_BAD_T;
    }

    /* Initialize the c_j^(k) and c_k^(k-1) */
    let mut cjk = [ZERO; MXORDP1];
    let mut cjk_1 = [ZERO; MXORDP1];

    let delt = t - ida_mem.ida_tn;

    let mut i: i32 = 0;
    while i <= k {
        let mut psij_1;
        if i == 0 {
            cjk[0] = 1.0;
            psij_1 = 0.0;
        } else {
            cjk[i as usize] = cjk[i as usize - 1] * i as f64 / ida_mem.ida_psi[i as usize - 1];
            psij_1 = ida_mem.ida_psi[i as usize - 1];
        }

        /* Update cjk based on the recurrence */
        let mut j: i32 = i + 1;
        while j <= ida_mem.ida_kused - k + i {
            cjk[j as usize] = (i as f64 * cjk_1[j as usize - 1]
                + cjk[j as usize - 1] * (delt + psij_1))
                / ida_mem.ida_psi[j as usize - 1];
            psij_1 = ida_mem.ida_psi[j as usize - 1];
            j += 1;
        }

        /* Update cjk_1 for the next step */
        let mut j: i32 = i + 1;
        while j <= ida_mem.ida_kused - k + i {
            cjk_1[j as usize] = cjk[j as usize];
            j += 1;
        }
        i += 1;
    }

    /* Compute sum (c_j(t) * phi(t))
       (C gathers Xvecs[j-k] = phiS[j][is] and calls the fused
       N_VLinearCombination; the Xvecs staging is dropped and the
       serial kernel is applied directly to phiS[j][is].) */
    let kk = k as usize;
    let kused = ida_mem.ida_kused as usize;
    let s = is as usize;
    for (d, p) in dkyS.data.iter_mut().zip(&ida_mem.ida_phiS[kk][s].data) {
        *d = cjk[kk] * *p;
    }
    for j in (kk + 1)..=kused {
        for (d, p) in dkyS.data.iter_mut().zip(&ida_mem.ida_phiS[j][s].data) {
            *d += cjk[j] * *p;
        }
    }

    IDA_SUCCESS
}

/*
 * IDAGetQuadSens
 *
 * This routine extracts quadrature sensitivity solution into yyQSout at the
 * time at which IDASolve returned the solution.
 * This is just a wrapper that calls IDAGetQuadSensDky1 with k=0 and
 * is=0, 1, ... ,NS-1.
 */
pub fn IDAGetQuadSens(ida_mem: &IDAMem, ptret: &mut f64, yyQSout: &mut [NVector]) -> i32 {
    /* Check the parameters */
    if yyQSout.is_empty() {
        IDAProcessError(Some(ida_mem), IDA_BAD_DKY, line!(), "IDAGetQuadSens", file!(), MSG_NULL_DKY);
        return IDA_BAD_DKY;
    }

    /* are sensitivities enabled? */
    /* (C-exact quirk: returns IDA_NO_SENS — not IDA_NO_QUADSENS —
       while printing MSG_NO_QUADSENSI) */
    if ida_mem.ida_quadr_sensi == SUNFALSE {
        IDAProcessError(Some(ida_mem), IDA_NO_SENS, line!(), "IDAGetQuadSens", file!(),
                        MSG_NO_QUADSENSI);
        return IDA_NO_SENS;
    }

    *ptret = ida_mem.ida_tretlast;

    let mut ierr = 0;
    for is in 0..ida_mem.ida_Ns {
        ierr = IDAGetQuadSensDky1(ida_mem, *ptret, 0, is, &mut yyQSout[is as usize]);
        if ierr != IDA_SUCCESS {
            break;
        }
    }

    ierr
}

/*
 * IDAGetQuadSensDky
 *
 * Computes the k-th derivative of all quadratures sensitivities of the y function at
 * time t. It repeatedly calls IDAGetQuadSensDky. The argument dkyS must be
 * a pointer to N_Vector and must be allocated by the user to hold at
 * least Ns vectors.
 */
pub fn IDAGetQuadSensDky(ida_mem: &IDAMem, t: f64, k: i32, dkyQSout: &mut [NVector]) -> i32 {
    if ida_mem.ida_sensi == SUNFALSE {
        IDAProcessError(Some(ida_mem), IDA_NO_SENS, line!(), "IDAGetQuadSensDky", file!(),
                        MSG_NO_SENSI);
        return IDA_NO_SENS;
    }

    if ida_mem.ida_quadr_sensi == SUNFALSE {
        IDAProcessError(Some(ida_mem), IDA_NO_QUADSENS, line!(), "IDAGetQuadSensDky", file!(),
                        MSG_NO_QUADSENSI);
        return IDA_NO_QUADSENS;
    }

    if dkyQSout.is_empty() {
        IDAProcessError(Some(ida_mem), IDA_BAD_DKY, line!(), "IDAGetQuadSensDky", file!(),
                        MSG_NULL_DKY);
        return IDA_BAD_DKY;
    }

    /* (C-exact: bound is ida_kk here) */
    if k < 0 || k > ida_mem.ida_kk {
        IDAProcessError(Some(ida_mem), IDA_BAD_K, line!(), "IDAGetQuadSensDky", file!(), MSG_BAD_K);
        return IDA_BAD_K;
    }

    let mut ier = 0;
    for is in 0..ida_mem.ida_Ns {
        ier = IDAGetQuadSensDky1(ida_mem, t, k, is, &mut dkyQSout[is as usize]);
        if ier != IDA_SUCCESS {
            break;
        }
    }

    ier
}

/*
 * IDAGetQuadSens1
 *
 * This routine extracts the is-th quadrature sensitivity solution into yQSout
 * at the time at which IDASolve returned the solution.
 * This is just a wrapper that calls IDASensDky1 with k=0.
 */
pub fn IDAGetQuadSens1(ida_mem: &IDAMem, ptret: &mut f64, is: i32, yyQSret: &mut NVector) -> i32 {
    if ida_mem.ida_sensi == SUNFALSE {
        IDAProcessError(Some(ida_mem), IDA_NO_SENS, line!(), "IDAGetQuadSens1", file!(),
                        MSG_NO_SENSI);
        return IDA_NO_SENS;
    }

    if ida_mem.ida_quadr_sensi == SUNFALSE {
        IDAProcessError(Some(ida_mem), IDA_NO_QUADSENS, line!(), "IDAGetQuadSens1", file!(),
                        MSG_NO_QUADSENSI);
        return IDA_NO_QUADSENS;
    }

    if yyQSret.is_empty() {
        IDAProcessError(Some(ida_mem), IDA_BAD_DKY, line!(), "IDAGetQuadSens1", file!(),
                        MSG_NULL_DKY);
        return IDA_BAD_DKY;
    }

    *ptret = ida_mem.ida_tretlast;

    IDAGetQuadSensDky1(ida_mem, *ptret, 0, is, yyQSret)
}

/*
 * IDAGetQuadSensDky1
 *
 * IDAGetQuadSensDky1 computes the kth derivative of the yS[is] function
 * at time t, where tn-hu <= t <= tn, tn denotes the current
 * internal time reached, and hu is the last internal step size
 * successfully used by the solver. The user may request
 * is=0, 1, ..., Ns-1 and k=0, 1, ..., kk, where kk is the current
 * order. The derivative vector is returned in dky. This vector
 * must be allocated by the caller. It is only legal to call this
 * function after a successful return from IDASolve with sensitivity
 * computation enabled.
 */
pub fn IDAGetQuadSensDky1(ida_mem: &IDAMem, t: f64, k: i32, is: i32, dkyQS: &mut NVector) -> i32 {
    if ida_mem.ida_sensi == SUNFALSE {
        IDAProcessError(Some(ida_mem), IDA_NO_SENS, line!(), "IDAGetQuadSensDky1", file!(),
                        MSG_NO_SENSI);
        return IDA_NO_SENS;
    }

    if ida_mem.ida_quadr_sensi == SUNFALSE {
        IDAProcessError(Some(ida_mem), IDA_NO_QUADSENS, line!(), "IDAGetQuadSensDky1", file!(),
                        MSG_NO_QUADSENSI);
        return IDA_NO_QUADSENS;
    }

    if dkyQS.is_empty() {
        IDAProcessError(Some(ida_mem), IDA_BAD_DKY, line!(), "IDAGetQuadSensDky1", file!(),
                        MSG_NULL_DKY);
        return IDA_BAD_DKY;
    }

    /* Is the requested sensitivity index valid */
    if is < 0 || is >= ida_mem.ida_Ns {
        IDAProcessError(Some(ida_mem), IDA_BAD_IS, line!(), "IDAGetQuadSensDky1", file!(),
                        MSG_BAD_IS);
        return IDA_BAD_IS;
    }

    /* Is the requested order valid? */
    if k < 0 || k > ida_mem.ida_kused {
        IDAProcessError(Some(ida_mem), IDA_BAD_K, line!(), "IDAGetQuadSensDky1", file!(), MSG_BAD_K);
        return IDA_BAD_K;
    }

    /* Check t for legality.  Here tn - hused is t_{n-1}. */

    let mut tfuzz = HUNDRED * ida_mem.ida_uround *
        (SUNRabs(ida_mem.ida_tn) + SUNRabs(ida_mem.ida_hh));
    if ida_mem.ida_hh < ZERO {
        tfuzz = -tfuzz;
    }
    let tp = ida_mem.ida_tn - ida_mem.ida_hused - tfuzz;
    if (t - tp) * ida_mem.ida_hh < ZERO {
        IDAProcessError(Some(ida_mem), IDA_BAD_T, line!(), "IDAGetQuadSensDky1", file!(),
                        &ida_msg_g(MSG_BAD_T,
                                   &[t, ida_mem.ida_tn - ida_mem.ida_hused, ida_mem.ida_tn]));
        return IDA_BAD_T;
    }

    /* Initialize the c_j^(k) and c_k^(k-1) */
    let mut cjk = [ZERO; MXORDP1];
    let mut cjk_1 = [ZERO; MXORDP1];

    let delt = t - ida_mem.ida_tn;

    let mut i: i32 = 0;
    while i <= k {
        let mut psij_1;
        if i == 0 {
            cjk[0] = 1.0;
            psij_1 = 0.0;
        } else {
            cjk[i as usize] = cjk[i as usize - 1] * i as f64 / ida_mem.ida_psi[i as usize - 1];
            psij_1 = ida_mem.ida_psi[i as usize - 1];
        }

        /* Update cjk based on the recurrence */
        let mut j: i32 = i + 1;
        while j <= ida_mem.ida_kused - k + i {
            cjk[j as usize] = (i as f64 * cjk_1[j as usize - 1]
                + cjk[j as usize - 1] * (delt + psij_1))
                / ida_mem.ida_psi[j as usize - 1];
            psij_1 = ida_mem.ida_psi[j as usize - 1];
            j += 1;
        }

        /* Update cjk_1 for the next step */
        let mut j: i32 = i + 1;
        while j <= ida_mem.ida_kused - k + i {
            cjk_1[j as usize] = cjk[j as usize];
            j += 1;
        }
        i += 1;
    }

    /* Compute sum (c_j(t) * phi(t))
       (Xvecs gather dropped; serial fused kernel applied directly to
       phiQS[j][is], as in IDAGetSensDky1.) */
    let kk = k as usize;
    let kused = ida_mem.ida_kused as usize;
    let s = is as usize;
    for (d, p) in dkyQS.data.iter_mut().zip(&ida_mem.ida_phiQS[kk][s].data) {
        *d = cjk[kk] * *p;
    }
    for j in (kk + 1)..=kused {
        for (d, p) in dkyQS.data.iter_mut().zip(&ida_mem.ida_phiQS[j][s].data) {
            *d += cjk[j] * *p;
        }
    }

    IDA_SUCCESS
}

/*
 * IDAComputeY
 *
 * Computes y based on the current prediction and given correction.
 */
pub fn IDAComputeY(ida_mem: &IDAMem, ycor: &NVector, y: &mut NVector) -> i32 {
    N_VLinearSum(ONE, &ida_mem.ida_yypredict, ONE, ycor, y);

    IDA_SUCCESS
}

/*
 * IDAComputeYp
 *
 * Computes y' based on the current prediction and given correction.
 */
pub fn IDAComputeYp(ida_mem: &IDAMem, ycor: &NVector, yp: &mut NVector) -> i32 {
    N_VLinearSum(ONE, &ida_mem.ida_yppredict, ida_mem.ida_cj, ycor, yp);

    IDA_SUCCESS
}

/*
 * IDAComputeYSens
 *
 * Computes yS based on the current prediction and given correction.
 */
pub fn IDAComputeYSens(ida_mem: &IDAMem, ycorS: &[NVector], yyS: &mut [NVector]) -> i32 {
    /* (C: N_VLinearSumVectorArray(Ns, ...); the serial fused kernel is
       a per-vector N_VLinearSum loop, replicated here; the cannot-fail
       retval branch vanishes) */
    for is in 0..(ida_mem.ida_Ns as usize) {
        N_VLinearSum(ONE, &ida_mem.ida_yySpredict[is], ONE, &ycorS[is], &mut yyS[is]);
    }

    IDA_SUCCESS
}

/*
 * IDAComputeYpSens
 *
 * Computes yS' based on the current prediction and given correction.
 */
pub fn IDAComputeYpSens(ida_mem: &IDAMem, ycorS: &[NVector], ypS: &mut [NVector]) -> i32 {
    for is in 0..(ida_mem.ida_Ns as usize) {
        N_VLinearSum(ONE, &ida_mem.ida_ypSpredict[is], ida_mem.ida_cj, &ycorS[is], &mut ypS[is]);
    }

    IDA_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Deallocation functions
 * -----------------------------------------------------------------
 */

/*
 * IDAFree
 *
 * This routine frees the problem memory allocated by IDAInit.
 * In C this releases the vectors (IDAFreeVectors), the quadrature,
 * sensitivity, quadrature-sensitivity and adjoint memory
 * (IDAQuadFree / IDASensFree / IDAQuadSensFree / IDAAdjFree), the
 * owned NLS objects (NLS, NLSsim, NLSstg), the linear solver memory
 * (lfree), the rootfinding arrays and the fused-op work arrays; here
 * all of that is RAII — dropping the Box releases everything.
 * (Donor ida.rs signature, which the examples call.)
 */
pub fn IDAFree(_ida_mem: Box<IDAMem>) {}

/*
 * IDAQuadFree
 *
 * IDAQuadFree frees the problem memory in ida_mem allocated
 * for quadrature integration. Its only argument is the pointer
 * ida_mem returned by IDACreate.
 *
 * (Unlike IDAFree this is a mid-lifetime operation in C — the
 * integrator remains usable afterwards — so the flag resets and the
 * IDAQuadFreeVectors lrw/liw bookkeeping are observable and are
 * ported for real.)
 */
pub fn IDAQuadFree(ida_mem: &mut IDAMem) {
    if ida_mem.ida_quadMallocDone {
        IDAQuadFreeVectors(ida_mem);
        ida_mem.ida_quadMallocDone = SUNFALSE;
        ida_mem.ida_quadr = SUNFALSE;
    }
}

/*
 * IDASensFree
 *
 * IDASensFree frees the problem memory in ida_mem allocated
 * for sensitivity analysis. Its only argument is the pointer
 * ida_mem returned by IDACreate.
 */
pub fn IDASensFree(ida_mem: &mut IDAMem) {
    if ida_mem.ida_sensMallocDone {
        IDASensFreeVectors(ida_mem);
        ida_mem.ida_sensMallocDone = SUNFALSE;
        ida_mem.ida_sensi = SUNFALSE;
    }

    /* free any vector wrappers */
    /* (C destroys the ypredictSim/ycorSim/ewtSim and ypredictStg/
       ycorStg/ewtStg senswrapper aliases; per the pinned idas_impl
       convention those wrappers are not stored as IDAMem fields — the
       sim/stg NLS state lives in the NLS objects — so only the flags
       reset here) */
    if ida_mem.simMallocDone {
        ida_mem.simMallocDone = SUNFALSE;
    }
    if ida_mem.stgMallocDone {
        ida_mem.stgMallocDone = SUNFALSE;
    }

    /* if IDA created the NLS object then free it */
    if ida_mem.ownNLSsim {
        ida_mem.NLSsim = None;
        ida_mem.ownNLSsim = SUNFALSE;
    }
    if ida_mem.ownNLSstg {
        ida_mem.NLSstg = None;
        ida_mem.ownNLSstg = SUNFALSE;
    }

    /* free min atol array if necessary */
    if !ida_mem.ida_atolSmin0.is_empty() {
        ida_mem.ida_atolSmin0 = Vec::new();
    }
}

/*
 * IDAQuadSensFree
 *
 * IDAQuadSensFree frees the problem memory in ida_mem allocated
 * for quadrature sensitivity analysis. Its only argument is the
 * pointer ida_mem returned by IDACreate.
 */
pub fn IDAQuadSensFree(ida_mem: &mut IDAMem) {
    if ida_mem.ida_quadSensMallocDone {
        IDAQuadSensFreeVectors(ida_mem);
        ida_mem.ida_quadSensMallocDone = SUNFALSE;
        ida_mem.ida_quadr_sensi = SUNFALSE;
    }

    /* free min atol array if necessary */
    if !ida_mem.ida_atolQSmin0.is_empty() {
        ida_mem.ida_atolQSmin0 = Vec::new();
    }
}

/*
 * =================================================================
 * PRIVATE FUNCTIONS
 * =================================================================
 */

/*
 * IDACheckNvector
 *
 * This routine checks if all required vector operations are present.
 * If any of them is missing it returns SUNFALSE.
 * (The serial NVector implements every required operation, so this
 * always returns SUNTRUE; kept for structural fidelity.)
 */
fn IDACheckNvector(_tmpl: &NVector) -> bool {
    SUNTRUE
}

/*
 * -----------------------------------------------------------------
 * Memory allocation/deallocation
 * -----------------------------------------------------------------
 */

/*
 * IDAAllocVectors
 *
 * This routine allocates the IDA vectors ewt, tempv1, tempv2, and
 * phi[0], ..., phi[maxord].
 * (Allocation is infallible in Rust, so the C SUNFALSE/rollback
 * branches vanish and no boolean is returned.)
 * This routine also sets the optional outputs lrw and liw, which are
 * (respectively) the lengths of the real and integer work spaces
 * allocated here.
 *
 * Adaptation (donor ida.rs): ida_yy and ida_yp — which in C are never
 * allocated because they alias the user's yret/ypret during IDASolve
 * (and IDACalcIC's yy0/yp0) — are allocated here as owned stand-ins.
 * They are deliberately NOT counted in lrw/liw, matching the C
 * workspace accounting.
 */
fn IDAAllocVectors(ida_mem: &mut IDAMem, tmpl: &NVector) {
    /* Allocate ewt, ee, delta, yypredict, yppredict, savres, tempv1, tempv2, tempv3 */

    ida_mem.ida_ewt = N_VClone(tmpl);
    ida_mem.ida_ee = N_VClone(tmpl);
    ida_mem.ida_delta = N_VClone(tmpl);
    ida_mem.ida_yypredict = N_VClone(tmpl);
    ida_mem.ida_yppredict = N_VClone(tmpl);
    ida_mem.ida_savres = N_VClone(tmpl);
    ida_mem.ida_tempv1 = N_VClone(tmpl);
    ida_mem.ida_tempv2 = N_VClone(tmpl);
    ida_mem.ida_tempv3 = N_VClone(tmpl);

    /* (adaptation — owned stand-ins for the C user-vector aliases) */
    ida_mem.ida_yy = N_VClone(tmpl);
    ida_mem.ida_yp = N_VClone(tmpl);

    /* Allocate phi[0] ... phi[maxord].  Make sure phi[2] and phi[3] are
       allocated (for use as temporary vectors), regardless of maxord.  */

    let maxcol = if ida_mem.ida_maxord > 3 { ida_mem.ida_maxord } else { 3 };
    ida_mem.ida_phi = (0..=maxcol).map(|_| N_VClone(tmpl)).collect();

    /* Update solver workspace lengths */
    ida_mem.ida_lrw += (maxcol as i64 + 10) * ida_mem.ida_lrw1;
    ida_mem.ida_liw += (maxcol as i64 + 10) * ida_mem.ida_liw1;

    /* Store the value of maxord used here */
    ida_mem.ida_maxord_alloc = ida_mem.ida_maxord;
}

/* IDAFreeVectors has no Rust counterpart: it is called only from
   IDAFree (and failed-init paths that vanish), the vectors are owned
   fields released by drop, and its lrw/liw bookkeeping is only
   observable through IDAGetWorkSpace before the memory is destroyed
   (donor ida.rs convention).  The Quad/Sens/QuadSens FreeVectors
   helpers below ARE real ports: their callers (IDAQuadFree,
   IDASensFree, IDAQuadSensFree) are mid-lifetime API. */

/*
 * IDAQuadAllocVectors
 *
 * NOTE: Space for ewtQ is allocated even when errconQ=SUNFALSE,
 * although in this case, ewtQ is never used. The reason for this
 * decision is to allow the user to re-initialize the quadrature
 * computation with errconQ=SUNTRUE, after an initialization with
 * errconQ=SUNFALSE, without new memory allocation within
 * IDAQuadReInit.
 */
fn IDAQuadAllocVectors(ida_mem: &mut IDAMem, tmpl: &NVector) {
    /* Allocate yyQ, ypQ, ewtQ, eeQ (infallible; the C rollback
       branches vanish) */
    ida_mem.ida_yyQ = N_VClone(tmpl);
    ida_mem.ida_ypQ = N_VClone(tmpl);
    ida_mem.ida_ewtQ = N_VClone(tmpl);
    ida_mem.ida_eeQ = N_VClone(tmpl);

    ida_mem.ida_phiQ = (0..=ida_mem.ida_maxord).map(|_| N_VClone(tmpl)).collect();

    ida_mem.ida_lrw += (ida_mem.ida_maxord as i64 + 4) * ida_mem.ida_lrw1Q;
    ida_mem.ida_liw += (ida_mem.ida_maxord as i64 + 4) * ida_mem.ida_liw1Q;
}

/*
 * IDAQuadFreeVectors
 *
 * This routine frees the IDAS vectors allocated in IDAQuadAllocVectors.
 */
fn IDAQuadFreeVectors(ida_mem: &mut IDAMem) {
    ida_mem.ida_yyQ = NVector::default();
    ida_mem.ida_ypQ = NVector::default();
    ida_mem.ida_ewtQ = NVector::default();
    ida_mem.ida_eeQ = NVector::default();
    ida_mem.ida_phiQ = Vec::new();

    /* (C-exact asymmetry, preserved: IDAQuadAllocVectors adds
       (maxord+4)*lrw1Q but this subtracts (maxord+5)*lrw1Q) */
    ida_mem.ida_lrw -= (ida_mem.ida_maxord as i64 + 5) * ida_mem.ida_lrw1Q;
    ida_mem.ida_liw -= (ida_mem.ida_maxord as i64 + 5) * ida_mem.ida_liw1Q;

    if ida_mem.ida_VatolQMallocDone {
        ida_mem.ida_VatolQ = NVector::default();
        ida_mem.ida_lrw -= ida_mem.ida_lrw1Q;
        ida_mem.ida_liw -= ida_mem.ida_liw1Q;
    }

    ida_mem.ida_VatolQMallocDone = SUNFALSE;
}

/*
 * IDASensAllocVectors
 *
 * Allocates space for the N_Vectors, plist, and pbar required for FSA.
 * (Allocation is infallible; the C rollback branches vanish and no
 * boolean is returned.)
 */
fn IDASensAllocVectors(ida_mem: &mut IDAMem, tmpl: &NVector) {
    /* (C: ida_tmpS1 = ida_tempv1; ida_tmpS2 = ida_tempv2 — plain
       aliases; per the pinned convention aliases are not stored as
       fields, and call sites use ida_tempv1/ida_tempv2 directly) */

    let ns = ida_mem.ida_Ns as usize;

    /* Allocate space for workspace vectors
       (N_VCloneVectorArray → map/collect) */

    ida_mem.ida_tmpS3 = N_VClone(tmpl);
    ida_mem.ida_ewtS = (0..ns).map(|_| N_VClone(tmpl)).collect();
    ida_mem.ida_eeS = (0..ns).map(|_| N_VClone(tmpl)).collect();
    ida_mem.ida_yyS = (0..ns).map(|_| N_VClone(tmpl)).collect();
    ida_mem.ida_ypS = (0..ns).map(|_| N_VClone(tmpl)).collect();
    ida_mem.ida_yySpredict = (0..ns).map(|_| N_VClone(tmpl)).collect();
    ida_mem.ida_ypSpredict = (0..ns).map(|_| N_VClone(tmpl)).collect();
    ida_mem.ida_deltaS = (0..ns).map(|_| N_VClone(tmpl)).collect();

    /* Update solver workspace lengths */
    /* (C-exact quirk, preserved: 7 arrays of Ns vectors plus tmpS3 are
       allocated above, but the C code counts only (5*Ns+1)*lrw1) */
    ida_mem.ida_lrw += (5 * ida_mem.ida_Ns as i64 + 1) * ida_mem.ida_lrw1;
    ida_mem.ida_liw += (5 * ida_mem.ida_Ns as i64 + 1) * ida_mem.ida_liw1;

    /* Allocate space for phiS */
    /*  Make sure phiS[2], phiS[3] and phiS[4] are
        allocated (for use as temporary vectors), regardless of maxord.*/

    let maxcol = if ida_mem.ida_maxord > 4 { ida_mem.ida_maxord } else { 4 };
    ida_mem.ida_phiS = (0..=maxcol)
        .map(|_| (0..ns).map(|_| N_VClone(tmpl)).collect())
        .collect();

    /* Update solver workspace lengths */
    /* (C-exact quirk, preserved: maxcol+1 rows are allocated above,
       but the C code counts only maxcol*Ns*lrw1) */
    ida_mem.ida_lrw += maxcol as i64 * ida_mem.ida_Ns as i64 * ida_mem.ida_lrw1;
    ida_mem.ida_liw += maxcol as i64 * ida_mem.ida_Ns as i64 * ida_mem.ida_liw1;

    /* Allocate space for pbar and plist */

    ida_mem.ida_pbar = vec![ZERO; ns];
    ida_mem.ida_plist = vec![0; ns];

    /* Update solver workspace lengths */
    ida_mem.ida_lrw += ida_mem.ida_Ns as i64;
    ida_mem.ida_liw += ida_mem.ida_Ns as i64;
}

/*
 * IDASensFreeVectors
 *
 * Frees memory allocated by IDASensAllocVectors.
 */
fn IDASensFreeVectors(ida_mem: &mut IDAMem) {
    ida_mem.ida_deltaS = Vec::new();
    ida_mem.ida_ypSpredict = Vec::new();
    ida_mem.ida_yySpredict = Vec::new();
    ida_mem.ida_ypS = Vec::new();
    ida_mem.ida_yyS = Vec::new();
    ida_mem.ida_eeS = Vec::new();
    ida_mem.ida_ewtS = Vec::new();
    ida_mem.ida_tmpS3 = NVector::default();

    /* (maxord_alloc here, unlike IDAQuadSensFreeVectors which uses
       maxord — C-exact) */
    let maxcol =
        (if ida_mem.ida_maxord_alloc > 4 { ida_mem.ida_maxord_alloc } else { 4 }) as i64;
    ida_mem.ida_phiS = Vec::new();

    ida_mem.ida_pbar = Vec::new();
    ida_mem.ida_plist = Vec::new();

    let ns = ida_mem.ida_Ns as i64;
    ida_mem.ida_lrw -= ((maxcol + 3) * ns + 1) * ida_mem.ida_lrw1 + ns;
    ida_mem.ida_liw -= ((maxcol + 3) * ns + 1) * ida_mem.ida_liw1 + ns;

    if ida_mem.ida_VatolSMallocDone {
        ida_mem.ida_VatolS = Vec::new();
        ida_mem.ida_lrw -= ns * ida_mem.ida_lrw1;
        ida_mem.ida_liw -= ns * ida_mem.ida_liw1;
        ida_mem.ida_VatolSMallocDone = SUNFALSE;
    }
    if ida_mem.ida_SatolSMallocDone {
        ida_mem.ida_SatolS = Vec::new();
        /* (C-exact: SatolS is counted in lrw only) */
        ida_mem.ida_lrw -= ns;
        ida_mem.ida_SatolSMallocDone = SUNFALSE;
    }
}

/*
 * IDAQuadSensAllocVectors
 *
 * Create (through duplication) N_Vectors used for quadrature sensitivity analysis,
 * using the N_Vector 'tmpl' as a template.
 */
fn IDAQuadSensAllocVectors(ida_mem: &mut IDAMem, tmpl: &NVector) {
    let ns = ida_mem.ida_Ns as usize;

    /* Allocate yQS, ewtQS, tempvQS, eeQS */
    ida_mem.ida_yyQS = (0..ns).map(|_| N_VClone(tmpl)).collect();
    ida_mem.ida_ewtQS = (0..ns).map(|_| N_VClone(tmpl)).collect();
    ida_mem.ida_tempvQS = (0..ns).map(|_| N_VClone(tmpl)).collect();
    ida_mem.ida_eeQS = (0..ns).map(|_| N_VClone(tmpl)).collect();

    /* (C bug, vanishes here: the savrhsQ failure branch destroys the
       arrays above but is missing its `return (SUNFALSE);` —
       allocation is infallible in Rust anyway) */
    ida_mem.ida_savrhsQ = N_VClone(tmpl);

    let maxcol = if ida_mem.ida_maxord > 4 { ida_mem.ida_maxord } else { 4 };
    /* Allocate phiQS */
    ida_mem.ida_phiQS = (0..=maxcol)
        .map(|_| (0..ns).map(|_| N_VClone(tmpl)).collect())
        .collect();

    /* Update solver workspace lengths */
    ida_mem.ida_lrw += (maxcol as i64 + 5) * ida_mem.ida_Ns as i64 * ida_mem.ida_lrw1Q;
    ida_mem.ida_liw += (maxcol as i64 + 5) * ida_mem.ida_Ns as i64 * ida_mem.ida_liw1Q;
}

/*
 * IDAQuadSensFreeVectors
 *
 * This routine frees the IDAS vectors allocated in IDAQuadSensAllocVectors.
 */
fn IDAQuadSensFreeVectors(ida_mem: &mut IDAMem) {
    /* (C-exact: maxord here, NOT maxord_alloc as in
       IDASensFreeVectors) */
    let maxcol = (if ida_mem.ida_maxord > 4 { ida_mem.ida_maxord } else { 4 }) as i64;

    ida_mem.ida_yyQS = Vec::new();
    ida_mem.ida_ewtQS = Vec::new();
    ida_mem.ida_eeQS = Vec::new();
    ida_mem.ida_tempvQS = Vec::new();
    ida_mem.ida_savrhsQ = NVector::default();

    ida_mem.ida_phiQS = Vec::new();

    let ns = ida_mem.ida_Ns as i64;
    ida_mem.ida_lrw -= (maxcol + 5) * ns * ida_mem.ida_lrw1Q;
    ida_mem.ida_liw -= (maxcol + 5) * ns * ida_mem.ida_liw1Q;

    if ida_mem.ida_VatolQSMallocDone {
        ida_mem.ida_VatolQS = Vec::new();
        ida_mem.ida_lrw -= ns * ida_mem.ida_lrw1Q;
        ida_mem.ida_liw -= ns * ida_mem.ida_liw1Q;
    }
    if ida_mem.ida_SatolQSMallocDone {
        ida_mem.ida_SatolQS = Vec::new();
        /* (C-exact: SatolQS is counted in lrw only) */
        ida_mem.ida_lrw -= ns;
    }
    /* (C-exact: both flags reset unconditionally at the end, unlike
       IDASensFreeVectors which resets them inside the conditionals) */
    ida_mem.ida_VatolQSMallocDone = SUNFALSE;
    ida_mem.ida_SatolQSMallocDone = SUNFALSE;
}

/*
 * -----------------------------------------------------------------
 * Initial setup
 * -----------------------------------------------------------------
 */

/*
 * IDAInitialSetup
 *
 * This routine is called by IDASolve once at the first step.
 * It performs all checks on optional inputs and inputs to
 * IDAInit/IDAReInit that could not be done before.
 *
 * If no error is encountered, IDAInitialSetup returns IDA_SUCCESS.
 * Otherwise, it returns an error flag and reported to the error
 * handler function.
 *
 * (Non-static in C: also called by IDACalcIC in idas_ic.c.)
 */
pub fn IDAInitialSetup(ida_mem: &mut IDAMem) -> i32 {
    /* Test for more vector operations, depending on options
       (nvwrmsnormmask is always implemented by the serial NVector,
       so the suppressalg MSG_BAD_NVECTOR branch vanishes) */

    /* Test id vector for legality */
    if ida_mem.ida_suppressalg && !ida_mem.ida_idMallocDone {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAInitialSetup", file!(),
                        MSG_MISSING_ID);
        return IDA_ILL_INPUT;
    }

    /* Did the user specify tolerances? */
    if ida_mem.ida_itol == IDA_NN {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAInitialSetup", file!(),
                        MSG_NO_TOLS);
        return IDA_ILL_INPUT;
    }

    /* Set data for efun: collapsed into ida_efun_dispatch (C sets
       ida_edata = user_data or IDA_mem here) */

    /* Initial error weight vector */
    let ier = ida_efun_apply_to_ewt(ida_mem);
    if ier != 0 {
        if ida_mem.ida_itol == IDA_WF {
            IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAInitialSetup", file!(),
                            MSG_FAIL_EWT);
        } else {
            IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAInitialSetup", file!(),
                            MSG_BAD_EWT);
        }
        return IDA_ILL_INPUT;
    }

    if ida_mem.ida_quadr {
        /* Evaluate quadrature rhs and set phiQ[1] */
        let ier = {
            let rhsQ = ida_mem.ida_rhsQ.unwrap();
            let tn = ida_mem.ida_tn;
            let IDAMem { ida_phi, ida_phiQ, ida_user_data, .. } = ida_mem;
            rhsQ(tn, &ida_phi[0], &ida_phi[1], &mut ida_phiQ[1], ida_user_data)
        };
        ida_mem.ida_nrQe += 1;
        if ier < 0 {
            IDAProcessError(Some(ida_mem), IDA_QRHS_FAIL, line!(), "IDAInitialSetup", file!(),
                            MSG_QRHSFUNC_FAILED);
            return IDA_QRHS_FAIL;
        } else if ier > 0 {
            IDAProcessError(Some(ida_mem), IDA_FIRST_QRHS_ERR, line!(), "IDAInitialSetup",
                            file!(), MSG_QRHSFUNC_FIRST);
            return IDA_FIRST_QRHS_ERR;
        }

        if ida_mem.ida_errconQ {
            /* Did the user specify tolerances? */
            if ida_mem.ida_itolQ == IDA_NN {
                IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAInitialSetup",
                                file!(), MSG_NO_TOLQ);
                return IDA_ILL_INPUT;
            }

            /* Load ewtQ */
            let mut wQ = std::mem::take(&mut ida_mem.ida_ewtQ);
            let ier = IDAQuadEwtSet(ida_mem, &ida_mem.ida_phiQ[0], &mut wQ);
            ida_mem.ida_ewtQ = wQ;
            if ier != 0 {
                IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAInitialSetup",
                                file!(), MSG_BAD_EWTQ);
                return IDA_ILL_INPUT;
            }
        }
    } else {
        ida_mem.ida_errconQ = SUNFALSE;
    }

    if ida_mem.ida_sensi {
        /* Did the user specify tolerances? */
        if ida_mem.ida_itolS == IDA_NN {
            IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAInitialSetup", file!(),
                            MSG_NO_TOLS);
            return IDA_ILL_INPUT;
        }

        /* Load ewtS */
        let mut wS = std::mem::take(&mut ida_mem.ida_ewtS);
        let yS0 = std::mem::take(&mut ida_mem.ida_phiS[0]);
        let ier = IDASensEwtSet(ida_mem, &yS0, &mut wS);
        ida_mem.ida_phiS[0] = yS0;
        ida_mem.ida_ewtS = wS;
        if ier != 0 {
            IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAInitialSetup", file!(),
                            MSG_BAD_EWTS);
            return IDA_ILL_INPUT;
        }
    } else {
        ida_mem.ida_errconS = SUNFALSE;
    }

    if ida_mem.ida_quadr_sensi {
        /* store the quadrature sensitivity residual. */
        /* (C call: rhsQS(Ns, tn, phi[0], phi[1], phiS[0], phiS[1],
           phiQ[1], phiQS[1], user_dataQS, tmpS1, tmpS2, tmpS3) with
           tmpS1/tmpS2 aliasing tempv1/tempv2.  The vectors are
           detached so the internal DQ path can borrow IDAMem mutably
           (it perturbs ida_p); user_dataQS self-pointer is not stored
           per the pinned convention — the rhsQSDQ flag selects the
           path and the user path receives ida_user_data directly.) */
        let ns = ida_mem.ida_Ns;
        let tn = ida_mem.ida_tn;
        let phi0 = std::mem::take(&mut ida_mem.ida_phi[0]);
        let phi1 = std::mem::take(&mut ida_mem.ida_phi[1]);
        let phiS0 = std::mem::take(&mut ida_mem.ida_phiS[0]);
        let phiS1 = std::mem::take(&mut ida_mem.ida_phiS[1]);
        let phiQ1 = std::mem::take(&mut ida_mem.ida_phiQ[1]);
        let mut phiQS1 = std::mem::take(&mut ida_mem.ida_phiQS[1]);
        let mut tempv1 = std::mem::take(&mut ida_mem.ida_tempv1);
        let mut tempv2 = std::mem::take(&mut ida_mem.ida_tempv2);
        let mut tmpS3 = std::mem::take(&mut ida_mem.ida_tmpS3);

        let ier = if ida_mem.ida_rhsQSDQ {
            IDAQuadSensRhsInternalDQ(ida_mem, ns, tn, &phi0, &phi1, &phiS0, &phiS1, &phiQ1,
                                     &mut phiQS1, &mut tempv1, &mut tempv2, &mut tmpS3)
        } else {
            let rhsQS = ida_mem.ida_rhsQS.unwrap();
            rhsQS(ns, tn, &phi0, &phi1, &phiS0, &phiS1, &phiQ1, &mut phiQS1,
                  &mut ida_mem.ida_user_data, &mut tempv1, &mut tempv2, &mut tmpS3)
        };

        ida_mem.ida_phi[0] = phi0;
        ida_mem.ida_phi[1] = phi1;
        ida_mem.ida_phiS[0] = phiS0;
        ida_mem.ida_phiS[1] = phiS1;
        ida_mem.ida_phiQ[1] = phiQ1;
        ida_mem.ida_phiQS[1] = phiQS1;
        ida_mem.ida_tempv1 = tempv1;
        ida_mem.ida_tempv2 = tempv2;
        ida_mem.ida_tmpS3 = tmpS3;

        ida_mem.ida_nrQSe += 1;
        if ier < 0 {
            /* (C-exact quirk, preserved: the error logged is
               IDA_QSRHS_FAIL but the value returned is
               IDA_QRHS_FAIL) */
            IDAProcessError(Some(ida_mem), IDA_QSRHS_FAIL, line!(), "IDAInitialSetup", file!(),
                            MSG_QSRHSFUNC_FAILED);
            return IDA_QRHS_FAIL;
        } else if ier > 0 {
            IDAProcessError(Some(ida_mem), IDA_FIRST_QSRHS_ERR, line!(), "IDAInitialSetup",
                            file!(), MSG_QSRHSFUNC_FIRST);
            return IDA_FIRST_QSRHS_ERR;
        }

        /* If using the internal DQ functions, we must have access to fQ
         * (i.e. quadrature integration must be enabled) and to the problem parameters */

        if ida_mem.ida_rhsQSDQ {
            /* Test if quadratures are defined, so we can use fQ */
            if !ida_mem.ida_quadr {
                IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAInitialSetup",
                                file!(), MSG_NULL_RHSQ);
                return IDA_ILL_INPUT;
            }

            /* Test if we have the problem parameters */
            if ida_mem.ida_p.is_empty() {
                IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAInitialSetup",
                                file!(), MSG_NULL_P);
                return IDA_ILL_INPUT;
            }
        }

        if ida_mem.ida_errconQS {
            /* Did the user specify tolerances? */
            if ida_mem.ida_itolQS == IDA_NN {
                IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAInitialSetup",
                                file!(), MSG_NO_TOLQS);
                return IDA_ILL_INPUT;
            }

            /* If needed, did the user provide quadrature tolerances? */
            if ida_mem.ida_itolQS == IDA_EE && ida_mem.ida_itolQ == IDA_NN {
                IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAInitialSetup",
                                file!(), MSG_NO_TOLQ);
                return IDA_ILL_INPUT;
            }

            /* Load ewtS */
            let mut wQS = std::mem::take(&mut ida_mem.ida_ewtQS);
            let ier = IDAQuadSensEwtSet(ida_mem, &ida_mem.ida_phiQS[0], &mut wQS);
            ida_mem.ida_ewtQS = wQS;
            if ier != 0 {
                IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAInitialSetup",
                                file!(), MSG_BAD_EWTQS);
                return IDA_ILL_INPUT;
            }
        }
    } else {
        ida_mem.ida_errconQS = SUNFALSE;
    }

    /* Check to see if y0 satisfies constraints. */
    if ida_mem.ida_constraintsSet {
        if ida_mem.ida_sensi && ida_mem.ida_ism == IDA_SIMULTANEOUS {
            IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAInitialSetup", file!(),
                            MSG_BAD_ISM_CONSTR);
            return IDA_ILL_INPUT;
        }

        let conOK = N_VConstrMask(&ida_mem.ida_constraints, &ida_mem.ida_phi[0],
                                  &mut ida_mem.ida_tempv2);
        if !conOK {
            IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAInitialSetup", file!(),
                            MSG_Y0_FAIL_CONSTR);
            return IDA_ILL_INPUT;
        }
    }

    /* Call linit function if it exists. */
    if !ida_mem.ida_lmem.is_none() {
        let ier = ida_linit_dispatch(ida_mem);
        if ier != 0 {
            IDAProcessError(Some(ida_mem), IDA_LINIT_FAIL, line!(), "IDAInitialSetup", file!(),
                            MSG_LINIT_FAIL);
            return IDA_LINIT_FAIL;
        }
    }

    /* Initialize the nonlinear solver (must occur after linear solver is initialize) so
     * that lsetup and lsolve pointers have been set */

    /* always initialize the DAE NLS in case the user disables sensitivities later */
    let ier = idaNlsInit(ida_mem);
    if ier != IDA_SUCCESS {
        IDAProcessError(Some(ida_mem), IDA_NLS_INIT_FAIL, line!(), "IDAInitialSetup", file!(),
                        MSG_NLS_INIT_FAIL);
        return IDA_NLS_INIT_FAIL;
    }

    if ida_mem.NLSsim.is_some() {
        let ier = idaNlsInitSensSim(ida_mem);
        if ier != IDA_SUCCESS {
            IDAProcessError(Some(ida_mem), IDA_NLS_INIT_FAIL, line!(), "IDAInitialSetup",
                            file!(), MSG_NLS_INIT_FAIL);
            return IDA_NLS_INIT_FAIL;
        }
    }

    if ida_mem.NLSstg.is_some() {
        let ier = idaNlsInitSensStg(ida_mem);
        if ier != IDA_SUCCESS {
            IDAProcessError(Some(ida_mem), IDA_NLS_INIT_FAIL, line!(), "IDAInitialSetup",
                            file!(), MSG_NLS_INIT_FAIL);
            return IDA_NLS_INIT_FAIL;
        }
    }

    IDA_SUCCESS
}

/*
 * IDAEwtSet
 *
 * This routine is responsible for loading the error weight vector
 * ewt, according to itol, as follows:
 * (1) ewt[i] = 1 / (rtol * SUNRabs(ycur[i]) + atol), i=0,...,Neq-1
 *     if itol = IDA_SS
 * (2) ewt[i] = 1 / (rtol * SUNRabs(ycur[i]) + atol[i]), i=0,...,Neq-1
 *     if itol = IDA_SV
 *
 *  IDAEwtSet returns 0 if ewt is successfully set as above to a
 *  positive vector and -1 otherwise. In the latter case, ewt is
 *  considered undefined.
 *
 * All the real work is done in the routines IDAEwtSetSS, IDAEwtSetSV.
 * (C's efun signature (ycur, weight, data) with data = IDA_mem is
 * realized by ida_efun_dispatch; itol == IDA_WF never reaches here.)
 */
pub fn IDAEwtSet(ida_mem: &IDAMem, ycur: &NVector, weight: &mut NVector) -> i32 {
    let mut flag = 0;

    match ida_mem.ida_itol {
        IDA_SS => flag = IDAEwtSetSS(ida_mem, ycur, weight),
        IDA_SV => flag = IDAEwtSetSV(ida_mem, ycur, weight),
        _ => {}
    }
    flag
}

/*
 * IDAEwtSetSS
 *
 * This routine sets ewt as described above in the case itol=IDA_SS.
 * If the absolute tolerance is zero, it tests for non-positive components
 * before inverting. IDAEwtSetSS returns 0 if ewt is successfully set to a
 * positive vector and -1 otherwise. In the latter case, ewt is considered
 * undefined.
 *
 * (C uses ida_tempv1 as scratch; a local is behaviorally identical
 * elementwise and keeps the &IDAMem signature.  The N_V calls are
 * reproduced as their serial kernels, one op per loop, preserving
 * FP order.)
 */
fn IDAEwtSetSS(ida_mem: &IDAMem, ycur: &NVector, weight: &mut NVector) -> i32 {
    let mut tempv = N_VClone(ycur);
    /* N_VAbs(ycur, tempv1) */
    for (t, y) in tempv.data.iter_mut().zip(&ycur.data) {
        *t = SUNRabs(*y);
    }
    /* N_VScale(rtol, tempv1, tempv1) */
    let rtol = ida_mem.ida_rtol;
    for t in tempv.data.iter_mut() {
        *t = rtol * *t;
    }
    /* N_VAddConst(tempv1, Satol, tempv1) */
    let satol = ida_mem.ida_Satol;
    for t in tempv.data.iter_mut() {
        *t += satol;
    }
    if ida_mem.ida_atolmin0 {
        if N_VMin(&tempv) <= ZERO {
            return -1;
        }
    }
    /* N_VInv(tempv1, weight) */
    for (w, t) in weight.data.iter_mut().zip(&tempv.data) {
        *w = ONE / *t;
    }
    0
}

/*
 * IDAEwtSetSV
 *
 * This routine sets ewt as described above in the case itol=IDA_SV.
 * If the absolute tolerance is zero, it tests for non-positive components
 * before inverting. IDAEwtSetSV returns 0 if ewt is successfully set to a
 * positive vector and -1 otherwise. In the latter case, ewt is considered
 * undefined.
 */
fn IDAEwtSetSV(ida_mem: &IDAMem, ycur: &NVector, weight: &mut NVector) -> i32 {
    let mut tempv = N_VClone(ycur);
    /* N_VAbs(ycur, tempv1) */
    for (t, y) in tempv.data.iter_mut().zip(&ycur.data) {
        *t = SUNRabs(*y);
    }
    /* N_VLinearSum(rtol, tempv1, ONE, Vatol, tempv1) — serial kernel
       VLin1: z[i] = (a*x[i]) + y[i] */
    let rtol = ida_mem.ida_rtol;
    for (t, a) in tempv.data.iter_mut().zip(&ida_mem.ida_Vatol.data) {
        *t = rtol * *t + *a;
    }
    if ida_mem.ida_atolmin0 {
        if N_VMin(&tempv) <= ZERO {
            return -1;
        }
    }
    /* N_VInv(tempv1, weight) */
    for (w, t) in weight.data.iter_mut().zip(&tempv.data) {
        *w = ONE / *t;
    }
    0
}

/*
 * IDAQuadEwtSet
 *
 */
fn IDAQuadEwtSet(ida_mem: &IDAMem, qcur: &NVector, weightQ: &mut NVector) -> i32 {
    let mut flag = 0;

    match ida_mem.ida_itolQ {
        IDA_SS => flag = IDAQuadEwtSetSS(ida_mem, qcur, weightQ),
        IDA_SV => flag = IDAQuadEwtSetSV(ida_mem, qcur, weightQ),
        _ => {}
    }

    flag
}

/*
 * IDAQuadEwtSetSS
 *
 * (C uses ida_ypQ as temporary storage; a local scratch vector is
 * behaviorally identical and keeps the &IDAMem signature.)
 */
fn IDAQuadEwtSetSS(ida_mem: &IDAMem, qcur: &NVector, weightQ: &mut NVector) -> i32 {
    let mut tempvQ = N_VClone(qcur);
    /* N_VAbs(qcur, tempvQ) */
    for (t, q) in tempvQ.data.iter_mut().zip(&qcur.data) {
        *t = SUNRabs(*q);
    }
    /* N_VScale(rtolQ, tempvQ, tempvQ) */
    let rtolQ = ida_mem.ida_rtolQ;
    for t in tempvQ.data.iter_mut() {
        *t = rtolQ * *t;
    }
    /* N_VAddConst(tempvQ, SatolQ, tempvQ) */
    let satolQ = ida_mem.ida_SatolQ;
    for t in tempvQ.data.iter_mut() {
        *t += satolQ;
    }
    if ida_mem.ida_atolQmin0 {
        if N_VMin(&tempvQ) <= ZERO {
            return -1;
        }
    }
    /* N_VInv(tempvQ, weightQ) */
    for (w, t) in weightQ.data.iter_mut().zip(&tempvQ.data) {
        *w = ONE / *t;
    }

    0
}

/*
 * IDAQuadEwtSetSV
 *
 * (C uses ida_ypQ as temporary storage; local scratch here.)
 */
fn IDAQuadEwtSetSV(ida_mem: &IDAMem, qcur: &NVector, weightQ: &mut NVector) -> i32 {
    let mut tempvQ = N_VClone(qcur);
    /* N_VAbs(qcur, tempvQ) */
    for (t, q) in tempvQ.data.iter_mut().zip(&qcur.data) {
        *t = SUNRabs(*q);
    }
    /* N_VLinearSum(rtolQ, tempvQ, ONE, VatolQ, tempvQ) */
    let rtolQ = ida_mem.ida_rtolQ;
    for (t, a) in tempvQ.data.iter_mut().zip(&ida_mem.ida_VatolQ.data) {
        *t = rtolQ * *t + *a;
    }
    if ida_mem.ida_atolQmin0 {
        if N_VMin(&tempvQ) <= ZERO {
            return -1;
        }
    }
    /* N_VInv(tempvQ, weightQ) */
    for (w, t) in weightQ.data.iter_mut().zip(&tempvQ.data) {
        *w = ONE / *t;
    }

    0
}

/*
 * IDASensEwtSet
 *
 * (Non-static in C: used in IC for sensitivities.  Takes &mut IDAMem
 * so IDASensEwtSetEE can invoke a user efun; callers detach weightS
 * and yScur — e.g. phiS[0] — from IDAMem before the call.)
 */
pub fn IDASensEwtSet(ida_mem: &mut IDAMem, yScur: &[NVector], weightS: &mut [NVector]) -> i32 {
    let mut flag = 0;

    match ida_mem.ida_itolS {
        IDA_EE => flag = IDASensEwtSetEE(ida_mem, yScur, weightS),
        IDA_SS => flag = IDASensEwtSetSS(ida_mem, yScur, weightS),
        IDA_SV => flag = IDASensEwtSetSV(ida_mem, yScur, weightS),
        _ => {}
    }

    flag
}

/*
 * IDASensEwtSetEE
 *
 * In this case, the error weight vector for the i-th sensitivity is set to
 *
 * ewtS_i = pbar_i * efun(pbar_i*yS_i)
 *
 * In other words, the scaled sensitivity pbar_i * yS_i has the same error
 * weight vector calculation as the solution vector.
 *
 * (C uses ida_tempv1 as scratch for the scaled sensitivity; a local
 * keeps the borrows disjoint.  The efun call reproduces the
 * ida_efun_dispatch logic: user efun with ida_user_data, or the
 * internal IDAEwtSet.)
 */
fn IDASensEwtSetEE(ida_mem: &mut IDAMem, yScur: &[NVector], weightS: &mut [NVector]) -> i32 {
    for is in 0..(ida_mem.ida_Ns as usize) {
        /* N_VScale(pbar[is], yScur[is], pyS) */
        let mut pyS = N_VClone(&yScur[is]);
        let pb = ida_mem.ida_pbar[is];
        for (p, y) in pyS.data.iter_mut().zip(&yScur[is].data) {
            *p = pb * *y;
        }

        let flag = if ida_mem.ida_user_efun {
            let efun = ida_mem.ida_efun.unwrap();
            efun(&pyS, &mut weightS[is], &mut ida_mem.ida_user_data)
        } else {
            IDAEwtSet(ida_mem, &pyS, &mut weightS[is])
        };
        if flag != 0 {
            return -1;
        }

        /* N_VScale(pbar[is], weightS[is], weightS[is]) */
        weightS[is].scale_inplace(pb);
    }

    0
}

/*
 * IDASensEwtSetSS
 *
 */
fn IDASensEwtSetSS(ida_mem: &mut IDAMem, yScur: &[NVector], weightS: &mut [NVector]) -> i32 {
    for is in 0..(ida_mem.ida_Ns as usize) {
        let mut tempv = N_VClone(&yScur[is]);
        /* N_VAbs(yScur[is], tempv1) */
        for (t, y) in tempv.data.iter_mut().zip(&yScur[is].data) {
            *t = SUNRabs(*y);
        }
        /* N_VScale(rtolS, tempv1, tempv1) */
        let rtolS = ida_mem.ida_rtolS;
        for t in tempv.data.iter_mut() {
            *t = rtolS * *t;
        }
        /* N_VAddConst(tempv1, SatolS[is], tempv1) */
        let satol = ida_mem.ida_SatolS[is];
        for t in tempv.data.iter_mut() {
            *t += satol;
        }
        if ida_mem.ida_atolSmin0[is] {
            if N_VMin(&tempv) <= ZERO {
                return -1;
            }
        }
        /* N_VInv(tempv1, weightS[is]) */
        for (w, t) in weightS[is].data.iter_mut().zip(&tempv.data) {
            *w = ONE / *t;
        }
    }
    0
}

/*
 * IDASensEwtSetSV
 *
 */
fn IDASensEwtSetSV(ida_mem: &mut IDAMem, yScur: &[NVector], weightS: &mut [NVector]) -> i32 {
    for is in 0..(ida_mem.ida_Ns as usize) {
        let mut tempv = N_VClone(&yScur[is]);
        /* N_VAbs(yScur[is], tempv1) */
        for (t, y) in tempv.data.iter_mut().zip(&yScur[is].data) {
            *t = SUNRabs(*y);
        }
        /* N_VLinearSum(rtolS, tempv1, ONE, VatolS[is], tempv1) */
        let rtolS = ida_mem.ida_rtolS;
        for (t, a) in tempv.data.iter_mut().zip(&ida_mem.ida_VatolS[is].data) {
            *t = rtolS * *t + *a;
        }
        if ida_mem.ida_atolSmin0[is] {
            if N_VMin(&tempv) <= ZERO {
                return -1;
            }
        }
        /* N_VInv(tempv1, weightS[is]) */
        for (w, t) in weightS[is].data.iter_mut().zip(&tempv.data) {
            *w = ONE / *t;
        }
    }

    0
}

/*
 * IDAQuadSensEwtSet
 *
 * (Non-static in C.  &IDAMem suffices: the EE variant dispatches to
 * IDAQuadEwtSet, never to a user efun.)
 */
pub fn IDAQuadSensEwtSet(ida_mem: &IDAMem, yQScur: &[NVector], weightQS: &mut [NVector]) -> i32 {
    let mut flag = 0;

    match ida_mem.ida_itolQS {
        IDA_EE => flag = IDAQuadSensEwtSetEE(ida_mem, yQScur, weightQS),
        IDA_SS => flag = IDAQuadSensEwtSetSS(ida_mem, yQScur, weightQS),
        IDA_SV => flag = IDAQuadSensEwtSetSV(ida_mem, yQScur, weightQS),
        _ => {}
    }

    flag
}

/*
 * IDAQuadSensEwtSetEE
 *
 * In this case, the error weight vector for the i-th quadrature sensitivity
 * is set to
 *
 * ewtQS_i = pbar_i * IDAQuadEwtSet(pbar_i*yQS_i)
 *
 * In other words, the scaled sensitivity pbar_i * yQS_i has the same error
 * weight vector calculation as the quadrature vector.
 *
 * (C uses tempvQS[0] as scratch; local here.)
 */
fn IDAQuadSensEwtSetEE(ida_mem: &IDAMem, yQScur: &[NVector], weightQS: &mut [NVector]) -> i32 {
    for is in 0..(ida_mem.ida_Ns as usize) {
        /* N_VScale(pbar[is], yQScur[is], pyS) */
        let mut pyS = N_VClone(&yQScur[is]);
        let pb = ida_mem.ida_pbar[is];
        for (p, y) in pyS.data.iter_mut().zip(&yQScur[is].data) {
            *p = pb * *y;
        }

        let flag = IDAQuadEwtSet(ida_mem, &pyS, &mut weightQS[is]);
        if flag != 0 {
            return -1;
        }

        /* N_VScale(pbar[is], weightQS[is], weightQS[is]) */
        weightQS[is].scale_inplace(pb);
    }

    0
}

fn IDAQuadSensEwtSetSS(ida_mem: &IDAMem, yQScur: &[NVector], weightQS: &mut [NVector]) -> i32 {
    /* (C uses ida_ypQ as temporary storage; local scratch here) */
    for is in 0..(ida_mem.ida_Ns as usize) {
        let mut tempvQ = N_VClone(&yQScur[is]);
        /* N_VAbs(yQScur[is], tempvQ) */
        for (t, y) in tempvQ.data.iter_mut().zip(&yQScur[is].data) {
            *t = SUNRabs(*y);
        }
        /* N_VScale(rtolQS, tempvQ, tempvQ) */
        let rtolQS = ida_mem.ida_rtolQS;
        for t in tempvQ.data.iter_mut() {
            *t = rtolQS * *t;
        }
        /* N_VAddConst(tempvQ, SatolQS[is], tempvQ) */
        let satol = ida_mem.ida_SatolQS[is];
        for t in tempvQ.data.iter_mut() {
            *t += satol;
        }
        if ida_mem.ida_atolQSmin0[is] {
            if N_VMin(&tempvQ) <= ZERO {
                return -1;
            }
        }
        /* N_VInv(tempvQ, weightQS[is]) */
        for (w, t) in weightQS[is].data.iter_mut().zip(&tempvQ.data) {
            *w = ONE / *t;
        }
    }

    0
}

fn IDAQuadSensEwtSetSV(ida_mem: &IDAMem, yQScur: &[NVector], weightQS: &mut [NVector]) -> i32 {
    /* (C uses ida_ypQ as temporary storage; local scratch here) */
    for is in 0..(ida_mem.ida_Ns as usize) {
        let mut tempvQ = N_VClone(&yQScur[is]);
        /* N_VAbs(yQScur[is], tempvQ) */
        for (t, y) in tempvQ.data.iter_mut().zip(&yQScur[is].data) {
            *t = SUNRabs(*y);
        }
        /* N_VLinearSum(rtolQS, tempvQ, ONE, VatolQS[is], tempvQ) */
        let rtolQS = ida_mem.ida_rtolQS;
        for (t, a) in tempvQ.data.iter_mut().zip(&ida_mem.ida_VatolQS[is].data) {
            *t = rtolQS * *t + *a;
        }
        if ida_mem.ida_atolQSmin0[is] {
            if N_VMin(&tempvQ) <= ZERO {
                return -1;
            }
        }
        /* N_VInv(tempvQ, weightQS[is]) */
        for (w, t) in weightQS[is].data.iter_mut().zip(&tempvQ.data) {
            *w = ONE / *t;
        }
    }

    0
}

/* (donor Part-2 pattern: imports needed from here on) */
use crate::sundials_errors::SUN_ERR_ARG_CORRUPT;

/*
 * -----------------------------------------------------------------
 * Stopping tests
 * -----------------------------------------------------------------
 */

/*
 * IDAStopTest1
 *
 * This routine tests for stop conditions before taking a step.
 * The tests depend on the value of itask.
 * The variable tretlast is the previously returned value of tret.
 *
 * The return values are:
 * CONTINUE_STEPS       if no stop conditions were found
 * IDA_SUCCESS          for a normal return to the user
 * IDA_TSTOP_RETURN     for a tstop-reached return to the user
 * IDA_ILL_INPUT        for an illegal-input return to the user
 *
 * In the tstop cases, this routine may adjust the stepsize hh to cause
 * the next step to reach tstop exactly.
 */
fn IDAStopTest1(ida_mem: &mut IDAMem, tout: f64, tret: &mut f64, yret: &mut NVector,
                ypret: &mut NVector, itask: i32) -> i32 {
    if ida_mem.ida_tstopset {
        /* Test for tn past tstop */
        if (ida_mem.ida_tn - ida_mem.ida_tstop) * ida_mem.ida_hh > ZERO {
            IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAStopTest1", file!(),
                            &ida_msg_g(MSG_BAD_TSTOP, &[ida_mem.ida_tstop, ida_mem.ida_tn]));
            return IDA_ILL_INPUT;
        }

        let troundoff = HUNDRED * ida_mem.ida_uround *
            (SUNRabs(ida_mem.ida_tn) + SUNRabs(ida_mem.ida_hh));

        /* Test for tn at tstop */
        if SUNRabs(ida_mem.ida_tn - ida_mem.ida_tstop) <= troundoff {
            /* Ensure tout >= tstop, otherwise check for tout return below */
            if (tout - ida_mem.ida_tstop) * ida_mem.ida_hh >= ZERO
                || SUNRabs(tout - ida_mem.ida_tstop) <= troundoff
            {
                let ier = IDAGetSolution(ida_mem, ida_mem.ida_tstop, yret, ypret);
                if ier != IDA_SUCCESS {
                    IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAStopTest1",
                                    file!(),
                                    &ida_msg_g(MSG_BAD_TSTOP,
                                               &[ida_mem.ida_tstop, ida_mem.ida_tn]));
                    return IDA_ILL_INPUT;
                }
                *tret = ida_mem.ida_tstop;
                ida_mem.ida_tretlast = ida_mem.ida_tstop;
                ida_mem.ida_tstopset = SUNFALSE;
                return IDA_TSTOP_RETURN;
            }
        }
        /* Test for tn approaching tstop */
        else if (ida_mem.ida_tn + ida_mem.ida_hh - ida_mem.ida_tstop) * ida_mem.ida_hh > ZERO {
            ida_mem.ida_hh = (ida_mem.ida_tstop - ida_mem.ida_tn)
                * (ONE - FOUR * ida_mem.ida_uround);
        }
    }

    match itask {
        IDA_NORMAL => {
            /* Test for tout = tretlast, and for tn past tout. */
            if tout == ida_mem.ida_tretlast {
                *tret = tout;
                ida_mem.ida_tretlast = tout;
                return IDA_SUCCESS;
            }
            if (ida_mem.ida_tn - tout) * ida_mem.ida_hh >= ZERO {
                let ier = IDAGetSolution(ida_mem, tout, yret, ypret);
                if ier != IDA_SUCCESS {
                    IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAStopTest1",
                                    file!(), &ida_msg_g(MSG_BAD_TOUT, &[tout]));
                    return IDA_ILL_INPUT;
                }
                *tret = tout;
                ida_mem.ida_tretlast = tout;
                return IDA_SUCCESS;
            }

            CONTINUE_STEPS
        }

        IDA_ONE_STEP => {
            /* Test for tn past tretlast. */
            if (ida_mem.ida_tn - ida_mem.ida_tretlast) * ida_mem.ida_hh > ZERO {
                let _ = IDAGetSolution(ida_mem, ida_mem.ida_tn, yret, ypret);
                *tret = ida_mem.ida_tn;
                ida_mem.ida_tretlast = ida_mem.ida_tn;
                return IDA_SUCCESS;
            }

            CONTINUE_STEPS
        }

        _ => IDA_ILL_INPUT, /* This return should never happen. */
    }
}

/*
 * IDAStopTest2
 *
 * This routine tests for stop conditions after taking a step.
 * The tests depend on the value of itask.
 *
 * The return values are:
 *  CONTINUE_STEPS     if no stop conditions were found
 *  IDA_SUCCESS        for a normal return to the user
 *  IDA_TSTOP_RETURN   for a tstop-reached return to the user
 *  IDA_ILL_INPUT      for an illegal-input return to the user
 *
 * In the two cases with tstop, this routine may reset the stepsize hh
 * to cause the next step to reach tstop exactly.
 *
 * In the two cases with ONE_STEP mode, no interpolation to tn is needed
 * because yret and ypret already contain the current y and y' values.
 *
 * Note: No test is made for an error return from IDAGetSolution here,
 * because the same test was made prior to the step.
 */
fn IDAStopTest2(ida_mem: &mut IDAMem, tout: f64, tret: &mut f64, yret: &mut NVector,
                ypret: &mut NVector, itask: i32) -> i32 {
    if ida_mem.ida_tstopset {
        let troundoff = HUNDRED * ida_mem.ida_uround *
            (SUNRabs(ida_mem.ida_tn) + SUNRabs(ida_mem.ida_hh));

        /* Test for tn at tstop */
        if SUNRabs(ida_mem.ida_tn - ida_mem.ida_tstop) <= troundoff {
            /* Ensure tout >= tstop, otherwise check for tout return below */
            if (tout - ida_mem.ida_tstop) * ida_mem.ida_hh >= ZERO
                || SUNRabs(tout - ida_mem.ida_tstop) <= troundoff
            {
                /* ier = */
                let _ = IDAGetSolution(ida_mem, ida_mem.ida_tstop, yret, ypret);
                *tret = ida_mem.ida_tstop;
                ida_mem.ida_tretlast = ida_mem.ida_tstop;
                ida_mem.ida_tstopset = SUNFALSE;
                return IDA_TSTOP_RETURN;
            }
        }
        /* Test for tn approaching tstop */
        else if (ida_mem.ida_tn + ida_mem.ida_hh - ida_mem.ida_tstop) * ida_mem.ida_hh > ZERO {
            ida_mem.ida_hh = (ida_mem.ida_tstop - ida_mem.ida_tn)
                * (ONE - FOUR * ida_mem.ida_uround);
        }
    }

    match itask {
        IDA_NORMAL => {
            /* Test for tn past tout. */
            if (ida_mem.ida_tn - tout) * ida_mem.ida_hh >= ZERO {
                /* ier = */
                let _ = IDAGetSolution(ida_mem, tout, yret, ypret);
                *tret = tout;
                ida_mem.ida_tretlast = tout;
                return IDA_SUCCESS;
            }

            CONTINUE_STEPS
        }

        IDA_ONE_STEP => {
            *tret = ida_mem.ida_tn;
            ida_mem.ida_tretlast = ida_mem.ida_tn;
            /* (C relies on yret == ida_yy / ypret == ida_yp holding
               the current solution; copy the owned vectors back) */
            yret.data.copy_from_slice(&ida_mem.ida_yy.data);
            ypret.data.copy_from_slice(&ida_mem.ida_yp.data);
            IDA_SUCCESS
        }

        _ => IDA_ILL_INPUT, /* This return should never happen. */
    }
}

/*
 * -----------------------------------------------------------------
 * Error handler
 * -----------------------------------------------------------------
 */

/*
 * IDAHandleFailure
 *
 * This routine prints error messages for all cases of failure by
 * IDAStep.  It returns to IDASolve the value that it is to return to
 * the user.
 *
 * (Messages whose format contains no %.15g placeholder ignore the
 * extra values passed to ida_msg_g, exactly as C printf ignores
 * surplus varargs — the C argument lists are mirrored verbatim.)
 */
fn IDAHandleFailure(ida_mem: &mut IDAMem, sflag: i32) -> i32 {
    /* Depending on sflag, print error message and return error flag */
    match sflag {
        IDA_ERR_FAIL => {
            IDAProcessError(Some(ida_mem), IDA_ERR_FAIL, line!(), "IDAHandleFailure", file!(),
                            &ida_msg_g(MSG_ERR_FAILS, &[ida_mem.ida_tn, ida_mem.ida_hh]));
            IDA_ERR_FAIL
        }

        IDA_CONV_FAIL => {
            IDAProcessError(Some(ida_mem), IDA_CONV_FAIL, line!(), "IDAHandleFailure", file!(),
                            &ida_msg_g(MSG_CONV_FAILS, &[ida_mem.ida_tn, ida_mem.ida_hh]));
            IDA_CONV_FAIL
        }

        IDA_LSETUP_FAIL => {
            IDAProcessError(Some(ida_mem), IDA_LSETUP_FAIL, line!(), "IDAHandleFailure", file!(),
                            &ida_msg_g(MSG_SETUP_FAILED, &[ida_mem.ida_tn]));
            IDA_LSETUP_FAIL
        }

        IDA_LSOLVE_FAIL => {
            IDAProcessError(Some(ida_mem), IDA_LSOLVE_FAIL, line!(), "IDAHandleFailure", file!(),
                            &ida_msg_g(MSG_SOLVE_FAILED, &[ida_mem.ida_tn]));
            IDA_LSOLVE_FAIL
        }

        IDA_REP_RES_ERR => {
            IDAProcessError(Some(ida_mem), IDA_REP_RES_ERR, line!(), "IDAHandleFailure", file!(),
                            &ida_msg_g(MSG_REP_RES_ERR, &[ida_mem.ida_tn]));
            IDA_REP_RES_ERR
        }

        IDA_RES_FAIL => {
            IDAProcessError(Some(ida_mem), IDA_RES_FAIL, line!(), "IDAHandleFailure", file!(),
                            &ida_msg_g(MSG_RES_NONRECOV, &[ida_mem.ida_tn]));
            IDA_RES_FAIL
        }

        IDA_REP_QRHS_ERR => {
            IDAProcessError(Some(ida_mem), IDA_REP_QRHS_ERR, line!(), "IDAHandleFailure",
                            file!(), &ida_msg_g(MSG_QRHSFUNC_REPTD, &[ida_mem.ida_tn]));
            IDA_REP_QRHS_ERR
        }

        IDA_QRHS_FAIL => {
            IDAProcessError(Some(ida_mem), IDA_QRHS_FAIL, line!(), "IDAHandleFailure", file!(),
                            &ida_msg_g(MSG_QRHSFUNC_FAILED, &[ida_mem.ida_tn]));
            IDA_QRHS_FAIL
        }

        IDA_REP_SRES_ERR => {
            IDAProcessError(Some(ida_mem), IDA_REP_SRES_ERR, line!(), "IDAHandleFailure",
                            file!(), &ida_msg_g(MSG_SRES_REPTD, &[ida_mem.ida_tn]));
            IDA_REP_SRES_ERR
        }

        IDA_SRES_FAIL => {
            IDAProcessError(Some(ida_mem), IDA_SRES_FAIL, line!(), "IDAHandleFailure", file!(),
                            &ida_msg_g(MSG_SRES_FAILED, &[ida_mem.ida_tn]));
            IDA_SRES_FAIL
        }

        IDA_REP_QSRHS_ERR => {
            IDAProcessError(Some(ida_mem), IDA_REP_QSRHS_ERR, line!(), "IDAHandleFailure",
                            file!(), &ida_msg_g(MSG_QSRHSFUNC_REPTD, &[ida_mem.ida_tn]));
            IDA_REP_QSRHS_ERR
        }

        IDA_QSRHS_FAIL => {
            IDAProcessError(Some(ida_mem), IDA_QSRHS_FAIL, line!(), "IDAHandleFailure", file!(),
                            &ida_msg_g(MSG_QSRHSFUNC_FAILED, &[ida_mem.ida_tn]));
            IDA_QSRHS_FAIL
        }

        IDA_CONSTR_FAIL => {
            IDAProcessError(Some(ida_mem), IDA_CONSTR_FAIL, line!(), "IDAHandleFailure", file!(),
                            &ida_msg_g(MSG_FAILED_CONSTR, &[ida_mem.ida_tn]));
            IDA_CONSTR_FAIL
        }

        IDA_MEM_NULL => {
            IDAProcessError(None, IDA_MEM_NULL, line!(), "IDAHandleFailure", file!(), MSG_NO_MEM);
            IDA_MEM_NULL
        }

        SUN_ERR_ARG_CORRUPT => {
            IDAProcessError(Some(ida_mem), IDA_MEM_NULL, line!(), "IDAHandleFailure", file!(),
                            &ida_msg_g(MSG_NLS_INPUT_NULL, &[ida_mem.ida_tn]));
            IDA_MEM_NULL
        }

        IDA_NLS_SETUP_FAIL => {
            IDAProcessError(Some(ida_mem), IDA_NLS_SETUP_FAIL, line!(), "IDAHandleFailure",
                            file!(), &ida_msg_g(MSG_NLS_SETUP_FAILED, &[ida_mem.ida_tn]));
            IDA_NLS_SETUP_FAIL
        }

        IDA_NLS_FAIL => {
            IDAProcessError(Some(ida_mem), IDA_NLS_FAIL, line!(), "IDAHandleFailure", file!(),
                            &ida_msg_g(MSG_NLS_FAIL, &[ida_mem.ida_tn]));
            IDA_NLS_FAIL
        }

        _ => {
            /* This return should never happen */
            IDAProcessError(Some(ida_mem), IDA_UNRECOGNIZED_ERROR, line!(), "IDAHandleFailure",
                            file!(),
                            "IDA encountered an unrecognized error. Please report this \
                             to the Sundials developers at sundials-users@llnl.gov");
            IDA_UNRECOGNIZED_ERROR
        }
    }
}

/* (donor Part-2 pattern: further imports needed from here on;
   IDANls/IDASensNls live in idas_nls.rs / idas_nls_stg.rs with the
   collapsed SUNNonlinSolSolve_Newton loops — IDAStep dispatches) */
use crate::idas_nls::IDANls;
use crate::idas_nls_stg::IDASensNls;
use crate::sundials_math::{SUNMAX, SUNMIN, SUNRdifferentsign, SUNRpowerR, SUNRsqrt};

/*
 * -----------------------------------------------------------------
 * Main IDAStep function
 * -----------------------------------------------------------------
 */

/*
 * IDAStep
 *
 * This routine performs one internal IDA step, from tn to tn + hh.
 * It calls other routines to do all the work.
 *
 * It solves a system of differential/algebraic equations of the form
 *       F(t,y,y') = 0, for one step. In IDA, tt is used for t,
 * yy is used for y, and yp is used for y'. The function F is supplied as 'res'
 * by the user.
 *
 * The methods used are modified divided difference, fixed leading
 * coefficient forms of backward differentiation formulas.
 * The code adjusts the stepsize and order to control the local error per step.
 *
 * The main operations done here are as follows:
 *  * initialize various quantities;
 *  * setting of multistep method coefficients;
 *  * solution of the nonlinear system for yy at t = tn + hh;
 *  * deciding on order reduction and testing the local error;
 *  * attempting to recover from failure in nonlinear solver or error test;
 *  * resetting stepsize and order for the next step.
 *  * updating phi and other state data if successful;
 *
 * On a failure in the nonlinear system solution or error test, the
 * step may be reattempted, depending on the nature of the failure.
 *
 *       Return values are:
 *       IDA_SUCCESS   IDA_RES_FAIL      LSETUP_ERROR_NONRECVR
 *                     IDA_LSOLVE_FAIL   IDA_ERR_FAIL
 *                     IDA_CONSTR_FAIL   IDA_CONV_FAIL
 *                     IDA_REP_RES_ERR
 *
 * (C passes &IDA_mem->ida_ncfn / &IDA_mem->ida_netf (or the Q
 * variants) as long* to IDAHandleNFlag alongside IDA_mem; the
 * scalars are taken out and restored around each call to keep the
 * borrows disjoint.  SUNLogInfo/SUNLogInfoIf instrumentation is
 * dropped per convention.)
 */
fn IDAStep(ida_mem: &mut IDAMem) -> i32 {
    /* Are we computing sensitivities with the staggered or simultaneous approach? */
    let sensi_stg = ida_mem.ida_sensi && ida_mem.ida_ism == IDA_STAGGERED;
    let sensi_sim = ida_mem.ida_sensi && ida_mem.ida_ism == IDA_SIMULTANEOUS;

    let saved_t = ida_mem.ida_tn;

    /* Initialize failure counters for this step attempt */

    let mut ncf: i32 = 0; /* corrector failures  */
    let mut nef: i32 = 0; /* error test failures */
    let mut step_constraint_fails: i32 = 0;

    if ida_mem.ida_nst == 0 {
        ida_mem.ida_kk = 1;
        ida_mem.ida_kused = 0;
        ida_mem.ida_hused = ZERO;
        ida_mem.ida_psi[0] = ida_mem.ida_hh;
        ida_mem.ida_cj = ONE / ida_mem.ida_hh;
        ida_mem.ida_phase = 0;
        ida_mem.ida_ns = 0;
    }

    /* To prevent 'uninitialized variable' warnings */
    let mut err_k = ZERO;
    let mut err_km1 = ZERO;
    let mut err_km2 = ZERO;

    let mut ck = ZERO;

    /* Looping point for attempts to take a step */

    loop {
        /*-----------------------
          Set method coefficients
          -----------------------*/

        IDASetCoeffs(ida_mem, &mut ck);

        /*----------------------------------------------------
          If tn is past tstop (by roundoff), reset it to tstop.
          -----------------------------------------------------*/

        ida_mem.ida_tn += ida_mem.ida_hh;
        if ida_mem.ida_tstopset {
            if (ida_mem.ida_tn - ida_mem.ida_tstop) * ida_mem.ida_hh > ZERO {
                ida_mem.ida_tn = ida_mem.ida_tstop;
            }
        }

        /*-----------------------
          Advance state variables
          -----------------------*/

        /* Compute predicted values for yy and yp */
        IDAPredict(ida_mem);

        /* Compute predicted values for yyS and ypS (if simultaneous approach) */
        if sensi_sim {
            let mut yySp = std::mem::take(&mut ida_mem.ida_yySpredict);
            let mut ypSp = std::mem::take(&mut ida_mem.ida_ypSpredict);
            IDASensPredict(ida_mem, &mut yySp, &mut ypSp);
            ida_mem.ida_yySpredict = yySp;
            ida_mem.ida_ypSpredict = ypSp;
        }

        /* Nonlinear system solution */
        let mut nflag = IDANls(ida_mem);

        /* Nonlinear solve was successful */
        if nflag == IDA_SUCCESS {
            /* Check and enforce inequality constraints */
            if ida_mem.ida_constraintsSet {
                nflag = IDACheckConstraints(ida_mem, saved_t, &mut step_constraint_fails);

                /* Constraint check failed, predict again */
                if nflag == PREDICT_AGAIN {
                    continue;
                }

                /* Exit on nonrecoverable failure */
                if nflag != IDA_SUCCESS {
                    return nflag;
                }
            }

            /* Perform error test */
            nflag = IDATestError(ida_mem, ck, &mut err_k, &mut err_km1, &mut err_km2);
        }

        /* Test for convergence or error test failures */
        if nflag != IDA_SUCCESS {
            /* restore and decide what to do */
            IDARestore(ida_mem, saved_t);
            let mut ncfn = ida_mem.ida_ncfn;
            let mut netf = ida_mem.ida_netf;
            let kflag = IDAHandleNFlag(ida_mem, nflag, err_k, err_km1, &mut ncfn, &mut ncf,
                                       &mut netf, &mut nef);
            ida_mem.ida_ncfn = ncfn;
            ida_mem.ida_netf = netf;

            /* exit on nonrecoverable failure */
            if kflag != PREDICT_AGAIN {
                return kflag;
            }

            /* recoverable error; predict again */
            if ida_mem.ida_nst == 0 {
                IDAReset(ida_mem);
            }
            continue;
        }

        /*----------------------------
          Advance quadrature variables
          ----------------------------*/
        if ida_mem.ida_quadr {
            nflag = IDAQuadNls(ida_mem);

            /* If NLS was successful, perform error test */
            if ida_mem.ida_errconQ && nflag == IDA_SUCCESS {
                nflag = IDAQuadTestError(ida_mem, ck, &mut err_k, &mut err_km1, &mut err_km2);
            }

            /* Test for convergence or error test failures */
            if nflag != IDA_SUCCESS {
                /* restore and decide what to do */
                IDARestore(ida_mem, saved_t);
                let mut ncfnQ = ida_mem.ida_ncfnQ;
                let mut netfQ = ida_mem.ida_netfQ;
                let kflag = IDAHandleNFlag(ida_mem, nflag, err_k, err_km1, &mut ncfnQ, &mut ncf,
                                           &mut netfQ, &mut nef);
                ida_mem.ida_ncfnQ = ncfnQ;
                ida_mem.ida_netfQ = netfQ;

                /* exit on nonrecoverable failure */
                if kflag != PREDICT_AGAIN {
                    return kflag;
                }

                /* recoverable error; predict again */
                if ida_mem.ida_nst == 0 {
                    IDAReset(ida_mem);
                }
                continue;
            }
        }

        /*--------------------------------------------------
          Advance sensitivity variables (Staggered approach)
          --------------------------------------------------*/
        if sensi_stg {
            /* Evaluate res at converged y, needed for future evaluations of sens. RHS
               If res() fails recoverably, treat it as a convergence failure and
               attempt the step again */

            let retval = {
                let res = ida_mem.ida_res.unwrap();
                let tn = ida_mem.ida_tn;
                let IDAMem { ida_yy, ida_yp, ida_delta, ida_user_data, .. } = ida_mem;
                res(tn, ida_yy, ida_yp, ida_delta, ida_user_data)
            };

            if retval < 0 {
                return IDA_RES_FAIL;
            }
            if retval > 0 {
                continue;
            }

            /* Compute predicted values for yyS and ypS */
            let mut yySp = std::mem::take(&mut ida_mem.ida_yySpredict);
            let mut ypSp = std::mem::take(&mut ida_mem.ida_ypSpredict);
            IDASensPredict(ida_mem, &mut yySp, &mut ypSp);
            ida_mem.ida_yySpredict = yySp;
            ida_mem.ida_ypSpredict = ypSp;

            /* Nonlinear system solution */
            nflag = IDASensNls(ida_mem);

            /* If NLS was successful, perform error test */
            if ida_mem.ida_errconS && nflag == IDA_SUCCESS {
                nflag = IDASensTestError(ida_mem, ck, &mut err_k, &mut err_km1, &mut err_km2);
            }

            /* Test for convergence or error test failures */
            if nflag != IDA_SUCCESS {
                /* restore and decide what to do */
                IDARestore(ida_mem, saved_t);
                /* (C-exact quirk, preserved: the ncfnQ/netfQ counters
                   are passed here, not ncfnS/netfS) */
                let mut ncfnQ = ida_mem.ida_ncfnQ;
                let mut netfQ = ida_mem.ida_netfQ;
                let kflag = IDAHandleNFlag(ida_mem, nflag, err_k, err_km1, &mut ncfnQ, &mut ncf,
                                           &mut netfQ, &mut nef);
                ida_mem.ida_ncfnQ = ncfnQ;
                ida_mem.ida_netfQ = netfQ;

                /* exit on nonrecoverable failure */
                if kflag != PREDICT_AGAIN {
                    return kflag;
                }

                /* recoverable error; predict again */
                if ida_mem.ida_nst == 0 {
                    IDAReset(ida_mem);
                }
                continue;
            }
        }

        /*-------------------------------------------
          Advance quadrature sensitivity variables
          -------------------------------------------*/
        if ida_mem.ida_quadr_sensi {
            nflag = IDAQuadSensNls(ida_mem);

            /* If NLS was successful, perform error test */
            if ida_mem.ida_errconQS && nflag == IDA_SUCCESS {
                nflag = IDAQuadSensTestError(ida_mem, ck, &mut err_k, &mut err_km1, &mut err_km2);
            }

            /* Test for convergence or error test failures */
            if nflag != IDA_SUCCESS {
                /* restore and decide what to do */
                IDARestore(ida_mem, saved_t);
                let mut ncfnQ = ida_mem.ida_ncfnQ;
                let mut netfQ = ida_mem.ida_netfQ;
                let kflag = IDAHandleNFlag(ida_mem, nflag, err_k, err_km1, &mut ncfnQ, &mut ncf,
                                           &mut netfQ, &mut nef);
                ida_mem.ida_ncfnQ = ncfnQ;
                ida_mem.ida_netfQ = netfQ;

                /* exit on nonrecoverable failure */
                if kflag != PREDICT_AGAIN {
                    return kflag;
                }

                /* recoverable error; predict again */
                if ida_mem.ida_nst == 0 {
                    IDAReset(ida_mem);
                }
                continue;
            }
        }

        /* kflag == IDA_SUCCESS */
        break;
    } /* end loop */

    /* Nonlinear system solve and error test were both successful;
       update data, and consider change of step and/or order */

    IDACompleteStep(ida_mem, err_k, err_km1);

    /*
       Rescale ee vector to be the estimated local error
       Notes:
         (1) altering the value of ee is permissible since
             it will be overwritten by
             IDASolve()->IDAStep()->IDANls()
             before it is needed again
         (2) the value of ee is only valid if IDAHandleNFlag()
             returns either PREDICT_AGAIN or IDA_SUCCESS
    */

    /* N_VScale(ck, ee, ee) */
    ida_mem.ida_ee.scale_inplace(ck);

    IDA_SUCCESS
}

/*
 * IDASetCoeffs
 *
 *  This routine computes the coefficients relevant to the current step.
 *  The counter ns counts the number of consecutive steps taken at
 *  constant stepsize h and order k, up to a maximum of k + 2.
 *  Then the first ns components of beta will be one, and on a step
 *  with ns = k + 2, the coefficients alpha, etc. need not be reset here.
 *  Also, IDACompleteStep prohibits an order increase until ns = k + 2.
 *
 * (C stages beta[i] into ida_cvals — and gathers phiS/phiQS columns
 * into ida_Xvecs — for fused N_VScaleVectorArray calls; the staging
 * is dropped and the serial per-vector scale kernels are applied
 * inline.)
 */
fn IDASetCoeffs(ida_mem: &mut IDAMem, ck: &mut f64) {
    /* Set coefficients for the current stepsize h */

    if ida_mem.ida_hh != ida_mem.ida_hused || ida_mem.ida_kk != ida_mem.ida_kused {
        ida_mem.ida_ns = 0;
    }
    /* ns = SUNMIN(ns + 1, kused + 2) — integer min written out */
    ida_mem.ida_ns = if ida_mem.ida_ns + 1 < ida_mem.ida_kused + 2 {
        ida_mem.ida_ns + 1
    } else {
        ida_mem.ida_kused + 2
    };
    if ida_mem.ida_kk + 1 >= ida_mem.ida_ns {
        ida_mem.ida_beta[0] = ONE;
        ida_mem.ida_alpha[0] = ONE;
        let mut temp1 = ida_mem.ida_hh;
        ida_mem.ida_gamma[0] = ZERO;
        ida_mem.ida_sigma[0] = ONE;
        let mut i: i32 = 1;
        while i <= ida_mem.ida_kk {
            let iu = i as usize;
            let temp2 = ida_mem.ida_psi[iu - 1];
            ida_mem.ida_psi[iu - 1] = temp1;
            ida_mem.ida_beta[iu] = ida_mem.ida_beta[iu - 1] * ida_mem.ida_psi[iu - 1] / temp2;
            temp1 = temp2 + ida_mem.ida_hh;
            ida_mem.ida_alpha[iu] = ida_mem.ida_hh / temp1;
            ida_mem.ida_sigma[iu] = i as f64 * ida_mem.ida_sigma[iu - 1] * ida_mem.ida_alpha[iu];
            ida_mem.ida_gamma[iu] = ida_mem.ida_gamma[iu - 1]
                + ida_mem.ida_alpha[iu - 1] / ida_mem.ida_hh;
            i += 1;
        }
        ida_mem.ida_psi[ida_mem.ida_kk as usize] = temp1;
    }
    /* compute alphas, alpha0 */
    let mut alphas = ZERO;
    let mut alpha0 = ZERO;
    let mut i: i32 = 0;
    while i < ida_mem.ida_kk {
        alphas -= ONE / (i + 1) as f64;
        alpha0 -= ida_mem.ida_alpha[i as usize];
        i += 1;
    }

    /* compute leading coefficient cj  */
    ida_mem.ida_cjlast = ida_mem.ida_cj;
    ida_mem.ida_cj = -alphas / ida_mem.ida_hh;

    /* compute variable stepsize error coefficient ck */

    *ck = SUNRabs(ida_mem.ida_alpha[ida_mem.ida_kk as usize] + alphas - alpha0);
    *ck = SUNMAX(*ck, ida_mem.ida_alpha[ida_mem.ida_kk as usize]);

    /* change phi to phi-star  */

    /* Scale i=IDA_mem->ida_ns to i<=IDA_mem->ida_kk */
    if ida_mem.ida_ns <= ida_mem.ida_kk {
        let ns_i = ida_mem.ida_ns as usize;
        let kk = ida_mem.ida_kk as usize;

        /* N_VScaleVectorArray(kk - ns + 1, beta[ns..], phi + ns, phi + ns) */
        for i in ns_i..=kk {
            let b = ida_mem.ida_beta[i];
            ida_mem.ida_phi[i].scale_inplace(b);
        }

        if ida_mem.ida_quadr {
            for i in ns_i..=kk {
                let b = ida_mem.ida_beta[i];
                ida_mem.ida_phiQ[i].scale_inplace(b);
            }
        }

        /* (C's sensi||quadr_sensi cvals-staging loop vanishes) */

        if ida_mem.ida_sensi {
            for i in ns_i..=kk {
                let b = ida_mem.ida_beta[i];
                for is in 0..(ida_mem.ida_Ns as usize) {
                    ida_mem.ida_phiS[i][is].scale_inplace(b);
                }
            }
        }

        if ida_mem.ida_quadr_sensi {
            for i in ns_i..=kk {
                let b = ida_mem.ida_beta[i];
                for is in 0..(ida_mem.ida_Ns as usize) {
                    ida_mem.ida_phiQS[i][is].scale_inplace(b);
                }
            }
        }
    }
}

/*
 * -----------------------------------------------------------------
 * Nonlinear solver functions
 * -----------------------------------------------------------------
 */

/* (IDANls sits here in idas.c; its body — with the collapsed
   SUNNonlinSolSolve_Newton loop — lives in idas_nls.rs per the donor
   decision, and IDAStep dispatches to it.) */

fn IDACheckConstraints(ida_mem: &mut IDAMem, saved_t: f64,
                       step_constraint_fails: &mut i32) -> i32 {
    /* N_Vector mm  = tempv2 (mask), tmp = tempv1 (workspace) — used
       directly as fields below; serial kernels inline */

    /* Get mask vector mm, 1 where constraints failed and 0 otherwise */
    let constraints_passed = {
        let IDAMem { ida_constraints, ida_yy, ida_tempv2, .. } = ida_mem;
        N_VConstrMask(ida_constraints, ida_yy, ida_tempv2)
    };
    if constraints_passed {
        return IDA_SUCCESS;
    }

    /* Constraints not met */

    /* Compute correction v such that y - v will satisfy the constraints
     *
     * 1. Create a mask array that is +1 where constraints are strictly greater
     *    than or less than zero (|c[i]| = 2) and 0 otherwise
     *
     * 2. Create a mask array that is +/- 2 where constraints are strictly greater
     *    than (+) or less than (-) zero and 0 otherwise
     *
     * 3. Use error weights to compute an adjustment vector for values with strict
     *    constraints, a[i] = +/- 2 * w[i] = +/- 2 * (atol * |y[i]| + rtol[i]),
     *    and is 0 otherwise
     *
     * 4. Save the adjustment vector for possible use later
     *
     * 5. Compute correction vector for all values, v[i] = y[i] - 0.1 * a[i] for
     *    strict constraints and v[i] = y[i] otherwise
     *
     * 6. Zero out entries where the constraints passed, v = mask * v
     */
    {
        let IDAMem { ida_constraints, ida_ewt, ida_yy, ida_tempv1, ida_tempv2,
                     ida_tempv3, .. } = ida_mem;
        let tmp = ida_tempv1;
        let mm = ida_tempv2;

        /* N_VCompare(ONEPT5, constraints, tmp) */
        for (t, c) in tmp.data.iter_mut().zip(&ida_constraints.data) {
            *t = if SUNRabs(*c) >= ONEPT5 { ONE } else { ZERO };
        }
        /* N_VProd(tmp, constraints, tmp) */
        for (t, c) in tmp.data.iter_mut().zip(&ida_constraints.data) {
            *t = *t * *c;
        }
        /* N_VDiv(tmp, ewt, tmp) */
        for (t, w) in tmp.data.iter_mut().zip(&ida_ewt.data) {
            *t = *t / *w;
        }
        /* N_VScale(-PT1, tmp, tempv3) */
        for (v, t) in ida_tempv3.data.iter_mut().zip(&tmp.data) {
            *v = -PT1 * *t;
        }
        /* N_VLinearSum(ONE, yy, -PT1, tmp, tmp) — serial kernel
           VLin1(b, y, x, z): z[i] = b*y[i] + x[i] */
        for (t, y) in tmp.data.iter_mut().zip(&ida_yy.data) {
            *t = -PT1 * *t + *y;
        }
        /* N_VProd(tmp, mm, tmp) */
        for (t, m) in tmp.data.iter_mut().zip(&mm.data) {
            *t = *t * *m;
        }
    }

    let vnorm = IDAWrmsNorm(ida_mem, &ida_mem.ida_tempv1, &ida_mem.ida_ewt,
                            SUNFALSE); /* ||v|| */

    /* If constraint correction vector is small in norm (satisfies the nonlinear
       solver convergence condition with R = 1), correct and accept this step */
    if vnorm <= ida_mem.ida_epsNewt {
        /* Update constraint correction count */
        ida_mem.constraint_corrections += 1;

        /* To reduce roundoff errors that can violate the constraints, split the
         * correction update, ee = ee - v, into three steps */
        let IDAMem { ida_ee, ida_yypredict, ida_tempv1, ida_tempv2, ida_tempv3, .. } = ida_mem;
        let tmp = ida_tempv1;
        let mm = ida_tempv2;

        /* Zero out the correction where any constraint failed */
        /* N_VProd(mm, ee, tmp); N_VLinearSum(ONE, ee, -ONE, tmp, ee) */
        for (t, (m, e)) in tmp.data.iter_mut().zip(mm.data.iter().zip(&ida_ee.data)) {
            *t = *m * *e;
        }
        for (e, t) in ida_ee.data.iter_mut().zip(&tmp.data) {
            *e = *e - *t;
        }

        /* Set correction to zero out the predictor where any constraint failed */
        /* N_VProd(mm, yypredict, tmp); N_VLinearSum(ONE, ee, -ONE, tmp, ee) */
        for (t, (m, y)) in tmp.data.iter_mut().zip(mm.data.iter().zip(&ida_yypredict.data)) {
            *t = *m * *y;
        }
        for (e, t) in ida_ee.data.iter_mut().zip(&tmp.data) {
            *e = *e - *t;
        }

        /* Update the correction where constraints failed and are strictly greater
           or less than zero to shift the state with the adjustment saved above */
        /* N_VProd(mm, tempv3, tempv3); N_VLinearSum(ONE, ee, -ONE, tempv3, ee) */
        for (v, m) in ida_tempv3.data.iter_mut().zip(&mm.data) {
            *v = *m * *v;
        }
        for (e, v) in ida_ee.data.iter_mut().zip(&ida_tempv3.data) {
            *e = *e - *v;
        }

        return IDA_SUCCESS;
    }

    /* update failure counts */
    *step_constraint_fails += 1;
    ida_mem.constraint_fails += 1;

    /* Return with error if |h| == hmin */
    if SUNRabs(ida_mem.ida_hh) <= ida_mem.ida_hmin * ONEPSM {
        return IDA_CONSTR_FAIL;
    }

    /* Return with error if max step attempt failures */
    if *step_constraint_fails == ida_mem.max_constraint_fails {
        return IDA_CONSTR_FAIL;
    }

    /* Constraints correction is too large, reduce h by computing rr = h'/h */
    {
        let IDAMem { ida_phi, ida_yy, ida_tempv1, ida_tempv2, .. } = ida_mem;
        let tmp = ida_tempv1;
        let mm = ida_tempv2;
        /* N_VLinearSum(ONE, phi[0], -ONE, yy, tmp) — VDiff */
        for (t, (p, y)) in tmp.data.iter_mut().zip(ida_phi[0].data.iter().zip(&ida_yy.data)) {
            *t = *p - *y;
        }
        /* N_VProd(mm, tmp, tmp) */
        for (t, m) in tmp.data.iter_mut().zip(&mm.data) {
            *t = *m * *t;
        }
    }
    ida_mem.ida_eta = PT9 * N_VMinQuotient(&ida_mem.ida_phi[0], &ida_mem.ida_tempv1);
    ida_mem.ida_eta = SUNMAX(ida_mem.ida_eta, PT1);
    ida_mem.ida_eta = SUNMAX(ida_mem.ida_eta,
                             ida_mem.ida_hmin / SUNRabs(ida_mem.ida_hh));

    /* Reattempt step with new step size */
    IDARestore(ida_mem, saved_t);
    ida_mem.ida_phase = 1;
    ida_mem.ida_hh *= ida_mem.ida_eta;
    if ida_mem.ida_nst == 0 {
        IDAReset(ida_mem);
    }

    PREDICT_AGAIN
}

/*
 * IDAPredict
 *
 * This routine predicts the new values for vectors yy and yp.
 *
 * (C: cvals[j] = ONE and fused N_VLinearCombination(kk+1, cvals, phi,
 * yypredict) / (kk, gamma+1, phi+1, yppredict); the coefficient
 * staging is dropped and the serial kernels — z = c[0]*X[0], then
 * z += c[j]*X[j] — are applied inline.  Multiplication by ONE is
 * exact, so the unit-coefficient combination is a copy plus
 * accumulations.)
 */
fn IDAPredict(ida_mem: &mut IDAMem) {
    let kk = ida_mem.ida_kk as usize;
    let IDAMem { ida_phi, ida_gamma, ida_yypredict, ida_yppredict, .. } = ida_mem;

    /* yypredict = sum_{j=0..kk} ONE * phi[j] */
    ida_yypredict.data.copy_from_slice(&ida_phi[0].data);
    for j in 1..=kk {
        for (y, p) in ida_yypredict.data.iter_mut().zip(&ida_phi[j].data) {
            *y += *p;
        }
    }

    /* yppredict = sum_{j=1..kk} gamma[j] * phi[j] */
    for (y, p) in ida_yppredict.data.iter_mut().zip(&ida_phi[1].data) {
        *y = ida_gamma[1] * *p;
    }
    for j in 2..=kk {
        for (y, p) in ida_yppredict.data.iter_mut().zip(&ida_phi[j].data) {
            *y += ida_gamma[j] * *p;
        }
    }
}

/*
 * IDAQuadNls
 *
 * This routine solves for the quadrature variables at the new step.
 * It does not solve a nonlinear system, but rather updates the
 * quadrature variables. The name for this function is just for
 * uniformity purposes.
 *
 */
fn IDAQuadNls(ida_mem: &mut IDAMem) -> i32 {
    /* Predict: load yyQ and ypQ */
    IDAQuadPredict(ida_mem);

    /* Compute correction eeQ */
    let retval = {
        let rhsQ = ida_mem.ida_rhsQ.unwrap();
        let tn = ida_mem.ida_tn;
        let IDAMem { ida_yy, ida_yp, ida_eeQ, ida_user_data, .. } = ida_mem;
        rhsQ(tn, ida_yy, ida_yp, ida_eeQ, ida_user_data)
    };
    ida_mem.ida_nrQe += 1;
    if retval < 0 {
        return IDA_QRHS_FAIL;
    } else if retval > 0 {
        return IDA_QRHS_RECVR;
    }

    if ida_mem.ida_quadr_sensi {
        /* N_VScale(ONE, eeQ, savrhsQ) */
        let IDAMem { ida_eeQ, ida_savrhsQ, .. } = ida_mem;
        ida_savrhsQ.data.copy_from_slice(&ida_eeQ.data);
    }

    {
        let cj = ida_mem.ida_cj;
        let IDAMem { ida_eeQ, ida_ypQ, ida_yyQ, .. } = ida_mem;

        /* N_VLinearSum(ONE, eeQ, -ONE, ypQ, eeQ) — VDiff */
        for (e, yp) in ida_eeQ.data.iter_mut().zip(&ida_ypQ.data) {
            *e = *e - *yp;
        }
        /* N_VScale(ONE/cj, eeQ, eeQ) */
        let rcj = ONE / cj;
        for e in ida_eeQ.data.iter_mut() {
            *e = rcj * *e;
        }

        /* Apply correction: yyQ = yyQ + eeQ */
        /* N_VLinearSum(ONE, yyQ, ONE, eeQ, yyQ) — VSum */
        for (y, e) in ida_yyQ.data.iter_mut().zip(&ida_eeQ.data) {
            *y = *y + *e;
        }
    }

    IDA_SUCCESS
}

/*
 * IDAQuadPredict
 *
 * This routine predicts the new value for vectors yyQ and ypQ
 */
fn IDAQuadPredict(ida_mem: &mut IDAMem) {
    let kk = ida_mem.ida_kk as usize;
    let IDAMem { ida_phiQ, ida_gamma, ida_yyQ, ida_ypQ, .. } = ida_mem;

    /* yyQ = sum_{j=0..kk} ONE * phiQ[j] */
    ida_yyQ.data.copy_from_slice(&ida_phiQ[0].data);
    for j in 1..=kk {
        for (y, p) in ida_yyQ.data.iter_mut().zip(&ida_phiQ[j].data) {
            *y += *p;
        }
    }

    /* ypQ = sum_{j=1..kk} gamma[j] * phiQ[j] */
    for (y, p) in ida_ypQ.data.iter_mut().zip(&ida_phiQ[1].data) {
        *y = ida_gamma[1] * *p;
    }
    for j in 2..=kk {
        for (y, p) in ida_ypQ.data.iter_mut().zip(&ida_phiQ[j].data) {
            *y += ida_gamma[j] * *p;
        }
    }
}

/* (IDASensNls sits here in idas.c; its body lives in idas_nls_stg.rs
   per the donor decision, and IDAStep dispatches to it.) */

/*
 * IDASensPredict
 *
 * This routine loads the predicted values for the is-th sensitivity
 * in the vectors yySens and ypSens.
 *
 * When ism=IDA_STAGGERED,  yySens = yyS[is] and ypSens = ypS[is]
 *
 * (Both C call sites pass ida_yySpredict/ida_ypSpredict; callers
 * detach those rows around the call.  Fused
 * N_VLinearCombinationVectorArray → per-is serial kernels inline.)
 */
fn IDASensPredict(ida_mem: &IDAMem, yySens: &mut [NVector], ypSens: &mut [NVector]) {
    let kk = ida_mem.ida_kk as usize;

    for is in 0..(ida_mem.ida_Ns as usize) {
        /* yySens[is] = sum_{j=0..kk} ONE * phiS[j][is] */
        yySens[is].data.copy_from_slice(&ida_mem.ida_phiS[0][is].data);
        for j in 1..=kk {
            for (y, p) in yySens[is].data.iter_mut().zip(&ida_mem.ida_phiS[j][is].data) {
                *y += *p;
            }
        }

        /* ypSens[is] = sum_{j=1..kk} gamma[j] * phiS[j][is] */
        for (y, p) in ypSens[is].data.iter_mut().zip(&ida_mem.ida_phiS[1][is].data) {
            *y = ida_mem.ida_gamma[1] * *p;
        }
        for j in 2..=kk {
            for (y, p) in ypSens[is].data.iter_mut().zip(&ida_mem.ida_phiS[j][is].data) {
                *y += ida_mem.ida_gamma[j] * *p;
            }
        }
    }
}

/*
 * IDAQuadSensNls
 *
 * This routine solves for the snesitivity quadrature variables at the
 * new step. It does not solve a nonlinear system, but rather updates
 * the sensitivity variables. The name for this function is just for
 * uniformity purposes.
 *
 */
fn IDAQuadSensNls(ida_mem: &mut IDAMem) -> i32 {
    /* Predict: load yyQS and ypQS for each sensitivity. Store
       1st order information in tempvQS. */
    /* (C: ypQS = ida_tempvQS alias; both rows detached for the call
       and for the rhsQS dispatch below) */
    let mut yyQS = std::mem::take(&mut ida_mem.ida_yyQS);
    let mut ypQS = std::mem::take(&mut ida_mem.ida_tempvQS);
    IDAQuadSensPredict(ida_mem, &mut yyQS, &mut ypQS);

    /* Compute correction eeQS */
    let ns = ida_mem.ida_Ns;
    let tn = ida_mem.ida_tn;
    let yy = std::mem::take(&mut ida_mem.ida_yy);
    let yp = std::mem::take(&mut ida_mem.ida_yp);
    let yyS = std::mem::take(&mut ida_mem.ida_yyS);
    let ypS = std::mem::take(&mut ida_mem.ida_ypS);
    let savrhsQ = std::mem::take(&mut ida_mem.ida_savrhsQ);
    let mut eeQS = std::mem::take(&mut ida_mem.ida_eeQS);
    let mut tempv1 = std::mem::take(&mut ida_mem.ida_tempv1);
    let mut tempv2 = std::mem::take(&mut ida_mem.ida_tempv2);
    let mut tmpS3 = std::mem::take(&mut ida_mem.ida_tmpS3);

    let retval = if ida_mem.ida_rhsQSDQ {
        IDAQuadSensRhsInternalDQ(ida_mem, ns, tn, &yy, &yp, &yyS, &ypS, &savrhsQ, &mut eeQS,
                                 &mut tempv1, &mut tempv2, &mut tmpS3)
    } else {
        let rhsQS = ida_mem.ida_rhsQS.unwrap();
        rhsQS(ns, tn, &yy, &yp, &yyS, &ypS, &savrhsQ, &mut eeQS, &mut ida_mem.ida_user_data,
              &mut tempv1, &mut tempv2, &mut tmpS3)
    };

    ida_mem.ida_yy = yy;
    ida_mem.ida_yp = yp;
    ida_mem.ida_yyS = yyS;
    ida_mem.ida_ypS = ypS;
    ida_mem.ida_savrhsQ = savrhsQ;
    ida_mem.ida_tempv1 = tempv1;
    ida_mem.ida_tempv2 = tempv2;
    ida_mem.ida_tmpS3 = tmpS3;

    ida_mem.ida_nrQSe += 1;

    let mut ret = IDA_SUCCESS;
    if retval < 0 {
        ret = IDA_QSRHS_FAIL;
    } else if retval > 0 {
        ret = IDA_QSRHS_RECVR;
    }

    if ret == IDA_SUCCESS {
        /* retval = N_VLinearSumVectorArray(Ns, ONE/cj, eeQS, -ONE/cj,
           ypQS, eeQS) — per-is general kernel z = a*x + b*y; the
           cannot-fail VECTOROP_ERR branch vanishes */
        let a = ONE / ida_mem.ida_cj;
        let b = -ONE / ida_mem.ida_cj;
        for is in 0..(ns as usize) {
            for (e, yp) in eeQS[is].data.iter_mut().zip(&ypQS[is].data) {
                *e = a * *e + b * *yp;
            }
        }

        /* Apply correction: yyQS[is] = yyQ[is] + eeQ[is] */
        /* N_VLinearSumVectorArray(Ns, ONE, yyQS, ONE, eeQS, yyQS) — VSum */
        for is in 0..(ns as usize) {
            for (y, e) in yyQS[is].data.iter_mut().zip(&eeQS[is].data) {
                *y = *y + *e;
            }
        }
    }

    ida_mem.ida_yyQS = yyQS;
    ida_mem.ida_tempvQS = ypQS;
    ida_mem.ida_eeQS = eeQS;

    ret
}

/*
 * IDAQuadSensPredict
 *
 * This routine predicts the new value for vectors yyQS and ypQS
 */
fn IDAQuadSensPredict(ida_mem: &IDAMem, yQS: &mut [NVector], ypQS: &mut [NVector]) {
    let kk = ida_mem.ida_kk as usize;

    for is in 0..(ida_mem.ida_Ns as usize) {
        /* yQS[is] = sum_{j=0..kk} ONE * phiQS[j][is] */
        yQS[is].data.copy_from_slice(&ida_mem.ida_phiQS[0][is].data);
        for j in 1..=kk {
            for (y, p) in yQS[is].data.iter_mut().zip(&ida_mem.ida_phiQS[j][is].data) {
                *y += *p;
            }
        }

        /* ypQS[is] = sum_{j=1..kk} gamma[j] * phiQS[j][is] */
        for (y, p) in ypQS[is].data.iter_mut().zip(&ida_mem.ida_phiQS[1][is].data) {
            *y = ida_mem.ida_gamma[1] * *p;
        }
        for j in 2..=kk {
            for (y, p) in ypQS[is].data.iter_mut().zip(&ida_mem.ida_phiQS[j][is].data) {
                *y += ida_mem.ida_gamma[j] * *p;
            }
        }
    }
}

/*
 * -----------------------------------------------------------------
 * Error test
 * -----------------------------------------------------------------
 */

/*
 * IDATestError
 *
 * This routine estimates errors at orders k, k-1, k-2, decides
 * whether or not to suggest an order decrease, and performs
 * the local error test.
 *
 * IDATestError returns either IDA_SUCCESS or ERROR_TEST_FAIL.
 */
fn IDATestError(ida_mem: &mut IDAMem, ck: f64, err_k: &mut f64, err_km1: &mut f64,
                err_km2: &mut f64) -> i32 {
    /* Compute error for order k. */
    let enorm_k = IDAWrmsNorm(ida_mem, &ida_mem.ida_ee, &ida_mem.ida_ewt,
                              ida_mem.ida_suppressalg);
    *err_k = ida_mem.ida_sigma[ida_mem.ida_kk as usize] * enorm_k;
    let terr_k = (ida_mem.ida_kk + 1) as f64 * *err_k;

    ida_mem.ida_knew = ida_mem.ida_kk;

    if ida_mem.ida_kk > 1 {
        /* Compute error at order k-1 */
        /* N_VLinearSum(ONE, phi[kk], ONE, ee, delta) — VSum */
        {
            let kk = ida_mem.ida_kk as usize;
            let IDAMem { ida_phi, ida_ee, ida_delta, .. } = ida_mem;
            for (d, (p, e)) in ida_delta.data.iter_mut()
                .zip(ida_phi[kk].data.iter().zip(&ida_ee.data))
            {
                *d = *p + *e;
            }
        }
        let enorm_km1 = IDAWrmsNorm(ida_mem, &ida_mem.ida_delta, &ida_mem.ida_ewt,
                                    ida_mem.ida_suppressalg);
        *err_km1 = ida_mem.ida_sigma[ida_mem.ida_kk as usize - 1] * enorm_km1;
        let terr_km1 = ida_mem.ida_kk as f64 * *err_km1;

        if ida_mem.ida_kk > 2 {
            /* Compute error at order k-2 */
            /* N_VLinearSum(ONE, phi[kk-1], ONE, delta, delta) — VSum,
               y aliases z: d = p + d */
            {
                let kk = ida_mem.ida_kk as usize;
                let IDAMem { ida_phi, ida_delta, .. } = ida_mem;
                for (d, p) in ida_delta.data.iter_mut().zip(&ida_phi[kk - 1].data) {
                    *d = *p + *d;
                }
            }
            let enorm_km2 = IDAWrmsNorm(ida_mem, &ida_mem.ida_delta, &ida_mem.ida_ewt,
                                        ida_mem.ida_suppressalg);
            *err_km2 = ida_mem.ida_sigma[ida_mem.ida_kk as usize - 2] * enorm_km2;
            let terr_km2 = (ida_mem.ida_kk - 1) as f64 * *err_km2;

            /* Decrease order if errors are reduced */
            if SUNMAX(terr_km1, terr_km2) <= terr_k {
                ida_mem.ida_knew = ida_mem.ida_kk - 1;
            }
        } else {
            /* Decrease order to 1 if errors are reduced by at least 1/2 */
            if terr_km1 <= HALF * terr_k {
                ida_mem.ida_knew = ida_mem.ida_kk - 1;
            }
        }
    }

    /* Perform error test */
    if ck * enorm_k > ONE {
        ERROR_TEST_FAIL
    } else {
        IDA_SUCCESS
    }
}

/*
 * IDAQuadTestError
 *
 * This routine estimates quadrature errors and updates errors at
 * orders k, k-1, k-2, decides whether or not to suggest an order reduction,
 * and performs the local error test.
 *
 * IDAQuadTestError returns the updated local error estimate at orders k,
 * k-1, and k-2. These are norms of type SUNMAX(|err|,|errQ|).
 *
 * The return flag can be either IDA_SUCCESS or ERROR_TEST_FAIL.
 *
 * (C renames ypQ as tempv; ida_ypQ is used directly as the scratch.)
 */
fn IDAQuadTestError(ida_mem: &mut IDAMem, ck: f64, err_k: &mut f64, err_km1: &mut f64,
                    err_km2: &mut f64) -> i32 {
    let mut check_for_reduction = SUNFALSE;

    /* Update error for order k. */
    let enormQ = N_VWrmsNorm(&ida_mem.ida_eeQ, &ida_mem.ida_ewtQ);
    let errQ_k = ida_mem.ida_sigma[ida_mem.ida_kk as usize] * enormQ;
    if errQ_k > *err_k {
        *err_k = errQ_k;
        check_for_reduction = SUNTRUE;
    }
    let terr_k = (ida_mem.ida_kk + 1) as f64 * *err_k;

    if ida_mem.ida_kk > 1 {
        /* Update error at order k-1 */
        /* N_VLinearSum(ONE, phiQ[kk], ONE, eeQ, tempv) — VSum */
        {
            let kk = ida_mem.ida_kk as usize;
            let IDAMem { ida_phiQ, ida_eeQ, ida_ypQ, .. } = ida_mem;
            for (t, (p, e)) in ida_ypQ.data.iter_mut()
                .zip(ida_phiQ[kk].data.iter().zip(&ida_eeQ.data))
            {
                *t = *p + *e;
            }
        }
        let errQ_km1 = ida_mem.ida_sigma[ida_mem.ida_kk as usize - 1]
            * N_VWrmsNorm(&ida_mem.ida_ypQ, &ida_mem.ida_ewtQ);
        if errQ_km1 > *err_km1 {
            *err_km1 = errQ_km1;
            check_for_reduction = SUNTRUE;
        }
        let terr_km1 = ida_mem.ida_kk as f64 * *err_km1;

        /* Has an order decrease already been decided in IDATestError? */
        if ida_mem.ida_knew != ida_mem.ida_kk {
            check_for_reduction = SUNFALSE;
        }

        if check_for_reduction {
            if ida_mem.ida_kk > 2 {
                /* Update error at order k-2 */
                /* N_VLinearSum(ONE, phiQ[kk-1], ONE, tempv, tempv) —
                   VSum, y aliases z */
                {
                    let kk = ida_mem.ida_kk as usize;
                    let IDAMem { ida_phiQ, ida_ypQ, .. } = ida_mem;
                    for (t, p) in ida_ypQ.data.iter_mut().zip(&ida_phiQ[kk - 1].data) {
                        *t = *p + *t;
                    }
                }
                let errQ_km2 = ida_mem.ida_sigma[ida_mem.ida_kk as usize - 2]
                    * N_VWrmsNorm(&ida_mem.ida_ypQ, &ida_mem.ida_ewtQ);
                if errQ_km2 > *err_km2 {
                    *err_km2 = errQ_km2;
                }
                let terr_km2 = (ida_mem.ida_kk - 1) as f64 * *err_km2;

                /* Decrease order if errors are reduced */
                if SUNMAX(terr_km1, terr_km2) <= terr_k {
                    ida_mem.ida_knew = ida_mem.ida_kk - 1;
                }
            } else {
                /* Decrease order to 1 if errors are reduced by at least 1/2 */
                if terr_km1 <= HALF * terr_k {
                    ida_mem.ida_knew = ida_mem.ida_kk - 1;
                }
            }
        }
    }

    /* Perform error test */
    if ck * enormQ > ONE {
        ERROR_TEST_FAIL
    } else {
        IDA_SUCCESS
    }
}

/*
 * IDASensTestError
 *
 * This routine estimates sensitivity errors and updates errors at
 * orders k, k-1, k-2, decides whether or not to suggest an order reduction,
 * and performs the local error test. (Used only in staggered approach).
 *
 * IDASensTestError returns the updated local error estimate at orders k,
 * k-1, and k-2. These are norms of type SUNMAX(|err|,|errQ|,|errS|).
 *
 * The return flag can be either IDA_SUCCESS or ERROR_TEST_FAIL.
 *
 * (C renames deltaS as tempv; ida_deltaS is used directly.  The
 * fused N_VLinearSumVectorArray calls become per-is VSum kernels and
 * their cannot-fail VECTOROP_ERR branches vanish.)
 */
fn IDASensTestError(ida_mem: &mut IDAMem, ck: f64, err_k: &mut f64, err_km1: &mut f64,
                    err_km2: &mut f64) -> i32 {
    let mut check_for_reduction = SUNFALSE;

    /* Update error for order k. */
    let enormS = IDASensWrmsNorm(ida_mem, &ida_mem.ida_eeS, &ida_mem.ida_ewtS,
                                 ida_mem.ida_suppressalg);
    let errS_k = ida_mem.ida_sigma[ida_mem.ida_kk as usize] * enormS;
    if errS_k > *err_k {
        *err_k = errS_k;
        check_for_reduction = SUNTRUE;
    }
    let terr_k = (ida_mem.ida_kk + 1) as f64 * *err_k;

    if ida_mem.ida_kk > 1 {
        /* Update error at order k-1 */
        {
            let kk = ida_mem.ida_kk as usize;
            let ns = ida_mem.ida_Ns as usize;
            let IDAMem { ida_phiS, ida_eeS, ida_deltaS, .. } = ida_mem;
            for is in 0..ns {
                for (d, (p, e)) in ida_deltaS[is].data.iter_mut()
                    .zip(ida_phiS[kk][is].data.iter().zip(&ida_eeS[is].data))
                {
                    *d = *p + *e;
                }
            }
        }
        let errS_km1 = ida_mem.ida_sigma[ida_mem.ida_kk as usize - 1]
            * IDASensWrmsNorm(ida_mem, &ida_mem.ida_deltaS, &ida_mem.ida_ewtS,
                              ida_mem.ida_suppressalg);
        if errS_km1 > *err_km1 {
            *err_km1 = errS_km1;
            check_for_reduction = SUNTRUE;
        }
        let terr_km1 = ida_mem.ida_kk as f64 * *err_km1;

        /* Has an order decrease already been decided in IDATestError? */
        if ida_mem.ida_knew != ida_mem.ida_kk {
            check_for_reduction = SUNFALSE;
        }

        if check_for_reduction {
            if ida_mem.ida_kk > 2 {
                /* Update error at order k-2 */
                {
                    let kk = ida_mem.ida_kk as usize;
                    let ns = ida_mem.ida_Ns as usize;
                    let IDAMem { ida_phiS, ida_deltaS, .. } = ida_mem;
                    for is in 0..ns {
                        for (d, p) in ida_deltaS[is].data.iter_mut()
                            .zip(&ida_phiS[kk - 1][is].data)
                        {
                            *d = *p + *d;
                        }
                    }
                }
                let errS_km2 = ida_mem.ida_sigma[ida_mem.ida_kk as usize - 2]
                    * IDASensWrmsNorm(ida_mem, &ida_mem.ida_deltaS, &ida_mem.ida_ewtS,
                                      ida_mem.ida_suppressalg);
                if errS_km2 > *err_km2 {
                    *err_km2 = errS_km2;
                }
                let terr_km2 = (ida_mem.ida_kk - 1) as f64 * *err_km2;

                /* Decrease order if errors are reduced */
                if SUNMAX(terr_km1, terr_km2) <= terr_k {
                    ida_mem.ida_knew = ida_mem.ida_kk - 1;
                }
            } else {
                /* Decrease order to 1 if errors are reduced by at least 1/2 */
                if terr_km1 <= HALF * terr_k {
                    ida_mem.ida_knew = ida_mem.ida_kk - 1;
                }
            }
        }
    }

    /* Perform error test */
    if ck * enormS > ONE {
        ERROR_TEST_FAIL
    } else {
        IDA_SUCCESS
    }
}

/*
 * IDAQuadSensTestError
 *
 * This routine estimates quadrature sensitivity errors and updates
 * errors at orders k, k-1, k-2, decides whether or not to suggest
 * an order reduction and performs the local error test. (Used
 * only in staggered approach).
 *
 * IDAQuadSensTestError returns the updated local error estimate at
 * orders k, k-1, and k-2. These are norms of type
 * SUNMAX(|err|,|errQ|,|errS|,|errQS|).
 *
 * The return flag can be either IDA_SUCCESS or ERROR_TEST_FAIL.
 *
 * (C renames yyQS as tempv; ida_yyQS is used directly.)
 */
fn IDAQuadSensTestError(ida_mem: &mut IDAMem, ck: f64, err_k: &mut f64, err_km1: &mut f64,
                        err_km2: &mut f64) -> i32 {
    let mut check_for_reduction = SUNFALSE;

    let enormQS = IDAQuadSensWrmsNorm(ida_mem, &ida_mem.ida_eeQS, &ida_mem.ida_ewtQS);
    let errQS_k = ida_mem.ida_sigma[ida_mem.ida_kk as usize] * enormQS;

    if errQS_k > *err_k {
        *err_k = errQS_k;
        check_for_reduction = SUNTRUE;
    }
    let terr_k = (ida_mem.ida_kk + 1) as f64 * *err_k;

    if ida_mem.ida_kk > 1 {
        /* Update error at order k-1 */
        {
            let kk = ida_mem.ida_kk as usize;
            let ns = ida_mem.ida_Ns as usize;
            let IDAMem { ida_phiQS, ida_eeQS, ida_yyQS, .. } = ida_mem;
            for is in 0..ns {
                for (t, (p, e)) in ida_yyQS[is].data.iter_mut()
                    .zip(ida_phiQS[kk][is].data.iter().zip(&ida_eeQS[is].data))
                {
                    *t = *p + *e;
                }
            }
        }
        let errQS_km1 = ida_mem.ida_sigma[ida_mem.ida_kk as usize - 1]
            * IDAQuadSensWrmsNorm(ida_mem, &ida_mem.ida_yyQS, &ida_mem.ida_ewtQS);
        if errQS_km1 > *err_km1 {
            *err_km1 = errQS_km1;
            check_for_reduction = SUNTRUE;
        }
        let terr_km1 = ida_mem.ida_kk as f64 * *err_km1;

        /* Has an order decrease already been decided in IDATestError? */
        if ida_mem.ida_knew != ida_mem.ida_kk {
            check_for_reduction = SUNFALSE;
        }

        if check_for_reduction {
            if ida_mem.ida_kk > 2 {
                /* Update error at order k-2 */
                {
                    let kk = ida_mem.ida_kk as usize;
                    let ns = ida_mem.ida_Ns as usize;
                    let IDAMem { ida_phiQS, ida_yyQS, .. } = ida_mem;
                    for is in 0..ns {
                        for (t, p) in ida_yyQS[is].data.iter_mut()
                            .zip(&ida_phiQS[kk - 1][is].data)
                        {
                            *t = *p + *t;
                        }
                    }
                }
                let errQS_km2 = ida_mem.ida_sigma[ida_mem.ida_kk as usize - 2]
                    * IDAQuadSensWrmsNorm(ida_mem, &ida_mem.ida_yyQS, &ida_mem.ida_ewtQS);
                if errQS_km2 > *err_km2 {
                    *err_km2 = errQS_km2;
                }
                let terr_km2 = (ida_mem.ida_kk - 1) as f64 * *err_km2;

                /* Decrease order if errors are reduced */
                if SUNMAX(terr_km1, terr_km2) <= terr_k {
                    ida_mem.ida_knew = ida_mem.ida_kk - 1;
                }
            } else {
                /* Decrease order to 1 if errors are reduced by at least 1/2 */
                if terr_km1 <= HALF * terr_k {
                    ida_mem.ida_knew = ida_mem.ida_kk - 1;
                }
            }
        }
    }

    /* Perform error test */
    if ck * enormQS > ONE {
        ERROR_TEST_FAIL
    } else {
        IDA_SUCCESS
    }
}

/*
 * IDARestore
 *
 * This routine restores tn, psi, and phi in the event of a failure.
 * It changes back phi-star to phi (changed in IDASetCoeffs)
 *
 * (C stages ONE/beta[i] into ida_cvals — and gathers phiS/phiQS into
 * ida_Xvecs — for fused N_VScaleVectorArray calls; staging dropped,
 * per-vector scale kernels inline.)
 */
fn IDARestore(ida_mem: &mut IDAMem, saved_t: f64) {
    ida_mem.ida_tn = saved_t;

    let mut i: i32 = 1;
    while i <= ida_mem.ida_kk {
        ida_mem.ida_psi[i as usize - 1] = ida_mem.ida_psi[i as usize] - ida_mem.ida_hh;
        i += 1;
    }

    if ida_mem.ida_ns <= ida_mem.ida_kk {
        let ns_i = ida_mem.ida_ns as usize;
        let kk = ida_mem.ida_kk as usize;

        /* N_VScaleVectorArray(kk - ns + 1, 1/beta[ns..], phi + ns, phi + ns) */
        for i in ns_i..=kk {
            let c = ONE / ida_mem.ida_beta[i];
            ida_mem.ida_phi[i].scale_inplace(c);
        }

        if ida_mem.ida_quadr {
            for i in ns_i..=kk {
                let c = ONE / ida_mem.ida_beta[i];
                ida_mem.ida_phiQ[i].scale_inplace(c);
            }
        }

        /* (C's sensi||quadr_sensi cvals-staging loop vanishes) */

        if ida_mem.ida_sensi {
            for i in ns_i..=kk {
                let c = ONE / ida_mem.ida_beta[i];
                for is in 0..(ida_mem.ida_Ns as usize) {
                    ida_mem.ida_phiS[i][is].scale_inplace(c);
                }
            }
        }

        if ida_mem.ida_quadr_sensi {
            for i in ns_i..=kk {
                let c = ONE / ida_mem.ida_beta[i];
                for is in 0..(ida_mem.ida_Ns as usize) {
                    ida_mem.ida_phiQS[i][is].scale_inplace(c);
                }
            }
        }
    }
}

/*
 * -----------------------------------------------------------------
 * Handler for convergence and/or error test failures
 * -----------------------------------------------------------------
 */

/*
 * IDAHandleNFlag
 *
 * This routine handles failures indicated by the input variable nflag.
 * Positive values indicate various recoverable failures while negative
 * values indicate nonrecoverable failures. This routine adjusts the
 * step size for recoverable failures.
 *
 *  Possible nflag values (input):
 *
 *   --convergence failures--
 *   IDA_RES_RECVR              > 0
 *   IDA_LSOLVE_RECVR           > 0
 *   IDA_CONSTR_RECVR           > 0
 *   SUN_NLS_CONV_RECVR         > 0
 *   IDA_QRHS_RECVR             > 0
 *   IDA_QSRHS_RECVR            > 0
 *   IDA_RES_FAIL               < 0
 *   IDA_LSOLVE_FAIL            < 0
 *   IDA_LSETUP_FAIL            < 0
 *   IDA_QRHS_FAIL              < 0
 *
 *   --error test failure--
 *   ERROR_TEST_FAIL            > 0
 *
 *  Possible kflag values (output):
 *
 *   --recoverable--
 *   PREDICT_AGAIN
 *
 *   --nonrecoverable--
 *   IDA_REP_RES_ERR
 *   IDA_ERR_FAIL
 *   IDA_CONV_FAIL
 *   IDA_RES_FAIL
 *   IDA_LSETUP_FAIL
 *   IDA_LSOLVE_FAIL
 *   IDA_QRHS_FAIL
 *   IDA_REP_QRHS_ERR
 *
 * (C signature takes long* ncfnPtr / netfPtr pointing at IDAMem
 * counter fields; callers take the scalars out and restore them,
 * so plain &mut i64 here.  SUNLogDebug instrumentation dropped.)
 */
fn IDAHandleNFlag(ida_mem: &mut IDAMem, nflag: i32, err_k: f64, err_km1: f64,
                  ncfn_ptr: &mut i64, ncf_ptr: &mut i32, netf_ptr: &mut i64,
                  nef_ptr: &mut i32) -> i32 {
    ida_mem.ida_phase = 1;

    if nflag != ERROR_TEST_FAIL {
        /*-----------------------
          Nonlinear solver failed
          -----------------------*/

        *ncf_ptr += 1; /* local counter for convergence failures */
        *ncfn_ptr += 1; /* global counter for convergence failures */

        if nflag < 0 {
            /* nonrecoverable failure */

            if nflag == IDA_LSOLVE_FAIL {
                return IDA_LSOLVE_FAIL;
            } else if nflag == IDA_LSETUP_FAIL {
                return IDA_LSETUP_FAIL;
            } else if nflag == IDA_RES_FAIL {
                return IDA_RES_FAIL;
            } else if nflag == IDA_QRHS_FAIL {
                return IDA_QRHS_FAIL;
            } else if nflag == IDA_SRES_FAIL {
                return IDA_SRES_FAIL;
            } else if nflag == IDA_QSRHS_FAIL {
                return IDA_QSRHS_FAIL;
            } else {
                return IDA_NLS_FAIL;
            }
        } else {
            /* recoverable failure    */

            /* Test if there were too many convergence failures or |h| = hmin */
            if *ncf_ptr == ida_mem.ida_maxncf
                || SUNRabs(ida_mem.ida_hh) <= ida_mem.ida_hmin * ONEPSM
            {
                if nflag == IDA_RES_RECVR {
                    return IDA_REP_RES_ERR;
                }
                if nflag == IDA_QRHS_RECVR {
                    return IDA_REP_QRHS_ERR;
                }
                if nflag == IDA_SRES_RECVR {
                    return IDA_REP_SRES_ERR;
                }
                if nflag == IDA_QSRHS_RECVR {
                    return IDA_REP_QSRHS_ERR;
                }
                return IDA_CONV_FAIL;
            }

            /* Reduce step size for a new prediction */
            ida_mem.ida_eta = SUNMAX(ida_mem.ida_eta_cf,
                                     ida_mem.ida_hmin / SUNRabs(ida_mem.ida_hh));
            ida_mem.ida_hh *= ida_mem.ida_eta;

            PREDICT_AGAIN
        }
    } else {
        /*-----------------
          Error Test failed
          -----------------*/

        *nef_ptr += 1; /* local counter for error test failures */
        *netf_ptr += 1; /* global counter for error test failures */

        if *nef_ptr == 1 {
            /* On first error test failure, keep current order or lower order by one.
               Compute new stepsize based on differences of the solution. */

            let err_knew = if ida_mem.ida_kk == ida_mem.ida_knew { err_k } else { err_km1 };

            ida_mem.ida_kk = ida_mem.ida_knew;
            ida_mem.ida_eta = PT9 * SUNRpowerR(TWO * err_knew + PT0001,
                                               -ONE / (ida_mem.ida_kk + 1) as f64);
            ida_mem.ida_eta = SUNMAX(ida_mem.ida_eta_min_ef,
                                     SUNMIN(ida_mem.ida_eta_low, ida_mem.ida_eta));
            ida_mem.ida_eta = SUNMAX(ida_mem.ida_eta,
                                     ida_mem.ida_hmin / SUNRabs(ida_mem.ida_hh));
            ida_mem.ida_hh *= ida_mem.ida_eta;

            PREDICT_AGAIN
        } else if *nef_ptr == 2 {
            /* On second error test failure, use current order or decrease order by one.
               Reduce stepsize by factor of 1/4. */

            ida_mem.ida_kk = ida_mem.ida_knew;
            ida_mem.ida_eta = SUNMAX(ida_mem.ida_eta_min_ef,
                                     ida_mem.ida_hmin / SUNRabs(ida_mem.ida_hh));
            ida_mem.ida_hh *= ida_mem.ida_eta;

            PREDICT_AGAIN
        } else if *nef_ptr < ida_mem.ida_maxnef {
            /* On third and subsequent error test failures, set order to 1.
               Reduce stepsize by factor of 1/4. */
            ida_mem.ida_kk = 1;
            ida_mem.ida_eta = SUNMAX(ida_mem.ida_eta_min_ef,
                                     ida_mem.ida_hmin / SUNRabs(ida_mem.ida_hh));
            ida_mem.ida_hh *= ida_mem.ida_eta;

            PREDICT_AGAIN
        } else {
            /* Too many error test failures */
            IDA_ERR_FAIL
        }
    }
}

/*
 * IDAReset
 *
 * This routine is called only if we need to predict again at the
 * very first step. In such a case, reset phi[1] and psi[0].
 *
 * (C cvals staging for the fused N_VScaleVectorArray dropped.)
 */
fn IDAReset(ida_mem: &mut IDAMem) {
    ida_mem.ida_psi[0] = ida_mem.ida_hh;

    let eta = ida_mem.ida_eta;

    /* N_VScale(eta, phi[1], phi[1]) */
    ida_mem.ida_phi[1].scale_inplace(eta);

    if ida_mem.ida_quadr {
        ida_mem.ida_phiQ[1].scale_inplace(eta);
    }

    if ida_mem.ida_sensi {
        for is in 0..(ida_mem.ida_Ns as usize) {
            ida_mem.ida_phiS[1][is].scale_inplace(eta);
        }
    }

    if ida_mem.ida_quadr_sensi {
        for is in 0..(ida_mem.ida_Ns as usize) {
            ida_mem.ida_phiQS[1][is].scale_inplace(eta);
        }
    }
}

/*
 * IDACompleteStep
 *
 * This routine completes a successful step.  It increments nst,
 * saves the stepsize and order used, makes the final selection of
 * stepsize and order for the next step, and updates the phi array.
 *
 * (C stages the phi update through ida_Xvecs/ida_Zvecs and one fused
 * N_VLinearSumVectorArray whose serial kernel processes the pairs IN
 * ARRAY ORDER — so pair j reads the phi column UPDATED by pair j-1.
 * That sequential-aliasing cascade is replicated by iterating
 * j = 0..=kused in order with split borrows.  The unit-coefficient
 * N_VScale / N_VScaleVectorArray saves are copies, and the
 * cannot-fail VECTOROP_ERR branches vanish.)
 */
fn IDACompleteStep(ida_mem: &mut IDAMem, err_k: f64, err_km1: f64) {
    ida_mem.ida_nst += 1;
    let kdiff = ida_mem.ida_kk - ida_mem.ida_kused;
    ida_mem.ida_kused = ida_mem.ida_kk;
    ida_mem.ida_hused = ida_mem.ida_hh;

    if ida_mem.ida_knew == ida_mem.ida_kk - 1 || ida_mem.ida_kk == ida_mem.ida_maxord {
        ida_mem.ida_phase = 1;
    }

    /* For the first few steps, until either a step fails, or the order is
       reduced, or the order reaches its maximum, we raise the order and double
       the stepsize. During these steps, phase = 0. Thereafter, phase = 1, and
       stepsize and order are set by the usual local error algorithm.

       Note that, after the first step, the order is not increased, as not all
       of the neccessary information is available yet. */

    if ida_mem.ida_phase == 0 {
        if ida_mem.ida_nst > 1 {
            ida_mem.ida_kk += 1;
            let mut hnew = TWO * ida_mem.ida_hh;
            let tmp = SUNRabs(hnew) * ida_mem.ida_hmax_inv;
            if tmp > ONE {
                hnew /= tmp;
            }
            ida_mem.ida_hh = hnew;
        }
    } else {
        /* action = UNSET; (definite-assignment analysis replaces the
           C sentinel — every branch below assigns exactly once) */
        let action;
        let mut err_kp1 = ZERO;

        /* Set action = LOWER/MAINTAIN/RAISE to specify order decision */

        if ida_mem.ida_knew == ida_mem.ida_kk - 1 {
            action = LOWER;
        } else if ida_mem.ida_kk == ida_mem.ida_maxord {
            action = MAINTAIN;
        } else if ida_mem.ida_kk + 1 >= ida_mem.ida_ns || kdiff == 1 {
            action = MAINTAIN;
        } else {
            /* Consider order raise */

            /* Estimate the error at order k+1 */
            /* N_VLinearSum(ONE, ee, -ONE, phi[kk+1], tempv1) — VDiff */
            {
                let kk = ida_mem.ida_kk as usize;
                let IDAMem { ida_ee, ida_phi, ida_tempv1, .. } = ida_mem;
                for (t, (e, p)) in ida_tempv1.data.iter_mut()
                    .zip(ida_ee.data.iter().zip(&ida_phi[kk + 1].data))
                {
                    *t = *e - *p;
                }
            }
            let mut enorm = IDAWrmsNorm(ida_mem, &ida_mem.ida_tempv1, &ida_mem.ida_ewt,
                                        ida_mem.ida_suppressalg);

            if ida_mem.ida_errconQ {
                /* (C renames ypQ as tempvQ) */
                /* N_VLinearSum(ONE, eeQ, -ONE, phiQ[kk+1], tempvQ) — VDiff */
                {
                    let kk = ida_mem.ida_kk as usize;
                    let IDAMem { ida_eeQ, ida_phiQ, ida_ypQ, .. } = ida_mem;
                    for (t, (e, p)) in ida_ypQ.data.iter_mut()
                        .zip(ida_eeQ.data.iter().zip(&ida_phiQ[kk + 1].data))
                    {
                        *t = *e - *p;
                    }
                }
                enorm = IDAQuadWrmsNormUpdate(ida_mem, enorm, &ida_mem.ida_ypQ,
                                              &ida_mem.ida_ewtQ);
            }

            if ida_mem.ida_errconS {
                /* (C renames ypS as tempvS; fused per-is VDiff) */
                {
                    let kk = ida_mem.ida_kk as usize;
                    let ns = ida_mem.ida_Ns as usize;
                    let IDAMem { ida_eeS, ida_phiS, ida_ypS, .. } = ida_mem;
                    for is in 0..ns {
                        for (t, (e, p)) in ida_ypS[is].data.iter_mut()
                            .zip(ida_eeS[is].data.iter().zip(&ida_phiS[kk + 1][is].data))
                        {
                            *t = *e - *p;
                        }
                    }
                }
                enorm = IDASensWrmsNormUpdate(ida_mem, enorm, &ida_mem.ida_ypS,
                                              &ida_mem.ida_ewtS, ida_mem.ida_suppressalg);
            }

            if ida_mem.ida_errconQS {
                /* (per-is VDiff into tempvQS) */
                {
                    let kk = ida_mem.ida_kk as usize;
                    let ns = ida_mem.ida_Ns as usize;
                    let IDAMem { ida_eeQS, ida_phiQS, ida_tempvQS, .. } = ida_mem;
                    for is in 0..ns {
                        for (t, (e, p)) in ida_tempvQS[is].data.iter_mut()
                            .zip(ida_eeQS[is].data.iter().zip(&ida_phiQS[kk + 1][is].data))
                        {
                            *t = *e - *p;
                        }
                    }
                }
                enorm = IDAQuadSensWrmsNormUpdate(ida_mem, enorm, &ida_mem.ida_tempvQS,
                                                  &ida_mem.ida_ewtQS);
            }
            err_kp1 = enorm / (ida_mem.ida_kk + 2) as f64;

            /* Choose among orders k-1, k, k+1 using local truncation error norms. */

            let terr_k = (ida_mem.ida_kk + 1) as f64 * err_k;
            let terr_kp1 = (ida_mem.ida_kk + 2) as f64 * err_kp1;

            if ida_mem.ida_kk == 1 {
                if terr_kp1 >= HALF * terr_k {
                    action = MAINTAIN;
                } else {
                    action = RAISE;
                }
            } else {
                let terr_km1 = ida_mem.ida_kk as f64 * err_km1;
                if terr_km1 <= SUNMIN(terr_k, terr_kp1) {
                    action = LOWER;
                } else if terr_kp1 >= terr_k {
                    action = MAINTAIN;
                } else {
                    action = RAISE;
                }
            }
        }

        /* Set the estimated error norm and, on change of order, reset kk. */
        let err_knew;
        if action == RAISE {
            ida_mem.ida_kk += 1;
            err_knew = err_kp1;
        } else if action == LOWER {
            ida_mem.ida_kk -= 1;
            err_knew = err_km1;
        } else {
            err_knew = err_k;
        }

        /* Compute tmp = tentative ratio hnew/hh from error norm estimate.
           1. If eta >= eta_max_fx (default = 2), the step size is doubled or, if
              eta_max was set to a value < 2, increased to at most eta_max with
              the maximum step size bound enforced.
           2. If eta <= eta_min_fx (default = 1), the step size is reduced to at
              most eta_low with the minimum eta and step size bounds enforced.
         */

        ida_mem.ida_eta = ONE;
        let tmp = SUNRpowerR(TWO * err_knew + PT0001, -ONE / (ida_mem.ida_kk + 1) as f64);

        if tmp >= ida_mem.ida_eta_max_fx {
            /* Enforce max growth factor bound and max step size */
            ida_mem.ida_eta = SUNMIN(tmp, ida_mem.ida_eta_max);
            ida_mem.ida_eta /= SUNMAX(ONE, ida_mem.ida_eta * SUNRabs(ida_mem.ida_hh)
                                           * ida_mem.ida_hmax_inv);
        } else if tmp <= ida_mem.ida_eta_min_fx {
            /* Enforce min eta bound, min step size */
            ida_mem.ida_eta = SUNMIN(tmp, ida_mem.ida_eta_low);
            ida_mem.ida_eta = SUNMAX(ida_mem.ida_eta, ida_mem.ida_eta_min);
            ida_mem.ida_eta = SUNMAX(ida_mem.ida_eta,
                                     ida_mem.ida_hmin / SUNRabs(ida_mem.ida_hh));
        }
        ida_mem.ida_hh *= ida_mem.ida_eta;
    } /* end of phase if block */

    /* Save ee for possible order increase on next step */
    if ida_mem.ida_kused < ida_mem.ida_maxord {
        /* N_VScale(ONE, ee, phi[kused+1]) — copy */
        {
            let ku = ida_mem.ida_kused as usize;
            let IDAMem { ida_ee, ida_phi, .. } = ida_mem;
            ida_phi[ku + 1].data.copy_from_slice(&ida_ee.data);
        }

        if ida_mem.ida_quadr {
            let ku = ida_mem.ida_kused as usize;
            let IDAMem { ida_eeQ, ida_phiQ, .. } = ida_mem;
            ida_phiQ[ku + 1].data.copy_from_slice(&ida_eeQ.data);
        }

        /* (C's sensi||quadr_sensi cvals-staging loop vanishes) */

        if ida_mem.ida_sensi {
            let ku = ida_mem.ida_kused as usize;
            let ns = ida_mem.ida_Ns as usize;
            let IDAMem { ida_eeS, ida_phiS, .. } = ida_mem;
            for is in 0..ns {
                ida_phiS[ku + 1][is].data.copy_from_slice(&ida_eeS[is].data);
            }
        }

        if ida_mem.ida_quadr_sensi {
            let ku = ida_mem.ida_kused as usize;
            let ns = ida_mem.ida_Ns as usize;
            let IDAMem { ida_eeQS, ida_phiQS, .. } = ida_mem;
            for is in 0..ns {
                ida_phiQS[ku + 1][is].data.copy_from_slice(&ida_eeQS[is].data);
            }
        }
    }

    /* Update phi arrays */

    /* To avoid scanning the phi array twice, the fused operation
       X = X + Z is used, with
         X = [ phi[kused], phi[kused-1], ..., phi[0] ]
         Z = [ ee,         phi[kused],   ..., phi[1] ]
       processed in order, so each pair reads the previous update. */
    {
        let ku = ida_mem.ida_kused as usize;

        /* j = 0: phi[kused] += ee */
        {
            let IDAMem { ida_phi, ida_ee, .. } = ida_mem;
            for (x, z) in ida_phi[ku].data.iter_mut().zip(&ida_ee.data) {
                *x = *x + *z;
            }
        }
        /* j = 1..=kused: phi[kused-j] += phi[kused-j+1] (the UPDATED column) */
        for j in 1..=ku {
            let (lo, hi) = ida_mem.ida_phi.split_at_mut(ku - j + 1);
            let x = &mut lo[ku - j];
            let z = &hi[0];
            for (xi, zi) in x.data.iter_mut().zip(&z.data) {
                *xi = *xi + *zi;
            }
        }
    }

    if ida_mem.ida_quadr {
        let ku = ida_mem.ida_kused as usize;

        /* j = 0: phiQ[kused] += eeQ */
        {
            let IDAMem { ida_phiQ, ida_eeQ, .. } = ida_mem;
            for (x, z) in ida_phiQ[ku].data.iter_mut().zip(&ida_eeQ.data) {
                *x = *x + *z;
            }
        }
        for j in 1..=ku {
            let (lo, hi) = ida_mem.ida_phiQ.split_at_mut(ku - j + 1);
            let x = &mut lo[ku - j];
            let z = &hi[0];
            for (xi, zi) in x.data.iter_mut().zip(&z.data) {
                *xi = *xi + *zi;
            }
        }
    }

    if ida_mem.ida_sensi {
        let ku = ida_mem.ida_kused as usize;
        let ns = ida_mem.ida_Ns as usize;

        /* (C interleaves the per-is cascades into one fused call over
           Ns*(kused+1) pairs; the cascades are independent across is
           and ordered within each is) */
        for is in 0..ns {
            /* j = 0: phiS[kused][is] += eeS[is] */
            {
                let IDAMem { ida_phiS, ida_eeS, .. } = ida_mem;
                for (x, z) in ida_phiS[ku][is].data.iter_mut().zip(&ida_eeS[is].data) {
                    *x = *x + *z;
                }
            }
            for j in 1..=ku {
                let (lo, hi) = ida_mem.ida_phiS.split_at_mut(ku - j + 1);
                let x = &mut lo[ku - j][is];
                let z = &hi[0][is];
                for (xi, zi) in x.data.iter_mut().zip(&z.data) {
                    *xi = *xi + *zi;
                }
            }
        }
    }

    if ida_mem.ida_quadr_sensi {
        let ku = ida_mem.ida_kused as usize;
        let ns = ida_mem.ida_Ns as usize;

        for is in 0..ns {
            /* j = 0: phiQS[kused][is] += eeQS[is] */
            {
                let IDAMem { ida_phiQS, ida_eeQS, .. } = ida_mem;
                for (x, z) in ida_phiQS[ku][is].data.iter_mut().zip(&ida_eeQS[is].data) {
                    *x = *x + *z;
                }
            }
            for j in 1..=ku {
                let (lo, hi) = ida_mem.ida_phiQS.split_at_mut(ku - j + 1);
                let x = &mut lo[ku - j][is];
                let z = &hi[0][is];
                for (xi, zi) in x.data.iter_mut().zip(&z.data) {
                    *xi = *xi + *zi;
                }
            }
        }
    }
}

/*
 * -----------------------------------------------------------------
 * Interpolated output
 * -----------------------------------------------------------------
 */

/*
 * IDAGetSolution
 *
 * This routine evaluates y(t) and y'(t) as the value and derivative of
 * the interpolating polynomial at the independent variable t, and stores
 * the results in the vectors yret and ypret.  It uses the current
 * independent variable value, tn, and the method order last used, kused.
 * This function is called by IDASolve with t = tout, t = tn, or t = tstop.
 *
 * If kused = 0 (no step has been taken), or if t = tn, then the order used
 * here is taken to be 1, giving yret = phi[0], ypret = phi[1]/psi[0].
 *
 * The return values are:
 *   IDA_SUCCESS  if t is legal, or
 *   IDA_BAD_T    if t is not within the interval of the last step taken.
 *   (the serial linear combinations below cannot fail, so the
 *    IDA_VECTOROP_ERR branches of the C code vanish)
 *
 * (C stages the coefficients in ida_cvals/ida_dvals and calls the
 * fused N_VLinearCombination twice; the work arrays are locals here —
 * the cvals/dvals/Xvecs/Zvecs IDAMem fields are not stored, per the
 * dropped-fused-ops convention — and the serial kernels are applied
 * inline.)
 */
pub fn IDAGetSolution(ida_mem: &IDAMem, t: f64, yret: &mut NVector, ypret: &mut NVector) -> i32 {
    /* Check t for legality.  Here tn - hused is t_{n-1}. */

    let mut tfuzz = HUNDRED * ida_mem.ida_uround *
        (SUNRabs(ida_mem.ida_tn) + SUNRabs(ida_mem.ida_hh));
    if ida_mem.ida_hh < ZERO {
        tfuzz = -tfuzz;
    }
    let tp = ida_mem.ida_tn - ida_mem.ida_hused - tfuzz;
    if (t - tp) * ida_mem.ida_hh < ZERO {
        IDAProcessError(Some(ida_mem), IDA_BAD_T, line!(), "IDAGetSolution", file!(),
                        &ida_msg_g(MSG_BAD_T,
                                   &[t, ida_mem.ida_tn - ida_mem.ida_hused, ida_mem.ida_tn]));
        return IDA_BAD_T;
    }

    /* Initialize kord = (kused or 1). */

    let mut kord = ida_mem.ida_kused;
    if ida_mem.ida_kused == 0 {
        kord = 1;
    }

    /* Accumulate multiples of columns phi[j] into yret and ypret. */

    let delt = t - ida_mem.ida_tn;
    let mut c = ONE;
    let mut d = ZERO;
    let mut gam = delt / ida_mem.ida_psi[0];

    let mut cvals = [ZERO; MXORDP1];
    let mut dvals = [ZERO; MXORDP1];

    cvals[0] = c;
    let mut j: i32 = 1;
    while j <= kord {
        d = d * gam + c / ida_mem.ida_psi[j as usize - 1];
        c = c * gam;
        gam = (delt + ida_mem.ida_psi[j as usize - 1]) / ida_mem.ida_psi[j as usize];

        cvals[j as usize] = c;
        dvals[j as usize - 1] = d;
        j += 1;
    }

    /* retval = N_VLinearCombination(kord + 1, cvals, phi, yret)
       (serial kernel: yret = cvals[0]*phi[0], then += cvals[j]*phi[j]) */
    let kord_u = kord as usize;
    for (y, p) in yret.data.iter_mut().zip(&ida_mem.ida_phi[0].data) {
        *y = cvals[0] * *p;
    }
    for j in 1..=kord_u {
        for (y, p) in yret.data.iter_mut().zip(&ida_mem.ida_phi[j].data) {
            *y += cvals[j] * *p;
        }
    }

    /* retval = N_VLinearCombination(kord, dvals, phi + 1, ypret) */
    for (y, p) in ypret.data.iter_mut().zip(&ida_mem.ida_phi[1].data) {
        *y = dvals[0] * *p;
    }
    for j in 2..=kord_u {
        for (y, p) in ypret.data.iter_mut().zip(&ida_mem.ida_phi[j].data) {
            *y += dvals[j - 1] * *p;
        }
    }

    IDA_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Norm functions
 * -----------------------------------------------------------------
 */

/*
 * IDAWrmsNorm
 *
 *  Returns the WRMS norm of vector x with weights w.
 *  If mask = SUNTRUE, the weight vector w is masked by id, i.e.,
 *      nrm = N_VWrmsNormMask(x,w,id);
 *  Otherwise,
 *      nrm = N_VWrmsNorm(x,w);
 *
 * mask = SUNFALSE       when the call is made from the nonlinear solver.
 * mask = suppressalg otherwise.
 */
pub fn IDAWrmsNorm(ida_mem: &IDAMem, x: &NVector, w: &NVector, mask: bool) -> f64 {
    if mask {
        N_VWrmsNormMask(x, w, &ida_mem.ida_id)
    } else {
        N_VWrmsNorm(x, w)
    }
}

/*
 * IDASensWrmsNorm
 *
 * This routine returns the maximum over the weighted root mean
 * square norm of xS with weight vectors wS:
 *
 *   max { wrms(xS[0],wS[0]) ... wrms(xS[Ns-1],wS[Ns-1]) }
 *
 * Called by IDASensUpdateNorm or directly in the IDA_STAGGERED approach
 * during the NLS solution and before the error test.
 *
 * Declared global for use in the computation of IC for sensitivities.
 *
 * (C computes the per-vector norms into ida_cvals via the fused
 * N_VWrmsNorm(Mask)VectorArray and then scans for the maximum; the
 * fused serial kernel is a per-vector norm, computed here into a
 * local buffer.)
 */
pub fn IDASensWrmsNorm(ida_mem: &IDAMem, xS: &[NVector], wS: &[NVector], mask: bool) -> f64 {
    let ns = ida_mem.ida_Ns as usize;
    let mut cvals = vec![ZERO; ns];

    if mask {
        for is in 0..ns {
            cvals[is] = N_VWrmsNormMask(&xS[is], &wS[is], &ida_mem.ida_id);
        }
    } else {
        for is in 0..ns {
            cvals[is] = N_VWrmsNorm(&xS[is], &wS[is]);
        }
    }

    let mut nrm = cvals[0];
    for is in 1..ns {
        if cvals[is] > nrm {
            nrm = cvals[is];
        }
    }

    nrm
}

/*
 * IDAQuadSensWrmsNorm
 *
 * This routine returns the maximum over the weighted root mean
 * square norm of xQS with weight vectors wQS:
 *
 *   max { wrms(xQS[0],wQS[0]) ... wrms(xQS[Ns-1],wQS[Ns-1]) }
 */
fn IDAQuadSensWrmsNorm(ida_mem: &IDAMem, xQS: &[NVector], wQS: &[NVector]) -> f64 {
    let ns = ida_mem.ida_Ns as usize;
    let mut cvals = vec![ZERO; ns];

    for is in 0..ns {
        cvals[is] = N_VWrmsNorm(&xQS[is], &wQS[is]);
    }

    let mut nrm = cvals[0];
    for is in 1..ns {
        if cvals[is] > nrm {
            nrm = cvals[is];
        }
    }

    nrm
}

/*
 * IDAQuadWrmsNormUpdate
 *
 * Updates the norm old_nrm to account for all quadratures.
 * (C marks IDA_mem SUNDIALS_MAYBE_UNUSED — underscore here.)
 */
fn IDAQuadWrmsNormUpdate(_ida_mem: &IDAMem, old_nrm: f64, xQ: &NVector, wQ: &NVector) -> f64 {
    let qnrm = N_VWrmsNorm(xQ, wQ);
    if old_nrm > qnrm {
        old_nrm
    } else {
        qnrm
    }
}

/*
 * IDASensWrmsNormUpdate
 *
 * Updates the norm old_nrm to account for all sensitivities.
 *
 * This function is declared global since it is used for finding
 * IC for sensitivities,
 */
pub fn IDASensWrmsNormUpdate(ida_mem: &IDAMem, old_nrm: f64, xS: &[NVector], wS: &[NVector],
                             mask: bool) -> f64 {
    let snrm = IDASensWrmsNorm(ida_mem, xS, wS, mask);
    if old_nrm > snrm {
        old_nrm
    } else {
        snrm
    }
}

fn IDAQuadSensWrmsNormUpdate(ida_mem: &IDAMem, old_nrm: f64, xQS: &[NVector],
                             wQS: &[NVector]) -> f64 {
    let qsnrm = IDAQuadSensWrmsNorm(ida_mem, xQS, wQS);
    if old_nrm > qsnrm {
        old_nrm
    } else {
        qsnrm
    }
}

/*
 * -----------------------------------------------------------------
 * Functions for rootfinding
 * -----------------------------------------------------------------
 */

/*
 * IDARcheck1
 *
 * This routine completes the initialization of rootfinding memory
 * information, and checks whether g has a zero both at and very near
 * the initial point of the IVP.
 *
 * This routine returns an int equal to:
 *  IDA_RTFUNC_FAIL < 0 if the g function failed, or
 *  IDA_SUCCESS     = 0 otherwise.
 */
fn IDARcheck1(ida_mem: &mut IDAMem) -> i32 {
    let nrtfn = ida_mem.ida_nrtfn as usize;

    for i in 0..nrtfn {
        ida_mem.ida_iroots[i] = 0;
    }
    ida_mem.ida_tlo = ida_mem.ida_tn;
    ida_mem.ida_ttol = (SUNRabs(ida_mem.ida_tn) + SUNRabs(ida_mem.ida_hh))
        * ida_mem.ida_uround * HUNDRED;

    /* Evaluate g at initial t and check for zero values. */
    let retval = {
        let gfun = ida_mem.ida_gfun.unwrap();
        let tlo = ida_mem.ida_tlo;
        let IDAMem { ida_phi, ida_glo, ida_user_data, .. } = ida_mem;
        gfun(tlo, &ida_phi[0], &ida_phi[1], ida_glo, ida_user_data)
    };
    ida_mem.ida_nge = 1;
    if retval != 0 {
        return IDA_RTFUNC_FAIL;
    }

    let mut zroot = SUNFALSE;
    for i in 0..nrtfn {
        if SUNRabs(ida_mem.ida_glo[i]) == ZERO {
            zroot = SUNTRUE;
            ida_mem.ida_gactive[i] = SUNFALSE;
        }
    }
    if !zroot {
        return IDA_SUCCESS;
    }

    /* Some g_i is zero at t0; look at g at t0+(small increment). */
    let hratio = SUNMAX(ida_mem.ida_ttol / SUNRabs(ida_mem.ida_hh), PT1);
    let smallh = hratio * ida_mem.ida_hh;
    let tplus = ida_mem.ida_tlo + smallh;
    /* N_VLinearSum(ONE, phi[0], smallh, phi[1], yy) — serial kernel
       VLin1(b, y, x, z): z = b*y + x (writes the owned ida_yy) */
    {
        let IDAMem { ida_phi, ida_yy, .. } = ida_mem;
        for (y, (p0, p1)) in ida_yy.data.iter_mut()
            .zip(ida_phi[0].data.iter().zip(&ida_phi[1].data))
        {
            *y = smallh * *p1 + *p0;
        }
    }
    let retval = {
        let gfun = ida_mem.ida_gfun.unwrap();
        let IDAMem { ida_yy, ida_phi, ida_ghi, ida_user_data, .. } = ida_mem;
        gfun(tplus, ida_yy, &ida_phi[1], ida_ghi, ida_user_data)
    };
    ida_mem.ida_nge += 1;
    if retval != 0 {
        return IDA_RTFUNC_FAIL;
    }

    /* We check now only the components of g which were exactly 0.0 at t0
     * to see if we can 'activate' them. */
    for i in 0..nrtfn {
        if !ida_mem.ida_gactive[i] && SUNRabs(ida_mem.ida_ghi[i]) != ZERO {
            ida_mem.ida_gactive[i] = SUNTRUE;
            ida_mem.ida_glo[i] = ida_mem.ida_ghi[i];
        }
    }
    IDA_SUCCESS
}

/*
 * IDARcheck2
 *
 * This routine checks for exact zeros of g at the last root found,
 * if the last return was a root.  It then checks for a close pair of
 * zeros (an error condition), and for a new root at a nearby point.
 * The array glo = g(tlo) at the left endpoint of the search interval
 * is adjusted if necessary to assure that all g_i are nonzero
 * there, before returning to do a root search in the interval.
 *
 * On entry, tlo = tretlast is the last value of tret returned by
 * IDASolve.  This may be the previous tn, the previous tout value,
 * or the last root location.
 *
 * This routine returns an int equal to:
 *     IDA_RTFUNC_FAIL < 0 if the g function failed, or
 *     CLOSERT         = 3 if a close pair of zeros was found, or
 *     RTFOUND         = 1 if a new zero of g was found near tlo, or
 *     IDA_SUCCESS     = 0 otherwise.
 */
fn IDARcheck2(ida_mem: &mut IDAMem) -> i32 {
    if ida_mem.ida_irfnd == 0 {
        return IDA_SUCCESS;
    }

    {
        let mut yy = std::mem::take(&mut ida_mem.ida_yy);
        let mut yp = std::mem::take(&mut ida_mem.ida_yp);
        let tlo = ida_mem.ida_tlo;
        let _ = IDAGetSolution(ida_mem, tlo, &mut yy, &mut yp);
        ida_mem.ida_yy = yy;
        ida_mem.ida_yp = yp;
    }
    let retval = {
        let gfun = ida_mem.ida_gfun.unwrap();
        let tlo = ida_mem.ida_tlo;
        let IDAMem { ida_yy, ida_yp, ida_glo, ida_user_data, .. } = ida_mem;
        gfun(tlo, ida_yy, ida_yp, ida_glo, ida_user_data)
    };
    ida_mem.ida_nge += 1;
    if retval != 0 {
        return IDA_RTFUNC_FAIL;
    }

    let nrtfn = ida_mem.ida_nrtfn as usize;
    let mut zroot = SUNFALSE;
    for i in 0..nrtfn {
        ida_mem.ida_iroots[i] = 0;
    }
    for i in 0..nrtfn {
        if !ida_mem.ida_gactive[i] {
            continue;
        }
        if SUNRabs(ida_mem.ida_glo[i]) == ZERO {
            zroot = SUNTRUE;
            ida_mem.ida_iroots[i] = 1;
        }
    }
    if !zroot {
        return IDA_SUCCESS;
    }

    /* One or more g_i has a zero at tlo.  Check g at tlo+smallh. */
    ida_mem.ida_ttol = (SUNRabs(ida_mem.ida_tn) + SUNRabs(ida_mem.ida_hh))
        * ida_mem.ida_uround * HUNDRED;
    let smallh = if ida_mem.ida_hh > ZERO { ida_mem.ida_ttol } else { -ida_mem.ida_ttol };
    let tplus = ida_mem.ida_tlo + smallh;
    if (tplus - ida_mem.ida_tn) * ida_mem.ida_hh >= ZERO {
        let hratio = smallh / ida_mem.ida_hh;
        /* N_VLinearSum(ONE, yy, hratio, phi[1], yy) — VLin1(b, y, x, z):
           z = b*y + x, with x aliasing z (the owned ida_yy) */
        let IDAMem { ida_yy, ida_phi, .. } = ida_mem;
        for (y, p1) in ida_yy.data.iter_mut().zip(&ida_phi[1].data) {
            *y = hratio * *p1 + *y;
        }
    } else {
        let mut yy = std::mem::take(&mut ida_mem.ida_yy);
        let mut yp = std::mem::take(&mut ida_mem.ida_yp);
        let _ = IDAGetSolution(ida_mem, tplus, &mut yy, &mut yp);
        ida_mem.ida_yy = yy;
        ida_mem.ida_yp = yp;
    }
    let retval = {
        let gfun = ida_mem.ida_gfun.unwrap();
        let IDAMem { ida_yy, ida_yp, ida_ghi, ida_user_data, .. } = ida_mem;
        gfun(tplus, ida_yy, ida_yp, ida_ghi, ida_user_data)
    };
    ida_mem.ida_nge += 1;
    if retval != 0 {
        return IDA_RTFUNC_FAIL;
    }

    /* Check for close roots (error return), for a new zero at tlo+smallh,
    and for a g_i that changed from zero to nonzero. */
    zroot = SUNFALSE;
    for i in 0..nrtfn {
        if !ida_mem.ida_gactive[i] {
            continue;
        }
        if SUNRabs(ida_mem.ida_ghi[i]) == ZERO {
            if ida_mem.ida_iroots[i] == 1 {
                return CLOSERT;
            }
            zroot = SUNTRUE;
            ida_mem.ida_iroots[i] = 1;
        } else {
            if ida_mem.ida_iroots[i] == 1 {
                ida_mem.ida_glo[i] = ida_mem.ida_ghi[i];
            }
        }
    }
    if zroot {
        return RTFOUND;
    }
    IDA_SUCCESS
}

/*
 * IDARcheck3
 *
 * This routine interfaces to IDARootfind to look for a root of g
 * between tlo and either tn or tout, whichever comes first.
 * Only roots beyond tlo in the direction of integration are sought.
 *
 * This routine returns an int equal to:
 *     IDA_RTFUNC_FAIL < 0 if the g function failed, or
 *     RTFOUND         = 1 if a root of g was found, or
 *     IDA_SUCCESS     = 0 otherwise.
 */
fn IDARcheck3(ida_mem: &mut IDAMem, tout: f64, itask: i32) -> i32 {
    /* Set thi = tn or tout, whichever comes first. */
    if itask == IDA_ONE_STEP {
        ida_mem.ida_thi = ida_mem.ida_tn;
    }
    if itask == IDA_NORMAL {
        ida_mem.ida_thi = if (tout - ida_mem.ida_tn) * ida_mem.ida_hh >= ZERO {
            ida_mem.ida_tn
        } else {
            tout
        };
    }

    /* Get y and y' at thi. */
    {
        let mut yy = std::mem::take(&mut ida_mem.ida_yy);
        let mut yp = std::mem::take(&mut ida_mem.ida_yp);
        let thi = ida_mem.ida_thi;
        let _ = IDAGetSolution(ida_mem, thi, &mut yy, &mut yp);
        ida_mem.ida_yy = yy;
        ida_mem.ida_yp = yp;
    }

    /* Set ghi = g(thi) and call IDARootfind to search (tlo,thi) for roots. */
    let retval = {
        let gfun = ida_mem.ida_gfun.unwrap();
        let thi = ida_mem.ida_thi;
        let IDAMem { ida_yy, ida_yp, ida_ghi, ida_user_data, .. } = ida_mem;
        gfun(thi, ida_yy, ida_yp, ida_ghi, ida_user_data)
    };
    ida_mem.ida_nge += 1;
    if retval != 0 {
        return IDA_RTFUNC_FAIL;
    }

    ida_mem.ida_ttol = (SUNRabs(ida_mem.ida_tn) + SUNRabs(ida_mem.ida_hh))
        * ida_mem.ida_uround * HUNDRED;
    let ier = IDARootfind(ida_mem);
    if ier == IDA_RTFUNC_FAIL {
        return IDA_RTFUNC_FAIL;
    }
    let nrtfn = ida_mem.ida_nrtfn as usize;
    for i in 0..nrtfn {
        if !ida_mem.ida_gactive[i] && ida_mem.ida_grout[i] != ZERO {
            ida_mem.ida_gactive[i] = SUNTRUE;
        }
    }
    ida_mem.ida_tlo = ida_mem.ida_trout;
    for i in 0..nrtfn {
        ida_mem.ida_glo[i] = ida_mem.ida_grout[i];
    }

    /* If no root found, return IDA_SUCCESS. */
    if ier == IDA_SUCCESS {
        return IDA_SUCCESS;
    }

    /* If a root was found, interpolate to get y(trout) and return.  */
    {
        let mut yy = std::mem::take(&mut ida_mem.ida_yy);
        let mut yp = std::mem::take(&mut ida_mem.ida_yp);
        let trout = ida_mem.ida_trout;
        let _ = IDAGetSolution(ida_mem, trout, &mut yy, &mut yp);
        ida_mem.ida_yy = yy;
        ida_mem.ida_yp = yp;
    }
    RTFOUND
}

/*
 * IDARootfind
 *
 * This routine solves for a root of g(t) between tlo and thi, if
 * one exists.  Only roots of odd multiplicity (i.e. with a change
 * of sign in one of the g_i), or exact zeros, are found.
 * Here the sign of tlo - thi is arbitrary, but if multiple roots
 * are found, the one closest to tlo is returned.
 *
 * The method used is the Illinois algorithm, a modified secant method.
 * Reference: Kathie L. Hiebert and Lawrence F. Shampine, Implicitly
 * Defined Output Points for Solutions of ODEs, Sandia National
 * Laboratory Report SAND80-0180, February 1980.
 *
 * This routine uses the following parameters for communication:
 *
 * nrtfn    = number of functions g_i, or number of components of
 *            the vector-valued function g(t).  Input only.
 *
 * gfun     = user-defined function for g(t).  Its form is
 *            (void) gfun(t, y, yp, gt, user_data)
 *
 * rootdir  = in array specifying the direction of zero-crossings.
 *            If rootdir[i] > 0, search for roots of g_i only if
 *            g_i is increasing; if rootdir[i] < 0, search for
 *            roots of g_i only if g_i is decreasing; otherwise
 *            always search for roots of g_i.
 *
 * gactive  = array specifying whether a component of g should
 *            or should not be monitored. gactive[i] is initially
 *            set to SUNTRUE for all i=0,...,nrtfn-1, but it may be
 *            reset to SUNFALSE if at the first step g[i] is 0.0
 *            both at the I.C. and at a small perturbation of them.
 *            gactive[i] is then set back on SUNTRUE only after the
 *            corresponding g function moves away from 0.0.
 *
 * nge      = cumulative counter for gfun calls.
 *
 * ttol     = a convergence tolerance for trout.  Input only.
 *            When a root at trout is found, it is located only to
 *            within a tolerance of ttol.  Typically, ttol should
 *            be set to a value on the order of
 *               100 * UROUND * max (SUNRabs(tlo), SUNRabs(thi))
 *            where UROUND is the unit roundoff of the machine.
 *
 * tlo, thi = endpoints of the interval in which roots are sought.
 *            On input, these must be distinct, but tlo - thi may
 *            be of either sign.  The direction of integration is
 *            assumed to be from tlo to thi.  On return, tlo and thi
 *            are the endpoints of the final relevant interval.
 *
 * glo, ghi = arrays of length nrtfn containing the vectors g(tlo)
 *            and g(thi) respectively.  Input and output.  On input,
 *            none of the glo[i] should be zero.
 *
 * trout    = root location, if a root was found, or thi if not.
 *            Output only.  If a root was found other than an exact
 *            zero of g, trout is the endpoint thi of the final
 *            interval bracketing the root, with size at most ttol.
 *
 * grout    = array of length nrtfn containing g(trout) on return.
 *
 * iroots   = int array of length nrtfn with root information.
 *            Output only.  If a root was found, iroots indicates
 *            which components g_i have a root at trout.  For
 *            i = 0, ..., nrtfn-1, iroots[i] = 1 if g_i has a root
 *            and g_i is increasing, iroots[i] = -1 if g_i has a
 *            root and g_i is decreasing, and iroots[i] = 0 if g_i
 *            has no roots or g_i varies in the direction opposite
 *            to that indicated by rootdir[i].
 *
 * This routine returns an int equal to:
 *      IDA_RTFUNC_FAIL < 0 if the g function failed, or
 *      RTFOUND         = 1 if a root of g was found, or
 *      IDA_SUCCESS     = 0 otherwise.
 *
 */
fn IDARootfind(ida_mem: &mut IDAMem) -> i32 {
    let nrtfn = ida_mem.ida_nrtfn as usize;
    let mut imax: usize = 0;

    /* First check for change in sign in ghi or for a zero in ghi. */
    let mut maxfrac = ZERO;
    let mut zroot = SUNFALSE;
    let mut sgnchg = SUNFALSE;
    for i in 0..nrtfn {
        if !ida_mem.ida_gactive[i] {
            continue;
        }
        if SUNRabs(ida_mem.ida_ghi[i]) == ZERO {
            if ida_mem.ida_rootdir[i] as f64 * ida_mem.ida_glo[i] <= ZERO {
                zroot = SUNTRUE;
            }
        } else {
            if SUNRdifferentsign(ida_mem.ida_glo[i], ida_mem.ida_ghi[i])
                && ida_mem.ida_rootdir[i] as f64 * ida_mem.ida_glo[i] <= ZERO
            {
                let gfrac = SUNRabs(ida_mem.ida_ghi[i]
                                    / (ida_mem.ida_ghi[i] - ida_mem.ida_glo[i]));
                if gfrac > maxfrac {
                    sgnchg = SUNTRUE;
                    maxfrac = gfrac;
                    imax = i;
                }
            }
        }
    }

    /* If no sign change was found, reset trout and grout.  Then return
       IDA_SUCCESS if no zero was found, or set iroots and return RTFOUND.  */
    if !sgnchg {
        ida_mem.ida_trout = ida_mem.ida_thi;
        for i in 0..nrtfn {
            ida_mem.ida_grout[i] = ida_mem.ida_ghi[i];
        }
        if !zroot {
            return IDA_SUCCESS;
        }
        for i in 0..nrtfn {
            ida_mem.ida_iroots[i] = 0;
            if !ida_mem.ida_gactive[i] {
                continue;
            }
            if SUNRabs(ida_mem.ida_ghi[i]) == ZERO
                && ida_mem.ida_rootdir[i] as f64 * ida_mem.ida_glo[i] <= ZERO
            {
                ida_mem.ida_iroots[i] = if ida_mem.ida_glo[i] > ZERO { -1 } else { 1 };
            }
        }
        return RTFOUND;
    }

    /* Initialize alph to avoid compiler warning */
    let mut alph = ONE;

    /* A sign change was found.  Loop to locate nearest root. */

    let mut side: i32 = 0;
    let mut sideprev: i32 = -1;
    loop {
        /* Looping point */

        /* If interval size is already less than tolerance ttol, break. */
        if SUNRabs(ida_mem.ida_thi - ida_mem.ida_tlo) <= ida_mem.ida_ttol {
            break;
        }

        /* Set weight alph.
           On the first two passes, set alph = 1.  Thereafter, reset alph
           according to the side (low vs high) of the subinterval in which
           the sign change was found in the previous two passes.
           If the sides were opposite, set alph = 1.
           If the sides were the same, then double alph (if high side),
           or halve alph (if low side).
           The next guess tmid is the secant method value if alph = 1, but
           is closer to tlo if alph < 1, and closer to thi if alph > 1.    */

        if sideprev == side {
            alph = if side == 2 { alph * TWO } else { alph * HALF };
        } else {
            alph = ONE;
        }

        /* Set next root approximation tmid and get g(tmid).
           If tmid is too close to tlo or thi, adjust it inward,
           by a fractional distance that is between 0.1 and 0.5.  */
        let mut tmid = ida_mem.ida_thi
            - (ida_mem.ida_thi - ida_mem.ida_tlo) * ida_mem.ida_ghi[imax]
                / (ida_mem.ida_ghi[imax] - alph * ida_mem.ida_glo[imax]);
        if SUNRabs(tmid - ida_mem.ida_tlo) < HALF * ida_mem.ida_ttol {
            let fracint = SUNRabs(ida_mem.ida_thi - ida_mem.ida_tlo) / ida_mem.ida_ttol;
            let fracsub = if fracint > FIVE { PT1 } else { HALF / fracint };
            tmid = ida_mem.ida_tlo + fracsub * (ida_mem.ida_thi - ida_mem.ida_tlo);
        }
        if SUNRabs(ida_mem.ida_thi - tmid) < HALF * ida_mem.ida_ttol {
            let fracint = SUNRabs(ida_mem.ida_thi - ida_mem.ida_tlo) / ida_mem.ida_ttol;
            let fracsub = if fracint > FIVE { PT1 } else { HALF / fracint };
            tmid = ida_mem.ida_thi - fracsub * (ida_mem.ida_thi - ida_mem.ida_tlo);
        }

        {
            let mut yy = std::mem::take(&mut ida_mem.ida_yy);
            let mut yp = std::mem::take(&mut ida_mem.ida_yp);
            let _ = IDAGetSolution(ida_mem, tmid, &mut yy, &mut yp);
            ida_mem.ida_yy = yy;
            ida_mem.ida_yp = yp;
        }
        let retval = {
            let gfun = ida_mem.ida_gfun.unwrap();
            let IDAMem { ida_yy, ida_yp, ida_grout, ida_user_data, .. } = ida_mem;
            gfun(tmid, ida_yy, ida_yp, ida_grout, ida_user_data)
        };
        ida_mem.ida_nge += 1;
        if retval != 0 {
            return IDA_RTFUNC_FAIL;
        }

        /* Check to see in which subinterval g changes sign, and reset imax.
           Set side = 1 if sign change is on low side, or 2 if on high side.  */
        maxfrac = ZERO;
        zroot = SUNFALSE;
        sgnchg = SUNFALSE;
        sideprev = side;
        for i in 0..nrtfn {
            if !ida_mem.ida_gactive[i] {
                continue;
            }
            if SUNRabs(ida_mem.ida_grout[i]) == ZERO {
                if ida_mem.ida_rootdir[i] as f64 * ida_mem.ida_glo[i] <= ZERO {
                    zroot = SUNTRUE;
                }
            } else {
                if SUNRdifferentsign(ida_mem.ida_glo[i], ida_mem.ida_grout[i])
                    && ida_mem.ida_rootdir[i] as f64 * ida_mem.ida_glo[i] <= ZERO
                {
                    let gfrac = SUNRabs(ida_mem.ida_grout[i]
                                        / (ida_mem.ida_grout[i] - ida_mem.ida_glo[i]));
                    if gfrac > maxfrac {
                        sgnchg = SUNTRUE;
                        maxfrac = gfrac;
                        imax = i;
                    }
                }
            }
        }
        if sgnchg {
            /* Sign change found in (tlo,tmid); replace thi with tmid. */
            ida_mem.ida_thi = tmid;
            for i in 0..nrtfn {
                ida_mem.ida_ghi[i] = ida_mem.ida_grout[i];
            }
            side = 1;
            /* Stop at root thi if converged; otherwise loop. */
            if SUNRabs(ida_mem.ida_thi - ida_mem.ida_tlo) <= ida_mem.ida_ttol {
                break;
            }
            continue; /* Return to looping point. */
        }

        if zroot {
            /* No sign change in (tlo,tmid), but g = 0 at tmid; return root tmid. */
            ida_mem.ida_thi = tmid;
            for i in 0..nrtfn {
                ida_mem.ida_ghi[i] = ida_mem.ida_grout[i];
            }
            break;
        }

        /* No sign change in (tlo,tmid), and no zero at tmid.
           Sign change must be in (tmid,thi).  Replace tlo with tmid. */
        ida_mem.ida_tlo = tmid;
        for i in 0..nrtfn {
            ida_mem.ida_glo[i] = ida_mem.ida_grout[i];
        }
        side = 2;
        /* Stop at root thi if converged; otherwise loop back. */
        if SUNRabs(ida_mem.ida_thi - ida_mem.ida_tlo) <= ida_mem.ida_ttol {
            break;
        }
    } /* End of root-search loop */

    /* Reset trout and grout, set iroots, and return RTFOUND. */
    ida_mem.ida_trout = ida_mem.ida_thi;
    for i in 0..nrtfn {
        ida_mem.ida_grout[i] = ida_mem.ida_ghi[i];
        ida_mem.ida_iroots[i] = 0;
        if !ida_mem.ida_gactive[i] {
            continue;
        }
        if SUNRabs(ida_mem.ida_ghi[i]) == ZERO
            && ida_mem.ida_rootdir[i] as f64 * ida_mem.ida_glo[i] <= ZERO
        {
            ida_mem.ida_iroots[i] = if ida_mem.ida_glo[i] > ZERO { -1 } else { 1 };
        }
        if SUNRdifferentsign(ida_mem.ida_glo[i], ida_mem.ida_ghi[i])
            && ida_mem.ida_rootdir[i] as f64 * ida_mem.ida_glo[i] <= ZERO
        {
            ida_mem.ida_iroots[i] = if ida_mem.ida_glo[i] > ZERO { -1 } else { 1 };
        }
    }
    RTFOUND
}

/*
 * =================================================================
 * Internal DQ approximations for sensitivity RHS
 * =================================================================
 */

/*
 * IDASensResDQ
 *
 * IDASensRhsDQ computes the residuals of the sensitivity equations
 * by finite differences. It is of type IDASensResFn.
 * Returns 0 if successful, <0 if an unrecoverable failure occurred,
 * >0 for a recoverable error.
 *
 * (C receives IDA_mem through the void* user_dataS slot; the Rust
 * signature moves it to the front, as with IDAQuadSensRhsInternalDQ.)
 */
pub fn IDASensResDQ(ida_mem: &mut IDAMem, Ns: i32, t: f64, yy: &NVector, yp: &NVector,
                    resval: &NVector, yyS: &[NVector], ypS: &[NVector],
                    resvalS: &mut [NVector], ytemp: &mut NVector, yptemp: &mut NVector,
                    restemp: &mut NVector) -> i32 {
    for is in 0..(Ns as usize) {
        let retval = IDASensRes1DQ(ida_mem, Ns, t, yy, yp, resval, is as i32, &yyS[is],
                                   &ypS[is], &mut resvalS[is], ytemp, yptemp, restemp);
        if retval != 0 {
            return retval;
        }
    }
    0
}

/*
 * IDASensRes1DQ
 *
 * IDASensRes1DQ computes the residual of the is-th sensitivity
 * equation by finite differences.
 *
 * Returns 0 if successful or the return value of res if res fails
 * (<0 if res fails unrecoverably, >0 if res has a recoverable error).
 *
 * (C locals Del/rDel/r2Del are delta/rdelta/r2delta here; Delp/Dely
 * families keep their names lowercased.  The scalar N_VLinearSum
 * serial kernels are: b==ONE → VLin1(a, x, y, z): z = a*x + y;
 * a==ONE → VLin1(b, y, x, z): z = b*y + x; a==-b → VScaleDiff(a, x,
 * y, z): z = a*(x - y); a==b==ONE → VSum.  Early error returns skip
 * the psave restore, exactly as in C.)
 */
fn IDASensRes1DQ(ida_mem: &mut IDAMem, _Ns: i32, t: f64, yy: &NVector, yp: &NVector,
                 resval: &NVector, is: i32, yyS: &NVector, ypS: &NVector,
                 resvalS: &mut NVector, ytemp: &mut NVector, yptemp: &mut NVector,
                 restemp: &mut NVector) -> i32 {
    /* Set base perturbation del */
    let del = SUNRsqrt(SUNMAX(ida_mem.ida_rtol, ida_mem.ida_uround));
    let rdel = ONE / del;

    let pbari = ida_mem.ida_pbar[is as usize];

    let which = ida_mem.ida_plist[is as usize] as usize;

    let psave = ida_mem.ida_p[which];

    let delp = pbari * del;
    let rdelp = ONE / delp;
    let norms = N_VWrmsNorm(yyS, &ida_mem.ida_ewt) * pbari;
    let rdely = SUNMAX(norms, rdel) / pbari;
    let dely = ONE / rdely;

    let method = if ida_mem.ida_DQrhomax == ZERO {
        /* No switching */
        if ida_mem.ida_DQtype == IDA_CENTERED { CENTERED1 } else { FORWARD1 }
    } else {
        /* switch between simultaneous/separate DQ */
        let ratio = dely * rdelp;
        if SUNMAX(ONE / ratio, ratio) <= ida_mem.ida_DQrhomax {
            if ida_mem.ida_DQtype == IDA_CENTERED { CENTERED1 } else { FORWARD1 }
        } else {
            if ida_mem.ida_DQtype == IDA_CENTERED { CENTERED2 } else { FORWARD2 }
        }
    };

    let res = ida_mem.ida_res.unwrap();

    match method {
        CENTERED1 => {
            let delta = SUNMIN(dely, delp);
            let r2delta = HALF / delta;

            /* Forward perturb y, y' and parameter */
            /* N_VLinearSum(Del, yyS, ONE, yy, ytemp) — VLin1(a, x, y): z = a*x + y */
            for (yt, (s, y)) in ytemp.data.iter_mut().zip(yyS.data.iter().zip(&yy.data)) {
                *yt = delta * *s + *y;
            }
            for (yt, (s, y)) in yptemp.data.iter_mut().zip(ypS.data.iter().zip(&yp.data)) {
                *yt = delta * *s + *y;
            }
            ida_mem.ida_p[which] = psave + delta;

            /* Save residual in resvalS */
            let retval = res(t, ytemp, yptemp, resvalS, &mut ida_mem.ida_user_data);
            ida_mem.ida_nreS += 1;
            if retval != 0 {
                return retval;
            }

            /* Backward perturb y, y' and parameter */
            for (yt, (s, y)) in ytemp.data.iter_mut().zip(yyS.data.iter().zip(&yy.data)) {
                *yt = -delta * *s + *y;
            }
            for (yt, (s, y)) in yptemp.data.iter_mut().zip(ypS.data.iter().zip(&yp.data)) {
                *yt = -delta * *s + *y;
            }
            ida_mem.ida_p[which] = psave - delta;

            /* Save residual in restemp */
            let retval = res(t, ytemp, yptemp, restemp, &mut ida_mem.ida_user_data);
            ida_mem.ida_nreS += 1;
            if retval != 0 {
                return retval;
            }

            /* Estimate the residual for the i-th sensitivity equation */
            /* N_VLinearSum(r2Del, resvalS, -r2Del, restemp, resvalS) —
               VScaleDiff: z = c*(x - y) */
            for (r, s) in resvalS.data.iter_mut().zip(&restemp.data) {
                *r = r2delta * (*r - *s);
            }
        }

        CENTERED2 => {
            let r2delp = HALF / delp;
            let r2dely = HALF / dely;

            /* Forward perturb y and y' */
            for (yt, (s, y)) in ytemp.data.iter_mut().zip(yyS.data.iter().zip(&yy.data)) {
                *yt = dely * *s + *y;
            }
            for (yt, (s, y)) in yptemp.data.iter_mut().zip(ypS.data.iter().zip(&yp.data)) {
                *yt = dely * *s + *y;
            }

            /* Save residual in resvalS */
            let retval = res(t, ytemp, yptemp, resvalS, &mut ida_mem.ida_user_data);
            ida_mem.ida_nreS += 1;
            if retval != 0 {
                return retval;
            }

            /* Backward perturb y and y' */
            for (yt, (s, y)) in ytemp.data.iter_mut().zip(yyS.data.iter().zip(&yy.data)) {
                *yt = -dely * *s + *y;
            }
            for (yt, (s, y)) in yptemp.data.iter_mut().zip(ypS.data.iter().zip(&yp.data)) {
                *yt = -dely * *s + *y;
            }

            /* Save residual in restemp */
            let retval = res(t, ytemp, yptemp, restemp, &mut ida_mem.ida_user_data);
            ida_mem.ida_nreS += 1;
            if retval != 0 {
                return retval;
            }

            /* Save the first difference quotient in resvalS */
            /* VScaleDiff */
            for (r, s) in resvalS.data.iter_mut().zip(&restemp.data) {
                *r = r2dely * (*r - *s);
            }

            /* Forward perturb parameter */
            ida_mem.ida_p[which] = psave + delp;

            /* Save residual in ytemp */
            let retval = res(t, yy, yp, ytemp, &mut ida_mem.ida_user_data);
            ida_mem.ida_nreS += 1;
            if retval != 0 {
                return retval;
            }

            /* Backward perturb parameter */
            ida_mem.ida_p[which] = psave - delp;

            /* Save residual in yptemp */
            let retval = res(t, yy, yp, yptemp, &mut ida_mem.ida_user_data);
            ida_mem.ida_nreS += 1;
            if retval != 0 {
                return retval;
            }

            /* Save the second difference quotient in restemp */
            /* N_VLinearSum(r2Delp, ytemp, -r2Delp, yptemp, restemp) — VScaleDiff */
            for (r, (a, b)) in restemp.data.iter_mut()
                .zip(ytemp.data.iter().zip(&yptemp.data))
            {
                *r = r2delp * (*a - *b);
            }

            /* Add the difference quotients for the sensitivity residual */
            /* N_VLinearSum(ONE, resvalS, ONE, restemp, resvalS) — VSum */
            for (r, s) in resvalS.data.iter_mut().zip(&restemp.data) {
                *r = *r + *s;
            }
        }

        FORWARD1 => {
            let delta = SUNMIN(dely, delp);
            let rdelta = ONE / delta;

            /* Forward perturb y, y' and parameter */
            for (yt, (s, y)) in ytemp.data.iter_mut().zip(yyS.data.iter().zip(&yy.data)) {
                *yt = delta * *s + *y;
            }
            for (yt, (s, y)) in yptemp.data.iter_mut().zip(ypS.data.iter().zip(&yp.data)) {
                *yt = delta * *s + *y;
            }
            ida_mem.ida_p[which] = psave + delta;

            /* Save residual in resvalS */
            let retval = res(t, ytemp, yptemp, resvalS, &mut ida_mem.ida_user_data);
            ida_mem.ida_nreS += 1;
            if retval != 0 {
                return retval;
            }

            /* Estimate the residual for the i-th sensitivity equation */
            /* N_VLinearSum(rDel, resvalS, -rDel, resval, resvalS) — VScaleDiff */
            for (r, v) in resvalS.data.iter_mut().zip(&resval.data) {
                *r = rdelta * (*r - *v);
            }
        }

        FORWARD2 => {
            /* Forward perturb y and y' */
            for (yt, (s, y)) in ytemp.data.iter_mut().zip(yyS.data.iter().zip(&yy.data)) {
                *yt = dely * *s + *y;
            }
            for (yt, (s, y)) in yptemp.data.iter_mut().zip(ypS.data.iter().zip(&yp.data)) {
                *yt = dely * *s + *y;
            }

            /* Save residual in resvalS */
            let retval = res(t, ytemp, yptemp, resvalS, &mut ida_mem.ida_user_data);
            ida_mem.ida_nreS += 1;
            if retval != 0 {
                return retval;
            }

            /* Save the first difference quotient in resvalS */
            /* VScaleDiff */
            for (r, v) in resvalS.data.iter_mut().zip(&resval.data) {
                *r = rdely * (*r - *v);
            }

            /* Forward perturb parameter */
            ida_mem.ida_p[which] = psave + delp;

            /* Save residual in restemp */
            let retval = res(t, yy, yp, restemp, &mut ida_mem.ida_user_data);
            ida_mem.ida_nreS += 1;
            if retval != 0 {
                return retval;
            }

            /* Save the second difference quotient in restemp */
            /* N_VLinearSum(rDelp, restemp, -rDelp, resval, restemp) — VScaleDiff */
            for (r, v) in restemp.data.iter_mut().zip(&resval.data) {
                *r = rdelp * (*r - *v);
            }

            /* Add the difference quotients for the sensitivity residual */
            /* VSum */
            for (r, s) in resvalS.data.iter_mut().zip(&restemp.data) {
                *r = *r + *s;
            }
        }

        _ => {}
    }

    /* Restore original value of parameter */
    ida_mem.ida_p[which] = psave;

    0
}

/* IDAQuadSensRhsInternalDQ   - internal IDAQuadSensRhsFn
 *
 * IDAQuadSensRhsInternalDQ computes right hand side of all quadrature
 * sensitivity equations by finite differences. All work is actually
 * done in IDAQuadSensRhs1InternalDQ.
 *
 * (C receives IDA_mem through the void* user-data slot; the pinned
 * Rust signature moves it to the front.)
 */
fn IDAQuadSensRhsInternalDQ(ida_mem: &mut IDAMem, Ns: i32, t: f64, yy: &NVector,
                            yp: &NVector, yyS: &[NVector], ypS: &[NVector], rrQ: &NVector,
                            resvalQS: &mut [NVector], yytmp: &mut NVector,
                            yptmp: &mut NVector, tmpQS: &mut NVector) -> i32 {
    for is in 0..(Ns as usize) {
        let retval = IDAQuadSensRhs1InternalDQ(ida_mem, is as i32, t, yy, yp, &yyS[is],
                                               &ypS[is], rrQ, &mut resvalQS[is], yytmp,
                                               yptmp, tmpQS);
        if retval != 0 {
            return retval;
        }
    }

    0
}

/*
 * IDAQuadSensRhs1InternalDQ
 *
 * (Early error returns skip both the psave restore and the nrQeS
 * increment, exactly as in C.  Note C's FORWARD1 arm reuses and
 * overwrites `rdel = ONE / Del`.)
 */
fn IDAQuadSensRhs1InternalDQ(ida_mem: &mut IDAMem, is: i32, t: f64, yy: &NVector,
                             yp: &NVector, yyS: &NVector, ypS: &NVector,
                             resvalQ: &NVector, resvalQS: &mut NVector,
                             yytmp: &mut NVector, yptmp: &mut NVector,
                             tmpQS: &mut NVector) -> i32 {
    let mut nfel: i64 = 0;

    let del = SUNRsqrt(SUNMAX(ida_mem.ida_rtol, ida_mem.ida_uround));
    let mut rdel = ONE / del;

    let pbari = ida_mem.ida_pbar[is as usize];

    let which = ida_mem.ida_plist[is as usize] as usize;

    let psave = ida_mem.ida_p[which];

    let delp = pbari * del;
    let norms = N_VWrmsNorm(yyS, &ida_mem.ida_ewt) * pbari;
    let rdely = SUNMAX(norms, rdel) / pbari;
    let dely = ONE / rdely;

    let method = if ida_mem.ida_DQtype == IDA_CENTERED { CENTERED1 } else { FORWARD1 };

    let rhsQ = ida_mem.ida_rhsQ.unwrap();

    match method {
        CENTERED1 => {
            let delta = SUNMIN(dely, delp);
            let r2delta = HALF / delta;

            /* N_VLinearSum(ONE, yy, Del, yyS, yytmp) — VLin1(b, y, x): z = b*y + x */
            for (yt, (y, s)) in yytmp.data.iter_mut().zip(yy.data.iter().zip(&yyS.data)) {
                *yt = delta * *s + *y;
            }
            for (yt, (y, s)) in yptmp.data.iter_mut().zip(yp.data.iter().zip(&ypS.data)) {
                *yt = delta * *s + *y;
            }
            ida_mem.ida_p[which] = psave + delta;

            let retval = rhsQ(t, yytmp, yptmp, resvalQS, &mut ida_mem.ida_user_data);
            nfel += 1;
            if retval != 0 {
                return retval;
            }

            /* N_VLinearSum(-Del, yyS, ONE, yy, yytmp) — VLin1(a, x, y): z = a*x + y */
            for (yt, (s, y)) in yytmp.data.iter_mut().zip(yyS.data.iter().zip(&yy.data)) {
                *yt = -delta * *s + *y;
            }
            for (yt, (s, y)) in yptmp.data.iter_mut().zip(ypS.data.iter().zip(&yp.data)) {
                *yt = -delta * *s + *y;
            }

            ida_mem.ida_p[which] = psave - delta;

            let retval = rhsQ(t, yytmp, yptmp, tmpQS, &mut ida_mem.ida_user_data);
            nfel += 1;
            if retval != 0 {
                return retval;
            }

            /* N_VLinearSum(r2Del, resvalQS, -r2Del, tmpQS, resvalQS) —
               VScaleDiff: z = c*(x - y) */
            for (r, s) in resvalQS.data.iter_mut().zip(&tmpQS.data) {
                *r = r2delta * (*r - *s);
            }
        }

        FORWARD1 => {
            let delta = SUNMIN(dely, delp);
            /* (C reuses/overwrites rdel here) */
            rdel = ONE / delta;

            /* N_VLinearSum(ONE, yy, Del, yyS, yytmp) — VLin1(b, y, x): z = b*y + x */
            for (yt, (y, s)) in yytmp.data.iter_mut().zip(yy.data.iter().zip(&yyS.data)) {
                *yt = delta * *s + *y;
            }
            for (yt, (y, s)) in yptmp.data.iter_mut().zip(yp.data.iter().zip(&ypS.data)) {
                *yt = delta * *s + *y;
            }
            ida_mem.ida_p[which] = psave + delta;

            let retval = rhsQ(t, yytmp, yptmp, resvalQS, &mut ida_mem.ida_user_data);
            nfel += 1;
            if retval != 0 {
                return retval;
            }

            /* N_VLinearSum(rdel, resvalQS, -rdel, resvalQ, resvalQS) — VScaleDiff */
            for (r, q) in resvalQS.data.iter_mut().zip(&resvalQ.data) {
                *r = rdel * (*r - *q);
            }
        }

        _ => {}
    }

    ida_mem.ida_p[which] = psave;
    /* Increment counter nrQeS */
    ida_mem.ida_nrQeS += nfel;

    0
}

/* (IDAProcessError, which closes idas.c, lives in idas_impl.rs.
   END of idas.c port.) */









