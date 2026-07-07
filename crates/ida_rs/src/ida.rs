/* -----------------------------------------------------------------
 * Translated from src/ida/ida.c (IDA 7.7.0) — PART 1 (ida.c:1-2074).
 * Main IDA integrator for DAE systems F(t, y, y') = 0:
 * creation/initialization/re-initialization, tolerances, rootfinding
 * init, the IDASolve driver, interpolated output (IDAGetDky,
 * IDAComputeY/Yp), deallocation, vector allocation helpers and
 * IDAInitialSetup.  The remaining internals (IDAStep and below) are
 * PART 2, appended after the trailer at the bottom of this file.
 * Conventions follow the donor cvode_rs/src/cvode.rs.
 * -----------------------------------------------------------------*/
use crate::ida_impl::*;
use crate::ida_ls::{idaLsInitialize, idaLsInitializeCounters, idaLsPerf};
use crate::ida_nls::{idaNlsInit, IDASetNonlinearSolver};
use crate::nvector_serial::*;
use crate::sundials_context::SUNContext;
use crate::sundials_math::SUNRabs;
use crate::sundials_types::*;
use crate::sundials_utils::fmt_g;
use crate::sunnonlinsol_newton::SUNNonlinSol_Newton;

/*
 * =================================================================
 * IDA PRIVATE CONSTANTS
 * =================================================================
 */

const ZERO: f64 = 0.0; /* real 0.0    */
const HALF: f64 = 0.5; /* real 0.5    */
/* TWOTHIRDS (0.667) is only used by IDACreate's steptol default,
   which lives in the ida_impl::IDAMem Default impl. */
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
 * IDA ROUTINE-SPECIFIC CONSTANTS
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

/* (The itol control constants IDA_NN/IDA_SS/IDA_SV/IDA_WF and the
   algorithmic constants MXNCF/MXNEF/MAXNH/MAXNJ/MAXNI/EPCON/MAXBACKS
   are defined in ida_impl.rs, where the IDAMem Default impl uses
   them.) */

/*
 * -----------------------------------------------------------------
 * Message rendering helper
 *
 * The C IDAProcessError is printf-style; the ida_impl.rs message
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
 * default-setting body of IDACreate is the ida_impl::IDAMem Default
 * impl, mirroring how the donor's CVodeCreate builds its literal.)
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
       Rust, so the MSG_MEM_FAIL branch of the C code vanishes */

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
    /* (C: ida_edata will be set to user_data in InitialSetup; the
       Rust dispatch passes ida_user_data directly.) */

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

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
 * ida_lsetup/ida_lsolve, whose only call sites are in ida_nls.c and
 * ida_ic.c and whose dispatch therefore lives there); the module is
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
   IDALsMem.iterative (see ida_impl.rs) */
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
           check for approach to tstop, and scale phi[1] by hh.
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
            let ypnorm = IDAWrmsNorm(ida_mem, &ida_mem.ida_phi[1], &ida_mem.ida_ewt,
                                     ida_mem.ida_suppressalg);
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

        /* Reset and check ewt (if not first call). */

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
        }

        /* Check for too much accuracy requested. */

        let nrm = IDAWrmsNorm(ida_mem, &ida_mem.ida_phi[0], &ida_mem.ida_ewt,
                              ida_mem.ida_suppressalg);
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
 * -----------------------------------------------------------------
 * Deallocation function
 * -----------------------------------------------------------------
 */

/*
 * IDAFree
 *
 * This routine frees the problem memory allocated by IDAInit.
 * In C this releases the vectors (IDAFreeVectors), the owned NLS
 * object, the linear solver memory (lfree) and the rootfinding
 * arrays; here all of that is RAII — dropping the Box releases
 * everything.
 */
pub fn IDAFree(_ida_mem: Box<IDAMem>) {}

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
 * Adaptation: ida_yy and ida_yp — which in C are never allocated
 * because they alias the user's yret/ypret during IDASolve (and
 * IDACalcIC's yy0/yp0) — are allocated here as owned stand-ins
 * (donor cv_y pattern).  They are deliberately NOT counted in
 * lrw/liw, matching the C workspace accounting.
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

/* IDAFreeVectors has no Rust counterpart: the vectors are owned
   fields released by drop, and the lrw/liw bookkeeping it performs
   is only observable through IDAGetWorkSpace before the memory is
   destroyed. */

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
 * (Non-static in C: also called by IDACalcIC in ida_ic.c.)
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

    /* Check to see if y0 satisfies constraints. */
    if ida_mem.ida_constraintsSet {
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
    let ier = idaNlsInit(ida_mem);
    if ier != IDA_SUCCESS {
        IDAProcessError(Some(ida_mem), IDA_NLS_INIT_FAIL, line!(), "IDAInitialSetup", file!(),
                        MSG_NLS_INIT_FAIL);
        return IDA_NLS_INIT_FAIL;
    }

    IDA_SUCCESS
}

// ===================== END PART 1 (ida.c:1-2074) =====================
// PART 2 (remaining internals of ida.c:2076-4101, appended below by the
// next agent):
//   IDAEwtSet / IDAEwtSetSS / IDAEwtSetSV        (ewt setters)
//   IDAStopTest1 / IDAStopTest2                  (stopping tests)
//   IDAHandleFailure                             (failure -> istate map)
//   IDAStep / IDASetCoeffs                       (main step driver)
//   IDANls / IDAPredict / IDACheckConstraints    (nonlinear drivers)
//   IDATestError / IDARestore / IDAHandleNFlag / IDAReset  (error tests)
//   IDACompleteStep                              (order/step selection)
//   IDAGetSolution                               (interpolation)
//   IDAWrmsNorm                                  (norms)
//   IDARcheck1 / IDARcheck2 / IDARcheck3 / IDARootfind     (rootfinding)
// (IDAProcessError/IDAErrHandler live in ida_impl.rs; the alloc/free
// helpers IDACheckNvector/IDAAllocVectors/IDAFreeVectors are above in
// PART 1.)
//
// Signatures Part 1 compiles against (must match exactly):
//   pub fn IDAEwtSet(ida_mem: &IDAMem, ycur: &NVector, weight: &mut NVector) -> i32
//   pub fn IDAWrmsNorm(ida_mem: &IDAMem, x: &NVector, w: &NVector, mask: bool) -> f64
//   pub fn IDAGetSolution(ida_mem: &IDAMem, t: f64, yret: &mut NVector,
//                         ypret: &mut NVector) -> i32
//   fn IDAStopTest1(ida_mem: &mut IDAMem, tout: f64, tret: &mut f64,
//                   yret: &mut NVector, ypret: &mut NVector, itask: i32) -> i32
//   fn IDAStopTest2(ida_mem: &mut IDAMem, tout: f64, tret: &mut f64,
//                   yret: &mut NVector, ypret: &mut NVector, itask: i32) -> i32
//     (aliasing: in the ONE_STEP/tstop branches where C relies on
//      yret==ida_yy, copy ida_yy -> yret and ida_yp -> ypret)
//   fn IDAHandleFailure(ida_mem: &mut IDAMem, sflag: i32) -> i32
//   fn IDAStep(ida_mem: &mut IDAMem) -> i32
//   fn IDARcheck1(ida_mem: &mut IDAMem) -> i32
//   fn IDARcheck2(ida_mem: &mut IDAMem) -> i32
//   fn IDARcheck3(ida_mem: &mut IDAMem, tout: f64, itask: i32) -> i32
//     (the Rcheck routines write the interpolated solution into the
//      owned ida_yy/ida_yp; IDASolve copies them back to yret/ypret)

/* -----------------------------------------------------------------
 * PART 2 (ida.c:2076-4101), appended below.
 * (IDANls itself lives in ida_nls.rs together with the collapsed
 * SUNNonlinSolSolve_Newton loop; IDAStep dispatches into it.
 * IDAProcessError lives in ida_impl.rs.)
 * -----------------------------------------------------------------*/
use crate::ida_nls::IDANls;
use crate::sundials_errors::SUN_ERR_ARG_CORRUPT;
use crate::sundials_math::{SUNMAX, SUNMIN, SUNRdifferentsign, SUNRpowerR};

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
 *
 * (In C the signature is IDAEwtSet(ycur, weight, data) with data
 * pointing to IDA_mem; here the memory block is the first argument,
 * donor cvEwtSet convention.)
 */
pub fn IDAEwtSet(ida_mem: &IDAMem, ycur: &NVector, weight: &mut NVector) -> i32 {
    match ida_mem.ida_itol {
        IDA_SS => IDAEwtSetSS(ida_mem, ycur, weight),
        IDA_SV => IDAEwtSetSV(ida_mem, ycur, weight),
        _ => 0,
    }
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
 * (The C version accumulates in ida_tempv1 and inverts into weight;
 * computing directly in `weight` performs the identical arithmetic —
 * donor cvEwtSetSS adaptation.)
 */
fn IDAEwtSetSS(ida_mem: &IDAMem, ycur: &NVector, weight: &mut NVector) -> i32 {
    N_VAbs(ycur, weight);
    weight.scale_inplace(ida_mem.ida_rtol);
    weight.add_const_inplace(ida_mem.ida_Satol);
    if ida_mem.ida_atolmin0 && N_VMin(weight) <= ZERO {
        return -1;
    }
    weight.invert_inplace();
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
    N_VAbs(ycur, weight);
    weight.linear_sum_with(ida_mem.ida_rtol, ONE, &ida_mem.ida_Vatol);
    if ida_mem.ida_atolmin0 && N_VMin(weight) <= ZERO {
        return -1;
    }
    weight.invert_inplace();
    0
}

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
fn IDAStopTest1(
    ida_mem: &mut IDAMem,
    tout: f64,
    tret: &mut f64,
    yret: &mut NVector,
    ypret: &mut NVector,
    itask: i32,
) -> i32 {
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
                    IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAStopTest1", file!(),
                                    &ida_msg_g(MSG_BAD_TSTOP,
                                               &[ida_mem.ida_tstop, ida_mem.ida_tn]));
                    return IDA_ILL_INPUT;
                }
                ida_mem.ida_tretlast = ida_mem.ida_tstop;
                *tret = ida_mem.ida_tstop;
                ida_mem.ida_tstopset = SUNFALSE;
                return IDA_TSTOP_RETURN;
            }
        }
        /* Test for tn approaching tstop */
        else if (ida_mem.ida_tn + ida_mem.ida_hh - ida_mem.ida_tstop) * ida_mem.ida_hh > ZERO {
            ida_mem.ida_hh =
                (ida_mem.ida_tstop - ida_mem.ida_tn) * (ONE - FOUR * ida_mem.ida_uround);
        }
    }

    match itask {
        IDA_NORMAL => {
            /* Test for tout = tretlast, and for tn past tout. */
            if tout == ida_mem.ida_tretlast {
                ida_mem.ida_tretlast = tout;
                *tret = tout;
                return IDA_SUCCESS;
            }
            if (ida_mem.ida_tn - tout) * ida_mem.ida_hh >= ZERO {
                let ier = IDAGetSolution(ida_mem, tout, yret, ypret);
                if ier != IDA_SUCCESS {
                    IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAStopTest1", file!(),
                                    &ida_msg_g(MSG_BAD_TOUT, &[tout]));
                    return IDA_ILL_INPUT;
                }
                ida_mem.ida_tretlast = tout;
                *tret = tout;
                return IDA_SUCCESS;
            }

            CONTINUE_STEPS
        }

        IDA_ONE_STEP => {
            /* Test for tn past tretlast. */
            if (ida_mem.ida_tn - ida_mem.ida_tretlast) * ida_mem.ida_hh > ZERO {
                let _ = IDAGetSolution(ida_mem, ida_mem.ida_tn, yret, ypret);
                ida_mem.ida_tretlast = ida_mem.ida_tn;
                *tret = ida_mem.ida_tn;
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
 * because yret and ypret already contain the current y and y' values
 * (in C through the ida_yy/ida_yp alias; here the owned vectors are
 * copied back — workspace rule 5).
 *
 * Note: No test is made for an error return from IDAGetSolution here,
 * because the same test was made prior to the step.
 */
fn IDAStopTest2(
    ida_mem: &mut IDAMem,
    tout: f64,
    tret: &mut f64,
    yret: &mut NVector,
    ypret: &mut NVector,
    itask: i32,
) -> i32 {
    if ida_mem.ida_tstopset {
        let troundoff = HUNDRED * ida_mem.ida_uround *
            (SUNRabs(ida_mem.ida_tn) + SUNRabs(ida_mem.ida_hh));

        /* Test for tn at tstop */
        if SUNRabs(ida_mem.ida_tn - ida_mem.ida_tstop) <= troundoff {
            /* Ensure tout >= tstop, otherwise check for tout return below */
            if (tout - ida_mem.ida_tstop) * ida_mem.ida_hh >= ZERO
                || SUNRabs(tout - ida_mem.ida_tstop) <= troundoff
            {
                let _ = IDAGetSolution(ida_mem, ida_mem.ida_tstop, yret, ypret);
                ida_mem.ida_tretlast = ida_mem.ida_tstop;
                *tret = ida_mem.ida_tstop;
                ida_mem.ida_tstopset = SUNFALSE;
                return IDA_TSTOP_RETURN;
            }
        }
        /* Test for tn approaching tstop */
        else if (ida_mem.ida_tn + ida_mem.ida_hh - ida_mem.ida_tstop) * ida_mem.ida_hh > ZERO {
            ida_mem.ida_hh =
                (ida_mem.ida_tstop - ida_mem.ida_tn) * (ONE - FOUR * ida_mem.ida_uround);
        }
    }

    match itask {
        IDA_NORMAL => {
            /* Test for tn past tout. */
            if (ida_mem.ida_tn - tout) * ida_mem.ida_hh >= ZERO {
                let _ = IDAGetSolution(ida_mem, tout, yret, ypret);
                ida_mem.ida_tretlast = tout;
                *tret = tout;
                return IDA_SUCCESS;
            }

            CONTINUE_STEPS
        }

        IDA_ONE_STEP => {
            ida_mem.ida_tretlast = ida_mem.ida_tn;
            *tret = ida_mem.ida_tn;
            /* C relies on yret == ida_yy / ypret == ida_yp here; copy
               the owned vectors back (Part-1 trailer contract) */
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
            IDAProcessError(Some(ida_mem), IDA_NLS_SETUP_FAIL, line!(), "IDAHandleFailure", file!(),
                            &ida_msg_g(MSG_NLS_SETUP_FAILED, &[ida_mem.ida_tn]));
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
                            "IDA encountered an unrecognized error. Please report this to the \
                             Sundials developers at sundials-users@llnl.gov");
            IDA_UNRECOGNIZED_ERROR
        }
    }
}

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
 * Variables or arrays (all in the IDAMem structure) used in IDAStep are:
 *
 * tt -- Independent variable.
 * yy -- Solution vector at tt.
 * yp -- Derivative of solution vector after successful stelp.
 * res -- User-supplied function to evaluate the residual. See the
 *        description given in file ida.h .
 * lsetup -- Routine to prepare for the linear solver call. It may either
 *        save or recalculate quantities used by lsolve. (Optional)
 * lsolve -- Routine to solve a linear system. A prior call to lsetup
 *        may be required.
 * hh  -- Appropriate step size for next step.
 * ewt -- Vector of weights used in all convergence tests.
 * phi -- Array of divided differences used by IDAStep. This array is composed
 *       of  (maxord+1) nvectors (each of size Neq). (maxord+1) is the maximum
 *       order for the problem, maxord, plus 1.
 *
 *       Return values are:
 *       IDA_SUCCESS   IDA_RES_FAIL      LSETUP_ERROR_NONRECVR
 *                     IDA_LSOLVE_FAIL   IDA_ERR_FAIL
 *                     IDA_CONSTR_FAIL   IDA_CONV_FAIL
 *                     IDA_REP_RES_ERR
 */
fn IDAStep(ida_mem: &mut IDAMem) -> i32 {
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

    /* To prevent 'unintialized variable' warnings */
    let mut err_k = ZERO;
    let mut err_km1 = ZERO;
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
        if ida_mem.ida_tstopset
            && (ida_mem.ida_tn - ida_mem.ida_tstop) * ida_mem.ida_hh > ZERO
        {
            ida_mem.ida_tn = ida_mem.ida_tstop;
        }

        /*-----------------------
          Advance state variables
          -----------------------*/

        /* Compute predicted values for yy and yp */
        IDAPredict(ida_mem);

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
            nflag = IDATestError(ida_mem, ck, &mut err_k, &mut err_km1);
        }

        /* Test for convergence or error test failures */
        if nflag != IDA_SUCCESS {
            /* restore and decide what to do */
            IDARestore(ida_mem, saved_t);
            let kflag = IDAHandleNFlag(ida_mem, nflag, err_k, err_km1, &mut ncf, &mut nef);

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
 */
fn IDASetCoeffs(ida_mem: &mut IDAMem, ck: &mut f64) {
    /* Set coefficients for the current stepsize h */

    if ida_mem.ida_hh != ida_mem.ida_hused || ida_mem.ida_kk != ida_mem.ida_kused {
        ida_mem.ida_ns = 0;
    }
    ida_mem.ida_ns = std::cmp::min(ida_mem.ida_ns + 1, ida_mem.ida_kused + 2);
    if ida_mem.ida_kk + 1 >= ida_mem.ida_ns {
        ida_mem.ida_beta[0] = ONE;
        ida_mem.ida_alpha[0] = ONE;
        let mut temp1 = ida_mem.ida_hh;
        ida_mem.ida_gamma[0] = ZERO;
        ida_mem.ida_sigma[0] = ONE;
        for i in 1..=(ida_mem.ida_kk as usize) {
            let temp2 = ida_mem.ida_psi[i - 1];
            ida_mem.ida_psi[i - 1] = temp1;
            ida_mem.ida_beta[i] = ida_mem.ida_beta[i - 1] * ida_mem.ida_psi[i - 1] / temp2;
            temp1 = temp2 + ida_mem.ida_hh;
            ida_mem.ida_alpha[i] = ida_mem.ida_hh / temp1;
            ida_mem.ida_sigma[i] = i as f64 * ida_mem.ida_sigma[i - 1] * ida_mem.ida_alpha[i];
            ida_mem.ida_gamma[i] =
                ida_mem.ida_gamma[i - 1] + ida_mem.ida_alpha[i - 1] / ida_mem.ida_hh;
        }
        ida_mem.ida_psi[ida_mem.ida_kk as usize] = temp1;
    }
    /* compute alphas, alpha0 */
    let mut alphas = ZERO;
    let mut alpha0 = ZERO;
    for i in 0..(ida_mem.ida_kk as usize) {
        alphas -= ONE / (i + 1) as f64;
        alpha0 -= ida_mem.ida_alpha[i];
    }

    /* compute leading coefficient cj  */
    ida_mem.ida_cjlast = ida_mem.ida_cj;
    ida_mem.ida_cj = -alphas / ida_mem.ida_hh;

    /* compute variable stepsize error coefficient ck */

    *ck = SUNRabs(ida_mem.ida_alpha[ida_mem.ida_kk as usize] + alphas - alpha0);
    *ck = SUNMAX(*ck, ida_mem.ida_alpha[ida_mem.ida_kk as usize]);

    /* change phi to phi-star  */

    /* Scale i=IDA_mem->ida_ns to i<=IDA_mem->ida_kk
       (C: fused N_VScaleVectorArray(kk-ns+1, beta+ns, phi+ns, phi+ns);
       the serial kernel scales each phi[i] in place by beta[i]) */
    if ida_mem.ida_ns <= ida_mem.ida_kk {
        for i in (ida_mem.ida_ns as usize)..=(ida_mem.ida_kk as usize) {
            let b = ida_mem.ida_beta[i];
            ida_mem.ida_phi[i].scale_inplace(b);
        }
    }
}

/*
 * -----------------------------------------------------------------
 * Nonlinear solver functions
 * -----------------------------------------------------------------
 */

/* (IDANls is translated in ida_nls.rs, collapsed with the Newton
   solve loop of sunnonlinsol_newton.c that the C code reaches
   through the SUNNonlinearSolver ops table.) */

/*
 * IDACheckConstraints
 *
 * Check and enforce the inequality constraints on the corrected
 * solution.  mm = ida_tempv2 (mask), tmp = ida_tempv1 (workspace).
 */
fn IDACheckConstraints(
    ida_mem: &mut IDAMem,
    saved_t: f64,
    step_constraint_fails: &mut i32,
) -> i32 {
    /* Get mask vector mm, 1 where constraints failed and 0 otherwise */
    let constraintsPassed = {
        let IDAMem { ida_constraints, ida_yy, ida_tempv2, .. } = ida_mem;
        N_VConstrMask(ida_constraints, ida_yy, ida_tempv2)
    };
    if constraintsPassed {
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
        let IDAMem {
            ida_tempv1: tmp,
            ida_tempv2: mm,
            ida_tempv3,
            ida_constraints,
            ida_ewt,
            ida_yy,
            ..
        } = ida_mem;
        N_VCompare(ONEPT5, ida_constraints, tmp);
        tmp.prod_with(ida_constraints); /* tmp = tmp * constraints */
        tmp.div_with(ida_ewt); /* tmp = tmp / ewt */
        N_VScale(-PT1, tmp, ida_tempv3); /* tempv3 = -0.1*tmp (saved adjustment) */
        tmp.linear_sum_with(-PT1, ONE, ida_yy); /* tmp = yy - 0.1*tmp */
        tmp.prod_with(mm); /* tmp = tmp * mm */
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
        {
            let IDAMem {
                ida_tempv1: tmp,
                ida_tempv2: mm,
                ida_tempv3,
                ida_ee,
                ida_yypredict,
                ..
            } = ida_mem;

            /* Zero out the correction where any constraint failed */
            N_VProd(mm, ida_ee, tmp);
            ida_ee.linear_sum_with(ONE, -ONE, tmp);

            /* Set correction to zero out the predictor where any constraint failed */
            N_VProd(mm, ida_yypredict, tmp);
            ida_ee.linear_sum_with(ONE, -ONE, tmp);

            /* Update the correction where constraints failed and are strictly greater
               or less than zero to shift the state with the adjustment saved above */
            ida_tempv3.prod_with(mm);
            ida_ee.linear_sum_with(ONE, -ONE, ida_tempv3);
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
        let IDAMem { ida_tempv1: tmp, ida_tempv2: mm, ida_phi, ida_yy, .. } = ida_mem;
        N_VLinearSum(ONE, &ida_phi[0], -ONE, ida_yy, tmp);
        tmp.prod_with(mm);
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
 */
fn IDAPredict(ida_mem: &mut IDAMem) {
    let kk = ida_mem.ida_kk as usize;

    /* (C loads ida_cvals[0..kk] = ONE and calls the fused
       N_VLinearCombination(kk+1, cvals, phi, yypredict); the serial
       kernel computes z = c0*X0 then accumulates z += cj*Xj,
       replicated here) */
    let IDAMem { ida_phi, ida_gamma, ida_yypredict, ida_yppredict, .. } = ida_mem;

    for (z, p) in ida_yypredict.data.iter_mut().zip(&ida_phi[0].data) {
        *z = ONE * *p;
    }
    for j in 1..=kk {
        for (z, p) in ida_yypredict.data.iter_mut().zip(&ida_phi[j].data) {
            *z += ONE * *p;
        }
    }

    /* N_VLinearCombination(kk, gamma+1, phi+1, yppredict) */
    for (z, p) in ida_yppredict.data.iter_mut().zip(&ida_phi[1].data) {
        *z = ida_gamma[1] * *p;
    }
    for j in 2..=kk {
        for (z, p) in ida_yppredict.data.iter_mut().zip(&ida_phi[j].data) {
            *z += ida_gamma[j] * *p;
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
fn IDATestError(ida_mem: &mut IDAMem, ck: f64, err_k: &mut f64, err_km1: &mut f64) -> i32 {
    /* Compute error for order k. */
    let enorm_k = IDAWrmsNorm(ida_mem, &ida_mem.ida_ee, &ida_mem.ida_ewt,
                              ida_mem.ida_suppressalg);
    *err_k = ida_mem.ida_sigma[ida_mem.ida_kk as usize] * enorm_k;
    let terr_k = (ida_mem.ida_kk + 1) as f64 * *err_k;

    ida_mem.ida_knew = ida_mem.ida_kk;

    if ida_mem.ida_kk > 1 {
        /* Compute error at order k-1 */
        {
            let kk = ida_mem.ida_kk as usize;
            let IDAMem { ida_phi, ida_ee, ida_delta, .. } = ida_mem;
            N_VLinearSum(ONE, &ida_phi[kk], ONE, ida_ee, ida_delta);
        }
        let enorm_km1 = IDAWrmsNorm(ida_mem, &ida_mem.ida_delta, &ida_mem.ida_ewt,
                                    ida_mem.ida_suppressalg);
        *err_km1 = ida_mem.ida_sigma[(ida_mem.ida_kk - 1) as usize] * enorm_km1;
        let terr_km1 = ida_mem.ida_kk as f64 * *err_km1;

        if ida_mem.ida_kk > 2 {
            /* Compute error at order k-2 */
            {
                let kk = ida_mem.ida_kk as usize;
                let IDAMem { ida_phi, ida_delta, .. } = ida_mem;
                /* N_VLinearSum(ONE, phi[kk-1], ONE, delta, delta) */
                ida_delta.linear_sum_with(ONE, ONE, &ida_phi[kk - 1]);
            }
            let enorm_km2 = IDAWrmsNorm(ida_mem, &ida_mem.ida_delta, &ida_mem.ida_ewt,
                                        ida_mem.ida_suppressalg);
            let err_km2 = ida_mem.ida_sigma[(ida_mem.ida_kk - 2) as usize] * enorm_km2;
            let terr_km2 = (ida_mem.ida_kk - 1) as f64 * err_km2;

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
 * IDARestore
 *
 * This routine restores tn, psi, and phi in the event of a failure.
 * It changes back phi-star to phi (changed in IDASetCoeffs)
 */
fn IDARestore(ida_mem: &mut IDAMem, saved_t: f64) {
    ida_mem.ida_tn = saved_t;

    for j in 1..=(ida_mem.ida_kk as usize) {
        ida_mem.ida_psi[j - 1] = ida_mem.ida_psi[j] - ida_mem.ida_hh;
    }

    if ida_mem.ida_ns <= ida_mem.ida_kk {
        /* (C: cvals[j-ns] = ONE/beta[j], then the fused
           N_VScaleVectorArray scales each phi[j] in place) */
        for j in (ida_mem.ida_ns as usize)..=(ida_mem.ida_kk as usize) {
            let c = ONE / ida_mem.ida_beta[j];
            ida_mem.ida_phi[j].scale_inplace(c);
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
 *   SUN_NLS_CONV_RECVR         > 0
 *   IDA_RES_FAIL               < 0
 *   IDA_LSOLVE_FAIL            < 0
 *   IDA_LSETUP_FAIL            < 0
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
 *
 * (Adaptation: the C signature also passes pointers to the global
 * counters ida_ncfn/ida_netf, which alias IDA_mem fields; here they
 * are incremented directly on ida_mem and only the per-step counters
 * ncf/nef are passed by reference.)
 */
fn IDAHandleNFlag(
    ida_mem: &mut IDAMem,
    nflag: i32,
    err_k: f64,
    err_km1: f64,
    ncfPtr: &mut i32,
    nefPtr: &mut i32,
) -> i32 {
    ida_mem.ida_phase = 1;

    if nflag != ERROR_TEST_FAIL {
        /*-----------------------
          Nonlinear solver failed
          -----------------------*/

        *ncfPtr += 1; /* local counter for convergence failures */
        ida_mem.ida_ncfn += 1; /* global counter for convergence failures */

        if nflag < 0 {
            /* nonrecoverable failure */

            if nflag == IDA_LSOLVE_FAIL {
                IDA_LSOLVE_FAIL
            } else if nflag == IDA_LSETUP_FAIL {
                IDA_LSETUP_FAIL
            } else if nflag == IDA_RES_FAIL {
                IDA_RES_FAIL
            } else {
                IDA_NLS_FAIL
            }
        } else {
            /* recoverable failure    */

            /* Test if there were too many convergence failures or |h| = hmin */
            if *ncfPtr == ida_mem.ida_maxncf
                || SUNRabs(ida_mem.ida_hh) <= ida_mem.ida_hmin * ONEPSM
            {
                if nflag == IDA_RES_RECVR {
                    return IDA_REP_RES_ERR;
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

        *nefPtr += 1; /* local counter for error test failures */
        ida_mem.ida_netf += 1; /* global counter for error test failures */

        if *nefPtr == 1 {
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
        } else if *nefPtr == 2 {
            /* On second error test failure, use current order or decrease order by one.
               Reduce stepsize by factor of 1/4. */

            ida_mem.ida_kk = ida_mem.ida_knew;
            ida_mem.ida_eta = SUNMAX(ida_mem.ida_eta_min_ef,
                                     ida_mem.ida_hmin / SUNRabs(ida_mem.ida_hh));
            ida_mem.ida_hh *= ida_mem.ida_eta;

            PREDICT_AGAIN
        } else if *nefPtr < ida_mem.ida_maxnef {
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
 */
fn IDAReset(ida_mem: &mut IDAMem) {
    ida_mem.ida_psi[0] = ida_mem.ida_hh;

    let eta = ida_mem.ida_eta;
    ida_mem.ida_phi[1].scale_inplace(eta);
}

/*
 * -----------------------------------------------------------------
 * Function called after a successful step
 * -----------------------------------------------------------------
 */

/*
 * IDACompleteStep
 *
 * This routine completes a successful step.  It increments nst,
 * saves the stepsize and order used, makes the final selection of
 * stepsize and order for the next step, and updates the phi array.
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
       of the necessary information is available yet. */

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
        /* err_kp1 is only read on the RAISE path, where the estimation
           branch below has computed it (C leaves it uninitialized) */
        let mut err_kp1 = ZERO;

        /* Set action = LOWER/MAINTAIN/RAISE to specify order decision */

        let action = if ida_mem.ida_knew == ida_mem.ida_kk - 1 {
            /* Already decided to reduce the order */
            LOWER
        } else if ida_mem.ida_kk == ida_mem.ida_maxord {
            /* Already using the maximum order */
            MAINTAIN
        } else if ida_mem.ida_kk + 1 >= ida_mem.ida_ns || kdiff == 1 {
            /* Step size has not been constant or the order was just raised */
            MAINTAIN
        } else {
            /* Estimate the error at order k+1 */

            {
                let kk = ida_mem.ida_kk as usize;
                let IDAMem { ida_phi, ida_ee, ida_tempv1, .. } = ida_mem;
                N_VLinearSum(ONE, ida_ee, -ONE, &ida_phi[kk + 1], ida_tempv1);
            }
            let enorm = IDAWrmsNorm(ida_mem, &ida_mem.ida_tempv1, &ida_mem.ida_ewt,
                                    ida_mem.ida_suppressalg);
            err_kp1 = enorm / (ida_mem.ida_kk + 2) as f64;

            /* Choose among orders k-1, k, k+1 using local truncation error norms. */

            let terr_k = (ida_mem.ida_kk + 1) as f64 * err_k;
            let terr_kp1 = (ida_mem.ida_kk + 2) as f64 * err_kp1;

            if ida_mem.ida_kk == 1 {
                if terr_kp1 >= HALF * terr_k {
                    MAINTAIN
                } else {
                    RAISE
                }
            } else {
                let terr_km1 = ida_mem.ida_kk as f64 * err_km1;
                if terr_km1 <= SUNMIN(terr_k, terr_kp1) {
                    LOWER
                } else if terr_kp1 >= terr_k {
                    MAINTAIN
                } else {
                    RAISE
                }
            }
        };

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
           1. If eta >= eta_max_fx (default = 2), increase hh to at most eta_max
              (default = 2) i.e., double the step size
           2. If eta <= eta_min_fx (default = 1), reduce hh to between eta_min
              (default 0.5) and eta_low (default 0.9),
           3. Otherwise leave hh as is i.e., eta = 1. */

        ida_mem.ida_eta = ONE;
        let tmp = SUNRpowerR(TWO * err_knew + PT0001, -ONE / (ida_mem.ida_kk + 1) as f64);

        if tmp >= ida_mem.ida_eta_max_fx {
            /* Enforce max growth factor bound and max step size */
            ida_mem.ida_eta = SUNMIN(tmp, ida_mem.ida_eta_max);
            ida_mem.ida_eta /= SUNMAX(ONE, ida_mem.ida_eta * SUNRabs(ida_mem.ida_hh) *
                                             ida_mem.ida_hmax_inv);
        } else if tmp <= ida_mem.ida_eta_min_fx {
            /* Enforce required reduction factor bound, min reduction bound, and min
               step size. Note if eta = eta_min_fx = 1 and ida_eta_low < 1 the step
               size is reduced. */
            ida_mem.ida_eta = SUNMIN(tmp, ida_mem.ida_eta_low);
            ida_mem.ida_eta = SUNMAX(ida_mem.ida_eta, ida_mem.ida_eta_min);
            ida_mem.ida_eta = SUNMAX(ida_mem.ida_eta,
                                     ida_mem.ida_hmin / SUNRabs(ida_mem.ida_hh));
        }
        ida_mem.ida_hh *= ida_mem.ida_eta;
    } /* end of phase if block */

    /* Save ee for possible order increase on next step */
    if ida_mem.ida_kused < ida_mem.ida_maxord {
        let kused = ida_mem.ida_kused as usize;
        let IDAMem { ida_phi, ida_ee, .. } = ida_mem;
        N_VScale(ONE, ida_ee, &mut ida_phi[kused + 1]);
    }

    /* Update phi arrays */

    /* To update phi arrays compute X += Z where                  */
    /* X = [ phi[kused], phi[kused-1], phi[kused-2], ... phi[1] ] */
    /* Z = [ ee,         phi[kused],   phi[kused-1], ... phi[0] ] */
    /* (the C fused N_VLinearSumVectorArray processes the pairs in
       order, so each sum after the first reads the vector already
       updated by the previous pair — a running accumulation from
       phi[kused] down to phi[0], replicated sequentially here) */
    {
        let kused = ida_mem.ida_kused as usize;
        let IDAMem { ida_phi, ida_ee, .. } = ida_mem;

        for (x, z) in ida_phi[kused].data.iter_mut().zip(&ida_ee.data) {
            *x = ONE * *x + ONE * *z;
        }
        for j in (0..kused).rev() {
            let (lo, hi) = ida_phi.split_at_mut(j + 1);
            for (x, z) in lo[j].data.iter_mut().zip(&hi[0].data) {
                *x = ONE * *x + ONE * *z;
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
    for j in 1..=(kord as usize) {
        d = d * gam + c / ida_mem.ida_psi[j - 1];
        c = c * gam;
        gam = (delt + ida_mem.ida_psi[j - 1]) / ida_mem.ida_psi[j];

        cvals[j] = c;
        dvals[j - 1] = d;
    }

    /* yret = N_VLinearCombination(kord+1, cvals, phi); the serial
       fused kernel computes z = c0*X0 then accumulates z += cj*Xj
       (donor IDAGetDky replication; it cannot fail, so the
       IDA_VECTOROP_ERR branch of the C code vanishes) */
    let kord = kord as usize;
    for (z, p) in yret.data.iter_mut().zip(&ida_mem.ida_phi[0].data) {
        *z = cvals[0] * *p;
    }
    for j in 1..=kord {
        for (z, p) in yret.data.iter_mut().zip(&ida_mem.ida_phi[j].data) {
            *z += cvals[j] * *p;
        }
    }

    /* ypret = N_VLinearCombination(kord, dvals, phi+1) */
    for (z, p) in ypret.data.iter_mut().zip(&ida_mem.ida_phi[1].data) {
        *z = dvals[0] * *p;
    }
    for j in 2..=kord {
        for (z, p) in ypret.data.iter_mut().zip(&ida_mem.ida_phi[j].data) {
            *z += dvals[j - 1] * *p;
        }
    }

    IDA_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Norm function
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
 * -----------------------------------------------------------------
 * Functions for rootfinding
 * -----------------------------------------------------------------
 */

/* The C call sites `IDAGetSolution(IDA_mem, t, IDA_mem->ida_yy,
   IDA_mem->ida_yp)` write the interpolant into the integrator's own
   (owned) yy/yp vectors; they are detached for the call (donor
   take() pattern) since IDAGetSolution borrows the memory block. */
fn ida_get_solution_into_yyyp(ida_mem: &mut IDAMem, t: f64) -> i32 {
    let mut yy = std::mem::take(&mut ida_mem.ida_yy);
    let mut yp = std::mem::take(&mut ida_mem.ida_yp);
    let flag = IDAGetSolution(ida_mem, t, &mut yy, &mut yp);
    ida_mem.ida_yy = yy;
    ida_mem.ida_yp = yp;
    flag
}

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
    for i in 0..(ida_mem.ida_nrtfn as usize) {
        ida_mem.ida_iroots[i] = 0;
    }
    ida_mem.ida_tlo = ida_mem.ida_tn;
    ida_mem.ida_ttol = (SUNRabs(ida_mem.ida_tn) + SUNRabs(ida_mem.ida_hh)) *
        ida_mem.ida_uround * HUNDRED;

    /* Evaluate g at initial t and check for zero values. */
    let gfun = ida_mem.ida_gfun.unwrap();
    let retval = {
        let IDAMem { ida_tlo, ida_phi, ida_glo, ida_user_data, .. } = ida_mem;
        gfun(*ida_tlo, &ida_phi[0], &ida_phi[1], ida_glo, ida_user_data)
    };
    ida_mem.ida_nge = 1;
    if retval != 0 {
        return IDA_RTFUNC_FAIL;
    }

    let mut zroot = SUNFALSE;
    for i in 0..(ida_mem.ida_nrtfn as usize) {
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
    {
        let IDAMem { ida_phi, ida_yy, .. } = ida_mem;
        N_VLinearSum(ONE, &ida_phi[0], smallh, &ida_phi[1], ida_yy);
    }
    let retval = {
        let IDAMem { ida_yy, ida_phi, ida_ghi, ida_user_data, .. } = ida_mem;
        gfun(tplus, ida_yy, &ida_phi[1], ida_ghi, ida_user_data)
    };
    ida_mem.ida_nge += 1;
    if retval != 0 {
        return IDA_RTFUNC_FAIL;
    }

    /* We check now only the components of g which were exactly 0.0 at t0
     * to see if we can 'activate' them. */
    for i in 0..(ida_mem.ida_nrtfn as usize) {
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

    let _ = ida_get_solution_into_yyyp(ida_mem, ida_mem.ida_tlo);
    let gfun = ida_mem.ida_gfun.unwrap();
    let retval = {
        let IDAMem { ida_tlo, ida_yy, ida_yp, ida_glo, ida_user_data, .. } = ida_mem;
        gfun(*ida_tlo, ida_yy, ida_yp, ida_glo, ida_user_data)
    };
    ida_mem.ida_nge += 1;
    if retval != 0 {
        return IDA_RTFUNC_FAIL;
    }

    let mut zroot = SUNFALSE;
    for i in 0..(ida_mem.ida_nrtfn as usize) {
        ida_mem.ida_iroots[i] = 0;
    }
    for i in 0..(ida_mem.ida_nrtfn as usize) {
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
    ida_mem.ida_ttol = (SUNRabs(ida_mem.ida_tn) + SUNRabs(ida_mem.ida_hh)) *
        ida_mem.ida_uround * HUNDRED;
    let smallh = if ida_mem.ida_hh > ZERO { ida_mem.ida_ttol } else { -ida_mem.ida_ttol };
    let tplus = ida_mem.ida_tlo + smallh;
    if (tplus - ida_mem.ida_tn) * ida_mem.ida_hh >= ZERO {
        let hratio = smallh / ida_mem.ida_hh;
        let IDAMem { ida_phi, ida_yy, .. } = ida_mem;
        /* N_VLinearSum(ONE, yy, hratio, phi[1], yy) */
        ida_yy.linear_sum_with(ONE, hratio, &ida_phi[1]);
    } else {
        let _ = ida_get_solution_into_yyyp(ida_mem, tplus);
    }
    let retval = {
        let IDAMem { ida_yy, ida_yp, ida_ghi, ida_user_data, .. } = ida_mem;
        gfun(tplus, ida_yy, ida_yp, ida_ghi, ida_user_data)
    };
    ida_mem.ida_nge += 1;
    if retval != 0 {
        return IDA_RTFUNC_FAIL;
    }

    /* Check for close roots (error return), for a new zero at tlo+smallh,
    and for a g_i that changed from zero to nonzero. */
    let mut zroot = SUNFALSE;
    for i in 0..(ida_mem.ida_nrtfn as usize) {
        if !ida_mem.ida_gactive[i] {
            continue;
        }
        if SUNRabs(ida_mem.ida_ghi[i]) == ZERO {
            if ida_mem.ida_iroots[i] == 1 {
                return CLOSERT;
            }
            zroot = SUNTRUE;
            ida_mem.ida_iroots[i] = 1;
        } else if ida_mem.ida_iroots[i] == 1 {
            ida_mem.ida_glo[i] = ida_mem.ida_ghi[i];
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
    let _ = ida_get_solution_into_yyyp(ida_mem, ida_mem.ida_thi);

    /* Set ghi = g(thi) and call IDARootfind to search (tlo,thi) for roots. */
    let gfun = ida_mem.ida_gfun.unwrap();
    let retval = {
        let IDAMem { ida_thi, ida_yy, ida_yp, ida_ghi, ida_user_data, .. } = ida_mem;
        gfun(*ida_thi, ida_yy, ida_yp, ida_ghi, ida_user_data)
    };
    ida_mem.ida_nge += 1;
    if retval != 0 {
        return IDA_RTFUNC_FAIL;
    }

    ida_mem.ida_ttol = (SUNRabs(ida_mem.ida_tn) + SUNRabs(ida_mem.ida_hh)) *
        ida_mem.ida_uround * HUNDRED;
    let ier = IDARootfind(ida_mem);
    if ier == IDA_RTFUNC_FAIL {
        return IDA_RTFUNC_FAIL;
    }
    for i in 0..(ida_mem.ida_nrtfn as usize) {
        if !ida_mem.ida_gactive[i] && ida_mem.ida_grout[i] != ZERO {
            ida_mem.ida_gactive[i] = SUNTRUE;
        }
    }
    ida_mem.ida_tlo = ida_mem.ida_trout;
    for i in 0..(ida_mem.ida_nrtfn as usize) {
        ida_mem.ida_glo[i] = ida_mem.ida_grout[i];
    }

    /* If no root found, return IDA_SUCCESS. */
    if ier == IDA_SUCCESS {
        return IDA_SUCCESS;
    }

    /* If a root was found, interpolate to get y(trout) and return.  */
    let _ = ida_get_solution_into_yyyp(ida_mem, ida_mem.ida_trout);
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
    let nrt = ida_mem.ida_nrtfn as usize;
    let mut imax = 0usize;

    /* First check for change in sign in ghi or for a zero in ghi. */
    let mut maxfrac = ZERO;
    let mut zroot = SUNFALSE;
    let mut sgnchg = SUNFALSE;
    for i in 0..nrt {
        if !ida_mem.ida_gactive[i] {
            continue;
        }
        if SUNRabs(ida_mem.ida_ghi[i]) == ZERO {
            if ida_mem.ida_rootdir[i] as f64 * ida_mem.ida_glo[i] <= ZERO {
                zroot = SUNTRUE;
            }
        } else if SUNRdifferentsign(ida_mem.ida_glo[i], ida_mem.ida_ghi[i])
            && ida_mem.ida_rootdir[i] as f64 * ida_mem.ida_glo[i] <= ZERO
        {
            let gfrac = SUNRabs(ida_mem.ida_ghi[i] /
                                (ida_mem.ida_ghi[i] - ida_mem.ida_glo[i]));
            if gfrac > maxfrac {
                sgnchg = SUNTRUE;
                maxfrac = gfrac;
                imax = i;
            }
        }
    }

    /* If no sign change was found, reset trout and grout.  Then return
       IDA_SUCCESS if no zero was found, or set iroots and return RTFOUND.  */
    if !sgnchg {
        ida_mem.ida_trout = ida_mem.ida_thi;
        for i in 0..nrt {
            ida_mem.ida_grout[i] = ida_mem.ida_ghi[i];
        }
        if !zroot {
            return IDA_SUCCESS;
        }
        for i in 0..nrt {
            ida_mem.ida_iroots[i] = 0;
            if !ida_mem.ida_gactive[i] {
                continue;
            }
            if SUNRabs(ida_mem.ida_ghi[i]) == ZERO
                && ida_mem.ida_rootdir[i] as f64 * ida_mem.ida_glo[i] <= ZERO
            {
                ida_mem.ida_iroots[i] = if ida_mem.ida_glo[i] > 0.0 { -1 } else { 1 };
            }
        }
        return RTFOUND;
    }

    /* Initialize alph to avoid compiler warning */
    let mut alph = ONE;

    /* A sign change was found.  Loop to locate nearest root. */

    let mut side = 0;
    let mut sideprev = -1;
    let gfun = ida_mem.ida_gfun.unwrap();
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
        let mut tmid = ida_mem.ida_thi -
            (ida_mem.ida_thi - ida_mem.ida_tlo) * ida_mem.ida_ghi[imax] /
                (ida_mem.ida_ghi[imax] - alph * ida_mem.ida_glo[imax]);
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

        let _ = ida_get_solution_into_yyyp(ida_mem, tmid);
        let retval = {
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
        for i in 0..nrt {
            if !ida_mem.ida_gactive[i] {
                continue;
            }
            if SUNRabs(ida_mem.ida_grout[i]) == ZERO {
                if ida_mem.ida_rootdir[i] as f64 * ida_mem.ida_glo[i] <= ZERO {
                    zroot = SUNTRUE;
                }
            } else if SUNRdifferentsign(ida_mem.ida_glo[i], ida_mem.ida_grout[i])
                && ida_mem.ida_rootdir[i] as f64 * ida_mem.ida_glo[i] <= ZERO
            {
                let gfrac = SUNRabs(ida_mem.ida_grout[i] /
                                    (ida_mem.ida_grout[i] - ida_mem.ida_glo[i]));
                if gfrac > maxfrac {
                    sgnchg = SUNTRUE;
                    maxfrac = gfrac;
                    imax = i;
                }
            }
        }
        if sgnchg {
            /* Sign change found in (tlo,tmid); replace thi with tmid. */
            ida_mem.ida_thi = tmid;
            for i in 0..nrt {
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
            for i in 0..nrt {
                ida_mem.ida_ghi[i] = ida_mem.ida_grout[i];
            }
            break;
        }

        /* No sign change in (tlo,tmid), and no zero at tmid.
           Sign change must be in (tmid,thi).  Replace tlo with tmid. */
        ida_mem.ida_tlo = tmid;
        for i in 0..nrt {
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
    for i in 0..nrt {
        ida_mem.ida_grout[i] = ida_mem.ida_ghi[i];
        ida_mem.ida_iroots[i] = 0;
        if !ida_mem.ida_gactive[i] {
            continue;
        }
        if SUNRabs(ida_mem.ida_ghi[i]) == ZERO
            && ida_mem.ida_rootdir[i] as f64 * ida_mem.ida_glo[i] <= ZERO
        {
            ida_mem.ida_iroots[i] = if ida_mem.ida_glo[i] > 0.0 { -1 } else { 1 };
        }
        if SUNRdifferentsign(ida_mem.ida_glo[i], ida_mem.ida_ghi[i])
            && ida_mem.ida_rootdir[i] as f64 * ida_mem.ida_glo[i] <= ZERO
        {
            ida_mem.ida_iroots[i] = if ida_mem.ida_glo[i] > 0.0 { -1 } else { 1 };
        }
    }
    RTFOUND
}

/* (IDAProcessError — the ida.c error message handler — lives in
   ida_impl.rs, together with the message constants.) */
