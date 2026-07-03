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
