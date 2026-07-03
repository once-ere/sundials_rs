/* -----------------------------------------------------------------
 * Translated from src/cvode/cvode.c (CVODE 7.7.0).
 * Main CVODE integrator: creation/initialization, tolerances,
 * rootfinding init, the CVode driver, cvStep and all its helper
 * routines (order/step-size control, error test, BDF stability
 * limit detection, rootfinding, error weights).
 * -----------------------------------------------------------------*/
use crate::cvode_impl::*;
use crate::cvode_ls::{cvLsInitialize, cvLsSetup, cvLsSolve};
use crate::cvode_diag::{cvDiagInit, cvDiagSetup, cvDiagSolve};
use crate::cvode_nls::cvNlsInit;
use crate::cvode_nls::cvNls;
use crate::cvode_proj::{cvDoProjection, cvProjInit};
use crate::nvector_serial::*;
use crate::sundials_context::SUNContext;
use crate::sundials_math::*;
use crate::sundials_nonlinearsolver::SUN_NLS_CONV_RECVR;
use crate::sundials_types::*;
use crate::sunnonlinsol_newton::SUNNonlinSol_Newton;

/*=================================================================*/
/* CVODE Private Constants                                         */
/*=================================================================*/

const ZERO: f64 = 0.0;
const TINY: f64 = 1.0e-10;
const PT1: f64 = 0.1;
const POINT2: f64 = 0.2;
const FOURTH: f64 = 0.25;
const HALF: f64 = 0.5;
const PT9: f64 = 0.9;
const ONE: f64 = 1.0;
const ONEPT5: f64 = 1.5;
const TWO: f64 = 2.0;
const THREE: f64 = 3.0;
const FOUR: f64 = 4.0;
const FIVE: f64 = 5.0;
const TWELVE: f64 = 12.0;
const HUNDRED: f64 = 100.0;

/* Control constants for lower-level rootfinding functions */
pub const RTFOUND: i32 = 1;
pub const CLOSERT: i32 = 3;

/* Algorithmic constants */
const FUZZ_FACTOR: f64 = 100.0;

const HLB_FACTOR: f64 = 100.0;
const HUB_FACTOR: f64 = 0.1;
const H_BIAS: f64 = HALF;
const MAX_ITERS: i32 = 4;

const CORTES: f64 = 0.1;

/*
 * =================================================================
 * Exported Functions Implementation
 * =================================================================
 */

/*
 * CVodeCreate
 *
 * Creates an internal memory block for a problem to be solved by
 * CVODE and returns it (C returns NULL on illegal lmm; the Rust
 * translation panics with the same diagnostic since a null-object
 * cannot exist in safe Rust).
 */
pub fn CVodeCreate(lmm: i32, sunctx: &SUNContext) -> Box<CVodeMem> {
    if lmm != CV_ADAMS && lmm != CV_BDF {
        cvProcessError(None, 0, line!(), "CVodeCreate", file!(), MSGCV_BAD_LMM);
        panic!("{}", MSGCV_BAD_LMM);
    }

    let maxord = if lmm == CV_ADAMS {
        ADAMS_Q_MAX as i32
    } else {
        BDF_Q_MAX as i32
    };

    Box::new(CVodeMem {
        cv_sunctx: sunctx.clone(),
        cv_uround: SUN_UNIT_ROUNDOFF,

        /* Problem specification data (defaults) */
        cv_f: None,
        cv_user_data: None,
        cv_lmm: lmm,
        cv_itol: CV_NN,
        cv_reltol: ZERO,
        cv_Sabstol: ZERO,
        cv_Vabstol: NVector::default(),
        cv_atolmin0: SUNTRUE,
        cv_user_efun: SUNFALSE,
        cv_efun: None,

        cv_zn: Vec::new(),

        cv_ewt: NVector::default(),
        cv_y: NVector::default(),
        cv_acor: NVector::default(),
        cv_tempv: NVector::default(),
        cv_ftemp: NVector::default(),
        cv_vtemp1: NVector::default(),
        cv_vtemp2: NVector::default(),
        cv_vtemp3: NVector::default(),

        cv_tstopset: SUNFALSE,
        cv_tstopinterp: SUNFALSE,
        cv_tstop: ZERO,

        cv_q: 0,
        cv_qprime: 0,
        cv_next_q: 0,
        cv_qwait: 0,
        cv_L: 0,

        cv_hin: ZERO,
        cv_h: ZERO,
        cv_hprime: ZERO,
        cv_next_h: ZERO,
        cv_eta: ZERO,
        cv_hscale: ZERO,
        cv_tn: ZERO,
        cv_tretlast: ZERO,

        cv_tau: [ZERO; L_MAX + 1],
        cv_tq: [ZERO; NUM_TESTS + 1],
        cv_l: [ZERO; L_MAX],

        cv_rl1: ZERO,
        cv_gamma: ZERO,
        cv_gammap: ZERO,
        cv_gamrat: ZERO,

        cv_crate: ZERO,
        cv_delp: ZERO,
        cv_acnrm: ZERO,
        cv_acnrmcur: SUNFALSE,
        cv_nlscoef: CORTES,

        cv_qmax: maxord,
        cv_mxstep: MXSTEP_DEFAULT,
        cv_mxhnil: MXHNIL_DEFAULT,
        cv_maxnef: MXNEF,
        cv_maxncf: MXNCF,

        cv_hmin: HMIN_DEFAULT,
        cv_hmax_inv: HMAX_INV_DEFAULT,
        cv_etamax: ZERO,
        cv_eta_min_fx: ETA_MIN_FX_DEFAULT,
        cv_eta_max_fx: ETA_MAX_FX_DEFAULT,
        cv_eta_max_fs: ETA_MAX_FS_DEFAULT,
        cv_eta_max_es: ETA_MAX_ES_DEFAULT,
        cv_eta_max_gs: ETA_MAX_GS_DEFAULT,
        cv_eta_min: ETA_MIN_DEFAULT,
        cv_eta_min_ef: ETA_MIN_EF_DEFAULT,
        cv_eta_max_ef: ETA_MAX_EF_DEFAULT,
        cv_eta_cf: ETA_CF_DEFAULT,

        cv_small_nst: SMALL_NST_DEFAULT,
        cv_small_nef: SMALL_NEF_DEFAULT,

        cv_nst: 0,
        cv_nfe: 0,
        cv_ncfn: 0,
        cv_nni: 0,
        cv_nnf: 0,
        cv_netf: 0,
        cv_nsetups: 0,
        cv_nhnil: 0,

        cv_etaqm1: ZERO,
        cv_etaq: ZERO,
        cv_etaqp1: ZERO,

        cv_lrw1: 0,
        cv_liw1: 0,
        cv_lrw: (58 + 2 * L_MAX + NUM_TESTS) as i64,
        cv_liw: 40,

        NLS: None,
        cv_nls_curiter: 0,
        ownNLS: SUNFALSE,
        nls_f: None,
        convfail: CV_NO_FAILURES,

        cv_lmem: LsModule::None,
        cv_msbp: MSBP_DEFAULT,
        cv_dgmax_lsetup: DGMAX_LSETUP_DEFAULT,

        cv_qu: 0,
        cv_nstlp: 0,
        cv_h0u: ZERO,
        cv_hu: ZERO,
        cv_saved_tq5: ZERO,
        cv_jcur: SUNFALSE,
        cv_tolsf: ZERO,
        cv_qmax_alloc: maxord,
        cv_indx_acor: 0,

        cv_VabstolMallocDone: SUNFALSE,
        cv_MallocDone: SUNFALSE,

        cv_sldeton: SUNFALSE,
        cv_ssdat: [[ZERO; 4]; 6],
        cv_nscon: 0,
        cv_nor: 0,

        cv_gfun: None,
        cv_nrtfn: 0,
        cv_iroots: Vec::new(),
        cv_rootdir: Vec::new(),
        cv_tlo: ZERO,
        cv_thi: ZERO,
        cv_trout: ZERO,
        cv_glo: Vec::new(),
        cv_ghi: Vec::new(),
        cv_grout: Vec::new(),
        cv_ttol: ZERO,
        cv_taskc: 0,
        cv_irfnd: 0,
        cv_nge: 0,
        cv_gactive: Vec::new(),
        cv_mxgnull: 1,

        cv_constraints: NVector::default(),
        cv_constraintsSet: SUNFALSE,
        constraint_corrections: 0,
        constraint_fails: 0,
        max_constraint_fails: MAX_CONSTRAINT_FAILS,

        proj_mem: None,
        proj_enabled: SUNFALSE,
        proj_applied: SUNFALSE,
        proj_p: [ZERO; L_MAX],

        cv_usefused: SUNFALSE,

        first_step_after_resize: SUNFALSE,
    })
}

/*
 * CVodeInit
 *
 * Allocates and initializes memory for a problem.
 */
pub fn CVodeInit(cv_mem: &mut CVodeMem, f: CVRhsFn, t0: f64, y0: &NVector) -> i32 {
    if y0.is_empty() {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeInit", file!(), MSGCV_NULL_Y0);
        return CV_ILL_INPUT;
    }

    /* Set space requirements for one N_Vector */
    let (lrw1, liw1) = N_VSpace(y0);
    cv_mem.cv_lrw1 = lrw1;
    cv_mem.cv_liw1 = liw1;

    /* Allocate the vectors (using y0 as a template) */
    cvAllocVectors(cv_mem, y0);

    /* Copy the input parameters into CVODE state */
    cv_mem.cv_f = Some(f);
    cv_mem.cv_tn = t0;

    /* Initialize zn[0] in the history array */
    cv_mem.cv_zn[0].data.copy_from_slice(&y0.data);

    /* create a Newton nonlinear solver object by default */
    let nls = SUNNonlinSol_Newton(y0, &cv_mem.cv_sunctx);

    /* attach the nonlinear solver to the CVODE memory */
    let retval = crate::cvode_nls::CVodeSetNonlinearSolver(cv_mem, nls);
    if retval != CV_SUCCESS {
        cvProcessError(Some(cv_mem), retval, line!(), "CVodeInit", file!(),
                       "Setting the nonlinear solver failed");
        return CV_MEM_FAIL;
    }

    /* set ownership flag */
    cv_mem.ownNLS = SUNTRUE;

    /* Set step parameters */
    cv_mem.cv_q = 1;
    cv_mem.cv_L = 2;
    cv_mem.cv_qwait = cv_mem.cv_L;
    cv_mem.cv_etamax = cv_mem.cv_eta_max_fs;

    cv_mem.cv_qu = 0;
    cv_mem.cv_hu = ZERO;
    cv_mem.cv_tolsf = ONE;

    /* Set the linear solver addresses to NULL */
    cv_mem.cv_lmem = LsModule::None;

    /* Initialize all the counters */
    cv_mem.cv_nst = 0;
    cv_mem.cv_nfe = 0;
    cv_mem.cv_ncfn = 0;
    cv_mem.cv_netf = 0;
    cv_mem.cv_nni = 0;
    cv_mem.cv_nnf = 0;
    cv_mem.cv_nsetups = 0;
    cv_mem.cv_nhnil = 0;
    cv_mem.cv_nstlp = 0;
    cv_mem.cv_nscon = 0;
    cv_mem.cv_nge = 0;

    cv_mem.cv_irfnd = 0;

    /* Initialize other integrator optional outputs */
    cv_mem.cv_h0u = ZERO;
    cv_mem.cv_next_h = ZERO;
    cv_mem.cv_next_q = 0;

    /* Initialize Stability Limit Detection data */
    cv_mem.cv_nor = 0;
    for i in 1..=5usize {
        for k in 1..=3usize {
            cv_mem.cv_ssdat[i - 1][k - 1] = ZERO;
        }
    }

    /* Problem has been successfully initialized */
    cv_mem.cv_MallocDone = SUNTRUE;

    CV_SUCCESS
}

/*
 * CVodeReInit
 *
 * Re-initializes CVODE's memory for a problem.
 */
pub fn CVodeReInit(cv_mem: &mut CVodeMem, t0: f64, y0: &NVector) -> i32 {
    if !cv_mem.cv_MallocDone {
        cvProcessError(Some(cv_mem), CV_NO_MALLOC, line!(), "CVodeReInit", file!(), MSGCV_NO_MALLOC);
        return CV_NO_MALLOC;
    }

    if y0.is_empty() {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeReInit", file!(), MSGCV_NULL_Y0);
        return CV_ILL_INPUT;
    }

    /* Copy the input parameters into CVODE state */
    cv_mem.cv_tn = t0;

    /* Set step parameters */
    cv_mem.cv_q = 1;
    cv_mem.cv_L = 2;
    cv_mem.cv_qwait = cv_mem.cv_L;
    cv_mem.cv_etamax = cv_mem.cv_eta_max_fs;

    cv_mem.cv_qu = 0;
    cv_mem.cv_hu = ZERO;
    cv_mem.cv_tolsf = ONE;

    /* Initialize zn[0] in the history array */
    cv_mem.cv_zn[0].data.copy_from_slice(&y0.data);

    /* Initialize all the counters */
    cv_mem.cv_nst = 0;
    cv_mem.cv_nfe = 0;
    cv_mem.cv_ncfn = 0;
    cv_mem.cv_netf = 0;
    cv_mem.cv_nni = 0;
    cv_mem.cv_nnf = 0;
    cv_mem.cv_nsetups = 0;
    cv_mem.cv_nhnil = 0;
    cv_mem.cv_nstlp = 0;
    cv_mem.cv_nscon = 0;
    cv_mem.cv_nge = 0;

    cv_mem.cv_irfnd = 0;

    cv_mem.constraint_corrections = 0;
    cv_mem.constraint_fails = 0;

    /* lreinit: reinitialize the attached linear solver interface */
    cv_lreinit_dispatch(cv_mem);

    /* Initialize other integrator optional outputs */
    cv_mem.cv_h0u = ZERO;
    cv_mem.cv_next_h = ZERO;
    cv_mem.cv_next_q = 0;

    /* Initialize Stability Limit Detection data */
    cv_mem.cv_nor = 0;
    for i in 1..=5usize {
        for k in 1..=3usize {
            cv_mem.cv_ssdat[i - 1][k - 1] = ZERO;
        }
    }

    CV_SUCCESS
}

/*
 * CVodeSStolerances / CVodeSVtolerances / CVodeWFtolerances
 */
pub fn CVodeSStolerances(cv_mem: &mut CVodeMem, reltol: f64, abstol: f64) -> i32 {
    if !cv_mem.cv_MallocDone {
        cvProcessError(Some(cv_mem), CV_NO_MALLOC, line!(), "CVodeSStolerances", file!(), MSGCV_NO_MALLOC);
        return CV_NO_MALLOC;
    }
    if reltol < ZERO {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSStolerances", file!(), MSGCV_BAD_RELTOL);
        return CV_ILL_INPUT;
    }
    if abstol < ZERO {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSStolerances", file!(), MSGCV_BAD_ABSTOL);
        return CV_ILL_INPUT;
    }

    cv_mem.cv_reltol = reltol;
    cv_mem.cv_Sabstol = abstol;
    cv_mem.cv_atolmin0 = abstol == ZERO;

    cv_mem.cv_itol = CV_SS;

    cv_mem.cv_user_efun = SUNFALSE;
    cv_mem.cv_efun = None; /* internal cvEwtSet dispatch is used */

    CV_SUCCESS
}

pub fn CVodeSVtolerances(cv_mem: &mut CVodeMem, reltol: f64, abstol: &NVector) -> i32 {
    if !cv_mem.cv_MallocDone {
        cvProcessError(Some(cv_mem), CV_NO_MALLOC, line!(), "CVodeSVtolerances", file!(), MSGCV_NO_MALLOC);
        return CV_NO_MALLOC;
    }
    if reltol < ZERO {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSVtolerances", file!(), MSGCV_BAD_RELTOL);
        return CV_ILL_INPUT;
    }
    let atolmin = N_VMin(abstol);
    if atolmin < ZERO {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSVtolerances", file!(), MSGCV_BAD_ABSTOL);
        return CV_ILL_INPUT;
    }

    if !cv_mem.cv_VabstolMallocDone {
        cv_mem.cv_Vabstol = N_VClone(&cv_mem.cv_ewt);
        cv_mem.cv_lrw += cv_mem.cv_lrw1;
        cv_mem.cv_liw += cv_mem.cv_liw1;
        cv_mem.cv_VabstolMallocDone = SUNTRUE;
    }

    cv_mem.cv_reltol = reltol;
    cv_mem.cv_Vabstol.data.copy_from_slice(&abstol.data);
    cv_mem.cv_atolmin0 = atolmin == ZERO;

    cv_mem.cv_itol = CV_SV;

    cv_mem.cv_user_efun = SUNFALSE;
    cv_mem.cv_efun = None;

    CV_SUCCESS
}

pub fn CVodeWFtolerances(cv_mem: &mut CVodeMem, efun: CVEwtFn) -> i32 {
    if !cv_mem.cv_MallocDone {
        cvProcessError(Some(cv_mem), CV_NO_MALLOC, line!(), "CVodeWFtolerances", file!(), MSGCV_NO_MALLOC);
        return CV_NO_MALLOC;
    }

    cv_mem.cv_itol = CV_WF;
    cv_mem.cv_user_efun = SUNTRUE;
    cv_mem.cv_efun = Some(efun);

    CV_SUCCESS
}

/*
 * CVodeRootInit
 */
pub fn CVodeRootInit(cv_mem: &mut CVodeMem, nrtfn: i32, g: Option<CVRootFn>) -> i32 {
    let nrt = if nrtfn < 0 { 0 } else { nrtfn };

    /* If rerunning with a different number of root functions, free memory */
    if nrt != cv_mem.cv_nrtfn && cv_mem.cv_nrtfn > 0 {
        cv_mem.cv_glo = Vec::new();
        cv_mem.cv_ghi = Vec::new();
        cv_mem.cv_grout = Vec::new();
        cv_mem.cv_iroots = Vec::new();
        cv_mem.cv_rootdir = Vec::new();
        cv_mem.cv_gactive = Vec::new();

        cv_mem.cv_lrw -= 3 * cv_mem.cv_nrtfn as i64;
        cv_mem.cv_liw -= 3 * cv_mem.cv_nrtfn as i64;
    }

    /* nrtfn == 0: disable rootfinding */
    if nrt == 0 {
        cv_mem.cv_nrtfn = nrt;
        cv_mem.cv_gfun = None;
        return CV_SUCCESS;
    }

    /* Same number of root functions: just (re)set g */
    if nrt == cv_mem.cv_nrtfn {
        match g {
            None => {
                cv_mem.cv_glo = Vec::new();
                cv_mem.cv_ghi = Vec::new();
                cv_mem.cv_grout = Vec::new();
                cv_mem.cv_iroots = Vec::new();
                cv_mem.cv_rootdir = Vec::new();
                cv_mem.cv_gactive = Vec::new();
                cv_mem.cv_lrw -= 3 * nrt as i64;
                cv_mem.cv_liw -= 3 * nrt as i64;
                cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeRootInit", file!(), MSGCV_NULL_G);
                return CV_ILL_INPUT;
            }
            Some(gf) => {
                cv_mem.cv_gfun = Some(gf);
                return CV_SUCCESS;
            }
        }
    }

    /* Set variable values in CVode memory block */
    cv_mem.cv_nrtfn = nrt;
    match g {
        None => {
            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeRootInit", file!(), MSGCV_NULL_G);
            return CV_ILL_INPUT;
        }
        Some(gf) => cv_mem.cv_gfun = Some(gf),
    }

    /* Allocate necessary memory and return */
    let n = nrt as usize;
    cv_mem.cv_glo = vec![ZERO; n];
    cv_mem.cv_ghi = vec![ZERO; n];
    cv_mem.cv_grout = vec![ZERO; n];
    cv_mem.cv_iroots = vec![0; n];
    cv_mem.cv_rootdir = vec![0; n]; /* both directions */
    cv_mem.cv_gactive = vec![SUNTRUE; n]; /* all active */

    cv_mem.cv_lrw += 3 * nrt as i64;
    cv_mem.cv_liw += 3 * nrt as i64;

    CV_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Linear solver module dispatch
 * (In C these are the cv_linit/cv_lsetup/cv_lsolve function
 *  pointers; the module is taken out of CVodeMem during the call so
 *  its routine can borrow the integrator memory mutably.)
 * -----------------------------------------------------------------
 */

pub(crate) fn cv_has_lsetup(cv_mem: &CVodeMem) -> bool {
    /* mirrors C's `cv_mem->cv_lsetup != NULL`: cvLsInitialize NULLs the
       lsetup pointer for matrix-free-without-preconditioner and
       matrix-embedded configurations (setup_disabled) */
    match &cv_mem.cv_lmem {
        LsModule::None => false,
        LsModule::Ls(ls) => !ls.setup_disabled,
        LsModule::Diag(_) => true,
    }
}

pub(crate) fn cv_linit_dispatch(cv_mem: &mut CVodeMem) -> i32 {
    let mut lmem = std::mem::take(&mut cv_mem.cv_lmem);
    let ier = match &mut lmem {
        LsModule::None => 0,
        LsModule::Ls(ls) => cvLsInitialize(cv_mem, ls),
        LsModule::Diag(dm) => cvDiagInit(cv_mem, dm),
    };
    cv_mem.cv_lmem = lmem;
    ier
}

fn cv_lreinit_dispatch(cv_mem: &mut CVodeMem) {
    /* cv_lreinit is only set by the CVLS interface; it resets counters */
    let mut lmem = std::mem::take(&mut cv_mem.cv_lmem);
    if let LsModule::Ls(ls) = &mut lmem {
        crate::cvode_ls::cvLsReInitialize(cv_mem, ls);
    }
    cv_mem.cv_lmem = lmem;
}

pub(crate) fn cv_lsetup_dispatch(
    cv_mem: &mut CVodeMem,
    convfail: i32,
    jcur_ptr: &mut bool,
) -> i32 {
    let mut lmem = std::mem::take(&mut cv_mem.cv_lmem);
    let ier = match &mut lmem {
        LsModule::None => 0,
        LsModule::Ls(ls) => cvLsSetup(cv_mem, ls, convfail, jcur_ptr),
        LsModule::Diag(dm) => cvDiagSetup(cv_mem, dm, convfail, jcur_ptr),
    };
    cv_mem.cv_lmem = lmem;
    ier
}

pub(crate) fn cv_lsolve_dispatch(cv_mem: &mut CVodeMem, b: &mut NVector) -> i32 {
    let mut lmem = std::mem::take(&mut cv_mem.cv_lmem);
    let ier = match &mut lmem {
        LsModule::None => {
            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "cv_lsolve", file!(), MSGCV_LSOLVE_NULL);
            -1
        }
        LsModule::Ls(ls) => cvLsSolve(cv_mem, ls, b),
        LsModule::Diag(dm) => cvDiagSolve(cv_mem, dm, b),
    };
    cv_mem.cv_lmem = lmem;
    ier
}

/*
 * -----------------------------------------------------------------
 * Main solver function
 * -----------------------------------------------------------------
 */

/*
 * CVode
 *
 * Main driver of the CVODE package: integrates over a time interval
 * by calling cvStep to do internal time steps.
 */
pub fn CVode(
    cv_mem: &mut CVodeMem,
    tout: f64,
    yout: &mut NVector,
    tret: &mut f64,
    itask: i32,
) -> i32 {
    /*
     * -------------------------------------
     * 1. Check and process inputs
     * -------------------------------------
     */

    if !cv_mem.cv_MallocDone {
        cvProcessError(Some(cv_mem), CV_NO_MALLOC, line!(), "CVode", file!(), MSGCV_NO_MALLOC);
        return CV_NO_MALLOC;
    }

    if yout.is_empty() {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVode", file!(), MSGCV_YOUT_NULL);
        return CV_ILL_INPUT;
    }

    if itask != CV_NORMAL && itask != CV_ONE_STEP {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVode", file!(), MSGCV_BAD_ITASK);
        return CV_ILL_INPUT;
    }

    cv_mem.cv_taskc = itask;

    /*
     * ----------------------------------------
     * 2. Initializations performed only at the first step (nst=0)
     * ----------------------------------------
     */

    if cv_mem.cv_nst == 0 {
        cv_mem.cv_tretlast = cv_mem.cv_tn;
        *tret = cv_mem.cv_tn;

        /* Check inputs for correctness */
        let ier = cvInitialSetup(cv_mem, tout);
        if ier != CV_SUCCESS {
            return ier;
        }

        /* Call f at (t0,y0), set zn[1] = y'(t0). */
        let f = cv_mem.cv_f.unwrap();
        let retval = {
            let (zn0, zn_rest) = cv_mem.cv_zn.split_at_mut(1);
            f(cv_mem.cv_tn, &zn0[0], &mut zn_rest[0], &mut cv_mem.cv_user_data)
        };
        cv_mem.cv_nfe += 1;
        if retval < 0 {
            cvProcessError(Some(cv_mem), CV_RHSFUNC_FAIL, line!(), "CVode", file!(),
                &format!("At t = {}, the right-hand side routine failed in an unrecoverable manner.", cv_mem.cv_tn));
            return CV_RHSFUNC_FAIL;
        }
        if retval > 0 {
            cvProcessError(Some(cv_mem), CV_FIRST_RHSFUNC_ERR, line!(), "CVode", file!(), MSGCV_RHSFUNC_FIRST);
            return CV_FIRST_RHSFUNC_ERR;
        }

        /* Test input tstop for legality. */
        if cv_mem.cv_tstopset
            && (cv_mem.cv_tstop - cv_mem.cv_tn) * (tout - cv_mem.cv_tn) <= ZERO
        {
            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVode", file!(),
                &format!("The value tstop = {} is behind current t = {} in the direction of integration.",
                         cv_mem.cv_tstop, cv_mem.cv_tn));
            return CV_ILL_INPUT;
        }

        /* Set initial h (from H0 or cvHin). */
        cv_mem.cv_h = cv_mem.cv_hin;
        if cv_mem.cv_h != ZERO && (tout - cv_mem.cv_tn) * cv_mem.cv_h < ZERO {
            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVode", file!(), MSGCV_BAD_H0);
            return CV_ILL_INPUT;
        }
        if cv_mem.cv_h == ZERO {
            let mut tout_hin = tout;
            if cv_mem.cv_tstopset && (tout - cv_mem.cv_tn) * (tout - cv_mem.cv_tstop) > ZERO {
                tout_hin = cv_mem.cv_tstop;
            }
            let hflag = cvHin(cv_mem, tout_hin);
            if hflag != CV_SUCCESS {
                return cvHandleFailure(cv_mem, hflag);
            }
        }

        /* Enforce hmax and hmin */
        let rh = SUNRabs(cv_mem.cv_h) * cv_mem.cv_hmax_inv;
        if rh > ONE {
            cv_mem.cv_h /= rh;
        }
        if SUNRabs(cv_mem.cv_h) < cv_mem.cv_hmin {
            cv_mem.cv_h *= cv_mem.cv_hmin / SUNRabs(cv_mem.cv_h);
        }

        /* Check for approach to tstop */
        if cv_mem.cv_tstopset
            && (cv_mem.cv_tn + cv_mem.cv_h - cv_mem.cv_tstop) * cv_mem.cv_h > ZERO
        {
            cv_mem.cv_h = (cv_mem.cv_tstop - cv_mem.cv_tn) * (ONE - FOUR * cv_mem.cv_uround);
        }

        /* Scale zn[1] by h. */
        cv_mem.cv_hscale = cv_mem.cv_h;
        cv_mem.cv_h0u = cv_mem.cv_h;
        cv_mem.cv_hprime = cv_mem.cv_h;

        let h = cv_mem.cv_h;
        cv_mem.cv_zn[1].scale_inplace(h);

        /* Check for zeros of root function g at and near t0. */
        if cv_mem.cv_nrtfn > 0 {
            let retval = cvRcheck1(cv_mem);
            if retval == CV_RTFUNC_FAIL {
                cvProcessError(Some(cv_mem), CV_RTFUNC_FAIL, line!(), "CVode", file!(),
                    &format!("At t = {}, the rootfinding routine failed in an unrecoverable manner.", cv_mem.cv_tn));
                return CV_RTFUNC_FAIL;
            }
        }
        /* end of first call block */
    } else if cv_mem.first_step_after_resize {
        /* Check if the resized y satisfies the constraints */
        if cv_mem.cv_constraintsSet {
            let conOK = {
                let (c, zn0, tempv) = (
                    &cv_mem.cv_constraints,
                    &cv_mem.cv_zn[0],
                    &mut cv_mem.cv_tempv,
                );
                N_VConstrMask(c, zn0, tempv)
            };
            if !conOK {
                cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVode", file!(),
                               "y does not satisfy the constraints");
                return CV_ILL_INPUT;
            }
        }

        /* Initialize the linear solver */
        if cv_has_lsetup(cv_mem) {
            let ier = cv_linit_dispatch(cv_mem);
            if ier != 0 {
                cvProcessError(Some(cv_mem), CV_LINIT_FAIL, line!(), "CVode", file!(), MSGCV_LINIT_FAIL);
                return CV_LINIT_FAIL;
            }
        }

        /* Initialize the nonlinear solver */
        let ier = cvNlsInit(cv_mem);
        if ier != 0 {
            cvProcessError(Some(cv_mem), CV_NLS_INIT_FAIL, line!(), "CVode", file!(), MSGCV_NLS_INIT_FAIL);
            return CV_NLS_INIT_FAIL;
        }
    }

    /*
     * ------------------------------------------------------
     * 3. At following steps, perform stop tests
     * -------------------------------------------------------
     */

    if cv_mem.cv_nst > 0 {
        let troundoff =
            FUZZ_FACTOR * cv_mem.cv_uround * (SUNRabs(cv_mem.cv_tn) + SUNRabs(cv_mem.cv_h));

        /* First, check for a root in the last step taken */
        if cv_mem.cv_nrtfn > 0 {
            let irfndp = cv_mem.cv_irfnd;

            let retval = cvRcheck2(cv_mem);

            if retval == CLOSERT {
                cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVode", file!(),
                    &format!("Root found at and very near t = {}.", cv_mem.cv_tlo));
                return CV_ILL_INPUT;
            } else if retval == CV_RTFUNC_FAIL {
                cvProcessError(Some(cv_mem), CV_RTFUNC_FAIL, line!(), "CVode", file!(),
                    &format!("At t = {}, the rootfinding routine failed in an unrecoverable manner.", cv_mem.cv_tlo));
                return CV_RTFUNC_FAIL;
            } else if retval == RTFOUND {
                cv_mem.cv_tretlast = cv_mem.cv_tlo;
                *tret = cv_mem.cv_tlo;
                /* cv_y (interpolated root solution) aliases yout in C */
                yout.data.copy_from_slice(&cv_mem.cv_y.data);
                return CV_ROOT_RETURN;
            }

            /* If tn is distinct from tretlast, check remaining interval for roots */
            if SUNRabs(cv_mem.cv_tn - cv_mem.cv_tretlast) > troundoff {
                let retval = cvRcheck3(cv_mem, tout, itask);

                if retval == CV_SUCCESS {
                    /* no root found */
                    cv_mem.cv_irfnd = 0;
                    if irfndp == 1 && itask == CV_ONE_STEP {
                        cv_mem.cv_tretlast = cv_mem.cv_tn;
                        *tret = cv_mem.cv_tn;
                        yout.data.copy_from_slice(&cv_mem.cv_zn[0].data);
                        return CV_SUCCESS;
                    }
                } else if retval == RTFOUND {
                    /* a new root was found */
                    cv_mem.cv_irfnd = 1;
                    cv_mem.cv_tretlast = cv_mem.cv_tlo;
                    *tret = cv_mem.cv_tlo;
                    /* cv_y (interpolated root solution) aliases yout in C */
                    yout.data.copy_from_slice(&cv_mem.cv_y.data);
                    return CV_ROOT_RETURN;
                } else if retval == CV_RTFUNC_FAIL {
                    cvProcessError(Some(cv_mem), CV_RTFUNC_FAIL, line!(), "CVode", file!(),
                        &format!("At t = {}, the rootfinding routine failed in an unrecoverable manner.", cv_mem.cv_tlo));
                    return CV_RTFUNC_FAIL;
                }
            }
        } /* end of root stop check */

        /* Test for tn at tstop or near tstop */
        if cv_mem.cv_tstopset {
            if SUNRabs(cv_mem.cv_tn - cv_mem.cv_tstop) <= troundoff {
                /* Ensure tout >= tstop, otherwise check for tout return below */
                if (tout - cv_mem.cv_tstop) * cv_mem.cv_h >= ZERO
                    || SUNRabs(tout - cv_mem.cv_tstop) <= troundoff
                {
                    if cv_mem.cv_tstopinterp {
                        let ier = cvGetDky_with_temp(cv_mem, cv_mem.cv_tstop, 0, yout);
                        if ier != CV_SUCCESS {
                            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVode", file!(),
                                &format!("The value tstop = {} is behind current t = {} in the direction of integration.",
                                         cv_mem.cv_tstop, cv_mem.cv_tn));
                            return CV_ILL_INPUT;
                        }
                    } else {
                        yout.data.copy_from_slice(&cv_mem.cv_zn[0].data);
                    }
                    cv_mem.cv_tretlast = cv_mem.cv_tstop;
                    *tret = cv_mem.cv_tstop;
                    cv_mem.cv_tstopset = SUNFALSE;
                    return CV_TSTOP_RETURN;
                }
            }
            /* If next step would overtake tstop, adjust stepsize */
            else if (cv_mem.cv_tn + cv_mem.cv_hprime - cv_mem.cv_tstop) * cv_mem.cv_h > ZERO {
                cv_mem.cv_hprime =
                    (cv_mem.cv_tstop - cv_mem.cv_tn) * (ONE - FOUR * cv_mem.cv_uround);
                cv_mem.cv_eta = cv_mem.cv_hprime / cv_mem.cv_h;
            }
        }

        /* In CV_NORMAL mode, test if tout was reached */
        if itask == CV_NORMAL && (cv_mem.cv_tn - tout) * cv_mem.cv_h >= ZERO {
            cv_mem.cv_tretlast = tout;
            *tret = tout;
            let ier = cvGetDky_with_temp(cv_mem, tout, 0, yout);
            if ier != CV_SUCCESS {
                cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVode", file!(),
                    &format!("Trouble interpolating at tout = {}. tout too far back in direction of integration", tout));
                return CV_ILL_INPUT;
            }
            return CV_SUCCESS;
        }

        /* In CV_ONE_STEP mode, test if tn was returned */
        if itask == CV_ONE_STEP && SUNRabs(cv_mem.cv_tn - cv_mem.cv_tretlast) > troundoff {
            cv_mem.cv_tretlast = cv_mem.cv_tn;
            *tret = cv_mem.cv_tn;
            yout.data.copy_from_slice(&cv_mem.cv_zn[0].data);
            return CV_SUCCESS;
        }
    } /* end stopping tests block */

    /*
     * --------------------------------------------------
     * 4. Looping point for internal steps
     * --------------------------------------------------
     */

    let mut nstloc: i64 = 0;
    let istate;
    loop {
        cv_mem.cv_next_h = cv_mem.cv_h;
        cv_mem.cv_next_q = cv_mem.cv_q;

        /* Reset and check ewt */
        if cv_mem.cv_nst > 0 {
            let ewtset_flag = cv_efun_apply_to_ewt(cv_mem);

            if ewtset_flag != 0 {
                if cv_mem.cv_itol == CV_WF {
                    cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVode", file!(),
                        &format!("At t = {}, the user-provide EwtSet function failed.", cv_mem.cv_tn));
                } else {
                    cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVode", file!(),
                        &format!("At t = {}, a component of ewt has become <= 0.", cv_mem.cv_tn));
                }
                istate = CV_ILL_INPUT;
                cv_mem.cv_tretlast = cv_mem.cv_tn;
                *tret = cv_mem.cv_tn;
                yout.data.copy_from_slice(&cv_mem.cv_zn[0].data);
                break;
            }
        }

        /* Check for too many steps */
        if cv_mem.cv_mxstep > 0 && nstloc >= cv_mem.cv_mxstep {
            cvProcessError(Some(cv_mem), CV_TOO_MUCH_WORK, line!(), "CVode", file!(),
                &format!("At t = {}, mxstep steps taken before reaching tout.", cv_mem.cv_tn));
            istate = CV_TOO_MUCH_WORK;
            cv_mem.cv_tretlast = cv_mem.cv_tn;
            *tret = cv_mem.cv_tn;
            yout.data.copy_from_slice(&cv_mem.cv_zn[0].data);
            break;
        }

        /* Check for too much accuracy requested */
        let nrm = N_VWrmsNorm(&cv_mem.cv_zn[0], &cv_mem.cv_ewt);
        cv_mem.cv_tolsf = cv_mem.cv_uround * nrm;
        if cv_mem.cv_tolsf > ONE {
            cvProcessError(Some(cv_mem), CV_TOO_MUCH_ACC, line!(), "CVode", file!(),
                &format!("At t = {}, too much accuracy requested.", cv_mem.cv_tn));
            istate = CV_TOO_MUCH_ACC;
            cv_mem.cv_tretlast = cv_mem.cv_tn;
            *tret = cv_mem.cv_tn;
            yout.data.copy_from_slice(&cv_mem.cv_zn[0].data);
            cv_mem.cv_tolsf *= TWO;
            break;
        } else {
            cv_mem.cv_tolsf = ONE;
        }

        /* Check for h below roundoff level in tn */
        if cv_mem.cv_tn + cv_mem.cv_h == cv_mem.cv_tn {
            cv_mem.cv_nhnil += 1;
            if cv_mem.cv_nhnil <= cv_mem.cv_mxhnil {
                cvProcessError(Some(cv_mem), CV_WARNING, line!(), "CVode", file!(),
                    &format!("Internal t = {} and h = {} are such that t + h = t on the next step. The solver will continue anyway.",
                             cv_mem.cv_tn, cv_mem.cv_h));
            }
            if cv_mem.cv_nhnil == cv_mem.cv_mxhnil {
                cvProcessError(Some(cv_mem), CV_WARNING, line!(), "CVode", file!(), MSGCV_HNIL_DONE);
            }
        }

        /* Call cvStep to take a step */
        let kflag = cvStep(cv_mem);

        /* Process failed step cases, and exit loop */
        if kflag != CV_SUCCESS {
            istate = cvHandleFailure(cv_mem, kflag);
            cv_mem.cv_tretlast = cv_mem.cv_tn;
            *tret = cv_mem.cv_tn;
            yout.data.copy_from_slice(&cv_mem.cv_zn[0].data);
            break;
        }

        nstloc += 1;

        /* If tstop is set and was reached, reset tn = tstop */
        if cv_mem.cv_tstopset {
            let troundoff =
                FUZZ_FACTOR * cv_mem.cv_uround * (SUNRabs(cv_mem.cv_tn) + SUNRabs(cv_mem.cv_h));
            if SUNRabs(cv_mem.cv_tn - cv_mem.cv_tstop) <= troundoff {
                cv_mem.cv_tn = cv_mem.cv_tstop;
            }
        }

        /* Check for root in last step taken. */
        if cv_mem.cv_nrtfn > 0 {
            let retval = cvRcheck3(cv_mem, tout, itask);

            if retval == RTFOUND {
                /* A new root was found */
                cv_mem.cv_irfnd = 1;
                istate = CV_ROOT_RETURN;
                cv_mem.cv_tretlast = cv_mem.cv_tlo;
                *tret = cv_mem.cv_tlo;
                /* cv_y (interpolated root solution) aliases yout in C */
                yout.data.copy_from_slice(&cv_mem.cv_y.data);
                break;
            } else if retval == CV_RTFUNC_FAIL {
                /* g failed */
                cvProcessError(Some(cv_mem), CV_RTFUNC_FAIL, line!(), "CVode", file!(),
                    &format!("At t = {}, the rootfinding routine failed in an unrecoverable manner.", cv_mem.cv_tlo));
                istate = CV_RTFUNC_FAIL;
                break;
            }

            /* Warn about inactive event functions at the end of the first step */
            if cv_mem.cv_nst == 1 {
                let inactive_roots = cv_mem
                    .cv_gactive
                    .iter()
                    .take(cv_mem.cv_nrtfn as usize)
                    .any(|&a| !a);
                if cv_mem.cv_mxgnull > 0 && inactive_roots {
                    cvProcessError(Some(cv_mem), CV_WARNING, line!(), "CVode", file!(), MSGCV_INACTIVE_ROOTS);
                }
            }
        }

        /* Check if tn is at tstop or near tstop */
        if cv_mem.cv_tstopset {
            let troundoff =
                FUZZ_FACTOR * cv_mem.cv_uround * (SUNRabs(cv_mem.cv_tn) + SUNRabs(cv_mem.cv_h));

            /* Test for tn at tstop */
            if SUNRabs(cv_mem.cv_tn - cv_mem.cv_tstop) <= troundoff {
                /* Ensure tout >= tstop, otherwise check for tout return below */
                if (tout - cv_mem.cv_tstop) * cv_mem.cv_h >= ZERO
                    || SUNRabs(tout - cv_mem.cv_tstop) <= troundoff
                {
                    if cv_mem.cv_tstopinterp {
                        let _ = cvGetDky_with_temp(cv_mem, cv_mem.cv_tstop, 0, yout);
                    } else {
                        yout.data.copy_from_slice(&cv_mem.cv_zn[0].data);
                    }
                    cv_mem.cv_tretlast = cv_mem.cv_tstop;
                    *tret = cv_mem.cv_tstop;
                    cv_mem.cv_tstopset = SUNFALSE;
                    istate = CV_TSTOP_RETURN;
                    break;
                }
            }
            /* If next step would overtake tstop, adjust stepsize */
            else if (cv_mem.cv_tn + cv_mem.cv_hprime - cv_mem.cv_tstop) * cv_mem.cv_h > ZERO {
                cv_mem.cv_hprime =
                    (cv_mem.cv_tstop - cv_mem.cv_tn) * (ONE - FOUR * cv_mem.cv_uround);
                cv_mem.cv_eta = cv_mem.cv_hprime / cv_mem.cv_h;
            }
        }

        /* In NORMAL mode, check if tout reached */
        if itask == CV_NORMAL && (cv_mem.cv_tn - tout) * cv_mem.cv_h >= ZERO {
            istate = CV_SUCCESS;
            cv_mem.cv_tretlast = tout;
            *tret = tout;
            let _ = cvGetDky_with_temp(cv_mem, tout, 0, yout);
            cv_mem.cv_next_q = cv_mem.cv_qprime;
            cv_mem.cv_next_h = cv_mem.cv_hprime;
            break;
        }

        /* In ONE_STEP mode, copy y and exit loop */
        if itask == CV_ONE_STEP {
            istate = CV_SUCCESS;
            cv_mem.cv_tretlast = cv_mem.cv_tn;
            *tret = cv_mem.cv_tn;
            yout.data.copy_from_slice(&cv_mem.cv_zn[0].data);
            cv_mem.cv_next_q = cv_mem.cv_qprime;
            cv_mem.cv_next_h = cv_mem.cv_hprime;
            break;
        }
    } /* end looping for internal steps */

    istate
}

/*
 * -----------------------------------------------------------------
 * Interpolated output and extraction functions
 * -----------------------------------------------------------------
 */

/*
 * CVodeGetDky
 *
 * Computes the k-th derivative of the interpolating polynomial at
 * the time t and stores the result in dky:
 *         q
 *  dky = SUM c(j,k) * (t - tn)^(j-k) * h^(-j) * zn[j]
 *        j=k
 */
pub fn CVodeGetDky(cv_mem: &CVodeMem, t: f64, k: i32, dky: &mut NVector) -> i32 {
    if dky.is_empty() {
        cvProcessError(Some(cv_mem), CV_BAD_DKY, line!(), "CVodeGetDky", file!(), MSGCV_NULL_DKY);
        return CV_BAD_DKY;
    }

    if k < 0 || k > cv_mem.cv_q {
        cvProcessError(Some(cv_mem), CV_BAD_K, line!(), "CVodeGetDky", file!(), MSGCV_BAD_K);
        return CV_BAD_K;
    }

    /* Allow for some slack */
    let mut tfuzz =
        FUZZ_FACTOR * cv_mem.cv_uround * (SUNRabs(cv_mem.cv_tn) + SUNRabs(cv_mem.cv_hu));
    if cv_mem.cv_hu < ZERO {
        tfuzz = -tfuzz;
    }
    let tp = cv_mem.cv_tn - cv_mem.cv_hu - tfuzz;
    let tn1 = cv_mem.cv_tn + tfuzz;
    if (t - tp) * (t - tn1) > ZERO {
        cvProcessError(Some(cv_mem), CV_BAD_T, line!(), "CVodeGetDky", file!(),
            &format!("Illegal value for t. t = {} is not between tcur - hold = {} and tcur = {}",
                     t, cv_mem.cv_tn - cv_mem.cv_hu, cv_mem.cv_tn));
        return CV_BAD_T;
    }

    /* Sum the differentiated interpolating polynomial (the C code uses
       the fused N_VLinearCombination; its serial kernel accumulates
       z = c0*X0 then z += ci*Xi, replicated here). */
    let s = (t - cv_mem.cv_tn) / cv_mem.cv_h;
    let mut first = true;
    let mut j = cv_mem.cv_q;
    while j >= k {
        let mut c = ONE;
        let mut i = j;
        while i >= j - k + 1 {
            c *= i as f64;
            i -= 1;
        }
        for _ in 0..(j - k) {
            c *= s;
        }
        let znj = &cv_mem.cv_zn[j as usize];
        if first {
            for (d, z) in dky.data.iter_mut().zip(&znj.data) {
                *d = c * *z;
            }
            first = false;
        } else {
            for (d, z) in dky.data.iter_mut().zip(&znj.data) {
                *d += c * *z;
            }
        }
        j -= 1;
    }

    if k == 0 {
        return CV_SUCCESS;
    }
    let r = SUNRpowerI(cv_mem.cv_h, -k);
    dky.scale_inplace(r);

    CV_SUCCESS
}

/* Internal helper: call CVodeGetDky writing into a vector while cv_mem
   is mutably held (dky is disjoint from cv_mem). */
pub(crate) fn cvGetDky_with_temp(cv_mem: &mut CVodeMem, t: f64, k: i32, dky: &mut NVector) -> i32 {
    CVodeGetDky(cv_mem, t, k, dky)
}

/* Internal helper: interpolate into cv_mem.cv_y (used by rootfinding).
   Takes cv_y out of the struct to satisfy the borrow checker; cv_y is
   never an input of the interpolation. */
pub(crate) fn cvGetDky_into_y(cv_mem: &mut CVodeMem, t: f64, k: i32) -> i32 {
    let mut y = std::mem::take(&mut cv_mem.cv_y);
    let ier = CVodeGetDky(cv_mem, t, k, &mut y);
    cv_mem.cv_y = y;
    ier
}

/*
 * CVodeComputeState
 *
 * Computes y based on the current prediction and given correction.
 */
pub fn CVodeComputeState(cv_mem: &CVodeMem, ycor: &NVector, y: &mut NVector) -> i32 {
    N_VLinearSum(ONE, &cv_mem.cv_zn[0], ONE, ycor, y);
    CV_SUCCESS
}

/*
 * CVodeFree
 *
 * Frees the problem memory allocated by CVodeInit (RAII: dropping
 * the Box releases everything the C code frees explicitly).
 */
pub fn CVodeFree(_cvode_mem: Box<CVodeMem>) {}

/*
 * =================================================================
 *  Private Functions Implementation
 * =================================================================
 */

/*
 * cvAllocVectors
 *
 * Allocates the CVODE vectors ewt, acor, tempv, ftemp, vtemp1-3 and
 * zn[0], ..., zn[qmax].
 */
fn cvAllocVectors(cv_mem: &mut CVodeMem, tmpl: &NVector) {
    cv_mem.cv_ewt = N_VClone(tmpl);
    cv_mem.cv_acor = N_VClone(tmpl);
    cv_mem.cv_tempv = N_VClone(tmpl);
    cv_mem.cv_ftemp = N_VClone(tmpl);
    cv_mem.cv_vtemp1 = N_VClone(tmpl);
    cv_mem.cv_vtemp2 = N_VClone(tmpl);
    cv_mem.cv_vtemp3 = N_VClone(tmpl);
    cv_mem.cv_y = N_VClone(tmpl);

    /* Allocate zn[0] ... zn[qmax] (the full L_MAX array is allocated so
       CVodeSetMaxOrd may raise qmax later exactly as in C, where the
       zn array has L_MAX slots) */
    cv_mem.cv_zn = (0..L_MAX).map(|_| N_VClone(tmpl)).collect();

    /* Update solver workspace lengths */
    cv_mem.cv_lrw += (cv_mem.cv_qmax as i64 + 8) * cv_mem.cv_lrw1;
    cv_mem.cv_liw += (cv_mem.cv_qmax as i64 + 8) * cv_mem.cv_liw1;

    /* Store the value of qmax used here */
    cv_mem.cv_qmax_alloc = cv_mem.cv_qmax;
}

/*
 * -----------------------------------------------------------------
 * Initial setup
 * -----------------------------------------------------------------
 */

fn cvInitialSetup(cv_mem: &mut CVodeMem, tout: f64) -> i32 {
    /* Is tout too close to tn? */
    let tdist = SUNRabs(tout - cv_mem.cv_tn);
    let tround = cv_mem.cv_uround * SUNMAX(SUNRabs(cv_mem.cv_tn), SUNRabs(tout));

    if tdist == ZERO || tdist < TWO * tround {
        cvProcessError(Some(cv_mem), CV_TOO_CLOSE, line!(), "cvInitialSetup", file!(), MSGCV_TOO_CLOSE);
        return CV_TOO_CLOSE;
    }

    /* Did the user specify tolerances? */
    if cv_mem.cv_itol == CV_NN {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "cvInitialSetup", file!(), MSGCV_NO_TOL);
        return CV_ILL_INPUT;
    }

    /* Check to see if y0 satisfies constraints */
    if cv_mem.cv_constraintsSet {
        let conOK = {
            let (c, zn0, tempv) = (
                &cv_mem.cv_constraints,
                &cv_mem.cv_zn[0],
                &mut cv_mem.cv_tempv,
            );
            N_VConstrMask(c, zn0, tempv)
        };
        if !conOK {
            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "cvInitialSetup", file!(), MSGCV_Y0_FAIL_CONSTR);
            return CV_ILL_INPUT;
        }
    }

    /* Load initial error weights */
    let ier = cv_efun_apply_to_ewt(cv_mem);
    if ier != 0 {
        if cv_mem.cv_itol == CV_WF {
            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "cvInitialSetup", file!(), MSGCV_EWT_FAIL);
        } else {
            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "cvInitialSetup", file!(), MSGCV_BAD_EWT);
        }
        return CV_ILL_INPUT;
    }

    /* Call linit function (if it exists) */
    if cv_has_lsetup(cv_mem) {
        let ier = cv_linit_dispatch(cv_mem);
        if ier != 0 {
            cvProcessError(Some(cv_mem), CV_LINIT_FAIL, line!(), "cvInitialSetup", file!(), MSGCV_LINIT_FAIL);
            return CV_LINIT_FAIL;
        }
    }

    /* Initialize the nonlinear solver */
    let ier = cvNlsInit(cv_mem);
    if ier != 0 {
        cvProcessError(Some(cv_mem), CV_NLS_INIT_FAIL, line!(), "cvInitialSetup", file!(), MSGCV_NLS_INIT_FAIL);
        return CV_NLS_INIT_FAIL;
    }

    /* Initialize projection data */
    if cv_mem.proj_enabled && cv_mem.proj_mem.is_none() {
        cvProcessError(Some(cv_mem), CV_PROJ_MEM_NULL, line!(), "cvInitialSetup", file!(),
                       "proj_mem = NULL illegal.");
        return CV_PROJ_MEM_NULL;
    }

    if let Some(pm) = cv_mem.proj_mem.as_deref_mut() {
        let ier = cvProjInit(pm);
        if ier != CV_SUCCESS {
            cvProcessError(Some(cv_mem), CV_MEM_FAIL, line!(), "cvInitialSetup", file!(), MSGCV_MEM_FAIL);
            return CV_MEM_FAIL;
        }
        cv_mem.proj_applied = SUNFALSE;
    }

    CV_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Initial stepsize calculation
 * -----------------------------------------------------------------
 */

/*
 * cvHin
 *
 * Computes a tentative initial step size h0 as an approximate
 * solution of (WRMS norm of (h0^2 ydd / 2)) = 1.
 */
fn cvHin(cv_mem: &mut CVodeMem, tout: f64) -> i32 {
    /* cvInitialSetup checks for tdiff = 0 or < 2 * troundoff */
    let tdiff = tout - cv_mem.cv_tn;
    let sign = if tdiff > ZERO { 1 } else { -1 };
    let tdist = SUNRabs(tdiff);
    let tround = cv_mem.cv_uround * SUNMAX(SUNRabs(cv_mem.cv_tn), SUNRabs(tout));

    /* Set lower and upper bounds on h0, take geometric mean as first trial */
    let hlb = HLB_FACTOR * tround;
    let hub = cvUpperBoundH0(cv_mem, tdist);

    let mut hg = SUNRsqrt(hlb * hub);

    if hub < hlb {
        if sign == -1 {
            cv_mem.cv_h = -hg;
        } else {
            cv_mem.cv_h = hg;
        }
        return CV_SUCCESS;
    }

    /* Outer loop */
    let mut hs = hg; /* safeguard against 'uninitialized variable' warning */
    let mut hnew = hs;

    let mut count1 = 1;
    while count1 <= MAX_ITERS {
        /* Attempts to estimate ydd */
        let mut hg_ok = SUNFALSE;
        let mut yddnrm = ZERO;

        for _count2 in 1..=MAX_ITERS {
            let hgs = hg * sign as f64;
            let retval = cvYddNorm(cv_mem, hgs, &mut yddnrm);
            /* If the RHS function failed unrecoverably, give up */
            if retval < 0 {
                return CV_RHSFUNC_FAIL;
            }
            /* If successful, we can use ydd */
            if retval == CV_SUCCESS {
                hg_ok = SUNTRUE;
                break;
            }
            /* f failed recoverably; cut step size and test again */
            hg *= POINT2;
        }

        /* If f failed recoverably MAX_ITERS times */
        if !hg_ok {
            /* Exit if this is the first or second pass. No recovery possible */
            if count1 <= 2 {
                return CV_REPTD_RHSFUNC_ERR;
            }
            /* We have a fall-back option: hs is a previous hnew which
               passed through f(). Use it and break */
            hnew = hs;
            break;
        }

        /* The proposed step size is feasible. Save it. */
        hs = hg;

        /* Propose new step size */
        hnew = if yddnrm * hub * hub > TWO {
            SUNRsqrt(TWO / yddnrm)
        } else {
            SUNRsqrt(hg * hub)
        };

        /* If last pass, stop now with hnew */
        if count1 == MAX_ITERS {
            break;
        }

        let hrat = hnew / hg;

        /* Accept hnew if it does not differ from hg by more than a factor of 2 */
        if hrat > HALF && hrat < TWO {
            break;
        }

        /* After one pass, if ydd seems to be bad, use fall-back value. */
        if count1 > 1 && hrat > TWO {
            hnew = hg;
            break;
        }

        /* Send this value back through f() */
        hg = hnew;
        count1 += 1;
    }

    /* Apply bounds, bias factor, and attach sign */
    let mut h0 = H_BIAS * hnew;
    if h0 < hlb {
        h0 = hlb;
    }
    if h0 > hub {
        h0 = hub;
    }
    if sign == -1 {
        h0 = -h0;
    }
    cv_mem.cv_h = h0;

    CV_SUCCESS
}

/*
 * cvUpperBoundH0
 *
 * Sets an upper bound on |h0| based on tdist = tn - t0 and y[i]/y'[i].
 */
fn cvUpperBoundH0(cv_mem: &mut CVodeMem, tdist: f64) -> f64 {
    /*
     * Bound based on |y0|/|y0'| -- allow at most an increase of
     * HUB_FACTOR in y0 (based on a forward Euler step). The weight
     * factor is used as a safeguard against zero components in y0.
     */
    /* temp1 = cv_tempv, temp2 = cv_acor in the C code */

    {
        let (zn0, acor) = (&cv_mem.cv_zn[0], &mut cv_mem.cv_acor);
        N_VAbs(zn0, acor);
    }
    /* efun writes the weight vector into tempv */
    {
        let mut w = std::mem::take(&mut cv_mem.cv_tempv);
        let _ = cv_efun_dispatch(cv_mem, &mut w);
        cv_mem.cv_tempv = w;
    }
    cv_mem.cv_tempv.invert_inplace();
    {
        let (acor, tempv) = (&cv_mem.cv_acor, &mut cv_mem.cv_tempv);
        /* temp1 = HUB_FACTOR*temp2 + temp1 */
        tempv.linear_sum_with(ONE, HUB_FACTOR, acor);
    }

    {
        let (zn1, acor) = (&cv_mem.cv_zn[1], &mut cv_mem.cv_acor);
        N_VAbs(zn1, acor);
    }

    {
        let (acor, tempv) = (&cv_mem.cv_acor, &mut cv_mem.cv_tempv);
        /* temp1 = temp2 / temp1 */
        for (t1, t2) in tempv.data.iter_mut().zip(&acor.data) {
            *t1 = *t2 / *t1;
        }
    }
    let hub_inv = N_VMaxNorm(&cv_mem.cv_tempv);

    /* bound based on tdist -- allow at most a step of HUB_FACTOR * tdist */
    let mut hub = HUB_FACTOR * tdist;

    /* Use the smaller of the two */
    if hub * hub_inv > ONE {
        hub = ONE / hub_inv;
    }

    hub
}

/*
 * cvYddNorm
 *
 * Computes an estimate of the second derivative of y using a
 * difference quotient, and returns its WRMS norm.
 */
fn cvYddNorm(cv_mem: &mut CVodeMem, hg: f64, yddnrm: &mut f64) -> i32 {
    {
        let (zn, y) = (&cv_mem.cv_zn, &mut cv_mem.cv_y);
        N_VLinearSum(hg, &zn[1], ONE, &zn[0], y);
    }
    let f = cv_mem.cv_f.unwrap();
    let retval = f(
        cv_mem.cv_tn + hg,
        &cv_mem.cv_y,
        &mut cv_mem.cv_tempv,
        &mut cv_mem.cv_user_data,
    );
    cv_mem.cv_nfe += 1;
    if retval < 0 {
        return CV_RHSFUNC_FAIL;
    }
    if retval > 0 {
        return RHSFUNC_RECVR;
    }

    {
        let (zn, tempv) = (&cv_mem.cv_zn, &mut cv_mem.cv_tempv);
        tempv.linear_sum_with(ONE / hg, -ONE / hg, &zn[1]);
    }

    *yddnrm = N_VWrmsNorm(&cv_mem.cv_tempv, &cv_mem.cv_ewt);

    CV_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Main cvStep function
 * -----------------------------------------------------------------
 */

/*
 * cvStep
 *
 * Performs one internal cvode step, from tn to tn + h.
 */
fn cvStep(cv_mem: &mut CVodeMem) -> i32 {
    /* Initialize failure counters for this step attempt */
    let mut ncf = 0; /* corrector failures  */
    let mut npf = 0; /* projection failures */
    let mut nef = 0; /* error test failures */
    let mut step_constraint_fails = 0;

    /* If the step size has changed, update the history array */
    if cv_mem.cv_nst > 0 && cv_mem.cv_hprime != cv_mem.cv_h {
        cvAdjustParams(cv_mem);
    }

    /* Check if this step should be projected */
    let mut do_projection = SUNFALSE;
    if cv_mem.proj_enabled {
        let pm = cv_mem.proj_mem.as_deref().unwrap();
        do_projection = pm.freq > 0
            && (cv_mem.cv_nst == 0 || cv_mem.cv_nst >= pm.nstlprj + pm.freq);
    }

    /* Looping point for attempts to take a step */
    let saved_t = cv_mem.cv_tn; /* tn is updated in cvPredict */
    let mut nflag = FIRST_CALL;
    let mut dsm = ZERO;

    loop {
        cvPredict(cv_mem);
        cvSet(cv_mem);

        nflag = cvNls(cv_mem, nflag);
        let kflag = cvHandleNFlag(cv_mem, &mut nflag, saved_t, &mut ncf);

        /* Go back in loop if we need to predict again (nflag=PREV_CONV_FAIL) */
        if kflag == PREDICT_AGAIN {
            continue;
        }

        /* Return if nonlinear solve failed and recovery is not possible. */
        if kflag != DO_ERROR_TEST {
            return kflag;
        }

        /* Check inequality constraints */
        if cv_mem.cv_constraintsSet {
            let cflag =
                cvCheckConstraints(cv_mem, &mut nflag, saved_t, &mut step_constraint_fails);

            /* Go back in loop if we need to predict again */
            if cflag == PREDICT_AGAIN {
                continue;
            }

            /* Return if the check failed and recovery is not possible. */
            if cflag != CV_SUCCESS {
                return cflag;
            }
        }

        /* Check if a projection needs to be performed */
        cv_mem.proj_applied = SUNFALSE;

        if do_projection {
            /* Perform projection (nflag=CV_SUCCESS) */
            let pflag = cvDoProjection(cv_mem, &mut nflag, saved_t, &mut npf);

            /* Go back in loop if we need to predict again (nflag=PREV_PROJ_FAIL) */
            if pflag == PREDICT_AGAIN {
                continue;
            }

            /* Return if projection failed and recovery is not possible */
            if pflag != CV_SUCCESS {
                return pflag;
            }
        }

        /* Perform error test (nflag=CV_SUCCESS) */
        let eflag = cvDoErrorTest(cv_mem, &mut nflag, saved_t, &mut nef, &mut dsm);

        /* Go back in loop if we need to predict again (nflag=PREV_ERR_FAIL) */
        if eflag == TRY_AGAIN {
            continue;
        }

        /* Return if error test failed and recovery is not possible. */
        if eflag != CV_SUCCESS {
            return eflag;
        }

        /* Error test passed (eflag=CV_SUCCESS), break from loop */
        break;
    }

    /* Nonlinear system solve and error test were both successful.
       Update data, and consider change of step and/or order. */

    cvCompleteStep(cv_mem);

    cvPrepareNextStep(cv_mem, dsm);

    /* If Stability Limit Detection is turned on, check stability */
    if cv_mem.cv_sldeton {
        cvBDFStab(cv_mem);
    }

    cv_mem.cv_etamax = if cv_mem.cv_nst <= cv_mem.cv_small_nst {
        cv_mem.cv_eta_max_es
    } else {
        cv_mem.cv_eta_max_gs
    };

    /* Finally, we rescale the acor array to be the
       estimated local error vector. */
    let tq2 = cv_mem.cv_tq[2];
    cv_mem.cv_acor.scale_inplace(tq2);

    CV_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function called at beginning of step
 * -----------------------------------------------------------------
 */

/*
 * cvAdjustParams
 */
fn cvAdjustParams(cv_mem: &mut CVodeMem) {
    if cv_mem.cv_qprime != cv_mem.cv_q {
        /* History adjustments for an order change were applied when resizing */
        if !cv_mem.first_step_after_resize {
            cvAdjustOrder(cv_mem, cv_mem.cv_qprime - cv_mem.cv_q);
        }
        cv_mem.cv_q = cv_mem.cv_qprime;
        cv_mem.cv_L = cv_mem.cv_q + 1;
        cv_mem.cv_qwait = cv_mem.cv_L;
    }
    cvRescale(cv_mem);
}

/*
 * cvAdjustOrder
 */
pub(crate) fn cvAdjustOrder(cv_mem: &mut CVodeMem, deltaq: i32) {
    if cv_mem.cv_q == 2 && deltaq != 1 {
        return;
    }

    match cv_mem.cv_lmm {
        CV_ADAMS => cvAdjustAdams(cv_mem, deltaq),
        CV_BDF => cvAdjustBDF(cv_mem, deltaq),
        _ => {}
    }
}

/*
 * cvAdjustAdams
 *
 * Adjusts the history array on a change of order q by deltaq for
 * lmm == CV_ADAMS.
 */
fn cvAdjustAdams(cv_mem: &mut CVodeMem, deltaq: i32) {
    /* On an order increase, set new column of zn to zero and return */
    if deltaq == 1 {
        let l = cv_mem.cv_L as usize;
        N_VConst(ZERO, &mut cv_mem.cv_zn[l]);
        return;
    }

    /* On an order decrease, each zn[j] is adjusted by a multiple of zn[q]. */
    for i in 0..=(cv_mem.cv_qmax as usize) {
        cv_mem.cv_l[i] = ZERO;
    }
    cv_mem.cv_l[1] = ONE;
    let mut hsum = ZERO;
    for j in 1..=(cv_mem.cv_q - 2) {
        hsum += cv_mem.cv_tau[j as usize];
        let xi = hsum / cv_mem.cv_hscale;
        let mut i = j + 1;
        while i >= 1 {
            cv_mem.cv_l[i as usize] = cv_mem.cv_l[i as usize] * xi + cv_mem.cv_l[(i - 1) as usize];
            i -= 1;
        }
    }

    for j in 1..=(cv_mem.cv_q - 2) {
        cv_mem.cv_l[(j + 1) as usize] = cv_mem.cv_q as f64 * (cv_mem.cv_l[j as usize] / (j + 1) as f64);
    }

    /* zn[j] += -l[j]*zn[q] for j = 2..q-1 (N_VScaleAddMulti in C) */
    if cv_mem.cv_q > 2 {
        let q = cv_mem.cv_q as usize;
        let (front, back) = cv_mem.cv_zn.split_at_mut(q);
        let znq = &back[0];
        for j in 2..q {
            let c = -cv_mem.cv_l[j];
            for (z, x) in front[j].data.iter_mut().zip(&znq.data) {
                *z += c * *x;
            }
        }
    }
}

/*
 * cvAdjustBDF
 */
fn cvAdjustBDF(cv_mem: &mut CVodeMem, deltaq: i32) {
    match deltaq {
        1 => cvIncreaseBDF(cv_mem),
        -1 => cvDecreaseBDF(cv_mem),
        _ => {}
    }
}

/*
 * cvIncreaseBDF
 *
 * Adjusts the history array on an increase in order q for CV_BDF.
 */
fn cvIncreaseBDF(cv_mem: &mut CVodeMem) {
    for i in 0..=(cv_mem.cv_qmax as usize) {
        cv_mem.cv_l[i] = ZERO;
    }
    cv_mem.cv_l[2] = ONE;
    let mut alpha1 = ONE;
    let mut prod = ONE;
    let mut xiold = ONE;
    let mut alpha0 = -ONE;
    let mut hsum = cv_mem.cv_hscale;
    if cv_mem.cv_q > 1 {
        for j in 1..cv_mem.cv_q {
            hsum += cv_mem.cv_tau[(j + 1) as usize];
            let xi = hsum / cv_mem.cv_hscale;
            prod *= xi;
            alpha0 -= ONE / (j + 1) as f64;
            alpha1 += ONE / xi;
            let mut i = j + 2;
            while i >= 2 {
                cv_mem.cv_l[i as usize] =
                    cv_mem.cv_l[i as usize] * xiold + cv_mem.cv_l[(i - 1) as usize];
                i -= 1;
            }
            xiold = xi;
        }
    }
    let a1 = (-alpha0 - alpha1) / prod;
    {
        let indx = cv_mem.cv_indx_acor as usize;
        let l = cv_mem.cv_L as usize;
        /* zn[L] = A1 * zn[indx_acor] */
        if indx == l {
            cv_mem.cv_zn[l].scale_inplace(a1);
        } else {
            let (lo, hi) = if indx < l { (indx, l) } else { (l, indx) };
            let (front, back) = cv_mem.cv_zn.split_at_mut(hi);
            let (src, dst) = if indx < l {
                (&front[lo], &mut back[0])
            } else {
                let tmp = &mut front[lo];
                (&back[0], tmp)
            };
            N_VScale(a1, src, dst);
        }
    }

    /* zn[j] += l[j]*zn[L] for j = 2..=q (N_VScaleAddMulti in C) */
    if cv_mem.cv_q > 1 {
        let l = cv_mem.cv_L as usize;
        let q = cv_mem.cv_q as usize;
        let (front, back) = cv_mem.cv_zn.split_at_mut(l);
        let znl = &back[0];
        for j in 2..=q {
            let c = cv_mem.cv_l[j];
            for (z, x) in front[j].data.iter_mut().zip(&znl.data) {
                *z += c * *x;
            }
        }
    }
}

/*
 * cvDecreaseBDF
 *
 * Adjusts the history array on a decrease in order q for CV_BDF.
 */
fn cvDecreaseBDF(cv_mem: &mut CVodeMem) {
    for i in 0..=(cv_mem.cv_qmax as usize) {
        cv_mem.cv_l[i] = ZERO;
    }
    cv_mem.cv_l[2] = ONE;
    let mut hsum = ZERO;
    for j in 1..=(cv_mem.cv_q - 2) {
        hsum += cv_mem.cv_tau[j as usize];
        let xi = hsum / cv_mem.cv_hscale;
        let mut i = j + 2;
        while i >= 2 {
            cv_mem.cv_l[i as usize] = cv_mem.cv_l[i as usize] * xi + cv_mem.cv_l[(i - 1) as usize];
            i -= 1;
        }
    }

    /* zn[j] += -l[j]*zn[q] for j = 2..q-1 */
    if cv_mem.cv_q > 2 {
        let q = cv_mem.cv_q as usize;
        let (front, back) = cv_mem.cv_zn.split_at_mut(q);
        let znq = &back[0];
        for j in 2..q {
            let c = -cv_mem.cv_l[j];
            for (z, x) in front[j].data.iter_mut().zip(&znq.data) {
                *z += c * *x;
            }
        }
    }
}

/*
 * cvRescale
 *
 * Rescales the Nordsieck array by multiplying the jth column zn[j]
 * by eta^j, j = 1, ..., q. Then h is rescaled by eta and hscale is
 * reset to h.
 */
pub fn cvRescale(cv_mem: &mut CVodeMem) {
    /* compute scaling factors and scale columns (N_VScaleVectorArray) */
    let mut factor = cv_mem.cv_eta;
    for j in 1..=(cv_mem.cv_q as usize) {
        cv_mem.cv_zn[j].scale_inplace(factor);
        factor *= cv_mem.cv_eta;
    }

    cv_mem.cv_h = cv_mem.cv_hscale * cv_mem.cv_eta;
    cv_mem.cv_next_h = cv_mem.cv_h;
    cv_mem.cv_hscale = cv_mem.cv_h;
    cv_mem.cv_nscon = 0;
}

/*
 * cvPredict
 *
 * Advances tn by the tentative step size h and computes the
 * predicted array z_n(0) by repeated additions.
 */
fn cvPredict(cv_mem: &mut CVodeMem) {
    cv_mem.cv_tn += cv_mem.cv_h;
    if cv_mem.cv_tstopset && (cv_mem.cv_tn - cv_mem.cv_tstop) * cv_mem.cv_h > ZERO {
        cv_mem.cv_tn = cv_mem.cv_tstop;
    }

    let q = cv_mem.cv_q;
    for k in 1..=q {
        let mut j = q;
        while j >= k {
            /* zn[j-1] += zn[j] */
            let (front, back) = cv_mem.cv_zn.split_at_mut(j as usize);
            let znj = &back[0];
            for (z, x) in front[(j - 1) as usize].data.iter_mut().zip(&znj.data) {
                *z += *x;
            }
            j -= 1;
        }
    }
}

/*
 * cvSet
 *
 * Sets the polynomial l, the test quantity array tq, and related
 * variables rl1, gamma, and gamrat.
 */
fn cvSet(cv_mem: &mut CVodeMem) {
    match cv_mem.cv_lmm {
        CV_ADAMS => cvSetAdams(cv_mem),
        CV_BDF => cvSetBDF(cv_mem),
        _ => {}
    }
    cv_mem.cv_rl1 = ONE / cv_mem.cv_l[1];
    cv_mem.cv_gamma = cv_mem.cv_h * cv_mem.cv_rl1;
    if cv_mem.cv_nst == 0 {
        cv_mem.cv_gammap = cv_mem.cv_gamma;
    }
    cv_mem.cv_gamrat = if cv_mem.cv_nst > 0 {
        cv_mem.cv_gamma / cv_mem.cv_gammap
    } else {
        ONE /* protect x / x != 1.0 */
    };
}

/*
 * cvSetAdams
 */
fn cvSetAdams(cv_mem: &mut CVodeMem) {
    let mut m = [ZERO; L_MAX];
    let mut big_m = [ZERO; 3];

    if cv_mem.cv_q == 1 {
        cv_mem.cv_l[0] = ONE;
        cv_mem.cv_l[1] = ONE;
        cv_mem.cv_tq[1] = ONE;
        cv_mem.cv_tq[5] = ONE;
        cv_mem.cv_tq[2] = HALF;
        cv_mem.cv_tq[3] = ONE / TWELVE;
        cv_mem.cv_tq[4] = cv_mem.cv_nlscoef / cv_mem.cv_tq[2]; /* = 0.1 / tq[2] */
        return;
    }

    let hsum = cvAdamsStart(cv_mem, &mut m);

    big_m[0] = cvAltSum(cv_mem.cv_q - 1, &m, 1);
    big_m[1] = cvAltSum(cv_mem.cv_q - 1, &m, 2);

    cvAdamsFinish(cv_mem, &mut m, &mut big_m, hsum);
}

/*
 * cvAdamsStart
 *
 * Generates in m[] the coefficients of the product polynomial
 * needed for the Adams l and tq coefficients for q > 1.
 */
fn cvAdamsStart(cv_mem: &mut CVodeMem, m: &mut [f64]) -> f64 {
    let mut hsum = cv_mem.cv_h;
    m[0] = ONE;
    for i in 1..=(cv_mem.cv_q as usize) {
        m[i] = ZERO;
    }
    for j in 1..cv_mem.cv_q {
        if j == cv_mem.cv_q - 1 && cv_mem.cv_qwait == 1 {
            let sum = cvAltSum(cv_mem.cv_q - 2, m, 2);
            cv_mem.cv_tq[1] = cv_mem.cv_q as f64 * sum / m[(cv_mem.cv_q - 2) as usize];
        }
        let xi_inv = cv_mem.cv_h / hsum;
        let mut i = j;
        while i >= 1 {
            m[i as usize] += m[(i - 1) as usize] * xi_inv;
            i -= 1;
        }
        hsum += cv_mem.cv_tau[j as usize];
        /* The m[i] are coefficients of product(1 to j) (1 + x/xi_i) */
    }
    hsum
}

/*
 * cvAdamsFinish
 *
 * Completes the calculation of the Adams l and tq.
 */
fn cvAdamsFinish(cv_mem: &mut CVodeMem, m: &mut [f64], big_m: &mut [f64], hsum: f64) {
    let m0_inv = ONE / big_m[0];

    cv_mem.cv_l[0] = ONE;
    for i in 1..=(cv_mem.cv_q as usize) {
        cv_mem.cv_l[i] = m0_inv * (m[i - 1] / i as f64);
    }
    let xi = hsum / cv_mem.cv_h;
    let xi_inv = ONE / xi;

    cv_mem.cv_tq[2] = big_m[1] * m0_inv / xi;
    cv_mem.cv_tq[5] = xi / cv_mem.cv_l[cv_mem.cv_q as usize];

    if cv_mem.cv_qwait == 1 {
        let mut i = cv_mem.cv_q;
        while i >= 1 {
            m[i as usize] += m[(i - 1) as usize] * xi_inv;
            i -= 1;
        }
        big_m[2] = cvAltSum(cv_mem.cv_q, m, 2);
        cv_mem.cv_tq[3] = big_m[2] * m0_inv / cv_mem.cv_L as f64;
    }

    cv_mem.cv_tq[4] = cv_mem.cv_nlscoef / cv_mem.cv_tq[2];
}

/*
 * cvAltSum
 *
 * Returns the value of the alternating sum
 *   sum (i = 0 ... iend) [ (-1)^i * (a[i] / (i + k)) ].
 */
fn cvAltSum(iend: i32, a: &[f64], k: i32) -> f64 {
    if iend < 0 {
        return ZERO;
    }

    let mut sum = ZERO;
    let mut sign = 1i32;
    for i in 0..=iend {
        sum += sign as f64 * (a[i as usize] / (i + k) as f64);
        sign = -sign;
    }
    sum
}

/*
 * cvSetBDF
 */
fn cvSetBDF(cv_mem: &mut CVodeMem) {
    cv_mem.cv_l[0] = ONE;
    cv_mem.cv_l[1] = ONE;
    let mut xi_inv = ONE;
    let mut xistar_inv = ONE;
    for i in 2..=(cv_mem.cv_q as usize) {
        cv_mem.cv_l[i] = ZERO;
    }
    let mut alpha0 = -ONE;
    let mut alpha0_hat = -ONE;
    let mut hsum = cv_mem.cv_h;

    if cv_mem.proj_enabled {
        for i in 0..=(cv_mem.cv_q as usize) {
            cv_mem.proj_p[i] = cv_mem.cv_l[i];
        }
    }

    if cv_mem.cv_q > 1 {
        for j in 2..cv_mem.cv_q {
            hsum += cv_mem.cv_tau[(j - 1) as usize];
            xi_inv = cv_mem.cv_h / hsum;
            alpha0 -= ONE / j as f64;
            let mut i = j;
            while i >= 1 {
                cv_mem.cv_l[i as usize] += cv_mem.cv_l[(i - 1) as usize] * xi_inv;
                i -= 1;
            }
            /* The l[i] are coefficients of product(1 to j) (1 + x/xi_i) */
        }

        /* j = q */
        alpha0 -= ONE / cv_mem.cv_q as f64;
        xistar_inv = -cv_mem.cv_l[1] - alpha0;
        hsum += cv_mem.cv_tau[(cv_mem.cv_q - 1) as usize];
        xi_inv = cv_mem.cv_h / hsum;
        alpha0_hat = -cv_mem.cv_l[1] - xi_inv;

        if cv_mem.proj_enabled {
            let mut i = cv_mem.cv_q;
            while i >= 1 {
                cv_mem.proj_p[i as usize] =
                    cv_mem.cv_l[i as usize] + cv_mem.proj_p[(i - 1) as usize] * xi_inv;
                i -= 1;
            }
        }

        let mut i = cv_mem.cv_q;
        while i >= 1 {
            cv_mem.cv_l[i as usize] += cv_mem.cv_l[(i - 1) as usize] * xistar_inv;
            i -= 1;
        }
    }

    cvSetTqBDF(cv_mem, hsum, alpha0, alpha0_hat, xi_inv, xistar_inv);
}

/*
 * cvSetTqBDF
 */
fn cvSetTqBDF(
    cv_mem: &mut CVodeMem,
    mut hsum: f64,
    alpha0: f64,
    alpha0_hat: f64,
    mut xi_inv: f64,
    xistar_inv: f64,
) {
    let a1 = ONE - alpha0_hat + alpha0;
    let a2 = ONE + cv_mem.cv_q as f64 * a1;
    cv_mem.cv_tq[2] = SUNRabs(a1 / (alpha0 * a2));
    cv_mem.cv_tq[5] = SUNRabs(a2 * xistar_inv / (cv_mem.cv_l[cv_mem.cv_q as usize] * xi_inv));
    if cv_mem.cv_qwait == 1 {
        if cv_mem.cv_q > 1 {
            let c = xistar_inv / cv_mem.cv_l[cv_mem.cv_q as usize];
            let a3 = alpha0 + ONE / cv_mem.cv_q as f64;
            let a4 = alpha0_hat + xi_inv;
            let cpinv = (ONE - a4 + a3) / a3;
            cv_mem.cv_tq[1] = SUNRabs(c * cpinv);
        } else {
            cv_mem.cv_tq[1] = ONE;
        }
        hsum += cv_mem.cv_tau[cv_mem.cv_q as usize];
        xi_inv = cv_mem.cv_h / hsum;
        let a5 = alpha0 - ONE / (cv_mem.cv_q + 1) as f64;
        let a6 = alpha0_hat - xi_inv;
        let cppinv = (ONE - a6 + a5) / a2;
        cv_mem.cv_tq[3] = SUNRabs(cppinv / (xi_inv * (cv_mem.cv_q + 2) as f64 * a5));
    }
    cv_mem.cv_tq[4] = cv_mem.cv_nlscoef / cv_mem.cv_tq[2];
}

/*
 * cvCheckConstraints
 *
 * Determines if the constraints of the problem are satisfied by
 * the proposed step.
 */
fn cvCheckConstraints(
    cv_mem: &mut CVodeMem,
    nflag_ptr: &mut i32,
    saved_t: f64,
    step_constraint_fails: &mut i32,
) -> i32 {
    /* mm = cv_ftemp (mask), tmp = cv_tempv (workspace) */

    /* Get mask vector mm, 1 where constraints failed and 0 otherwise */
    let constraints_passed = {
        let (c, y, mm) = (&cv_mem.cv_constraints, &cv_mem.cv_y, &mut cv_mem.cv_ftemp);
        N_VConstrMask(c, y, mm)
    };
    if constraints_passed {
        return CV_SUCCESS;
    }

    /* Constraints not met */

    /* Compute correction v such that y - v will satisfy the constraints */
    {
        let CVodeMem {
            cv_constraints,
            cv_ewt,
            cv_y,
            cv_tempv: tmp,
            cv_vtemp1,
            ..
        } = cv_mem;
        N_VCompare(ONEPT5, cv_constraints, tmp);
        tmp.prod_with(cv_constraints); /* tmp = tmp * constraints */
        tmp.div_with(cv_ewt); /* tmp = tmp / ewt */
        N_VScale(-PT1, tmp, cv_vtemp1); /* vtemp1 = -0.1*tmp (saved adjustment) */
        /* tmp = y - 0.1*tmp */
        tmp.linear_sum_with(-PT1, ONE, cv_y);
        /* tmp = tmp * mm */
        tmp.prod_with(&cv_mem.cv_ftemp);
    }

    let vnorm = N_VWrmsNorm(&cv_mem.cv_tempv, &cv_mem.cv_ewt); /* ||v|| */

    /* If the correction is small in norm, correct and accept this step */
    if vnorm <= cv_mem.cv_tq[4] {
        /* Update constraint correction count */
        cv_mem.constraint_corrections += 1;

        /* Split the correction update acor = acor - v into three steps */
        {
            let CVodeMem {
                cv_ftemp: mm,
                cv_tempv: tmp,
                cv_acor,
                cv_zn,
                cv_vtemp1,
                ..
            } = cv_mem;

            /* Zero out the correction where any constraint failed */
            N_VProd(mm, cv_acor, tmp);
            cv_acor.linear_sum_with(ONE, -ONE, tmp);

            /* Set correction to zero out the predictor where constraints failed */
            N_VProd(mm, &cv_zn[0], tmp);
            cv_acor.linear_sum_with(ONE, -ONE, tmp);

            /* Update the correction with the adjustment saved above */
            cv_vtemp1.prod_with(mm);
            cv_acor.linear_sum_with(ONE, -ONE, cv_vtemp1);
        }

        return CV_SUCCESS;
    }

    /* update failure counts */
    *step_constraint_fails += 1;
    cv_mem.constraint_fails += 1;

    /* restore zn */
    cvRestore(cv_mem, saved_t);

    /* Check for |h| == hmin */
    if SUNRabs(cv_mem.cv_h) <= cv_mem.cv_hmin * ONEPSM {
        return CV_CONSTR_FAIL;
    }

    /* Check for max step attempt failures */
    if *step_constraint_fails == cv_mem.max_constraint_fails {
        return CV_CONSTR_FAIL;
    }

    /* Constraint correction is too large, reduce h by computing eta = h'/h */
    {
        let CVodeMem {
            cv_ftemp: mm,
            cv_tempv: tmp,
            cv_zn,
            cv_y,
            ..
        } = cv_mem;
        N_VLinearSum(ONE, &cv_zn[0], -ONE, cv_y, tmp);
        tmp.prod_with(mm);
    }

    /* Reduce step size; return to reattempt the step */
    cv_mem.cv_eta = PT9 * N_VMinQuotient(&cv_mem.cv_zn[0], &cv_mem.cv_tempv);
    cv_mem.cv_eta = SUNMAX(cv_mem.cv_eta, PT1);
    cv_mem.cv_eta = SUNMAX(cv_mem.cv_eta, cv_mem.cv_hmin / SUNRabs(cv_mem.cv_h));
    cvRescale(cv_mem);
    *nflag_ptr = PREV_CONV_FAIL;

    PREDICT_AGAIN
}

/*
 * cvHandleNFlag
 *
 * Takes action on the return value nflag returned by cvNls.
 */
fn cvHandleNFlag(cv_mem: &mut CVodeMem, nflag_ptr: &mut i32, saved_t: f64, ncf_ptr: &mut i32) -> i32 {
    let nflag = *nflag_ptr;

    if nflag == CV_SUCCESS {
        return DO_ERROR_TEST;
    }

    /* The nonlinear soln. failed; increment ncfn and restore zn */
    cv_mem.cv_ncfn += 1;
    cvRestore(cv_mem, saved_t);

    /* Return if failed unrecoverably */
    if nflag < 0 {
        if nflag == CV_LSETUP_FAIL {
            return CV_LSETUP_FAIL;
        } else if nflag == CV_LSOLVE_FAIL {
            return CV_LSOLVE_FAIL;
        } else if nflag == CV_RHSFUNC_FAIL {
            return CV_RHSFUNC_FAIL;
        } else {
            return CV_NLS_FAIL;
        }
    }

    /* At this point, a recoverable error occurred. */
    *ncf_ptr += 1;
    cv_mem.cv_etamax = ONE;

    /* If we had maxncf failures or |h| = hmin, return failure. */
    if SUNRabs(cv_mem.cv_h) <= cv_mem.cv_hmin * ONEPSM || *ncf_ptr == cv_mem.cv_maxncf {
        if nflag == SUN_NLS_CONV_RECVR {
            return CV_CONV_FAILURE;
        }
        if nflag == RHSFUNC_RECVR {
            return CV_REPTD_RHSFUNC_ERR;
        }
    }

    /* Reduce step size; return to reattempt the step */
    cv_mem.cv_eta = SUNMAX(cv_mem.cv_eta_cf, cv_mem.cv_hmin / SUNRabs(cv_mem.cv_h));
    *nflag_ptr = PREV_CONV_FAIL;
    cvRescale(cv_mem);

    PREDICT_AGAIN
}

/*
 * cvRestore
 *
 * Restores tn to saved_t and undoes the prediction.
 */
pub fn cvRestore(cv_mem: &mut CVodeMem, saved_t: f64) {
    cv_mem.cv_tn = saved_t;
    let q = cv_mem.cv_q;
    for k in 1..=q {
        let mut j = q;
        while j >= k {
            /* zn[j-1] -= zn[j] */
            let (front, back) = cv_mem.cv_zn.split_at_mut(j as usize);
            let znj = &back[0];
            for (z, x) in front[(j - 1) as usize].data.iter_mut().zip(&znj.data) {
                *z -= *x;
            }
            j -= 1;
        }
    }
}

/*
 * -----------------------------------------------------------------
 * Error Test
 * -----------------------------------------------------------------
 */

/*
 * cvDoErrorTest
 *
 * Performs the local error test.
 */
fn cvDoErrorTest(
    cv_mem: &mut CVodeMem,
    nflag_ptr: &mut i32,
    saved_t: f64,
    nef_ptr: &mut i32,
    dsm_ptr: &mut f64,
) -> i32 {
    let dsm = cv_mem.cv_acnrm * cv_mem.cv_tq[2];

    /* If est. local error norm dsm passes test, return CV_SUCCESS */
    *dsm_ptr = dsm;
    if dsm <= ONE {
        return CV_SUCCESS;
    }

    /* Test failed; increment counters, set nflag, and restore zn array */
    *nef_ptr += 1;
    cv_mem.cv_netf += 1;
    *nflag_ptr = PREV_ERR_FAIL;
    cvRestore(cv_mem, saved_t);

    /* At maxnef failures or |h| = hmin, return CV_ERR_FAILURE */
    if SUNRabs(cv_mem.cv_h) <= cv_mem.cv_hmin * ONEPSM || *nef_ptr == cv_mem.cv_maxnef {
        return CV_ERR_FAILURE;
    }

    /* Set etamax = 1 to prevent step size increase at end of this step */
    cv_mem.cv_etamax = ONE;

    /* Set h ratio eta from dsm, rescale, and return for retry of step */
    if *nef_ptr <= MXNEF1 {
        cv_mem.cv_eta = ONE / (SUNRpowerR(BIAS2 * dsm, ONE / cv_mem.cv_L as f64) + ADDON);
        cv_mem.cv_eta = SUNMAX(
            cv_mem.cv_eta_min_ef,
            SUNMAX(cv_mem.cv_eta, cv_mem.cv_hmin / SUNRabs(cv_mem.cv_h)),
        );
        if *nef_ptr >= cv_mem.cv_small_nef {
            cv_mem.cv_eta = SUNMIN(cv_mem.cv_eta, cv_mem.cv_eta_max_ef);
        }

        cvRescale(cv_mem);
        return TRY_AGAIN;
    }

    /* After MXNEF1 failures, force an order reduction and retry step */
    if cv_mem.cv_q > 1 {
        cv_mem.cv_eta = SUNMAX(cv_mem.cv_eta_min_ef, cv_mem.cv_hmin / SUNRabs(cv_mem.cv_h));
        cvAdjustOrder(cv_mem, -1);
        cv_mem.cv_L = cv_mem.cv_q;
        cv_mem.cv_q -= 1;
        cv_mem.cv_qwait = cv_mem.cv_L;
        cvRescale(cv_mem);
        return TRY_AGAIN;
    }

    /* If already at order 1, restart: reload zn from scratch */
    cv_mem.cv_eta = SUNMAX(cv_mem.cv_eta_min_ef, cv_mem.cv_hmin / SUNRabs(cv_mem.cv_h));
    cv_mem.cv_h *= cv_mem.cv_eta;
    cv_mem.cv_next_h = cv_mem.cv_h;
    cv_mem.cv_hscale = cv_mem.cv_h;
    cv_mem.cv_qwait = LONG_WAIT;
    cv_mem.cv_nscon = 0;

    let f = cv_mem.cv_f.unwrap();
    let retval = f(
        cv_mem.cv_tn,
        &cv_mem.cv_zn[0],
        &mut cv_mem.cv_tempv,
        &mut cv_mem.cv_user_data,
    );
    cv_mem.cv_nfe += 1;
    if retval < 0 {
        return CV_RHSFUNC_FAIL;
    }
    if retval > 0 {
        return CV_UNREC_RHSFUNC_ERR;
    }

    {
        let (tempv, zn) = (&cv_mem.cv_tempv, &mut cv_mem.cv_zn);
        N_VScale(cv_mem.cv_h, tempv, &mut zn[1]);
    }

    TRY_AGAIN
}

/*
 * -----------------------------------------------------------------
 * Functions called after a successful step
 * -----------------------------------------------------------------
 */

/*
 * cvCompleteStep
 *
 * Performs update operations when the solution to the nonlinear
 * system has passed the local error test.
 */
fn cvCompleteStep(cv_mem: &mut CVodeMem) {
    cv_mem.cv_nst += 1;
    cv_mem.cv_nscon += 1;
    cv_mem.cv_hu = cv_mem.cv_h;
    cv_mem.cv_qu = cv_mem.cv_q;

    cv_mem.first_step_after_resize = SUNFALSE;

    let mut i = cv_mem.cv_q;
    while i >= 2 {
        cv_mem.cv_tau[i as usize] = cv_mem.cv_tau[(i - 1) as usize];
        i -= 1;
    }
    if cv_mem.cv_q == 1 && cv_mem.cv_nst > 1 {
        cv_mem.cv_tau[2] = cv_mem.cv_tau[1];
    }
    cv_mem.cv_tau[1] = cv_mem.cv_h;

    /* Apply correction to column j of zn: l_j * Delta_n */
    {
        let CVodeMem { cv_zn, cv_acor, cv_l, .. } = cv_mem;
        for j in 0..=(cv_mem.cv_q as usize) {
            let c = cv_l[j];
            for (z, a) in cv_zn[j].data.iter_mut().zip(&cv_acor.data) {
                *z += c * *a;
            }
        }
    }

    /* Apply the projection correction to column j of zn: p_j * Delta_n */
    if cv_mem.proj_applied {
        let CVodeMem { cv_zn, cv_tempv, proj_p, .. } = cv_mem;
        for j in 0..=(cv_mem.cv_q as usize) {
            let c = proj_p[j];
            for (z, a) in cv_zn[j].data.iter_mut().zip(&cv_tempv.data) {
                /* tempv = acorP */
                *z += c * *a;
            }
        }
    }

    cv_mem.cv_qwait -= 1;
    if cv_mem.cv_qwait == 1 && cv_mem.cv_q != cv_mem.cv_qmax {
        let qmax = cv_mem.cv_qmax as usize;
        let CVodeMem { cv_zn, cv_acor, .. } = cv_mem;
        cv_zn[qmax].data.copy_from_slice(&cv_acor.data);
        cv_mem.cv_saved_tq5 = cv_mem.cv_tq[5];
        cv_mem.cv_indx_acor = cv_mem.cv_qmax;
    }
}

/*
 * cvPrepareNextStep
 *
 * Handles the setting of stepsize and order for the next step.
 */
fn cvPrepareNextStep(cv_mem: &mut CVodeMem, dsm: f64) {
    /* If etamax = 1, defer step size or order changes */
    if cv_mem.cv_etamax == ONE {
        cv_mem.cv_qwait = cv_mem.cv_qwait.max(2);
        cv_mem.cv_qprime = cv_mem.cv_q;
        cv_mem.cv_hprime = cv_mem.cv_h;
        cv_mem.cv_eta = ONE;
    } else {
        /* etaq is the ratio of new to old h at the current order */
        cv_mem.cv_etaq = ONE / (SUNRpowerR(BIAS2 * dsm, ONE / cv_mem.cv_L as f64) + ADDON);

        /* If no order change, adjust eta and acor in cvSetEta and return */
        if cv_mem.cv_qwait != 0 {
            cv_mem.cv_eta = cv_mem.cv_etaq;
            cv_mem.cv_qprime = cv_mem.cv_q;
            cvSetEta(cv_mem);
        } else {
            /* If qwait = 0, consider an order change. */
            cv_mem.cv_qwait = 2;
            cv_mem.cv_etaqm1 = cvComputeEtaqm1(cv_mem);
            cv_mem.cv_etaqp1 = cvComputeEtaqp1(cv_mem);
            cvChooseEta(cv_mem);
            cvSetEta(cv_mem);
        }
    }
}

/*
 * cvSetEta
 *
 * Adjusts the value of eta according to the various heuristic
 * limits and the optional input hmax.
 */
fn cvSetEta(cv_mem: &mut CVodeMem) {
    if cv_mem.cv_eta > cv_mem.cv_eta_min_fx && cv_mem.cv_eta < cv_mem.cv_eta_max_fx {
        /* Eta is within the fixed step bounds, retain step size */
        cv_mem.cv_eta = ONE;
        cv_mem.cv_hprime = cv_mem.cv_h;
    } else {
        if cv_mem.cv_eta >= cv_mem.cv_eta_max_fx {
            /* Increase the step size, limit eta by etamax and hmax */
            cv_mem.cv_eta = SUNMIN(cv_mem.cv_eta, cv_mem.cv_etamax);
            cv_mem.cv_eta /= SUNMAX(
                ONE,
                SUNRabs(cv_mem.cv_h) * cv_mem.cv_hmax_inv * cv_mem.cv_eta,
            );
        } else {
            /* Reduce the step size, limit eta by etamin and hmin */
            cv_mem.cv_eta = SUNMAX(cv_mem.cv_eta, cv_mem.cv_eta_min);
            cv_mem.cv_eta = SUNMAX(cv_mem.cv_eta, cv_mem.cv_hmin / SUNRabs(cv_mem.cv_h));
        }
        /* Set hprime */
        cv_mem.cv_hprime = cv_mem.cv_h * cv_mem.cv_eta;
        if cv_mem.cv_qprime < cv_mem.cv_q {
            cv_mem.cv_nscon = 0;
        }
    }
}

/*
 * cvComputeEtaqm1
 */
fn cvComputeEtaqm1(cv_mem: &mut CVodeMem) -> f64 {
    cv_mem.cv_etaqm1 = ZERO;
    if cv_mem.cv_q > 1 {
        let ddn = N_VWrmsNorm(&cv_mem.cv_zn[cv_mem.cv_q as usize], &cv_mem.cv_ewt) * cv_mem.cv_tq[1];
        cv_mem.cv_etaqm1 = ONE / (SUNRpowerR(BIAS1 * ddn, ONE / cv_mem.cv_q as f64) + ADDON);
    }
    cv_mem.cv_etaqm1
}

/*
 * cvComputeEtaqp1
 */
fn cvComputeEtaqp1(cv_mem: &mut CVodeMem) -> f64 {
    cv_mem.cv_etaqp1 = ZERO;
    if cv_mem.cv_q != cv_mem.cv_qmax {
        if cv_mem.cv_saved_tq5 == ZERO {
            return cv_mem.cv_etaqp1;
        }
        let cquot = (cv_mem.cv_tq[5] / cv_mem.cv_saved_tq5)
            * SUNRpowerI(cv_mem.cv_h / cv_mem.cv_tau[2], cv_mem.cv_L);
        {
            let CVodeMem { cv_zn, cv_acor, cv_tempv, .. } = cv_mem;
            let qmax = cv_mem.cv_qmax as usize;
            N_VLinearSum(-cquot, &cv_zn[qmax], ONE, cv_acor, cv_tempv);
        }
        let dup = N_VWrmsNorm(&cv_mem.cv_tempv, &cv_mem.cv_ewt) * cv_mem.cv_tq[3];
        cv_mem.cv_etaqp1 = ONE / (SUNRpowerR(BIAS3 * dup, ONE / (cv_mem.cv_L + 1) as f64) + ADDON);
    }
    cv_mem.cv_etaqp1
}

/*
 * cvChooseEta
 *
 * Chooses the maximum eta value among etaqm1/etaq/etaqp1 and sets
 * qprime correspondingly.
 */
fn cvChooseEta(cv_mem: &mut CVodeMem) {
    let etam = SUNMAX(cv_mem.cv_etaqm1, SUNMAX(cv_mem.cv_etaq, cv_mem.cv_etaqp1));

    if etam > cv_mem.cv_eta_min_fx && etam < cv_mem.cv_eta_max_fx {
        cv_mem.cv_eta = ONE;
        cv_mem.cv_qprime = cv_mem.cv_q;
    } else if etam == cv_mem.cv_etaq {
        cv_mem.cv_eta = cv_mem.cv_etaq;
        cv_mem.cv_qprime = cv_mem.cv_q;
    } else if etam == cv_mem.cv_etaqm1 {
        cv_mem.cv_eta = cv_mem.cv_etaqm1;
        cv_mem.cv_qprime = cv_mem.cv_q - 1;
    } else {
        cv_mem.cv_eta = cv_mem.cv_etaqp1;
        cv_mem.cv_qprime = cv_mem.cv_q + 1;

        if cv_mem.cv_lmm == CV_BDF {
            /* Store Delta_n in zn[qmax] to be used in order increase */
            let qmax = cv_mem.cv_qmax as usize;
            let CVodeMem { cv_zn, cv_acor, .. } = cv_mem;
            cv_zn[qmax].data.copy_from_slice(&cv_acor.data);
        }
    }
}

/*
 * -----------------------------------------------------------------
 * Function to handle failures
 * -----------------------------------------------------------------
 */

/*
 * cvHandleFailure
 *
 * Prints error messages for all cases of failure by cvHin and
 * cvStep. Returns the value to be returned to the user.
 */
fn cvHandleFailure(cv_mem: &mut CVodeMem, flag: i32) -> i32 {
    match flag {
        CV_ERR_FAILURE => cvProcessError(Some(cv_mem), CV_ERR_FAILURE, line!(), "CVode", file!(),
            &format!("At t = {} and h = {}, the error test failed repeatedly or with |h| = hmin.",
                     cv_mem.cv_tn, cv_mem.cv_h)),
        CV_CONV_FAILURE => cvProcessError(Some(cv_mem), CV_CONV_FAILURE, line!(), "CVode", file!(),
            &format!("At t = {} and h = {}, the corrector convergence test failed repeatedly or with |h| = hmin.",
                     cv_mem.cv_tn, cv_mem.cv_h)),
        CV_LSETUP_FAIL => cvProcessError(Some(cv_mem), CV_LSETUP_FAIL, line!(), "CVode", file!(),
            &format!("At t = {}, the setup routine failed in an unrecoverable manner.", cv_mem.cv_tn)),
        CV_LSOLVE_FAIL => cvProcessError(Some(cv_mem), CV_LSOLVE_FAIL, line!(), "CVode", file!(),
            &format!("At t = {}, the solve routine failed in an unrecoverable manner.", cv_mem.cv_tn)),
        CV_RHSFUNC_FAIL => cvProcessError(Some(cv_mem), CV_RHSFUNC_FAIL, line!(), "CVode", file!(),
            &format!("At t = {}, the right-hand side routine failed in an unrecoverable manner.", cv_mem.cv_tn)),
        CV_UNREC_RHSFUNC_ERR => cvProcessError(Some(cv_mem), CV_UNREC_RHSFUNC_ERR, line!(), "CVode", file!(),
            &format!("At t = {}, the right-hand side failed in a recoverable manner, but no recovery is possible.",
                     cv_mem.cv_tn)),
        CV_REPTD_RHSFUNC_ERR => cvProcessError(Some(cv_mem), CV_REPTD_RHSFUNC_ERR, line!(), "CVode", file!(),
            &format!("At t = {} repeated recoverable right-hand side function errors.", cv_mem.cv_tn)),
        CV_RTFUNC_FAIL => cvProcessError(Some(cv_mem), CV_RTFUNC_FAIL, line!(), "CVode", file!(),
            &format!("At t = {}, the rootfinding routine failed in an unrecoverable manner.", cv_mem.cv_tn)),
        CV_TOO_CLOSE => cvProcessError(Some(cv_mem), CV_TOO_CLOSE, line!(), "CVode", file!(), MSGCV_TOO_CLOSE),
        CV_MEM_NULL => cvProcessError(None, CV_MEM_NULL, line!(), "CVode", file!(), MSGCV_NO_MEM),
        CV_NLS_SETUP_FAIL => cvProcessError(Some(cv_mem), CV_NLS_SETUP_FAIL, line!(), "CVode", file!(),
            &format!("At t = {}, the nonlinear solver setup failed unrecoverably.", cv_mem.cv_tn)),
        CV_CONSTR_FAIL => cvProcessError(Some(cv_mem), CV_CONSTR_FAIL, line!(), "CVode", file!(),
            &format!("At t = {}, unable to satisfy inequality constraints.", cv_mem.cv_tn)),
        CV_NLS_FAIL => cvProcessError(Some(cv_mem), CV_NLS_FAIL, line!(), "CVode", file!(),
            &format!("At t = {}, the nonlinear solver failed in an unrecoverable manner.", cv_mem.cv_tn)),
        CV_PROJ_MEM_NULL => cvProcessError(Some(cv_mem), CV_PROJ_MEM_NULL, line!(), "CVode", file!(),
            "proj_mem = NULL illegal."),
        CV_PROJFUNC_FAIL => cvProcessError(Some(cv_mem), CV_PROJFUNC_FAIL, line!(), "CVode", file!(),
            &format!("At t = {} the projection function failed with an unrecoverable error.", cv_mem.cv_tn)),
        CV_REPTD_PROJFUNC_ERR => cvProcessError(Some(cv_mem), CV_REPTD_PROJFUNC_ERR, line!(), "CVode", file!(),
            &format!("At t = {} the projection function had repeated recoverable errors.", cv_mem.cv_tn)),
        _ => {
            cvProcessError(Some(cv_mem), CV_UNRECOGNIZED_ERR, line!(), "CVode", file!(),
                "CVODE encountered an unrecognized error. Please report this to the SUNDIALS developers at sundials-users@llnl.gov");
            return CV_UNRECOGNIZED_ERR;
        }
    }

    flag
}

/*
 * -----------------------------------------------------------------
 * Functions for BDF Stability Limit Detection
 * -----------------------------------------------------------------
 */

/*
 * cvBDFStab
 *
 * Handles the BDF Stability Limit Detection Algorithm STALD.
 */
fn cvBDFStab(cv_mem: &mut CVodeMem) {
    /* If order is 3 or greater, then save scaled derivative data,
       push old data down in i, then add current values to top. */
    if cv_mem.cv_q >= 3 {
        for k in 1..=3usize {
            let mut i = 5usize;
            while i >= 2 {
                cv_mem.cv_ssdat[i][k] = cv_mem.cv_ssdat[i - 1][k];
                i -= 1;
            }
        }
        let mut factorial = 1i64;
        for i in 1..=(cv_mem.cv_q - 1) {
            factorial *= i as i64;
        }
        let sq = factorial as f64 * cv_mem.cv_q as f64 * (cv_mem.cv_q + 1) as f64
            * cv_mem.cv_acnrm
            / SUNMAX(cv_mem.cv_tq[5], TINY);
        let sqm1 = factorial as f64
            * cv_mem.cv_q as f64
            * N_VWrmsNorm(&cv_mem.cv_zn[cv_mem.cv_q as usize], &cv_mem.cv_ewt);
        let sqm2 = factorial as f64
            * N_VWrmsNorm(&cv_mem.cv_zn[(cv_mem.cv_q - 1) as usize], &cv_mem.cv_ewt);
        cv_mem.cv_ssdat[1][1] = sqm2 * sqm2;
        cv_mem.cv_ssdat[1][2] = sqm1 * sqm1;
        cv_mem.cv_ssdat[1][3] = sq * sq;
    }

    if cv_mem.cv_qprime >= cv_mem.cv_q {
        /* If order is 3 or greater, and enough ssdat has been saved,
           nscon >= q+5, then call stability limit detection routine. */
        if cv_mem.cv_q >= 3 && cv_mem.cv_nscon >= cv_mem.cv_q + 5 {
            let ldflag = cvSLdet(cv_mem);
            if ldflag > 3 {
                /* A stability limit violation is indicated by a return
                   flag of 4, 5, or 6. Reduce new order. */
                cv_mem.cv_qprime = cv_mem.cv_q - 1;
                cv_mem.cv_eta = cv_mem.cv_etaqm1;
                cv_mem.cv_eta = SUNMIN(cv_mem.cv_eta, cv_mem.cv_etamax);
                cv_mem.cv_eta = cv_mem.cv_eta
                    / SUNMAX(
                        ONE,
                        SUNRabs(cv_mem.cv_h) * cv_mem.cv_hmax_inv * cv_mem.cv_eta,
                    );
                cv_mem.cv_hprime = cv_mem.cv_h * cv_mem.cv_eta;
                cv_mem.cv_nor += 1;
            }
        }
    } else {
        /* Otherwise, let order increase happen and reset nscon. */
        cv_mem.cv_nscon = 0;
    }
}

/*
 * cvSLdet
 *
 * Detects stability limitation using stored scaled derivatives data.
 */
fn cvSLdet(cv_mem: &mut CVodeMem) -> i32 {
    let mut rat = [[ZERO; 4]; 5];
    let mut rav = [ZERO; 4];
    let mut qkr = [ZERO; 4];
    let mut sigsq = [ZERO; 4];
    let mut smax = [ZERO; 4];
    let mut ssmax = [ZERO; 4];
    let mut drr = [ZERO; 4];
    let mut rrc = [ZERO; 4];
    let mut sqmx = [ZERO; 4];
    let mut qjk = [[ZERO; 4]; 4];
    let mut vrat = [ZERO; 5];
    let mut qc = [[ZERO; 4]; 6];
    let mut qco = [[ZERO; 4]; 6];

    let mut kmin = 0usize;
    let mut kflag = 0;

    /* Cutoffs and tolerances used by this routine */
    let rrcut = 0.98;
    let vrrtol = 1.0e-4;
    let vrrt2 = 5.0e-4;
    let sqtol = 1.0e-3;
    let rrtol = 1.0e-2;

    let mut rr; /* (C initializes rr = ZERO; every reachable path assigns it) */

    /* get maxima, minima, and variances, and form quartic coefficients */
    for k in 1..=3usize {
        let mut smink = cv_mem.cv_ssdat[1][k];
        let mut smaxk = ZERO;

        for i in 1..=5usize {
            smink = SUNMIN(smink, cv_mem.cv_ssdat[i][k]);
            smaxk = SUNMAX(smaxk, cv_mem.cv_ssdat[i][k]);
        }

        if smink < TINY * smaxk {
            return -1;
        }
        smax[k] = smaxk;
        ssmax[k] = smaxk * smaxk;

        let mut sumrat = ZERO;
        let mut sumrsq = ZERO;
        for i in 1..=4usize {
            rat[i][k] = cv_mem.cv_ssdat[i][k] / cv_mem.cv_ssdat[i + 1][k];
            sumrat += rat[i][k];
            sumrsq += rat[i][k] * rat[i][k];
        }
        rav[k] = FOURTH * sumrat;
        vrat[k] = SUNRabs(FOURTH * sumrsq - rav[k] * rav[k]);

        qc[5][k] = cv_mem.cv_ssdat[1][k] * cv_mem.cv_ssdat[3][k]
            - cv_mem.cv_ssdat[2][k] * cv_mem.cv_ssdat[2][k];
        qc[4][k] = cv_mem.cv_ssdat[2][k] * cv_mem.cv_ssdat[3][k]
            - cv_mem.cv_ssdat[1][k] * cv_mem.cv_ssdat[4][k];
        qc[3][k] = ZERO;
        qc[2][k] = cv_mem.cv_ssdat[2][k] * cv_mem.cv_ssdat[5][k]
            - cv_mem.cv_ssdat[3][k] * cv_mem.cv_ssdat[4][k];
        qc[1][k] = cv_mem.cv_ssdat[4][k] * cv_mem.cv_ssdat[4][k]
            - cv_mem.cv_ssdat[3][k] * cv_mem.cv_ssdat[5][k];

        for i in 1..=5usize {
            qco[i][k] = qc[i][k];
        }
    } /* End of k loop */

    /* Isolate normal or nearly-normal matrix case. */
    let vmin = SUNMIN(vrat[1], SUNMIN(vrat[2], vrat[3]));
    let vmax = SUNMAX(vrat[1], SUNMAX(vrat[2], vrat[3]));

    if vmin < vrrtol * vrrtol {
        if vmax > vrrt2 * vrrt2 {
            return -2;
        }

        rr = (rav[1] + rav[2] + rav[3]) / THREE;
        let mut drrmax = ZERO;
        for k in 1..=3usize {
            let adrr = SUNRabs(rav[k] - rr);
            drrmax = SUNMAX(drrmax, adrr);
        }
        if drrmax > vrrt2 {
            return -3;
        }

        kflag = 1;
        /* can compute characteristic root, drop to next section */
    } else {
        /* use the quartics to get rr. */
        if SUNRabs(qco[1][1]) < TINY * ssmax[1] {
            return -4;
        }

        let mut tem = qco[1][2] / qco[1][1];
        for i in 2..=5usize {
            qco[i][2] -= tem * qco[i][1];
        }

        qco[1][2] = ZERO;
        tem = qco[1][3] / qco[1][1];
        for i in 2..=5usize {
            qco[i][3] -= tem * qco[i][1];
        }
        qco[1][3] = ZERO;

        if SUNRabs(qco[2][2]) < TINY * ssmax[2] {
            return -4;
        }

        tem = qco[2][3] / qco[2][2];
        for i in 3..=5usize {
            qco[i][3] -= tem * qco[i][2];
        }

        if SUNRabs(qco[4][3]) < TINY * ssmax[3] {
            return -4;
        }

        rr = -qco[5][3] / qco[4][3];

        if rr < TINY || rr > HUNDRED {
            return -5;
        }

        for k in 1..=3usize {
            qkr[k] = qc[5][k] + rr * (qc[4][k] + rr * rr * (qc[2][k] + rr * qc[1][k]));
        }

        let mut sqmax = ZERO;
        for k in 1..=3usize {
            let saqk = SUNRabs(qkr[k]) / ssmax[k];
            if saqk > sqmax {
                sqmax = saqk;
            }
        }

        if sqmax < sqtol {
            kflag = 2;
            /* can compute characteristic root, drop to "given rr,etc" */
        } else {
            /* do Newton corrections to improve rr. */
            let mut sqmin = ZERO;
            for _it in 1..=3 {
                for k in 1..=3usize {
                    let qp = qc[4][k] + rr * rr * (THREE * qc[2][k] + rr * FOUR * qc[1][k]);
                    drr[k] = ZERO;
                    if SUNRabs(qp) > TINY * ssmax[k] {
                        drr[k] = -qkr[k] / qp;
                    }
                    rrc[k] = rr + drr[k];
                }

                for k in 1..=3usize {
                    let s = rrc[k];
                    let mut sqmaxk = ZERO;
                    for j in 1..=3usize {
                        qjk[j][k] = qc[5][j] + s * (qc[4][j] + s * s * (qc[2][j] + s * qc[1][j]));
                        let saqj = SUNRabs(qjk[j][k]) / ssmax[j];
                        if saqj > sqmaxk {
                            sqmaxk = saqj;
                        }
                    }
                    sqmx[k] = sqmaxk;
                }

                sqmin = sqmx[1] + ONE;
                for k in 1..=3usize {
                    if sqmx[k] < sqmin {
                        kmin = k;
                        sqmin = sqmx[k];
                    }
                }
                rr = rrc[kmin];

                if sqmin < sqtol {
                    kflag = 3;
                    /* can compute characteristic root; break out */
                    break;
                } else {
                    for j in 1..=3usize {
                        qkr[j] = qjk[j][kmin];
                    }
                }
            } /* end of Newton correction loop */

            if sqmin > sqtol {
                return -6;
            }
        } /* end of if (sqmax < sqtol) else */
    } /* end of if (vmin < vrrtol*vrrtol) else */

    /* given rr, find sigsq[k] and verify rr. */
    /* All positive kflag drop to this section */
    for k in 1..=3usize {
        let rsa = cv_mem.cv_ssdat[1][k];
        let rsb = cv_mem.cv_ssdat[2][k] * rr;
        let rsc = cv_mem.cv_ssdat[3][k] * rr * rr;
        let rsd = cv_mem.cv_ssdat[4][k] * rr * rr * rr;
        let rd1a = rsa - rsb;
        let rd1b = rsb - rsc;
        let rd1c = rsc - rsd;
        let rd2a = rd1a - rd1b;
        let rd2b = rd1b - rd1c;
        let rd3a = rd2a - rd2b;

        if SUNRabs(rd1b) < TINY * smax[k] {
            return -7;
        }

        let cest1 = -rd3a / rd1b;
        if cest1 < TINY || cest1 > FOUR {
            return -7;
        }
        let corr1 = (rd2b / cest1) / (rr * rr);
        sigsq[k] = cv_mem.cv_ssdat[3][k] + corr1;
    }

    if sigsq[2] < TINY {
        return -8;
    }

    let ratp = sigsq[3] / sigsq[2];
    let ratm = sigsq[1] / sigsq[2];
    let qfac1 = FOURTH * (cv_mem.cv_q as f64 * cv_mem.cv_q as f64 - ONE);
    let qfac2 = TWO / (cv_mem.cv_q as f64 - ONE);
    let bb = ratp * ratm - ONE - qfac1 * ratp;
    let tem = ONE - qfac2 * bb;

    if SUNRabs(tem) < TINY {
        return -8;
    }

    let rrb = ONE / tem;

    if SUNRabs(rrb - rr) > rrtol {
        return -9;
    }

    /* Check to see if rr is above cutoff rrcut */
    if rr > rrcut {
        if kflag == 1 {
            kflag = 4;
        }
        if kflag == 2 {
            kflag = 5;
        }
        if kflag == 3 {
            kflag = 6;
        }
    }

    /* All positive kflag returned at this point */
    kflag
}

/*
 * -----------------------------------------------------------------
 * Functions for rootfinding
 * -----------------------------------------------------------------
 */

/*
 * cvRcheck1
 *
 * Completes the initialization of rootfinding memory information and
 * checks whether g has a zero both at and very near the initial
 * point of the IVP.
 */
fn cvRcheck1(cv_mem: &mut CVodeMem) -> i32 {
    for i in 0..(cv_mem.cv_nrtfn as usize) {
        cv_mem.cv_iroots[i] = 0;
    }
    cv_mem.cv_tlo = cv_mem.cv_tn;
    cv_mem.cv_ttol =
        (SUNRabs(cv_mem.cv_tn) + SUNRabs(cv_mem.cv_h)) * cv_mem.cv_uround * HUNDRED;

    /* Evaluate g at initial t and check for zero values. */
    let gfun = cv_mem.cv_gfun.unwrap();
    let retval = gfun(
        cv_mem.cv_tlo,
        &cv_mem.cv_zn[0],
        &mut cv_mem.cv_glo,
        &mut cv_mem.cv_user_data,
    );
    cv_mem.cv_nge = 1;
    if retval != 0 {
        return CV_RTFUNC_FAIL;
    }

    let mut zroot = SUNFALSE;
    for i in 0..(cv_mem.cv_nrtfn as usize) {
        if SUNRabs(cv_mem.cv_glo[i]) == ZERO {
            zroot = SUNTRUE;
            cv_mem.cv_gactive[i] = SUNFALSE;
        }
    }
    if !zroot {
        return CV_SUCCESS;
    }

    /* Some g_i is zero at t0; look at g at t0+(small increment). */
    let hratio = SUNMAX(cv_mem.cv_ttol / SUNRabs(cv_mem.cv_h), PT1);
    let smallh = hratio * cv_mem.cv_h;
    let tplus = cv_mem.cv_tlo + smallh;
    {
        let CVodeMem { cv_zn, cv_y, .. } = cv_mem;
        N_VLinearSum(ONE, &cv_zn[0], hratio, &cv_zn[1], cv_y);
    }
    let retval = gfun(
        tplus,
        &cv_mem.cv_y,
        &mut cv_mem.cv_ghi,
        &mut cv_mem.cv_user_data,
    );
    cv_mem.cv_nge += 1;
    if retval != 0 {
        return CV_RTFUNC_FAIL;
    }

    /* Check now only the components of g which were exactly 0.0 at t0
       to see if we can 'activate' them. */
    for i in 0..(cv_mem.cv_nrtfn as usize) {
        if !cv_mem.cv_gactive[i] && SUNRabs(cv_mem.cv_ghi[i]) != ZERO {
            cv_mem.cv_gactive[i] = SUNTRUE;
            cv_mem.cv_glo[i] = cv_mem.cv_ghi[i];
        }
    }
    CV_SUCCESS
}

/*
 * cvRcheck2
 *
 * Checks for exact zeros of g at the last root found, if the last
 * return was a root.
 */
fn cvRcheck2(cv_mem: &mut CVodeMem) -> i32 {
    if cv_mem.cv_irfnd == 0 {
        return CV_SUCCESS;
    }

    let _ = cvGetDky_into_y(cv_mem, cv_mem.cv_tlo, 0);
    let gfun = cv_mem.cv_gfun.unwrap();
    let retval = gfun(
        cv_mem.cv_tlo,
        &cv_mem.cv_y,
        &mut cv_mem.cv_glo,
        &mut cv_mem.cv_user_data,
    );
    cv_mem.cv_nge += 1;
    if retval != 0 {
        return CV_RTFUNC_FAIL;
    }

    let mut zroot = SUNFALSE;
    for i in 0..(cv_mem.cv_nrtfn as usize) {
        cv_mem.cv_iroots[i] = 0;
    }
    for i in 0..(cv_mem.cv_nrtfn as usize) {
        if !cv_mem.cv_gactive[i] {
            continue;
        }
        if SUNRabs(cv_mem.cv_glo[i]) == ZERO {
            zroot = SUNTRUE;
            cv_mem.cv_iroots[i] = 1;
        }
    }
    if !zroot {
        return CV_SUCCESS;
    }

    /* One or more g_i has a zero at tlo. Check g at tlo+smallh. */
    cv_mem.cv_ttol =
        (SUNRabs(cv_mem.cv_tn) + SUNRabs(cv_mem.cv_h)) * cv_mem.cv_uround * HUNDRED;
    let smallh = if cv_mem.cv_h > ZERO {
        cv_mem.cv_ttol
    } else {
        -cv_mem.cv_ttol
    };
    let tplus = cv_mem.cv_tlo + smallh;
    if (tplus - cv_mem.cv_tn) * cv_mem.cv_h >= ZERO {
        let hratio = smallh / cv_mem.cv_h;
        let CVodeMem { cv_zn, cv_y, .. } = cv_mem;
        cv_y.linear_sum_with(ONE, hratio, &cv_zn[1]);
    } else {
        let _ = cvGetDky_into_y(cv_mem, tplus, 0);
    }
    let retval = gfun(
        tplus,
        &cv_mem.cv_y,
        &mut cv_mem.cv_ghi,
        &mut cv_mem.cv_user_data,
    );
    cv_mem.cv_nge += 1;
    if retval != 0 {
        return CV_RTFUNC_FAIL;
    }

    /* Check for close roots, for a new zero at tlo+smallh, and for a
       g_i that changed from zero to nonzero. */
    let mut zroot = SUNFALSE;
    for i in 0..(cv_mem.cv_nrtfn as usize) {
        if !cv_mem.cv_gactive[i] {
            continue;
        }
        if SUNRabs(cv_mem.cv_ghi[i]) == ZERO {
            if cv_mem.cv_iroots[i] == 1 {
                return CLOSERT;
            }
            zroot = SUNTRUE;
            cv_mem.cv_iroots[i] = 1;
        } else if cv_mem.cv_iroots[i] == 1 {
            cv_mem.cv_glo[i] = cv_mem.cv_ghi[i];
        }
    }
    if zroot {
        return RTFOUND;
    }
    CV_SUCCESS
}

/*
 * cvRcheck3
 *
 * Interfaces to cvRootfind to look for a root of g between tlo and
 * either tn or tout, whichever comes first.
 */
fn cvRcheck3(cv_mem: &mut CVodeMem, tout: f64, itask: i32) -> i32 {
    /* Set thi = tn or tout, whichever comes first; set y = y(thi). */
    if itask == CV_ONE_STEP {
        cv_mem.cv_thi = cv_mem.cv_tn;
        let CVodeMem { cv_zn, cv_y, .. } = cv_mem;
        cv_y.data.copy_from_slice(&cv_zn[0].data);
    }
    if itask == CV_NORMAL {
        if (tout - cv_mem.cv_tn) * cv_mem.cv_h >= ZERO {
            cv_mem.cv_thi = cv_mem.cv_tn;
            let CVodeMem { cv_zn, cv_y, .. } = cv_mem;
            cv_y.data.copy_from_slice(&cv_zn[0].data);
        } else {
            cv_mem.cv_thi = tout;
            let _ = cvGetDky_into_y(cv_mem, cv_mem.cv_thi, 0);
        }
    }

    /* Set ghi = g(thi) and call cvRootfind. */
    let gfun = cv_mem.cv_gfun.unwrap();
    let retval = gfun(
        cv_mem.cv_thi,
        &cv_mem.cv_y,
        &mut cv_mem.cv_ghi,
        &mut cv_mem.cv_user_data,
    );
    cv_mem.cv_nge += 1;
    if retval != 0 {
        return CV_RTFUNC_FAIL;
    }

    cv_mem.cv_ttol =
        (SUNRabs(cv_mem.cv_tn) + SUNRabs(cv_mem.cv_h)) * cv_mem.cv_uround * HUNDRED;
    let ier = cvRootfind(cv_mem);
    if ier == CV_RTFUNC_FAIL {
        return CV_RTFUNC_FAIL;
    }
    for i in 0..(cv_mem.cv_nrtfn as usize) {
        if !cv_mem.cv_gactive[i] && cv_mem.cv_grout[i] != ZERO {
            cv_mem.cv_gactive[i] = SUNTRUE;
        }
    }
    cv_mem.cv_tlo = cv_mem.cv_trout;
    for i in 0..(cv_mem.cv_nrtfn as usize) {
        cv_mem.cv_glo[i] = cv_mem.cv_grout[i];
    }

    /* If no root found, return CV_SUCCESS. */
    if ier == CV_SUCCESS {
        return CV_SUCCESS;
    }

    /* If a root was found, interpolate to get y(trout) and return. */
    let _ = cvGetDky_into_y(cv_mem, cv_mem.cv_trout, 0);
    RTFOUND
}

/*
 * cvRootfind
 *
 * Solves for a root of g(t) between tlo and thi, if one exists,
 * using the Illinois algorithm (a modified secant method).
 */
fn cvRootfind(cv_mem: &mut CVodeMem) -> i32 {
    let nrt = cv_mem.cv_nrtfn as usize;
    let mut imax = 0usize;

    /* First check for change in sign in ghi or for a zero in ghi. */
    let mut maxfrac = ZERO;
    let mut zroot = SUNFALSE;
    let mut sgnchg = SUNFALSE;
    for i in 0..nrt {
        if !cv_mem.cv_gactive[i] {
            continue;
        }
        if SUNRabs(cv_mem.cv_ghi[i]) == ZERO {
            if cv_mem.cv_rootdir[i] as f64 * cv_mem.cv_glo[i] <= ZERO {
                zroot = SUNTRUE;
            }
        } else if SUNRdifferentsign(cv_mem.cv_glo[i], cv_mem.cv_ghi[i])
            && cv_mem.cv_rootdir[i] as f64 * cv_mem.cv_glo[i] <= ZERO
        {
            let gfrac = SUNRabs(cv_mem.cv_ghi[i] / (cv_mem.cv_ghi[i] - cv_mem.cv_glo[i]));
            if gfrac > maxfrac {
                sgnchg = SUNTRUE;
                maxfrac = gfrac;
                imax = i;
            }
        }
    }

    /* If no sign change was found, reset trout and grout. */
    if !sgnchg {
        cv_mem.cv_trout = cv_mem.cv_thi;
        for i in 0..nrt {
            cv_mem.cv_grout[i] = cv_mem.cv_ghi[i];
        }
        if !zroot {
            return CV_SUCCESS;
        }
        for i in 0..nrt {
            cv_mem.cv_iroots[i] = 0;
            if !cv_mem.cv_gactive[i] {
                continue;
            }
            if SUNRabs(cv_mem.cv_ghi[i]) == ZERO
                && cv_mem.cv_rootdir[i] as f64 * cv_mem.cv_glo[i] <= ZERO
            {
                cv_mem.cv_iroots[i] = if cv_mem.cv_glo[i] > 0.0 { -1 } else { 1 };
            }
        }
        return RTFOUND;
    }

    /* Initialize alph to avoid compiler warning */
    let mut alph = ONE;

    /* A sign change was found. Loop to locate nearest root. */
    let mut side = 0;
    let mut sideprev = -1;
    let gfun = cv_mem.cv_gfun.unwrap();
    loop {
        /* If interval size is already less than tolerance ttol, break. */
        if SUNRabs(cv_mem.cv_thi - cv_mem.cv_tlo) <= cv_mem.cv_ttol {
            break;
        }

        /* Set weight alph. */
        if sideprev == side {
            alph = if side == 2 { alph * TWO } else { alph * HALF };
        } else {
            alph = ONE;
        }

        /* Set next root approximation tmid and get g(tmid). */
        let mut tmid = cv_mem.cv_thi
            - (cv_mem.cv_thi - cv_mem.cv_tlo) * cv_mem.cv_ghi[imax]
                / (cv_mem.cv_ghi[imax] - alph * cv_mem.cv_glo[imax]);
        if SUNRabs(tmid - cv_mem.cv_tlo) < HALF * cv_mem.cv_ttol {
            let fracint = SUNRabs(cv_mem.cv_thi - cv_mem.cv_tlo) / cv_mem.cv_ttol;
            let fracsub = if fracint > FIVE { PT1 } else { HALF / fracint };
            tmid = cv_mem.cv_tlo + fracsub * (cv_mem.cv_thi - cv_mem.cv_tlo);
        }
        if SUNRabs(cv_mem.cv_thi - tmid) < HALF * cv_mem.cv_ttol {
            let fracint = SUNRabs(cv_mem.cv_thi - cv_mem.cv_tlo) / cv_mem.cv_ttol;
            let fracsub = if fracint > FIVE { PT1 } else { HALF / fracint };
            tmid = cv_mem.cv_thi - fracsub * (cv_mem.cv_thi - cv_mem.cv_tlo);
        }

        let _ = cvGetDky_into_y(cv_mem, tmid, 0);
        let retval = gfun(
            tmid,
            &cv_mem.cv_y,
            &mut cv_mem.cv_grout,
            &mut cv_mem.cv_user_data,
        );
        cv_mem.cv_nge += 1;
        if retval != 0 {
            return CV_RTFUNC_FAIL;
        }

        /* Check to see in which subinterval g changes sign, and reset imax. */
        maxfrac = ZERO;
        zroot = SUNFALSE;
        sgnchg = SUNFALSE;
        sideprev = side;
        for i in 0..nrt {
            if !cv_mem.cv_gactive[i] {
                continue;
            }
            if SUNRabs(cv_mem.cv_grout[i]) == ZERO {
                if cv_mem.cv_rootdir[i] as f64 * cv_mem.cv_glo[i] <= ZERO {
                    zroot = SUNTRUE;
                }
            } else if SUNRdifferentsign(cv_mem.cv_glo[i], cv_mem.cv_grout[i])
                && cv_mem.cv_rootdir[i] as f64 * cv_mem.cv_glo[i] <= ZERO
            {
                let gfrac = SUNRabs(cv_mem.cv_grout[i] / (cv_mem.cv_grout[i] - cv_mem.cv_glo[i]));
                if gfrac > maxfrac {
                    sgnchg = SUNTRUE;
                    maxfrac = gfrac;
                    imax = i;
                }
            }
        }
        if sgnchg {
            /* Sign change found in (tlo,tmid); replace thi with tmid. */
            cv_mem.cv_thi = tmid;
            for i in 0..nrt {
                cv_mem.cv_ghi[i] = cv_mem.cv_grout[i];
            }
            side = 1;
            /* Stop at root thi if converged; otherwise loop. */
            if SUNRabs(cv_mem.cv_thi - cv_mem.cv_tlo) <= cv_mem.cv_ttol {
                break;
            }
            continue;
        }

        if zroot {
            /* No sign change in (tlo,tmid), but g = 0 at tmid; root tmid. */
            cv_mem.cv_thi = tmid;
            for i in 0..nrt {
                cv_mem.cv_ghi[i] = cv_mem.cv_grout[i];
            }
            break;
        }

        /* No sign change, no zero at tmid: sign change must be in
           (tmid,thi). Replace tlo with tmid. */
        cv_mem.cv_tlo = tmid;
        for i in 0..nrt {
            cv_mem.cv_glo[i] = cv_mem.cv_grout[i];
        }
        side = 2;
        /* Stop at root thi if converged; otherwise loop back. */
        if SUNRabs(cv_mem.cv_thi - cv_mem.cv_tlo) <= cv_mem.cv_ttol {
            break;
        }
    } /* End of root-search loop */

    /* Reset trout and grout, set iroots, and return RTFOUND. */
    cv_mem.cv_trout = cv_mem.cv_thi;
    for i in 0..nrt {
        cv_mem.cv_grout[i] = cv_mem.cv_ghi[i];
        cv_mem.cv_iroots[i] = 0;
        if !cv_mem.cv_gactive[i] {
            continue;
        }
        if SUNRabs(cv_mem.cv_ghi[i]) == ZERO
            && cv_mem.cv_rootdir[i] as f64 * cv_mem.cv_glo[i] <= ZERO
        {
            cv_mem.cv_iroots[i] = if cv_mem.cv_glo[i] > 0.0 { -1 } else { 1 };
        }
        if SUNRdifferentsign(cv_mem.cv_glo[i], cv_mem.cv_ghi[i])
            && cv_mem.cv_rootdir[i] as f64 * cv_mem.cv_glo[i] <= ZERO
        {
            cv_mem.cv_iroots[i] = if cv_mem.cv_glo[i] > 0.0 { -1 } else { 1 };
        }
    }
    RTFOUND
}

/*
 * =================================================================
 * Internal EWT function
 * =================================================================
 */

/*
 * cvEwtSet core: sets the error weight vector ewt from zn[0]
 * according to itol (CV_SS or CV_SV).  Returns 0 on success, -1 if
 * ewt would have a non-positive component.
 *
 * (In C this is cvEwtSet(ycur, weight, data) with data = cv_mem;
 * the weight vector is detached from CVodeMem by the caller here so
 * ycur may be borrowed from the same struct. The C version uses
 * cv_tempv as scratch and inverts it into weight; computing directly
 * in `weight` performs the identical arithmetic.)
 */
pub fn cvEwtSet(cv_mem: &CVodeMem, ycur: &NVector, weight: &mut NVector) -> i32 {
    match cv_mem.cv_itol {
        CV_SS => cvEwtSetSS(cv_mem, ycur, weight),
        CV_SV => cvEwtSetSV(cv_mem, ycur, weight),
        _ => 0,
    }
}

fn cvEwtSetSS(cv_mem: &CVodeMem, ycur: &NVector, weight: &mut NVector) -> i32 {
    N_VAbs(ycur, weight);
    weight.scale_inplace(cv_mem.cv_reltol);
    weight.add_const_inplace(cv_mem.cv_Sabstol);
    if cv_mem.cv_atolmin0 && N_VMin(weight) <= ZERO {
        return -1;
    }
    weight.invert_inplace();
    0
}

fn cvEwtSetSV(cv_mem: &CVodeMem, ycur: &NVector, weight: &mut NVector) -> i32 {
    N_VAbs(ycur, weight);
    weight.linear_sum_with(cv_mem.cv_reltol, ONE, &cv_mem.cv_Vabstol);
    if cv_mem.cv_atolmin0 && N_VMin(weight) <= ZERO {
        return -1;
    }
    weight.invert_inplace();
    0
}

/* Dispatch the efun (user-supplied or internal) writing into `w`.
   ycur is always zn[0], matching every C call site. */
pub(crate) fn cv_efun_dispatch(cv_mem: &mut CVodeMem, w: &mut NVector) -> i32 {
    if cv_mem.cv_user_efun {
        let efun = cv_mem.cv_efun.unwrap();
        /* e_data = user_data for a user efun */
        let CVodeMem { cv_zn, cv_user_data, .. } = cv_mem;
        efun(&cv_zn[0], w, cv_user_data)
    } else {
        /* e_data = cv_mem for the internal efun */
        cvEwtSet(cv_mem, &cv_mem.cv_zn[0], w)
    }
}

/* Apply the efun to the integrator's own ewt vector. */
pub(crate) fn cv_efun_apply_to_ewt(cv_mem: &mut CVodeMem) -> i32 {
    let mut w = std::mem::take(&mut cv_mem.cv_ewt);
    let flag = cv_efun_dispatch(cv_mem, &mut w);
    cv_mem.cv_ewt = w;
    flag
}
