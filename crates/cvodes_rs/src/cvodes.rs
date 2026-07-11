/* -----------------------------------------------------------------
 * Translated from src/cvodes/cvodes.c (CVODES 7.7.0) — PART 1
 * (cvodes.c lines 1-2912): creation, allocation and
 * re-initialization functions (CVodeCreate, CVodeInit, CVodeReInit,
 * the tolerance families for states / quadratures / sensitivities /
 * quadrature sensitivities, CVodeQuadInit/ReInit, CVodeSensInit,
 * CVodeSensInit1, CVodeSensReInit, CVodeQuadSensInit/ReInit,
 * CVodeSensToggleOff, CVodeRootInit).
 *
 * The donor translation crates/cvode_rs/src/cvode.rs is reused
 * verbatim wherever the CVODES C matches the CVODE C; the
 * sensitivity / quadrature deltas follow the CVODES source.
 * Documented adaptations (all donor precedents):
 *  - CVodeCreate returns Box<CVodeMem> and panics on illegal lmm
 *    (a NULL object cannot exist in safe Rust); the sunctx == NULL
 *    check has no counterpart (&SUNContext is never null).
 *  - cvode_mem == NULL checks vanish (&mut CVodeMem is never null);
 *    N_Vector/array NULL checks map to .is_empty().
 *  - cvCheckNvector and the ops-table guards (nvspace/nvmin) are
 *    dropped: the concrete serial NVector implements every op.
 *  - malloc-failure branches vanish (Vec/Box allocation aborts);
 *    the cvAllocVectors/cvQuadAllocVectors/cvSensAllocVectors/
 *    cvQuadSensAllocVectors helpers (defined at cvodes.c:4631+,
 *    ported in PART 2 below) are therefore called as infallible.
 *  - The fused-op scratch arrays cv_cvals/cv_Xvecs/cv_Zvecs are
 *    omitted (pinned decision in cvodes_impl.rs); the
 *    N_VScaleVectorArray(.., ONE, ..) calls become per-vector
 *    copies and the CV_VECTOROP_ERR branches vanish.
 *  - cv_fS = cvSensRhsInternalDQ / cv_fS1 = cvSensRhs1InternalDQ /
 *    cv_fQS = cvQuadSensRhsInternalDQ cannot be stored through the
 *    public fn-pointer types (the internal DQ routines take
 *    &mut CVodeMem); per pinned decision 5 in cvodes_impl.rs the
 *    fields stay None and the cv_fSDQ/cv_fQSDQ flags steer dispatch
 *    to the internal routines (ported with the sensitivity RHS
 *    wrappers in a later part).
 *  - error paths keep the C control flow but drop the explicit
 *    free calls (cvSensFreeVectors / SUNNonlinSolFree): RAII.
 * -----------------------------------------------------------------*/
use crate::cvodes_impl::*;
use crate::nvector_serial::*;
use crate::sundials_context::SUNContext;
use crate::sundials_errors::SUN_ERR_ARG_CORRUPT;
use crate::sundials_nonlinearsolver::SUN_NLS_CONV_RECVR;
use crate::sundials_types::*;
use crate::sunnonlinsol_newton::SUNNonlinSol_Newton;

/*=================================================================*/
/* CVODE Private Constants                                         */
/*=================================================================*/

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;
/* (TINY, PT1, POINT2, FOURTH, HALF, PT9, ONEPT5, TWO, THREE, FOUR,
   FIVE, TWELVE, HUNDRED are first used by the CVode driver and its
   helpers; PART 2 defines them below the end-of-part marker.) */

/*=================================================================*/
/* CVODE Routine-Specific Constants                                */
/*=================================================================*/

/* Control constants for lower-level rootfinding functions */
pub const RTFOUND: i32 = 1;
pub const CLOSERT: i32 = 3;

/* Control constants for sensitivity DQ */
pub const CENTERED1: i32 = 1;
pub const CENTERED2: i32 = 2;
pub const FORWARD1: i32 = 3;
pub const FORWARD2: i32 = 4;

/* Control constants for type of sensitivity RHS */
pub const CV_ONESENS: i32 = 1;
pub const CV_ALLSENS: i32 = 2;

/* Control constants for tolerances (CV_NN/CV_SS/CV_SV/CV_WF live in
   cvodes_impl.rs; cvodes.c additionally defines CV_EE) */
pub const CV_EE: i32 = 4;

/* Algorithmic constants (CVodeCreate) */
const CORTES: f64 = 0.1;
/* (FUZZ_FACTOR, HLB_FACTOR, HUB_FACTOR, H_BIAS, MAX_ITERS belong to
   CVodeGetDky/cvHin — PART 2.) */

/*
 * =================================================================
 * Exported Functions Implementation
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Creation, allocation and re-initialization functions
 * -----------------------------------------------------------------
 */

/*
 * CVodeCreate
 *
 * CVodeCreate creates an internal memory block for a problem to
 * be solved by CVODES. (C returns NULL on illegal lmm; the Rust
 * translation panics with the same diagnostic since a null object
 * cannot exist in safe Rust.)
 */
pub fn CVodeCreate(lmm: i32, sunctx: &SUNContext) -> Box<CVodeMem> {
    /* Test inputs */
    if lmm != CV_ADAMS && lmm != CV_BDF {
        cvProcessError(None, 0, line!(), "CVodeCreate", file!(), MSGCV_BAD_LMM);
        panic!("{}", MSGCV_BAD_LMM);
    }

    let maxord = if lmm == CV_ADAMS {
        ADAMS_Q_MAX as i32
    } else {
        BDF_Q_MAX as i32
    };

    /* The C code memsets cv_mem to zero and then assigns the defaults
       below; fields the C leaves at the memset value get the matching
       zero/false/empty default in this struct literal. */
    Box::new(CVodeMem {
        /* Copy input parameters into cv_mem */
        cv_sunctx: sunctx.clone(),

        /* Set uround */
        cv_uround: SUN_UNIT_ROUNDOFF,

        /* Set default values for integrator optional inputs */
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

        /* Set default values for quad. optional inputs */
        cv_quadr: SUNFALSE,
        cv_fQ: None,
        cv_errconQ: SUNFALSE,
        cv_itolQ: CV_NN,
        cv_reltolQ: ZERO,
        cv_SabstolQ: ZERO,
        cv_VabstolQ: NVector::default(),
        cv_atolQmin0: SUNTRUE,

        /* Set default values for sensi. optional inputs
           (C: cv_fS = cvSensRhsInternalDQ, cv_fS1 = cvSensRhs1InternalDQ;
           here None + cv_fSDQ = SUNTRUE selects the internal DQ routines,
           see module header) */
        cv_sensi: SUNFALSE,
        cv_Ns: 0,
        cv_ism: 0,
        cv_fS: None,
        cv_fS1: None,
        cv_fSDQ: SUNTRUE,
        cv_ifS: CV_ONESENS,
        cv_p: Vec::new(),
        cv_pbar: Vec::new(),
        cv_plist: Vec::new(),
        cv_DQtype: CV_CENTERED,
        cv_DQrhomax: ZERO,
        cv_errconS: SUNFALSE,
        cv_itolS: CV_NN,
        cv_reltolS: ZERO,
        cv_SabstolS: Vec::new(),
        cv_VabstolS: Vec::new(),
        cv_atolSmin0: Vec::new(),

        /* Set default values for quad. sensi. optional inputs */
        cv_quadr_sensi: SUNFALSE,
        cv_fQS: None,
        cv_fQSDQ: SUNTRUE,
        cv_errconQS: SUNFALSE,
        cv_itolQS: CV_NN,
        cv_reltolQS: ZERO,
        cv_SabstolQS: Vec::new(),
        cv_VabstolQS: Vec::new(),
        cv_atolQSmin0: Vec::new(),

        /* Nordsieck history array and length-N work vectors
           (allocated by CVodeInit) */
        cv_zn: Vec::new(),
        cv_ewt: NVector::default(),
        cv_y: NVector::default(),
        cv_acor: NVector::default(),
        cv_tempv: NVector::default(),
        cv_ftemp: NVector::default(),
        cv_vtemp1: NVector::default(),
        cv_vtemp2: NVector::default(),
        cv_vtemp3: NVector::default(),

        /* Quadrature related vectors (allocated by CVodeQuadInit) */
        cv_znQ: Vec::new(),
        cv_ewtQ: NVector::default(),
        cv_yQ: NVector::default(),
        cv_acorQ: NVector::default(),
        cv_tempvQ: NVector::default(),

        /* Sensitivity related vectors (allocated by CVodeSensInit*) */
        cv_znS: std::array::from_fn(|_| Vec::new()),
        cv_ewtS: Vec::new(),
        cv_yS: Vec::new(),
        cv_acorS: Vec::new(),
        cv_tempvS: Vec::new(),
        cv_ftempS: Vec::new(),
        cv_stgr1alloc: SUNFALSE,

        /* Quadrature sensitivity related vectors */
        cv_znQS: std::array::from_fn(|_| Vec::new()),
        cv_ewtQS: Vec::new(),
        cv_yQS: Vec::new(),
        cv_acorQS: Vec::new(),
        cv_tempvQS: Vec::new(),
        cv_ftempQ: NVector::default(),

        /* Tstop information */
        cv_tstopset: SUNFALSE,
        cv_tstopinterp: SUNFALSE,
        cv_tstop: ZERO,

        /* Step data */
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
        cv_crateS: ZERO,
        cv_delp: ZERO,
        cv_acnrm: ZERO,
        cv_acnrmcur: SUNFALSE,
        cv_acnrmQ: ZERO,
        cv_acnrmS: ZERO,
        cv_acnrmScur: SUNFALSE,
        cv_acnrmQS: ZERO,
        cv_nlscoef: CORTES,
        cv_ncfS1: Vec::new(),

        /* Limits */
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

        /* Counters */
        cv_nst: 0,
        cv_nfe: 0,
        cv_nfQe: 0,
        cv_nfSe: 0,
        cv_nfeS: 0,
        cv_nfQSe: 0,
        cv_nfQeS: 0,
        cv_ncfn: 0,
        cv_ncfnS: 0,
        cv_ncfnS1: Vec::new(),
        cv_nni: 0,
        cv_nniS: 0,
        cv_nniS1: Vec::new(),
        cv_nnf: 0,
        cv_nnfS: 0,
        cv_nnfS1: Vec::new(),
        cv_netf: 0,
        cv_netfQ: 0,
        cv_netfS: 0,
        cv_netfQS: 0,
        cv_nsetups: 0,
        cv_nsetupsS: 0,
        cv_nhnil: 0,

        /* Step size ratios */
        cv_etaqm1: ZERO,
        cv_etaq: ZERO,
        cv_etaqp1: ZERO,

        /* Initialize lrw and liw */
        cv_lrw1: 0,
        cv_liw1: 0,
        cv_lrw1Q: 0,
        cv_liw1Q: 0,
        cv_lrw: (65 + 2 * L_MAX + NUM_TESTS) as i64,
        cv_liw: 52,

        /* Initialize nonlinear solver variables */
        NLS: None,
        cv_nls_curiter: 0,
        ownNLS: SUNFALSE,
        NLSsim: None,
        ownNLSsim: SUNFALSE,
        simMallocDone: SUNFALSE,
        NLSstg: None,
        ownNLSstg: SUNFALSE,
        stgMallocDone: SUNFALSE,
        NLSstg1: None,
        ownNLSstg1: SUNFALSE,
        sens_solve: SUNFALSE,
        sens_solve_idx: -1,
        nnip: 0,
        nls_f: None,
        convfail: CV_NO_FAILURES,

        /* Linear solver data (attached later) */
        cv_lmem: LsModule::None,
        cv_msbp: MSBP_DEFAULT,
        cv_dgmax_lsetup: DGMAX_LSETUP_DEFAULT,
        cv_forceSetup: SUNFALSE,

        /* Saved values */
        cv_qu: 0,
        cv_nstlp: 0,
        cv_h0u: ZERO,
        cv_hu: ZERO,
        cv_saved_tq5: ZERO,
        cv_jcur: SUNFALSE,
        cv_convfail: 0,
        cv_tolsf: ZERO,
        /* Set the saved values for qmax_alloc (cv_qmax_allocQS is only
           set by cvQuadSensAllocVectors; memset value 0 here) */
        cv_qmax_alloc: maxord,
        cv_qmax_allocQ: maxord,
        cv_qmax_allocS: maxord,
        cv_qmax_allocQS: 0,
        cv_indx_acor: 0,

        /* No mallocs have been done yet */
        cv_VabstolMallocDone: SUNFALSE,
        cv_MallocDone: SUNFALSE,
        cv_constraintsMallocDone: SUNFALSE,
        cv_VabstolQMallocDone: SUNFALSE,
        cv_QuadMallocDone: SUNFALSE,
        cv_VabstolSMallocDone: SUNFALSE,
        cv_SabstolSMallocDone: SUNFALSE,
        cv_SensMallocDone: SUNFALSE,
        cv_VabstolQSMallocDone: SUNFALSE,
        cv_SabstolQSMallocDone: SUNFALSE,
        cv_QuadSensMallocDone: SUNFALSE,

        /* Stability Limit Detection data */
        cv_sldeton: SUNFALSE,
        cv_ssdat: [[ZERO; 4]; 6],
        cv_nscon: 0,
        cv_nor: 0,

        /* Initialize root finding variables */
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
        cv_irfnd: 0,
        cv_nge: 0,
        cv_gactive: Vec::new(),
        cv_mxgnull: 1,

        /* Initialize inequality constraint variables */
        cv_constraints: NVector::default(),
        constraint_corrections: 0,
        constraint_fails: 0,
        max_constraint_fails: MAX_CONSTRAINT_FAILS,

        /* Initialize projection variables */
        proj_mem: None,
        proj_enabled: SUNFALSE,
        proj_applied: SUNFALSE,
        proj_p: [ZERO; L_MAX],

        /* Initialize resize variables */
        first_step_after_resize: SUNFALSE,

        /* Set default for ASA */
        cv_adj: SUNFALSE,
        cv_adj_mem: None,
        cv_adjMallocDone: SUNFALSE,
    })
}

/*-----------------------------------------------------------------*/

/*
 * CVodeInit
 *
 * CVodeInit allocates and initializes memory for a problem. All
 * problem inputs are checked for errors. If any error occurs during
 * initialization, an error flag is returned. Otherwise, it returns
 * CV_SUCCESS.
 */
pub fn CVodeInit(cv_mem: &mut CVodeMem, f: CVRhsFn, t0: f64, y0: &NVector) -> i32 {
    /* Check for legal input parameters */
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

    /* Input checks complete at this point and history array allocated */

    /* Copy the input parameters into CVODES state */
    cv_mem.cv_f = Some(f);
    cv_mem.cv_tn = t0;

    /* Initialize zn[0] in the history array */
    cv_mem.cv_zn[0].data.copy_from_slice(&y0.data);

    /* create a Newton nonlinear solver object by default */
    let nls = SUNNonlinSol_Newton(y0, &cv_mem.cv_sunctx);

    /* attach the nonlinear solver to the CVODE memory */
    let retval = crate::cvodes_nls::CVodeSetNonlinearSolver(cv_mem, nls);

    /* check that the nonlinear solver was successfully attached */
    if retval != CV_SUCCESS {
        cvProcessError(Some(cv_mem), retval, line!(), "CVodeInit", file!(),
                       "Setting the nonlinear solver failed");
        return CV_MEM_FAIL;
    }

    /* set ownership flag */
    cv_mem.ownNLS = SUNTRUE;

    /* All error checking is complete at this point */

    /* Set step parameters */
    cv_mem.cv_q = 1;
    cv_mem.cv_L = 2;
    cv_mem.cv_qwait = cv_mem.cv_L;
    cv_mem.cv_etamax = cv_mem.cv_eta_max_fs;

    cv_mem.cv_qu = 0;
    cv_mem.cv_hu = ZERO;
    cv_mem.cv_tolsf = ONE;

    /* Set the linear solver addresses to NULL.
       (We check != NULL later, in CVode) */
    cv_mem.cv_lmem = LsModule::None;

    /* Set forceSetup to SUNFALSE */
    cv_mem.cv_forceSetup = SUNFALSE;

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

    /* Initialize Stablilty Limit Detection data */
    /* NOTE: We do this even if stab lim det was not
       turned on yet. This way, the user can turn it
       on at any time */
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

/*-----------------------------------------------------------------*/

/*
 * CVodeReInit
 *
 * CVodeReInit re-initializes CVODES's memory for a problem, assuming
 * it has already been allocated in a prior CVodeInit call.
 */
pub fn CVodeReInit(cv_mem: &mut CVodeMem, t0: f64, y0: &NVector) -> i32 {
    /* Check if cvode_mem was allocated */
    if !cv_mem.cv_MallocDone {
        cvProcessError(Some(cv_mem), CV_NO_MALLOC, line!(), "CVodeReInit", file!(), MSGCV_NO_MALLOC);
        return CV_NO_MALLOC;
    }

    /* Check for legal input parameters */
    if y0.is_empty() {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeReInit", file!(), MSGCV_NULL_Y0);
        return CV_ILL_INPUT;
    }

    /* Copy the input parameters into CVODES state */
    cv_mem.cv_tn = t0;

    /* Set step parameters */
    cv_mem.cv_q = 1;
    cv_mem.cv_L = 2;
    cv_mem.cv_qwait = cv_mem.cv_L;
    cv_mem.cv_etamax = cv_mem.cv_eta_max_fs;

    cv_mem.cv_qu = 0;
    cv_mem.cv_hu = ZERO;
    cv_mem.cv_tolsf = ONE;

    /* Set forceSetup to SUNFALSE */
    cv_mem.cv_forceSetup = SUNFALSE;

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

    /* lreinit: reinitialize the attached linear solver interface
       (C: if (cv_mem->cv_lreinit) { cv_mem->cv_lreinit(cv_mem); }) */
    cv_lreinit_dispatch(cv_mem);

    /* Initialize other integrator optional outputs */
    cv_mem.cv_h0u = ZERO;
    cv_mem.cv_next_h = ZERO;
    cv_mem.cv_next_q = 0;

    /* Initialize Stablilty Limit Detection data */
    cv_mem.cv_nor = 0;
    for i in 1..=5usize {
        for k in 1..=3usize {
            cv_mem.cv_ssdat[i - 1][k - 1] = ZERO;
        }
    }

    /* Problem has been successfully re-initialized */
    CV_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Linear solver module dispatch (donor pattern: in C these are the
 * cv_linit/cv_lreinit/... function pointers; the module is taken
 * out of CVodeMem during the call so its routine can borrow the
 * integrator memory mutably. The remaining dispatch helpers land
 * with the CVode driver in PART 2.)
 * -----------------------------------------------------------------
 */

fn cv_lreinit_dispatch(cv_mem: &mut CVodeMem) {
    /* cv_lreinit is only set by the CVLS interface; it resets counters */
    let mut lmem = std::mem::take(&mut cv_mem.cv_lmem);
    if let LsModule::Ls(ls) = &mut lmem {
        crate::cvodes_ls::cvLsReInitialize(cv_mem, ls);
    }
    cv_mem.cv_lmem = lmem;
}

/*-----------------------------------------------------------------*/

/*
 * CVodeSStolerances
 * CVodeSVtolerances
 * CVodeWFtolerances
 *
 * These functions specify the integration tolerances. One of them
 * MUST be called before the first call to CVode.
 *
 * CVodeSStolerances specifies scalar relative and absolute tolerances.
 * CVodeSVtolerances specifies scalar relative tolerance and a vector
 *   absolute tolerance (a potentially different absolute tolerance
 *   for each vector component).
 * CVodeWFtolerances specifies a user-provides function (of type CVEwtFn)
 *   which will be called to set the error weight vector.
 */

pub fn CVodeSStolerances(cv_mem: &mut CVodeMem, reltol: f64, abstol: f64) -> i32 {
    if !cv_mem.cv_MallocDone {
        cvProcessError(Some(cv_mem), CV_NO_MALLOC, line!(), "CVodeSStolerances", file!(), MSGCV_NO_MALLOC);
        return CV_NO_MALLOC;
    }

    /* Check inputs */
    if reltol < ZERO {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSStolerances", file!(), MSGCV_BAD_RELTOL);
        return CV_ILL_INPUT;
    }
    if abstol < ZERO {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSStolerances", file!(), MSGCV_BAD_ABSTOL);
        return CV_ILL_INPUT;
    }

    /* Copy tolerances into memory */
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

    /* Check inputs */
    if reltol < ZERO {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSVtolerances", file!(), MSGCV_BAD_RELTOL);
        return CV_ILL_INPUT;
    }
    let atolmin = N_VMin(abstol);
    if atolmin < ZERO {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSVtolerances", file!(), MSGCV_BAD_ABSTOL);
        return CV_ILL_INPUT;
    }

    /* Copy tolerances into memory */
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

/*-----------------------------------------------------------------*/

/*
 * CVodeQuadInit
 *
 * CVodeQuadInit allocates and initializes quadrature related
 * memory for a problem. All problem specification inputs are
 * checked for errors.
 */
pub fn CVodeQuadInit(cv_mem: &mut CVodeMem, fQ: CVQuadRhsFn, yQ0: &NVector) -> i32 {
    /* Set space requirements for one N_Vector */
    let (lrw1Q, liw1Q) = N_VSpace(yQ0);
    cv_mem.cv_lrw1Q = lrw1Q;
    cv_mem.cv_liw1Q = liw1Q;

    /* Allocate the vectors (using yQ0 as a template) */
    cvQuadAllocVectors(cv_mem, yQ0);

    /* Initialize znQ[0] in the history array */
    cv_mem.cv_znQ[0].data.copy_from_slice(&yQ0.data);

    /* Copy the input parameters into CVODES state */
    cv_mem.cv_fQ = Some(fQ);

    /* Initialize counters */
    cv_mem.cv_nfQe = 0;
    cv_mem.cv_netfQ = 0;

    /* Quadrature integration turned ON */
    cv_mem.cv_quadr = SUNTRUE;
    cv_mem.cv_QuadMallocDone = SUNTRUE;

    /* Quadrature initialization was successful */
    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

/*
 * CVodeQuadReInit
 *
 * CVodeQuadReInit re-initializes CVODES's quadrature related memory
 * for a problem, assuming it has already been allocated in prior
 * calls to CVodeInit and CVodeQuadInit.
 */
pub fn CVodeQuadReInit(cv_mem: &mut CVodeMem, yQ0: &NVector) -> i32 {
    /* Check if quadrature was initialized? */
    if !cv_mem.cv_QuadMallocDone {
        cvProcessError(Some(cv_mem), CV_NO_QUAD, line!(), "CVodeQuadReInit", file!(), MSGCV_NO_QUAD);
        return CV_NO_QUAD;
    }

    /* Initialize znQ[0] in the history array */
    cv_mem.cv_znQ[0].data.copy_from_slice(&yQ0.data);

    /* Initialize counters */
    cv_mem.cv_nfQe = 0;
    cv_mem.cv_netfQ = 0;

    /* Quadrature integration turned ON */
    cv_mem.cv_quadr = SUNTRUE;

    /* Quadrature re-initialization was successful */
    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

/*
 * CVodeQuadSStolerances
 * CVodeQuadSVtolerances
 *
 * These functions specify the integration tolerances for quadrature
 * variables. One of them MUST be called before the first call to
 * CVode IF error control on the quadrature variables is enabled
 * (see CVodeSetQuadErrCon).
 */

pub fn CVodeQuadSStolerances(cv_mem: &mut CVodeMem, reltolQ: f64, abstolQ: f64) -> i32 {
    /* Check if quadrature was initialized? */
    if !cv_mem.cv_QuadMallocDone {
        cvProcessError(Some(cv_mem), CV_NO_QUAD, line!(), "CVodeQuadSStolerances", file!(), MSGCV_NO_QUAD);
        return CV_NO_QUAD;
    }

    /* Test user-supplied tolerances */
    if reltolQ < ZERO {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeQuadSStolerances", file!(), MSGCV_BAD_RELTOLQ);
        return CV_ILL_INPUT;
    }
    if abstolQ < ZERO {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeQuadSStolerances", file!(), MSGCV_BAD_ABSTOLQ);
        return CV_ILL_INPUT;
    }

    /* Copy tolerances into memory */
    cv_mem.cv_itolQ = CV_SS;

    cv_mem.cv_reltolQ = reltolQ;
    cv_mem.cv_SabstolQ = abstolQ;
    cv_mem.cv_atolQmin0 = abstolQ == ZERO;

    CV_SUCCESS
}

pub fn CVodeQuadSVtolerances(cv_mem: &mut CVodeMem, reltolQ: f64, abstolQ: &NVector) -> i32 {
    /* Check if quadrature was initialized? */
    if !cv_mem.cv_QuadMallocDone {
        cvProcessError(Some(cv_mem), CV_NO_QUAD, line!(), "CVodeQuadSVtolerances", file!(), MSGCV_NO_QUAD);
        return CV_NO_QUAD;
    }

    /* Test user-supplied tolerances */
    if reltolQ < ZERO {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeQuadSVtolerances", file!(), MSGCV_BAD_RELTOLQ);
        return CV_ILL_INPUT;
    }

    if abstolQ.is_empty() {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeQuadSVtolerances", file!(), MSGCV_NULL_ABSTOLQ);
        return CV_ILL_INPUT;
    }

    let atolmin = N_VMin(abstolQ);
    if atolmin < ZERO {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeQuadSVtolerances", file!(), MSGCV_BAD_ABSTOLQ);
        return CV_ILL_INPUT;
    }

    /* Copy tolerances into memory */
    cv_mem.cv_itolQ = CV_SV;

    cv_mem.cv_reltolQ = reltolQ;

    if !cv_mem.cv_VabstolQMallocDone {
        cv_mem.cv_VabstolQ = N_VClone(&cv_mem.cv_tempvQ);
        cv_mem.cv_lrw += cv_mem.cv_lrw1Q;
        cv_mem.cv_liw += cv_mem.cv_liw1Q;
        cv_mem.cv_VabstolQMallocDone = SUNTRUE;
    }

    cv_mem.cv_VabstolQ.data.copy_from_slice(&abstolQ.data);
    cv_mem.cv_atolQmin0 = atolmin == ZERO;

    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

/*
 * CVodeSensInit
 *
 * CVodeSensInit allocates and initializes sensitivity related
 * memory for a problem (using a sensitivity RHS function of type
 * CVSensRhsFn). All problem specification inputs are checked for
 * errors.
 */
pub fn CVodeSensInit(
    cv_mem: &mut CVodeMem,
    Ns: i32,
    ism: i32,
    fS: Option<CVSensRhsFn>,
    yS0: &[NVector],
) -> i32 {
    /* Check if CVodeSensInit or CVodeSensInit1 was already called */
    if cv_mem.cv_SensMallocDone {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSensInit", file!(), MSGCV_SENSINIT_2);
        return CV_ILL_INPUT;
    }

    /* Check if Ns is legal */
    if Ns <= 0 {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSensInit", file!(), MSGCV_BAD_NS);
        return CV_ILL_INPUT;
    }
    cv_mem.cv_Ns = Ns;

    /* Check if ism is compatible */
    if ism == CV_STAGGERED1 {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSensInit", file!(), MSGCV_BAD_ISM_IFS);
        return CV_ILL_INPUT;
    }

    /* Check if ism is legal */
    if ism != CV_SIMULTANEOUS && ism != CV_STAGGERED {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSensInit", file!(), MSGCV_BAD_ISM);
        return CV_ILL_INPUT;
    }
    cv_mem.cv_ism = ism;

    /* Check if yS0 is non-null */
    if yS0.is_empty() {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSensInit", file!(), MSGCV_NULL_YS0);
        return CV_ILL_INPUT;
    }

    /* Store sensitivity RHS-related data */
    cv_mem.cv_ifS = CV_ALLSENS;
    cv_mem.cv_fS1 = None;

    match fS {
        None => {
            /* C: cv_fS = cvSensRhsInternalDQ, cv_fS_data = cvode_mem;
               here cv_fS stays None and cv_fSDQ steers the dispatch to
               the internal DQ routine (module header) */
            cv_mem.cv_fSDQ = SUNTRUE;
            cv_mem.cv_fS = None;
        }
        Some(f) => {
            cv_mem.cv_fSDQ = SUNFALSE;
            cv_mem.cv_fS = Some(f);
        }
    }

    /* No memory allocation for STAGGERED1 */
    cv_mem.cv_stgr1alloc = SUNFALSE;

    /* Allocate the vectors (using yS0[0] as a template) */
    cvSensAllocVectors(cv_mem, &yS0[0]);

    /* (Fused-op work arrays cv_cvals/cv_Xvecs/cv_Zvecs omitted) */

    /*----------------------------------------------
      All error checking is complete at this point
      -----------------------------------------------*/

    /* Initialize znS[0] in the history array */
    for is in 0..Ns as usize {
        cv_mem.cv_znS[0][is].data.copy_from_slice(&yS0[is].data);
    }

    /* Initialize all sensitivity related counters */
    cv_mem.cv_nfSe = 0;
    cv_mem.cv_nfeS = 0;
    cv_mem.cv_ncfnS = 0;
    cv_mem.cv_netfS = 0;
    cv_mem.cv_nniS = 0;
    cv_mem.cv_nnfS = 0;
    cv_mem.cv_nsetupsS = 0;

    /* Set default values for plist and pbar */
    for is in 0..Ns as usize {
        cv_mem.cv_plist[is] = is as i32;
        cv_mem.cv_pbar[is] = ONE;
    }

    /* Sensitivities will be computed */
    cv_mem.cv_sensi = SUNTRUE;
    cv_mem.cv_SensMallocDone = SUNTRUE;

    /* create a Newton nonlinear solver object by default */
    let NLS = if ism == CV_SIMULTANEOUS {
        crate::sunnonlinsol_newton::SUNNonlinSol_NewtonSens(Ns + 1, &cv_mem.cv_acor,
                                                            &cv_mem.cv_sunctx)
    } else {
        crate::sunnonlinsol_newton::SUNNonlinSol_NewtonSens(Ns, &cv_mem.cv_acor,
                                                            &cv_mem.cv_sunctx)
    };

    /* attach the nonlinear solver to the CVODE memory */
    let retval = if ism == CV_SIMULTANEOUS {
        crate::cvodes_nls_sim::CVodeSetNonlinearSolverSensSim(cv_mem, NLS)
    } else {
        crate::cvodes_nls_stg::CVodeSetNonlinearSolverSensStg(cv_mem, NLS)
    };

    /* check that the nonlinear solver was successfully attached */
    if retval != CV_SUCCESS {
        cvProcessError(Some(cv_mem), retval, line!(), "CVodeSensInit", file!(),
                       "Setting the nonlinear solver failed");
        return CV_MEM_FAIL;
    }

    /* set ownership flag */
    if ism == CV_SIMULTANEOUS {
        cv_mem.ownNLSsim = SUNTRUE;
    } else {
        cv_mem.ownNLSstg = SUNTRUE;
    }

    /* Sensitivity initialization was successful */
    CV_SUCCESS
}

/*
 * CVodeSensInit1
 *
 * CVodeSensInit1 allocates and initializes sensitivity related
 * memory for a problem (using a sensitivity RHS function of type
 * CVSensRhs1Fn). All problem specification inputs are checked for
 * errors.
 */
pub fn CVodeSensInit1(
    cv_mem: &mut CVodeMem,
    Ns: i32,
    ism: i32,
    fS1: Option<CVSensRhs1Fn>,
    yS0: &[NVector],
) -> i32 {
    /* Check if CVodeSensInit or CVodeSensInit1 was already called */
    if cv_mem.cv_SensMallocDone {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSensInit1", file!(), MSGCV_SENSINIT_2);
        return CV_ILL_INPUT;
    }

    /* Check if Ns is legal */
    if Ns <= 0 {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSensInit1", file!(), MSGCV_BAD_NS);
        return CV_ILL_INPUT;
    }
    cv_mem.cv_Ns = Ns;

    /* Check if ism is legal */
    if ism != CV_SIMULTANEOUS && ism != CV_STAGGERED && ism != CV_STAGGERED1 {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSensInit1", file!(), MSGCV_BAD_ISM);
        return CV_ILL_INPUT;
    }
    cv_mem.cv_ism = ism;

    /* Check if yS0 is non-null */
    if yS0.is_empty() {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSensInit1", file!(), MSGCV_NULL_YS0);
        return CV_ILL_INPUT;
    }

    /* Store sensitivity RHS-related data */
    cv_mem.cv_ifS = CV_ONESENS;
    cv_mem.cv_fS = None;

    match fS1 {
        None => {
            /* C: cv_fS1 = cvSensRhs1InternalDQ, cv_fS_data = cvode_mem;
               here cv_fS1 stays None and cv_fSDQ steers the dispatch to
               the internal DQ routine (module header) */
            cv_mem.cv_fSDQ = SUNTRUE;
            cv_mem.cv_fS1 = None;
        }
        Some(f) => {
            cv_mem.cv_fSDQ = SUNFALSE;
            cv_mem.cv_fS1 = Some(f);
        }
    }

    /* Allocate ncfS1, ncfnS1, and nniS1 if needed */
    if ism == CV_STAGGERED1 {
        cv_mem.cv_stgr1alloc = SUNTRUE;
        cv_mem.cv_ncfS1 = vec![0; Ns as usize];
        cv_mem.cv_ncfnS1 = vec![0; Ns as usize];
        cv_mem.cv_nniS1 = vec![0; Ns as usize];
        cv_mem.cv_nnfS1 = vec![0; Ns as usize];
    } else {
        cv_mem.cv_stgr1alloc = SUNFALSE;
    }

    /* Allocate the vectors (using yS0[0] as a template) */
    cvSensAllocVectors(cv_mem, &yS0[0]);

    /* (Fused-op work arrays cv_cvals/cv_Xvecs/cv_Zvecs omitted) */

    /*----------------------------------------------
      All error checking is complete at this point
      -----------------------------------------------*/

    /* Initialize znS[0] in the history array */
    for is in 0..Ns as usize {
        cv_mem.cv_znS[0][is].data.copy_from_slice(&yS0[is].data);
    }

    /* Initialize all sensitivity related counters */
    cv_mem.cv_nfSe = 0;
    cv_mem.cv_nfeS = 0;
    cv_mem.cv_ncfnS = 0;
    cv_mem.cv_netfS = 0;
    cv_mem.cv_nniS = 0;
    cv_mem.cv_nnfS = 0;
    cv_mem.cv_nsetupsS = 0;
    if ism == CV_STAGGERED1 {
        for is in 0..Ns as usize {
            cv_mem.cv_ncfnS1[is] = 0;
            cv_mem.cv_nniS1[is] = 0;
            cv_mem.cv_nnfS1[is] = 0;
        }
    }

    /* Set default values for plist and pbar */
    for is in 0..Ns as usize {
        cv_mem.cv_plist[is] = is as i32;
        cv_mem.cv_pbar[is] = ONE;
    }

    /* Sensitivities will be computed */
    cv_mem.cv_sensi = SUNTRUE;
    cv_mem.cv_SensMallocDone = SUNTRUE;

    /* create a Newton nonlinear solver object by default */
    let NLS = if ism == CV_SIMULTANEOUS {
        crate::sunnonlinsol_newton::SUNNonlinSol_NewtonSens(Ns + 1, &cv_mem.cv_acor,
                                                            &cv_mem.cv_sunctx)
    } else if ism == CV_STAGGERED {
        crate::sunnonlinsol_newton::SUNNonlinSol_NewtonSens(Ns, &cv_mem.cv_acor,
                                                            &cv_mem.cv_sunctx)
    } else {
        SUNNonlinSol_Newton(&cv_mem.cv_acor, &cv_mem.cv_sunctx)
    };

    /* attach the nonlinear solver to the CVODE memory */
    let retval = if ism == CV_SIMULTANEOUS {
        crate::cvodes_nls_sim::CVodeSetNonlinearSolverSensSim(cv_mem, NLS)
    } else if ism == CV_STAGGERED {
        crate::cvodes_nls_stg::CVodeSetNonlinearSolverSensStg(cv_mem, NLS)
    } else {
        crate::cvodes_nls_stg1::CVodeSetNonlinearSolverSensStg1(cv_mem, NLS)
    };

    /* check that the nonlinear solver was successfully attached */
    if retval != CV_SUCCESS {
        cvProcessError(Some(cv_mem), retval, line!(), "CVodeSensInit1", file!(),
                       "Setting the nonlinear solver failed");
        return CV_MEM_FAIL;
    }

    /* set ownership flag */
    if ism == CV_SIMULTANEOUS {
        cv_mem.ownNLSsim = SUNTRUE;
    } else if ism == CV_STAGGERED {
        cv_mem.ownNLSstg = SUNTRUE;
    } else {
        cv_mem.ownNLSstg1 = SUNTRUE;
    }

    /* Sensitivity initialization was successful */
    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

/*
 * CVodeSensReInit
 *
 * CVodeSensReInit re-initializes CVODES's sensitivity related memory
 * for a problem, assuming it has already been allocated in prior
 * calls to CVodeInit and CVodeSensInit/CVodeSensInit1.
 * The number of sensitivities Ns is assumed to be unchanged since
 * the previous call to CVodeSensInit.
 */
pub fn CVodeSensReInit(cv_mem: &mut CVodeMem, ism: i32, yS0: &[NVector]) -> i32 {
    /* Was sensitivity initialized? */
    if !cv_mem.cv_SensMallocDone {
        cvProcessError(Some(cv_mem), CV_NO_SENS, line!(), "CVodeSensReInit", file!(), MSGCV_NO_SENSI);
        return CV_NO_SENS;
    }

    /* Check if ism is compatible */
    if cv_mem.cv_ifS == CV_ALLSENS && ism == CV_STAGGERED1 {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSensReInit", file!(), MSGCV_BAD_ISM_IFS);
        return CV_ILL_INPUT;
    }

    /* Check if ism is legal */
    if ism != CV_SIMULTANEOUS && ism != CV_STAGGERED && ism != CV_STAGGERED1 {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSensReInit", file!(), MSGCV_BAD_ISM);
        return CV_ILL_INPUT;
    }
    cv_mem.cv_ism = ism;

    /* Check if yS0 is non-null */
    if yS0.is_empty() {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSensReInit", file!(), MSGCV_NULL_YS0);
        return CV_ILL_INPUT;
    }

    /* Allocate ncfS1, ncfnS1, and nniS1 if needed */
    if ism == CV_STAGGERED1 && !cv_mem.cv_stgr1alloc {
        cv_mem.cv_stgr1alloc = SUNTRUE;
        cv_mem.cv_ncfS1 = vec![0; cv_mem.cv_Ns as usize];
        cv_mem.cv_ncfnS1 = vec![0; cv_mem.cv_Ns as usize];
        cv_mem.cv_nniS1 = vec![0; cv_mem.cv_Ns as usize];
        cv_mem.cv_nnfS1 = vec![0; cv_mem.cv_Ns as usize];
    }

    /*----------------------------------------------
      All error checking is complete at this point
      -----------------------------------------------*/

    /* Initialize znS[0] in the history array */
    for is in 0..cv_mem.cv_Ns as usize {
        cv_mem.cv_znS[0][is].data.copy_from_slice(&yS0[is].data);
    }

    /* Initialize all sensitivity related counters */
    cv_mem.cv_nfSe = 0;
    cv_mem.cv_nfeS = 0;
    cv_mem.cv_ncfnS = 0;
    cv_mem.cv_netfS = 0;
    cv_mem.cv_nniS = 0;
    cv_mem.cv_nnfS = 0;
    cv_mem.cv_nsetupsS = 0;
    if ism == CV_STAGGERED1 {
        for is in 0..cv_mem.cv_Ns as usize {
            cv_mem.cv_ncfnS1[is] = 0;
            cv_mem.cv_nniS1[is] = 0;
            cv_mem.cv_nnfS1[is] = 0;
        }
    }

    /* Problem has been successfully re-initialized */
    cv_mem.cv_sensi = SUNTRUE;

    /* Check if the NLS exists, create the default NLS if needed */
    if (ism == CV_SIMULTANEOUS && cv_mem.NLSsim.is_none())
        || (ism == CV_STAGGERED && cv_mem.NLSstg.is_none())
        || (ism == CV_STAGGERED1 && cv_mem.NLSstg1.is_none())
    {
        /* create a Newton nonlinear solver object by default */
        let NLS = if ism == CV_SIMULTANEOUS {
            crate::sunnonlinsol_newton::SUNNonlinSol_NewtonSens(cv_mem.cv_Ns + 1,
                                                                &cv_mem.cv_acor,
                                                                &cv_mem.cv_sunctx)
        } else if ism == CV_STAGGERED {
            crate::sunnonlinsol_newton::SUNNonlinSol_NewtonSens(cv_mem.cv_Ns,
                                                                &cv_mem.cv_acor,
                                                                &cv_mem.cv_sunctx)
        } else {
            SUNNonlinSol_Newton(&cv_mem.cv_acor, &cv_mem.cv_sunctx)
        };

        /* attach the nonlinear solver to the CVODES memory */
        let retval = if ism == CV_SIMULTANEOUS {
            crate::cvodes_nls_sim::CVodeSetNonlinearSolverSensSim(cv_mem, NLS)
        } else if ism == CV_STAGGERED {
            crate::cvodes_nls_stg::CVodeSetNonlinearSolverSensStg(cv_mem, NLS)
        } else {
            crate::cvodes_nls_stg1::CVodeSetNonlinearSolverSensStg1(cv_mem, NLS)
        };

        /* check that the nonlinear solver was successfully attached */
        if retval != CV_SUCCESS {
            cvProcessError(Some(cv_mem), retval, line!(), "CVodeSensReInit", file!(),
                           "Setting the nonlinear solver failed");
            return CV_MEM_FAIL;
        }

        /* set ownership flag */
        if ism == CV_SIMULTANEOUS {
            cv_mem.ownNLSsim = SUNTRUE;
        } else if ism == CV_STAGGERED {
            cv_mem.ownNLSstg = SUNTRUE;
        } else {
            cv_mem.ownNLSstg1 = SUNTRUE;
        }

        /* initialize the NLS object, this assumes that the linear solver has
           already been initialized in CVodeInit */
        let retval = if ism == CV_SIMULTANEOUS {
            crate::cvodes_nls_sim::cvNlsInitSensSim(cv_mem)
        } else if ism == CV_STAGGERED {
            crate::cvodes_nls_stg::cvNlsInitSensStg(cv_mem)
        } else {
            crate::cvodes_nls_stg1::cvNlsInitSensStg1(cv_mem)
        };

        if retval != CV_SUCCESS {
            cvProcessError(Some(cv_mem), CV_NLS_INIT_FAIL, line!(), "CVodeSensReInit", file!(),
                           MSGCV_NLS_INIT_FAIL);
            return CV_NLS_INIT_FAIL;
        }
    }

    /* Sensitivity re-initialization was successful */
    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

/*
 * CVodeSensSStolerances
 * CVodeSensSVtolerances
 * CVodeSensEEtolerances
 *
 * These functions specify the integration tolerances for sensitivity
 * variables. One of them MUST be called before the first call to CVode.
 *
 * CVodeSensSStolerances specifies scalar relative and absolute tolerances.
 * CVodeSensSVtolerances specifies scalar relative tolerance and a vector
 *   absolute tolerance for each sensitivity vector (a potentially different
 *   absolute tolerance for each vector component).
 * CVodeSensEEtolerances specifies that tolerances for sensitivity variables
 *   should be estimated from those provided for the state variables.
 */

pub fn CVodeSensSStolerances(cv_mem: &mut CVodeMem, reltolS: f64, abstolS: &[f64]) -> i32 {
    /* Was sensitivity initialized? */
    if !cv_mem.cv_SensMallocDone {
        cvProcessError(Some(cv_mem), CV_NO_SENS, line!(), "CVodeSensSStolerances", file!(), MSGCV_NO_SENSI);
        return CV_NO_SENS;
    }

    /* Test user-supplied tolerances */
    if reltolS < ZERO {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSensSStolerances", file!(), MSGCV_BAD_RELTOLS);
        return CV_ILL_INPUT;
    }

    if abstolS.is_empty() {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSensSStolerances", file!(), MSGCV_NULL_ABSTOLS);
        return CV_ILL_INPUT;
    }

    for is in 0..cv_mem.cv_Ns as usize {
        if abstolS[is] < ZERO {
            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSensSStolerances", file!(), MSGCV_BAD_ABSTOLS);
            return CV_ILL_INPUT;
        }
    }

    /* Copy tolerances into memory */
    cv_mem.cv_itolS = CV_SS;

    cv_mem.cv_reltolS = reltolS;

    if !cv_mem.cv_SabstolSMallocDone {
        cv_mem.cv_SabstolS = vec![ZERO; cv_mem.cv_Ns as usize];
        cv_mem.cv_atolSmin0 = vec![SUNFALSE; cv_mem.cv_Ns as usize];
        cv_mem.cv_lrw += cv_mem.cv_Ns as i64;
        cv_mem.cv_SabstolSMallocDone = SUNTRUE;
    }

    for is in 0..cv_mem.cv_Ns as usize {
        cv_mem.cv_SabstolS[is] = abstolS[is];
        cv_mem.cv_atolSmin0[is] = abstolS[is] == ZERO;
    }

    CV_SUCCESS
}

pub fn CVodeSensSVtolerances(cv_mem: &mut CVodeMem, reltolS: f64, abstolS: &[NVector]) -> i32 {
    /* Was sensitivity initialized? */
    if !cv_mem.cv_SensMallocDone {
        cvProcessError(Some(cv_mem), CV_NO_SENS, line!(), "CVodeSensSVtolerances", file!(), MSGCV_NO_SENSI);
        return CV_NO_SENS;
    }

    /* Test user-supplied tolerances */
    if reltolS < ZERO {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSensSVtolerances", file!(), MSGCV_BAD_RELTOLS);
        return CV_ILL_INPUT;
    }

    if abstolS.is_empty() {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSensSVtolerances", file!(), MSGCV_NULL_ABSTOLS);
        return CV_ILL_INPUT;
    }

    let mut atolmin = vec![ZERO; cv_mem.cv_Ns as usize];
    for is in 0..cv_mem.cv_Ns as usize {
        atolmin[is] = N_VMin(&abstolS[is]);
        if atolmin[is] < ZERO {
            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSensSVtolerances", file!(), MSGCV_BAD_ABSTOLS);
            return CV_ILL_INPUT;
        }
    }

    /* Copy tolerances into memory */
    cv_mem.cv_itolS = CV_SV;

    cv_mem.cv_reltolS = reltolS;

    if !cv_mem.cv_VabstolSMallocDone {
        cv_mem.cv_VabstolS = (0..cv_mem.cv_Ns as usize)
            .map(|_| N_VClone(&cv_mem.cv_tempv))
            .collect();
        cv_mem.cv_atolSmin0 = vec![SUNFALSE; cv_mem.cv_Ns as usize];
        cv_mem.cv_lrw += cv_mem.cv_Ns as i64 * cv_mem.cv_lrw1;
        cv_mem.cv_liw += cv_mem.cv_Ns as i64 * cv_mem.cv_liw1;
        cv_mem.cv_VabstolSMallocDone = SUNTRUE;
    }

    for is in 0..cv_mem.cv_Ns as usize {
        cv_mem.cv_atolSmin0[is] = atolmin[is] == ZERO;
    }

    for is in 0..cv_mem.cv_Ns as usize {
        cv_mem.cv_VabstolS[is].data.copy_from_slice(&abstolS[is].data);
    }

    CV_SUCCESS
}

pub fn CVodeSensEEtolerances(cv_mem: &mut CVodeMem) -> i32 {
    /* Was sensitivity initialized? */
    if !cv_mem.cv_SensMallocDone {
        cvProcessError(Some(cv_mem), CV_NO_SENS, line!(), "CVodeSensEEtolerances", file!(), MSGCV_NO_SENSI);
        return CV_NO_SENS;
    }

    cv_mem.cv_itolS = CV_EE;

    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

/*
 * CVodeQuadSensInit
 *
 */
pub fn CVodeQuadSensInit(
    cv_mem: &mut CVodeMem,
    fQS: Option<CVQuadSensRhsFn>,
    yQS0: &[NVector],
) -> i32 {
    /* Check if sensitivity analysis is active */
    if !cv_mem.cv_sensi {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeQuadSensInit", file!(), MSGCV_NO_SENSI);
        return CV_ILL_INPUT;
    }

    /* Check if yQS0 is non-null */
    if yQS0.is_empty() {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeQuadSensInit", file!(), MSGCV_NULL_YQS0);
        return CV_ILL_INPUT;
    }

    /* Allocate the vectors (using yQS0[0] as a template) */
    cvQuadSensAllocVectors(cv_mem, &yQS0[0]);

    /*----------------------------------------------
      All error checking is complete at this point
      -----------------------------------------------*/

    /* Set fQS */
    match fQS {
        None => {
            /* C: cv_fQS = cvQuadSensRhsInternalDQ, cv_fQS_data = cvode_mem;
               here cv_fQS stays None and cv_fQSDQ steers the dispatch to
               the internal DQ routine (module header) */
            cv_mem.cv_fQSDQ = SUNTRUE;
            cv_mem.cv_fQS = None;
        }
        Some(f) => {
            cv_mem.cv_fQSDQ = SUNFALSE;
            cv_mem.cv_fQS = Some(f);
        }
    }

    /* Initialize znQS[0] in the history array */
    for is in 0..cv_mem.cv_Ns as usize {
        cv_mem.cv_znQS[0][is].data.copy_from_slice(&yQS0[is].data);
    }

    /* Initialize all sensitivity related counters */
    cv_mem.cv_nfQSe = 0;
    cv_mem.cv_nfQeS = 0;
    cv_mem.cv_netfQS = 0;

    /* Quadrature sensitivities will be computed */
    cv_mem.cv_quadr_sensi = SUNTRUE;
    cv_mem.cv_QuadSensMallocDone = SUNTRUE;

    /* Sensitivity initialization was successful */
    CV_SUCCESS
}

/*
 * CVodeQuadSensReInit
 *
 */
pub fn CVodeQuadSensReInit(cv_mem: &mut CVodeMem, yQS0: &[NVector]) -> i32 {
    /* Check if sensitivity analysis is active */
    if !cv_mem.cv_sensi {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeQuadSensReInit", file!(), MSGCV_NO_SENSI);
        return CV_NO_SENS;
    }

    /* Was quadrature sensitivity initialized? */
    if !cv_mem.cv_QuadSensMallocDone {
        cvProcessError(Some(cv_mem), CV_NO_QUADSENS, line!(), "CVodeQuadSensReInit", file!(), MSGCV_NO_QUADSENSI);
        return CV_NO_QUADSENS;
    }

    /* Check if yQS0 is non-null */
    if yQS0.is_empty() {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeQuadSensReInit", file!(), MSGCV_NULL_YQS0);
        return CV_ILL_INPUT;
    }

    /*----------------------------------------------
      All error checking is complete at this point
      -----------------------------------------------*/

    /* Initialize znQS[0] in the history array */
    for is in 0..cv_mem.cv_Ns as usize {
        cv_mem.cv_znQS[0][is].data.copy_from_slice(&yQS0[is].data);
    }

    /* Initialize all sensitivity related counters */
    cv_mem.cv_nfQSe = 0;
    cv_mem.cv_nfQeS = 0;
    cv_mem.cv_netfQS = 0;

    /* Quadrature sensitivities will be computed */
    cv_mem.cv_quadr_sensi = SUNTRUE;

    /* Problem has been successfully re-initialized */
    CV_SUCCESS
}

/*
 * CVodeQuadSensSStolerances
 * CVodeQuadSensSVtolerances
 * CVodeQuadSensEEtolerances
 *
 * These functions specify the integration tolerances for quadrature
 * sensitivity variables. One of them MUST be called before the first
 * call to CVode IF these variables are included in the error test.
 *
 * CVodeQuadSensSStolerances specifies scalar relative and absolute tolerances.
 * CVodeQuadSensSVtolerances specifies scalar relative tolerance and a vector
 *   absolute tolerance for each quadrature sensitivity vector (a potentially
 *   different absolute tolerance for each vector component).
 * CVodeQuadSensEEtolerances specifies that tolerances for sensitivity variables
 *   should be estimated from those provided for the quadrature variables.
 *   In this case, tolerances for the quadrature variables must be
 *   specified through a call to one of CVodeQuad**tolerances.
 */

pub fn CVodeQuadSensSStolerances(cv_mem: &mut CVodeMem, reltolQS: f64, abstolQS: &[f64]) -> i32 {
    /* Check if sensitivity was initialized */
    if !cv_mem.cv_SensMallocDone {
        cvProcessError(Some(cv_mem), CV_NO_SENS, line!(), "CVodeQuadSensSStolerances", file!(), MSGCV_NO_SENSI);
        return CV_NO_SENS;
    }

    /* Check if quadrature sensitivity was initialized? */
    if !cv_mem.cv_QuadSensMallocDone {
        cvProcessError(Some(cv_mem), CV_NO_QUADSENS, line!(), "CVodeQuadSensSStolerances", file!(), MSGCV_NO_QUADSENSI);
        /* (C returns CV_NO_QUAD here, not CV_NO_QUADSENS) */
        return CV_NO_QUAD;
    }

    /* Test user-supplied tolerances */
    if reltolQS < ZERO {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeQuadSensSStolerances", file!(), MSGCV_BAD_RELTOLQS);
        return CV_ILL_INPUT;
    }

    if abstolQS.is_empty() {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeQuadSensSStolerances", file!(), MSGCV_NULL_ABSTOLQS);
        return CV_ILL_INPUT;
    }

    for is in 0..cv_mem.cv_Ns as usize {
        if abstolQS[is] < ZERO {
            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeQuadSensSStolerances", file!(), MSGCV_BAD_ABSTOLQS);
            return CV_ILL_INPUT;
        }
    }

    /* Copy tolerances into memory */
    cv_mem.cv_itolQS = CV_SS;

    cv_mem.cv_reltolQS = reltolQS;

    if !cv_mem.cv_SabstolQSMallocDone {
        cv_mem.cv_SabstolQS = vec![ZERO; cv_mem.cv_Ns as usize];
        cv_mem.cv_atolQSmin0 = vec![SUNFALSE; cv_mem.cv_Ns as usize];
        cv_mem.cv_lrw += cv_mem.cv_Ns as i64;
        cv_mem.cv_SabstolQSMallocDone = SUNTRUE;
    }

    for is in 0..cv_mem.cv_Ns as usize {
        cv_mem.cv_SabstolQS[is] = abstolQS[is];
        cv_mem.cv_atolQSmin0[is] = abstolQS[is] == ZERO;
    }

    CV_SUCCESS
}

pub fn CVodeQuadSensSVtolerances(
    cv_mem: &mut CVodeMem,
    reltolQS: f64,
    abstolQS: &[NVector],
) -> i32 {
    /* check if sensitivity was initialized */
    if !cv_mem.cv_SensMallocDone {
        cvProcessError(Some(cv_mem), CV_NO_SENS, line!(), "CVodeQuadSensSVtolerances", file!(), MSGCV_NO_SENSI);
        return CV_NO_SENS;
    }

    /* Check if quadrature sensitivity was initialized? */
    if !cv_mem.cv_QuadSensMallocDone {
        cvProcessError(Some(cv_mem), CV_NO_QUADSENS, line!(), "CVodeQuadSensSVtolerances", file!(), MSGCV_NO_QUADSENSI);
        /* (C returns CV_NO_QUAD here, not CV_NO_QUADSENS) */
        return CV_NO_QUAD;
    }

    /* Test user-supplied tolerances */
    if reltolQS < ZERO {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeQuadSensSVtolerances", file!(), MSGCV_BAD_RELTOLQS);
        return CV_ILL_INPUT;
    }

    if abstolQS.is_empty() {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeQuadSensSVtolerances", file!(), MSGCV_NULL_ABSTOLQS);
        return CV_ILL_INPUT;
    }

    let mut atolmin = vec![ZERO; cv_mem.cv_Ns as usize];
    for is in 0..cv_mem.cv_Ns as usize {
        atolmin[is] = N_VMin(&abstolQS[is]);
        if atolmin[is] < ZERO {
            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeQuadSensSVtolerances", file!(), MSGCV_BAD_ABSTOLQS);
            return CV_ILL_INPUT;
        }
    }

    /* Copy tolerances into memory */
    cv_mem.cv_itolQS = CV_SV;

    cv_mem.cv_reltolQS = reltolQS;

    if !cv_mem.cv_VabstolQSMallocDone {
        cv_mem.cv_VabstolQS = (0..cv_mem.cv_Ns as usize)
            .map(|_| N_VClone(&cv_mem.cv_tempvQ))
            .collect();
        cv_mem.cv_atolQSmin0 = vec![SUNFALSE; cv_mem.cv_Ns as usize];
        cv_mem.cv_lrw += cv_mem.cv_Ns as i64 * cv_mem.cv_lrw1Q;
        cv_mem.cv_liw += cv_mem.cv_Ns as i64 * cv_mem.cv_liw1Q;
        cv_mem.cv_VabstolQSMallocDone = SUNTRUE;
    }

    for is in 0..cv_mem.cv_Ns as usize {
        cv_mem.cv_atolQSmin0[is] = atolmin[is] == ZERO;
    }

    for is in 0..cv_mem.cv_Ns as usize {
        cv_mem.cv_VabstolQS[is].data.copy_from_slice(&abstolQS[is].data);
    }

    CV_SUCCESS
}

pub fn CVodeQuadSensEEtolerances(cv_mem: &mut CVodeMem) -> i32 {
    /* check if sensitivity was initialized */
    if !cv_mem.cv_SensMallocDone {
        cvProcessError(Some(cv_mem), CV_NO_SENS, line!(), "CVodeQuadSensEEtolerances", file!(), MSGCV_NO_SENSI);
        return CV_NO_SENS;
    }

    /* Check if quadrature sensitivity was initialized? */
    if !cv_mem.cv_QuadSensMallocDone {
        cvProcessError(Some(cv_mem), CV_NO_QUADSENS, line!(), "CVodeQuadSensEEtolerances", file!(), MSGCV_NO_QUADSENSI);
        /* (C returns CV_NO_QUAD here, not CV_NO_QUADSENS) */
        return CV_NO_QUAD;
    }

    cv_mem.cv_itolQS = CV_EE;

    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

/*
 * CVodeSensToggleOff
 *
 * CVodeSensToggleOff deactivates sensitivity calculations.
 * It does NOT deallocate sensitivity-related memory.
 */
pub fn CVodeSensToggleOff(cv_mem: &mut CVodeMem) -> i32 {
    /* Disable sensitivities */
    cv_mem.cv_sensi = SUNFALSE;
    cv_mem.cv_quadr_sensi = SUNFALSE;

    CV_SUCCESS
}

/*-----------------------------------------------------------------*/

/*
 * CVodeRootInit
 *
 * CVodeRootInit initializes a rootfinding problem to be solved
 * during the integration of the ODE system.  It loads the root
 * function pointer and the number of root functions, and allocates
 * workspace memory.  The return value is CV_SUCCESS = 0 if no errors
 * occurred, or a negative value otherwise.
 */
pub fn CVodeRootInit(cv_mem: &mut CVodeMem, nrtfn: i32, g: Option<CVRootFn>) -> i32 {
    let nrt = if nrtfn < 0 { 0 } else { nrtfn };

    /* If rerunning CVodeRootInit() with a different number of root
       functions (changing number of gfun components), then free
       currently held memory resources */
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

    /* If CVodeRootInit() was called with nrtfn == 0, then set cv_nrtfn to
       zero and cv_gfun to NULL before returning */
    if nrt == 0 {
        cv_mem.cv_nrtfn = nrt;
        cv_mem.cv_gfun = None;
        return CV_SUCCESS;
    }

    /* If rerunning CVodeRootInit() with the same number of root functions
       (not changing number of gfun components), then check if the root
       function argument has changed */
    /* If g != NULL then return as currently reserved memory resources
       will suffice */
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

    /* Set default values for rootdir (both directions) */
    cv_mem.cv_rootdir = vec![0; n];

    /* Set default values for gactive (all active) */
    cv_mem.cv_gactive = vec![SUNTRUE; n];

    cv_mem.cv_lrw += 3 * nrt as i64;
    cv_mem.cv_liw += 3 * nrt as i64;

    CV_SUCCESS
}

// ===================== END PART 1 (cvodes.c:1-2912) =====================
// PART 2 (CVode driver, GetDky/Get* extraction, CVodeFree, cvInitialSetup,
// cvHin, alloc/free helpers — cvodes.c:2913-6244) is appended below by the
// next agent in the chain.

/* -----------------------------------------------------------------
 * Translated from src/cvodes/cvodes.c (CVODES 7.7.0) — PART 2
 * (cvodes.c lines 2913-6244): the CVode main driver (with
 * quadrature / sensitivity / quadrature-sensitivity interleaving),
 * CVodeComputeState[Sens[1]], the CVodeGetDky / GetQuad* / GetSens*
 * / GetQuadSens* extraction family, the CVodeFree / CVodeQuadFree /
 * CVodeSensFree / CVodeQuadSensFree deallocation family, the
 * cv*AllocVectors / cv*FreeVectors helpers, cvInitialSetup and the
 * initial step-size machinery (cvHin, cvUpperBoundH0, cvYddNorm).
 * cvStep itself (cvodes.c:5874-6244) is step machinery and is
 * appended with PART 3 (see the end-of-part marker below).
 *
 * The donor translation crates/cvode_rs/src/cvode.rs is reused
 * verbatim wherever the CVODES C matches the CVODE C; the
 * quadrature / sensitivity deltas follow the CVODES source.
 * Documented adaptations (donor precedents unless noted):
 *  - In C, cv_mem->cv_y aliases the user's yout for the whole
 *    duration of CVode ("cv_mem->cv_y = yout" at entry). Here cv_y
 *    stays owned by CVodeMem and is copied into yout at EVERY
 *    return path that the C reaches with data in yout — including
 *    all three CV_ROOT_RETURN exits, where the C relies on the
 *    rootfinder having interpolated into cv_y==yout (CLAUDE.md
 *    rule 5).
 *  - cvode_mem == NULL and tret == NULL checks vanish (&mut refs).
 *  - cvCheckNvector is dropped entirely (donor precedent): the
 *    concrete serial NVector implements every required operation,
 *    so the ops-table guard is vacuously true.
 *  - malloc-failure branches vanish: the cv*AllocVectors helpers
 *    are infallible (Vec/Box allocation aborts on OOM) and the
 *    corresponding N_VDestroy error-unwind ladders disappear
 *    (RAII); the SUNTRUE/SUNFALSE returns become ().
 *  - Fused vector ops (N_VScaleVectorArray, N_VLinearSumVectorArray,
 *    N_VLinearCombination, N_VWrmsNormVectorArray) are expanded to
 *    per-vector loops; their CV_VECTOROP_ERR branches vanish and no
 *    cv_cvals/cv_Xvecs scratch arrays exist (pinned decision).
 *  - The cv_linit/cv_lsetup/cv_lsolve function-pointer dispatch
 *    becomes the LsModule take()-dispatch helpers below (donor
 *    pattern; the CVODES lsolve additionally forwards the
 *    staggered1 sensitivity index for the ewtS[is] weight).
 *  - cv_fQS dispatch: when cv_fQSDQ is SUNTRUE the C stores
 *    cvQuadSensRhsInternalDQ in cv_fQS with cv_fQS_data = cv_mem;
 *    here cv_fQS_dispatch branches on the flag (pinned decision 5
 *    in cvodes_impl.rs; the DQ routine itself is ported with the
 *    sensitivity RHS wrappers in PART 3).
 * -----------------------------------------------------------------*/
use crate::sundials_math::*;

/* Remaining CVODE private constants (see the PART 1 header note) */
const TINY: f64 = 1.0e-10;
const PT1: f64 = 0.1;
const POINT2: f64 = 0.2;
const FOURTH: f64 = 0.25;
const HALF: f64 = 0.5;
const PT9: f64 = 0.9;
const ONEPT5: f64 = 1.5;
const TWO: f64 = 2.0;
const THREE: f64 = 3.0;
const FOUR: f64 = 4.0;
const FIVE: f64 = 5.0;
const TWELVE: f64 = 12.0;
const HUNDRED: f64 = 100.0;

/* Algorithmic constants: CVodeGetDky and cvStep */
const FUZZ_FACTOR: f64 = 100.0;
/* cvHin */
const HLB_FACTOR: f64 = 100.0;
const HUB_FACTOR: f64 = 0.1;
const H_BIAS: f64 = HALF;
const MAX_ITERS: i32 = 4;

/*
 * -----------------------------------------------------------------
 * Linear solver module dispatch
 * (In C these are the cv_linit/cv_lsetup/cv_lsolve function
 *  pointers; the module is taken out of CVodeMem during the call so
 *  its routine can borrow the integrator memory mutably. Donor
 *  pattern. cv_lreinit_dispatch landed with PART 1; the shared
 *  cv_has_lsetup / cv_lsetup_dispatch / cv_lsolve_dispatch helpers
 *  live in cvodes_nls.rs next to their C callers (the CVODES lsolve
 *  additionally forwards the STAGGERED1 sensitivity index selecting
 *  the ewtS[is] weight, mirroring the C cv_lsolve weight argument);
 *  only the linit dispatch is needed here.)
 * -----------------------------------------------------------------
 */
use crate::cvodes_nls::cv_has_lsetup;

pub(crate) fn cv_linit_dispatch(cv_mem: &mut CVodeMem) -> i32 {
    let mut lmem = std::mem::take(&mut cv_mem.cv_lmem);
    let ier = match &mut lmem {
        LsModule::None => 0,
        LsModule::Ls(ls) => crate::cvodes_ls::cvLsInitialize(cv_mem, ls),
        LsModule::Diag(dm) => crate::cvodes_diag::cvDiagInit(cv_mem, dm),
    };
    cv_mem.cv_lmem = lmem;
    ier
}

/*
 * -----------------------------------------------------------------
 * Error weight dispatch helpers (donor pattern)
 *
 * The cvEwtSet / cvQuadEwtSet / cvSensEwtSet / cvQuadSensEwtSet
 * core routines are defined in PART 3 (cvodes.c:9260+); the
 * take()-based wrappers below detach the target weight vectors from
 * CVodeMem so those routines can borrow the integrator memory.
 * -----------------------------------------------------------------
 */

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

/* Apply cvQuadEwtSet to the integrator's own ewtQ vector
   (C: cvQuadEwtSet(cv_mem, cv_mem->cv_znQ[0], cv_mem->cv_ewtQ)). */
pub(crate) fn cvQuadEwtSet_apply_to_ewtQ(cv_mem: &mut CVodeMem) -> i32 {
    let mut w = std::mem::take(&mut cv_mem.cv_ewtQ);
    let flag = cvQuadEwtSet(cv_mem, &cv_mem.cv_znQ[0], &mut w);
    cv_mem.cv_ewtQ = w;
    flag
}

/* Apply cvSensEwtSet to the integrator's own ewtS vectors
   (C: cvSensEwtSet(cv_mem, cv_mem->cv_znS[0], cv_mem->cv_ewtS)). */
pub(crate) fn cvSensEwtSet_apply_to_ewtS(cv_mem: &mut CVodeMem) -> i32 {
    let yScur = std::mem::take(&mut cv_mem.cv_znS[0]);
    let mut wS = std::mem::take(&mut cv_mem.cv_ewtS);
    let flag = cvSensEwtSet(cv_mem, &yScur, &mut wS);
    cv_mem.cv_znS[0] = yScur;
    cv_mem.cv_ewtS = wS;
    flag
}

/* Apply cvQuadSensEwtSet to the integrator's own ewtQS vectors
   (C: cvQuadSensEwtSet(cv_mem, cv_mem->cv_znQS[0], cv_mem->cv_ewtQS)). */
pub(crate) fn cvQuadSensEwtSet_apply_to_ewtQS(cv_mem: &mut CVodeMem) -> i32 {
    let yQScur = std::mem::take(&mut cv_mem.cv_znQS[0]);
    let mut wQS = std::mem::take(&mut cv_mem.cv_ewtQS);
    let flag = cvQuadSensEwtSet(cv_mem, &yQScur, &mut wQS);
    cv_mem.cv_znQS[0] = yQScur;
    cv_mem.cv_ewtQS = wQS;
    flag
}

/*
 * cv_fQS_dispatch
 *
 * Dispatch of the quadrature-sensitivity RHS: the user's fQS or the
 * internal DQ routine cvQuadSensRhsInternalDQ (PART 3), selected by
 * cv_fQSDQ. In C the selection is done once (cv_fQS function
 * pointer, cv_fQS_data = cv_mem); the fQS counter increments stay
 * at the call sites exactly as in C.
 */
pub(crate) fn cv_fQS_dispatch(
    cv_mem: &mut CVodeMem,
    t: f64,
    y: &NVector,
    yS: &[NVector],
    yQdot: &NVector,
    yQSdot: &mut [NVector],
    tmp: &mut NVector,
    tmpQ: &mut NVector,
) -> i32 {
    let Ns = cv_mem.cv_Ns;
    if cv_mem.cv_fQSDQ {
        cvQuadSensRhsInternalDQ(cv_mem, Ns, t, y, yS, yQdot, yQSdot, tmp, tmpQ)
    } else {
        let fQS = cv_mem.cv_fQS.unwrap();
        fQS(Ns, t, y, yS, yQdot, yQSdot, &mut cv_mem.cv_user_data, tmp, tmpQ)
    }
}

/*
 * -----------------------------------------------------------------
 * Main solver function
 * -----------------------------------------------------------------
 */

/*
 * CVode
 *
 * This routine is the main driver of the CVODES package.
 *
 * It integrates over a time interval defined by the user, by calling
 * cvStep to do internal time steps.
 *
 * The first time that CVode is called for a successfully initialized
 * problem, it computes a tentative initial step size h.
 *
 * CVode supports two modes, specified by itask: CV_NORMAL, CV_ONE_STEP.
 * In the CV_NORMAL mode, the solver steps until it reaches or passes tout
 * and then interpolates to obtain y(tout).
 * In the CV_ONE_STEP mode, it takes one internal step and returns.
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

    /* Check if cvode_mem was allocated */
    if !cv_mem.cv_MallocDone {
        cvProcessError(Some(cv_mem), CV_NO_MALLOC, line!(), "CVode", file!(), MSGCV_NO_MALLOC);
        return CV_NO_MALLOC;
    }

    /* Check for yout != NULL */
    if yout.is_empty() {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVode", file!(), MSGCV_YOUT_NULL);
        return CV_ILL_INPUT;
    }

    /* Check for valid itask */
    if itask != CV_NORMAL && itask != CV_ONE_STEP {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVode", file!(), MSGCV_BAD_ITASK);
        return CV_ILL_INPUT;
    }

    /*
     * ----------------------------------------
     * 2. Initializations performed only at
     *    the first step (nst=0):
     *    - initial setup
     *    - initialize Nordsieck history array
     *    - compute initial step size
     *    - check for approach to tstop
     *    - check for approach to a root
     *    Or initializations performed after
     *    resizing the integrator
     *    - check constraints
     *    - initialize linear solver
     *    - initialize nonlinear solver
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

        /*
         * Call f at (t0,y0), set zn[1] = y'(t0).
         * If computing any quadratures, call fQ at (t0,y0), set znQ[1] = yQ'(t0)
         * If computing sensitivities, call fS at (t0,y0,yS0), set znS[1][is] = yS'(t0), is=1,...,Ns.
         * If computing quadr. sensi., call fQS at (t0,y0,yS0), set znQS[1][is] = yQS'(t0), is=1,...,Ns.
         */

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

        if cv_mem.cv_quadr {
            let fQ = cv_mem.cv_fQ.unwrap();
            let tn = cv_mem.cv_tn;
            let retval = {
                let CVodeMem { cv_zn, cv_znQ, cv_user_data, .. } = cv_mem;
                fQ(tn, &cv_zn[0], &mut cv_znQ[1], cv_user_data)
            };
            cv_mem.cv_nfQe += 1;
            if retval < 0 {
                cvProcessError(Some(cv_mem), CV_QRHSFUNC_FAIL, line!(), "CVode", file!(),
                    &format!("At t = {}, the quadrature right-hand side routine failed in an unrecoverable manner.", cv_mem.cv_tn));
                return CV_QRHSFUNC_FAIL;
            }
            if retval > 0 {
                cvProcessError(Some(cv_mem), CV_FIRST_QRHSFUNC_ERR, line!(), "CVode", file!(), MSGCV_QRHSFUNC_FIRST);
                return CV_FIRST_QRHSFUNC_ERR;
            }
        }

        if cv_mem.cv_sensi {
            let tn = cv_mem.cv_tn;
            let zn0 = std::mem::take(&mut cv_mem.cv_zn[0]);
            let zn1 = std::mem::take(&mut cv_mem.cv_zn[1]);
            let znS0 = std::mem::take(&mut cv_mem.cv_znS[0]);
            let mut znS1 = std::mem::take(&mut cv_mem.cv_znS[1]);
            let mut tempv = std::mem::take(&mut cv_mem.cv_tempv);
            let mut ftemp = std::mem::take(&mut cv_mem.cv_ftemp);
            let retval =
                cvSensRhsWrapper(cv_mem, tn, &zn0, &zn1, &znS0, &mut znS1, &mut tempv, &mut ftemp);
            cv_mem.cv_zn[0] = zn0;
            cv_mem.cv_zn[1] = zn1;
            cv_mem.cv_znS[0] = znS0;
            cv_mem.cv_znS[1] = znS1;
            cv_mem.cv_tempv = tempv;
            cv_mem.cv_ftemp = ftemp;
            if retval < 0 {
                cvProcessError(Some(cv_mem), CV_SRHSFUNC_FAIL, line!(), "CVode", file!(),
                    &format!("At t = {}, the sensitivity right-hand side routine failed in an unrecoverable manner.", cv_mem.cv_tn));
                return CV_SRHSFUNC_FAIL;
            }
            if retval > 0 {
                cvProcessError(Some(cv_mem), CV_FIRST_SRHSFUNC_ERR, line!(), "CVode", file!(), MSGCV_SRHSFUNC_FIRST);
                return CV_FIRST_SRHSFUNC_ERR;
            }
        }

        if cv_mem.cv_quadr_sensi {
            let tn = cv_mem.cv_tn;
            let zn0 = std::mem::take(&mut cv_mem.cv_zn[0]);
            let znS0 = std::mem::take(&mut cv_mem.cv_znS[0]);
            let znQ1 = std::mem::take(&mut cv_mem.cv_znQ[1]);
            let mut znQS1 = std::mem::take(&mut cv_mem.cv_znQS[1]);
            let mut tempv = std::mem::take(&mut cv_mem.cv_tempv);
            let mut tempvQ = std::mem::take(&mut cv_mem.cv_tempvQ);
            let retval =
                cv_fQS_dispatch(cv_mem, tn, &zn0, &znS0, &znQ1, &mut znQS1, &mut tempv, &mut tempvQ);
            cv_mem.cv_zn[0] = zn0;
            cv_mem.cv_znS[0] = znS0;
            cv_mem.cv_znQ[1] = znQ1;
            cv_mem.cv_znQS[1] = znQS1;
            cv_mem.cv_tempv = tempv;
            cv_mem.cv_tempvQ = tempvQ;
            cv_mem.cv_nfQSe += 1;
            if retval < 0 {
                cvProcessError(Some(cv_mem), CV_QSRHSFUNC_FAIL, line!(), "CVode", file!(),
                    &format!("At t = {}, the quadrature sensitivity right-hand side routine failed in an unrecoverable manner.", cv_mem.cv_tn));
                return CV_QSRHSFUNC_FAIL;
            }
            if retval > 0 {
                cvProcessError(Some(cv_mem), CV_FIRST_QSRHSFUNC_ERR, line!(), "CVode", file!(), MSGCV_QSRHSFUNC_FIRST);
                return CV_FIRST_QSRHSFUNC_ERR;
            }
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

        /*
         * Scale zn[1] by h.
         * If computing any quadratures, scale znQ[1] by h.
         * If computing sensitivities,  scale znS[1][is] by h.
         * If computing quadrature sensitivities,  scale znQS[1][is] by h.
         */

        cv_mem.cv_hscale = cv_mem.cv_h;
        cv_mem.cv_h0u = cv_mem.cv_h;
        cv_mem.cv_hprime = cv_mem.cv_h;

        let h = cv_mem.cv_h;
        cv_mem.cv_zn[1].scale_inplace(h);

        if cv_mem.cv_quadr {
            cv_mem.cv_znQ[1].scale_inplace(h);
        }

        if cv_mem.cv_sensi {
            for is in 0..cv_mem.cv_Ns as usize {
                cv_mem.cv_znS[1][is].scale_inplace(h);
            }
        }

        if cv_mem.cv_quadr_sensi {
            for is in 0..cv_mem.cv_Ns as usize {
                cv_mem.cv_znQS[1][is].scale_inplace(h);
            }
        }

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
        if !cv_mem.cv_constraints.is_empty() {
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

        /* Initialize the nonlinear solver (must occur after linear solver is
           initialized) so the lsetup and lsolve pointers have been set */
        let ier = crate::cvodes_nls::cvNlsInit(cv_mem);
        if ier != 0 {
            cvProcessError(Some(cv_mem), CV_NLS_INIT_FAIL, line!(), "CVode", file!(), MSGCV_NLS_INIT_FAIL);
            return CV_NLS_INIT_FAIL;
        }
    }

    /*
     * ------------------------------------------------------
     * 3. At following steps, perform stop tests:
     *    - check for root in last step
     *    - check if we passed tstop
     *    - check if we passed tout (NORMAL mode)
     *    - check if current tn was returned (ONE_STEP mode)
     *    - check if we are close to tstop
     *      (adjust step size if needed)
     * -------------------------------------------------------
     */

    if cv_mem.cv_nst > 0 {
        /* Estimate an infinitesimal time interval to be used as
           a roundoff for time quantities (based on current time
           and step size) */
        let troundoff =
            FUZZ_FACTOR * cv_mem.cv_uround * (SUNRabs(cv_mem.cv_tn) + SUNRabs(cv_mem.cv_h));

        /* First, check for a root in the last step taken, other than the
           last root found, if any.  If itask = CV_ONE_STEP and y(tn) was not
           returned because of an intervening root, return y(tn) now.     */
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

            /* If tn is distinct from tretlast (within roundoff),
               check remaining interval for roots */
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
                    /* g failed */
                    cvProcessError(Some(cv_mem), CV_RTFUNC_FAIL, line!(), "CVode", file!(),
                        &format!("At t = {}, the rootfinding routine failed in an unrecoverable manner.", cv_mem.cv_tlo));
                    return CV_RTFUNC_FAIL;
                }
            }
        } /* end of root stop check */

        /* Test for tn at tstop or near tstop */
        if cv_mem.cv_tstopset {
            /* Test for tn at tstop */
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
     *
     *    4.1. check for errors (too many steps, too much
     *         accuracy requested, step size too small)
     *    4.2. take a new step (call cvStep)
     *    4.3. stop on error
     *    4.4. perform stop tests:
     *         - check for root in last step
     *         - check if tout was passed
     *         - check if close to tstop
     *         - check if in ONE_STEP mode (must return)
     * --------------------------------------------------
     */

    let mut nstloc: i64 = 0;
    let istate;
    loop {
        cv_mem.cv_next_h = cv_mem.cv_h;
        cv_mem.cv_next_q = cv_mem.cv_q;

        /* Reset and check ewt, ewtQ, ewtS */
        if cv_mem.cv_nst > 0 {
            let ier = cv_efun_apply_to_ewt(cv_mem);
            if ier != 0 {
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

            if cv_mem.cv_quadr && cv_mem.cv_errconQ {
                let ier = cvQuadEwtSet_apply_to_ewtQ(cv_mem);
                if ier != 0 {
                    cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVode", file!(),
                        &format!("At t = {}, a component of ewtQ has become <= 0.", cv_mem.cv_tn));
                    istate = CV_ILL_INPUT;
                    cv_mem.cv_tretlast = cv_mem.cv_tn;
                    *tret = cv_mem.cv_tn;
                    yout.data.copy_from_slice(&cv_mem.cv_zn[0].data);
                    break;
                }
            }

            if cv_mem.cv_sensi {
                let ier = cvSensEwtSet_apply_to_ewtS(cv_mem);
                if ier != 0 {
                    cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVode", file!(),
                        &format!("At t = {}, a component of ewtS has become <= 0.", cv_mem.cv_tn));
                    istate = CV_ILL_INPUT;
                    cv_mem.cv_tretlast = cv_mem.cv_tn;
                    *tret = cv_mem.cv_tn;
                    yout.data.copy_from_slice(&cv_mem.cv_zn[0].data);
                    break;
                }
            }

            if cv_mem.cv_quadr_sensi && cv_mem.cv_errconQS {
                let ier = cvQuadSensEwtSet_apply_to_ewtQS(cv_mem);
                if ier != 0 {
                    cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVode", file!(),
                        &format!("At t = {}, a component of ewtQS has become <= 0.", cv_mem.cv_tn));
                    istate = CV_ILL_INPUT;
                    cv_mem.cv_tretlast = cv_mem.cv_tn;
                    *tret = cv_mem.cv_tn;
                    yout.data.copy_from_slice(&cv_mem.cv_zn[0].data);
                    break;
                }
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
        let mut nrm = N_VWrmsNorm(&cv_mem.cv_zn[0], &cv_mem.cv_ewt);
        if cv_mem.cv_quadr && cv_mem.cv_errconQ {
            nrm = cvQuadUpdateNorm(cv_mem, nrm, &cv_mem.cv_znQ[0], &cv_mem.cv_ewtQ);
        }
        if cv_mem.cv_sensi && cv_mem.cv_errconS {
            nrm = cvSensUpdateNorm(cv_mem, nrm, &cv_mem.cv_znS[0], &cv_mem.cv_ewtS);
        }
        if cv_mem.cv_quadr_sensi && cv_mem.cv_errconQS {
            nrm = cvQuadSensUpdateNorm(cv_mem, nrm, &cv_mem.cv_znQS[0], &cv_mem.cv_ewtQS);
        }
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

            /* If we are at the end of the first step and we still have
             * some event functions that are inactive, issue a warning
             * as this may indicate a user error in the implementation
             * of the root function. */

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

    /* Load optional output */
    if cv_mem.cv_sensi && cv_mem.cv_ism == CV_STAGGERED1 {
        cv_mem.cv_nniS = 0;
        cv_mem.cv_nnfS = 0;
        cv_mem.cv_ncfnS = 0;
        for is in 0..cv_mem.cv_Ns as usize {
            cv_mem.cv_nniS += cv_mem.cv_nniS1[is];
            cv_mem.cv_nnfS += cv_mem.cv_nnfS1[is];
            cv_mem.cv_ncfnS += cv_mem.cv_ncfnS1[is];
        }
    }

    istate
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
 * CVodeComputeStateSens
 *
 * Computes yS based on the current prediction and given correction.
 */
pub fn CVodeComputeStateSens(cv_mem: &CVodeMem, ycorS: &[NVector], yS: &mut [NVector]) -> i32 {
    for is in 0..cv_mem.cv_Ns as usize {
        N_VLinearSum(ONE, &cv_mem.cv_znS[0][is], ONE, &ycorS[is], &mut yS[is]);
    }
    CV_SUCCESS
}

/*
 * CVodeComputeStateSens1
 *
 * Computes yS[idx] based on the current prediction and given correction.
 */
pub fn CVodeComputeStateSens1(
    cv_mem: &CVodeMem,
    idx: i32,
    ycorS1: &NVector,
    yS1: &mut NVector,
) -> i32 {
    N_VLinearSum(ONE, &cv_mem.cv_znS[0][idx as usize], ONE, ycorS1, yS1);
    CV_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Interpolated output and extraction functions
 * -----------------------------------------------------------------
 */

/*
 * CVodeGetDky
 *
 * This routine computes the k-th derivative of the interpolating
 * polynomial at the time t and stores the result in the vector dky.
 * The formula is:
 *         q
 *  dky = SUM c(j,k) * (t - tn)^(j-k) * h^(-j) * zn[j] ,
 *        j=k
 * where c(j,k) = j*(j-1)*...*(j-k+1), q is the current order, and
 * zn[j] is the j-th column of the Nordsieck history array.
 *
 * This function is called by CVode with k = 0 and t = tout, but
 * may also be called directly by the user.
 */
pub fn CVodeGetDky(cv_mem: &CVodeMem, t: f64, k: i32, dky: &mut NVector) -> i32 {
    /* Check all inputs for legality */

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
 * CVodeGetQuad
 *
 * This routine extracts quadrature solution into yQout at the
 * time which CVode returned the solution.
 * This is just a wrapper that calls CVodeGetQuadDky with k=0.
 */
pub fn CVodeGetQuad(cv_mem: &CVodeMem, tret: &mut f64, yQout: &mut NVector) -> i32 {
    *tret = cv_mem.cv_tretlast;

    CVodeGetQuadDky(cv_mem, cv_mem.cv_tretlast, 0, yQout)
}

/*
 * CVodeGetQuadDky
 *
 * CVodeQuadDky computes the kth derivative of the yQ function at
 * time t, where tn-hu <= t <= tn, tn denotes the current
 * internal time reached, and hu is the last internal step size
 * successfully used by the solver. The user may request
 * k=0, 1, ..., qu, where qu is the current order.
 * The derivative vector is returned in dky. This vector
 * must be allocated by the caller. It is only legal to call this
 * function after a successful return from CVode with quadrature
 * computation enabled.
 */
pub fn CVodeGetQuadDky(cv_mem: &CVodeMem, t: f64, k: i32, dkyQ: &mut NVector) -> i32 {
    /* Check all inputs for legality */

    if cv_mem.cv_quadr != SUNTRUE {
        cvProcessError(Some(cv_mem), CV_NO_QUAD, line!(), "CVodeGetQuadDky", file!(), MSGCV_NO_QUAD);
        return CV_NO_QUAD;
    }

    if dkyQ.is_empty() {
        cvProcessError(Some(cv_mem), CV_BAD_DKY, line!(), "CVodeGetQuadDky", file!(), MSGCV_NULL_DKY);
        return CV_BAD_DKY;
    }

    if k < 0 || k > cv_mem.cv_q {
        cvProcessError(Some(cv_mem), CV_BAD_K, line!(), "CVodeGetQuadDky", file!(), MSGCV_BAD_K);
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
        cvProcessError(Some(cv_mem), CV_BAD_T, line!(), "CVodeGetQuadDky", file!(),
            &format!("Illegal value for t. t = {} is not between tcur - hold = {} and tcur = {}",
                     t, cv_mem.cv_tn - cv_mem.cv_hu, cv_mem.cv_tn));
        return CV_BAD_T;
    }

    /* Sum the differentiated interpolating polynomial */
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
        let znQj = &cv_mem.cv_znQ[j as usize];
        if first {
            for (d, z) in dkyQ.data.iter_mut().zip(&znQj.data) {
                *d = c * *z;
            }
            first = false;
        } else {
            for (d, z) in dkyQ.data.iter_mut().zip(&znQj.data) {
                *d += c * *z;
            }
        }
        j -= 1;
    }

    if k == 0 {
        return CV_SUCCESS;
    }
    let r = SUNRpowerI(cv_mem.cv_h, -k);
    dkyQ.scale_inplace(r);

    CV_SUCCESS
}

/*
 * CVodeGetSens
 *
 * This routine extracts sensitivity solution into ySout at the
 * time at which CVode returned the solution.
 * This is just a wrapper that calls CVodeSensDky with k=0.
 */
pub fn CVodeGetSens(cv_mem: &CVodeMem, tret: &mut f64, ySout: &mut [NVector]) -> i32 {
    *tret = cv_mem.cv_tretlast;

    CVodeGetSensDky(cv_mem, cv_mem.cv_tretlast, 0, ySout)
}

/*
 * CVodeGetSens1
 *
 * This routine extracts the is-th sensitivity solution into ySout
 * at the time at which CVode returned the solution.
 * This is just a wrapper that calls CVodeSensDky1 with k=0.
 */
pub fn CVodeGetSens1(cv_mem: &CVodeMem, tret: &mut f64, is: i32, ySout: &mut NVector) -> i32 {
    *tret = cv_mem.cv_tretlast;

    CVodeGetSensDky1(cv_mem, cv_mem.cv_tretlast, 0, is, ySout)
}

/*
 * CVodeGetSensDky
 *
 * If the user calls directly CVodeSensDky then s must be allocated
 * prior to this call. When CVodeSensDky is called by
 * CVodeGetSens, only ier=CV_SUCCESS, ier=CV_NO_SENS, or
 * ier=CV_BAD_T are possible.
 */
pub fn CVodeGetSensDky(cv_mem: &CVodeMem, t: f64, k: i32, dkyS: &mut [NVector]) -> i32 {
    let mut ier = CV_SUCCESS;

    if dkyS.is_empty() {
        cvProcessError(Some(cv_mem), CV_BAD_DKY, line!(), "CVodeGetSensDky", file!(), MSGCV_NULL_DKYA);
        return CV_BAD_DKY;
    }

    for is in 0..cv_mem.cv_Ns as usize {
        ier = CVodeGetSensDky1(cv_mem, t, k, is as i32, &mut dkyS[is]);
        if ier != CV_SUCCESS {
            break;
        }
    }

    ier
}

/*
 * CVodeGetSensDky1
 *
 * CVodeSensDky1 computes the kth derivative of the yS[is] function at
 * time t, where tn-hu <= t <= tn, tn denotes the current
 * internal time reached, and hu is the last internal step size
 * successfully used by the solver. The user may request
 * is=0, 1, ..., Ns-1 and k=0, 1, ..., qu, where qu is the current
 * order. The derivative vector is returned in dky. This vector
 * must be allocated by the caller. It is only legal to call this
 * function after a successful return from CVode with sensitivity
 * computation enabled.
 */
pub fn CVodeGetSensDky1(cv_mem: &CVodeMem, t: f64, k: i32, is: i32, dkyS: &mut NVector) -> i32 {
    /* Check all inputs for legality */

    if cv_mem.cv_sensi != SUNTRUE {
        cvProcessError(Some(cv_mem), CV_NO_SENS, line!(), "CVodeGetSensDky1", file!(), MSGCV_NO_SENSI);
        return CV_NO_SENS;
    }

    if dkyS.is_empty() {
        cvProcessError(Some(cv_mem), CV_BAD_DKY, line!(), "CVodeGetSensDky1", file!(), MSGCV_NULL_DKY);
        return CV_BAD_DKY;
    }

    if k < 0 || k > cv_mem.cv_q {
        cvProcessError(Some(cv_mem), CV_BAD_K, line!(), "CVodeGetSensDky1", file!(), MSGCV_BAD_K);
        return CV_BAD_K;
    }

    if is < 0 || is > cv_mem.cv_Ns - 1 {
        cvProcessError(Some(cv_mem), CV_BAD_IS, line!(), "CVodeGetSensDky1", file!(), MSGCV_BAD_IS);
        return CV_BAD_IS;
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
        cvProcessError(Some(cv_mem), CV_BAD_T, line!(), "CVodeGetSensDky1", file!(),
            &format!("Illegal value for t. t = {} is not between tcur - hold = {} and tcur = {}",
                     t, cv_mem.cv_tn - cv_mem.cv_hu, cv_mem.cv_tn));
        return CV_BAD_T;
    }

    /* Sum the differentiated interpolating polynomial */
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
        let znSj = &cv_mem.cv_znS[j as usize][is as usize];
        if first {
            for (d, z) in dkyS.data.iter_mut().zip(&znSj.data) {
                *d = c * *z;
            }
            first = false;
        } else {
            for (d, z) in dkyS.data.iter_mut().zip(&znSj.data) {
                *d += c * *z;
            }
        }
        j -= 1;
    }

    if k == 0 {
        return CV_SUCCESS;
    }
    let r = SUNRpowerI(cv_mem.cv_h, -k);
    dkyS.scale_inplace(r);

    CV_SUCCESS
}

/*
 * CVodeGetQuadSens and CVodeGetQuadSens1
 *
 * Extraction functions for all or only one of the quadrature sensitivity
 * vectors at the time at which CVode returned the ODE solution.
 */
pub fn CVodeGetQuadSens(cv_mem: &CVodeMem, tret: &mut f64, yQSout: &mut [NVector]) -> i32 {
    *tret = cv_mem.cv_tretlast;

    CVodeGetQuadSensDky(cv_mem, cv_mem.cv_tretlast, 0, yQSout)
}

pub fn CVodeGetQuadSens1(cv_mem: &CVodeMem, tret: &mut f64, is: i32, yQSout: &mut NVector) -> i32 {
    *tret = cv_mem.cv_tretlast;

    CVodeGetQuadSensDky1(cv_mem, cv_mem.cv_tretlast, 0, is, yQSout)
}

/*
 * CVodeGetQuadSensDky and CVodeGetQuadSensDky1
 *
 * Dense output functions for all or only one of the quadrature sensitivity
 * vectors (or derivative thereof).
 */
pub fn CVodeGetQuadSensDky(cv_mem: &CVodeMem, t: f64, k: i32, dkyQS_all: &mut [NVector]) -> i32 {
    let mut ier = CV_SUCCESS;

    if dkyQS_all.is_empty() {
        cvProcessError(Some(cv_mem), CV_BAD_DKY, line!(), "CVodeGetQuadSensDky", file!(), MSGCV_NULL_DKYA);
        return CV_BAD_DKY;
    }

    for is in 0..cv_mem.cv_Ns as usize {
        ier = CVodeGetQuadSensDky1(cv_mem, t, k, is as i32, &mut dkyQS_all[is]);
        if ier != CV_SUCCESS {
            break;
        }
    }

    ier
}

pub fn CVodeGetQuadSensDky1(
    cv_mem: &CVodeMem,
    t: f64,
    k: i32,
    is: i32,
    dkyQS: &mut NVector,
) -> i32 {
    /* Check all inputs for legality */

    if cv_mem.cv_quadr_sensi != SUNTRUE {
        cvProcessError(Some(cv_mem), CV_NO_QUADSENS, line!(), "CVodeGetQuadSensDky1", file!(),
                       MSGCV_NO_QUADSENSI);
        return CV_NO_QUADSENS;
    }

    if dkyQS.is_empty() {
        cvProcessError(Some(cv_mem), CV_BAD_DKY, line!(), "CVodeGetQuadSensDky1", file!(), MSGCV_NULL_DKY);
        return CV_BAD_DKY;
    }

    if k < 0 || k > cv_mem.cv_q {
        cvProcessError(Some(cv_mem), CV_BAD_K, line!(), "CVodeGetQuadSensDky1", file!(), MSGCV_BAD_K);
        return CV_BAD_K;
    }

    if is < 0 || is > cv_mem.cv_Ns - 1 {
        cvProcessError(Some(cv_mem), CV_BAD_IS, line!(), "CVodeGetQuadSensDky1", file!(), MSGCV_BAD_IS);
        return CV_BAD_IS;
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
        cvProcessError(Some(cv_mem), CV_BAD_T, line!(), "CVodeGetQuadSensDky1", file!(),
            &format!("Illegal value for t. t = {} is not between tcur - hold = {} and tcur = {}",
                     t, cv_mem.cv_tn - cv_mem.cv_hu, cv_mem.cv_tn));
        return CV_BAD_T;
    }

    /* Sum the differentiated interpolating polynomial */
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
        let znQSj = &cv_mem.cv_znQS[j as usize][is as usize];
        if first {
            for (d, z) in dkyQS.data.iter_mut().zip(&znQSj.data) {
                *d = c * *z;
            }
            first = false;
        } else {
            for (d, z) in dkyQS.data.iter_mut().zip(&znQSj.data) {
                *d += c * *z;
            }
        }
        j -= 1;
    }

    if k == 0 {
        return CV_SUCCESS;
    }
    let r = SUNRpowerI(cv_mem.cv_h, -k);
    dkyQS.scale_inplace(r);

    CV_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Deallocation functions
 * -----------------------------------------------------------------
 */

/*
 * CVodeFree
 *
 * This routine frees the problem memory allocated by CVodeInit.
 * Such memory includes all the vectors allocated by cvAllocVectors,
 * and the memory lmem for the linear solver (deallocated by a call
 * to lfree), as well as (if Ns!=0) all memory allocated for
 * sensitivity computations by CVodeSensInit.
 *
 * (RAII: dropping the Box releases everything the C code frees
 * explicitly. The C call sequence is kept so the workspace
 * bookkeeping and flag resets run exactly as in C before the drop;
 * cv_lfree, the rootfinding array frees, the fused-op scratch frees
 * and cvProjFree collapse into the final drop.)
 */
pub fn CVodeFree(mut cvode_mem: Box<CVodeMem>) {
    let cv_mem = &mut *cvode_mem;

    cvFreeVectors(cv_mem);

    /* if CVODE created the nonlinear solver object then free it */
    if cv_mem.ownNLS {
        cv_mem.NLS = None;
        cv_mem.ownNLS = SUNFALSE;
    }

    CVodeQuadFree(cv_mem);

    CVodeSensFree(cv_mem);

    CVodeQuadSensFree(cv_mem);

    crate::cvodea::CVodeAdjFree(cv_mem);

    /* cv_lfree: dropping the linear solver module frees its memory */
    cv_mem.cv_lmem = LsModule::None;

    /* rootfinding arrays, projection memory: dropped with the Box */
}

/*
 * CVodeQuadFree
 *
 * CVodeQuadFree frees the problem memory in cvode_mem allocated
 * for quadrature integration. Its only argument is the pointer
 * cvode_mem returned by CVodeCreate.
 */
pub fn CVodeQuadFree(cv_mem: &mut CVodeMem) {
    if cv_mem.cv_QuadMallocDone {
        cvQuadFreeVectors(cv_mem);
        cv_mem.cv_QuadMallocDone = SUNFALSE;
        cv_mem.cv_quadr = SUNFALSE;
    }
}

/*
 * CVodeSensFree
 *
 * CVodeSensFree frees the problem memory in cvode_mem allocated
 * for sensitivity analysis. Its only argument is the pointer
 * cvode_mem returned by CVodeCreate.
 */
pub fn CVodeSensFree(cv_mem: &mut CVodeMem) {
    if cv_mem.cv_SensMallocDone {
        if cv_mem.cv_stgr1alloc {
            cv_mem.cv_ncfS1 = Vec::new();
            cv_mem.cv_ncfnS1 = Vec::new();
            cv_mem.cv_nniS1 = Vec::new();
            cv_mem.cv_nnfS1 = Vec::new();
            cv_mem.cv_stgr1alloc = SUNFALSE;
        }
        cvSensFreeVectors(cv_mem);
        cv_mem.cv_SensMallocDone = SUNFALSE;
        cv_mem.cv_sensi = SUNFALSE;
    }

    /* free any vector wrappers (the C zn0Sim/ycorSim/ewtSim and
       zn0Stg/ycorStg/ewtStg senswrapper aliases are not stored in
       this port -- pinned decision 3 in cvodes_impl.rs -- so only
       the allocation flags are reset) */
    if cv_mem.simMallocDone {
        cv_mem.simMallocDone = SUNFALSE;
    }
    if cv_mem.stgMallocDone {
        cv_mem.stgMallocDone = SUNFALSE;
    }

    /* if CVODES created a NLS object then free it */
    if cv_mem.ownNLSsim {
        cv_mem.NLSsim = None;
        cv_mem.ownNLSsim = SUNFALSE;
    }
    if cv_mem.ownNLSstg {
        cv_mem.NLSstg = None;
        cv_mem.ownNLSstg = SUNFALSE;
    }
    if cv_mem.ownNLSstg1 {
        cv_mem.NLSstg1 = None;
        cv_mem.ownNLSstg1 = SUNFALSE;
    }

    /* free min atol array if necessary */
    if !cv_mem.cv_atolSmin0.is_empty() {
        cv_mem.cv_atolSmin0 = Vec::new();
    }
}

/*
 * CVodeQuadSensFree
 *
 * CVodeQuadSensFree frees the problem memory in cvode_mem allocated
 * for quadrature sensitivity analysis. Its only argument is the pointer
 * cvode_mem returned by CVodeCreate.
 */
pub fn CVodeQuadSensFree(cv_mem: &mut CVodeMem) {
    if cv_mem.cv_QuadSensMallocDone {
        cvQuadSensFreeVectors(cv_mem);
        cv_mem.cv_QuadSensMallocDone = SUNFALSE;
        cv_mem.cv_quadr_sensi = SUNFALSE;
    }

    /* free min atol array if necessary */
    if !cv_mem.cv_atolQSmin0.is_empty() {
        cv_mem.cv_atolQSmin0 = Vec::new();
    }
}

/*
 * =================================================================
 *  Private Functions Implementation
 * =================================================================
 */

/*
 * cvCheckNvector (cvodes.c:4631) is dropped entirely, mirroring the
 * donor port: it only verifies that the ops table of the template
 * N_Vector provides nvclone/nvdestroy/nvlinearsum/nvconst/nvprod/
 * nvdiv/nvscale/nvabs/nvinv/nvaddconst/nvmaxnorm/nvwrmsnorm, and the
 * concrete serial NVector implements every one of these, so the
 * check is vacuously SUNTRUE.
 */

/*
 * -----------------------------------------------------------------
 * Memory allocation/deallocation
 * -----------------------------------------------------------------
 */

/*
 * cvAllocVectors
 *
 * This routine allocates the CVODES vectors ewt, acor, tempv, ftemp, and
 * zn[0], ..., zn[maxord].
 * (Infallible in Rust: Vec/Box allocation aborts on OOM, so the
 * SUNFALSE unwind ladders vanish and the return type is ().)
 * This routine also sets the optional outputs lrw and liw, which are
 * (respectively) the lengths of the real and integer work spaces
 * allocated here.
 */
pub(crate) fn cvAllocVectors(cv_mem: &mut CVodeMem, tmpl: &NVector) {
    /* Allocate ewt, acor, tempv, ftemp */
    cv_mem.cv_ewt = N_VClone(tmpl);
    cv_mem.cv_acor = N_VClone(tmpl);
    cv_mem.cv_tempv = N_VClone(tmpl);
    cv_mem.cv_ftemp = N_VClone(tmpl);
    cv_mem.cv_vtemp1 = N_VClone(tmpl);
    cv_mem.cv_vtemp2 = N_VClone(tmpl);
    cv_mem.cv_vtemp3 = N_VClone(tmpl);
    /* cv_y is owned here (in C it aliases the user's yout and is not
       allocated); it is not counted in lrw/liw, as in C */
    cv_mem.cv_y = N_VClone(tmpl);

    /* Allocate zn[0] ... zn[qmax] (the full L_MAX array is allocated,
       donor precedent: in C the zn array has L_MAX slots of which
       qmax+1 are populated; lrw/liw count only qmax+8 as in C) */
    cv_mem.cv_zn = (0..L_MAX).map(|_| N_VClone(tmpl)).collect();

    /* Update solver workspace lengths  */
    cv_mem.cv_lrw += (cv_mem.cv_qmax as i64 + 8) * cv_mem.cv_lrw1;
    cv_mem.cv_liw += (cv_mem.cv_qmax as i64 + 8) * cv_mem.cv_liw1;

    /* Store the value of qmax used here */
    cv_mem.cv_qmax_alloc = cv_mem.cv_qmax;
}

/*
 * cvFreeVectors
 *
 * This routine frees the vectors allocated in cvAllocVectors.
 */
fn cvFreeVectors(cv_mem: &mut CVodeMem) {
    let maxord = cv_mem.cv_qmax_alloc;

    cv_mem.cv_ewt = NVector::default();
    cv_mem.cv_acor = NVector::default();
    cv_mem.cv_tempv = NVector::default();
    cv_mem.cv_ftemp = NVector::default();
    cv_mem.cv_vtemp1 = NVector::default();
    cv_mem.cv_vtemp2 = NVector::default();
    cv_mem.cv_vtemp3 = NVector::default();
    cv_mem.cv_zn = Vec::new();

    cv_mem.cv_lrw -= (maxord as i64 + 8) * cv_mem.cv_lrw1;
    cv_mem.cv_liw -= (maxord as i64 + 8) * cv_mem.cv_liw1;

    if cv_mem.cv_VabstolMallocDone {
        cv_mem.cv_Vabstol = NVector::default();
        cv_mem.cv_lrw -= cv_mem.cv_lrw1;
        cv_mem.cv_liw -= cv_mem.cv_liw1;
    }

    if !cv_mem.cv_constraints.is_empty() {
        cv_mem.cv_constraints = NVector::default();
        cv_mem.cv_lrw -= cv_mem.cv_lrw1;
        cv_mem.cv_liw -= cv_mem.cv_liw1;
    }
}

/*
 * CVodeQuadAllocVectors
 *
 * NOTE: Space for ewtQ is allocated even when errconQ=SUNFALSE,
 * although in this case, ewtQ is never used. The reason for this
 * decision is to allow the user to re-initialize the quadrature
 * computation with errconQ=SUNTRUE, after an initialization with
 * errconQ=SUNFALSE, without new memory allocation within
 * CVodeQuadReInit.
 */
pub(crate) fn cvQuadAllocVectors(cv_mem: &mut CVodeMem, tmpl: &NVector) {
    /* Allocate ewtQ */
    cv_mem.cv_ewtQ = N_VClone(tmpl);

    /* Allocate acorQ */
    cv_mem.cv_acorQ = N_VClone(tmpl);

    /* Allocate yQ */
    cv_mem.cv_yQ = N_VClone(tmpl);

    /* Allocate tempvQ */
    cv_mem.cv_tempvQ = N_VClone(tmpl);

    /* Allocate zQn[0] ... zQn[maxord] */
    cv_mem.cv_znQ = (0..=cv_mem.cv_qmax as usize).map(|_| N_VClone(tmpl)).collect();

    /* Store the value of qmax used here */
    cv_mem.cv_qmax_allocQ = cv_mem.cv_qmax;

    /* Update solver workspace lengths */
    cv_mem.cv_lrw += (cv_mem.cv_qmax as i64 + 5) * cv_mem.cv_lrw1Q;
    cv_mem.cv_liw += (cv_mem.cv_qmax as i64 + 5) * cv_mem.cv_liw1Q;
}

/*
 * cvQuadFreeVectors
 *
 * This routine frees the CVODES vectors allocated in cvQuadAllocVectors.
 */
fn cvQuadFreeVectors(cv_mem: &mut CVodeMem) {
    let maxord = cv_mem.cv_qmax_allocQ;

    cv_mem.cv_ewtQ = NVector::default();
    cv_mem.cv_acorQ = NVector::default();
    cv_mem.cv_yQ = NVector::default();
    cv_mem.cv_tempvQ = NVector::default();

    cv_mem.cv_znQ = Vec::new();

    cv_mem.cv_lrw -= (maxord as i64 + 5) * cv_mem.cv_lrw1Q;
    cv_mem.cv_liw -= (maxord as i64 + 5) * cv_mem.cv_liw1Q;

    if cv_mem.cv_VabstolQMallocDone {
        cv_mem.cv_VabstolQ = NVector::default();
        cv_mem.cv_lrw -= cv_mem.cv_lrw1Q;
        cv_mem.cv_liw -= cv_mem.cv_liw1Q;
    }

    cv_mem.cv_VabstolQMallocDone = SUNFALSE;
}

/*
 * cvSensAllocVectors
 *
 * Create (through duplication) N_Vectors used for sensitivity analysis,
 * using the N_Vector 'tmpl' as a template.
 */
pub(crate) fn cvSensAllocVectors(cv_mem: &mut CVodeMem, tmpl: &NVector) {
    let ns = cv_mem.cv_Ns as usize;

    /* Allocate yS */
    cv_mem.cv_yS = (0..ns).map(|_| N_VClone(tmpl)).collect();

    /* Allocate ewtS */
    cv_mem.cv_ewtS = (0..ns).map(|_| N_VClone(tmpl)).collect();

    /* Allocate acorS */
    cv_mem.cv_acorS = (0..ns).map(|_| N_VClone(tmpl)).collect();

    /* Allocate tempvS */
    cv_mem.cv_tempvS = (0..ns).map(|_| N_VClone(tmpl)).collect();

    /* Allocate ftempS */
    cv_mem.cv_ftempS = (0..ns).map(|_| N_VClone(tmpl)).collect();

    /* Allocate znS */
    for j in 0..=cv_mem.cv_qmax as usize {
        cv_mem.cv_znS[j] = (0..ns).map(|_| N_VClone(tmpl)).collect();
    }

    /* Allocate space for pbar and plist (C mallocs them here without
       initialization; CVodeSensInit/CVodeSensInit1 fill them in) */
    cv_mem.cv_pbar = vec![ZERO; ns];
    cv_mem.cv_plist = vec![0; ns];

    /* Update solver workspace lengths */
    cv_mem.cv_lrw +=
        (cv_mem.cv_qmax as i64 + 6) * cv_mem.cv_Ns as i64 * cv_mem.cv_lrw1 + cv_mem.cv_Ns as i64;
    cv_mem.cv_liw +=
        (cv_mem.cv_qmax as i64 + 6) * cv_mem.cv_Ns as i64 * cv_mem.cv_liw1 + cv_mem.cv_Ns as i64;

    /* Store the value of qmax used here */
    cv_mem.cv_qmax_allocS = cv_mem.cv_qmax;
}

/*
 * cvSensFreeVectors
 *
 * This routine frees the CVODES vectors allocated in cvSensAllocVectors.
 */
fn cvSensFreeVectors(cv_mem: &mut CVodeMem) {
    let maxord = cv_mem.cv_qmax_allocS;

    cv_mem.cv_yS = Vec::new();
    cv_mem.cv_ewtS = Vec::new();
    cv_mem.cv_acorS = Vec::new();
    cv_mem.cv_tempvS = Vec::new();
    cv_mem.cv_ftempS = Vec::new();

    for j in 0..=maxord as usize {
        cv_mem.cv_znS[j] = Vec::new();
    }

    cv_mem.cv_pbar = Vec::new();
    cv_mem.cv_plist = Vec::new();

    cv_mem.cv_lrw -=
        (maxord as i64 + 6) * cv_mem.cv_Ns as i64 * cv_mem.cv_lrw1 + cv_mem.cv_Ns as i64;
    cv_mem.cv_liw -=
        (maxord as i64 + 6) * cv_mem.cv_Ns as i64 * cv_mem.cv_liw1 + cv_mem.cv_Ns as i64;

    if cv_mem.cv_VabstolSMallocDone {
        cv_mem.cv_VabstolS = Vec::new();
        cv_mem.cv_lrw -= cv_mem.cv_Ns as i64 * cv_mem.cv_lrw1;
        cv_mem.cv_liw -= cv_mem.cv_Ns as i64 * cv_mem.cv_liw1;
    }
    if cv_mem.cv_SabstolSMallocDone {
        cv_mem.cv_SabstolS = Vec::new();
        cv_mem.cv_lrw -= cv_mem.cv_Ns as i64;
    }
    cv_mem.cv_VabstolSMallocDone = SUNFALSE;
    cv_mem.cv_SabstolSMallocDone = SUNFALSE;
}

/*
 * cvQuadSensAllocVectors
 *
 * Create (through duplication) N_Vectors used for quadrature sensitivity analysis,
 * using the N_Vector 'tmpl' as a template.
 */
pub(crate) fn cvQuadSensAllocVectors(cv_mem: &mut CVodeMem, tmpl: &NVector) {
    let ns = cv_mem.cv_Ns as usize;

    /* Allocate ftempQ */
    cv_mem.cv_ftempQ = N_VClone(tmpl);

    /* Allocate yQS */
    cv_mem.cv_yQS = (0..ns).map(|_| N_VClone(tmpl)).collect();

    /* Allocate ewtQS */
    cv_mem.cv_ewtQS = (0..ns).map(|_| N_VClone(tmpl)).collect();

    /* Allocate acorQS */
    cv_mem.cv_acorQS = (0..ns).map(|_| N_VClone(tmpl)).collect();

    /* Allocate tempvQS */
    cv_mem.cv_tempvQS = (0..ns).map(|_| N_VClone(tmpl)).collect();

    /* Allocate znQS */
    for j in 0..=cv_mem.cv_qmax as usize {
        cv_mem.cv_znQS[j] = (0..ns).map(|_| N_VClone(tmpl)).collect();
    }

    /* Update solver workspace lengths */
    cv_mem.cv_lrw += (cv_mem.cv_qmax as i64 + 5) * cv_mem.cv_Ns as i64 * cv_mem.cv_lrw1Q;
    cv_mem.cv_liw += (cv_mem.cv_qmax as i64 + 5) * cv_mem.cv_Ns as i64 * cv_mem.cv_liw1Q;

    /* Store the value of qmax used here */
    cv_mem.cv_qmax_allocQS = cv_mem.cv_qmax;
}

/*
 * cvQuadSensFreeVectors
 *
 * This routine frees the CVODES vectors allocated in cvQuadSensAllocVectors.
 */
fn cvQuadSensFreeVectors(cv_mem: &mut CVodeMem) {
    let maxord = cv_mem.cv_qmax_allocQS;

    cv_mem.cv_ftempQ = NVector::default();

    cv_mem.cv_yQS = Vec::new();
    cv_mem.cv_ewtQS = Vec::new();
    cv_mem.cv_acorQS = Vec::new();
    cv_mem.cv_tempvQS = Vec::new();

    for j in 0..=maxord as usize {
        cv_mem.cv_znQS[j] = Vec::new();
    }

    cv_mem.cv_lrw -= (maxord as i64 + 5) * cv_mem.cv_Ns as i64 * cv_mem.cv_lrw1Q;
    cv_mem.cv_liw -= (maxord as i64 + 5) * cv_mem.cv_Ns as i64 * cv_mem.cv_liw1Q;

    if cv_mem.cv_VabstolQSMallocDone {
        cv_mem.cv_VabstolQS = Vec::new();
        cv_mem.cv_lrw -= cv_mem.cv_Ns as i64 * cv_mem.cv_lrw1Q;
        cv_mem.cv_liw -= cv_mem.cv_Ns as i64 * cv_mem.cv_liw1Q;
    }
    if cv_mem.cv_SabstolQSMallocDone {
        cv_mem.cv_SabstolQS = Vec::new();
        cv_mem.cv_lrw -= cv_mem.cv_Ns as i64;
    }
    cv_mem.cv_VabstolQSMallocDone = SUNFALSE;
    cv_mem.cv_SabstolQSMallocDone = SUNFALSE;
}

/*
 * -----------------------------------------------------------------
 * Initial setup
 * -----------------------------------------------------------------
 */

/*
 * cvInitialSetup
 *
 * This routine performs input consistency checks at the first step.
 * If needed, it also checks the linear solver module and calls the
 * linear solver initialization routine.
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

    /* (The C check that N_VMin is available when abstol==0 is dropped:
       the concrete serial NVector always provides it. The cv_e_data
       assignment is omitted -- pinned decision 5 in cvodes_impl.rs.) */

    /* Check to see if y0 satisfies constraints */
    if !cv_mem.cv_constraints.is_empty() {
        if cv_mem.cv_sensi && cv_mem.cv_ism == CV_SIMULTANEOUS {
            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "cvInitialSetup", file!(),
                           MSGCV_BAD_ISM_CONSTR);
            return CV_ILL_INPUT;
        }

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

    /* Quadrature initial setup */

    if cv_mem.cv_quadr && cv_mem.cv_errconQ {
        /* Did the user specify tolerances? */
        if cv_mem.cv_itolQ == CV_NN {
            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "cvInitialSetup", file!(), MSGCV_NO_TOLQ);
            return CV_ILL_INPUT;
        }

        /* Load ewtQ */
        let ier = cvQuadEwtSet_apply_to_ewtQ(cv_mem);
        if ier != 0 {
            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "cvInitialSetup", file!(), MSGCV_BAD_EWTQ);
            return CV_ILL_INPUT;
        }
    }

    if !cv_mem.cv_quadr {
        cv_mem.cv_errconQ = SUNFALSE;
    }

    /* Forward sensitivity initial setup */

    if cv_mem.cv_sensi {
        /* Did the user specify tolerances? */
        if cv_mem.cv_itolS == CV_NN {
            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "cvInitialSetup", file!(), MSGCV_NO_TOLS);
            return CV_ILL_INPUT;
        }

        /* If using the internal DQ functions, we must have access to the problem parameters */
        if cv_mem.cv_fSDQ && cv_mem.cv_p.is_empty() {
            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "cvInitialSetup", file!(), MSGCV_NULL_P);
            return CV_ILL_INPUT;
        }

        /* (Rust-port guard, no C counterpart: the internal DQ perturbs
           p[which] through the user's own array in C.  Here that is
           only possible via the FSAUserData convention — without it
           the user RHS would silently see frozen parameters and every
           sensitivity would come out zero.) */
        if cv_mem.cv_fSDQ {
            let is_fsa = cv_mem
                .cv_user_data
                .as_ref()
                .map(|d| d.is::<FSAUserData>())
                .unwrap_or(false);
            if !is_fsa {
                cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "cvInitialSetup", file!(),
                               "Internal DQ sensitivities require the FSAUserData user-data \
                                convention (see sundials_types.rs).");
                return CV_ILL_INPUT;
            }
        }

        /* Load ewtS */
        let ier = cvSensEwtSet_apply_to_ewtS(cv_mem);
        if ier != 0 {
            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "cvInitialSetup", file!(), MSGCV_BAD_EWTS);
            return CV_ILL_INPUT;
        }
    }

    /* FSA of quadrature variables */

    if cv_mem.cv_quadr_sensi {
        /* If using the internal DQ functions, we must have access to fQ
         * (i.e. quadrature integration must be enabled) and to the problem parameters */

        if cv_mem.cv_fQSDQ {
            /* Test if quadratures are defined, so we can use fQ */
            if !cv_mem.cv_quadr {
                cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "cvInitialSetup", file!(), MSGCV_NULL_FQ);
                return CV_ILL_INPUT;
            }

            /* Test if we have the problem parameters */
            if cv_mem.cv_p.is_empty() {
                cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "cvInitialSetup", file!(), MSGCV_NULL_P);
                return CV_ILL_INPUT;
            }

            /* (Rust-port guard, no C counterpart — see the sensitivity
               FSAUserData guard above; the quad-sens DQ perturbs
               p[which] the same way.) */
            let is_fsa = cv_mem
                .cv_user_data
                .as_ref()
                .map(|d| d.is::<FSAUserData>())
                .unwrap_or(false);
            if !is_fsa {
                cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "cvInitialSetup", file!(),
                               "Internal DQ quadrature sensitivities require the FSAUserData \
                                user-data convention (see sundials_types.rs).");
                return CV_ILL_INPUT;
            }
        }

        if cv_mem.cv_errconQS {
            /* Did the user specify tolerances? */
            if cv_mem.cv_itolQS == CV_NN {
                cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "cvInitialSetup", file!(), MSGCV_NO_TOLQS);
                return CV_ILL_INPUT;
            }

            /* If needed, did the user provide quadrature tolerances? */
            if cv_mem.cv_itolQS == CV_EE && cv_mem.cv_itolQ == CV_NN {
                cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "cvInitialSetup", file!(), MSGCV_NO_TOLQ);
                return CV_ILL_INPUT;
            }

            /* Load ewtQS */
            let ier = cvQuadSensEwtSet_apply_to_ewtQS(cv_mem);
            if ier != 0 {
                cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "cvInitialSetup", file!(), MSGCV_BAD_EWTQS);
                return CV_ILL_INPUT;
            }
        }
    } else {
        cv_mem.cv_errconQS = SUNFALSE;
    }

    /* Call linit function (if it exists) */
    if cv_has_lsetup(cv_mem) {
        let ier = cv_linit_dispatch(cv_mem);
        if ier != 0 {
            cvProcessError(Some(cv_mem), CV_LINIT_FAIL, line!(), "cvInitialSetup", file!(), MSGCV_LINIT_FAIL);
            return CV_LINIT_FAIL;
        }
    }

    /* Initialize the nonlinear solver (must occur after linear solver is
       initialized) so that lsetup and lsolve pointer have been set */

    /* always initialize the ODE NLS in case the user disables sensitivities */
    let ier = crate::cvodes_nls::cvNlsInit(cv_mem);
    if ier != 0 {
        cvProcessError(Some(cv_mem), CV_NLS_INIT_FAIL, line!(), "cvInitialSetup", file!(), MSGCV_NLS_INIT_FAIL);
        return CV_NLS_INIT_FAIL;
    }

    if cv_mem.NLSsim.is_some() {
        let ier = crate::cvodes_nls_sim::cvNlsInitSensSim(cv_mem);
        if ier != 0 {
            cvProcessError(Some(cv_mem), CV_NLS_INIT_FAIL, line!(), "cvInitialSetup", file!(), MSGCV_NLS_INIT_FAIL);
            return CV_NLS_INIT_FAIL;
        }
    }

    if cv_mem.NLSstg.is_some() {
        let ier = crate::cvodes_nls_stg::cvNlsInitSensStg(cv_mem);
        if ier != 0 {
            cvProcessError(Some(cv_mem), CV_NLS_INIT_FAIL, line!(), "cvInitialSetup", file!(), MSGCV_NLS_INIT_FAIL);
            return CV_NLS_INIT_FAIL;
        }
    }

    if cv_mem.NLSstg1.is_some() {
        let ier = crate::cvodes_nls_stg1::cvNlsInitSensStg1(cv_mem);
        if ier != 0 {
            cvProcessError(Some(cv_mem), CV_NLS_INIT_FAIL, line!(), "cvInitialSetup", file!(), MSGCV_NLS_INIT_FAIL);
            return CV_NLS_INIT_FAIL;
        }
    }

    /* Initialize projection data */
    if cv_mem.proj_enabled && cv_mem.proj_mem.is_none() {
        cvProcessError(Some(cv_mem), CV_PROJ_MEM_NULL, line!(), "cvInitialSetup", file!(),
                       MSG_CV_PROJ_MEM_NULL);
        return CV_PROJ_MEM_NULL;
    }

    if let Some(pm) = cv_mem.proj_mem.as_deref_mut() {
        let ier = crate::cvodes_proj::cvProjInit(pm);
        if ier != CV_SUCCESS {
            cvProcessError(Some(cv_mem), CV_MEM_FAIL, line!(), "cvInitialSetup", file!(), MSGCV_MEM_FAIL);
            return CV_MEM_FAIL;
        }
        cv_mem.proj_applied = SUNFALSE;
    }

    /* Initial setup complete */
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
 * This routine computes a tentative initial step size h0. Note that here tout
 * is either the value passed to CVode at the first call or the value of tstop
 * (if tstop is enabled and it is closer to t0=tn than tout). If any RHS
 * function fails unrecoverably, cvHin returns CV_*RHSFUNC_FAIL. If any RHS
 * function fails recoverably too many times and recovery is not possible, cvHin
 * returns CV_REPTD_*RHSFUNC_ERR. Otherwise, cvHin sets h to the chosen value
 * h0 and returns CV_SUCCESS.
 *
 * The algorithm used seeks to find h0 as a solution of
 *       (WRMS norm of (h0^2 ydd / 2)) = 1,
 * where ydd = estimated second derivative of y. Here, y includes
 * all variables considered in the error test.
 *
 * We start with an initial estimate equal to the geometric mean of the
 * lower and upper bounds on the step size.
 *
 * Loop up to MAX_ITERS times to find h0.
 * Stop if new and previous values differ by a factor < 2.
 * Stop if hnew/hg > 2 after one iteration, as this probably means
 * that the ydd value is bad because of cancellation error.
 *
 * For each new proposed hg, we allow MAX_ITERS attempts to
 * resolve a possible recoverable failure from f() by reducing
 * the proposed stepsize by a factor of 0.2. If a legal stepsize
 * still cannot be found, fall back on a previous value if possible,
 * or else return CV_REPTD_RHSFUNC_ERR.
 *
 * Finally, we apply a bias (0.5) and verify that h0 is within bounds.
 */
fn cvHin(cv_mem: &mut CVodeMem, tout: f64) -> i32 {
    /* cvInitialSetup checks for tdiff = 0 or < 2 * troundoff */
    let tdiff = tout - cv_mem.cv_tn;
    let sign = if tdiff > ZERO { 1 } else { -1 };
    let tdist = SUNRabs(tdiff);
    let tround = cv_mem.cv_uround * SUNMAX(SUNRabs(cv_mem.cv_tn), SUNRabs(tout));

    /*
       Set lower and upper bounds on h0, and take geometric mean
       as first trial value.
       Exit with this value if the bounds cross each other.
    */

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
        let mut retval = 0;

        for _count2 in 1..=MAX_ITERS {
            let hgs = hg * sign as f64;
            retval = cvYddNorm(cv_mem, hgs, &mut yddnrm);
            /* If a RHS function failed unrecoverably, give up */
            if retval < 0 {
                return retval;
            }
            /* If successful, we can use ydd */
            if retval == CV_SUCCESS {
                hg_ok = SUNTRUE;
                break;
            }
            /* A RHS function failed recoverably; cut step size and test again */
            hg *= POINT2;
        }

        /* If a RHS function failed recoverably MAX_ITERS times */

        if !hg_ok {
            /* Exit if this is the first or second pass. No recovery possible */
            if count1 <= 2 {
                if retval == RHSFUNC_RECVR {
                    return CV_REPTD_RHSFUNC_ERR;
                }
                if retval == QRHSFUNC_RECVR {
                    return CV_REPTD_QRHSFUNC_ERR;
                }
                if retval == SRHSFUNC_RECVR {
                    return CV_REPTD_SRHSFUNC_ERR;
                }
            }
            /* We have a fall-back option. The value hs is a previous hnew which
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
 * This routine sets an upper bound on abs(h0) based on
 * tdist = tn - t0 and the values of y[i]/y'[i].
 */
fn cvUpperBoundH0(cv_mem: &mut CVodeMem, tdist: f64) -> f64 {
    /*
     * Bound based on |y|/|y'| -- allow at most an increase of
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
    let mut hub_inv = N_VMaxNorm(&cv_mem.cv_tempv);

    /* Bound based on |yQ|/|yQ'| */

    if cv_mem.cv_quadr && cv_mem.cv_errconQ {
        /* tempQ1 = cv_tempvQ, tempQ2 = cv_acorQ in the C code */

        {
            let (znQ0, acorQ) = (&cv_mem.cv_znQ[0], &mut cv_mem.cv_acorQ);
            N_VAbs(znQ0, acorQ);
        }
        {
            let mut wQ = std::mem::take(&mut cv_mem.cv_tempvQ);
            let _ = cvQuadEwtSet(cv_mem, &cv_mem.cv_znQ[0], &mut wQ);
            cv_mem.cv_tempvQ = wQ;
        }
        cv_mem.cv_tempvQ.invert_inplace();
        {
            let (acorQ, tempvQ) = (&cv_mem.cv_acorQ, &mut cv_mem.cv_tempvQ);
            tempvQ.linear_sum_with(ONE, HUB_FACTOR, acorQ);
        }

        {
            let (znQ1, acorQ) = (&cv_mem.cv_znQ[1], &mut cv_mem.cv_acorQ);
            N_VAbs(znQ1, acorQ);
        }

        {
            let (acorQ, tempvQ) = (&cv_mem.cv_acorQ, &mut cv_mem.cv_tempvQ);
            for (t1, t2) in tempvQ.data.iter_mut().zip(&acorQ.data) {
                *t1 = *t2 / *t1;
            }
        }
        let hubQ_inv = N_VMaxNorm(&cv_mem.cv_tempvQ);

        if hubQ_inv > hub_inv {
            hub_inv = hubQ_inv;
        }
    }

    /* Bound based on |yS|/|yS'| */

    if cv_mem.cv_sensi && cv_mem.cv_errconS {
        /* tempS1 = cv_acorS in the C code */
        {
            let yScur = std::mem::take(&mut cv_mem.cv_znS[0]);
            let mut tempS1 = std::mem::take(&mut cv_mem.cv_acorS);
            let _ = cvSensEwtSet(cv_mem, &yScur, &mut tempS1);
            cv_mem.cv_znS[0] = yScur;
            cv_mem.cv_acorS = tempS1;
        }

        for is in 0..cv_mem.cv_Ns as usize {
            {
                let (znS0is, acor) = (&cv_mem.cv_znS[0][is], &mut cv_mem.cv_acor);
                N_VAbs(znS0is, acor);
            }
            {
                let (acorSis, tempv) = (&cv_mem.cv_acorS[is], &mut cv_mem.cv_tempv);
                N_VInv(acorSis, tempv);
            }
            {
                let (acor, tempv) = (&cv_mem.cv_acor, &mut cv_mem.cv_tempv);
                tempv.linear_sum_with(ONE, HUB_FACTOR, acor);
            }

            {
                let (znS1is, acor) = (&cv_mem.cv_znS[1][is], &mut cv_mem.cv_acor);
                N_VAbs(znS1is, acor);
            }

            {
                let (acor, tempv) = (&cv_mem.cv_acor, &mut cv_mem.cv_tempv);
                for (t1, t2) in tempv.data.iter_mut().zip(&acor.data) {
                    *t1 = *t2 / *t1;
                }
            }
            let hubS_inv = N_VMaxNorm(&cv_mem.cv_tempv);

            if hubS_inv > hub_inv {
                hub_inv = hubS_inv;
            }
        }
    }

    /* Bound based on |yQS|/|yQS'| */

    if cv_mem.cv_quadr_sensi && cv_mem.cv_errconQS {
        /* tempQ1 = cv_tempvQ, tempQ2 = cv_acorQ, tempQS1 = cv_acorQS */
        {
            let yQScur = std::mem::take(&mut cv_mem.cv_znQS[0]);
            let mut tempQS1 = std::mem::take(&mut cv_mem.cv_acorQS);
            let _ = cvQuadSensEwtSet(cv_mem, &yQScur, &mut tempQS1);
            cv_mem.cv_znQS[0] = yQScur;
            cv_mem.cv_acorQS = tempQS1;
        }

        for is in 0..cv_mem.cv_Ns as usize {
            {
                let (znQS0is, acorQ) = (&cv_mem.cv_znQS[0][is], &mut cv_mem.cv_acorQ);
                N_VAbs(znQS0is, acorQ);
            }
            {
                let (acorQSis, tempvQ) = (&cv_mem.cv_acorQS[is], &mut cv_mem.cv_tempvQ);
                N_VInv(acorQSis, tempvQ);
            }
            {
                let (acorQ, tempvQ) = (&cv_mem.cv_acorQ, &mut cv_mem.cv_tempvQ);
                tempvQ.linear_sum_with(ONE, HUB_FACTOR, acorQ);
            }

            {
                let (znQS1is, acorQ) = (&cv_mem.cv_znQS[1][is], &mut cv_mem.cv_acorQ);
                N_VAbs(znQS1is, acorQ);
            }

            {
                let (acorQ, tempvQ) = (&cv_mem.cv_acorQ, &mut cv_mem.cv_tempvQ);
                for (t1, t2) in tempvQ.data.iter_mut().zip(&acorQ.data) {
                    *t1 = *t2 / *t1;
                }
            }
            let hubQS_inv = N_VMaxNorm(&cv_mem.cv_tempvQ);

            if hubQS_inv > hub_inv {
                hub_inv = hubQS_inv;
            }
        }
    }

    /*
     * bound based on tdist -- allow at most a step of magnitude
     * HUB_FACTOR * tdist
     */

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
 * This routine computes an estimate of the second derivative of Y
 * using a difference quotient, and returns its WRMS norm.
 *
 * Y contains all variables included in the error test.
 */
fn cvYddNorm(cv_mem: &mut CVodeMem, hg: f64, yddnrm: &mut f64) -> i32 {
    /* y <- h*y'(t) + y(t) */

    {
        let (zn, y) = (&cv_mem.cv_zn, &mut cv_mem.cv_y);
        N_VLinearSum(hg, &zn[1], ONE, &zn[0], y);
    }

    if cv_mem.cv_sensi && cv_mem.cv_errconS {
        let ns = cv_mem.cv_Ns as usize;
        let CVodeMem { cv_znS, cv_yS, .. } = cv_mem;
        for is in 0..ns {
            N_VLinearSum(hg, &cv_znS[1][is], ONE, &cv_znS[0][is], &mut cv_yS[is]);
        }
    }

    /* tempv <- f(t+h, h*y'(t)+y(t)) */

    let f = cv_mem.cv_f.unwrap();
    let t = cv_mem.cv_tn + hg;
    let retval = {
        let CVodeMem { cv_y, cv_tempv, cv_user_data, .. } = cv_mem;
        f(t, cv_y, cv_tempv, cv_user_data)
    };
    cv_mem.cv_nfe += 1;
    if retval < 0 {
        return CV_RHSFUNC_FAIL;
    }
    if retval > 0 {
        return RHSFUNC_RECVR;
    }

    if cv_mem.cv_quadr && cv_mem.cv_errconQ {
        let fQ = cv_mem.cv_fQ.unwrap();
        let retval = {
            let CVodeMem { cv_y, cv_tempvQ, cv_user_data, .. } = cv_mem;
            fQ(t, cv_y, cv_tempvQ, cv_user_data)
        };
        cv_mem.cv_nfQe += 1;
        if retval < 0 {
            return CV_QRHSFUNC_FAIL;
        }
        if retval > 0 {
            return QRHSFUNC_RECVR;
        }
    }

    if cv_mem.cv_sensi && cv_mem.cv_errconS {
        /* wrk1 = cv_ftemp, wrk2 = cv_acor in the C code */
        let y = std::mem::take(&mut cv_mem.cv_y);
        let tempv = std::mem::take(&mut cv_mem.cv_tempv);
        let yS = std::mem::take(&mut cv_mem.cv_yS);
        let mut tempvS = std::mem::take(&mut cv_mem.cv_tempvS);
        let mut wrk1 = std::mem::take(&mut cv_mem.cv_ftemp);
        let mut wrk2 = std::mem::take(&mut cv_mem.cv_acor);
        let retval =
            cvSensRhsWrapper(cv_mem, t, &y, &tempv, &yS, &mut tempvS, &mut wrk1, &mut wrk2);
        cv_mem.cv_y = y;
        cv_mem.cv_tempv = tempv;
        cv_mem.cv_yS = yS;
        cv_mem.cv_tempvS = tempvS;
        cv_mem.cv_ftemp = wrk1;
        cv_mem.cv_acor = wrk2;
        if retval < 0 {
            return CV_SRHSFUNC_FAIL;
        }
        if retval > 0 {
            return SRHSFUNC_RECVR;
        }
    }

    if cv_mem.cv_quadr_sensi && cv_mem.cv_errconQS {
        /* wrk1 = cv_ftemp, wrk2 = cv_acorQ in the C code */
        let y = std::mem::take(&mut cv_mem.cv_y);
        let yS = std::mem::take(&mut cv_mem.cv_yS);
        let tempvQ = std::mem::take(&mut cv_mem.cv_tempvQ);
        let mut tempvQS = std::mem::take(&mut cv_mem.cv_tempvQS);
        let mut wrk1 = std::mem::take(&mut cv_mem.cv_ftemp);
        let mut wrk2 = std::mem::take(&mut cv_mem.cv_acorQ);
        let retval =
            cv_fQS_dispatch(cv_mem, t, &y, &yS, &tempvQ, &mut tempvQS, &mut wrk1, &mut wrk2);
        cv_mem.cv_y = y;
        cv_mem.cv_yS = yS;
        cv_mem.cv_tempvQ = tempvQ;
        cv_mem.cv_tempvQS = tempvQS;
        cv_mem.cv_ftemp = wrk1;
        cv_mem.cv_acorQ = wrk2;

        cv_mem.cv_nfQSe += 1;
        if retval < 0 {
            return CV_QSRHSFUNC_FAIL;
        }
        if retval > 0 {
            return QSRHSFUNC_RECVR;
        }
    }

    /* Load estimate of ||y''|| into tempv:
     * tempv <-  (1/h) * f(t+h, h*y'(t)+y(t)) - y'(t) */

    {
        let (zn, tempv) = (&cv_mem.cv_zn, &mut cv_mem.cv_tempv);
        tempv.linear_sum_with(ONE / hg, -ONE / hg, &zn[1]);
    }

    *yddnrm = N_VWrmsNorm(&cv_mem.cv_tempv, &cv_mem.cv_ewt);

    if cv_mem.cv_quadr && cv_mem.cv_errconQ {
        {
            let (znQ, tempvQ) = (&cv_mem.cv_znQ, &mut cv_mem.cv_tempvQ);
            tempvQ.linear_sum_with(ONE / hg, -ONE / hg, &znQ[1]);
        }

        *yddnrm = cvQuadUpdateNorm(cv_mem, *yddnrm, &cv_mem.cv_tempvQ, &cv_mem.cv_ewtQ);
    }

    if cv_mem.cv_sensi && cv_mem.cv_errconS {
        {
            let ns = cv_mem.cv_Ns as usize;
            let CVodeMem { cv_znS, cv_tempvS, .. } = cv_mem;
            for is in 0..ns {
                cv_tempvS[is].linear_sum_with(ONE / hg, -ONE / hg, &cv_znS[1][is]);
            }
        }

        *yddnrm = cvSensUpdateNorm(cv_mem, *yddnrm, &cv_mem.cv_tempvS, &cv_mem.cv_ewtS);
    }

    if cv_mem.cv_quadr_sensi && cv_mem.cv_errconQS {
        {
            let ns = cv_mem.cv_Ns as usize;
            let CVodeMem { cv_znQS, cv_tempvQS, .. } = cv_mem;
            for is in 0..ns {
                cv_tempvQS[is].linear_sum_with(ONE / hg, -ONE / hg, &cv_znQS[1][is]);
            }
        }

        *yddnrm = cvQuadSensUpdateNorm(cv_mem, *yddnrm, &cv_mem.cv_tempvQS, &cv_mem.cv_ewtQS);
    }

    CV_SUCCESS
}

// ===================== END PART 2 (cvodes.c:2913-6244) =====================
// ===================== PART 3 (cvodes.c:5874-10126) ========================

/*
 * cvStep
 *
 * This routine performs one internal cvode step, from tn to tn + h.
 * It calls other routines to do all the work.
 *
 * The main operations done here are as follows:
 * - preliminary adjustments if a new step size was chosen;
 * - prediction of the Nordsieck history array zn at tn + h;
 * - setting of multistep method coefficients and test quantities;
 * - solution of the nonlinear system;
 * - testing the local error;
 * - updating zn and other state data if successful;
 * - resetting stepsize and order for the next step.
 * - if SLDET is on, check for stability, reduce order if necessary.
 * On a failure in the nonlinear system solution or error test, the
 * step may be reattempted, depending on the nature of the failure.
 */
fn cvStep(cv_mem: &mut CVodeMem) -> i32 {
    /* Are we computing sensitivities with a staggered approach? */
    let do_sensi_stg = cv_mem.cv_sensi && cv_mem.cv_ism == CV_STAGGERED;
    let do_sensi_stg1 = cv_mem.cv_sensi && cv_mem.cv_ism == CV_STAGGERED1;

    /* Initialize failure counters for this step attempt */
    let mut ncf = 0; /* corrector failures  */
    let mut npf = 0; /* projection failures */
    let mut nef = 0; /* error test failures */
    let mut step_constraint_fails = 0;

    let mut ncfS = 0; /* sensitivity corrector failures           */
    let mut nefS = 0; /* sensitivity error test fails             */
    let mut nefQ = 0; /* quadrature error test fails              */
    let mut nefQS = 0; /* quadrature sensitivity error test fails  */

    if do_sensi_stg1 {
        for is in 0..cv_mem.cv_Ns as usize {
            cv_mem.cv_ncfS1[is] = 0;
        }
    }

    /* If the step size has changed, update the history array */
    if cv_mem.cv_nst > 0 && cv_mem.cv_hprime != cv_mem.cv_h {
        cvAdjustParams(cv_mem);
    }

    /* Check if this step should be projected */
    let mut do_projection = SUNFALSE;
    if cv_mem.proj_enabled {
        let pm = cv_mem.proj_mem.as_deref().unwrap();
        do_projection =
            pm.freq > 0 && (cv_mem.cv_nst == 0 || cv_mem.cv_nst >= pm.nstlprj + pm.freq);
    }

    /* Looping point for attempts to take a step */

    let saved_t = cv_mem.cv_tn; /* tn is updated in cvPredict */
    let mut nflag = FIRST_CALL;
    let mut kflag;
    let mut dsm = ZERO;
    let mut dsmQ = ZERO;
    let mut dsmS = ZERO;
    let mut dsmQS = ZERO;

    loop {
        cvPredict(cv_mem);
        cvSet(cv_mem);

        /* ------ Correct state variables ------ */

        nflag = cvNls(cv_mem, nflag);
        {
            let mut ncfn_l = cv_mem.cv_ncfn;
            kflag = cvHandleNFlag(cv_mem, &mut nflag, saved_t, &mut ncf, &mut ncfn_l);
            cv_mem.cv_ncfn = ncfn_l;
        }

        /* Go back in loop if we need to predict again (nflag=PREV_CONV_FAIL) */
        if kflag == PREDICT_AGAIN {
            continue;
        }

        /* Return if nonlinear solve failed and recovery is not possible. */
        if kflag != DO_ERROR_TEST {
            return kflag;
        }

        /* Check inequality constraints */
        if !cv_mem.cv_constraints.is_empty() {
            let cflag =
                cvCheckConstraints(cv_mem, &mut nflag, saved_t, &mut step_constraint_fails);

            /* Go back in loop if we need to predict again (nflag=PREV_CONV_FAIL) */
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
            let pflag = crate::cvodes_proj::cvDoProjection(cv_mem, &mut nflag, saved_t, &mut npf);

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
        let eflag = {
            let acnrm = cv_mem.cv_acnrm;
            let mut netf_l = cv_mem.cv_netf;
            let e = cvDoErrorTest(cv_mem, &mut nflag, saved_t, acnrm, &mut nef, &mut netf_l, &mut dsm);
            cv_mem.cv_netf = netf_l;
            e
        };

        /* Go back in loop if we need to predict again (nflag=PREV_ERR_FAIL) */
        if eflag == TRY_AGAIN {
            continue;
        }

        /* Return if error test failed and recovery is not possible. */
        if eflag != CV_SUCCESS {
            return eflag;
        }

        /* Error test passed (eflag=CV_SUCCESS, nflag=CV_SUCCESS), go on */

        /* ------ Correct the quadrature variables ------ */

        if cv_mem.cv_quadr {
            ncf = 0;
            nef = 0; /* reset counters for states */

            nflag = cvQuadNls(cv_mem);
            {
                let mut ncfn_l = cv_mem.cv_ncfn;
                kflag = cvHandleNFlag(cv_mem, &mut nflag, saved_t, &mut ncf, &mut ncfn_l);
                cv_mem.cv_ncfn = ncfn_l;
            }

            if kflag == PREDICT_AGAIN {
                continue;
            }
            if kflag != DO_ERROR_TEST {
                return kflag;
            }

            /* Error test on quadratures */
            if cv_mem.cv_errconQ {
                cv_mem.cv_acnrmQ = N_VWrmsNorm(&cv_mem.cv_acorQ, &cv_mem.cv_ewtQ);
                let eflag = {
                    let acnrmQ = cv_mem.cv_acnrmQ;
                    let mut netfQ_l = cv_mem.cv_netfQ;
                    let e = cvDoErrorTest(cv_mem, &mut nflag, saved_t, acnrmQ, &mut nefQ, &mut netfQ_l, &mut dsmQ);
                    cv_mem.cv_netfQ = netfQ_l;
                    e
                };

                if eflag == TRY_AGAIN {
                    continue;
                }
                if eflag != CV_SUCCESS {
                    return eflag;
                }

                /* Set dsm = max(dsm, dsmQ) to be used in cvPrepareNextStep */
                if dsmQ > dsm {
                    dsm = dsmQ;
                }
            }
        }

        /* ------ Correct the sensitivity variables (STAGGERED or STAGGERED1) ------- */

        if do_sensi_stg || do_sensi_stg1 {
            ncf = 0;
            nef = 0; /* reset counters for states     */
            if cv_mem.cv_quadr {
                nefQ = 0; /* reset counter for quadratures */
            }

            /* Evaluate f at converged y, needed for future evaluations of sens. RHS
             * If f() fails recoverably, treat it as a convergence failure and
             * attempt the step again */

            let retval = {
                let f = cv_mem.cv_f.unwrap();
                f(
                    cv_mem.cv_tn,
                    &cv_mem.cv_y,
                    &mut cv_mem.cv_ftemp,
                    &mut cv_mem.cv_user_data,
                )
            };
            cv_mem.cv_nfe += 1;

            if retval < 0 {
                return CV_RHSFUNC_FAIL;
            }
            if retval > 0 {
                nflag = PREV_CONV_FAIL;
                continue;
            }

            if do_sensi_stg {
                /* Nonlinear solve for sensitivities (all-at-once) */
                nflag = cvStgrNls(cv_mem);
                let mut ncfnS_l = cv_mem.cv_ncfnS;
                kflag = cvHandleNFlag(cv_mem, &mut nflag, saved_t, &mut ncfS, &mut ncfnS_l);
                cv_mem.cv_ncfnS = ncfnS_l;
            } else {
                /* Nonlinear solve for sensitivities (one-by-one) */
                for is in 0..cv_mem.cv_Ns as usize {
                    cv_mem.sens_solve_idx = is as i32;

                    nflag = cvStgr1Nls(cv_mem, is as i32);
                    let mut ncf_l = cv_mem.cv_ncfS1[is] as i32;
                    let mut ncfn_l = cv_mem.cv_ncfnS1[is];
                    kflag = cvHandleNFlag(cv_mem, &mut nflag, saved_t, &mut ncf_l, &mut ncfn_l);
                    cv_mem.cv_ncfS1[is] = ncf_l as i64;
                    cv_mem.cv_ncfnS1[is] = ncfn_l;
                    if kflag != DO_ERROR_TEST {
                        break;
                    }
                }
            }

            if kflag == PREDICT_AGAIN {
                continue;
            }
            if kflag != DO_ERROR_TEST {
                return kflag;
            }

            /* Error test on sensitivities */
            if cv_mem.cv_errconS {
                if !cv_mem.cv_acnrmScur {
                    cv_mem.cv_acnrmS = cvSensNorm(cv_mem, &cv_mem.cv_acorS, &cv_mem.cv_ewtS);
                }

                let eflag = {
                    let acnrmS = cv_mem.cv_acnrmS;
                    let mut netfS_l = cv_mem.cv_netfS;
                    let e = cvDoErrorTest(cv_mem, &mut nflag, saved_t, acnrmS, &mut nefS, &mut netfS_l, &mut dsmS);
                    cv_mem.cv_netfS = netfS_l;
                    e
                };

                if eflag == TRY_AGAIN {
                    continue;
                }
                if eflag != CV_SUCCESS {
                    return eflag;
                }

                /* Set dsm = max(dsm, dsmS) to be used in cvPrepareNextStep */
                if dsmS > dsm {
                    dsm = dsmS;
                }
            }
        }

        /* ------ Correct the quadrature sensitivity variables ------ */

        if cv_mem.cv_quadr_sensi {
            /* Reset local convergence and error test failure counters */
            ncf = 0;
            nef = 0;
            if cv_mem.cv_quadr {
                nefQ = 0;
            }
            if do_sensi_stg {
                ncfS = 0;
                nefS = 0;
            }
            if do_sensi_stg1 {
                for is in 0..cv_mem.cv_Ns as usize {
                    cv_mem.cv_ncfS1[is] = 0;
                }
                nefS = 0;
            }

            /* Note that ftempQ contains yQdot evaluated at the converged y
             * (stored in cvQuadNls) and can be used in evaluating fQS */

            nflag = cvQuadSensNls(cv_mem);
            {
                let mut ncfn_l = cv_mem.cv_ncfn;
                kflag = cvHandleNFlag(cv_mem, &mut nflag, saved_t, &mut ncf, &mut ncfn_l);
                cv_mem.cv_ncfn = ncfn_l;
            }

            if kflag == PREDICT_AGAIN {
                continue;
            }
            if kflag != DO_ERROR_TEST {
                return kflag;
            }

            /* Error test on quadrature sensitivities */
            if cv_mem.cv_errconQS {
                cv_mem.cv_acnrmQS = cvQuadSensNorm(cv_mem, &cv_mem.cv_acorQS, &cv_mem.cv_ewtQS);
                let eflag = {
                    let acnrmQS = cv_mem.cv_acnrmQS;
                    let mut netfQS_l = cv_mem.cv_netfQS;
                    let e = cvDoErrorTest(cv_mem, &mut nflag, saved_t, acnrmQS, &mut nefQS, &mut netfQS_l, &mut dsmQS);
                    cv_mem.cv_netfQS = netfQS_l;
                    e
                };

                if eflag == TRY_AGAIN {
                    continue;
                }
                if eflag != CV_SUCCESS {
                    return eflag;
                }

                /* Set dsm = max(dsm, dsmQS) to be used in cvPrepareNextStep */
                if dsmQS > dsm {
                    dsm = dsmQS;
                }
            }
        }

        /* Error test passed (eflag=CV_SUCCESS), break from loop */
        break;
    }

    /* Nonlinear system solve and error test were both successful.
       Update data, and consider change of step and/or order.       */

    cvCompleteStep(cv_mem);

    cvPrepareNextStep(cv_mem, dsm);

    /* If Stablilty Limit Detection is turned on, call stability limit
       detection routine for possible order reduction. */

    if cv_mem.cv_sldeton {
        cvBDFStab(cv_mem);
    }

    cv_mem.cv_etamax = if cv_mem.cv_nst <= cv_mem.cv_small_nst {
        cv_mem.cv_eta_max_es
    } else {
        cv_mem.cv_eta_max_gs
    };

    /*  Finally, we rescale the acor array to be the
        estimated local error vector. */

    let tq2 = cv_mem.cv_tq[2];
    cv_mem.cv_acor.scale_inplace(tq2);

    if cv_mem.cv_quadr {
        cv_mem.cv_acorQ.scale_inplace(tq2);
    }

    if cv_mem.cv_sensi {
        for is in 0..cv_mem.cv_Ns as usize {
            cv_mem.cv_acorS[is].scale_inplace(tq2);
        }
    }

    if cv_mem.cv_quadr_sensi {
        for is in 0..cv_mem.cv_Ns as usize {
            cv_mem.cv_acorQS[is].scale_inplace(tq2);
        }
    }

    CV_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function called at beginning of step
 * -----------------------------------------------------------------
 */

/*
 * cvAdjustParams
 *
 * This routine is called when a change in step size was decided upon,
 * and it handles the required adjustments to the history array zn.
 * If there is to be a change in order, we call cvAdjustOrder and reset
 * q, L = q+1, and qwait.  Then in any case, we call cvRescale, which
 * resets h and rescales the Nordsieck array.
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
 *
 * This routine is a high level routine which handles an order
 * change by an amount deltaq (= +1 or -1). If a decrease in order
 * is requested and q==2, then the routine returns immediately.
 * Otherwise cvAdjustAdams or cvAdjustBDF is called to handle the
 * order change (depending on the value of lmm).
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
 * This routine adjusts the history array on a change of order q by
 * deltaq, in the case that lmm == CV_ADAMS.
 */
fn cvAdjustAdams(cv_mem: &mut CVodeMem, deltaq: i32) {
    /* On an order increase, set new column of zn to zero and return */
    if deltaq == 1 {
        let l = cv_mem.cv_L as usize;
        N_VConst(ZERO, &mut cv_mem.cv_zn[l]);
        if cv_mem.cv_quadr {
            N_VConst(ZERO, &mut cv_mem.cv_znQ[l]);
        }
        if cv_mem.cv_sensi {
            for is in 0..cv_mem.cv_Ns as usize {
                N_VConst(ZERO, &mut cv_mem.cv_znS[l][is]);
            }
        }
        return;
    }

    /*
     * On an order decrease, each zn[j] is adjusted by a multiple of zn[q].
     * The coeffs. in the adjustment are the coeffs. of the polynomial:
     *        x
     * q * INT { u * ( u + xi_1 ) * ... * ( u + xi_{q-2} ) } du
     *        0
     * where xi_j = [t_n - t_(n-j)]/h => xi_0 = 0
     */

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
        {
            let (front, back) = cv_mem.cv_zn.split_at_mut(q);
            let znq = &back[0];
            for j in 2..q {
                let c = -cv_mem.cv_l[j];
                for (z, x) in front[j].data.iter_mut().zip(&znq.data) {
                    *z += c * *x;
                }
            }
        }

        if cv_mem.cv_quadr {
            let (front, back) = cv_mem.cv_znQ.split_at_mut(q);
            let znQq = &back[0];
            for j in 2..q {
                let c = -cv_mem.cv_l[j];
                for (z, x) in front[j].data.iter_mut().zip(&znQq.data) {
                    *z += c * *x;
                }
            }
        }

        if cv_mem.cv_sensi {
            let ns = cv_mem.cv_Ns as usize;
            let (front, back) = cv_mem.cv_znS.split_at_mut(q);
            let znSq = &back[0];
            for j in 2..q {
                let c = -cv_mem.cv_l[j];
                for is in 0..ns {
                    for (z, x) in front[j][is].data.iter_mut().zip(&znSq[is].data) {
                        *z += c * *x;
                    }
                }
            }
        }
    }
}

/*
 * cvAdjustBDF
 *
 * This is a high level routine which handles adjustments to the
 * history array on a change of order by deltaq in the case that
 * lmm == CV_BDF.  cvAdjustBDF calls cvIncreaseBDF if deltaq = +1 and
 * cvDecreaseBDF if deltaq = -1 to do the actual work.
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
 * This routine adjusts the history array on an increase in the
 * order q in the case that lmm == CV_BDF.
 * A new column zn[q+1] is set equal to a multiple of the saved
 * vector (= acor) in zn[indx_acor].  Then each zn[j] is adjusted by
 * a multiple of zn[q+1].  The coefficients in the adjustment are the
 * coefficients of the polynomial x*x*(x+xi_1)*...*(x+xi_j),
 * where xi_j = [t_n - t_(n-j)]/h.
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

    /*
       zn[indx_acor] contains the value Delta_n = y_n - y_n(0)
       This value was stored there at the previous successful
       step (in cvCompleteStep)

       A1 contains dbar = (1/xi* - 1/xi_q)/prod(xi_j)
    */

    let indx = cv_mem.cv_indx_acor as usize;
    let l = cv_mem.cv_L as usize;
    let q = cv_mem.cv_q as usize;
    let ns = cv_mem.cv_Ns as usize;

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

    /* zn[j] += l[j]*zn[L] for j = 2..=q (N_VScaleAddMulti in C) */
    if cv_mem.cv_q > 1 {
        let (front, back) = cv_mem.cv_zn.split_at_mut(l);
        let znl = &back[0];
        for j in 2..=q {
            let c = cv_mem.cv_l[j];
            for (z, x) in front[j].data.iter_mut().zip(&znl.data) {
                *z += c * *x;
            }
        }
    }

    if cv_mem.cv_quadr {
        /* znQ[L] = A1 * znQ[indx_acor] */
        if indx == l {
            cv_mem.cv_znQ[l].scale_inplace(a1);
        } else {
            let (lo, hi) = if indx < l { (indx, l) } else { (l, indx) };
            let (front, back) = cv_mem.cv_znQ.split_at_mut(hi);
            let (src, dst) = if indx < l {
                (&front[lo], &mut back[0])
            } else {
                let tmp = &mut front[lo];
                (&back[0], tmp)
            };
            N_VScale(a1, src, dst);
        }

        /* znQ[j] += l[j]*znQ[L] for j = 2..=q */
        if cv_mem.cv_q > 1 {
            let (front, back) = cv_mem.cv_znQ.split_at_mut(l);
            let znQl = &back[0];
            for j in 2..=q {
                let c = cv_mem.cv_l[j];
                for (z, x) in front[j].data.iter_mut().zip(&znQl.data) {
                    *z += c * *x;
                }
            }
        }
    }

    if cv_mem.cv_sensi {
        /* znS[L][is] = A1 * znS[indx_acor][is] */
        if indx == l {
            for is in 0..ns {
                cv_mem.cv_znS[l][is].scale_inplace(a1);
            }
        } else {
            let (lo, hi) = if indx < l { (indx, l) } else { (l, indx) };
            let (front, back) = cv_mem.cv_znS.split_at_mut(hi);
            let (src_row, dst_row) = if indx < l {
                (&front[lo], &mut back[0])
            } else {
                let tmp = &mut front[lo];
                (&back[0], tmp)
            };
            for is in 0..ns {
                N_VScale(a1, &src_row[is], &mut dst_row[is]);
            }
        }

        /* znS[j][is] += l[j]*znS[L][is] for j = 2..=q */
        if cv_mem.cv_q > 1 {
            let (front, back) = cv_mem.cv_znS.split_at_mut(l);
            let znSl = &back[0];
            for j in 2..=q {
                let c = cv_mem.cv_l[j];
                for is in 0..ns {
                    for (z, x) in front[j][is].data.iter_mut().zip(&znSl[is].data) {
                        *z += c * *x;
                    }
                }
            }
        }
    }

    if cv_mem.cv_quadr_sensi {
        /* znQS[L][is] = A1 * znQS[indx_acor][is] */
        if indx == l {
            for is in 0..ns {
                cv_mem.cv_znQS[l][is].scale_inplace(a1);
            }
        } else {
            let (lo, hi) = if indx < l { (indx, l) } else { (l, indx) };
            let (front, back) = cv_mem.cv_znQS.split_at_mut(hi);
            let (src_row, dst_row) = if indx < l {
                (&front[lo], &mut back[0])
            } else {
                let tmp = &mut front[lo];
                (&back[0], tmp)
            };
            for is in 0..ns {
                N_VScale(a1, &src_row[is], &mut dst_row[is]);
            }
        }

        /* znQS[j][is] += l[j]*znQS[L][is] for j = 2..=q */
        if cv_mem.cv_q > 1 {
            let (front, back) = cv_mem.cv_znQS.split_at_mut(l);
            let znQSl = &back[0];
            for j in 2..=q {
                let c = cv_mem.cv_l[j];
                for is in 0..ns {
                    for (z, x) in front[j][is].data.iter_mut().zip(&znQSl[is].data) {
                        *z += c * *x;
                    }
                }
            }
        }
    }
}

/*
 * cvDecreaseBDF
 *
 * This routine adjusts the history array on a decrease in the
 * order q in the case that lmm == CV_BDF.
 * Each zn[j] is adjusted by a multiple of zn[q].  The coefficients
 * in the adjustment are the coefficients of the polynomial
 *   x*x*(x+xi_1)*...*(x+xi_j), where xi_j = [t_n - t_(n-j)]/h.
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
        {
            let (front, back) = cv_mem.cv_zn.split_at_mut(q);
            let znq = &back[0];
            for j in 2..q {
                let c = -cv_mem.cv_l[j];
                for (z, x) in front[j].data.iter_mut().zip(&znq.data) {
                    *z += c * *x;
                }
            }
        }

        if cv_mem.cv_quadr {
            let (front, back) = cv_mem.cv_znQ.split_at_mut(q);
            let znQq = &back[0];
            for j in 2..q {
                let c = -cv_mem.cv_l[j];
                for (z, x) in front[j].data.iter_mut().zip(&znQq.data) {
                    *z += c * *x;
                }
            }
        }

        if cv_mem.cv_sensi {
            let ns = cv_mem.cv_Ns as usize;
            let (front, back) = cv_mem.cv_znS.split_at_mut(q);
            let znSq = &back[0];
            for j in 2..q {
                let c = -cv_mem.cv_l[j];
                for is in 0..ns {
                    for (z, x) in front[j][is].data.iter_mut().zip(&znSq[is].data) {
                        *z += c * *x;
                    }
                }
            }
        }

        if cv_mem.cv_quadr_sensi {
            let ns = cv_mem.cv_Ns as usize;
            let (front, back) = cv_mem.cv_znQS.split_at_mut(q);
            let znQSq = &back[0];
            for j in 2..q {
                let c = -cv_mem.cv_l[j];
                for is in 0..ns {
                    for (z, x) in front[j][is].data.iter_mut().zip(&znQSq[is].data) {
                        *z += c * *x;
                    }
                }
            }
        }
    }
}

/*
 * cvRescale
 *
 * This routine rescales the Nordsieck array by multiplying the
 * jth column zn[j] by eta^j, j = 1, ..., q.  Then the value of
 * h is rescaled by eta, and hscale is reset to h.
 */
pub fn cvRescale(cv_mem: &mut CVodeMem) {
    /* compute scaling factors sequentially (cvals[j] = eta^j in C)
       and scale the columns (N_VScaleVectorArray) */
    let ns = cv_mem.cv_Ns as usize;
    let mut factor = cv_mem.cv_eta;
    for j in 1..=(cv_mem.cv_q as usize) {
        cv_mem.cv_zn[j].scale_inplace(factor);

        if cv_mem.cv_quadr {
            cv_mem.cv_znQ[j].scale_inplace(factor);
        }

        if cv_mem.cv_sensi {
            for is in 0..ns {
                cv_mem.cv_znS[j][is].scale_inplace(factor);
            }
        }

        if cv_mem.cv_quadr_sensi {
            for is in 0..ns {
                cv_mem.cv_znQS[j][is].scale_inplace(factor);
            }
        }

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
 * This routine advances tn by the tentative step size h, and computes
 * the predicted array z_n(0), which is overwritten on zn.  The
 * prediction of zn is done by repeated additions.
 * If tstop is enabled, it is possible for tn + h to be past tstop by roundoff,
 * and in that case, we reset tn (after incrementing by h) to tstop.
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

    if cv_mem.cv_quadr {
        for k in 1..=q {
            let mut j = q;
            while j >= k {
                /* znQ[j-1] += znQ[j] */
                let (front, back) = cv_mem.cv_znQ.split_at_mut(j as usize);
                let znQj = &back[0];
                for (z, x) in front[(j - 1) as usize].data.iter_mut().zip(&znQj.data) {
                    *z += *x;
                }
                j -= 1;
            }
        }
    }

    if cv_mem.cv_sensi {
        let ns = cv_mem.cv_Ns as usize;
        for k in 1..=q {
            let mut j = q;
            while j >= k {
                /* znS[j-1][is] += znS[j][is] */
                let (front, back) = cv_mem.cv_znS.split_at_mut(j as usize);
                let znSj = &back[0];
                for is in 0..ns {
                    for (z, x) in front[(j - 1) as usize][is].data.iter_mut().zip(&znSj[is].data) {
                        *z += *x;
                    }
                }
                j -= 1;
            }
        }
    }

    if cv_mem.cv_quadr_sensi {
        let ns = cv_mem.cv_Ns as usize;
        for k in 1..=q {
            let mut j = q;
            while j >= k {
                /* znQS[j-1][is] += znQS[j][is] */
                let (front, back) = cv_mem.cv_znQS.split_at_mut(j as usize);
                let znQSj = &back[0];
                for is in 0..ns {
                    for (z, x) in front[(j - 1) as usize][is].data.iter_mut().zip(&znQSj[is].data) {
                        *z += *x;
                    }
                }
                j -= 1;
            }
        }
    }
}

/*
 * cvSet
 *
 * This routine is a high level routine which calls cvSetAdams or
 * cvSetBDF to set the polynomial l, the test quantity array tq,
 * and the related variables  rl1, gamma, and gamrat.
 *
 * The array tq is loaded with constants used in the control of estimated
 * local errors and in the nonlinear convergence test.  Specifically, while
 * running at order q, the components of tq are as follows:
 *   tq[1] = a coefficient used to get the est. local error at order q-1
 *   tq[2] = a coefficient used to get the est. local error at order q
 *   tq[3] = a coefficient used to get the est. local error at order q+1
 *   tq[4] = constant used in nonlinear iteration convergence test
 *   tq[5] = coefficient used to get the order q+2 derivative vector used in
 *           the est. local error at order q+1
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
 *
 * This routine handles the computation of l and tq for the
 * case lmm == CV_ADAMS.
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
 * This routine generates in m[] the coefficients of the product
 * polynomial needed for the Adams l and tq coefficients for q > 1.
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
 * This routine completes the calculation of the Adams l and tq.
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
 * cvAltSum returns the value of the alternating sum
 *   sum (i= 0 ... iend) [ (-1)^i * (a[i] / (i + k)) ].
 * If iend < 0 then cvAltSum returns 0.
 * This operation is needed to compute the integral, from -1 to 0,
 * of a polynomial x^(k-1) M(x) given the coefficients of M(x).
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
 *
 * This routine computes the coefficients l and tq in the case
 * lmm == CV_BDF.  cvSetBDF calls cvSetTqBDF to set the test
 * quantity array tq.
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
 *
 * This routine sets the test quantity array tq in the case
 * lmm == CV_BDF.
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
 * -----------------------------------------------------------------
 * Nonlinear solver functions
 * -----------------------------------------------------------------
 */

/*
 * cvNls
 *
 * This routine attempts to solve the nonlinear system associated
 * with a single implicit step of the linear multistep method.
 */
fn cvNls(cv_mem: &mut CVodeMem, nflag: i32) -> i32 {
    let call_setup;

    /* Are we computing sensitivities with the CV_SIMULTANEOUS approach? */
    let do_sensi_sim = cv_mem.cv_sensi && cv_mem.cv_ism == CV_SIMULTANEOUS;

    /* Decide whether or not to call setup routine (if one exists) and */
    /* set flag convfail (input to lsetup for its evaluation decision) */
    if crate::cvodes_nls::cv_has_lsetup(cv_mem) {
        cv_mem.convfail = if nflag == FIRST_CALL || nflag == PREV_ERR_FAIL {
            CV_NO_FAILURES
        } else {
            CV_FAIL_OTHER
        };

        let mut cs = nflag == PREV_CONV_FAIL
            || nflag == PREV_ERR_FAIL
            || cv_mem.cv_nst == 0
            || cv_mem.first_step_after_resize
            || cv_mem.cv_nst >= cv_mem.cv_nstlp + cv_mem.cv_msbp
            || SUNRabs(cv_mem.cv_gamrat - ONE) > cv_mem.cv_dgmax_lsetup;

        /* Decide whether to force a call to setup */
        if cv_mem.cv_forceSetup {
            cs = SUNTRUE;
            cv_mem.convfail = CV_FAIL_OTHER;
        }

        call_setup = cs;
    } else {
        cv_mem.cv_crate = ONE;
        cv_mem.cv_crateS = ONE; /* if NO lsetup all conv. rates are set to ONE */
        call_setup = SUNFALSE;
    }

    /* initial guess for the correction to the predictor
       (for the simultaneous corrector, ycorSim = [acor, acorS]) */
    if do_sensi_sim {
        N_VConst(ZERO, &mut cv_mem.cv_acor);
        for is in 0..cv_mem.cv_Ns as usize {
            N_VConst(ZERO, &mut cv_mem.cv_acorS[is]);
        }
    } else {
        N_VConst(ZERO, &mut cv_mem.cv_acor);
    }

    /* (Newton and fixed-point solvers have no setup operation) */

    /* solve the nonlinear system */
    let tol = cv_mem.cv_tq[4];
    let flag;
    if do_sensi_sim {
        let mut nls = cv_mem
            .NLSsim
            .take()
            .expect("simultaneous-corrector nonlinear solver attached");
        flag = crate::cvodes_nls_sim::cvNlsSolveSensSim(cv_mem, &mut nls, tol, call_setup);

        /* increment counters */
        cv_mem.cv_nni += nls.get_num_iters();
        cv_mem.cv_nnf += nls.get_num_conv_fails();
        cv_mem.NLSsim = Some(nls);
    } else {
        let mut nls = cv_mem.NLS.take().expect("nonlinear solver attached");
        flag = crate::cvodes_nls::cvNlsSolve(cv_mem, &mut nls, tol, call_setup);

        /* increment counters */
        cv_mem.cv_nni += nls.get_num_iters();
        cv_mem.cv_nnf += nls.get_num_conv_fails();
        cv_mem.NLS = Some(nls);
    }

    /* if the solve failed return */
    if flag != 0 {
        return flag;
    }

    /* solve successful */

    /* update the state based on the final correction from the nonlinear solver */
    {
        let CVodeMem { cv_zn, cv_acor, cv_y, .. } = cv_mem;
        N_VLinearSum(ONE, &cv_zn[0], ONE, cv_acor, cv_y);
    }

    /* update the sensitivities based on the final correction from the nonlinear solver */
    if do_sensi_sim {
        let ns = cv_mem.cv_Ns as usize;
        let CVodeMem { cv_znS, cv_acorS, cv_yS, .. } = cv_mem;
        for is in 0..ns {
            N_VLinearSum(ONE, &cv_znS[0][is], ONE, &cv_acorS[is], &mut cv_yS[is]);
        }
    }

    /* compute acnrm if is was not already done by the nonlinear solver */
    if !cv_mem.cv_acnrmcur {
        if do_sensi_sim && cv_mem.cv_errconS {
            /* N_VWrmsNorm(ycorSim, ewtSim): the senswrapper WRMS norm is the
               max over the sub-vector norms */
            let del = N_VWrmsNorm(&cv_mem.cv_acor, &cv_mem.cv_ewt);
            cv_mem.cv_acnrm = cvSensUpdateNorm(cv_mem, del, &cv_mem.cv_acorS, &cv_mem.cv_ewtS);
        } else {
            cv_mem.cv_acnrm = N_VWrmsNorm(&cv_mem.cv_acor, &cv_mem.cv_ewt);
        }
    }

    /* update Jacobian status */
    cv_mem.cv_jcur = SUNFALSE;

    flag
}

/*
 * cvCheckConstraints
 *
 * This routine determines if the constraints of the problem
 * are satisfied by the proposed step
 *
 * Possible return values are:
 *
 *   CV_SUCCESS     ---> allows stepping forward
 *
 *   PREDICT_AGAIN  ---> values failed to satisfy constraints
 *
 *   CV_CONSTR_FAIL ---> values failed to satisfy constraints with hmin
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
 * cvQuadNls
 *
 * This routine solves for the quadrature variables at the new step.
 * It does not solve a nonlinear system, but rather updates the
 * quadrature variables. The name for this function is just for
 * uniformity purposes.
 *
 * Possible return values (interpreted by cvHandleNFlag)
 *
 *   CV_SUCCESS       -> continue with error test
 *   CV_QRHSFUNC_FAIL -> halt the integration
 *   QRHSFUNC_RECVR   -> predict again or stop if too many
 *
 */
fn cvQuadNls(cv_mem: &mut CVodeMem) -> i32 {
    /* Save quadrature correction in acorQ */
    let retval = {
        let fQ = cv_mem.cv_fQ.unwrap();
        fQ(
            cv_mem.cv_tn,
            &cv_mem.cv_y,
            &mut cv_mem.cv_acorQ,
            &mut cv_mem.cv_user_data,
        )
    };
    cv_mem.cv_nfQe += 1;
    if retval < 0 {
        return CV_QRHSFUNC_FAIL;
    }
    if retval > 0 {
        return QRHSFUNC_RECVR;
    }

    /* If needed, save the value of yQdot = fQ into ftempQ
     * for use in evaluating fQS */
    if cv_mem.cv_quadr_sensi {
        let CVodeMem { cv_acorQ, cv_ftempQ, .. } = cv_mem;
        N_VScale(ONE, cv_acorQ, cv_ftempQ);
    }

    {
        let h = cv_mem.cv_h;
        let rl1 = cv_mem.cv_rl1;
        let CVodeMem { cv_acorQ, cv_znQ, .. } = cv_mem;
        /* acorQ = h*acorQ - znQ[1]; acorQ *= rl1 */
        cv_acorQ.linear_sum_with(h, -ONE, &cv_znQ[1]);
        cv_acorQ.scale_inplace(rl1);
    }

    /* Apply correction to quadrature variables */
    {
        let CVodeMem { cv_znQ, cv_acorQ, cv_yQ, .. } = cv_mem;
        N_VLinearSum(ONE, &cv_znQ[0], ONE, cv_acorQ, cv_yQ);
    }

    CV_SUCCESS
}

/*
 * cvQuadSensNls
 *
 * This routine solves for the quadrature sensitivity variables
 * at the new step. It does not solve a nonlinear system, but
 * rather updates the quadrature variables. The name for this
 * function is just for uniformity purposes.
 *
 * Possible return values (interpreted by cvHandleNFlag)
 *
 *   CV_SUCCESS        -> continue with error test
 *   CV_QSRHSFUNC_FAIL -> halt the integration
 *   QSRHSFUNC_RECVR   -> predict again or stop if too many
 *
 */
fn cvQuadSensNls(cv_mem: &mut CVodeMem) -> i32 {
    /* Save quadrature sensitivity correction in acorQS */
    let y = std::mem::take(&mut cv_mem.cv_y);
    let yS = std::mem::take(&mut cv_mem.cv_yS);
    let ftempQ = std::mem::take(&mut cv_mem.cv_ftempQ);
    let mut acorQS = std::mem::take(&mut cv_mem.cv_acorQS);
    let mut tempv = std::mem::take(&mut cv_mem.cv_tempv);
    let mut tempvQ = std::mem::take(&mut cv_mem.cv_tempvQ);
    let tn = cv_mem.cv_tn;
    let retval = cv_fQS_dispatch(cv_mem, tn, &y, &yS, &ftempQ, &mut acorQS, &mut tempv, &mut tempvQ);
    cv_mem.cv_y = y;
    cv_mem.cv_yS = yS;
    cv_mem.cv_ftempQ = ftempQ;
    cv_mem.cv_acorQS = acorQS;
    cv_mem.cv_tempv = tempv;
    cv_mem.cv_tempvQ = tempvQ;
    cv_mem.cv_nfQSe += 1;
    if retval < 0 {
        return CV_QSRHSFUNC_FAIL;
    }
    if retval > 0 {
        return QSRHSFUNC_RECVR;
    }

    let h = cv_mem.cv_h;
    let rl1 = cv_mem.cv_rl1;
    let ns = cv_mem.cv_Ns as usize;
    let CVodeMem { cv_acorQS, cv_znQS, cv_yQS, .. } = cv_mem;
    for is in 0..ns {
        /* acorQS[is] = h*acorQS[is] - znQS[1][is]; acorQS[is] *= rl1 */
        cv_acorQS[is].linear_sum_with(h, -ONE, &cv_znQS[1][is]);
        cv_acorQS[is].scale_inplace(rl1);
        /* Apply correction to quadrature sensitivity variables */
        N_VLinearSum(ONE, &cv_znQS[0][is], ONE, &cv_acorQS[is], &mut cv_yQS[is]);
    }

    CV_SUCCESS
}

/*
 * cvStgrNls
 *
 * This is a high-level routine that attempts to solve the
 * sensitivity linear systems using the attached nonlinear solver
 * once the states y_n were obtained and passed the error test.
 */
fn cvStgrNls(cv_mem: &mut CVodeMem) -> i32 {
    let call_setup = SUNFALSE;
    if !crate::cvodes_nls::cv_has_lsetup(cv_mem) {
        cv_mem.cv_crateS = ONE;
    }

    /* initial guess for the correction to the predictor */
    for is in 0..cv_mem.cv_Ns as usize {
        N_VConst(ZERO, &mut cv_mem.cv_acorS[is]);
    }

    /* set sens solve flag */
    cv_mem.sens_solve = SUNTRUE;

    /* solve the nonlinear system */
    let tol = cv_mem.cv_tq[4];
    let mut nls = cv_mem
        .NLSstg
        .take()
        .expect("staggered-corrector nonlinear solver attached");
    let flag = crate::cvodes_nls_stg::cvNlsSolveSensStg(cv_mem, &mut nls, tol, call_setup);

    /* increment counters */
    cv_mem.cv_nniS += nls.get_num_iters();
    cv_mem.cv_nnfS += nls.get_num_conv_fails();
    cv_mem.NLSstg = Some(nls);

    /* reset sens solve flag */
    cv_mem.sens_solve = SUNFALSE;

    /* if the solve failed return */
    if flag != 0 {
        return flag;
    }

    /* solve successful */

    /* update the sensitivities based on the final correction from the nonlinear solver */
    {
        let ns = cv_mem.cv_Ns as usize;
        let CVodeMem { cv_znS, cv_acorS, cv_yS, .. } = cv_mem;
        for is in 0..ns {
            N_VLinearSum(ONE, &cv_znS[0][is], ONE, &cv_acorS[is], &mut cv_yS[is]);
        }
    }

    /* update Jacobian status */
    cv_mem.cv_jcur = SUNFALSE;

    flag
}

/*
 * cvStgr1Nls
 *
 * This is a high-level routine that attempts to solve the i-th
 * sensitivity linear system using the attached nonlinear solver
 * once the states y_n were obtained and passed the error test.
 */
fn cvStgr1Nls(cv_mem: &mut CVodeMem, is: i32) -> i32 {
    let isu = is as usize;
    let call_setup = SUNFALSE;
    if !crate::cvodes_nls::cv_has_lsetup(cv_mem) {
        cv_mem.cv_crateS = ONE;
    }

    /* initial guess for the correction to the predictor */
    N_VConst(ZERO, &mut cv_mem.cv_acorS[isu]);

    /* set sens solve flag */
    cv_mem.sens_solve = SUNTRUE;

    /* solve the nonlinear system */
    let tol = cv_mem.cv_tq[4];
    let mut nls = cv_mem
        .NLSstg1
        .take()
        .expect("staggered1-corrector nonlinear solver attached");
    let flag = crate::cvodes_nls_stg1::cvNlsSolveSensStg1(cv_mem, &mut nls, tol, call_setup);

    /* increment counters */
    cv_mem.cv_nniS1[isu] += nls.get_num_iters();
    cv_mem.cv_nnfS1[isu] += nls.get_num_conv_fails();
    cv_mem.NLSstg1 = Some(nls);

    /* reset sens solve flag */
    cv_mem.sens_solve = SUNFALSE;

    /* if the solve failed return */
    if flag != 0 {
        return flag;
    }

    /* solve successful */

    /* update the sensitivity with the final correction from the nonlinear solver */
    {
        let CVodeMem { cv_znS, cv_acorS, cv_yS, .. } = cv_mem;
        N_VLinearSum(ONE, &cv_znS[0][isu], ONE, &cv_acorS[isu], &mut cv_yS[isu]);
    }

    /* update Jacobian status */
    cv_mem.cv_jcur = SUNFALSE;

    flag
}

/*
 * cvHandleNFlag
 *
 * This routine takes action on the return value nflag = *nflagPtr
 * returned by cvNls, as follows:
 *
 * If cvNls succeeded in solving the nonlinear system, then
 * cvHandleNFlag returns the constant DO_ERROR_TEST, which tells cvStep
 * to perform the error test.
 *
 * If the nonlinear system was not solved successfully, then ncfn and
 * ncf = *ncfPtr are incremented and Nordsieck array zn is restored.
 *
 * If the solution of the nonlinear system failed due to an
 * unrecoverable failure by setup, we return the value CV_LSETUP_FAIL.
 *
 * If it failed due to an unrecoverable failure in solve, then we return
 * the value CV_LSOLVE_FAIL.
 *
 * If it failed due to an unrecoverable failure in rhs, then we return
 * the value CV_RHSFUNC_FAIL / CV_QRHSFUNC_FAIL / CV_SRHSFUNC_FAIL /
 * CV_QSRHSFUNC_FAIL.
 *
 * Otherwise, a recoverable failure occurred when solving the nonlinear
 * system (cvNls returned SUN_NLS_CONV_RECVR, RHSFUNC_RECVR, or
 * SRHSFUNC_RECVR).
 *
 * If ncf is now equal to maxncf or |h| = hmin, we return the value
 * CV_CONV_FAILURE (if SUN_NLS_CONV_RECVR),
 * CV_REPTD_RHSFUNC_ERR (if RHSFUNC_RECVR), or
 * CV_REPTD_SRHSFUNC_ERR (if SRHSFUNC_RECVR).
 * Otherwise, we set *nflagPtr = PREV_CONV_FAIL and return the value
 * PREDICT_AGAIN, telling cvStep to reattempt the step.
 *
 */
fn cvHandleNFlag(
    cv_mem: &mut CVodeMem,
    nflag_ptr: &mut i32,
    saved_t: f64,
    ncf_ptr: &mut i32,
    ncfn_ptr: &mut i64,
) -> i32 {
    let nflag = *nflag_ptr;

    if nflag == CV_SUCCESS {
        return DO_ERROR_TEST;
    }

    /* The nonlinear soln. failed; increment ncfn and restore zn */
    *ncfn_ptr += 1;
    cvRestore(cv_mem, saved_t);

    /* Return if failed unrecoverably */
    if nflag < 0 {
        if nflag == CV_LSETUP_FAIL {
            return CV_LSETUP_FAIL;
        } else if nflag == CV_LSOLVE_FAIL {
            return CV_LSOLVE_FAIL;
        } else if nflag == CV_RHSFUNC_FAIL {
            return CV_RHSFUNC_FAIL;
        } else if nflag == CV_QRHSFUNC_FAIL {
            return CV_QRHSFUNC_FAIL;
        } else if nflag == CV_SRHSFUNC_FAIL {
            return CV_SRHSFUNC_FAIL;
        } else if nflag == CV_QSRHSFUNC_FAIL {
            return CV_QSRHSFUNC_FAIL;
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
        if nflag == QRHSFUNC_RECVR {
            return CV_REPTD_QRHSFUNC_ERR;
        }
        if nflag == SRHSFUNC_RECVR {
            return CV_REPTD_SRHSFUNC_ERR;
        }
        if nflag == QSRHSFUNC_RECVR {
            return CV_REPTD_QSRHSFUNC_ERR;
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
 * This routine restores the value of tn to saved_t and undoes the
 * prediction.  After execution of cvRestore, the Nordsieck array zn has
 * the same values as before the call to cvPredict.
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

    if cv_mem.cv_quadr {
        for k in 1..=q {
            let mut j = q;
            while j >= k {
                /* znQ[j-1] -= znQ[j] */
                let (front, back) = cv_mem.cv_znQ.split_at_mut(j as usize);
                let znQj = &back[0];
                for (z, x) in front[(j - 1) as usize].data.iter_mut().zip(&znQj.data) {
                    *z -= *x;
                }
                j -= 1;
            }
        }
    }

    if cv_mem.cv_sensi {
        let ns = cv_mem.cv_Ns as usize;
        for k in 1..=q {
            let mut j = q;
            while j >= k {
                /* znS[j-1][is] -= znS[j][is] */
                let (front, back) = cv_mem.cv_znS.split_at_mut(j as usize);
                let znSj = &back[0];
                for is in 0..ns {
                    for (z, x) in front[(j - 1) as usize][is].data.iter_mut().zip(&znSj[is].data) {
                        *z -= *x;
                    }
                }
                j -= 1;
            }
        }
    }

    if cv_mem.cv_quadr_sensi {
        let ns = cv_mem.cv_Ns as usize;
        for k in 1..=q {
            let mut j = q;
            while j >= k {
                /* znQS[j-1][is] -= znQS[j][is] */
                let (front, back) = cv_mem.cv_znQS.split_at_mut(j as usize);
                let znQSj = &back[0];
                for is in 0..ns {
                    for (z, x) in front[(j - 1) as usize][is].data.iter_mut().zip(&znQSj[is].data) {
                        *z -= *x;
                    }
                }
                j -= 1;
            }
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
 * This routine performs the local error test, for the state, quadrature,
 * or sensitivity variables. Its last three arguments change depending
 * on which variables the error test is to be performed on.
 *
 * The weighted local error norm dsm is loaded into *dsmPtr, and
 * the test dsm ?<= 1 is made.
 *
 * If the test passes, cvDoErrorTest returns CV_SUCCESS.
 *
 * If the test fails, we undo the step just taken (call cvRestore) and
 *
 *   - if maxnef error test failures have occurred or if SUNRabs(h) = hmin,
 *     we return CV_ERR_FAILURE.
 *
 *   - if more than MXNEF1 error test failures have occurred, an order
 *     reduction is forced. If already at order 1, restart by reloading
 *     zn from scratch (also znQ and znS if appropriate).
 *     If f() fails, we return CV_RHSFUNC_FAIL or CV_UNREC_RHSFUNC_ERR;
 *     if fQ() fails, we return CV_QRHSFUNC_FAIL or CV_UNREC_QRHSFUNC_ERR;
 *     if cvSensRhsWrapper() fails, we return CV_SRHSFUNC_FAIL or
 *     CV_UNREC_SRHSFUNC_ERR; (no recovery is possible at this stage).
 *
 *   - otherwise, set *nflagPtr to PREV_ERR_FAIL, and return TRY_AGAIN.
 *
 */
fn cvDoErrorTest(
    cv_mem: &mut CVodeMem,
    nflag_ptr: &mut i32,
    saved_t: f64,
    acor_nrm: f64,
    nef_ptr: &mut i32,
    netf_ptr: &mut i64,
    dsm_ptr: &mut f64,
) -> i32 {
    let dsm = acor_nrm * cv_mem.cv_tq[2];

    /* If est. local error norm dsm passes test, return CV_SUCCESS */
    *dsm_ptr = dsm;
    if dsm <= ONE {
        return CV_SUCCESS;
    }

    /* Test failed; increment counters, set nflag, and restore zn array */
    *nef_ptr += 1;
    *netf_ptr += 1;
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

    /* If already at order 1, restart: reload zn, znQ, znS, znQS from scratch */

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
        let h = cv_mem.cv_h;
        let CVodeMem { cv_tempv, cv_zn, .. } = cv_mem;
        N_VScale(h, cv_tempv, &mut cv_zn[1]);
    }

    if cv_mem.cv_quadr {
        let fQ = cv_mem.cv_fQ.unwrap();
        let retval = fQ(
            cv_mem.cv_tn,
            &cv_mem.cv_zn[0],
            &mut cv_mem.cv_tempvQ,
            &mut cv_mem.cv_user_data,
        );
        cv_mem.cv_nfQe += 1;
        if retval < 0 {
            return CV_QRHSFUNC_FAIL;
        }
        if retval > 0 {
            return CV_UNREC_QRHSFUNC_ERR;
        }

        let h = cv_mem.cv_h;
        let CVodeMem { cv_tempvQ, cv_znQ, .. } = cv_mem;
        N_VScale(h, cv_tempvQ, &mut cv_znQ[1]);
    }

    if cv_mem.cv_sensi {
        /* wrk1 = ftemp, wrk2 = ftempS[0] */
        let tn = cv_mem.cv_tn;
        let zn0 = std::mem::take(&mut cv_mem.cv_zn[0]);
        let tempv = std::mem::take(&mut cv_mem.cv_tempv);
        let znS0 = std::mem::take(&mut cv_mem.cv_znS[0]);
        let mut tempvS = std::mem::take(&mut cv_mem.cv_tempvS);
        let mut wrk1 = std::mem::take(&mut cv_mem.cv_ftemp);
        let mut wrk2 = std::mem::take(&mut cv_mem.cv_ftempS[0]);
        let retval = cvSensRhsWrapper(cv_mem, tn, &zn0, &tempv, &znS0, &mut tempvS, &mut wrk1, &mut wrk2);
        cv_mem.cv_zn[0] = zn0;
        cv_mem.cv_tempv = tempv;
        cv_mem.cv_znS[0] = znS0;
        cv_mem.cv_tempvS = tempvS;
        cv_mem.cv_ftemp = wrk1;
        cv_mem.cv_ftempS[0] = wrk2;
        if retval < 0 {
            return CV_SRHSFUNC_FAIL;
        }
        if retval > 0 {
            return CV_UNREC_SRHSFUNC_ERR;
        }

        let ns = cv_mem.cv_Ns as usize;
        let h = cv_mem.cv_h;
        let CVodeMem { cv_tempvS, cv_znS, .. } = cv_mem;
        for is in 0..ns {
            N_VScale(h, &cv_tempvS[is], &mut cv_znS[1][is]);
        }
    }

    if cv_mem.cv_quadr_sensi {
        /* wrk1 = ftemp, wrk2 = ftempQ */
        let tn = cv_mem.cv_tn;
        let zn0 = std::mem::take(&mut cv_mem.cv_zn[0]);
        let znS0 = std::mem::take(&mut cv_mem.cv_znS[0]);
        let tempvQ = std::mem::take(&mut cv_mem.cv_tempvQ);
        let mut tempvQS = std::mem::take(&mut cv_mem.cv_tempvQS);
        let mut wrk1 = std::mem::take(&mut cv_mem.cv_ftemp);
        let mut wrk2 = std::mem::take(&mut cv_mem.cv_ftempQ);
        let retval = cv_fQS_dispatch(cv_mem, tn, &zn0, &znS0, &tempvQ, &mut tempvQS, &mut wrk1, &mut wrk2);
        cv_mem.cv_zn[0] = zn0;
        cv_mem.cv_znS[0] = znS0;
        cv_mem.cv_tempvQ = tempvQ;
        cv_mem.cv_tempvQS = tempvQS;
        cv_mem.cv_ftemp = wrk1;
        cv_mem.cv_ftempQ = wrk2;
        cv_mem.cv_nfQSe += 1;
        if retval < 0 {
            return CV_QSRHSFUNC_FAIL;
        }
        if retval > 0 {
            return CV_UNREC_QSRHSFUNC_ERR;
        }

        let ns = cv_mem.cv_Ns as usize;
        let h = cv_mem.cv_h;
        let CVodeMem { cv_tempvQS, cv_znQS, .. } = cv_mem;
        for is in 0..ns {
            N_VScale(h, &cv_tempvQS[is], &mut cv_znQS[1][is]);
        }
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
 * This routine performs various update operations when the solution
 * to the nonlinear system has passed the local error test.
 * We increment the step counter nst, record the values hu and qu,
 * update the tau array, and apply the corrections to the zn array.
 * The tau[i] are the last q values of h, with tau[1] the most recent.
 * The counter qwait is decremented, and if qwait == 1 (and q < qmax)
 * we save acor and tq[5] for a possible order increase.
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

    let q = cv_mem.cv_q as usize;
    let ns = cv_mem.cv_Ns as usize;

    /* Apply correction to column j of zn: l_j * Delta_n */
    {
        let CVodeMem { cv_zn, cv_acor, cv_l, .. } = cv_mem;
        for j in 0..=q {
            let c = cv_l[j];
            for (z, a) in cv_zn[j].data.iter_mut().zip(&cv_acor.data) {
                *z += c * *a;
            }
        }
    }

    /* Apply the projection correction to column j of zn: p_j * Delta_n */
    if cv_mem.proj_applied {
        let CVodeMem { cv_zn, cv_tempv, proj_p, .. } = cv_mem;
        for j in 0..=q {
            let c = proj_p[j];
            for (z, a) in cv_zn[j].data.iter_mut().zip(&cv_tempv.data) {
                /* tempv = acorP */
                *z += c * *a;
            }
        }
    }

    if cv_mem.cv_quadr {
        let CVodeMem { cv_znQ, cv_acorQ, cv_l, .. } = cv_mem;
        for j in 0..=q {
            let c = cv_l[j];
            for (z, a) in cv_znQ[j].data.iter_mut().zip(&cv_acorQ.data) {
                *z += c * *a;
            }
        }
    }

    if cv_mem.cv_sensi {
        let CVodeMem { cv_znS, cv_acorS, cv_l, .. } = cv_mem;
        for j in 0..=q {
            let c = cv_l[j];
            for is in 0..ns {
                for (z, a) in cv_znS[j][is].data.iter_mut().zip(&cv_acorS[is].data) {
                    *z += c * *a;
                }
            }
        }
    }

    if cv_mem.cv_quadr_sensi {
        let CVodeMem { cv_znQS, cv_acorQS, cv_l, .. } = cv_mem;
        for j in 0..=q {
            let c = cv_l[j];
            for is in 0..ns {
                for (z, a) in cv_znQS[j][is].data.iter_mut().zip(&cv_acorQS[is].data) {
                    *z += c * *a;
                }
            }
        }
    }

    /* If necessary, store Delta_n in zn[qmax] to be used in order increase.
     * This actually will be Delta_{n-1} in the ELTE at q+1 since it happens at
     * the next to last step of order q before a possible one at order q+1
     */

    cv_mem.cv_qwait -= 1;
    if cv_mem.cv_qwait == 1 && cv_mem.cv_q != cv_mem.cv_qmax {
        let qmax = cv_mem.cv_qmax as usize;

        {
            let CVodeMem { cv_zn, cv_acor, .. } = cv_mem;
            cv_zn[qmax].data.copy_from_slice(&cv_acor.data);
        }

        if cv_mem.cv_quadr {
            let CVodeMem { cv_znQ, cv_acorQ, .. } = cv_mem;
            cv_znQ[qmax].data.copy_from_slice(&cv_acorQ.data);
        }

        if cv_mem.cv_sensi {
            let CVodeMem { cv_znS, cv_acorS, .. } = cv_mem;
            for is in 0..ns {
                cv_znS[qmax][is].data.copy_from_slice(&cv_acorS[is].data);
            }
        }

        if cv_mem.cv_quadr_sensi {
            let CVodeMem { cv_znQS, cv_acorQS, .. } = cv_mem;
            for is in 0..ns {
                cv_znQS[qmax][is].data.copy_from_slice(&cv_acorQS[is].data);
            }
        }

        cv_mem.cv_saved_tq5 = cv_mem.cv_tq[5];
        cv_mem.cv_indx_acor = cv_mem.cv_qmax;
    }
}

/*
 * cvPrepareNextStep
 *
 * This routine handles the setting of stepsize and order for the
 * next step -- hprime and qprime.  Along with hprime, it sets the
 * ratio eta = hprime/h.  It also updates other state variables
 * related to a change of step size or order.
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
            /* If qwait = 0, consider an order change.   etaqm1 and etaqp1 are
              the ratios of new to old h at orders q-1 and q+1, respectively.
              cvChooseEta selects the largest; cvSetEta adjusts eta and acor */
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
 * This routine adjusts the value of eta according to the various
 * heuristic limits and the optional input hmax.
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
 *
 * This routine computes and returns the value of etaqm1 for a
 * possible decrease in order by 1.
 */
fn cvComputeEtaqm1(cv_mem: &mut CVodeMem) -> f64 {
    cv_mem.cv_etaqm1 = ZERO;
    if cv_mem.cv_q > 1 {
        let q = cv_mem.cv_q as usize;
        let mut ddn = N_VWrmsNorm(&cv_mem.cv_zn[q], &cv_mem.cv_ewt);

        if cv_mem.cv_quadr && cv_mem.cv_errconQ {
            ddn = cvQuadUpdateNorm(cv_mem, ddn, &cv_mem.cv_znQ[q], &cv_mem.cv_ewtQ);
        }

        if cv_mem.cv_sensi && cv_mem.cv_errconS {
            ddn = cvSensUpdateNorm(cv_mem, ddn, &cv_mem.cv_znS[q], &cv_mem.cv_ewtS);
        }

        if cv_mem.cv_quadr_sensi && cv_mem.cv_errconQS {
            ddn = cvQuadSensUpdateNorm(cv_mem, ddn, &cv_mem.cv_znQS[q], &cv_mem.cv_ewtQS);
        }

        ddn *= cv_mem.cv_tq[1];
        cv_mem.cv_etaqm1 = ONE / (SUNRpowerR(BIAS1 * ddn, ONE / cv_mem.cv_q as f64) + ADDON);
    }
    cv_mem.cv_etaqm1
}

/*
 * cvComputeEtaqp1
 *
 * This routine computes and returns the value of etaqp1 for a
 * possible increase in order by 1.
 */
fn cvComputeEtaqp1(cv_mem: &mut CVodeMem) -> f64 {
    cv_mem.cv_etaqp1 = ZERO;
    if cv_mem.cv_q != cv_mem.cv_qmax {
        if cv_mem.cv_saved_tq5 == ZERO {
            return cv_mem.cv_etaqp1;
        }
        let cquot = (cv_mem.cv_tq[5] / cv_mem.cv_saved_tq5)
            * SUNRpowerI(cv_mem.cv_h / cv_mem.cv_tau[2], cv_mem.cv_L);
        let qmax = cv_mem.cv_qmax as usize;
        {
            let CVodeMem { cv_zn, cv_acor, cv_tempv, .. } = cv_mem;
            N_VLinearSum(-cquot, &cv_zn[qmax], ONE, cv_acor, cv_tempv);
        }
        let mut dup = N_VWrmsNorm(&cv_mem.cv_tempv, &cv_mem.cv_ewt);

        if cv_mem.cv_quadr && cv_mem.cv_errconQ {
            {
                let CVodeMem { cv_znQ, cv_acorQ, cv_tempvQ, .. } = cv_mem;
                N_VLinearSum(-cquot, &cv_znQ[qmax], ONE, cv_acorQ, cv_tempvQ);
            }
            dup = cvQuadUpdateNorm(cv_mem, dup, &cv_mem.cv_tempvQ, &cv_mem.cv_ewtQ);
        }

        if cv_mem.cv_sensi && cv_mem.cv_errconS {
            {
                let ns = cv_mem.cv_Ns as usize;
                let CVodeMem { cv_znS, cv_acorS, cv_tempvS, .. } = cv_mem;
                for is in 0..ns {
                    N_VLinearSum(-cquot, &cv_znS[qmax][is], ONE, &cv_acorS[is], &mut cv_tempvS[is]);
                }
            }
            dup = cvSensUpdateNorm(cv_mem, dup, &cv_mem.cv_tempvS, &cv_mem.cv_ewtS);
        }

        if cv_mem.cv_quadr_sensi && cv_mem.cv_errconQS {
            {
                let ns = cv_mem.cv_Ns as usize;
                let CVodeMem { cv_znQS, cv_acorQS, cv_tempvQS, .. } = cv_mem;
                for is in 0..ns {
                    N_VLinearSum(-cquot, &cv_znQS[qmax][is], ONE, &cv_acorQS[is], &mut cv_tempvQS[is]);
                }
            }
            /* (C calls cvSensUpdateNorm here, not cvQuadSensUpdateNorm) */
            dup = cvSensUpdateNorm(cv_mem, dup, &cv_mem.cv_tempvQS, &cv_mem.cv_ewtQS);
        }

        dup *= cv_mem.cv_tq[3];
        cv_mem.cv_etaqp1 = ONE / (SUNRpowerR(BIAS3 * dup, ONE / (cv_mem.cv_L + 1) as f64) + ADDON);
    }
    cv_mem.cv_etaqp1
}

/*
 * cvChooseEta
 * Given etaqm1, etaq, etaqp1 (the values of eta for qprime =
 * q - 1, q, or q + 1, respectively), this routine chooses the
 * maximum eta value, sets eta to that value, and sets qprime to the
 * corresponding value of q.  If there is a tie, the preference
 * order is to (1) keep the same order, then (2) decrease the order,
 * and finally (3) increase the order.  If the maximum eta value
 * is within the fixed step bounds, the order is kept unchanged and
 * eta is set to 1.
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
            /*
             * Store Delta_n in zn[qmax] to be used in order increase
             *
             * This happens at the last step of order q before an increase
             * to order q+1, so it represents Delta_n in the ELTE at q+1
             */

            let qmax = cv_mem.cv_qmax as usize;
            let ns = cv_mem.cv_Ns as usize;

            {
                let CVodeMem { cv_zn, cv_acor, .. } = cv_mem;
                cv_zn[qmax].data.copy_from_slice(&cv_acor.data);
            }

            if cv_mem.cv_quadr && cv_mem.cv_errconQ {
                let CVodeMem { cv_znQ, cv_acorQ, .. } = cv_mem;
                cv_znQ[qmax].data.copy_from_slice(&cv_acorQ.data);
            }

            if cv_mem.cv_sensi && cv_mem.cv_errconS {
                let CVodeMem { cv_znS, cv_acorS, .. } = cv_mem;
                for is in 0..ns {
                    cv_znS[qmax][is].data.copy_from_slice(&cv_acorS[is].data);
                }
            }

            if cv_mem.cv_quadr_sensi && cv_mem.cv_errconQS {
                let CVodeMem { cv_znQS, cv_acorQS, .. } = cv_mem;
                for is in 0..ns {
                    cv_znQS[qmax][is].data.copy_from_slice(&cv_acorQS[is].data);
                }
            }
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
 * This routine prints error messages for all cases of failure by
 * cvHin and cvStep.
 * It returns to CVode the value that CVode is to return to the user.
 */
fn cvHandleFailure(cv_mem: &mut CVodeMem, flag: i32) -> i32 {
    /* Depending on flag, print error message and return error flag */
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
        CV_QRHSFUNC_FAIL => cvProcessError(Some(cv_mem), CV_QRHSFUNC_FAIL, line!(), "CVode", file!(),
            &format!("At t = {}, the quadrature right-hand side routine failed in an unrecoverable manner.",
                     cv_mem.cv_tn)),
        CV_UNREC_QRHSFUNC_ERR => cvProcessError(Some(cv_mem), CV_UNREC_QRHSFUNC_ERR, line!(), "CVode", file!(),
            &format!("At t = {}, the quadrature right-hand side failed in a recoverable manner, but no recovery is possible.",
                     cv_mem.cv_tn)),
        CV_REPTD_QRHSFUNC_ERR => cvProcessError(Some(cv_mem), CV_REPTD_QRHSFUNC_ERR, line!(), "CVode", file!(),
            &format!("At t = {} repeated recoverable quadrature right-hand side function errors.", cv_mem.cv_tn)),
        CV_SRHSFUNC_FAIL => cvProcessError(Some(cv_mem), CV_SRHSFUNC_FAIL, line!(), "CVode", file!(),
            &format!("At t = {}, the sensitivity right-hand side routine failed in an unrecoverable manner.",
                     cv_mem.cv_tn)),
        CV_UNREC_SRHSFUNC_ERR => cvProcessError(Some(cv_mem), CV_UNREC_SRHSFUNC_ERR, line!(), "CVode", file!(),
            &format!("At t = {}, the sensitivity right-hand side failed in a recoverable manner, but no recovery is possible.",
                     cv_mem.cv_tn)),
        CV_REPTD_SRHSFUNC_ERR => cvProcessError(Some(cv_mem), CV_REPTD_SRHSFUNC_ERR, line!(), "CVode", file!(),
            &format!("At t = {} repeated recoverable sensitivity right-hand side function errors.", cv_mem.cv_tn)),
        CV_QSRHSFUNC_FAIL => cvProcessError(Some(cv_mem), CV_QSRHSFUNC_FAIL, line!(), "CVode", file!(),
            &format!("At t = {}, the quadrature sensitivity right-hand side routine failed in an unrecoverable manner.",
                     cv_mem.cv_tn)),
        CV_UNREC_QSRHSFUNC_ERR => cvProcessError(Some(cv_mem), CV_UNREC_QSRHSFUNC_ERR, line!(), "CVode", file!(),
            &format!("At t = {}, the quadrature sensitivity right-hand side failed in a recoverable manner, but no recovery is possible.",
                     cv_mem.cv_tn)),
        CV_REPTD_QSRHSFUNC_ERR => cvProcessError(Some(cv_mem), CV_REPTD_QSRHSFUNC_ERR, line!(), "CVode", file!(),
            &format!("At t = {} repeated recoverable quadrature sensitivity right-hand side function errors.",
                     cv_mem.cv_tn)),
        CV_TOO_CLOSE => cvProcessError(Some(cv_mem), CV_TOO_CLOSE, line!(), "CVode", file!(), MSGCV_TOO_CLOSE),
        CV_MEM_NULL => cvProcessError(None, CV_MEM_NULL, line!(), "CVode", file!(), MSGCV_NO_MEM),
        SUN_ERR_ARG_CORRUPT => {
            cvProcessError(Some(cv_mem), CV_MEM_NULL, line!(), "CVode", file!(),
                &format!("At t = {}, the nonlinear solver was passed a NULL input.", cv_mem.cv_tn));
        }
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
            /* This return should never happen */
            cvProcessError(Some(cv_mem), CV_UNRECOGNIZED_ERR, line!(), "CVode", file!(),
                "CVODES encountered an unrecognized error. Please report this to the SUNDIALS developers at sundials-users@llnl.gov");
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
 * This routine handles the BDF Stability Limit Detection Algorithm
 * STALD.  It is called if lmm = CV_BDF and the SLDET option is on.
 * If the order is 3 or more, the required norm data is saved.
 * If a decision to reduce order has not already been made, and
 * enough data has been saved, cvSLdet is called.  If it signals
 * a stability limit violation, the order is reduced, and the step
 * size is reset accordingly.
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
        /* Otherwise, let order increase happen, and
           reset stability limit counter, nscon. */
        cv_mem.cv_nscon = 0;
    }
}

/*
 * cvSLdet
 *
 * This routine detects stability limitation using stored scaled
 * derivatives data. cvSLdet returns the magnitude of the
 * dominate characteristic root, rr. The presence of a stability
 * limit is indicated by rr > "something a little less then 1.0",
 * and a positive kflag. This routine should only be called if
 * order is greater than or equal to 3, and data has been collected
 * for 5 time steps.
 *
 * Returned values:
 *    kflag = 1 -> Found stable characteristic root, normal matrix case
 *    kflag = 2 -> Found stable characteristic root, quartic solution
 *    kflag = 3 -> Found stable characteristic root, quartic solution,
 *                 with Newton correction
 *    kflag = 4 -> Found stability violation, normal matrix case
 *    kflag = 5 -> Found stability violation, quartic solution
 *    kflag = 6 -> Found stability violation, quartic solution,
 *                 with Newton correction
 *
 *    kflag < 0 -> No stability limitation,
 *                 or could not compute limitation.
 *
 *    kflag = -1 -> Min/max ratio of ssdat too small.
 *    kflag = -2 -> For normal matrix case, vmax > vrrt2*vrrt2
 *    kflag = -3 -> For normal matrix case, The three ratios
 *                  are inconsistent.
 *    kflag = -4 -> Small coefficient prevents elimination of quartics.
 *    kflag = -5 -> R value from quartics not consistent.
 *    kflag = -6 -> No corrected root passes test on qk values
 *    kflag = -7 -> Trouble solving for sigsq.
 *    kflag = -8 -> Trouble solving for B, or R via B.
 *    kflag = -9 -> R via sigsq[k] disagrees with R from data.
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

    /* The following are cutoffs and tolerances used by this routine */
    let rrcut = 0.98;
    let vrrtol = 1.0e-4;
    let vrrt2 = 5.0e-4;
    let sqtol = 1.0e-3;
    let rrtol = 1.0e-2;

    let mut rr; /* (C initializes rr = ZERO; every reachable path assigns it) */

    /*  Index k corresponds to the degree of the interpolating polynomial. */
    /*      k = 1 -> q-1          */
    /*      k = 2 -> q            */
    /*      k = 3 -> q+1          */

    /*  Index i is a backward-in-time index, i = 1 -> current time, */
    /*      i = 2 -> previous step, etc    */

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
 * This routine completes the initialization of rootfinding memory
 * information, and checks whether g has a zero both at and very near
 * the initial point of the IVP.
 *
 * This routine returns an int equal to:
 *  CV_RTFUNC_FAIL < 0 if the g function failed, or
 *  CV_SUCCESS     = 0 otherwise.
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

    /* We check now only the components of g which were exactly 0.0 at t0
     * to see if we can 'activate' them. */
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
 * This routine checks for exact zeros of g at the last root found,
 * if the last return was a root.  It then checks for a close pair of
 * zeros (an error condition), and for a new root at a nearby point.
 * The array glo = g(tlo) at the left endpoint of the search interval
 * is adjusted if necessary to assure that all g_i are nonzero
 * there, before returning to do a root search in the interval.
 *
 * On entry, tlo = tretlast is the last value of tret returned by
 * CVode.  This may be the previous tn, the previous tout value,
 * or the last root location.
 *
 * This routine returns an int equal to:
 *  CV_RTFUNC_FAIL  < 0 if the g function failed, or
 *  CLOSERT         = 3 if a close pair of zeros was found, or
 *  RTFOUND         = 1 if a new zero of g was found near tlo, or
 *  CV_SUCCESS      = 0 otherwise.
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

    /* Check for close roots (error return), for a new zero at tlo+smallh,
       and for a g_i that changed from zero to nonzero. */
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
 * This routine interfaces to cvRootfind to look for a root of g
 * between tlo and either tn or tout, whichever comes first.
 * Only roots beyond tlo in the direction of integration are sought.
 *
 * This routine returns an int equal to:
 *  CV_RTFUNC_FAIL  < 0 if the g function failed, or
 *  RTFOUND         = 1 if a root of g was found, or
 *  CV_SUCCESS      = 0 otherwise.
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

    /* Set ghi = g(thi) and call cvRootfind to search (tlo,thi) for roots. */
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
 * This routine solves for a root of g(t) between tlo and thi, if
 * one exists.  Only roots of odd multiplicity (i.e. with a change
 * of sign in one of the g_i), or exact zeros, are found.
 * Here the sign of tlo - thi is arbitrary, but if multiple roots
 * are found, the one closest to tlo is returned.
 *
 * The method used is the Illinois algorithm, a modified secant method.
 *
 * This routine returns an int equal to:
 *  CV_RTFUNC_FAIL  < 0 if the g function failed, or
 *  RTFOUND         = 1 if a root of g was found, or
 *  CV_SUCCESS      = 0 otherwise.
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

    /* If no sign change was found, reset trout and grout. Then return
       CV_SUCCESS if no zero was found, or set iroots and return RTFOUND. */
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

    /* A sign change was found.  Loop to locate nearest root. */
    let mut side = 0;
    let mut sideprev = -1;
    let gfun = cv_mem.cv_gfun.unwrap();
    loop {
        /* If interval size is already less than tolerance ttol, break. */
        if SUNRabs(cv_mem.cv_thi - cv_mem.cv_tlo) <= cv_mem.cv_ttol {
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

        /* Check to see in which subinterval g changes sign, and reset imax.
           Set side = 1 if sign change is on low side, or 2 if on high side. */
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
            continue; /* Return to looping point. */
        }

        if zroot {
            /* No sign change in (tlo,tmid), but g = 0 at tmid; return root tmid. */
            cv_mem.cv_thi = tmid;
            for i in 0..nrt {
                cv_mem.cv_ghi[i] = cv_mem.cv_grout[i];
            }
            break;
        }

        /* No sign change in (tlo,tmid), and no zero at tmid.
           Sign change must be in (tmid,thi).  Replace tlo with tmid. */
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
 * cvEwtSet
 *
 * This routine is responsible for setting the error weight vector ewt,
 * according to tol_type, as follows:
 *
 * (1) ewt[i] = 1 / (reltol * SUNRabs(ycur[i]) + abstol), i=0,...,neq-1
 *     if tol_type = CV_SS
 * (2) ewt[i] = 1 / (reltol * SUNRabs(ycur[i]) + abstol[i]), i=0,...,neq-1
 *     if tol_type = CV_SV
 *
 * cvEwtSet returns 0 if ewt is successfully set as above to a
 * positive vector and -1 otherwise. In the latter case, ewt is
 * considered undefined.
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

/*
 * cvEwtSetSS
 *
 * This routine sets ewt as described above in the case tol_type = CV_SS.
 * If the absolute tolerance is zero, it tests for non-positive components
 * before inverting. cvEwtSetSS returns 0 if ewt is successfully set to a
 * positive vector and -1 otherwise. In the latter case, ewt is considered
 * undefined.
 */
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

/*
 * cvEwtSetSV
 *
 * This routine sets ewt as described above in the case tol_type = CV_SV.
 * If any absolute tolerance is zero, it tests for non-positive components
 * before inverting. cvEwtSetSV returns 0 if ewt is successfully set to a
 * positive vector and -1 otherwise. In the latter case, ewt is considered
 * undefined.
 */
fn cvEwtSetSV(cv_mem: &CVodeMem, ycur: &NVector, weight: &mut NVector) -> i32 {
    N_VAbs(ycur, weight);
    weight.linear_sum_with(cv_mem.cv_reltol, ONE, &cv_mem.cv_Vabstol);
    if cv_mem.cv_atolmin0 && N_VMin(weight) <= ZERO {
        return -1;
    }
    weight.invert_inplace();
    0
}

/*
 * cvQuadEwtSet
 *
 */
fn cvQuadEwtSet(cv_mem: &CVodeMem, qcur: &NVector, weightQ: &mut NVector) -> i32 {
    match cv_mem.cv_itolQ {
        CV_SS => cvQuadEwtSetSS(cv_mem, qcur, weightQ),
        CV_SV => cvQuadEwtSetSV(cv_mem, qcur, weightQ),
        _ => 0,
    }
}

/*
 * cvQuadEwtSetSS
 *
 * (The C version computes in cv_tempvQ and inverts into weightQ;
 * computing directly in weightQ performs the identical arithmetic.)
 */
fn cvQuadEwtSetSS(cv_mem: &CVodeMem, qcur: &NVector, weightQ: &mut NVector) -> i32 {
    N_VAbs(qcur, weightQ);
    weightQ.scale_inplace(cv_mem.cv_reltolQ);
    weightQ.add_const_inplace(cv_mem.cv_SabstolQ);
    if cv_mem.cv_atolQmin0 && N_VMin(weightQ) <= ZERO {
        return -1;
    }
    weightQ.invert_inplace();
    0
}

/*
 * cvQuadEwtSetSV
 *
 */
fn cvQuadEwtSetSV(cv_mem: &CVodeMem, qcur: &NVector, weightQ: &mut NVector) -> i32 {
    N_VAbs(qcur, weightQ);
    weightQ.linear_sum_with(cv_mem.cv_reltolQ, ONE, &cv_mem.cv_VabstolQ);
    if cv_mem.cv_atolQmin0 && N_VMin(weightQ) <= ZERO {
        return -1;
    }
    weightQ.invert_inplace();
    0
}

/*
 * cvSensEwtSet
 *
 */
fn cvSensEwtSet(cv_mem: &mut CVodeMem, yScur: &[NVector], weightS: &mut [NVector]) -> i32 {
    match cv_mem.cv_itolS {
        CV_EE => cvSensEwtSetEE(cv_mem, yScur, weightS),
        CV_SS => cvSensEwtSetSS(cv_mem, yScur, weightS),
        CV_SV => cvSensEwtSetSV(cv_mem, yScur, weightS),
        _ => 0,
    }
}

/*
 * cvSensEwtSetEE
 *
 * In this case, the error weight vector for the i-th sensitivity is set to
 *
 * ewtS_i = pbar_i * efun(pbar_i*yS_i)
 *
 * In other words, the scaled sensitivity pbar_i * yS_i has the same error
 * weight vector calculation as the solution vector.
 *
 */
fn cvSensEwtSetEE(cv_mem: &mut CVodeMem, yScur: &[NVector], weightS: &mut [NVector]) -> i32 {
    /* Use tempvS[0] as temporary storage for the scaled sensitivity pyS */
    for is in 0..cv_mem.cv_Ns as usize {
        let pbari = cv_mem.cv_pbar[is];
        {
            let CVodeMem { cv_tempvS, .. } = cv_mem;
            N_VScale(pbari, &yScur[is], &mut cv_tempvS[0]);
        }
        let flag = if cv_mem.cv_user_efun {
            let efun = cv_mem.cv_efun.unwrap();
            let CVodeMem { cv_tempvS, cv_user_data, .. } = cv_mem;
            efun(&cv_tempvS[0], &mut weightS[is], cv_user_data)
        } else {
            cvEwtSet(cv_mem, &cv_mem.cv_tempvS[0], &mut weightS[is])
        };
        if flag != 0 {
            return -1;
        }
        weightS[is].scale_inplace(pbari);
    }
    0
}

/*
 * cvSensEwtSetSS
 *
 */
fn cvSensEwtSetSS(cv_mem: &mut CVodeMem, yScur: &[NVector], weightS: &mut [NVector]) -> i32 {
    for is in 0..cv_mem.cv_Ns as usize {
        N_VAbs(&yScur[is], &mut cv_mem.cv_tempv);
        let reltolS = cv_mem.cv_reltolS;
        cv_mem.cv_tempv.scale_inplace(reltolS);
        let sabstolS = cv_mem.cv_SabstolS[is];
        cv_mem.cv_tempv.add_const_inplace(sabstolS);
        if cv_mem.cv_atolSmin0[is] && N_VMin(&cv_mem.cv_tempv) <= ZERO {
            return -1;
        }
        N_VInv(&cv_mem.cv_tempv, &mut weightS[is]);
    }
    0
}

/*
 * cvSensEwtSetSV
 *
 */
fn cvSensEwtSetSV(cv_mem: &mut CVodeMem, yScur: &[NVector], weightS: &mut [NVector]) -> i32 {
    for is in 0..cv_mem.cv_Ns as usize {
        N_VAbs(&yScur[is], &mut cv_mem.cv_tempv);
        {
            let reltolS = cv_mem.cv_reltolS;
            let CVodeMem { cv_tempv, cv_VabstolS, .. } = cv_mem;
            cv_tempv.linear_sum_with(reltolS, ONE, &cv_VabstolS[is]);
        }
        if cv_mem.cv_atolSmin0[is] && N_VMin(&cv_mem.cv_tempv) <= ZERO {
            return -1;
        }
        N_VInv(&cv_mem.cv_tempv, &mut weightS[is]);
    }
    0
}

/*
 * cvQuadSensEwtSet
 *
 */
fn cvQuadSensEwtSet(cv_mem: &mut CVodeMem, yQScur: &[NVector], weightQS: &mut [NVector]) -> i32 {
    match cv_mem.cv_itolQS {
        CV_EE => cvQuadSensEwtSetEE(cv_mem, yQScur, weightQS),
        CV_SS => cvQuadSensEwtSetSS(cv_mem, yQScur, weightQS),
        CV_SV => cvQuadSensEwtSetSV(cv_mem, yQScur, weightQS),
        _ => 0,
    }
}

/*
 * cvQuadSensEwtSetEE
 *
 * In this case, the error weight vector for the i-th quadrature sensitivity
 * is set to
 *
 * ewtQS_i = pbar_i * cvQuadEwtSet(pbar_i*yQS_i)
 *
 * In other words, the scaled sensitivity pbar_i * yQS_i has the same error
 * weight vector calculation as the quadrature vector.
 *
 */
fn cvQuadSensEwtSetEE(cv_mem: &mut CVodeMem, yQScur: &[NVector], weightQS: &mut [NVector]) -> i32 {
    /* Use tempvQS[0] as temporary storage for the scaled sensitivity pyS */
    for is in 0..cv_mem.cv_Ns as usize {
        let pbari = cv_mem.cv_pbar[is];
        {
            let CVodeMem { cv_tempvQS, .. } = cv_mem;
            N_VScale(pbari, &yQScur[is], &mut cv_tempvQS[0]);
        }
        let flag = cvQuadEwtSet(cv_mem, &cv_mem.cv_tempvQS[0], &mut weightQS[is]);
        if flag != 0 {
            return -1;
        }
        weightQS[is].scale_inplace(pbari);
    }
    0
}

/*
 * cvQuadSensEwtSetSS
 *
 */
fn cvQuadSensEwtSetSS(cv_mem: &mut CVodeMem, yQScur: &[NVector], weightQS: &mut [NVector]) -> i32 {
    for is in 0..cv_mem.cv_Ns as usize {
        N_VAbs(&yQScur[is], &mut cv_mem.cv_tempvQ);
        let reltolQS = cv_mem.cv_reltolQS;
        cv_mem.cv_tempvQ.scale_inplace(reltolQS);
        let sabstolQS = cv_mem.cv_SabstolQS[is];
        cv_mem.cv_tempvQ.add_const_inplace(sabstolQS);
        if cv_mem.cv_atolQSmin0[is] && N_VMin(&cv_mem.cv_tempvQ) <= ZERO {
            return -1;
        }
        N_VInv(&cv_mem.cv_tempvQ, &mut weightQS[is]);
    }
    0
}

/*
 * cvQuadSensEwtSetSV
 *
 */
fn cvQuadSensEwtSetSV(cv_mem: &mut CVodeMem, yQScur: &[NVector], weightQS: &mut [NVector]) -> i32 {
    for is in 0..cv_mem.cv_Ns as usize {
        N_VAbs(&yQScur[is], &mut cv_mem.cv_tempvQ);
        {
            let reltolQS = cv_mem.cv_reltolQS;
            let CVodeMem { cv_tempvQ, cv_VabstolQS, .. } = cv_mem;
            cv_tempvQ.linear_sum_with(reltolQS, ONE, &cv_VabstolQS[is]);
        }
        if cv_mem.cv_atolQSmin0[is] && N_VMin(&cv_mem.cv_tempvQ) <= ZERO {
            return -1;
        }
        N_VInv(&cv_mem.cv_tempvQ, &mut weightQS[is]);
    }
    0
}

/*
 * -----------------------------------------------------------------
 * Functions for combined norms
 * -----------------------------------------------------------------
 */

/*
 * cvQuadUpdateNorm
 *
 * Updates the norm old_nrm to account for all quadratures.
 */
pub(crate) fn cvQuadUpdateNorm(_cv_mem: &CVodeMem, old_nrm: f64, xQ: &NVector, wQ: &NVector) -> f64 {
    let qnrm = N_VWrmsNorm(xQ, wQ);
    if old_nrm > qnrm {
        old_nrm
    } else {
        qnrm
    }
}

/*
 * cvSensNorm
 *
 * This routine returns the maximum over the weighted root mean
 * square norm of xS with weight vectors wS:
 *
 *  max { wrms(xS[0],wS[0]) ... wrms(xS[Ns-1],wS[Ns-1]) }
 *
 * Called by cvSensUpdateNorm or directly in the CV_STAGGERED approach
 * during the NLS solution and before the error test.
 */
pub(crate) fn cvSensNorm(cv_mem: &CVodeMem, xS: &[NVector], wS: &[NVector]) -> f64 {
    /* (N_VWrmsNormVectorArray + max reduction, translated per-is) */
    let mut nrm = N_VWrmsNorm(&xS[0], &wS[0]);
    for is in 1..cv_mem.cv_Ns as usize {
        let tmp = N_VWrmsNorm(&xS[is], &wS[is]);
        if tmp > nrm {
            nrm = tmp;
        }
    }
    nrm
}

/*
 * cvSensUpdateNorm
 *
 * Updates the norm old_nrm to account for all sensitivities.
 */
pub(crate) fn cvSensUpdateNorm(cv_mem: &CVodeMem, old_nrm: f64, xS: &[NVector], wS: &[NVector]) -> f64 {
    let snrm = cvSensNorm(cv_mem, xS, wS);
    if old_nrm > snrm {
        old_nrm
    } else {
        snrm
    }
}

/*
 * cvQuadSensNorm
 *
 * This routine returns the maximum over the weighted root mean
 * square norm of xQS with weight vectors wQS:
 *
 *  max { wrms(xQS[0],wS[0]) ... wrms(xQS[Ns-1],wS[Ns-1]) }
 *
 * Called by cvQuadSensUpdateNorm.
 */
fn cvQuadSensNorm(cv_mem: &CVodeMem, xQS: &[NVector], wQS: &[NVector]) -> f64 {
    /* (N_VWrmsNormVectorArray + max reduction, translated per-is) */
    let mut nrm = N_VWrmsNorm(&xQS[0], &wQS[0]);
    for is in 1..cv_mem.cv_Ns as usize {
        let tmp = N_VWrmsNorm(&xQS[is], &wQS[is]);
        if tmp > nrm {
            nrm = tmp;
        }
    }
    nrm
}

/*
 * cvQuadSensUpdateNorm
 *
 * Updates the norm old_nrm to account for all quadrature sensitivities.
 */
pub(crate) fn cvQuadSensUpdateNorm(cv_mem: &CVodeMem, old_nrm: f64, xQS: &[NVector], wQS: &[NVector]) -> f64 {
    let snrm = cvQuadSensNorm(cv_mem, xQS, wQS);
    if old_nrm > snrm {
        old_nrm
    } else {
        snrm
    }
}

/*
 * -----------------------------------------------------------------
 * Wrappers for sensitivity RHS
 * -----------------------------------------------------------------
 */

/*
 * cvSensRhsWrapper
 *
 * CVSensRhs is a high level routine that returns right hand side
 * of sensitivity equations. Depending on the 'ifS' flag, it either
 * calls directly the fS routine (ifS=CV_ALLSENS) or (if ifS=CV_ONESENS)
 * calls the fS1 routine in a loop over all sensitivities.
 *
 * CVSensRhs is called:
 *  (*) by CVode at the first step
 *  (*) by cvYddNorm if errcon=SUNTRUE
 *  (*) by the nonlinear solver if ism=CV_SIMULTANEOUS
 *  (*) by cvDoErrorTest when restarting from scratch
 *  (*) in the corrector loop if ism=CV_STAGGERED
 *  (*) by cvStgrDoErrorTest when restarting from scratch
 *
 * The return value is that of the sensitivity RHS function fS,
 *
 */
pub(crate) fn cvSensRhsWrapper(
    cv_mem: &mut CVodeMem,
    time: f64,
    ycur: &NVector,
    fcur: &NVector,
    yScur: &[NVector],
    fScur: &mut [NVector],
    temp1: &mut NVector,
    temp2: &mut NVector,
) -> i32 {
    let ns = cv_mem.cv_Ns;

    if cv_mem.cv_ifS == CV_ALLSENS {
        let retval = if cv_mem.cv_fSDQ {
            /* C: cv_fS = cvSensRhsInternalDQ with fS_data = cvode_mem */
            cvSensRhsInternalDQ(cv_mem, ns, time, ycur, fcur, yScur, fScur, temp1, temp2)
        } else {
            let fS = cv_mem.cv_fS.unwrap();
            fS(ns, time, ycur, fcur, yScur, fScur, &mut cv_mem.cv_user_data, temp1, temp2)
        };
        cv_mem.cv_nfSe += 1;
        retval
    } else {
        let mut retval = 0;
        for is in 0..ns as usize {
            retval = if cv_mem.cv_fSDQ {
                /* C: cv_fS1 = cvSensRhs1InternalDQ with fS_data = cvode_mem */
                cvSensRhs1InternalDQ(cv_mem, ns, time, ycur, fcur, is as i32, &yScur[is], &mut fScur[is], temp1, temp2)
            } else {
                let fS1 = cv_mem.cv_fS1.unwrap();
                fS1(ns, time, ycur, fcur, is as i32, &yScur[is], &mut fScur[is], &mut cv_mem.cv_user_data, temp1, temp2)
            };
            cv_mem.cv_nfSe += 1;
            if retval != 0 {
                break;
            }
        }
        retval
    }
}

/*
 * cvSensRhs1Wrapper
 *
 * cvSensRhs1Wrapper is a high level routine that returns right-hand
 * side of the is-th sensitivity equation.
 *
 * cvSensRhs1Wrapper is called only during the CV_STAGGERED1 corrector loop
 * (ifS must be CV_ONESENS, otherwise CVodeSensInit would have
 * issued an error message).
 *
 * The return value is that of the sensitivity RHS function fS1,
 */
pub(crate) fn cvSensRhs1Wrapper(
    cv_mem: &mut CVodeMem,
    time: f64,
    ycur: &NVector,
    fcur: &NVector,
    is: i32,
    yScur: &NVector,
    fScur: &mut NVector,
    temp1: &mut NVector,
    temp2: &mut NVector,
) -> i32 {
    let ns = cv_mem.cv_Ns;
    let retval = if cv_mem.cv_fSDQ {
        cvSensRhs1InternalDQ(cv_mem, ns, time, ycur, fcur, is, yScur, fScur, temp1, temp2)
    } else {
        let fS1 = cv_mem.cv_fS1.unwrap();
        fS1(ns, time, ycur, fcur, is, yScur, fScur, &mut cv_mem.cv_user_data, temp1, temp2)
    };
    cv_mem.cv_nfSe += 1;
    retval
}

/*
 * -----------------------------------------------------------------
 * Internal DQ approximations for sensitivity RHS
 * -----------------------------------------------------------------
 */

/*
 * cvSensRhsInternalDQ   - internal CVSensRhsFn
 *
 * cvSensRhsInternalDQ computes right hand side of all sensitivity equations
 * by finite differences
 */
pub(crate) fn cvSensRhsInternalDQ(
    cv_mem: &mut CVodeMem,
    Ns: i32,
    t: f64,
    y: &NVector,
    ydot: &NVector,
    yS: &[NVector],
    ySdot: &mut [NVector],
    ytemp: &mut NVector,
    ftemp: &mut NVector,
) -> i32 {
    for is in 0..Ns as usize {
        let retval =
            cvSensRhs1InternalDQ(cv_mem, Ns, t, y, ydot, is as i32, &yS[is], &mut ySdot[is], ytemp, ftemp);
        if retval != 0 {
            return retval;
        }
    }

    0
}

/* C's sens-DQ routines perturb the parameter THROUGH the user's own p
   array (CVodeSetSensParams stores the pointer, so the user RHS sees
   the perturbed value).  cv_p is an owned copy here; the perturbation
   is mirrored into the user data through the FSAUserData convention
   (sundials_types.rs) and restored the same way. */
fn cv_dq_set_p(cv_mem: &mut CVodeMem, which: usize, value: f64) {
    cv_mem.cv_p[which] = value;
    if let Some(d) = cv_mem.cv_user_data.as_mut() {
        if let Some(f) = d.downcast_mut::<FSAUserData>() {
            f.p[which] = value;
        }
    }
}

/*
 * cvSensRhs1InternalDQ   - internal CVSensRhs1Fn
 *
 * cvSensRhs1InternalDQ computes the right hand side of the is-th sensitivity
 * equation by finite differences
 *
 * cvSensRhs1InternalDQ returns 0 if successful. Otherwise it returns the
 * non-zero return value from f().
 */
fn cvSensRhs1InternalDQ(
    cv_mem: &mut CVodeMem,
    _Ns: i32,
    t: f64,
    y: &NVector,
    ydot: &NVector,
    is: i32,
    yS: &NVector,
    ySdot: &mut NVector,
    ytemp: &mut NVector,
    ftemp: &mut NVector,
) -> i32 {
    let mut nfel: i64 = 0;

    let delta = SUNRsqrt(SUNMAX(cv_mem.cv_reltol, cv_mem.cv_uround));
    let rdelta = ONE / delta;

    let pbari = cv_mem.cv_pbar[is as usize];

    let which = cv_mem.cv_plist[is as usize] as usize;

    let psave = cv_mem.cv_p[which];

    let Deltap = pbari * delta;
    let rDeltap = ONE / Deltap;
    let norms = N_VWrmsNorm(yS, &cv_mem.cv_ewt) * pbari;
    let rDeltay = SUNMAX(norms, rdelta) / pbari;
    let Deltay = ONE / rDeltay;

    let method = if cv_mem.cv_DQrhomax == ZERO {
        /* No switching */
        if cv_mem.cv_DQtype == CV_CENTERED { CENTERED1 } else { FORWARD1 }
    } else {
        /* switch between simultaneous/separate DQ */
        let ratio = Deltay * rDeltap;
        if SUNMAX(ONE / ratio, ratio) <= cv_mem.cv_DQrhomax {
            if cv_mem.cv_DQtype == CV_CENTERED { CENTERED1 } else { FORWARD1 }
        } else if cv_mem.cv_DQtype == CV_CENTERED {
            CENTERED2
        } else {
            FORWARD2
        }
    };

    /* (fn pointer is Copy; capture before mutating cv_p) */
    let f = cv_mem.cv_f.unwrap();

    match method {
        CENTERED1 => {
            let Delta = SUNMIN(Deltay, Deltap);
            let r2Delta = HALF / Delta;

            N_VLinearSum(ONE, y, Delta, yS, ytemp);
            cv_dq_set_p(cv_mem, which, psave + Delta);

            let retval = f(t, ytemp, ySdot, &mut cv_mem.cv_user_data);
            nfel += 1;
            if retval != 0 {
                return retval;
            }

            N_VLinearSum(ONE, y, -Delta, yS, ytemp);
            cv_dq_set_p(cv_mem, which, psave - Delta);

            let retval = f(t, ytemp, ftemp, &mut cv_mem.cv_user_data);
            nfel += 1;
            if retval != 0 {
                return retval;
            }

            /* ySdot = r2Delta*ySdot - r2Delta*ftemp (aliased) */
            ySdot.linear_sum_with(r2Delta, -r2Delta, ftemp);
        }

        CENTERED2 => {
            let r2Deltap = HALF / Deltap;
            let r2Deltay = HALF / Deltay;

            N_VLinearSum(ONE, y, Deltay, yS, ytemp);

            let retval = f(t, ytemp, ySdot, &mut cv_mem.cv_user_data);
            nfel += 1;
            if retval != 0 {
                return retval;
            }

            N_VLinearSum(ONE, y, -Deltay, yS, ytemp);

            let retval = f(t, ytemp, ftemp, &mut cv_mem.cv_user_data);
            nfel += 1;
            if retval != 0 {
                return retval;
            }

            /* ySdot = r2Deltay*ySdot - r2Deltay*ftemp (aliased) */
            ySdot.linear_sum_with(r2Deltay, -r2Deltay, ftemp);

            cv_dq_set_p(cv_mem, which, psave + Deltap);
            let retval = f(t, y, ytemp, &mut cv_mem.cv_user_data);
            nfel += 1;
            if retval != 0 {
                return retval;
            }

            cv_dq_set_p(cv_mem, which, psave - Deltap);
            let retval = f(t, y, ftemp, &mut cv_mem.cv_user_data);
            nfel += 1;
            if retval != 0 {
                return retval;
            }

            /* ySdot = ySdot + r2Deltap * ytemp - r2Deltap * ftemp
               (N_VLinearCombination(3, {1, r2Deltap, -r2Deltap},
                {ySdot, ytemp, ftemp}, ySdot)) */
            for j in 0..ySdot.data.len() {
                let a = ySdot.data[j];
                let b = ytemp.data[j];
                let c = ftemp.data[j];
                ySdot.data[j] = a + r2Deltap * b + (-r2Deltap) * c;
            }
        }

        FORWARD1 => {
            let Delta = SUNMIN(Deltay, Deltap);
            let rDelta = ONE / Delta;

            N_VLinearSum(ONE, y, Delta, yS, ytemp);
            cv_dq_set_p(cv_mem, which, psave + Delta);

            let retval = f(t, ytemp, ySdot, &mut cv_mem.cv_user_data);
            nfel += 1;
            if retval != 0 {
                return retval;
            }

            /* ySdot = rDelta*ySdot - rDelta*ydot (aliased) */
            ySdot.linear_sum_with(rDelta, -rDelta, ydot);
        }

        FORWARD2 => {
            N_VLinearSum(ONE, y, Deltay, yS, ytemp);

            let retval = f(t, ytemp, ySdot, &mut cv_mem.cv_user_data);
            nfel += 1;
            if retval != 0 {
                return retval;
            }

            /* ySdot = rDeltay*ySdot - rDeltay*ydot (aliased) */
            ySdot.linear_sum_with(rDeltay, -rDeltay, ydot);

            cv_dq_set_p(cv_mem, which, psave + Deltap);
            let retval = f(t, y, ytemp, &mut cv_mem.cv_user_data);
            nfel += 1;
            if retval != 0 {
                return retval;
            }

            /* ySdot = ySdot + rDeltap * ytemp - rDeltap * ydot
               (N_VLinearCombination(3, {1, rDeltap, -rDeltap},
                {ySdot, ytemp, ydot}, ySdot)) */
            for j in 0..ySdot.data.len() {
                let a = ySdot.data[j];
                let b = ytemp.data[j];
                let c = ydot.data[j];
                ySdot.data[j] = a + rDeltap * b + (-rDeltap) * c;
            }
        }

        _ => {}
    }

    cv_dq_set_p(cv_mem, which, psave);

    /* Increment counter nfeS */
    cv_mem.cv_nfeS += nfel;

    0
}

/*
 * cvQuadSensRhsInternalDQ   - internal CVQuadSensRhsFn
 *
 * cvQuadSensRhsInternalDQ computes right hand side of all quadrature
 * sensitivity equations by finite differences. All work is actually
 * done in cvQuadSensRhs1InternalDQ.
 */
pub(crate) fn cvQuadSensRhsInternalDQ(
    cv_mem: &mut CVodeMem,
    Ns: i32,
    t: f64,
    y: &NVector,
    yS: &[NVector],
    yQdot: &NVector,
    yQSdot: &mut [NVector],
    tmp: &mut NVector,
    tmpQ: &mut NVector,
) -> i32 {
    for is in 0..Ns as usize {
        let retval =
            cvQuadSensRhs1InternalDQ(cv_mem, is as i32, t, y, &yS[is], yQdot, &mut yQSdot[is], tmp, tmpQ);
        if retval != 0 {
            return retval;
        }
    }

    0
}

fn cvQuadSensRhs1InternalDQ(
    cv_mem: &mut CVodeMem,
    is: i32,
    t: f64,
    y: &NVector,
    yS: &NVector,
    yQdot: &NVector,
    yQSdot: &mut NVector,
    tmp: &mut NVector,
    tmpQ: &mut NVector,
) -> i32 {
    let mut nfel: i64 = 0;

    let delta = SUNRsqrt(SUNMAX(cv_mem.cv_reltol, cv_mem.cv_uround));
    let rdelta = ONE / delta;

    let pbari = cv_mem.cv_pbar[is as usize];

    let which = cv_mem.cv_plist[is as usize] as usize;

    let psave = cv_mem.cv_p[which];

    let Deltap = pbari * delta;
    let norms = N_VWrmsNorm(yS, &cv_mem.cv_ewt) * pbari;
    let rDeltay = SUNMAX(norms, rdelta) / pbari;
    let Deltay = ONE / rDeltay;

    let method = if cv_mem.cv_DQtype == CV_CENTERED { CENTERED1 } else { FORWARD1 };

    /* (fn pointer is Copy; capture before mutating cv_p) */
    let fQ = cv_mem.cv_fQ.unwrap();

    match method {
        CENTERED1 => {
            let Delta = SUNMIN(Deltay, Deltap);
            let r2Delta = HALF / Delta;

            N_VLinearSum(ONE, y, Delta, yS, tmp);
            cv_dq_set_p(cv_mem, which, psave + Delta);

            let retval = fQ(t, tmp, yQSdot, &mut cv_mem.cv_user_data);
            nfel += 1;
            if retval != 0 {
                return retval;
            }

            N_VLinearSum(ONE, y, -Delta, yS, tmp);
            cv_dq_set_p(cv_mem, which, psave - Delta);

            let retval = fQ(t, tmp, tmpQ, &mut cv_mem.cv_user_data);
            nfel += 1;
            if retval != 0 {
                return retval;
            }

            /* yQSdot = r2Delta*yQSdot - r2Delta*tmpQ (aliased) */
            yQSdot.linear_sum_with(r2Delta, -r2Delta, tmpQ);
        }

        FORWARD1 => {
            let Delta = SUNMIN(Deltay, Deltap);
            let rDelta = ONE / Delta;

            N_VLinearSum(ONE, y, Delta, yS, tmp);
            cv_dq_set_p(cv_mem, which, psave + Delta);

            let retval = fQ(t, tmp, yQSdot, &mut cv_mem.cv_user_data);
            nfel += 1;
            if retval != 0 {
                return retval;
            }

            /* yQSdot = rDelta*yQSdot - rDelta*yQdot (aliased) */
            yQSdot.linear_sum_with(rDelta, -rDelta, yQdot);
        }

        _ => {}
    }

    cv_dq_set_p(cv_mem, which, psave);

    /* Increment counter nfQeS */
    cv_mem.cv_nfQeS += nfel;

    0
}

// ===================== END PART 3 (cvodes.c:5874-10126) =====================
// (cvProcessError is defined in cvodes_impl.rs)
