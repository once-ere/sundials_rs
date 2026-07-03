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
