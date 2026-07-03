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
 *  pattern; cv_lreinit_dispatch already landed with PART 1.)
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
        LsModule::Ls(ls) => crate::cvodes_ls::cvLsInitialize(cv_mem, ls),
        LsModule::Diag(dm) => crate::cvodes_diag::cvDiagInit(cv_mem, dm),
    };
    cv_mem.cv_lmem = lmem;
    ier
}

pub(crate) fn cv_lsetup_dispatch(
    cv_mem: &mut CVodeMem,
    convfail: i32,
    jcur_ptr: &mut bool,
) -> i32 {
    let mut lmem = std::mem::take(&mut cv_mem.cv_lmem);
    let ier = match &mut lmem {
        LsModule::None => 0,
        LsModule::Ls(ls) => crate::cvodes_ls::cvLsSetup(cv_mem, ls, convfail, jcur_ptr),
        LsModule::Diag(dm) => crate::cvodes_diag::cvDiagSetup(cv_mem, dm, convfail, jcur_ptr),
    };
    cv_mem.cv_lmem = lmem;
    ier
}

/* The CVODES lsolve receives the weight vector to use: cv_ewt for
   state solves (ewtS_is = None) or cv_ewtS[is] for STAGGERED1
   sensitivity solves (ewtS_is = Some(is)); this mirrors the C
   cv_lsolve(cv_mem, b, weight, ycur, fcur) weight argument. */
pub(crate) fn cv_lsolve_dispatch(
    cv_mem: &mut CVodeMem,
    b: &mut NVector,
    ewtS_is: Option<usize>,
) -> i32 {
    let mut lmem = std::mem::take(&mut cv_mem.cv_lmem);
    let ier = match &mut lmem {
        LsModule::None => {
            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "cv_lsolve", file!(), MSGCV_LSOLVE_NULL);
            -1
        }
        LsModule::Ls(ls) => crate::cvodes_ls::cvLsSolve(cv_mem, ls, b, ewtS_is),
        LsModule::Diag(dm) => crate::cvodes_diag::cvDiagSolve(cv_mem, dm, b),
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
// PART 3 (step machinery: cvAdjustParams..cvBDFStab, nls drivers cvNls/
// cvStgrNls/cvStgr1Nls, error/convergence handling, ewt/sens norms, sens RHS
// wrappers + DQ, rootfinding cvRcheck/cvRootfind, cvProcessError —
// cvodes.c:6245-10126) is appended below by the next agent.
