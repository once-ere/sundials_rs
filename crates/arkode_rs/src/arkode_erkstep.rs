/* -----------------------------------------------------------------
 * Translation of sundials-7.7.0/src/arkode/arkode_erkstep.c
 * (PART I: everything except the discrete-adjoint machinery —
 * erkStep_TakeStep_Adjoint / erkStep_fe_Adj /
 * ERKStepCreateAdjointStepper need the ManyVector module and are
 * deferred — and erkStep_RelaxDeltaE, deferred with
 * arkode_relaxation.c).
 *
 * step_mem access: every op takes ark_mem.step_mem out (Box<dyn
 * Any> downcast to ARKodeERKStepMem), runs an _inner worker, and
 * puts it back (Addendum C.1). erkStep_TakeStep re-enters
 * erkStep_FullRHS through the step_fullrhs pointer, so the wrapper
 * releases step_mem around that call (the C statements executed
 * in-between do not touch step_mem).
 *
 * The C cvals/Xvecs scratch arrays feed N_VLinearCombination; here
 * the coefficients/operands are assembled in locals at each call
 * site (identical element-wise arithmetic); cvals/nfusedopvecs are
 * kept for the C workspace accounting.
 * -----------------------------------------------------------------*/

use crate::arkode::{arkAllocVecArray, arkFreeVecArray, arkInit, arkResizeVec};
use crate::arkode_butcher::{
    ARKodeButcherTable_IsStifflyAccurate, ARKodeButcherTable_Space, ARKodeButcherTable_Write,
};
use crate::arkode_butcher_erk::ARKodeButcherTable_LoadERK;
use crate::arkode_erkstep_impl::{
    ARKodeERKStepMem, ERKSTEP_DEFAULT_1, ERKSTEP_DEFAULT_2, ERKSTEP_DEFAULT_3,
    ERKSTEP_DEFAULT_4, ERKSTEP_DEFAULT_5, ERKSTEP_DEFAULT_6, ERKSTEP_DEFAULT_7,
    ERKSTEP_DEFAULT_8, ERKSTEP_DEFAULT_9, MSG_ERKSTEP_NO_MEM,
};
use crate::arkode_impl::{
    arkProcessError, ARKRhsFn, ARKVecResizeFn, ARKodeMem, ARK_ACCUMERROR_NONE,
    ARK_FULLRHS_END, ARK_FULLRHS_OTHER, ARK_FULLRHS_START, ARK_ILL_INPUT, ARK_INVALID_TABLE,
    ARK_MEM_NULL, ARK_NO_MALLOC, ARK_POSTPROCESS_STAGE_FAIL, ARK_POSTPROCESS_STEP_FAIL,
    ARK_PRERHSFN_FAIL, ARK_RHSFUNC_FAIL, ARK_SUCCESS, ARK_UNREC_RHSFUNC_ERR, ARK_VECTOROP_ERR,
    ARK_WARNING, FIRST_INIT, ONE, RESET_INIT, RESIZE_INIT, ZERO,
};
use crate::nvector_serial::{N_VLinearCombination, N_VScale, N_VWrmsNorm, NVector};
use crate::sundials_math::SUNRabs;
use crate::sundials_types::UserData;
use crate::sundials_utils::fmt_g;

/*---------------------------------------------------------------
  erkStep_AccessStepMem: takes the ERKStep memory out of ark_mem
  (callers must put it back).
  ---------------------------------------------------------------*/
pub(crate) fn erkStep_AccessStepMem(
    ark_mem: &mut ARKodeMem,
    fname: &str,
) -> Option<Box<ARKodeERKStepMem>> {
    match ark_mem.step_mem.take() {
        None => {
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_NULL,
                line!(),
                fname,
                file!(),
                MSG_ERKSTEP_NO_MEM,
            );
            None
        }
        Some(b) => match b.downcast::<ARKodeERKStepMem>() {
            Ok(sm) => Some(sm),
            Err(b) => {
                ark_mem.step_mem = Some(b);
                arkProcessError(
                    Some(ark_mem),
                    ARK_MEM_NULL,
                    line!(),
                    fname,
                    file!(),
                    MSG_ERKSTEP_NO_MEM,
                );
                None
            }
        },
    }
}

/*===============================================================
  Exported functions
  ===============================================================*/

/// C ERKStepCreate(f, t0, y0, sunctx) returns the ARKODE memory (the
/// C NULL-argument and allocation failures cannot occur here).
pub fn ERKStepCreate(
    f: ARKRhsFn,
    t0: f64,
    y0: &NVector,
    sunctx: &crate::sundials_context::SUNContext,
) -> Box<ARKodeMem> {
    /* Create ark_mem structure and set default values */
    let mut ark_mem = crate::arkode::arkCreate(sunctx);

    /* Allocate ARKodeERKStepMem structure, and initialize to zero */
    let mut step_mem = Box::new(ARKodeERKStepMem {
        f: None,
        F: Vec::new(),
        q: 0,
        p: 0,
        istage: 0,
        stages: 0,
        B: None,
        nfe: 0,
        cvals: Vec::new(),
        nfusedopvecs: 0,
        tshift: 0.0,
        tscale: 0.0,
        forcing: Vec::new(),
        nforcing: 0,
        stage_times: Vec::new(),
        stage_coefs: Vec::new(),
    });

    /* Copy the input parameters into ARKODE state */
    step_mem.f = Some(f);

    /* Attach step_mem structure and function pointers to ark_mem */
    ark_mem.step_init = Some(erkStep_Init);
    ark_mem.step_fullrhs = Some(erkStep_FullRHS);
    ark_mem.step = Some(erkStep_TakeStep);
    ark_mem.step_printallstats = Some(crate::arkode_erkstep_io::erkStep_PrintAllStats);
    ark_mem.step_writeparameters = Some(crate::arkode_erkstep_io::erkStep_WriteParameters);
    ark_mem.step_setusecompensatedsums = None;
    ark_mem.step_resize = Some(erkStep_Resize);
    ark_mem.step_free = Some(erkStep_Free);
    ark_mem.step_printmem = Some(erkStep_PrintMem);
    ark_mem.step_setoptions = Some(crate::arkode_erkstep_io::erkStep_SetOptions);
    ark_mem.step_setdefaults = Some(crate::arkode_erkstep_io::erkStep_SetDefaults);
    ark_mem.step_setrelaxfn = None; /* erkStep_SetRelaxFn: relaxation module pending */
    ark_mem.step_setorder = Some(crate::arkode_erkstep_io::erkStep_SetOrder);
    ark_mem.step_getnumrhsevals = Some(crate::arkode_erkstep_io::erkStep_GetNumRhsEvals);
    ark_mem.step_getestlocalerrors = Some(crate::arkode_erkstep_io::erkStep_GetEstLocalErrors);
    ark_mem.step_setforcing = Some(erkStep_SetInnerForcing);
    ark_mem.step_getstageindex = Some(crate::arkode_erkstep_io::erkStep_GetStageIndex);
    ark_mem.step_supports_adaptive = true;
    ark_mem.step_supports_relaxation = true;
    ark_mem.step_mem = Some(step_mem);

    /* Set default values for optional inputs */
    let retval = crate::arkode_erkstep_io::erkStep_SetDefaults(&mut ark_mem);
    debug_assert_eq!(retval, ARK_SUCCESS);

    /* NOTE: F, cvals and Xvecs will be allocated later on
    (based on the number of ERK stages) */

    /* Update the ARKODE workspace requirements -- UPDATE */
    ark_mem.liw += 41; /* fcn/data ptr, int, long int, sunindextype, sunbooleantype */
    ark_mem.lrw += 10;

    /* (counters, fused-op workspace and forcing data already zeroed) */

    /* Initialize main ARKODE infrastructure */
    let retval = arkInit(&mut ark_mem, t0, y0, FIRST_INIT);
    debug_assert_eq!(retval, ARK_SUCCESS);

    ark_mem
}

/*---------------------------------------------------------------
  ERKStepReInit:

  This routine re-initializes the ERKStep module to solve a new
  problem of the same size as was previously solved. Note all
  internal counters are set to 0 on re-initialization.
  ---------------------------------------------------------------*/
pub fn ERKStepReInit(ark_mem: &mut ARKodeMem, f: ARKRhsFn, t0: f64, y0: &NVector) -> i32 {
    /* access ARKodeERKStepMem structure */
    let mut step_mem = match erkStep_AccessStepMem(ark_mem, "ERKStepReInit") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* Check if ark_mem was allocated */
    if !ark_mem.MallocDone {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_NO_MALLOC,
            line!(),
            "ERKStepReInit",
            file!(),
            "Attempt to call before ARKODE initialized.",
        );
        return ARK_NO_MALLOC;
    }

    /* Copy the input parameters into ARKODE state */
    step_mem.f = Some(f);
    step_mem.nfe = 0;
    ark_mem.step_mem = Some(step_mem);

    /* Initialize main ARKODE infrastructure */
    let retval = arkInit(ark_mem, t0, y0, FIRST_INIT);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            retval,
            line!(),
            "ERKStepReInit",
            file!(),
            "Unable to initialize main ARKODE infrastructure",
        );
        return retval;
    }

    ARK_SUCCESS
}

/*===============================================================
  Interface routines supplied to ARKODE
  ===============================================================*/

/*---------------------------------------------------------------
  erkStep_Resize:

  This routine resizes the memory within the ERKStep module.
  ---------------------------------------------------------------*/
pub fn erkStep_Resize(
    ark_mem: &mut ARKodeMem,
    y0: &NVector,
    _hscale: f64,
    _t0: f64,
    resize: Option<ARKVecResizeFn>,
    resize_data: &mut UserData,
) -> i32 {
    use crate::arkode_impl::ARK_MEM_FAIL;

    /* access ARKodeERKStepMem structure */
    let mut step_mem = match erkStep_AccessStepMem(ark_mem, "erkStep_Resize") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* Determine change in vector sizes */
    let lrw1 = y0.data.len() as i64;
    let liw1 = 1i64;
    let lrw_diff = lrw1 - ark_mem.lrw1;
    let liw_diff = liw1 - ark_mem.liw1;
    ark_mem.lrw1 = lrw1;
    ark_mem.liw1 = liw1;

    /* Resize the RHS vectors */
    for i in 0..step_mem.stages as usize {
        if !arkResizeVec(
            ark_mem,
            resize,
            resize_data,
            lrw_diff,
            liw_diff,
            y0,
            &mut step_mem.F[i],
        ) {
            ark_mem.step_mem = Some(step_mem);
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_FAIL,
                line!(),
                "erkStep_Resize",
                file!(),
                "Unable to resize vector",
            );
            return ARK_MEM_FAIL;
        }
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  erkStep_Free frees all ERKStep memory (the drops themselves are
  handled by ownership; this keeps the C workspace accounting).
  ---------------------------------------------------------------*/
pub fn erkStep_Free(ark_mem: &mut ARKodeMem) {
    /* conditional frees on non-NULL ERKStep module */
    if let Some(b) = ark_mem.step_mem.take() {
        if let Ok(mut step_mem) = b.downcast::<ARKodeERKStepMem>() {
            /* free the Butcher table */
            if let Some(bt) = step_mem.B.take() {
                let (mut bliw, mut blrw) = (0i64, 0i64);
                ARKodeButcherTable_Space(&bt, &mut bliw, &mut blrw);
                ark_mem.liw -= bliw;
                ark_mem.lrw -= blrw;
            }

            /* free the RHS vectors */
            {
                let stages = step_mem.stages;
                let (lrw1, liw1) = (ark_mem.lrw1, ark_mem.liw1);
                let (mut lrw, mut liw) = (ark_mem.lrw, ark_mem.liw);
                arkFreeVecArray(stages, &mut step_mem.F, lrw1, &mut lrw, liw1, &mut liw);
                ark_mem.lrw = lrw;
                ark_mem.liw = liw;
            }

            /* free the reusable arrays for fused vector interface */
            if !step_mem.cvals.is_empty() {
                step_mem.cvals = Vec::new();
                ark_mem.lrw -= step_mem.nfusedopvecs as i64;
            }
            /* (Xvecs pointer array: liw accounting only) */
            if step_mem.nfusedopvecs > 0 {
                ark_mem.liw -= step_mem.nfusedopvecs as i64;
            }
            step_mem.nfusedopvecs = 0;

            /* free work arrays for MRI forcing */
            if !step_mem.stage_times.is_empty() {
                step_mem.stage_times = Vec::new();
                ark_mem.lrw -= step_mem.stages as i64;
            }
            if !step_mem.stage_coefs.is_empty() {
                step_mem.stage_coefs = Vec::new();
                ark_mem.lrw -= step_mem.stages as i64;
            }

            /* free the time stepper module itself: drop */
            drop(step_mem);
        }
    }
}

/*---------------------------------------------------------------
  erkStep_PrintMem:

  This routine outputs the memory from the ERKStep structure to
  a specified file pointer (useful when debugging).
  ---------------------------------------------------------------*/
pub fn erkStep_PrintMem(ark_mem: &mut ARKodeMem, outfile: &mut dyn std::io::Write) {
    /* access ARKodeERKStepMem structure */
    let step_mem = match erkStep_AccessStepMem(ark_mem, "erkStep_PrintMem") {
        None => return,
        Some(sm) => sm,
    };

    /* output integer quantities */
    let _ = write!(outfile, "ERKStep: q = {}\n", step_mem.q);
    let _ = write!(outfile, "ERKStep: p = {}\n", step_mem.p);
    let _ = write!(outfile, "ERKStep: istage = {}\n", step_mem.istage);
    let _ = write!(outfile, "ERKStep: stages = {}\n", step_mem.stages);

    /* output long integer quantities */
    let _ = write!(outfile, "ERKStep: nfe = {}\n", step_mem.nfe);

    /* output sunrealtype quantities */
    let _ = write!(outfile, "ERKStep: Butcher table:\n");
    if let Some(bt) = step_mem.B.as_ref() {
        ARKodeButcherTable_Write(bt, outfile);
    }

    ark_mem.step_mem = Some(step_mem);
}

/*---------------------------------------------------------------
  erkStep_Init:

  This routine is called just prior to performing internal time
  steps (after all user "set" routines have been called) from
  within arkInitialSetup.

  With initialization types FIRST_INIT this routine:
  - sets/checks the ARK Butcher tables to be used
  - allocates any memory that depends on the number of
    stages, method order, or solver options
  - sets the call_fullrhs flag

  With other initialization types, this routine does nothing.
  ---------------------------------------------------------------*/
pub fn erkStep_Init(ark_mem: &mut ARKodeMem, init_type: i32) -> i32 {
    /* access ARKodeERKStepMem structure */
    let mut step_mem = match erkStep_AccessStepMem(ark_mem, "erkStep_Init") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };
    let ret = erkStep_Init_inner(ark_mem, &mut step_mem, init_type);
    ark_mem.step_mem = Some(step_mem);
    ret
}

fn erkStep_Init_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut ARKodeERKStepMem,
    init_type: i32,
) -> i32 {
    /* immediately return if resize or reset */
    if init_type == RESIZE_INIT || init_type == RESET_INIT {
        return ARK_SUCCESS;
    }

    /* enforce use of arkEwtSmallReal if using a fixed step size,
    an internal error weight function, and not performing accumulated
    temporal error estimation */
    let mut reset_efun = true;
    if !ark_mem.fixedstep {
        reset_efun = false;
    }
    if ark_mem.user_efun {
        reset_efun = false;
    }
    if ark_mem.AccumErrorType != ARK_ACCUMERROR_NONE {
        reset_efun = false;
    }
    if reset_efun {
        /* internal SmallReal efun (C: efun = arkEwtSetSmallReal with
        e_data = ark_mem); the dispatch helper calls a non-user Some(efun)
        directly */
        ark_mem.user_efun = false;
        ark_mem.efun = Some(erk_ewt_small_real);
        ark_mem.e_data = None;
    }

    /* Create Butcher table (if not already set) */
    let retval = erkStep_SetButcherTable_inner(ark_mem, step_mem);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "erkStep_Init",
            file!(),
            "Could not create Butcher table",
        );
        return ARK_ILL_INPUT;
    }

    /* Check that Butcher table are OK */
    let retval = erkStep_CheckButcherTable_inner(ark_mem, step_mem);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "erkStep_Init",
            file!(),
            "Error in Butcher table",
        );
        return ARK_ILL_INPUT;
    }

    /* Retrieve/store method and embedding orders now that table is
    finalized */
    step_mem.q = step_mem.B.as_ref().unwrap().q;
    step_mem.p = step_mem.B.as_ref().unwrap().p;
    {
        let hadapt_mem = ark_mem.hadapt_mem.as_mut().unwrap();
        hadapt_mem.q = step_mem.q;
        hadapt_mem.p = step_mem.p;
    }

    /* Ensure that if adaptivity or error accumulation is enabled, then
    method includes embedding coefficients */
    if (!ark_mem.fixedstep || (ark_mem.AccumErrorType != ARK_ACCUMERROR_NONE))
        && (step_mem.p == 0)
    {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "erkStep_Init",
            file!(),
            "Temporal error estimation cannot be performed without embedding coefficients",
        );
        return ARK_ILL_INPUT;
    }

    /* Allocate RHS vector memory, update storage requirements */
    /*   Allocate F[0] ... F[stages-1] if needed */
    {
        let stages = step_mem.stages;
        let tmpl_len = ark_mem.ewt.data.len();
        let (lrw1, liw1) = (ark_mem.lrw1, ark_mem.liw1);
        let (mut lrw, mut liw) = (ark_mem.lrw, ark_mem.liw);
        arkAllocVecArray(stages, tmpl_len, &mut step_mem.F, lrw1, &mut lrw, liw1, &mut liw);
        ark_mem.lrw = lrw;
        ark_mem.liw = liw;
    }

    /* Allocate reusable arrays for fused vector interface */
    step_mem.nfusedopvecs = 2 * step_mem.stages + 2 + step_mem.nforcing;
    if step_mem.cvals.is_empty() {
        step_mem.cvals = vec![0.0; step_mem.nfusedopvecs as usize];
        ark_mem.lrw += step_mem.nfusedopvecs as i64;
        /* (Xvecs pointer array: operands are assembled at the call sites;
        keep the C liw accounting — C allocates Xvecs alongside cvals and
        only then adds liw, so a ReInit does not re-count it) */
        ark_mem.liw += step_mem.nfusedopvecs as i64;
    }

    /* Allocate workspace for MRI forcing -- need to allocate here as the
    number of stages may not be set before this point */
    if step_mem.stage_times.is_empty() {
        step_mem.stage_times = vec![0.0; step_mem.stages as usize];
        ark_mem.lrw += step_mem.stages as i64;
    }
    if step_mem.stage_coefs.is_empty() {
        step_mem.stage_coefs = vec![0.0; step_mem.stages as usize];
        ark_mem.lrw += step_mem.stages as i64;
    }

    /* Override the interpolant degree (if needed), used in
    arkInitialSetup */
    if step_mem.q > 1 && ark_mem.interp_degree > (step_mem.q - 1) {
        /* Limit max degree to at most one less than the method global
        order */
        ark_mem.interp_degree = step_mem.q - 1;
    } else if step_mem.q == 1 && ark_mem.interp_degree > 1 {
        /* Allow for linear interpolant with first order methods to ensure
        solution values are returned at the time interval end points */
        ark_mem.interp_degree = 1;
    }

    /* set appropriate TakeStep routine based on problem configuration
    (do_adjoint can only be set by the deferred adjoint machinery) */
    ark_mem.step = Some(erkStep_TakeStep);

    /* Signal to shared arkode module that full RHS evaluations are
    required */
    ark_mem.call_fullrhs = true;

    ARK_SUCCESS
}

/// C arkEwtSetSmallReal installed as an efun (see erkStep_Init).
fn erk_ewt_small_real(_ycur: &NVector, weight: &mut NVector, _e_data: &mut UserData) -> i32 {
    crate::nvector_serial::N_VConst(crate::sundials_types::SUN_SMALL_REAL, weight);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  erkStep_FullRHS:

  This is just a wrapper to call the user-supplied RHS function,
  f(t,y). See the C source for the full description of the three
  ARK_FULLRHS_* modes.
  ---------------------------------------------------------------*/
pub fn erkStep_FullRHS(ark_mem: &mut ARKodeMem, t: f64, y: &NVector, f: &mut NVector, mode: i32) -> i32 {
    /* access ARKodeERKStepMem structure */
    let mut step_mem = match erkStep_AccessStepMem(ark_mem, "erkStep_FullRHS") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };
    let ret = erkStep_FullRHS_inner(ark_mem, &mut step_mem, t, y, f, mode);
    ark_mem.step_mem = Some(step_mem);
    ret
}

fn erkStep_FullRHS_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut ARKodeERKStepMem,
    t: f64,
    y: &NVector,
    f: &mut NVector,
    mode: i32,
) -> i32 {
    /* perform RHS functions contingent on 'mode' argument */
    match mode {
        m if m == ARK_FULLRHS_START => {
            /* compute the RHS if needed */
            if !ark_mem.fn_is_current {
                /* call the user-supplied pre-RHS function (if supplied) */
                if let Some(pre_rhs_fn) = ark_mem.PreRhsFn {
                    let retval = pre_rhs_fn(t, y, &mut ark_mem.user_data);
                    if retval != 0 {
                        return ARK_PRERHSFN_FAIL;
                    }
                }

                /* call f */
                let ff = step_mem.f.unwrap();
                let retval = ff(t, y, &mut step_mem.F[0], &mut ark_mem.user_data);
                step_mem.nfe += 1;
                if retval != 0 {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_RHSFUNC_FAIL,
                        line!(),
                        "erkStep_FullRHS",
                        file!(),
                        &format!(
                            "At t = {}, the right-hand side routine failed in an unrecoverable manner.",
                            fmt_g(t, 0, 15)
                        ),
                    );
                    return ARK_RHSFUNC_FAIL;
                }
            }

            /* copy RHS into output */
            N_VScale(ONE, &step_mem.F[0], f);

            /* apply external polynomial forcing */
            if step_mem.nforcing > 0 {
                let vals = erkStep_ApplyForcing_coeffs(step_mem, &[t], &[ONE], 1);
                erk_accumulate_forcing(step_mem, &vals, f);
            }
        }

        m if m == ARK_FULLRHS_END => {
            /* determine if RHS function needs to be recomputed */
            if !ark_mem.fn_is_current {
                let mut recompute_rhs =
                    !ARKodeButcherTable_IsStifflyAccurate(step_mem.B.as_ref().unwrap());

                /* First Same As Last methods are not FSAL when relaxation is
                enabled */
                if ark_mem.relax_enabled {
                    recompute_rhs = true;
                }

                /* base RHS call on recomputeRHS argument */
                if recompute_rhs {
                    /* call the user-supplied pre-RHS function (if supplied) */
                    if let Some(pre_rhs_fn) = ark_mem.PreRhsFn {
                        let retval = pre_rhs_fn(t, y, &mut ark_mem.user_data);
                        if retval != 0 {
                            return ARK_PRERHSFN_FAIL;
                        }
                    }

                    /* call f */
                    let ff = step_mem.f.unwrap();
                    let retval = ff(t, y, &mut step_mem.F[0], &mut ark_mem.user_data);
                    step_mem.nfe += 1;
                    if retval != 0 {
                        arkProcessError(
                            Some(ark_mem),
                            ARK_RHSFUNC_FAIL,
                            line!(),
                            "erkStep_FullRHS",
                            file!(),
                            &format!(
                                "At t = {}, the right-hand side routine failed in an unrecoverable manner.",
                                fmt_g(t, 0, 15)
                            ),
                        );
                        return ARK_RHSFUNC_FAIL;
                    }
                } else {
                    /* N_VScale(ONE, F[stages-1], F[0]) */
                    let stages = step_mem.stages as usize;
                    let (first, rest) = step_mem.F.split_at_mut(1);
                    first[0].data.copy_from_slice(&rest[stages - 2].data);
                }

                /* copy RHS vector into output */
                N_VScale(ONE, &step_mem.F[0], f);

                /* apply external polynomial forcing */
                if step_mem.nforcing > 0 {
                    let vals = erkStep_ApplyForcing_coeffs(step_mem, &[t], &[ONE], 1);
                    erk_accumulate_forcing(step_mem, &vals, f);
                }
            }
        }

        m if m == ARK_FULLRHS_OTHER => {
            /* call the user-supplied pre-RHS function (if supplied) */
            if let Some(pre_rhs_fn) = ark_mem.PreRhsFn {
                let retval = pre_rhs_fn(t, y, &mut ark_mem.user_data);
                if retval != 0 {
                    return ARK_PRERHSFN_FAIL;
                }
            }

            /* call f */
            let ff = step_mem.f.unwrap();
            let retval = ff(t, y, f, &mut ark_mem.user_data);
            step_mem.nfe += 1;
            if retval != 0 {
                arkProcessError(
                    Some(ark_mem),
                    ARK_RHSFUNC_FAIL,
                    line!(),
                    "erkStep_FullRHS",
                    file!(),
                    &format!(
                        "At t = {}, the right-hand side routine failed in an unrecoverable manner.",
                        fmt_g(t, 0, 15)
                    ),
                );
                return ARK_RHSFUNC_FAIL;
            }
            /* apply external polynomial forcing */
            if step_mem.nforcing > 0 {
                let vals = erkStep_ApplyForcing_coeffs(step_mem, &[t], &[ONE], 1);
                erk_accumulate_forcing(step_mem, &vals, f);
            }
        }

        _ => {
            /* return with RHS failure if unknown mode is passed */
            arkProcessError(
                Some(ark_mem),
                ARK_RHSFUNC_FAIL,
                line!(),
                "erkStep_FullRHS",
                file!(),
                "Unknown full RHS mode",
            );
            return ARK_RHSFUNC_FAIL;
        }
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  erkStep_TakeStep:

  This routine serves the primary purpose of the ERKStep module:
  it performs a single ERK step (with embedding, if possible).
  See the C source for the output/return-value conventions.
  ---------------------------------------------------------------*/
pub fn erkStep_TakeStep(ark_mem: &mut ARKodeMem, dsmPtr: &mut f64, nflagPtr: &mut i32) -> i32 {
    /* initialize algebraic solver convergence flag to success */
    *nflagPtr = ARK_SUCCESS;

    /* access ARKodeERKStepMem structure; initialize the current stage
    index and determine the FSAL property, then release step_mem for the
    full-RHS re-entry below */
    let fsal;
    {
        let mut step_mem = match erkStep_AccessStepMem(ark_mem, "erkStep_TakeStep") {
            None => return ARK_MEM_NULL,
            Some(sm) => sm,
        };
        fsal = ARKodeButcherTable_IsStifflyAccurate(step_mem.B.as_ref().unwrap());
        step_mem.istage = 0;
        ark_mem.step_mem = Some(step_mem);
    }

    /* Call the full RHS if needed. If this is the first step then we may
    need to evaluate or copy the RHS values from an earlier evaluation
    (e.g., to compute h0). For subsequent steps treat this RHS evaluation
    as an evaluation at the end of the just completed step to potentially
    reuse (FSAL methods) RHS evaluations from the end of the last step. */
    if !ark_mem.fn_is_current {
        let mode = if ark_mem.initsetup {
            ARK_FULLRHS_START
        } else {
            ARK_FULLRHS_END
        };
        let retval = crate::arkode_impl::ark_step_fullrhs_yn_fn(ark_mem, ark_mem.tn, mode);
        if retval != 0 {
            return ARK_RHSFUNC_FAIL;
        }
        ark_mem.fn_is_current = true;
    }

    /* checkpoint the step state (if necessary) */
    if ark_mem.checkpoint_scheme.is_some() {
        let retval = erk_checkpoint(ark_mem, 0, ark_mem.tn, CheckpointVec::Yn);
        if retval != ARK_SUCCESS {
            return retval;
        }
    }

    /* re-take step_mem for the stage loop and solution computation */
    let mut step_mem = match erkStep_AccessStepMem(ark_mem, "erkStep_TakeStep") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };
    let ret = erkStep_TakeStep_inner(ark_mem, &mut step_mem, dsmPtr, fsal);
    ark_mem.step_mem = Some(step_mem);
    ret
}

fn erkStep_TakeStep_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut ARKodeERKStepMem,
    dsmPtr: &mut f64,
    fsal: bool,
) -> i32 {
    /* Loop over internal stages to the step; since the method is explicit
    the first stage RHS is just the full RHS from the start of the step */
    for is in 1..step_mem.stages as usize {
        /* Set current stage time and index */
        ark_mem.tcur =
            ark_mem.tn + step_mem.B.as_ref().unwrap().c[is] * ark_mem.h;
        step_mem.istage = is as i32;

        /* Set ycur to current stage solution (C: fused linear
        combination over [h*A[is][js]*F[js] ..., 1*yn, forcing...]) */
        {
            let b = step_mem.B.as_ref().unwrap();
            let mut cv: Vec<f64> = Vec::with_capacity(is + 2 + step_mem.nforcing as usize);
            for js in 0..is {
                cv.push(ark_mem.h * b.A[is][js]);
            }
            cv.push(ONE);

            /* apply external polynomial forcing */
            if step_mem.nforcing > 0 {
                for js in 0..is {
                    step_mem.stage_times[js] = ark_mem.tn + b.c[js] * ark_mem.h;
                    step_mem.stage_coefs[js] = ark_mem.h * b.A[is][js];
                }
                let fvals = erkStep_ApplyForcing_coeffs(
                    step_mem,
                    &step_mem.stage_times[..is],
                    &step_mem.stage_coefs[..is],
                    is,
                );
                cv.extend_from_slice(&fvals);
            }

            let mut xrefs: Vec<&NVector> = Vec::with_capacity(cv.len());
            for js in 0..is {
                xrefs.push(&step_mem.F[js]);
            }
            xrefs.push(&ark_mem.yn);
            for k in 0..step_mem.nforcing as usize {
                xrefs.push(&step_mem.forcing[k]);
            }

            let retval =
                N_VLinearCombination(cv.len() as i32, &cv, &xrefs, &mut ark_mem.ycur);
            if retval != 0 {
                return ARK_VECTOROP_ERR;
            }
        }

        /* apply user-supplied stage postprocessing function (if supplied)
        unless this is the last stage of a FSAL method, then apply the
        user-supplied step postprocessing function (if supplied) */
        let last_fsal_stage = is as i32 == step_mem.stages - 1 && fsal;
        if last_fsal_stage && ark_mem.PostProcessStepFn.is_some() {
            if let Some(post) = ark_mem.PostProcessStepFn {
                let retval = post(ark_mem.tcur, &ark_mem.ycur, &mut ark_mem.user_data);
                if retval != 0 {
                    return ARK_POSTPROCESS_STEP_FAIL;
                }
            }
        } else if let Some(post) = ark_mem.PostProcessStageFn {
            let retval = post(ark_mem.tcur, &ark_mem.ycur, &mut ark_mem.user_data);
            if retval != 0 {
                return ARK_POSTPROCESS_STAGE_FAIL;
            }
        }

        /* call the user-supplied pre-RHS function (if supplied) */
        if let Some(pre_rhs_fn) = ark_mem.PreRhsFn {
            let retval = pre_rhs_fn(ark_mem.tcur, &ark_mem.ycur, &mut ark_mem.user_data);
            if retval != 0 {
                return ARK_PRERHSFN_FAIL;
            }
        }

        /* compute updated RHS */
        let ff = step_mem.f.unwrap();
        let retval = ff(
            ark_mem.tcur,
            &ark_mem.ycur,
            &mut step_mem.F[is],
            &mut ark_mem.user_data,
        );
        step_mem.nfe += 1;

        if retval < 0 {
            return ARK_RHSFUNC_FAIL;
        }
        if retval > 0 {
            return ARK_UNREC_RHSFUNC_ERR;
        }

        /* checkpoint stage for adjoint (if necessary) */
        if ark_mem.checkpoint_scheme.is_some() {
            let retval = erk_checkpoint(ark_mem, is as i64, ark_mem.tcur, CheckpointVec::Ycur);
            if retval != ARK_SUCCESS {
                return retval;
            }
        }
    } /* loop over stages */

    /* compute time-evolved solution (in ark_ycur), error estimate (in
    dsm) */
    ark_mem.tcur = ark_mem.tn + ark_mem.h;

    let retval = erkStep_ComputeSolutions_inner(ark_mem, step_mem, dsmPtr);
    if retval < 0 {
        return retval;
    }

    /* checkpoint the step solution (if necessary) */
    if ark_mem.checkpoint_scheme.is_some() {
        let retval = erk_checkpoint(
            ark_mem,
            step_mem.B.as_ref().unwrap().stages as i64,
            ark_mem.tn + ark_mem.h,
            CheckpointVec::Ycur,
        );
        if retval != ARK_SUCCESS {
            return retval;
        }
    }

    ARK_SUCCESS
}

/// Which ark_mem vector a checkpoint insertion snapshots.
enum CheckpointVec {
    Yn,
    Ycur,
}

/// The C checkpoint blocks in erkStep_TakeStep (NeedsSaving followed
/// by InsertVector, with the shared error reporting).
fn erk_checkpoint(ark_mem: &mut ARKodeMem, stage: i64, t: f64, which: CheckpointVec) -> i32 {
    use crate::arkode_impl::ARK_ADJ_CHECKPOINT_FAIL;
    use crate::sundials_adjointcheckpointscheme::{
        SUNAdjointCheckpointScheme_InsertVector, SUNAdjointCheckpointScheme_NeedsSaving,
    };

    let mut scheme = ark_mem.checkpoint_scheme.take().unwrap();

    let mut do_save = false;
    let errcode = SUNAdjointCheckpointScheme_NeedsSaving(
        &mut scheme,
        ark_mem.checkpoint_step_idx,
        stage,
        t,
        &mut do_save,
    );
    if errcode != 0 {
        ark_mem.checkpoint_scheme = Some(scheme);
        arkProcessError(
            Some(ark_mem),
            ARK_ADJ_CHECKPOINT_FAIL,
            line!(),
            "erkStep_TakeStep",
            file!(),
            &format!("SUNAdjointCheckpointScheme_NeedsSaving returned {}", errcode),
        );
        return ARK_ADJ_CHECKPOINT_FAIL;
    }

    if do_save {
        let y = match which {
            CheckpointVec::Yn => &ark_mem.yn,
            CheckpointVec::Ycur => &ark_mem.ycur,
        };
        let errcode = SUNAdjointCheckpointScheme_InsertVector(
            &mut scheme,
            ark_mem.checkpoint_step_idx,
            stage,
            t,
            y,
        );
        if errcode != 0 {
            ark_mem.checkpoint_scheme = Some(scheme);
            arkProcessError(
                Some(ark_mem),
                ARK_ADJ_CHECKPOINT_FAIL,
                line!(),
                "erkStep_TakeStep",
                file!(),
                &format!("SUNAdjointCheckpointScheme_InsertVector returned {}", errcode),
            );
            return ARK_ADJ_CHECKPOINT_FAIL;
        }
    }

    ark_mem.checkpoint_scheme = Some(scheme);
    ARK_SUCCESS
}

/*===============================================================
  Internal utility routines
  ===============================================================*/

/*---------------------------------------------------------------
  erkStep_SetButcherTable

  This routine determines the ERK method to use, based on the
  desired accuracy.
  ---------------------------------------------------------------*/
pub fn erkStep_SetButcherTable(ark_mem: &mut ARKodeMem) -> i32 {
    let mut step_mem = match erkStep_AccessStepMem(ark_mem, "erkStep_SetButcherTable") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };
    let ret = erkStep_SetButcherTable_inner(ark_mem, &mut step_mem);
    ark_mem.step_mem = Some(step_mem);
    ret
}

fn erkStep_SetButcherTable_inner(ark_mem: &mut ARKodeMem, step_mem: &mut ARKodeERKStepMem) -> i32 {
    /* if table has already been specified, just return */
    if step_mem.B.is_some() {
        return ARK_SUCCESS;
    }

    /* select method based on order */
    let etable = match step_mem.q {
        1 => ERKSTEP_DEFAULT_1,
        2 => ERKSTEP_DEFAULT_2,
        3 => ERKSTEP_DEFAULT_3,
        4 => ERKSTEP_DEFAULT_4,
        5 => ERKSTEP_DEFAULT_5,
        6 => ERKSTEP_DEFAULT_6,
        7 => ERKSTEP_DEFAULT_7,
        8 => ERKSTEP_DEFAULT_8,
        9 => ERKSTEP_DEFAULT_9,
        _ => {
            /* no available method, set default */
            arkProcessError(
                Some(ark_mem),
                ARK_WARNING,
                line!(),
                "erkStep_SetButcherTable",
                file!(),
                "No explicit method at requested order, using q=9.",
            );
            ERKSTEP_DEFAULT_9
        }
    };

    step_mem.B = ARKodeButcherTable_LoadERK(etable);

    /* note Butcher table space requirements */
    let (mut bliw, mut blrw) = (0i64, 0i64);
    ARKodeButcherTable_Space(step_mem.B.as_ref().unwrap(), &mut bliw, &mut blrw);
    ark_mem.liw += bliw;
    ark_mem.lrw += blrw;

    /* set [redundant] stored values for stage numbers and method orders */
    if let Some(b) = step_mem.B.as_ref() {
        step_mem.stages = b.stages;
        step_mem.q = b.q;
        step_mem.p = b.p;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  erkStep_CheckButcherTable

  This routine runs through the explicit Butcher table to ensure
  that it meets all necessary requirements, including:
    strictly lower-triangular (ERK)
    method order q > 0 (all)
    embedding order q > 0 (all -- if adaptive time-stepping enabled)
    stages > 0 (all)

  Returns ARK_SUCCESS if tables pass, ARK_INVALID_TABLE otherwise.
  ---------------------------------------------------------------*/
pub fn erkStep_CheckButcherTable(ark_mem: &mut ARKodeMem) -> i32 {
    let mut step_mem = match erkStep_AccessStepMem(ark_mem, "erkStep_CheckButcherTable") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };
    let ret = erkStep_CheckButcherTable_inner(ark_mem, &mut step_mem);
    ark_mem.step_mem = Some(step_mem);
    ret
}

fn erkStep_CheckButcherTable_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut ARKodeERKStepMem,
) -> i32 {
    let tol: f64 = 1.0e-12;

    /* check that stages > 0 */
    if step_mem.stages < 1 {
        arkProcessError(
            Some(ark_mem),
            ARK_INVALID_TABLE,
            line!(),
            "erkStep_CheckButcherTable",
            file!(),
            "stages < 1!",
        );
        return ARK_INVALID_TABLE;
    }

    /* check that method order q > 0 */
    if step_mem.q < 1 {
        arkProcessError(
            Some(ark_mem),
            ARK_INVALID_TABLE,
            line!(),
            "erkStep_CheckButcherTable",
            file!(),
            "method order < 1!",
        );
        return ARK_INVALID_TABLE;
    }

    /* check that embedding order p > 0 */
    if (step_mem.p < 1) && !ark_mem.fixedstep {
        arkProcessError(
            Some(ark_mem),
            ARK_INVALID_TABLE,
            line!(),
            "erkStep_CheckButcherTable",
            file!(),
            "embedding order < 1!",
        );
        return ARK_INVALID_TABLE;
    }

    /* check that embedding exists */
    if (step_mem.p > 0) && !ark_mem.fixedstep {
        if step_mem.B.as_ref().unwrap().d.is_none() {
            arkProcessError(
                Some(ark_mem),
                ARK_INVALID_TABLE,
                line!(),
                "erkStep_CheckButcherTable",
                file!(),
                "no embedding!",
            );
            return ARK_INVALID_TABLE;
        }
    }

    /* check that ERK table is strictly lower triangular */
    let mut okay = true;
    {
        let b = step_mem.B.as_ref().unwrap();
        for i in 0..step_mem.stages as usize {
            for j in i..step_mem.stages as usize {
                if SUNRabs(b.A[i][j]) > tol {
                    okay = false;
                }
            }
        }
    }
    if !okay {
        arkProcessError(
            Some(ark_mem),
            ARK_INVALID_TABLE,
            line!(),
            "erkStep_CheckButcherTable",
            file!(),
            "Ae Butcher table is implicit!",
        );
        return ARK_INVALID_TABLE;
    }

    /* check if all b values are positive for relaxation */
    if ark_mem.relax_enabled {
        if step_mem.q < 2 {
            arkProcessError(
                Some(ark_mem),
                ARK_INVALID_TABLE,
                line!(),
                "erkStep_CheckButcherTable",
                file!(),
                "The Butcher table must be at least second order!",
            );
            return ARK_INVALID_TABLE;
        }

        let b = step_mem.B.as_ref().unwrap();
        for i in 0..step_mem.stages as usize {
            if b.b[i] < ZERO {
                arkProcessError(
                    Some(ark_mem),
                    ARK_INVALID_TABLE,
                    line!(),
                    "erkStep_CheckButcherTable",
                    file!(),
                    "The Butcher table has a negative b value!",
                );
                return ARK_INVALID_TABLE;
            }
        }
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  erkStep_ComputeSolutions

  This routine calculates the final RK solution using the existing
  data.  This solution is placed directly in ark_ycur.  This
  routine also computes the error estimate ||y-ytilde||_WRMS, where
  ytilde is the embedded solution, and the norm weights come from
  ark_ewt.  This norm value is returned.  The vector form of this
  estimated error (y-ytilde) is stored in ark_tempv1, in case the
  calling routine wishes to examine the error locations.
  ---------------------------------------------------------------*/
pub(crate) fn erkStep_ComputeSolutions_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut ARKodeERKStepMem,
    dsmPtr: &mut f64,
) -> i32 {
    /* initialize output */
    *dsmPtr = ZERO;

    /* determine if method has fsal property */
    let fsal = ARKodeButcherTable_IsStifflyAccurate(step_mem.B.as_ref().unwrap());

    /* Compute time step solution. For FSAL methods, ycur already contains
    the new solution. */
    if !fsal {
        let b = step_mem.B.as_ref().unwrap();
        let stages = step_mem.stages as usize;
        let mut cv: Vec<f64> = Vec::with_capacity(stages + 1 + step_mem.nforcing as usize);
        for j in 0..stages {
            cv.push(ark_mem.h * b.b[j]);
        }
        cv.push(ONE);

        /* apply external polynomial forcing */
        if step_mem.nforcing > 0 {
            for j in 0..stages {
                step_mem.stage_times[j] = ark_mem.tn + b.c[j] * ark_mem.h;
                step_mem.stage_coefs[j] = ark_mem.h * b.b[j];
            }
            let fvals = erkStep_ApplyForcing_coeffs(
                step_mem,
                &step_mem.stage_times[..stages],
                &step_mem.stage_coefs[..stages],
                stages,
            );
            cv.extend_from_slice(&fvals);
        }

        let mut xrefs: Vec<&NVector> = Vec::with_capacity(cv.len());
        for j in 0..stages {
            xrefs.push(&step_mem.F[j]);
        }
        xrefs.push(&ark_mem.yn);
        for k in 0..step_mem.nforcing as usize {
            xrefs.push(&step_mem.forcing[k]);
        }

        let retval = N_VLinearCombination(cv.len() as i32, &cv, &xrefs, &mut ark_mem.ycur);
        if retval != 0 {
            return ARK_VECTOROP_ERR;
        }

        /* apply user-supplied step postprocessing function (if supplied) */
        if let Some(post) = ark_mem.PostProcessStepFn {
            let retval = post(ark_mem.tcur, &ark_mem.ycur, &mut ark_mem.user_data);
            if retval != 0 {
                return ARK_POSTPROCESS_STEP_FAIL;
            }
        }
    }

    /* Compute yerr (if step adaptivity or error accumulation enabled) */
    if !ark_mem.fixedstep || (ark_mem.AccumErrorType != ARK_ACCUMERROR_NONE) {
        let b = step_mem.B.as_ref().unwrap();
        let stages = step_mem.stages as usize;
        let d = b.d.as_ref().unwrap();
        let mut cv: Vec<f64> = Vec::with_capacity(stages + step_mem.nforcing as usize);
        for j in 0..stages {
            cv.push(ark_mem.h * (b.b[j] - d[j]));
        }

        /* apply external polynomial forcing */
        if step_mem.nforcing > 0 {
            for j in 0..stages {
                step_mem.stage_times[j] = ark_mem.tn + b.c[j] * ark_mem.h;
                step_mem.stage_coefs[j] = ark_mem.h * (b.b[j] - d[j]);
            }
            let fvals = erkStep_ApplyForcing_coeffs(
                step_mem,
                &step_mem.stage_times[..stages],
                &step_mem.stage_coefs[..stages],
                stages,
            );
            cv.extend_from_slice(&fvals);
        }

        let mut xrefs: Vec<&NVector> = Vec::with_capacity(cv.len());
        for j in 0..stages {
            xrefs.push(&step_mem.F[j]);
        }
        for k in 0..step_mem.nforcing as usize {
            xrefs.push(&step_mem.forcing[k]);
        }

        let retval = N_VLinearCombination(cv.len() as i32, &cv, &xrefs, &mut ark_mem.tempv1);
        if retval != 0 {
            return ARK_VECTOROP_ERR;
        }

        /* fill error norm */
        *dsmPtr = N_VWrmsNorm(&ark_mem.tempv1, &ark_mem.ewt);
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  erkStep_ApplyForcing_coeffs

  The coefficient half of C's erkStep_ApplyForcing: computes the
  polynomial-forcing coefficients (C fills cvals[offset..]; the
  vectors are appended by the call sites).
  ---------------------------------------------------------------*/
fn erkStep_ApplyForcing_coeffs(
    step_mem: &ARKodeERKStepMem,
    stage_times: &[f64],
    stage_coefs: &[f64],
    jmax: usize,
) -> Vec<f64> {
    let nforcing = step_mem.nforcing as usize;
    let mut vals = vec![ZERO; nforcing];

    for j in 0..jmax {
        let tau = (stage_times[j] - step_mem.tshift) / step_mem.tscale;
        let mut taui = ONE;

        for k in 0..nforcing {
            vals[k] += stage_coefs[j] * taui;
            taui *= tau;
        }
    }

    vals
}

/// z += sum_k vals[k] * forcing[k] — the z == X[0], c[0] == 1 branch
/// of the C N_VLinearCombination kernel used when applying forcing
/// to a full-RHS output.
fn erk_accumulate_forcing(step_mem: &ARKodeERKStepMem, vals: &[f64], f: &mut NVector) {
    for (k, val) in vals.iter().enumerate() {
        for e in 0..f.data.len() {
            f.data[e] += val * step_mem.forcing[k].data[e];
        }
    }
}

/*---------------------------------------------------------------
  erkStep_SetInnerForcing

  Sets an array of coefficient vectors for a time-dependent
  external polynomial forcing term in the ODE RHS. Primarily for
  use with MRIStep; the C code stores the caller's vector-array
  pointer, the Rust port stores owned copies (MRIStep re-sets the
  forcing before each fast integration).
  ---------------------------------------------------------------*/
pub fn erkStep_SetInnerForcing(
    ark_mem: &mut ARKodeMem,
    tshift: f64,
    tscale: f64,
    forcing: &[NVector],
    nvecs: i32,
) -> i32 {
    /* access ARKodeERKStepMem structure */
    let mut step_mem = match erkStep_AccessStepMem(ark_mem, "erkStep_SetInnerForcing") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    if nvecs > 0 {
        /* store forcing inputs */
        step_mem.tshift = tshift;
        step_mem.tscale = tscale;
        step_mem.forcing = forcing.to_vec();
        step_mem.nforcing = nvecs;

        /* check if there are enough reusable arrays for fused operations
        (only applies after erkStep_Init has allocated them) */
        if !step_mem.cvals.is_empty() && (step_mem.nfusedopvecs - nvecs) < (step_mem.stages + 1)
        {
            /* free current work space */
            ark_mem.lrw -= step_mem.nfusedopvecs as i64;
            ark_mem.liw -= step_mem.nfusedopvecs as i64;

            /* allocate reusable arrays for fused vector operations */
            step_mem.nfusedopvecs = step_mem.stages + 1 + nvecs;
            step_mem.cvals = vec![0.0; step_mem.nfusedopvecs as usize];
            ark_mem.lrw += step_mem.nfusedopvecs as i64;
            ark_mem.liw += step_mem.nfusedopvecs as i64;
        }
    } else {
        /* disable forcing */
        step_mem.tshift = ZERO;
        step_mem.tscale = ONE;
        step_mem.forcing = Vec::new();
        step_mem.nforcing = 0;
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arkode::{ARKodeEvolve, ARKodeSStolerances};
    use crate::arkode_erkstep_io::{ERKStepGetTimestepperStats, ERKStepSetTableNum};
    use crate::arkode_impl::{ARK_NORMAL, ARK_ONE_STEP};
    use crate::arkode_io::ARKodeSetFixedStep;
    use crate::sundials_context::SUNContext_Create;

    /* ark_analytic problem: ydot = lambda*y + 1/(1+t^2) - lambda*atan(t),
    y(0) = 0, exact solution y(t) = atan(t); lambda = -100 */
    fn analytic_rhs(t: f64, y: &NVector, ydot: &mut NVector, _ud: &mut UserData) -> i32 {
        let lamda = -100.0;
        ydot.data[0] = lamda * y.data[0] + 1.0 / (1.0 + t * t) - lamda * t.atan();
        0
    }

    #[test]
    fn erkstep_integrates_ark_analytic() {
        let ctx = SUNContext_Create();
        let mut y = NVector::new(1);
        y.data[0] = 0.0;
        let mut ark_mem = ERKStepCreate(analytic_rhs, 0.0, &y, &ctx);
        assert_eq!(ARKodeSStolerances(&mut ark_mem, 1.0e-6, 1.0e-10), ARK_SUCCESS);

        let mut tret = 0.0;
        let istate = ARKodeEvolve(&mut ark_mem, 1.0, &mut y, &mut tret, ARK_NORMAL);
        assert_eq!(istate, ARK_SUCCESS);
        assert_eq!(tret, 1.0);
        let err = (y.data[0] - 1.0_f64.atan()).abs();
        assert!(err < 1.0e-5, "err = {:e}", err);

        /* statistics are consistent: default order 4 table is
        Sofroniou-Spaletta (5 stages, FSAL) */
        let (mut exps, mut accs, mut att, mut nfe, mut netf) = (0, 0, 0, 0, 0);
        assert_eq!(
            ERKStepGetTimestepperStats(
                &mut ark_mem, &mut exps, &mut accs, &mut att, &mut nfe, &mut netf
            ),
            ARK_SUCCESS
        );
        assert!(ark_mem.nst > 10, "nst = {}", ark_mem.nst);
        assert!(att >= ark_mem.nst);
        assert!(nfe > 4 * ark_mem.nst, "nfe = {}, nst = {}", nfe, ark_mem.nst);
        assert_eq!(accs, ark_mem.nst_attempts);

        /* continue integrating to t = 2 */
        let istate = ARKodeEvolve(&mut ark_mem, 2.0, &mut y, &mut tret, ARK_NORMAL);
        assert_eq!(istate, ARK_SUCCESS);
        let err = (y.data[0] - 2.0_f64.atan()).abs();
        assert!(err < 1.0e-5, "err = {:e}", err);
    }

    /* y' = 2t (y = t^2): second-order methods reproduce it exactly */
    fn t2_rhs(t: f64, _y: &NVector, ydot: &mut NVector, _ud: &mut UserData) -> i32 {
        ydot.data[0] = 2.0 * t;
        0
    }

    #[test]
    fn erkstep_fixed_step_heun_exact_on_quadratic() {
        use crate::arkode_butcher_erk::ARKODE_HEUN_EULER_2_1_2;
        let ctx = SUNContext_Create();
        let mut y = NVector::new(1);
        y.data[0] = 0.0;
        let mut ark_mem = ERKStepCreate(t2_rhs, 0.0, &y, &ctx);
        assert_eq!(
            ERKStepSetTableNum(&mut ark_mem, ARKODE_HEUN_EULER_2_1_2),
            ARK_SUCCESS
        );
        assert_eq!(ARKodeSetFixedStep(&mut ark_mem, 0.125), ARK_SUCCESS);

        let mut tret = 0.0;
        for _ in 0..8 {
            let istate = ARKodeEvolve(&mut ark_mem, 1.0, &mut y, &mut tret, ARK_ONE_STEP);
            assert_eq!(istate, ARK_SUCCESS);
        }
        assert!((tret - 1.0).abs() < 1e-14, "tret = {}", tret);
        assert!((y.data[0] - 1.0).abs() < 1e-14, "y = {}", y.data[0]);
        assert_eq!(ark_mem.nst, 8);
    }

    #[test]
    fn erkstep_set_table_name_and_rejection() {
        let ctx = SUNContext_Create();
        let y = NVector::new(1);
        let mut ark_mem = ERKStepCreate(t2_rhs, 0.0, &y, &ctx);
        assert_eq!(
            crate::arkode_erkstep_io::ERKStepSetTableName(&mut ark_mem, "ARKODE_HEUN_EULER_2_1_2"),
            ARK_SUCCESS
        );
        /* unknown name maps to ARKODE_ERK_NONE and is rejected */
        assert_eq!(
            crate::arkode_erkstep_io::ERKStepSetTableName(&mut ark_mem, "NOT_A_TABLE"),
            ARK_ILL_INPUT
        );
    }

    #[test]
    fn erkstep_reinit_resets_counters() {
        let ctx = SUNContext_Create();
        let mut y = NVector::new(1);
        let mut ark_mem = ERKStepCreate(analytic_rhs, 0.0, &y, &ctx);
        assert_eq!(ARKodeSStolerances(&mut ark_mem, 1.0e-6, 1.0e-10), ARK_SUCCESS);
        let mut tret = 0.0;
        assert_eq!(
            ARKodeEvolve(&mut ark_mem, 0.5, &mut y, &mut tret, ARK_NORMAL),
            ARK_SUCCESS
        );
        assert!(ark_mem.nst > 0);

        y.data[0] = 0.0;
        assert_eq!(ERKStepReInit(&mut ark_mem, analytic_rhs, 0.0, &y), ARK_SUCCESS);
        assert_eq!(ark_mem.nst, 0);
        assert_eq!(
            ARKodeEvolve(&mut ark_mem, 1.0, &mut y, &mut tret, ARK_NORMAL),
            ARK_SUCCESS
        );
        let err = (y.data[0] - 1.0_f64.atan()).abs();
        assert!(err < 1.0e-5, "err = {:e}", err);
    }
}
