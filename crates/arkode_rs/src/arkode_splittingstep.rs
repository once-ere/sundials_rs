/* -----------------------------------------------------------------
 * Translation of src/arkode/arkode_splittingstep.c (SUNDIALS 7.7.0).
 *
 * This is the implementation file for ARKODE's operator splitting
 * module.
 *
 * step_mem follows the crate's take/downcast/put-back Access pattern;
 * TakeStep/FullRHS additionally mem::take the ark_mem vectors they
 * pass to the inner SUNSteppers (ycur/tempv1) and restore them on
 * every return path.  C NULL checks on `steppers`/`y0` arguments have
 * no Rust equivalent (slices/references) and are noted at their
 * sites; indexing steppers[0..partitions] with partitions >
 * steppers.len() panics where C reads out of bounds.
 * -----------------------------------------------------------------*/

use crate::arkode::{arkCreate, arkInit};
use crate::arkode_impl::{
    arkProcessError, ARKodeMem, ARK_ILL_INPUT, ARK_INTERP_HERMITE, ARK_INTERP_LAGRANGE,
    ARK_MEM_FAIL, ARK_MEM_NULL, ARK_NO_MALLOC, ARK_RHSFUNC_FAIL, ARK_SUCCESS, ARK_SUNSTEPPER_ERR,
    ARK_WARNING, FIRST_INIT, RESET_INIT, RESIZE_INIT, MSG_ARK_NO_MALLOC, ZERO,
};
use crate::sundials_types::SUNOutputFormat;
use crate::arkode_io::ARKodeSetInterpolantType;
use crate::arkode_splittingstep_coefficients::{
    SplittingStepCoefficientsMem, SplittingStepCoefficients_Copy,
    SplittingStepCoefficients_Destroy, SplittingStepCoefficients_LieTrotter,
    SplittingStepCoefficients_LoadCoefficientsByName, SplittingStepCoefficients_ThirdOrderSuzuki,
    SplittingStepCoefficients_TripleJump, SplittingStepCoefficients_Write,
};
use crate::arkode_splittingstep_impl::ARKodeSplittingStepMem;
use crate::nvector_serial::{NVector, N_VScale};
use crate::sundials_errors::SUN_SUCCESS;
use crate::sundials_stepper::{
    SUNStepper, SUNStepper_Evolve, SUNStepper_FullRhs, SUNStepper_Reset,
    SUNStepper_SetStepDirection, SUNStepper_SetStopTime, SUN_FULLRHS_OTHER,
};
use crate::sundials_utils::fmt_g;

const ONE: f64 = 1.0;
/* value of SUN_TABLE_WIDTH used by the sunfprintf_* printers */
const SUN_TABLE_WIDTH: usize = 29;

/*------------------------------------------------------------------------------
  Shortcut routine to unpack step_mem structure from ark_mem. If missing it
  returns ARK_MEM_NULL.  (Take semantics: caller must put the box back.)
  ----------------------------------------------------------------------------*/
pub(crate) fn splittingStep_AccessStepMem(
    ark_mem: &mut ARKodeMem,
    fname: &str,
) -> Option<Box<ARKodeSplittingStepMem>> {
    let taken = ark_mem.step_mem.take();
    match taken {
        Some(b) => match b.downcast::<ARKodeSplittingStepMem>() {
            Ok(sm) => Some(sm),
            Err(other) => {
                ark_mem.step_mem = Some(other);
                arkProcessError(
                    Some(ark_mem),
                    ARK_MEM_NULL,
                    line!(),
                    fname,
                    file!(),
                    "Time step module memory is NULL.",
                );
                None
            }
        },
        None => {
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_NULL,
                line!(),
                fname,
                file!(),
                "Time step module memory is NULL.",
            );
            None
        }
    }
}

/*------------------------------------------------------------------------------
  This routine determines the splitting coefficients to use based on the desired
  accuracy.
  ----------------------------------------------------------------------------*/
fn splittingStep_SetCoefficients(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut ARKodeSplittingStepMem,
) -> i32 {
    if step_mem.coefficients.is_some() {
        return ARK_SUCCESS;
    }

    if step_mem.order <= 1 {
        /* Lie-Trotter is the default (order < 1) */
        step_mem.coefficients = SplittingStepCoefficients_LieTrotter(step_mem.partitions);
    } else if step_mem.order == 3 {
        step_mem.coefficients = SplittingStepCoefficients_ThirdOrderSuzuki(step_mem.partitions);
    } else if step_mem.order % 2 == 0 {
        /* Triple jump only works for even order */
        step_mem.coefficients =
            SplittingStepCoefficients_TripleJump(step_mem.partitions, step_mem.order);
    } else {
        /* Bump the order up to be even but with a warning */
        let new_order = step_mem.order + 1;
        arkProcessError(
            Some(ark_mem),
            ARK_WARNING,
            line!(),
            "splittingStep_SetCoefficients",
            file!(),
            &format!(
                "No splitting method at requested order, using q={}.",
                new_order
            ),
        );
        step_mem.coefficients =
            SplittingStepCoefficients_TripleJump(step_mem.partitions, new_order);
    }

    if step_mem.coefficients.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_FAIL,
            line!(),
            "splittingStep_SetCoefficients",
            file!(),
            "Failed to allocate splitting coefficients",
        );
        return ARK_MEM_FAIL;
    }

    ARK_SUCCESS
}

/*-----------------------------------------------------------------------------
  This routine is called just prior to performing internal time steps (after all
  user "set" routines have been called) from within arkInitialSetup.

  With initialization types FIRST_INIT this routine:
  - sets/checks the splitting coefficients to be used

  With other initialization types, this routine does nothing.
  ----------------------------------------------------------------------------*/
fn splittingStep_Init(ark_mem: &mut ARKodeMem, init_type: i32) -> i32 {
    let mut step_mem = match splittingStep_AccessStepMem(ark_mem, "splittingStep_Init") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };
    let retval = splittingStep_Init_inner(ark_mem, &mut step_mem, init_type);
    ark_mem.step_mem = Some(step_mem);
    retval
}

fn splittingStep_Init_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut ARKodeSplittingStepMem,
    init_type: i32,
) -> i32 {
    if ark_mem.interp_type == ARK_INTERP_HERMITE {
        for i in 0..step_mem.partitions as usize {
            if step_mem.steppers[i].ops.fullrhs.is_none() {
                arkProcessError(
                    Some(ark_mem),
                    ARK_ILL_INPUT,
                    line!(),
                    "splittingStep_Init",
                    file!(),
                    &format!(
                        "steppers[{}] must implement SUNStepper_FullRhs when using Hermite interpolation",
                        i
                    ),
                );
                return ARK_ILL_INPUT;
            }
        }
    }

    /* inform arkode to ensure that ycur==yn upon entry to TakeStep function */
    ark_mem.ensure_ycur = true;

    /* immediately return if resize or reset */
    if init_type == RESIZE_INIT || init_type == RESET_INIT {
        return ARK_SUCCESS;
    }

    /* assume fixed step size */
    if !ark_mem.fixedstep {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "splittingStep_Init",
            file!(),
            "Adaptive outer time stepping is not currently supported",
        );
        return ARK_ILL_INPUT;
    }

    let retval = splittingStep_SetCoefficients(ark_mem, step_mem);
    if retval != ARK_SUCCESS {
        return retval;
    }

    ark_mem.interp_degree = std::cmp::max(
        1,
        std::cmp::min(
            step_mem.coefficients.as_ref().unwrap().order - 1,
            ark_mem.interp_degree,
        ),
    );

    ARK_SUCCESS
}

/*------------------------------------------------------------------------------
  This is just a wrapper to call the user-supplied RHS function,
  f^1(t,y) + f^2(t,y) + ... + f^P(t,y).

  This will be called in one of three 'modes':

     ARK_FULLRHS_START -> called at the beginning of a simulation i.e., at
                          (tn, yn) = (t0, y0) or (tR, yR)

     ARK_FULLRHS_END   -> called at the end of a successful step i.e, at
                          (tcur, ycur) or the start of the subsequent step i.e.,
                          at (tn, yn) = (tcur, ycur) from the end of the last
                          step

     ARK_FULLRHS_OTHER -> called elsewhere (e.g. for dense output)

  In SplittingStep, we accumulate the RHS functions in ARK_FULLRHS_OTHER mode.
  Generally, inner steppers will not have the correct yn when this function is
  called and will not be able to reuse a function evaluation since their state
  resets at the next SUNStepper_Evolve call.
  ----------------------------------------------------------------------------*/
fn splittingStep_FullRHS(
    ark_mem: &mut ARKodeMem,
    t: f64,
    y: &NVector,
    f: &mut NVector,
    _mode: i32,
) -> i32 {
    let mut step_mem = match splittingStep_AccessStepMem(ark_mem, "splittingStep_FullRHS") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };
    let mut tempv1 = std::mem::take(&mut ark_mem.tempv1);
    let retval = splittingStep_FullRHS_inner(ark_mem, &mut step_mem, t, y, f, &mut tempv1);
    ark_mem.tempv1 = tempv1;
    ark_mem.step_mem = Some(step_mem);
    retval
}

fn splittingStep_FullRHS_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut ARKodeSplittingStepMem,
    t: f64,
    y: &NVector,
    f: &mut NVector,
    tempv1: &mut NVector,
) -> i32 {
    for i in 0..step_mem.partitions as usize {
        let err = if i == 0 {
            SUNStepper_FullRhs(&mut step_mem.steppers[i], t, y, f, SUN_FULLRHS_OTHER)
        } else {
            SUNStepper_FullRhs(&mut step_mem.steppers[i], t, y, tempv1, SUN_FULLRHS_OTHER)
        };
        if err != SUN_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_RHSFUNC_FAIL,
                line!(),
                "splittingStep_FullRHS",
                file!(),
                &format!(
                    "At t = {}, the right-hand side routine failed in an unrecoverable manner.",
                    fmt_g(t, 0, 15)
                ),
            );
            return ARK_RHSFUNC_FAIL;
        }
        if i > 0 {
            /* C N_VLinearSum(ONE, f, ONE, tempv1, f) aliases z with x */
            f.linear_sum_with(ONE, ONE, tempv1);
        }
    }

    ARK_SUCCESS
}

/*------------------------------------------------------------------------------
  This routine performs a sequential operator splitting method
  ----------------------------------------------------------------------------*/
fn splittingStep_SequentialMethod(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut ARKodeSplittingStepMem,
    i: usize,
    y: &mut NVector,
) -> i32 {
    let stages = step_mem.coefficients.as_ref().unwrap().stages;
    let partitions = step_mem.coefficients.as_ref().unwrap().partitions;

    for j in 0..stages as usize {
        for k in 0..partitions as usize {
            let (beta_start, beta_end) = {
                let coefficients = step_mem.coefficients.as_ref().unwrap();
                (coefficients.beta[i][j][k], coefficients.beta[i][j + 1][k])
            };

            if beta_start == beta_end {
                continue;
            }

            let t_start = ark_mem.tn + beta_start * ark_mem.h;
            let t_end = ark_mem.tn + beta_end * ark_mem.h;

            let stepper = &mut step_mem.steppers[k];
            /* TODO(SBR): A potential future optimization is removing this reset and
             * a call to SUNStepper_SetStopTime later for methods that start a step
             * evolving the same partition the last step ended with (essentially a
             * FSAL property). Care is needed when a reset occurs, the step direction
             * changes, the coefficients change, etc. */
            let err = SUNStepper_Reset(stepper, t_start, y);
            if err != SUN_SUCCESS {
                return ARK_SUNSTEPPER_ERR;
            }

            let err = SUNStepper_SetStepDirection(stepper, t_end - t_start);
            if err != SUN_SUCCESS {
                return ARK_SUNSTEPPER_ERR;
            }

            let err = SUNStepper_SetStopTime(stepper, t_end);
            if err != SUN_SUCCESS {
                return ARK_SUNSTEPPER_ERR;
            }

            let mut tret = ZERO;
            let err = SUNStepper_Evolve(stepper, t_end, y, &mut tret);
            if err != SUN_SUCCESS {
                return ARK_SUNSTEPPER_ERR;
            }
            step_mem.n_stepper_evolves[k] += 1;
        }
    }

    ARK_SUCCESS
}

/*------------------------------------------------------------------------------
  This routine performs a single step of the splitting method.
  ----------------------------------------------------------------------------*/
fn splittingStep_TakeStep(ark_mem: &mut ARKodeMem, dsmPtr: &mut f64, nflagPtr: &mut i32) -> i32 {
    let mut step_mem = match splittingStep_AccessStepMem(ark_mem, "splittingStep_TakeStep") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };
    let mut ycur = std::mem::take(&mut ark_mem.ycur);
    let mut tempv1 = std::mem::take(&mut ark_mem.tempv1);
    let retval =
        splittingStep_TakeStep_inner(ark_mem, &mut step_mem, &mut ycur, &mut tempv1, dsmPtr, nflagPtr);
    ark_mem.ycur = ycur;
    ark_mem.tempv1 = tempv1;
    ark_mem.step_mem = Some(step_mem);
    retval
}

fn splittingStep_TakeStep_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut ARKodeSplittingStepMem,
    ycur: &mut NVector,
    tempv1: &mut NVector,
    dsmPtr: &mut f64,
    nflagPtr: &mut i32,
) -> i32 {
    *nflagPtr = ARK_SUCCESS; /* No algebraic solver */
    *dsmPtr = ZERO; /* No error estimate */

    step_mem.istage = 0;
    let retval = splittingStep_SequentialMethod(ark_mem, step_mem, 0, ycur);
    if retval != ARK_SUCCESS {
        return retval;
    }

    let alpha0 = step_mem.coefficients.as_ref().unwrap().alpha[0];
    if alpha0 != ONE {
        /* C N_VScale(alpha[0], ycur, ycur) operates in place */
        ycur.scale_inplace(alpha0);
    }

    let sequential_methods = step_mem.coefficients.as_ref().unwrap().sequential_methods;
    for i in 1..sequential_methods as usize {
        step_mem.istage = i as i32;

        N_VScale(ONE, &ark_mem.yn, tempv1);
        let retval = splittingStep_SequentialMethod(ark_mem, step_mem, i, tempv1);
        if retval != ARK_SUCCESS {
            return retval;
        }
        /* C N_VLinearSum(ONE, ycur, alpha[i], tempv1, ycur) aliases z with x */
        let alpha_i = step_mem.coefficients.as_ref().unwrap().alpha[i];
        ycur.linear_sum_with(ONE, alpha_i, tempv1);
    }

    ARK_SUCCESS
}

/*------------------------------------------------------------------------------
  Prints integrator statistics
  ----------------------------------------------------------------------------*/
fn splittingStep_PrintAllStats(
    ark_mem: &mut ARKodeMem,
    outfile: &mut dyn std::io::Write,
    fmt: SUNOutputFormat,
) -> i32 {
    let step_mem = match splittingStep_AccessStepMem(ark_mem, "splittingStep_PrintAllStats") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    for k in 0..step_mem.partitions as usize {
        /* C snprintf(name_buf, SUN_TABLE_WIDTH, ...) truncates to width-1 */
        let mut name_buf = format!("Partition {} evolves", k + 1);
        name_buf.truncate(SUN_TABLE_WIDTH - 1);
        crate::arkode_io::sunfprintf_long(outfile, fmt, false, &name_buf, step_mem.n_stepper_evolves[k]);
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*------------------------------------------------------------------------------
  Outputs all solver parameters to the provided file pointer.
  ----------------------------------------------------------------------------*/
fn splittingStep_WriteParameters(ark_mem: &mut ARKodeMem, fp: &mut dyn std::io::Write) -> i32 {
    let step_mem = match splittingStep_AccessStepMem(ark_mem, "splittingStep_WriteParameters") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    let _ = write!(
        fp,
        "SplittingStep time step module parameters:\n  Method order {}\n\n",
        step_mem.order
    );

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*------------------------------------------------------------------------------
  Frees all SplittingStep memory.
  ----------------------------------------------------------------------------*/
fn splittingStep_Free(ark_mem: &mut ARKodeMem) {
    /* ownership drop releases the steppers, counters and coefficients */
    ark_mem.step_mem = None;
}

/*------------------------------------------------------------------------------
  This routine outputs the memory from the SplittingStep structure to a
  specified file pointer (useful when debugging).
  ----------------------------------------------------------------------------*/
fn splittingStep_PrintMem(ark_mem: &mut ARKodeMem, outfile: &mut dyn std::io::Write) {
    let step_mem = match splittingStep_AccessStepMem(ark_mem, "splittingStep_PrintMem") {
        Some(sm) => sm,
        None => return,
    };

    /* output integer quantities */
    let _ = writeln!(outfile, "SplittingStep: istage = {}", step_mem.istage);
    let _ = writeln!(outfile, "SplittingStep: partitions = {}", step_mem.partitions);
    let _ = writeln!(outfile, "SplittingStep: order = {}", step_mem.order);

    /* output long integer quantities */
    for k in 0..step_mem.partitions as usize {
        let _ = writeln!(
            outfile,
            "SplittingStep: partition {}: n_stepper_evolves = {}",
            k, step_mem.n_stepper_evolves[k]
        );
    }

    /* output sunrealtype quantities */
    let _ = writeln!(outfile, "SplittingStep: Coefficients:");
    if let Some(coefficients) = step_mem.coefficients.as_ref() {
        SplittingStepCoefficients_Write(coefficients, outfile);
    }

    ark_mem.step_mem = Some(step_mem);
}

/*------------------------------------------------------------------------------
  Specifies the method order
  ----------------------------------------------------------------------------*/
fn splittingStep_SetOrder(ark_mem: &mut ARKodeMem, order: i32) -> i32 {
    let mut step_mem = match splittingStep_AccessStepMem(ark_mem, "splittingStep_SetOrder") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* set user-provided value, or default, depending on argument */
    step_mem.order = std::cmp::max(1, order);

    SplittingStepCoefficients_Destroy(&mut step_mem.coefficients);

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  Returns the current stage index and number of stages
  ---------------------------------------------------------------*/
fn splittingStep_GetStageIndex(ark_mem: &mut ARKodeMem, istage: &mut i32, num_stages: &mut i32) -> i32 {
    let step_mem = match splittingStep_AccessStepMem(ark_mem, "splittingStep_GetStageIndex") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };
    let retval = ARK_SUCCESS;

    /* if coefficients structure is not yet available, return defaults */
    if let Some(coefficients) = step_mem.coefficients.as_ref() {
        *istage = step_mem.istage;
        *num_stages = coefficients.sequential_methods;
    } else {
        *istage = -1;
        *num_stages = -1;
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_NULL,
            line!(),
            "splittingStep_GetStageIndex",
            file!(),
            "coefficient table not allocated",
        );
        /* C returns retval (== ARK_SUCCESS at this point) */
        return retval;
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*------------------------------------------------------------------------------
  Routine to set SplittingStep options
  ----------------------------------------------------------------------------*/
fn splittingStep_SetOptions(
    ark_mem: &mut ARKodeMem,
    argidx: &mut usize,
    argv: &[String],
    offset: usize,
    arg_used: &mut bool,
) -> i32 {
    /* The only SplittingStep-specific "Set" routine takes a custom set of
       coefficients; however, these may be specified by name, so here we'll support
       a key to specify the SplittingStepCoefficients by name,
       create the coefficients with that name, attach it to SplittingStep (who copies its
       values), and then frees the coefficients. */
    if &argv[*argidx][offset..] == "splitting_coefficients_name" {
        *argidx += 1;
        let mut Coefficients = SplittingStepCoefficients_LoadCoefficientsByName(&argv[*argidx]);
        if Coefficients.is_none() {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "splittingStep_SetOptions",
                file!(),
                &format!(
                    "error setting key {} {} (invalid coefficients name)",
                    argv[*argidx - 1],
                    argv[*argidx]
                ),
            );
            return ARK_ILL_INPUT;
        }
        let retval = SplittingStepSetCoefficients(ark_mem, Coefficients.as_deref().unwrap());
        SplittingStepCoefficients_Destroy(&mut Coefficients);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                retval,
                line!(),
                "splittingStep_SetOptions",
                file!(),
                &format!(
                    "error setting key {} {} (SetCoefficients failed)",
                    argv[*argidx - 1],
                    argv[*argidx]
                ),
            );
            return retval;
        }
        *arg_used = true;
        return ARK_SUCCESS;
    }

    ARK_SUCCESS
}

/*------------------------------------------------------------------------------
  Resets all SplittingStep optional inputs to their default values. Does not
  change problem-defining function pointers or user_data pointer.
  ----------------------------------------------------------------------------*/
fn splittingStep_SetDefaults(ark_mem: &mut ARKodeMem) -> i32 {
    /* C does an AccessStepMem null check before delegating; the delegate
    repeats it (take semantics forbid holding both) */
    if ark_mem.step_mem.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_NULL,
            line!(),
            "splittingStep_SetDefaults",
            file!(),
            "Time step module memory is NULL.",
        );
        return ARK_MEM_NULL;
    }

    splittingStep_SetOrder(ark_mem, 0)
}

/*------------------------------------------------------------------------------
  This routine checks if all required SUNStepper operations are present. If any
  of them are missing it return SUNFALSE.
  ----------------------------------------------------------------------------*/
fn splittingStep_CheckSUNStepper(stepper: &SUNStepper) -> bool {
    let ops = &stepper.ops;
    ops.evolve.is_some()
        && ops.reset.is_some()
        && ops.setstoptime.is_some()
        && ops.setstepdirection.is_some()
}

/*------------------------------------------------------------------------------
  This routine validates arguments when (re)initializing a SplittingStep
  integrator.  (C `steppers == NULL`, `steppers[i] == NULL` and `y0 == NULL`
  checks have no Rust equivalent.)
  ----------------------------------------------------------------------------*/
fn splittingStep_CheckArgs(
    ark_mem: Option<&ARKodeMem>,
    steppers: &[SUNStepper],
    partitions: i32,
) -> i32 {
    if partitions <= 1 {
        arkProcessError(
            ark_mem,
            ARK_ILL_INPUT,
            line!(),
            "splittingStep_CheckArgs",
            file!(),
            "The number of partitions must be greater than one",
        );
        return ARK_ILL_INPUT;
    }

    for (i, stepper) in steppers.iter().enumerate().take(partitions as usize) {
        if !splittingStep_CheckSUNStepper(stepper) {
            arkProcessError(
                ark_mem,
                ARK_ILL_INPUT,
                line!(),
                "splittingStep_CheckArgs",
                file!(),
                &format!("stepper[{}] does not implement the required operations.", i),
            );
            return ARK_ILL_INPUT;
        }
    }

    ARK_SUCCESS
}

/*------------------------------------------------------------------------------
  This routine initializes the step memory and resets the statistics
  ----------------------------------------------------------------------------*/
fn splittingStep_InitStepMem(
    _ark_mem: &mut ARKodeMem,
    step_mem: &mut ARKodeSplittingStepMem,
    steppers: Vec<SUNStepper>,
    partitions: i32,
) -> i32 {
    step_mem.steppers = steppers;
    step_mem.n_stepper_evolves = vec![0i64; partitions as usize];

    /* If the number of partitions changed, the coefficients are no longer
     * compatible and must be cleared. If a user previously called ARKodeSetOrder
     * that will still be respected at the next call to ARKodeEvolve */
    if step_mem.partitions != partitions {
        SplittingStepCoefficients_Destroy(&mut step_mem.coefficients);
    }
    step_mem.partitions = partitions;

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  Creates the SplittingStep integrator
  ---------------------------------------------------------------*/
pub fn SplittingStepCreate(
    steppers: Vec<SUNStepper>,
    partitions: i32,
    t0: f64,
    y0: &NVector,
    sunctx: &crate::sundials_context::SUNContext,
) -> Option<Box<ARKodeMem>> {
    let retval = splittingStep_CheckArgs(None, &steppers, partitions);
    if retval != ARK_SUCCESS {
        return None;
    }

    /* C checks sunctx == NULL (inexpressible here) */

    /* Create ark_mem structure and set default values */
    let mut ark_mem = arkCreate(sunctx);

    let mut step_mem = Box::new(ARKodeSplittingStepMem {
        steppers: Vec::new(),
        coefficients: None,
        n_stepper_evolves: Vec::new(),
        istage: 0,
        partitions,
        order: 0,
    });
    let retval = splittingStep_InitStepMem(&mut ark_mem, &mut step_mem, steppers, partitions);
    if retval != ARK_SUCCESS {
        return None;
    }

    /* Attach step_mem structure and function pointers to ark_mem */
    ark_mem.step_init = Some(splittingStep_Init);
    ark_mem.step_fullrhs = Some(splittingStep_FullRHS);
    ark_mem.step = Some(splittingStep_TakeStep);
    ark_mem.step_printallstats = Some(splittingStep_PrintAllStats);
    ark_mem.step_writeparameters = Some(splittingStep_WriteParameters);
    ark_mem.step_free = Some(splittingStep_Free);
    ark_mem.step_printmem = Some(splittingStep_PrintMem);
    ark_mem.step_setoptions = Some(splittingStep_SetOptions);
    ark_mem.step_setdefaults = Some(splittingStep_SetDefaults);
    ark_mem.step_setorder = Some(splittingStep_SetOrder);
    ark_mem.step_getstageindex = Some(splittingStep_GetStageIndex);
    ark_mem.step_mem = Some(step_mem);

    /* Set default values for ARKStep optional inputs */
    let retval = splittingStep_SetDefaults(&mut ark_mem);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(&ark_mem),
            retval,
            line!(),
            "SplittingStepCreate",
            file!(),
            "Error setting default solver options",
        );
        return None;
    }

    /* Initialize main ARKODE infrastructure */
    let retval = arkInit(&mut ark_mem, t0, y0, FIRST_INIT);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(&ark_mem),
            retval,
            line!(),
            "SplittingStepCreate",
            file!(),
            "Unable to initialize main ARKODE infrastructure",
        );
        return None;
    }

    ARKodeSetInterpolantType(&mut ark_mem, ARK_INTERP_LAGRANGE);

    Some(ark_mem)
}

/*------------------------------------------------------------------------------
  This routine re-initializes the SplittingStep module to solve a new problem of
  the same size as was previously solved. This routine should also be called
  when the problem dynamics or desired solvers have changed dramatically, so
  that the problem integration should resume as if started from scratch.

  Note all internal counters are set to 0 on re-initialization.
  ----------------------------------------------------------------------------*/
pub fn SplittingStepReInit(
    ark_mem: &mut ARKodeMem,
    steppers: Vec<SUNStepper>,
    partitions: i32,
    t0: f64,
    y0: &NVector,
) -> i32 {
    let mut step_mem = match splittingStep_AccessStepMem(ark_mem, "SplittingStepReInit") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* Check if ark_mem was allocated */
    if !ark_mem.MallocDone {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_NO_MALLOC,
            line!(),
            "SplittingStepReInit",
            file!(),
            MSG_ARK_NO_MALLOC,
        );
        return ARK_NO_MALLOC;
    }

    let retval = splittingStep_CheckArgs(Some(ark_mem), &steppers, partitions);
    if retval != ARK_SUCCESS {
        ark_mem.step_mem = Some(step_mem);
        return retval;
    }

    splittingStep_InitStepMem(ark_mem, &mut step_mem, steppers, partitions);
    ark_mem.step_mem = Some(step_mem);

    /* Initialize main ARKODE infrastructure */
    let retval = arkInit(ark_mem, t0, y0, FIRST_INIT);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            retval,
            line!(),
            "SplittingStepReInit",
            file!(),
            "Unable to initialize main ARKODE infrastructure",
        );
        return retval;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  Sets the SplittingStep coefficients.
  ---------------------------------------------------------------*/
pub fn SplittingStepSetCoefficients(
    ark_mem: &mut ARKodeMem,
    coefficients: &SplittingStepCoefficientsMem,
) -> i32 {
    let mut step_mem = match splittingStep_AccessStepMem(ark_mem, "SplittingStepSetCoefficients") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* C checks coefficients == NULL (inexpressible here) */

    if step_mem.partitions != coefficients.partitions {
        let msg = format!(
            "The splitting method has {} partitions but the coefficients have {}.",
            step_mem.partitions, coefficients.partitions
        );
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "SplittingStepSetCoefficients",
            file!(),
            &msg,
        );
        return ARK_ILL_INPUT;
    }

    SplittingStepCoefficients_Destroy(&mut step_mem.coefficients);
    step_mem.coefficients = SplittingStepCoefficients_Copy(coefficients);
    if step_mem.coefficients.is_none() {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_FAIL,
            line!(),
            "SplittingStepSetCoefficients",
            file!(),
            "Failed to copy splitting coefficients",
        );
        return ARK_MEM_NULL;
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*------------------------------------------------------------------------------
  Accesses the number of times a given partition was evolved
  ----------------------------------------------------------------------------*/
pub fn SplittingStepGetNumEvolves(ark_mem: &mut ARKodeMem, partition: i32, evolves: &mut i64) -> i32 {
    let step_mem = match splittingStep_AccessStepMem(ark_mem, "SplittingStepGetNumEvolves") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    if partition >= step_mem.partitions {
        let msg = format!(
            "The partition index is {} but there are only {} partitions",
            partition, step_mem.partitions
        );
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "SplittingStepGetNumEvolves",
            file!(),
            &msg,
        );
        return ARK_ILL_INPUT;
    }

    if partition < 0 {
        *evolves = 0;
        for k in 0..step_mem.partitions as usize {
            *evolves += step_mem.n_stepper_evolves[k];
        }
    } else {
        *evolves = step_mem.n_stepper_evolves[partition as usize];
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}
