/* -----------------------------------------------------------------
 * Translation of src/arkode/arkode_forcingstep.c (SUNDIALS 7.7.0).
 *
 * This is the implementation file for ARKODE's forcing method.
 *
 * step_mem follows the crate's take/downcast/put-back Access pattern;
 * TakeStep additionally mem::takes the ark_mem vectors it passes to
 * the inner SUNSteppers (ycur/tempv1) and restores them on every
 * return path.  C NULL checks on stepper/y0 arguments have no Rust
 * equivalent and are noted at their sites.
 * -----------------------------------------------------------------*/

use crate::arkode::{arkCreate, arkInit};
use crate::arkode_forcingstep_impl::{ARKodeForcingStepMem, NUM_PARTITIONS};
use crate::arkode_impl::{
    arkProcessError, ARKodeMem, ARK_FULLRHS_END, ARK_ILL_INPUT, ARK_INTERP_HERMITE,
    ARK_INTERP_LAGRANGE, ARK_MEM_NULL, ARK_NO_MALLOC, ARK_RHSFUNC_FAIL, ARK_SUCCESS,
    ARK_SUNSTEPPER_ERR, FIRST_INIT, RESET_INIT, RESIZE_INIT, MSG_ARK_NO_MALLOC,
    ZERO,
};
use crate::sundials_types::SUNOutputFormat;
use crate::arkode_io::ARKodeSetInterpolantType;
use crate::nvector_serial::{NVector, N_VLinearSum};
use crate::sundials_errors::SUN_SUCCESS;
use crate::sundials_stepper::{
    SUNStepper, SUNStepper_Evolve, SUNStepper_FullRhs, SUNStepper_Reset, SUNStepper_SetForcing,
    SUNStepper_SetStepDirection, SUNStepper_SetStopTime, SUN_FULLRHS_END, SUN_FULLRHS_OTHER,
};
use crate::sundials_utils::fmt_g;

const ONE: f64 = 1.0;

/*------------------------------------------------------------------------------
  Shortcut routine to unpack step_mem structure from ark_mem. If missing it
  returns ARK_MEM_NULL.  (Take semantics: caller must put the box back.)
  ----------------------------------------------------------------------------*/
pub(crate) fn forcingStep_AccessStepMem(
    ark_mem: &mut ARKodeMem,
    fname: &str,
) -> Option<Box<ARKodeForcingStepMem>> {
    let taken = ark_mem.step_mem.take();
    match taken {
        Some(b) => match b.downcast::<ARKodeForcingStepMem>() {
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
  This routine is called just prior to performing internal time steps (after
  all user "set" routines have been called) from within arkInitialSetup.
  ----------------------------------------------------------------------------*/
fn forcingStep_Init(ark_mem: &mut ARKodeMem, init_type: i32) -> i32 {
    let mut step_mem = match forcingStep_AccessStepMem(ark_mem, "forcingStep_Init") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };
    let retval = forcingStep_Init_inner(ark_mem, &mut step_mem, init_type);
    ark_mem.step_mem = Some(step_mem);
    retval
}

fn forcingStep_Init_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut ARKodeForcingStepMem,
    init_type: i32,
) -> i32 {
    /* assume fixed outer step size */
    if !ark_mem.fixedstep {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "forcingStep_Init",
            file!(),
            "Adaptive outer time stepping is not currently supported",
        );
        return ARK_ILL_INPUT;
    }

    if ark_mem.interp_type == ARK_INTERP_HERMITE
        && (step_mem.stepper[0].ops.fullrhs.is_none() || step_mem.stepper[1].ops.fullrhs.is_none())
    {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "forcingStep_Init",
            file!(),
            "The SUNSteppers must implement SUNStepper_FullRhs when using Hermite interpolation",
        );
        return ARK_ILL_INPUT;
    }

    /* immediately return if resize or reset */
    if init_type == RESIZE_INIT || init_type == RESET_INIT {
        return ARK_SUCCESS;
    }

    /* On first initialization, make the SUNStepper consistent with the current
     * state in case a user provided a different initial condition for the
     * ForcingStep integrator and SUNStepper. */
    let err = SUNStepper_Reset(&mut step_mem.stepper[1], ark_mem.tn, &ark_mem.yn);
    if err != SUN_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            ARK_SUNSTEPPER_ERR,
            line!(),
            "forcingStep_Init",
            file!(),
            "Resetting the second partition SUNStepper failed",
        );
        return ARK_SUNSTEPPER_ERR;
    }

    ark_mem.interp_degree = 1;

    ARK_SUCCESS
}

/*------------------------------------------------------------------------------
  This routine resets the ForcingStep integrator by resetting the partition
  integrators
  ----------------------------------------------------------------------------*/
fn forcingStep_Reset(ark_mem: &mut ARKodeMem, tR: f64, yR: &NVector) -> i32 {
    let mut step_mem = match forcingStep_AccessStepMem(ark_mem, "forcingStep_Reset") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };
    let retval = forcingStep_Reset_inner(ark_mem, &mut step_mem, tR, yR);
    ark_mem.step_mem = Some(step_mem);
    retval
}

fn forcingStep_Reset_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut ARKodeForcingStepMem,
    tR: f64,
    yR: &NVector,
) -> i32 {
    let err = SUNStepper_Reset(&mut step_mem.stepper[0], tR, yR);
    if err != SUN_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            ARK_SUNSTEPPER_ERR,
            line!(),
            "forcingStep_Reset",
            file!(),
            "Resetting the first partition SUNStepper failed",
        );
        return ARK_SUNSTEPPER_ERR;
    }

    let err = SUNStepper_Reset(&mut step_mem.stepper[1], tR, yR);
    if err != SUN_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            ARK_SUNSTEPPER_ERR,
            line!(),
            "forcingStep_Reset",
            file!(),
            "Resetting the second partition SUNStepper failed",
        );
        return ARK_SUNSTEPPER_ERR;
    }

    ARK_SUCCESS
}

/*------------------------------------------------------------------------------
  This routine sets the step direction of the partition integrators and is
  called once the ForcingStep integrator has updated its step direction.
  ----------------------------------------------------------------------------*/
fn forcingStep_SetStepDirection(ark_mem: &mut ARKodeMem, stepdir: f64) -> i32 {
    let mut step_mem = match forcingStep_AccessStepMem(ark_mem, "forcingStep_SetStepDirection") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };
    let retval = forcingStep_SetStepDirection_inner(ark_mem, &mut step_mem, stepdir);
    ark_mem.step_mem = Some(step_mem);
    retval
}

fn forcingStep_SetStepDirection_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut ARKodeForcingStepMem,
    stepdir: f64,
) -> i32 {
    let err = SUNStepper_SetStepDirection(&mut step_mem.stepper[0], stepdir);
    if err != SUN_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            ARK_SUNSTEPPER_ERR,
            line!(),
            "forcingStep_SetStepDirection",
            file!(),
            "Setting the step direction for the first partition SUNStepper failed",
        );
        return ARK_SUNSTEPPER_ERR;
    }

    let err = SUNStepper_SetStepDirection(&mut step_mem.stepper[1], stepdir);
    if err != SUN_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            ARK_SUNSTEPPER_ERR,
            line!(),
            "forcingStep_SetStepDirection",
            file!(),
            "Setting the step direction for the second partition SUNStepper failed",
        );
        return ARK_SUNSTEPPER_ERR;
    }

    ARK_SUCCESS
}

/*------------------------------------------------------------------------------
  This is just a wrapper to call the user-supplied RHS function,
  f^1(t,y) + f^2(t,y).

  The stepper for partition 1 has a state that is inconsistent with the
  ForcingStep integrator, so we cannot pass it the SUN_FULLRHS_END option. For
  partition 2, the state should be consistent, and we can use SUN_FULLRHS_END.
  ----------------------------------------------------------------------------*/
fn forcingStep_FullRHS(ark_mem: &mut ARKodeMem, t: f64, y: &NVector, f: &mut NVector, mode: i32) -> i32 {
    let mut step_mem = match forcingStep_AccessStepMem(ark_mem, "forcingStep_FullRHS") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };
    let mut tempv1 = std::mem::take(&mut ark_mem.tempv1);
    let retval = forcingStep_FullRHS_inner(ark_mem, &mut step_mem, t, y, f, &mut tempv1, mode);
    ark_mem.tempv1 = tempv1;
    ark_mem.step_mem = Some(step_mem);
    retval
}

fn forcingStep_FullRHS_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut ARKodeForcingStepMem,
    t: f64,
    y: &NVector,
    f: &mut NVector,
    tempv1: &mut NVector,
    mode: i32,
) -> i32 {
    /* TODO(SBR): Possible optimization in FULLRHS_START mode. Currently that
     * mode is not forwarded to the SUNSteppers */
    let err = SUNStepper_FullRhs(&mut step_mem.stepper[0], t, y, tempv1, SUN_FULLRHS_OTHER);
    if err != SUN_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            ARK_RHSFUNC_FAIL,
            line!(),
            "forcingStep_FullRHS",
            file!(),
            &format!(
                "At t = {}, the right-hand side routine failed in an unrecoverable manner.",
                fmt_g(t, 0, 15)
            ),
        );
        return ARK_RHSFUNC_FAIL;
    }

    let err = SUNStepper_FullRhs(
        &mut step_mem.stepper[1],
        t,
        y,
        f,
        if mode == ARK_FULLRHS_END {
            SUN_FULLRHS_END
        } else {
            SUN_FULLRHS_OTHER
        },
    );
    if err != SUN_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            ARK_RHSFUNC_FAIL,
            line!(),
            "forcingStep_FullRHS",
            file!(),
            &format!(
                "At t = {}, the right-hand side routine failed in an unrecoverable manner.",
                fmt_g(t, 0, 15)
            ),
        );
        return ARK_RHSFUNC_FAIL;
    }
    /* C N_VLinearSum(ONE, f, ONE, tempv1, f) aliases z with x */
    f.linear_sum_with(ONE, ONE, tempv1);

    ARK_SUCCESS
}

/*------------------------------------------------------------------------------
  This routine performs a single step of the forcing method.
  ----------------------------------------------------------------------------*/
fn forcingStep_TakeStep(ark_mem: &mut ARKodeMem, dsmPtr: &mut f64, nflagPtr: &mut i32) -> i32 {
    let mut step_mem = match forcingStep_AccessStepMem(ark_mem, "forcingStep_TakeStep") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };
    let mut ycur = std::mem::take(&mut ark_mem.ycur);
    let mut tempv1 = std::mem::take(&mut ark_mem.tempv1);
    let retval =
        forcingStep_TakeStep_inner(ark_mem, &mut step_mem, &mut ycur, &mut tempv1, dsmPtr, nflagPtr);
    ark_mem.ycur = ycur;
    ark_mem.tempv1 = tempv1;
    ark_mem.step_mem = Some(step_mem);
    retval
}

fn forcingStep_TakeStep_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut ARKodeForcingStepMem,
    ycur: &mut NVector,
    tempv1: &mut NVector,
    dsmPtr: &mut f64,
    nflagPtr: &mut i32,
) -> i32 {
    *nflagPtr = ARK_SUCCESS; /* No algebraic solver */
    *dsmPtr = ZERO; /* No error estimate */

    let tout = ark_mem.tn + ark_mem.h;
    let mut tret = ZERO;

    /* Evolve stepper 0 on its own */
    let s0 = &mut step_mem.stepper[0];
    let err = SUNStepper_Reset(s0, ark_mem.tn, &ark_mem.yn);
    if err != SUN_SUCCESS {
        return ARK_SUNSTEPPER_ERR;
    }

    let err = SUNStepper_SetStopTime(s0, tout);
    if err != SUN_SUCCESS {
        return ARK_SUNSTEPPER_ERR;
    }

    let err = SUNStepper_Evolve(s0, tout, ycur, &mut tret);
    if err != SUN_SUCCESS {
        return ARK_SUNSTEPPER_ERR;
    }
    step_mem.n_stepper_evolves[0] += 1;

    let s1 = &mut step_mem.stepper[1];
    /* A reset is not needed because steeper 1's state is consistent with the
     * forcing method */
    let err = SUNStepper_SetStopTime(s1, tout);
    if err != SUN_SUCCESS {
        return ARK_SUNSTEPPER_ERR;
    }

    /* Write tendency (ycur - yn)/h into stepper 1 forcing */
    let hinv = ONE / ark_mem.h;
    N_VLinearSum(hinv, ycur, -hinv, &ark_mem.yn, tempv1);
    let err = SUNStepper_SetForcing(s1, ZERO, ZERO, std::slice::from_ref(tempv1), 1);
    if err != SUN_SUCCESS {
        return ARK_SUNSTEPPER_ERR;
    }

    /* Evolve stepper 1 with the forcing */
    let err = SUNStepper_Evolve(s1, tout, ycur, &mut tret);
    if err != SUN_SUCCESS {
        return ARK_SUNSTEPPER_ERR;
    }
    step_mem.n_stepper_evolves[1] += 1;

    /* Clear the forcing so it doesn't get included in a fullRhs call */
    let err = SUNStepper_SetForcing(s1, ZERO, ZERO, &[], 0);
    if err != SUN_SUCCESS {
        return ARK_SUNSTEPPER_ERR;
    }

    ARK_SUCCESS
}

/*------------------------------------------------------------------------------
  Prints integrator statistics
  ----------------------------------------------------------------------------*/
fn forcingStep_PrintAllStats(
    ark_mem: &mut ARKodeMem,
    outfile: &mut dyn std::io::Write,
    fmt: SUNOutputFormat,
) -> i32 {
    let step_mem = match forcingStep_AccessStepMem(ark_mem, "forcingStep_PrintAllStats") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    crate::arkode_io::sunfprintf_long(
        outfile,
        fmt,
        false,
        "Partition 1 evolves",
        step_mem.n_stepper_evolves[0],
    );
    crate::arkode_io::sunfprintf_long(
        outfile,
        fmt,
        false,
        "Partition 2 evolves",
        step_mem.n_stepper_evolves[1],
    );

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*------------------------------------------------------------------------------
  Frees all ForcingStep memory.
  ----------------------------------------------------------------------------*/
fn forcingStep_Free(ark_mem: &mut ARKodeMem) {
    ark_mem.step_mem = None;
}

/*------------------------------------------------------------------------------
  This routine outputs the memory from the ForcingStep structure to a specified
  file pointer (useful when debugging).
  ----------------------------------------------------------------------------*/
fn forcingStep_PrintMem(ark_mem: &mut ARKodeMem, outfile: &mut dyn std::io::Write) {
    let step_mem = match forcingStep_AccessStepMem(ark_mem, "forcingStep_PrintMem") {
        Some(sm) => sm,
        None => return,
    };

    /* output long integer quantities */
    for k in 0..NUM_PARTITIONS {
        let _ = writeln!(
            outfile,
            "ForcingStep: partition {}: n_stepper_evolves = {}",
            k, step_mem.n_stepper_evolves[k]
        );
    }

    ark_mem.step_mem = Some(step_mem);
}

/*------------------------------------------------------------------------------
  This routine checks if all required SUNStepper operations are present. If any
  of them are missing it return SUNFALSE.
  ----------------------------------------------------------------------------*/
fn forcingStep_CheckSUNStepper(stepper: &SUNStepper, needs_forcing: bool) -> bool {
    let ops = &stepper.ops;
    ops.evolve.is_some()
        && ops.reset.is_some()
        && ops.setstoptime.is_some()
        && (!needs_forcing || ops.setforcing.is_some())
}

/*------------------------------------------------------------------------------
  This routine validates arguments when (re)initializing a ForcingStep
  integrator.  (C `stepper1 == NULL`, `stepper2 == NULL` and `y0 == NULL`
  checks have no Rust equivalent.)
  ----------------------------------------------------------------------------*/
fn forcingStep_CheckArgs(
    ark_mem: Option<&ARKodeMem>,
    stepper1: &SUNStepper,
    stepper2: &SUNStepper,
) -> i32 {
    if !forcingStep_CheckSUNStepper(stepper1, false) {
        arkProcessError(
            ark_mem,
            ARK_ILL_INPUT,
            line!(),
            "forcingStep_CheckArgs",
            file!(),
            "stepper1 does not implement the required operations.",
        );
        return ARK_ILL_INPUT;
    }

    if !forcingStep_CheckSUNStepper(stepper2, true) {
        arkProcessError(
            ark_mem,
            ARK_ILL_INPUT,
            line!(),
            "forcingStep_CheckArgs",
            file!(),
            "stepper2 does not implement the required operations.",
        );
        return ARK_ILL_INPUT;
    }

    ARK_SUCCESS
}

/*------------------------------------------------------------------------------
  This routine initializes the step memory and resets the statistics
  ----------------------------------------------------------------------------*/
fn forcingStep_InitStepMem(
    step_mem: &mut ARKodeForcingStepMem,
    stepper1: SUNStepper,
    stepper2: SUNStepper,
) {
    step_mem.stepper[0] = stepper1;
    step_mem.stepper[1] = stepper2;
    step_mem.n_stepper_evolves[0] = 0;
    step_mem.n_stepper_evolves[1] = 0;
}

/*------------------------------------------------------------------------------
  Creates the ForcingStep integrator
  ----------------------------------------------------------------------------*/
pub fn ForcingStepCreate(
    stepper1: SUNStepper,
    stepper2: SUNStepper,
    t0: f64,
    y0: &NVector,
    sunctx: &crate::sundials_context::SUNContext,
) -> Option<Box<ARKodeMem>> {
    let retval = forcingStep_CheckArgs(None, &stepper1, &stepper2);
    if retval != ARK_SUCCESS {
        return None;
    }

    /* C checks sunctx == NULL (inexpressible here) */

    /* Create ark_mem structure and set default values */
    let mut ark_mem = arkCreate(sunctx);

    let step_mem = Box::new(ARKodeForcingStepMem {
        stepper: [stepper1, stepper2],
        n_stepper_evolves: [0, 0],
    });

    /* Attach step_mem structure and function pointers to ark_mem */
    ark_mem.step_init = Some(forcingStep_Init);
    ark_mem.step_fullrhs = Some(forcingStep_FullRHS);
    ark_mem.step_reset = Some(forcingStep_Reset);
    ark_mem.step_setstepdirection = Some(forcingStep_SetStepDirection);
    ark_mem.step = Some(forcingStep_TakeStep);
    ark_mem.step_printallstats = Some(forcingStep_PrintAllStats);
    ark_mem.step_free = Some(forcingStep_Free);
    ark_mem.step_printmem = Some(forcingStep_PrintMem);
    ark_mem.step_mem = Some(step_mem);

    /* Initialize main ARKODE infrastructure */
    let retval = arkInit(&mut ark_mem, t0, y0, FIRST_INIT);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(&ark_mem),
            retval,
            line!(),
            "ForcingStepCreate",
            file!(),
            "Unable to initialize main ARKODE infrastructure",
        );
        return None;
    }

    ARKodeSetInterpolantType(&mut ark_mem, ARK_INTERP_LAGRANGE);

    Some(ark_mem)
}

/*------------------------------------------------------------------------------
  This routine re-initializes the ForcingStep module to solve a new problem of
  the same size as was previously solved. This routine should also be called
  when the problem dynamics or desired solvers have changed dramatically, so
  that the problem integration should resume as if started from scratch.

  Note all internal counters are set to 0 on re-initialization.
  ----------------------------------------------------------------------------*/
pub fn ForcingStepReInit(
    ark_mem: &mut ARKodeMem,
    stepper1: SUNStepper,
    stepper2: SUNStepper,
    t0: f64,
    y0: &NVector,
) -> i32 {
    let mut step_mem = match forcingStep_AccessStepMem(ark_mem, "ForcingStepReInit") {
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
            "ForcingStepReInit",
            file!(),
            MSG_ARK_NO_MALLOC,
        );
        return ARK_NO_MALLOC;
    }

    let retval = forcingStep_CheckArgs(Some(ark_mem), &stepper1, &stepper2);
    if retval != ARK_SUCCESS {
        ark_mem.step_mem = Some(step_mem);
        return retval;
    }

    forcingStep_InitStepMem(&mut step_mem, stepper1, stepper2);
    ark_mem.step_mem = Some(step_mem);

    /* Initialize main ARKODE infrastructure */
    let retval = arkInit(ark_mem, t0, y0, FIRST_INIT);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            retval,
            line!(),
            "ForcingStepReInit",
            file!(),
            "Unable to initialize main ARKODE infrastructure",
        );
        return retval;
    }

    ARK_SUCCESS
}

/*------------------------------------------------------------------------------
  Accesses the number of times a given partition was evolved
  ----------------------------------------------------------------------------*/
pub fn ForcingStepGetNumEvolves(ark_mem: &mut ARKodeMem, partition: i32, evolves: &mut i64) -> i32 {
    let step_mem = match forcingStep_AccessStepMem(ark_mem, "ForcingStepGetNumEvolves") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    if partition >= NUM_PARTITIONS as i32 {
        let msg = format!(
            "The partition index is {} but there are only 2 partitions",
            partition
        );
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "ForcingStepGetNumEvolves",
            file!(),
            &msg,
        );
        return ARK_ILL_INPUT;
    }

    if partition < 0 {
        *evolves = step_mem.n_stepper_evolves[0] + step_mem.n_stepper_evolves[1];
    } else {
        *evolves = step_mem.n_stepper_evolves[partition as usize];
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}
