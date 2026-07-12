/* -----------------------------------------------------------------
 * Translation of sundials-7.7.0/src/arkode/arkode_erkstep_io.c
 * (PART I: the erkStep_* option/output routines attached to ARKODE
 * and the ERKStepSetTable family; the deprecated ERKStep* wrapper
 * aliases simply forward to the generic ARKode* routines and are
 * added as needed by the examples).
 * -----------------------------------------------------------------*/

use crate::arkode_butcher::{
    ARKodeButcherTable, ARKodeButcherTable_Copy, ARKodeButcherTable_Space,
};
use crate::arkode_butcher_erk::{
    arkButcherTableERKNameToID, ARKODE_ERKTableID, ARKodeButcherTable_LoadERK,
    ARKODE_MAX_ERK_NUM, ARKODE_MIN_ERK_NUM,
};
use crate::arkode_erkstep::erkStep_AccessStepMem;
use crate::arkode_impl::{
    arkProcessError, ARKodeMem, ARK_ACCUMERROR_NONE, ARK_ILL_INPUT, ARK_MEM_NULL,
    ARK_STEPPER_UNSUPPORTED, ARK_SUCCESS, ONE, Q_DEFAULT,
};
use crate::nvector_serial::{N_VScale, NVector};
use crate::sundials_types::{SUNOutputFormat, SUN_OUTPUTFORMAT_TABLE};

/*===============================================================
  Exported optional input functions.
  ===============================================================*/

/*---------------------------------------------------------------
  ERKStepSetTable:

  Specifies to use a customized Butcher table for the explicit
  portion of the system.
  ---------------------------------------------------------------*/
pub fn ERKStepSetTable(ark_mem: &mut ARKodeMem, b: &ARKodeButcherTable) -> i32 {
    /* access ARKodeERKStepMem structure */
    let mut step_mem = match erkStep_AccessStepMem(ark_mem, "ERKStepSetTable") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* clear any existing parameters and Butcher tables */
    step_mem.stages = 0;
    step_mem.q = 0;
    step_mem.p = 0;

    let (mut bliw, mut blrw) = (0i64, 0i64);
    if let Some(old) = step_mem.B.take() {
        ARKodeButcherTable_Space(&old, &mut bliw, &mut blrw);
    }
    ark_mem.liw -= bliw;
    ark_mem.lrw -= blrw;

    /* set the relevant parameters */
    step_mem.stages = b.stages;
    step_mem.q = b.q;
    step_mem.p = b.p;

    /* copy the table into step memory */
    step_mem.B = ARKodeButcherTable_Copy(b);

    let (mut bliw, mut blrw) = (0i64, 0i64);
    ARKodeButcherTable_Space(step_mem.B.as_ref().unwrap(), &mut bliw, &mut blrw);
    ark_mem.liw += bliw;
    ark_mem.lrw += blrw;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ERKStepSetTableNum:

  Specifies to use a pre-existing Butcher table for the problem,
  based on the integer flag passed to ARKodeButcherTable_LoadERK()
  within the file arkode_butcher_erk.c.
  ---------------------------------------------------------------*/
pub fn ERKStepSetTableNum(ark_mem: &mut ARKodeMem, etable: ARKODE_ERKTableID) -> i32 {
    /* access ARKodeERKStepMem structure */
    let mut step_mem = match erkStep_AccessStepMem(ark_mem, "ERKStepSetTableNum") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* check that argument specifies an explicit table */
    if !(ARKODE_MIN_ERK_NUM..=ARKODE_MAX_ERK_NUM).contains(&etable) {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_NULL,
            line!(),
            "ERKStepSetTableNum",
            file!(),
            "Illegal ERK table number",
        );
        return ARK_ILL_INPUT;
    }

    /* clear any existing parameters and Butcher tables */
    step_mem.stages = 0;
    step_mem.q = 0;
    step_mem.p = 0;

    let (mut bliw, mut blrw) = (0i64, 0i64);
    if let Some(old) = step_mem.B.take() {
        ARKodeButcherTable_Space(&old, &mut bliw, &mut blrw);
    }
    ark_mem.liw -= bliw;
    ark_mem.lrw -= blrw;

    /* fill in table based on argument */
    step_mem.B = ARKodeButcherTable_LoadERK(etable);
    if step_mem.B.is_none() {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_NULL,
            line!(),
            "ERKStepSetTableNum",
            file!(),
            "Error setting table with that index",
        );
        return ARK_ILL_INPUT;
    }
    {
        let b = step_mem.B.as_ref().unwrap();
        step_mem.stages = b.stages;
        step_mem.q = b.q;
        step_mem.p = b.p;
    }

    let (mut bliw, mut blrw) = (0i64, 0i64);
    ARKodeButcherTable_Space(step_mem.B.as_ref().unwrap(), &mut bliw, &mut blrw);
    ark_mem.liw += bliw;
    ark_mem.lrw += blrw;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ERKStepSetTableName:

  Specifies to use a pre-existing Butcher table for the problem,
  based on the string name of the table.
  ---------------------------------------------------------------*/
pub fn ERKStepSetTableName(ark_mem: &mut ARKodeMem, etable: &str) -> i32 {
    ERKStepSetTableNum(ark_mem, arkButcherTableERKNameToID(etable))
}

/*===============================================================
  Exported optional output functions.
  ===============================================================*/

/*---------------------------------------------------------------
  erkStep_GetNumRhsEvals:

  Returns the current number of RHS calls
  ---------------------------------------------------------------*/
pub fn erkStep_GetNumRhsEvals(
    ark_mem: &mut ARKodeMem,
    partition_index: i32,
    rhs_evals: &mut i64,
) -> i32 {
    /* access ARKodeERKStepMem structure */
    let step_mem = match erkStep_AccessStepMem(ark_mem, "erkStep_GetNumRhsEvals") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    if partition_index > 0 {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "erkStep_GetNumRhsEvals",
            file!(),
            "Invalid partition index",
        );
        return ARK_ILL_INPUT;
    }

    *rhs_evals = step_mem.nfe;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ERKStepGetCurrentButcherTable:

  Returns a copy of the Butcher table currently in use (C hands
  out the internal pointer).
  ---------------------------------------------------------------*/
pub fn ERKStepGetCurrentButcherTable(ark_mem: &mut ARKodeMem) -> Option<ARKodeButcherTable> {
    let step_mem = erkStep_AccessStepMem(ark_mem, "ERKStepGetCurrentButcherTable")?;
    let b = step_mem
        .B
        .as_ref()
        .and_then(ARKodeButcherTable_Copy);
    ark_mem.step_mem = Some(step_mem);
    b
}

/*---------------------------------------------------------------
  ERKStepGetTimestepperStats:

  Returns integrator statistics
  ---------------------------------------------------------------*/
pub fn ERKStepGetTimestepperStats(
    ark_mem: &mut ARKodeMem,
    expsteps: &mut i64,
    accsteps: &mut i64,
    attempts: &mut i64,
    fevals: &mut i64,
    netfails: &mut i64,
) -> i32 {
    let step_mem = match erkStep_AccessStepMem(ark_mem, "ERKStepGetTimestepperStats") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* set expsteps and accsteps from adaptivity structure */
    let hadapt_mem = ark_mem.hadapt_mem.as_ref().unwrap();
    *expsteps = hadapt_mem.nst_exp;
    *accsteps = hadapt_mem.nst_acc;

    /* set remaining outputs */
    *attempts = ark_mem.nst_attempts;
    *fevals = step_mem.nfe;
    *netfails = ark_mem.netf;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*===============================================================
  Private functions attached to ARKODE
  ===============================================================*/

/*---------------------------------------------------------------
  ERKStepSetAdaptivityMethod: user-callable deprecated wrapper
  around arkSetAdaptivityMethod.
  ---------------------------------------------------------------*/
pub fn ERKStepSetAdaptivityMethod(
    ark_mem: &mut ARKodeMem,
    imethod: i32,
    idefault: i32,
    pq: i32,
    adapt_params: Option<&[f64; 3]>,
) -> i32 {
    crate::arkode_io::arkSetAdaptivityMethod(ark_mem, imethod, idefault, pq, adapt_params)
}

/*---------------------------------------------------------------
  erkStep_SetRelaxFn:

  Sets up the relaxation module using ERKStep's utility routines.
  ---------------------------------------------------------------*/
pub fn erkStep_SetRelaxFn(
    ark_mem: &mut ARKodeMem,
    rfn: Option<crate::arkode_impl::ARKRelaxFn>,
    rjac: Option<crate::arkode_impl::ARKRelaxJacFn>,
) -> i32 {
    crate::arkode_relaxation::arkRelaxCreate(
        ark_mem,
        rfn,
        rjac,
        Some(crate::arkode_erkstep::erkStep_RelaxDeltaE),
        Some(crate::arkode_erkstep::erkStep_GetOrder),
    )
}

/*---------------------------------------------------------------
  erkStep_SetDefaults:

  Resets all ERKStep optional inputs to their default values.
  Does not change problem-defining function pointers or
  user_data pointer.
  ---------------------------------------------------------------*/
pub fn erkStep_SetDefaults(ark_mem: &mut ARKodeMem) -> i32 {
    /* access ARKodeERKStepMem structure */
    let mut step_mem = match erkStep_AccessStepMem(ark_mem, "erkStep_SetDefaults") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* Set default values for integrator optional inputs */
    step_mem.q = Q_DEFAULT; /* method order */
    step_mem.p = 0; /* embedding order */
    step_mem.stages = 0; /* no stages */
    {
        let hadapt_mem = ark_mem.hadapt_mem.as_mut().unwrap();
        hadapt_mem.etamxf = 0.3; /* max change on error-failed step */
        hadapt_mem.safety = 0.99; /* step adaptivity safety factor  */
        hadapt_mem.growth = 25.0; /* step adaptivity growth factor */
    }

    /* Remove pre-existing Butcher table */
    if let Some(old) = step_mem.B.take() {
        let (mut bliw, mut blrw) = (0i64, 0i64);
        ARKodeButcherTable_Space(&old, &mut bliw, &mut blrw);
        ark_mem.liw -= bliw;
        ark_mem.lrw -= blrw;
    }
    step_mem.B = None;

    ark_mem.step_mem = Some(step_mem);

    /* Load the default SUNAdaptController */
    let retval = crate::arkode_io::arkReplaceAdaptController(ark_mem, None, true);
    if retval != 0 {
        return retval;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  erkStep_SetOrder:

  Specifies the method order
  ---------------------------------------------------------------*/
pub fn erkStep_SetOrder(ark_mem: &mut ARKodeMem, ord: i32) -> i32 {
    /* access ARKodeERKStepMem structure */
    let mut step_mem = match erkStep_AccessStepMem(ark_mem, "erkStep_SetOrder") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* set user-provided value, or default, depending on argument */
    if ord <= 0 {
        step_mem.q = Q_DEFAULT;
    } else {
        step_mem.q = ord;
    }

    /* clear Butcher tables, since user is requesting a change in method
    or a reset to defaults.  Tables will be set in ARKInitialSetup. */
    step_mem.stages = 0;
    step_mem.p = 0;

    let (mut bliw, mut blrw) = (0i64, 0i64);
    if let Some(old) = step_mem.B.take() {
        ARKodeButcherTable_Space(&old, &mut bliw, &mut blrw);
    }
    step_mem.B = None;
    ark_mem.liw -= bliw;
    ark_mem.lrw -= blrw;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  erkStep_GetEstLocalErrors: Returns the current local truncation
  error estimate vector
  ---------------------------------------------------------------*/
pub fn erkStep_GetEstLocalErrors(ark_mem: &mut ARKodeMem, ele: &mut NVector) -> i32 {
    let step_mem = match erkStep_AccessStepMem(ark_mem, "erkStep_GetEstLocalErrors") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* return an error if local truncation error is not computed */
    if (ark_mem.fixedstep && (ark_mem.AccumErrorType == ARK_ACCUMERROR_NONE))
        || (step_mem.p <= 0)
    {
        ark_mem.step_mem = Some(step_mem);
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* otherwise, copy local truncation error vector to output */
    N_VScale(ONE, &ark_mem.tempv1, ele);
    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  erkStep_GetStageIndex:

  Returns the current stage index and number of stages
  ---------------------------------------------------------------*/
pub fn erkStep_GetStageIndex(ark_mem: &mut ARKodeMem, stage: &mut i32, max_stages: &mut i32) -> i32 {
    let step_mem = match erkStep_AccessStepMem(ark_mem, "erkStep_GetStageIndex") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    *stage = step_mem.istage;
    *max_stages = step_mem.stages;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  erkStep_PrintAllStats:

  Prints integrator statistics
  ---------------------------------------------------------------*/
pub fn erkStep_PrintAllStats(
    ark_mem: &mut ARKodeMem,
    outfile: &mut dyn std::io::Write,
    fmt: SUNOutputFormat,
) -> i32 {
    let step_mem = match erkStep_AccessStepMem(ark_mem, "erkStep_PrintAllStats") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* sunfprintf_long(outfile, fmt, SUNFALSE, "RHS fn evals", nfe) */
    if fmt == SUN_OUTPUTFORMAT_TABLE {
        let _ = write!(outfile, "{:<width$} = {}\n", "RHS fn evals", step_mem.nfe, width = 29);
    } else {
        let _ = write!(outfile, ",RHS fn evals,{}", step_mem.nfe);
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  erkStep_WriteParameters:

  Outputs all solver parameters to the provided file pointer.
  ---------------------------------------------------------------*/
pub fn erkStep_WriteParameters(ark_mem: &mut ARKodeMem, fp: &mut dyn std::io::Write) -> i32 {
    let step_mem = match erkStep_AccessStepMem(ark_mem, "erkStep_WriteParameters") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* print integrator parameters to file */
    let _ = write!(fp, "ERKStep time step module parameters:\n");
    let _ = write!(fp, "  Method order {}\n", step_mem.q);
    let _ = write!(fp, "\n");

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  erkStep_SetOptions:

  Provides command-line control over ERKStep-specific "set"
  routines (arkode_erkstep_io.c).
  ---------------------------------------------------------------*/
pub fn erkStep_SetOptions(
    ark_mem: &mut ARKodeMem,
    argidx: &mut usize,
    argv: &[String],
    offset: usize,
    arg_used: &mut bool,
) -> i32 {
    use crate::sundials_cli::{sunCheckAndSetCharArgs, sunKeyCharPair};

    /* Set lists of keys, and the corresponding set routines */
    let char_pairs: [sunKeyCharPair<ARKodeMem>; 1] = [sunKeyCharPair {
        key: "table_name",
        set: ERKStepSetTableName,
    }];

    /* check all "char" keys */
    let mut j: usize = 0;
    let retval =
        sunCheckAndSetCharArgs(ark_mem, argidx, argv, offset, &char_pairs, arg_used, &mut j);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            retval,
            line!(),
            "erkStep_SetOptions",
            file!(),
            &format!("error setting key: {}", char_pairs[j].key),
        );
        return retval;
    }

    ARK_SUCCESS
}
