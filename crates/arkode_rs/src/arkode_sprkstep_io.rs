/* -----------------------------------------------------------------
 * Translation of src/arkode/arkode_sprkstep_io.c (ARKODE 7.7.0).
 * SPRKStep optional input/output functions.
 *
 * Not ported: sprkStep_SetOptions (CLI module) and the deprecated
 * SPRKStep* aliases for shared ARKODE routines (call the ARKode*
 * versions directly).
 * -----------------------------------------------------------------*/

use crate::arkode::arkAllocVec;
use crate::arkode_impl::*;
use crate::arkode_sprk::{
    ARKodeSPRKTable, ARKodeSPRKTable_Copy, ARKodeSPRKTable_Free, ARKodeSPRKTable_LoadByName,
};
use crate::arkode_sprkstep::{
    sprkStep_AccessStepMem, sprkStep_TakeStep, sprkStep_TakeStep_Compensated,
};
use crate::arkode_io::sunfprintf_long;
use crate::nvector_serial::N_VConst;
use crate::sundials_types::{SUNOutputFormat, SUNFALSE};

/*===============================================================
  Exported optional input functions.
  ===============================================================*/

/*---------------------------------------------------------------
  SPRKStepSetUseCompensatedSums:

  Turns on/off compensated summation in SPRKStep and ARKODE.
  ---------------------------------------------------------------*/
pub fn SPRKStepSetUseCompensatedSums(ark_mem: &mut ARKodeMem, onoff: bool) -> i32 {
    ark_mem.use_compensated_sums = onoff;
    sprkStep_SetUseCompensatedSums(ark_mem, onoff)
}

/*---------------------------------------------------------------
  SPRKStepSetMethod:

  Specifies the SPRK method

  ** Note in documentation that this should not be called along
  with ARKodeSetOrder. **
  ---------------------------------------------------------------*/
pub fn SPRKStepSetMethod(ark_mem: &mut ARKodeMem, sprk_storage: &ARKodeSPRKTable) -> i32 {
    let mut step_mem = match sprkStep_AccessStepMem(ark_mem, "SPRKStepSetMethod") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    if let Some(method) = step_mem.method.take() {
        ARKodeSPRKTable_Free(method);
    }

    step_mem.method = Some(ARKodeSPRKTable_Copy(sprk_storage));

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  SPRKStepSetMethodName:

  Specifies the SPRK method.
  ---------------------------------------------------------------*/
pub fn SPRKStepSetMethodName(ark_mem: &mut ARKodeMem, method: &str) -> i32 {
    let mut step_mem = match sprkStep_AccessStepMem(ark_mem, "SPRKStepSetMethodName") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    if let Some(old) = step_mem.method.take() {
        ARKodeSPRKTable_Free(old);
    }

    step_mem.method = ARKodeSPRKTable_LoadByName(method);

    let ok = step_mem.method.is_some();
    ark_mem.step_mem = Some(step_mem);
    if ok {
        ARK_SUCCESS
    } else {
        ARK_ILL_INPUT
    }
}

/*===============================================================
  Exported optional output functions.
  ===============================================================*/

/*---------------------------------------------------------------
  SPRKStepGetCurrentMethod:

  Returns a copy of the stepper method structure (C hands out the
  internal pointer).
  ---------------------------------------------------------------*/
pub fn SPRKStepGetCurrentMethod(
    ark_mem: &mut ARKodeMem,
    sprk_storage: &mut Option<ARKodeSPRKTable>,
) -> i32 {
    let step_mem = match sprkStep_AccessStepMem(ark_mem, "SPRKStepGetCurrentMethod") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    *sprk_storage = step_mem.method.as_ref().map(ARKodeSPRKTable_Copy);

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  sprkStep_GetNumRhsEvals:

  Returns the current number of RHS calls
  ---------------------------------------------------------------*/
pub fn sprkStep_GetNumRhsEvals(
    ark_mem: &mut ARKodeMem,
    partition_index: i32,
    rhs_evals: &mut i64,
) -> i32 {
    let step_mem = match sprkStep_AccessStepMem(ark_mem, "sprkStep_GetNumRhsEvals") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    if partition_index > 1 {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "sprkStep_GetNumRhsEvals",
            file!(),
            "Invalid partition index",
        );
        return ARK_ILL_INPUT;
    }

    match partition_index {
        0 => *rhs_evals = step_mem.nf1,
        1 => *rhs_evals = step_mem.nf2,
        _ => *rhs_evals = step_mem.nf1 + step_mem.nf2,
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*===============================================================
  Private functions attached to ARKODE
  ===============================================================*/

/*---------------------------------------------------------------
  sprkStep_SetDefaults:

  Resets all SPRKStep optional inputs to their default values.
  Does not change problem-defining function pointers or
  user_data pointer.
  ---------------------------------------------------------------*/
pub fn sprkStep_SetDefaults(ark_mem: &mut ARKodeMem) -> i32 {
    /* use the default method order */
    sprkStep_SetOrder(ark_mem, 0)
}

/*---------------------------------------------------------------
  sprkStep_SetOrder:

  Specifies the method order
  ---------------------------------------------------------------*/
pub fn sprkStep_SetOrder(ark_mem: &mut ARKodeMem, ord: i32) -> i32 {
    let mut step_mem = match sprkStep_AccessStepMem(ark_mem, "sprkStep_SetOrder") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* Invalid orders result in the default order being used. */
    let mut ord = ord;
    if ord == 7 || ord == 9 || ord > 10 {
        ord = -1;
    }

    /* set user-provided value, or default, depending on argument */
    if ord <= 0 {
        step_mem.q = 4;
    } else {
        step_mem.q = ord;
    }

    if let Some(method) = step_mem.method.take() {
        ARKodeSPRKTable_Free(method);
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  sprkStep_GetStageIndex:

  Returns the current stage index and number of stages
  ---------------------------------------------------------------*/
pub fn sprkStep_GetStageIndex(ark_mem: &mut ARKodeMem, stage: &mut i32, max_stages: &mut i32) -> i32 {
    let step_mem = match sprkStep_AccessStepMem(ark_mem, "sprkStep_GetStageIndex") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* if table is not yet set, return defaults */
    if let Some(method) = &step_mem.method {
        *stage = step_mem.istage;
        *max_stages = method.stages;
    } else {
        *stage = -1;
        *max_stages = -1;
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_NULL,
            line!(),
            "sprkStep_GetStageIndex",
            file!(),
            "method structure not allocated",
        );
        /* (C returns the earlier successful access retval here) */
        return ARK_SUCCESS;
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  sprkStep_PrintAllStats:

  Prints integrator statistics
  ---------------------------------------------------------------*/
pub fn sprkStep_PrintAllStats(
    ark_mem: &mut ARKodeMem,
    outfile: &mut dyn std::io::Write,
    fmt: SUNOutputFormat,
) -> i32 {
    let step_mem = match sprkStep_AccessStepMem(ark_mem, "sprkStep_PrintAllStats") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    sunfprintf_long(outfile, fmt, SUNFALSE, "f1 RHS fn evals", step_mem.nf1);
    sunfprintf_long(outfile, fmt, SUNFALSE, "f2 RHS fn evals", step_mem.nf2);

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  sprkStep_WriteParameters:

  Outputs all solver parameters to the provided file pointer.
  ---------------------------------------------------------------*/
pub fn sprkStep_WriteParameters(ark_mem: &mut ARKodeMem, fp: &mut dyn std::io::Write) -> i32 {
    let step_mem = match sprkStep_AccessStepMem(ark_mem, "sprkStep_WriteParameters") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* print integrator parameters to file */
    let _ = write!(fp, "SPRKStep time step module parameters:\n");
    let _ = write!(fp, "  Method order {}\n", step_mem.method.as_ref().unwrap().q);
    let _ = write!(
        fp,
        "  Method stages {}\n",
        step_mem.method.as_ref().unwrap().stages
    );

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

pub fn sprkStep_SetUseCompensatedSums(ark_mem: &mut ARKodeMem, onoff: bool) -> i32 {
    let mut step_mem = match sprkStep_AccessStepMem(ark_mem, "sprkStep_SetUseCompensatedSums") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    if onoff {
        ark_mem.step = Some(sprkStep_TakeStep_Compensated);
        if step_mem.yerr.data.is_empty() {
            let tmpl_len = ark_mem.yn.data.len();
            arkAllocVec(ark_mem, tmpl_len, &mut step_mem.yerr);
            /* Zero yerr for compensated summation */
            N_VConst(ZERO, &mut step_mem.yerr);
        }
    } else {
        ark_mem.step = Some(sprkStep_TakeStep);
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}
