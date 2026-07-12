/* -----------------------------------------------------------------
 * Translation of src/arkode/arkode_forcingstep_impl.h
 * (SUNDIALS 7.7.0).
 *
 * Implementation header for ARKODE's forcing method.
 *
 * The two SUNSteppers are owned by the step memory (see
 * arkode_sunstepper.rs for the ownership adaptation).
 * -----------------------------------------------------------------*/

use crate::sundials_stepper::SUNStepper;

pub const NUM_PARTITIONS: usize = 2;

/// struct ARKodeForcingStepMemRec
pub struct ARKodeForcingStepMem {
    pub stepper: [SUNStepper; NUM_PARTITIONS],
    pub n_stepper_evolves: [i64; NUM_PARTITIONS],
}
