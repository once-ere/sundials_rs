/* -----------------------------------------------------------------
 * Translation of src/arkode/arkode_splittingstep_impl.h
 * (SUNDIALS 7.7.0).
 *
 * Implementation header for ARKODE's operator splitting module.
 *
 * C's `SUNStepper* steppers` (caller-owned pointer array copied by
 * value) and `long int* n_stepper_evolves` become owning Vecs; see
 * arkode_sunstepper.rs for the SUNStepper ownership adaptation (each
 * stepper owns its wrapped inner integrator).
 * -----------------------------------------------------------------*/

use crate::arkode_splittingstep_coefficients::SplittingStepCoefficients;
use crate::sundials_stepper::SUNStepper;

/// struct ARKodeSplittingStepMemRec
pub struct ARKodeSplittingStepMem {
    pub steppers: Vec<SUNStepper>,
    pub coefficients: Option<SplittingStepCoefficients>,
    pub n_stepper_evolves: Vec<i64>,

    pub istage: i32,
    pub partitions: i32,
    pub order: i32,
}
