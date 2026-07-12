/* -----------------------------------------------------------------
 * Translation of src/arkode/arkode_sprkstep_impl.h (+ the
 * SPRKSTEP_DEFAULT_* constants of include/arkode/arkode_sprkstep.h).
 * SPRKStep time-stepper module memory structure.
 * -----------------------------------------------------------------*/

use crate::arkode_impl::ARKRhsFn;
use crate::arkode_sprk::{
    ARKodeSPRKTable, ARKODE_SPRKMethodID, ARKODE_SPRK_EULER_1_1, ARKODE_SPRK_LEAPFROG_2_2,
    ARKODE_SPRK_MCLACHLAN_3_3, ARKODE_SPRK_MCLACHLAN_4_4, ARKODE_SPRK_MCLACHLAN_5_6,
    ARKODE_SPRK_SOFRONIOU_10_36, ARKODE_SPRK_SUZUKI_UMENO_8_16, ARKODE_SPRK_YOSHIDA_6_8,
};
use crate::nvector_serial::NVector;

/* Default SPRK tables per order (arkode_sprkstep.h) */
pub const SPRKSTEP_DEFAULT_1: ARKODE_SPRKMethodID = ARKODE_SPRK_EULER_1_1;
pub const SPRKSTEP_DEFAULT_2: ARKODE_SPRKMethodID = ARKODE_SPRK_LEAPFROG_2_2;
pub const SPRKSTEP_DEFAULT_3: ARKODE_SPRKMethodID = ARKODE_SPRK_MCLACHLAN_3_3;
pub const SPRKSTEP_DEFAULT_4: ARKODE_SPRKMethodID = ARKODE_SPRK_MCLACHLAN_4_4;
pub const SPRKSTEP_DEFAULT_5: ARKODE_SPRKMethodID = ARKODE_SPRK_MCLACHLAN_5_6;
pub const SPRKSTEP_DEFAULT_6: ARKODE_SPRKMethodID = ARKODE_SPRK_YOSHIDA_6_8;
pub const SPRKSTEP_DEFAULT_8: ARKODE_SPRKMethodID = ARKODE_SPRK_SUZUKI_UMENO_8_16;
pub const SPRKSTEP_DEFAULT_10: ARKODE_SPRKMethodID = ARKODE_SPRK_SOFRONIOU_10_36;

/// struct ARKodeSPRKStepMemRec (arkode_sprkstep_impl.h)
#[derive(Default)]
pub struct ARKodeSPRKStepMem {
    /* SPRK method and storage */
    pub method: Option<ARKodeSPRKTable>, /* method spec  */
    pub q: i32,                          /* method order */
    pub sdata: NVector,                  /* persisted stage data */
    pub yerr: NVector,                   /* error vector for compensated
                                         summation (empty = C NULL) */

    /* SPRK problem specification */
    pub f1: Option<ARKRhsFn>, /* p' = f1(t,q) = - dV(t,q)/dq  */
    pub f2: Option<ARKRhsFn>, /* q' = f2(t,p) =   dT(t,p)/dp  */

    /* Counters */
    pub nf1: i64, /* number of calls to f1        */
    pub nf2: i64, /* number of calls to f2        */
    pub istage: i32,
}

/* Initialization and I/O error messages */
pub const MSG_SPRKSTEP_NO_MEM: &str = "Time step module memory is NULL.";
