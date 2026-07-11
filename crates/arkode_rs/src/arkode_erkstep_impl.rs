/* -----------------------------------------------------------------
 * Translation of sundials-7.7.0/src/arkode/arkode_erkstep_impl.h
 * (+ the ERKSTEP_DEFAULT_* constants of include/arkode/
 * arkode_erkstep.h).
 *
 * Modeling notes:
 *  - `N_Vector* F` -> Vec<NVector> (empty = unallocated).
 *  - `Xvecs` (a scratch array of N_Vector POINTERS for the fused
 *    ops) cannot be stored in safe Rust; the operand list is
 *    assembled at each linear-combination call site instead. The
 *    liw accounting for it is kept (nfusedopvecs).
 *  - `forcing` in C aliases the caller's vectors (MRIStep); here it
 *    is an owned copy refreshed on every erkStep_SetInnerForcing
 *    call (MRIStep re-sets the forcing before each fast integration).
 *  - `adj_f` (SUNAdjRhsFn, discrete-adjoint RHS) is deferred with
 *    the erkStep_TakeStep_Adjoint machinery (needs the ManyVector
 *    module).
 * -----------------------------------------------------------------*/

use crate::arkode_butcher::ARKodeButcherTable;
use crate::arkode_butcher_erk::{
    ARKODE_BOGACKI_SHAMPINE_4_2_3, ARKODE_FORWARD_EULER_1_1, ARKODE_RALSTON_3_1_2,
    ARKODE_SOFRONIOU_SPALETTA_5_3_4, ARKODE_TSITOURAS_7_4_5, ARKODE_VERNER_10_6_7,
    ARKODE_VERNER_13_7_8, ARKODE_VERNER_16_8_9, ARKODE_VERNER_9_5_6, ARKODE_ERKTableID,
};
use crate::arkode_impl::ARKRhsFn;
use crate::nvector_serial::NVector;

/* Default Butcher tables per order (arkode_erkstep.h) */
pub const ERKSTEP_DEFAULT_1: ARKODE_ERKTableID = ARKODE_FORWARD_EULER_1_1;
pub const ERKSTEP_DEFAULT_2: ARKODE_ERKTableID = ARKODE_RALSTON_3_1_2;
pub const ERKSTEP_DEFAULT_3: ARKODE_ERKTableID = ARKODE_BOGACKI_SHAMPINE_4_2_3;
pub const ERKSTEP_DEFAULT_4: ARKODE_ERKTableID = ARKODE_SOFRONIOU_SPALETTA_5_3_4;
pub const ERKSTEP_DEFAULT_5: ARKODE_ERKTableID = ARKODE_TSITOURAS_7_4_5;
pub const ERKSTEP_DEFAULT_6: ARKODE_ERKTableID = ARKODE_VERNER_9_5_6;
pub const ERKSTEP_DEFAULT_7: ARKODE_ERKTableID = ARKODE_VERNER_10_6_7;
pub const ERKSTEP_DEFAULT_8: ARKODE_ERKTableID = ARKODE_VERNER_13_7_8;
pub const ERKSTEP_DEFAULT_9: ARKODE_ERKTableID = ARKODE_VERNER_16_8_9;

pub const MSG_ERKSTEP_NO_MEM: &str = "Time step module memory is NULL.";

/// struct ARKodeERKStepMemRec (arkode_erkstep_impl.h)
pub struct ARKodeERKStepMem {
    /* ERK problem specification */
    pub f: Option<ARKRhsFn>, /* y' = f(t,y) */

    /* ARK method storage and parameters */
    pub F: Vec<NVector>, /* explicit RHS at each stage */
    pub q: i32,          /* method order               */
    pub p: i32,          /* embedding order            */
    pub istage: i32,     /* current stage              */
    pub stages: i32,     /* number of stages           */
    pub B: Option<ARKodeButcherTable>, /* ERK Butcher table */

    /* Counters */
    pub nfe: i64, /* num fe calls               */

    /* Reusable arrays for fused vector operations */
    pub cvals: Vec<f64>,
    pub nfusedopvecs: i32, /* length of cvals and Xvecs arrays */

    /* Data for using ERKStep with external polynomial forcing */
    pub tshift: f64,           /* time normalization shift       */
    pub tscale: f64,           /* time normalization scaling     */
    pub forcing: Vec<NVector>, /* array of forcing vectors       */
    pub nforcing: i32,         /* number of forcing vectors      */
    pub stage_times: Vec<f64>, /* workspace for applying forcing */
    pub stage_coefs: Vec<f64>, /* workspace for applying forcing */
}
