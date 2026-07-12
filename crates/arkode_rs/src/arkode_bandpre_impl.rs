/* -----------------------------------------------------------------
 * Translated from src/arkode/arkode_bandpre_impl.h (ARKODE 7.7.0).
 * ARKBANDPRE banded preconditioner module data.
 *
 * C keeps a `void* arkode_mem` back-pointer in the pdata block; the
 * Rust setup/solve routines receive `&mut ARKodeMem` as an argument
 * instead (the module lives inside ARKLsMem.prec_module, which is
 * detached from ARKodeMem for the duration of every ARKLS call).
 * -----------------------------------------------------------------*/
use crate::nvector_serial::NVector;
use crate::sundials_linearsolver::LinearSolver;
use crate::sundials_matrix::SUNMatrix;

/* -----------------------------------------------------------------
   Type: ARKBandPrecData
   -----------------------------------------------------------------*/
pub struct ARKBandPrecData {
    /* Data set by user in ARKBandPrecInit */
    pub N: i64,
    pub ml: i64,
    pub mu: i64,

    /* Data set by ARKBandPrecSetup */
    pub savedJ: SUNMatrix,
    pub savedP: SUNMatrix,
    pub LS: LinearSolver,
    pub tmp1: NVector,
    pub tmp2: NVector,

    /* Rhs calls */
    pub nfeBP: i64,
}

/* -----------------------------------------------------------------
   ARKBANDPRE error messages
   -----------------------------------------------------------------*/
pub const MSG_BP_MEM_FAIL: &str = "A memory request failed.";
pub const MSG_BP_BAD_NVECTOR: &str = "A required vector operation is not implemented.";
pub const MSG_BP_SUNMAT_FAIL: &str = "An error arose from a SUNBandMatrix routine.";
pub const MSG_BP_SUNLS_FAIL: &str = "An error arose from a SUNBandLinearSolver routine.";
pub const MSG_BP_PMEM_NULL: &str =
    "Band preconditioner memory is NULL. ARKBandPrecInit must be called.";
pub const MSG_BP_RHSFUNC_FAILED: &str =
    "The right-hand side routine failed in an unrecoverable manner.";
