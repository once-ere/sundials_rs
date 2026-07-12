/* -----------------------------------------------------------------
 * Translated from src/arkode/arkode_bbdpre_impl.h and the type
 * definitions of include/arkode/arkode_bbdpre.h (ARKODE 7.7.0).
 * ARKBBDPRE band-block-diagonal preconditioner module data
 * (serial reduction: a single block of dimension n_local).
 *
 * C keeps a `void* arkode_mem` back-pointer in the pdata block; the
 * Rust setup/solve routines receive `&mut ARKodeMem` as an argument
 * instead (the module lives inside ARKLsMem.prec_module, which is
 * detached from ARKodeMem for the duration of every ARKLS call).
 * -----------------------------------------------------------------*/
use crate::nvector_serial::NVector;
use crate::sundials_linearsolver::LinearSolver;
use crate::sundials_matrix::SUNMatrix;
use crate::sundials_types::UserData;

/* ARKLocalFn (arkode_bbdpre.h): computes g(t,y) approximating the
   right-hand side function f */
pub type ARKLocalFn =
    fn(nlocal: i64, t: f64, y: &NVector, g: &mut NVector, user_data: &mut UserData) -> i32;

/* ARKCommFn (arkode_bbdpre.h): inter-process communication — a
   serial build never communicates, but the type is preserved. */
pub type ARKCommFn = fn(nlocal: i64, t: f64, y: &NVector, user_data: &mut UserData) -> i32;

/* -----------------------------------------------------------------
   Type: ARKBBDPrecData
   -----------------------------------------------------------------*/
pub struct ARKBBDPrecData {
    /* passed by user to ARKBBDPrecAlloc and used by PrecSetup/PrecSolve */
    pub mudq: i64,
    pub mldq: i64,
    pub mukeep: i64,
    pub mlkeep: i64,
    pub dqrely: f64,
    pub gloc: Option<ARKLocalFn>,
    pub cfn: Option<ARKCommFn>,

    /* set by ARKBBDPrecSetup and used by ARKBBDPrecSolve */
    pub savedJ: SUNMatrix,
    pub savedP: SUNMatrix,
    pub LS: LinearSolver,
    pub tmp1: NVector,
    pub tmp2: NVector,
    pub tmp3: NVector,
    pub zlocal: NVector,
    pub rlocal: NVector,

    /* set by ARKBBDPrecAlloc and used by ARKBBDPrecSetup */
    pub n_local: i64,

    /* available for optional output */
    pub rpwsize: i64,
    pub ipwsize: i64,
    pub nge: i64,
}

/* -----------------------------------------------------------------
   ARKBBDPRE error messages
   -----------------------------------------------------------------*/
pub const MSG_BBD_MEM_FAIL: &str = "A memory request failed.";
pub const MSG_BBD_BAD_NVECTOR: &str = "A required vector operation is not implemented.";
pub const MSG_BBD_SUNMAT_FAIL: &str = "An error arose from a SUNBandMatrix routine.";
pub const MSG_BBD_SUNLS_FAIL: &str = "An error arose from a SUNBandLinearSolver routine.";
pub const MSG_BBD_PMEM_NULL: &str =
    "BBD peconditioner memory is NULL. ARKBBDPrecInit must be called.";
pub const MSG_BBD_FUNC_FAILED: &str =
    "The gloc or cfn routine failed in an unrecoverable manner.";
