/* -----------------------------------------------------------------
 * Translated from src/kinsol/kinsol_bbdpre_impl.h and the type
 * definitions of include/kinsol/kinsol_bbdpre.h (KINSOL 7.7.0).
 * KINBBDPRE band-block-diagonal preconditioner module data
 * (serial reduction: a single block). Conventions follow the donor
 * cvode_bbdpre_impl.rs.
 * -----------------------------------------------------------------*/
use crate::nvector_serial::NVector;
use crate::sundials_linearsolver::LinearSolver;
use crate::sundials_matrix::SUNMatrix;
use crate::sundials_types::UserData;

/* KINBBDCommFn (kinsol_bbdpre.h): inter-process communication — a
   serial build never communicates, but the type is preserved. */
pub type KINBBDCommFn = fn(Nlocal: i64, u: &NVector, user_data: &mut UserData) -> i32;

/* KINBBDLocalFn (kinsol_bbdpre.h): computes g(u) approximating f */
pub type KINBBDLocalFn =
    fn(Nlocal: i64, uu: &NVector, gval: &mut NVector, user_data: &mut UserData) -> i32;

/*------------------------------------------------------------------
  Definition of KBBDData
  ------------------------------------------------------------------*/
pub struct KBBDPrecData {
    /* passed by user to KINBBDPrecAlloc, used by pset/psolve functions */
    pub mudq: i64,
    pub mldq: i64,
    pub mukeep: i64,
    pub mlkeep: i64,
    pub rel_uu: f64, /* relative error for the Jacobian DQ routine */
    pub gloc: Option<KINBBDLocalFn>,
    pub gcomm: Option<KINBBDCommFn>,

    /* set by KINBBDPrecSetup and used by KINBBDPrecSetup and
       KINBBDPrecSolve functions */
    pub n_local: i64,
    pub PP: SUNMatrix,
    pub LS: LinearSolver,
    pub rlocal: NVector,
    pub zlocal: NVector,
    pub tempv1: NVector,
    pub tempv2: NVector,
    pub tempv3: NVector,

    /* available for optional output */
    pub rpwsize: i64,
    pub ipwsize: i64,
    pub nge: i64,
    /* (C `void* kin_mem` back-pointer: the solver memory is passed to
       the module routines as a &mut KINMem argument instead.) */
}
