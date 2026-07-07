/* -----------------------------------------------------------------
 * Translated from src/ida/ida_bbdpre_impl.h and the type
 * definitions of include/ida/ida_bbdpre.h (IDA 7.7.0).
 * IDABBDPRE band-block-diagonal preconditioner module data
 * (serial reduction: a single block). Conventions follow the landed
 * kinsol_bbdpre_impl.rs.
 * -----------------------------------------------------------------*/
use crate::nvector_serial::NVector;
use crate::sundials_linearsolver::LinearSolver;
use crate::sundials_matrix::SUNMatrix;
use crate::sundials_types::UserData;

/* IDABBDLocalFn (ida_bbdpre.h): computes G(t,y,y') approximating the
   residual function F */
pub type IDABBDLocalFn = fn(
    Nlocal: i64,
    tt: f64,
    yy: &NVector,
    yp: &NVector,
    gval: &mut NVector,
    user_data: &mut UserData,
) -> i32;

/* IDABBDCommFn (ida_bbdpre.h): inter-process communication — a
   serial build never communicates, but the type is preserved. */
pub type IDABBDCommFn =
    fn(Nlocal: i64, tt: f64, yy: &NVector, yp: &NVector, user_data: &mut UserData) -> i32;

/*------------------------------------------------------------------
  Definition of IBBDPrecData
  ------------------------------------------------------------------*/
pub struct IBBDPrecData {
    /* passed by user to IDABBDPrecAlloc and used by
       IDABBDPrecSetup/IDABBDPrecSolve functions */
    pub mudq: i64,
    pub mldq: i64,
    pub mukeep: i64,
    pub mlkeep: i64,
    pub rel_yy: f64,
    pub glocal: Option<IDABBDLocalFn>,
    pub gcomm: Option<IDABBDCommFn>,

    /* set by IDABBDPrecSetup and used by IDABBDPrecSetup and
       IDABBDPrecSolve functions */
    pub n_local: i64,
    pub PP: SUNMatrix,
    pub LS: LinearSolver,
    pub zlocal: NVector,
    pub rlocal: NVector,
    pub tempv1: NVector,
    pub tempv2: NVector,
    pub tempv3: NVector,
    pub tempv4: NVector,

    /* available for optional output */
    pub rpwsize: i64,
    pub ipwsize: i64,
    pub nge: i64,
    /* (C `void* ida_mem` back-pointer: the solver memory is passed to
       the module routines as a &mut IDAMem argument instead.) */
}
