/* -----------------------------------------------------------------
 * Translated from src/idas/idas_bbdpre_impl.h and the type
 * definitions of include/idas/idas_bbdpre.h (IDAS 7.7.0).
 * IDASBBDPRE band-block-diagonal preconditioner module data
 * (serial reduction: a single block). Forward-problem conventions
 * follow the verified ida_bbdpre_impl.rs donor; the backward-problem
 * additions (IDABBDPrecDataB + the *FnB callback types) mirror the
 * idas_ls_impl.rs PART II modeling.
 * -----------------------------------------------------------------*/
use crate::nvector_serial::NVector;
use crate::sundials_linearsolver::LinearSolver;
use crate::sundials_matrix::SUNMatrix;
use crate::sundials_types::UserData;

/* IDABBDLocalFn (idas_bbdpre.h): computes G(t,y,y') approximating the
   residual function F */
pub type IDABBDLocalFn = fn(
    Nlocal: i64,
    tt: f64,
    yy: &NVector,
    yp: &NVector,
    gval: &mut NVector,
    user_data: &mut UserData,
) -> i32;

/* IDABBDCommFn (idas_bbdpre.h): inter-process communication — a
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

/*------------------------------------------------------------------
  Backward-problem types (idas_bbdpre.h BACKWARD PROBLEMS section)
  ------------------------------------------------------------------*/

/* IDABBDLocalFnB: backward G^B(t, y, y', yB, yB') local residual
   approximation; receives the interpolated forward solution (yy, yp). */
pub type IDABBDLocalFnB = fn(
    NlocalB: i64,
    tt: f64,
    yy: &NVector,
    yp: &NVector,
    yyB: &NVector,
    ypB: &NVector,
    gvalB: &mut NVector,
    user_dataB: &mut UserData,
) -> i32;

/* IDABBDCommFnB: backward inter-process communication (serial no-op,
   type preserved). */
pub type IDABBDCommFnB = fn(
    NlocalB: i64,
    tt: f64,
    yy: &NVector,
    yp: &NVector,
    yyB: &NVector,
    ypB: &NVector,
    user_dataB: &mut UserData,
) -> i32;

/*------------------------------------------------------------------
  Type: IDABBDPrecDataB (idas_bbdpre_impl.h)

  Stored behind IDAB_mem.ida_pmem (Option<Box<dyn Any>>, downcast at
  the wrapper call sites); the C ida_pfree hook is Rust Drop.
  ------------------------------------------------------------------*/
pub struct IDABBDPrecDataB {
    /* BBD user functions (glocB and cfnB) for backward run */
    pub glocalB: Option<IDABBDLocalFnB>,
    pub gcommB: Option<IDABBDCommFnB>,
}

/* (C idas_bbdpre_impl.h also defines MSGBBD_AMEM_NULL and
   MSGBBD_PDATAB_NULL; they are unused in the 7.7.0 sources and are
   not carried.  The used messages live in idas_bbdpre.rs, donor
   convention.) */
