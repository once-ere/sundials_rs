/* -----------------------------------------------------------------
 * Translated from src/cvode/cvode_bbdpre_impl.h (CVODE 7.7.0).
 * CVBBDPRE band-block-diagonal preconditioner module data
 * (serial reduction: a single block).
 * -----------------------------------------------------------------*/
use crate::nvector_serial::NVector;
use crate::sundials_linearsolver::LinearSolver;
use crate::sundials_matrix::SUNMatrix;
use crate::sundials_types::UserData;

/* CVLocalFn (cvode_bbdpre.h): computes g(t,y) approximating f */
pub type CVLocalFn =
    fn(nlocal: i64, t: f64, y: &NVector, g: &mut NVector, user_data: &mut UserData) -> i32;

/* CVCommFn (cvode_bbdpre.h): inter-process communication — a serial
   build never communicates, but the type is preserved. */
pub type CVCommFn = fn(nlocal: i64, t: f64, y: &NVector, user_data: &mut UserData) -> i32;

/* -----------------------------------------------------------------
   Type: CVBBDPrecData
   -----------------------------------------------------------------*/
pub struct CVBBDPrecData {
    /* passed by user to CVBBDPrecInit and used by PrecSetup/PrecSolve */
    pub mudq: i64,
    pub mldq: i64,
    pub mukeep: i64,
    pub mlkeep: i64,
    pub dqrely: f64,
    pub gloc: Option<CVLocalFn>,
    pub cfn: Option<CVCommFn>,

    /* set by CVBBDPrecSetup and used by CVBBDPrecSolve */
    pub savedJ: SUNMatrix,
    pub savedP: SUNMatrix,
    pub LS: LinearSolver,
    pub tmp1: NVector,
    pub tmp2: NVector,
    pub tmp3: NVector,
    pub zlocal: NVector,
    pub rlocal: NVector,

    /* set by CVBBDPrecInit and used by CVBBDPrecSetup */
    pub n_local: i64,

    /* available for optional output */
    pub rpwsize: i64,
    pub ipwsize: i64,
    pub nge: i64,
}
