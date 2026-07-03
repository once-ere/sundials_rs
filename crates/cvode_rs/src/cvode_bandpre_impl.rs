/* -----------------------------------------------------------------
 * Translated from src/cvode/cvode_bandpre_impl.h (CVODE 7.7.0).
 * CVBANDPRE banded preconditioner module data.
 * -----------------------------------------------------------------*/
use crate::nvector_serial::NVector;
use crate::sundials_linearsolver::LinearSolver;
use crate::sundials_matrix::SUNMatrix;

/* -----------------------------------------------------------------
   Type: CVBandPrecData
   -----------------------------------------------------------------*/
pub struct CVBandPrecData {
    /* Data set by user in CVBandPrecInit */
    pub N: i64,
    pub ml: i64,
    pub mu: i64,

    /* Data set by CVBandPrecSetup */
    pub savedJ: SUNMatrix,
    pub savedP: SUNMatrix,
    pub LS: LinearSolver,
    pub tmp1: NVector,
    pub tmp2: NVector,

    /* Rhs calls */
    pub nfeBP: i64,
}
