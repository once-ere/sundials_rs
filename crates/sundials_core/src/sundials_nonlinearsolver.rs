/* -----------------------------------------------------------------
 * Translated from src/sundials/sundials_nonlinearsolver.c and
 * include/sundials/sundials_nonlinearsolver.h (SUNDIALS 7.7.0).
 *
 * The generic SUNNonlinearSolver ops-table becomes an enum over the
 * two concrete solvers CVODE can use (Newton for CV_BDF/CV_ADAMS with
 * a linear solver; accelerated fixed point otherwise). The solve
 * loops themselves are driven from cvode_nls.rs, where the C code
 * installs its Sys/LSetup/LSolve/CTest callbacks.
 * -----------------------------------------------------------------*/
use crate::sunnonlinsol_fixedpoint::FixedPointSolver;
use crate::sunnonlinsol_newton::NewtonSolver;

/* SUNNonlinearSolver_Type */
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SUNNonlinearSolver_Type {
    SUNNONLINEARSOLVER_ROOTFIND,
    SUNNONLINEARSOLVER_FIXEDPOINT,
}
pub use SUNNonlinearSolver_Type::*;

/* Return flags (sundials_nonlinearsolver.h) */
pub const SUN_NLS_CONTINUE: i32 = 901; /* not converged, keep iterating      */
pub const SUN_NLS_CONV_RECVR: i32 = 902; /* convergence failure, try to recover */

pub enum NonlinearSolver {
    Newton(NewtonSolver),
    FixedPoint(FixedPointSolver),
}

impl NonlinearSolver {
    /// SUNNonlinSolGetType
    pub fn nls_type(&self) -> SUNNonlinearSolver_Type {
        match self {
            NonlinearSolver::Newton(_) => SUNNONLINEARSOLVER_ROOTFIND,
            NonlinearSolver::FixedPoint(_) => SUNNONLINEARSOLVER_FIXEDPOINT,
        }
    }

    /// SUNNonlinSolInitialize
    pub fn initialize(&mut self) -> i32 {
        match self {
            NonlinearSolver::Newton(s) => {
                s.niters = 0;
                s.nconvfails = 0;
                s.jcur = false;
                0
            }
            NonlinearSolver::FixedPoint(s) => {
                s.niters = 0;
                s.nconvfails = 0;
                0
            }
        }
    }

    /// SUNNonlinSolSetMaxIters
    pub fn set_max_iters(&mut self, maxiters: i32) -> i32 {
        if maxiters < 1 {
            return crate::sundials_errors::SUN_ERR_ARG_OUTOFRANGE;
        }
        match self {
            NonlinearSolver::Newton(s) => s.maxiters = maxiters,
            NonlinearSolver::FixedPoint(s) => s.maxiters = maxiters,
        }
        0
    }

    /// SUNNonlinSolGetNumIters
    pub fn get_num_iters(&self) -> i64 {
        match self {
            NonlinearSolver::Newton(s) => s.niters,
            NonlinearSolver::FixedPoint(s) => s.niters,
        }
    }

    /// SUNNonlinSolGetCurIter
    pub fn get_cur_iter(&self) -> i32 {
        match self {
            NonlinearSolver::Newton(s) => s.curiter,
            NonlinearSolver::FixedPoint(s) => s.curiter,
        }
    }

    /// SUNNonlinSolGetNumConvFails
    pub fn get_num_conv_fails(&self) -> i64 {
        match self {
            NonlinearSolver::Newton(s) => s.nconvfails,
            NonlinearSolver::FixedPoint(s) => s.nconvfails,
        }
    }

    /// SUNNonlinSolSetDamping (fixed point only)
    pub fn set_damping(&mut self, beta: f64) -> i32 {
        match self {
            NonlinearSolver::FixedPoint(s) => s.set_damping(beta),
            _ => crate::sundials_errors::SUN_ERR_NOT_IMPLEMENTED,
        }
    }
}
