/* -----------------------------------------------------------------
 * Translated from src/sunnonlinsol/newton/sunnonlinsol_newton.c
 * (SUNDIALS 7.7.0). State/counters of the Newton nonlinear solver;
 * the solve loop (SUNNonlinSolSolve_Newton) is driven from
 * cvode_nls.rs with the CVODE-specific Sys/LSetup/LSolve/CTest
 * callbacks inlined, exactly mirroring the C control flow.
 * -----------------------------------------------------------------*/
use crate::nvector_serial::NVector;
use crate::sundials_context::SUNContext;
use crate::sundials_nonlinearsolver::NonlinearSolver;

pub struct NewtonSolver {
    /// is Jacobian information current?
    pub jcur: bool,
    /// current iteration count within this solve
    pub curiter: i32,
    /// maximum Newton iterations per solve attempt (default 3)
    pub maxiters: i32,
    /// total iterations in the last solve
    pub niters: i64,
    /// total convergence failures across the solver lifetime
    pub nconvfails: i64,
    /// Newton update workspace
    pub delta: NVector,
}

/// SUNNonlinSol_Newton
pub fn SUNNonlinSol_Newton(y: &NVector, _sunctx: &SUNContext) -> NonlinearSolver {
    NonlinearSolver::Newton(NewtonSolver {
        jcur: false,
        curiter: 0,
        maxiters: 3,
        niters: 0,
        nconvfails: 0,
        delta: NVector::new(y.len()),
    })
}
