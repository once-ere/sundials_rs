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
use crate::sundials_nvector_senswrapper::{N_VNew_SensWrapper, NVectorSensWrapper};

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
    /// Newton update workspace of a solver built by SUNNonlinSol_NewtonSens.
    /// In C the Sens constructor clones a SensWrapper template, so the ONE
    /// `delta` workspace *is* a senswrapper holding `count` sub-vectors; in
    /// this port that wrapper lives here and `delta` stays empty (the plain
    /// constructor leaves `deltaS` empty instead). The CVODES sensitivity
    /// correctors (cvodes_nls_sim / cvodes_nls_stg) drive the iteration on
    /// `deltaS`.
    pub deltaS: NVectorSensWrapper,
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
        deltaS: NVectorSensWrapper::default(),
    })
}

/// SUNNonlinSol_NewtonSens (sunnonlinsol_newton.c): constructor wrapper to
/// create a new Newton solver for the CVODES/IDAS sensitivity correctors.
/// In C this builds a temporary senswrapper w = N_VNew_SensWrapper(count, y)
/// (i.e. `count` clones of the template `y`: Ns+1 sub-vectors for the
/// SIMULTANEOUS corrector — state + Ns sensitivities — or Ns for STAGGERED),
/// calls SUNNonlinSol_Newton(w, sunctx) — whose N_VClone(w) makes the solver's
/// `delta` workspace a senswrapper of the same shape — and destroys w. Here
/// the cloned wrapper is created directly in `deltaS` and `delta` stays
/// empty. Panics if count < 1 (C would return NULL from N_VNew_SensWrapper
/// and crash in N_VClone).
pub fn SUNNonlinSol_NewtonSens(count: i32, y: &NVector, _sunctx: &SUNContext) -> NonlinearSolver {
    NonlinearSolver::Newton(NewtonSolver {
        jcur: false,
        curiter: 0,
        maxiters: 3,
        niters: 0,
        nconvfails: 0,
        delta: NVector::default(),
        deltaS: N_VNew_SensWrapper(count, y).expect("SUNNonlinSol_NewtonSens: count < 1"),
    })
}
