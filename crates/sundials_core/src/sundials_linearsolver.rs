/* -----------------------------------------------------------------
 * Translated from src/sundials/sundials_linearsolver.c and
 * include/sundials/sundials_linearsolver.h (SUNDIALS 7.7.0).
 *
 * The C generic SUNLinearSolver (ops-table over void* content)
 * becomes an enum over the concrete solver implementations; the
 * generic SUNLinSol* wrappers become the methods below.
 * -----------------------------------------------------------------*/
use crate::nvector_serial::NVector;
use crate::sundials_errors::*;
use crate::sundials_matrix::SUNMatrix;
use crate::sundials_types::UserData;
use crate::sunlinsol_band::BandLS;
use crate::sunlinsol_dense::DenseLS;
use crate::sunlinsol_pcg::PcgLS;
use crate::sunlinsol_spbcgs::SpbcgsLS;
use crate::sunlinsol_spfgmr::SpfgmrLS;
use crate::sunlinsol_spgmr::SpgmrLS;
use crate::sunlinsol_sptfqmr::SptfqmrLS;

/* SUNLinearSolver_Type */
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SUNLinearSolver_Type {
    SUNLINEARSOLVER_DIRECT,
    SUNLINEARSOLVER_ITERATIVE,
    SUNLINEARSOLVER_MATRIX_ITERATIVE,
    SUNLINEARSOLVER_MATRIX_EMBEDDED,
}
pub use SUNLinearSolver_Type::*;

/* Preconditioning and Gram-Schmidt types (sundials_iterative.h) */
pub const SUN_PREC_NONE: i32 = 0;
pub const SUN_PREC_LEFT: i32 = 1;
pub const SUN_PREC_RIGHT: i32 = 2;
pub const SUN_PREC_BOTH: i32 = 3;

pub const SUN_MODIFIED_GS: i32 = 1;
pub const SUN_CLASSICAL_GS: i32 = 2;

/* SUNLinearSolver return flags (sundials_linearsolver.h) */
pub const SUNLS_ATIMES_NULL: i32 = -804;
pub const SUNLS_ATIMES_FAIL_UNREC: i32 = -805;
pub const SUNLS_PSET_FAIL_UNREC: i32 = -806;
pub const SUNLS_PSOLVE_NULL: i32 = -807;
pub const SUNLS_PSOLVE_FAIL_UNREC: i32 = -808;
pub const SUNLS_GS_FAIL: i32 = -810;
pub const SUNLS_QRSOL_FAIL: i32 = -811;

pub const SUNLS_RECOV_FAILURE: i32 = 800;
pub const SUNLS_RES_REDUCED: i32 = 801;
pub const SUNLS_CONV_FAIL: i32 = 802;
pub const SUNLS_ATIMES_FAIL_REC: i32 = 803;
pub const SUNLS_PSET_FAIL_REC: i32 = 804;
pub const SUNLS_PSOLVE_FAIL_REC: i32 = 805;
pub const SUNLS_PACKAGE_FAIL_REC: i32 = 806;
pub const SUNLS_QRFACT_FAIL: i32 = 807;
pub const SUNLS_LUFACT_FAIL: i32 = 808;

/// ATimes callback: jv = J*v (SUNATimesFn).
pub type ATimesFn<'a> = dyn FnMut(&NVector, &mut NVector) -> i32 + 'a;
/// PSolve callback: solve P z = r with tolerance tol; lr=1 left, 2 right
/// (SUNPSolveFn).
pub type PSolveFn<'a> = dyn FnMut(&NVector, &mut NVector, f64, i32) -> i32 + 'a;

/// User-defined MATRIX_EMBEDDED solvers (SUNLinSolNewEmpty pattern,
/// e.g. cvAnalytic_mels). The integrator passes the current (t, gamma)
/// and user data so the solve can form and invert I - gamma*J itself.
pub trait CustomLinSol {
    fn ls_type(&self) -> SUNLinearSolver_Type {
        SUNLINEARSOLVER_MATRIX_EMBEDDED
    }
    fn solve(
        &mut self,
        x: &mut NVector,
        b: &NVector,
        tol: f64,
        t: f64,
        gamma: f64,
        user_data: &mut UserData,
    ) -> i32;
    fn last_flag(&self) -> i64 {
        0
    }
}

/// Generic SUNLinearSolver.
pub enum LinearSolver {
    Dense(DenseLS),
    Band(BandLS),
    Spgmr(SpgmrLS),
    Spfgmr(SpfgmrLS),
    Spbcgs(SpbcgsLS),
    Sptfqmr(SptfqmrLS),
    Pcg(PcgLS),
    Custom(Box<dyn CustomLinSol>),
}

impl LinearSolver {
    /// SUNLinSolGetType
    pub fn ls_type(&self) -> SUNLinearSolver_Type {
        match self {
            LinearSolver::Dense(_) | LinearSolver::Band(_) => SUNLINEARSOLVER_DIRECT,
            LinearSolver::Spgmr(_)
            | LinearSolver::Spfgmr(_)
            | LinearSolver::Spbcgs(_)
            | LinearSolver::Sptfqmr(_)
            | LinearSolver::Pcg(_) => SUNLINEARSOLVER_ITERATIVE,
            LinearSolver::Custom(c) => c.ls_type(),
        }
    }

    /// SUNLinSolInitialize
    pub fn initialize(&mut self) -> i32 {
        match self {
            LinearSolver::Dense(_) | LinearSolver::Band(_) | LinearSolver::Custom(_) => {
                SUN_SUCCESS
            }
            LinearSolver::Spgmr(s) => s.initialize(),
            LinearSolver::Spfgmr(s) => s.initialize(),
            LinearSolver::Spbcgs(s) => s.initialize(),
            LinearSolver::Sptfqmr(s) => s.initialize(),
            LinearSolver::Pcg(s) => s.initialize(),
        }
    }

    /// SUNLinSolSetup (direct solvers factor A; iterative ones no-op —
    /// their preconditioner setup happens through cvLsPSetup).
    pub fn setup(&mut self, a: Option<&mut SUNMatrix>) -> i32 {
        match self {
            LinearSolver::Dense(s) => match a {
                Some(SUNMatrix::Dense(am)) => s.setup(am),
                _ => SUN_ERR_ARG_INCOMPATIBLE,
            },
            LinearSolver::Band(s) => match a {
                Some(SUNMatrix::Band(am)) => s.setup(am),
                _ => SUN_ERR_ARG_INCOMPATIBLE,
            },
            _ => SUN_SUCCESS,
        }
    }

    /// SUNLinSolSolve for matrix-based and iterative solvers.
    /// Custom (matrix-embedded) solvers are dispatched separately in
    /// cvode_ls.rs because they need (t, gamma, user_data).
    pub fn solve(
        &mut self,
        a: Option<&mut SUNMatrix>,
        x: &mut NVector,
        b: &NVector,
        tol: f64,
        atimes: &mut ATimesFn,
        psolve: Option<&mut PSolveFn>,
        s1: Option<&NVector>,
        s2: Option<&NVector>,
    ) -> i32 {
        match self {
            LinearSolver::Dense(s) => match a {
                Some(SUNMatrix::Dense(am)) => s.solve(am, x, b),
                _ => SUN_ERR_ARG_INCOMPATIBLE,
            },
            LinearSolver::Band(s) => match a {
                Some(SUNMatrix::Band(am)) => s.solve(am, x, b),
                _ => SUN_ERR_ARG_INCOMPATIBLE,
            },
            LinearSolver::Spgmr(s) => s.solve(x, b, tol, atimes, psolve, s1, s2),
            LinearSolver::Spfgmr(s) => s.solve(x, b, tol, atimes, psolve, s1, s2),
            LinearSolver::Spbcgs(s) => s.solve(x, b, tol, atimes, psolve, s1, s2),
            LinearSolver::Sptfqmr(s) => s.solve(x, b, tol, atimes, psolve, s1, s2),
            LinearSolver::Pcg(s) => s.solve(x, b, tol, atimes, psolve, s1, s2),
            LinearSolver::Custom(_) => SUN_ERR_NOT_IMPLEMENTED,
        }
    }

    /// SUNLinSolSetPrecType
    pub fn set_prec_type(&mut self, pretype: i32) -> i32 {
        match self {
            LinearSolver::Spgmr(s) => s.set_prec_type(pretype),
            LinearSolver::Spfgmr(s) => s.set_prec_type(pretype),
            LinearSolver::Spbcgs(s) => s.set_prec_type(pretype),
            LinearSolver::Sptfqmr(s) => s.set_prec_type(pretype),
            LinearSolver::Pcg(s) => s.set_prec_type(pretype),
            _ => SUN_ERR_NOT_IMPLEMENTED,
        }
    }

    /// SUNLinSolSetGSType (SPGMR/SPFGMR only)
    pub fn set_gs_type(&mut self, gstype: i32) -> i32 {
        match self {
            LinearSolver::Spgmr(s) => s.set_gs_type(gstype),
            LinearSolver::Spfgmr(s) => s.set_gs_type(gstype),
            _ => SUN_ERR_NOT_IMPLEMENTED,
        }
    }

    /// SUNLinSol_*SetMaxl (SPBCGS/SPTFQMR/PCG)
    pub fn set_maxl(&mut self, maxl: i32) -> i32 {
        match self {
            LinearSolver::Spbcgs(s) => s.set_maxl(maxl),
            LinearSolver::Sptfqmr(s) => s.set_maxl(maxl),
            LinearSolver::Pcg(s) => s.set_maxl(maxl),
            _ => SUN_ERR_NOT_IMPLEMENTED,
        }
    }

    /// SUNLinSolSetZeroGuess
    pub fn set_zero_guess(&mut self, onoff: bool) {
        match self {
            LinearSolver::Spgmr(s) => s.zeroguess = onoff,
            LinearSolver::Spfgmr(s) => s.zeroguess = onoff,
            LinearSolver::Spbcgs(s) => s.zeroguess = onoff,
            LinearSolver::Sptfqmr(s) => s.zeroguess = onoff,
            LinearSolver::Pcg(s) => s.zeroguess = onoff,
            _ => {}
        }
    }

    /// SUNLinSolNumIters
    pub fn num_iters(&self) -> i32 {
        match self {
            LinearSolver::Spgmr(s) => s.numiters,
            LinearSolver::Spfgmr(s) => s.numiters,
            LinearSolver::Spbcgs(s) => s.numiters,
            LinearSolver::Sptfqmr(s) => s.numiters,
            LinearSolver::Pcg(s) => s.numiters,
            _ => 0,
        }
    }

    /// SUNLinSolResNorm
    pub fn res_norm(&self) -> f64 {
        match self {
            LinearSolver::Spgmr(s) => s.resnorm,
            LinearSolver::Spfgmr(s) => s.resnorm,
            LinearSolver::Spbcgs(s) => s.resnorm,
            LinearSolver::Sptfqmr(s) => s.resnorm,
            LinearSolver::Pcg(s) => s.resnorm,
            _ => 0.0,
        }
    }

    /// SUNLinSolSpace: (lenrwLS, leniwLS). Direct solvers report
    /// 2 + N integer words (pivots + counters), 0 real words, exactly
    /// like SUNLinSolSpace_Dense/_Band; iterative solvers report their
    /// workspace-vector formulas.
    pub fn space(&self) -> (i64, i64) {
        match self {
            LinearSolver::Dense(s) => (0, 2 + s.pivots.len() as i64),
            LinearSolver::Band(s) => (0, 2 + s.pivots.len() as i64),
            LinearSolver::Spgmr(s) => s.space(),
            LinearSolver::Spfgmr(s) => s.space(),
            LinearSolver::Spbcgs(s) => s.space(),
            LinearSolver::Sptfqmr(s) => s.space(),
            LinearSolver::Pcg(s) => s.space(),
            LinearSolver::Custom(_) => (0, 0),
        }
    }

    /// SUNLinSolResid: the residual N_Vector each iterative solver exposes
    /// (spgmr/spfgmr vtemp, spbcgs/pcg r, sptfqmr vtemp1). idaLsSolve/
    /// cvLsSolve copy this into b when the solve converges in 0 iterations
    /// (the preconditioned residual of the zero initial guess). Not defined
    /// for direct/matrix-embedded solvers.
    pub fn resid(&self) -> &NVector {
        match self {
            LinearSolver::Spgmr(s) => &s.vtemp,
            LinearSolver::Spfgmr(s) => &s.vtemp,
            LinearSolver::Spbcgs(s) => &s.r,
            LinearSolver::Sptfqmr(s) => &s.vtemp1,
            LinearSolver::Pcg(s) => &s.r,
            LinearSolver::Dense(_) | LinearSolver::Band(_) | LinearSolver::Custom(_) => {
                unreachable!("SUNLinSolResid is only defined for iterative solvers")
            }
        }
    }

    /// SUNLinSolLastFlag
    pub fn last_flag(&self) -> i64 {
        match self {
            LinearSolver::Dense(s) => s.last_flag,
            LinearSolver::Band(s) => s.last_flag,
            LinearSolver::Spgmr(s) => s.last_flag,
            LinearSolver::Spfgmr(s) => s.last_flag,
            LinearSolver::Spbcgs(s) => s.last_flag,
            LinearSolver::Sptfqmr(s) => s.last_flag,
            LinearSolver::Pcg(s) => s.last_flag,
            LinearSolver::Custom(c) => c.last_flag(),
        }
    }
}
