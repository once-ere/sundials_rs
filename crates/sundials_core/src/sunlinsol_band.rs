/* -----------------------------------------------------------------
 * Translated from src/sunlinsol/band/sunlinsol_band.c
 * (SUNDIALS/CVODE 7.7.0).
 *
 * The C SUNLinearSolverContent_Band (N, pivots, last_flag) becomes
 * `BandLS` (N is pivots.len()); the ops-table entries become the
 * methods below, dispatched from the `LinearSolver` enum in
 * sundials_linearsolver.rs.
 * -----------------------------------------------------------------*/
use crate::nvector_serial::{NVector, N_VScale};
use crate::sundials_band::{SUNDlsMat_bandGBTRF, SUNDlsMat_bandGBTRS};
use crate::sundials_context::SUNContext;
use crate::sundials_errors::{SUN_ERR_ARG_INCOMPATIBLE, SUN_SUCCESS};
use crate::sundials_linearsolver::{LinearSolver, SUNLS_LUFACT_FAIL};
use crate::sundials_matrix::SUNMatrix;
use crate::sunmatrix_band::BandMatrix;

const ONE: f64 = 1.0;

/// Content of the band linear solver (SUNLinearSolverContent_Band).
pub struct BandLS {
    pub pivots: Vec<i64>,
    pub last_flag: i64,
}

/// SUNLinSol_Band: create a band linear solver for the band matrix
/// `a` and template vector `y`. A must have storage upper bandwidth
/// s_mu >= min(n-1, mu+ml) so the factorization has room for fill-in.
/// Panics on the argument errors for which the C version returns NULL.
pub fn SUNLinSol_Band(y: &NVector, a: &SUNMatrix, _sunctx: &SUNContext) -> LinearSolver {
    let ab = match a {
        SUNMatrix::Band(ab) => ab,
        _ => panic!("SUNLinSol_Band: SUN_ERR_ARG_WRONGTYPE (matrix is not band)"),
    };

    /* Check that A has appropriate storage upper bandwidth for
    factorization */
    let matrix_rows = ab.n;
    assert!(
        ab.s_mu >= i64::min(matrix_rows - 1, ab.ml + ab.mu),
        "SUNLinSol_Band: SUN_ERR_ARG_INCOMPATIBLE (insufficient storage upper bandwidth)"
    );
    assert!(
        matrix_rows == y.len() as i64,
        "SUNLinSol_Band: SUN_ERR_ARG_DIMSMISMATCH (matrix and vector sizes differ)"
    );

    LinearSolver::Band(BandLS {
        pivots: vec![0; matrix_rows as usize],
        last_flag: 0,
    })
}

impl BandLS {
    /// SUNLinSolInitialize_Band: all solver-specific memory has
    /// already been allocated.
    pub fn initialize(&mut self) -> i32 {
        self.last_flag = SUN_SUCCESS as i64;
        SUN_SUCCESS
    }

    /// SUNLinSolSetup_Band: perform the band LU factorization of `a`
    /// in place. On a zero pivot, `last_flag` holds the (1-based)
    /// index returned by GBTRF and SUNLS_LUFACT_FAIL (recoverable,
    /// positive) is returned; on success last_flag = 0 / SUN_SUCCESS.
    pub fn setup(&mut self, a: &mut BandMatrix) -> i32 {
        /* ensure that storage upper bandwidth is sufficient for
        fill-in */
        if a.s_mu < i64::min(a.n - 1, a.mu + a.ml) {
            return SUN_ERR_ARG_INCOMPATIBLE;
        }

        /* perform LU factorization of input matrix */
        self.last_flag = SUNDlsMat_bandGBTRF(a, &mut self.pivots);

        /* store error flag (if nonzero, that row encountered a
        zero-valued pivot) */
        if self.last_flag > 0 {
            return SUNLS_LUFACT_FAIL;
        }
        SUN_SUCCESS
    }

    /// SUNLinSolSolve_Band: x = A^{-1} b using the LU factors in `a`.
    pub fn solve(&mut self, a: &mut BandMatrix, x: &mut NVector, b: &NVector) -> i32 {
        /* copy b into x */
        N_VScale(ONE, b, x);

        /* solve using LU factors */
        SUNDlsMat_bandGBTRS(a, &self.pivots, &mut x.data);
        self.last_flag = SUN_SUCCESS as i64;
        SUN_SUCCESS
    }

    /// SUNLinSolLastFlag_Band.
    pub fn last_flag(&self) -> i64 {
        self.last_flag
    }

    /// SUNLinSolSpace_Band: (lenrwLS, leniwLS).
    pub fn space(&self) -> (i64, i64) {
        (0, 2 + self.pivots.len() as i64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sundials_context::SUNContext_Create;
    use crate::sunmatrix_band::SUNBandMatrix;

    #[test]
    fn setup_and_solve_band_system() {
        let ctx = SUNContext_Create();
        /* tridiagonal: 4 on diagonal, -1 off-diagonals; n = 5 */
        let mut amat = SUNBandMatrix(5, 1, 1, &ctx);
        if let SUNMatrix::Band(ab) = &mut amat {
            for j in 0..5i64 {
                ab.set(j, j, 4.0);
                if j > 0 {
                    ab.set(j, j - 1, -1.0);
                    ab.set(j - 1, j, -1.0);
                }
            }
        }
        let y = NVector::new(5);
        let mut ls = match SUNLinSol_Band(&y, &amat, &ctx) {
            LinearSolver::Band(b) => b,
            _ => unreachable!(),
        };
        assert_eq!(ls.initialize(), SUN_SUCCESS);

        let ab = match &mut amat {
            SUNMatrix::Band(ab) => ab,
            _ => unreachable!(),
        };
        assert_eq!(ls.setup(ab), SUN_SUCCESS);
        assert_eq!(ls.last_flag(), 0);

        /* b = A * [1, 2, 3, 4, 5] */
        let b = NVector::from_slice(&[2.0, 4.0, 6.0, 8.0, 16.0]);
        let mut x = NVector::new(5);
        assert_eq!(ls.solve(ab, &mut x, &b), SUN_SUCCESS);
        for (xi, ei) in x.data.iter().zip([1.0, 2.0, 3.0, 4.0, 5.0].iter()) {
            assert!((xi - ei).abs() < 1e-12, "got {xi}, want {ei}");
        }
    }

    #[test]
    fn setup_reports_singular_matrix() {
        let ctx = SUNContext_Create();
        let mut amat = SUNBandMatrix(3, 1, 1, &ctx);
        /* left all-zero: zero pivot on the first elimination step */
        let y = NVector::new(3);
        let mut ls = match SUNLinSol_Band(&y, &amat, &ctx) {
            LinearSolver::Band(b) => b,
            _ => unreachable!(),
        };
        let ab = match &mut amat {
            SUNMatrix::Band(ab) => ab,
            _ => unreachable!(),
        };
        assert_eq!(ls.setup(ab), SUNLS_LUFACT_FAIL);
        assert_eq!(ls.last_flag(), 1); /* GBTRF returned k+1 = 1 */
    }
}
