/* -----------------------------------------------------------------
 * Translated from src/sunlinsol/dense/sunlinsol_dense.c
 * (SUNDIALS/CVODE 7.7.0).
 *
 * The C SUNLinearSolverContent_Dense (N, pivots, last_flag) becomes
 * `DenseLS` (N is pivots.len()); the ops-table entries become the
 * methods below, dispatched from the `LinearSolver` enum in
 * sundials_linearsolver.rs.
 * -----------------------------------------------------------------*/
use crate::nvector_serial::{NVector, N_VScale};
use crate::sundials_context::SUNContext;
use crate::sundials_dense::{SUNDlsMat_denseGETRF, SUNDlsMat_denseGETRS};
use crate::sundials_errors::SUN_SUCCESS;
use crate::sundials_linearsolver::{LinearSolver, SUNLS_LUFACT_FAIL};
use crate::sundials_matrix::SUNMatrix;
use crate::sunmatrix_dense::DenseMatrix;

const ONE: f64 = 1.0;

/// Content of the dense linear solver (SUNLinearSolverContent_Dense).
pub struct DenseLS {
    pub pivots: Vec<i64>,
    pub last_flag: i64,
}

/// SUNLinSol_Dense: create a dense linear solver for the square dense
/// matrix `a` and template vector `y`. Panics on the argument errors
/// for which the C version returns NULL.
pub fn SUNLinSol_Dense(y: &NVector, a: &SUNMatrix, _sunctx: &SUNContext) -> LinearSolver {
    let ad = match a {
        SUNMatrix::Dense(ad) => ad,
        _ => panic!("SUNLinSol_Dense: SUN_ERR_ARG_WRONGTYPE (matrix is not dense)"),
    };
    assert!(
        ad.m == ad.n,
        "SUNLinSol_Dense: SUN_ERR_ARG_DIMSMISMATCH (matrix not square)"
    );

    let matrix_rows = ad.m;
    assert!(
        matrix_rows == y.len() as i64,
        "SUNLinSol_Dense: SUN_ERR_ARG_DIMSMISMATCH (matrix and vector sizes differ)"
    );

    LinearSolver::Dense(DenseLS {
        pivots: vec![0; matrix_rows as usize],
        last_flag: 0,
    })
}

impl DenseLS {
    /// SUNLinSolInitialize_Dense: all solver-specific memory has
    /// already been allocated.
    pub fn initialize(&mut self) -> i32 {
        self.last_flag = SUN_SUCCESS as i64;
        SUN_SUCCESS
    }

    /// SUNLinSolSetup_Dense: perform the LU factorization of `a` in
    /// place. On a zero pivot, `last_flag` holds the (1-based) column
    /// index returned by GETRF and SUNLS_LUFACT_FAIL (recoverable,
    /// positive) is returned; on success last_flag = 0 / SUN_SUCCESS.
    pub fn setup(&mut self, a: &mut DenseMatrix) -> i32 {
        /* perform LU factorization of input matrix */
        self.last_flag = SUNDlsMat_denseGETRF(a, &mut self.pivots);

        /* store error flag (if nonzero, this row encountered a
        zero-valued pivot) */
        if self.last_flag > 0 {
            return SUNLS_LUFACT_FAIL;
        }
        SUN_SUCCESS
    }

    /// SUNLinSolSolve_Dense: x = A^{-1} b using the LU factors in `a`.
    pub fn solve(&mut self, a: &mut DenseMatrix, x: &mut NVector, b: &NVector) -> i32 {
        /* copy b into x */
        N_VScale(ONE, b, x);

        /* solve using LU factors */
        SUNDlsMat_denseGETRS(a, &self.pivots, &mut x.data);
        self.last_flag = SUN_SUCCESS as i64;
        SUN_SUCCESS
    }

    /// SUNLinSolLastFlag_Dense.
    pub fn last_flag(&self) -> i64 {
        self.last_flag
    }

    /// SUNLinSolSpace_Dense: (lenrwLS, leniwLS).
    pub fn space(&self) -> (i64, i64) {
        (0, 2 + self.pivots.len() as i64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sundials_context::SUNContext_Create;
    use crate::sunmatrix_dense::SUNDenseMatrix;

    #[test]
    fn setup_and_solve_dense_system() {
        let ctx = SUNContext_Create();
        let mut amat = SUNDenseMatrix(3, 3, &ctx);
        if let SUNMatrix::Dense(ad) = &mut amat {
            let rows = [
                [2.0, 1.0, 1.0],
                [4.0, -6.0, 0.0],
                [-2.0, 7.0, 2.0],
            ];
            for (i, r) in rows.iter().enumerate() {
                for (j, &v) in r.iter().enumerate() {
                    ad.set(i as i64, j as i64, v);
                }
            }
        }
        let y = NVector::new(3);
        let mut ls = match SUNLinSol_Dense(&y, &amat, &ctx) {
            LinearSolver::Dense(d) => d,
            _ => unreachable!(),
        };
        assert_eq!(ls.initialize(), SUN_SUCCESS);

        let ad = match &mut amat {
            SUNMatrix::Dense(ad) => ad,
            _ => unreachable!(),
        };
        assert_eq!(ls.setup(ad), SUN_SUCCESS);
        assert_eq!(ls.last_flag(), 0);

        /* b = A * [1, 2, 3] */
        let b = NVector::from_slice(&[7.0, -8.0, 18.0]);
        let mut x = NVector::new(3);
        assert_eq!(ls.solve(ad, &mut x, &b), SUN_SUCCESS);
        for (xi, ei) in x.data.iter().zip([1.0, 2.0, 3.0].iter()) {
            assert!((xi - ei).abs() < 1e-12, "got {xi}, want {ei}");
        }
    }

    #[test]
    fn setup_reports_singular_matrix() {
        let ctx = SUNContext_Create();
        let mut amat = SUNDenseMatrix(2, 2, &ctx);
        /* all-zero matrix: zero pivot in the first column */
        let y = NVector::new(2);
        let mut ls = match SUNLinSol_Dense(&y, &amat, &ctx) {
            LinearSolver::Dense(d) => d,
            _ => unreachable!(),
        };
        let ad = match &mut amat {
            SUNMatrix::Dense(ad) => ad,
            _ => unreachable!(),
        };
        assert_eq!(ls.setup(ad), SUNLS_LUFACT_FAIL);
        assert_eq!(ls.last_flag(), 1); /* GETRF returned k+1 = 1 */
    }
}
