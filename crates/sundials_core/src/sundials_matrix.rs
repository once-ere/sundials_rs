/* -----------------------------------------------------------------
 * Translated from src/sundials/sundials_matrix.c and
 * include/sundials/sundials_matrix.h (SUNDIALS/CVODE 7.7.0).
 *
 * The C generic SUNMatrix (ops-table over void* content) becomes an
 * enum over the concrete matrix implementations; the generic
 * SUNMat* wrappers become the methods / free functions below. The
 * NULL / ops-compatibility checks of the C code are covered by the
 * type system; shape checks that produced error returns are kept.
 * -----------------------------------------------------------------*/
use crate::nvector_serial::NVector;
use crate::sundials_errors::{SUNErrCode, SUN_ERR_ARG_WRONGTYPE};
use crate::sunmatrix_band::{self, BandMatrix};
use crate::sunmatrix_dense::{self, DenseMatrix};
use crate::sunmatrix_sparse::{self, SparseMatrix};

/// SUNMatrix_ID (sundials_matrix.h)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SUNMatrix_ID {
    SUNMATRIX_DENSE,
    SUNMATRIX_BAND,
    SUNMATRIX_SPARSE,
}
pub use SUNMatrix_ID::*;

/// Generic SUNMatrix.
#[derive(Clone, Debug)]
pub enum SUNMatrix {
    Dense(DenseMatrix),
    Band(BandMatrix),
    Sparse(SparseMatrix),
}

impl SUNMatrix {
    /// SUNMatGetID
    pub fn get_id(&self) -> SUNMatrix_ID {
        match self {
            SUNMatrix::Dense(_) => SUNMATRIX_DENSE,
            SUNMatrix::Band(_) => SUNMATRIX_BAND,
            SUNMatrix::Sparse(_) => SUNMATRIX_SPARSE,
        }
    }

    /// Number of rows.
    pub fn rows(&self) -> i64 {
        match self {
            SUNMatrix::Dense(a) => a.m,
            SUNMatrix::Band(a) => a.n,
            SUNMatrix::Sparse(a) => a.m,
        }
    }

    /// Number of columns.
    pub fn cols(&self) -> i64 {
        match self {
            SUNMatrix::Dense(a) => a.n,
            SUNMatrix::Band(a) => a.n,
            SUNMatrix::Sparse(a) => a.n,
        }
    }

    /// SUNMatZero: A_ij = 0.
    pub fn zero(&mut self) -> SUNErrCode {
        match self {
            SUNMatrix::Dense(a) => sunmatrix_dense::SUNMatZero_Dense(a),
            SUNMatrix::Band(a) => sunmatrix_band::SUNMatZero_Band(a),
            SUNMatrix::Sparse(a) => sunmatrix_sparse::SUNMatZero_Sparse(a),
        }
    }

    /// SUNMatCopy(A, B): B gets A's contents (B may grow: band
    /// bandwidth, sparse storage).
    pub fn copy_to(&self, b: &mut SUNMatrix) -> SUNErrCode {
        match (self, b) {
            (SUNMatrix::Dense(a), SUNMatrix::Dense(bm)) => {
                sunmatrix_dense::SUNMatCopy_Dense(a, bm)
            }
            (SUNMatrix::Band(a), SUNMatrix::Band(bm)) => {
                sunmatrix_band::SUNMatCopy_Band(a, bm)
            }
            (SUNMatrix::Sparse(a), SUNMatrix::Sparse(bm)) => {
                sunmatrix_sparse::SUNMatCopy_Sparse(a, bm)
            }
            _ => SUN_ERR_ARG_WRONGTYPE,
        }
    }

    /// SUNMatScaleAdd: A = c*A + B.
    pub fn scale_add(&mut self, c: f64, b: &SUNMatrix) -> SUNErrCode {
        match (self, b) {
            (SUNMatrix::Dense(a), SUNMatrix::Dense(bm)) => {
                sunmatrix_dense::SUNMatScaleAdd_Dense(c, a, bm)
            }
            (SUNMatrix::Band(a), SUNMatrix::Band(bm)) => {
                sunmatrix_band::SUNMatScaleAdd_Band(c, a, bm)
            }
            (SUNMatrix::Sparse(a), SUNMatrix::Sparse(bm)) => {
                sunmatrix_sparse::SUNMatScaleAdd_Sparse(c, a, bm)
            }
            _ => SUN_ERR_ARG_WRONGTYPE,
        }
    }

    /// SUNMatScaleAddI: A = c*A + I.
    pub fn scale_addi(&mut self, c: f64) -> SUNErrCode {
        match self {
            SUNMatrix::Dense(a) => sunmatrix_dense::SUNMatScaleAddI_Dense(c, a),
            SUNMatrix::Band(a) => sunmatrix_band::SUNMatScaleAddI_Band(c, a),
            SUNMatrix::Sparse(a) => sunmatrix_sparse::SUNMatScaleAddI_Sparse(c, a),
        }
    }

    /// SUNMatMatvec: y = A*x.
    pub fn matvec(&self, x: &NVector, y: &mut NVector) -> SUNErrCode {
        match self {
            SUNMatrix::Dense(a) => sunmatrix_dense::SUNMatMatvec_Dense(a, x, y),
            SUNMatrix::Band(a) => sunmatrix_band::SUNMatMatvec_Band(a, x, y),
            SUNMatrix::Sparse(a) => sunmatrix_sparse::SUNMatMatvec_Sparse(a, x, y),
        }
    }

    /// SUNMatHermitianTransposeVec: y = A^T x.
    pub fn hermitian_transpose_vec(&self, x: &NVector, y: &mut NVector) -> SUNErrCode {
        match self {
            SUNMatrix::Dense(a) => {
                sunmatrix_dense::SUNMatHermitianTransposeVec_Dense(a, x, y)
            }
            SUNMatrix::Band(a) => sunmatrix_band::SUNMatHermitianTransposeVec_Band(a, x, y),
            SUNMatrix::Sparse(a) => {
                sunmatrix_sparse::SUNMatHermitianTransposeVec_Sparse(a, x, y)
            }
        }
    }

    /// SUNMatClone: new matrix of the same shape/storage, zeroed.
    pub fn clone_empty(&self) -> SUNMatrix {
        match self {
            SUNMatrix::Dense(a) => SUNMatrix::Dense(sunmatrix_dense::SUNMatClone_Dense(a)),
            SUNMatrix::Band(a) => SUNMatrix::Band(sunmatrix_band::SUNMatClone_Band(a)),
            SUNMatrix::Sparse(a) => {
                SUNMatrix::Sparse(sunmatrix_sparse::SUNMatClone_Sparse(a))
            }
        }
    }

    /// SUNMatSpace: (lenrw, leniw).
    pub fn space(&self) -> (i64, i64) {
        match self {
            SUNMatrix::Dense(a) => sunmatrix_dense::SUNMatSpace_Dense(a),
            SUNMatrix::Band(a) => sunmatrix_band::SUNMatSpace_Band(a),
            SUNMatrix::Sparse(a) => sunmatrix_sparse::SUNMatSpace_Sparse(a),
        }
    }
}

/* -----------------------------------------------------------------
 * Free wrappers keeping the C names (sundials_matrix.c)
 * -----------------------------------------------------------------*/

/// SUNMatGetID
pub fn SUNMatGetID(a: &SUNMatrix) -> SUNMatrix_ID {
    a.get_id()
}

/// SUNMatClone
pub fn SUNMatClone(a: &SUNMatrix) -> SUNMatrix {
    a.clone_empty()
}

/// SUNMatZero
pub fn SUNMatZero(a: &mut SUNMatrix) -> SUNErrCode {
    a.zero()
}

/// SUNMatCopy: B gets A.
pub fn SUNMatCopy(a: &SUNMatrix, b: &mut SUNMatrix) -> SUNErrCode {
    a.copy_to(b)
}

/// SUNMatScaleAdd: A = c*A + B.
pub fn SUNMatScaleAdd(c: f64, a: &mut SUNMatrix, b: &SUNMatrix) -> SUNErrCode {
    a.scale_add(c, b)
}

/// SUNMatScaleAddI: A = c*A + I.
pub fn SUNMatScaleAddI(c: f64, a: &mut SUNMatrix) -> SUNErrCode {
    a.scale_addi(c)
}

/// SUNMatMatvec: y = A*x.
pub fn SUNMatMatvec(a: &SUNMatrix, x: &NVector, y: &mut NVector) -> SUNErrCode {
    a.matvec(x, y)
}

/// SUNMatHermitianTransposeVec: y = A^T x.
pub fn SUNMatHermitianTransposeVec(a: &SUNMatrix, x: &NVector, y: &mut NVector) -> SUNErrCode {
    a.hermitian_transpose_vec(x, y)
}

/// SUNMatSpace
pub fn SUNMatSpace(a: &SUNMatrix, lenrw: &mut i64, leniw: &mut i64) -> SUNErrCode {
    let (rw, iw) = a.space();
    *lenrw = rw;
    *leniw = iw;
    crate::sundials_errors::SUN_SUCCESS
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sundials_context::SUNContext_Create;
    use crate::sundials_errors::SUN_SUCCESS;
    use crate::sunmatrix_dense::SUNDenseMatrix;

    #[test]
    fn generic_ops_on_dense() {
        let ctx = SUNContext_Create();
        let mut a = SUNDenseMatrix(2, 2, &ctx);
        if let SUNMatrix::Dense(ad) = &mut a {
            ad.set(0, 0, 1.0);
            ad.set(0, 1, 2.0);
            ad.set(1, 0, 3.0);
            ad.set(1, 1, 4.0);
        }
        assert_eq!(a.get_id(), SUNMATRIX_DENSE);
        assert_eq!((a.rows(), a.cols()), (2, 2));

        let mut b = SUNMatClone(&a);
        assert_eq!(SUNMatCopy(&a, &mut b), SUN_SUCCESS);
        assert_eq!(SUNMatScaleAddI(2.0, &mut a), SUN_SUCCESS); /* A = 2A + I */
        assert_eq!(SUNMatScaleAdd(1.0, &mut a, &b), SUN_SUCCESS); /* A += B */

        /* A = 3*orig + I */
        if let SUNMatrix::Dense(ad) = &a {
            assert!((ad.get(0, 0) - 4.0).abs() < 1e-15);
            assert!((ad.get(0, 1) - 6.0).abs() < 1e-15);
            assert!((ad.get(1, 0) - 9.0).abs() < 1e-15);
            assert!((ad.get(1, 1) - 13.0).abs() < 1e-15);
        }

        assert_eq!(SUNMatZero(&mut a), SUN_SUCCESS);
        if let SUNMatrix::Dense(ad) = &a {
            assert!(ad.data.iter().all(|&v| v == 0.0));
        }

        /* cross-type dispatch is rejected */
        let mut band = crate::sunmatrix_band::SUNBandMatrix(2, 0, 0, &ctx);
        assert_eq!(SUNMatCopy(&b, &mut band), SUN_ERR_ARG_WRONGTYPE);
    }

    #[test]
    fn generic_matvec_band_and_space() {
        let ctx = SUNContext_Create();
        let mut a = crate::sunmatrix_band::SUNBandMatrix(3, 1, 1, &ctx);
        if let SUNMatrix::Band(ab) = &mut a {
            for j in 0..3i64 {
                ab.set(j, j, 2.0);
                if j > 0 {
                    ab.set(j, j - 1, -1.0);
                    ab.set(j - 1, j, -1.0);
                }
            }
        }
        let x = NVector::from_slice(&[1.0, 1.0, 1.0]);
        let mut y = NVector::new(3);
        assert_eq!(SUNMatMatvec(&a, &x, &mut y), SUN_SUCCESS);
        for (yi, ei) in y.data.iter().zip([1.0, 0.0, 1.0].iter()) {
            assert!((yi - ei).abs() < 1e-15, "got {yi}, want {ei}");
        }

        let (mut lenrw, mut leniw) = (0i64, 0i64);
        assert_eq!(SUNMatSpace(&a, &mut lenrw, &mut leniw), SUN_SUCCESS);
        /* smu = min(2, 2) = 2 -> lenrw = 3*(2+1+1) = 12, leniw = 10 */
        assert_eq!(lenrw, 12);
        assert_eq!(leniw, 10);
    }
}
