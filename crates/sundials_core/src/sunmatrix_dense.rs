/* -----------------------------------------------------------------
 * Translated from src/sunmatrix/dense/sunmatrix_dense.c
 * (SUNDIALS/CVODE 7.7.0).
 *
 * The C SUNMatrixContent_Dense (M, N, ldata, data, cols) becomes
 * `DenseMatrix`; the `cols` pointer array is replaced by index
 * arithmetic into the column-major `data` Vec (element (i,j) at
 * data[j*m + i], matching SM_ELEMENT_D). The ops-table entries
 * (SUNMatZero_Dense, ...) are free functions here, dispatched from
 * the `SUNMatrix` enum in sundials_matrix.rs.
 * -----------------------------------------------------------------*/
use crate::nvector_serial::NVector;
use crate::sundials_context::SUNContext;
use crate::sundials_errors::{SUNErrCode, SUN_ERR_ARG_DIMSMISMATCH, SUN_SUCCESS};
use crate::sundials_matrix::SUNMatrix;

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/// Dense matrix, column-major storage: element (i,j) at data[j*m + i].
#[derive(Clone, Debug)]
pub struct DenseMatrix {
    pub m: i64, /* number of rows    (SM_ROWS_D)    */
    pub n: i64, /* number of columns (SM_COLUMNS_D) */
    pub data: Vec<f64>,
}

impl DenseMatrix {
    /// Zero-filled m-by-n matrix (calloc in the C constructor).
    pub fn new(m: i64, n: i64) -> Self {
        DenseMatrix {
            m,
            n,
            data: vec![ZERO; (m * n) as usize],
        }
    }

    /// SM_ELEMENT_D(A, i, j)
    #[inline]
    pub fn get(&self, i: i64, j: i64) -> f64 {
        self.data[(j * self.m + i) as usize]
    }

    /// SM_ELEMENT_D(A, i, j) = v
    #[inline]
    pub fn set(&mut self, i: i64, j: i64, v: f64) {
        self.data[(j * self.m + i) as usize] = v;
    }

    /// SM_COLUMN_D(A, j) (read-only)
    #[inline]
    pub fn col(&self, j: i64) -> &[f64] {
        let start = (j * self.m) as usize;
        &self.data[start..start + self.m as usize]
    }

    /// SM_COLUMN_D(A, j)
    #[inline]
    pub fn col_mut(&mut self, j: i64) -> &mut [f64] {
        let start = (j * self.m) as usize;
        let m = self.m as usize;
        &mut self.data[start..start + m]
    }

    /// SM_LDATA_D(A)
    #[inline]
    pub fn ldata(&self) -> i64 {
        self.m * self.n
    }
}

/// SUNDenseMatrix: create a new dense matrix (zero-filled), wrapped in
/// the generic `SUNMatrix` enum. Panics on illegal dimensions (the C
/// version returns NULL under SUN_ERR_ARG_OUTOFRANGE).
pub fn SUNDenseMatrix(m: i64, n: i64, _sunctx: &SUNContext) -> SUNMatrix {
    assert!(n > 0 && m > 0, "SUNDenseMatrix: SUN_ERR_ARG_OUTOFRANGE");
    SUNMatrix::Dense(DenseMatrix::new(m, n))
}

/// SUNMatClone_Dense: new dense matrix of the same shape (zeroed).
pub fn SUNMatClone_Dense(a: &DenseMatrix) -> DenseMatrix {
    DenseMatrix::new(a.m, a.n)
}

/// SUNMatZero_Dense: A_ij = 0.
pub fn SUNMatZero_Dense(a: &mut DenseMatrix) -> SUNErrCode {
    for v in a.data.iter_mut() {
        *v = ZERO;
    }
    SUN_SUCCESS
}

/// SUNMatCopy_Dense: B_ij = A_ij.
pub fn SUNMatCopy_Dense(a: &DenseMatrix, b: &mut DenseMatrix) -> SUNErrCode {
    /* both matrices must have the same shape */
    if a.m != b.m || a.n != b.n {
        return SUN_ERR_ARG_DIMSMISMATCH;
    }

    for j in 0..a.n {
        for i in 0..a.m {
            b.set(i, j, a.get(i, j));
        }
    }
    SUN_SUCCESS
}

/// SUNMatScaleAddI_Dense: A = c*A + I.
pub fn SUNMatScaleAddI_Dense(c: f64, a: &mut DenseMatrix) -> SUNErrCode {
    for j in 0..a.n {
        for i in 0..a.m {
            let idx = (j * a.m + i) as usize;
            a.data[idx] *= c;
            if i == j {
                a.data[idx] += ONE;
            }
        }
    }
    SUN_SUCCESS
}

/// SUNMatScaleAdd_Dense: A = c*A + B.
pub fn SUNMatScaleAdd_Dense(c: f64, a: &mut DenseMatrix, b: &DenseMatrix) -> SUNErrCode {
    if a.m != b.m || a.n != b.n {
        return SUN_ERR_ARG_DIMSMISMATCH;
    }

    for j in 0..a.n {
        for i in 0..a.m {
            let idx = (j * a.m + i) as usize;
            a.data[idx] = c * a.data[idx] + b.data[idx];
        }
    }
    SUN_SUCCESS
}

/// SUNMatMatvec_Dense: y = A*x.
pub fn SUNMatMatvec_Dense(a: &DenseMatrix, x: &NVector, y: &mut NVector) -> SUNErrCode {
    if x.len() as i64 != a.n || y.len() as i64 != a.m {
        return SUN_ERR_ARG_DIMSMISMATCH;
    }

    let xd = &x.data;
    let yd = &mut y.data;
    let m = a.m as usize;
    let n = a.n as usize;

    for yi in yd.iter_mut().take(m) {
        *yi = ZERO;
    }
    for j in 0..n {
        let cj = j * m;
        for i in 0..m {
            yd[i] += a.data[cj + i] * xd[j];
        }
    }
    SUN_SUCCESS
}

/// SUNMatHermitianTransposeVec_Dense: y = A^T x.
pub fn SUNMatHermitianTransposeVec_Dense(
    a: &DenseMatrix,
    x: &NVector,
    y: &mut NVector,
) -> SUNErrCode {
    if x.len() as i64 != a.m || y.len() as i64 != a.n {
        return SUN_ERR_ARG_DIMSMISMATCH;
    }

    let xd = &x.data;
    let yd = &mut y.data;
    let m = a.m as usize;
    let n = a.n as usize;

    for yi in yd.iter_mut().take(n) {
        *yi = ZERO;
    }
    for i in 0..n {
        let row_i = i * m; /* column i of A = row i of A^T */
        for j in 0..m {
            yd[i] += a.data[row_i + j] * xd[j];
        }
    }
    SUN_SUCCESS
}

/// SUNMatSpace_Dense: (lenrw, leniw).
pub fn SUNMatSpace_Dense(a: &DenseMatrix) -> (i64, i64) {
    (a.ldata(), 3 + a.n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matvec_and_transpose_vec() {
        let mut a = DenseMatrix::new(2, 3);
        /* A = [1 2 3; 4 5 6] */
        for (j, col) in [[1.0, 4.0], [2.0, 5.0], [3.0, 6.0]].iter().enumerate() {
            a.set(0, j as i64, col[0]);
            a.set(1, j as i64, col[1]);
        }
        let x = NVector::from_slice(&[1.0, -1.0, 2.0]);
        let mut y = NVector::new(2);
        assert_eq!(SUNMatMatvec_Dense(&a, &x, &mut y), SUN_SUCCESS);
        assert!((y.ith(0) - 5.0).abs() < 1e-15);
        assert!((y.ith(1) - 11.0).abs() < 1e-15);

        let xt = NVector::from_slice(&[1.0, 2.0]);
        let mut yt = NVector::new(3);
        assert_eq!(SUNMatHermitianTransposeVec_Dense(&a, &xt, &mut yt), SUN_SUCCESS);
        assert!((yt.ith(0) - 9.0).abs() < 1e-15);
        assert!((yt.ith(1) - 12.0).abs() < 1e-15);
        assert!((yt.ith(2) - 15.0).abs() < 1e-15);

        /* dimension mismatch is reported, not computed */
        let bad = NVector::new(5);
        let mut ybad = NVector::new(2);
        assert_eq!(
            SUNMatMatvec_Dense(&a, &bad, &mut ybad),
            SUN_ERR_ARG_DIMSMISMATCH
        );
    }

    #[test]
    fn scale_add_and_scale_addi() {
        let mut a = DenseMatrix::new(2, 2);
        a.set(0, 0, 1.0);
        a.set(0, 1, 2.0);
        a.set(1, 0, 3.0);
        a.set(1, 1, 4.0);
        let mut b = DenseMatrix::new(2, 2);
        assert_eq!(SUNMatCopy_Dense(&a, &mut b), SUN_SUCCESS);

        /* A = 2*A + B = 3*orig */
        assert_eq!(SUNMatScaleAdd_Dense(2.0, &mut a, &b), SUN_SUCCESS);
        assert!((a.get(1, 0) - 9.0).abs() < 1e-15);

        /* B = -1*B + I */
        assert_eq!(SUNMatScaleAddI_Dense(-1.0, &mut b), SUN_SUCCESS);
        assert!((b.get(0, 0) - 0.0).abs() < 1e-15);
        assert!((b.get(0, 1) + 2.0).abs() < 1e-15);
        assert!((b.get(1, 1) + 3.0).abs() < 1e-15);
    }
}
