/* -----------------------------------------------------------------
 * Translated from src/sunmatrix/band/sunmatrix_band.c
 * (SUNDIALS/CVODE 7.7.0).
 *
 * The C SUNMatrixContent_Band (M, N, mu, ml, s_mu, ldim, ldata,
 * data, cols) becomes `BandMatrix` (square: M == N == n). Storage:
 * column j occupies data[j*ldim .. j*ldim+ldim] and element (i,j)
 * lives at data[j*ldim + s_mu + i - j] (SM_ELEMENT_B). The
 * ops-table entries (SUNMatZero_Band, ...) are free functions,
 * dispatched from the `SUNMatrix` enum in sundials_matrix.rs.
 * -----------------------------------------------------------------*/
use crate::nvector_serial::NVector;
use crate::sundials_context::SUNContext;
use crate::sundials_errors::{SUNErrCode, SUN_ERR_ARG_DIMSMISMATCH, SUN_SUCCESS};
use crate::sundials_matrix::SUNMatrix;

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/// Band matrix (square, n-by-n), column-major banded storage.
#[derive(Clone, Debug)]
pub struct BandMatrix {
    pub n: i64,    /* dimension (SM_ROWS_B == SM_COLUMNS_B)      */
    pub mu: i64,   /* upper bandwidth (SM_UBAND_B)               */
    pub ml: i64,   /* lower bandwidth (SM_LBAND_B)               */
    pub s_mu: i64, /* stored upper bandwidth (SM_SUBAND_B)       */
    pub ldim: i64, /* leading dimension = s_mu + ml + 1          */
    pub data: Vec<f64>,
}

impl BandMatrix {
    /// Zero-filled band matrix with the given stored upper bandwidth.
    pub fn new(n: i64, mu: i64, ml: i64, smu: i64) -> Self {
        let ldim = smu + ml + 1;
        BandMatrix {
            n,
            mu,
            ml,
            s_mu: smu,
            ldim,
            data: vec![ZERO; (n * ldim) as usize],
        }
    }

    /// SM_ELEMENT_B(A, i, j) — valid for j-mu <= i <= j+ml
    /// (and up to j - s_mu during factorizations).
    #[inline]
    pub fn get(&self, i: i64, j: i64) -> f64 {
        self.data[(j * self.ldim + self.s_mu + i - j) as usize]
    }

    /// SM_ELEMENT_B(A, i, j) = v
    #[inline]
    pub fn set(&mut self, i: i64, j: i64, v: f64) {
        self.data[(j * self.ldim + self.s_mu + i - j) as usize] = v;
    }

    /// Whole stored column j (length ldim). The diagonal element of
    /// column j sits at offset s_mu; SM_COLUMN_B-style offset
    /// addressing (s_mu + i - j) is left to the caller.
    #[inline]
    pub fn col_mut(&mut self, j: i64) -> &mut [f64] {
        let start = (j * self.ldim) as usize;
        let ldim = self.ldim as usize;
        &mut self.data[start..start + ldim]
    }

    /// SM_LDATA_B(A)
    #[inline]
    pub fn ldata(&self) -> i64 {
        self.n * self.ldim
    }
}

/// SUNBandMatrix: create a band matrix with default stored upper
/// bandwidth s_mu = min(n-1, mu+ml) — enough for GBTRF fill-in.
pub fn SUNBandMatrix(n: i64, mu: i64, ml: i64, sunctx: &SUNContext) -> SUNMatrix {
    SUNBandMatrixStorage(n, mu, ml, i64::min(n - 1, mu + ml), sunctx)
}

/// SUNBandMatrixStorage: create a band matrix with specified storage
/// upper bandwidth. Panics on illegal dimension input (C returns NULL
/// under SUN_ERR_ARG_OUTOFRANGE).
pub fn SUNBandMatrixStorage(
    n: i64,
    mu: i64,
    ml: i64,
    smu: i64,
    _sunctx: &SUNContext,
) -> SUNMatrix {
    assert!(n > 0, "SUNBandMatrixStorage: SUN_ERR_ARG_OUTOFRANGE");
    assert!(smu >= 0, "SUNBandMatrixStorage: SUN_ERR_ARG_OUTOFRANGE");
    assert!(ml >= 0, "SUNBandMatrixStorage: SUN_ERR_ARG_OUTOFRANGE");
    SUNMatrix::Band(BandMatrix::new(n, mu, ml, smu))
}

/// SUNMatClone_Band: new band matrix of identical shape/bandwidths
/// (zeroed).
pub fn SUNMatClone_Band(a: &BandMatrix) -> BandMatrix {
    BandMatrix::new(a.n, a.mu, a.ml, a.s_mu)
}

/// SUNMatZero_Band.
pub fn SUNMatZero_Band(a: &mut BandMatrix) -> SUNErrCode {
    for v in a.data.iter_mut() {
        *v = ZERO;
    }
    SUN_SUCCESS
}

/// SUNMatCopy_Band: B = A. Grows B's stored bandwidth if A's
/// bandwidth is larger (as the C code reallocates B's data).
pub fn SUNMatCopy_Band(a: &BandMatrix, b: &mut BandMatrix) -> SUNErrCode {
    /* both matrices must have the same dimension */
    if a.n != b.n {
        return SUN_ERR_ARG_DIMSMISMATCH;
    }

    /* Grow B if A's bandwidth is larger */
    if a.mu > b.mu || a.ml > b.ml {
        let ml = i64::max(b.ml, a.ml);
        let mu = i64::max(b.mu, a.mu);
        let smu = i64::max(b.s_mu, a.s_mu);
        let col_size = smu + ml + 1;
        b.mu = mu;
        b.ml = ml;
        b.s_mu = smu;
        b.ldim = col_size;
        b.data = vec![ZERO; (b.n * col_size) as usize];
    }

    /* Perform operation */
    SUNMatZero_Band(b);
    for j in 0..b.n {
        let a_off = j * a.ldim + a.s_mu;
        let b_off = j * b.ldim + b.s_mu;
        for i in -a.mu..=a.ml {
            b.data[(b_off + i) as usize] = a.data[(a_off + i) as usize];
        }
    }
    SUN_SUCCESS
}

/// SUNMatScaleAddI_Band: A = c*A + I.
pub fn SUNMatScaleAddI_Band(c: f64, a: &mut BandMatrix) -> SUNErrCode {
    for j in 0..a.n {
        let off = j * a.ldim + a.s_mu;
        for i in -a.mu..=a.ml {
            a.data[(off + i) as usize] *= c;
        }
        a.data[off as usize] += ONE;
    }
    SUN_SUCCESS
}

/// SUNMatScaleAdd_Band: A = c*A + B. If B has larger bandwidth(s)
/// than A, A's storage grows to hold the result (SMScaleAddNew_Band).
pub fn SUNMatScaleAdd_Band(c: f64, a: &mut BandMatrix, b: &BandMatrix) -> SUNErrCode {
    if a.n != b.n {
        return SUN_ERR_ARG_DIMSMISMATCH;
    }

    /* Call separate routine if B has larger bandwidth(s) than A */
    if b.mu > a.mu || b.ml > a.ml {
        return SMScaleAddNew_Band(c, a, b);
    }

    /* Otherwise, perform operation in-place */
    for j in 0..a.n {
        let a_off = j * a.ldim + a.s_mu;
        let b_off = j * b.ldim + b.s_mu;
        for i in -b.mu..=b.ml {
            let ai = (a_off + i) as usize;
            a.data[ai] = c * a.data[ai] + b.data[(b_off + i) as usize];
        }
    }
    SUN_SUCCESS
}

/// SMScaleAddNew_Band (private in C): A = c*A + B computed in a new
/// matrix large enough to hold both, which then replaces A's content.
fn SMScaleAddNew_Band(c: f64, a: &mut BandMatrix, b: &BandMatrix) -> SUNErrCode {
    /* create new matrix large enough to hold both A and B */
    let ml = i64::max(a.ml, b.ml);
    let mu = i64::max(a.mu, b.mu);
    let smu = i64::min(a.n - 1, mu + ml);
    let mut cm = BandMatrix::new(a.n, mu, ml, smu);

    /* scale/add c*A into new matrix */
    for j in 0..a.n {
        let a_off = j * a.ldim + a.s_mu;
        let c_off = j * cm.ldim + cm.s_mu;
        for i in -a.mu..=a.ml {
            cm.data[(c_off + i) as usize] = c * a.data[(a_off + i) as usize];
        }
    }

    /* add B into new matrix */
    for j in 0..b.n {
        let b_off = j * b.ldim + b.s_mu;
        let c_off = j * cm.ldim + cm.s_mu;
        for i in -b.mu..=b.ml {
            cm.data[(c_off + i) as usize] += b.data[(b_off + i) as usize];
        }
    }

    /* replace A contents with C contents */
    *a = cm;
    SUN_SUCCESS
}

/// SUNMatMatvec_Band: y = A*x.
pub fn SUNMatMatvec_Band(a: &BandMatrix, x: &NVector, y: &mut NVector) -> SUNErrCode {
    if x.len() as i64 != a.n || y.len() as i64 != a.n {
        return SUN_ERR_ARG_DIMSMISMATCH;
    }

    let xd = &x.data;
    let yd = &mut y.data;
    let n = a.n;

    for yi in yd.iter_mut().take(n as usize) {
        *yi = ZERO;
    }
    for j in 0..n {
        let off = j * a.ldim + a.s_mu - j; /* col_j with SM_COLUMN_B offset */
        let is = i64::max(0, j - a.mu);
        let ie = i64::min(n - 1, j + a.ml);
        for i in is..=ie {
            yd[i as usize] += a.data[(off + i) as usize] * xd[j as usize];
        }
    }
    SUN_SUCCESS
}

/// SUNMatHermitianTransposeVec_Band: y = A^T x.
pub fn SUNMatHermitianTransposeVec_Band(
    a: &BandMatrix,
    x: &NVector,
    y: &mut NVector,
) -> SUNErrCode {
    if x.len() as i64 != a.n || y.len() as i64 != a.n {
        return SUN_ERR_ARG_DIMSMISMATCH;
    }

    let xd = &x.data;
    let yd = &mut y.data;
    let n = a.n;

    for yi in yd.iter_mut().take(n as usize) {
        *yi = ZERO;
    }
    for j in 0..n {
        let off = j * a.ldim + a.s_mu - j;
        let is = i64::max(0, j - a.mu);
        let ie = i64::min(n - 1, j + a.ml);
        for i in is..=ie {
            yd[j as usize] += a.data[(off + i) as usize] * xd[i as usize];
        }
    }
    SUN_SUCCESS
}

/// SUNMatSpace_Band: (lenrw, leniw).
pub fn SUNMatSpace_Band(a: &BandMatrix) -> (i64, i64) {
    (a.n * (a.s_mu + a.ml + 1), 7 + a.n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matvec_tridiagonal() {
        let mut a = BandMatrix::new(4, 1, 1, 2);
        for j in 0..4i64 {
            a.set(j, j, 2.0);
            if j > 0 {
                a.set(j, j - 1, 1.0);
                a.set(j - 1, j, -1.0);
            }
        }
        let x = NVector::from_slice(&[1.0, 2.0, 3.0, 4.0]);
        let mut y = NVector::new(4);
        assert_eq!(SUNMatMatvec_Band(&a, &x, &mut y), SUN_SUCCESS);
        /* rows: 2*1 - 2 = 0 ; 1 + 4 - 3 = 2 ; 2 + 6 - 4 = 4 ; 3 + 8 = 11 */
        for (yi, ei) in y.data.iter().zip([0.0, 2.0, 4.0, 11.0].iter()) {
            assert!((yi - ei).abs() < 1e-15, "got {yi}, want {ei}");
        }
    }

    #[test]
    fn copy_grows_bandwidth() {
        let mut a = BandMatrix::new(4, 2, 1, 3);
        a.set(0, 2, 7.0); /* needs mu = 2 */
        a.set(3, 2, -3.0);
        let mut b = BandMatrix::new(4, 0, 0, 0);
        assert_eq!(SUNMatCopy_Band(&a, &mut b), SUN_SUCCESS);
        assert_eq!(b.mu, 2);
        assert_eq!(b.ml, 1);
        assert!((b.get(0, 2) - 7.0).abs() < 1e-15);
        assert!((b.get(3, 2) + 3.0).abs() < 1e-15);
    }

    #[test]
    fn scale_add_grows_when_b_is_wider() {
        /* A diagonal, B tridiagonal: A = 2*A + B must grow A */
        let mut a = BandMatrix::new(3, 0, 0, 0);
        for j in 0..3i64 {
            a.set(j, j, 1.0);
        }
        let mut b = BandMatrix::new(3, 1, 1, 1);
        for j in 0..3i64 {
            b.set(j, j, 10.0);
            if j > 0 {
                b.set(j, j - 1, 5.0);
                b.set(j - 1, j, -5.0);
            }
        }
        assert_eq!(SUNMatScaleAdd_Band(2.0, &mut a, &b), SUN_SUCCESS);
        assert_eq!(a.mu, 1);
        assert_eq!(a.ml, 1);
        assert_eq!(a.s_mu, 2); /* min(n-1, mu+ml) = 2 */
        assert!((a.get(0, 0) - 12.0).abs() < 1e-15);
        assert!((a.get(1, 0) - 5.0).abs() < 1e-15);
        assert!((a.get(0, 1) + 5.0).abs() < 1e-15);

        /* in-place path: same bandwidths */
        let mut a2 = BandMatrix::new(3, 1, 1, 2);
        for j in 0..3i64 {
            a2.set(j, j, 1.0);
        }
        assert_eq!(SUNMatScaleAdd_Band(3.0, &mut a2, &b), SUN_SUCCESS);
        assert!((a2.get(0, 0) - 13.0).abs() < 1e-15);
        assert!((a2.get(2, 1) - 5.0).abs() < 1e-15);
    }
}
