/* -----------------------------------------------------------------
 * Translated from src/sundials/sundials_dense.c (SUNDIALS/CVODE
 * 7.7.0), with the column-major storage macros of
 * include/sundials/sundials_direct.h realized by `DenseMatrix`
 * (element (i,j) lives at data[j*m + i]).
 *
 * Kernels: LU with partial pivoting (denseGETRF/denseGETRS),
 * Cholesky (densePOTRF/densePOTRS), Householder QR
 * (denseGEQRF/denseORMQR) plus the small helpers denseCopy,
 * denseScale, denseAddIdentity and denseMatvec.
 *
 * The C `SUNDlsMat_DenseGETRF`-style wrappers (operating on the
 * legacy SUNDlsMat struct) collapse into the lowercase kernels here,
 * which take `DenseMatrix` / slices directly.
 * -----------------------------------------------------------------*/
use crate::sundials_math::{SUNRabs, SUNRsqrt};
use crate::sundials_types::sunindextype;
use crate::sunmatrix_dense::DenseMatrix;

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;

/// SUNDlsMat_denseGETRF: LU factorization with partial (row) pivoting of
/// the m-by-n matrix `a` (m >= n), done in place. `p` (length >= n)
/// receives the pivot row chosen at each elimination step.
///
/// Returns 0 on success. If a zero pivot is encountered at elimination
/// step k (0-based), returns k+1 exactly as the C source does
/// (`return (k + 1);`).
pub fn SUNDlsMat_denseGETRF(a: &mut DenseMatrix, p: &mut [sunindextype]) -> sunindextype {
    let m = a.m as usize;
    let n = a.n as usize;
    let d = &mut a.data;

    /* k-th elimination step number */
    for k in 0..n {
        let ck = k * m; /* start of column k */

        /* find l = pivot row number */
        let mut l = k;
        for i in (k + 1)..m {
            if SUNRabs(d[ck + i]) > SUNRabs(d[ck + l]) {
                l = i;
            }
        }
        p[k] = l as sunindextype;

        /* check for zero pivot element */
        if d[ck + l] == ZERO {
            return (k + 1) as sunindextype;
        }

        /* swap a(k,1:n) and a(l,1:n) if necessary */
        if l != k {
            for i in 0..n {
                d.swap(i * m + l, i * m + k);
            }
        }

        /* Scale the elements below the diagonal in column k by
         * 1.0/a(k,k). After the above swap a(k,k) holds the pivot
         * element. This scaling stores the pivot row multipliers
         * a(i,k)/a(k,k) in a(i,k), i=k+1, ..., m-1. */
        let mult = ONE / d[ck + k];
        for i in (k + 1)..m {
            d[ck + i] *= mult;
        }

        /* row_i = row_i - [a(i,k)/a(k,k)] row_k, i=k+1, ..., m-1 */
        /* row k is the pivot row after swapping with row l.       */
        /* The computation is done one column at a time,           */
        /* column j=k+1, ..., n-1.                                 */
        for j in (k + 1)..n {
            let cj = j * m;
            let a_kj = d[cj + k];

            /* a(i,j) = a(i,j) - [a(i,k)/a(k,k)]*a(k,j)  */
            /* a_kj = a(k,j), col_k[i] = - a(i,k)/a(k,k) */
            if a_kj != ZERO {
                for i in (k + 1)..m {
                    d[cj + i] -= a_kj * d[ck + i];
                }
            }
        }
    }

    /* return 0 to indicate success */
    0
}

/// SUNDlsMat_denseGETRS: solve A x = b using the LU factors and pivot
/// array produced by [`SUNDlsMat_denseGETRF`] (square system, n = a.n).
/// The solution overwrites `b`.
pub fn SUNDlsMat_denseGETRS(a: &DenseMatrix, p: &[sunindextype], b: &mut [f64]) {
    let m = a.m as usize; /* column stride */
    let n = a.n as usize;
    let d = &a.data;

    /* Permute b, based on pivot information in p */
    for k in 0..n {
        let pk = p[k] as usize;
        if pk != k {
            b.swap(k, pk);
        }
    }

    /* Solve Ly = b, store solution y in b */
    for k in 0..n.saturating_sub(1) {
        let ck = k * m;
        for i in (k + 1)..n {
            b[i] -= d[ck + i] * b[k];
        }
    }

    /* Solve Ux = y, store solution x in b */
    for k in (1..n).rev() {
        let ck = k * m;
        b[k] /= d[ck + k];
        for i in 0..k {
            b[i] -= d[ck + i] * b[k];
        }
    }
    b[0] /= d[0];
}

/// SUNDlsMat_densePOTRF: Cholesky decomposition of a symmetric
/// positive-definite matrix A = C*C^T (gaxpy version). Only the lower
/// triangle of A is accessed and it is overwritten with the lower
/// triangle of C. Returns 0 on success, j+1 if the (0-based) j-th
/// diagonal element is not positive.
pub fn SUNDlsMat_densePOTRF(a: &mut DenseMatrix) -> sunindextype {
    let m = a.m as usize;
    let stride = m; /* column stride equals row count */
    let d = &mut a.data;

    for j in 0..m {
        let cj = j * stride;

        if j > 0 {
            for i in j..m {
                for k in 0..j {
                    let ck = k * stride;
                    d[cj + i] -= d[ck + i] * d[ck + j];
                }
            }
        }

        let mut a_diag = d[cj + j];
        if a_diag <= ZERO {
            return (j + 1) as sunindextype;
        }
        a_diag = SUNRsqrt(a_diag);

        for i in j..m {
            d[cj + i] /= a_diag;
        }
    }

    0
}

/// SUNDlsMat_densePOTRS: solution of Ax=b, with A s.p.d., based on the
/// Cholesky decomposition obtained with densePOTRF (A = C*C^T, C lower
/// triangular). The solution overwrites `b`.
pub fn SUNDlsMat_densePOTRS(a: &DenseMatrix, b: &mut [f64]) {
    let m = a.m as usize;
    let stride = m;
    let d = &a.data;

    /* Solve C y = b, forward substitution - column version.
    Store solution y in b */
    for j in 0..m.saturating_sub(1) {
        let cj = j * stride;
        b[j] /= d[cj + j];
        for i in (j + 1)..m {
            b[i] -= b[j] * d[cj + i];
        }
    }
    b[m - 1] /= d[(m - 1) * stride + (m - 1)];

    /* Solve C^T x = y, backward substitution - row version.
    Store solution x in b */
    b[m - 1] /= d[(m - 1) * stride + (m - 1)];
    for i in (0..m - 1).rev() {
        let ci = i * stride;
        for j in (i + 1)..m {
            b[i] -= d[ci + j] * b[j];
        }
        b[i] /= d[ci + i];
    }
}

/// SUNDlsMat_denseGEQRF: QR factorization of a rectangular m-by-n
/// matrix A (m >= n) using Householder reflections. On exit the
/// elements on and above the diagonal contain R; the elements below
/// the diagonal, with the array `beta` (length n), represent the
/// orthogonal matrix Q as a product of elementary reflectors.
/// `v` (length m) must be provided as workspace. Returns 0.
pub fn SUNDlsMat_denseGEQRF(a: &mut DenseMatrix, beta: &mut [f64], v: &mut [f64]) -> i32 {
    let m = a.m as usize;
    let n = a.n as usize;
    let d = &mut a.data;

    /* For each column...*/
    for j in 0..n {
        let cj = j * m;
        let ajj = d[cj + j];

        /* Compute the j-th Householder vector (of length m-j) */
        v[0] = ONE;
        let mut s = ZERO;
        for i in 1..(m - j) {
            v[i] = d[cj + i + j];
            s += v[i] * v[i];
        }

        if s != ZERO {
            let mu = SUNRsqrt(ajj * ajj + s);
            let v1 = if ajj <= ZERO { ajj - mu } else { -s / (ajj + mu) };
            let v1_2 = v1 * v1;
            beta[j] = TWO * v1_2 / (s + v1_2);
            for i in 1..(m - j) {
                v[i] /= v1;
            }
        } else {
            beta[j] = ZERO;
        }

        /* Update upper triangle of A (load R) */
        for k in j..n {
            let ck = k * m;
            let mut s = ZERO;
            for i in 0..(m - j) {
                s += d[ck + i + j] * v[i];
            }
            s *= beta[j];
            for i in 0..(m - j) {
                d[ck + i + j] -= s * v[i];
            }
        }

        /* Update A (load Householder vector) */
        if j < m - 1 {
            for i in 1..(m - j) {
                d[cj + i + j] = v[i];
            }
        }
    }

    0
}

/// SUNDlsMat_denseORMQR: computes vm = Q * vn, where the orthogonal
/// matrix Q is stored as elementary reflectors in the m-by-n matrix A
/// and in the vector `beta` (as produced by denseGEQRF). `vn` has
/// length n, `vm` has length m (m >= n), and `v` (length m) is
/// workspace. Returns 0.
pub fn SUNDlsMat_denseORMQR(
    a: &DenseMatrix,
    beta: &[f64],
    vn: &[f64],
    vm: &mut [f64],
    v: &mut [f64],
) -> i32 {
    let m = a.m as usize;
    let n = a.n as usize;
    let d = &a.data;

    /* Initialize vm */
    for i in 0..n {
        vm[i] = vn[i];
    }
    for i in n..m {
        vm[i] = ZERO;
    }

    /* Accumulate (backwards) corrections into vm */
    for j in (0..n).rev() {
        let cj = j * m;

        v[0] = ONE;
        let mut s = vm[j];
        for i in 1..(m - j) {
            v[i] = d[cj + i + j];
            s += v[i] * vm[i + j];
        }
        s *= beta[j];

        for i in 0..(m - j) {
            vm[i + j] -= s * v[i];
        }
    }

    0
}

/// SUNDlsMat_denseCopy: B = A (copies a.m x a.n block).
pub fn SUNDlsMat_denseCopy(a: &DenseMatrix, b: &mut DenseMatrix) {
    let m = a.m as usize;
    let n = a.n as usize;
    let bm = b.m as usize;
    for j in 0..n {
        for i in 0..m {
            b.data[j * bm + i] = a.data[j * m + i];
        }
    }
}

/// SUNDlsMat_denseScale: A = c*A.
pub fn SUNDlsMat_denseScale(c: f64, a: &mut DenseMatrix) {
    let m = a.m as usize;
    let n = a.n as usize;
    for j in 0..n {
        for i in 0..m {
            a.data[j * m + i] *= c;
        }
    }
}

/// SUNDlsMat_denseAddIdentity: A += I (over the a.n leading diagonal
/// entries, as the C version is called with n = N on square matrices).
pub fn SUNDlsMat_denseAddIdentity(a: &mut DenseMatrix) {
    let m = a.m as usize;
    let n = a.n as usize;
    for i in 0..n {
        a.data[i * m + i] += ONE;
    }
}

/// SUNDlsMat_denseMatvec: y = A*x (x length n, y length m).
pub fn SUNDlsMat_denseMatvec(a: &DenseMatrix, x: &[f64], y: &mut [f64]) {
    let m = a.m as usize;
    let n = a.n as usize;

    for yi in y.iter_mut().take(m) {
        *yi = ZERO;
    }

    for j in 0..n {
        let cj = j * m;
        for i in 0..m {
            y[i] += a.data[cj + i] * x[j];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dense_from_rows(rows: &[&[f64]]) -> DenseMatrix {
        let m = rows.len() as i64;
        let n = rows[0].len() as i64;
        let mut a = DenseMatrix::new(m, n);
        for (i, r) in rows.iter().enumerate() {
            for (j, &v) in r.iter().enumerate() {
                a.set(i as i64, j as i64, v);
            }
        }
        a
    }

    #[test]
    fn getrf_getrs_solves_3x3() {
        let mut a = dense_from_rows(&[
            &[2.0, 1.0, 1.0],
            &[4.0, -6.0, 0.0],
            &[-2.0, 7.0, 2.0],
        ]);
        /* b = A * [1, 2, 3] */
        let mut b = [7.0, -8.0, 18.0];
        let mut p = [0i64; 3];
        let ret = SUNDlsMat_denseGETRF(&mut a, &mut p);
        assert_eq!(ret, 0);
        SUNDlsMat_denseGETRS(&a, &p, &mut b);
        for (bi, xi) in b.iter().zip([1.0, 2.0, 3.0].iter()) {
            assert!((bi - xi).abs() < 1e-12, "got {bi}, want {xi}");
        }
    }

    #[test]
    fn getrf_zero_pivot_returns_k_plus_one() {
        /* Singular: second elimination step hits a zero pivot -> returns 2 */
        let mut a = dense_from_rows(&[&[1.0, 2.0], &[2.0, 4.0]]);
        let mut p = [0i64; 2];
        assert_eq!(SUNDlsMat_denseGETRF(&mut a, &mut p), 2);

        /* Zero first column -> zero pivot at step 0 -> returns 1 */
        let mut a0 = dense_from_rows(&[&[0.0, 1.0], &[0.0, 1.0]]);
        let mut p0 = [0i64; 2];
        assert_eq!(SUNDlsMat_denseGETRF(&mut a0, &mut p0), 1);
    }

    #[test]
    fn potrf_potrs_solves_spd_3x3() {
        let mut a = dense_from_rows(&[
            &[4.0, 2.0, 0.0],
            &[2.0, 3.0, 1.0],
            &[0.0, 1.0, 2.0],
        ]);
        /* b = A * [1, -1, 2] */
        let mut b = [2.0, 1.0, 3.0];
        assert_eq!(SUNDlsMat_densePOTRF(&mut a), 0);
        SUNDlsMat_densePOTRS(&a, &mut b);
        for (bi, xi) in b.iter().zip([1.0, -1.0, 2.0].iter()) {
            assert!((bi - xi).abs() < 1e-12, "got {bi}, want {xi}");
        }
        /* non-positive-definite matrix reports 1-based failure column */
        let mut bad = dense_from_rows(&[&[-1.0, 0.0], &[0.0, 1.0]]);
        assert_eq!(SUNDlsMat_densePOTRF(&mut bad), 1);
    }

    #[test]
    fn geqrf_ormqr_reproduce_a_times_x() {
        /* A is 3x2; y = A*x must equal Q*(R*x) */
        let a0 = dense_from_rows(&[&[1.0, 2.0], &[3.0, 4.0], &[5.0, 6.0]]);
        let x = [1.0, 2.0];
        let mut y = [0.0; 3];
        SUNDlsMat_denseMatvec(&a0, &x, &mut y);

        let mut a = DenseMatrix::new(3, 2);
        SUNDlsMat_denseCopy(&a0, &mut a);
        let mut beta = [0.0; 2];
        let mut v = [0.0; 3];
        assert_eq!(SUNDlsMat_denseGEQRF(&mut a, &mut beta, &mut v), 0);

        /* R*x from the factored upper triangle */
        let rx = [
            a.get(0, 0) * x[0] + a.get(0, 1) * x[1],
            a.get(1, 1) * x[1],
        ];
        let mut vm = [0.0; 3];
        assert_eq!(SUNDlsMat_denseORMQR(&a, &beta, &rx, &mut vm, &mut v), 0);
        for (vi, yi) in vm.iter().zip(y.iter()) {
            assert!((vi - yi).abs() < 1e-12, "got {vi}, want {yi}");
        }
    }
}
