/* -----------------------------------------------------------------
 * Translated from src/sundials/sundials_band.c (SUNDIALS/CVODE
 * 7.7.0). Band storage follows sundials_direct.h / BandMatrix:
 * column j occupies data[j*ldim .. j*ldim+ldim] and element (i,j)
 * lives at data[j*ldim + s_mu + i - j]  (the C macro
 * ROW(i,j,smu) = i - j + smu inside column j).
 *
 * Kernels: band LU with partial pivoting (bandGBTRF/bandGBTRS) plus
 * the helpers bandCopy, bandScale, bandAddIdentity and bandMatvec.
 * The legacy `SUNDlsMat_BandGBTRF`-style wrappers over SUNDlsMat
 * collapse into the lowercase kernels here.
 * -----------------------------------------------------------------*/
use crate::sundials_types::sunindextype;
use crate::sundials_math::SUNRabs;
use crate::sunmatrix_band::BandMatrix;

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/// SUNDlsMat_bandGBTRF: in-place LU factorization with partial
/// pivoting of the n-by-n band matrix `a` (bandwidths mu/ml, stored
/// upper bandwidth s_mu >= min(n-1, mu+ml) to leave room for fill-in).
/// `p` (length >= n) receives the pivot rows.
///
/// Returns 0 on success; if a zero pivot is found at elimination step
/// k (0-based, k < n-1) it returns k+1, and if the last diagonal
/// element is zero it returns n — exactly as the C source.
pub fn SUNDlsMat_bandGBTRF(a: &mut BandMatrix, p: &mut [sunindextype]) -> sunindextype {
    let n = a.n as usize;
    let mu = a.mu as usize;
    let ml = a.ml as usize;
    let smu = a.s_mu as usize;
    let ldim = a.ldim as usize;
    let d = &mut a.data;

    /* zero out the first smu - mu rows of the rectangular array a */
    let num_rows = smu.saturating_sub(mu);
    if num_rows > 0 {
        for c in 0..n {
            let base = c * ldim;
            for r in 0..num_rows {
                d[base + r] = ZERO;
            }
        }
    }

    /* k = elimination step number */
    for k in 0..(n - 1) {
        let ck = k * ldim; /* start of column k */
        let diag_k = ck + smu; /* a(k,k) */
        let last_row_k = usize::min(n - 1, k + ml);

        /* find l = pivot row number */
        let mut l = k;
        let mut max = SUNRabs(d[diag_k]);
        for i in (k + 1)..=last_row_k {
            let v = SUNRabs(d[ck + smu + i - k]);
            if v > max {
                l = i;
                max = v;
            }
        }
        let storage_l = ck + smu + l - k; /* ROW(l,k,smu) in column k */
        p[k] = l as sunindextype;

        /* check for zero pivot element */
        if d[storage_l] == ZERO {
            return (k + 1) as sunindextype;
        }

        /* swap a(l,k) and a(k,k) if necessary */
        let swap = l != k;
        if swap {
            d.swap(storage_l, diag_k);
        }

        /* Scale the elements below the diagonal in         */
        /* column k by -1.0 / a(k,k). After the above swap, */
        /* a(k,k) holds the pivot element. This scaling     */
        /* stores the pivot row multipliers -a(i,k)/a(k,k)  */
        /* in a(i,k), i=k+1, ..., min(n-1,k+ml).            */
        let mult = -ONE / d[diag_k];
        for i in (k + 1)..=last_row_k {
            d[ck + smu + i - k] *= mult;
        }

        /* row_i = row_i - [a(i,k)/a(k,k)] row_k, i=k+1, ..., min(n-1,k+ml) */
        /* row k is the pivot row after swapping with row l.                */
        /* The computation is done one column at a time,                    */
        /* column j=k+1, ..., min(k+smu,n-1).                               */
        let last_col_k = usize::min(k + smu, n - 1);
        for j in (k + 1)..=last_col_k {
            let cj = j * ldim;
            let storage_l = cj + smu + l - j; /* ROW(l,j,smu) */
            let storage_k = cj + smu + k - j; /* ROW(k,j,smu) */
            let a_kj = d[storage_l];

            /* Swap the elements a(k,j) and a(k,l) if l!=k. */
            if swap {
                d[storage_l] = d[storage_k];
                d[storage_k] = a_kj;
            }

            /* a(i,j) = a(i,j) - [a(i,k)/a(k,k)]*a(k,j) */
            /* a_kj = a(k,j), *kptr = - a(i,k)/a(k,k), *jptr = a(i,j) */
            if a_kj != ZERO {
                for i in (k + 1)..=last_row_k {
                    d[cj + smu + i - j] += a_kj * d[ck + smu + i - k];
                }
            }
        }
    }

    /* set the last pivot row to be n-1 and check for a zero pivot */
    p[n - 1] = (n - 1) as sunindextype;
    if d[(n - 1) * ldim + smu] == ZERO {
        return n as sunindextype;
    }

    /* return 0 to indicate success */
    0
}

/// SUNDlsMat_bandGBTRS: solve A x = b using the band LU factors and
/// pivot array produced by [`SUNDlsMat_bandGBTRF`]. The solution
/// overwrites `b`.
pub fn SUNDlsMat_bandGBTRS(a: &BandMatrix, p: &[sunindextype], b: &mut [f64]) {
    let n = a.n as usize;
    let ml = a.ml as usize;
    let smu = a.s_mu as usize;
    let ldim = a.ldim as usize;
    let d = &a.data;

    /* Solve Ly = Pb, store solution y in b */
    for k in 0..(n - 1) {
        let l = p[k] as usize;
        let mult = b[l];
        if l != k {
            b[l] = b[k];
            b[k] = mult;
        }
        let diag_k = k * ldim + smu;
        let last_row_k = usize::min(n - 1, k + ml);
        for i in (k + 1)..=last_row_k {
            b[i] += mult * d[diag_k + i - k];
        }
    }

    /* Solve Ux = y, store solution x in b */
    for k in (0..n).rev() {
        let diag_k = k * ldim + smu;
        let first_row_k = k.saturating_sub(smu);
        b[k] /= d[diag_k];
        let mult = -b[k];
        for i in first_row_k..k {
            b[i] += mult * d[diag_k + i - k];
        }
    }
}

/// SUNDlsMat_bandCopy: B = A over the band (copymu, copyml).
pub fn SUNDlsMat_bandCopy(
    a: &BandMatrix,
    b: &mut BandMatrix,
    copymu: sunindextype,
    copyml: sunindextype,
) {
    let n = a.n;
    let copy_size = copymu + copyml + 1;

    for j in 0..n {
        let a_off = j * a.ldim + a.s_mu - copymu;
        let b_off = j * b.ldim + b.s_mu - copymu;
        for i in 0..copy_size {
            b.data[(b_off + i) as usize] = a.data[(a_off + i) as usize];
        }
    }
}

/// SUNDlsMat_bandScale: A = c*A (over the stored band mu..ml).
pub fn SUNDlsMat_bandScale(c: f64, a: &mut BandMatrix) {
    let n = a.n;
    let col_size = a.mu + a.ml + 1;

    for j in 0..n {
        let off = j * a.ldim + a.s_mu - a.mu;
        for i in 0..col_size {
            a.data[(off + i) as usize] *= c;
        }
    }
}

/// SUNDlsMat_bandAddIdentity: A += I.
pub fn SUNDlsMat_bandAddIdentity(a: &mut BandMatrix) {
    let n = a.n;
    for j in 0..n {
        a.data[(j * a.ldim + a.s_mu) as usize] += ONE;
    }
}

/// SUNDlsMat_bandMatvec: y = A*x (both of length n).
pub fn SUNDlsMat_bandMatvec(a: &BandMatrix, x: &[f64], y: &mut [f64]) {
    let n = a.n;

    for yi in y.iter_mut().take(n as usize) {
        *yi = ZERO;
    }

    for j in 0..n {
        let off = j * a.ldim + a.s_mu - a.mu;
        let is = if 0 > j - a.mu { 0 } else { j - a.mu };
        let ie = if n - 1 < j + a.ml { n - 1 } else { j + a.ml };
        for i in is..=ie {
            y[i as usize] += a.data[(off + i - j + a.mu) as usize] * x[j as usize];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tridiagonal test matrix: 4 on the diagonal, -1 on the sub- and
    /// super-diagonals. n=5, mu=ml=1, smu=min(n-1, mu+ml)=2.
    fn tridiag() -> BandMatrix {
        let mut a = BandMatrix::new(5, 1, 1, 2);
        for j in 0..5i64 {
            a.set(j, j, 4.0);
            if j > 0 {
                a.set(j, j - 1, -1.0);
                a.set(j - 1, j, -1.0);
            }
        }
        a
    }

    #[test]
    fn gbtrf_gbtrs_solves_tridiagonal() {
        let mut a = tridiag();
        /* b = A * [1, 2, 3, 4, 5] */
        let mut b = [2.0, 4.0, 6.0, 8.0, 16.0];
        let mut p = [0i64; 5];
        assert_eq!(SUNDlsMat_bandGBTRF(&mut a, &mut p), 0);
        SUNDlsMat_bandGBTRS(&a, &p, &mut b);
        for (bi, xi) in b.iter().zip([1.0, 2.0, 3.0, 4.0, 5.0].iter()) {
            assert!((bi - xi).abs() < 1e-12, "got {bi}, want {xi}");
        }
    }

    #[test]
    fn gbtrf_reports_zero_pivot() {
        /* all-zero matrix: zero pivot on the very first step -> 1 */
        let mut a = BandMatrix::new(3, 1, 1, 2);
        let mut p = [0i64; 3];
        assert_eq!(SUNDlsMat_bandGBTRF(&mut a, &mut p), 1);

        /* singular in the last diagonal position -> returns n */
        let mut a2 = BandMatrix::new(2, 0, 0, 0);
        a2.set(0, 0, 3.0);
        /* a2(1,1) left zero */
        let mut p2 = [0i64; 2];
        assert_eq!(SUNDlsMat_bandGBTRF(&mut a2, &mut p2), 2);
    }

    #[test]
    fn band_matvec_matches_dense_product() {
        let a = tridiag();
        let x = [1.0, -1.0, 2.0, 0.5, -2.0];
        let mut y = [0.0; 5];
        SUNDlsMat_bandMatvec(&a, &x, &mut y);
        /* rows: 4*1 -1*(-1) = 5 ; -1+(-4)+(-2)... compute directly */
        let expect = [
            4.0 * 1.0 - (-1.0),
            -1.0 + 4.0 * (-1.0) - 2.0,
            1.0 + 4.0 * 2.0 - 0.5,
            -2.0 + 4.0 * 0.5 + 2.0,
            -0.5 + 4.0 * (-2.0),
        ];
        for (yi, ei) in y.iter().zip(expect.iter()) {
            assert!((yi - ei).abs() < 1e-12, "got {yi}, want {ei}");
        }
    }
}
