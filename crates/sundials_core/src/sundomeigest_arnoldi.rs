/* -----------------------------------------------------------------
 * Translation of
 * sundials-7.7.0/src/sundomeigest/arnoldi/sundomeigest_arnoldi.c
 * (+ include/sundomeigest/sundomeigest_arnoldi.h).
 *
 * Arnoldi Iteration implementation of the SUNDomEigEst package. The
 * C content struct is the Arnoldi variant payload of the
 * SUNDomEigEstimator enum; implementation ops called on an estimator
 * of another variant return SUN_ERR_ARG_INCOMPATIBLE.
 *
 * LAPACK adaptation: the C module has no non-LAPACK path — it packs
 * the (kry_dim x kry_dim) upper-Hessenberg matrix produced by the
 * Arnoldi/SUNModifiedGS loop into LAPACK_A (column-major) and calls
 * dgeev (jobvl = jobvr = 'N', eigenvalues only). Since LAPACK/FFI is
 * excluded from this workspace, the eigenvalue computation is done by
 * a faithful pure-Rust Francis double-shift QR iteration for real
 * upper-Hessenberg matrices (the classical EISPACK HQR algorithm —
 * the same algorithm dgeev's dhseqr applies after reducing a general
 * matrix to Hessenberg form; here the input is already Hessenberg, so
 * no reduction step is needed). As in HQR/dhseqr, the iteration for
 * any single eigenvalue fails after 30 shifts without deflation; that
 * failure maps to LAPACK's info > 0 and is reported as
 * SUN_ERR_EXT_FAIL exactly like the C code. Complex conjugate pairs
 * are stored with the positive imaginary part first, matching the
 * dgeev output convention. Eigenvalues may differ from dgeev results
 * in the last floating-point digits (dgeev balances the matrix
 * first); the estimator's users only consume the dominant eigenvalue
 * to modest accuracy.
 *
 * Consequently the C content fields LAPACK_work / LAPACK_lwork (the
 * dgeev workspace and its size query) have no Rust counterpart. The
 * remaining LAPACK_* fields (packed matrix, eigenvalue real and
 * imaginary parts, sort array) are kept with the C names and used by
 * the same steps as the C code.
 *
 * Note (fidelity): SUNModifiedGS only fills the upper-Hessenberg part
 * of Hes; the entries below the first subdiagonal are never written.
 * The C code mallocs Hes and packs those uninitialized entries into
 * LAPACK_A (in practice zero pages). Here Hes is zero-initialized, so
 * the packed matrix is deterministically upper Hessenberg, and the
 * QR iteration never reads below the first subdiagonal anyway.
 *
 * Other adaptations as in sundomeigest_power.rs: ATimes is supplied
 * at Estimate time as a closure argument (pinned adaptation; the C
 * stores it via SetATimes); it is no longer a stored
 * boxed closure (C ATdata captured); constructor SUNAssertNull
 * argument checks are debug-only in C and drop out; the SUNAssert
 * guards for a missing ATimes / not-yet-Initialized workspace are
 * kept as early returns; Destroy_Arnoldi is ownership drop (see the
 * base module).
 * -----------------------------------------------------------------*/

use crate::nvector_serial::{NVector, N_VClone, N_VDotProd, N_VScale};
use crate::sundials_context::SUNContext;
use crate::sundials_domeigestimator::SUNDomEigEstimator;
use crate::sundials_errors::{
    SUNErrCode, SUN_ERR_ARG_CORRUPT, SUN_ERR_ARG_INCOMPATIBLE, SUN_ERR_EXT_FAIL,
    SUN_ERR_USER_FCN_FAIL, SUN_SUCCESS,
};
use crate::sundials_iterative::SUNModifiedGS;
use crate::sundials_linearsolver::ATimesFn;
use crate::sundials_math::{SUNRabs, SUNRsqrt};

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/* Default estimator parameters */
const DEE_NUM_OF_WARMUPS_ARNOLDI_DEFAULT: i32 = 100;

/* Default Arnoldi Iteration parameters */
const DEE_KRYLOV_DIM_DEFAULT: i32 = 3;

/// C struct SUNDomEigEstimatorContent_Arnoldi_
pub struct SUNDomEigEstimatorContent_Arnoldi {
    /// Krylov subspace vectors
    pub V: Vec<NVector>,
    pub q: NVector,
    /// Krylov subspace dimension
    pub kry_dim: i32,
    /// Number of preprocessing iterations
    pub num_warmups: i32,
    /// Number of iterations in last Estimate call
    pub num_iters: i64,
    /// Number of ATimes calls
    pub num_ATimes: i64,
    /// The Hessenberg matrix packed in column-major order (dgeev input)
    pub LAPACK_A: Vec<f64>,
    /// Real parts of eigenvalues
    pub LAPACK_wr: Vec<f64>,
    /// Imaginary parts of eigenvalues
    pub LAPACK_wi: Vec<f64>,
    /// an array to sort eigenvalues
    pub LAPACK_arr: Vec<[f64; 2]>,
    /// Hessenberg matrix Hes
    pub Hes: Vec<Vec<f64>>,
}

fn content(
    dee: &SUNDomEigEstimator,
) -> Result<&SUNDomEigEstimatorContent_Arnoldi, SUNErrCode> {
    match dee {
        SUNDomEigEstimator::Arnoldi(c) => Ok(c),
        _ => Err(SUN_ERR_ARG_INCOMPATIBLE),
    }
}

fn content_mut(
    dee: &mut SUNDomEigEstimator,
) -> Result<&mut SUNDomEigEstimatorContent_Arnoldi, SUNErrCode> {
    match dee {
        SUNDomEigEstimator::Arnoldi(c) => Ok(c),
        _ => Err(SUN_ERR_ARG_INCOMPATIBLE),
    }
}

/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/// Function to create a new Arnoldi estimator
pub fn SUNDomEigEstimator_Arnoldi(
    q: &NVector,
    kry_dim: i32,
    _sunctx: &SUNContext,
) -> SUNDomEigEstimator {
    /* Check if kry_dim >= 2 */
    let kry_dim = if kry_dim < 3 {
        DEE_KRYLOV_DIM_DEFAULT
    } else {
        kry_dim
    };

    /* Allocate content: content->q = N_VClone(q); N_VScale(ONE, q, content->q);
    content->V = N_VCloneVectorArray(kry_dim + 1, q) */
    let mut cq = N_VClone(q);
    N_VScale(ONE, q, &mut cq);
    let cv: Vec<NVector> = (0..(kry_dim + 1)).map(|_| N_VClone(q)).collect();

    SUNDomEigEstimator::Arnoldi(SUNDomEigEstimatorContent_Arnoldi {
        V: cv,
        q: cq,
        kry_dim,
        num_warmups: DEE_NUM_OF_WARMUPS_ARNOLDI_DEFAULT,
        num_iters: 0,
        num_ATimes: 0,
        LAPACK_A: Vec::new(),
        LAPACK_wr: Vec::new(),
        LAPACK_wi: Vec::new(),
        LAPACK_arr: Vec::new(),
        Hes: Vec::new(),
    })
}

/* -----------------------------------------------------------------
 * implementation of dominant eigenvalue estimator operations
 * ----------------------------------------------------------------- */

pub fn SUNDomEigEstimator_Initialize_Arnoldi(dee: &mut SUNDomEigEstimator) -> SUNErrCode {
    let c = match content_mut(dee) {
        Ok(c) => c,
        Err(e) => return e,
    };

    if c.kry_dim < 2 {
        c.kry_dim = DEE_KRYLOV_DIM_DEFAULT;
    }
    if c.num_warmups < 0 {
        c.num_warmups = DEE_NUM_OF_WARMUPS_ARNOLDI_DEFAULT;
    }

    /* (ATimes presence check dropped: it arrives at Estimate time) */

    let kd = c.kry_dim as usize;
    if c.LAPACK_A.is_empty() {
        c.LAPACK_A = vec![ZERO; kd * kd];
    }
    if c.LAPACK_wr.is_empty() {
        c.LAPACK_wr = vec![ZERO; kd];
    }
    if c.LAPACK_wi.is_empty() {
        c.LAPACK_wi = vec![ZERO; kd];
    }

    /* (C queries the dgeev workspace size here and allocates
    LAPACK_work; the pure-Rust QR iteration needs no workspace.) */

    /* LAPACK array */
    c.LAPACK_arr = vec![[ZERO; 2]; kd];

    /* Hessenberg matrix Hes */
    c.Hes = vec![vec![ZERO; kd]; kd + 1];

    let mut normq = N_VDotProd(&c.q, &c.q);
    normq = SUNRsqrt(normq);

    N_VScale(ONE / normq, &c.q, &mut c.V[0]);

    SUN_SUCCESS
}

pub fn SUNDomEigEstimator_SetNumPreprocessIters_Arnoldi(
    dee: &mut SUNDomEigEstimator,
    num_iters: i32,
) -> SUNErrCode {
    let c = match content_mut(dee) {
        Ok(c) => c,
        Err(e) => return e,
    };
    /* Check if num_iters >= 0 */
    let num_iters = if num_iters < 0 {
        DEE_NUM_OF_WARMUPS_ARNOLDI_DEFAULT
    } else {
        num_iters
    };
    /* set the number of warmups */
    c.num_warmups = num_iters;
    SUN_SUCCESS
}

pub fn SUNDomEigEstimator_SetInitialGuess_Arnoldi(
    dee: &mut SUNDomEigEstimator,
    q: &NVector,
) -> SUNErrCode {
    let c = match content_mut(dee) {
        Ok(c) => c,
        Err(e) => return e,
    };

    let mut normq = N_VDotProd(q, q);
    normq = SUNRsqrt(normq);

    /* set the initial guess */
    N_VScale(ONE / normq, q, &mut c.V[0]);

    SUN_SUCCESS
}

pub fn SUNDomEigEstimator_Estimate_Arnoldi(
    dee: &mut SUNDomEigEstimator,
    atimes: &mut ATimesFn,
    lambdaR: &mut f64,
    lambdaI: &mut f64,
) -> SUNErrCode {
    let c = match content_mut(dee) {
        Ok(c) => c,
        Err(e) => return e,
    };

    /* C: SUNAssert(Hes) — the workspace exists only after Initialize */
    if c.Hes.is_empty() {
        return SUN_ERR_ARG_CORRUPT;
    }

    let n = c.kry_dim as usize;
    let mut normq: f64;
    c.num_ATimes = 0;
    c.num_iters = 0;

    /* Set the initial q = A^{num_warmups}q/||A^{num_warmups}q|| */
    for _i in 0..c.num_warmups {
        let retval = atimes(&c.V[0], &mut c.q);
        c.num_ATimes += 1;
        c.num_iters += 1;
        if retval != 0 {
            return SUN_ERR_USER_FCN_FAIL;
        }

        normq = N_VDotProd(&c.q, &c.q);
        normq = SUNRsqrt(normq);
        N_VScale(ONE / normq, &c.q, &mut c.V[0]);
    }

    for i in 0..n {
        /* Compute the next Krylov vector */
        {
            let (vh, vt) = c.V.split_at_mut(i + 1);
            let retval = atimes(&vh[i], &mut vt[0]);
            c.num_ATimes += 1;
            c.num_iters += 1;
            if retval != 0 {
                return SUN_ERR_USER_FCN_FAIL;
            }
        }

        /* C: SUNModifiedGS(V, Hes, i+1, n, &(Hes[i+1][i])) — the norm
        destination aliases an entry of Hes, which MGS neither reads
        nor writes for k = i+1; assign it after the call. */
        let mut new_vk_norm = ZERO;
        let retval = SUNModifiedGS(&mut c.V, &mut c.Hes, (i + 1) as i32, n as i32, &mut new_vk_norm);
        if retval != SUN_SUCCESS {
            return retval;
        }
        c.Hes[i + 1][i] = new_vk_norm;

        /* Unitize the computed orthogonal vector (aliased N_VScale) */
        c.V[i + 1].scale_inplace(ONE / c.Hes[i + 1][i]);
    }

    /* Pack the Hessenberg matrix in column-major order for the (LAPACK
    dgeev in C) eigenvalue computation */
    let mut k = 0;
    for j in 0..n {
        for i in 0..n {
            c.LAPACK_A[k] = c.Hes[i][j];
            k += 1;
        }
    }

    /* Compute all eigenvalues of the packed upper-Hessenberg matrix
    (C: dgeev with jobvl = jobvr = 'N'). info semantics match LAPACK:
    0 = success, > 0 = the QR iteration failed to converge. */
    let mut a: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n).map(|j| c.LAPACK_A[j * n + i]).collect())
        .collect();
    let info = sundomeigest_hqr(n, &mut a, &mut c.LAPACK_wr, &mut c.LAPACK_wi);

    if info != 0 {
        return SUN_ERR_EXT_FAIL;
    }

    /* order the eigenvalues by their magnitude */
    for i in 0..n {
        c.LAPACK_arr[i][0] = c.LAPACK_wr[i];
        c.LAPACK_arr[i][1] = c.LAPACK_wi[i];
    }

    /* Sort the array (C: qsort with sundomeigest_Compare) */
    c.LAPACK_arr.sort_by(sundomeigest_Compare);

    /* Substitute the ordered eigenvalues back in LAPACK_w* */
    for i in 0..n {
        c.LAPACK_wr[i] = c.LAPACK_arr[i][0];
        c.LAPACK_wi[i] = c.LAPACK_arr[i][1];
    }

    /* Copy the dominant eigenvalue */
    *lambdaR = c.LAPACK_wr[0];
    *lambdaI = c.LAPACK_wi[0];

    SUN_SUCCESS
}

pub fn SUNDomEigEstimator_GetNumIters_Arnoldi(
    dee: &SUNDomEigEstimator,
    num_iters: &mut i64,
) -> SUNErrCode {
    let c = match content(dee) {
        Ok(c) => c,
        Err(e) => return e,
    };
    *num_iters = c.num_iters;
    SUN_SUCCESS
}

pub fn SUNDomEigEstimator_GetNumATimesCalls_Arnoldi(
    dee: &SUNDomEigEstimator,
    num_ATimes: &mut i64,
) -> SUNErrCode {
    let c = match content(dee) {
        Ok(c) => c,
        Err(e) => return e,
    };
    *num_ATimes = c.num_ATimes;
    SUN_SUCCESS
}

pub fn SUNDomEigEstimator_Write_Arnoldi(
    dee: &SUNDomEigEstimator,
    outfile: &mut dyn std::io::Write,
) -> SUNErrCode {
    let c = match content(dee) {
        Ok(c) => c,
        Err(e) => return e,
    };
    let r = (|| -> std::io::Result<()> {
        write!(outfile, "\nArnoldi Iteration SUNDomEigEstimator:\n")?;
        write!(outfile, "Krylov dimension         = {}\n", c.kry_dim)?;
        write!(outfile, "Num. preprocessing iters = {}\n", c.num_warmups)?;
        write!(outfile, "Num. iters               = {}\n", c.num_iters)?;
        write!(outfile, "Num. ATimes calls        = {}\n\n", c.num_ATimes)?;
        Ok(())
    })();
    if r.is_err() {
        return SUN_ERR_ARG_CORRUPT;
    }
    SUN_SUCCESS
}

/* -----------------------------------------------------------------
 * Arnoldi module private functions
 * ----------------------------------------------------------------- */

/// Comparison function for the eigenvalue sort (C sundomeigest_Compare
/// for qsort): descending order of complex magnitude.
fn sundomeigest_Compare(a: &[f64; 2], b: &[f64; 2]) -> std::cmp::Ordering {
    let mag_a = SUNRsqrt(a[0] * a[0] + a[1] * a[1]);
    let mag_b = SUNRsqrt(b[0] * b[0] + b[1] * b[1]);
    let c = (mag_b > mag_a) as i32 - ((mag_b < mag_a) as i32); /* Descending order */
    c.cmp(&0)
}

/// Eigenvalues of a real upper-Hessenberg matrix `a` (n x n, row
/// indexed `a[i][j]`) by the Francis double-shift QR iteration — the
/// classical EISPACK HQR algorithm (the eigenvalues-only kernel that
/// LAPACK's dgeev, called by the C module, applies to the Hessenberg
/// form). `a` is destroyed. Real parts land in `wr`, imaginary parts
/// in `wi`; complex conjugate pairs are stored consecutively with the
/// positive imaginary part first (the dgeev output convention).
///
/// Returns 0 on success; if the iteration for an eigenvalue fails to
/// deflate after 30 shifts, returns a positive count (mirroring
/// LAPACK's info > 0 "failed to compute all the eigenvalues").
fn sundomeigest_hqr(n: usize, a: &mut [Vec<f64>], wr: &mut [f64], wi: &mut [f64]) -> i32 {
    /* Fortran SIGN(a, b): |a| with the sign of b */
    fn sign(a: f64, b: f64) -> f64 {
        if b >= 0.0 {
            a.abs()
        } else {
            -a.abs()
        }
    }

    let ni = n as isize;

    /* Compute the norm of the Hessenberg part of the matrix */
    let mut anorm = ZERO;
    for i in 0..ni {
        let j0 = if i - 1 > 0 { i - 1 } else { 0 };
        for j in j0..ni {
            anorm += SUNRabs(a[i as usize][j as usize]);
        }
    }

    let mut nn: isize = ni - 1;
    let mut t = ZERO;
    /* Search for the next eigenvalue(s) */
    while nn >= 0 {
        let mut its = 0;
        loop {
            /* Look for a single small subdiagonal element to split the
            matrix; a[l][l-1] is declared negligible when adding it to
            the magnitude of its diagonal neighbors changes nothing. */
            let mut l: isize = 0;
            let mut ll = nn;
            while ll >= 1 {
                let mut s =
                    SUNRabs(a[(ll - 1) as usize][(ll - 1) as usize]) + SUNRabs(a[ll as usize][ll as usize]);
                if s == ZERO {
                    s = anorm;
                }
                if SUNRabs(a[ll as usize][(ll - 1) as usize]) + s == s {
                    a[ll as usize][(ll - 1) as usize] = ZERO;
                    l = ll;
                    break;
                }
                ll -= 1;
            }

            let mut x = a[nn as usize][nn as usize];
            if l == nn {
                /* One (real) root found */
                wr[nn as usize] = x + t;
                wi[nn as usize] = ZERO;
                nn -= 1;
            } else {
                let mut y = a[(nn - 1) as usize][(nn - 1) as usize];
                let mut w = a[nn as usize][(nn - 1) as usize] * a[(nn - 1) as usize][nn as usize];
                if l == nn - 1 {
                    /* Two roots found: eigenvalues of the trailing 2x2 block */
                    let p = 0.5 * (y - x);
                    let q = p * p + w;
                    let mut z = SUNRsqrt(SUNRabs(q));
                    x += t;
                    if q >= ZERO {
                        /* ...a real pair */
                        z = p + sign(z, p);
                        wr[(nn - 1) as usize] = x + z;
                        wr[nn as usize] = x + z;
                        if z != ZERO {
                            wr[nn as usize] = x - w / z;
                        }
                        wi[(nn - 1) as usize] = ZERO;
                        wi[nn as usize] = ZERO;
                    } else {
                        /* ...a complex pair (positive imaginary part first) */
                        wr[(nn - 1) as usize] = x + p;
                        wr[nn as usize] = x + p;
                        wi[(nn - 1) as usize] = z;
                        wi[nn as usize] = -z;
                    }
                    nn -= 2;
                } else {
                    /* No roots found yet; continue the iteration */
                    if its == 30 {
                        return (nn + 1) as i32;
                    }
                    if its == 10 || its == 20 {
                        /* Form an exceptional shift */
                        t += x;
                        for i in 0..=nn {
                            a[i as usize][i as usize] -= x;
                        }
                        let s = SUNRabs(a[nn as usize][(nn - 1) as usize])
                            + SUNRabs(a[(nn - 1) as usize][(nn - 2) as usize]);
                        y = 0.75 * s;
                        x = 0.75 * s;
                        w = -0.4375 * s * s;
                    }
                    its += 1;

                    /* Form the (double) shift and look for two consecutive
                    small subdiagonal elements to start the QR sweep at. */
                    let mut m = nn - 2;
                    let mut p = ZERO;
                    let mut q = ZERO;
                    let mut r = ZERO;
                    while m >= l {
                        let z = a[m as usize][m as usize];
                        let rr = x - z;
                        let ss = y - z;
                        p = (rr * ss - w) / a[(m + 1) as usize][m as usize]
                            + a[m as usize][(m + 1) as usize];
                        q = a[(m + 1) as usize][(m + 1) as usize] - z - rr - ss;
                        r = a[(m + 2) as usize][(m + 1) as usize];
                        let s = SUNRabs(p) + SUNRabs(q) + SUNRabs(r);
                        p /= s;
                        q /= s;
                        r /= s;
                        if m == l {
                            break;
                        }
                        let u = SUNRabs(a[m as usize][(m - 1) as usize]) * (SUNRabs(q) + SUNRabs(r));
                        let v = SUNRabs(p)
                            * (SUNRabs(a[(m - 1) as usize][(m - 1) as usize])
                                + SUNRabs(z)
                                + SUNRabs(a[(m + 1) as usize][(m + 1) as usize]));
                        if u + v == v {
                            break;
                        }
                        m -= 1;
                    }

                    for i in (m + 2)..=nn {
                        a[i as usize][(i - 2) as usize] = ZERO;
                        if i != m + 2 {
                            a[i as usize][(i - 3) as usize] = ZERO;
                        }
                    }

                    /* Double QR step on rows l..nn and columns m..nn */
                    for kk in m..=(nn - 1) {
                        if kk != m {
                            p = a[kk as usize][(kk - 1) as usize];
                            q = a[(kk + 1) as usize][(kk - 1) as usize];
                            r = ZERO;
                            if kk != nn - 1 {
                                r = a[(kk + 2) as usize][(kk - 1) as usize];
                            }
                            x = SUNRabs(p) + SUNRabs(q) + SUNRabs(r);
                            if x != ZERO {
                                p /= x;
                                q /= x;
                                r /= x;
                            }
                        }
                        let s = sign(SUNRsqrt(p * p + q * q + r * r), p);
                        if s != ZERO {
                            if kk == m {
                                if l != m {
                                    a[kk as usize][(kk - 1) as usize] =
                                        -a[kk as usize][(kk - 1) as usize];
                                }
                            } else {
                                a[kk as usize][(kk - 1) as usize] = -s * x;
                            }
                            p += s;
                            x = p / s;
                            y = q / s;
                            let z = r / s;
                            q /= p;
                            r /= p;
                            /* Row modification */
                            for j in kk..=nn {
                                let mut pp = a[kk as usize][j as usize]
                                    + q * a[(kk + 1) as usize][j as usize];
                                if kk != nn - 1 {
                                    pp += r * a[(kk + 2) as usize][j as usize];
                                    a[(kk + 2) as usize][j as usize] -= pp * z;
                                }
                                a[(kk + 1) as usize][j as usize] -= pp * y;
                                a[kk as usize][j as usize] -= pp * x;
                            }
                            /* Column modification */
                            let mmin = if nn < kk + 3 { nn } else { kk + 3 };
                            for i in l..=mmin {
                                let mut pp = x * a[i as usize][kk as usize]
                                    + y * a[i as usize][(kk + 1) as usize];
                                if kk != nn - 1 {
                                    pp += z * a[i as usize][(kk + 2) as usize];
                                    a[i as usize][(kk + 2) as usize] -= pp * r;
                                }
                                a[i as usize][(kk + 1) as usize] -= pp * q;
                                a[i as usize][kk as usize] -= pp;
                            }
                        }
                    }
                }
            }
            if l >= nn - 1 {
                break;
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sundials_context::SUNContext;
    use crate::sundials_domeigestimator::*;

    /// Dense matrix-vector product ATimes closure (owns its matrix).
    fn dense_atimes(a: Vec<Vec<f64>>) -> Box<ATimesFn<'static>> {
        Box::new(move |v: &NVector, av: &mut NVector| -> i32 {
            let n = a.len();
            for i in 0..n {
                let mut s = 0.0;
                for j in 0..n {
                    s += a[i][j] * v.data[j];
                }
                av.data[i] = s;
            }
            0
        })
    }

    fn sorted(mut v: Vec<f64>) -> Vec<f64> {
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        v
    }

    #[test]
    fn hqr_companion_matrix_real_roots() {
        /* Companion matrix of (x-1)(x-2)(x-3)(x-4)
           = x^4 - 10x^3 + 35x^2 - 50x + 24 (upper Hessenberg) */
        let mut a = vec![
            vec![10.0, -35.0, 50.0, -24.0],
            vec![1.0, 0.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0, 0.0],
            vec![0.0, 0.0, 1.0, 0.0],
        ];
        let mut wr = vec![0.0; 4];
        let mut wi = vec![0.0; 4];
        assert_eq!(sundomeigest_hqr(4, &mut a, &mut wr, &mut wi), 0);
        let wr = sorted(wr);
        for (got, want) in wr.iter().zip([1.0, 2.0, 3.0, 4.0]) {
            assert!((got - want).abs() < 1e-10, "wr = {:?}", wr);
        }
        for x in wi {
            assert!(x.abs() < 1e-10);
        }
    }

    #[test]
    fn hqr_complex_pair() {
        /* Companion matrix of x^2 + 1: eigenvalues +/- i,
           positive imaginary part first */
        let mut a = vec![vec![0.0, -1.0], vec![1.0, 0.0]];
        let mut wr = vec![0.0; 2];
        let mut wi = vec![0.0; 2];
        assert_eq!(sundomeigest_hqr(2, &mut a, &mut wr, &mut wi), 0);
        assert!(wr[0].abs() < 1e-14 && wr[1].abs() < 1e-14);
        assert!((wi[0] - 1.0).abs() < 1e-14);
        assert!((wi[1] + 1.0).abs() < 1e-14);
    }

    #[test]
    fn constructor_defaults() {
        let sunctx = SUNContext::default();
        let q = NVector::from_slice(&[1.0, 2.0, 3.0]);
        /* kry_dim < 3 falls back onto the default */
        let dee = SUNDomEigEstimator_Arnoldi(&q, 2, &sunctx);
        match &dee {
            SUNDomEigEstimator::Arnoldi(c) => {
                assert_eq!(c.kry_dim, 3);
                assert_eq!(c.num_warmups, 100);
                assert_eq!(c.V.len(), 4); /* kry_dim + 1 clones */
                assert_eq!(c.q.data, vec![1.0, 2.0, 3.0]);
                        assert!(c.Hes.is_empty()); /* allocated by Initialize */
            }
            _ => panic!(),
        }
        assert_eq!(SUNDomEigEstimator_Destroy(dee), SUN_SUCCESS);
    }

    #[test]
    fn estimate_diag_1_2_10_via_dispatch() {
        let sunctx = SUNContext::default();
        let q = NVector::from_slice(&[1.0, 1.0, 1.0]);
        let mut dee = SUNDomEigEstimator_Arnoldi(&q, 3, &sunctx);
        let a = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 2.0, 0.0],
            vec![0.0, 0.0, 10.0],
        ];
        let mut atimes = dense_atimes(a);
        /* keep the Krylov basis exact: no preprocessing iterations */
        assert_eq!(SUNDomEigEstimator_SetNumPreprocessIters(&mut dee, 0), SUN_SUCCESS);
        assert_eq!(SUNDomEigEstimator_Initialize(&mut dee), SUN_SUCCESS);

        let (mut lr, mut li) = (0.0, -1.0);
        assert_eq!(SUNDomEigEstimator_Estimate(&mut dee, &mut atimes, &mut lr, &mut li), SUN_SUCCESS);
        /* kry_dim = n = 3, so the Hessenberg matrix is similar to A:
        the dominant eigenvalue is reproduced (up to roundoff) */
        assert!((lr - 10.0).abs() < 1e-8, "lambdaR = {}", lr);
        assert!(li.abs() < 1e-8, "lambdaI = {}", li);

        /* one ATimes call per Krylov vector, none for warmups */
        let (mut ni, mut na) = (0i64, 0i64);
        assert_eq!(SUNDomEigEstimator_GetNumIters(&dee, &mut ni), SUN_SUCCESS);
        assert_eq!(SUNDomEigEstimator_GetNumATimesCalls(&dee, &mut na), SUN_SUCCESS);
        assert_eq!(ni, 3);
        assert_eq!(na, 3);
    }

    #[test]
    fn estimate_complex_dominant_pair() {
        /* 2D rotation-plus-scaling embedded in 3D:
        eigenvalues {2i, -2i, 1}, dominant pair magnitude 2 */
        let sunctx = SUNContext::default();
        let q = NVector::from_slice(&[1.0, 1.0, 1.0]);
        let mut dee = SUNDomEigEstimator_Arnoldi(&q, 3, &sunctx);
        let a = vec![
            vec![0.0, -2.0, 0.0],
            vec![2.0, 0.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let mut atimes = dense_atimes(a);
        assert_eq!(SUNDomEigEstimator_SetNumPreprocessIters(&mut dee, 0), SUN_SUCCESS);
        assert_eq!(SUNDomEigEstimator_Initialize(&mut dee), SUN_SUCCESS);

        let (mut lr, mut li) = (0.0, 0.0);
        assert_eq!(SUNDomEigEstimator_Estimate(&mut dee, &mut atimes, &mut lr, &mut li), SUN_SUCCESS);
        assert!(lr.abs() < 1e-8, "lambdaR = {}", lr);
        assert!((li - 2.0).abs() < 1e-8, "lambdaI = {}", li);

        /* GetRes falls back onto the base-class default (op absent) */
        let mut res = -1.0;
        assert_eq!(SUNDomEigEstimator_GetRes(&dee, &mut res), SUN_SUCCESS);
        assert_eq!(res, 0.0);
    }

    #[test]
    fn warmups_count_and_reestimate() {
        let sunctx = SUNContext::default();
        let q = NVector::from_slice(&[1.0, 1.0, 1.0, 1.0]);
        let mut dee = SUNDomEigEstimator_Arnoldi(&q, 3, &sunctx);
        let a = vec![
            vec![4.0, 1.0, 0.0, 0.0],
            vec![1.0, 3.0, 0.0, 0.0],
            vec![0.0, 0.0, 2.0, 0.0],
            vec![0.0, 0.0, 0.0, 1.0],
        ];
        let mut atimes = dense_atimes(a);
        assert_eq!(SUNDomEigEstimator_SetNumPreprocessIters(&mut dee, 2), SUN_SUCCESS);
        assert_eq!(SUNDomEigEstimator_Initialize(&mut dee), SUN_SUCCESS);

        let (mut lr, mut li) = (0.0, 0.0);
        assert_eq!(SUNDomEigEstimator_Estimate(&mut dee, &mut atimes, &mut lr, &mut li), SUN_SUCCESS);
        /* dominant eigenvalue of [[4,1],[1,3]] block: (7+sqrt(5))/2.
        kry_dim (3) < n (4), so the Ritz value is an approximation */
        let want = (7.0 + 5.0f64.sqrt()) / 2.0;
        assert!((lr - want).abs() < 1e-3, "lambdaR = {} want {}", lr, want);

        /* counters: 2 warmups + kry_dim Krylov steps, reset per call */
        let (mut ni, mut na) = (0i64, 0i64);
        assert_eq!(SUNDomEigEstimator_GetNumIters(&dee, &mut ni), SUN_SUCCESS);
        assert_eq!(SUNDomEigEstimator_GetNumATimesCalls(&dee, &mut na), SUN_SUCCESS);
        assert_eq!(ni, 5);
        assert_eq!(na, 5);

        /* a second Estimate call resets the counters */
        assert_eq!(SUNDomEigEstimator_Estimate(&mut dee, &mut atimes, &mut lr, &mut li), SUN_SUCCESS);
        assert_eq!(SUNDomEigEstimator_GetNumIters(&dee, &mut ni), SUN_SUCCESS);
        assert_eq!(ni, 5);
    }

    #[test]
    fn estimate_before_initialize_fails() {
        let sunctx = SUNContext::default();
        let q = NVector::from_slice(&[1.0, 1.0, 1.0]);
        let mut dee = SUNDomEigEstimator_Arnoldi(&q, 3, &sunctx);
        let a = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 2.0, 0.0],
            vec![0.0, 0.0, 3.0],
        ];
        let mut atimes = dense_atimes(a);
        let (mut lr, mut li) = (0.0, 0.0);
        assert_eq!(
            SUNDomEigEstimator_Estimate(&mut dee, &mut atimes, &mut lr, &mut li),
            SUN_ERR_ARG_CORRUPT
        );
    }

    #[test]
    fn write_output() {
        let sunctx = SUNContext::default();
        let q = NVector::from_slice(&[1.0, 1.0, 1.0]);
        let mut dee = SUNDomEigEstimator_Arnoldi(&q, 3, &sunctx);
        let a = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 2.0, 0.0],
            vec![0.0, 0.0, 10.0],
        ];
        let mut atimes = dense_atimes(a);
        assert_eq!(SUNDomEigEstimator_SetNumPreprocessIters(&mut dee, 2), SUN_SUCCESS);
        assert_eq!(SUNDomEigEstimator_Initialize(&mut dee), SUN_SUCCESS);
        let (mut lr, mut li) = (0.0, 0.0);
        assert_eq!(SUNDomEigEstimator_Estimate(&mut dee, &mut atimes, &mut lr, &mut li), SUN_SUCCESS);

        let mut buf: Vec<u8> = Vec::new();
        assert_eq!(SUNDomEigEstimator_Write(&dee, &mut buf), SUN_SUCCESS);
        let s = String::from_utf8(buf).unwrap();
        assert!(s.starts_with("\nArnoldi Iteration SUNDomEigEstimator:\n"));
        assert!(s.contains("Krylov dimension         = 3\n"));
        assert!(s.contains("Num. preprocessing iters = 2\n"));
        assert!(s.contains("Num. iters               = 5\n"));
        assert!(s.ends_with("Num. ATimes calls        = 5\n\n"));
    }

    #[test]
    fn wrong_variant_returns_incompatible() {
        let sunctx = SUNContext::default();
        let q = NVector::from_slice(&[1.0, 1.0]);
        let mut power = crate::sundomeigest_power::SUNDomEigEstimator_Power(&q, 10, 0.1, &sunctx);
        let mut ni = 0i64;
        assert_eq!(
            SUNDomEigEstimator_GetNumIters_Arnoldi(&power, &mut ni),
            SUN_ERR_ARG_INCOMPATIBLE
        );
        assert_eq!(
            SUNDomEigEstimator_Initialize_Arnoldi(&mut power),
            SUN_ERR_ARG_INCOMPATIBLE
        );
    }
}
