/* -----------------------------------------------------------------
 * Translated from src/sundials/sundials_iterative.c and
 * src/sundials/sundials_iterative_impl.h (SUNDIALS 7.7.0).
 *
 * Helper routines shared by the iterative linear solvers:
 * SUNModifiedGS / SUNClassicalGS (Gram-Schmidt orthogonalization
 * with conditional reorthogonalization), SUNQRfact (Givens-rotation
 * QR update of a Hessenberg matrix) and SUNQRsol (least-squares
 * solve with the factored matrix).
 *
 * Signature adaptations for safe Rust:
 *  - the Krylov basis `v` is a `&mut [NVector]` (C: `N_Vector*`);
 *  - the Hessenberg matrix `h` is row-major `Vec<Vec<f64>>` indexed
 *    `h[i][j]` exactly like the C `h[i][j]`;
 *  - SUNClassicalGS keeps the scalar workspace `stemp` but drops the
 *    C `vtemp` argument: that array only stored N_Vector *aliases*
 *    feeding the fused N_VLinearCombination kernel, whose serial
 *    implementation is reproduced inline here with identical
 *    floating-point operation order.
 *
 * The SUNQRAddFn type and SUNQRData workspace struct (used by
 * KINSOL's Anderson-acceleration orthogonalization options) are
 * defined at the end of this module; the SUNQRAdd_* routine
 * implementations land together with the kinsol_rs port.
 * -----------------------------------------------------------------*/
use crate::nvector_serial::{NVector, N_VDotProd};
use crate::sundials_errors::{SUNErrCode, SUN_SUCCESS};
use crate::sundials_math::{SUNRsqrt, SUNSQR};

const FACTOR: f64 = 1000.0;
const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/* -----------------------------------------------------------------
 * Function : SUNModifiedGS
 *
 * Orthogonalize v[k] against v[max(k-p,0)] .. v[k-1] with modified
 * Gram-Schmidt, filling column k-1 of h and returning the new norm
 * of v[k] in new_vk_norm; reorthogonalizes when severe cancellation
 * is detected.
 * -----------------------------------------------------------------*/
pub fn SUNModifiedGS(
    v: &mut [NVector],
    h: &mut Vec<Vec<f64>>,
    k: i32,
    p: i32,
    new_vk_norm: &mut f64,
) -> SUNErrCode {
    let ku = k as usize;
    let k_minus_1 = ku - 1;
    let i0 = if k - p > 0 { (k - p) as usize } else { 0 };

    let (head, tail) = v.split_at_mut(ku);
    let vk = &mut tail[0];

    let mut vk_norm = N_VDotProd(vk, vk);
    vk_norm = SUNRsqrt(vk_norm);

    /* Perform modified Gram-Schmidt */

    for i in i0..ku {
        h[i][k_minus_1] = N_VDotProd(&head[i], vk);
        /* v[k] = v[k] - h[i][k-1]*v[i] */
        vk.linear_sum_with(ONE, -h[i][k_minus_1], &head[i]);
    }

    /* Compute the norm of the new vector at v[k] */

    *new_vk_norm = N_VDotProd(vk, vk);
    *new_vk_norm = SUNRsqrt(*new_vk_norm);

    /* If the norm of the new vector at v[k] is less than
       FACTOR (== 1000) times unit roundoff times the norm of the
       input vector v[k], then the vector will be reorthogonalized
       in order to ensure that nonorthogonality is not being masked
       by a very small vector length. */

    let temp = FACTOR * vk_norm;
    if temp + *new_vk_norm != temp {
        return SUN_SUCCESS;
    }

    let mut new_norm_2 = ZERO;

    for i in i0..ku {
        let new_product = N_VDotProd(&head[i], vk);
        let temp = FACTOR * h[i][k_minus_1];
        if temp + new_product == temp {
            continue;
        }
        h[i][k_minus_1] += new_product;
        vk.linear_sum_with(ONE, -new_product, &head[i]);
        new_norm_2 += SUNSQR(new_product);
    }

    if new_norm_2 != ZERO {
        let new_product = SUNSQR(*new_vk_norm) - new_norm_2;
        *new_vk_norm = if new_product > ZERO {
            SUNRsqrt(new_product)
        } else {
            ZERO
        };
    }

    SUN_SUCCESS
}

/* -----------------------------------------------------------------
 * Function : SUNClassicalGS
 *
 * Classical Gram-Schmidt orthogonalization of v[k] against
 * v[max(k-p,0)] .. v[k-1], with one conditional reorthogonalization
 * pass. `stemp` is scalar workspace of length >= k - max(k-p,0) + 1.
 * -----------------------------------------------------------------*/
pub fn SUNClassicalGS(
    v: &mut [NVector],
    h: &mut Vec<Vec<f64>>,
    k: i32,
    p: i32,
    new_vk_norm: &mut f64,
    stemp: &mut [f64],
) -> SUNErrCode {
    let ku = k as usize;
    let k_minus_1 = ku - 1;
    let i0 = if k - p > 0 { (k - p) as usize } else { 0 };
    let kmi0 = ku - i0;

    let (head, tail) = v.split_at_mut(ku);
    let vk = &mut tail[0];

    /* Perform Classical Gram-Schmidt */

    /* N_VDotProdMulti(k - i0 + 1, v[k], v + i0, stemp) */
    for j in 0..kmi0 {
        stemp[j] = N_VDotProd(vk, &head[i0 + j]);
    }
    stemp[kmi0] = N_VDotProd(vk, vk);

    let vk_norm = SUNRsqrt(stemp[kmi0]);
    for i in (0..kmi0).rev() {
        h[i][k_minus_1] = stemp[i];
        stemp[i + 1] = -stemp[i];
        /* C: vtemp[i + 1] = v[i] (alias only, handled below) */
    }
    stemp[0] = ONE;
    /* C: vtemp[0] = v[k] */

    /* N_VLinearCombination(k - i0 + 1, stemp, vtemp, v[k]) with
       vtemp = [v[k], v[0], ..., v[k-i0-1]]; the serial kernel takes
       the (X[0] == z && c[0] == ONE) path: z += sum_j c[j]*X[j],
       accumulated one vector at a time. */
    for j in 1..=kmi0 {
        let c = stemp[j];
        for (zi, xi) in vk.data.iter_mut().zip(&head[j - 1].data) {
            *zi += c * *xi;
        }
    }

    /* Compute the norm of the new vector at v[k] */

    *new_vk_norm = SUNRsqrt(N_VDotProd(vk, vk));

    /* Reorthogonalize if necessary */

    if (FACTOR * *new_vk_norm) < vk_norm {
        /* N_VDotProdMulti(k - i0, v[k], v + i0, stemp + 1) */
        for j in 0..kmi0 {
            stemp[j + 1] = N_VDotProd(vk, &head[i0 + j]);
        }

        stemp[0] = ONE;
        /* C: vtemp[0] = v[k] */
        for i in i0..ku {
            h[i][k_minus_1] += stemp[i - i0 + 1];
            stemp[i - i0 + 1] = -stemp[i - i0 + 1];
            /* C: vtemp[i - i0 + 1] = v[i - i0] */
        }

        /* N_VLinearCombination(k + 1, stemp, vtemp, v[k]) with
           vtemp = [v[k], v[0], ..., v[k-1]] (i0 == 0 whenever this
           routine is called from SPGMR/SPFGMR, since p = maxl >= k) */
        for j in 1..=ku {
            let c = stemp[j];
            for (zi, xi) in vk.data.iter_mut().zip(&head[j - 1].data) {
                *zi += c * *xi;
            }
        }

        *new_vk_norm = SUNRsqrt(N_VDotProd(vk, vk));
    }

    SUN_SUCCESS
}

/* Compute the Givens rotation components c and s from the pair
   (temp1, temp2) exactly as both branches of the C SUNQRfact do. */
#[inline]
fn givens_cs(temp1: f64, temp2: f64) -> (f64, f64) {
    let c;
    let s;
    if temp2 == ZERO {
        c = ONE;
        s = ZERO;
    } else if temp2.abs() >= temp1.abs() {
        let temp3 = temp1 / temp2;
        s = -ONE / SUNRsqrt(ONE + SUNSQR(temp3));
        c = -s * temp3;
    } else {
        let temp3 = temp2 / temp1;
        c = ONE / SUNRsqrt(ONE + SUNSQR(temp3));
        s = -c * temp3;
    }
    (c, s)
}

/* -----------------------------------------------------------------
 * Function : SUNQRfact
 *
 * QR factorization (job == 0) or single-column QR update (job != 0)
 * of the (n+1) x n Hessenberg matrix h via Givens rotations stored
 * in q (c,s pairs). Returns 0 on success, k+1 if a zero appears on
 * the diagonal of R at position k.
 * -----------------------------------------------------------------*/
pub fn SUNQRfact(n: i32, h: &mut Vec<Vec<f64>>, q: &mut [f64], job: i32) -> i32 {
    let mut code = 0;

    match job {
        0 => {
            /* Compute a new factorization of H */

            for k in 0..n {
                let ku = k as usize;

                /* Multiply column k by the previous k-1 Givens rotations */

                for j in 0..(k - 1) {
                    let ju = j as usize;
                    let i = 2 * ju;
                    let temp1 = h[ju][ku];
                    let temp2 = h[ju + 1][ku];
                    let c = q[i];
                    let s = q[i + 1];
                    h[ju][ku] = c * temp1 - s * temp2;
                    h[ju + 1][ku] = s * temp1 + c * temp2;
                }

                /* Compute the Givens rotation components c and s */

                let q_ptr = 2 * ku;
                let temp1 = h[ku][ku];
                let temp2 = h[ku + 1][ku];
                let (c, s) = givens_cs(temp1, temp2);
                q[q_ptr] = c;
                q[q_ptr + 1] = s;
                h[ku][ku] = c * temp1 - s * temp2;
                if h[ku][ku] == ZERO {
                    code = k + 1;
                }
            }
        }
        _ => {
            /* Update the factored H to which a new column has been added */

            let n_minus_1 = (n - 1) as usize;

            /* Multiply the new column by the previous n-1 Givens rotations */

            for k in 0..n_minus_1 {
                let i = 2 * k;
                let temp1 = h[k][n_minus_1];
                let temp2 = h[k + 1][n_minus_1];
                let c = q[i];
                let s = q[i + 1];
                h[k][n_minus_1] = c * temp1 - s * temp2;
                h[k + 1][n_minus_1] = s * temp1 + c * temp2;
            }

            /* Compute new Givens rotation and multiply it times the last two
               entries in the new column of H.  Note that the second entry of
               this product will be 0, so it is not necessary to compute it. */

            let temp1 = h[n_minus_1][n_minus_1];
            let temp2 = h[n as usize][n_minus_1];
            let (c, s) = givens_cs(temp1, temp2);
            let q_ptr = 2 * n_minus_1;
            q[q_ptr] = c;
            q[q_ptr + 1] = s;
            h[n_minus_1][n_minus_1] = c * temp1 - s * temp2;
            if h[n_minus_1][n_minus_1] == ZERO {
                code = n;
            }
        }
    }

    code
}

/* -----------------------------------------------------------------
 * Function : SUNQRsol
 *
 * Solve the least-squares problem min ||b - H*y|| with the QR
 * factors produced by SUNQRfact; the solution overwrites b[0..n].
 * Returns 0 on success, k+1 if R has a zero diagonal entry at k.
 * -----------------------------------------------------------------*/
pub fn SUNQRsol(n: i32, h: &[Vec<f64>], q: &[f64], b: &mut [f64]) -> i32 {
    let mut code = 0;
    let nu = n as usize;

    /* Compute Q*b */

    for k in 0..nu {
        let q_ptr = 2 * k;
        let c = q[q_ptr];
        let s = q[q_ptr + 1];
        let temp1 = b[k];
        let temp2 = b[k + 1];
        b[k] = c * temp1 - s * temp2;
        b[k + 1] = s * temp1 + c * temp2;
    }

    /* Solve  R*x = Q*b */

    for k in (0..nu).rev() {
        if h[k][k] == ZERO {
            code = (k + 1) as i32;
            break;
        }
        b[k] /= h[k][k];
        for i in 0..k {
            b[i] -= b[k] * h[i][k];
        }
    }

    code
}

/* -----------------------------------------------------------------
 * Type : SUNQRData (src/sundials/sundials_iterative_impl.h)
 *
 * Holds temporary workspace vectors and a sunrealtype array for a
 * SUNQRAddFn. In C the N_Vectors and array are created by the
 * routine calling a SUNQRAdd function; here they are owned.
 * -----------------------------------------------------------------*/
#[derive(Debug, Clone, Default)]
pub struct SUNQRData {
    pub vtemp: NVector,
    pub vtemp2: NVector,
    pub temp_array: Vec<f64>,
}

/* -----------------------------------------------------------------
 * Type : SUNQRAddFn (include/sundials/sundials_iterative.h)
 *
 * Updates the QR factorization (Q, R) with the input vector df.
 *   m    : number of vectors already in the QR factorization
 *   mMax : maximum number of vectors in the QR factorization
 * C: int (*SUNQRAddFn)(N_Vector* Q_1d, sunrealtype* R_1d, N_Vector f,
 *                      int m, int mMax, void* QR_data)
 * The void* QR_data becomes &mut SUNQRData (workspace convention).
 * -----------------------------------------------------------------*/
pub type SUNQRAddFn = fn(
    Q: &mut [NVector],
    R: &mut [f64],
    df: &NVector,
    m: i32,
    mMax: i32,
    qr_data: &mut SUNQRData,
) -> i32;

#[cfg(test)]
mod tests {
    use super::*;

    /// Orthogonalize a second vector against the first with MGS and check
    /// orthogonality plus the returned norm.
    #[test]
    fn modified_gs_orthogonalizes() {
        let mut v = vec![
            NVector::from_slice(&[1.0, 0.0, 0.0]),
            NVector::from_slice(&[1.0, 1.0, 0.0]),
        ];
        let mut h = vec![vec![0.0; 2]; 3];
        let mut norm = 0.0;
        let flag = SUNModifiedGS(&mut v, &mut h, 1, 2, &mut norm);
        assert_eq!(flag, SUN_SUCCESS);
        assert!((h[0][0] - 1.0).abs() < 1e-15);
        assert!(N_VDotProd(&v[0], &v[1]).abs() < 1e-15);
        assert!((norm - 1.0).abs() < 1e-15);
    }

    /// QR-factor a small 3x2 Hessenberg matrix column by column and solve;
    /// verify the least-squares solution reproduces the projected rhs.
    #[test]
    fn qrfact_qrsol_roundtrip() {
        // H = [[2, 1], [1, 3], [0, 1]] built up one column at a time,
        // exactly as SPGMR does (job = l).
        let mut h = vec![vec![0.0; 2], vec![0.0; 2], vec![0.0; 2]];
        let mut q = vec![0.0; 4];
        h[0][0] = 2.0;
        h[1][0] = 1.0;
        assert_eq!(SUNQRfact(1, &mut h, &mut q, 0), 0);
        h[0][1] = 1.0;
        h[1][1] = 3.0;
        h[2][1] = 1.0;
        assert_eq!(SUNQRfact(2, &mut h, &mut q, 1), 0);

        // Solve min ||g - H y|| for g = (1, 0, 0): the residual of the
        // normal equations H^T (g - H y) must vanish.
        let mut b = [1.0, 0.0, 0.0];
        assert_eq!(SUNQRsol(2, &h, &q, &mut b), 0);
        let y = [b[0], b[1]];
        let hy = [
            2.0 * y[0] + 1.0 * y[1],
            1.0 * y[0] + 3.0 * y[1],
            1.0 * y[1],
        ];
        let r = [1.0 - hy[0], -hy[1], -hy[2]];
        let nt_r0 = 2.0 * r[0] + 1.0 * r[1];
        let nt_r1 = 1.0 * r[0] + 3.0 * r[1] + 1.0 * r[2];
        assert!(nt_r0.abs() < 1e-14 && nt_r1.abs() < 1e-14);
    }
}
