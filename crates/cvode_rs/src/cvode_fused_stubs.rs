/* -----------------------------------------------------------------
 * Translated from src/cvode/cvode_fused_stubs.c (CVODE 7.7.0).
 * This file implements fused stub kernels for CVODE.
 *
 * They are never called in this build (cv_usefused is always false,
 * matching a C build without SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS),
 * but are kept as correct, safe-Rust translations. Aliased C vector
 * operations (output == input) use the in-place NVector methods.
 * -----------------------------------------------------------------*/
use crate::nvector_serial::*;

const ZERO: f64 = 0.0;
const PT1: f64 = 0.1;
const FRACT: f64 = 0.1;
const ONEPT5: f64 = 1.50;
const ONE: f64 = 1.0;

/*
 * -----------------------------------------------------------------
 * Compute the ewt vector when the tol type is CV_SS.
 * -----------------------------------------------------------------
 */

pub fn cvEwtSetSS_fused(
    atolmin0: bool,
    reltol: f64,
    Sabstol: f64,
    ycur: &NVector,
    tempv: &mut NVector,
    weight: &mut NVector,
) -> i32 {
    N_VAbs(ycur, tempv);
    tempv.scale_inplace(reltol); /* N_VScale(reltol, tempv, tempv) */
    tempv.add_const_inplace(Sabstol); /* N_VAddConst(tempv, Sabstol, tempv) */
    if atolmin0 {
        if N_VMin(tempv) <= ZERO {
            return -1;
        }
    }
    N_VInv(tempv, weight);
    0
}

/*
 * -----------------------------------------------------------------
 * Compute the ewt vector when the tol type is CV_SV.
 * -----------------------------------------------------------------
 */

pub fn cvEwtSetSV_fused(
    atolmin0: bool,
    reltol: f64,
    Vabstol: &NVector,
    ycur: &NVector,
    tempv: &mut NVector,
    weight: &mut NVector,
) -> i32 {
    N_VAbs(ycur, tempv);
    /* N_VLinearSum(reltol, tempv, ONE, Vabstol, tempv) */
    tempv.linear_sum_with(reltol, ONE, Vabstol);
    if atolmin0 {
        if N_VMin(tempv) <= ZERO {
            return -1;
        }
    }
    N_VInv(tempv, weight);
    0
}

/*
 * -----------------------------------------------------------------
 * Determine if the constraints of the problem are satisfied by
 * the proposed step.
 * -----------------------------------------------------------------
 */

pub fn cvCheckConstraints_fused(
    c: &NVector,
    ewt: &NVector,
    y: &NVector,
    mm: &NVector,
    tmp: &mut NVector,
    save: &mut NVector,
) -> i32 {
    N_VCompare(ONEPT5, c, tmp); /* a[i]=1 when |c[i]|=2  */
    tmp.prod_with(c); /* a * c                 */
    tmp.div_with(ewt); /* a * c * wt            */
    N_VScale(-PT1, tmp, save);
    tmp.linear_sum_with(-PT1, ONE, y); /* y - 0.1 * a * c * wt  */
    tmp.prod_with(mm); /* v = mm*(y-0.1*a*c*wt) */
    0
}

/*
 * -----------------------------------------------------------------
 * Compute the nonlinear residual.
 * -----------------------------------------------------------------
 */

pub fn cvNlsResid_fused(
    rl1: f64,
    ngamma: f64,
    zn1: &NVector,
    ycor: &NVector,
    ftemp: &NVector,
    res: &mut NVector,
) -> i32 {
    N_VLinearSum(rl1, zn1, ONE, ycor, res);
    /* N_VLinearSum(ngamma, ftemp, ONE, res, res) */
    res.linear_sum_with(ONE, ngamma, ftemp);
    0
}

/*
 * -----------------------------------------------------------------
 * Form y with perturbation = FRACT*(func. iter. correction)
 * -----------------------------------------------------------------
 */

pub fn cvDiagSetup_formY(
    h: f64,
    r: f64,
    fpred: &NVector,
    zn1: &NVector,
    ypred: &NVector,
    ftemp: &mut NVector,
    y: &mut NVector,
) -> i32 {
    N_VLinearSum(h, fpred, -ONE, zn1, ftemp);
    N_VLinearSum(r, ftemp, ONE, ypred, y);
    0
}

/*
 * -----------------------------------------------------------------
 * Construct M = I - gamma*J with J = diag(deltaf_i/deltay_i)
 * protecting against deltay_i being at roundoff level.
 * -----------------------------------------------------------------
 */

pub fn cvDiagSetup_buildM(
    _fract: f64,
    uround: f64,
    h: f64,
    ftemp: &NVector,
    fpred: &NVector,
    ewt: &NVector,
    bit: &mut NVector,
    bitcomp: &mut NVector,
    y: &mut NVector,
    M: &mut NVector,
) -> i32 {
    /* N_VLinearSum(ONE, M, -ONE, fpred, M) */
    M.linear_sum_with(ONE, -ONE, fpred);
    /* N_VLinearSum(FRACT, ftemp, -h, M, M) */
    M.linear_sum_with(-h, FRACT, ftemp);
    N_VProd(ftemp, ewt, y);
    /* Protect against deltay_i being at roundoff level */
    N_VCompare(uround, y, bit);
    N_VAddConst(bit, -ONE, bitcomp);
    N_VProd(ftemp, bit, y);
    /* N_VLinearSum(FRACT, y, -ONE, bitcomp, y) */
    y.linear_sum_with(FRACT, -ONE, bitcomp);
    M.div_with(y); /* N_VDiv(M, y, M) */
    M.prod_with(bit); /* N_VProd(M, bit, M) */
    /* N_VLinearSum(ONE, M, -ONE, bitcomp, M) */
    M.linear_sum_with(ONE, -ONE, bitcomp);
    0
}

/*
 * -----------------------------------------------------------------
 *  Update M with changed gamma so that M = I - gamma*J.
 * -----------------------------------------------------------------
 */

pub fn cvDiagSolve_updateM(r: f64, M: &mut NVector) -> i32 {
    M.invert_inplace(); /* N_VInv(M, M) */
    M.add_const_inplace(-ONE); /* N_VAddConst(M, -ONE, M) */
    M.scale_inplace(r); /* N_VScale(r, M, M) */
    M.add_const_inplace(ONE); /* N_VAddConst(M, ONE, M) */
    0
}
