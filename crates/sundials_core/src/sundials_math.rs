/* -----------------------------------------------------------------
 * Translated from src/sundials/sundials_math.c and
 * include/sundials/sundials_math.h (SUNDIALS 7.7.0).
 * -----------------------------------------------------------------*/
use crate::sundials_types::*;

#[inline]
pub fn SUNMIN(a: f64, b: f64) -> f64 {
    if a < b { a } else { b }
}
#[inline]
pub fn SUNMAX(a: f64, b: f64) -> f64 {
    if a > b { a } else { b }
}
#[inline]
pub fn SUNSQR(a: f64) -> f64 {
    a * a
}
#[inline]
pub fn SUNRabs(x: f64) -> f64 {
    x.abs()
}
#[inline]
pub fn SUNRsqrt(x: f64) -> f64 {
    if x <= 0.0 { 0.0 } else { x.sqrt() }
}
#[inline]
pub fn SUNRexp(x: f64) -> f64 {
    x.exp()
}
#[inline]
pub fn SUNRlog(x: f64) -> f64 {
    x.ln()
}
#[inline]
pub fn SUNRceil(x: f64) -> f64 {
    x.ceil()
}
#[inline]
pub fn SUNRround(x: f64) -> f64 {
    x.round()
}

/// SUNRpowerI: base^exponent via repeated multiplication (kept loop-exact
/// with the C source so results match bit-for-bit).
pub fn SUNRpowerI(base: f64, exponent: i32) -> f64 {
    let mut prod: f64 = 1.0;
    let expt = exponent.unsigned_abs();
    for _ in 1..=expt {
        prod *= base;
    }
    if exponent < 0 {
        prod = 1.0 / prod;
    }
    prod
}

/// SUNRpowerR: base^exponent for real exponent; 0 for negative base.
pub fn SUNRpowerR(base: f64, exponent: f64) -> f64 {
    if base <= 0.0 {
        return 0.0;
    }
    base.powf(exponent)
}

/// SUNRsamesign(x, y): true if x and y share the same sign bit.
#[inline]
pub fn SUNRsamesign(x: f64, y: f64) -> bool {
    x.is_sign_negative() == y.is_sign_negative()
}

/// SUNRdifferentsign(x, y): true if x and y have different sign bits.
#[inline]
pub fn SUNRdifferentsign(x: f64, y: f64) -> bool {
    !SUNRsamesign(x, y)
}

/// SUNRCompare: returns SUNTRUE if a and b are NOT equal to within
/// 10*unit-roundoff relative tolerance (note the inverted convention).
pub fn SUNRCompare(a: f64, b: f64) -> bool {
    SUNRCompareTol(a, b, 10.0 * SUN_UNIT_ROUNDOFF)
}

pub fn SUNRCompareTol(a: f64, b: f64, tol: f64) -> bool {
    if a == b {
        return SUNFALSE;
    }
    let diff = SUNRabs(a - b);
    let norm = SUNMIN(SUNRabs(a + b), SUN_BIG_REAL);
    // C uses !isless(diff, max(10*uround, tol*norm)) to be NaN-safe.
    !(diff < SUNMAX(10.0 * SUN_UNIT_ROUNDOFF, tol * norm))
}

/// SUNStrToReal (sundials_math.c): strtod equivalent.
pub fn SUNStrToReal(s: &str) -> f64 {
    s.trim().parse::<f64>().unwrap_or(0.0)
}
