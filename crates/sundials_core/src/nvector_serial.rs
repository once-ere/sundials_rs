/* -----------------------------------------------------------------
 * Translated from src/nvector/serial/nvector_serial.c and
 * include/nvector/nvector_serial.h (SUNDIALS 7.7.0).
 *
 * The C serial N_Vector is a length + data pointer plus an ops table;
 * here it is a plain owned buffer. Ops used by CVODE are provided as
 * free functions with the original names (distinct-operand forms) and
 * as in-place methods for the call sites where the C code aliases an
 * output with an input (illegal under Rust borrows).
 * -----------------------------------------------------------------*/
use crate::sundials_context::SUNContext;

#[derive(Debug, Clone, Default, PartialEq)]
pub struct NVector {
    pub data: Vec<f64>,
}

impl NVector {
    /// N_VNew_Serial: zero-initialized vector of length n.
    pub fn new(n: usize) -> Self {
        NVector { data: vec![0.0; n] }
    }

    /// N_VMake_Serial: build from existing data.
    pub fn from_slice(s: &[f64]) -> Self {
        NVector { data: s.to_vec() }
    }

    /// N_VGetLength_Serial
    #[inline]
    pub fn len(&self) -> usize {
        self.data.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// NV_Ith_S(v, i) — 0-based like the C macro.
    #[inline]
    pub fn ith(&self, i: usize) -> f64 {
        self.data[i]
    }

    #[inline]
    pub fn set_ith(&mut self, i: usize, v: f64) {
        self.data[i] = v;
    }

    /* ---------- in-place forms for aliased C call sites ---------- */

    /// z = c*z  (N_VScale(c, z, z))
    pub fn scale_inplace(&mut self, c: f64) {
        for zi in &mut self.data {
            *zi *= c;
        }
    }

    /// z = a*z + b*y  (N_VLinearSum(a, z, b, y, z))
    ///
    /// Reproduces the C serial kernel's special cases bit-for-bit: the
    /// C N_VLinearSum_Serial dispatches to VScaleSum (z = a*(x+y)) when
    /// a == b and VScaleDiff (z = a*(x-y)) when a == -b, which are NOT
    /// bitwise-equal to the generic two-multiply form. The ±1 cases
    /// (Vaxpy/VSum/VDiff/VLin1/VLin2) are bitwise-identical to the
    /// generic formula (multiplication by ±1.0 is exact) and need no
    /// special handling.
    pub fn linear_sum_with(&mut self, a: f64, b: f64, y: &NVector) {
        if a == b {
            for (zi, yi) in self.data.iter_mut().zip(&y.data) {
                *zi = a * (*zi + *yi);
            }
        } else if a == -b {
            for (zi, yi) in self.data.iter_mut().zip(&y.data) {
                *zi = a * (*zi - *yi);
            }
        } else {
            for (zi, yi) in self.data.iter_mut().zip(&y.data) {
                *zi = a * *zi + b * *yi;
            }
        }
    }

    /// z = z + b  (N_VAddConst(z, b, z))
    pub fn add_const_inplace(&mut self, b: f64) {
        for zi in &mut self.data {
            *zi += b;
        }
    }

    /// z = z .* x  (N_VProd(z, x, z))
    pub fn prod_with(&mut self, x: &NVector) {
        for (zi, xi) in self.data.iter_mut().zip(&x.data) {
            *zi *= *xi;
        }
    }

    /// z = z ./ x  (N_VDiv(z, x, z))
    pub fn div_with(&mut self, x: &NVector) {
        for (zi, xi) in self.data.iter_mut().zip(&x.data) {
            *zi /= *xi;
        }
    }

    /// z = 1 ./ z  (N_VInv(z, z))
    pub fn invert_inplace(&mut self) {
        for zi in &mut self.data {
            *zi = 1.0 / *zi;
        }
    }

    /// z = |z|  (N_VAbs(z, z))
    pub fn abs_inplace(&mut self) {
        for zi in &mut self.data {
            *zi = zi.abs();
        }
    }
}

/// N_VNew_Serial with C-style signature (context unused in serial build).
pub fn N_VNew_Serial(n: i64, _sunctx: &SUNContext) -> NVector {
    NVector::new(n as usize)
}

/// N_VClone
pub fn N_VClone(v: &NVector) -> NVector {
    NVector::new(v.len())
}

/// N_VGetArrayPointer / NV_DATA_S
pub fn N_VGetArrayPointer(v: &mut NVector) -> &mut [f64] {
    &mut v.data
}

/// N_VGetLength
pub fn N_VGetLength(v: &NVector) -> i64 {
    v.len() as i64
}

/* ------------------- standard vector operations ------------------- */

/// z = a*x + b*y
///
/// Matches the C serial kernel's arithmetic exactly: N_VLinearSum_Serial
/// dispatches to VScaleSum (z = a*(x+y)) for a == b and VScaleDiff
/// (z = a*(x-y)) for a == -b — one multiply instead of two, which is not
/// bitwise-equal to the generic form. All the ±1 special cases in C
/// (VSum/VDiff/VLin1/VLin2) are bitwise-identical to the generic formula
/// and need no separate branch.
pub fn N_VLinearSum(a: f64, x: &NVector, b: f64, y: &NVector, z: &mut NVector) {
    if a == b {
        for i in 0..z.data.len() {
            z.data[i] = a * (x.data[i] + y.data[i]);
        }
    } else if a == -b {
        for i in 0..z.data.len() {
            z.data[i] = a * (x.data[i] - y.data[i]);
        }
    } else {
        for i in 0..z.data.len() {
            z.data[i] = a * x.data[i] + b * y.data[i];
        }
    }
}

/// z = c
pub fn N_VConst(c: f64, z: &mut NVector) {
    for zi in &mut z.data {
        *zi = c;
    }
}

/// Print vector to stdout (C N_VPrint_Serial: one element per line,
/// SUN_FORMAT_E = "% .15e" — space flag pads non-negative values)
pub fn N_VPrint(x: &NVector) {
    for i in 0..x.data.len() {
        let s = crate::sundials_utils::fmt_e(x.data[i], 0, 15);
        if s.starts_with('-') {
            println!("{}", s);
        } else {
            println!(" {}", s);
        }
    }
}

/// z = x .* y
pub fn N_VProd(x: &NVector, y: &NVector, z: &mut NVector) {
    for i in 0..z.data.len() {
        z.data[i] = x.data[i] * y.data[i];
    }
}

/// z = x ./ y
pub fn N_VDiv(x: &NVector, y: &NVector, z: &mut NVector) {
    for i in 0..z.data.len() {
        z.data[i] = x.data[i] / y.data[i];
    }
}

/// z = c*x
pub fn N_VScale(c: f64, x: &NVector, z: &mut NVector) {
    for i in 0..z.data.len() {
        z.data[i] = c * x.data[i];
    }
}

/// z = |x|
pub fn N_VAbs(x: &NVector, z: &mut NVector) {
    for i in 0..z.data.len() {
        z.data[i] = x.data[i].abs();
    }
}

/// z = 1 ./ x
pub fn N_VInv(x: &NVector, z: &mut NVector) {
    for i in 0..z.data.len() {
        z.data[i] = 1.0 / x.data[i];
    }
}

/// z = x + b
pub fn N_VAddConst(x: &NVector, b: f64, z: &mut NVector) {
    for i in 0..z.data.len() {
        z.data[i] = x.data[i] + b;
    }
}

/// dot product
pub fn N_VDotProd(x: &NVector, y: &NVector) -> f64 {
    let mut sum = 0.0;
    for i in 0..x.data.len() {
        sum += x.data[i] * y.data[i];
    }
    sum
}

/// max norm
pub fn N_VMaxNorm(x: &NVector) -> f64 {
    let mut m = 0.0;
    for &xi in &x.data {
        if xi.abs() > m {
            m = xi.abs();
        }
    }
    m
}

/// WRMS norm: sqrt( sum((x_i * w_i)^2) / n )
pub fn N_VWrmsNorm(x: &NVector, w: &NVector) -> f64 {
    let n = x.data.len();
    if n == 0 {
        return 0.0;
    }
    let mut sum = 0.0;
    for i in 0..n {
        let p = x.data[i] * w.data[i];
        sum += p * p;
    }
    (sum / n as f64).sqrt()
}

/// masked WRMS norm (only entries with id_i > 0 contribute)
pub fn N_VWrmsNormMask(x: &NVector, w: &NVector, id: &NVector) -> f64 {
    let n = x.data.len();
    if n == 0 {
        return 0.0;
    }
    let mut sum = 0.0;
    for i in 0..n {
        if id.data[i] > 0.0 {
            let p = x.data[i] * w.data[i];
            sum += p * p;
        }
    }
    (sum / n as f64).sqrt()
}

/// minimum entry
pub fn N_VMin(x: &NVector) -> f64 {
    if x.data.is_empty() {
        return 0.0;
    }
    let mut m = x.data[0];
    for &xi in &x.data[1..] {
        if xi < m {
            m = xi;
        }
    }
    m
}

/// weighted Euclidean l2 norm: sqrt( sum((x_i * w_i)^2) )
pub fn N_VWL2Norm(x: &NVector, w: &NVector) -> f64 {
    let mut sum = 0.0;
    for i in 0..x.data.len() {
        let p = x.data[i] * w.data[i];
        sum += p * p;
    }
    sum.sqrt()
}

/// l1 norm
pub fn N_VL1Norm(x: &NVector) -> f64 {
    let mut sum = 0.0;
    for &xi in &x.data {
        sum += xi.abs();
    }
    sum
}

/// z_i = 1 if |x_i| >= c, else 0
pub fn N_VCompare(c: f64, x: &NVector, z: &mut NVector) {
    for i in 0..z.data.len() {
        z.data[i] = if x.data[i].abs() >= c { 1.0 } else { 0.0 };
    }
}

/// z_i = 1/x_i with test; returns false if any x_i == 0
pub fn N_VInvTest(x: &NVector, z: &mut NVector) -> bool {
    let mut no_zero_found = true;
    for i in 0..z.data.len() {
        if x.data[i] == 0.0 {
            no_zero_found = false;
        } else {
            z.data[i] = 1.0 / x.data[i];
        }
    }
    no_zero_found
}

/// Constraint test (c_i in {-2,-1,0,1,2}); m_i = 1 where violated.
/// Returns true if all constraints satisfied.
pub fn N_VConstrMask(c: &NVector, x: &NVector, m: &mut NVector) -> bool {
    let mut temp = 0.0;
    for i in 0..x.data.len() {
        m.data[i] = 0.0;
        if c.data[i] == 0.0 {
            continue;
        }
        if c.data[i] > 1.5 || c.data[i] < -1.5 {
            // c = +-2 : x must be strictly positive/negative
            if x.data[i] * c.data[i] <= 0.0 {
                temp = 1.0;
                m.data[i] = 1.0;
            }
            continue;
        }
        // c = +-1 : x must be non-negative/non-positive
        if x.data[i] * c.data[i] < 0.0 {
            temp = 1.0;
            m.data[i] = 1.0;
        }
    }
    temp != 1.0
}

/// min over i with denom_i != 0 of num_i/denom_i (BIG_REAL if none)
pub fn N_VMinQuotient(num: &NVector, denom: &NVector) -> f64 {
    let mut min = f64::MAX;
    let mut notevenonce = true;
    for i in 0..num.data.len() {
        if denom.data[i] == 0.0 {
            continue;
        }
        let q = num.data[i] / denom.data[i];
        if notevenonce {
            min = q;
            notevenonce = false;
        } else if q < min {
            min = q;
        }
    }
    min
}

/// N_VSpace_Serial
pub fn N_VSpace(v: &NVector) -> (i64, i64) {
    (v.len() as i64, 1)
}

/// N_VLinearCombination_Serial for a destination DISTINCT from every
/// operand (the C kernel's general branch; the z==X[0] in-place
/// branches are separate in-place calls under workspace rule 5).
/// Per element the accumulation order is c[0]*x0, += c[1]*x1, ... —
/// bit-identical to the C special cases for nvec = 1 and 2 as well.
pub fn N_VLinearCombination(
    nvec: i32,
    c: &[f64],
    x: &[&NVector],
    z: &mut NVector,
) -> i32 {
    /* invalid number of vectors */
    if nvec < 1 {
        return -1;
    }
    let nvec = nvec as usize;
    for k in 0..z.data.len() {
        let mut acc = c[0] * x[0].data[k];
        for i in 1..nvec {
            acc += c[i] * x[i].data[k];
        }
        z.data[k] = acc;
    }
    0
}

/// N_VBufSize_Serial: buffer size in bytes (length * sizeof(sunrealtype)).
pub fn N_VBufSize(x: &NVector, size: &mut i64) -> crate::sundials_errors::SUNErrCode {
    *size = (x.data.len() as i64) * (std::mem::size_of::<f64>() as i64);
    crate::sundials_errors::SUN_SUCCESS
}

/// N_VBufPack_Serial: copy the vector data into a byte buffer (C copies
/// sunrealtype-by-sunrealtype through a void* buffer; native-endian bytes
/// here so Pack/Unpack round-trip bit-exactly).
pub fn N_VBufPack(x: &NVector, buf: &mut [u8]) -> crate::sundials_errors::SUNErrCode {
    let n = x.data.len();
    for i in 0..n {
        buf[8 * i..8 * i + 8].copy_from_slice(&x.data[i].to_ne_bytes());
    }
    crate::sundials_errors::SUN_SUCCESS
}

/// N_VBufUnpack_Serial: copy a byte buffer back into the vector data.
pub fn N_VBufUnpack(x: &mut NVector, buf: &[u8]) -> crate::sundials_errors::SUNErrCode {
    let n = x.data.len();
    for i in 0..n {
        let mut b = [0u8; 8];
        b.copy_from_slice(&buf[8 * i..8 * i + 8]);
        x.data[i] = f64::from_ne_bytes(b);
    }
    crate::sundials_errors::SUN_SUCCESS
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wrmsnorm() {
        let x = NVector::from_slice(&[3.0, 4.0]);
        let w = NVector::from_slice(&[1.0, 1.0]);
        let n = N_VWrmsNorm(&x, &w);
        assert!((n - (12.5f64).sqrt()).abs() < 1e-15);
    }

    #[test]
    fn constrmask() {
        let c = NVector::from_slice(&[2.0, 1.0, 0.0, -1.0, -2.0]);
        let x = NVector::from_slice(&[1.0, 0.0, -5.0, 0.0, -1.0]);
        let mut m = NVector::new(5);
        assert!(N_VConstrMask(&c, &x, &mut m));
        let x2 = NVector::from_slice(&[0.0, 0.0, -5.0, 0.0, -1.0]);
        assert!(!N_VConstrMask(&c, &x2, &mut m));
        assert_eq!(m.data[0], 1.0);
    }
}
