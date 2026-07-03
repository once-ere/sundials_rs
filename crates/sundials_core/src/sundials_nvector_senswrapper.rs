/* -----------------------------------------------------------------
 * Translated from src/sundials/sundials_nvector_senswrapper.c and
 * include/sundials/sundials_nvector_senswrapper.h (SUNDIALS 7.7.0).
 *
 * The C SensWrapper is a "fake" N_Vector whose content is an array of
 * `nvecs` real N_Vectors; every attached op loops the corresponding
 * sub-vector op over that array. It is used by the CVODES/IDAS
 * sensitivity nonlinear solvers so a Newton iteration can treat the
 * whole set of sensitivity vectors as a single vector.
 *
 * Rust design (mirrors nvector_serial.rs): a plain struct holding an
 * owned `Vec<NVector>` plus the C `own_vecs` flag, with the C function
 * names preserved as free functions and in-place method variants for
 * the call sites where the C code aliases an output with an input.
 *
 * Ownership: in C, `own_vecs == SUNFALSE` marks a wrapper whose slots
 * borrow vectors owned by someone else (the empty-constructor path,
 * where CVODES stores raw pointers to its own sensitivity vectors).
 * Safe Rust cannot hold those borrows in an owning struct, so — as
 * with the SUNMemoryHelper_Alias port — borrow semantics collapse to
 * ownership: the `vecs` are always owned by the wrapper and callers
 * clone data in (clone-on-wrap) and copy results back out. The
 * `own_vecs` flag is kept for structural fidelity and still tracks
 * which C constructor path produced the wrapper (false for the
 * *Empty* constructors, true once N_VNew/N_VClone populate the slots),
 * but it no longer changes destruction behavior (RAII frees
 * everything).
 *
 * Reduction semantics (verified against the C loops — these matter
 * for the CVODES sensitivity convergence tests):
 *   - N_VDotProd_SensWrapper       SUMS the per-sub-vector dot products
 *   - N_VMaxNorm_SensWrapper       MAX over sub-vectors (init 0)
 *   - N_VWrmsNorm_SensWrapper      MAX of per-sub-vector WRMS norms
 *   - N_VWrmsNormMask_SensWrapper  MAX of per-sub-vector masked norms
 *   - N_VWL2Norm_SensWrapper       MAX of per-sub-vector WL2 norms
 *                                  (NOT the root-sum-of-squares)
 *   - N_VL1Norm_SensWrapper        MAX of per-sub-vector L1 norms
 *                                  (NOT the sum)
 *   - N_VMin_SensWrapper           MIN over sub-vectors, initialized
 *                                  from sub-vector 0
 *   - N_VMinQuotient_SensWrapper   MIN over sub-vectors, initialized
 *                                  from sub-vector 0
 *   - N_VInvTest_SensWrapper       logical AND over sub-vectors, all
 *                                  sub-vectors always processed
 *   - N_VConstrMask_SensWrapper    logical AND over sub-vectors, all
 *                                  sub-vectors always processed; the
 *                                  constraint vector `c` is passed to
 *                                  every sub-vector op WHOLE (it is a
 *                                  plain vector, not a wrapper)
 * -----------------------------------------------------------------*/
use crate::nvector_serial::{
    N_VAbs, N_VAddConst, N_VClone, N_VCompare, N_VConst, N_VConstrMask, N_VDiv, N_VDotProd,
    N_VInv, N_VInvTest, N_VL1Norm, N_VLinearSum, N_VMaxNorm, N_VMin, N_VMinQuotient, N_VProd,
    N_VScale, N_VWL2Norm, N_VWrmsNorm, N_VWrmsNormMask, NVector,
};

/// struct _N_VectorContent_SensWrapper (vector object and content are
/// merged into one struct, as everywhere else in this port).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct NVectorSensWrapper {
    /// NV_VECS_SW: the wrapped sub-vectors. In C this is an array of
    /// `nvecs` (possibly NULL) N_Vector pointers; here the slots are
    /// always present and empty slots are zero-length NVectors.
    pub vecs: Vec<NVector>,
    /// NV_OWN_VECS_SW: kept for structural fidelity; see module header
    /// (in safe Rust the wrapper always owns its sub-vectors).
    pub own_vecs: bool,
}

impl NVectorSensWrapper {
    /// NV_NVECS_SW(v)
    #[inline]
    pub fn nvecs(&self) -> usize {
        self.vecs.len()
    }

    /// NV_VEC_SW(v, i)
    #[inline]
    pub fn vec(&self, i: usize) -> &NVector {
        &self.vecs[i]
    }

    /// NV_VEC_SW(v, i) as assignable lvalue.
    #[inline]
    pub fn vec_mut(&mut self, i: usize) -> &mut NVector {
        &mut self.vecs[i]
    }

    /* ---------- in-place forms for aliased C call sites ----------
     * These mirror the nvector_serial.rs convention: when a C call
     * site aliases an operand with the output (e.g.
     * N_VLinearSum(a, x, b, y, y)), Rust borrows forbid the free
     * function, so the aliased operand becomes `self`. Each method
     * delegates to the sub-vector in-place method, so the per-element
     * arithmetic (including the a == ±b single-multiply kernels) is
     * bit-identical to the free-function path. */

    /// z = c*z  (N_VScale_SensWrapper(c, z, z))
    pub fn scale_inplace(&mut self, c: f64) {
        for i in 0..self.nvecs() {
            self.vecs[i].scale_inplace(c);
        }
    }

    /// z = a*z + b*y  (N_VLinearSum_SensWrapper(a, z, b, y, z))
    pub fn linear_sum_with(&mut self, a: f64, b: f64, y: &NVectorSensWrapper) {
        for i in 0..self.nvecs() {
            self.vecs[i].linear_sum_with(a, b, y.vec(i));
        }
    }

    /// z = z + b  (N_VAddConst_SensWrapper(z, b, z))
    pub fn add_const_inplace(&mut self, b: f64) {
        for i in 0..self.nvecs() {
            self.vecs[i].add_const_inplace(b);
        }
    }

    /// z = z .* x  (N_VProd_SensWrapper(z, x, z))
    pub fn prod_with(&mut self, x: &NVectorSensWrapper) {
        for i in 0..self.nvecs() {
            self.vecs[i].prod_with(x.vec(i));
        }
    }

    /// z = z ./ x  (N_VDiv_SensWrapper(z, x, z))
    pub fn div_with(&mut self, x: &NVectorSensWrapper) {
        for i in 0..self.nvecs() {
            self.vecs[i].div_with(x.vec(i));
        }
    }

    /// z = 1 ./ z  (N_VInv_SensWrapper(z, z))
    pub fn invert_inplace(&mut self) {
        for i in 0..self.nvecs() {
            self.vecs[i].invert_inplace();
        }
    }

    /// z = |z|  (N_VAbs_SensWrapper(z, z))
    pub fn abs_inplace(&mut self) {
        for i in 0..self.nvecs() {
            self.vecs[i].abs_inplace();
        }
    }
}

/*==============================================================================
  Constructors
  ============================================================================*/

/// N_VNewEmpty_SensWrapper: create a new empty vector wrapper with space
/// for `nvecs` vectors. Returns None if `nvecs < 1` (C returns NULL).
/// The C slots start as NULL pointers with own_vecs = SUNFALSE; here they
/// are zero-length NVectors the caller fills via `vec_mut` (clone-on-wrap
/// replaces the C borrow of caller-owned vectors). The SUNContext
/// argument is dropped (serial build; see sundials_context.rs).
pub fn N_VNewEmpty_SensWrapper(nvecs: i32) -> Option<NVectorSensWrapper> {
    /* return if wrapper is empty */
    if nvecs < 1 {
        return None;
    }

    Some(NVectorSensWrapper {
        vecs: vec![NVector::default(); nvecs as usize],
        own_vecs: false,
    })
}

/// N_VNew_SensWrapper: create a wrapper holding `count` clones of the
/// template vector `w` (N_VClone semantics: same length, zero data).
/// Returns None if `count < 1` (C returns NULL).
pub fn N_VNew_SensWrapper(count: i32, w: &NVector) -> Option<NVectorSensWrapper> {
    let mut v = N_VNewEmpty_SensWrapper(count)?;

    for i in 0..v.nvecs() {
        *v.vec_mut(i) = N_VClone(w);
    }

    /* update own vectors status */
    v.own_vecs = true;

    Some(v)
}

/*==============================================================================
  Clone operations
  ============================================================================*/

/// N_VCloneEmpty_SensWrapper: create an empty clone of the wrapper `w`
/// (same nvecs, slots unallocated, own_vecs = SUNFALSE). Returns None if
/// `w` has fewer than 1 sub-vector (C returns NULL). The ops-table copy
/// in the C code has no Rust counterpart (dispatch is static).
pub fn N_VCloneEmpty_SensWrapper(w: &NVectorSensWrapper) -> Option<NVectorSensWrapper> {
    if w.nvecs() < 1 {
        return None;
    }

    Some(NVectorSensWrapper {
        vecs: vec![NVector::default(); w.nvecs()],
        own_vecs: false,
    })
}

/// N_VClone_SensWrapper: create a clone of the wrapper `w`. Like the C
/// (which calls N_VClone on every sub-vector), the clone has sub-vectors
/// of matching lengths but does NOT copy their data — use Rust's derived
/// `.clone()` for a deep copy. Returns None if `w` is empty.
pub fn N_VClone_SensWrapper(w: &NVectorSensWrapper) -> Option<NVectorSensWrapper> {
    /* create empty wrapper */
    let mut v = N_VCloneEmpty_SensWrapper(w)?;

    /* update own vectors status */
    v.own_vecs = true;

    /* allocate arrays */
    for i in 0..v.nvecs() {
        *v.vec_mut(i) = N_VClone(w.vec(i));
    }

    Some(v)
}

/*==============================================================================
  Destructor
  ============================================================================*/

/// N_VDestroy_SensWrapper — RAII frees the sub-vectors regardless of the
/// own_vecs flag (in C, borrowed slots were skipped; borrows collapse to
/// ownership here). Kept for API parity.
pub fn N_VDestroy_SensWrapper(v: NVectorSensWrapper) {
    drop(v);
}

/*==============================================================================
  Standard vector operations
  ============================================================================*/

/// z = a*x + b*y, sub-vector-wise. Aliased form: `linear_sum_with`.
pub fn N_VLinearSum_SensWrapper(
    a: f64,
    x: &NVectorSensWrapper,
    b: f64,
    y: &NVectorSensWrapper,
    z: &mut NVectorSensWrapper,
) {
    for i in 0..x.nvecs() {
        N_VLinearSum(a, x.vec(i), b, y.vec(i), z.vec_mut(i));
    }
}

/// z = c, sub-vector-wise.
pub fn N_VConst_SensWrapper(c: f64, z: &mut NVectorSensWrapper) {
    for i in 0..z.nvecs() {
        N_VConst(c, z.vec_mut(i));
    }
}

/// z = x .* y, sub-vector-wise. Aliased form: `prod_with`.
pub fn N_VProd_SensWrapper(x: &NVectorSensWrapper, y: &NVectorSensWrapper, z: &mut NVectorSensWrapper) {
    for i in 0..x.nvecs() {
        N_VProd(x.vec(i), y.vec(i), z.vec_mut(i));
    }
}

/// z = x ./ y, sub-vector-wise. Aliased form: `div_with`.
pub fn N_VDiv_SensWrapper(x: &NVectorSensWrapper, y: &NVectorSensWrapper, z: &mut NVectorSensWrapper) {
    for i in 0..x.nvecs() {
        N_VDiv(x.vec(i), y.vec(i), z.vec_mut(i));
    }
}

/// z = c*x, sub-vector-wise. Aliased form: `scale_inplace`.
pub fn N_VScale_SensWrapper(c: f64, x: &NVectorSensWrapper, z: &mut NVectorSensWrapper) {
    for i in 0..x.nvecs() {
        N_VScale(c, x.vec(i), z.vec_mut(i));
    }
}

/// z = |x|, sub-vector-wise. Aliased form: `abs_inplace`.
pub fn N_VAbs_SensWrapper(x: &NVectorSensWrapper, z: &mut NVectorSensWrapper) {
    for i in 0..x.nvecs() {
        N_VAbs(x.vec(i), z.vec_mut(i));
    }
}

/// z = 1 ./ x, sub-vector-wise. Aliased form: `invert_inplace`.
pub fn N_VInv_SensWrapper(x: &NVectorSensWrapper, z: &mut NVectorSensWrapper) {
    for i in 0..x.nvecs() {
        N_VInv(x.vec(i), z.vec_mut(i));
    }
}

/// z = x + b, sub-vector-wise. Aliased form: `add_const_inplace`.
pub fn N_VAddConst_SensWrapper(x: &NVectorSensWrapper, b: f64, z: &mut NVectorSensWrapper) {
    for i in 0..x.nvecs() {
        N_VAddConst(x.vec(i), b, z.vec_mut(i));
    }
}

/// Dot product: SUM of the per-sub-vector dot products.
pub fn N_VDotProd_SensWrapper(x: &NVectorSensWrapper, y: &NVectorSensWrapper) -> f64 {
    let mut sum = 0.0;

    for i in 0..x.nvecs() {
        sum += N_VDotProd(x.vec(i), y.vec(i));
    }

    sum
}

/// Max norm: MAX of the per-sub-vector max norms (init 0).
pub fn N_VMaxNorm_SensWrapper(x: &NVectorSensWrapper) -> f64 {
    let mut max = 0.0;

    for i in 0..x.nvecs() {
        let tmp = N_VMaxNorm(x.vec(i));
        if tmp > max {
            max = tmp;
        }
    }

    max
}

/// WRMS norm: MAX of the per-sub-vector WRMS norms (init 0) — the
/// sub-vector norms are NOT combined into one RMS.
pub fn N_VWrmsNorm_SensWrapper(x: &NVectorSensWrapper, w: &NVectorSensWrapper) -> f64 {
    let mut nrm = 0.0;

    for i in 0..x.nvecs() {
        let tmp = N_VWrmsNorm(x.vec(i), w.vec(i));
        if tmp > nrm {
            nrm = tmp;
        }
    }

    nrm
}

/// Masked WRMS norm: MAX of the per-sub-vector masked WRMS norms (init 0).
pub fn N_VWrmsNormMask_SensWrapper(
    x: &NVectorSensWrapper,
    w: &NVectorSensWrapper,
    id: &NVectorSensWrapper,
) -> f64 {
    let mut nrm = 0.0;

    for i in 0..x.nvecs() {
        let tmp = N_VWrmsNormMask(x.vec(i), w.vec(i), id.vec(i));
        if tmp > nrm {
            nrm = tmp;
        }
    }

    nrm
}

/// Min: MIN of the per-sub-vector minima, initialized from sub-vector 0.
pub fn N_VMin_SensWrapper(x: &NVectorSensWrapper) -> f64 {
    let mut min = N_VMin(x.vec(0));

    for i in 1..x.nvecs() {
        let tmp = N_VMin(x.vec(i));
        if tmp < min {
            min = tmp;
        }
    }

    min
}

/// Weighted L2 norm: MAX of the per-sub-vector WL2 norms (init 0) — NOT
/// the root-sum-of-squares over all sub-vectors.
pub fn N_VWL2Norm_SensWrapper(x: &NVectorSensWrapper, w: &NVectorSensWrapper) -> f64 {
    let mut nrm = 0.0;

    for i in 0..x.nvecs() {
        let tmp = N_VWL2Norm(x.vec(i), w.vec(i));
        if tmp > nrm {
            nrm = tmp;
        }
    }

    nrm
}

/// L1 norm: MAX of the per-sub-vector L1 norms (init 0) — NOT the sum.
pub fn N_VL1Norm_SensWrapper(x: &NVectorSensWrapper) -> f64 {
    let mut nrm = 0.0;

    for i in 0..x.nvecs() {
        let tmp = N_VL1Norm(x.vec(i));
        if tmp > nrm {
            nrm = tmp;
        }
    }

    nrm
}

/// z_i = 1 if |x_i| >= c else 0, sub-vector-wise.
pub fn N_VCompare_SensWrapper(c: f64, x: &NVectorSensWrapper, z: &mut NVectorSensWrapper) {
    for i in 0..x.nvecs() {
        N_VCompare(c, x.vec(i), z.vec_mut(i));
    }
}

/// Inversion with zero test: logical AND of the per-sub-vector results;
/// every sub-vector is processed (no short-circuit) so all valid entries
/// of z are written even when a zero is found.
pub fn N_VInvTest_SensWrapper(x: &NVectorSensWrapper, z: &mut NVectorSensWrapper) -> bool {
    let mut no_zero_found = true;

    for i in 0..x.nvecs() {
        let tmp = N_VInvTest(x.vec(i), z.vec_mut(i));
        if !tmp {
            no_zero_found = false;
        }
    }

    no_zero_found
}

/// Constraint mask: logical AND of the per-sub-vector tests; every
/// sub-vector is processed (no short-circuit). NOTE: exactly as in the C,
/// the constraint vector `c` is a PLAIN vector applied whole against each
/// sub-vector of x (the C passes `c` un-unwrapped to N_VConstrMask).
pub fn N_VConstrMask_SensWrapper(
    c: &NVector,
    x: &NVectorSensWrapper,
    m: &mut NVectorSensWrapper,
) -> bool {
    let mut test = true;

    for i in 0..x.nvecs() {
        let tmp = N_VConstrMask(c, x.vec(i), m.vec_mut(i));
        if !tmp {
            test = false;
        }
    }

    test
}

/// Min quotient: MIN of the per-sub-vector min quotients, initialized
/// from sub-vector 0.
pub fn N_VMinQuotient_SensWrapper(num: &NVectorSensWrapper, denom: &NVectorSensWrapper) -> f64 {
    let mut min = N_VMinQuotient(num.vec(0), denom.vec(0));

    for i in 1..num.nvecs() {
        let tmp = N_VMinQuotient(num.vec(i), denom.vec(i));
        if tmp < min {
            min = tmp;
        }
    }

    min
}

#[cfg(test)]
mod tests {
    use super::*;

    /// x = [[1, 2], [3, -4]]
    fn make_x() -> NVectorSensWrapper {
        NVectorSensWrapper {
            vecs: vec![
                NVector::from_slice(&[1.0, 2.0]),
                NVector::from_slice(&[3.0, -4.0]),
            ],
            own_vecs: true,
        }
    }

    /// y = [[5, 6], [7, 8]]
    fn make_y() -> NVectorSensWrapper {
        NVectorSensWrapper {
            vecs: vec![
                NVector::from_slice(&[5.0, 6.0]),
                NVector::from_slice(&[7.0, 8.0]),
            ],
            own_vecs: true,
        }
    }

    /// all-ones weight wrapper
    fn make_ones() -> NVectorSensWrapper {
        NVectorSensWrapper {
            vecs: vec![
                NVector::from_slice(&[1.0, 1.0]),
                NVector::from_slice(&[1.0, 1.0]),
            ],
            own_vecs: true,
        }
    }

    fn make_z() -> NVectorSensWrapper {
        N_VNew_SensWrapper(2, &NVector::new(2)).unwrap()
    }

    fn flat(v: &NVectorSensWrapper) -> Vec<f64> {
        v.vecs.iter().flat_map(|s| s.data.clone()).collect()
    }

    /* ------------------------- construction ------------------------- */

    #[test]
    fn new_empty() {
        let v = N_VNewEmpty_SensWrapper(3).unwrap();
        assert_eq!(v.nvecs(), 3);
        assert!(!v.own_vecs);
        for i in 0..3 {
            assert!(v.vec(i).is_empty());
        }
        /* C returns NULL for nvecs < 1 */
        assert!(N_VNewEmpty_SensWrapper(0).is_none());
        assert!(N_VNewEmpty_SensWrapper(-1).is_none());
    }

    #[test]
    fn new_from_template() {
        let w = NVector::from_slice(&[1.0, 2.0, 3.0]);
        let v = N_VNew_SensWrapper(2, &w).unwrap();
        assert_eq!(v.nvecs(), 2);
        assert!(v.own_vecs);
        /* N_VClone allocates same length, zero data (no copy) */
        for i in 0..2 {
            assert_eq!(v.vec(i).len(), 3);
            assert_eq!(v.vec(i).data, vec![0.0; 3]);
        }
        assert!(N_VNew_SensWrapper(0, &w).is_none());
    }

    #[test]
    fn clone_ops() {
        let x = make_x();
        let e = N_VCloneEmpty_SensWrapper(&x).unwrap();
        assert_eq!(e.nvecs(), 2);
        assert!(!e.own_vecs);
        assert!(e.vec(0).is_empty() && e.vec(1).is_empty());

        let c = N_VClone_SensWrapper(&x).unwrap();
        assert_eq!(c.nvecs(), 2);
        assert!(c.own_vecs);
        /* clone matches lengths but does not copy data (C N_VClone) */
        assert_eq!(c.vec(0).data, vec![0.0, 0.0]);
        assert_eq!(c.vec(1).data, vec![0.0, 0.0]);

        let empty = NVectorSensWrapper::default();
        assert!(N_VCloneEmpty_SensWrapper(&empty).is_none());
        assert!(N_VClone_SensWrapper(&empty).is_none());

        N_VDestroy_SensWrapper(c);
    }

    /* ---------------------- arithmetic (2 sub-vectors) ---------------------- */

    #[test]
    fn linear_sum() {
        /* z = 2*x + 3*y = [[2+15, 4+18], [6+21, -8+24]] */
        let (x, y, mut z) = (make_x(), make_y(), make_z());
        N_VLinearSum_SensWrapper(2.0, &x, 3.0, &y, &mut z);
        assert_eq!(flat(&z), vec![17.0, 22.0, 27.0, 16.0]);

        /* aliased form: y := 2*y + 3*x */
        let mut ya = make_y();
        ya.linear_sum_with(2.0, 3.0, &x);
        assert_eq!(flat(&ya), vec![13.0, 18.0, 23.0, 4.0]);
    }

    #[test]
    fn const_op() {
        let mut z = make_z();
        N_VConst_SensWrapper(4.5, &mut z);
        assert_eq!(flat(&z), vec![4.5, 4.5, 4.5, 4.5]);
    }

    #[test]
    fn prod() {
        let (x, y, mut z) = (make_x(), make_y(), make_z());
        N_VProd_SensWrapper(&x, &y, &mut z);
        assert_eq!(flat(&z), vec![5.0, 12.0, 21.0, -32.0]);

        let mut xa = make_x();
        xa.prod_with(&y);
        assert_eq!(flat(&xa), vec![5.0, 12.0, 21.0, -32.0]);
    }

    #[test]
    fn div() {
        let (x, y, mut z) = (make_x(), make_y(), make_z());
        N_VDiv_SensWrapper(&x, &y, &mut z);
        assert_eq!(flat(&z), vec![1.0 / 5.0, 2.0 / 6.0, 3.0 / 7.0, -0.5]);

        let mut xa = make_x();
        xa.div_with(&y);
        assert_eq!(flat(&xa), vec![1.0 / 5.0, 2.0 / 6.0, 3.0 / 7.0, -0.5]);
    }

    #[test]
    fn scale() {
        let (x, mut z) = (make_x(), make_z());
        N_VScale_SensWrapper(2.0, &x, &mut z);
        assert_eq!(flat(&z), vec![2.0, 4.0, 6.0, -8.0]);

        let mut xa = make_x();
        xa.scale_inplace(2.0);
        assert_eq!(flat(&xa), vec![2.0, 4.0, 6.0, -8.0]);
    }

    #[test]
    fn abs() {
        let (x, mut z) = (make_x(), make_z());
        N_VAbs_SensWrapper(&x, &mut z);
        assert_eq!(flat(&z), vec![1.0, 2.0, 3.0, 4.0]);

        let mut xa = make_x();
        xa.abs_inplace();
        assert_eq!(flat(&xa), vec![1.0, 2.0, 3.0, 4.0]);
    }

    #[test]
    fn inv() {
        let (x, mut z) = (make_x(), make_z());
        N_VInv_SensWrapper(&x, &mut z);
        assert_eq!(flat(&z), vec![1.0, 0.5, 1.0 / 3.0, -0.25]);

        let mut xa = make_x();
        xa.invert_inplace();
        assert_eq!(flat(&xa), vec![1.0, 0.5, 1.0 / 3.0, -0.25]);
    }

    #[test]
    fn add_const() {
        let (x, mut z) = (make_x(), make_z());
        N_VAddConst_SensWrapper(&x, 1.0, &mut z);
        assert_eq!(flat(&z), vec![2.0, 3.0, 4.0, -3.0]);

        let mut xa = make_x();
        xa.add_const_inplace(1.0);
        assert_eq!(flat(&xa), vec![2.0, 3.0, 4.0, -3.0]);
    }

    #[test]
    fn compare() {
        /* z_i = |x_i| >= 2 → [[0,1],[1,1]] */
        let (x, mut z) = (make_x(), make_z());
        N_VCompare_SensWrapper(2.0, &x, &mut z);
        assert_eq!(flat(&z), vec![0.0, 1.0, 1.0, 1.0]);
    }

    /* -------------------- reductions (cross-sub-vector) --------------------
     * Expected values written from the C loops in
     * sundials_nvector_senswrapper.c, not from this implementation. */

    #[test]
    fn dotprod_sums_across_subvectors() {
        /* C: sum += N_VDotProd(x_i, y_i)
         * sub0: 1*5 + 2*6 = 17;  sub1: 3*7 + (-4)*8 = -11;  sum = 6 */
        let (x, y) = (make_x(), make_y());
        assert_eq!(N_VDotProd_SensWrapper(&x, &y), 6.0);
    }

    #[test]
    fn maxnorm_max_across_subvectors() {
        /* C: max of per-sub-vector max norms; sub0: 2, sub1: 4 → 4 */
        let x = make_x();
        assert_eq!(N_VMaxNorm_SensWrapper(&x), 4.0);
    }

    #[test]
    fn wrmsnorm_max_across_subvectors() {
        /* C: max of per-sub-vector WRMS norms, NOT a combined RMS.
         * sub0: sqrt((1+4)/2) = sqrt(2.5); sub1: sqrt((9+16)/2) = sqrt(12.5)
         * → sqrt(12.5) */
        let (x, w) = (make_x(), make_ones());
        assert_eq!(N_VWrmsNorm_SensWrapper(&x, &w), (12.5f64).sqrt());
    }

    #[test]
    fn wrmsnormmask_max_across_subvectors() {
        /* id = [[1,0],[0,1]]: sub0 uses x=1 → sqrt(1/2);
         * sub1 uses x=-4 → sqrt(16/2) = sqrt(8) → max sqrt(8) */
        let (x, w) = (make_x(), make_ones());
        let id = NVectorSensWrapper {
            vecs: vec![
                NVector::from_slice(&[1.0, 0.0]),
                NVector::from_slice(&[0.0, 1.0]),
            ],
            own_vecs: true,
        };
        assert_eq!(N_VWrmsNormMask_SensWrapper(&x, &w, &id), (8.0f64).sqrt());
    }

    #[test]
    fn min_across_subvectors() {
        /* C: min initialized from sub-vector 0 (min 1), then sub1 min -4 → -4 */
        let x = make_x();
        assert_eq!(N_VMin_SensWrapper(&x), -4.0);
        /* single sub-vector: init path only */
        let x0 = NVectorSensWrapper {
            vecs: vec![NVector::from_slice(&[7.0, 5.0])],
            own_vecs: true,
        };
        assert_eq!(N_VMin_SensWrapper(&x0), 5.0);
    }

    #[test]
    fn wl2norm_max_across_subvectors() {
        /* C: MAX of per-sub-vector WL2 norms, not root-sum-square.
         * sub0: sqrt(1+4) = sqrt(5); sub1: sqrt(9+16) = 5 → 5 */
        let (x, w) = (make_x(), make_ones());
        assert_eq!(N_VWL2Norm_SensWrapper(&x, &w), 5.0);
    }

    #[test]
    fn l1norm_max_across_subvectors() {
        /* C: MAX of per-sub-vector L1 norms, NOT the sum.
         * sub0: 1+2 = 3; sub1: 3+4 = 7 → 7 (a summing port would give 10) */
        let x = make_x();
        assert_eq!(N_VL1Norm_SensWrapper(&x), 7.0);
    }

    #[test]
    fn invtest_and_across_subvectors() {
        /* x = [[1,2],[3,0]]: sub0 passes, sub1 fails → overall false;
         * both sub-vectors still processed, zero entry left untouched */
        let x = NVectorSensWrapper {
            vecs: vec![
                NVector::from_slice(&[1.0, 2.0]),
                NVector::from_slice(&[3.0, 0.0]),
            ],
            own_vecs: true,
        };
        let mut z = NVectorSensWrapper {
            vecs: vec![
                NVector::from_slice(&[99.0, 99.0]),
                NVector::from_slice(&[99.0, 99.0]),
            ],
            own_vecs: true,
        };
        assert!(!N_VInvTest_SensWrapper(&x, &mut z));
        assert_eq!(z.vec(0).data, vec![1.0, 0.5]);
        assert_eq!(z.vec(1).data[0], 1.0 / 3.0);
        assert_eq!(z.vec(1).data[1], 99.0); /* untouched, as in serial C */

        /* all nonzero → true */
        let ok = make_x();
        let mut z2 = make_z();
        assert!(N_VInvTest_SensWrapper(&ok, &mut z2));
    }

    #[test]
    fn constrmask_plain_c_and_across_subvectors() {
        /* C passes the constraint vector whole to each sub-vector:
         * c = [2, -1]; x = [[1,-1],[1,1]]
         * sub0: 1 > 0 ok, -1 <= 0 ok → pass
         * sub1: 1 > 0 ok,  1 <= 0 violated → fail → overall false */
        let c = NVector::from_slice(&[2.0, -1.0]);
        let x = NVectorSensWrapper {
            vecs: vec![
                NVector::from_slice(&[1.0, -1.0]),
                NVector::from_slice(&[1.0, 1.0]),
            ],
            own_vecs: true,
        };
        let mut m = make_z();
        assert!(!N_VConstrMask_SensWrapper(&c, &x, &mut m));
        assert_eq!(m.vec(0).data, vec![0.0, 0.0]);
        assert_eq!(m.vec(1).data, vec![0.0, 1.0]);

        /* both satisfy → true */
        let x_ok = NVectorSensWrapper {
            vecs: vec![
                NVector::from_slice(&[1.0, -1.0]),
                NVector::from_slice(&[2.0, 0.0]),
            ],
            own_vecs: true,
        };
        let mut m2 = make_z();
        assert!(N_VConstrMask_SensWrapper(&c, &x_ok, &mut m2));
    }

    #[test]
    fn minquotient_min_across_subvectors() {
        /* C: init from sub-vector 0, then min.
         * num = [[1,2],[6,3]], denom = [[2,0],[3,1]]
         * sub0: only 1/2 (denom 0 skipped) → 0.5
         * sub1: min(6/3, 3/1) = 2 → overall 0.5 */
        let num = NVectorSensWrapper {
            vecs: vec![
                NVector::from_slice(&[1.0, 2.0]),
                NVector::from_slice(&[6.0, 3.0]),
            ],
            own_vecs: true,
        };
        let denom = NVectorSensWrapper {
            vecs: vec![
                NVector::from_slice(&[2.0, 0.0]),
                NVector::from_slice(&[3.0, 1.0]),
            ],
            own_vecs: true,
        };
        assert_eq!(N_VMinQuotient_SensWrapper(&num, &denom), 0.5);
    }
}
