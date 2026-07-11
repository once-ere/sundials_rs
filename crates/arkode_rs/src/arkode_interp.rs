/* -----------------------------------------------------------------
 * Translation of sundials-7.7.0/src/arkode/arkode_interp.c.
 *
 * Section I generic dispatchers: C passes (ark_mem, interp) where
 * interp is ALWAYS ark_mem->interp; the Rust dispatchers take only
 * ark_mem and use take/put-back on ark_mem.interp so the
 * implementations can borrow ark_mem and the interp content
 * disjointly (pinned convention, see ARCHITECTURE Addendum C).
 * A None interp returns ARK_SUCCESS exactly like C's NULL check.
 *
 * The Hermite quartic/quintic bootstrap recursion goes through C's
 * generic arkInterpEvaluate; with the content already borrowed the
 * Rust port recurses directly on arkInterpEvaluate_Hermite (the
 * dispatch would land there anyway).
 *
 * N_Vector NULL-ness follows the crate convention (empty = NULL);
 * allocation failures cannot occur, so C's arkInterpFree-on-failure
 * paths in Init have no translation.
 * -----------------------------------------------------------------*/

use crate::arkode::{arkAllocVec, arkFreeVec, arkResizeVec};
use crate::arkode_impl::{
    ark_step_fullrhs_yn_fn, arkProcessError, ARKInterp, ARKVecResizeFn, ARKodeMem,
    ARK_FULLRHS_END, ARK_FULLRHS_OTHER, ARK_FULLRHS_START, ARK_ILL_INPUT, ARK_INTERP_FAIL,
    ARK_INTERP_MAX_DEGREE, ARK_MEM_FAIL, ARK_RHSFUNC_FAIL, ARK_SUCCESS, FOUR, FUZZ_FACTOR,
    HALF, ONE, THREE as THREE_ARK, TWO, ZERO,
};
use crate::arkode_interp_impl::{
    ARKInterpContent_Hermite, ARKInterpContent_Lagrange, FOURTH, SIX, THREE, TWELVE,
};
use crate::nvector_serial::{
    N_VConst, N_VLinearCombination, N_VLinearSum, N_VScale, NVector,
};
use crate::sundials_math::{SUNMIN, SUNRabs};
use crate::sundials_types::UserData;
use crate::sundials_utils::fmt_g;

/* integer min/max (C SUNMIN/SUNMAX on ints) */
fn imin(a: i32, b: i32) -> i32 {
    if a < b { a } else { b }
}
fn imax(a: i32, b: i32) -> i32 {
    if a > b { a } else { b }
}

/*---------------------------------------------------------------
  Section I: generic ARKInterp functions provided by all
  interpolation modules
  ---------------------------------------------------------------*/

pub fn arkInterpResize(
    ark_mem: &mut ARKodeMem,
    resize: Option<ARKVecResizeFn>,
    resize_data: &mut UserData,
    lrw_diff: i64,
    liw_diff: i64,
    tmpl: &NVector,
) -> i32 {
    let mut interp = match ark_mem.interp.take() {
        None => return ARK_SUCCESS,
        Some(i) => i,
    };
    let ret = match &mut interp {
        ARKInterp::Hermite(hc) => {
            arkInterpResize_Hermite(ark_mem, hc, resize, resize_data, lrw_diff, liw_diff, tmpl)
        }
        ARKInterp::Lagrange(lc) => {
            arkInterpResize_Lagrange(ark_mem, lc, resize, resize_data, lrw_diff, liw_diff, tmpl)
        }
    };
    ark_mem.interp = Some(interp);
    ret
}

/// C arkInterpFree: runs the implementation's free (bookkeeping) and
/// drops the structure; ark_mem.interp is left as None.
pub fn arkInterpFree(ark_mem: &mut ARKodeMem) {
    let mut interp = match ark_mem.interp.take() {
        None => return,
        Some(i) => i,
    };
    match &mut interp {
        ARKInterp::Hermite(hc) => arkInterpFree_Hermite(ark_mem, hc),
        ARKInterp::Lagrange(lc) => arkInterpFree_Lagrange(ark_mem, lc),
    }
    drop(interp);
}

pub fn arkInterpPrintMem(interp: Option<&ARKInterp>, outfile: &mut dyn std::io::Write) {
    match interp {
        None => {}
        Some(ARKInterp::Hermite(hc)) => arkInterpPrintMem_Hermite(hc, outfile),
        Some(ARKInterp::Lagrange(lc)) => arkInterpPrintMem_Lagrange(lc, outfile),
    }
}

pub fn arkInterpSetDegree(ark_mem: &mut ARKodeMem, degree: i32) -> i32 {
    let mut interp = match ark_mem.interp.take() {
        None => return ARK_SUCCESS,
        Some(i) => i,
    };
    let ret = match &mut interp {
        ARKInterp::Hermite(hc) => arkInterpSetDegree_Hermite(ark_mem, hc, degree),
        ARKInterp::Lagrange(lc) => arkInterpSetDegree_Lagrange(ark_mem, lc, degree),
    };
    ark_mem.interp = Some(interp);
    ret
}

pub fn arkInterpInit(ark_mem: &mut ARKodeMem, tnew: f64) -> i32 {
    let mut interp = match ark_mem.interp.take() {
        None => return ARK_SUCCESS,
        Some(i) => i,
    };
    let ret = match &mut interp {
        ARKInterp::Hermite(hc) => arkInterpInit_Hermite(ark_mem, hc, tnew),
        ARKInterp::Lagrange(lc) => arkInterpInit_Lagrange(ark_mem, lc, tnew),
    };
    ark_mem.interp = Some(interp);
    ret
}

pub fn arkInterpUpdate(ark_mem: &mut ARKodeMem, tnew: f64) -> i32 {
    let mut interp = match ark_mem.interp.take() {
        None => return ARK_SUCCESS,
        Some(i) => i,
    };
    let ret = match &mut interp {
        ARKInterp::Hermite(hc) => arkInterpUpdate_Hermite(ark_mem, hc, tnew),
        ARKInterp::Lagrange(lc) => arkInterpUpdate_Lagrange(ark_mem, lc, tnew),
    };
    ark_mem.interp = Some(interp);
    ret
}

pub fn arkInterpEvaluate(
    ark_mem: &mut ARKodeMem,
    tau: f64,
    d: i32,
    order: i32,
    yout: &mut NVector,
) -> i32 {
    let mut interp = match ark_mem.interp.take() {
        None => return ARK_SUCCESS,
        Some(i) => i,
    };
    let ret = match &mut interp {
        ARKInterp::Hermite(hc) => arkInterpEvaluate_Hermite(ark_mem, hc, tau, d, order, yout),
        ARKInterp::Lagrange(lc) => arkInterpEvaluate_Lagrange(ark_mem, lc, tau, d, order, yout),
    };
    ark_mem.interp = Some(interp);
    ret
}

/*---------------------------------------------------------------
  Section II: Hermite interpolation module implementation
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  arkInterpCreate_Hermite:

  This routine creates an ARKInterp structure. This returns a
  non-NULL structure if no errors occurred, or a NULL value
  otherwise.
  ---------------------------------------------------------------*/
pub fn arkInterpCreate_Hermite(ark_mem: &mut ARKodeMem, degree: i32) -> Option<ARKInterp> {
    /* check for valid degree */
    if degree < 0 || degree > ARK_INTERP_MAX_DEGREE {
        return None;
    }

    /* create content, and initialize everything to zero/NULL */
    let mut content = ARKInterpContent_Hermite {
        degree: 0,
        fold: NVector::new(0),
        yold: NVector::new(0),
        fa: NVector::new(0),
        fb: NVector::new(0),
        told: 0.0,
        tnew: 0.0,
        h: 0.0,
    };

    /* set maximum interpolant degree */
    content.degree = imin(ARK_INTERP_MAX_DEGREE, degree);

    /* update workspace sizes */
    ark_mem.lrw += 2;
    ark_mem.liw += 5;

    /* initialize time values */
    content.told = ark_mem.tcur;
    content.tnew = ark_mem.tcur;
    content.h = 0.0;

    Some(ARKInterp::Hermite(content))
}

/*---------------------------------------------------------------
  arkInterpResize_Hermite:

  This routine resizes the internal vectors.
  ---------------------------------------------------------------*/
fn arkInterpResize_Hermite(
    ark_mem: &mut ARKodeMem,
    hc: &mut ARKInterpContent_Hermite,
    resize: Option<ARKVecResizeFn>,
    resize_data: &mut UserData,
    lrw_diff: i64,
    liw_diff: i64,
    y0: &NVector,
) -> i32 {
    /* resize vectors */
    if !arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, y0, &mut hc.fold) {
        return ARK_MEM_FAIL;
    }
    if !arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, y0, &mut hc.yold) {
        return ARK_MEM_FAIL;
    }
    if !arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, y0, &mut hc.fa) {
        return ARK_MEM_FAIL;
    }
    if !arkResizeVec(ark_mem, resize, resize_data, lrw_diff, liw_diff, y0, &mut hc.fb) {
        return ARK_MEM_FAIL;
    }

    /* reinitialize time values */
    hc.told = ark_mem.tcur;
    hc.tnew = ark_mem.tcur;
    hc.h = 0.0;

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkInterpFree_Hermite:

  Workspace bookkeeping for freeing the Hermite structure (the
  memory itself is dropped by the caller).
  ---------------------------------------------------------------*/
fn arkInterpFree_Hermite(ark_mem: &mut ARKodeMem, hc: &mut ARKInterpContent_Hermite) {
    if !hc.fold.data.is_empty() {
        arkFreeVec(ark_mem, &mut hc.fold);
    }
    if !hc.yold.data.is_empty() {
        arkFreeVec(ark_mem, &mut hc.yold);
    }
    if !hc.fa.data.is_empty() {
        arkFreeVec(ark_mem, &mut hc.fa);
    }
    if !hc.fb.data.is_empty() {
        arkFreeVec(ark_mem, &mut hc.fb);
    }

    /* update work space sizes */
    ark_mem.lrw -= 2;
    ark_mem.liw -= 5;
}

/*---------------------------------------------------------------
  arkInterpPrintMem_Hermite
  ---------------------------------------------------------------*/
fn arkInterpPrintMem_Hermite(hc: &ARKInterpContent_Hermite, outfile: &mut dyn std::io::Write) {
    let _ = write!(outfile, "arkode_interp (Hermite): degree = {}\n", hc.degree);
    let _ = write!(
        outfile,
        "arkode_interp (Hermite): told = {}\n",
        fmt_g(hc.told, 0, 15)
    );
    let _ = write!(
        outfile,
        "arkode_interp (Hermite): tnew = {}\n",
        fmt_g(hc.tnew, 0, 15)
    );
    let _ = write!(
        outfile,
        "arkode_interp (Hermite): h = {}\n",
        fmt_g(hc.h, 0, 15)
    );
}

/*---------------------------------------------------------------
  arkInterpSetDegree_Hermite
  ---------------------------------------------------------------*/
fn arkInterpSetDegree_Hermite(
    ark_mem: &mut ARKodeMem,
    hc: &mut ARKInterpContent_Hermite,
    degree: i32,
) -> i32 {
    if degree > ARK_INTERP_MAX_DEGREE || degree < 0 {
        arkProcessError(
            Some(ark_mem),
            ARK_INTERP_FAIL,
            line!(),
            "arkInterpSetDegree_Hermite",
            file!(),
            "Illegal degree specified.",
        );
        return ARK_ILL_INPUT;
    }

    hc.degree = degree;

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkInterpInit_Hermite

  This routine performs the following steps:
  1. Sets tnew and told to the input time
  2. Allocates any missing/needed N_Vector storage (for reinit)
  3. Signals that full RHS data is required for interpolation
  ---------------------------------------------------------------*/
fn arkInterpInit_Hermite(
    ark_mem: &mut ARKodeMem,
    hc: &mut ARKInterpContent_Hermite,
    tnew: f64,
) -> i32 {
    /* initialize time values */
    hc.told = tnew;
    hc.tnew = tnew;
    hc.h = 0.0;

    /* allocate vectors based on interpolant degree */
    let yn_len = ark_mem.yn.data.len();
    if hc.fold.data.is_empty() {
        arkAllocVec(ark_mem, yn_len, &mut hc.fold);
    }
    if hc.yold.data.is_empty() {
        arkAllocVec(ark_mem, yn_len, &mut hc.yold);
    }
    if (hc.degree > 3) && hc.fa.data.is_empty() {
        arkAllocVec(ark_mem, yn_len, &mut hc.fa);
    }
    if (hc.degree > 4) && hc.fb.data.is_empty() {
        arkAllocVec(ark_mem, yn_len, &mut hc.fb);
    }

    /* signal that a full RHS data is required for interpolation */
    ark_mem.call_fullrhs = true;

    /* return with success */
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkInterpUpdate_Hermite

  This routine copies ynew into yold, and fnew into fold, so that
  yold and fold contain the previous values.
  ---------------------------------------------------------------*/
fn arkInterpUpdate_Hermite(
    ark_mem: &mut ARKodeMem,
    hc: &mut ARKInterpContent_Hermite,
    tnew: f64,
) -> i32 {
    /* call full RHS if needed -- called just BEFORE the end of a step, so yn
    has NOT been updated to ycur yet */
    if !ark_mem.fn_is_current {
        let retval = ark_step_fullrhs_yn_fn(ark_mem, ark_mem.tn, ARK_FULLRHS_START);
        if retval != 0 {
            return ARK_RHSFUNC_FAIL;
        }
        ark_mem.fn_is_current = true;
    }

    /* copy ynew and fnew into yold and fold, respectively */
    N_VScale(ONE, &ark_mem.yn, &mut hc.yold);
    N_VScale(ONE, &ark_mem.fn_, &mut hc.fold);

    /* update time values */
    hc.told = hc.tnew;
    hc.tnew = tnew;
    hc.h = ark_mem.h;

    /* return with success */
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkInterpEvaluate_Hermite

  This routine evaluates a temporal interpolation/extrapolation
  based on the data in the interpolation structure (see the C
  source for the full description). The input 'tau' is defined
  over the most-recently-completed interval [told,tnew]:
               t = tnew + tau*(tnew-told).
  ---------------------------------------------------------------*/
fn arkInterpEvaluate_Hermite(
    ark_mem: &mut ARKodeMem,
    hc: &mut ARKInterpContent_Hermite,
    tau: f64,
    d: i32,
    order: i32,
    yout: &mut NVector,
) -> i32 {
    /* local variables */
    let mut a = [0.0_f64; 6];

    /* set constants */
    let tau2 = tau * tau;
    let tau3 = tau * tau2;
    let tau4 = tau * tau3;
    let tau5 = tau * tau4;

    let h = hc.h;
    let h2 = h * h;
    let h3 = h * h2;
    let h4 = h * h3;
    let h5 = h * h4;

    /* determine polynomial order q */
    let mut q = imax(order, 0); /* respect lower bound  */
    q = imin(q, hc.degree); /* respect max possible */

    /* call full RHS if needed -- called just AFTER the end of a step, so yn
    has been updated to ycur */
    if !ark_mem.fn_is_current {
        let retval = ark_step_fullrhs_yn_fn(ark_mem, ark_mem.tn, ARK_FULLRHS_END);
        if retval != 0 {
            return ARK_RHSFUNC_FAIL;
        }
        ark_mem.fn_is_current = true;
    }

    /* error on illegal d */
    if d < 0 {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "arkInterpEvaluate_Hermite",
            file!(),
            "Requested illegal derivative.",
        );
        return ARK_ILL_INPUT;
    }

    /* if d is too high, just return zeros */
    if d > q {
        N_VConst(ZERO, yout);
        return ARK_SUCCESS;
    }

    /* build polynomial based on order */
    match q {
        0 => {
            /* constant interpolant, yout = 0.5*(yn+yp) */
            N_VLinearSum(HALF, &hc.yold, HALF, &ark_mem.yn, yout);
        }

        1 => {
            /* linear interpolant */
            let (a0, a1);
            if d == 0 {
                a0 = -tau;
                a1 = ONE + tau;
            } else {
                /* d=1 */
                a0 = -ONE / h;
                a1 = ONE / h;
            }
            N_VLinearSum(a0, &hc.yold, a1, &ark_mem.yn, yout);
        }

        2 => {
            /* quadratic interpolant */
            if d == 0 {
                a[0] = tau2;
                a[1] = ONE - tau2;
                a[2] = h * (tau2 + tau);
            } else if d == 1 {
                a[0] = TWO * tau / h;
                a[1] = -TWO * tau / h;
                a[2] = ONE + TWO * tau;
            } else {
                /* d == 2 */
                a[0] = TWO / h / h;
                a[1] = -TWO / h / h;
                a[2] = TWO / h;
            }
            let x = [&hc.yold, &ark_mem.yn, &ark_mem.fn_];
            N_VLinearCombination(3, &a, &x, yout);
        }

        3 => {
            /* cubic interpolant */
            if d == 0 {
                a[0] = THREE_ARK * tau2 + TWO * tau3;
                a[1] = ONE - THREE_ARK * tau2 - TWO * tau3;
                a[2] = h * (tau2 + tau3);
                a[3] = h * (tau + TWO * tau2 + tau3);
            } else if d == 1 {
                a[0] = SIX * (tau + tau2) / h;
                a[1] = -SIX * (tau + tau2) / h;
                a[2] = TWO * tau + THREE * tau2;
                a[3] = ONE + FOUR * tau + THREE * tau2;
            } else if d == 2 {
                a[0] = SIX * (ONE + TWO * tau) / h2;
                a[1] = -SIX * (ONE + TWO * tau) / h2;
                a[2] = (TWO + SIX * tau) / h;
                a[3] = (FOUR + SIX * tau) / h;
            } else {
                /* d == 3 */
                a[0] = TWELVE / h3;
                a[1] = -TWELVE / h3;
                a[2] = SIX / h2;
                a[3] = SIX / h2;
            }
            let x = [&hc.yold, &ark_mem.yn, &hc.fold, &ark_mem.fn_];
            N_VLinearCombination(4, &a, &x, yout);
        }

        4 => {
            /* quartic interpolant */

            /* first, evaluate cubic interpolant at tau=-1/3 */
            let mut tval = -ONE / THREE_ARK;
            let retval = arkInterpEvaluate_Hermite(ark_mem, hc, tval, 0, 3, yout);
            if retval != 0 {
                return ARK_RHSFUNC_FAIL;
            }

            /* second, evaluate RHS at tau=-1/3, storing the result in fa */
            tval = hc.tnew - h / THREE_ARK;
            let fullrhs = ark_mem.step_fullrhs.unwrap();
            let retval = fullrhs(ark_mem, tval, yout, &mut hc.fa, ARK_FULLRHS_OTHER);
            if retval != 0 {
                return ARK_RHSFUNC_FAIL;
            }

            /* evaluate desired function */
            if d == 0 {
                a[0] = -SIX * tau2 - 16.0 * tau3 - 9.0 * tau4;
                a[1] = ONE + SIX * tau2 + 16.0 * tau3 + 9.0 * tau4;
                a[2] = h * FOURTH * (-FIVE_C * tau2 - 14.0 * tau3 - 9.0 * tau4);
                a[3] = h * (tau + TWO * tau2 + tau3);
                a[4] = h * 27.0 * FOURTH * (-tau4 - TWO * tau3 - tau2);
            } else if d == 1 {
                a[0] = (-TWELVE * tau - 48.0 * tau2 - 36.0 * tau3) / h;
                a[1] = (TWELVE * tau + 48.0 * tau2 + 36.0 * tau3) / h;
                a[2] = HALF * (-FIVE_C * tau - 21.0 * tau2 - 18.0 * tau3);
                a[3] = ONE + FOUR * tau + THREE * tau2;
                a[4] = -27.0 * HALF * (TWO * tau3 + THREE * tau2 + tau);
            } else if d == 2 {
                a[0] = (-TWELVE - 96.0 * tau - 108.0 * tau2) / h2;
                a[1] = (TWELVE + 96.0 * tau + 108.0 * tau2) / h2;
                a[2] = (-FIVE_C * HALF - 21.0 * tau - 27.0 * tau2) / h;
                a[3] = (FOUR + SIX * tau) / h;
                a[4] = (-27.0 * HALF - 81.0 * tau - 81.0 * tau2) / h;
            } else if d == 3 {
                a[0] = (-96.0 - 216.0 * tau) / h3;
                a[1] = (96.0 + 216.0 * tau) / h3;
                a[2] = (-21.0 - 54.0 * tau) / h2;
                a[3] = SIX / h2;
                a[4] = (-81.0 - 162.0 * tau) / h2;
            } else {
                /* d == 4 */
                a[0] = -216.0 / h4;
                a[1] = 216.0 / h4;
                a[2] = -54.0 / h3;
                a[3] = ZERO;
                a[4] = -162.0 / h3;
            }
            let x = [&hc.yold, &ark_mem.yn, &hc.fold, &ark_mem.fn_, &hc.fa];
            N_VLinearCombination(5, &a, &x, yout);
        }

        5 => {
            /* quintic interpolant */

            /* first, evaluate quartic interpolant at tau=-1/3 */
            let mut tval = -ONE / THREE_ARK;
            let retval = arkInterpEvaluate_Hermite(ark_mem, hc, tval, 0, 4, yout);
            if retval != 0 {
                return ARK_RHSFUNC_FAIL;
            }

            /* second, evaluate RHS at tau=-1/3, storing the result in fa */
            tval = hc.tnew - h / THREE_ARK;
            let fullrhs = ark_mem.step_fullrhs.unwrap();
            let retval = fullrhs(ark_mem, tval, yout, &mut hc.fa, ARK_FULLRHS_OTHER);
            if retval != 0 {
                return ARK_RHSFUNC_FAIL;
            }

            /* third, evaluate quartic interpolant at tau=-2/3 */
            tval = -TWO / THREE_ARK;
            let retval = arkInterpEvaluate_Hermite(ark_mem, hc, tval, 0, 4, yout);
            if retval != 0 {
                return ARK_RHSFUNC_FAIL;
            }

            /* fourth, evaluate RHS at tau=-2/3, storing the result in fb */
            tval = hc.tnew - h * TWO / THREE_ARK;
            let fullrhs = ark_mem.step_fullrhs.unwrap();
            let retval = fullrhs(ark_mem, tval, yout, &mut hc.fb, ARK_FULLRHS_OTHER);
            if retval != 0 {
                return ARK_RHSFUNC_FAIL;
            }

            /* evaluate desired function */
            if d == 0 {
                a[0] = 54.0 * tau5 + 135.0 * tau4 + 110.0 * tau3 + 30.0 * tau2;
                a[1] = ONE - a[0];
                a[2] = h / FOUR * (27.0 * tau5 + 63.0 * tau4 + 49.0 * tau3 + 13.0 * tau2);
                a[3] = h / FOUR
                    * (27.0 * tau5 + 72.0 * tau4 + 67.0 * tau3 + 26.0 * tau2 + FOUR * tau);
                a[4] = h / FOUR * (81.0 * tau5 + 189.0 * tau4 + 135.0 * tau3 + 27.0 * tau2);
                a[5] = h / FOUR * (81.0 * tau5 + 216.0 * tau4 + 189.0 * tau3 + 54.0 * tau2);
            } else if d == 1 {
                a[0] = (270.0 * tau4 + 540.0 * tau3 + 330.0 * tau2 + 60.0 * tau) / h;
                a[1] = -a[0];
                a[2] = (135.0 * tau4 + 252.0 * tau3 + 147.0 * tau2 + 26.0 * tau) / FOUR;
                a[3] = (135.0 * tau4 + 288.0 * tau3 + 201.0 * tau2 + 52.0 * tau + FOUR) / FOUR;
                a[4] = (405.0 * tau4 + 4.0 * 189.0 * tau3 + 405.0 * tau2 + 54.0 * tau) / FOUR;
                a[5] = (405.0 * tau4 + 864.0 * tau3 + 567.0 * tau2 + 108.0 * tau) / FOUR;
            } else if d == 2 {
                a[0] = (1080.0 * tau3 + 1620.0 * tau2 + 660.0 * tau + 60.0) / h2;
                a[1] = -a[0];
                a[2] = (270.0 * tau3 + 378.0 * tau2 + 147.0 * tau + 13.0) / (TWO * h);
                a[3] = (270.0 * tau3 + 432.0 * tau2 + 201.0 * tau + 26.0) / (TWO * h);
                a[4] = (810.0 * tau3 + 1134.0 * tau2 + 405.0 * tau + 27.0) / (TWO * h);
                a[5] = (810.0 * tau3 + 1296.0 * tau2 + 567.0 * tau + 54.0) / (TWO * h);
            } else if d == 3 {
                a[0] = (3240.0 * tau2 + 3240.0 * tau + 660.0) / h3;
                a[1] = -a[0];
                a[2] = (810.0 * tau2 + 756.0 * tau + 147.0) / (TWO * h2);
                a[3] = (810.0 * tau2 + 864.0 * tau + 201.0) / (TWO * h2);
                a[4] = (2430.0 * tau2 + 2268.0 * tau + 405.0) / (TWO * h2);
                a[5] = (2430.0 * tau2 + 2592.0 * tau + 567.0) / (TWO * h2);
            } else if d == 4 {
                a[0] = (6480.0 * tau + 3240.0) / h4;
                a[1] = -a[0];
                a[2] = (810.0 * tau + 378.0) / h3;
                a[3] = (810.0 * tau + 432.0) / h3;
                a[4] = (2430.0 * tau + 1134.0) / h3;
                a[5] = (2430.0 * tau + 1296.0) / h3;
            } else {
                /* d == 5 */
                a[0] = 6480.0 / h5;
                a[1] = -a[0];
                a[2] = 810.0 / h4;
                a[3] = a[2];
                a[4] = 2430.0 / h4;
                a[5] = a[4];
            }
            let x = [
                &hc.yold,
                &ark_mem.yn,
                &hc.fold,
                &ark_mem.fn_,
                &hc.fa,
                &hc.fb,
            ];
            N_VLinearCombination(6, &a, &x, yout);
        }

        _ => {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "arkInterpEvaluate_Hermite",
                file!(),
                "Illegal polynomial order",
            );
            return ARK_ILL_INPUT;
        }
    }

    ARK_SUCCESS
}

/* C FIVE constant from arkode_impl.h (renamed locally to avoid the
   unused-import lint dance with the interp-impl constants) */
const FIVE_C: f64 = crate::arkode_impl::FIVE;

/*---------------------------------------------------------------
  Section III: Lagrange interpolation module implementation
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  arkInterpCreate_Lagrange:

  This routine creates an ARKInterp structure. This returns a
  non-NULL structure if no errors occurred, or a NULL value
  otherwise.
  ---------------------------------------------------------------*/
pub fn arkInterpCreate_Lagrange(ark_mem: &mut ARKodeMem, degree: i32) -> Option<ARKInterp> {
    /* check for valid degree */
    if degree < 0 || degree > ARK_INTERP_MAX_DEGREE {
        return None;
    }

    /* create content, and initialize everything to zero/NULL */
    let mut content = ARKInterpContent_Lagrange {
        nmax: 0,
        nmaxalloc: 0,
        yhist: Vec::new(),
        thist: Vec::new(),
        nhist: 0,
        tround: 0.0,
    };

    /* maximum/current history length */
    content.nmax = imin(degree + 1, ARK_INTERP_MAX_DEGREE + 1); /* respect maximum possible */
    content.nmaxalloc = 0;
    content.nhist = 0;

    /* initial t roundoff value */
    content.tround = FUZZ_FACTOR * ark_mem.uround;

    /* update workspace sizes */
    ark_mem.lrw += (content.nmax + 1) as i64;
    ark_mem.liw += (content.nmax + 2) as i64;

    Some(ARKInterp::Lagrange(content))
}

/*---------------------------------------------------------------
  arkInterpResize_Lagrange:

  This routine resizes the internal vectors.
  ---------------------------------------------------------------*/
fn arkInterpResize_Lagrange(
    ark_mem: &mut ARKodeMem,
    lc: &mut ARKInterpContent_Lagrange,
    resize: Option<ARKVecResizeFn>,
    resize_data: &mut UserData,
    lrw_diff: i64,
    liw_diff: i64,
    y0: &NVector,
) -> i32 {
    /* resize vectors */
    if !lc.yhist.is_empty() {
        for i in 0..lc.nmaxalloc as usize {
            if !arkResizeVec(
                ark_mem,
                resize,
                resize_data,
                lrw_diff,
                liw_diff,
                y0,
                &mut lc.yhist[i],
            ) {
                return ARK_MEM_FAIL;
            }
        }
    }

    /* reset active history length */
    lc.nhist = 0;

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkInterpFree_Lagrange:

  Workspace bookkeeping for freeing the Lagrange structure (the
  memory itself is dropped by the caller).
  ---------------------------------------------------------------*/
fn arkInterpFree_Lagrange(ark_mem: &mut ARKodeMem, lc: &mut ARKInterpContent_Lagrange) {
    if !lc.yhist.is_empty() {
        for i in 0..lc.nmaxalloc as usize {
            if !lc.yhist[i].data.is_empty() {
                arkFreeVec(ark_mem, &mut lc.yhist[i]);
            }
        }
        lc.yhist = Vec::new();
    }
    if !lc.thist.is_empty() {
        lc.thist = Vec::new();
    }

    /* update work space sizes */
    ark_mem.lrw -= (lc.nmax + 1) as i64;
    ark_mem.liw -= (lc.nmax + 2) as i64;
}

/*---------------------------------------------------------------
  arkInterpPrintMem_Lagrange
  ---------------------------------------------------------------*/
fn arkInterpPrintMem_Lagrange(lc: &ARKInterpContent_Lagrange, outfile: &mut dyn std::io::Write) {
    let _ = write!(outfile, "arkode_interp (Lagrange): nmax = {}\n", lc.nmax);
    let _ = write!(outfile, "arkode_interp (Lagrange): nhist = {}\n", lc.nhist);
    if !lc.thist.is_empty() {
        let _ = write!(outfile, "arkode_interp (Lagrange): thist =");
        for i in 0..lc.nmax as usize {
            let _ = write!(outfile, "  {}", fmt_g(lc.thist[i], 0, 15));
        }
        let _ = write!(outfile, "\n");
    }
    if !lc.yhist.is_empty() {
        let _ = write!(outfile, "arkode_interp (Lagrange): yhist ptrs =");
        for i in 0..lc.nmax as usize {
            /* C prints the heap pointers */
            let _ = write!(outfile, "  {:p}", &lc.yhist[i]);
        }
        let _ = write!(outfile, "\n");
    }
}

/*---------------------------------------------------------------
  arkInterpSetDegree_Lagrange
  ---------------------------------------------------------------*/
fn arkInterpSetDegree_Lagrange(
    ark_mem: &mut ARKodeMem,
    lc: &mut ARKInterpContent_Lagrange,
    degree: i32,
) -> i32 {
    if degree > ARK_INTERP_MAX_DEGREE || degree < 0 {
        arkProcessError(
            Some(ark_mem),
            ARK_INTERP_FAIL,
            line!(),
            "arkInterpSetDegree_Lagrange",
            file!(),
            "Illegal degree specified.",
        );
        return ARK_ILL_INPUT;
    }

    lc.nmax = degree + 1;

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkInterpInit_Lagrange

  This routine performs the following steps:
  1. allocates any missing/needed (t,y) history arrays
  2. zeros out stored (t,y) history
  3. copies current (t,y) from main ARKODE memory into history
  4. updates the 'active' history counter to 1
  ---------------------------------------------------------------*/
fn arkInterpInit_Lagrange(
    ark_mem: &mut ARKodeMem,
    lc: &mut ARKInterpContent_Lagrange,
    tnew: f64,
) -> i32 {
    /* check if storage has increased since the last init */
    if lc.nmax > lc.nmaxalloc {
        if !lc.thist.is_empty() {
            lc.thist = Vec::new();
        }
        if !lc.yhist.is_empty() {
            for i in 0..lc.nmaxalloc as usize {
                if !lc.yhist[i].data.is_empty() {
                    arkFreeVec(ark_mem, &mut lc.yhist[i]);
                }
            }
            lc.yhist = Vec::new();
        }
    }

    /* allocate storage for time and solution histories */
    if lc.thist.is_empty() {
        lc.thist = vec![0.0; lc.nmax as usize];
    }

    /* solution history allocation */
    if lc.yhist.is_empty() {
        let yn_len = ark_mem.yn.data.len();
        lc.yhist = (0..lc.nmax as usize).map(|_| NVector::new(0)).collect();
        for i in 0..lc.nmax as usize {
            arkAllocVec(ark_mem, yn_len, &mut lc.yhist[i]);
        }
    }

    /* update allocated size if necessary */
    if lc.nmax > lc.nmaxalloc {
        lc.nmaxalloc = lc.nmax;
    }

    /* zero out history (to be safe) */
    for i in 0..lc.nmaxalloc as usize {
        lc.thist[i] = 0.0;
    }
    for i in 0..lc.nmaxalloc as usize {
        /* N_VConstVectorArray(nmaxalloc, 0.0, yhist) */
        N_VConst(0.0, &mut lc.yhist[i]);
    }

    /* set current time and state as first entries of (t,y) history, update
    counter */
    lc.thist[0] = tnew;
    N_VScale(ONE, &ark_mem.yn, &mut lc.yhist[0]);
    lc.nhist = 1;

    /* return with success */
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkInterpUpdate_Lagrange

  If the current time is 'different enough' from the stored
  values, shifts the (t,y) history and prepends (tnew, ycur);
  otherwise it just returns with success.
  ---------------------------------------------------------------*/
fn arkInterpUpdate_Lagrange(
    ark_mem: &mut ARKodeMem,
    lc: &mut ARKInterpContent_Lagrange,
    tnew: f64,
) -> i32 {
    /* set readability shortcuts */
    let nhist = lc.nhist;
    let nmax = lc.nmax as usize;

    /* update t roundoff value */
    lc.tround = FUZZ_FACTOR * ark_mem.uround * (SUNRabs(ark_mem.tcur) + SUNRabs(ark_mem.h));

    /* determine if tnew differs sufficiently from stored values */
    let mut tdiff = SUNRabs(tnew - lc.thist[0]);
    for i in 1..nhist as usize {
        tdiff = SUNMIN(tdiff, SUNRabs(tnew - lc.thist[i]));
    }
    if tdiff <= lc.tround {
        return ARK_SUCCESS;
    }

    /* shift (t,y) history arrays by one (the C code rotates the
    y-vector pointer array: last slot's storage moves to the front) */
    lc.yhist[..nmax].rotate_right(1);
    for i in (1..nmax).rev() {
        lc.thist[i] = lc.thist[i - 1];
    }

    /* copy tnew and ycur into first entry of history arrays */
    lc.thist[0] = tnew;
    N_VScale(ONE, &ark_mem.ycur, &mut lc.yhist[0]);

    /* update 'nhist' (first few steps) */
    lc.nhist = imin(nhist + 1, nmax as i32);

    /* return with success */
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkInterpEvaluate_Lagrange

  This routine evaluates a temporal interpolation/extrapolation
  based on the stored solution data. The input 'tau' is defined
  over the most-recently-completed interval [t1,t0]:
               t = t0 + tau*(t0-t1).
  ---------------------------------------------------------------*/
fn arkInterpEvaluate_Lagrange(
    ark_mem: &mut ARKodeMem,
    lc: &mut ARKInterpContent_Lagrange,
    tau: f64,
    deriv: i32,
    degree: i32,
    yout: &mut NVector,
) -> i32 {
    /* local variables */
    let mut a = [0.0_f64; 6];

    /* set readability shortcuts */
    let nhist = lc.nhist;

    /* determine polynomial degree q */
    let mut q = imax(degree, 0); /* respect lower bound */
    q = imin(q, nhist - 1); /* respect max possible */

    /* error on illegal deriv */
    if !(0..=3).contains(&deriv) {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "arkInterpEvaluate_Lagrange",
            file!(),
            "Requested illegal derivative.",
        );
        return ARK_ILL_INPUT;
    }

    /* if deriv is too high, just return zeros */
    if deriv > q {
        N_VConst(ZERO, yout);
        return ARK_SUCCESS;
    }

    /* if constant interpolant is requested, just return ynew */
    if q == 0 {
        N_VScale(ONE, &lc.yhist[0], yout);
        return ARK_SUCCESS;
    }

    /* convert from tau back to t (both tnew and told are valid since
    q > 0 => NHIST > 1) */
    let tval = lc.thist[0] + tau * (lc.thist[0] - lc.thist[1]);

    /* linear interpolant */
    if q == 1 {
        if deriv == 0 {
            a[0] = LBasis(lc, 0, tval);
            a[1] = LBasis(lc, 1, tval);
        } else {
            /* deriv == 1 */
            a[0] = LBasisD(lc, 0, tval);
            a[1] = LBasisD(lc, 1, tval);
        }
        N_VLinearSum(a[0], &lc.yhist[0], a[1], &lc.yhist[1], yout);
        return ARK_SUCCESS;
    }

    /* higher-degree interpolant */
    /*    construct linear combination coefficients based on derivative
    requested */
    match deriv {
        0 => {
            /* p(t) */
            for (j, aj) in a.iter_mut().enumerate().take(q as usize + 1) {
                *aj = LBasis(lc, j as i32, tval);
            }
        }
        1 => {
            /* p'(t) */
            for (j, aj) in a.iter_mut().enumerate().take(q as usize + 1) {
                *aj = LBasisD(lc, j as i32, tval);
            }
        }
        2 => {
            /* p''(t) */
            for (j, aj) in a.iter_mut().enumerate().take(q as usize + 1) {
                *aj = LBasisD2(lc, j as i32, tval);
            }
        }
        3 => {
            /* p'''(t) */
            for (j, aj) in a.iter_mut().enumerate().take(q as usize + 1) {
                *aj = LBasisD3(lc, j as i32, tval);
            }
        }
        _ => {}
    }

    /*    evaluate the linear combination and return */
    let x: Vec<&NVector> = lc.yhist[..q as usize + 1].iter().collect();
    N_VLinearCombination(q + 1, &a, &x, yout);

    ARK_SUCCESS
}

/* Lagrange utility routines (basis functions and their derivatives) */
#[allow(non_snake_case)]
fn LBasis(lc: &ARKInterpContent_Lagrange, j: i32, t: f64) -> f64 {
    let mut p = ONE;
    for k in 0..lc.nhist {
        if k == j {
            continue;
        }
        p *= (t - lc.thist[k as usize]) / (lc.thist[j as usize] - lc.thist[k as usize]);
    }
    p
}

#[allow(non_snake_case)]
fn LBasisD(lc: &ARKInterpContent_Lagrange, j: i32, t: f64) -> f64 {
    let mut p = ZERO;
    for i in 0..lc.nhist {
        if i == j {
            continue;
        }
        let mut q = ONE;
        for k in 0..lc.nhist {
            if k == j {
                continue;
            }
            if k == i {
                continue;
            }
            q *= (t - lc.thist[k as usize]) / (lc.thist[j as usize] - lc.thist[k as usize]);
        }
        p += q / (lc.thist[j as usize] - lc.thist[i as usize]);
    }

    p
}

#[allow(non_snake_case)]
fn LBasisD2(lc: &ARKInterpContent_Lagrange, j: i32, t: f64) -> f64 {
    let mut p = ZERO;
    for l in 0..lc.nhist {
        if l == j {
            continue;
        }
        let mut q = ZERO;
        for i in 0..lc.nhist {
            if i == j {
                continue;
            }
            if i == l {
                continue;
            }
            let mut r = ONE;
            for k in 0..lc.nhist {
                if k == j {
                    continue;
                }
                if k == i {
                    continue;
                }
                if k == l {
                    continue;
                }
                r *= (t - lc.thist[k as usize]) / (lc.thist[j as usize] - lc.thist[k as usize]);
            }
            q += r / (lc.thist[j as usize] - lc.thist[i as usize]);
        }
        p += q / (lc.thist[j as usize] - lc.thist[l as usize]);
    }

    p
}

#[allow(non_snake_case)]
fn LBasisD3(lc: &ARKInterpContent_Lagrange, j: i32, t: f64) -> f64 {
    let mut p = ZERO;
    for m in 0..lc.nhist {
        if m == j {
            continue;
        }
        let mut q = ZERO;
        for l in 0..lc.nhist {
            if l == j {
                continue;
            }
            if l == m {
                continue;
            }
            let mut r = ZERO;
            for i in 0..lc.nhist {
                if i == j {
                    continue;
                }
                if i == m {
                    continue;
                }
                if i == l {
                    continue;
                }
                let mut s = ONE;
                for k in 0..lc.nhist {
                    if k == j {
                        continue;
                    }
                    if k == m {
                        continue;
                    }
                    if k == l {
                        continue;
                    }
                    if k == i {
                        continue;
                    }
                    s *= (t - lc.thist[k as usize])
                        / (lc.thist[j as usize] - lc.thist[k as usize]);
                }
                r += s / (lc.thist[j as usize] - lc.thist[i as usize]);
            }
            q += r / (lc.thist[j as usize] - lc.thist[l as usize]);
        }
        p += q / (lc.thist[j as usize] - lc.thist[m as usize]);
    }

    p
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arkode_impl::ARK_INTERP_HERMITE;

    /* full RHS for y' = 3t^2 (y = t^3): independent of y, so the
    Hermite bootstrap stages reproduce the exact polynomial */
    fn cubic_rhs(_ark_mem: &mut ARKodeMem, t: f64, _y: &NVector, f: &mut NVector, _m: i32) -> i32 {
        f.data[0] = 3.0 * t * t;
        0
    }

    fn hermite_mem(degree: i32) -> ARKodeMem {
        let mut ark_mem = ARKodeMem::default();
        ark_mem.uround = crate::sundials_types::SUN_UNIT_ROUNDOFF;
        ark_mem.lrw1 = 1;
        ark_mem.liw1 = 2;
        ark_mem.interp_type = ARK_INTERP_HERMITE;
        ark_mem.step_fullrhs = Some(cubic_rhs);
        /* state at t = 1 for y = t^3 */
        ark_mem.tn = 1.0;
        ark_mem.tcur = 1.0;
        ark_mem.yn = NVector::new(1);
        ark_mem.yn.data[0] = 1.0;
        ark_mem.fn_ = NVector::new(1);
        ark_mem.fn_is_current = false;
        ark_mem.h = 1.0;
        ark_mem.interp = arkInterpCreate_Hermite(&mut ark_mem, degree);
        assert!(ark_mem.interp.is_some());
        assert_eq!(arkInterpInit(&mut ark_mem, 1.0), ARK_SUCCESS);
        /* take the step from t=1 to t=2: Update runs just before the
        end of the step (yn still holds y(1)) */
        assert_eq!(arkInterpUpdate(&mut ark_mem, 2.0), ARK_SUCCESS);
        /* complete the step: yn <- y(2) = 8, fn stale */
        ark_mem.tn = 2.0;
        ark_mem.tcur = 2.0;
        ark_mem.yn.data[0] = 8.0;
        ark_mem.fn_is_current = false;
        ark_mem
    }

    #[test]
    fn hermite_cubic_reproduces_t3() {
        let mut ark_mem = hermite_mem(3);
        let mut yout = NVector::new(1);
        /* t = tnew + tau*h: tau = -0.5 -> t = 1.5 */
        assert_eq!(arkInterpEvaluate(&mut ark_mem, -0.5, 0, 3, &mut yout), ARK_SUCCESS);
        assert!((yout.data[0] - 1.5_f64.powi(3)).abs() < 1e-12, "{}", yout.data[0]);
        /* first derivative: 3 t^2 at t = 1.5 */
        assert_eq!(arkInterpEvaluate(&mut ark_mem, -0.5, 1, 3, &mut yout), ARK_SUCCESS);
        assert!((yout.data[0] - 6.75).abs() < 1e-12, "{}", yout.data[0]);
        /* second derivative: 6 t at t = 1.5 */
        assert_eq!(arkInterpEvaluate(&mut ark_mem, -0.5, 2, 3, &mut yout), ARK_SUCCESS);
        assert!((yout.data[0] - 9.0).abs() < 1e-12, "{}", yout.data[0]);
        /* d > q returns zeros */
        assert_eq!(arkInterpEvaluate(&mut ark_mem, -0.5, 4, 3, &mut yout), ARK_SUCCESS);
        assert_eq!(yout.data[0], 0.0);
        /* illegal derivative */
        assert_eq!(arkInterpEvaluate(&mut ark_mem, -0.5, -1, 3, &mut yout), ARK_ILL_INPUT);
    }

    #[test]
    fn hermite_quintic_bootstrap_reproduces_t3() {
        let mut ark_mem = hermite_mem(5);
        let mut yout = NVector::new(1);
        /* quartic (uses fa via one RHS bootstrap) */
        assert_eq!(arkInterpEvaluate(&mut ark_mem, -0.25, 0, 4, &mut yout), ARK_SUCCESS);
        assert!((yout.data[0] - 1.75_f64.powi(3)).abs() < 1e-10, "{}", yout.data[0]);
        /* quintic (uses fa and fb) */
        assert_eq!(arkInterpEvaluate(&mut ark_mem, -0.25, 0, 5, &mut yout), ARK_SUCCESS);
        assert!((yout.data[0] - 1.75_f64.powi(3)).abs() < 1e-10, "{}", yout.data[0]);
        /* free: workspace bookkeeping returns to the pre-create state */
        let lrw_before = ark_mem.lrw;
        arkInterpFree(&mut ark_mem);
        assert!(ark_mem.interp.is_none());
        assert!(ark_mem.lrw < lrw_before);
    }

    #[test]
    fn lagrange_quadratic_reproduces_t2() {
        let mut ark_mem = ARKodeMem::default();
        ark_mem.uround = crate::sundials_types::SUN_UNIT_ROUNDOFF;
        ark_mem.lrw1 = 1;
        ark_mem.liw1 = 2;
        /* y = t^2 sampled at t = 0 */
        ark_mem.tcur = 0.0;
        ark_mem.yn = NVector::new(1);
        ark_mem.ycur = NVector::new(1);
        ark_mem.h = 1.0;
        ark_mem.interp = arkInterpCreate_Lagrange(&mut ark_mem, 2);
        assert!(ark_mem.interp.is_some());
        assert_eq!(arkInterpInit(&mut ark_mem, 0.0), ARK_SUCCESS);
        /* two updates: (t, y) = (1, 1), (2, 4) */
        ark_mem.tcur = 1.0;
        ark_mem.ycur.data[0] = 1.0;
        assert_eq!(arkInterpUpdate(&mut ark_mem, 1.0), ARK_SUCCESS);
        ark_mem.tcur = 2.0;
        ark_mem.ycur.data[0] = 4.0;
        assert_eq!(arkInterpUpdate(&mut ark_mem, 2.0), ARK_SUCCESS);

        let mut yout = NVector::new(1);
        /* t = t0 + tau*(t0 - t1): tau = -0.5 -> t = 1.5 */
        assert_eq!(arkInterpEvaluate(&mut ark_mem, -0.5, 0, 2, &mut yout), ARK_SUCCESS);
        assert!((yout.data[0] - 2.25).abs() < 1e-12, "{}", yout.data[0]);
        /* p'(1.5) = 3, p''(1.5) = 2, p''' of quadratic = 0 (deriv > q) */
        assert_eq!(arkInterpEvaluate(&mut ark_mem, -0.5, 1, 2, &mut yout), ARK_SUCCESS);
        assert!((yout.data[0] - 3.0).abs() < 1e-12, "{}", yout.data[0]);
        assert_eq!(arkInterpEvaluate(&mut ark_mem, -0.5, 2, 2, &mut yout), ARK_SUCCESS);
        assert!((yout.data[0] - 2.0).abs() < 1e-12, "{}", yout.data[0]);
        assert_eq!(arkInterpEvaluate(&mut ark_mem, -0.5, 3, 2, &mut yout), ARK_SUCCESS);
        assert_eq!(yout.data[0], 0.0);

        /* a repeated update within tround is skipped (nhist unchanged) */
        let nhist_before = match ark_mem.interp.as_ref().unwrap() {
            ARKInterp::Lagrange(lc) => lc.nhist,
            _ => unreachable!(),
        };
        assert_eq!(arkInterpUpdate(&mut ark_mem, 2.0), ARK_SUCCESS);
        let nhist_after = match ark_mem.interp.as_ref().unwrap() {
            ARKInterp::Lagrange(lc) => lc.nhist,
            _ => unreachable!(),
        };
        assert_eq!(nhist_before, nhist_after);

        /* q = 1 with 3 stored points: C's LBasis spans the FULL nhist
        history, so the two coefficients are 3-point basis values —
        0.375*4 + 0.75*1 = 2.25 (not the 2-point linear 2.5) */
        assert_eq!(arkInterpEvaluate(&mut ark_mem, -0.5, 0, 1, &mut yout), ARK_SUCCESS);
        assert!((yout.data[0] - 2.25).abs() < 1e-12, "{}", yout.data[0]);
    }
}
