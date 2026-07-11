/* -----------------------------------------------------------------
 * Translation of sundials-7.7.0/src/arkode/arkode.c — PART I:
 * vector-management utilities (arkAllocVec / arkFreeVec /
 * arkResizeVec). The main infrastructure (arkCreate, ARKodeEvolve,
 * arkInitialSetup, arkHin, arkCompleteStep, rootfinding driver,
 * ...) follows in later parts.
 *
 * Vector nullability convention (pinned for the crate): C's
 * `N_Vector v == NULL` is a ZERO-LENGTH NVector (`v.data.is_empty()`)
 * — ARKODE never allocates zero-length problem vectors. N_VClone of
 * a template only needs its length in the serial build, so
 * arkAllocVec takes `tmpl_len` where C takes the template vector
 * (callers pass ark_mem.yn.data.len() etc., avoiding a second
 * ark_mem borrow).
 * -----------------------------------------------------------------*/

use crate::arkode_impl::{
    arkProcessError, ARKVecResizeFn, ARKodeMem, ARK_BAD_T, ARK_INTERP_MAX_DEGREE, ARK_MEM_FAIL,
    ARK_MEM_NULL, ARK_SUCCESS, FUZZ_FACTOR, ZERO,
};
use crate::arkode_interp::arkInterpEvaluate;
use crate::nvector_serial::NVector;
use crate::sundials_math::SUNRabs;
use crate::sundials_types::UserData;
use crate::sundials_utils::fmt_g;

pub const MSG_ARK_RESIZE_FAIL: &str = "Error in user-supplied resize() function.";

/*---------------------------------------------------------------
  ARKodeGetDky:

  This routine computes the k-th derivative of the interpolating
  polynomial at the time t and stores the result in the vector
  dky. This routine internally calls arkInterpEvaluate to perform
  the interpolation.
  ---------------------------------------------------------------*/
pub fn ARKodeGetDky(ark_mem: &mut ARKodeMem, t: f64, k: i32, dky: &mut NVector) -> i32 {
    /* Check all inputs for legality (dky NULL unrepresentable) */
    if ark_mem.interp.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_NULL,
            line!(),
            "ARKodeGetDky",
            file!(),
            "Missing interpolation structure",
        );
        return ARK_MEM_NULL;
    }

    /* Allow for some slack */
    let mut tfuzz =
        FUZZ_FACTOR * ark_mem.uround * (SUNRabs(ark_mem.tcur) + SUNRabs(ark_mem.hold));
    if ark_mem.hold < ZERO {
        tfuzz = -tfuzz;
    }
    let tp = ark_mem.tcur - ark_mem.hold - tfuzz;
    let tn1 = ark_mem.tcur + tfuzz;
    if (t - tp) * (t - tn1) > ZERO {
        arkProcessError(
            Some(ark_mem),
            ARK_BAD_T,
            line!(),
            "ARKodeGetDky",
            file!(),
            &format!(
                "Illegal value for t. t = {} is not between tcur - hold = {} and tcur = {}",
                fmt_g(t, 0, 15),
                fmt_g(ark_mem.tcur - ark_mem.hold, 0, 15),
                fmt_g(ark_mem.tcur, 0, 15)
            ),
        );
        return ARK_BAD_T;
    }

    /* call arkInterpEvaluate to evaluate result */
    let s = (t - ark_mem.tcur) / ark_mem.h;
    let retval = arkInterpEvaluate(ark_mem, s, k, ARK_INTERP_MAX_DEGREE, dky);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            retval,
            line!(),
            "ARKodeGetDky",
            file!(),
            "Error calling arkInterpEvaluate",
        );
        return retval;
    }
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkAllocVec:

  Allocates a single vector (by cloning the template) if the
  target is NULL, updating lrw/liw. SUNTRUE is returned if the
  allocation is successful (or if the target vector already
  exists). (The C clone-failure path cannot occur here.)
  ---------------------------------------------------------------*/
pub fn arkAllocVec(ark_mem: &mut ARKodeMem, tmpl_len: usize, v: &mut NVector) -> bool {
    /* allocate the new vector if necessary */
    if v.data.is_empty() {
        *v = NVector::new(tmpl_len); /* N_VClone(tmpl) */
        ark_mem.lrw += ark_mem.lrw1;
        ark_mem.liw += ark_mem.liw1;
    }
    true
}

/*---------------------------------------------------------------
  arkFreeVec:

  Frees a single vector if non-NULL, updating lrw/liw.
  ---------------------------------------------------------------*/
pub fn arkFreeVec(ark_mem: &mut ARKodeMem, v: &mut NVector) {
    if !v.data.is_empty() {
        *v = NVector::new(0); /* N_VDestroy; *v = NULL */
        ark_mem.lrw -= ark_mem.lrw1;
        ark_mem.liw -= ark_mem.liw1;
    }
}

/*---------------------------------------------------------------
  arkResizeVec:

  Resizes a single vector based on a template vector: with a
  user ARKVecResizeFn the user routine performs the resize,
  otherwise the vector is re-cloned from the template. Updates
  lrw/liw by the given differences. SUNTRUE on success.
  ---------------------------------------------------------------*/
pub fn arkResizeVec(
    ark_mem: &mut ARKodeMem,
    resize: Option<ARKVecResizeFn>,
    resize_data: &mut UserData,
    lrw_diff: i64,
    liw_diff: i64,
    tmpl: &NVector,
    v: &mut NVector,
) -> bool {
    if !v.data.is_empty() {
        match resize {
            None => {
                /* N_VDestroy + N_VClone(tmpl) */
                *v = NVector::new(tmpl.data.len());
            }
            Some(resize) => {
                if resize(v, tmpl, resize_data) != 0 {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_MEM_FAIL,
                        line!(),
                        "arkResizeVec",
                        file!(),
                        MSG_ARK_RESIZE_FAIL,
                    );
                    return false;
                }
            }
        }
        ark_mem.lrw += lrw_diff;
        ark_mem.liw += liw_diff;
    }
    true
}
