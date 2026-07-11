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

/*---------------------------------------------------------------
  arkCreate:

  Create and set default values in the ARKodeMem structure. The C
  routine returns NULL on allocation failure; allocation cannot
  fail here (and ARKodeSetDefaults cannot fail with no stepper
  attached), so the memory is returned directly.
  ---------------------------------------------------------------*/
pub fn arkCreate(_sunctx: &crate::sundials_context::SUNContext) -> Box<ARKodeMem> {
    /* C: malloc + memset(0) = ARKodeMem::default() */
    let mut ark_mem = Box::new(ARKodeMem::default());

    /* Set uround */
    ark_mem.uround = crate::sundials_types::SUN_UNIT_ROUNDOFF;

    /* The time step module table, rootfinding, constraints and
    relaxation fields are already NULL/false from Default */

    /* Initialize lrw and liw */
    ark_mem.lrw = 18;
    ark_mem.liw = 53; /* fcn/data ptr, int, long int, sunindextype, sunbooleantype */

    /* Allocate step adaptivity structure and note storage */
    ark_mem.hadapt_mem = Some(crate::arkode_adapt::arkAdaptInit());
    ark_mem.lrw += crate::arkode_adapt_impl::ARK_ADAPT_LRW;
    ark_mem.liw += crate::arkode_adapt_impl::ARK_ADAPT_LIW;

    /* Initialize the interpolation structure to NULL */
    ark_mem.interp = None;
    ark_mem.interp_type = crate::arkode_impl::ARK_INTERP_HERMITE;
    ark_mem.interp_degree = ARK_INTERP_MAX_DEGREE;

    /* Initially, rwt should point to ewt */
    ark_mem.rwt_is_ewt = true;

    /* Indicate that calling the full RHS function is not required, this flag
    is updated to SUNTRUE by the interpolation module initialization function
    and/or the stepper initialization function in arkInitialSetup */
    ark_mem.call_fullrhs = false;

    /* Indicate that the problem needs to be initialized */
    ark_mem.initsetup = true;
    ark_mem.init_type = crate::arkode_impl::FIRST_INIT;
    ark_mem.firststage = true;
    ark_mem.initialized = false;

    /* Initial step size has not been determined yet */
    ark_mem.h = ZERO;
    ark_mem.h0u = ZERO;

    /* Accumulated error estimation strategy */
    ark_mem.AccumErrorType = crate::arkode_impl::ARK_ACCUMERROR_NONE;
    ark_mem.AccumError = ZERO;

    /* Default to having stepper initialize ycur during evolution */
    ark_mem.ensure_ycur = false;

    /* Set default values for integrator and stepper optional inputs
    (cannot fail: no stepper is attached yet) */
    let _ = crate::arkode_io::ARKodeSetDefaults(&mut ark_mem);

    ark_mem.load_checkpoint_fail = false;
    ark_mem.do_adjoint = false;

    ark_mem
}

/*---------------------------------------------------------------
  arkRwtSet

  This routine is responsible for setting the residual weight
  vector rwt (C prototype: (y, weight, data) with data = ark_mem).
  ---------------------------------------------------------------*/
pub fn arkRwtSet(ark_mem: &mut ARKodeMem, y: &NVector, weight: &mut NVector) -> i32 {
    /* return if rwt is just ewt */
    if ark_mem.rwt_is_ewt {
        return 0;
    }

    /* put M*y into ark_tempv1 */
    if let Some(mmult) = ark_mem.step_mmult {
        let mut my = std::mem::replace(&mut ark_mem.tempv1, NVector::new(0));
        let flag = mmult(ark_mem, y, &mut my);
        ark_mem.tempv1 = my;
        if flag != ARK_SUCCESS {
            return crate::arkode_impl::ARK_MASSMULT_FAIL;
        }
    } else {
        /* this condition should not apply, but just in case */
        crate::nvector_serial::N_VScale(1.0, y, &mut ark_mem.tempv1);
    }

    /* call appropriate routine to fill rwt */
    let mut flag = 0;
    let my = std::mem::replace(&mut ark_mem.tempv1, NVector::new(0));
    match ark_mem.ritol {
        crate::arkode_impl::ARK_SS => flag = arkRwtSetSS(ark_mem, &my, weight),
        crate::arkode_impl::ARK_SV => flag = arkRwtSetSV(ark_mem, &my, weight),
        _ => {}
    }
    ark_mem.tempv1 = my;

    flag
}

/*---------------------------------------------------------------
  arkEwtSetSS / arkEwtSetSV / arkEwtSetSmallReal

  Error weight vector routines (C prototype: (ycur, weight,
  arkode_mem) with arkode_mem = ark_mem for the internal
  functions). Following the cvode donor idiom, the weight is built
  in place in `weight` (C stages through tempv1; the element-wise
  operations are identical).
  ---------------------------------------------------------------*/
pub fn arkEwtSetSS(ark_mem: &ARKodeMem, ycur: &NVector, weight: &mut NVector) -> i32 {
    crate::nvector_serial::N_VAbs(ycur, weight);
    weight.scale_inplace(ark_mem.reltol);
    weight.add_const_inplace(ark_mem.Sabstol);
    if ark_mem.atolmin0 && crate::nvector_serial::N_VMin(weight) <= ZERO {
        return -1;
    }
    weight.invert_inplace();
    0
}

pub fn arkEwtSetSV(ark_mem: &ARKodeMem, ycur: &NVector, weight: &mut NVector) -> i32 {
    crate::nvector_serial::N_VAbs(ycur, weight);
    weight.linear_sum_with(ark_mem.reltol, 1.0, ark_mem.Vabstol.as_ref().unwrap());
    if ark_mem.atolmin0 && crate::nvector_serial::N_VMin(weight) <= ZERO {
        return -1;
    }
    weight.invert_inplace();
    0
}

pub fn arkEwtSetSmallReal(_ycur: &NVector, weight: &mut NVector) -> i32 {
    crate::nvector_serial::N_VConst(crate::sundials_types::SUN_SMALL_REAL, weight);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkRwtSetSS / arkRwtSetSV
  ---------------------------------------------------------------*/
pub fn arkRwtSetSS(ark_mem: &ARKodeMem, my: &NVector, weight: &mut NVector) -> i32 {
    crate::nvector_serial::N_VAbs(my, weight);
    weight.scale_inplace(ark_mem.reltol);
    weight.add_const_inplace(ark_mem.SRabstol);
    if ark_mem.Ratolmin0 && crate::nvector_serial::N_VMin(weight) <= ZERO {
        return -1;
    }
    weight.invert_inplace();
    0
}

pub fn arkRwtSetSV(ark_mem: &ARKodeMem, my: &NVector, weight: &mut NVector) -> i32 {
    crate::nvector_serial::N_VAbs(my, weight);
    weight.linear_sum_with(ark_mem.reltol, 1.0, ark_mem.VRabstol.as_ref().unwrap());
    if ark_mem.Ratolmin0 && crate::nvector_serial::N_VMin(weight) <= ZERO {
        return -1;
    }
    weight.invert_inplace();
    0
}
