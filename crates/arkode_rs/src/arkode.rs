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
    ARK_MEM_NULL, ARK_SUCCESS, FUZZ_FACTOR, HALF, ONE, TWO, ZERO,
};
use crate::arkode_interp::arkInterpEvaluate;
use crate::nvector_serial::{NVector, N_VLinearCombination};
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

/*---------------------------------------------------------------
  arkInit:

  arkInit allocates and initializes memory for a problem. All
  inputs are checked for errors. If any error occurs during
  initialization, an error flag is returned. Otherwise, it
  returns ARK_SUCCESS.  This routine should only be called by
  ARKodeReset (RESET_INIT) or a timestepper module (re-)
  initialization routine (FIRST_INIT); never by the user.
  ---------------------------------------------------------------*/
pub fn arkInit(ark_mem: &mut ARKodeMem, t0: f64, y0: &NVector, init_type: i32) -> i32 {
    use crate::arkode_impl::{ARK_CONTROLLER_ERR, ARK_ILL_INPUT, FIRST_INIT, ONE, RESET_INIT};

    let mut init_type = init_type;

    /* Check for legal input parameters (y0 NULL unrepresentable) */

    /* Check if reset was called before the first Evolve call */
    if init_type == RESET_INIT && !ark_mem.initialized {
        init_type = FIRST_INIT;
    }

    /* Check if allocations have been done i.e., is this first init call */
    if !ark_mem.MallocDone {
        /* Test if all required time stepper operations are implemented */
        if !arkCheckTimestepper(ark_mem) {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "arkInit",
                file!(),
                "Time stepper module is missing required functionality",
            );
            return ARK_ILL_INPUT;
        }

        /* Required vector operations are always present in the serial
        build (arkCheckNvectorRequired) */

        /* Set space requirements for one N_Vector (serial N_VSpace) */
        ark_mem.lrw1 = y0.data.len() as i64;
        ark_mem.liw1 = 1;

        /* Allocate the solver vectors (using y0 as a template) */
        arkAllocVectors(ark_mem, y0.data.len());

        /* All allocations are complete */
        ark_mem.MallocDone = true;
    }

    /* All allocation and error checking is complete at this point */

    /* Copy the input parameters into ARKODE state */
    ark_mem.tcur = t0;
    ark_mem.tn = t0;

    /* Initialize yn */
    crate::nvector_serial::N_VScale(ONE, y0, &mut ark_mem.yn);
    ark_mem.fn_is_current = false;

    /* Clear any previous 'tstop' */
    ark_mem.tstopset = false;

    /* Initializations on (re-)initialization call, skip on reset */
    if init_type == FIRST_INIT {
        /* Counters */
        ark_mem.nst_attempts = 0;
        ark_mem.nst = 0;
        ark_mem.nhnil = 0;
        ark_mem.ncfn = 0;
        ark_mem.netf = 0;
        ark_mem.nconstrfails = 0;

        /* Initial, old, and next step sizes */
        ark_mem.h0u = ZERO;
        ark_mem.hold = ZERO;
        ark_mem.next_h = ZERO;

        /* Tolerance scale factor */
        ark_mem.tolsf = ONE;

        /* Reset error controller object */
        let hadapt_mem = ark_mem.hadapt_mem.as_mut().unwrap();
        if let Some(hcontroller) = hadapt_mem.hcontroller.as_mut() {
            let retval = crate::sundials_adaptcontroller::SUNAdaptController_Reset(hcontroller);
            if retval != crate::sundials_errors::SUN_SUCCESS {
                arkProcessError(
                    Some(ark_mem),
                    ARK_CONTROLLER_ERR,
                    line!(),
                    "arkInit",
                    file!(),
                    "Unable to reset error controller object",
                );
                return ARK_CONTROLLER_ERR;
            }
        }

        /* Adaptivity counters */
        let hadapt_mem = ark_mem.hadapt_mem.as_mut().unwrap();
        hadapt_mem.nst_acc = 0;
        hadapt_mem.nst_exp = 0;

        /* Accumulated error estimate */
        ark_mem.AccumError = ZERO;

        /* Indicate that calling the full RHS function is not required, this
        flag is updated to SUNTRUE by the interpolation module initialization
        function and/or the stepper initialization function in
        arkInitialSetup */
        ark_mem.call_fullrhs = false;

        /* Adjoint related */
        ark_mem.checkpoint_step_idx = 0;

        /* Indicate that initialization has not been done before */
        ark_mem.initialized = false;
    }

    /* Indicate initialization is needed */
    ark_mem.initsetup = true;
    ark_mem.init_type = init_type;
    ark_mem.firststage = true;

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkCheckTimestepper:

  This routine checks if all required time stepper function
  pointers have been supplied.  If any of them is missing it
  returns SUNFALSE.
  ---------------------------------------------------------------*/
pub fn arkCheckTimestepper(ark_mem: &ARKodeMem) -> bool {
    !(ark_mem.step_init.is_none() || ark_mem.step.is_none() || ark_mem.step_mem.is_none())
}

/*---------------------------------------------------------------
  arkCheckNvectorRequired / arkCheckNvectorOptional:

  The serial NVector provides every operation these routines
  test for, so both always hold in this port.
  ---------------------------------------------------------------*/
pub fn arkCheckNvectorRequired(_tmpl: &NVector) -> bool {
    true
}

pub fn arkCheckNvectorOptional(_ark_mem: &ARKodeMem) -> bool {
    true
}

/*---------------------------------------------------------------
  arkAllocVectors:

  This routine allocates the ARKODE vectors ewt, yn, tempv* and
  ftemp. If any of these vectors already exist, they are left
  alone. It also sets the optional outputs lrw and liw, which are
  (respectively) the lengths of the real and integer work spaces.

  rwt aliasing note: when rwt_is_ewt is SUNTRUE the C code sets
  the rwt pointer equal to ewt; here rwt simply stays unallocated
  and readers go through the ewt vector (Addendum C.1).
  ---------------------------------------------------------------*/
pub fn arkAllocVectors(ark_mem: &mut ARKodeMem, tmpl_len: usize) -> bool {
    /* Allocate ewt if needed */
    let mut v = std::mem::replace(&mut ark_mem.ewt, NVector::new(0));
    arkAllocVec(ark_mem, tmpl_len, &mut v);
    ark_mem.ewt = v;

    /* Set rwt to point at ewt: represented by leaving rwt empty */

    /* Allocate yn if needed */
    let mut v = std::mem::replace(&mut ark_mem.yn, NVector::new(0));
    arkAllocVec(ark_mem, tmpl_len, &mut v);
    ark_mem.yn = v;

    /* Allocate tempv1 if needed */
    let mut v = std::mem::replace(&mut ark_mem.tempv1, NVector::new(0));
    arkAllocVec(ark_mem, tmpl_len, &mut v);
    ark_mem.tempv1 = v;

    /* Allocate tempv2 if needed */
    let mut v = std::mem::replace(&mut ark_mem.tempv2, NVector::new(0));
    arkAllocVec(ark_mem, tmpl_len, &mut v);
    ark_mem.tempv2 = v;

    /* Allocate tempv3 if needed */
    let mut v = std::mem::replace(&mut ark_mem.tempv3, NVector::new(0));
    arkAllocVec(ark_mem, tmpl_len, &mut v);
    ark_mem.tempv3 = v;

    /* Allocate tempv4 if needed */
    let mut v = std::mem::replace(&mut ark_mem.tempv4, NVector::new(0));
    arkAllocVec(ark_mem, tmpl_len, &mut v);
    ark_mem.tempv4 = v;

    true
}

/*---------------------------------------------------------------
  arkFreeVectors

  This routine frees the ARKODE vectors allocated in both
  arkAllocVectors and arkAllocVec.
  ---------------------------------------------------------------*/
pub fn arkFreeVectors(ark_mem: &mut ARKodeMem) {
    macro_rules! free_field {
        ($field:ident) => {
            let mut v = std::mem::replace(&mut ark_mem.$field, NVector::new(0));
            arkFreeVec(ark_mem, &mut v);
            ark_mem.$field = v;
        };
    }
    free_field!(ewt);
    if !ark_mem.rwt_is_ewt {
        free_field!(rwt);
    }
    free_field!(tempv1);
    free_field!(tempv2);
    free_field!(tempv3);
    free_field!(tempv4);
    free_field!(tempv5);
    free_field!(yn);
    free_field!(fn_);
    /* Vabstol / constraints are Option<NVector>: drop + bookkeeping */
    if let Some(v) = ark_mem.Vabstol.take() {
        drop(v);
        ark_mem.lrw -= ark_mem.lrw1;
        ark_mem.liw -= ark_mem.liw1;
    }
    if let Some(v) = ark_mem.constraints.take() {
        drop(v);
        ark_mem.lrw -= ark_mem.lrw1;
        ark_mem.liw -= ark_mem.liw1;
    }
}

/*---------------------------------------------------------------
  ark_efun_apply_yn:

  Dispatch point for C calls `ark_mem->efun(ark_mem->yn, w,
  ark_mem->e_data)`: a user efun receives user_data (C sets
  e_data = user_data), the internal path branches on itol
  (C installs arkEwtSetSS/SV with e_data = ark_mem).
  ---------------------------------------------------------------*/
pub(crate) fn ark_efun_apply_yn(ark_mem: &mut ARKodeMem, w: &mut NVector) -> i32 {
    if ark_mem.user_efun {
        let efun = ark_mem.efun.unwrap();
        efun(&ark_mem.yn, w, &mut ark_mem.user_data)
    } else {
        match ark_mem.itol {
            crate::arkode_impl::ARK_SS => arkEwtSetSS(ark_mem, &ark_mem.yn, w),
            crate::arkode_impl::ARK_SV => arkEwtSetSV(ark_mem, &ark_mem.yn, w),
            _ => -1,
        }
    }
}

/*---------------------------------------------------------------
  arkHin

  This routine computes a tentative initial step size h0.  If
  tout is too close to tn (= t0), then arkHin returns
  ARK_TOO_CLOSE and h remains uninitialized (the caller checks
  the distance first).  Otherwise it sets ark_mem->h and returns
  ARK_SUCCESS (see the C source for the full algorithm notes).
  ---------------------------------------------------------------*/
pub fn arkHin(ark_mem: &mut ARKodeMem, tout: f64) -> i32 {
    use crate::arkode_impl::{
        ARK_FULLRHS_START, ARK_REPTD_RHSFUNC_ERR, ARK_RHSFUNC_FAIL, H0_BIAS, H0_ITERS,
        H0_LBFACTOR, HALF, TWO,
    };
    use crate::sundials_math::{SUNMAX, SUNRsqrt};

    /* arkInitialSetup checks for tdiff = 0 or < 2 * troundoff */
    let tdiff = tout - ark_mem.tcur;
    let sign = if tdiff > ZERO { 1 } else { -1 };
    let tdist = SUNRabs(tdiff);
    let tround = ark_mem.uround * SUNMAX(SUNRabs(ark_mem.tcur), SUNRabs(tout));

    /* call full RHS if needed */
    if !ark_mem.fn_is_current {
        /* NOTE: The step size (h) is used in setting the tolerance in a
        potential mass matrix solve when computing the full RHS. Before
        calling arkHin, h is set to |tout - tcur| or 1 and so we do not need
        to guard against h == 0 here before calling the full RHS. */
        let retval =
            crate::arkode_impl::ark_step_fullrhs_yn_fn(ark_mem, ark_mem.tn, ARK_FULLRHS_START);
        if retval != 0 {
            return ARK_RHSFUNC_FAIL;
        }
        ark_mem.fn_is_current = true;
    }

    /* Set lower and upper bounds on h0, and take geometric mean
    as first trial value.
    Exit with this value if the bounds cross each other. */
    let hlb = H0_LBFACTOR * tround;
    let hub = arkUpperBoundH0(ark_mem, tdist);

    let mut hg = SUNRsqrt(hlb * hub);

    if hub < hlb {
        if sign == -1 {
            ark_mem.h = -hg;
        } else {
            ark_mem.h = hg;
        }
        return ARK_SUCCESS;
    }

    /* Outer loop */
    let mut hs = hg; /* safeguard against 'uninitialized variable' warning */
    let mut hnew = hs;
    let mut yddnrm = 0.0;
    let mut count1 = 1;
    while count1 <= H0_ITERS {
        /* Attempts to estimate ydd */
        let mut hg_ok = false;

        for _count2 in 1..=H0_ITERS {
            let hgs = hg * sign as f64;
            let retval = arkYddNorm(ark_mem, hgs, &mut yddnrm);
            /* If f() failed unrecoverably, give up */
            if retval < 0 {
                return ARK_RHSFUNC_FAIL;
            }
            /* If successful, we can use ydd */
            if retval == ARK_SUCCESS {
                hg_ok = true;
                break;
            }
            /* f() failed recoverably; cut step size and test it again */
            hg *= 0.2;
        }

        /* If f() failed recoverably H0_ITERS times */
        if !hg_ok {
            /* Exit if this is the first or second pass. No recovery possible */
            if count1 <= 2 {
                return ARK_REPTD_RHSFUNC_ERR;
            }
            /* We have a fall-back option. The value hs is a previous hnew
            which passed through f(). Use it and break */
            hnew = hs;
            break;
        }

        /* The proposed step size is feasible. Save it. */
        hs = hg;

        /* Propose new step size */
        hnew = if yddnrm * hub * hub > TWO {
            SUNRsqrt(TWO / yddnrm)
        } else {
            SUNRsqrt(hg * hub)
        };

        /* If last pass, stop now with hnew */
        if count1 == H0_ITERS {
            break;
        }

        let hrat = hnew / hg;

        /* Accept hnew if it does not differ from hg by more than a factor of
        2 */
        if (hrat > HALF) && (hrat < TWO) {
            break;
        }

        /* After one pass, if ydd seems to be bad, use fall-back value. */
        if (count1 > 1) && (hrat > TWO) {
            hnew = hg;
            break;
        }

        /* Send this value back through f() */
        hg = hnew;
        count1 += 1;
    }

    /* Apply bounds, bias factor, and attach sign */
    let mut h0 = H0_BIAS * hnew;
    if h0 < hlb {
        h0 = hlb;
    }
    if h0 > hub {
        h0 = hub;
    }
    if sign == -1 {
        h0 = -h0;
    }
    ark_mem.h = h0;

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkUpperBoundH0

  This routine sets an upper bound on abs(h0) based on
  tdist = tn - t0 and the values of y[i]/y'[i].
  ---------------------------------------------------------------*/
pub fn arkUpperBoundH0(ark_mem: &mut ARKodeMem, tdist: f64) -> f64 {
    use crate::arkode_impl::{H0_UBFACTOR, ONE};
    use crate::nvector_serial::{N_VAbs, N_VDiv, N_VMaxNorm};

    /* Bound based on |y0|/|y0'| -- allow at most an increase of
     * H0_UBFACTOR in y0 (based on a forward Euler step). The weight
     * factor is used as a safeguard against zero components in y0. */
    let mut temp1 = std::mem::replace(&mut ark_mem.tempv1, NVector::new(0));
    let mut temp2 = std::mem::replace(&mut ark_mem.tempv2, NVector::new(0));

    N_VAbs(&ark_mem.yn, &mut temp2);
    let _ = ark_efun_apply_yn(ark_mem, &mut temp1);
    temp1.invert_inplace(); /* N_VInv(temp1, temp1) */
    temp1.linear_sum_with(ONE, H0_UBFACTOR, &temp2); /* temp1 = H0_UBFACTOR*temp2 + temp1 */

    N_VAbs(&ark_mem.fn_, &mut temp2);

    let mut quot = NVector::new(temp1.data.len());
    N_VDiv(&temp2, &temp1, &mut quot);
    let hub_inv = N_VMaxNorm(&quot);
    /* (N_VDiv(temp2, temp1, temp1) in C: quotient built in a fresh
    vector here, same element-wise values) */
    temp1 = quot;

    ark_mem.tempv1 = temp1;
    ark_mem.tempv2 = temp2;

    /* bound based on tdist -- allow at most a step of magnitude
     * H0_UBFACTOR * tdist */
    let mut hub = H0_UBFACTOR * tdist;

    /* Use the smaller of the two */
    if hub * hub_inv > ONE {
        hub = ONE / hub_inv;
    }

    hub
}

/*---------------------------------------------------------------
  arkYddNorm

  This routine computes an estimate of the second derivative of y
  using a difference quotient, and returns its WRMS norm.
  ---------------------------------------------------------------*/
pub fn arkYddNorm(ark_mem: &mut ARKodeMem, hg: f64, yddnrm: &mut f64) -> i32 {
    use crate::arkode_impl::{ARK_FULLRHS_OTHER, ARK_RHSFUNC_FAIL, ONE};
    use crate::nvector_serial::{N_VLinearSum, N_VScale, N_VWrmsNorm};

    /* increment y with a multiple of f */
    N_VLinearSum(hg, &ark_mem.fn_, ONE, &ark_mem.yn, &mut ark_mem.ycur);

    /* compute y', via the ODE RHS routine */
    let fullrhs = ark_mem.step_fullrhs.unwrap();
    let ycur = std::mem::replace(&mut ark_mem.ycur, NVector::new(0));
    let mut tempv1 = std::mem::replace(&mut ark_mem.tempv1, NVector::new(0));
    let retval = fullrhs(ark_mem, ark_mem.tcur + hg, &ycur, &mut tempv1, ARK_FULLRHS_OTHER);
    ark_mem.ycur = ycur;
    ark_mem.tempv1 = tempv1;
    if retval != 0 {
        return ARK_RHSFUNC_FAIL;
    }

    /* difference new f and original f to estimate y'' (C output
    aliases the first operand: in-place method family) */
    let fn_ = std::mem::replace(&mut ark_mem.fn_, NVector::new(0));
    ark_mem.tempv1.linear_sum_with(ONE / hg, -ONE / hg, &fn_);
    ark_mem.fn_ = fn_;

    /* reset ycur to equal yn (unnecessary?) */
    let mut ycur = std::mem::replace(&mut ark_mem.ycur, NVector::new(0));
    N_VScale(ONE, &ark_mem.yn, &mut ycur);
    ark_mem.ycur = ycur;

    /* compute norm of y'' */
    *yddnrm = N_VWrmsNorm(&ark_mem.tempv1, &ark_mem.ewt);

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ark_rfun_apply_yn:

  Dispatch point for C calls `ark_mem->rfun(ark_mem->yn, w,
  ark_mem->r_data)` (only reached when !rwt_is_ewt): a user rfun
  receives user_data, the internal path is arkRwtSet.
  ---------------------------------------------------------------*/
pub(crate) fn ark_rfun_apply_yn(ark_mem: &mut ARKodeMem, w: &mut NVector) -> i32 {
    if ark_mem.user_rfun {
        let rfun = ark_mem.rfun.unwrap();
        rfun(&ark_mem.yn, w, &mut ark_mem.user_data)
    } else {
        let yn = std::mem::replace(&mut ark_mem.yn, NVector::new(0));
        let flag = arkRwtSet(ark_mem, &yn, w);
        ark_mem.yn = yn;
        flag
    }
}

/*---------------------------------------------------------------
  arkCompleteStep

  This routine performs various update operations when the step
  solution is complete.  It is assumed that the timestepper
  module has stored the time-evolved solution in ark_mem->ycur,
  and the step that gave rise to this solution in ark_mem->h.
  We update the current time (tn), the current solution (yn),
  increment the overall step counter nst, record the values hold
  and tnew, allow for user-provided postprocessing, and update
  the interpolation structure.
  ---------------------------------------------------------------*/
pub fn arkCompleteStep(ark_mem: &mut ARKodeMem, dsm: f64) -> i32 {
    use crate::arkode_impl::{
        ARK_ACCUMERROR_MAX, ARK_ACCUMERROR_NONE, ARK_ACCUMERROR_SUM, ARK_CONTROLLER_ERR,
        ARK_POSTSTEPFN_FAIL, ONE,
    };
    use crate::sundials_math::SUNMAX;
    use crate::sundials_utils::sunCompensatedSum;

    /* Set current time to the end of the step (in case the last stage time
    does not coincide with the step solution time). If tstop is enabled, it
    is possible for tn + h to be past tstop by roundoff, and in that case,
    we reset tn (after incrementing by h) to tstop. */

    /* During long-time integration, roundoff can creep into tcur.
    Compensated summation fixes this but with increased cost, so it is
    optional. */
    if ark_mem.use_compensated_sums {
        let (tn, h) = (ark_mem.tn, ark_mem.h);
        sunCompensatedSum(tn, h, &mut ark_mem.tcur, &mut ark_mem.terr);
    } else {
        ark_mem.tcur = ark_mem.tn + ark_mem.h;
    }

    if ark_mem.tstopset {
        let troundoff =
            FUZZ_FACTOR * ark_mem.uround * (SUNRabs(ark_mem.tcur) + SUNRabs(ark_mem.h));
        if SUNRabs(ark_mem.tcur - ark_mem.tstop) <= troundoff {
            ark_mem.tcur = ark_mem.tstop;
        }
    }

    /* store this step's contribution to accumulated temporal error */
    if ark_mem.AccumErrorType != ARK_ACCUMERROR_NONE {
        if ark_mem.AccumErrorType == ARK_ACCUMERROR_MAX {
            ark_mem.AccumError = SUNMAX(dsm, ark_mem.AccumError);
        } else if ark_mem.AccumErrorType == ARK_ACCUMERROR_SUM {
            ark_mem.AccumError += dsm;
        } else {
            /* ARK_ACCUMERROR_AVG */
            ark_mem.AccumError += dsm * ark_mem.h;
        }
    }

    /* call the user-supplied post-step function (if supplied) */
    if let Some(post_step_fn) = ark_mem.PostStepFn {
        let retval = post_step_fn(
            ark_mem.tcur,
            &ark_mem.ycur,
            ark_mem.nst,
            &mut ark_mem.user_data,
        );
        if retval != 0 {
            return ARK_POSTSTEPFN_FAIL;
        }
    }

    /* update interpolation structure

    NOTE: This must be called before updating yn with ycur as the
    interpolation module may need to save tn, yn from the start of this
    step. */
    if ark_mem.interp.is_some() {
        let retval = crate::arkode_interp::arkInterpUpdate(ark_mem, ark_mem.tcur);
        if retval != ARK_SUCCESS {
            return retval;
        }
    }

    /* update yn to current solution */
    {
        let ycur = std::mem::replace(&mut ark_mem.ycur, NVector::new(0));
        crate::nvector_serial::N_VScale(ONE, &ycur, &mut ark_mem.yn);
        ark_mem.ycur = ycur;
    }
    ark_mem.fn_is_current = false;

    /* Notify time step controller object of successful step */
    let (h, _eta) = (ark_mem.h, ark_mem.eta);
    if ark_mem.hadapt_mem.as_ref().unwrap().hcontroller.is_some() {
        let hc = ark_mem
            .hadapt_mem
            .as_mut()
            .unwrap()
            .hcontroller
            .as_mut()
            .unwrap();
        let retval = crate::sundials_adaptcontroller::SUNAdaptController_UpdateH(hc, h, dsm);
        if retval != crate::sundials_errors::SUN_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_CONTROLLER_ERR,
                line!(),
                "arkCompleteStep",
                file!(),
                "Failure updating controller object",
            );
            return ARK_CONTROLLER_ERR;
        }
    }

    /* update scalar quantities */
    ark_mem.nst += 1;
    ark_mem.checkpoint_step_idx += 1;
    ark_mem.hold = ark_mem.h;
    ark_mem.tn = ark_mem.tcur;
    ark_mem.hprime = ark_mem.h * ark_mem.eta;

    /* Reset growth factor for subsequent time step */
    {
        let hadapt_mem = ark_mem.hadapt_mem.as_mut().unwrap();
        hadapt_mem.etamax = hadapt_mem.growth;
    }

    /* Turn off flag indicating initial step and first stage */
    ark_mem.initsetup = false;
    ark_mem.firststage = false;

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkHandleFailure

  This routine prints error messages for all cases of failure by
  arkHin and ark_step. It returns to ARKODE the value that ARKODE
  is to return to the user.
  ---------------------------------------------------------------*/
pub fn arkHandleFailure(ark_mem: &mut ARKodeMem, flag: i32) -> i32 {
    use crate::arkode_impl::*;

    let t = fmt_g(ark_mem.tcur, 0, 15);
    let h = fmt_g(ark_mem.h, 0, 15);

    /* Depending on flag, print error message and return error flag */
    let msg: String = match flag {
        ARK_ERR_FAILURE => format!(
            "At t = {} and h = {}, the error test failed repeatedly or with |h| = hmin.",
            t, h
        ),
        ARK_CONV_FAILURE => format!(
            "At t = {} and h = {}, the solver convergence test failed repeatedly or with |h| = hmin.",
            t, h
        ),
        ARK_LSETUP_FAIL => {
            format!("At t = {}, the setup routine failed in an unrecoverable manner.", t)
        }
        ARK_LSOLVE_FAIL => {
            format!("At t = {}, the solve routine failed in an unrecoverable manner.", t)
        }
        ARK_RHSFUNC_FAIL => format!(
            "At t = {}, the right-hand side routine failed in an unrecoverable manner.",
            t
        ),
        ARK_UNREC_RHSFUNC_ERR => format!(
            "At t = {}, the right-hand side failed in a recoverable manner, but no recovery is possible.",
            t
        ),
        ARK_REPTD_RHSFUNC_ERR => {
            format!("At t = {} repeated recoverable right-hand side function errors.", t)
        }
        ARK_RTFUNC_FAIL => format!(
            "At t = {}, the rootfinding routine failed in an unrecoverable manner.",
            t
        ),
        ARK_TOO_CLOSE => "tout too close to t0 to start integration.".to_string(),
        ARK_CONSTR_FAIL => {
            format!("At t = {}, unable to satisfy inequality constraints.", t)
        }
        ARK_MASSSOLVE_FAIL => "The mass matrix solver failed.".to_string(),
        ARK_NLS_SETUP_FAIL => {
            format!("At t = {} the nonlinear solver setup failed unrecoverably", t)
        }
        ARK_VECTOROP_ERR => format!("At t = {}, a vector operation failed.", t),
        ARK_INNERSTEP_FAIL => {
            format!("At t = {}, the inner stepper failed in an unrecoverable manner.", t)
        }
        ARK_NLS_OP_ERR => {
            format!("At t = {} the nonlinear solver failed in an unrecoverable manner.", t)
        }
        ARK_USER_PREDICT_FAIL => format!(
            "At t = {} the user-supplied predictor failed in an unrecoverable manner.",
            t
        ),
        ARK_POSTPROCESS_STEP_FAIL => format!(
            "At t = {}, the step postprocessing routine failed in an unrecoverable manner.",
            t
        ),
        ARK_POSTPROCESS_STAGE_FAIL => format!(
            "At t = {}, the stage postprocessing routine failed in an unrecoverable manner.",
            t
        ),
        ARK_PRESTEPFN_FAIL => format!(
            "At t = {}, the pre-step function failed in an unrecoverable manner.",
            t
        ),
        ARK_POSTSTEPFN_FAIL => format!(
            "At t = {}, the post-step function failed in an unrecoverable manner.",
            t
        ),
        ARK_PRERHSFN_FAIL => format!(
            "At t = {}, the pre-RHS function failed in an unrecoverable manner.",
            t
        ),
        ARK_INTERP_FAIL => {
            format!("At t = {} the interpolation module failed unrecoverably", t)
        }
        ARK_INVALID_TABLE => "ARKODE was provided an invalid method table".to_string(),
        ARK_RELAX_FAIL => format!("At t = {} the relaxation module failed", t),
        ARK_RELAX_MEM_NULL => "The ARKODE relaxation module memory is NULL".to_string(),
        ARK_RELAX_FUNC_FAIL => "The relaxation function failed unrecoverably".to_string(),
        ARK_RELAX_JAC_FAIL => "The relaxation Jacobian failed unrecoverably".to_string(),
        ARK_ADJ_RECOMPUTE_FAIL => {
            "The forward recomputation of step failed unrecoverably".to_string()
        }
        ARK_ADJ_CHECKPOINT_FAIL => "A checkpoint operation failed unrecoverably".to_string(),
        ARK_SUNADJSTEPPER_ERR => {
            "A SUNAdjStepper operation failed unrecoverably".to_string()
        }
        ARK_DOMEIG_FAIL => "The dominant eigenvalue function failed unrecoverably".to_string(),
        ARK_MAX_STAGE_LIMIT_FAIL => "The max stage limit failed unrecoverably".to_string(),
        ARK_SUNSTEPPER_ERR => "An inner SUNStepper error occurred".to_string(),
        _ => {
            /* This return should never happen */
            arkProcessError(
                Some(ark_mem),
                ARK_UNRECOGNIZED_ERROR,
                line!(),
                "arkHandleFailure",
                file!(),
                "ARKODE encountered an unrecognized error. Please report this to the \
                 Sundials developers at sundials-users@llnl.gov",
            );
            return ARK_UNRECOGNIZED_ERROR;
        }
    };
    arkProcessError(Some(ark_mem), flag, line!(), "arkHandleFailure", file!(), &msg);

    flag
}

/*---------------------------------------------------------------
  arkPredict_MaximumOrder

  This routine predicts the nonlinear implicit stage solution
  using the ARKode interpolation module.  This uses the
  highest-degree interpolant supported by the module (stored
  in the interpolation module).
  ---------------------------------------------------------------*/
pub fn arkPredict_MaximumOrder(ark_mem: &mut ARKodeMem, tau: f64, yguess: &mut NVector) -> i32 {
    /* verify that the interpolation structure is provided */
    if ark_mem.interp.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_NULL,
            line!(),
            "arkPredict_MaximumOrder",
            file!(),
            "ARKodeInterpMem structure is NULL",
        );
        return ARK_MEM_NULL;
    }

    /* call the interpolation module to do the work */
    crate::arkode_interp::arkInterpEvaluate(ark_mem, tau, 0, ARK_INTERP_MAX_DEGREE, yguess)
}

/*---------------------------------------------------------------
  arkPredict_VariableOrder

  This routine predicts the nonlinear implicit stage solution
  using the ARKODE interpolation module.  The degree of the
  interpolant is based on the level of extrapolation outside the
  preceding time step.
  ---------------------------------------------------------------*/
pub fn arkPredict_VariableOrder(ark_mem: &mut ARKodeMem, tau: f64, yguess: &mut NVector) -> i32 {
    let tau_tol: f64 = HALF;
    let tau_tol2: f64 = 0.75;

    /* verify that the interpolation structure is provided */
    if ark_mem.interp.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_NULL,
            line!(),
            "arkPredict_VariableOrder",
            file!(),
            "ARKodeInterpMem structure is NULL",
        );
        return ARK_MEM_NULL;
    }

    /* set the polynomial order based on tau input */
    let ord = if tau <= tau_tol {
        3
    } else if tau <= tau_tol2 {
        2
    } else {
        1
    };

    /* call the interpolation module to do the work */
    crate::arkode_interp::arkInterpEvaluate(ark_mem, tau, 0, ord, yguess)
}

/*---------------------------------------------------------------
  arkPredict_CutoffOrder

  This routine predicts the nonlinear implicit stage solution
  using the ARKODE interpolation module.  If the level of
  extrapolation is small enough, it uses the maximum degree
  polynomial available (stored in the interpolation module
  structure); otherwise it uses a linear polynomial.
  ---------------------------------------------------------------*/
pub fn arkPredict_CutoffOrder(ark_mem: &mut ARKodeMem, tau: f64, yguess: &mut NVector) -> i32 {
    let tau_tol: f64 = HALF;

    /* verify that the interpolation structure is provided */
    if ark_mem.interp.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_NULL,
            line!(),
            "arkPredict_CutoffOrder",
            file!(),
            "ARKodeInterpMem structure is NULL",
        );
        return ARK_MEM_NULL;
    }

    /* set the polynomial order based on tau input */
    let ord = if tau <= tau_tol { ARK_INTERP_MAX_DEGREE } else { 1 };

    /* call the interpolation module to do the work */
    crate::arkode_interp::arkInterpEvaluate(ark_mem, tau, 0, ord, yguess)
}

/*---------------------------------------------------------------
  arkPredict_Bootstrap

  This routine predicts the nonlinear implicit stage solution
  using a quadratic Hermite interpolating polynomial, based on
  the data {y_n, f(t_n,y_n), f(t_n+hj,z_j)}.

  Note: we assume that ftemp = f(t_n+hj,z_j) can be computed via
     N_VLinearCombination(nvec, cvals, Xvecs, ftemp),
  i.e. the inputs cvals[0:nvec-1] and Xvecs[0:nvec-1] may be
  combined to form f(t_n+hj,z_j).  Here the operand list arrives
  as owned coefficient/vector-reference slices assembled by the
  caller (Xvecs replaced by call-site assembly, Addendum C).
  ---------------------------------------------------------------*/
pub fn arkPredict_Bootstrap(
    ark_mem: &mut ARKodeMem,
    hj: f64,
    tau: f64,
    nvec: usize,
    cvals: &[f64],
    xvecs: &[&NVector],
    yguess: &mut NVector,
) -> i32 {
    /* verify that the interpolation structure is provided */
    if ark_mem.interp.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_NULL,
            line!(),
            "arkPredict_Bootstrap",
            file!(),
            "ARKodeInterpMem structure is NULL",
        );
        return ARK_MEM_NULL;
    }

    /* set coefficients for Hermite interpolant */
    let a0 = ONE;
    let a2 = tau * tau / TWO / hj;
    let a1 = tau - a2;

    /* set arrays for fused vector operation; shift inputs for
       f(t_n+hj,z_j) to end of queue */
    let mut cv: Vec<f64> = Vec::with_capacity(nvec + 2);
    let mut xr: Vec<&NVector> = Vec::with_capacity(nvec + 2);
    cv.push(a0);
    xr.push(&ark_mem.yn);
    cv.push(a1);
    xr.push(&ark_mem.fn_);
    for i in 0..nvec {
        cv.push(a2 * cvals[i]);
        xr.push(xvecs[i]);
    }

    /* call fused vector operation to compute prediction */
    N_VLinearCombination((nvec + 2) as i32, &cv, &xr, yguess);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkCheckConvergence

  This routine checks the return flag from the time-stepper module
  and handles solver convergence failures (see the C source for
  the full description).
  ---------------------------------------------------------------*/
pub fn arkCheckConvergence(ark_mem: &mut ARKodeMem, nflagPtr: &mut i32, ncfPtr: &mut i32) -> i32 {
    use crate::arkode_impl::{
        ARK_CONV_FAILURE, ARK_LSETUP_FAIL, ARK_LSOLVE_FAIL, ARK_MEM_NULL,
        ARK_REPTD_RHSFUNC_ERR, ARK_RETRY_STEP, ARK_RHSFUNC_FAIL, CONV_FAIL, MSG_ARKADAPT_NO_MEM,
        ONE, ONEPSM, PREDICT_AGAIN, PREV_CONV_FAIL, RHSFUNC_RECVR,
    };

    /* If nonlinear solver succeeded, return with ARK_SUCCESS */
    if *nflagPtr == ARK_SUCCESS {
        return ARK_SUCCESS;
    }
    /* Returns with an ARK_RETRY_STEP flag occur at a stage well before any
    algebraic solvers are involved. On the other hand, the
    arkCheckConvergence function handles the results from algebraic solvers,
    which never take place with an ARK_RETRY_STEP flag. Therefore, we
    immediately return from arkCheckConvergence, as it is irrelevant in the
    case of an ARK_RETRY_STEP */
    if *nflagPtr == ARK_RETRY_STEP {
        return ARK_RETRY_STEP;
    }

    /* The nonlinear soln. failed; increment ncfn */
    ark_mem.ncfn += 1;

    /* If fixed time stepping, then return with convergence failure */
    if ark_mem.fixedstep {
        return ARK_CONV_FAILURE;
    }

    /* Otherwise, access adaptivity structure */
    if ark_mem.hadapt_mem.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_NULL,
            line!(),
            "arkCheckConvergence",
            file!(),
            MSG_ARKADAPT_NO_MEM,
        );
        return ARK_MEM_NULL;
    }

    /* Return if lsetup, lsolve, or rhs failed unrecoverably */
    if *nflagPtr < 0 {
        if *nflagPtr == ARK_LSETUP_FAIL {
            return ARK_LSETUP_FAIL;
        } else if *nflagPtr == ARK_LSOLVE_FAIL {
            return ARK_LSOLVE_FAIL;
        } else if *nflagPtr == ARK_RHSFUNC_FAIL {
            return ARK_RHSFUNC_FAIL;
        } else {
            return crate::arkode_impl::ARK_NLS_OP_ERR;
        }
    }

    /* At this point, nflag = CONV_FAIL or RHSFUNC_RECVR; increment ncf */
    *ncfPtr += 1;
    ark_mem.hadapt_mem.as_mut().unwrap().etamax = ONE;

    /* If we had maxncf failures, or if |h| = hmin,
    return ARK_CONV_FAILURE or ARK_REPTD_RHSFUNC_ERR. */
    if (*ncfPtr == ark_mem.maxncf) || (SUNRabs(ark_mem.h) <= ark_mem.hmin * ONEPSM) {
        if *nflagPtr == CONV_FAIL {
            return ARK_CONV_FAILURE;
        }
        if *nflagPtr == RHSFUNC_RECVR {
            return ARK_REPTD_RHSFUNC_ERR;
        }
    }

    /* Reduce step size due to convergence failure */
    ark_mem.eta = ark_mem.hadapt_mem.as_ref().unwrap().etacf;

    /* Signal for Jacobian/preconditioner setup */
    *nflagPtr = PREV_CONV_FAIL;

    /* Return to reattempt the step */
    PREDICT_AGAIN
}

/*---------------------------------------------------------------
  arkCheckConstraints

  This routine determines if the constraints of the problem
  are satisfied by the proposed step

  Returns ARK_SUCCESS if successful, otherwise CONSTR_RECVR
  --------------------------------------------------------------*/
pub fn arkCheckConstraints(ark_mem: &mut ARKodeMem, constrfails: &mut i32, nflag: &mut i32) -> i32 {
    use crate::arkode_impl::{
        ARK_CONSTR_FAIL, CONSTR_RECVR, ONE, ONEPSM, PREV_CONV_FAIL, TENTH,
    };
    use crate::nvector_serial::{N_VConstrMask, N_VLinearSum, N_VMinQuotient};
    use crate::sundials_math::SUNMAX;

    /* Check constraints and get mask vector mm (tempv4) for where
    constraints failed */
    let constraints_passed = N_VConstrMask(
        ark_mem.constraints.as_ref().unwrap(),
        &ark_mem.ycur,
        &mut ark_mem.tempv4,
    );
    if constraints_passed {
        return ARK_SUCCESS;
    }

    /* Constraints not met */

    /* Update total fails and fails in current step */
    ark_mem.nconstrfails += 1;
    *constrfails += 1;

    /* Return with error if reached max fails in a step */
    if *constrfails == ark_mem.maxconstrfails {
        return ARK_CONSTR_FAIL;
    }

    /* Return with error if using fixed step sizes */
    if ark_mem.fixedstep {
        return ARK_CONSTR_FAIL;
    }

    /* Return with error if |h| == hmin */
    if SUNRabs(ark_mem.h) <= ark_mem.hmin * ONEPSM {
        return ARK_CONSTR_FAIL;
    }

    /* Reduce h by computing eta = h'/h */
    N_VLinearSum(ONE, &ark_mem.yn, -ONE, &ark_mem.ycur, &mut ark_mem.tempv3);
    /* N_VProd(mm, tmp, tmp): output aliases the second operand;
    z[i] = x[i]*y[i] with y == z (multiplication commutes bit-exactly) */
    for k in 0..ark_mem.tempv3.data.len() {
        ark_mem.tempv3.data[k] *= ark_mem.tempv4.data[k];
    }
    ark_mem.eta = 0.9 * N_VMinQuotient(&ark_mem.yn, &ark_mem.tempv3);
    ark_mem.eta = SUNMAX(ark_mem.eta, TENTH);

    /* Signal for Jacobian/preconditioner setup */
    *nflag = PREV_CONV_FAIL;

    /* Return to reattempt the step */
    CONSTR_RECVR
}

/*---------------------------------------------------------------
  arkCheckTemporalError

  This routine performs the local error test for the method (see
  the C source for the full description).
  --------------------------------------------------------------*/
pub fn arkCheckTemporalError(
    ark_mem: &mut ARKodeMem,
    nflagPtr: &mut i32,
    nefPtr: &mut i32,
    dsm: f64,
) -> i32 {
    use crate::arkode_impl::{
        ARK_ERR_FAILURE, ARK_MEM_NULL, MSG_ARKADAPT_NO_MEM, ONE, PREV_ERR_FAIL, TRY_AGAIN,
    };
    use crate::sundials_math::{SUNMAX, SUNMIN};

    /* Access hadapt_mem structure */
    if ark_mem.hadapt_mem.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_NULL,
            line!(),
            "arkCheckTemporalError",
            file!(),
            MSG_ARKADAPT_NO_MEM,
        );
        return ARK_MEM_NULL;
    }

    /* consider change of step size for next step attempt (may be
    larger/smaller than current step, depending on dsm) */
    let ttmp = if dsm <= ONE {
        ark_mem.tn + ark_mem.h
    } else {
        ark_mem.tn
    };
    let (h, hmin, hmax_inv) = (ark_mem.h, ark_mem.hmin, ark_mem.hmax_inv);
    let retval = crate::arkode_adapt::arkAdapt(ark_mem, ttmp, h, dsm);
    if retval != ARK_SUCCESS {
        return ARK_ERR_FAILURE;
    }

    /* if we've made it here then no nonrecoverable failures occurred;
    someone above has recommended an 'eta' value for the next step --
    enforce bounds on that value and set upcoming step size */
    ark_mem.eta = SUNMIN(ark_mem.eta, ark_mem.hadapt_mem.as_ref().unwrap().etamax);
    ark_mem.eta = SUNMAX(ark_mem.eta, hmin / SUNRabs(h));
    ark_mem.eta /= SUNMAX(ONE, SUNRabs(h) * hmax_inv * ark_mem.eta);

    /* If est. local error norm dsm passes test, return ARK_SUCCESS */
    if dsm <= ONE {
        return ARK_SUCCESS;
    }

    /* Test failed; increment counters, set nflag */
    *nefPtr += 1;
    ark_mem.netf += 1;
    *nflagPtr = PREV_ERR_FAIL;

    /* At maxnef failures, return ARK_ERR_FAILURE */
    if *nefPtr == ark_mem.maxnef {
        return ARK_ERR_FAILURE;
    }

    /* Set etamax=1 to prevent step size increase at end of this step */
    ark_mem.hadapt_mem.as_mut().unwrap().etamax = ONE;

    /* Enforce failure bounds on eta */
    if *nefPtr >= ark_mem.hadapt_mem.as_ref().unwrap().small_nef {
        ark_mem.eta = SUNMIN(ark_mem.eta, ark_mem.hadapt_mem.as_ref().unwrap().etamxf);
    }

    /* Enforce min/max step bounds once again due to adjustments above */
    ark_mem.eta = SUNMIN(ark_mem.eta, ark_mem.hadapt_mem.as_ref().unwrap().etamax);
    ark_mem.eta = SUNMAX(ark_mem.eta, hmin / SUNRabs(h));
    ark_mem.eta /= SUNMAX(ONE, SUNRabs(h) * hmax_inv * ark_mem.eta);

    TRY_AGAIN
}

/*---------------------------------------------------------------
  arkInitialSetup

  This routine performs all necessary items to prepare ARKODE for
  the first internal step after initialization, reinitialization,
  a reset() call, or a resize() call, including:
  - input consistency checks
  - (re)initializes the stepper
  - computes error and residual weights
  - (re)initialize the interpolation structure
  - checks for valid initial step input or estimates first step
  - checks for approach to tstop
  - checks for root near t0
  ---------------------------------------------------------------*/
pub fn arkInitialSetup(ark_mem: &mut ARKodeMem, tout: f64) -> i32 {
    use crate::arkode_impl::{
        ARK_ILL_INPUT, ARK_INTERP_LAGRANGE, ARK_INTERP_NONE, ARK_STEP_H0_FAIL,
        ARK_TOO_CLOSE, ARK_WF, FOUR, MSG_ARK_MISSING_FULLRHS, ONE, TWO,
    };
    use crate::nvector_serial::N_VConstrMask;
    use crate::sundials_math::SUNMAX;

    /* Is tout too close to tn? */
    let tdist = SUNRabs(tout - ark_mem.tcur);
    let tround = ark_mem.uround * SUNMAX(SUNRabs(ark_mem.tcur), SUNRabs(tout));

    if tdist == ZERO || tdist < TWO * tround {
        arkProcessError(
            Some(ark_mem),
            ARK_TOO_CLOSE,
            line!(),
            "arkInitialSetup",
            file!(),
            "tout too close to t0 to start integration.",
        );
        return ARK_TOO_CLOSE;
    }

    /* Check that user has supplied an initial step size if fixedstep mode is
    on */
    if ark_mem.fixedstep && ark_mem.hin == ZERO {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "arkInitialSetup",
            file!(),
            "Fixed step mode enabled, but no step size set",
        );
        return ARK_ILL_INPUT;
    }

    /* Optional N_Vector checks always pass in the serial build */

    /* Test input tstop for legality (correct direction of integration) */
    if ark_mem.tstopset {
        let htmp = if ark_mem.h == ZERO {
            tout - ark_mem.tcur
        } else {
            ark_mem.h
        };
        if (ark_mem.tstop - ark_mem.tcur) * htmp <= ZERO {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "arkInitialSetup",
                file!(),
                &format!(
                    "The value tstop = {} is behind current t = {} in the direction of integration.",
                    fmt_g(ark_mem.tstop, 0, 15),
                    fmt_g(ark_mem.tcur, 0, 15)
                ),
            );
            return ARK_ILL_INPUT;
        }
    }

    /* Check to see if y0 satisfies constraints */
    if let Some(constraints) = ark_mem.constraints.as_ref() {
        let con_ok = N_VConstrMask(constraints, &ark_mem.yn, &mut ark_mem.tempv1);
        if !con_ok {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "arkInitialSetup",
                file!(),
                "y0 fails to satisfy constraints.",
            );
            return ARK_ILL_INPUT;
        }
    }

    /* Load initial error weights */
    {
        let mut ewt = std::mem::replace(&mut ark_mem.ewt, NVector::new(0));
        let retval = ark_efun_apply_yn(ark_mem, &mut ewt);
        ark_mem.ewt = ewt;
        if retval != 0 {
            if ark_mem.itol == ARK_WF {
                arkProcessError(
                    Some(ark_mem),
                    ARK_ILL_INPUT,
                    line!(),
                    "arkInitialSetup",
                    file!(),
                    "The user-provide EwtSet function failed.",
                );
            } else {
                arkProcessError(
                    Some(ark_mem),
                    ARK_ILL_INPUT,
                    line!(),
                    "arkInitialSetup",
                    file!(),
                    "Initial ewt has component(s) equal to zero (illegal).",
                );
            }
            return ARK_ILL_INPUT;
        }
    }

    /* Set up the time stepper module if not done so already */
    if !ark_mem.preallocated {
        if ark_mem.step_init.is_none() {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "arkInitialSetup",
                file!(),
                "Time stepper module is missing",
            );
            return ARK_ILL_INPUT;
        }
        let step_init = ark_mem.step_init.unwrap();
        let retval = step_init(ark_mem, ark_mem.init_type);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                retval,
                line!(),
                "arkInitialSetup",
                file!(),
                "Error in initialization of time stepper module",
            );
            return retval;
        }
    }

    /* Load initial residual weights */
    if ark_mem.rwt_is_ewt {
        /* C updates the rwt pointer to ewt; readers dispatch on the flag */
    } else {
        let mut rwt = std::mem::replace(&mut ark_mem.rwt, NVector::new(0));
        let retval = ark_rfun_apply_yn(ark_mem, &mut rwt);
        ark_mem.rwt = rwt;
        if retval != 0 {
            if ark_mem.itol == ARK_WF {
                arkProcessError(
                    Some(ark_mem),
                    ARK_ILL_INPUT,
                    line!(),
                    "arkInitialSetup",
                    file!(),
                    "The user-provide RwtSet function failed.",
                );
            } else {
                arkProcessError(
                    Some(ark_mem),
                    ARK_ILL_INPUT,
                    line!(),
                    "arkInitialSetup",
                    file!(),
                    "Initial rwt has component(s) equal to zero (illegal).",
                );
            }
            return ARK_ILL_INPUT;
        }
    }

    /* Create default interpolation module (if needed) */
    if ark_mem.interp_type != ARK_INTERP_NONE && ark_mem.interp.is_none() {
        ark_mem.interp = if ark_mem.interp_type == ARK_INTERP_LAGRANGE {
            crate::arkode_interp::arkInterpCreate_Lagrange(ark_mem, ark_mem.interp_degree)
        } else {
            crate::arkode_interp::arkInterpCreate_Hermite(ark_mem, ark_mem.interp_degree)
        };
        if ark_mem.interp.is_none() {
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_FAIL,
                line!(),
                "arkInitialSetup",
                file!(),
                "Unable to allocate interpolation module",
            );
            return ARK_MEM_FAIL;
        }
    }

    /* Fill initial interpolation data (if needed) */
    if ark_mem.interp.is_some() {
        /* Stepper init may have limited the interpolation degree */
        if crate::arkode_interp::arkInterpSetDegree(ark_mem, ark_mem.interp_degree) != 0 {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "arkInitialSetup",
                file!(),
                "Unable to update interpolation polynomial degree",
            );
            return ARK_ILL_INPUT;
        }

        if crate::arkode_interp::arkInterpInit(ark_mem, ark_mem.tcur) != 0 {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "arkInitialSetup",
                file!(),
                "Unable to initialize interpolation module",
            );
            return ARK_ILL_INPUT;
        }
    }

    /* Check if the configuration requires interpolation */
    if ark_mem.root_mem.is_some() && ark_mem.interp.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "arkInitialSetup",
            file!(),
            "Rootfinding requires an interpolation module",
        );
        return ARK_ILL_INPUT;
    }

    if ark_mem.tstopinterp && ark_mem.interp.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "arkInitialSetup",
            file!(),
            "Stop time interpolation requires an interpolation module",
        );
        return ARK_ILL_INPUT;
    }

    /* Call stepper-provided initial step size estimation routine to fill
    ark_mem->hin, if applicable. */
    if ark_mem.h0u == ZERO && ark_mem.hin == ZERO && !ark_mem.fixedstep
        && ark_mem.step_H0.is_some()
    {
        let step_h0 = ark_mem.step_H0.unwrap();
        let mut hin = ark_mem.hin;
        let retval = step_h0(ark_mem, tout, &mut hin);
        ark_mem.hin = hin;
        if retval != 0 {
            arkProcessError(
                Some(ark_mem),
                ARK_STEP_H0_FAIL,
                line!(),
                "arkInitialSetup",
                file!(),
                "Failure in timestepping module h0 calculation",
            );
            return ARK_STEP_H0_FAIL;
        }
    }

    /* If fullrhs will be called (to estimate initial step, explicit
    steppers, Hermite interpolation module, and possibly (but not always)
    arkRootCheck1), then ensure that it is provided, and space is allocated
    for fn.  Otherwise, we should free ark_mem->fn if it is allocated. */
    if ark_mem.call_fullrhs || (ark_mem.h0u == ZERO && ark_mem.hin == ZERO)
        || ark_mem.root_mem.is_some()
    {
        if ark_mem.step_fullrhs.is_none() {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "arkInitialSetup",
                file!(),
                MSG_ARK_MISSING_FULLRHS,
            );
            return ARK_ILL_INPUT;
        }

        let yn_len = ark_mem.yn.data.len();
        let mut fn_ = std::mem::replace(&mut ark_mem.fn_, NVector::new(0));
        arkAllocVec(ark_mem, yn_len, &mut fn_);
        ark_mem.fn_ = fn_;
    } else if !ark_mem.fn_.data.is_empty() {
        let mut fn_ = std::mem::replace(&mut ark_mem.fn_, NVector::new(0));
        arkFreeVec(ark_mem, &mut fn_);
        ark_mem.fn_ = fn_;
    }

    /* initialization complete */
    ark_mem.initialized = true;

    /* Set initial step size */
    if ark_mem.h0u == ZERO {
        /* Check input h for validity */
        ark_mem.h = ark_mem.hin;
        if (ark_mem.h != ZERO) && ((tout - ark_mem.tcur) * ark_mem.h < ZERO) {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "arkInitialSetup",
                file!(),
                "h0 and tout - t0 inconsistent.",
            );
            return ARK_ILL_INPUT;
        }

        /* Estimate initial h if not set */
        if ark_mem.h == ZERO {
            /* If necessary, temporarily set h as it is used to compute the
            tolerance in a potential mass matrix solve when computing the
            full rhs */
            ark_mem.h = SUNRabs(tout - ark_mem.tcur);
            if ark_mem.h == ZERO {
                ark_mem.h = ONE;
            }

            /* Estimate the first step size */
            let mut tout_hin = tout;
            if ark_mem.tstopset && (tout - ark_mem.tcur) * (tout - ark_mem.tstop) > ZERO {
                tout_hin = ark_mem.tstop;
            }
            let hflag = arkHin(ark_mem, tout_hin);
            if hflag != ARK_SUCCESS {
                return arkHandleFailure(ark_mem, hflag);
            }

            /* Use first step growth factor for estimated h */
            let hadapt_mem = ark_mem.hadapt_mem.as_mut().unwrap();
            hadapt_mem.etamax = hadapt_mem.etamx1;
        } else if ark_mem.nst == 0 {
            /* Use first step growth factor for user defined h */
            let hadapt_mem = ark_mem.hadapt_mem.as_mut().unwrap();
            hadapt_mem.etamax = hadapt_mem.etamx1;
        } else {
            /* Use standard growth factor (e.g., for reset) */
            let hadapt_mem = ark_mem.hadapt_mem.as_mut().unwrap();
            hadapt_mem.etamax = hadapt_mem.growth;
        }

        /* Enforce step size bounds */
        let rh = SUNRabs(ark_mem.h) * ark_mem.hmax_inv;
        if rh > ONE {
            ark_mem.h /= rh;
        }
        if SUNRabs(ark_mem.h) < ark_mem.hmin {
            ark_mem.h *= ark_mem.hmin / SUNRabs(ark_mem.h);
        }

        /* Check for approach to tstop */
        if ark_mem.tstopset {
            if (ark_mem.tcur + ark_mem.h - ark_mem.tstop) * ark_mem.h > ZERO {
                ark_mem.h =
                    (ark_mem.tstop - ark_mem.tcur) * (ONE - FOUR * ark_mem.uround);
            }
        }

        /* Set initial time step factors */
        ark_mem.h0u = ark_mem.h;
        ark_mem.eta = ONE;
        ark_mem.hprime = ark_mem.h;
    } else {
        /* If next step would overtake tstop, adjust stepsize */
        if ark_mem.tstopset {
            if (ark_mem.tcur + ark_mem.hprime - ark_mem.tstop) * ark_mem.h > ZERO {
                ark_mem.hprime =
                    (ark_mem.tstop - ark_mem.tcur) * (ONE - FOUR * ark_mem.uround);
                ark_mem.eta = ark_mem.hprime / ark_mem.h;
            }
        }
    }

    /* Check for zeros of root function g at and near t0. */
    if ark_mem.root_mem.is_some() {
        if ark_mem.root_mem.as_ref().unwrap().nrtfn > 0 {
            let retval = crate::arkode_root::arkRootCheck1(ark_mem);
            if retval != ARK_SUCCESS {
                return retval;
            }
        }
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStopTests

  This routine performs relevant stopping tests:
  - check for root in last step
  - check if we passed tstop
  - check if we passed tout (NORMAL mode)
  - check if current tn was returned (ONE_STEP mode)
  - check if we are close to tstop (adjust step size if needed)

  Returns 1 (with *ier set) if the Evolve loop should return to
  the user, 0 to continue.
  ---------------------------------------------------------------*/
pub fn arkStopTests(
    ark_mem: &mut ARKodeMem,
    tout: f64,
    yout: &mut NVector,
    tret: &mut f64,
    itask: i32,
    ier: &mut i32,
) -> i32 {
    use crate::arkode_impl::{
        ARK_FULLRHS_END, ARK_ILL_INPUT, ARK_NORMAL, ARK_ONE_STEP, ARK_RHSFUNC_FAIL,
        ARK_ROOT_RETURN, ARK_RTFUNC_FAIL, ARK_TSTOP_RETURN, CLOSERT, FOUR, ONE, RTFOUND,
    };
    use crate::nvector_serial::N_VScale;

    /* Estimate an infinitesimal time interval to be used as a roundoff for
    time quantities (based on current time and step size) */
    let troundoff =
        FUZZ_FACTOR * ark_mem.uround * (SUNRabs(ark_mem.tcur) + SUNRabs(ark_mem.h));

    /* First, check for a root in the last step taken, other than the last
    root found, if any.  If itask = ARK_ONE_STEP and y(tn) was not returned
    because of an intervening root, return y(tn) now.     */
    if ark_mem.root_mem.is_some() {
        if ark_mem.root_mem.as_ref().unwrap().nrtfn > 0 {
            /* Shortcut to roots found in previous step */
            let irfndp = ark_mem.root_mem.as_ref().unwrap().irfnd;

            /* If the full RHS was not computed in the last call to
            arkCompleteStep and roots were found in the previous step, then
            compute the full rhs for possible use in arkRootCheck2 (not
            always necessary) */
            if !ark_mem.fn_is_current && irfndp != 0 {
                let retval = crate::arkode_impl::ark_step_fullrhs_yn_fn(
                    ark_mem,
                    ark_mem.tn,
                    ARK_FULLRHS_END,
                );
                if retval != 0 {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_RHSFUNC_FAIL,
                        line!(),
                        "arkStopTests",
                        file!(),
                        "The right-hand side routine failed in an unrecoverable manner.",
                    );
                    *ier = ARK_RHSFUNC_FAIL;
                    return 1;
                }
                ark_mem.fn_is_current = true;
            }

            let retval = crate::arkode_root::arkRootCheck2(ark_mem);

            if retval == CLOSERT {
                arkProcessError(
                    Some(ark_mem),
                    ARK_ILL_INPUT,
                    line!(),
                    "arkStopTests",
                    file!(),
                    &format!(
                        "Root found at and very near t = {}.",
                        fmt_g(ark_mem.root_mem.as_ref().unwrap().tlo, 0, 15)
                    ),
                );
                *ier = ARK_ILL_INPUT;
                return 1;
            } else if retval == ARK_RTFUNC_FAIL {
                arkProcessError(
                    Some(ark_mem),
                    ARK_RTFUNC_FAIL,
                    line!(),
                    "arkStopTests",
                    file!(),
                    &format!(
                        "At t = {}, the rootfinding routine failed in an unrecoverable manner.",
                        fmt_g(ark_mem.root_mem.as_ref().unwrap().tlo, 0, 15)
                    ),
                );
                *ier = ARK_RTFUNC_FAIL;
                return 1;
            } else if retval == RTFOUND {
                let tlo = ark_mem.root_mem.as_ref().unwrap().tlo;
                ark_mem.tretlast = tlo;
                *tret = tlo;
                /* C: yout aliases ycur, which holds the Check2 probe state */
                N_VScale(ONE, &ark_mem.ycur, yout);
                *ier = ARK_ROOT_RETURN;
                return 1;
            }

            /* If tn is distinct from tretlast (within roundoff), check
            remaining interval for roots */
            if SUNRabs(ark_mem.tcur - ark_mem.tretlast) > troundoff {
                let retval = crate::arkode_root::arkRootCheck3(ark_mem, tout, itask);

                if retval == ARK_SUCCESS {
                    /* no root found */
                    ark_mem.root_mem.as_mut().unwrap().irfnd = 0;
                    if (irfndp == 1) && (itask == ARK_ONE_STEP) {
                        ark_mem.tretlast = ark_mem.tcur;
                        *tret = ark_mem.tcur;
                        N_VScale(ONE, &ark_mem.yn, yout);
                        *ier = ARK_SUCCESS;
                        return 1;
                    }
                } else if retval == RTFOUND {
                    /* a new root was found */
                    ark_mem.root_mem.as_mut().unwrap().irfnd = 1;
                    let tlo = ark_mem.root_mem.as_ref().unwrap().tlo;
                    ark_mem.tretlast = tlo;
                    *tret = tlo;
                    /* C: yout aliases ycur = y(trout) from arkRootCheck3 */
                    N_VScale(ONE, &ark_mem.ycur, yout);
                    *ier = ARK_ROOT_RETURN;
                    return 1;
                } else if retval == ARK_RTFUNC_FAIL {
                    /* g failed */
                    arkProcessError(
                        Some(ark_mem),
                        ARK_RTFUNC_FAIL,
                        line!(),
                        "arkStopTests",
                        file!(),
                        &format!(
                            "At t = {}, the rootfinding routine failed in an unrecoverable manner.",
                            fmt_g(ark_mem.root_mem.as_ref().unwrap().tlo, 0, 15)
                        ),
                    );
                    *ier = ARK_RTFUNC_FAIL;
                    return 1;
                }
            }
        } /* end of root stop check */
    }

    /* Test for tn at tstop or near tstop */
    if ark_mem.tstopset {
        if SUNRabs(ark_mem.tcur - ark_mem.tstop) <= troundoff {
            /* Ensure tout >= tstop, otherwise check for tout return below */
            if (tout - ark_mem.tstop) * ark_mem.h >= ZERO
                || SUNRabs(tout - ark_mem.tstop) <= troundoff
            {
                if ark_mem.tstopinterp && ark_mem.interp.is_some() {
                    *ier = ARKodeGetDky(ark_mem, ark_mem.tstop, 0, yout);
                    if *ier != ARK_SUCCESS {
                        arkProcessError(
                            Some(ark_mem),
                            ARK_ILL_INPUT,
                            line!(),
                            "arkStopTests",
                            file!(),
                            &format!(
                                "The value tstop = {} is behind current t = {} in the direction of integration.",
                                fmt_g(ark_mem.tstop, 0, 15),
                                fmt_g(ark_mem.tcur, 0, 15)
                            ),
                        );
                        *ier = ARK_ILL_INPUT;
                        return 1;
                    }
                } else {
                    N_VScale(ONE, &ark_mem.yn, yout);
                }
                ark_mem.tretlast = ark_mem.tstop;
                *tret = ark_mem.tstop;
                ark_mem.tstopset = false;
                *ier = ARK_TSTOP_RETURN;
                return 1;
            }
        }
        /* If next step would overtake tstop, adjust stepsize */
        else if (ark_mem.tcur + ark_mem.hprime - ark_mem.tstop) * ark_mem.h > ZERO {
            ark_mem.hprime =
                (ark_mem.tstop - ark_mem.tcur) * (ONE - FOUR * ark_mem.uround);
            ark_mem.eta = ark_mem.hprime / ark_mem.h;
        }
    }

    /* In ARK_NORMAL mode, test if tout was reached */
    if (itask == ARK_NORMAL) && ((ark_mem.tcur - tout) * ark_mem.h >= ZERO) {
        if ark_mem.interp.is_some() {
            *ier = ARKodeGetDky(ark_mem, tout, 0, yout);
            if *ier != ARK_SUCCESS {
                arkProcessError(
                    Some(ark_mem),
                    ARK_ILL_INPUT,
                    line!(),
                    "arkStopTests",
                    file!(),
                    &format!(
                        "Trouble interpolating at tout = {}. tout too far back in direction of integration",
                        fmt_g(tout, 0, 15)
                    ),
                );
                *ier = ARK_ILL_INPUT;
                return 1;
            }
            ark_mem.tretlast = tout;
            *tret = tout;
        } else {
            N_VScale(ONE, &ark_mem.yn, yout);
            ark_mem.tretlast = ark_mem.tcur;
            *tret = ark_mem.tcur;
        }
        *ier = ARK_SUCCESS;
        return 1;
    }

    /* In ARK_ONE_STEP mode, test if tn was returned */
    if itask == ARK_ONE_STEP && SUNRabs(ark_mem.tcur - ark_mem.tretlast) > troundoff {
        ark_mem.tretlast = ark_mem.tcur;
        *tret = ark_mem.tcur;
        N_VScale(ONE, &ark_mem.yn, yout);
        *ier = ARK_SUCCESS;
        return 1;
    }

    0
}

/*---------------------------------------------------------------
  ARKodeEvolve:

  This routine is the main driver of ARKODE-based integrators.

  It integrates over a time interval defined by the user, by
  calling the time step module to do internal time steps.

  The first time that ARKodeEvolve is called for a successfully
  initialized problem, it computes a tentative initial step size.

  ARKodeEvolve supports two modes as specified by itask: ARK_NORMAL and
  ARK_ONE_STEP.  In the ARK_NORMAL mode, the solver steps until
  it reaches or passes tout and then interpolates to obtain
  y(tout).  In the ARK_ONE_STEP mode, it takes one internal step
  and returns.  The behavior of both modes can be over-ridden
  through user-specification of ark_tstop (through the
  ARKodeSetStopTime function), in which case if a solver step
  would pass tstop, the step is shortened so that it stops at
  exactly the specified stop time, and hence interpolation of
  y(tout) is not required.

  ycur/yout aliasing note (Addendum C.1): C sets ycur = yout for
  the whole call; here ycur owns storage and the return paths that
  rely on the alias copy ycur into yout explicitly.
  ---------------------------------------------------------------*/
pub fn ARKodeEvolve(
    ark_mem: &mut ARKodeMem,
    tout: f64,
    yout: &mut NVector,
    tret: &mut f64,
    itask: i32,
) -> i32 {
    use crate::arkode_impl::{
        ARK_ERR_FAILURE, ARK_ILL_INPUT, ARK_NORMAL, ARK_NO_MALLOC,
        ARK_ONE_STEP, ARK_PRESTEPFN_FAIL, ARK_RELAX_MEM_NULL, ARK_RETRY_STEP, ARK_ROOT_RETURN,
        ARK_RTFUNC_FAIL, ARK_TOO_MUCH_ACC, ARK_TOO_MUCH_WORK, ARK_TSTOP_RETURN, ARK_WARNING,
        ARK_WF, FIRST_CALL, FOUR, ONE, ONEPSM, RTFOUND, TWO,
    };
    use crate::nvector_serial::{N_VScale, N_VWrmsNorm};

    /* C leaves istate uninitialized (every break sets it) */
    #[allow(unused_assignments)]
    let mut istate = ARK_SUCCESS;

    /* Check and process inputs */

    /* Check if ark_mem was allocated */
    if !ark_mem.MallocDone {
        arkProcessError(
            Some(ark_mem),
            ARK_NO_MALLOC,
            line!(),
            "ARKodeEvolve",
            file!(),
            "Attempt to call before ARKODE initialized.",
        );
        return ARK_NO_MALLOC;
    }

    /* Check for yout != NULL: C sets ark_mem->ycur = yout; here ycur
    owns matching storage instead (see the aliasing note above) */
    if ark_mem.ycur.data.len() != ark_mem.yn.data.len() {
        ark_mem.ycur = NVector::new(ark_mem.yn.data.len());
    }

    /* Check for valid itask */
    if (itask != ARK_NORMAL) && (itask != ARK_ONE_STEP) {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "ARKodeEvolve",
            file!(),
            "Illegal value for itask.",
        );
        return ARK_ILL_INPUT;
    }

    /* perform first-step-specific initializations:
    - initialize tret values to initialization time
    - perform initial integrator setup  */
    if ark_mem.initsetup {
        ark_mem.tretlast = ark_mem.tcur;
        *tret = ark_mem.tcur;
        let retval = arkInitialSetup(ark_mem, tout);
        if retval != ARK_SUCCESS {
            return retval;
        }
    }

    /* perform stopping tests */
    if !ark_mem.initsetup {
        let mut retval = ARK_SUCCESS;
        if arkStopTests(ark_mem, tout, yout, tret, itask, &mut retval) != 0 {
            return retval;
        }
    }

    /* fill current independent variable (and optionally ycur with yn) */
    ark_mem.tcur = ark_mem.tn;
    if ark_mem.ensure_ycur {
        let mut ycur = std::mem::replace(&mut ark_mem.ycur, NVector::new(0));
        N_VScale(ONE, &ark_mem.yn, &mut ycur);
        ark_mem.ycur = ycur;
    }

    /*--------------------------------------------------
    Looping point for successful internal steps (see the C source
    for the full description of the loop stages)
    --------------------------------------------------*/
    let mut nstloc: i64 = 0;
    loop {
        ark_mem.next_h = ark_mem.h;

        /* Reset and check ewt and rwt */
        if !ark_mem.initsetup {
            let mut ewt = std::mem::replace(&mut ark_mem.ewt, NVector::new(0));
            let ewtset_ok = ark_efun_apply_yn(ark_mem, &mut ewt);
            ark_mem.ewt = ewt;
            if ewtset_ok != 0 {
                if ark_mem.itol == ARK_WF {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_ILL_INPUT,
                        line!(),
                        "ARKodeEvolve",
                        file!(),
                        &format!(
                            "At t = {}, the user-provide EwtSet function failed.",
                            fmt_g(ark_mem.tcur, 0, 15)
                        ),
                    );
                } else {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_ILL_INPUT,
                        line!(),
                        "ARKodeEvolve",
                        file!(),
                        &format!(
                            "At t = {}, a component of ewt has become <= 0.",
                            fmt_g(ark_mem.tcur, 0, 15)
                        ),
                    );
                }

                istate = ARK_ILL_INPUT;
                ark_mem.tretlast = ark_mem.tcur;
                *tret = ark_mem.tcur;
                N_VScale(ONE, &ark_mem.yn, yout);
                break;
            }

            if !ark_mem.rwt_is_ewt {
                let mut rwt = std::mem::replace(&mut ark_mem.rwt, NVector::new(0));
                let ewtset_ok = ark_rfun_apply_yn(ark_mem, &mut rwt);
                ark_mem.rwt = rwt;
                if ewtset_ok != 0 {
                    if ark_mem.itol == ARK_WF {
                        arkProcessError(
                            Some(ark_mem),
                            ARK_ILL_INPUT,
                            line!(),
                            "ARKodeEvolve",
                            file!(),
                            &format!(
                                "At t = {}, the user-provide RwtSet function failed.",
                                fmt_g(ark_mem.tcur, 0, 15)
                            ),
                        );
                    } else {
                        arkProcessError(
                            Some(ark_mem),
                            ARK_ILL_INPUT,
                            line!(),
                            "ARKodeEvolve",
                            file!(),
                            &format!(
                                "At t = {}, a component of rwt has become <= 0.",
                                fmt_g(ark_mem.tcur, 0, 15)
                            ),
                        );
                    }

                    istate = ARK_ILL_INPUT;
                    ark_mem.tretlast = ark_mem.tcur;
                    *tret = ark_mem.tcur;
                    N_VScale(ONE, &ark_mem.yn, yout);
                    break;
                }
            }
        }

        /* Check for too many steps */
        if (ark_mem.mxstep > 0) && (nstloc >= ark_mem.mxstep) {
            arkProcessError(
                Some(ark_mem),
                ARK_TOO_MUCH_WORK,
                line!(),
                "ARKodeEvolve",
                file!(),
                &format!(
                    "At t = {}, mxstep steps taken before reaching tout.",
                    fmt_g(ark_mem.tcur, 0, 15)
                ),
            );
            istate = ARK_TOO_MUCH_WORK;
            ark_mem.tretlast = ark_mem.tcur;
            *tret = ark_mem.tcur;
            N_VScale(ONE, &ark_mem.yn, yout);
            break;
        }

        /* Check for too much accuracy requested */
        let nrm = N_VWrmsNorm(&ark_mem.yn, &ark_mem.ewt);
        ark_mem.tolsf = ark_mem.uround * nrm;
        if ark_mem.tolsf > ONE && !ark_mem.fixedstep {
            arkProcessError(
                Some(ark_mem),
                ARK_TOO_MUCH_ACC,
                line!(),
                "ARKodeEvolve",
                file!(),
                &format!(
                    "At t = {}, too much accuracy requested.",
                    fmt_g(ark_mem.tcur, 0, 15)
                ),
            );
            istate = ARK_TOO_MUCH_ACC;
            ark_mem.tretlast = ark_mem.tcur;
            *tret = ark_mem.tcur;
            N_VScale(ONE, &ark_mem.yn, yout);
            ark_mem.tolsf *= TWO;
            break;
        } else {
            ark_mem.tolsf = ONE;
        }

        /* Check for h below roundoff level in tn */
        if ark_mem.tcur + ark_mem.h == ark_mem.tcur {
            ark_mem.nhnil += 1;
            if ark_mem.nhnil <= ark_mem.mxhnil {
                arkProcessError(
                    Some(ark_mem),
                    ARK_WARNING,
                    line!(),
                    "ARKodeEvolve",
                    file!(),
                    &format!(
                        "Internal t = {} and h = {} are such that t + h = t on the next step. The solver will continue anyway.",
                        fmt_g(ark_mem.tcur, 0, 15),
                        fmt_g(ark_mem.h, 0, 15)
                    ),
                );
            }
            if ark_mem.nhnil == ark_mem.mxhnil {
                arkProcessError(
                    Some(ark_mem),
                    ARK_WARNING,
                    line!(),
                    "ARKodeEvolve",
                    file!(),
                    "The above warning has been issued mxhnil times and will not be issued again for this problem.",
                );
            }
        }

        /* Update parameter for upcoming step size */
        if ark_mem.hprime != ark_mem.h {
            ark_mem.h *= ark_mem.eta;
            ark_mem.next_h = ark_mem.h;
        }
        if ark_mem.fixedstep {
            ark_mem.h = ark_mem.hin;
            ark_mem.next_h = ark_mem.h;

            /* patch for 'fixedstep' + 'tstop' use case:
            limit fixed step size if step would overtake tstop */
            if ark_mem.tstopset {
                if (ark_mem.tcur + ark_mem.h - ark_mem.tstop) * ark_mem.h > ZERO {
                    ark_mem.h =
                        (ark_mem.tstop - ark_mem.tcur) * (ONE - FOUR * ark_mem.uround);
                }
            }
        }

        /* Looping point for step attempts */
        let mut dsm = ZERO;
        let mut kflag = ARK_SUCCESS;
        let relax_fails = 0; /* arkRelax bookkeeping (module pending) */
        let mut nflag = FIRST_CALL;
        let mut attempts = 0;
        let _ = relax_fails;
        let mut ncf = 0;
        let mut nef = 0;
        let mut constrfails = 0;
        ark_mem.last_kflag = 0;
        loop {
            /* increment attempt counters
            Note: kflag can only equal ARK_RETRY_STEP if the stepper rejected
            the current step size before performing calculations. Thus, we
            do not include those when keeping track of step "attempts". */
            if kflag != ARK_RETRY_STEP {
                attempts += 1;
                ark_mem.nst_attempts += 1;
            }

            /* fill tcur with the last accepted step time */
            ark_mem.tcur = ark_mem.tn;

            /* call the user-supplied pre-step function (if it exists) */
            if let Some(pre_step_fn) = ark_mem.PreStepFn {
                let retval = if ark_mem.ensure_ycur {
                    pre_step_fn(
                        ark_mem.tcur,
                        &ark_mem.ycur,
                        ark_mem.nst,
                        attempts,
                        &mut ark_mem.user_data,
                    )
                } else {
                    pre_step_fn(
                        ark_mem.tcur,
                        &ark_mem.yn,
                        ark_mem.nst,
                        attempts,
                        &mut ark_mem.user_data,
                    )
                };
                if retval != 0 {
                    return ARK_PRESTEPFN_FAIL;
                }
            }

            /* Call time stepper module to attempt a step:
            0 => step completed successfully
            >0 => step encountered recoverable failure; reduce step if possible
            <0 => step encountered unrecoverable failure */
            let step = ark_mem.step.unwrap();
            kflag = step(ark_mem, &mut dsm, &mut nflag);
            if kflag < 0 {
                break;
            }

            /* handle solver convergence failures */
            kflag = arkCheckConvergence(ark_mem, &mut nflag, &mut ncf);

            if kflag < 0 {
                break;
            }

            /* Perform relaxation (arkode_relaxation.c not yet ported;
            relax_enabled can only be set by ARKodeSetRelaxFn, which is
            also pending — fail loudly if ever reached) */
            if ark_mem.relax_enabled && (kflag == ARK_SUCCESS) {
                kflag = ARK_RELAX_MEM_NULL;
                break;
            }

            /* perform constraint-handling (if selected, and if solver check
            passed) */
            if ark_mem.constraints.is_some() && (kflag == ARK_SUCCESS) {
                kflag = arkCheckConstraints(ark_mem, &mut constrfails, &mut nflag);

                if kflag < 0 {
                    break;
                }
            }

            /* when fixed time-stepping is enabled, 'success' == successful
            stage solves (checked in previous block), so just enforce no step
            size change */
            if ark_mem.fixedstep {
                ark_mem.eta = ONE;
                break;
            }

            /* check temporal error (if checks above passed) */
            if kflag == ARK_SUCCESS {
                kflag = arkCheckTemporalError(ark_mem, &mut nflag, &mut nef, dsm);

                if kflag < 0 {
                    break;
                }
            }

            /* if ignoring temporal error test result (XBraid) force step to
            pass */
            if ark_mem.force_pass {
                ark_mem.last_kflag = kflag;
                kflag = ARK_SUCCESS;
                break;
            }

            /* break attempt loop on successful step */
            if kflag == ARK_SUCCESS {
                break;
            }

            /* unsuccessful step, if |h| = hmin, return ARK_ERR_FAILURE */
            if SUNRabs(ark_mem.h) <= ark_mem.hmin * ONEPSM {
                return ARK_ERR_FAILURE;
            }

            /* update h, hprime and next_h for next iteration */
            ark_mem.h *= ark_mem.eta;
            ark_mem.hprime = ark_mem.h;
            ark_mem.next_h = ark_mem.h;

            /* reset tcur to last saved internal time before reattempting step
            (and optionally ycur to yn ) */
            ark_mem.tcur = ark_mem.tn;
            if ark_mem.ensure_ycur {
                let mut ycur = std::mem::replace(&mut ark_mem.ycur, NVector::new(0));
                N_VScale(ONE, &ark_mem.yn, &mut ycur);
                ark_mem.ycur = ycur;
            }
        } /* end looping for step attempts */

        /* If step attempt loop succeeded, complete step (update current
        time, solution, error stepsize history arrays; call user-supplied
        step postprocessing function) */
        if kflag == ARK_SUCCESS {
            kflag = arkCompleteStep(ark_mem, dsm);
        }

        /* If step attempt loop failed, process flag and return to user */
        if kflag != ARK_SUCCESS {
            istate = arkHandleFailure(ark_mem, kflag);
            ark_mem.tretlast = ark_mem.tcur;
            *tret = ark_mem.tcur;
            N_VScale(ONE, &ark_mem.yn, yout);
            break;
        }

        nstloc += 1;

        /* Check for root in last step taken. */
        if ark_mem.root_mem.is_some() {
            if ark_mem.root_mem.as_ref().unwrap().nrtfn > 0 {
                let retval = crate::arkode_root::arkRootCheck3(ark_mem, tout, itask);
                if retval == RTFOUND {
                    /* A new root was found */
                    ark_mem.root_mem.as_mut().unwrap().irfnd = 1;
                    istate = ARK_ROOT_RETURN;
                    let tlo = ark_mem.root_mem.as_ref().unwrap().tlo;
                    ark_mem.tretlast = tlo;
                    *tret = tlo;
                    /* C: yout aliases ycur = y(trout) from arkRootCheck3 */
                    N_VScale(ONE, &ark_mem.ycur, yout);
                    break;
                } else if retval == ARK_RTFUNC_FAIL {
                    /* g failed */
                    arkProcessError(
                        Some(ark_mem),
                        ARK_RTFUNC_FAIL,
                        line!(),
                        "ARKodeEvolve",
                        file!(),
                        &format!(
                            "At t = {}, the rootfinding routine failed in an unrecoverable manner.",
                            fmt_g(ark_mem.root_mem.as_ref().unwrap().tlo, 0, 15)
                        ),
                    );
                    istate = ARK_RTFUNC_FAIL;
                    break;
                }

                /* If we are at the end of the first step and we still have
                some event functions that are inactive, issue a warning as
                this may indicate a user error in the implementation of the
                root function. */
                if ark_mem.nst == 1 {
                    let rootmem = ark_mem.root_mem.as_ref().unwrap();
                    let mut inactive_roots = false;
                    for ir in 0..rootmem.nrtfn as usize {
                        if !rootmem.gactive[ir] {
                            inactive_roots = true;
                            break;
                        }
                    }
                    let mxgnull = rootmem.mxgnull;
                    if (mxgnull > 0) && inactive_roots {
                        arkProcessError(
                            Some(ark_mem),
                            ARK_WARNING,
                            line!(),
                            "ARKodeEvolve",
                            file!(),
                            "At the end of the first step, there are still some root functions identically 0. This warning will not be issued again.",
                        );
                    }
                }
            }
        }

        /* Check if tn is at tstop or near tstop */
        if ark_mem.tstopset {
            let troundoff =
                FUZZ_FACTOR * ark_mem.uround * (SUNRabs(ark_mem.tcur) + SUNRabs(ark_mem.h));

            if SUNRabs(ark_mem.tcur - ark_mem.tstop) <= troundoff {
                /* Ensure tout >= tstop, otherwise check for tout return
                below */
                if (tout - ark_mem.tstop) * ark_mem.h >= ZERO
                    || SUNRabs(tout - ark_mem.tstop) <= troundoff
                {
                    if ark_mem.tstopinterp && ark_mem.interp.is_some() {
                        let retval = ARKodeGetDky(ark_mem, ark_mem.tstop, 0, yout);
                        if retval != ARK_SUCCESS {
                            arkProcessError(
                                Some(ark_mem),
                                retval,
                                line!(),
                                "ARKodeEvolve",
                                file!(),
                                &format!(
                                    "At t = {}, interpolating the solution failed.",
                                    fmt_g(ark_mem.tstop, 0, 15)
                                ),
                            );
                            istate = retval;
                            break;
                        }
                    } else {
                        N_VScale(ONE, &ark_mem.yn, yout);
                    }
                    ark_mem.tretlast = ark_mem.tstop;
                    *tret = ark_mem.tstop;
                    ark_mem.tstopset = false;
                    istate = ARK_TSTOP_RETURN;
                    break;
                }
            }
            /* limit upcoming step if it will overcome tstop */
            else if (ark_mem.tcur + ark_mem.hprime - ark_mem.tstop) * ark_mem.h > ZERO {
                ark_mem.hprime =
                    (ark_mem.tstop - ark_mem.tcur) * (ONE - FOUR * ark_mem.uround);
                ark_mem.eta = ark_mem.hprime / ark_mem.h;
            }
        }

        /* In NORMAL mode, check if tout reached */
        if (itask == ARK_NORMAL) && (ark_mem.tcur - tout) * ark_mem.h >= ZERO {
            if ark_mem.interp.is_some() {
                let retval = ARKodeGetDky(ark_mem, tout, 0, yout);
                if retval != ARK_SUCCESS {
                    arkProcessError(
                        Some(ark_mem),
                        retval,
                        line!(),
                        "ARKodeEvolve",
                        file!(),
                        &format!(
                            "At t = {}, interpolating the solution failed.",
                            fmt_g(tout, 0, 15)
                        ),
                    );
                    istate = retval;
                    break;
                }
                ark_mem.tretlast = tout;
                *tret = tout;
            } else {
                N_VScale(ONE, &ark_mem.yn, yout);
                ark_mem.tretlast = ark_mem.tcur;
                *tret = ark_mem.tcur;
            }
            ark_mem.next_h = ark_mem.hprime;
            istate = ARK_SUCCESS;
            break;
        }

        /* In ONE_STEP mode, exit loop (arkCompleteStep already copied ycur
        to yn; C relies on ycur being an alias of yout) */
        if itask == ARK_ONE_STEP {
            istate = ARK_SUCCESS;
            ark_mem.tretlast = ark_mem.tcur;
            *tret = ark_mem.tcur;
            ark_mem.next_h = ark_mem.hprime;
            N_VScale(ONE, &ark_mem.ycur, yout);
            break;
        }
    } /* end looping for internal steps */

    istate
}

/*---------------------------------------------------------------
  ARKodeReset:

  This routine resets an ARKode module to solve the same
  problem from the given time with the input state (all counter
  values are retained).
  ---------------------------------------------------------------*/
pub fn ARKodeReset(ark_mem: &mut ARKodeMem, tR: f64, yR: &NVector) -> i32 {
    /* Reset main ARKODE infrastructure */
    let retval = arkInit(ark_mem, tR, yR, crate::arkode_impl::RESET_INIT);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            retval,
            line!(),
            "ARKodeReset",
            file!(),
            "ARKode reset failure",
        );
        return retval;
    }

    /* Call stepper routine to perform remaining reset operations (if provided) */
    if let Some(step_reset) = ark_mem.step_reset {
        return step_reset(ark_mem, tR, yR);
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSStolerances, ARKodeSVtolerances, ARKodeWFtolerances:

  These functions specify the integration tolerances. One of them
  SHOULD be called before the first call to ARKodeEvolve; otherwise
  default values of reltol=1e-4 and abstol=1e-9 will be used, which
  may be entirely incorrect for a specific problem.
  ---------------------------------------------------------------*/
pub fn ARKodeSStolerances(ark_mem: &mut ARKodeMem, reltol: f64, abstol: f64) -> i32 {
    use crate::arkode_impl::{ARK_ILL_INPUT, ARK_NO_MALLOC, ARK_SS};

    /* Check inputs */
    if !ark_mem.MallocDone {
        arkProcessError(
            Some(ark_mem),
            ARK_NO_MALLOC,
            line!(),
            "ARKodeSStolerances",
            file!(),
            "Attempt to call before ARKODE initialized.",
        );
        return ARK_NO_MALLOC;
    }
    if reltol < ZERO {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "ARKodeSStolerances",
            file!(),
            "reltol < 0 illegal.",
        );
        return ARK_ILL_INPUT;
    }
    if abstol < ZERO {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "ARKodeSStolerances",
            file!(),
            "abstol has negative component(s) (illegal).",
        );
        return ARK_ILL_INPUT;
    }

    /* (N_VAddConst is always provided by the serial vector) */

    /* Set flag indicating whether abstol == 0 */
    ark_mem.atolmin0 = abstol == ZERO;

    /* Copy tolerances into memory */
    ark_mem.reltol = reltol;
    ark_mem.Sabstol = abstol;
    ark_mem.itol = ARK_SS;

    /* enforce use of arkEwtSetSS (internal dispatch, Addendum C.1) */
    ark_mem.user_efun = false;
    ark_mem.efun = None;
    ark_mem.e_data = None;

    ARK_SUCCESS
}

pub fn ARKodeSVtolerances(ark_mem: &mut ARKodeMem, reltol: f64, abstol: &NVector) -> i32 {
    use crate::arkode_impl::{ARK_ILL_INPUT, ARK_NO_MALLOC, ARK_SV};
    use crate::nvector_serial::{N_VMin, N_VScale};

    /* Check inputs */
    if !ark_mem.MallocDone {
        arkProcessError(
            Some(ark_mem),
            ARK_NO_MALLOC,
            line!(),
            "ARKodeSVtolerances",
            file!(),
            "Attempt to call before ARKODE initialized.",
        );
        return ARK_NO_MALLOC;
    }
    if reltol < ZERO {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "ARKodeSVtolerances",
            file!(),
            "reltol < 0 illegal.",
        );
        return ARK_ILL_INPUT;
    }
    let abstolmin = N_VMin(abstol);
    if abstolmin < ZERO {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "ARKodeSVtolerances",
            file!(),
            "abstol has negative component(s) (illegal).",
        );
        return ARK_ILL_INPUT;
    }

    /* Set flag indicating whether min(abstol) == 0 */
    ark_mem.atolmin0 = abstolmin == ZERO;

    /* Copy tolerances into memory */
    if !ark_mem.VabstolMallocDone {
        /* arkAllocVec(ark_mem, ewt, &Vabstol) */
        ark_mem.Vabstol = Some(NVector::new(ark_mem.ewt.data.len()));
        ark_mem.lrw += ark_mem.lrw1;
        ark_mem.liw += ark_mem.liw1;
        ark_mem.VabstolMallocDone = true;
    }
    N_VScale(1.0, abstol, ark_mem.Vabstol.as_mut().unwrap());
    ark_mem.reltol = reltol;
    ark_mem.itol = ARK_SV;

    /* enforce use of arkEwtSetSV (internal dispatch, Addendum C.1) */
    ark_mem.user_efun = false;
    ark_mem.efun = None;
    ark_mem.e_data = None;

    ARK_SUCCESS
}

pub fn ARKodeWFtolerances(ark_mem: &mut ARKodeMem, efun: crate::arkode_impl::ARKEwtFn) -> i32 {
    use crate::arkode_impl::{ARK_NO_MALLOC, ARK_WF};

    if !ark_mem.MallocDone {
        arkProcessError(
            Some(ark_mem),
            ARK_NO_MALLOC,
            line!(),
            "ARKodeWFtolerances",
            file!(),
            "Attempt to call before ARKODE initialized.",
        );
        return ARK_NO_MALLOC;
    }

    /* Copy tolerance data into memory (C: e_data = user_data; the user
    efun receives ark_mem.user_data through the dispatch helper) */
    ark_mem.itol = ARK_WF;
    ark_mem.user_efun = true;
    ark_mem.efun = Some(efun);
    ark_mem.e_data = None;

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeFree:

  This routine frees the ARKODE infrastructure memory (all the
  C free() calls become drops).
  ---------------------------------------------------------------*/
pub fn ARKodeFree(arkode_mem: &mut Option<Box<ARKodeMem>>) {
    let mut ark_mem = match arkode_mem.take() {
        None => return,
        Some(m) => m,
    };

    /* free the time-stepper module memory (if provided) */
    if let Some(step_free) = ark_mem.step_free {
        step_free(&mut ark_mem);
    }

    /* free vector storage */
    arkFreeVectors(&mut ark_mem);

    /* free the time step adaptivity module (the controller and the
    structure drop with it) */
    ark_mem.hadapt_mem = None;

    /* free the interpolation module */
    crate::arkode_interp::arkInterpFree(&mut ark_mem);

    /* free the root-finding module */
    if ark_mem.root_mem.is_some() {
        let _ = crate::arkode_root::arkRootFree(&mut ark_mem);
    }

    /* free the relaxation module, constraints, step memory: drop */
    drop(ark_mem);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arkode_impl::{ARK_ILL_INPUT, FIRST_INIT};
    use crate::sundials_context::SUNContext_Create;

    fn stub_init(_a: &mut ARKodeMem, _t: i32) -> i32 {
        ARK_SUCCESS
    }
    fn stub_step(_a: &mut ARKodeMem, _d: &mut f64, _n: &mut i32) -> i32 {
        ARK_SUCCESS
    }

    #[test]
    fn init_rejects_missing_stepper_and_allocates() {
        let ctx = SUNContext_Create();
        let mut ark_mem = arkCreate(&ctx);
        let y0 = {
            let mut y = NVector::new(3);
            y.data.fill(2.0);
            y
        };

        /* no stepper attached: ARK_ILL_INPUT */
        assert_eq!(arkInit(&mut ark_mem, 0.5, &y0, FIRST_INIT), ARK_ILL_INPUT);

        /* attach a minimal stepper */
        ark_mem.step_init = Some(stub_init);
        ark_mem.step = Some(stub_step);
        ark_mem.step_mem = Some(Box::new(0i32));

        let lrw_before = ark_mem.lrw;
        assert_eq!(arkInit(&mut ark_mem, 0.5, &y0, FIRST_INIT), ARK_SUCCESS);
        assert!(ark_mem.MallocDone);
        assert_eq!(ark_mem.lrw1, 3);
        assert_eq!(ark_mem.liw1, 1);
        /* ewt + yn + tempv1..4 = 6 vectors (rwt aliases ewt) */
        assert_eq!(ark_mem.lrw, lrw_before + 6 * 3);
        assert_eq!(ark_mem.ewt.data.len(), 3);
        assert_eq!(ark_mem.tempv4.data.len(), 3);
        assert!(ark_mem.rwt.data.is_empty()); /* rwt_is_ewt */
        assert_eq!(ark_mem.yn.data, vec![2.0; 3]);
        assert_eq!((ark_mem.tcur, ark_mem.tn), (0.5, 0.5));
        assert_eq!(ark_mem.tolsf, 1.0);
        assert!(ark_mem.initsetup);
        assert!(ark_mem.firststage);

        /* free restores the workspace bookkeeping */
        arkFreeVectors(&mut ark_mem);
        assert_eq!(ark_mem.lrw, lrw_before);
    }

    fn const_rhs(_a: &mut ARKodeMem, _t: f64, _y: &NVector, f: &mut NVector, _m: i32) -> i32 {
        f.data[0] = 1.0;
        0
    }

    #[test]
    fn hin_constant_rhs_matches_algorithm() {
        use crate::arkode_impl::{H0_BIAS, H0_LBFACTOR, H0_UBFACTOR};
        use crate::sundials_math::SUNRsqrt;
        let ctx = SUNContext_Create();
        let mut ark_mem = arkCreate(&ctx);
        ark_mem.step_init = Some(stub_init);
        ark_mem.step = Some(stub_step);
        ark_mem.step_mem = Some(Box::new(0i32));
        ark_mem.step_fullrhs = Some(const_rhs);
        let mut y0 = NVector::new(1);
        y0.data[0] = 1.0;
        assert_eq!(arkInit(&mut ark_mem, 0.0, &y0, FIRST_INIT), ARK_SUCCESS);
        /* fn and ycur as arkInitialSetup / ARKodeEvolve would provide them */
        ark_mem.fn_ = NVector::new(1);
        ark_mem.ycur = NVector::new(1);
        /* ewt as arkInitialSetup would set it */
        let mut ewt = std::mem::replace(&mut ark_mem.ewt, NVector::new(0));
        assert_eq!(ark_efun_apply_yn(&mut ark_mem, &mut ewt), 0);
        ark_mem.ewt = ewt;

        assert_eq!(arkHin(&mut ark_mem, 1.0), ARK_SUCCESS);

        /* y' = 1 (constant): yddnrm = 0 every pass, so
           hg1   = sqrt(hlb*hub),
           hnew1 = sqrt(hg1*hub),
           pass 2 has hrat > 2 with count1 > 1 -> fall back to hnew1,
           h     = H0_BIAS * hnew1 (within [hlb, hub]).
           hub = H0_UBFACTOR * tdist = 0.1 here (the |y|/|y'| bound is
           larger). */
        let tround = ark_mem.uround * 1.0;
        let hlb = H0_LBFACTOR * tround;
        let hub = H0_UBFACTOR * 1.0;
        let expected = H0_BIAS * SUNRsqrt(SUNRsqrt(hlb * hub) * hub);
        assert_eq!(ark_mem.h, expected);

        /* integrating backwards flips the sign */
        ark_mem.fn_is_current = true;
        assert_eq!(arkHin(&mut ark_mem, -1.0), ARK_SUCCESS);
        assert!(ark_mem.h < 0.0);
    }

    /* ------------- ARKodeEvolve end-to-end machinery tests -------------
    A manufactured stepper that advances y = t^2 exactly: each step fills
    ycur with the exact solution at tn + h and reports dsm = 0.5 (always
    passes the error test). This exercises InitialSetup, the internal
    step loop, ewt refresh, CompleteStep, the Hermite interpolation
    chain, GetDky output interpolation, tstop and root returns. */

    fn t2_rhs(_a: &mut ARKodeMem, t: f64, _y: &NVector, f: &mut NVector, _m: i32) -> i32 {
        f.data[0] = 2.0 * t;
        0
    }

    fn t2_step(a: &mut ARKodeMem, dsm: &mut f64, nflag: &mut i32) -> i32 {
        let t1 = a.tcur + a.h;
        a.ycur.data[0] = t1 * t1;
        *dsm = 0.5;
        *nflag = ARK_SUCCESS;
        ARK_SUCCESS
    }

    fn t2_mem() -> Box<ARKodeMem> {
        let ctx = SUNContext_Create();
        let mut ark_mem = arkCreate(&ctx);
        ark_mem.step_init = Some(stub_init);
        ark_mem.step = Some(t2_step);
        ark_mem.step_fullrhs = Some(t2_rhs);
        ark_mem.step_mem = Some(Box::new(0i32));
        let mut y0 = NVector::new(1);
        y0.data[0] = 1.0; /* y(1) = 1 */
        assert_eq!(arkInit(&mut ark_mem, 1.0, &y0, FIRST_INIT), ARK_SUCCESS);
        ark_mem.hin = 0.1; /* user-supplied initial step */
        ark_mem
    }

    #[test]
    fn evolve_normal_mode_interpolates_tout() {
        use crate::arkode_impl::ARK_NORMAL;
        let mut ark_mem = t2_mem();
        let mut yout = NVector::new(1);
        let mut tret = 0.0;
        let istate = ARKodeEvolve(&mut ark_mem, 2.0, &mut yout, &mut tret, ARK_NORMAL);
        assert_eq!(istate, ARK_SUCCESS);
        assert_eq!(tret, 2.0);
        assert!((yout.data[0] - 4.0).abs() < 1e-9, "yout = {}", yout.data[0]);
        assert!(ark_mem.nst >= 10, "nst = {}", ark_mem.nst);
        /* continuation call also works (stop tests path) */
        let istate = ARKodeEvolve(&mut ark_mem, 3.0, &mut yout, &mut tret, ARK_NORMAL);
        assert_eq!(istate, ARK_SUCCESS);
        assert_eq!(tret, 3.0);
        assert!((yout.data[0] - 9.0).abs() < 1e-8, "yout = {}", yout.data[0]);
    }

    #[test]
    fn evolve_one_step_and_fixedstep_modes() {
        use crate::arkode_impl::ARK_ONE_STEP;
        let mut ark_mem = t2_mem();
        ark_mem.fixedstep = true;
        ark_mem.hin = 0.25;
        let mut yout = NVector::new(1);
        let mut tret = 0.0;
        let istate = ARKodeEvolve(&mut ark_mem, 2.0, &mut yout, &mut tret, ARK_ONE_STEP);
        assert_eq!(istate, ARK_SUCCESS);
        assert!((tret - 1.25).abs() < 1e-12, "tret = {}", tret);
        assert!((yout.data[0] - 1.25 * 1.25).abs() < 1e-12);
        assert_eq!(ark_mem.nst, 1);
        /* second single step */
        let istate = ARKodeEvolve(&mut ark_mem, 2.0, &mut yout, &mut tret, ARK_ONE_STEP);
        assert_eq!(istate, ARK_SUCCESS);
        assert!((tret - 1.5).abs() < 1e-12, "tret = {}", tret);
        assert_eq!(ark_mem.nst, 2);
    }

    #[test]
    fn evolve_tstop_return() {
        use crate::arkode_impl::{ARK_NORMAL, ARK_TSTOP_RETURN};
        let mut ark_mem = t2_mem();
        ark_mem.tstopset = true;
        ark_mem.tstop = 1.5;
        let mut yout = NVector::new(1);
        let mut tret = 0.0;
        let istate = ARKodeEvolve(&mut ark_mem, 2.0, &mut yout, &mut tret, ARK_NORMAL);
        assert_eq!(istate, ARK_TSTOP_RETURN);
        assert_eq!(tret, 1.5);
        /* tstop return copies yn = y(tstop) (steps were limited to end
        exactly at tstop) */
        assert!((yout.data[0] - 2.25).abs() < 1e-9, "yout = {}", yout.data[0]);
    }

    #[test]
    fn evolve_root_return() {
        use crate::arkode_impl::{ARK_NORMAL, ARK_ROOT_RETURN};
        use crate::arkode_root::ARKodeRootInit;

        fn g_at_1p5(_t: f64, y: &NVector, gout: &mut [f64], _ud: &mut UserData) -> i32 {
            gout[0] = y.data[0] - 2.25; /* zero at t = 1.5 on y = t^2 */
            0
        }

        let mut ark_mem = t2_mem();
        assert_eq!(ARKodeRootInit(&mut ark_mem, 1, Some(g_at_1p5)), ARK_SUCCESS);
        let mut yout = NVector::new(1);
        let mut tret = 0.0;
        let istate = ARKodeEvolve(&mut ark_mem, 2.0, &mut yout, &mut tret, ARK_NORMAL);
        assert_eq!(istate, ARK_ROOT_RETURN);
        assert!((tret - 1.5).abs() < 1e-6, "tret = {}", tret);
        assert!((yout.data[0] - 2.25).abs() < 1e-6, "yout = {}", yout.data[0]);
        assert_eq!(ark_mem.root_mem.as_ref().unwrap().iroots[0], 1);

        /* continuing past the root reaches tout */
        let istate = ARKodeEvolve(&mut ark_mem, 2.0, &mut yout, &mut tret, ARK_NORMAL);
        assert_eq!(istate, ARK_SUCCESS);
        assert_eq!(tret, 2.0);
        assert!((yout.data[0] - 4.0).abs() < 1e-8, "yout = {}", yout.data[0]);
    }
}

/*---------------------------------------------------------------
  arkAllocVecArray / arkFreeVecArray:

  Allocate or free a vector array (C keeps the workspace pointers
  lrw/liw separate from ark_mem exactly as here).
  ---------------------------------------------------------------*/
#[allow(clippy::too_many_arguments)]
pub fn arkAllocVecArray(
    count: i32,
    tmpl_len: usize,
    v: &mut Vec<NVector>,
    lrw1: i64,
    lrw: &mut i64,
    liw1: i64,
    liw: &mut i64,
) -> bool {
    /* allocate the new vector array if necessary */
    if v.is_empty() {
        *v = (0..count).map(|_| NVector::new(tmpl_len)).collect();
        *lrw += count as i64 * lrw1;
        *liw += count as i64 * liw1;
    }
    true
}

pub fn arkFreeVecArray(
    count: i32,
    v: &mut Vec<NVector>,
    lrw1: i64,
    lrw: &mut i64,
    liw1: i64,
    liw: &mut i64,
) {
    if !v.is_empty() {
        *v = Vec::new();
        *lrw -= count as i64 * lrw1;
        *liw -= count as i64 * liw1;
    }
}
