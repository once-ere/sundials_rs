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
}
