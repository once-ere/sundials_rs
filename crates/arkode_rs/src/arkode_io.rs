/* -----------------------------------------------------------------
 * Translation of sundials-7.7.0/src/arkode/arkode_io.c — PART I:
 * ARKodeSetDefaults. The remaining Set / Get families follow in
 * later parts.
 *
 * The internal weight functions in C are installed as fn pointers
 * with e_data/r_data = ark_mem (self-referential); per the cvode
 * donor convention, efun = None means "internal dispatch" (the ewt
 * call sites branch on user_efun and call arkEwtSetSS/SV or
 * arkRwtSet directly), so SetDefaults leaves efun/rfun None here.
 * -----------------------------------------------------------------*/

use crate::arkode_adapt_impl::{
    ADJUST, CFLFAC, ETACF, ETAMIN, ETAMX1, ETAMXF, GROWTH, HFIXED_LB, HFIXED_UB, PQ, SAFETY,
    SMALL_NEF,
};
use crate::arkode_impl::{
    arkProcessError, ARKodeMem, ARK_SS, ARK_SUCCESS, MAXCONSTRFAILS, MAXNCF, MAXNEF, MXHNIL,
    MXSTEP_DEFAULT, ZERO,
};

/*---------------------------------------------------------------
  ARKodeSetDefaults:

  Resets all optional inputs to ARKODE default values.  Does not
  change problem-defining function pointers fe and fi or
  user_data pointer.  Also leaves alone any data
  structures/options related to root-finding (those can be reset
  using ARKodeRootInit) or post-processing a step (ProcessStep).
  ---------------------------------------------------------------*/
pub fn ARKodeSetDefaults(ark_mem: &mut ARKodeMem) -> i32 {
    /* Set default values for integrator optional inputs */
    ark_mem.use_compensated_sums = false;
    ark_mem.fixedstep = false; /* default to use adaptive steps */
    ark_mem.reltol = 1.0e-4; /* relative tolerance */
    ark_mem.itol = ARK_SS; /* scalar-scalar solution tolerances */
    ark_mem.ritol = ARK_SS; /* scalar-scalar residual tolerances */
    ark_mem.Sabstol = 1.0e-9; /* solution absolute tolerance */
    ark_mem.atolmin0 = false; /* min(abstol) > 0 */
    ark_mem.SRabstol = 1.0e-9; /* residual absolute tolerance */
    ark_mem.Ratolmin0 = false; /* min(Rabstol) > 0 */
    ark_mem.user_efun = false; /* no user-supplied ewt function */
    ark_mem.efun = None; /* internal arkEwtSetSS dispatch is used */
    ark_mem.e_data = None; /* (C: e_data = ark_mem) */
    ark_mem.user_rfun = false; /* no user-supplied rwt function */
    ark_mem.rfun = None; /* internal arkRwtSet dispatch is used */
    ark_mem.r_data = None; /* (C: r_data = ark_mem) */
    ark_mem.mxstep = MXSTEP_DEFAULT; /* max number of steps */
    ark_mem.mxhnil = MXHNIL; /* max warns of t+h==t */
    ark_mem.maxnef = MAXNEF; /* max error test fails */
    ark_mem.maxncf = MAXNCF; /* max convergence fails */
    ark_mem.maxconstrfails = MAXCONSTRFAILS; /* max number of constraint fails */
    ark_mem.preallocated = false; /* data was not preallocated */
    ark_mem.hin = ZERO; /* determine initial step on-the-fly */
    ark_mem.hmin = ZERO; /* no minimum step size */
    ark_mem.hmax_inv = ZERO; /* no maximum step size */
    ark_mem.tstopset = false; /* no stop time set */
    ark_mem.tstopinterp = false; /* copy at stop time */
    ark_mem.tstop = ZERO; /* no fixed stop time */
    {
        let hadapt_mem = ark_mem.hadapt_mem.as_mut().unwrap();
        hadapt_mem.etamx1 = ETAMX1; /* max change on first step */
        hadapt_mem.etamxf = ETAMXF; /* max change on error-failed step */
        hadapt_mem.etamin = ETAMIN; /* min bound on time step reduction */
        hadapt_mem.small_nef = SMALL_NEF; /* num error fails before ETAMXF enforced */
        hadapt_mem.etacf = ETACF; /* max change on convergence failure */
        hadapt_mem.cfl = CFLFAC; /* explicit stability factor */
        hadapt_mem.safety = SAFETY; /* step adaptivity safety factor  */
        hadapt_mem.growth = GROWTH; /* step adaptivity growth factor */
        hadapt_mem.lbound = HFIXED_LB; /* step adaptivity no-change lower bound */
        hadapt_mem.ubound = HFIXED_UB; /* step adaptivity no-change upper bound */
        hadapt_mem.expstab = None; /* no explicit stability fn */
        hadapt_mem.estab_data = None; /* no explicit stability fn data */
        hadapt_mem.pq = PQ; /* embedding order */
        hadapt_mem.p = 0; /* no default embedding order */
        hadapt_mem.q = 0; /* no default method order */
        hadapt_mem.adjust = ADJUST; /* controller order adjustment */
    }

    /* Set stepper defaults (if provided) */
    if let Some(step_setdefaults) = ark_mem.step_setdefaults {
        let retval = step_setdefaults(ark_mem);
        if retval != ARK_SUCCESS {
            return retval;
        }
    }

    ARK_SUCCESS
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arkode::{arkCreate, arkEwtSetSS, arkEwtSetSV, arkEwtSetSmallReal};
    use crate::arkode_adapt_impl::{ARK_ADAPT_LIW, ARK_ADAPT_LRW};
    use crate::arkode_impl::{ARK_INTERP_HERMITE, ARK_INTERP_MAX_DEGREE, FIRST_INIT};
    use crate::nvector_serial::NVector;
    use crate::sundials_context::SUNContext_Create;
    use crate::sundials_types::SUN_SMALL_REAL;

    #[test]
    fn create_sets_c_defaults() {
        let ctx = SUNContext_Create();
        let ark_mem = arkCreate(&ctx);

        /* arkCreate-set fields */
        assert_eq!(ark_mem.lrw, 18 + ARK_ADAPT_LRW);
        assert_eq!(ark_mem.liw, 53 + ARK_ADAPT_LIW);
        assert!(ark_mem.hadapt_mem.is_some());
        assert!(ark_mem.interp.is_none());
        assert_eq!(ark_mem.interp_type, ARK_INTERP_HERMITE);
        assert_eq!(ark_mem.interp_degree, ARK_INTERP_MAX_DEGREE);
        assert!(ark_mem.rwt_is_ewt);
        assert!(ark_mem.initsetup);
        assert_eq!(ark_mem.init_type, FIRST_INIT);
        assert!(ark_mem.firststage);
        assert!(!ark_mem.initialized);

        /* ARKodeSetDefaults-set fields */
        assert_eq!(ark_mem.reltol, 1.0e-4);
        assert_eq!(ark_mem.Sabstol, 1.0e-9);
        assert_eq!(ark_mem.itol, ARK_SS);
        assert_eq!(ark_mem.mxstep, MXSTEP_DEFAULT);
        assert_eq!(ark_mem.mxhnil, MXHNIL);
        assert_eq!(ark_mem.maxnef, MAXNEF);
        assert_eq!(ark_mem.maxncf, MAXNCF);
        let hm = ark_mem.hadapt_mem.as_ref().unwrap();
        assert_eq!(hm.etamx1, 10000.0);
        assert_eq!(hm.etamxf, 0.3);
        assert_eq!(hm.safety, 0.9);
        assert_eq!(hm.growth, 20.0);
        assert_eq!(hm.small_nef, 2);
    }

    #[test]
    fn internal_weight_functions() {
        let ctx = SUNContext_Create();
        let mut ark_mem = arkCreate(&ctx);
        ark_mem.reltol = 1.0e-2;
        ark_mem.Sabstol = 1.0e-3;

        let mut y = NVector::new(2);
        y.data[0] = -2.0;
        y.data[1] = 0.5;
        let mut w = NVector::new(2);

        /* SS: w_i = 1/(reltol*|y_i| + abstol) */
        assert_eq!(arkEwtSetSS(&ark_mem, &y, &mut w), 0);
        assert_eq!(w.data[0], 1.0 / (1.0e-2 * 2.0 + 1.0e-3));
        assert_eq!(w.data[1], 1.0 / (1.0e-2 * 0.5 + 1.0e-3));

        /* SV with a zero abstol component and atolmin0: negative test */
        let mut vab = NVector::new(2);
        vab.data[0] = 0.0;
        vab.data[1] = 1.0e-3;
        ark_mem.Vabstol = Some(vab);
        assert_eq!(arkEwtSetSV(&ark_mem, &y, &mut w), 0);
        assert_eq!(w.data[0], 1.0 / (1.0e-2 * 2.0));
        ark_mem.atolmin0 = true;
        let mut y0 = NVector::new(2);
        y0.data[1] = 1.0; /* y0[0] = 0 with abstol[0] = 0 -> failure */
        assert_eq!(arkEwtSetSV(&ark_mem, &y0, &mut w), -1);

        /* SmallReal fills with SUN_SMALL_REAL */
        assert_eq!(arkEwtSetSmallReal(&y, &mut w), ARK_SUCCESS);
        assert_eq!(w.data[0], SUN_SMALL_REAL);
    }
}

/*---------------------------------------------------------------
  ARKodeSetUserData:

  Specifies the user data pointer for f
  ---------------------------------------------------------------*/
pub fn ARKodeSetUserData(ark_mem: &mut ARKodeMem, user_data: crate::sundials_types::UserData) -> i32 {
    ark_mem.user_data = user_data;

    /* efun/rfun data and root_data follow user_data automatically in
    this port (the dispatch helpers read ark_mem.user_data) */

    /* Set user data into stepper (if provided; the stepper op reads
    ark_mem.user_data itself) */
    if let Some(step_setuserdata) = ark_mem.step_setuserdata {
        return step_setuserdata(ark_mem);
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetInitStep:

  Specifies the initial step size to be attempted.  Passing 0
  sets the default, otherwise use input.
  ---------------------------------------------------------------*/
pub fn ARKodeSetInitStep(ark_mem: &mut ARKodeMem, hin: f64) -> i32 {
    use crate::arkode_impl::{arkProcessError, ARK_CONTROLLER_ERR, ARK_STEPPER_UNSUPPORTED};

    /* Guard against hin==0 for non-adaptive time stepper modules */
    if !ark_mem.step_supports_adaptive && hin == ZERO {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!(),
            "ARKodeSetInitStep",
            file!(),
            "time-stepping module does not support temporal adaptivity",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* Passing hin=0 sets the default, otherwise use input. */
    if hin == ZERO {
        ark_mem.hin = ZERO;
    } else {
        ark_mem.hin = hin;
    }

    /* Clear previous initial step */
    ark_mem.h0u = ZERO;

    /* Reset error controller (e.g., error and step size history) */
    if let Some(uc) = ark_mem.hadapt_mem.as_mut().unwrap().usercontrol.as_mut() {
        let _ = crate::arkode_user_controller::SUNAdaptController_Reset_ARKUserControl(uc);
    }
    if let Some(hcontroller) = ark_mem.hadapt_mem.as_mut().unwrap().hcontroller.as_mut() {
        let retval = crate::sundials_adaptcontroller::SUNAdaptController_Reset(hcontroller);
        if retval != crate::sundials_errors::SUN_SUCCESS {
            return ARK_CONTROLLER_ERR;
        }
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetFixedStep:

  Specifies to use a fixed time step size instead of performing
  any form of temporal adaptivity.  ARKODE will use this step size
  for all steps (unless tstop is set, in which case it may need to
  modify that last step approaching tstop.  If any solver failure
  occurs in the timestepping module, ARKODE will typically
  immediately return with an error message indicating that the
  selected step size cannot be used.

  Any nonzero argument will result in the use of that fixed step
  size; an argument of 0 will re-enable temporal adaptivity.
  ---------------------------------------------------------------*/
pub fn ARKodeSetFixedStep(ark_mem: &mut ARKodeMem, hfixed: f64) -> i32 {
    use crate::arkode_impl::{arkProcessError, ARK_STEPPER_UNSUPPORTED, ARK_SV};

    /* ensure that when hfixed=0, the time step module supports adaptivity */
    if hfixed == ZERO && !ark_mem.step_supports_adaptive {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!(),
            "ARKodeSetFixedStep",
            file!(),
            "temporal adaptivity is not supported by this time step module",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* re-attach internal error weight functions if necessary */
    if hfixed == ZERO && !ark_mem.user_efun {
        let retval = if ark_mem.itol == ARK_SV && ark_mem.Vabstol.is_some() {
            let vabstol = ark_mem.Vabstol.take().unwrap();
            let r = crate::arkode::ARKodeSVtolerances(ark_mem, ark_mem.reltol, &vabstol);
            ark_mem.Vabstol = Some(vabstol);
            r
        } else {
            crate::arkode::ARKodeSStolerances(ark_mem, ark_mem.reltol, ark_mem.Sabstol)
        };
        if retval != ARK_SUCCESS {
            return retval;
        }
    }

    /* set ark_mem "fixedstep" entry */
    if hfixed != ZERO {
        ark_mem.fixedstep = true;
        ark_mem.hin = hfixed;
    } else {
        ark_mem.fixedstep = false;
    }

    /* Notify ARKODE to use hfixed as the initial step size, and return */
    ARKodeSetInitStep(ark_mem, hfixed)
}

/*---------------------------------------------------------------
  ARKodeSetStepDirection:

  Specifies the direction of integration (forward or backward)
  based on the sign of stepdir. If 0, the direction will remain
  unchanged. Note that if a fixed step size was previously set,
  this function can change the sign of that.

  This should only be called after ARKodeReset, or between
  creating a stepper and ARKodeEvolve.
  ---------------------------------------------------------------*/
pub fn ARKodeSetStepDirection(ark_mem: &mut ARKodeMem, stepdir: f64) -> i32 {
    /* stepdir is a sunrealtype because the direction typically comes from a time
     * step h or tend-tstart which are sunrealtypes. If stepdir was in int,
     * conversions would be required which can cause undefined behavior when
     * greater than MAX_INT */
    use crate::arkode_impl::{arkProcessError, ARK_CONTROLLER_ERR, ARK_STEP_DIRECTION_ERR};
    use crate::sundials_math::SUNRcopysign;

    /* do not change direction once the module has been initialized i.e., after calling
       ARKodeEvolve unless ReInit or Reset are called. */
    if !ark_mem.initsetup {
        arkProcessError(
            Some(ark_mem),
            ARK_STEP_DIRECTION_ERR,
            line!(),
            "ARKodeSetStepDirection",
            file!(),
            "Step direction cannot be specified after module initialization.",
        );
        return ARK_STEP_DIRECTION_ERR;
    }

    if stepdir != ZERO {
        let mut h = ZERO;
        let retval = ARKodeGetStepDirection(ark_mem, &mut h);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                retval,
                line!(),
                "ARKodeSetStepDirection",
                file!(),
                "Unable to access step direction",
            );
            return retval;
        }

        if h != SUNRcopysign(h, stepdir) {
            /* Reverse the sign of h. If adaptive, h will be overwritten anyway by the
             * initial step estimation since ARKodeReset must be called before this.
             * However, the sign of h will be used to check if the integration
             * direction and stop time are consistent, e.g., in ARKodeSetStopTime, so
             * we should not set h = 0. */
            ark_mem.h = -h;
            /* Clear previous initial step and force an initial step recomputation.
             * Normally, this would not occur after a reset, but it is necessary here
             * because the timestep used in one direction may not be suitable for the
             * other */
            ark_mem.h0u = ZERO;
            /* Reverse the step if in fixed mode. If adaptive, reset to 0 to clear any
             * old value from a call to ARKodeSetInit */
            ark_mem.hin = if ark_mem.fixedstep { -h } else { ZERO };

            /* Reset error controller (e.g., error and step size history) */
            if let Some(hadapt_mem) = ark_mem.hadapt_mem.as_mut() {
                if let Some(uc) = hadapt_mem.usercontrol.as_mut() {
                    let _ =
                        crate::arkode_user_controller::SUNAdaptController_Reset_ARKUserControl(uc);
                }
                if let Some(hcontroller) = hadapt_mem.hcontroller.as_mut() {
                    let err =
                        crate::sundials_adaptcontroller::SUNAdaptController_Reset(hcontroller);
                    if err != crate::sundials_errors::SUN_SUCCESS {
                        arkProcessError(
                            Some(ark_mem),
                            ARK_CONTROLLER_ERR,
                            line!(),
                            "ARKodeSetStepDirection",
                            file!(),
                            "Unable to reset error controller object",
                        );
                        return ARK_CONTROLLER_ERR;
                    }
                }
            }
        }
    }

    if let Some(step_setstepdirection) = ark_mem.step_setstepdirection {
        return step_setstepdirection(ark_mem, stepdir);
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeGetStepDirection:

  Gets the direction of integration (forward or backward) based
  on the sign of stepdir. A value of 0 indicates integration can
  proceed in either direction.
  ---------------------------------------------------------------*/
pub fn ARKodeGetStepDirection(ark_mem: &mut ARKodeMem, stepdir: &mut f64) -> i32 {
    *stepdir = if ark_mem.fixedstep || ark_mem.h == ZERO {
        ark_mem.hin
    } else {
        ark_mem.h
    };
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetMaxNumSteps:

  Specifies the maximum number of integration steps
  ---------------------------------------------------------------*/
pub fn ARKodeSetMaxNumSteps(ark_mem: &mut ARKodeMem, mxsteps: i64) -> i32 {
    /* Passing mxsteps=0 sets the default. Passing mxsteps<0 disables the
    test. */
    if mxsteps == 0 {
        ark_mem.mxstep = MXSTEP_DEFAULT;
    } else {
        ark_mem.mxstep = mxsteps;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetStopTime:

  Specifies the time beyond which the integration is not to proceed.
  ---------------------------------------------------------------*/
pub fn ARKodeSetStopTime(ark_mem: &mut ARKodeMem, tstop: f64) -> i32 {
    use crate::arkode_impl::{arkProcessError, ARK_ILL_INPUT};
    use crate::sundials_utils::fmt_g;

    /* If ARKODE was called at least once, test if tstop is legal (i.e. if
    it was not already passed). If ARKodeSetStopTime is called before the
    first call to ARKODE, tstop will be checked in ARKODE. */
    if ark_mem.nst > 0 {
        if (tstop - ark_mem.tcur) * ark_mem.h < ZERO {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "ARKodeSetStopTime",
                file!(),
                &format!(
                    "The value tstop = {} is behind current t = {} in the direction of integration.",
                    fmt_g(tstop, 0, 15),
                    fmt_g(ark_mem.tcur, 0, 15)
                ),
            );
            return ARK_ILL_INPUT;
        }
    }

    ark_mem.tstop = tstop;
    ark_mem.tstopset = true;

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetInterpolantDegree:

  Specifies the polynomial degree for the dense output
  interpolation module.
  ---------------------------------------------------------------*/
pub fn ARKodeSetInterpolantDegree(ark_mem: &mut ARKodeMem, degree: i32) -> i32 {
    use crate::arkode_impl::{arkProcessError, ARK_ILL_INPUT, ARK_INTERP_FAIL, ARK_INTERP_MAX_DEGREE};

    /* do not change degree once the module has been initialized */
    if ark_mem.initialized {
        arkProcessError(
            Some(ark_mem),
            ARK_INTERP_FAIL,
            line!(),
            "ARKodeSetInterpolantDegree",
            file!(),
            "Degree cannot be specified after module initialization.",
        );
        return ARK_ILL_INPUT;
    }

    if degree > ARK_INTERP_MAX_DEGREE {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "ARKodeSetInterpolantDegree",
            file!(),
            "Illegal degree specified.",
        );
        return ARK_ILL_INPUT;
    } else if degree < 0 {
        ark_mem.interp_degree = ARK_INTERP_MAX_DEGREE;
    } else {
        ark_mem.interp_degree = degree;
    }

    /* Set the degree now if possible otherwise it will be used when
    creating the interpolation module */
    if ark_mem.interp.is_some() {
        return crate::arkode_interp::arkInterpSetDegree(ark_mem, ark_mem.interp_degree);
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkReplaceAdaptController:

  Replaces the current SUNAdaptController time step controller
  object. On NULL-valued input the default (I) controller is
  created.
  ---------------------------------------------------------------*/
pub fn arkReplaceAdaptController(
    ark_mem: &mut ARKodeMem,
    c: Option<crate::sundials_adaptcontroller::SUNAdaptController>,
    take_ownership: bool,
) -> i32 {
    use crate::sundials_adaptcontroller::SUNAdaptController_Space;

    let mut lenrw: i64 = 0;
    let mut leniw: i64 = 0;

    /* Remove current SUNAdaptController object
    (delete if owned, and then nullify pointer); an owned ARKUserControl
    wrapper lives in the usercontrol slot */
    if let Some(uc) = ark_mem.hadapt_mem.as_mut().unwrap().usercontrol.take() {
        let _ = crate::arkode_user_controller::SUNAdaptController_Space_ARKUserControl(
            &uc, &mut lenrw, &mut leniw,
        );
        ark_mem.liw -= leniw;
        ark_mem.lrw -= lenrw;
    }
    let hadapt_mem = ark_mem.hadapt_mem.as_mut().unwrap();
    if hadapt_mem.owncontroller {
        if let Some(hc) = hadapt_mem.hcontroller.as_ref() {
            let retval = SUNAdaptController_Space(hc, &mut lenrw, &mut leniw);
            if retval == crate::sundials_errors::SUN_SUCCESS {
                ark_mem.liw -= leniw;
                ark_mem.lrw -= lenrw;
            }
        }

        /* SUNAdaptController_Destroy = drop */
        let hadapt_mem = ark_mem.hadapt_mem.as_mut().unwrap();
        hadapt_mem.hcontroller = None;
        hadapt_mem.owncontroller = false;
    }
    ark_mem.hadapt_mem.as_mut().unwrap().hcontroller = None;

    /* On NULL-valued input, create default SUNAdaptController object */
    let c = match c {
        None => {
            ark_mem.hadapt_mem.as_mut().unwrap().owncontroller = true;
            crate::sunadaptcontroller_soderlind::SUNAdaptController_I()
        }
        Some(c) => {
            ark_mem.hadapt_mem.as_mut().unwrap().owncontroller = take_ownership;
            c
        }
    };

    /* Attach new SUNAdaptController object */
    let retval = SUNAdaptController_Space(&c, &mut lenrw, &mut leniw);
    if retval == crate::sundials_errors::SUN_SUCCESS {
        ark_mem.liw += leniw;
        ark_mem.lrw += lenrw;
    }
    ark_mem.hadapt_mem.as_mut().unwrap().hcontroller = Some(c);

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkSetAdaptivityMethod:

  Specifies the built-in time step adaptivity algorithm (and
  optionally, its associated parameters) to use.  All parameters
  will be checked for validity when used by the solver.

  Users should transition to constructing non-default SUNAdaptController
  objects directly, and providing those directly to the integrator
  via the time-stepping module *SetController routines.
  ---------------------------------------------------------------*/
pub fn arkSetAdaptivityMethod(
    ark_mem: &mut ARKodeMem,
    imethod: i32,
    idefault: i32,
    pq: i32,
    adapt_params: Option<&[f64; 3]>,
) -> i32 {
    use crate::sunadaptcontroller_imexgus::{
        SUNAdaptController_ImExGus, SUNAdaptController_SetParams_ImExGus,
    };
    use crate::sunadaptcontroller_soderlind::{
        SUNAdaptController_ExpGus, SUNAdaptController_I, SUNAdaptController_ImpGus,
        SUNAdaptController_PI, SUNAdaptController_PID, SUNAdaptController_SetParams_ExpGus,
        SUNAdaptController_SetParams_I, SUNAdaptController_SetParams_ImpGus,
        SUNAdaptController_SetParams_PI, SUNAdaptController_SetParams_PID,
    };
    use crate::arkode_impl::{ARK_CONTROLLER_ERR, ARK_ILL_INPUT};
    use crate::sundials_adaptcontroller::SUNAdaptController_Space;
    use crate::sundials_errors::SUN_SUCCESS;

    /* the ARK_ADAPT_* method constants (arkode.h) */
    const ARK_ADAPT_PID: i32 = 0;
    const ARK_ADAPT_PI: i32 = 1;
    const ARK_ADAPT_I: i32 = 2;
    const ARK_ADAPT_EXP_GUS: i32 = 3;
    const ARK_ADAPT_IMP_GUS: i32 = 4;
    const ARK_ADAPT_IMEX_GUS: i32 = 5;

    /* Check for illegal inputs */
    if idefault != 1 && adapt_params.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "arkSetAdaptivityMethod",
            file!(),
            "NULL-valued adapt_params provided",
        );
        return ARK_ILL_INPUT;
    }

    /* Remove current SUNAdaptController object
       (delete if owned, and then nullify pointer); an owned
       ARKUserControl wrapper lives in the usercontrol slot */
    {
        let mut lenrw: i64 = 0;
        let mut leniw: i64 = 0;
        if let Some(uc) = ark_mem.hadapt_mem.as_mut().unwrap().usercontrol.take() {
            let _ = crate::arkode_user_controller::SUNAdaptController_Space_ARKUserControl(
                &uc, &mut lenrw, &mut leniw,
            );
            ark_mem.liw -= leniw;
            ark_mem.lrw -= lenrw;
        }
        let hadapt_mem = ark_mem.hadapt_mem.as_mut().unwrap();
        if hadapt_mem.owncontroller {
            if let Some(hc) = hadapt_mem.hcontroller.as_ref() {
                let retval = SUNAdaptController_Space(hc, &mut lenrw, &mut leniw);
                if retval == SUN_SUCCESS {
                    ark_mem.liw -= leniw;
                    ark_mem.lrw -= lenrw;
                }
            }
            /* SUNAdaptController_Destroy = drop */
            let hadapt_mem = ark_mem.hadapt_mem.as_mut().unwrap();
            hadapt_mem.owncontroller = false;
        }
        ark_mem.hadapt_mem.as_mut().unwrap().hcontroller = None;
    }

    /* set adaptivity parameters from inputs */
    let mut k1 = ZERO;
    let mut k2 = ZERO;
    let mut k3 = ZERO;
    if idefault != 1 {
        let p = adapt_params.unwrap();
        k1 = p[0];
        k2 = p[1];
        k3 = p[2];
    }
    ark_mem.hadapt_mem.as_mut().unwrap().pq = pq;

    /* Create new SUNAdaptController object based on "imethod" input,
       optionally setting the specified controller parameters */
    macro_rules! params_or_controller_err {
        ($retval:expr, $msg:expr) => {
            if $retval != SUN_SUCCESS {
                arkProcessError(
                    Some(ark_mem),
                    ARK_CONTROLLER_ERR,
                    line!(),
                    "arkSetAdaptivityMethod",
                    file!(),
                    $msg,
                );
                return ARK_CONTROLLER_ERR;
            }
        };
    }
    let c = match imethod {
        ARK_ADAPT_PID => {
            let mut c = SUNAdaptController_PID();
            if idefault != 1 {
                let retval = SUNAdaptController_SetParams_PID(&mut c, k1, -k2, k3);
                params_or_controller_err!(retval, "SUNAdaptController_SetParams_PID failure");
            }
            c
        }
        ARK_ADAPT_PI => {
            let mut c = SUNAdaptController_PI();
            if idefault != 1 {
                let retval = SUNAdaptController_SetParams_PI(&mut c, k1, -k2);
                params_or_controller_err!(retval, "SUNAdaptController_SetParams_PI failure");
            }
            c
        }
        ARK_ADAPT_I => {
            let mut c = SUNAdaptController_I();
            if idefault != 1 {
                let retval = SUNAdaptController_SetParams_I(&mut c, k1);
                params_or_controller_err!(retval, "SUNAdaptController_SetParams_I failure");
            }
            c
        }
        ARK_ADAPT_EXP_GUS => {
            let mut c = SUNAdaptController_ExpGus();
            if idefault != 1 {
                let retval = SUNAdaptController_SetParams_ExpGus(&mut c, k1, k2);
                params_or_controller_err!(retval, "SUNAdaptController_SetParams_ExpGus failure");
            }
            c
        }
        ARK_ADAPT_IMP_GUS => {
            let mut c = SUNAdaptController_ImpGus();
            if idefault != 1 {
                let retval = SUNAdaptController_SetParams_ImpGus(&mut c, k1, k2);
                params_or_controller_err!(retval, "SUNAdaptController_SetParams_ImpGus failure");
            }
            c
        }
        ARK_ADAPT_IMEX_GUS => {
            let mut c = SUNAdaptController_ImExGus();
            if idefault != 1 {
                let retval = SUNAdaptController_SetParams_ImExGus(&mut c, k1, k2, k3, k3);
                params_or_controller_err!(retval, "SUNAdaptController_SetParams_ImExGus failure");
            }
            c
        }
        _ => {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "arkSetAdaptivityMethod",
                file!(),
                "Illegal imethod",
            );
            return ARK_ILL_INPUT;
        }
    };

    /* Attach new SUNAdaptController object */
    {
        let mut lenrw: i64 = 0;
        let mut leniw: i64 = 0;
        let retval = SUNAdaptController_Space(&c, &mut lenrw, &mut leniw);
        if retval == SUN_SUCCESS {
            ark_mem.liw += leniw;
            ark_mem.lrw += lenrw;
        }
        let hadapt_mem = ark_mem.hadapt_mem.as_mut().unwrap();
        hadapt_mem.hcontroller = Some(c);
        hadapt_mem.owncontroller = true;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkSetAdaptivityFn:

  Specifies the user-provided time step adaptivity function to use.
  If 'hfun' is NULL-valued, then the default I controller will
  be used instead.

  Users should transition to constructing a custom SUNAdaptController
  object, and providing this directly to the integrator
  via the time-stepping module *SetController routines.
  ---------------------------------------------------------------*/
pub fn arkSetAdaptivityFn(
    ark_mem: &mut ARKodeMem,
    hfun: Option<crate::arkode_impl::ARKAdaptFn>,
    h_data: crate::sundials_types::UserData,
) -> i32 {
    use crate::sundials_adaptcontroller::SUNAdaptController_Space;
    use crate::sundials_errors::SUN_SUCCESS;

    /* Remove current SUNAdaptController object
       (delete if owned, and then nullify pointer); an owned
       ARKUserControl wrapper lives in the usercontrol slot */
    {
        let mut lenrw: i64 = 0;
        let mut leniw: i64 = 0;
        if let Some(uc) = ark_mem.hadapt_mem.as_mut().unwrap().usercontrol.take() {
            let _ = crate::arkode_user_controller::SUNAdaptController_Space_ARKUserControl(
                &uc, &mut lenrw, &mut leniw,
            );
            ark_mem.liw -= leniw;
            ark_mem.lrw -= lenrw;
        }
        let hadapt_mem = ark_mem.hadapt_mem.as_mut().unwrap();
        if hadapt_mem.owncontroller {
            if let Some(hc) = hadapt_mem.hcontroller.as_ref() {
                let retval = SUNAdaptController_Space(hc, &mut lenrw, &mut leniw);
                if retval == SUN_SUCCESS {
                    ark_mem.liw -= leniw;
                    ark_mem.lrw -= lenrw;
                }
            }
            /* SUNAdaptController_Destroy = drop */
            let hadapt_mem = ark_mem.hadapt_mem.as_mut().unwrap();
            hadapt_mem.owncontroller = false;
        }
        ark_mem.hadapt_mem.as_mut().unwrap().hcontroller = None;
    }

    /* Create new SUNAdaptController object depending on NULL-ity of 'hfun' */
    match hfun {
        None => {
            let c = crate::sunadaptcontroller_soderlind::SUNAdaptController_I();
            /* Attach new SUNAdaptController object */
            let mut lenrw: i64 = 0;
            let mut leniw: i64 = 0;
            let retval = SUNAdaptController_Space(&c, &mut lenrw, &mut leniw);
            if retval == SUN_SUCCESS {
                ark_mem.liw += leniw;
                ark_mem.lrw += lenrw;
            }
            let hadapt_mem = ark_mem.hadapt_mem.as_mut().unwrap();
            hadapt_mem.hcontroller = Some(c);
            hadapt_mem.owncontroller = true;
        }
        Some(hfun) => {
            let uc = crate::arkode_user_controller::ARKUserControl(hfun, h_data);
            /* Attach new SUNAdaptController object (the ARKUserControl
               wrapper occupies the usercontrol slot; ownership implied) */
            let mut lenrw: i64 = 0;
            let mut leniw: i64 = 0;
            let _ = crate::arkode_user_controller::SUNAdaptController_Space_ARKUserControl(
                &uc, &mut lenrw, &mut leniw,
            );
            ark_mem.liw += leniw;
            ark_mem.lrw += lenrw;
            ark_mem.hadapt_mem.as_mut().unwrap().usercontrol = Some(uc);
        }
    }

    ARK_SUCCESS
}

/* -----------------------------------------------------------------
 * Counterparts of sunfprintf_real / sunfprintf_long
 * (src/sundials/sundials_utils.h). SUN_FORMAT_G is "%.15g" and
 * SUN_FORMAT_E is "% .15e" for double precision.
 * -----------------------------------------------------------------*/

const SUN_TABLE_WIDTH: usize = 29;

pub(crate) fn sunfprintf_real(
    outfile: &mut dyn std::io::Write,
    fmt: crate::sundials_types::SUNOutputFormat,
    start: bool,
    name: &str,
    value: f64,
) {
    use crate::sundials_types::SUN_OUTPUTFORMAT_TABLE;
    use crate::sundials_utils::{fmt_e, fmt_g};
    if fmt == SUN_OUTPUTFORMAT_TABLE {
        let _ = writeln!(outfile, "{:<width$} = {}", name, fmt_g(value, 0, 15),
                         width = SUN_TABLE_WIDTH);
    } else {
        if !start {
            let _ = write!(outfile, ",");
        }
        /* C "% .15e": a space is printed in place of a plus sign */
        let e = fmt_e(value, 0, 15);
        let e = if e.starts_with('-') { e } else { format!(" {}", e) };
        let _ = write!(outfile, "{},{}", name, e);
    }
}

pub(crate) fn sunfprintf_long(
    outfile: &mut dyn std::io::Write,
    fmt: crate::sundials_types::SUNOutputFormat,
    start: bool,
    name: &str,
    value: i64,
) {
    use crate::sundials_types::SUN_OUTPUTFORMAT_TABLE;
    if fmt == SUN_OUTPUTFORMAT_TABLE {
        let _ = writeln!(outfile, "{:<width$} = {}", name, value, width = SUN_TABLE_WIDTH);
    } else {
        if !start {
            let _ = write!(outfile, ",");
        }
        let _ = write!(outfile, "{},{}", name, value);
    }
}

/*---------------------------------------------------------------
  ARKodePrintAllStats:

  Prints the current value of all statistics
  ---------------------------------------------------------------*/
pub fn ARKodePrintAllStats(
    ark_mem: &mut ARKodeMem,
    outfile: &mut dyn std::io::Write,
    fmt: crate::sundials_types::SUNOutputFormat,
) -> i32 {
    use crate::sundials_types::{SUNFALSE, SUNTRUE};

    /* (invalid formatting options are unrepresentable in the enum) */

    sunfprintf_real(outfile, fmt, SUNTRUE, "Current time", ark_mem.tcur);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Steps", ark_mem.nst);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Step attempts", ark_mem.nst_attempts);
    sunfprintf_long(
        outfile,
        fmt,
        SUNFALSE,
        "Stability limited steps",
        ark_mem.hadapt_mem.as_ref().unwrap().nst_exp,
    );
    sunfprintf_long(
        outfile,
        fmt,
        SUNFALSE,
        "Accuracy limited steps",
        ark_mem.hadapt_mem.as_ref().unwrap().nst_acc,
    );
    sunfprintf_long(outfile, fmt, SUNFALSE, "Error test fails", ark_mem.netf);
    sunfprintf_long(outfile, fmt, SUNFALSE, "NLS step fails", ark_mem.ncfn);
    sunfprintf_long(
        outfile,
        fmt,
        SUNFALSE,
        "Inequality constraint fails",
        ark_mem.nconstrfails,
    );
    sunfprintf_real(outfile, fmt, SUNFALSE, "Initial step size", ark_mem.h0u);
    sunfprintf_real(outfile, fmt, SUNFALSE, "Last step size", ark_mem.hold);
    sunfprintf_real(outfile, fmt, SUNFALSE, "Current step size", ark_mem.next_h);
    if let Some(root_mem) = ark_mem.root_mem.as_ref() {
        sunfprintf_long(outfile, fmt, SUNFALSE, "Root fn evals", root_mem.nge);
    }

    /* Print relaxation stats */
    if ark_mem.relax_enabled {
        let retval = crate::arkode_relaxation::arkRelaxPrintAllStats(ark_mem, outfile, fmt);
        if retval != ARK_SUCCESS {
            return retval;
        }
    }

    /* Print stepper stats (if provided) */
    if let Some(step_printallstats) = ark_mem.step_printallstats {
        return step_printallstats(ark_mem, outfile, fmt);
    }

    ARK_SUCCESS
}

/*===============================================================
  Stepper-op dispatch wrappers (arkode_io.c) — implicit-solver
  option/stat families, added with the ARKStep port.
  ===============================================================*/

/* Set routines guarded by step_supports_implicit; fallback error
   matches C ("time-stepping module does not support this function"). */
macro_rules! ark_io_dispatch_implicit {
    ($name:ident, $op:ident $(, $arg:ident : $ty:ty)*) => {
        pub fn $name(ark_mem: &mut ARKodeMem $(, $arg: $ty)*) -> i32 {
            /* Guard against use for time steppers that do not need an
               algebraic solver */
            if !ark_mem.step_supports_implicit {
                arkProcessError(
                    Some(ark_mem),
                    crate::arkode_impl::ARK_STEPPER_UNSUPPORTED,
                    line!(),
                    stringify!($name),
                    file!(),
                    "time-stepping module does not require an algebraic solver",
                );
                return crate::arkode_impl::ARK_STEPPER_UNSUPPORTED;
            }

            /* Call stepper routine (if provided) */
            if let Some(op) = ark_mem.$op {
                return op(ark_mem $(, $arg)*);
            }
            arkProcessError(
                Some(ark_mem),
                crate::arkode_impl::ARK_STEPPER_UNSUPPORTED,
                line!(),
                stringify!($name),
                file!(),
                "time-stepping module does not support this function",
            );
            crate::arkode_impl::ARK_STEPPER_UNSUPPORTED
        }
    };
}

/* Stat getters with a zero fallback when the stepper lacks the op. */
macro_rules! ark_io_stat_zero_fallback {
    ($name:ident, $op:ident $(, $arg:ident : $ty:ty)*) => {
        pub fn $name(ark_mem: &mut ARKodeMem $(, $arg: $ty)*) -> i32 {
            /* Call stepper routine (if provided) */
            if let Some(op) = ark_mem.$op {
                return op(ark_mem $(, $arg)*);
            }
            $( *$arg = 0; )*
            ARK_SUCCESS
        }
    };
}

ark_io_dispatch_implicit!(
    ARKodeSetNonlinearSolver,
    step_setnonlinearsolver,
    nls: crate::sundials_nonlinearsolver::NonlinearSolver
);
ark_io_dispatch_implicit!(ARKodeSetLinear, step_setlinear, timedepend: i32);
ark_io_dispatch_implicit!(ARKodeSetNonlinear, step_setnonlinear);
ark_io_dispatch_implicit!(ARKodeSetAutonomous, step_setautonomous, autonomous: bool);
ark_io_dispatch_implicit!(
    ARKodeSetNlsRhsFn,
    step_setnlsrhsfn,
    nls_fi: Option<crate::arkode_impl::ARKRhsFn>
);
ark_io_dispatch_implicit!(
    ARKodeSetDeduceImplicitRhs,
    step_setdeduceimplicitrhs,
    deduce: bool
);
/* ARKodeSetPredictorMethod: full version below (adds C's ARK_INTERP_NONE check) */
ark_io_dispatch_implicit!(ARKodeSetMaxNonlinIters, step_setmaxnonliniters, maxcor: i32);
ark_io_dispatch_implicit!(ARKodeSetNonlinConvCoef, step_setnonlinconvcoef, nlscoef: f64);
ark_io_dispatch_implicit!(ARKodeSetNonlinCRDown, step_setnonlincrdown, crdown: f64);
ark_io_dispatch_implicit!(ARKodeSetNonlinRDiv, step_setnonlinrdiv, rdiv: f64);
ark_io_dispatch_implicit!(ARKodeSetDeltaGammaMax, step_setdeltagammamax, dgmax: f64);
ark_io_dispatch_implicit!(ARKodeSetLSetupFrequency, step_setlsetupfrequency, msbp: i32);
ark_io_dispatch_implicit!(
    ARKodeSetStagePredictFn,
    step_setstagepredictfn,
    predict_stage: Option<crate::arkode_impl::ARKStagePredictFn>
);
ark_io_dispatch_implicit!(ARKodeGetCurrentGamma, step_getcurrentgamma, gamma: &mut f64);

ark_io_stat_zero_fallback!(
    ARKodeGetNumLinSolvSetups,
    step_getnumlinsolvsetups,
    nlinsetups: &mut i64
);
ark_io_stat_zero_fallback!(
    ARKodeGetNumNonlinSolvIters,
    step_getnumnonlinsolviters,
    nniters: &mut i64
);
ark_io_stat_zero_fallback!(
    ARKodeGetNumNonlinSolvConvFails,
    step_getnumnonlinsolvconvfails,
    nnfails: &mut i64
);
ark_io_stat_zero_fallback!(
    ARKodeGetNonlinSolvStats,
    step_getnonlinsolvstats,
    nniters: &mut i64,
    nnfails: &mut i64
);

/*---------------------------------------------------------------
  ARKodeGetNumRhsEvals: dispatches to the stepper (no zero
  fallback in C).
  ---------------------------------------------------------------*/
pub fn ARKodeGetNumRhsEvals(
    ark_mem: &mut ARKodeMem,
    partition_index: i32,
    num_rhs_evals: &mut i64,
) -> i32 {
    /* Call stepper routine (if provided) */
    if let Some(op) = ark_mem.step_getnumrhsevals {
        return op(ark_mem, partition_index, num_rhs_evals);
    }
    arkProcessError(
        Some(ark_mem),
        crate::arkode_impl::ARK_STEPPER_UNSUPPORTED,
        line!(),
        "ARKodeGetNumRhsEvals",
        file!(),
        "time-stepping module does not support this function",
    );
    crate::arkode_impl::ARK_STEPPER_UNSUPPORTED
}

/*---------------------------------------------------------------
  ARKodeGetEstLocalErrors: Returns the current local truncation
  error estimate vector.
  ---------------------------------------------------------------*/
pub fn ARKodeGetEstLocalErrors(
    ark_mem: &mut ARKodeMem,
    ele: &mut crate::nvector_serial::NVector,
) -> i32 {
    /* Call stepper-specific routine (if provided); otherwise return an error */
    if let Some(op) = ark_mem.step_getestlocalerrors {
        return op(ark_mem, ele);
    }
    arkProcessError(
        Some(ark_mem),
        crate::arkode_impl::ARK_STEPPER_UNSUPPORTED,
        line!(),
        "ARKodeGetEstLocalErrors",
        file!(),
        "time-stepping module does provide a temporal error estimate",
    );
    crate::arkode_impl::ARK_STEPPER_UNSUPPORTED
}

/*---------------------------------------------------------------
  ARKodeGetNumSteps:

  Returns the current number of integration steps
  ---------------------------------------------------------------*/
pub fn ARKodeGetNumSteps(ark_mem: &mut ARKodeMem, nsteps: &mut i64) -> i32 {
    *nsteps = ark_mem.nst;
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeGetNumStepAttempts:

  Returns the current number of steps attempted by the solver
  ---------------------------------------------------------------*/
pub fn ARKodeGetNumStepAttempts(ark_mem: &mut ARKodeMem, nstep_attempts: &mut i64) -> i32 {
    *nstep_attempts = ark_mem.nst_attempts;
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeGetNumErrTestFails:

  Returns the current number of error test failures
  ---------------------------------------------------------------*/
pub fn ARKodeGetNumErrTestFails(ark_mem: &mut ARKodeMem, netfails: &mut i64) -> i32 {
    *netfails = ark_mem.netf;
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeWriteParameters:

  Outputs all solver parameters to the provided file pointer.
  ---------------------------------------------------------------*/
pub fn ARKodeWriteParameters(ark_mem: &mut ARKodeMem, fp: &mut dyn std::io::Write) -> i32 {
    use crate::arkode_impl::{ARK_WF, ONE};
    use crate::sundials_utils::fmt_g;

    /* print integrator parameters to file */
    let _ = write!(fp, "ARKODE solver parameters:\n");
    if ark_mem.hmin != ZERO {
        let _ = write!(fp, "  Minimum step size = {}\n", fmt_g(ark_mem.hmin, 0, 15));
    }
    if ark_mem.hmax_inv != ZERO {
        let _ = write!(fp, "  Maximum step size = {}\n", fmt_g(ONE / ark_mem.hmax_inv, 0, 15));
    }
    if ark_mem.fixedstep {
        let _ = write!(fp, "  Fixed time-stepping enabled\n");
    }
    if ark_mem.itol == ARK_WF {
        let _ = write!(fp, "  User provided error weight function\n");
    } else {
        let _ = write!(fp, "  Solver relative tolerance = {}\n", fmt_g(ark_mem.reltol, 0, 15));
        if ark_mem.itol == ARK_SS {
            let _ = write!(fp, "  Solver absolute tolerance = {}\n", fmt_g(ark_mem.Sabstol, 0, 15));
        } else {
            let _ = write!(fp, "  Vector-valued solver absolute tolerance\n");
        }
    }
    if !ark_mem.rwt_is_ewt {
        if ark_mem.ritol == ARK_WF {
            let _ = write!(fp, "  User provided residual weight function\n");
        } else {
            if ark_mem.ritol == ARK_SS {
                let _ = write!(
                    fp,
                    "  Absolute residual tolerance = {}\n",
                    fmt_g(ark_mem.SRabstol, 0, 15)
                );
            } else {
                let _ = write!(fp, "  Vector-valued residual absolute tolerance\n");
            }
        }
    }
    if ark_mem.hin != ZERO {
        let _ = write!(fp, "  Initial step size = {}\n", fmt_g(ark_mem.hin, 0, 15));
    }
    let _ = write!(fp, "\n");
    {
        let ha = ark_mem.hadapt_mem.as_mut().unwrap();
        let _ = write!(
            fp,
            "  Maximum step increase (first step) = {}\n",
            fmt_g(ha.etamx1, 0, 15)
        );
        let _ = write!(
            fp,
            "  Step reduction factor on multiple error fails = {}\n",
            fmt_g(ha.etamxf, 0, 15)
        );
        let _ = write!(
            fp,
            "  Minimum error fails before above factor is used = {}\n",
            ha.small_nef
        );
        let _ = write!(
            fp,
            "  Step reduction factor on nonlinear convergence failure = {}\n",
            fmt_g(ha.etacf, 0, 15)
        );
        let _ = write!(fp, "  Explicit safety factor = {}\n", fmt_g(ha.cfl, 0, 15));
        let _ = write!(fp, "  Safety factor = {}\n", fmt_g(ha.safety, 0, 15));
        let _ = write!(fp, "  Growth factor = {}\n", fmt_g(ha.growth, 0, 15));
        let _ = write!(fp, "  Step growth lower bound = {}\n", fmt_g(ha.lbound, 0, 15));
        let _ = write!(fp, "  Step growth upper bound = {}\n", fmt_g(ha.ubound, 0, 15));
        if ha.expstab.is_none() {
            let _ = write!(fp, "  No explicit stability function supplied\n");
        } else {
            let _ = write!(fp, "  User provided explicit stability function\n");
        }
        if let Some(hc) = ha.hcontroller.as_mut() {
            let _ = crate::sundials_adaptcontroller::SUNAdaptController_Write(hc, fp);
        }
        if let Some(uc) = ark_mem.hadapt_mem.as_ref().unwrap().usercontrol.as_ref() {
            let _ = crate::arkode_user_controller::SUNAdaptController_Write_ARKUserControl(uc, fp);
        }
    }

    let _ = write!(fp, "  Maximum number of error test failures = {}\n", ark_mem.maxnef);
    let _ = write!(
        fp,
        "  Maximum number of convergence test failures = {}\n",
        ark_mem.maxncf
    );

    /* Call stepper routine (if provided) */
    if let Some(op) = ark_mem.step_writeparameters {
        return op(ark_mem, fp);
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetInterpolantType:

  Specifies the interpolation module (Hermite / Lagrange / none)
  to use for dense output and predictors.  May not be called after
  module initialization.
  ---------------------------------------------------------------*/
pub fn ARKodeSetInterpolantType(ark_mem: &mut ARKodeMem, itype: i32) -> i32 {
    use crate::arkode_impl::{
        ARK_ILL_INPUT, ARK_INTERP_FAIL, ARK_INTERP_HERMITE, ARK_INTERP_LAGRANGE, ARK_INTERP_NONE,
    };
    use crate::arkode_interp::{arkInterpCreate_Hermite, arkInterpCreate_Lagrange, arkInterpFree};

    /* check for legal itype input */
    if itype != ARK_INTERP_HERMITE && itype != ARK_INTERP_LAGRANGE && itype != ARK_INTERP_NONE {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "ARKodeSetInterpolantType",
            file!(),
            "Illegal interpolation type input.",
        );
        return ARK_ILL_INPUT;
    }

    /* do not change type once the module has been initialized */
    if ark_mem.initialized {
        arkProcessError(
            Some(ark_mem),
            ARK_INTERP_FAIL,
            line!(),
            "ARKodeSetInterpolantType",
            file!(),
            "Type cannot be specified after module initialization.",
        );
        return ARK_ILL_INPUT;
    }

    /* delete any existing interpolation module */
    if ark_mem.interp.is_some() {
        arkInterpFree(ark_mem);
    }

    /* create requested interpolation module, initially specifying
       the maximum possible interpolant degree (the C NULL-return
       allocation-failure paths cannot occur). */
    if itype == ARK_INTERP_HERMITE {
        ark_mem.interp = arkInterpCreate_Hermite(ark_mem, ark_mem.interp_degree);
        ark_mem.interp_type = ARK_INTERP_HERMITE;
    } else if itype == ARK_INTERP_LAGRANGE {
        ark_mem.interp = arkInterpCreate_Lagrange(ark_mem, ark_mem.interp_degree);
        ark_mem.interp_type = ARK_INTERP_LAGRANGE;
    } else {
        ark_mem.interp = None;
        ark_mem.interp_type = ARK_INTERP_NONE;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeGetNumStepSolveFails:

  Returns the current number of failed steps due to an algebraic
  solver convergence failure.
  ---------------------------------------------------------------*/
pub fn ARKodeGetNumStepSolveFails(ark_mem: &mut ARKodeMem, nncfails: &mut i64) -> i32 {
    *nncfails = ark_mem.ncfn;
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetMaxErrTestFails:

  Specifies the maximum number of error test failures during one
  step try.  A non-positive input implies a reset to the default.
  ---------------------------------------------------------------*/
pub fn ARKodeSetMaxErrTestFails(ark_mem: &mut ARKodeMem, maxnef: i32) -> i32 {
    use crate::arkode_impl::ARK_STEPPER_UNSUPPORTED;

    /* Guard against use for non-adaptive time stepper modules */
    if !ark_mem.step_supports_adaptive {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!(),
            "ARKodeSetMaxErrTestFails",
            file!(),
            "time-stepping module does not support temporal adaptivity",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* argument <= 0 sets default, otherwise set input */
    if maxnef <= 0 {
        ark_mem.maxnef = MAXNEF;
    } else {
        ark_mem.maxnef = maxnef;
    }
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeGetNumGEvals:

  Returns the current number of calls to the root function g
  ---------------------------------------------------------------*/
pub fn ARKodeGetNumGEvals(ark_mem: &mut ARKodeMem, ngevals: &mut i64) -> i32 {
    use crate::arkode_impl::ARK_MEM_NULL;

    let root_mem = match ark_mem.root_mem.as_ref() {
        Some(rm) => rm,
        None => {
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_NULL,
                line!(),
                "ARKodeGetNumGEvals",
                file!(),
                "arkode_mem = NULL illegal.",
            );
            return ARK_MEM_NULL;
        }
    };
    *ngevals = root_mem.nge;
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeGetRootInfo:

  Returns pointer to array rootsfound showing roots found
  ---------------------------------------------------------------*/
pub fn ARKodeGetRootInfo(ark_mem: &mut ARKodeMem, rootsfound: &mut [i32]) -> i32 {
    use crate::arkode_impl::ARK_MEM_NULL;

    let root_mem = match ark_mem.root_mem.as_ref() {
        Some(rm) => rm,
        None => {
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_NULL,
                line!(),
                "ARKodeGetRootInfo",
                file!(),
                "arkode_mem = NULL illegal.",
            );
            return ARK_MEM_NULL;
        }
    };
    for i in 0..root_mem.nrtfn as usize {
        rootsfound[i] = root_mem.iroots[i];
    }
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeGetWorkSpace:

  Returns integrator work space requirements
  ---------------------------------------------------------------*/
pub fn ARKodeGetWorkSpace(ark_mem: &mut ARKodeMem, lenrw: &mut i64, leniw: &mut i64) -> i32 {
    *leniw = ark_mem.liw;
    *lenrw = ark_mem.lrw;
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeGetLastStep:

  Returns the step size used on the last successful step
  ---------------------------------------------------------------*/
pub fn ARKodeGetLastStep(ark_mem: &mut ARKodeMem, hlast: &mut f64) -> i32 {
    *hlast = ark_mem.hold;
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeGetCurrentStep:

  Returns the step size to be attempted on the next step
  ---------------------------------------------------------------*/
pub fn ARKodeGetCurrentStep(ark_mem: &mut ARKodeMem, hcur: &mut f64) -> i32 {
    *hcur = ark_mem.next_h;
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetUseCompensatedSums:

  Turns compensated summation on/off in the shared and stepper
  modules.
  ---------------------------------------------------------------*/
pub fn ARKodeSetUseCompensatedSums(ark_mem: &mut ARKodeMem, onoff: bool) -> i32 {
    ark_mem.use_compensated_sums = onoff;

    /* Call stepper routine (if provided) */
    if let Some(set) = ark_mem.step_setusecompensatedsums {
        return set(ark_mem, onoff);
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeGetAccumulatedError:

  This routine returns the accumulated temporal error estimate.
  ---------------------------------------------------------------*/
#[allow(clippy::if_same_then_else)] /* C's MAX/SUM branches kept verbatim */
pub fn ARKodeGetAccumulatedError(ark_mem: &mut ARKodeMem, accum_error: &mut f64) -> i32 {
    use crate::arkode_impl::{
        arkProcessError, ARK_ACCUMERROR_AVG, ARK_ACCUMERROR_MAX, ARK_ACCUMERROR_SUM,
        ARK_STEPPER_UNSUPPORTED, ARK_WARNING,
    };

    /* Return an error if the stepper cannot accumulate temporal error */
    if !ark_mem.step_supports_adaptive {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!(),
            "ARKodeGetAccumulatedError",
            file!(),
            "time-stepping module does not support accumulated error estimation",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* Get time since last accumulated error reset */
    let time_interval = ark_mem.tcur - ark_mem.AccumErrorStart;

    /* Fill output based on error accumulation type */
    if ark_mem.AccumErrorType == ARK_ACCUMERROR_MAX {
        *accum_error = ark_mem.AccumError * ark_mem.reltol;
    } else if ark_mem.AccumErrorType == ARK_ACCUMERROR_SUM {
        *accum_error = ark_mem.AccumError * ark_mem.reltol;
    } else if ark_mem.AccumErrorType == ARK_ACCUMERROR_AVG {
        *accum_error = ark_mem.AccumError * ark_mem.reltol / time_interval;
    } else {
        arkProcessError(
            Some(ark_mem),
            ARK_WARNING,
            line!(),
            "ARKodeGetAccumulatedError",
            file!(),
            "temporal error accumulation is currently disabled",
        );
        return ARK_WARNING;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeResetAccumulatedError:

  This routine resets the accumulated temporal error estimate.
  ---------------------------------------------------------------*/
pub fn ARKodeResetAccumulatedError(ark_mem: &mut ARKodeMem) -> i32 {
    use crate::arkode_impl::{arkProcessError, ARK_STEPPER_UNSUPPORTED};

    /* Guard against use for non-adaptive time stepper modules */
    if !ark_mem.step_supports_adaptive {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!(),
            "ARKodeResetAccumulatedError",
            file!(),
            "time-stepping module does not support temporal adaptivity",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* Reset value and counter, and return */
    ark_mem.AccumErrorStart = ark_mem.tn;
    ark_mem.AccumError = ZERO;

    ARK_SUCCESS
}

/*===============================================================
  Optional input functions required by the ARKODE CLI module
  (arkode_io.c PART VI)
  ===============================================================*/

/*---------------------------------------------------------------
  ARKodeSetOrder: Specifies the method order
  ---------------------------------------------------------------*/
pub fn ARKodeSetOrder(ark_mem: &mut ARKodeMem, ord: i32) -> i32 {
    use crate::arkode_impl::{arkProcessError, ARK_STEPPER_UNSUPPORTED};

    /* Call stepper routine (if provided) */
    if let Some(step_setorder) = ark_mem.step_setorder {
        step_setorder(ark_mem, ord)
    } else {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!(),
            "ARKodeSetOrder",
            file!(),
            "time-stepping module does not support this function",
        );
        ARK_STEPPER_UNSUPPORTED
    }
}

/*---------------------------------------------------------------
  ARKodeSetPredictorMethod: implicit-solution predictor method.
  ---------------------------------------------------------------*/
pub fn ARKodeSetPredictorMethod(ark_mem: &mut ARKodeMem, pred_method: i32) -> i32 {
    use crate::arkode_impl::{arkProcessError, ARK_ILL_INPUT, ARK_INTERP_NONE,
        ARK_STEPPER_UNSUPPORTED};

    /* Guard against use for time steppers that do not need an algebraic solver */
    if !ark_mem.step_supports_implicit {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!(),
            "ARKodeSetPredictorMethod",
            file!(),
            "time-stepping module does not require an algebraic solver",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* Higher-order predictors require interpolation */
    if ark_mem.interp_type == ARK_INTERP_NONE && pred_method != 0 {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "ARKodeSetPredictorMethod",
            file!(),
            "Non-trival predictors require an interpolation module",
        );
        return ARK_ILL_INPUT;
    }

    /* Call stepper routine (if provided) */
    if let Some(step_setpredictormethod) = ark_mem.step_setpredictormethod {
        step_setpredictormethod(ark_mem, pred_method)
    } else {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!(),
            "ARKodeSetPredictorMethod",
            file!(),
            "time-stepping module does not support this function",
        );
        ARK_STEPPER_UNSUPPORTED
    }
}

/*---------------------------------------------------------------
  ARKodeSetMaxHnilWarns: max warnings of t+h==t.
  ---------------------------------------------------------------*/
pub fn ARKodeSetMaxHnilWarns(ark_mem: &mut ARKodeMem, mxhnil: i32) -> i32 {
    use crate::arkode_impl::{arkProcessError, ARK_STEPPER_UNSUPPORTED};

    /* Guard against use for non-adaptive time stepper modules */
    if !ark_mem.step_supports_adaptive {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!(),
            "ARKodeSetMaxHnilWarns",
            file!(),
            "time-stepping module does not support temporal adaptivity",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* Passing mxhnil=0 sets the default, otherwise use input. */
    if mxhnil == 0 {
        ark_mem.mxhnil = 10;
    } else {
        ark_mem.mxhnil = mxhnil;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetInterpolateStopTime: interpolate at tstop.
  ---------------------------------------------------------------*/
pub fn ARKodeSetInterpolateStopTime(ark_mem: &mut ARKodeMem, interp: bool) -> i32 {
    ark_mem.tstopinterp = interp;
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetMaxNumConstrFails: max constraint failures per step.
  ---------------------------------------------------------------*/
pub fn ARKodeSetMaxNumConstrFails(ark_mem: &mut ARKodeMem, maxfails: i32) -> i32 {
    use crate::arkode_impl::{arkProcessError, ARK_STEPPER_UNSUPPORTED, MAXCONSTRFAILS};

    /* Guard against use for non-adaptive time stepper modules */
    if !ark_mem.step_supports_adaptive {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!(),
            "ARKodeSetMaxNumConstrFails",
            file!(),
            "time-stepping module does not support temporal adaptivity",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* Passing maxfails = 0 sets the default, otherwise set to input */
    if maxfails <= 0 {
        ark_mem.maxconstrfails = MAXCONSTRFAILS;
    } else {
        ark_mem.maxconstrfails = maxfails;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetAdaptivityAdjustment: controller order adjustment.
  ---------------------------------------------------------------*/
pub fn ARKodeSetAdaptivityAdjustment(ark_mem: &mut ARKodeMem, adjust: i32) -> i32 {
    use crate::arkode_impl::{arkProcessError, ARK_STEPPER_UNSUPPORTED};

    /* Guard against use for non-adaptive time stepper modules */
    if !ark_mem.step_supports_adaptive {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!(),
            "ARKodeSetAdaptivityAdjustment",
            file!(),
            "time-stepping module does not support temporal adaptivity",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* store requested adjustment */
    ark_mem.hadapt_mem.as_mut().unwrap().adjust = adjust;
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetSmallNumEFails: num error fails before ETAMXF enforced.
  ---------------------------------------------------------------*/
pub fn ARKodeSetSmallNumEFails(ark_mem: &mut ARKodeMem, small_nef: i32) -> i32 {
    use crate::arkode_adapt_impl::SMALL_NEF;
    use crate::arkode_impl::{arkProcessError, ARK_STEPPER_UNSUPPORTED};

    /* Guard against use for non-adaptive time stepper modules */
    if !ark_mem.step_supports_adaptive {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!(),
            "ARKodeSetSmallNumEFails",
            file!(),
            "time-stepping module does not support temporal adaptivity",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* if argument legal set it, otherwise set default */
    let hadapt_mem = ark_mem.hadapt_mem.as_mut().unwrap();
    if small_nef <= 0 {
        hadapt_mem.small_nef = SMALL_NEF;
    } else {
        hadapt_mem.small_nef = small_nef;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetMaxConvFails: max convergence failures per step.
  ---------------------------------------------------------------*/
pub fn ARKodeSetMaxConvFails(ark_mem: &mut ARKodeMem, maxncf: i32) -> i32 {
    use crate::arkode_impl::{arkProcessError, ARK_STEPPER_UNSUPPORTED, MAXNCF};

    /* Guard against use for time steppers that do not need an algebraic solver */
    if !ark_mem.step_supports_implicit {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!(),
            "ARKodeSetMaxConvFails",
            file!(),
            "time-stepping module does not require an algebraic solver",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* argument <= 0 sets default, otherwise set input */
    if maxncf <= 0 {
        ark_mem.maxncf = MAXNCF;
    } else {
        ark_mem.maxncf = maxncf;
    }
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetMinStep: minimum absolute step size.
  ---------------------------------------------------------------*/
pub fn ARKodeSetMinStep(ark_mem: &mut ARKodeMem, hmin: f64) -> i32 {
    use crate::arkode_impl::{arkProcessError, ARK_ILL_INPUT, ARK_STEPPER_UNSUPPORTED,
        MSG_ARK_BAD_HMIN_HMAX, ONE};

    /* Guard against use for non-adaptive time stepper modules */
    if !ark_mem.step_supports_adaptive {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!(),
            "ARKodeSetMinStep",
            file!(),
            "time-stepping module does not support temporal adaptivity",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* Passing a value <= 0 sets hmin = 0 */
    if hmin <= ZERO {
        ark_mem.hmin = ZERO;
        return ARK_SUCCESS;
    }

    /* check that hmin and hmax are agreeable */
    if hmin * ark_mem.hmax_inv > ONE {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "ARKodeSetMinStep",
            file!(),
            MSG_ARK_BAD_HMIN_HMAX,
        );
        return ARK_ILL_INPUT;
    }

    /* set the value */
    ark_mem.hmin = hmin;

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetMaxStep: maximum absolute step size.
  ---------------------------------------------------------------*/
pub fn ARKodeSetMaxStep(ark_mem: &mut ARKodeMem, hmax: f64) -> i32 {
    use crate::arkode_impl::{arkProcessError, ARK_ILL_INPUT, ARK_STEPPER_UNSUPPORTED,
        MSG_ARK_BAD_HMIN_HMAX, ONE};

    /* Guard against use for non-adaptive time stepper modules */
    if !ark_mem.step_supports_adaptive {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!(),
            "ARKodeSetMaxStep",
            file!(),
            "time-stepping module does not support temporal adaptivity",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* Passing a value <= 0 sets hmax = infinity */
    if hmax <= ZERO {
        ark_mem.hmax_inv = ZERO;
        return ARK_SUCCESS;
    }

    /* check that hmax and hmin are agreeable */
    let hmax_inv = ONE / hmax;
    if hmax_inv * ark_mem.hmin > ONE {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "ARKodeSetMaxStep",
            file!(),
            MSG_ARK_BAD_HMIN_HMAX,
        );
        return ARK_ILL_INPUT;
    }

    /* set the value */
    ark_mem.hmax_inv = hmax_inv;

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetCFLFraction: explicit stability safety factor.
  ---------------------------------------------------------------*/
pub fn ARKodeSetCFLFraction(ark_mem: &mut ARKodeMem, cfl_frac: f64) -> i32 {
    use crate::arkode_adapt_impl::CFLFAC;
    use crate::arkode_impl::{arkProcessError, ARK_STEPPER_UNSUPPORTED};

    /* Guard against use for non-adaptive time stepper modules */
    if !ark_mem.step_supports_adaptive {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!(),
            "ARKodeSetCFLFraction",
            file!(),
            "time-stepping module does not support temporal adaptivity",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* set positive-valued parameters, otherwise set default */
    let hadapt_mem = ark_mem.hadapt_mem.as_mut().unwrap();
    if cfl_frac <= ZERO {
        hadapt_mem.cfl = CFLFAC;
    } else {
        hadapt_mem.cfl = cfl_frac;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetSafetyFactor: step adaptivity safety factor.
  ---------------------------------------------------------------*/
pub fn ARKodeSetSafetyFactor(ark_mem: &mut ARKodeMem, safety: f64) -> i32 {
    use crate::arkode_adapt_impl::SAFETY;
    use crate::arkode_impl::{arkProcessError, ARK_ILL_INPUT, ARK_STEPPER_UNSUPPORTED, ONE};

    /* Guard against use for non-adaptive time stepper modules */
    if !ark_mem.step_supports_adaptive {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!(),
            "ARKodeSetSafetyFactor",
            file!(),
            "time-stepping module does not support temporal adaptivity",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* check for allowable parameters */
    if safety > ONE {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "ARKodeSetSafetyFactor",
            file!(),
            "Illegal safety factor",
        );
        return ARK_ILL_INPUT;
    }

    /* set positive-valued parameters, otherwise set default */
    let hadapt_mem = ark_mem.hadapt_mem.as_mut().unwrap();
    if safety <= ZERO {
        hadapt_mem.safety = SAFETY;
    } else {
        hadapt_mem.safety = safety;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetErrorBias: temporal error bias factor.
  ---------------------------------------------------------------*/
pub fn ARKodeSetErrorBias(ark_mem: &mut ARKodeMem, bias: f64) -> i32 {
    use crate::arkode_impl::{arkProcessError, ARK_CONTROLLER_ERR, ARK_MEM_NULL,
        ARK_STEPPER_UNSUPPORTED, ONE};
    use crate::sundials_adaptcontroller::SUNAdaptController_SetErrorBias;

    /* Guard against use for non-adaptive time stepper modules */
    if !ark_mem.step_supports_adaptive {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!(),
            "ARKodeSetErrorBias",
            file!(),
            "time-stepping module does not support temporal adaptivity",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* Return an error if there is not a current SUNAdaptController */
    if ark_mem.hadapt_mem.as_ref().unwrap().hcontroller.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_NULL,
            line!(),
            "ARKodeSetErrorBias",
            file!(),
            "SUNAdaptController NULL -- must be set before setting the error bias",
        );
        return ARK_MEM_NULL;
    }

    /* set allowed value, otherwise set default */
    let hcontroller = ark_mem
        .hadapt_mem
        .as_mut()
        .unwrap()
        .hcontroller
        .as_mut()
        .unwrap();
    let retval = if bias < ONE {
        SUNAdaptController_SetErrorBias(hcontroller, -ONE)
    } else {
        SUNAdaptController_SetErrorBias(hcontroller, bias)
    };
    if retval != crate::sundials_errors::SUN_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            ARK_CONTROLLER_ERR,
            line!(),
            "ARKodeSetErrorBias",
            file!(),
            "SUNAdaptController_SetErrorBias failure",
        );
        return ARK_CONTROLLER_ERR;
    }
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetMaxGrowth: max step growth factor.
  ---------------------------------------------------------------*/
pub fn ARKodeSetMaxGrowth(ark_mem: &mut ARKodeMem, mx_growth: f64) -> i32 {
    use crate::arkode_adapt_impl::GROWTH;
    use crate::arkode_impl::{arkProcessError, ARK_STEPPER_UNSUPPORTED, ONE};

    /* Guard against use for non-adaptive time stepper modules */
    if !ark_mem.step_supports_adaptive {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!(),
            "ARKodeSetMaxGrowth",
            file!(),
            "time-stepping module does not support temporal adaptivity",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* set allowed value, otherwise set default */
    let hadapt_mem = ark_mem.hadapt_mem.as_mut().unwrap();
    if mx_growth <= ONE {
        hadapt_mem.growth = GROWTH;
    } else {
        hadapt_mem.growth = mx_growth;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetMinReduction: min step reduction factor.
  ---------------------------------------------------------------*/
pub fn ARKodeSetMinReduction(ark_mem: &mut ARKodeMem, eta_min: f64) -> i32 {
    use crate::arkode_adapt_impl::ETAMIN;
    use crate::arkode_impl::{arkProcessError, ARK_STEPPER_UNSUPPORTED, ONE};

    /* Guard against use for non-adaptive time stepper modules */
    if !ark_mem.step_supports_adaptive {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!(),
            "ARKodeSetMinReduction",
            file!(),
            "time-stepping module does not support temporal adaptivity",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* set allowed value, otherwise set default */
    let hadapt_mem = ark_mem.hadapt_mem.as_mut().unwrap();
    if eta_min >= ONE || eta_min <= ZERO {
        hadapt_mem.etamin = ETAMIN;
    } else {
        hadapt_mem.etamin = eta_min;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetMaxFirstGrowth: max first-step growth factor.
  ---------------------------------------------------------------*/
pub fn ARKodeSetMaxFirstGrowth(ark_mem: &mut ARKodeMem, etamx1: f64) -> i32 {
    use crate::arkode_adapt_impl::ETAMX1;
    use crate::arkode_impl::{arkProcessError, ARK_STEPPER_UNSUPPORTED, ONE};

    /* Guard against use for non-adaptive time stepper modules */
    if !ark_mem.step_supports_adaptive {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!(),
            "ARKodeSetMaxFirstGrowth",
            file!(),
            "time-stepping module does not support temporal adaptivity",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* if argument legal set it, otherwise set default */
    let hadapt_mem = ark_mem.hadapt_mem.as_mut().unwrap();
    if etamx1 <= ONE {
        hadapt_mem.etamx1 = ETAMX1;
    } else {
        hadapt_mem.etamx1 = etamx1;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetMaxEFailGrowth: max error-fail growth factor.
  ---------------------------------------------------------------*/
pub fn ARKodeSetMaxEFailGrowth(ark_mem: &mut ARKodeMem, etamxf: f64) -> i32 {
    use crate::arkode_adapt_impl::ETAMXF;
    use crate::arkode_impl::{arkProcessError, ARK_STEPPER_UNSUPPORTED, ONE};

    /* Guard against use for non-adaptive time stepper modules */
    if !ark_mem.step_supports_adaptive {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!(),
            "ARKodeSetMaxEFailGrowth",
            file!(),
            "time-stepping module does not support temporal adaptivity",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* if argument legal set it, otherwise set default */
    let hadapt_mem = ark_mem.hadapt_mem.as_mut().unwrap();
    if etamxf <= ZERO || etamxf > ONE {
        hadapt_mem.etamxf = ETAMXF;
    } else {
        hadapt_mem.etamxf = etamxf;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetMaxCFailGrowth: max convergence-fail growth factor.
  ---------------------------------------------------------------*/
pub fn ARKodeSetMaxCFailGrowth(ark_mem: &mut ARKodeMem, etacf: f64) -> i32 {
    use crate::arkode_adapt_impl::ETACF;
    use crate::arkode_impl::{arkProcessError, ARK_STEPPER_UNSUPPORTED, ONE};

    /* Guard against use for non-adaptive time stepper modules */
    if !ark_mem.step_supports_adaptive {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!(),
            "ARKodeSetMaxCFailGrowth",
            file!(),
            "time-stepping module does not support temporal adaptivity",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* if argument legal set it, otherwise set default */
    let hadapt_mem = ark_mem.hadapt_mem.as_mut().unwrap();
    if etacf <= ZERO || etacf > ONE {
        hadapt_mem.etacf = ETACF;
    } else {
        hadapt_mem.etacf = etacf;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetFixedStepBounds: no-change step growth interval.
  ---------------------------------------------------------------*/
pub fn ARKodeSetFixedStepBounds(ark_mem: &mut ARKodeMem, lb: f64, ub: f64) -> i32 {
    use crate::arkode_adapt_impl::{HFIXED_LB, HFIXED_UB};
    use crate::arkode_impl::{arkProcessError, ARK_STEPPER_UNSUPPORTED, ONE};

    /* Guard against use for non-adaptive time stepper modules */
    if !ark_mem.step_supports_adaptive {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!(),
            "ARKodeSetFixedStepBounds",
            file!(),
            "time-stepping module does not support temporal adaptivity",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* set allowable interval, otherwise set defaults */
    let hadapt_mem = ark_mem.hadapt_mem.as_mut().unwrap();
    if lb <= ONE && ub >= ONE {
        hadapt_mem.lbound = lb;
        hadapt_mem.ubound = ub;
    } else {
        hadapt_mem.lbound = HFIXED_LB;
        hadapt_mem.ubound = HFIXED_UB;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeClearStopTime: disables the stop time.
  ---------------------------------------------------------------*/
pub fn ARKodeClearStopTime(ark_mem: &mut ARKodeMem) -> i32 {
    ark_mem.tstopset = false;
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetNoInactiveRootWarn: disables the inactive-root warning.
  ---------------------------------------------------------------*/
pub fn ARKodeSetNoInactiveRootWarn(ark_mem: &mut ARKodeMem) -> i32 {
    use crate::arkode_impl::{arkProcessError, ARK_MEM_NULL, MSG_ARK_NO_ROOT};

    if ark_mem.root_mem.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_NULL,
            line!(),
            "ARKodeSetNoInactiveRootWarn",
            file!(),
            MSG_ARK_NO_ROOT,
        );
        return ARK_MEM_NULL;
    }
    ark_mem.root_mem.as_mut().unwrap().mxgnull = 0;
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetAccumulatedErrorType: enables accumulated-error
  estimation.
  ---------------------------------------------------------------*/
pub fn ARKodeSetAccumulatedErrorType(
    ark_mem: &mut ARKodeMem,
    accum_type: crate::arkode_impl::ARKAccumError,
) -> i32 {
    let retval = ARKodeResetAccumulatedError(ark_mem);
    if retval != ARK_SUCCESS {
        return retval;
    }
    ark_mem.AccumErrorType = accum_type;
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeSetConstraints:

  Activates or deactivates inequality constraint checking.
  ---------------------------------------------------------------*/
pub fn ARKodeSetConstraints(
    ark_mem: &mut ARKodeMem,
    constraints: Option<&crate::nvector_serial::NVector>,
) -> i32 {
    use crate::arkode_impl::{arkProcessError, ARK_ILL_INPUT, ARK_STEPPER_UNSUPPORTED,
        HALF, MSG_ARK_BAD_CONSTR, ONE};
    use crate::nvector_serial::{N_VMaxNorm, N_VScale};

    /* Guard against use for non-adaptive time stepper modules */
    if !ark_mem.step_supports_adaptive && constraints.is_some() {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!(),
            "ARKodeSetConstraints",
            file!(),
            "time-stepping module does not support temporal adaptivity",
        );
        return ARK_STEPPER_UNSUPPORTED;
    }

    /* If there are no constraints, destroy data structures */
    let constraints = match constraints {
        None => {
            if ark_mem.constraints.take().is_some() {
                ark_mem.lrw -= ark_mem.lrw1;
                ark_mem.liw -= ark_mem.liw1;
            }
            return ARK_SUCCESS;
        }
        Some(c) => c,
    };

    /* (C also tests that the required vector ops are defined; the
    serial NVector provides them all) */

    /* Check the constraints vector */
    let temptest = N_VMaxNorm(constraints);
    if temptest > 2.5 || temptest < HALF {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "ARKodeSetConstraints",
            file!(),
            MSG_ARK_BAD_CONSTR,
        );
        return ARK_ILL_INPUT;
    }

    /* Allocate the internal constraints vector (if necessary) */
    if ark_mem.constraints.is_none() {
        ark_mem.constraints = Some(crate::nvector_serial::NVector::new(constraints.data.len()));
        ark_mem.lrw += ark_mem.lrw1;
        ark_mem.liw += ark_mem.liw1;
    }

    /* Load the constraints vector */
    N_VScale(ONE, constraints, ark_mem.constraints.as_mut().unwrap());

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  ARKodeGetNumConstrFails:

  Returns the current number of constraint fails
  ---------------------------------------------------------------*/
pub fn ARKodeGetNumConstrFails(ark_mem: &mut ARKodeMem, nconstrfails: &mut i64) -> i32 {
    *nconstrfails = ark_mem.nconstrfails;
    ARK_SUCCESS
}
