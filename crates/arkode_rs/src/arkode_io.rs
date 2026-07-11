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
    ARKodeMem, ARK_SS, ARK_SUCCESS, MAXCONSTRFAILS, MAXNCF, MAXNEF, MXHNIL, MXSTEP_DEFAULT, ZERO,
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
    (delete if owned, and then nullify pointer) */
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
