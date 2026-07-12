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

    /* (relaxation stats: module pending; relax_enabled cannot be set) */

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
ark_io_dispatch_implicit!(ARKodeSetPredictorMethod, step_setpredictormethod, method: i32);
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
