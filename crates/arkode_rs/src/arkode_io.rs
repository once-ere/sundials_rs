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
