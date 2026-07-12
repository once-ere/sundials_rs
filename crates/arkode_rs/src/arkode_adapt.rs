/* -----------------------------------------------------------------
 * Translation of sundials-7.7.0/src/arkode/arkode_adapt.c.
 *
 * arkAdapt in C receives (ark_mem, hadapt_mem, ycur, ...) where the
 * latter two are always ark_mem->hadapt_mem and ark_mem->ycur (one
 * caller: arkCheckTemporalError); the Rust signature drops them and
 * takes/puts back internally per the crate aliasing conventions
 * (ARCHITECTURE Addendum C.1). The SUNLogDebug lines are compiled
 * out as in the rest of the workspace.
 * -----------------------------------------------------------------*/

use crate::arkode_adapt_impl::ARKodeHAdaptMem;
use crate::arkode_impl::{
    arkProcessError, ARKodeMem, ARK_CONTROLLER_ERR, ARK_ILL_INPUT, ARK_SUCCESS, ONE, ONEMSM,
    ONEPSM, ZERO,
};
use crate::sundials_adaptcontroller::{
    SUNAdaptController_EstimateStep, SUNAdaptController_Write,
};
use crate::sundials_errors::SUN_SUCCESS;
use crate::sundials_math::{SUNMAX, SUNMIN, SUNRabs, SUNRcopysign};
use crate::sundials_utils::fmt_g;

/*---------------------------------------------------------------
  arkAdaptInit:

  This routine creates and sets default values in an
  ARKodeHAdaptMem structure.  This returns a non-NULL structure
  if no errors occurred, or a NULL value otherwise.
  ---------------------------------------------------------------*/
pub fn arkAdaptInit() -> ARKodeHAdaptMem {
    /* initialize values (default parameters are set in
    ARKodeSetDefaults); C memsets the structure to zero */
    ARKodeHAdaptMem {
        etamax: 0.0,
        etamx1: 0.0,
        etamxf: 0.0,
        etamin: 0.0,
        small_nef: 0,
        etacf: 0.0,
        cfl: 0.0,
        safety: 0.0,
        growth: 0.0,
        lbound: 0.0,
        ubound: 0.0,
        p: 0,
        q: 0,
        pq: 0,
        adjust: 0,
        hcontroller: None,
        owncontroller: false,
        expstab: None,
        estab_data: None,
        nst_acc: 0,
        nst_exp: 0,
    }
}

/*---------------------------------------------------------------
  arkPrintAdaptMem

  This routine outputs the time step adaptivity memory structure
  to a specified file pointer.
  ---------------------------------------------------------------*/
pub fn arkPrintAdaptMem(hadapt_mem: &ARKodeHAdaptMem, outfile: &mut dyn std::io::Write) {
    let _ = write!(outfile, "ark_hadapt: etamax = {}\n", fmt_g(hadapt_mem.etamax, 0, 15));
    let _ = write!(outfile, "ark_hadapt: etamx1 = {}\n", fmt_g(hadapt_mem.etamx1, 0, 15));
    let _ = write!(outfile, "ark_hadapt: etamxf = {}\n", fmt_g(hadapt_mem.etamxf, 0, 15));
    let _ = write!(outfile, "ark_hadapt: etamin = {}\n", fmt_g(hadapt_mem.etamin, 0, 15));
    let _ = write!(outfile, "ark_hadapt: small_nef = {}\n", hadapt_mem.small_nef);
    let _ = write!(outfile, "ark_hadapt: etacf = {}\n", fmt_g(hadapt_mem.etacf, 0, 15));
    let _ = write!(outfile, "ark_hadapt: cfl = {}\n", fmt_g(hadapt_mem.cfl, 0, 15));
    let _ = write!(outfile, "ark_hadapt: safety = {}\n", fmt_g(hadapt_mem.safety, 0, 15));
    let _ = write!(outfile, "ark_hadapt: growth = {}\n", fmt_g(hadapt_mem.growth, 0, 15));
    let _ = write!(outfile, "ark_hadapt: lbound = {}\n", fmt_g(hadapt_mem.lbound, 0, 15));
    let _ = write!(outfile, "ark_hadapt: ubound = {}\n", fmt_g(hadapt_mem.ubound, 0, 15));
    let _ = write!(outfile, "ark_hadapt: nst_acc = {}\n", hadapt_mem.nst_acc);
    let _ = write!(outfile, "ark_hadapt: nst_exp = {}\n", hadapt_mem.nst_exp);
    let _ = write!(outfile, "ark_hadapt: pq = {}\n", hadapt_mem.pq);
    let _ = write!(outfile, "ark_hadapt: p = {}\n", hadapt_mem.p);
    let _ = write!(outfile, "ark_hadapt: q = {}\n", hadapt_mem.q);
    let _ = write!(outfile, "ark_hadapt: adjust = {}\n", hadapt_mem.adjust);
    if hadapt_mem.expstab.is_none() {
        let _ = write!(outfile, "  ark_hadapt: No explicit stability function supplied\n");
    } else {
        let _ = write!(outfile, "  ark_hadapt: User provided explicit stability function\n");
        let _ = write!(
            outfile,
            "  ark_hadapt: stability function data pointer = {:p}\n",
            &hadapt_mem.estab_data
        );
    }
    if let Some(c) = hadapt_mem.hcontroller.as_ref() {
        let _ = SUNAdaptController_Write(c, outfile);
    }
}

/*---------------------------------------------------------------
  arkAdapt is the time step adaptivity wrapper function.  This
  computes and sets the value of ark_eta inside of the ARKodeMem
  data structure.
  ---------------------------------------------------------------*/
pub fn arkAdapt(ark_mem: &mut ARKodeMem, tcur: f64, hcur: f64, dsm: f64) -> i32 {
    let mut hadapt_mem = ark_mem.hadapt_mem.take().unwrap();
    let ret = arkAdapt_inner(ark_mem, &mut hadapt_mem, tcur, hcur, dsm);
    ark_mem.hadapt_mem = Some(hadapt_mem);
    ret
}

fn arkAdapt_inner(
    ark_mem: &mut ARKodeMem,
    hadapt_mem: &mut ARKodeHAdaptMem,
    tcur: f64,
    hcur: f64,
    dsm: f64,
) -> i32 {
    let mut h_acc = 0.0;

    /* Return with no stepsize adjustment if the controller is NULL */
    if hadapt_mem.hcontroller.is_none() {
        ark_mem.eta = ONE;
        return ARK_SUCCESS;
    }

    /* Request error-based step size from adaptivity controller */
    let controller_order = if hadapt_mem.pq == 0 {
        hadapt_mem.p + hadapt_mem.adjust
    } else if hadapt_mem.pq == 1 {
        hadapt_mem.q + hadapt_mem.adjust
    } else {
        std::cmp::min(hadapt_mem.p, hadapt_mem.q) + hadapt_mem.adjust
    };
    /* (an MRI-H-TOL controller carries C's SUNAdaptController_MRIStep
       wrapper semantics: dispatch through the MRIStep step memory) */
    let hc = hadapt_mem.hcontroller.as_mut().unwrap();
    let mut retval = if crate::sundials_adaptcontroller::SUNAdaptController_GetType(hc)
        == crate::sundials_adaptcontroller::SUNAdaptController_Type::SUN_ADAPTCONTROLLER_MRI_H_TOL
    {
        crate::arkode_mristep_controller::SUNAdaptController_EstimateStep_MRIStep(
            ark_mem,
            hc,
            hcur,
            controller_order,
            dsm,
            &mut h_acc,
        )
    } else {
        SUNAdaptController_EstimateStep(hc, hcur, controller_order, dsm, &mut h_acc)
    };
    if retval != SUN_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            ARK_CONTROLLER_ERR,
            line!(),
            "arkAdapt",
            file!(),
            "SUNAdaptController_EstimateStep failure.",
        );
        return ARK_CONTROLLER_ERR;
    }

    /* enforce safety factors */
    h_acc *= hadapt_mem.safety;

    /* enforce maximum bound on time step growth */
    h_acc = SUNMIN(SUNRabs(h_acc), SUNRabs(hadapt_mem.etamax * hcur));

    /* enforce minimum bound time step reduction */
    h_acc = SUNMAX(h_acc, SUNRabs(hadapt_mem.etamin * hcur));

    if let Some(expstab) = hadapt_mem.expstab {
        let mut h_cfl = ZERO;
        /* C passes ark_mem->ycur (only caller: arkCheckTemporalError) */
        let ycur = std::mem::replace(&mut ark_mem.ycur, crate::nvector_serial::NVector::new(0));
        retval = expstab(&ycur, tcur, &mut h_cfl, &mut hadapt_mem.estab_data);
        ark_mem.ycur = ycur;
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "arkAdapt",
                file!(),
                "Error in explicit stability function.",
            );
            return ARK_ILL_INPUT;
        }

        h_cfl *= hadapt_mem.cfl;

        if h_cfl > ZERO && h_cfl < h_acc {
            hadapt_mem.nst_exp += 1;
            h_acc = h_cfl;
        } else {
            hadapt_mem.nst_acc += 1;
        }
    } else {
        hadapt_mem.nst_acc += 1;
    }

    /* enforce adaptivity bounds to retain Jacobian/preconditioner accuracy */
    if dsm <= ONE {
        if (h_acc > SUNRabs(hcur * hadapt_mem.lbound * ONEMSM))
            && (h_acc < SUNRabs(hcur * hadapt_mem.ubound * ONEPSM))
        {
            h_acc = hcur;
        }
    }
    h_acc = SUNRcopysign(h_acc, hcur);

    /* set basic value of ark_eta */
    ark_mem.eta = h_acc / hcur;

    /* enforce minimum time step size */
    ark_mem.eta = SUNMAX(ark_mem.eta, ark_mem.hmin / SUNRabs(hcur));

    /* enforce maximum time step size */
    ark_mem.eta /= SUNMAX(ONE, SUNRabs(hcur) * ark_mem.hmax_inv * ark_mem.eta);

    retval
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arkode_adapt_impl::{ETAMIN, GROWTH, SAFETY};
    use crate::nvector_serial::NVector;
    use crate::sunadaptcontroller_soderlind::SUNAdaptController_PID;
    use crate::sundials_context::SUNContext_Create;

    #[test]
    fn adapt_no_controller_sets_eta_one() {
        let mut ark_mem = ARKodeMem::default();
        ark_mem.hadapt_mem = Some(arkAdaptInit());
        assert_eq!(arkAdapt(&mut ark_mem, 0.0, 0.1, 0.5), ARK_SUCCESS);
        assert_eq!(ark_mem.eta, ONE);
    }

    #[test]
    fn adapt_with_pid_controller() {
        let _ctx = SUNContext_Create();
        let mut ark_mem = ARKodeMem::default();
        ark_mem.ycur = NVector::new(1);
        ark_mem.hmin = 0.0;
        ark_mem.hmax_inv = 0.0;
        let mut hm = arkAdaptInit();
        hm.etamax = GROWTH;
        hm.etamin = ETAMIN;
        hm.safety = SAFETY;
        hm.lbound = 1.0;
        hm.ubound = 1.5;
        hm.p = 2;
        hm.q = 3;
        hm.pq = 0;
        hm.hcontroller = Some(SUNAdaptController_PID());
        ark_mem.hadapt_mem = Some(hm);

        /* dsm > 1 (failed step): controller shrinks the step */
        assert_eq!(arkAdapt(&mut ark_mem, 0.0, 0.1, 4.0), ARK_SUCCESS);
        assert!(ark_mem.eta < ONE, "eta = {}", ark_mem.eta);
        assert!(ark_mem.eta >= ETAMIN);
        let hm = ark_mem.hadapt_mem.as_ref().unwrap();
        assert_eq!(hm.nst_acc, 1);
        assert_eq!(hm.nst_exp, 0);
    }
}
