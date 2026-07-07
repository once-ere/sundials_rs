/* -----------------------------------------------------------------
 * Translated from src/ida/ida_ic.c (IDA 7.7.0).
 * IC (consistent initial condition) calculation for IDA; independent
 * of the linear solver in use.
 *
 * The C code aliases four IDACalcIC work "vectors" onto existing
 * IDAMem fields (see ida_impl.rs): ynew = tempv2, ypnew = ee,
 * delnew = phi[2], dtemp = phi[3], and mc = ee (used before ypnew).
 * This port references those underlying fields directly. Vectors are
 * detached with std::mem::take only where the linear-solver dispatch
 * needs &mut IDAMem alongside borrowed vector arguments.
 * -----------------------------------------------------------------*/
use crate::ida::{IDAEwtSet, IDAInitialSetup, IDAWrmsNorm};
use crate::ida_impl::*;
use crate::ida_ls::{idaLsSetup, idaLsSolve};
use crate::ida_nls::ida_has_lsetup;
use crate::nvector_serial::{
    NVector, N_VConstrMask, N_VLinearSum, N_VMin, N_VMinQuotient, N_VProd, N_VScale,
};
use crate::sundials_math::SUNRabs;
use crate::sundials_types::*;

/* Private Constants */
const ZERO: f64 = 0.0;
const HALF: f64 = 0.5;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;
const PT99: f64 = 0.99;
const PT1: f64 = 0.1;
const PT001: f64 = 0.001;

/* IDACalcIC control constants */
const ICRATEMAX: f64 = 0.9; /* max. Newton conv. rate         */
const ALPHALS: f64 = 0.0001; /* alpha in linesearch conv. test */

/* Return values for lower level routines used by IDACalcIC */
const IC_FAIL_RECOV: i32 = 1;
const IC_CONSTR_FAILED: i32 = 2;
const IC_LINESRCH_FAILED: i32 = 3;
const IC_CONV_FAIL: i32 = 4;
const IC_SLOW_CONVRG: i32 = 5;

/*
 * =================================================================
 * EXPORTED FUNCTIONS IMPLEMENTATION
 * =================================================================
 */

/*
 * IDACalcIC computes consistent initial conditions, given the user's
 * initial guess for unknown components of yy0 and/or yp0.
 */
pub fn IDACalcIC(ida_mem: &mut IDAMem, icopt: i32, tout1: f64) -> i32 {
    /* Check if problem was malloc'ed */
    if !ida_mem.ida_MallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_MALLOC, line!(), "IDACalcIC", file!(),
                        MSG_NO_MALLOC);
        return IDA_NO_MALLOC;
    }

    /* Check inputs to IDA for correctness and consistency */
    let ier = IDAInitialSetup(ida_mem);
    if ier != IDA_SUCCESS {
        return IDA_ILL_INPUT;
    }
    ida_mem.ida_SetupDone = SUNTRUE;

    /* Check legality of input arguments, and set IDA memory copies. */
    if icopt != IDA_YA_YDP_INIT && icopt != IDA_Y_INIT {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDACalcIC", file!(),
                        MSG_IC_BAD_ICOPT);
        return IDA_ILL_INPUT;
    }
    ida_mem.ida_icopt = icopt;

    if icopt == IDA_YA_YDP_INIT && !ida_mem.ida_idMallocDone {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDACalcIC", file!(),
                        MSG_IC_MISSING_ID);
        return IDA_ILL_INPUT;
    }

    let tdist = SUNRabs(tout1 - ida_mem.ida_tn);
    let troundoff = TWO * ida_mem.ida_uround * (SUNRabs(ida_mem.ida_tn) + SUNRabs(tout1));
    if tdist < troundoff {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDACalcIC", file!(),
                        MSG_IC_TOO_CLOSE);
        return IDA_ILL_INPUT;
    }

    /* Allocate space and initialize temporary vectors */
    ida_mem.ida_yy0 = NVector::new(ida_mem.ida_ee.len());
    ida_mem.ida_yp0 = NVector::new(ida_mem.ida_ee.len());
    ida_mem.ida_t0 = ida_mem.ida_tn;
    N_VScale(ONE, &ida_mem.ida_phi[0], &mut ida_mem.ida_yy0);
    N_VScale(ONE, &ida_mem.ida_phi[1], &mut ida_mem.ida_yp0);

    /* For use in the IDA_YA_YP_INIT case, set sysindex and tscale. */
    ida_mem.ida_sysindex = 1;
    ida_mem.ida_tscale = tdist;
    if icopt == IDA_YA_YDP_INIT {
        let minid = N_VMin(&ida_mem.ida_id);
        if minid < ZERO {
            IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDACalcIC", file!(),
                            MSG_IC_BAD_ID);
            return IDA_ILL_INPUT;
        }
        if minid > HALF {
            ida_mem.ida_sysindex = 0;
        }
    }

    /* Set the test constant in the Newton convergence test */
    ida_mem.ida_epsNewt = ida_mem.ida_epiccon;

    /* Initializations: cjratio = 1 (for use in direct linear solvers);
       set nbacktr = 0; */
    ida_mem.ida_cjratio = ONE;
    ida_mem.ida_nbacktr = 0;

    /* Set hic, hh, cj, and mxnh. */
    let mut hic = PT001 * tdist;
    let ypnorm = IDAWrmsNorm(ida_mem, &ida_mem.ida_yp0, &ida_mem.ida_ewt, ida_mem.ida_suppressalg);
    if ypnorm > HALF / hic {
        hic = HALF / ypnorm;
    }
    if tout1 < ida_mem.ida_tn {
        hic = -hic;
    }
    ida_mem.ida_hh = hic;
    let mxnh;
    if icopt == IDA_YA_YDP_INIT {
        ida_mem.ida_cj = ONE / hic;
        mxnh = ida_mem.ida_maxnh;
    } else {
        ida_mem.ida_cj = ZERO;
        mxnh = 1;
    }

    let mut retval = 0;

    /* Loop over nwt = number of evaluations of ewt vector. */
    'nwt: for _nwt in 1..=2 {
        /* Loop over nh = number of h values. */
        for nh in 1..=mxnh {
            /* Call the IC nonlinear solver function. */
            retval = IDANlsIC(ida_mem);

            /* Cut h and loop on recoverable IDA_YA_YDP_INIT failure; else break. */
            if retval == IDA_SUCCESS {
                break;
            }
            ida_mem.ida_ncfn += 1;
            if retval < 0 {
                break;
            }
            if nh == mxnh {
                break;
            }
            /* If looping to try again, reset yy0 and yp0 if not converging. */
            if retval != IC_SLOW_CONVRG {
                N_VScale(ONE, &ida_mem.ida_phi[0], &mut ida_mem.ida_yy0);
                N_VScale(ONE, &ida_mem.ida_phi[1], &mut ida_mem.ida_yp0);
            }
            hic *= PT1;
            ida_mem.ida_cj = ONE / hic;
            ida_mem.ida_hh = hic;
        } /* End of nh loop */

        /* Break on failure; else reset ewt, save yy0, yp0 in phi, and loop. */
        if retval != IDA_SUCCESS {
            break 'nwt;
        }
        let ewtset_ok = ic_efun(ida_mem);
        if ewtset_ok != 0 {
            retval = IDA_BAD_EWT;
            break 'nwt;
        }
        N_VScale(ONE, &ida_mem.ida_yy0, &mut ida_mem.ida_phi[0]);
        N_VScale(ONE, &ida_mem.ida_yp0, &mut ida_mem.ida_phi[1]);
    } /* End of nwt loop */

    /* Free temporary space */
    ida_mem.ida_yy0 = NVector::default();
    ida_mem.ida_yp0 = NVector::default();

    /* Load the optional outputs. */
    if icopt == IDA_YA_YDP_INIT {
        ida_mem.ida_hused = hic;
    }

    /* On any failure, print message and return proper flag. */
    if retval != IDA_SUCCESS {
        return IDAICFailFlag(ida_mem, retval);
    }

    /* Otherwise return success flag. */
    IDA_SUCCESS
}

/* Error-weight evaluation at yy0 (C: ida_efun(yy0, ewt, edata)). */
fn ic_efun(ida_mem: &mut IDAMem) -> i32 {
    let mut ewt = std::mem::take(&mut ida_mem.ida_ewt);
    let r = if ida_mem.ida_user_efun {
        let efun = ida_mem.ida_efun.unwrap();
        let IDAMem { ida_yy0, ida_user_data, .. } = ida_mem;
        efun(ida_yy0, &mut ewt, ida_user_data)
    } else {
        let yy0 = std::mem::take(&mut ida_mem.ida_yy0);
        let r = IDAEwtSet(ida_mem, &yy0, &mut ewt);
        ida_mem.ida_yy0 = yy0;
        r
    };
    ida_mem.ida_ewt = ewt;
    r
}

/*
 * =================================================================
 * PRIVATE FUNCTIONS IMPLEMENTATION
 * =================================================================
 */

/*
 * IDANlsIC solves a nonlinear system for consistent initial
 * conditions. It calls IDANewtonIC to do most of the work.
 */
fn IDANlsIC(ida_mem: &mut IDAMem) -> i32 {
    /* Evaluate RHS. (writes ida_delta from yy0, yp0) */
    let res = ida_mem.ida_res.unwrap();
    let retval = {
        let IDAMem { ida_t0, ida_yy0, ida_yp0, ida_delta, ida_user_data, .. } = &mut *ida_mem;
        res(*ida_t0, ida_yy0, ida_yp0, ida_delta, ida_user_data)
    };
    ida_mem.ida_nre += 1;
    if retval < 0 {
        return IDA_RES_FAIL;
    }
    if retval > 0 {
        return IDA_FIRST_RES_FAIL;
    }

    /* Save the residual. */
    N_VScale(ONE, &ida_mem.ida_delta, &mut ida_mem.ida_savres);

    let mut retval = IC_SLOW_CONVRG;

    /* Loop over nj = number of linear solve Jacobian setups. */
    for _nj in 1..=ida_mem.ida_maxnj {
        /* If there is a setup routine, call it.
           (tv1 = ee, tv2 = tempv2, tv3 = phi[2] in C are the lsetup
           scratch vectors; idaLsSetup pulls its own temporaries.) */
        if ida_has_lsetup(ida_mem) {
            ida_mem.ida_nsetups += 1;
            let mut lmem = std::mem::take(&mut ida_mem.ida_lmem);
            let r = match &mut lmem {
                LsModule::Ls(idals_mem) => {
                    let yy0 = std::mem::take(&mut ida_mem.ida_yy0);
                    let yp0 = std::mem::take(&mut ida_mem.ida_yp0);
                    let delta = std::mem::take(&mut ida_mem.ida_delta);
                    let r = idaLsSetup(ida_mem, idals_mem, &yy0, &yp0, &delta);
                    ida_mem.ida_yy0 = yy0;
                    ida_mem.ida_yp0 = yp0;
                    ida_mem.ida_delta = delta;
                    r
                }
                LsModule::None => 0,
            };
            ida_mem.ida_lmem = lmem;
            if r < 0 {
                return IDA_LSETUP_FAIL;
            }
            if r > 0 {
                return IC_FAIL_RECOV;
            }
        }

        /* Call the Newton iteration routine, and return if successful. */
        retval = IDANewtonIC(ida_mem);
        if retval == IDA_SUCCESS {
            return IDA_SUCCESS;
        }

        /* If converging slowly and lsetup is nontrivial, retry. */
        if retval == IC_SLOW_CONVRG && ida_has_lsetup(ida_mem) {
            N_VScale(ONE, &ida_mem.ida_savres, &mut ida_mem.ida_delta);
            continue;
        } else {
            return retval;
        }
    } /* End of nj loop */

    /* No convergence after maxnj tries; return with retval = IC_SLOW_CONVRG */
    retval
}

/*
 * IDANewtonIC performs the Newton iteration to solve for consistent
 * initial conditions. It calls IDALineSrch within each iteration.
 * On return, savres contains the current residual vector.
 */
fn IDANewtonIC(ida_mem: &mut IDAMem) -> i32 {
    /* (delnew = phi[2].) Call the linear solve function to get the
       Newton step, delta. */
    let mut lmem = std::mem::take(&mut ida_mem.ida_lmem);
    let retval = match &mut lmem {
        LsModule::Ls(idals_mem) => {
            let mut delta = std::mem::take(&mut ida_mem.ida_delta);
            let ewt = std::mem::take(&mut ida_mem.ida_ewt);
            let yy0 = std::mem::take(&mut ida_mem.ida_yy0);
            let yp0 = std::mem::take(&mut ida_mem.ida_yp0);
            let savres = std::mem::take(&mut ida_mem.ida_savres);
            let r = idaLsSolve(ida_mem, idals_mem, &mut delta, &ewt, &yy0, &yp0, &savres);
            ida_mem.ida_delta = delta;
            ida_mem.ida_ewt = ewt;
            ida_mem.ida_yy0 = yy0;
            ida_mem.ida_yp0 = yp0;
            ida_mem.ida_savres = savres;
            r
        }
        LsModule::None => -1,
    };
    ida_mem.ida_lmem = lmem;
    if retval < 0 {
        return IDA_LSOLVE_FAIL;
    }
    if retval > 0 {
        return IC_FAIL_RECOV;
    }

    /* Compute the norm of the step; return now if this is small. */
    let mut fnorm = IDAWrmsNorm(ida_mem, &ida_mem.ida_delta, &ida_mem.ida_ewt, SUNFALSE);
    if ida_mem.ida_sysindex == 0 {
        fnorm *= ida_mem.ida_tscale * SUNRabs(ida_mem.ida_cj);
    }
    if fnorm <= ida_mem.ida_epsNewt {
        return IDA_SUCCESS;
    }
    let fnorm0 = fnorm;

    /* Initialize rate to avoid compiler warning message */
    let mut rate = ZERO;

    /* Newton iteration loop */
    for _mnewt in 0..ida_mem.ida_maxnit {
        ida_mem.ida_nni += 1;
        let mut delnorm = fnorm;
        let oldfnrm = fnorm;

        /* Call the Linesearch function and return if it failed. */
        let retval = IDALineSrch(ida_mem, &mut delnorm, &mut fnorm);
        if retval != IDA_SUCCESS {
            return retval;
        }

        /* Set the observed convergence rate and test for convergence. */
        rate = fnorm / oldfnrm;
        if fnorm <= ida_mem.ida_epsNewt {
            return IDA_SUCCESS;
        }

        /* If not converged, copy new step vector, and loop. (delnew = phi[2]) */
        N_VScale(ONE, &ida_mem.ida_phi[2], &mut ida_mem.ida_delta);
    } /* End of Newton iteration loop */

    /* Return either IC_SLOW_CONVRG or recoverable fail flag. */
    if rate <= ICRATEMAX || fnorm < PT1 * fnorm0 {
        return IC_SLOW_CONVRG;
    }
    IC_CONV_FAIL
}

/*
 * IDALineSrch performs the Linesearch algorithm with the calculation
 * of consistent initial conditions.
 *
 * On a successful return, yy0, yp0, and savres have been updated,
 * delnew contains the current value of J-inverse F, and fnorm is
 * WRMS-norm(delnew).
 */
fn IDALineSrch(ida_mem: &mut IDAMem, delnorm: &mut f64, fnorm: &mut f64) -> i32 {
    /* Work space aliases: mc = ee, dtemp = phi[3], ynew = tempv2,
       ypnew = ee (use of mc does not conflict with ypnew). */
    let f1norm = (*fnorm) * (*fnorm) * HALF;
    let mut ratio = ONE;

    /* If there are constraints, check and reduce step if necessary. */
    if ida_mem.ida_constraintsSet {
        /* Update y and check constraints. */
        IDANewy(ida_mem);
        let con_ok = N_VConstrMask(&ida_mem.ida_constraints, &ida_mem.ida_tempv2,
                                   &mut ida_mem.ida_ee);

        if !con_ok {
            /* Not satisfied. Compute scaled step to satisfy constraints. */
            N_VProd(&ida_mem.ida_ee, &ida_mem.ida_delta, &mut ida_mem.ida_phi[3]);
            ratio = PT99 * N_VMinQuotient(&ida_mem.ida_yy0, &ida_mem.ida_phi[3]);
            *delnorm *= ratio;
            if *delnorm <= ida_mem.ida_steptol {
                return IC_CONSTR_FAILED;
            }
            ida_mem.ida_delta.scale_inplace(ratio);
        }
    } /* End of constraints check */

    let slpi = -TWO * f1norm * ratio;
    let minlam = ida_mem.ida_steptol / (*delnorm);
    let mut lambda = ONE;
    let mut nbacks = 0;

    /* In IDA_Y_INIT case, set ypnew = yp0 (fixed) for linesearch. */
    if ida_mem.ida_icopt == IDA_Y_INIT {
        N_VScale(ONE, &ida_mem.ida_yp0, &mut ida_mem.ida_ee);
    }

    /* Loop on linesearch variable lambda. */
    let mut fnormp = ZERO;
    loop {
        if nbacks == ida_mem.ida_maxbacks {
            return IC_LINESRCH_FAILED;
        }
        /* Get new (y,y') = (ynew,ypnew) and norm of new function value. */
        IDANewyyp(ida_mem, lambda);
        let retval = IDAfnorm(ida_mem, &mut fnormp);
        if retval != IDA_SUCCESS {
            return retval;
        }

        /* If lsoff option is on, break out. */
        if ida_mem.ida_lsoff {
            break;
        }

        /* Do alpha-condition test. */
        let f1normp = fnormp * fnormp * HALF;
        if f1normp <= f1norm + ALPHALS * slpi * lambda {
            break;
        }
        if lambda < minlam {
            return IC_LINESRCH_FAILED;
        }
        lambda /= TWO;
        ida_mem.ida_nbacktr += 1;
        nbacks += 1;
    } /* End of breakout linesearch loop */

    /* Update yy0, yp0, and fnorm, then return. (ynew = tempv2, ypnew = ee) */
    N_VScale(ONE, &ida_mem.ida_tempv2, &mut ida_mem.ida_yy0);
    if ida_mem.ida_icopt == IDA_YA_YDP_INIT {
        N_VScale(ONE, &ida_mem.ida_ee, &mut ida_mem.ida_yp0);
    }
    *fnorm = fnormp;
    IDA_SUCCESS
}

/*
 * IDAfnorm computes the norm of the current function value, by
 * evaluating the DAE residual function, calling the linear system
 * solver, and computing a WRMS-norm.
 *
 * On return, savres contains the current residual vector F, and
 * delnew contains J-inverse F.
 */
fn IDAfnorm(ida_mem: &mut IDAMem, fnorm: &mut f64) -> i32 {
    /* Get residual vector F, return if failed, and save F in savres.
       (ynew = tempv2, ypnew = ee, delnew = phi[2]) */
    let res = ida_mem.ida_res.unwrap();
    let retval = {
        let IDAMem { ida_t0, ida_tempv2, ida_ee, ida_phi, ida_user_data, .. } = &mut *ida_mem;
        res(*ida_t0, ida_tempv2, ida_ee, &mut ida_phi[2], ida_user_data)
    };
    ida_mem.ida_nre += 1;
    if retval < 0 {
        return IDA_RES_FAIL;
    }
    if retval > 0 {
        return IC_FAIL_RECOV;
    }

    N_VScale(ONE, &ida_mem.ida_phi[2], &mut ida_mem.ida_savres);

    /* Call the linear solve function to get J-inverse F; return if failed. */
    let mut lmem = std::mem::take(&mut ida_mem.ida_lmem);
    let retval = match &mut lmem {
        LsModule::Ls(idals_mem) => {
            let mut delnew = std::mem::take(&mut ida_mem.ida_phi[2]);
            let ewt = std::mem::take(&mut ida_mem.ida_ewt);
            let ynew = std::mem::take(&mut ida_mem.ida_tempv2);
            let ypnew = std::mem::take(&mut ida_mem.ida_ee);
            let savres = std::mem::take(&mut ida_mem.ida_savres);
            let r = idaLsSolve(ida_mem, idals_mem, &mut delnew, &ewt, &ynew, &ypnew, &savres);
            ida_mem.ida_phi[2] = delnew;
            ida_mem.ida_ewt = ewt;
            ida_mem.ida_tempv2 = ynew;
            ida_mem.ida_ee = ypnew;
            ida_mem.ida_savres = savres;
            r
        }
        LsModule::None => -1,
    };
    ida_mem.ida_lmem = lmem;
    if retval < 0 {
        return IDA_LSOLVE_FAIL;
    }
    if retval > 0 {
        return IC_FAIL_RECOV;
    }

    /* Compute the WRMS-norm; rescale if index = 0. */
    *fnorm = IDAWrmsNorm(ida_mem, &ida_mem.ida_phi[2], &ida_mem.ida_ewt, SUNFALSE);
    if ida_mem.ida_sysindex == 0 {
        *fnorm *= ida_mem.ida_tscale * SUNRabs(ida_mem.ida_cj);
    }

    IDA_SUCCESS
}

/*
 * IDANewyyp updates the vectors ynew and ypnew from yy0 and yp0,
 * using the current step vector lambda*delta, in a manner depending
 * on icopt and the input id vector. (ynew = tempv2, ypnew = ee,
 * dtemp = phi[3].)
 */
fn IDANewyyp(ida_mem: &mut IDAMem, lambda: f64) -> i32 {
    /* IDA_YA_YDP_INIT case: ynew  = yy0 - lambda*delta    where id_i = 0
                             ypnew = yp0 - cj*lambda*delta where id_i = 1. */
    if ida_mem.ida_icopt == IDA_YA_YDP_INIT {
        N_VProd(&ida_mem.ida_id, &ida_mem.ida_delta, &mut ida_mem.ida_phi[3]);
        {
            let neg = -ida_mem.ida_cj * lambda;
            N_VLinearSum(ONE, &ida_mem.ida_yp0, neg, &ida_mem.ida_phi[3], &mut ida_mem.ida_ee);
        }
        /* dtemp = delta - dtemp (C: N_VLinearSum(1, delta, -1, dtemp, dtemp)) */
        {
            let IDAMem { ida_phi, ida_delta, .. } = &mut *ida_mem;
            ida_phi[3].linear_sum_with(-ONE, ONE, ida_delta);
        }
        N_VLinearSum(ONE, &ida_mem.ida_yy0, -lambda, &ida_mem.ida_phi[3], &mut ida_mem.ida_tempv2);
        return IDA_SUCCESS;
    }

    /* IDA_Y_INIT case: ynew = yy0 - lambda*delta. (ypnew = yp0 preset.) */
    N_VLinearSum(ONE, &ida_mem.ida_yy0, -lambda, &ida_mem.ida_delta, &mut ida_mem.ida_tempv2);
    IDA_SUCCESS
}

/*
 * IDANewy updates the vector ynew from yy0, using the current step
 * vector delta, in a manner depending on icopt and the input id
 * vector. (ynew = tempv2, dtemp = phi[3].)
 */
fn IDANewy(ida_mem: &mut IDAMem) -> i32 {
    /* IDA_YA_YDP_INIT case: ynew = yy0 - delta    where id_i = 0. */
    if ida_mem.ida_icopt == IDA_YA_YDP_INIT {
        N_VProd(&ida_mem.ida_id, &ida_mem.ida_delta, &mut ida_mem.ida_phi[3]);
        /* dtemp = delta - dtemp (C: N_VLinearSum(1, delta, -1, dtemp, dtemp)) */
        {
            let IDAMem { ida_phi, ida_delta, .. } = &mut *ida_mem;
            ida_phi[3].linear_sum_with(-ONE, ONE, ida_delta);
        }
        N_VLinearSum(ONE, &ida_mem.ida_yy0, -ONE, &ida_mem.ida_phi[3], &mut ida_mem.ida_tempv2);
        return IDA_SUCCESS;
    }

    /* IDA_Y_INIT case: ynew = yy0 - delta. */
    N_VLinearSum(ONE, &ida_mem.ida_yy0, -ONE, &ida_mem.ida_delta, &mut ida_mem.ida_tempv2);
    IDA_SUCCESS
}

/*
 * IDAICFailFlag prints a message and sets the IDACalcIC return value
 * appropriate to the flag retval returned by IDANlsIC.
 */
fn IDAICFailFlag(ida_mem: &mut IDAMem, retval: i32) -> i32 {
    match retval {
        IDA_RES_FAIL => {
            IDAProcessError(Some(ida_mem), IDA_RES_FAIL, line!(), "IDACalcIC", file!(),
                            MSG_IC_RES_NONREC);
            IDA_RES_FAIL
        }
        IDA_FIRST_RES_FAIL => {
            IDAProcessError(Some(ida_mem), IDA_FIRST_RES_FAIL, line!(), "IDACalcIC", file!(),
                            MSG_IC_RES_FAIL);
            IDA_FIRST_RES_FAIL
        }
        IDA_LSETUP_FAIL => {
            IDAProcessError(Some(ida_mem), IDA_LSETUP_FAIL, line!(), "IDACalcIC", file!(),
                            MSG_IC_SETUP_FAIL);
            IDA_LSETUP_FAIL
        }
        IDA_LSOLVE_FAIL => {
            IDAProcessError(Some(ida_mem), IDA_LSOLVE_FAIL, line!(), "IDACalcIC", file!(),
                            MSG_IC_SOLVE_FAIL);
            IDA_LSOLVE_FAIL
        }
        IC_FAIL_RECOV => {
            IDAProcessError(Some(ida_mem), IDA_NO_RECOVERY, line!(), "IDACalcIC", file!(),
                            MSG_IC_NO_RECOVERY);
            IDA_NO_RECOVERY
        }
        IC_CONSTR_FAILED => {
            IDAProcessError(Some(ida_mem), IDA_CONSTR_FAIL, line!(), "IDACalcIC", file!(),
                            MSG_IC_FAIL_CONSTR);
            IDA_CONSTR_FAIL
        }
        IC_LINESRCH_FAILED => {
            IDAProcessError(Some(ida_mem), IDA_LINESEARCH_FAIL, line!(), "IDACalcIC", file!(),
                            MSG_IC_FAILED_LINS);
            IDA_LINESEARCH_FAIL
        }
        IC_CONV_FAIL => {
            IDAProcessError(Some(ida_mem), IDA_CONV_FAIL, line!(), "IDACalcIC", file!(),
                            MSG_IC_CONV_FAILED);
            IDA_CONV_FAIL
        }
        IC_SLOW_CONVRG => {
            IDAProcessError(Some(ida_mem), IDA_CONV_FAIL, line!(), "IDACalcIC", file!(),
                            MSG_IC_CONV_FAILED);
            IDA_CONV_FAIL
        }
        IDA_BAD_EWT => {
            IDAProcessError(Some(ida_mem), IDA_BAD_EWT, line!(), "IDACalcIC", file!(),
                            MSG_IC_BAD_EWT);
            IDA_BAD_EWT
        }
        _ => -99,
    }
}
