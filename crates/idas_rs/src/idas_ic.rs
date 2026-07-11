/* -----------------------------------------------------------------
 * Translated from src/idas/idas_ic.c (IDAS 7.7.0).
 * IC (consistent initial condition) calculation for IDAS; independent
 * of the linear solver in use.
 *
 * The C code aliases IDACalcIC work "vectors" onto existing IDAMem
 * fields: ynew = tempv2, ypnew = ee, delnew = phi[2], dtemp = phi[3],
 * mc = ee (used before ypnew), and for sensitivities savresS =
 * phiS[2], delnewS = phiS[3], yyS0new = phiS[4], ypS0new = eeS,
 * tmpS1 = tempv1, tmpS2 = tempv2.  This port references the
 * underlying fields directly (aliases not stored, pinned convention).
 * Vectors are detached with std::mem::take only where the
 * linear-solver or DQ-residual dispatch needs &mut IDAMem alongside
 * borrowed vector arguments.
 *
 * NOT registered in lib.rs yet; registers together with idas.rs once
 * idas_nls / idas_nls_sim / idas_nls_stg and idas_ls land.
 * -----------------------------------------------------------------*/
use crate::idas::{IDAEwtSet, IDAInitialSetup, IDASensEwtSet, IDASensResDQ, IDASensWrmsNorm,
                  IDASensWrmsNormUpdate, IDAWrmsNorm};
use crate::idas_impl::*;
use crate::idas_ls::{idaLsSetup, idaLsSolve};
use crate::idas_nls::ida_has_lsetup;
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
 * -----------------------------------------------------------------
 * IDACalcIC
 * -----------------------------------------------------------------
 * IDACalcIC computes consistent initial conditions, given the
 * user's initial guess for unknown components of yy0 and/or yp0.
 *
 * The return value is IDA_SUCCESS = 0 if no error occurred.
 *
 * The error return values (fully described in ida.h) are:
 *   IDA_MEM_NULL        ida_mem is NULL
 *   IDA_NO_MALLOC       ida_mem was not allocated
 *   IDA_ILL_INPUT       bad value for icopt, tout1, or id
 *   IDA_LINIT_FAIL      the linear solver linit routine failed
 *   IDA_BAD_EWT         zero value of some component of ewt
 *   IDA_RES_FAIL        res had a non-recoverable error
 *   IDA_FIRST_RES_FAIL  res failed recoverably on the first call
 *   IDA_LSETUP_FAIL     lsetup had a non-recoverable error
 *   IDA_LSOLVE_FAIL     lsolve had a non-recoverable error
 *   IDA_NO_RECOVERY     res, lsetup, or lsolve had a recoverable
 *                       error, but IDACalcIC could not recover
 *   IDA_CONSTR_FAIL     the inequality constraints could not be met
 *   IDA_LINESEARCH_FAIL the linesearch failed (either on steptol test
 *                       or on the maxbacks test)
 *   IDA_CONV_FAIL       the Newton iterations failed to converge
 * -----------------------------------------------------------------
 */
pub fn IDACalcIC(ida_mem: &mut IDAMem, icopt: i32, tout1: f64) -> i32 {
    /* (C's ida_mem NULL check vanishes: &mut IDAMem) */

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

    /* Are we computing sensitivities? */
    let sensi_stg = ida_mem.ida_sensi && ida_mem.ida_ism == IDA_STAGGERED;
    let sensi_sim = ida_mem.ida_sensi && ida_mem.ida_ism == IDA_SIMULTANEOUS;

    /* Allocate space and initialize temporary vectors */
    ida_mem.ida_yy0 = NVector::new(ida_mem.ida_ee.len());
    ida_mem.ida_yp0 = NVector::new(ida_mem.ida_ee.len());
    ida_mem.ida_t0 = ida_mem.ida_tn;
    N_VScale(ONE, &ida_mem.ida_phi[0], &mut ida_mem.ida_yy0);
    N_VScale(ONE, &ida_mem.ida_phi[1], &mut ida_mem.ida_yp0);

    if ida_mem.ida_sensi {
        /* Allocate temporary space required for sensitivity IC: yyS0 and ypS0. */
        let n = ida_mem.ida_ee.len();
        let ns = ida_mem.ida_Ns as usize;
        ida_mem.ida_yyS0 = (0..ns).map(|_| NVector::new(n)).collect();
        ida_mem.ida_ypS0 = (0..ns).map(|_| NVector::new(n)).collect();

        /* Initialize sensitivity vector. */
        for is in 0..ns {
            N_VScale(ONE, &ida_mem.ida_phiS[0][is], &mut ida_mem.ida_yyS0[is]);
            N_VScale(ONE, &ida_mem.ida_phiS[1][is], &mut ida_mem.ida_ypS0[is]);
        }

        /* Initialize work space vectors needed for sensitivities:
           savresS = phiS[2], delnewS = phiS[3], yyS0new = phiS[4],
           ypS0new = eeS.  (Pointer aliases in C; the underlying
           fields are referenced directly at use sites here.) */
    }

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

    /* Initializations:
       cjratio = 1 (for use in direct linear solvers);
       set nbacktr = 0; */
    ida_mem.ida_cjratio = ONE;
    ida_mem.ida_nbacktr = 0;

    /* Set hic, hh, cj, and mxnh. */
    let mut hic = PT001 * tdist;
    let mut ypnorm = IDAWrmsNorm(ida_mem, &ida_mem.ida_yp0, &ida_mem.ida_ewt,
                                 ida_mem.ida_suppressalg);

    if sensi_sim {
        ypnorm = IDASensWrmsNormUpdate(ida_mem, ypnorm, &ida_mem.ida_ypS0,
                                       &ida_mem.ida_ewtS, SUNFALSE);
    }

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
                if sensi_sim {
                    /* Reset yyS0 and ypS0. */
                    /* Copy phiS[0] and phiS[1] into yyS0 and ypS0. */
                    for is in 0..(ida_mem.ida_Ns as usize) {
                        N_VScale(ONE, &ida_mem.ida_phiS[0][is], &mut ida_mem.ida_yyS0[is]);
                        N_VScale(ONE, &ida_mem.ida_phiS[1][is], &mut ida_mem.ida_ypS0[is]);
                    }
                }
            }
            hic *= PT1;
            ida_mem.ida_cj = ONE / hic;
            ida_mem.ida_hh = hic;
        } /* End of nh loop */

        /* Break on failure */
        if retval != IDA_SUCCESS {
            break 'nwt;
        }

        /* Reset ewt, save yy0, yp0 in phi, and loop. */
        let ewtset_ok = ic_efun(ida_mem);
        if ewtset_ok != 0 {
            retval = IDA_BAD_EWT;
            break 'nwt;
        }
        N_VScale(ONE, &ida_mem.ida_yy0, &mut ida_mem.ida_phi[0]);
        N_VScale(ONE, &ida_mem.ida_yp0, &mut ida_mem.ida_phi[1]);

        if sensi_sim {
            /* Reevaluate ewtS. */
            let yyS0 = std::mem::take(&mut ida_mem.ida_yyS0);
            let mut ewtS = std::mem::take(&mut ida_mem.ida_ewtS);
            let ewtset_ok = IDASensEwtSet(ida_mem, &yyS0, &mut ewtS);
            ida_mem.ida_yyS0 = yyS0;
            ida_mem.ida_ewtS = ewtS;
            if ewtset_ok != 0 {
                retval = IDA_BAD_EWT;
                break 'nwt;
            }

            /* Save yyS0 and ypS0. */
            for is in 0..(ida_mem.ida_Ns as usize) {
                N_VScale(ONE, &ida_mem.ida_yyS0[is], &mut ida_mem.ida_phiS[0][is]);
                N_VScale(ONE, &ida_mem.ida_ypS0[is], &mut ida_mem.ida_phiS[1][is]);
            }
        }
    } /* End of nwt loop */

    /* Load the optional outputs. */
    if icopt == IDA_YA_YDP_INIT {
        ida_mem.ida_hused = hic;
    }

    /* On any failure, free memory, print error message and return */
    if retval != IDA_SUCCESS {
        ida_mem.ida_yy0 = NVector::default();
        ida_mem.ida_yp0 = NVector::default();

        if ida_mem.ida_sensi {
            ida_mem.ida_yyS0 = Vec::new();
            ida_mem.ida_ypS0 = Vec::new();
        }

        return IDAICFailFlag(ida_mem, retval);
    }

    /* Unless using the STAGGERED approach for sensitivities, return now */
    if !sensi_stg {
        ida_mem.ida_yy0 = NVector::default();
        ida_mem.ida_yp0 = NVector::default();

        if ida_mem.ida_sensi {
            ida_mem.ida_yyS0 = Vec::new();
            ida_mem.ida_ypS0 = Vec::new();
        }

        return IDA_SUCCESS;
    }

    /* Find consistent I.C. for sensitivities using a staggered approach */

    /* Evaluate res at converged y, needed for future evaluations of sens. RHS
       If res() fails recoverably, treat it as a convergence failure and
       attempt the step again */

    let res = ida_mem.ida_res.unwrap();
    retval = {
        let IDAMem { ida_t0, ida_yy0, ida_yp0, ida_delta, ida_user_data, .. } = &mut *ida_mem;
        res(*ida_t0, ida_yy0, ida_yp0, ida_delta, ida_user_data)
    };
    ida_mem.ida_nre += 1;
    if retval < 0 {
        /* res function failed unrecoverably.
           (C returns without freeing yy0/yp0/yyS0/ypS0; the owned
           fields simply retain the clones — harmless.) */
        return IDA_RES_FAIL;
    }
    if retval > 0 {
        /* res function failed recoverably but no recovery possible. */
        return IDA_FIRST_RES_FAIL;
    }

    /* Loop over nwt = number of evaluations of ewt vector. */
    'nwt2: for _nwt in 1..=2 {
        /* Loop over nh = number of h values. */
        for nh in 1..=mxnh {
            retval = IDASensNlsIC(ida_mem);
            if retval == IDA_SUCCESS {
                break;
            }

            /* Increment the number of the sensitivity related corrector convergence failures. */
            ida_mem.ida_ncfnS += 1;

            if retval < 0 {
                break;
            }
            if nh == mxnh {
                break;
            }

            /* If looping to try again, reset yyS0 and ypS0 if not converging. */
            if retval != IC_SLOW_CONVRG {
                for is in 0..(ida_mem.ida_Ns as usize) {
                    N_VScale(ONE, &ida_mem.ida_phiS[0][is], &mut ida_mem.ida_yyS0[is]);
                    N_VScale(ONE, &ida_mem.ida_phiS[1][is], &mut ida_mem.ida_ypS0[is]);
                }
            }
            hic *= PT1;
            ida_mem.ida_cj = ONE / hic;
            ida_mem.ida_hh = hic;
        } /* End of nh loop */

        /* Break on failure */
        if retval != IDA_SUCCESS {
            break 'nwt2;
        }

        /* Since it was successful, reevaluate ewtS with the new values of yyS0, save
           yyS0 and ypS0 in phiS[0] and phiS[1] and loop one more time to check and
           maybe correct the  new sensitivities IC with respect to the new weights. */

        /* Reevaluate ewtS. */
        let yyS0 = std::mem::take(&mut ida_mem.ida_yyS0);
        let mut ewtS = std::mem::take(&mut ida_mem.ida_ewtS);
        let ewtset_ok = IDASensEwtSet(ida_mem, &yyS0, &mut ewtS);
        ida_mem.ida_yyS0 = yyS0;
        ida_mem.ida_ewtS = ewtS;
        if ewtset_ok != 0 {
            retval = IDA_BAD_EWT;
            break 'nwt2;
        }

        /* Save yyS0 and ypS0. */
        for is in 0..(ida_mem.ida_Ns as usize) {
            N_VScale(ONE, &ida_mem.ida_yyS0[is], &mut ida_mem.ida_phiS[0][is]);
            N_VScale(ONE, &ida_mem.ida_ypS0[is], &mut ida_mem.ida_phiS[1][is]);
        }
    } /* End of nwt loop */

    /* Load the optional outputs. */
    if icopt == IDA_YA_YDP_INIT {
        ida_mem.ida_hused = hic;
    }

    /* Free temporary space */
    ida_mem.ida_yy0 = NVector::default();
    ida_mem.ida_yp0 = NVector::default();

    /* Here sensi is SUNTRUE, so deallocate sensitivity temporary vectors. */
    ida_mem.ida_yyS0 = Vec::new();
    ida_mem.ida_ypS0 = Vec::new();

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
 * -----------------------------------------------------------------
 * IDANlsIC
 * -----------------------------------------------------------------
 * IDANlsIC solves a nonlinear system for consistent initial
 * conditions.  It calls IDANewtonIC to do most of the work.
 *
 * The return value is IDA_SUCCESS = 0 if no error occurred.
 * The error return values (positive) considered recoverable are:
 *  IC_FAIL_RECOV      if res, lsetup, or lsolve failed recoverably
 *  IC_CONSTR_FAILED   if the constraints could not be met
 *  IC_LINESRCH_FAILED if the linesearch failed (either on steptol test
 *                     or on maxbacks test)
 *  IC_CONV_FAIL       if the Newton iterations failed to converge
 *  IC_SLOW_CONVRG     if the iterations are converging slowly
 *                     (failed the convergence test, but showed
 *                     norm reduction or convergence rate < 1)
 * The error return values (negative) considered non-recoverable are:
 *  IDA_RES_FAIL       if res had a non-recoverable error
 *  IDA_FIRST_RES_FAIL if res failed recoverably on the first call
 *  IDA_LSETUP_FAIL    if lsetup had a non-recoverable error
 *  IDA_LSOLVE_FAIL    if lsolve had a non-recoverable error
 * -----------------------------------------------------------------
 */
fn IDANlsIC(ida_mem: &mut IDAMem) -> i32 {
    /* Are we computing sensitivities with the IDA_SIMULTANEOUS approach? */
    let sensi_sim = ida_mem.ida_sensi && ida_mem.ida_ism == IDA_SIMULTANEOUS;

    /* (tv1 = ee, tv2 = tempv2, tv3 = phi[2] are C's lsetup scratch
       vectors; idaLsSetup pulls its own temporaries.) */

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

    if sensi_sim {
        /* Evaluate sensitivity RHS and save it in savresS (= phiS[2]). */
        let ns = ida_mem.ida_Ns;
        let t0 = ida_mem.ida_t0;
        let yy0 = std::mem::take(&mut ida_mem.ida_yy0);
        let yp0 = std::mem::take(&mut ida_mem.ida_yp0);
        let delta = std::mem::take(&mut ida_mem.ida_delta);
        let yyS0 = std::mem::take(&mut ida_mem.ida_yyS0);
        let ypS0 = std::mem::take(&mut ida_mem.ida_ypS0);
        let mut deltaS = std::mem::take(&mut ida_mem.ida_deltaS);
        /* tmpS1 = tempv1, tmpS2 = tempv2 (aliases), tmpS3 real */
        let mut tmp1 = std::mem::take(&mut ida_mem.ida_tempv1);
        let mut tmp2 = std::mem::take(&mut ida_mem.ida_tempv2);
        let mut tmp3 = std::mem::take(&mut ida_mem.ida_tmpS3);

        let retval = if ida_mem.ida_resSDQ {
            IDASensResDQ(ida_mem, ns, t0, &yy0, &yp0, &delta, &yyS0, &ypS0, &mut deltaS,
                         &mut tmp1, &mut tmp2, &mut tmp3)
        } else {
            let resS = ida_mem.ida_resS.unwrap();
            resS(ns, t0, &yy0, &yp0, &delta, &yyS0, &ypS0, &mut deltaS,
                 &mut ida_mem.ida_user_data, &mut tmp1, &mut tmp2, &mut tmp3)
        };

        ida_mem.ida_yy0 = yy0;
        ida_mem.ida_yp0 = yp0;
        ida_mem.ida_delta = delta;
        ida_mem.ida_yyS0 = yyS0;
        ida_mem.ida_ypS0 = ypS0;
        ida_mem.ida_deltaS = deltaS;
        ida_mem.ida_tempv1 = tmp1;
        ida_mem.ida_tempv2 = tmp2;
        ida_mem.ida_tmpS3 = tmp3;

        ida_mem.ida_nrSe += 1;
        if retval < 0 {
            return IDA_RES_FAIL;
        }
        if retval > 0 {
            return IDA_FIRST_RES_FAIL;
        }

        for is in 0..(ida_mem.ida_Ns as usize) {
            N_VScale(ONE, &ida_mem.ida_deltaS[is], &mut ida_mem.ida_phiS[2][is]);
        }
    }

    let mut retval = IC_SLOW_CONVRG;

    /* Loop over nj = number of linear solve Jacobian setups. */
    for _nj in 1..=ida_mem.ida_maxnj {
        /* If there is a setup routine, call it. */
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

            if sensi_sim {
                for is in 0..(ida_mem.ida_Ns as usize) {
                    N_VScale(ONE, &ida_mem.ida_phiS[2][is], &mut ida_mem.ida_deltaS[is]);
                }
            }

            continue;
        } else {
            return retval;
        }
    } /* End of nj loop */

    /* No convergence after maxnj tries; return with retval=IC_SLOW_CONVRG */
    retval
}

/*
 * -----------------------------------------------------------------
 * IDANewtonIC
 * -----------------------------------------------------------------
 * IDANewtonIC performs the Newton iteration to solve for consistent
 * initial conditions.  It calls IDALineSrch within each iteration.
 * On return, savres contains the current residual vector.
 *
 * The return value is IDA_SUCCESS = 0 if no error occurred.
 * The error return values (positive) considered recoverable are:
 *  IC_FAIL_RECOV      if res or lsolve failed recoverably
 *  IC_CONSTR_FAILED   if the constraints could not be met
 *  IC_LINESRCH_FAILED if the linesearch failed (either on steptol test
 *                     or on maxbacks test)
 *  IC_CONV_FAIL       if the Newton iterations failed to converge
 *  IC_SLOW_CONVRG     if the iterations appear to be converging slowly.
 *                     They failed the convergence test, but showed
 *                     an overall norm reduction (by a factor of < 0.1)
 *                     or a convergence rate <= ICRATEMAX).
 * The error return values (negative) considered non-recoverable are:
 *  IDA_RES_FAIL   if res had a non-recoverable error
 *  IDA_LSOLVE_FAIL      if lsolve had a non-recoverable error
 * -----------------------------------------------------------------
 */
fn IDANewtonIC(ida_mem: &mut IDAMem) -> i32 {
    /* Are we computing sensitivities with the IDA_SIMULTANEOUS approach? */
    let sensi_sim = ida_mem.ida_sensi && ida_mem.ida_ism == IDA_SIMULTANEOUS;

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

    /* Compute the norm of the step. */
    let mut fnorm = IDAWrmsNorm(ida_mem, &ida_mem.ida_delta, &ida_mem.ida_ewt, SUNFALSE);

    /* Call the lsolve function to get correction vectors deltaS. */
    if sensi_sim {
        let ns = ida_mem.ida_Ns as usize;
        let mut lmem = std::mem::take(&mut ida_mem.ida_lmem);
        let retval = match &mut lmem {
            LsModule::Ls(idals_mem) => {
                let mut deltaS = std::mem::take(&mut ida_mem.ida_deltaS);
                let ewtS = std::mem::take(&mut ida_mem.ida_ewtS);
                let yy0 = std::mem::take(&mut ida_mem.ida_yy0);
                let yp0 = std::mem::take(&mut ida_mem.ida_yp0);
                let savres = std::mem::take(&mut ida_mem.ida_savres);
                let mut r = 0;
                for is in 0..ns {
                    r = idaLsSolve(ida_mem, idals_mem, &mut deltaS[is], &ewtS[is], &yy0,
                                   &yp0, &savres);
                    if r != 0 {
                        break;
                    }
                }
                ida_mem.ida_deltaS = deltaS;
                ida_mem.ida_ewtS = ewtS;
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

        /* Update the norm of delta. */
        fnorm = IDASensWrmsNormUpdate(ida_mem, fnorm, &ida_mem.ida_deltaS,
                                      &ida_mem.ida_ewtS, SUNFALSE);
    }

    /* Test for convergence. Return now if the norm is small. */
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

        if sensi_sim {
            /* Update the iteration's step for sensitivities. (delnewS = phiS[3]) */
            for is in 0..(ida_mem.ida_Ns as usize) {
                N_VScale(ONE, &ida_mem.ida_phiS[3][is], &mut ida_mem.ida_deltaS[is]);
            }
        }
    } /* End of Newton iteration loop */

    /* Return either IC_SLOW_CONVRG or recoverable fail flag. */
    if rate <= ICRATEMAX || fnorm < PT1 * fnorm0 {
        return IC_SLOW_CONVRG;
    }
    IC_CONV_FAIL
}

/*
 * -----------------------------------------------------------------
 * IDALineSrch
 * -----------------------------------------------------------------
 * IDALineSrch performs the Linesearch algorithm with the
 * calculation of consistent initial conditions.
 *
 * On entry, yy0 and yp0 are the current values of y and y', the
 * Newton step is delta, the current residual vector F is savres,
 * delnorm is WRMS-norm(delta), and fnorm is the norm of the vector
 * J-inverse F.
 *
 * On a successful return, yy0, yp0, and savres have been updated,
 * delnew contains the current value of J-inverse F, and fnorm is
 * WRMS-norm(delnew).
 *
 * The return value is IDA_SUCCESS = 0 if no error occurred.
 * The error return values (positive) considered recoverable are:
 *  IC_FAIL_RECOV      if res or lsolve failed recoverably
 *  IC_CONSTR_FAILED   if the constraints could not be met
 *  IC_LINESRCH_FAILED if the linesearch failed (either on steptol test
 *                     or on maxbacks test)
 * The error return values (negative) considered non-recoverable are:
 *  IDA_RES_FAIL   if res had a non-recoverable error
 *  IDA_LSOLVE_FAIL      if lsolve had a non-recoverable error
 * -----------------------------------------------------------------
 */
fn IDALineSrch(ida_mem: &mut IDAMem, delnorm: &mut f64, fnorm: &mut f64) -> i32 {
    /* Initialize work space pointers, f1norm, ratio.
       (Use of mc in constraint check does not conflict with ypnew.)
       Work space aliases: mc = ee, dtemp = phi[3], ynew = tempv2,
       ypnew = ee. */
    let f1norm = (*fnorm) * (*fnorm) * HALF;
    let mut ratio = ONE;

    /* If there are constraints, check and reduce step if necessary. */
    if ida_mem.ida_constraintsSet {
        /* Update y and check constraints. */
        IDANewy(ida_mem);
        let con_ok = N_VConstrMask(&ida_mem.ida_constraints, &ida_mem.ida_tempv2,
                                   &mut ida_mem.ida_ee);

        if !con_ok {
            /* Not satisfied.  Compute scaled step to satisfy constraints. */
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

    /* Are we computing sensitivities with the IDA_SIMULTANEOUS approach? */
    let sensi_sim = ida_mem.ida_sensi && ida_mem.ida_ism == IDA_SIMULTANEOUS;

    /* In IDA_Y_INIT case, set ypnew = yp0 (fixed) for linesearch. (ypnew = ee) */
    if ida_mem.ida_icopt == IDA_Y_INIT {
        N_VScale(ONE, &ida_mem.ida_yp0, &mut ida_mem.ida_ee);

        /* do the same for sensitivities. (ypS0new = eeS) */
        if sensi_sim {
            for is in 0..(ida_mem.ida_Ns as usize) {
                N_VScale(ONE, &ida_mem.ida_ypS0[is], &mut ida_mem.ida_eeS[is]);
            }
        }
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

    /* Update yy0, yp0. (ynew = tempv2) */
    N_VScale(ONE, &ida_mem.ida_tempv2, &mut ida_mem.ida_yy0);

    if sensi_sim {
        /* Update yyS0 and ypS0. (yyS0new = phiS[4]) */
        for is in 0..(ida_mem.ida_Ns as usize) {
            N_VScale(ONE, &ida_mem.ida_phiS[4][is], &mut ida_mem.ida_yyS0[is]);
        }
    }

    if ida_mem.ida_icopt == IDA_YA_YDP_INIT {
        N_VScale(ONE, &ida_mem.ida_ee, &mut ida_mem.ida_yp0);

        if sensi_sim {
            /* (ypS0new = eeS) */
            for is in 0..(ida_mem.ida_Ns as usize) {
                N_VScale(ONE, &ida_mem.ida_eeS[is], &mut ida_mem.ida_ypS0[is]);
            }
        }
    }
    /* Update fnorm, then return. */
    *fnorm = fnormp;
    IDA_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * IDAfnorm
 * -----------------------------------------------------------------
 * IDAfnorm computes the norm of the current function value, by
 * evaluating the DAE residual function, calling the linear
 * system solver, and computing a WRMS-norm.
 *
 * On return, savres contains the current residual vector F, and
 * delnew contains J-inverse F.
 *
 * The return value is IDA_SUCCESS = 0 if no error occurred, or
 *  IC_FAIL_RECOV    if res or lsolve failed recoverably, or
 *  IDA_RES_FAIL     if res had a non-recoverable error, or
 *  IDA_LSOLVE_FAIL  if lsolve had a non-recoverable error.
 * -----------------------------------------------------------------
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

    /* Compute the WRMS-norm. */
    *fnorm = IDAWrmsNorm(ida_mem, &ida_mem.ida_phi[2], &ida_mem.ida_ewt, SUNFALSE);

    /* Are we computing SENSITIVITIES with the IDA_SIMULTANEOUS approach? */
    if ida_mem.ida_sensi && ida_mem.ida_ism == IDA_SIMULTANEOUS {
        /* Evaluate the residual for sensitivities.
           (ynew = tempv2, ypnew = ee, yyS0new = phiS[4],
           ypS0new = eeS, delnewS = phiS[3]) */
        let ns = ida_mem.ida_Ns;
        let t0 = ida_mem.ida_t0;
        let ynew = std::mem::take(&mut ida_mem.ida_tempv2);
        let ypnew = std::mem::take(&mut ida_mem.ida_ee);
        let savres = std::mem::take(&mut ida_mem.ida_savres);
        let eeS = std::mem::take(&mut ida_mem.ida_eeS);
        let mut phiS = std::mem::take(&mut ida_mem.ida_phiS);
        let mut tmp1 = std::mem::take(&mut ida_mem.ida_tempv1);
        let mut tmp3 = std::mem::take(&mut ida_mem.ida_tmpS3);
        /* UPSTREAM ALIASING DEFECT (documented deviation): C passes
           tempv2 BOTH as ynew (the yy argument) and as tmpS2 (the
           yptemp scratch).  With the internal DQ residual (CENTERED1,
           the default) the scratch writes clobber ynew mid-evaluation
           in C.  That aliasing is inexpressible under Rust borrow
           rules; this port uses a fresh scratch vector, i.e. the
           intended semantics.  Revisit only if a sim-sensitivity
           IDACalcIC example ever diffs against the C reference. */
        let mut tmp2 = NVector::new(tmp1.len());

        let retval = {
            let (lo, hi) = phiS.split_at_mut(4);
            let yyS0new = &hi[0];
            let delnewS = &mut lo[3];
            if ida_mem.ida_resSDQ {
                IDASensResDQ(ida_mem, ns, t0, &ynew, &ypnew, &savres, yyS0new, &eeS,
                             delnewS, &mut tmp1, &mut tmp2, &mut tmp3)
            } else {
                let resS = ida_mem.ida_resS.unwrap();
                resS(ns, t0, &ynew, &ypnew, &savres, yyS0new, &eeS, delnewS,
                     &mut ida_mem.ida_user_data, &mut tmp1, &mut tmp2, &mut tmp3)
            }
        };

        ida_mem.ida_tempv2 = ynew;
        ida_mem.ida_ee = ypnew;
        ida_mem.ida_savres = savres;
        ida_mem.ida_eeS = eeS;
        ida_mem.ida_phiS = phiS;
        ida_mem.ida_tempv1 = tmp1;
        ida_mem.ida_tmpS3 = tmp3;

        ida_mem.ida_nrSe += 1;
        if retval < 0 {
            return IDA_RES_FAIL;
        }
        if retval > 0 {
            return IC_FAIL_RECOV;
        }

        /* Save delnewS in savresS. (N_VScale(ONE, x, z) = copy;
           delnewS = phiS[3], savresS = phiS[2]) */
        {
            let (lo, hi) = ida_mem.ida_phiS.split_at_mut(3);
            for is in 0..(ns as usize) {
                lo[2][is].data.copy_from_slice(&hi[0][is].data);
            }
        }

        /* Call the linear solve function to get J-inverse deltaS. */
        let mut lmem = std::mem::take(&mut ida_mem.ida_lmem);
        let retval = match &mut lmem {
            LsModule::Ls(idals_mem) => {
                let mut phiS = std::mem::take(&mut ida_mem.ida_phiS);
                let ewtS = std::mem::take(&mut ida_mem.ida_ewtS);
                let ynew = std::mem::take(&mut ida_mem.ida_tempv2);
                let ypnew = std::mem::take(&mut ida_mem.ida_ee);
                let savres = std::mem::take(&mut ida_mem.ida_savres);
                let mut r = 0;
                for is in 0..(ns as usize) {
                    r = idaLsSolve(ida_mem, idals_mem, &mut phiS[3][is], &ewtS[is],
                                   &ynew, &ypnew, &savres);
                    if r != 0 {
                        break;
                    }
                }
                ida_mem.ida_phiS = phiS;
                ida_mem.ida_ewtS = ewtS;
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

        /* Include sensitivities in norm. */
        *fnorm = IDASensWrmsNormUpdate(ida_mem, *fnorm, &ida_mem.ida_phiS[3],
                                       &ida_mem.ida_ewtS, SUNFALSE);
    }

    /* Rescale norm if index = 0. */
    if ida_mem.ida_sysindex == 0 {
        *fnorm *= ida_mem.ida_tscale * SUNRabs(ida_mem.ida_cj);
    }

    IDA_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * IDANewyyp
 * -----------------------------------------------------------------
 * IDANewyyp updates the vectors ynew and ypnew from yy0 and yp0,
 * using the current step vector lambda*delta, in a manner
 * depending on icopt and the input id vector.
 *
 * The return value is always IDA_SUCCESS = 0.
 * (ynew = tempv2, ypnew = ee, dtemp = phi[3].)
 * -----------------------------------------------------------------
 */
fn IDANewyyp(ida_mem: &mut IDAMem, lambda: f64) -> i32 {
    let mut retval = IDA_SUCCESS;

    /* IDA_YA_YDP_INIT case: ynew  = yy0 - lambda*delta    where id_i = 0
                             ypnew = yp0 - cj*lambda*delta where id_i = 1. */
    if ida_mem.ida_icopt == IDA_YA_YDP_INIT {
        N_VProd(&ida_mem.ida_id, &ida_mem.ida_delta, &mut ida_mem.ida_phi[3]);
        {
            let neg = -ida_mem.ida_cj * lambda;
            N_VLinearSum(ONE, &ida_mem.ida_yp0, neg, &ida_mem.ida_phi[3],
                         &mut ida_mem.ida_ee);
        }
        /* dtemp = delta - dtemp (C: N_VLinearSum(1, delta, -1, dtemp, dtemp)) */
        {
            let IDAMem { ida_phi, ida_delta, .. } = &mut *ida_mem;
            ida_phi[3].linear_sum_with(-ONE, ONE, ida_delta);
        }
        N_VLinearSum(ONE, &ida_mem.ida_yy0, -lambda, &ida_mem.ida_phi[3],
                     &mut ida_mem.ida_tempv2);
    } else if ida_mem.ida_icopt == IDA_Y_INIT {
        /* IDA_Y_INIT case: ynew = yy0 - lambda*delta. (ypnew = yp0 preset.) */
        N_VLinearSum(ONE, &ida_mem.ida_yy0, -lambda, &ida_mem.ida_delta,
                     &mut ida_mem.ida_tempv2);
    }

    if ida_mem.ida_sensi && ida_mem.ida_ism == IDA_SIMULTANEOUS {
        retval = IDASensNewyyp(ida_mem, lambda);
    }

    retval
}

/*
 * -----------------------------------------------------------------
 * IDANewy
 * -----------------------------------------------------------------
 * IDANewy updates the vector ynew from yy0,
 * using the current step vector delta, in a manner
 * depending on icopt and the input id vector.
 *
 * The return value is always IDA_SUCCESS = 0.
 * (ynew = tempv2, dtemp = phi[3].)
 * -----------------------------------------------------------------
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
        N_VLinearSum(ONE, &ida_mem.ida_yy0, -ONE, &ida_mem.ida_phi[3],
                     &mut ida_mem.ida_tempv2);
        return IDA_SUCCESS;
    }

    /* IDA_Y_INIT case: ynew = yy0 - delta. */
    N_VLinearSum(ONE, &ida_mem.ida_yy0, -ONE, &ida_mem.ida_delta,
                 &mut ida_mem.ida_tempv2);
    IDA_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Sensitivity I.C. functions
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * IDASensNlsIC
 * -----------------------------------------------------------------
 * IDASensNlsIC solves nonlinear systems for sensitivities consistent
 * initial conditions.  It mainly relies on IDASensNewtonIC.
 *
 * The return value is IDA_SUCCESS = 0 if no error occurred.
 * The error return values (positive) considered recoverable are:
 *  IC_FAIL_RECOV      if res, lsetup, or lsolve failed recoverably
 *  IC_CONSTR_FAILED   if the constraints could not be met
 *  IC_LINESRCH_FAILED if the linesearch failed (either on steptol test
 *                     or on maxbacks test)
 *  IC_CONV_FAIL       if the Newton iterations failed to converge
 *  IC_SLOW_CONVRG     if the iterations are converging slowly
 *                     (failed the convergence test, but showed
 *                     norm reduction or convergence rate < 1)
 * The error return values (negative) considered non-recoverable are:
 *  IDA_RES_FAIL       if res had a non-recoverable error
 *  IDA_FIRST_RES_FAIL if res failed recoverably on the first call
 *  IDA_LSETUP_FAIL    if lsetup had a non-recoverable error
 *  IDA_LSOLVE_FAIL    if lsolve had a non-recoverable error
 * -----------------------------------------------------------------
 */
fn IDASensNlsIC(ida_mem: &mut IDAMem) -> i32 {
    let ns = ida_mem.ida_Ns;
    let t0 = ida_mem.ida_t0;
    let yy0 = std::mem::take(&mut ida_mem.ida_yy0);
    let yp0 = std::mem::take(&mut ida_mem.ida_yp0);
    let delta = std::mem::take(&mut ida_mem.ida_delta);
    let yyS0 = std::mem::take(&mut ida_mem.ida_yyS0);
    let ypS0 = std::mem::take(&mut ida_mem.ida_ypS0);
    let mut deltaS = std::mem::take(&mut ida_mem.ida_deltaS);
    /* tmpS1 = tempv1, tmpS2 = tempv2 (aliases), tmpS3 real */
    let mut tmp1 = std::mem::take(&mut ida_mem.ida_tempv1);
    let mut tmp2 = std::mem::take(&mut ida_mem.ida_tempv2);
    let mut tmp3 = std::mem::take(&mut ida_mem.ida_tmpS3);

    let mut retval = if ida_mem.ida_resSDQ {
        IDASensResDQ(ida_mem, ns, t0, &yy0, &yp0, &delta, &yyS0, &ypS0, &mut deltaS,
                     &mut tmp1, &mut tmp2, &mut tmp3)
    } else {
        let resS = ida_mem.ida_resS.unwrap();
        resS(ns, t0, &yy0, &yp0, &delta, &yyS0, &ypS0, &mut deltaS,
             &mut ida_mem.ida_user_data, &mut tmp1, &mut tmp2, &mut tmp3)
    };

    ida_mem.ida_yy0 = yy0;
    ida_mem.ida_yp0 = yp0;
    ida_mem.ida_delta = delta;
    ida_mem.ida_yyS0 = yyS0;
    ida_mem.ida_ypS0 = ypS0;
    ida_mem.ida_deltaS = deltaS;
    ida_mem.ida_tempv1 = tmp1;
    ida_mem.ida_tempv2 = tmp2;
    ida_mem.ida_tmpS3 = tmp3;

    ida_mem.ida_nrSe += 1;
    if retval < 0 {
        return IDA_RES_FAIL;
    }
    if retval > 0 {
        return IDA_FIRST_RES_FAIL;
    }

    /* Save deltaS (savresS = phiS[2]) */
    for is in 0..(ns as usize) {
        N_VScale(ONE, &ida_mem.ida_deltaS[is], &mut ida_mem.ida_phiS[2][is]);
    }

    /* Loop over nj = number of linear solve Jacobian setups. */
    for nj in 1..=2 {
        /* Call the Newton iteration routine */
        retval = IDASensNewtonIC(ida_mem);
        if retval == IDA_SUCCESS {
            return IDA_SUCCESS;
        }

        /* If converging slowly and lsetup is nontrivial and this is the first pass,
           update Jacobian and retry. */
        if retval == IC_SLOW_CONVRG && ida_has_lsetup(ida_mem) && nj == 1 {
            /* Restore deltaS. (savresS = phiS[2]) */
            for is in 0..(ns as usize) {
                N_VScale(ONE, &ida_mem.ida_phiS[2][is], &mut ida_mem.ida_deltaS[is]);
            }

            ida_mem.ida_nsetupsS += 1;
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

            continue;
        } else {
            return retval;
        }
    }

    IDA_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * IDASensNewtonIC
 * -----------------------------------------------------------------
 * IDANewtonIC performs the Newton iteration to solve for
 * sensitivities consistent initial conditions.  It calls
 * IDASensLineSrch within each iteration.
 * On return, savresS contains the current residual vectors.
 *
 * The return value is IDA_SUCCESS = 0 if no error occurred.
 * The error return values (positive) considered recoverable are:
 *  IC_FAIL_RECOV      if res or lsolve failed recoverably
 *  IC_CONSTR_FAILED   if the constraints could not be met
 *  IC_LINESRCH_FAILED if the linesearch failed (either on steptol test
 *                     or on maxbacks test)
 *  IC_CONV_FAIL       if the Newton iterations failed to converge
 *  IC_SLOW_CONVRG     if the iterations appear to be converging slowly.
 *                     They failed the convergence test, but showed
 *                     an overall norm reduction (by a factor of < 0.1)
 *                     or a convergence rate <= ICRATEMAX).
 * The error return values (negative) considered non-recoverable are:
 *  IDA_RES_FAIL   if res had a non-recoverable error
 *  IDA_LSOLVE_FAIL      if lsolve had a non-recoverable error
 * -----------------------------------------------------------------
 */
fn IDASensNewtonIC(ida_mem: &mut IDAMem) -> i32 {
    let ns = ida_mem.ida_Ns as usize;

    /* Call the linear solve function to get the Newton step, delta.
       (NOTE: C passes ida_delta as the residual argument here.) */
    let mut lmem = std::mem::take(&mut ida_mem.ida_lmem);
    let retval = match &mut lmem {
        LsModule::Ls(idals_mem) => {
            let mut deltaS = std::mem::take(&mut ida_mem.ida_deltaS);
            let ewtS = std::mem::take(&mut ida_mem.ida_ewtS);
            let yy0 = std::mem::take(&mut ida_mem.ida_yy0);
            let yp0 = std::mem::take(&mut ida_mem.ida_yp0);
            let delta = std::mem::take(&mut ida_mem.ida_delta);
            let mut r = 0;
            for is in 0..ns {
                r = idaLsSolve(ida_mem, idals_mem, &mut deltaS[is], &ewtS[is], &yy0,
                               &yp0, &delta);
                if r != 0 {
                    break;
                }
            }
            ida_mem.ida_deltaS = deltaS;
            ida_mem.ida_ewtS = ewtS;
            ida_mem.ida_yy0 = yy0;
            ida_mem.ida_yp0 = yp0;
            ida_mem.ida_delta = delta;
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

    /* Compute the norm of the step and return if it is small enough */
    let mut fnorm = IDASensWrmsNorm(ida_mem, &ida_mem.ida_deltaS, &ida_mem.ida_ewtS,
                                    SUNFALSE);
    if ida_mem.ida_sysindex == 0 {
        fnorm *= ida_mem.ida_tscale * SUNRabs(ida_mem.ida_cj);
    }
    if fnorm <= ida_mem.ida_epsNewt {
        return IDA_SUCCESS;
    }
    let fnorm0 = fnorm;

    let mut rate = ZERO;

    /* Newton iteration loop */
    for _mnewt in 0..ida_mem.ida_maxnit {
        ida_mem.ida_nniS += 1;
        let mut delnorm = fnorm;
        let oldfnrm = fnorm;

        /* Call the Linesearch function and return if it failed. */
        let retval = IDASensLineSrch(ida_mem, &mut delnorm, &mut fnorm);
        if retval != IDA_SUCCESS {
            return retval;
        }

        /* Set the observed convergence rate and test for convergence. */
        rate = fnorm / oldfnrm;
        if fnorm <= ida_mem.ida_epsNewt {
            return IDA_SUCCESS;
        }

        /* If not converged, copy new step vectors, and loop. (delnewS = phiS[3]) */
        for is in 0..ns {
            N_VScale(ONE, &ida_mem.ida_phiS[3][is], &mut ida_mem.ida_deltaS[is]);
        }
    } /* End of Newton iteration loop */

    /* Return either IC_SLOW_CONVRG or recoverable fail flag. */
    if rate <= ICRATEMAX || fnorm < PT1 * fnorm0 {
        return IC_SLOW_CONVRG;
    }
    IC_CONV_FAIL
}

/*
 * -----------------------------------------------------------------
 * IDASensLineSrch
 * -----------------------------------------------------------------
 * IDASensLineSrch performs the Linesearch algorithm with the
 * calculation of consistent initial conditions for sensitivities
 * systems.
 *
 * On entry, yyS0 and ypS0 contain the current values, the Newton
 * steps are contained in deltaS, the current residual vectors FS are
 * savresS, delnorm is sens-WRMS-norm(deltaS), and fnorm is
 * max { WRMS-norm( J-inverse FS[is] ) : is=1,2,...,Ns }
 *
 * On a successful return, yy0, yp0, and savres have been updated,
 * delnew contains the current values of J-inverse FS, and fnorm is
 * max { WRMS-norm(delnewS[is]) : is = 1,2,...Ns }
 *
 * The return value is IDA_SUCCESS = 0 if no error occurred.
 * The error return values (positive) considered recoverable are:
 *  IC_FAIL_RECOV      if res or lsolve failed recoverably
 *  IC_CONSTR_FAILED   if the constraints could not be met
 *  IC_LINESRCH_FAILED if the linesearch failed (either on steptol test
 *                     or on maxbacks test)
 * The error return values (negative) considered non-recoverable are:
 *  IDA_RES_FAIL   if res had a non-recoverable error
 *  IDA_LSOLVE_FAIL      if lsolve had a non-recoverable error
 * -----------------------------------------------------------------
 */
fn IDASensLineSrch(ida_mem: &mut IDAMem, delnorm: &mut f64, fnorm: &mut f64) -> i32 {
    /* Set work space pointer. (dtemp = phi[3]) */

    let f1norm = (*fnorm) * (*fnorm) * HALF;

    /* Initialize local variables. */
    let ratio = ONE;
    let slpi = -TWO * f1norm * ratio;
    let minlam = ida_mem.ida_steptol / (*delnorm);
    let mut lambda = ONE;
    let mut nbacks = 0;

    let mut fnormp = ZERO;
    loop {
        if nbacks == ida_mem.ida_maxbacks {
            return IC_LINESRCH_FAILED;
        }
        /* Get new iteration in (ySnew, ypSnew). */
        IDASensNewyyp(ida_mem, lambda);

        /* Get the norm of new function value. */
        let retval = IDASensfnorm(ida_mem, &mut fnormp);
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
    }

    /* Update yyS0, ypS0 and fnorm and return. (yyS0new = phiS[4]) */
    for is in 0..(ida_mem.ida_Ns as usize) {
        N_VScale(ONE, &ida_mem.ida_phiS[4][is], &mut ida_mem.ida_yyS0[is]);
    }

    if ida_mem.ida_icopt == IDA_YA_YDP_INIT {
        /* (ypS0new = eeS) */
        for is in 0..(ida_mem.ida_Ns as usize) {
            N_VScale(ONE, &ida_mem.ida_eeS[is], &mut ida_mem.ida_ypS0[is]);
        }
    }

    *fnorm = fnormp;
    IDA_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * IDASensfnorm
 * -----------------------------------------------------------------
 * IDASensfnorm computes the norm of the current function value, by
 * evaluating the sensitivity residual function, calling the linear
 * system solver, and computing a WRMS-norm.
 *
 * On return, savresS contains the current residual vectors FS, and
 * delnewS contains J-inverse FS.
 *
 * The return value is IDA_SUCCESS = 0 if no error occurred, or
 *  IC_FAIL_RECOV    if res or lsolve failed recoverably, or
 *  IDA_RES_FAIL     if res had a non-recoverable error, or
 *  IDA_LSOLVE_FAIL  if lsolve had a non-recoverable error.
 * -----------------------------------------------------------------
 */
fn IDASensfnorm(ida_mem: &mut IDAMem, fnorm: &mut f64) -> i32 {
    /* Get sensitivity residual.
       (yyS0new = phiS[4], ypS0new = eeS, delnewS = phiS[3].
       Unlike IDAfnorm, there is no yy/scratch aliasing here: the
       state arguments are yy0/yp0, so tmpS2 = tempv2 detaches
       normally.) */
    let ns = ida_mem.ida_Ns;
    let t0 = ida_mem.ida_t0;
    let yy0 = std::mem::take(&mut ida_mem.ida_yy0);
    let yp0 = std::mem::take(&mut ida_mem.ida_yp0);
    let delta = std::mem::take(&mut ida_mem.ida_delta);
    let eeS = std::mem::take(&mut ida_mem.ida_eeS);
    let mut phiS = std::mem::take(&mut ida_mem.ida_phiS);
    let mut tmp1 = std::mem::take(&mut ida_mem.ida_tempv1);
    let mut tmp2 = std::mem::take(&mut ida_mem.ida_tempv2);
    let mut tmp3 = std::mem::take(&mut ida_mem.ida_tmpS3);

    let retval = {
        let (lo, hi) = phiS.split_at_mut(4);
        let yyS0new = &hi[0];
        let delnewS = &mut lo[3];
        if ida_mem.ida_resSDQ {
            IDASensResDQ(ida_mem, ns, t0, &yy0, &yp0, &delta, yyS0new, &eeS, delnewS,
                         &mut tmp1, &mut tmp2, &mut tmp3)
        } else {
            let resS = ida_mem.ida_resS.unwrap();
            resS(ns, t0, &yy0, &yp0, &delta, yyS0new, &eeS, delnewS,
                 &mut ida_mem.ida_user_data, &mut tmp1, &mut tmp2, &mut tmp3)
        }
    };

    ida_mem.ida_yy0 = yy0;
    ida_mem.ida_yp0 = yp0;
    ida_mem.ida_delta = delta;
    ida_mem.ida_eeS = eeS;
    ida_mem.ida_phiS = phiS;
    ida_mem.ida_tempv1 = tmp1;
    ida_mem.ida_tempv2 = tmp2;
    ida_mem.ida_tmpS3 = tmp3;

    ida_mem.ida_nrSe += 1;
    if retval < 0 {
        return IDA_RES_FAIL;
    }
    if retval > 0 {
        return IC_FAIL_RECOV;
    }

    /* Save delnewS in savresS. (N_VScale(ONE, x, z) = copy;
       delnewS = phiS[3], savresS = phiS[2]) */
    {
        let (lo, hi) = ida_mem.ida_phiS.split_at_mut(3);
        for is in 0..(ns as usize) {
            lo[2][is].data.copy_from_slice(&hi[0][is].data);
        }
    }

    /* Call linear solve function.
       (NOTE: C passes ida_delta as the residual argument here.) */
    let mut lmem = std::mem::take(&mut ida_mem.ida_lmem);
    let retval = match &mut lmem {
        LsModule::Ls(idals_mem) => {
            let mut phiS = std::mem::take(&mut ida_mem.ida_phiS);
            let ewtS = std::mem::take(&mut ida_mem.ida_ewtS);
            let yy0 = std::mem::take(&mut ida_mem.ida_yy0);
            let yp0 = std::mem::take(&mut ida_mem.ida_yp0);
            let delta = std::mem::take(&mut ida_mem.ida_delta);
            let mut r = 0;
            for is in 0..(ns as usize) {
                r = idaLsSolve(ida_mem, idals_mem, &mut phiS[3][is], &ewtS[is], &yy0,
                               &yp0, &delta);
                if r != 0 {
                    break;
                }
            }
            ida_mem.ida_phiS = phiS;
            ida_mem.ida_ewtS = ewtS;
            ida_mem.ida_yy0 = yy0;
            ida_mem.ida_yp0 = yp0;
            ida_mem.ida_delta = delta;
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

    /* Compute the WRMS-norm; rescale if index = 0. (delnewS = phiS[3]) */
    *fnorm = IDASensWrmsNorm(ida_mem, &ida_mem.ida_phiS[3], &ida_mem.ida_ewtS, SUNFALSE);
    if ida_mem.ida_sysindex == 0 {
        *fnorm *= ida_mem.ida_tscale * SUNRabs(ida_mem.ida_cj);
    }

    IDA_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * IDASensNewyyp
 * -----------------------------------------------------------------
 * IDASensNewyyp computes the Newton updates for each of the
 * sensitivities systems using the current step vector lambda*delta,
 * in a manner depending on icopt and the input id vector.
 *
 * The return value is always IDA_SUCCESS = 0.
 * (dtemp = phi[3], yyS0new = phiS[4], ypS0new = eeS.)
 * -----------------------------------------------------------------
 */
fn IDASensNewyyp(ida_mem: &mut IDAMem, lambda: f64) -> i32 {
    let ns = ida_mem.ida_Ns as usize;

    if ida_mem.ida_icopt == IDA_YA_YDP_INIT {
        /* IDA_YA_YDP_INIT case:
         - ySnew  = yS0  - lambda*deltaS    where id_i = 0
         - ypSnew = ypS0 - cj*lambda*delta  where id_i = 1. */

        for is in 0..ns {
            /* It is ok to use dtemp as temporary vector here. */
            N_VProd(&ida_mem.ida_id, &ida_mem.ida_deltaS[is], &mut ida_mem.ida_phi[3]);
            {
                let neg = -ida_mem.ida_cj * lambda;
                N_VLinearSum(ONE, &ida_mem.ida_ypS0[is], neg, &ida_mem.ida_phi[3],
                             &mut ida_mem.ida_eeS[is]);
            }
            /* dtemp = deltaS[is] - dtemp
               (C: N_VLinearSum(1, deltaS[is], -1, dtemp, dtemp)) */
            {
                let IDAMem { ida_phi, ida_deltaS, .. } = &mut *ida_mem;
                ida_phi[3].linear_sum_with(-ONE, ONE, &ida_deltaS[is]);
            }
            N_VLinearSum(ONE, &ida_mem.ida_yyS0[is], -lambda, &ida_mem.ida_phi[3],
                         &mut ida_mem.ida_phiS[4][is]);
        } /* end loop is */
    } else {
        /* IDA_Y_INIT case:
           - ySnew = yS0 - lambda*deltaS. (ypnew = yp0 preset.) */

        for is in 0..ns {
            N_VLinearSum(ONE, &ida_mem.ida_yyS0[is], -lambda, &ida_mem.ida_deltaS[is],
                         &mut ida_mem.ida_phiS[4][is]);
        }
    } /* end loop is */
    IDA_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * IDAICFailFlag
 * -----------------------------------------------------------------
 * IDAICFailFlag prints a message and sets the IDACalcIC return
 * value appropriate to the flag retval returned by IDANlsIC.
 * -----------------------------------------------------------------
 */
fn IDAICFailFlag(ida_mem: &mut IDAMem, retval: i32) -> i32 {
    /* Depending on retval, print error message and return error flag. */
    match retval {
        IDA_RES_FAIL => {
            IDAProcessError(Some(ida_mem), IDA_RES_FAIL, line!(), "IDAICFailFlag", file!(),
                            MSG_IC_RES_NONREC);
            IDA_RES_FAIL
        }
        IDA_FIRST_RES_FAIL => {
            IDAProcessError(Some(ida_mem), IDA_FIRST_RES_FAIL, line!(), "IDAICFailFlag",
                            file!(), MSG_IC_RES_FAIL);
            IDA_FIRST_RES_FAIL
        }
        IDA_LSETUP_FAIL => {
            IDAProcessError(Some(ida_mem), IDA_LSETUP_FAIL, line!(), "IDAICFailFlag",
                            file!(), MSG_IC_SETUP_FAIL);
            IDA_LSETUP_FAIL
        }
        IDA_LSOLVE_FAIL => {
            IDAProcessError(Some(ida_mem), IDA_LSOLVE_FAIL, line!(), "IDAICFailFlag",
                            file!(), MSG_IC_SOLVE_FAIL);
            IDA_LSOLVE_FAIL
        }
        IC_FAIL_RECOV => {
            IDAProcessError(Some(ida_mem), IDA_NO_RECOVERY, line!(), "IDAICFailFlag",
                            file!(), MSG_IC_NO_RECOVERY);
            IDA_NO_RECOVERY
        }
        IC_CONSTR_FAILED => {
            IDAProcessError(Some(ida_mem), IDA_CONSTR_FAIL, line!(), "IDAICFailFlag",
                            file!(), MSG_IC_FAIL_CONSTR);
            IDA_CONSTR_FAIL
        }
        IC_LINESRCH_FAILED => {
            IDAProcessError(Some(ida_mem), IDA_LINESEARCH_FAIL, line!(), "IDAICFailFlag",
                            file!(), MSG_IC_FAILED_LINS);
            IDA_LINESEARCH_FAIL
        }
        IC_CONV_FAIL => {
            IDAProcessError(Some(ida_mem), IDA_CONV_FAIL, line!(), "IDAICFailFlag",
                            file!(), MSG_IC_CONV_FAILED);
            IDA_CONV_FAIL
        }
        IC_SLOW_CONVRG => {
            IDAProcessError(Some(ida_mem), IDA_CONV_FAIL, line!(), "IDAICFailFlag",
                            file!(), MSG_IC_CONV_FAILED);
            IDA_CONV_FAIL
        }
        IDA_BAD_EWT => {
            IDAProcessError(Some(ida_mem), IDA_BAD_EWT, line!(), "IDAICFailFlag", file!(),
                            MSG_IC_BAD_EWT);
            IDA_BAD_EWT
        }
        _ => -99,
    }
}

/* END of idas_ic.c port. */
