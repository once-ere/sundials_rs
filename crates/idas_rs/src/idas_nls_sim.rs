/* -----------------------------------------------------------------
 * Translated from src/idas/idas_nls_sim.c (IDAS 7.7.0), together
 * with the Newton solve driver of sunnonlinsol_newton.c that the C
 * code reaches through the SUNNonlinearSolver ops table: the
 * SIMULTANEOUS corrector, which iterates on the composite vector
 * [ee, eeS] (state + all Ns sensitivity corrections at once).
 *
 * In C the composite vectors are N_Vector senswrappers of Ns+1
 * sub-vectors whose slots alias IDAMem storage:
 *   ypredictSim = [ida_yypredict, ida_yySpredict[0..Ns]]
 *   ycorSim     = [ida_ee,        ida_eeS[0..Ns]]
 *   ewtSim      = [ida_ewt,       ida_ewtS[0..Ns]]
 * Per the pinned senswrapper decision (idas_impl.rs, mirroring
 * cvodes_nls_sim.rs) those aliases are NOT stored: this module
 * operates directly on the IDAMem fields, and only the solver-owned
 * workspace (NewtonSolver.deltaS — created by SUNNonlinSol_NewtonSens
 * with Ns+1 sub-vectors: sub-vector 0 is the state part, sub-vectors
 * 1..=Ns the sensitivities) remains a real senswrapper.  The
 * cross-sub-vector reduction semantics of the C senswrapper are
 * reproduced exactly where composite norms are needed:
 * N_VWrmsNorm(wrapper) = MAX of the per-sub-vector WRMS norms
 * (init 0, state sub-vector first).
 *
 * Exported entry points for idas_nls.rs / idas.rs:
 *   IDASetNonlinearSolverSensSim / IDAGetNonlinearSystemDataSens
 *   idaNlsInitSensSim(ida_mem: &mut IDAMem) -> i32
 *   idaNlsSolveSensSim(ida_mem, nls, tol, callLSetup) -> i32, the
 *       specialization of SUNNonlinSolSolve(NLSsim, ypredictSim,
 *       ycorSim, ewtSim, tol, callLSetup, IDA_mem) invoked by
 *       IDANls's sensi_sim branch (the caller detaches the solver
 *       via NLSsim.take() and reads the niters/nconvfails counters).
 * -----------------------------------------------------------------*/
use crate::idas_impl::*;
use crate::idas_ls::{idaLsSetup, idaLsSolve};
use crate::idas_nls::ida_has_lsetup;
use crate::nvector_serial::*;
use crate::sundials_errors::SUN_SUCCESS;
use crate::sundials_math::*;
use crate::sundials_nonlinearsolver::*;
use crate::sundials_types::*;
use crate::sunnonlinsol_newton::NewtonSolver;
use sundials_core::sundials_nvector_senswrapper::NVectorSensWrapper;

/* constant macros */
const PT0001: f64 = 0.0001; /* real 0.0001 */
const ZERO: f64 = 0.0; /* real 0.0 (sunnonlinsol_newton.c) */
const ONE: f64 = 1.0; /* real 1.0    */
const TWENTY: f64 = 20.0; /* real 20.0   */

/* nonlinear solver parameters */
const MAXIT: i32 = 4; /* default max number of nonlinear iterations    */
const RATEMAX: f64 = 0.9; /* max convergence rate used in divergence check */

/* -----------------------------------------------------------------------------
 * Exported functions
 * ---------------------------------------------------------------------------*/

pub fn IDASetNonlinearSolverSensSim(ida_mem: &mut IDAMem, NLS: NonlinearSolver) -> i32 {
    /* (The C NULL-input and missing-ops checks — gettype/solve/setsysfn
       on the NLS — cannot fail here: the workspace enum implements
       them.) */

    /* check for allowed nonlinear solver types */
    if NLS.nls_type() != SUNNONLINEARSOLVER_ROOTFIND {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASetNonlinearSolverSensSim",
                        file!(), "NLS type must be SUNNONLINEARSOLVER_ROOTFIND");
        return IDA_ILL_INPUT;
    }

    /* check that sensitivities were initialized */
    if !ida_mem.ida_sensi {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASetNonlinearSolverSensSim",
                        file!(), MSG_NO_SENSI);
        return IDA_ILL_INPUT;
    }

    /* check that the simultaneous corrector was selected */
    if ida_mem.ida_ism != IDA_SIMULTANEOUS {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASetNonlinearSolverSensSim",
                        file!(), "Sensitivity solution method is not IDA_SIMULTANEOUS");
        return IDA_ILL_INPUT;
    }

    /* free any existing nonlinear solver: RAII (drop on overwrite) */

    /* set SUNNonlinearSolver pointer */
    ida_mem.NLSsim = Some(NLS);

    /* Set NLS ownership flag. If this function was called to attach the default
       NLS, IDA will set the flag to SUNTRUE after this function returns. */
    ida_mem.ownNLSsim = SUNFALSE;

    /* the nonlinear residual and convergence test functions
       (idaNlsResidualSensSim / idaNlsConvTestSensSim) are wired
       statically into the Newton solve driver below; the C
       SUNNonlinSolSetSysFn / SUNNonlinSolSetConvTestFn registrations
       have no Rust counterpart */

    /* set max allowed nonlinear iterations */
    let retval = ida_mem.NLSsim.as_mut().unwrap().set_max_iters(MAXIT);
    if retval != IDA_SUCCESS {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASetNonlinearSolverSensSim",
                        file!(), "Setting maximum number of nonlinear iterations failed");
        return IDA_ILL_INPUT;
    }

    /* create vector wrappers if necessary: per the pinned senswrapper
       decision the ypredictSim/ycorSim/ewtSim aliases are not stored —
       this module reads ida_yypredict/ida_yySpredict, ida_ee/ida_eeS
       and ida_ewt/ida_ewtS directly — so only the allocation flag
       survives (the C attach-loop below it vanishes with them) */
    if !ida_mem.simMallocDone {
        ida_mem.simMallocDone = SUNTRUE;
    }

    /* Set the nonlinear system RES function */
    if ida_mem.ida_res.is_none() {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASetNonlinearSolverSensSim",
                        file!(), "The DAE residual function is NULL");
        return IDA_ILL_INPUT;
    }
    ida_mem.nls_res = ida_mem.ida_res;

    IDA_SUCCESS
}

/*---------------------------------------------------------------
  IDAGetNonlinearSystemDataSens:

  This routine provides access to the relevant data needed to
  compute the nonlinear system function (out-pointers become a
  returned tuple (tcur, cj), donor convention; the vector arrays —
  yySpredict, ypSpredict, yyS, ypS — and user_data remain
  accessible as IDAMem fields).
  ---------------------------------------------------------------*/
pub fn IDAGetNonlinearSystemDataSens(ida_mem: &IDAMem) -> (f64, f64) {
    (ida_mem.ida_tn, ida_mem.ida_cj)
}

/* -----------------------------------------------------------------------------
 * Private functions
 * ---------------------------------------------------------------------------*/

pub fn idaNlsInitSensSim(ida_mem: &mut IDAMem) -> i32 {
    /* In C this wires the idaNlsLSetupSensSim/idaNlsLSolveSensSim
       wrapper functions into the NLS depending on whether
       ida_lsetup/ida_lsolve exist (infallible registrations); here the
       dispatch is dynamic through ida_has_lsetup / the LsModule enum. */

    /* initialize nonlinear solver */
    let retval = match ida_mem.NLSsim.as_mut() {
        Some(nls) => nls.initialize(),
        None => -1,
    };

    if retval != IDA_SUCCESS {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "idaNlsInitSensSim", file!(),
                        MSG_NLS_INIT_FAIL);
        return IDA_NLS_INIT_FAIL;
    }

    IDA_SUCCESS
}

/* idaNlsLSetupSensSim (idas_nls_sim.c): wrapper around the lsetup
   dispatch.  In C the lsetup call receives (yy, yp, savres,
   tempv1..3); the vectors alias IDAMem fields, so they are detached
   for the call (the three tmps come from the IDAMem fields inside
   idaLsSetup). */
fn idaNlsLSetupSensSim(ida_mem: &mut IDAMem, _jbad: bool, jcur: &mut bool) -> i32 {
    ida_mem.ida_nsetups += 1;
    ida_mem.ida_forceSetup = SUNFALSE;

    let mut lmem = std::mem::take(&mut ida_mem.ida_lmem);
    let retval = match &mut lmem {
        LsModule::Ls(idals_mem) => {
            let yy = std::mem::take(&mut ida_mem.ida_yy);
            let yp = std::mem::take(&mut ida_mem.ida_yp);
            let savres = std::mem::take(&mut ida_mem.ida_savres);
            let r = idaLsSetup(ida_mem, idals_mem, &yy, &yp, &savres);
            ida_mem.ida_yy = yy;
            ida_mem.ida_yp = yp;
            ida_mem.ida_savres = savres;
            r
        }
        /* cannot arise: dispatch is guarded by ida_has_lsetup */
        LsModule::None => 0,
    };
    ida_mem.ida_lmem = lmem;

    /* update Jacobian status */
    *jcur = SUNTRUE;

    /* update convergence test constants */
    ida_mem.ida_cjold = ida_mem.ida_cj;
    ida_mem.ida_cjratio = ONE;
    ida_mem.ida_ss = TWENTY;
    ida_mem.ida_ssS = TWENTY;

    if retval < 0 {
        return IDA_LSETUP_FAIL;
    }
    if retval > 0 {
        return IDA_LSETUP_RECVR;
    }

    IDA_SUCCESS
}

/* idaNlsLSolveSensSim (idas_nls_sim.c): solves the state linear
   system (weight = ida_ewt) followed by the Ns sensitivity linear
   systems (weight = ida_ewtS[is]); deltaSim sub-vector 0 = state
   delta, sub-vectors 1..=Ns = sensitivity deltas.  In C each lsolve
   call receives (delta, w, yy, yp, savres); the weight and
   current-state vectors alias IDAMem fields, so they are detached
   for the calls.  (C maps every failing call's retval identically,
   so mapping the last nonzero retval after the restore block is
   behaviorally exact.) */
fn idaNlsLSolveSensSim(ida_mem: &mut IDAMem, deltaSim: &mut NVectorSensWrapper) -> i32 {
    let mut lmem = std::mem::take(&mut ida_mem.ida_lmem);
    let retval = match &mut lmem {
        LsModule::Ls(idals_mem) => {
            let ewt = std::mem::take(&mut ida_mem.ida_ewt);
            let ewtS = std::mem::take(&mut ida_mem.ida_ewtS);
            let yy = std::mem::take(&mut ida_mem.ida_yy);
            let yp = std::mem::take(&mut ida_mem.ida_yp);
            let savres = std::mem::take(&mut ida_mem.ida_savres);

            /* solve the state linear system */
            let mut r = idaLsSolve(ida_mem, idals_mem, &mut deltaSim.vecs[0], &ewt, &yy, &yp,
                                   &savres);

            /* solve the sensitivity linear systems */
            if r == 0 {
                for is in 0..ida_mem.ida_Ns as usize {
                    r = idaLsSolve(ida_mem, idals_mem, &mut deltaSim.vecs[is + 1], &ewtS[is],
                                   &yy, &yp, &savres);
                    if r != 0 {
                        break;
                    }
                }
            }

            ida_mem.ida_ewt = ewt;
            ida_mem.ida_ewtS = ewtS;
            ida_mem.ida_yy = yy;
            ida_mem.ida_yp = yp;
            ida_mem.ida_savres = savres;
            r
        }
        /* cannot arise: IDAInitialSetup rejects a missing lsolve */
        LsModule::None => -1,
    };
    ida_mem.ida_lmem = lmem;

    if retval < 0 {
        return IDA_LSOLVE_FAIL;
    }
    if retval > 0 {
        return IDA_LSOLVE_RECVR;
    }

    IDA_SUCCESS
}

/* idaNlsResidualSensSim (idas_nls_sim.c): updates yy/yp and yyS/ypS
   from the current composite correction ycorSim = [ida_ee, ida_eeS]
   (read directly from IDAMem — aliases not stored) and evaluates the
   DAE residual (into resSim sub-vector 0, saved in savres) followed
   by the Ns sensitivity residuals (into sub-vectors 1..=Ns). */
fn idaNlsResidualSensSim(ida_mem: &mut IDAMem, resSim: &mut NVectorSensWrapper) -> i32 {
    /* update yy and yp based on the current correction */
    {
        let IDAMem { ida_yypredict, ida_ee, ida_yy, .. } = ida_mem;
        N_VLinearSum(ONE, ida_yypredict, ONE, ida_ee, ida_yy);
    }
    {
        let IDAMem { ida_yppredict, ida_ee, ida_yp, ida_cj, .. } = ida_mem;
        N_VLinearSum(ONE, ida_yppredict, *ida_cj, ida_ee, ida_yp);
    }

    /* evaluate residual */
    let res_fn = ida_mem.nls_res.unwrap();
    let retval = {
        let IDAMem { ida_tn, ida_yy, ida_yp, ida_user_data, .. } = ida_mem;
        res_fn(*ida_tn, ida_yy, ida_yp, &mut resSim.vecs[0], ida_user_data)
    };

    /* increment the number of residual evaluations */
    ida_mem.ida_nre += 1;

    /* save a copy of the residual vector in savres */
    N_VScale(ONE, &resSim.vecs[0], &mut ida_mem.ida_savres);

    if retval < 0 {
        return IDA_RES_FAIL;
    }
    if retval > 0 {
        return IDA_RES_RECVR;
    }

    /* update yS and ypS based on the current correction
       (N_VLinearSumVectorArray expanded to per-vector loops) */
    {
        let ns = ida_mem.ida_Ns as usize;
        let IDAMem { ida_yySpredict, ida_eeS, ida_yyS, .. } = ida_mem;
        for is in 0..ns {
            N_VLinearSum(ONE, &ida_yySpredict[is], ONE, &ida_eeS[is], &mut ida_yyS[is]);
        }
    }
    {
        let ns = ida_mem.ida_Ns as usize;
        let IDAMem { ida_ypSpredict, ida_eeS, ida_ypS, ida_cj, .. } = ida_mem;
        for is in 0..ns {
            N_VLinearSum(ONE, &ida_ypSpredict[is], *ida_cj, &ida_eeS[is], &mut ida_ypS[is]);
        }
    }

    /* evaluate sens residual (C: ida_resS(Ns, tn, yy, yp, res, yyS, ypS,
       resS, user_dataS, tmpS1, tmpS2, tmpS3); pinned dispatch — resSDQ
       selects the internal DQ routine, which receives IDA_mem in place
       of the user_dataS self-pointer; tmpS1/tmpS2 alias tempv1/tempv2.
       The argument vectors are IDAMem fields, taken out for the
       duration of the call — donor take()/restore pattern.) */
    let ns = ida_mem.ida_Ns;
    let tn = ida_mem.ida_tn;
    let yy = std::mem::take(&mut ida_mem.ida_yy);
    let yp = std::mem::take(&mut ida_mem.ida_yp);
    let yyS = std::mem::take(&mut ida_mem.ida_yyS);
    let ypS = std::mem::take(&mut ida_mem.ida_ypS);
    let mut tmp1 = std::mem::take(&mut ida_mem.ida_tempv1);
    let mut tmp2 = std::mem::take(&mut ida_mem.ida_tempv2);
    let mut tmp3 = std::mem::take(&mut ida_mem.ida_tmpS3);

    let retval = {
        let (res, resS) = resSim.vecs.split_at_mut(1);
        if ida_mem.ida_resSDQ {
            crate::idas::IDASensResDQ(ida_mem, ns, tn, &yy, &yp, &res[0], &yyS, &ypS, resS,
                                      &mut tmp1, &mut tmp2, &mut tmp3)
        } else {
            let resS_fn = ida_mem.ida_resS.unwrap();
            resS_fn(ns, tn, &yy, &yp, &res[0], &yyS, &ypS, resS, &mut ida_mem.ida_user_data,
                    &mut tmp1, &mut tmp2, &mut tmp3)
        }
    };

    ida_mem.ida_yy = yy;
    ida_mem.ida_yp = yp;
    ida_mem.ida_yyS = yyS;
    ida_mem.ida_ypS = ypS;
    ida_mem.ida_tempv1 = tmp1;
    ida_mem.ida_tempv2 = tmp2;
    ida_mem.ida_tmpS3 = tmp3;

    /* increment the number of sens residual evaluations */
    ida_mem.ida_nrSe += 1;

    if retval < 0 {
        return IDA_SRES_FAIL;
    }
    if retval > 0 {
        return IDA_SRES_RECVR;
    }

    IDA_SUCCESS
}

/* idaNlsConvTestSensSim (idas_nls_sim.c).  m is the current nonlinear
   iteration count (C SUNNonlinSolGetCurIter); del/ewt are the C
   senswrapper composites — delSim is the solver-owned workspace and
   ewtSim = [ida_ewt, ida_ewtS] is read directly from IDAMem; the C
   ycor input is unused. */
fn idaNlsConvTestSensSim(
    ida_mem: &mut IDAMem,
    delSim: &NVectorSensWrapper,
    tol: f64,
    m: i32,
) -> i32 {
    /* compute the norm of the correction: N_VWrmsNorm(del, ewt) on
       senswrappers = MAX of the per-sub-vector WRMS norms (init 0,
       state sub-vector first) */
    let mut delnrm = ZERO;
    let tmp = N_VWrmsNorm(&delSim.vecs[0], &ida_mem.ida_ewt);
    if tmp > delnrm {
        delnrm = tmp;
    }
    for is in 0..ida_mem.ida_Ns as usize {
        let tmp = N_VWrmsNorm(&delSim.vecs[is + 1], &ida_mem.ida_ewtS[is]);
        if tmp > delnrm {
            delnrm = tmp;
        }
    }

    /* test for convergence, first directly, then with rate estimate. */
    if m == 0 {
        ida_mem.ida_oldnrm = delnrm;
        if delnrm <= PT0001 * ida_mem.ida_toldel {
            return SUN_SUCCESS;
        }
    } else {
        let rate = SUNRpowerR(delnrm / ida_mem.ida_oldnrm, ONE / (m as f64));
        if rate > RATEMAX {
            return SUN_NLS_CONV_RECVR;
        }
        ida_mem.ida_ss = rate / (ONE - rate);
    }

    if ida_mem.ida_ss * delnrm <= tol {
        return SUN_SUCCESS;
    }

    /* not yet converged */
    SUN_NLS_CONTINUE
}

/*
 * idaNlsSolveSensSim — SUNNonlinSolSolve(IDA_mem->NLSsim,
 * IDA_mem->ypredictSim, IDA_mem->ycorSim, IDA_mem->ewtSim, tol,
 * callLSetup, IDA_mem) as invoked by idas.c's IDANls (idas_nls.rs).
 * The caller detaches the solver from IDAMem (NLSsim.take()) and
 * reattaches it afterwards, reading the niters/nconvfails counters
 * for the nni/nnf updates.
 */
pub fn idaNlsSolveSensSim(
    ida_mem: &mut IDAMem,
    nls: &mut NonlinearSolver,
    tol: f64,
    callLSetup: bool,
) -> i32 {
    match nls {
        NonlinearSolver::Newton(ns) => idaNlsSolveNewtonSensSim(ida_mem, ns, tol, callLSetup),
        /* cannot arise: IDASetNonlinearSolverSensSim rejects non-ROOTFIND types */
        NonlinearSolver::FixedPoint(_) => IDA_NLS_FAIL,
    }
}

/* zero the composite correction: N_VConst(ZERO, ycorSim) */
fn idaNlsZeroYcorSim(ida_mem: &mut IDAMem) {
    N_VConst(ZERO, &mut ida_mem.ida_ee);
    for is in 0..ida_mem.ida_Ns as usize {
        N_VConst(ZERO, &mut ida_mem.ida_eeS[is]);
    }
}

/*
 * SUNNonlinSolSolve_Newton (sunnonlinsol_newton.c), specialized to
 * the IDAS simultaneous-corrector callbacks.  ycor = [ida_ee,
 * ida_eeS], w = [ida_ewt, ida_ewtS] (read directly from IDAMem); the
 * Newton update workspace is the senswrapper ns.deltaS (Ns+1
 * sub-vectors, from SUNNonlinSol_NewtonSens).
 */
fn idaNlsSolveNewtonSensSim(
    ida_mem: &mut IDAMem,
    ns: &mut NewtonSolver,
    tol: f64,
    mut call_lsetup: bool,
) -> i32 {
    /* assume the Jacobian is good */
    let mut jbad = SUNFALSE;

    /* initialize iteration and convergence fail counters for this solve */
    ns.niters = 0;
    ns.nconvfails = 0;

    /* looping point for attempts at solution of the nonlinear system:
       Evaluate the nonlinear residual function (store in delta)
       Setup the linear solver if necessary
       Perform Newton iteration.

       NOTE on break levels (matches the C exactly): a failure of the
       *initial* residual evaluation or of lsetup breaks out of the whole
       setup loop and returns; only failures arising *inside* the Newton
       iteration reach the bad-Jacobian retry below. */
    let mut retval: i32;
    'outer: loop {
        /* initialize current iteration counter for this solve attempt */
        ns.curiter = 0;

        /* compute the nonlinear residual, store in delta */
        retval = {
            let mut deltaS = std::mem::take(&mut ns.deltaS);
            let r = idaNlsResidualSensSim(ida_mem, &mut deltaS);
            ns.deltaS = deltaS;
            r
        };
        if retval != 0 {
            break 'outer;
        }

        /* if indicated, setup the linear system */
        if call_lsetup {
            let mut jcur = ns.jcur;
            retval = idaNlsLSetupSensSim(ida_mem, jbad, &mut jcur);
            ns.jcur = jcur;
            if retval != 0 {
                break 'outer;
            }
        }

        /* looping point for Newton iteration. Break out on any error. */
        loop {
            /* increment nonlinear solver iteration counter */
            ns.niters += 1;

            /* compute the negative of the residual for the linear
               system rhs */
            ns.deltaS.scale_inplace(-ONE);

            /* solve the linear systems to get the Newton update delta */
            retval = idaNlsLSolveSensSim(ida_mem, &mut ns.deltaS);
            if retval != 0 {
                break;
            }

            /* update the Newton iterate: ycorSim = ycorSim + deltaSim */
            ida_mem.ida_ee.linear_sum_with(ONE, ONE, &ns.deltaS.vecs[0]);
            for is in 0..ida_mem.ida_Ns as usize {
                ida_mem.ida_eeS[is].linear_sum_with(ONE, ONE, &ns.deltaS.vecs[is + 1]);
            }

            /* test for convergence */
            retval = {
                let m = ns.curiter;
                let deltaS = std::mem::take(&mut ns.deltaS);
                let r = idaNlsConvTestSensSim(ida_mem, &deltaS, tol, m);
                ns.deltaS = deltaS;
                r
            };

            ns.curiter += 1;

            /* if successful update Jacobian status and return */
            if retval == SUN_SUCCESS {
                ns.jcur = SUNFALSE;
                return SUN_SUCCESS;
            }

            /* check if the iteration should continue; otherwise exit
               Newton loop */
            if retval != SUN_NLS_CONTINUE {
                break;
            }

            /* not yet converged, test for max allowed iterations. */
            if ns.curiter >= ns.maxiters {
                retval = SUN_NLS_CONV_RECVR;
                break;
            }

            /* compute the nonlinear residual, store in delta */
            retval = {
                let mut deltaS = std::mem::take(&mut ns.deltaS);
                let r = idaNlsResidualSensSim(ida_mem, &mut deltaS);
                ns.deltaS = deltaS;
                r
            };
            if retval != 0 {
                break;
            }
        } /* end of Newton iteration loop */

        /* all errors from the Newton iteration go here */

        /* If there is a recoverable convergence failure and the
           Jacobian-related data appears not to be current, increment the
           convergence failure count, reset the initial correction to
           zero, and loop again with a call to lsetup in which jbad is
           TRUE. Otherwise break out and return. */
        if retval > 0 && !ns.jcur && ida_has_lsetup(ida_mem) {
            ns.nconvfails += 1;
            call_lsetup = SUNTRUE;
            jbad = SUNTRUE;
            idaNlsZeroYcorSim(ida_mem);
            continue 'outer;
        } else {
            break 'outer;
        }
    } /* end of setup loop */

    /* increment number of convergence failures */
    ns.nconvfails += 1;

    /* all error returns exit here */
    retval
}

/*===============================================================
  Tests
  ===============================================================*/
#[cfg(test)]
mod tests {
    use super::*;
    use crate::idas_ls::IDASetLinearSolver;
    use crate::idas_nls::IDANls;
    use crate::sundials_context::SUNContext;
    use crate::sunlinsol_dense::SUNLinSol_Dense;
    use crate::sunmatrix_dense::SUNDenseMatrix;
    use crate::sunnonlinsol_newton::SUNNonlinSol_NewtonSens;

    /* decoupled linear DAE: F_i = yp_i + 2*y_i - 4 with one sensitivity
       system FS_i = ypS_i + 2*yS_i - 1 (a user-supplied resS standing in
       for dF/dp = -1).  The Newton matrix for both is
       J = dF/dy + cj*dF/dyp = (2 + cj) I, so the exact corrections are
       e_i  = (4 - yppredict_i - 2*yypredict_i)/(2 + cj) and
       eS_i = (1 - ypSpredict_i - 2*yySpredict_i)/(2 + cj). */
    fn resfn(_t: f64, yy: &NVector, yp: &NVector, rr: &mut NVector, _ud: &mut UserData) -> i32 {
        for i in 0..yy.len() {
            rr.data[i] = yp.data[i] + 2.0 * yy.data[i] - 4.0;
        }
        0
    }

    #[allow(clippy::too_many_arguments)]
    fn resSfn(_Ns: i32, _t: f64, _yy: &NVector, _yp: &NVector, _resval: &NVector,
              yyS: &[NVector], ypS: &[NVector], resvalS: &mut [NVector],
              _user_data: &mut UserData, _tmp1: &mut NVector, _tmp2: &mut NVector,
              _tmp3: &mut NVector) -> i32 {
        for is in 0..yyS.len() {
            for i in 0..yyS[is].len() {
                resvalS[is].data[i] = ypS[is].data[i] + 2.0 * yyS[is].data[i] - 1.0;
            }
        }
        0
    }

    fn make_ida_mem_sens(n: usize, ns: usize) -> IDAMem {
        let mut ida_mem = IDAMem::default();
        ida_mem.ida_res = Some(resfn);
        ida_mem.ida_ewt = NVector::from_slice(&vec![1.0; n]);
        ida_mem.ida_yy = NVector::new(n);
        ida_mem.ida_yp = NVector::new(n);
        ida_mem.ida_ee = NVector::new(n);
        ida_mem.ida_savres = NVector::new(n);
        ida_mem.ida_yypredict = NVector::new(n);
        ida_mem.ida_yppredict = NVector::new(n);
        ida_mem.ida_tempv1 = NVector::new(n);
        ida_mem.ida_tempv2 = NVector::new(n);
        ida_mem.ida_tempv3 = NVector::new(n);
        ida_mem.ida_hh = 1.0e-2;
        /* FSA (simultaneous corrector) state */
        ida_mem.ida_sensi = SUNTRUE;
        ida_mem.ida_ism = IDA_SIMULTANEOUS;
        ida_mem.ida_Ns = ns as i32;
        ida_mem.ida_resS = Some(resSfn);
        ida_mem.ida_resSDQ = SUNFALSE;
        ida_mem.ida_ewtS = (0..ns).map(|_| NVector::from_slice(&vec![1.0; n])).collect();
        ida_mem.ida_eeS = (0..ns).map(|_| NVector::new(n)).collect();
        ida_mem.ida_yyS = (0..ns).map(|_| NVector::new(n)).collect();
        ida_mem.ida_ypS = (0..ns).map(|_| NVector::new(n)).collect();
        ida_mem.ida_yySpredict = (0..ns).map(|_| NVector::new(n)).collect();
        ida_mem.ida_ypSpredict = (0..ns).map(|_| NVector::new(n)).collect();
        ida_mem.ida_tmpS3 = NVector::new(n);
        ida_mem
    }

    /* IDASetNonlinearSolverSensSim stores the solver, clears ownership,
       sets maxiters = MAXIT = 4, marks simMallocDone and wires nls_res
       (idas_nls_sim.c lines 47-200); sensitivity-less memory and a
       non-SIMULTANEOUS ism are rejected. */
    #[test]
    fn idasetnonlinearsolversenssim_defaults_and_rejections() {
        let sunctx = SUNContext::default();
        let mut ida_mem = make_ida_mem_sens(2, 1);

        /* rejection: sensitivities not initialized */
        ida_mem.ida_sensi = SUNFALSE;
        let nls = SUNNonlinSol_NewtonSens(2, &ida_mem.ida_yy, &sunctx);
        assert_eq!(IDASetNonlinearSolverSensSim(&mut ida_mem, nls), IDA_ILL_INPUT);
        assert!(ida_mem.NLSsim.is_none());
        ida_mem.ida_sensi = SUNTRUE;

        /* rejection: staggered corrector selected */
        ida_mem.ida_ism = IDA_STAGGERED;
        let nls = SUNNonlinSol_NewtonSens(2, &ida_mem.ida_yy, &sunctx);
        assert_eq!(IDASetNonlinearSolverSensSim(&mut ida_mem, nls), IDA_ILL_INPUT);
        assert!(ida_mem.NLSsim.is_none());
        ida_mem.ida_ism = IDA_SIMULTANEOUS;

        /* acceptance */
        let nls = SUNNonlinSol_NewtonSens(2, &ida_mem.ida_yy, &sunctx);
        assert_eq!(IDASetNonlinearSolverSensSim(&mut ida_mem, nls), IDA_SUCCESS);
        assert!(!ida_mem.ownNLSsim);
        assert!(ida_mem.simMallocDone);
        assert!(ida_mem.nls_res.is_some());
        match ida_mem.NLSsim.as_ref().unwrap() {
            NonlinearSolver::Newton(ns) => {
                assert_eq!(ns.maxiters, MAXIT);
                assert_eq!(ns.deltaS.vecs.len(), 2);
            }
            _ => unreachable!(),
        }

        /* NLS init succeeds with the solver attached */
        assert_eq!(idaNlsInitSensSim(&mut ida_mem), IDA_SUCCESS);
    }

    /* IDANls (idas_nls.rs) with sensi_sim performs the composite Newton
       solve on [ee, eeS] through idaNlsResidualSensSim /
       idaNlsLSetupSensSim(idaLsSetup) / idaNlsLSolveSensSim(idaLsSolve)
       / idaNlsConvTestSensSim: with linear residuals the first
       correction is exact and the second (zero) correction certifies
       convergence, so niters = 2, nsetups = 1 (first step forces
       lsetup), and yy/yp/yyS/ypS carry the corrected states
       (idas.c IDANls sensi_sim branch). */
    #[test]
    fn idanls_newton_sens_sim_linear_dae() {
        let sunctx = SUNContext::default();
        let mut ida_mem = make_ida_mem_sens(2, 1);

        /* predictor and step data */
        ida_mem.ida_yypredict = NVector::from_slice(&[1.0, 0.5]);
        ida_mem.ida_yppredict = NVector::from_slice(&[0.0, 0.0]);
        ida_mem.ida_cj = 1.0;
        ida_mem.ida_cjlast = 1.0;
        ida_mem.ida_epsNewt = 1.0e-8;
        ida_mem.ida_nst = 0; /* first step: forces callLSetup */

        /* attach dense LS + Newton NLSsim (Ns+1 = 2 sub-vectors) */
        let a = SUNDenseMatrix(2, 2, &sunctx);
        let ls = SUNLinSol_Dense(&ida_mem.ida_tempv1, &a, &sunctx);
        assert_eq!(IDASetLinearSolver(&mut ida_mem, ls, Some(a)), 0);
        {
            /* linit */
            let mut lmem = std::mem::take(&mut ida_mem.ida_lmem);
            if let LsModule::Ls(m) = &mut lmem {
                assert_eq!(crate::idas_ls::idaLsInit(&mut ida_mem, m), 0);
            }
            ida_mem.ida_lmem = lmem;
        }
        let nls = SUNNonlinSol_NewtonSens(2, &ida_mem.ida_yy, &sunctx);
        assert_eq!(IDASetNonlinearSolverSensSim(&mut ida_mem, nls), IDA_SUCCESS);
        assert_eq!(idaNlsInitSensSim(&mut ida_mem), IDA_SUCCESS);

        let retval = IDANls(&mut ida_mem);
        assert_eq!(retval, IDA_SUCCESS);

        /* exact corrections: e = (4 - yppredict - 2*yypredict)/(2 + cj),
           eS = (1 - ypSpredict - 2*yySpredict)/(2 + cj) */
        let e0 = (4.0 - 0.0 - 2.0 * 1.0) / 3.0; /* 2/3 */
        let e1 = (4.0 - 0.0 - 2.0 * 0.5) / 3.0; /* 1   */
        let es = (1.0 - 0.0 - 2.0 * 0.0) / 3.0; /* 1/3 */
        assert!((ida_mem.ida_ee.data[0] - e0).abs() < 1.0e-9);
        assert!((ida_mem.ida_ee.data[1] - e1).abs() < 1.0e-9);
        assert!((ida_mem.ida_eeS[0].data[0] - es).abs() < 1.0e-9);
        assert!((ida_mem.ida_eeS[0].data[1] - es).abs() < 1.0e-9);

        /* yy = yypredict + ee, yp = yppredict + cj*ee; likewise yyS/ypS */
        assert!((ida_mem.ida_yy.data[0] - (1.0 + e0)).abs() < 1.0e-9);
        assert!((ida_mem.ida_yy.data[1] - (0.5 + e1)).abs() < 1.0e-9);
        assert!((ida_mem.ida_yp.data[0] - e0).abs() < 1.0e-9);
        assert!((ida_mem.ida_yp.data[1] - e1).abs() < 1.0e-9);
        assert!((ida_mem.ida_yyS[0].data[0] - es).abs() < 1.0e-9);
        assert!((ida_mem.ida_ypS[0].data[0] - es).abs() < 1.0e-9);

        /* state residual at the solution vanishes and was saved in savres */
        assert!(ida_mem.ida_savres.data[0].abs() < 1.0e-8);
        assert!(ida_mem.ida_savres.data[1].abs() < 1.0e-8);

        /* counters: 2 composite Newton iterations (exact + certifying),
           1 lsetup, 2 state + 2 sens residual evaluations, no failures */
        assert_eq!(ida_mem.ida_nni, 2);
        assert_eq!(ida_mem.ida_nnf, 0);
        assert_eq!(ida_mem.ida_nsetups, 1);
        assert_eq!(ida_mem.ida_nre, 2);
        assert_eq!(ida_mem.ida_nrSe, 2);

        /* lsetup refreshed the convergence-test constants (incl. ssS) */
        assert_eq!(ida_mem.ida_cjold, 1.0);
        assert_eq!(ida_mem.ida_cjratio, 1.0);
        assert_eq!(ida_mem.ida_ssS, 20.0);
        assert!(!ida_mem.ida_forceSetup);
    }
}
