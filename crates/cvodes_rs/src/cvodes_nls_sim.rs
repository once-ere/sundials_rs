/* -----------------------------------------------------------------
 * Translated from src/cvodes/cvodes_nls_sim.c (CVODES 7.7.0),
 * together with the solve drivers of sunnonlinsol_newton.c and
 * sunnonlinsol_fixedpoint.c that the C code reaches through the
 * SUNNonlinearSolver ops table: the SIMULTANEOUS corrector, which
 * iterates on the composite vector [zn[0], znS[0]] (state + all Ns
 * sensitivities at once).
 *
 * In C the composite vectors are N_Vector senswrappers of Ns+1
 * sub-vectors whose slots alias CVodeMem storage:
 *   zn0Sim  = [cv_zn[0],  cv_znS[0][0..Ns]]
 *   ycorSim = [cv_acor,   cv_acorS[0..Ns]]
 *   ewtSim  = [cv_ewt,    cv_ewtS[0..Ns]]
 * Per pinned decision 3 (cvodes_impl.rs) those aliases are NOT
 * stored: this module operates directly on the CVodeMem fields, and
 * only the solver-owned workspaces (NewtonSolver.deltaS /
 * FixedPointSolver.*S — created by SUNNonlinSol_NewtonSens /
 * SUNNonlinSol_FixedPointSens with Ns+1 sub-vectors: sub-vector 0 is
 * the state part, sub-vectors 1..=Ns the sensitivities) remain real
 * senswrappers. Cross-sub-vector reduction semantics of the C
 * senswrapper are reproduced exactly where composite norms are
 * needed: N_VWrmsNorm(wrapper) = MAX of the per-sub-vector WRMS
 * norms (init 0, state sub-vector first), N_VDotProd = SUM (inside
 * anderson_accelerate_sens).
 *
 * Exported solve entry (called by cvodes.c's cvNls, ported in a
 * later part): cvNlsSolveSensSim(cv_mem, nls, tol, callLSetup) is
 * the specialization of
 *   SUNNonlinSolSolve(NLSsim, zn0Sim, ycorSim, ewtSim, tol,
 *                     callLSetup, cv_mem).
 * The caller detaches the solver (cv_mem.NLSsim.take()) for the
 * duration of the call.
 *
 * Forward references to cvodes.c symbols ported in later parts
 * (assumed signatures, pinned for the Part 2/3 briefs):
 *   crate::cvodes::cvSensRhsWrapper(cv_mem: &mut CVodeMem, time: f64,
 *       ycur: &NVector, fcur: &NVector, yScur: &[NVector],
 *       fScur: &mut [NVector], temp1: &mut NVector,
 *       temp2: &mut NVector) -> i32   (increments cv_nfSe itself)
 *   crate::cvodes::cvSensUpdateNorm(cv_mem: &CVodeMem, old_nrm: f64,
 *       xS: &[NVector], wS: &[NVector]) -> f64
 * -----------------------------------------------------------------*/
use crate::cvodes_impl::*;
use crate::cvodes_nls::{cv_has_lsetup, cv_lsetup_dispatch, cv_lsolve_dispatch};
use crate::nvector_serial::*;
use crate::sundials_math::*;
use crate::sundials_nonlinearsolver::*;
use crate::sundials_types::*;
use crate::sunnonlinsol_fixedpoint::FixedPointSolver;
use crate::sunnonlinsol_newton::NewtonSolver;
use sundials_core::sundials_nvector_senswrapper::NVectorSensWrapper;

/* constant macros */
const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/* -----------------------------------------------------------------------------
 * Exported functions
 * ---------------------------------------------------------------------------*/

pub fn CVodeSetNonlinearSolverSensSim(cv_mem: &mut CVodeMem, nls: NonlinearSolver) -> i32 {
    /* check that sensitivities were initialized */
    if !cv_mem.cv_sensi {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSetNonlinearSolverSensSim",
                       file!(), MSGCV_NO_SENSI);
        return CV_ILL_INPUT;
    }

    /* check that simultaneous corrector was selected */
    if cv_mem.cv_ism != CV_SIMULTANEOUS {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSetNonlinearSolverSensSim",
                       file!(), "Sensitivity solution method is not CV_SIMULTANEOUS");
        return CV_ILL_INPUT;
    }

    /* free any existing nonlinear solver: RAII (drop on overwrite) */

    /* set SUNNonlinearSolver pointer; the Sys function is selected by the
       solver type at solve time (cvNlsResidualSensSim /
       cvNlsFPFunctionSensSim), and the convergence test is inlined in the
       solve drivers below */
    cv_mem.NLSsim = Some(nls);

    /* Set NLS ownership flag. If this function was called to attach the
       default NLS, CVODES will set the flag to SUNTRUE after this function
       returns. */
    cv_mem.ownNLSsim = SUNFALSE;

    /* set max allowed nonlinear iterations */
    if let Some(s) = cv_mem.NLSsim.as_mut() {
        let retval = s.set_max_iters(NLS_MAXCOR);
        if retval != CV_SUCCESS {
            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSetNonlinearSolverSensSim",
                           file!(),
                           "Setting maximum number of nonlinear iterations failed");
            return CV_ILL_INPUT;
        }
    }

    /* create vector wrappers if necessary: per pinned decision 3 the
       zn0Sim/ycorSim/ewtSim senswrapper aliases are not stored — this
       module reads cv_zn[0]/cv_znS[0], cv_acor/cv_acorS and cv_ewt/cv_ewtS
       directly — so only the allocation flag survives */
    if !cv_mem.simMallocDone {
        cv_mem.simMallocDone = SUNTRUE;
    }

    /* Reset the acnrmcur flag to SUNFALSE */
    cv_mem.cv_acnrmcur = SUNFALSE;

    /* Set the nonlinear system RHS function */
    if cv_mem.cv_f.is_none() {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSetNonlinearSolverSensSim",
                       file!(), "The ODE RHS function is NULL");
        return CV_ILL_INPUT;
    }
    cv_mem.nls_f = cv_mem.cv_f;

    CV_SUCCESS
}

/*---------------------------------------------------------------
  CVodeGetNonlinearSystemDataSens:

  This routine provides access to the relevant data needed to
  compute the nonlinear system function (out-pointers become a
  returned tuple: (tcur, gamma, rl1); the vector arrays
  ySpred = cv_znS[0], ySn = cv_yS, znS1 = cv_znS[1] and user_data
  remain accessible as CVodeMem fields — donor adaptation).
  ---------------------------------------------------------------*/
pub fn CVodeGetNonlinearSystemDataSens(cv_mem: &CVodeMem) -> (f64, f64, f64) {
    (cv_mem.cv_tn, cv_mem.cv_gamma, cv_mem.cv_rl1)
}

/* -----------------------------------------------------------------------------
 * Private functions
 * ---------------------------------------------------------------------------*/

pub fn cvNlsInitSensSim(cv_mem: &mut CVodeMem) -> i32 {
    /* In C this wires cvNlsLSetupSensSim/cvNlsLSolveSensSim into the NLS
       depending on whether cv_lsetup/cv_lsolve exist; here the dispatch is
       dynamic. A Newton solver without an attached linear solver cannot
       work (SUNNonlinSolSolve_Newton requires an LSolve function). */
    if let Some(nls) = cv_mem.NLSsim.as_ref() {
        if nls.nls_type() == SUNNONLINEARSOLVER_ROOTFIND && cv_mem.cv_lmem.is_none() {
            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "cvNlsInitSensSim", file!(),
                           MSGCV_LSOLVE_NULL);
            return CV_NLS_INIT_FAIL;
        }
    }

    /* initialize nonlinear solver */
    let retval = match cv_mem.NLSsim.as_mut() {
        Some(nls) => nls.initialize(),
        None => -1,
    };

    if retval != CV_SUCCESS {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "cvNlsInitSensSim", file!(),
                       MSGCV_NLS_INIT_FAIL);
        return CV_NLS_INIT_FAIL;
    }

    CV_SUCCESS
}

/* cvNlsLSetupSensSim (cvodes_nls_sim.c): wrapper around the lsetup
   dispatch */
fn cvNlsLSetupSensSim(cv_mem: &mut CVodeMem, jbad: bool, jcur: &mut bool) -> i32 {
    /* if the nonlinear solver marked the Jacobian as bad update convfail */
    if jbad {
        cv_mem.convfail = CV_FAIL_BAD_J;
    }

    /* setup the linear solver */
    let convfail = cv_mem.convfail;
    let mut cv_jcur = cv_mem.cv_jcur;
    let retval = cv_lsetup_dispatch(cv_mem, convfail, &mut cv_jcur);
    cv_mem.cv_jcur = cv_jcur;
    cv_mem.cv_nsetups += 1;

    /* update Jacobian status */
    *jcur = cv_mem.cv_jcur;

    cv_mem.cv_forceSetup = SUNFALSE;
    cv_mem.cv_gamrat = ONE;
    cv_mem.cv_gammap = cv_mem.cv_gamma;
    cv_mem.cv_crate = ONE;
    cv_mem.cv_crateS = ONE;
    cv_mem.cv_nstlp = cv_mem.cv_nst;

    if retval < 0 {
        return CV_LSETUP_FAIL;
    }
    if retval > 0 {
        return SUN_NLS_CONV_RECVR;
    }

    CV_SUCCESS
}

/* cvNlsLSolveSensSim (cvodes_nls_sim.c): solves the state linear system
   (weight = cv_ewt) followed by the Ns sensitivity linear systems
   (weight = cv_ewtS[is]); deltaSim sub-vector 0 = state delta,
   sub-vectors 1..=Ns = sensitivity deltas */
fn cvNlsLSolveSensSim(cv_mem: &mut CVodeMem, deltaSim: &mut NVectorSensWrapper) -> i32 {
    /* solve the state linear system */
    let retval = cv_lsolve_dispatch(cv_mem, &mut deltaSim.vecs[0], None);

    if retval < 0 {
        return CV_LSOLVE_FAIL;
    }
    if retval > 0 {
        return SUN_NLS_CONV_RECVR;
    }

    /* solve the sensitivity linear systems */
    for is in 0..cv_mem.cv_Ns as usize {
        let retval = cv_lsolve_dispatch(cv_mem, &mut deltaSim.vecs[is + 1], Some(is));

        if retval < 0 {
            return CV_LSOLVE_FAIL;
        }
        if retval > 0 {
            return SUN_NLS_CONV_RECVR;
        }
    }

    CV_SUCCESS
}

/* cvNlsConvTestSensSim (cvodes_nls_sim.c); ycorSim = [cv_acor, cv_acorS],
   ewtSim = [cv_ewt, cv_ewtS] (read directly from cv_mem), m = current
   iteration count */
fn cvNlsConvTestSensSim(
    cv_mem: &mut CVodeMem,
    deltaSim: &NVectorSensWrapper,
    tol: f64,
    m: i32,
) -> i32 {
    /* compute the norm of the state and sensitivity corrections */
    let del = N_VWrmsNorm(&deltaSim.vecs[0], &cv_mem.cv_ewt);
    let delS = crate::cvodes::cvSensUpdateNorm(cv_mem, del, &deltaSim.vecs[1..], &cv_mem.cv_ewtS);

    /* norm used in error test */
    let Del = delS;

    /* Test for convergence. If m > 0, an estimate of the convergence
       rate constant is stored in crate, and used in the test.

       Recall that, even when errconS=SUNFALSE, all variables are used in the
       convergence test. Hence, we use Del (and not del). However, acnrm is used
       in the error test and thus it has different forms depending on errconS
       (and this explains why we have to carry around del and delS).
    */
    if m > 0 {
        cv_mem.cv_crate = SUNMAX(CRDOWN * cv_mem.cv_crate, Del / cv_mem.cv_delp);
    }
    let dcon = Del * SUNMIN(ONE, cv_mem.cv_crate) / tol;

    /* check if nonlinear system was solved successfully */
    if dcon <= ONE {
        if m == 0 {
            cv_mem.cv_acnrm = if cv_mem.cv_errconS { delS } else { del };
        } else {
            cv_mem.cv_acnrm = if cv_mem.cv_errconS {
                /* N_VWrmsNorm(ycorSim, ewtSim): senswrapper WRMS norm = MAX
                   of the per-sub-vector WRMS norms (init 0), state
                   sub-vector first */
                let mut nrm = ZERO;
                let tmp = N_VWrmsNorm(&cv_mem.cv_acor, &cv_mem.cv_ewt);
                if tmp > nrm {
                    nrm = tmp;
                }
                for is in 0..cv_mem.cv_Ns as usize {
                    let tmp = N_VWrmsNorm(&cv_mem.cv_acorS[is], &cv_mem.cv_ewtS[is]);
                    if tmp > nrm {
                        nrm = tmp;
                    }
                }
                nrm
            } else {
                N_VWrmsNorm(&cv_mem.cv_acor, &cv_mem.cv_ewt)
            };
        }
        cv_mem.cv_acnrmcur = SUNTRUE;
        return CV_SUCCESS;
    }

    /* check if the iteration seems to be diverging */
    if m >= 1 && Del > RDIV * cv_mem.cv_delp {
        return SUN_NLS_CONV_RECVR;
    }

    /* Save norm of correction and loop again */
    cv_mem.cv_delp = Del;

    /* Not yet converged */
    SUN_NLS_CONTINUE
}

/* cvNlsResidualSensSim (cvodes_nls_sim.c): evaluates the state residual
   (into resSim sub-vector 0) and the Ns sensitivity residuals (into
   sub-vectors 1..=Ns); ycorSim = [cv_acor, cv_acorS]. */
fn cvNlsResidualSensSim(cv_mem: &mut CVodeMem, resSim: &mut NVectorSensWrapper) -> i32 {
    /* update the state based on the current correction */
    {
        let CVodeMem { cv_zn, cv_acor, cv_y, .. } = cv_mem;
        N_VLinearSum(ONE, &cv_zn[0], ONE, cv_acor, cv_y);
    }

    /* evaluate the rhs function */
    let f = cv_mem.nls_f.unwrap();
    let retval = f(
        cv_mem.cv_tn,
        &cv_mem.cv_y,
        &mut cv_mem.cv_ftemp,
        &mut cv_mem.cv_user_data,
    );
    cv_mem.cv_nfe += 1;
    if retval < 0 {
        return CV_RHSFUNC_FAIL;
    }
    if retval > 0 {
        return RHSFUNC_RECVR;
    }

    /* compute the residual */
    {
        let CVodeMem { cv_zn, cv_acor, cv_ftemp, cv_rl1, cv_gamma, .. } = cv_mem;
        let res = &mut resSim.vecs[0];
        N_VLinearSum(*cv_rl1, &cv_zn[1], ONE, cv_acor, res);
        /* res = -gamma*ftemp + res */
        res.linear_sum_with(ONE, -*cv_gamma, cv_ftemp);
    }

    /* update sensitivities based on the current correction
       (N_VLinearSumVectorArray expanded to a per-vector loop) */
    {
        let CVodeMem { cv_znS, cv_acorS, cv_yS, cv_Ns, .. } = cv_mem;
        for is in 0..*cv_Ns as usize {
            N_VLinearSum(ONE, &cv_znS[0][is], ONE, &cv_acorS[is], &mut cv_yS[is]);
        }
    }

    /* evaluate the sensitivity rhs function (cvSensRhsWrapper increments
       cv_nfSe itself; the argument vectors are CVodeMem fields, taken out
       for the duration of the call — donor take()/restore pattern) */
    let cv_y = std::mem::take(&mut cv_mem.cv_y);
    let cv_ftemp = std::mem::take(&mut cv_mem.cv_ftemp);
    let cv_yS = std::mem::take(&mut cv_mem.cv_yS);
    let mut cv_ftempS = std::mem::take(&mut cv_mem.cv_ftempS);
    let mut vtemp1 = std::mem::take(&mut cv_mem.cv_vtemp1);
    let mut vtemp2 = std::mem::take(&mut cv_mem.cv_vtemp2);
    let tn = cv_mem.cv_tn;
    let retval = crate::cvodes::cvSensRhsWrapper(
        cv_mem,
        tn,
        &cv_y,
        &cv_ftemp,
        &cv_yS,
        &mut cv_ftempS,
        &mut vtemp1,
        &mut vtemp2,
    );
    cv_mem.cv_y = cv_y;
    cv_mem.cv_ftemp = cv_ftemp;
    cv_mem.cv_yS = cv_yS;
    cv_mem.cv_ftempS = cv_ftempS;
    cv_mem.cv_vtemp1 = vtemp1;
    cv_mem.cv_vtemp2 = vtemp2;

    if retval < 0 {
        return CV_SRHSFUNC_FAIL;
    }
    if retval > 0 {
        return SRHSFUNC_RECVR;
    }

    /* compute the sensitivity residual: resS = rl1*znS[1] + ycorS
       - gamma*ftempS. The C fused call
       N_VLinearCombinationVectorArray(Ns, 3, {rl1, 1, -gamma},
       {znS[1], ycorS, ftempS}, resS) expands (serial fallback:
       N_VScale + N_VLinearSum accumulation) to the bit-identical
       two-step form below. */
    {
        let CVodeMem { cv_znS, cv_acorS, cv_ftempS, cv_rl1, cv_gamma, cv_Ns, .. } = cv_mem;
        for is in 0..*cv_Ns as usize {
            let resS = &mut resSim.vecs[is + 1];
            N_VLinearSum(*cv_rl1, &cv_znS[1][is], ONE, &cv_acorS[is], resS);
            /* resS = -gamma*ftempS + resS */
            resS.linear_sum_with(ONE, -*cv_gamma, &cv_ftempS[is]);
        }
    }

    CV_SUCCESS
}

/* cvNlsFPFunctionSensSim (cvodes_nls_sim.c): fixed-point function for the
   composite system; resSim sub-vector 0 = state, 1..=Ns = sensitivities.
   NOTE (as in C): the state fixed-point value res — already transformed to
   rl1*(h*f - zn[1]) — is what gets passed as `fcur` to cvSensRhsWrapper. */
fn cvNlsFPFunctionSensSim(cv_mem: &mut CVodeMem, resSim: &mut NVectorSensWrapper) -> i32 {
    /* update the state based on the current correction */
    {
        let CVodeMem { cv_zn, cv_acor, cv_y, .. } = cv_mem;
        N_VLinearSum(ONE, &cv_zn[0], ONE, cv_acor, cv_y);
    }

    /* evaluate the rhs function */
    let f = cv_mem.nls_f.unwrap();
    let retval = f(
        cv_mem.cv_tn,
        &cv_mem.cv_y,
        &mut resSim.vecs[0],
        &mut cv_mem.cv_user_data,
    );
    cv_mem.cv_nfe += 1;
    if retval < 0 {
        return CV_RHSFUNC_FAIL;
    }
    if retval > 0 {
        return RHSFUNC_RECVR;
    }

    /* evaluate fixed point function */
    {
        let CVodeMem { cv_zn, cv_h, cv_rl1, .. } = cv_mem;
        let res = &mut resSim.vecs[0];
        /* res = h*res - zn[1] */
        res.linear_sum_with(*cv_h, -ONE, &cv_zn[1]);
        res.scale_inplace(*cv_rl1);
    }

    /* update the sensitivities based on the current correction */
    {
        let CVodeMem { cv_znS, cv_acorS, cv_yS, cv_Ns, .. } = cv_mem;
        for is in 0..*cv_Ns as usize {
            N_VLinearSum(ONE, &cv_znS[0][is], ONE, &cv_acorS[is], &mut cv_yS[is]);
        }
    }

    /* evaluate the sensitivity rhs function (fcur = the state fixed-point
       value res, fScur = the sensitivity sub-vectors of resSim) */
    let cv_y = std::mem::take(&mut cv_mem.cv_y);
    let cv_yS = std::mem::take(&mut cv_mem.cv_yS);
    let mut vtemp1 = std::mem::take(&mut cv_mem.cv_vtemp1);
    let mut vtemp2 = std::mem::take(&mut cv_mem.cv_vtemp2);
    let tn = cv_mem.cv_tn;
    let retval = {
        let (res, resS) = resSim.vecs.split_at_mut(1);
        crate::cvodes::cvSensRhsWrapper(
            cv_mem,
            tn,
            &cv_y,
            &res[0],
            &cv_yS,
            resS,
            &mut vtemp1,
            &mut vtemp2,
        )
    };
    cv_mem.cv_y = cv_y;
    cv_mem.cv_yS = cv_yS;
    cv_mem.cv_vtemp1 = vtemp1;
    cv_mem.cv_vtemp2 = vtemp2;

    if retval < 0 {
        return CV_SRHSFUNC_FAIL;
    }
    if retval > 0 {
        return SRHSFUNC_RECVR;
    }

    /* evaluate sensitivity fixed point function */
    {
        let CVodeMem { cv_znS, cv_h, cv_rl1, cv_Ns, .. } = cv_mem;
        for is in 0..*cv_Ns as usize {
            let resS = &mut resSim.vecs[is + 1];
            /* resS = h*resS - znS[1] */
            resS.linear_sum_with(*cv_h, -ONE, &cv_znS[1][is]);
            resS.scale_inplace(*cv_rl1);
        }
    }

    CV_SUCCESS
}

/*
 * cvNlsSolveSensSim — SUNNonlinSolSolve(cv_mem->NLSsim, cv_mem->zn0Sim,
 * cv_mem->ycorSim, cv_mem->ewtSim, tol, callLSetup, cv_mem) as invoked
 * by cvodes.c's cvNls (ported in a later part). The caller detaches the
 * solver from cv_mem (NLSsim.take()) and reattaches it afterwards,
 * reading the niters/nconvfails counters for the nniS/nnfS updates.
 */
pub fn cvNlsSolveSensSim(
    cv_mem: &mut CVodeMem,
    nls: &mut NonlinearSolver,
    tol: f64,
    callLSetup: bool,
) -> i32 {
    match nls {
        NonlinearSolver::Newton(ns) => cvNlsSolveNewtonSensSim(cv_mem, ns, tol, callLSetup),
        NonlinearSolver::FixedPoint(fps) => cvNlsSolveFixedPointSensSim(cv_mem, fps, tol),
    }
}

/* zero the composite correction: N_VConst(ZERO, ycorSim) */
fn cvNlsZeroYcorSim(cv_mem: &mut CVodeMem) {
    N_VConst(ZERO, &mut cv_mem.cv_acor);
    for is in 0..cv_mem.cv_Ns as usize {
        N_VConst(ZERO, &mut cv_mem.cv_acorS[is]);
    }
}

/*
 * SUNNonlinSolSolve_Newton (sunnonlinsol_newton.c), specialized to the
 * CVODES simultaneous-corrector callbacks. ycor = [cv_acor, cv_acorS],
 * w = [cv_ewt, cv_ewtS]; the Newton update workspace is the senswrapper
 * ns.deltaS (Ns+1 sub-vectors, from SUNNonlinSol_NewtonSens).
 */
fn cvNlsSolveNewtonSensSim(
    cv_mem: &mut CVodeMem,
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
       Perform Newton iteration (see the break-level NOTE in
       cvodes_nls.rs::cvNlsSolveNewton — the structure is identical). */
    let mut retval: i32;
    'outer: loop {
        /* initialize current iteration counter for this solve attempt */
        ns.curiter = 0;

        /* compute the nonlinear residual, store in delta */
        retval = {
            let mut deltaS = std::mem::take(&mut ns.deltaS);
            let r = cvNlsResidualSensSim(cv_mem, &mut deltaS);
            ns.deltaS = deltaS;
            r
        };
        if retval != 0 {
            break 'outer;
        }

        /* if indicated, setup the linear system */
        if call_lsetup {
            let mut jcur = ns.jcur;
            retval = cvNlsLSetupSensSim(cv_mem, jbad, &mut jcur);
            ns.jcur = jcur;
            if retval != 0 {
                break 'outer;
            }
        }

        /* looping point for Newton iteration. Break out on any error. */
        loop {
            /* increment nonlinear solver iteration counter */
            ns.niters += 1;
            cv_mem.cv_nls_curiter = ns.curiter;

            /* compute the negative of the residual for the linear
               system rhs */
            ns.deltaS.scale_inplace(-ONE);

            /* solve the linear systems to get the Newton update delta */
            retval = cvNlsLSolveSensSim(cv_mem, &mut ns.deltaS);
            if retval != 0 {
                break;
            }

            /* update the Newton iterate: ycorSim = ycorSim + deltaSim */
            cv_mem.cv_acor.linear_sum_with(ONE, ONE, &ns.deltaS.vecs[0]);
            for is in 0..cv_mem.cv_Ns as usize {
                cv_mem.cv_acorS[is].linear_sum_with(ONE, ONE, &ns.deltaS.vecs[is + 1]);
            }

            /* test for convergence */
            retval = cvNlsConvTestSensSim(cv_mem, &ns.deltaS, tol, ns.curiter);

            ns.curiter += 1;

            /* if successful update Jacobian status and return */
            if retval == CV_SUCCESS {
                ns.jcur = SUNFALSE;
                return 0;
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
                let r = cvNlsResidualSensSim(cv_mem, &mut deltaS);
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
        if retval > 0 && !ns.jcur && cv_has_lsetup(cv_mem) {
            ns.nconvfails += 1;
            call_lsetup = SUNTRUE;
            jbad = SUNTRUE;
            cvNlsZeroYcorSim(cv_mem);
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

/*
 * SUNNonlinSolSolve_FixedPoint (sunnonlinsol_fixedpoint.c), specialized
 * to the CVODES simultaneous-corrector callbacks. ycor = [cv_acor,
 * cv_acorS]; the workspaces are the senswrappers fps.yprevS/gyS/deltaS
 * (Ns+1 sub-vectors, from SUNNonlinSol_FixedPointSens).
 */
fn cvNlsSolveFixedPointSensSim(
    cv_mem: &mut CVodeMem,
    fps: &mut FixedPointSolver,
    tol: f64,
) -> i32 {
    /* initialize iteration and convergence fail counters for this solve */
    fps.niters = 0;
    fps.nconvfails = 0;

    /* Looping point for attempts at solution of the nonlinear system:
       Evaluate fixed-point function (store in gy).
       Perform the accelerated fixed-point iteration.
       Perform stopping tests. */
    fps.curiter = 0;
    while fps.curiter < fps.maxiters {
        cv_mem.cv_nls_curiter = fps.curiter;
        /* update previous solution guess: yprev = ycorSim */
        fps.yprevS.vecs[0].data.copy_from_slice(&cv_mem.cv_acor.data);
        for is in 0..cv_mem.cv_Ns as usize {
            fps.yprevS.vecs[is + 1]
                .data
                .copy_from_slice(&cv_mem.cv_acorS[is].data);
        }

        /* compute fixed-point iteration function, store in gy */
        let retval = cvNlsFPFunctionSensSim(cv_mem, &mut fps.gyS);
        if retval != 0 {
            return retval;
        }

        /* perform fixed point update, based on choice of acceleration or not */
        if fps.m == 0 {
            /* basic fixed-point solver: ycorSim = gy */
            cv_mem.cv_acor.data.copy_from_slice(&fps.gyS.vecs[0].data);
            for is in 0..cv_mem.cv_Ns as usize {
                cv_mem.cv_acorS[is]
                    .data
                    .copy_from_slice(&fps.gyS.vecs[is + 1].data);
            }
        } else {
            /* Anderson-accelerated solver: assemble the composite iterate
               x = [cv_acor, cv_acorS] by moving the CVodeMem vectors into
               a senswrapper (no data copies), run the composite kernel,
               and move them back */
            let mut x = NVectorSensWrapper {
                vecs: Vec::with_capacity(cv_mem.cv_Ns as usize + 1),
                own_vecs: false,
            };
            x.vecs.push(std::mem::take(&mut cv_mem.cv_acor));
            x.vecs.append(&mut cv_mem.cv_acorS);
            let iter = fps.curiter;
            fps.anderson_accelerate_sens(&mut x, iter);
            let mut it = x.vecs.into_iter();
            cv_mem.cv_acor = it.next().unwrap();
            cv_mem.cv_acorS = it.collect();
        }

        /* increment nonlinear solver iteration counter */
        fps.niters += 1;

        /* compute change in solution: delta = ycorSim - yprev */
        {
            let FixedPointSolver { yprevS, deltaS, .. } = fps;
            N_VLinearSum(ONE, &cv_mem.cv_acor, -ONE, &yprevS.vecs[0], &mut deltaS.vecs[0]);
            for is in 0..cv_mem.cv_Ns as usize {
                N_VLinearSum(
                    ONE,
                    &cv_mem.cv_acorS[is],
                    -ONE,
                    &yprevS.vecs[is + 1],
                    &mut deltaS.vecs[is + 1],
                );
            }
        }

        /* test for convergence */
        let retval = {
            let m = fps.curiter;
            cvNlsConvTestSensSim(cv_mem, &fps.deltaS, tol, m)
        };

        /* return if successful */
        if retval == 0 {
            return 0;
        }

        /* check if the iterations should continue; otherwise increment the
           convergence failure count and return error flag */
        if retval != SUN_NLS_CONTINUE {
            fps.nconvfails += 1;
            return retval;
        }

        fps.curiter += 1;
    }

    /* if we've reached this point, then we exhausted the iteration limit;
       increment the convergence failure count and return */
    fps.nconvfails += 1;
    SUN_NLS_CONV_RECVR
}
