/* -----------------------------------------------------------------
 * Translated from src/cvodes/cvodes_nls_stg.c (CVODES 7.7.0),
 * together with the solve drivers of sunnonlinsol_newton.c and
 * sunnonlinsol_fixedpoint.c that the C code reaches through the
 * SUNNonlinearSolver ops table: the STAGGERED corrector, which
 * iterates on all Ns sensitivities at once AFTER the state
 * nonlinear system has converged.
 *
 * In C the composite vectors are N_Vector senswrappers of Ns
 * sub-vectors whose slots alias CVodeMem storage:
 *   zn0Stg  = cv_znS[0][0..Ns]
 *   ycorStg = cv_acorS[0..Ns]
 *   ewtStg  = cv_ewtS[0..Ns]
 * Per pinned decision 3 (cvodes_impl.rs) those aliases are NOT
 * stored: this module operates directly on the CVodeMem fields, and
 * only the solver-owned workspaces (NewtonSolver.deltaS /
 * FixedPointSolver.*S — created by SUNNonlinSol_NewtonSens /
 * SUNNonlinSol_FixedPointSens with Ns sub-vectors: sub-vector is =
 * sensitivity is) remain real senswrappers. Composite norms use the
 * C senswrapper reduction semantics: WRMS norm = MAX of the
 * per-sub-vector WRMS norms (cvSensNorm), DotProd = SUM (inside
 * anderson_accelerate_sens).
 *
 * Exported solve entry (called by cvodes.c's cvStgrNls, ported in a
 * later part): cvNlsSolveSensStg(cv_mem, nls, tol, callLSetup) is
 * the specialization of
 *   SUNNonlinSolSolve(NLSstg, zn0Stg, ycorStg, ewtStg, tol,
 *                     callLSetup, cv_mem).
 * The caller detaches the solver (cv_mem.NLSstg.take()) for the
 * duration of the call.
 *
 * Forward references to cvodes.c symbols ported in later parts
 * (assumed signatures, pinned for the Part 2/3 briefs):
 *   crate::cvodes::cvSensRhsWrapper(cv_mem: &mut CVodeMem, time: f64,
 *       ycur: &NVector, fcur: &NVector, yScur: &[NVector],
 *       fScur: &mut [NVector], temp1: &mut NVector,
 *       temp2: &mut NVector) -> i32   (increments cv_nfSe itself)
 *   crate::cvodes::cvSensNorm(cv_mem: &CVodeMem, xS: &[NVector],
 *       wS: &[NVector]) -> f64
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

pub fn CVodeSetNonlinearSolverSensStg(cv_mem: &mut CVodeMem, nls: NonlinearSolver) -> i32 {
    /* check that sensitivities were initialized */
    if !cv_mem.cv_sensi {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSetNonlinearSolverSensStg",
                       file!(), MSGCV_NO_SENSI);
        return CV_ILL_INPUT;
    }

    /* check that staggered corrector was selected */
    if cv_mem.cv_ism != CV_STAGGERED {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSetNonlinearSolverSensStg",
                       file!(), "Sensitivity solution method is not CV_STAGGERED");
        return CV_ILL_INPUT;
    }

    /* free any existing nonlinear solver: RAII (drop on overwrite) */

    /* set SUNNonlinearSolver pointer; the Sys function is selected by the
       solver type at solve time (cvNlsResidualSensStg /
       cvNlsFPFunctionSensStg), and the convergence test is inlined in the
       solve drivers below */
    cv_mem.NLSstg = Some(nls);

    /* Set NLS ownership flag. If this function was called to attach the
       default NLS, CVODES will set the flag to SUNTRUE after this function
       returns. */
    cv_mem.ownNLSstg = SUNFALSE;

    /* set max allowed nonlinear iterations */
    if let Some(s) = cv_mem.NLSstg.as_mut() {
        let retval = s.set_max_iters(NLS_MAXCOR);
        if retval != CV_SUCCESS {
            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSetNonlinearSolverSensStg",
                           file!(),
                           "Setting maximum number of nonlinear iterations failed");
            return CV_ILL_INPUT;
        }
    }

    /* create vector wrappers if necessary: per pinned decision 3 the
       zn0Stg/ycorStg/ewtStg senswrapper aliases are not stored — this
       module reads cv_znS[0], cv_acorS and cv_ewtS directly — so only the
       allocation flag survives */
    if !cv_mem.stgMallocDone {
        cv_mem.stgMallocDone = SUNTRUE;
    }

    /* Reset the acnrmScur flag to SUNFALSE */
    cv_mem.cv_acnrmScur = SUNFALSE;

    CV_SUCCESS
}

/* -----------------------------------------------------------------------------
 * Private functions
 * ---------------------------------------------------------------------------*/

pub fn cvNlsInitSensStg(cv_mem: &mut CVodeMem) -> i32 {
    /* In C this wires cvNlsLSetupSensStg/cvNlsLSolveSensStg into the NLS
       depending on whether cv_lsetup/cv_lsolve exist; here the dispatch is
       dynamic. A Newton solver without an attached linear solver cannot
       work (SUNNonlinSolSolve_Newton requires an LSolve function). */
    if let Some(nls) = cv_mem.NLSstg.as_ref() {
        if nls.nls_type() == SUNNONLINEARSOLVER_ROOTFIND && cv_mem.cv_lmem.is_none() {
            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "cvNlsInitSensStg", file!(),
                           MSGCV_LSOLVE_NULL);
            return CV_NLS_INIT_FAIL;
        }
    }

    /* initialize nonlinear solver */
    let retval = match cv_mem.NLSstg.as_mut() {
        Some(nls) => nls.initialize(),
        None => -1,
    };

    if retval != CV_SUCCESS {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "cvNlsInitSensStg", file!(),
                       MSGCV_NLS_INIT_FAIL);
        return CV_NLS_INIT_FAIL;
    }

    CV_SUCCESS
}

/* cvNlsLSetupSensStg (cvodes_nls_stg.c): wrapper around the lsetup
   dispatch (unlike the state/simultaneous variants this also counts the
   setup against the sensitivities and does not touch cv_forceSetup) */
fn cvNlsLSetupSensStg(cv_mem: &mut CVodeMem, jbad: bool, jcur: &mut bool) -> i32 {
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
    cv_mem.cv_nsetupsS += 1;

    /* update Jacobian status */
    *jcur = cv_mem.cv_jcur;

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

/* cvNlsLSolveSensStg (cvodes_nls_stg.c): solves the Ns sensitivity linear
   systems (weight = cv_ewtS[is]); deltaStg sub-vector is = sensitivity is */
fn cvNlsLSolveSensStg(cv_mem: &mut CVodeMem, deltaStg: &mut NVectorSensWrapper) -> i32 {
    /* solve the sensitivity linear systems */
    for is in 0..cv_mem.cv_Ns as usize {
        let retval = cv_lsolve_dispatch(cv_mem, &mut deltaStg.vecs[is], Some(is));

        if retval < 0 {
            return CV_LSOLVE_FAIL;
        }
        if retval > 0 {
            return SUN_NLS_CONV_RECVR;
        }
    }

    CV_SUCCESS
}

/* cvNlsConvTestSensStg (cvodes_nls_stg.c); ycorStg = cv_acorS,
   ewtStg = cv_ewtS (read directly from cv_mem), m = current iteration
   count */
fn cvNlsConvTestSensStg(
    cv_mem: &mut CVodeMem,
    deltaStg: &NVectorSensWrapper,
    tol: f64,
    m: i32,
) -> i32 {
    /* compute the norm of the state and sensitivity corrections */
    let Del = crate::cvodes::cvSensNorm(cv_mem, &deltaStg.vecs, &cv_mem.cv_ewtS);

    /* Test for convergence. If m > 0, an estimate of the convergence
       rate constant is stored in crate, and used in the test.

       Recall that, even when errconS=SUNFALSE, all variables are used in the
       convergence test. Hence, we use Del (and not del). However, acnrm is used
       in the error test and thus it has different forms depending on errconS
       (and this explains why we have to carry around del and delS).
    */
    if m > 0 {
        cv_mem.cv_crateS = SUNMAX(CRDOWN * cv_mem.cv_crateS, Del / cv_mem.cv_delp);
    }
    let dcon = Del * SUNMIN(ONE, cv_mem.cv_crateS) / tol;

    /* check if nonlinear system was solved successfully */
    if dcon <= ONE {
        if cv_mem.cv_errconS {
            cv_mem.cv_acnrmS = if m == 0 {
                Del
            } else {
                crate::cvodes::cvSensNorm(cv_mem, &cv_mem.cv_acorS, &cv_mem.cv_ewtS)
            };
            cv_mem.cv_acnrmScur = SUNTRUE;
        }
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

/* cvNlsResidualSensStg (cvodes_nls_stg.c): evaluates the Ns sensitivity
   residuals into resStg; ycorStg = cv_acorS. The state values cv_y and
   cv_ftemp hold the already-converged state solution. */
fn cvNlsResidualSensStg(cv_mem: &mut CVodeMem, resStg: &mut NVectorSensWrapper) -> i32 {
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
            let resS = &mut resStg.vecs[is];
            N_VLinearSum(*cv_rl1, &cv_znS[1][is], ONE, &cv_acorS[is], resS);
            /* resS = -gamma*ftempS + resS */
            resS.linear_sum_with(ONE, -*cv_gamma, &cv_ftempS[is]);
        }
    }

    CV_SUCCESS
}

/* cvNlsFPFunctionSensStg (cvodes_nls_stg.c): fixed-point function for the
   staggered sensitivity system; resStg sub-vector is = sensitivity is. */
fn cvNlsFPFunctionSensStg(cv_mem: &mut CVodeMem, resStg: &mut NVectorSensWrapper) -> i32 {
    /* update the sensitivities based on the current correction */
    {
        let CVodeMem { cv_znS, cv_acorS, cv_yS, cv_Ns, .. } = cv_mem;
        for is in 0..*cv_Ns as usize {
            N_VLinearSum(ONE, &cv_znS[0][is], ONE, &cv_acorS[is], &mut cv_yS[is]);
        }
    }

    /* evaluate the sensitivity rhs function (fScur = resStg) */
    let cv_y = std::mem::take(&mut cv_mem.cv_y);
    let cv_ftemp = std::mem::take(&mut cv_mem.cv_ftemp);
    let cv_yS = std::mem::take(&mut cv_mem.cv_yS);
    let mut vtemp1 = std::mem::take(&mut cv_mem.cv_vtemp1);
    let mut vtemp2 = std::mem::take(&mut cv_mem.cv_vtemp2);
    let tn = cv_mem.cv_tn;
    let retval = crate::cvodes::cvSensRhsWrapper(
        cv_mem,
        tn,
        &cv_y,
        &cv_ftemp,
        &cv_yS,
        &mut resStg.vecs,
        &mut vtemp1,
        &mut vtemp2,
    );
    cv_mem.cv_y = cv_y;
    cv_mem.cv_ftemp = cv_ftemp;
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
            let resS = &mut resStg.vecs[is];
            /* resS = h*resS - znS[1] */
            resS.linear_sum_with(*cv_h, -ONE, &cv_znS[1][is]);
            resS.scale_inplace(*cv_rl1);
        }
    }

    CV_SUCCESS
}

/*
 * cvNlsSolveSensStg — SUNNonlinSolSolve(cv_mem->NLSstg, cv_mem->zn0Stg,
 * cv_mem->ycorStg, cv_mem->ewtStg, tol, callLSetup, cv_mem) as invoked
 * by cvodes.c's cvStgrNls (ported in a later part). The caller detaches
 * the solver from cv_mem (NLSstg.take()) and reattaches it afterwards,
 * reading the niters/nconvfails counters for the nniS/nnfS updates.
 */
pub fn cvNlsSolveSensStg(
    cv_mem: &mut CVodeMem,
    nls: &mut NonlinearSolver,
    tol: f64,
    callLSetup: bool,
) -> i32 {
    match nls {
        NonlinearSolver::Newton(ns) => cvNlsSolveNewtonSensStg(cv_mem, ns, tol, callLSetup),
        NonlinearSolver::FixedPoint(fps) => cvNlsSolveFixedPointSensStg(cv_mem, fps, tol),
    }
}

/* zero the composite correction: N_VConst(ZERO, ycorStg) */
fn cvNlsZeroYcorStg(cv_mem: &mut CVodeMem) {
    for is in 0..cv_mem.cv_Ns as usize {
        N_VConst(ZERO, &mut cv_mem.cv_acorS[is]);
    }
}

/*
 * SUNNonlinSolSolve_Newton (sunnonlinsol_newton.c), specialized to the
 * CVODES staggered-corrector callbacks. ycor = cv_acorS, w = cv_ewtS;
 * the Newton update workspace is the senswrapper ns.deltaS (Ns
 * sub-vectors, from SUNNonlinSol_NewtonSens).
 */
fn cvNlsSolveNewtonSensStg(
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
            let r = cvNlsResidualSensStg(cv_mem, &mut deltaS);
            ns.deltaS = deltaS;
            r
        };
        if retval != 0 {
            break 'outer;
        }

        /* if indicated, setup the linear system */
        if call_lsetup {
            let mut jcur = ns.jcur;
            retval = cvNlsLSetupSensStg(cv_mem, jbad, &mut jcur);
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
            retval = cvNlsLSolveSensStg(cv_mem, &mut ns.deltaS);
            if retval != 0 {
                break;
            }

            /* update the Newton iterate: ycorStg = ycorStg + deltaStg */
            for is in 0..cv_mem.cv_Ns as usize {
                cv_mem.cv_acorS[is].linear_sum_with(ONE, ONE, &ns.deltaS.vecs[is]);
            }

            /* test for convergence */
            retval = cvNlsConvTestSensStg(cv_mem, &ns.deltaS, tol, ns.curiter);

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
                let r = cvNlsResidualSensStg(cv_mem, &mut deltaS);
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
            cvNlsZeroYcorStg(cv_mem);
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
 * to the CVODES staggered-corrector callbacks. ycor = cv_acorS; the
 * workspaces are the senswrappers fps.yprevS/gyS/deltaS (Ns
 * sub-vectors, from SUNNonlinSol_FixedPointSens).
 */
fn cvNlsSolveFixedPointSensStg(
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
        /* update previous solution guess: yprev = ycorStg */
        for is in 0..cv_mem.cv_Ns as usize {
            fps.yprevS.vecs[is]
                .data
                .copy_from_slice(&cv_mem.cv_acorS[is].data);
        }

        /* compute fixed-point iteration function, store in gy */
        let retval = cvNlsFPFunctionSensStg(cv_mem, &mut fps.gyS);
        if retval != 0 {
            return retval;
        }

        /* perform fixed point update, based on choice of acceleration or not */
        if fps.m == 0 {
            /* basic fixed-point solver: ycorStg = gy */
            for is in 0..cv_mem.cv_Ns as usize {
                cv_mem.cv_acorS[is]
                    .data
                    .copy_from_slice(&fps.gyS.vecs[is].data);
            }
        } else {
            /* Anderson-accelerated solver: assemble the composite iterate
               x = cv_acorS by moving the CVodeMem vectors into a
               senswrapper (no data copies), run the composite kernel, and
               move them back */
            let mut x = NVectorSensWrapper {
                vecs: std::mem::take(&mut cv_mem.cv_acorS),
                own_vecs: false,
            };
            let iter = fps.curiter;
            fps.anderson_accelerate_sens(&mut x, iter);
            cv_mem.cv_acorS = x.vecs;
        }

        /* increment nonlinear solver iteration counter */
        fps.niters += 1;

        /* compute change in solution: delta = ycorStg - yprev */
        {
            let FixedPointSolver { yprevS, deltaS, .. } = fps;
            for is in 0..cv_mem.cv_Ns as usize {
                N_VLinearSum(
                    ONE,
                    &cv_mem.cv_acorS[is],
                    -ONE,
                    &yprevS.vecs[is],
                    &mut deltaS.vecs[is],
                );
            }
        }

        /* test for convergence */
        let retval = {
            let m = fps.curiter;
            cvNlsConvTestSensStg(cv_mem, &fps.deltaS, tol, m)
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
