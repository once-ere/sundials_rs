/* -----------------------------------------------------------------
 * Translated from src/cvodes/cvodes_nls_stg1.c (CVODES 7.7.0),
 * together with the solve drivers of sunnonlinsol_newton.c and
 * sunnonlinsol_fixedpoint.c that the C code reaches through the
 * SUNNonlinearSolver ops table: the STAGGERED1 corrector, which
 * solves the Ns sensitivity systems one parameter at a time (after
 * the state system has converged), the current parameter index
 * living in cv_mem.sens_solve_idx.
 *
 * No senswrappers are involved (the iteration vectors are single
 * N_Vectors): in C the solver is a plain SUNNonlinSol_Newton /
 * SUNNonlinSol_FixedPoint built on the state template, and cvodes.c
 * passes ycor = cv_acorS[is], w = cv_ewtS[is] for each solve. Per
 * pinned decision 3 (cvodes_impl.rs) this module reads
 * cv_znS[0]/cv_znS[1], cv_acorS, cv_ewtS at index is =
 * cv_mem.sens_solve_idx directly.
 *
 * Exported solve entry (called by cvodes.c's cvStgr1Nls, ported in a
 * later part): cvNlsSolveSensStg1(cv_mem, nls, tol, callLSetup) is
 * the specialization of
 *   SUNNonlinSolSolve(NLSstg1, znS[0][is], acorS[is], ewtS[is], tol,
 *                     callLSetup, cv_mem)
 * — the caller must set cv_mem.sens_solve_idx = is beforehand (as
 * cvStgr1Nls does in C) and detaches the solver
 * (cv_mem.NLSstg1.take()) for the duration of the call.
 *
 * Forward reference to a cvodes.c symbol ported in a later part
 * (assumed signature, pinned for the Part 2/3 briefs):
 *   crate::cvodes::cvSensRhs1Wrapper(cv_mem: &mut CVodeMem, time: f64,
 *       ycur: &NVector, fcur: &NVector, is: i32, yScur: &NVector,
 *       fScur: &mut NVector, temp1: &mut NVector,
 *       temp2: &mut NVector) -> i32   (increments cv_nfSe itself)
 * -----------------------------------------------------------------*/
use crate::cvodes_impl::*;
use crate::cvodes_nls::{cv_has_lsetup, cv_lsetup_dispatch, cv_lsolve_dispatch};
use crate::nvector_serial::*;
use crate::sundials_math::*;
use crate::sundials_nonlinearsolver::*;
use crate::sundials_types::*;
use crate::sunnonlinsol_fixedpoint::FixedPointSolver;
use crate::sunnonlinsol_newton::NewtonSolver;

/* constant macros */
const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/* -----------------------------------------------------------------------------
 * Exported functions
 * ---------------------------------------------------------------------------*/

pub fn CVodeSetNonlinearSolverSensStg1(cv_mem: &mut CVodeMem, nls: NonlinearSolver) -> i32 {
    /* check that sensitivities were initialized */
    if !cv_mem.cv_sensi {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSetNonlinearSolverSensStg1",
                       file!(), MSGCV_NO_SENSI);
        return CV_ILL_INPUT;
    }

    /* check that staggered corrector was selected */
    if cv_mem.cv_ism != CV_STAGGERED1 {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSetNonlinearSolverSensStg1",
                       file!(), "Sensitivity solution method is not CV_STAGGERED1");
        return CV_ILL_INPUT;
    }

    /* free any existing nonlinear solver: RAII (drop on overwrite) */

    /* set SUNNonlinearSolver pointer; the Sys function is selected by the
       solver type at solve time (cvNlsResidualSensStg1 /
       cvNlsFPFunctionSensStg1), and the convergence test is inlined in the
       solve drivers below */
    cv_mem.NLSstg1 = Some(nls);

    /* Set NLS ownership flag. If this function was called to attach the
       default NLS, CVODES will set the flag to SUNTRUE after this function
       returns. */
    cv_mem.ownNLSstg1 = SUNFALSE;

    /* set max allowed nonlinear iterations */
    if let Some(s) = cv_mem.NLSstg1.as_mut() {
        let retval = s.set_max_iters(NLS_MAXCOR);
        if retval != CV_SUCCESS {
            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(),
                           "CVodeSetNonlinearSolverSensStg1", file!(),
                           "Setting maximum number of nonlinear iterations failed");
            return CV_ILL_INPUT;
        }
    }

    /* Reset the acnrmScur flag to SUNFALSE (always false for stg1) */
    cv_mem.cv_acnrmScur = SUNFALSE;

    CV_SUCCESS
}

/* -----------------------------------------------------------------------------
 * Private functions
 * ---------------------------------------------------------------------------*/

pub fn cvNlsInitSensStg1(cv_mem: &mut CVodeMem) -> i32 {
    /* In C this wires cvNlsLSetupSensStg1/cvNlsLSolveSensStg1 into the NLS
       depending on whether cv_lsetup/cv_lsolve exist; here the dispatch is
       dynamic. A Newton solver without an attached linear solver cannot
       work (SUNNonlinSolSolve_Newton requires an LSolve function). */
    if let Some(nls) = cv_mem.NLSstg1.as_ref() {
        if nls.nls_type() == SUNNONLINEARSOLVER_ROOTFIND && cv_mem.cv_lmem.is_none() {
            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "cvNlsInitSensStg1", file!(),
                           MSGCV_LSOLVE_NULL);
            return CV_NLS_INIT_FAIL;
        }
    }

    /* initialize nonlinear solver */
    let retval = match cv_mem.NLSstg1.as_mut() {
        Some(nls) => nls.initialize(),
        None => -1,
    };

    if retval != CV_SUCCESS {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "cvNlsInitSensStg1", file!(),
                       MSGCV_NLS_INIT_FAIL);
        return CV_NLS_INIT_FAIL;
    }

    /* reset previous iteration count for updating nniS1 */
    cv_mem.nnip = 0;

    CV_SUCCESS
}

/* cvNlsLSetupSensStg1 (cvodes_nls_stg1.c): wrapper around the lsetup
   dispatch (also counts the setup against the sensitivities; does not
   touch cv_forceSetup) */
fn cvNlsLSetupSensStg1(cv_mem: &mut CVodeMem, jbad: bool, jcur: &mut bool) -> i32 {
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

/* cvNlsLSolveSensStg1 (cvodes_nls_stg1.c): solves the linear system for
   the current sensitivity (weight = cv_ewtS[sens_solve_idx]) */
fn cvNlsLSolveSensStg1(cv_mem: &mut CVodeMem, delta: &mut NVector) -> i32 {
    /* get index of current sensitivity solve */
    let is = cv_mem.sens_solve_idx as usize;

    /* solve the sensitivity linear system */
    let retval = cv_lsolve_dispatch(cv_mem, delta, Some(is));

    if retval < 0 {
        return CV_LSOLVE_FAIL;
    }
    if retval > 0 {
        return SUN_NLS_CONV_RECVR;
    }

    CV_SUCCESS
}

/* cvNlsConvTestSensStg1 (cvodes_nls_stg1.c); ewt = cv_ewtS[sens_solve_idx],
   m = current iteration count. The ycor argument is unused in C
   (SUNDIALS_MAYBE_UNUSED) and dropped here; no acnrm update is made. */
fn cvNlsConvTestSensStg1(cv_mem: &mut CVodeMem, delta: &NVector, tol: f64, m: i32) -> i32 {
    /* compute the norm of the state and sensitivity corrections */
    let is = cv_mem.sens_solve_idx as usize;
    let del = N_VWrmsNorm(delta, &cv_mem.cv_ewtS[is]);

    /* Test for convergence. If m > 0, an estimate of the convergence
       rate constant is stored in crate, and used in the test.
    */
    if m > 0 {
        cv_mem.cv_crateS = SUNMAX(CRDOWN * cv_mem.cv_crateS, del / cv_mem.cv_delp);
    }
    let dcon = del * SUNMIN(ONE, cv_mem.cv_crateS) / tol;

    /* check if nonlinear system was solved successfully */
    if dcon <= ONE {
        return CV_SUCCESS;
    }

    /* check if the iteration seems to be diverging */
    if m >= 1 && del > RDIV * cv_mem.cv_delp {
        return SUN_NLS_CONV_RECVR;
    }

    /* Save norm of correction and loop again */
    cv_mem.cv_delp = del;

    /* Not yet converged */
    SUN_NLS_CONTINUE
}

/* cvNlsResidualSensStg1 (cvodes_nls_stg1.c): evaluates the residual of
   the current (is = sens_solve_idx) sensitivity system;
   ycor = cv_acorS[is]. The state values cv_y and cv_ftemp hold the
   already-converged state solution. */
fn cvNlsResidualSensStg1(cv_mem: &mut CVodeMem, res: &mut NVector) -> i32 {
    /* get index of current sensitivity solve */
    let is = cv_mem.sens_solve_idx as usize;

    /* update sensitivity based on the current correction */
    {
        let CVodeMem { cv_znS, cv_acorS, cv_yS, .. } = cv_mem;
        N_VLinearSum(ONE, &cv_znS[0][is], ONE, &cv_acorS[is], &mut cv_yS[is]);
    }

    /* evaluate the sensitivity rhs function (cvSensRhs1Wrapper increments
       cv_nfSe itself; the argument vectors are CVodeMem fields, taken out
       for the duration of the call — donor take()/restore pattern) */
    let cv_y = std::mem::take(&mut cv_mem.cv_y);
    let cv_ftemp = std::mem::take(&mut cv_mem.cv_ftemp);
    let yS_is = std::mem::take(&mut cv_mem.cv_yS[is]);
    let mut ftempS_is = std::mem::take(&mut cv_mem.cv_ftempS[is]);
    let mut vtemp1 = std::mem::take(&mut cv_mem.cv_vtemp1);
    let mut vtemp2 = std::mem::take(&mut cv_mem.cv_vtemp2);
    let tn = cv_mem.cv_tn;
    let retval = crate::cvodes::cvSensRhs1Wrapper(
        cv_mem,
        tn,
        &cv_y,
        &cv_ftemp,
        is as i32,
        &yS_is,
        &mut ftempS_is,
        &mut vtemp1,
        &mut vtemp2,
    );
    cv_mem.cv_y = cv_y;
    cv_mem.cv_ftemp = cv_ftemp;
    cv_mem.cv_yS[is] = yS_is;
    cv_mem.cv_ftempS[is] = ftempS_is;
    cv_mem.cv_vtemp1 = vtemp1;
    cv_mem.cv_vtemp2 = vtemp2;

    if retval < 0 {
        return CV_SRHSFUNC_FAIL;
    }
    if retval > 0 {
        return SRHSFUNC_RECVR;
    }

    /* compute the sensitivity residual */
    {
        let CVodeMem { cv_znS, cv_acorS, cv_ftempS, cv_rl1, cv_gamma, .. } = cv_mem;
        N_VLinearSum(*cv_rl1, &cv_znS[1][is], ONE, &cv_acorS[is], res);
        /* res = -gamma*ftempS + res */
        res.linear_sum_with(ONE, -*cv_gamma, &cv_ftempS[is]);
    }

    CV_SUCCESS
}

/* cvNlsFPFunctionSensStg1 (cvodes_nls_stg1.c): fixed-point function for
   the current (is = sens_solve_idx) sensitivity system;
   ycor = cv_acorS[is]. */
fn cvNlsFPFunctionSensStg1(cv_mem: &mut CVodeMem, res: &mut NVector) -> i32 {
    /* get index of current sensitivity solve */
    let is = cv_mem.sens_solve_idx as usize;

    /* update the sensitivities based on the current correction */
    {
        let CVodeMem { cv_znS, cv_acorS, cv_yS, .. } = cv_mem;
        N_VLinearSum(ONE, &cv_znS[0][is], ONE, &cv_acorS[is], &mut cv_yS[is]);
    }

    /* evaluate the sensitivity rhs function (fScur = res) */
    let cv_y = std::mem::take(&mut cv_mem.cv_y);
    let cv_ftemp = std::mem::take(&mut cv_mem.cv_ftemp);
    let yS_is = std::mem::take(&mut cv_mem.cv_yS[is]);
    let mut vtemp1 = std::mem::take(&mut cv_mem.cv_vtemp1);
    let mut vtemp2 = std::mem::take(&mut cv_mem.cv_vtemp2);
    let tn = cv_mem.cv_tn;
    let retval = crate::cvodes::cvSensRhs1Wrapper(
        cv_mem,
        tn,
        &cv_y,
        &cv_ftemp,
        is as i32,
        &yS_is,
        res,
        &mut vtemp1,
        &mut vtemp2,
    );
    cv_mem.cv_y = cv_y;
    cv_mem.cv_ftemp = cv_ftemp;
    cv_mem.cv_yS[is] = yS_is;
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
        let CVodeMem { cv_znS, cv_h, cv_rl1, .. } = cv_mem;
        /* res = h*res - znS[1] */
        res.linear_sum_with(*cv_h, -ONE, &cv_znS[1][is]);
        res.scale_inplace(*cv_rl1);
    }

    CV_SUCCESS
}

/*
 * cvNlsSolveSensStg1 — SUNNonlinSolSolve(cv_mem->NLSstg1,
 * cv_mem->cv_znS[0][is], cv_mem->cv_acorS[is], cv_mem->cv_ewtS[is], tol,
 * callLSetup, cv_mem) as invoked by cvodes.c's cvStgr1Nls (ported in a
 * later part) with is = cv_mem.sens_solve_idx (set by the caller before
 * the call). The caller detaches the solver from cv_mem
 * (NLSstg1.take()) and reattaches it afterwards, reading the
 * niters/nconvfails counters (together with cv_mem.nnip) for the
 * nniS1/nnfS1 updates.
 */
pub fn cvNlsSolveSensStg1(
    cv_mem: &mut CVodeMem,
    nls: &mut NonlinearSolver,
    tol: f64,
    callLSetup: bool,
) -> i32 {
    match nls {
        NonlinearSolver::Newton(ns) => cvNlsSolveNewtonSensStg1(cv_mem, ns, tol, callLSetup),
        NonlinearSolver::FixedPoint(fps) => cvNlsSolveFixedPointSensStg1(cv_mem, fps, tol),
    }
}

/*
 * SUNNonlinSolSolve_Newton (sunnonlinsol_newton.c), specialized to the
 * CVODES staggered1-corrector callbacks. ycor = cv_acorS[is],
 * w = cv_ewtS[is] with is = cv_mem.sens_solve_idx; the solver is a plain
 * SUNNonlinSol_Newton (state-length delta workspace).
 */
fn cvNlsSolveNewtonSensStg1(
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
        retval = cvNlsResidualSensStg1(cv_mem, &mut ns.delta);
        if retval != 0 {
            break 'outer;
        }

        /* if indicated, setup the linear system */
        if call_lsetup {
            let mut jcur = ns.jcur;
            retval = cvNlsLSetupSensStg1(cv_mem, jbad, &mut jcur);
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
            ns.delta.scale_inplace(-ONE);

            /* solve the linear system to get Newton update delta */
            retval = cvNlsLSolveSensStg1(cv_mem, &mut ns.delta);
            if retval != 0 {
                break;
            }

            /* update the Newton iterate: ycor = ycor + delta */
            {
                let is = cv_mem.sens_solve_idx as usize;
                cv_mem.cv_acorS[is].linear_sum_with(ONE, ONE, &ns.delta);
            }

            /* test for convergence */
            retval = cvNlsConvTestSensStg1(cv_mem, &ns.delta, tol, ns.curiter);

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
            retval = cvNlsResidualSensStg1(cv_mem, &mut ns.delta);
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
            let is = cv_mem.sens_solve_idx as usize;
            N_VConst(ZERO, &mut cv_mem.cv_acorS[is]);
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
 * to the CVODES staggered1-corrector callbacks. ycor = cv_acorS[is] with
 * is = cv_mem.sens_solve_idx; the solver is a plain
 * SUNNonlinSol_FixedPoint (state-length workspaces).
 */
fn cvNlsSolveFixedPointSensStg1(
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
        /* update previous solution guess: yprev = ycor */
        {
            let is = cv_mem.sens_solve_idx as usize;
            fps.yprev.data.copy_from_slice(&cv_mem.cv_acorS[is].data);
        }

        /* compute fixed-point iteration function, store in gy */
        let retval = cvNlsFPFunctionSensStg1(cv_mem, &mut fps.gy);
        if retval != 0 {
            return retval;
        }

        /* perform fixed point update, based on choice of acceleration or not */
        let is = cv_mem.sens_solve_idx as usize;
        if fps.m == 0 {
            /* basic fixed-point solver: ycor = gy */
            cv_mem.cv_acorS[is].data.copy_from_slice(&fps.gy.data);
        } else {
            /* Anderson-accelerated solver */
            let mut acorS_is = std::mem::take(&mut cv_mem.cv_acorS[is]);
            let iter = fps.curiter;
            fps.anderson_accelerate(&mut acorS_is, iter);
            cv_mem.cv_acorS[is] = acorS_is;
        }

        /* increment nonlinear solver iteration counter */
        fps.niters += 1;

        /* compute change in solution: delta = ycor - yprev */
        {
            let FixedPointSolver { yprev, delta, .. } = fps;
            N_VLinearSum(ONE, &cv_mem.cv_acorS[is], -ONE, yprev, delta);
        }

        /* test for convergence */
        let retval = {
            let m = fps.curiter;
            cvNlsConvTestSensStg1(cv_mem, &fps.delta, tol, m)
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
