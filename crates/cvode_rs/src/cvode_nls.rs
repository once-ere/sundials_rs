/* -----------------------------------------------------------------
 * Translated from src/cvode/cvode_nls.c (CVODE 7.7.0), together
 * with the solve drivers of sunnonlinsol_newton.c and
 * sunnonlinsol_fixedpoint.c that the C code reaches through the
 * SUNNonlinearSolver ops table.  The CVODE-specific callbacks
 * (cvNlsResidual, cvNlsFPFunction, cvNlsLSetup, cvNlsLSolve,
 * cvNlsConvTest) are inlined into the drivers; control flow and
 * arithmetic order match the C sources statement for statement.
 * -----------------------------------------------------------------*/
use crate::cvode::{cv_has_lsetup, cv_lsetup_dispatch, cv_lsolve_dispatch};
use crate::cvode_impl::*;
use crate::nvector_serial::*;
use crate::sundials_math::*;
use crate::sundials_nonlinearsolver::*;
use crate::sundials_types::*;
use crate::sunnonlinsol_fixedpoint::FixedPointSolver;
use crate::sunnonlinsol_newton::NewtonSolver;

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/* nonlinear solver constants
   NLS_MAXCOR  maximum no. of corrector iterations for the nonlinear solver
   CRDOWN      constant used in the estimation of the convergence rate (crate)
   RDIV        declare divergence if ratio del/delp > RDIV */
const NLS_MAXCOR: i32 = 3;
const CRDOWN: f64 = 0.3;
const RDIV: f64 = 2.0;

/* -----------------------------------------------------------------------------
 * Exported functions
 * ---------------------------------------------------------------------------*/

pub fn CVodeSetNonlinearSolver(cv_mem: &mut CVodeMem, nls: NonlinearSolver) -> i32 {
    /* free any existing nonlinear solver: RAII (drop on overwrite) */

    /* set SUNNonlinearSolver pointer; the Sys function is selected by
       the solver type at solve time (cvNlsResidual / cvNlsFPFunction) */
    cv_mem.NLS = Some(nls);

    /* set ownership flag (CVODE resets to true when attaching its default) */
    cv_mem.ownNLS = SUNFALSE;

    /* set max allowed nonlinear iterations */
    if let Some(s) = cv_mem.NLS.as_mut() {
        let retval = s.set_max_iters(NLS_MAXCOR);
        if retval != 0 {
            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSetNonlinearSolver", file!(),
                           "Setting maximum number of nonlinear iterations failed");
            return CV_ILL_INPUT;
        }
    }

    /* Reset the acnrmcur flag to SUNFALSE */
    cv_mem.cv_acnrmcur = SUNFALSE;

    /* Set the nonlinear system RHS function */
    if cv_mem.cv_f.is_none() {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSetNonlinearSolver", file!(),
                       "The ODE RHS function is NULL");
        return CV_ILL_INPUT;
    }
    cv_mem.nls_f = cv_mem.cv_f;

    CV_SUCCESS
}

/*---------------------------------------------------------------
  CVodeSetNlsRhsFn:

  Sets an alternative user-supplied ODE right-hand side function
  to use in the evaluation of nonlinear system functions.
  ---------------------------------------------------------------*/
pub fn CVodeSetNlsRhsFn(cv_mem: &mut CVodeMem, f: Option<CVRhsFn>) -> i32 {
    match f {
        Some(func) => cv_mem.nls_f = Some(func),
        None => cv_mem.nls_f = cv_mem.cv_f,
    }
    CV_SUCCESS
}

/*---------------------------------------------------------------
  CVodeGetNonlinearSystemData:

  Provides access to the relevant data needed to compute the
  nonlinear system function (out-pointers become a returned tuple:
  (tcur, gamma, rl1); the vectors remain accessible as fields).
  ---------------------------------------------------------------*/
pub fn CVodeGetNonlinearSystemData(cv_mem: &CVodeMem) -> (f64, f64, f64) {
    (cv_mem.cv_tn, cv_mem.cv_gamma, cv_mem.cv_rl1)
}

/* -----------------------------------------------------------------------------
 * Private functions
 * ---------------------------------------------------------------------------*/

pub fn cvNlsInit(cv_mem: &mut CVodeMem) -> i32 {
    /* In C this wires cvNlsLSetup/cvNlsLSolve into the NLS depending on
       whether cv_lsetup/cv_lsolve exist; here the dispatch is dynamic.
       A Newton solver without an attached linear solver cannot work. */
    if let Some(nls) = cv_mem.NLS.as_ref() {
        if nls.nls_type() == SUNNONLINEARSOLVER_ROOTFIND && cv_mem.cv_lmem.is_none() {
            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "cvNlsInit", file!(), MSGCV_LSOLVE_NULL);
            return CV_NLS_INIT_FAIL;
        }
    }

    /* initialize nonlinear solver */
    let retval = match cv_mem.NLS.as_mut() {
        Some(nls) => nls.initialize(),
        None => -1,
    };

    if retval != CV_SUCCESS {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "cvNlsInit", file!(), MSGCV_NLS_INIT_FAIL);
        return CV_NLS_INIT_FAIL;
    }

    CV_SUCCESS
}

/* cvNlsLSetup (cvode_nls.c): wrapper around the lsetup dispatch */
fn cvNlsLSetup(cv_mem: &mut CVodeMem, jbad: bool, jcur: &mut bool) -> i32 {
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

    cv_mem.cv_gamrat = ONE;
    cv_mem.cv_gammap = cv_mem.cv_gamma;
    cv_mem.cv_crate = ONE;
    cv_mem.cv_nstlp = cv_mem.cv_nst;

    if retval < 0 {
        return CV_LSETUP_FAIL;
    }
    if retval > 0 {
        return SUN_NLS_CONV_RECVR;
    }

    CV_SUCCESS
}

/* cvNlsLSolve (cvode_nls.c): wrapper around the lsolve dispatch */
fn cvNlsLSolve(cv_mem: &mut CVodeMem, delta: &mut NVector) -> i32 {
    let retval = cv_lsolve_dispatch(cv_mem, delta);

    if retval < 0 {
        return CV_LSOLVE_FAIL;
    }
    if retval > 0 {
        return SUN_NLS_CONV_RECVR;
    }

    CV_SUCCESS
}

/* cvNlsConvTest (cvode_nls.c) */
fn cvNlsConvTest(
    cv_mem: &mut CVodeMem,
    ycor: &NVector,
    delta: &NVector,
    tol: f64,
    m: i32,
) -> i32 {
    /* compute the norm of the correction */
    let del = N_VWrmsNorm(delta, &cv_mem.cv_ewt);

    /* Test for convergence. If m > 0, an estimate of the convergence
       rate constant is stored in crate, and used in the test. */
    if m > 0 {
        cv_mem.cv_crate = SUNMAX(CRDOWN * cv_mem.cv_crate, del / cv_mem.cv_delp);
    }
    let dcon = del * SUNMIN(ONE, cv_mem.cv_crate) / tol;

    if dcon <= ONE {
        cv_mem.cv_acnrm = if m == 0 {
            del
        } else {
            N_VWrmsNorm(ycor, &cv_mem.cv_ewt)
        };
        cv_mem.cv_acnrmcur = SUNTRUE;
        return CV_SUCCESS; /* Nonlinear system was solved successfully */
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

/* cvNlsResidual (cvode_nls.c): evaluates the nonlinear residual
   res = ycor - gamma*f(t, zn[0]+ycor) + rl1*zn[1] with ycor = cv_acor. */
fn cvNlsResidual(cv_mem: &mut CVodeMem, res: &mut NVector) -> i32 {
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

    {
        let CVodeMem { cv_zn, cv_acor, cv_ftemp, cv_rl1, cv_gamma, .. } = cv_mem;
        N_VLinearSum(*cv_rl1, &cv_zn[1], ONE, cv_acor, res);
        /* res = -gamma*ftemp + res */
        res.linear_sum_with(ONE, -*cv_gamma, cv_ftemp);
    }

    CV_SUCCESS
}

/* cvNlsFPFunction (cvode_nls.c): fixed-point function evaluation
   res = rl1 * (h*f(t, zn[0]+ycor) - zn[1]) with ycor = cv_acor. */
fn cvNlsFPFunction(cv_mem: &mut CVodeMem, res: &mut NVector) -> i32 {
    /* update the state based on the current correction */
    {
        let CVodeMem { cv_zn, cv_acor, cv_y, .. } = cv_mem;
        N_VLinearSum(ONE, &cv_zn[0], ONE, cv_acor, cv_y);
    }

    /* evaluate the rhs function */
    let f = cv_mem.nls_f.unwrap();
    let retval = f(cv_mem.cv_tn, &cv_mem.cv_y, res, &mut cv_mem.cv_user_data);
    cv_mem.cv_nfe += 1;
    if retval < 0 {
        return CV_RHSFUNC_FAIL;
    }
    if retval > 0 {
        return RHSFUNC_RECVR;
    }

    {
        let CVodeMem { cv_zn, cv_h, cv_rl1, .. } = cv_mem;
        /* res = h*res - zn[1] */
        res.linear_sum_with(*cv_h, -ONE, &cv_zn[1]);
        res.scale_inplace(*cv_rl1);
    }

    CV_SUCCESS
}

/*
 * cvNls (cvode.c)
 *
 * Attempts to solve the nonlinear system associated with a single
 * implicit step of the linear multistep method.
 */
pub fn cvNls(cv_mem: &mut CVodeMem, nflag: i32) -> i32 {
    let call_setup;

    /* Decide whether or not to call setup routine (if one exists) and
       set flag convfail (input to lsetup for its evaluation decision) */
    if cv_has_lsetup(cv_mem) {
        cv_mem.convfail = if nflag == FIRST_CALL || nflag == PREV_ERR_FAIL {
            CV_NO_FAILURES
        } else {
            CV_FAIL_OTHER
        };

        call_setup = nflag == PREV_CONV_FAIL
            || nflag == PREV_ERR_FAIL
            || cv_mem.cv_nst == 0
            || cv_mem.first_step_after_resize
            || cv_mem.cv_nst >= cv_mem.cv_nstlp + cv_mem.cv_msbp
            || SUNRabs(cv_mem.cv_gamrat - ONE) > cv_mem.cv_dgmax_lsetup;
    } else {
        cv_mem.cv_crate = ONE;
        call_setup = SUNFALSE;
    }

    /* initial guess for the correction to the predictor */
    N_VConst(ZERO, &mut cv_mem.cv_acor);

    /* (Newton and fixed-point solvers have no setup operation) */

    /* solve the nonlinear system */
    let mut nls = cv_mem.NLS.take().expect("nonlinear solver attached");
    let tol = cv_mem.cv_tq[4];
    let flag = match &mut nls {
        NonlinearSolver::Newton(ns) => cvNlsSolveNewton(cv_mem, ns, tol, call_setup),
        NonlinearSolver::FixedPoint(fps) => cvNlsSolveFixedPoint(cv_mem, fps, tol),
    };

    /* increment counters */
    cv_mem.cv_nni += nls.get_num_iters();
    cv_mem.cv_nnf += nls.get_num_conv_fails();
    cv_mem.NLS = Some(nls);

    /* if the solve failed return */
    if flag != 0 {
        return flag;
    }

    /* solve successful */

    /* update the state based on the final correction from the nonlinear solver */
    {
        let CVodeMem { cv_zn, cv_acor, cv_y, .. } = cv_mem;
        N_VLinearSum(ONE, &cv_zn[0], ONE, cv_acor, cv_y);
    }

    /* compute acnrm if it was not already done by the nonlinear solver */
    if !cv_mem.cv_acnrmcur {
        cv_mem.cv_acnrm = N_VWrmsNorm(&cv_mem.cv_acor, &cv_mem.cv_ewt);
    }

    /* update Jacobian status */
    cv_mem.cv_jcur = SUNFALSE;

    flag
}

/*
 * SUNNonlinSolSolve_Newton (sunnonlinsol_newton.c), specialized to
 * the CVODE callbacks. ycor = cv_mem.cv_acor, w = cv_mem.cv_ewt.
 */
fn cvNlsSolveNewton(
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
        retval = cvNlsResidual(cv_mem, &mut ns.delta);
        if retval != 0 {
            break 'outer;
        }

        /* if indicated, setup the linear system */
        if call_lsetup {
            let mut jcur = ns.jcur;
            retval = cvNlsLSetup(cv_mem, jbad, &mut jcur);
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
            retval = cvNlsLSolve(cv_mem, &mut ns.delta);
            if retval != 0 {
                break;
            }

            /* update the Newton iterate */
            cv_mem.cv_acor.linear_sum_with(ONE, ONE, &ns.delta);

            /* test for convergence */
            retval = {
                let CVodeMem { cv_acor, .. } = cv_mem;
                let acor = std::mem::take(cv_acor);
                let r = cvNlsConvTest(cv_mem, &acor, &ns.delta, tol, ns.curiter);
                cv_mem.cv_acor = acor;
                r
            };

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
            retval = cvNlsResidual(cv_mem, &mut ns.delta);
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
            N_VConst(ZERO, &mut cv_mem.cv_acor);
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
 * SUNNonlinSolSolve_FixedPoint (sunnonlinsol_fixedpoint.c),
 * specialized to the CVODE callbacks.  ycor = cv_mem.cv_acor.
 */
fn cvNlsSolveFixedPoint(cv_mem: &mut CVodeMem, fps: &mut FixedPointSolver, tol: f64) -> i32 {
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
        /* update previous solution guess */
        fps.yprev.data.copy_from_slice(&cv_mem.cv_acor.data);

        /* compute fixed-point iteration function, store in gy */
        {
            let mut gy = std::mem::take(&mut fps.gy);
            let retval = cvNlsFPFunction(cv_mem, &mut gy);
            fps.gy = gy;
            if retval != 0 {
                return retval;
            }
        }

        /* perform fixed point update, based on choice of acceleration or not */
        if fps.m == 0 {
            /* basic fixed-point solver */
            cv_mem.cv_acor.data.copy_from_slice(&fps.gy.data);
        } else {
            /* Anderson-accelerated solver */
            let mut acor = std::mem::take(&mut cv_mem.cv_acor);
            let iter = fps.curiter;
            fps.anderson_accelerate(&mut acor, iter);
            cv_mem.cv_acor = acor;
        }

        /* increment nonlinear solver iteration counter */
        fps.niters += 1;

        /* compute change in solution, and call the convergence test function */
        {
            let CVodeMem { cv_acor, .. } = cv_mem;
            N_VLinearSum(ONE, cv_acor, -ONE, &fps.yprev, &mut fps.delta);
        }

        /* test for convergence */
        let retval = {
            let acor = std::mem::take(&mut cv_mem.cv_acor);
            let m = fps.curiter;
            let r = cvNlsConvTest(cv_mem, &acor, &fps.delta, tol, m);
            cv_mem.cv_acor = acor;
            r
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
