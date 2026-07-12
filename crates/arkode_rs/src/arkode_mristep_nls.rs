/* -----------------------------------------------------------------
 * Translation of src/arkode/arkode_mristep_nls.c (SUNDIALS 7.7.0).
 *
 * This is the interface between MRIStep and the SUNNonlinearSolver
 * object.  Following the crate convention (see
 * arkode_arkstep_nls.rs), the generic SUNNonlinSolSolve loops are
 * INLINED here (Newton + accelerated fixed-point), with the core
 * NewtonSolver/FixedPointSolver structs holding only state and
 * counters.
 * -----------------------------------------------------------------*/

use crate::arkode_impl::{
    arkProcessError, ARKRhsFn, ARKodeMem, ARK_FAIL_BAD_J, ARK_FAIL_OTHER, ARK_ILL_INPUT,
    ARK_LSETUP_FAIL, ARK_LSOLVE_FAIL, ARK_MEM_NULL, ARK_NLS_INIT_FAIL, ARK_NO_FAILURES,
    ARK_PRERHSFN_FAIL, ARK_RHSFUNC_FAIL, ARK_SUCCESS, CONV_FAIL, FIRST_CALL, PREV_CONV_FAIL,
    PREV_ERR_FAIL, RHSFUNC_RECVR, ZERO,
};
use crate::arkode_mristep::{mriStep_AccessStepMem, with_step_mem_installed};
use crate::arkode_mristep_impl::{ARKodeMRIStepMem, MSG_NLS_INIT_FAIL};
use crate::nvector_serial::{NVector, N_VConst, N_VLinearCombination, N_VLinearSum, N_VWrmsNorm};
use crate::sundials_errors::SUN_SUCCESS;
use crate::sundials_math::{SUNMAX, SUNMIN, SUNRabs};
use crate::sundials_nonlinearsolver::{NonlinearSolver, SUN_NLS_CONTINUE, SUN_NLS_CONV_RECVR};
use crate::sunnonlinsol_fixedpoint::FixedPointSolver;
use crate::sunnonlinsol_newton::NewtonSolver;

const ONE: f64 = 1.0;
const SUNFALSE: bool = false;
const SUNTRUE: bool = true;

/*===============================================================
  Interface routines supplied to ARKODE
  ===============================================================*/

/*---------------------------------------------------------------
  mriStep_SetNonlinearSolver:

  This routine attaches a SUNNonlinearSolver object to the MRIStep
  module.
  ---------------------------------------------------------------*/
pub fn mriStep_SetNonlinearSolver(ark_mem: &mut ARKodeMem, NLS: NonlinearSolver) -> i32 {
    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_SetNonlinearSolver") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    /* (C checks the NLS ops table for gettype/solve/setsysfn; the Rust
       NonlinearSolver enum always provides them.  The nonlinear
       residual/fixed-point function and convergence test are dispatched
       directly by the inlined solve loop.) */

    /* free any existing nonlinear solver; set SUNNonlinearSolver pointer */
    step_mem.NLS = Some(NLS);
    step_mem.ownNLS = false;

    /* set default nonlinear iterations */
    let maxcor = step_mem.maxcor;
    step_mem.NLS.as_mut().unwrap().set_max_iters(maxcor);

    /* set the nonlinear system RHS function */
    step_mem.nls_fsi = None;

    if step_mem.implicit_rhs {
        if step_mem.fsi.is_none() {
            ark_mem.step_mem = Some(step_mem);
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "mriStep_SetNonlinearSolver",
                file!(),
                "The implicit slow ODE RHS function is NULL",
            );
            return ARK_ILL_INPUT;
        }
        step_mem.nls_fsi = step_mem.fsi;
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_SetNlsRhsFn:

  This routine sets an alternative user-supplied slow ODE
  right-hand side function to use in the evaluation of nonlinear
  system functions.
  ---------------------------------------------------------------*/
pub fn mriStep_SetNlsRhsFn(ark_mem: &mut ARKodeMem, nls_fsi: Option<ARKRhsFn>) -> i32 {
    let mut step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_SetNlsRhsFn") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    if nls_fsi.is_some() {
        step_mem.nls_fsi = nls_fsi;
    } else {
        step_mem.nls_fsi = step_mem.fsi;
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_GetNonlinearSystemData:

  This routine provides access to the relevant data needed to
  compute the nonlinear system function.  (C hands out internal
  pointers; the vector out-parameters receive CLONES -- see the
  ARKTimestepGetNonlinearSystemData op note in arkode_impl.rs.)
  ---------------------------------------------------------------*/
pub fn mriStep_GetNonlinearSystemData(
    ark_mem: &mut ARKodeMem,
    tcur: &mut f64,
    zpred: &mut NVector,
    z: &mut NVector,
    F: &mut NVector,
    gamma: &mut f64,
    sdata: &mut NVector,
) -> i32 {
    let step_mem = match mriStep_AccessStepMem(ark_mem, "mriStep_GetNonlinearSystemData") {
        Some(sm) => sm,
        None => return ARK_MEM_NULL,
    };

    *tcur = ark_mem.tcur;
    zpred.data.clone_from(&step_mem.zpred.data);
    z.data.clone_from(&ark_mem.ycur.data);
    {
        let smap = step_mem.stage_map[step_mem.istage as usize] as usize;
        let fsi_arr = if step_mem.unify_Fs {
            &step_mem.Fse
        } else {
            &step_mem.Fsi
        };
        F.data.clone_from(&fsi_arr[smap].data);
    }
    *gamma = step_mem.gamma;
    sdata.data.clone_from(&step_mem.sdata.data);

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*===============================================================
  Utility routines called by MRIStep
  ===============================================================*/

/*---------------------------------------------------------------
  mriStep_NlsInit:

  This routine initializes the nonlinear solver object.  This
  should only be called at the start of a simulation, after a
  re-init, or after a re-size.  (The lsetup/lsolve wrapper
  functions are dispatched by the inlined solve loop on
  step_mem.lsetup / step_mem.lsolve presence.)
  ---------------------------------------------------------------*/
pub(crate) fn mriStep_NlsInit(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeMRIStepMem>,
) -> i32 {
    /* reset counters */
    step_mem.nls_iters = 0;
    step_mem.nls_fails = 0;

    /* initialize nonlinear solver */
    let retval = step_mem.NLS.as_mut().map(|n| n.initialize()).unwrap_or(0);
    if retval != 0 {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "mriStep_NlsInit",
            file!(),
            MSG_NLS_INIT_FAIL,
        );
        return ARK_NLS_INIT_FAIL;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_Nls

  This routine attempts to solve the nonlinear system associated
  with a single solve-decoupled implicit stage.

  Upon entry, the predicted solution is held in step_mem->zpred,
  which is never changed throughout this routine.

  Upon a successful solve, the solution is held in ark_mem->ycur.
  ---------------------------------------------------------------*/
pub(crate) fn mriStep_Nls(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeMRIStepMem>,
    nflag: i32,
) -> i32 {
    /* If a linear solver 'setup' is supplied, set various flags for
       determining whether it should be called */
    let call_lsetup;
    if step_mem.lsetup.is_some() {
        /* Set interface 'convfail' flag for use inside lsetup */
        if step_mem.linear {
            step_mem.convfail = if nflag == FIRST_CALL {
                ARK_NO_FAILURES
            } else {
                ARK_FAIL_OTHER
            };
        } else {
            step_mem.convfail = if nflag == FIRST_CALL || nflag == PREV_ERR_FAIL {
                ARK_NO_FAILURES
            } else {
                ARK_FAIL_OTHER
            };
        }

        /* Decide whether to recommend call to lsetup within nonlinear solver */
        let mut cls = ark_mem.firststage
            || step_mem.msbp < 0
            || SUNRabs(step_mem.gamrat - ONE) > step_mem.dgmax;
        if step_mem.linear {
            /* linearly-implicit problem */
            cls = cls || step_mem.linear_timedep;
        } else {
            /* nonlinearly-implicit problem */
            cls = cls
                || nflag == PREV_CONV_FAIL
                || nflag == PREV_ERR_FAIL
                || ark_mem.nst >= step_mem.nstlp + (step_mem.msbp.abs() as i64);
        }
        call_lsetup = cls;
    } else {
        step_mem.crate_ = ONE;
        call_lsetup = false;
    }

    /* set a zero guess for correction */
    N_VConst(ZERO, &mut step_mem.zcor);

    /* Reset the stored residual norm (for iterative linear solvers) */
    step_mem.eRNrm = 0.1 * step_mem.nlscoef;

    /* solve the nonlinear system for the actual correction; take the
       nonlinear solver out of step_mem for the iteration */
    let nls = match step_mem.NLS.take() {
        Some(n) => n,
        None => return ARK_MEM_NULL,
    };
    let (retval, nls) = match nls {
        NonlinearSolver::Newton(mut ns) => {
            let tol = step_mem.nlscoef;
            let ret = mriStep_NlsSolveNewton(ark_mem, step_mem, &mut ns, tol, call_lsetup);
            (ret, NonlinearSolver::Newton(ns))
        }
        NonlinearSolver::FixedPoint(mut fps) => {
            let tol = step_mem.nlscoef;
            let ret = mriStep_NlsSolveFixedPoint(ark_mem, step_mem, &mut fps, tol);
            (ret, NonlinearSolver::FixedPoint(fps))
        }
    };

    /* increment counters */
    let (iters_inc, fails_inc) = (nls.get_num_iters(), nls.get_num_conv_fails());
    step_mem.NLS = Some(nls);
    step_mem.nls_iters += iters_inc;
    step_mem.nls_fails += fails_inc;

    /* successful solve -- reset the jcur flag and apply correction */
    if retval == SUN_SUCCESS {
        step_mem.jcur = false;
        let ARKodeMRIStepMem { zcor, zpred, .. } = &mut **step_mem;
        N_VLinearSum(ONE, zcor, ONE, zpred, &mut ark_mem.ycur);
        return ARK_SUCCESS;
    }

    /* check for recoverable failure, return ARKODE::CONV_FAIL */
    if retval == SUN_NLS_CONV_RECVR {
        return CONV_FAIL;
    }

    retval
}

/*---------------------------------------------------------------
  mriStep_NlsSolveNewton: the generic SUNNonlinSolSolve_Newton
  loop (sunnonlinsol_newton.c), specialized to the MRIStep
  callbacks.  ycor = step_mem->zcor, w = ark_mem->ewt.
  ---------------------------------------------------------------*/
fn mriStep_NlsSolveNewton(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeMRIStepMem>,
    ns: &mut NewtonSolver,
    tol: f64,
    mut call_lsetup: bool,
) -> i32 {
    /* assume the Jacobian is good */
    let mut jbad = SUNFALSE;

    /* initialize iteration and convergence fail counters for this solve */
    ns.niters = 0;
    ns.nconvfails = 0;

    let mut retval: i32;
    'outer: loop {
        /* initialize current iteration counter for this solve attempt */
        ns.curiter = 0;

        /* compute the nonlinear residual, store in delta */
        retval = mriStep_NlsResidual(ark_mem, step_mem, &mut ns.delta);
        if retval != 0 {
            break 'outer;
        }

        /* if indicated, setup the linear system */
        if call_lsetup {
            let mut jcur = ns.jcur;
            retval = mriStep_NlsLSetup(ark_mem, step_mem, jbad, &mut jcur);
            ns.jcur = jcur;
            if retval != 0 {
                break 'outer;
            }
        }

        /* looping point for Newton iteration. Break out on any error. */
        loop {
            /* increment nonlinear solver iteration counter */
            ns.niters += 1;

            /* compute the negative of the residual for the linear system rhs */
            ns.delta.scale_inplace(-ONE);

            /* solve the linear system to get Newton update delta */
            retval = mriStep_NlsLSolve(ark_mem, step_mem, ns.curiter, &mut ns.delta);
            if retval != 0 {
                break;
            }

            /* update the Newton iterate */
            step_mem.zcor.linear_sum_with(ONE, ONE, &ns.delta);

            /* test for convergence */
            retval = mriStep_NlsConvTest(ark_mem, step_mem, &ns.delta, tol, ns.curiter);

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
            ns.curiter += 1;
            if ns.curiter >= ns.maxiters {
                retval = SUN_NLS_CONV_RECVR;
                break;
            }

            /* compute the nonlinear residual, store in delta */
            retval = mriStep_NlsResidual(ark_mem, step_mem, &mut ns.delta);
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
        if retval > 0 && !ns.jcur && step_mem.lsetup.is_some() {
            ns.nconvfails += 1;
            call_lsetup = SUNTRUE;
            jbad = SUNTRUE;
            N_VConst(ZERO, &mut step_mem.zcor);
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

/*---------------------------------------------------------------
  mriStep_NlsSolveFixedPoint: the generic
  SUNNonlinSolSolve_FixedPoint loop (sunnonlinsol_fixedpoint.c),
  specialized to the MRIStep callbacks.
  ---------------------------------------------------------------*/
fn mriStep_NlsSolveFixedPoint(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeMRIStepMem>,
    fps: &mut FixedPointSolver,
    tol: f64,
) -> i32 {
    /* initialize iteration and convergence fail counters for this solve */
    fps.niters = 0;
    fps.nconvfails = 0;

    fps.curiter = 0;
    while fps.curiter < fps.maxiters {
        /* update previous solution guess */
        fps.yprev.data.copy_from_slice(&step_mem.zcor.data);

        /* compute fixed-point iteration function, store in gy */
        {
            let mut gy = std::mem::take(&mut fps.gy);
            let retval = mriStep_NlsFPFunction(ark_mem, step_mem, &mut gy);
            fps.gy = gy;
            if retval != 0 {
                return retval;
            }
        }

        /* perform fixed point update, based on choice of acceleration or not */
        if fps.m == 0 {
            /* basic fixed-point solver */
            step_mem.zcor.data.copy_from_slice(&fps.gy.data);
        } else {
            /* Anderson-accelerated solver */
            let mut zcor = std::mem::take(&mut step_mem.zcor);
            let iter = fps.curiter;
            fps.anderson_accelerate(&mut zcor, iter);
            step_mem.zcor = zcor;
        }

        /* increment nonlinear solver iteration counter */
        fps.niters += 1;

        /* compute change in solution */
        {
            let ARKodeMRIStepMem { zcor, .. } = &mut **step_mem;
            N_VLinearSum(ONE, zcor, -ONE, &fps.yprev, &mut fps.delta);
        }

        /* test for convergence */
        let retval = {
            let delta = std::mem::take(&mut fps.delta);
            let r = mriStep_NlsConvTest(ark_mem, step_mem, &delta, tol, fps.curiter);
            fps.delta = delta;
            r
        };

        /* return if successful */
        if retval == SUN_SUCCESS {
            return SUN_SUCCESS;
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

/*===============================================================
  Interface routines supplied to the SUNNonlinearSolver module
  ===============================================================*/

/*---------------------------------------------------------------
  mriStep_NlsLSetup:

  This routine wraps the ARKODE linear solver interface 'setup'
  routine for use by the nonlinear solver object.
  ---------------------------------------------------------------*/
fn mriStep_NlsLSetup(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeMRIStepMem>,
    jbad: bool,
    jcur: &mut bool,
) -> i32 {
    /* update convfail based on jbad flag */
    if jbad {
        step_mem.convfail = ARK_FAIL_BAD_J;
    }

    /* Use ARKODE's tempv1, tempv2 and tempv3 as
       temporary vectors for the linear solver setup routine */
    step_mem.nsetups += 1;

    let lsetup = step_mem.lsetup.unwrap();
    let convfail = step_mem.convfail;
    let smap = step_mem.stage_map[step_mem.istage as usize] as usize;
    let tcur = ark_mem.tcur;
    /* detach the aliased vectors for the call */
    let fpred = {
        let ARKodeMRIStepMem { Fse, Fsi, unify_Fs, .. } = &mut **step_mem;
        let fsi_arr = if *unify_Fs { Fse } else { Fsi };
        std::mem::take(&mut fsi_arr[smap])
    };
    let mut ypred = std::mem::take(&mut ark_mem.ycur);
    let mut t1 = std::mem::take(&mut ark_mem.tempv1);
    let mut t2 = std::mem::take(&mut ark_mem.tempv2);
    let mut t3 = std::mem::take(&mut ark_mem.tempv3);
    let mut jcur_local = step_mem.jcur;

    let retval = with_step_mem_installed(ark_mem, step_mem, |am| {
        lsetup(
            am,
            convfail,
            tcur,
            &mut ypred,
            &fpred,
            &mut jcur_local,
            &mut t1,
            &mut t2,
            &mut t3,
        )
    });

    /* reattach the vectors */
    ark_mem.ycur = ypred;
    ark_mem.tempv1 = t1;
    ark_mem.tempv2 = t2;
    ark_mem.tempv3 = t3;
    {
        let ARKodeMRIStepMem { Fse, Fsi, unify_Fs, .. } = &mut **step_mem;
        let fsi_arr = if *unify_Fs { Fse } else { Fsi };
        fsi_arr[smap] = fpred;
    }

    /* update Jacobian status */
    step_mem.jcur = jcur_local;
    *jcur = step_mem.jcur;

    /* update flags and 'gamma' values for last lsetup call */
    ark_mem.firststage = false;
    step_mem.gamrat = ONE;
    step_mem.crate_ = ONE;
    step_mem.gammap = step_mem.gamma;
    step_mem.nstlp = ark_mem.nst;

    if retval < 0 {
        return ARK_LSETUP_FAIL;
    }
    if retval > 0 {
        return CONV_FAIL;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_NlsLSolve:

  This routine wraps the ARKODE linear solver interface 'solve'
  routine for use by the nonlinear solver object.
  ---------------------------------------------------------------*/
fn mriStep_NlsLSolve(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeMRIStepMem>,
    nonlin_iter: i32,
    b: &mut NVector,
) -> i32 {
    let lsolve = step_mem.lsolve.unwrap();
    let smap = step_mem.stage_map[step_mem.istage as usize] as usize;
    let ernrm = step_mem.eRNrm;
    let tcur = ark_mem.tcur;
    /* detach the aliased vectors for the call */
    let fcur = {
        let ARKodeMRIStepMem { Fse, Fsi, unify_Fs, .. } = &mut **step_mem;
        let fsi_arr = if *unify_Fs { Fse } else { Fsi };
        std::mem::take(&mut fsi_arr[smap])
    };
    let ycur = std::mem::take(&mut ark_mem.ycur);

    /* call linear solver interface, and handle return value */
    let retval = with_step_mem_installed(ark_mem, step_mem, |am| {
        lsolve(am, b, tcur, &ycur, &fcur, ernrm, nonlin_iter)
    });

    ark_mem.ycur = ycur;
    {
        let ARKodeMRIStepMem { Fse, Fsi, unify_Fs, .. } = &mut **step_mem;
        let fsi_arr = if *unify_Fs { Fse } else { Fsi };
        fsi_arr[smap] = fcur;
    }

    if retval < 0 {
        return ARK_LSOLVE_FAIL;
    }
    if retval > 0 {
        return CONV_FAIL;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_NlsResidual:

  This routine evaluates the nonlinear residual for this
  solve-decoupled implicit MRI stage:
     z = zp + zc (stored in ark_mem->ycur)
     Fsi(z) (stored step_mem->Fsi[step_mem->stage_map[istage]])
     r = zc - gamma*Fsi(z) - step_mem->sdata
  ---------------------------------------------------------------*/
fn mriStep_NlsResidual(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeMRIStepMem>,
    r: &mut NVector,
) -> i32 {
    let smap = step_mem.stage_map[step_mem.istage as usize] as usize;

    /* update 'ycur' value as stored predictor + current corrector */
    {
        let ARKodeMRIStepMem { zpred, zcor, .. } = &mut **step_mem;
        N_VLinearSum(ONE, zpred, ONE, zcor, &mut ark_mem.ycur);
    }

    /* call the user-supplied pre-RHS function (if supplied), then call RHS */
    if let Some(pre_rhs) = ark_mem.PreRhsFn {
        let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
        let retval = pre_rhs(*tcur, ycur, user_data);
        if retval != 0 {
            return ARK_PRERHSFN_FAIL;
        }
    }
    let nls_fsi = step_mem.nls_fsi.unwrap();
    let retval = {
        let ARKodeMRIStepMem { Fse, Fsi, unify_Fs, .. } = &mut **step_mem;
        let fsi_arr = if *unify_Fs { Fse } else { Fsi };
        let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
        nls_fsi(*tcur, ycur, &mut fsi_arr[smap], user_data)
    };
    step_mem.nfsi += 1;
    if retval < 0 {
        return ARK_RHSFUNC_FAIL;
    }
    if retval > 0 {
        return RHSFUNC_RECVR;
    }

    /* compute residual: zcor - gamma*Fsi - sdata */
    {
        let gamma = step_mem.gamma;
        let ARKodeMRIStepMem { Fse, Fsi, unify_Fs, zcor, sdata, .. } = &mut **step_mem;
        let fsi_arr = if *unify_Fs { Fse } else { Fsi };
        let c = [ONE, -ONE, -gamma];
        let x: [&NVector; 3] = [zcor, sdata, &fsi_arr[smap]];
        N_VLinearCombination(3, &c, &x, r);
    }
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_NlsFPFunction:

  This routine evaluates the fixed point iteration function for
  this solve-decoupled implicit MRI stage:
     z = zp + zc (stored in ark_mem->ycur)
     Fsi(z) (store in step_mem->Fsi[step_mem->istage])
     g = gamma*Fsi(z) + step_mem->sdata
  ---------------------------------------------------------------*/
fn mriStep_NlsFPFunction(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeMRIStepMem>,
    g: &mut NVector,
) -> i32 {
    let smap = step_mem.stage_map[step_mem.istage as usize] as usize;

    /* update 'ycur' value as stored predictor + current corrector */
    {
        let ARKodeMRIStepMem { zpred, zcor, .. } = &mut **step_mem;
        N_VLinearSum(ONE, zpred, ONE, zcor, &mut ark_mem.ycur);
    }

    /* call the user-supplied pre-RHS function (if supplied), then call RHS */
    if let Some(pre_rhs) = ark_mem.PreRhsFn {
        let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
        let retval = pre_rhs(*tcur, ycur, user_data);
        if retval != 0 {
            return ARK_PRERHSFN_FAIL;
        }
    }
    let nls_fsi = step_mem.nls_fsi.unwrap();
    let retval = {
        let ARKodeMRIStepMem { Fse, Fsi, unify_Fs, .. } = &mut **step_mem;
        let fsi_arr = if *unify_Fs { Fse } else { Fsi };
        let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
        nls_fsi(*tcur, ycur, &mut fsi_arr[smap], user_data)
    };
    step_mem.nfsi += 1;
    if retval < 0 {
        return ARK_RHSFUNC_FAIL;
    }
    if retval > 0 {
        return RHSFUNC_RECVR;
    }

    /* combine parts:  g = gamma*Fsi(z) + sdata */
    {
        let gamma = step_mem.gamma;
        let ARKodeMRIStepMem { Fse, Fsi, unify_Fs, sdata, .. } = &mut **step_mem;
        let fsi_arr = if *unify_Fs { Fse } else { Fsi };
        N_VLinearSum(gamma, &fsi_arr[smap], ONE, sdata, g);
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  mriStep_NlsConvTest:

  This routine provides the nonlinear solver convergence test for
  this solve-decoupled implicit MRI stage.
  ---------------------------------------------------------------*/
fn mriStep_NlsConvTest(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeMRIStepMem>,
    del: &NVector,
    tol: f64,
    m: i32,
) -> i32 {
    /* if the problem is linearly implicit, just return success */
    if step_mem.linear {
        return SUN_SUCCESS;
    }

    /* compute the norm of the correction */
    let delnrm = N_VWrmsNorm(del, &ark_mem.ewt);

    /* update the stored estimate of the convergence rate (assumes linear
       convergence) */
    if m > 0 {
        step_mem.crate_ = SUNMAX(step_mem.crdown * step_mem.crate_, delnrm / step_mem.delp);
    }

    /* compute our scaled error norm for testing convergence */
    let dcon = SUNMIN(step_mem.crate_, ONE) * delnrm / tol;

    /* check for convergence; if so return with success */
    if dcon <= ONE {
        return SUN_SUCCESS;
    }

    /* check for divergence */
    if m >= 1 && delnrm > step_mem.rdiv * step_mem.delp {
        return SUN_NLS_CONV_RECVR;
    }

    /* save norm of correction for next iteration */
    step_mem.delp = delnrm;

    /* return with flag that there is more work to do */
    SUN_NLS_CONTINUE
}
