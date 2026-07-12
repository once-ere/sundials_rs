/* -----------------------------------------------------------------
 * Translated from src/arkode/arkode_arkstep_nls.c (ARKODE 7.7.0).
 * ARKStep nonlinear-solver interface.
 *
 * The generic C Newton iteration (SUNNonlinSolSolve_Newton) is
 * inlined here with the ARKStep Sys/LSetup/LSolve/CTest callbacks,
 * exactly mirroring the C control flow (donor: cvode_nls.rs).
 * C installs the callbacks into the NLS object via SetSysFn etc.;
 * here the residual variant is selected by the same flags C uses
 * when (re)installing (mass_type / predictor / autonomous), so the
 * dispatch is evaluated at the call, with identical results.
 *
 *
 * lsetup/lsolve op re-entries need step_mem installed in ark_mem
 * (arkLsSetup/arkLsSolve call step_getgammas etc.), so the wrappers
 * temporarily swap step_mem back in around the call (Addendum C.2).
 * -----------------------------------------------------------------*/
use crate::arkode_arkstep::arkStep_AccessStepMem;
use crate::arkode_arkstep_impl::*;
use crate::arkode_impl::*;
use crate::nvector_serial::*;
use crate::sundials_errors::SUN_SUCCESS;
use crate::sundials_math::*;
use crate::sundials_nonlinearsolver::*;
use crate::sundials_types::*;
use crate::sunnonlinsol_fixedpoint::FixedPointSolver;
use crate::sunnonlinsol_newton::NewtonSolver;

/*===============================================================
  Exported functions
  ===============================================================*/

/*---------------------------------------------------------------
  arkStep_SetNonlinearSolver:

  This routine attaches a SUNNonlinearSolver object to the
  ARKStep module.
  ---------------------------------------------------------------*/
pub fn arkStep_SetNonlinearSolver(ark_mem: &mut ARKodeMem, NLS: NonlinearSolver) -> i32 {
    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_SetNonlinearSolver") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* (C checks the NLS ops table for gettype/solve/setsysfn; the Rust
       NonlinearSolver enum always provides them) */

    /* free any existing nonlinear solver; set SUNNonlinearSolver pointer */
    step_mem.NLS = Some(NLS);
    step_mem.ownNLS = false;

    /* (default convergence test function arkStep_NlsConvTest is invoked
       directly by the inlined solve loop) */

    /* set default nonlinear iterations */
    if let Some(nls) = step_mem.NLS.as_mut() {
        let retval = nls.set_max_iters(step_mem.maxcor);
        if retval != 0 {
            ark_mem.step_mem = Some(step_mem);
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "arkStep_SetNonlinearSolver",
                file!(),
                "Setting maximum number of nonlinear iterations failed",
            );
            return ARK_ILL_INPUT;
        }
    }

    /* set the nonlinear system RHS function */
    if step_mem.fi.is_none() {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "arkStep_SetNonlinearSolver",
            file!(),
            "The implicit ODE RHS function is NULL",
        );
        return ARK_ILL_INPUT;
    }
    step_mem.nls_fi = step_mem.fi;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_SetNlsRhsFn:

  This routine sets an alternative user-supplied implicit ODE
  right-hand side function to use in the evaluation of nonlinear
  system functions.
  ---------------------------------------------------------------*/
pub fn arkStep_SetNlsRhsFn(ark_mem: &mut ARKodeMem, nls_fi: Option<ARKRhsFn>) -> i32 {
    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_SetNlsRhsFn") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    if nls_fi.is_some() {
        step_mem.nls_fi = nls_fi;
    } else {
        step_mem.nls_fi = step_mem.fi;
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_SetNlsSysFn:

  This routine sets the appropriate version of the nonlinear
  system function based on the current settings.  In C this
  installs a function pointer into the NLS object; the Rust solve
  loop derives the same choice from the flags at the call, so this
  routine performs the C validation only.
  ---------------------------------------------------------------*/
pub fn arkStep_SetNlsSysFn(ark_mem: &mut ARKodeMem) -> i32 {
    let step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_SetNlsSysFn") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* validate the mass matrix type as C does when selecting the
       residual/fixed-point function */
    if step_mem.mass_type != MASS_IDENTITY
        && step_mem.mass_type != MASS_FIXED
        && step_mem.mass_type != MASS_TIMEDEP
    {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "arkStep_SetNlsSysFn",
            file!(),
            "Invalid mass matrix type",
        );
        return ARK_ILL_INPUT;
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_GetNonlinearSystemData:

  This routine provides access to the relevant data needed to
  compute the nonlinear system function (Rust: scalar outputs and
  clones of the data vectors; C hands out raw pointers).
  ---------------------------------------------------------------*/
#[allow(clippy::too_many_arguments)]
pub fn arkStep_GetNonlinearSystemData(
    ark_mem: &mut ARKodeMem,
    tcur: &mut f64,
    zpred: &mut NVector,
    z: &mut NVector,
    Fi: &mut NVector,
    gamma: &mut f64,
    sdata: &mut NVector,
) -> i32 {
    let step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_GetNonlinearSystemData") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    *tcur = ark_mem.tcur;
    *zpred = step_mem.zpred.clone();
    *z = ark_mem.ycur.clone();
    *Fi = step_mem.Fi[step_mem.istage as usize].clone();
    *gamma = step_mem.gamma;
    *sdata = step_mem.sdata.clone();

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*===============================================================
  Utility routines called by ARKStep
  ===============================================================*/

/*---------------------------------------------------------------
  arkStep_NlsInit:

  This routine attaches the linear solver 'setup' and 'solve'
  routines to the nonlinear solver object (implicitly, via the
  inlined solve loop), and then initializes the nonlinear solver
  object itself.  This should only be called at the start of a
  simulation, after a re-init, or after a re-size.
  ---------------------------------------------------------------*/
pub fn arkStep_NlsInit(ark_mem: &mut ARKodeMem) -> i32 {
    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_NlsInit") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* reset counters */
    step_mem.nls_iters = 0;
    step_mem.nls_fails = 0;

    /* (the lsetup/lsolve wrapper functions are dispatched by the solve
       loop on step_mem.lsetup / step_mem.lsolve presence) */

    ark_mem.step_mem = Some(step_mem);

    let retval = arkStep_SetNlsSysFn(ark_mem);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "arkStep_NlsInit",
            file!(),
            "Setting nonlinear system function failed",
        );
        return ARK_ILL_INPUT;
    }

    /* initialize nonlinear solver */
    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_NlsInit") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };
    let retval = step_mem.NLS.as_mut().map(|n| n.initialize()).unwrap_or(0);
    ark_mem.step_mem = Some(step_mem);
    if retval != 0 {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "arkStep_NlsInit",
            file!(),
            MSG_NLS_INIT_FAIL,
        );
        return ARK_NLS_INIT_FAIL;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_Nls

  This routine attempts to solve the nonlinear system associated
  with a single implicit stage.  It calls the supplied
  SUNNonlinearSolver object (here: the inlined Newton loop) to
  perform the solve.

  Upon entry, the predicted solution is held in step_mem->zpred,
  which is never changed throughout this routine.  If an initial
  attempt at solving the nonlinear system fails (e.g. due to a
  stale Jacobian), this allows for new attempts at the solution.

  Upon a successful solve, the solution is held in ark_mem->ycur.
  ---------------------------------------------------------------*/
pub fn arkStep_Nls(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeARKStepMem>,
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
       Newton solver out of step_mem for the iteration */
    let nls = match step_mem.NLS.take() {
        Some(n) => n,
        None => return ARK_MEM_NULL,
    };
    let (retval, nls) = match nls {
        NonlinearSolver::Newton(mut ns) => {
            let tol = step_mem.nlscoef;
            let ret = arkStep_NlsSolveNewton(ark_mem, step_mem, &mut ns, tol, call_lsetup);
            (ret, NonlinearSolver::Newton(ns))
        }
        NonlinearSolver::FixedPoint(mut fps) => {
            let tol = step_mem.nlscoef;
            let ret = arkStep_NlsSolveFixedPoint(ark_mem, step_mem, &mut fps, tol);
            (ret, NonlinearSolver::FixedPoint(fps))
        }
    };

    /* increment counters */
    let (iters_inc, fails_inc) = (nls.get_num_iters(), nls.get_num_conv_fails());
    step_mem.NLS = Some(nls);
    step_mem.nls_iters += iters_inc;
    step_mem.nls_fails += fails_inc;

    /* successful solve -- reset jcur flag and apply correction */
    if retval == SUN_SUCCESS {
        step_mem.jcur = false;
        let ARKodeARKStepMem { zcor, zpred, .. } = &mut **step_mem;
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
  arkStep_NlsSolveNewton: the generic SUNNonlinSolSolve_Newton
  loop (sunnonlinsol_newton.c), specialized to the ARKStep
  callbacks.  ycor = step_mem->zcor, w = ark_mem->ewt.
  ---------------------------------------------------------------*/
fn arkStep_NlsSolveNewton(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeARKStepMem>,
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
        retval = arkStep_NlsResidual(ark_mem, step_mem, ns.curiter, &mut ns.delta);
        if retval != 0 {
            break 'outer;
        }

        /* if indicated, setup the linear system */
        if call_lsetup {
            let mut jcur = ns.jcur;
            retval = arkStep_NlsLSetup(ark_mem, step_mem, jbad, &mut jcur);
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
            retval = arkStep_NlsLSolve(ark_mem, step_mem, ns.curiter, &mut ns.delta);
            if retval != 0 {
                break;
            }

            /* update the Newton iterate */
            step_mem.zcor.linear_sum_with(ONE, ONE, &ns.delta);

            /* test for convergence */
            retval = arkStep_NlsConvTest(ark_mem, step_mem, &ns.delta, tol, ns.curiter);

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
            retval = arkStep_NlsResidual(ark_mem, step_mem, ns.curiter, &mut ns.delta);
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

/*===============================================================
  Interface routines supplied to the SUNNonlinearSolver module
  ===============================================================*/

/* helper: temporarily re-install step_mem into ark_mem around an
   ARKLS lsetup/lsolve call (those re-enter the step_getgammas /
   step_getimplicitrhs / step_setjcur ops, which need step_mem in
   place — Addendum C.2). */
fn with_step_mem_installed<R>(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeARKStepMem>,
    f: impl FnOnce(&mut ARKodeMem) -> R,
) -> R {
    let owned = std::mem::replace(step_mem, Box::new(ARKodeARKStepMem::default()));
    ark_mem.step_mem = Some(owned);
    let r = f(ark_mem);
    let owned = ark_mem
        .step_mem
        .take()
        .unwrap()
        .downcast::<ARKodeARKStepMem>()
        .unwrap();
    *step_mem = owned;
    r
}

/*---------------------------------------------------------------
  arkStep_NlsLSetup:

  This routine wraps the ARKODE linear solver interface 'setup'
  routine for use by the nonlinear solver object.
  ---------------------------------------------------------------*/
fn arkStep_NlsLSetup(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeARKStepMem>,
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
    let istage = step_mem.istage as usize;
    let tcur = ark_mem.tcur;
    /* detach the aliased vectors for the call (C passes pointers into
       ark_mem/step_mem storage) */
    let fpred = std::mem::take(&mut step_mem.Fi[istage]);
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
    step_mem.Fi[istage] = fpred;

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
  arkStep_NlsLSolve:

  This routine wraps the ARKODE linear solver interface 'solve'
  routine for use by the nonlinear solver object.
  ---------------------------------------------------------------*/
fn arkStep_NlsLSolve(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeARKStepMem>,
    nonlin_iter: i32,
    b: &mut NVector,
) -> i32 {
    let lsolve = step_mem.lsolve.unwrap();
    let istage = step_mem.istage as usize;
    let ernrm = step_mem.eRNrm;
    let tcur = ark_mem.tcur;
    /* detach the aliased vectors for the call */
    let fcur = std::mem::take(&mut step_mem.Fi[istage]);
    let ycur = std::mem::take(&mut ark_mem.ycur);

    /* call linear solver interface, and handle return value */
    let retval = with_step_mem_installed(ark_mem, step_mem, |am| {
        lsolve(am, b, tcur, &ycur, &fcur, ernrm, nonlin_iter)
    });

    ark_mem.ycur = ycur;
    step_mem.Fi[istage] = fcur;

    if retval < 0 {
        return ARK_LSOLVE_FAIL;
    }
    if retval > 0 {
        return CONV_FAIL;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_NlsResidual (MassIdent and the TrivialPredAutonomous
  variant, selected by the same flags C uses in SetNlsSysFn)

  This routine evaluates the nonlinear residual for the additive
  Runge-Kutta method.  It assumes that any data from previous
  time steps/stages is contained in step_mem, and merely combines
  this old data with the current implicit ODE RHS vector to
  compute the nonlinear residual r:
     z = zp + zc (stored in ark_mem->ycur)
     Fi(z) (stored step_mem->Fi[step_mem->istage])
     r = zc - gamma*Fi(z) - step_mem->sdata

  The "TrivialPredAutonomous" version reuses the implicit RHS
  evaluation at the beginning of the step in the initial residual
  evaluation.
  ---------------------------------------------------------------*/
fn arkStep_NlsResidual(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeARKStepMem>,
    nls_iter: i32,
    r: &mut NVector,
) -> i32 {
    let istage = step_mem.istage as usize;

    /* update 'ycur' value as stored predictor + current corrector */
    {
        let ARKodeARKStepMem { zpred, zcor, .. } = &mut **step_mem;
        N_VLinearSum(ONE, zpred, ONE, zcor, &mut ark_mem.ycur);
    }

    /* MassTDep variant: r = M(t)*(zcor - sdata) - gamma*Fi(z) */
    if step_mem.mass_type == MASS_TIMEDEP {
        /* put M*(zcor - sdata) in r (use Fi[is] as temporary storage) */
        {
            let ARKodeARKStepMem { zcor, sdata, Fi, .. } = &mut **step_mem;
            N_VLinearSum(ONE, zcor, -ONE, sdata, &mut Fi[istage]);
        }
        let mmult = step_mem.mmult.unwrap();
        let retval = mmult(ark_mem, &step_mem.Fi[istage], r);
        if retval != ARK_SUCCESS {
            return ARK_MASSMULT_FAIL;
        }

        /* call the user-supplied pre-RHS function (if supplied), then call RHS */
        if let Some(pre_rhs) = ark_mem.PreRhsFn {
            let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
            let retval = pre_rhs(*tcur, ycur, user_data);
            if retval != 0 {
                return ARK_PRERHSFN_FAIL;
            }
        }
        let nls_fi = step_mem.nls_fi.unwrap();
        let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
        let retval = nls_fi(*tcur, ycur, &mut step_mem.Fi[istage], user_data);
        step_mem.nfi += 1;
        if retval < 0 {
            return ARK_RHSFUNC_FAIL;
        }
        if retval > 0 {
            return RHSFUNC_RECVR;
        }

        /* compute residual via linear sum */
        r.linear_sum_with(ONE, -step_mem.gamma, &step_mem.Fi[istage]);
        return ARK_SUCCESS;
    }

    /* TrivialPredAutonomous variant (MassIdent / MassFixed): reuse the
       saved implicit RHS evaluation on the first iteration */
    let tpa = step_mem.predictor == 0
        && step_mem.autonomous
        && step_mem.fn_implicit != FnImplicitAlias::None;
    if tpa && nls_iter == 0 {
        match step_mem.fn_implicit {
            FnImplicitAlias::Fi0 => {
                let (head, tail) = step_mem.Fi.split_at_mut(istage);
                tail[0].data.copy_from_slice(&head[0].data);
            }
            FnImplicitAlias::Tempv5 => {
                step_mem.Fi[istage].data.copy_from_slice(&ark_mem.tempv5.data);
            }
            FnImplicitAlias::ArkFn => {
                step_mem.Fi[istage].data.copy_from_slice(&ark_mem.fn_.data);
            }
            FnImplicitAlias::None => unreachable!(),
        }
    } else {
        /* call the user-supplied pre-RHS function (if supplied), then call RHS */
        if let Some(pre_rhs) = ark_mem.PreRhsFn {
            let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
            let retval = pre_rhs(*tcur, ycur, user_data);
            if retval != 0 {
                return ARK_PRERHSFN_FAIL;
            }
        }
        let nls_fi = step_mem.nls_fi.unwrap();
        let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
        let retval = nls_fi(*tcur, ycur, &mut step_mem.Fi[istage], user_data);
        step_mem.nfi += 1;
        if retval < 0 {
            return ARK_RHSFUNC_FAIL;
        }
        if retval > 0 {
            return RHSFUNC_RECVR;
        }
    }

    if step_mem.mass_type == MASS_FIXED {
        /* put M*zcor in r, then r = r - sdata - gamma*Fi */
        let mmult = step_mem.mmult.unwrap();
        let retval = mmult(ark_mem, &step_mem.zcor, r);
        if retval != ARK_SUCCESS {
            return ARK_MASSMULT_FAIL;
        }
        /* z == X[0] with c0 == 1: in-place accumulate form of the
           3-term N_VLinearCombination */
        let ARKodeARKStepMem { sdata, Fi, gamma, .. } = &mut **step_mem;
        for e in 0..r.data.len() {
            r.data[e] += -ONE * sdata.data[e];
        }
        for e in 0..r.data.len() {
            r.data[e] += -*gamma * Fi[istage].data[e];
        }
    } else {
        /* compute residual via linear combination */
        let c = [ONE, -ONE, -step_mem.gamma];
        let ARKodeARKStepMem { zcor, sdata, Fi, .. } = &mut **step_mem;
        let x: [&NVector; 3] = [zcor, sdata, &Fi[istage]];
        N_VLinearCombination(3, &c, &x, r);
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_NlsConvTest:

  This routine provides the nonlinear solver convergence test for
  the additive Runge-Kutta method.
  ---------------------------------------------------------------*/
fn arkStep_NlsConvTest(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeARKStepMem>,
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

/*---------------------------------------------------------------
  arkStep_NlsSolveFixedPoint: the generic
  SUNNonlinSolSolve_FixedPoint loop (sunnonlinsol_fixedpoint.c),
  specialized to the ARKStep callbacks (donor cvode_nls pattern).
  ycor = step_mem->zcor, w = ark_mem->ewt.
  ---------------------------------------------------------------*/
fn arkStep_NlsSolveFixedPoint(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeARKStepMem>,
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
        /* update previous solution guess */
        fps.yprev.data.copy_from_slice(&step_mem.zcor.data);

        /* compute fixed-point iteration function, store in gy */
        {
            let mut gy = std::mem::take(&mut fps.gy);
            let retval = arkStep_NlsFPFunction(ark_mem, step_mem, fps.curiter, &mut gy);
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
            let ARKodeARKStepMem { zcor, .. } = &mut **step_mem;
            N_VLinearSum(ONE, zcor, -ONE, &fps.yprev, &mut fps.delta);
        }

        /* test for convergence */
        let retval = {
            let delta = std::mem::take(&mut fps.delta);
            let r = arkStep_NlsConvTest(ark_mem, step_mem, &delta, tol, fps.curiter);
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

/*---------------------------------------------------------------
  arkStep_NlsFPFunction (MassIdent and the TrivialPredAutonomous
  variant, selected by the same flags C uses in SetNlsSysFn)

  This routine evaluates the fixed point iteration function for
  the additive Runge-Kutta method (identity mass matrix):
     Fi(z) (store in step_mem->Fi[step_mem->istage])
     g = gamma*Fi(z) + step_mem->sdata

  The "TrivialPredAutonomous" version reuses the implicit RHS
  evaluation at the beginning of the step in the initial FP
  function evaluation.
  ---------------------------------------------------------------*/
fn arkStep_NlsFPFunction(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeARKStepMem>,
    nls_iter: i32,
    g: &mut NVector,
) -> i32 {
    let istage = step_mem.istage as usize;

    /* update 'ycur' value as stored predictor + current corrector */
    {
        let ARKodeARKStepMem { zpred, zcor, .. } = &mut **step_mem;
        N_VLinearSum(ONE, zpred, ONE, zcor, &mut ark_mem.ycur);
    }

    /* TrivialPredAutonomous variant (MassIdent / MassFixed; the
       MassTDep fixed-point function has no TPA form): reuse the saved
       implicit RHS evaluation on the first iteration */
    let tpa = step_mem.predictor == 0
        && step_mem.autonomous
        && step_mem.mass_type != MASS_TIMEDEP
        && step_mem.fn_implicit != FnImplicitAlias::None;
    if tpa && nls_iter == 0 {
        match step_mem.fn_implicit {
            FnImplicitAlias::Fi0 => {
                let (head, tail) = step_mem.Fi.split_at_mut(istage);
                tail[0].data.copy_from_slice(&head[0].data);
            }
            FnImplicitAlias::Tempv5 => {
                step_mem.Fi[istage].data.copy_from_slice(&ark_mem.tempv5.data);
            }
            FnImplicitAlias::ArkFn => {
                step_mem.Fi[istage].data.copy_from_slice(&ark_mem.fn_.data);
            }
            FnImplicitAlias::None => unreachable!(),
        }
    } else {
        /* call the user-supplied pre-RHS function (if supplied), then call RHS */
        if let Some(pre_rhs) = ark_mem.PreRhsFn {
            let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
            let retval = pre_rhs(*tcur, ycur, user_data);
            if retval != 0 {
                return ARK_PRERHSFN_FAIL;
            }
        }
        let nls_fi = step_mem.nls_fi.unwrap();
        let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
        let retval = nls_fi(*tcur, ycur, &mut step_mem.Fi[istage], user_data);
        step_mem.nfi += 1;
        if retval < 0 {
            return ARK_RHSFUNC_FAIL;
        }
        if retval > 0 {
            return RHSFUNC_RECVR;
        }
    }

    if step_mem.mass_type == MASS_TIMEDEP {
        /* copy gamma*Fi into g, perform mass matrix solve, then
           combine parts:  g = g + sdata */
        {
            let ARKodeARKStepMem { Fi, gamma, .. } = &mut **step_mem;
            N_VScale(*gamma, &Fi[istage], g);
        }
        let msolve = step_mem.msolve.unwrap();
        let tol = step_mem.nlscoef;
        let retval = msolve(ark_mem, g, tol);
        if retval < 0 {
            return ARK_RHSFUNC_FAIL;
        }
        if retval > 0 {
            return RHSFUNC_RECVR;
        }
        g.linear_sum_with(ONE, ONE, &step_mem.sdata);
    } else {
        /* combine parts:  g = gamma*Fi(z) + sdata */
        {
            let ARKodeARKStepMem { sdata, Fi, gamma, .. } = &mut **step_mem;
            N_VLinearSum(*gamma, &Fi[istage], ONE, sdata, g);
        }

        /* perform mass matrix solve (fixed mass matrix) */
        if step_mem.mass_type == MASS_FIXED {
            let msolve = step_mem.msolve.unwrap();
            let tol = step_mem.nlscoef;
            let retval = msolve(ark_mem, g, tol);
            if retval < 0 {
                return ARK_RHSFUNC_FAIL;
            }
            if retval > 0 {
                return RHSFUNC_RECVR;
            }
        }
    }

    ARK_SUCCESS
}
