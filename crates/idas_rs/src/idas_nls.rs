/* -----------------------------------------------------------------
 * Translated from src/idas/idas_nls.c (IDAS 7.7.0), together with
 * idas.c's IDANls driver and the Newton solve loop of
 * sunnonlinsol_newton.c that the C code reaches through the
 * SUNNonlinearSolver ops table.  The IDAS-specific callbacks
 * (idaNlsResidual, idaNlsLSetup, idaNlsLSolve, idaNlsConvTest) are
 * collapsed donor-style (ida_nls.rs) into direct functions over
 * &mut IDAMem; control flow and arithmetic order match the C
 * sources statement for statement.
 *
 * IDAS deltas vs the verified ida_nls.rs donor:
 *   - idaNlsLSetup additionally clears ida_forceSetup and resets
 *     ida_ssS = TWENTY (idas_nls.c lines 253/266);
 *   - IDANls (idas.c) adds the IDA_SIMULTANEOUS corrector dispatch:
 *     when sensi_sim the composite system [ee, eeS] is solved through
 *     NLSsim (idas_nls_sim.rs::idaNlsSolveSensSim), the ida_forceSetup
 *     flag participates in the lsetup decision, ida_ssS mirrors the
 *     ida_ss resets, and the final correction also updates yyS/ypS.
 *
 * Exported entry points for idas.rs:
 *   idaNlsInit(ida_mem: &mut IDAMem) -> i32      (C idaNlsInit)
 *   IDANls(ida_mem: &mut IDAMem) -> i32          (C idas.c IDANls,
 *       the nonlinear solve dispatched from IDAStep)
 * -----------------------------------------------------------------*/
use crate::idas_impl::*;
use crate::idas_ls::{idaLsSetup, idaLsSolve};
use crate::nvector_serial::*;
use crate::sundials_errors::SUN_SUCCESS;
use crate::sundials_math::*;
use crate::sundials_nonlinearsolver::*;
use crate::sundials_types::*;
use crate::sunnonlinsol_newton::NewtonSolver;

/* constant macros */
const PT0001: f64 = 0.0001; /* real 0.0001 */
const ZERO: f64 = 0.0; /* real 0.0 (idas.c) */
const ONE: f64 = 1.0; /* real 1.0    */
const TWENTY: f64 = 20.0; /* real 20.0   */
const HUNDRED: f64 = 100.0; /* real 100.0 (idas.c) */

/* nonlinear solver parameters */
const MAXIT: i32 = 4; /* default max number of nonlinear iterations    */
const RATEMAX: f64 = 0.9; /* max convergence rate used in divergence check */

/* -----------------------------------------------------------------------------
 * Exported functions
 * ---------------------------------------------------------------------------*/

pub fn IDASetNonlinearSolver(ida_mem: &mut IDAMem, NLS: NonlinearSolver) -> i32 {
    /* (The C NULL-input and missing-ops checks — gettype/solve/setsysfn
       on the NLS — cannot fail here: the workspace enum implements
       them.) */

    /* check for allowed nonlinear solver types */
    if NLS.nls_type() != SUNNONLINEARSOLVER_ROOTFIND {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASetNonlinearSolver", file!(),
                        "NLS type must be SUNNONLINEARSOLVER_ROOTFIND");
        return IDA_ILL_INPUT;
    }

    /* free any existing nonlinear solver: RAII (drop on overwrite) */

    /* set SUNNonlinearSolver pointer */
    ida_mem.NLS = Some(NLS);

    /* Set NLS ownership flag. If this function was called to attach the default
       NLS, IDA will set the flag to SUNTRUE after this function returns. */
    ida_mem.ownNLS = SUNFALSE;

    /* the nonlinear residual and convergence test functions
       (idaNlsResidual / idaNlsConvTest) are wired statically into the
       Newton solve driver below; the C SUNNonlinSolSetSysFn /
       SUNNonlinSolSetConvTestFn registrations have no Rust
       counterpart */

    /* set max allowed nonlinear iterations */
    let retval = ida_mem.NLS.as_mut().unwrap().set_max_iters(MAXIT);
    if retval != IDA_SUCCESS {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASetNonlinearSolver", file!(),
                        "Setting maximum number of nonlinear iterations failed");
        return IDA_ILL_INPUT;
    }

    /* Set the nonlinear system RES function */
    if ida_mem.ida_res.is_none() {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASetNonlinearSolver", file!(),
                        "The DAE residual function is NULL");
        return IDA_ILL_INPUT;
    }
    ida_mem.nls_res = ida_mem.ida_res;

    IDA_SUCCESS
}

/*---------------------------------------------------------------
  IDASetNlsResFn:

  This routine sets an alternative user-supplied DAE residual
  function to use in the evaluation of nonlinear system functions.
  ---------------------------------------------------------------*/
pub fn IDASetNlsResFn(ida_mem: &mut IDAMem, res: Option<IDAResFn>) -> i32 {
    match res {
        Some(f) => ida_mem.nls_res = Some(f),
        None => ida_mem.nls_res = ida_mem.ida_res,
    }

    IDA_SUCCESS
}

/*---------------------------------------------------------------
  IDAGetNonlinearSystemData:

  This routine provides access to the relevant data needed to
  compute the nonlinear system function (out-pointers become a
  returned tuple (tcur, cj), donor convention; the vectors —
  yypredict, yppredict, yy, yp, savres — and user_data remain
  accessible as IDAMem fields).
  ---------------------------------------------------------------*/
pub fn IDAGetNonlinearSystemData(ida_mem: &IDAMem) -> (f64, f64) {
    (ida_mem.ida_tn, ida_mem.ida_cj)
}

/* -----------------------------------------------------------------------------
 * Private functions
 * ---------------------------------------------------------------------------*/

/* `IDA_mem->ida_lsetup != NULL` guard: the lsetup hook exists when an
   IDALS module is attached and idaLsInitialize has not disabled it
   (idas_ls_impl.rs setup_disabled contract). */
pub(crate) fn ida_has_lsetup(ida_mem: &IDAMem) -> bool {
    match &ida_mem.ida_lmem {
        LsModule::Ls(idals_mem) => !idals_mem.setup_disabled,
        LsModule::None => false,
    }
}

pub fn idaNlsInit(ida_mem: &mut IDAMem) -> i32 {
    /* In C this wires the idaNlsLSetup/idaNlsLSolve wrapper functions
       into the NLS depending on whether ida_lsetup/ida_lsolve exist
       (infallible registrations); here the dispatch is dynamic through
       ida_has_lsetup / the LsModule enum. */

    /* initialize nonlinear solver */
    let retval = match ida_mem.NLS.as_mut() {
        Some(nls) => nls.initialize(),
        None => -1,
    };

    if retval != IDA_SUCCESS {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "idaNlsInit", file!(),
                        MSG_NLS_INIT_FAIL);
        return IDA_NLS_INIT_FAIL;
    }

    IDA_SUCCESS
}

/* idaNlsLSetup (idas_nls.c): wrapper around the lsetup dispatch.
   In C the lsetup call receives (yy, yp, savres, tempv1..3); the
   vectors alias IDAMem fields, so they are detached for the call
   (the three tmps come from the IDAMem fields inside idaLsSetup). */
fn idaNlsLSetup(ida_mem: &mut IDAMem, _jbad: bool, jcur: &mut bool) -> i32 {
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

/* idaNlsLSolve (idas_nls.c): wrapper around the lsolve dispatch.
   In C the lsolve call receives (delta, ewt, yy, yp, savres); the
   weight and current-state vectors alias IDAMem fields, so they are
   detached for the call. */
fn idaNlsLSolve(ida_mem: &mut IDAMem, delta: &mut NVector) -> i32 {
    let mut lmem = std::mem::take(&mut ida_mem.ida_lmem);
    let retval = match &mut lmem {
        LsModule::Ls(idals_mem) => {
            let ewt = std::mem::take(&mut ida_mem.ida_ewt);
            let yy = std::mem::take(&mut ida_mem.ida_yy);
            let yp = std::mem::take(&mut ida_mem.ida_yp);
            let savres = std::mem::take(&mut ida_mem.ida_savres);
            let r = idaLsSolve(ida_mem, idals_mem, delta, &ewt, &yy, &yp, &savres);
            ida_mem.ida_ewt = ewt;
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

/* idaNlsResidual (idas_nls.c): updates yy/yp from the current
   correction and evaluates the DAE residual; saves a copy in
   savres.  ycor = the Newton correction vector (ida_ee, detached
   by the caller). */
fn idaNlsResidual(ida_mem: &mut IDAMem, ycor: &NVector, res: &mut NVector) -> i32 {
    /* update yy and yp based on the current correction */
    {
        let IDAMem { ida_yypredict, ida_yy, .. } = ida_mem;
        N_VLinearSum(ONE, ida_yypredict, ONE, ycor, ida_yy);
    }
    {
        let IDAMem { ida_yppredict, ida_yp, ida_cj, .. } = ida_mem;
        N_VLinearSum(ONE, ida_yppredict, *ida_cj, ycor, ida_yp);
    }

    /* evaluate residual */
    let res_fn = ida_mem.nls_res.unwrap();
    let retval = {
        let IDAMem { ida_tn, ida_yy, ida_yp, ida_user_data, .. } = ida_mem;
        res_fn(*ida_tn, ida_yy, ida_yp, res, ida_user_data)
    };

    /* increment the number of residual evaluations */
    ida_mem.ida_nre += 1;

    /* save a copy of the residual vector in savres */
    N_VScale(ONE, res, &mut ida_mem.ida_savres);

    if retval < 0 {
        return IDA_RES_FAIL;
    }
    if retval > 0 {
        return IDA_RES_RECVR;
    }

    IDA_SUCCESS
}

/* idaNlsConvTest (idas_nls.c).  m is the current nonlinear iteration
   count (C SUNNonlinSolGetCurIter); ewt is detached from IDAMem by
   the caller; the C ycor input is unused. */
fn idaNlsConvTest(ida_mem: &mut IDAMem, del: &NVector, tol: f64, ewt: &NVector, m: i32) -> i32 {
    /* compute the norm of the correction */
    let delnrm = N_VWrmsNorm(del, ewt);

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
 * IDANls (idas.c)
 *
 * This routine attempts to solve the nonlinear system using the linear
 * solver specified.
 *
 *  Possible return values:
 *
 *  IDA_SUCCESS
 *
 *  IDA_RES_RECVR       IDA_RES_FAIL
 *  IDA_SRES_RECVR      IDA_SRES_FAIL
 *  IDA_LSETUP_RECVR    IDA_LSETUP_FAIL
 *  IDA_LSOLVE_RECVR    IDA_LSOLVE_FAIL
 *
 *  SUN_NLS_CONV_RECVR
 *  IDA_MEM_NULL
 */
pub fn IDANls(ida_mem: &mut IDAMem) -> i32 {
    let mut callLSetup = SUNFALSE;

    /* Are we computing sensitivities with the IDA_SIMULTANEOUS approach? */
    let sensi_sim = ida_mem.ida_sensi && ida_mem.ida_ism == IDA_SIMULTANEOUS;

    /* Initialize if the first time called */

    if ida_mem.ida_nst == 0 {
        ida_mem.ida_cjold = ida_mem.ida_cj;
        ida_mem.ida_ss = TWENTY;
        ida_mem.ida_ssS = TWENTY;
        if ida_has_lsetup(ida_mem) {
            callLSetup = SUNTRUE;
        }
    }

    /* Decide if lsetup is to be called */

    if ida_has_lsetup(ida_mem) {
        ida_mem.ida_cjratio = ida_mem.ida_cj / ida_mem.ida_cjold;
        let temp1 = (ONE - ida_mem.ida_dcj) / (ONE + ida_mem.ida_dcj);
        let temp2 = ONE / temp1;
        if ida_mem.ida_cjratio < temp1 || ida_mem.ida_cjratio > temp2 {
            callLSetup = SUNTRUE;
        }
        if ida_mem.ida_forceSetup {
            callLSetup = SUNTRUE;
        }
        if ida_mem.ida_cj != ida_mem.ida_cjlast {
            ida_mem.ida_ss = HUNDRED;
            ida_mem.ida_ssS = HUNDRED;
        }
    }

    /* solve the nonlinear system */
    let tol = ida_mem.ida_epsNewt;
    let retval;
    if sensi_sim {
        /* initial guess for the correction to the predictor: the C
           N_VConst(ZERO, ycorSim) zeroes the composite [ee, eeS]
           (senswrapper aliases are not stored — pinned convention) */
        N_VConst(ZERO, &mut ida_mem.ida_ee);
        for is in 0..ida_mem.ida_Ns as usize {
            N_VConst(ZERO, &mut ida_mem.ida_eeS[is]);
        }

        /* call nonlinear solver setup if it exists: the workspace Newton
           solver has no setup operation, so the C guard is a no-op */

        let mut nls = match ida_mem.NLSsim.take() {
            Some(nls) => nls,
            None => return IDA_MEM_NULL,
        };
        retval = crate::idas_nls_sim::idaNlsSolveSensSim(ida_mem, &mut nls, tol, callLSetup);

        /* increment counters */
        ida_mem.ida_nni += nls.get_num_iters();
        ida_mem.ida_nnf += nls.get_num_conv_fails();
        ida_mem.NLSsim = Some(nls);
    } else {
        /* initial guess for the correction to the predictor
           (ee is detached from IDAMem for the duration of the solve: it
           is the ycor vector the C code hands to SUNNonlinSolSolve) */
        let mut ee = std::mem::take(&mut ida_mem.ida_ee);
        N_VConst(ZERO, &mut ee);

        /* call nonlinear solver setup if it exists: the workspace Newton
           solver has no setup operation, so the C guard is a no-op */

        let mut nls = match ida_mem.NLS.take() {
            Some(nls) => nls,
            None => {
                ida_mem.ida_ee = ee;
                return IDA_MEM_NULL;
            }
        };
        retval = match &mut nls {
            NonlinearSolver::Newton(ns) => idaNlsSolveNewton(ida_mem, ns, &mut ee, tol, callLSetup),
            /* cannot arise: IDASetNonlinearSolver rejects non-ROOTFIND types */
            NonlinearSolver::FixedPoint(_) => IDA_NLS_FAIL,
        };

        /* increment counters */
        ida_mem.ida_nni += nls.get_num_iters();
        ida_mem.ida_nnf += nls.get_num_conv_fails();
        ida_mem.NLS = Some(nls);
        ida_mem.ida_ee = ee;
    }

    /* return if nonlinear solver failed */
    if retval != SUN_SUCCESS {
        return retval;
    }

    /* update yy and yp based on the final correction from the nonlinear solver */
    {
        let IDAMem { ida_yypredict, ida_ee, ida_yy, .. } = ida_mem;
        N_VLinearSum(ONE, ida_yypredict, ONE, ida_ee, ida_yy);
    }
    {
        let IDAMem { ida_yppredict, ida_ee, ida_yp, ida_cj, .. } = ida_mem;
        N_VLinearSum(ONE, ida_yppredict, *ida_cj, ida_ee, ida_yp);
    }

    /* update the sensitivities based on the final correction from the nonlinear
       solver (N_VLinearSumVectorArray expanded to per-vector loops) */
    if sensi_sim {
        let ns = ida_mem.ida_Ns as usize;
        {
            let IDAMem { ida_yySpredict, ida_eeS, ida_yyS, .. } = ida_mem;
            for is in 0..ns {
                N_VLinearSum(ONE, &ida_yySpredict[is], ONE, &ida_eeS[is], &mut ida_yyS[is]);
            }
        }
        {
            let IDAMem { ida_ypSpredict, ida_eeS, ida_ypS, ida_cj, .. } = ida_mem;
            for is in 0..ns {
                N_VLinearSum(ONE, &ida_ypSpredict[is], *ida_cj, &ida_eeS[is], &mut ida_ypS[is]);
            }
        }
    }

    IDA_SUCCESS
}

/*
 * SUNNonlinSolSolve_Newton (sunnonlinsol_newton.c), specialized to
 * the IDAS callbacks.  ycor = ida_ee (detached by IDANls),
 * w = ida_ewt (detached transiently around the convergence test).
 */
fn idaNlsSolveNewton(
    ida_mem: &mut IDAMem,
    ns: &mut NewtonSolver,
    ycor: &mut NVector,
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
        retval = idaNlsResidual(ida_mem, ycor, &mut ns.delta);
        if retval != 0 {
            break 'outer;
        }

        /* if indicated, setup the linear system */
        if call_lsetup {
            let mut jcur = ns.jcur;
            retval = idaNlsLSetup(ida_mem, jbad, &mut jcur);
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
            ns.delta.scale_inplace(-ONE);

            /* solve the linear system to get Newton update delta */
            retval = idaNlsLSolve(ida_mem, &mut ns.delta);
            if retval != 0 {
                break;
            }

            /* update the Newton iterate */
            ycor.linear_sum_with(ONE, ONE, &ns.delta);

            /* test for convergence */
            retval = {
                let ewt = std::mem::take(&mut ida_mem.ida_ewt);
                let r = idaNlsConvTest(ida_mem, &ns.delta, tol, &ewt, ns.curiter);
                ida_mem.ida_ewt = ewt;
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
            retval = idaNlsResidual(ida_mem, ycor, &mut ns.delta);
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
            N_VConst(ZERO, ycor);
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
    use crate::idas_ls::{idaLs_AccessLMem, IDASetLinearSolver};
    use crate::sundials_context::SUNContext;
    use crate::sundials_linearsolver::LinearSolver;
    use crate::sundials_matrix::SUNMatrix;
    use crate::sunlinsol_dense::SUNLinSol_Dense;
    use crate::sunmatrix_dense::SUNDenseMatrix;
    use crate::sunnonlinsol_newton::SUNNonlinSol_Newton;

    /* decoupled linear DAE: F_i = yp_i + 2*y_i - 4, so the Newton
       matrix is J = dF/dy + cj*dF/dyp = (2 + cj) I and the nonlinear
       system F(ypredict + e, yppredict + cj*e) = 0 has the exact
       solution e_i = (4 - yppredict_i - 2*yypredict_i)/(2 + cj). */
    fn resfn(_t: f64, yy: &NVector, yp: &NVector, rr: &mut NVector, _ud: &mut UserData) -> i32 {
        for i in 0..yy.len() {
            rr.data[i] = yp.data[i] + 2.0 * yy.data[i] - 4.0;
        }
        0
    }

    fn make_ida_mem(n: usize) -> IDAMem {
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
        ida_mem
    }

    /* IDASetNonlinearSolver stores the solver, clears ownership, sets
       maxiters = MAXIT = 4 and wires nls_res (idas_nls.c lines 44-132);
       a fixed-point (non-ROOTFIND) solver is rejected. */
    #[test]
    fn idasetnonlinearsolver_defaults() {
        let sunctx = SUNContext::default();
        let mut ida_mem = make_ida_mem(2);

        let nls = SUNNonlinSol_Newton(&ida_mem.ida_yy, &sunctx);
        assert_eq!(IDASetNonlinearSolver(&mut ida_mem, nls), IDA_SUCCESS);
        assert!(!ida_mem.ownNLS);
        assert!(ida_mem.nls_res.is_some());
        match ida_mem.NLS.as_ref().unwrap() {
            NonlinearSolver::Newton(ns) => assert_eq!(ns.maxiters, MAXIT),
            _ => unreachable!(),
        }

        /* NLS init succeeds with the solver attached */
        assert_eq!(idaNlsInit(&mut ida_mem), IDA_SUCCESS);
    }

    /* IDANls performs the full Newton solve on a linear DAE through
       idaNlsResidual/idaNlsLSetup(idaLsSetup)/idaNlsLSolve(idaLsSolve)/
       idaNlsConvTest: with a linear residual the first correction is
       exact and the second (zero) correction certifies convergence, so
       niters = 2, nsetups = 1 (first step forces lsetup), and yy/yp
       carry the corrected states (idas.c IDANls).  The IDAS deltas are
       asserted too: idaNlsLSetup clears ida_forceSetup and resets
       ida_ssS = TWENTY. */
    #[test]
    fn idanls_newton_linear_dae() {
        let sunctx = SUNContext::default();
        let mut ida_mem = make_ida_mem(2);

        /* predictor and step data */
        ida_mem.ida_yypredict = NVector::from_slice(&[1.0, 0.5]);
        ida_mem.ida_yppredict = NVector::from_slice(&[0.0, 0.0]);
        ida_mem.ida_cj = 1.0;
        ida_mem.ida_cjlast = 1.0;
        ida_mem.ida_epsNewt = 1.0e-8;
        ida_mem.ida_nst = 0; /* first step: forces callLSetup */
        ida_mem.ida_forceSetup = SUNTRUE; /* cleared by idaNlsLSetup (IDAS delta) */

        /* attach dense LS + Newton NLS */
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
        let nls = SUNNonlinSol_Newton(&ida_mem.ida_yy, &sunctx);
        assert_eq!(IDASetNonlinearSolver(&mut ida_mem, nls), IDA_SUCCESS);
        assert_eq!(idaNlsInit(&mut ida_mem), IDA_SUCCESS);

        let retval = IDANls(&mut ida_mem);
        assert_eq!(retval, IDA_SUCCESS);

        /* exact corrections: e = (4 - yppredict - 2*yypredict)/(2 + cj) */
        let e0 = (4.0 - 0.0 - 2.0 * 1.0) / 3.0; /* 2/3 */
        let e1 = (4.0 - 0.0 - 2.0 * 0.5) / 3.0; /* 1   */
        assert!((ida_mem.ida_ee.data[0] - e0).abs() < 1.0e-9);
        assert!((ida_mem.ida_ee.data[1] - e1).abs() < 1.0e-9);

        /* yy = yypredict + ee, yp = yppredict + cj*ee */
        assert!((ida_mem.ida_yy.data[0] - (1.0 + e0)).abs() < 1.0e-9);
        assert!((ida_mem.ida_yy.data[1] - (0.5 + e1)).abs() < 1.0e-9);
        assert!((ida_mem.ida_yp.data[0] - e0).abs() < 1.0e-9);
        assert!((ida_mem.ida_yp.data[1] - e1).abs() < 1.0e-9);

        /* residual at the solution vanishes and was saved in savres */
        assert!(ida_mem.ida_savres.data[0].abs() < 1.0e-8);
        assert!(ida_mem.ida_savres.data[1].abs() < 1.0e-8);

        /* counters: 2 Newton iterations (exact + certifying), 1 lsetup,
           2 residual evaluations (initial + after the first update; the
           certifying iteration converges before another res call), no
           failures */
        assert_eq!(ida_mem.ida_nni, 2);
        assert_eq!(ida_mem.ida_nnf, 0);
        assert_eq!(ida_mem.ida_nsetups, 1);
        assert_eq!(ida_mem.ida_nre, 2);

        /* lsetup refreshed the convergence-test constants */
        assert_eq!(ida_mem.ida_cjold, 1.0);
        assert_eq!(ida_mem.ida_cjratio, 1.0);

        /* IDAS deltas: forceSetup cleared, ssS reset alongside ss */
        assert!(!ida_mem.ida_forceSetup);
        assert_eq!(ida_mem.ida_ssS, 20.0);

        /* one Jacobian evaluation through the DQ path */
        let idals_mem = idaLs_AccessLMem(&mut ida_mem, "test").unwrap();
        assert_eq!(idals_mem.nje, 1);
        assert!(matches!(idals_mem.J, Some(SUNMatrix::Dense(_))));
        assert!(matches!(idals_mem.LS, LinearSolver::Dense(_)));
    }

    /* a fixed-point NLS is rejected: IDA requires ROOTFIND type */
    #[test]
    fn idasetnonlinearsolver_rejects_fixedpoint() {
        let sunctx = SUNContext::default();
        let mut ida_mem = make_ida_mem(2);
        let nls = crate::sunnonlinsol_fixedpoint::SUNNonlinSol_FixedPoint(&ida_mem.ida_yy, 0,
                                                                          &sunctx);
        assert_eq!(IDASetNonlinearSolver(&mut ida_mem, nls), IDA_ILL_INPUT);
        assert!(ida_mem.NLS.is_none());
    }
}
