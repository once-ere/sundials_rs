/* -----------------------------------------------------------------
 * Translated from src/idas/idas_nls_stg.c (IDAS 7.7.0), together
 * with idas.c's IDASensNls driver and the Newton solve loop of
 * sunnonlinsol_newton.c that the C code reaches through the
 * SUNNonlinearSolver ops table: the STAGGERED corrector, which —
 * after the state nonlinear system has converged — iterates on the
 * composite vector [eeS[0..Ns]] (all Ns sensitivity corrections at
 * once, no state part).
 *
 * In C the composite vectors are N_Vector senswrappers of Ns
 * sub-vectors whose slots alias IDAMem storage:
 *   ypredictStg = [ida_yySpredict[0..Ns]]
 *   ycorStg     = [ida_eeS[0..Ns]]
 *   ewtStg      = [ida_ewtS[0..Ns]]
 * Per the pinned senswrapper decision (idas_impl.rs) those aliases
 * are NOT stored: this module operates directly on the IDAMem
 * fields, and only the solver-owned workspace (NewtonSolver.deltaS —
 * created by SUNNonlinSol_NewtonSens with Ns sub-vectors, one per
 * sensitivity) remains a real senswrapper.  Composite norms use the
 * C senswrapper WRMS semantics: MAX of the per-sub-vector norms
 * (init 0).
 *
 * Staggered-specific deltas vs idas_nls_sim.rs:
 *   - lsetup counts ida_nsetupsS (not ida_nsetups) and does NOT
 *     clear ida_forceSetup; its residual/temp arguments are
 *     ida_delta and tmpS1..3 (tmpS1/tmpS2 alias tempv1/tempv2);
 *   - lsolve/residual work on the Ns sensitivity systems only, with
 *     rescur/resval = ida_delta (the saved state residual);
 *   - the m == 0 direct convergence test uses toldel itself (no
 *     PT0001 factor) and the rate estimate reads/updates ida_ssS.
 *
 * Exported entry points for idas.rs:
 *   IDASetNonlinearSolverSensStg / idaNlsInitSensStg
 *   IDASensNls(ida_mem: &mut IDAMem) -> i32   (C idas.c IDASensNls,
 *       the staggered sensitivity solve dispatched from IDAStep)
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
const ZERO: f64 = 0.0; /* real 0.0 (idas.c) */
const ONE: f64 = 1.0; /* real 1.0    */
const TWENTY: f64 = 20.0; /* real 20.0   */

/* nonlinear solver parameters */
const MAXIT: i32 = 4; /* default max number of nonlinear iterations    */
const RATEMAX: f64 = 0.9; /* max convergence rate used in divergence check */

/* -----------------------------------------------------------------------------
 * Exported functions
 * ---------------------------------------------------------------------------*/

pub fn IDASetNonlinearSolverSensStg(ida_mem: &mut IDAMem, NLS: NonlinearSolver) -> i32 {
    /* (The C NULL-input and missing-ops checks — gettype/solve/setsysfn
       on the NLS — cannot fail here: the workspace enum implements
       them.) */

    /* check for allowed nonlinear solver types */
    if NLS.nls_type() != SUNNONLINEARSOLVER_ROOTFIND {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASetNonlinearSolverSensStg",
                        file!(), "NLS type must be SUNNONLINEARSOLVER_ROOTFIND");
        return IDA_ILL_INPUT;
    }

    /* check that sensitivities were initialized */
    if !ida_mem.ida_sensi {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASetNonlinearSolverSensStg",
                        file!(), MSG_NO_SENSI);
        return IDA_ILL_INPUT;
    }

    /* check that the staggered corrector was selected */
    if ida_mem.ida_ism != IDA_STAGGERED {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASetNonlinearSolverSensStg",
                        file!(), "Sensitivity solution method is not IDA_STAGGERED");
        return IDA_ILL_INPUT;
    }

    /* free any existing nonlinear solver: RAII (drop on overwrite) */

    /* set SUNNonlinearSolver pointer */
    ida_mem.NLSstg = Some(NLS);

    /* Set NLS ownership flag. If this function was called to attach the default
       NLS, IDA will set the flag to SUNTRUE after this function returns. */
    ida_mem.ownNLSstg = SUNFALSE;

    /* the nonlinear residual and convergence test functions
       (idaNlsResidualSensStg / idaNlsConvTestSensStg) are wired
       statically into the Newton solve driver below; the C
       SUNNonlinSolSetSysFn / SUNNonlinSolSetConvTestFn registrations
       have no Rust counterpart */

    /* set max allowed nonlinear iterations */
    let retval = ida_mem.NLSstg.as_mut().unwrap().set_max_iters(MAXIT);
    if retval != IDA_SUCCESS {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASetNonlinearSolverSensStg",
                        file!(), "Setting maximum number of nonlinear iterations failed");
        return IDA_ILL_INPUT;
    }

    /* create vector wrappers if necessary: per the pinned senswrapper
       decision the ypredictStg/ycorStg/ewtStg aliases are not stored —
       this module reads ida_yySpredict, ida_eeS and ida_ewtS directly —
       so only the allocation flag survives (the C attach-loop below it
       vanishes with them) */
    if !ida_mem.stgMallocDone {
        ida_mem.stgMallocDone = SUNTRUE;
    }

    IDA_SUCCESS
}

/* -----------------------------------------------------------------------------
 * Private functions
 * ---------------------------------------------------------------------------*/

pub fn idaNlsInitSensStg(ida_mem: &mut IDAMem) -> i32 {
    /* In C this wires the idaNlsLSetupSensStg/idaNlsLSolveSensStg
       wrapper functions into the NLS depending on whether
       ida_lsetup/ida_lsolve exist (infallible registrations); here the
       dispatch is dynamic through ida_has_lsetup / the LsModule enum. */

    /* initialize nonlinear solver */
    let retval = match ida_mem.NLSstg.as_mut() {
        Some(nls) => nls.initialize(),
        None => -1,
    };

    if retval != IDA_SUCCESS {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "idaNlsInitSensStg", file!(),
                        MSG_NLS_INIT_FAIL);
        return IDA_NLS_INIT_FAIL;
    }

    IDA_SUCCESS
}

/* idaNlsLSetupSensStg (idas_nls_stg.c): wrapper around the lsetup
   dispatch.  In C the lsetup call receives (yy, yp, delta,
   tmpS1..3); the vectors alias IDAMem fields (tmpS1/tmpS2 alias
   tempv1/tempv2 — pinned convention), so they are detached for the
   call.  NOTE: idaLsSetup internally uses ida_tempv1/tempv2/tempv3
   as its three temporaries, while the C staggered wrapper passes
   tmpS1/tmpS2/tmpS3 = tempv1/tempv2/tmpS3 — i.e. the band-DQ
   Jacobian's yptemp scratch lands in tempv3 here instead of tmpS3.
   All three are dead scratch between calls (every consumer rewrites
   them before reading), and the assembled Jacobian is bit-identical,
   so routing through idaLsSetup is behaviorally exact. */
fn idaNlsLSetupSensStg(ida_mem: &mut IDAMem, _jbad: bool, jcur: &mut bool) -> i32 {
    ida_mem.ida_nsetupsS += 1;

    let mut lmem = std::mem::take(&mut ida_mem.ida_lmem);
    let retval = match &mut lmem {
        LsModule::Ls(idals_mem) => {
            let yy = std::mem::take(&mut ida_mem.ida_yy);
            let yp = std::mem::take(&mut ida_mem.ida_yp);
            let delta = std::mem::take(&mut ida_mem.ida_delta);
            let r = idaLsSetup(ida_mem, idals_mem, &yy, &yp, &delta);
            ida_mem.ida_yy = yy;
            ida_mem.ida_yp = yp;
            ida_mem.ida_delta = delta;
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

/* idaNlsLSolveSensStg (idas_nls_stg.c): solves the Ns sensitivity
   linear systems (weight = ida_ewtS[is], rescur = ida_delta);
   deltaStg sub-vector is = the is-th sensitivity delta.  In C each
   lsolve call receives (delta, w, yy, yp, delta-residual); the
   weight and current-state vectors alias IDAMem fields, so they are
   detached for the calls.  (C maps every failing call's retval
   identically, so mapping the last nonzero retval after the restore
   block is behaviorally exact.) */
fn idaNlsLSolveSensStg(ida_mem: &mut IDAMem, deltaStg: &mut NVectorSensWrapper) -> i32 {
    let mut lmem = std::mem::take(&mut ida_mem.ida_lmem);
    let retval = match &mut lmem {
        LsModule::Ls(idals_mem) => {
            let ewtS = std::mem::take(&mut ida_mem.ida_ewtS);
            let yy = std::mem::take(&mut ida_mem.ida_yy);
            let yp = std::mem::take(&mut ida_mem.ida_yp);
            let delta = std::mem::take(&mut ida_mem.ida_delta);

            let mut r = 0;
            for is in 0..ida_mem.ida_Ns as usize {
                r = idaLsSolve(ida_mem, idals_mem, &mut deltaStg.vecs[is], &ewtS[is], &yy, &yp,
                               &delta);
                if r != 0 {
                    break;
                }
            }

            ida_mem.ida_ewtS = ewtS;
            ida_mem.ida_yy = yy;
            ida_mem.ida_yp = yp;
            ida_mem.ida_delta = delta;
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

/* idaNlsResidualSensStg (idas_nls_stg.c): updates yyS/ypS from the
   current correction ycorStg = [ida_eeS] (read directly from IDAMem —
   aliases not stored) and evaluates the Ns sensitivity residuals
   (into resStg sub-vectors), with resval = ida_delta (the saved
   state residual). */
fn idaNlsResidualSensStg(ida_mem: &mut IDAMem, resStg: &mut NVectorSensWrapper) -> i32 {
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

    /* evaluate sens residual (C: ida_resS(Ns, tn, yy, yp, delta, yyS,
       ypS, resS, user_dataS, tmpS1, tmpS2, tmpS3); pinned dispatch —
       resSDQ selects the internal DQ routine, which receives IDA_mem in
       place of the user_dataS self-pointer; tmpS1/tmpS2 alias
       tempv1/tempv2.  The argument vectors are IDAMem fields, taken out
       for the duration of the call — donor take()/restore pattern.) */
    let ns = ida_mem.ida_Ns;
    let tn = ida_mem.ida_tn;
    let yy = std::mem::take(&mut ida_mem.ida_yy);
    let yp = std::mem::take(&mut ida_mem.ida_yp);
    let delta = std::mem::take(&mut ida_mem.ida_delta);
    let yyS = std::mem::take(&mut ida_mem.ida_yyS);
    let ypS = std::mem::take(&mut ida_mem.ida_ypS);
    let mut tmp1 = std::mem::take(&mut ida_mem.ida_tempv1);
    let mut tmp2 = std::mem::take(&mut ida_mem.ida_tempv2);
    let mut tmp3 = std::mem::take(&mut ida_mem.ida_tmpS3);

    let retval = if ida_mem.ida_resSDQ {
        crate::idas::IDASensResDQ(ida_mem, ns, tn, &yy, &yp, &delta, &yyS, &ypS,
                                  &mut resStg.vecs, &mut tmp1, &mut tmp2, &mut tmp3)
    } else {
        let resS_fn = ida_mem.ida_resS.unwrap();
        resS_fn(ns, tn, &yy, &yp, &delta, &yyS, &ypS, &mut resStg.vecs,
                &mut ida_mem.ida_user_data, &mut tmp1, &mut tmp2, &mut tmp3)
    };

    ida_mem.ida_yy = yy;
    ida_mem.ida_yp = yp;
    ida_mem.ida_delta = delta;
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

/* idaNlsConvTestSensStg (idas_nls_stg.c).  m is the current nonlinear
   iteration count (C SUNNonlinSolGetCurIter); del/ewt are the C
   senswrapper composites — delStg is the solver-owned workspace and
   ewtStg = [ida_ewtS] is read directly from IDAMem; the C ycor input
   is unused.  Staggered deltas: the m == 0 direct test compares
   against toldel itself (no PT0001 factor) and the rate estimate
   reads/updates ida_ssS. */
fn idaNlsConvTestSensStg(
    ida_mem: &mut IDAMem,
    delStg: &NVectorSensWrapper,
    tol: f64,
    m: i32,
) -> i32 {
    /* compute the norm of the correction: N_VWrmsNorm(del, ewt) on
       senswrappers = MAX of the per-sub-vector WRMS norms (init 0) */
    let mut delnrm = ZERO;
    for is in 0..ida_mem.ida_Ns as usize {
        let tmp = N_VWrmsNorm(&delStg.vecs[is], &ida_mem.ida_ewtS[is]);
        if tmp > delnrm {
            delnrm = tmp;
        }
    }

    /* test for convergence, first directly, then with rate estimate. */
    if m == 0 {
        ida_mem.ida_oldnrm = delnrm;
        if delnrm <= ida_mem.ida_toldel {
            return SUN_SUCCESS;
        }
    } else {
        let rate = SUNRpowerR(delnrm / ida_mem.ida_oldnrm, ONE / (m as f64));
        if rate > RATEMAX {
            return SUN_NLS_CONV_RECVR;
        }
        ida_mem.ida_ssS = rate / (ONE - rate);
    }

    if ida_mem.ida_ssS * delnrm <= tol {
        return SUN_SUCCESS;
    }

    /* not yet converged */
    SUN_NLS_CONTINUE
}

/*
 * IDASensNls (idas.c)
 *
 * This routine attempts to solve, one by one, all the sensitivity
 * linear systems using nonlinear iterations and the linear solver
 * specified (Staggered approach).
 */
pub fn IDASensNls(ida_mem: &mut IDAMem) -> i32 {
    let callLSetup = SUNFALSE;

    /* initial guess for the correction to the predictor: the C
       N_VConst(ZERO, ycorStg) zeroes the composite [eeS]
       (senswrapper aliases are not stored — pinned convention) */
    for is in 0..ida_mem.ida_Ns as usize {
        N_VConst(ZERO, &mut ida_mem.ida_eeS[is]);
    }

    /* solve the nonlinear system */
    let mut nls = match ida_mem.NLSstg.take() {
        Some(nls) => nls,
        None => return IDA_MEM_NULL,
    };
    let tol = ida_mem.ida_epsNewt;
    let retval = match &mut nls {
        NonlinearSolver::Newton(ns) => idaNlsSolveNewtonSensStg(ida_mem, ns, tol, callLSetup),
        /* cannot arise: IDASetNonlinearSolverSensStg rejects non-ROOTFIND types */
        NonlinearSolver::FixedPoint(_) => IDA_NLS_FAIL,
    };

    /* increment counters */
    ida_mem.ida_nniS += nls.get_num_iters();
    ida_mem.ida_nnfS += nls.get_num_conv_fails();
    ida_mem.NLSstg = Some(nls);

    if retval != SUN_SUCCESS {
        ida_mem.ida_ncfnS += 1;
        return retval;
    }

    /* update using the final correction from the nonlinear solver
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

    retval
}

/* zero the composite correction: N_VConst(ZERO, ycorStg) */
fn idaNlsZeroYcorStg(ida_mem: &mut IDAMem) {
    for is in 0..ida_mem.ida_Ns as usize {
        N_VConst(ZERO, &mut ida_mem.ida_eeS[is]);
    }
}

/*
 * SUNNonlinSolSolve_Newton (sunnonlinsol_newton.c), specialized to
 * the IDAS staggered-corrector callbacks.  ycor = [ida_eeS],
 * w = [ida_ewtS] (read directly from IDAMem); the Newton update
 * workspace is the senswrapper ns.deltaS (Ns sub-vectors, from
 * SUNNonlinSol_NewtonSens).
 */
fn idaNlsSolveNewtonSensStg(
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
            let r = idaNlsResidualSensStg(ida_mem, &mut deltaS);
            ns.deltaS = deltaS;
            r
        };
        if retval != 0 {
            break 'outer;
        }

        /* if indicated, setup the linear system */
        if call_lsetup {
            let mut jcur = ns.jcur;
            retval = idaNlsLSetupSensStg(ida_mem, jbad, &mut jcur);
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
            retval = idaNlsLSolveSensStg(ida_mem, &mut ns.deltaS);
            if retval != 0 {
                break;
            }

            /* update the Newton iterate: ycorStg = ycorStg + deltaStg */
            for is in 0..ida_mem.ida_Ns as usize {
                ida_mem.ida_eeS[is].linear_sum_with(ONE, ONE, &ns.deltaS.vecs[is]);
            }

            /* test for convergence */
            retval = {
                let m = ns.curiter;
                let deltaS = std::mem::take(&mut ns.deltaS);
                let r = idaNlsConvTestSensStg(ida_mem, &deltaS, tol, m);
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
                let r = idaNlsResidualSensStg(ida_mem, &mut deltaS);
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
            idaNlsZeroYcorStg(ida_mem);
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
    use crate::idas_nls::{idaNlsInit, IDANls, IDASetNonlinearSolver};
    use crate::sundials_context::SUNContext;
    use crate::sunlinsol_dense::SUNLinSol_Dense;
    use crate::sunmatrix_dense::SUNDenseMatrix;
    use crate::sunnonlinsol_newton::{SUNNonlinSol_Newton, SUNNonlinSol_NewtonSens};

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

    fn make_ida_mem_sens_stg(n: usize, ns: usize) -> IDAMem {
        let mut ida_mem = IDAMem::default();
        ida_mem.ida_res = Some(resfn);
        ida_mem.ida_ewt = NVector::from_slice(&vec![1.0; n]);
        ida_mem.ida_yy = NVector::new(n);
        ida_mem.ida_yp = NVector::new(n);
        ida_mem.ida_ee = NVector::new(n);
        ida_mem.ida_savres = NVector::new(n);
        ida_mem.ida_delta = NVector::new(n);
        ida_mem.ida_yypredict = NVector::new(n);
        ida_mem.ida_yppredict = NVector::new(n);
        ida_mem.ida_tempv1 = NVector::new(n);
        ida_mem.ida_tempv2 = NVector::new(n);
        ida_mem.ida_tempv3 = NVector::new(n);
        ida_mem.ida_hh = 1.0e-2;
        /* FSA (staggered corrector) state */
        ida_mem.ida_sensi = SUNTRUE;
        ida_mem.ida_ism = IDA_STAGGERED;
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

    /* IDASetNonlinearSolverSensStg stores the solver, clears ownership,
       sets maxiters = MAXIT = 4 and marks stgMallocDone (idas_nls_stg.c
       lines 46-186); sensitivity-less memory and a non-STAGGERED ism
       are rejected.  (Unlike the sim variant it does NOT touch
       nls_res.) */
    #[test]
    fn idasetnonlinearsolversensstg_defaults_and_rejections() {
        let sunctx = SUNContext::default();
        let mut ida_mem = make_ida_mem_sens_stg(2, 1);

        /* rejection: sensitivities not initialized */
        ida_mem.ida_sensi = SUNFALSE;
        let nls = SUNNonlinSol_NewtonSens(1, &ida_mem.ida_yy, &sunctx);
        assert_eq!(IDASetNonlinearSolverSensStg(&mut ida_mem, nls), IDA_ILL_INPUT);
        assert!(ida_mem.NLSstg.is_none());
        ida_mem.ida_sensi = SUNTRUE;

        /* rejection: simultaneous corrector selected */
        ida_mem.ida_ism = IDA_SIMULTANEOUS;
        let nls = SUNNonlinSol_NewtonSens(1, &ida_mem.ida_yy, &sunctx);
        assert_eq!(IDASetNonlinearSolverSensStg(&mut ida_mem, nls), IDA_ILL_INPUT);
        assert!(ida_mem.NLSstg.is_none());
        ida_mem.ida_ism = IDA_STAGGERED;

        /* acceptance */
        let nls = SUNNonlinSol_NewtonSens(1, &ida_mem.ida_yy, &sunctx);
        assert_eq!(IDASetNonlinearSolverSensStg(&mut ida_mem, nls), IDA_SUCCESS);
        assert!(!ida_mem.ownNLSstg);
        assert!(ida_mem.stgMallocDone);
        match ida_mem.NLSstg.as_ref().unwrap() {
            NonlinearSolver::Newton(ns) => {
                assert_eq!(ns.maxiters, MAXIT);
                assert_eq!(ns.deltaS.vecs.len(), 1);
            }
            _ => unreachable!(),
        }

        /* NLS init succeeds with the solver attached */
        assert_eq!(idaNlsInitSensStg(&mut ida_mem), IDA_SUCCESS);
    }

    /* Staggered flow: IDANls (idas_nls.rs, plain path — ism is
       STAGGERED so sensi_sim is false) converges the state system
       first (calling lsetup on the first step), then IDASensNls
       iterates the Ns sensitivity systems on the lagged Jacobian
       without any further lsetup (idas.c IDASensNls: callLSetup =
       SUNFALSE), so nniS = 2, nrSe = 2, nsetupsS = 0. */
    #[test]
    fn idasensnls_newton_staggered_linear_dae() {
        let sunctx = SUNContext::default();
        let mut ida_mem = make_ida_mem_sens_stg(2, 1);

        /* predictor and step data */
        ida_mem.ida_yypredict = NVector::from_slice(&[1.0, 0.5]);
        ida_mem.ida_yppredict = NVector::from_slice(&[0.0, 0.0]);
        ida_mem.ida_cj = 1.0;
        ida_mem.ida_cjlast = 1.0;
        ida_mem.ida_epsNewt = 1.0e-8;
        ida_mem.ida_nst = 0; /* first step: forces callLSetup in IDANls */

        /* attach dense LS + state Newton NLS + staggered NLSstg */
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
        let nlss = SUNNonlinSol_NewtonSens(1, &ida_mem.ida_yy, &sunctx);
        assert_eq!(IDASetNonlinearSolverSensStg(&mut ida_mem, nlss), IDA_SUCCESS);
        assert_eq!(idaNlsInitSensStg(&mut ida_mem), IDA_SUCCESS);

        /* state solve, then staggered sensitivity solve */
        assert_eq!(IDANls(&mut ida_mem), IDA_SUCCESS);
        assert_eq!(IDASensNls(&mut ida_mem), IDA_SUCCESS);

        /* exact corrections: e = (4 - yppredict - 2*yypredict)/(2 + cj),
           eS = (1 - ypSpredict - 2*yySpredict)/(2 + cj) */
        let e0 = (4.0 - 0.0 - 2.0 * 1.0) / 3.0; /* 2/3 */
        let e1 = (4.0 - 0.0 - 2.0 * 0.5) / 3.0; /* 1   */
        let es = (1.0 - 0.0 - 2.0 * 0.0) / 3.0; /* 1/3 */
        assert!((ida_mem.ida_ee.data[0] - e0).abs() < 1.0e-9);
        assert!((ida_mem.ida_ee.data[1] - e1).abs() < 1.0e-9);
        assert!((ida_mem.ida_eeS[0].data[0] - es).abs() < 1.0e-9);
        assert!((ida_mem.ida_eeS[0].data[1] - es).abs() < 1.0e-9);

        /* yyS = yySpredict + eeS, ypS = ypSpredict + cj*eeS */
        assert!((ida_mem.ida_yyS[0].data[0] - es).abs() < 1.0e-9);
        assert!((ida_mem.ida_yyS[0].data[1] - es).abs() < 1.0e-9);
        assert!((ida_mem.ida_ypS[0].data[0] - es).abs() < 1.0e-9);
        assert!((ida_mem.ida_ypS[0].data[1] - es).abs() < 1.0e-9);

        /* counters: the state solve took 2 iterations/1 setup (donor
           behavior); the staggered solve adds 2 sens iterations, 2 sens
           residuals, no sens setups and no failures */
        assert_eq!(ida_mem.ida_nni, 2);
        assert_eq!(ida_mem.ida_nsetups, 1);
        assert_eq!(ida_mem.ida_nniS, 2);
        assert_eq!(ida_mem.ida_nnfS, 0);
        assert_eq!(ida_mem.ida_ncfnS, 0);
        assert_eq!(ida_mem.ida_nrSe, 2);
        assert_eq!(ida_mem.ida_nsetupsS, 0);
    }
}
