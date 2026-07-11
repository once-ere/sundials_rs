/* -----------------------------------------------------------------
 * Translated from src/idas/idas_io.c (IDAS 7.7.0).
 * Optional input and output functions for the IDAS solver.
 *
 * The C functions take `void* ida_mem` and start with a NULL check;
 * here the memory is `&mut IDAMem`, which cannot be null, so those
 * checks vanish (donor ida_io.rs / cvode_io.rs convention). All
 * other checks, defaults and messages are translated line-for-line.
 *
 * NOT registered in lib.rs yet; registers together with idas.rs once
 * idas_nls / idas_nls_sim / idas_nls_stg and idas_ls land
 * (IDAPrintAllStats dispatches on the LsModule::Ls variant pinned
 * for idas_ls_impl.h).
 * -----------------------------------------------------------------*/
use crate::idas::ida_msg_g;
use crate::idas_impl::*;
use crate::nvector_serial::{NVector, N_VMaxNorm, N_VScale};
use crate::sundials_math::SUNRabs;
use crate::sundials_types::*;
use crate::sundials_utils::{fmt_e, fmt_g};

const ZERO: f64 = 0.0;
const HALF: f64 = 0.5;
const ONE: f64 = 1.0;
const TWOPT5: f64 = 2.5;

/*
 * =================================================================
 * IDA optional input functions
 * =================================================================
 */

pub fn IDASetDeltaCjLSetup(ida_mem: &mut IDAMem, dcj: f64) -> i32 {
    if dcj < ZERO || dcj >= ONE {
        ida_mem.ida_dcj = DCJ_DEFAULT;
    } else {
        ida_mem.ida_dcj = dcj;
    }

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDASetUserData(ida_mem: &mut IDAMem, user_data: UserData) -> i32 {
    ida_mem.ida_user_data = user_data;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDASetEtaFixedStepBounds(ida_mem: &mut IDAMem, eta_min_fx: f64, eta_max_fx: f64) -> i32 {
    /* set allowed value or use default */
    if eta_min_fx >= ZERO && eta_min_fx <= ONE {
        ida_mem.ida_eta_min_fx = eta_min_fx;
    } else {
        ida_mem.ida_eta_min_fx = ETA_MIN_FX_DEFAULT;
    }

    if eta_max_fx >= ONE {
        ida_mem.ida_eta_max_fx = eta_max_fx;
    } else {
        ida_mem.ida_eta_max_fx = ETA_MAX_FX_DEFAULT;
    }

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDASetEtaMax(ida_mem: &mut IDAMem, eta_max: f64) -> i32 {
    /* set allowed value or use default */
    if eta_max <= ONE {
        ida_mem.ida_eta_max = ETA_MAX_DEFAULT;
    } else {
        ida_mem.ida_eta_max = eta_max;
    }

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDASetEtaMin(ida_mem: &mut IDAMem, eta_min: f64) -> i32 {
    /* set allowed value or use default */
    if eta_min <= ZERO || eta_min >= ONE {
        ida_mem.ida_eta_min = ETA_MIN_DEFAULT;
    } else {
        ida_mem.ida_eta_min = eta_min;
    }

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDASetEtaLow(ida_mem: &mut IDAMem, eta_low: f64) -> i32 {
    /* set allowed value or use default */
    if eta_low <= ZERO || eta_low >= ONE {
        ida_mem.ida_eta_low = ETA_LOW_DEFAULT;
    } else {
        ida_mem.ida_eta_low = eta_low;
    }

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDASetEtaMinErrFail(ida_mem: &mut IDAMem, eta_min_ef: f64) -> i32 {
    /* set allowed value or use default */
    if eta_min_ef <= ZERO || eta_min_ef >= ONE {
        ida_mem.ida_eta_min_ef = ETA_MIN_EF_DEFAULT;
    } else {
        ida_mem.ida_eta_min_ef = eta_min_ef;
    }

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDASetEtaConvFail(ida_mem: &mut IDAMem, eta_cf: f64) -> i32 {
    /* set allowed value or use default */
    if eta_cf <= ZERO || eta_cf >= ONE {
        ida_mem.ida_eta_cf = ETA_CF_DEFAULT;
    } else {
        ida_mem.ida_eta_cf = eta_cf;
    }

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDASetMaxOrd(ida_mem: &mut IDAMem, maxord: i32) -> i32 {
    if maxord <= 0 {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASetMaxOrd", file!(),
                        MSG_NEG_MAXORD);
        return IDA_ILL_INPUT;
    }

    /* Cannot increase maximum order beyond the value that
    was used when allocating memory */
    let maxord_alloc = ida_mem.ida_maxord_alloc;

    if maxord > maxord_alloc {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASetMaxOrd", file!(),
                        MSG_BAD_MAXORD);
        return IDA_ILL_INPUT;
    }

    /* (C: SUNMIN(maxord, MAXORD_DEFAULT)) */
    ida_mem.ida_maxord = maxord.min(MAXORD_DEFAULT);

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDASetMaxNumSteps(ida_mem: &mut IDAMem, mxsteps: i64) -> i32 {
    /* Passing mxsteps=0 sets the default. Passing mxsteps<0 disables the test. */

    if mxsteps == 0 {
        ida_mem.ida_mxstep = MXSTEP_DEFAULT;
    } else {
        ida_mem.ida_mxstep = mxsteps;
    }

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDASetInitStep(ida_mem: &mut IDAMem, hin: f64) -> i32 {
    ida_mem.ida_hin = hin;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDASetMaxStep(ida_mem: &mut IDAMem, hmax: f64) -> i32 {
    if hmax < ZERO {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASetMaxStep", file!(),
                        MSG_NEG_HMAX);
        return IDA_ILL_INPUT;
    }

    /* Passing 0 sets hmax = infinity */
    if hmax == ZERO {
        ida_mem.ida_hmax_inv = HMAX_INV_DEFAULT;
        return IDA_SUCCESS;
    }

    ida_mem.ida_hmax_inv = ONE / hmax;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDASetMinStep(ida_mem: &mut IDAMem, hmin: f64) -> i32 {
    if hmin < ZERO {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASetMinStep", file!(),
                        MSG_NEG_HMIN);
        return IDA_ILL_INPUT;
    }

    /* Passing 0 sets hmin = zero */
    if hmin == ZERO {
        ida_mem.ida_hmin = HMIN_DEFAULT;
        return IDA_SUCCESS;
    }

    ida_mem.ida_hmin = hmin;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDASetStopTime(ida_mem: &mut IDAMem, tstop: f64) -> i32 {
    /* If IDASolve was called at least once, test if tstop is legal
     * (i.e. if it was not already passed).
     * If IDASetStopTime is called before the first call to IDASolve,
     * tstop will be checked in IDASolve. */
    if ida_mem.ida_nst > 0 {
        if (tstop - ida_mem.ida_tn) * ida_mem.ida_hh < ZERO {
            IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASetStopTime", file!(),
                            &ida_msg_g(MSG_BAD_TSTOP, &[tstop, ida_mem.ida_tn]));
            return IDA_ILL_INPUT;
        }
    }

    ida_mem.ida_tstop = tstop;
    ida_mem.ida_tstopset = SUNTRUE;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAClearStopTime(ida_mem: &mut IDAMem) -> i32 {
    ida_mem.ida_tstopset = SUNFALSE;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDASetNonlinConvCoef(ida_mem: &mut IDAMem, epcon: f64) -> i32 {
    if epcon <= ZERO {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASetNonlinConvCoef", file!(),
                        MSG_NEG_EPCON);
        return IDA_ILL_INPUT;
    }

    ida_mem.ida_epcon = epcon;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDASetMaxErrTestFails(ida_mem: &mut IDAMem, maxnef: i32) -> i32 {
    ida_mem.ida_maxnef = maxnef;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDASetMaxConvFails(ida_mem: &mut IDAMem, maxncf: i32) -> i32 {
    ida_mem.ida_maxncf = maxncf;

    IDA_SUCCESS
}

pub fn IDASetMaxNonlinIters(ida_mem: &mut IDAMem, maxcor: i32) -> i32 {
    /* Are we computing sensitivities with the simultaneous approach? */
    let sensi_sim = ida_mem.ida_sensi && ida_mem.ida_ism == IDA_SIMULTANEOUS;

    if sensi_sim {
        /* check that the NLS is non-NULL */
        if ida_mem.NLSsim.is_none() {
            IDAProcessError(None, IDA_MEM_FAIL, line!(), "IDASetMaxNonlinIters", file!(),
                            MSG_MEM_FAIL);
            return IDA_MEM_FAIL;
        }

        /* (C: SUNNonlinSolSetMaxIters(IDA_mem->NLSsim, maxcor)) */
        ida_mem.NLSsim.as_mut().unwrap().set_max_iters(maxcor)
    } else {
        /* check that the NLS is non-NULL */
        if ida_mem.NLS.is_none() {
            IDAProcessError(None, IDA_MEM_FAIL, line!(), "IDASetMaxNonlinIters", file!(),
                            MSG_MEM_FAIL);
            return IDA_MEM_FAIL;
        }

        /* (C: SUNNonlinSolSetMaxIters(IDA_mem->NLS, maxcor)) */
        ida_mem.NLS.as_mut().unwrap().set_max_iters(maxcor)
    }
}

/*-----------------------------------------------------------------*/

pub fn IDASetSuppressAlg(ida_mem: &mut IDAMem, suppressalg: bool) -> i32 {
    ida_mem.ida_suppressalg = suppressalg;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDASetId(ida_mem: &mut IDAMem, id: Option<&NVector>) -> i32 {
    /* Disable id */
    let id = match id {
        None => {
            if ida_mem.ida_idMallocDone {
                ida_mem.ida_id = NVector::default();
                ida_mem.ida_lrw -= ida_mem.ida_lrw1;
                ida_mem.ida_liw -= ida_mem.ida_liw1;
            }
            ida_mem.ida_idMallocDone = SUNFALSE;
            return IDA_SUCCESS;
        }
        Some(v) => v,
    };

    if !ida_mem.ida_idMallocDone {
        ida_mem.ida_id = NVector::new(id.len());
        ida_mem.ida_lrw += ida_mem.ida_lrw1;
        ida_mem.ida_liw += ida_mem.ida_liw1;
        ida_mem.ida_idMallocDone = SUNTRUE;
    }

    /* Load the id vector */
    N_VScale(ONE, id, &mut ida_mem.ida_id);

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDASetConstraints(ida_mem: &mut IDAMem, constraints: Option<&NVector>) -> i32 {
    /* Disable constraints */
    let constraints = match constraints {
        None => {
            if ida_mem.ida_constraintsSet {
                ida_mem.ida_constraints = NVector::default();
                ida_mem.ida_constraintsSet = SUNFALSE;
                ida_mem.ida_lrw -= ida_mem.ida_lrw1;
                ida_mem.ida_liw -= ida_mem.ida_liw1;
            }
            return IDA_SUCCESS;
        }
        Some(c) => c,
    };

    /* (The C code tests here that the required vector ops are defined;
       the serial NVector implements them all.) */

    /*  Check the constraints vector */
    let temptest = N_VMaxNorm(constraints);
    if (temptest > TWOPT5) || (temptest < HALF) {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASetConstraints", file!(),
                        MSG_BAD_CONSTR);
        return IDA_ILL_INPUT;
    }

    if !ida_mem.ida_constraintsSet {
        ida_mem.ida_constraints = NVector::new(constraints.len());
        ida_mem.ida_constraintsSet = SUNTRUE;
        ida_mem.ida_lrw += ida_mem.ida_lrw1;
        ida_mem.ida_liw += ida_mem.ida_liw1;
    }

    /* Load the constraints vector */
    N_VScale(ONE, constraints, &mut ida_mem.ida_constraints);

    IDA_SUCCESS
}

/*
 * IDASetMaxNumConstraintFails
 *
 * Set the maximum number of constraint failure allowed in a step
 */

pub fn IDASetMaxNumConstraintFails(ida_mem: &mut IDAMem, max_fails: i32) -> i32 {
    if max_fails <= 0 {
        ida_mem.max_constraint_fails = MAX_CONSTRAINT_FAILS;
    } else {
        ida_mem.max_constraint_fails = max_fails;
    }

    IDA_SUCCESS
}

/*
 * IDAGetNumConstraintFails
 *
 * Get the number of failed steps due to constraint violation
 */

pub fn IDAGetNumConstraintFails(ida_mem: &mut IDAMem, num_fails_out: &mut i64) -> i32 {
    *num_fails_out = ida_mem.constraint_fails;

    IDA_SUCCESS
}

/*
 * IDAGetNumConstraintCorrections
 *
 * Get the number of constraint corrections
 */

pub fn IDAGetNumConstraintCorrections(ida_mem: &mut IDAMem, num_corrections_out: &mut i64) -> i32 {
    *num_corrections_out = ida_mem.constraint_corrections;

    IDA_SUCCESS
}

/*
 * IDASetRootDirection
 *
 * Specifies the direction of zero-crossings to be monitored.
 * The default is to monitor both crossings.
 */

pub fn IDASetRootDirection(ida_mem: &mut IDAMem, rootdir: &[i32]) -> i32 {
    let nrt = ida_mem.ida_nrtfn;
    if nrt == 0 {
        IDAProcessError(None, IDA_ILL_INPUT, line!(), "IDASetRootDirection", file!(),
                        MSG_NO_ROOT);
        return IDA_ILL_INPUT;
    }

    for i in 0..nrt as usize {
        ida_mem.ida_rootdir[i] = rootdir[i];
    }

    IDA_SUCCESS
}

/*
 * IDASetNoInactiveRootWarn
 *
 * Disables issuing a warning if some root function appears
 * to be identically zero at the beginning of the integration
 */

pub fn IDASetNoInactiveRootWarn(ida_mem: &mut IDAMem) -> i32 {
    ida_mem.ida_mxgnull = 0;

    IDA_SUCCESS
}

/*
 * =================================================================
 * IDA IC optional input functions
 * =================================================================
 */

pub fn IDASetNonlinConvCoefIC(ida_mem: &mut IDAMem, epiccon: f64) -> i32 {
    if epiccon <= ZERO {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASetNonlinConvCoefIC", file!(),
                        MSG_BAD_EPICCON);
        return IDA_ILL_INPUT;
    }

    ida_mem.ida_epiccon = epiccon;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDASetMaxNumStepsIC(ida_mem: &mut IDAMem, maxnh: i32) -> i32 {
    if maxnh <= 0 {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASetMaxNumStepsIC", file!(),
                        MSG_BAD_MAXNH);
        return IDA_ILL_INPUT;
    }

    ida_mem.ida_maxnh = maxnh;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDASetMaxNumJacsIC(ida_mem: &mut IDAMem, maxnj: i32) -> i32 {
    if maxnj <= 0 {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASetMaxNumJacsIC", file!(),
                        MSG_BAD_MAXNJ);
        return IDA_ILL_INPUT;
    }

    ida_mem.ida_maxnj = maxnj;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDASetMaxNumItersIC(ida_mem: &mut IDAMem, maxnit: i32) -> i32 {
    if maxnit <= 0 {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASetMaxNumItersIC", file!(),
                        MSG_BAD_MAXNIT);
        return IDA_ILL_INPUT;
    }

    ida_mem.ida_maxnit = maxnit;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDASetMaxBacksIC(ida_mem: &mut IDAMem, maxbacks: i32) -> i32 {
    if maxbacks <= 0 {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASetMaxBacksIC", file!(),
                        MSG_IC_BAD_MAXBACKS);
        return IDA_ILL_INPUT;
    }

    ida_mem.ida_maxbacks = maxbacks;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDASetLineSearchOffIC(ida_mem: &mut IDAMem, lsoff: bool) -> i32 {
    ida_mem.ida_lsoff = lsoff;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDASetStepToleranceIC(ida_mem: &mut IDAMem, steptol: f64) -> i32 {
    if steptol <= ZERO {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASetStepToleranceIC", file!(),
                        MSG_BAD_STEPTOL);
        return IDA_ILL_INPUT;
    }

    ida_mem.ida_steptol = steptol;

    IDA_SUCCESS
}

/*
 * =================================================================
 * Quadrature optional input functions
 * =================================================================
 */

/*-----------------------------------------------------------------*/

pub fn IDASetQuadErrCon(ida_mem: &mut IDAMem, errconQ: bool) -> i32 {
    if !ida_mem.ida_quadMallocDone {
        /* (C passes NULL as the memory argument here) */
        IDAProcessError(None, IDA_NO_QUAD, line!(), "IDASetQuadErrCon", file!(),
                        MSG_NO_QUAD);
        return IDA_NO_QUAD;
    }

    ida_mem.ida_errconQ = errconQ;

    IDA_SUCCESS
}

/*
 * =================================================================
 * FSA optional input functions
 * =================================================================
 */

pub fn IDASetSensDQMethod(ida_mem: &mut IDAMem, DQtype: i32, DQrhomax: f64) -> i32 {
    if DQtype != IDA_CENTERED && DQtype != IDA_FORWARD {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASetSensDQMethod", file!(),
                        MSG_BAD_DQTYPE);
        return IDA_ILL_INPUT;
    }

    if DQrhomax < ZERO {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASetSensDQMethod", file!(),
                        MSG_BAD_DQRHO);
        return IDA_ILL_INPUT;
    }

    ida_mem.ida_DQtype = DQtype;
    ida_mem.ida_DQrhomax = DQrhomax;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDASetSensErrCon(ida_mem: &mut IDAMem, errconS: bool) -> i32 {
    ida_mem.ida_errconS = errconS;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDASetSensMaxNonlinIters(ida_mem: &mut IDAMem, maxcorS: i32) -> i32 {
    /* check that the NLS is non-NULL */
    if ida_mem.NLSstg.is_none() {
        IDAProcessError(None, IDA_MEM_FAIL, line!(), "IDASetSensMaxNonlinIters", file!(),
                        MSG_MEM_FAIL);
        return IDA_MEM_FAIL;
    }

    /* (C: SUNNonlinSolSetMaxIters(IDA_mem->NLSstg, maxcorS)) */
    ida_mem.NLSstg.as_mut().unwrap().set_max_iters(maxcorS)
}

/*-----------------------------------------------------------------*/

pub fn IDASetSensParams(ida_mem: &mut IDAMem, p: Option<&[f64]>, pbar: Option<&[f64]>,
                        plist: Option<&[i32]>) -> i32 {
    /* Was sensitivity initialized? */
    if !ida_mem.ida_sensMallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_SENS, line!(), "IDASetSensParams", file!(),
                        MSG_NO_SENSI);
        return IDA_NO_SENS;
    }

    let ns = ida_mem.ida_Ns as usize;

    /* Parameters */
    /* DEVIATION NOTE: C stores the USER'S pointer (ida_p = p), so the
       internal DQ residual's perturbations of ida_p[which] are visible
       to the user's res through its own parameter block.  This port
       copies into the owned Vec; the DQ path (IDASensRes1DQ /
       IDAQuadSensRhs1InternalDQ in idas.rs) reads and perturbs
       ida_mem.ida_p.  Rust example ports must therefore route their
       res parameter reads so the perturbation is seen (open point
       shared with cvodes_rs, to be settled at FSA example
       verification). */
    ida_mem.ida_p = match p {
        Some(s) => s.to_vec(),
        None => Vec::new(),
    };

    /* pbar */
    if let Some(pbar) = pbar {
        for is in 0..ns {
            if pbar[is] == ZERO {
                IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASetSensParams",
                                file!(), MSG_BAD_PBAR);
                return IDA_ILL_INPUT;
            }
            ida_mem.ida_pbar[is] = SUNRabs(pbar[is]);
        }
    } else {
        for is in 0..ns {
            ida_mem.ida_pbar[is] = ONE;
        }
    }

    /* plist */
    if let Some(plist) = plist {
        for is in 0..ns {
            if plist[is] < 0 {
                IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASetSensParams",
                                file!(), MSG_BAD_PLIST);
                return IDA_ILL_INPUT;
            }
            ida_mem.ida_plist[is] = plist[is];
        }
    } else {
        for is in 0..ns {
            ida_mem.ida_plist[is] = is as i32;
        }
    }

    IDA_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function: IDASetQuadSensErrCon
 * -----------------------------------------------------------------
 * IDASetQuadSensErrCon specifies if quadrature sensitivity variables
 * are considered or not in the error control.
 * -----------------------------------------------------------------
 */
pub fn IDASetQuadSensErrCon(ida_mem: &mut IDAMem, errconQS: bool) -> i32 {
    /* Was sensitivity initialized? */
    if !ida_mem.ida_sensMallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_SENS, line!(), "IDASetQuadSensErrCon", file!(),
                        MSG_NO_SENSI);
        return IDA_NO_SENS;
    }

    /* Was quadrature sensitivity initialized?
       (C quirk: this arm also uses MSG_NO_SENSI, not MSG_NO_QUADSENSI) */
    if !ida_mem.ida_quadSensMallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_QUADSENS, line!(), "IDASetQuadSensErrCon",
                        file!(), MSG_NO_SENSI);
        return IDA_NO_QUADSENS;
    }

    ida_mem.ida_errconQS = errconQS;

    IDA_SUCCESS
}

/*
 * =================================================================
 * IDA optional output functions
 * =================================================================
 */

pub fn IDAGetNumSteps(ida_mem: &mut IDAMem, nsteps: &mut i64) -> i32 {
    *nsteps = ida_mem.ida_nst;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetNumResEvals(ida_mem: &mut IDAMem, nrevals: &mut i64) -> i32 {
    *nrevals = ida_mem.ida_nre;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetNumLinSolvSetups(ida_mem: &mut IDAMem, nlinsetups: &mut i64) -> i32 {
    *nlinsetups = ida_mem.ida_nsetups;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetNumErrTestFails(ida_mem: &mut IDAMem, netfails: &mut i64) -> i32 {
    *netfails = ida_mem.ida_netf;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetNumBacktrackOps(ida_mem: &mut IDAMem, nbacktracks: &mut i64) -> i32 {
    *nbacktracks = ida_mem.ida_nbacktr as i64;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetConsistentIC(
    ida_mem: &mut IDAMem,
    yy0: Option<&mut NVector>,
    yp0: Option<&mut NVector>,
) -> i32 {
    if ida_mem.ida_kused != 0 {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAGetConsistentIC", file!(),
                        MSG_TOO_LATE);
        return IDA_ILL_INPUT;
    }

    if let Some(yy0) = yy0 {
        N_VScale(ONE, &ida_mem.ida_phi[0], yy0);
    }
    if let Some(yp0) = yp0 {
        N_VScale(ONE, &ida_mem.ida_phi[1], yp0);
    }

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetLastOrder(ida_mem: &mut IDAMem, klast: &mut i32) -> i32 {
    *klast = ida_mem.ida_kused;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetCurrentOrder(ida_mem: &mut IDAMem, kcur: &mut i32) -> i32 {
    *kcur = ida_mem.ida_kk;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetCurrentCj(ida_mem: &mut IDAMem, cj: &mut f64) -> i32 {
    *cj = ida_mem.ida_cj;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

/* (C returns a pointer to the internal work vector; the port returns a
   borrow of it — always IDA_SUCCESS in C.) */
pub fn IDAGetCurrentY(ida_mem: &IDAMem) -> &NVector {
    &ida_mem.ida_yy
}

/*-----------------------------------------------------------------*/

pub fn IDAGetCurrentYSens(ida_mem: &IDAMem) -> &[NVector] {
    &ida_mem.ida_yyS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetCurrentYp(ida_mem: &IDAMem) -> &NVector {
    &ida_mem.ida_yp
}

/*-----------------------------------------------------------------*/

pub fn IDAGetCurrentYpSens(ida_mem: &IDAMem) -> &[NVector] {
    &ida_mem.ida_ypS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetActualInitStep(ida_mem: &mut IDAMem, hinused: &mut f64) -> i32 {
    *hinused = ida_mem.ida_h0u;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetLastStep(ida_mem: &mut IDAMem, hlast: &mut f64) -> i32 {
    *hlast = ida_mem.ida_hused;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetCurrentStep(ida_mem: &mut IDAMem, hcur: &mut f64) -> i32 {
    *hcur = ida_mem.ida_hh;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetCurrentTime(ida_mem: &mut IDAMem, tcur: &mut f64) -> i32 {
    *tcur = ida_mem.ida_tn;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetTolScaleFactor(ida_mem: &mut IDAMem, tolsfact: &mut f64) -> i32 {
    *tolsfact = ida_mem.ida_tolsf;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetErrWeights(ida_mem: &mut IDAMem, eweight: &mut NVector) -> i32 {
    N_VScale(ONE, &ida_mem.ida_ewt, eweight);

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetEstLocalErrors(ida_mem: &mut IDAMem, ele: &mut NVector) -> i32 {
    N_VScale(ONE, &ida_mem.ida_ee, ele);

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetWorkSpace(ida_mem: &mut IDAMem, lenrw: &mut i64, leniw: &mut i64) -> i32 {
    *leniw = ida_mem.ida_liw;
    *lenrw = ida_mem.ida_lrw;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

#[allow(clippy::too_many_arguments)]
pub fn IDAGetIntegratorStats(
    ida_mem: &mut IDAMem,
    nsteps: &mut i64,
    nrevals: &mut i64,
    nlinsetups: &mut i64,
    netfails: &mut i64,
    klast: &mut i32,
    kcur: &mut i32,
    hinused: &mut f64,
    hlast: &mut f64,
    hcur: &mut f64,
    tcur: &mut f64,
) -> i32 {
    *nsteps = ida_mem.ida_nst;
    *nrevals = ida_mem.ida_nre;
    *nlinsetups = ida_mem.ida_nsetups;
    *netfails = ida_mem.ida_netf;
    *klast = ida_mem.ida_kused;
    *kcur = ida_mem.ida_kk;
    *hinused = ida_mem.ida_h0u;
    *hlast = ida_mem.ida_hused;
    *hcur = ida_mem.ida_hh;
    *tcur = ida_mem.ida_tn;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetNumGEvals(ida_mem: &mut IDAMem, ngevals: &mut i64) -> i32 {
    *ngevals = ida_mem.ida_nge;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetRootInfo(ida_mem: &mut IDAMem, rootsfound: &mut [i32]) -> i32 {
    let nrt = ida_mem.ida_nrtfn;

    for i in 0..nrt as usize {
        rootsfound[i] = ida_mem.ida_iroots[i];
    }

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetNumNonlinSolvIters(ida_mem: &mut IDAMem, nniters: &mut i64) -> i32 {
    *nniters = ida_mem.ida_nni;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetNumNonlinSolvConvFails(ida_mem: &mut IDAMem, nnfails: &mut i64) -> i32 {
    *nnfails = ida_mem.ida_nnf;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetNonlinSolvStats(ida_mem: &mut IDAMem, nniters: &mut i64, nnfails: &mut i64) -> i32 {
    *nniters = ida_mem.ida_nni;
    *nnfails = ida_mem.ida_nnf;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetNumStepSolveFails(ida_mem: &mut IDAMem, nncfails: &mut i64) -> i32 {
    *nncfails = ida_mem.ida_ncfn;

    IDA_SUCCESS
}

/*
 * =================================================================
 * Quadrature optional output functions
 * =================================================================
 */

/*-----------------------------------------------------------------*/

pub fn IDAGetQuadNumRhsEvals(ida_mem: &mut IDAMem, nrQevals: &mut i64) -> i32 {
    if !ida_mem.ida_quadr {
        IDAProcessError(Some(ida_mem), IDA_NO_QUAD, line!(), "IDAGetQuadNumRhsEvals", file!(),
                        MSG_NO_QUAD);
        return IDA_NO_QUAD;
    }

    *nrQevals = ida_mem.ida_nrQe;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetQuadNumErrTestFails(ida_mem: &mut IDAMem, nQetfails: &mut i64) -> i32 {
    if !ida_mem.ida_quadr {
        IDAProcessError(Some(ida_mem), IDA_NO_QUAD, line!(), "IDAGetQuadNumErrTestFails",
                        file!(), MSG_NO_QUAD);
        return IDA_NO_QUAD;
    }

    *nQetfails = ida_mem.ida_netfQ;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetQuadErrWeights(ida_mem: &mut IDAMem, eQweight: &mut NVector) -> i32 {
    if !ida_mem.ida_quadr {
        IDAProcessError(Some(ida_mem), IDA_NO_QUAD, line!(), "IDAGetQuadErrWeights", file!(),
                        MSG_NO_QUAD);
        return IDA_NO_QUAD;
    }

    if ida_mem.ida_errconQ {
        N_VScale(ONE, &ida_mem.ida_ewtQ, eQweight);
    }

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetQuadStats(ida_mem: &mut IDAMem, nrQevals: &mut i64, nQetfails: &mut i64) -> i32 {
    if !ida_mem.ida_quadr {
        IDAProcessError(Some(ida_mem), IDA_NO_QUAD, line!(), "IDAGetQuadStats", file!(),
                        MSG_NO_QUAD);
        return IDA_NO_QUAD;
    }

    *nrQevals = ida_mem.ida_nrQe;
    *nQetfails = ida_mem.ida_netfQ;

    IDA_SUCCESS
}

/*
 * =================================================================
 * Quadrature FSA optional output functions
 * =================================================================
 */

/*-----------------------------------------------------------------*/

pub fn IDAGetQuadSensNumRhsEvals(ida_mem: &mut IDAMem, nrhsQSevals: &mut i64) -> i32 {
    if !ida_mem.ida_quadr_sensi {
        IDAProcessError(Some(ida_mem), IDA_NO_QUADSENS, line!(), "IDAGetQuadSensNumRhsEvals",
                        file!(), MSG_NO_QUADSENSI);
        return IDA_NO_QUADSENS;
    }

    *nrhsQSevals = ida_mem.ida_nrQSe;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetQuadSensNumErrTestFails(ida_mem: &mut IDAMem, nQSetfails: &mut i64) -> i32 {
    if !ida_mem.ida_quadr_sensi {
        IDAProcessError(Some(ida_mem), IDA_NO_QUADSENS, line!(),
                        "IDAGetQuadSensNumErrTestFails", file!(), MSG_NO_QUADSENSI);
        return IDA_NO_QUADSENS;
    }

    *nQSetfails = ida_mem.ida_netfQS;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetQuadSensErrWeights(ida_mem: &mut IDAMem, eQSweight: &mut [NVector]) -> i32 {
    if !ida_mem.ida_quadr_sensi {
        IDAProcessError(Some(ida_mem), IDA_NO_QUADSENS, line!(), "IDAGetQuadSensErrWeights",
                        file!(), MSG_NO_QUADSENSI);
        return IDA_NO_QUADSENS;
    }
    let ns = ida_mem.ida_Ns as usize;

    if ida_mem.ida_errconQS {
        for is in 0..ns {
            N_VScale(ONE, &ida_mem.ida_ewtQS[is], &mut eQSweight[is]);
        }
    }

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetQuadSensStats(ida_mem: &mut IDAMem, nrhsQSevals: &mut i64,
                           nQSetfails: &mut i64) -> i32 {
    if !ida_mem.ida_quadr_sensi {
        IDAProcessError(Some(ida_mem), IDA_NO_QUADSENS, line!(), "IDAGetQuadSensStats",
                        file!(), MSG_NO_QUADSENSI);
        return IDA_NO_QUADSENS;
    }

    *nrhsQSevals = ida_mem.ida_nrQSe;
    *nQSetfails = ida_mem.ida_netfQS;

    IDA_SUCCESS
}

/*
 * =================================================================
 * FSA optional output functions
 * =================================================================
 */

/*-----------------------------------------------------------------*/

pub fn IDAGetSensConsistentIC(
    ida_mem: &mut IDAMem,
    yyS0: Option<&mut [NVector]>,
    ypS0: Option<&mut [NVector]>,
) -> i32 {
    if !ida_mem.ida_sensi {
        IDAProcessError(Some(ida_mem), IDA_NO_SENS, line!(), "IDAGetSensConsistentIC", file!(),
                        MSG_NO_SENSI);
        return IDA_NO_SENS;
    }

    if ida_mem.ida_kused != 0 {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAGetSensConsistentIC",
                        file!(), MSG_TOO_LATE);
        return IDA_ILL_INPUT;
    }

    if let Some(yyS0) = yyS0 {
        for is in 0..ida_mem.ida_Ns as usize {
            N_VScale(ONE, &ida_mem.ida_phiS[0][is], &mut yyS0[is]);
        }
    }

    if let Some(ypS0) = ypS0 {
        for is in 0..ida_mem.ida_Ns as usize {
            N_VScale(ONE, &ida_mem.ida_phiS[1][is], &mut ypS0[is]);
        }
    }

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetSensNumResEvals(ida_mem: &mut IDAMem, nrSevals: &mut i64) -> i32 {
    if !ida_mem.ida_sensi {
        IDAProcessError(Some(ida_mem), IDA_NO_SENS, line!(), "IDAGetSensNumResEvals", file!(),
                        MSG_NO_SENSI);
        return IDA_NO_SENS;
    }

    *nrSevals = ida_mem.ida_nrSe;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetNumResEvalsSens(ida_mem: &mut IDAMem, nrevalsS: &mut i64) -> i32 {
    if !ida_mem.ida_sensi {
        IDAProcessError(Some(ida_mem), IDA_NO_SENS, line!(), "IDAGetNumResEvalsSens", file!(),
                        MSG_NO_SENSI);
        return IDA_NO_SENS;
    }

    *nrevalsS = ida_mem.ida_nreS;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetSensNumErrTestFails(ida_mem: &mut IDAMem, nSetfails: &mut i64) -> i32 {
    if !ida_mem.ida_sensi {
        IDAProcessError(Some(ida_mem), IDA_NO_SENS, line!(), "IDAGetSensNumErrTestFails",
                        file!(), MSG_NO_SENSI);
        return IDA_NO_SENS;
    }

    *nSetfails = ida_mem.ida_netfS;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetSensNumLinSolvSetups(ida_mem: &mut IDAMem, nlinsetupsS: &mut i64) -> i32 {
    if !ida_mem.ida_sensi {
        IDAProcessError(Some(ida_mem), IDA_NO_SENS, line!(), "IDAGetSensNumLinSolvSetups",
                        file!(), MSG_NO_SENSI);
        return IDA_NO_SENS;
    }

    *nlinsetupsS = ida_mem.ida_nsetupsS;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetSensErrWeights(ida_mem: &mut IDAMem, eSweight: &mut [NVector]) -> i32 {
    if !ida_mem.ida_sensi {
        IDAProcessError(Some(ida_mem), IDA_NO_SENS, line!(), "IDAGetSensErrWeights", file!(),
                        MSG_NO_SENSI);
        return IDA_NO_SENS;
    }

    /* (unconditional per-is copies — unlike the quadrature weight
       getters, C does not gate this on errconS) */
    for is in 0..ida_mem.ida_Ns as usize {
        N_VScale(ONE, &ida_mem.ida_ewtS[is], &mut eSweight[is]);
    }

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetSensStats(ida_mem: &mut IDAMem, nrSevals: &mut i64, nrevalsS: &mut i64,
                       nSetfails: &mut i64, nlinsetupsS: &mut i64) -> i32 {
    if !ida_mem.ida_sensi {
        IDAProcessError(Some(ida_mem), IDA_NO_SENS, line!(), "IDAGetSensStats", file!(),
                        MSG_NO_SENSI);
        return IDA_NO_SENS;
    }

    *nrSevals = ida_mem.ida_nrSe;
    *nrevalsS = ida_mem.ida_nreS;
    *nSetfails = ida_mem.ida_netfS;
    *nlinsetupsS = ida_mem.ida_nsetupsS;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetSensNumNonlinSolvIters(ida_mem: &mut IDAMem, nSniters: &mut i64) -> i32 {
    if !ida_mem.ida_sensi {
        IDAProcessError(Some(ida_mem), IDA_NO_SENS, line!(), "IDAGetSensNumNonlinSolvIters",
                        file!(), MSG_NO_SENSI);
        return IDA_NO_SENS;
    }

    *nSniters = ida_mem.ida_nniS;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetSensNumNonlinSolvConvFails(ida_mem: &mut IDAMem, nSnfails: &mut i64) -> i32 {
    if !ida_mem.ida_sensi {
        IDAProcessError(Some(ida_mem), IDA_NO_SENS, line!(),
                        "IDAGetSensNumNonlinSolvConvFails", file!(), MSG_NO_SENSI);
        return IDA_NO_SENS;
    }

    *nSnfails = ida_mem.ida_nnfS;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetSensNonlinSolvStats(ida_mem: &mut IDAMem, nSniters: &mut i64,
                                 nSnfails: &mut i64) -> i32 {
    if !ida_mem.ida_sensi {
        IDAProcessError(Some(ida_mem), IDA_NO_SENS, line!(), "IDAGetSensNonlinSolvStats",
                        file!(), MSG_NO_SENSI);
        return IDA_NO_SENS;
    }

    *nSniters = ida_mem.ida_nniS;
    *nSnfails = ida_mem.ida_nnfS;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

pub fn IDAGetNumStepSensSolveFails(ida_mem: &mut IDAMem, nSncfails: &mut i64) -> i32 {
    /* (C quirk: no sensi gate, and this reads ida_ncfn — the DAE step
       convergence-failure counter — NOT ida_ncfnS; preserved as-is) */
    *nSncfails = ida_mem.ida_ncfn;

    IDA_SUCCESS
}

/*-----------------------------------------------------------------*/

/* (C writes the user_data pointer through an out-param; the port
   returns a borrow of the owned UserData, per donor convention.) */
pub fn IDAGetUserData(ida_mem: &mut IDAMem) -> &mut UserData {
    &mut ida_mem.ida_user_data
}

/*-----------------------------------------------------------------*/

pub fn IDAPrintAllStats(
    ida_mem: &mut IDAMem,
    outfile: &mut dyn std::io::Write,
    fmt: SUNOutputFormat,
) -> i32 {
    if fmt != SUN_OUTPUTFORMAT_TABLE && fmt != SUN_OUTPUTFORMAT_CSV {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAPrintAllStats", file!(),
                        "Invalid formatting option.");
        return IDA_ILL_INPUT;
    }

    /* step and method stats */
    sunfprintf_real(outfile, fmt, SUNTRUE, "Current time", ida_mem.ida_tn);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Steps", ida_mem.ida_nst);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Error test fails", ida_mem.ida_netf);
    sunfprintf_long(outfile, fmt, SUNFALSE, "NLS step fails", ida_mem.ida_ncfn);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Constraint fails", ida_mem.constraint_fails);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Constraint corrections",
                    ida_mem.constraint_corrections);
    sunfprintf_real(outfile, fmt, SUNFALSE, "Initial step size", ida_mem.ida_h0u);
    sunfprintf_real(outfile, fmt, SUNFALSE, "Last step size", ida_mem.ida_hused);
    sunfprintf_real(outfile, fmt, SUNFALSE, "Current step size", ida_mem.ida_hh);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Last method order", ida_mem.ida_kused as i64);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Current method order", ida_mem.ida_kk as i64);

    /* function evaluations */
    sunfprintf_long(outfile, fmt, SUNFALSE, "Residual fn evals", ida_mem.ida_nre);

    /* IC calculation stats */
    sunfprintf_long(outfile, fmt, SUNFALSE, "IC linesearch backtrack ops",
                    ida_mem.ida_nbacktr as i64);

    /* nonlinear solver stats */
    sunfprintf_long(outfile, fmt, SUNFALSE, "NLS iters", ida_mem.ida_nni);
    sunfprintf_long(outfile, fmt, SUNFALSE, "NLS fails", ida_mem.ida_nnf);
    if ida_mem.ida_nst > 0 {
        sunfprintf_real(outfile, fmt, SUNFALSE, "NLS iters per step",
                        ida_mem.ida_nre as f64 / ida_mem.ida_nst as f64);
    }

    /* linear solver stats */
    sunfprintf_long(outfile, fmt, SUNFALSE, "LS setups", ida_mem.ida_nsetups);
    if let LsModule::Ls(idals_mem) = &ida_mem.ida_lmem {
        sunfprintf_long(outfile, fmt, SUNFALSE, "Jac fn evals", idals_mem.nje);
        sunfprintf_long(outfile, fmt, SUNFALSE, "LS residual fn evals", idals_mem.nreDQ);
        sunfprintf_long(outfile, fmt, SUNFALSE, "Prec setup evals", idals_mem.npe);
        sunfprintf_long(outfile, fmt, SUNFALSE, "Prec solves", idals_mem.nps);
        sunfprintf_long(outfile, fmt, SUNFALSE, "LS iters", idals_mem.nli);
        sunfprintf_long(outfile, fmt, SUNFALSE, "LS fails", idals_mem.ncfl);
        sunfprintf_long(outfile, fmt, SUNFALSE, "Jac-times setups", idals_mem.njtsetup);
        sunfprintf_long(outfile, fmt, SUNFALSE, "Jac-times evals", idals_mem.njtimes);
        if ida_mem.ida_nni > 0 {
            sunfprintf_real(outfile, fmt, SUNFALSE, "LS iters per NLS iter",
                            idals_mem.nli as f64 / ida_mem.ida_nni as f64);
            sunfprintf_real(outfile, fmt, SUNFALSE, "Jac evals per NLS iter",
                            idals_mem.nje as f64 / ida_mem.ida_nni as f64);
            sunfprintf_real(outfile, fmt, SUNFALSE, "Prec evals per NLS iter",
                            idals_mem.npe as f64 / ida_mem.ida_nni as f64);
        }
    }

    /* rootfinding stats */
    sunfprintf_long(outfile, fmt, SUNFALSE, "Root fn evals", ida_mem.ida_nge);

    /* quadrature stats */
    if ida_mem.ida_quadr {
        sunfprintf_long(outfile, fmt, SUNFALSE, "Quad fn evals", ida_mem.ida_nrQe);
        sunfprintf_long(outfile, fmt, SUNFALSE, "Quad error test fails", ida_mem.ida_netfQ);
    }

    /* sensitivity stats */
    if ida_mem.ida_sensi {
        sunfprintf_long(outfile, fmt, SUNFALSE, "Sens fn evals", ida_mem.ida_nrSe);
        sunfprintf_long(outfile, fmt, SUNFALSE, "Sens residual fn evals", ida_mem.ida_nreS);
        sunfprintf_long(outfile, fmt, SUNFALSE, "Sens error test fails", ida_mem.ida_netfS);
        if ida_mem.ida_ism == IDA_STAGGERED {
            sunfprintf_long(outfile, fmt, SUNFALSE, "Sens NLS iters", ida_mem.ida_nniS);
            sunfprintf_long(outfile, fmt, SUNFALSE, "Sens NLS fails", ida_mem.ida_nnfS);
            sunfprintf_long(outfile, fmt, SUNFALSE, "Sens NLS step fails", ida_mem.ida_ncfnS);
        }
        sunfprintf_long(outfile, fmt, SUNFALSE, "Sens LS setups", ida_mem.ida_nsetupsS);
    }

    /* quadrature-sensitivity stats */
    if ida_mem.ida_quadr_sensi {
        sunfprintf_long(outfile, fmt, SUNFALSE, "QuadSens residual evals", ida_mem.ida_nrQSe);
        sunfprintf_long(outfile, fmt, SUNFALSE, "QuadSens error test fails",
                        ida_mem.ida_netfQS);
    }

    IDA_SUCCESS
}

/*
 * =================================================================
 * IDAGetReturnFlagName
 * =================================================================
 */

pub fn IDAGetReturnFlagName(flag: i64) -> String {
    let name = match flag as i32 {
        IDA_SUCCESS => "IDA_SUCCESS",
        IDA_TSTOP_RETURN => "IDA_TSTOP_RETURN",
        IDA_ROOT_RETURN => "IDA_ROOT_RETURN",
        IDA_TOO_MUCH_WORK => "IDA_TOO_MUCH_WORK",
        IDA_TOO_MUCH_ACC => "IDA_TOO_MUCH_ACC",
        IDA_ERR_FAIL => "IDA_ERR_FAIL",
        IDA_CONV_FAIL => "IDA_CONV_FAIL",
        IDA_LINIT_FAIL => "IDA_LINIT_FAIL",
        IDA_LSETUP_FAIL => "IDA_LSETUP_FAIL",
        IDA_LSOLVE_FAIL => "IDA_LSOLVE_FAIL",
        IDA_CONSTR_FAIL => "IDA_CONSTR_FAIL",
        IDA_RES_FAIL => "IDA_RES_FAIL",
        IDA_FIRST_RES_FAIL => "IDA_FIRST_RES_FAIL",
        IDA_REP_RES_ERR => "IDA_REP_RES_ERR",
        IDA_RTFUNC_FAIL => "IDA_RTFUNC_FAIL",
        IDA_MEM_FAIL => "IDA_MEM_FAIL",
        IDA_MEM_NULL => "IDA_MEM_NULL",
        IDA_ILL_INPUT => "IDA_ILL_INPUT",
        IDA_NO_MALLOC => "IDA_NO_MALLOC",
        IDA_BAD_T => "IDA_BAD_T",
        IDA_BAD_K => "IDA_BAD_K",
        IDA_BAD_DKY => "IDA_BAD_DKY",
        IDA_BAD_EWT => "IDA_BAD_EWT",
        IDA_NO_RECOVERY => "IDA_NO_RECOVERY",
        IDA_LINESEARCH_FAIL => "IDA_LINESEARCH_FAIL",
        IDA_NO_SENS => "IDA_NO_SENS",
        IDA_SRES_FAIL => "IDA_SRES_FAIL",
        IDA_REP_SRES_ERR => "IDA_REP_SRES_ERR",
        IDA_BAD_IS => "IDA_BAD_IS",
        IDA_NO_QUAD => "IDA_NO_QUAD",
        IDA_NO_QUADSENS => "IDA_NO_QUADSENS",
        IDA_QRHS_FAIL => "IDA_QRHS_FAIL",
        IDA_REP_QRHS_ERR => "IDA_REP_QRHS_ERR",
        IDA_QSRHS_FAIL => "IDA_QSRHS_FAIL",
        IDA_REP_QSRHS_ERR => "IDA_REP_QSRHS_ERR",

        /* IDAA flags follow below. */
        IDA_NO_ADJ => "IDA_NO_ADJ",
        IDA_BAD_TB0 => "IDA_BAD_TB0",
        IDA_REIFWD_FAIL => "IDA_REIFWD_FAIL",
        IDA_FWD_FAIL => "IDA_FWD_FAIL",
        IDA_GETY_BADT => "IDA_GETY_BADT",
        IDA_NO_BCK => "IDA_NO_BCK",
        IDA_NO_FWD => "IDA_NO_FWD",
        IDA_NLS_SETUP_FAIL => "IDA_NLS_SETUP_FAIL",
        IDA_NLS_FAIL => "IDA_NLS_FAIL",
        _ => "NONE",
    };
    name.to_string()
}

/* -----------------------------------------------------------------
 * Counterparts of sunfprintf_real / sunfprintf_long
 * (src/sundials/sundials_utils.h). SUN_FORMAT_G is "%.15g" and
 * SUN_FORMAT_E is "% .15e" for double precision.
 * -----------------------------------------------------------------*/

const SUN_TABLE_WIDTH: usize = 29;

fn sunfprintf_real(
    outfile: &mut dyn std::io::Write,
    fmt: SUNOutputFormat,
    start: bool,
    name: &str,
    value: f64,
) {
    if fmt == SUN_OUTPUTFORMAT_TABLE {
        let _ = writeln!(outfile, "{:<width$} = {}", name, fmt_g(value, 0, 15),
                         width = SUN_TABLE_WIDTH);
    } else {
        if !start {
            let _ = write!(outfile, ",");
        }
        /* C "% .15e": a space is printed in place of a plus sign */
        let e = fmt_e(value, 0, 15);
        let e = if e.starts_with('-') { e } else { format!(" {}", e) };
        let _ = write!(outfile, "{},{}", name, e);
    }
}

fn sunfprintf_long(
    outfile: &mut dyn std::io::Write,
    fmt: SUNOutputFormat,
    start: bool,
    name: &str,
    value: i64,
) {
    if fmt == SUN_OUTPUTFORMAT_TABLE {
        let _ = writeln!(outfile, "{:<width$} = {}", name, value, width = SUN_TABLE_WIDTH);
    } else {
        if !start {
            let _ = write!(outfile, ",");
        }
        let _ = write!(outfile, "{},{}", name, value);
    }
}

/* END of idas_io.c port. */
