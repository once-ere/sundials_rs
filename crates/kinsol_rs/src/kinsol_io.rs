/* -----------------------------------------------------------------
 * Translated from src/kinsol/kinsol_io.c (KINSOL 7.7.0).
 * Optional input and output functions for the KINSOL solver.
 *
 * The C functions take `void* kinmem` and start with a NULL check;
 * here the memory is `&mut KINMem`, which cannot be null, so those
 * checks vanish (donor cvode_io.rs convention). All other checks,
 * defaults and messages are translated line-for-line.
 * -----------------------------------------------------------------*/
use crate::kinsol_impl::*;
use crate::nvector_serial::{NVector, N_VClone, N_VMaxNorm, N_VScale};
use crate::sundials_math::{SUNRpowerR, SUNRsqrt};
use crate::sundials_types::*;
use crate::sundials_utils::{fmt_e, fmt_g};

const ZERO: f64 = 0.0;
const POINT1: f64 = 0.1;
const ONETHIRD: f64 = 0.3333333333333333;
const TWOTHIRDS: f64 = 0.6666666666666667;
const POINT9: f64 = 0.9;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;
const TWOPT5: f64 = 2.5;

/*
 * =================================================================
 * KINSOL optional input functions
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Function : KINSetUserData
 * -----------------------------------------------------------------
 */

pub fn KINSetUserData(kin_mem: &mut KINMem, user_data: UserData) -> i32 {
    kin_mem.kin_user_data = user_data;

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetDamping
 * -----------------------------------------------------------------
 */

pub fn KINSetDamping(kin_mem: &mut KINMem, beta: f64) -> i32 {
    /* check for illegal input value */
    if beta <= ZERO {
        KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSetDamping", file!(),
                        "beta <= 0 illegal");
        return KIN_ILL_INPUT;
    }

    if beta < ONE {
        /* enable damping */
        kin_mem.kin_beta = beta;
        kin_mem.kin_damping = SUNTRUE;
    } else {
        /* disable damping */
        kin_mem.kin_beta = ONE;
        kin_mem.kin_damping = SUNFALSE;
    }

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetMAA
 * -----------------------------------------------------------------
 */

pub fn KINSetMAA(kin_mem: &mut KINMem, maa: i64) -> i32 {
    if maa < 0 {
        KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSetMAA", file!(),
                        MSG_BAD_MAA);
        return KIN_ILL_INPUT;
    }

    // To allow for setting the depth and max number of iterations in any order we
    // do not limit maa here and instead enforce maa < mxiter in the AA
    // initialization function (KINInitAA)
    kin_mem.kin_m_aa = maa;

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetDelayAA
 * -----------------------------------------------------------------
 */

pub fn KINSetDelayAA(kin_mem: &mut KINMem, delay: i64) -> i32 {
    /* check for illegal input value */
    if delay < 0 {
        KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSetDelayAA", file!(),
                        "delay < 0 illegal");
        return KIN_ILL_INPUT;
    }

    kin_mem.kin_delay_aa = delay;

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetOrthAA
 * -----------------------------------------------------------------
 */

pub fn KINSetOrthAA(kin_mem: &mut KINMem, orthaa: i32) -> i32 {
    if (orthaa < KIN_ORTH_MGS) || (orthaa > KIN_ORTH_DCGS2) {
        KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSetOrthAA", file!(),
                        MSG_BAD_ORTHAA);
        return KIN_ILL_INPUT;
    }

    kin_mem.kin_orth_aa = orthaa;

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetDampingAA
 * -----------------------------------------------------------------
 */

pub fn KINSetDampingAA(kin_mem: &mut KINMem, beta: f64) -> i32 {
    /* check for illegal input value */
    if beta <= ZERO {
        KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSetDampingAA", file!(),
                        "beta <= 0 illegal");
        return KIN_ILL_INPUT;
    }

    if beta < ONE {
        /* enable damping */
        kin_mem.kin_beta_aa = beta;
        kin_mem.kin_damping_aa = SUNTRUE;
    } else {
        /* disable damping */
        kin_mem.kin_beta_aa = ONE;
        kin_mem.kin_damping_aa = SUNFALSE;
    }

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetDampingFn
 * -----------------------------------------------------------------
 */

pub fn KINSetDampingFn(kin_mem: &mut KINMem, damping_fn: Option<KINDampingFn>) -> i32 {
    kin_mem.kin_damping_fn = damping_fn;

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetDepthFn
 * -----------------------------------------------------------------
 */

pub fn KINSetDepthFn(kin_mem: &mut KINMem, depth_fn: Option<KINDepthFn>) -> i32 {
    kin_mem.kin_depth_fn = depth_fn;

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetReturnNewest
 * -----------------------------------------------------------------
 */

pub fn KINSetReturnNewest(kin_mem: &mut KINMem, ret_newest: bool) -> i32 {
    kin_mem.kin_ret_newest = ret_newest;

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetNumMaxIters
 * -----------------------------------------------------------------
 */

pub fn KINSetNumMaxIters(kin_mem: &mut KINMem, mxiter: i64) -> i32 {
    if mxiter < 0 {
        KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSetNumMaxIters", file!(),
                        MSG_BAD_MXITER);
        return KIN_ILL_INPUT;
    }

    if mxiter == 0 {
        kin_mem.kin_mxiter = MXITER_DEFAULT;
    } else {
        kin_mem.kin_mxiter = mxiter;
    }

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetNoInitSetup
 * -----------------------------------------------------------------
 */

pub fn KINSetNoInitSetup(kin_mem: &mut KINMem, noInitSetup: bool) -> i32 {
    kin_mem.kin_noInitSetup = noInitSetup;

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetNoResMon
 * -----------------------------------------------------------------
 */

pub fn KINSetNoResMon(kin_mem: &mut KINMem, noResMon: bool) -> i32 {
    kin_mem.kin_noResMon = noResMon;

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetMaxSetupCalls
 * -----------------------------------------------------------------
 */

pub fn KINSetMaxSetupCalls(kin_mem: &mut KINMem, msbset: i64) -> i32 {
    if msbset < 0 {
        KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSetMaxSetupCalls", file!(),
                        MSG_BAD_MSBSET);
        return KIN_ILL_INPUT;
    }

    if msbset == 0 {
        kin_mem.kin_msbset = MSBSET_DEFAULT;
    } else {
        kin_mem.kin_msbset = msbset;
    }

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetMaxSubSetupCalls
 * -----------------------------------------------------------------
 */

pub fn KINSetMaxSubSetupCalls(kin_mem: &mut KINMem, msbsetsub: i64) -> i32 {
    if msbsetsub < 0 {
        KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSetMaxSubSetupCalls",
                        file!(), MSG_BAD_MSBSETSUB);
        return KIN_ILL_INPUT;
    }

    if msbsetsub == 0 {
        kin_mem.kin_msbset_sub = MSBSET_SUB_DEFAULT;
    } else {
        kin_mem.kin_msbset_sub = msbsetsub;
    }

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetEtaForm
 * -----------------------------------------------------------------
 */

pub fn KINSetEtaForm(kin_mem: &mut KINMem, etachoice: i32) -> i32 {
    if (etachoice != KIN_ETACONSTANT)
        && (etachoice != KIN_ETACHOICE1)
        && (etachoice != KIN_ETACHOICE2)
    {
        KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSetEtaForm", file!(),
                        MSG_BAD_ETACHOICE);
        return KIN_ILL_INPUT;
    }

    kin_mem.kin_etaflag = etachoice;

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetEtaConstValue
 * -----------------------------------------------------------------
 */

pub fn KINSetEtaConstValue(kin_mem: &mut KINMem, eta: f64) -> i32 {
    if (eta < ZERO) || (eta > ONE) {
        KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSetEtaConstValue", file!(),
                        MSG_BAD_ETACONST);
        return KIN_ILL_INPUT;
    }

    if eta == ZERO {
        kin_mem.kin_eta = POINT1;
    } else {
        kin_mem.kin_eta = eta;
    }

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetEtaParams
 * -----------------------------------------------------------------
 */

pub fn KINSetEtaParams(kin_mem: &mut KINMem, egamma: f64, ealpha: f64) -> i32 {
    if (ealpha <= ONE) || (ealpha > TWO) {
        if ealpha != ZERO {
            KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSetEtaParams", file!(),
                            MSG_BAD_ALPHA);
            return KIN_ILL_INPUT;
        }
    }

    if ealpha == ZERO {
        kin_mem.kin_eta_alpha = TWO;
    } else {
        kin_mem.kin_eta_alpha = ealpha;
    }

    if (egamma <= ZERO) || (egamma > ONE) {
        if egamma != ZERO {
            KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSetEtaParams", file!(),
                            MSG_BAD_GAMMA);
            return KIN_ILL_INPUT;
        }
    }

    if egamma == ZERO {
        kin_mem.kin_eta_gamma = POINT9;
    } else {
        kin_mem.kin_eta_gamma = egamma;
    }

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetResMonParams
 * -----------------------------------------------------------------
 */

pub fn KINSetResMonParams(kin_mem: &mut KINMem, omegamin: f64, omegamax: f64) -> i32 {
    /* check omegamin */

    if omegamin < ZERO {
        KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSetResMonParams", file!(),
                        MSG_BAD_OMEGA);
        return KIN_ILL_INPUT;
    }

    if omegamin == ZERO {
        kin_mem.kin_omega_min = OMEGA_MIN;
    } else {
        kin_mem.kin_omega_min = omegamin;
    }

    /* check omegamax */

    if omegamax < ZERO {
        KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSetResMonParams", file!(),
                        MSG_BAD_OMEGA);
        return KIN_ILL_INPUT;
    }

    if omegamax == ZERO {
        if kin_mem.kin_omega_min > OMEGA_MAX {
            KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSetResMonParams",
                            file!(), MSG_BAD_OMEGA);
            return KIN_ILL_INPUT;
        } else {
            kin_mem.kin_omega_max = OMEGA_MAX;
        }
    } else {
        if kin_mem.kin_omega_min > omegamax {
            KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSetResMonParams",
                            file!(), MSG_BAD_OMEGA);
            return KIN_ILL_INPUT;
        } else {
            kin_mem.kin_omega_max = omegamax;
        }
    }

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetResMonConstValue
 * -----------------------------------------------------------------
 */

pub fn KINSetResMonConstValue(kin_mem: &mut KINMem, omegaconst: f64) -> i32 {
    /* check omegaconst */

    if omegaconst < ZERO {
        KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSetResMonConstValue",
                        file!(), MSG_BAD_OMEGA);
        return KIN_ILL_INPUT;
    }

    /* Load omega value. A value of 0 will force using omega_min and omega_max */
    kin_mem.kin_omega = omegaconst;

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetNoMinEps
 * -----------------------------------------------------------------
 */

pub fn KINSetNoMinEps(kin_mem: &mut KINMem, noMinEps: bool) -> i32 {
    kin_mem.kin_noMinEps = noMinEps;

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetMaxNewtonStep
 * -----------------------------------------------------------------
 */

pub fn KINSetMaxNewtonStep(kin_mem: &mut KINMem, mxnewtstep: f64) -> i32 {
    if mxnewtstep < ZERO {
        KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSetMaxNewtonStep", file!(),
                        MSG_BAD_MXNEWTSTEP);
        return KIN_ILL_INPUT;
    }

    /* Note: passing a value of 0.0 will use the default
       value (computed in KINSolInit) */

    kin_mem.kin_mxnstepin = mxnewtstep;

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetMaxBetaFails
 * -----------------------------------------------------------------
 */

pub fn KINSetMaxBetaFails(kin_mem: &mut KINMem, mxnbcf: i64) -> i32 {
    if mxnbcf < 0 {
        KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSetMaxBetaFails", file!(),
                        MSG_BAD_MXNBCF);
        return KIN_ILL_INPUT;
    }

    if mxnbcf == 0 {
        kin_mem.kin_mxnbcf = MXNBCF_DEFAULT;
    } else {
        kin_mem.kin_mxnbcf = mxnbcf;
    }

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetRelErrFunc
 * -----------------------------------------------------------------
 */

pub fn KINSetRelErrFunc(kin_mem: &mut KINMem, relfunc: f64) -> i32 {
    let uround: f64;

    if relfunc < ZERO {
        KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSetRelErrFunc", file!(),
                        MSG_BAD_RELFUNC);
        return KIN_ILL_INPUT;
    }

    if relfunc == ZERO {
        uround = kin_mem.kin_uround;
        kin_mem.kin_sqrt_relfunc = SUNRsqrt(uround);
    } else {
        kin_mem.kin_sqrt_relfunc = SUNRsqrt(relfunc);
    }

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetFuncNormTol
 * -----------------------------------------------------------------
 */

pub fn KINSetFuncNormTol(kin_mem: &mut KINMem, fnormtol: f64) -> i32 {
    let uround: f64;

    if fnormtol < ZERO {
        KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSetFuncNormTol", file!(),
                        MSG_BAD_FNORMTOL);
        return KIN_ILL_INPUT;
    }

    if fnormtol == ZERO {
        uround = kin_mem.kin_uround;
        kin_mem.kin_fnormtol = SUNRpowerR(uround, ONETHIRD);
    } else {
        kin_mem.kin_fnormtol = fnormtol;
    }

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetScaledStepTol
 * -----------------------------------------------------------------
 */

pub fn KINSetScaledStepTol(kin_mem: &mut KINMem, scsteptol: f64) -> i32 {
    let uround: f64;

    if scsteptol < ZERO {
        KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSetScaledStepTol", file!(),
                        MSG_BAD_SCSTEPTOL);
        return KIN_ILL_INPUT;
    }

    if scsteptol == ZERO {
        uround = kin_mem.kin_uround;
        kin_mem.kin_scsteptol = SUNRpowerR(uround, TWOTHIRDS);
    } else {
        kin_mem.kin_scsteptol = scsteptol;
    }

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetConstraints
 * -----------------------------------------------------------------
 */

pub fn KINSetConstraints(kin_mem: &mut KINMem, constraints: Option<&NVector>) -> i32 {
    let constraints = match constraints {
        None => {
            if kin_mem.kin_constraintsSet {
                /* (C: N_VDestroy(kin_mem->kin_constraints)) */
                kin_mem.kin_constraints = NVector::default();
                kin_mem.kin_lrw -= kin_mem.kin_lrw1;
                kin_mem.kin_liw -= kin_mem.kin_liw1;
            }
            kin_mem.kin_constraintsSet = SUNFALSE;
            return KIN_SUCCESS;
        }
        Some(c) => c,
    };

    /* Check the constraints vector */

    let temptest = N_VMaxNorm(constraints);
    if temptest > TWOPT5 {
        KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSetConstraints", file!(),
                        MSG_BAD_CONSTRAINTS);
        return KIN_ILL_INPUT;
    }

    if !kin_mem.kin_constraintsSet {
        kin_mem.kin_constraints = N_VClone(constraints);
        kin_mem.kin_lrw += kin_mem.kin_lrw1;
        kin_mem.kin_liw += kin_mem.kin_liw1;
        kin_mem.kin_constraintsSet = SUNTRUE;
    }

    /* Load the constraint vector */

    N_VScale(ONE, constraints, &mut kin_mem.kin_constraints);

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINSetSysFunc
 * -----------------------------------------------------------------
 */

pub fn KINSetSysFunc(kin_mem: &mut KINMem, func: Option<KINSysFn>) -> i32 {
    if func.is_none() {
        KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSetSysFunc", file!(),
                        MSG_FUNC_NULL);
        return KIN_ILL_INPUT;
    }

    kin_mem.kin_func = func;

    KIN_SUCCESS
}

/*
 * =================================================================
 * KINSOL optional output functions
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Function : KINGetWorkSpace
 * -----------------------------------------------------------------
 */

pub fn KINGetWorkSpace(kin_mem: &mut KINMem, lenrw: &mut i64, leniw: &mut i64) -> i32 {
    *lenrw = kin_mem.kin_lrw;
    *leniw = kin_mem.kin_liw;

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetNumNonlinSolvIters
 * -----------------------------------------------------------------
 */

pub fn KINGetNumNonlinSolvIters(kin_mem: &mut KINMem, nniters: &mut i64) -> i32 {
    *nniters = kin_mem.kin_nni;

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetNumFuncEvals
 * -----------------------------------------------------------------
 */

pub fn KINGetNumFuncEvals(kin_mem: &mut KINMem, nfevals: &mut i64) -> i32 {
    *nfevals = kin_mem.kin_nfe;

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetNumBetaCondFails
 * -----------------------------------------------------------------
 */

pub fn KINGetNumBetaCondFails(kin_mem: &mut KINMem, nbcfails: &mut i64) -> i32 {
    *nbcfails = kin_mem.kin_nbcf;

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetNumBacktrackOps
 * -----------------------------------------------------------------
 */

pub fn KINGetNumBacktrackOps(kin_mem: &mut KINMem, nbacktr: &mut i64) -> i32 {
    *nbacktr = kin_mem.kin_nbktrk;

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetFuncNorm
 * -----------------------------------------------------------------
 */

pub fn KINGetFuncNorm(kin_mem: &mut KINMem, funcnorm: &mut f64) -> i32 {
    *funcnorm = kin_mem.kin_fnorm;

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetStepLength
 * -----------------------------------------------------------------
 */

pub fn KINGetStepLength(kin_mem: &mut KINMem, steplength: &mut f64) -> i32 {
    *steplength = kin_mem.kin_stepl;

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetUserData
 * -----------------------------------------------------------------
 */

pub fn KINGetUserData(kin_mem: &mut KINMem) -> &mut UserData {
    &mut kin_mem.kin_user_data
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

/*
 * -----------------------------------------------------------------
 * Function : KINPrintAllStats
 * -----------------------------------------------------------------
 */

pub fn KINPrintAllStats(
    kin_mem: &mut KINMem,
    outfile: &mut dyn std::io::Write,
    fmt: SUNOutputFormat,
) -> i32 {
    if fmt != SUN_OUTPUTFORMAT_TABLE && fmt != SUN_OUTPUTFORMAT_CSV {
        KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINPrintAllStats", file!(),
                        "Invalid formatting option.");
        return KIN_ILL_INPUT;
    }

    sunfprintf_long(outfile, fmt, SUNTRUE, "Nonlinear iters", kin_mem.kin_nni);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Nonlinear fn evals", kin_mem.kin_nfe);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Beta condition fails", kin_mem.kin_nbcf);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Backtrack operations", kin_mem.kin_nbktrk);
    sunfprintf_real(outfile, fmt, SUNFALSE, "Nonlinear fn norm", kin_mem.kin_fnorm);
    sunfprintf_real(outfile, fmt, SUNFALSE, "Step length", kin_mem.kin_stepl);

    /* linear solver stats */
    if let LsModule::Ls(kinls_mem) = &kin_mem.kin_lmem {
        sunfprintf_long(outfile, fmt, SUNFALSE, "Jac fn evals", kinls_mem.nje);
        sunfprintf_long(outfile, fmt, SUNFALSE, "LS Nonlinear fn evals", kinls_mem.nfeDQ);
        sunfprintf_long(outfile, fmt, SUNFALSE, "Prec setup evals", kinls_mem.npe);
        sunfprintf_long(outfile, fmt, SUNFALSE, "Prec solves", kinls_mem.nps);
        sunfprintf_long(outfile, fmt, SUNFALSE, "LS iters", kinls_mem.nli);
        sunfprintf_long(outfile, fmt, SUNFALSE, "LS fails", kinls_mem.ncfl);
        sunfprintf_long(outfile, fmt, SUNFALSE, "Jac-times evals", kinls_mem.njtimes);
        if kin_mem.kin_nni > 0 {
            sunfprintf_real(outfile, fmt, SUNFALSE, "LS iters per NLS iter",
                            kinls_mem.nli as f64 / kin_mem.kin_nni as f64);
            sunfprintf_real(outfile, fmt, SUNFALSE, "Jac evals per NLS iter",
                            kinls_mem.nje as f64 / kin_mem.kin_nni as f64);
            sunfprintf_real(outfile, fmt, SUNFALSE, "Prec evals per NLS iter",
                            kinls_mem.npe as f64 / kin_mem.kin_nni as f64);
        }
    }

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Function : KINGetReturnFlagName
 * -----------------------------------------------------------------
 */

pub fn KINGetReturnFlagName(flag: i64) -> String {
    let flag_i32 = i32::try_from(flag).unwrap_or(i32::MIN); /* out-of-range -> "NONE" */
    let name = match flag_i32 {
        KIN_SUCCESS => "KIN_SUCCESS",
        KIN_INITIAL_GUESS_OK => "KIN_INITIAL_GUESS_OK",
        KIN_STEP_LT_STPTOL => "KIN_STEP_LT_STPTOL",
        KIN_WARNING => "KIN_WARNING",
        KIN_MEM_NULL => "KIN_MEM_NULL",
        KIN_ILL_INPUT => "KIN_ILL_INPUT",
        KIN_NO_MALLOC => "KIN_NO_MALLOC",
        KIN_MEM_FAIL => "KIN_MEM_FAIL",
        KIN_LINESEARCH_NONCONV => "KIN_LINESEARCH_NONCONV",
        KIN_MAXITER_REACHED => "KIN_MAXITER_REACHED",
        KIN_MXNEWT_5X_EXCEEDED => "KIN_MXNEWT_5X_EXCEEDED",
        KIN_LINESEARCH_BCFAIL => "KIN_LINESEARCH_BCFAIL",
        KIN_LINSOLV_NO_RECOVERY => "KIN_LINSOLV_NO_RECOVERY",
        KIN_LINIT_FAIL => "KIN_LINIT_FAIL",
        KIN_LSETUP_FAIL => "KIN_LSETUP_FAIL",
        KIN_LSOLVE_FAIL => "KIN_LSOLVE_FAIL",
        KIN_SYSFUNC_FAIL => "KIN_SYSFUNC_FAIL",
        KIN_FIRST_SYSFUNC_ERR => "KIN_FIRST_SYSFUNC_ERR",
        KIN_REPTD_SYSFUNC_ERR => "KIN_REPTD_SYSFUNC_ERR",
        KIN_VECTOROP_ERR => "KIN_VECTOROP_ERR",
        KIN_CONTEXT_ERR => "KIN_CONTEXT_ERR",
        KIN_DAMPING_FN_ERR => "KIN_DAMPING_FN_ERR",
        /* (KIN_DEPTH_FN_ERR is not in the C switch — falls to "NONE") */
        _ => "NONE",
    };

    name.to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    /*
     * Optional-input validation: bounds checks, error returns and
     * default-restoring semantics (kinsol_io.c behavior).
     */

    #[test]
    fn set_damping_validation_and_defaults() {
        let mut kin_mem = KINMem::default();

        /* beta <= 0 illegal; state untouched */
        assert_eq!(KINSetDamping(&mut kin_mem, -1.0), KIN_ILL_INPUT);
        assert_eq!(KINSetDamping(&mut kin_mem, 0.0), KIN_ILL_INPUT);
        assert_eq!(kin_mem.kin_beta, 1.0);
        assert!(!kin_mem.kin_damping);

        /* beta < 1 enables damping */
        assert_eq!(KINSetDamping(&mut kin_mem, 0.5), KIN_SUCCESS);
        assert_eq!(kin_mem.kin_beta, 0.5);
        assert!(kin_mem.kin_damping);

        /* beta >= 1 disables damping */
        assert_eq!(KINSetDamping(&mut kin_mem, 1.5), KIN_SUCCESS);
        assert_eq!(kin_mem.kin_beta, 1.0);
        assert!(!kin_mem.kin_damping);
    }

    #[test]
    fn set_num_max_iters_validation_and_defaults() {
        let mut kin_mem = KINMem::default();

        assert_eq!(KINSetNumMaxIters(&mut kin_mem, -1), KIN_ILL_INPUT);
        assert_eq!(kin_mem.kin_mxiter, MXITER_DEFAULT);

        assert_eq!(KINSetNumMaxIters(&mut kin_mem, 50), KIN_SUCCESS);
        assert_eq!(kin_mem.kin_mxiter, 50);

        /* 0 restores the default */
        assert_eq!(KINSetNumMaxIters(&mut kin_mem, 0), KIN_SUCCESS);
        assert_eq!(kin_mem.kin_mxiter, MXITER_DEFAULT);
    }

    #[test]
    fn set_maa_orth_delay_validation() {
        let mut kin_mem = KINMem::default();

        assert_eq!(KINSetMAA(&mut kin_mem, -1), KIN_ILL_INPUT);
        assert_eq!(kin_mem.kin_m_aa, 0);
        /* maa is not limited by mxiter here (enforced in KINInitAA) */
        assert_eq!(KINSetMAA(&mut kin_mem, 5000), KIN_SUCCESS);
        assert_eq!(kin_mem.kin_m_aa, 5000);

        assert_eq!(KINSetOrthAA(&mut kin_mem, -1), KIN_ILL_INPUT);
        assert_eq!(KINSetOrthAA(&mut kin_mem, KIN_ORTH_DCGS2 + 1), KIN_ILL_INPUT);
        assert_eq!(kin_mem.kin_orth_aa, KIN_ORTH_MGS);
        assert_eq!(KINSetOrthAA(&mut kin_mem, KIN_ORTH_ICWY), KIN_SUCCESS);
        assert_eq!(kin_mem.kin_orth_aa, KIN_ORTH_ICWY);

        assert_eq!(KINSetDelayAA(&mut kin_mem, -1), KIN_ILL_INPUT);
        assert_eq!(KINSetDelayAA(&mut kin_mem, 3), KIN_SUCCESS);
        assert_eq!(kin_mem.kin_delay_aa, 3);
    }

    #[test]
    fn set_eta_params_validation_and_defaults() {
        let mut kin_mem = KINMem::default();

        /* etachoice must be one of the three enumerations */
        assert_eq!(KINSetEtaForm(&mut kin_mem, 0), KIN_ILL_INPUT);
        assert_eq!(kin_mem.kin_etaflag, KIN_ETACHOICE1);
        assert_eq!(KINSetEtaForm(&mut kin_mem, KIN_ETACONSTANT), KIN_SUCCESS);
        assert_eq!(kin_mem.kin_etaflag, KIN_ETACONSTANT);

        /* eta out of range */
        assert_eq!(KINSetEtaConstValue(&mut kin_mem, -0.1), KIN_ILL_INPUT);
        assert_eq!(KINSetEtaConstValue(&mut kin_mem, 1.1), KIN_ILL_INPUT);
        /* 0 restores the default POINT1 */
        assert_eq!(KINSetEtaConstValue(&mut kin_mem, 0.0), KIN_SUCCESS);
        assert_eq!(kin_mem.kin_eta, 0.1);
        assert_eq!(KINSetEtaConstValue(&mut kin_mem, 0.4), KIN_SUCCESS);
        assert_eq!(kin_mem.kin_eta, 0.4);

        /* ealpha must satisfy 1 < ealpha <= 2 (or 0 for default) */
        assert_eq!(KINSetEtaParams(&mut kin_mem, 0.5, 3.0), KIN_ILL_INPUT);
        assert_eq!(kin_mem.kin_eta_alpha, 2.0);
        assert_eq!(kin_mem.kin_eta_gamma, 0.9);
        /* egamma error path leaves the (already validated) alpha set,
           as in the C source */
        assert_eq!(KINSetEtaParams(&mut kin_mem, 1.5, 1.5), KIN_ILL_INPUT);
        assert_eq!(kin_mem.kin_eta_alpha, 1.5);
        assert_eq!(kin_mem.kin_eta_gamma, 0.9);
        /* 0/0 restores both defaults */
        assert_eq!(KINSetEtaParams(&mut kin_mem, 0.0, 0.0), KIN_SUCCESS);
        assert_eq!(kin_mem.kin_eta_alpha, 2.0);
        assert_eq!(kin_mem.kin_eta_gamma, 0.9);
    }

    #[test]
    fn set_res_mon_params_validation_and_defaults() {
        let mut kin_mem = KINMem::default();

        assert_eq!(KINSetResMonParams(&mut kin_mem, -1.0, 0.5), KIN_ILL_INPUT);
        assert_eq!(KINSetResMonParams(&mut kin_mem, 0.1, -1.0), KIN_ILL_INPUT);
        /* the failed omegamax check above already stored omegamin (C behavior) */
        assert_eq!(kin_mem.kin_omega_min, 0.1);

        /* omegamin > omegamax illegal */
        assert_eq!(KINSetResMonParams(&mut kin_mem, 0.5, 0.1), KIN_ILL_INPUT);
        /* omegamin (0.5) stored before the failure, omega_max untouched */
        assert_eq!(kin_mem.kin_omega_min, 0.5);
        assert_eq!(kin_mem.kin_omega_max, OMEGA_MAX);

        /* omegamin = 0 restores OMEGA_MIN; then omegamax = 0 restores OMEGA_MAX */
        assert_eq!(KINSetResMonParams(&mut kin_mem, 0.0, 0.0), KIN_SUCCESS);
        assert_eq!(kin_mem.kin_omega_min, OMEGA_MIN);
        assert_eq!(kin_mem.kin_omega_max, OMEGA_MAX);

        assert_eq!(KINSetResMonConstValue(&mut kin_mem, -1.0), KIN_ILL_INPUT);
        assert_eq!(KINSetResMonConstValue(&mut kin_mem, 0.7), KIN_SUCCESS);
        assert_eq!(kin_mem.kin_omega, 0.7);
    }

    #[test]
    fn set_tolerances_zero_restores_defaults() {
        let mut kin_mem = KINMem::default();

        assert_eq!(KINSetFuncNormTol(&mut kin_mem, -1.0), KIN_ILL_INPUT);
        assert_eq!(KINSetFuncNormTol(&mut kin_mem, 1e-6), KIN_SUCCESS);
        assert_eq!(kin_mem.kin_fnormtol, 1e-6);
        assert_eq!(KINSetFuncNormTol(&mut kin_mem, 0.0), KIN_SUCCESS);
        assert_eq!(kin_mem.kin_fnormtol, SUNRpowerR(kin_mem.kin_uround, ONETHIRD));

        assert_eq!(KINSetScaledStepTol(&mut kin_mem, -1.0), KIN_ILL_INPUT);
        assert_eq!(KINSetScaledStepTol(&mut kin_mem, 1e-8), KIN_SUCCESS);
        assert_eq!(kin_mem.kin_scsteptol, 1e-8);
        assert_eq!(KINSetScaledStepTol(&mut kin_mem, 0.0), KIN_SUCCESS);
        assert_eq!(kin_mem.kin_scsteptol, SUNRpowerR(kin_mem.kin_uround, TWOTHIRDS));

        assert_eq!(KINSetRelErrFunc(&mut kin_mem, -1.0), KIN_ILL_INPUT);
        assert_eq!(KINSetRelErrFunc(&mut kin_mem, 1e-10), KIN_SUCCESS);
        assert_eq!(kin_mem.kin_sqrt_relfunc, SUNRsqrt(1e-10));
        assert_eq!(KINSetRelErrFunc(&mut kin_mem, 0.0), KIN_SUCCESS);
        assert_eq!(kin_mem.kin_sqrt_relfunc, SUNRsqrt(kin_mem.kin_uround));
    }

    #[test]
    fn set_constraints_check_and_workspace() {
        let mut kin_mem = KINMem::default();
        kin_mem.kin_lrw1 = 4;
        kin_mem.kin_liw1 = 2;
        let (lrw0, liw0) = (kin_mem.kin_lrw, kin_mem.kin_liw);

        /* illegal values (|c| > 2.5) in the constraints vector */
        let mut c = NVector::from_slice(&[1.0, -1.0, 2.0, 3.0]);
        assert_eq!(KINSetConstraints(&mut kin_mem, Some(&c)), KIN_ILL_INPUT);
        assert!(!kin_mem.kin_constraintsSet);

        /* legal constraints are cloned in and the workspace grows */
        c.data.copy_from_slice(&[1.0, -1.0, 2.0, 0.0]);
        assert_eq!(KINSetConstraints(&mut kin_mem, Some(&c)), KIN_SUCCESS);
        assert!(kin_mem.kin_constraintsSet);
        assert_eq!(kin_mem.kin_lrw, lrw0 + 4);
        assert_eq!(kin_mem.kin_liw, liw0 + 2);
        assert_eq!(kin_mem.kin_constraints.data, vec![1.0, -1.0, 2.0, 0.0]);

        /* None disables constraints and releases the workspace */
        assert_eq!(KINSetConstraints(&mut kin_mem, None), KIN_SUCCESS);
        assert!(!kin_mem.kin_constraintsSet);
        assert_eq!(kin_mem.kin_lrw, lrw0);
        assert_eq!(kin_mem.kin_liw, liw0);
    }

    #[test]
    fn getters_report_counters() {
        let mut kin_mem = KINMem::default();
        kin_mem.kin_nni = 5;
        kin_mem.kin_nfe = 8;
        kin_mem.kin_nbcf = 1;
        kin_mem.kin_nbktrk = 2;
        kin_mem.kin_fnorm = 0.25;
        kin_mem.kin_stepl = 1.5;

        let (mut l, mut v) = (0i64, 0i64);
        assert_eq!(KINGetWorkSpace(&mut kin_mem, &mut l, &mut v), KIN_SUCCESS);
        assert_eq!((l, v), (17, 22));
        assert_eq!(KINGetNumNonlinSolvIters(&mut kin_mem, &mut l), KIN_SUCCESS);
        assert_eq!(l, 5);
        assert_eq!(KINGetNumFuncEvals(&mut kin_mem, &mut l), KIN_SUCCESS);
        assert_eq!(l, 8);
        assert_eq!(KINGetNumBetaCondFails(&mut kin_mem, &mut l), KIN_SUCCESS);
        assert_eq!(l, 1);
        assert_eq!(KINGetNumBacktrackOps(&mut kin_mem, &mut l), KIN_SUCCESS);
        assert_eq!(l, 2);
        let mut r = 0.0;
        assert_eq!(KINGetFuncNorm(&mut kin_mem, &mut r), KIN_SUCCESS);
        assert_eq!(r, 0.25);
        assert_eq!(KINGetStepLength(&mut kin_mem, &mut r), KIN_SUCCESS);
        assert_eq!(r, 1.5);
    }

    #[test]
    fn print_all_stats_table_format() {
        let mut kin_mem = KINMem::default();
        kin_mem.kin_nni = 5;
        kin_mem.kin_nfe = 8;
        kin_mem.kin_nbcf = 1;
        kin_mem.kin_nbktrk = 2;
        kin_mem.kin_fnorm = 0.25;
        kin_mem.kin_stepl = 1.5;

        let mut out: Vec<u8> = Vec::new();
        assert_eq!(
            KINPrintAllStats(&mut kin_mem, &mut out, SUN_OUTPUTFORMAT_TABLE),
            KIN_SUCCESS
        );
        let expected = "Nonlinear iters               = 5\n\
                        Nonlinear fn evals            = 8\n\
                        Beta condition fails          = 1\n\
                        Backtrack operations          = 2\n\
                        Nonlinear fn norm             = 0.25\n\
                        Step length                   = 1.5\n";
        assert_eq!(String::from_utf8(out).unwrap(), expected);
    }

    #[test]
    fn print_all_stats_csv_format() {
        let mut kin_mem = KINMem::default();
        kin_mem.kin_nni = 5;
        kin_mem.kin_nfe = 8;
        kin_mem.kin_nbcf = 1;
        kin_mem.kin_nbktrk = 2;
        kin_mem.kin_fnorm = 0.25;
        kin_mem.kin_stepl = -1.5;

        let mut out: Vec<u8> = Vec::new();
        assert_eq!(
            KINPrintAllStats(&mut kin_mem, &mut out, SUN_OUTPUTFORMAT_CSV),
            KIN_SUCCESS
        );
        let expected = "Nonlinear iters,5,\
                        Nonlinear fn evals,8,\
                        Beta condition fails,1,\
                        Backtrack operations,2,\
                        Nonlinear fn norm, 2.500000000000000e-01,\
                        Step length,-1.500000000000000e+00";
        assert_eq!(String::from_utf8(out).unwrap(), expected);
    }

    #[test]
    fn return_flag_names() {
        assert_eq!(KINGetReturnFlagName(0), "KIN_SUCCESS");
        assert_eq!(KINGetReturnFlagName(2), "KIN_STEP_LT_STPTOL");
        assert_eq!(KINGetReturnFlagName(99), "KIN_WARNING");
        assert_eq!(KINGetReturnFlagName(-2), "KIN_ILL_INPUT");
        assert_eq!(KINGetReturnFlagName(-18), "KIN_DAMPING_FN_ERR");
        /* KIN_DEPTH_FN_ERR (-19) is absent from the C switch */
        assert_eq!(KINGetReturnFlagName(-19), "NONE");
        assert_eq!(KINGetReturnFlagName(12345), "NONE");
    }
}
