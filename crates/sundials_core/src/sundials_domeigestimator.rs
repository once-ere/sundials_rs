/* -----------------------------------------------------------------
 * Translation of sundials-7.7.0/src/sundials/sundials_domeigestimator.c
 * (+ include/sundials/sundials_domeigestimator.h).
 *
 * The C generic SUNDomEigEstimator is a base "class" holding an ops
 * table; per the workspace architecture, the ops table becomes enum
 * dispatch. The C-object-model plumbing NewEmpty/FreeEmpty has no
 * Rust counterpart (enum variants are constructed directly by the
 * implementation constructors in sundomeigest_*.rs); Destroy is
 * ownership drop. FILE* becomes &mut dyn std::io::Write.
 *
 * The C SUNATimesFn pair (void* A_data + function pointer) becomes a
 * boxed closure of the workspace-wide ATimesFn type from
 * sundials_linearsolver.rs; A_data is captured by the closure, so the
 * separate void* argument of SetATimes drops out.
 *
 * Ops provided by each implementation (mirroring the C ops tables):
 *   Power:   setatimes, setmaxiters, setnumpreprocessiters, setreltol,
 *            setinitialguess, initialize, estimate, getres,
 *            getnumiters, getnumatimescalls, write, destroy
 *   Arnoldi: setatimes, setnumpreprocessiters, setinitialguess,
 *            initialize, estimate, getnumiters, getnumatimescalls,
 *            write, destroy
 * Neither implementation provides setoptions. Where an op is absent
 * the base-class default from the C dispatch applies (SUN_SUCCESS
 * no-op, or zero outputs for the getters).
 * -----------------------------------------------------------------*/

use crate::nvector_serial::NVector;
use crate::sundials_errors::{SUNErrCode, SUN_ERR_ARG_INCOMPATIBLE, SUN_SUCCESS};
use crate::sundials_linearsolver::ATimesFn;
use crate::sundials_math::SUNStrToReal;
use crate::sundomeigest_arnoldi::SUNDomEigEstimatorContent_Arnoldi;
use crate::sundomeigest_power::SUNDomEigEstimatorContent_Power;

/// The generic estimator: C's ops-table polymorphism as enum dispatch.
pub enum SUNDomEigEstimator<'a> {
    Power(SUNDomEigEstimatorContent_Power<'a>),
    Arnoldi(SUNDomEigEstimatorContent_Arnoldi<'a>),
}

/* -----------------------------------------------------------------
 * internal utility routines
 * ----------------------------------------------------------------- */

/// C `atol`: leading whitespace, optional sign, digits (0 if none).
fn atol(s: &str) -> i64 {
    let t = s.trim_start();
    let (neg, t) = match t.strip_prefix('-') {
        Some(r) => (true, r),
        None => (false, t.strip_prefix('+').unwrap_or(t)),
    };
    let mut v: i64 = 0;
    for ch in t.chars() {
        match ch.to_digit(10) {
            Some(d) => v = v.wrapping_mul(10).wrapping_add(d as i64),
            None => break,
        }
    }
    if neg {
        -v
    } else {
        v
    }
}

/// C `atoi` (same as atol, truncated to int).
fn atoi(s: &str) -> i32 {
    atol(s) as i32
}

/// C sunDEESetFromCommandLine: process base-class options. `args`
/// corresponds to argv (args[0] is the program name and is skipped).
fn sunDEESetFromCommandLine(
    dee: &mut SUNDomEigEstimator<'_>,
    did: Option<&str>,
    args: &[String],
) -> SUNErrCode {
    /* Prefix for options to set */
    let default_id = "sundomeigestimator";
    let prefix = match did {
        Some(id) if !id.is_empty() => format!("{}.", id),
        _ => format!("{}.", default_id),
    };
    let offset = prefix.len();

    let mut idx = 1;
    while idx < args.len() {
        let arg = &args[idx];
        /* skip command-line arguments that do not begin with correct prefix */
        if !arg.starts_with(&prefix) {
            idx += 1;
            continue;
        }
        let key = &arg[offset..];

        /* control over SetMaxIters function */
        if key == "max_iters" {
            idx += 1;
            let large = atol(&args[idx]);
            let retval = SUNDomEigEstimator_SetMaxIters(dee, large);
            if retval != SUN_SUCCESS {
                return retval;
            }
            idx += 1;
            continue;
        }

        /* control over SetNumPreprocessIters function */
        if key == "num_preprocess_iters" {
            idx += 1;
            let iarg = atoi(&args[idx]);
            let retval = SUNDomEigEstimator_SetNumPreprocessIters(dee, iarg);
            if retval != SUN_SUCCESS {
                return retval;
            }
            idx += 1;
            continue;
        }

        /* control over SetRelTol function */
        if key == "rel_tol" {
            idx += 1;
            let rarg = SUNStrToReal(&args[idx]);
            let retval = SUNDomEigEstimator_SetRelTol(dee, rarg);
            if retval != SUN_SUCCESS {
                return retval;
            }
            idx += 1;
            continue;
        }

        idx += 1;
    }
    SUN_SUCCESS
}

/* -----------------------------------------------------------------
 * Functions in the 'ops' structure
 * ----------------------------------------------------------------- */

pub fn SUNDomEigEstimator_SetATimes<'a>(
    dee: &mut SUNDomEigEstimator<'a>,
    a_times: Box<ATimesFn<'a>>,
) -> SUNErrCode {
    match dee {
        SUNDomEigEstimator::Power(_) => {
            crate::sundomeigest_power::SUNDomEigEstimator_SetATimes_Power(dee, a_times)
        }
        SUNDomEigEstimator::Arnoldi(_) => {
            crate::sundomeigest_arnoldi::SUNDomEigEstimator_SetATimes_Arnoldi(dee, a_times)
        }
    }
}

/// C SUNDomEigEstimator_SetOptions(DEE, Did, file_name, argc, argv).
pub fn SUNDomEigEstimator_SetOptions(
    dee: &mut SUNDomEigEstimator<'_>,
    did: Option<&str>,
    file_name: Option<&str>,
    args: &[String],
) -> SUNErrCode {
    /* File-based option control is currently unimplemented */
    if let Some(f) = file_name {
        if !f.is_empty() {
            return SUN_ERR_ARG_INCOMPATIBLE;
        }
    }

    /* First, process all base-class options */
    if !args.is_empty() {
        let retval = sunDEESetFromCommandLine(dee, did, args);
        if retval != SUN_SUCCESS {
            return retval;
        }
    }

    /* Second, ask the implementation to process any remaining options:
    neither the Power nor the Arnoldi ops table provides setoptions. */
    SUN_SUCCESS
}

pub fn SUNDomEigEstimator_SetMaxIters(
    dee: &mut SUNDomEigEstimator<'_>,
    max_iters: i64,
) -> SUNErrCode {
    match dee {
        SUNDomEigEstimator::Power(_) => {
            crate::sundomeigest_power::SUNDomEigEstimator_SetMaxIters_Power(dee, max_iters)
        }
        SUNDomEigEstimator::Arnoldi(_) => SUN_SUCCESS, /* op not provided */
    }
}

pub fn SUNDomEigEstimator_SetNumPreprocessIters(
    dee: &mut SUNDomEigEstimator<'_>,
    num_iters: i32,
) -> SUNErrCode {
    match dee {
        SUNDomEigEstimator::Power(_) => {
            crate::sundomeigest_power::SUNDomEigEstimator_SetNumPreprocessIters_Power(
                dee, num_iters,
            )
        }
        SUNDomEigEstimator::Arnoldi(_) => {
            crate::sundomeigest_arnoldi::SUNDomEigEstimator_SetNumPreprocessIters_Arnoldi(
                dee, num_iters,
            )
        }
    }
}

pub fn SUNDomEigEstimator_SetRelTol(dee: &mut SUNDomEigEstimator<'_>, rel_tol: f64) -> SUNErrCode {
    match dee {
        SUNDomEigEstimator::Power(_) => {
            crate::sundomeigest_power::SUNDomEigEstimator_SetRelTol_Power(dee, rel_tol)
        }
        SUNDomEigEstimator::Arnoldi(_) => SUN_SUCCESS, /* op not provided */
    }
}

pub fn SUNDomEigEstimator_SetInitialGuess(
    dee: &mut SUNDomEigEstimator<'_>,
    q: &NVector,
) -> SUNErrCode {
    match dee {
        SUNDomEigEstimator::Power(_) => {
            crate::sundomeigest_power::SUNDomEigEstimator_SetInitialGuess_Power(dee, q)
        }
        SUNDomEigEstimator::Arnoldi(_) => {
            crate::sundomeigest_arnoldi::SUNDomEigEstimator_SetInitialGuess_Arnoldi(dee, q)
        }
    }
}

pub fn SUNDomEigEstimator_Initialize(dee: &mut SUNDomEigEstimator<'_>) -> SUNErrCode {
    match dee {
        SUNDomEigEstimator::Power(_) => {
            crate::sundomeigest_power::SUNDomEigEstimator_Initialize_Power(dee)
        }
        SUNDomEigEstimator::Arnoldi(_) => {
            crate::sundomeigest_arnoldi::SUNDomEigEstimator_Initialize_Arnoldi(dee)
        }
    }
}

pub fn SUNDomEigEstimator_Estimate(
    dee: &mut SUNDomEigEstimator<'_>,
    lambda_r: &mut f64,
    lambda_i: &mut f64,
) -> SUNErrCode {
    /* C: SUN_ERR_NOT_IMPLEMENTED when the op is absent; both provide it */
    match dee {
        SUNDomEigEstimator::Power(_) => {
            crate::sundomeigest_power::SUNDomEigEstimator_Estimate_Power(dee, lambda_r, lambda_i)
        }
        SUNDomEigEstimator::Arnoldi(_) => {
            crate::sundomeigest_arnoldi::SUNDomEigEstimator_Estimate_Arnoldi(
                dee, lambda_r, lambda_i,
            )
        }
    }
}

pub fn SUNDomEigEstimator_GetRes(dee: &SUNDomEigEstimator<'_>, res: &mut f64) -> SUNErrCode {
    match dee {
        SUNDomEigEstimator::Power(_) => {
            crate::sundomeigest_power::SUNDomEigEstimator_GetRes_Power(dee, res)
        }
        SUNDomEigEstimator::Arnoldi(_) => {
            /* op not provided */
            *res = 0.0;
            SUN_SUCCESS
        }
    }
}

pub fn SUNDomEigEstimator_GetNumIters(
    dee: &SUNDomEigEstimator<'_>,
    num_iters: &mut i64,
) -> SUNErrCode {
    match dee {
        SUNDomEigEstimator::Power(_) => {
            crate::sundomeigest_power::SUNDomEigEstimator_GetNumIters_Power(dee, num_iters)
        }
        SUNDomEigEstimator::Arnoldi(_) => {
            crate::sundomeigest_arnoldi::SUNDomEigEstimator_GetNumIters_Arnoldi(dee, num_iters)
        }
    }
}

pub fn SUNDomEigEstimator_GetNumATimesCalls(
    dee: &SUNDomEigEstimator<'_>,
    num_atimes: &mut i64,
) -> SUNErrCode {
    match dee {
        SUNDomEigEstimator::Power(_) => {
            crate::sundomeigest_power::SUNDomEigEstimator_GetNumATimesCalls_Power(dee, num_atimes)
        }
        SUNDomEigEstimator::Arnoldi(_) => {
            crate::sundomeigest_arnoldi::SUNDomEigEstimator_GetNumATimesCalls_Arnoldi(
                dee, num_atimes,
            )
        }
    }
}

pub fn SUNDomEigEstimator_Write(
    dee: &SUNDomEigEstimator<'_>,
    outfile: &mut dyn std::io::Write,
) -> SUNErrCode {
    match dee {
        SUNDomEigEstimator::Power(_) => {
            crate::sundomeigest_power::SUNDomEigEstimator_Write_Power(dee, outfile)
        }
        SUNDomEigEstimator::Arnoldi(_) => {
            crate::sundomeigest_arnoldi::SUNDomEigEstimator_Write_Arnoldi(dee, outfile)
        }
    }
}

/// C SUNDomEigEstimator_Destroy: the per-implementation destroy ops
/// (SUNDomEigEstimator_Destroy_Power / _Arnoldi) free the workspace
/// vectors and the content/ops structures; in Rust all of that is
/// ownership drop.
pub fn SUNDomEigEstimator_Destroy(dee: SUNDomEigEstimator<'_>) -> SUNErrCode {
    drop(dee);
    SUN_SUCCESS
}
