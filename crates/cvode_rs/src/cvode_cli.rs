/* ---------------------------------------------------------------
 * Translated from src/cvode/cvode_cli.c (CVODE 7.7.0), together
 * with the sunCheckAndSet* helpers it uses from
 * src/sundials/sundials_cli.c.
 *
 * Provides command-line control over optional inputs to CVODE:
 * arguments of the form "<cvid>.<key>" followed by value
 * argument(s) are mapped onto the corresponding CVodeSet* calls.
 * The key names are identical to the C dispatch tables.
 * ---------------------------------------------------------------*/
use crate::cvode::CVodeSStolerances;
use crate::cvode_impl::*;
use crate::cvode_io::*;
use crate::cvode_ls::{
    CVodeSetDeltaGammaMaxBadJac, CVodeSetEpsLin, CVodeSetJacEvalFrequency,
    CVodeSetLSNormFactor, CVodeSetLinearSolutionScaling,
};
use crate::cvode_proj::{
    CVodeSetEpsProj, CVodeSetMaxNumProjFails, CVodeSetProjErrEst, CVodeSetProjFailEta,
    CVodeSetProjFrequency,
};

/* Set-routine signatures used in the dispatch tables (sundials_cli.h) */
type SunIntSetFn = fn(&mut CVodeMem, i32) -> i32;
type SunLongSetFn = fn(&mut CVodeMem, i64) -> i32;
type SunRealSetFn = fn(&mut CVodeMem, f64) -> i32;
type SunTwoRealSetFn = fn(&mut CVodeMem, f64, f64) -> i32;
type SunActionSetFn = fn(&mut CVodeMem) -> i32;

/*---------------------------------------------------------------
  CVodeSetOptions:

  Sets CVODE options using strings.
  ---------------------------------------------------------------*/

pub fn CVodeSetOptions(cv_mem: &mut CVodeMem, cvid: &str, file_name: &str, args: &[String]) -> i32 {
    if !file_name.is_empty() {
        let retval = CV_ILL_INPUT;
        cvProcessError(Some(cv_mem), retval, line!(), "CVodeSetOptions", file!(),
                       "file-based options are not currently supported.");
        return retval;
    }

    if !args.is_empty() {
        let retval = cvSetFromCommandLine(cv_mem, cvid, args);
        if retval != CV_SUCCESS {
            return retval;
        }
    }

    CV_SUCCESS
}

/* -----------------------------------------------------------------
 * Function to control CVODE options from the command line
 */

fn cvSetFromCommandLine(cv_mem: &mut CVodeMem, cvid: &str, argv: &[String]) -> i32 {
    /* Set lists of command-line arguments, and the corresponding set routines */
    let int_pairs: &[(&str, SunIntSetFn)] = &[
        ("max_conv_fails", CVodeSetMaxConvFails),
        ("max_err_test_fails", CVodeSetMaxErrTestFails),
        ("max_hnil_warns", CVodeSetMaxHnilWarns),
        ("max_nonlin_iters", CVodeSetMaxNonlinIters),
        ("max_order", CVodeSetMaxOrd),
        ("stab_lim_det", |m, v| CVodeSetStabLimDet(m, v != 0)),
        ("interpolate_stop_time", |m, v| CVodeSetInterpolateStopTime(m, v != 0)),
        ("use_integrator_fused_kernels", |m, v| {
            CVodeSetUseIntegratorFusedKernels(m, v != 0)
        }),
        ("num_fails_eta_max_err_fail", CVodeSetNumFailsEtaMaxErrFail),
        ("linear_solution_scaling", |m, v| CVodeSetLinearSolutionScaling(m, v != 0)),
        ("proj_err_est", |m, v| CVodeSetProjErrEst(m, v != 0)),
        ("max_num_proj_fails", CVodeSetMaxNumProjFails),
        ("max_num_constraint_fails", CVodeSetMaxNumConstraintFails),
    ];

    let long_pairs: &[(&str, SunLongSetFn)] = &[
        ("lsetup_frequency", CVodeSetLSetupFrequency),
        ("max_num_steps", CVodeSetMaxNumSteps),
        ("monitor_frequency", CVodeSetMonitorFrequency),
        ("num_steps_eta_max_early_step", CVodeSetNumStepsEtaMaxEarlyStep),
        ("jac_eval_frequency", CVodeSetJacEvalFrequency),
        ("proj_frequency", CVodeSetProjFrequency),
    ];

    let real_pairs: &[(&str, SunRealSetFn)] = &[
        ("delta_gamma_max_lsetup", CVodeSetDeltaGammaMaxLSetup),
        ("init_step", CVodeSetInitStep),
        ("max_step", CVodeSetMaxStep),
        ("min_step", CVodeSetMinStep),
        ("stop_time", CVodeSetStopTime),
        ("nonlin_conv_coef", CVodeSetNonlinConvCoef),
        ("eta_max_first_step", CVodeSetEtaMaxFirstStep),
        ("eta_max_early_step", CVodeSetEtaMaxEarlyStep),
        ("eta_max", CVodeSetEtaMax),
        ("eta_min", CVodeSetEtaMin),
        ("eta_min_err_fail", CVodeSetEtaMinErrFail),
        ("eta_max_err_fail", CVodeSetEtaMaxErrFail),
        ("eta_conv_fail", CVodeSetEtaConvFail),
        ("delta_gamma_max_bad_jac", CVodeSetDeltaGammaMaxBadJac),
        ("eps_lin", CVodeSetEpsLin),
        ("ls_norm_factor", CVodeSetLSNormFactor),
        ("eps_proj", CVodeSetEpsProj),
        ("proj_fail_eta", CVodeSetProjFailEta),
    ];

    let tworeal_pairs: &[(&str, SunTwoRealSetFn)] = &[
        ("eta_fixed_step_bounds", CVodeSetEtaFixedStepBounds),
        ("scalar_tolerances", CVodeSStolerances),
    ];

    let action_pairs: &[(&str, SunActionSetFn)] = &[
        ("clear_stop_time", CVodeClearStopTime),
        ("no_inactive_root_warn", CVodeSetNoInactiveRootWarn),
    ];

    /* Prefix for options to set */
    let default_id = "cvode";
    let prefix = if !cvid.is_empty() {
        format!("{}.", cvid)
    } else {
        format!("{}.", default_id)
    };
    let offset = prefix.len();

    let mut idx: usize = 1;
    while idx < argv.len() {
        'this_arg: {
            /* skip command-line arguments that do not begin with correct prefix */
            if !argv[idx].starts_with(&prefix) {
                break 'this_arg;
            }

            let mut arg_used = false;
            let mut j: usize = 0;

            /* check all "int" command-line options */
            let retval = sunCheckAndSetIntArgs(cv_mem, &mut idx, argv, offset, int_pairs,
                                               &mut arg_used, &mut j);
            if retval != CV_SUCCESS {
                cvProcessError(Some(cv_mem), retval, line!(), "cvSetFromCommandLine", file!(),
                               &format!("error setting key: {}", int_pairs[j].0));
                return retval;
            }
            if arg_used {
                break 'this_arg;
            }

            /* check all long int command-line options */
            let retval = sunCheckAndSetLongArgs(cv_mem, &mut idx, argv, offset, long_pairs,
                                                &mut arg_used, &mut j);
            if retval != CV_SUCCESS {
                cvProcessError(Some(cv_mem), retval, line!(), "cvSetFromCommandLine", file!(),
                               &format!("error setting key: {}", long_pairs[j].0));
                return retval;
            }
            if arg_used {
                break 'this_arg;
            }

            /* check all real command-line options */
            let retval = sunCheckAndSetRealArgs(cv_mem, &mut idx, argv, offset, real_pairs,
                                                &mut arg_used, &mut j);
            if retval != CV_SUCCESS {
                cvProcessError(Some(cv_mem), retval, line!(), "cvSetFromCommandLine", file!(),
                               &format!("error setting key: {}", real_pairs[j].0));
                return retval;
            }
            if arg_used {
                break 'this_arg;
            }

            /* check all pair-of-real command-line options */
            let retval = sunCheckAndSetTwoRealArgs(cv_mem, &mut idx, argv, offset,
                                                   tworeal_pairs, &mut arg_used, &mut j);
            if retval != CV_SUCCESS {
                cvProcessError(Some(cv_mem), retval, line!(), "cvSetFromCommandLine", file!(),
                               &format!("error setting key: {}", tworeal_pairs[j].0));
                return retval;
            }
            if arg_used {
                break 'this_arg;
            }

            /* check all action command-line options */
            let retval = sunCheckAndSetActionArgs(cv_mem, &mut idx, argv, offset, action_pairs,
                                                  &mut arg_used, &mut j);
            if retval != CV_SUCCESS {
                cvProcessError(Some(cv_mem), retval, line!(), "cvSetFromCommandLine", file!(),
                               &format!("error setting key: {}", action_pairs[j].0));
                return retval;
            }
            if arg_used {
                break 'this_arg;
            }

            /* warn for uninterpreted cvid.X arguments */
            cvProcessError(Some(cv_mem), CV_WARNING, line!(), "cvSetFromCommandLine", file!(),
                           &format!("WARNING: key {} was not handled\n", argv[idx]));
        }
        idx += 1;
    }

    CV_SUCCESS
}

/*===============================================================
  Command-line input utility routines (sundials_cli.c)
  ===============================================================*/

/* atoi/atol: parse a leading (optionally signed) integer, 0 on failure */
fn sun_atol(s: &str) -> i64 {
    let t = s.trim_start();
    let bytes = t.as_bytes();
    let mut i = 0;
    if i < bytes.len() && (bytes[i] == b'+' || bytes[i] == b'-') {
        i += 1;
    }
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        i += 1;
    }
    t[..i].parse().unwrap_or(0)
}

fn sun_atoi(s: &str) -> i32 {
    sun_atol(s) as i32
}

/* SUNStrToReal (strtod): parse a real number, 0.0 on failure */
fn sun_strtod(s: &str) -> f64 {
    s.trim().parse().unwrap_or(0.0)
}

fn next_arg<'a>(argv: &'a [String], argidx: &mut usize) -> &'a str {
    *argidx += 1;
    argv.get(*argidx).map(String::as_str).unwrap_or("")
}

fn sunCheckAndSetIntArgs(
    cv_mem: &mut CVodeMem,
    argidx: &mut usize,
    argv: &[String],
    offset: usize,
    testpairs: &[(&str, SunIntSetFn)],
    arg_used: &mut bool,
    failedarg: &mut usize,
) -> i32 {
    for (j, (key, set)) in testpairs.iter().enumerate() {
        *arg_used = false;
        if &argv[*argidx][offset..] == *key {
            let iarg = sun_atoi(next_arg(argv, argidx));
            let retval = set(cv_mem, iarg);
            if retval != CV_SUCCESS {
                *failedarg = j;
                return retval;
            }
            *arg_used = true;
            return CV_SUCCESS;
        }
    }
    CV_SUCCESS
}

fn sunCheckAndSetLongArgs(
    cv_mem: &mut CVodeMem,
    argidx: &mut usize,
    argv: &[String],
    offset: usize,
    testpairs: &[(&str, SunLongSetFn)],
    arg_used: &mut bool,
    failedarg: &mut usize,
) -> i32 {
    for (j, (key, set)) in testpairs.iter().enumerate() {
        *arg_used = false;
        if &argv[*argidx][offset..] == *key {
            let iarg = sun_atol(next_arg(argv, argidx));
            let retval = set(cv_mem, iarg);
            if retval != CV_SUCCESS {
                *failedarg = j;
                return retval;
            }
            *arg_used = true;
            return CV_SUCCESS;
        }
    }
    CV_SUCCESS
}

fn sunCheckAndSetRealArgs(
    cv_mem: &mut CVodeMem,
    argidx: &mut usize,
    argv: &[String],
    offset: usize,
    testpairs: &[(&str, SunRealSetFn)],
    arg_used: &mut bool,
    failedarg: &mut usize,
) -> i32 {
    for (j, (key, set)) in testpairs.iter().enumerate() {
        *arg_used = false;
        if &argv[*argidx][offset..] == *key {
            let rarg = sun_strtod(next_arg(argv, argidx));
            let retval = set(cv_mem, rarg);
            if retval != CV_SUCCESS {
                *failedarg = j;
                return retval;
            }
            *arg_used = true;
            return CV_SUCCESS;
        }
    }
    CV_SUCCESS
}

fn sunCheckAndSetTwoRealArgs(
    cv_mem: &mut CVodeMem,
    argidx: &mut usize,
    argv: &[String],
    offset: usize,
    testpairs: &[(&str, SunTwoRealSetFn)],
    arg_used: &mut bool,
    failedarg: &mut usize,
) -> i32 {
    for (j, (key, set)) in testpairs.iter().enumerate() {
        *arg_used = false;
        if &argv[*argidx][offset..] == *key {
            let rarg1 = sun_strtod(next_arg(argv, argidx));
            let rarg2 = sun_strtod(next_arg(argv, argidx));
            let retval = set(cv_mem, rarg1, rarg2);
            if retval != CV_SUCCESS {
                *failedarg = j;
                return retval;
            }
            *arg_used = true;
            return CV_SUCCESS;
        }
    }
    CV_SUCCESS
}

fn sunCheckAndSetActionArgs(
    cv_mem: &mut CVodeMem,
    argidx: &mut usize,
    argv: &[String],
    offset: usize,
    testpairs: &[(&str, SunActionSetFn)],
    arg_used: &mut bool,
    failedarg: &mut usize,
) -> i32 {
    for (j, (key, set)) in testpairs.iter().enumerate() {
        *arg_used = false;
        if &argv[*argidx][offset..] == *key {
            let retval = set(cv_mem);
            if retval != CV_SUCCESS {
                *failedarg = j;
                return retval;
            }
            *arg_used = true;
            return CV_SUCCESS;
        }
    }
    CV_SUCCESS
}
