/* -----------------------------------------------------------------
 * Translation of src/arkode/arkode_cli.c (SUNDIALS 7.7.0).
 *
 * This file provides command-line control over optional inputs
 * to ARKODE.
 *
 * C's sunbooleantype-taking setters are listed in the int tables
 * (implicit int conversion); the Rust setters take bool, so those
 * entries go through local `cli_*` shims mapping nonzero -> true.
 * -----------------------------------------------------------------*/

use crate::arkode_impl::{
    arkProcessError, ARKodeMem, ARK_ACCUMERROR_AVG, ARK_ACCUMERROR_MAX, ARK_ACCUMERROR_NONE,
    ARK_ACCUMERROR_SUM, ARK_ILL_INPUT, ARK_INTERP_HERMITE, ARK_INTERP_LAGRANGE, ARK_INTERP_NONE,
    ARK_SUCCESS, ARK_WARNING,
};
use crate::arkode_io::{
    ARKodeClearStopTime, ARKodeResetAccumulatedError,
    ARKodeSetAccumulatedErrorType, ARKodeSetAdaptivityAdjustment, ARKodeSetAutonomous,
    ARKodeSetCFLFraction, ARKodeSetDeduceImplicitRhs, ARKodeSetDeltaGammaMax, ARKodeSetErrorBias,
    ARKodeSetFixedStep, ARKodeSetFixedStepBounds, ARKodeSetInitStep, ARKodeSetInterpolantDegree,
    ARKodeSetInterpolantType, ARKodeSetInterpolateStopTime, ARKodeSetLSetupFrequency,
    ARKodeSetLinear, ARKodeSetMaxCFailGrowth, ARKodeSetMaxConvFails, ARKodeSetMaxEFailGrowth,
    ARKodeSetMaxErrTestFails, ARKodeSetMaxFirstGrowth, ARKodeSetMaxGrowth, ARKodeSetMaxHnilWarns,
    ARKodeSetMaxNonlinIters, ARKodeSetMaxNumConstrFails, ARKodeSetMaxNumSteps, ARKodeSetMaxStep,
    ARKodeSetMinReduction, ARKodeSetMinStep, ARKodeSetNoInactiveRootWarn, ARKodeSetNonlinCRDown,
    ARKodeSetNonlinConvCoef, ARKodeSetNonlinRDiv, ARKodeSetNonlinear, ARKodeSetOrder,
    ARKodeSetPredictorMethod, ARKodeSetSafetyFactor, ARKodeSetSmallNumEFails,
    ARKodeSetStepDirection, ARKodeSetStopTime, ARKodeSetUseCompensatedSums,
    ARKodeWriteParameters,
};
use crate::arkode_ls::{
    ARKodeSetEpsLin, ARKodeSetJacEvalFrequency, ARKodeSetLSNormFactor,
    ARKodeSetLinearSolutionScaling, ARKodeSetMassEpsLin, ARKodeSetMassLSNormFactor,
};
use crate::sundials_cli::{
    sunCheckAndSetActionArgs, sunCheckAndSetIntArgs, sunCheckAndSetLongArgs,
    sunCheckAndSetRealArgs, sunCheckAndSetTwoRealArgs, sunKeyActionPair, sunKeyIntPair,
    sunKeyLongPair, sunKeyRealPair, sunKeyTwoRealPair,
};

/* sunbooleantype-taking setters wrapped for the int tables */
fn cli_set_autonomous(ark_mem: &mut ARKodeMem, arg: i32) -> i32 {
    ARKodeSetAutonomous(ark_mem, arg != 0)
}
fn cli_set_deduce_implicit_rhs(ark_mem: &mut ARKodeMem, arg: i32) -> i32 {
    ARKodeSetDeduceImplicitRhs(ark_mem, arg != 0)
}
fn cli_set_interpolate_stop_time(ark_mem: &mut ARKodeMem, arg: i32) -> i32 {
    ARKodeSetInterpolateStopTime(ark_mem, arg != 0)
}
fn cli_set_linear_solution_scaling(ark_mem: &mut ARKodeMem, arg: i32) -> i32 {
    ARKodeSetLinearSolutionScaling(ark_mem, arg != 0)
}
fn cli_set_use_compensated_sums(ark_mem: &mut ARKodeMem, arg: i32) -> i32 {
    ARKodeSetUseCompensatedSums(ark_mem, arg != 0)
}

/*---------------------------------------------------------------
  ARKodeSetOptions:

  Sets ARKODE options using strings.
  ---------------------------------------------------------------*/
pub fn ARKodeSetOptions(
    ark_mem: &mut ARKodeMem,
    arkid: Option<&str>,
    file_name: Option<&str>,
    argv: &[String],
) -> i32 {
    if let Some(fname) = file_name {
        if !fname.is_empty() {
            let retval = ARK_ILL_INPUT;
            arkProcessError(
                Some(ark_mem),
                retval,
                line!(),
                "ARKodeSetOptions",
                file!(),
                "file-based options are not currently supported.",
            );
            return retval;
        }
    }

    if !argv.is_empty() {
        let retval = arkSetFromCommandLine(ark_mem, arkid, argv);
        if retval != ARK_SUCCESS {
            return retval;
        }
    }

    ARK_SUCCESS
}

/* -----------------------------------------------------------------
 * Function to control ARKODE options from the command line
 */
fn arkSetFromCommandLine(ark_mem: &mut ARKodeMem, arkid: Option<&str>, argv: &[String]) -> i32 {
    /* Set lists of command-line arguments, and the corresponding set routines */
    let int_pairs: [sunKeyIntPair<ARKodeMem>; 17] = [
        sunKeyIntPair { key: "order", set: ARKodeSetOrder },
        sunKeyIntPair { key: "interpolant_degree", set: ARKodeSetInterpolantDegree },
        sunKeyIntPair { key: "linear", set: ARKodeSetLinear },
        sunKeyIntPair { key: "autonomous", set: cli_set_autonomous },
        sunKeyIntPair { key: "deduce_implicit_rhs", set: cli_set_deduce_implicit_rhs },
        sunKeyIntPair { key: "lsetup_frequency", set: ARKodeSetLSetupFrequency },
        sunKeyIntPair { key: "predictor_method", set: ARKodeSetPredictorMethod },
        sunKeyIntPair { key: "max_nonlin_iters", set: ARKodeSetMaxNonlinIters },
        sunKeyIntPair { key: "max_hnil_warns", set: ARKodeSetMaxHnilWarns },
        sunKeyIntPair { key: "interpolate_stop_time", set: cli_set_interpolate_stop_time },
        sunKeyIntPair { key: "max_num_constr_fails", set: ARKodeSetMaxNumConstrFails },
        sunKeyIntPair { key: "adaptivity_adjustment", set: ARKodeSetAdaptivityAdjustment },
        sunKeyIntPair { key: "small_num_efails", set: ARKodeSetSmallNumEFails },
        sunKeyIntPair { key: "max_err_test_fails", set: ARKodeSetMaxErrTestFails },
        sunKeyIntPair { key: "max_conv_fails", set: ARKodeSetMaxConvFails },
        sunKeyIntPair { key: "linear_solution_scaling", set: cli_set_linear_solution_scaling },
        sunKeyIntPair { key: "use_compensated_sums", set: cli_set_use_compensated_sums },
    ];

    let long_pairs: [sunKeyLongPair<ARKodeMem>; 2] = [
        sunKeyLongPair { key: "max_num_steps", set: ARKodeSetMaxNumSteps },
        sunKeyLongPair { key: "jac_eval_frequency", set: ARKodeSetJacEvalFrequency },
    ];

    let real_pairs: [sunKeyRealPair<ARKodeMem>; 22] = [
        sunKeyRealPair { key: "nonlin_crdown", set: ARKodeSetNonlinCRDown },
        sunKeyRealPair { key: "nonlin_rdiv", set: ARKodeSetNonlinRDiv },
        sunKeyRealPair { key: "delta_gamma_max", set: ARKodeSetDeltaGammaMax },
        sunKeyRealPair { key: "nonlin_conv_coef", set: ARKodeSetNonlinConvCoef },
        sunKeyRealPair { key: "init_step", set: ARKodeSetInitStep },
        sunKeyRealPair { key: "min_step", set: ARKodeSetMinStep },
        sunKeyRealPair { key: "max_step", set: ARKodeSetMaxStep },
        sunKeyRealPair { key: "stop_time", set: ARKodeSetStopTime },
        sunKeyRealPair { key: "fixed_step", set: ARKodeSetFixedStep },
        sunKeyRealPair { key: "step_direction", set: ARKodeSetStepDirection },
        sunKeyRealPair { key: "cfl_fraction", set: ARKodeSetCFLFraction },
        sunKeyRealPair { key: "safety_factor", set: ARKodeSetSafetyFactor },
        sunKeyRealPair { key: "error_bias", set: ARKodeSetErrorBias },
        sunKeyRealPair { key: "max_growth", set: ARKodeSetMaxGrowth },
        sunKeyRealPair { key: "min_reduction", set: ARKodeSetMinReduction },
        sunKeyRealPair { key: "max_first_growth", set: ARKodeSetMaxFirstGrowth },
        sunKeyRealPair { key: "max_efail_growth", set: ARKodeSetMaxEFailGrowth },
        sunKeyRealPair { key: "max_cfail_growth", set: ARKodeSetMaxCFailGrowth },
        sunKeyRealPair { key: "eps_lin", set: ARKodeSetEpsLin },
        sunKeyRealPair { key: "mass_eps_lin", set: ARKodeSetMassEpsLin },
        sunKeyRealPair { key: "ls_norm_factor", set: ARKodeSetLSNormFactor },
        sunKeyRealPair { key: "mass_ls_norm_factor", set: ARKodeSetMassLSNormFactor },
    ];

    let tworeal_pairs: [sunKeyTwoRealPair<ARKodeMem>; 2] = [
        sunKeyTwoRealPair { key: "scalar_tolerances", set: crate::arkode::ARKodeSStolerances },
        sunKeyTwoRealPair { key: "fixed_step_bounds", set: ARKodeSetFixedStepBounds },
    ];

    let action_pairs: [sunKeyActionPair<ARKodeMem>; 4] = [
        sunKeyActionPair { key: "nonlinear", set: ARKodeSetNonlinear },
        sunKeyActionPair { key: "clear_stop_time", set: ARKodeClearStopTime },
        sunKeyActionPair { key: "no_inactive_root_warn", set: ARKodeSetNoInactiveRootWarn },
        sunKeyActionPair { key: "reset_accumulated_error", set: ARKodeResetAccumulatedError },
    ];

    /* Prefix for options to set */
    let default_id = "arkode";
    let prefix = match arkid {
        Some(id) if !id.is_empty() => format!("{}.", id),
        _ => format!("{}.", default_id),
    };
    let offset = prefix.len();

    let mut write_parameters = false;
    let mut idx: usize = 1;
    while idx < argv.len() {
        let mut arg_used = false;
        let mut j: usize = 0;

        /* skip command-line arguments that do not begin with correct prefix */
        if !argv[idx].starts_with(&prefix) {
            idx += 1;
            continue;
        }

        /* check all "int" command-line options */
        let retval = sunCheckAndSetIntArgs(
            ark_mem,
            &mut idx,
            argv,
            offset,
            &int_pairs,
            &mut arg_used,
            &mut j,
        );
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                retval,
                line!(),
                "arkSetFromCommandLine",
                file!(),
                &format!("error setting key: {}", int_pairs[j].key),
            );
            return retval;
        }
        if arg_used {
            idx += 1;
            continue;
        }

        /* check all long int command-line options */
        let retval = sunCheckAndSetLongArgs(
            ark_mem,
            &mut idx,
            argv,
            offset,
            &long_pairs,
            &mut arg_used,
            &mut j,
        );
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                retval,
                line!(),
                "arkSetFromCommandLine",
                file!(),
                &format!("error setting key: {}", long_pairs[j].key),
            );
            return retval;
        }
        if arg_used {
            idx += 1;
            continue;
        }

        /* check all real command-line options */
        let retval = sunCheckAndSetRealArgs(
            ark_mem,
            &mut idx,
            argv,
            offset,
            &real_pairs,
            &mut arg_used,
            &mut j,
        );
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                retval,
                line!(),
                "arkSetFromCommandLine",
                file!(),
                &format!("error setting key: {}", real_pairs[j].key),
            );
            return retval;
        }
        if arg_used {
            idx += 1;
            continue;
        }

        /* check all pair-of-real command-line options */
        let retval = sunCheckAndSetTwoRealArgs(
            ark_mem,
            &mut idx,
            argv,
            offset,
            &tworeal_pairs,
            &mut arg_used,
            &mut j,
        );
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                retval,
                line!(),
                "arkSetFromCommandLine",
                file!(),
                &format!("error setting key: {}", tworeal_pairs[j].key),
            );
            return retval;
        }
        if arg_used {
            idx += 1;
            continue;
        }

        /* check all action command-line options */
        let retval = sunCheckAndSetActionArgs(
            ark_mem,
            &mut idx,
            argv,
            offset,
            &action_pairs,
            &mut arg_used,
            &mut j,
        );
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                retval,
                line!(),
                "arkSetFromCommandLine",
                file!(),
                &format!("error setting key: {}", action_pairs[j].key),
            );
            return retval;
        }
        if arg_used {
            idx += 1;
            continue;
        }

        /*** handle all remaining command-line options ***/

        if &argv[idx][offset..] == "interpolant_type" {
            idx += 1;
            let mut retval = ARK_ILL_INPUT;
            if argv[idx] == "ARK_INTERP_HERMITE" {
                retval = ARKodeSetInterpolantType(ark_mem, ARK_INTERP_HERMITE);
            } else if argv[idx] == "ARK_INTERP_LAGRANGE" {
                retval = ARKodeSetInterpolantType(ark_mem, ARK_INTERP_LAGRANGE);
            } else if argv[idx] == "ARK_INTERP_NONE" {
                retval = ARKodeSetInterpolantType(ark_mem, ARK_INTERP_NONE);
            }
            if retval != ARK_SUCCESS {
                arkProcessError(
                    Some(ark_mem),
                    retval,
                    line!(),
                    "arkSetFromCommandLine",
                    file!(),
                    &format!("error setting key: {} {}", argv[idx - 1], argv[idx]),
                );
                return retval;
            }
            idx += 1;
            continue;
        }

        if &argv[idx][offset..] == "accumulated_error_type" {
            idx += 1;
            let mut retval = ARK_ILL_INPUT;
            if argv[idx] == "ARK_ACCUMERROR_NONE" {
                retval = ARKodeSetAccumulatedErrorType(ark_mem, ARK_ACCUMERROR_NONE);
            } else if argv[idx] == "ARK_ACCUMERROR_MAX" {
                retval = ARKodeSetAccumulatedErrorType(ark_mem, ARK_ACCUMERROR_MAX);
            } else if argv[idx] == "ARK_ACCUMERROR_SUM" {
                retval = ARKodeSetAccumulatedErrorType(ark_mem, ARK_ACCUMERROR_SUM);
            } else if argv[idx] == "ARK_ACCUMERROR_AVG" {
                retval = ARKodeSetAccumulatedErrorType(ark_mem, ARK_ACCUMERROR_AVG);
            }
            if retval != ARK_SUCCESS {
                arkProcessError(
                    Some(ark_mem),
                    retval,
                    line!(),
                    "arkSetFromCommandLine",
                    file!(),
                    &format!("error setting key: {} {}", argv[idx - 1], argv[idx]),
                );
                return retval;
            }
            idx += 1;
            continue;
        }

        if &argv[idx][offset..] == "write_parameters" {
            write_parameters = true;
            idx += 1;
            continue;
        }

        /* Call stepper-specific SetFromCommandLine routine (if supplied) to
           process this command-line argument */
        if let Some(step_setoptions) = ark_mem.step_setoptions {
            let retval = step_setoptions(ark_mem, &mut idx, argv, offset, &mut arg_used);
            if retval != ARK_SUCCESS {
                return retval;
            }
            if arg_used {
                idx += 1;
                continue;
            }
        }

        /* warn for uninterpreted arkid.X arguments */
        arkProcessError(
            Some(ark_mem),
            ARK_WARNING,
            line!(),
            "arkSetFromCommandLine",
            file!(),
            &format!("WARNING: key {} was not handled\n", argv[idx]),
        );
        idx += 1;
    }

    /* Call ARKodeWriteParameters (if requested) now that all
       command-line options have been set */
    if write_parameters {
        let mut stdout = std::io::stdout();
        let retval = ARKodeWriteParameters(ark_mem, &mut stdout);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                retval,
                line!(),
                "arkSetFromCommandLine",
                file!(),
                "error writing parameters to stdout",
            );
            return retval;
        }
    }

    ARK_SUCCESS
}
