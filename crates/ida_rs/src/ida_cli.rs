/* ---------------------------------------------------------------
 * Translated from src/ida/ida_cli.c (IDA 7.7.0).
 *
 * Provides command-line control over optional inputs to IDA:
 * arguments of the form "<idaid>.<key>" followed by value
 * argument(s) are mapped onto the corresponding IDASet* calls.
 * The key names are identical to the C dispatch tables; the
 * sunCheckAndSet* helpers come from sundials_core::sundials_cli
 * (the shared translation of src/sundials/sundials_cli.c).
 * ---------------------------------------------------------------*/
use crate::ida::IDASStolerances;
use crate::ida_impl::*;
use crate::ida_io::*;
use crate::ida_ls::{
    IDASetEpsLin, IDASetIncrementFactor, IDASetLSNormFactor, IDASetLinearSolutionScaling,
};
use sundials_core::sundials_cli::{
    sunCheckAndSetActionArgs, sunCheckAndSetIntArgs, sunCheckAndSetLongArgs, sunCheckAndSetRealArgs,
    sunCheckAndSetTwoRealArgs, sunKeyActionPair, sunKeyIntPair, sunKeyLongPair, sunKeyRealPair,
    sunKeyTwoRealPair,
};

/*---------------------------------------------------------------
  IDASetOptions:

  Sets IDA options using strings.
  ---------------------------------------------------------------*/
pub fn IDASetOptions(ida_mem: &mut IDAMem, idaid: &str, file_name: &str, args: &[String]) -> i32 {
    if !file_name.is_empty() {
        let retval = IDA_ILL_INPUT;
        IDAProcessError(Some(ida_mem), retval, line!(), "IDASetOptions", file!(),
                        "file-based options are not currently supported.");
        return retval;
    }

    if !args.is_empty() {
        let retval = idaSetFromCommandLine(ida_mem, idaid, args);
        if retval != IDA_SUCCESS {
            return retval;
        }
    }

    IDA_SUCCESS
}

/* -----------------------------------------------------------------
 * Function to control IDA options from the command line
 */
fn idaSetFromCommandLine(ida_mem: &mut IDAMem, idaid: &str, argv: &[String]) -> i32 {
    /* Set lists of command-line arguments, and the corresponding set routines */
    let int_pairs: &[sunKeyIntPair<IDAMem>] = &[
        sunKeyIntPair { key: "max_num_steps_ic", set: IDASetMaxNumStepsIC },
        sunKeyIntPair { key: "max_num_jacs_ic", set: IDASetMaxNumJacsIC },
        sunKeyIntPair { key: "max_num_iters_ic", set: IDASetMaxNumItersIC },
        sunKeyIntPair { key: "line_search_off_ic", set: |m, v| IDASetLineSearchOffIC(m, v != 0) },
        sunKeyIntPair { key: "max_backs_ic", set: IDASetMaxBacksIC },
        sunKeyIntPair { key: "max_order", set: IDASetMaxOrd },
        sunKeyIntPair { key: "max_err_test_fails", set: IDASetMaxErrTestFails },
        sunKeyIntPair { key: "suppress_alg", set: |m, v| IDASetSuppressAlg(m, v != 0) },
        sunKeyIntPair { key: "max_conv_fails", set: IDASetMaxConvFails },
        sunKeyIntPair { key: "max_nonlin_iters", set: IDASetMaxNonlinIters },
        sunKeyIntPair {
            key: "linear_solution_scaling",
            set: |m, v| IDASetLinearSolutionScaling(m, v != 0),
        },
        sunKeyIntPair { key: "max_num_constraint_fails", set: IDASetMaxNumConstraintFails },
    ];

    let long_pairs: &[sunKeyLongPair<IDAMem>] =
        &[sunKeyLongPair { key: "max_num_steps", set: IDASetMaxNumSteps }];

    let real_pairs: &[sunKeyRealPair<IDAMem>] = &[
        sunKeyRealPair { key: "nonlin_conv_coef_ic", set: IDASetNonlinConvCoefIC },
        sunKeyRealPair { key: "step_tolerance_ic", set: IDASetStepToleranceIC },
        sunKeyRealPair { key: "delta_cj_lsetup", set: IDASetDeltaCjLSetup },
        sunKeyRealPair { key: "init_step", set: IDASetInitStep },
        sunKeyRealPair { key: "max_step", set: IDASetMaxStep },
        sunKeyRealPair { key: "min_step", set: IDASetMinStep },
        sunKeyRealPair { key: "stop_time", set: IDASetStopTime },
        sunKeyRealPair { key: "eta_min", set: IDASetEtaMin },
        sunKeyRealPair { key: "eta_max", set: IDASetEtaMax },
        sunKeyRealPair { key: "eta_low", set: IDASetEtaLow },
        sunKeyRealPair { key: "eta_min_err_fail", set: IDASetEtaMinErrFail },
        sunKeyRealPair { key: "eta_conv_fail", set: IDASetEtaConvFail },
        sunKeyRealPair { key: "nonlin_conv_coef", set: IDASetNonlinConvCoef },
        sunKeyRealPair { key: "eps_lin", set: IDASetEpsLin },
        sunKeyRealPair { key: "ls_norm_factor", set: IDASetLSNormFactor },
        sunKeyRealPair { key: "increment_factor", set: IDASetIncrementFactor },
    ];

    let tworeal_pairs: &[sunKeyTwoRealPair<IDAMem>] = &[
        sunKeyTwoRealPair { key: "eta_fixed_step_bounds", set: IDASetEtaFixedStepBounds },
        sunKeyTwoRealPair { key: "scalar_tolerances", set: IDASStolerances },
    ];

    let action_pairs: &[sunKeyActionPair<IDAMem>] = &[
        sunKeyActionPair { key: "clear_stop_time", set: IDAClearStopTime },
        sunKeyActionPair { key: "no_inactive_root_warn", set: IDASetNoInactiveRootWarn },
    ];

    /* Prefix for options to set */
    let default_id = "ida";
    let prefix = if !idaid.is_empty() {
        format!("{}.", idaid)
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
            let retval = sunCheckAndSetIntArgs(ida_mem, &mut idx, argv, offset, int_pairs,
                                               &mut arg_used, &mut j);
            if retval != IDA_SUCCESS {
                IDAProcessError(Some(ida_mem), retval, line!(), "idaSetFromCommandLine", file!(),
                                &format!("error setting key: {}", int_pairs[j].key));
                return retval;
            }
            if arg_used {
                break 'this_arg;
            }

            /* check all long int command-line options */
            let retval = sunCheckAndSetLongArgs(ida_mem, &mut idx, argv, offset, long_pairs,
                                                &mut arg_used, &mut j);
            if retval != IDA_SUCCESS {
                IDAProcessError(Some(ida_mem), retval, line!(), "idaSetFromCommandLine", file!(),
                                &format!("error setting key: {}", long_pairs[j].key));
                return retval;
            }
            if arg_used {
                break 'this_arg;
            }

            /* check all real command-line options */
            let retval = sunCheckAndSetRealArgs(ida_mem, &mut idx, argv, offset, real_pairs,
                                                &mut arg_used, &mut j);
            if retval != IDA_SUCCESS {
                IDAProcessError(Some(ida_mem), retval, line!(), "idaSetFromCommandLine", file!(),
                                &format!("error setting key: {}", real_pairs[j].key));
                return retval;
            }
            if arg_used {
                break 'this_arg;
            }

            /* check all pair-of-real command-line options */
            let retval = sunCheckAndSetTwoRealArgs(ida_mem, &mut idx, argv, offset, tworeal_pairs,
                                                   &mut arg_used, &mut j);
            if retval != IDA_SUCCESS {
                IDAProcessError(Some(ida_mem), retval, line!(), "idaSetFromCommandLine", file!(),
                                &format!("error setting key: {}", tworeal_pairs[j].key));
                return retval;
            }
            if arg_used {
                break 'this_arg;
            }

            /* check all action command-line options */
            let retval = sunCheckAndSetActionArgs(ida_mem, &mut idx, argv, offset, action_pairs,
                                                  &mut arg_used, &mut j);
            if retval != IDA_SUCCESS {
                IDAProcessError(Some(ida_mem), retval, line!(), "idaSetFromCommandLine", file!(),
                                &format!("error setting key: {}", action_pairs[j].key));
                return retval;
            }
            if arg_used {
                break 'this_arg;
            }

            /* warn for uninterpreted idaid.X arguments */
            IDAProcessError(Some(ida_mem), IDA_WARNING, line!(), "idaSetFromCommandLine", file!(),
                            &format!("WARNING: key {} was not handled\n", argv[idx]));
        }
        idx += 1;
    }

    IDA_SUCCESS
}

/*==================================================================
  Tests
  ==================================================================*/
#[cfg(test)]
mod tests {
    use super::*;

    fn argv(items: &[&str]) -> Vec<String> {
        items.iter().map(|s| s.to_string()).collect()
    }

    /* options with the default "ida." prefix reach the IDASet* routines;
       foreign prefixes are skipped (argv[0] is the program name). */
    #[test]
    fn idasetoptions_parses_int_long_real_and_tworeal_keys() {
        let mut ida_mem = IDAMem::default();
        let args = argv(&[
            "prog",
            "ida.max_num_steps", "1234",
            "other.max_num_steps", "5",
            "ida.max_order", "3",
            "ida.suppress_alg", "1",
            "ida.init_step", "0.25",
            "ida.eta_fixed_step_bounds", "0.5", "1.5",
        ]);
        assert_eq!(IDASetOptions(&mut ida_mem, "", "", &args), IDA_SUCCESS);
        assert_eq!(ida_mem.ida_mxstep, 1234); /* "other." prefix skipped */
        assert_eq!(ida_mem.ida_maxord, 3);
        assert!(ida_mem.ida_suppressalg);
        assert_eq!(ida_mem.ida_hin, 0.25);
        assert_eq!(ida_mem.ida_eta_min_fx, 0.5);
        assert_eq!(ida_mem.ida_eta_max_fx, 1.5);
    }

    /* a custom idaid replaces the default prefix */
    #[test]
    fn idasetoptions_honors_custom_idaid() {
        let mut ida_mem = IDAMem::default();
        let args = argv(&[
            "prog",
            "mysolver.max_conv_fails", "9",
            "ida.max_conv_fails", "40", /* wrong prefix now: skipped */
        ]);
        assert_eq!(IDASetOptions(&mut ida_mem, "mysolver", "", &args), IDA_SUCCESS);
        assert_eq!(ida_mem.ida_maxncf, 9);
    }

    /* action keys take no value argument */
    #[test]
    fn idasetoptions_handles_action_keys() {
        let mut ida_mem = IDAMem::default();
        ida_mem.ida_tstopset = true;
        let args = argv(&["prog", "ida.clear_stop_time"]);
        assert_eq!(IDASetOptions(&mut ida_mem, "", "", &args), IDA_SUCCESS);
        assert!(!ida_mem.ida_tstopset);
    }

    /* file-based options are rejected (IDA_ILL_INPUT) */
    #[test]
    fn idasetoptions_rejects_file_input() {
        let mut ida_mem = IDAMem::default();
        assert_eq!(IDASetOptions(&mut ida_mem, "", "opts.txt", &[]), IDA_ILL_INPUT);
    }

    /* unhandled "ida.*" keys only warn; later options still apply */
    #[test]
    fn idasetoptions_warns_and_continues_on_unknown_key() {
        let mut ida_mem = IDAMem::default();
        let args = argv(&["prog", "ida.bogus_key", "ida.max_conv_fails", "12"]);
        assert_eq!(IDASetOptions(&mut ida_mem, "", "", &args), IDA_SUCCESS);
        assert_eq!(ida_mem.ida_maxncf, 12);
    }
}
