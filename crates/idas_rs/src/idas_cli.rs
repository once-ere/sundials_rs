/* ---------------------------------------------------------------
 * Translated from src/idas/idas_cli.c (IDAS 7.7.0).
 *
 * Provides command-line control over optional inputs to IDAS:
 * arguments of the form "<idaid>.<key>" followed by value
 * argument(s) are mapped onto the corresponding IDASet* calls.
 * The key names are identical to the C dispatch tables; the
 * sunCheckAndSet* helpers come from sundials_core::sundials_cli
 * (the shared translation of src/sundials/sundials_cli.c).
 * Relative to the IDA donor (ida_cli.rs) this adds the quadrature /
 * sensitivity keys and the backward-problem (int-prefixed) tables.
 * ---------------------------------------------------------------*/
use crate::idaa::{IDAQuadSStolerancesB, IDASStolerancesB};
use crate::idaa_io::{
    IDAAdjSetNoSensi, IDASetInitStepB, IDASetMaxNumStepsB, IDASetMaxOrdB, IDASetMaxStepB,
    IDASetQuadErrConB, IDASetSuppressAlgB,
};
use crate::idas::{IDAQuadSStolerances, IDASStolerances, IDASensToggleOff};
use crate::idas_impl::*;
use crate::idas_io::*;
use crate::idas_ls::{
    IDASetEpsLin, IDASetEpsLinB, IDASetIncrementFactor, IDASetIncrementFactorB,
    IDASetLSNormFactor, IDASetLSNormFactorB, IDASetLinearSolutionScaling,
    IDASetLinearSolutionScalingB,
};
use sundials_core::sundials_cli::{
    sunCheckAndSetActionArgs, sunCheckAndSetIntArgs, sunCheckAndSetIntLongArgs,
    sunCheckAndSetIntRealArgs, sunCheckAndSetIntRealRealArgs, sunCheckAndSetLongArgs,
    sunCheckAndSetRealArgs, sunCheckAndSetTwoIntArgs, sunCheckAndSetTwoRealArgs,
    sunKeyActionPair, sunKeyIntLongPair, sunKeyIntPair, sunKeyIntRealPair, sunKeyIntRealRealPair,
    sunKeyLongPair, sunKeyRealPair, sunKeyTwoIntPair, sunKeyTwoRealPair,
};

/*---------------------------------------------------------------
  IDASetOptions:

  Sets IDAS options using strings.
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
 * Function to control IDAS options from the command line
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
        sunKeyIntPair { key: "quad_err_con", set: |m, v| IDASetQuadErrCon(m, v != 0) },
        sunKeyIntPair { key: "sens_err_con", set: |m, v| IDASetSensErrCon(m, v != 0) },
        sunKeyIntPair { key: "sens_max_nonlin_iters", set: IDASetSensMaxNonlinIters },
        sunKeyIntPair { key: "quad_sens_err_con", set: |m, v| IDASetQuadSensErrCon(m, v != 0) },
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
        sunKeyTwoRealPair { key: "quad_scalar_tolerances", set: IDAQuadSStolerances },
    ];

    let twoint_pairs: &[sunKeyTwoIntPair<IDAMem>] = &[
        sunKeyTwoIntPair { key: "max_order_b", set: IDASetMaxOrdB },
        sunKeyTwoIntPair { key: "suppress_alg_b", set: |m, w, v| IDASetSuppressAlgB(m, w, v != 0) },
        sunKeyTwoIntPair { key: "quad_err_con_b", set: |m, w, v| IDASetQuadErrConB(m, w, v != 0) },
        sunKeyTwoIntPair {
            key: "linear_solution_scaling_b",
            set: |m, w, v| IDASetLinearSolutionScalingB(m, w, v != 0),
        },
    ];

    let action_pairs: &[sunKeyActionPair<IDAMem>] = &[
        sunKeyActionPair { key: "clear_stop_time", set: IDAClearStopTime },
        sunKeyActionPair { key: "no_inactive_root_warn", set: IDASetNoInactiveRootWarn },
        sunKeyActionPair { key: "sens_toggle_off", set: IDASensToggleOff },
        sunKeyActionPair { key: "adj_no_sensi", set: IDAAdjSetNoSensi },
    ];

    let int_real_pairs: &[sunKeyIntRealPair<IDAMem>] = &[
        sunKeyIntRealPair { key: "sens_dq_method", set: IDASetSensDQMethod },
        sunKeyIntRealPair { key: "init_step_b", set: IDASetInitStepB },
        sunKeyIntRealPair { key: "max_step_b", set: IDASetMaxStepB },
        sunKeyIntRealPair { key: "eps_lin_b", set: IDASetEpsLinB },
        sunKeyIntRealPair { key: "ls_norm_factor_b", set: IDASetLSNormFactorB },
        sunKeyIntRealPair { key: "increment_factor_b", set: IDASetIncrementFactorB },
    ];

    let int_real_real_pairs: &[sunKeyIntRealRealPair<IDAMem>] = &[
        sunKeyIntRealRealPair { key: "scalar_tolerances_b", set: IDASStolerancesB },
        sunKeyIntRealRealPair { key: "quad_scalar_tolerances_b", set: IDAQuadSStolerancesB },
    ];

    let int_long_pairs: &[sunKeyIntLongPair<IDAMem>] =
        &[sunKeyIntLongPair { key: "max_num_steps_b", set: IDASetMaxNumStepsB }];

    /* Prefix for options to set */
    let default_id = "idas";
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

            /* check all pair-of-int command-line options */
            let retval = sunCheckAndSetTwoIntArgs(ida_mem, &mut idx, argv, offset, twoint_pairs,
                                                  &mut arg_used, &mut j);
            if retval != IDA_SUCCESS {
                IDAProcessError(Some(ida_mem), retval, line!(), "idaSetFromCommandLine", file!(),
                                &format!("error setting key: {}", twoint_pairs[j].key));
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

            /* check all int+real command-line options */
            let retval = sunCheckAndSetIntRealArgs(ida_mem, &mut idx, argv, offset,
                                                   int_real_pairs, &mut arg_used, &mut j);
            if retval != IDA_SUCCESS {
                IDAProcessError(Some(ida_mem), retval, line!(), "idaSetFromCommandLine", file!(),
                                &format!("error setting key: {}", int_real_pairs[j].key));
                return retval;
            }
            if arg_used {
                break 'this_arg;
            }

            /* check all int+long command-line options */
            let retval = sunCheckAndSetIntLongArgs(ida_mem, &mut idx, argv, offset,
                                                   int_long_pairs, &mut arg_used, &mut j);
            if retval != IDA_SUCCESS {
                IDAProcessError(Some(ida_mem), retval, line!(), "idaSetFromCommandLine", file!(),
                                &format!("error setting key: {}", int_long_pairs[j].key));
                return retval;
            }
            if arg_used {
                break 'this_arg;
            }

            /* check all int+real+real command-line options */
            let retval = sunCheckAndSetIntRealRealArgs(ida_mem, &mut idx, argv, offset,
                                                       int_real_real_pairs, &mut arg_used,
                                                       &mut j);
            if retval != IDA_SUCCESS {
                IDAProcessError(Some(ida_mem), retval, line!(), "idaSetFromCommandLine", file!(),
                                &format!("error setting key: {}", int_real_real_pairs[j].key));
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
    use crate::idaa::{IDAAdjInit, IDACreateB, IDAInitB};
    use crate::nvector_serial::NVector;
    use crate::sundials_types::UserData;

    fn argv(items: &[&str]) -> Vec<String> {
        items.iter().map(|s| s.to_string()).collect()
    }

    /* options with the default "idas." prefix reach the IDASet* routines;
       foreign prefixes are skipped (argv[0] is the program name). */
    #[test]
    fn idasetoptions_parses_int_long_real_and_tworeal_keys() {
        let mut ida_mem = IDAMem::default();
        /* quad keys require initialized quadratures (IDA_NO_QUAD guard) */
        ida_mem.ida_quadMallocDone = true;
        let args = argv(&[
            "prog",
            "idas.max_num_steps", "1234",
            "other.max_num_steps", "5",
            "idas.max_order", "3",
            "idas.suppress_alg", "1",
            "idas.quad_err_con", "1",
            "idas.sens_err_con", "1",
            "idas.init_step", "0.25",
            "idas.eta_fixed_step_bounds", "0.5", "1.5",
            "idas.quad_scalar_tolerances", "1e-4", "1e-8",
        ]);
        assert_eq!(IDASetOptions(&mut ida_mem, "", "", &args), IDA_SUCCESS);
        assert_eq!(ida_mem.ida_mxstep, 1234); /* "other." prefix skipped */
        assert_eq!(ida_mem.ida_maxord, 3);
        assert!(ida_mem.ida_suppressalg);
        assert!(ida_mem.ida_errconQ);
        assert!(ida_mem.ida_errconS);
        assert_eq!(ida_mem.ida_hin, 0.25);
        assert_eq!(ida_mem.ida_eta_min_fx, 0.5);
        assert_eq!(ida_mem.ida_eta_max_fx, 1.5);
        assert_eq!(ida_mem.ida_rtolQ, 1e-4);
        assert_eq!(ida_mem.ida_SatolQ, 1e-8);
    }

    /* backward-problem keys (which-prefixed tables) reach the nested
       solver through the ***B wrappers */
    #[test]
    fn idasetoptions_parses_backward_problem_keys() {
        fn resB(_tt: f64, _yy: &NVector, _yp: &NVector, _yyB: &NVector, _ypB: &NVector,
                _rrB: &mut NVector, _ud: &mut UserData) -> i32 {
            0
        }

        let mut ida_mem = IDAMem::default();
        ida_mem.ida_tempv1 = NVector::new(2);
        assert_eq!(IDAAdjInit(&mut ida_mem, 5, IDA_HERMITE), IDA_SUCCESS);
        let mut which = -1;
        assert_eq!(IDACreateB(&mut ida_mem, &mut which), IDA_SUCCESS);
        let yyB0 = NVector::from_slice(&[1.0, 2.0]);
        let ypB0 = NVector::from_slice(&[0.0, 0.0]);
        assert_eq!(IDAInitB(&mut ida_mem, which, resB, 0.0, &yyB0, &ypB0), IDA_SUCCESS);

        let args = argv(&[
            "prog",
            "idas.max_order_b", "0", "2",
            "idas.max_num_steps_b", "0", "321",
            "idas.init_step_b", "0", "0.125",
            "idas.scalar_tolerances_b", "0", "1e-5", "1e-9",
            "idas.adj_no_sensi",
        ]);
        assert_eq!(IDASetOptions(&mut ida_mem, "", "", &args), IDA_SUCCESS);
        {
            let idaadj_mem = ida_mem.ida_adj_mem.as_ref().unwrap();
            assert!(!idaadj_mem.ia_storeSensi); /* adj_no_sensi */
            let nested = &idaadj_mem.IDAB_mem[0].IDA_mem;
            assert_eq!(nested.ida_maxord, 2);
            assert_eq!(nested.ida_mxstep, 321);
            assert_eq!(nested.ida_hin, 0.125);
            assert_eq!(nested.ida_rtol, 1e-5);
            assert_eq!(nested.ida_Satol, 1e-9);
        }
    }

    /* a custom idaid replaces the default prefix */
    #[test]
    fn idasetoptions_honors_custom_idaid() {
        let mut ida_mem = IDAMem::default();
        let args = argv(&[
            "prog",
            "mysolver.max_conv_fails", "9",
            "idas.max_conv_fails", "40", /* wrong prefix now: skipped */
        ]);
        assert_eq!(IDASetOptions(&mut ida_mem, "mysolver", "", &args), IDA_SUCCESS);
        assert_eq!(ida_mem.ida_maxncf, 9);
    }

    /* file-based options are rejected (IDA_ILL_INPUT) */
    #[test]
    fn idasetoptions_rejects_file_input() {
        let mut ida_mem = IDAMem::default();
        assert_eq!(IDASetOptions(&mut ida_mem, "", "opts.txt", &[]), IDA_ILL_INPUT);
    }

    /* unhandled "idas.*" keys only warn; later options still apply */
    #[test]
    fn idasetoptions_warns_and_continues_on_unknown_key() {
        let mut ida_mem = IDAMem::default();
        let args = argv(&["prog", "idas.bogus_key", "idas.max_conv_fails", "12"]);
        assert_eq!(IDASetOptions(&mut ida_mem, "", "", &args), IDA_SUCCESS);
        assert_eq!(ida_mem.ida_maxncf, 12);
    }
}
