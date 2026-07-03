/* ---------------------------------------------------------------
 * Translated from src/kinsol/kinsol_cli.c (KINSOL 7.7.0).
 *
 * Provides command-line control over optional inputs to KINSOL:
 * arguments of the form "<kinid>.<key>" followed by value
 * argument(s) are mapped onto the corresponding KINSet* calls.
 * The key names are identical to the C dispatch tables; the
 * sunCheckAndSet* helpers come from sundials_core::sundials_cli
 * (the shared translation of src/sundials/sundials_cli.c).
 * ---------------------------------------------------------------*/
use crate::kinsol_impl::*;
use crate::kinsol_io::*;
use sundials_core::sundials_cli::{
    sunCheckAndSetIntArgs, sunCheckAndSetLongArgs, sunCheckAndSetRealArgs,
    sunCheckAndSetTwoRealArgs, sunKeyIntPair, sunKeyLongPair, sunKeyRealPair, sunKeyTwoRealPair,
};

/*---------------------------------------------------------------
  KINSetOptions:

  Sets KINSOL options using strings.
  ---------------------------------------------------------------*/

pub fn KINSetOptions(kin_mem: &mut KINMem, kinid: &str, file_name: &str, args: &[String]) -> i32 {
    if !file_name.is_empty() {
        let retval = KIN_ILL_INPUT;
        KINProcessError(Some(kin_mem), retval, line!(), "KINSetOptions", file!(),
                        "file-based options are not currently supported.");
        return retval;
    }

    if !args.is_empty() {
        let retval = kinSetFromCommandLine(kin_mem, kinid, args);
        if retval != KIN_SUCCESS {
            return retval;
        }
    }

    KIN_SUCCESS
}

/* -----------------------------------------------------------------
 * Function to control KINSOL options from the command line
 */

fn kinSetFromCommandLine(kin_mem: &mut KINMem, kinid: &str, argv: &[String]) -> i32 {
    /* (The C kinmem NULL check with MSG_NO_MEM cannot arise: kinmem is
       a &mut KINMem here.) */

    /* Set lists of command-line arguments, and the corresponding set routines */
    let int_pairs: &[sunKeyIntPair<KINMem>] = &[
        sunKeyIntPair { key: "orth_aa", set: KINSetOrthAA },
        sunKeyIntPair { key: "return_newest", set: |m, v| KINSetReturnNewest(m, v != 0) },
        sunKeyIntPair { key: "no_init_setup", set: |m, v| KINSetNoInitSetup(m, v != 0) },
        sunKeyIntPair { key: "no_res_mon", set: |m, v| KINSetNoResMon(m, v != 0) },
        sunKeyIntPair { key: "eta_form", set: KINSetEtaForm },
        sunKeyIntPair { key: "no_min_eps", set: |m, v| KINSetNoMinEps(m, v != 0) },
    ];

    let long_pairs: &[sunKeyLongPair<KINMem>] = &[
        sunKeyLongPair { key: "m_aa", set: KINSetMAA },
        sunKeyLongPair { key: "delay_aa", set: KINSetDelayAA },
        sunKeyLongPair { key: "num_max_iters", set: KINSetNumMaxIters },
        sunKeyLongPair { key: "max_setup_calls", set: KINSetMaxSetupCalls },
        sunKeyLongPair { key: "max_sub_setup_calls", set: KINSetMaxSubSetupCalls },
        sunKeyLongPair { key: "max_beta_fails", set: KINSetMaxBetaFails },
    ];

    let real_pairs: &[sunKeyRealPair<KINMem>] = &[
        sunKeyRealPair { key: "damping", set: KINSetDamping },
        sunKeyRealPair { key: "damping_aa", set: KINSetDampingAA },
        sunKeyRealPair { key: "eta_const_value", set: KINSetEtaConstValue },
        sunKeyRealPair { key: "res_mon_const_value", set: KINSetResMonConstValue },
        sunKeyRealPair { key: "max_newton_step", set: KINSetMaxNewtonStep },
        sunKeyRealPair { key: "rel_err_func", set: KINSetRelErrFunc },
        sunKeyRealPair { key: "func_norm_tol", set: KINSetFuncNormTol },
        sunKeyRealPair { key: "scaled_step_tol", set: KINSetScaledStepTol },
    ];

    let tworeal_pairs: &[sunKeyTwoRealPair<KINMem>] = &[
        sunKeyTwoRealPair { key: "eta_params", set: KINSetEtaParams },
        sunKeyTwoRealPair { key: "res_mon_params", set: KINSetResMonParams },
    ];

    /* Prefix for options to set */
    let default_id = "kinsol";
    let prefix = if !kinid.is_empty() {
        format!("{}.", kinid)
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
            let retval = sunCheckAndSetIntArgs(kin_mem, &mut idx, argv, offset, int_pairs,
                                               &mut arg_used, &mut j);
            if retval != KIN_SUCCESS {
                KINProcessError(Some(kin_mem), retval, line!(), "kinSetFromCommandLine", file!(),
                                &format!("error setting key: {}", int_pairs[j].key));
                return retval;
            }
            if arg_used {
                break 'this_arg;
            }

            /* check all long int command-line options */
            let retval = sunCheckAndSetLongArgs(kin_mem, &mut idx, argv, offset, long_pairs,
                                                &mut arg_used, &mut j);
            if retval != KIN_SUCCESS {
                KINProcessError(Some(kin_mem), retval, line!(), "kinSetFromCommandLine", file!(),
                                &format!("error setting key: {}", long_pairs[j].key));
                return retval;
            }
            if arg_used {
                break 'this_arg;
            }

            /* check all real command-line options */
            let retval = sunCheckAndSetRealArgs(kin_mem, &mut idx, argv, offset, real_pairs,
                                                &mut arg_used, &mut j);
            if retval != KIN_SUCCESS {
                KINProcessError(Some(kin_mem), retval, line!(), "kinSetFromCommandLine", file!(),
                                &format!("error setting key: {}", real_pairs[j].key));
                return retval;
            }
            if arg_used {
                break 'this_arg;
            }

            /* check all pair-of-real command-line options */
            let retval = sunCheckAndSetTwoRealArgs(kin_mem, &mut idx, argv, offset,
                                                   tworeal_pairs, &mut arg_used, &mut j);
            if retval != KIN_SUCCESS {
                KINProcessError(Some(kin_mem), retval, line!(), "kinSetFromCommandLine", file!(),
                                &format!("error setting key: {}", tworeal_pairs[j].key));
                return retval;
            }
            if arg_used {
                break 'this_arg;
            }

            /* warn for uninterpreted kinid.X arguments */
            KINProcessError(Some(kin_mem), KIN_WARNING, line!(), "kinSetFromCommandLine", file!(),
                            &format!("WARNING: key {} was not handled\n", argv[idx]));
        }
        idx += 1;
    }

    KIN_SUCCESS
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

    /* options with the default "kinsol." prefix reach the KINSet*
       routines; foreign prefixes are skipped (argv[0] is the program
       name, exactly as in the C argc/argv loop) */
    #[test]
    fn kinsetoptions_parses_int_long_real_and_tworeal_keys() {
        let mut kin_mem = KINMem::default();
        let args = argv(&[
            "prog",
            "kinsol.num_max_iters", "77",
            "other.num_max_iters", "5",
            "kinsol.no_min_eps", "1",
            "kinsol.eta_form", "3",
            "kinsol.damping", "0.5",
            "kinsol.func_norm_tol", "1.5e-5",
            "kinsol.eta_params", "0.7", "1.5",
        ]);
        assert_eq!(KINSetOptions(&mut kin_mem, "", "", &args), KIN_SUCCESS);
        assert_eq!(kin_mem.kin_mxiter, 77); /* "other." prefix skipped */
        assert!(kin_mem.kin_noMinEps);
        assert_eq!(kin_mem.kin_etaflag, KIN_ETACONSTANT);
        assert!(kin_mem.kin_damping);
        assert_eq!(kin_mem.kin_beta, 0.5);
        assert_eq!(kin_mem.kin_fnormtol, 1.5e-5);
        assert_eq!(kin_mem.kin_eta_gamma, 0.7);
        assert_eq!(kin_mem.kin_eta_alpha, 1.5);
    }

    /* a custom kinid replaces the default prefix */
    #[test]
    fn kinsetoptions_honors_custom_kinid() {
        let mut kin_mem = KINMem::default();
        let args = argv(&[
            "prog",
            "mysolver.max_setup_calls", "25",
            "kinsol.max_setup_calls", "40", /* wrong prefix now: skipped */
        ]);
        assert_eq!(KINSetOptions(&mut kin_mem, "mysolver", "", &args), KIN_SUCCESS);
        assert_eq!(kin_mem.kin_msbset, 25);
    }

    /* a failing set routine propagates its error code */
    #[test]
    fn kinsetoptions_propagates_set_errors() {
        let mut kin_mem = KINMem::default();
        let args = argv(&["prog", "kinsol.eta_form", "7"]); /* invalid etachoice */
        assert_eq!(KINSetOptions(&mut kin_mem, "", "", &args), KIN_ILL_INPUT);
    }

    /* file-based options are rejected (KIN_ILL_INPUT) */
    #[test]
    fn kinsetoptions_rejects_file_input() {
        let mut kin_mem = KINMem::default();
        assert_eq!(KINSetOptions(&mut kin_mem, "", "opts.txt", &[]), KIN_ILL_INPUT);
    }

    /* unhandled "kinsol.*" keys only warn; later options still apply */
    #[test]
    fn kinsetoptions_warns_and_continues_on_unknown_key() {
        let mut kin_mem = KINMem::default();
        let args = argv(&["prog", "kinsol.bogus_key", "kinsol.max_beta_fails", "12"]);
        assert_eq!(KINSetOptions(&mut kin_mem, "", "", &args), KIN_SUCCESS);
        assert_eq!(kin_mem.kin_mxnbcf, 12);
    }
}
