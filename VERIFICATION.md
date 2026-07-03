# VERIFICATION — example outputs vs upstream references

Run `tools/verify_examples.sh [crate|all]`, then Read `logs/summary.txt`.
Reference files: `../sundials-7.7.0/examples/<solver>/<serial-dir>/<name>.out`.
Statuses: `identical` | `last-digit(reason)` | `local-C(reason)` (matches a
locally-built C binary byte-for-byte; shipped .out from a foreign libm) |
`noref` | `todo` | `excluded(reason)`.

## cvode_rs (donor; re-verified in workspace 2026-07-03)

| example | status |
|---|---|
| cvRoberts_dns | identical |
| cvRoberts_dnsL | last-digit(LAPACK ref vs native LU; donor-documented) |
| cvRoberts_dns_uw | identical |
| cvRoberts_dns_constraints | identical |
| cvRoberts_dns_negsol | local-C(shipped .out has stale stats-line spacing; port matches current C source) |
| cvAdvDiff_bnd | identical |
| cvAdvDiff_bndL | identical |
| cvDiurnal_kry | local-C(shipped .out from foreign libm; donor-documented) |
| cvDiurnal_kry_bp | todo (donor stub — translate) |
| cvDirectDemo_ls | local-C(foreign-libm ref; donor-documented) |
| cvKrylovDemo_ls | todo (donor stub — translate) |
| cvKrylovDemo_prec | todo (donor stub — translate) |
| cvDisc_dns | identical |
| cvAnalytic_mels | identical |
| cvParticle_dns | todo (donor stub — translate) |
| cvPendulum_dns | todo (donor stub — translate) |
| cvRocket_dns | identical |
| solar_system | noref (new example; energy-conservation self-check) |
| cvRoberts_klu / cvRoberts_block_klu / cvRoberts_sps | excluded(KLU/SuperLU) |

## cvodes_rs (Phase 2)

| example | status |
|---|---|
| cvsRoberts_dns | todo |
| cvsRoberts_dnsL | todo |
| cvsRoberts_dns_uw | todo |
| cvsRoberts_dns_constraints | todo |
| cvsAdvDiff_bnd | todo |
| cvsAdvDiff_bndL | todo |
| cvsDiurnal_kry | todo |
| cvsDiurnal_kry_bp | todo |
| cvsDirectDemo_ls | todo |
| cvsKrylovDemo_ls | todo |
| cvsKrylovDemo_prec | todo |
| cvsAnalytic_mels | todo |
| cvsParticle_dns | todo |
| cvsPendulum_dns | todo |
| cvsAdvDiff_FSA_non | todo |
| cvsDiurnal_FSA_kry | todo |
| cvsRoberts_FSA_dns | todo |
| cvsRoberts_FSA_dns_Switch | todo |
| cvsRoberts_FSA_dns_constraints | todo |
| cvsAdvDiff_ASAi_bnd | todo |
| cvsFoodWeb_ASAi_kry | todo |
| cvsFoodWeb_ASAp_kry | todo |
| cvsHessian_ASA_FSA | todo |
| cvsLotkaVolterra_ASA | todo |
| cvsRoberts_ASAi_dns | todo |
| cvsRoberts_ASAi_dns_constraints | todo |
| cvsRoberts_klu / _sps / _ASAi_klu / _ASAi_sps / _FSA_klu / _FSA_sps | excluded(KLU/SuperLU) |

## kinsol_rs (Phase 3)

| example | status |
|---|---|
| kinAnalytic_fp | todo |
| kinFerTron_dns | todo |
| kinFoodWeb_kry | todo |
| kinKrylovDemo_ls | todo |
| kinLaplace_bnd | todo |
| kinLaplace_picard_bnd | todo |
| kinLaplace_picard_kry | todo |
| kinRoberts_fp | todo |
| kinRoboKin_dns | todo |
| kinFerTron_klu / kinRoboKin_slu | excluded(KLU/SuperLU) |

## ida_rs (Phase 4)

| example | status |
|---|---|
| idaAnalytic_mels | todo |
| idaFoodWeb_bnd | todo |
| idaFoodWeb_kry | todo |
| idaHeat2D_bnd | todo |
| idaHeat2D_kry | todo |
| idaKrylovDemo_ls | todo |
| idaRoberts_dns | todo |
| idaSlCrank_dns | todo |
| idaHeat2D_klu / idaRoberts_klu / idaRoberts_sps | excluded(KLU/SuperLU) |

## idas_rs (Phase 5)

| example | status |
|---|---|
| idasAkzoNob_ASAi_dns | todo |
| idasAkzoNob_dns | todo |
| idasAnalytic_mels | todo |
| idasFoodWeb_bnd | todo |
| idasHeat2D_bnd | todo |
| idasHeat2D_kry | todo |
| idasHessian_ASA_FSA | todo |
| idasKrylovDemo_ls | todo |
| idasRoberts_ASAi_dns | todo |
| idasRoberts_FSA_dns | todo |
| idasRoberts_dns | todo |
| idasSlCrank_FSA_dns | todo |
| idasSlCrank_dns | todo |
| idasRoberts_klu / _sps / _ASAi_klu / _ASAi_sps / _FSA_klu / _FSA_sps | excluded(KLU/SuperLU) |

## arkode_rs (Phase 6)

| example | status |
|---|---|
| ark_KrylovDemo_prec | todo |
| ark_advection_diffusion_reaction_splitting | todo |
| ark_analytic | todo |
| ark_analytic_lsrk | todo |
| ark_analytic_lsrk_domeigest | todo |
| ark_analytic_lsrk_varjac | todo |
| ark_analytic_mels | todo |
| ark_analytic_nonlin | todo |
| ark_analytic_partitioned | todo |
| ark_analytic_ssprk | todo |
| ark_brusselator | todo |
| ark_brusselator1D | todo |
| ark_brusselator1D_imexmri | todo |
| ark_brusselator_1D_mri | todo |
| ark_brusselator_fp | todo |
| ark_brusselator_mri | todo |
| ark_conserved_exp_entropy_ark | todo |
| ark_conserved_exp_entropy_erk | todo |
| ark_damped_harmonic_symplectic | todo |
| ark_dissipated_exp_entropy | todo |
| ark_harmonic_symplectic | todo |
| ark_heat1D | todo |
| ark_heat1D_adapt | todo |
| ark_kepler | todo |
| ark_kpr_mri | todo |
| ark_lotka_volterra_ASA | todo |
| ark_onewaycouple_mri | todo |
| ark_reaction_diffusion_mri | todo |
| ark_robertson | todo |
| ark_robertson_constraints | todo |
| ark_robertson_root | todo |
| ark_twowaycouple_mri | todo |
