# VERIFICATION — example outputs vs upstream references

Run `tools/verify_examples.sh [crate|all]`, then Read `logs/summary.txt`.
Reference files: `../sundials-7.7.0/examples/<solver>/<serial-dir>/<name>.out`.
Statuses: `identical` | `last-digit(reason)` | `local-C(reason)` (matches a
locally-built C binary byte-for-byte; shipped .out from a foreign libm) |
`noref` | `todo` | `excluded(reason)`.

2026-07-03: C SUNDIALS 7.7.0 was built locally (Release, `-ffp-contract=off`,
no LAPACK/KLU) in the session scratchpad; every serial example's local output
is saved under `<scratchpad>/localref/`. All `local-C` statuses below were
re-confirmed byte-for-byte against those binaries on this machine.

## cvode_rs (donor; re-verified in workspace 2026-07-03)

| example | status |
|---|---|
| cvRoberts_dns | identical |
| cvRoberts_dnsL | last-digit(LAPACK ref vs native LU; donor-documented) |
| cvRoberts_dns_uw | identical |
| cvRoberts_dns_constraints | identical |
| cvRoberts_dns_negsol | local-C(byte-identical to local C build; shipped .out stale stats-line spacing) |
| cvAdvDiff_bnd | identical |
| cvAdvDiff_bndL | identical |
| cvDiurnal_kry | local-C(byte-identical to local C build) |
| cvDiurnal_kry_bp | local-C(byte-identical to local C build; shipped .out foreign-libm, drift onset t=2.88e4) |
| cvDirectDemo_ls | local-C(byte-identical to local C build) |
| cvKrylovDemo_ls | local-C(byte-identical to local C build; shipped .out foreign-libm) |
| cvKrylovDemo_prec | identical |
| cvDisc_dns | identical |
| cvAnalytic_mels | identical |
| cvParticle_dns | local-C(byte-identical to local C -ffp-contract=off build; shipped .out foreign platform, 100-orbit chaotic amplification) |
| cvPendulum_dns | local-C(shipped .out has stale one-digit atol exponent unproducible by current C source; all numeric rows byte-identical) |
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
| kinAnalytic_fp | identical (incl. all 10 shipped arg-variant .out files) |
| kinFerTron_dns | identical |
| kinFoodWeb_kry | identical |
| kinKrylovDemo_ls | identical (4 solver passes) |
| kinLaplace_bnd | identical |
| kinLaplace_picard_bnd | identical |
| kinLaplace_picard_kry | local-C(shipped .out adds trailing newline absent from C printf; content identical) |
| kinRoberts_fp | identical (incl. kinsol.m_aa_1 KINSetOptions variant) |
| kinRoboKin_dns | local-C(byte-identical stdout+csv to local C build; shipped .out stale SUN_TABLE_WIDTH 28 vs 29 in 7.7.0, shipped .csv stale %g-style reals vs SUN_FORMAT_E) |
| kinFerTron_klu / kinRoboKin_slu | excluded(KLU/SuperLU) |

## ida_rs (Phase 4)

| example | status |
|---|---|
| idaRoberts_dns | identical (stdout; shipped .csv is stale pre-7.7.0 stats format, our CSV matches 7.7.0 ida_io.c) |
| idaAnalytic_mels | identical (custom matrix-embedded LS) |
| idaFoodWeb_bnd | local-C(byte-identical stdout to local C 7.7.0 -ffp-contract=off build; shipped .out foreign-libm, last-digit h drift at t≥0.7 after sin() usage) |
| idaFoodWeb_kry | solution-correct; UNRESOLVED nst/nli divergence (see note) |
| idaHeat2D_bnd | identical (band DQ Jacobian + IDACalcIC + constraints) |
| idaHeat2D_kry | todo |
| idaKrylovDemo_ls | todo |
| idaSlCrank_dns | todo |
| idaHeat2D_klu / idaRoberts_klu / idaRoberts_sps | excluded(KLU/SuperLU) |

### idaFoodWeb_kry note (open issue)

Concentration columns (t, bottom-left, top-right) match the shipped .out
(== local C 7.7.0) to displayed precision, but nst/k/h and the final-stats
counters (nli 641 vs 1034, nps, nst 188 vs 165) diverge. Root cause NOT
found despite deep instrumented C-vs-Rust comparison. Established, on this
machine (so not libm — idaFoodWeb_bnd is RUST==LOCALC):
- First Precond call byte-identical (hh, ewt, cxy, cpxy, rates, inc, P block).
- Preconditioner `rates` always consistent with cc; delta/tol identical every
  solve; solve boundaries identical; SPGMR default gstype = MODIFIED_GS (same).
- Jv formula (idaLsDQJtimes) + linear_sum_with a==-b VScaleDiff special case
  verified correct; nrmfac/eplifac/dqincfac identical; restart path not taken.
- Shared SPGMR core is byte-identical in cvKrylovDemo_prec (both GS types),
  and IDA/CVODE both pass s1=s2=weight, so the scaled-GMRES path is exercised.
- Divergence localises to the Krylov vector V input to Jv call ~13 (v-norm
  0.577 vs 0.799 — gross, not sub-ULP) after Jv calls 0-12 matched exactly;
  i.e. a GMRES solve converged in a different iteration count. The first ~11
  Jv calls have near-constant v-norm ≈1.00002 (likely the IDACalcIC phase).
Next lead to try: instrument the SPGMR core (rho/givens/Hessenberg per inner
iteration) rebuilt for C, and IDA's IC-phase linear solves specifically; or
element-wise (not checksum) compare Jv call 12 output for a hidden sub-ULP
seed. Scratch has kry_dbg.c/ida_ls_dbg.o harness for instrumented C builds.

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

NOTE for the Phase 6 harness: many arkode references encode command-line
arguments in the filename (e.g. `ark_brusselator1D_imexmri_2_0.001.out` ←
argv "2 0.001"; `ark_kepler_--stepper_ERK_--step-mode_adapt.out`). There are
76 reference outputs for 32 examples. Verification must run each example
once per reference file with the args decoded from its suffix, and local-C
references must be regenerated the same way (C binaries live in the
scratchpad sunbuild). One VERIFICATION line per reference file, not per
example.

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
