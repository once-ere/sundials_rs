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

2026-07-11: that scratchpad did not survive; the C library was rebuilt with
the same configuration and the local outputs now live in the COMMITTED
`localref/` tree at the workspace root (currently idas/serial only —
regenerate other solvers there if a new local-C comparison is needed).
`tools/verify_examples.sh` defaults to it (`SUNDIALS_LOCALREF` overrides).

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
| idaFoodWeb_kry | identical (SPGMR + user block-diag precond + IDACalcIC; fixed by nli==0 residual copy) |
| idaHeat2D_bnd | identical (band DQ Jacobian + IDACalcIC + constraints) |
| idaHeat2D_kry | identical (SPGMR + diagonal precond, 2 GS-type cases) |
| idaKrylovDemo_ls | identical (loops SPGMR/SPBCGS/SPTFQMR + diagonal precond) |
| idaSlCrank_dns | identical (index-2 GGL DAE, dense DQ Jacobian, IDASetSuppressAlg) |
| idaHeat2D_klu / idaRoberts_klu / idaRoberts_sps | excluded(KLU/SuperLU) |

### idaFoodWeb_kry / idaHeat2D_kry — RESOLVED (was: nst/nli divergence)

Root cause: idaLsSolve's iterative finish path returned x (=0) as the Newton
correction whenever GMRES converged in 0 iterations (the preconditioner alone
met the linear tolerance). C copies SUNLinSolResid(LS) — the left-
preconditioned residual of the zero guess — instead. Returning 0 gave a
spurious zero correction that derailed the step trajectory. Fixed by exposing
each iterative solver's residual vector (LinearSolver::resid()) and copying it
into b when nli_inc==0. NOTE: cvode_ls.rs has the same latent shortcut but no
cvode example currently hits nli==0; worth fixing for parity/robustness.

## idas_rs (Phase 5)

| example | status |
|---|---|
| idasAkzoNob_ASAi_dns | local-C(byte-identical to local C build; shipped .out stale trailing-space/blank-line vs current C printf) |
| idasAkzoNob_dns | identical |
| idasAnalytic_mels | identical (incl. idas.init_step_1e-5 CLI variant) |
| idasFoodWeb_bnd | local-C(byte-identical stdout to local C 7.7.0 -ffp-contract=off build; shipped .out foreign-libm, last-digit h drift at t≥0.7 — same signature as idaFoodWeb_bnd) |
| idasHeat2D_bnd | identical |
| idasHeat2D_kry | identical |
| idasHessian_ASA_FSA | identical |
| idasKrylovDemo_ls | identical (incl. _1/_2 nrmfactor arg variants) |
| idasRoberts_ASAi_dns | identical |
| idasRoberts_FSA_dns | identical (shipped stg_t ref; sim_t/sim_f/stg_f/-nosensi byte-identical to local C build, localref committed) |
| idasRoberts_dns | identical |
| idasSlCrank_FSA_dns | tolerance-level(stats/G-tail/one dG/dp digit vs local C; FSA dG/dp agrees with the example's own FD checks, FD sections byte-identical to C. The C binary is flag-unstable here: 232 vs 263 steps between -fmath-errno/default builds of the same source (sincos fusion in the DQ-perturbed residual). Exposed+fixed the IDASetSensParams p-copy defect via the pinned FSAUserData convention — see ARCHITECTURE.md §3.6) |
| idasSlCrank_dns | last-digit(only the G quadrature line, 13th digit vs local C; all 44 trajectory/stats lines byte-identical. Root-caused by bit-level instrumentation: library-side values are bit-identical through the first divergent step — the deviation enters inside the EXAMPLE's own residual, where apple-clang fuses adjacent sin/cos into __sincos_stret (rare 1-ulp differences vs separate calls, first hit at step 74). The C binary itself is flag-fragile at 1e-11 here: -fmath-errno moves G to ...3317697 vs the fused build's ...3378475; ours ...3381925 is 3.5e-13 from the fused build) |
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

Module-only units without any serial C driver are validated by unit test
instead of a reference diff: ARKBANDPRE and ARKBBDPRE (C drivers are
MPI-only) are covered by `arkode_bandpre::tests` / `arkode_bbdpre::tests`
— implicit 1D heat equation with ARKStep + SPGMR(left) + the module
preconditioner vs the exact semi-discrete solution, the internally
generated DQ Jacobian vs the analytic tridiagonal Jacobian, band-vs-BBD
trajectory cross-agreement (<1e-8), and counter/workspace/ReInit/
error-path checks.

| example | status |
|---|---|
| ark_KrylovDemo_prec | IDENTICAL incl. _1/_2 nrmfactor arg variants (shipped refs; matrix-free SPGMR + user Precond/PSolve + DQ Jtimes + ARKStepReInit + workspace accounting, all four jpre/gstype runs byte-exact) |
| ark_advection_diffusion_reaction_splitting | todo |
| ark_analytic | IDENTICAL x2 (shipped refs; first implicit ARKStep example — DIRK + Newton + dense ARKLS + SetLinear; the scalar_tolerances/table_names CLI arg variant passes via ARKodeSetOptions + arkStep_SetOptions with an ESDIRK547L2SA table override) |
| ark_analytic_lsrk (+_varjac, _domeigest, _domeigest arg variant, _ssprk) | LOCAL-C x5 (byte-identical to fresh 7.7.0 C build incl. all stats; shipped refs stale; covers RKL/RKC STS with user dom-eig fn, power-iteration SUNDomEigEstimator with warmup/max_iters options, and the SSP(9,3) path) |
| ark_analytic_mels | IDENTICAL (shipped ref; custom matrix-embedded SUNLinearSolver via the core CustomLinSol trait — solve callback receives (t, gamma, user_data) per the pinned cvAnalytic_mels adaptation; exposed a latent bug where linit ran with step_mem detached so the MELS lsetup-disable no-op'd) |
| ark_analytic_nonlin | IDENTICAL (shipped ref; the FIRST arkode example, verified byte-exact) |
| ark_analytic_partitioned | LOCAL-C x5 (splitting default, forcing, and BEST_2_2_2/RUTH_3_3_2/YOSHIDA_8_6_2 named-coefficient variants all byte-identical to a fresh 7.7.0 C build; shipped refs predate the SUN_TABLE_WIDTH 28->29 change; exercises SplittingStep + ForcingStep over ERKStep/ARKStep-backed SUNSteppers incl. ARKStep inner forcing) |
| ark_brusselator | IDENTICAL (shipped ref; exercises deduce_rhs + SetAutonomous/TrivialPredAutonomous residual + Lagrange interpolant + Newton failure-retry paths, verified byte-exact) |
| ark_brusselator1D | todo |
| ark_brusselator1D_imexmri | todo |
| ark_brusselator_1D_mri | IDENTICAL (shipped ref; explicit-slow MRIStep over an implicit adaptive ARKStep inner (ARK324L2SA DIRK) with band ARKLS + user Jacobian, verified byte-exact on first run) |
| ark_brusselator_fp | IDENTICAL (shipped ref; ImEx ARK pair + Anderson-accelerated fixed-point NLS + autonomous TPA reuse, verified byte-exact; the _fp_1 monitor-arg reference is byte-identical to the no-arg one and also passes) |
| ark_brusselator_mri | IDENTICAL (shipped ref; the first MRIStep example — MIS_KW3 slow coupling over an explicit KNOTH_WOLKE_3_3 ARKStep inner via ARKodeCreateMRIStepInnerStepper, verified byte-exact on first run) |
| ark_conserved_exp_entropy_ark | LOCAL-C x2 (1_0 explicit ERK relax + 1_1 DIRK relax with dense ARKLS+Newton, byte-identical to a fresh 7.7.0 C build; shipped refs stale at 1e-14 delta-e/Newton-iteration level; exercises arkStep_RelaxDeltaE stored-stage (implicit) and reconstructed-stage (explicit) paths + the arkRelax Newton solve) |
| ark_conserved_exp_entropy_erk | LOCAL-C (1-arg relax-on run byte-identical to a fresh 7.7.0 C build; shipped ref stale; exercises erkStep_RelaxDeltaE + the FSAL-disabled-under-relaxation FullRHS path) |
| ark_damped_harmonic_symplectic | todo |
| ark_dissipated_exp_entropy | LOCAL-C x2 (1_0 explicit + 1_1 DIRK relaxed runs byte-identical to a fresh 7.7.0 C build; shipped refs stale; scalar entropy-dissipation counterpart of the conserved problem) |
| ark_harmonic_symplectic | todo |
| ark_heat1D | todo |
| ark_heat1D_adapt | todo |
| ark_kepler | 1 IDENTICAL + 12 LOCAL-C (all 13 arg variants byte-match a freshly built 7.7.0 C binary; shipped .out files predate the SUN_TABLE_WIDTH 28->29 change; covers SPRK standard + compensated-sum steps, 9 SPRK tables, SPRK rootfinding, explicit ARKStep ERK fixed/adaptive, tstop, check-order) |
| ark_kpr_mri | todo |
| ark_lotka_volterra_ASA | todo |
| ark_onewaycouple_mri | todo |
| ark_reaction_diffusion_mri | todo |
| ark_robertson | IDENTICAL (shipped ref; stiff DIRK + arkPredict_MaximumOrder (predictor 1) + full ARKodePrintAllStats TABLE output to stdout, verified byte-exact) |
| ark_robertson_constraints | IDENTICAL (shipped ref; exercises ARKodeSetConstraints + the arkCheckConstraints path and ARKodeGetNumConstrFails, verified byte-exact on first run) |
| ark_robertson_root | IDENTICAL (shipped ref; exercises ARKodeRootInit/GetRootInfo end-to-end — Illinois rootfinder, SVtolerances, arkHin auto initial step — verified byte-exact) |
| ark_twowaycouple_mri | todo |
