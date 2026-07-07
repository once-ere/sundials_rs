# PROGRESS — SUNDIALS 7.7.0 → pure Rust

Status legend: `done(donor)` = inherited verified from cvode_rs donor;
`todo` → `ported` → `committed`. Excluded files are listed once with reason.
Update this file as each unit lands. Resume point after compaction:
CLAUDE.md + this file + `git log`.

## Phase 0 — bootstrap: COMPLETE (see git tag phase0-cvode-green)

## crates/sundials_core

### inherited from donor — committed
nvector_serial.c; sundials_band.c; sundials_context.c; sundials_dense.c;
sundials_errors.c; sundials_iterative.c; sundials_linearsolver.c;
sundials_math.c; sundials_matrix.c; sundials_nonlinearsolver.c;
sundials_nvector.c (fused into nvector_serial per ARCHITECTURE);
sundials_types.h; sundials_utils (priv);
sunmatrix_{dense,band,sparse}.c; sunlinsol_{dense,band,spgmr,spfgmr,spbcgs,sptfqmr,pcg}.c;
sunnonlinsol_{newton,fixedpoint}.c

### Phase 1 — todo
- [x] sundials/sundials_version.c — committed
- [x] sundials/sundials_hashmap.c — committed
- [x] sundials/sundials_logger.c — committed
- [x] sundials/sundials_profiler.c — committed (std::time::Instant timers)
- [x] sundials/sundials_memory.c — committed
- [x] sunmemory/system/sundials_system_memory.c — committed
- [x] sundials/sundials_direct.c — committed
- [x] sundials/sundials_futils.c — committed (std::fs)
- [x] sundials/sundials_cli.c — committed
- [x] sundials/sundials_adaptcontroller.c — committed
- [x] sunadaptcontroller/soderlind/sunadaptcontroller_soderlind.c — committed
- [x] sunadaptcontroller/imexgus/sunadaptcontroller_imexgus.c — committed
- [x] sunadaptcontroller/mrihtol/sunadaptcontroller_mrihtol.c — committed
- [x] sundials/sundials_domeigestimator.c — committed
- [x] sundomeigest/power/sundomeigest_power.c — committed
- [x] sundomeigest/arnoldi/sundomeigest_arnoldi.c — committed (dgeev → EISPACK HQR, documented)
- [x] sundials/sundials_nvector_senswrapper.c — committed
- [ ] sundials/sundials_stepper.c — DEFERRED to Phase 6 start (only arkode adjoint needs it)
- [ ] sundials/sundials_adjointstepper.c — DEFERRED to Phase 6 start (only arkode adjoint needs it)
- [ ] sundials/sundials_adjointcheckpointscheme.c — DEFERRED to Phase 6 start (only arkode adjoint needs it)
- [ ] sunadjointcheckpointscheme/fixed/sunadjointcheckpointscheme_fixed.c — DEFERRED to Phase 6 start
- [ ] sundials/sundials_datanode.c — DEFERRED to Phase 6 start (only arkode adjoint needs it)
- [ ] sundials/sundatanode/sundatanode_inmem.c — DEFERRED to Phase 6 start (stl/sunstl_vector.h → Vec)

### excluded (core)
sundials_mpi_errors.c (MPI); sundials_xbraid.c (XBraid);
sundials_cuda.h/hip.h/sycl.h etc. (GPU headers); fmod_* (Fortran)

## crates/cvode_rs — library committed (donor); examples see VERIFICATION.md
- [x] donor example stubs translated: cvDiurnal_kry_bp, cvKrylovDemo_ls,
      cvParticle_dns, cvPendulum_dns — committed (all LOCAL-C verified)
- [x] cvKrylovDemo_prec — committed (IDENTICAL to shipped .out; exposed+fixed cv_jcur aliasing bug in cvode_ls.rs)
### excluded (cvode)
cvode_fused_gpu.cpp (GPU; stubs ported); examples cvRoberts_klu,
cvRoberts_block_klu, cvRoberts_sps (KLU/SuperLU)

## crates/cvodes_rs — Phase 2
- [x] cvodes_impl.h — committed
- [~] cvodes.c — PARTS 1+2/3 written (lines 1-6244: init/tolerances + CVode driver, Get* extraction, Free family, alloc/free vectors, cvInitialSetup, cvHin family; part 2 type-checked in an isolated harness against pinned PART 3 stubs); part 3 (6245-end, incl. cvStep 5874-6244) pending
- [ ] cvodes_io.c — todo
- [x] cvodes_ls_impl.h — committed
- [x] cvodes_ls.c — committed
- [ ] cvodes_nls.c — todo
- [ ] cvodes_nls_sim.c — todo
- [ ] cvodes_nls_stg.c — todo
- [ ] cvodes_nls_stg1.c — todo
- [x] cvodes_diag_impl.h — committed
- [x] cvodes_diag.c — committed
- [x] cvodes_proj_impl.h — committed
- [x] cvodes_proj.c — committed
- [ ] cvodes_bandpre_impl.h — todo
- [ ] cvodes_bandpre.c — todo
- [ ] cvodes_bbdpre_impl.h — todo
- [ ] cvodes_bbdpre.c — todo
- [ ] cvodes_resize.c — todo
- [ ] cvodes_cli.c — todo
- [ ] cvodea.c — todo
- [ ] cvodea_io.c — todo
### excluded (cvodes) examples
cvsRoberts_klu, cvsRoberts_sps, cvsRoberts_ASAi_klu, cvsRoberts_ASAi_sps,
cvsRoberts_FSA_klu, cvsRoberts_FSA_sps (KLU/SuperLU)

## crates/kinsol_rs — Phase 3
- [x] kinsol_impl.h — ported (kinsol_impl.rs: KINMem, kinsol.h constants/typedefs, KINProcessError; uncommitted)
- [x] kinsol.c — committed
- [x] kinsol_io.c — committed
- [x] kinsol_ls_impl.h — ported (kinsol_ls_impl.rs: KINLsMem, kinsol_ls.h types/codes; uncommitted)
- [x] kinsol_ls.c — committed
- [x] kinsol_aa.c — committed
- [x] kinsol_orth.c — committed
- [x] kinsol_bbdpre_impl.h — ported (kinsol_bbdpre_impl.rs: KBBDPrecData, KINBBDLocalFn/KINBBDCommFn; uncommitted)
- [x] kinsol_bbdpre.c — ported (kinsol_bbdpre.rs; PrecModule::BBDPre variant + kinLsPSetup/psolve/setup_disabled dispatch in kinsol_ls(_impl).rs; uncommitted)
- [x] kinsol_cli.c — ported (kinsol_cli.rs over sundials_core::sundials_cli generics; uncommitted)
### excluded (kinsol) examples
kinFerTron_klu (KLU), kinRoboKin_slu (SuperLU)

## crates/ida_rs — Phase 4: LIBRARY COMPLETE (every ida/*.c ported)
- [x] ida_impl.h — committed (IDAMem)
- [x] ida.c — committed (PART 1 init/tolerances/RootInit/IDASolve/GetDky/alloc/
      InitialSetup; PART 2 IDAStep/Predict/TestError/HandleNFlag/CompleteStep/
      Restore/Reset + rootfinding Rcheck1-3/Rootfind + GetSolution)
- [x] ida_ic.c — committed (IDACalcIC + IDANlsIC/NewtonIC/LineSrch/fnorm/
      Newyyp/Newy/ICFailFlag; workspace aliasing delnew=phi[2]/dtemp=phi[3]/
      ynew=tempv2/ypnew=ee)
- [x] ida_io.c — committed (62/62 fns: all Set/Get families, IC setters,
      IDAPrintAllStats byte-format helpers, IDAGetReturnFlagName)
- [x] ida_ls_impl.h — committed (IDALsMem; PrecModule +BBDPre variant)
- [x] ida_ls.c — committed (full LS interface; BBDPre psetup/psolve dispatch)
- [x] ida_nls.c — committed (IDANls Newton solve collapsed onto IDAMem)
- [x] ida_bbdpre_impl.h — committed (IBBDPrecData)
- [x] ida_bbdpre.c — committed (serial band reduction; DQ Jacobian + band LU
      round-trip test green)
- [x] ida_cli.c — committed (IDASetOptions over sundials_core::sundials_cli)
- 17/17 ida_rs unit tests green; crate builds warning-free.
### ida examples — COMPLETE (Phase 4 verification): all 8 serial examples ported + verified
  idaRoberts_dns, idaAnalytic_mels, idaHeat2D_bnd, idaHeat2D_kry,
  idaKrylovDemo_ls, idaSlCrank_dns = IDENTICAL to shipped; idaFoodWeb_kry =
  IDENTICAL; idaFoodWeb_bnd = local-C (foreign-libm last-digit h drift). See
  VERIFICATION.md. Library fix landed: idaLsSolve returned x(=0) instead of the
  preconditioned residual on 0-iteration Krylov solves (LinearSolver::resid()
  now exposes it) — was making every IDA Krylov example diverge.
  NEXT: Phase 5 (idas_rs). cvode_ls.rs has the same latent nli==0 bug (task
  chip spawned) — fix before cvsKrylov* examples.
### excluded (ida) examples
idaHeat2D_klu, idaRoberts_klu, idaRoberts_sps (KLU/SuperLU)

## crates/idas_rs — Phase 5
- [ ] idas_impl.h — todo
- [ ] idas.c — todo
- [ ] idas_ic.c — todo
- [ ] idas_io.c — todo
- [ ] idas_ls_impl.h — todo
- [ ] idas_ls.c — todo
- [ ] idas_nls.c — todo
- [ ] idas_nls_sim.c — todo
- [ ] idas_nls_stg.c — todo
- [ ] idas_bbdpre_impl.h — todo
- [ ] idas_bbdpre.c — todo
- [ ] idas_cli.c — todo
- [ ] idaa.c — todo
- [ ] idaa_io.c — todo
### excluded (idas) examples
idasRoberts_klu, idasRoberts_sps, idasRoberts_ASAi_klu, idasRoberts_ASAi_sps,
idasRoberts_FSA_klu, idasRoberts_FSA_sps (KLU/SuperLU)

## crates/arkode_rs — Phase 6
- [ ] arkode_types_impl.h — todo
- [ ] arkode_impl.h — todo
- [ ] arkode_butcher.c — todo
- [ ] arkode_butcher_erk.c (+ .def) — todo
- [ ] arkode_butcher_dirk.c (+ .def) — todo
- [ ] arkode_interp_impl.h — todo
- [ ] arkode_interp.c — todo
- [ ] arkode.c — todo
- [ ] arkode_io.c — todo
- [ ] arkode_adapt_impl.h — todo
- [ ] arkode_adapt.c — todo
- [ ] arkode_root_impl.h — todo
- [ ] arkode_root.c — todo
- [ ] arkode_erkstep_impl.h — todo
- [ ] arkode_erkstep.c — todo
- [ ] arkode_erkstep_io.c — todo
- [ ] arkode_ls_impl.h — todo
- [ ] arkode_ls.c — todo
- [ ] arkode_arkstep_impl.h — todo
- [ ] arkode_arkstep.c — todo
- [ ] arkode_arkstep_io.c — todo
- [ ] arkode_arkstep_nls.c — todo
- [ ] arkode_sprk.c — todo
- [ ] arkode_sprkstep_impl.h — todo
- [ ] arkode_sprkstep.c — todo
- [ ] arkode_sprkstep_io.c — todo
- [ ] arkode_lsrkstep_impl.h — todo
- [ ] arkode_lsrkstep.c — todo
- [ ] arkode_lsrkstep_io.c — todo
- [ ] arkode_mri_tables.c (+ .def) — todo
- [ ] arkode_mristep_impl.h — todo
- [ ] arkode_mristep.c — todo
- [ ] arkode_mristep_io.c — todo
- [ ] arkode_mristep_nls.c — todo
- [ ] arkode_mristep_controller.c — todo
- [ ] arkode_splittingstep_coefficients.c (+ .def) — todo
- [ ] arkode_splittingstep_impl.h — todo
- [ ] arkode_splittingstep.c — todo
- [ ] arkode_forcingstep_impl.h — todo
- [ ] arkode_forcingstep.c — todo
- [ ] arkode_bandpre_impl.h — todo
- [ ] arkode_bandpre.c — todo
- [ ] arkode_bbdpre_impl.h — todo
- [ ] arkode_bbdpre.c — todo
- [ ] arkode_relaxation_impl.h — todo
- [ ] arkode_relaxation.c — todo
- [ ] arkode_user_controller.c/.h — todo
- [ ] arkode_sunstepper.c — todo
- [ ] arkode_cli.c — todo
### excluded (arkode)
xbraid/ (XBraid); fmod_* (Fortran)

## Phase 7 — docs
- [ ] sundials.md — todo
- [ ] final CLAUDE.md / ARCHITECTURE.md pass — todo
