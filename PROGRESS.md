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
- [x] sundials/sundials_stepper.c — committed (Phase 6 start; ops stay an Option<fn> table since integrators register them one at a time via Set*Fn; content = UserData; C's uninitialized reinit/resetcheckpointindex/getnumsteps ops noted)
- [x] sundials/sundials_adjointstepper.c — committed (owning Option slots for fwd/adj SUNSteppers + checkpoint scheme; own_* flags honored in Destroy; ReInit discards ReInit retvals as C does)
- [x] sundials/sundials_adjointcheckpointscheme.c — committed (Option<fn> registration table like SUNStepper; content = UserData)
- [x] sunadjointcheckpointscheme/fixed/sunadjointcheckpointscheme_fixed.c — committed (cached node aliases re-derived by key lookup from step_num fields; C UAF corner on re-loading a consumed step returns CHECKPOINT_NOT_FOUND instead)
- [x] sundials/sundials_datanode.c — committed (enum dispatch; tree owns children, Get* lend &mut borrows, Remove* move the child out; Destroy = drop)
- [x] sundials/sundatanode/sundatanode_inmem.c — committed (SUNStlVector → Vec, hashmap children own nodes, write-only parent ptr dropped, leaf = SUNMemory bytes [t|BufPack(v)]; N_VBufSize/Pack/Unpack added to nvector_serial)

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

## crates/cvodes_rs — Phase 2: LIBRARY COMPLETE (every cvodes/*.c ported)
- [x] cvodes_impl.h — committed
- [x] cvodes.c — committed (PART 1 init/tolerances, PART 2 CVode driver +
      extraction + Free + cvInitialSetup/cvHin, PART 3 cvStep machinery, NLS
      drivers, error/eta handling, rootfinding, ewt/norms, sens DQ RHS)
- [x] cvodes_io.c — committed
- [x] cvodes_ls_impl.h — committed
- [x] cvodes_ls.c — committed
- [x] cvodes_nls.c — committed
- [x] cvodes_nls_sim.c — committed
- [x] cvodes_nls_stg.c — committed
- [x] cvodes_nls_stg1.c — committed
- [x] cvodes_diag_impl.h — committed
- [x] cvodes_diag.c — committed
- [x] cvodes_proj_impl.h — committed
- [x] cvodes_proj.c — committed
- [x] cvodes_bandpre_impl.h — committed
- [x] cvodes_bandpre.c — committed
- [x] cvodes_bbdpre_impl.h — committed
- [x] cvodes_bbdpre.c — committed
- [x] cvodes_resize.c — committed
- [x] cvodes_cli.c — committed
- [x] cvodea.c — committed (adjoint: checkpointing, backward problems, Hermite/polynomial interpolation)
- [x] cvodea_io.c — committed
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
- [x] idas_impl.h — ported (idas_impl.rs: IDAMem + quad/sens/quadSens/adjoint
      extensions, idas.h constants + callback types incl. backward-problem
      types, IDAProcessError, all IDAS/IDAA messages, adjoint blocks
      IDAckpntMem/DtpntContent/IDAdtpntMem/IDABMem/IDAadjMem; modeling mirrors
      cvodes_impl.rs pinned decisions — Vecs for linked lists, Option<usize>
      indices, DtpntContent enum, interp dispatch on ia_interpType, senswrapper
      aliases not stored; LsModule has None variant only until idas_ls_impl.h
      lands; verified `cargo build -p idas_rs` clean + zero warnings and
      donor-carried default test on Linux sandbox 2026-07-08; round-trip diff
      vs disk IDENTICAL; uncommitted)
- [x] idas.c — PORTED 2026-07-08 (Parts 1–3 complete, ≈ 7100 lines;
      NOT registered in lib.rs / NOT compile-verified — pending
      idas_nls.rs, idas_nls_sim.rs, idas_nls_stg.rs, idas_ls.rs.
      Port log follows.) (idas.rs Part 1 on disk 2026-07-08, covers
      idas.c ≈ lines 1–3280: ALL init/tolerance exported functions
      (constants incl. CENTERED1/2+FORWARD1/2, ida_msg_g, IDACreate,
      IDAInit, IDAReInit, all 14 tolerance setters, IDAQuad/Sens/QuadSens
      Init+ReInit, IDASensToggleOff, IDARootInit) PLUS the dispatch
      helpers (ida_linit/lperf/efun — donor-verbatim) and the FULL
      IDASolve driver: donor ida.rs yy/yp copy-back aliasing pattern
      carried at every return path; IDAS deltas ported from idas.c —
      entry DQ checks (resSDQ/rhsQSDQ → MSG_NULL_P; user_dataS/QS self-
      pointers not stored per pinned convention), first-call quad/sens
      ypnorm updates + phiQ[1]/phiS[1]/phiQS[1] × hh scaling (fused ops
      inline), step-loop ewtQ/ewtS/ewtQS resets via take/restore weight
      detachment (EWTQ/EWTS/EWTQS_NOW_BAD paths) and quad/sens norm
      updates in the too-much-accuracy test.
      CAUTION: file is NOT registered in lib.rs and NOT compile-verified —
      forward refs: IDACheckNvector / IDAAllocVectors / IDAQuad-Sens-
      QuadSensAllocVectors / IDAInitialSetup / IDAEwtSet + Quad/Sens/
      QuadSensEwtSet / IDAWrmsNorm + norm-update family / IDARcheck1-3 /
      IDAStopTest1-2 / IDAHandleFailure / IDAGetSolution / IDAStep (later
      sections of this file) and crate::idas_nls / idas_nls_sim /
      idas_nls_stg / idas_ls (later units). Pinned signatures for later
      sections: IDAGetSolution(&IDAMem,f64,&mut NVector,&mut NVector);
      IDAWrmsNorm(&IDAMem,&NVector,&NVector,bool)->f64; norm updates take
      &IDAMem + shared refs; Quad/Sens/QuadSensEwtSet(&IDAMem, cur,
      &mut w)->i32 with weights detached at call sites. Do NOT add
      `pub mod idas;` yet. DONE 2026-07-08 (cont.): full extraction family
      — IDAGetDky (donor-verbatim), IDAGetQuad/QuadDky (tfuzz WITHOUT
      SUNRabs/negation preserved as C quirk; k bound is ida_kk),
      IDAGetSens/SensDky/Sens1/SensDky1, IDAGetQuadSens/QuadSensDky/
      QuadSens1/QuadSensDky1 (IDAGetQuadSens returns IDA_NO_SENS with
      MSG_NO_QUADSENSI — C-exact quirk; Dky1 k bound is ida_kused;
      Xvecs gathers dropped, fused kernels inline over phiS/phiQS[j][is]);
      IDAComputeY/Yp (donor-verbatim, signatures verified against donor
      ida.rs: (&IDAMem, ycor: &NVector, out: &mut NVector)) +
      IDAComputeYSens/YpSens ((&IDAMem, &[NVector], &mut [NVector]),
      LinearSumVectorArray → per-vector N_VLinearSum loop); IDAFree
      (donor-verbatim: pub fn IDAFree(_: Box<IDAMem>) {} — RAII);
      IDAQuadFree (mid-lifetime: flag resets + forward-ref
      IDAQuadFreeVectors for lrw/liw bookkeeping).
      NEXT: read idas.c tail (≥ line 4200) for IDASensFree tail +
      IDAQuadSensFree, then port them (NLSsim/NLSstg = None +
      ownNLSsim/stg + sim/stgMallocDone flag resets; senswrapper
      aliases not stored per pinned convention) and the private
      alloc/dealloc section (IDACheckNvector + IDAAllocVectors donor
      text in hand incl. ida_yy/ida_yp owned stand-ins not counted in
      lrw/liw; Quad/Sens/QuadSensAllocVectors + FreeVectors with
      lrw/liw accounting), then IDAInitialSetup (donor base text in
      hand; IDAS adds quad/sens tolerance checks + ewtQ/S/QS init —
      C text pending same read), ewt-set family, stopping tests,
      norms. DONE 2026-07-08 (cont. 2): idas.c fully banked in
      transcript (tail=4720 read, ≈ 4184–EOF); IDASensFree tail
      confirmed = atolSmin0 only (atolQSmin0 lives in IDAQuadSensFree)
      — both ported (NLSsim/NLSstg = None, ownNLS*/sim/stgMallocDone
      flag resets, wrappers-not-stored comment); full private
      alloc/dealloc section: IDACheckNvector, IDAAllocVectors
      (donor + yy/yp stand-ins), IDAFreeVectors dropped w/ rationale,
      IDAQuadAlloc/FreeVectors (+4/−5 lrw asymmetry preserved),
      IDASensAlloc/FreeVectors ((5Ns+1) and maxcol*Ns alloc-count
      quirks; free uses maxord_alloc; Satol* lrw-only), IDAQuadSens
      Alloc/FreeVectors (savrhsQ missing-return C bug vanishes; free
      uses maxord NOT maxord_alloc; both *QSMallocDone reset
      unconditionally); IDAInitialSetup (full IDAS version:
      rhsQ→phiQ[1] eval, QSRHS logs IDA_QSRHS_FAIL but returns
      IDA_QRHS_FAIL quirk preserved, rhsQSDQ NULL_RHSQ/NULL_P checks,
      errconQ/S/QS=false else-branches, ISM_CONSTR check, idaNlsInit
      always + NLSsim/NLSstg init; rhsQS dispatch detaches phi[0..1]/
      phiS[0..1]/phiQ[1]/phiQS[1]/tempv1/tempv2/tmpS3 then DQ-vs-user
      split); import line patched to add idaNlsInit; full ewt-set
      family (12 fns): IDAEwtSet/SS/SV + IDAQuadEwtSet/SS/SV +
      IDAQuadSensEwtSet/EE/SS/SV take &IDAMem with LOCAL scratch
      replacing C's tempv1/ypQ/tempvQS[0] (behaviorally identical,
      serial kernels inlined one-op-per-loop); IDASensEwtSet/EE/SS/SV
      take &mut IDAMem (EE dispatches user efun needing &mut
      user_data) — IDASolve call site patched to also detach phiS[0],
      InitialSetup site written likewise. PINNED forward signature
      (Part 3 must match): IDAQuadSensRhsInternalDQ(&mut IDAMem, Ns:
      i32, t: f64, yy/yp: &NVector, yyS/ypS: &[NVector], rrQ:
      &NVector, resvalQS: &mut [NVector], yytmp/yptmp/tmpQS: &mut
      NVector) -> i32. File ≈ 3790 lines, covers idas.c ≈ 1–5570.
      NEXT: IDAStopTest1/2 (aliasing copy-backs per donor pins),
      IDAHandleFailure (all IDAS cases + SUN_ERR_ARG_CORRUPT →
      MSG_NLS_INPUT_NULL; needs sundials_errors::SUN_ERR_ARG_CORRUPT
      import), IDAGetSolution (cvals+dvals recurrence — VERIFY
      ida_dvals field exists in idas_impl.rs before writing), norm
      family (IDAWrmsNorm/IDASensWrmsNorm pub, IDAQuadSensWrmsNorm,
      three NormUpdates — max-of semantics; QuadWrmsNormUpdate C is
      MAYBE_UNUSED ida_mem), then Part 2: IDAStep/IDASetCoeffs/
      IDACheckConstraints/predictors/IDAQuadNls/IDAQuadSensNls/
      TestError quartet/IDARestore/IDAHandleNFlag/IDAReset/
      IDACompleteStep (donor decision: IDANls/IDASensNls collapsed
      Newton loops live in idas_nls.rs/idas_nls_stg.rs, IDAStep
      dispatches), then Part 3: Rcheck1-3/IDARootfind/IDASensResDQ/
      Res1DQ/QuadSensRhs(1)InternalDQ. Donor Part-2 pin list + imports
      (IDANls, SUNMAX/SUNMIN/SUNRdifferentsign/SUNRpowerR,
      SUN_ERR_ARG_CORRUPT) captured from donor ida.rs ≈ 1400.
      DONE 2026-07-08 (cont. 3) — PART 1 COMPLETE: IDAStopTest1
      (tstop past/at/approaching; NORMAL tout==tretlast shortcut +
      explicit GetSolution paths; ONE_STEP GetSolution; match with
      `_ => IDA_ILL_INPUT` for C's unreachable default) and
      IDAStopTest2 (same tstop block; NORMAL explicit GetSolution;
      ONE_STEP is the alias-reliant branch → copies ida_yy→yret,
      ida_yp→ypret per the owned-stand-in convention);
      IDAHandleFailure(&mut IDAMem, sflag) with all 17 IDAS arms
      (ERR/CONV two-arg messages; QRHS/QSRHS_FAIL mirror C's surplus
      tn arg — harmless under ida_msg_g's substitute-if-present
      semantics; MEM_NULL → IDAProcessError(None,..); 
      SUN_ERR_ARG_CORRUPT → logs MSG_NLS_INPUT_NULL, returns
      IDA_MEM_NULL; default → UNRECOGNIZED_ERROR with
      report-to-sundials-users string); mid-file `use
      crate::sundials_errors::SUN_ERR_ARG_CORRUPT;` added (donor
      Part-2 pattern); IDAGetSolution(&IDAMem, t, yret, ypret) with
      LOCAL cvals/dvals ([ZERO; MXORDP1]) — the "verify ida_dvals
      field" question resolved: donor IDAInit allocates no
      cvals/dvals/Xvecs/Zvecs, they are not IDAMem fields, locals per
      dropped-fused-ops convention; VECTOROP_ERR branches vanish;
      norm family: IDAWrmsNorm (pub, N_VWrmsNormMask/N_VWrmsNorm),
      IDASensWrmsNorm (pub, local cvals buffer + max scan),
      IDAQuadSensWrmsNorm (priv), IDAQuadWrmsNormUpdate (priv,
      _ida_mem for C MAYBE_UNUSED), IDASensWrmsNormUpdate (pub),
      IDAQuadSensWrmsNormUpdate (priv). File ≈ 4350 lines. Two slot
      markers remain: PART 2 slot (above IDAGetSolution): IDAStep,
      IDASetCoeffs, predictors, IDAQuadNls/IDAQuadSensNls,
      IDACheckConstraints, TestError quartet, IDARestore,
      IDAHandleNFlag, IDAReset, IDACompleteStep (idas.c ≈ 6100–7050);
      PART 3 slot (end of file): Rcheck1-3, IDARootfind, IDASensResDQ,
      IDASensRes1DQ, IDAQuadSensRhs(1)InternalDQ (idas.c ≈ 7310–8800;
      IDAQuadSensRhsInternalDQ signature pinned in the marker).
      Part-2 additional needs: `use crate::sundials_math::{SUNMAX,
      SUNMIN, SUNRdifferentsign, SUNRpowerR};` (SUNRsqrt too for Part
      3); IDANls/IDASensNls dispatch into idas_nls.rs/idas_nls_stg.rs
      (donor decision); IDAHandleNFlag ncfn/netf field-pointer args →
      take/restore scalar pattern at call sites. All remaining C text
      banked in transcript.
      DONE 2026-07-08 (cont. 4) — PART 2 COMPLETE (idas.c ≈ 6100–7050
      inserted into the Part-2 slot above IDAGetSolution, C order):
      mid-file imports (idas_nls::IDANls, idas_nls_stg::IDASensNls,
      sundials_math::{SUNMAX, SUNMIN, SUNRdifferentsign, SUNRpowerR,
      SUNRsqrt} — last two ahead of Part 3); IDAStep (sensi_stg/
      sensi_sim; nst==0 init; step loop with SetCoeffs → tstop clamp
      → Predict (+SensPredict via yySpredict/ypSpredict take/restore)
      → IDANls → constraints → TestError; failure paths IDARestore +
      IDAHandleNFlag with ncfn/netf take/restore — quad, staggered-
      sens, and quad-sens paths all pass ncfnQ/netfQ, sens quirk
      commented C-exact; staggered res eval via unwrap+destructure
      into delta; nst==0→IDAReset; final IDACompleteStep +
      ee.scale_inplace(ck)); IDASetCoeffs (integer SUNMIN written
      out; psi/beta/alpha/sigma/gamma recurrence C-literal; alphas −=
      ONE/(i+1) as f64; ck via SUNRabs+SUNMAX; phi→phi-star scale
      loops incl. phiQ/phiS/phiQS); IDACheckConstraints (ConstrMask
      via destructure; correction kernels inline — Compare ONEPT5,
      Prod, Div, Scale −PT1→tempv3, VLin1 z=b*y+x order, Prod mm;
      vnorm<=epsNewt → constraint_corrections+=1 + three-step ee
      update; else step_constraint_fails/constraint_fails, hmin*ONEPSM
      and max_constraint_fails → CONSTR_FAIL; eta=PT9*N_VMinQuotient
      clamps, IDARestore, phase=1, hh*=eta, nst==0→IDAReset,
      PREDICT_AGAIN; non-prefixed C field names used);
      IDAPredict/IDAQuadPredict (&mut IDAMem, destructured fields,
      unit-coefficient combination = copy + accumulations, gamma
      combination) and IDASensPredict/IDAQuadSensPredict (&IDAMem +
      &mut [NVector] targets — both C call sites pass the *predict
      rows, detached by callers); IDAQuadNls (rhsQ destructure
      dispatch, nrQe, savrhsQ copy under quadr_sensi, VDiff/scale
      ONE/cj/VSum kernels); IDAQuadSensNls (yyQS+tempvQS(=ypQS)
      take → QuadSensPredict → nine-field take for rhsQSDQ split
      (DQ arm IDAQuadSensRhsInternalDQ, user arm rhsQS w/ user_data),
      restore, nrQSe, ret staging so all paths restore, corrections
      a=ONE/cj b=−ONE/cj per-is + VSum); TestError quartet (sigma-
      scaled enorms; delta/ypQ/deltaS/yyQS as C's scratch renames;
      aliased VSum d=p+d for km2; check_for_reduction flow C-exact
      incl. errQ_km2/errS_km2/errQS_km2 not setting the flag and the
      knew!=kk cancellation; ck*enorm>ONE→ERROR_TEST_FAIL);
      IDARestore (psi rollback while-loop, ONE/beta rescale of
      phi/phiQ/phiS/phiQS); IDAHandleNFlag (&mut i64/&mut i32 counter
      params; nonrecoverable 6-way chain else NLS_FAIL; recoverable
      maxncf/hmin gate → REP_*_ERR else CONV_FAIL, eta_cf clamp;
      ERROR_TEST_FAIL arms nef==1 (err_knew select, kk=knew, PT9*
      SUNRpowerR(TWO*err_knew+PT0001, −ONE/(kk+1)), eta_min_ef/
      eta_low/hmin clamps), nef==2, nef<maxnef→kk=1, else
      IDA_ERR_FAIL); IDAReset (psi[0]=hh, eta scale of phi[1]/
      phiQ[1]/phiS[1][is]/phiQS[1][is]); IDACompleteStep (nst/kdiff/
      kused/hused; phase-0 doubling w/ hmax_inv clamp; `let action;`
      definite-assignment replaces UNSET sentinel (avoids
      unused_assignments warning); order-raise err_kp1 estimate w/
      VDiff into tempv1/ypQ/ypS/tempvQS + the three *WrmsNormUpdate
      calls; terr chain → LOWER/MAINTAIN/RAISE; eta_max_fx/eta_min_fx
      clamps C-exact; ee→phi[kused+1] saves (copies) for
      phi/phiQ/phiS/phiQS; X+=Z cascade with SEQUENTIAL-ALIASING
      semantics — j=0 phi[ku]+=ee then j=1..=ku
      phi[ku−j]+=phi[ku−j+1] via split_at_mut so each pair reads the
      previous update, per-is for phiS/phiQS). File ≈ 5900 lines.
      REMAINING in idas.c: Part 3 only (EOF slot): IDARcheck1/2/3,
      IDARootfind, IDASensResDQ, IDASensRes1DQ,
      IDAQuadSensRhs(1)InternalDQ (idas.c ≈ 7310–8800; C banked;
      IDAQuadSensRhsInternalDQ signature pinned in the marker).
      DONE 2026-07-08 (cont. 5) — PART 3 COMPLETE, idas.c FULLY
      PORTED. C text re-read fresh from ../sundials-7.7.0 (tail 1600)
      rather than transcript. Rootfinding: IDARcheck1 (iroots zero,
      tlo/ttol, gfun eval via unwrap+destructure — glo/ghi pass as
      &mut Vec through deref coercion; zroot→gactive[i]=false;
      smallh=SUNMAX(ttol/|hh|,PT1)*hh; yy = smallh*phi[1]+phi[0]
      (VLin1 kernel, writes owned ida_yy); reactivation scan);
      IDARcheck2 (irfnd gate; IDAGetSolution via yy/yp take/restore
      — GetSolution reads only phi/psi so defaults-in-place safe;
      smallh sign select; (tplus−tn)*hh>=0 → aliased VLin1
      yy=hratio*phi[1]+yy else GetSolution(tplus); CLOSERT on
      double-zero, glo refresh on zero→nonzero, RTFOUND);
      IDARcheck3(&mut, tout, itask) (thi select; GetSolution
      take/restore; Rootfind; gactive reactivation from grout;
      tlo=trout + glo←grout; root → GetSolution(trout)+RTFOUND);
      IDARootfind (Illinois: full doc comment carried; rootdir[i] as
      f64 * glo[i] for C's int*double; imax: usize; alph=ONE init
      comment carried; side/sideprev; tmid secant + two
      fracint/fracsub inward clamps (FIVE/PT1/HALF); per-iteration
      GetSolution take/restore + gfun→grout; sgnchg→side=1 thi=tmid,
      zroot→break, else side=2 tlo=tmid; final trout/grout/iroots
      block with BOTH zero-test and sign-change-test assignments,
      C-verbatim). DQ residuals: IDASensResDQ (pub; ida_mem moved to
      front replacing void* user_dataS, loops Res1DQ); IDASensRes1DQ
      (_Ns MAYBE_UNUSED; del=SUNRsqrt(SUNMAX(rtol,uround));
      pbar/plist/psave; Delp/Dely families lowercased,
      Del/rDel/r2Del → delta/rdelta/r2delta; DQrhomax==0 vs ratio
      switching → CENTERED1/2 FORWARD1/2 via match with `_ => {}`;
      res bound once (fn ptr Copy) then called with param temps +
      &mut ida_mem.ida_user_data; KERNEL-EXACT scalar N_VLinearSum
      dispatch honored — b==ONE→VLin1 z=a*x+y, a==−b→VScaleDiff
      z=c*(x−y) (NOT a*x+b*y — FP differs), VSum; early error
      returns skip psave restore, C-exact);
      IDAQuadSensRhsInternalDQ (private — C static; pinned signature
      honored; loops Rhs1); IDAQuadSensRhs1InternalDQ (nfel: i64 →
      nrQeS += nfel only on success path — early returns skip BOTH
      psave restore and counter, C-exact; CENTERED1/FORWARD1 only;
      FORWARD1 reuses/overwrites `rdel = ONE/Del` — `let mut rdel`;
      VLin1/VScaleDiff kernels commented). IDAProcessError closes
      idas.c and lives in idas_impl.rs — EOF note in file. NOTE for
      the nls/ls units: the fused N_V*VectorArray serial kernels have
      NO special-case dispatch (plain z=a*x+b*y), unlike scalar
      N_VLinearSum — conventions already applied throughout idas.rs.
      NEXT: idas_ic.c.)
- [x] idas_ic.c — PORTED 2026-07-08 (idas_ic.rs ≈ 1660 lines; NOT
      registered in lib.rs / NOT compile-verified — registers together
      with idas.rs once idas_nls*/idas_ls land. Structural donor:
      ida_rs/src/ida_ic.rs (verified Phase 4). Workspace aliases
      referenced directly, never stored: scalar mc=ee, dtemp=phi[3],
      ynew=tempv2, ypnew=ee, delnew=phi[2]; sens savresS=phiS[2],
      delnewS=phiS[3], yyS0new=phiS[4], ypS0new=eeS, tmpS1=tempv1,
      tmpS2=tempv2. IDACalcIC single fn per one-fn-per-C-fn (both nwt
      loops: scalar 'nwt + staggered 'nwt2 with IDASensNlsIC/ncfnS;
      yy0/yp0 = NVector::new + phi copies; yyS0/ypS0 =
      (0..ns).map(NVector::new).collect + phiS[0..1] copies; sensi_sim
      ypnorm IDASensWrmsNormUpdate; IDASensEwtSet call sites detach
      yyS0+ewtS; frees = NVector::default()/Vec::new(); staggered res
      eval C-leaks temporaries — fields just retain clones, comment).
      ic_efun donor-verbatim. Four resS dispatch shapes, all
      take/detach + resSDQ split (DQ arm IDASensResDQ w/ leading &mut
      IDAMem, user arm resS w/ &mut ida_user_data): (a) IDANlsIC and
      (c) IDASensNlsIC on yy0/yp0/delta→deltaS; (b) IDAfnorm on
      ynew/ypnew/savres→delnewS via phiS.split_at_mut(4)
      (yyS0new=&hi[0], delnewS=&mut lo[3]); (d) IDASensfnorm on
      yy0/yp0/delta→delnewS same split, tempv2 detaches normally.
      *** UPSTREAM ALIASING DEFECT (documented deviation): C's
      IDAfnorm passes tempv2 BOTH as ynew(yy) and tmpS2(yptemp
      scratch); internal-DQ CENTERED1 scratch writes clobber ynew
      mid-eval in C. Inexpressible under borrow rules — Rust uses
      fresh scratch NVector::new(tmp1.len()) = intended semantics.
      Prominent comment in file; revisit only if a sim-sensitivity
      IDACalcIC example diffs. *** phiS[3]→phiS[2] savresS copies via
      split_at_mut(3) + copy_from_slice. lsetup/lsolve donor dispatch
      (LsModule::Ls take/restore; lsolve-missing → -1, lsetup-missing
      → 0); per-is lsolve loops break on r!=0 then map (C-identical);
      IDASensNewtonIC/IDASensfnorm pass ida_delta as rescur (C-exact,
      noted). IDANewtonIC sens lsolve+SensWrmsNormUpdate placed BEFORE
      sysindex rescale per C. IDANewyyp `let mut retval` + trailing
      sensi_sim IDASensNewyyp; IDASensNewyyp per-is
      N_VProd/VLin1/linear_sum_with(-ONE,ONE,..)/VLin1 kernel-exact.
      IDASensNlsIC nj 1..=2 w/ lsetup retry gated nj==1 + trailing
      IDA_SUCCESS per C. IDAICFailFlag func string "IDAICFailFlag"
      per C __func__ (donor slipped w/ "IDACalcIC"); `_ => -99`.
      PINS for coming units: idas_ls_impl.h must add
      Ls(Box<IDALsMem>) variant to idas_impl's LsModule; idas_ls.rs
      must provide idaLsSetup(ida_mem, idals_mem, &yy, &yp, &res)->i32
      and idaLsSolve(ida_mem, idals_mem, &mut b, &ewt, &ycur, &ypcur,
      &rescur)->i32 (donor ida_ls signatures); idas_nls.rs must
      provide ida_has_lsetup(&IDAMem)->bool. NEXT: idas_io.c.)
- [x] idas_io.c — PORTED 2026-07-08 (idas_io.rs ≈ 1510 lines, all ≈88
      C functions in C order; NOT registered in lib.rs / NOT
      compile-verified — registers together with idas.rs /
      idas_ic.rs once idas_nls*/idas_ls land. Structural donor:
      ida_rs/src/ida_io.rs (verified Phase 4); shared IDA functions
      donor-verbatim. Conventions: void*-NULL checks vanish; Get
      out-params &mut i64/i32/f64/NVector; N_Vector-NULL semantics →
      Option (SetId/SetConstraints/GetConsistentIC Option<&mut
      NVector>, GetSensConsistentIC Option<&mut [NVector]>);
      IDAGetCurrentY/Yp/YSens/YpSens + IDAGetUserData are
      borrow-returns. IDAS deltas vs donor: IDASetMaxNonlinIters
      dispatches sensi_sim → NLSsim else NLS (is_none →
      IDAProcessError(None, IDA_MEM_FAIL) per C NULL mem, else
      .as_mut().unwrap().set_max_iters); IDASetSensMaxNonlinIters
      same on NLSstg; IDASetQuadErrCon (quadMallocDone gate, C
      passes NULL mem — None preserved); IDASetSensDQMethod /
      IDASetSensErrCon; IDASetQuadSensErrCon (C quirk preserved:
      quadSensMallocDone arm uses MSG_NO_SENSI, not
      MSG_NO_QUADSENSI). *** IDASetSensParams(p/pbar/plist all
      Option<&[..]>) — DEVIATION + FSA WATCH ITEM: C stores the
      USER'S pointer (ida_p = p) so DQ perturbations of
      ida_p[which] are visible to the user's res through its own
      parameter block; port copies into owned Vec
      (p.to_vec()/Vec::new()). Rust example ports must route res
      parameter reads so the perturbation is seen — open point
      shared with cvodes_rs, settle at FSA example verification.
      *** pbar None→ONE fill / Some→ZERO check BAD_PBAR + SUNRabs;
      plist None→identity / Some→<0 check BAD_PLIST. Quad/QuadSens/
      Sens output getters gate on quadr/quadr_sensi/sensi
      (NO_QUAD/NO_QUADSENS/NO_SENS); IDAGetQuadErrWeights copies
      only if errconQ, QuadSens weights only if errconQS, but
      IDAGetSensErrWeights copies UNCONDITIONALLY per C (comment).
      C quirk preserved + commented: IDAGetNumStepSensSolveFails
      has NO sensi gate and reads ida_ncfn (NOT ncfnS).
      IDAGetNumBacktrackOps `nbacktr as i64`. IDAPrintAllStats
      (outfile: &mut dyn std::io::Write, fmt: SUNOutputFormat) =
      donor base (LS block gated on the pinned LsModule::Ls variant
      — forward reference to idas_ls_impl.h — with
      nje/nreDQ/npe/nps/nli/ncfl/njtsetup/njtimes + three
      per-NLS-iter ratios; NLS-iters-per-step keeps C's nre/nst)
      PLUS IDAS blocks after root stats in C order: quadr →
      nrQe/netfQ; sensi → nrSe/nreS/netfS + (ism==STAGGERED →
      nniS/nnfS/ncfnS) + nsetupsS; quadr_sensi → nrQSe/netfQS.
      IDAGetReturnFlagName(flag: i64) -> String, donor structure +
      IDAS flags NO_SENS..REP_QSRHS_ERR and IDAA flags
      NO_ADJ/BAD_TB0/REIFWD_FAIL/FWD_FAIL/GETY_BADT/NO_BCK/NO_FWD in
      C order, then NLS_SETUP_FAIL/NLS_FAIL, `_ => "NONE"`. Private
      sunfprintf_real/long donor-verbatim (SUN_TABLE_WIDTH=29;
      CSV "% .15e" space-for-plus via fmt_e(v,0,15); TABLE
      fmt_g(v,0,15)). NEXT: idas_ls_impl.h.)
- [x] idas_ls_impl.h — PORTED 2026-07-08 (idas_ls_impl.rs ≈ 455
      lines; REGISTERED in lib.rs (pub mod + flat re-export) — depends
      only on already-compiled modules (idas_impl::IDAResFn + core),
      so the crate should stay green; compile check on Nash's side
      pending. Structural donor: ida_rs/src/ida_ls_impl.rs (verified
      Phase 4); forward half donor-verbatim: IDALS_* return codes,
      IDALsJacFn/PrecSetupFn/PrecSolveFn/JacTimesSetupFn/JacTimesVecFn
      (PrecSetupFn keeps the donor's extra ewt/hh args for the
      pure-Rust re-borrow problem, rationale comment carried),
      IDALsMem struct (ycur/ypcur/rcur borrowed-pointer fields
      dropped — passed as arguments per pinned idaLsSetup/idaLsSolve
      signatures; J_data/jt_data/pdata dropped — user_data from
      IDAMem; pfree → RAII; setup_disabled flag carries C's
      lsetup-NULLing), MSG_LS_* messages with SUN_FORMAT_G expanded
      to %.15g. PrecModule enum has None/User only; the
      BBDPre(Box<idas_bbdpre_impl::IBBDPrecData>) variant lands with
      the idas_bbdpre units (placeholder comment, mirrors the old
      LsModule::Ls precedent). IDAS additions: IDALS_NO_ADJ=-101 /
      IDALS_LMEMB_NULL=-102; PART II backward function types verified
      against include/idas/idas_ls.h — IDALsJacFnB/BS,
      IDALsPrecSetupFnB/BS, IDALsPrecSolveFnB/BS,
      IDALsJacTimesSetupFnB/BS, IDALsJacTimesVecFnB/BS (N_Vector*
      yS_1d/ypS_1d → &[NVector]; user_dataB from IDABMem; backward
      pset types keep C shapes — no ewt/hh extension unless a
      backward-pset example needs it); #[derive(Default)] IDALsMemB
      with ten Option callbacks, P_dataB dropped; PART II messages
      MSG_LS_CAMEM_NULL/LMEMB_NULL/BAD_T/BAD_WHICH/NO_ADJ. ALSO:
      idas_impl.rs LsModule placeholder comment replaced by the
      pinned Ls(Box<crate::idas_ls_impl::IDALsMem>) variant — the
      forward references in idas.rs/idas_ic.rs/idas_io.rs
      (LsModule::Ls dispatch, IDAPrintAllStats LS block) now resolve
      once those register. NEXT: idas_ls.c.)
- [x] idas_ls.c — done → `crates/idas_rs/src/idas_ls.rs` (≈2930 lines,
  registered in lib.rs + flat re-export; compile + clippy `-D warnings`
  + tests green, 8/8 idas_rs).
  PART I (forward): donor-verbatim from ida_rs/ida_ls.rs (verified Phase 4),
  incl. the nli_inc==0 resid-copy fix in idaLsSolve and the pinned hook
  signatures (idaLsInit/Setup/Solve/Perf carry IDALsMem explicitly; callers
  detach y/yp/r/weight via take; tmp vectors = ida_tempv1/2/3).
  IDAS deltas vs donor:
  * IDAGetLinWorkSpace: `*leniwLS = 34` (IDAS C; donor IDA counts 33).
  * Three ida_bbdpre dispatch points (PIN SATISFIED by the idas_bbdpre
    unit 2026-07-11): (1)/(2) idaLsInitialize setup_disabled and the
    idaLsSetup psetup branch were written variant-agnostically
    (`matches!(prec_module, None | User)` / negation) and needed no
    edit; (3) idaLsPSetup carries the IDABBDPrecSetup arm and
    idaLsSolveIterative the IDABBDPrecSolve psolve arm + has_psolve
    contribution (both donor-verbatim from ida_ls.rs).
  PART II (backward): full port following the pinned cvodes_ls.rs design —
  wrappers are forward-callback-typed fns whose &mut UserData downcasts to
  the OUTER IDAMem (idaLs_AccessIDAMem); IDAB_mem.ida_lmem holds
  Box<IDALsMemB> via dyn Any (idaLsB_downcast); idaLs_AccessLMemB/BCur
  return Result<usize, i32> Vec indices; ia_yyTmp/ia_ypTmp (+S) taken as
  owned locals around user-callback calls; IDA-specific ia_noInterp gate
  preserved in all 10 wrappers; BS wrappers pass the (possibly stale)
  yyS/ypS tmps when !interpSensi, as C does. All 12 exported B routines
  (SetLinearSolverB, JacFnB/BS, EpsLinB, LSNormFactorB,
  LinearSolutionScalingB, IncrementFactorB, PreconditionerB/BS,
  JacTimesB/BS, JacTimesResFnB) + idaLsFreeB.
  ***PIN (idaa.c unit): idaLsGetY is a forward-reference bridge whose body
  is `unreachable!` — statically unreachable until idaa.rs can construct an
  IDAadjMem — and MUST be replaced with the ia_interpType dispatch to
  IDAAhermiteGetY / IDAApolynomialGetY, mirroring cvodes_ls.rs cvLsIMget.
  (Now pub(crate): the idas_bbdpre.rs IDAAglocal/IDAAgcomm wrappers call
  the same bridge, so the one replacement fixes both modules.)***
  Tests: donor's 6 carried (defaults, direct-needs-matrix, dense DQ
  Jacobian, direct solve round-trip incl. cjratio scaling, DQ Jtimes vs
  analytic, flag names) + new backward_set_routines_require_adjoint
  (IDALS_NO_ADJ guards). Crate-level clippy-1.94 allows added to idas_rs
  lib.rs with rationale (unnecessary_unwrap / needless_borrow /
  explicit_auto_deref / ptr_arg / field_reassign_with_default — all on
  C-faithful constructs).
- [x] idas_nls.c — done → `crates/idas_rs/src/idas_nls.rs` (donor
      ida_nls.rs verbatim + IDAS deltas: idaNlsLSetup clears ida_forceSetup
      and resets ida_ssS = TWENTY; IDANls carries the idas.c driver —
      sensi_sim = sensi && ism==IDA_SIMULTANEOUS dispatches the composite
      solve to idas_nls_sim::idaNlsSolveSensSim on NLSsim, ida_forceSetup
      joins the lsetup decision, cj!=cjlast also sets ssS=HUNDRED, and the
      final correction updates yyS/ypS via per-vector N_VLinearSum loops.
      ida_forceSetup field ADDED to idas_impl.rs (was missing; C
      idas_impl.h:479) — IDAInit/IDAReInit already cleared it in the banked
      idas.rs. 3 donor tests carried (+ forceSetup/ssS delta assertions).)
- [x] idas_nls_sim.c — done → `crates/idas_rs/src/idas_nls_sim.rs`
      (SIMULTANEOUS corrector on [ee, eeS] following the cvodes_nls_sim.rs
      pinned senswrapper design: ypredictSim/ycorSim/ewtSim aliases NOT
      stored — module reads yypredict/yySpredict, ee/eeS, ewt/ewtS
      directly; only NewtonSolver.deltaS (Ns+1 sub-vectors, sub-vector 0 =
      state) is a real senswrapper; composite WRMS = MAX over sub-norms,
      init 0, state first. idaNlsResidualSensSim uses the pinned resS
      dispatch (resSDQ → crate::idas::IDASensResDQ with IDA_mem replacing
      the user_dataS self-pointer; tmpS1/tmpS2 alias tempv1/tempv2, tmpS3
      real); nrSe incremented at the call site; lsolve solves state system
      with ewt then Ns sens systems with ewtS[is] (last nonzero retval
      mapped after the restore block — C maps every call identically).
      2 tests: setter rejections/defaults + full IDANls sensi_sim solve.)
- [x] idas_nls_stg.c — done → `crates/idas_rs/src/idas_nls_stg.rs`
      (STAGGERED corrector on [eeS] only (Ns sub-vectors, no state slot) +
      idas.c IDASensNls driver: callLSetup always SUNFALSE, failure bumps
      ncfnS, success updates yyS/ypS. Staggered deltas vs sim: lsetup
      counts nsetupsS (not nsetups), does NOT clear forceSetup, and passes
      ida_delta as the residual + tmpS1..3 as temps (NOTE recorded: the
      band-DQ yptemp scratch lands in tempv3 instead of C's tmpS3 — dead
      scratch either way, Jacobian bit-identical); lsolve/residual use
      rescur/resval = ida_delta; the m==0 direct conv test compares
      delnrm <= toldel itself (no PT0001 factor) and the rate estimate
      reads/updates ida_ssS. 2 tests: setter rejections/defaults + state
      IDANls → IDASensNls staggered flow (nsetupsS stays 0).)
      REGISTRATION (this unit): lib.rs now registers idas, idas_ic,
      idas_io, idas_nls, idas_nls_sim, idas_nls_stg + flat re-exports —
      the banked idas.rs/idas_ic.rs/idas_io.rs are compile-verified for
      the first time. Fallout fixed: idas_impl.rs ida_phiS/ida_phiQS
      changed [Vec<NVector>; MXORDP1] → Vec<Vec<NVector>> (matches the
      ida_phi Vec modeling and the maxcol+1/maxord+1 live-row allocation
      in idas.rs; only indexing uses elsewhere); two
      N_VCloneVectorArray call sites in idas.rs (IDASensSVtolerances /
      IDAQuadSensSVtolerances) → map/N_VClone/collect idiom; crate-level
      clippy allows extended with rationale (assign_op_pattern,
      collapsible_if, needless_return, manual_memcpy, if_same_then_else —
      the idas.c order-decision ladder's twin MAINTAIN arms,
      manual_range_contains). Full gate green 2026-07-11: check + clippy
      -D warnings (core+idas) + workspace build + workspace tests (15
      idas_rs tests: 8 prior + 7 new NLS).
- [x] idas_bbdpre_impl.h — done → `crates/idas_rs/src/idas_bbdpre_impl.rs`
      (donor ida_bbdpre_impl.rs + backward additions: IDABBDLocalFnB /
      IDABBDCommFnB callback types (idas_bbdpre.h) and IDABBDPrecDataB
      {glocalB, gcommB} stored behind IDAB_mem.ida_pmem via dyn Any —
      the C ida_pfree hook is Rust Drop.  MSGBBD_AMEM_NULL /
      MSGBBD_PDATAB_NULL are unused in the C 7.7.0 sources and not
      carried.)
- [x] idas_bbdpre.c — done → `crates/idas_rs/src/idas_bbdpre.rs`
      (PART I donor-verbatim from verified ida_bbdpre.rs — serial
      single-block reduction, zlocal/rlocal copy in/out, PrecModule::
      BBDPre(Box<IBBDPrecData>) installed by IDABBDPrecInit.  PART II
      per the idas_ls.rs pinned design: IDABBDPrecInitB/ReInitB check
      adjMallocDone/which then delegate to the forward routines on
      IDAB_mem.IDA_mem with IDAAglocal/IDAAgcomm installed as the inner
      glocal/gcomm (C malloc-failure branch vanishes); the wrappers are
      forward-callback-typed, downcast &mut UserData to the OUTER IDAMem
      (idaLs_AccessIDAMem), read ia_bckpbCrt + pmem downcast, and take
      ia_yyTmp/ia_ypTmp as owned locals around the user call under the
      ia_noInterp gate (interp via the shared idaLsGetY bridge, now
      pub(crate) — see the idaa.c PIN).  idas_ls integration landed per
      PIN: PrecModule::BBDPre variant in idas_ls_impl.rs, IDABBDPrecSetup
      arm in idaLsPSetup, IDABBDPrecSolve arm + has_psolve in
      idaLsSolveIterative.  4 tests: donor's needs-lmem + DQ-Jacobian/
      solve round-trip, new ReInit reset check, new NO_ADJ guards for
      the B entry points.  Gate green 2026-07-11 (19 idas_rs tests).)
- [x] idaa.c — done → `crates/idas_rs/src/idaa.rs` (~2900 lines; ported
      2026-07-11 following the cvodea.rs pinned modeling.  BOTH standing
      adjoint pins SATISFIED: (1) idas_ls.rs idaLsGetY unreachable!
      replaced by crate::idaa::IDAAgetY — the ia_interpType dispatch to
      IDAAhermiteGetY/IDAApolynomialGetY that also serves the
      idas_bbdpre wrappers; (2) IDASolveB/IDACalcICB(S) install the
      OUTER IDAMem as the nested backward problem's ida_user_data
      transiently around each IDASolve/IDACalcIC call (ownership dance,
      idaa_integrate_backward / idaa_calc_ic_backward), which the
      idas_ls.rs and idas_bbdpre.rs PART II wrappers downcast back out.
      Modeling: ck list → Vec (index 0 = initial ckpnt, last = most
      recent, ck_next walk = idx-1); backward problems → Vec in creation
      order (ida_index == position; C-reverse traversals run forward —
      problems independent); ia_malloc/free/storePnt/getY fn pointers →
      interpType dispatch (frees are RAII: ckpnt delete, bckpb delete,
      dataFree, hermite/polynomial free all vanish into Drop);
      ia_Y/ia_YS = OWNED scratch allocated in IDASolveF (C aliases
      phi/phiS rows; every use overwrites before reading — pinned as
      cvodea.rs ca_Y).  Exported: IDAAdjInit/ReInit/Free, IDASolveF
      (root-return caching, tstop replay, ckpnt every nsteps +
      forceSetup), IDACreateB (nested IDACreate + ida_Ns copy),
      IDAInitB/BS, IDAReInitB, IDASStolerancesB/SV, IDAQuadSS/SV-
      tolerancesB, IDAQuadInitB/BS/ReInitB, IDACalcICB/BS (noInterp
      gate + yyTmp/ypTmp(+S) preload), IDASolveB (ckpnt scan +
      IDAAdataStore second forward pass + per-problem active logic),
      IDAGetB, IDAGetQuadB (nst==0 → phiQ[0]), IDAGetAdjY.  Private:
      ckpntInit/New/AllocVectors/CopyVectors (fused ScaleVectorArray
      copies → N_VScale loops), dataStore (yyTmp/ypTmp detached around
      the ONE_STEP replay), ckpntGet (adjoint detached; idx==0 →
      SetInitStep(h0u)+ReInits, else scalar/array copy + forceSetup),
      hermite + polynomial malloc/storePnt/getY (fused
      LinearCombination(+VectorArray) kernels expanded to elementwise
      accumulation loops per the IDAGetSolution idiom; polynomial DD
      update via split_at_mut in-place linear_sum_with), GettnSolutionYp
      /YpS, findIndex (ilast i64), IDAAres/IDAArhsQ (forward-callback-
      typed; interpSensi branches; yyTmp/ypTmp(+S) taken as owned locals
      around user calls; MSGAM_BAD_TINTERP via ida_msg_g).  Deviation
      notes: IDASolveF earlyret bool+flag → Option<i32> (C discards the
      root-branch IDAGetSolution flag — preserved); MSGAM_BAD_TB0/
      BAD_TBOUT C varargs index has no %d conversion — not printed,
      plain message kept; the IDASolveB back-error message formats the
      problem index inline (cvodea.rs precedent).  3 tests: AdjInit
      validation/state + ReInit/Free, CreateB indices + which checks +
      NO_MALLOC-before-InitB + tolerance propagation + NO_FWD, InitB
      BAD_TB0 + workspace init.  Gate green 2026-07-11 (22 idas_rs
      tests).  ASA numerics verification lands with the idas serial
      examples (idasRoberts_ASAi_dns etc.).)
- [x] idaa_io.c — done → `crates/idas_rs/src/idaa_io.rs` (cvodea_io.rs
      conventions: shared idaa_io_which_index preamble; IDAAdjSetNoSensi;
      11 ***B optional-input wrappers delegating to the nested solver
      (SetNonlinearSolverB via idas_nls, UserDataB stored on IDABMem,
      MaxOrdB/MaxNumStepsB/InitStepB/MaxStepB/SuppressAlgB/IdB/
      ConstraintsB/QuadErrConB via idas_io); outputs: IDAGetAdjIDABmem
      (borrow of nested Box<IDAMem>, None on C NULL paths),
      IDAadjCheckPointRec with my_addr/next_addr as Option<usize> ck_mem
      indices (list walked newest-first = Vec reversed),
      IDAGetAdjCheckPointsInfo, IDAGetConsistentICB, IDAGetUserDataB
      (borrow-back), IDAGetAdjDataPointHermite/Polynomial (C NULL
      outputs → Option), IDAGetAdjCurrentCheckPoint (Option<usize>).
      2 tests: setter delegation + user-data round trip + nested access;
      checkpoint info/current + wrong-interp guard.  Gate green
      2026-07-11 (24 idas_rs tests).)
- [x] idas_cli.c — done → `crates/idas_rs/src/idas_cli.rs` (donor
      ida_cli.rs + IDAS additions, C table order preserved: int keys
      + quad_err_con/sens_err_con/sens_max_nonlin_iters/
      quad_sens_err_con; tworeal + quad_scalar_tolerances; NEW twoint
      table (max_order_b, suppress_alg_b, quad_err_con_b,
      linear_solution_scaling_b); action + sens_toggle_off/adj_no_sensi;
      NEW int-real table (sens_dq_method, init_step_b, max_step_b,
      eps_lin_b, ls_norm_factor_b, increment_factor_b); NEW
      int-real-real table (scalar_tolerances_b,
      quad_scalar_tolerances_b); NEW int-long table (max_num_steps_b);
      default_id "idas"; dispatch order = C (int, long, real, twoint,
      tworeal, action, int+real, int+long, int+real+real); struct-pair
      helpers from sundials_core::sundials_cli per the ida_cli.rs donor
      convention.  6 tests incl. backward-key delegation to the nested
      solver.  Gate green 2026-07-11 (29 idas_rs tests).
      *** PHASE 5 LIBRARY COMPLETE: all idas/*.c and idaa*.c ported and
      compiled; remaining Phase 5 work = idas serial examples
      (idasRoberts_dns etc., incl. FSA/ASA verification and the
      FSA-verification watch items). ***
### idas examples (Phase 5 verification)
- [x] idasRoberts_dns — IDENTICAL (donor idaRoberts_dns.rs; deltas: csv
      filename, two header lines.  First idas example run — validates
      the idas.c driver + idas_nls + idas_ls + rootfinding end to end.
      lib.rs prelude gained the sunnonlinsol_newton/fixedpoint
      re-exports the examples need, matching ida_rs.)
- [ ] idasAnalytic_mels — todo (donor idaAnalytic_mels.rs; IDAS adds
      reltol/abstol header lines + IDAGetActualInitStep h0 output)
- [x] idasHeat2D_bnd — IDENTICAL (donor + header-indent deltas)
- [x] idasHeat2D_kry — IDENTICAL (donor + header-indent deltas)
- [x] idasFoodWeb_bnd — LOCAL-C (byte-identical to the local C build;
      shipped .out foreign-libm last-digit h drift at t>=0.7, same
      signature as the Phase 4 idaFoodWeb_bnd case.  The local C
      reference tree was REBUILT 2026-07-11 and now lives COMMITTED at
      localref/ — verify_examples.sh defaults to it.)
- [x] idasKrylovDemo_ls — IDENTICAL (donor + header-indent deltas)
- [ ] idasSlCrank_dns — todo (donor idaSlCrank_dns.rs)
- [x] idasSlCrank_dns — last-digit (first quadrature-carrying example:
      IDAQuadInit + IDAQuadSStolerances + IDASetQuadErrCon + IDAGetQuad
      all exercised; 44/45 output lines byte-identical to local C; the
      G line differs in the 13th digit — root-caused to the C example
      binary's __sincos_stret fusion, see VERIFICATION.md.  The quad
      error-control path (IDAQuadTestError feeding err_k into eta) was
      instrumented bit-level against a patched local C build: all
      library-side values identical.)
- [x] idasAkzoNob_dns — IDENTICAL (fresh port; 2nd quadrature example,
      G matches shipped .out byte-for-byte incl. the 16-decimal print)
- [x] idasRoberts_FSA_dns — IDENTICAL (first FSA example in the
      workspace.  Shipped ref `-sensi stg t`: byte-identical.  All
      other modes verified byte-identical against the local C build
      (localref variants committed): -sensi sim t / sim f / stg f and
      -nosensi.  Validates end-to-end: IDASensInit, SIMULTANEOUS
      (idas_nls_sim) and STAGGERED (idas_nls_stg) correctors with and
      without sens error control, IDASensEEtolerances,
      IDASetSensParams (p-copy WATCH ITEM SETTLED — byte-exact),
      IDACalcIC with sensitivities incl. IDAGetSensConsistentIC,
      IDAQuadInit + IDAQuadSensInit(None → internal QuadSens DQ),
      IDAGetSens/IDAGetQuad/IDAGetQuadSens, sens stats in
      IDAPrintAllStats, and idas_cli via the Analytic_mels
      idas.init_step variant.  (The IDAfnorm tempv2-aliasing deviation
      is NOT yet exercised: this example supplies an analytic resS, so
      the internal sens-DQ IC path stays cold — revisit with a
      resSDQ-based IC example.)  Harness: verify_examples.sh now runs
      argument-encoded reference variants (name_argv.out, decode
      overrides in tools/verify_args.map) for idas_rs/arkode_rs — also
      auto-verifies idasKrylovDemo_ls_1/_2 and the Analytic_mels
      variant, all IDENTICAL.
- [x] idasSlCrank_FSA_dns — tolerance-level (fresh port; the first
      INTERNAL-DQ sensitivity example, and it EXPOSED AND SETTLED the
      pinned IDASetSensParams p-copy watch item: the DQ perturbs the
      owned ida_p copy, which the user residual never sees, so all
      sensitivities silently came out zero (run bit-matched the
      no-sensi SlCrank_dns).  FIX (pinned, ARCHITECTURE.md §3.6):
      FSAUserData { p, user } wrapper in sundials_types.rs; the DQ
      routines mirror every p[which] perturbation into the user data
      via ida_dq_set_p; IDASolve rejects internal-DQ sensi without the
      wrapper (ILL_INPUT) instead of silently zeroing.  cvodes_rs has
      the same latent defect (task chip spawned).  With the fix the
      example self-validates: FSA dG/dp (3.3346e-01, -3.6375e-01)
      agrees with its own finite-difference checks; the FD sections are
      byte-identical to local C.  Remaining diffs vs local C are stats/
      G-tail/one dG/dp last digit and are tolerance-level only: this C
      example is provably compiler-flag-unstable (232 vs 263 steps
      between -fmath-errno and default builds of the SAME source —
      sincos fusion in the sin/cos-heavy residual, amplified by the DQ's
      4 extra residual calls per step through the near-tie MIN(dely,
      delp) pick).  See VERIFICATION.md.
- [x] idasRoberts_ASAi_dns — IDENTICAL (first adjoint example:
      byte-exact end-to-end validation of idaa.rs + idas_ls PART II —
      IDAAdjInit(HERMITE), two-phase IDASolveF checkpointing over
      t=4e10, IDACreateB/InitB/SStolerancesB/UserDataB/MaxNumStepsB,
      IDASetLinearSolverB + user JacB through the B-wrappers,
      IDAQuadInitB + backward quadrature error control, IDASolveB with
      checkpoint replay (IDAAdataStore hot restarts + ownership dance +
      Hermite interpolation into resB/JacB/rhsQB), IDAGetB/IDAGetQuadB,
      then IDAReInitB + IDACalcICB + IDAGetConsistentICB and a second
      backward solve from TB1.  Also exercises IDAWFtolerances (user
      ewt) and user forward Jac.)
- [x] idasAkzoNob_ASAi_dns — LOCAL-C (byte-identical to the local C
      build; shipped .out is stale — missing the trailing space the
      current C `%24.16f \n` prints and a blank line.  Second adjoint
      validation: dG/dy0 w.r.t. initial conditions, stiff 6-eq DAE.)
- [x] idasHessian_ASA_FSA — IDENTICAL (2nd-order adjoint byte-exact:
      SIMULTANEOUS FSA forward with analytic resS/rhsQS + quad-sens
      error control, then TWO sensitivity-dependent backward problems
      via IDAInitBS/IDAQuadInitBS — exercises the BS wrappers and the
      Hermite SENSITIVITY interpolation (ia_interpSensi path,
      yySTmp/ypSTmp) — plus the FD verification reruns.
      *** PHASE 5 COMPLETE: library + all 13 serial examples verified
      (9 identical incl. arg variants, 2 local-C, 2 last-digit/
      tolerance-level — both root-caused to C-binary sincos fusion).
      Next: Phase 6 (arkode), starting with the deferred Phase-1 core
      files (sundials_stepper etc., see the Phase 1 section). ***
### excluded (idas) examples
idasRoberts_klu, idasRoberts_sps, idasRoberts_ASAi_klu, idasRoberts_ASAi_sps,
idasRoberts_FSA_klu, idasRoberts_FSA_sps (KLU/SuperLU)

## crates/arkode_rs — Phase 6
- [x] arkode_types_impl.h — committed (forward decls subsumed by arkode_impl.rs)
- [x] arkode_impl.h — committed (arkode_impl.rs: arkode.h constants/fn types + ARKodeMem with Option<fn> step_* table, step_mem = Box<dyn Any>, ARKInterp enum, `fn` field renamed fn_)
- [x] arkode_butcher.c — committed (PART I table object; PART II CheckOrder/CheckARKOrder + mv/vv/vp/dot, rowsum/order1..order6s, __ButcherSimplifyingAssumptions; C's duplicated method/embedding blocks factored into byte-identical-output helpers; CheckARKOrder d[1]=B1->d quirk preserved)
- [x] arkode_butcher_erk.c (+ .def) — committed (X-macro -> match; all 27 tables transcribed; validated against ark_test_butcher.out expectations)
- [x] arkode_butcher_dirk.c (+ .def) — committed (27 tables; DIRK+ERK ARK-pair CheckARKOrder validated against ark_test_butcher.out)
- [x] arkode_interp_impl.h — committed (Hermite/Lagrange content structs; generic ARKInterp enum lives in arkode_impl.rs)
- [x] arkode_interp.c — committed (dispatchers take ark_mem only, take/put-back on ark_mem.interp; Hermite quartic/quintic bootstrap recurses on the impl; LBasis family spans full nhist as in C; polynomial-exactness tests)
- [ ] arkode.c — PARTS I-II committed (vector utils, ARKodeGetDky, arkCreate, arkEwtSetSS/SV/SmallReal + arkRwtSet family with donor in-place idiom; efun=None means internal dispatch). PART III adds arkInit + arkCheckTimestepper/Nvector + arkAllocVectors/arkFreeVectors (rwt_is_ewt = rwt left unallocated). PART V adds the full Evolve chain: ARKodeEvolve, arkInitialSetup, arkStopTests, arkCompleteStep, arkHandleFailure, arkCheckConvergence/Constraints/TemporalError, ark_rfun_apply_yn (relaxation call site stubs to ARK_RELAX_MEM_NULL until arkode_relaxation.c lands). PART VI adds ARKodeSStolerances/SVtolerances/WFtolerances + ARKodeFree; arkode_io.rs PART II adds SetUserData/SetInitStep/SetFixedStep/SetMaxNumSteps/SetStopTime/SetInterpolantDegree. Remaining: Reset/Resize/predictors/remaining io Set-Get families
- [ ] arkode_io.c — PARTS I-III committed (SetDefaults, essential setters, arkReplaceAdaptController, ARKodePrintAllStats + sunfprintf helpers). Remaining: remaining Set/Get families, WriteParameters, SetOptions CLI
- [x] arkode_adapt_impl.h — committed (ARKodeHAdaptMem struct + adaptivity constants)
- [x] arkode_adapt.c — committed (arkAdapt signature drops the always-ark_mem-derived hadapt_mem/ycur args per Addendum C.1; SUNRcopysign added to sundials_math)
- [x] arkode_root_impl.h — committed (ARKodeRootMem struct, C int*/realtype* arrays -> Vec)
- [x] arkode_root.c — committed (RootInit/Free/PrintMem/Check1-3/Rootfind; root_mem take/put-back wrappers; root_data alias collapsed to ark_mem.user_data; ARKodeGetDky ported into arkode.rs PART I)
- [x] arkode_erkstep_impl.h — committed (ARKodeERKStepMem; Xvecs pointer array replaced by call-site operand assembly, liw accounting kept; adj_f deferred with the adjoint machinery)
- [ ] arkode_erkstep.c — PART I committed (ERKStepCreate/ReInit, erkStep Init/FullRHS/TakeStep/Resize/Free/PrintMem, SetButcherTable/CheckButcherTable/ComputeSolutions, ApplyForcing/SetInnerForcing; step_mem take/put-back with release around the FullRHS re-entry). Deferred: TakeStep_Adjoint/fe_Adj/CreateAdjointStepper (ManyVector), RelaxDeltaE (relaxation)
- [ ] arkode_erkstep_io.c — PART I committed (erkStep SetDefaults/SetOrder/GetNumRhsEvals/GetEstLocalErrors/GetStageIndex/PrintAllStats/WriteParameters + ERKStepSetTable/Num/Name, GetCurrentButcherTable (copy), GetTimestepperStats). Remaining: deprecated ERKStep* wrapper aliases, SetOptions/SetRelaxFn
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
- [x] arkode_relaxation_impl.h — committed (ARKodeRelaxMem struct + delta-E/get-order fn types)
- [ ] arkode_relaxation.c — todo
- [ ] arkode_user_controller.c/.h — todo
- [ ] arkode_sunstepper.c — todo
- [ ] arkode_cli.c — todo
### excluded (arkode)
xbraid/ (XBraid); fmod_* (Fortran)

## Phase 7 — docs
- [ ] sundials.md — todo
- [ ] final CLAUDE.md / ARCHITECTURE.md pass — todo
