/* -----------------------------------------------------------------
 * Translated from src/arkode/arkode_arkstep.c (ARKODE 7.7.0).
 * ARKStep: additive (IMEX) Runge-Kutta time stepper.
 *
 * PART I: creation/attach ops, arkStep_Init, arkStep_FullRHS,
 * arkStep_TakeStep_Z, Butcher table set/check, predictors,
 * stage setup and solution computation (identity mass matrix
 * paths; the fixed/time-dependent mass branches follow the C
 * shape but the mass-solver ops are not yet ported and mass_type
 * stays MASS_IDENTITY).
 *
 * Deferred with their subsystems:
 *  - arkStep_TakeStep_ERK_Adjoint / arkStep_fe_Adj /
 *    ARKStepCreateAdjointStepper (ManyVector + adjoint modules)
 *  - arkStep_AttachMasssol / arkStep_DisableMSetup /
 *    arkStep_GetMassMem (ARKLS mass half)
 *  - arkStep_RelaxDeltaE / arkStep_SetRelaxFn / arkStep_GetOrder
 *    (relaxation module)
 *  - arkStep_Resize / arkStep_Reset (with arkode.c Resize/Reset)
 *  - arkStep_SetInnerForcing / arkStep_ApplyForcing (MRIStep;
 *    forcing flags stay false until then)
 *
 * step_mem access follows the erkstep take/put-back convention;
 * the TakeStep wrapper releases step_mem around the step_fullrhs
 * op re-entry.  The ARKLS lmem lives on ark_mem.lmem (Addendum
 * C.2), so lsetup/lsolve calls need no step_mem release.
 * -----------------------------------------------------------------*/
use crate::arkode::{arkAllocVec, arkAllocVecArray, arkCreate, arkFreeVec, arkInit};
use crate::arkode_arkstep_impl::*;
use crate::arkode_butcher::{ARKodeButcherTable_IsStifflyAccurate, ARKodeButcherTable_Space};
use crate::arkode_butcher_dirk::ARKodeButcherTable_LoadDIRK;
use crate::arkode_butcher_erk::ARKodeButcherTable_LoadERK;
use crate::arkode_impl::*;
use crate::nvector_serial::*;
use crate::sundials_math::*;
use crate::sundials_types::*;

/*===============================================================
  Exported functions
  ===============================================================*/

/// C ARKStepCreate(fe, fi, t0, y0, sunctx); returns None on the
/// fe == fi == NULL input error (allocation failures cannot occur).
pub fn ARKStepCreate(
    fe: Option<ARKRhsFn>,
    fi: Option<ARKRhsFn>,
    t0: f64,
    y0: &NVector,
    sunctx: &crate::sundials_context::SUNContext,
) -> Option<Box<ARKodeMem>> {
    /* Check that at least one of fe, fi is supplied and is to be used */
    if fe.is_none() && fi.is_none() {
        arkProcessError(
            None,
            ARK_ILL_INPUT,
            line!(),
            "ARKStepCreate",
            file!(),
            "Must specify at least one of fe, fi (both NULL).",
        );
        return None;
    }

    /* Create ark_mem structure and set default values */
    let mut ark_mem = arkCreate(sunctx);

    /* Allocate ARKodeARKStepMem structure, and initialize to zero */
    let step_mem = Box::new(ARKodeARKStepMem::default());

    /* Attach step_mem structure and function pointers to ark_mem */
    ark_mem.step_attachlinsol = Some(arkStep_AttachLinsol);
    ark_mem.step_attachmasssol = Some(arkStep_AttachMasssol);
    ark_mem.step_disablelsetup = Some(arkStep_DisableLSetup);
    ark_mem.step_disablemsetup = Some(arkStep_DisableMSetup);
    ark_mem.step_getlinmem = Some(arkStep_GetLmem);
    ark_mem.step_getmassmem = Some(arkStep_GetMassMem);
    ark_mem.step_setjcur = Some(arkStep_SetJcur);
    ark_mem.step_getimplicitrhs = Some(arkStep_GetImplicitRHS);
    ark_mem.step_getgammas = Some(arkStep_GetGammas);
    ark_mem.step_init = Some(arkStep_Init);
    ark_mem.step_fullrhs = Some(arkStep_FullRHS);
    ark_mem.step = Some(arkStep_TakeStep_Z);
    ark_mem.step_setuserdata = Some(crate::arkode_arkstep_io::arkStep_SetUserData);
    ark_mem.step_printallstats = Some(crate::arkode_arkstep_io::arkStep_PrintAllStats);
    ark_mem.step_writeparameters = Some(crate::arkode_arkstep_io::arkStep_WriteParameters);
    ark_mem.step_setusecompensatedsums = None;
    ark_mem.step_resize = None; /* arkStep_Resize: with arkode.c Resize */
    ark_mem.step_free = Some(arkStep_Free);
    ark_mem.step_printmem = Some(arkStep_PrintMem);
    ark_mem.step_setdefaults = Some(crate::arkode_arkstep_io::arkStep_SetDefaults);
    ark_mem.step_computestate = Some(arkStep_ComputeState);
    ark_mem.step_setoptions = Some(crate::arkode_arkstep_io::arkStep_SetOptions);
    ark_mem.step_setrelaxfn = None; /* arkStep_SetRelaxFn: relaxation module pending */
    ark_mem.step_setorder = Some(crate::arkode_arkstep_io::arkStep_SetOrder);
    ark_mem.step_setnonlinearsolver =
        Some(crate::arkode_arkstep_nls::arkStep_SetNonlinearSolver);
    ark_mem.step_setlinear = Some(crate::arkode_arkstep_io::arkStep_SetLinear);
    ark_mem.step_setnonlinear = Some(crate::arkode_arkstep_io::arkStep_SetNonlinear);
    ark_mem.step_setautonomous = Some(crate::arkode_arkstep_io::arkStep_SetAutonomous);
    ark_mem.step_setnlsrhsfn = Some(crate::arkode_arkstep_nls::arkStep_SetNlsRhsFn);
    ark_mem.step_setdeduceimplicitrhs =
        Some(crate::arkode_arkstep_io::arkStep_SetDeduceImplicitRhs);
    ark_mem.step_setnonlincrdown = Some(crate::arkode_arkstep_io::arkStep_SetNonlinCRDown);
    ark_mem.step_setnonlinrdiv = Some(crate::arkode_arkstep_io::arkStep_SetNonlinRDiv);
    ark_mem.step_setdeltagammamax = Some(crate::arkode_arkstep_io::arkStep_SetDeltaGammaMax);
    ark_mem.step_setlsetupfrequency =
        Some(crate::arkode_arkstep_io::arkStep_SetLSetupFrequency);
    ark_mem.step_setpredictormethod =
        Some(crate::arkode_arkstep_io::arkStep_SetPredictorMethod);
    ark_mem.step_setmaxnonliniters = Some(crate::arkode_arkstep_io::arkStep_SetMaxNonlinIters);
    ark_mem.step_setnonlinconvcoef = Some(crate::arkode_arkstep_io::arkStep_SetNonlinConvCoef);
    ark_mem.step_setstagepredictfn = Some(crate::arkode_arkstep_io::arkStep_SetStagePredictFn);
    ark_mem.step_getnumrhsevals = Some(crate::arkode_arkstep_io::arkStep_GetNumRhsEvals);
    ark_mem.step_getnumlinsolvsetups =
        Some(crate::arkode_arkstep_io::arkStep_GetNumLinSolvSetups);
    ark_mem.step_getcurrentgamma = Some(crate::arkode_arkstep_io::arkStep_GetCurrentGamma);
    ark_mem.step_getestlocalerrors = Some(crate::arkode_arkstep_io::arkStep_GetEstLocalErrors);
    ark_mem.step_getnumnonlinsolviters =
        Some(crate::arkode_arkstep_io::arkStep_GetNumNonlinSolvIters);
    ark_mem.step_getnumnonlinsolvconvfails =
        Some(crate::arkode_arkstep_io::arkStep_GetNumNonlinSolvConvFails);
    ark_mem.step_getnonlinsolvstats = Some(crate::arkode_arkstep_io::arkStep_GetNonlinSolvStats);
    ark_mem.step_setforcing = Some(arkStep_SetInnerForcing);
    ark_mem.step_getstageindex = Some(crate::arkode_arkstep_io::arkStep_GetStageIndex);
    ark_mem.step_supports_adaptive = true;
    ark_mem.step_supports_implicit = true;
    ark_mem.step_supports_massmatrix = true;
    ark_mem.step_supports_relaxation = false; /* SUNTRUE in C; relaxation pending */
    ark_mem.step_mem = Some(step_mem);

    /* Set default values for optional inputs */
    let retval = crate::arkode_arkstep_io::arkStep_SetDefaults(&mut ark_mem);
    debug_assert_eq!(retval, ARK_SUCCESS);

    /* re-take step_mem for the remaining initialization */
    let mut step_mem = arkStep_AccessStepMem(&mut ark_mem, "ARKStepCreate").unwrap();

    /* Set implicit/explicit problem based on function pointers */
    step_mem.explicit = fe.is_some();
    step_mem.implicit = fi.is_some();

    /* Allocate the general ARK stepper vectors using y0 as a template */
    /* NOTE: Fe, Fi, cvals and Xvecs will be allocated later on
       (based on the number of ARK stages) */

    /* Clone the input vector to create sdata, zpred and zcor */
    let tmpl_len = y0.data.len();
    arkAllocVec(&mut ark_mem, tmpl_len, &mut step_mem.sdata);
    arkAllocVec(&mut ark_mem, tmpl_len, &mut step_mem.zpred);
    arkAllocVec(&mut ark_mem, tmpl_len, &mut step_mem.zcor);

    /* Copy the input parameters into ARKODE state */
    step_mem.fe = fe;
    step_mem.fi = fi;

    /* Update the ARKODE workspace requirements */
    ark_mem.liw += 41; /* fcn/data ptr, int, long int, sunindextype, sunbooleantype */
    ark_mem.lrw += 10;

    /* If an implicit component is to be solved, create default Newton NLS object */
    step_mem.ownNLS = false;
    let implicit = step_mem.implicit;

    /* Initialize initial error norm  */
    step_mem.eRNrm = ONE;

    /* (counters, fused-op workspace, forcing data and the fn_implicit
       alias were zeroed at construction) */

    ark_mem.step_mem = Some(step_mem);

    if implicit {
        let NLS = crate::sunnonlinsol_newton::SUNNonlinSol_Newton(y0, sunctx);
        let retval = crate::arkode_io::ARKodeSetNonlinearSolver(&mut ark_mem, NLS);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(&ark_mem),
                ARK_MEM_FAIL,
                line!(),
                "ARKStepCreate",
                file!(),
                "Error attaching default Newton solver",
            );
            return None;
        }
        if let Some(mut sm) = arkStep_AccessStepMem(&mut ark_mem, "ARKStepCreate") {
            sm.ownNLS = true;
            ark_mem.step_mem = Some(sm);
        }
    }

    /* (linear/mass solver addresses already NULL from construction) */

    /* Initialize main ARKODE infrastructure */
    let retval = arkInit(&mut ark_mem, t0, y0, FIRST_INIT);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(&ark_mem),
            retval,
            line!(),
            "ARKStepCreate",
            file!(),
            "Unable to initialize main ARKODE infrastructure",
        );
        return None;
    }

    Some(ark_mem)
}

/*---------------------------------------------------------------
  ARKStepReInit:

  This routine re-initializes the ARKStep module to solve a new
  problem of the same size as was previously solved. This routine
  should also be called when the problem dynamics or desired solvers
  have changed dramatically, so that the problem integration should
  resume as if started from scratch.

  Note all internal counters are set to 0 on re-initialization.
  ---------------------------------------------------------------*/
pub fn ARKStepReInit(
    ark_mem: &mut ARKodeMem,
    fe: Option<ARKRhsFn>,
    fi: Option<ARKRhsFn>,
    t0: f64,
    y0: &NVector,
) -> i32 {
    /* access ARKodeARKStepMem structure */
    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "ARKStepReInit") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* Check if ark_mem was allocated */
    if !ark_mem.MallocDone {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_NO_MALLOC,
            line!(),
            "ARKStepReInit",
            file!(),
            "Attempt to call before ARKODE initialized.",
        );
        return ARK_NO_MALLOC;
    }

    /* Check that at least one of fe, fi is supplied and is to be used */
    if fe.is_none() && fi.is_none() {
        ark_mem.step_mem = Some(step_mem);
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "ARKStepReInit",
            file!(),
            "Must specify at least one of fe, fi (both NULL).",
        );
        return ARK_ILL_INPUT;
    }

    /* Set implicit/explicit problem based on function pointers */
    step_mem.explicit = fe.is_some();
    step_mem.implicit = fi.is_some();

    /* Copy the input parameters into ARKODE state */
    step_mem.fe = fe;
    step_mem.fi = fi;

    /* Initialize initial error norm  */
    step_mem.eRNrm = ONE;

    ark_mem.step_mem = Some(step_mem);

    /* Initialize main ARKODE infrastructure */
    let retval = arkInit(ark_mem, t0, y0, FIRST_INIT);
    if retval != ARK_SUCCESS {
        arkProcessError(
            Some(ark_mem),
            retval,
            line!(),
            "ARKStepReInit",
            file!(),
            "Unable to reinitialize main ARKODE infrastructure",
        );
        return retval;
    }

    /* Initialize all the counters */
    let mut step_mem = arkStep_AccessStepMem(ark_mem, "ARKStepReInit").unwrap();
    step_mem.nfe = 0;
    step_mem.nfi = 0;
    step_mem.nsetups = 0;
    step_mem.nstlp = 0;
    ark_mem.step_mem = Some(step_mem);

    if let Some(mut lmem) = ark_mem.lmem.take() {
        crate::arkode_ls::arkLsInitializeCounters(&mut lmem);
        ark_mem.lmem = Some(lmem);
    }

    ARK_SUCCESS
}

/*===============================================================
  Interface routines supplied to ARKODE
  ===============================================================*/

/*---------------------------------------------------------------
  arkStep_ComputeState:

  Computes y based on the current prediction and a given
  correction.
  ---------------------------------------------------------------*/
pub fn arkStep_ComputeState(ark_mem: &mut ARKodeMem, zcor: &NVector, z: &mut NVector) -> i32 {
    let step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_ComputeState") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    N_VLinearSum(ONE, &step_mem.zpred, ONE, zcor, z);

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_Free frees all ARKStep memory (drops handled by
  ownership; this keeps the C workspace accounting).
  ---------------------------------------------------------------*/
pub fn arkStep_Free(ark_mem: &mut ARKodeMem) {
    /* conditional frees on non-NULL ARKStep module */
    if let Some(b) = ark_mem.step_mem.take() {
        if let Ok(mut step_mem) = b.downcast::<ARKodeARKStepMem>() {
            /* free the Butcher tables */
            if let Some(bt) = step_mem.Be.take() {
                let (mut bliw, mut blrw) = (0i64, 0i64);
                ARKodeButcherTable_Space(&bt, &mut bliw, &mut blrw);
                ark_mem.liw -= bliw;
                ark_mem.lrw -= blrw;
            }
            if let Some(bt) = step_mem.Bi.take() {
                let (mut bliw, mut blrw) = (0i64, 0i64);
                ARKodeButcherTable_Space(&bt, &mut bliw, &mut blrw);
                ark_mem.liw -= bliw;
                ark_mem.lrw -= blrw;
            }

            /* free the nonlinear solver memory (if applicable) */
            step_mem.NLS = None;
            step_mem.ownNLS = false;

            /* free the linear solver memory */
            if let Some(lfree) = step_mem.lfree {
                lfree(ark_mem);
            }

            /* free the mass matrix solver memory */
            if let Some(mfree) = step_mem.mfree {
                mfree(ark_mem);
            }

            /* free the sdata, zpred and zcor vectors */
            arkFreeVec(ark_mem, &mut step_mem.sdata);
            arkFreeVec(ark_mem, &mut step_mem.zpred);
            arkFreeVec(ark_mem, &mut step_mem.zcor);

            /* free the RHS and stage vectors */
            if !step_mem.Fe.is_empty() {
                for v in step_mem.Fe.iter_mut() {
                    arkFreeVec(ark_mem, v);
                }
                step_mem.Fe.clear();
                ark_mem.liw -= step_mem.stages as i64;
            }
            if !step_mem.Fi.is_empty() {
                for v in step_mem.Fi.iter_mut() {
                    arkFreeVec(ark_mem, v);
                }
                step_mem.Fi.clear();
                ark_mem.liw -= step_mem.stages as i64;
            }
            if !step_mem.z.is_empty() {
                for v in step_mem.z.iter_mut() {
                    arkFreeVec(ark_mem, v);
                }
                step_mem.z.clear();
                ark_mem.liw -= step_mem.stages as i64;
            }

            /* free the reusable arrays for fused vector interface */
            if !step_mem.cvals.is_empty() {
                step_mem.cvals.clear();
                ark_mem.lrw -= step_mem.nfusedopvecs as i64;
            }
            /* (Xvecs assembled at call sites; keep its liw accounting) */
            if step_mem.nfusedopvecs > 0 {
                ark_mem.liw -= step_mem.nfusedopvecs as i64;
            }
            step_mem.nfusedopvecs = 0;

            /* free work arrays for MRI forcing */
            if !step_mem.stage_times.is_empty() {
                step_mem.stage_times.clear();
                ark_mem.lrw -= step_mem.stages as i64;
            }
            if !step_mem.stage_coefs.is_empty() {
                step_mem.stage_coefs.clear();
                ark_mem.lrw -= step_mem.stages as i64;
            }

            /* the time stepper module itself drops here */
        }
    }
}

/*---------------------------------------------------------------
  arkStep_PrintMem:

  This routine outputs the memory from the ARKStep structure to
  a specified stream (for debugging; extra vector output elided
  as in the erkstep port).
  ---------------------------------------------------------------*/
pub fn arkStep_PrintMem(ark_mem: &mut ARKodeMem, outfile: &mut dyn std::io::Write) {
    let step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_PrintMem") {
        None => return,
        Some(sm) => sm,
    };

    /* output integer quantities */
    let _ = writeln!(outfile, "ARKStep: q = {}", step_mem.q);
    let _ = writeln!(outfile, "ARKStep: p = {}", step_mem.p);
    let _ = writeln!(outfile, "ARKStep: istage = {}", step_mem.istage);
    let _ = writeln!(outfile, "ARKStep: stages = {}", step_mem.stages);
    let _ = writeln!(outfile, "ARKStep: maxcor = {}", step_mem.maxcor);
    let _ = writeln!(outfile, "ARKStep: msbp = {}", step_mem.msbp);
    let _ = writeln!(outfile, "ARKStep: predictor = {}", step_mem.predictor);

    /* output long integer quantities */
    let _ = writeln!(outfile, "ARKStep: nstlp = {}", step_mem.nstlp);
    let _ = writeln!(outfile, "ARKStep: nfe = {}", step_mem.nfe);
    let _ = writeln!(outfile, "ARKStep: nfi = {}", step_mem.nfi);
    let _ = writeln!(outfile, "ARKStep: nsetups = {}", step_mem.nsetups);
    let _ = writeln!(outfile, "ARKStep: nls_iters = {}", step_mem.nls_iters);
    let _ = writeln!(outfile, "ARKStep: nls_fails = {}", step_mem.nls_fails);

    ark_mem.step_mem = Some(step_mem);
}

/*---------------------------------------------------------------
  arkStep_AttachLinsol:

  This routine attaches the linear solver interface routines and
  solve type to the ARKStep module; the ARKLsMem box goes to
  ark_mem.lmem (Addendum C.2).
  ---------------------------------------------------------------*/
pub fn arkStep_AttachLinsol(
    ark_mem: &mut ARKodeMem,
    linit: Option<ARKLinsolInitFn>,
    lsetup: Option<ARKLinsolSetupFn>,
    lsolve: Option<ARKLinsolSolveFn>,
    lfree: Option<ARKLinsolFreeFn>,
    lsolve_type: crate::sundials_linearsolver::SUNLinearSolver_Type,
    lmem: Box<crate::arkode_ls_impl::ARKLsMem>,
) -> i32 {
    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_AttachLinsol") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* free any existing system solver */
    if let Some(old_lfree) = step_mem.lfree {
        old_lfree(ark_mem);
    }

    /* Attach the provided routines, data structure and solve type */
    step_mem.linit = linit;
    step_mem.lsetup = lsetup;
    step_mem.lsolve = lsolve;
    step_mem.lfree = lfree;
    step_mem.lsolve_type = lsolve_type;
    ark_mem.lmem = Some(lmem);

    /* Reset all linear solver counters */
    step_mem.nsetups = 0;
    step_mem.nstlp = 0;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_DisableLSetup:

  This routine NULLifies the lsetup function pointer in the
  ARKStep module.
  ---------------------------------------------------------------*/
pub fn arkStep_DisableLSetup(ark_mem: &mut ARKodeMem) {
    if let Some(b) = ark_mem.step_mem.as_mut() {
        if let Some(step_mem) = b.downcast_mut::<ARKodeARKStepMem>() {
            step_mem.lsetup = None;
        }
    }
}

/*---------------------------------------------------------------
  arkStep_AttachMasssol:

  This routine attaches the mass matrix linear solver interface
  routines and solver type to the ARKStep module; the
  ARKLsMassMem box goes to ark_mem.mass_mem (Addendum C.2).
  ---------------------------------------------------------------*/
#[allow(clippy::too_many_arguments)]
pub fn arkStep_AttachMasssol(
    ark_mem: &mut ARKodeMem,
    minit: Option<ARKMassInitFn>,
    msetup: Option<ARKMassSetupFn>,
    mmult: Option<ARKMassMultFn>,
    msolve: Option<ARKMassSolveFn>,
    mfree: Option<ARKMassFreeFn>,
    time_dep: bool,
    msolve_type: crate::sundials_linearsolver::SUNLinearSolver_Type,
    mass_mem: Box<crate::arkode_ls_impl::ARKLsMassMem>,
) -> i32 {
    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_AttachMasssol") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* free any existing mass matrix solver */
    if let Some(old_mfree) = step_mem.mfree {
        old_mfree(ark_mem);
    }

    /* Attach the provided routines, data structure and solve type */
    step_mem.minit = minit;
    step_mem.msetup = msetup;
    step_mem.mmult = mmult;
    step_mem.msolve = msolve;
    step_mem.mfree = mfree;
    step_mem.mass_type = if time_dep { MASS_TIMEDEP } else { MASS_FIXED };
    step_mem.msolve_type = msolve_type;
    ark_mem.mass_mem = Some(mass_mem);

    /* Attach mmult function pointer to ark_mem as well */
    ark_mem.step_mmult = mmult;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_DisableMSetup:

  This routine NULLifies the msetup function pointer in the
  ARKStep module.
  ---------------------------------------------------------------*/
pub fn arkStep_DisableMSetup(ark_mem: &mut ARKodeMem) {
    if let Some(b) = ark_mem.step_mem.as_mut() {
        if let Some(step_mem) = b.downcast_mut::<ARKodeARKStepMem>() {
            step_mem.msetup = None;
        }
    }
}

/*---------------------------------------------------------------
  arkStep_GetLmem:

  This routine returns the system linear solver interface memory
  (take semantics; put-back writes ark_mem.lmem — Addendum C.2).
  ---------------------------------------------------------------*/
pub fn arkStep_GetLmem(ark_mem: &mut ARKodeMem) -> Option<Box<crate::arkode_ls_impl::ARKLsMem>> {
    ark_mem.lmem.take()
}

/*---------------------------------------------------------------
  arkStep_GetMassMem:

  This routine returns the mass matrix solver interface memory
  (take semantics; put-back writes ark_mem.mass_mem).
  ---------------------------------------------------------------*/
pub fn arkStep_GetMassMem(
    ark_mem: &mut ARKodeMem,
) -> Option<Box<crate::arkode_ls_impl::ARKLsMassMem>> {
    ark_mem.mass_mem.take()
}

/*---------------------------------------------------------------
  arkStep_GetImplicitRHS:

  This routine returns the implicit RHS function pointer, fi.
  ---------------------------------------------------------------*/
pub fn arkStep_GetImplicitRHS(ark_mem: &mut ARKodeMem) -> Option<ARKRhsFn> {
    let step_mem = arkStep_AccessStepMem(ark_mem, "arkStep_GetImplicitRHS")?;
    let fi = step_mem.fi;
    ark_mem.step_mem = Some(step_mem);
    fi
}

/*---------------------------------------------------------------
  arkStep_GetGammas:

  This routine fills the current value of gamma, and states
  whether the gamma ratio fails the dgmax criteria.  (C hands out
  the ADDRESS of jcur; the Rust op copies it — write-back goes
  through arkStep_SetJcur.)
  ---------------------------------------------------------------*/
pub fn arkStep_GetGammas(
    ark_mem: &mut ARKodeMem,
    gamma: &mut f64,
    gamrat: &mut f64,
    jcur: &mut bool,
    dgamma_fail: &mut bool,
) -> i32 {
    let step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_GetGammas") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* set outputs */
    *gamma = step_mem.gamma;
    *gamrat = step_mem.gamrat;
    *jcur = step_mem.jcur;
    *dgamma_fail = SUNRabs(*gamrat - ONE) >= step_mem.dgmax;

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/// Rust-only companion of arkStep_GetGammas (Addendum C.2): carries
/// the ARKLS preconditioner-setup jcur write-back into the stepper.
pub fn arkStep_SetJcur(ark_mem: &mut ARKodeMem, jcur: bool) {
    if let Some(b) = ark_mem.step_mem.as_mut() {
        if let Some(step_mem) = b.downcast_mut::<ARKodeARKStepMem>() {
            step_mem.jcur = jcur;
        }
    }
}

/*---------------------------------------------------------------
  arkStep_Init:

  This routine is called just prior to performing internal time
  steps (after all user "set" routines have been called) from
  within arkInitialSetup (see the C source for the full
  init_type-dependent description).
  ---------------------------------------------------------------*/
pub fn arkStep_Init(ark_mem: &mut ARKodeMem, init_type: i32) -> i32 {
    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_Init") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };
    let ret = arkStep_Init_inner(ark_mem, &mut step_mem, init_type);
    ark_mem.step_mem = Some(step_mem);

    /* Initialize the nonlinear solver object (if it exists); done after
       putting step_mem back since arkStep_NlsInit re-accesses it */
    if ret == ARK_SUCCESS && init_type != RESET_INIT {
        let has_nls = {
            let sm = arkStep_AccessStepMem(ark_mem, "arkStep_Init").unwrap();
            let h = sm.NLS.is_some();
            ark_mem.step_mem = Some(sm);
            h
        };
        if has_nls {
            let retval = crate::arkode_arkstep_nls::arkStep_NlsInit(ark_mem);
            if retval != ARK_SUCCESS {
                arkProcessError(
                    Some(ark_mem),
                    ARK_NLS_INIT_FAIL,
                    line!(),
                    "arkStep_Init",
                    file!(),
                    "Unable to initialize SUNNonlinearSolver object",
                );
                return ARK_NLS_INIT_FAIL;
            }
        }

        /* Signal to shared arkode module that full RHS evaluations are required */
        ark_mem.call_fullrhs = true;
    }

    ret
}

fn arkStep_Init_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut ARKodeARKStepMem,
    init_type: i32,
) -> i32 {
    /* immediately return if reset */
    if init_type == RESET_INIT {
        return ARK_SUCCESS;
    }

    /* initializations/checks for (re-)initialization call */
    if init_type == FIRST_INIT {
        /* enforce use of arkEwtSmallReal if using a fixed step size for
           an explicit method, an internal error weight function, not
           using an iterative mass matrix solver with rwt=ewt, and not
           performing accumulated temporal error estimation */
        let mut reset_efun = true;
        if step_mem.implicit {
            reset_efun = false;
        }
        if !ark_mem.fixedstep {
            reset_efun = false;
        }
        if ark_mem.user_efun {
            reset_efun = false;
        }
        if ark_mem.AccumErrorType != ARK_ACCUMERROR_NONE {
            reset_efun = false;
        }
        /* (iterative mass-solver rwt checks: mass half not ported) */
        if reset_efun {
            ark_mem.user_efun = false;
            ark_mem.efun = Some(ark_ewt_small_real);
        }

        /* Create Butcher tables (if not already set) */
        let retval = arkStep_SetButcherTables(ark_mem, step_mem);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "arkStep_Init",
                file!(),
                "Could not create Butcher table(s)",
            );
            return ARK_ILL_INPUT;
        }

        /* Check that Butcher tables are OK */
        let retval = arkStep_CheckButcherTables(ark_mem, step_mem);
        if retval != ARK_SUCCESS {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "arkStep_Init",
                file!(),
                "Error in Butcher table(s)",
            );
            return ARK_ILL_INPUT;
        }

        /* Retrieve/store method and embedding orders now that tables are finalized */
        if let Some(bi) = &step_mem.Bi {
            step_mem.q = bi.q;
            step_mem.p = bi.p;
        } else {
            step_mem.q = step_mem.Be.as_ref().unwrap().q;
            step_mem.p = step_mem.Be.as_ref().unwrap().p;
        }
        if let Some(ha) = ark_mem.hadapt_mem.as_mut() {
            ha.q = step_mem.q;
            ha.p = step_mem.p;
        }

        /* Ensure that if adaptivity or error accumulation is enabled, then
           method includes embedding coefficients */
        if (!ark_mem.fixedstep || ark_mem.AccumErrorType != ARK_ACCUMERROR_NONE)
            && step_mem.p <= 0
        {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "arkStep_Init",
                file!(),
                "Temporal error estimation cannot be performed without embedding coefficients",
            );
            return ARK_ILL_INPUT;
        }

        /* Relaxation is incompatible with implicit RHS deduction */
        if ark_mem.relax_enabled && step_mem.implicit && step_mem.deduce_rhs {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "arkStep_Init",
                file!(),
                "Relaxation cannot be performed when deducing implicit RHS values",
            );
            return ARK_ILL_INPUT;
        }

        /* Allocate ARK RHS vector memory, update storage requirements */
        let tmpl_len = ark_mem.ewt.data.len();
        let (lrw1, liw1) = (ark_mem.lrw1, ark_mem.liw1);
        /*   Allocate Fe[0] ... Fe[stages-1] if needed */
        if step_mem.explicit {
            let ARKodeMem { lrw, liw, .. } = ark_mem;
            arkAllocVecArray(step_mem.stages, tmpl_len, &mut step_mem.Fe, lrw1, lrw, liw1, liw);
        }

        /*   Allocate Fi[0] ... Fi[stages-1] if needed */
        if step_mem.implicit {
            let ARKodeMem { lrw, liw, .. } = ark_mem;
            arkAllocVecArray(step_mem.stages, tmpl_len, &mut step_mem.Fi, lrw1, lrw, liw1, liw);
        }

        /* Allocate stage storage for relaxation with implicit/IMEX methods or
           if a fixed mass matrix is present */
        if ark_mem.relax_enabled
            && (step_mem.implicit || step_mem.mass_type == MASS_FIXED)
        {
            let ARKodeMem { lrw, liw, .. } = ark_mem;
            arkAllocVecArray(step_mem.stages, tmpl_len, &mut step_mem.z, lrw1, lrw, liw1, liw);
        }

        /* Allocate reusable arrays for fused vector operations */
        step_mem.nfusedopvecs = 2 * step_mem.stages + 2 + step_mem.nforcing;
        if step_mem.cvals.is_empty() {
            step_mem.cvals = vec![ZERO; step_mem.nfusedopvecs as usize];
            ark_mem.lrw += step_mem.nfusedopvecs as i64;
            /* (Xvecs assembled at call sites; keep the C liw accounting —
               C allocates Xvecs alongside cvals and only then adds liw,
               so a ReInit does not re-count it) */
            ark_mem.liw += step_mem.nfusedopvecs as i64;
        }

        /* Allocate workspace for MRI forcing */
        if step_mem.stage_times.is_empty() {
            step_mem.stage_times = vec![ZERO; step_mem.stages as usize];
            ark_mem.lrw += step_mem.stages as i64;
        }
        if step_mem.stage_coefs.is_empty() {
            step_mem.stage_coefs = vec![ZERO; step_mem.stages as usize];
            ark_mem.lrw += step_mem.stages as i64;
        }

        /* Override the interpolant degree (if needed), used in arkInitialSetup */
        if step_mem.q > 1 && ark_mem.interp_degree > (step_mem.q - 1) {
            /* Limit max degree to at most one less than the method global order */
            ark_mem.interp_degree = step_mem.q - 1;
        } else if step_mem.q == 1 && ark_mem.interp_degree > 1 {
            /* Allow for linear interpolant with first order methods to ensure
               solution values are returned at the time interval end points */
            ark_mem.interp_degree = 1;
        }

        /* Higher-order predictors require interpolation */
        if ark_mem.interp_type == ARK_INTERP_NONE && step_mem.predictor != 0 {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "arkStep_Init",
                file!(),
                "Non-trival predictors require an interpolation module",
            );
            return ARK_ILL_INPUT;
        }
    }

    /* set appropriate TakeStep routine based on problem configuration
       (adjoint TakeStep pending the adjoint machinery) */
    ark_mem.step = Some(arkStep_TakeStep_Z);

    /* Check for consistency between mass system and system linear system
       modules (e.g., if lsolve is direct, msolve needs to match) */
    if step_mem.mass_type != MASS_IDENTITY && ark_mem.lmem.is_some() {
        if step_mem.lsolve_type != step_mem.msolve_type {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "arkStep_Init",
                file!(),
                "Incompatible linear and mass matrix solvers",
            );
            return ARK_ILL_INPUT;
        }
    }

    /* Perform mass matrix solver initialization and setup (if applicable) */
    if step_mem.mass_type != MASS_IDENTITY {
        /* Call minit (if it exists) */
        if let Some(minit) = step_mem.minit {
            let retval = minit(ark_mem);
            if retval != 0 {
                arkProcessError(
                    Some(ark_mem),
                    ARK_MASSINIT_FAIL,
                    line!(),
                    "arkStep_Init",
                    file!(),
                    "The mass matrix solver's init routine failed.",
                );
                return ARK_MASSINIT_FAIL;
            }
        }

        /* Call msetup (if it exists) */
        if let Some(msetup) = step_mem.msetup {
            let tcur = ark_mem.tcur;
            let mut t1 = std::mem::take(&mut ark_mem.tempv1);
            let mut t2 = std::mem::take(&mut ark_mem.tempv2);
            let mut t3 = std::mem::take(&mut ark_mem.tempv3);
            let retval = msetup(ark_mem, tcur, &mut t1, &mut t2, &mut t3);
            ark_mem.tempv1 = t1;
            ark_mem.tempv2 = t2;
            ark_mem.tempv3 = t3;
            if retval != 0 {
                arkProcessError(
                    Some(ark_mem),
                    ARK_MASSSETUP_FAIL,
                    line!(),
                    "arkStep_Init",
                    file!(),
                    "The mass matrix solver's setup routine failed.",
                );
                return ARK_MASSSETUP_FAIL;
            }
        }
    }

    /* Call linit (if it exists) */
    if let Some(linit) = step_mem.linit {
        let retval = linit(ark_mem);
        if retval != 0 {
            arkProcessError(
                Some(ark_mem),
                ARK_LINIT_FAIL,
                line!(),
                "arkStep_Init",
                file!(),
                "The linear solver's init routine failed.",
            );
            return ARK_LINIT_FAIL;
        }
    }

    /* (NLS initialization and the call_fullrhs signal happen in the
       arkStep_Init wrapper, after step_mem is put back) */

    ARK_SUCCESS
}

/// C arkEwtSetSmallReal installed as an efun (see arkStep_Init).
fn ark_ewt_small_real(_ycur: &NVector, weight: &mut NVector, _e_data: &mut UserData) -> i32 {
    N_VConst(crate::sundials_types::SUN_SMALL_REAL, weight);
    ARK_SUCCESS
}

/*------------------------------------------------------------------------------
  arkStep_FullRHS:

  Rewriting the problem
    My' = fe(t,y) + fi(t,y)
  in the form
    y' = M^{-1}*[ fe(t,y) + fi(t,y) ],
  this routine computes the full right-hand side vector,
    f = M^{-1}*[ fe(t,y) + fi(t,y) ]

  (see the C source for the full mode description; the mass-matrix
  and MRI-forcing branches are unreachable until those halves land)
  ----------------------------------------------------------------------------*/
pub fn arkStep_FullRHS(ark_mem: &mut ARKodeMem, t: f64, y: &NVector, f: &mut NVector, mode: i32) -> i32 {
    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_FullRHS") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };
    let ret = arkStep_FullRHS_inner(ark_mem, &mut step_mem, t, y, f, mode);
    ark_mem.step_mem = Some(step_mem);
    ret
}

fn arkStep_FullRHS_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut ARKodeARKStepMem,
    t: f64,
    y: &NVector,
    f: &mut NVector,
    mode: i32,
) -> i32 {
    /* setup mass-matrix if required (use output f as a temporary) */
    if let (MASS_TIMEDEP, Some(msetup)) = (step_mem.mass_type, step_mem.msetup) {
        let mut t2 = std::mem::take(&mut ark_mem.tempv2);
        let mut t3 = std::mem::take(&mut ark_mem.tempv3);
        let retval = msetup(ark_mem, t, f, &mut t2, &mut t3);
        ark_mem.tempv2 = t2;
        ark_mem.tempv3 = t3;
        if retval != ARK_SUCCESS {
            return ARK_MASSSETUP_FAIL;
        }
    }

    /* perform RHS functions contingent on 'mode' argument */
    match mode {
        ARK_FULLRHS_START | ARK_FULLRHS_END => {
            /* compute the full RHS */
            if !ark_mem.fn_is_current {
                /* for ARK_FULLRHS_END, determine if RHS functions need to be
                   recomputed (they may be copied from the last stage) */
                let mut recompute_rhs = mode == ARK_FULLRHS_START;
                if mode == ARK_FULLRHS_END {
                    if step_mem.explicit
                        && !ARKodeButcherTable_IsStifflyAccurate(step_mem.Be.as_ref().unwrap())
                    {
                        recompute_rhs = true;
                    }
                    if step_mem.implicit
                        && !ARKodeButcherTable_IsStifflyAccurate(step_mem.Bi.as_ref().unwrap())
                    {
                        recompute_rhs = true;
                    }
                    /* Stiffly Accurate methods are not SA when relaxation is enabled */
                    if ark_mem.relax_enabled {
                        recompute_rhs = true;
                    }
                }

                if recompute_rhs {
                    /* call the user-supplied pre-RHS function (if supplied) */
                    if let Some(pre_rhs) = ark_mem.PreRhsFn {
                        let retval = pre_rhs(t, y, &mut ark_mem.user_data);
                        if retval != 0 {
                            arkProcessError(
                                Some(ark_mem),
                                ARK_PRERHSFN_FAIL,
                                line!(),
                                "arkStep_FullRHS",
                                file!(),
                                &format!(
                                    "At t = {}, the pre-RHS function failed in an unrecoverable manner.",
                                    t
                                ),
                            );
                            return ARK_PRERHSFN_FAIL;
                        }
                    }

                    /* compute the implicit component */
                    if step_mem.implicit {
                        let fi = step_mem.fi.unwrap();
                        let retval = fi(t, y, &mut step_mem.Fi[0], &mut ark_mem.user_data);
                        step_mem.nfi += 1;
                        if retval != 0 {
                            arkProcessError(
                                Some(ark_mem),
                                ARK_RHSFUNC_FAIL,
                                line!(),
                                "arkStep_FullRHS",
                                file!(),
                                &format!(
                                    "At t = {}, the right-hand side routine failed in an unrecoverable manner.",
                                    t
                                ),
                            );
                            return ARK_RHSFUNC_FAIL;
                        }

                        /* compute and store M(t)^{-1} fi */
                        if step_mem.mass_type == MASS_TIMEDEP {
                            let msolve = step_mem.msolve.unwrap();
                            let tol = step_mem.nlscoef / ark_mem.h;
                            let retval = msolve(ark_mem, &mut step_mem.Fi[0], tol);
                            if retval != 0 {
                                arkProcessError(
                                    Some(ark_mem),
                                    ARK_MASSSOLVE_FAIL,
                                    line!(),
                                    "arkStep_FullRHS",
                                    file!(),
                                    "Mass matrix solver failure",
                                );
                                return ARK_MASSSOLVE_FAIL;
                            }
                        }
                    }

                    /* compute the explicit component */
                    if step_mem.explicit {
                        let fe = step_mem.fe.unwrap();
                        let retval = fe(t, y, &mut step_mem.Fe[0], &mut ark_mem.user_data);
                        step_mem.nfe += 1;
                        if retval != 0 {
                            arkProcessError(
                                Some(ark_mem),
                                ARK_RHSFUNC_FAIL,
                                line!(),
                                "arkStep_FullRHS",
                                file!(),
                                &format!(
                                    "At t = {}, the right-hand side routine failed in an unrecoverable manner.",
                                    t
                                ),
                            );
                            return ARK_RHSFUNC_FAIL;
                        }

                        /* compute and store M(t)^{-1} fe */
                        if step_mem.mass_type == MASS_TIMEDEP {
                            let msolve = step_mem.msolve.unwrap();
                            let tol = step_mem.nlscoef / ark_mem.h;
                            let retval = msolve(ark_mem, &mut step_mem.Fe[0], tol);
                            if retval != 0 {
                                arkProcessError(
                                    Some(ark_mem),
                                    ARK_MASSSOLVE_FAIL,
                                    line!(),
                                    "arkStep_FullRHS",
                                    file!(),
                                    "Mass matrix solver failure",
                                );
                                return ARK_MASSSOLVE_FAIL;
                            }
                        }
                    }
                } else {
                    /* copy the RHS values from the last stage (FSAL) */
                    let s = step_mem.stages as usize;
                    if step_mem.explicit {
                        let (head, tail) = step_mem.Fe.split_at_mut(s - 1);
                        if s > 1 {
                            head[0].data.copy_from_slice(&tail[0].data);
                        }
                    }
                    if step_mem.implicit {
                        let (head, tail) = step_mem.Fi.split_at_mut(s - 1);
                        if s > 1 {
                            head[0].data.copy_from_slice(&tail[0].data);
                        }
                    }
                }
            }

            /* combine RHS vector(s) into output */
            if step_mem.explicit && step_mem.implicit {
                /* ImEx */
                N_VLinearSum(ONE, &step_mem.Fi[0], ONE, &step_mem.Fe[0], f);
            } else if step_mem.implicit {
                /* implicit */
                N_VScale(ONE, &step_mem.Fi[0], f);
            } else {
                /* explicit */
                N_VScale(ONE, &step_mem.Fe[0], f);
            }

            /* compute M^{-1} f for output but do not store */
            if step_mem.mass_type == MASS_FIXED {
                let msolve = step_mem.msolve.unwrap();
                let tol = step_mem.nlscoef / ark_mem.h;
                let retval = msolve(ark_mem, f, tol);
                if retval != 0 {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_MASSSOLVE_FAIL,
                        line!(),
                        "arkStep_FullRHS",
                        file!(),
                        "Mass matrix solver failure",
                    );
                    return ARK_MASSSOLVE_FAIL;
                }
            }

            /* apply external polynomial (MRI) forcing (M = I required) */
            if step_mem.expforcing || step_mem.impforcing {
                let vals = arkStep_ApplyForcing_coeffs(step_mem, &[t], &[ONE], 1);
                ark_accumulate_forcing(step_mem, &vals, f);
            }
        }

        ARK_FULLRHS_OTHER => {
            /* call the user-supplied pre-RHS function (if supplied) */
            if let Some(pre_rhs) = ark_mem.PreRhsFn {
                let retval = pre_rhs(t, y, &mut ark_mem.user_data);
                if retval != 0 {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_PRERHSFN_FAIL,
                        line!(),
                        "arkStep_FullRHS",
                        file!(),
                        &format!(
                            "At t = {}, the pre-RHS function failed in an unrecoverable manner.",
                            t
                        ),
                    );
                    return ARK_PRERHSFN_FAIL;
                }
            }

            /* compute the implicit component and store in sdata */
            if step_mem.implicit {
                let fi = step_mem.fi.unwrap();
                let retval = fi(t, y, &mut step_mem.sdata, &mut ark_mem.user_data);
                step_mem.nfi += 1;
                if retval != 0 {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_RHSFUNC_FAIL,
                        line!(),
                        "arkStep_FullRHS",
                        file!(),
                        &format!(
                            "At t = {}, the right-hand side routine failed in an unrecoverable manner.",
                            t
                        ),
                    );
                    return ARK_RHSFUNC_FAIL;
                }
            }

            /* compute the explicit component and store in ark_tempv2 */
            if step_mem.explicit {
                let fe = step_mem.fe.unwrap();
                let retval = fe(t, y, &mut ark_mem.tempv2, &mut ark_mem.user_data);
                step_mem.nfe += 1;
                if retval != 0 {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_RHSFUNC_FAIL,
                        line!(),
                        "arkStep_FullRHS",
                        file!(),
                        &format!(
                            "At t = {}, the right-hand side routine failed in an unrecoverable manner.",
                            t
                        ),
                    );
                    return ARK_RHSFUNC_FAIL;
                }
            }

            /* combine RHS vector(s) into output */
            if step_mem.explicit && step_mem.implicit {
                /* ImEx */
                N_VLinearSum(ONE, &step_mem.sdata, ONE, &ark_mem.tempv2, f);
            } else if step_mem.implicit {
                /* implicit */
                N_VScale(ONE, &step_mem.sdata, f);
            } else {
                /* explicit */
                N_VScale(ONE, &ark_mem.tempv2, f);
            }

            /* compute M^{-1} f for output but do not store */
            if step_mem.mass_type != MASS_IDENTITY {
                let msolve = step_mem.msolve.unwrap();
                let tol = step_mem.nlscoef / ark_mem.h;
                let retval = msolve(ark_mem, f, tol);
                if retval != 0 {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_MASSSOLVE_FAIL,
                        line!(),
                        "arkStep_FullRHS",
                        file!(),
                        "Mass matrix solver failure",
                    );
                    return ARK_MASSSOLVE_FAIL;
                }
            }

            /* apply external polynomial (MRI) forcing (M = I required) */
            if step_mem.expforcing || step_mem.impforcing {
                let vals = arkStep_ApplyForcing_coeffs(step_mem, &[t], &[ONE], 1);
                ark_accumulate_forcing(step_mem, &vals, f);
            }
        }

        _ => {
            /* return with RHS failure if unknown mode is passed */
            arkProcessError(
                Some(ark_mem),
                ARK_RHSFUNC_FAIL,
                line!(),
                "arkStep_FullRHS",
                file!(),
                "Unknown full RHS mode",
            );
            return ARK_RHSFUNC_FAIL;
        }
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_TakeStep_Z:

  This routine serves the primary purpose of the ARKStep module:
  it performs a single ARK step (with embedding, if possible).
  This version solves for each ARK stage vector, z_i.

  See the C source for the dsmPtr/nflagPtr conventions.
  ---------------------------------------------------------------*/
pub fn arkStep_TakeStep_Z(ark_mem: &mut ARKodeMem, dsmPtr: &mut f64, nflagPtr: &mut i32) -> i32 {
    /* access ARKodeARKStepMem structure; compute the step-level flags,
       then release step_mem for the full-RHS re-entry below */
    let (implicit_stage0, is_start, save_stages, imex_method, stiffly_accurate);
    let (save_fn_for_interp, save_fn_for_residual, eval_rhs);
    {
        let mut step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_TakeStep_Z") {
            None => return ARK_MEM_NULL,
            Some(sm) => sm,
        };

        /* if problem will involve no algebraic solvers, initialize nflagPtr
           to success */
        if !step_mem.implicit && step_mem.mass_type == MASS_IDENTITY {
            *nflagPtr = ARK_SUCCESS;
        }

        /* initialize the current stage index */
        step_mem.istage = 0;

        /* (nonlinear solver setup: the Newton and fixed-point modules have
           no setup routine, matching the C ops tables) */

        /* check if we need to store stage values */
        save_stages = ark_mem.relax_enabled
            && (step_mem.implicit || step_mem.mass_type == MASS_FIXED);

        /* check for an ImEx method */
        imex_method = step_mem.implicit && step_mem.explicit;

        /* check for implicit method with an explicit first stage */
        let mut ims = false;
        let mut isr = 1;
        if step_mem.implicit {
            if SUNRabs(step_mem.Bi.as_ref().unwrap().A[0][0]) > TINY {
                ims = true;
                isr = 0;
            }
        }
        implicit_stage0 = ims;
        is_start = isr;

        /* check if the method is Stiffly Accurate (SA) */
        let mut sa = true;
        if step_mem.explicit
            && !ARKodeButcherTable_IsStifflyAccurate(step_mem.Be.as_ref().unwrap())
        {
            sa = false;
        }
        if step_mem.implicit
            && !ARKodeButcherTable_IsStifflyAccurate(step_mem.Bi.as_ref().unwrap())
        {
            sa = false;
        }
        stiffly_accurate = sa;

        /* Save f(tn, yn) for Hermite interpolation */
        save_fn_for_interp = ark_mem.interp_type == ARK_INTERP_HERMITE;

        /* For an implicit or ImEx method using the trivial predictor with an
           autonomous problem with an identity or fixed mass matrix, save
           fi(tn, yn) for reuse in the first residual evaluation of each
           stage solve */
        save_fn_for_residual = step_mem.implicit
            && step_mem.predictor == 0
            && step_mem.autonomous
            && step_mem.mass_type != MASS_TIMEDEP;

        /* Call the RHS if needed. */
        eval_rhs = !implicit_stage0 || save_fn_for_interp || save_fn_for_residual;

        ark_mem.step_mem = Some(step_mem);
    }

    /* explicit first stage -- store stage if necessary for relaxation or
       checkpointing */
    if is_start == 1 {
        if save_stages {
            let mut step_mem = arkStep_AccessStepMem(ark_mem, "arkStep_TakeStep_Z").unwrap();
            let ARKodeMem { yn, .. } = ark_mem;
            step_mem.z[0].data.copy_from_slice(&yn.data);
            ark_mem.step_mem = Some(step_mem);
        }
        if ark_mem.checkpoint_scheme.is_some() {
            let retval = arkstep_checkpoint(ark_mem, 0, ark_mem.tn, CheckpointVec::Yn);
            if retval != ARK_SUCCESS {
                return retval;
            }
        }
    }

    if !ark_mem.fn_is_current && eval_rhs {
        /* If saving the RHS evaluation for reuse in the residual, call the
           full RHS for all implicit methods or for ImEx methods with an
           explicit first stage. */
        let res_full_rhs = save_fn_for_residual && implicit_stage0 && !imex_method;

        if !implicit_stage0 || save_fn_for_interp || res_full_rhs {
            /* Need full RHS evaluation (step_mem released around the
               step_fullrhs op re-entry) */
            let mode = if ark_mem.initsetup {
                ARK_FULLRHS_START
            } else {
                ARK_FULLRHS_END
            };
            let retval = ark_step_fullrhs_yn_fn(ark_mem, ark_mem.tn, mode);
            if retval != 0 {
                return ARK_RHSFUNC_FAIL;
            }
            ark_mem.fn_is_current = true;
        } else {
            /* For an ImEx method with implicit first stage and an
               interpolation method that does not need fn (e.g., Lagrange),
               only evaluate fi (if necessary) for reuse in the residual */
            let mut step_mem = arkStep_AccessStepMem(ark_mem, "arkStep_TakeStep_Z").unwrap();
            if stiffly_accurate {
                let s = step_mem.stages as usize;
                let (head, tail) = step_mem.Fi.split_at_mut(s - 1);
                if s > 1 {
                    head[0].data.copy_from_slice(&tail[0].data);
                }
            } else {
                /* call the user-supplied pre-RHS function (if supplied) */
                if let Some(pre_rhs) = ark_mem.PreRhsFn {
                    let retval = pre_rhs(ark_mem.tn, &ark_mem.yn, &mut ark_mem.user_data);
                    if retval != 0 {
                        ark_mem.step_mem = Some(step_mem);
                        return ARK_PRERHSFN_FAIL;
                    }
                }
                let fi = step_mem.fi.unwrap();
                let ARKodeMem { yn, user_data, tn, .. } = ark_mem;
                let retval = fi(*tn, yn, &mut step_mem.Fi[0], user_data);
                step_mem.nfi += 1;
                if retval < 0 {
                    ark_mem.step_mem = Some(step_mem);
                    return ARK_RHSFUNC_FAIL;
                }
                if retval > 0 {
                    ark_mem.step_mem = Some(step_mem);
                    return ARK_UNREC_RHSFUNC_ERR;
                }
            }
            ark_mem.step_mem = Some(step_mem);
        }
    }

    /* re-take step_mem for the stage loop and solution computation */
    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_TakeStep_Z") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    /* Set alias to implicit RHS evaluation for reuse in residual */
    step_mem.fn_implicit = FnImplicitAlias::None;
    if save_fn_for_residual {
        if !implicit_stage0 {
            /* Explicit first stage -- Fi[0] will be retained */
            step_mem.fn_implicit = FnImplicitAlias::Fi0;
        } else if imex_method || step_mem.mass_type == MASS_FIXED {
            /* Implicit first stage -- Fi[0] will be overwritten;
               copy from Fi[0] as fn includes fe or M^{-1} */
            ark_mem.tempv5.data.copy_from_slice(&step_mem.Fi[0].data);
            step_mem.fn_implicit = FnImplicitAlias::Tempv5;
        } else {
            /* fn is the same as Fi[0] but will not be overwritten */
            step_mem.fn_implicit = FnImplicitAlias::ArkFn;
        }
    }

    let ret = arkStep_TakeStep_Z_inner(
        ark_mem,
        &mut step_mem,
        dsmPtr,
        nflagPtr,
        is_start,
        save_stages,
        stiffly_accurate,
    );
    ark_mem.step_mem = Some(step_mem);
    ret
}

#[allow(clippy::too_many_arguments)]
fn arkStep_TakeStep_Z_inner(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut Box<ARKodeARKStepMem>,
    dsmPtr: &mut f64,
    nflagPtr: &mut i32,
    is_start: i32,
    save_stages: bool,
    stiffly_accurate: bool,
) -> i32 {
    /* loop over internal stages to the step */
    for is in is_start..step_mem.stages {
        /* store current stage index */
        step_mem.istage = is;

        /* determine whether implicit solve is required */
        let mut implicit_stage = false;
        if step_mem.implicit {
            if SUNRabs(step_mem.Bi.as_ref().unwrap().A[is as usize][is as usize]) > TINY {
                implicit_stage = true;
            }
        }

        /* determine if the stage RHS will be deduced from the implicit solve */
        let deduce_stage = step_mem.deduce_rhs && implicit_stage;

        /* set current stage time(s) */
        if step_mem.implicit {
            ark_mem.tcur = ark_mem.tn + step_mem.Bi.as_ref().unwrap().c[is as usize] * ark_mem.h;
        } else {
            ark_mem.tcur = ark_mem.tn + step_mem.Be.as_ref().unwrap().c[is as usize] * ark_mem.h;
        }

        /* setup time-dependent mass matrix */
        if let (MASS_TIMEDEP, Some(msetup)) = (step_mem.mass_type, step_mem.msetup) {
            let tcur = ark_mem.tcur;
            let mut t1 = std::mem::take(&mut ark_mem.tempv1);
            let mut t2 = std::mem::take(&mut ark_mem.tempv2);
            let mut t3 = std::mem::take(&mut ark_mem.tempv3);
            let retval = msetup(ark_mem, tcur, &mut t1, &mut t2, &mut t3);
            ark_mem.tempv1 = t1;
            ark_mem.tempv2 = t2;
            ark_mem.tempv3 = t3;
            if retval != ARK_SUCCESS {
                return ARK_MASSSETUP_FAIL;
            }
        }

        /* if implicit, call built-in and user-supplied predictors
           (results placed in zpred) */
        if implicit_stage {
            let retval = arkStep_Predict(ark_mem, step_mem, is);
            if retval != ARK_SUCCESS {
                return retval;
            }

            /* if a user-supplied predictor routine is provided, call that here.
               Note that arkStep_Predict is *still* called, so this user-supplied
               routine can just 'clean up' the built-in prediction, if desired. */
            if let Some(stage_predict) = step_mem.stage_predict {
                let retval =
                    stage_predict(ark_mem.tcur, &mut step_mem.zpred, &mut ark_mem.user_data);
                if retval < 0 {
                    return ARK_USER_PREDICT_FAIL;
                }
                if retval > 0 {
                    return TRY_AGAIN;
                }
            }
        }

        /* set up explicit data for evaluation of ARK stage (store in sdata) */
        let retval = arkStep_StageSetup(ark_mem, step_mem, implicit_stage);
        if retval != ARK_SUCCESS {
            return retval;
        }

        /* perform implicit solve if required */
        if implicit_stage {
            /* implicit solve result is stored in ark_mem->ycur;
               return with positive value on anything but success */
            *nflagPtr = crate::arkode_arkstep_nls::arkStep_Nls(ark_mem, step_mem, *nflagPtr);
            if *nflagPtr != ARK_SUCCESS {
                return TRY_AGAIN;
            }

            /* otherwise no implicit solve is needed */
        } else {
            /* if M is fixed, solve with it to compute update (place back
               in sdata) */
            if step_mem.mass_type == MASS_FIXED {
                /* perform solve; return with positive value on anything
                   but success */
                let msolve = step_mem.msolve.unwrap();
                let tol = step_mem.nlscoef;
                *nflagPtr = msolve(ark_mem, &mut step_mem.sdata, tol);
                if *nflagPtr != ARK_SUCCESS {
                    return TRY_AGAIN;
                }
            }

            /* set y to be yn + sdata (either computed in arkStep_StageSetup,
               or updated in prev. block) */
            let ARKodeMem { yn, ycur, .. } = ark_mem;
            N_VLinearSum(ONE, yn, ONE, &step_mem.sdata, ycur);
        }

        /* apply user-supplied stage postprocessing function (if supplied)
           unless this is the last stage of a FSAL method, then apply the
           user-supplied step postprocessing function instead (if supplied) */
        let post_step = ark_mem.PostProcessStepFn;
        let post_stage = ark_mem.PostProcessStageFn;
        let last_sa_stage = is == step_mem.stages - 1 && stiffly_accurate;
        if let (true, Some(post)) = (last_sa_stage, post_step) {
            let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
            let retval = post(*tcur, ycur, user_data);
            if retval != 0 {
                return ARK_POSTPROCESS_STEP_FAIL;
            }
        } else if let Some(post) = post_stage {
            let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
            let retval = post(*tcur, ycur, user_data);
            if retval != 0 {
                return ARK_POSTPROCESS_STAGE_FAIL;
            }
        }

        /* successful stage solve */

        /*    store stage (if necessary for relaxation) */
        if save_stages {
            let ARKodeMem { ycur, .. } = ark_mem;
            step_mem.z[is as usize].data.copy_from_slice(&ycur.data);
        }

        /*    checkpoint stage for adjoint (if necessary) */
        if ark_mem.checkpoint_scheme.is_some() {
            let retval =
                arkstep_checkpoint(ark_mem, is as i64, ark_mem.tcur, CheckpointVec::Ycur);
            if retval != ARK_SUCCESS {
                return retval;
            }
        }

        /* call the user-supplied pre-RHS function (if supplied) */
        if let Some(pre_rhs) = ark_mem.PreRhsFn {
            if (step_mem.implicit && !deduce_stage) || step_mem.explicit {
                let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
                let retval = pre_rhs(*tcur, ycur, user_data);
                if retval != 0 {
                    return ARK_PRERHSFN_FAIL;
                }
            }
        }

        /*    store implicit RHS (value in Fi[is] is from preceding nonlinear
              iteration) */
        if step_mem.implicit {
            if !deduce_stage {
                let fi = step_mem.fi.unwrap();
                let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
                let retval = fi(*tcur, ycur, &mut step_mem.Fi[is as usize], user_data);
                step_mem.nfi += 1;
                if retval < 0 {
                    return ARK_RHSFUNC_FAIL;
                }
                if retval > 0 {
                    return ARK_UNREC_RHSFUNC_ERR;
                }
            } else if step_mem.mass_type == MASS_FIXED {
                let mmult = step_mem.mmult.unwrap();
                let mut t1 = std::mem::take(&mut ark_mem.tempv1);
                let retval = mmult(ark_mem, &step_mem.zcor, &mut t1);
                if retval != ARK_SUCCESS {
                    ark_mem.tempv1 = t1;
                    return ARK_MASSMULT_FAIL;
                }
                let g = step_mem.gamma;
                let ARKodeARKStepMem { sdata, Fi, .. } = &mut **step_mem;
                N_VLinearSum(ONE / g, &t1, -ONE / g, sdata, &mut Fi[is as usize]);
                ark_mem.tempv1 = t1;
            } else {
                let g = step_mem.gamma;
                let ARKodeARKStepMem { zcor, sdata, Fi, .. } = &mut **step_mem;
                N_VLinearSum(ONE / g, zcor, -ONE / g, sdata, &mut Fi[is as usize]);
            }
        }

        /*    store explicit RHS */
        if step_mem.explicit {
            let fe = step_mem.fe.unwrap();
            let t_e = ark_mem.tn + step_mem.Be.as_ref().unwrap().c[is as usize] * ark_mem.h;
            let ARKodeMem { ycur, user_data, .. } = ark_mem;
            let retval = fe(t_e, ycur, &mut step_mem.Fe[is as usize], user_data);
            step_mem.nfe += 1;
            if retval < 0 {
                return ARK_RHSFUNC_FAIL;
            }
            if retval > 0 {
                return ARK_UNREC_RHSFUNC_ERR;
            }
        }

        /* if using a time-dependent mass matrix, update Fe[is] and/or
           Fi[is] with M(t)^{-1} */
        if step_mem.mass_type == MASS_TIMEDEP {
            /* If the implicit stage was deduced, it already includes
               M(t)^{-1} */
            if step_mem.implicit && !deduce_stage {
                let msolve = step_mem.msolve.unwrap();
                let tol = step_mem.nlscoef;
                *nflagPtr = msolve(ark_mem, &mut step_mem.Fi[is as usize], tol);
                if *nflagPtr != ARK_SUCCESS {
                    return TRY_AGAIN;
                }
            }
            if step_mem.explicit {
                let msolve = step_mem.msolve.unwrap();
                let tol = step_mem.nlscoef;
                *nflagPtr = msolve(ark_mem, &mut step_mem.Fe[is as usize], tol);
                if *nflagPtr != ARK_SUCCESS {
                    return TRY_AGAIN;
                }
            }
        }
    } /* loop over stages */

    /* compute time-evolved solution (in ark_ycur), error estimate (in dsm).
       This can fail recoverably due to nonconvergence of the mass matrix
       solve, so handle that appropriately. */
    ark_mem.tcur = ark_mem.tn + ark_mem.h;

    *nflagPtr = if step_mem.mass_type == MASS_FIXED {
        arkStep_ComputeSolutions_MassFixed(ark_mem, step_mem, dsmPtr)
    } else {
        arkStep_ComputeSolutions(ark_mem, step_mem, dsmPtr)
    };
    if *nflagPtr < 0 {
        return *nflagPtr;
    }
    if *nflagPtr > 0 {
        return TRY_AGAIN;
    }

    /* checkpoint the step solution (if necessary) */
    if ark_mem.checkpoint_scheme.is_some() {
        let retval = arkstep_checkpoint(
            ark_mem,
            step_mem.Be.as_ref().map(|b| b.stages).unwrap_or(step_mem.stages) as i64,
            ark_mem.tn + ark_mem.h,
            CheckpointVec::Ycur,
        );
        if retval != ARK_SUCCESS {
            return retval;
        }
    }

    ARK_SUCCESS
}

/* which state vector to checkpoint */
enum CheckpointVec {
    Yn,
    Ycur,
}

fn arkstep_checkpoint(ark_mem: &mut ARKodeMem, stage: i64, t: f64, which: CheckpointVec) -> i32 {
    use crate::sundials_adjointcheckpointscheme::{
        SUNAdjointCheckpointScheme_InsertVector, SUNAdjointCheckpointScheme_NeedsSaving,
    };

    let mut scheme = ark_mem.checkpoint_scheme.take().unwrap();

    let mut do_save = false;
    let errcode = SUNAdjointCheckpointScheme_NeedsSaving(
        &mut scheme,
        ark_mem.checkpoint_step_idx,
        stage,
        t,
        &mut do_save,
    );
    if errcode != 0 {
        ark_mem.checkpoint_scheme = Some(scheme);
        arkProcessError(
            Some(ark_mem),
            ARK_ADJ_CHECKPOINT_FAIL,
            line!(),
            "arkStep_TakeStep_Z",
            file!(),
            &format!("SUNAdjointCheckpointScheme_NeedsSaving returned {}", errcode),
        );
        return ARK_ADJ_CHECKPOINT_FAIL;
    }

    if do_save {
        let y = match which {
            CheckpointVec::Yn => &ark_mem.yn,
            CheckpointVec::Ycur => &ark_mem.ycur,
        };
        let errcode = SUNAdjointCheckpointScheme_InsertVector(
            &mut scheme,
            ark_mem.checkpoint_step_idx,
            stage,
            t,
            y,
        );
        if errcode != 0 {
            ark_mem.checkpoint_scheme = Some(scheme);
            arkProcessError(
                Some(ark_mem),
                ARK_ADJ_CHECKPOINT_FAIL,
                line!(),
                "arkStep_TakeStep_Z",
                file!(),
                &format!("SUNAdjointCheckpointScheme_InsertVector returned {}", errcode),
            );
            return ARK_ADJ_CHECKPOINT_FAIL;
        }
    }

    ark_mem.checkpoint_scheme = Some(scheme);
    ARK_SUCCESS
}

/*===============================================================
  Internal utility routines
  ===============================================================*/

/*---------------------------------------------------------------
  arkStep_AccessStepMem: takes the ARKStep memory out of ark_mem
  (callers must put it back).  C's arkStep_AccessARKODEStepMem
  (void* entry) collapses onto the same helper.
  ---------------------------------------------------------------*/
pub(crate) fn arkStep_AccessStepMem(
    ark_mem: &mut ARKodeMem,
    fname: &str,
) -> Option<Box<ARKodeARKStepMem>> {
    match ark_mem.step_mem.take() {
        None => {
            arkProcessError(
                Some(ark_mem),
                ARK_MEM_NULL,
                line!(),
                fname,
                file!(),
                MSG_ARKSTEP_NO_MEM,
            );
            None
        }
        Some(b) => match b.downcast::<ARKodeARKStepMem>() {
            Ok(sm) => Some(sm),
            Err(b) => {
                ark_mem.step_mem = Some(b);
                arkProcessError(
                    Some(ark_mem),
                    ARK_MEM_NULL,
                    line!(),
                    fname,
                    file!(),
                    MSG_ARKSTEP_NO_MEM,
                );
                None
            }
        },
    }
}

/*---------------------------------------------------------------
  arkStep_SetButcherTables

  This routine determines the ERK/DIRK/ARK method to use, based
  on the desired accuracy and information on whether the problem
  is explicit, implicit or imex.
  ---------------------------------------------------------------*/
pub(crate) fn arkStep_SetButcherTables(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut ARKodeARKStepMem,
) -> i32 {
    /* if tables have already been specified, just return */
    if step_mem.Be.is_some() || step_mem.Bi.is_some() {
        return ARK_SUCCESS;
    }

    /* initialize table numbers to illegal values */
    let mut etable = -1;
    let mut itable = -1;

    /**** ImEx methods ****/
    if step_mem.explicit && step_mem.implicit {
        match step_mem.q {
            2 => {
                etable = ARKSTEP_DEFAULT_ARK_ETABLE_2;
                itable = ARKSTEP_DEFAULT_ARK_ITABLE_2;
            }
            3 => {
                etable = ARKSTEP_DEFAULT_ARK_ETABLE_3;
                itable = ARKSTEP_DEFAULT_ARK_ITABLE_3;
            }
            4 => {
                etable = ARKSTEP_DEFAULT_ARK_ETABLE_4;
                itable = ARKSTEP_DEFAULT_ARK_ITABLE_4;
            }
            5 => {
                etable = ARKSTEP_DEFAULT_ARK_ETABLE_5;
                itable = ARKSTEP_DEFAULT_ARK_ITABLE_5;
            }
            _ => {
                /* no available method, set default */
                arkProcessError(
                    Some(ark_mem),
                    ARK_WARNING,
                    line!(),
                    "arkStep_SetButcherTables",
                    file!(),
                    "No ImEx method at requested order, using q=5.",
                );
                etable = ARKSTEP_DEFAULT_ARK_ETABLE_5;
                itable = ARKSTEP_DEFAULT_ARK_ITABLE_5;
            }
        }

        /**** implicit methods ****/
    } else if step_mem.implicit {
        match step_mem.q {
            1 => itable = ARKSTEP_DEFAULT_DIRK_1,
            2 => itable = ARKSTEP_DEFAULT_DIRK_2,
            3 => itable = ARKSTEP_DEFAULT_DIRK_3,
            4 => itable = ARKSTEP_DEFAULT_DIRK_4,
            5 => itable = ARKSTEP_DEFAULT_DIRK_5,
            _ => {
                /* no available method, set default */
                arkProcessError(
                    Some(ark_mem),
                    ARK_WARNING,
                    line!(),
                    "arkStep_SetButcherTables",
                    file!(),
                    "No implicit method at requested order, using q=5.",
                );
                itable = ARKSTEP_DEFAULT_DIRK_5;
            }
        }

        /**** explicit methods ****/
    } else {
        match step_mem.q {
            1 => etable = ARKSTEP_DEFAULT_ERK_1,
            2 => etable = ARKSTEP_DEFAULT_ERK_2,
            3 => etable = ARKSTEP_DEFAULT_ERK_3,
            4 => etable = ARKSTEP_DEFAULT_ERK_4,
            5 => etable = ARKSTEP_DEFAULT_ERK_5,
            6 => etable = ARKSTEP_DEFAULT_ERK_6,
            7 => etable = ARKSTEP_DEFAULT_ERK_7,
            8 => etable = ARKSTEP_DEFAULT_ERK_8,
            9 => etable = ARKSTEP_DEFAULT_ERK_9,
            _ => {
                /* no available method, set default */
                arkProcessError(
                    Some(ark_mem),
                    ARK_WARNING,
                    line!(),
                    "arkStep_SetButcherTables",
                    file!(),
                    "No explicit method at requested order, using q=9.",
                );
                etable = ARKSTEP_DEFAULT_ERK_9;
            }
        }
    }

    if etable > -1 {
        step_mem.Be = ARKodeButcherTable_LoadERK(etable);
    }
    if itable > -1 {
        step_mem.Bi = ARKodeButcherTable_LoadDIRK(itable);
    }

    /* note Butcher table space requirements */
    let (mut bliw, mut blrw) = (0i64, 0i64);
    if let Some(be) = &step_mem.Be {
        ARKodeButcherTable_Space(be, &mut bliw, &mut blrw);
        ark_mem.liw += bliw;
        ark_mem.lrw += blrw;
    }
    if let Some(bi) = &step_mem.Bi {
        ARKodeButcherTable_Space(bi, &mut bliw, &mut blrw);
        ark_mem.liw += bliw;
        ark_mem.lrw += blrw;
    }

    /* set [redundant] ARK stored values for stage numbers and method orders */
    if let Some(be) = &step_mem.Be {
        step_mem.stages = be.stages;
        step_mem.q = be.q;
        step_mem.p = be.p;
    }
    if let Some(bi) = &step_mem.Bi {
        step_mem.stages = bi.stages;
        step_mem.q = bi.q;
        step_mem.p = bi.p;
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_CheckButcherTables

  This routine runs through the explicit and/or implicit Butcher
  tables to ensure that they meet all necessary requirements
  (see the C source for the list).  Returns ARK_SUCCESS if tables
  pass, ARK_INVALID_TABLE otherwise.
  ---------------------------------------------------------------*/
pub(crate) fn arkStep_CheckButcherTables(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut ARKodeARKStepMem,
) -> i32 {
    let tol: f64 = 100.0 * SUN_UNIT_ROUNDOFF;

    /* check that the expected tables are set */
    if step_mem.explicit && step_mem.Be.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_INVALID_TABLE,
            line!(),
            "arkStep_CheckButcherTables",
            file!(),
            "explicit table is NULL!",
        );
        return ARK_INVALID_TABLE;
    }

    if step_mem.implicit && step_mem.Bi.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_INVALID_TABLE,
            line!(),
            "arkStep_CheckButcherTables",
            file!(),
            "implicit table is NULL!",
        );
        return ARK_INVALID_TABLE;
    }

    /* check that stages > 0 */
    if step_mem.stages < 1 {
        arkProcessError(
            Some(ark_mem),
            ARK_INVALID_TABLE,
            line!(),
            "arkStep_CheckButcherTables",
            file!(),
            "stages < 1!",
        );
        return ARK_INVALID_TABLE;
    }

    /* check that method order q > 0 */
    if step_mem.q < 1 {
        arkProcessError(
            Some(ark_mem),
            ARK_INVALID_TABLE,
            line!(),
            "arkStep_CheckButcherTables",
            file!(),
            "method order < 1!",
        );
        return ARK_INVALID_TABLE;
    }

    /* check that embedding order p > 0 */
    if step_mem.p < 1 && !ark_mem.fixedstep {
        arkProcessError(
            Some(ark_mem),
            ARK_INVALID_TABLE,
            line!(),
            "arkStep_CheckButcherTables",
            file!(),
            "embedding order < 1!",
        );
        return ARK_INVALID_TABLE;
    }

    /* check that embedding exists */
    if step_mem.p > 0 && !ark_mem.fixedstep {
        if step_mem.implicit && step_mem.Bi.as_ref().unwrap().d.is_none() {
            arkProcessError(
                Some(ark_mem),
                ARK_INVALID_TABLE,
                line!(),
                "arkStep_CheckButcherTables",
                file!(),
                "no implicit embedding!",
            );
            return ARK_INVALID_TABLE;
        }
        if step_mem.explicit && step_mem.Be.as_ref().unwrap().d.is_none() {
            arkProcessError(
                Some(ark_mem),
                ARK_INVALID_TABLE,
                line!(),
                "arkStep_CheckButcherTables",
                file!(),
                "no explicit embedding!",
            );
            return ARK_INVALID_TABLE;
        }
    }

    /* check that ERK table is strictly lower triangular */
    if step_mem.explicit {
        let mut okay = true;
        let be = step_mem.Be.as_ref().unwrap();
        for i in 0..step_mem.stages as usize {
            for j in i..step_mem.stages as usize {
                if SUNRabs(be.A[i][j]) > tol {
                    okay = false;
                }
            }
        }
        if !okay {
            arkProcessError(
                Some(ark_mem),
                ARK_INVALID_TABLE,
                line!(),
                "arkStep_CheckButcherTables",
                file!(),
                "Ae Butcher table is implicit!",
            );
            return ARK_INVALID_TABLE;
        }
    }

    /* check that IRK table is implicit and lower triangular */
    if step_mem.implicit {
        let bi = step_mem.Bi.as_ref().unwrap();
        let mut okay = false;
        for i in 0..step_mem.stages as usize {
            if SUNRabs(bi.A[i][i]) > tol {
                okay = true;
            }
        }
        if !okay {
            arkProcessError(
                Some(ark_mem),
                ARK_INVALID_TABLE,
                line!(),
                "arkStep_CheckButcherTables",
                file!(),
                "Ai Butcher table is explicit!",
            );
            return ARK_INVALID_TABLE;
        }

        let mut okay = true;
        for i in 0..step_mem.stages as usize {
            for j in (i + 1)..step_mem.stages as usize {
                if SUNRabs(bi.A[i][j]) > tol {
                    okay = false;
                }
            }
        }
        if !okay {
            arkProcessError(
                Some(ark_mem),
                ARK_INVALID_TABLE,
                line!(),
                "arkStep_CheckButcherTables",
                file!(),
                "Ai Butcher table has entries above diagonal!",
            );
            return ARK_INVALID_TABLE;
        }
    }

    /* Check if the method is compatible with relaxation */
    if ark_mem.relax_enabled {
        if step_mem.q < 2 {
            arkProcessError(
                Some(ark_mem),
                ARK_INVALID_TABLE,
                line!(),
                "arkStep_CheckButcherTables",
                file!(),
                "The Butcher table(s) must be at least second order!",
            );
            return ARK_INVALID_TABLE;
        }

        if step_mem.explicit {
            /* Check if all b values are positive */
            let be = step_mem.Be.as_ref().unwrap();
            for i in 0..step_mem.stages as usize {
                if be.b[i] < ZERO {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_INVALID_TABLE,
                        line!(),
                        "arkStep_CheckButcherTables",
                        file!(),
                        "The explicit Butcher table has a negative b value!",
                    );
                    return ARK_INVALID_TABLE;
                }
            }
        }

        if step_mem.implicit {
            /* Check if all b values are positive */
            let bi = step_mem.Bi.as_ref().unwrap();
            for i in 0..step_mem.stages as usize {
                if bi.b[i] < ZERO {
                    arkProcessError(
                        Some(ark_mem),
                        ARK_INVALID_TABLE,
                        line!(),
                        "arkStep_CheckButcherTables",
                        file!(),
                        "The implicit Butcher table has a negative b value!",
                    );
                    return ARK_INVALID_TABLE;
                }
            }
        }
    }

    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_Predict

  This routine computes the prediction for a specific internal
  stage solution, storing the result in step_mem->zpred.  The
  prediction is done using the interpolation structure in
  extrapolation mode, hence stages "far" from the previous time
  interval are predicted using lower order polynomials than the
  "nearby" stages.
  ---------------------------------------------------------------*/
pub(crate) fn arkStep_Predict(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut ARKodeARKStepMem,
    istage: i32,
) -> i32 {
    /* verify that interpolation structure is provided */
    if ark_mem.interp.is_none() && step_mem.predictor > 0 && step_mem.predictor < 4 {
        arkProcessError(
            Some(ark_mem),
            ARK_MEM_NULL,
            line!(),
            "arkStep_Predict",
            file!(),
            "Interpolation structure is NULL",
        );
        return ARK_MEM_NULL;
    }

    /* if the first step, use initial condition as guess */
    if ark_mem.initsetup {
        step_mem.zpred.data.copy_from_slice(&ark_mem.yn.data);
        return ARK_SUCCESS;
    }

    /* set evaluation time tau as relative shift from previous successful time */
    let tau =
        step_mem.Bi.as_ref().unwrap().c[istage as usize] * ark_mem.h / ark_mem.hold;

    /* use requested predictor formula */
    match step_mem.predictor {
        1 => {
            /***** Interpolatory Predictor 1 -- all to max order *****/
            let retval =
                crate::arkode::arkPredict_MaximumOrder(ark_mem, tau, &mut step_mem.zpred);
            if retval != ARK_ILL_INPUT {
                return retval;
            }
        }
        2 => {
            /***** Interpolatory Predictor 2 -- decrease order w/ increasing
                   level of extrapolation *****/
            let retval =
                crate::arkode::arkPredict_VariableOrder(ark_mem, tau, &mut step_mem.zpred);
            if retval != ARK_ILL_INPUT {
                return retval;
            }
        }
        3 => {
            /***** Cutoff predictor: max order interpolatory output for stages
                   "close" to previous step, first-order predictor for
                   subsequent stages *****/
            let retval =
                crate::arkode::arkPredict_CutoffOrder(ark_mem, tau, &mut step_mem.zpred);
            if retval != ARK_ILL_INPUT {
                return retval;
            }
        }
        4 => {
            /***** Bootstrap predictor: if any previous stage in step has
                   nonzero c_i, construct a quadratic Hermite interpolant for
                   prediction; otherwise use the trivial predictor. *****/
            let bi_c: Vec<f64> = step_mem.Bi.as_ref().unwrap().c.clone();

            /* determine if any previous stages in step meet criteria */
            let mut jstage: i32 = -1;
            for i in 0..istage {
                jstage = if bi_c[i as usize] != ZERO { i } else { jstage };
            }

            /* if using the trivial predictor, fall through */
            if jstage != -1 {
                /* find the "optimal" previous stage to use */
                for i in 0..istage {
                    if bi_c[i as usize] > bi_c[jstage as usize] && bi_c[i as usize] != ZERO
                    {
                        jstage = i;
                    }
                }

                /* set stage time, stage RHS and interpolation values */
                let h = ark_mem.h * bi_c[jstage as usize];
                let tau_b = ark_mem.h * bi_c[istage as usize];
                let mut cv: Vec<f64> = Vec::new();
                let mut xr: Vec<&NVector> = Vec::new();
                if step_mem.implicit {
                    /* Implicit piece */
                    cv.push(ONE);
                    xr.push(&step_mem.Fi[jstage as usize]);
                }
                if step_mem.explicit {
                    /* Explicit piece */
                    cv.push(ONE);
                    xr.push(&step_mem.Fe[jstage as usize]);
                }

                /* call predictor routine */
                let nvec = cv.len();
                let retval = crate::arkode::arkPredict_Bootstrap(
                    ark_mem,
                    h,
                    tau_b,
                    nvec,
                    &cv,
                    &xr,
                    &mut step_mem.zpred,
                );
                if retval != ARK_ILL_INPUT {
                    return retval;
                }
            }
        }
        5 => {
            /***** Minimal correction predictor: use all previous stage
                   information in this step *****/
            let mut cv: Vec<f64> = Vec::new();
            let mut xr: Vec<&NVector> = Vec::new();
            if step_mem.explicit {
                /* Explicit pieces */
                let be = step_mem.Be.as_ref().unwrap();
                for jstage in 0..istage as usize {
                    cv.push(ark_mem.h * be.A[istage as usize][jstage]);
                    xr.push(&step_mem.Fe[jstage]);
                }
            }
            if step_mem.implicit {
                /* Implicit pieces */
                let bi = step_mem.Bi.as_ref().unwrap();
                for jstage in 0..istage as usize {
                    cv.push(ark_mem.h * bi.A[istage as usize][jstage]);
                    xr.push(&step_mem.Fi[jstage]);
                }
            }
            cv.push(ONE);
            xr.push(&ark_mem.yn);

            /* compute predictor */
            N_VLinearCombination(cv.len() as i32, &cv, &xr, &mut step_mem.zpred);
            return ARK_SUCCESS;
        }
        _ => {}
    }

    /* if we made it here, use the trivial predictor (previous step solution) */
    step_mem.zpred.data.copy_from_slice(&ark_mem.yn.data);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_StageSetup

  This routine sets up the stage data for computing the RK
  residual, along with the step- and method-related factors
  gamma, gammap and gamrat.  (See the C source for the full
  explicit/implicit x mass-matrix mode description; only the
  identity-mass modes are reachable until the mass half lands:

  Explicit:            sdata = h*sum_{j<i} (Ae(i,j)*Fe(j) + Ai(i,j)*Fi(j))
  Implicit (M == I):   sdata = yn - zp + h*sum_{j<i} (...)
  ---------------------------------------------------------------*/
pub(crate) fn arkStep_StageSetup(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut ARKodeARKStepMem,
    implicit: bool,
) -> i32 {
    /* Set shortcut to current stage index */
    let i = step_mem.istage as usize;

    /* Update gamma if stage is implicit */
    if implicit {
        step_mem.gamma = ark_mem.h * step_mem.Bi.as_ref().unwrap().A[i][i];
        if ark_mem.firststage {
            step_mem.gammap = step_mem.gamma;
        }
        step_mem.gamrat = if ark_mem.firststage {
            ONE
        } else {
            step_mem.gamma / step_mem.gammap /* protect x/x != 1.0 */
        };
    }

    /* If implicit, initialize sdata to yn - zpred (here: zpred = zp), and
       set first entries for eventual N_VLinearCombination call */
    if implicit {
        let ARKodeARKStepMem { zpred, sdata, .. } = step_mem;
        N_VLinearSum(ONE, &ark_mem.yn, -ONE, zpred, sdata);
    }

    /* If implicit with fixed M!=I, update sdata with M*sdata */
    if implicit && step_mem.mass_type == MASS_FIXED {
        let mmult = step_mem.mmult.unwrap();
        let mut t1 = std::mem::take(&mut ark_mem.tempv1);
        t1.data.copy_from_slice(&step_mem.sdata.data);
        let retval = mmult(ark_mem, &t1, &mut step_mem.sdata);
        ark_mem.tempv1 = t1;
        if retval != ARK_SUCCESS {
            return ARK_MASSMULT_FAIL;
        }
    }

    /* Update sdata with prior stage information: assemble the fused-op
       operand list (Xvecs at call site; in-place accumulate variant when
       sdata leads the list) */
    let mut cv: Vec<f64> = Vec::new();
    {
        if implicit {
            cv.push(ONE);
        }
        if step_mem.explicit {
            let be = step_mem.Be.as_ref().unwrap();
            for j in 0..i {
                cv.push(ark_mem.h * be.A[i][j]);
            }
        }
        if step_mem.implicit {
            let bi = step_mem.Bi.as_ref().unwrap();
            for j in 0..i {
                cv.push(ark_mem.h * bi.A[i][j]);
            }
        }
    }

    /* apply external polynomial (MRI) forcing (M = I required) */
    if step_mem.expforcing || step_mem.impforcing {
        let (jmax, is_exp) = if step_mem.expforcing { (i, true) } else { (i + 1, false) };
        let mut stage_times: Vec<f64> = Vec::with_capacity(jmax);
        let mut stage_coefs: Vec<f64> = Vec::with_capacity(jmax);
        {
            let b = if is_exp {
                step_mem.Be.as_ref().unwrap()
            } else {
                step_mem.Bi.as_ref().unwrap()
            };
            for j in 0..jmax {
                stage_times.push(ark_mem.tn + b.c[j] * ark_mem.h);
                stage_coefs.push(ark_mem.h * b.A[i][j]);
            }
        }
        let vals = arkStep_ApplyForcing_coeffs(step_mem, &stage_times, &stage_coefs, jmax);
        cv.extend_from_slice(&vals);
    }

    /* call fused vector operation to do the work */
    if implicit {
        /* z == X[0] with c0 == 1: in-place accumulate form */
        let ARKodeARKStepMem { sdata, Fe, Fi, forcing, .. } = step_mem;
        let mut xr: Vec<&NVector> = Vec::new();
        if step_mem.explicit {
            for j in 0..i {
                xr.push(&Fe[j]);
            }
        }
        if step_mem.implicit {
            for j in 0..i {
                xr.push(&Fi[j]);
            }
        }
        if step_mem.expforcing || step_mem.impforcing {
            for v in forcing.iter() {
                xr.push(v);
            }
        }
        ark_lincomb_accumulate(&cv[1..], &xr, sdata);
    } else {
        let ARKodeARKStepMem { sdata, Fe, Fi, forcing, .. } = step_mem;
        let mut xr: Vec<&NVector> = Vec::new();
        if step_mem.explicit {
            for j in 0..i {
                xr.push(&Fe[j]);
            }
        }
        if step_mem.implicit {
            for j in 0..i {
                xr.push(&Fi[j]);
            }
        }
        if step_mem.expforcing || step_mem.impforcing {
            for v in forcing.iter() {
                xr.push(v);
            }
        }
        N_VLinearCombination(cv.len() as i32, &cv, &xr, sdata);
    }

    /* return with success */
    ARK_SUCCESS
}

/// z += sum_k c[k] * x[k] — the z == X[0], c[0] == 1 branch of the
/// C N_VLinearCombination kernel (outer loop over vectors, matching
/// the serial kernel's accumulation order bit-for-bit).
fn ark_lincomb_accumulate(cvals: &[f64], xvecs: &[&NVector], z: &mut NVector) {
    for (k, val) in cvals.iter().enumerate() {
        for e in 0..z.data.len() {
            z.data[e] += val * xvecs[k].data[e];
        }
    }
}

/*------------------------------------------------------------------------------
  arkStep_ApplyForcing

  Determines the scaling values necessary for the MRI polynomial
  forcing terms.  C appends the values and N_Vector pointers to the
  cvals/Xvecs arrays; the Rust port returns the scaling values (the
  forcing vectors are appended to the operand list at the call site).
  ----------------------------------------------------------------------------*/
fn arkStep_ApplyForcing_coeffs(
    step_mem: &ARKodeARKStepMem,
    stage_times: &[f64],
    stage_coefs: &[f64],
    jmax: usize,
) -> Vec<f64> {
    let nforcing = step_mem.nforcing as usize;
    let mut vals = vec![ZERO; nforcing];

    for j in 0..jmax {
        let tau = (stage_times[j] - step_mem.tshift) / step_mem.tscale;
        let mut taui = ONE;

        for k in 0..nforcing {
            vals[k] += stage_coefs[j] * taui;
            taui *= tau;
        }
    }

    vals
}

/// z += sum_k vals[k] * forcing[k] — the z == X[0], c[0] == 1 branch
/// of the C N_VLinearCombination kernel used when applying forcing
/// to a full-RHS output.
fn ark_accumulate_forcing(step_mem: &ARKodeARKStepMem, vals: &[f64], f: &mut NVector) {
    for (k, val) in vals.iter().enumerate() {
        for e in 0..f.data.len() {
            f.data[e] += val * step_mem.forcing[k].data[e];
        }
    }
}

/*------------------------------------------------------------------------------
  arkStep_SetInnerForcing

  Sets an array of coefficient vectors for a time-dependent external polynomial
  forcing term in the ODE RHS i.e., y' = fe(t,y) + fi(t,y) + p(t). This
  function is primarily intended for use with multirate integration methods
  (e.g., MRIStep) where ARKStep is used to solve a modified ODE at a fast time
  scale. The polynomial is of the form

  p(t) = sum_{i = 0}^{nvecs - 1} forcing[i] * ((t - tshift) / (tscale))^i

  where tshift and tscale are used to normalize the time t (e.g., with MRIGARK
  methods).  The C code stores the caller's vector-array pointer, the
  Rust port stores owned copies.
  ----------------------------------------------------------------------------*/
pub fn arkStep_SetInnerForcing(
    ark_mem: &mut ARKodeMem,
    tshift: f64,
    tscale: f64,
    forcing: &[NVector],
    nvecs: i32,
) -> i32 {
    /* access ARKodeARKStepMem structure */
    let mut step_mem = match arkStep_AccessStepMem(ark_mem, "arkStep_SetInnerForcing") {
        None => return ARK_MEM_NULL,
        Some(sm) => sm,
    };

    if nvecs > 0 {
        /* enable forcing */
        if step_mem.explicit {
            step_mem.expforcing = true;
            step_mem.impforcing = false;
        } else {
            step_mem.expforcing = false;
            step_mem.impforcing = true;
        }
        step_mem.tshift = tshift;
        step_mem.tscale = tscale;
        step_mem.forcing = forcing.to_vec();
        step_mem.nforcing = nvecs;

        /* If cvals and Xvecs are not allocated then arkStep_Init has not been
           called and the number of stages has not been set yet. These arrays will
           be allocated in arkStep_Init and take into account the value of nforcing.
           On subsequent calls will check if enough space has allocated in case
           nforcing has increased since the original allocation. */
        if !step_mem.cvals.is_empty()
            && (step_mem.nfusedopvecs - nvecs) < (2 * step_mem.stages + 2)
        {
            /* free current work space */
            ark_mem.lrw -= step_mem.nfusedopvecs as i64;
            ark_mem.liw -= step_mem.nfusedopvecs as i64;

            /* allocate reusable arrays for fused vector operations */
            step_mem.nfusedopvecs = 2 * step_mem.stages + 2 + nvecs;
            step_mem.cvals = vec![ZERO; step_mem.nfusedopvecs as usize];
            ark_mem.lrw += step_mem.nfusedopvecs as i64;
            ark_mem.liw += step_mem.nfusedopvecs as i64;
        }
    } else {
        /* disable forcing */
        step_mem.expforcing = false;
        step_mem.impforcing = false;
        step_mem.tshift = ZERO;
        step_mem.tscale = ONE;
        step_mem.forcing = Vec::new();
        step_mem.nforcing = 0;
    }

    ark_mem.step_mem = Some(step_mem);
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkStep_ComputeSolutions

  This routine calculates the final RK solution using the existing
  data.  This solution is placed directly in ark_ycur.  This
  routine also computes the error estimate ||y-ytilde||_WRMS,
  where ytilde is the embedded solution, and the norm weights come
  from ark_ewt.  This norm value is returned.  The vector form of
  this estimated error (y-ytilde) is stored in ark_mem->tempv1, in
  case the calling routine wishes to examine the error locations.

  This version assumes either an identity or time-dependent mass
  matrix (identical steps).
  ---------------------------------------------------------------*/
pub(crate) fn arkStep_ComputeSolutions(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut ARKodeARKStepMem,
    dsmPtr: &mut f64,
) -> i32 {
    /* initialize output */
    *dsmPtr = ZERO;

    /* check if the method is stiffly accurate */
    let mut stiffly_accurate = true;
    if step_mem.explicit
        && !ARKodeButcherTable_IsStifflyAccurate(step_mem.Be.as_ref().unwrap())
    {
        stiffly_accurate = false;
    }
    if step_mem.implicit
        && !ARKodeButcherTable_IsStifflyAccurate(step_mem.Bi.as_ref().unwrap())
    {
        stiffly_accurate = false;
    }

    /* If the method is stiffly accurate, ycur is already the new solution */

    if !stiffly_accurate {
        /* Compute time step solution (if necessary) */
        /*   set arrays for fused vector operation */
        let mut cv: Vec<f64> = vec![ONE];
        for j in 0..step_mem.stages as usize {
            if step_mem.explicit {
                /* Explicit pieces */
                cv.push(ark_mem.h * step_mem.Be.as_ref().unwrap().b[j]);
            }
            if step_mem.implicit {
                /* Implicit pieces */
                cv.push(ark_mem.h * step_mem.Bi.as_ref().unwrap().b[j]);
            }
        }

        /* apply external polynomial (MRI) forcing (M = I required) */
        if step_mem.expforcing || step_mem.impforcing {
            let stages = step_mem.stages as usize;
            let mut stage_times: Vec<f64> = Vec::with_capacity(stages);
            let mut stage_coefs: Vec<f64> = Vec::with_capacity(stages);
            {
                let b = if step_mem.expforcing {
                    step_mem.Be.as_ref().unwrap()
                } else {
                    step_mem.Bi.as_ref().unwrap()
                };
                for j in 0..stages {
                    stage_times.push(ark_mem.tn + b.c[j] * ark_mem.h);
                    stage_coefs.push(ark_mem.h * b.b[j]);
                }
            }
            let vals = arkStep_ApplyForcing_coeffs(step_mem, &stage_times, &stage_coefs, stages);
            cv.extend_from_slice(&vals);
        }

        /*   call fused vector operation to do the work */
        {
            let ARKodeMem { yn, ycur, .. } = ark_mem;
            let mut xr: Vec<&NVector> = vec![yn];
            for j in 0..step_mem.stages as usize {
                if step_mem.explicit {
                    xr.push(&step_mem.Fe[j]);
                }
                if step_mem.implicit {
                    xr.push(&step_mem.Fi[j]);
                }
            }
            if step_mem.expforcing || step_mem.impforcing {
                for v in step_mem.forcing.iter() {
                    xr.push(v);
                }
            }
            N_VLinearCombination(cv.len() as i32, &cv, &xr, ycur);
        }

        if let Some(post) = ark_mem.PostProcessStepFn {
            let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
            let retval = post(*tcur, ycur, user_data);
            if retval != 0 {
                return ARK_POSTPROCESS_STEP_FAIL;
            }
        }
    }

    /* Compute yerr (if temporal error estimation is enabled). */
    if !ark_mem.fixedstep || ark_mem.AccumErrorType != ARK_ACCUMERROR_NONE {
        /* set arrays for fused vector operation */
        let mut cv: Vec<f64> = Vec::new();
        for j in 0..step_mem.stages as usize {
            if step_mem.explicit {
                /* Explicit pieces */
                let be = step_mem.Be.as_ref().unwrap();
                let d = be.d.as_ref().unwrap();
                cv.push(ark_mem.h * (be.b[j] - d[j]));
            }
            if step_mem.implicit {
                /* Implicit pieces */
                let bi = step_mem.Bi.as_ref().unwrap();
                let d = bi.d.as_ref().unwrap();
                cv.push(ark_mem.h * (bi.b[j] - d[j]));
            }
        }

        /* apply external polynomial (MRI) forcing (M = I required) */
        if step_mem.expforcing || step_mem.impforcing {
            let stages = step_mem.stages as usize;
            let mut stage_times: Vec<f64> = Vec::with_capacity(stages);
            let mut stage_coefs: Vec<f64> = Vec::with_capacity(stages);
            {
                let b = if step_mem.expforcing {
                    step_mem.Be.as_ref().unwrap()
                } else {
                    step_mem.Bi.as_ref().unwrap()
                };
                let d = b.d.as_ref().unwrap();
                for j in 0..stages {
                    stage_times.push(ark_mem.tn + b.c[j] * ark_mem.h);
                    stage_coefs.push(ark_mem.h * (b.b[j] - d[j]));
                }
            }
            let vals = arkStep_ApplyForcing_coeffs(step_mem, &stage_times, &stage_coefs, stages);
            cv.extend_from_slice(&vals);
        }

        /* call fused vector operation to do the work */
        {
            let ARKodeMem { tempv1, .. } = ark_mem;
            let mut xr: Vec<&NVector> = Vec::new();
            for j in 0..step_mem.stages as usize {
                if step_mem.explicit {
                    xr.push(&step_mem.Fe[j]);
                }
                if step_mem.implicit {
                    xr.push(&step_mem.Fi[j]);
                }
            }
            if step_mem.expforcing || step_mem.impforcing {
                for v in step_mem.forcing.iter() {
                    xr.push(v);
                }
            }
            N_VLinearCombination(cv.len() as i32, &cv, &xr, tempv1);
        }

        /* fill error norm */
        *dsmPtr = N_VWrmsNorm(&ark_mem.tempv1, &ark_mem.ewt);
    }

    ARK_SUCCESS
}

/*===============================================================
  Tests
  ===============================================================*/
#[cfg(test)]
mod tests {
    use super::*;
    use crate::arkode::{ARKodeEvolve, ARKodeSStolerances};
    use crate::arkode_arkstep_io::{ARKStepSetTableName, ARKStepSetTables};
    use crate::arkode_butcher_dirk::ARKodeButcherTable_LoadDIRK;
    use crate::sundials_context::SUNContext_Create;

    /* dy/dt = lambda*(y - atan(t)) + 1/(1+t^2), y(0)=0, y(t)=atan(t) */
    fn fi_stiff(t: f64, y: &NVector, ydot: &mut NVector, _ud: &mut UserData) -> i32 {
        let lambda = -10.0;
        ydot.data[0] = lambda * (y.data[0] - t.atan()) + 1.0 / (1.0 + t * t);
        0
    }

    /* ARKStepSetTableName with an implicit-only selection replaces the
       default DIRK table; the integration then uses the requested
       method (checked via the stored q/p and a short solve). */
    #[test]
    fn arkstep_set_table_name_dirk() {
        let ctx = SUNContext_Create();
        let mut y = NVector::new(1);
        y.data[0] = 0.0;
        let mut ark = ARKStepCreate(None, Some(fi_stiff), 0.0, &y, &ctx).unwrap();
        assert_eq!(ARKodeSStolerances(&mut ark, 1.0e-4, 1.0e-8), ARK_SUCCESS);

        let flag = ARKStepSetTableName(&mut ark, "ARKODE_SDIRK_2_1_2", "ARKODE_ERK_NONE");
        assert_eq!(flag, ARK_SUCCESS);
        {
            let sm = arkStep_AccessStepMem(&mut ark, "test").unwrap();
            assert_eq!((sm.q, sm.p, sm.stages), (2, 1, 2));
            assert!(sm.implicit && !sm.explicit);
            ark.step_mem = Some(sm);
        }

        /* attach dense linear solver and integrate to t=1 */
        let a = crate::sunmatrix_dense::SUNDenseMatrix(1, 1, &ctx);
        let ls = crate::sunlinsol_dense::SUNLinSol_Dense(&y, &a, &ctx);
        assert_eq!(
            crate::arkode_ls::ARKodeSetLinearSolver(&mut ark, ls, Some(a)),
            crate::arkode_ls_impl::ARKLS_SUCCESS
        );
        let mut t = 0.0;
        let flag = ARKodeEvolve(&mut ark, 1.0, &mut y, &mut t, ARK_NORMAL);
        assert!(flag >= 0, "ARKodeEvolve flag = {}", flag);
        assert!(SUNRabs(y.data[0] - 1.0f64.atan()) < 1.0e-3);
    }

    /* Fixed non-identity mass matrix: solve M*y' = fi(t,y) with
       M = 2*I and fi = -2*y, i.e. y' = -y with y(0)=1, so
       y(t) = exp(-t).  Exercises ARKodeSetMassLinearSolver/SetMassFn
       (direct dense), the MASS_FIXED stage/residual/solution paths
       (A = M - gamma*J in arkLsLinSys, M*(yn-zp) in StageSetup,
       NlsResidual_MassFixed, ComputeSolutions_MassFixed) and the
       mass statistics. */
    #[test]
    fn arkstep_fixed_mass_matrix() {
        use crate::arkode_ls::{
            ARKodeGetNumMassSetups, ARKodeGetNumMassSolves, ARKodeSetJacFn,
            ARKodeSetLinearSolver, ARKodeSetMassFn, ARKodeSetMassLinearSolver,
        };

        fn fi_mass(_t: f64, y: &NVector, ydot: &mut NVector, _ud: &mut UserData) -> i32 {
            ydot.data[0] = -2.0 * y.data[0];
            0
        }
        fn jac_mass(
            _t: f64,
            _y: &NVector,
            _fy: &NVector,
            j: &mut crate::sundials_matrix::SUNMatrix,
            _ud: &mut UserData,
            _t1: &mut NVector,
            _t2: &mut NVector,
            _t3: &mut NVector,
        ) -> i32 {
            if let crate::sundials_matrix::SUNMatrix::Dense(dm) = j {
                dm.data[0] = -2.0;
            }
            0
        }
        fn mass_fn(
            _t: f64,
            m: &mut crate::sundials_matrix::SUNMatrix,
            _ud: &mut UserData,
            _t1: &mut NVector,
            _t2: &mut NVector,
            _t3: &mut NVector,
        ) -> i32 {
            if let crate::sundials_matrix::SUNMatrix::Dense(dm) = m {
                dm.data[0] = 2.0;
            }
            0
        }

        let ctx = SUNContext_Create();
        let mut y = NVector::new(1);
        y.data[0] = 1.0;
        let mut ark = ARKStepCreate(None, Some(fi_mass), 0.0, &y, &ctx).unwrap();
        assert_eq!(ARKodeSStolerances(&mut ark, 1.0e-6, 1.0e-10), ARK_SUCCESS);

        /* system solver (dense) */
        let a = crate::sunmatrix_dense::SUNDenseMatrix(1, 1, &ctx);
        let ls = crate::sunlinsol_dense::SUNLinSol_Dense(&y, &a, &ctx);
        assert_eq!(
            ARKodeSetLinearSolver(&mut ark, ls, Some(a)),
            crate::arkode_ls_impl::ARKLS_SUCCESS
        );
        assert_eq!(ARKodeSetJacFn(&mut ark, Some(jac_mass)), 0);

        /* mass solver (dense, time-independent) */
        let m = crate::sunmatrix_dense::SUNDenseMatrix(1, 1, &ctx);
        let mls = crate::sunlinsol_dense::SUNLinSol_Dense(&y, &m, &ctx);
        assert_eq!(
            ARKodeSetMassLinearSolver(&mut ark, mls, Some(m), false),
            crate::arkode_ls_impl::ARKLS_SUCCESS
        );
        assert_eq!(ARKodeSetMassFn(&mut ark, Some(mass_fn)), 0);

        let mut t = 0.0;
        let flag = ARKodeEvolve(&mut ark, 1.0, &mut y, &mut t, ARK_NORMAL);
        assert!(flag >= 0, "ARKodeEvolve flag = {}", flag);
        assert!(
            SUNRabs(y.data[0] - (-1.0f64).exp()) < 1.0e-5,
            "y(1) = {} vs exp(-1) = {}",
            y.data[0],
            (-1.0f64).exp()
        );

        /* mass statistics were exercised */
        let (mut nmsetups, mut nmsolves) = (0i64, 0i64);
        assert_eq!(ARKodeGetNumMassSetups(&mut ark, &mut nmsetups), 0);
        assert_eq!(ARKodeGetNumMassSolves(&mut ark, &mut nmsolves), 0);
        assert!(nmsetups > 0, "nmsetups = {}", nmsetups);
        assert!(nmsolves > 0, "nmsolves = {}", nmsolves);
    }

    /* ARKStepSetTables (implicit-only) copies a user table into step
       memory and switches to purely implicit mode. */
    #[test]
    fn arkstep_set_tables_copy() {
        let ctx = SUNContext_Create();
        let mut y = NVector::new(1);
        y.data[0] = 0.0;
        let mut ark = ARKStepCreate(None, Some(fi_stiff), 0.0, &y, &ctx).unwrap();

        let bi = ARKodeButcherTable_LoadDIRK(
            crate::arkode_butcher_dirk::ARKODE_BACKWARD_EULER_1_1,
        )
        .unwrap();
        let flag = ARKStepSetTables(&mut ark, 1, 0, Some(&bi), None);
        assert_eq!(flag, ARK_SUCCESS);
        let sm = arkStep_AccessStepMem(&mut ark, "test").unwrap();
        assert_eq!((sm.q, sm.stages), (1, 1));
        assert!(sm.Bi.is_some() && sm.Be.is_none());
        ark.step_mem = Some(sm);
    }
}

/*---------------------------------------------------------------
  arkStep_ComputeSolutions_MassFixed

  This routine calculates the final RK solution using the existing
  data.  This solution is placed directly in ark_ycur.  This
  routine also computes the error estimate ||y-ytilde||_WRMS,
  where ytilde is the embedded solution, and the norm weights come
  from ark_ewt.  This norm value is returned.  The vector form of
  this estimated error (y-ytilde) is stored in ark_mem->tempv1, in
  case the calling routine wishes to examine the error locations.

  This version assumes a fixed mass matrix.
  ---------------------------------------------------------------*/
pub(crate) fn arkStep_ComputeSolutions_MassFixed(
    ark_mem: &mut ARKodeMem,
    step_mem: &mut ARKodeARKStepMem,
    dsmPtr: &mut f64,
) -> i32 {
    /* initialize output */
    *dsmPtr = ZERO;

    /* check if the method is stiffly accurate */
    let mut stiffly_accurate = true;
    if step_mem.explicit
        && !ARKodeButcherTable_IsStifflyAccurate(step_mem.Be.as_ref().unwrap())
    {
        stiffly_accurate = false;
    }
    if step_mem.implicit
        && !ARKodeButcherTable_IsStifflyAccurate(step_mem.Bi.as_ref().unwrap())
    {
        stiffly_accurate = false;
    }

    /* If the method is stiffly accurate, ycur is already the new solution */

    if !stiffly_accurate {
        /* compute y RHS (store in y) */
        /*   set arrays for fused vector operation */
        let mut cv: Vec<f64> = Vec::new();
        for j in 0..step_mem.stages as usize {
            if step_mem.explicit {
                /* Explicit pieces */
                cv.push(ark_mem.h * step_mem.Be.as_ref().unwrap().b[j]);
            }
            if step_mem.implicit {
                /* Implicit pieces */
                cv.push(ark_mem.h * step_mem.Bi.as_ref().unwrap().b[j]);
            }
        }

        /*   call fused vector operation to compute RHS */
        {
            let ARKodeMem { ycur, .. } = ark_mem;
            let mut xr: Vec<&NVector> = Vec::new();
            for j in 0..step_mem.stages as usize {
                if step_mem.explicit {
                    xr.push(&step_mem.Fe[j]);
                }
                if step_mem.implicit {
                    xr.push(&step_mem.Fi[j]);
                }
            }
            N_VLinearCombination(cv.len() as i32, &cv, &xr, ycur);
        }

        /* solve for y update (stored in y) */
        {
            let msolve = step_mem.msolve.unwrap();
            let tol = step_mem.nlscoef;
            let mut y = std::mem::take(&mut ark_mem.ycur);
            let retval = msolve(ark_mem, &mut y, tol);
            ark_mem.ycur = y;
            if retval < 0 {
                *dsmPtr = 2.0; /* indicate too much error, step with smaller step */
                let ARKodeMem { yn, ycur, .. } = ark_mem;
                ycur.data.copy_from_slice(&yn.data); /* place old solution into y */
                return CONV_FAIL;
            }
        }

        /* compute y = yn + update */
        {
            let ARKodeMem { yn, ycur, .. } = ark_mem;
            ycur.linear_sum_with(ONE, ONE, yn);
        }

        if let Some(post) = ark_mem.PostProcessStepFn {
            let ARKodeMem { ycur, user_data, tcur, .. } = ark_mem;
            let retval = post(*tcur, ycur, user_data);
            if retval != 0 {
                return ARK_POSTPROCESS_STEP_FAIL;
            }
        }
    }

    /* compute yerr (if step adaptivity enabled) */
    if !ark_mem.fixedstep {
        /* compute yerr RHS vector */
        /*   set arrays for fused vector operation */
        let mut cv: Vec<f64> = Vec::new();
        for j in 0..step_mem.stages as usize {
            if step_mem.explicit {
                /* Explicit pieces */
                let be = step_mem.Be.as_ref().unwrap();
                let d = be.d.as_ref().unwrap();
                cv.push(ark_mem.h * (be.b[j] - d[j]));
            }
            if step_mem.implicit {
                /* Implicit pieces */
                let bi = step_mem.Bi.as_ref().unwrap();
                let d = bi.d.as_ref().unwrap();
                cv.push(ark_mem.h * (bi.b[j] - d[j]));
            }
        }

        /*   call fused vector operation to compute yerr RHS */
        {
            let ARKodeMem { tempv1, .. } = ark_mem;
            let mut xr: Vec<&NVector> = Vec::new();
            for j in 0..step_mem.stages as usize {
                if step_mem.explicit {
                    xr.push(&step_mem.Fe[j]);
                }
                if step_mem.implicit {
                    xr.push(&step_mem.Fi[j]);
                }
            }
            N_VLinearCombination(cv.len() as i32, &cv, &xr, tempv1);
        }

        /* solve for yerr */
        {
            let msolve = step_mem.msolve.unwrap();
            let tol = step_mem.nlscoef;
            let mut yerr = std::mem::take(&mut ark_mem.tempv1);
            let retval = msolve(ark_mem, &mut yerr, tol);
            ark_mem.tempv1 = yerr;
            if retval < 0 {
                *dsmPtr = 2.0; /* next attempt will reduce step by 'etacf';
                               insert dsmPtr placeholder here */
                return CONV_FAIL;
            }
        }
        /* fill error norm */
        *dsmPtr = N_VWrmsNorm(&ark_mem.tempv1, &ark_mem.ewt);
    }

    ARK_SUCCESS
}
