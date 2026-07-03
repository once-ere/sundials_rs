/* -----------------------------------------------------------------
 * Translated from src/kinsol/kinsol.c (KINSOL 7.7.0).
 * Main KINSol solver: creation/initialization (KINCreate, KINInit),
 * the KINSol driver with its step strategies (KINLinSolDrv,
 * KINFullNewton, KINLineSearch, KINConstraint), stopping tests
 * (KINStop, KINForcingTerm), scaled norms (KINScFNorm, KINScSNorm),
 * the Picard and fixed-point drivers (KINPicardAA, KINPicardFcnEval,
 * KINFP) with Anderson acceleration (AndersonAcc,
 * AndersonAccQRDelete), verbose output (KINPrintInfo) and
 * deallocation (KINFree).
 *
 * Conventions (donor cvode.rs):
 *  - The C `kin_linit/kin_lsetup/kin_lsolve` function pointers are
 *    the LsModule enum (kinsol_impl.rs), take()n out of KINMem for
 *    the duration of each call; `kin_lsetup != NULL` maps to
 *    kin_has_lsetup() (kinLsInitialize records the C NULLing of the
 *    pointer in KINLsMem.setup_disabled).
 *  - In C, KINSol aliases kin_uu with the caller's `u` and stores
 *    the caller's scaling vectors; here `u` is copied in on entry
 *    and the final kin_uu is copied back on every return path
 *    (rule 5), while u_scale/f_scale are cloned into KINMem.
 *  - All KINPrintInfo call sites in kinsol.c sit inside
 *    `#if SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_INFO`, which
 *    the reference (default) build excludes — they are compiled out
 *    here as well (marked in place). KINPrintInfo itself is exported
 *    because kinsol_ls.c calls it unconditionally.
 *  - The C fused N_Vector ops (N_VLinearCombination, N_VDotProdMulti
 *    and the kin_cv/kin_Xv scratch arrays) are reproduced inline
 *    with identical floating-point operation order.
 * -----------------------------------------------------------------*/
use crate::kinsol_aa::KINInitAA;
use crate::kinsol_impl::*;
use crate::kinsol_ls::{kinLsInit, kinLsSetup, kinLsSolve};
use crate::kinsol_orth::KINInitOrth;
use crate::nvector_serial::*;
use crate::sundials_context::SUNContext;
use crate::sundials_math::*;
use crate::sundials_types::*;

/*=================================================================*/
/* KINSOL Private Constants                                        */
/*=================================================================*/

const HALF: f64 = 0.5;
const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;
const ONEPT5: f64 = 1.5;
const TWO: f64 = 2.0;
const THREE: f64 = 3.0;
const FIVE: f64 = 5.0;
const TWELVE: f64 = 12.0;
const POINT1: f64 = 0.1;
const POINT01: f64 = 0.01;
const POINT99: f64 = 0.99;
const THOUSAND: f64 = 1000.0;
const POINT9: f64 = 0.9;
const POINT0001: f64 = 0.0001;

/*=================================================================*/
/* KINSOL Routine-Specific Constants                               */
/*=================================================================*/

/*
 * Control constants for lower-level functions used by KINSol
 * ----------------------------------------------------------
 *
 * KINStop return value requesting more iterations
 *    RETRY_ITERATION
 *    CONTINUE_ITERATIONS
 *
 * KINFullNewton, KINLineSearch, KINFP, and KINPicardAA return values:
 *    KIN_SUCCESS
 *    KIN_SYSFUNC_FAIL
 *    STEP_TOO_SMALL
 *
 * KINConstraint return values:
 *    KIN_SUCCESS
 *    CONSTR_VIOLATED
 */

pub const RETRY_ITERATION: i32 = -998;
pub const CONTINUE_ITERATIONS: i32 = -999;
pub const STEP_TOO_SMALL: i32 = -997;
pub const CONSTR_VIOLATED: i32 = -996;

/*
 * Algorithmic constants
 * ---------------------
 *
 * MAX_RECVR   max. no. of attempts to correct a recoverable func error
 */

const MAX_RECVR: i32 = 5;

/*
 * Keys for KINPrintInfo (kinsol.c defines keys 2..13 only when
 * SUNDIALS_LOGGING_LEVEL >= SUNDIALS_LOGGING_INFO; they are kept
 * here unconditionally for API completeness — PRNT_NLI/PRNT_EPS
 * live in kinsol_ls_impl.rs)
 */

pub const PRNT_RETVAL: i32 = 1;
pub const PRNT_NNI: i32 = 2;
pub const PRNT_TOL: i32 = 3;
pub const PRNT_FMAX: i32 = 4;
pub const PRNT_PNORM: i32 = 5;
pub const PRNT_PNORM1: i32 = 6;
pub const PRNT_FNORM: i32 = 7;
pub const PRNT_LAM: i32 = 8;
pub const PRNT_ALPHA: i32 = 9;
pub const PRNT_BETA: i32 = 10;
pub const PRNT_ALPHABETA: i32 = 11;
pub const PRNT_ADJ: i32 = 12;
pub const PRNT_OTHER: i32 = 13;

/*
 * =================================================================
 * Exported Functions Implementation
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Creation and allocation functions
 * -----------------------------------------------------------------
 */

/*
 * Function : KINCreate
 *
 * KINCreate creates an internal memory block for a problem to
 * be solved by KINSOL. (C returns NULL on a NULL sunctx or a failed
 * malloc; neither can occur in safe Rust — the &SUNContext cannot be
 * null and allocation is infallible.) The default optional-input
 * values set here in C live in KINMem::default() (kinsol_impl.rs).
 */
pub fn KINCreate(sunctx: &SUNContext) -> Box<KINMem> {
    Box::new(KINMem {
        kin_sunctx: sunctx.clone(),
        ..KINMem::default()
    })
}

/*
 * Function : KINInit
 *
 * KINInit allocates memory for a problem or execution of KINSol.
 * If memory is successfully allocated, KIN_SUCCESS is returned.
 * Otherwise, an error message is printed and an error flag
 * returned.
 *
 * (C also rejects kinmem == NULL / func == NULL; a &mut KINMem and a
 * plain fn pointer cannot be null, so MSG_NO_MEM / MSG_FUNC_NULL are
 * unreachable.)
 */
pub fn KINInit(kin_mem: &mut KINMem, func: KINSysFn, tmpl: &NVector) -> i32 {
    /* check if all required vector operations are implemented */
    let nvectorOK = KINCheckNvector(tmpl);
    if !nvectorOK {
        KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINInit", file!(),
                        MSG_BAD_NVECTOR);
        return KIN_ILL_INPUT;
    }

    /* set space requirements for one N_Vector */
    let (lrw1, liw1) = N_VSpace(tmpl);
    kin_mem.kin_lrw1 = lrw1;
    kin_mem.kin_liw1 = liw1;

    /* allocate necessary vectors */
    let allocOK = KINAllocVectors(kin_mem, tmpl);
    if !allocOK {
        KINProcessError(Some(kin_mem), KIN_MEM_FAIL, line!(), "KINInit", file!(),
                        MSG_MEM_FAIL);
        return KIN_MEM_FAIL;
    }

    /* copy the input parameter into KINSol state */
    kin_mem.kin_func = Some(func);

    /* set the linear solver addresses to NULL */
    kin_mem.kin_lmem = LsModule::None;

    /* problem memory has been successfully allocated */
    kin_mem.kin_MallocDone = SUNTRUE;

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Linear solver module dispatch (C: kin_linit/kin_lsetup/kin_lsolve
 * function pointers). The module is take()n out of KINMem so the
 * KINLS routines can borrow the solver memory mutably (donor
 * cvode.rs pattern).
 * -----------------------------------------------------------------
 */

/// Mirrors C's `kin_mem->kin_lsetup != NULL`: kinLsInitialize NULLs the
/// lsetup pointer for matrix-free configurations without a full
/// preconditioner (KINLsMem.setup_disabled carries that state).
pub fn kin_has_lsetup(kin_mem: &KINMem) -> bool {
    match &kin_mem.kin_lmem {
        LsModule::None => false,
        LsModule::Ls(ls) => !ls.setup_disabled,
    }
}

fn kin_linit_dispatch(kin_mem: &mut KINMem) -> i32 {
    let mut lmem = std::mem::take(&mut kin_mem.kin_lmem);
    let retval = match &mut lmem {
        LsModule::None => 0,
        LsModule::Ls(ls) => kinLsInit(kin_mem, ls),
    };
    kin_mem.kin_lmem = lmem;
    retval
}

fn kin_lsetup_dispatch(kin_mem: &mut KINMem) -> i32 {
    let mut lmem = std::mem::take(&mut kin_mem.kin_lmem);
    let retval = match &mut lmem {
        LsModule::None => 0,
        LsModule::Ls(ls) => kinLsSetup(kin_mem, ls),
    };
    kin_mem.kin_lmem = lmem;
    retval
}

/// C: `kin_lsolve(kin_mem, xx, bb, &sJpnorm, &sFdotJp)`. The scalar
/// outputs are threaded through locals seeded from and stored back to
/// the KINMem fields the C code passes by address.
fn kin_lsolve_dispatch(kin_mem: &mut KINMem, xx: &mut NVector, bb: &mut NVector) -> i32 {
    let mut lmem = std::mem::take(&mut kin_mem.kin_lmem);
    let mut sJpnorm = kin_mem.kin_sJpnorm;
    let mut sFdotJp = kin_mem.kin_sFdotJp;
    let retval = match &mut lmem {
        LsModule::None => {
            /* C would dereference a NULL kin_lsolve; report and fail */
            KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINLinSolDrv", file!(),
                            MSG_LSOLV_NO_MEM);
            -1
        }
        LsModule::Ls(ls) => kinLsSolve(kin_mem, ls, xx, bb, &mut sJpnorm, &mut sFdotJp),
    };
    kin_mem.kin_lmem = lmem;
    kin_mem.kin_sJpnorm = sJpnorm;
    kin_mem.kin_sFdotJp = sFdotJp;
    retval
}

/*
 * -----------------------------------------------------------------
 * Main solver function
 * -----------------------------------------------------------------
 */

/*
 * Function : KINSol
 *
 * KINSol (main KINSOL driver routine) manages the computational
 * process of computing an approximate solution of the nonlinear
 * system F(uu) = 0. The KINSol routine calls the following
 * subroutines:
 *
 *  KINSolInit    checks if initial guess satisfies user-supplied
 *                constraints and initializes linear solver
 *
 *  KINLinSolDrv  interfaces with linear solver to find a
 *                solution of the system J(uu)*x = b (calculate
 *                Newton step)
 *
 *  KINFullNewton/KINLineSearch  implement the global strategy
 *
 *  KINForcingTerm  computes the forcing term (eta)
 *
 *  KINStop  determines if an approximate solution has been found
 */
pub fn KINSol(
    kin_mem: &mut KINMem,
    u: &mut NVector,
    strategy_in: i32,
    u_scale: &NVector,
    f_scale: &NVector,
) -> i32 {
    /* check for kinmem non-NULL: guaranteed by &mut KINMem */

    if kin_mem.kin_MallocDone == SUNFALSE {
        KINProcessError(Some(kin_mem), KIN_NO_MALLOC, line!(), "KINSol", file!(),
                        MSG_NO_MALLOC);
        return KIN_NO_MALLOC;
    }

    /* load input arguments */
    kin_mem.kin_uu = u.clone();
    kin_mem.kin_uscale = u_scale.clone();
    kin_mem.kin_fscale = f_scale.clone();
    kin_mem.kin_globalstrategy = strategy_in;

    let ret = KINSol_body(kin_mem);

    /* In C kin_uu aliases the caller's u for the whole solve; copy the
       final iterate back on every return path past this point (rule 5) */
    *u = kin_mem.kin_uu.clone();

    ret
}

/* Body of KINSol after the input arguments are loaded (single exit so
   the kin_uu -> u copy-back above covers every C return path). */
fn KINSol_body(kin_mem: &mut KINMem) -> i32 {
    /* initialize to avoid compiler warning messages */
    let mut maxStepTaken = SUNFALSE;
    let mut fnormp = -ONE;
    let mut f1normp = -ONE;

    /* initialize epsmin to avoid compiler warning message */
    let mut epsmin = ZERO;

    /* Setup Anderson acceleration for FP or Picard */
    if (kin_mem.kin_globalstrategy == KIN_FP || kin_mem.kin_globalstrategy == KIN_PICARD)
        && kin_mem.kin_m_aa != 0
    {
        /* Initialize Anderson acceleration workspace */
        let ret = KINInitAA(kin_mem);
        if ret != 0 {
            KINProcessError(Some(kin_mem), ret, line!(), "KINSol", file!(),
                            "Initializing Anderson acceleration failed");
            return ret;
        }

        /* Initialize orthogonalization workspace */
        let ret = KINInitOrth(kin_mem);
        if ret != 0 {
            KINProcessError(Some(kin_mem), ret, line!(), "KINSol", file!(),
                            "Initializing the orthogonalization method failed");
            return ret;
        }
    }

    /* CSW:
       Call fixed point solver if requested.  Note that this should probably
       be forked off to a FPSOL solver instead of kinsol in the future. */
    if kin_mem.kin_globalstrategy == KIN_FP {
        if kin_mem.kin_uu.is_empty() {
            KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSol", file!(),
                            MSG_UU_NULL);
            return KIN_ILL_INPUT;
        }

        if kin_mem.kin_constraintsSet != SUNFALSE {
            KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSol", file!(),
                            MSG_CONSTRAINTS_NOTOK);
            return KIN_ILL_INPUT;
        }

        /* (KINPrintInfo PRNT_TOL: SUNDIALS_LOGGING_LEVEL >= INFO only) */

        kin_mem.kin_nfe = 0;
        kin_mem.kin_nnilset = 0;
        kin_mem.kin_nnilset_sub = 0;
        kin_mem.kin_nni = 0;
        kin_mem.kin_nbcf = 0;
        kin_mem.kin_nbktrk = 0;

        let ret = KINFP(kin_mem);

        match ret {
            KIN_SYSFUNC_FAIL => {
                KINProcessError(Some(kin_mem), KIN_SYSFUNC_FAIL, line!(), "KINSol", file!(),
                                MSG_SYSFUNC_FAILED);
            }
            KIN_MAXITER_REACHED => {
                KINProcessError(Some(kin_mem), KIN_MAXITER_REACHED, line!(), "KINSol",
                                file!(), MSG_MAXITER_REACHED);
            }
            _ => {}
        }

        return ret;
    }

    /* initialize solver */
    let ret = KINSolInit(kin_mem);
    if ret != KIN_SUCCESS {
        return ret;
    }

    kin_mem.kin_ncscmx = 0;

    /* Note: The following logic allows the choice of whether or not
       to force a call to the linear solver setup upon a given call to
       KINSol */

    if kin_mem.kin_noInitSetup {
        kin_mem.kin_sthrsh = ONE;
    } else {
        kin_mem.kin_sthrsh = TWO;
    }

    /* if eps is to be bounded from below, set the bound */

    if kin_mem.kin_inexact_ls && !kin_mem.kin_noMinEps {
        epsmin = POINT01 * kin_mem.kin_fnormtol;
    }

    /* if omega is zero at this point, make sure it will be evaluated
       at each iteration based on the provided min/max bounds and the
       current function norm. */
    if kin_mem.kin_omega == ZERO {
        kin_mem.kin_eval_omega = SUNTRUE;
    } else {
        kin_mem.kin_eval_omega = SUNFALSE;
    }

    /* CSW:
       Call fixed point solver for Picard method if requested.
       Note that this should probably be forked off to a part of an
       FPSOL solver instead of kinsol in the future. */
    if kin_mem.kin_globalstrategy == KIN_PICARD {
        if kin_mem.kin_gval.is_empty() {
            kin_mem.kin_gval = N_VClone(&kin_mem.kin_unew);
            kin_mem.kin_liw += kin_mem.kin_liw1;
            kin_mem.kin_lrw += kin_mem.kin_lrw1;
        }
        return KINPicardAA(kin_mem);
    }

    let mut ret = KIN_SUCCESS;

    'outer: loop {
        kin_mem.kin_retry_nni = SUNFALSE;

        kin_mem.kin_nni += 1;

        /* calculate the epsilon (stopping criteria for iterative linear solver)
           for this iteration based on eta from the routine KINForcingTerm */

        if kin_mem.kin_inexact_ls {
            kin_mem.kin_eps = (kin_mem.kin_eta + kin_mem.kin_uround) * kin_mem.kin_fnorm;
            if !kin_mem.kin_noMinEps {
                kin_mem.kin_eps = SUNMAX(epsmin, kin_mem.kin_eps);
            }
        }

        /* repeat_nni: (C goto target) */
        loop {
            /* call the appropriate routine to calculate an acceptable step pp */

            let mut sflag = 0;

            if kin_mem.kin_globalstrategy == KIN_NONE {
                /* Full Newton Step */

                /* call KINLinSolDrv to calculate the (approximate) Newton step, pp */
                ret = KINLinSolDrv(kin_mem);
                if ret != KIN_SUCCESS {
                    break 'outer;
                }

                sflag = KINFullNewton(kin_mem, &mut fnormp, &mut f1normp, &mut maxStepTaken);

                /* if sysfunc failed unrecoverably, stop */
                if sflag == KIN_SYSFUNC_FAIL || sflag == KIN_REPTD_SYSFUNC_ERR {
                    ret = sflag;
                    break 'outer;
                }
            } else if kin_mem.kin_globalstrategy == KIN_LINESEARCH {
                /* Line Search */

                /* call KINLinSolDrv to calculate the (approximate) Newton step, pp */
                ret = KINLinSolDrv(kin_mem);
                if ret != KIN_SUCCESS {
                    break 'outer;
                }

                sflag = KINLineSearch(kin_mem, &mut fnormp, &mut f1normp, &mut maxStepTaken);

                /* if sysfunc failed unrecoverably, stop */
                if sflag == KIN_SYSFUNC_FAIL || sflag == KIN_REPTD_SYSFUNC_ERR {
                    ret = sflag;
                    break 'outer;
                }

                /* if too many beta condition failures, then stop iteration */
                if kin_mem.kin_nbcf > kin_mem.kin_mxnbcf {
                    ret = KIN_LINESEARCH_BCFAIL;
                    break 'outer;
                }
            }

            if kin_mem.kin_globalstrategy != KIN_PICARD
                && kin_mem.kin_globalstrategy != KIN_FP
            {
                /* evaluate eta by calling the forcing term routine */
                if kin_mem.kin_callForcingTerm {
                    KINForcingTerm(kin_mem, fnormp);
                }

                kin_mem.kin_fnorm = fnormp;

                /* call KINStop to check if tolerances where met by this iteration */
                ret = KINStop(kin_mem, maxStepTaken, sflag);

                if ret == RETRY_ITERATION {
                    kin_mem.kin_retry_nni = SUNTRUE;
                    continue; /* goto repeat_nni */
                }
            }
            break;
        }

        /* update uu after the iteration */
        N_VScale(ONE, &kin_mem.kin_unew, &mut kin_mem.kin_uu);

        kin_mem.kin_f1norm = f1normp;

        /* print the current nni, fnorm, and nfe values */
        /* (KINPrintInfo PRNT_NNI: SUNDIALS_LOGGING_LEVEL >= INFO only) */

        if ret != CONTINUE_ITERATIONS {
            break;
        }
    } /* end of loop; return */

    /* (KINPrintInfo PRNT_RETVAL: SUNDIALS_LOGGING_LEVEL >= INFO only) */

    match ret {
        KIN_SYSFUNC_FAIL => {
            KINProcessError(Some(kin_mem), KIN_SYSFUNC_FAIL, line!(), "KINSol", file!(),
                            MSG_SYSFUNC_FAILED);
        }
        KIN_REPTD_SYSFUNC_ERR => {
            KINProcessError(Some(kin_mem), KIN_REPTD_SYSFUNC_ERR, line!(), "KINSol", file!(),
                            MSG_SYSFUNC_REPTD);
        }
        KIN_LSETUP_FAIL => {
            KINProcessError(Some(kin_mem), KIN_LSETUP_FAIL, line!(), "KINSol", file!(),
                            MSG_LSETUP_FAILED);
        }
        KIN_LSOLVE_FAIL => {
            KINProcessError(Some(kin_mem), KIN_LSOLVE_FAIL, line!(), "KINSol", file!(),
                            MSG_LSOLVE_FAILED);
        }
        KIN_LINSOLV_NO_RECOVERY => {
            KINProcessError(Some(kin_mem), KIN_LINSOLV_NO_RECOVERY, line!(), "KINSol",
                            file!(), MSG_LINSOLV_NO_RECOVERY);
        }
        KIN_LINESEARCH_NONCONV => {
            KINProcessError(Some(kin_mem), KIN_LINESEARCH_NONCONV, line!(), "KINSol",
                            file!(), MSG_LINESEARCH_NONCONV);
        }
        KIN_LINESEARCH_BCFAIL => {
            KINProcessError(Some(kin_mem), KIN_LINESEARCH_BCFAIL, line!(), "KINSol",
                            file!(), MSG_LINESEARCH_BCFAIL);
        }
        KIN_MAXITER_REACHED => {
            KINProcessError(Some(kin_mem), KIN_MAXITER_REACHED, line!(), "KINSol", file!(),
                            MSG_MAXITER_REACHED);
        }
        KIN_MXNEWT_5X_EXCEEDED => {
            KINProcessError(Some(kin_mem), KIN_MXNEWT_5X_EXCEEDED, line!(), "KINSol",
                            file!(), MSG_MXNEWT_5X_EXCEEDED);
        }
        _ => {}
    }

    ret
}

/*
 * -----------------------------------------------------------------
 * Deallocation function
 * -----------------------------------------------------------------
 */

/*
 * Function : KINFree
 *
 * This routine frees the problem memory allocated by KINInit
 * (KINFreeVectors, lfree, KINFreeAA, KINFreeOrth and free itself
 * all collapse to RAII drop; the SUNDIALS_ENABLE_PYTHON teardown is
 * excluded with the foreign-runtime backends).
 */
pub fn KINFree(_kinmem: Box<KINMem>) {}

/*
 * =================================================================
 * Private Functions
 * =================================================================
 */

/*
 * Function : KINCheckNvector
 *
 * This routine checks if all required vector operations are
 * implemented (excluding those required by KINConstraint). The
 * pure-Rust serial NVector implements every op the C code tests
 * (nvclone, nvdestroy, nvlinearsum, nvprod, nvdiv, nvscale, nvabs,
 * nvinv, nvmaxnorm, nvmin, nvwl2norm), so this always returns
 * SUNTRUE.
 */
fn KINCheckNvector(_tmpl: &NVector) -> bool {
    SUNTRUE
}

/*
 * -----------------------------------------------------------------
 * Memory allocation/deallocation
 * -----------------------------------------------------------------
 */

/*
 * Function : KINAllocVectors
 *
 * This routine allocates the KINSol vectors (unew, fval, pp, vtemp1
 * and vtemp2; df/dg/q for Anderson acceleration are allocated by
 * KINInitAA). Allocation cannot fail in safe Rust, so the C
 * partial-unwind error paths are unreachable and SUNTRUE is always
 * returned.
 */
fn KINAllocVectors(kin_mem: &mut KINMem, tmpl: &NVector) -> bool {
    if kin_mem.kin_unew.is_empty() {
        kin_mem.kin_unew = N_VClone(tmpl);
        kin_mem.kin_liw += kin_mem.kin_liw1;
        kin_mem.kin_lrw += kin_mem.kin_lrw1;
    }

    if kin_mem.kin_fval.is_empty() {
        kin_mem.kin_fval = N_VClone(tmpl);
        kin_mem.kin_liw += kin_mem.kin_liw1;
        kin_mem.kin_lrw += kin_mem.kin_lrw1;
    }

    if kin_mem.kin_pp.is_empty() {
        kin_mem.kin_pp = N_VClone(tmpl);
        kin_mem.kin_liw += kin_mem.kin_liw1;
        kin_mem.kin_lrw += kin_mem.kin_lrw1;
    }

    if kin_mem.kin_vtemp1.is_empty() {
        kin_mem.kin_vtemp1 = N_VClone(tmpl);
        kin_mem.kin_liw += kin_mem.kin_liw1;
        kin_mem.kin_lrw += kin_mem.kin_lrw1;
    }

    if kin_mem.kin_vtemp2.is_empty() {
        kin_mem.kin_vtemp2 = N_VClone(tmpl);
        kin_mem.kin_liw += kin_mem.kin_liw1;
        kin_mem.kin_lrw += kin_mem.kin_lrw1;
    }

    SUNTRUE
}

/* (KINFreeVectors — called only from KINFree in C — collapses to the
   RAII drop of KINMem and is not ported.) */

/*
 * -----------------------------------------------------------------
 * Initial setup
 * -----------------------------------------------------------------
 */

/*
 * KINSolInit
 *
 * KINSolInit initializes the problem for the specific input
 * received in this call to KINSol (which calls KINSolInit). All
 * problem specification inputs are checked for errors.
 *
 * The possible return values for KINSolInit are:
 *   KIN_SUCCESS : indicates a normal initialization
 *
 *   KIN_ILL_INPUT : indicates that an input error has been found
 *
 *   KIN_INITIAL_GUESS_OK : indicates that the guess uu
 *                          satisfied the system func(uu) = 0
 *                          within the tolerances specified
 */
fn KINSolInit(kin_mem: &mut KINMem) -> i32 {
    /* check for illegal input parameters */

    if kin_mem.kin_uu.is_empty() {
        KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSolInit", file!(),
                        MSG_UU_NULL);
        return KIN_ILL_INPUT;
    }

    /* check for valid strategy */

    if kin_mem.kin_globalstrategy != KIN_NONE
        && kin_mem.kin_globalstrategy != KIN_LINESEARCH
        && kin_mem.kin_globalstrategy != KIN_PICARD
        && kin_mem.kin_globalstrategy != KIN_FP
    {
        KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSolInit", file!(),
                        MSG_BAD_GLSTRAT);
        return KIN_ILL_INPUT;
    }

    if kin_mem.kin_uscale.is_empty() {
        KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSolInit", file!(),
                        MSG_BAD_USCALE);
        return KIN_ILL_INPUT;
    }

    if N_VMin(&kin_mem.kin_uscale) <= ZERO {
        KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSolInit", file!(),
                        MSG_USCALE_NONPOSITIVE);
        return KIN_ILL_INPUT;
    }

    if kin_mem.kin_fscale.is_empty() {
        KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSolInit", file!(),
                        MSG_BAD_FSCALE);
        return KIN_ILL_INPUT;
    }

    if N_VMin(&kin_mem.kin_fscale) <= ZERO {
        KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSolInit", file!(),
                        MSG_FSCALE_NONPOSITIVE);
        return KIN_ILL_INPUT;
    }

    if !kin_mem.kin_constraints.is_empty()
        && (kin_mem.kin_globalstrategy == KIN_PICARD
            || kin_mem.kin_globalstrategy == KIN_FP)
    {
        KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSolInit", file!(),
                        MSG_CONSTRAINTS_NOTOK);
        return KIN_ILL_INPUT;
    }

    /* set the constraints flag */

    if kin_mem.kin_constraints.is_empty() {
        kin_mem.kin_constraintsSet = SUNFALSE;
    } else {
        kin_mem.kin_constraintsSet = SUNTRUE;
        /* (the serial NVector implements nvconstrmask and nvminquotient,
           so the C MSG_BAD_NVECTOR rejection cannot fire) */
    }

    /* check the initial guess uu against the constraints */

    if kin_mem.kin_constraintsSet {
        let KINMem { kin_constraints, kin_uu, kin_vtemp1, .. } = kin_mem;
        if !N_VConstrMask(kin_constraints, kin_uu, kin_vtemp1) {
            KINProcessError(Some(kin_mem), KIN_ILL_INPUT, line!(), "KINSolInit", file!(),
                            MSG_INITIAL_CNSTRNT);
            return KIN_ILL_INPUT;
        }
    }

    /* all error checking is complete at this point */
    /* (KINPrintInfo PRNT_TOL: SUNDIALS_LOGGING_LEVEL >= INFO only) */

    /* calculate the default value for mxnewtstep (maximum Newton step) */

    if kin_mem.kin_mxnstepin == ZERO {
        kin_mem.kin_mxnewtstep =
            THOUSAND * N_VWL2Norm(&kin_mem.kin_uu, &kin_mem.kin_uscale);
    } else {
        kin_mem.kin_mxnewtstep = kin_mem.kin_mxnstepin;
    }

    if kin_mem.kin_mxnewtstep < ONE {
        kin_mem.kin_mxnewtstep = ONE;
    }

    /* additional set-up for inexact linear solvers */

    if kin_mem.kin_inexact_ls {
        /* set up the coefficients for the eta calculation */

        kin_mem.kin_callForcingTerm = kin_mem.kin_etaflag != KIN_ETACONSTANT;

        /* this value is always used for choice #1 */

        if kin_mem.kin_etaflag == KIN_ETACHOICE1 {
            kin_mem.kin_eta_alpha = (ONE + SUNRsqrt(FIVE)) * HALF;
        }

        /* initial value for eta set to 0.5 for other than the
           KIN_ETACONSTANT option */

        if kin_mem.kin_etaflag != KIN_ETACONSTANT {
            kin_mem.kin_eta = HALF;
        }

        /* disable residual monitoring if using an inexact linear solver */

        kin_mem.kin_noResMon = SUNTRUE;
    } else {
        kin_mem.kin_callForcingTerm = SUNFALSE;
    }

    /* initialize counters */

    kin_mem.kin_nfe = 0;
    kin_mem.kin_nnilset = 0;
    kin_mem.kin_nnilset_sub = 0;
    kin_mem.kin_nni = 0;
    kin_mem.kin_nbcf = 0;
    kin_mem.kin_nbktrk = 0;

    /* see if the initial guess uu satisfies the nonlinear system */
    let func = kin_mem.kin_func.unwrap();
    let retval = func(&kin_mem.kin_uu, &mut kin_mem.kin_fval, &mut kin_mem.kin_user_data);
    kin_mem.kin_nfe += 1;

    if retval < 0 {
        KINProcessError(Some(kin_mem), KIN_SYSFUNC_FAIL, line!(), "KINSolInit", file!(),
                        MSG_SYSFUNC_FAILED);
        return KIN_SYSFUNC_FAIL;
    } else if retval > 0 {
        KINProcessError(Some(kin_mem), KIN_FIRST_SYSFUNC_ERR, line!(), "KINSolInit",
                        file!(), MSG_SYSFUNC_FIRST);
        return KIN_FIRST_SYSFUNC_ERR;
    }

    let fmax = {
        let KINMem { kin_vtemp1, kin_fval, kin_fscale, .. } = kin_mem;
        KINScFNorm(kin_vtemp1, kin_fval, kin_fscale)
    };
    if fmax <= POINT01 * kin_mem.kin_fnormtol {
        kin_mem.kin_fnorm = N_VWL2Norm(&kin_mem.kin_fval, &kin_mem.kin_fscale);
        return KIN_INITIAL_GUESS_OK;
    }

    /* (KINPrintInfo PRNT_FMAX: SUNDIALS_LOGGING_LEVEL >= INFO only) */

    /* initialize the linear solver if linit != NULL */

    if !kin_mem.kin_lmem.is_none() {
        let retval = kin_linit_dispatch(kin_mem);
        if retval != 0 {
            KINProcessError(Some(kin_mem), KIN_LINIT_FAIL, line!(), "KINSolInit", file!(),
                            MSG_LINIT_FAIL);
            return KIN_LINIT_FAIL;
        }
    }

    /* initialize the L2 (Euclidean) norms of f for the linear iteration steps */

    kin_mem.kin_fnorm = N_VWL2Norm(&kin_mem.kin_fval, &kin_mem.kin_fscale);
    kin_mem.kin_f1norm = HALF * kin_mem.kin_fnorm * kin_mem.kin_fnorm;
    kin_mem.kin_fnorm_sub = kin_mem.kin_fnorm;
    /* (KINPrintInfo PRNT_NNI: SUNDIALS_LOGGING_LEVEL >= INFO only) */

    /* problem has now been successfully initialized */

    KIN_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Step functions
 * -----------------------------------------------------------------
 */

/*
 * KINLinSolDrv
 *
 * This routine handles the process of solving for the approximate
 * solution of the Newton equations in the Newton iteration.
 * Subsequent routines handle the nonlinear aspects of its
 * application.
 */
fn KINLinSolDrv(kin_mem: &mut KINMem) -> i32 {
    if (kin_mem.kin_nni - kin_mem.kin_nnilset) >= kin_mem.kin_msbset {
        kin_mem.kin_sthrsh = TWO;
        kin_mem.kin_update_fnorm_sub = SUNTRUE;
    }

    loop {
        kin_mem.kin_jacCurrent = SUNFALSE;

        if kin_mem.kin_sthrsh > ONEPT5 && kin_has_lsetup(kin_mem) {
            let retval = kin_lsetup_dispatch(kin_mem);
            kin_mem.kin_jacCurrent = SUNTRUE;
            kin_mem.kin_nnilset = kin_mem.kin_nni;
            kin_mem.kin_nnilset_sub = kin_mem.kin_nni;
            if retval != 0 {
                return KIN_LSETUP_FAIL;
            }
        }

        /* rename vectors for readability */

        let mut b = std::mem::take(&mut kin_mem.kin_unew);
        let mut x = std::mem::take(&mut kin_mem.kin_pp);

        /* load b with the current value of -fval */

        N_VScale(-ONE, &kin_mem.kin_fval, &mut b);

        /* call the generic 'lsolve' routine to solve the system Jx = b */

        let retval = kin_lsolve_dispatch(kin_mem, &mut x, &mut b);

        kin_mem.kin_unew = b;
        kin_mem.kin_pp = x;

        if retval == 0 {
            return KIN_SUCCESS;
        } else if retval < 0 {
            return KIN_LSOLVE_FAIL;
        } else if !kin_has_lsetup(kin_mem) || kin_mem.kin_jacCurrent {
            return KIN_LINSOLV_NO_RECOVERY;
        }

        /* loop back only if the linear solver setup is in use
           and Jacobian information is not current */

        kin_mem.kin_sthrsh = TWO;
    }
}

/*
 * KINFullNewton
 *
 * This routine is the main driver for the Full Newton
 * algorithm. Its purpose is to compute unew = uu + pp in the
 * direction pp from uu, taking the full Newton step. The
 * step may be constrained if the constraint conditions are
 * violated, or if the norm of pp is greater than mxnewtstep.
 */
fn KINFullNewton(
    kin_mem: &mut KINMem,
    fnormp: &mut f64,
    f1normp: &mut f64,
    maxStepTaken: &mut bool,
) -> i32 {
    *maxStepTaken = SUNFALSE;
    let mut pnorm = N_VWL2Norm(&kin_mem.kin_pp, &kin_mem.kin_uscale);
    let mut ratio = ONE;
    if pnorm > kin_mem.kin_mxnewtstep {
        ratio = kin_mem.kin_mxnewtstep / pnorm;
        /* N_VScale(ratio, pp, pp) — output aliases the input */
        kin_mem.kin_pp.scale_inplace(ratio);
        pnorm = kin_mem.kin_mxnewtstep;
    }
    /* (KINPrintInfo PRNT_PNORM: SUNDIALS_LOGGING_LEVEL >= INFO only) */

    /* If constraints are active, then constrain the step accordingly */

    kin_mem.kin_stepl = pnorm;
    kin_mem.kin_stepmul = ONE;
    if kin_mem.kin_constraintsSet {
        let retval = KINConstraint(kin_mem);
        if retval == CONSTR_VIOLATED {
            /* Apply stepmul set in KINConstraint */
            let stepmul = kin_mem.kin_stepmul;
            ratio *= stepmul;
            kin_mem.kin_pp.scale_inplace(stepmul);
            pnorm *= stepmul;
            kin_mem.kin_stepl = pnorm;
            /* (KINPrintInfo PRNT_PNORM: SUNDIALS_LOGGING_LEVEL >= INFO only) */
            if pnorm <= kin_mem.kin_scsteptol {
                N_VLinearSum(ONE, &kin_mem.kin_uu, ONE, &kin_mem.kin_pp,
                             &mut kin_mem.kin_unew);
                return STEP_TOO_SMALL;
            }
        }
    }

    /* Attempt (at most MAX_RECVR times) to evaluate function at the new iterate */

    let mut fOK = SUNFALSE;
    let func = kin_mem.kin_func.unwrap();

    for _ircvr in 1..=MAX_RECVR {
        /* compute the iterate unew = uu + pp */
        N_VLinearSum(ONE, &kin_mem.kin_uu, ONE, &kin_mem.kin_pp, &mut kin_mem.kin_unew);

        /* evaluate func(unew) and its norm, and return */
        let retval =
            func(&kin_mem.kin_unew, &mut kin_mem.kin_fval, &mut kin_mem.kin_user_data);
        kin_mem.kin_nfe += 1;

        /* if func was successful, accept pp */
        if retval == 0 {
            fOK = SUNTRUE;
            break;
        }
        /* if func failed unrecoverably, give up */
        else if retval < 0 {
            return KIN_SYSFUNC_FAIL;
        }

        /* func failed recoverably; cut step in half and try again */
        ratio *= HALF;
        kin_mem.kin_pp.scale_inplace(HALF);
        pnorm *= HALF;
        kin_mem.kin_stepl = pnorm;
    }

    /* If func() failed recoverably MAX_RECVR times, give up */

    if !fOK {
        return KIN_REPTD_SYSFUNC_ERR;
    }

    /* Evaluate function norms */

    *fnormp = N_VWL2Norm(&kin_mem.kin_fval, &kin_mem.kin_fscale);
    *f1normp = HALF * (*fnormp) * (*fnormp);

    /* scale sFdotJp and sJpnorm by ratio for later use in KINForcingTerm */

    kin_mem.kin_sFdotJp *= ratio;
    kin_mem.kin_sJpnorm *= ratio;

    /* (KINPrintInfo PRNT_FNORM: SUNDIALS_LOGGING_LEVEL >= INFO only) */

    if pnorm > POINT99 * kin_mem.kin_mxnewtstep {
        *maxStepTaken = SUNTRUE;
    }

    KIN_SUCCESS
}

/*
 * KINLineSearch
 *
 * The routine KINLineSearch implements the LineSearch algorithm.
 * Its purpose is to find unew = uu + rl * pp in the direction pp
 * from uu so that:
 *                                    t
 *  func(unew) <= func(uu) + alpha * g  (unew - uu) (alpha = 1.e-4)
 *
 *    and
 *                                   t
 *  func(unew) >= func(uu) + beta * g  (unew - uu) (beta = 0.9)
 *
 * where 0 < rlmin <= rl <= rlmax.
 *
 * Note:
 *             mxnewtstep
 *  rlmax = ----------------   if uu+pp is feasible
 *          ||uscale*pp||_L2
 *
 *  rlmax = 1   otherwise
 *
 *    and
 *
 *                 scsteptol
 *  rlmin = --------------------------
 *          ||           pp         ||
 *          || -------------------- ||_L-infinity
 *          || (1/uscale + SUNRabs(uu)) ||
 *
 *
 * If the system function fails unrecoverably at any time, KINLineSearch
 * returns KIN_SYSFUNC_FAIL which will halt the solver.
 *
 * We attempt to correct recoverable system function failures only before
 * the alpha-condition loop; i.e. when the solution is updated with the
 * full Newton step (possibly reduced due to constraint violations).
 * Once we find a feasible pp, we assume that any update up to pp is
 * feasible.
 *
 * If the step size is limited due to constraint violations and/or
 * recoverable system function failures, we set rlmax=1 to ensure
 * that the update remains feasible during the attempts to enforce
 * the beta-condition (this is not an issue while enforcing the alpha
 * condition, as rl can only decrease from 1 at that stage)
 */
#[allow(unused_assignments)] /* C-faithful dead stores (f1nprv bookkeeping) */
fn KINLineSearch(
    kin_mem: &mut KINMem,
    fnormp: &mut f64,
    f1normp: &mut f64,
    maxStepTaken: &mut bool,
) -> i32 {
    /* Initializations */

    let mut nbktrk_l: i64 = 0; /* local backtracking counter */
    let mut ratio = ONE; /* step change ratio          */
    let alpha = POINT0001;
    let beta = POINT9;

    let mut firstBacktrack = SUNTRUE;
    *maxStepTaken = SUNFALSE;

    let mut rlprev = ZERO;
    let mut f1nprv = ZERO;

    /* Compute length of Newton step */

    let mut pnorm = N_VWL2Norm(&kin_mem.kin_pp, &kin_mem.kin_uscale);
    let mut rlmax = kin_mem.kin_mxnewtstep / pnorm;
    kin_mem.kin_stepl = pnorm;

    /* If the full Newton step is too large, set it to the maximum allowable value */

    if pnorm > kin_mem.kin_mxnewtstep {
        ratio = kin_mem.kin_mxnewtstep / pnorm;
        kin_mem.kin_pp.scale_inplace(ratio); /* N_VScale(ratio, pp, pp) */
        pnorm = kin_mem.kin_mxnewtstep;
        rlmax = ONE;
        kin_mem.kin_stepl = pnorm;
    }

    /* If constraint checking is activated, check and correct violations */

    kin_mem.kin_stepmul = ONE;

    if kin_mem.kin_constraintsSet {
        let retval = KINConstraint(kin_mem);
        if retval == CONSTR_VIOLATED {
            /* Apply stepmul set in KINConstraint */
            let stepmul = kin_mem.kin_stepmul;
            kin_mem.kin_pp.scale_inplace(stepmul);
            ratio *= stepmul;
            pnorm *= stepmul;
            rlmax = ONE;
            kin_mem.kin_stepl = pnorm;
            /* (KINPrintInfo PRNT_PNORM1: SUNDIALS_LOGGING_LEVEL >= INFO only) */
            if pnorm <= kin_mem.kin_scsteptol {
                N_VLinearSum(ONE, &kin_mem.kin_uu, ONE, &kin_mem.kin_pp,
                             &mut kin_mem.kin_unew);
                return STEP_TOO_SMALL;
            }
        }
    }

    /* Attempt (at most MAX_RECVR times) to evaluate function at the new iterate */

    let mut fOK = SUNFALSE;
    let func = kin_mem.kin_func.unwrap();

    for _ircvr in 1..=MAX_RECVR {
        /* compute the iterate unew = uu + pp */
        N_VLinearSum(ONE, &kin_mem.kin_uu, ONE, &kin_mem.kin_pp, &mut kin_mem.kin_unew);

        /* evaluate func(unew) and its norm, and return */
        let retval =
            func(&kin_mem.kin_unew, &mut kin_mem.kin_fval, &mut kin_mem.kin_user_data);
        kin_mem.kin_nfe += 1;

        /* if func was successful, accept pp */
        if retval == 0 {
            fOK = SUNTRUE;
            break;
        }
        /* if func failed unrecoverably, give up */
        else if retval < 0 {
            return KIN_SYSFUNC_FAIL;
        }

        /* func failed recoverably; cut step in half and try again */
        kin_mem.kin_pp.scale_inplace(HALF);
        ratio *= HALF;
        pnorm *= HALF;
        rlmax = ONE;
        kin_mem.kin_stepl = pnorm;
    }

    /* If func() failed recoverably MAX_RECVR times, give up */

    if !fOK {
        return KIN_REPTD_SYSFUNC_ERR;
    }

    /* Evaluate function norms */

    *fnormp = N_VWL2Norm(&kin_mem.kin_fval, &kin_mem.kin_fscale);
    *f1normp = HALF * (*fnormp) * (*fnormp);

    /* Estimate the line search value rl (lambda) to satisfy both ALPHA and BETA conditions */

    let slpi = kin_mem.kin_sFdotJp * ratio;
    let rlength = {
        let KINMem { kin_vtemp1, kin_vtemp2, kin_uscale, kin_pp, kin_uu, .. } = kin_mem;
        KINScSNorm(kin_vtemp1, kin_vtemp2, kin_uscale, kin_pp, kin_uu)
    };
    let rlmin = kin_mem.kin_scsteptol / rlength;
    let mut rl = ONE;

    /* (KINPrintInfo PRNT_LAM: SUNDIALS_LOGGING_LEVEL >= INFO only) */

    /* Loop until the ALPHA condition is satisfied. Terminate if rl becomes too small */

    let mut alpha_cond;
    loop {
        /* Evaluate test quantity */

        alpha_cond = kin_mem.kin_f1norm + (alpha * slpi * rl);

        /* (KINPrintInfo PRNT_ALPHA: SUNDIALS_LOGGING_LEVEL >= INFO only) */

        /* If ALPHA condition is satisfied, break out from loop */

        if *f1normp <= alpha_cond {
            break;
        }

        /* Backtracking. Use quadratic fit the first time and cubic fit afterwards. */

        let mut rltmp;
        if firstBacktrack {
            rltmp = -slpi / (TWO * ((*f1normp) - kin_mem.kin_f1norm - slpi));
            firstBacktrack = SUNFALSE;
        } else {
            let mut tmp1 = (*f1normp) - kin_mem.kin_f1norm - (rl * slpi);
            let tmp2 = f1nprv - kin_mem.kin_f1norm - (rlprev * slpi);
            let mut rl_a = ((ONE / (rl * rl)) * tmp1) - ((ONE / (rlprev * rlprev)) * tmp2);
            let mut rl_b =
                ((-rlprev / (rl * rl)) * tmp1) + ((rl / (rlprev * rlprev)) * tmp2);
            tmp1 = ONE / (rl - rlprev);
            rl_a *= tmp1;
            rl_b *= tmp1;
            let disc = (rl_b * rl_b) - (THREE * rl_a * slpi);

            if SUNRabs(rl_a) < kin_mem.kin_uround {
                /* cubic is actually just a quadratic (rl_a ~ 0) */
                rltmp = -slpi / (TWO * rl_b);
            } else {
                /* real cubic */
                rltmp = (-rl_b + SUNRsqrt(disc)) / (THREE * rl_a);
            }
        }
        if rltmp > HALF * rl {
            rltmp = HALF * rl;
        }

        /* Set new rl (do not allow a reduction by a factor larger than 10) */

        rlprev = rl;
        f1nprv = *f1normp;
        let pt1trl = POINT1 * rl;
        rl = SUNMAX(pt1trl, rltmp);
        nbktrk_l += 1;

        /* Update unew and re-evaluate function */

        N_VLinearSum(ONE, &kin_mem.kin_uu, rl, &kin_mem.kin_pp, &mut kin_mem.kin_unew);

        let retval =
            func(&kin_mem.kin_unew, &mut kin_mem.kin_fval, &mut kin_mem.kin_user_data);
        kin_mem.kin_nfe += 1;
        if retval != 0 {
            return KIN_SYSFUNC_FAIL;
        }

        *fnormp = N_VWL2Norm(&kin_mem.kin_fval, &kin_mem.kin_fscale);
        *f1normp = HALF * (*fnormp) * (*fnormp);

        /* Check if rl (lambda) is too small */

        if rl < rlmin {
            /* unew sufficiently distinct from uu cannot be found.
               copy uu into unew (step remains unchanged) and
               return STEP_TOO_SMALL */
            N_VScale(ONE, &kin_mem.kin_uu, &mut kin_mem.kin_unew);
            return STEP_TOO_SMALL;
        }
    } /* end ALPHA condition loop */

    /* ALPHA condition is satisfied. Now check the BETA condition */

    let mut beta_cond = kin_mem.kin_f1norm + (beta * slpi * rl);

    if *f1normp < beta_cond {
        /* BETA condition not satisfied */

        if rl == ONE && pnorm < kin_mem.kin_mxnewtstep {
            loop {
                rlprev = rl;
                f1nprv = *f1normp;
                rl = SUNMIN(TWO * rl, rlmax);
                nbktrk_l += 1;

                N_VLinearSum(ONE, &kin_mem.kin_uu, rl, &kin_mem.kin_pp,
                             &mut kin_mem.kin_unew);
                let retval = func(&kin_mem.kin_unew, &mut kin_mem.kin_fval,
                                  &mut kin_mem.kin_user_data);
                kin_mem.kin_nfe += 1;
                if retval != 0 {
                    return KIN_SYSFUNC_FAIL;
                }
                *fnormp = N_VWL2Norm(&kin_mem.kin_fval, &kin_mem.kin_fscale);
                *f1normp = HALF * (*fnormp) * (*fnormp);

                alpha_cond = kin_mem.kin_f1norm + (alpha * slpi * rl);
                beta_cond = kin_mem.kin_f1norm + (beta * slpi * rl);

                /* (KINPrintInfo PRNT_BETA: SUNDIALS_LOGGING_LEVEL >= INFO only) */

                if !((*f1normp <= alpha_cond) && (*f1normp < beta_cond) && (rl < rlmax)) {
                    break;
                }
            }
        } /* end if (rl == ONE) block */

        if (rl < ONE) || ((rl > ONE) && (*f1normp > alpha_cond)) {
            let mut rllo = SUNMIN(rl, rlprev);
            let mut rldiff = SUNRabs(rlprev - rl);

            loop {
                let rlinc = HALF * rldiff;
                rl = rllo + rlinc;
                nbktrk_l += 1;

                N_VLinearSum(ONE, &kin_mem.kin_uu, rl, &kin_mem.kin_pp,
                             &mut kin_mem.kin_unew);
                let retval = func(&kin_mem.kin_unew, &mut kin_mem.kin_fval,
                                  &mut kin_mem.kin_user_data);
                kin_mem.kin_nfe += 1;
                if retval != 0 {
                    return KIN_SYSFUNC_FAIL;
                }
                *fnormp = N_VWL2Norm(&kin_mem.kin_fval, &kin_mem.kin_fscale);
                *f1normp = HALF * (*fnormp) * (*fnormp);

                alpha_cond = kin_mem.kin_f1norm + (alpha * slpi * rl);
                beta_cond = kin_mem.kin_f1norm + (beta * slpi * rl);

                /* (KINPrintInfo PRNT_ALPHABETA: SUNDIALS_LOGGING_LEVEL >= INFO only) */

                if *f1normp > alpha_cond {
                    rldiff = rlinc;
                } else if *f1normp < beta_cond {
                    rllo = rl;
                    rldiff -= rlinc;
                }

                if !((*f1normp > alpha_cond)
                    || ((*f1normp < beta_cond) && (rldiff >= rlmin)))
                {
                    break;
                }
            }

            if (*f1normp < beta_cond) || ((rldiff < rlmin) && (*f1normp > alpha_cond)) {
                /* beta condition could not be satisfied or rldiff too small
                   and alpha_cond not satisfied, so set unew to last u value
                   that satisfied the alpha condition and continue */

                N_VLinearSum(ONE, &kin_mem.kin_uu, rllo, &kin_mem.kin_pp,
                             &mut kin_mem.kin_unew);
                let retval = func(&kin_mem.kin_unew, &mut kin_mem.kin_fval,
                                  &mut kin_mem.kin_user_data);
                kin_mem.kin_nfe += 1;
                if retval != 0 {
                    return KIN_SYSFUNC_FAIL;
                }
                *fnormp = N_VWL2Norm(&kin_mem.kin_fval, &kin_mem.kin_fscale);
                *f1normp = HALF * (*fnormp) * (*fnormp);

                /* increment beta-condition failures counter */

                kin_mem.kin_nbcf += 1;
            }
        } /* end of if (rl < ONE) block */
    } /* end of if (f1normp < beta_cond) block */

    /* Update number of backtracking operations */

    kin_mem.kin_nbktrk += nbktrk_l;

    /* (KINPrintInfo PRNT_ADJ: SUNDIALS_LOGGING_LEVEL >= INFO only) */

    /* scale sFdotJp and sJpnorm by rl * ratio for later use in KINForcingTerm */

    kin_mem.kin_sFdotJp = kin_mem.kin_sFdotJp * rl * ratio;
    kin_mem.kin_sJpnorm = kin_mem.kin_sJpnorm * rl * ratio;

    if (rl * pnorm) > (POINT99 * kin_mem.kin_mxnewtstep) {
        *maxStepTaken = SUNTRUE;
    }

    KIN_SUCCESS
}

/*
 * Function : KINConstraint
 *
 * This routine checks if the proposed solution vector uu + pp
 * violates any constraints. If a constraint is violated, then the
 * scalar stepmul is determined such that uu + stepmul * pp does
 * not violate any constraints.
 *
 * Note: This routine is called by the functions
 *       KINLineSearch and KINFullNewton.
 */
fn KINConstraint(kin_mem: &mut KINMem) -> i32 {
    N_VLinearSum(ONE, &kin_mem.kin_uu, ONE, &kin_mem.kin_pp, &mut kin_mem.kin_vtemp1);

    /* if vtemp1[i] violates constraint[i] then vtemp2[i] = 1
       else vtemp2[i] = 0 (vtemp2 is the mask vector) */

    let constr_ok = {
        let KINMem { kin_constraints, kin_vtemp1, kin_vtemp2, .. } = kin_mem;
        N_VConstrMask(kin_constraints, kin_vtemp1, kin_vtemp2)
    };
    if constr_ok {
        return KIN_SUCCESS;
    }

    /* vtemp1[i] = SUNRabs(pp[i]) */

    {
        let KINMem { kin_pp, kin_vtemp1, .. } = kin_mem;
        N_VAbs(kin_pp, kin_vtemp1);
    }

    /* consider vtemp1[i] only if vtemp2[i] = 1 (constraint violated) */
    /* N_VProd(vtemp2, vtemp1, vtemp1) — output aliases the second operand */

    {
        let KINMem { kin_vtemp1, kin_vtemp2, .. } = kin_mem;
        kin_vtemp1.prod_with(kin_vtemp2);
    }

    {
        let KINMem { kin_uu, kin_vtemp2, .. } = kin_mem;
        N_VAbs(kin_uu, kin_vtemp2);
    }
    kin_mem.kin_stepmul =
        POINT9 * N_VMinQuotient(&kin_mem.kin_vtemp2, &kin_mem.kin_vtemp1);

    CONSTR_VIOLATED
}

/*
 * -----------------------------------------------------------------
 * Stopping tests
 * -----------------------------------------------------------------
 */

/*
 * KINStop
 *
 * This routine checks the current iterate unew to see if the
 * system func(unew) = 0 is satisfied by a variety of tests.
 *
 * strategy is one of KIN_NONE or KIN_LINESEARCH
 * sflag    is one of KIN_SUCCESS, STEP_TOO_SMALL
 */
fn KINStop(kin_mem: &mut KINMem, maxStepTaken: bool, sflag: i32) -> i32 {
    /* Check for too small a step */

    if sflag == STEP_TOO_SMALL {
        if kin_has_lsetup(kin_mem) && !kin_mem.kin_jacCurrent {
            /* If the Jacobian is out of date, update it and retry */
            kin_mem.kin_sthrsh = TWO;
            return RETRY_ITERATION;
        } else {
            /* Give up */
            if kin_mem.kin_globalstrategy == KIN_NONE {
                return KIN_STEP_LT_STPTOL;
            } else {
                return KIN_LINESEARCH_NONCONV;
            }
        }
    }

    /* Check tolerance on scaled function norm at the current iterate */

    let fmax = {
        let KINMem { kin_vtemp1, kin_fval, kin_fscale, .. } = kin_mem;
        KINScFNorm(kin_vtemp1, kin_fval, kin_fscale)
    };

    /* (KINPrintInfo PRNT_FMAX: SUNDIALS_LOGGING_LEVEL >= INFO only) */

    if fmax <= kin_mem.kin_fnormtol {
        return KIN_SUCCESS;
    }

    /* Check if the scaled distance between the last two steps is too small */
    /* NOTE: pp used as work space to store this distance */

    /* delta = kin_pp */
    N_VLinearSum(ONE, &kin_mem.kin_unew, -ONE, &kin_mem.kin_uu, &mut kin_mem.kin_pp);
    let rlength = {
        let KINMem { kin_vtemp1, kin_vtemp2, kin_uscale, kin_pp, kin_unew, .. } = kin_mem;
        KINScSNorm(kin_vtemp1, kin_vtemp2, kin_uscale, kin_pp, kin_unew)
    };

    if rlength <= kin_mem.kin_scsteptol {
        if kin_has_lsetup(kin_mem) && !kin_mem.kin_jacCurrent {
            /* If the Jacobian is out of date, update it and retry */
            kin_mem.kin_sthrsh = TWO;
            return CONTINUE_ITERATIONS;
        } else {
            /* give up */
            return KIN_STEP_LT_STPTOL;
        }
    }

    /* Check if the maximum number of iterations is reached */

    if kin_mem.kin_nni >= kin_mem.kin_mxiter {
        return KIN_MAXITER_REACHED;
    }

    /* Check for consecutive number of steps taken of size mxnewtstep
       and if not maxStepTaken, then set ncscmx to 0 */

    if maxStepTaken {
        kin_mem.kin_ncscmx += 1;
    } else {
        kin_mem.kin_ncscmx = 0;
    }

    if kin_mem.kin_ncscmx == 5 {
        return KIN_MXNEWT_5X_EXCEEDED;
    }

    /* Proceed according to the type of linear solver used */

    if kin_mem.kin_inexact_ls {
        /* We're doing inexact Newton.
           Load threshold for reevaluating the Jacobian. */

        kin_mem.kin_sthrsh = rlength;
    } else if !kin_mem.kin_noResMon {
        /* We're doing modified Newton and the user did not disable residual monitoring.
           Check if it is time to monitor residual. */

        if (kin_mem.kin_nni - kin_mem.kin_nnilset_sub) >= kin_mem.kin_msbset_sub {
            /* Residual monitoring needed */

            kin_mem.kin_nnilset_sub = kin_mem.kin_nni;

            /* If indicated, estimate new OMEGA value */
            if kin_mem.kin_eval_omega {
                let omexp = SUNMAX(ZERO, (kin_mem.kin_fnorm / kin_mem.kin_fnormtol) - ONE);
                kin_mem.kin_omega = if omexp > TWELVE {
                    kin_mem.kin_omega_max
                } else {
                    SUNMIN(kin_mem.kin_omega_min * SUNRexp(omexp), kin_mem.kin_omega_max)
                };
            }
            /* Check if making satisfactory progress */

            if kin_mem.kin_fnorm > kin_mem.kin_omega * kin_mem.kin_fnorm_sub {
                /* Insufficient progress */
                if kin_has_lsetup(kin_mem) && !kin_mem.kin_jacCurrent {
                    /* If the Jacobian is out of date, update it and retry */
                    kin_mem.kin_sthrsh = TWO;
                    return CONTINUE_ITERATIONS;
                }
                /* Otherwise, we cannot do anything, so just return. */
            } else {
                /* Sufficient progress */
                kin_mem.kin_fnorm_sub = kin_mem.kin_fnorm;
                kin_mem.kin_sthrsh = ONE;
            }
        } else {
            /* Residual monitoring not needed */

            /* Reset sthrsh */
            if kin_mem.kin_retry_nni || kin_mem.kin_update_fnorm_sub {
                kin_mem.kin_fnorm_sub = kin_mem.kin_fnorm;
            }
            if kin_mem.kin_update_fnorm_sub {
                kin_mem.kin_update_fnorm_sub = SUNFALSE;
            }
            kin_mem.kin_sthrsh = ONE;
        }
    }

    /* if made it to here, then the iteration process is not finished
       so return CONTINUE_ITERATIONS flag */

    CONTINUE_ITERATIONS
}

/*
 * KINForcingTerm
 *
 * This routine computes eta, the scaling factor in the linear
 * convergence stopping tolerance eps when choice #1 or choice #2
 * forcing terms are used. Eta is computed here for all but the
 * first iterative step, which is set to the default in routine
 * KINSolInit.
 *
 * This routine was written by Homer Walker of Utah State
 * University with subsequent modifications by Allan Taylor @ LLNL.
 *
 * It is based on the concepts of the paper 'Choosing the forcing
 * terms in an inexact Newton method', SIAM J Sci Comput, 17
 * (1996), pp 16 - 32, or Utah State University Research Report
 * 6/94/75 of the same title.
 */
fn KINForcingTerm(kin_mem: &mut KINMem, fnormp: f64) {
    let eta_max = POINT9;
    let eta_min = POINT0001;
    let mut eta_safe = HALF;

    /* choice #1 forcing term */

    if kin_mem.kin_etaflag == KIN_ETACHOICE1 {
        /* compute the norm of f + Jp , scaled L2 norm */

        let linmodel_norm = SUNRsqrt(
            (kin_mem.kin_fnorm * kin_mem.kin_fnorm) + (TWO * kin_mem.kin_sFdotJp)
                + (kin_mem.kin_sJpnorm * kin_mem.kin_sJpnorm),
        );

        /* form the safeguarded for choice #1 */

        eta_safe = SUNRpowerR(kin_mem.kin_eta, kin_mem.kin_eta_alpha);
        kin_mem.kin_eta = SUNRabs(fnormp - linmodel_norm) / kin_mem.kin_fnorm;
    }

    /* choice #2 forcing term */

    if kin_mem.kin_etaflag == KIN_ETACHOICE2 {
        eta_safe = kin_mem.kin_eta_gamma
            * SUNRpowerR(kin_mem.kin_eta, kin_mem.kin_eta_alpha);

        kin_mem.kin_eta = kin_mem.kin_eta_gamma
            * SUNRpowerR(fnormp / kin_mem.kin_fnorm, kin_mem.kin_eta_alpha);
    }

    /* apply safeguards */

    if eta_safe < POINT1 {
        eta_safe = ZERO;
    }
    kin_mem.kin_eta = SUNMAX(kin_mem.kin_eta, eta_safe);
    kin_mem.kin_eta = SUNMAX(kin_mem.kin_eta, eta_min);
    kin_mem.kin_eta = SUNMIN(kin_mem.kin_eta, eta_max);
}

/*
 * -----------------------------------------------------------------
 * Norm functions
 * -----------------------------------------------------------------
 */

/*
 * Function : KINScFNorm
 *
 * This routine computes the max norm for scaled vectors. The
 * scaling vector is scale, and the vector of which the norm is to
 * be determined is vv. The returned value, fnormval, is the
 * resulting scaled vector norm. This routine uses N_Vector
 * functions from the vector module.
 *
 * (C signature takes kin_mem; only kin_vtemp1 is used, so the
 * workspace vector is passed explicitly to permit disjoint field
 * borrows. The one call site where v aliases vtemp1 — KINFP — is
 * inlined there.)
 */
fn KINScFNorm(vtemp1: &mut NVector, v: &NVector, scale: &NVector) -> f64 {
    N_VProd(scale, v, vtemp1);
    N_VMaxNorm(vtemp1)
}

/*
 * Function : KINScSNorm
 *
 * This routine computes the max norm of the scaled steplength, ss.
 * Here ucur is the current step and usc is the u scale factor.
 *
 * (C signature takes kin_mem; kin_vtemp1/kin_vtemp2/kin_uscale are
 * passed explicitly to permit disjoint field borrows.)
 */
fn KINScSNorm(
    vtemp1: &mut NVector,
    vtemp2: &mut NVector,
    uscale: &NVector,
    v: &NVector,
    u: &NVector,
) -> f64 {
    N_VInv(uscale, vtemp1);
    N_VAbs(u, vtemp2);
    /* N_VLinearSum(ONE, vtemp1, ONE, vtemp2, vtemp1) — output aliases
       the first operand */
    vtemp1.linear_sum_with(ONE, ONE, vtemp2);
    /* N_VDiv(v, vtemp1, vtemp1) — output aliases the denominator; the
       serial kernel computes zd[i] = xd[i]/yd[i] */
    for i in 0..vtemp1.data.len() {
        vtemp1.data[i] = v.data[i] / vtemp1.data[i];
    }

    N_VMaxNorm(vtemp1)
}

/*
 * =================================================================
 * KINSOL Verbose output functions
 * =================================================================
 */

/*
 * KINPrintInfo
 *
 * KINPrintInfo is a high level error handling function.
 * Based on the value info_code, it composes the info message and
 * passes it to the info handler function.
 *
 * Adaptations: the printf varargs are rendered at each call site
 * (sundials_utils::fmt_* helpers), so `msg` arrives fully formatted;
 * for PRNT_RETVAL the numeric value (C: va_arg) is recovered from
 * the rendered INFO_RETVAL text ("Return value: %d"). The composed
 * message is queued in C on the SUNContext logger at
 * SUN_LOGLEVEL_INFO (SUNLogger_QueueMsg); the reference build
 * compiles with SUNDIALS_LOGGING_LEVEL < INFO, so the message is
 * discarded — the workspace SUNContext carries no logger
 * (ARCHITECTURE.md Addendum A) and this port discards it likewise.
 */
pub fn KINPrintInfo(_kin_mem: &KINMem, info_code: i32, _module: &str, _fname: &str,
                    msg: &str) {
    let composed: String;
    let msg_final: &str = if info_code == PRNT_RETVAL {
        /* If info_code = PRNT_RETVAL, decode the numeric value */
        let ret: i32 = msg
            .rsplit(' ')
            .next()
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);
        let retstr = match ret {
            KIN_SUCCESS => "KIN_SUCCESS",
            KIN_SYSFUNC_FAIL => "KIN_SYSFUNC_FAIL",
            KIN_REPTD_SYSFUNC_ERR => "KIN_REPTD_SYSFUNC_ERR",
            KIN_STEP_LT_STPTOL => "KIN_STEP_LT_STPTOL",
            KIN_LINESEARCH_NONCONV => "KIN_LINESEARCH_NONCONV",
            KIN_LINESEARCH_BCFAIL => "KIN_LINESEARCH_BCFAIL",
            KIN_MAXITER_REACHED => "KIN_MAXITER_REACHED",
            KIN_MXNEWT_5X_EXCEEDED => "KIN_MXNEWT_5X_EXCEEDED",
            KIN_LINSOLV_NO_RECOVERY => "KIN_LINSOLV_NO_RECOVERY",
            KIN_LSETUP_FAIL => "KIN_PRECONDSET_FAILURE",
            KIN_LSOLVE_FAIL => "KIN_PRECONDSOLVE_FAILURE",
            _ => "",
        };
        /* Compose the message: sprintf(msg, "%s (%s)", msg1, retstr) */
        composed = format!("{msg} ({retstr})");
        &composed
    } else {
        msg
    };

    /* SUNLogger_QueueMsg(KIN_LOGGER, SUN_LOGLEVEL_INFO, fname, "KINSOL",
       "%s", msg): no logger in the default build — message discarded */
    let _ = msg_final;
}

/*
 * =================================================================
 * KINSOL Error Handling functions
 * =================================================================
 *
 * KINProcessError is defined in kinsol_impl.rs (donor convention:
 * the C varargs/error-handler-stack plumbing routes to stderr).
 * The KINInfoHandler prototype in kinsol_impl.h has no definition
 * in the 7.7.0 sources and is not ported.
 */

/*
 * =======================================================================
 * Picard and fixed point solvers
 * =======================================================================
 */

/*
 * KINPicardAA
 *
 * This routine is the main driver for the Picard iteration with
 * accelerated fixed point.
 */
fn KINPicardAA(kin_mem: &mut KINMem) -> i32 {
    let mut ret; /* iteration status            */
    let mut epsmin = ZERO;

    /* initialize iteration count */
    kin_mem.kin_nni = 0;

    /* if eps is to be bounded from below, set the bound */
    if kin_mem.kin_inexact_ls && !kin_mem.kin_noMinEps {
        epsmin = POINT01 * kin_mem.kin_fnormtol;
    }

    ret = CONTINUE_ITERATIONS;

    while ret == CONTINUE_ITERATIONS {
        /* update iteration count */
        kin_mem.kin_nni += 1;

        /* Update the forcing term for the inexact linear solves */
        if kin_mem.kin_inexact_ls {
            kin_mem.kin_eps = (kin_mem.kin_eta + kin_mem.kin_uround) * kin_mem.kin_fnorm;
            if !kin_mem.kin_noMinEps {
                kin_mem.kin_eps = SUNMAX(epsmin, kin_mem.kin_eps);
            }
        }

        /* evaluate g = uu - L^{-1}func(uu) and return if failed.
           For Picard, assume that the fval vector has been filled
           with an eval of the nonlinear residual prior to this call. */
        let retval = KINPicardFcnEval(kin_mem);

        if retval < 0 {
            ret = KIN_SYSFUNC_FAIL;
            break;
        }

        /* compute new solution */
        if kin_mem.kin_m_aa == 0 || kin_mem.kin_nni - 1 < kin_mem.kin_delay_aa {
            if kin_mem.kin_damping || kin_mem.kin_damping_fn.is_some() {
                if let Some(damping_fn) = kin_mem.kin_damping_fn {
                    let retval = {
                        let KINMem { kin_nni, kin_uu, kin_gval, kin_user_data, kin_beta, .. } =
                            kin_mem;
                        damping_fn(*kin_nni, kin_uu, kin_gval, &[], 0, kin_user_data,
                                   kin_beta)
                    };
                    if retval != 0 {
                        KINProcessError(Some(kin_mem), KIN_DAMPING_FN_ERR, line!(),
                                        "KINPicardAA", file!(),
                                        "The damping function failed.");
                        ret = KIN_DAMPING_FN_ERR;
                        break;
                    }
                    if kin_mem.kin_beta <= ZERO || kin_mem.kin_beta > ONE {
                        KINProcessError(Some(kin_mem), KIN_DAMPING_FN_ERR, line!(),
                                        "KINPicardAA", file!(),
                                        "The damping parameter is outside of the range (0, 1].");
                        ret = KIN_DAMPING_FN_ERR;
                        break;
                    }
                }

                /* damped fixed point */
                N_VLinearSum(ONE - kin_mem.kin_beta, &kin_mem.kin_uu, kin_mem.kin_beta,
                             &kin_mem.kin_gval, &mut kin_mem.kin_unew);
            } else {
                /* standard fixed point */
                N_VScale(ONE, &kin_mem.kin_gval, &mut kin_mem.kin_unew);
            }
        } else {
            /* compute iteration count for Anderson acceleration */
            let iter_aa = if kin_mem.kin_delay_aa > 0 {
                kin_mem.kin_nni - 1 - kin_mem.kin_delay_aa
            } else {
                kin_mem.kin_nni - 1
            };

            /* C: AndersonAcc(kin_mem, gval, delta(=vtemp1), unew, uu,
               iter_aa, R_aa, gamma_aa) — the vector/array arguments are
               KINMem fields, take()n out around the call */
            let gval = std::mem::take(&mut kin_mem.kin_gval);
            let mut fv = std::mem::take(&mut kin_mem.kin_vtemp1); /* delta (temp) */
            let mut x = std::mem::take(&mut kin_mem.kin_unew);
            let xold = std::mem::take(&mut kin_mem.kin_uu);
            let mut r = std::mem::take(&mut kin_mem.kin_R_aa);
            let mut gamma = std::mem::take(&mut kin_mem.kin_gamma_aa);
            let retval = AndersonAcc(kin_mem, /* kinsol memory            */
                                     &gval,   /* G(u_cur)       in        */
                                     &mut fv, /* F(u_cur)       in (temp) */
                                     &mut x,  /* u_new output   out       */
                                     &xold,   /* u_cur input    in        */
                                     iter_aa, /* AA iteration   in        */
                                     &mut r,  /* R matrix       in/out    */
                                     &mut gamma); /* gamma vector in (temp) */
            kin_mem.kin_gval = gval;
            kin_mem.kin_vtemp1 = fv;
            kin_mem.kin_unew = x;
            kin_mem.kin_uu = xold;
            kin_mem.kin_R_aa = r;
            kin_mem.kin_gamma_aa = gamma;
            if retval != 0 {
                ret = retval;
                break;
            }
        }

        /* Fill the Newton residual based on the new solution iterate */
        let func = kin_mem.kin_func.unwrap();
        let retval =
            func(&kin_mem.kin_unew, &mut kin_mem.kin_fval, &mut kin_mem.kin_user_data);
        kin_mem.kin_nfe += 1;

        if retval < 0 {
            ret = KIN_SYSFUNC_FAIL;
            break;
        }

        /* Measure || F(x) ||_max */
        kin_mem.kin_fnorm = {
            let KINMem { kin_vtemp1, kin_fval, kin_fscale, .. } = kin_mem;
            KINScFNorm(kin_vtemp1, kin_fval, kin_fscale)
        };

        /* (KINPrintInfo PRNT_FMAX / PRNT_NNI: SUNDIALS_LOGGING_LEVEL >= INFO only) */

        /* Check if the maximum number of iterations is reached */
        if kin_mem.kin_nni >= kin_mem.kin_mxiter {
            ret = KIN_MAXITER_REACHED;
        }
        if kin_mem.kin_fnorm <= kin_mem.kin_fnormtol {
            ret = KIN_SUCCESS;
        }

        /* Update the solution. Always return the newest iteration. Note this is
           also consistent with last function evaluation. */
        N_VScale(ONE, &kin_mem.kin_unew, &mut kin_mem.kin_uu);

        if ret == CONTINUE_ITERATIONS && kin_mem.kin_callForcingTerm {
            /* evaluate eta by calling the forcing term routine */
            let fnormp = N_VWL2Norm(&kin_mem.kin_fval, &kin_mem.kin_fscale);
            KINForcingTerm(kin_mem, fnormp);
        }
    } /* end of loop; return */

    /* (KINPrintInfo PRNT_RETVAL: SUNDIALS_LOGGING_LEVEL >= INFO only) */

    ret
}

/*
 * KINPicardFcnEval
 *
 * This routine evaluates the Picard fixed point function
 * using the linear solver, gval = u - L^{-1}F(u).
 * The function assumes the user has defined L either through
 * a user-supplied matvec if using a SPILS solver or through
 * a supplied matrix if using a dense solver.  This assumption is
 * tested by a check on the strategy and the requisite functionality
 * within the linear solve routines.
 *
 * This routine fills gval = uu - L^{-1}F(uu) given uu and fval = F(uu).
 *
 * (C signature: KINPicardFcnEval(kin_mem, gval, uval, fval1), always
 * called with gval = kin_gval, uval = kin_uu, fval1 = kin_fval; the
 * indirection is collapsed onto the KINMem fields, donor cvode_nls
 * convention. In C the lsolve rhs `bb` *aliases* kin_fval; here a
 * stand-in copy is passed and written back after the call, so the
 * KINMem state kinLsSolve observes — and the value kin_fval holds
 * afterwards — match the C aliasing.)
 */
fn KINPicardFcnEval(kin_mem: &mut KINMem) -> i32 {
    if (kin_mem.kin_nni - kin_mem.kin_nnilset) >= kin_mem.kin_msbset {
        kin_mem.kin_sthrsh = TWO;
        kin_mem.kin_update_fnorm_sub = SUNTRUE;
    }

    loop {
        kin_mem.kin_jacCurrent = SUNFALSE;

        if kin_mem.kin_sthrsh > ONEPT5 && kin_has_lsetup(kin_mem) {
            let retval = kin_lsetup_dispatch(kin_mem);
            kin_mem.kin_jacCurrent = SUNTRUE;
            kin_mem.kin_nnilset = kin_mem.kin_nni;
            kin_mem.kin_nnilset_sub = kin_mem.kin_nni;
            if retval != 0 {
                return KIN_LSETUP_FAIL;
            }
        }

        /* call the generic 'lsolve' routine to solve the system Lx = -fval
           Note that we are using gval to hold x. */
        /* N_VScale(-ONE, fval1, fval1) — output aliases the input */
        kin_mem.kin_fval.scale_inplace(-ONE);

        let mut xx = std::mem::take(&mut kin_mem.kin_gval);
        let mut bb = kin_mem.kin_fval.clone(); /* bb aliases kin_fval in C */
        let retval = kin_lsolve_dispatch(kin_mem, &mut xx, &mut bb);
        kin_mem.kin_fval = bb; /* preserve the C alias side effects */
        kin_mem.kin_gval = xx;

        if retval == 0 {
            /* Update gval = uval + gval since gval = -L^{-1}F(uu)  */
            /* N_VLinearSum(ONE, uval, ONE, gval, gval) — output aliases
               the second operand */
            {
                let KINMem { kin_uu, kin_gval, .. } = kin_mem;
                kin_gval.linear_sum_with(ONE, ONE, kin_uu);
            }
            return KIN_SUCCESS;
        } else if retval < 0 {
            return KIN_LSOLVE_FAIL;
        } else if !kin_has_lsetup(kin_mem) || kin_mem.kin_jacCurrent {
            return KIN_LINSOLV_NO_RECOVERY;
        }

        /* loop back only if the linear solver setup is in use
           and matrix information is not current */

        kin_mem.kin_sthrsh = TWO;
    }
}

/*
 * KINFP
 *
 * This routine is the main driver for the fixed point iteration with
 * Anderson Acceleration.
 */
#[allow(unused_assignments)] /* C-faithful dead store (tolfac init) */
fn KINFP(kin_mem: &mut KINMem) -> i32 {
    let mut ret; /* iteration status            */
    let mut tolfac; /* tolerance adjustment factor */

    ret = CONTINUE_ITERATIONS;
    tolfac = ONE;

    /* initialize iteration count */
    kin_mem.kin_nni = 0;

    while ret == CONTINUE_ITERATIONS {
        /* update iteration count */
        kin_mem.kin_nni += 1;

        /* evaluate func(uu) and return if failed */
        let func = kin_mem.kin_func.unwrap();
        let retval =
            func(&kin_mem.kin_uu, &mut kin_mem.kin_fval, &mut kin_mem.kin_user_data);
        kin_mem.kin_nfe += 1;

        if retval < 0 {
            ret = KIN_SYSFUNC_FAIL;
            break;
        }

        /* compute new solution */
        if kin_mem.kin_m_aa == 0 || kin_mem.kin_nni - 1 < kin_mem.kin_delay_aa {
            if kin_mem.kin_damping || kin_mem.kin_damping_fn.is_some() {
                if let Some(damping_fn) = kin_mem.kin_damping_fn {
                    let retval = {
                        let KINMem { kin_nni, kin_uu, kin_fval, kin_user_data, kin_beta, .. } =
                            kin_mem;
                        damping_fn(*kin_nni, kin_uu, kin_fval, &[], 0, kin_user_data,
                                   kin_beta)
                    };
                    if retval != 0 {
                        KINProcessError(Some(kin_mem), KIN_DAMPING_FN_ERR, line!(), "KINFP",
                                        file!(), "The damping function failed.");
                        ret = KIN_DAMPING_FN_ERR;
                        break;
                    }
                    if kin_mem.kin_beta <= ZERO || kin_mem.kin_beta > ONE {
                        KINProcessError(Some(kin_mem), KIN_DAMPING_FN_ERR, line!(), "KINFP",
                                        file!(),
                                        "The damping parameter is outside of the range (0, 1].");
                        ret = KIN_DAMPING_FN_ERR;
                        break;
                    }
                }

                /* damped fixed point */
                N_VLinearSum(ONE - kin_mem.kin_beta, &kin_mem.kin_uu, kin_mem.kin_beta,
                             &kin_mem.kin_fval, &mut kin_mem.kin_unew);

                /* tolerance adjustment */
                tolfac = kin_mem.kin_beta;
            } else {
                /* standard fixed point */
                N_VScale(ONE, &kin_mem.kin_fval, &mut kin_mem.kin_unew);

                /* tolerance adjustment */
                tolfac = ONE;
            }
        } else {
            /* compute iteration count for Anderson acceleration */
            let iter_aa = if kin_mem.kin_delay_aa > 0 {
                kin_mem.kin_nni - 1 - kin_mem.kin_delay_aa
            } else {
                kin_mem.kin_nni - 1
            };

            /* apply Anderson acceleration */
            /* C: AndersonAcc(kin_mem, fval, delta(=vtemp1), unew, uu,
               iter_aa, R_aa, gamma_aa) — fields take()n out around the
               call */
            let gval = std::mem::take(&mut kin_mem.kin_fval);
            let mut fv = std::mem::take(&mut kin_mem.kin_vtemp1); /* delta (temp) */
            let mut x = std::mem::take(&mut kin_mem.kin_unew);
            let xold = std::mem::take(&mut kin_mem.kin_uu);
            let mut r = std::mem::take(&mut kin_mem.kin_R_aa);
            let mut gamma = std::mem::take(&mut kin_mem.kin_gamma_aa);
            let retval =
                AndersonAcc(kin_mem, &gval, &mut fv, &mut x, &xold, iter_aa, &mut r,
                            &mut gamma);
            kin_mem.kin_fval = gval;
            kin_mem.kin_vtemp1 = fv;
            kin_mem.kin_unew = x;
            kin_mem.kin_uu = xold;
            kin_mem.kin_R_aa = r;
            kin_mem.kin_gamma_aa = gamma;
            if retval != 0 {
                ret = retval;
                break;
            }

            /* tolerance adjustment (first iteration is standard fixed point) */
            if iter_aa == 0
                && (kin_mem.kin_damping_aa || kin_mem.kin_damping_fn.is_some())
            {
                tolfac = kin_mem.kin_beta;
            } else {
                tolfac = ONE;
            }
        }

        /* compute change between iterations */
        /* delta = kin_vtemp1 */
        N_VLinearSum(ONE, &kin_mem.kin_unew, -ONE, &kin_mem.kin_uu,
                     &mut kin_mem.kin_vtemp1);

        /* measure || g(x) - x || */
        /* KINScFNorm(kin_mem, delta, fscale) with v = delta = vtemp1:
           N_VProd(fscale, delta, vtemp1) — output aliases the input */
        {
            let KINMem { kin_vtemp1, kin_fscale, .. } = kin_mem;
            kin_vtemp1.prod_with(kin_fscale);
        }
        kin_mem.kin_fnorm = N_VMaxNorm(&kin_mem.kin_vtemp1);

        /* (KINPrintInfo PRNT_FMAX / PRNT_NNI: SUNDIALS_LOGGING_LEVEL >= INFO only) */

        /* Check if the maximum number of iterations is reached */
        if kin_mem.kin_nni >= kin_mem.kin_mxiter {
            ret = KIN_MAXITER_REACHED;
        }
        if kin_mem.kin_fnorm <= tolfac * kin_mem.kin_fnormtol {
            ret = KIN_SUCCESS;
        }

        /* Update the solution if taking another iteration or returning the newest
           iterate. Otherwise return the solution consistent with the last function
           evaluation. */
        if ret == CONTINUE_ITERATIONS || kin_mem.kin_ret_newest {
            N_VScale(ONE, &kin_mem.kin_unew, &mut kin_mem.kin_uu);
        }
    } /* end of loop; return */

    /* (KINPrintInfo PRNT_RETVAL: SUNDIALS_LOGGING_LEVEL >= INFO only) */

    ret
}

/*
 * ========================================================================
 * Anderson Acceleration
 * ========================================================================
 */

fn AndersonAccQRDelete(kin_mem: &mut KINMem, Q: &mut [NVector], R: &mut [f64],
                       depth: i32) -> i32 {
    /* Delete left-most column vector from QR factorization */
    let depth = depth as usize;

    for i in 0..depth - 1 {
        let a = R[(i + 1) * depth + i];
        let b = R[(i + 1) * depth + i + 1];
        let temp = SUNRsqrt(a * a + b * b);
        let c = a / temp;
        let s = b / temp;
        R[(i + 1) * depth + i] = temp;
        R[(i + 1) * depth + i + 1] = ZERO;
        /* OK to reuse temp */
        if i < depth - 1 {
            for j in i + 2..depth {
                let a = R[j * depth + i];
                let b = R[j * depth + i + 1];
                let temp = c * a + s * b;
                R[j * depth + i + 1] = -s * a + c * b;
                R[j * depth + i] = temp;
            }
        }
        {
            let (q_head, q_tail) = Q.split_at_mut(i + 1);
            let qi = &mut q_head[i];
            let qip1 = &mut q_tail[0];
            N_VLinearSum(c, qi, s, qip1, &mut kin_mem.kin_vtemp2);
            /* N_VLinearSum(-s, Q[i], c, Q[i+1], Q[i+1]) — output aliases
               the second operand; the serial kernel's a==b / a==-b
               dispatch is reproduced (a = -s, b = c) */
            if -s == c {
                for k in 0..qip1.data.len() {
                    qip1.data[k] = -s * (qi.data[k] + qip1.data[k]);
                }
            } else if s == c {
                for k in 0..qip1.data.len() {
                    qip1.data[k] = -s * (qi.data[k] - qip1.data[k]);
                }
            } else {
                for k in 0..qip1.data.len() {
                    qip1.data[k] = -s * qi.data[k] + c * qip1.data[k];
                }
            }
            N_VScale(ONE, &kin_mem.kin_vtemp2, qi);
        }
    }

    /* Shift R to the left by one. */
    for i in 1..depth {
        for j in 0..depth - 1 {
            R[(i - 1) * depth + j] = R[i * depth + j];
        }
    }

    /* If ICWY orthogonalization, then update T */
    if kin_mem.kin_orth_aa == KIN_ORTH_ICWY {
        /* kin_dot_prod_sb is always SUNFALSE for the serial N_Vector (no
           nvdotprodmultiallreduce op; see KINInitOrth), so only the
           standard branch is reachable. */
        kin_mem.kin_T_aa[0] = ONE;
        for i in 2..depth {
            /* N_VDotProdMulti(i-1, Q[i-1], Q, T_aa + (i-1)*depth) */
            for j in 0..i - 1 {
                kin_mem.kin_T_aa[(i - 1) * depth + j] = N_VDotProd(&Q[i - 1], &Q[j]);
            }
            kin_mem.kin_T_aa[(i - 1) * depth + (i - 1)] = ONE;
        }
    }

    KIN_SUCCESS
}

/* AndersonAcc
 *
 * (C signature: AndersonAcc(kin_mem, gval, fv, x, xold, iter, R,
 * gamma); the vector/array arguments are KINMem fields — kin_gval or
 * kin_fval, kin_vtemp1, kin_unew, kin_uu, kin_R_aa, kin_gamma_aa —
 * take()n out by the callers KINFP / KINPicardAA. The kin_cv/kin_Xv
 * fused-operation scratch arrays become locals; the serial fused
 * N_VLinearCombination kernel is reproduced inline.)
 */
fn AndersonAcc(kin_mem: &mut KINMem, gval: &NVector, fv: &mut NVector, x: &mut NVector,
               xold: &NVector, iter: i64, R: &mut [f64], gamma: &mut [f64]) -> i32 {
    /* Compute residual F(x) = G(x_old) - x_old */
    N_VLinearSum(ONE, gval, -ONE, xold, fv);

    if iter > 0 {
        /* If we've filled the acceleration subspace, start recycling */
        if kin_mem.kin_current_depth == kin_mem.kin_m_aa {
            /* Move the left-most column vector (oldest value) to the end so it gets
               overwritten with the newest value below. */
            kin_mem.kin_dg_aa.rotate_left(1);
            kin_mem.kin_df_aa.rotate_left(1);

            /* Delete left-most column vector from QR factorization */
            let mut q = std::mem::take(&mut kin_mem.kin_q_aa);
            let retval =
                AndersonAccQRDelete(kin_mem, &mut q, R, kin_mem.kin_m_aa as i32);
            kin_mem.kin_q_aa = q;
            if retval != 0 {
                return retval;
            }

            kin_mem.kin_current_depth -= 1;
        }

        /* compute dg_new = gval - gval_old */
        {
            let depth = kin_mem.kin_current_depth as usize;
            let KINMem { kin_gold_aa, kin_dg_aa, .. } = kin_mem;
            N_VLinearSum(ONE, gval, -ONE, kin_gold_aa, &mut kin_dg_aa[depth]);
        }

        /* compute df_new = fval - fval_old */
        {
            let depth = kin_mem.kin_current_depth as usize;
            let KINMem { kin_fold_aa, kin_df_aa, .. } = kin_mem;
            N_VLinearSum(ONE, fv, -ONE, kin_fold_aa, &mut kin_df_aa[depth]);
        }

        kin_mem.kin_current_depth += 1;
    }

    /* (KINPrintInfo PRNT_OTHER "current_depth": logging-level INFO only) */

    N_VScale(ONE, gval, &mut kin_mem.kin_gold_aa);
    N_VScale(ONE, fv, &mut kin_mem.kin_fold_aa);

    /* on first iteration, do fixed point update */
    if kin_mem.kin_current_depth == 0 {
        if kin_mem.kin_damping_aa || kin_mem.kin_damping_fn.is_some() {
            if let Some(damping_fn) = kin_mem.kin_damping_fn {
                let retval = {
                    let KINMem { kin_nni, kin_user_data, kin_beta_aa, .. } = kin_mem;
                    damping_fn(*kin_nni, xold, gval, &[], 0, kin_user_data, kin_beta_aa)
                };
                if retval != 0 {
                    KINProcessError(Some(kin_mem), KIN_DAMPING_FN_ERR, line!(),
                                    "AndersonAcc", file!(),
                                    "The damping function failed.");
                    return KIN_DAMPING_FN_ERR;
                }
                if kin_mem.kin_beta_aa <= ZERO || kin_mem.kin_beta_aa > ONE {
                    KINProcessError(Some(kin_mem), KIN_DAMPING_FN_ERR, line!(),
                                    "AndersonAcc", file!(),
                                    "The damping parameter is outside of the range (0, 1].");
                    return KIN_DAMPING_FN_ERR;
                }
            }

            /* damped fixed point */
            N_VLinearSum(ONE - kin_mem.kin_beta_aa, xold, kin_mem.kin_beta_aa, gval, x);
        } else {
            /* standard fixed point */
            N_VScale(ONE, gval, x);
        }

        return KIN_SUCCESS;
    }

    /* Add a column to the QR factorization */

    if kin_mem.kin_current_depth == 1 {
        R[0] = SUNRsqrt(N_VDotProd(&kin_mem.kin_df_aa[0], &kin_mem.kin_df_aa[0]));
        let alfa = ONE / R[0];
        {
            let KINMem { kin_df_aa, kin_q_aa, .. } = kin_mem;
            N_VScale(alfa, &kin_df_aa[0], &mut kin_q_aa[0]);
        }
    } else {
        /* C: kin_qr_func(q_aa, R, df_aa[depth-1], depth-1, m_aa, qr_data).
           For ICWY the C qr_data->temp_array aliases kin_T_aa (also
           updated by AndersonAccQRDelete): the canonical kin_T_aa is
           swapped into the owned SUNQRData around the call. */
        let qr_func = kin_mem.kin_qr_func.unwrap();
        let mut qr_data = kin_mem.kin_qr_data.take().unwrap();
        if kin_mem.kin_orth_aa == KIN_ORTH_ICWY {
            std::mem::swap(&mut qr_data.temp_array, &mut kin_mem.kin_T_aa);
        }
        let mut q = std::mem::take(&mut kin_mem.kin_q_aa);
        /* (C discards the SUNQRAdd return value) */
        qr_func(&mut q, R,
                &kin_mem.kin_df_aa[(kin_mem.kin_current_depth - 1) as usize],
                (kin_mem.kin_current_depth - 1) as i32, kin_mem.kin_m_aa as i32,
                &mut qr_data);
        kin_mem.kin_q_aa = q;
        if kin_mem.kin_orth_aa == KIN_ORTH_ICWY {
            std::mem::swap(&mut qr_data.temp_array, &mut kin_mem.kin_T_aa);
        }
        kin_mem.kin_qr_data = Some(qr_data);
    }

    /* Adjust the depth */
    if let Some(depth_fn) = kin_mem.kin_depth_fn {
        let mut new_depth = kin_mem.kin_current_depth;

        let retval = {
            let KINMem { kin_nni, kin_df_aa, kin_current_depth, kin_user_data, .. } =
                kin_mem;
            depth_fn(*kin_nni, xold, gval, fv, kin_df_aa, R, *kin_current_depth,
                     kin_user_data, &mut new_depth, &mut [])
        };
        if retval != 0 {
            KINProcessError(Some(kin_mem), KIN_DEPTH_FN_ERR, line!(), "AndersonAcc",
                            file!(), "The depth function failed.");
            return KIN_DEPTH_FN_ERR;
        }

        new_depth = std::cmp::min(new_depth, kin_mem.kin_current_depth);
        new_depth = std::cmp::max(new_depth, 0);

        /* (KINPrintInfo PRNT_OTHER "new_depth": logging-level INFO only) */

        if new_depth == 0 {
            kin_mem.kin_current_depth = new_depth;

            /* do fixed point update */
            if kin_mem.kin_damping_aa || kin_mem.kin_damping_fn.is_some() {
                if let Some(damping_fn) = kin_mem.kin_damping_fn {
                    let retval = {
                        let KINMem { kin_nni, kin_user_data, kin_beta_aa, .. } = kin_mem;
                        damping_fn(*kin_nni, xold, gval, &[], 0, kin_user_data,
                                   kin_beta_aa)
                    };
                    if retval != 0 {
                        KINProcessError(Some(kin_mem), KIN_DAMPING_FN_ERR, line!(),
                                        "AndersonAcc", file!(),
                                        "The damping function failed.");
                        return KIN_DAMPING_FN_ERR;
                    }
                    if kin_mem.kin_beta_aa <= ZERO || kin_mem.kin_beta_aa > ONE {
                        KINProcessError(Some(kin_mem), KIN_DAMPING_FN_ERR, line!(),
                                        "AndersonAcc", file!(),
                                        "The damping parameter is outside of the range (0, 1].");
                        return KIN_DAMPING_FN_ERR;
                    }
                }

                /* damped fixed point */
                N_VLinearSum(ONE - kin_mem.kin_beta_aa, xold, kin_mem.kin_beta_aa, gval,
                             x);
            } else {
                /* standard fixed point */
                N_VScale(ONE, gval, x);
            }

            return KIN_SUCCESS;
        }

        /* TODO(DJG): In the future, update QRDelete to support removing arbitrary
           columns from the factorization */
        if new_depth < kin_mem.kin_current_depth {
            /* Remove columns from the left one at a time */
            for _j in 0..(kin_mem.kin_current_depth - new_depth) {
                let depth = kin_mem.kin_current_depth as usize;
                kin_mem.kin_dg_aa[..depth].rotate_left(1);
                kin_mem.kin_df_aa[..depth].rotate_left(1);

                let mut q = std::mem::take(&mut kin_mem.kin_q_aa);
                let retval = AndersonAccQRDelete(kin_mem, &mut q, R,
                                                 kin_mem.kin_current_depth as i32);
                kin_mem.kin_q_aa = q;
                if retval != 0 {
                    return retval;
                }

                kin_mem.kin_current_depth -= 1;
            }
        }
    }

    /* Solve least squares problem and update solution */
    let lAA = kin_mem.kin_current_depth;

    /* Compute Q^T fv (fused N_VDotProdMulti, serial kernel inline;
       the serial op cannot fail, so C's KIN_VECTOROP_ERR path is
       unreachable) */
    for j in 0..lAA as usize {
        gamma[j] = N_VDotProd(fv, &kin_mem.kin_q_aa[j]);
    }

    /* Compute the damping factor before overwriting gamma below so we can pass
       gamma = Q^T fv (just computed above) to the damping function as it can be
       used to compute the acceleration gain = sqrt(1 - ||Q^T fv||^2/||fv||^2). */
    if let Some(damping_fn) = kin_mem.kin_damping_fn {
        let retval = {
            let KINMem { kin_nni, kin_user_data, kin_beta_aa, .. } = kin_mem;
            damping_fn(*kin_nni, xold, gval, &gamma[..lAA as usize], lAA, kin_user_data,
                       kin_beta_aa)
        };
        if retval != 0 {
            KINProcessError(Some(kin_mem), KIN_DAMPING_FN_ERR, line!(), "AndersonAcc",
                            file!(), "The damping function failed.");
            return KIN_DAMPING_FN_ERR;
        }
        if kin_mem.kin_beta_aa <= ZERO || kin_mem.kin_beta_aa > ONE {
            KINProcessError(Some(kin_mem), KIN_DAMPING_FN_ERR, line!(), "AndersonAcc",
                            file!(),
                            "The damping parameter is outside of the range (0, 1].");
            return KIN_DAMPING_FN_ERR;
        }
    }

    /* set arrays for fused vector operation (C kin_cv/kin_Xv scratch) */
    let mut cv: Vec<f64> = Vec::with_capacity(2 * (kin_mem.kin_m_aa as usize + 1));
    let mut Xv: Vec<&NVector> = Vec::with_capacity(2 * (kin_mem.kin_m_aa as usize + 1));
    cv.push(ONE);
    Xv.push(gval);

    /* Solve the upper triangular system R gamma = Q^T fv */
    for i in (0..lAA).rev() {
        for j in i + 1..lAA {
            gamma[i as usize] -=
                R[(j * kin_mem.kin_m_aa + i) as usize] * gamma[j as usize];
        }
        gamma[i as usize] /= R[(i * kin_mem.kin_m_aa + i) as usize];

        cv.push(-gamma[i as usize]);
        Xv.push(&kin_mem.kin_dg_aa[i as usize]);
    }

    /* if enabled, apply damping */
    if kin_mem.kin_damping_aa || kin_mem.kin_damping_fn.is_some() {
        let onembeta = ONE - kin_mem.kin_beta_aa;
        cv.push(-onembeta);
        Xv.push(&*fv);
        for i in (0..lAA).rev() {
            cv.push(onembeta * gamma[i as usize]);
            Xv.push(&kin_mem.kin_df_aa[i as usize]);
        }
    }

    /* update solution: N_VLinearCombination(nvec, cv, Xv, x) — serial
       fused kernel reproduced inline (nvec >= 2 always: cv[0] plus at
       least one gamma entry; Xv[0] = gval never aliases x, so the
       general branch matches the serial VLinearCombination path; the
       serial op cannot fail, so C's KIN_VECTOROP_ERR is unreachable) */
    let nvec = cv.len();
    if nvec == 2 {
        N_VLinearSum(cv[0], Xv[0], cv[1], Xv[1], x);
    } else {
        for j in 0..x.data.len() {
            x.data[j] = cv[0] * Xv[0].data[j];
        }
        for i in 1..nvec {
            for j in 0..x.data.len() {
                x.data[j] += cv[i] * Xv[i].data[j];
            }
        }
    }

    KIN_SUCCESS
}
