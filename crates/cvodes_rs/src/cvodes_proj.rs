/* ---------------------------------------------------------------------------
 * Translated from src/cvodes/cvodes_proj.c (CVODES 7.7.0).
 * Implementation file for projections in CVODE.
 *
 * The C functions take `void* cvode_mem` and start with a NULL check;
 * here the memory is `&mut CVodeMem`, which cannot be null, so those
 * checks vanish. cvProjFree is not needed (Box is dropped by RAII).
 * ---------------------------------------------------------------------------*/
use crate::cvodes::{cvRescale, cvRestore};
use crate::cvodes_impl::*;
use crate::cvodes_proj_impl::*;
use crate::nvector_serial::{N_VScale, N_VWrmsNorm};
use crate::sundials_math::{SUNMAX, SUNRabs};
use crate::sundials_types::*;

/* Private constants */
const ZERO: f64 = 0.0; /* real 0.0 */
const ONE: f64 = 1.0; /* real 1.0 */

/* (ONEPSM = 1.000001 comes from cvodes_impl) */

/* ===========================================================================
 * Exported Functions - projection initialization
 * ===========================================================================*/

/* -----------------------------------------------------------------------------
 * CVodeSetProjFn sets a user defined projection function
 * ---------------------------------------------------------------------------*/
pub fn CVodeSetProjFn(cv_mem: &mut CVodeMem, pfun: CVProjFn) -> i32 {
    /* (The projection function cannot be NULL in Rust: `pfun` is a fn item.) */

    /* Check for compatible method */
    if cv_mem.cv_lmm != CV_BDF {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeSetProjFn", file!(),
                       "Projection is only supported with BDF methods.");
        return CV_ILL_INPUT;
    }

    /* Create the projection memory (if necessary) */
    cvProjCreate(&mut cv_mem.proj_mem);

    /* Shortcut to projection memory */
    let proj_mem = cv_mem.proj_mem.as_deref_mut().unwrap();

    /* User-defined projection */
    proj_mem.internal_proj = SUNFALSE;

    /* Set the projection function */
    proj_mem.pfun = Some(pfun);

    /* Enable projection */
    cv_mem.proj_enabled = SUNTRUE;

    CV_SUCCESS
}

/* ===========================================================================
 * Exported Functions - projection set function
 * ===========================================================================*/

pub fn CVodeSetProjErrEst(cv_mem: &mut CVodeMem, onoff: bool) -> i32 {
    /* Access memory structures */
    let proj_mem = match cvAccessProjMem(cv_mem, "CVodeSetProjErrEst") {
        Ok(pm) => pm,
        Err(retval) => return retval,
    };

    /* Set projection error flag */
    proj_mem.err_proj = onoff;

    CV_SUCCESS
}

pub fn CVodeSetProjFrequency(cv_mem: &mut CVodeMem, freq: i64) -> i32 {
    /* Access memory structures */
    if let Err(retval) = cvAccessProjMem(cv_mem, "CVodeSetProjFrequency") {
        return retval;
    }

    /* Set projection frequency */
    if freq < 0 {
        /* Restore default */
        cv_mem.proj_mem.as_deref_mut().unwrap().freq = 1;
        cv_mem.proj_enabled = SUNTRUE;
    } else if freq == 0 {
        /* Disable projection */
        cv_mem.proj_mem.as_deref_mut().unwrap().freq = 0;
        cv_mem.proj_enabled = SUNFALSE;
    } else {
        /* Enable projection at given frequency */
        cv_mem.proj_mem.as_deref_mut().unwrap().freq = freq;
        cv_mem.proj_enabled = SUNTRUE;
    }

    CV_SUCCESS
}

pub fn CVodeSetMaxNumProjFails(cv_mem: &mut CVodeMem, max_fails: i32) -> i32 {
    /* Access memory structures */
    let proj_mem = match cvAccessProjMem(cv_mem, "CVodeSetMaxNumProjFails") {
        Ok(pm) => pm,
        Err(retval) => return retval,
    };

    /* Set maximum number of projection failures in a step attempt */
    if max_fails < 1 {
        /* Restore default */
        proj_mem.max_fails = PROJ_MAX_FAILS;
    } else {
        /* Update max number of fails */
        proj_mem.max_fails = max_fails;
    }

    CV_SUCCESS
}

pub fn CVodeSetEpsProj(cv_mem: &mut CVodeMem, eps: f64) -> i32 {
    /* Access memory structures */
    let proj_mem = match cvAccessProjMem(cv_mem, "CVodeSetEpsProj") {
        Ok(pm) => pm,
        Err(retval) => return retval,
    };

    /* Set the projection tolerance */
    if eps <= ZERO {
        /* Restore default */
        proj_mem.eps_proj = PROJ_EPS;
    } else {
        /* Update projection tolerance */
        proj_mem.eps_proj = eps;
    }

    CV_SUCCESS
}

pub fn CVodeSetProjFailEta(cv_mem: &mut CVodeMem, eta: f64) -> i32 {
    /* Access memory structures */
    let proj_mem = match cvAccessProjMem(cv_mem, "CVodeSetProjFailEta") {
        Ok(pm) => pm,
        Err(retval) => return retval,
    };

    /* Set the step size reduction factor for a projection failure */
    if (eta <= ZERO) || (eta > ONE) {
        /* Restore default */
        proj_mem.eta_pfail = PROJ_FAIL_ETA;
    } else {
        /* Update the eta value */
        proj_mem.eta_pfail = eta;
    }

    CV_SUCCESS
}

/* ===========================================================================
 * Exported Functions - projection get functions
 * ===========================================================================*/

pub fn CVodeGetNumProjEvals(cv_mem: &mut CVodeMem, nproj: &mut i64) -> i32 {
    /* Access memory structures */
    let proj_mem = match cvAccessProjMem(cv_mem, "CVodeGetNumProjEvals") {
        Ok(pm) => pm,
        Err(retval) => return retval,
    };

    /* Get number of projection evaluations */
    *nproj = proj_mem.nproj;

    CV_SUCCESS
}

pub fn CVodeGetNumProjFails(cv_mem: &mut CVodeMem, npfails: &mut i64) -> i32 {
    /* Access memory structures */
    let proj_mem = match cvAccessProjMem(cv_mem, "CVodeGetNumProjFails") {
        Ok(pm) => pm,
        Err(retval) => return retval,
    };

    /* Get number of projection fails */
    *npfails = proj_mem.npfails;

    CV_SUCCESS
}

/* ===========================================================================
 * Internal Functions
 * ===========================================================================*/

/*
 * cvProjection
 *
 * For user supplied projection function, use ftemp as temporary storage
 * for the current error estimate (acor) and use tempv to store the
 * accumulated correction due to projection, acorP (tempv is not touched
 * until it is potentially used in cvCompleteStep).
 */

pub fn cvDoProjection(
    cv_mem: &mut CVodeMem,
    nflagPtr: &mut i32,
    saved_t: f64,
    npfailPtr: &mut i32,
) -> i32 {
    /* Access projection memory */
    let Some(mut proj_mem) = cv_mem.proj_mem.take() else {
        cvProcessError(Some(cv_mem), CV_PROJ_MEM_NULL, line!(), "cvDoProjection", file!(),
                       MSG_CV_PROJ_MEM_NULL);
        return CV_PROJ_MEM_NULL;
    };

    /* proj_mem is taken out of cv_mem so the integrator memory can be
       borrowed freely while calling the user projection function */
    let retval = cvDoProjectionImpl(cv_mem, &mut proj_mem, nflagPtr, saved_t, npfailPtr);

    cv_mem.proj_mem = Some(proj_mem);

    retval
}

fn cvDoProjectionImpl(
    cv_mem: &mut CVodeMem,
    proj_mem: &mut CVodeProjMem,
    nflagPtr: &mut i32,
    saved_t: f64,
    npfailPtr: &mut i32,
) -> i32 {
    /* Use tempv to store acorP and, if projecting the error, ftemp to store
       errP (recall that in this case we did not allocate vectors to for
       acorP and errP). */

    /* Copy acor into errP (if projecting the error) */
    if proj_mem.err_proj {
        let CVodeMem { cv_acor, cv_ftemp, .. } = cv_mem;
        N_VScale(ONE, cv_acor, cv_ftemp);
    }

    /* Call the user projection function */
    let pfun = proj_mem.pfun.expect("projection function is NULL");
    let retval = {
        let CVodeMem { cv_tn, cv_y, cv_tempv, cv_ftemp, cv_user_data, .. } = cv_mem;
        let acorP = cv_tempv;
        let errP = if proj_mem.err_proj { Some(&mut *cv_ftemp) } else { None };
        pfun(*cv_tn, cv_y, acorP, proj_mem.eps_proj, errP, cv_user_data)
    };
    proj_mem.nproj += 1;

    /* This is not the first projection anymore */
    proj_mem.first_proj = SUNFALSE;

    /* Check the return value */
    if retval == CV_SUCCESS {
        /* Recompute acnrm to be used in error test (if projecting the error) */
        if proj_mem.err_proj {
            cv_mem.cv_acnrm = N_VWrmsNorm(&cv_mem.cv_ftemp, &cv_mem.cv_ewt);
        }

        /* The projection was successful, return now */
        cv_mem.proj_applied = SUNTRUE;
        return CV_SUCCESS;
    }

    /* The projection failed, update the return value */
    let retval = if retval < 0 { CV_PROJFUNC_FAIL } else { PROJFUNC_RECVR };

    /* Increment cumulative failure count and restore zn */
    proj_mem.npfails += 1;
    cvRestore(cv_mem, saved_t);

    /* Return if failed unrecoverably */
    if retval == CV_PROJFUNC_FAIL {
        return CV_PROJFUNC_FAIL;
    }

    /* Recoverable failure, increment failure count for this step attempt */
    *npfailPtr += 1;
    cv_mem.cv_etamax = ONE;

    /* Check for maximum number of failures or |h| = hmin */
    if (SUNRabs(cv_mem.cv_h) <= cv_mem.cv_hmin * ONEPSM) || (*npfailPtr == proj_mem.max_fails) {
        if retval == PROJFUNC_RECVR {
            return CV_REPTD_PROJFUNC_ERR;
        }
    }

    /* Reduce step size; return to reattempt the step */
    cv_mem.cv_eta = SUNMAX(proj_mem.eta_pfail, cv_mem.cv_hmin / SUNRabs(cv_mem.cv_h));
    *nflagPtr = PREV_PROJ_FAIL;
    cvRescale(cv_mem);

    PREDICT_AGAIN
}

pub fn cvProjInit(proj_mem: &mut CVodeProjMem) -> i32 {
    /* reset flags and counters */
    proj_mem.first_proj = SUNTRUE;
    proj_mem.nstlprj = 0;
    proj_mem.nproj = 0;
    proj_mem.npfails = 0;

    CV_SUCCESS
}

/* ===========================================================================
 * Utility Functions
 * ===========================================================================*/

fn cvProjCreate(proj_mem: &mut Option<Box<CVodeProjMem>>) {
    /* Allocate projection memory if necessary, otherwise return success */
    if proj_mem.is_none() {
        /* Initialize projection variables (cvProjSetDefaults) */
        *proj_mem = Some(Box::new(CVodeProjMem {
            internal_proj: SUNTRUE,
            err_proj: SUNTRUE,
            first_proj: SUNTRUE,

            freq: 1,
            nstlprj: 0,

            max_fails: PROJ_MAX_FAILS,

            pfun: None,

            eps_proj: PROJ_EPS,
            eta_pfail: PROJ_FAIL_ETA,

            nproj: 0,
            npfails: 0,
        }));
    }
}

fn cvAccessProjMem<'a>(
    cv_mem: &'a mut CVodeMem,
    fname: &str,
) -> Result<&'a mut CVodeProjMem, i32> {
    /* Access projection memory */
    if cv_mem.proj_mem.is_none() {
        cvProcessError(Some(cv_mem), CV_PROJ_MEM_NULL, line!(), fname, file!(),
                       MSG_CV_PROJ_MEM_NULL);
        return Err(CV_PROJ_MEM_NULL);
    }
    Ok(cv_mem.proj_mem.as_deref_mut().unwrap())
}
