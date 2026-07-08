/* -----------------------------------------------------------------
 * Translated from src/cvodes/cvodea_io.c (CVODES 7.7.0).
 *
 * Optional input and output functions for the adjoint module in
 * the CVODES solver.
 *
 * Port notes:
 *  - Every ***B function shares the C preamble (check ASA
 *    initialized, check `which`, find the CVodeBMem entry): the
 *    helper cva_io_which_index performs it and returns the Vec
 *    index of the backward problem (cv_index == position).
 *  - The C void* addresses of check points (CVadjCheckPointRec
 *    my_addr/next_addr, CVodeGetAdjCurrentCheckPoint) become
 *    indices into the ck_mem Vec (Option<usize>; None = NULL),
 *    ordered as the C linked list is walked (newest first).
 * -----------------------------------------------------------------*/

use crate::cvodes_impl::*;
use crate::cvodes_io::{
    CVodeSetConstraints, CVodeSetInitStep, CVodeSetMaxNumSteps, CVodeSetMaxOrd, CVodeSetMaxStep,
    CVodeSetMinStep, CVodeSetQuadErrCon, CVodeSetStabLimDet,
};
use crate::nvector_serial::*;
use crate::sundials_nonlinearsolver::NonlinearSolver;
use crate::sundials_types::*;

const ONE: f64 = 1.0;

/* Shared preamble of every ***B function (see cvodea.rs). */
fn cva_io_which_index(cv_mem: &mut CVodeMem, which: i32, fname: &str) -> Result<usize, i32> {
    if !cv_mem.cv_adjMallocDone {
        cvProcessError(Some(cv_mem), CV_NO_ADJ, line!(), fname, file!(), MSGCV_NO_ADJ);
        return Err(CV_NO_ADJ);
    }

    if which >= cv_mem.cv_adj_mem.as_ref().unwrap().ca_nbckpbs {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), fname, file!(), MSGCV_BAD_WHICH);
        return Err(CV_ILL_INPUT);
    }

    /* Find the CVodeBMem entry corresponding to which */
    Ok(cv_mem
        .cv_adj_mem
        .as_ref()
        .unwrap()
        .cvB_mem
        .iter()
        .position(|b| b.cv_index == which)
        .unwrap())
}

/*
 * -----------------------------------------------------------------
 * Optional input functions for ASA
 * -----------------------------------------------------------------
 */

pub fn CVodeSetAdjNoSensi(cv_mem: &mut CVodeMem) -> i32 {
    /* Was ASA initialized? */
    if !cv_mem.cv_adjMallocDone {
        cvProcessError(Some(cv_mem), CV_NO_ADJ, line!(), "CVodeSetAdjNoSensi", file!(),
                       MSGCV_NO_ADJ);
        return CV_NO_ADJ;
    }

    cv_mem.cv_adj_mem.as_mut().unwrap().ca_IMstoreSensi = SUNFALSE;

    CV_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Optional input functions for backward integration
 * -----------------------------------------------------------------
 */

pub fn CVodeSetNonlinearSolverB(cv_mem: &mut CVodeMem, which: i32, nls: NonlinearSolver) -> i32 {
    let idx = match cva_io_which_index(cv_mem, which, "CVodeSetNonlinearSolverB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    crate::cvodes_nls::CVodeSetNonlinearSolver(
        &mut cv_mem.cv_adj_mem.as_mut().unwrap().cvB_mem[idx].cv_mem,
        nls,
    )
}

pub fn CVodeSetUserDataB(cv_mem: &mut CVodeMem, which: i32, user_dataB: UserData) -> i32 {
    let idx = match cva_io_which_index(cv_mem, which, "CVodeSetUserDataB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    cv_mem.cv_adj_mem.as_mut().unwrap().cvB_mem[idx].cv_user_data = user_dataB;

    CV_SUCCESS
}

/* (C: int CVodeGetUserDataB(cvode_mem, which, void** user_dataB); the
   pointer out-parameter becomes a borrow of the stored user data, like
   the forward CVodeGetUserData.) */
pub fn CVodeGetUserDataB(cv_mem: &mut CVodeMem, which: i32) -> Result<&mut UserData, i32> {
    let idx = match cva_io_which_index(cv_mem, which, "CVodeGetUserDataB") {
        Ok(i) => i,
        Err(e) => return Err(e),
    };

    Ok(&mut cv_mem.cv_adj_mem.as_mut().unwrap().cvB_mem[idx].cv_user_data)
}

pub fn CVodeSetMaxOrdB(cv_mem: &mut CVodeMem, which: i32, maxordB: i32) -> i32 {
    let idx = match cva_io_which_index(cv_mem, which, "CVodeSetMaxOrdB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    CVodeSetMaxOrd(&mut cv_mem.cv_adj_mem.as_mut().unwrap().cvB_mem[idx].cv_mem, maxordB)
}

pub fn CVodeSetMaxNumStepsB(cv_mem: &mut CVodeMem, which: i32, mxstepsB: i64) -> i32 {
    let idx = match cva_io_which_index(cv_mem, which, "CVodeSetMaxNumStepsB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    CVodeSetMaxNumSteps(&mut cv_mem.cv_adj_mem.as_mut().unwrap().cvB_mem[idx].cv_mem, mxstepsB)
}

pub fn CVodeSetStabLimDetB(cv_mem: &mut CVodeMem, which: i32, stldetB: bool) -> i32 {
    let idx = match cva_io_which_index(cv_mem, which, "CVodeSetStabLimDetB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    CVodeSetStabLimDet(&mut cv_mem.cv_adj_mem.as_mut().unwrap().cvB_mem[idx].cv_mem, stldetB)
}

pub fn CVodeSetInitStepB(cv_mem: &mut CVodeMem, which: i32, hinB: f64) -> i32 {
    let idx = match cva_io_which_index(cv_mem, which, "CVodeSetInitStepB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    CVodeSetInitStep(&mut cv_mem.cv_adj_mem.as_mut().unwrap().cvB_mem[idx].cv_mem, hinB)
}

pub fn CVodeSetMinStepB(cv_mem: &mut CVodeMem, which: i32, hminB: f64) -> i32 {
    let idx = match cva_io_which_index(cv_mem, which, "CVodeSetMinStepB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    CVodeSetMinStep(&mut cv_mem.cv_adj_mem.as_mut().unwrap().cvB_mem[idx].cv_mem, hminB)
}

pub fn CVodeSetMaxStepB(cv_mem: &mut CVodeMem, which: i32, hmaxB: f64) -> i32 {
    let idx = match cva_io_which_index(cv_mem, which, "CVodeSetMaxStepB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    CVodeSetMaxStep(&mut cv_mem.cv_adj_mem.as_mut().unwrap().cvB_mem[idx].cv_mem, hmaxB)
}

pub fn CVodeSetConstraintsB(
    cv_mem: &mut CVodeMem,
    which: i32,
    constraintsB: Option<&NVector>,
) -> i32 {
    /* (The C prints MSGCV_NO_ADJ but is missing the early return here -
       reading through the NULL adjoint memory is undefined behavior; the
       port returns CV_NO_ADJ like every other ***B function.) */
    let idx = match cva_io_which_index(cv_mem, which, "CVodeSetConstraintsB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    CVodeSetConstraints(
        &mut cv_mem.cv_adj_mem.as_mut().unwrap().cvB_mem[idx].cv_mem,
        constraintsB,
    )
}

/*
 * CVodeSetQuad*B
 *
 * Wrappers for the backward phase around the corresponding
 * CVODES quadrature optional input functions
 */

pub fn CVodeSetQuadErrConB(cv_mem: &mut CVodeMem, which: i32, errconQB: bool) -> i32 {
    let idx = match cva_io_which_index(cv_mem, which, "CVodeSetQuadErrConB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    CVodeSetQuadErrCon(&mut cv_mem.cv_adj_mem.as_mut().unwrap().cvB_mem[idx].cv_mem, errconQB)
}

/*
 * -----------------------------------------------------------------
 * Optional output functions for backward integration
 * -----------------------------------------------------------------
 */

/*
 * CVodeGetAdjCVodeBmem
 *
 * This function returns a pointer to the CVODES memory allocated for
 * the backward problem. This pointer can then be used to call any of
 * the CVodeGet* CVODES routines to extract optional output for the
 * backward integration phase.
 * (C returns void*, NULL on error; here Option<&mut CVodeMem>.)
 */
pub fn CVodeGetAdjCVodeBmem(cv_mem: &mut CVodeMem, which: i32) -> Option<&mut CVodeMem> {
    if !cv_mem.cv_adjMallocDone {
        cvProcessError(Some(cv_mem), 0, line!(), "CVodeGetAdjCVodeBmem", file!(), MSGCV_NO_ADJ);
        return None;
    }

    if which >= cv_mem.cv_adj_mem.as_ref().unwrap().ca_nbckpbs {
        cvProcessError(Some(cv_mem), 0, line!(), "CVodeGetAdjCVodeBmem", file!(), MSGCV_BAD_WHICH);
        return None;
    }

    let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
    let idx = ca_mem.cvB_mem.iter().position(|b| b.cv_index == which).unwrap();
    Some(&mut ca_mem.cvB_mem[idx].cv_mem)
}

/*
 * CVadjCheckPointRec (cvodes.h)
 *
 * The C my_addr/next_addr void* fields become indices into the
 * ck_mem Vec (None = the C NULL next pointer of the initial check
 * point).
 */
#[derive(Debug, Clone, Copy, Default)]
pub struct CVadjCheckPointRec {
    pub my_addr: usize,
    pub next_addr: Option<usize>,
    pub t0: f64,
    pub t1: f64,
    pub nstep: i64,
    pub order: i32,
    pub step: f64,
}

/*
 * CVodeGetAdjCheckPointsInfo
 *
 * This routine loads an array of nckpnts structures of type
 * CVadjCheckPointRec. The user must allocate space for ckpnt.
 * (Walked newest check point first, like the C linked list.)
 */
pub fn CVodeGetAdjCheckPointsInfo(cv_mem: &mut CVodeMem, ckpnt: &mut [CVadjCheckPointRec]) -> i32 {
    /* Was ASA initialized? */
    if !cv_mem.cv_adjMallocDone {
        cvProcessError(Some(cv_mem), CV_NO_ADJ, line!(), "CVodeGetAdjCheckPointsInfo", file!(),
                       MSGCV_NO_ADJ);
        return CV_NO_ADJ;
    }

    let ca_mem = cv_mem.cv_adj_mem.as_ref().unwrap();

    let mut i = 0usize;
    for idx in (0..ca_mem.ck_mem.len()).rev() {
        let ck_mem = &ca_mem.ck_mem[idx];
        ckpnt[i].my_addr = idx;
        ckpnt[i].next_addr = if idx > 0 { Some(idx - 1) } else { None };
        ckpnt[i].t0 = ck_mem.ck_t0;
        ckpnt[i].t1 = ck_mem.ck_t1;
        ckpnt[i].nstep = ck_mem.ck_nst;
        ckpnt[i].order = ck_mem.ck_q;
        ckpnt[i].step = ck_mem.ck_h;
        i += 1;
    }

    CV_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Undocumented Development User-Callable Functions
 * -----------------------------------------------------------------
 */

/*
 * CVodeGetAdjDataPointHermite
 *
 * This routine returns the solution stored in the data structure
 * at the 'which' data point. Cubic Hermite interpolation.
 */
pub fn CVodeGetAdjDataPointHermite(
    cv_mem: &mut CVodeMem,
    which: i32,
    t: &mut f64,
    y: Option<&mut NVector>,
    yd: Option<&mut NVector>,
) -> i32 {
    /* Was ASA initialized? */
    if !cv_mem.cv_adjMallocDone {
        cvProcessError(Some(cv_mem), CV_NO_ADJ, line!(), "CVodeGetAdjDataPointHermite", file!(),
                       MSGCV_NO_ADJ);
        return CV_NO_ADJ;
    }

    if cv_mem.cv_adj_mem.as_ref().unwrap().ca_IMtype != CV_HERMITE {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeGetAdjDataPointHermite",
                       file!(), MSGCV_WRONG_INTERP);
        return CV_ILL_INPUT;
    }

    let ca_mem = cv_mem.cv_adj_mem.as_ref().unwrap();
    let d = &ca_mem.dt_mem[which as usize];

    *t = d.t;

    if let DtpntContent::Hermite { y: cy, yd: cyd, .. } = &d.content {
        if let Some(y) = y {
            N_VScale(ONE, cy, y);
        }
        if let Some(yd) = yd {
            N_VScale(ONE, cyd, yd);
        }
    }

    CV_SUCCESS
}

/*
 * CVodeGetAdjDataPointPolynomial
 *
 * This routine returns the solution stored in the data structure
 * at the 'which' data point. Polynomial interpolation.
 */
pub fn CVodeGetAdjDataPointPolynomial(
    cv_mem: &mut CVodeMem,
    which: i32,
    t: &mut f64,
    order: &mut i32,
    y: Option<&mut NVector>,
) -> i32 {
    /* Was ASA initialized? */
    if !cv_mem.cv_adjMallocDone {
        cvProcessError(Some(cv_mem), CV_NO_ADJ, line!(), "CVodeGetAdjDataPointPolynomial",
                       file!(), MSGCV_NO_ADJ);
        return CV_NO_ADJ;
    }

    if cv_mem.cv_adj_mem.as_ref().unwrap().ca_IMtype != CV_POLYNOMIAL {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeGetAdjDataPointPolynomial",
                       file!(), MSGCV_WRONG_INTERP);
        return CV_ILL_INPUT;
    }

    let ca_mem = cv_mem.cv_adj_mem.as_ref().unwrap();
    let d = &ca_mem.dt_mem[which as usize];

    *t = d.t;

    if let DtpntContent::Polynomial { y: cy, order: corder, .. } = &d.content {
        if let Some(y) = y {
            N_VScale(ONE, cy, y);
        }
        *order = *corder;
    }

    CV_SUCCESS
}

/*
 * CVodeGetAdjCurrentCheckPoint
 *
 * Returns the address of the 'active' check point.
 * (C writes a void*; here the index into ck_mem, None = NULL.)
 */
pub fn CVodeGetAdjCurrentCheckPoint(cv_mem: &mut CVodeMem, addr: &mut Option<usize>) -> i32 {
    /* Was ASA initialized? */
    if !cv_mem.cv_adjMallocDone {
        cvProcessError(Some(cv_mem), CV_NO_ADJ, line!(), "CVodeGetAdjCurrentCheckPoint", file!(),
                       MSGCV_NO_ADJ);
        return CV_NO_ADJ;
    }

    *addr = cv_mem.cv_adj_mem.as_ref().unwrap().ca_ckpntData;

    CV_SUCCESS
}
