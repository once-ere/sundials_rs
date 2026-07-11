/* -----------------------------------------------------------------
 * Translated from src/idas/idaa_io.c (IDAS 7.7.0).
 *
 * Optional input and output functions for the adjoint module in
 * the IDAS solver.
 *
 * Port notes (cvodea_io.rs conventions):
 *  - Every ***B function shares the C preamble (check ASA
 *    initialized, check `which`, find the IDABMem entry): the helper
 *    idaa_io_which_index performs it and returns the Vec index of
 *    the backward problem (ida_index == position).
 *  - The C void* addresses of check points (IDAadjCheckPointRec
 *    my_addr/next_addr, IDAGetAdjCurrentCheckPoint) become indices
 *    into the ck_mem Vec (Option<usize>; None = NULL), ordered as
 *    the C linked list is walked (newest first).
 *  - C IDAGetAdjIDABmem returns the nested solver memory as void*;
 *    here it is a borrow of the Box<IDAMem> (None on the C NULL
 *    error returns).
 * -----------------------------------------------------------------*/

use crate::idas_impl::*;
use crate::idas_io::{
    IDAGetConsistentIC, IDASetConstraints, IDASetId, IDASetInitStep, IDASetMaxNumSteps,
    IDASetMaxOrd, IDASetMaxStep, IDASetQuadErrCon, IDASetSuppressAlg,
};
use crate::idas_nls::IDASetNonlinearSolver;
use crate::nvector_serial::*;
use crate::sundials_nonlinearsolver::NonlinearSolver;
use crate::sundials_types::*;

const ONE: f64 = 1.0;

/* Shared preamble of every ***B function (see idaa.rs). */
fn idaa_io_which_index(ida_mem: &mut IDAMem, which: i32, fname: &str) -> Result<usize, i32> {
    /* Is ASA initialized? */
    if !ida_mem.ida_adjMallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_ADJ, line!(), fname, file!(), MSGAM_NO_ADJ);
        return Err(IDA_NO_ADJ);
    }

    /* Check the value of which */
    if which >= ida_mem.ida_adj_mem.as_ref().unwrap().ia_nbckpbs {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), fname, file!(), MSGAM_BAD_WHICH);
        return Err(IDA_ILL_INPUT);
    }

    /* Find the IDABMem entry in the linked list corresponding to which */
    Ok(ida_mem
        .ida_adj_mem
        .as_ref()
        .unwrap()
        .IDAB_mem
        .iter()
        .position(|b| b.ida_index == which)
        .unwrap())
}

/*
 * -----------------------------------------------------------------
 * Optional input functions for ASA
 * -----------------------------------------------------------------
 */

pub fn IDAAdjSetNoSensi(ida_mem: &mut IDAMem) -> i32 {
    /* Is ASA initialized? */
    if !ida_mem.ida_adjMallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_ADJ, line!(), "IDAAdjSetNoSensi", file!(),
                        MSGAM_NO_ADJ);
        return IDA_NO_ADJ;
    }

    ida_mem.ida_adj_mem.as_mut().unwrap().ia_storeSensi = SUNFALSE;

    IDA_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Optional input functions for backward integration
 * -----------------------------------------------------------------
 */

pub fn IDASetNonlinearSolverB(ida_mem: &mut IDAMem, which: i32, NLS: NonlinearSolver) -> i32 {
    let idx = match idaa_io_which_index(ida_mem, which, "IDASetNonlinearSolverB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    IDASetNonlinearSolver(
        &mut ida_mem.ida_adj_mem.as_mut().unwrap().IDAB_mem[idx].IDA_mem,
        NLS,
    )
}

pub fn IDASetUserDataB(ida_mem: &mut IDAMem, which: i32, user_dataB: UserData) -> i32 {
    let idx = match idaa_io_which_index(ida_mem, which, "IDASetUserDataB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    /* Set user data for this backward problem. */
    ida_mem.ida_adj_mem.as_mut().unwrap().IDAB_mem[idx].ida_user_data = user_dataB;

    IDA_SUCCESS
}

pub fn IDASetMaxOrdB(ida_mem: &mut IDAMem, which: i32, maxordB: i32) -> i32 {
    let idx = match idaa_io_which_index(ida_mem, which, "IDASetMaxOrdB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    IDASetMaxOrd(&mut ida_mem.ida_adj_mem.as_mut().unwrap().IDAB_mem[idx].IDA_mem, maxordB)
}

pub fn IDASetMaxNumStepsB(ida_mem: &mut IDAMem, which: i32, mxstepsB: i64) -> i32 {
    let idx = match idaa_io_which_index(ida_mem, which, "IDASetMaxNumStepsB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    IDASetMaxNumSteps(&mut ida_mem.ida_adj_mem.as_mut().unwrap().IDAB_mem[idx].IDA_mem, mxstepsB)
}

pub fn IDASetInitStepB(ida_mem: &mut IDAMem, which: i32, hinB: f64) -> i32 {
    let idx = match idaa_io_which_index(ida_mem, which, "IDASetInitStepB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    IDASetInitStep(&mut ida_mem.ida_adj_mem.as_mut().unwrap().IDAB_mem[idx].IDA_mem, hinB)
}

pub fn IDASetMaxStepB(ida_mem: &mut IDAMem, which: i32, hmaxB: f64) -> i32 {
    let idx = match idaa_io_which_index(ida_mem, which, "IDASetMaxStepB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    IDASetMaxStep(&mut ida_mem.ida_adj_mem.as_mut().unwrap().IDAB_mem[idx].IDA_mem, hmaxB)
}

pub fn IDASetSuppressAlgB(ida_mem: &mut IDAMem, which: i32, suppressalgB: bool) -> i32 {
    let idx = match idaa_io_which_index(ida_mem, which, "IDASetSuppressAlgB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    IDASetSuppressAlg(
        &mut ida_mem.ida_adj_mem.as_mut().unwrap().IDAB_mem[idx].IDA_mem,
        suppressalgB,
    )
}

pub fn IDASetIdB(ida_mem: &mut IDAMem, which: i32, idB: Option<&NVector>) -> i32 {
    let idx = match idaa_io_which_index(ida_mem, which, "IDASetIdB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    IDASetId(&mut ida_mem.ida_adj_mem.as_mut().unwrap().IDAB_mem[idx].IDA_mem, idB)
}

pub fn IDASetConstraintsB(ida_mem: &mut IDAMem, which: i32, constraintsB: Option<&NVector>) -> i32 {
    let idx = match idaa_io_which_index(ida_mem, which, "IDASetConstraintsB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    IDASetConstraints(
        &mut ida_mem.ida_adj_mem.as_mut().unwrap().IDAB_mem[idx].IDA_mem,
        constraintsB,
    )
}

pub fn IDASetQuadErrConB(ida_mem: &mut IDAMem, which: i32, errconQB: bool) -> i32 {
    let idx = match idaa_io_which_index(ida_mem, which, "IDASetQuadErrConB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    IDASetQuadErrCon(&mut ida_mem.ida_adj_mem.as_mut().unwrap().IDAB_mem[idx].IDA_mem, errconQB)
}

/*
 * -----------------------------------------------------------------
 * Optional output functions for backward integration
 * -----------------------------------------------------------------
 */

/*
 * IDAGetAdjIDABmem
 *
 * This function returns a pointer to the IDAS memory allocated for
 * the backward problem. This pointer can then be used to call any of
 * the IDAGet* IDAS routines to extract optional output for the
 * backward integration phase.  (C returns a void pointer or NULL;
 * here a borrow of the nested solver memory / None.)
 */
pub fn IDAGetAdjIDABmem(ida_mem: &mut IDAMem, which: i32) -> Option<&mut IDAMem> {
    /* Is ASA initialized? */
    if !ida_mem.ida_adjMallocDone {
        IDAProcessError(Some(ida_mem), 0, line!(), "IDAGetAdjIDABmem", file!(), MSGAM_NO_ADJ);
        return None;
    }

    /* Check the value of which */
    if which >= ida_mem.ida_adj_mem.as_ref().unwrap().ia_nbckpbs {
        IDAProcessError(Some(ida_mem), 0, line!(), "IDAGetAdjIDABmem", file!(), MSGAM_BAD_WHICH);
        return None;
    }

    /* Find the IDABMem entry in the linked list corresponding to 'which'. */
    let idaadj_mem = ida_mem.ida_adj_mem.as_mut().unwrap();
    let idx = idaadj_mem.IDAB_mem.iter().position(|b| b.ida_index == which).unwrap();
    Some(&mut idaadj_mem.IDAB_mem[idx].IDA_mem)
}

/*
 * IDAadjCheckPointRec (idas.h)
 *
 * my_addr/next_addr: the C void* checkpoint addresses become indices
 * into IDAadjMem.ck_mem (None = NULL).
 */
#[derive(Default, Clone, Copy)]
pub struct IDAadjCheckPointRec {
    pub my_addr: Option<usize>,
    pub next_addr: Option<usize>,
    pub t0: f64,
    pub t1: f64,
    pub nstep: i64,
    pub order: i32,
    pub step: f64,
}

/*
 * IDAGetAdjCheckPointsInfo
 *
 * Loads an array of nckpnts structures of type IDAadjCheckPointRec.
 * The user must allocate space for ckpnt (ncheck+1).
 * (The C walks the linked list newest-first; the Vec is walked in
 * reverse to match: ckpnt[0] = most recent check point.)
 */
pub fn IDAGetAdjCheckPointsInfo(ida_mem: &mut IDAMem, ckpnt: &mut [IDAadjCheckPointRec]) -> i32 {
    /* Is ASA initialized? */
    if !ida_mem.ida_adjMallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_ADJ, line!(), "IDAGetAdjCheckPointsInfo", file!(),
                        MSGAM_NO_ADJ);
        return IDA_NO_ADJ;
    }
    let idaadj_mem = ida_mem.ida_adj_mem.as_ref().unwrap();

    for (i, ck_idx) in (0..idaadj_mem.ck_mem.len()).rev().enumerate() {
        let ck_mem = &idaadj_mem.ck_mem[ck_idx];
        ckpnt[i].my_addr = Some(ck_idx);
        ckpnt[i].next_addr = if ck_idx == 0 { None } else { Some(ck_idx - 1) };
        ckpnt[i].t0 = ck_mem.ck_t0;
        ckpnt[i].t1 = ck_mem.ck_t1;
        ckpnt[i].nstep = ck_mem.ck_nst;
        ckpnt[i].order = ck_mem.ck_kk;
        ckpnt[i].step = ck_mem.ck_hh;
    }

    IDA_SUCCESS
}

/* IDAGetConsistentICB
 *
 * Returns the consistent initial conditions computed by IDACalcICB or
 * IDACalcICBS
 *
 * It must be preceded by a successful call to IDACalcICB or
 * IDACalcICBS for 'which' backward problem.
 */
pub fn IDAGetConsistentICB(
    ida_mem: &mut IDAMem,
    which: i32,
    yyB0_mod: Option<&mut NVector>,
    ypB0_mod: Option<&mut NVector>,
) -> i32 {
    let idx = match idaa_io_which_index(ida_mem, which, "IDAGetConsistentICB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    IDAGetConsistentIC(
        &mut ida_mem.ida_adj_mem.as_mut().unwrap().IDAB_mem[idx].IDA_mem,
        yyB0_mod,
        ypB0_mod,
    )
}

/* (C: int IDAGetUserDataB(ida_mem, which, void** user_dataB); the
   pointer out-parameter becomes a borrow of the stored user data, like
   the cvodea_io.rs precedent.) */
pub fn IDAGetUserDataB(ida_mem: &mut IDAMem, which: i32) -> Result<&mut UserData, i32> {
    let idx = idaa_io_which_index(ida_mem, which, "IDAGetUserDataB")?;

    Ok(&mut ida_mem.ida_adj_mem.as_mut().unwrap().IDAB_mem[idx].ida_user_data)
}

/*
 * -----------------------------------------------------------------
 * Undocumented development user-callable functions
 * -----------------------------------------------------------------
 */

/*
 * IDAGetAdjDataPointHermite
 *
 * Returns the 2 vectors stored for cubic Hermite interpolation at
 * the data point 'which'. The user must allocate space for yy and
 * yd.  (C NULL yy/yd outputs become None.)
 *
 * Returns IDA_MEM_NULL if ida_mem is NULL, IDA_ILL_INPUT if the
 * interpolation type previously specified is not IDA_HERMITE or
 * IDA_SUCCESS otherwise.
 */
pub fn IDAGetAdjDataPointHermite(
    ida_mem: &mut IDAMem,
    which: i32,
    t: &mut f64,
    yy: Option<&mut NVector>,
    yd: Option<&mut NVector>,
) -> i32 {
    /* Is ASA initialized? */
    if !ida_mem.ida_adjMallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_ADJ, line!(), "IDAGetAdjDataPointHermite", file!(),
                        MSGAM_NO_ADJ);
        return IDA_NO_ADJ;
    }
    let idaadj_mem = ida_mem.ida_adj_mem.as_ref().unwrap();

    if idaadj_mem.ia_interpType != IDA_HERMITE {
        IDAProcessError(None, IDA_ILL_INPUT, line!(), "IDAGetAdjDataPointHermite", file!(),
                        MSGAM_WRONG_INTERP);
        return IDA_ILL_INPUT;
    }

    let dt = &idaadj_mem.dt_mem[which as usize];
    *t = dt.t;

    if let DtpntContent::Hermite { y, yd: ydc, .. } = &dt.content {
        if let Some(yy) = yy {
            N_VScale(ONE, y, yy);
        }
        if let Some(yd) = yd {
            N_VScale(ONE, ydc, yd);
        }
    }

    IDA_SUCCESS
}

/*
 * IDAGetAdjDataPointPolynomial
 *
 * Returns the vector stored for polynomial interpolation at the
 * data point 'which'. The user must allocate space for y.
 *
 * Returns IDA_MEM_NULL if ida_mem is NULL, IDA_ILL_INPUT if the
 * interpolation type previously specified is not IDA_POLYNOMIAL or
 * IDA_SUCCESS otherwise.
 */
pub fn IDAGetAdjDataPointPolynomial(
    ida_mem: &mut IDAMem,
    which: i32,
    t: &mut f64,
    order: &mut i32,
    y: Option<&mut NVector>,
) -> i32 {
    /* Is ASA initialized? */
    if !ida_mem.ida_adjMallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_ADJ, line!(), "IDAGetAdjDataPointPolynomial",
                        file!(), MSGAM_NO_ADJ);
        return IDA_NO_ADJ;
    }
    let idaadj_mem = ida_mem.ida_adj_mem.as_ref().unwrap();

    if idaadj_mem.ia_interpType != IDA_POLYNOMIAL {
        IDAProcessError(None, IDA_ILL_INPUT, line!(), "IDAGetAdjDataPointPolynomial", file!(),
                        MSGAM_WRONG_INTERP);
        return IDA_ILL_INPUT;
    }

    let dt = &idaadj_mem.dt_mem[which as usize];
    *t = dt.t;

    if let DtpntContent::Polynomial { y: yc, order: ord, .. } = &dt.content {
        if let Some(y) = y {
            N_VScale(ONE, yc, y);
        }
        *order = *ord;
    }

    IDA_SUCCESS
}

/*
 * IDAGetAdjCurrentCheckPoint
 *
 * Returns the address of the 'active' check point.  (C void* out
 * becomes the Vec index Option.)
 */
pub fn IDAGetAdjCurrentCheckPoint(ida_mem: &mut IDAMem, addr: &mut Option<usize>) -> i32 {
    /* Is ASA initialized? */
    if !ida_mem.ida_adjMallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_ADJ, line!(), "IDAGetAdjCurrentCheckPoint", file!(),
                        MSGAM_NO_ADJ);
        return IDA_NO_ADJ;
    }

    *addr = ida_mem.ida_adj_mem.as_ref().unwrap().ia_ckpntData;

    IDA_SUCCESS
}

/*===============================================================
  Tests
  ===============================================================*/
#[cfg(test)]
mod tests {
    use super::*;
    use crate::idaa::{IDAAdjInit, IDACreateB};

    fn make_ida_mem(n: usize) -> IDAMem {
        let mut ida_mem = IDAMem::default();
        ida_mem.ida_tempv1 = NVector::new(n);
        ida_mem
    }

    /* All ***B optional-input wrappers require ASA init and a valid
       which, then delegate to the nested solver's setters (idaa_io.c
       lines 60-561). */
    #[test]
    fn adjoint_setters_delegate_to_nested_solver() {
        let mut ida_mem = make_ida_mem(2);

        assert_eq!(IDAAdjSetNoSensi(&mut ida_mem), IDA_NO_ADJ);
        assert_eq!(IDASetMaxOrdB(&mut ida_mem, 0, 3), IDA_NO_ADJ);

        assert_eq!(IDAAdjInit(&mut ida_mem, 4, IDA_POLYNOMIAL), IDA_SUCCESS);
        let mut which = -1;
        assert_eq!(IDACreateB(&mut ida_mem, &mut which), IDA_SUCCESS);

        assert_eq!(IDAAdjSetNoSensi(&mut ida_mem), IDA_SUCCESS);
        assert!(!ida_mem.ida_adj_mem.as_ref().unwrap().ia_storeSensi);

        assert_eq!(IDASetMaxOrdB(&mut ida_mem, 5, 3), IDA_ILL_INPUT); /* bad which */

        assert_eq!(IDASetMaxOrdB(&mut ida_mem, which, 3), IDA_SUCCESS);
        assert_eq!(IDASetMaxNumStepsB(&mut ida_mem, which, 250), IDA_SUCCESS);
        assert_eq!(IDASetInitStepB(&mut ida_mem, which, 1.0e-3), IDA_SUCCESS);
        assert_eq!(IDASetMaxStepB(&mut ida_mem, which, 10.0), IDA_SUCCESS);
        assert_eq!(IDASetSuppressAlgB(&mut ida_mem, which, true), IDA_SUCCESS);
        {
            let nested = &ida_mem.ida_adj_mem.as_ref().unwrap().IDAB_mem[0].IDA_mem;
            assert_eq!(nested.ida_maxord, 3);
            assert_eq!(nested.ida_mxstep, 250);
            assert_eq!(nested.ida_hin, 1.0e-3);
            assert!(nested.ida_suppressalg);
        }

        /* user data attach + borrow-back */
        assert_eq!(IDASetUserDataB(&mut ida_mem, which, Some(Box::new(7i32))), IDA_SUCCESS);
        let ud = IDAGetUserDataB(&mut ida_mem, which).unwrap();
        assert_eq!(*ud.as_ref().unwrap().downcast_ref::<i32>().unwrap(), 7);

        /* nested solver memory access */
        assert!(IDAGetAdjIDABmem(&mut ida_mem, which).is_some());
        assert!(IDAGetAdjIDABmem(&mut ida_mem, 9).is_none());
    }

    /* Check point info: index-based my_addr/next_addr, newest first;
       current check point mirrors ia_ckpntData. */
    #[test]
    fn checkpoint_info_and_current() {
        let mut ida_mem = make_ida_mem(2);
        assert_eq!(IDAAdjInit(&mut ida_mem, 4, IDA_HERMITE), IDA_SUCCESS);

        let mut addr = Some(99);
        assert_eq!(IDAGetAdjCurrentCheckPoint(&mut ida_mem, &mut addr), IDA_SUCCESS);
        assert_eq!(addr, None);

        /* wrong-interpolation guard for the data-point getters */
        let mut t = 0.0;
        let mut order = 0;
        assert_eq!(IDAGetAdjDataPointPolynomial(&mut ida_mem, 0, &mut t, &mut order, None),
                   IDA_ILL_INPUT);

        let mut ckpnt = [IDAadjCheckPointRec::default(); 1];
        assert_eq!(IDAGetAdjCheckPointsInfo(&mut ida_mem, &mut ckpnt), IDA_SUCCESS);
    }
}
