/* -----------------------------------------------------------------
 * Translated from src/cvodes/cvodea.c (CVODES 7.7.0).
 *
 * Adjoint sensitivity analysis: checkpointing (CVodeF), backward
 * problems (CVodeCreateB/CVodeInitB/CVodeB, ...), and the cubic
 * Hermite / variable-order polynomial interpolation modules.
 *
 * Port notes (pinned decisions, see cvodes_impl.rs):
 *  - The C check point linked list (head = MOST RECENT check point,
 *    ck_next walks back in time, tail = initial check point) becomes
 *    Vec<CVckpntMem> with index 0 = initial check point and
 *    last() = most recent; "ck_next == NULL" <=> index == 0, and
 *    "ck_mem = ck_mem->ck_next" <=> index -= 1.
 *  - The C backward-problem list (C prepends: head = newest) becomes
 *    Vec<CVodeBMem> in creation order, so cv_index == position.
 *    Loops "over all backward problems" therefore run in creation
 *    order instead of the C's reverse order; the problems are
 *    mutually independent so the results are unchanged.
 *  - The C interpolation-module function pointers ca_IMmalloc /
 *    ca_IMfree / ca_IMstore / ca_IMget are replaced by dispatch on
 *    ca_IMtype (cvaIMstore_dispatch / cvaIMget_dispatch). The free
 *    routines are absorbed by RAII (CVodeAdjFree drops CVadjMem).
 *  - In C, ca_Y[i]/ca_YS[i] are aliases of zn[i]/znS[i] used as
 *    interpolation scratch; here they are owned workspace vectors
 *    allocated in CVodeF (every use overwrites them before reading,
 *    so copying semantics are identical).
 *  - In C, the backward problem's CVODES memory has user_data ==
 *    the FORWARD cvode_mem (set once in CVodeCreateB). Here that
 *    self-reference is created transiently: CVodeB moves the outer
 *    CVodeMem into the nested problem's cv_user_data around each
 *    CVode() call on the backward problem (and CVArhs/CVArhsQ and
 *    the cvodes_ls.rs B-wrappers downcast it back out).
 *  - An empty Vec<NVector> plays the role of the C NULL N_Vector*
 *    argument of the IMget routines.
 * -----------------------------------------------------------------*/

use crate::cvodes::{
    CVode, CVodeCreate, CVodeGetDky, CVodeGetQuad, CVodeInit, CVodeQuadInit, CVodeQuadReInit,
    CVodeQuadSStolerances, CVodeQuadSVtolerances, CVodeQuadSensReInit, CVodeReInit,
    CVodeSStolerances, CVodeSVtolerances, CVodeSensReInit, cvSensRhsWrapper,
};
use crate::cvodes_impl::*;
use crate::cvodes_io::{CVodeSetInitStep, CVodeSetMaxHnilWarns, CVodeSetStopTime};
use crate::nvector_serial::*;
use crate::sundials_math::*;
use crate::sundials_types::*;

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;
const HUNDRED: f64 = 100.0;

/* fuzz factor for IMget */
const FUZZ_FACTOR: f64 = 1000000.0;

/* Helper: mutable access to the (attached) adjoint memory. */
fn adj(cv_mem: &mut CVodeMem) -> &mut CVadjMem {
    cv_mem.cv_adj_mem.as_mut().unwrap()
}

/*
 * =================================================================
 * EXPORTED FUNCTIONS IMPLEMENTATION
 * =================================================================
 */

/*
 * CVodeAdjInit
 *
 * This routine initializes ASA and allocates space for the adjoint
 * memory structure.
 */
pub fn CVodeAdjInit(cv_mem: &mut CVodeMem, steps: i64, interp: i32) -> i32 {
    if steps <= 0 {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeAdjInit", file!(),
                       MSGCV_BAD_STEPS);
        return CV_ILL_INPUT;
    }

    if interp != CV_HERMITE && interp != CV_POLYNOMIAL {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeAdjInit", file!(),
                       MSGCV_BAD_INTERP);
        return CV_ILL_INPUT;
    }

    /* ----------------------------
     * Allocate CVODEA memory block
     * ---------------------------- */

    /* Allocate space for the array of Data Point structures */
    let mut dt_mem = Vec::with_capacity((steps + 1) as usize);
    for _ in 0..=steps {
        dt_mem.push(CVdtpntMem { t: ZERO, content: DtpntContent::None });
    }

    let ca_mem = CVadjMem {
        /* Initialization of check points */
        ck_mem: Vec::new(),
        ca_nckpnts: 0,
        ca_ckpntData: None,

        /* Initialization of interpolation data */
        ca_IMtype: interp,
        ca_nsteps: steps,
        ca_ilast: -1, /* invalid value */
        dt_mem,
        ca_np: 0,

        /* (The interpolation-module function pointers of the C become
           dispatch on ca_IMtype; see cvaIMstore_dispatch/cvaIMget_dispatch.) */

        /* The interpolation module has not been initialized yet */
        ca_IMmallocDone: SUNFALSE,

        /* By default we will store but not interpolate sensitivities
         *  - IMstoreSensi will be set in CVodeF to SUNFALSE if FSA is not
         *    enabled or if the user can force this through CVodeSetAdjNoSensi
         *  - IMinterpSensi will be set in CVodeB to SUNTRUE if IMstoreSensi
         *    is SUNTRUE and if at least one backward problem requires
         *    sensitivities */
        ca_IMstoreSensi: SUNTRUE,
        ca_IMinterpSensi: SUNFALSE,
        ca_IMnewData: SUNFALSE,

        /* Initialize list of backward problems */
        cvB_mem: Vec::new(),
        ca_bckpbCrt: None,
        ca_nbckpbs: 0,

        /* CVodeF and CVodeB not called yet */
        ca_firstCVodeFcall: SUNTRUE,
        ca_tstopCVodeFcall: SUNFALSE,
        ca_tstopCVodeF: ZERO,
        ca_firstCVodeBcall: SUNTRUE,
        ca_rootret: SUNFALSE,
        ca_troot: ZERO,

        ca_tinitial: ZERO,
        ca_tfinal: ZERO,

        /* Interpolation workspace (allocated in CVodeF) */
        ca_Y: Vec::new(),
        ca_YS: std::array::from_fn(|_| Vec::new()),
        ca_T: [ZERO; L_MAX],

        /* Workspace for wrapper functions (allocated by IMmalloc) */
        ca_ytmp: NVector::default(),
        ca_yStmp: Vec::new(),
    };

    /* Attach ca_mem to CVodeMem structure */
    cv_mem.cv_adj_mem = Some(Box::new(ca_mem));

    /* ASA initialized and allocated */
    cv_mem.cv_adj = SUNTRUE;
    cv_mem.cv_adjMallocDone = SUNTRUE;

    CV_SUCCESS
}

/* CVodeAdjReInit
 *
 * This routine reinitializes the CVODEA memory structure assuming that the
 * the number of steps between check points and the type of interpolation
 * remain unchanged.
 * The list of check points (and associated memory) is deleted.
 * The list of backward problems is kept (however, new backward problems can
 * be added to this list by calling CVodeCreateB).
 * The CVODES memory for the forward and backward problems can be reinitialized
 * separately by calling CVodeReInit and CVodeReInitB, respectively.
 * NOTE: if a completely new list of backward problems is also needed, then
 *       simply free the adjoint memory (by calling CVodeAdjFree) and
 *       reinitialize ASA with CVodeAdjInit.
 */
pub fn CVodeAdjReInit(cv_mem: &mut CVodeMem) -> i32 {
    /* Was ASA initialized? */
    if !cv_mem.cv_adjMallocDone {
        cvProcessError(Some(cv_mem), CV_NO_ADJ, line!(), "CVodeAdjReInit", file!(), MSGCV_NO_ADJ);
        return CV_NO_ADJ;
    }

    let ca_mem = adj(cv_mem);

    /* Free current list of Check Points */
    ca_mem.ck_mem.clear();

    /* Initialization of check points */
    ca_mem.ca_nckpnts = 0;
    ca_mem.ca_ckpntData = None;

    /* CVodeF and CVodeB not called yet */
    ca_mem.ca_firstCVodeFcall = SUNTRUE;
    ca_mem.ca_tstopCVodeFcall = SUNFALSE;
    ca_mem.ca_firstCVodeBcall = SUNTRUE;

    CV_SUCCESS
}

/*
 * CVodeAdjFree
 *
 * This routine frees the memory allocated by CVodeAdjInit.
 * (Check points, interpolation data, and backward problems are all
 * owned by CVadjMem, so dropping it frees everything - the explicit
 * CVAckpntDelete / IMfree / CVAbckpbDelete loops of the C collapse
 * into the drop.)
 */
pub fn CVodeAdjFree(cv_mem: &mut CVodeMem) {
    if cv_mem.cv_adjMallocDone {
        cv_mem.cv_adj_mem = None;
        cv_mem.cv_adjMallocDone = SUNFALSE;
    }
}

/*
 * CVodeF
 *
 * This routine integrates to tout and returns solution into yout.
 * In the same time, it stores check point data every 'steps' steps.
 *
 * CVodeF can be called repeatedly by the user.
 *
 * ncheckPtr points to the number of check points stored so far.
 */
pub fn CVodeF(
    cv_mem: &mut CVodeMem,
    tout: f64,
    yout: &mut NVector,
    tret: &mut f64,
    itask: i32,
    ncheckPtr: &mut i32,
) -> i32 {
    /* Was ASA initialized? */
    if !cv_mem.cv_adjMallocDone {
        cvProcessError(Some(cv_mem), CV_NO_ADJ, line!(), "CVodeF", file!(), MSGCV_NO_ADJ);
        return CV_NO_ADJ;
    }

    /* (yout and tret cannot be NULL in the Rust port) */

    /* Check for valid itask */
    if itask != CV_NORMAL && itask != CV_ONE_STEP {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeF", file!(), MSGCV_BAD_ITASK);
        return CV_ILL_INPUT;
    }

    /* All error checking done */

    /* If tstop is enabled, store some info */
    if cv_mem.cv_tstopset {
        let tstop = cv_mem.cv_tstop;
        let ca_mem = adj(cv_mem);
        ca_mem.ca_tstopCVodeFcall = SUNTRUE;
        ca_mem.ca_tstopCVodeF = tstop;
    }

    /* On the first step:
     *   - set tinitial
     *   - initialize list of check points
     *   - if needed, initialize the interpolation module
     *   - load dt_mem[0]
     * On subsequent steps, test if taking a new step is necessary.
     */
    let mut flag;
    if adj(cv_mem).ca_firstCVodeFcall {
        adj(cv_mem).ca_tinitial = cv_mem.cv_tn;

        let ck = CVAckpntInit(cv_mem);
        {
            let ca_mem = adj(cv_mem);
            ca_mem.ck_mem.push(ck);
        }

        if !adj(cv_mem).ca_IMmallocDone {
            /* Do we need to store sensitivities? */
            if !cv_mem.cv_sensi {
                adj(cv_mem).ca_IMstoreSensi = SUNFALSE;
            }

            /* Allocate space for interpolation data */
            match adj(cv_mem).ca_IMtype {
                CV_HERMITE => cvaHermiteMalloc(cv_mem),
                _ => cvaPolynomialMalloc(cv_mem),
            }

            /* Rename zn and, if needed, znS for use in interpolation:
               here, allocate the owned scratch vectors ca_Y / ca_YS
               (the C aliases zn[i]/znS[i]; every use overwrites them
               before reading, so owned scratch is equivalent). */
            {
                let tmpl = cv_mem.cv_tempv.clone();
                let ns = cv_mem.cv_Ns as usize;
                let store_sensi = adj(cv_mem).ca_IMstoreSensi;
                let ca_mem = adj(cv_mem);
                ca_mem.ca_Y = (0..L_MAX).map(|_| tmpl.clone()).collect();
                if store_sensi {
                    for i in 0..L_MAX {
                        ca_mem.ca_YS[i] = (0..ns).map(|_| tmpl.clone()).collect();
                    }
                }
            }

            adj(cv_mem).ca_IMmallocDone = SUNTRUE;
        }

        {
            let ca_mem = adj(cv_mem);
            ca_mem.dt_mem[0].t = ca_mem.ck_mem.last().unwrap().ck_t0;
        }
        cvaIMstore_dispatch(cv_mem, 0);

        adj(cv_mem).ca_firstCVodeFcall = SUNFALSE;
    } else if itask == CV_NORMAL {
        /* When in normal mode, check if tout was passed or if a previous root
           was not reported and return an interpolated solution. No changes to
           ck_mem or dt_mem are needed. */

        /* flag to signal if an early return is needed */
        let mut earlyret = SUNFALSE;
        flag = CV_SUCCESS;

        /* if a root needs to be reported compare tout to troot otherwise
           compare to the current time tn */
        let ttest = if adj(cv_mem).ca_rootret {
            adj(cv_mem).ca_troot
        } else {
            cv_mem.cv_tn
        };

        if (ttest - tout) * cv_mem.cv_h >= ZERO {
            /* ttest is after tout, interpolate to tout */
            *tret = tout;
            flag = CVodeGetDky(cv_mem, tout, 0, yout);
            earlyret = SUNTRUE;
        } else if adj(cv_mem).ca_rootret {
            /* tout is after troot, interpolate to troot */
            let troot = adj(cv_mem).ca_troot;
            *tret = troot;
            let _ = CVodeGetDky(cv_mem, troot, 0, yout);
            flag = CV_ROOT_RETURN;
            adj(cv_mem).ca_rootret = SUNFALSE;
            earlyret = SUNTRUE;
        }

        /* return if necessary */
        if earlyret {
            let nst = cv_mem.cv_nst;
            let ca_mem = adj(cv_mem);
            *ncheckPtr = ca_mem.ca_nckpnts;
            ca_mem.ca_IMnewData = SUNTRUE;
            ca_mem.ca_ckpntData = Some(ca_mem.ck_mem.len() - 1);
            ca_mem.ca_np = nst % ca_mem.ca_nsteps + 1;
            return flag;
        }
    }

    /* Integrate to tout (in CV_ONE_STEP mode) while loading check points */
    let mut nstloc: i64 = 0;
    loop {
        /* Check for too many steps */
        if cv_mem.cv_mxstep > 0 && nstloc >= cv_mem.cv_mxstep {
            cvProcessError(Some(cv_mem), CV_TOO_MUCH_WORK, line!(), "CVodeF", file!(),
                &format!("At t = {}, mxstep steps taken before reaching tout.", cv_mem.cv_tn));
            flag = CV_TOO_MUCH_WORK;
            break;
        }

        /* Perform one step of the integration */
        flag = CVode(cv_mem, tout, yout, tret, CV_ONE_STEP);
        if flag < 0 {
            break;
        }

        nstloc += 1;

        /* Test if a new check point is needed */
        if cv_mem.cv_nst % adj(cv_mem).ca_nsteps == 0 {
            let tn = cv_mem.cv_tn;
            adj(cv_mem).ck_mem.last_mut().unwrap().ck_t1 = tn;

            /* Create a new check point, load it, and append it to the list
               (C prepends to the linked list; Vec push = newest last) */
            let tmp = CVAckpntNew(cv_mem);
            {
                let ca_mem = adj(cv_mem);
                ca_mem.ck_mem.push(tmp);
                ca_mem.ca_nckpnts += 1;
            }
            cv_mem.cv_forceSetup = SUNTRUE;

            /* Reset i=0 and load dt_mem[0] */
            {
                let ca_mem = adj(cv_mem);
                ca_mem.dt_mem[0].t = ca_mem.ck_mem.last().unwrap().ck_t0;
            }
            cvaIMstore_dispatch(cv_mem, 0);
        } else {
            /* Load next point in dt_mem */
            let idx = (cv_mem.cv_nst % adj(cv_mem).ca_nsteps) as usize;
            let tn = cv_mem.cv_tn;
            adj(cv_mem).dt_mem[idx].t = tn;
            cvaIMstore_dispatch(cv_mem, idx);
        }

        /* Set t1 field of the current check point structure
           for the case in which there will be no future
           check points */
        let tn = cv_mem.cv_tn;
        adj(cv_mem).ck_mem.last_mut().unwrap().ck_t1 = tn;

        /* tfinal is now set to tn */
        adj(cv_mem).ca_tfinal = tn;

        /* Return if in CV_ONE_STEP mode */
        if itask == CV_ONE_STEP {
            break;
        }

        /* CV_NORMAL_STEP returns */

        /* Return if tout reached */
        if (*tret - tout) * cv_mem.cv_h >= ZERO {
            /* If this was a root return, save the root time to return later */
            if flag == CV_ROOT_RETURN {
                let ca_mem = adj(cv_mem);
                ca_mem.ca_rootret = SUNTRUE;
                ca_mem.ca_troot = *tret;
            }

            /* Get solution value at tout to return now */
            *tret = tout;
            flag = CVodeGetDky(cv_mem, tout, 0, yout);

            /* Reset tretlast in cv_mem so that CVodeGetQuad and CVodeGetSens
             * evaluate quadratures and/or sensitivities at the proper time */
            cv_mem.cv_tretlast = tout;

            break;
        }

        /* Return if tstop or a root was found */
        if flag == CV_TSTOP_RETURN || flag == CV_ROOT_RETURN {
            break;
        }
    } /* end of loop() */

    /* Get ncheck from ca_mem */
    let nst = cv_mem.cv_nst;
    let ca_mem = adj(cv_mem);
    *ncheckPtr = ca_mem.ca_nckpnts;

    /* Data is available for the last interval */
    ca_mem.ca_IMnewData = SUNTRUE;
    ca_mem.ca_ckpntData = Some(ca_mem.ck_mem.len() - 1);
    ca_mem.ca_np = nst % ca_mem.ca_nsteps + 1;

    flag
}

/*
 * =================================================================
 * FUNCTIONS FOR BACKWARD PROBLEMS
 * =================================================================
 */

/* Shared preamble of every ***B function: check ASA initialized and
   `which`, then return the index of the matching backward problem
   (C: linked-list search on cv_index). */
fn cva_which_index(cv_mem: &mut CVodeMem, which: i32, fname: &str) -> Result<usize, i32> {
    if !cv_mem.cv_adjMallocDone {
        cvProcessError(Some(cv_mem), CV_NO_ADJ, line!(), fname, file!(), MSGCV_NO_ADJ);
        return Err(CV_NO_ADJ);
    }

    if which >= cv_mem.cv_adj_mem.as_ref().unwrap().ca_nbckpbs {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), fname, file!(), MSGCV_BAD_WHICH);
        return Err(CV_ILL_INPUT);
    }

    /* Find the CVodeBMem entry corresponding to which */
    let ca_mem = cv_mem.cv_adj_mem.as_ref().unwrap();
    Ok(ca_mem
        .cvB_mem
        .iter()
        .position(|b| b.cv_index == which)
        .unwrap())
}

pub fn CVodeCreateB(cv_mem: &mut CVodeMem, lmmB: i32, which: &mut i32) -> i32 {
    /* Was ASA initialized? */
    if !cv_mem.cv_adjMallocDone {
        cvProcessError(Some(cv_mem), CV_NO_ADJ, line!(), "CVodeCreateB", file!(), MSGCV_NO_ADJ);
        return CV_NO_ADJ;
    }

    /* Create and set a new CVODES object for the backward problem
       (C returns CV_MEM_FAIL if CVodeCreate fails on an illegal lmmB;
       the Rust CVodeCreate panics instead, so pre-check here). */
    if lmmB != CV_ADAMS && lmmB != CV_BDF {
        cvProcessError(Some(cv_mem), CV_MEM_FAIL, line!(), "CVodeCreateB", file!(), MSGCV_MEM_FAIL);
        return CV_MEM_FAIL;
    }

    let sunctx = cv_mem.cv_sunctx.clone();
    let mut cvodeB_mem = CVodeCreate(lmmB, &sunctx);

    /* We need to ensure Ns is set in the new CVODES object so that Ns is
       accessible in callbacks which only have access to cvodeB_mem */
    cvodeB_mem.cv_Ns = cv_mem.cv_Ns;

    /* C: CVodeSetUserData(cvodeB_mem, cvode_mem) - the forward memory is
       installed in the backward problem's user_data transiently by CVodeB
       in this port (see the ownership note at the top of this file). */

    CVodeSetMaxHnilWarns(&mut cvodeB_mem, -1);

    /* Set/initialize fields in the new CVodeBMem object, new_cvB_mem */
    let ca_mem = adj(cv_mem);

    let new_cvB_mem = CVodeBMem {
        cv_index: ca_mem.ca_nbckpbs,
        cv_t0: ZERO,
        cv_mem: cvodeB_mem,
        cv_f: None,
        cv_fs: None,
        cv_fQ: None,
        cv_fQs: None,
        cv_user_data: None,
        cv_lmem: None,
        cv_pmem: None,
        cv_y: NVector::default(),
        cv_tout: ZERO,
        cv_f_withSensi: SUNFALSE,
        cv_fQ_withSensi: SUNFALSE,
    };

    /* Attach the new object to the list cvB_mem (C prepends; the Vec
       appends, keeping cv_index == position) */
    ca_mem.cvB_mem.push(new_cvB_mem);

    /* Return the index of the newly created CVodeBMem object.
     * This must be passed to CVodeInitB and to other ***B
     * functions to set optional inputs for this backward problem */
    *which = ca_mem.ca_nbckpbs;

    ca_mem.ca_nbckpbs += 1;

    CV_SUCCESS
}

pub fn CVodeInitB(cv_mem: &mut CVodeMem, which: i32, fB: CVRhsFnB, tB0: f64, yB0: &NVector) -> i32 {
    let idx = match cva_which_index(cv_mem, which, "CVodeInitB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    /* Allocate and set the CVODES object */
    let cvB_mem = &mut adj(cv_mem).cvB_mem[idx];
    let flag = CVodeInit(&mut cvB_mem.cv_mem, CVArhs, tB0, yB0);

    if flag != CV_SUCCESS {
        return flag;
    }

    /* Copy fB function in cvB_mem */
    cvB_mem.cv_f_withSensi = SUNFALSE;
    cvB_mem.cv_f = Some(fB);

    /* Allocate space and initialize the y Nvector in cvB_mem */
    cvB_mem.cv_t0 = tB0;
    cvB_mem.cv_y = yB0.clone();

    CV_SUCCESS
}

pub fn CVodeInitBS(
    cv_mem: &mut CVodeMem,
    which: i32,
    fBs: CVRhsFnBS,
    tB0: f64,
    yB0: &NVector,
) -> i32 {
    let idx = match cva_which_index(cv_mem, which, "CVodeInitBS") {
        Ok(i) => i,
        Err(e) => return e,
    };

    /* Allocate and set the CVODES object */
    let cvB_mem = &mut adj(cv_mem).cvB_mem[idx];
    let flag = CVodeInit(&mut cvB_mem.cv_mem, CVArhs, tB0, yB0);

    if flag != CV_SUCCESS {
        return flag;
    }

    /* Copy fBs function in cvB_mem */
    cvB_mem.cv_f_withSensi = SUNTRUE;
    cvB_mem.cv_fs = Some(fBs);

    /* Allocate space and initialize the y Nvector in cvB_mem */
    cvB_mem.cv_t0 = tB0;
    cvB_mem.cv_y = yB0.clone();

    CV_SUCCESS
}

pub fn CVodeReInitB(cv_mem: &mut CVodeMem, which: i32, tB0: f64, yB0: &NVector) -> i32 {
    let idx = match cva_which_index(cv_mem, which, "CVodeReInitB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    /* Reinitialize CVODES object */
    CVodeReInit(&mut adj(cv_mem).cvB_mem[idx].cv_mem, tB0, yB0)
}

pub fn CVodeSStolerancesB(cv_mem: &mut CVodeMem, which: i32, reltolB: f64, abstolB: f64) -> i32 {
    let idx = match cva_which_index(cv_mem, which, "CVodeSStolerancesB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    /* Set tolerances */
    CVodeSStolerances(&mut adj(cv_mem).cvB_mem[idx].cv_mem, reltolB, abstolB)
}

pub fn CVodeSVtolerancesB(
    cv_mem: &mut CVodeMem,
    which: i32,
    reltolB: f64,
    abstolB: &NVector,
) -> i32 {
    let idx = match cva_which_index(cv_mem, which, "CVodeSVtolerancesB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    /* Set tolerances */
    CVodeSVtolerances(&mut adj(cv_mem).cvB_mem[idx].cv_mem, reltolB, abstolB)
}

pub fn CVodeQuadInitB(cv_mem: &mut CVodeMem, which: i32, fQB: CVQuadRhsFnB, yQB0: &NVector) -> i32 {
    let idx = match cva_which_index(cv_mem, which, "CVodeQuadInitB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    let cvB_mem = &mut adj(cv_mem).cvB_mem[idx];
    let flag = CVodeQuadInit(&mut cvB_mem.cv_mem, CVArhsQ, yQB0);
    if flag != CV_SUCCESS {
        return flag;
    }

    cvB_mem.cv_fQ_withSensi = SUNFALSE;
    cvB_mem.cv_fQ = Some(fQB);

    CV_SUCCESS
}

pub fn CVodeQuadInitBS(
    cv_mem: &mut CVodeMem,
    which: i32,
    fQBs: CVQuadRhsFnBS,
    yQB0: &NVector,
) -> i32 {
    let idx = match cva_which_index(cv_mem, which, "CVodeQuadInitBS") {
        Ok(i) => i,
        Err(e) => return e,
    };

    let cvB_mem = &mut adj(cv_mem).cvB_mem[idx];
    let flag = CVodeQuadInit(&mut cvB_mem.cv_mem, CVArhsQ, yQB0);
    if flag != CV_SUCCESS {
        return flag;
    }

    cvB_mem.cv_fQ_withSensi = SUNTRUE;
    cvB_mem.cv_fQs = Some(fQBs);

    CV_SUCCESS
}

pub fn CVodeQuadReInitB(cv_mem: &mut CVodeMem, which: i32, yQB0: &NVector) -> i32 {
    let idx = match cva_which_index(cv_mem, which, "CVodeQuadReInitB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    let flag = CVodeQuadReInit(&mut adj(cv_mem).cvB_mem[idx].cv_mem, yQB0);
    if flag != CV_SUCCESS {
        return flag;
    }

    CV_SUCCESS
}

pub fn CVodeQuadSStolerancesB(
    cv_mem: &mut CVodeMem,
    which: i32,
    reltolQB: f64,
    abstolQB: f64,
) -> i32 {
    let idx = match cva_which_index(cv_mem, which, "CVodeQuadSStolerancesB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    CVodeQuadSStolerances(&mut adj(cv_mem).cvB_mem[idx].cv_mem, reltolQB, abstolQB)
}

pub fn CVodeQuadSVtolerancesB(
    cv_mem: &mut CVodeMem,
    which: i32,
    reltolQB: f64,
    abstolQB: &NVector,
) -> i32 {
    let idx = match cva_which_index(cv_mem, which, "CVodeQuadSVtolerancesB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    CVodeQuadSVtolerances(&mut adj(cv_mem).cvB_mem[idx].cv_mem, reltolQB, abstolQB)
}

/*
 * CVodeB
 *
 * This routine performs the backward integration towards tBout
 * of all backward problems that were defined.
 * When necessary, it performs a forward integration between two
 * consecutive check points to update interpolation data.
 *
 * On a successful return, CVodeB returns CV_SUCCESS.
 *
 * NOTE that CVodeB DOES NOT return the solution for the backward
 * problem(s). Use CVodeGetB to extract the solution at tBret
 * for any given backward problem.
 *
 * If there are multiple backward problems and multiple check points,
 * CVodeB may not succeed in getting all problems to take one step
 * when called in ONE_STEP mode.
 */
pub fn CVodeB(cv_mem: &mut CVodeMem, mut tBout: f64, itaskB: i32) -> i32 {
    /* Was ASA initialized? */
    if !cv_mem.cv_adjMallocDone {
        cvProcessError(Some(cv_mem), CV_NO_ADJ, line!(), "CVodeB", file!(), MSGCV_NO_ADJ);
        return CV_NO_ADJ;
    }

    /* Check if any backward problem has been defined */
    if cv_mem.cv_adj_mem.as_ref().unwrap().ca_nbckpbs == 0 {
        cvProcessError(Some(cv_mem), CV_NO_BCK, line!(), "CVodeB", file!(), MSGCV_NO_BCK);
        return CV_NO_BCK;
    }

    /* Check whether CVodeF has been called */
    if cv_mem.cv_adj_mem.as_ref().unwrap().ca_firstCVodeFcall {
        cvProcessError(Some(cv_mem), CV_NO_FWD, line!(), "CVodeB", file!(), MSGCV_NO_FWD);
        return CV_NO_FWD;
    }

    let (tinitial, tfinal) = {
        let ca_mem = cv_mem.cv_adj_mem.as_ref().unwrap();
        (ca_mem.ca_tinitial, ca_mem.ca_tfinal)
    };
    let sign: f64 = if tfinal - tinitial > ZERO { 1.0 } else { -1.0 };

    /* If this is the first call, loop over all backward problems and
     *   - check that tB0 is valid
     *   - check that tBout is ahead of tB0 in the backward direction
     *   - check whether we need to interpolate forward sensitivities
     */
    if cv_mem.cv_adj_mem.as_ref().unwrap().ca_firstCVodeBcall {
        let nb = cv_mem.cv_adj_mem.as_ref().unwrap().cvB_mem.len();
        for i in 0..nb {
            let (tBn, index, with_sensi) = {
                let b = &cv_mem.cv_adj_mem.as_ref().unwrap().cvB_mem[i];
                (b.cv_mem.cv_tn, b.cv_index, b.cv_f_withSensi || b.cv_fQ_withSensi)
            };

            if sign * (tBn - tinitial) < ZERO || sign * (tfinal - tBn) < ZERO {
                cvProcessError(Some(cv_mem), CV_BAD_TB0, line!(), "CVodeB", file!(),
                    &format!("The initial time tB0 for problem {} is outside the interval over which the forward problem was solved.", index));
                return CV_BAD_TB0;
            }

            if sign * (tBn - tBout) <= ZERO {
                cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeB", file!(),
                    MSGCV_BAD_TBOUT);
                return CV_ILL_INPUT;
            }

            if with_sensi {
                adj(cv_mem).ca_IMinterpSensi = SUNTRUE;
            }
        }

        {
            let ca_mem = cv_mem.cv_adj_mem.as_ref().unwrap();
            if ca_mem.ca_IMinterpSensi && !ca_mem.ca_IMstoreSensi {
                cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeB", file!(),
                    MSGCV_BAD_SENSI);
                return CV_ILL_INPUT;
            }
        }

        adj(cv_mem).ca_firstCVodeBcall = SUNFALSE;
    }

    /* Check if itaskB is legal */
    if itaskB != CV_NORMAL && itaskB != CV_ONE_STEP {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeB", file!(), MSGCV_BAD_ITASKB);
        return CV_ILL_INPUT;
    }

    /* Check if tBout is legal */
    if sign * (tBout - tinitial) < ZERO || sign * (tfinal - tBout) < ZERO {
        let tfuzz = HUNDRED * cv_mem.cv_uround * (SUNRabs(tinitial) + SUNRabs(tfinal));
        if sign * (tBout - tinitial) < ZERO && SUNRabs(tBout - tinitial) < tfuzz {
            tBout = tinitial;
        } else {
            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeB", file!(),
                MSGCV_BAD_TBOUT);
            return CV_ILL_INPUT;
        }
    }

    /* Loop through the check points and stop as soon as a backward
     * problem has its tn value behind the current check point's t0_
     * value (in the backward direction).
     * (C walks the linked list from newest to oldest; the Vec is walked
     * from last index down to 0.) */
    let mut ck_idx = cv_mem.cv_adj_mem.as_ref().unwrap().ck_mem.len() - 1;

    'search: loop {
        let ca_mem = cv_mem.cv_adj_mem.as_ref().unwrap();
        let ck_t0 = ca_mem.ck_mem[ck_idx].ck_t0;

        for b in ca_mem.cvB_mem.iter() {
            let tBn = b.cv_mem.cv_tn;

            if sign * (tBn - ck_t0) > ZERO {
                break 'search;
            }

            if itaskB == CV_NORMAL && tBn == ck_t0 && sign * (tBout - ck_t0) >= ZERO {
                break 'search;
            }
        }

        if ck_idx == 0 {
            /* ck_next == NULL */
            break;
        }

        ck_idx -= 1;
    }

    /* Starting with the current check point from above, loop over check
       points while propagating backward problems */
    let mut flag = 0;
    loop {
        /* Store interpolation data if not available.
           This is the 2nd forward integration pass */
        if cv_mem.cv_adj_mem.as_ref().unwrap().ca_ckpntData != Some(ck_idx) {
            flag = CVAdataStore(cv_mem, ck_idx);
            if flag != CV_SUCCESS {
                break;
            }
        }

        /* Loop through all backward problems and, if needed,
         * propagate their solution towards tBout */
        let nb = cv_mem.cv_adj_mem.as_ref().unwrap().cvB_mem.len();
        let mut err_index = 0;
        for i in 0..nb {
            let (tBn, ck_t0) = {
                let ca_mem = cv_mem.cv_adj_mem.as_ref().unwrap();
                (ca_mem.cvB_mem[i].cv_mem.cv_tn, ca_mem.ck_mem[ck_idx].ck_t0)
            };

            /* Decide if current backward problem is "active" in this check point */
            let mut is_active = SUNTRUE;

            if tBn == ck_t0 && sign * (tBout - ck_t0) < ZERO {
                is_active = SUNFALSE;
            }
            if tBn == ck_t0 && itaskB == CV_ONE_STEP {
                is_active = SUNFALSE;
            }
            if sign * (tBn - ck_t0) < ZERO {
                is_active = SUNFALSE;
            }

            if is_active {
                /* Store the index of the current backward problem memory
                 * in ca_mem to be used in the wrapper functions */
                adj(cv_mem).ca_bckpbCrt = Some(i);

                /* Integrate current backward problem */
                CVodeSetStopTime(&mut adj(cv_mem).cvB_mem[i].cv_mem, ck_t0);

                let mut tBret = ZERO;
                flag = cva_integrate_backward(cv_mem, i, tBout, itaskB, &mut tBret);

                /* Set the time at which we will report solution and/or
                   quadratures */
                adj(cv_mem).cvB_mem[i].cv_tout = tBret;

                /* If an error occurred, exit the loop */
                if flag < 0 {
                    err_index = adj(cv_mem).cvB_mem[i].cv_index;
                    break;
                }
            } else {
                flag = CV_SUCCESS;
                adj(cv_mem).cvB_mem[i].cv_tout = tBn;
            }
        }

        /* If an error occurred, return now */
        if flag < 0 {
            cvProcessError(Some(cv_mem), flag, line!(), "CVodeB", file!(),
                &format!("Error occurred while integrating backward problem # {}", err_index));
            return flag;
        }

        /* If in CV_ONE_STEP mode, return now (flag = CV_SUCCESS) */
        if itaskB == CV_ONE_STEP {
            break;
        }

        /* If all backward problems have successfully reached tBout, return now */
        let mut reached_tBout = SUNTRUE;
        {
            let ca_mem = cv_mem.cv_adj_mem.as_ref().unwrap();
            for b in ca_mem.cvB_mem.iter() {
                if sign * (b.cv_tout - tBout) > ZERO {
                    reached_tBout = SUNFALSE;
                    break;
                }
            }
        }

        if reached_tBout {
            break;
        }

        /* Move check point in linked list to next one */
        ck_idx -= 1;
    }

    flag
}

/* Integrate backward problem `idx` towards tBout: the ownership dance
   that realizes the C's permanent user_data == forward cvode_mem link
   (see the note at the top of this file). */
fn cva_integrate_backward(
    cv_mem: &mut CVodeMem,
    idx: usize,
    tBout: f64,
    itaskB: i32,
    tBret: &mut f64,
) -> i32 {
    let sunctx = cv_mem.cv_sunctx.clone();

    /* detach the backward problem's solver memory and workspace vector */
    let (mut nested, mut yB) = {
        let placeholder = CVodeCreate(CV_BDF, &sunctx);
        let ca_mem = adj(cv_mem);
        let nested = std::mem::replace(&mut ca_mem.cvB_mem[idx].cv_mem, placeholder);
        let y = std::mem::take(&mut ca_mem.cvB_mem[idx].cv_y);
        (nested, y)
    };

    /* move the forward memory into the backward problem's user_data */
    let shell = CVodeCreate(CV_BDF, &sunctx);
    let outer = std::mem::replace(cv_mem, *shell);
    nested.cv_user_data = Some(Box::new(outer));

    let flag = CVode(&mut nested, tBout, &mut yB, tBret, itaskB);

    /* restore the forward memory and reattach the backward problem */
    let outer = nested
        .cv_user_data
        .take()
        .unwrap()
        .downcast::<CVodeMem>()
        .unwrap();
    *cv_mem = *outer;

    let ca_mem = adj(cv_mem);
    ca_mem.cvB_mem[idx].cv_y = yB;
    ca_mem.cvB_mem[idx].cv_mem = nested;

    flag
}

pub fn CVodeGetB(cv_mem: &mut CVodeMem, which: i32, tret: &mut f64, yB: &mut NVector) -> i32 {
    let idx = match cva_which_index(cv_mem, which, "CVodeGetB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    let cvB_mem = &adj(cv_mem).cvB_mem[idx];
    N_VScale(ONE, &cvB_mem.cv_y, yB);
    *tret = cvB_mem.cv_tout;

    CV_SUCCESS
}

/*
 * CVodeGetQuadB
 */
pub fn CVodeGetQuadB(cv_mem: &mut CVodeMem, which: i32, tret: &mut f64, qB: &mut NVector) -> i32 {
    let idx = match cva_which_index(cv_mem, which, "CVodeGetQuadB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    let cvB_mem = &adj(cv_mem).cvB_mem[idx];

    /* If the integration for this backward problem has not started yet,
     * simply return the current value of qB (i.e. the final conditions) */
    let nstB = cvB_mem.cv_mem.cv_nst; /* C: CVodeGetNumSteps */

    if nstB == 0 {
        N_VScale(ONE, &cvB_mem.cv_mem.cv_znQ[0], qB);
        *tret = cvB_mem.cv_tout;
        CV_SUCCESS
    } else {
        CVodeGetQuad(&cvB_mem.cv_mem, tret, qB)
    }
}

/*
 * =================================================================
 * PRIVATE FUNCTIONS FOR CHECK POINTS
 * =================================================================
 */

/*
 * CVAckpntInit
 *
 * This routine initializes the check point linked list with
 * information from the initial time.
 */
fn CVAckpntInit(cv_mem: &mut CVodeMem) -> CVckpntMem {
    let ns = cv_mem.cv_Ns as usize;

    let mut ck_mem = CVckpntMem {
        ck_t0: cv_mem.cv_tn,
        ck_t1: ZERO,
        ck_zn: std::iter::repeat_with(NVector::default).take(L_MAX).collect(),
        ck_quadr: SUNFALSE,
        ck_znQ: std::iter::repeat_with(NVector::default).take(L_MAX).collect(),
        ck_sensi: SUNFALSE,
        ck_Ns: 0,
        ck_znS: std::array::from_fn(|_| Vec::new()),
        ck_quadr_sensi: SUNFALSE,
        ck_znQS: std::array::from_fn(|_| Vec::new()),
        /* ck_mem->ck_zn[qmax] was not allocated */
        ck_zqm: 0,
        ck_nst: 0,
        ck_tretlast: ZERO,
        ck_q: 1,
        ck_qprime: 0,
        ck_qwait: 0,
        ck_L: 0,
        ck_gammap: ZERO,
        ck_h: ZERO,
        ck_hprime: ZERO,
        ck_hscale: ZERO,
        ck_eta: ZERO,
        ck_etamax: ZERO,
        ck_tau: [ZERO; L_MAX + 1],
        ck_tq: [ZERO; NUM_TESTS + 1],
        ck_l: [ZERO; L_MAX],
        ck_saved_tq5: ZERO,
    };

    /* Load ckdata from cv_mem (zn[1] is allocated in C but only zn[0] is
       loaded; the Rust CVAckpntGet reinitializes from zn[0] only) */
    ck_mem.ck_zn[0] = cv_mem.cv_zn[0].clone();
    ck_mem.ck_zn[1] = NVector { data: vec![ZERO; cv_mem.cv_tempv.data.len()] };

    /* Do we need to carry quadratures */
    ck_mem.ck_quadr = cv_mem.cv_quadr && cv_mem.cv_errconQ;
    if ck_mem.ck_quadr {
        ck_mem.ck_znQ[0] = cv_mem.cv_znQ[0].clone();
    }

    /* Do we need to carry sensitivities? */
    ck_mem.ck_sensi = cv_mem.cv_sensi;
    if ck_mem.ck_sensi {
        ck_mem.ck_Ns = cv_mem.cv_Ns;
        ck_mem.ck_znS[0] = (0..ns).map(|is| cv_mem.cv_znS[0][is].clone()).collect();
    }

    /* Do we need to carry quadrature sensitivities? */
    ck_mem.ck_quadr_sensi = cv_mem.cv_quadr_sensi && cv_mem.cv_errconQS;
    if ck_mem.ck_quadr_sensi {
        ck_mem.ck_znQS[0] = (0..ns).map(|is| cv_mem.cv_znQS[0][is].clone()).collect();
    }

    ck_mem
}

/*
 * CVAckpntNew
 *
 * This routine allocates space for a new check point and sets
 * its data from current values in cv_mem.
 */
fn CVAckpntNew(cv_mem: &mut CVodeMem) -> CVckpntMem {
    let q = cv_mem.cv_q as usize;
    let qmax = cv_mem.cv_qmax as usize;
    let ns = cv_mem.cv_Ns as usize;

    let mut ck_mem = CVckpntMem {
        ck_t0: cv_mem.cv_tn,
        ck_t1: ZERO,
        ck_zn: std::iter::repeat_with(NVector::default).take(L_MAX).collect(),
        ck_quadr: cv_mem.cv_quadr && cv_mem.cv_errconQ,
        ck_znQ: std::iter::repeat_with(NVector::default).take(L_MAX).collect(),
        ck_sensi: cv_mem.cv_sensi,
        ck_Ns: cv_mem.cv_Ns,
        ck_znS: std::array::from_fn(|_| Vec::new()),
        ck_quadr_sensi: cv_mem.cv_quadr_sensi && cv_mem.cv_errconQS,
        ck_znQS: std::array::from_fn(|_| Vec::new()),
        /* NOTE: zn(qmax) may be needed for a hot restart, if an order
         * increase is deemed necessary at the first step after a check
         * point */
        ck_zqm: if cv_mem.cv_q < cv_mem.cv_qmax { cv_mem.cv_qmax } else { 0 },
        ck_nst: cv_mem.cv_nst,
        ck_tretlast: cv_mem.cv_tretlast,
        ck_q: cv_mem.cv_q,
        ck_qprime: cv_mem.cv_qprime,
        ck_qwait: cv_mem.cv_qwait,
        ck_L: cv_mem.cv_L,
        ck_gammap: cv_mem.cv_gammap,
        ck_h: cv_mem.cv_h,
        ck_hprime: cv_mem.cv_hprime,
        ck_hscale: cv_mem.cv_hscale,
        ck_eta: cv_mem.cv_eta,
        ck_etamax: cv_mem.cv_etamax,
        ck_tau: cv_mem.cv_tau,
        ck_tq: cv_mem.cv_tq,
        ck_l: [ZERO; L_MAX],
        ck_saved_tq5: cv_mem.cv_saved_tq5,
    };

    /* Load check point data from cv_mem */

    for j in 0..=q {
        ck_mem.ck_zn[j] = cv_mem.cv_zn[j].clone();
    }
    if q < qmax {
        ck_mem.ck_zn[qmax] = cv_mem.cv_zn[qmax].clone();
    }

    if ck_mem.ck_quadr {
        for j in 0..=q {
            ck_mem.ck_znQ[j] = cv_mem.cv_znQ[j].clone();
        }
        if q < qmax {
            ck_mem.ck_znQ[qmax] = cv_mem.cv_znQ[qmax].clone();
        }
    }

    if ck_mem.ck_sensi {
        for j in 0..=q {
            ck_mem.ck_znS[j] = (0..ns).map(|is| cv_mem.cv_znS[j][is].clone()).collect();
        }
        if q < qmax {
            ck_mem.ck_znS[qmax] = (0..ns).map(|is| cv_mem.cv_znS[qmax][is].clone()).collect();
        }
    }

    if ck_mem.ck_quadr_sensi {
        for j in 0..=q {
            ck_mem.ck_znQS[j] = (0..ns).map(|is| cv_mem.cv_znQS[j][is].clone()).collect();
        }
        if q < qmax {
            ck_mem.ck_znQS[qmax] = (0..ns).map(|is| cv_mem.cv_znQS[qmax][is].clone()).collect();
        }
    }

    for j in 0..=q {
        ck_mem.ck_l[j] = cv_mem.cv_l[j];
    }

    ck_mem
}

/* (CVAckpntDelete and CVAbckpbDelete of the C are absorbed by RAII:
   Vec::clear / dropping CVadjMem frees the same memory.) */

/*
 * =================================================================
 * PRIVATE FUNCTIONS FOR INTERPOLATION
 * =================================================================
 */

/* Dispatch on ca_IMtype (replaces the C ca_IMstore fn pointer). */
fn cvaIMstore_dispatch(cv_mem: &mut CVodeMem, idx: usize) -> i32 {
    if cv_mem.cv_adj_mem.as_ref().unwrap().ca_IMtype == CV_HERMITE {
        cvaHermiteStorePnt(cv_mem, idx)
    } else {
        cvaPolynomialStorePnt(cv_mem, idx)
    }
}

/* Dispatch on ca_IMtype (replaces the C ca_IMget fn pointer).
   An empty yS Vec plays the role of the C NULL argument. */
fn cvaIMget_dispatch(cv_mem: &mut CVodeMem, t: f64, y: &mut NVector, yS: &mut Vec<NVector>) -> i32 {
    if cv_mem.cv_adj_mem.as_ref().unwrap().ca_IMtype == CV_HERMITE {
        cvaHermiteGetY(cv_mem, t, y, yS)
    } else {
        cvaPolynomialGetY(cv_mem, t, y, yS)
    }
}

/*
 * CVAdataStore
 *
 * This routine integrates the forward model starting at the check
 * point ck_mem and stores y and yprime at all intermediate steps.
 *
 * Return values:
 * CV_SUCCESS
 * CV_REIFWD_FAIL
 * CV_FWD_FAIL
 */
fn CVAdataStore(cv_mem: &mut CVodeMem, ck_idx: usize) -> i32 {
    /* Initialize cv_mem with data from ck_mem */
    let flag = CVAckpntGet(cv_mem, ck_idx);
    if flag != CV_SUCCESS {
        return CV_REIFWD_FAIL;
    }

    /* Set first structure in dt_mem[0] */
    {
        let ca_mem = adj(cv_mem);
        ca_mem.dt_mem[0].t = ca_mem.ck_mem[ck_idx].ck_t0;
    }
    cvaIMstore_dispatch(cv_mem, 0);

    /* Decide whether TSTOP must be activated */
    if adj(cv_mem).ca_tstopCVodeFcall {
        let tstop = adj(cv_mem).ca_tstopCVodeF;
        CVodeSetStopTime(cv_mem, tstop);
    }

    let (tinitial, tfinal, ck_t1) = {
        let ca_mem = cv_mem.cv_adj_mem.as_ref().unwrap();
        (ca_mem.ca_tinitial, ca_mem.ca_tfinal, ca_mem.ck_mem[ck_idx].ck_t1)
    };
    let sign: f64 = if tfinal - tinitial > ZERO { 1.0 } else { -1.0 };

    /* Run CVode to set following structures in dt_mem[i] */
    let mut i: i64 = 1;
    loop {
        let mut ytmp = std::mem::take(&mut adj(cv_mem).ca_ytmp);
        let mut t = ZERO;
        let flag = CVode(cv_mem, ck_t1, &mut ytmp, &mut t, CV_ONE_STEP);
        adj(cv_mem).ca_ytmp = ytmp;
        if flag < 0 {
            return CV_FWD_FAIL;
        }

        adj(cv_mem).dt_mem[i as usize].t = t;
        cvaIMstore_dispatch(cv_mem, i as usize);
        i += 1;

        if sign * (ck_t1 - t) <= ZERO {
            break;
        }
    }

    let ca_mem = adj(cv_mem);
    ca_mem.ca_IMnewData = SUNTRUE; /* New data is now available    */
    ca_mem.ca_ckpntData = Some(ck_idx); /* starting at this check point */
    ca_mem.ca_np = i; /* and we have this many points */

    CV_SUCCESS
}

/*
 * CVAckpntGet
 *
 * This routine prepares CVODES for a hot restart from
 * the check point ck_mem
 */
fn CVAckpntGet(cv_mem: &mut CVodeMem, ck_idx: usize) -> i32 {
    /* Detach the adjoint memory so the check point data can be borrowed
       while cv_mem is reinitialized */
    let ca_mem = cv_mem.cv_adj_mem.take().unwrap();
    let flag = cva_ckpnt_get_inner(cv_mem, &ca_mem, ck_idx);
    cv_mem.cv_adj_mem = Some(ca_mem);
    flag
}

fn cva_ckpnt_get_inner(cv_mem: &mut CVodeMem, ca_mem: &CVadjMem, ck_idx: usize) -> i32 {
    let ck_mem = &ca_mem.ck_mem[ck_idx];

    if ck_idx == 0 {
        /* ck_next == NULL: this is the check point at the initial time.
         * In this case, we just call the reinitialization routine,
         * but make sure we use the same initial stepsize as on
         * the first run. */

        CVodeSetInitStep(cv_mem, cv_mem.cv_h0u);

        let flag = CVodeReInit(cv_mem, ck_mem.ck_t0, &ck_mem.ck_zn[0]);
        if flag != CV_SUCCESS {
            return flag;
        }

        if ck_mem.ck_quadr {
            let flag = CVodeQuadReInit(cv_mem, &ck_mem.ck_znQ[0]);
            if flag != CV_SUCCESS {
                return flag;
            }
        }

        if ck_mem.ck_sensi {
            let flag = CVodeSensReInit(cv_mem, cv_mem.cv_ism, &ck_mem.ck_znS[0]);
            if flag != CV_SUCCESS {
                return flag;
            }
        }

        if ck_mem.ck_quadr_sensi {
            let flag = CVodeQuadSensReInit(cv_mem, &ck_mem.ck_znQS[0]);
            if flag != CV_SUCCESS {
                return flag;
            }
        }
    } else {
        let qmax = cv_mem.cv_qmax as usize;
        let ns = cv_mem.cv_Ns as usize;

        /* Copy parameters from check point data structure */

        cv_mem.cv_nst = ck_mem.ck_nst;
        cv_mem.cv_tretlast = ck_mem.ck_tretlast;
        cv_mem.cv_q = ck_mem.ck_q;
        cv_mem.cv_qprime = ck_mem.ck_qprime;
        cv_mem.cv_qwait = ck_mem.ck_qwait;
        cv_mem.cv_L = ck_mem.ck_L;
        cv_mem.cv_gammap = ck_mem.ck_gammap;
        cv_mem.cv_h = ck_mem.ck_h;
        cv_mem.cv_hprime = ck_mem.ck_hprime;
        cv_mem.cv_hscale = ck_mem.ck_hscale;
        cv_mem.cv_eta = ck_mem.ck_eta;
        cv_mem.cv_etamax = ck_mem.ck_etamax;
        cv_mem.cv_tn = ck_mem.ck_t0;
        cv_mem.cv_saved_tq5 = ck_mem.ck_saved_tq5;

        let q = cv_mem.cv_q as usize;

        /* Copy the arrays from check point data structure */

        for j in 0..=q {
            cv_mem.cv_zn[j].data.copy_from_slice(&ck_mem.ck_zn[j].data);
        }
        if q < qmax {
            cv_mem.cv_zn[qmax].data.copy_from_slice(&ck_mem.ck_zn[qmax].data);
        }

        if ck_mem.ck_quadr {
            for j in 0..=q {
                cv_mem.cv_znQ[j].data.copy_from_slice(&ck_mem.ck_znQ[j].data);
            }
            if q < qmax {
                cv_mem.cv_znQ[qmax].data.copy_from_slice(&ck_mem.ck_znQ[qmax].data);
            }
        }

        if ck_mem.ck_sensi {
            for j in 0..=q {
                for is in 0..ns {
                    cv_mem.cv_znS[j][is].data.copy_from_slice(&ck_mem.ck_znS[j][is].data);
                }
            }
            if q < qmax {
                for is in 0..ns {
                    cv_mem.cv_znS[qmax][is].data.copy_from_slice(&ck_mem.ck_znS[qmax][is].data);
                }
            }
        }

        if ck_mem.ck_quadr_sensi {
            for j in 0..=q {
                for is in 0..ns {
                    cv_mem.cv_znQS[j][is].data.copy_from_slice(&ck_mem.ck_znQS[j][is].data);
                }
            }
            if q < qmax {
                for is in 0..ns {
                    cv_mem.cv_znQS[qmax][is].data.copy_from_slice(&ck_mem.ck_znQS[qmax][is].data);
                }
            }
        }

        cv_mem.cv_tau = ck_mem.ck_tau;
        cv_mem.cv_tq = ck_mem.ck_tq;
        for j in 0..=q {
            cv_mem.cv_l[j] = ck_mem.ck_l[j];
        }

        /* Force a call to setup */
        cv_mem.cv_forceSetup = SUNTRUE;
    }

    CV_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Functions for interpolation
 * -----------------------------------------------------------------
 */

/*
 * CVAfindIndex
 *
 * Finds the index in the array of data point structures such that
 *     dt_mem[index-1].t <= t < dt_mem[index].t
 * If index is changed from the previous invocation, then newpoint = SUNTRUE
 *
 * If t is beyond the leftmost limit, but close enough, index=0.
 *
 * Returns CV_SUCCESS if successful and CV_GETY_BADT if unable to
 * find index (t is too far beyond limits).
 */
fn CVAfindIndex(cv_mem: &mut CVodeMem, t: f64, index: &mut i64, newpoint: &mut bool) -> i32 {
    let uround = cv_mem.cv_uround;
    let ca_mem = adj(cv_mem);

    *newpoint = SUNFALSE;

    /* Find the direction of integration */
    let sign: f64 = if ca_mem.ca_tfinal - ca_mem.ca_tinitial > ZERO { 1.0 } else { -1.0 };

    /* If this is the first time we use new data */
    if ca_mem.ca_IMnewData {
        ca_mem.ca_ilast = ca_mem.ca_np - 1;
        *newpoint = SUNTRUE;
        ca_mem.ca_IMnewData = SUNFALSE;
    }

    /* Search for index starting from ilast */
    let dt_mem = &ca_mem.dt_mem;
    /* (guard ilast == 0: the C indexes dt_mem[-1] here, whose garbage
       value can only route through the equivalent index == 0 path) */
    let to_left = ca_mem.ca_ilast >= 1
        && sign * (t - dt_mem[(ca_mem.ca_ilast - 1) as usize].t) < ZERO;
    let to_right = sign * (t - dt_mem[ca_mem.ca_ilast as usize].t) > ZERO;

    if to_left {
        /* look for a new index to the left */

        *newpoint = SUNTRUE;

        *index = ca_mem.ca_ilast;
        loop {
            if *index == 0 {
                break;
            }
            if sign * (t - dt_mem[(*index - 1) as usize].t) <= ZERO {
                *index -= 1;
            } else {
                break;
            }
        }

        if *index == 0 {
            ca_mem.ca_ilast = 1;
        } else {
            ca_mem.ca_ilast = *index;
        }

        if *index == 0 {
            /* t is beyond leftmost limit. Is it too far? */
            if SUNRabs(t - dt_mem[0].t) > FUZZ_FACTOR * uround {
                return CV_GETY_BADT;
            }
        }
    } else if to_right {
        /* look for a new index to the right */

        *newpoint = SUNTRUE;

        *index = ca_mem.ca_ilast;
        loop {
            if sign * (t - dt_mem[*index as usize].t) > ZERO {
                *index += 1;
            } else {
                break;
            }
        }

        ca_mem.ca_ilast = *index;
    } else {
        /* ilast is still OK */
        *index = ca_mem.ca_ilast;
    }

    CV_SUCCESS
}

/*
 * CVodeGetAdjY
 *
 * This routine returns the interpolated forward solution at time t.
 * The user must allocate space for y.
 */
pub fn CVodeGetAdjY(cv_mem: &mut CVodeMem, t: f64, y: &mut NVector) -> i32 {
    let mut no_yS: Vec<NVector> = Vec::new();
    cvaIMget_dispatch(cv_mem, t, y, &mut no_yS)
}

/*
 * -----------------------------------------------------------------
 * Functions specific to cubic Hermite interpolation
 * -----------------------------------------------------------------
 */

/*
 * cvaHermiteMalloc (C: CVAhermiteMalloc)
 *
 * This routine allocates memory for storing information at all
 * intermediate points between two consecutive check points.
 * This data is then used to interpolate the forward solution
 * at any other time.
 */
fn cvaHermiteMalloc(cv_mem: &mut CVodeMem) {
    let tmpl = cv_mem.cv_tempv.clone();
    let ns = cv_mem.cv_Ns as usize;
    let ca_mem = adj(cv_mem);

    /* Allocate space for the vectors ytmp and yStmp */
    ca_mem.ca_ytmp = tmpl.clone();
    if ca_mem.ca_IMstoreSensi {
        ca_mem.ca_yStmp = (0..ns).map(|_| tmpl.clone()).collect();
    }

    /* Allocate space for the content field of the dt structures */
    let store_sensi = ca_mem.ca_IMstoreSensi;
    for i in 0..=(ca_mem.ca_nsteps as usize) {
        ca_mem.dt_mem[i].content = DtpntContent::Hermite {
            y: tmpl.clone(),
            yd: tmpl.clone(),
            yS: if store_sensi { (0..ns).map(|_| tmpl.clone()).collect() } else { Vec::new() },
            ySd: if store_sensi { (0..ns).map(|_| tmpl.clone()).collect() } else { Vec::new() },
        };
    }
}

/* (cvaHermiteFree is absorbed by RAII.) */

/*
 * cvaHermiteStorePnt (C: CVAhermiteStorePnt -> IMstore)
 *
 * This routine stores a new point (y,yd) in the structure d for use
 * in the cubic Hermite interpolation.
 * Note that the time is already stored.
 */
fn cvaHermiteStorePnt(cv_mem: &mut CVodeMem, idx: usize) -> i32 {
    /* detach the content so cv_mem can be borrowed while filling it */
    let mut content = std::mem::replace(&mut adj(cv_mem).dt_mem[idx].content, DtpntContent::None);
    let store_sensi = adj(cv_mem).ca_IMstoreSensi;
    let ns = cv_mem.cv_Ns as usize;

    if let DtpntContent::Hermite { y, yd, yS, ySd } = &mut content {
        /* Load solution */
        y.data.copy_from_slice(&cv_mem.cv_zn[0].data);

        if store_sensi {
            for is in 0..ns {
                yS[is].data.copy_from_slice(&cv_mem.cv_znS[0][is].data);
            }
        }

        /* Load derivative */
        if cv_mem.cv_nst == 0 {
            let f = cv_mem.cv_f.unwrap();
            /* retval = */
            let _ = f(cv_mem.cv_tn, y, yd, &mut cv_mem.cv_user_data);

            if store_sensi {
                let tn = cv_mem.cv_tn;
                let mut wrk1 = std::mem::take(&mut cv_mem.cv_tempv);
                let mut wrk2 = std::mem::take(&mut cv_mem.cv_ftemp);
                /* retval = */
                let _ = cvSensRhsWrapper(cv_mem, tn, y, yd, yS, ySd, &mut wrk1, &mut wrk2);
                cv_mem.cv_tempv = wrk1;
                cv_mem.cv_ftemp = wrk2;
            }
        } else {
            let hinv = ONE / cv_mem.cv_h;
            N_VScale(hinv, &cv_mem.cv_zn[1], yd);

            if store_sensi {
                for is in 0..ns {
                    N_VScale(hinv, &cv_mem.cv_znS[1][is], &mut ySd[is]);
                }
            }
        }
    }

    adj(cv_mem).dt_mem[idx].content = content;

    0
}

/*
 * cvaHermiteGetY (C: CVAhermiteGetY -> IMget)
 *
 * This routine uses cubic piece-wise Hermite interpolation for
 * the forward solution vector.
 * It is typically called by the wrapper routines before calling
 * user provided routines (fB, djacB, bjacB, jtimesB, psolB) but
 * can be directly called by the user through CVodeGetAdjY
 *
 * (An empty yS Vec plays the role of the C NULL argument.)
 */
pub(crate) fn cvaHermiteGetY(
    cv_mem: &mut CVodeMem,
    t: f64,
    y: &mut NVector,
    yS: &mut Vec<NVector>,
) -> i32 {
    /* Get the index in dt_mem */
    let mut index: i64 = 0;
    let mut newpoint = SUNFALSE;
    let flag = CVAfindIndex(cv_mem, t, &mut index, &mut newpoint);
    if flag != CV_SUCCESS {
        return flag;
    }

    let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();

    /* Local value of Ns */
    let ns = if ca_mem.ca_IMinterpSensi && !yS.is_empty() { yS.len() } else { 0 };

    let CVadjMem { dt_mem, ca_Y, ca_YS, .. } = &mut **ca_mem;

    /* If we are beyond the left limit but close enough,
       then return y at the left limit. */
    if index == 0 {
        if let DtpntContent::Hermite { y: y0, yS: yS0, .. } = &dt_mem[0].content {
            y.data.copy_from_slice(&y0.data);
            for is in 0..ns {
                yS[is].data.copy_from_slice(&yS0[is].data);
            }
        }
        return CV_SUCCESS;
    }

    /* Extract stuff from the appropriate data points */
    let idx = index as usize;
    let t0 = dt_mem[idx - 1].t;
    let t1 = dt_mem[idx].t;
    let delta = t1 - t0;

    let (y0, yd0, ys0, ysd0) = match &dt_mem[idx - 1].content {
        DtpntContent::Hermite { y, yd, yS, ySd } => (y, yd, yS, ySd),
        _ => unreachable!(),
    };

    if newpoint {
        /* Recompute Y0 and Y1 */
        let (y1, yd1, ys1, ysd1) = match &dt_mem[idx].content {
            DtpntContent::Hermite { y, yd, yS, ySd } => (y, yd, yS, ySd),
            _ => unreachable!(),
        };

        /* Y1 = delta (yd1 + yd0) - 2 (y1 - y0)
           (N_VLinearCombination(4, {-2, 2, delta, delta},
            {y1, y0, yd1, yd0}, Y[1])) */
        {
            let (yfront, yback) = ca_Y.split_at_mut(1);
            let y_1 = &mut yback[0];
            let _ = &yfront; /* Y[0] untouched here */
            for k in 0..y_1.data.len() {
                y_1.data[k] =
                    -TWO * y1.data[k] + TWO * y0.data[k] + delta * yd1.data[k] + delta * yd0.data[k];
            }
        }

        /* Y0 = y1 - y0 - delta * yd0
           (N_VLinearCombination(3, {1, -1, -delta}, {y1, y0, yd0}, Y[0])) */
        {
            let y_0 = &mut ca_Y[0];
            for k in 0..y_0.data.len() {
                y_0.data[k] = ONE * y1.data[k] + (-ONE) * y0.data[k] + (-delta) * yd0.data[k];
            }
        }

        /* Recompute YS0 and YS1, if needed */
        if ns > 0 {
            /* YS1 = delta (ySd1 + ySd0) - 2 (yS1 - yS0) */
            {
                let (ysfront, ysback) = ca_YS.split_at_mut(1);
                let ys_1 = &mut ysback[0];
                let _ = &ysfront;
                for is in 0..ns {
                    for k in 0..ys_1[is].data.len() {
                        ys_1[is].data[k] = -TWO * ys1[is].data[k]
                            + TWO * ys0[is].data[k]
                            + delta * ysd1[is].data[k]
                            + delta * ysd0[is].data[k];
                    }
                }
            }

            /* YS0 = yS1 - yS0 - delta * ySd0 */
            {
                let ys_0 = &mut ca_YS[0];
                for is in 0..ns {
                    for k in 0..ys_0[is].data.len() {
                        ys_0[is].data[k] = ONE * ys1[is].data[k]
                            + (-ONE) * ys0[is].data[k]
                            + (-delta) * ysd0[is].data[k];
                    }
                }
            }
        }
    }

    /* Perform the actual interpolation. */

    let factor1 = t - t0;

    let mut factor2 = factor1 / delta;
    factor2 = factor2 * factor2;

    let factor3 = factor2 * (t - t1) / delta;

    /* y = y0 + factor1 yd0 + factor2 * Y[0] + factor3 Y[1]
       (N_VLinearCombination(4, {1, factor1, factor2, factor3},
        {y0, yd0, Y[0], Y[1]}, y)) */
    for k in 0..y.data.len() {
        y.data[k] = ONE * y0.data[k]
            + factor1 * yd0.data[k]
            + factor2 * ca_Y[0].data[k]
            + factor3 * ca_Y[1].data[k];
    }

    /* yS = yS0 + factor1 ySd0 + factor2 * YS[0] + factor3 YS[1], if needed */
    for is in 0..ns {
        for k in 0..yS[is].data.len() {
            yS[is].data[k] = ONE * ys0[is].data[k]
                + factor1 * ysd0[is].data[k]
                + factor2 * ca_YS[0][is].data[k]
                + factor3 * ca_YS[1][is].data[k];
        }
    }

    CV_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Functions specific to Polynomial interpolation
 * -----------------------------------------------------------------
 */

/*
 * cvaPolynomialMalloc (C: CVApolynomialMalloc)
 *
 * This routine allocates memory for storing information at all
 * intermediate points between two consecutive check points.
 * This data is then used to interpolate the forward solution
 * at any other time.
 */
fn cvaPolynomialMalloc(cv_mem: &mut CVodeMem) {
    let tmpl = cv_mem.cv_tempv.clone();
    let ns = cv_mem.cv_Ns as usize;
    let ca_mem = adj(cv_mem);

    /* Allocate space for the vectors ytmp and yStmp */
    ca_mem.ca_ytmp = tmpl.clone();
    if ca_mem.ca_IMstoreSensi {
        ca_mem.ca_yStmp = (0..ns).map(|_| tmpl.clone()).collect();
    }

    /* Allocate space for the content field of the dt structures */
    let store_sensi = ca_mem.ca_IMstoreSensi;
    for i in 0..=(ca_mem.ca_nsteps as usize) {
        ca_mem.dt_mem[i].content = DtpntContent::Polynomial {
            y: tmpl.clone(),
            yS: if store_sensi { (0..ns).map(|_| tmpl.clone()).collect() } else { Vec::new() },
            order: 0,
        };
    }
}

/* (cvaPolynomialFree is absorbed by RAII.) */

/*
 * cvaPolynomialStorePnt (C: CVApolynomialStorePnt -> IMstore)
 *
 * This routine stores a new point y in the structure d for use
 * in the Polynomial interpolation.
 * Note that the time is already stored.
 */
fn cvaPolynomialStorePnt(cv_mem: &mut CVodeMem, idx: usize) -> i32 {
    let mut content = std::mem::replace(&mut adj(cv_mem).dt_mem[idx].content, DtpntContent::None);
    let store_sensi = adj(cv_mem).ca_IMstoreSensi;
    let ns = cv_mem.cv_Ns as usize;

    if let DtpntContent::Polynomial { y, yS, order } = &mut content {
        y.data.copy_from_slice(&cv_mem.cv_zn[0].data);

        if store_sensi {
            for is in 0..ns {
                yS[is].data.copy_from_slice(&cv_mem.cv_znS[0][is].data);
            }
        }

        *order = cv_mem.cv_qu;
    }

    adj(cv_mem).dt_mem[idx].content = content;

    0
}

/*
 * cvaPolynomialGetY (C: CVApolynomialGetY -> IMget)
 *
 * This routine uses polynomial interpolation for the forward solution
 * vector. It is typically called by the wrapper routines before calling
 * user provided routines (fB, djacB, bjacB, jtimesB, psolB) but
 * can be directly called by the user through CVodeGetAdjY.
 *
 * (An empty yS Vec plays the role of the C NULL argument.)
 */
pub(crate) fn cvaPolynomialGetY(
    cv_mem: &mut CVodeMem,
    t: f64,
    y: &mut NVector,
    yS: &mut Vec<NVector>,
) -> i32 {
    /* Get the index in dt_mem */
    let mut index: i64 = 0;
    let mut newpoint = SUNFALSE;
    let flag = CVAfindIndex(cv_mem, t, &mut index, &mut newpoint);
    if flag != CV_SUCCESS {
        return flag;
    }

    let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();

    /* Local value of Ns */
    let ns = if ca_mem.ca_IMinterpSensi && !yS.is_empty() { yS.len() } else { 0 };

    /* If we are beyond the left limit but close enough,
       then return y at the left limit. */
    if index == 0 {
        if let DtpntContent::Polynomial { y: y0, yS: yS0, .. } = &ca_mem.dt_mem[0].content {
            y.data.copy_from_slice(&y0.data);
            for is in 0..ns {
                yS[is].data.copy_from_slice(&yS0[is].data);
            }
        }
        return CV_SUCCESS;
    }

    let idx = index as usize;

    /* Scaling factor */
    let dt = SUNRabs(ca_mem.dt_mem[idx].t - ca_mem.dt_mem[idx - 1].t);

    /* Find the direction of the forward integration */
    let dir: i32 = if ca_mem.ca_tfinal - ca_mem.ca_tinitial > ZERO { 1 } else { -1 };

    /* Establish the base point depending on the integration direction.
       Modify the base if there are not enough points for the current order */
    let (base, order): (i64, i32);
    if dir == 1 {
        let mut b = index;
        let ord = match &ca_mem.dt_mem[b as usize].content {
            DtpntContent::Polynomial { order, .. } => *order,
            _ => unreachable!(),
        };
        if index < ord as i64 {
            b += ord as i64 - index;
        }
        base = b;
        order = ord;
    } else {
        let mut b = index - 1;
        let ord = match &ca_mem.dt_mem[b as usize].content {
            DtpntContent::Polynomial { order, .. } => *order,
            _ => unreachable!(),
        };
        if ca_mem.ca_np - index > ord as i64 {
            b -= index + ord as i64 - ca_mem.ca_np;
        }
        base = b;
        order = ord;
    }

    let ord = order as usize;

    /* Recompute Y (divided differences for Newton polynomial) if needed */
    if newpoint {
        let CVadjMem { dt_mem, ca_Y, ca_YS, ca_T, .. } = &mut **ca_mem;

        /* Store 0-th order DD */
        if dir == 1 {
            for j in 0..=ord {
                let src = (base - j as i64) as usize;
                ca_T[j] = dt_mem[src].t;
                if let DtpntContent::Polynomial { y: cy, yS: cyS, .. } = &dt_mem[src].content {
                    ca_Y[j].data.copy_from_slice(&cy.data);
                    for is in 0..ns {
                        ca_YS[j][is].data.copy_from_slice(&cyS[is].data);
                    }
                }
            }
        } else {
            for j in 0..=ord {
                let src = (base - 1 + j as i64) as usize;
                ca_T[j] = dt_mem[src].t;
                if let DtpntContent::Polynomial { y: cy, yS: cyS, .. } = &dt_mem[src].content {
                    ca_Y[j].data.copy_from_slice(&cy.data);
                    for is in 0..ns {
                        ca_YS[j][is].data.copy_from_slice(&cyS[is].data);
                    }
                }
            }
        }

        /* Compute higher-order DD */
        for i in 1..=ord {
            let mut j = ord;
            while j >= i {
                let factor = dt / (ca_T[j] - ca_T[j - i]);
                /* Y[j] = factor*Y[j] - factor*Y[j-1] (aliased) */
                {
                    let (front, back) = ca_Y.split_at_mut(j);
                    back[0].linear_sum_with(factor, -factor, &front[j - 1]);
                }
                if ns > 0 {
                    let (front, back) = ca_YS.split_at_mut(j);
                    for is in 0..ns {
                        back[0][is].linear_sum_with(factor, -factor, &front[j - 1][is]);
                    }
                }
                j -= 1;
            }
        }
    }

    /* Perform the actual interpolation using nested multiplications */
    let mut cvals = [ZERO; L_MAX + 1];
    cvals[0] = ONE;
    for i in 0..ord {
        cvals[i + 1] = cvals[i] * (t - ca_mem.ca_T[i]) / dt;
    }

    /* y = sum_j cvals[j] * Y[j] (N_VLinearCombination) */
    for k in 0..y.data.len() {
        let mut acc = cvals[0] * ca_mem.ca_Y[0].data[k];
        for j in 1..=ord {
            acc += cvals[j] * ca_mem.ca_Y[j].data[k];
        }
        y.data[k] = acc;
    }

    for is in 0..ns {
        for k in 0..yS[is].data.len() {
            let mut acc = cvals[0] * ca_mem.ca_YS[0][is].data[k];
            for j in 1..=ord {
                acc += cvals[j] * ca_mem.ca_YS[j][is].data[k];
            }
            yS[is].data[k] = acc;
        }
    }

    CV_SUCCESS
}

/*
 * =================================================================
 * WRAPPERS FOR ADJOINT SYSTEM
 * =================================================================
 */

/*
 * CVArhs
 *
 * This routine interfaces to the CVRhsFnB (or CVRhsFnBS) routine
 * provided by the user.
 *
 * (Registered as the CVRhsFn of the backward problem; the user_data
 * of the backward problem holds the FORWARD CVodeMem, installed by
 * CVodeB for the duration of the backward integration.)
 */
fn CVArhs(t: f64, yB: &NVector, yBdot: &mut NVector, cvode_mem: &mut UserData) -> i32 {
    let cv_mem = cvode_mem
        .as_mut()
        .unwrap()
        .downcast_mut::<CVodeMem>()
        .unwrap();

    /* Get forward solution from interpolation */
    let interp_sensi = adj(cv_mem).ca_IMinterpSensi;
    let mut ytmp = std::mem::take(&mut adj(cv_mem).ca_ytmp);
    let mut yStmp = if interp_sensi {
        std::mem::take(&mut adj(cv_mem).ca_yStmp)
    } else {
        Vec::new()
    };

    let flag = cvaIMget_dispatch(cv_mem, t, &mut ytmp, &mut yStmp);

    let ca_mem = adj(cv_mem);
    ca_mem.ca_ytmp = ytmp;
    if interp_sensi {
        ca_mem.ca_yStmp = yStmp;
    }

    if flag != CV_SUCCESS {
        cvProcessError(Some(cv_mem), -1, line!(), "CVArhs", file!(),
            &format!("Bad t = {} for interpolation.", t));
        return -1;
    }

    /* Call the user's RHS function */
    let ca_mem = adj(cv_mem);
    let CVadjMem { ca_ytmp, ca_yStmp, cvB_mem, ca_bckpbCrt, .. } = ca_mem;
    let cvB = &mut cvB_mem[ca_bckpbCrt.unwrap()];

    if cvB.cv_f_withSensi {
        (cvB.cv_fs.unwrap())(t, &*ca_ytmp, &*ca_yStmp, yB, yBdot, &mut cvB.cv_user_data)
    } else {
        (cvB.cv_f.unwrap())(t, &*ca_ytmp, yB, yBdot, &mut cvB.cv_user_data)
    }
}

/*
 * CVArhsQ
 *
 * This routine interfaces to the CVQuadRhsFnB (or CVQuadRhsFnBS) routine
 * provided by the user.
 */
fn CVArhsQ(t: f64, yB: &NVector, qBdot: &mut NVector, cvode_mem: &mut UserData) -> i32 {
    let cv_mem = cvode_mem
        .as_mut()
        .unwrap()
        .downcast_mut::<CVodeMem>()
        .unwrap();

    /* Get forward solution from interpolation */
    let interp_sensi = adj(cv_mem).ca_IMinterpSensi;
    let mut ytmp = std::mem::take(&mut adj(cv_mem).ca_ytmp);
    let mut yStmp = if interp_sensi {
        std::mem::take(&mut adj(cv_mem).ca_yStmp)
    } else {
        Vec::new()
    };

    /* flag = */
    let _ = cvaIMget_dispatch(cv_mem, t, &mut ytmp, &mut yStmp);

    let ca_mem = adj(cv_mem);
    ca_mem.ca_ytmp = ytmp;
    if interp_sensi {
        ca_mem.ca_yStmp = yStmp;
    }

    /* Call the user's RHS function */
    let ca_mem = adj(cv_mem);
    let CVadjMem { ca_ytmp, ca_yStmp, cvB_mem, ca_bckpbCrt, .. } = ca_mem;
    let cvB = &mut cvB_mem[ca_bckpbCrt.unwrap()];

    if cvB.cv_fQ_withSensi {
        (cvB.cv_fQs.unwrap())(t, &*ca_ytmp, &*ca_yStmp, yB, qBdot, &mut cvB.cv_user_data)
    } else {
        (cvB.cv_fQ.unwrap())(t, &*ca_ytmp, yB, qBdot, &mut cvB.cv_user_data)
    }
}
