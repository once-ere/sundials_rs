/* -----------------------------------------------------------------
 * Translated from src/idas/idaa.c (IDAS 7.7.0).
 *
 * Adjoint sensitivity analysis: checkpointing (IDASolveF), backward
 * problems (IDACreateB/IDAInitB/IDASolveB, ...), and the cubic
 * Hermite / variable-order polynomial interpolation modules.
 *
 * Port notes (pinned decisions, mirroring cvodea.rs / cvodes_impl.rs;
 * see idas_impl.rs):
 *  - The C check point linked list (head = MOST RECENT check point,
 *    ck_next walks back in time, tail = initial check point) becomes
 *    Vec<IDAckpntMem> with index 0 = initial check point and
 *    last() = most recent; "ck_next == NULL" <=> index == 0, and
 *    "ck_mem = ck_mem->ck_next" <=> index -= 1.
 *  - The C backward-problem list (C prepends: head = newest) becomes
 *    Vec<IDABMem> in creation order, so ida_index == position.
 *    Loops "over all backward problems" therefore run in creation
 *    order instead of the C's reverse order; the problems are
 *    mutually independent so the results are unchanged.
 *  - The C interpolation-module function pointers ia_malloc/ia_free/
 *    ia_storePnt/ia_getY are replaced by dispatch on ia_interpType
 *    (idaa_storePnt_dispatch / IDAAgetY).  IDAAgetY is also the body
 *    of the idas_ls.rs idaLsGetY bridge (one dispatch point serves
 *    the LS and BBDPre backward wrappers, per the recorded PIN).
 *    The free routines are absorbed by RAII (IDAAdjFree drops
 *    IDAadjMem).
 *  - In C, ia_Y[i]/ia_YS[i] are aliases of phi[i]/phiS[i] used as
 *    interpolation scratch; here they are owned workspace vectors
 *    allocated in IDASolveF (every use overwrites them before
 *    reading, so copying semantics are identical).
 *  - In C, each backward problem's IDAS memory has user_data == the
 *    FORWARD ida_mem (set once in IDACreateB).  Here that
 *    self-reference is created transiently: IDASolveB / IDACalcICB(S)
 *    move the outer IDAMem into the nested problem's ida_user_data
 *    around each IDASolve/IDACalcIC call on the backward problem
 *    (and IDAAres/IDAArhsQ and the idas_ls.rs / idas_bbdpre.rs
 *    B-wrappers downcast it back out).
 *  - An empty Vec<NVector> plays the role of the C NULL N_Vector*
 *    argument of the getY routines.
 *  - The C SUNDIALS_MARK_FUNCTION profiler brackets vanish (the
 *    workspace profiler is a no-op shell).
 * -----------------------------------------------------------------*/

use crate::idas::{
    ida_msg_g, IDACreate, IDAGetQuad, IDAGetSolution, IDAInit, IDAQuadInit, IDAQuadReInit,
    IDAQuadSStolerances, IDAQuadSVtolerances, IDAQuadSensReInit, IDAReInit, IDASStolerances,
    IDASVtolerances, IDASensReInit, IDASolve,
};
use crate::idas_ic::IDACalcIC;
use crate::idas_impl::*;
use crate::idas_io::{IDASetInitStep, IDASetStopTime};
use crate::nvector_serial::*;
use crate::sundials_math::*;
use crate::sundials_types::*;

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;
const HUNDRED: f64 = 100.0;

/* fuzz factor for IDAAgetY */
const FUZZ_FACTOR: f64 = 1000000.0;

/* Helper: mutable access to the (attached) adjoint memory. */
fn adj(ida_mem: &mut IDAMem) -> &mut IDAadjMem {
    ida_mem.ida_adj_mem.as_mut().unwrap()
}

/*=================================================================*/
/*                  Exported Functions                             */
/*=================================================================*/

/*
 * IDAAdjInit
 *
 * This routine allocates space for the global IDAA memory
 * structure.
 */
pub fn IDAAdjInit(ida_mem: &mut IDAMem, steps: i64, interp: i32) -> i32 {
    /* Check arguments */

    if steps <= 0 {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAAdjInit", file!(),
                        MSGAM_BAD_STEPS);
        return IDA_ILL_INPUT;
    }

    if interp != IDA_HERMITE && interp != IDA_POLYNOMIAL {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAAdjInit", file!(),
                        MSGAM_BAD_INTERP);
        return IDA_ILL_INPUT;
    }

    /* Allocate memory block for IDAadjMem (infallible; the C
       malloc-failure branch vanishes).  IDAAdataMalloc: allocate
       space for the array of Data Point structures. */
    let mut dt_mem = Vec::with_capacity((steps + 1) as usize);
    for _ in 0..=steps {
        dt_mem.push(IDAdtpntMem { t: ZERO, content: DtpntContent::None });
    }

    let idaadj_mem = IDAadjMem {
        /* Forward-problem data (set by IDASolveF) */
        ia_tinitial: ZERO,
        ia_tfinal: ZERO,

        /* Initialization of check points. */
        ck_mem: Vec::new(),
        ia_nckpnts: 0,
        ia_ckpntData: None,

        /* Initialization of interpolation data. */
        ia_interpType: interp,
        ia_nsteps: steps,

        /* Last index used in IDAAfindIndex, initialize to invalid value */
        ia_ilast: -1,

        dt_mem,
        ia_np: 0,

        /* (The interpolation-module function pointers of the C become
           dispatch on ia_interpType; see idaa_storePnt_dispatch /
           IDAAgetY.) */

        /* The interpolation module has not been initialized yet */
        ia_mallocDone: SUNFALSE,
        ia_newData: SUNFALSE,

        /* By default we will store but not interpolate sensitivities
         *  - storeSensi will be set in IDASolveF to SUNFALSE if FSA is not
         *    enabled or if the user forced this through IDAAdjSetNoSensi
         *  - interpSensi will be set in IDASolveB to SUNTRUE if storeSensi
         *    is SUNTRUE and if at least one backward problem requires
         *    sensitivities
         *  - noInterp will be set in IDACalcICB to SUNTRUE before the call
         *    to IDACalcIC and SUNFALSE after. */
        ia_storeSensi: SUNTRUE,
        ia_interpSensi: SUNFALSE,
        ia_noInterp: SUNFALSE,

        /* Workspaces (allocated in IDASolveF / the IM malloc routines) */
        ia_Y: Vec::new(),
        ia_YS: std::array::from_fn(|_| Vec::new()),
        ia_T: [ZERO; MXORDP1],
        ia_yyTmp: NVector::default(),
        ia_ypTmp: NVector::default(),
        ia_yySTmp: Vec::new(),
        ia_ypSTmp: Vec::new(),

        /* Initialize backward problems. */
        IDAB_mem: Vec::new(),
        ia_bckpbCrt: None,
        ia_nbckpbs: 0,

        /* IDASolveF and IDASolveB not called yet. */
        ia_firstIDAFcall: SUNTRUE,
        ia_tstopIDAFcall: SUNFALSE,
        ia_tstopIDAF: ZERO,

        ia_firstIDABcall: SUNTRUE,

        ia_rootret: SUNFALSE,
        ia_troot: ZERO,
    };

    /* Attach IDAS memory for forward runs */
    ida_mem.ida_adj_mem = Some(Box::new(idaadj_mem));

    /* Adjoint module initialized and allocated. */
    ida_mem.ida_adj = SUNTRUE;
    ida_mem.ida_adjMallocDone = SUNTRUE;

    IDA_SUCCESS
}

/*
 * IDAAdjReInit
 *
 * IDAAdjReInit reinitializes the IDAS memory structure for ASA
 */
pub fn IDAAdjReInit(ida_mem: &mut IDAMem) -> i32 {
    /* Was ASA previously initialized? */
    if !ida_mem.ida_adjMallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_ADJ, line!(), "IDAAdjReInit", file!(), MSGAM_NO_ADJ);
        return IDA_NO_ADJ;
    }

    let idaadj_mem = adj(ida_mem);

    /* Free all stored checkpoints (RAII drop). */
    idaadj_mem.ck_mem.clear();
    idaadj_mem.ia_nckpnts = 0;
    idaadj_mem.ia_ckpntData = None;

    /* Flags for tracking the first calls to IDASolveF and IDASolveB. */
    idaadj_mem.ia_firstIDAFcall = SUNTRUE;
    idaadj_mem.ia_tstopIDAFcall = SUNFALSE;
    idaadj_mem.ia_firstIDABcall = SUNTRUE;

    IDA_SUCCESS
}

/*
 * IDAAdjFree
 *
 * IDAAdjFree routine frees the memory allocated by IDAAdjInit.
 * (Check points, interpolation data and all backward problems are
 * released by RAII when the IDAadjMem Box is dropped.)
 */
pub fn IDAAdjFree(ida_mem: &mut IDAMem) {
    if ida_mem.ida_adjMallocDone {
        ida_mem.ida_adj_mem = None;
        ida_mem.ida_adjMallocDone = SUNFALSE;
    }
}

/* (IDAAbckpbDelete: the per-problem teardown — IDAFree on the nested
   solver, lfree/pfree hooks, workspace vectors — is RAII here: dropping
   an IDABMem drops its nested Box<IDAMem>, its dyn-Any lmem/pmem
   attachments and its yy/yp workspaces.) */

/*=================================================================*/
/*                    Wrappers for IDAA                            */
/*=================================================================*/

/*
 *                      IDASolveF
 *
 * This routine integrates to tout and returns solution into yout.
 * In the same time, it stores check point data every 'steps' steps.
 *
 * IDASolveF can be called repeatedly by the user. The last tout
 *  will be used as the starting time for the backward integration.
 *
 *  ncheckPtr points to the number of check points stored so far.
 */
pub fn IDASolveF(
    ida_mem: &mut IDAMem,
    tout: f64,
    tret: &mut f64,
    yret: &mut NVector,
    ypret: &mut NVector,
    itask: i32,
    ncheckPtr: &mut i32,
) -> i32 {
    /* Is ASA initialized ? */
    if !ida_mem.ida_adjMallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_ADJ, line!(), "IDASolveF", file!(), MSGAM_NO_ADJ);
        return IDA_NO_ADJ;
    }

    /* Check for yret != NULL */
    if yret.is_empty() {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASolveF", file!(),
                        MSG_YRET_NULL);
        return IDA_ILL_INPUT;
    }
    /* Check for ypret != NULL */
    if ypret.is_empty() {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASolveF", file!(),
                        MSG_YPRET_NULL);
        return IDA_ILL_INPUT;
    }
    /* (tret != NULL check vanishes: &mut receiver) */

    /* Check for valid itask */
    if itask != IDA_NORMAL && itask != IDA_ONE_STEP {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASolveF", file!(),
                        MSG_BAD_ITASK);
        return IDA_ILL_INPUT;
    }

    /* All memory checks done, proceed ... */

    /* If tstop is enabled, store some info */
    if ida_mem.ida_tstopset {
        let tstop = ida_mem.ida_tstop;
        let idaadj_mem = adj(ida_mem);
        idaadj_mem.ia_tstopIDAFcall = SUNTRUE;
        idaadj_mem.ia_tstopIDAF = tstop;
    }

    /* On the first step:
     *   - set tinitial
     *   - initialize list of check points
     *   - if needed, initialize the interpolation module
     *   - load dt_mem[0]
     * On subsequent steps, test if taking a new step is necessary.
     */
    if adj(ida_mem).ia_firstIDAFcall {
        let tn = ida_mem.ida_tn;
        adj(ida_mem).ia_tinitial = tn;

        /* initialize the check point list (allocation is infallible;
           the C NULL/IDA_MEM_FAIL branch vanishes) */
        let ck = IDAAckpntInit(ida_mem);
        adj(ida_mem).ck_mem.push(ck);

        if !adj(ida_mem).ia_mallocDone {
            /* Do we need to store sensitivities? */
            if !ida_mem.ida_sensi {
                adj(ida_mem).ia_storeSensi = SUNFALSE;
            }

            /* Allocate space for interpolation data (infallible) */
            match adj(ida_mem).ia_interpType {
                IDA_HERMITE => IDAAhermiteMalloc(ida_mem),
                _ => IDAApolynomialMalloc(ida_mem),
            }

            /* Rename phi and, if needed, phiS for use in interpolation:
               in C ia_Y[i]/ia_YS[i] alias phi[i]/phiS[i]; here they are
               owned scratch rows (pinned decision, see file header) */
            let ns = ida_mem.ida_Ns as usize;
            let store_sensi = adj(ida_mem).ia_storeSensi;
            let tmpl = std::mem::take(&mut ida_mem.ida_tempv1);
            {
                let idaadj_mem = adj(ida_mem);
                idaadj_mem.ia_Y = (0..MXORDP1).map(|_| N_VClone(&tmpl)).collect();
                if store_sensi {
                    for row in idaadj_mem.ia_YS.iter_mut() {
                        *row = (0..ns).map(|_| N_VClone(&tmpl)).collect();
                    }
                }
            }
            ida_mem.ida_tempv1 = tmpl;

            adj(ida_mem).ia_mallocDone = SUNTRUE;
        }

        let t0 = adj(ida_mem).ck_mem.last().unwrap().ck_t0;
        adj(ida_mem).dt_mem[0].t = t0;
        idaa_storePnt_dispatch(ida_mem, 0);

        adj(ida_mem).ia_firstIDAFcall = SUNFALSE;
    } else if itask == IDA_NORMAL {
        /* When in normal mode, check if tout was passed or if a previous
           root was not reported and return an interpolated solution. No
           changes to ck_mem or dt_mem are needed. */

        /* flag to signal if an early return is needed (holds the return
           flag when set — the C uses a separate earlyret boolean) */
        let mut earlyret: Option<i32> = None;

        /* if a root needs to be reported compare tout to troot otherwise
           compare to the current time tn */
        let ttest = if adj(ida_mem).ia_rootret {
            adj(ida_mem).ia_troot
        } else {
            ida_mem.ida_tn
        };

        if (ttest - tout) * ida_mem.ida_hh >= ZERO {
            /* ttest is after tout, interpolate to tout */
            *tret = tout;
            earlyret = Some(IDAGetSolution(ida_mem, tout, yret, ypret));
        } else if adj(ida_mem).ia_rootret {
            /* tout is after troot, interpolate to troot
               (C discards the IDAGetSolution flag and reports
               IDA_ROOT_RETURN) */
            let troot = adj(ida_mem).ia_troot;
            *tret = troot;
            let _ = IDAGetSolution(ida_mem, troot, yret, ypret);
            adj(ida_mem).ia_rootret = SUNFALSE;
            earlyret = Some(IDA_ROOT_RETURN);
        }

        /* return if necessary */
        if let Some(flag) = earlyret {
            let nst = ida_mem.ida_nst;
            let idaadj_mem = adj(ida_mem);
            *ncheckPtr = idaadj_mem.ia_nckpnts;
            idaadj_mem.ia_newData = SUNTRUE;
            idaadj_mem.ia_ckpntData = Some(idaadj_mem.ck_mem.len() - 1);
            idaadj_mem.ia_np = nst % idaadj_mem.ia_nsteps + 1;
            return flag;
        }
    }

    /* Integrate to tout (in IDA_ONE_STEP mode) while loading check points */
    let mut nstloc: i64 = 0;
    let mut flag: i32;
    loop {
        /* Check for too many steps */

        if ida_mem.ida_mxstep > 0 && nstloc >= ida_mem.ida_mxstep {
            let tn = ida_mem.ida_tn;
            IDAProcessError(Some(ida_mem), IDA_TOO_MUCH_WORK, line!(), "IDASolveF", file!(),
                            &ida_msg_g(MSG_MAX_STEPS, &[tn]));
            flag = IDA_TOO_MUCH_WORK;
            break;
        }

        /* Perform one step of the integration */

        flag = IDASolve(ida_mem, tout, tret, yret, ypret, IDA_ONE_STEP);
        if flag < 0 {
            break;
        }

        nstloc += 1;

        /* Test if a new check point is needed */

        if ida_mem.ida_nst % adj(ida_mem).ia_nsteps == 0 {
            let tn = ida_mem.ida_tn;
            adj(ida_mem).ck_mem.last_mut().unwrap().ck_t1 = tn;

            /* Create a new check point, load it, and append it to the list
               (allocation is infallible; the C NULL/IDA_MEM_FAIL branch
               vanishes; C prepends — the Vec push keeps most-recent last) */
            let tmp = IDAAckpntNew(ida_mem);
            adj(ida_mem).ck_mem.push(tmp);
            adj(ida_mem).ia_nckpnts += 1;

            ida_mem.ida_forceSetup = SUNTRUE;

            /* Reset i=0 and load dt_mem[0] */
            let t0 = adj(ida_mem).ck_mem.last().unwrap().ck_t0;
            adj(ida_mem).dt_mem[0].t = t0;
            idaa_storePnt_dispatch(ida_mem, 0);
        } else {
            /* Load next point in dt_mem */
            let idx = (ida_mem.ida_nst % adj(ida_mem).ia_nsteps) as usize;
            let tn = ida_mem.ida_tn;
            adj(ida_mem).dt_mem[idx].t = tn;
            idaa_storePnt_dispatch(ida_mem, idx);
        }

        /* Set t1 field of the current check point structure
           for the case in which there will be no future
           check points */
        let tn = ida_mem.ida_tn;
        adj(ida_mem).ck_mem.last_mut().unwrap().ck_t1 = tn;

        /* tfinal is now set to tn */
        adj(ida_mem).ia_tfinal = tn;

        /* Return if in IDA_ONE_STEP mode */
        if itask == IDA_ONE_STEP {
            break;
        }

        /* IDA_NORMAL_STEP returns */

        /* Return if tout reached */
        if (*tret - tout) * ida_mem.ida_hh >= ZERO {
            /* If this was a root return, save the root time to return later */
            if flag == IDA_ROOT_RETURN {
                adj(ida_mem).ia_rootret = SUNTRUE;
                adj(ida_mem).ia_troot = *tret;
            }

            /* Get solution value at tout to return now */
            *tret = tout;
            flag = IDAGetSolution(ida_mem, tout, yret, ypret);

            /* Reset tretlast in IDA_mem so that IDAGetQuad and IDAGetSens
             * evaluate quadratures and/or sensitivities at the proper time */
            ida_mem.ida_tretlast = tout;

            break;
        }

        /* Return if tstop or a root was found */
        if flag == IDA_TSTOP_RETURN || flag == IDA_ROOT_RETURN {
            break;
        }
    } /* end of for(;;) */

    /* Get ncheck from IDAADJ_mem */
    let nst = ida_mem.ida_nst;
    let idaadj_mem = adj(ida_mem);
    *ncheckPtr = idaadj_mem.ia_nckpnts;

    /* Data is available for the last interval */
    idaadj_mem.ia_newData = SUNTRUE;
    idaadj_mem.ia_ckpntData = Some(idaadj_mem.ck_mem.len() - 1);
    idaadj_mem.ia_np = nst % idaadj_mem.ia_nsteps + 1;

    flag
}

/*
 * =================================================================
 * FUNCTIONS FOR BACKWARD PROBLEMS
 * =================================================================
 */

/* Common preamble for the ***B functions: ASA-initialized check and
   `which` validation, returning the Vec index of the backward problem
   (ida_index == position, creation order — see file header). */
fn idaa_which_index(ida_mem: &mut IDAMem, which: i32, fname: &str) -> Result<usize, i32> {
    /* Is ASA initialized ? */
    if !ida_mem.ida_adjMallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_ADJ, line!(), fname, file!(), MSGAM_NO_ADJ);
        return Err(IDA_NO_ADJ);
    }
    let idaadj_mem = ida_mem.ida_adj_mem.as_ref().unwrap();

    /* Check the value of which */
    if which >= idaadj_mem.ia_nbckpbs {
        IDAProcessError(None, IDA_ILL_INPUT, line!(), fname, file!(), MSGAM_BAD_WHICH);
        return Err(IDA_ILL_INPUT);
    }

    /* Find the IDABMem entry in the linked list corresponding to 'which'. */
    Ok(idaadj_mem.IDAB_mem.iter().position(|b| b.ida_index == which).unwrap())
}

pub fn IDACreateB(ida_mem: &mut IDAMem, which: &mut i32) -> i32 {
    /* Is ASA initialized ? */
    if !ida_mem.ida_adjMallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_ADJ, line!(), "IDACreateB", file!(), MSGAM_NO_ADJ);
        return IDA_NO_ADJ;
    }

    /* Allocate the IDAMem struct needed by this backward problem
       (infallible; the C malloc/NULL branches vanish). */
    let sunctx = ida_mem.ida_sunctx.clone();
    let mut ida_memB = IDACreate(&sunctx);

    /* We need to ensure Ns is set in the new IDAS object so that Ns is
       accessible in callbacks which only have access to ida_memB, not
       the original ida_mem */
    ida_memB.ida_Ns = ida_mem.ida_Ns;

    /* (C: IDASetUserData(ida_memB, ida_mem) — the outer-memory
       self-reference is created transiently by IDASolveB/IDACalcICB(S)
       around each nested call, per the pinned ownership dance.) */

    let idaadj_mem = adj(ida_mem);

    /* Initialize fields in the IDABMem struct and attach it (C prepends
       to the linked list; the Vec push keeps creation order). */
    let new_IDAB_mem = IDABMem {
        ida_index: idaadj_mem.ia_nbckpbs,
        ida_t0: ZERO,
        IDA_mem: ida_memB,
        ida_res_withSensi: SUNFALSE,
        ida_rhsQ_withSensi: SUNFALSE,
        ida_res: None,
        ida_resS: None,
        ida_rhsQ: None,
        ida_rhsQS: None,
        ida_user_data: None,
        ida_lmem: None,
        ida_pmem: None,
        ida_tout: ZERO,
        ida_yy: NVector::default(),
        ida_yp: NVector::default(),
    };
    idaadj_mem.IDAB_mem.push(new_IDAB_mem);

    /* Return the assigned index. This id is used as identificator and
       has to be passed to IDAInitB and other ***B functions that set the
       optional inputs for this backward problem. */
    *which = idaadj_mem.ia_nbckpbs;

    /* Increase the counter of the backward problems stored. */
    idaadj_mem.ia_nbckpbs += 1;

    IDA_SUCCESS
}

pub fn IDAInitB(
    ida_mem: &mut IDAMem,
    which: i32,
    resB: IDAResFnB,
    tB0: f64,
    yyB0: &NVector,
    ypB0: &NVector,
) -> i32 {
    /* Is ASA initialized ? */
    if !ida_mem.ida_adjMallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_ADJ, line!(), "IDAInitB", file!(), MSGAM_NO_ADJ);
        return IDA_NO_ADJ;
    }

    /* Check the initial time for this backward problem against the
       adjoint data. */
    {
        let idaadj_mem = ida_mem.ida_adj_mem.as_ref().unwrap();
        if tB0 < idaadj_mem.ia_tinitial || tB0 > idaadj_mem.ia_tfinal {
            IDAProcessError(Some(ida_mem), IDA_BAD_TB0, line!(), "IDAInitB", file!(),
                            MSGAM_BAD_TB0);
            return IDA_BAD_TB0;
        }
    }

    /* Check the value of which + find the IDABMem entry */
    let idx = match idaa_which_index(ida_mem, which, "IDAInitB") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let idaB_mem = &mut adj(ida_mem).IDAB_mem[idx];

    /* Call the IDAInit for this backward problem. */
    let flag = IDAInit(&mut idaB_mem.IDA_mem, IDAAres, tB0, yyB0, ypB0);
    if flag != IDA_SUCCESS {
        return flag;
    }

    /* Copy residual function in IDAB_mem. */
    idaB_mem.ida_res = Some(resB);
    idaB_mem.ida_res_withSensi = SUNFALSE;

    /* Initialized the initial time field. */
    idaB_mem.ida_t0 = tB0;

    /* Allocate and initialize space workspace vectors. */
    idaB_mem.ida_yy = N_VClone(yyB0);
    idaB_mem.ida_yp = N_VClone(yyB0);
    N_VScale(ONE, yyB0, &mut idaB_mem.ida_yy);
    N_VScale(ONE, ypB0, &mut idaB_mem.ida_yp);

    flag
}

pub fn IDAInitBS(
    ida_mem: &mut IDAMem,
    which: i32,
    resS: IDAResFnBS,
    tB0: f64,
    yyB0: &NVector,
    ypB0: &NVector,
) -> i32 {
    /* Is ASA initialized ? */
    if !ida_mem.ida_adjMallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_ADJ, line!(), "IDAInitBS", file!(), MSGAM_NO_ADJ);
        return IDA_NO_ADJ;
    }

    /* Check the initial time for this backward problem against the
       adjoint data. */
    {
        let idaadj_mem = ida_mem.ida_adj_mem.as_ref().unwrap();
        if tB0 < idaadj_mem.ia_tinitial || tB0 > idaadj_mem.ia_tfinal {
            IDAProcessError(Some(ida_mem), IDA_BAD_TB0, line!(), "IDAInitBS", file!(),
                            MSGAM_BAD_TB0);
            return IDA_BAD_TB0;
        }

        /* Were sensitivities active during the forward integration? */
        if !idaadj_mem.ia_storeSensi {
            IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDAInitBS", file!(),
                            MSGAM_BAD_SENSI);
            return IDA_ILL_INPUT;
        }
    }

    /* Check the value of which + find the IDABMem entry */
    let idx = match idaa_which_index(ida_mem, which, "IDAInitBS") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let idaB_mem = &mut adj(ida_mem).IDAB_mem[idx];

    /* Allocate and set the IDAS object */
    let flag = IDAInit(&mut idaB_mem.IDA_mem, IDAAres, tB0, yyB0, ypB0);
    if flag != IDA_SUCCESS {
        return flag;
    }

    /* Copy residual function pointer in IDAB_mem. */
    idaB_mem.ida_res_withSensi = SUNTRUE;
    idaB_mem.ida_resS = Some(resS);

    /* Allocate space and initialize the yy and yp vectors. */
    idaB_mem.ida_t0 = tB0;
    idaB_mem.ida_yy = N_VClone(yyB0);
    idaB_mem.ida_yp = N_VClone(ypB0);
    N_VScale(ONE, yyB0, &mut idaB_mem.ida_yy);
    N_VScale(ONE, ypB0, &mut idaB_mem.ida_yp);

    IDA_SUCCESS
}

pub fn IDAReInitB(
    ida_mem: &mut IDAMem,
    which: i32,
    tB0: f64,
    yyB0: &NVector,
    ypB0: &NVector,
) -> i32 {
    /* Is ASA initialized ? */
    if !ida_mem.ida_adjMallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_ADJ, line!(), "IDAReInitB", file!(), MSGAM_NO_ADJ);
        return IDA_NO_ADJ;
    }

    /* Check the initial time for this backward problem against the
       adjoint data. */
    {
        let idaadj_mem = ida_mem.ida_adj_mem.as_ref().unwrap();
        if tB0 < idaadj_mem.ia_tinitial || tB0 > idaadj_mem.ia_tfinal {
            IDAProcessError(Some(ida_mem), IDA_BAD_TB0, line!(), "IDAReInitB", file!(),
                            MSGAM_BAD_TB0);
            return IDA_BAD_TB0;
        }
    }

    /* Check the value of which + find the IDABMem entry */
    let idx = match idaa_which_index(ida_mem, which, "IDAReInitB") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let idaB_mem = &mut adj(ida_mem).IDAB_mem[idx];

    /* Call the IDAReInit for this backward problem. */
    IDAReInit(&mut idaB_mem.IDA_mem, tB0, yyB0, ypB0)
}

pub fn IDASStolerancesB(ida_mem: &mut IDAMem, which: i32, relTolB: f64, absTolB: f64) -> i32 {
    let idx = match idaa_which_index(ida_mem, which, "IDASStolerancesB") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let idaB_mem = &mut adj(ida_mem).IDAB_mem[idx];

    /* Set tolerances and return. */
    IDASStolerances(&mut idaB_mem.IDA_mem, relTolB, absTolB)
}

pub fn IDASVtolerancesB(ida_mem: &mut IDAMem, which: i32, relTolB: f64, absTolB: &NVector) -> i32 {
    let idx = match idaa_which_index(ida_mem, which, "IDASVtolerancesB") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let idaB_mem = &mut adj(ida_mem).IDAB_mem[idx];

    /* Set tolerances and return. */
    IDASVtolerances(&mut idaB_mem.IDA_mem, relTolB, absTolB)
}

pub fn IDAQuadSStolerancesB(ida_mem: &mut IDAMem, which: i32, reltolQB: f64, abstolQB: f64) -> i32 {
    let idx = match idaa_which_index(ida_mem, which, "IDAQuadSStolerancesB") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let idaB_mem = &mut adj(ida_mem).IDAB_mem[idx];

    IDAQuadSStolerances(&mut idaB_mem.IDA_mem, reltolQB, abstolQB)
}

pub fn IDAQuadSVtolerancesB(
    ida_mem: &mut IDAMem,
    which: i32,
    reltolQB: f64,
    abstolQB: &NVector,
) -> i32 {
    let idx = match idaa_which_index(ida_mem, which, "IDAQuadSVtolerancesB") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let idaB_mem = &mut adj(ida_mem).IDAB_mem[idx];

    IDAQuadSVtolerances(&mut idaB_mem.IDA_mem, reltolQB, abstolQB)
}

pub fn IDAQuadInitB(ida_mem: &mut IDAMem, which: i32, rhsQB: IDAQuadRhsFnB, yQB0: &NVector) -> i32 {
    let idx = match idaa_which_index(ida_mem, which, "IDAQuadInitB") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let idaB_mem = &mut adj(ida_mem).IDAB_mem[idx];

    let flag = IDAQuadInit(&mut idaB_mem.IDA_mem, IDAArhsQ, yQB0);
    if flag != IDA_SUCCESS {
        return flag;
    }

    idaB_mem.ida_rhsQ_withSensi = SUNFALSE;
    idaB_mem.ida_rhsQ = Some(rhsQB);

    flag
}

pub fn IDAQuadInitBS(
    ida_mem: &mut IDAMem,
    which: i32,
    rhsQS: IDAQuadRhsFnBS,
    yQB0: &NVector,
) -> i32 {
    let idx = match idaa_which_index(ida_mem, which, "IDAQuadInitBS") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let idaB_mem = &mut adj(ida_mem).IDAB_mem[idx];

    /* Allocate and set the IDAS object */
    let flag = IDAQuadInit(&mut idaB_mem.IDA_mem, IDAArhsQ, yQB0);
    if flag != IDA_SUCCESS {
        return flag;
    }

    /* Copy RHS function pointer in IDAB_mem and enable quad
       sensitivities. */
    idaB_mem.ida_rhsQ_withSensi = SUNTRUE;
    idaB_mem.ida_rhsQS = Some(rhsQS);

    IDA_SUCCESS
}

pub fn IDAQuadReInitB(ida_mem: &mut IDAMem, which: i32, yQB0: &NVector) -> i32 {
    let idx = match idaa_which_index(ida_mem, which, "IDAQuadReInitB") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let idaB_mem = &mut adj(ida_mem).IDAB_mem[idx];

    IDAQuadReInit(&mut idaB_mem.IDA_mem, yQB0)
}

/*
 * ----------------------------------------------------------------
 * Function : IDACalcICB
 * ----------------------------------------------------------------
 * IDACalcIC calculates corrected initial conditions for a DAE
 * backward system (index-one in semi-implicit form).
 * It uses Newton iteration combined with a Linesearch algorithm.
 * Calling IDACalcICB is optional. It is only necessary when the
 * initial conditions do not solve the given system.  I.e., if
 * yB0 and ypB0 are known to satisfy the backward problem, then
 * a call to IDACalcIC is NOT necessary (for index-one problems).
 */
pub fn IDACalcICB(
    ida_mem: &mut IDAMem,
    which: i32,
    tout1: f64,
    yy0: &NVector,
    yp0: &NVector,
) -> i32 {
    let idx = match idaa_which_index(ida_mem, which, "IDACalcICB") {
        Ok(i) => i,
        Err(e) => return e,
    };

    {
        let idaadj_mem = adj(ida_mem);

        /* The wrapper for user supplied res function requires ia_bckpbCrt
           from IDAAdjMem to be set to current problem. */
        idaadj_mem.ia_bckpbCrt = Some(idx);

        /* Save (y, y') in yyTmp and ypTmp for use in the res wrapper. */
        /* yyTmp and ypTmp workspaces are safe to use if IDAADataStore is
           not called. */
        N_VScale(ONE, yy0, &mut idaadj_mem.ia_yyTmp);
        N_VScale(ONE, yp0, &mut idaadj_mem.ia_ypTmp);

        /* Set noInterp flag to SUNTRUE, so IDAARes will use user provided
           values for y and y' and will not call the interpolation
           routine(s). */
        idaadj_mem.ia_noInterp = SUNTRUE;
    }

    let flag = idaa_calc_ic_backward(ida_mem, idx, tout1);

    /* Set interpolation on in IDAARes. */
    adj(ida_mem).ia_noInterp = SUNFALSE;

    flag
}

/*
 * ----------------------------------------------------------------
 * Function : IDACalcICBS
 * ----------------------------------------------------------------
 * IDACalcIC calculates corrected initial conditions for a DAE
 * backward system (index-one in semi-implicit form) that also
 * dependes on the sensivities.
 *
 * It calls IDACalcIC for the 'which' backward problem.
 */
pub fn IDACalcICBS(
    ida_mem: &mut IDAMem,
    which: i32,
    tout1: f64,
    yy0: &NVector,
    yp0: &NVector,
    yyS0: &[NVector],
    ypS0: &[NVector],
) -> i32 {
    /* Is ASA initialized? */
    if !ida_mem.ida_adjMallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_ADJ, line!(), "IDACalcICBS", file!(), MSGAM_NO_ADJ);
        return IDA_NO_ADJ;
    }

    /* Were sensitivities active during the forward integration? */
    if !adj(ida_mem).ia_storeSensi {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDACalcICBS", file!(),
                        MSGAM_BAD_SENSI);
        return IDA_ILL_INPUT;
    }

    let idx = match idaa_which_index(ida_mem, which, "IDACalcICBS") {
        Ok(i) => i,
        Err(e) => return e,
    };

    /* Was InitBS called for this problem? */
    if !adj(ida_mem).IDAB_mem[idx].ida_res_withSensi {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDACalcICBS", file!(),
                        MSGAM_NO_INITBS);
        return IDA_ILL_INPUT;
    }

    {
        let ns = ida_mem.ida_Ns as usize;
        let idaadj_mem = ida_mem.ida_adj_mem.as_mut().unwrap();

        /* The wrapper for user supplied res function requires ia_bckpbCrt
           from IDAAdjMem to be set to current problem. */
        idaadj_mem.ia_bckpbCrt = Some(idx);

        /* Save (y, y') and (y_p, y'_p) in yyTmp, ypTmp and yySTmp, ypSTmp.
           The wrapper for residual will use these values instead of
           calling interpolation routine. */
        /* The four workspaces variables are safe to use if IDAADataStore
           is not called. */
        N_VScale(ONE, yy0, &mut idaadj_mem.ia_yyTmp);
        N_VScale(ONE, yp0, &mut idaadj_mem.ia_ypTmp);

        /* (C: cvals[is]=ONE + N_VScaleVectorArray(Ns, cvals, yyS0/ypS0,
           yySTmp/ypSTmp); serial fused kernel reproduced inline; the
           cannot-fail IDA_VECTOROP_ERR branches vanish.) */
        for is in 0..ns {
            N_VScale(ONE, &yyS0[is], &mut idaadj_mem.ia_yySTmp[is]);
        }
        for is in 0..ns {
            N_VScale(ONE, &ypS0[is], &mut idaadj_mem.ia_ypSTmp[is]);
        }

        /* Set noInterp flag to SUNTRUE, so IDAARes will use user provided
           values for y and y' and will not call the interpolation
           routine(s). */
        idaadj_mem.ia_noInterp = SUNTRUE;
    }

    let flag = idaa_calc_ic_backward(ida_mem, idx, tout1);

    /* Set interpolation on in IDAARes. */
    adj(ida_mem).ia_noInterp = SUNFALSE;

    flag
}

/*
 * IDASolveB
 *
 * This routine performs the backward integration from tB0
 * to tinitial through a sequence of forward-backward runs in
 * between consecutive check points. It returns the values of
 * the adjoint variables and any existing quadrature variables
 * at tinitial.
 *
 * On a successful return, IDASolveB returns IDA_SUCCESS.
 *
 * NOTE that IDASolveB DOES NOT return the solution for the
 * backward problem(s). Use IDAGetB to extract the solution
 * for any given backward problem.
 *
 * If there are multiple backward problems and multiple check points,
 * IDASolveB may not succeed in getting all problems to take one step
 * when called in ONE_STEP mode.
 */
pub fn IDASolveB(ida_mem: &mut IDAMem, mut tBout: f64, itaskB: i32) -> i32 {
    /* Is ASA initialized ? */
    if !ida_mem.ida_adjMallocDone {
        IDAProcessError(Some(ida_mem), IDA_NO_ADJ, line!(), "IDASolveB", file!(), MSGAM_NO_ADJ);
        return IDA_NO_ADJ;
    }

    if adj(ida_mem).ia_nbckpbs == 0 {
        IDAProcessError(Some(ida_mem), IDA_NO_BCK, line!(), "IDASolveB", file!(), MSGAM_NO_BCK);
        return IDA_NO_BCK;
    }

    /* Check whether IDASolveF has been called */
    if adj(ida_mem).ia_firstIDAFcall {
        IDAProcessError(Some(ida_mem), IDA_NO_FWD, line!(), "IDASolveB", file!(), MSGAM_NO_FWD);
        return IDA_NO_FWD;
    }
    let sign: f64 = {
        let idaadj_mem = ida_mem.ida_adj_mem.as_ref().unwrap();
        if idaadj_mem.ia_tfinal - idaadj_mem.ia_tinitial > ZERO { 1.0 } else { -1.0 }
    };

    /* If this is the first call, loop over all backward problems and
     *   - check that tB0 is valid
     *   - check that tBout is ahead of tB0 in the backward direction
     *   - check whether we need to interpolate forward sensitivities
     */
    if adj(ida_mem).ia_firstIDABcall {
        let (tinitial, tfinal) = {
            let idaadj_mem = ida_mem.ida_adj_mem.as_ref().unwrap();
            (idaadj_mem.ia_tinitial, idaadj_mem.ia_tfinal)
        };
        let nbck = adj(ida_mem).IDAB_mem.len();
        let mut interp_sensi = adj(ida_mem).ia_interpSensi;
        let mut bad: Option<(i32, i32)> = None; /* (error code, problem index) */
        {
            let idaadj_mem = ida_mem.ida_adj_mem.as_ref().unwrap();
            for idx in 0..nbck {
                let b = &idaadj_mem.IDAB_mem[idx];
                let tBn = b.IDA_mem.ida_tn;

                if sign * (tBn - tinitial) < ZERO || sign * (tfinal - tBn) < ZERO {
                    bad = Some((IDA_BAD_TB0, b.ida_index));
                    break;
                }

                if sign * (tBn - tBout) <= ZERO {
                    bad = Some((IDA_ILL_INPUT, b.ida_index));
                    break;
                }

                if b.ida_res_withSensi || b.ida_rhsQ_withSensi {
                    interp_sensi = SUNTRUE;
                }
            }
        }
        if let Some((code, _index)) = bad {
            /* (C passes the problem index as a varargs extra; the message
               strings carry no %d conversion, so it is not printed) */
            let msg = if code == IDA_BAD_TB0 { MSGAM_BAD_TB0 } else { MSGAM_BAD_TBOUT };
            IDAProcessError(Some(ida_mem), code, line!(), "IDASolveB", file!(), msg);
            return code;
        }
        adj(ida_mem).ia_interpSensi = interp_sensi;

        if adj(ida_mem).ia_interpSensi && !adj(ida_mem).ia_storeSensi {
            IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASolveB", file!(),
                            MSGAM_BAD_SENSI);
            return IDA_ILL_INPUT;
        }

        adj(ida_mem).ia_firstIDABcall = SUNFALSE;
    }

    /* Check for valid itask */
    if itaskB != IDA_NORMAL && itaskB != IDA_ONE_STEP {
        IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASolveB", file!(),
                        MSG_BAD_ITASK);
        return IDA_ILL_INPUT;
    }

    /* Check if tBout is legal */
    let (tinitial, tfinal) = {
        let idaadj_mem = ida_mem.ida_adj_mem.as_ref().unwrap();
        (idaadj_mem.ia_tinitial, idaadj_mem.ia_tfinal)
    };
    if sign * (tBout - tinitial) < ZERO || sign * (tfinal - tBout) < ZERO {
        let tfuzz = HUNDRED * ida_mem.ida_uround * (SUNRabs(tinitial) + SUNRabs(tfinal));
        if sign * (tBout - tinitial) < ZERO && SUNRabs(tBout - tinitial) < tfuzz {
            tBout = tinitial;
        } else {
            IDAProcessError(Some(ida_mem), IDA_ILL_INPUT, line!(), "IDASolveB", file!(),
                            MSGAM_BAD_TBOUT);
            return IDA_ILL_INPUT;
        }
    }

    /* Loop through the check points and stop as soon as a backward
     * problem has its tn value behind the current check point's t0_
     * value (in the backward direction) */

    let mut ck_idx = adj(ida_mem).ck_mem.len() - 1; /* head = most recent */

    {
        let idaadj_mem = ida_mem.ida_adj_mem.as_ref().unwrap();
        let mut gotCkpnt = SUNFALSE;
        loop {
            let ck_t0 = idaadj_mem.ck_mem[ck_idx].ck_t0;
            for b in idaadj_mem.IDAB_mem.iter() {
                let tBn = b.IDA_mem.ida_tn;

                if sign * (tBn - ck_t0) > ZERO {
                    gotCkpnt = SUNTRUE;
                    break;
                }

                if itaskB == IDA_NORMAL && tBn == ck_t0 && sign * (tBout - ck_t0) >= ZERO {
                    gotCkpnt = SUNTRUE;
                    break;
                }
            }

            if gotCkpnt {
                break;
            }

            if ck_idx == 0 {
                /* ck_next == NULL */
                break;
            }

            ck_idx -= 1;
        }
    }

    /* Loop while propagating backward problems */
    let mut flag = 0;
    loop {
        /* Store interpolation data if not available.
           This is the 2nd forward integration pass */
        if adj(ida_mem).ia_ckpntData != Some(ck_idx) {
            flag = IDAAdataStore(ida_mem, ck_idx);
            if flag != IDA_SUCCESS {
                break;
            }
        }

        /* Starting with the current check point from above, loop over
           check points while propagating backward problems */

        let nbck = adj(ida_mem).IDAB_mem.len();
        let mut err_index = 0;
        for idx in 0..nbck {
            /* Decide if current backward problem is "active" in this
               check point */
            let mut isActive = SUNTRUE;

            let ck_t0 = adj(ida_mem).ck_mem[ck_idx].ck_t0;
            let tBn = adj(ida_mem).IDAB_mem[idx].IDA_mem.ida_tn;

            if tBn == ck_t0 && sign * (tBout - ck_t0) < ZERO {
                isActive = SUNFALSE;
            }
            if tBn == ck_t0 && itaskB == IDA_ONE_STEP {
                isActive = SUNFALSE;
            }
            if sign * (tBn - ck_t0) < ZERO {
                isActive = SUNFALSE;
            }

            if isActive {
                /* Store the address of current backward problem memory
                 * in IDAADJ_mem to be used in the wrapper functions */
                adj(ida_mem).ia_bckpbCrt = Some(idx);

                /* Integrate current backward problem */
                IDASetStopTime(&mut adj(ida_mem).IDAB_mem[idx].IDA_mem, ck_t0);
                let mut tBret = ZERO;
                flag = idaa_integrate_backward(ida_mem, idx, tBout, itaskB, &mut tBret);

                /* Set the time at which we will report solution and/or
                   quadratures */
                adj(ida_mem).IDAB_mem[idx].ida_tout = tBret;

                /* If an error occurred, exit while loop */
                if flag < 0 {
                    err_index = adj(ida_mem).IDAB_mem[idx].ida_index;
                    break;
                }
            } else {
                flag = IDA_SUCCESS;
                adj(ida_mem).IDAB_mem[idx].ida_tout = tBn;
            }
        } /* End of while: iteration through backward problems. */

        /* If an error occurred, return now */
        if flag < 0 {
            IDAProcessError(Some(ida_mem), flag, line!(), "IDASolveB", file!(),
                            &format!("Error occurred while integrating backward problem # {}",
                                     err_index));
            return flag;
        }

        /* If in IDA_ONE_STEP mode, return now (flag = IDA_SUCCESS) */
        if itaskB == IDA_ONE_STEP {
            break;
        }

        /* If all backward problems have successfully reached tBout,
           return now */
        let mut reachedTBout = SUNTRUE;
        {
            let idaadj_mem = ida_mem.ida_adj_mem.as_ref().unwrap();
            for b in idaadj_mem.IDAB_mem.iter() {
                if sign * (b.ida_tout - tBout) > ZERO {
                    reachedTBout = SUNFALSE;
                    break;
                }
            }
        }

        if reachedTBout {
            break;
        }

        /* Move check point in linked list to next one */
        ck_idx -= 1;
    } /* End of loop. */

    flag
}

/* Integrate backward problem `idx` towards tBout: the ownership dance
   that realizes the C's permanent user_data == forward ida_mem link
   (see the note at the top of this file). */
fn idaa_integrate_backward(
    ida_mem: &mut IDAMem,
    idx: usize,
    tBout: f64,
    itaskB: i32,
    tBret: &mut f64,
) -> i32 {
    let sunctx = ida_mem.ida_sunctx.clone();

    /* detach the backward problem's solver memory and workspace vectors */
    let (mut nested, mut yy, mut yp) = {
        let idaadj_mem = adj(ida_mem);
        let b = &mut idaadj_mem.IDAB_mem[idx];
        let nested = std::mem::replace(&mut b.IDA_mem, IDACreate(&sunctx));
        let yy = std::mem::take(&mut b.ida_yy);
        let yp = std::mem::take(&mut b.ida_yp);
        (nested, yy, yp)
    };

    /* move the forward memory into the backward problem's user_data */
    let outer = std::mem::replace(ida_mem, *IDACreate(&sunctx));
    nested.ida_user_data = Some(Box::new(outer));

    let flag = IDASolve(&mut nested, tBout, tBret, &mut yy, &mut yp, itaskB);

    /* restore the forward memory and reattach the backward problem */
    let outer = nested.ida_user_data.take().unwrap().downcast::<IDAMem>().unwrap();
    *ida_mem = *outer;

    let idaadj_mem = adj(ida_mem);
    let b = &mut idaadj_mem.IDAB_mem[idx];
    b.IDA_mem = nested;
    b.ida_yy = yy;
    b.ida_yp = yp;

    flag
}

/* IDACalcIC on backward problem `idx`: same ownership dance as
   idaa_integrate_backward (IDAAres reaches the forward memory through
   the nested problem's user_data). */
fn idaa_calc_ic_backward(ida_mem: &mut IDAMem, idx: usize, tout1: f64) -> i32 {
    let sunctx = ida_mem.ida_sunctx.clone();

    /* detach the backward problem's solver memory */
    let mut nested = {
        let idaadj_mem = adj(ida_mem);
        std::mem::replace(&mut idaadj_mem.IDAB_mem[idx].IDA_mem, IDACreate(&sunctx))
    };

    /* move the forward memory into the backward problem's user_data */
    let outer = std::mem::replace(ida_mem, *IDACreate(&sunctx));
    nested.ida_user_data = Some(Box::new(outer));

    let flag = IDACalcIC(&mut nested, IDA_YA_YDP_INIT, tout1);

    /* restore the forward memory and reattach the backward problem */
    let outer = nested.ida_user_data.take().unwrap().downcast::<IDAMem>().unwrap();
    *ida_mem = *outer;
    adj(ida_mem).IDAB_mem[idx].IDA_mem = nested;

    flag
}

/*
 * IDAGetB
 *
 * IDAGetB returns the state variables at the same time (also returned
 * in tret) as that at which IDASolveB returned the solution.
 */
pub fn IDAGetB(
    ida_mem: &mut IDAMem,
    which: i32,
    tret: &mut f64,
    yy: &mut NVector,
    yp: &mut NVector,
) -> i32 {
    let idx = match idaa_which_index(ida_mem, which, "IDAGetB") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let idaB_mem = &adj(ida_mem).IDAB_mem[idx];

    N_VScale(ONE, &idaB_mem.ida_yy, yy);
    N_VScale(ONE, &idaB_mem.ida_yp, yp);
    *tret = idaB_mem.ida_tout;

    IDA_SUCCESS
}

/*
 * IDAGetQuadB
 *
 * IDAGetQuadB returns the quadrature variables at the same
 * time (also returned in tret) as that at which IDASolveB
 * returned the solution.
 */
pub fn IDAGetQuadB(ida_mem: &mut IDAMem, which: i32, tret: &mut f64, qB: &mut NVector) -> i32 {
    let idx = match idaa_which_index(ida_mem, which, "IDAGetQuadB") {
        Ok(i) => i,
        Err(e) => return e,
    };
    let idaB_mem = &adj(ida_mem).IDAB_mem[idx];

    /* If the integration for this backward problem has not started yet,
     * simply return the current value of qB (i.e. the final conditions) */

    let nstB = idaB_mem.IDA_mem.ida_nst; /* C: IDAGetNumSteps */

    if nstB == 0 {
        N_VScale(ONE, &idaB_mem.IDA_mem.ida_phiQ[0], qB);
        *tret = idaB_mem.ida_tout;
        IDA_SUCCESS
    } else {
        IDAGetQuad(&idaB_mem.IDA_mem, tret, qB)
    }
}

/*=================================================================*/
/*                Private Functions Implementation                 */
/*=================================================================*/

/*
 * IDAAckpntInit
 *
 * This routine initializes the check point linked list with
 * information from the initial time.
 */
fn IDAAckpntInit(ida_mem: &mut IDAMem) -> IDAckpntMem {
    let mut ck_mem = IDAckpntMem {
        ck_t0: ida_mem.ida_tn,
        ck_t1: ZERO,
        ck_nst: 0,
        ck_kk: 1,
        ck_hh: ZERO,

        /* Test if we need to carry quadratures */
        ck_quadr: ida_mem.ida_quadr && ida_mem.ida_errconQ,

        /* Test if we need to carry sensitivities */
        ck_sensi: ida_mem.ida_sensi,
        ck_Ns: if ida_mem.ida_sensi { ida_mem.ida_Ns } else { 0 },

        /* Test if we need to carry quadrature sensitivities */
        ck_quadr_sensi: ida_mem.ida_quadr_sensi && ida_mem.ida_errconQS,

        /* Alloc 3: current order, i.e. 1,  +   2. */
        ck_phi_alloc: 3,

        /* (remaining step data is loaded by IDAAckpntNew only; the C
           leaves these fields uninitialized for the initial check
           point) */
        ck_tretlast: ZERO,
        ck_ns: 0,
        ck_kused: 0,
        ck_knew: 0,
        ck_phase: 0,
        ck_hused: ZERO,
        ck_eta: ZERO,
        ck_cj: ZERO,
        ck_cjlast: ZERO,
        ck_cjold: ZERO,
        ck_cjratio: ZERO,
        ck_ss: ZERO,
        ck_ssS: ZERO,
        ck_psi: [ZERO; MXORDP1],
        ck_alpha: [ZERO; MXORDP1],
        ck_beta: [ZERO; MXORDP1],
        ck_sigma: [ZERO; MXORDP1],
        ck_gamma: [ZERO; MXORDP1],

        ck_phi: Vec::new(),
        ck_phiQ: Vec::new(),
        ck_phiS: std::array::from_fn(|_| Vec::new()),
        ck_phiQS: std::array::from_fn(|_| Vec::new()),
    };

    /* (IDAAckpntAllocVectors is infallible here; the C rollback/NULL
       branches vanish) */
    IDAAckpntAllocVectors(ida_mem, &mut ck_mem);
    /* Save phi* vectors from IDA_mem to ck_mem. */
    IDAAckpntCopyVectors(ida_mem, &mut ck_mem);

    ck_mem
}

/*
 * IDAAckpntNew
 *
 * This routine allocates space for a new check point and sets
 * its data from current values in IDA_mem.
 */
fn IDAAckpntNew(ida_mem: &mut IDAMem) -> IDAckpntMem {
    let mut ck_mem = IDAckpntMem {
        ck_nst: ida_mem.ida_nst,
        ck_tretlast: ida_mem.ida_tretlast,
        ck_kk: ida_mem.ida_kk,
        ck_kused: ida_mem.ida_kused,
        ck_knew: ida_mem.ida_knew,
        ck_phase: ida_mem.ida_phase,
        ck_ns: ida_mem.ida_ns,
        ck_hh: ida_mem.ida_hh,
        ck_hused: ida_mem.ida_hused,
        ck_eta: ida_mem.ida_eta,
        ck_cj: ida_mem.ida_cj,
        ck_cjlast: ida_mem.ida_cjlast,
        ck_cjold: ida_mem.ida_cjold,
        ck_cjratio: ida_mem.ida_cjratio,
        ck_ss: ida_mem.ida_ss,
        ck_ssS: ida_mem.ida_ssS,
        ck_t0: ida_mem.ida_tn,
        ck_t1: ZERO,

        ck_psi: ida_mem.ida_psi,
        ck_alpha: ida_mem.ida_alpha,
        ck_beta: ida_mem.ida_beta,
        ck_sigma: ida_mem.ida_sigma,
        ck_gamma: ida_mem.ida_gamma,

        /* Test if we need to carry quadratures */
        ck_quadr: ida_mem.ida_quadr && ida_mem.ida_errconQ,

        /* Test if we need to carry sensitivities */
        ck_sensi: ida_mem.ida_sensi,
        ck_Ns: if ida_mem.ida_sensi { ida_mem.ida_Ns } else { 0 },

        /* Test if we need to carry quadrature sensitivities */
        ck_quadr_sensi: ida_mem.ida_quadr_sensi && ida_mem.ida_errconQS,

        ck_phi_alloc: if ida_mem.ida_kk + 2 < MXORDP1 as i32 {
            ida_mem.ida_kk + 2
        } else {
            MXORDP1 as i32
        },

        ck_phi: Vec::new(),
        ck_phiQ: Vec::new(),
        ck_phiS: std::array::from_fn(|_| Vec::new()),
        ck_phiQS: std::array::from_fn(|_| Vec::new()),
    };

    IDAAckpntAllocVectors(ida_mem, &mut ck_mem);
    /* Save phi* vectors from IDA_mem to ck_mem. */
    IDAAckpntCopyVectors(ida_mem, &mut ck_mem);

    ck_mem
}

/* (IDAAckpntDelete: RAII — dropping an IDAckpntMem releases its phi,
   phiQ, phiS and phiQS vectors.) */

/*
 * IDAAckpntAllocVectors
 *
 * Allocate checkpoint's phi, phiQ, phiS, phiQS vectors needed to save
 * current state of IDAMem.  (Allocation is infallible; the C rollback
 * branches vanish and no boolean is returned.)
 */
fn IDAAckpntAllocVectors(ida_mem: &mut IDAMem, ck_mem: &mut IDAckpntMem) {
    let alloc = ck_mem.ck_phi_alloc as usize;
    let ns = ida_mem.ida_Ns as usize;

    ck_mem.ck_phi = (0..alloc).map(|_| N_VClone(&ida_mem.ida_tempv1)).collect();

    /* Do we need to carry quadratures? */
    if ck_mem.ck_quadr {
        ck_mem.ck_phiQ = (0..alloc).map(|_| N_VClone(&ida_mem.ida_eeQ)).collect();
    }

    /* Do we need to carry sensitivities? */
    if ck_mem.ck_sensi {
        for j in 0..alloc {
            ck_mem.ck_phiS[j] = (0..ns).map(|_| N_VClone(&ida_mem.ida_tempv1)).collect();
        }
    }

    /* Do we need to carry quadrature sensitivities? */
    if ck_mem.ck_quadr_sensi {
        for j in 0..alloc {
            ck_mem.ck_phiQS[j] = (0..ns).map(|_| N_VClone(&ida_mem.ida_eeQ)).collect();
        }
    }
}

/*
 * IDAAckpntCopyVectors
 *
 * Copy phi* vectors from IDAMem in the corresponding vectors from
 * checkpoint.  (C stages cvals/Xvecs/Zvecs for fused
 * N_VScaleVectorArray copies; the serial fused kernel with unit
 * coefficients is a plain copy, reproduced as N_VScale loops.)
 */
fn IDAAckpntCopyVectors(ida_mem: &mut IDAMem, ck_mem: &mut IDAckpntMem) {
    let alloc = ck_mem.ck_phi_alloc as usize;
    let ns = ida_mem.ida_Ns as usize;

    /* Save phi* arrays from IDA_mem */

    for j in 0..alloc {
        N_VScale(ONE, &ida_mem.ida_phi[j], &mut ck_mem.ck_phi[j]);
    }

    if ck_mem.ck_quadr {
        for j in 0..alloc {
            N_VScale(ONE, &ida_mem.ida_phiQ[j], &mut ck_mem.ck_phiQ[j]);
        }
    }

    if ck_mem.ck_sensi {
        for j in 0..alloc {
            for is in 0..ns {
                N_VScale(ONE, &ida_mem.ida_phiS[j][is], &mut ck_mem.ck_phiS[j][is]);
            }
        }
    }

    if ck_mem.ck_quadr_sensi {
        for j in 0..alloc {
            for is in 0..ns {
                N_VScale(ONE, &ida_mem.ida_phiQS[j][is], &mut ck_mem.ck_phiQS[j][is]);
            }
        }
    }
}

/* (IDAAdataMalloc is folded into IDAAdjInit — the dt_mem Vec is built
   there; IDAAdataFree is RAII: dropping IDAadjMem releases dt_mem and
   the interpolation content.) */

/*
 * IDAAdataStore
 *
 * This routine integrates the forward model starting at the check
 * point ck_mem and stores y and yprime at all intermediate
 * steps.
 *
 * Return values:
 *   - the flag that IDASolve may return on error
 *   - IDA_REIFWD_FAIL if no check point is available for this hot start
 *   - IDA_SUCCESS
 */
fn IDAAdataStore(ida_mem: &mut IDAMem, ck_idx: usize) -> i32 {
    /* Initialize IDA_mem with data from ck_mem. */
    let flag = IDAAckpntGet(ida_mem, ck_idx);
    if flag != IDA_SUCCESS {
        return IDA_REIFWD_FAIL;
    }

    /* Set first structure in dt_mem[0] */
    {
        let idaadj_mem = adj(ida_mem);
        idaadj_mem.dt_mem[0].t = idaadj_mem.ck_mem[ck_idx].ck_t0;
    }
    idaa_storePnt_dispatch(ida_mem, 0);

    /* Decide whether TSTOP must be activated */
    if adj(ida_mem).ia_tstopIDAFcall {
        let tstop = adj(ida_mem).ia_tstopIDAF;
        IDASetStopTime(ida_mem, tstop);
    }

    let (tinitial, tfinal, ck_t1) = {
        let idaadj_mem = ida_mem.ida_adj_mem.as_ref().unwrap();
        (idaadj_mem.ia_tinitial, idaadj_mem.ia_tfinal, idaadj_mem.ck_mem[ck_idx].ck_t1)
    };
    let sign: f64 = if tfinal - tinitial > ZERO { 1.0 } else { -1.0 };

    /* Run IDASolve in IDA_ONE_STEP mode to set following structures in
       dt_mem[i]. */
    let mut i: i64 = 1;
    loop {
        let mut yyTmp = std::mem::take(&mut adj(ida_mem).ia_yyTmp);
        let mut ypTmp = std::mem::take(&mut adj(ida_mem).ia_ypTmp);
        let mut t = ZERO;
        let flag = IDASolve(ida_mem, ck_t1, &mut t, &mut yyTmp, &mut ypTmp, IDA_ONE_STEP);
        {
            let idaadj_mem = adj(ida_mem);
            idaadj_mem.ia_yyTmp = yyTmp;
            idaadj_mem.ia_ypTmp = ypTmp;
        }
        if flag < 0 {
            return IDA_FWD_FAIL;
        }

        adj(ida_mem).dt_mem[i as usize].t = t;
        idaa_storePnt_dispatch(ida_mem, i as usize);

        i += 1;

        if sign * (ck_t1 - t) <= ZERO {
            break;
        }
    }

    /* New data is now available. */
    let idaadj_mem = adj(ida_mem);
    idaadj_mem.ia_ckpntData = Some(ck_idx);
    idaadj_mem.ia_newData = SUNTRUE;
    idaadj_mem.ia_np = i;

    IDA_SUCCESS
}

/*
 * IDAAckpntGet  (C comment header says CVAckpntGet)
 *
 * This routine prepares IDAS for a hot restart from
 * the check point ck_mem
 */
fn IDAAckpntGet(ida_mem: &mut IDAMem, ck_idx: usize) -> i32 {
    /* Detach the adjoint memory so the check point data can be borrowed
       while ida_mem is reinitialized */
    let idaadj_mem = ida_mem.ida_adj_mem.take().unwrap();
    let flag = idaa_ckpnt_get_inner(ida_mem, &idaadj_mem, ck_idx);
    ida_mem.ida_adj_mem = Some(idaadj_mem);
    flag
}

fn idaa_ckpnt_get_inner(ida_mem: &mut IDAMem, idaadj_mem: &IDAadjMem, ck_idx: usize) -> i32 {
    let ck_mem = &idaadj_mem.ck_mem[ck_idx];

    if ck_idx == 0 {
        /* ck_next == NULL: this is the check point at the initial time.
         * In this case, we just call the reinitialization routine,
         * but make sure we use the same initial stepsize as on
         * the first run. */

        let h0u = ida_mem.ida_h0u;
        IDASetInitStep(ida_mem, h0u);

        let flag = IDAReInit(ida_mem, ck_mem.ck_t0, &ck_mem.ck_phi[0], &ck_mem.ck_phi[1]);
        if flag != IDA_SUCCESS {
            return flag;
        }

        if ck_mem.ck_quadr {
            let flag = IDAQuadReInit(ida_mem, &ck_mem.ck_phiQ[0]);
            if flag != IDA_SUCCESS {
                return flag;
            }
        }

        if ck_mem.ck_sensi {
            let ism = ida_mem.ida_ism;
            let flag = IDASensReInit(ida_mem, ism, &ck_mem.ck_phiS[0], &ck_mem.ck_phiS[1]);
            if flag != IDA_SUCCESS {
                return flag;
            }
        }

        if ck_mem.ck_quadr_sensi {
            let flag = IDAQuadSensReInit(ida_mem, &ck_mem.ck_phiQS[0]);
            if flag != IDA_SUCCESS {
                return flag;
            }
        }
    } else {
        /* Copy parameters from check point data structure */
        ida_mem.ida_nst = ck_mem.ck_nst;
        ida_mem.ida_tretlast = ck_mem.ck_tretlast;
        ida_mem.ida_kk = ck_mem.ck_kk;
        ida_mem.ida_kused = ck_mem.ck_kused;
        ida_mem.ida_knew = ck_mem.ck_knew;
        ida_mem.ida_phase = ck_mem.ck_phase;
        ida_mem.ida_ns = ck_mem.ck_ns;
        ida_mem.ida_hh = ck_mem.ck_hh;
        ida_mem.ida_hused = ck_mem.ck_hused;
        ida_mem.ida_eta = ck_mem.ck_eta;
        ida_mem.ida_cj = ck_mem.ck_cj;
        ida_mem.ida_cjlast = ck_mem.ck_cjlast;
        ida_mem.ida_cjold = ck_mem.ck_cjold;
        ida_mem.ida_cjratio = ck_mem.ck_cjratio;
        ida_mem.ida_tn = ck_mem.ck_t0;
        ida_mem.ida_ss = ck_mem.ck_ss;
        ida_mem.ida_ssS = ck_mem.ck_ssS;

        /* Copy the arrays from check point data structure */
        for j in 0..ck_mem.ck_phi_alloc as usize {
            N_VScale(ONE, &ck_mem.ck_phi[j], &mut ida_mem.ida_phi[j]);
        }

        if ck_mem.ck_quadr {
            for j in 0..ck_mem.ck_phi_alloc as usize {
                N_VScale(ONE, &ck_mem.ck_phiQ[j], &mut ida_mem.ida_phiQ[j]);
            }
        }

        if ck_mem.ck_sensi {
            for is in 0..ida_mem.ida_Ns as usize {
                for j in 0..ck_mem.ck_phi_alloc as usize {
                    N_VScale(ONE, &ck_mem.ck_phiS[j][is], &mut ida_mem.ida_phiS[j][is]);
                }
            }
        }

        if ck_mem.ck_quadr_sensi {
            for is in 0..ida_mem.ida_Ns as usize {
                for j in 0..ck_mem.ck_phi_alloc as usize {
                    N_VScale(ONE, &ck_mem.ck_phiQS[j][is], &mut ida_mem.ida_phiQS[j][is]);
                }
            }
        }

        for j in 0..MXORDP1 {
            ida_mem.ida_psi[j] = ck_mem.ck_psi[j];
            ida_mem.ida_alpha[j] = ck_mem.ck_alpha[j];
            ida_mem.ida_beta[j] = ck_mem.ck_beta[j];
            ida_mem.ida_sigma[j] = ck_mem.ck_sigma[j];
            ida_mem.ida_gamma[j] = ck_mem.ck_gamma[j];
        }

        /* Force a call to setup */
        ida_mem.ida_forceSetup = SUNTRUE;
    }

    IDA_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Functions specific to cubic Hermite interpolation
 * -----------------------------------------------------------------
 */

/*
 * IDAAhermiteMalloc
 *
 * This routine allocates memory for storing information at all
 * intermediate points between two consecutive check points.
 * This data is then used to interpolate the forward solution
 * at any other time.  (Allocation is infallible; the C rollback/
 * SUNFALSE branches vanish and no boolean is returned.)
 */
fn IDAAhermiteMalloc(ida_mem: &mut IDAMem) {
    let ns = ida_mem.ida_Ns as usize;
    let tmpl = std::mem::take(&mut ida_mem.ida_tempv1);

    let idaadj_mem = adj(ida_mem);

    /* Allocate space for the vectors yyTmp and ypTmp. */
    idaadj_mem.ia_yyTmp = N_VClone(&tmpl);
    idaadj_mem.ia_ypTmp = N_VClone(&tmpl);

    /* Allocate space for sensitivities temporary vectors. */
    if idaadj_mem.ia_storeSensi {
        idaadj_mem.ia_yySTmp = (0..ns).map(|_| N_VClone(&tmpl)).collect();
        idaadj_mem.ia_ypSTmp = (0..ns).map(|_| N_VClone(&tmpl)).collect();
    }

    /* Allocate space for the content field of the dt structures */
    let store_sensi = idaadj_mem.ia_storeSensi;
    for i in 0..=(idaadj_mem.ia_nsteps as usize) {
        idaadj_mem.dt_mem[i].content = DtpntContent::Hermite {
            y: N_VClone(&tmpl),
            yd: N_VClone(&tmpl),
            yS: if store_sensi { (0..ns).map(|_| N_VClone(&tmpl)).collect() } else { Vec::new() },
            ySd: if store_sensi { (0..ns).map(|_| N_VClone(&tmpl)).collect() } else { Vec::new() },
        };
    }

    ida_mem.ida_tempv1 = tmpl;
}

/* (IDAAhermiteFree: RAII — replaced by dropping the DtpntContent and
   workspace vectors with the IDAadjMem.) */

/*
 * IDAAhermiteStorePnt
 *
 * This routine stores a new point (y,yd) in the structure d for use
 * in the cubic Hermite interpolation.
 * Note that the time is already stored.
 */
fn IDAAhermiteStorePnt(ida_mem: &mut IDAMem, idx: usize) -> i32 {
    /* Detach the adjoint memory so the dt point can be filled from
       ida_mem (borrow discipline) */
    let mut idaadj_mem = ida_mem.ida_adj_mem.take().unwrap();
    let store_sensi = idaadj_mem.ia_storeSensi;

    if let DtpntContent::Hermite { y, yd, yS, ySd } = &mut idaadj_mem.dt_mem[idx].content {
        /* Load solution(s) */
        N_VScale(ONE, &ida_mem.ida_phi[0], y);

        if store_sensi {
            /* (C: cvals[is]=ONE + N_VScaleVectorArray — plain copy) */
            for is in 0..ida_mem.ida_Ns as usize {
                N_VScale(ONE, &ida_mem.ida_phiS[0][is], &mut yS[is]);
            }
        }

        /* Load derivative(s). */
        IDAAGettnSolutionYp(ida_mem, yd);

        if store_sensi {
            IDAAGettnSolutionYpS(ida_mem, ySd);
        }
    }

    ida_mem.ida_adj_mem = Some(idaadj_mem);
    0
}

/*
 * IDAAhermiteGetY
 *
 * This routine uses cubic piece-wise Hermite interpolation for
 * the forward solution vector.
 * It is typically called by the wrapper routines before calling
 * user provided routines (fB, djacB, bjacB, jtimesB, psolB) but
 * can be directly called by the user through IDAGetAdjY
 */
fn IDAAhermiteGetY(
    ida_mem: &mut IDAMem,
    t: f64,
    yy: &mut NVector,
    yp: &mut NVector,
    yyS: &mut Vec<NVector>,
    ypS: &mut Vec<NVector>,
) -> i32 {
    /* Local value of Ns */
    let NS = if adj(ida_mem).ia_interpSensi && !yyS.is_empty() {
        ida_mem.ida_Ns as usize
    } else {
        0
    };

    /* Get the index in dt_mem */
    let mut index: i64 = 0;
    let mut newpoint = SUNFALSE;
    let flag = IDAAfindIndex(ida_mem, t, &mut index, &mut newpoint);
    if flag != IDA_SUCCESS {
        return flag;
    }

    /* Detach the adjoint memory: the interpolation reads dt_mem and
       writes the ia_Y/ia_YS scratch (disjoint fields, split borrow) */
    let mut idaadj_mem = ida_mem.ida_adj_mem.take().unwrap();
    let flag = {
        let IDAadjMem { dt_mem, ia_Y, ia_YS, .. } = &mut *idaadj_mem;

        idaa_hermite_get_y_inner(dt_mem, ia_Y, ia_YS, t, index, newpoint, NS, yy, yp, yyS, ypS)
    };
    ida_mem.ida_adj_mem = Some(idaadj_mem);
    flag
}

#[allow(clippy::too_many_arguments)]
fn idaa_hermite_get_y_inner(
    dt_mem: &[IDAdtpntMem],
    ia_Y: &mut [NVector],
    ia_YS: &mut [Vec<NVector>; MXORDP1],
    t: f64,
    index: i64,
    newpoint: bool,
    NS: usize,
    yy: &mut NVector,
    yp: &mut NVector,
    yyS: &mut Vec<NVector>,
    ypS: &mut Vec<NVector>,
) -> i32 {
    /* If we are beyond the left limit but close enough,
       then return y at the left limit. */

    if index == 0 {
        if let DtpntContent::Hermite { y, yd, yS, ySd } = &dt_mem[0].content {
            N_VScale(ONE, y, yy);
            N_VScale(ONE, yd, yp);

            if NS > 0 {
                /* (C: cvals[is]=ONE + N_VScaleVectorArray — plain copy) */
                for is in 0..NS {
                    N_VScale(ONE, &yS[is], &mut yyS[is]);
                }
                for is in 0..NS {
                    N_VScale(ONE, &ySd[is], &mut ypS[is]);
                }
            }
        }
        return IDA_SUCCESS;
    }

    /* Extract stuff from the appropriate data points */
    let iu = index as usize;
    let t0 = dt_mem[iu - 1].t;
    let t1 = dt_mem[iu].t;
    let delta = t1 - t0;

    let (y0, yd0, yS0, ySd0) = match &dt_mem[iu - 1].content {
        DtpntContent::Hermite { y, yd, yS, ySd } => (y, yd, yS, ySd),
        _ => unreachable!("hermite content"),
    };

    if newpoint {
        /* Recompute Y0 and Y1 */
        let (y1, yd1, yS1, ySd1) = match &dt_mem[iu].content {
            DtpntContent::Hermite { y, yd, yS, ySd } => (y, yd, yS, ySd),
            _ => unreachable!("hermite content"),
        };

        /* Y1 = delta (yd1 + yd0) - 2 (y1 - y0)
           (C: N_VLinearCombination(4, {-2, 2, delta, delta},
           {y1, y0, yd1, yd0}, Y[1]); serial kernel: z = c0*x0 then
           z += ci*xi) */
        for (z, x) in ia_Y[1].data.iter_mut().zip(&y1.data) {
            *z = -TWO * *x;
        }
        for (z, x) in ia_Y[1].data.iter_mut().zip(&y0.data) {
            *z += TWO * *x;
        }
        for (z, x) in ia_Y[1].data.iter_mut().zip(&yd1.data) {
            *z += delta * *x;
        }
        for (z, x) in ia_Y[1].data.iter_mut().zip(&yd0.data) {
            *z += delta * *x;
        }

        /* Y0 = y1 - y0 - delta * yd0 */
        for (z, x) in ia_Y[0].data.iter_mut().zip(&y1.data) {
            *z = *x;
        }
        for (z, x) in ia_Y[0].data.iter_mut().zip(&y0.data) {
            *z += -ONE * *x;
        }
        for (z, x) in ia_Y[0].data.iter_mut().zip(&yd0.data) {
            *z += -delta * *x;
        }

        /* Recompute YS0 and YS1, if needed */

        if NS > 0 {
            /* YS1 = delta (ySd1 + ySd0) - 2 (yS1 - yS0) */
            for is in 0..NS {
                for (z, x) in ia_YS[1][is].data.iter_mut().zip(&yS1[is].data) {
                    *z = -TWO * *x;
                }
                for (z, x) in ia_YS[1][is].data.iter_mut().zip(&yS0[is].data) {
                    *z += TWO * *x;
                }
                for (z, x) in ia_YS[1][is].data.iter_mut().zip(&ySd1[is].data) {
                    *z += delta * *x;
                }
                for (z, x) in ia_YS[1][is].data.iter_mut().zip(&ySd0[is].data) {
                    *z += delta * *x;
                }
            }

            /* YS0 = yS1 - yS0 - delta * ySd0 */
            for is in 0..NS {
                for (z, x) in ia_YS[0][is].data.iter_mut().zip(&yS1[is].data) {
                    *z = *x;
                }
                for (z, x) in ia_YS[0][is].data.iter_mut().zip(&yS0[is].data) {
                    *z += -ONE * *x;
                }
                for (z, x) in ia_YS[0][is].data.iter_mut().zip(&ySd0[is].data) {
                    *z += -delta * *x;
                }
            }
        }
    }

    /* Perform the actual interpolation. */

    /* For y. */
    let factor1 = t - t0;

    let mut factor2 = factor1 / delta;
    factor2 = factor2 * factor2;

    let factor3 = factor2 * (t - t1) / delta;

    /* y = y0 + factor1 yd0 + factor2 * Y[0] + factor3 Y[1] */
    for (z, x) in yy.data.iter_mut().zip(&y0.data) {
        *z = *x;
    }
    for (z, x) in yy.data.iter_mut().zip(&yd0.data) {
        *z += factor1 * *x;
    }
    for (z, x) in yy.data.iter_mut().zip(&ia_Y[0].data) {
        *z += factor2 * *x;
    }
    for (z, x) in yy.data.iter_mut().zip(&ia_Y[1].data) {
        *z += factor3 * *x;
    }

    /* Sensi Interpolation. */

    /* yS = yS0 + factor1 ySd0 + factor2 * YS[0] + factor3 YS[1], if needed */
    if NS > 0 {
        for is in 0..NS {
            for (z, x) in yyS[is].data.iter_mut().zip(&yS0[is].data) {
                *z = *x;
            }
            for (z, x) in yyS[is].data.iter_mut().zip(&ySd0[is].data) {
                *z += factor1 * *x;
            }
            for (z, x) in yyS[is].data.iter_mut().zip(&ia_YS[0][is].data) {
                *z += factor2 * *x;
            }
            for (z, x) in yyS[is].data.iter_mut().zip(&ia_YS[1][is].data) {
                *z += factor3 * *x;
            }
        }
    }

    /* For y'. */
    let mut factor1 = factor1 / delta / delta; /* factor1 = 2(t-t0)/(t1-t0)^2             */
    let factor2 = factor1 * ((3.0 * t - TWO * t1 - t0) /
                             delta); /* factor2 = (t-t0)(3*t-2*t1-t0)/(t1-t0)^3 */
    factor1 *= 2.0;

    /* yp = yd0 + factor1 Y[0] + factor 2 Y[1] */
    for (z, x) in yp.data.iter_mut().zip(&yd0.data) {
        *z = *x;
    }
    for (z, x) in yp.data.iter_mut().zip(&ia_Y[0].data) {
        *z += factor1 * *x;
    }
    for (z, x) in yp.data.iter_mut().zip(&ia_Y[1].data) {
        *z += factor2 * *x;
    }

    /* Sensi interpolation for 1st derivative. */

    /* ypS = ySd0 + factor1 YS[0] + factor 2 YS[1], if needed */
    if NS > 0 {
        for is in 0..NS {
            for (z, x) in ypS[is].data.iter_mut().zip(&ySd0[is].data) {
                *z = *x;
            }
            for (z, x) in ypS[is].data.iter_mut().zip(&ia_YS[0][is].data) {
                *z += factor1 * *x;
            }
            for (z, x) in ypS[is].data.iter_mut().zip(&ia_YS[1][is].data) {
                *z += factor2 * *x;
            }
        }
    }

    IDA_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * Functions specific to Polynomial interpolation
 * -----------------------------------------------------------------
 */

/*
 * IDAApolynomialMalloc
 *
 * This routine allocates memory for storing information at all
 * intermediate points between two consecutive check points.
 * This data is then used to interpolate the forward solution
 * at any other time.  (Allocation is infallible; the C rollback/
 * SUNFALSE branches vanish and no boolean is returned.)
 *
 * Information about the first derivative is stored only for the first
 * data point.
 */
fn IDAApolynomialMalloc(ida_mem: &mut IDAMem) {
    let ns = ida_mem.ida_Ns as usize;
    let tmpl = std::mem::take(&mut ida_mem.ida_tempv1);

    let idaadj_mem = adj(ida_mem);

    /* Allocate space for the vectors yyTmp and ypTmp */
    idaadj_mem.ia_yyTmp = N_VClone(&tmpl);
    idaadj_mem.ia_ypTmp = N_VClone(&tmpl);

    if idaadj_mem.ia_storeSensi {
        idaadj_mem.ia_yySTmp = (0..ns).map(|_| N_VClone(&tmpl)).collect();
        idaadj_mem.ia_ypSTmp = (0..ns).map(|_| N_VClone(&tmpl)).collect();
    }

    /* Allocate space for the content field of the dt structures:
       yd/ySd only for the first data point (i == 0). */
    let store_sensi = idaadj_mem.ia_storeSensi;
    for i in 0..=(idaadj_mem.ia_nsteps as usize) {
        idaadj_mem.dt_mem[i].content = DtpntContent::Polynomial {
            y: N_VClone(&tmpl),
            yd: if i == 0 { Some(N_VClone(&tmpl)) } else { None },
            yS: if store_sensi { (0..ns).map(|_| N_VClone(&tmpl)).collect() } else { Vec::new() },
            ySd: if store_sensi && i == 0 {
                (0..ns).map(|_| N_VClone(&tmpl)).collect()
            } else {
                Vec::new()
            },
            order: 0,
        };
    }

    ida_mem.ida_tempv1 = tmpl;
}

/* (IDAApolynomialFree: RAII, as IDAAhermiteFree.) */

/*
 * IDAApolynomialStorePnt
 *
 * This routine stores a new point y in the structure d for use
 * in the Polynomial interpolation.
 *
 * Note that the time is already stored. Information about the
 * first derivative is available only for the first data point,
 * in which case content->yp is non-null.
 */
fn IDAApolynomialStorePnt(ida_mem: &mut IDAMem, idx: usize) -> i32 {
    let mut idaadj_mem = ida_mem.ida_adj_mem.take().unwrap();
    let store_sensi = idaadj_mem.ia_storeSensi;

    if let DtpntContent::Polynomial { y, yd, yS, ySd, order } =
        &mut idaadj_mem.dt_mem[idx].content
    {
        N_VScale(ONE, &ida_mem.ida_phi[0], y);

        /* copy also the derivative for the first data point (in this case
           content->yp is non-null). */
        if let Some(yd) = yd {
            IDAAGettnSolutionYp(ida_mem, yd);
        }

        if store_sensi {
            /* (C: cvals[is]=ONE + N_VScaleVectorArray — plain copy) */
            for is in 0..ida_mem.ida_Ns as usize {
                N_VScale(ONE, &ida_mem.ida_phiS[0][is], &mut yS[is]);
            }

            /* store the derivative if it is the first data point. */
            if !ySd.is_empty() {
                IDAAGettnSolutionYpS(ida_mem, ySd);
            }
        }

        *order = ida_mem.ida_kused;
    }

    ida_mem.ida_adj_mem = Some(idaadj_mem);
    0
}

/*
 * IDAApolynomialGetY
 *
 * This routine uses polynomial interpolation for the forward solution
 * vector.  It is typically called by the wrapper routines before
 * calling user provided routines (fB, djacB, bjacB, jtimesB, psolB))
 * but can be directly called by the user through CVodeGetAdjY.
 */
fn IDAApolynomialGetY(
    ida_mem: &mut IDAMem,
    t: f64,
    yy: &mut NVector,
    yp: &mut NVector,
    yyS: &mut Vec<NVector>,
    ypS: &mut Vec<NVector>,
) -> i32 {
    /* Local value of Ns */
    let NS = if adj(ida_mem).ia_interpSensi && !yyS.is_empty() {
        ida_mem.ida_Ns as usize
    } else {
        0
    };

    /* Get the index in dt_mem */
    let mut index: i64 = 0;
    let mut newpoint = SUNFALSE;
    let flag = IDAAfindIndex(ida_mem, t, &mut index, &mut newpoint);
    if flag != IDA_SUCCESS {
        return flag;
    }

    let (tinitial, tfinal, ia_np) = {
        let idaadj_mem = ida_mem.ida_adj_mem.as_ref().unwrap();
        (idaadj_mem.ia_tinitial, idaadj_mem.ia_tfinal, idaadj_mem.ia_np)
    };

    let mut idaadj_mem = ida_mem.ida_adj_mem.take().unwrap();
    let flag = {
        let IDAadjMem { dt_mem, ia_Y, ia_YS, ia_T, .. } = &mut *idaadj_mem;
        idaa_polynomial_get_y_inner(dt_mem, ia_Y, ia_YS, ia_T, t, index, newpoint, NS, tinitial,
                                    tfinal, ia_np, yy, yp, yyS, ypS)
    };
    ida_mem.ida_adj_mem = Some(idaadj_mem);
    flag
}

#[allow(clippy::too_many_arguments)]
fn idaa_polynomial_get_y_inner(
    dt_mem: &[IDAdtpntMem],
    ia_Y: &mut [NVector],
    ia_YS: &mut [Vec<NVector>; MXORDP1],
    ia_T: &mut [f64; MXORDP1],
    t: f64,
    index: i64,
    newpoint: bool,
    NS: usize,
    tinitial: f64,
    tfinal: f64,
    ia_np: i64,
    yy: &mut NVector,
    yp: &mut NVector,
    yyS: &mut Vec<NVector>,
    ypS: &mut Vec<NVector>,
) -> i32 {
    /* If we are beyond the left limit but close enough,
       then return y at the left limit. */

    if index == 0 {
        if let DtpntContent::Polynomial { y, yd, yS, ySd, .. } = &dt_mem[0].content {
            N_VScale(ONE, y, yy);
            N_VScale(ONE, yd.as_ref().unwrap(), yp);

            if NS > 0 {
                for is in 0..NS {
                    N_VScale(ONE, &yS[is], &mut yyS[is]);
                }
                for is in 0..NS {
                    N_VScale(ONE, &ySd[is], &mut ypS[is]);
                }
            }
        }
        return IDA_SUCCESS;
    }

    /* Scaling factor */
    let iu = index as usize;
    let delt = SUNRabs(dt_mem[iu].t - dt_mem[iu - 1].t);

    /* Find the direction of the forward integration */
    let dir: i32 = if tfinal - tinitial > ZERO { 1 } else { -1 };

    /* Establish the base point depending on the integration direction.
       Modify the base if there are not enough points for the current
       order */

    let poly_order = |i: usize| -> i32 {
        match &dt_mem[i].content {
            DtpntContent::Polynomial { order, .. } => *order,
            _ => unreachable!("polynomial content"),
        }
    };

    let mut base: i64;
    let order: i32;
    if dir == 1 {
        base = index;
        order = poly_order(base as usize);
        if index < order as i64 {
            base += order as i64 - index;
        }
    } else {
        base = index - 1;
        order = poly_order(base as usize);
        if ia_np - index > order as i64 {
            base -= index + order as i64 - ia_np;
        }
    }
    let order_u = order as usize;

    /* Recompute Y (divided differences for Newton polynomial) if needed */

    if newpoint {
        /* Store 0-th order DD */
        if dir == 1 {
            for j in 0..=order_u {
                let di = (base - j as i64) as usize;
                ia_T[j] = dt_mem[di].t;
                if let DtpntContent::Polynomial { y, yS, .. } = &dt_mem[di].content {
                    N_VScale(ONE, y, &mut ia_Y[j]);

                    if NS > 0 {
                        for is in 0..NS {
                            N_VScale(ONE, &yS[is], &mut ia_YS[j][is]);
                        }
                    }
                }
            }
        } else {
            for j in 0..=order_u {
                let di = (base - 1 + j as i64) as usize;
                ia_T[j] = dt_mem[di].t;
                if let DtpntContent::Polynomial { y, yS, .. } = &dt_mem[di].content {
                    N_VScale(ONE, y, &mut ia_Y[j]);

                    if NS > 0 {
                        for is in 0..NS {
                            N_VScale(ONE, &yS[is], &mut ia_YS[j][is]);
                        }
                    }
                }
            }
        }

        /* Compute higher-order DD */
        for i in 1..=order_u {
            for j in (i..=order_u).rev() {
                let factor = delt / (ia_T[j] - ia_T[j - i]);
                /* Y[j] = factor*Y[j] - factor*Y[j-1] (aliased in-place) */
                {
                    let (lo, hi) = ia_Y.split_at_mut(j);
                    hi[0].linear_sum_with(factor, -factor, &lo[j - 1]);
                }

                for is in 0..NS {
                    let (lo, hi) = ia_YS.split_at_mut(j);
                    hi[0][is].linear_sum_with(factor, -factor, &lo[j - 1][is]);
                }
            }
        }
    }

    /* Perform the actual interpolation for yy using nested
       multiplications */

    let mut cvals = [ZERO; MXORDP1];
    cvals[0] = ONE;
    for i in 0..order_u {
        cvals[i + 1] = cvals[i] * (t - ia_T[i]) / delt;
    }

    /* yy = N_VLinearCombination(order + 1, cvals, Y, yy) */
    for (z, x) in yy.data.iter_mut().zip(&ia_Y[0].data) {
        *z = cvals[0] * *x;
    }
    for j in 1..=order_u {
        for (z, x) in yy.data.iter_mut().zip(&ia_Y[j].data) {
            *z += cvals[j] * *x;
        }
    }

    if NS > 0 {
        for is in 0..NS {
            for (z, x) in yyS[is].data.iter_mut().zip(&ia_YS[0][is].data) {
                *z = cvals[0] * *x;
            }
            for j in 1..=order_u {
                for (z, x) in yyS[is].data.iter_mut().zip(&ia_YS[j][is].data) {
                    *z += cvals[j] * *x;
                }
            }
        }
    }

    /* Perform the actual interpolation for yp.

       Writing p(t) = y0 + (t-t0)*f[t0,t1] + ... + (t-t0)(t-t1)...(t-tn)*f[t0,t1,...tn],
       denote psi_k(t) = (t-t0)(t-t1)...(t-tk).

       The formula used for p'(t) is:
         - p'(t) = f[t0,t1] + psi_1'(t)*f[t0,t1,t2] + ... + psi_n'(t)*f[t0,t1,...,tn]

       We recursively compute psi_k'(t) from:
         - psi_k'(t) = (t-tk)*psi_{k-1}'(t) + psi_{k-1}

       psi_k is rescaled with 1/delt each time is computed, because the
       Newton DDs from Y were scaled with delt.
    */

    let mut Psi = ONE;
    let mut Psiprime = ZERO;

    for i in 1..=order_u {
        let factor = (t - ia_T[i - 1]) / delt;

        Psiprime = Psi / delt + factor * Psiprime;
        Psi = Psi * factor;

        cvals[i - 1] = Psiprime;
    }

    /* yp = N_VLinearCombination(order, cvals, Y + 1, yp) */
    for (z, x) in yp.data.iter_mut().zip(&ia_Y[1].data) {
        *z = cvals[0] * *x;
    }
    for j in 2..=order_u {
        for (z, x) in yp.data.iter_mut().zip(&ia_Y[j].data) {
            *z += cvals[j - 1] * *x;
        }
    }

    if NS > 0 {
        for is in 0..NS {
            for (z, x) in ypS[is].data.iter_mut().zip(&ia_YS[1][is].data) {
                *z = cvals[0] * *x;
            }
            for j in 2..=order_u {
                for (z, x) in ypS[is].data.iter_mut().zip(&ia_YS[j][is].data) {
                    *z += cvals[j - 1] * *x;
                }
            }
        }
    }

    IDA_SUCCESS
}

/*
 * IDAAGettnSolutionYp
 *
 * Evaluates the first derivative of the solution at the last time
 * returned by IDASolve (tretlast).
 *
 * The function implements the same algorithm as in IDAGetSolution but
 * in the particular case when t=tn (i.e. delta=0).
 *
 * This function was implemented to avoid calls to IDAGetSolution which
 * computes y by doing a loop that is not necessary for this particular
 * situation.
 */
fn IDAAGettnSolutionYp(ida_mem: &IDAMem, yp: &mut NVector) -> i32 {
    if ida_mem.ida_nst == 0 {
        /* If no integration was done, return the yp supplied by user.*/
        N_VScale(ONE, &ida_mem.ida_phi[1], yp);

        return 0;
    }

    /* Compute yp as in IDAGetSolution for this particular case when t=tn. */

    let mut kord = ida_mem.ida_kused;
    if ida_mem.ida_kused == 0 {
        kord = 1;
    }

    let mut c = ONE;
    let mut d = ZERO;
    let mut gam = ZERO;
    let mut dvals = [ZERO; MXORDP1];
    for j in 1..=(kord as usize) {
        d = d * gam + c / ida_mem.ida_psi[j - 1];
        c = c * gam;
        gam = ida_mem.ida_psi[j - 1] / ida_mem.ida_psi[j];

        dvals[j - 1] = d;
    }

    /* retval = N_VLinearCombination(kord, dvals, phi + 1, yp) */
    let kord_u = kord as usize;
    for (z, x) in yp.data.iter_mut().zip(&ida_mem.ida_phi[1].data) {
        *z = dvals[0] * *x;
    }
    for j in 2..=kord_u {
        for (z, x) in yp.data.iter_mut().zip(&ida_mem.ida_phi[j].data) {
            *z += dvals[j - 1] * *x;
        }
    }

    0
}

/*
 * IDAAGettnSolutionYpS
 *
 * Same as IDAAGettnSolutionYp, but for first derivative of the
 * sensitivities.
 */
fn IDAAGettnSolutionYpS(ida_mem: &IDAMem, ypS: &mut [NVector]) -> i32 {
    if ida_mem.ida_nst == 0 {
        /* If no integration was done, return the ypS supplied by user.*/
        /* (C: cvals[is]=ONE + N_VScaleVectorArray — plain copy) */
        for is in 0..ida_mem.ida_Ns as usize {
            N_VScale(ONE, &ida_mem.ida_phiS[1][is], &mut ypS[is]);
        }

        return 0;
    }

    let mut kord = ida_mem.ida_kused;
    if ida_mem.ida_kused == 0 {
        kord = 1;
    }

    let mut c = ONE;
    let mut d = ZERO;
    let mut gam = ZERO;
    let mut dvals = [ZERO; MXORDP1];
    for j in 1..=(kord as usize) {
        d = d * gam + c / ida_mem.ida_psi[j - 1];
        c = c * gam;
        gam = ida_mem.ida_psi[j - 1] / ida_mem.ida_psi[j];

        dvals[j - 1] = d;
    }

    /* retval = N_VLinearCombinationVectorArray(Ns, kord, dvals,
       phiS + 1, ypS) */
    let kord_u = kord as usize;
    for is in 0..ida_mem.ida_Ns as usize {
        for (z, x) in ypS[is].data.iter_mut().zip(&ida_mem.ida_phiS[1][is].data) {
            *z = dvals[0] * *x;
        }
        for j in 2..=kord_u {
            for (z, x) in ypS[is].data.iter_mut().zip(&ida_mem.ida_phiS[j][is].data) {
                *z += dvals[j - 1] * *x;
            }
        }
    }

    0
}

/*
 * IDAAfindIndex
 *
 * Finds the index in the array of data point structures such that
 *     dt_mem[index-1].t <= t < dt_mem[index].t
 * If index is changed from the previous invocation, then newpoint =
 * SUNTRUE
 *
 * If t is beyond the leftmost limit, but close enough, index=0.
 *
 * Returns IDA_SUCCESS if successful and IDA_GETY_BADT if unable to
 * find index (t is too far beyond limits).
 */
fn IDAAfindIndex(
    ida_mem: &mut IDAMem,
    t: f64,
    index: &mut i64,
    newpoint: &mut bool,
) -> i32 {
    let uround = ida_mem.ida_uround;
    let idaadj_mem = adj(ida_mem);

    *newpoint = SUNFALSE;

    /* Find the direction of integration */
    let sign: f64 = if idaadj_mem.ia_tfinal - idaadj_mem.ia_tinitial > ZERO { 1.0 } else { -1.0 };

    /* If this is the first time we use new data */
    if idaadj_mem.ia_newData {
        idaadj_mem.ia_ilast = idaadj_mem.ia_np - 1;
        *newpoint = SUNTRUE;
        idaadj_mem.ia_newData = SUNFALSE;
    }

    /* Search for index starting from ilast */
    let dt_mem = &idaadj_mem.dt_mem;
    let ilast = idaadj_mem.ia_ilast;
    let to_left = sign * (t - dt_mem[(ilast - 1) as usize].t) < ZERO;
    let to_right = sign * (t - dt_mem[ilast as usize].t) > ZERO;

    if to_left {
        /* look for a new index to the left */

        *newpoint = SUNTRUE;

        *index = ilast;
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
            idaadj_mem.ia_ilast = 1;
        } else {
            idaadj_mem.ia_ilast = *index;
        }

        if *index == 0 {
            /* t is beyond leftmost limit. Is it too far? */
            if SUNRabs(t - idaadj_mem.dt_mem[0].t) > FUZZ_FACTOR * uround {
                return IDA_GETY_BADT;
            }
        }
    } else if to_right {
        /* look for a new index to the right */

        *newpoint = SUNTRUE;

        *index = ilast;
        loop {
            if sign * (t - dt_mem[*index as usize].t) > ZERO {
                *index += 1;
            } else {
                break;
            }
        }

        idaadj_mem.ia_ilast = *index;
    } else {
        /* ilast is still OK */

        *index = ilast;
    }
    IDA_SUCCESS
}

/*
 * IDAGetAdjY
 *
 * This routine returns the interpolated forward solution at time t.
 * The user must allocate space for y.
 */
pub fn IDAGetAdjY(ida_mem: &mut IDAMem, t: f64, yy: &mut NVector, yp: &mut NVector) -> i32 {
    /* (C dereferences ida_adj_mem unconditionally after the NULL
       ida_mem check; here the unwrap in IDAAgetY plays that role) */
    let mut noS: Vec<NVector> = Vec::new();
    let mut noSp: Vec<NVector> = Vec::new();
    IDAAgetY(ida_mem, t, yy, yp, &mut noS, &mut noSp)
}

/* Dispatch on ia_interpType (replaces the C ia_storePnt fn pointer). */
fn idaa_storePnt_dispatch(ida_mem: &mut IDAMem, idx: usize) -> i32 {
    if ida_mem.ida_adj_mem.as_ref().unwrap().ia_interpType == IDA_HERMITE {
        IDAAhermiteStorePnt(ida_mem, idx)
    } else {
        IDAApolynomialStorePnt(ida_mem, idx)
    }
}

/* Dispatch on ia_interpType (replaces the C ia_getY fn pointer).
   Empty yyS/ypS Vecs play the role of the C NULL arguments.  This is
   also the body of the idas_ls.rs idaLsGetY bridge (see the PIN
   recorded there). */
pub(crate) fn IDAAgetY(
    ida_mem: &mut IDAMem,
    t: f64,
    yy: &mut NVector,
    yp: &mut NVector,
    yyS: &mut Vec<NVector>,
    ypS: &mut Vec<NVector>,
) -> i32 {
    if ida_mem.ida_adj_mem.as_ref().unwrap().ia_interpType == IDA_HERMITE {
        IDAAhermiteGetY(ida_mem, t, yy, yp, yyS, ypS)
    } else {
        IDAApolynomialGetY(ida_mem, t, yy, yp, yyS, ypS)
    }
}

/*=================================================================*/
/*             Wrappers for adjoint system                         */
/*=================================================================*/

/*
 * IDAAres
 *
 * This routine interfaces to the RhsFnB routine provided by
 * the user.  (Forward-callback-typed: the backward problem's
 * UserData holds the OUTER IDAMem, installed transiently by
 * IDASolveB/IDACalcICB(S) — the pinned ownership dance.)
 */
pub(crate) fn IDAAres(
    tt: f64,
    yyB: &NVector,
    ypB: &NVector,
    rrB: &mut NVector,
    ida_mem: &mut UserData,
) -> i32 {
    /* (C: IDA_mem = (IDAMem) ida_mem) */
    let ida_mem = ida_mem.as_mut().unwrap().downcast_mut::<IDAMem>().unwrap();

    /* Get the current backward problem. */
    let which = ida_mem.ida_adj_mem.as_ref().unwrap().ia_bckpbCrt.unwrap();

    /* Get forward solution from interpolation. */
    if !adj(ida_mem).ia_noInterp {
        let interp_sensi = adj(ida_mem).ia_interpSensi;
        let (mut yyTmp, mut ypTmp, mut yySTmp, mut ypSTmp) = {
            let idaadj_mem = adj(ida_mem);
            (std::mem::take(&mut idaadj_mem.ia_yyTmp),
             std::mem::take(&mut idaadj_mem.ia_ypTmp),
             std::mem::take(&mut idaadj_mem.ia_yySTmp),
             std::mem::take(&mut idaadj_mem.ia_ypSTmp))
        };
        let flag = if interp_sensi {
            IDAAgetY(ida_mem, tt, &mut yyTmp, &mut ypTmp, &mut yySTmp, &mut ypSTmp)
        } else {
            let mut noS: Vec<NVector> = Vec::new();
            let mut noSp: Vec<NVector> = Vec::new();
            IDAAgetY(ida_mem, tt, &mut yyTmp, &mut ypTmp, &mut noS, &mut noSp)
        };
        {
            let idaadj_mem = adj(ida_mem);
            idaadj_mem.ia_yyTmp = yyTmp;
            idaadj_mem.ia_ypTmp = ypTmp;
            idaadj_mem.ia_yySTmp = yySTmp;
            idaadj_mem.ia_ypSTmp = ypSTmp;
        }

        if flag != IDA_SUCCESS {
            IDAProcessError(Some(ida_mem), -1, line!(), "IDAAres", file!(),
                            &ida_msg_g(MSGAM_BAD_TINTERP, &[tt]));
            return -1;
        }
    }

    /* Call the user supplied residual. */
    let (yyTmp, ypTmp, yySTmp, ypSTmp) = {
        let idaadj_mem = adj(ida_mem);
        (std::mem::take(&mut idaadj_mem.ia_yyTmp),
         std::mem::take(&mut idaadj_mem.ia_ypTmp),
         std::mem::take(&mut idaadj_mem.ia_yySTmp),
         std::mem::take(&mut idaadj_mem.ia_ypSTmp))
    };
    let idaB_mem = &mut adj(ida_mem).IDAB_mem[which];
    let retval = if idaB_mem.ida_res_withSensi {
        let resS = idaB_mem.ida_resS.unwrap();
        resS(tt, &yyTmp, &ypTmp, &yySTmp, &ypSTmp, yyB, ypB, rrB, &mut idaB_mem.ida_user_data)
    } else {
        let res = idaB_mem.ida_res.unwrap();
        res(tt, &yyTmp, &ypTmp, yyB, ypB, rrB, &mut idaB_mem.ida_user_data)
    };
    {
        let idaadj_mem = adj(ida_mem);
        idaadj_mem.ia_yyTmp = yyTmp;
        idaadj_mem.ia_ypTmp = ypTmp;
        idaadj_mem.ia_yySTmp = yySTmp;
        idaadj_mem.ia_ypSTmp = ypSTmp;
    }
    retval
}

/*
 * IDAArhsQ
 *
 * This routine interfaces to the IDAQuadRhsFnB routine provided by
 * the user.
 *
 * It is passed to IDAQuadInit calls for backward problem, so it must
 * be of IDAQuadRhsFn type.
 */
pub(crate) fn IDAArhsQ(
    tt: f64,
    yyB: &NVector,
    ypB: &NVector,
    resvalQB: &mut NVector,
    ida_mem: &mut UserData,
) -> i32 {
    /* (C: IDA_mem = (IDAMem) ida_mem) */
    let ida_mem = ida_mem.as_mut().unwrap().downcast_mut::<IDAMem>().unwrap();

    /* Get current backward problem. */
    let which = ida_mem.ida_adj_mem.as_ref().unwrap().ia_bckpbCrt.unwrap();

    /* Get forward solution from interpolation. */
    if !adj(ida_mem).ia_noInterp {
        let interp_sensi = adj(ida_mem).ia_interpSensi;
        let (mut yyTmp, mut ypTmp, mut yySTmp, mut ypSTmp) = {
            let idaadj_mem = adj(ida_mem);
            (std::mem::take(&mut idaadj_mem.ia_yyTmp),
             std::mem::take(&mut idaadj_mem.ia_ypTmp),
             std::mem::take(&mut idaadj_mem.ia_yySTmp),
             std::mem::take(&mut idaadj_mem.ia_ypSTmp))
        };
        let flag = if interp_sensi {
            IDAAgetY(ida_mem, tt, &mut yyTmp, &mut ypTmp, &mut yySTmp, &mut ypSTmp)
        } else {
            let mut noS: Vec<NVector> = Vec::new();
            let mut noSp: Vec<NVector> = Vec::new();
            IDAAgetY(ida_mem, tt, &mut yyTmp, &mut ypTmp, &mut noS, &mut noSp)
        };
        {
            let idaadj_mem = adj(ida_mem);
            idaadj_mem.ia_yyTmp = yyTmp;
            idaadj_mem.ia_ypTmp = ypTmp;
            idaadj_mem.ia_yySTmp = yySTmp;
            idaadj_mem.ia_ypSTmp = ypSTmp;
        }

        if flag != IDA_SUCCESS {
            IDAProcessError(Some(ida_mem), -1, line!(), "IDAArhsQ", file!(),
                            &ida_msg_g(MSGAM_BAD_TINTERP, &[tt]));
            return -1;
        }
    }

    /* Call user's adjoint quadrature RHS routine */
    let (yyTmp, ypTmp, yySTmp, ypSTmp) = {
        let idaadj_mem = adj(ida_mem);
        (std::mem::take(&mut idaadj_mem.ia_yyTmp),
         std::mem::take(&mut idaadj_mem.ia_ypTmp),
         std::mem::take(&mut idaadj_mem.ia_yySTmp),
         std::mem::take(&mut idaadj_mem.ia_ypSTmp))
    };
    let idaB_mem = &mut adj(ida_mem).IDAB_mem[which];
    let retval = if idaB_mem.ida_rhsQ_withSensi {
        let rhsQS = idaB_mem.ida_rhsQS.unwrap();
        rhsQS(tt, &yyTmp, &ypTmp, &yySTmp, &ypSTmp, yyB, ypB, resvalQB,
              &mut idaB_mem.ida_user_data)
    } else {
        let rhsQ = idaB_mem.ida_rhsQ.unwrap();
        rhsQ(tt, &yyTmp, &ypTmp, yyB, ypB, resvalQB, &mut idaB_mem.ida_user_data)
    };
    {
        let idaadj_mem = adj(ida_mem);
        idaadj_mem.ia_yyTmp = yyTmp;
        idaadj_mem.ia_ypTmp = ypTmp;
        idaadj_mem.ia_yySTmp = yySTmp;
        idaadj_mem.ia_ypSTmp = ypSTmp;
    }
    retval
}

/*===============================================================
  Tests
  ===============================================================*/
#[cfg(test)]
mod tests {
    use super::*;

    fn make_ida_mem(n: usize) -> IDAMem {
        let mut ida_mem = IDAMem::default();
        ida_mem.ida_tempv1 = NVector::new(n);
        ida_mem
    }

    /* IDAAdjInit validates steps/interp and initializes the adjoint
       block (idaa.c lines 106-228); IDAAdjReInit requires prior init
       and clears the check point list; IDAAdjFree drops everything. */
    #[test]
    fn idaadjinit_validation_and_state() {
        let mut ida_mem = make_ida_mem(2);

        assert_eq!(IDAAdjInit(&mut ida_mem, 0, IDA_HERMITE), IDA_ILL_INPUT);
        assert_eq!(IDAAdjInit(&mut ida_mem, -5, IDA_POLYNOMIAL), IDA_ILL_INPUT);
        assert_eq!(IDAAdjInit(&mut ida_mem, 10, 0), IDA_ILL_INPUT);
        assert_eq!(IDAAdjReInit(&mut ida_mem), IDA_NO_ADJ);

        assert_eq!(IDAAdjInit(&mut ida_mem, 10, IDA_HERMITE), IDA_SUCCESS);
        assert!(ida_mem.ida_adj);
        assert!(ida_mem.ida_adjMallocDone);
        {
            let idaadj_mem = ida_mem.ida_adj_mem.as_ref().unwrap();
            assert_eq!(idaadj_mem.ia_nsteps, 10);
            assert_eq!(idaadj_mem.dt_mem.len(), 11);
            assert_eq!(idaadj_mem.ia_interpType, IDA_HERMITE);
            assert!(idaadj_mem.ia_storeSensi);
            assert!(!idaadj_mem.ia_interpSensi);
            assert!(!idaadj_mem.ia_noInterp);
            assert!(idaadj_mem.ia_firstIDAFcall);
            assert_eq!(idaadj_mem.ia_ilast, -1);
        }

        assert_eq!(IDAAdjReInit(&mut ida_mem), IDA_SUCCESS);

        IDAAdjFree(&mut ida_mem);
        assert!(ida_mem.ida_adj_mem.is_none());
        assert!(!ida_mem.ida_adjMallocDone);
    }

    /* IDACreateB assigns increasing indices and stores the backward
       problems in creation order; the ***B entry points validate
       `which` and the adjoint-initialized state. */
    #[test]
    fn idacreateb_indices_and_which_checks() {
        let mut ida_mem = make_ida_mem(2);

        let mut which = -1;
        assert_eq!(IDACreateB(&mut ida_mem, &mut which), IDA_NO_ADJ);

        assert_eq!(IDAAdjInit(&mut ida_mem, 5, IDA_POLYNOMIAL), IDA_SUCCESS);

        assert_eq!(IDACreateB(&mut ida_mem, &mut which), IDA_SUCCESS);
        assert_eq!(which, 0);
        assert_eq!(IDACreateB(&mut ida_mem, &mut which), IDA_SUCCESS);
        assert_eq!(which, 1);
        {
            let idaadj_mem = ida_mem.ida_adj_mem.as_ref().unwrap();
            assert_eq!(idaadj_mem.ia_nbckpbs, 2);
            assert_eq!(idaadj_mem.IDAB_mem[0].ida_index, 0);
            assert_eq!(idaadj_mem.IDAB_mem[1].ida_index, 1);
        }

        /* which out of range */
        assert_eq!(IDASStolerancesB(&mut ida_mem, 2, 1.0e-4, 1.0e-6), IDA_ILL_INPUT);

        /* before IDAInitB the nested solver rejects tolerances (as C:
           IDASStolerances requires MallocDone) */
        assert_eq!(IDASStolerancesB(&mut ida_mem, 0, 1.0e-4, 1.0e-6), IDA_NO_MALLOC);

        /* tolerances propagate to the nested solver after IDAInitB */
        fn resB(_tt: f64, _yy: &NVector, _yp: &NVector, _yyB: &NVector, _ypB: &NVector,
                _rrB: &mut NVector, _ud: &mut UserData) -> i32 {
            0
        }
        let yyB0 = NVector::from_slice(&[1.0, 2.0]);
        let ypB0 = NVector::from_slice(&[0.0, 0.0]);
        assert_eq!(IDAInitB(&mut ida_mem, 0, resB, 0.0, &yyB0, &ypB0), IDA_SUCCESS);
        assert_eq!(IDASStolerancesB(&mut ida_mem, 0, 1.0e-4, 1.0e-6), IDA_SUCCESS);
        {
            let idaadj_mem = ida_mem.ida_adj_mem.as_ref().unwrap();
            assert_eq!(idaadj_mem.IDAB_mem[0].IDA_mem.ida_rtol, 1.0e-4);
        }

        /* IDASolveB without IDASolveF: no backward integration possible */
        assert_eq!(IDASolveB(&mut ida_mem, 0.0, IDA_NORMAL), IDA_NO_FWD);
    }

    /* IDAInitB rejects a tB0 outside [tinitial, tfinal] (both 0 before
       IDASolveF ran) and otherwise initializes the nested solver plus
       the yy/yp workspaces. */
    #[test]
    fn idainitb_time_check_and_workspaces() {
        fn resB(_tt: f64, _yy: &NVector, _yp: &NVector, _yyB: &NVector, _ypB: &NVector,
                _rrB: &mut NVector, _ud: &mut UserData) -> i32 {
            0
        }

        let mut ida_mem = make_ida_mem(2);
        assert_eq!(IDAAdjInit(&mut ida_mem, 5, IDA_HERMITE), IDA_SUCCESS);
        let mut which = -1;
        assert_eq!(IDACreateB(&mut ida_mem, &mut which), IDA_SUCCESS);

        let yyB0 = NVector::from_slice(&[1.0, 2.0]);
        let ypB0 = NVector::from_slice(&[0.5, -0.5]);

        /* tB0 outside [ia_tinitial, ia_tfinal] = [0, 0] */
        assert_eq!(IDAInitB(&mut ida_mem, which, resB, 1.0, &yyB0, &ypB0), IDA_BAD_TB0);

        assert_eq!(IDAInitB(&mut ida_mem, which, resB, 0.0, &yyB0, &ypB0), IDA_SUCCESS);
        {
            let b = &ida_mem.ida_adj_mem.as_ref().unwrap().IDAB_mem[0];
            assert!(b.ida_res.is_some());
            assert!(!b.ida_res_withSensi);
            assert_eq!(b.ida_t0, 0.0);
            assert_eq!(b.ida_yy.data, vec![1.0, 2.0]);
            assert_eq!(b.ida_yp.data, vec![0.5, -0.5]);
            assert!(b.IDA_mem.ida_MallocDone);
        }
    }
}
