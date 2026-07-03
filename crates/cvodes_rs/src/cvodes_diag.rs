/* -----------------------------------------------------------------
 * Translated from src/cvodes/cvodes_diag.c (CVODES 7.7.0).
 * CVDIAG diagonal linear solver module: approximates the Newton
 * matrix M = I - gamma*J by a diagonal matrix built from difference
 * quotients, and applies its inverse.
 *
 * In C the module installs cv_linit/cv_lsetup/cv_lsolve/cv_lfree
 * function pointers plus a void* cv_lmem; here CVDiag attaches
 * LsModule::Diag(Box<CVDiagMem>) and cvodes.rs dispatches to
 * cvDiagInit/cvDiagSetup/cvDiagSolve (the CVDiagMem is detached
 * from CVodeMem for the duration of each call).
 * -----------------------------------------------------------------*/
use crate::cvodes_diag_impl::*;
use crate::cvodes_impl::{cvProcessError, CVodeMem, LsModule};
use crate::nvector_serial::{
    N_VAddConst, N_VClone, N_VCompare, N_VLinearSum, N_VProd, NVector,
};
use crate::sundials_types::{SUNFALSE, SUNTRUE};

/* Other Constants */

const FRACT: f64 = 0.1;
const ONE: f64 = 1.0;

/* Error messages (cvodes_diag_impl.h) */
const MSGDG_LMEM_NULL: &str = "CVDIAG memory is NULL.";
const MSGDG_RHSFUNC_FAILED: &str =
    "The right-hand side routine failed in an unrecoverable manner.";
const MSGDG_NO_ADJ: &str = "Illegal attempt to call before calling CVodeAdjMalloc.";
const MSGDG_BAD_WHICH: &str = "Illegal value for which.";

/* N_VInvTest(v, v): the C call sites alias x with z, which the free
   function form cannot express under Rust borrows. Semantics match
   nvector_serial::N_VInvTest with x == z (zero entries untouched). */
fn N_VInvTest_inplace(z: &mut NVector) -> bool {
    let mut no_zero_found = true;
    for zi in &mut z.data {
        if *zi == 0.0 {
            no_zero_found = false;
        } else {
            *zi = ONE / *zi;
        }
    }
    no_zero_found
}

/*
 * ================================================================
 *
 *                   PART I - forward problems
 *
 * ================================================================
 */

/*
 * -----------------------------------------------------------------
 * CVDiag
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the diagonal linear solver module.  It allocates
 * memory for a structure of type CVDiagMemRec and sets the cv_lmem
 * field to this structure.  Finally, it allocates memory for M, bit,
 * and bitcomp.
 * -----------------------------------------------------------------
 */
pub fn CVDiag(cv_mem: &mut CVodeMem) -> i32 {
    /* Return immediately if cvode_mem is NULL: impossible under &mut. */

    /* Check if N_VCompare and N_VInvTest are present: always true for
       the serial NVector (covered by the type system). */

    /* if (cv_mem->cv_lfree != NULL) cv_mem->cv_lfree(cv_mem):
       overwriting cv_lmem below drops any previous module (RAII). */

    /* Set four main function fields in cv_mem: handled by the
       LsModule::Diag dispatch in cvodes.rs. */

    /* Get memory for CVDiagMemRec; allocate memory for M, bit, and
       bitcomp (allocation failures abort in Rust). */
    let cvdiag_mem = Box::new(CVDiagMem {
        di_gammasv: 0.0, /* set by cvDiagSetup before first use */
        di_M: N_VClone(&cv_mem.cv_tempv),
        di_bit: N_VClone(&cv_mem.cv_tempv),
        di_bitcomp: N_VClone(&cv_mem.cv_tempv),
        di_nfeDI: 0,
        di_last_flag: CVDIAG_SUCCESS as i64,
    });

    /* Attach linear solver memory to integrator memory */
    cv_mem.cv_lmem = LsModule::Diag(cvdiag_mem);

    CVDIAG_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * CVDiagGetWorkSpace
 * -----------------------------------------------------------------
 */
pub fn CVDiagGetWorkSpace(cv_mem: &mut CVodeMem, lenrwLS: &mut i64, leniwLS: &mut i64) -> i32 {
    *lenrwLS = 3 * cv_mem.cv_lrw1;
    *leniwLS = 3 * cv_mem.cv_liw1;

    CVDIAG_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * CVDiagGetNumRhsEvals
 * -----------------------------------------------------------------
 */
pub fn CVDiagGetNumRhsEvals(cv_mem: &mut CVodeMem, nfevalsLS: &mut i64) -> i32 {
    let cvdiag_mem = match &cv_mem.cv_lmem {
        LsModule::Diag(dm) => dm,
        _ => {
            cvProcessError(None, CVDIAG_LMEM_NULL, line!(), "CVDiagGetNumRhsEvals", file!(),
                           MSGDG_LMEM_NULL);
            return CVDIAG_LMEM_NULL;
        }
    };

    *nfevalsLS = cvdiag_mem.di_nfeDI;

    CVDIAG_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * CVDiagGetLastFlag
 * -----------------------------------------------------------------
 */
pub fn CVDiagGetLastFlag(cv_mem: &mut CVodeMem, flag: &mut i64) -> i32 {
    let cvdiag_mem = match &cv_mem.cv_lmem {
        LsModule::Diag(dm) => dm,
        _ => {
            cvProcessError(None, CVDIAG_LMEM_NULL, line!(), "CVDiagGetLastFlag", file!(),
                           MSGDG_LMEM_NULL);
            return CVDIAG_LMEM_NULL;
        }
    };

    *flag = cvdiag_mem.di_last_flag;

    CVDIAG_SUCCESS
}

/*
 * -----------------------------------------------------------------
 * CVDiagGetReturnFlagName
 * -----------------------------------------------------------------
 */
pub fn CVDiagGetReturnFlagName(flag: i64) -> String {
    let name = if flag == CVDIAG_SUCCESS as i64 {
        "CVDIAG_SUCCESS"
    } else if flag == CVDIAG_MEM_NULL as i64 {
        "CVDIAG_MEM_NULL"
    } else if flag == CVDIAG_LMEM_NULL as i64 {
        "CVDIAG_LMEM_NULL"
    } else if flag == CVDIAG_ILL_INPUT as i64 {
        "CVDIAG_ILL_INPUT"
    } else if flag == CVDIAG_MEM_FAIL as i64 {
        "CVDIAG_MEM_FAIL"
    } else if flag == CVDIAG_INV_FAIL as i64 {
        "CVDIAG_INV_FAIL"
    } else if flag == CVDIAG_RHSFUNC_UNRECVR as i64 {
        "CVDIAG_RHSFUNC_UNRECVR"
    } else if flag == CVDIAG_RHSFUNC_RECVR as i64 {
        "CVDIAG_RHSFUNC_RECVR"
    } else if flag == CVDIAG_NO_ADJ as i64 {
        "CVDIAG_NO_ADJ"
    } else {
        "NONE"
    };

    name.to_string()
}

/*
 * -----------------------------------------------------------------
 * CVDiagInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the diagonal
 * linear solver.
 * -----------------------------------------------------------------
 */
pub fn cvDiagInit(_cv_mem: &mut CVodeMem, dm: &mut CVDiagMem) -> i32 {
    dm.di_nfeDI = 0;

    dm.di_last_flag = CVDIAG_SUCCESS as i64;
    0
}

/*
 * -----------------------------------------------------------------
 * CVDiagSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the diagonal linear
 * solver.  It constructs a diagonal approximation to the Newton matrix
 * M = I - gamma*J, updates counters, and inverts M.
 *
 * C arguments map to CVodeMem fields: ypred = cv_y, fpred = cv_ftemp,
 * vtemp1 = cv_vtemp1, vtemp2 = cv_vtemp2 (vtemp3 unused).
 * -----------------------------------------------------------------
 */
pub fn cvDiagSetup(
    cv_mem: &mut CVodeMem,
    dm: &mut CVDiagMem,
    _convfail: i32,
    jcur_ptr: &mut bool,
) -> i32 {
    /* Rename work vectors for use as temporary values of y and f:
       ftemp = vtemp1 = cv_vtemp1, y = vtemp2 = cv_vtemp2 */

    /* Form y with perturbation = FRACT*(func. iter. correction) */
    let r = FRACT * cv_mem.cv_rl1;
    N_VLinearSum(
        cv_mem.cv_h,
        &cv_mem.cv_ftemp,
        -ONE,
        &cv_mem.cv_zn[1],
        &mut cv_mem.cv_vtemp1,
    );
    N_VLinearSum(r, &cv_mem.cv_vtemp1, ONE, &cv_mem.cv_y, &mut cv_mem.cv_vtemp2);

    /* Evaluate f at perturbed y */
    let f = cv_mem.cv_f.unwrap();
    let retval = f(
        cv_mem.cv_tn,
        &cv_mem.cv_vtemp2,
        &mut dm.di_M,
        &mut cv_mem.cv_user_data,
    );
    dm.di_nfeDI += 1;
    if retval < 0 {
        cvProcessError(Some(cv_mem), CVDIAG_RHSFUNC_UNRECVR, line!(), "CVDiagSetup", file!(),
                       MSGDG_RHSFUNC_FAILED);
        dm.di_last_flag = CVDIAG_RHSFUNC_UNRECVR as i64;
        return -1;
    }
    if retval > 0 {
        dm.di_last_flag = CVDIAG_RHSFUNC_RECVR as i64;
        return 1;
    }

    /* Construct M = I - gamma*J with J = diag(deltaf_i/deltay_i) */
    /* M = M - fpred */
    dm.di_M.linear_sum_with(ONE, -ONE, &cv_mem.cv_ftemp);
    /* M = FRACT*ftemp - h*M */
    dm.di_M.linear_sum_with(-(cv_mem.cv_h), FRACT, &cv_mem.cv_vtemp1);
    N_VProd(&cv_mem.cv_vtemp1, &cv_mem.cv_ewt, &mut cv_mem.cv_vtemp2);
    /* Protect against deltay_i being at roundoff level */
    N_VCompare(cv_mem.cv_uround, &cv_mem.cv_vtemp2, &mut dm.di_bit);
    N_VAddConst(&dm.di_bit, -ONE, &mut dm.di_bitcomp);
    N_VProd(&cv_mem.cv_vtemp1, &dm.di_bit, &mut cv_mem.cv_vtemp2);
    /* y = FRACT*y - bitcomp */
    cv_mem.cv_vtemp2.linear_sum_with(FRACT, -ONE, &dm.di_bitcomp);
    /* M = M / y */
    dm.di_M.div_with(&cv_mem.cv_vtemp2);
    /* M = M .* bit */
    dm.di_M.prod_with(&dm.di_bit);
    /* M = M - bitcomp */
    dm.di_M.linear_sum_with(ONE, -ONE, &dm.di_bitcomp);

    /* Invert M with test for zero components */
    let invOK = N_VInvTest_inplace(&mut dm.di_M);
    if !invOK {
        dm.di_last_flag = CVDIAG_INV_FAIL as i64;
        return 1;
    }

    /* Set jcur = SUNTRUE, save gamma in gammasv, and return */
    *jcur_ptr = SUNTRUE;
    dm.di_gammasv = cv_mem.cv_gamma;
    dm.di_last_flag = CVDIAG_SUCCESS as i64;
    0
}

/*
 * -----------------------------------------------------------------
 * CVDiagSolve
 * -----------------------------------------------------------------
 * This routine performs the solve operation for the diagonal linear
 * solver.  If necessary it first updates gamma in M = I - gamma*J.
 * -----------------------------------------------------------------
 */
pub fn cvDiagSolve(cv_mem: &mut CVodeMem, dm: &mut CVDiagMem, b: &mut NVector) -> i32 {
    /* If gamma has changed, update factor in M, and save gamma value */

    if dm.di_gammasv != cv_mem.cv_gamma {
        let r = cv_mem.cv_gamma / dm.di_gammasv;
        dm.di_M.invert_inplace();
        dm.di_M.add_const_inplace(-ONE);
        dm.di_M.scale_inplace(r);
        dm.di_M.add_const_inplace(ONE);
        let invOK = N_VInvTest_inplace(&mut dm.di_M);
        if !invOK {
            dm.di_last_flag = CVDIAG_INV_FAIL as i64;
            return 1;
        }
        dm.di_gammasv = cv_mem.cv_gamma;
    }

    /* Apply M-inverse to b */
    b.prod_with(&dm.di_M);

    dm.di_last_flag = CVDIAG_SUCCESS as i64;
    0
}

/*
 * -----------------------------------------------------------------
 * CVDiagFree
 * -----------------------------------------------------------------
 * In C this frees M, bit, and bitcomp and the memory record; here
 * dropping the LsModule::Diag value (e.g. on overwrite of cv_lmem or
 * in CVodeFree) releases everything.
 * -----------------------------------------------------------------
 */

/*
 * ================================================================
 *
 *                   PART II - backward problems
 *
 * ================================================================
 */

/*
 * CVDiagB
 *
 * Wrappers for the backward phase around the corresponding
 * CVODES functions
 */

pub fn CVDiagB(cv_mem: &mut CVodeMem, which: i32) -> i32 {
    /* Check if cvode_mem exists: impossible under &mut. */

    /* Was ASA initialized? */
    if cv_mem.cv_adjMallocDone == SUNFALSE {
        cvProcessError(Some(cv_mem), CVDIAG_NO_ADJ, line!(), "CVDiagB", file!(), MSGDG_NO_ADJ);
        return CVDIAG_NO_ADJ;
    }

    /* Check which */
    if which >= cv_mem.cv_adj_mem.as_ref().unwrap().ca_nbckpbs {
        cvProcessError(Some(cv_mem), CVDIAG_ILL_INPUT, line!(), "CVDiagB", file!(),
                       MSGDG_BAD_WHICH);
        return CVDIAG_ILL_INPUT;
    }
    let ca_mem = cv_mem.cv_adj_mem.as_deref_mut().unwrap();

    /* Find the CVodeBMem entry in the linked list corresponding to which
       (C: walk the cv_next chain; here: search the owning Vec) */
    let cvB_mem = ca_mem
        .cvB_mem
        .iter_mut()
        .find(|cvB| cvB.cv_index == which)
        .expect("backward problem for `which` not found");

    CVDiag(&mut cvB_mem.cv_mem)
}
