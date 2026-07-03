/* -----------------------------------------------------------------
 * Translated from src/cvodes/cvodes_bandpre.c (CVODES 7.7.0).
 * Banded difference-quotient Jacobian-based preconditioner and
 * solver routines for use with the CVSLS linear solver interface.
 *
 * In C the module stores its data behind cvls_mem->P_data and
 * installs cvBandPrecSetup/cvBandPrecSolve via
 * CVodeSetPreconditioner; here CVBandPrecInit installs
 * PrecModule::BandPre(Box<CVBandPrecData>) in the CVLS memory and
 * cvLsPSetup / the psolve closure in cvodes_ls.rs dispatch to
 * CVBandPrecSetup / CVBandPrecSolve.
 * -----------------------------------------------------------------*/
use crate::cvodes_bandpre_impl::CVBandPrecData;
use crate::cvodes_impl::{cvProcessError, CVodeMem, LsModule};
use crate::cvodes_ls_impl::{
    PrecModule, CVLS_ILL_INPUT, CVLS_LMEM_NULL, CVLS_NO_ADJ, CVLS_PMEM_NULL, CVLS_SUCCESS,
    CVLS_SUNLS_FAIL,
};
use crate::nvector_serial::{N_VClone, N_VScale, N_VSpace, N_VWrmsNorm, NVector};
use crate::sundials_errors::SUN_SUCCESS;
use crate::sundials_linearsolver::LinearSolver;
use crate::sundials_math::{SUNMAX, SUNRabs, SUNRsqrt};
use crate::sundials_matrix::{SUNMatCopy, SUNMatScaleAddI, SUNMatZero, SUNMatrix};
use crate::sundials_types::{SUNFALSE, SUNTRUE};
use crate::sunlinsol_band::SUNLinSol_Band;
use crate::sunmatrix_band::SUNBandMatrixStorage;

const MIN_INC_MULT: f64 = 1000.0;
const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;

/* Error messages (cvodes_bandpre_impl.h) */
const MSGBP_LMEM_NULL: &str =
    "Linear solver memory is NULL. One of the SPILS linear solvers must be attached.";
const MSGBP_SUNMAT_FAIL: &str = "An error arose from a SUNBandMatrix routine.";
const MSGBP_SUNLS_FAIL: &str = "An error arose from a SUNBandLinearSolver routine.";
const MSGBP_PMEM_NULL: &str =
    "Band preconditioner memory is NULL. CVBandPrecInit must be called.";
const MSGBP_RHSFUNC_FAILED: &str =
    "The right-hand side routine failed in an unrecoverable manner.";
const MSGBP_NO_ADJ: &str = "Illegal attempt to call before calling CVodeAdjInit.";
const MSGBP_BAD_WHICH: &str = "Illegal value for parameter which.";

/*================================================================
  PART I - Forward Problems
  ================================================================*/

/*-----------------------------------------------------------------
  Initialization, Free, and Get Functions
  NOTE: The band linear solver assumes a serial implementation of
        the NVECTOR package (always true here).
  -----------------------------------------------------------------*/
pub fn CVBandPrecInit(cv_mem: &mut CVodeMem, n: i64, mu: i64, ml: i64) -> i32 {
    /* Test if the CVSLS linear solver interface has been attached */
    if !matches!(cv_mem.cv_lmem, LsModule::Ls(_)) {
        cvProcessError(Some(cv_mem), CVLS_LMEM_NULL, line!(), "CVBandPrecInit", file!(),
                       MSGBP_LMEM_NULL);
        return CVLS_LMEM_NULL;
    }

    /* Test compatibility of NVECTOR package with the BAND preconditioner:
       N_VGetArrayPointer always exists for the serial NVector. */

    /* Load pointers and bandwidths into pdata block. */
    let N = n;
    let mup = i64::min(N - 1, i64::max(0, mu));
    let mlp = i64::min(N - 1, i64::max(0, ml));

    /* Allocate memory for saved banded Jacobian approximation. */
    let savedJ = SUNBandMatrixStorage(N, mup, mlp, mup, &cv_mem.cv_sunctx);

    /* Allocate memory for banded preconditioner. */
    let storagemu = i64::min(N - 1, mup + mlp);
    let savedP = SUNBandMatrixStorage(N, mup, mlp, storagemu, &cv_mem.cv_sunctx);

    /* Allocate memory for banded linear solver */
    let mut LS = SUNLinSol_Band(&cv_mem.cv_tempv, &savedP, &cv_mem.cv_sunctx);

    /* allocate memory for temporary N_Vectors */
    let tmp1 = N_VClone(&cv_mem.cv_tempv);
    let tmp2 = N_VClone(&cv_mem.cv_tempv);

    /* initialize band linear solver object */
    let flag = LS.initialize();
    if flag != SUN_SUCCESS {
        cvProcessError(Some(cv_mem), CVLS_SUNLS_FAIL, line!(), "CVBandPrecInit", file!(),
                       MSGBP_SUNLS_FAIL);
        return CVLS_SUNLS_FAIL;
    }

    let pdata = Box::new(CVBandPrecData {
        N,
        ml: mlp,
        mu: mup,
        savedJ,
        savedP,
        LS,
        tmp1,
        tmp2,
        /* Initialize nfeBP counter */
        nfeBP: 0,
    });

    /* make sure P_data is free from any previous allocations (RAII on
       overwrite), point to the new P_data field in the LS memory, and
       attach the preconditioner setup and solve functions: in C this is
       CVodeSetPreconditioner(cvode_mem, cvBandPrecSetup, cvBandPrecSolve);
       here the prec_module drives the cvLsPSetup/cvLsPSolve dispatch. */
    let mut ill_input = false;
    if let LsModule::Ls(cvls_mem) = &mut cv_mem.cv_lmem {
        cvls_mem.pset = None;
        cvls_mem.psolve = None;
        cvls_mem.prec_module = PrecModule::BandPre(pdata);

        /* CVodeSetPreconditioner issues an error if the LS object does
           not allow user-supplied preconditioning */
        ill_input = !matches!(
            cvls_mem.LS,
            LinearSolver::Spgmr(_)
                | LinearSolver::Spfgmr(_)
                | LinearSolver::Spbcgs(_)
                | LinearSolver::Sptfqmr(_)
                | LinearSolver::Pcg(_)
        );
    }
    if ill_input {
        cvProcessError(Some(cv_mem), CVLS_ILL_INPUT, line!(), "CVBandPrecInit", file!(),
                       "SUNLinearSolver object does not support user-supplied preconditioning");
        return CVLS_ILL_INPUT;
    }

    CVLS_SUCCESS
}

pub fn CVBandPrecGetWorkSpace(cv_mem: &mut CVodeMem, lenrwBP: &mut i64, leniwBP: &mut i64) -> i32 {
    let cvls_mem = match &cv_mem.cv_lmem {
        LsModule::Ls(ls) => ls,
        _ => {
            cvProcessError(None, CVLS_LMEM_NULL, line!(), "CVBandPrecGetWorkSpace", file!(),
                           MSGBP_LMEM_NULL);
            return CVLS_LMEM_NULL;
        }
    };
    let pdata = match &cvls_mem.prec_module {
        PrecModule::BandPre(bp) => bp,
        _ => {
            cvProcessError(None, CVLS_PMEM_NULL, line!(), "CVBandPrecGetWorkSpace", file!(),
                           MSGBP_PMEM_NULL);
            return CVLS_PMEM_NULL;
        }
    };

    /* sum space requirements for all objects in pdata */
    *leniwBP = 4;
    *lenrwBP = 0;
    {
        let (lrw1, liw1) = N_VSpace(&cv_mem.cv_tempv);
        *leniwBP += 2 * liw1;
        *lenrwBP += 2 * lrw1;
    }
    {
        let (lrw, liw) = pdata.savedJ.space();
        *leniwBP += liw;
        *lenrwBP += lrw;
    }
    {
        let (lrw, liw) = pdata.savedP.space();
        *leniwBP += liw;
        *lenrwBP += lrw;
    }
    if let LinearSolver::Band(bls) = &pdata.LS {
        let (lrw, liw) = bls.space();
        *leniwBP += liw;
        *lenrwBP += lrw;
    }

    CVLS_SUCCESS
}

pub fn CVBandPrecGetNumRhsEvals(cv_mem: &mut CVodeMem, nfevalsBP: &mut i64) -> i32 {
    let cvls_mem = match &cv_mem.cv_lmem {
        LsModule::Ls(ls) => ls,
        _ => {
            cvProcessError(None, CVLS_LMEM_NULL, line!(), "CVBandPrecGetNumRhsEvals", file!(),
                           MSGBP_LMEM_NULL);
            return CVLS_LMEM_NULL;
        }
    };
    let pdata = match &cvls_mem.prec_module {
        PrecModule::BandPre(bp) => bp,
        _ => {
            cvProcessError(None, CVLS_PMEM_NULL, line!(), "CVBandPrecGetNumRhsEvals", file!(),
                           MSGBP_PMEM_NULL);
            return CVLS_PMEM_NULL;
        }
    };

    *nfevalsBP = pdata.nfeBP;

    CVLS_SUCCESS
}

/*-----------------------------------------------------------------
  cvBandPrecSetup
  -----------------------------------------------------------------
  Together cvBandPrecSetup and cvBandPrecSolve use a banded
  difference quotient Jacobian to create a preconditioner.
  cvBandPrecSetup calculates a new J, if necessary, then
  calculates P = I - gamma*J, and does an LU factorization of P.

  C arguments map to CVodeMem fields: t = cv_tn, y = cv_y,
  fy = cv_ftemp, gamma = cv_gamma, and *jcurPtr = cv_jcur (as passed
  by cvLsPSetup).

  The value returned is
    0  if successful, or
    1  if the band factorization failed (recoverable),
   <0  on an unrecoverable error.
  -----------------------------------------------------------------*/
pub fn CVBandPrecSetup(cv_mem: &mut CVodeMem, bp: &mut CVBandPrecData, jok: bool) -> i32 {
    /* Assume matrix and lpivots have already been allocated. */

    if jok {
        /* If jok = SUNTRUE, use saved copy of J. */
        cv_mem.cv_jcur = SUNFALSE;
        let retval = SUNMatCopy(&bp.savedJ, &mut bp.savedP);
        if retval < 0 {
            cvProcessError(Some(cv_mem), -1, line!(), "cvBandPrecSetup", file!(),
                           MSGBP_SUNMAT_FAIL);
            return -1;
        }
        if retval > 0 {
            return 1;
        }
    } else {
        /* If jok = SUNFALSE, call CVBandPDQJac for new J value. */
        cv_mem.cv_jcur = SUNTRUE;
        let retval = SUNMatZero(&mut bp.savedJ);
        if retval < 0 {
            cvProcessError(Some(cv_mem), -1, line!(), "cvBandPrecSetup", file!(),
                           MSGBP_SUNMAT_FAIL);
            return -1;
        }
        if retval > 0 {
            return 1;
        }

        let retval = cvBandPrecDQJac(cv_mem, bp);
        if retval < 0 {
            cvProcessError(Some(cv_mem), -1, line!(), "cvBandPrecSetup", file!(),
                           MSGBP_RHSFUNC_FAILED);
            return -1;
        }
        if retval > 0 {
            return 1;
        }

        let retval = SUNMatCopy(&bp.savedJ, &mut bp.savedP);
        if retval < 0 {
            cvProcessError(Some(cv_mem), -1, line!(), "cvBandPrecSetup", file!(),
                           MSGBP_SUNMAT_FAIL);
            return -1;
        }
        if retval > 0 {
            return 1;
        }
    }

    /* Scale and add identity to get savedP = I - gamma*J. */
    let retval = SUNMatScaleAddI(-cv_mem.cv_gamma, &mut bp.savedP);
    if retval != 0 {
        cvProcessError(Some(cv_mem), -1, line!(), "cvBandPrecSetup", file!(),
                       MSGBP_SUNMAT_FAIL);
        return -1;
    }

    /* Do LU factorization of matrix and return error flag */
    let CVBandPrecData { LS, savedP, .. } = bp;
    LS.setup(Some(savedP))
}

/*-----------------------------------------------------------------
  cvBandPrecSolve
  -----------------------------------------------------------------
  cvBandPrecSolve solves a linear system P z = r, where P is the
  matrix computed by cvBandPrecSetup.

  The value returned is 0 on success (as returned by the band
  SUNLinSolSolve).
  -----------------------------------------------------------------*/
pub fn CVBandPrecSolve(
    _cv_mem: &mut CVodeMem,
    bp: &mut CVBandPrecData,
    r: &NVector,
    z: &mut NVector,
) -> i32 {
    /* Call banded solver object to do the work */
    let CVBandPrecData { LS, savedP, .. } = bp;
    let mut atimes = |_v: &NVector, _jv: &mut NVector| -> i32 { 0 };
    LS.solve(Some(savedP), z, r, ZERO, &mut atimes, None, None, None)
}

/* cvBandPrecFree: the C routine frees LS, savedP, savedJ, tmp1, tmp2;
   here dropping the PrecModule::BandPre value releases everything. */

/*-----------------------------------------------------------------
  cvBandPrecDQJac
  -----------------------------------------------------------------
  This routine generates a banded difference quotient approximation
  to the Jacobian of f(t,y). It assumes that a band SUNMatrix is
  stored column-wise, and that elements within each column are
  contiguous.

  C arguments map as: t = cv_tn, y = cv_y, fy = cv_ftemp,
  ftemp = pdata.tmp1, ytemp = pdata.tmp2.
  (CVODES tests constraint presence through the cv_constraints
  vector itself: empty = C NULL.)
  -----------------------------------------------------------------*/
fn cvBandPrecDQJac(cv_mem: &mut CVodeMem, pdata: &mut CVBandPrecData) -> i32 {
    let t = cv_mem.cv_tn;
    let f = cv_mem.cv_f.unwrap();
    let constraints_set = !cv_mem.cv_constraints.data.is_empty();

    /* Load ytemp with y = predicted y vector. */
    N_VScale(ONE, &cv_mem.cv_y, &mut pdata.tmp2);

    /* Set minimum increment based on uround and norm of f. */
    let srur = SUNRsqrt(cv_mem.cv_uround);
    let fnorm = N_VWrmsNorm(&cv_mem.cv_ftemp, &cv_mem.cv_ewt);
    let minInc = if fnorm != ZERO {
        MIN_INC_MULT * SUNRabs(cv_mem.cv_h) * cv_mem.cv_uround * pdata.N as f64 * fnorm
    } else {
        ONE
    };

    /* Set bandwidth and number of column groups for band differencing. */
    let width = pdata.ml + pdata.mu + 1;
    let ngroups = i64::min(width, pdata.N);

    for group in 1..=ngroups {
        /* Increment all y_j in group. */
        let mut j = group - 1;
        while j < pdata.N {
            let ju = j as usize;
            let mut inc = SUNMAX(
                srur * SUNRabs(cv_mem.cv_y.data[ju]),
                minInc / cv_mem.cv_ewt.data[ju],
            );
            let yj = cv_mem.cv_y.data[ju];

            /* Adjust sign(inc) again if yj has an inequality constraint. */
            if constraints_set {
                let conj = cv_mem.cv_constraints.data[ju];
                if SUNRabs(conj) == ONE {
                    if (yj + inc) * conj < ZERO {
                        inc = -inc;
                    }
                } else if SUNRabs(conj) == TWO && (yj + inc) * conj <= ZERO {
                    inc = -inc;
                }
            }

            pdata.tmp2.data[ju] += inc;
            j += width;
        }

        /* Evaluate f with incremented y. */
        let retval = f(t, &pdata.tmp2, &mut pdata.tmp1, &mut cv_mem.cv_user_data);
        pdata.nfeBP += 1;
        if retval != 0 {
            return retval;
        }

        /* Restore ytemp, then form and load difference quotients. */
        let jmat = match &mut pdata.savedJ {
            SUNMatrix::Band(m) => m,
            /* unreachable: savedJ is band by construction */
            _ => return -1,
        };
        let mut j = group - 1;
        while j < pdata.N {
            let ju = j as usize;
            let yj = cv_mem.cv_y.data[ju];
            pdata.tmp2.data[ju] = cv_mem.cv_y.data[ju];
            let mut inc = SUNMAX(
                srur * SUNRabs(cv_mem.cv_y.data[ju]),
                minInc / cv_mem.cv_ewt.data[ju],
            );

            /* Adjust sign(inc) as before. */
            if constraints_set {
                let conj = cv_mem.cv_constraints.data[ju];
                if SUNRabs(conj) == ONE {
                    if (yj + inc) * conj < ZERO {
                        inc = -inc;
                    }
                } else if SUNRabs(conj) == TWO && (yj + inc) * conj <= ZERO {
                    inc = -inc;
                }
            }

            let inc_inv = ONE / inc;
            let i1 = i64::max(0, j - pdata.mu);
            let i2 = i64::min(j + pdata.ml, pdata.N - 1);
            for i in i1..=i2 {
                /* SM_COLUMN_ELEMENT_B(col_j, i, j) */
                let val =
                    inc_inv * (pdata.tmp1.data[i as usize] - cv_mem.cv_ftemp.data[i as usize]);
                jmat.set(i, j, val);
            }
            j += width;
        }
    }

    0
}

/*================================================================
  PART II - Backward Problems
  ================================================================*/

/*---------------------------------------------------------------
  User-Callable initialization function: wrapper for the backward
  phase around the corresponding CVODES functions
  ---------------------------------------------------------------*/
pub fn CVBandPrecInitB(cv_mem: &mut CVodeMem, which: i32, nB: i64, muB: i64, mlB: i64) -> i32 {
    /* Was ASA initialized? */
    if !cv_mem.cv_adjMallocDone {
        cvProcessError(Some(cv_mem), CVLS_NO_ADJ, line!(), "CVBandPrecInitB", file!(),
                       MSGBP_NO_ADJ);
        return CVLS_NO_ADJ;
    }

    /* Check which */
    if which >= cv_mem.cv_adj_mem.as_ref().unwrap().ca_nbckpbs {
        cvProcessError(Some(cv_mem), CVLS_ILL_INPUT, line!(), "CVBandPrecInitB", file!(),
                       MSGBP_BAD_WHICH);
        return CVLS_ILL_INPUT;
    }
    let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();

    /* Find the CVodeBMem entry corresponding to which (C: linked-list
       search by cv_index; backward problems are created with
       cv_index = 0..nbckpbs-1, so the entry always exists) */
    let idx = ca_mem.cvB_mem.iter().position(|b| b.cv_index == which).unwrap();
    let cvB_mem = &mut ca_mem.cvB_mem[idx];

    /* Set pfree: the C code nulls cvB_mem->cv_pfree; the cv_pfree
       pointer is dropped in this port (RAII) and any previously
       attached preconditioner pmem is released here instead. */
    cvB_mem.cv_pmem = None;

    /* Initialize the band preconditioner for this backward problem. */
    CVBandPrecInit(&mut cvB_mem.cv_mem, nB, muB, mlB)
}
