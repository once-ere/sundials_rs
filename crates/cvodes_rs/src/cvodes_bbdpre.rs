/* -----------------------------------------------------------------
 * Translated from src/cvodes/cvodes_bbdpre.c (CVODES 7.7.0).
 * Band-block-diagonal preconditioner, i.e. a block-diagonal matrix
 * with banded blocks, for use with CVODE and the CVSLS linear
 * solver interface. This pure-Rust build is the serial reduction:
 * a single block of dimension n_local.
 *
 * In C the module stores its data behind cvls_mem->P_data and
 * installs cvBBDPrecSetup/cvBBDPrecSolve via CVodeSetPreconditioner;
 * here CVBBDPrecInit installs PrecModule::BBDPre(Box<CVBBDPrecData>)
 * in the CVLS memory and cvLsPSetup / the psolve closure in
 * cvodes_ls.rs dispatch to CVBBDPrecSetup / CVBBDPrecSolve.
 * -----------------------------------------------------------------*/
use crate::cvodes_bbdpre_impl::{CVBBDPrecData, CVBBDPrecDataB, CVCommFn, CVCommFnB, CVLocalFn,
                                CVLocalFnB};
use crate::cvodes_impl::{cvProcessError, CVodeMem, LsModule, CV_HERMITE, CV_SUCCESS};
use crate::cvodes_ls_impl::{
    PrecModule, CVLS_ILL_INPUT, CVLS_LMEM_NULL, CVLS_MEM_NULL, CVLS_NO_ADJ, CVLS_PMEM_NULL,
    CVLS_SUCCESS, CVLS_SUNLS_FAIL,
};
use crate::nvector_serial::{N_VClone, N_VScale, N_VSpace, N_VWrmsNorm, NVector};
use crate::sundials_errors::SUN_SUCCESS;
use crate::sundials_linearsolver::LinearSolver;
use crate::sundials_math::{SUNMAX, SUNRabs, SUNRsqrt};
use crate::sundials_matrix::{SUNMatCopy, SUNMatScaleAddI, SUNMatZero, SUNMatrix};
use crate::sundials_types::{SUNFALSE, SUNTRUE, UserData};
use crate::sunlinsol_band::SUNLinSol_Band;
use crate::sunmatrix_band::SUNBandMatrixStorage;

const MIN_INC_MULT: f64 = 1000.0;
const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;

/* Error messages (cvodes_bbdpre_impl.h) */
const MSGBBD_MEM_NULL: &str = "Integrator memory is NULL.";
const MSGBBD_LMEM_NULL: &str =
    "Linear solver memory is NULL. One of the SPILS linear solvers must be attached.";
const MSGBBD_SUNMAT_FAIL: &str = "An error arose from a SUNBandMatrix routine.";
const MSGBBD_SUNLS_FAIL: &str = "An error arose from a SUNBandLinearSolver routine.";
const MSGBBD_PMEM_NULL: &str =
    "BBD peconditioner memory is NULL. CVBBDPrecInit must be called.";
const MSGBBD_FUNC_FAILED: &str = "The gloc or cfn routine failed in an unrecoverable manner.";
const MSGBBD_NO_ADJ: &str = "Illegal attempt to call before calling CVodeAdjInit.";
const MSGBBD_BAD_WHICH: &str = "Illegal value for the which parameter.";
const MSGBBD_BAD_TINTERP: &str = "Bad t for interpolation.";

/*================================================================
  PART I - forward problems
  ================================================================*/

/*-----------------------------------------------------------------
  User-Callable Functions: initialization, reinit and free
  -----------------------------------------------------------------*/
#[allow(clippy::too_many_arguments)]
pub fn CVBBDPrecInit(
    cv_mem: &mut CVodeMem,
    nlocal: i64,
    mudq: i64,
    mldq: i64,
    mukeep: i64,
    mlkeep: i64,
    dqrely: f64,
    gloc: Option<CVLocalFn>,
    cfn: Option<CVCommFn>,
) -> i32 {
    /* Test if the CVSLS linear solver interface has been created */
    if !matches!(cv_mem.cv_lmem, LsModule::Ls(_)) {
        cvProcessError(Some(cv_mem), CVLS_LMEM_NULL, line!(), "CVBBDPrecInit", file!(),
                       MSGBBD_LMEM_NULL);
        return CVLS_LMEM_NULL;
    }

    /* Test compatibility of NVECTOR package with the BBD preconditioner:
       N_VGetArrayPointer always exists for the serial NVector. */

    /* Set pointers to gloc and cfn; load half-bandwidths */
    let mudq_c = i64::min(nlocal - 1, i64::max(0, mudq));
    let mldq_c = i64::min(nlocal - 1, i64::max(0, mldq));
    let muk = i64::min(nlocal - 1, i64::max(0, mukeep));
    let mlk = i64::min(nlocal - 1, i64::max(0, mlkeep));

    /* Allocate memory for saved Jacobian */
    let savedJ = SUNBandMatrixStorage(nlocal, muk, mlk, muk, &cv_mem.cv_sunctx);

    /* Allocate memory for preconditioner matrix */
    let storage_mu = i64::min(nlocal - 1, muk + mlk);
    let savedP = SUNBandMatrixStorage(nlocal, muk, mlk, storage_mu, &cv_mem.cv_sunctx);

    /* Allocate memory for temporary N_Vectors: in C zlocal/rlocal are
       empty serial wrappers that borrow the data of z/r during the
       solve; here they own their storage and CVBBDPrecSolve copies. */
    let zlocal = NVector::new(nlocal as usize);
    let rlocal = NVector::new(nlocal as usize);
    let tmp1 = N_VClone(&cv_mem.cv_tempv);
    let tmp2 = N_VClone(&cv_mem.cv_tempv);
    let tmp3 = N_VClone(&cv_mem.cv_tempv);

    /* Allocate memory for banded linear solver */
    let mut LS = SUNLinSol_Band(&rlocal, &savedP, &cv_mem.cv_sunctx);

    /* initialize band linear solver object */
    let flag = LS.initialize();
    if flag != SUN_SUCCESS {
        cvProcessError(Some(cv_mem), CVLS_SUNLS_FAIL, line!(), "CVBBDPrecInit", file!(),
                       MSGBBD_SUNLS_FAIL);
        return CVLS_SUNLS_FAIL;
    }

    /* Set dqrely based on input dqrely (0 implies default). */
    let dqrely = if dqrely > ZERO {
        dqrely
    } else {
        SUNRsqrt(cv_mem.cv_uround)
    };

    /* Set work space sizes and initialize nge */
    let mut rpwsize: i64 = 0;
    let mut ipwsize: i64 = 0;
    {
        let (lrw1, liw1) = N_VSpace(&cv_mem.cv_tempv);
        rpwsize += 3 * lrw1;
        ipwsize += 3 * liw1;
    }
    {
        let (lrw1, liw1) = N_VSpace(&rlocal);
        rpwsize += 2 * lrw1;
        ipwsize += 2 * liw1;
    }
    {
        let (lrw, liw) = savedJ.space();
        rpwsize += lrw;
        ipwsize += liw;
    }
    {
        let (lrw, liw) = savedP.space();
        rpwsize += lrw;
        ipwsize += liw;
    }
    if let LinearSolver::Band(bls) = &LS {
        let (lrw, liw) = bls.space();
        rpwsize += lrw;
        ipwsize += liw;
    }

    let pdata = Box::new(CVBBDPrecData {
        mudq: mudq_c,
        mldq: mldq_c,
        mukeep: muk,
        mlkeep: mlk,
        dqrely,
        gloc,
        cfn,
        savedJ,
        savedP,
        LS,
        tmp1,
        tmp2,
        tmp3,
        zlocal,
        rlocal,
        /* Store Nlocal to be used in CVBBDPrecSetup */
        n_local: nlocal,
        rpwsize,
        ipwsize,
        nge: 0,
    });

    /* make sure P_data is free from any previous allocations (RAII on
       overwrite), point to the new P_data field in the LS memory, and
       attach the preconditioner setup and solve functions: in C this is
       CVodeSetPreconditioner(cvode_mem, cvBBDPrecSetup, cvBBDPrecSolve);
       here the prec_module drives the cvLsPSetup/cvLsPSolve dispatch. */
    let mut ill_input = false;
    if let LsModule::Ls(cvls_mem) = &mut cv_mem.cv_lmem {
        cvls_mem.pset = None;
        cvls_mem.psolve = None;
        cvls_mem.prec_module = PrecModule::BBDPre(pdata);

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
        cvProcessError(Some(cv_mem), CVLS_ILL_INPUT, line!(), "CVBBDPrecInit", file!(),
                       "SUNLinearSolver object does not support user-supplied preconditioning");
        return CVLS_ILL_INPUT;
    }

    CVLS_SUCCESS
}

pub fn CVBBDPrecReInit(cv_mem: &mut CVodeMem, mudq: i64, mldq: i64, dqrely: f64) -> i32 {
    let uround = cv_mem.cv_uround;

    /* Test if the LS linear solver interface has been created */
    let cvls_mem = match &mut cv_mem.cv_lmem {
        LsModule::Ls(ls) => ls,
        _ => {
            cvProcessError(None, CVLS_LMEM_NULL, line!(), "CVBBDPrecReInit", file!(),
                           MSGBBD_LMEM_NULL);
            return CVLS_LMEM_NULL;
        }
    };

    /* Test if the preconditioner data is non-NULL */
    let pdata = match &mut cvls_mem.prec_module {
        PrecModule::BBDPre(p) => p,
        _ => {
            cvProcessError(None, CVLS_PMEM_NULL, line!(), "CVBBDPrecReInit", file!(),
                           MSGBBD_PMEM_NULL);
            return CVLS_PMEM_NULL;
        }
    };

    /* Load half-bandwidths */
    let nlocal = pdata.n_local;
    pdata.mudq = i64::min(nlocal - 1, i64::max(0, mudq));
    pdata.mldq = i64::min(nlocal - 1, i64::max(0, mldq));

    /* Set pdata->dqrely based on input dqrely (0 implies default). */
    pdata.dqrely = if dqrely > ZERO { dqrely } else { SUNRsqrt(uround) };

    /* Re-initialize nge */
    pdata.nge = 0;

    CVLS_SUCCESS
}

pub fn CVBBDPrecGetWorkSpace(
    cv_mem: &mut CVodeMem,
    lenrwBBDP: &mut i64,
    leniwBBDP: &mut i64,
) -> i32 {
    let cvls_mem = match &cv_mem.cv_lmem {
        LsModule::Ls(ls) => ls,
        _ => {
            cvProcessError(None, CVLS_LMEM_NULL, line!(), "CVBBDPrecGetWorkSpace", file!(),
                           MSGBBD_LMEM_NULL);
            return CVLS_LMEM_NULL;
        }
    };
    let pdata = match &cvls_mem.prec_module {
        PrecModule::BBDPre(p) => p,
        _ => {
            cvProcessError(None, CVLS_PMEM_NULL, line!(), "CVBBDPrecGetWorkSpace", file!(),
                           MSGBBD_PMEM_NULL);
            return CVLS_PMEM_NULL;
        }
    };

    *lenrwBBDP = pdata.rpwsize;
    *leniwBBDP = pdata.ipwsize;

    CVLS_SUCCESS
}

pub fn CVBBDPrecGetNumGfnEvals(cv_mem: &mut CVodeMem, ngevalsBBDP: &mut i64) -> i32 {
    let cvls_mem = match &cv_mem.cv_lmem {
        LsModule::Ls(ls) => ls,
        _ => {
            cvProcessError(None, CVLS_LMEM_NULL, line!(), "CVBBDPrecGetNumGfnEvals", file!(),
                           MSGBBD_LMEM_NULL);
            return CVLS_LMEM_NULL;
        }
    };
    let pdata = match &cvls_mem.prec_module {
        PrecModule::BBDPre(p) => p,
        _ => {
            cvProcessError(None, CVLS_PMEM_NULL, line!(), "CVBBDPrecGetNumGfnEvals", file!(),
                           MSGBBD_PMEM_NULL);
            return CVLS_PMEM_NULL;
        }
    };

    *ngevalsBBDP = pdata.nge;

    CVLS_SUCCESS
}

/*-----------------------------------------------------------------
  Function : cvBBDPrecSetup
  -----------------------------------------------------------------
  CVBBDPrecSetup generates and factors a banded block of the
  preconditioner matrix, via calls to the user-supplied gloc and cfn
  functions. It uses difference quotient approximations to the
  Jacobian elements.

  CVBBDPrecSetup calculates a new J, if necessary, then calculates
  P = I - gamma*J, and does an LU factorization of P.

  C arguments map to CVodeMem fields: t = cv_tn, y = cv_y,
  fy = cv_ftemp (unused), gamma = cv_gamma, and *jcurPtr = cv_jcur
  (as passed by cvLsPSetup).

  Return value:
    0  if successful,
    1  for a recoverable error (step will be retried),
   <0  on an unrecoverable error.
  -----------------------------------------------------------------*/
pub fn CVBBDPrecSetup(cv_mem: &mut CVodeMem, bbd: &mut CVBBDPrecData, jok: bool) -> i32 {
    /* If jok = SUNTRUE, use saved copy of J */
    if jok {
        cv_mem.cv_jcur = SUNFALSE;
        let retval = SUNMatCopy(&bbd.savedJ, &mut bbd.savedP);
        if retval < 0 {
            cvProcessError(Some(cv_mem), -1, line!(), "cvBBDPrecSetup", file!(),
                           MSGBBD_SUNMAT_FAIL);
            return -1;
        }
        if retval > 0 {
            return 1;
        }

        /* Otherwise call cvBBDDQJac for new J value */
    } else {
        cv_mem.cv_jcur = SUNTRUE;
        let retval = SUNMatZero(&mut bbd.savedJ);
        if retval < 0 {
            cvProcessError(Some(cv_mem), -1, line!(), "cvBBDPrecSetup", file!(),
                           MSGBBD_SUNMAT_FAIL);
            return -1;
        }
        if retval > 0 {
            return 1;
        }

        let retval = cvBBDDQJac(cv_mem, bbd);
        if retval < 0 {
            cvProcessError(Some(cv_mem), -1, line!(), "cvBBDPrecSetup", file!(),
                           MSGBBD_FUNC_FAILED);
            return -1;
        }
        if retval > 0 {
            return 1;
        }

        let retval = SUNMatCopy(&bbd.savedJ, &mut bbd.savedP);
        if retval < 0 {
            cvProcessError(Some(cv_mem), -1, line!(), "cvBBDPrecSetup", file!(),
                           MSGBBD_SUNMAT_FAIL);
            return -1;
        }
        if retval > 0 {
            return 1;
        }
    }

    /* Scale and add I to get P = I - gamma*J */
    let retval = SUNMatScaleAddI(-cv_mem.cv_gamma, &mut bbd.savedP);
    if retval != 0 {
        cvProcessError(Some(cv_mem), -1, line!(), "cvBBDPrecSetup", file!(),
                       MSGBBD_SUNMAT_FAIL);
        return -1;
    }

    /* Do LU factorization of matrix and return error flag */
    let CVBBDPrecData { LS, savedP, .. } = bbd;
    LS.setup(Some(savedP))
}

/*-----------------------------------------------------------------
  Function : cvBBDPrecSolve
  -----------------------------------------------------------------
  CVBBDPrecSolve solves a linear system P z = r, with the
  band-block-diagonal preconditioner matrix P generated and factored
  by CVBBDPrecSetup.

  The value returned is 0 on success (as returned by the band
  SUNLinSolSolve).
  -----------------------------------------------------------------*/
pub fn CVBBDPrecSolve(
    _cv_mem: &mut CVodeMem,
    bbd: &mut CVBBDPrecData,
    r: &NVector,
    z: &mut NVector,
) -> i32 {
    /* Attach local data arrays for r and z to rlocal and zlocal:
       the C code aliases the data pointers; here the serial
       single-block reduction copies in and out instead. */
    N_VScale(ONE, r, &mut bbd.rlocal);

    /* Call banded solver object to do the work */
    let retval = {
        let CVBBDPrecData {
            LS,
            savedP,
            zlocal,
            rlocal,
            ..
        } = bbd;
        let mut atimes = |_v: &NVector, _jv: &mut NVector| -> i32 { 0 };
        LS.solve(Some(savedP), zlocal, rlocal, ZERO, &mut atimes, None, None, None)
    };

    /* Detach local data arrays from rlocal and zlocal */
    N_VScale(ONE, &bbd.zlocal, z);

    retval
}

/* cvBBDPrecFree: the C routine frees LS, tmp1-3, zlocal, rlocal,
   savedP, savedJ; here dropping the PrecModule::BBDPre value releases
   everything. */

/*-----------------------------------------------------------------
  Function : cvBBDDQJac
  -----------------------------------------------------------------
  This routine generates a banded difference quotient approximation
  to the local block of the Jacobian of g(t,y). It assumes that a
  band SUNMatrix is stored columnwise, and that elements within each
  column are contiguous. All matrix elements are generated as
  difference quotients, by way of calls to the user routine gloc.
  By virtue of the band structure, the number of these calls is
  bandwidth + 1, where bandwidth = mldq + mudq + 1.
  But the band matrix kept has bandwidth = mlkeep + mukeep + 1.

  C arguments map as: t = cv_tn, y = cv_y, gy = pdata.tmp1,
  ytemp = pdata.tmp2, gtemp = pdata.tmp3.
  (CVODES tests constraint presence through the cv_constraints
  vector itself: empty = C NULL.)
  -----------------------------------------------------------------*/
fn cvBBDDQJac(cv_mem: &mut CVodeMem, pdata: &mut CVBBDPrecData) -> i32 {
    let t = cv_mem.cv_tn;
    let gloc = pdata.gloc.unwrap();
    let constraints_set = !cv_mem.cv_constraints.data.is_empty();

    /* Load ytemp with y = predicted solution vector */
    N_VScale(ONE, &cv_mem.cv_y, &mut pdata.tmp2);

    /* Call cfn and gloc to get base value of g(t,y) */
    if let Some(cfn) = pdata.cfn {
        let retval = cfn(pdata.n_local, t, &cv_mem.cv_y, &mut cv_mem.cv_user_data);
        if retval != 0 {
            return retval;
        }
    }

    let retval = gloc(
        pdata.n_local,
        t,
        &pdata.tmp2,
        &mut pdata.tmp1,
        &mut cv_mem.cv_user_data,
    );
    pdata.nge += 1;
    if retval != 0 {
        return retval;
    }

    /* Set minimum increment based on uround and norm of g */
    let gnorm = N_VWrmsNorm(&pdata.tmp1, &cv_mem.cv_ewt);
    let minInc = if gnorm != ZERO {
        MIN_INC_MULT * SUNRabs(cv_mem.cv_h) * cv_mem.cv_uround * pdata.n_local as f64 * gnorm
    } else {
        ONE
    };

    /* Set bandwidth and number of column groups for band differencing */
    let width = pdata.mldq + pdata.mudq + 1;
    let ngroups = i64::min(width, pdata.n_local);

    /* Loop over groups */
    for group in 1..=ngroups {
        /* Increment all y_j in group */
        let mut j = group - 1;
        while j < pdata.n_local {
            let ju = j as usize;
            let mut inc = SUNMAX(
                pdata.dqrely * SUNRabs(cv_mem.cv_y.data[ju]),
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

        /* Evaluate g with incremented y */
        let retval = gloc(
            pdata.n_local,
            t,
            &pdata.tmp2,
            &mut pdata.tmp3,
            &mut cv_mem.cv_user_data,
        );
        pdata.nge += 1;
        if retval != 0 {
            return retval;
        }

        /* Restore ytemp, then form and load difference quotients */
        let jmat = match &mut pdata.savedJ {
            SUNMatrix::Band(m) => m,
            /* unreachable: savedJ is band by construction */
            _ => return -1,
        };
        let mut j = group - 1;
        while j < pdata.n_local {
            let ju = j as usize;
            let yj = cv_mem.cv_y.data[ju];
            pdata.tmp2.data[ju] = cv_mem.cv_y.data[ju];
            let mut inc = SUNMAX(
                pdata.dqrely * SUNRabs(cv_mem.cv_y.data[ju]),
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
            let i1 = i64::max(0, j - pdata.mukeep);
            let i2 = i64::min(j + pdata.mlkeep, pdata.n_local - 1);
            for i in i1..=i2 {
                /* SM_COLUMN_ELEMENT_B(col_j, i, j) */
                let val =
                    inc_inv * (pdata.tmp3.data[i as usize] - pdata.tmp1.data[i as usize]);
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
  User-Callable Functions: initialization, reinit and free
  ---------------------------------------------------------------*/
#[allow(clippy::too_many_arguments)]
pub fn CVBBDPrecInitB(
    cv_mem: &mut CVodeMem,
    which: i32,
    NlocalB: i64,
    mudqB: i64,
    mldqB: i64,
    mukeepB: i64,
    mlkeepB: i64,
    dqrelyB: f64,
    glocB: Option<CVLocalFnB>,
    cfnB: Option<CVCommFnB>,
) -> i32 {
    /* Was ASA initialized? */
    if !cv_mem.cv_adjMallocDone {
        cvProcessError(Some(cv_mem), CVLS_NO_ADJ, line!(), "CVBBDPrecInitB", file!(),
                       MSGBBD_NO_ADJ);
        return CVLS_NO_ADJ;
    }

    /* Check which */
    if which >= cv_mem.cv_adj_mem.as_ref().unwrap().ca_nbckpbs {
        cvProcessError(Some(cv_mem), CVLS_ILL_INPUT, line!(), "CVBBDPrecInitB", file!(),
                       MSGBBD_BAD_WHICH);
        return CVLS_ILL_INPUT;
    }

    /* Find the CVodeBMem entry corresponding to which (C: linked-list
       search by cv_index; backward problems are created with
       cv_index = 0..nbckpbs-1, so the entry always exists) */
    let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
    let idx = ca_mem.cvB_mem.iter().position(|b| b.cv_index == which).unwrap();

    /* Initialize the BBD preconditioner for this backward problem. */
    let flag = CVBBDPrecInit(
        &mut ca_mem.cvB_mem[idx].cv_mem,
        NlocalB,
        mudqB,
        mldqB,
        mukeepB,
        mlkeepB,
        dqrelyB,
        Some(cvGlocWrapper),
        Some(cvCfnWrapper),
    );
    if flag != CV_SUCCESS {
        return flag;
    }

    /* Allocate memory for CVBBDPrecDataB to store the user-provided
       functions which will be called from the wrappers */
    let cvbbdB_mem = Box::new(CVBBDPrecDataB {
        /* set pointers to user-provided functions */
        glocB,
        cfnB,
    });

    /* Attach pmem (the C pfree attachment, CVBBDPrecFreeB, is
       subsumed by Drop) */
    let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
    ca_mem.cvB_mem[idx].cv_pmem = Some(cvbbdB_mem);

    CVLS_SUCCESS
}

pub fn CVBBDPrecReInitB(
    cv_mem: &mut CVodeMem,
    which: i32,
    mudqB: i64,
    mldqB: i64,
    dqrelyB: f64,
) -> i32 {
    /* Was ASA initialized? */
    if !cv_mem.cv_adjMallocDone {
        cvProcessError(Some(cv_mem), CVLS_NO_ADJ, line!(), "CVBBDPrecReInitB", file!(),
                       MSGBBD_NO_ADJ);
        return CVLS_NO_ADJ;
    }

    /* Check which */
    if which >= cv_mem.cv_adj_mem.as_ref().unwrap().ca_nbckpbs {
        cvProcessError(Some(cv_mem), CVLS_ILL_INPUT, line!(), "CVBBDPrecReInitB", file!(),
                       MSGBBD_BAD_WHICH);
        return CVLS_ILL_INPUT;
    }

    /* Find the CVodeBMem entry corresponding to which */
    let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
    let idx = ca_mem.cvB_mem.iter().position(|b| b.cv_index == which).unwrap();

    /* ReInitialize the BBD preconditioner for this backward problem. */
    CVBBDPrecReInit(&mut ca_mem.cvB_mem[idx].cv_mem, mudqB, mldqB, dqrelyB)
}

/* CVBBDPrecFreeB: the C routine frees cvB_mem->cv_pmem; here dropping
   the CVodeBMem (or overwriting cv_pmem) releases it. */

/*----------------------------------------------------------------
  Wrapper functions
  ----------------------------------------------------------------*/

/* The C wrappers receive the *forward* integrator memory through the
   void* user_data of the backward problem (cvodea.c: CVodeCreateB
   does CVodeSetUserData(cvodeB_mem, cvode_mem)); same cross-module
   contract as the cvodes_ls.rs backward wrappers. */
fn cvBBD_AccessCVodeMem<'a>(
    cvode_mem: &'a mut UserData,
    fname: &str,
) -> Result<&'a mut CVodeMem, i32> {
    match cvode_mem.as_mut().and_then(|d| d.downcast_mut::<CVodeMem>()) {
        Some(m) => Ok(m),
        None => {
            cvProcessError(None, CVLS_MEM_NULL, line!(), fname, file!(), MSGBBD_MEM_NULL);
            Err(CVLS_MEM_NULL)
        }
    }
}

/* ca_mem->ca_IMget(cv_mem, t, y, NULL): dispatch on ca_IMtype per the
   pinned interpolation-module design (cvodes_impl.rs decision 4). An
   empty yS Vec plays the role of the C NULL argument.
   FORWARD REFERENCE: cvaHermiteGetY / cvaPolynomialGetY are
   implemented by the cvodea.c port (cvodea.rs). */
fn cvBBDIMget(cv_mem: &mut CVodeMem, t: f64, y: &mut NVector, yS: &mut Vec<NVector>) -> i32 {
    if cv_mem.cv_adj_mem.as_ref().unwrap().ca_IMtype == CV_HERMITE {
        crate::cvodea::cvaHermiteGetY(cv_mem, t, y, yS)
    } else {
        crate::cvodea::cvaPolynomialGetY(cv_mem, t, y, yS)
    }
}

/* cvGlocWrapper interfaces to the CVLocalFnB routine provided by the
   user. cvGlocWrapper is of type CVLocalFn. */
fn cvGlocWrapper(
    NlocalB: i64,
    t: f64,
    yB: &NVector,
    gB: &mut NVector,
    cvode_mem: &mut UserData,
) -> i32 {
    /* access relevant memory structures */
    let cv_mem = match cvBBD_AccessCVodeMem(cvode_mem, "cvGlocWrapper") {
        Ok(m) => m,
        Err(e) => return e,
    };
    let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
    /* C: cvB_mem = ca_mem->ca_bckpbCrt (set for the duration of the
       backward run; the C code dereferences it unconditionally) */
    let which = ca_mem.ca_bckpbCrt.unwrap();

    /* Get forward solution from interpolation */
    let mut ytmp = std::mem::take(&mut ca_mem.ca_ytmp);
    let mut noS: Vec<NVector> = Vec::new();
    let flag = cvBBDIMget(cv_mem, t, &mut ytmp, &mut noS);
    if flag != CV_SUCCESS {
        cv_mem.cv_adj_mem.as_mut().unwrap().ca_ytmp = ytmp;
        cvProcessError(Some(cv_mem), -1, line!(), "cvGlocWrapper", file!(), MSGBBD_BAD_TINTERP);
        return -1;
    }

    /* Call user's adjoint glocB routine */
    let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
    let cvB_mem = &mut ca_mem.cvB_mem[which];
    let glocB = cvB_mem
        .cv_pmem
        .as_mut()
        .unwrap()
        .downcast_mut::<CVBBDPrecDataB>()
        .unwrap()
        .glocB
        .unwrap();
    let retval = glocB(NlocalB, t, &ytmp, yB, gB, &mut cvB_mem.cv_user_data);
    ca_mem.ca_ytmp = ytmp;
    retval
}

/* cvCfnWrapper interfaces to the CVCommFnB routine provided by the
   user. cvCfnWrapper is of type CVCommFn. */
fn cvCfnWrapper(NlocalB: i64, t: f64, yB: &NVector, cvode_mem: &mut UserData) -> i32 {
    /* access relevant memory structures */
    let cv_mem = match cvBBD_AccessCVodeMem(cvode_mem, "cvCfnWrapper") {
        Ok(m) => m,
        Err(e) => return e,
    };
    let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
    let which = ca_mem.ca_bckpbCrt.unwrap();
    let cfnB = ca_mem.cvB_mem[which]
        .cv_pmem
        .as_ref()
        .unwrap()
        .downcast_ref::<CVBBDPrecDataB>()
        .unwrap()
        .cfnB;
    let cfnB = match cfnB {
        Some(f) => f,
        None => return 0,
    };

    /* Get forward solution from interpolation */
    let mut ytmp = std::mem::take(&mut ca_mem.ca_ytmp);
    let mut noS: Vec<NVector> = Vec::new();
    let flag = cvBBDIMget(cv_mem, t, &mut ytmp, &mut noS);
    if flag != CV_SUCCESS {
        cv_mem.cv_adj_mem.as_mut().unwrap().ca_ytmp = ytmp;
        cvProcessError(Some(cv_mem), -1, line!(), "cvCfnWrapper", file!(), MSGBBD_BAD_TINTERP);
        return -1;
    }

    /* Call user's adjoint cfnB routine */
    let ca_mem = cv_mem.cv_adj_mem.as_mut().unwrap();
    let cvB_mem = &mut ca_mem.cvB_mem[which];
    let retval = cfnB(NlocalB, t, &ytmp, yB, &mut cvB_mem.cv_user_data);
    ca_mem.ca_ytmp = ytmp;
    retval
}
