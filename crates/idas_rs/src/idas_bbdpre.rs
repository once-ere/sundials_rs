/* -----------------------------------------------------------------
 * Translated from src/idas/idas_bbdpre.c (IDAS 7.7.0).
 * Band-block-diagonal preconditioner (block-diagonal matrix with
 * banded blocks) for use with IDAS and the IDASLS linear solver
 * interface. This pure-Rust build is the serial reduction: a single
 * banded block of dimension n_local.
 *
 * Note: With only one process, a banded matrix results rather than a
 * b-b-d matrix with banded blocks. Diagonal blocking occurs at the
 * process level.
 *
 * PART I (forward problems) is donor-verbatim from the verified
 * ida_bbdpre.rs: the module stores its data behind
 * PrecModule::BBDPre(Box<IBBDPrecData>) in the IDALS memory and
 * idaLsPSetup / the psolve closure in idas_ls.rs dispatch to
 * IDABBDPrecSetup / IDABBDPrecSolve (C: IDASetPreconditioner with
 * pset = IDABBDPrecSetup, psolve = IDABBDPrecSolve).
 *
 * PART II (backward problems) follows the idas_ls.rs PART II pinned
 * design: IDAAglocal/IDAAgcomm are forward-callback-typed wrappers
 * whose &mut UserData downcasts to the OUTER IDAMem
 * (idaLs_AccessIDAMem); the user's glocalB/gcommB live in
 * IDABBDPrecDataB behind IDAB_mem.ida_pmem (Box<dyn Any>, the C
 * ida_pfree hook is Rust Drop); the forward solution comes from
 * idaLsGetY under the ia_noInterp gate with ia_yyTmp/ia_ypTmp taken
 * as owned locals around the user call.
 * -----------------------------------------------------------------*/
use crate::idas_bbdpre_impl::{
    IBBDPrecData, IDABBDCommFn, IDABBDCommFnB, IDABBDLocalFn, IDABBDLocalFnB, IDABBDPrecDataB,
};
use crate::idas_impl::{IDAMem, IDAProcessError, LsModule, IDA_SUCCESS};
use crate::idas_ls::{idaLsGetY, idaLs_AccessIDAMem};
use crate::idas_ls_impl::{
    PrecModule, IDALS_ILL_INPUT, IDALS_LMEM_NULL, IDALS_NO_ADJ, IDALS_PMEM_NULL, IDALS_SUCCESS,
    IDALS_SUNLS_FAIL, MSG_LS_BAD_WHICH, MSG_LS_NO_ADJ,
};
use crate::nvector_serial::{N_VClone, N_VScale, N_VSpace, NVector};
use crate::sundials_errors::SUN_SUCCESS;
use crate::sundials_linearsolver::LinearSolver;
use crate::sundials_math::{SUNMAX, SUNRabs, SUNRsqrt};
use crate::sundials_matrix::{SUNMatZero, SUNMatrix};
use crate::sundials_types::UserData;
use crate::sunlinsol_band::SUNLinSol_Band;
use crate::sunmatrix_band::SUNBandMatrixStorage;

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;

/* Error messages (idas_bbdpre_impl.h) */
const MSGBBD_LMEM_NULL: &str =
    "Linear solver memory is NULL. One of the SPILS linear solvers must be attached.";
const MSGBBD_SUNLS_FAIL: &str = "An error arose from a SUNBandLinearSolver routine.";
const MSGBBD_PMEM_NULL: &str =
    "BBD peconditioner memory is NULL. IDABBDPrecInit must be called.";
const MSGBBD_FUNC_FAILED: &str =
    "The Glocal or Gcomm routine failed in an unrecoverable manner.";
const MSGBBD_BAD_T: &str = "Bad t for interpolation.";

/*================================================================
  PART I - forward problems
  ================================================================*/

/*---------------------------------------------------------------
  User-Callable Functions: initialization, reinit and free
  ---------------------------------------------------------------*/
#[allow(clippy::too_many_arguments)]
pub fn IDABBDPrecInit(
    ida_mem: &mut IDAMem,
    Nlocal: i64,
    mudq: i64,
    mldq: i64,
    mukeep: i64,
    mlkeep: i64,
    dq_rel_yy: f64,
    Gres: Option<IDABBDLocalFn>,
    Gcomm: Option<IDABBDCommFn>,
) -> i32 {
    /* Test if the LS linear solver interface has been created */
    if !matches!(ida_mem.ida_lmem, LsModule::Ls(_)) {
        IDAProcessError(Some(ida_mem), IDALS_LMEM_NULL, line!(), "IDABBDPrecInit", file!(),
                        MSGBBD_LMEM_NULL);
        return IDALS_LMEM_NULL;
    }

    /* Test compatibility of NVECTOR package with the BBD preconditioner:
       N_VGetArrayPointer always exists for the serial NVector. */

    /* Set pointers to glocal and gcomm; load half-bandwidths. */
    let mudq_c = i64::min(Nlocal - 1, i64::max(0, mudq));
    let mldq_c = i64::min(Nlocal - 1, i64::max(0, mldq));
    let muk = i64::min(Nlocal - 1, i64::max(0, mukeep));
    let mlk = i64::min(Nlocal - 1, i64::max(0, mlkeep));

    /* Set extended upper half-bandwidth for PP (required for pivoting). */
    let storage_mu = i64::min(Nlocal - 1, muk + mlk);

    /* Allocate memory for preconditioner matrix. */
    let PP = SUNBandMatrixStorage(Nlocal, muk, mlk, storage_mu, &ida_mem.ida_sunctx);

    /* Allocate memory for temporary N_Vectors: in C zlocal/rlocal are
       empty serial wrappers that borrow the data of zvec/rvec during
       the solve; here they own their storage and IDABBDPrecSolve copies
       in and out. */
    let zlocal = NVector::new(Nlocal as usize);
    let rlocal = NVector::new(Nlocal as usize);
    let tempv1 = N_VClone(&ida_mem.ida_tempv1);
    let tempv2 = N_VClone(&ida_mem.ida_tempv1);
    let tempv3 = N_VClone(&ida_mem.ida_tempv1);
    let tempv4 = N_VClone(&ida_mem.ida_tempv1);

    /* Allocate memory for banded linear solver */
    let mut LS = SUNLinSol_Band(&rlocal, &PP, &ida_mem.ida_sunctx);

    /* initialize band linear solver object */
    let flag = LS.initialize();
    if flag != SUN_SUCCESS {
        IDAProcessError(Some(ida_mem), IDALS_SUNLS_FAIL, line!(), "IDABBDPrecInit", file!(),
                        MSGBBD_SUNLS_FAIL);
        return IDALS_SUNLS_FAIL;
    }

    /* Set rel_yy based on input value dq_rel_yy (0 implies default). */
    let rel_yy = if dq_rel_yy > ZERO {
        dq_rel_yy
    } else {
        SUNRsqrt(ida_mem.ida_uround)
    };

    /* Set work space sizes and initialize nge. */
    let mut rpwsize: i64 = 0;
    let mut ipwsize: i64 = 0;
    {
        let (lrw1, liw1) = N_VSpace(&ida_mem.ida_tempv1);
        rpwsize += 4 * lrw1;
        ipwsize += 4 * liw1;
    }
    {
        let (lrw1, liw1) = N_VSpace(&rlocal);
        rpwsize += 2 * lrw1;
        ipwsize += 2 * liw1;
    }
    {
        let (lrw, liw) = PP.space();
        rpwsize += lrw;
        ipwsize += liw;
    }
    if let LinearSolver::Band(bls) = &LS {
        let (lrw, liw) = bls.space();
        rpwsize += lrw;
        ipwsize += liw;
    }

    let pdata = Box::new(IBBDPrecData {
        mudq: mudq_c,
        mldq: mldq_c,
        mukeep: muk,
        mlkeep: mlk,
        rel_yy,
        glocal: Gres,
        gcomm: Gcomm,
        /* Store Nlocal to be used in IDABBDPrecSetup */
        n_local: Nlocal,
        PP,
        LS,
        zlocal,
        rlocal,
        tempv1,
        tempv2,
        tempv3,
        tempv4,
        rpwsize,
        ipwsize,
        nge: 0,
    });

    /* make sure pdata is free from any previous allocations (RAII on
       overwrite), point to the new pdata in the LS memory, attach the
       pfree function (RAII drop) and the preconditioner solve/setup
       functions: in C this is
       IDASetPreconditioner(idamem, IDABBDPrecSetup, IDABBDPrecSolve);
       here the prec_module drives the idaLsPSetup/psolve dispatch. */
    if let LsModule::Ls(idals_mem) = &mut ida_mem.ida_lmem {
        idals_mem.pset = None;
        idals_mem.psolve = None;
        idals_mem.prec_module = PrecModule::BBDPre(pdata);
    }

    IDALS_SUCCESS
}

/*-------------------------------------------------------------*/
pub fn IDABBDPrecReInit(ida_mem: &mut IDAMem, mudq: i64, mldq: i64, dq_rel_yy: f64) -> i32 {
    let uround = ida_mem.ida_uround;
    let idals_mem = match &mut ida_mem.ida_lmem {
        LsModule::Ls(m) => m,
        LsModule::None => {
            IDAProcessError(None, IDALS_LMEM_NULL, line!(), "IDABBDPrecReInit", file!(),
                            MSGBBD_LMEM_NULL);
            return IDALS_LMEM_NULL;
        }
    };
    let pdata = match &mut idals_mem.prec_module {
        PrecModule::BBDPre(p) => p,
        _ => {
            IDAProcessError(None, IDALS_PMEM_NULL, line!(), "IDABBDPrecReInit", file!(),
                            MSGBBD_PMEM_NULL);
            return IDALS_PMEM_NULL;
        }
    };

    /* Load half-bandwidths. */
    let Nlocal = pdata.n_local;
    pdata.mudq = i64::min(Nlocal - 1, i64::max(0, mudq));
    pdata.mldq = i64::min(Nlocal - 1, i64::max(0, mldq));

    /* Set rel_yy based on input value dq_rel_yy (0 implies default). */
    pdata.rel_yy = if dq_rel_yy > ZERO {
        dq_rel_yy
    } else {
        SUNRsqrt(uround)
    };

    /* Re-initialize nge */
    pdata.nge = 0;

    IDALS_SUCCESS
}

/*-------------------------------------------------------------*/
pub fn IDABBDPrecGetWorkSpace(
    ida_mem: &mut IDAMem,
    lenrwBBDP: &mut i64,
    leniwBBDP: &mut i64,
) -> i32 {
    let idals_mem = match &ida_mem.ida_lmem {
        LsModule::Ls(m) => m,
        LsModule::None => {
            IDAProcessError(None, IDALS_LMEM_NULL, line!(), "IDABBDPrecGetWorkSpace", file!(),
                            MSGBBD_LMEM_NULL);
            return IDALS_LMEM_NULL;
        }
    };
    let pdata = match &idals_mem.prec_module {
        PrecModule::BBDPre(p) => p,
        _ => {
            IDAProcessError(None, IDALS_PMEM_NULL, line!(), "IDABBDPrecGetWorkSpace", file!(),
                            MSGBBD_PMEM_NULL);
            return IDALS_PMEM_NULL;
        }
    };

    *lenrwBBDP = pdata.rpwsize;
    *leniwBBDP = pdata.ipwsize;

    IDALS_SUCCESS
}

/*-------------------------------------------------------------*/
pub fn IDABBDPrecGetNumGfnEvals(ida_mem: &mut IDAMem, ngevalsBBDP: &mut i64) -> i32 {
    let idals_mem = match &ida_mem.ida_lmem {
        LsModule::Ls(m) => m,
        LsModule::None => {
            IDAProcessError(None, IDALS_LMEM_NULL, line!(), "IDABBDPrecGetNumGfnEvals", file!(),
                            MSGBBD_LMEM_NULL);
            return IDALS_LMEM_NULL;
        }
    };
    let pdata = match &idals_mem.prec_module {
        PrecModule::BBDPre(p) => p,
        _ => {
            IDAProcessError(None, IDALS_PMEM_NULL, line!(), "IDABBDPrecGetNumGfnEvals", file!(),
                            MSGBBD_PMEM_NULL);
            return IDALS_PMEM_NULL;
        }
    };

    *ngevalsBBDP = pdata.nge;

    IDALS_SUCCESS
}

/*---------------------------------------------------------------
  IDABBDPrecSetup:

  Generates a band-block-diagonal preconditioner matrix (the local
  block is a band matrix computed by difference quotients via glocal,
  gcomm), then does an LU factorization in place in PP.

  The C arguments tt, yy, yp, c_j are supplied by idaLsPSetup as
  ida_tn, ycur, ypcur, ida_cj; rr is unused. bbd_data is the
  PrecModule::BBDPre payload.

  Return: 0 success, > 0 recoverable error, < 0 nonrecoverable.
 ----------------------------------------------------------------*/
pub fn IDABBDPrecSetup(
    ida_mem: &mut IDAMem,
    pdata: &mut IBBDPrecData,
    tt: f64,
    c_j: f64,
    yy: &NVector,
    yp: &NVector,
) -> i32 {
    /* Call IBBDDQJac for a new Jacobian calculation and store in PP. */
    SUNMatZero(&mut pdata.PP);
    let retval = IBBDDQJac(ida_mem, pdata, tt, c_j, yy, yp);
    if retval < 0 {
        IDAProcessError(Some(ida_mem), -1, line!(), "IDABBDPrecSetup", file!(),
                        MSGBBD_FUNC_FAILED);
        return -1;
    }
    if retval > 0 {
        return 1;
    }

    /* Do LU factorization of matrix and return error flag */
    let IBBDPrecData { LS, PP, .. } = pdata;
    LS.setup(Some(PP))
}

/*---------------------------------------------------------------
  IDABBDPrecSolve

  Computes a solution to P z = r, with the banded preconditioner
  matrix P generated and factored by IDABBDPrecSetup. r comes in as
  rvec, z is returned in zvec. The C arguments tt, yy, yp, rr, c_j
  and delta are unused. Returns the linear solver's flag.
  ---------------------------------------------------------------*/
pub fn IDABBDPrecSolve(pdata: &mut IBBDPrecData, rvec: &NVector, zvec: &mut NVector) -> i32 {
    /* Attach local data array for rvec to rlocal: the C code aliases the
       data pointers; here the serial single-block reduction copies. */
    N_VScale(ONE, rvec, &mut pdata.rlocal);

    /* Call banded solver object to do the work */
    let retval = {
        let IBBDPrecData { LS, PP, zlocal, rlocal, .. } = pdata;
        let mut atimes = |_v: &NVector, _jv: &mut NVector| -> i32 { 0 };
        LS.solve(Some(PP), zlocal, rlocal, ZERO, &mut atimes, None, None, None)
    };

    /* Copy result into zvec */
    for i in 0..(pdata.n_local as usize) {
        zvec.data[i] = pdata.zlocal.data[i];
    }

    retval
}

/* IDABBDPrecFree: the C routine frees LS, rlocal, zlocal, tempv1-4 and
   PP; here dropping the PrecModule::BBDPre value releases everything. */

/*---------------------------------------------------------------
  IBBDDQJac

  Generates a banded difference quotient approximation to the local
  block of the Jacobian of G(t,y,y'). All matrix elements are
  generated as difference quotients via calls to glocal; the number
  of calls is bandwidth + 1, where bandwidth = mldq + mudq + 1.

  gref = tempv1, ytemp = tempv2, yptemp = tempv3, gtemp = tempv4
  (as passed by IDABBDPrecSetup in C).

  Returns: 0 (success), > 0 (recoverable), < 0 (nonrecoverable).
  ----------------------------------------------------------------*/
fn IBBDDQJac(
    ida_mem: &mut IDAMem,
    pdata: &mut IBBDPrecData,
    tt: f64,
    cj: f64,
    yy: &NVector,
    yp: &NVector,
) -> i32 {
    let glocal = pdata.glocal.unwrap();

    /* Initialize ytemp and yptemp. */
    N_VScale(ONE, yy, &mut pdata.tempv2);
    N_VScale(ONE, yp, &mut pdata.tempv3);

    let hh = ida_mem.ida_hh;
    let have_constraints = ida_mem.ida_constraintsSet;

    /* Call gcomm and glocal to get base value of G(t,y,y'). */
    if let Some(gcomm) = pdata.gcomm {
        let retval = gcomm(pdata.n_local, tt, yy, yp, &mut ida_mem.ida_user_data);
        if retval != 0 {
            return retval;
        }
    }

    let retval = glocal(pdata.n_local, tt, yy, yp, &mut pdata.tempv1, &mut ida_mem.ida_user_data);
    pdata.nge += 1;
    if retval != 0 {
        return retval;
    }

    /* Set bandwidth and number of column groups for band differencing. */
    let width = pdata.mldq + pdata.mudq + 1;
    let ngroups = i64::min(width, pdata.n_local);

    /* Loop over groups. */
    for group in 1..=ngroups {
        /* Loop over the components in this group; increment yj and ypj. */
        let mut j = group - 1;
        while j < pdata.n_local {
            let ju = j as usize;
            let yj = yy.data[ju];
            let ypj = yp.data[ju];
            let ewtj = ida_mem.ida_ewt.data[ju];

            /* Set increment inc to yj based on rel_yy*abs(yj), with
               adjustments using ypj and ewtj if this is small, and a
               further adjustment to give it the same sign as hh*ypj. */
            let mut inc =
                pdata.rel_yy * SUNMAX(SUNRabs(yj), SUNMAX(SUNRabs(hh * ypj), ONE / ewtj));
            if hh * ypj < ZERO {
                inc = -inc;
            }
            inc = (yj + inc) - yj;

            /* Adjust sign(inc) again if yj has an inequality constraint. */
            if have_constraints {
                let conj = ida_mem.ida_constraints.data[ju];
                if SUNRabs(conj) == ONE {
                    if (yj + inc) * conj < ZERO {
                        inc = -inc;
                    }
                } else if SUNRabs(conj) == TWO && (yj + inc) * conj <= ZERO {
                    inc = -inc;
                }
            }

            /* Increment yj and ypj. */
            pdata.tempv2.data[ju] += inc;
            pdata.tempv3.data[ju] += cj * inc;
            j += width;
        }

        /* Evaluate G with incremented y and yp arguments. */
        let retval = {
            let IBBDPrecData { n_local, tempv2, tempv3, tempv4, .. } = &mut *pdata;
            glocal(*n_local, tt, tempv2, tempv3, tempv4, &mut ida_mem.ida_user_data)
        };
        pdata.nge += 1;
        if retval != 0 {
            return retval;
        }

        /* Loop over components of the group again; restore ytemp and
           yptemp, then form and load difference quotients into PP. */
        let pmat = match &mut pdata.PP {
            SUNMatrix::Band(m) => m,
            /* unreachable: PP is band by construction */
            _ => return -1,
        };
        let mut j = group - 1;
        while j < pdata.n_local {
            let ju = j as usize;
            let yj = yy.data[ju];
            let ypj = yp.data[ju];
            pdata.tempv2.data[ju] = yj;
            pdata.tempv3.data[ju] = ypj;
            let ewtj = ida_mem.ida_ewt.data[ju];

            /* Set increment inc as before. */
            let mut inc =
                pdata.rel_yy * SUNMAX(SUNRabs(yj), SUNMAX(SUNRabs(hh * ypj), ONE / ewtj));
            if hh * ypj < ZERO {
                inc = -inc;
            }
            inc = (yj + inc) - yj;
            if have_constraints {
                let conj = ida_mem.ida_constraints.data[ju];
                if SUNRabs(conj) == ONE {
                    if (yj + inc) * conj < ZERO {
                        inc = -inc;
                    }
                } else if SUNRabs(conj) == TWO && (yj + inc) * conj <= ZERO {
                    inc = -inc;
                }
            }

            /* Form difference quotients and load into PP. */
            let inc_inv = ONE / inc;
            let i1 = i64::max(0, j - pdata.mukeep);
            let i2 = i64::min(j + pdata.mlkeep, pdata.n_local - 1);
            for i in i1..=i2 {
                let val =
                    inc_inv * (pdata.tempv4.data[i as usize] - pdata.tempv1.data[i as usize]);
                pmat.set(i, j, val);
            }
            j += width;
        }
    }

    0
}

/*================================================================
  PART II - backward problems
  ================================================================*/

/*---------------------------------------------------------------
  User-Callable Functions: initialization, reinit and free
  ---------------------------------------------------------------*/
#[allow(clippy::too_many_arguments)]
pub fn IDABBDPrecInitB(
    ida_mem: &mut IDAMem,
    which: i32,
    NlocalB: i64,
    mudqB: i64,
    mldqB: i64,
    mukeepB: i64,
    mlkeepB: i64,
    dq_rel_yyB: f64,
    glocalB: Option<IDABBDLocalFnB>,
    gcommB: Option<IDABBDCommFnB>,
) -> i32 {
    /* (C NULL ida_mem check vanishes: &mut receiver) */

    /* Is ASA initialized? */
    if !ida_mem.ida_adjMallocDone {
        IDAProcessError(Some(ida_mem), IDALS_NO_ADJ, line!(), "IDABBDPrecInitB", file!(),
                        MSG_LS_NO_ADJ);
        return IDALS_NO_ADJ;
    }
    let idaadj_mem = ida_mem.ida_adj_mem.as_mut().unwrap();

    /* Check the value of which */
    if which >= idaadj_mem.ia_nbckpbs {
        IDAProcessError(None, IDALS_ILL_INPUT, line!(), "IDABBDPrecInitB", file!(),
                        MSG_LS_BAD_WHICH);
        return IDALS_ILL_INPUT;
    }

    /* Find the IDABMem entry in the linked list corresponding to 'which'. */
    let idx = idaadj_mem.IDAB_mem.iter().position(|b| b.ida_index == which).unwrap();
    let idaB_mem = &mut idaadj_mem.IDAB_mem[idx];

    /* Initialize the BBD preconditioner for this backward problem. */
    let flag = IDABBDPrecInit(&mut idaB_mem.IDA_mem, NlocalB, mudqB, mldqB, mukeepB, mlkeepB,
                              dq_rel_yyB, Some(IDAAglocal), Some(IDAAgcomm));
    if flag != IDA_SUCCESS {
        return flag;
    }

    /* Allocate memory for IDABBDPrecDataB to store the user-provided
       functions which will be called from the wrappers (the C malloc
       failure branch vanishes: Box::new is infallible) */
    /* set pointers to user-provided functions; attach pmem and pfree
       (pfree = IDABBDPrecFreeB in C; Rust Drop frees the Box) */
    idaB_mem.ida_pmem = Some(Box::new(IDABBDPrecDataB { glocalB, gcommB }));

    IDALS_SUCCESS
}

/*-------------------------------------------------------------*/
pub fn IDABBDPrecReInitB(
    ida_mem: &mut IDAMem,
    which: i32,
    mudqB: i64,
    mldqB: i64,
    dq_rel_yyB: f64,
) -> i32 {
    /* (C NULL ida_mem check vanishes: &mut receiver) */

    /* Is ASA initialized? */
    if !ida_mem.ida_adjMallocDone {
        IDAProcessError(Some(ida_mem), IDALS_NO_ADJ, line!(), "IDABBDPrecReInitB", file!(),
                        MSG_LS_NO_ADJ);
        return IDALS_NO_ADJ;
    }
    let idaadj_mem = ida_mem.ida_adj_mem.as_mut().unwrap();

    /* Check the value of which */
    if which >= idaadj_mem.ia_nbckpbs {
        IDAProcessError(None, IDALS_ILL_INPUT, line!(), "IDABBDPrecReInitB", file!(),
                        MSG_LS_BAD_WHICH);
        return IDALS_ILL_INPUT;
    }

    /* Find the IDABMem entry in the linked list corresponding to 'which'. */
    let idx = idaadj_mem.IDAB_mem.iter().position(|b| b.ida_index == which).unwrap();
    let idaB_mem = &mut idaadj_mem.IDAB_mem[idx];

    /* ReInitialize the BBD preconditioner for this backward problem. */
    IDABBDPrecReInit(&mut idaB_mem.IDA_mem, mudqB, mldqB, dq_rel_yyB)
}

/* (C IDABBDPrecFreeB frees IDAB_mem->ida_pmem; Rust Drop releases the
   IDABBDPrecDataB Box when the IDABMem entry is dropped or the pmem
   slot is overwritten.) */

/*----------------------------------------------------------------
  Wrapper functions
  ----------------------------------------------------------------*/

/*----------------------------------------------------------------
  IDAAglocal

  This routine interfaces to the IDALocalFnB routine
  provided by the user.  (Forward-callback-typed: the inner backward
  problem's UserData holds the OUTER IDAMem, installed by idaa.rs.)
  ----------------------------------------------------------------*/
fn IDAAglocal(
    NlocalB: i64,
    tt: f64,
    yyB: &NVector,
    ypB: &NVector,
    gvalB: &mut NVector,
    ida_mem: &mut UserData,
) -> i32 {
    let ida_mem = match idaLs_AccessIDAMem(ida_mem, "IDAAglocal") {
        Ok(m) => m,
        Err(e) => return e,
    };

    /* Get current backward problem. */
    let which = ida_mem.ida_adj_mem.as_ref().unwrap().ia_bckpbCrt.unwrap();

    /* Get forward solution from interpolation. */
    let (mut yyTmp, mut ypTmp, no_interp) = {
        let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
        (std::mem::take(&mut adj.ia_yyTmp), std::mem::take(&mut adj.ia_ypTmp), adj.ia_noInterp)
    };
    if !no_interp {
        let mut noS: Vec<NVector> = Vec::new();
        let mut noSp: Vec<NVector> = Vec::new();
        let flag = idaLsGetY(ida_mem, tt, &mut yyTmp, &mut ypTmp, &mut noS, &mut noSp);
        if flag != IDA_SUCCESS {
            let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
            adj.ia_yyTmp = yyTmp;
            adj.ia_ypTmp = ypTmp;
            IDAProcessError(Some(ida_mem), -1, line!(), "IDAAglocal", file!(), MSGBBD_BAD_T);
            return -1;
        }
    }

    /* Get the preconditioner's memory and call user's adjoint LocalFnB
       function. */
    let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
    let idaB_mem = &mut adj.IDAB_mem[which];
    let glocalB = idaB_mem
        .ida_pmem
        .as_mut()
        .unwrap()
        .downcast_mut::<IDABBDPrecDataB>()
        .unwrap()
        .glocalB
        .unwrap();
    let retval = glocalB(NlocalB, tt, &yyTmp, &ypTmp, yyB, ypB, gvalB,
                         &mut idaB_mem.ida_user_data);
    adj.ia_yyTmp = yyTmp;
    adj.ia_ypTmp = ypTmp;
    retval
}

/*----------------------------------------------------------------
  IDAAgcomm

  This routine interfaces to the IDACommFnB routine
  provided by the user.  (Forward-callback-typed, as IDAAglocal.)
  ----------------------------------------------------------------*/
fn IDAAgcomm(NlocalB: i64, tt: f64, yyB: &NVector, ypB: &NVector, ida_mem: &mut UserData) -> i32 {
    let ida_mem = match idaLs_AccessIDAMem(ida_mem, "IDAAgcomm") {
        Ok(m) => m,
        Err(e) => return e,
    };

    /* Get current backward problem and the preconditioner's memory. */
    let which = ida_mem.ida_adj_mem.as_ref().unwrap().ia_bckpbCrt.unwrap();
    let gcommB = {
        let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
        adj.IDAB_mem[which]
            .ida_pmem
            .as_mut()
            .unwrap()
            .downcast_mut::<IDABBDPrecDataB>()
            .unwrap()
            .gcommB
    };
    let gcommB = match gcommB {
        Some(g) => g,
        None => return 0,
    };

    /* Get forward solution from interpolation. */
    let (mut yyTmp, mut ypTmp, no_interp) = {
        let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
        (std::mem::take(&mut adj.ia_yyTmp), std::mem::take(&mut adj.ia_ypTmp), adj.ia_noInterp)
    };
    if !no_interp {
        let mut noS: Vec<NVector> = Vec::new();
        let mut noSp: Vec<NVector> = Vec::new();
        let flag = idaLsGetY(ida_mem, tt, &mut yyTmp, &mut ypTmp, &mut noS, &mut noSp);
        if flag != IDA_SUCCESS {
            let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
            adj.ia_yyTmp = yyTmp;
            adj.ia_ypTmp = ypTmp;
            IDAProcessError(Some(ida_mem), -1, line!(), "IDAAgcomm", file!(), MSGBBD_BAD_T);
            return -1;
        }
    }

    /* Call user's adjoint CommFnB routine */
    let adj = ida_mem.ida_adj_mem.as_mut().unwrap();
    let idaB_mem = &mut adj.IDAB_mem[which];
    let retval = gcommB(NlocalB, tt, &yyTmp, &ypTmp, yyB, ypB, &mut idaB_mem.ida_user_data);
    adj.ia_yyTmp = yyTmp;
    adj.ia_ypTmp = ypTmp;
    retval
}

/*==================================================================
  Tests
  ==================================================================*/
#[cfg(test)]
mod tests {
    use super::*;
    use crate::idas_ls::{idaLsPSetup, IDASetLinearSolver};
    use crate::sundials_linearsolver::SUN_PREC_RIGHT;
    use crate::sundials_types::UserData;
    use crate::sunlinsol_spgmr::SUNLinSol_SPGMR;

    /* analytic tridiagonal residual in y only: G_i = 4 y_i - y_{i-1}
       - y_{i+1}, so dG/dy = tridiag(-1, 4, -1) and dG/dy' = 0, hence
       the (y,y') DQ block is tridiag(-1,4,-1) independent of cj. */
    fn glocal_tridiag(
        Nlocal: i64,
        _tt: f64,
        yy: &NVector,
        _yp: &NVector,
        gval: &mut NVector,
        _ud: &mut UserData,
    ) -> i32 {
        let n = Nlocal as usize;
        for i in 0..n {
            let mut v = 4.0 * yy.data[i];
            if i > 0 {
                v -= yy.data[i - 1];
            }
            if i + 1 < n {
                v -= yy.data[i + 1];
            }
            gval.data[i] = v;
        }
        0
    }

    fn gcomm_count(
        _Nlocal: i64,
        _tt: f64,
        _yy: &NVector,
        _yp: &NVector,
        user_data: &mut UserData,
    ) -> i32 {
        if let Some(d) = user_data {
            if let Some(c) = d.downcast_mut::<i64>() {
                *c += 1;
            }
        }
        0
    }

    fn make_ida_mem(n: usize) -> IDAMem {
        let mut ida_mem = IDAMem::default();
        ida_mem.ida_tempv1 = NVector::new(n);
        ida_mem.ida_ewt = NVector::from_slice(&vec![1.0; n]);
        ida_mem.ida_hh = 1.0;
        ida_mem.ida_cj = 1.0;
        ida_mem.ida_tn = 0.0;
        ida_mem
    }

    /* IDABBDPrecInit requires the IDALS interface (IDALS_LMEM_NULL) */
    #[test]
    fn idabbdprecinit_needs_lmem() {
        let mut ida_mem = make_ida_mem(4);
        let flag = IDABBDPrecInit(&mut ida_mem, 4, 1, 1, 1, 1, 0.0, Some(glocal_tridiag), None);
        assert_eq!(flag, IDALS_LMEM_NULL);
    }

    /* banded DQ Jacobian against the analytic tridiagonal glocal, then
       the band LU solve must round-trip: solve P z = A x and get x */
    #[test]
    fn idabbdprec_dq_jacobian_and_solve_roundtrip() {
        const N: usize = 5;
        let mut ida_mem = make_ida_mem(N);
        let sunctx = ida_mem.ida_sunctx.clone();

        /* iterative LS (SPGMR) so the preconditioner is legal */
        let ls = SUNLinSol_SPGMR(&ida_mem.ida_tempv1, SUN_PREC_RIGHT, 0, &sunctx);
        assert_eq!(IDASetLinearSolver(&mut ida_mem, ls, None), IDALS_SUCCESS);

        /* count gcomm calls through user_data */
        ida_mem.ida_user_data = Some(Box::new(0i64));

        let flag = IDABBDPrecInit(&mut ida_mem, N as i64, 1, 1, 1, 1, 0.0,
                                  Some(glocal_tridiag), Some(gcomm_count));
        assert_eq!(flag, IDALS_SUCCESS);

        /* run the setup through the IDALS dispatch (idaLsPSetup) */
        let yy = NVector::from_slice(&(1..=N).map(|i| 0.5 + i as f64).collect::<Vec<_>>());
        let yp = NVector::new(N);
        let rr = NVector::new(N);
        let mut lmem = std::mem::take(&mut ida_mem.ida_lmem);
        let idals_mem = match &mut lmem {
            LsModule::Ls(m) => m,
            _ => unreachable!(),
        };
        assert_eq!(idaLsPSetup(&mut ida_mem, idals_mem, &yy, &yp, &rr), 0);
        assert_eq!(idals_mem.npe, 1);

        let pdata = match &mut idals_mem.prec_module {
            PrecModule::BBDPre(p) => p,
            _ => unreachable!(),
        };

        /* base G eval + one eval per column group (width = mldq+mudq+1 = 3) */
        assert_eq!(pdata.nge, 1 + 3);

        /* band LU solve round-trip: r = A*x with A = tridiag(-1,4,-1) */
        let x = [1.0, -2.0, 3.0, -4.0, 5.0];
        let mut rvec = NVector::new(N);
        let mut zvec = NVector::new(N);
        for i in 0..N {
            let mut v = 4.0 * x[i];
            if i > 0 {
                v -= x[i - 1];
            }
            if i + 1 < N {
                v -= x[i + 1];
            }
            rvec.data[i] = v;
        }
        assert_eq!(IDABBDPrecSolve(pdata, &rvec, &mut zvec), 0);
        for i in 0..N {
            assert!((zvec.data[i] - x[i]).abs() < 1.0e-6, "z[{}] = {} != {}", i, zvec.data[i], x[i]);
        }

        ida_mem.ida_lmem = lmem;

        /* gcomm was called once, before the base glocal evaluation */
        let ncomm = *ida_mem.ida_user_data.as_ref().unwrap().downcast_ref::<i64>().unwrap();
        assert_eq!(ncomm, 1);

        /* optional output getters */
        let mut nge = -1;
        assert_eq!(IDABBDPrecGetNumGfnEvals(&mut ida_mem, &mut nge), IDALS_SUCCESS);
        assert_eq!(nge, 4);
        let (mut lenrw, mut leniw) = (-1, -1);
        assert_eq!(IDABBDPrecGetWorkSpace(&mut ida_mem, &mut lenrw, &mut leniw), IDALS_SUCCESS);
        assert!(lenrw > 0);
    }

    /* fully-decoupled residual G_i = 2 y_i: the DQ block is 2 I, so
       P z = r gives z = r/2 — exercises IDABBDPrecReInit too */
    fn glocal_diag(
        Nlocal: i64,
        _tt: f64,
        yy: &NVector,
        _yp: &NVector,
        gval: &mut NVector,
        _ud: &mut UserData,
    ) -> i32 {
        for i in 0..Nlocal as usize {
            gval.data[i] = 2.0 * yy.data[i];
        }
        0
    }

    #[test]
    fn idabbdprecreinit_resets_counters_and_bandwidths() {
        const N: usize = 4;
        let mut ida_mem = make_ida_mem(N);
        let sunctx = ida_mem.ida_sunctx.clone();
        let ls = SUNLinSol_SPGMR(&ida_mem.ida_tempv1, SUN_PREC_RIGHT, 0, &sunctx);
        assert_eq!(IDASetLinearSolver(&mut ida_mem, ls, None), IDALS_SUCCESS);
        assert_eq!(IDABBDPrecInit(&mut ida_mem, N as i64, 2, 2, 0, 0, 0.0,
                                  Some(glocal_diag), None),
                   IDALS_SUCCESS);

        /* reinit narrows the DQ bandwidths and resets nge */
        assert_eq!(IDABBDPrecReInit(&mut ida_mem, 0, 0, 0.5), IDALS_SUCCESS);
        if let LsModule::Ls(m) = &mut ida_mem.ida_lmem {
            if let PrecModule::BBDPre(p) = &mut m.prec_module {
                assert_eq!(p.mudq, 0);
                assert_eq!(p.mldq, 0);
                assert_eq!(p.rel_yy, 0.5);
                assert_eq!(p.nge, 0);
            } else {
                unreachable!();
            }
        }
    }

    /* backward-problem entry points demand an initialized adjoint
       module (IDALS_NO_ADJ), per idas_bbdpre.c PART II */
    fn glocalB_stub(
        _NlocalB: i64,
        _tt: f64,
        _yy: &NVector,
        _yp: &NVector,
        _yyB: &NVector,
        _ypB: &NVector,
        _gvalB: &mut NVector,
        _ud: &mut UserData,
    ) -> i32 {
        0
    }

    #[test]
    fn backward_entry_points_require_adjoint() {
        let mut ida_mem = make_ida_mem(4);
        assert_eq!(IDABBDPrecInitB(&mut ida_mem, 0, 4, 1, 1, 1, 1, 0.0,
                                   Some(glocalB_stub), None),
                   IDALS_NO_ADJ);
        assert_eq!(IDABBDPrecReInitB(&mut ida_mem, 0, 1, 1, 0.0), IDALS_NO_ADJ);
    }
}
