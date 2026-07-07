/* -----------------------------------------------------------------
 * Translated from src/ida/ida_bbdpre.c (IDA 7.7.0).
 * Band-block-diagonal preconditioner (block-diagonal matrix with
 * banded blocks) for use with IDA and the IDALS linear solver
 * interface. This pure-Rust build is the serial reduction: a single
 * banded block of dimension n_local.
 *
 * Note: With only one process, a banded matrix results rather than a
 * b-b-d matrix with banded blocks. Diagonal blocking occurs at the
 * process level.
 *
 * In C the module stores its data behind idals_mem->pdata and installs
 * IDABBDPrecSetup/IDABBDPrecSolve via IDASetPreconditioner; here
 * IDABBDPrecInit installs PrecModule::BBDPre(Box<IBBDPrecData>) in the
 * IDALS memory and idaLsPSetup / the psolve closure in ida_ls.rs
 * dispatch to IDABBDPrecSetup / IDABBDPrecSolve.
 * -----------------------------------------------------------------*/
use crate::ida_bbdpre_impl::{IBBDPrecData, IDABBDCommFn, IDABBDLocalFn};
use crate::ida_impl::{IDAMem, IDAProcessError, LsModule};
use crate::ida_ls_impl::{
    PrecModule, IDALS_LMEM_NULL, IDALS_PMEM_NULL, IDALS_SUCCESS, IDALS_SUNLS_FAIL,
};
use crate::nvector_serial::{N_VClone, N_VScale, N_VSpace, NVector};
use crate::sundials_errors::SUN_SUCCESS;
use crate::sundials_linearsolver::LinearSolver;
use crate::sundials_math::{SUNMAX, SUNRabs, SUNRsqrt};
use crate::sundials_matrix::{SUNMatZero, SUNMatrix};
use crate::sunlinsol_band::SUNLinSol_Band;
use crate::sunmatrix_band::SUNBandMatrixStorage;

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;

/* Error messages (ida_bbdpre_impl.h) */
const MSGBBD_LMEM_NULL: &str =
    "Linear solver memory is NULL. One of the SPILS linear solvers must be attached.";
const MSGBBD_SUNLS_FAIL: &str = "An error arose from a SUNBandLinearSolver routine.";
const MSGBBD_PMEM_NULL: &str =
    "BBD peconditioner memory is NULL. IDABBDPrecInit must be called.";
const MSGBBD_FUNC_FAILED: &str =
    "The Glocal or Gcomm routine failed in an unrecoverable manner.";

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

/*==================================================================
  Tests
  ==================================================================*/
#[cfg(test)]
mod tests {
    use super::*;
    use crate::ida_ls::{idaLsPSetup, IDASetLinearSolver};
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
}
