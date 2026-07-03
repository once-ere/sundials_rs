/* -----------------------------------------------------------------
 * Translated from src/kinsol/kinsol_bbdpre.c (KINSOL 7.7.0).
 * Band-block-diagonal preconditioner, i.e. a block-diagonal matrix
 * with banded blocks, for use with KINSol and the KINLS linear
 * solver interface. This pure-Rust build is the serial reduction:
 * a single block of dimension n_local.
 *
 * Note: With only one process, a banded matrix results rather than
 * a b-b-d matrix with banded blocks. Diagonal blocking occurs at
 * the process level.
 *
 * In C the module stores its data behind kinls_mem->pdata and
 * installs KINBBDPrecSetup/KINBBDPrecSolve via KINSetPreconditioner;
 * here KINBBDPrecInit installs PrecModule::BBDPre(Box<KBBDPrecData>)
 * in the KINLS memory and kinLsPSetup / the psolve closure in
 * kinsol_ls.rs dispatch to KINBBDPrecSetup / KINBBDPrecSolve.
 * -----------------------------------------------------------------*/
use crate::kinsol_bbdpre_impl::{KBBDPrecData, KINBBDCommFn, KINBBDLocalFn};
use crate::kinsol_impl::{KINMem, KINProcessError, LsModule};
use crate::kinsol_ls_impl::{
    PrecModule, KINLS_ILL_INPUT, KINLS_LMEM_NULL, KINLS_PMEM_NULL, KINLS_SUCCESS, KINLS_SUNLS_FAIL,
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

/* Error messages (kinsol_bbdpre_impl.h) */
const MSGBBD_LMEM_NULL: &str =
    "Linear solver memory is NULL. One of the SPILS linear solvers must be attached.";
const MSGBBD_SUNMAT_FAIL: &str = "An error arose from a SUNBandMatrix routine.";
const MSGBBD_SUNLS_FAIL: &str = "An error arose from a SUNBandLinearSolver routine.";
const MSGBBD_PMEM_NULL: &str =
    "BBD peconditioner memory is NULL. IDABBDPrecInit must be called.";
const MSGBBD_FUNC_FAILED: &str = "The gloc or gcomm routine failed in an unrecoverable manner.";

/*------------------------------------------------------------------
  user-callable functions
  ------------------------------------------------------------------*/

/*------------------------------------------------------------------
  KINBBDPrecInit
  ------------------------------------------------------------------*/
pub fn KINBBDPrecInit(
    kin_mem: &mut KINMem,
    Nlocal: i64,
    mudq: i64,
    mldq: i64,
    mukeep: i64,
    mlkeep: i64,
    dq_rel_uu: f64,
    gloc: Option<KINBBDLocalFn>,
    gcomm: Option<KINBBDCommFn>,
) -> i32 {
    /* Test if the LS linear solver interface has been created */
    if !matches!(kin_mem.kin_lmem, LsModule::Ls(_)) {
        KINProcessError(Some(kin_mem), KINLS_LMEM_NULL, line!(), "KINBBDPrecInit", file!(),
                        MSGBBD_LMEM_NULL);
        return KINLS_LMEM_NULL;
    }

    /* Test compatibility of NVECTOR package with the BBD preconditioner:
       N_VGetArrayPointer always exists for the serial NVector. */

    /* Set pointers to gloc and gcomm; load half-bandwidths */
    let mudq_c = i64::min(Nlocal - 1, i64::max(0, mudq));
    let mldq_c = i64::min(Nlocal - 1, i64::max(0, mldq));
    let muk = i64::min(Nlocal - 1, i64::max(0, mukeep));
    let mlk = i64::min(Nlocal - 1, i64::max(0, mlkeep));

    /* Set extended upper half-bandwidth for PP (required for pivoting) */
    let storage_mu = i64::min(Nlocal - 1, muk + mlk);

    /* Allocate memory for preconditioner matrix */
    let PP = SUNBandMatrixStorage(Nlocal, muk, mlk, storage_mu, &kin_mem.kin_sunctx);

    /* Allocate memory for temporary N_Vectors: in C rlocal is an empty
       serial wrapper (N_VNewEmpty_Serial) that borrows the data of vv
       during the solve; here it owns its storage and KINBBDPrecSolve
       copies in and out. */
    let zlocal = NVector::new(Nlocal as usize);
    let rlocal = NVector::new(Nlocal as usize);
    let tempv1 = N_VClone(&kin_mem.kin_vtemp1);
    let tempv2 = N_VClone(&kin_mem.kin_vtemp1);
    let tempv3 = N_VClone(&kin_mem.kin_vtemp1);

    /* Allocate memory for banded linear solver */
    let mut LS = SUNLinSol_Band(&zlocal, &PP, &kin_mem.kin_sunctx);

    /* initialize band linear solver object */
    let flag = LS.initialize();
    if flag != SUN_SUCCESS {
        KINProcessError(Some(kin_mem), KINLS_SUNLS_FAIL, line!(), "KINBBDPrecInit", file!(),
                        MSGBBD_SUNLS_FAIL);
        return KINLS_SUNLS_FAIL;
    }

    /* Set rel_uu based on input value dq_rel_uu (0 implies default) */
    let rel_uu = if dq_rel_uu > ZERO {
        dq_rel_uu
    } else {
        SUNRsqrt(kin_mem.kin_uround)
    };

    /* Set work space sizes and initialize nge */
    let mut rpwsize: i64 = 0;
    let mut ipwsize: i64 = 0;
    {
        let (lrw1, liw1) = N_VSpace(&kin_mem.kin_vtemp1);
        rpwsize += 3 * lrw1;
        ipwsize += 3 * liw1;
    }
    {
        let (lrw1, liw1) = N_VSpace(&zlocal);
        rpwsize += lrw1;
        ipwsize += liw1;
    }
    {
        let (lrw1, liw1) = N_VSpace(&rlocal);
        rpwsize += lrw1;
        ipwsize += liw1;
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

    let pdata = Box::new(KBBDPrecData {
        gloc,
        gcomm,
        mudq: mudq_c,
        mldq: mldq_c,
        mukeep: muk,
        mlkeep: mlk,
        rel_uu,
        /* Store Nlocal to be used in KINBBDPrecSetup */
        n_local: Nlocal,
        PP,
        LS,
        rlocal,
        zlocal,
        tempv1,
        tempv2,
        tempv3,
        rpwsize,
        ipwsize,
        nge: 0,
    });

    /* make sure pdata is free from any previous allocations (RAII on
       overwrite), point to the new pdata field in the LS memory, attach
       the pfree function (RAII drop), and attach the preconditioner
       solve and setup functions: in C this is
       KINSetPreconditioner(kinmem, KINBBDPrecSetup, KINBBDPrecSolve);
       here the prec_module drives the kinLsPSetup/psolve dispatch. */
    let mut ill_input = false;
    if let LsModule::Ls(kinls_mem) = &mut kin_mem.kin_lmem {
        kinls_mem.pset = None;
        kinls_mem.psolve = None;
        kinls_mem.prec_module = PrecModule::BBDPre(pdata);

        /* KINSetPreconditioner issues an error if the LS object does
           not allow user-supplied preconditioning */
        ill_input = matches!(
            kinls_mem.LS,
            LinearSolver::Dense(_) | LinearSolver::Band(_) | LinearSolver::Custom(_)
        );
    }
    if ill_input {
        KINProcessError(Some(kin_mem), KINLS_ILL_INPUT, line!(), "KINBBDPrecInit", file!(),
                        "SUNLinearSolver object does not support user-supplied preconditioning");
        return KINLS_ILL_INPUT;
    }

    KINLS_SUCCESS
}

/*------------------------------------------------------------------
  KINBBDPrecGetWorkSpace
  ------------------------------------------------------------------*/
pub fn KINBBDPrecGetWorkSpace(
    kin_mem: &mut KINMem,
    lenrwBBDP: &mut i64,
    leniwBBDP: &mut i64,
) -> i32 {
    let kinls_mem = match &kin_mem.kin_lmem {
        LsModule::Ls(ls) => ls,
        _ => {
            KINProcessError(None, KINLS_LMEM_NULL, line!(), "KINBBDPrecGetWorkSpace", file!(),
                            MSGBBD_LMEM_NULL);
            return KINLS_LMEM_NULL;
        }
    };
    let pdata = match &kinls_mem.prec_module {
        PrecModule::BBDPre(p) => p,
        _ => {
            KINProcessError(None, KINLS_PMEM_NULL, line!(), "KINBBDPrecGetWorkSpace", file!(),
                            MSGBBD_PMEM_NULL);
            return KINLS_PMEM_NULL;
        }
    };

    *lenrwBBDP = pdata.rpwsize;
    *leniwBBDP = pdata.ipwsize;

    KINLS_SUCCESS
}

/*------------------------------------------------------------------
 KINBBDPrecGetNumGfnEvals
 -------------------------------------------------------------------*/
pub fn KINBBDPrecGetNumGfnEvals(kin_mem: &mut KINMem, ngevalsBBDP: &mut i64) -> i32 {
    let kinls_mem = match &kin_mem.kin_lmem {
        LsModule::Ls(ls) => ls,
        _ => {
            KINProcessError(None, KINLS_LMEM_NULL, line!(), "KINBBDPrecGetNumGfnEvals", file!(),
                            MSGBBD_LMEM_NULL);
            return KINLS_LMEM_NULL;
        }
    };
    let pdata = match &kinls_mem.prec_module {
        PrecModule::BBDPre(p) => p,
        _ => {
            KINProcessError(None, KINLS_PMEM_NULL, line!(), "KINBBDPrecGetNumGfnEvals", file!(),
                            MSGBBD_PMEM_NULL);
            return KINLS_PMEM_NULL;
        }
    };

    *ngevalsBBDP = pdata.nge;

    KINLS_SUCCESS
}

/*------------------------------------------------------------------
  KINBBDPrecSetup

  KINBBDPrecSetup generates and factors a banded block of the
  preconditioner matrix on each processor, via calls to the
  user-supplied gloc and gcomm functions. It uses difference
  quotient approximations to the Jacobian elements.

  KINBBDPrecSetup calculates a new Jacobian, stored in banded
  matrix PP and does an LU factorization of P in place in PP.

  The C arguments map to KINMem fields: uu = kin_uu (the current
  value of the dependent variable vector), uscale = kin_uscale, and
  fval/fscale are unused (as passed by kinLsPSetup); bbd_data is
  the PrecModule::BBDPre payload.

  Note: The value to be returned by the KINBBDPrecSetup function
  is a flag indicating whether it was successful. This value is:
    0 if successful,
    > 0 for a recoverable error - step will be retried.
  ------------------------------------------------------------------*/
pub fn KINBBDPrecSetup(kin_mem: &mut KINMem, pdata: &mut KBBDPrecData) -> i32 {
    /* Call KBBDDQJac for a new Jacobian calculation and store in PP */
    let retval = SUNMatZero(&mut pdata.PP);
    if retval != 0 {
        KINProcessError(Some(kin_mem), -1, line!(), "KINBBDPrecSetup", file!(),
                        MSGBBD_SUNMAT_FAIL);
        return -1;
    }

    /* gu = pdata.tempv1, gtemp = pdata.tempv2, utemp = pdata.tempv3 */
    let retval = KBBDDQJac(kin_mem, pdata);
    if retval != 0 {
        KINProcessError(Some(kin_mem), -1, line!(), "KINBBDPrecSetup", file!(),
                        MSGBBD_FUNC_FAILED);
        return -1;
    }

    /* Do LU factorization of P and return error flag */
    let KBBDPrecData { LS, PP, .. } = pdata;
    LS.setup(Some(PP))
}

/*------------------------------------------------------------------
  KINBBDPrecSolve

  KINBBDPrecSolve solves a linear system P z = r, with the
  banded blocked preconditioner matrix P generated and factored
  by KINBBDPrecSetup. Here, r comes in as vv and z is
  returned in vv as well.

  The C arguments uu/uscale/fval/fscale are all unused; vv is the
  vector initially set to the right-hand side vector r, but which
  upon return contains a solution of the linear system P*z = r.

  Note: The value returned by the KINBBDPrecSolve function is a
  flag returned from the linear solver object.
  ------------------------------------------------------------------*/
pub fn KINBBDPrecSolve(_kin_mem: &mut KINMem, pdata: &mut KBBDPrecData, vv: &mut NVector) -> i32 {
    /* Attach local data array for vv to rlocal: the C code aliases the
       data pointers; here the serial single-block reduction copies. */
    N_VScale(ONE, vv, &mut pdata.rlocal);

    /* Call banded solver object to do the work */
    let retval = {
        let KBBDPrecData { LS, PP, zlocal, rlocal, .. } = pdata;
        let mut atimes = |_v: &NVector, _jv: &mut NVector| -> i32 { 0 };
        LS.solve(Some(PP), zlocal, rlocal, ZERO, &mut atimes, None, None, None)
    };

    /* Copy result into vv */
    for i in 0..(pdata.n_local as usize) {
        vv.data[i] = pdata.zlocal.data[i];
    }

    retval
}

/* KINBBDPrecFree: the C routine frees LS, zlocal, rlocal, tempv1-3
   and PP; here dropping the PrecModule::BBDPre value releases
   everything. */

/*------------------------------------------------------------------
  KBBDDQJac

  This routine generates a banded difference quotient
  approximation to the Jacobian of f(u). It assumes that a band
  matrix of type SUNMatrix is stored column-wise, and that elements
  within each column are contiguous. All matrix elements are
  generated as difference quotients, by way of calls to the user
  routine gloc. By virtue of the band structure, the number of
  these calls is bandwidth + 1, where bandwidth = ml + mu + 1.
  This routine also assumes that the local elements of a vector
  are stored contiguously.

  C arguments map as: uu = kin_mem.kin_uu, uscale = kin_mem.kin_uscale,
  gu = pdata.tempv1, gtemp = pdata.tempv2, utemp = pdata.tempv3
  (as passed by KINBBDPrecSetup).
  ------------------------------------------------------------------*/
fn KBBDDQJac(kin_mem: &mut KINMem, pdata: &mut KBBDPrecData) -> i32 {
    let gloc = pdata.gloc.unwrap();

    /* load utemp with uu = predicted solution vector */
    N_VScale(ONE, &kin_mem.kin_uu, &mut pdata.tempv3);

    /* Call gcomm and gloc to get base value of g(uu) */
    if let Some(gcomm) = pdata.gcomm {
        let retval = gcomm(pdata.n_local, &kin_mem.kin_uu, &mut kin_mem.kin_user_data);
        if retval != 0 {
            return retval;
        }
    }

    let retval = {
        let KINMem { kin_uu, kin_user_data, .. } = kin_mem;
        gloc(pdata.n_local, kin_uu, &mut pdata.tempv1, kin_user_data)
    };
    pdata.nge += 1;
    if retval != 0 {
        return retval;
    }

    /* Set bandwidth and number of column groups for band differencing */
    let width = pdata.mldq + pdata.mudq + 1;
    let ngroups = i64::min(width, pdata.n_local);

    /* Loop over groups */
    for group in 1..=ngroups {
        /* increment all u_j in group */
        let mut j = group - 1;
        while j < pdata.n_local {
            let ju = j as usize;
            let inc = pdata.rel_uu
                * SUNMAX(
                    SUNRabs(kin_mem.kin_uu.data[ju]),
                    ONE / kin_mem.kin_uscale.data[ju],
                );
            pdata.tempv3.data[ju] += inc;
            j += width;
        }

        /* Evaluate g with incremented u */
        let retval = gloc(
            pdata.n_local,
            &pdata.tempv3,
            &mut pdata.tempv2,
            &mut kin_mem.kin_user_data,
        );
        pdata.nge += 1;
        if retval != 0 {
            return retval;
        }

        /* restore utemp, then form and load difference quotients */
        let pmat = match &mut pdata.PP {
            SUNMatrix::Band(m) => m,
            /* unreachable: PP is band by construction */
            _ => return -1,
        };
        let mut j = group - 1;
        while j < pdata.n_local {
            let ju = j as usize;
            pdata.tempv3.data[ju] = kin_mem.kin_uu.data[ju];
            let inc = pdata.rel_uu
                * SUNMAX(
                    SUNRabs(kin_mem.kin_uu.data[ju]),
                    ONE / kin_mem.kin_uscale.data[ju],
                );
            let inc_inv = ONE / inc;
            let i1 = i64::max(0, j - pdata.mukeep);
            let i2 = i64::min(j + pdata.mlkeep, pdata.n_local - 1);
            for i in i1..=i2 {
                /* SM_COLUMN_ELEMENT_B(col_j, i, j) */
                let val =
                    inc_inv * (pdata.tempv2.data[i as usize] - pdata.tempv1.data[i as usize]);
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
    use crate::kinsol_ls::{kinLsPSetup, KINSetLinearSolver};
    use crate::kinsol_ls_impl::KINLS_SUCCESS;
    use crate::nvector_serial::NVector;
    use crate::sundials_linearsolver::SUN_PREC_RIGHT;
    use crate::sundials_types::UserData;
    use crate::sunlinsol_spgmr::SUNLinSol_SPGMR;

    /* analytic tridiagonal gloc: g_i = 4 u_i - u_{i-1} - u_{i+1},
       so J = tridiag(-1, 4, -1) exactly (linear => DQ is exact up
       to roundoff in the increments) */
    fn gloc_tridiag(Nlocal: i64, uu: &NVector, gval: &mut NVector, _ud: &mut UserData) -> i32 {
        let n = Nlocal as usize;
        for i in 0..n {
            let mut v = 4.0 * uu.data[i];
            if i > 0 {
                v -= uu.data[i - 1];
            }
            if i + 1 < n {
                v -= uu.data[i + 1];
            }
            gval.data[i] = v;
        }
        0
    }

    fn gcomm_count(_Nlocal: i64, _u: &NVector, user_data: &mut UserData) -> i32 {
        if let Some(d) = user_data {
            if let Some(c) = d.downcast_mut::<i64>() {
                *c += 1;
            }
        }
        0
    }

    fn make_kin_mem(n: usize) -> KINMem {
        let mut kin_mem = KINMem::default();
        kin_mem.kin_uu = NVector::from_slice(&(1..=n).map(|i| 0.5 + i as f64).collect::<Vec<_>>());
        kin_mem.kin_fval = NVector::new(n);
        kin_mem.kin_uscale = NVector::from_slice(&vec![1.0; n]);
        kin_mem.kin_fscale = NVector::from_slice(&vec![1.0; n]);
        kin_mem.kin_vtemp1 = NVector::new(n);
        kin_mem.kin_vtemp2 = NVector::new(n);
        kin_mem
    }

    /* KINBBDPrecInit requires the KINLS interface (KINLS_LMEM_NULL) */
    #[test]
    fn kinbbdprecinit_needs_lmem() {
        let mut kin_mem = make_kin_mem(4);
        let flag = KINBBDPrecInit(&mut kin_mem, 4, 1, 1, 1, 1, 0.0,
                                  Some(gloc_tridiag), None);
        assert_eq!(flag, KINLS_LMEM_NULL);
    }

    /* banded DQ Jacobian against the analytic tridiagonal gloc, then
       the band LU solve must round-trip: solve P z = A x and get x */
    #[test]
    fn kinbbdprec_dq_jacobian_and_solve_roundtrip() {
        const N: usize = 5;
        let mut kin_mem = make_kin_mem(N);
        let sunctx = kin_mem.kin_sunctx.clone();

        /* iterative LS (SPGMR) so the preconditioner is legal */
        let ls = SUNLinSol_SPGMR(&kin_mem.kin_vtemp1, SUN_PREC_RIGHT, 0, &sunctx);
        assert_eq!(KINSetLinearSolver(&mut kin_mem, ls, None), KINLS_SUCCESS);

        /* count gcomm calls through user_data */
        kin_mem.kin_user_data = Some(Box::new(0i64));

        let flag = KINBBDPrecInit(&mut kin_mem, N as i64, 1, 1, 1, 1, 0.0,
                                  Some(gloc_tridiag), Some(gcomm_count));
        assert_eq!(flag, KINLS_SUCCESS);

        /* run the setup through the KINLS dispatch (kinLsPSetup) */
        let mut lmem = std::mem::take(&mut kin_mem.kin_lmem);
        let kinls_mem = match &mut lmem {
            LsModule::Ls(m) => m,
            _ => unreachable!(),
        };
        assert_eq!(kinLsPSetup(&mut kin_mem, kinls_mem), 0);
        assert_eq!(kinls_mem.npe, 1);

        let pdata = match &mut kinls_mem.prec_module {
            PrecModule::BBDPre(p) => p,
            _ => unreachable!(),
        };

        /* base g eval + one eval per column group (width = mldq+mudq+1 = 3) */
        assert_eq!(pdata.nge, 1 + 3);

        /* band LU solve round-trip: r = A*x with A = tridiag(-1,4,-1) */
        let x = [1.0, -2.0, 3.0, -4.0, 5.0];
        let mut vv = NVector::new(N);
        for i in 0..N {
            let mut v = 4.0 * x[i];
            if i > 0 {
                v -= x[i - 1];
            }
            if i + 1 < N {
                v -= x[i + 1];
            }
            vv.data[i] = v;
        }
        assert_eq!(KINBBDPrecSolve(&mut kin_mem, pdata, &mut vv), 0);
        for i in 0..N {
            assert!(
                (vv.data[i] - x[i]).abs() < 1.0e-6,
                "z[{}] = {} != {}",
                i,
                vv.data[i],
                x[i]
            );
        }

        kin_mem.kin_lmem = lmem;

        /* gcomm was called once, before the base gloc evaluation */
        let ncomm = *kin_mem
            .kin_user_data
            .as_ref()
            .unwrap()
            .downcast_ref::<i64>()
            .unwrap();
        assert_eq!(ncomm, 1);

        /* optional output getters */
        let mut nge = -1;
        assert_eq!(KINBBDPrecGetNumGfnEvals(&mut kin_mem, &mut nge), KINLS_SUCCESS);
        assert_eq!(nge, 4);
        let (mut lenrw, mut leniw) = (-1, -1);
        assert_eq!(KINBBDPrecGetWorkSpace(&mut kin_mem, &mut lenrw, &mut leniw), KINLS_SUCCESS);
        assert!(lenrw > 0);
    }
}
