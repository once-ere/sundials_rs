/* -----------------------------------------------------------------
 * Translated from src/arkode/arkode_bbdpre.c (ARKODE 7.7.0).
 * Band-block-diagonal preconditioner, i.e. a block-diagonal matrix
 * with banded blocks, for use with the ARKLS linear solver interface
 * and (in C) the MPI-parallel implementation of NVECTOR.  This
 * pure-Rust build is the serial reduction: a single block of
 * dimension n_local.
 *
 * In C the module stores its data behind arkls_mem->P_data and
 * installs ARKBBDPrecSetup/ARKBBDPrecSolve via
 * ARKodeSetPreconditioner; here ARKBBDPrecInit installs
 * PrecModule::BBDPre(Box<ARKBBDPrecData>) in the ARKLS memory and
 * arkLsSetup / the psolve closure in arkode_ls.rs dispatch to
 * ARKBBDPrecSetup / ARKBBDPrecSolve.
 * -----------------------------------------------------------------*/
use crate::arkode::{arkAllocVec, arkFreeVec};
use crate::arkode_bbdpre_impl::{
    ARKBBDPrecData, ARKCommFn, ARKLocalFn, MSG_BBD_FUNC_FAILED, MSG_BBD_PMEM_NULL,
    MSG_BBD_SUNLS_FAIL, MSG_BBD_SUNMAT_FAIL,
};
use crate::arkode_impl::{arkProcessError, ARKodeMem};
use crate::arkode_ls::{arkLs_AccessLMem, ark_rwt};
use crate::arkode_ls_impl::{
    PrecModule, ARKLS_ILL_INPUT, ARKLS_PMEM_NULL, ARKLS_SUCCESS, ARKLS_SUNLS_FAIL,
};
use crate::nvector_serial::{N_VScale, N_VSpace, N_VWrmsNorm, NVector};
use crate::sundials_context::SUNContext_Create;
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

/*---------------------------------------------------------------
 User-Callable Functions: initialization, reinit and free
---------------------------------------------------------------*/
#[allow(clippy::too_many_arguments)]
pub fn ARKBBDPrecInit(
    ark_mem: &mut ARKodeMem,
    Nlocal: i64,
    mudq: i64,
    mldq: i64,
    mukeep: i64,
    mlkeep: i64,
    dqrely: f64,
    gloc: Option<ARKLocalFn>,
    cfn: Option<ARKCommFn>,
) -> i32 {
    /* access ARKodeMem and ARKLsMem structure */
    let mut arkls_mem = match arkLs_AccessLMem(ark_mem, "ARKBBDPrecInit") {
        Ok(l) => l,
        Err(e) => return e,
    };

    /* Test compatibility of NVECTOR package with the BBD preconditioner:
       N_VGetArrayPointer/N_VSetArrayPointer always exist for the serial
       NVector. */

    /* Set pointers to gloc and cfn; load half-bandwidths */
    let mudq_c = i64::min(Nlocal - 1, i64::max(0, mudq));
    let mldq_c = i64::min(Nlocal - 1, i64::max(0, mldq));
    let muk = i64::min(Nlocal - 1, i64::max(0, mukeep));
    let mlk = i64::min(Nlocal - 1, i64::max(0, mlkeep));

    /* Allocate memory for saved Jacobian */
    let savedJ = SUNBandMatrixStorage(Nlocal, muk, mlk, muk, &SUNContext_Create());

    /* Allocate memory for preconditioner matrix */
    let storage_mu = i64::min(Nlocal - 1, muk + mlk);
    let savedP = SUNBandMatrixStorage(Nlocal, muk, mlk, storage_mu, &SUNContext_Create());

    /* Allocate memory for temporary N_Vectors: in C zlocal/rlocal are
       empty serial wrappers that borrow the data of z/r during the
       solve; here they own their storage and ARKBBDPrecSolve copies. */
    let zlocal = NVector::new(Nlocal as usize);
    let rlocal = NVector::new(Nlocal as usize);
    let tmpl_len = ark_mem.tempv1.len();
    let mut tmp1 = NVector::new(0);
    arkAllocVec(ark_mem, tmpl_len, &mut tmp1);
    let mut tmp2 = NVector::new(0);
    arkAllocVec(ark_mem, tmpl_len, &mut tmp2);
    let mut tmp3 = NVector::new(0);
    arkAllocVec(ark_mem, tmpl_len, &mut tmp3);

    /* Allocate memory for banded linear solver */
    let mut LS = SUNLinSol_Band(&rlocal, &savedP, &SUNContext_Create());

    /* initialize band linear solver object */
    let retval = LS.initialize();
    if retval != SUN_SUCCESS {
        arkFreeVec(ark_mem, &mut tmp1);
        arkFreeVec(ark_mem, &mut tmp2);
        arkFreeVec(ark_mem, &mut tmp3);
        arkProcessError(
            Some(ark_mem),
            ARKLS_SUNLS_FAIL,
            line!(),
            "ARKBBDPrecInit",
            file!(),
            MSG_BBD_SUNLS_FAIL,
        );
        ark_mem.lmem = Some(arkls_mem);
        return ARKLS_SUNLS_FAIL;
    }

    /* Set dqrely based on input dqrely (0 implies default). */
    let dqrely = if dqrely > ZERO {
        dqrely
    } else {
        SUNRsqrt(ark_mem.uround)
    };

    /* Set work space sizes and initialize nge */
    let mut rpwsize: i64 = 0;
    let mut ipwsize: i64 = 0;
    {
        let (lrw1, liw1) = N_VSpace(&ark_mem.tempv1);
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

    let pdata = Box::new(ARKBBDPrecData {
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
        /* Store Nlocal to be used in ARKBBDPrecSetup */
        n_local: Nlocal,
        rpwsize,
        ipwsize,
        nge: 0,
    });

    /* make sure P_data is free from any previous allocations */
    match std::mem::replace(&mut arkls_mem.prec_module, PrecModule::None) {
        PrecModule::BandPre(old) => crate::arkode_bandpre::ARKBandPrecFree(ark_mem, old),
        PrecModule::BBDPre(old) => ARKBBDPrecFree(ark_mem, old),
        PrecModule::None => {}
    }

    /* Point to the new P_data field in the LS memory (the pfree
       function is the PrecModule drop path) and attach the
       preconditioner solve and setup functions: in C this is
       ARKodeSetPreconditioner(arkode_mem, ARKBBDPrecSetup,
       ARKBBDPrecSolve); here the prec_module drives the dispatch in
       arkLsSetup / arkLsSolveIterative. */
    arkls_mem.pset = None;
    arkls_mem.psolve = None;
    arkls_mem.prec_module = PrecModule::BBDPre(pdata);

    /* ARKodeSetPreconditioner issues an error if the LS object does
       not allow user-supplied preconditioning */
    if !arkls_mem.iterative {
        arkProcessError(
            Some(ark_mem),
            ARKLS_ILL_INPUT,
            line!(),
            "ARKBBDPrecInit",
            file!(),
            "SUNLinearSolver object does not support user-supplied preconditioning",
        );
        ark_mem.lmem = Some(arkls_mem);
        return ARKLS_ILL_INPUT;
    }

    ark_mem.lmem = Some(arkls_mem);
    ARKLS_SUCCESS
}

/*-------------------------------------------------------------*/
pub fn ARKBBDPrecReInit(ark_mem: &mut ARKodeMem, mudq: i64, mldq: i64, dqrely: f64) -> i32 {
    /* access ARKodeMem and ARKLsMem structure */
    let mut arkls_mem = match arkLs_AccessLMem(ark_mem, "ARKBBDPrecReInit") {
        Ok(l) => l,
        Err(e) => return e,
    };

    /* Return immediately ARKBBDPrecData is NULL */
    let pdata = match &mut arkls_mem.prec_module {
        PrecModule::BBDPre(p) => p,
        _ => {
            arkProcessError(
                Some(ark_mem),
                ARKLS_PMEM_NULL,
                line!(),
                "ARKBBDPrecReInit",
                file!(),
                MSG_BBD_PMEM_NULL,
            );
            ark_mem.lmem = Some(arkls_mem);
            return ARKLS_PMEM_NULL;
        }
    };

    /* Load half-bandwidths */
    let Nlocal = pdata.n_local;
    pdata.mudq = i64::min(Nlocal - 1, i64::max(0, mudq));
    pdata.mldq = i64::min(Nlocal - 1, i64::max(0, mldq));

    /* Set dqrely based on input dqrely (0 implies default). */
    pdata.dqrely = if dqrely > ZERO {
        dqrely
    } else {
        SUNRsqrt(ark_mem.uround)
    };

    /* Re-initialize nge */
    pdata.nge = 0;

    ark_mem.lmem = Some(arkls_mem);
    ARKLS_SUCCESS
}

/*-------------------------------------------------------------*/
pub fn ARKBBDPrecGetWorkSpace(
    ark_mem: &mut ARKodeMem,
    lenrwBBDP: &mut i64,
    leniwBBDP: &mut i64,
) -> i32 {
    /* access ARKodeMem and ARKLsMem structure */
    let arkls_mem = match arkLs_AccessLMem(ark_mem, "ARKBBDPrecGetWorkSpace") {
        Ok(l) => l,
        Err(e) => return e,
    };

    /* Return immediately ARKBBDPrecData is NULL */
    let pdata = match &arkls_mem.prec_module {
        PrecModule::BBDPre(p) => p,
        _ => {
            arkProcessError(
                Some(ark_mem),
                ARKLS_PMEM_NULL,
                line!(),
                "ARKBBDPrecGetWorkSpace",
                file!(),
                MSG_BBD_PMEM_NULL,
            );
            ark_mem.lmem = Some(arkls_mem);
            return ARKLS_PMEM_NULL;
        }
    };

    /* set outputs */
    *lenrwBBDP = pdata.rpwsize;
    *leniwBBDP = pdata.ipwsize;

    ark_mem.lmem = Some(arkls_mem);
    ARKLS_SUCCESS
}

/*-------------------------------------------------------------*/
pub fn ARKBBDPrecGetNumGfnEvals(ark_mem: &mut ARKodeMem, ngevalsBBDP: &mut i64) -> i32 {
    /* access ARKodeMem and ARKLsMem structure */
    let arkls_mem = match arkLs_AccessLMem(ark_mem, "ARKBBDPrecGetNumGfnEvals") {
        Ok(l) => l,
        Err(e) => return e,
    };

    /* Return immediately if ARKBBDPrecData is NULL */
    let pdata = match &arkls_mem.prec_module {
        PrecModule::BBDPre(p) => p,
        _ => {
            arkProcessError(
                Some(ark_mem),
                ARKLS_PMEM_NULL,
                line!(),
                "ARKBBDPrecGetNumGfnEvals",
                file!(),
                MSG_BBD_PMEM_NULL,
            );
            ark_mem.lmem = Some(arkls_mem);
            return ARKLS_PMEM_NULL;
        }
    };

    /* set output */
    *ngevalsBBDP = pdata.nge;

    ark_mem.lmem = Some(arkls_mem);
    ARKLS_SUCCESS
}

/*---------------------------------------------------------------
 ARKBBDPrecSetup:

 ARKBBDPrecSetup generates and factors a banded block of the
 preconditioner matrix on each processor, via calls to the
 user-supplied gloc and cfn functions. It uses difference
 quotient approximations to the Jacobian elements.

 ARKBBDPrecSetup calculates a new J, if necessary, then
 calculates P = M - gamma*J, and does an LU factorization of P.

 Return value:
   0  if successful,
   1  for a recoverable error (step will be retried).
---------------------------------------------------------------*/
#[allow(clippy::too_many_arguments)]
pub fn ARKBBDPrecSetup(
    ark_mem: &mut ARKodeMem,
    pdata: &mut ARKBBDPrecData,
    t: f64,
    y: &NVector,
    _fy: &NVector,
    jok: bool,
    jcurPtr: &mut bool,
    gamma: f64,
) -> i32 {
    /* If jok = SUNTRUE, use saved copy of J */
    if jok {
        *jcurPtr = SUNFALSE;
        let retval = SUNMatCopy(&pdata.savedJ, &mut pdata.savedP);
        if retval < 0 {
            arkProcessError(
                Some(ark_mem),
                -1,
                line!(),
                "ARKBBDPrecSetup",
                file!(),
                MSG_BBD_SUNMAT_FAIL,
            );
            return -1;
        }
        if retval > 0 {
            return 1;
        }

        /* Otherwise call ARKBBDDQJac for new J value */
    } else {
        *jcurPtr = SUNTRUE;
        let retval = SUNMatZero(&mut pdata.savedJ);
        if retval < 0 {
            arkProcessError(
                Some(ark_mem),
                -1,
                line!(),
                "ARKBBDPrecSetup",
                file!(),
                MSG_BBD_SUNMAT_FAIL,
            );
            return -1;
        }
        if retval > 0 {
            return 1;
        }

        /* (gy = pdata.tmp1, ytemp = pdata.tmp2, gtemp = pdata.tmp3) */
        let retval = ARKBBDDQJac(ark_mem, pdata, t, y);
        if retval < 0 {
            arkProcessError(
                Some(ark_mem),
                -1,
                line!(),
                "ARKBBDPrecSetup",
                file!(),
                MSG_BBD_FUNC_FAILED,
            );
            return -1;
        }
        if retval > 0 {
            return 1;
        }

        let retval = SUNMatCopy(&pdata.savedJ, &mut pdata.savedP);
        if retval < 0 {
            arkProcessError(
                Some(ark_mem),
                -1,
                line!(),
                "ARKBBDPrecSetup",
                file!(),
                MSG_BBD_SUNMAT_FAIL,
            );
            return -1;
        }
        if retval > 0 {
            return 1;
        }
    }

    /* Scale and add I to get P = I - gamma*J */
    let retval = SUNMatScaleAddI(-gamma, &mut pdata.savedP);
    if retval != 0 {
        arkProcessError(
            Some(ark_mem),
            -1,
            line!(),
            "ARKBBDPrecSetup",
            file!(),
            MSG_BBD_SUNMAT_FAIL,
        );
        return -1;
    }

    /* Do LU factorization of matrix and return error flag */
    let ARKBBDPrecData { LS, savedP, .. } = pdata;
    LS.setup(Some(savedP))
}

/*---------------------------------------------------------------
 ARKBBDPrecSolve:

 ARKBBDPrecSolve solves a linear system P z = r, with the
 band-block-diagonal preconditioner matrix P generated and
 factored by ARKBBDPrecSetup.

 The value returned by the ARKBBDPrecSolve function is the same
 as the value returned from the linear solver object.
---------------------------------------------------------------*/
pub fn ARKBBDPrecSolve(pdata: &mut ARKBBDPrecData, r: &NVector, z: &mut NVector) -> i32 {
    /* Attach local data arrays for r and z to rlocal and zlocal:
       the C code aliases the data pointers; here the serial
       single-block reduction copies in and out instead. */
    N_VScale(ONE, r, &mut pdata.rlocal);

    /* Call banded solver object to do the work */
    let retval = {
        let ARKBBDPrecData {
            LS,
            savedP,
            zlocal,
            rlocal,
            ..
        } = pdata;
        let mut atimes = |_v: &NVector, _jv: &mut NVector| -> i32 { 0 };
        LS.solve(Some(savedP), zlocal, rlocal, ZERO, &mut atimes, None, None, None)
    };

    /* Detach local data arrays from rlocal and zlocal */
    N_VScale(ONE, &pdata.zlocal, z);

    retval
}

/*---------------------------------------------------------------
 ARKBBDPrecFree:

 Frees data associated with the ARKBBD preconditioner.  (The C
 routine also frees LS/zlocal/rlocal/savedP/savedJ; here dropping
 the pdata box releases them, while arkFreeVec keeps the lrw/liw
 accounting for tmp1-3.)
---------------------------------------------------------------*/
pub(crate) fn ARKBBDPrecFree(ark_mem: &mut ARKodeMem, mut pdata: Box<ARKBBDPrecData>) {
    arkFreeVec(ark_mem, &mut pdata.tmp1);
    arkFreeVec(ark_mem, &mut pdata.tmp2);
    arkFreeVec(ark_mem, &mut pdata.tmp3);
}

/*---------------------------------------------------------------
 ARKBBDDQJac:

 This routine generates a banded difference quotient approximation
 to the local block of the Jacobian of g(t,y). It assumes that a
 band matrix of type SUNMatrix is stored columnwise, and that
 elements within each column are contiguous. All matrix elements
 are generated as difference quotients, by way of calls to the
 user routine gloc.  By virtue of the band structure, the number
 of these calls is bandwidth + 1, where bandwidth = mldq + mudq + 1.
 But the band matrix kept has bandwidth = mlkeep + mukeep + 1.

 C arguments map as: gy = pdata.tmp1, ytemp = pdata.tmp2,
 gtemp = pdata.tmp3.
---------------------------------------------------------------*/
fn ARKBBDDQJac(ark_mem: &mut ARKodeMem, pdata: &mut ARKBBDPrecData, t: f64, y: &NVector) -> i32 {
    let gloc = pdata.gloc.unwrap();

    /* Load ytemp with y = predicted solution vector */
    N_VScale(ONE, y, &mut pdata.tmp2);

    /* Call cfn and gloc to get base value of g(t,y) */
    if let Some(cfn) = pdata.cfn {
        let retval = cfn(pdata.n_local, t, y, &mut ark_mem.user_data);
        if retval != 0 {
            return retval;
        }
    }

    let retval = gloc(
        pdata.n_local,
        t,
        &pdata.tmp2,
        &mut pdata.tmp1,
        &mut ark_mem.user_data,
    );
    pdata.nge += 1;
    if retval != 0 {
        return retval;
    }

    /* Set minimum increment based on uround and norm of g */
    let gnorm = N_VWrmsNorm(&pdata.tmp1, ark_rwt(ark_mem));
    let minInc = if gnorm != ZERO {
        MIN_INC_MULT * SUNRabs(ark_mem.h) * ark_mem.uround * pdata.n_local as f64 * gnorm
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
                pdata.dqrely * SUNRabs(y.data[ju]),
                minInc / ark_mem.ewt.data[ju],
            );
            let yj = y.data[ju];

            /* Adjust sign(inc) again if yj has an inequality constraint. */
            if let Some(constraints) = &ark_mem.constraints {
                let conj = constraints.data[ju];
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
            &mut ark_mem.user_data,
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
            let yj = y.data[ju];
            pdata.tmp2.data[ju] = y.data[ju];
            let mut inc = SUNMAX(
                pdata.dqrely * SUNRabs(y.data[ju]),
                minInc / ark_mem.ewt.data[ju],
            );

            if let Some(constraints) = &ark_mem.constraints {
                let conj = constraints.data[ju];
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
                let val = inc_inv * (pdata.tmp3.data[i as usize] - pdata.tmp1.data[i as usize]);
                jmat.set(i, j, val);
            }
            j += width;
        }
    }

    0
}

/* -----------------------------------------------------------------
   Unit tests: the C ARKBBDPRE drivers are MPI-parallel (excluded
   backend), so the serial reduction is validated directly — the same
   implicit 1D heat equation as the ARKBANDPRE tests, solved with
   ARKStep + SPGMR + ARKBBDPrecInit (gloc = the full RHS, cfn = None),
   checked against the exact semi-discrete solution and against a
   run preconditioned with ARKBANDPRE (identical bandwidths and the
   default dqrely = sqrt(uround) make the two DQ Jacobians agree).
   -----------------------------------------------------------------*/
#[cfg(test)]
mod tests {
    use super::*;
    use crate::arkode::{ARKodeEvolve, ARKodeFree, ARKodeSStolerances};
    use crate::arkode_arkstep::ARKStepCreate;
    use crate::arkode_bandpre::ARKBandPrecInit;
    use crate::arkode_impl::{ARKodeMem, ARK_NORMAL};
    use crate::arkode_ls::ARKodeSetLinearSolver;
    use crate::nvector_serial::N_VNew_Serial;
    use crate::sundials_context::SUNContext;
    use crate::sundials_linearsolver::SUN_PREC_LEFT;
    use crate::sundials_types::UserData;
    use crate::sunlinsol_spgmr::SUNLinSol_SPGMR;
    use std::f64::consts::PI;

    const NEQ: i64 = 32;

    /* 1D heat equation, second-order centered differences, zero
       Dirichlet boundaries: udot_j = (u_{j-1} - 2 u_j + u_{j+1})/dx^2 */
    fn heat_rhs(_t: f64, y: &NVector, ydot: &mut NVector, _user_data: &mut UserData) -> i32 {
        let n = y.data.len();
        let dx = 1.0 / (n as f64 + 1.0);
        let coef = 1.0 / (dx * dx);
        for j in 0..n {
            let ul = if j == 0 { 0.0 } else { y.data[j - 1] };
            let ur = if j == n - 1 { 0.0 } else { y.data[j + 1] };
            ydot.data[j] = coef * (ul - 2.0 * y.data[j] + ur);
        }
        0
    }

    /* gloc: the local approximation g(t,y) is the full RHS here */
    fn heat_gloc(_nlocal: i64, t: f64, y: &NVector, g: &mut NVector, user_data: &mut UserData) -> i32 {
        heat_rhs(t, y, g, user_data)
    }

    /* build the integrator with SPGMR attached and y0 = sin(pi x) */
    fn setup(ctx: &SUNContext) -> (Box<ARKodeMem>, NVector) {
        let dx = 1.0 / (NEQ as f64 + 1.0);
        let mut y = N_VNew_Serial(NEQ, ctx);
        for j in 0..NEQ as usize {
            y.data[j] = (PI * dx * (j as f64 + 1.0)).sin();
        }
        let mut arkode_mem = ARKStepCreate(None, Some(heat_rhs), 0.0, &y, ctx).expect("ARKStepCreate");
        assert_eq!(ARKodeSStolerances(&mut arkode_mem, 1.0e-6, 1.0e-10), 0);
        let ls = SUNLinSol_SPGMR(&y, SUN_PREC_LEFT, 0, ctx);
        assert_eq!(ARKodeSetLinearSolver(&mut arkode_mem, ls, None), ARKLS_SUCCESS);
        (arkode_mem, y)
    }

    #[test]
    fn bbdpre_heat1d_spgmr() {
        let ctx = SUNContext_Create();
        let dx = 1.0 / (NEQ as f64 + 1.0);

        /* BBD-preconditioned run */
        let (mut arkode_mem, mut y) = setup(&ctx);
        assert_eq!(
            ARKBBDPrecInit(&mut arkode_mem, NEQ, 1, 1, 1, 1, 0.0, Some(heat_gloc), None),
            ARKLS_SUCCESS
        );
        let mut t = 0.0;
        let retval = ARKodeEvolve(&mut arkode_mem, 0.1, &mut y, &mut t, ARK_NORMAL);
        assert!(retval >= 0, "ARKodeEvolve failed: {}", retval);

        /* solution accuracy vs the exact semi-discrete solution */
        let lambda = (2.0 * (PI * dx).cos() - 2.0) / (dx * dx);
        let decay = (lambda * t).exp();
        for j in 0..NEQ as usize {
            let exact = decay * (PI * dx * (j as f64 + 1.0)).sin();
            assert!(
                (y.data[j] - exact).abs() < 1.0e-4,
                "solution error at j={}: {} vs {}",
                j,
                y.data[j],
                exact
            );
        }

        /* gloc evaluations: each fresh DQ Jacobian costs 1 (base) +
           ngroups = mldq + mudq + 1 = 3 evaluations */
        let mut nge: i64 = 0;
        assert_eq!(ARKBBDPrecGetNumGfnEvals(&mut arkode_mem, &mut nge), ARKLS_SUCCESS);
        assert!(nge > 0 && nge % 4 == 0, "nge = {}", nge);

        /* workspace query */
        let mut lenrw: i64 = 0;
        let mut leniw: i64 = 0;
        assert_eq!(
            ARKBBDPrecGetWorkSpace(&mut arkode_mem, &mut lenrw, &mut leniw),
            ARKLS_SUCCESS
        );
        assert!(lenrw > 0 && leniw > 0, "lenrw = {}, leniw = {}", lenrw, leniw);

        /* the saved DQ Jacobian of the linear gloc must match the
           analytic tridiagonal Jacobian to difference-quotient
           rounding accuracy */
        let coef = 1.0 / (dx * dx);
        {
            let arkls_mem = arkode_mem.lmem.as_ref().unwrap();
            let bbd = match &arkls_mem.prec_module {
                PrecModule::BBDPre(p) => p,
                _ => panic!("prec_module is not BBDPre"),
            };
            let jm = match &bbd.savedJ {
                SUNMatrix::Band(m) => m,
                _ => panic!("savedJ is not band"),
            };
            for j in 0..NEQ {
                for i in i64::max(0, j - 1)..=i64::min(j + 1, NEQ - 1) {
                    let exact = if i == j { -2.0 * coef } else { coef };
                    assert!(
                        (jm.get(i, j) - exact).abs() < 1.0e-4 * coef,
                        "savedJ({},{}) = {} vs {}",
                        i,
                        j,
                        jm.get(i, j),
                        exact
                    );
                }
            }
        }

        /* ARKBBDPrecReInit reloads bandwidths/dqrely and resets nge */
        assert_eq!(ARKBBDPrecReInit(&mut arkode_mem, 1, 1, 0.0), ARKLS_SUCCESS);
        let mut nge2: i64 = -1;
        assert_eq!(ARKBBDPrecGetNumGfnEvals(&mut arkode_mem, &mut nge2), ARKLS_SUCCESS);
        assert_eq!(nge2, 0);

        /* BAND-preconditioned reference run: with identical bandwidths
           and the default dqrely = sqrt(uround) increments, both
           modules build the same preconditioner, so the trajectories
           agree far below the integration tolerances */
        let (mut arkode_mem2, mut y2) = setup(&ctx);
        assert_eq!(ARKBandPrecInit(&mut arkode_mem2, NEQ, 1, 1), ARKLS_SUCCESS);
        let mut t2 = 0.0;
        let retval = ARKodeEvolve(&mut arkode_mem2, 0.1, &mut y2, &mut t2, ARK_NORMAL);
        assert!(retval >= 0, "ARKodeEvolve (bandpre) failed: {}", retval);
        assert_eq!(t, t2);
        for j in 0..NEQ as usize {
            assert!(
                (y.data[j] - y2.data[j]).abs() < 1.0e-8,
                "bbd vs band mismatch at j={}: {} vs {}",
                j,
                y.data[j],
                y2.data[j]
            );
        }

        /* module mismatch: the ARKBANDPRE getters see BBD data (and
           vice versa) as ARKLS_PMEM_NULL */
        let mut dummy: i64 = 0;
        assert_eq!(
            crate::arkode_bandpre::ARKBandPrecGetNumRhsEvals(&mut arkode_mem, &mut dummy),
            ARKLS_PMEM_NULL
        );
        assert_eq!(
            ARKBBDPrecGetNumGfnEvals(&mut arkode_mem2, &mut dummy),
            ARKLS_PMEM_NULL
        );

        /* free path exercises ARKBBDPrecFree via arkLsFree */
        let mut slot = Some(arkode_mem);
        ARKodeFree(&mut slot);
        let mut slot2 = Some(arkode_mem2);
        ARKodeFree(&mut slot2);
    }
}
