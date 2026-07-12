/* -----------------------------------------------------------------
 * Translated from src/arkode/arkode_bandpre.c (ARKODE 7.7.0).
 * Banded difference quotient Jacobian-based preconditioner and
 * solver routines for use with the ARKLS linear solver interface.
 *
 * In C the module stores its data behind arkls_mem->P_data and
 * installs ARKBandPrecSetup/ARKBandPrecSolve via
 * ARKodeSetPreconditioner; here ARKBandPrecInit installs
 * PrecModule::BandPre(Box<ARKBandPrecData>) in the ARKLS memory and
 * arkLsSetup / the psolve closure in arkode_ls.rs dispatch to
 * ARKBandPrecSetup / ARKBandPrecSolve.
 * -----------------------------------------------------------------*/
use crate::arkode::{arkAllocVec, arkFreeVec};
use crate::arkode_bandpre_impl::{
    ARKBandPrecData, MSG_BP_PMEM_NULL, MSG_BP_RHSFUNC_FAILED, MSG_BP_SUNLS_FAIL,
    MSG_BP_SUNMAT_FAIL,
};
use crate::arkode_impl::{arkProcessError, ARKodeMem, ARK_PRERHSFN_FAIL};
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
 Initialization, Free, and Get Functions
 NOTE: The band linear solver assumes a serial implementation
       of the NVECTOR package (always true here, so the
       N_VGetArrayPointer compatibility test is a no-op).
---------------------------------------------------------------*/
pub fn ARKBandPrecInit(ark_mem: &mut ARKodeMem, N: i64, mu: i64, ml: i64) -> i32 {
    /* access ARKodeMem and ARKLsMem structures */
    let mut arkls_mem = match arkLs_AccessLMem(ark_mem, "ARKBandPrecInit") {
        Ok(l) => l,
        Err(e) => return e,
    };

    /* Test compatibility of NVECTOR package with the BAND preconditioner:
       N_VGetArrayPointer always exists for the serial NVector. */

    /* Load pointers and bandwidths into pdata block. */
    let mup = i64::min(N - 1, i64::max(0, mu));
    let mlp = i64::min(N - 1, i64::max(0, ml));

    /* Allocate memory for saved banded Jacobian approximation. */
    let savedJ = SUNBandMatrixStorage(N, mup, mlp, mup, &SUNContext_Create());

    /* Allocate memory for banded preconditioner. */
    let storagemu = i64::min(N - 1, mup + mlp);
    let savedP = SUNBandMatrixStorage(N, mup, mlp, storagemu, &SUNContext_Create());

    /* Allocate memory for banded linear solver */
    let mut LS = SUNLinSol_Band(&ark_mem.tempv1, &savedP, &SUNContext_Create());

    /* allocate memory for temporary N_Vectors */
    let tmpl_len = ark_mem.tempv1.len();
    let mut tmp1 = NVector::new(0);
    arkAllocVec(ark_mem, tmpl_len, &mut tmp1);
    let mut tmp2 = NVector::new(0);
    arkAllocVec(ark_mem, tmpl_len, &mut tmp2);

    /* initialize band linear solver object */
    let retval = LS.initialize();
    if retval != SUN_SUCCESS {
        arkFreeVec(ark_mem, &mut tmp1);
        arkFreeVec(ark_mem, &mut tmp2);
        arkProcessError(
            Some(ark_mem),
            ARKLS_SUNLS_FAIL,
            line!(),
            "ARKBandPrecInit",
            file!(),
            MSG_BP_SUNLS_FAIL,
        );
        ark_mem.lmem = Some(arkls_mem);
        return ARKLS_SUNLS_FAIL;
    }

    let pdata = Box::new(ARKBandPrecData {
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

    /* make sure s_P_data is free from any previous allocations */
    match std::mem::replace(&mut arkls_mem.prec_module, PrecModule::None) {
        PrecModule::BandPre(old) => ARKBandPrecFree(ark_mem, old),
        PrecModule::BBDPre(old) => crate::arkode_bbdpre::ARKBBDPrecFree(ark_mem, old),
        PrecModule::None => {}
    }

    /* Point to the new P_data field in the LS memory (the pfree
       function is the PrecModule drop path) and attach the
       preconditioner solve and setup functions: in C this is
       ARKodeSetPreconditioner(arkode_mem, ARKBandPrecSetup,
       ARKBandPrecSolve); here the prec_module drives the dispatch in
       arkLsSetup / arkLsSolveIterative. */
    arkls_mem.pset = None;
    arkls_mem.psolve = None;
    arkls_mem.prec_module = PrecModule::BandPre(pdata);

    /* ARKodeSetPreconditioner issues an error if the LS object does
       not allow user-supplied preconditioning */
    if !arkls_mem.iterative {
        arkProcessError(
            Some(ark_mem),
            ARKLS_ILL_INPUT,
            line!(),
            "ARKBandPrecInit",
            file!(),
            "SUNLinearSolver object does not support user-supplied preconditioning",
        );
        ark_mem.lmem = Some(arkls_mem);
        return ARKLS_ILL_INPUT;
    }

    ark_mem.lmem = Some(arkls_mem);
    ARKLS_SUCCESS
}

pub fn ARKBandPrecGetWorkSpace(
    ark_mem: &mut ARKodeMem,
    lenrwBP: &mut i64,
    leniwBP: &mut i64,
) -> i32 {
    /* access ARKodeMem and ARKLsMem structures */
    let arkls_mem = match arkLs_AccessLMem(ark_mem, "ARKBandPrecGetWorkSpace") {
        Ok(l) => l,
        Err(e) => return e,
    };

    /* Return immediately if ARKBandPrecData is NULL */
    let pdata = match &arkls_mem.prec_module {
        PrecModule::BandPre(bp) => bp,
        _ => {
            arkProcessError(
                Some(ark_mem),
                ARKLS_PMEM_NULL,
                line!(),
                "ARKBandPrecGetWorkSpace",
                file!(),
                MSG_BP_PMEM_NULL,
            );
            ark_mem.lmem = Some(arkls_mem);
            return ARKLS_PMEM_NULL;
        }
    };

    /* sum space requirements for all objects in pdata */
    *leniwBP = 4;
    *lenrwBP = 0;
    {
        let (lrw1, liw1) = N_VSpace(&ark_mem.tempv1);
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

    ark_mem.lmem = Some(arkls_mem);
    ARKLS_SUCCESS
}

pub fn ARKBandPrecGetNumRhsEvals(ark_mem: &mut ARKodeMem, nfevalsBP: &mut i64) -> i32 {
    /* access ARKodeMem and ARKLsMem structures */
    let arkls_mem = match arkLs_AccessLMem(ark_mem, "ARKBandPrecGetNumRhsEvals") {
        Ok(l) => l,
        Err(e) => return e,
    };

    /* Return immediately if ARKBandPrecData is NULL */
    let pdata = match &arkls_mem.prec_module {
        PrecModule::BandPre(bp) => bp,
        _ => {
            arkProcessError(
                Some(ark_mem),
                ARKLS_PMEM_NULL,
                line!(),
                "ARKBandPrecGetNumRhsEvals",
                file!(),
                MSG_BP_PMEM_NULL,
            );
            ark_mem.lmem = Some(arkls_mem);
            return ARKLS_PMEM_NULL;
        }
    };

    /* set output */
    *nfevalsBP = pdata.nfeBP;

    ark_mem.lmem = Some(arkls_mem);
    ARKLS_SUCCESS
}

/*---------------------------------------------------------------
 ARKBandPrecSetup:

 Together ARKBandPrecSetup and ARKBandPrecSolve use a banded
 difference quotient Jacobian to create a preconditioner.
 ARKBandPrecSetup calculates a new J, if necessary, then
 calculates P = I - gamma*J, and does an LU factorization of P.

 jok == SUNFALSE means recompute Jacobian-related data from
 scratch; jok == SUNTRUE means that Jacobian data from the
 previous PrecSetup call will be reused (with the current value
 of gamma).  *jcurPtr is set to SUNTRUE if Jacobian data was
 recomputed, SUNFALSE if saved data was reused.

 The value returned is
   0  if successful, or
   1  if the band factorization failed.
---------------------------------------------------------------*/
#[allow(clippy::too_many_arguments)]
pub fn ARKBandPrecSetup(
    ark_mem: &mut ARKodeMem,
    pdata: &mut ARKBandPrecData,
    t: f64,
    y: &NVector,
    fy: &NVector,
    jok: bool,
    jcurPtr: &mut bool,
    gamma: f64,
) -> i32 {
    /* Assume matrix and lpivots have already been allocated. */

    if jok {
        /* If jok = SUNTRUE, use saved copy of J. */
        *jcurPtr = SUNFALSE;
        let retval = SUNMatCopy(&pdata.savedJ, &mut pdata.savedP);
        if retval < 0 {
            arkProcessError(
                Some(ark_mem),
                -1,
                line!(),
                "ARKBandPrecSetup",
                file!(),
                MSG_BP_SUNMAT_FAIL,
            );
            return -1;
        }
        if retval > 0 {
            return 1;
        }
    } else {
        /* If jok = SUNFALSE, call ARKBandPDQJac for new J value. */
        *jcurPtr = SUNTRUE;
        let retval = SUNMatZero(&mut pdata.savedJ);
        if retval < 0 {
            arkProcessError(
                Some(ark_mem),
                -1,
                line!(),
                "ARKBandPrecSetup",
                file!(),
                MSG_BP_SUNMAT_FAIL,
            );
            return -1;
        }
        if retval > 0 {
            return 1;
        }

        /* (ftemp = pdata.tmp1, ytemp = pdata.tmp2) */
        let retval = ARKBandPDQJac(ark_mem, pdata, t, y, fy);
        if retval < 0 {
            arkProcessError(
                Some(ark_mem),
                -1,
                line!(),
                "ARKBandPrecSetup",
                file!(),
                MSG_BP_RHSFUNC_FAILED,
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
                "ARKBandPrecSetup",
                file!(),
                MSG_BP_SUNMAT_FAIL,
            );
            return -1;
        }
        if retval > 0 {
            return 1;
        }
    }

    /* Scale and add identity to get savedP = I - gamma*J. */
    let retval = SUNMatScaleAddI(-gamma, &mut pdata.savedP);
    if retval != 0 {
        arkProcessError(
            Some(ark_mem),
            -1,
            line!(),
            "ARKBandPrecSetup",
            file!(),
            MSG_BP_SUNMAT_FAIL,
        );
        return -1;
    }

    /* Do LU factorization of matrix and return error flag */
    let ARKBandPrecData { LS, savedP, .. } = pdata;
    LS.setup(Some(savedP))
}

/*---------------------------------------------------------------
 ARKBandPrecSolve:

 ARKBandPrecSolve solves a linear system P z = r, where P is the
 matrix computed by ARKBandPrecond.

 The value returned by the ARKBandPrecSolve function is always 0,
 indicating success.
---------------------------------------------------------------*/
pub fn ARKBandPrecSolve(pdata: &mut ARKBandPrecData, r: &NVector, z: &mut NVector) -> i32 {
    /* Assume matrix and linear solver have already been allocated. */

    /* Call banded solver object to do the work */
    let ARKBandPrecData { LS, savedP, .. } = pdata;
    let mut atimes = |_v: &NVector, _jv: &mut NVector| -> i32 { 0 };
    LS.solve(Some(savedP), z, r, ZERO, &mut atimes, None, None, None)
}

/*---------------------------------------------------------------
 ARKBandPrecFree:

 Frees data associated with the ARKBand preconditioner.  (The C
 routine also frees LS/savedP/savedJ; here dropping the pdata box
 releases them, while arkFreeVec keeps the lrw/liw accounting.)
---------------------------------------------------------------*/
pub(crate) fn ARKBandPrecFree(ark_mem: &mut ARKodeMem, mut pdata: Box<ARKBandPrecData>) {
    arkFreeVec(ark_mem, &mut pdata.tmp1);
    arkFreeVec(ark_mem, &mut pdata.tmp2);
}

/*---------------------------------------------------------------
 ARKBandPDQJac:

 This routine generates a banded difference quotient approximation to
 the Jacobian of f(t,y). It assumes that a band matrix of type
 SUNMatrix is stored column-wise, and that elements within each column
 are contiguous.

 C arguments map as: ftemp = pdata.tmp1, ytemp = pdata.tmp2.
---------------------------------------------------------------*/
fn ARKBandPDQJac(
    ark_mem: &mut ARKodeMem,
    pdata: &mut ARKBandPrecData,
    t: f64,
    y: &NVector,
    fy: &NVector,
) -> i32 {
    /* Access implicit RHS function */
    let fi = match ark_mem.step_getimplicitrhs {
        Some(get) => get(ark_mem),
        None => None,
    };
    let fi = match fi {
        Some(f) => f,
        None => return -1,
    };

    /* Load ytemp with y = predicted y vector. */
    N_VScale(ONE, y, &mut pdata.tmp2);

    /* Set minimum increment based on uround and norm of f. */
    let srur = SUNRsqrt(ark_mem.uround);
    let fnorm = N_VWrmsNorm(fy, ark_rwt(ark_mem));
    let minInc = if fnorm != ZERO {
        MIN_INC_MULT * SUNRabs(ark_mem.h) * ark_mem.uround * pdata.N as f64 * fnorm
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
            let mut inc = SUNMAX(srur * SUNRabs(y.data[ju]), minInc / ark_mem.ewt.data[ju]);
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

        /* call the user-supplied pre-RHS function (if supplied), then call RHS */
        if let Some(pre_rhs) = ark_mem.PreRhsFn {
            let retval = pre_rhs(t, &pdata.tmp2, &mut ark_mem.user_data);
            if retval != 0 {
                return ARK_PRERHSFN_FAIL;
            }
        }
        let retval = fi(t, &pdata.tmp2, &mut pdata.tmp1, &mut ark_mem.user_data);
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
            let yj = y.data[ju];
            pdata.tmp2.data[ju] = y.data[ju];
            let mut inc = SUNMAX(srur * SUNRabs(y.data[ju]), minInc / ark_mem.ewt.data[ju]);

            /* Adjust sign(inc) as before. */
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
            let i1 = i64::max(0, j - pdata.mu);
            let i2 = i64::min(j + pdata.ml, pdata.N - 1);
            for i in i1..=i2 {
                /* SM_COLUMN_ELEMENT_B(col_j, i, j) */
                let val = inc_inv * (pdata.tmp1.data[i as usize] - fy.data[i as usize]);
                jmat.set(i, j, val);
            }
            j += width;
        }
    }

    0
}

/* -----------------------------------------------------------------
   Unit tests: no serial C example exercises ARKBANDPRE (the C usage
   is documented for large banded problems and the shipped drivers
   are parallel), so the module is validated directly — an implicit
   1D heat equation solved with ARKStep + SPGMR + ARKBandPrecInit is
   checked against the exact solution of the semi-discrete system,
   and the internally generated difference-quotient Jacobian is
   compared against the analytic tridiagonal Jacobian.
   -----------------------------------------------------------------*/
#[cfg(test)]
mod tests {
    use super::*;
    use crate::arkode::{ARKodeEvolve, ARKodeFree, ARKodeSStolerances};
    use crate::arkode_arkstep::ARKStepCreate;
    use crate::arkode_impl::ARK_NORMAL;
    use crate::arkode_ls::ARKodeSetLinearSolver;
    use crate::arkode_ls_impl::ARKLS_LMEM_NULL;
    use crate::nvector_serial::N_VNew_Serial;
    use crate::sundials_linearsolver::SUN_PREC_LEFT;
    use crate::sundials_types::UserData;
    use crate::sunlinsol_dense::SUNLinSol_Dense;
    use crate::sunlinsol_spgmr::SUNLinSol_SPGMR;
    use crate::sunmatrix_dense::SUNDenseMatrix;
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

    #[test]
    fn bandpre_heat1d_spgmr() {
        let ctx = SUNContext_Create();
        let dx = 1.0 / (NEQ as f64 + 1.0);

        /* y0 = sin(pi x): eigenvector of the discrete Laplacian, so the
           exact solution of the ODE system is exp(lambda t) sin(pi x)
           with lambda = (2 cos(pi dx) - 2)/dx^2 */
        let mut y = N_VNew_Serial(NEQ, &ctx);
        for j in 0..NEQ as usize {
            y.data[j] = (PI * dx * (j as f64 + 1.0)).sin();
        }

        let mut arkode_mem =
            ARKStepCreate(None, Some(heat_rhs), 0.0, &y, &ctx).expect("ARKStepCreate");
        assert_eq!(ARKodeSStolerances(&mut arkode_mem, 1.0e-6, 1.0e-10), 0);

        /* attach SPGMR + banded preconditioner (mu = ml = 1) */
        let ls = SUNLinSol_SPGMR(&y, SUN_PREC_LEFT, 0, &ctx);
        assert_eq!(ARKodeSetLinearSolver(&mut arkode_mem, ls, None), ARKLS_SUCCESS);
        assert_eq!(ARKBandPrecInit(&mut arkode_mem, NEQ, 1, 1), ARKLS_SUCCESS);

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

        /* preconditioner RHS evaluations: each fresh DQ Jacobian costs
           ngroups = ml + mu + 1 = 3 evaluations */
        let mut nfeBP: i64 = 0;
        assert_eq!(ARKBandPrecGetNumRhsEvals(&mut arkode_mem, &mut nfeBP), ARKLS_SUCCESS);
        assert!(nfeBP > 0 && nfeBP % 3 == 0, "nfeBP = {}", nfeBP);

        /* workspace query */
        let mut lenrw: i64 = 0;
        let mut leniw: i64 = 0;
        assert_eq!(
            ARKBandPrecGetWorkSpace(&mut arkode_mem, &mut lenrw, &mut leniw),
            ARKLS_SUCCESS
        );
        assert!(lenrw > 0 && leniw >= 4, "lenrw = {}, leniw = {}", lenrw, leniw);

        /* the saved DQ Jacobian of the linear RHS must match the
           analytic tridiagonal Jacobian to difference-quotient
           rounding accuracy */
        let coef = 1.0 / (dx * dx);
        let arkls_mem = arkode_mem.lmem.as_ref().unwrap();
        let bp = match &arkls_mem.prec_module {
            PrecModule::BandPre(bp) => bp,
            _ => panic!("prec_module is not BandPre"),
        };
        let jm = match &bp.savedJ {
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

        /* free path exercises ARKBandPrecFree via arkLsFree */
        let mut slot = Some(arkode_mem);
        ARKodeFree(&mut slot);
    }

    #[test]
    fn bandpre_error_returns() {
        let ctx = SUNContext_Create();
        let y = N_VNew_Serial(NEQ, &ctx);
        let mut arkode_mem =
            ARKStepCreate(None, Some(heat_rhs), 0.0, &y, &ctx).expect("ARKStepCreate");

        /* no linear solver attached: ARKLS_LMEM_NULL */
        assert_eq!(ARKBandPrecInit(&mut arkode_mem, NEQ, 1, 1), ARKLS_LMEM_NULL);

        /* iterative LS attached but no preconditioner module:
           ARKLS_PMEM_NULL from the get routines */
        let ls = SUNLinSol_SPGMR(&y, SUN_PREC_LEFT, 0, &ctx);
        assert_eq!(ARKodeSetLinearSolver(&mut arkode_mem, ls, None), ARKLS_SUCCESS);
        let mut nfeBP: i64 = 0;
        assert_eq!(
            ARKBandPrecGetNumRhsEvals(&mut arkode_mem, &mut nfeBP),
            ARKLS_PMEM_NULL
        );

        /* direct (non-iterative) LS: ARKLS_ILL_INPUT, matching the C
           ARKodeSetPreconditioner error inside ARKBandPrecInit */
        let a = SUNDenseMatrix(NEQ, NEQ, &ctx);
        let dls = SUNLinSol_Dense(&y, &a, &ctx);
        assert_eq!(ARKodeSetLinearSolver(&mut arkode_mem, dls, Some(a)), ARKLS_SUCCESS);
        assert_eq!(ARKBandPrecInit(&mut arkode_mem, NEQ, 1, 1), ARKLS_ILL_INPUT);
    }
}
