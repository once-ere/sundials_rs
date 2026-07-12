/* -----------------------------------------------------------------
 * Translation of
 * sundials-7.7.0/src/sundomeigest/power/sundomeigest_power.c
 * (+ include/sundomeigest/sundomeigest_power.h).
 *
 * Power Iteration (PI) implementation of the SUNDomEigEst package.
 * The C content struct is the Power variant payload of the
 * SUNDomEigEstimator enum; implementation ops called on an estimator
 * of another variant return SUN_ERR_ARG_INCOMPATIBLE (in C this would
 * be undefined behavior through a mismatched content cast).
 *
 * Adaptations (see also sundials_domeigestimator.rs):
 *  - ATimes is supplied at Estimate time as a closure argument
 *    (the C SetATimes-stored callback cannot live inside the
 *    integrator that owns the estimator; same pinned adaptation as
 *    LinearSolver.solve); the C
 *    void* ATdata is captured by the closure and drops out.
 *  - The constructor returns the enum directly; the C SUNAssertNull
 *    argument checks (q/ops non-NULL, required ops present, q != 0)
 *    are debug-only in C (SUNDIALS_ENABLE_ERROR_CHECKS) and their
 *    NULL-return paths have no Rust counterpart. The SUNAssert
 *    guards on a missing ATimes in Initialize/Estimate are kept as
 *    early returns (the un-set state is representable here).
 *  - Destroy_Power is ownership drop (see the base module).
 * -----------------------------------------------------------------*/

use crate::nvector_serial::{NVector, N_VClone, N_VDotProd, N_VScale};
use crate::sundials_context::SUNContext;
use crate::sundials_domeigestimator::SUNDomEigEstimator;
use crate::sundials_errors::{
    SUNErrCode, SUN_ERR_ARG_CORRUPT, SUN_ERR_ARG_INCOMPATIBLE, SUN_ERR_USER_FCN_FAIL, SUN_SUCCESS,
};
use crate::sundials_linearsolver::ATimesFn;
use crate::sundials_math::{SUNRabs, SUNRsqrt};
use crate::sundials_types::SUN_SMALL_REAL;
use crate::sundials_utils::fmt_g;

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/* Default estimator parameters */
const DEE_NUM_OF_WARMUPS_PI_DEFAULT: i32 = 100;

/* Default Power Iteration parameters */
const DEE_TOL_DEFAULT: f64 = 0.005;
const DEE_MAX_ITER_DEFAULT: i64 = 100;

/// C struct SUNDomEigEstimatorContent_Power_
pub struct SUNDomEigEstimatorContent_Power {
    /// workspace vectors
    pub V: NVector,
    pub q: NVector,
    /// Number of preprocessing iterations
    pub num_warmups: i32,
    /// Maximum number of power iterations
    pub max_iters: i64,
    /// Number of iterations in last Estimate call
    pub num_iters: i64,
    /// Number of ATimes calls
    pub num_ATimes: i64,
    /// Convergence criteria for the power iteration
    pub rel_tol: f64,
    /// Residual from the last Estimate call
    pub res: f64,
}

fn content(
    dee: &SUNDomEigEstimator,
) -> Result<&SUNDomEigEstimatorContent_Power, SUNErrCode> {
    match dee {
        SUNDomEigEstimator::Power(c) => Ok(c),
        _ => Err(SUN_ERR_ARG_INCOMPATIBLE),
    }
}

fn content_mut(
    dee: &mut SUNDomEigEstimator,
) -> Result<&mut SUNDomEigEstimatorContent_Power, SUNErrCode> {
    match dee {
        SUNDomEigEstimator::Power(c) => Ok(c),
        _ => Err(SUN_ERR_ARG_INCOMPATIBLE),
    }
}

/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/// Function to create a new PI estimator
pub fn SUNDomEigEstimator_Power(
    q: &NVector,
    max_iters: i64,
    rel_tol: f64,
    _sunctx: &SUNContext,
) -> SUNDomEigEstimator {
    /* check for max_iters values; if illegal use defaults */
    let max_iters = if max_iters <= 0 {
        DEE_MAX_ITER_DEFAULT
    } else {
        max_iters
    };

    /* Check if rel_tol > 0 */
    let rel_tol = if rel_tol < SUN_SMALL_REAL {
        DEE_TOL_DEFAULT
    } else {
        rel_tol
    };

    /* Allocate content: content->q = N_VClone(q); N_VScale(ONE, q, content->q) */
    let mut cq = N_VClone(q);
    N_VScale(ONE, q, &mut cq);
    let cv = N_VClone(q);

    SUNDomEigEstimator::Power(SUNDomEigEstimatorContent_Power {
        V: cv,
        q: cq,
        max_iters,
        num_warmups: DEE_NUM_OF_WARMUPS_PI_DEFAULT,
        rel_tol,
        res: ZERO,
        num_iters: 0,
        num_ATimes: 0,
    })
}

/* -----------------------------------------------------------------
 * implementation of dominant eigenvalue estimator operations
 * ----------------------------------------------------------------- */

pub fn SUNDomEigEstimator_Initialize_Power(dee: &mut SUNDomEigEstimator) -> SUNErrCode {
    let c = match content_mut(dee) {
        Ok(c) => c,
        Err(e) => return e,
    };

    if c.rel_tol < SUN_SMALL_REAL {
        c.rel_tol = DEE_TOL_DEFAULT;
    }
    if c.num_warmups < 0 {
        c.num_warmups = DEE_NUM_OF_WARMUPS_PI_DEFAULT;
    }
    if c.max_iters <= 0 {
        c.max_iters = DEE_MAX_ITER_DEFAULT;
    }

    /* (ATimes presence check dropped: it arrives at Estimate time) */

    /* Initialize the vector V */
    let mut normq = N_VDotProd(&c.q, &c.q);
    normq = SUNRsqrt(normq);

    N_VScale(ONE / normq, &c.q, &mut c.V);

    SUN_SUCCESS
}

pub fn SUNDomEigEstimator_SetNumPreprocessIters_Power(
    dee: &mut SUNDomEigEstimator,
    num_iters: i32,
) -> SUNErrCode {
    let c = match content_mut(dee) {
        Ok(c) => c,
        Err(e) => return e,
    };
    /* Check if num_iters >= 0 */
    let num_iters = if num_iters < 0 {
        DEE_NUM_OF_WARMUPS_PI_DEFAULT
    } else {
        num_iters
    };
    /* set the number of warmups */
    c.num_warmups = num_iters;
    SUN_SUCCESS
}

pub fn SUNDomEigEstimator_SetRelTol_Power(
    dee: &mut SUNDomEigEstimator,
    rel_tol: f64,
) -> SUNErrCode {
    let c = match content_mut(dee) {
        Ok(c) => c,
        Err(e) => return e,
    };
    /* Check if rel_tol > 0 */
    let rel_tol = if rel_tol < SUN_SMALL_REAL {
        DEE_TOL_DEFAULT
    } else {
        rel_tol
    };
    /* set the tolerance */
    c.rel_tol = rel_tol;
    SUN_SUCCESS
}

pub fn SUNDomEigEstimator_SetMaxIters_Power(
    dee: &mut SUNDomEigEstimator,
    max_iters: i64,
) -> SUNErrCode {
    let c = match content_mut(dee) {
        Ok(c) => c,
        Err(e) => return e,
    };
    /* Check for legal number of iters */
    let max_iters = if max_iters <= 0 {
        DEE_MAX_ITER_DEFAULT
    } else {
        max_iters
    };
    /* Set max iters */
    c.max_iters = max_iters;
    SUN_SUCCESS
}

pub fn SUNDomEigEstimator_SetInitialGuess_Power(
    dee: &mut SUNDomEigEstimator,
    q: &NVector,
) -> SUNErrCode {
    let c = match content_mut(dee) {
        Ok(c) => c,
        Err(e) => return e,
    };

    let mut normq = N_VDotProd(q, q);
    normq = SUNRsqrt(normq);

    /* set the initial guess */
    N_VScale(ONE / normq, q, &mut c.V);

    SUN_SUCCESS
}

pub fn SUNDomEigEstimator_Estimate_Power(
    dee: &mut SUNDomEigEstimator,
    atimes: &mut ATimesFn,
    lambdaR: &mut f64,
    lambdaI: &mut f64,
) -> SUNErrCode {
    let c = match content_mut(dee) {
        Ok(c) => c,
        Err(e) => return e,
    };

    let mut newlambdaR = ZERO;
    let mut oldlambdaR = ZERO;

    let mut normq: f64;
    c.num_ATimes = 0;
    c.num_iters = 0;

    for _i in 0..c.num_warmups {
        let retval = atimes(&c.V, &mut c.q);
        c.num_ATimes += 1;
        c.num_iters += 1;
        if retval != 0 {
            return SUN_ERR_USER_FCN_FAIL;
        }

        normq = N_VDotProd(&c.q, &c.q);
        normq = SUNRsqrt(normq);
        N_VScale(ONE / normq, &c.q, &mut c.V);
    }

    for _k in 0..c.max_iters {
        let retval = atimes(&c.V, &mut c.q);
        c.num_ATimes += 1;
        c.num_iters += 1;
        if retval != 0 {
            return SUN_ERR_USER_FCN_FAIL;
        }

        newlambdaR = N_VDotProd(&c.V, &c.q); /* Rayleigh quotient */

        c.res = SUNRabs(newlambdaR - oldlambdaR) / SUNRabs(newlambdaR);

        if c.res < c.rel_tol {
            break;
        }

        normq = N_VDotProd(&c.q, &c.q);
        normq = SUNRsqrt(normq);
        N_VScale(ONE / normq, &c.q, &mut c.V);

        oldlambdaR = newlambdaR;
    }

    *lambdaI = ZERO;
    *lambdaR = newlambdaR;

    SUN_SUCCESS
}

pub fn SUNDomEigEstimator_GetRes_Power(
    dee: &SUNDomEigEstimator,
    res: &mut f64,
) -> SUNErrCode {
    let c = match content(dee) {
        Ok(c) => c,
        Err(e) => return e,
    };
    *res = c.res;
    SUN_SUCCESS
}

pub fn SUNDomEigEstimator_GetNumIters_Power(
    dee: &SUNDomEigEstimator,
    num_iters: &mut i64,
) -> SUNErrCode {
    let c = match content(dee) {
        Ok(c) => c,
        Err(e) => return e,
    };
    *num_iters = c.num_iters;
    SUN_SUCCESS
}

pub fn SUNDomEigEstimator_GetNumATimesCalls_Power(
    dee: &SUNDomEigEstimator,
    num_ATimes: &mut i64,
) -> SUNErrCode {
    let c = match content(dee) {
        Ok(c) => c,
        Err(e) => return e,
    };
    *num_ATimes = c.num_ATimes;
    SUN_SUCCESS
}

pub fn SUNDomEigEstimator_Write_Power(
    dee: &SUNDomEigEstimator,
    outfile: &mut dyn std::io::Write,
) -> SUNErrCode {
    let c = match content(dee) {
        Ok(c) => c,
        Err(e) => return e,
    };
    let r = (|| -> std::io::Result<()> {
        write!(outfile, "\nPower Iteration SUNDomEigEstimator:\n")?;
        write!(outfile, "Max. iters               = {}\n", c.max_iters)?;
        write!(outfile, "Num. preprocessing iters = {}\n", c.num_warmups)?;
        write!(
            outfile,
            "Relative tolerance       = {}\n",
            fmt_g(c.rel_tol, 0, 15)
        )?;
        write!(
            outfile,
            "Residual                 = {}\n",
            fmt_g(c.res, 0, 15)
        )?;
        write!(outfile, "Num. iters               = {}\n", c.num_iters)?;
        write!(outfile, "Num. ATimes calls        = {}\n\n", c.num_ATimes)?;
        Ok(())
    })();
    if r.is_err() {
        return SUN_ERR_ARG_CORRUPT;
    }
    SUN_SUCCESS
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sundials_context::SUNContext;
    use crate::sundials_domeigestimator::*;

    /// Dense matrix-vector product ATimes closure (owns its matrix).
    fn dense_atimes(a: Vec<Vec<f64>>) -> Box<ATimesFn<'static>> {
        Box::new(move |v: &NVector, av: &mut NVector| -> i32 {
            let n = a.len();
            for i in 0..n {
                let mut s = 0.0;
                for j in 0..n {
                    s += a[i][j] * v.data[j];
                }
                av.data[i] = s;
            }
            0
        })
    }

    #[test]
    fn constructor_defaults() {
        let sunctx = SUNContext::default();
        let q = NVector::from_slice(&[1.0, 1.0, 1.0]);
        /* illegal max_iters/rel_tol fall back onto the defaults */
        let dee = SUNDomEigEstimator_Power(&q, 0, -1.0, &sunctx);
        match &dee {
            SUNDomEigEstimator::Power(c) => {
                assert_eq!(c.max_iters, 100);
                assert_eq!(c.rel_tol, 0.005);
                assert_eq!(c.num_warmups, 100);
                assert_eq!(c.num_iters, 0);
                assert_eq!(c.num_ATimes, 0);
                assert_eq!(c.res, 0.0);
                assert_eq!(c.q.data, vec![1.0, 1.0, 1.0]);
                assert_eq!(c.V.len(), 3);
                    }
            _ => panic!(),
        }
        assert_eq!(SUNDomEigEstimator_Destroy(dee), SUN_SUCCESS);
    }

    #[test]
    fn estimate_diag_1_2_10_via_dispatch() {
        let sunctx = SUNContext::default();
        let q = NVector::from_slice(&[1.0, 1.0, 1.0]);
        let mut dee = SUNDomEigEstimator_Power(&q, 1000, 1e-10, &sunctx);

        let a = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 2.0, 0.0],
            vec![0.0, 0.0, 10.0],
        ];
        let mut atimes = dense_atimes(a);
        assert_eq!(SUNDomEigEstimator_SetNumPreprocessIters(&mut dee, 5), SUN_SUCCESS);
        assert_eq!(SUNDomEigEstimator_Initialize(&mut dee), SUN_SUCCESS);

        let (mut lr, mut li) = (0.0, -1.0);
        assert_eq!(SUNDomEigEstimator_Estimate(&mut dee, &mut atimes, &mut lr, &mut li), SUN_SUCCESS);
        assert!((lr - 10.0).abs() < 1e-6, "lambdaR = {}", lr);
        assert_eq!(li, 0.0);

        /* residual satisfies the convergence criterion */
        let mut res = -1.0;
        assert_eq!(SUNDomEigEstimator_GetRes(&dee, &mut res), SUN_SUCCESS);
        assert!(res < 1e-10);

        /* every iteration (warmup + power) performs exactly one ATimes call */
        let (mut ni, mut na) = (0i64, 0i64);
        assert_eq!(SUNDomEigEstimator_GetNumIters(&dee, &mut ni), SUN_SUCCESS);
        assert_eq!(SUNDomEigEstimator_GetNumATimesCalls(&dee, &mut na), SUN_SUCCESS);
        assert_eq!(ni, na);
        assert!(ni > 5); /* 5 warmups plus at least one power iteration */
        assert!(ni <= 5 + 1000);
    }

    #[test]
    fn set_initial_guess_and_reestimate() {
        let sunctx = SUNContext::default();
        let q = NVector::from_slice(&[1.0, 1.0]);
        let mut dee = SUNDomEigEstimator_Power(&q, 500, 1e-12, &sunctx);
        let a = vec![vec![3.0, 1.0], vec![1.0, 3.0]]; /* eigenvalues 4, 2 */
        let mut atimes = dense_atimes(a);
        assert_eq!(SUNDomEigEstimator_SetNumPreprocessIters(&mut dee, 0), SUN_SUCCESS);
        assert_eq!(SUNDomEigEstimator_Initialize(&mut dee), SUN_SUCCESS);
        /* a non-normalized guess must be normalized into V */
        let g = NVector::from_slice(&[3.0, -1.0]);
        assert_eq!(SUNDomEigEstimator_SetInitialGuess(&mut dee, &g), SUN_SUCCESS);
        match &dee {
            SUNDomEigEstimator::Power(c) => {
                let norm = (c.V.data[0] * c.V.data[0] + c.V.data[1] * c.V.data[1]).sqrt();
                assert!((norm - 1.0).abs() < 1e-15);
            }
            _ => panic!(),
        }
        let (mut lr, mut li) = (0.0, 0.0);
        assert_eq!(SUNDomEigEstimator_Estimate(&mut dee, &mut atimes, &mut lr, &mut li), SUN_SUCCESS);
        assert!((lr - 4.0).abs() < 1e-8, "lambdaR = {}", lr);
        assert_eq!(li, 0.0);
    }

    #[test]
    fn atimes_failure_is_reported() {
        let sunctx = SUNContext::default();
        let q = NVector::from_slice(&[1.0, 1.0]);
        let mut dee = SUNDomEigEstimator_Power(&q, 10, 0.1, &sunctx);
        let mut atimes: Box<ATimesFn<'static>> =
            Box::new(|_v: &NVector, _av: &mut NVector| -> i32 { 1 });
        assert_eq!(SUNDomEigEstimator_Initialize(&mut dee), SUN_SUCCESS);
        let (mut lr, mut li) = (0.0, 0.0);
        assert_eq!(
            SUNDomEigEstimator_Estimate(&mut dee, &mut atimes, &mut lr, &mut li),
            SUN_ERR_USER_FCN_FAIL
        );
    }

    #[test]
    fn setters_enforce_defaults_and_options_parse() {
        let sunctx = SUNContext::default();
        let q = NVector::from_slice(&[1.0, 1.0]);
        let mut dee = SUNDomEigEstimator_Power(&q, 10, 0.1, &sunctx);

        /* illegal values fall back onto defaults */
        assert_eq!(SUNDomEigEstimator_SetMaxIters(&mut dee, -3), SUN_SUCCESS);
        assert_eq!(SUNDomEigEstimator_SetRelTol(&mut dee, 0.0), SUN_SUCCESS);
        assert_eq!(SUNDomEigEstimator_SetNumPreprocessIters(&mut dee, -1), SUN_SUCCESS);
        match &dee {
            SUNDomEigEstimator::Power(c) => {
                assert_eq!(c.max_iters, 100);
                assert_eq!(c.rel_tol, 0.005);
                assert_eq!(c.num_warmups, 100);
            }
            _ => panic!(),
        }

        /* base-class command-line option processing */
        let args: Vec<String> = [
            "prog",
            "sundomeigestimator.max_iters",
            "250",
            "other.max_iters", /* wrong prefix: skipped */
            "999",
            "sundomeigestimator.num_preprocess_iters",
            "7",
            "sundomeigestimator.rel_tol",
            "0.001",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
        assert_eq!(
            SUNDomEigEstimator_SetOptions(&mut dee, None, None, &args),
            SUN_SUCCESS
        );
        match &dee {
            SUNDomEigEstimator::Power(c) => {
                assert_eq!(c.max_iters, 250);
                assert_eq!(c.num_warmups, 7);
                assert_eq!(c.rel_tol, 0.001);
            }
            _ => panic!(),
        }

        /* custom identifier prefix */
        let args: Vec<String> = ["prog", "myid.max_iters", "42"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert_eq!(
            SUNDomEigEstimator_SetOptions(&mut dee, Some("myid"), None, &args),
            SUN_SUCCESS
        );
        match &dee {
            SUNDomEigEstimator::Power(c) => assert_eq!(c.max_iters, 42),
            _ => panic!(),
        }

        /* file-based option control is unimplemented */
        assert_eq!(
            SUNDomEigEstimator_SetOptions(&mut dee, None, Some("opts.txt"), &args),
            SUN_ERR_ARG_INCOMPATIBLE
        );
    }

    #[test]
    fn write_output() {
        let sunctx = SUNContext::default();
        let q = NVector::from_slice(&[1.0, 1.0, 1.0]);
        let mut dee = SUNDomEigEstimator_Power(&q, 1000, 1e-10, &sunctx);
        let a = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 2.0, 0.0],
            vec![0.0, 0.0, 10.0],
        ];
        let mut atimes = dense_atimes(a);
        assert_eq!(SUNDomEigEstimator_SetNumPreprocessIters(&mut dee, 5), SUN_SUCCESS);
        assert_eq!(SUNDomEigEstimator_Initialize(&mut dee), SUN_SUCCESS);
        let (mut lr, mut li) = (0.0, 0.0);
        assert_eq!(SUNDomEigEstimator_Estimate(&mut dee, &mut atimes, &mut lr, &mut li), SUN_SUCCESS);

        let (mut ni, mut na) = (0i64, 0i64);
        assert_eq!(SUNDomEigEstimator_GetNumIters(&dee, &mut ni), SUN_SUCCESS);
        assert_eq!(SUNDomEigEstimator_GetNumATimesCalls(&dee, &mut na), SUN_SUCCESS);

        let mut buf: Vec<u8> = Vec::new();
        assert_eq!(SUNDomEigEstimator_Write(&dee, &mut buf), SUN_SUCCESS);
        let s = String::from_utf8(buf).unwrap();
        assert!(s.starts_with("\nPower Iteration SUNDomEigEstimator:\n"));
        assert!(s.contains("Max. iters               = 1000\n"));
        assert!(s.contains("Num. preprocessing iters = 5\n"));
        assert!(s.contains("Relative tolerance       = 1e-10\n"));
        assert!(s.contains("Residual                 = "));
        assert!(s.contains(&format!("Num. iters               = {}\n", ni)));
        assert!(s.ends_with(&format!("Num. ATimes calls        = {}\n\n", na)));
    }

    #[test]
    fn wrong_variant_returns_incompatible() {
        let sunctx = SUNContext::default();
        let q = NVector::from_slice(&[1.0, 1.0]);
        let mut arnoldi = crate::sundomeigest_arnoldi::SUNDomEigEstimator_Arnoldi(&q, 3, &sunctx);
        let mut res = 0.0;
        assert_eq!(
            SUNDomEigEstimator_GetRes_Power(&arnoldi, &mut res),
            SUN_ERR_ARG_INCOMPATIBLE
        );
        assert_eq!(
            SUNDomEigEstimator_SetMaxIters_Power(&mut arnoldi, 10),
            SUN_ERR_ARG_INCOMPATIBLE
        );
        /* base-class dispatch handles the missing Arnoldi ops instead */
        assert_eq!(SUNDomEigEstimator_SetMaxIters(&mut arnoldi, 10), SUN_SUCCESS);
        assert_eq!(SUNDomEigEstimator_SetRelTol(&mut arnoldi, 0.5), SUN_SUCCESS);
        let mut r = -1.0;
        assert_eq!(SUNDomEigEstimator_GetRes(&arnoldi, &mut r), SUN_SUCCESS);
        assert_eq!(r, 0.0);
    }
}
