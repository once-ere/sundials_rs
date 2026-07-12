/* -----------------------------------------------------------------
 * Translated from src/arkode/arkode_relaxation.c (ARKODE 7.7.0).
 * ARKODE's relaxation (in time) functionality.
 *
 * Temporary vectors utilized in the functions below:
 *   tempv2 - holds delta_y, the update direction vector
 *   tempv3 - holds y_relax, the relaxed solution vector
 *   tempv4 - holds J_relax, the Jacobian of the relaxation function
 *
 * Storage adaptation: C keeps ARKodeRelaxMem behind a pointer on
 * ARKodeMem; here it is Option<ARKodeRelaxMem> in place.  The
 * option/stat routines borrow it in place; arkRelax detaches it for
 * the duration of the solve (the residual/Jacobian/solver internals
 * need relax_mem and the ark_mem temp vectors simultaneously) and
 * puts it back before returning.
 * -----------------------------------------------------------------*/
use crate::arkode_impl::*;
use crate::arkode_relaxation_impl::*;
use crate::nvector_serial::{N_VDotProd, N_VLinearSum};
use crate::sundials_math::{SUNRabs, SUNRcopysign, SUNRpowerI, SUNRsamesign, SUNMIN};
use crate::sundials_types::SUNOutputFormat;

const ZERO: f64 = 0.0;
const HALF: f64 = 0.5;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;
const THREE: f64 = 3.0;

/* =============================================================================
 * Private Functions
 * ===========================================================================*/

/* C arkRelaxAccessMem: the relax_mem NULL check for the option/stat
   routines (the ark_mem NULL half is subsumed by &mut ARKodeMem). */
macro_rules! relax_access {
    ($ark_mem:expr, $fname:expr) => {
        if $ark_mem.relax_mem.is_none() {
            arkProcessError(
                Some($ark_mem),
                ARK_RELAX_MEM_NULL,
                line!(),
                $fname,
                file!(),
                MSG_RELAX_MEM_NULL,
            );
            return ARK_RELAX_MEM_NULL;
        }
    };
}

/* C's step_supports_relaxation guard shared by every option/stat
   routine below. */
macro_rules! relax_supported {
    ($ark_mem:expr, $fname:expr) => {
        if !$ark_mem.step_supports_relaxation {
            arkProcessError(
                Some($ark_mem),
                ARK_STEPPER_UNSUPPORTED,
                line!(),
                $fname,
                file!(),
                "time-stepping module does not support relaxation",
            );
            return ARK_STEPPER_UNSUPPORTED;
        }
    };
}

/* Evaluates the relaxation residual function */
fn arkRelaxResidual(
    relax_param: f64,
    relax_res: &mut f64,
    ark_mem: &mut ARKodeMem,
    relax_mem: &mut ARKodeRelaxMem,
) -> i32 {
    let e_old = relax_mem.e_old;
    let delta_e = relax_mem.delta_e;

    /* y_relax = y_n + r * delta_y  (delta_y = tempv2, y_relax = tempv3) */
    {
        let ARKodeMem { yn, tempv2, tempv3, .. } = ark_mem;
        N_VLinearSum(ONE, yn, relax_param, tempv2, tempv3);
    }

    /* Evaluate entropy function */
    let relax_fn = relax_mem.relax_fn.unwrap();
    let retval = relax_fn(&ark_mem.tempv3, relax_res, &mut ark_mem.user_data);
    relax_mem.num_relax_fn_evals += 1;
    if retval < 0 {
        return ARK_RELAX_FUNC_FAIL;
    }
    if retval > 0 {
        return ARK_RELAX_FUNC_RECV;
    }

    /* Compute relaxation residual */
    *relax_res = *relax_res - e_old - relax_param * delta_e;

    ARK_SUCCESS
}

/* Evaluates the Jacobian of the relaxation residual function */
fn arkRelaxResidualJacobian(
    relax_param: f64,
    relax_jac: &mut f64,
    ark_mem: &mut ARKodeMem,
    relax_mem: &mut ARKodeRelaxMem,
) -> i32 {
    let delta_e = relax_mem.delta_e;

    /* y_relax = y_n + r * delta_y  (delta_y = tempv2, y_relax = tempv3) */
    {
        let ARKodeMem { yn, tempv2, tempv3, .. } = ark_mem;
        N_VLinearSum(ONE, yn, relax_param, tempv2, tempv3);
    }

    /* Evaluate Jacobian of entropy functions (J_relax = tempv4) */
    let relax_jac_fn = relax_mem.relax_jac_fn.unwrap();
    let retval = {
        let ARKodeMem { tempv3, tempv4, user_data, .. } = ark_mem;
        relax_jac_fn(tempv3, tempv4, user_data)
    };
    relax_mem.num_relax_jac_evals += 1;
    if retval < 0 {
        return ARK_RELAX_JAC_FAIL;
    }
    if retval > 0 {
        return ARK_RELAX_JAC_RECV;
    }

    /* Compute relaxation residual Jacobian */
    *relax_jac = N_VDotProd(&ark_mem.tempv2, &ark_mem.tempv4);
    *relax_jac -= delta_e;

    ARK_SUCCESS
}

/* Solve the relaxation residual equation using Newton's method */
fn arkRelaxNewtonSolve(ark_mem: &mut ARKodeMem, relax_mem: &mut ARKodeRelaxMem) -> i32 {
    for _i in 0..relax_mem.max_iters {
        /* Compute the current residual */
        let mut res = relax_mem.res;
        let retval = arkRelaxResidual(relax_mem.relax_param, &mut res, ark_mem, relax_mem);
        relax_mem.res = res;
        if retval != 0 {
            return retval;
        }

        /* Check for convergence */
        if SUNRabs(relax_mem.res) < relax_mem.res_tol {
            return ARK_SUCCESS;
        }

        /* Compute Jacobian */
        let mut jac = relax_mem.jac;
        let retval = arkRelaxResidualJacobian(relax_mem.relax_param, &mut jac, ark_mem, relax_mem);
        relax_mem.jac = jac;
        if retval != 0 {
            return retval;
        }

        /* Update step length tolerance and solution */
        let tol = relax_mem.rel_tol * SUNRabs(relax_mem.relax_param) + relax_mem.abs_tol;

        let delta = relax_mem.res / relax_mem.jac;
        relax_mem.relax_param -= delta;

        /* Update cumulative iteration count */
        relax_mem.nls_iters += 1;

        /* Check for small update */
        if SUNRabs(delta) < tol {
            return ARK_SUCCESS;
        }
    }

    ARK_RELAX_SOLVE_RECV
}

/* Solve the relaxation residual equation using Brent's method */
fn arkRelaxBrentSolve(ark_mem: &mut ARKodeMem, relax_mem: &mut ARKodeRelaxMem) -> i32 {
    let mut xa: f64; /* previous solution and function value */
    let mut fa: f64 = ZERO;
    let mut xb: f64; /* current solution and function value  */
    let mut fb: f64 = ZERO;
    let mut xc: f64; /* together brac and curr bracket zero  */
    let mut fc: f64;
    let mut xm: f64; /* midpoint between brac and curr       */
    let mut old_update: f64; /* previous iteration update      */
    let mut new_update: f64; /* new iteration update           */
    let mut tol: f64; /* iteration tolerance                  */
    let mut pt: f64; /* temporary values                      */
    let mut qt: f64;
    let mut rt: f64;
    let mut st: f64;

    /* Compute interval that brackets the root */
    xa = 0.9 * relax_mem.relax_param;
    xb = 1.1 * relax_mem.relax_param;

    for _i in 0..10 {
        /* Compute relaxation residual */
        let retval = arkRelaxResidual(xa, &mut fa, ark_mem, relax_mem);
        relax_mem.num_relax_fn_evals += 1;
        if retval < 0 {
            return ARK_RELAX_FUNC_FAIL;
        }
        if retval > 0 {
            return ARK_RELAX_FUNC_RECV;
        }

        /* Check if we got lucky */
        if SUNRabs(fa) < relax_mem.res_tol {
            relax_mem.res = fa;
            relax_mem.relax_param = xa;
            return ARK_SUCCESS;
        }

        if fa < ZERO {
            break;
        }

        fb = fa;
        xb = xa;
        xa *= 0.9;
    }
    if fa > ZERO {
        return ARK_RELAX_SOLVE_RECV;
    }

    for _i in 0..10 {
        /* Compute relaxation residual */
        let retval = arkRelaxResidual(xb, &mut fb, ark_mem, relax_mem);
        relax_mem.num_relax_fn_evals += 1;
        if retval < 0 {
            return ARK_RELAX_FUNC_FAIL;
        }
        if retval > 0 {
            return ARK_RELAX_FUNC_RECV;
        }

        /* Check if we got lucky */
        if SUNRabs(fb) < relax_mem.res_tol {
            relax_mem.res = fb;
            relax_mem.relax_param = xb;
            return ARK_SUCCESS;
        }

        if fb > ZERO {
            break;
        }

        fa = fb;
        xa = xb;
        xb *= 1.1;
    }
    if fb < ZERO {
        return ARK_RELAX_SOLVE_RECV;
    }

    /* Initialize values bracketing values to lower bound and updates */
    xc = xa;
    fc = fa;

    old_update = ZERO;
    new_update = ZERO;

    /* Find root */
    for _i in 0..relax_mem.max_iters {
        /* Ensure xc and xb bracket zero */
        if SUNRsamesign(fc, fb) {
            xc = xa;
            fc = fa;
            new_update = xb - xa;
            old_update = new_update;
        }

        /* Ensure xb is closer to zero than xc */
        if SUNRabs(fb) > SUNRabs(fc) {
            xa = xb;
            xb = xc;
            xc = xa;

            fa = fb;
            fb = fc;
            fc = fa;
        }

        /* Update tolerance */
        tol = relax_mem.rel_tol * SUNRabs(xb) + HALF * relax_mem.abs_tol;

        /* Compute midpoint for bisection */
        xm = HALF * (xc - xb);

        /* Check for convergence */
        if SUNRabs(xm) < tol || SUNRabs(fb) < relax_mem.res_tol {
            relax_mem.res = fb;
            relax_mem.relax_param = xb;
            return ARK_SUCCESS;
        }

        /* Compute iteration update */
        if SUNRabs(old_update) >= tol && SUNRabs(fb) < SUNRabs(fa) {
            /* Converging sufficiently fast, interpolate solution */
            st = fb / fa;

            if xa == xc {
                /* Two unique values available, try linear interpolant (secant) */
                pt = TWO * xm * st;
                qt = ONE - st;
            } else {
                /* Three unique values available, try inverse quadratic interpolant */
                qt = fa / fc;
                rt = fb / fc;
                pt = st * (TWO * xm * qt * (qt - rt) - (xb - xa) * (rt - ONE));
                qt = (qt - ONE) * (rt - ONE) * (st - ONE);
            }

            /* Ensure updates produce values within [xc, xb] or [xb, xc] */
            if pt > ZERO {
                qt = -qt;
            } else {
                pt = -pt;
            }

            /* Check if interpolant is acceptable, otherwise use bisection */
            st = THREE * xm * qt - SUNRabs(tol * qt);
            rt = SUNRabs(old_update * qt);

            if TWO * pt < SUNMIN(st, rt) {
                old_update = new_update;
                new_update = pt / qt;
            } else {
                new_update = xm;
                old_update = xm;
            }
        } else {
            /* Converging too slowly, use bisection */
            new_update = xm;
            old_update = xm;
        }

        /* Update solution */
        xa = xb;
        fa = fb;

        /* If update is small, use tolerance in bisection direction */
        if SUNRabs(new_update) > tol {
            xb += new_update;
        } else {
            xb += SUNRcopysign(tol, xm);
        }

        /* Compute relaxation residual */
        let retval = arkRelaxResidual(xb, &mut fb, ark_mem, relax_mem);
        relax_mem.num_relax_fn_evals += 1;
        if retval < 0 {
            return ARK_RELAX_FUNC_FAIL;
        }
        if retval > 0 {
            return ARK_RELAX_FUNC_RECV;
        }
    }

    ARK_RELAX_SOLVE_RECV
}

/* Compute and apply relaxation parameter */
fn arkRelaxSolve(
    ark_mem: &mut ARKodeMem,
    relax_mem: &mut ARKodeRelaxMem,
    relax_val_out: &mut f64,
) -> i32 {
    /* Get the change in entropy (uses temp vectors 2 and 3) */
    let delta_e_fn = relax_mem.delta_e_fn.unwrap();
    let retval = delta_e_fn(
        ark_mem,
        relax_mem.relax_jac_fn.unwrap(),
        &mut relax_mem.num_relax_jac_evals,
        &mut relax_mem.delta_e,
    );
    if retval != 0 {
        return retval;
    }

    /* Get the change in state (delta_y = tempv2) */
    {
        let ARKodeMem { ycur, yn, tempv2, .. } = ark_mem;
        N_VLinearSum(ONE, ycur, -ONE, yn, tempv2);
    }

    /* Store the current relaxation function value */
    let relax_fn = relax_mem.relax_fn.unwrap();
    let retval = relax_fn(&ark_mem.yn, &mut relax_mem.e_old, &mut ark_mem.user_data);
    relax_mem.num_relax_fn_evals += 1;
    if retval < 0 {
        return ARK_RELAX_FUNC_FAIL;
    }
    if retval > 0 {
        return ARK_RELAX_FUNC_RECV;
    }

    /* Initial guess for relaxation parameter */
    relax_mem.relax_param = relax_mem.relax_param_prev;

    let retval = match relax_mem.solver {
        ARK_RELAX_BRENT => arkRelaxBrentSolve(ark_mem, relax_mem),
        ARK_RELAX_NEWTON => arkRelaxNewtonSolve(ark_mem, relax_mem),
    };

    /* Check for solver failure */
    if retval != 0 {
        relax_mem.nls_fails += 1;
        return retval;
    }

    /* Check for bad relaxation value */
    if relax_mem.relax_param < relax_mem.lower_bound
        || relax_mem.relax_param > relax_mem.upper_bound
    {
        relax_mem.bound_fails += 1;
        return ARK_RELAX_SOLVE_RECV;
    }

    /* Save parameter for next initial guess */
    relax_mem.relax_param_prev = relax_mem.relax_param;

    /* Return relaxation value */
    *relax_val_out = relax_mem.relax_param;

    ARK_SUCCESS
}

/* =============================================================================
 * User Functions
 * ===========================================================================*/

/* -----------------------------------------------------------------------------
 * Set functions
 * ---------------------------------------------------------------------------*/

pub fn ARKodeSetRelaxFn(
    ark_mem: &mut ARKodeMem,
    rfn: Option<ARKRelaxFn>,
    rjac: Option<ARKRelaxJacFn>,
) -> i32 {
    /* Ensure that the current N_Vector supports N_VDotProd: always
       true for the serial NVector. */

    /* Call stepper-specific routine (if it exists) */
    if let Some(setrelaxfn) = ark_mem.step_setrelaxfn {
        setrelaxfn(ark_mem, rfn, rjac)
    } else {
        arkProcessError(
            Some(ark_mem),
            ARK_STEPPER_UNSUPPORTED,
            line!(),
            "ARKodeSetRelaxFn",
            file!(),
            "time-stepping module does not support relaxation",
        );
        ARK_STEPPER_UNSUPPORTED
    }
}

pub fn ARKodeSetRelaxEtaFail(ark_mem: &mut ARKodeMem, eta_fail: f64) -> i32 {
    relax_access!(ark_mem, "ARKodeSetRelaxEtaFail");
    relax_supported!(ark_mem, "ARKodeSetRelaxEtaFail");
    let relax_mem = ark_mem.relax_mem.as_mut().unwrap();

    if eta_fail > ZERO && eta_fail < ONE {
        relax_mem.eta_fail = eta_fail;
    } else {
        relax_mem.eta_fail = ARK_RELAX_DEFAULT_ETA_FAIL;
    }

    ARK_SUCCESS
}

pub fn ARKodeSetRelaxLowerBound(ark_mem: &mut ARKodeMem, lower: f64) -> i32 {
    relax_access!(ark_mem, "ARKodeSetRelaxLowerBound");
    relax_supported!(ark_mem, "ARKodeSetRelaxLowerBound");
    let relax_mem = ark_mem.relax_mem.as_mut().unwrap();

    if lower > ZERO && lower < ONE {
        relax_mem.lower_bound = lower;
    } else {
        relax_mem.lower_bound = ARK_RELAX_DEFAULT_LOWER_BOUND;
    }

    ARK_SUCCESS
}

pub fn ARKodeSetRelaxMaxFails(ark_mem: &mut ARKodeMem, max_fails: i32) -> i32 {
    relax_access!(ark_mem, "ARKodeSetRelaxMaxFails");
    relax_supported!(ark_mem, "ARKodeSetRelaxMaxFails");
    let relax_mem = ark_mem.relax_mem.as_mut().unwrap();

    if max_fails > 0 {
        relax_mem.max_fails = max_fails;
    } else {
        relax_mem.max_fails = ARK_RELAX_DEFAULT_MAX_FAILS;
    }

    ARK_SUCCESS
}

pub fn ARKodeSetRelaxMaxIters(ark_mem: &mut ARKodeMem, max_iters: i32) -> i32 {
    relax_access!(ark_mem, "ARKodeSetRelaxMaxIters");
    relax_supported!(ark_mem, "ARKodeSetRelaxMaxIters");
    let relax_mem = ark_mem.relax_mem.as_mut().unwrap();

    if max_iters > 0 {
        relax_mem.max_iters = max_iters;
    } else {
        relax_mem.max_iters = ARK_RELAX_DEFAULT_MAX_ITERS;
    }

    ARK_SUCCESS
}

pub fn ARKodeSetRelaxSolver(ark_mem: &mut ARKodeMem, solver: ARKRelaxSolver) -> i32 {
    relax_access!(ark_mem, "ARKodeSetRelaxSolver");
    relax_supported!(ark_mem, "ARKodeSetRelaxSolver");

    /* (the C invalid-enum check collapses: ARKRelaxSolver only has the
       BRENT and NEWTON variants) */
    let relax_mem = ark_mem.relax_mem.as_mut().unwrap();
    relax_mem.solver = solver;

    ARK_SUCCESS
}

pub fn ARKodeSetRelaxResTol(ark_mem: &mut ARKodeMem, res_tol: f64) -> i32 {
    relax_access!(ark_mem, "ARKodeSetRelaxResTol");
    relax_supported!(ark_mem, "ARKodeSetRelaxResTol");
    let relax_mem = ark_mem.relax_mem.as_mut().unwrap();

    if res_tol > ZERO {
        relax_mem.res_tol = res_tol;
    } else {
        relax_mem.res_tol = ARK_RELAX_DEFAULT_RES_TOL;
    }

    ARK_SUCCESS
}

pub fn ARKodeSetRelaxTol(ark_mem: &mut ARKodeMem, rel_tol: f64, abs_tol: f64) -> i32 {
    relax_access!(ark_mem, "ARKodeSetRelaxTol");
    relax_supported!(ark_mem, "ARKodeSetRelaxTol");
    let relax_mem = ark_mem.relax_mem.as_mut().unwrap();

    if rel_tol > ZERO {
        relax_mem.rel_tol = rel_tol;
    } else {
        relax_mem.rel_tol = ARK_RELAX_DEFAULT_REL_TOL;
    }

    if abs_tol > ZERO {
        relax_mem.abs_tol = abs_tol;
    } else {
        relax_mem.abs_tol = ARK_RELAX_DEFAULT_ABS_TOL;
    }

    ARK_SUCCESS
}

pub fn ARKodeSetRelaxUpperBound(ark_mem: &mut ARKodeMem, upper: f64) -> i32 {
    relax_access!(ark_mem, "ARKodeSetRelaxUpperBound");
    relax_supported!(ark_mem, "ARKodeSetRelaxUpperBound");
    let relax_mem = ark_mem.relax_mem.as_mut().unwrap();

    if upper > ONE {
        relax_mem.upper_bound = upper;
    } else {
        relax_mem.upper_bound = ARK_RELAX_DEFAULT_UPPER_BOUND;
    }

    ARK_SUCCESS
}

/* -----------------------------------------------------------------------------
 * Get functions
 * ---------------------------------------------------------------------------*/

pub fn ARKodeGetNumRelaxFnEvals(ark_mem: &mut ARKodeMem, r_evals: &mut i64) -> i32 {
    relax_access!(ark_mem, "ARKodeGetNumRelaxFnEvals");
    relax_supported!(ark_mem, "ARKodeGetNumRelaxFnEvals");

    *r_evals = ark_mem.relax_mem.as_ref().unwrap().num_relax_fn_evals;

    ARK_SUCCESS
}

pub fn ARKodeGetNumRelaxJacEvals(ark_mem: &mut ARKodeMem, J_evals: &mut i64) -> i32 {
    relax_access!(ark_mem, "ARKodeGetNumRelaxJacEvals");
    relax_supported!(ark_mem, "ARKodeGetNumRelaxJacEvals");

    *J_evals = ark_mem.relax_mem.as_ref().unwrap().num_relax_jac_evals;

    ARK_SUCCESS
}

pub fn ARKodeGetNumRelaxFails(ark_mem: &mut ARKodeMem, relax_fails: &mut i64) -> i32 {
    relax_access!(ark_mem, "ARKodeGetNumRelaxFails");
    relax_supported!(ark_mem, "ARKodeGetNumRelaxFails");

    *relax_fails = ark_mem.relax_mem.as_ref().unwrap().num_fails;

    ARK_SUCCESS
}

pub fn ARKodeGetNumRelaxSolveFails(ark_mem: &mut ARKodeMem, fails: &mut i64) -> i32 {
    relax_access!(ark_mem, "ARKodeGetNumRelaxSolveFails");
    relax_supported!(ark_mem, "ARKodeGetNumRelaxSolveFails");

    *fails = ark_mem.relax_mem.as_ref().unwrap().nls_fails;

    ARK_SUCCESS
}

pub fn ARKodeGetNumRelaxBoundFails(ark_mem: &mut ARKodeMem, fails: &mut i64) -> i32 {
    relax_access!(ark_mem, "ARKodeGetNumRelaxBoundFails");
    relax_supported!(ark_mem, "ARKodeGetNumRelaxBoundFails");

    *fails = ark_mem.relax_mem.as_ref().unwrap().bound_fails;

    ARK_SUCCESS
}

pub fn ARKodeGetNumRelaxSolveIters(ark_mem: &mut ARKodeMem, iters: &mut i64) -> i32 {
    relax_access!(ark_mem, "ARKodeGetNumRelaxSolveIters");
    relax_supported!(ark_mem, "ARKodeGetNumRelaxSolveIters");

    *iters = ark_mem.relax_mem.as_ref().unwrap().nls_iters;

    ARK_SUCCESS
}

/* =============================================================================
 * Driver and Stepper Functions
 * ===========================================================================*/

/* Constructor called by stepper */
pub fn arkRelaxCreate(
    ark_mem: &mut ARKodeMem,
    relax_fn: Option<ARKRelaxFn>,
    relax_jac_fn: Option<ARKRelaxJacFn>,
    delta_e_fn: Option<ARKRelaxDeltaEFn>,
    get_order_fn: Option<ARKRelaxGetOrderFn>,
) -> i32 {
    /* Disable relaxation if both user inputs are NULL */
    if relax_fn.is_none() && relax_jac_fn.is_none() {
        ark_mem.relax_enabled = false;
        return ARK_SUCCESS;
    }

    /* Ensure both the relaxation function and Jacobian are provided */
    if relax_fn.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "arkRelaxCreate",
            file!(),
            "The relaxation function is NULL.",
        );
        return ARK_ILL_INPUT;
    }

    if relax_jac_fn.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "arkRelaxCreate",
            file!(),
            "The relaxation Jacobian function is NULL.",
        );
        return ARK_ILL_INPUT;
    }

    /* Ensure stepper supplied inputs are provided */
    if delta_e_fn.is_none() || get_order_fn.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "arkRelaxCreate",
            file!(),
            "The Delta y, Delta e, or get order function is NULL.",
        );
        return ARK_ILL_INPUT;
    }

    /* Allocate and initialize relaxation memory structure */
    if ark_mem.relax_mem.is_none() {
        ark_mem.relax_mem = Some(ARKodeRelaxMem {
            /* (memset-0 for all counters/values) */
            relax_fn: None,
            relax_jac_fn: None,
            delta_e_fn: None,
            get_order_fn: None,
            /* Set defaults */
            max_fails: ARK_RELAX_DEFAULT_MAX_FAILS,
            num_relax_fn_evals: 0,
            num_relax_jac_evals: 0,
            num_fails: 0,
            e_old: ZERO,
            delta_e: ZERO,
            res: ZERO,
            jac: ZERO,
            relax_param: ZERO,
            /* Initialize values */
            relax_param_prev: ONE,
            lower_bound: ARK_RELAX_DEFAULT_LOWER_BOUND,
            upper_bound: ARK_RELAX_DEFAULT_UPPER_BOUND,
            eta_fail: ARK_RELAX_DEFAULT_ETA_FAIL,
            solver: ARK_RELAX_NEWTON,
            res_tol: ARK_RELAX_DEFAULT_RES_TOL,
            rel_tol: ARK_RELAX_DEFAULT_REL_TOL,
            abs_tol: ARK_RELAX_DEFAULT_ABS_TOL,
            max_iters: ARK_RELAX_DEFAULT_MAX_ITERS,
            nls_iters: 0,
            nls_fails: 0,
            bound_fails: 0,
        });

        /* Update workspace sizes */
        ark_mem.lrw += 12;
        ark_mem.liw += 14;
    }

    /* Set function pointers */
    let relax_mem = ark_mem.relax_mem.as_mut().unwrap();
    relax_mem.relax_fn = relax_fn;
    relax_mem.relax_jac_fn = relax_jac_fn;
    relax_mem.delta_e_fn = delta_e_fn;
    relax_mem.get_order_fn = get_order_fn;

    /* Enable relaxation */
    ark_mem.relax_enabled = true;

    ARK_SUCCESS
}

/* arkRelaxDestroy: dropping ARKodeMem.relax_mem releases the C
   free()d structure. */

/* Compute and apply relaxation, called by driver */
pub fn arkRelax(ark_mem: &mut ARKodeMem, relax_fails: &mut i32, dsm_inout: &mut f64) -> i32 {
    /* Get the relaxation memory structure */
    let mut relax_mem = match ark_mem.relax_mem.take() {
        Some(rm) => rm,
        None => {
            arkProcessError(
                Some(ark_mem),
                ARK_RELAX_MEM_NULL,
                line!(),
                "arkRelax",
                file!(),
                MSG_RELAX_MEM_NULL,
            );
            return ARK_RELAX_MEM_NULL;
        }
    };
    let retval = arkRelax_inner(ark_mem, &mut relax_mem, relax_fails, dsm_inout);
    ark_mem.relax_mem = Some(relax_mem);
    retval
}

fn arkRelax_inner(
    ark_mem: &mut ARKodeMem,
    relax_mem: &mut ARKodeRelaxMem,
    relax_fails: &mut i32,
    dsm_inout: &mut f64,
) -> i32 {
    /* Compute the relaxation parameter */
    let mut relax_val = ZERO;
    let retval = arkRelaxSolve(ark_mem, relax_mem, &mut relax_val);
    if retval < 0 {
        return retval;
    }
    if retval > 0 {
        /* Update failure counts */
        relax_mem.num_fails += 1;
        *relax_fails += 1;

        /* Check for max fails in a step */
        if *relax_fails == relax_mem.max_fails {
            return ARK_RELAX_FAIL;
        }

        /* Return with an error if |h| == hmin */
        if SUNRabs(ark_mem.h) <= ark_mem.hmin * ONEPSM {
            return ARK_RELAX_FAIL;
        }

        /* Return with error if using fixed step sizes */
        if ark_mem.fixedstep {
            return ARK_RELAX_FAIL;
        }

        /* Cut step size and try again */
        ark_mem.eta = relax_mem.eta_fail;

        return TRY_AGAIN;
    }

    /* Update step size and error estimate */
    ark_mem.h *= relax_val;
    *dsm_inout *= SUNRpowerI(relax_val, relax_mem.get_order_fn.unwrap()(ark_mem));

    /* Relax solution: ycur = relax_val*ycur + (1-relax_val)*yn (the C
       N_VLinearSum aliases its output with ycur -> in-place form) */
    {
        let ARKodeMem { ycur, yn, .. } = ark_mem;
        ycur.linear_sum_with(relax_val, ONE - relax_val, yn);
    }

    ARK_SUCCESS
}

/* Print relaxation solver statistics, called by ARKODE */
pub fn arkRelaxPrintAllStats(
    ark_mem: &mut ARKodeMem,
    outfile: &mut dyn std::io::Write,
    fmt: SUNOutputFormat,
) -> i32 {
    relax_access!(ark_mem, "arkRelaxPrintAllStats");
    let relax_mem = ark_mem.relax_mem.as_ref().unwrap();

    crate::arkode_io::sunfprintf_long(
        outfile,
        fmt,
        false,
        "Relax fn evals",
        relax_mem.num_relax_fn_evals,
    );
    crate::arkode_io::sunfprintf_long(
        outfile,
        fmt,
        false,
        "Relax Jac evals",
        relax_mem.num_relax_jac_evals,
    );
    crate::arkode_io::sunfprintf_long(outfile, fmt, false, "Relax fails", relax_mem.num_fails);
    crate::arkode_io::sunfprintf_long(
        outfile,
        fmt,
        false,
        "Relax bound fails",
        relax_mem.bound_fails,
    );
    crate::arkode_io::sunfprintf_long(outfile, fmt, false, "Relax NLS iters", relax_mem.nls_iters);
    crate::arkode_io::sunfprintf_long(outfile, fmt, false, "Relax NLS fails", relax_mem.nls_fails);

    ARK_SUCCESS
}
