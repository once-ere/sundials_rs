/* -----------------------------------------------------------------
 * Translation of sundials-7.7.0/src/arkode/arkode_relaxation_impl.h.
 * Implementation header for ARKODE's relaxation (entropy) support;
 * the routines (arkRelaxCreate / arkRelax / ...) live in
 * arkode_relaxation.rs (from arkode_relaxation.c).
 * -----------------------------------------------------------------*/

use crate::arkode_impl::{ARKRelaxFn, ARKRelaxJacFn, ARKRelaxSolver, ARKodeMem};
use crate::sundials_types::SUN_UNIT_ROUNDOFF;

/* Relaxation Constants */
pub const ARK_RELAX_DEFAULT_MAX_FAILS: i32 = 10;
pub const ARK_RELAX_DEFAULT_RES_TOL: f64 = 10.0 * SUN_UNIT_ROUNDOFF;
pub const ARK_RELAX_DEFAULT_REL_TOL: f64 = 4.0 * SUN_UNIT_ROUNDOFF;
pub const ARK_RELAX_DEFAULT_ABS_TOL: f64 = 1.0e-14;
pub const ARK_RELAX_DEFAULT_MAX_ITERS: i32 = 10;
pub const ARK_RELAX_DEFAULT_LOWER_BOUND: f64 = 0.8;
pub const ARK_RELAX_DEFAULT_UPPER_BOUND: f64 = 1.2;
pub const ARK_RELAX_DEFAULT_ETA_FAIL: f64 = 0.25;

/* Relaxation Private Return Values (see arkode_impl.rs for public values) */
pub const ARK_RELAX_FUNC_RECV: i32 = 1;
pub const ARK_RELAX_JAC_RECV: i32 = 2;
pub const ARK_RELAX_SOLVE_RECV: i32 = 3;

/* Stepper Supplied Relaxation Functions */

/// Compute the estimated change in entropy for this step delta_e
pub type ARKRelaxDeltaEFn = fn(
    ark_mem: &mut ARKodeMem,
    relax_jac_fn: ARKRelaxJacFn,
    evals_out: &mut i64,
    delta_e_out: &mut f64,
) -> i32;

/// Get the method order
pub type ARKRelaxGetOrderFn = fn(ark_mem: &mut ARKodeMem) -> i32;

/// struct ARKodeRelaxMemRec (arkode_relaxation_impl.h)
pub struct ARKodeRelaxMem {
    /* user-supplied and stepper supplied functions */
    pub relax_fn: Option<ARKRelaxFn>, /* user relaxation function ("entropy") */
    pub relax_jac_fn: Option<ARKRelaxJacFn>, /* user relaxation Jacobian     */
    pub delta_e_fn: Option<ARKRelaxDeltaEFn>, /* get delta entropy from stepper */
    pub get_order_fn: Option<ARKRelaxGetOrderFn>, /* get the method order     */

    /* relaxation variables */
    pub max_fails: i32,           /* max allowed relax fails in a step   */
    pub num_relax_fn_evals: i64,  /* counter for total function evals    */
    pub num_relax_jac_evals: i64, /* counter for total jacobian evals    */
    pub num_fails: i64,           /* counter for total relaxation fails  */
    pub e_old: f64,               /* entropy at start of step y(t_{n-1}) */
    pub delta_e: f64,             /* change in entropy                   */
    pub res: f64,                 /* relaxation residual value           */
    pub jac: f64,                 /* relaxation Jacobian value           */
    pub relax_param: f64,         /* current relaxation parameter value  */
    pub relax_param_prev: f64,    /* previous relaxation parameter value */
    pub lower_bound: f64,         /* smallest allowed relaxation value   */
    pub upper_bound: f64,         /* largest allowed relaxation value    */
    pub eta_fail: f64,            /* failed relaxation step size factor  */

    /* nonlinear solver settings */
    pub solver: ARKRelaxSolver, /* choice of relaxation solver          */
    pub res_tol: f64,           /* nonlinear residual solve tolerance   */
    pub rel_tol: f64,           /* nonlinear iterate relative tolerance */
    pub abs_tol: f64,           /* nonlinear iterate absolute tolerance */
    pub max_iters: i32,         /* nonlinear solve max iterations       */
    pub nls_iters: i64,         /* total nonlinear iterations           */
    pub nls_fails: i64,         /* number of nonlinear solver fails     */
    pub bound_fails: i64,       /* number of relax param bound fails    */
}

/* Error Messages */
pub const MSG_RELAX_MEM_NULL: &str = "Relaxation memory is NULL.";
