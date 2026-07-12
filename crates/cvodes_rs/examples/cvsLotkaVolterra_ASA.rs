/* -----------------------------------------------------------------
 * Translation of examples/cvodes/serial/cvsLotkaVolterra_ASA.c
 * (SUNDIALS 7.7.0).
 *
 * This example solves the Lotka-Volterra ODE with four parameters,
 *
 *     u = [dx/dt] = [ p_0*x - p_1*x*y  ]
 *         [dy/dt]   [ -p_2*y + p_3*x*y ].
 *
 * The initial condition is u(t_0) = 1.0 and we use the parameters
 * p = [1.5, 1.0, 3.0, 1.0]. We compute the sensitivities of the
 * scalar cost function
 *     g(u_f, p) = 0.5*||1 - u(t_f, u_0, p)||^2
 * with respect to the parameters via adjoint sensitivity analysis:
 * backward integration of the adjoint ODE with terminal condition
 * dg/du(t_f), plus quadratures mu^T (df/dp) for dg/dp.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]

use cvodes_rs::cvodea::{
    CVodeAdjInit, CVodeB, CVodeCreateB, CVodeF, CVodeGetB, CVodeGetQuadB, CVodeInitB,
    CVodeQuadInitB, CVodeQuadSStolerancesB, CVodeSStolerancesB,
};
use cvodes_rs::cvodea_io::{CVodeSetQuadErrConB, CVodeSetUserDataB};
use cvodes_rs::cvodes::{CVodeCreate, CVodeFree, CVodeInit, CVodeSStolerances};
use cvodes_rs::cvodes_io::{CVodeSetMaxNumSteps, CVodeSetUserData};
use cvodes_rs::cvodes_ls::CVodeSetLinearSolver;
use cvodes_rs::sundials_utils::fmt_g;
use cvodes_rs::*;

/* Problem Constants */
const NEQ: i64 = 2; /* number of equations  */
const NP: i64 = 4; /* number of params     */
const T0: f64 = 0.0; /* initial time         */
const TF: f64 = 10.0; /* final time           */
const RTOL: f64 = 1.0e-10; /* relative tolerance   */
const ATOL: f64 = 1.0e-14; /* absolute tolerance   */
const STEPS: i64 = 5; /* checkpoint interval  */

const PARAMS: [f64; 4] = [1.5, 1.0, 3.0, 1.0];

/* Function to compute the ODE right-hand side */
fn lotka_volterra(_t: f64, uvec: &NVector, udotvec: &mut NVector, user_data: &mut UserData) -> i32 {
    let p = user_data.as_ref().unwrap().downcast_ref::<[f64; 4]>().unwrap();
    let u = &uvec.data;
    let udot = &mut udotvec.data;

    udot[0] = p[0] * u[0] - p[1] * u[0] * u[1];
    udot[1] = -p[2] * u[1] + p[3] * u[0] * u[1];

    0
}

/* Function to compute v^T (df/du) */
fn vjp(vvec: &NVector, Jvvec: &mut NVector, _t: f64, uvec: &NVector, user_data: &mut UserData) -> i32 {
    let p = user_data.as_ref().unwrap().downcast_ref::<[f64; 4]>().unwrap();
    let u = &uvec.data;
    let v = &vvec.data;
    let Jv = &mut Jvvec.data;

    Jv[0] = (p[0] - p[1] * u[1]) * v[0] + p[3] * u[1] * v[1];
    Jv[1] = -p[1] * u[0] * v[0] + (-p[2] + p[3] * u[0]) * v[1];

    0
}

/* Function to compute v^T (df/dp) */
fn parameter_vjp(
    vvec: &NVector,
    Jvvec: &mut NVector,
    _t: f64,
    uvec: &NVector,
    user_data: &mut UserData,
) -> i32 {
    /* C checks the user_data pointer is the params array */
    if user_data
        .as_ref()
        .and_then(|d| d.downcast_ref::<[f64; 4]>())
        .is_none()
    {
        return -1;
    }

    let u = &uvec.data;
    let v = &vvec.data;
    let Jv = &mut Jvvec.data;

    Jv[0] = u[0] * v[0];
    Jv[1] = -u[0] * u[1] * v[0];
    Jv[2] = -u[1] * v[1];
    Jv[3] = u[0] * u[1] * v[1];

    0
}

/* Gradient of the cost function w.r.t to u.
   The gradient w.r.t to p is zero since the cost function
   does not depend on the parameters. */
fn dgdu(uvec: &NVector, dgvec: &mut NVector) {
    let u = &uvec.data;
    let dg = &mut dgvec.data;

    dg[0] = -1.0 + u[0];
    dg[1] = -1.0 + u[1];
}

/* Function to compute the adjoint ODE right-hand side:
    -mu^T (df/du)
 */
fn adjoint_rhs(t: f64, uvec: &NVector, lvec: &NVector, ldotvec: &mut NVector, user_data: &mut UserData) -> i32 {
    vjp(lvec, ldotvec, t, uvec, user_data);
    /* C: N_VScale(-1.0, ldotvec, ldotvec) — aliased operands */
    ldotvec.scale_inplace(-1.0);
    0
}

/* Function to compute the quadrature right-hand side:
    mu^T (df/dp)
 */
fn quad_rhs(t: f64, uvec: &NVector, muvec: &NVector, qBdotvec: &mut NVector, user_dataB: &mut UserData) -> i32 {
    parameter_vjp(muvec, qBdotvec, t, uvec, user_dataB);
    0
}

/* Check if a SUNDIALS function returned a negative value */
fn check_retval(retval: i32, funcname: &str) -> bool {
    if retval < 0 {
        eprintln!("\nSUNDIALS_ERROR: {}() failed with retval = {}\n", funcname, retval);
        return true;
    }
    false
}

fn main() {
    let sunctx = SUNContext_Create();

    /* Allocate memory for the solution vector */
    let mut u = N_VNew_Serial(NEQ, &sunctx);

    /* Initialize the solution vector */
    N_VConst(1.0, &mut u);

    /* Set the tolerances */
    let reltol = RTOL;
    let abstol = ATOL;

    /* Create the CVODES object */
    let mut cvode_mem = CVodeCreate(CV_BDF, &sunctx);

    /* Initialize the CVODES solver */
    let retval = CVodeInit(&mut cvode_mem, lotka_volterra, T0, &u);
    if check_retval(retval, "CVodeInit") {
        return;
    }

    /* Set the user data */
    let retval = CVodeSetUserData(&mut cvode_mem, Some(Box::new(PARAMS)));
    if check_retval(retval, "CVodeSetUserData") {
        return;
    }

    /* Set the tolerances */
    let retval = CVodeSStolerances(&mut cvode_mem, reltol, abstol);
    if check_retval(retval, "CVodeSStolerances") {
        return;
    }

    let ls = SUNLinSol_SPGMR(&u, SUN_PREC_NONE, 3, &sunctx);

    let retval = CVodeSetLinearSolver(&mut cvode_mem, ls, None);
    if check_retval(retval, "CVodeSetLinearSolver") {
        return;
    }

    let retval = CVodeSetMaxNumSteps(&mut cvode_mem, 100000);
    if check_retval(retval, "CVodeSetMaxNumSteps") {
        return;
    }

    /* Initialize ASA */
    let retval = CVodeAdjInit(&mut cvode_mem, STEPS, CV_HERMITE);
    if check_retval(retval, "CVodeAdjInit") {
        return;
    }

    /* Integrate the ODE */
    let tout = TF;
    let mut t = 0.0;
    let mut ncheck: i32 = 0;
    let retval = CVodeF(&mut cvode_mem, tout, &mut u, &mut t, CV_NORMAL, &mut ncheck);
    if check_retval(retval, "CVode") {
        return;
    }

    /* Print the final solution */
    println!("Forward Solution at t = {}:", fmt_g(t, 0, 6));
    N_VPrint(&u);

    /* Allocate memory for the adjoint solution vector */
    let mut uB = N_VNew_Serial(NEQ, &sunctx);

    /* Allocate memory for the quadrature equations and initialize it to zero */
    let mut qB = N_VNew_Serial(NP, &sunctx);
    N_VConst(0.0, &mut qB);

    /* Initialize the adjoint solution vector */
    dgdu(&u, &mut uB);

    println!("Adjoint terminal condition:");
    N_VPrint(&uB);
    N_VPrint(&qB);

    /* Create the CVODES object for the backward problem */
    let mut which: i32 = 0;
    let _retval = CVodeCreateB(&mut cvode_mem, CV_BDF, &mut which);

    /* Initialize the CVODES solver for the backward problem */
    let retval = CVodeInitB(&mut cvode_mem, which, adjoint_rhs, TF, &uB);
    if check_retval(retval, "CVodeInitB") {
        return;
    }

    /* Set the user data for the backward problem */
    let retval = CVodeSetUserDataB(&mut cvode_mem, which, Some(Box::new(PARAMS)));
    if check_retval(retval, "CVodeSetUserDataB") {
        return;
    }

    /* Set the tolerances for the backward problem */
    let retval = CVodeSStolerancesB(&mut cvode_mem, which, reltol, abstol);
    if check_retval(retval, "CVodeSStolerancesB") {
        return;
    }

    /* Create the linear solver for the backward problem */
    let lsB = SUNLinSol_SPGMR(&uB, SUN_PREC_NONE, 3, &sunctx);

    let retval = CVodeSetLinearSolverB(&mut cvode_mem, which, lsB, None);
    if check_retval(retval, "CVodeSetLinearSolver") {
        return;
    }

    /* Call CVodeQuadInitB to allocate internal memory and initialize backward
       quadrature integration. This gives the sensitivities w.r.t. the
       parameters. */
    let retval = CVodeQuadInitB(&mut cvode_mem, which, quad_rhs, &qB);
    if check_retval(retval, "CVodeQuadInitB") {
        return;
    }

    /* Include the quadrature variables in the error control mechanism */
    let retval = CVodeSetQuadErrConB(&mut cvode_mem, which, SUNTRUE);
    if check_retval(retval, "CVodeSetQuadErrConB") {
        return;
    }

    /* Call CVodeQuadSStolerancesB to specify the scalar relative and
       absolute tolerances for the backward problem. */
    let retval = CVodeQuadSStolerancesB(&mut cvode_mem, which, reltol, abstol);
    if check_retval(retval, "CVodeQuadSStolerancesB") {
        return;
    }

    /* Integrate the adjoint ODE */
    let retval = CVodeB(&mut cvode_mem, T0, CV_NORMAL);
    if check_retval(retval, "CVodeB") {
        return;
    }

    /* Get the final adjoint solution */
    let retval = CVodeGetB(&mut cvode_mem, which, &mut t, &mut uB);
    if check_retval(retval, "CVodeGetB") {
        return;
    }

    /* Call CVodeGetQuadB to get the quadrature solution vector after a
       successful return from CVodeB. */
    let retval = CVodeGetQuadB(&mut cvode_mem, which, &mut t, &mut qB);
    if check_retval(retval, "CVodeGetQuadB") {
        return;
    }

    /* dg/dp = -qB   (C: N_VScale(-1.0, qB, qB) — aliased operands) */
    qB.scale_inplace(-1.0);

    /* Print the final adjoint solution */
    println!("Adjoint Solution at t = {}:", fmt_g(t, 0, 6));
    N_VPrint(&uB);
    N_VPrint(&qB);

    /* Free memory */
    drop(u);
    drop(uB);
    drop(qB);
    CVodeFree(cvode_mem);
}
