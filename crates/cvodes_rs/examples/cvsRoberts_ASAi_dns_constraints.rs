/* -----------------------------------------------------------------
 * Translation of examples/cvodes/serial/cvsRoberts_ASAi_dns_constraints.c
 * (SUNDIALS 7.7.0).
 *
 * Adjoint sensitivity example problem (chemical kinetics):
 *
 * The following is a simple example problem, with the coding
 * needed for its solution by CVODES. The problem is from chemical
 * kinetics, and consists of the following three rate equations:
 *    dy1/dt = -p1*y1 + p2*y2*y3
 *    dy2/dt =  p1*y1 - p2*y2*y3 - p3*(y2)^2
 *    dy3/dt =  p3*(y2)^2
 * on the interval from t = 0.0 to t = 4.e7, with initial
 * conditions: y1 = 1.0, y2 = y3 = 0. The reaction rates are:
 * p1=0.04, p2=1e4, and p3=3e7. The problem is stiff.
 * The constraint y_i >= 0 is posed for all components (forward and
 * backward problems).
 *
 * Additionally, CVODES computes the sensitivities of
 *   G = int_t0^tB0 g(t,p,y) dt,  g(t,p,y) = y3
 * with respect to the problem parameters through backward
 * (adjoint) integration.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use cvodes_rs::cvodea::{
    CVodeAdjInit, CVodeB, CVodeCreateB, CVodeF, CVodeGetAdjY, CVodeGetB, CVodeGetQuadB,
    CVodeInitB, CVodeQuadInitB, CVodeQuadReInitB, CVodeQuadSStolerancesB, CVodeReInitB,
    CVodeSStolerancesB,
};
use cvodes_rs::cvodea_io::{
    CVodeGetAdjCVodeBmem, CVodeSetConstraintsB, CVodeSetQuadErrConB, CVodeSetUserDataB,
};
use cvodes_rs::cvodes::{
    CVodeCreate, CVodeFree, CVodeGetQuad, CVodeInit, CVodeQuadInit, CVodeQuadSStolerances,
    CVodeWFtolerances,
};
use cvodes_rs::cvodes_io::{
    CVodeGetNumSteps, CVodeSetConstraints, CVodeSetQuadErrCon, CVodeSetUserData,
};
use cvodes_rs::cvodes_ls::{
    CVodeSetJacFn, CVodeSetJacFnB, CVodeSetLinearSolver, CVodeSetLinearSolverB,
};
use cvodes_rs::sundials_utils::fmt_e;
use cvodes_rs::*;

/* Problem Constants */

const NEQ: i64 = 3; /* number of equations                  */

const RTOL: f64 = 1e-4; /* scalar relative tolerance            */

const ATOL1: f64 = 1e-4; /* vector absolute tolerance components */
const ATOL2: f64 = 1e-8;
const ATOL3: f64 = 1e-4;

const ATOLl: f64 = 1e-8; /* absolute tolerance for adjoint vars. */
const ATOLq: f64 = 1e-6; /* absolute tolerance for quadratures   */

const T0: f64 = 0.0; /* initial time                         */
const TOUT: f64 = 4e7; /* final time                           */

const TB1: f64 = 4e7; /* starting point for adjoint problem   */
const TB2: f64 = 50.0; /* starting point for adjoint problem   */
const TBout1: f64 = 40.0; /* intermediate t for adjoint problem   */

const STEPS: i64 = 150; /* number of steps between check points */

const NP: i64 = 3; /* number of problem parameters         */

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/* Type : UserData */
#[derive(Clone)]
struct RobData {
    p: [f64; 3],
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY CVODES
 *--------------------------------------------------------------------
 */

/*
 * f routine. Compute f(t,y).
 */
fn f(_t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let data = user_data.as_ref().unwrap().downcast_ref::<RobData>().unwrap();
    let y1 = y.data[0];
    let y2 = y.data[1];
    let y3 = y.data[2];
    let p1 = data.p[0];
    let p2 = data.p[1];
    let p3 = data.p[2];

    let yd1 = -p1 * y1 + p2 * y2 * y3;
    ydot.data[0] = yd1;
    let yd3 = p3 * y2 * y2;
    ydot.data[2] = yd3;
    ydot.data[1] = -yd1 - yd3;

    0
}

/*
 * Jacobian routine. Compute J(t,y).
 */
#[allow(clippy::too_many_arguments)]
fn jac(
    _t: f64,
    y: &NVector,
    _fy: &NVector,
    j: &mut SUNMatrix,
    user_data: &mut UserData,
    _tmp1: &mut NVector,
    _tmp2: &mut NVector,
    _tmp3: &mut NVector,
) -> i32 {
    let data = user_data.as_ref().unwrap().downcast_ref::<RobData>().unwrap();
    let y2 = y.data[1];
    let y3 = y.data[2];
    let p1 = data.p[0];
    let p2 = data.p[1];
    let p3 = data.p[2];

    if let SUNMatrix::Dense(dm) = j {
        /* IJth(J,i,j) = column-major dense storage */
        let m = 3usize;
        dm.data[0] = -p1; /* (1,1) */
        dm.data[m] = p2 * y3; /* (1,2) */
        dm.data[2 * m] = p2 * y2; /* (1,3) */
        dm.data[1] = p1; /* (2,1) */
        dm.data[1 + m] = -p2 * y3 - 2.0 * p3 * y2; /* (2,2) */
        dm.data[1 + 2 * m] = -p2 * y2; /* (2,3) */
        dm.data[2] = ZERO; /* (3,1) */
        dm.data[2 + m] = 2.0 * p3 * y2; /* (3,2) */
        dm.data[2 + 2 * m] = ZERO; /* (3,3) */
    }

    0
}

/*
 * fQ routine. Compute fQ(t,y).
 */
fn fQ(_t: f64, y: &NVector, qdot: &mut NVector, _user_data: &mut UserData) -> i32 {
    qdot.data[0] = y.data[2];

    0
}

/*
 * EwtSet function. Computes the error weights at the current solution.
 */
fn ewt(y: &NVector, w: &mut NVector, _user_data: &mut UserData) -> i32 {
    let rtol = RTOL;
    let atol = [ATOL1, ATOL2, ATOL3];

    for i in 0..3usize {
        let yy = y.data[i];
        let ww = rtol * yy.abs() + atol[i];
        if ww <= 0.0 {
            return -1;
        }
        w.data[i] = 1.0 / ww;
    }

    0
}

/*
 * fB routine. Compute fB(t,y,yB).
 */
fn fB(_t: f64, y: &NVector, yB: &NVector, yBdot: &mut NVector, user_dataB: &mut UserData) -> i32 {
    let data = user_dataB.as_ref().unwrap().downcast_ref::<RobData>().unwrap();

    /* The p vector */
    let p1 = data.p[0];
    let p2 = data.p[1];
    let p3 = data.p[2];

    /* The y vector */
    let y2 = y.data[1];
    let y3 = y.data[2];

    /* The lambda vector */
    let l1 = yB.data[0];
    let l2 = yB.data[1];
    let l3 = yB.data[2];

    /* Temporary variables */
    let l21 = l2 - l1;
    let l32 = l3 - l2;

    /* Load yBdot */
    yBdot.data[0] = -p1 * l21;
    yBdot.data[1] = p2 * y3 * l21 - 2.0 * p3 * y2 * l32;
    yBdot.data[2] = p2 * y2 * l21 - 1.0;

    0
}

/*
 * JacB routine. Compute JB(t,y,yB).
 */
#[allow(clippy::too_many_arguments)]
fn jacB(
    _t: f64,
    y: &NVector,
    _yB: &NVector,
    _fyB: &NVector,
    jB: &mut SUNMatrix,
    user_dataB: &mut UserData,
    _tmp1B: &mut NVector,
    _tmp2B: &mut NVector,
    _tmp3B: &mut NVector,
) -> i32 {
    let data = user_dataB.as_ref().unwrap().downcast_ref::<RobData>().unwrap();

    /* The p vector */
    let p1 = data.p[0];
    let p2 = data.p[1];
    let p3 = data.p[2];

    /* The y vector */
    let y2 = y.data[1];
    let y3 = y.data[2];

    /* Load JB */
    if let SUNMatrix::Dense(dm) = jB {
        let m = 3usize;
        dm.data[0] = p1; /* (1,1) */
        dm.data[m] = -p1; /* (1,2) */
        dm.data[2 * m] = ZERO; /* (1,3) */
        dm.data[1] = -p2 * y3; /* (2,1) */
        dm.data[1 + m] = p2 * y3 + 2.0 * p3 * y2; /* (2,2) */
        dm.data[1 + 2 * m] = -2.0 * p3 * y2; /* (2,3) */
        dm.data[2] = -p2 * y2; /* (3,1) */
        dm.data[2 + m] = p2 * y2; /* (3,2) */
        dm.data[2 + 2 * m] = ZERO; /* (3,3) */
    }

    0
}

/*
 * fQB routine. Compute integrand for quadratures
 */
fn fQB(_t: f64, y: &NVector, yB: &NVector, qBdot: &mut NVector, _user_dataB: &mut UserData) -> i32 {
    /* The y vector */
    let y1 = y.data[0];
    let y2 = y.data[1];
    let y3 = y.data[2];

    /* The lambda vector */
    let l1 = yB.data[0];
    let l2 = yB.data[1];
    let l3 = yB.data[2];

    /* Temporary variables */
    let l21 = l2 - l1;
    let l32 = l3 - l2;
    let y23 = y2 * y3;

    qBdot.data[0] = y1 * l21;
    qBdot.data[1] = -y23 * l21;
    qBdot.data[2] = y2 * y2 * l32;

    0
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * Print heading for backward integration
 */
fn print_head(tB0: f64) {
    println!("Backward integration from tB0 = {}\n", fmt_e(tB0, 12, 4));
}

/*
 * Print intermediate results during backward integration
 */
fn print_output1(time: f64, t: f64, y: &NVector, yB: &NVector) {
    println!("--------------------------------------------------------");
    println!("returned t: {}", fmt_e(time, 12, 4));
    println!("tout:       {}", fmt_e(t, 12, 4));
    println!(
        "lambda(t):  {} {} {}",
        fmt_e(yB.data[0], 12, 4),
        fmt_e(yB.data[1], 12, 4),
        fmt_e(yB.data[2], 12, 4)
    );
    println!(
        "y(t):       {} {} {}",
        fmt_e(y.data[0], 12, 4),
        fmt_e(y.data[1], 12, 4),
        fmt_e(y.data[2], 12, 4)
    );
    println!("--------------------------------------------------------\n");
}

/*
 * Print final results of backward integration
 */
fn print_output(tfinal: f64, y: &NVector, yB: &NVector, qB: &NVector) {
    println!("--------------------------------------------------------");
    println!("returned t: {}", fmt_e(tfinal, 12, 4));
    println!(
        "lambda(t0): {} {} {}",
        fmt_e(yB.data[0], 12, 4),
        fmt_e(yB.data[1], 12, 4),
        fmt_e(yB.data[2], 12, 4)
    );
    println!(
        "y(t0):      {} {} {}",
        fmt_e(y.data[0], 12, 4),
        fmt_e(y.data[1], 12, 4),
        fmt_e(y.data[2], 12, 4)
    );
    println!(
        "dG/dp:      {} {} {}",
        fmt_e(-qB.data[0], 12, 4),
        fmt_e(-qB.data[1], 12, 4),
        fmt_e(-qB.data[2], 12, 4)
    );
    println!("--------------------------------------------------------\n");
}

/* Check if a SUNDIALS function returned a negative value */
fn check_retval(retval: i32, funcname: &str) -> bool {
    if retval < 0 {
        eprintln!("\nSUNDIALS_ERROR: {}() failed with retval = {}\n", funcname, retval);
        return true;
    }
    false
}

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */
fn main() {
    /* Print problem description */
    println!("\nAdjoint Sensitivity Example for Chemical Kinetics");
    println!("-------------------------------------------------\n");
    println!("ODE: dy1/dt = -p1*y1 + p2*y2*y3");
    println!("     dy2/dt =  p1*y1 - p2*y2*y3 - p3*(y2)^2");
    println!("     dy3/dt =  p3*(y2)^2\n");
    println!("Find dG/dp for");
    println!("     G = int_t0^tB0 g(t,p,y) dt");
    println!("     g(t,p,y) = y3\n\n");

    /* User data structure */
    let data = RobData {
        p: [0.04, 1.0e4, 3.0e7],
    };

    /* Create the SUNDIALS simulation context that all SUNDIALS objects require */
    let sunctx = SUNContext_Create();

    /* Initialize y */
    let mut y = N_VNew_Serial(NEQ, &sunctx);
    y.data[0] = 1.0;
    y.data[1] = ZERO;
    y.data[2] = ZERO;

    /* Set constraints to all 1's for nonnegative solution values. */
    let mut constraints = N_VNew_Serial(NEQ, &sunctx);
    N_VConst(ONE, &mut constraints);

    /* Initialize q */
    let mut q = N_VNew_Serial(1, &sunctx);
    q.data[0] = ZERO;

    /* Set the scalar relative and absolute tolerances reltolQ and abstolQ */
    let reltolQ = RTOL;
    let abstolQ = ATOLq;

    /* Create and allocate CVODES memory for forward run */
    println!("Create and allocate CVODES memory for forward runs");

    /* Call CVodeCreate to create the solver memory and specify the
       Backward Differentiation Formula */
    let mut cvode_mem = CVodeCreate(CV_BDF, &sunctx);

    /* Call CVodeInit to initialize the integrator memory */
    let retval = CVodeInit(&mut cvode_mem, f, T0, &y);
    if check_retval(retval, "CVodeInit") {
        return;
    }

    /* Call CVodeWFtolerances to specify a user-supplied function ewt */
    let retval = CVodeWFtolerances(&mut cvode_mem, ewt);
    if check_retval(retval, "CVodeWFtolerances") {
        return;
    }

    /* Attach user data */
    let retval = CVodeSetUserData(&mut cvode_mem, Some(Box::new(data.clone())));
    if check_retval(retval, "CVodeSetUserData") {
        return;
    }

    /* Call CVodeSetConstraints to initialize constraints */
    let retval = CVodeSetConstraints(&mut cvode_mem, Some(&constraints));
    if check_retval(retval, "CVODESetConstraints") {
        return;
    }
    drop(constraints);

    /* Create dense SUNMatrix and dense SUNLinearSolver object */
    let a_mat = SUNDenseMatrix(NEQ, NEQ, &sunctx);
    let ls = SUNLinSol_Dense(&y, &a_mat, &sunctx);

    /* Attach the matrix and linear solver */
    let retval = CVodeSetLinearSolver(&mut cvode_mem, ls, Some(a_mat));
    if check_retval(retval, "CVodeSetLinearSolver") {
        return;
    }

    /* Set the user-supplied Jacobian routine Jac */
    let retval = CVodeSetJacFn(&mut cvode_mem, Some(jac));
    if check_retval(retval, "CVodeSetJacFn") {
        return;
    }

    /* Call CVodeQuadInit to allocate internal memory and initialize
       quadrature integration */
    let retval = CVodeQuadInit(&mut cvode_mem, fQ, &q);
    if check_retval(retval, "CVodeQuadInit") {
        return;
    }

    /* Call CVodeSetQuadErrCon to include the quadrature variables in the
       step size control mechanism */
    let retval = CVodeSetQuadErrCon(&mut cvode_mem, SUNTRUE);
    if check_retval(retval, "CVodeSetQuadErrCon") {
        return;
    }

    /* Call CVodeQuadSStolerances to specify scalar tolerances */
    let retval = CVodeQuadSStolerances(&mut cvode_mem, reltolQ, abstolQ);
    if check_retval(retval, "CVodeQuadSStolerances") {
        return;
    }

    /* Allocate global memory */

    /* Call CVodeAdjInit to update CVODES memory block by allocating the
       internal memory needed for backward integration. */
    let steps = STEPS; /* no. of integration steps between two consecutive checkpoints */
    let retval = CVodeAdjInit(&mut cvode_mem, steps, CV_HERMITE);
    if check_retval(retval, "CVodeAdjInit") {
        return;
    }

    /* Perform forward run */
    print!("Forward integration ... ");

    /* Call CVodeF to integrate the forward problem and save checkpointing data */
    let mut time = 0.0;
    let mut ncheck: i32 = 0;
    let retval = CVodeF(&mut cvode_mem, TOUT, &mut y, &mut time, CV_NORMAL, &mut ncheck);
    if check_retval(retval, "CVodeF") {
        return;
    }
    let mut nst: i64 = 0;
    let retval = CVodeGetNumSteps(&mut cvode_mem, &mut nst);
    if check_retval(retval, "CVodeGetNumSteps") {
        return;
    }

    println!("done ( nst = {} )", nst);
    println!("\nncheck = {}\n", ncheck);

    let retval = CVodeGetQuad(&cvode_mem, &mut time, &mut q);
    if check_retval(retval, "CVodeGetQuad") {
        return;
    }

    println!("--------------------------------------------------------");
    println!("G:          {} ", fmt_e(q.data[0], 12, 4));
    println!("--------------------------------------------------------\n");

    /* Initialize yB */
    let mut yB = N_VNew_Serial(NEQ, &sunctx);
    yB.data[0] = ZERO;
    yB.data[1] = ZERO;
    yB.data[2] = ZERO;

    /* Initialize qB */
    let mut qB = N_VNew_Serial(NP, &sunctx);
    qB.data[0] = ZERO;
    qB.data[1] = ZERO;
    qB.data[2] = ZERO;

    /* Set the scalar relative tolerance reltolB */
    let reltolB = RTOL;

    /* Set the scalar absolute tolerance abstolB */
    let abstolB = ATOLl;

    /* Set the scalar absolute tolerance abstolQB */
    let abstolQB = ATOLq;

    /* Set constraints to all 1's for nonnegative solution values. */
    let mut constraintsB = N_VNew_Serial(NEQ, &sunctx);
    N_VConst(ONE, &mut constraintsB);

    /* Create and allocate CVODES memory for backward run */
    println!("Create and allocate CVODES memory for backward run");

    /* Call CVodeCreateB to specify the solution method for the backward
       problem. */
    let mut indexB: i32 = 0;
    let retval = CVodeCreateB(&mut cvode_mem, CV_BDF, &mut indexB);
    if check_retval(retval, "CVodeCreateB") {
        return;
    }

    /* Call CVodeInitB to allocate internal memory and initialize the
       backward problem. */
    let retval = CVodeInitB(&mut cvode_mem, indexB, fB, TB1, &yB);
    if check_retval(retval, "CVodeInitB") {
        return;
    }

    /* Set the scalar relative and absolute tolerances. */
    let retval = CVodeSStolerancesB(&mut cvode_mem, indexB, reltolB, abstolB);
    if check_retval(retval, "CVodeSStolerancesB") {
        return;
    }

    /* Attach the user data for backward problem. */
    let retval = CVodeSetUserDataB(&mut cvode_mem, indexB, Some(Box::new(data.clone())));
    if check_retval(retval, "CVodeSetUserDataB") {
        return;
    }

    /* Call CVodeSetConstraintsB to initialize constraints */
    let retval = CVodeSetConstraintsB(&mut cvode_mem, indexB, Some(&constraintsB));
    if check_retval(retval, "CVodeSetConstraintsB") {
        return;
    }
    drop(constraintsB);

    /* Create dense SUNMatrix and dense SUNLinearSolver for the backward
       problem */
    let aB_mat = SUNDenseMatrix(NEQ, NEQ, &sunctx);
    let lsB = SUNLinSol_Dense(&yB, &aB_mat, &sunctx);

    /* Attach the matrix and linear solver */
    let retval = CVodeSetLinearSolverB(&mut cvode_mem, indexB, lsB, Some(aB_mat));
    if check_retval(retval, "CVodeSetLinearSolverB") {
        return;
    }

    /* Set the user-supplied Jacobian routine JacB */
    let retval = CVodeSetJacFnB(&mut cvode_mem, indexB, Some(jacB));
    if check_retval(retval, "CVodeSetJacFnB") {
        return;
    }

    /* Call CVodeQuadInitB to allocate internal memory and initialize
       backward quadrature integration. */
    let retval = CVodeQuadInitB(&mut cvode_mem, indexB, fQB, &qB);
    if check_retval(retval, "CVodeQuadInitB") {
        return;
    }

    /* Include the quadrature variables in the error control mechanism */
    let retval = CVodeSetQuadErrConB(&mut cvode_mem, indexB, SUNTRUE);
    if check_retval(retval, "CVodeSetQuadErrConB") {
        return;
    }

    /* Call CVodeQuadSStolerancesB to specify the scalar tolerances for the
       backward problem. */
    let retval = CVodeQuadSStolerancesB(&mut cvode_mem, indexB, reltolB, abstolQB);
    if check_retval(retval, "CVodeQuadSStolerancesB") {
        return;
    }

    /* Backward Integration */

    print_head(TB1);

    /* First get results at t = TBout1 */

    /* Call CVodeB to integrate the backward ODE problem. */
    let retval = CVodeB(&mut cvode_mem, TBout1, CV_NORMAL);
    if check_retval(retval, "CVodeB") {
        return;
    }

    /* Call CVodeGetB to get yB of the backward ODE problem. */
    let retval = CVodeGetB(&mut cvode_mem, indexB, &mut time, &mut yB);
    if check_retval(retval, "CVodeGetB") {
        return;
    }

    /* Call CVodeGetAdjY to get the interpolated value of the forward
       solution y during a backward integration. */
    let retval = CVodeGetAdjY(&mut cvode_mem, TBout1, &mut y);
    if check_retval(retval, "CVodeGetAdjY") {
        return;
    }

    print_output1(time, TBout1, &y, &yB);

    /* Then at t = T0 */

    let retval = CVodeB(&mut cvode_mem, T0, CV_NORMAL);
    if check_retval(retval, "CVodeB") {
        return;
    }
    let mut nstB: i64 = 0;
    if let Some(bmem) = CVodeGetAdjCVodeBmem(&mut cvode_mem, indexB) {
        CVodeGetNumSteps(bmem, &mut nstB);
    }
    println!("Done ( nst = {} )", nstB);

    let retval = CVodeGetB(&mut cvode_mem, indexB, &mut time, &mut yB);
    if check_retval(retval, "CVodeGetB") {
        return;
    }

    /* Call CVodeGetQuadB to get the quadrature solution vector after a
       successful return from CVodeB. */
    let retval = CVodeGetQuadB(&mut cvode_mem, indexB, &mut time, &mut qB);
    if check_retval(retval, "CVodeGetQuadB") {
        return;
    }

    let retval = CVodeGetAdjY(&mut cvode_mem, T0, &mut y);
    if check_retval(retval, "CVodeGetAdjY") {
        return;
    }

    print_output(time, &y, &yB, &qB);

    /* Reinitialize backward phase (new tB0) */

    yB.data[0] = ZERO;
    yB.data[1] = ZERO;
    yB.data[2] = ZERO;

    qB.data[0] = ZERO;
    qB.data[1] = ZERO;
    qB.data[2] = ZERO;

    println!("Re-initialize CVODES memory for backward run");

    let retval = CVodeReInitB(&mut cvode_mem, indexB, TB2, &yB);
    if check_retval(retval, "CVodeReInitB") {
        return;
    }

    let retval = CVodeQuadReInitB(&mut cvode_mem, indexB, &qB);
    if check_retval(retval, "CVodeQuadReInitB") {
        return;
    }

    print_head(TB2);

    /* First get results at t = TBout1 */

    let retval = CVodeB(&mut cvode_mem, TBout1, CV_NORMAL);
    if check_retval(retval, "CVodeB") {
        return;
    }

    let retval = CVodeGetB(&mut cvode_mem, indexB, &mut time, &mut yB);
    if check_retval(retval, "CVodeGetB") {
        return;
    }

    let retval = CVodeGetAdjY(&mut cvode_mem, TBout1, &mut y);
    if check_retval(retval, "CVodeGetAdjY") {
        return;
    }

    print_output1(time, TBout1, &y, &yB);

    /* Then at t = T0 */

    let retval = CVodeB(&mut cvode_mem, T0, CV_NORMAL);
    if check_retval(retval, "CVodeB") {
        return;
    }
    if let Some(bmem) = CVodeGetAdjCVodeBmem(&mut cvode_mem, indexB) {
        CVodeGetNumSteps(bmem, &mut nstB);
    }
    println!("Done ( nst = {} )", nstB);

    let retval = CVodeGetB(&mut cvode_mem, indexB, &mut time, &mut yB);
    if check_retval(retval, "CVodeGetB") {
        return;
    }

    let retval = CVodeGetQuadB(&mut cvode_mem, indexB, &mut time, &mut qB);
    if check_retval(retval, "CVodeGetQuadB") {
        return;
    }

    let retval = CVodeGetAdjY(&mut cvode_mem, T0, &mut y);
    if check_retval(retval, "CVodeGetAdjY") {
        return;
    }

    print_output(time, &y, &yB, &qB);

    /* Free memory */
    println!("Free memory\n");

    CVodeFree(cvode_mem);
    drop(y);
    drop(q);
    drop(yB);
    drop(qB);
}
