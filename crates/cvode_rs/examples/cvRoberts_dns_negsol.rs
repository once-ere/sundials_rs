/* -----------------------------------------------------------------
 * Translated from examples/cvode/serial/cvRoberts_dns_negsol.c
 * (CVODE 7.7.0)
 *
 * Modification of the CVODE example cvRoberts_dns to illustrate
 * the treatment of unphysical solution components through the RHS
 * function return retval.
 *
 * Note that, to make possible negative solution components, the
 * absolute tolerances had to be loosened a bit from their values
 * in cvRoberts_dns.
 *
 * The following is a simple example problem, with the coding
 * needed for its solution by CVODE. The problem is from
 * chemical kinetics, and consists of the following three rate
 * equations:
 *    dy1/dt = -.04*y1 + 1.e4*y2*y3
 *    dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*(y2)^2
 *    dy3/dt = 3.e7*(y2)^2
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1.0, y2 = y3 = 0. The problem is stiff.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use cvode_rs::sundials_utils::fmt_e;
use cvode_rs::*;

/* Problem Constants */
const NEQ: usize = 3; /* number of equations  */
const Y1: f64 = 1.0; /* initial y components */
const Y2: f64 = 0.0;
const Y3: f64 = 0.0;
const RTOL: f64 = 1.0e-4; /* scalar relative tolerance            */
const ATOL1: f64 = 1.0e-7; /* vector absolute tolerance components */
const ATOL2: f64 = 1.0e-13;
const ATOL3: f64 = 1.0e-5;
const T0: f64 = 0.0; /* initial time           */
const T1: f64 = 0.4; /* first output time      */
const TMULT: f64 = 10.0; /* output time factor     */
const NOUT: i32 = 14; /* number of output times */

/*
 * f routine. Compute function f(t,y).
 */
fn f(_t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let check_negative = *user_data
        .as_mut()
        .unwrap()
        .downcast_mut::<bool>()
        .unwrap();

    let y1 = y.data[0];
    let y2 = y.data[1];
    let y3 = y.data[2];

    if check_negative && (y1 < 0.0 || y2 < 0.0 || y3 < 0.0) {
        return 1;
    }

    let yd1 = -0.04 * y1 + 1.0e4 * y2 * y3;
    ydot.data[0] = yd1;
    let yd3 = 3.0e7 * y2 * y2;
    ydot.data[2] = yd3;
    ydot.data[1] = -yd1 - yd3;

    0
}

/*
 * Private helper functions
 */
fn PrintOutput(t: f64, y1: f64, y2: f64, y3: f64) {
    println!(
        "At t = {}      y ={}  {}  {}",
        fmt_e(t, 0, 4),
        fmt_e(y1, 14, 6),
        fmt_e(y2, 14, 6),
        fmt_e(y3, 14, 6)
    );
}

fn PrintFinalStats(cvode_mem: &mut CVodeMem) {
    let (mut nst, mut nfe, mut nsetups, mut nje, mut nfeLS) = (0i64, 0i64, 0i64, 0i64, 0i64);
    let (mut nni, mut nnf, mut ncfn, mut netf) = (0i64, 0i64, 0i64, 0i64);

    let mut retval = CVodeGetNumSteps(cvode_mem, &mut nst);
    check_retval(retval, "CVodeGetNumSteps");
    retval = CVodeGetNumRhsEvals(cvode_mem, &mut nfe);
    check_retval(retval, "CVodeGetNumRhsEvals");
    retval = CVodeGetNumLinSolvSetups(cvode_mem, &mut nsetups);
    check_retval(retval, "CVodeGetNumLinSolvSetups");
    retval = CVodeGetNumErrTestFails(cvode_mem, &mut netf);
    check_retval(retval, "CVodeGetNumErrTestFails");
    retval = CVodeGetNumNonlinSolvIters(cvode_mem, &mut nni);
    check_retval(retval, "CVodeGetNumNonlinSolvIters");
    retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &mut nnf);
    check_retval(retval, "CVodeGetNumNonlinSolvConvFails");
    retval = CVodeGetNumStepSolveFails(cvode_mem, &mut ncfn);
    check_retval(retval, "CVodeGetNumStepSolveFails");

    retval = CVodeGetNumJacEvals(cvode_mem, &mut nje);
    check_retval(retval, "CVodeGetNumJacEvals");
    retval = CVodeGetNumLinRhsEvals(cvode_mem, &mut nfeLS);
    check_retval(retval, "CVodeGetNumLinRhsEvals");

    println!("\nFinal Statistics:");
    println!(
        "nst = {:<6} nfe = {:<6} nsetups = {:<6} nfeLS = {:<6} nje = {}",
        nst, nfe, nsetups, nfeLS, nje
    );
    println!(
        "nni = {:<6} nnf = {:<6} netf = {:<6}    ncfn = {:<6}\n",
        nni, nnf, netf, ncfn
    );
}

fn check_retval(retval: i32, funcname: &str) -> bool {
    if retval < 0 {
        eprintln!("\nSUNDIALS_ERROR: {}() failed with retval = {}\n", funcname, retval);
        return true;
    }
    false
}

fn main() {
    /* Create the SUNDIALS context */
    let sunctx = SUNContext_Create();

    /* Initial conditions */
    let mut y = N_VNew_Serial(NEQ as i64, &sunctx);

    /* Initialize y */
    y.data[0] = Y1;
    y.data[1] = Y2;
    y.data[2] = Y3;

    /* Set the vector absolute tolerance */
    let mut abstol = N_VNew_Serial(NEQ as i64, &sunctx);
    abstol.data[0] = ATOL1;
    abstol.data[1] = ATOL2;
    abstol.data[2] = ATOL3;

    /* Call CVodeCreate to create the solver memory and specify the
     * Backward Differentiation Formula */
    let mut cvode_mem = CVodeCreate(CV_BDF, &sunctx);

    /* Call CVodeInit to initialize the integrator memory and specify the
     * user's right hand side function in y'=f(t,y), the initial time T0, and
     * the initial dependent variable vector y. */
    let mut retval = CVodeInit(&mut cvode_mem, f, T0, &y);
    if check_retval(retval, "CVodeInit") {
        std::process::exit(1);
    }

    /* Call CVodeSVtolerances to specify the scalar relative tolerance
     * and vector absolute tolerances */
    retval = CVodeSVtolerances(&mut cvode_mem, RTOL, &abstol);
    if check_retval(retval, "CVodeSVtolerances") {
        std::process::exit(1);
    }

    /* Call CVodeSetUserData to pass the check negative retval as user data */
    retval = CVodeSetUserData(&mut cvode_mem, Some(Box::new(SUNFALSE)));
    if check_retval(retval, "CVodeSetUserData") {
        std::process::exit(1);
    }

    /* Create dense SUNMatrix for use in linear solves */
    let a = SUNDenseMatrix(NEQ as i64, NEQ as i64, &sunctx);

    /* Create dense SUNLinearSolver object for use by CVode */
    let ls = SUNLinSol_Dense(&y, &a, &sunctx);

    /* Attach the matrix and linear solver */
    retval = CVodeSetLinearSolver(&mut cvode_mem, ls, Some(a));
    if check_retval(retval, "CVodeSetLinearSolver") {
        std::process::exit(1);
    }

    /* Case 1: ignore negative solution components */
    println!("Ignore negative solution components\n");
    *CVodeGetUserData(&mut cvode_mem)
        .as_mut()
        .unwrap()
        .downcast_mut::<bool>()
        .unwrap() = SUNFALSE;
    /* In loop, call CVode in CV_NORMAL mode */
    let mut iout = 0;
    let mut tout = T1;
    let mut t = 0.0;
    loop {
        let _ = CVode(&mut cvode_mem, tout, &mut y, &mut t, CV_NORMAL);
        PrintOutput(t, y.data[0], y.data[1], y.data[2]);
        iout += 1;
        tout *= TMULT;
        if iout == NOUT {
            break;
        }
    }
    /* Print some final statistics */
    PrintFinalStats(&mut cvode_mem);

    /* Case 2: intercept negative solution components */
    println!("Intercept negative solution components\n");
    *CVodeGetUserData(&mut cvode_mem)
        .as_mut()
        .unwrap()
        .downcast_mut::<bool>()
        .unwrap() = SUNTRUE;
    /* Reinitialize solver */
    y.data[0] = Y1;
    y.data[1] = Y2;
    y.data[2] = Y3;
    let _ = CVodeReInit(&mut cvode_mem, T0, &y);
    /* In loop, call CVode in CV_NORMAL mode */
    iout = 0;
    tout = T1;
    loop {
        let _ = CVode(&mut cvode_mem, tout, &mut y, &mut t, CV_NORMAL);
        PrintOutput(t, y.data[0], y.data[1], y.data[2]);
        iout += 1;
        tout *= TMULT;
        if iout == NOUT {
            break;
        }
    }
    /* Print some final statistics */
    PrintFinalStats(&mut cvode_mem);

    /* Free memory (RAII) */
    CVodeFree(cvode_mem);
}
