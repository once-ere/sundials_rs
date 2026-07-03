/* -----------------------------------------------------------------
 * Translated from examples/cvode/serial/cvRocket_dns.c (CVODE 7.7.0)
 *
 * Example problem:
 *
 * The following is a simple example problem, with the coding needed
 * for its solution by CVODE. The problem is a simpliflied model of a
 * rocket, ascending vertically, with mass decreasing over time. The
 * system (of size 2) is given by
 *    y_1 = rocket height H, y_1(0) = 0,
 *    y_2 = rocket velocity v, y_2(0) = 0,
 *    dH/dt = v,
 *    dv/dt = a(t,v).
 * The upward acceleration a(t,v) is given by
 *    a(t,v) = F/(M_r + M_f) - Dv - g,
 * where F = engine thrust force (constant) M_r = rocket mass without
 * fuel, M_f = fuel mass = M_f0 - r*t, r = fuel burn rate,
 * D = drag coefficient, g = gravitational acceleration.
 * The engine force is reset to 0 when the fuel mass reaches 0, or
 * when H reaches a preset height H_c, whichever happens first.
 * Rootfinding is used to locate the time at which M_f = 0 or H =
 * H_c, and also the time at which the rocket reaches its maximum
 * height, given by the condition v = 0, t > 0.
 *
 * The problem is solved with the BDF method and Dense linear solver.
 *
 * Run statistics (optional outputs) are printed at the end.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use cvode_rs::sundials_utils::fmt_e;
use cvode_rs::*;

/* Problem Constants */

const NEQ: usize = 2; /* number of equations  */

const Force: f64 = 2200.0; /* engine force */
const massr: f64 = 10.0; /* rocket mass (empty) */
const massf0: f64 = 1.0; /* initial fuel mass */
const brate: f64 = 0.1; /* fuel burn rate */
const Drag: f64 = 0.3; /* Drag coefficient */
const grav: f64 = 32.0; /* acceleration due to gravity */
const Hcut: f64 = 4000.0; /* height of engine cutoff */

const Y1: f64 = 0.0; /* initial y components */
const Y2: f64 = 0.0;
const RTOL: f64 = 1.0e-5; /* scalar relative tolerance            */
const ATOL1: f64 = 1.0e-2; /* vector absolute tolerance components */
const ATOL2: f64 = 1.0e-1;
const T0: f64 = 0.0; /* initial time           */
const T1: f64 = 1.0; /* first output time      */
const TINC: f64 = 1.0; /* output time increment  */
const NOUT: i32 = 70; /* number of output times */

const ZERO: f64 = 0.0;

/*
 * f routine. Compute function f(t,y).
 */
fn f(t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let engine_on = *user_data
        .as_mut()
        .unwrap()
        .downcast_mut::<bool>()
        .unwrap();

    let v = y.data[1];
    ydot.data[0] = v;

    let acc = if engine_on {
        Force / (massr + massf0 - brate * t)
    } else {
        ZERO
    };

    ydot.data[1] = acc - Drag * v - grav;

    0
}

/*
 * Jacobian routine. Compute J(t,y) = df/dy.
 */
fn Jac(
    _t: f64,
    _y: &NVector,
    _fy: &NVector,
    j: &mut SUNMatrix,
    _user_data: &mut UserData,
    _tmp1: &mut NVector,
    _tmp2: &mut NVector,
    _tmp3: &mut NVector,
) -> i32 {
    let jm = match j {
        SUNMatrix::Dense(m) => m,
        _ => return -1,
    };
    /* C writes column-major Jdata[1] and Jdata[3] */
    jm.set(1, 0, 1.0);
    jm.set(1, 1, -Drag);

    0
}

/*
 * g routine. Compute functions g_i(t,y).
 */
fn g(t: f64, y: &NVector, gout: &mut [f64], user_data: &mut UserData) -> i32 {
    let engine_on = *user_data
        .as_mut()
        .unwrap()
        .downcast_mut::<bool>()
        .unwrap();

    if engine_on {
        gout[0] = massf0 - brate * t;
        let h = y.data[0];
        gout[1] = h - Hcut;
    } else {
        let v = y.data[1];
        gout[0] = v;
    }

    0
}

/*
 * Private helper functions
 */
fn PrintOutput(t: f64, y1: f64, y2: f64) {
    println!(
        "At t = {}      y ={}  {}",
        fmt_e(t, 0, 4),
        fmt_e(y1, 14, 6),
        fmt_e(y2, 14, 6)
    );
}

fn PrintRootInfo(root_f1: i32, root_f2: i32, numroot: i32) {
    if numroot == 2 {
        println!("    rootsfound[] = {:3} {:3}", root_f1, root_f2);
    }
    if numroot == 1 {
        println!("    rootsfound[] = {:3}", root_f1);
    }
}

/*
 * Get and print some final statistics
 */
fn PrintFinalStats(cvode_mem: &mut CVodeMem) {
    let (mut nst, mut nfe, mut nsetups, mut nje, mut nfeLS) = (0i64, 0i64, 0i64, 0i64, 0i64);
    let (mut nni, mut ncfn, mut netf, mut nge) = (0i64, 0i64, 0i64, 0i64);

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
    retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &mut ncfn);
    check_retval(retval, "CVodeGetNumNonlinSolvConvFails");

    retval = CVodeGetNumJacEvals(cvode_mem, &mut nje);
    check_retval(retval, "CVodeGetNumJacEvals");
    retval = CVodeGetNumLinRhsEvals(cvode_mem, &mut nfeLS);
    check_retval(retval, "CVodeGetNumLinRhsEvals");

    retval = CVodeGetNumGEvals(cvode_mem, &mut nge);
    check_retval(retval, "CVodeGetNumGEvals");

    println!("\nFinal Statistics:");
    /* C: "nje = % ld" — the space flag prints a blank before positive values */
    println!(
        "nst = {:<6} nfe  = {:<6} nsetups = {:<6} nfeLS = {:<6} nje = {}{}",
        nst,
        nfe,
        nsetups,
        nfeLS,
        if nje >= 0 { " " } else { "" },
        nje
    );
    /* C: "... nge = %ld\n \n" */
    println!("nni = {:<6} ncfn = {:<6} netf = {:<6} nge = {}", nni, ncfn, netf, nge);
    println!(" ");
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

    /* Create serial vector of length NEQ for I.C. and abstol */
    let mut y = N_VNew_Serial(NEQ as i64, &sunctx);
    let mut abstol = N_VNew_Serial(NEQ as i64, &sunctx);

    /* Initialize y */
    y.data[0] = Y1;
    y.data[1] = Y2;

    /* Set the scalar relative tolerance */
    let reltol = RTOL;
    /* Set the vector absolute tolerance */
    abstol.data[0] = ATOL1;
    abstol.data[1] = ATOL2;

    /* Call CVodeCreate to create the solver memory and specify the
     * Backward Differentiation Formula */
    let mut cvode_mem = CVodeCreate(CV_BDF, &sunctx);

    /* Call CVodeInit to initialize the integrator memory and specify the
     * right-hand side function in y'=f(t,y), the initial time T0, and the
     * initial dependent variable vector y. */
    let mut retval = CVodeInit(&mut cvode_mem, f, T0, &y);
    if check_retval(retval, "CVodeInit") {
        std::process::exit(1);
    }

    /* Call CVodeSVtolerances to specify the scalar relative tolerance and
     * vector absolute tolerances */
    retval = CVodeSVtolerances(&mut cvode_mem, reltol, &abstol);
    if check_retval(retval, "CVodeSVtolerances") {
        std::process::exit(1);
    }

    /* Provide engine_on as user data for use in f and g routines
       (set SUNTRUE before the first CVode call, as in the C main). */
    retval = CVodeSetUserData(&mut cvode_mem, Some(Box::new(SUNTRUE)));
    if check_retval(retval, "CVodeSetUserData") {
        std::process::exit(1);
    }

    /* Call CVodeRootInit to specify the root function g with 2 components */
    retval = CVodeRootInit(&mut cvode_mem, 2, Some(g));
    if check_retval(retval, "CVodeRootInit") {
        std::process::exit(1);
    }

    /* Create dense SUNMatrix for use in linear solves */
    let a = SUNDenseMatrix(NEQ as i64, NEQ as i64, &sunctx);

    /* Create dense SUNLinearSolver object for use by CVode */
    let ls = SUNLinSol_Dense(&y, &a, &sunctx);

    /* Call CVodeSetLinearSolver to attach the matrix and linear solver */
    retval = CVodeSetLinearSolver(&mut cvode_mem, ls, Some(a));
    if check_retval(retval, "CVodeSetLinearSolver") {
        std::process::exit(1);
    }

    /* Set the user-supplied Jacobian routine Jac */
    retval = CVodeSetJacFn(&mut cvode_mem, Some(Jac));
    if check_retval(retval, "CVodeSetJacFn") {
        std::process::exit(1);
    }

    /* In loop, call CVode, print results, check for root stops, and test for
       error.  On the first root return, restart with engine turned off. Break
       out of loop when NOUT preset output times have been reached, or when the
       returned value of H is negative.  */
    println!(" \nAccelerating rocket problem\n");

    let mut iout = 0;
    let mut tout = T1;
    let mut engine_on = SUNTRUE;
    let mut numroot = 2;
    let mut t = 0.0;
    let mut rootsfound = [0i32; 2];
    loop {
        retval = CVode(&mut cvode_mem, tout, &mut y, &mut t, CV_NORMAL);
        if check_retval(retval, "CVode") {
            std::process::exit(1);
        }

        PrintOutput(t, y.data[0], y.data[1]);

        if engine_on && (retval == CV_ROOT_RETURN) {
            /* engine cutoff */
            let retvalr = CVodeGetRootInfo(&mut cvode_mem, &mut rootsfound);
            if check_retval(retvalr, "CVodeGetRootInfo") {
                std::process::exit(1);
            }
            PrintRootInfo(rootsfound[0], rootsfound[1], numroot);
            engine_on = SUNFALSE;
            *CVodeGetUserData(&mut cvode_mem)
                .as_mut()
                .unwrap()
                .downcast_mut::<bool>()
                .unwrap() = SUNFALSE;
            numroot = 1;
            /* Call CVodeRootInit to specify the root function g with 1
               component.  NOTE: the C example reuses the `retval` variable
               here, so a successful RootInit/ReInit makes the loop's
               `retval == CV_SUCCESS` test advance iout/tout on this same
               iteration — replicated exactly (no shadowing). */
            retval = CVodeRootInit(&mut cvode_mem, 1, Some(g));
            if check_retval(retval, "CVodeRootInit") {
                std::process::exit(1);
            }
            /* Reinitialize the solver with current t and y values. */
            retval = CVodeReInit(&mut cvode_mem, t, &y);
            if check_retval(retval, "CVodeReInit") {
                std::process::exit(1);
            }
        } else if !engine_on && (retval == CV_ROOT_RETURN) {
            /* max.  height */
            let retvalr = CVodeGetRootInfo(&mut cvode_mem, &mut rootsfound);
            if check_retval(retvalr, "CVodeGetRootInfo") {
                std::process::exit(1);
            }
            PrintRootInfo(rootsfound[0], rootsfound[1], numroot);
        }

        if retval == CV_SUCCESS {
            iout += 1;
            tout += TINC;
        }

        if iout == NOUT {
            break;
        }
        if y.data[0] < ZERO {
            break;
        }
    }

    /* Print some final statistics */
    PrintFinalStats(&mut cvode_mem);

    /* Free memory (RAII) */
    CVodeFree(cvode_mem);

    std::process::exit(retval);
}
