/* -----------------------------------------------------------------
 * Translated from examples/cvode/serial/cvDisc_dns.c (CVODE 7.7.0)
 *
 * Simple 1D example to illustrate integrating over discontinuities:
 *
 * A) Discontinuity in solution
 *       y' = -y   ; y(0) = 1    ; t = [0,1]
 *       y' = -y   ; y(1) = 1    ; t = [1,2]
 *
 * B) Discontinuity in RHS (y')
 *       y' = -y   ; y(0) = 1    ; t = [0,1]
 *       z' = -5*z ; z(1) = y(1) ; t = [1,2]
 *    This case is solved twice, first by explicitly treating the
 *    discontinuity point and secondly by letting the integrator
 *    deal with the discontinuity.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use cvode_rs::sundials_utils::fmt_e;
use cvode_rs::*;

/* Problem Constants */
const NEQ: usize = 1; /* number of equations */

const RHS1: i32 = 1;
const RHS2: i32 = 2;

/*
 * RHS function
 * The form of the RHS function is controlled by the flag passed as f_data:
 *   flag = RHS1 -> y' = -y
 *   flag = RHS2 -> y' = -5*y
 */
fn f(_t: f64, y: &NVector, ydot: &mut NVector, f_data: &mut UserData) -> i32 {
    let flag = *f_data.as_mut().unwrap().downcast_mut::<i32>().unwrap();

    match flag {
        RHS1 => ydot.data[0] = -y.data[0],
        RHS2 => ydot.data[0] = -5.0 * y.data[0],
        _ => {}
    }

    0
}

/* Update the RHS flag stored in the integrator's user data
   (the C example shares `flag` with the solver by address). */
fn set_flag(cvode_mem: &mut CVodeMem, flag: i32) {
    *CVodeGetUserData(cvode_mem)
        .as_mut()
        .unwrap()
        .downcast_mut::<i32>()
        .unwrap() = flag;
}

fn check_retval(retval: i32, funcname: &str) -> bool {
    if retval < 0 {
        eprintln!("\nSUNDIALS_ERROR: {}() failed with flag = {}\n", funcname, retval);
        return true;
    }
    false
}

fn PrintSolution(t: f64, y0: f64) {
    println!("{}  {}", fmt_e(t, 12, 8), fmt_e(y0, 12, 8));
}

fn main() {
    let reltol = 1.0e-3;
    let abstol = 1.0e-4;

    let t0 = 0.0;
    let t1 = 1.0;
    let t2 = 2.0;

    /* Create the SUNDIALS context */
    let sunctx = SUNContext_Create();

    /* Allocate the vector of initial conditions */
    let mut y = N_VNew_Serial(NEQ as i64, &sunctx);

    /* Set initial condition */
    y.data[0] = 1.0;

    /*
     * ------------------------------------------------------------
     *  Shared initialization and setup
     * ------------------------------------------------------------
     */

    /* Call CVodeCreate to create CVODE memory block and specify the
     * Backward Differentiation Formula */
    let mut cvode_mem = CVodeCreate(CV_BDF, &sunctx);

    /* Call CVodeInit to initialize integrator memory and specify the
     * user's right hand side function y'=f(t,y), the initial time T0
     * and the initial condition vector y. */
    let mut retval = CVodeInit(&mut cvode_mem, f, t0, &y);
    if check_retval(retval, "CVodeInit") {
        std::process::exit(1);
    }

    /* Call CVodeSStolerances to specify integration tolerances,
     * specifically the scalar relative and absolute tolerance. */
    retval = CVodeSStolerances(&mut cvode_mem, reltol, abstol);
    if check_retval(retval, "CVodeSStolerances") {
        std::process::exit(1);
    }

    /* Provide RHS flag as user data which can be accessed in user
     * provided routines */
    retval = CVodeSetUserData(&mut cvode_mem, Some(Box::new(RHS1)));
    if check_retval(retval, "CVodeSetUserData") {
        std::process::exit(1);
    }

    /* Create dense SUNMatrix for use in linear solver */
    let a = SUNDenseMatrix(NEQ as i64, NEQ as i64, &sunctx);

    /* Create dense linear solver for use by CVode */
    let ls = SUNLinSol_Dense(&y, &a, &sunctx);

    /* Attach the linear solver and matrix to CVode by calling
     * CVodeSetLinearSolver */
    retval = CVodeSetLinearSolver(&mut cvode_mem, ls, Some(a));
    if check_retval(retval, "CVodeSetLinearSolver") {
        std::process::exit(1);
    }

    /*
     * ---------------------------------------------------------------
     * Discontinuity in the solution
     *
     * 1) Integrate to the discontinuity
     * 2) Integrate from the discontinuity
     * ---------------------------------------------------------------
     */

    /* ---- Integrate to the discontinuity */

    println!("\nDiscontinuity in solution\n");

    /* set TSTOP (max time solution proceeds to) - this is not required */
    retval = CVodeSetStopTime(&mut cvode_mem, t1);
    if check_retval(retval, "CVodeSetStopTime") {
        std::process::exit(1);
    }

    set_flag(&mut cvode_mem, RHS1); /* use -y for RHS */
    let mut t = t0; /* set the integrator start time */

    PrintSolution(t, y.data[0]);
    while t < t1 {
        /* advance solver just one internal step */
        retval = CVode(&mut cvode_mem, t1, &mut y, &mut t, CV_ONE_STEP);
        if check_retval(retval, "CVode") {
            std::process::exit(1);
        }
        PrintSolution(t, y.data[0]);
    }
    /* Get the number of steps the solver took to get to the discont. */
    let mut nst1 = 0i64;
    retval = CVodeGetNumSteps(&mut cvode_mem, &mut nst1);
    if check_retval(retval, "CvodeGetNumSteps") {
        std::process::exit(1);
    }

    /* ---- Integrate from the discontinuity */

    /* Include discontinuity */
    y.data[0] = 1.0;

    /* Reinitialize the solver */
    retval = CVodeReInit(&mut cvode_mem, t1, &y);
    if check_retval(retval, "CVodeReInit") {
        std::process::exit(1);
    }

    /* set TSTOP (max time solution proceeds to) - this is not required */
    retval = CVodeSetStopTime(&mut cvode_mem, t2);
    if check_retval(retval, "CVodeSetStopTime") {
        std::process::exit(1);
    }

    set_flag(&mut cvode_mem, RHS1); /* use -y for RHS */
    t = t1; /* set the integrator start time */

    PrintSolution(t, y.data[0]);

    while t < t2 {
        /* advance solver just one internal step */
        retval = CVode(&mut cvode_mem, t2, &mut y, &mut t, CV_ONE_STEP);
        if check_retval(retval, "CVode") {
            std::process::exit(1);
        }
        PrintSolution(t, y.data[0]);
    }

    /* Get the number of steps the solver took after the discont. */
    let mut nst2 = 0i64;
    retval = CVodeGetNumSteps(&mut cvode_mem, &mut nst2);
    if check_retval(retval, "CvodeGetNumSteps") {
        std::process::exit(1);
    }

    /* Print statistics */
    let mut nst = nst1 + nst2;
    println!("\nNumber of steps: {} + {} = {}", nst1, nst2, nst);

    /*
     * ---------------------------------------------------------------
     * Discontinuity in RHS: Case 1 - explicit treatment
     * Note that it is not required to set TSTOP, but without it
     * we would have to find y(t1) to reinitialize the solver.
     * ---------------------------------------------------------------
     */

    println!("\nDiscontinuity in RHS: Case 1 - explicit treatment\n");

    /* Set initial condition */
    y.data[0] = 1.0;

    /* Reinitialize the solver. CVodeReInit does not reallocate memory
     * so it can only be used when the new problem size is the same as
     * the problem size when CVodeCreate was called. */
    retval = CVodeReInit(&mut cvode_mem, t0, &y);
    if check_retval(retval, "CVodeReInit") {
        std::process::exit(1);
    }

    /* ---- Integrate to the discontinuity */

    /* Set TSTOP (max time solution proceeds to) to location of discont. */
    retval = CVodeSetStopTime(&mut cvode_mem, t1);
    if check_retval(retval, "CVodeSetStopTime") {
        std::process::exit(1);
    }

    set_flag(&mut cvode_mem, RHS1); /* use -y for RHS */
    t = t0; /* set the integrator start time */

    PrintSolution(t, y.data[0]);
    while t < t1 {
        /* advance solver just one internal step */
        retval = CVode(&mut cvode_mem, t1, &mut y, &mut t, CV_ONE_STEP);
        if check_retval(retval, "CVode") {
            std::process::exit(1);
        }
        PrintSolution(t, y.data[0]);
    }

    /* Get the number of steps the solver took to get to the discont. */
    retval = CVodeGetNumSteps(&mut cvode_mem, &mut nst1);
    if check_retval(retval, "CvodeGetNumSteps") {
        std::process::exit(1);
    }

    /* If TSTOP was not set, we'd need to find y(t1): */
    /* CVodeGetDky(cvode_mem, t1, 0, y); */

    /* ---- Integrate from the discontinuity */

    /* Reinitialize solver */
    retval = CVodeReInit(&mut cvode_mem, t1, &y);

    /* set TSTOP (max time solution proceeds to) - this is not required */
    retval = CVodeSetStopTime(&mut cvode_mem, t2);
    if check_retval(retval, "CVodeSetStopTime") {
        std::process::exit(1);
    }

    set_flag(&mut cvode_mem, RHS2); /* use -5y for RHS */
    t = t1; /* set the integrator start time */

    PrintSolution(t, y.data[0]);

    while t < t2 {
        /* advance solver just one internal step */
        retval = CVode(&mut cvode_mem, t2, &mut y, &mut t, CV_ONE_STEP);
        if check_retval(retval, "CVode") {
            std::process::exit(1);
        }
        PrintSolution(t, y.data[0]);
    }

    /* Get the number of steps the solver took after the discont. */
    retval = CVodeGetNumSteps(&mut cvode_mem, &mut nst2);
    if check_retval(retval, "CvodeGetNumSteps") {
        std::process::exit(1);
    }

    /* Print statistics */
    nst = nst1 + nst2;
    println!("\nNumber of steps: {} + {} = {}", nst1, nst2, nst);

    /*
     * ---------------------------------------------------------------
     * Discontinuity in RHS: Case 2 - let CVODE deal with it
     * Note that here we MUST set TSTOP to ensure that the
     * change in the RHS happens at the appropriate time
     * ---------------------------------------------------------------
     */

    println!("\nDiscontinuity in RHS: Case 2 - let CVODE deal with it\n");

    /* Set initial condition */
    y.data[0] = 1.0;

    /* Reinitialize the solver. CVodeReInit does not reallocate memory
     * so it can only be used when the new problem size is the same as
     * the problem size when CVodeCreate was called. */
    retval = CVodeReInit(&mut cvode_mem, t0, &y);
    if check_retval(retval, "CVodeReInit") {
        std::process::exit(1);
    }

    /* ---- Integrate to the discontinuity */

    /* Set TSTOP (max time solution proceeds to) to location of discont. */
    retval = CVodeSetStopTime(&mut cvode_mem, t1);
    if check_retval(retval, "CVodeSetStopTime") {
        std::process::exit(1);
    }

    set_flag(&mut cvode_mem, RHS1); /* use -y for RHS */
    t = t0; /* set the integrator start time */

    PrintSolution(t, y.data[0]);
    while t < t1 {
        /* advance solver just one internal step */
        retval = CVode(&mut cvode_mem, t1, &mut y, &mut t, CV_ONE_STEP);
        if check_retval(retval, "CVode") {
            std::process::exit(1);
        }
        PrintSolution(t, y.data[0]);
    }

    /* Get the number of steps the solver took to get to the discont. */
    retval = CVodeGetNumSteps(&mut cvode_mem, &mut nst1);
    if check_retval(retval, "CvodeGetNumSteps") {
        std::process::exit(1);
    }

    /* ---- Integrate from the discontinuity */

    /* set TSTOP (max time solution proceeds to) - this is not required */
    retval = CVodeSetStopTime(&mut cvode_mem, t2);
    if check_retval(retval, "CVodeSetStopTime") {
        std::process::exit(1);
    }

    set_flag(&mut cvode_mem, RHS2); /* use -5y for RHS */
    t = t1; /* set the integrator start time */

    PrintSolution(t, y.data[0]);

    while t < t2 {
        /* advance solver just one internal step */
        retval = CVode(&mut cvode_mem, t2, &mut y, &mut t, CV_ONE_STEP);
        if check_retval(retval, "CVode") {
            std::process::exit(1);
        }
        PrintSolution(t, y.data[0]);
    }

    /* Get the number of steps the solver took after the discont. */
    retval = CVodeGetNumSteps(&mut cvode_mem, &mut nst);
    if check_retval(retval, "CvodeGetNumSteps") {
        std::process::exit(1);
    }

    /* Print statistics */
    nst2 = nst - nst1;
    println!("\nNumber of steps: {} + {} = {}", nst1, nst2, nst);

    /* Free memory (RAII) */
    CVodeFree(cvode_mem);
}
