/* -----------------------------------------------------------------
 * Translated from examples/idas/serial/idasAnalytic_mels.c (IDAS 7.7.0)
 * Programmer(s): Daniel R. Reynolds @ UMBC
 *
 * The following is a simple example problem with analytical solution
 * adapted from example 10.2 of Ascher & Petzold, "Computer Methods
 * for Ordinary Differential Equations and Differential-Algebraic
 * Equations," SIAM, 1998, page 267:
 *    x1'(t) = (1-alpha)/(t-2)*x1 - x1 + (alpha-1)*x2 + 2*exp(t)
 *         0 = (t+2)*x1 - (t+2)*exp(t)
 * for t in the interval [0.0, 1.0], with initial condition:
 *    x1(0) = 1   and   x2(0) = -1/2.
 * The problem has true solution
 *    x1(t) = exp(t)  and  x2(t) = exp(t)/(t-2)
 *
 * This program solves the problem with IDAS using a custom
 * 'matrix-embedded' SUNLinearSolver. Output is printed every 0.1
 * units of time (10 total). Run statistics (optional outputs) are
 * printed at the end.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use idas_rs::sundials_utils::{fmt_e, fmt_f, fmt_g};
use idas_rs::*;

/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

/* System residual function:
      0 = (1-alpha)/(t-2)*x1 - x1 + (alpha-1)*x2 + 2*exp(t) - x1'(t)
      0 = (t+2)*x1 - (t+2)*exp(t)
*/
fn fres(t: f64, yy: &NVector, yp: &NVector, rr: &mut NVector, user_data: &mut UserData) -> i32 {
    let alpha = *user_data.as_mut().unwrap().downcast_mut::<f64>().unwrap();
    let x1 = yy.data[0];
    let x2 = yy.data[1];
    let x1p = yp.data[0];
    let ONE = 1.0;
    let TWO = 2.0;

    rr.data[0] =
        (ONE - alpha) / (t - TWO) * x1 - x1 + (alpha - ONE) * x2 + TWO * t.exp() - x1p;
    rr.data[1] = (t + TWO) * x1 - (t + TWO) * t.exp();

    0
}

/*-------------------------------------
 * Custom matrix-embedded linear solver
 *-------------------------------------*/

/* In C this is built with SUNLinSolNewEmpty + ops-table overrides
   (MatrixEmbeddedLS / MatrixEmbeddedLSType / MatrixEmbeddedLSSolve /
   MatrixEmbeddedLSFree); here it is a CustomLinSol implementation.
   The integrator hands solve() the current (tcur, cj) and user data,
   replacing the C call to IDAGetNonlinearSystemData (IDA's cj rides in
   the shared CustomLinSol `gamma` slot). */
struct MatrixEmbeddedLS;

impl CustomLinSol for MatrixEmbeddedLS {
    /* linear solve routine */
    fn solve(
        &mut self,
        x: &mut NVector,
        b: &NVector,
        _tol: f64,
        tcur: f64,
        cj: f64,
        user_data: &mut UserData,
    ) -> i32 {
        /* extract stiffness parameter from user_data */
        let alpha = *user_data.as_mut().unwrap().downcast_mut::<f64>().unwrap();
        let ONE = 1.0;
        let TWO = 2.0;

        /* perform linear solve: A*x=b
               A = df/dy + cj*df/dyp
            =>
               A = [ - cj - (alpha - 1)/(t - 2) - 1, alpha - 1]
                   [                          t + 2,         0]
        */
        let a11 = -cj - (alpha - ONE) / (tcur - TWO) - ONE;
        let a12 = alpha - ONE;
        let a21 = tcur + TWO;
        let b1 = b.data[0];
        let b2 = b.data[1];
        x.data[0] = b2 / a21;
        x.data[1] = -(a11 * b2 - a21 * b1) / (a12 * a21);

        /* return with success */
        SUN_SUCCESS
    }
}

/*-------------------------------
 * Private helper functions
 *-------------------------------*/

fn check_retval(retval: i32, funcname: &str) -> bool {
    if retval < 0 {
        eprintln!("\nSUNDIALS_ERROR: {}() failed with flag = {}\n", funcname, retval);
        return true;
    }
    false
}

/* routine to fill analytical solution and its derivative */
fn analytical_solution(t: f64, y: &mut NVector, yp: &mut NVector) {
    y.data[0] = t.exp();
    y.data[1] = t.exp() / (t - 2.0);
    yp.data[0] = t.exp();
    yp.data[1] = t.exp() / (t - 2.0) - t.exp() / (t - 2.0) / (t - 2.0);
}

/* check the computed solution */
fn check_ans(y: &NVector, t: f64, rtol: f64, atol: f64) -> i32 {
    let ONE = 1.0;

    /* create solution and error weight vectors */
    let mut ytrue = N_VClone(y);
    let mut ewt = N_VClone(y);
    let mut abstol = N_VClone(y);

    /* set the solution data */
    analytical_solution(t, &mut ytrue, &mut abstol);

    /* compute the error weight vector, loosen atol */
    N_VConst(atol, &mut abstol);
    N_VAbs(&ytrue, &mut ewt);
    ewt.linear_sum_with(rtol, 10.0, &abstol);
    if N_VMin(&ewt) <= 0.0 {
        eprintln!("\nSUNDIALS_ERROR: check_ans failed - ewt <= 0\n");
        return -1;
    }
    ewt.invert_inplace();

    /* compute the solution error */
    ytrue.linear_sum_with(-ONE, ONE, y);
    let err = N_VWrmsNorm(&ytrue, &ewt);

    /* is the solution within the tolerances? */
    let passfail = if err < ONE { 0 } else { 1 };

    if passfail != 0 {
        print!("\nSUNDIALS_WARNING: check_ans error={}\n\n", fmt_g(err, 0, 6));
    }

    passfail
}

/* Main Program */
fn main() {
    /* general problem parameters */
    let T0: f64 = 0.0; /* initial time */
    let Tf: f64 = 1.0; /* final time */
    let dTout: f64 = 0.1; /* time between outputs */
    let NEQ: i64 = 2; /* number of dependent vars. */
    let reltol: f64 = 1.0e-4; /* tolerances */
    let abstol: f64 = 1.0e-9;
    let alpha: f64 = 10.0; /* stiffness parameter */

    /* Initial diagnostics output */
    println!("\nAnalytical DAE test problem:");
    println!("    alpha = {}", fmt_g(alpha, 0, 6));
    println!("   reltol = {}", fmt_e(reltol, 0, 1));
    println!("   abstol = {}\n", fmt_e(abstol, 0, 1));

    /* Create the SUNDIALS context object for this simulation */
    let sunctx = SUNContext_Create();

    /* Initialize data structures */
    let mut yy = N_VNew_Serial(NEQ, &sunctx); /* Create serial vector for solution */
    let mut yp = N_VClone(&yy); /* Create serial vector for solution derivative */
    analytical_solution(T0, &mut yy, &mut yp); /* Specify initial conditions */

    /* Call IDACreate and IDAInit to initialize IDA memory */
    let mut ida_mem = IDACreate(&sunctx);
    let mut retval = IDAInit(&mut ida_mem, fres, T0, &yy, &yp);
    if check_retval(retval, "IDAInit") {
        std::process::exit(1);
    }

    /* Set routines */
    retval = IDASetUserData(&mut ida_mem, Some(Box::new(alpha)));
    if check_retval(retval, "IDASetUserData") {
        std::process::exit(1);
    }
    retval = IDASStolerances(&mut ida_mem, reltol, abstol);
    if check_retval(retval, "IDASStolerances") {
        std::process::exit(1);
    }

    /* Create custom matrix-embedded linear solver and attach it (NULL matrix) */
    retval = IDASetLinearSolver(
        &mut ida_mem,
        LinearSolver::Custom(Box::new(MatrixEmbeddedLS)),
        None,
    );
    if check_retval(retval, "IDASetLinearSolver") {
        std::process::exit(1);
    }

    /* Override any current settings with command-line options */
    let args: Vec<String> = std::env::args().collect();
    retval = IDASetOptions(&mut ida_mem, "", "", &args);
    if check_retval(retval, "IDASetOptions") {
        std::process::exit(1);
    }

    /* In loop, call IDASolve, print results, and test for error.
       Stops when the final time has been reached. */
    let mut t = T0;
    let mut tout = T0 + dTout;
    println!("        t          x1         x2");
    println!("   ----------------------------------");
    while Tf - t > 1.0e-15 {
        retval = IDASolve(&mut ida_mem, tout, &mut t, &mut yy, &mut yp, IDA_NORMAL);
        if check_retval(retval, "IDASolve") {
            std::process::exit(1);
        }
        println!(
            "  {}  {}  {}",
            fmt_f(t, 10, 6),
            fmt_f(yy.data[0], 10, 6),
            fmt_f(yy.data[1], 10, 6)
        );
        if retval >= 0 {
            /* successful solve: update time */
            tout += dTout;
            tout = if tout > Tf { Tf } else { tout };
        } else {
            /* unsuccessful solve: break */
            eprintln!("Solver failure, stopping integration");
            break;
        }
    }
    println!("   ----------------------------------");

    /* Get/print some final statistics on how the solve progressed */
    let mut h0 = 0.0f64;
    let (mut nst, mut nre, mut nni, mut netf, mut ncfn, mut nreLS) =
        (0i64, 0i64, 0i64, 0i64, 0i64, 0i64);
    retval = IDAGetActualInitStep(&mut ida_mem, &mut h0);
    check_retval(retval, "IDAGetActualInitStep");
    retval = IDAGetNumSteps(&mut ida_mem, &mut nst);
    check_retval(retval, "IDAGetNumSteps");
    retval = IDAGetNumResEvals(&mut ida_mem, &mut nre);
    check_retval(retval, "IDAGetNumResEvals");
    retval = IDAGetNumNonlinSolvIters(&mut ida_mem, &mut nni);
    check_retval(retval, "IDAGetNumNonlinSolvIters");
    retval = IDAGetNumErrTestFails(&mut ida_mem, &mut netf);
    check_retval(retval, "IDAGetNumErrTestFails");
    retval = IDAGetNumNonlinSolvConvFails(&mut ida_mem, &mut ncfn);
    check_retval(retval, "IDAGetNumNonlinSolvConvFails");
    retval = IDAGetNumLinResEvals(&mut ida_mem, &mut nreLS);
    check_retval(retval, "IDAGetNumLinResEvals");

    println!("\nFinal Solver Statistics: \n");
    println!("Initial time step                  = {}", fmt_f(h0, 8, 6));
    println!("Number of steps                    = {}", nst);
    println!("Number of residual evaluations     = {}", nre + nreLS);
    println!("Number of nonlinear iterations     = {}", nni);
    println!("Number of error test failures      = {}", netf);
    println!("Number of nonlinear conv. failures = {}", ncfn);

    /* check the solution error */
    let retval = check_ans(&yy, t, reltol, abstol);

    /* Clean up and return (RAII) */
    IDAFree(ida_mem);

    std::process::exit(retval);
}
