/* -----------------------------------------------------------------
 * Translated from examples/cvode/serial/cvsAnalytic_mels.c (CVODE 7.7.0)
 *
 * Example problem:
 *
 * The following is a simple example problem with analytical
 * solution,
 *    dy/dt = lambda*y + 1/(1+t^2) - lambda*atan(t)
 * for t in the interval [0.0, 10.0], with initial condition: y=0.
 *
 * The stiffness of the problem is directly proportional to the
 * value of "lambda".  The value of lambda should be negative to
 * result in a well-posed ODE; for values with magnitude larger
 * than 100 the problem becomes quite stiff.
 *
 * This program solves the problem with the BDF method, Newton
 * iteration, and a custom 'matrix-embedded' SUNLinearSolver. Output
 * is printed every 1.0 units of time (10 total).  Run statistics
 * (optional outputs) are printed at the end.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use cvodes_rs::cvodes_cli::CVodeSetOptions;
use cvodes_rs::sundials_utils::{fmt_e, fmt_f, fmt_g};
use cvodes_rs::*;

/* f routine to compute the ODE RHS function f(t,y). */
fn f(t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    /* set shortcut for stiffness parameter */
    let lamda = *user_data.as_mut().unwrap().downcast_mut::<f64>().unwrap();
    let u = y.data[0]; /* access current solution value */

    /* fill in the RHS function */
    ydot.data[0] = lamda * u + 1.0 / (1.0 + t * t) - lamda * t.atan();

    0 /* return with success */
}

/*-------------------------------------
 * Custom matrix-embedded linear solver
 *-------------------------------------*/

/* In C this is built with SUNLinSolNewEmpty + ops-table overrides
   (MatrixEmbeddedLS / MatrixEmbeddedLSType / MatrixEmbeddedLSSolve /
   MatrixEmbeddedLSFree); here it is a CustomLinSol implementation.
   The integrator hands solve() the current (t, gamma) and user data,
   replacing the C call to CVodeGetNonlinearSystemData. */
struct MatrixEmbeddedLS;

impl CustomLinSol for MatrixEmbeddedLS {
    /* linear solve routine */
    fn solve(
        &mut self,
        x: &mut NVector,
        b: &NVector,
        _tol: f64,
        _t: f64,
        gamma: f64,
        user_data: &mut UserData,
    ) -> i32 {
        /* extract stiffness parameter from user_data */
        let lamda = *user_data.as_mut().unwrap().downcast_mut::<f64>().unwrap();

        /* perform linear solve: (1-gamma*lamda)*x = b */
        x.data[0] = b.data[0] / (1.0 - gamma * lamda);

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

/* check the computed solution */
fn check_ans(y: &NVector, t: f64, rtol: f64, atol: f64) -> i32 {
    /* compute solution error */
    let ans = t.atan();
    let ewt = 1.0 / (rtol * ans.abs() + atol);
    let err = ewt * (y.data[0] - ans).abs();

    /* is the solution within the tolerances? */
    let passfail = if err < 1.0 { 0 } else { 1 };

    if passfail != 0 {
        println!("\nSUNDIALS_WARNING: check_ans error={}\n", fmt_g(err, 0, 6));
    }

    passfail
}

/* Main Program */
fn main() {
    /* general problem parameters */
    let T0: f64 = 0.0; /* initial time */
    let Tf: f64 = 10.0; /* final time */
    let dTout: f64 = 1.0; /* time between outputs */
    let NEQ: i64 = 1; /* number of dependent vars. */
    let reltol: f64 = 1.0e-6; /* tolerances */
    let abstol: f64 = 1.0e-10;
    let lamda: f64 = -100.0; /* stiffness parameter */

    /* Create the SUNDIALS context */
    let sunctx = SUNContext_Create();

    /* Initial diagnostics output */
    println!("\nAnalytical ODE test problem:");
    println!("   lambda = {}", fmt_g(lamda, 0, 6));
    println!("   reltol = {}", fmt_e(reltol, 0, 1));
    println!("   abstol = {}\n", fmt_e(abstol, 0, 1));

    /* Initialize data structures */
    let mut y = N_VNew_Serial(NEQ, &sunctx); /* Create serial vector for solution */
    N_VConst(0.0, &mut y); /* Specify initial condition */

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

    /* Call CVodeSetUserData to specify the stiffness factor */
    retval = CVodeSetUserData(&mut cvode_mem, Some(Box::new(lamda)));
    if check_retval(retval, "CVodeSetUserData") {
        std::process::exit(1);
    }

    /* Call CVodeSStolerances to specify the scalar relative and absolute
     * tolerances */
    retval = CVodeSStolerances(&mut cvode_mem, reltol, abstol);
    if check_retval(retval, "CVodeSStolerances") {
        std::process::exit(1);
    }

    /* Create custom matrix-embedded linear solver and attach it to CVode */
    retval = CVodeSetLinearSolver(
        &mut cvode_mem,
        LinearSolver::Custom(Box::new(MatrixEmbeddedLS)),
        None,
    );
    if check_retval(retval, "CVodeSetLinearSolver") {
        std::process::exit(1);
    }

    /* Override any current settings with command-line options */
    let args: Vec<String> = std::env::args().collect();
    retval = CVodeSetOptions(&mut cvode_mem, "", "", &args);
    if check_retval(retval, "CVodeSetOptions") {
        std::process::exit(1);
    }

    /* In loop, call CVode, print results, and test for error.
    Break out of loop when NOUT preset output times have been reached.  */
    let mut t = T0;
    let mut tout = T0 + dTout;
    println!("        t           u");
    println!("   ---------------------");
    while Tf - t > 1.0e-15 {
        retval = CVode(&mut cvode_mem, tout, &mut y, &mut t, CV_NORMAL); /* call integrator */
        if check_retval(retval, "CVode") {
            break;
        }
        /* access/print solution */
        println!("  {}  {}", fmt_f(t, 10, 6), fmt_f(y.data[0], 10, 6));
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
    println!("   ---------------------");

    /* Get/print some final statistics on how the solve progressed */
    let (mut nst, mut nfe, mut nsetups, mut nje, mut nfeLS) = (0i64, 0i64, 0i64, 0i64, 0i64);
    let (mut nni, mut ncfn, mut netf) = (0i64, 0i64, 0i64);
    retval = CVodeGetNumSteps(&mut cvode_mem, &mut nst);
    check_retval(retval, "CVodeGetNumSteps");
    retval = CVodeGetNumRhsEvals(&mut cvode_mem, &mut nfe);
    check_retval(retval, "CVodeGetNumRhsEvals");
    retval = CVodeGetNumLinSolvSetups(&mut cvode_mem, &mut nsetups);
    check_retval(retval, "CVodeGetNumLinSolvSetups");
    retval = CVodeGetNumErrTestFails(&mut cvode_mem, &mut netf);
    check_retval(retval, "CVodeGetNumErrTestFails");
    retval = CVodeGetNumNonlinSolvIters(&mut cvode_mem, &mut nni);
    check_retval(retval, "CVodeGetNumNonlinSolvIters");
    retval = CVodeGetNumNonlinSolvConvFails(&mut cvode_mem, &mut ncfn);
    check_retval(retval, "CVodeGetNumNonlinSolvConvFails");
    retval = CVodeGetNumJacEvals(&mut cvode_mem, &mut nje);
    check_retval(retval, "CVodeGetNumJacEvals");
    retval = CVodeGetNumLinRhsEvals(&mut cvode_mem, &mut nfeLS);
    check_retval(retval, "CVodeGetNumLinRhsEvals");

    println!("\nFinal Solver Statistics:");
    println!("   Internal solver steps = {}", nst);
    println!("   Total RHS evals = {}", nfe);
    println!("   Total linear solver setups = {}", nsetups);
    println!("   Total RHS evals for setting up the linear system = {}", nfeLS);
    println!("   Total number of Jacobian evaluations = {}", nje);
    println!("   Total number of Newton iterations = {}", nni);
    println!("   Total number of linear solver convergence failures = {}", ncfn);
    println!("   Total number of error test failures = {}\n", netf);

    /* check the solution error */
    let retval = check_ans(&y, t, reltol, abstol);

    /* Clean up and return (RAII) */
    CVodeFree(cvode_mem);

    std::process::exit(retval);
}
