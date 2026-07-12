/* -----------------------------------------------------------------
 * Translation of examples/arkode/C_serial/ark_analytic_mels.c
 * (SUNDIALS 7.7.0).
 *
 * Example problem: the following is a simple example problem with
 * analytical solution,
 *    dy/dt = lambda*y + 1/(1+t^2) - lambda*atan(t)
 * for t in the interval [0.0, 10.0], with initial condition: y=0.
 *
 * The stiffness of the problem is directly proportional to the
 * value of "lambda".  The value of lambda should be negative to
 * result in a well-posed ODE; for values with magnitude larger
 * than 100 the problem becomes quite stiff.
 *
 * This program solves the problem with the DIRK method and a
 * custom 'matrix-embedded' SUNLinearSolver. Output is printed
 * every 1.0 units of time (10 total).  Run statistics (optional
 * outputs) are printed at the end.
 *
 * Translation note: C's solver content is the arkode_mem pointer and
 * its solve callback queries ARKodeGetNonlinearSystemData; the Rust
 * CustomLinSol adaptation instead receives (t, gamma, user_data)
 * directly from the ARKLS solve path (same pinned pattern as
 * cvAnalytic_mels).
 * -----------------------------------------------------------------*/

use arkode_rs::arkode::{ARKodeEvolve, ARKodeFree, ARKodeSStolerances};
use arkode_rs::arkode_arkstep::ARKStepCreate;
use arkode_rs::arkode_io::{
    ARKodeGetNumErrTestFails, ARKodeGetNumLinSolvSetups, ARKodeGetNumNonlinSolvConvFails,
    ARKodeGetNumNonlinSolvIters, ARKodeGetNumRhsEvals, ARKodeGetNumStepAttempts,
    ARKodeGetNumSteps, ARKodeSetLinear, ARKodeSetUserData,
};
use arkode_rs::arkode_ls::{ARKodeGetNumJacEvals, ARKodeGetNumLinRhsEvals, ARKodeSetLinearSolver};
use arkode_rs::sundials_linearsolver::CustomLinSol;
use arkode_rs::sundials_utils::{fmt_e, fmt_f, fmt_g};
use arkode_rs::*;

/* f routine to compute the ODE RHS function f(t,y). */
fn f(t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let rdata = user_data.as_mut().unwrap().downcast_mut::<[f64; 1]>().unwrap();
    let lambda = rdata[0]; /* set shortcut for stiffness parameter */
    let u = y.data[0]; /* access current solution value */

    /* fill in the RHS function */
    ydot.data[0] = lambda * u + 1.0 / (1.0 + t * t) - lambda * t.atan();

    0 /* return with success */
}

/*-------------------------------------
 * Custom matrix-embedded linear solver
 *-------------------------------------*/

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
        let rdata = user_data.as_mut().unwrap().downcast_mut::<[f64; 1]>().unwrap();
        let lambda = rdata[0];

        /* perform linear solve: (1-gamma*lambda)*x = b */
        x.data[0] = b.data[0] / (1.0 - gamma * lambda);

        /* return with success */
        SUN_SUCCESS
    }
}

/* check the computed solution */
fn check_ans(y: &NVector, t: f64, rtol: f64, atol: f64) -> i32 {
    /* compute solution error */
    let ans = t.atan();
    let ewt = 1.0 / (rtol * SUNRabs(ans) + atol);
    let err = ewt * SUNRabs(y.data[0] - ans);

    /* The local errors accumulate from step to step so that the global error
     * is not quite within the local error tolerances. This factor accounts
     * for this. */
    let global_bound: f64 = 1.5;

    /* is the solution within the tolerances? */
    let passfail = if err < global_bound { 0 } else { 1 };

    if passfail != 0 {
        println!("\nSUNDIALS_WARNING: check_ans error={}\n", fmt_g(err, 0, 6));
    }

    passfail
}

/* Main Program */
fn main() {
    /* general problem parameters */
    let t0: f64 = 0.0; /* initial time */
    let tf: f64 = 10.0; /* final time */
    let dtout: f64 = 1.0; /* time between outputs */
    let neq: i64 = 1; /* number of dependent vars. */
    let reltol: f64 = 1.0e-5; /* tolerances */
    let abstol: f64 = 1.0e-10;
    let lambda: f64 = -100.0; /* stiffness parameter */

    /* Create the SUNDIALS context object for this simulation */
    let ctx = SUNContext_Create();

    /* Initial diagnostics output */
    println!("\nAnalytical ODE test problem:");
    println!("   lambda = {}", fmt_g(lambda, 0, 6));
    println!("   reltol = {}", fmt_e(reltol, 0, 1));
    println!("   abstol = {}\n", fmt_e(abstol, 0, 1));

    /* Initialize data structures */
    let mut y = N_VNew_Serial(neq, &ctx); /* Create serial vector for solution */
    N_VConst(0.0, &mut y); /* Specify initial condition */

    /* Call ARKStepCreate to initialize the ARK timestepper module and
    specify the right-hand side function in y'=f(t,y), the initial time
    T0, and the initial dependent variable vector y.  Note: since this
    problem is fully implicit, we set f_E to NULL and f_I to f. */
    let mut arkode_mem = ARKStepCreate(None, Some(f), t0, &y, &ctx).expect("ARKStepCreate");

    /* Set routines */
    let retval = ARKodeSetUserData(&mut arkode_mem, Some(Box::new([lambda])));
    assert!(retval >= 0, "ARKodeSetUserData failed with retval = {}", retval);
    let retval = ARKodeSStolerances(&mut arkode_mem, reltol, abstol);
    assert!(retval >= 0, "ARKodeSStolerances failed with retval = {}", retval);

    /* Initialize custom matrix-embedded linear solver */
    let ls = LinearSolver::Custom(Box::new(MatrixEmbeddedLS));
    let retval = ARKodeSetLinearSolver(&mut arkode_mem, ls, None); /* Attach linear solver */
    assert!(retval >= 0, "ARKodeSetLinearSolver failed with retval = {}", retval);

    /* Specify linearly implicit RHS, with non-time-dependent Jacobian */
    let retval = ARKodeSetLinear(&mut arkode_mem, 0);
    assert!(retval >= 0, "ARKodeSetLinear failed with retval = {}", retval);

    /* Main time-stepping loop: calls ARKodeEvolve to perform the
    integration, then prints results.  Stops when the final time has been
    reached. */
    let mut t = t0;
    let mut tout = t0 + dtout;
    println!("        t           u");
    println!("   ---------------------");
    while tf - t > 1.0e-15 {
        let retval = ARKodeEvolve(&mut arkode_mem, tout, &mut y, &mut t, ARK_NORMAL);
        if retval < 0 {
            eprintln!("\nSUNDIALS_ERROR: ARKodeEvolve() failed with retval = {}\n", retval);
            break;
        }
        println!("  {}  {}", fmt_f(t, 10, 6), fmt_f(y.data[0], 10, 6));
        if retval >= 0 {
            /* successful solve: update time */
            tout += dtout;
            tout = if tout > tf { tf } else { tout };
        } else {
            /* unsuccessful solve: break */
            eprintln!("Solver failure, stopping integration");
            break;
        }
    }
    println!("   ---------------------");

    /* Get/print some final statistics on how the solve progressed */
    let mut nst: i64 = 0;
    let mut nst_a: i64 = 0;
    let mut nfe: i64 = 0;
    let mut nfi: i64 = 0;
    let mut nsetups: i64 = 0;
    let mut nje: i64 = 0;
    let mut nfe_ls: i64 = 0;
    let mut nni: i64 = 0;
    let mut ncfn: i64 = 0;
    let mut netf: i64 = 0;
    ARKodeGetNumSteps(&mut arkode_mem, &mut nst);
    ARKodeGetNumStepAttempts(&mut arkode_mem, &mut nst_a);
    ARKodeGetNumRhsEvals(&mut arkode_mem, 0, &mut nfe);
    ARKodeGetNumRhsEvals(&mut arkode_mem, 1, &mut nfi);
    ARKodeGetNumLinSolvSetups(&mut arkode_mem, &mut nsetups);
    ARKodeGetNumErrTestFails(&mut arkode_mem, &mut netf);
    ARKodeGetNumNonlinSolvIters(&mut arkode_mem, &mut nni);
    ARKodeGetNumNonlinSolvConvFails(&mut arkode_mem, &mut ncfn);
    ARKodeGetNumJacEvals(&mut arkode_mem, &mut nje);
    ARKodeGetNumLinRhsEvals(&mut arkode_mem, &mut nfe_ls);

    println!("\nFinal Solver Statistics:");
    println!("   Internal solver steps = {} (attempted = {})", nst, nst_a);
    println!("   Total RHS evals:  Fe = {},  Fi = {}", nfe, nfi);
    println!("   Total linear solver setups = {}", nsetups);
    println!("   Total RHS evals for setting up the linear system = {}", nfe_ls);
    println!("   Total number of Jacobian evaluations = {}", nje);
    println!("   Total number of Newton iterations = {}", nni);
    println!("   Total number of linear solver convergence failures = {}", ncfn);
    println!("   Total number of error test failures = {}\n", netf);

    /* check the solution error */
    let retval = check_ans(&y, t, reltol, abstol);

    /* Clean up and return */
    drop(y); /* Free y vector */
    let mut slot = Some(arkode_mem);
    ARKodeFree(&mut slot); /* Free integrator memory */

    std::process::exit(retval);
}
