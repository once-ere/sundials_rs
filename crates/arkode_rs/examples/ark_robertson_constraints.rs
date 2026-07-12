/* -----------------------------------------------------------------
 * Translation of examples/arkode/C_serial/ark_robertson_constraints.c
 * (SUNDIALS 7.7.0).
 *
 * Example problem: the following test simulates the Robertson
 * problem, corresponding to the kinetics of an autocatalytic
 * reaction.  This is an ODE system with 3 components, Y = [u,v,w],
 * satisfying the equations,
 *    du/dt = -0.04*u + 1e4*v*w
 *    dv/dt = 0.04*u - 1e4*v*w - 3e7*v^2
 *    dw/dt = 3e7*v^2
 * for t in the interval [0.0, 1e11], with initial conditions
 * Y0 = [1,0,0].
 *
 * This program solves the problem with one of the solvers, ERK,
 * DIRK or ARK.  For DIRK and ARK, implicit subsystems are solved
 * using a Newton iteration with the dense SUNLinearSolver, and a
 * user-supplied Jacobian routine. The constraint y_i >= 0 is
 * posed for all components.
 *
 * 100 outputs are printed at equal intervals, and run statistics
 * are printed at the end.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]

use arkode_rs::arkode::{ARKodeEvolve, ARKodeFree, ARKodeSStolerances};
use arkode_rs::arkode_arkstep::ARKStepCreate;
use arkode_rs::arkode_io::{
    ARKodeGetNumConstrFails, ARKodeGetNumErrTestFails, ARKodeGetNumLinSolvSetups,
    ARKodeGetNumNonlinSolvConvFails, ARKodeGetNumNonlinSolvIters, ARKodeGetNumRhsEvals,
    ARKodeGetNumStepAttempts, ARKodeGetNumStepSolveFails, ARKodeGetNumSteps, ARKodeSetConstraints,
    ARKodeSetInitStep, ARKodeSetMaxErrTestFails, ARKodeSetMaxNonlinIters, ARKodeSetMaxNumSteps,
    ARKodeSetNonlinConvCoef, ARKodeSetPredictorMethod,
};
use arkode_rs::arkode_ls::{
    ARKodeGetNumJacEvals, ARKodeGetNumLinRhsEvals, ARKodeSetJacFn, ARKodeSetLinearSolver,
};
use arkode_rs::sundials_utils::{fmt_e, fmt_g};
use arkode_rs::*;
use std::io::Write;

/* f routine to compute the ODE RHS function f(t,y). */
fn f(_t: f64, y: &NVector, ydot: &mut NVector, _user_data: &mut UserData) -> i32 {
    let u = y.data[0]; /* access current solution */
    let v = y.data[1];
    let w = y.data[2];

    /* Fill in ODE RHS function */
    ydot.data[0] = -0.04 * u + 1.0e4 * v * w;
    ydot.data[1] = 0.04 * u - 1.0e4 * v * w - 3.0e7 * v * v;
    ydot.data[2] = 3.0e7 * v * v;

    0 /* Return with success */
}

/* Jacobian routine to compute J(t,y) = df/dy. */
#[allow(clippy::too_many_arguments)] /* fixed ARKLsJacFn signature */
fn jac(
    _t: f64,
    y: &NVector,
    _fy: &NVector,
    j: &mut SUNMatrix,
    _user_data: &mut UserData,
    _tmp1: &mut NVector,
    _tmp2: &mut NVector,
    _tmp3: &mut NVector,
) -> i32 {
    let v = y.data[1]; /* access current solution */
    let w = y.data[2];
    SUNMatZero(j); /* initialize Jacobian to zero */

    /* Fill in the Jacobian of the ODE RHS function (column-major) */
    if let SUNMatrix::Dense(dm) = j {
        let m = 3usize;
        dm.data[0] = -0.04; /* (0,0) */
        dm.data[m] = 1.0e4 * w; /* (0,1) */
        dm.data[2 * m] = 1.0e4 * v; /* (0,2) */

        dm.data[1] = 0.04; /* (1,0) */
        dm.data[m + 1] = -1.0e4 * w - 6.0e7 * v; /* (1,1) */
        dm.data[2 * m + 1] = -1.0e4 * v; /* (1,2) */

        dm.data[m + 2] = 6.0e7 * v; /* (2,1) */
    }

    0 /* Return with success */
}

/* compare the solution at the final time 1e11s to a reference solution
   computed using a relative tolerance of 1e-8 and absolute tolerance of
   1e-14 */
#[allow(clippy::excessive_precision)] /* C source's full decimal text kept */
fn check_ans(y: &NVector, _t: f64, rtol: f64, atol: f64) -> i32 {
    /* create reference solution and error weight vectors */
    let mut refv = N_VClone(y);
    let mut ewt = N_VClone(y);

    /* set the reference solution data */
    refv.data[0] = 2.0833403356917897e-08;
    refv.data[1] = 8.1470714598028223e-14;
    refv.data[2] = 9.9999997916651040e-01;

    /* compute the error weight vector */
    N_VAbs(&refv.clone(), &mut ewt);
    ewt.scale_inplace(rtol);
    ewt.add_const_inplace(atol);
    if N_VMin(&ewt) <= 0.0 {
        eprintln!("\nSUNDIALS_ERROR: check_ans failed - ewt <= 0\n");
        return -1;
    }
    ewt.invert_inplace();

    /* compute the solution error */
    refv.linear_sum_with(-1.0, 1.0, y);
    let err = N_VWrmsNorm(&refv, &ewt);

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
    let t0: f64 = 0.0; /* initial time */
    let tf: f64 = 1.0e11; /* final time */
    let dtout: f64 = (tf - t0) / 100.0; /* time between outputs */
    let nt = (tf / dtout).ceil() as i32; /* number of output times */
    let neq: i64 = 3; /* number of dependent vars. */

    /* set up the initial conditions, tolerances, initial time step size */
    let u0: f64 = 1.0;
    let v0: f64 = 0.0;
    let w0: f64 = 0.0;
    let reltol: f64 = 1.0e-3;
    let abstol: f64 = 1.0e-7;
    let h0: f64 = 1.0e-4 * reltol;

    /* Create the SUNDIALS context object for this simulation */
    let ctx = SUNContext_Create();

    /* Initial problem output */
    println!("\nRobertson ODE test problem:");
    println!(
        "    initial conditions:  u0 = {},  v0 = {},  w0 = {}",
        fmt_g(u0, 0, 6),
        fmt_g(v0, 0, 6),
        fmt_g(w0, 0, 6)
    );

    /* Initialize data structures */
    let mut y = N_VNew_Serial(neq, &ctx); /* Create serial vector for solution */
    y.data[0] = u0; /* Set initial conditions into y */
    y.data[1] = v0;
    y.data[2] = w0;

    let mut constraints = N_VClone(&y);
    /* Set constraints to all 1's for nonnegative solution values. */
    N_VConst(1.0, &mut constraints);

    /* Call ARKStepCreate to initialize the ARK timestepper module and
    specify the right-hand side function in y'=f(t,y), the initial time
    T0, and the initial dependent variable vector y.  Note: since this
    problem is fully implicit, we set f_E to NULL and f_I to f. */
    let mut arkode_mem = ARKStepCreate(None, Some(f), t0, &y, &ctx).expect("ARKStepCreate");

    /* Set routines */
    let flag = ARKodeSetInitStep(&mut arkode_mem, h0); /* Set custom initial step */
    assert!(flag >= 0, "ARKodeSetInitStep failed with flag = {}", flag);
    let flag = ARKodeSetMaxErrTestFails(&mut arkode_mem, 20); /* Increase max error test fails */
    assert!(flag >= 0, "ARKodeSetMaxErrTestFails failed with flag = {}", flag);
    let flag = ARKodeSetMaxNonlinIters(&mut arkode_mem, 8); /* Increase max nonlin iters  */
    assert!(flag >= 0, "ARKodeSetMaxNonlinIters failed with flag = {}", flag);
    let flag = ARKodeSetNonlinConvCoef(&mut arkode_mem, 1.0e-7); /* Set nonlinear convergence coeff. */
    assert!(flag >= 0, "ARKodeSetNonlinConvCoef failed with flag = {}", flag);
    let flag = ARKodeSetMaxNumSteps(&mut arkode_mem, 100000); /* Increase max num steps */
    assert!(flag >= 0, "ARKodeSetMaxNumSteps failed with flag = {}", flag);
    let flag = ARKodeSetPredictorMethod(&mut arkode_mem, 1); /* Specify maximum-order predictor */
    assert!(flag >= 0, "ARKodeSetPredictorMethod failed with flag = {}", flag);
    let flag = ARKodeSStolerances(&mut arkode_mem, reltol, abstol); /* Specify tolerances */
    assert!(flag >= 0, "ARKodeSStolerances failed with flag = {}", flag);
    let flag = ARKodeSetConstraints(&mut arkode_mem, Some(&constraints)); /* Set constraints */
    assert!(flag >= 0, "ARKodeSetConstraints failed with flag = {}", flag);

    /* Initialize dense matrix data structure and solver */
    let a = SUNDenseMatrix(neq, neq, &ctx);
    let ls = SUNLinSol_Dense(&y, &a, &ctx);

    /* Linear solver interface */
    let flag = ARKodeSetLinearSolver(&mut arkode_mem, ls, Some(a));
    assert!(flag >= 0, "ARKodeSetLinearSolver failed with flag = {}", flag);
    let flag = ARKodeSetJacFn(&mut arkode_mem, Some(jac));
    assert!(flag >= 0, "ARKodeSetJacFn failed with flag = {}", flag);

    /* Open output stream for results, output comment line */
    let mut ufid = std::fs::File::create("solution.txt").unwrap();
    let _ = writeln!(ufid, "# t u v w");

    /* output initial condition to disk */
    let _ = writeln!(
        ufid,
        " {} {} {} {}",
        fmt_e(t0, 0, 16),
        fmt_e(y.data[0], 0, 16),
        fmt_e(y.data[1], 0, 16),
        fmt_e(y.data[2], 0, 16)
    );

    /* Main time-stepping loop: calls ARKodeEvolve to perform the
    integration, then prints results.  Stops when the final time has
    been reached */
    let mut t = t0;
    let mut tout = t0 + dtout;
    println!("        t           u           v           w");
    println!("   --------------------------------------------------");
    println!(
        "  {}  {}  {}  {}",
        fmt_e(t, 10, 3),
        fmt_e(y.data[0], 12, 5),
        fmt_e(y.data[1], 12, 5),
        fmt_e(y.data[2], 12, 5)
    );
    for _iout in 0..nt {
        let flag = ARKodeEvolve(&mut arkode_mem, tout, &mut y, &mut t, ARK_NORMAL);
        if flag < 0 {
            break;
        }
        println!(
            "  {}  {}  {}  {}",
            fmt_e(t, 10, 3),
            fmt_e(y.data[0], 12, 5),
            fmt_e(y.data[1], 12, 5),
            fmt_e(y.data[2], 12, 5)
        );
        let _ = writeln!(
            ufid,
            " {} {} {} {}",
            fmt_e(t, 0, 16),
            fmt_e(y.data[0], 0, 16),
            fmt_e(y.data[1], 0, 16),
            fmt_e(y.data[2], 0, 16)
        );
        if flag >= 0 {
            /* successful solve: update time */
            tout += dtout;
            tout = if tout > tf { tf } else { tout };
        } else {
            /* unsuccessful solve: break */
            eprintln!("Solver failure, stopping integration");
            break;
        }
    }
    println!("   --------------------------------------------------");
    drop(ufid);

    /* Print some final statistics */
    let mut nst: i64 = 0;
    let mut nst_a: i64 = 0;
    let mut nfe: i64 = 0;
    let mut nfi: i64 = 0;
    let mut nsetups: i64 = 0;
    let mut nje: i64 = 0;
    let mut nfeLS: i64 = 0;
    let mut nni: i64 = 0;
    let mut nnf: i64 = 0;
    let mut ncfn: i64 = 0;
    let mut netf: i64 = 0;
    let mut nctf: i64 = 0;
    ARKodeGetNumSteps(&mut arkode_mem, &mut nst);
    ARKodeGetNumStepAttempts(&mut arkode_mem, &mut nst_a);
    ARKodeGetNumRhsEvals(&mut arkode_mem, 0, &mut nfe);
    ARKodeGetNumRhsEvals(&mut arkode_mem, 1, &mut nfi);
    ARKodeGetNumLinSolvSetups(&mut arkode_mem, &mut nsetups);
    ARKodeGetNumErrTestFails(&mut arkode_mem, &mut netf);
    ARKodeGetNumStepSolveFails(&mut arkode_mem, &mut ncfn);
    ARKodeGetNumNonlinSolvIters(&mut arkode_mem, &mut nni);
    ARKodeGetNumNonlinSolvConvFails(&mut arkode_mem, &mut nnf);
    ARKodeGetNumJacEvals(&mut arkode_mem, &mut nje);
    ARKodeGetNumLinRhsEvals(&mut arkode_mem, &mut nfeLS);
    ARKodeGetNumConstrFails(&mut arkode_mem, &mut nctf);

    println!("\nFinal Solver Statistics:");
    println!("   Internal solver steps = {} (attempted = {})", nst, nst_a);
    println!("   Total RHS evals:  Fe = {},  Fi = {}", nfe, nfi);
    println!("   Total linear solver setups = {}", nsetups);
    println!("   Total RHS evals for setting up the linear system = {}", nfeLS);
    println!("   Total number of Jacobian evaluations = {}", nje);
    println!("   Total number of Newton iterations = {}", nni);
    println!("   Total number of nonlinear solver convergence failures = {}", nnf);
    println!("   Total number of error test failures = {}", netf);
    println!("   Total number of constraint test failures = {}", nctf);
    println!("   Total number of failed steps from solver failure = {}", ncfn);

    /* check the solution error */
    let flag = check_ans(&y, t, reltol, abstol);

    /* Clean up and return with successful completion */
    drop(y); /* Free y vector */
    drop(constraints); /* Free constraints vector */
    let mut slot = Some(arkode_mem);
    ARKodeFree(&mut slot); /* Free integrator memory */

    std::process::exit(flag);
}
