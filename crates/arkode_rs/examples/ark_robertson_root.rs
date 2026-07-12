/* -----------------------------------------------------------------
 * Translation of examples/arkode/C_serial/ark_robertson_root.c
 * (SUNDIALS 7.7.0).
 *
 * Example problem: the following test simulates the Robertson
 * problem, corresponding to the kinetics of an autocatalytic
 * reaction, with added rootfinding.  This is an ODE system with 3
 * components, Y = [u,v,w], satisfying the equations,
 *    du/dt = -0.04*u + 1e4*v*w
 *    dv/dt = 0.04*u - 1e4*v*w - 3e7*v^2
 *    dw/dt = 3e7*v^2
 * for t in the interval [0.0, 1e11], with initial conditions
 * Y0 = [1,0,0].
 *
 * While integrating the system, we use the rootfinding feature
 * to find the times at which either u=1e-4 or w=1e-2.
 *
 * This program solves the problem with the DIRK method, using a
 * Newton iteration with the dense SUNLinearSolver, and a
 * user-supplied Jacobian routine.
 *
 * 100 outputs are printed at equal intervals, and run statistics
 * are printed at the end.
 * -----------------------------------------------------------------*/

use arkode_rs::arkode::{ARKodeEvolve, ARKodeFree, ARKodeSVtolerances};
use arkode_rs::arkode_io::{
    ARKodeGetNumErrTestFails, ARKodeGetNumGEvals, ARKodeGetNumLinSolvSetups,
    ARKodeGetNumNonlinSolvConvFails, ARKodeGetNumNonlinSolvIters, ARKodeGetNumRhsEvals,
    ARKodeGetNumStepAttempts, ARKodeGetNumStepSolveFails, ARKodeGetNumSteps, ARKodeGetRootInfo,
    ARKodeSetMaxErrTestFails, ARKodeSetMaxNonlinIters, ARKodeSetMaxNumSteps,
    ARKodeSetNonlinConvCoef, ARKodeSetPredictorMethod,
};
use arkode_rs::arkode_arkstep::ARKStepCreate;
use arkode_rs::arkode_ls::{
    ARKodeGetNumJacEvals, ARKodeGetNumLinRhsEvals, ARKodeSetJacFn, ARKodeSetLinearSolver,
};
use arkode_rs::arkode_root::ARKodeRootInit;
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

/* g routine to compute the root-finding function g(t,y). */
fn g(_t: f64, y: &NVector, gout: &mut [f64], _user_data: &mut UserData) -> i32 {
    let u = y.data[0]; /* access current solution */
    let w = y.data[2];

    gout[0] = u - 0.0001; /* check for u == 1e-4 */
    gout[1] = w - 0.01; /* check for w == 1e-2 */

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

fn main() {
    /* general problem parameters */
    let t0: f64 = 0.0; /* initial time */
    let t1: f64 = 0.4; /* first output time */
    let tmult: f64 = 10.0; /* output time multiplication factor */
    let nt = 12; /* total number of output times */
    let neq: i64 = 3; /* number of dependent vars. */

    /* set up the initial conditions */
    let u0: f64 = 1.0;
    let v0: f64 = 0.0;
    let w0: f64 = 0.0;
    let reltol: f64 = 1.0e-4;

    /* Create the SUNDIALS context object for this simulation */
    let ctx = SUNContext_Create();

    /* Initial problem output */
    println!("\nRobertson ODE test problem (with rootfinding):");
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

    let mut atols = N_VClone(&y); /* Create serial vector absolute tolerances */

    /* Set absolute tolerances */
    atols.data[0] = 1.0e-8;
    atols.data[1] = 1.0e-11;
    atols.data[2] = 1.0e-8;

    /* Call ARKStepCreate to initialize the ARK timestepper module and
    specify the right-hand side function in y'=f(t,y), the initial time
    T0, and the initial dependent variable vector y.  Note: since this
    problem is fully implicit, we set f_E to NULL and f_I to f. */
    let mut arkode_mem = ARKStepCreate(None, Some(f), t0, &y, &ctx).expect("ARKStepCreate");

    /* Set routines */
    let flag = ARKodeSetMaxErrTestFails(&mut arkode_mem, 20);
    assert!(flag >= 0, "ARKodeSetMaxErrTestFails failed with flag = {}", flag);
    let flag = ARKodeSetMaxNonlinIters(&mut arkode_mem, 8);
    assert!(flag >= 0, "ARKodeSetMaxNonlinIters failed with flag = {}", flag);
    let flag = ARKodeSetNonlinConvCoef(&mut arkode_mem, 1.0e-7);
    assert!(flag >= 0, "ARKodeSetNonlinConvCoef failed with flag = {}", flag);
    let flag = ARKodeSetMaxNumSteps(&mut arkode_mem, 100000);
    assert!(flag >= 0, "ARKodeSetMaxNumSteps failed with flag = {}", flag);
    let flag = ARKodeSetPredictorMethod(&mut arkode_mem, 1);
    assert!(flag >= 0, "ARKodeSetPredictorMethod failed with flag = {}", flag);
    let flag = ARKodeSVtolerances(&mut arkode_mem, reltol, &atols);
    assert!(flag >= 0, "ARKodeSStolerances failed with flag = {}", flag);

    /* Specify the root-finding function, having 2 equations */
    let flag = ARKodeRootInit(&mut arkode_mem, 2, Some(g));
    assert!(flag >= 0, "ARKodeRootInit failed with flag = {}", flag);

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
    println!("        t             u             v             w");
    println!("   -----------------------------------------------------");
    println!(
        "  {}  {}  {}  {}",
        fmt_e(t, 12, 5),
        fmt_e(y.data[0], 12, 5),
        fmt_e(y.data[1], 12, 5),
        fmt_e(y.data[2], 12, 5)
    );
    let mut tout = t1;
    let mut iout = 0;
    let mut rootsfound = [0i32; 2];
    loop {
        let flag = ARKodeEvolve(&mut arkode_mem, tout, &mut y, &mut t, ARK_NORMAL);
        if flag < 0 {
            eprintln!("Solver failure, stopping integration");
            break;
        }
        println!(
            "  {}  {}  {}  {}",
            fmt_e(t, 12, 5),
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
        if flag == ARK_ROOT_RETURN {
            /* check if a root was found */
            let rtflag = ARKodeGetRootInfo(&mut arkode_mem, &mut rootsfound);
            assert!(rtflag >= 0, "ARKodeGetRootInfo failed with flag = {}", rtflag);
            println!("      rootsfound[] = {:3} {:3}", rootsfound[0], rootsfound[1]);
        }
        /* successful solve: update output time */
        iout += 1;
        tout *= tmult;
        if iout == nt {
            break; /* stop after enough outputs */
        }
    }
    println!("   -----------------------------------------------------");
    drop(ufid);

    /* Print some final statistics */
    let (mut nst, mut nst_a, mut nfe, mut nfi) = (0i64, 0i64, 0i64, 0i64);
    let (mut nsetups, mut nje, mut nfe_ls, mut nni) = (0i64, 0i64, 0i64, 0i64);
    let (mut nnf, mut ncfn, mut netf, mut nge) = (0i64, 0i64, 0i64, 0i64);
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
    ARKodeGetNumLinRhsEvals(&mut arkode_mem, &mut nfe_ls);
    ARKodeGetNumGEvals(&mut arkode_mem, &mut nge);

    println!("\nFinal Solver Statistics:");
    println!("   Internal solver steps = {} (attempted = {})", nst, nst_a);
    println!("   Total RHS evals:  Fe = {},  Fi = {}", nfe, nfi);
    println!("   Total linear solver setups = {}", nsetups);
    println!("   Total RHS evals for setting up the linear system = {}", nfe_ls);
    println!("   Total number of Jacobian evaluations = {}", nje);
    println!("   Total number of Newton iterations = {}", nni);
    println!("   Total root-function g evals = {}", nge);
    println!("   Total number of nonlinear solver convergence failures = {}", nnf);
    println!("   Total number of error test failures = {}", netf);
    println!("   Total number of failed steps from solver failure = {}", ncfn);

    /* Clean up and return with successful completion */
    drop(y); /* Free y vector */
    drop(atols); /* Free atols vector */
    let mut slot = Some(arkode_mem);
    ARKodeFree(&mut slot); /* Free integrator memory */
}
