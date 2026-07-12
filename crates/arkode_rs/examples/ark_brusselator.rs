/* -----------------------------------------------------------------
 * Translation of examples/arkode/C_serial/ark_brusselator.c
 * (SUNDIALS 7.7.0).
 *
 * Example problem: the following test simulates a brusselator
 * problem from chemical kinetics.  This is an ODE system with 3
 * components, Y = [u,v,w], satisfying the equations,
 *    du/dt = a - (w+1)*u + v*u^2
 *    dv/dt = w*u - v*u^2
 *    dw/dt = (b-w)/ep - w*u
 * for t in the interval [0.0, 10.0], with initial conditions
 * Y0 = [u0,v0,w0].  We have 3 different testing scenarios
 * (test = 2 selected below):
 *
 * Test 1:  u0=3.9,  v0=1.1,  w0=2.8,  a=1.2,  b=2.5,  ep=1.0e-5
 * Test 2:  u0=1.2,  v0=3.1,  w0=3,    a=1,    b=3.5,  ep=5.0e-6
 * Test 3:  u0=3,    v0=3,    w0=3.5,  a=0.5,  b=3,    ep=5.0e-4
 *
 * This program solves the problem with the DIRK method, using a
 * Newton iteration with the dense SUNLinearSolver, and a
 * user-supplied Jacobian routine.
 *
 * 100 outputs are printed at equal intervals, and run statistics
 * are printed at the end.
 * -----------------------------------------------------------------*/

use arkode_rs::arkode::{ARKodeEvolve, ARKodeFree, ARKodeSStolerances};
use arkode_rs::arkode_arkstep::ARKStepCreate;
use arkode_rs::arkode_io::{
    ARKodeGetNumErrTestFails, ARKodeGetNumLinSolvSetups, ARKodeGetNumNonlinSolvConvFails,
    ARKodeGetNumNonlinSolvIters, ARKodeGetNumRhsEvals, ARKodeGetNumStepAttempts,
    ARKodeGetNumStepSolveFails, ARKodeGetNumSteps, ARKodeSetAutonomous,
    ARKodeSetDeduceImplicitRhs, ARKodeSetInterpolantType, ARKodeSetUserData,
};
use arkode_rs::arkode_ls::{
    ARKodeGetNumJacEvals, ARKodeGetNumLinRhsEvals, ARKodeSetJacFn, ARKodeSetLinearSolver,
};
use arkode_rs::sundials_utils::{fmt_e, fmt_f, fmt_g};
use arkode_rs::*;
use std::io::Write;

/* f routine to compute the ODE RHS function f(t,y). */
fn f(_t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let rdata = user_data.as_mut().unwrap().downcast_mut::<[f64; 3]>().unwrap();
    let a = rdata[0]; /* access data entries */
    let b = rdata[1];
    let ep = rdata[2];
    let u = y.data[0]; /* access solution values */
    let v = y.data[1];
    let w = y.data[2];

    /* fill in the RHS function */
    ydot.data[0] = a - (w + 1.0) * u + v * u * u;
    ydot.data[1] = w * u - v * u * u;
    ydot.data[2] = (b - w) / ep - w * u;

    0 /* Return with success */
}

/* Jacobian routine to compute J(t,y) = df/dy. */
#[allow(clippy::too_many_arguments)] /* fixed ARKLsJacFn signature */
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
    let rdata = user_data.as_mut().unwrap().downcast_mut::<[f64; 3]>().unwrap();
    let ep = rdata[2]; /* access data entries */
    let u = y.data[0]; /* access solution values */
    let v = y.data[1];
    let w = y.data[2];

    /* fill in the Jacobian (SM_ELEMENT_D(J, i, jcol): column-major) */
    if let SUNMatrix::Dense(dm) = j {
        let m = 3usize;
        dm.data[m] = u * u; /* (0,1) */
        dm.data[2 * m] = -u; /* (0,2) */
        dm.data[0] = -(w + 1.0) + 2.0 * u * v; /* (0,0) */

        dm.data[1] = w - 2.0 * u * v; /* (1,0) */
        dm.data[m + 1] = -u * u; /* (1,1) */
        dm.data[2 * m + 1] = u; /* (1,2) */

        dm.data[2] = -w; /* (2,0) */
        dm.data[m + 2] = 0.0; /* (2,1) */
        dm.data[2 * m + 2] = -1.0 / ep - u; /* (2,2) */
    }

    0 /* Return with success */
}

fn main() {
    /* general problem parameters */
    let t0: f64 = 0.0; /* initial time */
    let tf: f64 = 10.0; /* final time */
    let dtout: f64 = 1.0; /* time between outputs */
    let neq: i64 = 3; /* number of dependent vars. */
    let nt = (tf / dtout).ceil() as i32; /* number of output times */
    let test = 2; /* test problem to run */
    let reltol: f64 = 1.0e-6; /* tolerances */
    let abstol: f64 = 1.0e-10;

    /* Create the SUNDIALS context object for this simulation */
    let ctx = SUNContext_Create();

    /* set up the test problem according to the desired test */
    let (u0, v0, w0, a, b, ep): (f64, f64, f64, f64, f64, f64) = if test == 1 {
        (3.9, 1.1, 2.8, 1.2, 2.5, 1.0e-5)
    } else if test == 3 {
        (3.0, 3.0, 3.5, 0.5, 3.0, 5.0e-4)
    } else {
        (1.2, 3.1, 3.0, 1.0, 3.5, 5.0e-6)
    };

    /* Initial problem output */
    println!("\nBrusselator ODE test problem:");
    println!(
        "    initial conditions:  u0 = {},  v0 = {},  w0 = {}",
        fmt_g(u0, 0, 6),
        fmt_g(v0, 0, 6),
        fmt_g(w0, 0, 6)
    );
    println!(
        "    problem parameters:  a = {},  b = {},  ep = {}",
        fmt_g(a, 0, 6),
        fmt_g(b, 0, 6),
        fmt_g(ep, 0, 6)
    );
    println!(
        "    reltol = {},  abstol = {}\n",
        fmt_e(reltol, 0, 1),
        fmt_e(abstol, 0, 1)
    );

    /* Initialize data structures */
    let rdata: [f64; 3] = [a, b, ep]; /* set user data */
    let mut y = N_VNew_Serial(neq, &ctx); /* Create serial vector for solution */
    y.data[0] = u0; /* Set initial conditions into y */
    y.data[1] = v0;
    y.data[2] = w0;

    /* Call ARKStepCreate to initialize the ARK timestepper module and
    specify the right-hand side function in y'=f(t,y), the initial time
    T0, and the initial dependent variable vector y.  Note: since this
    problem is fully implicit, we set f_E to NULL and f_I to f. */
    let mut arkode_mem = ARKStepCreate(None, Some(f), t0, &y, &ctx).expect("ARKStepCreate");

    /* Set routines */
    let flag = ARKodeSetUserData(&mut arkode_mem, Some(Box::new(rdata)));
    assert!(flag >= 0, "ARKodeSetUserData failed with flag = {}", flag);
    let flag = ARKodeSStolerances(&mut arkode_mem, reltol, abstol);
    assert!(flag >= 0, "ARKodeSStolerances failed with flag = {}", flag);
    let flag = ARKodeSetInterpolantType(&mut arkode_mem, ARK_INTERP_LAGRANGE);
    assert!(flag >= 0, "ARKodeSetInterpolantType failed with flag = {}", flag);
    let flag = ARKodeSetDeduceImplicitRhs(&mut arkode_mem, true);
    assert!(flag >= 0, "ARKodeSetDeduceImplicitRhs failed with flag = {}", flag);

    /* Initialize dense matrix data structure and solver */
    let a_mat = SUNDenseMatrix(neq, neq, &ctx);
    let ls = SUNLinSol_Dense(&y, &a_mat, &ctx);

    /* Linear solver interface */
    let flag = ARKodeSetLinearSolver(&mut arkode_mem, ls, Some(a_mat));
    assert!(flag >= 0, "ARKodeSetLinearSolver failed with flag = {}", flag);
    let flag = ARKodeSetJacFn(&mut arkode_mem, Some(jac));
    assert!(flag >= 0, "ARKodeSetJacFn failed with flag = {}", flag);

    /* Signal that the problem is autonomous i.e., the RHS does not explicitly
    depend on time.  This is combined with ARKodeSetDeduceImplicitRhs because
    the implicit method is stiffly accurate, so we can reuse the implicit RHS
    evaluation from the beginning of the step in the residual. */
    let flag = ARKodeSetAutonomous(&mut arkode_mem, true);
    assert!(flag >= 0, "ARKodeSetAutonomous failed with flag = {}", flag);

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
    println!("   -------------------------------------------");
    println!(
        "  {}  {}  {}  {}",
        fmt_f(t, 10, 6),
        fmt_f(y.data[0], 10, 6),
        fmt_f(y.data[1], 10, 6),
        fmt_f(y.data[2], 10, 6)
    );

    for _iout in 0..nt {
        let flag = ARKodeEvolve(&mut arkode_mem, tout, &mut y, &mut t, ARK_NORMAL);
        if flag < 0 {
            eprintln!("Solver failure, stopping integration");
            break;
        }
        println!(
            "  {}  {}  {}  {}",
            fmt_f(t, 10, 6),
            fmt_f(y.data[0], 10, 6),
            fmt_f(y.data[1], 10, 6),
            fmt_f(y.data[2], 10, 6)
        );
        let _ = writeln!(
            ufid,
            " {} {} {} {}",
            fmt_e(t, 0, 16),
            fmt_e(y.data[0], 0, 16),
            fmt_e(y.data[1], 0, 16),
            fmt_e(y.data[2], 0, 16)
        );
        /* successful solve: update time */
        tout += dtout;
        tout = if tout > tf { tf } else { tout };
    }
    println!("   -------------------------------------------");
    drop(ufid);

    /* Print some final statistics */
    let (mut nst, mut nst_a, mut nfe, mut nfi) = (0i64, 0i64, 0i64, 0i64);
    let (mut nsetups, mut nje, mut nfe_ls, mut nni) = (0i64, 0i64, 0i64, 0i64);
    let (mut nnf, mut ncfn, mut netf) = (0i64, 0i64, 0i64);
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

    println!("\nFinal Solver Statistics:");
    println!("   Internal solver steps = {} (attempted = {})", nst, nst_a);
    println!("   Total RHS evals:  Fe = {},  Fi = {}", nfe, nfi);
    println!("   Total linear solver setups = {}", nsetups);
    println!("   Total RHS evals for setting up the linear system = {}", nfe_ls);
    println!("   Total number of Jacobian evaluations = {}", nje);
    println!("   Total number of Newton iterations = {}", nni);
    println!("   Total number of nonlinear solver convergence failures = {}", nnf);
    println!("   Total number of error test failures = {}", netf);
    println!("   Total number of failed steps from solver failure = {}", ncfn);

    /* Clean up and return with successful completion */
    drop(y); /* Free y vector */
    let mut slot = Some(arkode_mem);
    ARKodeFree(&mut slot); /* Free integrator memory */
}
