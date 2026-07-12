/*---------------------------------------------------------------
 * Translation of examples/arkode/C_serial/ark_heat1D.c
 * (SUNDIALS 7.7.0).
 *
 * Example problem:
 *
 * The following test simulates a simple 1D heat equation,
 *    u_t = k*u_xx + f
 * for t in [0, 10], x in [0, 1], with initial conditions
 *    u(0,x) =  0
 * Dirichlet boundary conditions, i.e.
 *    u_t(t,0) = u_t(t,1) = 0,
 * and a point-source heating term,
 *    f = 0.01 for x=0.5.
 *
 * The spatial derivatives are computed using second-order
 * centered differences, with the data distributed over N points
 * on a uniform spatial grid.
 *
 * This program solves the problem with either an ERK or DIRK
 * method.  For the DIRK method, we use a Newton iteration with
 * the SUNLinSol_PCG linear solver, and a user-supplied
 * Jacobian-vector product routine.
 *
 * 100 outputs are printed at equal intervals, and run statistics
 * are printed at the end.
 *
 * Translation note: C's trailing SUNLinSolSetOptions(LS, ...) call
 * is omitted — the Rust PCG solver carries no CLI option plumbing
 * (and the reference run passes no arguments).
 *---------------------------------------------------------------*/
#![allow(non_snake_case)]

use std::io::Write;

use arkode_rs::arkode::{ARKodeEvolve, ARKodeFree, ARKodeSStolerances};
use arkode_rs::arkode_arkstep::ARKStepCreate;
use arkode_rs::arkode_cli::ARKodeSetOptions;
use arkode_rs::arkode_io::{
    ARKodeGetNumErrTestFails, ARKodeGetNumLinSolvSetups, ARKodeGetNumNonlinSolvConvFails,
    ARKodeGetNumNonlinSolvIters, ARKodeGetNumRhsEvals, ARKodeGetNumStepAttempts,
    ARKodeGetNumSteps, ARKodeSetLinear, ARKodeSetMaxNumSteps, ARKodeSetPredictorMethod, ARKodeSetUserData,
};
use arkode_rs::arkode_ls::{
    ARKodeGetNumJtimesEvals, ARKodeGetNumLinConvFails, ARKodeGetNumLinIters,
    ARKodeSetJacTimes, ARKodeSetLinearSolver,
};
use arkode_rs::sundials_utils::{fmt_e, fmt_f, fmt_g};
use arkode_rs::*;

/* user data structure */
struct HeatData {
    N: i64,  /* number of intervals   */
    dx: f64, /* mesh spacing          */
    k: f64,  /* diffusion coefficient */
}

/*--------------------------------
 * Functions called by the solver
 *--------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
fn f(_t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let udata = user_data.as_ref().unwrap().downcast_ref::<HeatData>().unwrap();
    let N = udata.N; /* set variable shortcuts */
    let k = udata.k;
    let dx = udata.dx;

    N_VConst(0.0, ydot); /* Initialize ydot to zero */

    /* iterate over domain, computing all equations */
    let c1 = k / dx / dx;
    let c2 = -2.0 * k / dx / dx;
    let isource = (N / 2) as usize;
    ydot.data[0] = 0.0; /* left boundary condition */
    for i in 1..(N - 1) as usize {
        ydot.data[i] = c1 * y.data[i - 1] + c2 * y.data[i] + c1 * y.data[i + 1];
    }
    ydot.data[(N - 1) as usize] = 0.0; /* right boundary condition */
    ydot.data[isource] += 0.01 / dx; /* source term */

    0 /* Return with success */
}

/* Jacobian routine to compute J(t,y) = df/dy. */
fn jac(
    v: &NVector,
    jv: &mut NVector,
    _t: f64,
    _y: &NVector,
    _fy: &NVector,
    user_data: &mut UserData,
    _tmp: &mut NVector,
) -> i32 {
    let udata = user_data.as_ref().unwrap().downcast_ref::<HeatData>().unwrap();
    let N = udata.N; /* variable shortcuts */
    let k = udata.k;
    let dx = udata.dx;

    N_VConst(0.0, jv); /* initialize Jv product to zero */

    /* iterate over domain, computing all Jacobian-vector products */
    let c1 = k / dx / dx;
    let c2 = -2.0 * k / dx / dx;
    jv.data[0] = 0.0;
    for i in 1..(N - 1) as usize {
        jv.data[i] = c1 * v.data[i - 1] + c2 * v.data[i] + c1 * v.data[i + 1];
    }
    jv.data[(N - 1) as usize] = 0.0;

    0 /* Return with success */
}

/*-------------------------------
 * Private helper functions
 *-------------------------------*/

/* Check if a SUNDIALS function returned a negative flag */
fn check_flag(flag: i32, funcname: &str) -> bool {
    if flag < 0 {
        eprintln!("\nSUNDIALS_ERROR: {}() failed with flag = {}\n", funcname, flag);
        return true;
    }
    false
}

/* Main Program */
fn main() {
    /* general problem parameters */
    let T0: f64 = 0.0; /* initial time */
    let Tf: f64 = 1.0; /* final time */
    let Nt: i32 = 10; /* total number of output times */
    let rtol: f64 = 1.0e-6; /* relative tolerance */
    let atol: f64 = 1.0e-10; /* absolute tolerance */
    let N: i64 = 201; /* spatial mesh size */
    let k: f64 = 0.5; /* heat conductivity */

    /* Create the SUNDIALS context object for this simulation */
    let ctx = SUNContext_Create();

    /* allocate and fill udata structure */
    let udata = HeatData {
        N,
        k,
        dx: 1.0 / (N - 1) as f64, /* mesh spacing */
    };
    let (u_n, u_k, u_dx) = (udata.N, udata.k, udata.dx);

    /* Initial problem output */
    println!("\n1D Heat PDE test problem:");
    println!("  N = {}", u_n);
    println!("  diffusion coefficient:  k = {}", fmt_g(u_k, 0, 6));

    /* Initialize data structures */
    let mut y = N_VNew_Serial(N, &ctx); /* Create serial vector for solution */
    N_VConst(0.0, &mut y); /* Set initial conditions */

    /* Call ARKStepCreate to initialize the ARK timestepper module and
       specify the right-hand side function in y'=f(t,y), the initial time
       T0, and the initial dependent variable vector y.  Note: since this
       problem is fully implicit, we set f_E to NULL and f_I to f. */
    let mut arkode_mem = ARKStepCreate(None, Some(f), T0, &y, &ctx).expect("ARKStepCreate");

    /* Set routines */
    let flag = ARKodeSetUserData(&mut arkode_mem, Some(Box::new(udata)));
    if check_flag(flag, "ARKodeSetUserData") {
        return;
    }
    let flag = ARKodeSetMaxNumSteps(&mut arkode_mem, 10000); /* Increase max num steps */
    if check_flag(flag, "ARKodeSetMaxNumSteps") {
        return;
    }
    let flag = ARKodeSetPredictorMethod(&mut arkode_mem, 1); /* Specify maximum-order predictor */
    if check_flag(flag, "ARKodeSetPredictorMethod") {
        return;
    }
    let flag = ARKodeSStolerances(&mut arkode_mem, rtol, atol); /* Specify tolerances */
    if check_flag(flag, "ARKodeSStolerances") {
        return;
    }

    /* Initialize PCG solver -- no preconditioning, with up to N iterations */
    let ls = SUNLinSol_PCG(&y, SUN_PREC_NONE, N as i32, &ctx);

    /* Linear solver interface -- set user-supplied J*v routine (no 'jtsetup'
       required) */
    let flag = ARKodeSetLinearSolver(&mut arkode_mem, ls, None); /* Attach linear solver */
    if check_flag(flag, "ARKodeSetLinearSolver") {
        return;
    }
    let flag = ARKodeSetJacTimes(&mut arkode_mem, None, Some(jac)); /* Set the Jacobian routine */
    if check_flag(flag, "ARKodeSetJacTimes") {
        return;
    }

    /* Specify linearly implicit RHS, with non-time-dependent Jacobian */
    let flag = ARKodeSetLinear(&mut arkode_mem, 0);
    if check_flag(flag, "ARKodeSetLinear") {
        return;
    }

    /* Override any current settings with command-line options */
    let args: Vec<String> = std::env::args().collect();
    let flag = ARKodeSetOptions(&mut arkode_mem, None, None, &args);
    if check_flag(flag, "ARKodeSetOptions") {
        return;
    }
    /* (SUNLinSolSetOptions: no PCG CLI option plumbing in this build) */

    /* output mesh to disk */
    let mut fid = std::fs::File::create("heat_mesh.txt").expect("fopen");
    for i in 0..N {
        let _ = writeln!(fid, "  {}", fmt_e(u_dx * i as f64, 0, 16));
    }
    drop(fid);

    /* Open output stream for results */
    let mut ufid = std::fs::File::create("heat1D.txt").expect("fopen");

    /* output initial condition to disk */
    for i in 0..N as usize {
        let _ = write!(ufid, " {}", fmt_e(y.data[i], 0, 16));
    }
    let _ = writeln!(ufid);

    /* Main time-stepping loop: calls ARKodeEvolve to perform the integration,
       then prints results.  Stops when the final time has been reached */
    let mut t = T0;
    let dTout = (Tf - T0) / Nt as f64;
    let mut tout = T0 + dTout;
    println!("        t      ||u||_rms");
    println!("   -------------------------");
    println!(
        "  {}  {}",
        fmt_f(t, 10, 6),
        fmt_f((N_VDotProd(&y, &y) / N as f64).sqrt(), 10, 6)
    );
    for _iout in 0..Nt {
        let flag = ARKodeEvolve(&mut arkode_mem, tout, &mut y, &mut t, ARK_NORMAL); /* call integrator */
        if check_flag(flag, "ARKodeEvolve") {
            break;
        }
        println!(
            "  {}  {}",
            fmt_f(t, 10, 6),
            fmt_f((N_VDotProd(&y, &y) / N as f64).sqrt(), 10, 6)
        ); /* print solution stats */
        if flag >= 0 {
            /* successful solve: update output time */
            tout += dTout;
            tout = if tout > Tf { Tf } else { tout };
        } else {
            /* unsuccessful solve: break */
            eprintln!("Solver failure, stopping integration");
            break;
        }

        /* output results to disk */
        for i in 0..N as usize {
            let _ = write!(ufid, " {}", fmt_e(y.data[i], 0, 16));
        }
        let _ = writeln!(ufid);
    }
    println!("   -------------------------");
    drop(ufid);

    /* Print some final statistics */
    let mut nst: i64 = 0;
    let mut nst_a: i64 = 0;
    let mut nfe: i64 = 0;
    let mut nfi: i64 = 0;
    let mut nsetups: i64 = 0;
    let mut nli: i64 = 0;
    let mut nJv: i64 = 0;
    let mut nlcf: i64 = 0;
    let mut nni: i64 = 0;
    let mut ncfn: i64 = 0;
    let mut netf: i64 = 0;

    let flag = ARKodeGetNumSteps(&mut arkode_mem, &mut nst);
    check_flag(flag, "ARKodeGetNumSteps");
    let flag = ARKodeGetNumStepAttempts(&mut arkode_mem, &mut nst_a);
    check_flag(flag, "ARKodeGetNumStepAttempts");
    let flag = ARKodeGetNumRhsEvals(&mut arkode_mem, 0, &mut nfe);
    check_flag(flag, "ARKodeGetNumRhsEvals");
    let flag = ARKodeGetNumRhsEvals(&mut arkode_mem, 1, &mut nfi);
    check_flag(flag, "ARKodeGetNumRhsEvals");
    let flag = ARKodeGetNumLinSolvSetups(&mut arkode_mem, &mut nsetups);
    check_flag(flag, "ARKodeGetNumLinSolvSetups");
    let flag = ARKodeGetNumErrTestFails(&mut arkode_mem, &mut netf);
    check_flag(flag, "ARKodeGetNumErrTestFails");
    let flag = ARKodeGetNumNonlinSolvIters(&mut arkode_mem, &mut nni);
    check_flag(flag, "ARKodeGetNumNonlinSolvIters");
    let flag = ARKodeGetNumNonlinSolvConvFails(&mut arkode_mem, &mut ncfn);
    check_flag(flag, "ARKodeGetNumNonlinSolvConvFails");
    let flag = ARKodeGetNumLinIters(&mut arkode_mem, &mut nli);
    check_flag(flag, "ARKodeGetNumLinIters");
    let flag = ARKodeGetNumJtimesEvals(&mut arkode_mem, &mut nJv);
    check_flag(flag, "ARKodeGetNumJtimesEvals");
    let flag = ARKodeGetNumLinConvFails(&mut arkode_mem, &mut nlcf);
    check_flag(flag, "ARKodeGetNumLinConvFails");

    println!("\nFinal Solver Statistics:");
    println!("   Internal solver steps = {} (attempted = {})", nst, nst_a);
    println!("   Total RHS evals:  Fe = {},  Fi = {}", nfe, nfi);
    println!("   Total linear solver setups = {}", nsetups);
    println!("   Total linear iterations = {}", nli);
    println!("   Total number of Jacobian-vector products = {}", nJv);
    println!("   Total number of linear solver convergence failures = {}", nlcf);
    println!("   Total number of Newton iterations = {}", nni);
    println!("   Total number of nonlinear solver convergence failures = {}", ncfn);
    println!("   Total number of error test failures = {}", netf);

    /* Clean up and return with successful completion */
    drop(y); /* Free vectors */
    let mut slot = Some(arkode_mem);
    ARKodeFree(&mut slot); /* Free integrator memory */
}
