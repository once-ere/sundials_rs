/* ----------------------------------------------------------------
 * Translation of examples/arkode/C_serial/ark_twowaycouple_mri.c
 * (SUNDIALS 7.7.0).
 *
 * Based a linear example program by Rujeko Chinomona @ UMBC.
 *
 * Example problem:
 *
 * This example simulates an ODE system with 3 components,
 * Y = [u,v,w], given by the equations,
 *
 *   du/dt =  100v+w
 *   dv/dt = -100u
 *   dw/dt = -w+u
 *
 * for t in the interval [0.0, 2.0] with initial conditions
 * u(0)=9001/10001, v(0)=-1e-5/10001, and w(0)=1000. In this problem
 * the slow (w) and fast (u and v) components depend on one another.
 *
 * This program solves the problem with the MRI stepper. Outputs are
 * printed at equal intervals of 0.1 and run statistics are printed
 * at the end.
 *
 * Note on ownership: the stepper owns the wrapped inner integrator
 * (and the outer integrator owns the stepper), so the final fast
 * statistics are read by borrowing the inner integrator back out of
 * the outer step memory.
 * ----------------------------------------------------------------*/

use std::io::Write;

use arkode_rs::arkode::{ARKodeCreateMRIStepInnerStepper, ARKodeEvolve, ARKodeFree};
use arkode_rs::arkode_arkstep::ARKStepCreate;
use arkode_rs::arkode_arkstep_io::ARKStepSetTableNum;
use arkode_rs::arkode_butcher_erk::ARKODE_KNOTH_WOLKE_3_3;
use arkode_rs::arkode_io::{ARKodeGetNumRhsEvals, ARKodeGetNumSteps, ARKodeSetFixedStep};
use arkode_rs::arkode_mristep::MRIStepCreate;
use arkode_rs::arkode_mristep_impl::ARKodeMRIStepMem;
use arkode_rs::sundials_utils::{fmt_e, fmt_f, fmt_g};
use arkode_rs::*;

/* ------------------------------
 * Functions called by the solver
 * ------------------------------*/

/* ff routine to compute the fast portion of the ODE RHS. */
fn ff(_t: f64, y: &NVector, ydot: &mut NVector, _user_data: &mut UserData) -> i32 {
    let c1: f64 = 100.0; /* problem constant */
    let u = y.data[0]; /* access solution values */
    let v = y.data[1];

    /* fill in the RHS function */
    ydot.data[0] = c1 * v;
    ydot.data[1] = -c1 * u;
    ydot.data[2] = u;

    /* Return with success */
    0
}

/* fs routine to compute the slow portion of the ODE RHS. */
fn fs(_t: f64, y: &NVector, ydot: &mut NVector, _user_data: &mut UserData) -> i32 {
    let w = y.data[2]; /* access solution values */

    /* fill in the RHS function */
    ydot.data[0] = w;
    ydot.data[1] = 0.0;
    ydot.data[2] = -w;

    /* Return with success */
    0
}

/* ------------------------------
 * Private helper functions
 * ------------------------------*/

/* Check if a SUNDIALS function returned a negative flag */
fn check_retval(retval: i32, funcname: &str) -> bool {
    if retval < 0 {
        eprintln!("\nSUNDIALS_ERROR: {}() failed with retval = {}\n", funcname, retval);
        return true;
    }
    false
}

/* Main Program */
fn main() {
    /* general problem parameters */
    let t0: f64 = 0.0; /* initial time */
    let tf: f64 = 2.0; /* final time */
    let d_tout: f64 = 0.1; /* time between outputs */
    let neq: i64 = 3; /* number of dependent vars. */
    let nt: i32 = (tf / d_tout).ceil() as i32; /* number of output times */
    let hs: f64 = 0.001; /* slow step size */
    let hf: f64 = 0.00002; /* fast step size */

    /* Create the SUNDIALS context object for this simulation */
    let ctx = SUNContext_Create();

    /*
     * Initialization
     */

    /* Set the initial contions */
    let u0: f64 = 9001.0 / 10001.0;
    let v0: f64 = -1.0e5 / 10001.0;
    let w0: f64 = 1000.0;

    /* Initial problem output */
    println!("\nTwo way coupling ODE test problem:");
    println!(
        "    initial conditions:  u0 = {},  v0 = {},  w0 = {}",
        fmt_g(u0, 0, 6),
        fmt_g(v0, 0, 6),
        fmt_g(w0, 0, 6)
    );
    println!("    hs = {},  hf = {}\n", fmt_g(hs, 0, 6), fmt_g(hf, 0, 6));

    /* Create and initialize serial vector for the solution */
    let mut y = N_VNew_Serial(neq, &ctx);
    y.data[0] = u0;
    y.data[1] = v0;
    y.data[2] = w0;

    /*
     * Create the fast integrator and set options
     */

    /* Initialize the fast integrator. Specify the explicit fast right-hand
       side function in y'=fe(t,y)+fi(t,y)+ff(t,y), the initial time T0, and
       the initial dependent variable vector y. */
    let mut inner_arkode_mem = ARKStepCreate(Some(ff), None, t0, &y, &ctx).expect("ARKStepCreate");

    /* Set the fast method */
    let retval = ARKStepSetTableNum(&mut inner_arkode_mem, -1, ARKODE_KNOTH_WOLKE_3_3);
    if check_retval(retval, "ARKStepSetTableNum") {
        return;
    }

    /* Set the fast step size */
    let retval = ARKodeSetFixedStep(&mut inner_arkode_mem, hf);
    if check_retval(retval, "ARKodeSetFixedStep") {
        return;
    }

    /* Create inner stepper */
    let inner_stepper =
        ARKodeCreateMRIStepInnerStepper(inner_arkode_mem).expect("ARKodeCreateMRIStepInnerStepper");

    /*
     * Create the slow integrator and set options
     */

    /* Initialize the slow integrator. Specify the explicit slow right-hand
       side function in y'=fe(t,y)+fi(t,y)+ff(t,y), the initial time T0, the
       initial dependent variable vector y, and the fast integrator. */
    let mut arkode_mem =
        MRIStepCreate(Some(fs), None, t0, &y, inner_stepper, &ctx).expect("MRIStepCreate");

    /* Set the slow step size */
    let retval = ARKodeSetFixedStep(&mut arkode_mem, hs);
    if check_retval(retval, "ARKodeSetFixedStep") {
        return;
    }

    /*
     * Integrate ODE
     */

    /* Open output stream for results, output comment line */
    let mut ufid = std::fs::File::create("ark_twowaycouple_mri_solution.txt").expect("fopen");
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
       integration, then prints results. Stops when the final time
       has been reached */
    let mut t = t0;
    let mut tout = t0 + d_tout;
    println!("        t           u           v           w");
    println!("   -----------------------------------------------");
    println!(
        "  {}  {}  {}  {}",
        fmt_f(t, 10, 6),
        fmt_f(y.data[0], 10, 6),
        fmt_f(y.data[1], 10, 6),
        fmt_f(y.data[2], 10, 6)
    );

    for _iout in 0..nt {
        /* call integrator */
        let retval = ARKodeEvolve(&mut arkode_mem, tout, &mut y, &mut t, ARK_NORMAL);
        if check_retval(retval, "ARKodeEvolve") {
            break;
        }

        /* access/print solution */
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
        tout += d_tout;
        tout = if tout > tf { tf } else { tout };
    }
    println!("   -----------------------------------------------");
    drop(ufid);

    /*
     * Finalize
     */

    /* Get some slow integrator statistics */
    let mut nsts: i64 = 0;
    let mut nfse: i64 = 0;
    ARKodeGetNumSteps(&mut arkode_mem, &mut nsts);
    ARKodeGetNumRhsEvals(&mut arkode_mem, 0, &mut nfse);

    /* Get some fast integrator statistics (borrow the inner integrator
       back out of the outer step memory) */
    let mut nstf: i64 = 0;
    let mut nff: i64 = 0;
    {
        let step_mem = arkode_mem
            .step_mem
            .as_mut()
            .unwrap()
            .downcast_mut::<ARKodeMRIStepMem>()
            .unwrap();
        let inner = step_mem
            .stepper
            .content
            .as_mut()
            .unwrap()
            .downcast_mut::<ARKodeMem>()
            .unwrap();
        ARKodeGetNumSteps(inner, &mut nstf);
        ARKodeGetNumRhsEvals(inner, 0, &mut nff);
    }

    /* Print some final statistics */
    println!("\nFinal Solver Statistics:");
    println!("   Steps: nsts = {}, nstf = {}", nsts, nstf);
    println!("   Total RHS evals:  Fs = {},  Ff = {}", nfse, nff);

    /* Clean up and return */
    drop(y); /* Free y vector */
    let mut slot = Some(arkode_mem);
    ARKodeFree(&mut slot); /* Free integrator memory */
}
