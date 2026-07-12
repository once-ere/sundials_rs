/* -----------------------------------------------------------------
 * Translation of examples/arkode/C_serial/ark_brusselator_mri.c
 * (SUNDIALS 7.7.0).
 *
 * Example problem:
 *
 * The following test simulates a brusselator problem from chemical
 * kinetics. This is an ODE system with 3 components, Y = [u,v,w],
 * satisfying the equations,
 *
 *    du/dt = a - (w+1)*u + v*u^2
 *    dv/dt = w*u - v*u^2
 *    dw/dt = (b-w)/ep - w*u
 *
 * for t in the interval [0.0, 2.0], with parameter values a=1,
 * b=3.5, and ep=1.0e-2. The initial conditions Y0 = [u0,v0,w0] are
 * u0=1.2, v0=3.1, and w0=3.
 *
 * This program solves the problem with the MRI stepper. Outputs are
 * printed at equal intervals of 0.1 and run statistics are printed
 * at the end.
 *
 * Note on ownership: the C example keeps its own inner_arkode_mem
 * pointer alongside the inner stepper; in this port the stepper
 * owns the wrapped inner integrator (and the outer integrator owns
 * the stepper), so the fast-integrator statistics are read by
 * borrowing the inner integrator back out of the outer step memory.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]

use arkode_rs::arkode::{ARKodeCreateMRIStepInnerStepper, ARKodeEvolve, ARKodeFree};
use arkode_rs::arkode_arkstep::ARKStepCreate;
use arkode_rs::arkode_arkstep_io::ARKStepSetTableNum;
use arkode_rs::arkode_butcher_erk::ARKODE_KNOTH_WOLKE_3_3;
use arkode_rs::arkode_io::{
    ARKodeGetNumRhsEvals, ARKodeGetNumSteps, ARKodeSetFixedStep, ARKodeSetUserData,
};
use arkode_rs::arkode_mristep::MRIStepCreate;
use arkode_rs::arkode_mristep_impl::ARKodeMRIStepMem;
use arkode_rs::sundials_utils::{fmt_f, fmt_g};
use arkode_rs::*;
use std::io::Write;

/* ff routine to compute the fast portion of the ODE RHS. */
fn ff(_t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let rdata = user_data.as_mut().unwrap().downcast_mut::<[f64; 3]>().unwrap();
    let b = rdata[1]; /* access data entries */
    let ep = rdata[2];
    let w = y.data[2]; /* access solution values */

    /* fill in the RHS function */
    ydot.data[0] = 0.0;
    ydot.data[1] = 0.0;
    ydot.data[2] = (b - w) / ep;

    /* Return with success */
    0
}

/* fs routine to compute the slow portion of the ODE RHS. */
fn fs(_t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let rdata = user_data.as_mut().unwrap().downcast_mut::<[f64; 3]>().unwrap();
    let a = rdata[0]; /* access data entries */
    let u = y.data[0]; /* access solution values */
    let v = y.data[1];
    let w = y.data[2];

    /* fill in the RHS function */
    ydot.data[0] = a - (w + 1.0) * u + v * u * u;
    ydot.data[1] = w * u - v * u * u;
    ydot.data[2] = -w * u;

    /* Return with success */
    0
}

/* Main Program */
fn main() {
    /* general problem parameters */
    let T0: f64 = 0.0; /* initial time */
    let Tf: f64 = 2.0; /* final time */
    let dTout: f64 = 0.1; /* time between outputs */
    let NEQ: i64 = 3; /* number of dependent vars. */
    let Nt: i32 = (Tf / dTout).ceil() as i32; /* number of output times */
    let hs: f64 = 0.025; /* slow step size */
    let hf: f64 = 0.001; /* fast step size */

    /*
     * Initialization
     */

    /* Create the SUNDIALS context object for this simulation */
    let ctx = SUNContext_Create();

    /* Set up the test problem parameters */
    let a: f64 = 1.0;
    let b: f64 = 3.5;
    let ep: f64 = 1.0e-2;

    /* Set the initial contions */
    let u0: f64 = 1.2;
    let v0: f64 = 3.1;
    let w0: f64 = 3.0;

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
    println!("    hs = {},  hf = {}\n", fmt_g(hs, 0, 6), fmt_g(hf, 0, 6));

    /* Set parameters in user data */
    let rdata: [f64; 3] = [a, b, ep];

    /* Create and initialize serial vector for the solution */
    let mut y = N_VNew_Serial(NEQ, &ctx);
    y.data[0] = u0;
    y.data[1] = v0;
    y.data[2] = w0;

    /*
     * Create the fast integrator and set options
     */

    /* Initialize the fast integrator. */
    let mut inner_arkode_mem = ARKStepCreate(Some(ff), None, T0, &y, &ctx).expect("ARKStepCreate");

    /* Attach user data to fast integrator */
    let retval = ARKodeSetUserData(&mut inner_arkode_mem, Some(Box::new(rdata)));
    assert!(retval >= 0, "ARKodeSetUserData failed with retval = {}", retval);

    /* Set the fast method */
    let retval = ARKStepSetTableNum(&mut inner_arkode_mem, -1, ARKODE_KNOTH_WOLKE_3_3);
    assert!(retval >= 0, "ARKStepSetTableNum failed with retval = {}", retval);

    /* Set the fast step size */
    let retval = ARKodeSetFixedStep(&mut inner_arkode_mem, hf);
    assert!(retval >= 0, "ARKodeSetFixedStep failed with retval = {}", retval);

    /* Create inner stepper */
    let inner_stepper =
        ARKodeCreateMRIStepInnerStepper(inner_arkode_mem).expect("ARKodeCreateMRIStepInnerStepper");

    /*
     * Create the slow integrator and set options
     */

    /* Initialize the slow integrator. */
    let mut arkode_mem =
        MRIStepCreate(Some(fs), None, T0, &y, inner_stepper, &ctx).expect("MRIStepCreate");

    /* Pass rdata to user functions */
    let retval = ARKodeSetUserData(&mut arkode_mem, Some(Box::new(rdata)));
    assert!(retval >= 0, "ARKodeSetUserData failed with retval = {}", retval);

    /* Set the slow step size */
    let retval = ARKodeSetFixedStep(&mut arkode_mem, hs);
    assert!(retval >= 0, "ARKodeSetFixedStep failed with retval = {}", retval);

    /*
     * Integrate ODE
     */

    /* Open output stream for results, output comment line */
    let mut UFID = std::fs::File::create("ark_brusselator_mri_solution.txt").unwrap();
    let _ = writeln!(UFID, "# t u v w");

    /* output initial condition to disk */
    let _ = writeln!(
        UFID,
        " {} {} {} {}",
        fmt_e16(T0),
        fmt_e16(y.data[0]),
        fmt_e16(y.data[1]),
        fmt_e16(y.data[2])
    );

    /* Main time-stepping loop */
    let mut t = T0;
    let mut tout = T0 + dTout;
    println!("        t           u           v           w");
    println!("   ----------------------------------------------");
    println!(
        "  {}  {}  {}  {}",
        fmt_f(t, 10, 6),
        fmt_f(y.data[0], 10, 6),
        fmt_f(y.data[1], 10, 6),
        fmt_f(y.data[2], 10, 6)
    );

    for _iout in 0..Nt {
        /* call integrator */
        let retval = ARKodeEvolve(&mut arkode_mem, tout, &mut y, &mut t, ARK_NORMAL);
        if retval < 0 {
            eprintln!("\nSUNDIALS_ERROR: ARKodeEvolve() failed with retval = {}\n", retval);
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
            UFID,
            " {} {} {} {}",
            fmt_e16(t),
            fmt_e16(y.data[0]),
            fmt_e16(y.data[1]),
            fmt_e16(y.data[2])
        );

        /* successful solve: update time */
        tout += dTout;
        tout = if tout > Tf { Tf } else { tout };
    }
    println!("   ----------------------------------------------");
    drop(UFID);

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
        let inner_arkode_mem = step_mem
            .stepper
            .content
            .as_mut()
            .unwrap()
            .downcast_mut::<ARKodeMem>()
            .unwrap();
        ARKodeGetNumSteps(inner_arkode_mem, &mut nstf);
        ARKodeGetNumRhsEvals(inner_arkode_mem, 0, &mut nff);
    }

    /* Print some final statistics */
    println!("\nFinal Solver Statistics:");
    println!("   Steps: nsts = {}, nstf = {}", nsts, nstf);
    println!("   Total RHS evals:  Fs = {},  Ff = {}", nfse, nff);

    /* Clean up and return */
    drop(y);
    let mut slot = Some(arkode_mem);
    ARKodeFree(&mut slot);
}

/* C "%.16e" */
fn fmt_e16(x: f64) -> String {
    arkode_rs::sundials_utils::fmt_e(x, 0, 16)
}
