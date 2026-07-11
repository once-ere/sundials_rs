/* -----------------------------------------------------------------
 * Translation of examples/arkode/C_serial/ark_analytic_nonlin.c
 * (SUNDIALS 7.7.0).
 *
 * Example problem: the following is a simple example problem with
 * analytical solution,
 *     dy/dt = (t+1)*exp(-y)
 * for t in the interval [0.0, 10.0], with initial condition: y=0.
 * This has analytical solution y(t) = log(0.5*t^2 + t + 1).
 *
 * This program solves the problem with the ERK method.
 * Output is printed every 1.0 units of time (10 total).
 * Run statistics (optional outputs) are printed at the end.
 * -----------------------------------------------------------------*/

use arkode_rs::arkode::{ARKodeEvolve, ARKodeFree, ARKodeSStolerances};
use arkode_rs::arkode_erkstep::ERKStepCreate;
use arkode_rs::arkode_io::ARKodePrintAllStats;
use arkode_rs::sundials_utils::{fmt_e, fmt_f};
use arkode_rs::*;
use std::io::Write;

/* f routine to compute the ODE RHS function f(t,y). */
fn f(t: f64, y: &NVector, ydot: &mut NVector, _user_data: &mut UserData) -> i32 {
    ydot.data[0] = (t + 1.0) * SUNRexp(-y.data[0]);
    0
}

fn main() {
    /* general problem parameters */
    let t0: f64 = 0.0; /* initial time */
    let tf: f64 = 10.0; /* final time */
    let dtout: f64 = 1.0; /* time between outputs */
    let neq: i64 = 1; /* number of dependent vars. */
    let reltol: f64 = 1.0e-6; /* tolerances */
    let abstol: f64 = 1.0e-10;

    /* Create the SUNDIALS context object for this simulation */
    let ctx = SUNContext_Create();

    /* Initial problem output */
    println!("\nAnalytical ODE test problem:");
    println!("   reltol = {}", fmt_e(reltol, 0, 1));
    println!("   abstol = {}\n", fmt_e(abstol, 0, 1));

    /* Initialize data structures */
    let mut y = N_VNew_Serial(neq, &ctx); /* Create serial vector for solution */
    y.data[0] = 0.0; /* Specify initial condition */

    /* Call ERKStepCreate to initialize the ERK timestepper module and
    specify the right-hand side function in y'=f(t,y), the initial time
    T0, and the initial dependent variable vector y. */
    let mut arkode_mem = ERKStepCreate(f, t0, &y, &ctx);

    /* Specify tolerances */
    let flag = ARKodeSStolerances(&mut arkode_mem, reltol, abstol);
    assert!(flag >= 0, "ARKodeSStolerances failed with flag = {}", flag);

    /* Open output stream for results, output comment line */
    let mut ufid = std::fs::File::create("solution.txt").unwrap();
    let _ = writeln!(ufid, "# t u");

    /* output initial condition to disk */
    let _ = writeln!(ufid, " {} {}", fmt_e(t0, 0, 16), fmt_e(y.data[0], 0, 16));

    /* Main time-stepping loop: calls ARKodeEvolve to perform the
    integration, then prints results.  Stops when the final time has been
    reached */
    let mut t = t0;
    let mut tout = t0 + dtout;
    println!("        t           u");
    println!("   ---------------------");
    while tf - t > 1.0e-15 {
        let flag = ARKodeEvolve(&mut arkode_mem, tout, &mut y, &mut t, ARK_NORMAL);
        if flag < 0 {
            eprintln!("\nSUNDIALS_ERROR: ARKodeEvolve() failed with flag = {}\n", flag);
            break;
        }
        /* access/print solution */
        println!("  {}  {}", fmt_f(t, 10, 6), fmt_f(y.data[0], 10, 6));
        let _ = writeln!(ufid, " {} {}", fmt_e(t, 0, 16), fmt_e(y.data[0], 0, 16));
        if flag >= 0 {
            /* successful solve: update time */
            tout += dtout;
            tout = if tout > tf { tf } else { tout };
        }
    }
    println!("   ---------------------");
    drop(ufid);

    /* Print final statistics */
    println!("\nFinal Statistics:");
    let mut stdout = std::io::stdout();
    let _ = ARKodePrintAllStats(&mut arkode_mem, &mut stdout, SUN_OUTPUTFORMAT_TABLE);

    /* Print final statistics to a file in CSV format */
    let mut fid = std::fs::File::create("ark_analytic_nonlin_stats.csv").unwrap();
    let _ = ARKodePrintAllStats(&mut arkode_mem, &mut fid, SUN_OUTPUTFORMAT_CSV);
    drop(fid);

    /* Clean up and return with successful completion */
    drop(y); /* Free y vector */
    let mut slot = Some(arkode_mem);
    ARKodeFree(&mut slot); /* Free integrator memory */
}
