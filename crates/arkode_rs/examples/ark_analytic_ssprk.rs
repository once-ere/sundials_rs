/* -----------------------------------------------------------------
 * Translation of examples/arkode/C_serial/ark_analytic_ssprk.c
 * (SUNDIALS 7.7.0).
 *
 * Example problem: the following is a simple example problem with
 * analytical solution,
 *     dy/dt = lambda*y + 1/(1+t^2) - lambda*atan(t)
 * for t in the interval [0.0, 10.0], with initial condition: y=0.
 *
 * The stiffness of the problem is directly proportional to the
 * value of "lambda".  The value of lambda should be negative to
 * result in a well-posed ODE; for values with magnitude larger
 * than 100 the problem becomes quite stiff.
 *
 * This program solves the problem with the SSPRK method.  Output is
 * printed every 1.0 units of time (10 total).  Run statistics
 * (optional outputs) are printed at the end.
 * -----------------------------------------------------------------*/

use arkode_rs::arkode::{ARKodeEvolve, ARKodeFree, ARKodeSStolerances};
use arkode_rs::arkode_io::{ARKodePrintAllStats, ARKodeSetUserData};
use arkode_rs::arkode_lsrkstep::LSRKStepCreateSSP;
use arkode_rs::arkode_lsrkstep_impl::ARKODE_LSRK_SSP_S_3;
use arkode_rs::arkode_lsrkstep_io::{LSRKStepSetNumSSPStages, LSRKStepSetSSPMethod};
use arkode_rs::sundials_utils::{fmt_e, fmt_f, fmt_g};
use arkode_rs::*;
use std::io::Write;

/* f routine to compute the ODE RHS function f(t,y). */
fn f(t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let rdata = user_data.as_mut().unwrap().downcast_mut::<[f64; 1]>().unwrap();
    let lambda = rdata[0]; /* set shortcut for stiffness parameter */
    let u = y.data[0]; /* access current solution value */

    /* fill in the RHS function */
    ydot.data[0] = lambda * u + 1.0 / (1.0 + t * t) - lambda * t.atan();

    0 /* return with success */
}

/* check the error */
fn compute_error(y: &NVector, t: f64) -> i32 {
    /* compute solution error */
    let ans = t.atan();
    let err = SUNRabs(y.data[0] - ans);

    println!("\nACCURACY at the final time = {}", fmt_g(err, 0, 6));
    0
}

fn main() {
    /* general problem parameters */
    let t0: f64 = 0.0; /* initial time */
    let tf: f64 = 10.0; /* final time */
    let dtout: f64 = 1.0; /* time between outputs */
    let neq: i64 = 1; /* number of dependent vars. */

    let reltol: f64 = 1.0e-8; /* tolerances */
    let abstol: f64 = 1.0e-8;
    let lambda: f64 = -10.0; /* stiffness parameter */

    /* Create the SUNDIALS context object for this simulation */
    let ctx = SUNContext_Create();

    /* Initial diagnostics output */
    println!("\nAnalytical ODE test problem:");
    println!("    lambda = {}", fmt_g(lambda, 0, 6));
    println!("   reltol = {}", fmt_e(reltol, 0, 1));
    println!("   abstol = {}\n", fmt_e(abstol, 0, 1));

    /* Initialize data structures */
    let mut y = N_VNew_Serial(neq, &ctx); /* Create serial vector for solution */
    N_VConst(0.0, &mut y); /* Specify initial condition */

    /* Call LSRKStepCreateSSP to initialize the ARK timestepper module
    and specify the right-hand side function in y'=f(t,y), the initial
    time T0, and the initial dependent variable vector y. */
    let mut arkode_mem = LSRKStepCreateSSP(Some(f), t0, &y, &ctx).expect("LSRKStepCreateSSP");

    /* Set routines */
    let flag = ARKodeSetUserData(&mut arkode_mem, Some(Box::new([lambda])));
    assert!(flag >= 0, "ARKodeSetUserData failed with flag = {}", flag);

    /* Specify tolerances */
    let flag = ARKodeSStolerances(&mut arkode_mem, reltol, abstol);
    assert!(flag >= 0, "ARKodeSStolerances failed with flag = {}", flag);

    /* Specify the SSP method */
    let flag = LSRKStepSetSSPMethod(&mut arkode_mem, ARKODE_LSRK_SSP_S_3);
    assert!(flag >= 0, "LSRKStepSetSSPMethod failed with flag = {}", flag);

    /* Specify the number of SSP stages */
    let flag = LSRKStepSetNumSSPStages(&mut arkode_mem, 9);
    assert!(flag >= 0, "LSRKStepSetNumSSPStages failed with flag = {}", flag);

    /* Open output stream for results, output comment line */
    let mut ufid = std::fs::File::create("solution.txt").unwrap();
    let _ = writeln!(ufid, "# t u");

    /* output initial condition to disk */
    let _ = writeln!(ufid, " {} {}", fmt_e(t0, 0, 16), fmt_e(y.data[0], 0, 16));

    /* Main time-stepping loop: calls ARKodeEvolve to perform the
    integration, then prints results.  Stops when the final time has
    been reached */
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
        println!("  {}  {}", fmt_f(t, 10, 6), fmt_f(y.data[0], 10, 6));
        let _ = writeln!(ufid, " {} {}", fmt_e(t, 0, 16), fmt_e(y.data[0], 0, 16));
        /* successful solve: update time */
        tout += dtout;
        tout = if tout > tf { tf } else { tout };
    }
    println!("   ---------------------");
    drop(ufid);

    /* Print final statistics */
    println!("\nFinal Statistics:");
    let mut stdout = std::io::stdout();
    ARKodePrintAllStats(&mut arkode_mem, &mut stdout, SUN_OUTPUTFORMAT_TABLE);

    /* Print final statistics to a file in CSV format */
    let mut fid = std::fs::File::create("ark_analytic_nonlin_stats.csv").unwrap();
    ARKodePrintAllStats(&mut arkode_mem, &mut fid, SUN_OUTPUTFORMAT_CSV);
    drop(fid);

    /* check the solution error */
    let flag = compute_error(&y, t);

    /* Clean up and return */
    drop(y);
    let mut slot = Some(arkode_mem);
    ARKodeFree(&mut slot);

    std::process::exit(flag);
}
