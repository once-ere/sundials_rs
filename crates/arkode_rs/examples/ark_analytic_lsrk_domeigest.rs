/* -----------------------------------------------------------------
 * Translation of examples/arkode/C_serial/ark_analytic_lsrk_domeigest.c
 * (SUNDIALS 7.7.0).
 *
 * Uses a SUNDomEigEstimator (power iteration) instead of a
 * user-supplied dominant eigenvalue function.
 *
 * Translation notes:
 * - The C initial eigenvector is filled with rand()/RAND_MAX values;
 *   for NEQ = 1 the power iteration normalizes any positive entry to
 *   exactly 1, so a fixed positive value is used here.
 * - C calls SUNDomEigEstimator_SetOptions on the DEE after attaching
 *   it to LSRKStep; here the options are applied before the attach
 *   (the DEE moves into the integrator).  The observable state at the
 *   first Estimate call is identical.
 * - C also calls ARKodeSetOptions(arkode_mem, NULL, ...): its option
 *   prefix is "arkode.", so the "arkid."-prefixed argument of the
 *   reference variant is skipped (by C as well).
 *
 * Example problem: the following is a simple example problem with
 * analytical solution,
 *     dy/dt = (lambda - alpha*cos((10 - t)/10*pi))*y
 *             + 1/(1+t^2) - (lambda - alpha*cos((10 - t)/10*pi))*atan(t)
 * for t in the interval [0.0, 10.0], with initial condition: y=0.
 *
 * The stiffness of the problem is directly proportional to
 * "lambda - alpha*cos((10 - t)/10*pi)", which varies in time.
 *
 * This program solves the problem with the LSRK method.  Output is
 * printed every 1.0 units of time (10 total).  Run statistics
 * (optional outputs) are printed at the end.
 * -----------------------------------------------------------------*/

use arkode_rs::arkode::{ARKodeEvolve, ARKodeFree, ARKodeSStolerances};
use arkode_rs::arkode_io::{ARKodePrintAllStats, ARKodeSetMaxNumSteps, ARKodeSetUserData};
use arkode_rs::arkode_lsrkstep::LSRKStepCreateSTS;
use arkode_rs::arkode_lsrkstep_io::{
    LSRKStepSetDomEigEstimator, LSRKStepSetDomEigFrequency, LSRKStepSetDomEigSafetyFactor,
    LSRKStepSetMaxNumStages, LSRKStepSetNumDomEigEstInitPreprocessIters,
    LSRKStepSetNumDomEigEstPreprocessIters, LSRKStepSetSTSMethodByName,
};
use arkode_rs::arkode_cli::ARKodeSetOptions;
use arkode_rs::sundials_domeigestimator::SUNDomEigEstimator_SetOptions;
use arkode_rs::sundomeigest_power::SUNDomEigEstimator_Power;
use arkode_rs::sundials_utils::{fmt_e, fmt_f, fmt_g};
use arkode_rs::*;
use std::io::Write;

/* f routine to compute the ODE RHS function f(t,y). */
fn f(t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let rdata = user_data.as_mut().unwrap().downcast_mut::<[f64; 2]>().unwrap();
    let lambda = rdata[0]; /* set shortcut for stiffness parameter 1 */
    let alpha = rdata[1]; /* set shortcut for stiffness parameter 2 */
    let u = y.data[0]; /* access current solution value */

    /* fill in the RHS function */
    ydot.data[0] = (lambda - alpha * ((10.0 - t) / 10.0 * (-1.0f64).acos()).cos()) * u
        + 1.0 / (1.0 + t * t)
        - (lambda - alpha * ((10.0 - t) / 10.0 * (-1.0f64).acos()).cos()) * t.atan();

    0 /* return with success */
}

/* check the computed solution */
fn check_ans(y: &NVector, t: f64, rtol: f64, atol: f64) -> i32 {
    /* compute solution error */
    let ans = t.atan();
    let ewt = 1.0 / (rtol * SUNRabs(ans) + atol);
    let err = ewt * SUNRabs(y.data[0] - ans);

    /* is the solution within the tolerances? */
    let passfail = if err < 1.0 { 0 } else { 1 };

    if passfail != 0 {
        println!("\nSUNDIALS_WARNING: check_ans error={}\n", fmt_g(err, 0, 6));
    }

    passfail
}

/* check the error */
fn compute_error(y: &NVector, t: f64) -> i32 {
    /* compute solution error */
    let ans = t.atan();
    let err = SUNRabs(y.data[0] - ans);

    println!("\nACCURACY at the final time   = {}", fmt_g(err, 0, 6));
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
    let lambda: f64 = -1.0e+6; /* stiffness parameter 1 */
    let alpha: f64 = 1.0e+2; /* stiffness parameter 2 */

    let user_data_arr: [f64; 2] = [lambda, alpha];

    /* Create the SUNDIALS context object for this simulation */
    let ctx = SUNContext_Create();

    /* Initial diagnostics output */
    /* Retrieve the command-line options (forwarded to the DEE) */
    let argv: Vec<String> = std::env::args().collect();

    print!("\nAnalytical ODE test problem with a variable Jacobian:");
    print!("\nThe stiffness of the problem is directly proportional to");
    print!("\n\"lambda - alpha*cos((10 - t)/10*pi)\"\n\n");
    println!("    lambda = {}", fmt_g(lambda, 0, 6));
    println!("     alpha = {}", fmt_g(alpha, 0, 6));
    println!("    reltol = {}", fmt_e(reltol, 0, 1));
    println!("    abstol = {}\n", fmt_e(abstol, 0, 1));

    /* Initialize data structures */
    let mut y = N_VNew_Serial(neq, &ctx); /* Create serial vector for solution */
    N_VConst(0.0, &mut y); /* Specify initial condition */

    /* Call LSRKStepCreateSTS to initialize the ARK timestepper module
    and specify the right-hand side function in y'=f(t,y), the initial
    time T0, and the initial dependent variable vector y. */
    let mut arkode_mem = LSRKStepCreateSTS(Some(f), t0, &y, &ctx).expect("LSRKStepCreateSTS");

    /* Set routines */
    let flag = ARKodeSetUserData(&mut arkode_mem, Some(Box::new(user_data_arr)));
    assert!(flag >= 0, "ARKodeSetUserData failed with flag = {}", flag);

    /* Specify tolerances */
    let flag = ARKodeSStolerances(&mut arkode_mem, reltol, abstol);
    assert!(flag >= 0, "ARKodeSStolerances failed with flag = {}", flag);

    /* Dominant Eigenvalue Estimator (DEE) variables */
    let max_iters: i64 = 100; /* max number of power iterations (PI) */
    let numwarmup: i32 = 10; /* number of preprocessing warmups */
    let rel_tol: f64 = 5.0e-3; /* relative error for PI */

    /* Set the initial eigenvector for the DEE (C: rand()/RAND_MAX per
       entry; normalization makes any positive value equivalent for
       NEQ = 1) */
    let mut q = N_VClone(&y);
    q.data[0] = 0.5;

    /* Create power iteration dominant eigenvalue estimator (DEE) */
    let mut dee = SUNDomEigEstimator_Power(&q, max_iters, rel_tol, &ctx);
    drop(q);

    /* Override any current settings with command-line options (C does
       this after the attach; see the header note) */
    let flag = SUNDomEigEstimator_SetOptions(&mut dee, None, None, &argv);
    assert!(flag == 0, "SUNDomEigEstimator_SetOptions failed with flag = {}", flag);

    /* Attach the DEE to the LSRKStep module. */
    let flag = LSRKStepSetDomEigEstimator(&mut arkode_mem, Some(dee));
    assert!(flag >= 0, "LSRKStepSetDomEigEstimator failed with flag = {}", flag);

    /* Set the number of preprocessing warmups for the first estimate */
    let flag = LSRKStepSetNumDomEigEstInitPreprocessIters(&mut arkode_mem, numwarmup);
    assert!(
        flag >= 0,
        "LSRKStepSetNumDomEigEstInitPreprocessIters failed with flag = {}",
        flag
    );

    /* Specify after how many successful steps dom_eig is recomputed */
    let flag = LSRKStepSetDomEigFrequency(&mut arkode_mem, 25);
    assert!(flag >= 0, "LSRKStepSetDomEigFrequency failed with flag = {}", flag);

    /* Specify max number of stages allowed */
    let flag = LSRKStepSetMaxNumStages(&mut arkode_mem, 200);
    assert!(flag >= 0, "LSRKStepSetMaxNumStages failed with flag = {}", flag);

    /* Specify max number of steps allowed */
    let flag = ARKodeSetMaxNumSteps(&mut arkode_mem, 2000);
    assert!(flag >= 0, "ARKodeSetMaxNumSteps failed with flag = {}", flag);

    /* Specify safety factor for user provided dom_eig */
    let flag = LSRKStepSetDomEigSafetyFactor(&mut arkode_mem, 1.01);
    assert!(flag >= 0, "LSRKStepSetDomEigSafetyFactor failed with flag = {}", flag);

    /* Specify the number of preprocessing warmups before each estimate
       call succeeding the very first estimate call. */
    let flag = LSRKStepSetNumDomEigEstPreprocessIters(&mut arkode_mem, 0);
    assert!(
        flag >= 0,
        "LSRKStepSetNumDomEigEstPreprocessIters failed with flag = {}",
        flag
    );

    /* Specify the Runge--Kutta--Chebyshev LSRK method by name */
    let flag = LSRKStepSetSTSMethodByName(&mut arkode_mem, "ARKODE_LSRK_RKC_2");
    assert!(flag >= 0, "LSRKStepSetSTSMethodByName failed with flag = {}", flag);

    /* Override any current settings with command-line options */
    let flag = ARKodeSetOptions(&mut arkode_mem, None, None, &argv);
    assert!(flag >= 0, "ARKodeSetOptions failed with flag = {}", flag);

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
    let flag = check_ans(&y, t, reltol, abstol);
    compute_error(&y, t);

    /* Clean up and return */
    drop(y);
    let mut slot = Some(arkode_mem);
    ARKodeFree(&mut slot);

    std::process::exit(flag);
}
