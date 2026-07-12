/* -----------------------------------------------------------------------------
 * Translated from examples/cvode/serial/cvsPendulum_dns.c (CVODE 7.7.0)
 * Programmer(s): Radu Serban and David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * This example solves a simple pendulum equation in Cartesian coordinates where
 * the pendulum bob has mass 1 and is suspended from the origin with a rod of
 * length 1. The governing equations are
 *
 * x'  = vx
 * y'  = vy
 * vx' = -x * T
 * vy' = -y * T - g
 *
 * with the constraints
 *
 * x^2 + y^2 - 1 = 0
 * x * vx + y * vy = 0
 *
 * where x and y are the pendulum bob position, vx and vy are the bob velocity
 * in the x and y directions respectively, T is the tension in the rod, and
 * g is acceleration due to gravity chosen such that the pendulum has period 2.
 * The initial condition at t = 0 is x = 1, y = 0, vx = 0, and vy = 0.
 *
 * A reference solution is computed using the pendulum equation in terms of the
 * angle between the x-axis and the pendulum rod i.e., theta in [0, -pi]. The
 * governing equations are
 *
 * theta'  = vtheta
 * vtheta' = -g * cos(theta)
 *
 * where theta is the angle from the x-axis, vtheta is the angular velocity, and
 * g the same acceleration due to gravity from above. The initial condition at
 * t = 0 is theta = 0 and vtheta = 0.
 *
 * The Cartesian formulation is run to a final time tf (default 30) with and
 * without projection for various integration tolerances. The error in the
 * position and velocity at tf compared to the reference solution, the error in
 * the position constraint equation, and various integrator statistics are
 * printed to the screen for each run.
 *
 * When projection is enabled a user-supplied function is used to project the
 * position, velocity, and error to the constraint manifold.
 *
 * Optional command line inputs may be used to change the final simulation time
 * (default 30), the initial tolerance (default 1e-5), the number of outputs
 * (default 1), or disable error projection. Use the option --help for a list
 * of the command line flags.
 * ---------------------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use std::io::Write;

use cvodes_rs::sundials_utils::fmt_e;
use cvodes_rs::*;

/* Problem Constants */
const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;
const GRAV: f64 = 13.750371636040745654980191559621114395801712;

/* -----------------------------------------------------------------------------
 * Main Program
 * ---------------------------------------------------------------------------*/

fn main() {
    let mut nout: i32 = 1; /* number of outputs       */
    let mut rtol: f64 = 1.0e-5; /* base relative tolerance */
    let mut atol: f64 = 1.0e-5; /* base absolute tolerance */
    let mut tf: f64 = 30.0; /* final integration time  */
    let mut projerr: bool = SUNTRUE; /* enable error projection */

    /* Create the SUNDIALS context */
    let sunctx = SUNContext_Create();

    /* Read command line inputs */
    let args: Vec<String> = std::env::args().collect();
    let mut retval = ReadInputs(&args, &mut rtol, &mut atol, &mut tf, &mut nout, &mut projerr);
    if check_retval(retval, "ReadInputs") {
        std::process::exit(1);
    }

    /* Compute reference solution */
    let mut yref = N_VNew_Serial(4, &sunctx);

    retval = RefSol(tf, &mut yref, nout, &sunctx);
    if check_retval(retval, "RefSol") {
        std::process::exit(1);
    }

    /* Create serial vector to store the initial condition */
    let mut yy0 = N_VNew_Serial(4, &sunctx);

    /* Set the initial condition values */
    yy0.data[0] = ONE; /* x  */
    yy0.data[1] = ZERO; /* y  */
    yy0.data[2] = ZERO; /* xd */
    yy0.data[3] = ZERO; /* yd */

    /* Create CVODE memory */
    let mut cvode_mem = CVodeCreate(CV_BDF, &sunctx);

    /* Initialize CVODE */
    retval = CVodeInit(&mut cvode_mem, f, ZERO, &yy0);
    if check_retval(retval, "CVodeInit") {
        std::process::exit(1);
    }

    /* Create dense SUNMatrix for use in linear solves */
    let A = SUNDenseMatrix(4, 4, &sunctx);

    /* Create dense SUNLinearSolver object */
    let LS = SUNLinSol_Dense(&yy0, &A, &sunctx);

    /* Attach the matrix and linear solver to CVODE */
    retval = CVodeSetLinearSolver(&mut cvode_mem, LS, Some(A));
    if check_retval(retval, "CVodeSetLinearSolver") {
        std::process::exit(1);
    }

    /* Set a user-supplied projection function */
    retval = CVodeSetProjFn(&mut cvode_mem, proj);
    if check_retval(retval, "CVodeSetProjFn") {
        std::process::exit(1);
    }

    /* Set maximum number of steps between outputs */
    retval = CVodeSetMaxNumSteps(&mut cvode_mem, 50000);
    if check_retval(retval, "CVodeSetMaxNumSteps") {
        std::process::exit(1);
    }

    /* Compute the solution with various tolerances */
    for _i in 0..5 {
        /* Output tolerance and output header for this run */
        println!("\n\nrtol = {}, atol = {}", fmt_e(rtol, 8, 2), fmt_e(atol, 8, 2));
        print!("Project    x         y");
        print!("         x'        y'     |     g      |    ");
        print!("nst     rhs eval    setups (J eval)  |   cf   ef\n");

        /* Compute solution with projection */
        retval = GetSol(&mut cvode_mem, &yy0, rtol, atol, tf, nout, SUNTRUE, projerr, &yref, &sunctx);
        if check_retval(retval, "GetSol") {
            std::process::exit(1);
        }

        /* Compute solution without projection */
        retval = GetSol(&mut cvode_mem, &yy0, rtol, atol, tf, nout, SUNFALSE, SUNFALSE, &yref, &sunctx);
        if check_retval(retval, "GetSol") {
            std::process::exit(1);
        }

        /* Reduce tolerance for next run */
        rtol /= 10.0;
        atol /= 10.0;
    }

    /* Free memory (RAII for vectors, matrix, and linear solver) */
    CVodeFree(cvode_mem);
}

/* -----------------------------------------------------------------------------
 * Functions to integrate the Cartesian and reference systems
 * ---------------------------------------------------------------------------*/

/* Compute the Cartesian system solution */
fn GetSol(
    cvode_mem: &mut CVodeMem,
    yy0: &NVector,
    rtol: f64,
    atol: f64,
    tf: f64,
    nout: i32,
    proj: bool,
    projerr: bool,
    yref: &NVector,
    sunctx: &SUNContext,
) -> i32 {
    let mut retval: i32; /* reusable return flag */

    /* Enable or disable projection */
    if proj {
        print!("  YES   ");
        retval = CVodeSetProjFrequency(cvode_mem, 1);
        if check_retval(retval, "CVodeSetProjFrequency") {
            return 1;
        }

        /* Enable or disable error projection */
        retval = CVodeSetProjErrEst(cvode_mem, projerr);
        if check_retval(retval, "CVodeSetProjErrEst") {
            return 1;
        }
    } else {
        retval = CVodeSetProjFrequency(cvode_mem, 0);
        if check_retval(retval, "CVodeSetProjFrequency") {
            return 1;
        }
        print!("  NO    ");
    }

    /* Create vector to store the solution */
    let mut yy = N_VNew_Serial(4, sunctx);

    /* Copy initial condition into solution vector */
    N_VScale(ONE, yy0, &mut yy);

    /* Reinitialize CVODE for this run */
    retval = CVodeReInit(cvode_mem, ZERO, yy0);
    if check_retval(retval, "CVodeReInit") {
        return retval;
    }

    /* Set integration tolerances for this run */
    retval = CVodeSStolerances(cvode_mem, rtol, atol);
    if check_retval(retval, "CVodeSStolerances") {
        return retval;
    }

    /* Open output file */
    let outname = if proj {
        format!(
            "cvsPendulum_dns_rtol_{}_atol_{}_proj.txt",
            fmt_e(rtol, 3, 2),
            fmt_e(atol, 3, 2)
        )
    } else {
        format!(
            "cvsPendulum_dns_rtol_{}_atol_{}.txt",
            fmt_e(rtol, 3, 2),
            fmt_e(atol, 3, 2)
        )
    };
    let mut FID = std::fs::File::create(&outname).expect("fopen");

    /* Output initial condition */
    let _ = writeln!(
        FID,
        "{} {} {} {} {}",
        fmt_e(ZERO, 24, 16),
        fmt_e(yy.data[0], 24, 16),
        fmt_e(yy.data[1], 24, 16),
        fmt_e(yy.data[2], 24, 16),
        fmt_e(yy.data[3], 24, 16)
    );

    /* Integrate to tf and peridoically output the solution */
    let dtout = tf / nout as f64; /* output frequency */
    let mut tout = dtout; /* output time      */
    let mut t = ZERO; /* return time      */

    for out in 0..nout {
        /* Set stop time (do not interpolate output) */
        retval = CVodeSetStopTime(cvode_mem, tout);
        if check_retval(retval, "CVodeSetStopTime") {
            return retval;
        }

        /* Integrate to tout */
        retval = CVode(cvode_mem, tout, &mut yy, &mut t, CV_NORMAL);
        if check_retval(retval, "CVode") {
            return retval;
        }

        /* Write output */
        let _ = writeln!(
            FID,
            "{} {} {} {} {}",
            fmt_e(t, 24, 16),
            fmt_e(yy.data[0], 24, 16),
            fmt_e(yy.data[1], 24, 16),
            fmt_e(yy.data[2], 24, 16),
            fmt_e(yy.data[3], 24, 16)
        );

        /* Update output time */
        if out < nout - 1 {
            tout += dtout;
        } else {
            tout = tf;
        }
    }

    /* Close output file */
    drop(FID);

    /* Compute the constraint violation */
    let mut x = yy.data[0]; /* position values  */
    let mut y = yy.data[1];
    let g = (x * x + y * y - ONE).abs(); /* constraint value */

    /* Compute the absolute error compared to the reference solution
       (C: N_VLinearSum(ONE, yy, -ONE, yref, yy); N_VAbs(yy, yy) — aliased) */
    yy.linear_sum_with(ONE, -ONE, yref);
    yy.abs_inplace();

    x = yy.data[0];
    y = yy.data[1];
    let xd = yy.data[2]; /* velocity values */
    let yd = yy.data[3];

    /* Output errors */
    print!(
        "{}  {}  {}  {}  |  {}  |",
        fmt_e(x, 8, 2),
        fmt_e(y, 8, 2),
        fmt_e(xd, 8, 2),
        fmt_e(yd, 8, 2),
        fmt_e(g, 8, 2)
    );

    /* Free solution vector */
    drop(yy);

    /* Integrator stats */
    let (mut nst, mut nfe, mut nsetups, mut nje, mut nfeLS, mut ncfn, mut netf) =
        (0i64, 0i64, 0i64, 0i64, 0i64, 0i64, 0i64);

    /* Get integrator stats */
    retval = CVodeGetNumSteps(cvode_mem, &mut nst);
    if check_retval(retval, "CVodeGetNumSteps") {
        return retval;
    }

    retval = CVodeGetNumRhsEvals(cvode_mem, &mut nfe);
    if check_retval(retval, "CVodeGetNumFctEvals") {
        return retval;
    }

    retval = CVodeGetNumLinSolvSetups(cvode_mem, &mut nsetups);
    if check_retval(retval, "CVodeGetNumLinSolvSetups") {
        return retval;
    }

    retval = CVodeGetNumErrTestFails(cvode_mem, &mut netf);
    if check_retval(retval, "CVodeGetNumErrTestFails") {
        return retval;
    }

    retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &mut ncfn);
    if check_retval(retval, "CVodeGetNumNonlinSolvConvFails") {
        return retval;
    }

    retval = CVodeGetNumJacEvals(cvode_mem, &mut nje);
    if check_retval(retval, "CVodeGetNumJacEvals") {
        return retval;
    }

    retval = CVodeGetNumLinRhsEvals(cvode_mem, &mut nfeLS);
    if check_retval(retval, "CVodeGetNumLinRhsEvals") {
        return retval;
    }

    /* Output stats */
    println!(
        " {:6}   {:6}+{:<4}     {:4} ({:3})     |  {:3}  {:3}",
        nst, nfe, nfeLS, nsetups, nje, ncfn, netf
    );

    0
}

/* Compute the reference system solution */
fn RefSol(tf: f64, yref: &mut NVector, nout: i32, sunctx: &SUNContext) -> i32 {
    let tol: f64 = 1.0e-14; /* integration tolerance */

    /* Create the solution vector */
    let mut yy = N_VNew_Serial(2, sunctx);

    /* Set the initial condition */
    yy.data[0] = ZERO; /* theta  */
    yy.data[1] = ZERO; /* theta' */

    /* Create CVODE memory */
    let mut cvode_mem = CVodeCreate(CV_BDF, sunctx);

    /* Initialize CVODE */
    let mut retval = CVodeInit(&mut cvode_mem, fref, ZERO, &yy);
    if check_retval(retval, "CVodeInit") {
        return 1;
    }

    /* Set integration tolerances */
    retval = CVodeSStolerances(&mut cvode_mem, tol, tol);
    if check_retval(retval, "CVodeSStolerances") {
        return 1;
    }

    /* Create dense SUNMatrix for use in linear solves */
    let A = SUNDenseMatrix(2, 2, sunctx);

    /* Create dense SUNLinearSolver object */
    let LS = SUNLinSol_Dense(&yy, &A, sunctx);

    /* Attach the matrix and linear solver to CVODE */
    retval = CVodeSetLinearSolver(&mut cvode_mem, LS, Some(A));
    if check_retval(retval, "CVodeSetLinearSolver") {
        return 1;
    }

    /* Set CVODE optional inputs */
    retval = CVodeSetMaxNumSteps(&mut cvode_mem, 100000);
    if check_retval(retval, "CVodeSetMaxNumSteps") {
        return 1;
    }

    retval = CVodeSetStopTime(&mut cvode_mem, tf);
    if check_retval(retval, "CVodeSetStopTime") {
        return 1;
    }

    /* Open output file */
    let mut FID = std::fs::File::create("cvsPendulum_dns_ref.txt").expect("fopen");

    /* Output initial condition */
    let mut th = yy.data[0]; /* theta     */
    let mut thd = yy.data[1]; /* theta dot */
    let _ = writeln!(
        FID,
        "{} {} {} {} {}",
        fmt_e(ZERO, 24, 16),
        fmt_e(th.cos(), 24, 16),
        fmt_e(th.sin(), 24, 16),
        fmt_e(-thd * th.sin(), 24, 16),
        fmt_e(thd * th.cos(), 24, 16)
    );

    /* Integrate to tf and periodically output the solution */
    let dtout = tf / nout as f64; /* output frequency */
    let mut tout = dtout; /* output time      */
    let mut t = ZERO; /* return time      */

    for out in 0..nout {
        /* Set stop time (do not interpolate output) */
        retval = CVodeSetStopTime(&mut cvode_mem, tout);
        if check_retval(retval, "CVodeSetStopTime") {
            return retval;
        }

        /* Integrate to tout */
        retval = CVode(&mut cvode_mem, tf, &mut yy, &mut t, CV_NORMAL);
        if check_retval(retval, "CVode") {
            return retval;
        }

        /* Write output */
        th = yy.data[0];
        thd = yy.data[1];
        let _ = writeln!(
            FID,
            "{} {} {} {} {}",
            fmt_e(t, 24, 16),
            fmt_e(th.cos(), 24, 16),
            fmt_e(th.sin(), 24, 16),
            fmt_e(-thd * th.sin(), 24, 16),
            fmt_e(thd * th.cos(), 24, 16)
        );

        /* Update output time */
        if out < nout - 1 {
            tout += dtout;
        } else {
            tout = tf;
        }
    }

    /* Close output file */
    drop(FID);

    /* Get solution components */
    th = yy.data[0];
    thd = yy.data[1];

    /* Convert to Cartesian reference solution */
    yref.data[0] = th.cos();
    yref.data[1] = th.sin();
    yref.data[2] = -thd * th.sin();
    yref.data[3] = thd * th.cos();

    /* Free memory (RAII for yy, A, LS) */
    CVodeFree(cvode_mem);

    0
}

/* -----------------------------------------------------------------------------
 * Functions provided to CVODE
 * ---------------------------------------------------------------------------*/

/* ODE RHS function for the reference system */
fn fref(_t: f64, yy: &NVector, fy: &mut NVector, _f_data: &mut UserData) -> i32 {
    fy.data[0] = yy.data[1]; /* theta'          */
    fy.data[1] = -GRAV * (yy.data[0]).cos(); /* -g * cos(theta) */
    0
}

/* ODE RHS function for the Cartesian system */
fn f(_t: f64, yy: &NVector, fy: &mut NVector, _f_data: &mut UserData) -> i32 {
    /* Get vector components */
    let x = yy.data[0]; /* positions  */
    let y = yy.data[1];
    let xd = yy.data[2]; /* velocities */
    let yd = yy.data[3];

    /* Compute tension */
    let tmp = xd * xd + yd * yd - GRAV * y;

    /* Compute RHS */
    fy.data[0] = xd;
    fy.data[1] = yd;
    fy.data[2] = -x * tmp;
    fy.data[3] = -y * tmp - GRAV;

    0
}

/* Projection function */
fn proj(
    _t: f64,
    yy: &NVector,
    corr: &mut NVector,
    _epsProj: f64,
    err: Option<&mut NVector>,
    _pdata: &mut UserData,
) -> i32 {
    /* Extract current solution */

    let x = yy.data[0]; /* positions  */
    let y = yy.data[1];
    let xd = yy.data[2]; /* velocities */
    let yd = yy.data[3];

    /* Project positions */

    let R = (x * x + y * y).sqrt();

    let x_new = x / R;
    let y_new = y / R;

    /* Project velocities
     *
     *        +-            -+  +-    -+
     *        |  y*y    -x*y |  |  xd  |
     *  P v = |              |  |      |
     *        | -x*y     x*x |  |  yd  |
     *        +-            -+  +-    -+
     */

    let xd_new = xd * y_new * y_new - yd * x_new * y_new;
    let yd_new = -xd * x_new * y_new + yd * x_new * x_new;

    /* Return position and velocity corrections */

    corr.data[0] = x_new - x;
    corr.data[1] = y_new - y;
    corr.data[2] = xd_new - xd;
    corr.data[3] = yd_new - yd;

    /* Project error P * err */
    if let Some(err) = err {
        let e1 = err.data[0];
        let e2 = err.data[1];
        let e3 = err.data[2];
        let e4 = err.data[3];

        let e1_new = y_new * y_new * e1 - x_new * y_new * e2;
        let e2_new = -x_new * y_new * e1 + x_new * x_new * e2;

        let e3_new = y_new * y_new * e3 - x_new * y_new * e4;
        let e4_new = -x_new * y_new * e3 + x_new * x_new * e4;

        err.data[0] = e1_new;
        err.data[1] = e2_new;
        err.data[2] = e3_new;
        err.data[3] = e4_new;
    }

    0
}

/* -----------------------------------------------------------------------------
 * Private helper functions
 * ---------------------------------------------------------------------------*/

/* C atof: skip leading whitespace, convert the longest valid leading numeric
   prefix, return 0.0 if no conversion can be performed */
fn atof(s: &str) -> f64 {
    let t = s.trim_start();
    let b = t.as_bytes();
    let mut end = 0;
    if end < b.len() && (b[end] == b'+' || b[end] == b'-') {
        end += 1;
    }
    let mut seen_digit = false;
    while end < b.len() && b[end].is_ascii_digit() {
        end += 1;
        seen_digit = true;
    }
    if end < b.len() && b[end] == b'.' {
        end += 1;
        while end < b.len() && b[end].is_ascii_digit() {
            end += 1;
            seen_digit = true;
        }
    }
    if !seen_digit {
        return 0.0;
    }
    if end < b.len() && (b[end] == b'e' || b[end] == b'E') {
        let mut e = end + 1;
        if e < b.len() && (b[e] == b'+' || b[e] == b'-') {
            e += 1;
        }
        let mut exp_digit = false;
        while e < b.len() && b[e].is_ascii_digit() {
            e += 1;
            exp_digit = true;
        }
        if exp_digit {
            end = e;
        }
    }
    t[..end].parse().unwrap_or(0.0)
}

/* C atoi: skip leading whitespace, convert the longest valid leading integer
   prefix, return 0 if no conversion can be performed */
fn atoi(s: &str) -> i32 {
    let t = s.trim_start();
    let b = t.as_bytes();
    let mut end = 0;
    if end < b.len() && (b[end] == b'+' || b[end] == b'-') {
        end += 1;
    }
    let mut seen_digit = false;
    while end < b.len() && b[end].is_ascii_digit() {
        end += 1;
        seen_digit = true;
    }
    if !seen_digit {
        return 0;
    }
    t[..end].parse().unwrap_or(0)
}

/* Read command line unputs */
fn ReadInputs(
    args: &[String],
    rtol: &mut f64,
    atol: &mut f64,
    tf: &mut f64,
    nout: &mut i32,
    projerr: &mut bool,
) -> i32 {
    let mut arg_idx = 1;

    /* check for input args */
    while arg_idx < args.len() {
        if args[arg_idx] == "--tol" {
            arg_idx += 1;
            *rtol = atof(&args[arg_idx]);
            arg_idx += 1;
            *atol = atof(&args[arg_idx]);
            arg_idx += 1;
        } else if args[arg_idx] == "--tf" {
            arg_idx += 1;
            *tf = atof(&args[arg_idx]);
            arg_idx += 1;
        } else if args[arg_idx] == "--nout" {
            arg_idx += 1;
            *nout = atoi(&args[arg_idx]);
            arg_idx += 1;
        } else if args[arg_idx] == "--noerrproj" {
            arg_idx += 1;
            *projerr = SUNFALSE;
        } else if args[arg_idx] == "--help" {
            InputHelp();
            return -1;
        } else {
            eprint!("ERROR: Invalid input {}", args[arg_idx]);
            InputHelp();
            return -1;
        }
    }

    0
}

/* Print command line options */
fn InputHelp() {
    println!("\nCommand line options:");
    println!("  --tol <rtol> <atol> : relative and absolute tolerance");
    println!("  --tf <time>         : final simulation time");
    println!("  --nout <outputs>    : number of outputs");
    println!("  --noerrproj         : disable error projection");
}

/* Check function return value (C opt 1: retval < 0 is an error; the C opt 0
   NULL-pointer checks have no Rust counterpart) */
fn check_retval(retval: i32, funcname: &str) -> bool {
    if retval < 0 {
        eprint!("\nERROR: {}() returned = {}\n\n", funcname, retval);
        return true;
    }
    false
}
