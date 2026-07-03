/* -----------------------------------------------------------------------------
 * Translated from examples/cvode/serial/cvParticle_dns.c (CVODE 7.7.0)
 * Programmer(s): David J. Gardner @ LLNL
 * Based on an example from Jean-Luc Fattebert @ ORNL
 * -----------------------------------------------------------------------------
 * This example solves the equation for a particle moving conterclockwise with
 * velocity alpha on the unit circle in the xy-plane. The ODE system is given by
 *
 *   x' = -alpha * y
 *   y' =  alpha * x
 *
 * where x and y are subject to the constraint
 *
 *   x^2 + y^2 - 1 = 0
 *
 * with initial condition x = 1 and y = 0 at t = 0. The system has the analytic
 * solution
 *
 *  x(t) = cos(alpha * t)
 *  y(t) = sin(alpha * t)
 *
 * For a description of the command line options for this example run the
 * program with the --help flag.
 * ---------------------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use std::io::Write;

use cvode_rs::sundials_utils::{fmt_e, fmt_g};
use cvode_rs::*;

/* Problem Constants */
const PI: f64 = 3.141592653589793238462643383279502884197169;
const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;

/* User-defined data structure */
#[derive(Clone)]
struct UserDataStruct {
    alpha: f64, /* particle velocity */

    orbits: i32, /* number of orbits */
    torbit: f64, /* orbit time       */

    rtol: f64, /* integration tolerances */
    atol: f64,

    proj: i32,    /* enable/disable solution projection */
    projerr: i32, /* enable/disable error projection */

    tstop: i32, /* use tstop mode */
    nout: i32,  /* number of outputs per orbit */
}

/* C atof/atoi (numeric args as used by this example; non-numeric input -> 0) */
fn atof(s: &str) -> f64 {
    s.trim_start().parse().unwrap_or(0.0)
}

fn atoi(s: &str) -> i32 {
    s.trim_start().parse().unwrap_or(0)
}

/* -----------------------------------------------------------------------------
 * Functions provided to CVODE
 * ---------------------------------------------------------------------------*/

/* Compute the right-hand side function, y' = f(t,y) */
fn f(_t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let udata = user_data
        .as_mut()
        .unwrap()
        .downcast_mut::<UserDataStruct>()
        .unwrap();
    let ydata = &y.data;
    let fdata = &mut ydot.data;

    fdata[0] = -(udata.alpha) * ydata[1];
    fdata[1] = (udata.alpha) * ydata[0];

    0
}

/* Compute the Jacobian of the right-hand side function, J(t,y) = df/dy */
fn Jac(
    _t: f64,
    _y: &NVector,
    _fy: &NVector,
    j: &mut SUNMatrix,
    user_data: &mut UserData,
    _tmp1: &mut NVector,
    _tmp2: &mut NVector,
    _tmp3: &mut NVector,
) -> i32 {
    let udata = user_data
        .as_mut()
        .unwrap()
        .downcast_mut::<UserDataStruct>()
        .unwrap();
    let jm = match j {
        SUNMatrix::Dense(m) => m,
        _ => return -1,
    };

    /* C writes the column-major Jdata array: Jdata[k] -> (row k%2, col k/2) */
    jm.set(0, 0, ZERO);
    jm.set(1, 0, -(udata.alpha));
    jm.set(0, 1, udata.alpha);
    jm.set(1, 1, ZERO);

    0
}

/* Project the solution onto the constraint manifold */
fn Proj(
    _t: f64,
    ycur: &NVector,
    corr: &mut NVector,
    _epsProj: f64,
    err: Option<&mut NVector>,
    _user_data: &mut UserData,
) -> i32 {
    let ydata = &ycur.data;
    let cdata = &mut corr.data;
    let x = ydata[0];
    let y = ydata[1];

    /* project onto the unit circle */
    let r = (x * x + y * y).sqrt();

    let xp = x / r;
    let yp = y / r;

    /* correction to the unprojected solution */
    cdata[0] = xp - x;
    cdata[1] = yp - y;

    /* project the error */
    if let Some(e) = err {
        let edata = &mut e.data;

        let errxp = edata[0] * yp * yp - edata[1] * xp * yp;
        let erryp = -edata[0] * xp * yp + edata[1] * xp * xp;

        edata[0] = errxp;
        edata[1] = erryp;
    }

    0
}

/* -----------------------------------------------------------------------------
 * Private helper functions
 * ---------------------------------------------------------------------------*/

fn InitUserData(argv: &[String], udata: &mut UserDataStruct) -> i32 {
    let mut arg_idx = 1;

    /* set default values */
    udata.alpha = ONE;

    udata.orbits = 100;
    udata.torbit = (TWO * PI) / udata.alpha;

    udata.rtol = 1.0e-4;
    udata.atol = 1.0e-9;

    udata.proj = 1;
    udata.projerr = 0;

    udata.tstop = 0;
    udata.nout = 0;

    /* check for input args */
    while arg_idx < argv.len() {
        if argv[arg_idx] == "--alpha" {
            arg_idx += 1;
            udata.alpha = atof(&argv[arg_idx]);
            arg_idx += 1;
            udata.torbit = (TWO * PI) / udata.alpha;
        } else if argv[arg_idx] == "--orbits" {
            arg_idx += 1;
            udata.orbits = atoi(&argv[arg_idx]);
            arg_idx += 1;
        } else if argv[arg_idx] == "--rtol" {
            arg_idx += 1;
            udata.rtol = atof(&argv[arg_idx]);
            arg_idx += 1;
        } else if argv[arg_idx] == "--atol" {
            arg_idx += 1;
            udata.atol = atof(&argv[arg_idx]);
            arg_idx += 1;
        } else if argv[arg_idx] == "--proj" {
            arg_idx += 1;
            udata.proj = atoi(&argv[arg_idx]);
            arg_idx += 1;
        } else if argv[arg_idx] == "--projerr" {
            arg_idx += 1;
            udata.projerr = atoi(&argv[arg_idx]);
            arg_idx += 1;
        } else if argv[arg_idx] == "--nout" {
            arg_idx += 1;
            udata.nout = atoi(&argv[arg_idx]);
            arg_idx += 1;
        } else if argv[arg_idx] == "--tstop" {
            arg_idx += 1;
            udata.tstop = 1;
        } else if argv[arg_idx] == "--help" {
            InputHelp();
            return -1;
        } else {
            eprint!("ERROR: Invalid input {}", argv[arg_idx]);
            InputHelp();
            return -1;
        }
    }

    /* If projection is disabled then disable error projection */
    if udata.proj == 0 {
        udata.projerr = 0;
    }

    0
}

fn PrintUserData(udata: &UserDataStruct) -> i32 {
    println!("\nParticle traveling on the unit circle example");
    println!("---------------------------------------------");
    println!("alpha      = {}", fmt_e(udata.alpha, 0, 4));
    println!("num orbits = {}", udata.orbits);
    println!("---------------------------------------------");
    println!("rtol       = {}", fmt_g(udata.rtol, 0, 6));
    println!("atol       = {}", fmt_g(udata.atol, 0, 6));
    println!("proj sol   = {}", udata.proj);
    println!("proj err   = {}", udata.projerr);
    println!("nout       = {}", udata.nout);
    println!("tstop      = {}", udata.tstop);
    println!("---------------------------------------------");

    0
}

/* Print command line options */
fn InputHelp() {
    println!("\nCommand line options:");
    println!("  --alpha <vel>      : particle velocity");
    println!("  --orbits <orbits>  : number of orbits to perform");
    println!("  --rtol <rtol>      : relative tolerance");
    println!("  --atol <atol>      : absolute tolerance");
    println!("  --proj <1 or 0>    : enable (1) / disable (0) projection");
    println!("  --projerr <1 or 0> : enable (1) / disable (0) error projection");
    println!("  --nout <nout>      : outputs per period");
    println!("  --tstop            : stop at output time (do not interpolate)");
}

/* Compute the analytical solution */
fn ComputeSolution(t: f64, y: &mut NVector, udata: &UserDataStruct) -> i32 {
    let ydata = &mut y.data;

    ydata[0] = ((udata.alpha) * t).cos();
    ydata[1] = ((udata.alpha) * t).sin();

    0
}

/* Compute the error in the solution and constraint */
fn ComputeError(t: f64, y: &NVector, e: &mut NVector, ec: &mut f64, udata: &UserDataStruct) -> i32 {
    let ydata = &y.data;

    /* solution error */
    let retval = ComputeSolution(t, e, udata);
    if check_retval(retval, "ComputeSolution") {
        return 1;
    }
    /* N_VLinearSum(ONE, y, -ONE, e, e): e is aliased -> in-place method
       (e = -(e - y) is bitwise equal to y - e) */
    e.linear_sum_with(-ONE, ONE, y);

    /* constraint error */
    *ec = ydata[0] * ydata[0] + ydata[1] * ydata[1] - ONE;

    0
}

/* Output the solution to the screen or disk */
fn WriteOutput(
    t: f64,
    y: &NVector,
    e: &NVector,
    ec: f64,
    screenfile: i32,
    YFID: Option<&mut std::fs::File>,
    EFID: Option<&mut std::fs::File>,
) -> i32 {
    let ydata = &y.data;
    let edata = &e.data;

    if screenfile == 0 {
        /* output solution and error to screen */
        println!(
            "{} {} {} {} {} {}",
            fmt_e(t, 0, 4),
            fmt_e(ydata[0], 14, 6),
            fmt_e(ydata[1], 14, 6),
            fmt_e(edata[0], 14, 6),
            fmt_e(edata[1], 14, 6),
            fmt_e(ec, 14, 6)
        );
    } else {
        /* check file pointers */
        let (yfid, efid) = match (YFID, EFID) {
            (Some(yfid), Some(efid)) => (yfid, efid),
            _ => return 1,
        };

        /* output solution to disk */
        let _ = writeln!(
            yfid,
            "{} {} {}",
            fmt_e(t, 24, 16),
            fmt_e(ydata[0], 24, 16),
            fmt_e(ydata[1], 24, 16)
        );

        /* output error to disk */
        let _ = writeln!(
            efid,
            "{} {} {} {}",
            fmt_e(t, 24, 16),
            fmt_e(edata[0], 24, 16),
            fmt_e(edata[1], 24, 16),
            fmt_e(ec, 24, 16)
        );
    }

    0
}

/* Print final statistics */
fn PrintStats(cvode_mem: &mut CVodeMem) -> i32 {
    let (mut nst, mut nfe, mut nsetups, mut nje) = (0i64, 0i64, 0i64, 0i64);
    let (mut nni, mut ncfn, mut netf) = (0i64, 0i64, 0i64);

    let mut retval = CVodeGetNumSteps(cvode_mem, &mut nst);
    check_retval(retval, "CVodeGetNumSteps");
    retval = CVodeGetNumRhsEvals(cvode_mem, &mut nfe);
    check_retval(retval, "CVodeGetNumRhsEvals");
    retval = CVodeGetNumLinSolvSetups(cvode_mem, &mut nsetups);
    check_retval(retval, "CVodeGetNumLinSolvSetups");
    retval = CVodeGetNumErrTestFails(cvode_mem, &mut netf);
    check_retval(retval, "CVodeGetNumErrTestFails");
    retval = CVodeGetNumNonlinSolvIters(cvode_mem, &mut nni);
    check_retval(retval, "CVodeGetNumNonlinSolvIters");
    retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &mut ncfn);
    check_retval(retval, "CVodeGetNumNonlinSolvConvFails");

    retval = CVodeGetNumJacEvals(cvode_mem, &mut nje);
    check_retval(retval, "CVodeGetNumJacEvals");

    println!("\nIntegration Statistics:");

    println!("Number of steps taken = {:<6}", nst);
    println!("Number of function evaluations = {:<6}", nfe);

    println!("Number of linear solver setups = {:<6}", nsetups);
    println!("Number of Jacobian evaluations = {:<6}", nje);

    println!("Number of nonlinear solver iterations = {:<6}", nni);
    println!("Number of convergence failures = {:<6}", ncfn);
    println!("Number of error test failures = {:<6}", netf);

    0
}

/* Check function return value (C opt 1; opt-0 NULL checks vanish in Rust) */
fn check_retval(retval: i32, funcname: &str) -> bool {
    if retval < 0 {
        eprintln!("\nERROR: {}() returned = {}\n", funcname, retval);
        return true;
    }
    false
}

/* -----------------------------------------------------------------------------
 * Main Program
 * ---------------------------------------------------------------------------*/

fn main() {
    let mut t = ZERO; /* current integration time   */
    let dtout; /* output spacing             */
    let mut tout; /* next output time           */
    let mut ec = ZERO; /* constraint error           */
    let totalout; /* output counter             */

    /* Create the SUNDIALS context */
    let sunctx = SUNContext_Create();

    /* Allocate and initialize user data structure */
    let mut udata = UserDataStruct {
        alpha: ZERO,
        orbits: 0,
        torbit: ZERO,
        rtol: ZERO,
        atol: ZERO,
        proj: 0,
        projerr: 0,
        tstop: 0,
        nout: 0,
    };

    let argv: Vec<String> = std::env::args().collect();
    let mut retval = InitUserData(&argv, &mut udata);
    if check_retval(retval, "InitUserData") {
        std::process::exit(1);
    }

    /* Create serial vector to store the solution */
    let mut y = N_VNew_Serial(2, &sunctx);

    /* Set initial contion */
    y.data[0] = ONE;
    y.data[1] = ZERO;

    /* Create serial vector to store the solution error */
    let mut e = N_VClone(&y);

    /* Set initial error */
    N_VConst(ZERO, &mut e);

    /* Create CVODE memory */
    let mut cvode_mem = CVodeCreate(CV_BDF, &sunctx);

    /* Initialize CVODE */
    retval = CVodeInit(&mut cvode_mem, f, t, &y);
    if check_retval(retval, "CVodeInit") {
        std::process::exit(1);
    }

    /* Attach user-defined data structure to CVODE (the C code shares one
       struct; it is read-only after InitUserData, so a clone is identical) */
    retval = CVodeSetUserData(&mut cvode_mem, Some(Box::new(udata.clone())));
    if check_retval(retval, "CVodeSetUserData") {
        std::process::exit(1);
    }

    /* Set integration tolerances */
    retval = CVodeSStolerances(&mut cvode_mem, udata.rtol, udata.atol);
    if check_retval(retval, "CVodeSStolerances") {
        std::process::exit(1);
    }

    /* Create dense SUNMatrix for use in linear solves */
    let a = SUNDenseMatrix(2, 2, &sunctx);

    /* Create dense SUNLinearSolver object */
    let ls = SUNLinSol_Dense(&y, &a, &sunctx);

    /* Attach the matrix and linear solver to CVODE */
    retval = CVodeSetLinearSolver(&mut cvode_mem, ls, Some(a));
    if check_retval(retval, "CVodeSetLinearSolver") {
        std::process::exit(1);
    }

    /* Set a user-supplied Jacobian function */
    retval = CVodeSetJacFn(&mut cvode_mem, Some(Jac));
    if check_retval(retval, "CVodeSetJacFn") {
        std::process::exit(1);
    }

    /* Set a user-supplied projection function */
    if udata.proj != 0 {
        retval = CVodeSetProjFn(&mut cvode_mem, Proj);
        if check_retval(retval, "CVodeSetProjFn") {
            std::process::exit(1);
        }

        retval = CVodeSetProjErrEst(&mut cvode_mem, udata.projerr != 0);
        if check_retval(retval, "CVodeSetProjErrEst") {
            std::process::exit(1);
        }
    }

    /* Set max steps between outputs */
    retval = CVodeSetMaxNumSteps(&mut cvode_mem, 100000);
    if check_retval(retval, "CVodeSetMaxNumSteps") {
        std::process::exit(1);
    }

    /* Output problem setup */
    retval = PrintUserData(&udata);
    if check_retval(retval, "PrintUserData") {
        std::process::exit(1);
    }

    /* Output initial condition */
    println!("\n     t            x              y             err x          err y       err constr");
    WriteOutput(t, &y, &e, ec, 0, None, None);

    let mut YFID: Option<std::fs::File> = None; /* solution output file */
    let mut EFID: Option<std::fs::File> = None; /* error output file    */
    if udata.nout > 0 {
        YFID = std::fs::File::create("cvParticle_solution.txt").ok();
        EFID = std::fs::File::create("cvParticle_error.txt").ok();
        WriteOutput(t, &y, &e, ec, 1, YFID.as_mut(), EFID.as_mut());
    }

    /* Integrate in time and periodically output the solution and error */
    if udata.nout > 0 {
        totalout = udata.orbits * udata.nout;
        dtout = udata.torbit / (udata.nout as f64);
    } else {
        totalout = 1;
        dtout = udata.torbit * (udata.orbits as f64);
    }
    tout = dtout;

    for out in 0..totalout {
        /* Stop at output time (do not interpolate output) */
        if udata.tstop != 0 || udata.nout == 0 {
            retval = CVodeSetStopTime(&mut cvode_mem, tout);
            if check_retval(retval, "CVodeSetStopTime") {
                std::process::exit(1);
            }
        }

        /* Advance in time */
        retval = CVode(&mut cvode_mem, tout, &mut y, &mut t, CV_NORMAL);
        if check_retval(retval, "CVode") {
            break;
        }

        /* Output solution and error */
        if udata.nout > 0 {
            retval = ComputeError(t, &y, &mut e, &mut ec, &udata);
            if check_retval(retval, "ComputeError") {
                break;
            }

            WriteOutput(t, &y, &e, ec, 1, YFID.as_mut(), EFID.as_mut());
            if check_retval(retval, "WriteOutput") {
                break;
            }
        }

        /* Update output time */
        if out < totalout - 1 {
            tout += dtout;
        } else {
            tout = udata.torbit * (udata.orbits as f64);
        }
    }

    /* Close output files */
    if udata.nout > 0 {
        drop(YFID);
        drop(EFID);
    }

    /* Output final solution and error to screen
       (the C code checks the stale `retval` here, not the calls' results) */
    ComputeError(t, &y, &mut e, &mut ec, &udata);
    if check_retval(retval, "ComputeError") {
        std::process::exit(1);
    }

    WriteOutput(t, &y, &e, ec, 0, None, None);
    if check_retval(retval, "WriteOutput") {
        std::process::exit(1);
    }

    /* Print some final statistics */
    PrintStats(&mut cvode_mem);

    /* Free memory (RAII) */
    CVodeFree(cvode_mem);
}
