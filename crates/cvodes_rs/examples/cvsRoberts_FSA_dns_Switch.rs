/* -----------------------------------------------------------------
 * Translation of examples/cvodes/serial/cvsRoberts_FSA_dns_Switch.c
 * (SUNDIALS 7.7.0).
 *
 * Modification of the cvsRoberts_FSA_dns to illustrate switching
 * on and off sensitivity computations.
 *
 * Example problem (chemical kinetics):
 *    dy1/dt = -p1*y1 + p2*y2*y3
 *    dy2/dt =  p1*y1 - p2*y2*y3 - p3*(y2)^2
 *    dy3/dt =  p3*(y2)^2
 *
 * The problem is solved five times: with sensitivities enabled
 * (user fS), disabled via CVodeSensToggleOff, re-enabled with the
 * internal DQ sensitivity RHS, re-enabled staggered with user fS
 * and partial error control, and finally with the sensitivity
 * memory freed.
 *
 * Translation notes: C stores p in the user data and hands the SAME
 * array pointer to CVodeSetSensParams, so later parameter changes
 * alias into the solver's cv_p; here the user data is the pinned
 * FSAUserData wrapper (required by the internal-DQ run) and the
 * parameter changes are mirrored into both the wrapper and cv_p.
 * The C {sensi, errconS, fsDQ, meth} bookkeeping fields live in a
 * main-side struct used only by the print helpers.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]

use cvodes_rs::cvodes::{
    CVode, CVodeCreate, CVodeFree, CVodeInit, CVodeReInit, CVodeSVtolerances,
    CVodeSensEEtolerances, CVodeSensFree, CVodeSensInit1, CVodeSensToggleOff,
};
use cvodes_rs::cvodes_io::{
    CVodeGetNumErrTestFails, CVodeGetNumNonlinSolvConvFails, CVodeGetNumNonlinSolvIters,
    CVodeGetNumRhsEvals, CVodeGetNumRhsEvalsSens, CVodeGetNumLinSolvSetups, CVodeGetNumSteps,
    CVodeGetSensNumErrTestFails, CVodeGetSensNumLinSolvSetups,
    CVodeGetSensNumNonlinSolvConvFails, CVodeGetSensNumNonlinSolvIters,
    CVodeGetSensNumRhsEvals, CVodeSetMaxNumSteps, CVodeSetSensErrCon, CVodeSetSensParams,
    CVodeSetUserData,
};
use cvodes_rs::cvodes_ls::{
    CVodeGetNumJacEvals, CVodeGetNumLinRhsEvals, CVodeSetJacFn, CVodeSetLinearSolver,
};
use cvodes_rs::sundials_utils::fmt_e;
use cvodes_rs::*;

/* Problem Constants */
const MXSTEPS: i64 = 2000; /* max number of steps */
const NEQ: i64 = 3; /* number of equations */
const T0: f64 = 0.0; /* initial time        */
const T1: f64 = 4.0e10; /* first output time   */

const ZERO: f64 = 0.0;

/* main-side bookkeeping (the C UserData control fields) */
struct RunConfig {
    sensi: bool,   /* turn on (T) or off (F) sensitivity analysis    */
    errconS: bool, /* full (T) or partial error control (F)          */
    fsDQ: bool,    /* internal DQ r.h.s sensitivity analysis (T/F)   */
    meth: i32,     /* sensitivity method                             */
}

fn get_p(user_data: &UserData) -> [f64; 3] {
    let fsa = user_data.as_ref().unwrap().downcast_ref::<FSAUserData>().unwrap();
    [fsa.p[0], fsa.p[1], fsa.p[2]]
}

/*
 * f routine. Compute f(t,y).
 */
fn f(_t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let [p1, p2, p3] = get_p(user_data);
    let y1 = y.data[0];
    let y2 = y.data[1];
    let y3 = y.data[2];

    let yd1 = -p1 * y1 + p2 * y2 * y3;
    ydot.data[0] = yd1;
    let yd3 = p3 * y2 * y2;
    ydot.data[2] = yd3;
    ydot.data[1] = -yd1 - yd3;

    0
}

/*
 * Jacobian routine. Compute J(t,y).
 */
#[allow(clippy::too_many_arguments)]
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
    let [p1, p2, p3] = get_p(user_data);
    let y2 = y.data[1];
    let y3 = y.data[2];

    if let SUNMatrix::Dense(dm) = j {
        let m = 3usize;
        dm.data[0] = -p1;
        dm.data[m] = p2 * y3;
        dm.data[2 * m] = p2 * y2;
        dm.data[1] = p1;
        dm.data[1 + m] = -p2 * y3 - 2.0 * p3 * y2;
        dm.data[1 + 2 * m] = -p2 * y2;
        dm.data[2 + m] = 2.0 * p3 * y2;
    }

    0
}

/*
 * fS routine. Compute sensitivity r.h.s.
 */
#[allow(clippy::too_many_arguments)]
fn fS(
    _Ns: i32,
    _t: f64,
    y: &NVector,
    _ydot: &NVector,
    iS: i32,
    yS: &NVector,
    ySdot: &mut NVector,
    user_data: &mut UserData,
    _tmp1: &mut NVector,
    _tmp2: &mut NVector,
) -> i32 {
    let [p1, p2, p3] = get_p(user_data);
    let y1 = y.data[0];
    let y2 = y.data[1];
    let y3 = y.data[2];
    let s1 = yS.data[0];
    let s2 = yS.data[1];
    let s3 = yS.data[2];

    let mut sd1 = -p1 * s1 + p2 * y3 * s2 + p2 * y2 * s3;
    let mut sd3 = 2.0 * p3 * y2 * s2;
    let mut sd2 = -sd1 - sd3;

    match iS {
        0 => {
            sd1 += -y1;
            sd2 += y1;
        }
        1 => {
            sd1 += y2 * y3;
            sd2 += -y2 * y3;
        }
        2 => {
            sd2 += -y2 * y2;
            sd3 += y2 * y2;
        }
        _ => {}
    }

    ySdot.data[0] = sd1;
    ySdot.data[1] = sd2;
    ySdot.data[2] = sd3;

    0
}

/* set the problem parameters: C's data->p aliases the solver's cv_p
   pointer, so both Rust owned copies are updated together */
fn set_params(cvode_mem: &mut CVodeMem, vals: [f64; 3]) {
    if let Some(d) = cvode_mem.cv_user_data.as_mut() {
        if let Some(fsa) = d.downcast_mut::<FSAUserData>() {
            fsa.p = vals.to_vec();
        }
    }
    if !cvode_mem.cv_p.is_empty() {
        cvode_mem.cv_p = vals.to_vec();
    }
}

/*
 * Runs integrator and prints final statistics when complete.
 */
fn run_cvode(cvode_mem: &mut CVodeMem, y: &mut NVector, config: &RunConfig) -> i32 {
    /* Print header for current run */
    print_header(cvode_mem, config);

    /* Call CVode in CV_NORMAL mode */
    let mut t = 0.0;
    let retval = CVode(cvode_mem, T1, y, &mut t, CV_NORMAL);
    if retval != 0 {
        return retval;
    }

    /* Print final statistics */
    let retval = print_final_stats(cvode_mem, config);
    println!();

    retval
}

fn print_header(cvode_mem: &mut CVodeMem, config: &RunConfig) {
    /* Print sensitivity control retvals */
    print!("Sensitivity: ");
    if config.sensi {
        print!("YES (");
        match config.meth {
            CV_SIMULTANEOUS => print!("SIMULTANEOUS + "),
            CV_STAGGERED => print!("STAGGERED + "),
            CV_STAGGERED1 => print!("STAGGERED-1 + "),
            _ => {}
        }
        if config.errconS {
            print!("FULL ERROR CONTROL + ");
        } else {
            print!("PARTIAL ERROR CONTROL + ");
        }
        if config.fsDQ {
            println!("DQ sensitivity RHS)");
        } else {
            println!("user-provided sensitivity RHS)");
        }
    } else {
        println!("NO");
    }

    /* Print current problem parameters */
    let p = get_p(&cvode_mem.cv_user_data);
    println!(
        "Parameters: [{}  {}  {}]",
        fmt_e(p[0], 8, 4),
        fmt_e(p[1], 8, 4),
        fmt_e(p[2], 8, 4)
    );
}

/*
 * Print some final statistics from the CVODES memory.
 */
fn print_final_stats(cvode_mem: &mut CVodeMem, config: &RunConfig) -> i32 {
    let mut nst: i64 = 0;
    let mut nfe: i64 = 0;
    let mut nsetups: i64 = 0;
    let mut nni: i64 = 0;
    let mut ncfn: i64 = 0;
    let mut netf: i64 = 0;
    let mut nfSe: i64 = 0;
    let mut nfeS: i64 = 0;
    let mut nsetupsS: i64 = 0;
    let mut nniS: i64 = 0;
    let mut ncfnS: i64 = 0;
    let mut netfS: i64 = 0;
    let mut njeD: i64 = 0;
    let mut nfeD: i64 = 0;

    CVodeGetNumSteps(cvode_mem, &mut nst);
    CVodeGetNumRhsEvals(cvode_mem, &mut nfe);
    CVodeGetNumLinSolvSetups(cvode_mem, &mut nsetups);
    CVodeGetNumErrTestFails(cvode_mem, &mut netf);
    CVodeGetNumNonlinSolvIters(cvode_mem, &mut nni);
    CVodeGetNumNonlinSolvConvFails(cvode_mem, &mut ncfn);

    if config.sensi {
        CVodeGetSensNumRhsEvals(cvode_mem, &mut nfSe);
        CVodeGetNumRhsEvalsSens(cvode_mem, &mut nfeS);
        CVodeGetSensNumLinSolvSetups(cvode_mem, &mut nsetupsS);
        if config.errconS {
            CVodeGetSensNumErrTestFails(cvode_mem, &mut netfS);
        } else {
            netfS = 0;
        }
        if config.meth == CV_STAGGERED {
            CVodeGetSensNumNonlinSolvIters(cvode_mem, &mut nniS);
            CVodeGetSensNumNonlinSolvConvFails(cvode_mem, &mut ncfnS);
        } else {
            nniS = 0;
            ncfnS = 0;
        }
    }

    let retval = CVodeGetNumJacEvals(cvode_mem, &mut njeD);
    let _ = retval;
    let retval = CVodeGetNumLinRhsEvals(cvode_mem, &mut nfeD);

    println!("Run statistics:");

    println!("   nst     = {:5}", nst);
    println!("   nfe     = {:5}", nfe);
    println!("   netf    = {:5}    nsetups  = {:5}", netf, nsetups);
    println!("   nni     = {:5}    ncfn     = {:5}", nni, ncfn);

    println!("   njeD    = {:5}    nfeD     = {:5}", njeD, nfeD);

    if config.sensi {
        println!("   -----------------------------------");
        println!("   nfSe    = {:5}    nfeS     = {:5}", nfSe, nfeS);
        println!("   netfs   = {:5}    nsetupsS = {:5}", netfS, nsetupsS);
        println!("   nniS    = {:5}    ncfnS    = {:5}", nniS, ncfnS);
    }

    retval
}

/* Check if a SUNDIALS function returned a negative value */
fn check_retval(retval: i32, funcname: &str) -> bool {
    if retval < 0 {
        eprintln!("\nSUNDIALS_ERROR: {}() failed with retval = {}\n", funcname, retval);
        return true;
    }
    false
}

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */
fn main() {
    /* Create the SUNDIALS context object for this simulation */
    let sunctx = SUNContext_Create();

    /* Initialize sensitivity variables (reaction rates for this problem) */
    let p0 = [0.04, 1.0e4, 3.0e7];

    /* Allocate initial condition vector and set context */
    let mut y0 = N_VNew_Serial(NEQ, &sunctx);

    /* Create solution and absolute tolerance vectors */
    let mut y = N_VClone(&y0);
    let mut abstol = N_VClone(&y0);

    /* Set initial conditions */
    y0.data[0] = 1.0;
    y0.data[1] = 0.0;
    y0.data[2] = 0.0;

    /* Set integration tolerances */
    let reltol = 1e-6;
    abstol.data[0] = 1e-8;
    abstol.data[1] = 1e-14;
    abstol.data[2] = 1e-6;

    /* Call CVodeCreate to create the solver memory (BDF) */
    let mut cvode_mem = CVodeCreate(CV_BDF, &sunctx);

    /* Call CVodeInit to initialize the integrator memory */
    let retval = CVodeInit(&mut cvode_mem, f, T0, &y0);
    if check_retval(retval, "CVodeInit") {
        return;
    }

    /* Call CVodeSVtolerances */
    let retval = CVodeSVtolerances(&mut cvode_mem, reltol, &abstol);
    if check_retval(retval, "CVodeSVtolerances") {
        return;
    }

    /* Call CVodeSetUserData so the sensitivity params can be accessed
     * from user provided routines (FSAUserData wrapper for the DQ run) */
    let retval = CVodeSetUserData(
        &mut cvode_mem,
        Some(Box::new(FSAUserData {
            p: p0.to_vec(),
            user: Box::new(()),
        })),
    );
    if check_retval(retval, "CVodeSetUserData") {
        return;
    }

    /* Call CVodeSetMaxNumSteps */
    let retval = CVodeSetMaxNumSteps(&mut cvode_mem, MXSTEPS);
    if check_retval(retval, "CVodeSetMaxNumSteps") {
        return;
    }

    /* Create dense SUNMatrix and SUNLinearSolver; attach */
    let a_mat = SUNDenseMatrix(NEQ, NEQ, &sunctx);
    let ls = SUNLinSol_Dense(&y, &a_mat, &sunctx);
    let retval = CVodeSetLinearSolver(&mut cvode_mem, ls, Some(a_mat));
    if check_retval(retval, "CVodeSetLinearSolver") {
        return;
    }

    /* Specify the Jacobian approximation routine to be used */
    let retval = CVodeSetJacFn(&mut cvode_mem, Some(jac));
    if check_retval(retval, "CVodeSetJacFn") {
        return;
    }

    /* Sensitivity-related settings */
    let mut config = RunConfig {
        sensi: true,             /* sensitivity ON                */
        meth: CV_SIMULTANEOUS,   /* simultaneous corrector method */
        errconS: true,           /* full error control            */
        fsDQ: false,             /* user-provided sensitivity RHS */
    };

    let ns = 3;

    let pbar = [p0[0], p0[1], p0[2]];
    let plist: Vec<i32> = (0..ns).collect();

    let mut yS0: Vec<NVector> = (0..ns).map(|_| N_VClone(&y)).collect();
    for ys in yS0.iter_mut() {
        N_VConst(ZERO, ys);
    }

    let retval = CVodeSensInit1(&mut cvode_mem, ns, config.meth, Some(fS), &yS0);
    if check_retval(retval, "CVodeSensInit1") {
        return;
    }

    let retval = CVodeSetSensParams(&mut cvode_mem, Some(&p0), Some(&pbar), Some(&plist));
    if check_retval(retval, "CVodeSetSensParams") {
        return;
    }

    /*
      Sensitivities are enabled
      Set full error control
      Set user-provided sensitivity RHS
      Run CVODES
    */
    let retval = CVodeSensEEtolerances(&mut cvode_mem);
    if check_retval(retval, "CVodeSensEEtolerances") {
        return;
    }

    let retval = CVodeSetSensErrCon(&mut cvode_mem, config.errconS);
    if check_retval(retval, "CVodeSetSensErrCon") {
        return;
    }

    let retval = run_cvode(&mut cvode_mem, &mut y, &config);
    if check_retval(retval, "runCVode") {
        return;
    }

    /*
      Change parameters
      Toggle sensitivities OFF
      Reinitialize and run CVODES
    */
    set_params(&mut cvode_mem, [0.05, 2.0e4, 2.9e7]);
    config.sensi = false;

    let retval = CVodeReInit(&mut cvode_mem, T0, &y0);
    if check_retval(retval, "CVodeReInit") {
        return;
    }

    let retval = CVodeSensToggleOff(&mut cvode_mem);
    if check_retval(retval, "CVodeSensToggleOff") {
        return;
    }

    let retval = run_cvode(&mut cvode_mem, &mut y, &config);
    if check_retval(retval, "runCVode") {
        return;
    }

    /*
      Change parameters
      Switch to internal DQ sensitivity RHS function
      Toggle sensitivities ON (reinitialize sensitivities)
      Reinitialize and run CVODES
    */
    set_params(&mut cvode_mem, [0.06, 3.0e4, 2.8e7]);
    config.sensi = true;
    config.fsDQ = true;

    let retval = CVodeReInit(&mut cvode_mem, T0, &y0);
    if check_retval(retval, "CVodeReInit") {
        return;
    }

    CVodeSensFree(&mut cvode_mem);
    let retval = CVodeSensInit1(&mut cvode_mem, ns, config.meth, None, &yS0);
    if check_retval(retval, "CVodeSensInit1") {
        return;
    }

    let retval = run_cvode(&mut cvode_mem, &mut y, &config);
    if check_retval(retval, "runCVode") {
        return;
    }

    /*
      Switch to partial error control
      Switch back to user-provided sensitivity RHS
      Toggle sensitivities ON (reinitialize sensitivities)
      Change method to staggered
      Reinitialize and run CVODES
    */
    config.sensi = true;
    config.errconS = false;
    config.fsDQ = false;
    config.meth = CV_STAGGERED;

    let retval = CVodeReInit(&mut cvode_mem, T0, &y0);
    if check_retval(retval, "CVodeReInit") {
        return;
    }

    let retval = CVodeSetSensErrCon(&mut cvode_mem, config.errconS);
    if check_retval(retval, "CVodeSetSensErrCon") {
        return;
    }

    CVodeSensFree(&mut cvode_mem);
    let retval = CVodeSensInit1(&mut cvode_mem, ns, config.meth, Some(fS), &yS0);
    if check_retval(retval, "CVodeSensInit1") {
        return;
    }

    let retval = run_cvode(&mut cvode_mem, &mut y, &config);
    if check_retval(retval, "runCVode") {
        return;
    }

    /*
      Free sensitivity-related memory
      (CVodeSensToggle is not needed, as CVodeSensFree toggles
      sensitivities OFF)
      Reinitialize and run CVODES
    */
    config.sensi = false;

    CVodeSensFree(&mut cvode_mem);

    let retval = CVodeReInit(&mut cvode_mem, T0, &y0);
    if check_retval(retval, "CVodeReInit") {
        return;
    }

    let retval = run_cvode(&mut cvode_mem, &mut y, &config);
    if check_retval(retval, "runCVode") {
        return;
    }

    /* Free memory */
    drop(y0);
    drop(y);
    drop(abstol);
    drop(yS0);
    CVodeFree(cvode_mem);
}
