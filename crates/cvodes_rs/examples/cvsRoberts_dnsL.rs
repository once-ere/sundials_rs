/* -----------------------------------------------------------------
 * Translated from examples/cvode/serial/cvsRoberts_dnsL.c (CVODE 7.7.0)
 *
 * Chemical kinetics 3-species problem (see cvRoberts_dns), solved
 * with the BDF method, Newton iteration and the *LAPACK* dense
 * linear solver in the original. The pure-Rust build has no LAPACK;
 * SUNLinSol_LapackDense maps onto the native dense LU
 * (SUNLinSol_Dense), which performs the same factorization up to
 * floating-point ordering inside the elimination.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use cvodes_rs::sundials_utils::fmt_e;
use cvodes_rs::*;

const NEQ: usize = 3;
const Y1: f64 = 1.0;
const Y2: f64 = 0.0;
const Y3: f64 = 0.0;
const RTOL: f64 = 1.0e-4;
const ATOL1: f64 = 1.0e-8;
const ATOL2: f64 = 1.0e-14;
const ATOL3: f64 = 1.0e-6;
const T0: f64 = 0.0;
const T1: f64 = 0.4;
const TMULT: f64 = 10.0;
const NOUT: i32 = 12;

const ZERO: f64 = 0.0;

macro_rules! Ith {
    ($v:expr, $i:expr) => {
        $v.data[$i - 1]
    };
}

fn f(_t: f64, y: &NVector, ydot: &mut NVector, _user_data: &mut UserData) -> i32 {
    let y1 = Ith!(y, 1);
    let y2 = Ith!(y, 2);
    let y3 = Ith!(y, 3);

    let yd1 = -0.04 * y1 + 1.0e4 * y2 * y3;
    Ith!(ydot, 1) = yd1;
    let yd3 = 3.0e7 * y2 * y2;
    Ith!(ydot, 3) = yd3;
    Ith!(ydot, 2) = -yd1 - yd3;

    0
}

fn g(_t: f64, y: &NVector, gout: &mut [f64], _user_data: &mut UserData) -> i32 {
    gout[0] = Ith!(y, 1) - 0.0001;
    gout[1] = Ith!(y, 3) - 0.01;
    0
}

fn Jac(
    _t: f64,
    y: &NVector,
    _fy: &NVector,
    j: &mut SUNMatrix,
    _user_data: &mut UserData,
    _tmp1: &mut NVector,
    _tmp2: &mut NVector,
    _tmp3: &mut NVector,
) -> i32 {
    let y2 = Ith!(y, 2);
    let y3 = Ith!(y, 3);

    let jm = match j {
        SUNMatrix::Dense(m) => m,
        _ => return -1,
    };
    jm.set(0, 0, -0.04);
    jm.set(0, 1, 1.0e4 * y3);
    jm.set(0, 2, 1.0e4 * y2);
    jm.set(1, 0, 0.04);
    jm.set(1, 1, -1.0e4 * y3 - 6.0e7 * y2);
    jm.set(1, 2, -1.0e4 * y2);
    jm.set(2, 0, ZERO);
    jm.set(2, 1, 6.0e7 * y2);
    jm.set(2, 2, ZERO);

    0
}

fn PrintOutput(t: f64, y1: f64, y2: f64, y3: f64) {
    println!(
        "At t = {}      y ={}  {}  {}",
        fmt_e(t, 0, 4),
        fmt_e(y1, 14, 6),
        fmt_e(y2, 14, 6),
        fmt_e(y3, 14, 6)
    );
}

fn PrintRootInfo(root_f1: i32, root_f2: i32) {
    println!("    rootsfound[] = {:3} {:3}", root_f1, root_f2);
}

fn PrintFinalStats(cvode_mem: &mut CVodeMem) {
    let (mut nst, mut nfe, mut nsetups, mut nje, mut nfeLS) = (0i64, 0i64, 0i64, 0i64, 0i64);
    let (mut nni, mut nnf, mut ncfn, mut netf, mut nge) = (0i64, 0i64, 0i64, 0i64, 0i64);

    CVodeGetNumSteps(cvode_mem, &mut nst);
    CVodeGetNumRhsEvals(cvode_mem, &mut nfe);
    CVodeGetNumLinSolvSetups(cvode_mem, &mut nsetups);
    CVodeGetNumErrTestFails(cvode_mem, &mut netf);
    CVodeGetNumNonlinSolvIters(cvode_mem, &mut nni);
    CVodeGetNumNonlinSolvConvFails(cvode_mem, &mut nnf);
    CVodeGetNumStepSolveFails(cvode_mem, &mut ncfn);
    CVodeGetNumJacEvals(cvode_mem, &mut nje);
    CVodeGetNumLinRhsEvals(cvode_mem, &mut nfeLS);
    CVodeGetNumGEvals(cvode_mem, &mut nge);

    println!("\nFinal Statistics:");
    println!(
        "nst = {:<6} nfe = {:<6} nsetups = {:<6} nfeLS = {:<6} nje = {}",
        nst, nfe, nsetups, nfeLS, nje
    );
    println!(
        "nni = {:<6} nnf = {:<6} netf = {:<6}    ncfn = {:<6}  nge = {}\n",
        nni, nnf, netf, ncfn, nge
    );
}

fn check_retval(retval: i32, funcname: &str) -> bool {
    if retval < 0 {
        eprintln!("\nSUNDIALS_ERROR: {}() failed with retval = {}\n", funcname, retval);
        return true;
    }
    false
}

fn main() {
    let sunctx = SUNContext_Create();

    let mut y = N_VNew_Serial(NEQ as i64, &sunctx);
    Ith!(y, 1) = Y1;
    Ith!(y, 2) = Y2;
    Ith!(y, 3) = Y3;

    let mut abstol = N_VNew_Serial(NEQ as i64, &sunctx);
    Ith!(abstol, 1) = ATOL1;
    Ith!(abstol, 2) = ATOL2;
    Ith!(abstol, 3) = ATOL3;

    let mut cvode_mem = CVodeCreate(CV_BDF, &sunctx);

    let mut retval = CVodeInit(&mut cvode_mem, f, T0, &y);
    if check_retval(retval, "CVodeInit") {
        std::process::exit(1);
    }

    retval = CVodeSVtolerances(&mut cvode_mem, RTOL, &abstol);
    if check_retval(retval, "CVodeSVtolerances") {
        std::process::exit(1);
    }

    retval = CVodeRootInit(&mut cvode_mem, 2, Some(g));
    if check_retval(retval, "CVodeRootInit") {
        std::process::exit(1);
    }

    /* Create dense SUNMatrix for use in linear solves */
    let a = SUNDenseMatrix(NEQ as i64, NEQ as i64, &sunctx);

    /* Create (Lapack)Dense solver object — native dense LU here */
    let ls = SUNLinSol_Dense(&y, &a, &sunctx);

    retval = CVodeSetLinearSolver(&mut cvode_mem, ls, Some(a));
    if check_retval(retval, "CVodeSetLinearSolver") {
        std::process::exit(1);
    }

    retval = CVodeSetJacFn(&mut cvode_mem, Some(Jac));
    if check_retval(retval, "CVodeSetJacFn") {
        std::process::exit(1);
    }

    println!(" \n3-species kinetics problem\n");

    let mut iout = 0;
    let mut tout = T1;
    let mut t = 0.0;
    let mut rootsfound = [0i32; 2];
    loop {
        retval = CVode(&mut cvode_mem, tout, &mut y, &mut t, CV_NORMAL);
        PrintOutput(t, Ith!(y, 1), Ith!(y, 2), Ith!(y, 3));

        if retval == CV_ROOT_RETURN {
            let retvalr = CVodeGetRootInfo(&mut cvode_mem, &mut rootsfound);
            if check_retval(retvalr, "CVodeGetRootInfo") {
                std::process::exit(1);
            }
            PrintRootInfo(rootsfound[0], rootsfound[1]);
        }

        if check_retval(retval, "CVode") {
            break;
        }
        if retval == CV_SUCCESS {
            iout += 1;
            tout *= TMULT;
        }

        if iout == NOUT {
            break;
        }
    }

    PrintFinalStats(&mut cvode_mem);

    CVodeFree(cvode_mem);
}
