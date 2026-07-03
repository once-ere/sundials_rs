/* -----------------------------------------------------------------
 * Translated from examples/cvode/serial/cvRoberts_dns.c (CVODE 7.7.0)
 *
 * Example problem:
 *
 * The following is a simple example problem, with the coding
 * needed for its solution by CVODE. The problem is from
 * chemical kinetics, and consists of the following three rate
 * equations:
 *    dy1/dt = -.04*y1 + 1.e4*y2*y3
 *    dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*(y2)^2
 *    dy3/dt = 3.e7*(y2)^2
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1.0, y2 = y3 = 0. The problem is stiff.
 * While integrating the system, we also use the rootfinding
 * feature to find the points at which y1 = 1e-4 or at which
 * y3 = 0.01. This program solves the problem with the BDF method,
 * Newton iteration with the dense linear solver, and a
 * user-supplied Jacobian routine.
 * It uses a scalar relative tolerance and a vector absolute
 * tolerance. Output is printed in decades from t = .4 to t = 4.e10.
 * Run statistics (optional outputs) are printed at the end.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use cvode_rs::sundials_utils::fmt_e;
use cvode_rs::*;

/* Problem Constants */
const NEQ: usize = 3; /* number of equations  */
const Y1: f64 = 1.0; /* initial y components */
const Y2: f64 = 0.0;
const Y3: f64 = 0.0;
const RTOL: f64 = 1.0e-4; /* scalar relative tolerance            */
const ATOL1: f64 = 1.0e-8; /* vector absolute tolerance components */
const ATOL2: f64 = 1.0e-14;
const ATOL3: f64 = 1.0e-6;
const T0: f64 = 0.0; /* initial time       */
const T1: f64 = 0.4; /* first output time  */
const TMULT: f64 = 10.0; /* output time factor */
const NOUT: i32 = 12; /* number of output times */

const ZERO: f64 = 0.0;

/* Ith numbers components 1..NEQ (C macro) */
macro_rules! Ith {
    ($v:expr, $i:expr) => {
        $v.data[$i - 1]
    };
}

/*
 * f routine. Compute function f(t,y).
 */
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

/*
 * g routine. Compute functions g_i(t,y) for i = 0,1.
 */
fn g(_t: f64, y: &NVector, gout: &mut [f64], _user_data: &mut UserData) -> i32 {
    let y1 = Ith!(y, 1);
    let y3 = Ith!(y, 3);
    gout[0] = y1 - 0.0001;
    gout[1] = y3 - 0.01;

    0
}

/*
 * Jacobian routine. Compute J(t,y) = df/dy.
 */
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
    /* IJth(A,i,j) is 1-based in the C example */
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

/*
 * Private helper functions
 */
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

fn check_retval(retval: i32, funcname: &str) -> bool {
    if retval < 0 {
        eprintln!("\nSUNDIALS_ERROR: {}() failed with retval = {}\n", funcname, retval);
        return true;
    }
    false
}

/* compare the solution at the final time 4e10s to a reference solution computed
   using a relative tolerance of 1e-8 and absolute tolerance of 1e-14 */
fn check_ans(y: &NVector, _t: f64, rtol: f64, atol: &NVector) -> i32 {
    let mut refv = N_VClone(y); /* reference solution vector */
    let mut ewt = N_VClone(y); /* error weight vector       */

    /* set the reference solution data */
    refv.data[0] = 5.2083495894337328e-08;
    refv.data[1] = 2.0833399429795671e-13;
    refv.data[2] = 9.9999994791629776e-01;

    /* compute the error weight vector, loosen atol */
    N_VAbs(&refv, &mut ewt);
    ewt.linear_sum_with(rtol, 10.0, atol);
    if N_VMin(&ewt) <= ZERO {
        eprintln!("\nSUNDIALS_ERROR: check_ans failed - ewt <= 0\n");
        return -1;
    }
    ewt.invert_inplace();

    /* compute the solution error: ref = y - ref */
    refv.linear_sum_with(-1.0, 1.0, y);
    let err = N_VWrmsNorm(&refv, &ewt);

    /* is the solution within the tolerances? */
    let passfail = if err < 1.0 { 0 } else { 1 };
    if passfail != 0 {
        eprintln!("\nSUNDIALS_WARNING: check_ans error={} \n", err);
    }
    passfail
}

fn main() {
    /* Create the SUNDIALS context */
    let sunctx = SUNContext_Create();

    /* Initial conditions */
    let mut y = N_VNew_Serial(NEQ as i64, &sunctx);
    Ith!(y, 1) = Y1;
    Ith!(y, 2) = Y2;
    Ith!(y, 3) = Y3;

    /* Set the vector absolute tolerance */
    let mut abstol = N_VNew_Serial(NEQ as i64, &sunctx);
    Ith!(abstol, 1) = ATOL1;
    Ith!(abstol, 2) = ATOL2;
    Ith!(abstol, 3) = ATOL3;

    /* Create the solver memory with the BDF method */
    let mut cvode_mem = CVodeCreate(CV_BDF, &sunctx);

    /* Initialize the integrator memory */
    let mut retval = CVodeInit(&mut cvode_mem, f, T0, &y);
    if check_retval(retval, "CVodeInit") {
        std::process::exit(1);
    }

    /* Specify the scalar relative tolerance and vector absolute tolerances */
    retval = CVodeSVtolerances(&mut cvode_mem, RTOL, &abstol);
    if check_retval(retval, "CVodeSVtolerances") {
        std::process::exit(1);
    }

    /* Specify the root function g with 2 components */
    retval = CVodeRootInit(&mut cvode_mem, 2, Some(g));
    if check_retval(retval, "CVodeRootInit") {
        std::process::exit(1);
    }

    /* Create dense SUNMatrix for use in linear solves */
    let a = SUNDenseMatrix(NEQ as i64, NEQ as i64, &sunctx);

    /* Create dense SUNLinearSolver object for use by CVode */
    let ls = SUNLinSol_Dense(&y, &a, &sunctx);

    /* Attach the matrix and linear solver */
    retval = CVodeSetLinearSolver(&mut cvode_mem, ls, Some(a));
    if check_retval(retval, "CVodeSetLinearSolver") {
        std::process::exit(1);
    }

    /* Set the user-supplied Jacobian routine Jac */
    retval = CVodeSetJacFn(&mut cvode_mem, Some(Jac));
    if check_retval(retval, "CVodeSetJacFn") {
        std::process::exit(1);
    }

    /* In loop, call CVode, print results, and test for error. */
    println!(" \n3-species kinetics problem\n");

    /* Open file for printing statistics */
    let mut fid = std::fs::File::create("cvRoberts_dns_stats.csv").expect("create csv");

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

        let _ = CVodePrintAllStats(&mut cvode_mem, &mut fid, SUN_OUTPUTFORMAT_CSV);

        if iout == NOUT {
            break;
        }
    }
    drop(fid);

    /* Print final statistics to the screen */
    println!("\nFinal Statistics:");
    let mut stdout = std::io::stdout();
    let _ = CVodePrintAllStats(&mut cvode_mem, &mut stdout, SUN_OUTPUTFORMAT_TABLE);

    /* check the solution error */
    let retval = check_ans(&y, t, RTOL, &abstol);

    /* Free memory (RAII) */
    CVodeFree(cvode_mem);

    std::process::exit(retval);
}
