/* -----------------------------------------------------------------
 * Translated from examples/kinsol/serial/kinRoberts_fp.c (KINSOL 7.7.0)
 *
 * Example problem:
 *
 * The following is a simple example problem, with the coding
 * needed for its solution by the accelerated fixed point solver in
 * KINSOL.
 * The problem is from chemical kinetics, and consists of solving
 * the first time step in a Backward Euler solution for the
 * following three rate equations:
 *    dy1/dt = -.04*y1 + 1.e4*y2*y3
 *    dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e2*(y2)^2
 *    dy3/dt = 3.e2*(y2)^2
 * on the interval from t = 0.0 to t = 0.1, with initial
 * conditions: y1 = 1.0, y2 = y3 = 0. The problem is stiff.
 * Run statistics (optional outputs) are printed at the end.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use kinsol_rs::sundials_utils::{fmt_e, fmt_g};
use kinsol_rs::*;

/* Problem Constants */

const NEQ: i64 = 3; /* number of equations  */
const Y10: f64 = 1.0; /* initial y components */
const Y20: f64 = 0.0;
const Y30: f64 = 0.0;
const TOL: f64 = 1.0e-10; /* function tolerance */
const DSTEP: f64 = 0.1; /* Size of the single time step used */

const PRIORS: i64 = 2;

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/* Ith numbers components 1..NEQ (C macro) */
macro_rules! Ith {
    ($v:expr, $i:expr) => {
        $v.data[$i - 1]
    };
}

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

fn main() {
    let mut fnorm: f64 = 0.0;

    /* -------------------------
     * Print problem description
     * ------------------------- */

    println!("Example problem from chemical kinetics solving");
    println!("the first time step in a Backward Euler solution for the");
    println!("following three rate equations:");
    println!("    dy1/dt = -.04*y1 + 1.e4*y2*y3");
    println!("    dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e2*(y2)^2");
    println!("    dy3/dt = 3.e2*(y2)^2");
    println!("on the interval from t = 0.0 to t = 0.1, with initial");
    println!("conditions: y1 = 1.0, y2 = y3 = 0.");
    println!("Solution method: Anderson accelerated fixed point iteration.");

    /* Create the SUNDIALS context that all SUNDIALS objects require */
    let sunctx = SUNContext_Create();

    /* --------------------------------------
     * Create vectors for solution and scales
     * -------------------------------------- */

    let mut y = N_VNew_Serial(NEQ, &sunctx);

    let mut scale = N_VNew_Serial(NEQ, &sunctx);

    /* -----------------------------------------
     * Initialize and allocate memory for KINSOL
     * ----------------------------------------- */

    let mut kmem = KINCreate(&sunctx);

    /* y is used as a template */

    /* Set number of prior residuals used in Anderson acceleration */
    let _ = KINSetMAA(&mut kmem, PRIORS);

    let mut retval = KINInit(&mut kmem, funcRoberts, &y);
    if check_retval(retval, "KINInit") {
        std::process::exit(1);
    }

    /* -------------------
     * Set optional inputs
     * ------------------- */

    /* Specify stopping tolerance based on residual */

    let fnormtol = TOL;
    retval = KINSetFuncNormTol(&mut kmem, fnormtol);
    if check_retval(retval, "KINSetFuncNormTol") {
        std::process::exit(1);
    }

    /* Override any current settings with command-line options */

    let args: Vec<String> = std::env::args().collect();
    retval = KINSetOptions(&mut kmem, "", "", &args);
    if check_retval(retval, "KINSetOptions") {
        std::process::exit(1);
    }

    /* -------------
     * Initial guess
     * ------------- */

    N_VConst(ZERO, &mut y);
    Ith!(y, 1) = ONE;

    /* ----------------------------
     * Call KINSol to solve problem
     * ---------------------------- */

    /* No scaling used */
    N_VConst(ONE, &mut scale);

    /* Call main solver */
    retval = KINSol(&mut kmem, /* KINSol memory block */
                    &mut y,    /* initial guess on input; solution vector */
                    KIN_FP,    /* global strategy choice */
                    &scale,    /* scaling vector, for the variable cc */
                    &scale);   /* scaling vector for function values fval */
    if check_retval(retval, "KINSol") {
        std::process::exit(1);
    }

    /* ------------------------------------
     * Print solution and solver statistics
     * ------------------------------------ */

    /* Get scaled norm of the system function */

    retval = KINGetFuncNorm(&mut kmem, &mut fnorm);
    if check_retval(retval, "KINGetfuncNorm") {
        std::process::exit(1);
    }

    println!("\nComputed solution (||F|| = {}):\n", fmt_g(fnorm, 0, 6));
    PrintOutput(&y);

    PrintFinalStats(&mut kmem);

    /* check the solution error */
    let retval = check_ans(&y, 1e-4, 1e-6);

    /* -----------
     * Free memory
     * ----------- */

    KINFree(kmem);

    std::process::exit(retval);
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * System function
 */

fn funcRoberts(y: &NVector, g: &mut NVector, _user_data: &mut UserData) -> i32 {
    let y1 = Ith!(y, 1);
    let y2 = Ith!(y, 2);
    let y3 = Ith!(y, 3);

    let yd1 = DSTEP * (-0.04 * y1 + 1.0e4 * y2 * y3);
    let yd3 = DSTEP * 3.0e2 * y2 * y2;

    Ith!(g, 1) = yd1 + Y10;
    Ith!(g, 2) = -yd1 - yd3 + Y20;
    Ith!(g, 3) = yd3 + Y30;

    0
}

/*
 * Print solution at selected points
 */

fn PrintOutput(y: &NVector) {
    let y1 = Ith!(y, 1);
    let y2 = Ith!(y, 2);
    let y3 = Ith!(y, 3);

    println!(
        "y ={}  {}  {}",
        fmt_e(y1, 14, 6),
        fmt_e(y2, 14, 6),
        fmt_e(y3, 14, 6)
    );
}

/*
 * Print final statistics
 */

fn PrintFinalStats(kmem: &mut KINMem) {
    let mut nni: i64 = 0;
    let mut nfe: i64 = 0;

    /* Main solver statistics */

    let mut retval = KINGetNumNonlinSolvIters(kmem, &mut nni);
    check_retval(retval, "KINGetNumNonlinSolvIters");
    retval = KINGetNumFuncEvals(kmem, &mut nfe);
    check_retval(retval, "KINGetNumFuncEvals");

    println!("\nFinal Statistics.. \n");
    println!("nni      = {:6}    nfe     = {:6} ", nni, nfe);
}

/*
 * Check function return value...
 *    opt == 1 means SUNDIALS function returns a retval so check if
 *             retval >= 0
 *    (the opt == 0 / opt == 2 NULL-pointer cases do not arise in the
 *     Rust port)
 */

fn check_retval(retval: i32, funcname: &str) -> bool {
    if retval < 0 {
        eprintln!("\nSUNDIALS_ERROR: {}() failed with retval = {}\n", funcname, retval);
        return true;
    }
    false
}

/* compare the solution to a reference solution computed with a
   tolerance of 1e-14 */
fn check_ans(u: &NVector, rtol: f64, atol: f64) -> i32 {
    /* create reference solution and error weight vectors */
    let mut refv = N_VClone(u);
    let mut ewt = N_VClone(u);

    /* set the reference solution data */
    refv.data[0] = 9.9678538655358029e-01;
    refv.data[1] = 2.9530060962800345e-03;
    refv.data[2] = 2.6160735013975683e-04;

    /* compute the error weight vector */
    N_VAbs(&refv, &mut ewt);
    ewt.scale_inplace(rtol);
    ewt.add_const_inplace(atol);
    if N_VMin(&ewt) <= ZERO {
        eprintln!("\nSUNDIALS_ERROR: check_ans failed - ewt <= 0\n");
        return -1;
    }
    ewt.invert_inplace();

    /* compute the solution error */
    refv.linear_sum_with(-ONE, ONE, u);
    let err = N_VWrmsNorm(&refv, &ewt);

    /* is the solution within the tolerances? */
    let passfail = if err < ONE { 0 } else { 1 };

    if passfail != 0 {
        println!("\nSUNDIALS_WARNING: check_ans error={}\n", fmt_g(err, 0, 6));
    }

    passfail
}
