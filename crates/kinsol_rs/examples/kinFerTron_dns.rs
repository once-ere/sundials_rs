/* -----------------------------------------------------------------
 * Translated from examples/kinsol/serial/kinFerTron_dns.c (KINSOL 7.7.0)
 *
 * Example (serial):
 *
 * This example solves a nonlinear system from.
 *
 * Source: "Handbook of Test Problems in Local and Global Optimization",
 *             C.A. Floudas, P.M. Pardalos et al.
 *             Kluwer Academic Publishers, 1999.
 * Test problem 4 from Section 14.1, Chapter 14: Ferraris and Tronconi
 *
 * This problem involves a blend of trigonometric and exponential terms.
 *    0.5 sin(x1 x2) - 0.25 x2/pi - 0.5 x1 = 0
 *    (1-0.25/pi) ( exp(2 x1)-e ) + e x2 / pi - 2 e x1 = 0
 * such that
 *    0.25 <= x1 <=1.0
 *    1.5 <= x2 <= 2 pi
 *
 * The treatment of the bound constraints on x1 and x2 is done using
 * the additional variables
 *    l1 = x1 - x1_min >= 0
 *    L1 = x1 - x1_max <= 0
 *    l2 = x2 - x2_min >= 0
 *    L2 = x2 - x2_max >= 0
 *
 * and using the constraint feature in KINSOL to impose
 *    l1 >= 0    l2 >= 0
 *    L1 <= 0    L2 <= 0
 *
 * The Ferraris-Tronconi test problem has two known solutions.
 * The nonlinear system is solved by KINSOL using different
 * combinations of globalization and Jacobian update strategies
 * and with different initial guesses (leading to one or the other
 * of the known solutions).
 *
 * Constraints are imposed to make all components of the solution
 * positive.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use kinsol_rs::sundials_utils::fmt_g;
use kinsol_rs::*;

/* Problem Constants */

const NVAR: i64 = 2;
const NEQ: i64 = 3 * NVAR;

const FTOL: f64 = 1.0e-5; /* function tolerance */
const STOL: f64 = 1.0e-5; /* step tolerance     */

const ZERO: f64 = 0.0;
const PT25: f64 = 0.25;
const PT5: f64 = 0.5;
const ONE: f64 = 1.0;
const ONEPT5: f64 = 1.5;
const TWO: f64 = 2.0;

const PI: f64 = 3.1415926;
const E: f64 = 2.7182818;

struct UserDataStruct {
    lb: [f64; NVAR as usize],
    ub: [f64; NVAR as usize],
}

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

fn main() {
    /* Create the SUNDIALS context that all SUNDIALS objects require */
    let sunctx = SUNContext_Create();

    /* User data */
    let data = UserDataStruct {
        lb: [PT25, ONEPT5],
        ub: [ONE, TWO * PI],
    };

    /* Create serial vectors of length NEQ */
    let mut u1 = N_VNew_Serial(NEQ, &sunctx);

    let mut u2 = N_VNew_Serial(NEQ, &sunctx);

    let mut u = N_VNew_Serial(NEQ, &sunctx);

    let mut s = N_VNew_Serial(NEQ, &sunctx);

    let mut c = N_VNew_Serial(NEQ, &sunctx);

    SetInitialGuess1(&mut u1, &data);
    SetInitialGuess2(&mut u2, &data);

    N_VConst(ONE, &mut s); /* no scaling */

    let cdata = N_VGetArrayPointer(&mut c);
    cdata[0] = ZERO; /* no constraint on x1 */
    cdata[1] = ZERO; /* no constraint on x2 */
    cdata[2] = ONE; /* l1 = x1 - x1_min >= 0 */
    cdata[3] = -ONE; /* L1 = x1 - x1_max <= 0 */
    cdata[4] = ONE; /* l2 = x2 - x2_min >= 0 */
    cdata[5] = -ONE; /* L2 = x2 - x22_min <= 0 */

    let fnormtol = FTOL; /* residual tolerance    */
    let scsteptol = STOL; /* scaled step tolerance */

    let mut kmem = KINCreate(&sunctx);

    let mut retval = KINSetUserData(&mut kmem, Some(Box::new(data)));
    if check_retval(retval, "KINSetUserData") {
        std::process::exit(1);
    }

    retval = KINSetConstraints(&mut kmem, Some(&c));
    if check_retval(retval, "KINSetConstraints") {
        std::process::exit(1);
    }

    retval = KINSetFuncNormTol(&mut kmem, fnormtol);
    if check_retval(retval, "KINSetFuncNormTol") {
        std::process::exit(1);
    }

    retval = KINSetScaledStepTol(&mut kmem, scsteptol);
    if check_retval(retval, "KINSetScaledStepTol") {
        std::process::exit(1);
    }

    retval = KINInit(&mut kmem, func, &u);
    if check_retval(retval, "KINInit") {
        std::process::exit(1);
    }

    /* Create dense SUNMatrix */
    let J = SUNDenseMatrix(NEQ, NEQ, &sunctx);

    /* Create dense SUNLinearSolver object */
    let LS = SUNLinSol_Dense(&u, &J, &sunctx);

    /* Attach the matrix and linear solver to KINSOL */
    retval = KINSetLinearSolver(&mut kmem, LS, Some(J));
    if check_retval(retval, "KINSetLinearSolver") {
        std::process::exit(1);
    }

    /* Print out the problem size, solution parameters, initial guess. */
    PrintHeader(fnormtol, scsteptol);

    /* --------------------------- */

    let mut glstr: i32; /* KINSOL globalization strategy flag */
    let mut mset: i64; /* KINSOL method selection flag */

    println!("\n------------------------------------------");
    println!("\nInitial guess on lower bounds");
    print!("  [x1,x2] = ");
    PrintOutput(&u1);

    N_VScale(ONE, &u1, &mut u);
    glstr = KIN_NONE;
    mset = 1;
    SolveIt(&mut kmem, &mut u, &s, glstr, mset);

    /* --------------------------- */

    N_VScale(ONE, &u1, &mut u);
    glstr = KIN_LINESEARCH;
    mset = 1;
    SolveIt(&mut kmem, &mut u, &s, glstr, mset);

    /* --------------------------- */

    N_VScale(ONE, &u1, &mut u);
    glstr = KIN_NONE;
    mset = 0;
    SolveIt(&mut kmem, &mut u, &s, glstr, mset);

    /* --------------------------- */

    N_VScale(ONE, &u1, &mut u);
    glstr = KIN_LINESEARCH;
    mset = 0;
    SolveIt(&mut kmem, &mut u, &s, glstr, mset);

    /* --------------------------- */

    println!("\n------------------------------------------");
    println!("\nInitial guess in middle of feasible region");
    print!("  [x1,x2] = ");
    PrintOutput(&u2);

    N_VScale(ONE, &u2, &mut u);
    glstr = KIN_NONE;
    mset = 1;
    SolveIt(&mut kmem, &mut u, &s, glstr, mset);

    /* --------------------------- */

    N_VScale(ONE, &u2, &mut u);
    glstr = KIN_LINESEARCH;
    mset = 1;
    SolveIt(&mut kmem, &mut u, &s, glstr, mset);

    /* --------------------------- */

    N_VScale(ONE, &u2, &mut u);
    glstr = KIN_NONE;
    mset = 0;
    SolveIt(&mut kmem, &mut u, &s, glstr, mset);

    /* --------------------------- */

    N_VScale(ONE, &u2, &mut u);
    glstr = KIN_LINESEARCH;
    mset = 0;
    SolveIt(&mut kmem, &mut u, &s, glstr, mset);

    /* Free memory */

    KINFree(kmem);

    std::process::exit(0);
}

fn SolveIt(kmem: &mut KINMem, u: &mut NVector, s: &NVector, glstr: i32, mset: i64) -> i32 {
    println!();

    if mset == 1 {
        print!("Exact Newton");
    } else {
        print!("Modified Newton");
    }

    if glstr == KIN_NONE {
        println!();
    } else {
        println!(" with line search");
    }

    let mut retval = KINSetMaxSetupCalls(kmem, mset);
    if check_retval(retval, "KINSetMaxSetupCalls") {
        return 1;
    }

    retval = KINSol(kmem, u, glstr, s, s);
    if check_retval(retval, "KINSol") {
        return 1;
    }

    print!("Solution:\n  [x1,x2] = ");
    PrintOutput(u);

    PrintFinalStats(kmem);

    0
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY KINSOL
 *--------------------------------------------------------------------
 */

/*
 * System function for predator-prey system
 */

fn func(u: &NVector, f: &mut NVector, user_data: &mut UserData) -> i32 {
    let data = user_data
        .as_mut()
        .unwrap()
        .downcast_mut::<UserDataStruct>()
        .unwrap();
    let lb = &data.lb;
    let ub = &data.ub;

    let udata = &u.data;
    let fdata = &mut f.data;

    let x1 = udata[0];
    let x2 = udata[1];
    let l1 = udata[2];
    let L1 = udata[3];
    let l2 = udata[4];
    let L2 = udata[5];

    fdata[0] = PT5 * (x1 * x2).sin() - PT25 * x2 / PI - PT5 * x1;
    fdata[1] = (ONE - PT25 / PI) * ((TWO * x1).exp() - E) + E * x2 / PI - TWO * E * x1;
    fdata[2] = l1 - x1 + lb[0];
    fdata[3] = L1 - x1 + ub[0];
    fdata[4] = l2 - x2 + lb[1];
    fdata[5] = L2 - x2 + ub[1];

    0
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * Initial guesses
 */

fn SetInitialGuess1(u: &mut NVector, data: &UserDataStruct) {
    let udata = &mut u.data;

    let lb = &data.lb;
    let ub = &data.ub;

    /* There are two known solutions for this problem */

    /* this init. guess should take us to (0.29945; 2.83693) */
    let x1 = lb[0];
    let x2 = lb[1];

    udata[0] = x1;
    udata[1] = x2;
    udata[2] = x1 - lb[0];
    udata[3] = x1 - ub[0];
    udata[4] = x2 - lb[1];
    udata[5] = x2 - ub[1];
}

fn SetInitialGuess2(u: &mut NVector, data: &UserDataStruct) {
    let udata = &mut u.data;

    let lb = &data.lb;
    let ub = &data.ub;

    /* There are two known solutions for this problem */

    /* this init. guess should take us to (0.5; 3.1415926) */
    let x1 = PT5 * (lb[0] + ub[0]);
    let x2 = PT5 * (lb[1] + ub[1]);

    udata[0] = x1;
    udata[1] = x2;
    udata[2] = x1 - lb[0];
    udata[3] = x1 - ub[0];
    udata[4] = x2 - lb[1];
    udata[5] = x2 - ub[1];
}

/*
 * Print first lines of output (problem description)
 */

fn PrintHeader(fnormtol: f64, scsteptol: f64) {
    println!("\nFerraris and Tronconi test problem");
    println!("Tolerance parameters:");
    println!(
        "  fnormtol  = {}\n  scsteptol = {}",
        fmt_g(fnormtol, 10, 6),
        fmt_g(scsteptol, 10, 6)
    );
}

/*
 * Print solution
 */

fn PrintOutput(u: &NVector) {
    let udata = &u.data;
    println!(" {}  {}", fmt_g(udata[0], 8, 6), fmt_g(udata[1], 8, 6));
}

/*
 * Print final statistics contained in iopt
 */

fn PrintFinalStats(kmem: &mut KINMem) {
    let mut nni: i64 = 0;
    let mut nfe: i64 = 0;
    let mut nje: i64 = 0;
    let mut nfeD: i64 = 0;

    let mut retval = KINGetNumNonlinSolvIters(kmem, &mut nni);
    check_retval(retval, "KINGetNumNonlinSolvIters");
    retval = KINGetNumFuncEvals(kmem, &mut nfe);
    check_retval(retval, "KINGetNumFuncEvals");
    retval = KINGetNumJacEvals(kmem, &mut nje);
    check_retval(retval, "KINGetNumJacEvals");
    retval = KINGetNumLinFuncEvals(kmem, &mut nfeD);
    check_retval(retval, "KINGetNumLinFuncEvals");

    println!("Final Statistics:");
    println!("  nni = {:5}    nfe  = {:5} ", nni, nfe);
    println!("  nje = {:5}    nfeD = {:5} ", nje, nfeD);
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
