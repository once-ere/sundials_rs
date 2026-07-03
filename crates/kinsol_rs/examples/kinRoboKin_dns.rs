/* -----------------------------------------------------------------
 * Translated from examples/kinsol/serial/kinRoboKin_dns.c (KINSOL 7.7.0)
 *
 * This example solves a nonlinear system from robot kinematics.
 *
 * Source: "Handbook of Test Problems in Local and Global Optimization",
 *             C.A. Floudas, P.M. Pardalos et al.
 *             Kluwer Academic Publishers, 1999.
 * Test problem 6 from Section 14.1, Chapter 14
 *
 * The nonlinear system is solved by KINSOL using the DENSE linear
 * solver.
 *
 * Constraints are imposed to make all components of the solution
 * be within [-1,1].
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use kinsol_rs::sundials_utils::fmt_g;
use kinsol_rs::*;

/* Problem Constants */

const NVAR: usize = 8; /* variables */
const NEQ: usize = 3 * NVAR; /* equations + bounds */

const FTOL: f64 = 1.0e-5; /* function tolerance */
const STOL: f64 = 1.0e-5; /* step tolerance */

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;

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
    println!("\nRobot Kinematics Example");
    println!("8 variables; -1 <= x_i <= 1");
    println!("KINSOL problem size: 8 + 2*8 = 24 \n");

    /* Create the SUNDIALS context that all SUNDIALS objects require */
    let sunctx = SUNContext_Create();

    /* Create vectors for solution, scales, and constraints */

    let mut y = N_VNew_Serial(NEQ as i64, &sunctx);

    let mut scale = N_VNew_Serial(NEQ as i64, &sunctx);

    let mut constraints = N_VNew_Serial(NEQ as i64, &sunctx);

    /* Initialize and allocate memory for KINSOL */

    let mut kmem = KINCreate(&sunctx);

    let mut retval = KINInit(&mut kmem, func, &y); /* y passed as a template */
    if check_retval(retval, "KINInit") {
        std::process::exit(1);
    }

    /* Set optional inputs */

    N_VConst(ZERO, &mut constraints);
    for i in (NVAR + 1)..=NEQ {
        Ith!(constraints, i) = ONE;
    }

    retval = KINSetConstraints(&mut kmem, Some(&constraints));
    if check_retval(retval, "KINSetConstraints") {
        std::process::exit(1);
    }

    let fnormtol = FTOL;
    retval = KINSetFuncNormTol(&mut kmem, fnormtol);
    if check_retval(retval, "KINSetFuncNormTol") {
        std::process::exit(1);
    }

    let scsteptol = STOL;
    retval = KINSetScaledStepTol(&mut kmem, scsteptol);
    if check_retval(retval, "KINSetScaledStepTol") {
        std::process::exit(1);
    }

    /* Create dense SUNMatrix */
    let J = SUNDenseMatrix(NEQ as i64, NEQ as i64, &sunctx);

    /* Create dense SUNLinearSolver object */
    let LS = SUNLinSol_Dense(&y, &J, &sunctx);

    /* Attach the matrix and linear solver to KINSOL */
    retval = KINSetLinearSolver(&mut kmem, LS, Some(J));
    if check_retval(retval, "KINSetLinearSolver") {
        std::process::exit(1);
    }

    /* Set the Jacobian function */
    retval = KINSetJacFn(&mut kmem, Some(jac));
    if check_retval(retval, "KINSetJacFn") {
        std::process::exit(1);
    }

    /* Indicate exact Newton */

    let mset: i64 = 1;
    retval = KINSetMaxSetupCalls(&mut kmem, mset);
    if check_retval(retval, "KINSetMaxSetupCalls") {
        std::process::exit(1);
    }

    /* Initial guess */

    N_VConst(ONE, &mut y);
    for i in 1..=NVAR {
        Ith!(y, i) = SUNRsqrt(TWO) / TWO;
    }

    println!("Initial guess:");
    PrintOutput(&y);

    /* Call KINSol to solve problem */

    N_VConst(ONE, &mut scale);
    retval = KINSol(&mut kmem,      /* KINSol memory block */
                    &mut y,         /* initial guess on input; solution vector */
                    KIN_LINESEARCH, /* global strategy choice */
                    &scale,         /* scaling vector, for the variable cc */
                    &scale);        /* scaling vector for function values fval */
    if check_retval(retval, "KINSol") {
        std::process::exit(1);
    }

    println!("\nComputed solution:");
    PrintOutput(&y);

    /* Print final statistics to screen and file */

    println!("\nFinal statsistics:");
    let mut stdout = std::io::stdout();
    let _retval = KINPrintAllStats(&mut kmem, &mut stdout, SUN_OUTPUTFORMAT_TABLE);

    let mut fid = std::fs::File::create("kinRoboKin_dns_stats.csv").expect("create csv");
    let _retval = KINPrintAllStats(&mut kmem, &mut fid, SUN_OUTPUTFORMAT_CSV);
    drop(fid);

    /* free memory */
    KINFree(kmem);

    std::process::exit(0);
}

/*
 * System function
 */

fn func(y: &NVector, f: &mut NVector, _user_data: &mut UserData) -> i32 {
    let yd = &y.data;
    let fd = &mut f.data;

    let x1 = yd[0];
    let l1 = yd[8];
    let u1 = yd[16];
    let x2 = yd[1];
    let l2 = yd[9];
    let u2 = yd[17];
    let x3 = yd[2];
    let l3 = yd[10];
    let u3 = yd[18];
    let x4 = yd[3];
    let l4 = yd[11];
    let u4 = yd[19];
    let x5 = yd[4];
    let l5 = yd[12];
    let u5 = yd[20];
    let x6 = yd[5];
    let l6 = yd[13];
    let u6 = yd[21];
    let x7 = yd[6];
    let l7 = yd[14];
    let u7 = yd[22];
    let x8 = yd[7];
    let l8 = yd[15];
    let u8 = yd[23];

    /* Nonlinear equations */

    let eq1 = -0.1238 * x1 + x7 - 0.001637 * x2 - 0.9338 * x4 + 0.004731 * x1 * x3
        - 0.3578 * x2 * x3
        - 0.3571;
    let eq2 = 0.2638 * x1 - x7 - 0.07745 * x2 - 0.6734 * x4 + 0.2238 * x1 * x3
        + 0.7623 * x2 * x3
        - 0.6022;
    let eq3 = 0.3578 * x1 + 0.004731 * x2 + x6 * x8;
    let eq4 = -0.7623 * x1 + 0.2238 * x2 + 0.3461;
    let eq5 = x1 * x1 + x2 * x2 - 1.0;
    let eq6 = x3 * x3 + x4 * x4 - 1.0;
    let eq7 = x5 * x5 + x6 * x6 - 1.0;
    let eq8 = x7 * x7 + x8 * x8 - 1.0;

    /* Lower bounds ( l_i = 1 + x_i >= 0)*/

    let lb1 = l1 - 1.0 - x1;
    let lb2 = l2 - 1.0 - x2;
    let lb3 = l3 - 1.0 - x3;
    let lb4 = l4 - 1.0 - x4;
    let lb5 = l5 - 1.0 - x5;
    let lb6 = l6 - 1.0 - x6;
    let lb7 = l7 - 1.0 - x7;
    let lb8 = l8 - 1.0 - x8;

    /* Upper bounds ( u_i = 1 - x_i >= 0)*/

    let ub1 = u1 - 1.0 + x1;
    let ub2 = u2 - 1.0 + x2;
    let ub3 = u3 - 1.0 + x3;
    let ub4 = u4 - 1.0 + x4;
    let ub5 = u5 - 1.0 + x5;
    let ub6 = u6 - 1.0 + x6;
    let ub7 = u7 - 1.0 + x7;
    let ub8 = u8 - 1.0 + x8;

    fd[0] = eq1;
    fd[8] = lb1;
    fd[16] = ub1;
    fd[1] = eq2;
    fd[9] = lb2;
    fd[17] = ub2;
    fd[2] = eq3;
    fd[10] = lb3;
    fd[18] = ub3;
    fd[3] = eq4;
    fd[11] = lb4;
    fd[19] = ub4;
    fd[4] = eq5;
    fd[12] = lb5;
    fd[20] = ub5;
    fd[5] = eq6;
    fd[13] = lb6;
    fd[21] = ub6;
    fd[6] = eq7;
    fd[14] = lb7;
    fd[22] = ub7;
    fd[7] = eq8;
    fd[15] = lb8;
    fd[23] = ub8;

    0
}

/*
 * System Jacobian
 */

fn jac(
    y: &NVector,
    _f: &NVector,
    J: &mut SUNMatrix,
    _user_data: &mut UserData,
    _tmp1: &mut NVector,
    _tmp2: &mut NVector,
) -> i32 {
    let yd = &y.data;

    let x1 = yd[0];
    let x2 = yd[1];
    let x3 = yd[2];
    let x4 = yd[3];
    let x5 = yd[4];
    let x6 = yd[5];
    let x7 = yd[6];
    let x8 = yd[7];

    let jm = match J {
        SUNMatrix::Dense(m) => m,
        _ => return -1,
    };

    /* IJth(J,i,j) is 1-based in the C example */
    macro_rules! IJth {
        ($i:expr, $j:expr, $v:expr) => {
            jm.set($i - 1, $j - 1, $v)
        };
    }

    /* Nonlinear equations */

    /*
       - 0.1238*x1 + x7 - 0.001637*x2
       - 0.9338*x4 + 0.004731*x1*x3 - 0.3578*x2*x3 - 0.3571
    */
    IJth!(1, 1, -0.1238 + 0.004731 * x3);
    IJth!(1, 2, -0.001637 - 0.3578 * x3);
    IJth!(1, 3, 0.004731 * x1 - 0.3578 * x2);
    IJth!(1, 4, -0.9338);
    IJth!(1, 7, 1.0);

    /*
      0.2638*x1 - x7 - 0.07745*x2
      - 0.6734*x4 + 0.2238*x1*x3 + 0.7623*x2*x3 - 0.6022
    */
    IJth!(2, 1, 0.2638 + 0.2238 * x3);
    IJth!(2, 2, -0.07745 + 0.7623 * x3);
    IJth!(2, 3, 0.2238 * x1 + 0.7623 * x2);
    IJth!(2, 4, -0.6734);
    IJth!(2, 7, -1.0);

    /*
      0.3578*x1 + 0.004731*x2 + x6*x8
    */
    IJth!(3, 1, 0.3578);
    IJth!(3, 2, 0.004731);
    IJth!(3, 6, x8);
    IJth!(3, 8, x6);

    /*
      - 0.7623*x1 + 0.2238*x2 + 0.3461
    */
    IJth!(4, 1, -0.7623);
    IJth!(4, 2, 0.2238);

    /*
      x1*x1 + x2*x2 - 1
    */
    IJth!(5, 1, 2.0 * x1);
    IJth!(5, 2, 2.0 * x2);

    /*
      x3*x3 + x4*x4 - 1
    */
    IJth!(6, 3, 2.0 * x3);
    IJth!(6, 4, 2.0 * x4);

    /*
      x5*x5 + x6*x6 - 1
    */
    IJth!(7, 5, 2.0 * x5);
    IJth!(7, 6, 2.0 * x6);

    /*
      x7*x7 + x8*x8 - 1
    */
    IJth!(8, 7, 2.0 * x7);
    IJth!(8, 8, 2.0 * x8);

    /*
      Lower bounds ( l_i = 1 + x_i >= 0)
      l_i - 1.0 - x_i
     */

    for i in 1..=8 {
        IJth!(8 + i, i, -1.0);
        IJth!(8 + i, 8 + i, 1.0);
    }

    /*
      Upper bounds ( u_i = 1 - x_i >= 0)
      u_i - 1.0 + x_i
     */

    for i in 1..=8 {
        IJth!(16 + i, i, 1.0);
        IJth!(16 + i, 16 + i, 1.0);
    }

    0
}

/*
 * Print solution
 */

fn PrintOutput(y: &NVector) {
    println!("     l=x+1          x         u=1-x");
    println!("   ----------------------------------");

    for i in 1..=NVAR {
        println!(
            " {}   {}   {}",
            fmt_g(Ith!(y, i + NVAR), 10, 6),
            fmt_g(Ith!(y, i), 10, 6),
            fmt_g(Ith!(y, i + 2 * NVAR), 10, 6)
        );
    }
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
