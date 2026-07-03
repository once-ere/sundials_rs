/* -----------------------------------------------------------------
 * Translated from examples/kinsol/serial/kinLaplace_bnd.c (KINSOL 7.7.0)
 *
 * This example solves a 2D elliptic PDE
 *
 *    d^2 u / dx^2 + d^2 u / dy^2 = u^3 - u - 2.0
 *
 * subject to homogeneous Dirichlet boundary conditions.
 * The PDE is discretized on a uniform NX+2 by NY+2 grid with
 * central differencing, and with boundary values eliminated,
 * leaving a system of size NEQ = NX*NY.
 * The nonlinear system is solved by KINSOL using the SUNBAND linear
 * solver.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use kinsol_rs::sundials_utils::{fmt_f, fmt_g};
use kinsol_rs::*;

/* Problem Constants */

const NX: i64 = 31; /* no. of points in x direction */
const NY: i64 = 31; /* no. of points in y direction */
const NEQ: i64 = NX * NY; /* problem dimension */

const SKIP: i64 = 3; /* no. of points skipped for printing */

const FTOL: f64 = 1.0e-12; /* function tolerance */

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;

/* IJth is defined in order to isolate the translation from the
   mathematical 2-dimensional structure of the dependent variable vector
   to the underlying 1-dimensional storage.
   IJth(vdata,i,j) references the element in the vdata array for
   u at mesh point (i,j), where 1 <= i <= NX, 1 <= j <= NY.
   The variables are ordered by the y index j, then by the x index i. */

macro_rules! IJth {
    ($vdata:expr, $i:expr, $j:expr) => {
        $vdata[(($j - 1) + ($i - 1) * NY) as usize]
    };
}

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

fn main() {
    /* -------------------------
     * Print problem description
     * ------------------------- */

    println!("\n2D elliptic PDE on unit square");
    println!("   d^2 u / dx^2 + d^2 u / dy^2 = u^3 - u + 2.0");
    print!(" + homogeneous Dirichlet boundary conditions\n\n");
    println!("Solution method: Modified Newton with band linear solver");
    println!("Problem size: {:2} x {:2} = {:4}", NX, NY, NEQ);

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

    let mut retval = KINInit(&mut kmem, func, &y);
    if check_retval(retval, "KINInit") {
        std::process::exit(1);
    }

    /* -------------------
     * Set optional inputs
     * ------------------- */

    /* Specify stopping tolerance based on residual */

    let fnormtol = FTOL;
    retval = KINSetFuncNormTol(&mut kmem, fnormtol);
    if check_retval(retval, "KINSetFuncNormTol") {
        std::process::exit(1);
    }

    /* -------------------------
     * Create band SUNMatrix
     * ------------------------- */

    let J = SUNBandMatrix(NEQ, NX, NX, &sunctx);

    /* ---------------------------
     * Create band SUNLinearSolver
     * --------------------------- */

    let LS = SUNLinSol_Band(&y, &J, &sunctx);

    /* -------------------------
     * Attach band linear solver
     * ------------------------- */

    retval = KINSetLinearSolver(&mut kmem, LS, Some(J));
    if check_retval(retval, "KINSetLinearSolver") {
        std::process::exit(1);
    }

    /* ------------------------------
     * Parameters for Modified Newton
     * ------------------------------ */

    /* Force a Jacobian re-evaluation every mset iterations */
    let mset = 100;
    retval = KINSetMaxSetupCalls(&mut kmem, mset);
    if check_retval(retval, "KINSetMaxSetupCalls") {
        std::process::exit(1);
    }

    /* Every msubset iterations, test if a Jacobian evaluation
    is necessary */
    let msubset = 1;
    retval = KINSetMaxSubSetupCalls(&mut kmem, msubset);
    if check_retval(retval, "KINSetMaxSubSetupCalls") {
        std::process::exit(1);
    }

    /* -------------
     * Initial guess
     * ------------- */

    N_VConst(ZERO, &mut y);

    /* ----------------------------
     * Call KINSol to solve problem
     * ---------------------------- */

    /* No scaling used */
    N_VConst(ONE, &mut scale);

    /* Call main solver */
    retval = KINSol(
        &mut kmem,      /* KINSol memory block */
        &mut y,         /* initial guess on input; solution vector */
        KIN_LINESEARCH, /* global strategy choice */
        &scale,         /* scaling vector, for the variable cc */
        &scale,         /* scaling vector for function values fval */
    );
    if check_retval(retval, "KINSol") {
        std::process::exit(1);
    }

    /* ------------------------------------
     * Print solution and solver statistics
     * ------------------------------------ */

    /* Get scaled norm of the system function */

    let mut fnorm = ZERO;
    retval = KINGetFuncNorm(&mut kmem, &mut fnorm);
    if check_retval(retval, "KINGetfuncNorm") {
        std::process::exit(1);
    }

    print!("\nComputed solution (||F|| = {}):\n\n", fmt_g(fnorm, 0, 6));
    PrintOutput(&y);

    PrintFinalStats(&mut kmem);

    /* -----------
     * Free memory
     * ----------- */

    KINFree(kmem);
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * System function
 */

fn func(u: &NVector, f: &mut NVector, _user_data: &mut UserData) -> i32 {
    let dx = ONE / (NX + 1) as f64;
    let dy = ONE / (NY + 1) as f64;
    let hdc = ONE / (dx * dx);
    let vdc = ONE / (dy * dy);

    let udata = &u.data;
    let fdata = &mut f.data;

    for j in 1..=NY {
        for i in 1..=NX {
            /* Extract u at x_i, y_j and four neighboring points */

            let uij = IJth!(udata, i, j);
            let udn = if j == 1 { ZERO } else { IJth!(udata, i, j - 1) };
            let uup = if j == NY { ZERO } else { IJth!(udata, i, j + 1) };
            let ult = if i == 1 { ZERO } else { IJth!(udata, i - 1, j) };
            let urt = if i == NX { ZERO } else { IJth!(udata, i + 1, j) };

            /* Evaluate diffusion components */

            let hdiff = hdc * (ult - TWO * uij + urt);
            let vdiff = vdc * (uup - TWO * uij + udn);

            /* Set residual at x_i, y_j */

            IJth!(fdata, i, j) = hdiff + vdiff + uij - uij * uij * uij + 2.0;
        }
    }

    0
}

/*
 * Print solution at selected points
 */

fn PrintOutput(u: &NVector) {
    let dx = ONE / (NX + 1) as f64;
    let dy = ONE / (NY + 1) as f64;

    let udata = &u.data;

    print!("            ");
    let mut i = 1;
    while i <= NX {
        let x = i as f64 * dx;
        print!("{:<8} ", fmt_f(x, 0, 5)); /* %-8.5f */
        i += SKIP;
    }
    print!("\n\n");

    let mut j = 1;
    while j <= NY {
        let y = j as f64 * dy;
        print!("{:<8}    ", fmt_f(y, 0, 5)); /* %-8.5f */
        let mut i = 1;
        while i <= NX {
            print!("{:<8} ", fmt_f(IJth!(udata, i, j), 0, 5)); /* %-8.5f */
            i += SKIP;
        }
        println!();
        j += SKIP;
    }
}

/*
 * Print final statistics
 */

fn PrintFinalStats(kmem: &mut KINMem) {
    let (mut nni, mut nfe, mut nje, mut nfeD) = (0i64, 0i64, 0i64, 0i64);
    let (mut lenrw, mut leniw, mut lenrwB, mut leniwB) = (0i64, 0i64, 0i64, 0i64);
    let (mut nbcfails, mut nbacktr) = (0i64, 0i64);

    /* Main solver statistics */

    let mut retval = KINGetNumNonlinSolvIters(kmem, &mut nni);
    check_retval(retval, "KINGetNumNonlinSolvIters");
    retval = KINGetNumFuncEvals(kmem, &mut nfe);
    check_retval(retval, "KINGetNumFuncEvals");

    /* Linesearch statistics */

    retval = KINGetNumBetaCondFails(kmem, &mut nbcfails);
    check_retval(retval, "KINGetNumBetacondFails");
    retval = KINGetNumBacktrackOps(kmem, &mut nbacktr);
    check_retval(retval, "KINGetNumBacktrackOps");

    /* Main solver workspace size */

    retval = KINGetWorkSpace(kmem, &mut lenrw, &mut leniw);
    check_retval(retval, "KINGetWorkSpace");

    /* Band linear solver statistics */

    retval = KINGetNumJacEvals(kmem, &mut nje);
    check_retval(retval, "KINGetNumJacEvals");
    retval = KINGetNumLinFuncEvals(kmem, &mut nfeD);
    check_retval(retval, "KINGetNumLinFuncEvals");

    /* Band linear solver workspace size */

    retval = KINGetLinWorkSpace(kmem, &mut lenrwB, &mut leniwB);
    check_retval(retval, "KINGetLinWorkSpace");

    print!("\nFinal Statistics.. \n\n");
    println!("nni      = {:6}    nfe     = {:6} ", nni, nfe);
    println!("nbcfails = {:6}    nbacktr = {:6} ", nbcfails, nbacktr);
    println!("nje      = {:6}    nfeB    = {:6} ", nje, nfeD);
    println!();
    println!("lenrw    = {:6}    leniw   = {:6} ", lenrw, leniw);
    println!("lenrwB   = {:6}    leniwB  = {:6} ", lenrwB, leniwB);
}

/*
 * Check function return value (opt == 1 case of the C check_retval)
 */

fn check_retval(retval: i32, funcname: &str) -> bool {
    if retval < 0 {
        eprintln!(
            "\nSUNDIALS_ERROR: {}() failed with retval = {}\n",
            funcname, retval
        );
        return true;
    }
    false
}
