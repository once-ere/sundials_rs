/* -----------------------------------------------------------------
 * Translated from examples/kinsol/serial/kinLaplace_picard_kry.c
 * (KINSOL 7.7.0)
 *
 * This example solves a 2D elliptic PDE
 *
 *    d^2 u / dx^2 + d^2 u / dy^2 = u^3 - u - 2.0
 *
 * subject to homogeneous Dirichlet boundary conditions.
 * The PDE is discretized on a uniform NX+2 by NY+2 grid with
 * central differencing, and with boundary values eliminated,
 * leaving a system of size NEQ = NX*NY.
 * The nonlinear system is solved by KINSOL using the Picard
 * iteration and the SPGMR linear solver.
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

/* function tolerance */
const FTOL: f64 = 1.0e-12;

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
    println!("Solution method: Anderson accelerated Picard iteration with SPGMR linear solver.");
    println!("Problem size: {:2} x {:2} = {:4}", NX, NY, NEQ);

    /* Create the SUNDIALS context that all SUNDIALS objects require */
    let sunctx = SUNContext_Create();

    /* --------------------------------------
     * Create vectors for solution and scales
     * -------------------------------------- */

    let mut y = N_VNew_Serial(NEQ, &sunctx);
    let mut scale = N_VNew_Serial(NEQ, &sunctx);

    /* ----------------------------------------------------------------------------------
     * Initialize and allocate memory for KINSOL, set parameters for Anderson acceleration
     * ---------------------------------------------------------------------------------- */

    let mut kmem = KINCreate(&sunctx);

    /* y is used as a template */

    let mut retval = KINInit(&mut kmem, func, &y);
    if check_retval(retval, "KINInit") {
        std::process::exit(1);
    }

    /* -------------------
     * Set optional inputs
     * ------------------- */

    /* Use acceleration with up to 3 prior residuals */

    retval = KINSetMAA(&mut kmem, 3);
    if check_retval(retval, "KINSetMAA") {
        std::process::exit(1);
    }

    /* Specify stopping tolerance based on residual */

    let fnormtol = FTOL;
    retval = KINSetFuncNormTol(&mut kmem, fnormtol);
    if check_retval(retval, "KINSetFuncNormTol") {
        std::process::exit(1);
    }

    /* ----------------------
     * Create SUNLinearSolver
     * ---------------------- */

    let LS = SUNLinSol_SPGMR(&y, SUN_PREC_NONE, 10, &sunctx);

    /* --------------------
     * Attach linear solver
     * -------------------- */

    retval = KINSetLinearSolver(&mut kmem, LS, None);
    if check_retval(retval, "KINSetLinearSolver") {
        std::process::exit(1);
    }

    /* ------------------------------------
     * Set Jacobian vector product function
     * ------------------------------------ */

    retval = KINSetJacTimesVecFn(&mut kmem, Some(jactimes));
    if check_retval(retval, "KINSetJacTimesVecFn") {
        std::process::exit(1);
    }

    /* -------------
     * Initial guess
     * ------------- */

    N_VConst(ZERO, &mut y);
    IJth!(y.data, 2, 2) = ONE;

    /* ----------------------------
     * Call KINSol to solve problem
     * ---------------------------- */

    /* No scaling used */
    N_VConst(ONE, &mut scale);

    /* Call main solver */
    retval = KINSol(
        &mut kmem,  /* KINSol memory block */
        &mut y,     /* initial guess on input; solution vector */
        KIN_PICARD, /* global strategy choice */
        &scale,     /* scaling vector, for the variable cc */
        &scale,     /* scaling vector for function values fval */
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
 * Jacobian vector product function
 */

fn jactimes(
    v: &NVector,
    Jv: &mut NVector,
    _u: &NVector,
    _new_u: &mut bool,
    _user_data: &mut UserData,
) -> i32 {
    let dx = ONE / (NX + 1) as f64;
    let dy = ONE / (NY + 1) as f64;
    let hdc = ONE / (dx * dx);
    let vdc = ONE / (dy * dy);

    let vdata = &v.data;
    let Jvdata = &mut Jv.data;

    for j in 1..=NY {
        for i in 1..=NX {
            /* Extract v at x_i, y_j and four neighboring points */

            let vij = IJth!(vdata, i, j);
            let vdn = if j == 1 { ZERO } else { IJth!(vdata, i, j - 1) };
            let vup = if j == NY { ZERO } else { IJth!(vdata, i, j + 1) };
            let vlt = if i == 1 { ZERO } else { IJth!(vdata, i - 1, j) };
            let vrt = if i == NX { ZERO } else { IJth!(vdata, i + 1, j) };

            /* Evaluate diffusion components */

            let hdiff = hdc * (vlt - TWO * vij + vrt);
            let vdiff = vdc * (vup - TWO * vij + vdn);

            /* Set Jv at x_i, y_j */

            IJth!(Jvdata, i, j) = hdiff + vdiff;
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
    let (mut nni, mut nfe, mut nli, mut npe, mut nps, mut ncfl, mut nfeLS, mut njvevals) =
        (0i64, 0i64, 0i64, 0i64, 0i64, 0i64, 0i64, 0i64);

    /* Main solver statistics */

    let mut retval = KINGetNumNonlinSolvIters(kmem, &mut nni);
    check_retval(retval, "KINGetNumNonlinSolvIters");
    retval = KINGetNumFuncEvals(kmem, &mut nfe);
    check_retval(retval, "KINGetNumFuncEvals");

    /* Linear solver statistics */

    retval = KINGetNumLinIters(kmem, &mut nli);
    check_retval(retval, "KINGetNumLinIters");
    retval = KINGetNumLinFuncEvals(kmem, &mut nfeLS);
    check_retval(retval, "KINGetNumLinFuncEvals");
    retval = KINGetNumLinConvFails(kmem, &mut ncfl);
    check_retval(retval, "KINGetNumLinConvFails");
    retval = KINGetNumJtimesEvals(kmem, &mut njvevals);
    check_retval(retval, "KINGetNumJtimesEvals");
    retval = KINGetNumPrecEvals(kmem, &mut npe);
    check_retval(retval, "KINGetNumPrecEvals");
    retval = KINGetNumPrecSolves(kmem, &mut nps);
    check_retval(retval, "KINGetNumPrecSolves");

    print!("\nFinal Statistics.. \n\n");
    println!("nni = {:6}  nli   = {:6}  ncfl = {:6}", nni, nli, ncfl);
    println!("nfe = {:6}  nfeLS = {:6}  njt  = {:6}", nfe, nfeLS, njvevals);
    print!("npe = {:6}  nps   = {:6}", npe, nps);
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
