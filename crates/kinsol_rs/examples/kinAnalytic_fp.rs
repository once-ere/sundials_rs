/* -----------------------------------------------------------------------------
 * Translated from examples/kinsol/serial/kinAnalytic_fp.c (KINSOL 7.7.0)
 *
 * This example solves the nonlinear system
 *
 * 3x - cos((y-1)z) - 1/2 = 0
 * x^2 - 81(y-0.9)^2 + sin(z) + 1.06 = 0
 * exp(-x(y-1)) + 20z + (10 pi - 3)/3 = 0
 *
 * using the accelerated fixed pointer solver in KINSOL. The nonlinear fixed
 * point function is
 *
 * g1(x,y,z) = 1/3 cos((y-1)yz) + 1/6
 * g2(x,y,z) = 1/9 sqrt(x^2 + sin(z) + 1.06) + 0.9
 * g3(x,y,z) = -1/20 exp(-x(y-1)) - (10 pi - 3) / 60
 *
 * This system has the analytic solution x = 1/2, y = 1, z = -pi/6.
 * ---------------------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use kinsol_rs::sundials_utils::fmt_g;
use kinsol_rs::*;

/* problem constants */
const NEQ: i64 = 3; /* number of equations */

const ZERO: f64 = 0.0; /* real 0.0  */
const PTONE: f64 = 0.1; /* real 0.1  */
const HALF: f64 = 0.5; /* real 0.5  */
const PTNINE: f64 = 0.9; /* real 0.9  */
const ONE: f64 = 1.0; /* real 1.0  */
const ONEPTZEROSIX: f64 = 1.06; /* real 1.06 */
const THREE: f64 = 3.0; /* real 3.0  */
const SIX: f64 = 6.0; /* real 6.0  */
const NINE: f64 = 9.0; /* real 9.0  */
const TEN: f64 = 10.0; /* real 10.0 */
const TWENTY: f64 = 20.0; /* real 20.0 */
const SIXTY: f64 = 60.0; /* real 60.0 */
const PI: f64 = 3.1415926535898; /* real pi   */

/* analytic solution */
const XTRUE: f64 = HALF;
const YTRUE: f64 = ONE;
const ZTRUE: f64 = -PI / SIX;

/* problem options */
struct UserOpt {
    tol: f64,             /* solve tolerance                  */
    maxiter: i64,         /* max number of iterations         */
    m_aa: i64,            /* number of acceleration vectors   */
    delay_aa: i64,        /* number of iterations to delay AA */
    orth_aa: i32,         /* orthogonalization method         */
    damping_fp: f64,      /* damping parameter for FP         */
    damping_aa: f64,      /* damping parameter for AA         */
    use_damping_fn: bool, /* damping function                 */
    use_depth_fn: bool,   /* depth function                   */
}

/* C atoi semantics: leading (optionally signed) integer, 0 on failure */
fn atoi(s: &str) -> i64 {
    let t = s.trim_start();
    let bytes = t.as_bytes();
    let mut i = 0;
    if i < bytes.len() && (bytes[i] == b'+' || bytes[i] == b'-') {
        i += 1;
    }
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        i += 1;
    }
    t[..i].parse().unwrap_or(0)
}

/* -----------------------------------------------------------------------------
 * Main program
 * ---------------------------------------------------------------------------*/
fn main() {
    /* Set default options */
    let mut uopt = SetDefaults();

    let args: Vec<String> = std::env::args().collect();
    let mut retval = ReadInputs(&args, &mut uopt);
    if check_retval(retval, "ReadInputs") {
        std::process::exit(1);
    }

    /* -------------------------
     * Print problem description
     * ------------------------- */

    println!("Solve the nonlinear system:");
    println!("    3x - cos((y-1)z) - 1/2 = 0");
    println!("    x^2 - 81(y-0.9)^2 + sin(z) + 1.06 = 0");
    println!("    exp(-x(y-1)) + 20z + (10 pi - 3)/3 = 0");
    println!("Analytic solution:");
    println!("    x = {}", fmt_g(XTRUE, 0, 6));
    println!("    y = {}", fmt_g(YTRUE, 0, 6));
    println!("    z = {}", fmt_g(ZTRUE, 0, 6));
    println!("Solution method: Anderson accelerated fixed point iteration.");
    println!("    tolerance    = {}", fmt_g(uopt.tol, 0, 6));
    println!("    max iters    = {}", uopt.maxiter);
    println!("    m_aa         = {}", uopt.m_aa);
    println!("    delay_aa     = {}", uopt.delay_aa);
    println!("    damping_aa   = {}", fmt_g(uopt.damping_aa, 0, 6));
    println!("    damping_fp   = {}", fmt_g(uopt.damping_fp, 0, 6));
    if uopt.use_damping_fn {
        println!("    damping_fn   = ON");
    } else {
        println!("    damping_fn   = OFF");
    }
    if uopt.use_depth_fn {
        println!("    depth_fn     = ON");
    } else {
        println!("    depth_fn     = OFF");
    }
    println!("    orth routine = {}", uopt.orth_aa);

    /* Create the SUNDIALS context that all SUNDIALS objects require */
    let sunctx = SUNContext_Create();

    /* --------------------------------------
     * Create vectors for solution and scales
     * -------------------------------------- */

    let mut u = N_VNew_Serial(NEQ, &sunctx);

    let mut scale = N_VClone(&u);

    /* -----------------------------------------
     * Initialize and allocate memory for KINSOL
     * ----------------------------------------- */

    let mut kmem = KINCreate(&sunctx);

    /* Set number of prior residuals used in Anderson acceleration */
    let _ = KINSetMAA(&mut kmem, uopt.m_aa);

    /* Set orthogonalization routine used in Anderson acceleration */
    retval = KINSetOrthAA(&mut kmem, uopt.orth_aa);
    if check_retval(retval, "KINSetOrthAA") {
        std::process::exit(1);
    }

    retval = KINInit(&mut kmem, FPFunction, &u);
    if check_retval(retval, "KINInit") {
        std::process::exit(1);
    }

    /* -------------------
     * Set optional inputs
     * ------------------- */

    /* Specify stopping tolerance based on residual */
    retval = KINSetFuncNormTol(&mut kmem, uopt.tol);
    if check_retval(retval, "KINSetFuncNormTol") {
        std::process::exit(1);
    }

    /* Set maximum number of iterations */
    retval = KINSetNumMaxIters(&mut kmem, uopt.maxiter);
    if check_retval(retval, "KINSetNumMaxItersFuncNormTol") {
        std::process::exit(1);
    }

    /* Set Fixed point damping parameter */
    if uopt.m_aa == 0 {
        let _ = KINSetDamping(&mut kmem, uopt.damping_fp);
    }

    /* Set Anderson acceleration options */
    if uopt.m_aa > 0 {
        /* Set damping parameter */
        retval = KINSetDampingAA(&mut kmem, uopt.damping_aa);
        if check_retval(retval, "KINSetDampingAA") {
            std::process::exit(1);
        }

        /* Set acceleration delay */
        retval = KINSetDelayAA(&mut kmem, uopt.delay_aa);
        if check_retval(retval, "KINSetDelayAA") {
            std::process::exit(1);
        }
    }

    if uopt.use_damping_fn {
        /* Attach user defined damping function */
        retval = KINSetDampingFn(&mut kmem, Some(DampingFn));
        if check_retval(retval, "KINSetDampingFn") {
            std::process::exit(1);
        }
    }

    if uopt.use_depth_fn {
        /* Attach user defined depth function */
        retval = KINSetDepthFn(&mut kmem, Some(DepthFn));
        if check_retval(retval, "KINSetDepthFn") {
            std::process::exit(1);
        }
    }

    /* -------------
     * Initial guess
     * ------------- */

    /* Get vector data array */
    let data = N_VGetArrayPointer(&mut u);

    data[0] = PTONE;
    data[1] = PTONE;
    data[2] = -PTONE;

    /* ----------------------------
     * Call KINSol to solve problem
     * ---------------------------- */

    /* No scaling used */
    N_VConst(ONE, &mut scale);

    /* Call main solver */
    retval = KINSol(&mut kmem, /* KINSol memory block */
                    &mut u,    /* initial guess on input; solution vector */
                    KIN_FP,    /* global strategy choice */
                    &scale,    /* scaling vector, for the variable cc */
                    &scale);   /* scaling vector for function values fval */
    if check_retval(retval, "KINSol") {
        std::process::exit(1);
    }

    /* ------------------------------------
     * Get solver statistics
     * ------------------------------------ */

    /* get solver stats */
    let mut nni: i64 = 0;
    let mut nfe: i64 = 0;
    retval = KINGetNumNonlinSolvIters(&mut kmem, &mut nni);
    check_retval(retval, "KINGetNumNonlinSolvIters");

    retval = KINGetNumFuncEvals(&mut kmem, &mut nfe);
    check_retval(retval, "KINGetNumFuncEvals");

    println!("\nFinal Statistics:");
    println!("Number of nonlinear iterations: {:6}", nni);
    println!("Number of function evaluations: {:6}", nfe);

    /* ------------------------------------
     * Print solution and check error
     * ------------------------------------ */

    /* check solution */
    let retval = check_ans(&u, uopt.tol);

    /* -----------
     * Free memory
     * ----------- */

    KINFree(kmem);

    std::process::exit(retval);
}

/* -----------------------------------------------------------------------------
 * Nonlinear system
 *
 * 3x - cos((y-1)z) - 1/2 = 0
 * x^2 - 81(y-0.9)^2 + sin(z) + 1.06 = 0
 * exp(-x(y-1)) + 20z + (10 pi - 3)/3 = 0
 *
 * Nonlinear fixed point function
 *
 * g1(x,y,z) = 1/3 cos((y-1)z) + 1/6
 * g2(x,y,z) = 1/9 sqrt(x^2 + sin(z) + 1.06) + 0.9
 * g3(x,y,z) = -1/20 exp(-x(y-1)) - (10 pi - 3) / 60
 *
 * ---------------------------------------------------------------------------*/
fn FPFunction(u: &NVector, g: &mut NVector, _user_data: &mut UserData) -> i32 {
    /* Get vector data arrays */
    let udata = &u.data;
    let gdata = &mut g.data;

    let x = udata[0];
    let y = udata[1];
    let z = udata[2];

    gdata[0] = (ONE / THREE) * ((y - ONE) * z).cos() + (ONE / SIX);
    gdata[1] = (ONE / NINE) * (x * x + z.sin() + ONEPTZEROSIX).sqrt() + PTNINE;
    gdata[2] = -(ONE / TWENTY) * (-x * (y - ONE)).exp() - (TEN * PI - THREE) / SIXTY;

    0
}

fn DampingFn(
    _iter: i64,
    u_val: &NVector,
    g_val: &NVector,
    qt_fn: &[f64],
    depth: i64,
    _user_data: &mut UserData,
    damping_factor: &mut f64,
) -> i32 {
    if depth == 0 {
        *damping_factor = 0.5;
    } else {
        /* Compute ||Q^T fn||^2 */
        let mut qt_fn_norm_sqr = ZERO;
        for i in 0..depth as usize {
            qt_fn_norm_sqr += qt_fn[i] * qt_fn[i];
        }

        /* Compute ||fn||^2 = ||G(u_n) - u_n||^2 */
        let g_data = &g_val.data;
        let u_data = &u_val.data;
        let mut f_n = [0.0; 3];
        for i in 0..3 {
            f_n[i] = g_data[i] - u_data[i];
        }
        let mut fn_norm_sqr = ZERO;
        for i in 0..3 {
            fn_norm_sqr += f_n[i] * f_n[i];
        }

        /* Compute the gain = sqrt(1 - ||Q^T fn||^2 / ||fn||^2) */
        let gain = SUNRsqrt(ONE - qt_fn_norm_sqr / fn_norm_sqr);

        *damping_factor = 0.9 - 0.5 * gain;
    }

    0
}

fn DepthFn(
    iter: i64,
    _u_val: &NVector,
    _g_val: &NVector,
    _f_val: &NVector,
    _df: &[NVector],
    _r_mat: &[f64],
    depth: i64,
    _user_data: &mut UserData,
    new_depth: &mut i64,
    _remove_indices: &mut [bool],
) -> i32 {
    if iter < 2 {
        *new_depth = 1;
    } else {
        *new_depth = depth;
    }

    0
}

/* -----------------------------------------------------------------------------
 * Check the solution of the nonlinear system and return PASS or FAIL
 * ---------------------------------------------------------------------------*/
fn check_ans(u: &NVector, mut tol: f64) -> i32 {
    /* Get vector data array */
    let data = &u.data;

    /* print the solution */
    println!("Computed solution:");
    println!("    x = {}", fmt_g(data[0], 0, 6));
    println!("    y = {}", fmt_g(data[1], 0, 6));
    println!("    z = {}", fmt_g(data[2], 0, 6));

    /* solution error */
    let ex = (data[0] - XTRUE).abs();
    let ey = (data[1] - YTRUE).abs();
    let ez = (data[2] - ZTRUE).abs();

    /* print the solution error */
    println!("Solution error:");
    println!("    ex = {}", fmt_g(ex, 0, 6));
    println!("    ey = {}", fmt_g(ey, 0, 6));
    println!("    ez = {}", fmt_g(ez, 0, 6));

    tol *= TEN;
    if ex > tol || ey > tol || ez > tol {
        println!("FAIL");
        return 1;
    }

    println!("PASS");
    0
}

/* -----------------------------------------------------------------------------
 * Set default options
 * ---------------------------------------------------------------------------*/
fn SetDefaults() -> UserOpt {
    UserOpt {
        tol: 100.0 * SUNRsqrt(SUN_UNIT_ROUNDOFF),
        maxiter: 30,
        m_aa: 0,               /* no acceleration */
        delay_aa: 0,           /* no delay        */
        orth_aa: 0,            /* MGS             */
        damping_fp: 1.0,       /* no FP dampig    */
        damping_aa: 1.0,       /* no AA damping   */
        use_damping_fn: false, /* no damping fn   */
        use_depth_fn: false,   /* no depth fn     */
    }
}

/* -----------------------------------------------------------------------------
 * Read command line inputs
 * ---------------------------------------------------------------------------*/
fn ReadInputs(argv: &[String], uopt: &mut UserOpt) -> i32 {
    let mut arg_index = 1;

    let arg_at = |i: usize| -> &str { argv.get(i).map(String::as_str).unwrap_or("") };

    while arg_index < argv.len() {
        if argv[arg_index] == "--tol" {
            arg_index += 1;
            uopt.tol = SUNStrToReal(arg_at(arg_index));
            arg_index += 1;
        } else if argv[arg_index] == "--maxiter" {
            arg_index += 1;
            uopt.maxiter = atoi(arg_at(arg_index));
            arg_index += 1;
        } else if argv[arg_index] == "--m_aa" {
            arg_index += 1;
            uopt.m_aa = atoi(arg_at(arg_index));
            arg_index += 1;
        } else if argv[arg_index] == "--delay_aa" {
            arg_index += 1;
            uopt.delay_aa = atoi(arg_at(arg_index));
            arg_index += 1;
        } else if argv[arg_index] == "--damping_fp" {
            arg_index += 1;
            uopt.damping_fp = SUNStrToReal(arg_at(arg_index));
            arg_index += 1;
        } else if argv[arg_index] == "--damping_aa" {
            arg_index += 1;
            uopt.damping_aa = SUNStrToReal(arg_at(arg_index));
            arg_index += 1;
        } else if argv[arg_index] == "--damping_fn" {
            arg_index += 1;
            uopt.use_damping_fn = true;
        } else if argv[arg_index] == "--depth_fn" {
            arg_index += 1;
            uopt.use_depth_fn = true;
        } else if argv[arg_index] == "--orth_aa" {
            arg_index += 1;
            uopt.orth_aa = atoi(arg_at(arg_index)) as i32;
            arg_index += 1;
        } else if argv[arg_index] == "--help" {
            InputHelp();
            return -1;
        } else {
            println!("Error: Invalid command line parameter {}", argv[arg_index]);
            InputHelp();
            return -1;
        }
    }

    0
}

/* -----------------------------------------------------------------------------
 * Print command line options
 * ---------------------------------------------------------------------------*/
fn InputHelp() {
    println!();
    println!(" Command line options:");
    println!("   --tol        : nonlinear solver tolerance");
    println!("   --maxiter    : max number of nonlinear iterations");
    println!("   --m_aa       : number of Anderson acceleration vectors");
    println!("   --delay_aa   : Anderson acceleration delay");
    println!("   --damping_fp : fixed point damping parameter");
    println!("   --damping_aa : Anderson acceleration damping parameter");
    println!("   --orth_aa    : Anderson acceleration orthogonalization method");
    println!("   --damping_fn : user defined damping function");
    println!("   --depth_fn   : user defined depth function");
}

/* -----------------------------------------------------------------------------
 * Check function return value (C opt == 1 case: non-zero value is an error;
 * the opt == 0 NULL-pointer case does not arise in the Rust port)
 * ---------------------------------------------------------------------------*/
fn check_retval(retval: i32, funcname: &str) -> bool {
    if retval != 0 {
        eprintln!("\nERROR: {}() failed -- returned {}\n", funcname, retval);
        return true;
    }
    false
}
