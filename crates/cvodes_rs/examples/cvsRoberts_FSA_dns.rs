/* -----------------------------------------------------------------
 * Translation of examples/cvodes/serial/cvsRoberts_FSA_dns.c
 * (SUNDIALS 7.7.0).
 *
 * Example problem (chemical kinetics with sensitivity analysis):
 *
 * The following is a simple example problem, with the coding
 * needed for its solution by CVODES for Forward Sensitivity
 * Analysis. The problem is from chemical kinetics, and consists
 * of the following three rate equations:
 *    dy1/dt = -p1*y1 + p2*y2*y3
 *    dy2/dt =  p1*y1 - p2*y2*y3 - p3*(y2)^2
 *    dy3/dt =  p3*(y2)^2
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions y1 = 1.0, y2 = y3 = 0. The reaction rates are:
 * p1=0.04, p2=1e4, and p3=3e7. The problem is stiff.
 * This program solves the problem with the BDF method, Newton
 * iteration with the dense linear solver, and a user-supplied
 * Jacobian routine. It uses a scalar relative tolerance and a
 * vector absolute tolerance (supplied through the user ewt fn).
 * Solution sensitivities with respect to p1/p2/p3 use the
 * user-supplied fS routine (of type SensRhs1Fn).
 *
 * Execution:
 *    % cvsRoberts_FSA_dns -nosensi
 *    % cvsRoberts_FSA_dns -sensi sensi_meth err_con
 * where sensi_meth is one of {sim, stg, stg1} and err_con is one
 * of {t, f}.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]

use cvodes_rs::cvodes::{
    CVode, CVodeCreate, CVodeFree, CVodeInit, CVodeSensEEtolerances, CVodeSensInit1,
    CVodeGetSens, CVodeWFtolerances,
};
use cvodes_rs::cvodes_io::{
    CVodeGetLastOrder, CVodeGetLastStep, CVodeGetNumSteps, CVodePrintAllStats,
    CVodeSetSensErrCon, CVodeSetSensParams, CVodeSetUserData,
};
use cvodes_rs::cvodes_ls::{CVodeSetJacFn, CVodeSetLinearSolver};
use cvodes_rs::sundials_utils::fmt_e;
use cvodes_rs::*;

/* Problem Constants */
const NEQ: i64 = 3; /* number of equations  */
const Y1: f64 = 1.0; /* initial y components */
const Y2: f64 = 0.0;
const Y3: f64 = 0.0;
const RTOL: f64 = 1.0e-4; /* scalar relative tolerance            */
const ATOL1: f64 = 1.0e-8; /* vector absolute tolerance components */
const ATOL2: f64 = 1.0e-14;
const ATOL3: f64 = 1.0e-6;
const T0: f64 = 0.0; /* initial time           */
const T1: f64 = 0.4; /* first output time      */
const TMULT: f64 = 10.0; /* output time factor     */
const NOUT: i32 = 12; /* number of output times */

const NS: usize = 3; /* number of sensitivities computed */

const ZERO: f64 = 0.0;

/* Type : UserData */
struct RobData {
    p: [f64; 3], /* problem parameters */
}

/*
 * f routine. Compute function f(t,y).
 */
fn f(_t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let data = user_data.as_ref().unwrap().downcast_ref::<RobData>().unwrap();
    let y1 = y.data[0];
    let y2 = y.data[1];
    let y3 = y.data[2];
    let p1 = data.p[0];
    let p2 = data.p[1];
    let p3 = data.p[2];

    let yd1 = -p1 * y1 + p2 * y2 * y3;
    ydot.data[0] = yd1;
    let yd3 = p3 * y2 * y2;
    ydot.data[2] = yd3;
    ydot.data[1] = -yd1 - yd3;

    0
}

/*
 * Jacobian routine. Compute J(t,y) = df/dy.
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
    let data = user_data.as_ref().unwrap().downcast_ref::<RobData>().unwrap();
    let y2 = y.data[1];
    let y3 = y.data[2];
    let p1 = data.p[0];
    let p2 = data.p[1];
    let p3 = data.p[2];

    if let SUNMatrix::Dense(dm) = j {
        /* IJth(J,i,j) = column-major dense storage */
        let m = 3usize;
        dm.data[0] = -p1; /* (1,1) */
        dm.data[m] = p2 * y3; /* (1,2) */
        dm.data[2 * m] = p2 * y2; /* (1,3) */
        dm.data[1] = p1; /* (2,1) */
        dm.data[1 + m] = -p2 * y3 - 2.0 * p3 * y2; /* (2,2) */
        dm.data[1 + 2 * m] = -p2 * y2; /* (2,3) */
        dm.data[2 + m] = 2.0 * p3 * y2; /* (3,2) */
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
    let data = user_data.as_ref().unwrap().downcast_ref::<RobData>().unwrap();
    let p1 = data.p[0];
    let p2 = data.p[1];
    let p3 = data.p[2];

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

/*
 * EwtSet function. Computes the error weights at the current solution.
 */
fn ewt(y: &NVector, w: &mut NVector, _user_data: &mut UserData) -> i32 {
    let rtol = RTOL;
    let atol = [ATOL1, ATOL2, ATOL3];

    for i in 0..3usize {
        let yy = y.data[i];
        let ww = rtol * yy.abs() + atol[i];
        if ww <= 0.0 {
            return -1;
        }
        w.data[i] = 1.0 / ww;
    }

    0
}

/*
 * Process and verify arguments to cvsfwddenx.
 */
fn process_args(argv: &[String]) -> (bool, i32, bool) {
    let sensi;
    let mut sensi_meth: i32 = -1;
    let mut err_con = false;

    if argv.len() < 2 {
        wrong_args(&argv[0]);
    }

    if argv[1] == "-nosensi" {
        sensi = false;
    } else if argv[1] == "-sensi" {
        sensi = true;
    } else {
        wrong_args(&argv[0]);
    }

    if sensi {
        if argv.len() != 4 {
            wrong_args(&argv[0]);
        }

        if argv[2] == "sim" {
            sensi_meth = CV_SIMULTANEOUS;
        } else if argv[2] == "stg" {
            sensi_meth = CV_STAGGERED;
        } else if argv[2] == "stg1" {
            sensi_meth = CV_STAGGERED1;
        } else {
            wrong_args(&argv[0]);
        }

        if argv[3] == "t" {
            err_con = true;
        } else if argv[3] == "f" {
            err_con = false;
        } else {
            wrong_args(&argv[0]);
        }
    }

    (sensi, sensi_meth, err_con)
}

fn wrong_args(name: &str) -> ! {
    println!("\nUsage: {} [-nosensi] [-sensi sensi_meth err_con]", name);
    println!("         sensi_meth = sim, stg, or stg1");
    println!("         err_con    = t or f");
    std::process::exit(0);
}

/*
 * Print current t, step count, order, stepsize, and solution.
 */
fn print_output(cvode_mem: &mut CVodeMem, t: f64, u: &NVector) {
    let mut nst: i64 = 0;
    let mut qu: i32 = 0;
    let mut hu: f64 = 0.0;

    CVodeGetNumSteps(cvode_mem, &mut nst);
    CVodeGetLastOrder(cvode_mem, &mut qu);
    CVodeGetLastStep(cvode_mem, &mut hu);

    println!("{} {:2}  {} {:5}", fmt_e(t, 8, 3), qu, fmt_e(hu, 8, 3), nst);

    print!("                  Solution       ");
    println!(
        "{} {} {} ",
        fmt_e(u.data[0], 12, 4),
        fmt_e(u.data[1], 12, 4),
        fmt_e(u.data[2], 12, 4)
    );
}

/*
 * Print sensitivities.
 */
fn print_output_s(uS: &[NVector]) {
    print!("                  Sensitivity 1  ");
    println!(
        "{} {} {} ",
        fmt_e(uS[0].data[0], 12, 4),
        fmt_e(uS[0].data[1], 12, 4),
        fmt_e(uS[0].data[2], 12, 4)
    );

    print!("                  Sensitivity 2  ");
    println!(
        "{} {} {} ",
        fmt_e(uS[1].data[0], 12, 4),
        fmt_e(uS[1].data[1], 12, 4),
        fmt_e(uS[1].data[2], 12, 4)
    );

    print!("                  Sensitivity 3  ");
    println!(
        "{} {} {} ",
        fmt_e(uS[2].data[0], 12, 4),
        fmt_e(uS[2].data[1], 12, 4),
        fmt_e(uS[2].data[2], 12, 4)
    );
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
 *-------------------------------
 * Main Program
 *-------------------------------
 */
fn main() {
    /* Process arguments */
    let argv: Vec<String> = std::env::args().collect();
    let (sensi, sensi_meth, err_con) = process_args(&argv);

    /* User data structure: initialize sensitivity variables (reaction
       rates for this problem) */
    let data = RobData {
        p: [0.04, 1.0e4, 3.0e7],
    };
    let p = data.p;

    /* Create the SUNDIALS context that all SUNDIALS objects require */
    let sunctx = SUNContext_Create();

    /* Initial conditions */
    let mut y = N_VNew_Serial(NEQ, &sunctx);
    y.data[0] = Y1;
    y.data[1] = Y2;
    y.data[2] = Y3;

    /* Call CVodeCreate to create the solver memory and specify the
     * Backward Differentiation Formula */
    let mut cvode_mem = CVodeCreate(CV_BDF, &sunctx);

    /* Call CVodeInit to initialize the integrator memory */
    let retval = CVodeInit(&mut cvode_mem, f, T0, &y);
    if check_retval(retval, "CVodeInit") {
        return;
    }

    /* Call CVodeWFtolerances to specify a user-supplied function ewt that
     * sets the multiplicative error weights w_i */
    let retval = CVodeWFtolerances(&mut cvode_mem, ewt);
    if check_retval(retval, "CVodeWFtolerances") {
        return;
    }

    /* Attach user data */
    let retval = CVodeSetUserData(&mut cvode_mem, Some(Box::new(data)));
    if check_retval(retval, "CVodeSetUserData") {
        return;
    }

    /* Create dense SUNMatrix and dense SUNLinearSolver object */
    let a_mat = SUNDenseMatrix(NEQ, NEQ, &sunctx);
    let ls = SUNLinSol_Dense(&y, &a_mat, &sunctx);

    /* Attach the matrix and linear solver */
    let retval = CVodeSetLinearSolver(&mut cvode_mem, ls, Some(a_mat));
    if check_retval(retval, "CVodeSetLinearSolver") {
        return;
    }

    /* Set the user-supplied Jacobian routine Jac */
    let retval = CVodeSetJacFn(&mut cvode_mem, Some(jac));
    if check_retval(retval, "CVodeSetJacFn") {
        return;
    }

    println!(" \n3-species kinetics problem");

    /* Sensitivity-related settings */
    let mut yS: Vec<NVector> = Vec::new();
    if sensi {
        /* Set parameter scaling factor */
        let pbar = [p[0], p[1], p[2]];

        /* Set sensitivity initial conditions */
        yS = (0..NS).map(|_| N_VClone(&y)).collect();
        for ys in yS.iter_mut() {
            N_VConst(ZERO, ys);
        }

        /* Call CVodeSensInit1 to activate forward sensitivity computations */
        let retval = CVodeSensInit1(&mut cvode_mem, NS as i32, sensi_meth, Some(fS), &yS);
        if check_retval(retval, "CVodeSensInit") {
            return;
        }

        /* Call CVodeSensEEtolerances to estimate tolerances for sensitivity
         * variables based on the tolerances supplied for states variables
         * and the scaling factor pbar */
        let retval = CVodeSensEEtolerances(&mut cvode_mem);
        if check_retval(retval, "CVodeSensEEtolerances") {
            return;
        }

        /* Set sensitivity analysis optional inputs */
        let retval = CVodeSetSensErrCon(&mut cvode_mem, err_con);
        if check_retval(retval, "CVodeSetSensErrCon") {
            return;
        }

        /* Call CVodeSetSensParams to specify problem parameter information
         * for sensitivity calculations */
        let retval = CVodeSetSensParams(&mut cvode_mem, None, Some(&pbar), None);
        if check_retval(retval, "CVodeSetSensParams") {
            return;
        }

        print!("Sensitivity: YES ");
        if sensi_meth == CV_SIMULTANEOUS {
            print!("( SIMULTANEOUS +");
        } else if sensi_meth == CV_STAGGERED {
            print!("( STAGGERED +");
        } else {
            print!("( STAGGERED1 +");
        }
        if err_con {
            print!(" FULL ERROR CONTROL )");
        } else {
            print!(" PARTIAL ERROR CONTROL )");
        }
    } else {
        print!("Sensitivity: NO ");
    }

    /* In loop, call CVode, print results, and test for error. */
    println!("\n");
    print!("===========================================");
    println!("============================");
    print!("     T     Q       H      NST           y1");
    println!("           y2           y3    ");
    print!("===========================================");
    println!("============================");

    let mut t = T0;
    let mut tout = T1;
    for _iout in 1..=NOUT {
        let retval = CVode(&mut cvode_mem, tout, &mut y, &mut t, CV_NORMAL);
        if check_retval(retval, "CVode") {
            break;
        }

        print_output(&mut cvode_mem, t, &y);

        /* Call CVodeGetSens to get the sensitivity solution vector after a
         * successful return from CVode */
        if sensi {
            let retval = CVodeGetSens(&cvode_mem, &mut t, &mut yS);
            if check_retval(retval, "CVodeGetSens") {
                break;
            }
            print_output_s(&yS);
        }
        print!("-----------------------------------------");
        println!("------------------------------");

        tout *= TMULT;
    }

    /* Print final statistics to the screen */
    println!("\nFinal Statistics:");
    let mut stdout = std::io::stdout();
    let _ = CVodePrintAllStats(&mut cvode_mem, &mut stdout, SUN_OUTPUTFORMAT_TABLE);

    /* Print final statistics to a file in CSV format */
    let mut fname = String::from("cvsRoberts_FSA_dns_stats");
    if sensi {
        if sensi_meth == CV_SIMULTANEOUS {
            fname.push_str("_-sensi_sim");
        } else if sensi_meth == CV_STAGGERED {
            fname.push_str("_-sensi_stg");
        } else {
            fname.push_str("_-sensi_stg1");
        }
        if err_con {
            fname.push_str("_t");
        } else {
            fname.push_str("_f");
        }
    }
    fname.push_str(".csv");
    let mut fid = std::fs::File::create(&fname).expect("fopen");
    let _ = CVodePrintAllStats(&mut cvode_mem, &mut fid, SUN_OUTPUTFORMAT_CSV);
    drop(fid);

    /* Free memory */
    drop(y);
    drop(yS);
    CVodeFree(cvode_mem);
}
