/* -----------------------------------------------------------------
 * Translation of examples/cvodes/serial/cvsAdvDiff_FSA_non.c
 * (SUNDIALS 7.7.0).
 *
 * Example problem:
 *
 * The following is a simple example problem, with the program for
 * its solution by CVODES. The problem is the semi-discrete form of
 * the advection-diffusion equation in 1-D:
 *   du/dt = q1 * d^2 u / dx^2 + q2 * du/dx
 * on the interval 0 <= x <= 2, and the time interval 0 <= t <= 5.
 * Homogeneous Dirichlet boundary conditions are posed, and the
 * initial condition is:
 *   u(x,y,t=0) = x(2-x)exp(2x).
 * The PDE is discretized on a uniform grid of size MX+2 with
 * central differencing, and with boundary values eliminated,
 * leaving an ODE system of size NEQ = MX.
 * This program solves the problem with the option for nonstiff
 * systems: ADAMS method and fixed-point iteration.
 * It uses scalar relative and absolute tolerances.
 * Output is printed at t = .5, 1.0, ..., 5.
 *
 * Optionally, CVODES can compute sensitivities with respect to the
 * problem parameters q1 and q2 (internal difference-quotient
 * sensitivity RHS -> the user data uses the pinned FSAUserData
 * wrapper; ARCHITECTURE.md section 3.6).
 * Any of three sensitivity methods (SIMULTANEOUS, STAGGERED, and
 * STAGGERED1) can be used and sensitivities may be included in the
 * error test or not (error control set on SUNTRUE or SUNFALSE,
 * respectively).
 *
 * Execution:
 *    % cvsAdvDiff_FSA_non -nosensi
 *    % cvsAdvDiff_FSA_non -sensi sensi_meth err_con
 * where sensi_meth is one of {sim, stg, stg1} and err_con is one of
 * {t, f}.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]

use cvodes_rs::cvodes::{
    CVode, CVodeCreate, CVodeFree, CVodeGetSens, CVodeInit, CVodeSStolerances,
    CVodeSensEEtolerances, CVodeSensInit1,
};
use cvodes_rs::cvodes_io::{
    CVodeGetLastOrder, CVodeGetLastStep, CVodeGetNumErrTestFails, CVodeGetNumNonlinSolvConvFails,
    CVodeGetNumNonlinSolvIters, CVodeGetNumRhsEvals, CVodeGetNumRhsEvalsSens,
    CVodeGetNumLinSolvSetups, CVodeGetNumSteps, CVodeGetSensNumErrTestFails,
    CVodeGetSensNumLinSolvSetups, CVodeGetSensNumNonlinSolvConvFails,
    CVodeGetSensNumNonlinSolvIters, CVodeGetSensNumRhsEvals, CVodeSetSensDQMethod,
    CVodeSetSensErrCon, CVodeSetSensParams, CVodeSetUserData,
};
use cvodes_rs::cvodes_nls::CVodeSetNonlinearSolver;
use cvodes_rs::cvodes_nls_sim::CVodeSetNonlinearSolverSensSim;
use cvodes_rs::cvodes_nls_stg::CVodeSetNonlinearSolverSensStg;
use cvodes_rs::cvodes_nls_stg1::CVodeSetNonlinearSolverSensStg1;
use cvodes_rs::sundials_utils::fmt_e;
use cvodes_rs::*;

/* Problem Constants */
const XMAX: f64 = 2.0; /* domain boundary           */
const MX: i64 = 10; /* mesh dimension            */
const NEQ: i64 = MX; /* number of equations       */
const ATOL: f64 = 1.0e-5; /* scalar absolute tolerance */
const T0: f64 = 0.0; /* initial time              */
const T1: f64 = 0.5; /* first output time         */
const DTOUT: f64 = 0.5; /* output time increment     */
const NOUT: i32 = 10; /* number of output times    */

const NS: usize = 2;

const ZERO: f64 = 0.0;

/* Type : UserData — C carries {p, dx}; the internal-DQ sensitivity
   convention wraps p in FSAUserData with dx in the .user slot. */

/*
 * f routine. Compute f(t,u).
 */
fn f(_t: f64, u: &NVector, udot: &mut NVector, user_data: &mut UserData) -> i32 {
    /* Extract needed problem constants from data */
    let fsa = user_data.as_ref().unwrap().downcast_ref::<FSAUserData>().unwrap();
    let dx = *fsa.user.downcast_ref::<f64>().unwrap();
    let hordc = fsa.p[0] / (dx * dx);
    let horac = fsa.p[1] / (2.0 * dx);

    /* Loop over all grid points. */
    for i in 0..NEQ as usize {
        /* Extract u at x_i and two neighboring points */
        let ui = u.data[i];
        let ult = if i != 0 { u.data[i - 1] } else { ZERO };
        let urt = if i != (NEQ - 1) as usize { u.data[i + 1] } else { ZERO };

        /* Set diffusion and advection terms and load into udot */
        let hdiff = hordc * (ult - 2.0 * ui + urt);
        let hadv = horac * (urt - ult);
        udot.data[i] = hdiff + hadv;
    }

    0
}

/*
 * Process and verify arguments.
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
 * Set initial conditions in u vector.
 */
fn set_ic(u: &mut NVector, dx: f64) {
    /* Load initial profile into u vector */
    for i in 0..NEQ as usize {
        let x = (i + 1) as f64 * dx;
        u.data[i] = x * (XMAX - x) * (2.0 * x).exp();
    }
}

/*
 * Print current t, step count, order, stepsize, and max norm of solution
 */
fn print_output(cvode_mem: &mut CVodeMem, t: f64, u: &NVector) {
    let mut nst: i64 = 0;
    let mut qu: i32 = 0;
    let mut hu: f64 = 0.0;

    CVodeGetNumSteps(cvode_mem, &mut nst);
    CVodeGetLastOrder(cvode_mem, &mut qu);
    CVodeGetLastStep(cvode_mem, &mut hu);

    println!("{} {:2}  {} {:5}", fmt_e(t, 8, 3), qu, fmt_e(hu, 8, 3), nst);

    print!("                                Solution       ");
    println!("{} ", fmt_e(N_VMaxNorm(u), 12, 4));
}

/*
 * Print max norm of sensitivities
 */
fn print_output_s(uS: &[NVector]) {
    print!("                                Sensitivity 1  ");
    println!("{} ", fmt_e(N_VMaxNorm(&uS[0]), 12, 4));

    print!("                                Sensitivity 2  ");
    println!("{} ", fmt_e(N_VMaxNorm(&uS[1]), 12, 4));
}

/*
 * Print some final statistics located in the CVODES memory
 */
fn print_final_stats(cvode_mem: &mut CVodeMem, sensi: bool, err_con: bool, sensi_meth: i32) {
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

    CVodeGetNumSteps(cvode_mem, &mut nst);
    CVodeGetNumRhsEvals(cvode_mem, &mut nfe);
    CVodeGetNumLinSolvSetups(cvode_mem, &mut nsetups);
    CVodeGetNumErrTestFails(cvode_mem, &mut netf);
    CVodeGetNumNonlinSolvIters(cvode_mem, &mut nni);
    CVodeGetNumNonlinSolvConvFails(cvode_mem, &mut ncfn);

    if sensi {
        CVodeGetSensNumRhsEvals(cvode_mem, &mut nfSe);
        CVodeGetNumRhsEvalsSens(cvode_mem, &mut nfeS);
        CVodeGetSensNumLinSolvSetups(cvode_mem, &mut nsetupsS);
        if err_con {
            CVodeGetSensNumErrTestFails(cvode_mem, &mut netfS);
        } else {
            netfS = 0;
        }
        if sensi_meth == CV_STAGGERED || sensi_meth == CV_STAGGERED1 {
            CVodeGetSensNumNonlinSolvIters(cvode_mem, &mut nniS);
            CVodeGetSensNumNonlinSolvConvFails(cvode_mem, &mut ncfnS);
        } else {
            nniS = 0;
            ncfnS = 0;
        }
    }

    println!("\nFinal Statistics\n");
    println!("nst     = {:5}\n", nst);
    println!("nfe     = {:5}", nfe);
    println!("netf    = {:5}    nsetups  = {:5}", netf, nsetups);
    println!("nni     = {:5}    ncfn     = {:5}", nni, ncfn);

    if sensi {
        println!();
        println!("nfSe    = {:5}    nfeS     = {:5}", nfSe, nfeS);
        println!("netfs   = {:5}    nsetupsS = {:5}", netfS, nsetupsS);
        println!("nniS    = {:5}    ncfnS    = {:5}", nniS, ncfnS);
    }
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
    /* Process arguments */
    let argv: Vec<String> = std::env::args().collect();
    let (sensi, sensi_meth, err_con) = process_args(&argv);

    /* Create SUNDIALS context */
    let sunctx = SUNContext_Create();

    /* Set user data: internal DQ sensitivities use the pinned
       FSAUserData wrapper ({p} + dx in the user slot) */
    let dx = XMAX / ((MX + 1) as f64);
    let p = vec![1.0, 0.5];

    /* Allocate and set initial states */
    let mut u = N_VNew_Serial(NEQ, &sunctx);
    set_ic(&mut u, dx);

    /* Set integration tolerances */
    let reltol = ZERO;
    let abstol = ATOL;

    /* Create CVODES object */
    let mut cvode_mem = CVodeCreate(CV_ADAMS, &sunctx);

    let retval = CVodeSetUserData(
        &mut cvode_mem,
        Some(Box::new(FSAUserData {
            p: p.clone(),
            user: Box::new(dx),
        })),
    );
    if check_retval(retval, "CVodeSetUserData") {
        return;
    }

    /* Allocate CVODES memory */
    let retval = CVodeInit(&mut cvode_mem, f, T0, &u);
    if check_retval(retval, "CVodeInit") {
        return;
    }

    let retval = CVodeSStolerances(&mut cvode_mem, reltol, abstol);
    if check_retval(retval, "CVodeSStolerances") {
        return;
    }

    /* create fixed point nonlinear solver object and attach to CVode */
    let nls = SUNNonlinSol_FixedPoint(&u, 0, &sunctx);
    let retval = CVodeSetNonlinearSolver(&mut cvode_mem, nls);
    if check_retval(retval, "CVodeSetNonlinearSolver") {
        return;
    }

    println!("\n1-D advection-diffusion equation, mesh size ={:3}", MX);

    /* Sensitivity-related settings */
    let mut uS: Vec<NVector> = Vec::new();
    if sensi {
        let plist: Vec<i32> = (0..NS as i32).collect();
        let pbar: Vec<f64> = plist.iter().map(|&i| p[i as usize]).collect();

        uS = (0..NS).map(|_| N_VClone(&u)).collect();
        for us in uS.iter_mut() {
            N_VConst(ZERO, us);
        }

        let retval = CVodeSensInit1(&mut cvode_mem, NS as i32, sensi_meth, None, &uS);
        if check_retval(retval, "CVodeSensInit1") {
            return;
        }

        let retval = CVodeSensEEtolerances(&mut cvode_mem);
        if check_retval(retval, "CVodeSensEEtolerances") {
            return;
        }

        let retval = CVodeSetSensErrCon(&mut cvode_mem, err_con);
        if check_retval(retval, "CVodeSetSensErrCon") {
            return;
        }

        let retval = CVodeSetSensDQMethod(&mut cvode_mem, CV_CENTERED, ZERO);
        if check_retval(retval, "CVodeSetSensDQMethod") {
            return;
        }

        let retval = CVodeSetSensParams(&mut cvode_mem, Some(&p), Some(&pbar), Some(&plist));
        if check_retval(retval, "CVodeSetSensParams") {
            return;
        }

        /* create sensitivity fixed point nonlinear solver object and
           attach to CVode */
        if sensi_meth == CV_SIMULTANEOUS {
            let nls_sens = SUNNonlinSol_FixedPointSens((NS + 1) as i32, &u, 0, &sunctx);
            let retval = CVodeSetNonlinearSolverSensSim(&mut cvode_mem, nls_sens);
            if check_retval(retval, "CVodeSetNonlinearSolver") {
                return;
            }
        } else if sensi_meth == CV_STAGGERED {
            let nls_sens = SUNNonlinSol_FixedPointSens(NS as i32, &u, 0, &sunctx);
            let retval = CVodeSetNonlinearSolverSensStg(&mut cvode_mem, nls_sens);
            if check_retval(retval, "CVodeSetNonlinearSolver") {
                return;
            }
        } else {
            let nls_sens = SUNNonlinSol_FixedPoint(&u, 0, &sunctx);
            let retval = CVodeSetNonlinearSolverSensStg1(&mut cvode_mem, nls_sens);
            if check_retval(retval, "CVodeSetNonlinearSolver") {
                return;
            }
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

    /* In loop over output points, call CVode, print results, test for
       error */
    println!("\n");
    println!("============================================================");
    println!("     T     Q       H      NST                    Max norm   ");
    println!("============================================================");

    let mut t = T0;
    let mut tout = T1;
    for _iout in 1..=NOUT {
        let retval = CVode(&mut cvode_mem, tout, &mut u, &mut t, CV_NORMAL);
        if check_retval(retval, "CVode") {
            break;
        }
        print_output(&mut cvode_mem, t, &u);
        if sensi {
            let retval = CVodeGetSens(&cvode_mem, &mut t, &mut uS);
            if check_retval(retval, "CVodeGetSens") {
                break;
            }
            print_output_s(&uS);
        }
        println!("------------------------------------------------------------");

        tout += DTOUT;
    }

    /* Print final statistics */
    print_final_stats(&mut cvode_mem, sensi, err_con, sensi_meth);

    /* Free memory */
    drop(u);
    drop(uS);
    CVodeFree(cvode_mem);
}
