/* -----------------------------------------------------------------------------
 * Translation of examples/arkode/C_serial/ark_conserved_exp_entropy_ark.c
 * (SUNDIALS 7.7.0).
 *
 * This example problem is adapted from:
 *
 * H. Ranocha, M. Sayyari, L. Dalcin, M. Parsani, and D.I. Ketcheson,
 * "Relaxation Runge-Kutta Methods: Fully-Discrete Explicit Entropy-Stable
 * Schemes for the Compressible Euler and Navier-Stokes Equations," SIAM Journal
 * on Scientific Computing, 42(2), 2020, https://doi.org/10.1137/19M1263480.
 * -----------------------------------------------------------------------------
 * This example evolves system
 *
 *   du/dt = -exp(v)
 *   dv/dt =  exp(u)
 *
 * for t in the interval [0, 5] with the initial condition
 *
 *   u(0) = 1.0
 *   v(0) = 0.5
 *
 * The system has the analytic solution
 *
 *   u = log(e + e^(3/2)) - log(b)
 *   v = log(a * e^(a * t)) - log(b)
 *
 * where log is the natural logarithm, a = sqrt(e) + e, and
 * b = sqrt(e) + e^(a * t).
 *
 * The conserved exponential entropy for the system is given by
 * ent(u,v) = exp(u) + exp(v) with the Jacobian
 * ent'(u,v) = [ de/du de/dv ]^T = [ exp(u) exp(v) ]^T.
 *
 * The problem is advanced in time with an explicit or implicit relaxed
 * Runge-Kutta method from ARKStep to ensure conservation of the entropy.
 * ---------------------------------------------------------------------------*/
#![allow(non_snake_case)]

use std::io::Write;

use arkode_rs::arkode::{ARKodeEvolve, ARKodeFree, ARKodeSStolerances};
use arkode_rs::arkode_arkstep::ARKStepCreate;
use arkode_rs::arkode_io::{
    ARKodeGetNumErrTestFails, ARKodeGetNumLinSolvSetups, ARKodeGetNumNonlinSolvConvFails,
    ARKodeGetNumNonlinSolvIters, ARKodeGetNumRhsEvals, ARKodeGetNumStepAttempts,
    ARKodeGetNumSteps, ARKodeSetFixedStep, ARKodeSetNonlinConvCoef, ARKodeSetOrder,
};
use arkode_rs::arkode_ls::{
    ARKodeGetNumJacEvals, ARKodeGetNumLinRhsEvals, ARKodeSetJacFn, ARKodeSetLinearSolver,
};
use arkode_rs::arkode_relaxation::{
    ARKodeGetNumRelaxBoundFails, ARKodeGetNumRelaxFails, ARKodeGetNumRelaxFnEvals,
    ARKodeGetNumRelaxJacEvals, ARKodeGetNumRelaxSolveFails, ARKodeGetNumRelaxSolveIters,
    ARKodeSetRelaxFn,
};
use arkode_rs::sundials_utils::fmt_e;
use arkode_rs::*;

/* Value of the natural number e */
#[allow(clippy::approx_constant)]
const EVAL: f64 = 2.718281828459045235360287471352662497757247093699959574966;

/* ----------------------- *
 * User-supplied functions *
 * ----------------------- */

/* ODE RHS function f(t,y). */
fn f(_t: f64, y: &NVector, ydot: &mut NVector, _user_data: &mut UserData) -> i32 {
    ydot.data[0] = -(y.data[1]).exp();
    ydot.data[1] = (y.data[0]).exp();

    0
}

/* ODE RHS Jacobian function J(t,y) = df/dy. */
#[allow(clippy::too_many_arguments)]
fn jac(
    _t: f64,
    y: &NVector,
    _fy: &NVector,
    j: &mut SUNMatrix,
    _user_data: &mut UserData,
    _tmp1: &mut NVector,
    _tmp2: &mut NVector,
    _tmp3: &mut NVector,
) -> i32 {
    if let SUNMatrix::Dense(dm) = j {
        /* column 0 */
        dm.data[0] = 0.0;
        dm.data[1] = (y.data[0]).exp();

        /* column 1 */
        dm.data[2] = -(y.data[1]).exp();
        dm.data[3] = 0.0;
    }

    0
}

/* Entropy function e(y) */
fn ent(y: &NVector, e: &mut f64, _user_data: &mut UserData) -> i32 {
    *e = (y.data[0]).exp() + (y.data[1]).exp();

    0
}

/* Entropy function Jacobian Je(y) = de/dy */
fn jac_ent(y: &NVector, j: &mut NVector, _user_data: &mut UserData) -> i32 {
    j.data[0] = (y.data[0]).exp();
    j.data[1] = (y.data[1]).exp();

    0
}

/* ----------------- *
 * Utility functions *
 * ----------------- */

/* Analytic solution */
fn ans(t: f64, y: &mut NVector) -> i32 {
    let a = EVAL.sqrt() + EVAL;
    let b = EVAL.sqrt() + (a * t).exp();

    y.data[0] = (EVAL + (1.5f64).exp()).ln() - b.ln();
    y.data[1] = (a * (a * t).exp()).ln() - b.ln();

    0
}

/* Check for an unrecoverable (negative) return flag from a SUNDIALS function */
fn check_flag(flag: i32, funcname: &str) -> bool {
    if flag < 0 {
        eprintln!("ERROR: {}() returned {}", funcname, flag);
        return true;
    }
    false
}

/* ------------ *
 * Main Program *
 * ------------ */

fn main() {
    /* Initial and final times */
    let t0: f64 = 0.0;
    let tf: f64 = 5.0;

    /* Relative and absolute tolerances */
    let reltol: f64 = 1.0e-6;
    let abstol: f64 = 1.0e-10;

    /* Command line options */
    let argv: Vec<String> = std::env::args().collect();
    let relax: i32 = if argv.len() > 1 { argv[1].parse().unwrap_or(0) } else { 1 };
    let implicit: i32 = if argv.len() > 2 { argv[2].parse().unwrap_or(0) } else { 1 };
    let fixed_h: f64 = if argv.len() > 3 { argv[3].parse().unwrap_or(0.0) } else { 0.0 };

    /* -------------------- *
     * Output Problem Setup *
     * -------------------- */

    println!("\nConserved Exponential Entropy problem:");
    if implicit != 0 {
        println!("   method     = DIRK");
    } else {
        println!("   method     = ERK");
    }
    println!("   reltol     = {}", fmt_e(reltol, 0, 1));
    println!("   abstol     = {}", fmt_e(abstol, 0, 1));
    if fixed_h > 0.0 {
        println!("   fixed h    = {}", fmt_e(fixed_h, 0, 1));
    }
    if relax != 0 {
        println!("   relaxation = ON");
    } else {
        println!("   relaxation = OFF");
    }
    println!();

    /* ------------ *
     * Setup ARKODE *
     * ------------ */

    /* Create the SUNDIALS context object for this simulation */
    let ctx = SUNContext_Create();

    /* Create serial vector and set the initial condition values */
    let mut y = N_VNew_Serial(2, &ctx);
    y.data[0] = 1.0;
    y.data[1] = 0.5;

    let mut ytrue = N_VClone(&y);

    /* Initialize ARKStep */
    let mut arkode_mem = if implicit != 0 {
        ARKStepCreate(None, Some(f), t0, &y, &ctx).expect("ARKStepCreate")
    } else {
        ARKStepCreate(Some(f), None, t0, &y, &ctx).expect("ARKStepCreate")
    };

    /* (no user data is attached; direct entropy calls below pass a
       standalone None, matching C's NULL) */
    let mut no_user_data: UserData = None;

    /* Set order */
    let retval = ARKodeSetOrder(&mut arkode_mem, 2);
    if check_flag(retval, "ARKodeSetOrder") {
        return;
    }

    /* Specify tolerances */
    let retval = ARKodeSStolerances(&mut arkode_mem, reltol, abstol);
    if check_flag(retval, "ARKodeSStolerances") {
        return;
    }

    if relax != 0 {
        /* Enable relaxation methods */
        let retval = ARKodeSetRelaxFn(&mut arkode_mem, Some(ent), Some(jac_ent));
        if check_flag(retval, "ARKodeSetRelaxFn") {
            return;
        }
    }

    if implicit != 0 {
        /* Create dense matrix and linear solver */
        let a_mat = SUNDenseMatrix(2, 2, &ctx);
        let ls = SUNLinSol_Dense(&y, &a_mat, &ctx);

        /* Attach the matrix and linear solver */
        let retval = ARKodeSetLinearSolver(&mut arkode_mem, ls, Some(a_mat));
        if check_flag(retval, "ARKodeSetLinearSolver") {
            return;
        }

        /* Set Jacobian routine */
        let retval = ARKodeSetJacFn(&mut arkode_mem, Some(jac));
        if check_flag(retval, "ARKodeSetJacFn") {
            return;
        }

        /* Tighten nonlinear solver tolerance */
        let retval = ARKodeSetNonlinConvCoef(&mut arkode_mem, 0.01);
        if check_flag(retval, "ARKodeSetNonlinConvCoef") {
            return;
        }
    }

    if fixed_h > 0.0 {
        let retval = ARKodeSetFixedStep(&mut arkode_mem, fixed_h);
        if check_flag(retval, "ARKodeSetFixedStep") {
            return;
        }
    }

    /* Open output stream for results, output comment line */
    let mut ufid = std::fs::File::create("ark_conserved_exp_entropy_ark.txt").expect("fopen");
    let _ = writeln!(ufid, "# vars: t u v entropy u_err v_err entropy_error");

    /* --------------- *
     * Advance in Time *
     * --------------- */

    /* Initial time */
    let mut t = t0;

    /* Output the initial condition and entropy */
    let mut ent0 = 0.0;
    let retval = ent(&y, &mut ent0, &mut no_user_data);
    if check_flag(retval, "Ent") {
        return;
    }

    let _ = writeln!(
        ufid,
        "{} {} {} {} {} {} {}",
        fmt_e(t0, 23, 16),
        fmt_e(y.data[0], 23, 16),
        fmt_e(y.data[1], 23, 16),
        fmt_e(ent0, 23, 16),
        fmt_e(0.0, 23, 16),
        fmt_e(0.0, 23, 16),
        fmt_e(0.0, 23, 16)
    );

    println!(
        " step   t              u              v              e              delta e"
    );
    println!(
        " -------------------------------------------------------------------------------"
    );
    println!(
        "{:5} {} {} {} {} {}",
        0,
        fmt_e(t, 14, 6),
        fmt_e(y.data[0], 14, 6),
        fmt_e(y.data[1], 14, 6),
        fmt_e(ent0, 14, 6),
        fmt_e(0.0, 14, 6)
    );

    let mut flag = 0;
    let mut nst: i64 = 0;
    while t < tf {
        /* Evolve in time */
        flag = ARKodeEvolve(&mut arkode_mem, tf, &mut y, &mut t, ARK_ONE_STEP);
        if check_flag(flag, "ARKodeEvolve") {
            break;
        }

        /* Output solution and errors */
        let mut ent_t = 0.0;
        let retval = ent(&y, &mut ent_t, &mut no_user_data);
        if check_flag(retval, "Ent") {
            return;
        }

        let retval = ans(t, &mut ytrue);
        if check_flag(retval, "ans") {
            return;
        }

        let ent_err = ent_t - ent0;
        let u_err = y.data[0] - ytrue.data[0];
        let v_err = y.data[1] - ytrue.data[1];

        /* Output to the screen periodically */
        let retval = ARKodeGetNumSteps(&mut arkode_mem, &mut nst);
        check_flag(retval, "ARKodeGetNumSteps");

        if nst % 40 == 0 {
            println!(
                "{:5} {} {} {} {} {}",
                nst,
                fmt_e(t, 14, 6),
                fmt_e(y.data[0], 14, 6),
                fmt_e(y.data[1], 14, 6),
                fmt_e(ent_t, 14, 6),
                fmt_e(ent_err, 14, 6)
            );
        }

        /* Write all steps to file */
        let _ = writeln!(
            ufid,
            "{} {} {} {} {} {} {}",
            fmt_e(t, 23, 16),
            fmt_e(y.data[0], 23, 16),
            fmt_e(y.data[1], 23, 16),
            fmt_e(ent_t, 23, 16),
            fmt_e(u_err, 23, 16),
            fmt_e(v_err, 23, 16),
            fmt_e(ent_err, 23, 16)
        );
    }

    println!(
        " -------------------------------------------------------------------------------"
    );
    drop(ufid);

    /* ------------ *
     * Output Stats *
     * ------------ */

    /* Get final statistics on how the solve progressed */
    let mut nst_a: i64 = 0;
    let mut netf: i64 = 0;
    let mut nfe: i64 = 0;
    let mut nfi: i64 = 0;

    let retval = ARKodeGetNumSteps(&mut arkode_mem, &mut nst);
    check_flag(retval, "ARKodeGetNumSteps");

    let retval = ARKodeGetNumStepAttempts(&mut arkode_mem, &mut nst_a);
    check_flag(retval, "ARKodeGetNumStepAttempts");

    let retval = ARKodeGetNumErrTestFails(&mut arkode_mem, &mut netf);
    check_flag(retval, "ARKodeGetNumErrTestFails");

    let retval = ARKodeGetNumRhsEvals(&mut arkode_mem, 0, &mut nfe);
    check_flag(retval, "ARKodeGetNumRhsEvals");

    let retval = ARKodeGetNumRhsEvals(&mut arkode_mem, 1, &mut nfi);
    check_flag(retval, "ARKodeGetNumRhsEvals");

    println!("\nFinal Solver Statistics:");
    println!("   Internal solver steps = {} (attempted = {})", nst, nst_a);
    println!("   Total number of error test failures = {}", netf);
    println!("   Total RHS evals:  Fe = {},  Fi = {}", nfe, nfi);

    if implicit != 0 {
        let mut nni: i64 = 0;
        let mut ncfn: i64 = 0;
        let mut nsetups: i64 = 0;
        let mut nje: i64 = 0;
        let mut nfeLS: i64 = 0;

        let retval = ARKodeGetNumNonlinSolvIters(&mut arkode_mem, &mut nni);
        check_flag(retval, "ARKodeGetNumNonlinSolvIters");

        let retval = ARKodeGetNumNonlinSolvConvFails(&mut arkode_mem, &mut ncfn);
        check_flag(retval, "ARKodeGetNumNonlinSolvConvFails");

        let retval = ARKodeGetNumLinSolvSetups(&mut arkode_mem, &mut nsetups);
        check_flag(retval, "ARKodeGetNumLinSolvSetups");

        let retval = ARKodeGetNumJacEvals(&mut arkode_mem, &mut nje);
        check_flag(retval, "ARKodeGetNumJacEvals");

        let retval = ARKodeGetNumLinRhsEvals(&mut arkode_mem, &mut nfeLS);
        check_flag(retval, "ARKodeGetNumLinRhsEvals");

        println!("   Total number of Newton iterations = {}", nni);
        println!("   Total number of linear solver convergence failures = {}", ncfn);
        println!("   Total linear solver setups = {}", nsetups);
        println!("   Total number of Jacobian evaluations = {}", nje);
        println!("   Total RHS evals for setting up the linear system = {}", nfeLS);
    }

    if relax != 0 {
        let mut nre: i64 = 0;
        let mut nrje: i64 = 0;
        let mut nrf: i64 = 0;
        let mut nrbf: i64 = 0;
        let mut nrnlsf: i64 = 0;
        let mut nrnlsi: i64 = 0;

        let retval = ARKodeGetNumRelaxFnEvals(&mut arkode_mem, &mut nre);
        check_flag(retval, "ARKodeGetNumRelaxFnEvals");

        let retval = ARKodeGetNumRelaxJacEvals(&mut arkode_mem, &mut nrje);
        check_flag(retval, "ARKodeGetNumRelaxJacEvals");

        let retval = ARKodeGetNumRelaxFails(&mut arkode_mem, &mut nrf);
        check_flag(retval, "ARKodeGetNumRelaxFails");

        let retval = ARKodeGetNumRelaxBoundFails(&mut arkode_mem, &mut nrbf);
        check_flag(retval, "ARKodeGetNumRelaxBoundFails");

        let retval = ARKodeGetNumRelaxSolveFails(&mut arkode_mem, &mut nrnlsf);
        check_flag(retval, "ARKodeGetNumRelaxSolveFails");

        let retval = ARKodeGetNumRelaxSolveIters(&mut arkode_mem, &mut nrnlsi);
        check_flag(retval, "ARKodeGetNumRelaxSolveIters");

        println!("   Total Relaxation Fn evals    = {}", nre);
        println!("   Total Relaxation Jac evals   = {}", nrje);
        println!("   Total Relaxation fails       = {}", nrf);
        println!("   Total Relaxation bound fails = {}", nrbf);
        println!("   Total Relaxation NLS fails   = {}", nrnlsf);
        println!("   Total Relaxation NLS iters   = {}", nrnlsi);
    }
    println!();

    /* -------- *
     * Clean up *
     * -------- */

    /* Free ARKode integrator and SUNDIALS objects */
    let mut slot = Some(arkode_mem);
    ARKodeFree(&mut slot);
    drop(y);
    drop(ytrue);

    std::process::exit(if flag < 0 { flag } else { 0 });
}
