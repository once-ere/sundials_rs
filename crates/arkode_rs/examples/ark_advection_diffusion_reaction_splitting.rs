/*---------------------------------------------------------------
 * Translation of
 * examples/arkode/C_serial/ark_advection_diffusion_reaction_splitting.c
 * (SUNDIALS 7.7.0).
 *
 * Example problem:
 *
 * The following test simulates a simple 1D advection-diffusion-
 * reaction equation,
 *    u_t = (a/2)*(u^2)_x + b*u_xx + c*(u - u^3)
 * for t in [0, 1], x in [0, 1], with initial conditions
 *    u(0,x) = u_0
 * and Dirichlet boundary conditions at x=0 and x=1
 *    u(0,t) = u(1,t) = u_0
 *
 * The spatial derivatives are computed using second-order
 * centered differences, with the data distributed over N
 * points (excluding boundary points) on a uniform spatial grid.
 *
 * This program solves the problem with an operator splitting
 * method where advection is treated with a strong stability
 * preserving ERK method, diffusion is treated with a DIRK
 * method, and reaction is treated with a different ERK method.
 *
 * Outputs are printed at equal intervals, and run statistics are
 * printed at the end.
 *
 * Note on ownership: the C example keeps the partition integrator
 * pointers alongside their SUNSteppers; in this port each stepper
 * owns its wrapped integrator (and the splitting integrator owns
 * the steppers), so the final per-partition statistics are printed
 * by borrowing each integrator back out of the outer step memory.
 * The shared C udata pointer becomes one identical UserData clone
 * per integrator (the problem data is immutable during the run).
 *---------------------------------------------------------------*/
#![allow(non_snake_case)]

use arkode_rs::arkode::{ARKodeEvolve, ARKodeFree};
use arkode_rs::arkode_arkstep::ARKStepCreate;
use arkode_rs::arkode_erkstep::ERKStepCreate;
use arkode_rs::arkode_erkstep_io::ERKStepSetTableNum;
use arkode_rs::arkode_butcher_erk::ARKODE_SHU_OSHER_3_2_3;
use arkode_rs::arkode_io::{
    ARKodePrintAllStats, ARKodeSetFixedStep, ARKodeSetLinear, ARKodeSetOrder,
    ARKodeSetStopTime, ARKodeSetUserData,
};
use arkode_rs::arkode_ls::{ARKodeSetJacFn, ARKodeSetLinearSolver};
use arkode_rs::arkode_splittingstep::SplittingStepCreate;
use arkode_rs::arkode_splittingstep_impl::ARKodeSplittingStepMem;
use arkode_rs::arkode_sunstepper::ARKodeCreateSUNStepper;
use arkode_rs::sundials_utils::{fmt_f, fmt_g};
use arkode_rs::*;

/* user data structure */
#[derive(Clone)]
struct AdrData {
    N: i64,  /* number of grid points (excluding boundaries) */
    dx: f64, /* mesh spacing */
    a: f64,  /* advection coefficient */
    b: f64,  /* diffusion coefficient */
    c: f64,  /* reaction coefficient */
    u0: f64, /* initial and boundary values */
}

/* f routine to compute the advection RHS function. */
fn f_advection(_t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let udata = user_data.as_ref().unwrap().downcast_ref::<AdrData>().unwrap();

    let coeff = udata.a / (4.0 * udata.dx);
    let u0_sqr = udata.u0 * udata.u0;
    let n = udata.N as usize;

    /* Left boundary */
    ydot.data[0] = coeff * (y.data[1] * y.data[1] - u0_sqr);
    /* Interior */
    for i in 1..n - 1 {
        ydot.data[i] = coeff * (y.data[i + 1] * y.data[i + 1] - y.data[i - 1] * y.data[i - 1]);
    }
    /* Right boundary */
    ydot.data[n - 1] = coeff * (u0_sqr - y.data[n - 1] * y.data[n - 1]);

    0
}

/* f routine to compute the diffusion RHS function. */
fn f_diffusion(_t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let udata = user_data.as_ref().unwrap().downcast_ref::<AdrData>().unwrap();

    let coeff = udata.b / (udata.dx * udata.dx);
    let n = udata.N as usize;

    /* Left boundary */
    ydot.data[0] = coeff * (udata.u0 - 2.0 * y.data[0] + y.data[1]);
    /* Interior */
    for i in 1..n - 1 {
        ydot.data[i] = coeff * (y.data[i + 1] - 2.0 * y.data[i] + y.data[i - 1]);
    }
    /* Right boundary */
    ydot.data[n - 1] = coeff * (y.data[n - 2] - 2.0 * y.data[n - 1] + udata.u0);

    0
}

/* Routine to compute the diffusion Jacobian function. */
#[allow(clippy::too_many_arguments)]
fn jac_diffusion(
    _t: f64,
    _y: &NVector,
    _fy: &NVector,
    jac: &mut SUNMatrix,
    user_data: &mut UserData,
    _tmp1: &mut NVector,
    _tmp2: &mut NVector,
    _tmp3: &mut NVector,
) -> i32 {
    let udata = user_data.as_ref().unwrap().downcast_ref::<AdrData>().unwrap();
    let coeff = udata.b / (udata.dx * udata.dx);

    let jm = match jac {
        SUNMatrix::Band(m) => m,
        _ => return 1,
    };

    jm.set(0, 0, -2.0 * coeff);
    for i in 1..udata.N {
        jm.set(i - 1, i, coeff);
        jm.set(i, i, -2.0 * coeff);
        jm.set(i, i - 1, coeff);
    }

    0
}

/* f routine to compute the reaction RHS function. */
fn f_reaction(_t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let udata = user_data.as_ref().unwrap().downcast_ref::<AdrData>().unwrap();

    for i in 0..udata.N as usize {
        ydot.data[i] = udata.c * y.data[i] * (1.0 - y.data[i] * y.data[i]);
    }

    0
}

/* Check if a SUNDIALS function returned a negative flag */
fn check_flag(flag: i32, funcname: &str) -> bool {
    if flag < 0 {
        eprintln!("\nSUNDIALS_ERROR: {}() failed with flag = {}\n", funcname, flag);
        return true;
    }
    false
}

fn main() {
    /* Problem parameters */
    let T0: f64 = 0.0;
    let Tf: f64 = 1.0;
    let DT: f64 = 0.06;
    let N: i64 = 128;
    let udata = AdrData {
        N,
        dx: 1.0 / (N + 1) as f64,
        a: 1.0,
        b: 0.125,
        c: 4.0,
        u0: 0.1,
    };

    println!("\n1D Advection-Diffusion-Reaction PDE test problem:");
    println!("  N = {}", udata.N);
    println!("  advection coefficient = {}", fmt_g(udata.a, 0, 6));
    println!("  diffusion coefficient = {}", fmt_g(udata.b, 0, 6));
    println!("  reaction coefficient  = {}\n", fmt_g(udata.c, 0, 6));

    /* Create the SUNDIALS context object for this simulation */
    let ctx = SUNContext_Create();

    /* Initialize vector with initial condition */
    let mut y = N_VNew_Serial(udata.N, &ctx);
    N_VConst(udata.u0, &mut y);

    /* Create advection integrator */
    let mut advection_mem = ERKStepCreate(f_advection, T0, &y, &ctx);

    let flag = ARKodeSetUserData(&mut advection_mem, Some(Box::new(udata.clone())));
    if check_flag(flag, "ARKodeSetUserData") {
        return;
    }

    /* Choose a strong stability preserving method for advecton */
    let flag = ERKStepSetTableNum(&mut advection_mem, ARKODE_SHU_OSHER_3_2_3);
    if check_flag(flag, "ERKStepSetTableNum") {
        return;
    }

    let advection_stepper = ARKodeCreateSUNStepper(advection_mem);

    /* Create diffusion integrator */
    let mut diffusion_mem =
        ARKStepCreate(None, Some(f_diffusion), T0, &y, &ctx).expect("ARKStepCreate");

    let flag = ARKodeSetUserData(&mut diffusion_mem, Some(Box::new(udata.clone())));
    if check_flag(flag, "ARKodeSetUserData") {
        return;
    }

    let flag = ARKodeSetOrder(&mut diffusion_mem, 3);
    if check_flag(flag, "ARKStepSetOrder") {
        return;
    }

    let jac_mat = SUNBandMatrix(udata.N, 1, 1, &ctx);
    let ls = SUNLinSol_Band(&y, &jac_mat, &ctx);

    let flag = ARKodeSetLinearSolver(&mut diffusion_mem, ls, Some(jac_mat));
    if check_flag(flag, "ARKStepSetOrder") {
        return;
    }

    let flag = ARKodeSetJacFn(&mut diffusion_mem, Some(jac_diffusion));
    if check_flag(flag, "ARKodeSetJacFn") {
        return;
    }

    let flag = ARKodeSetLinear(&mut diffusion_mem, 0);
    if check_flag(flag, "ARKodeSetLinear") {
        return;
    }

    let diffusion_stepper = ARKodeCreateSUNStepper(diffusion_mem);

    /* Create reaction integrator */
    let mut reaction_mem = ERKStepCreate(f_reaction, T0, &y, &ctx);

    let flag = ARKodeSetUserData(&mut reaction_mem, Some(Box::new(udata.clone())));
    if check_flag(flag, "ARKodeSetUserData") {
        return;
    }

    let flag = ARKodeSetOrder(&mut reaction_mem, 3);
    if check_flag(flag, "ARKodeSetOrder") {
        return;
    }

    let reaction_stepper = ARKodeCreateSUNStepper(reaction_mem);

    /* Create operator splitting integrator */
    let steppers = vec![advection_stepper, diffusion_stepper, reaction_stepper];
    let mut arkode_mem =
        SplittingStepCreate(steppers, 3, T0, &y, &ctx).expect("SplittingStepCreate");

    let flag = ARKodeSetFixedStep(&mut arkode_mem, DT);
    if check_flag(flag, "ARKodeSetFixedStep") {
        return;
    }

    let flag = ARKodeSetStopTime(&mut arkode_mem, Tf);
    if check_flag(flag, "ARKodeSetStopTime") {
        return;
    }

    /* Evolve solution in time */
    let mut tret = T0;
    println!("        t      ||u||_rms");
    println!("   ----------------------");
    println!(
        "  {}  {}",
        fmt_f(tret, 10, 6),
        fmt_f((N_VDotProd(&y, &y) / udata.N as f64).sqrt(), 10, 6)
    );
    while tret < Tf {
        let flag = ARKodeEvolve(&mut arkode_mem, Tf, &mut y, &mut tret, ARK_ONE_STEP);
        if check_flag(flag, "ARKodeEvolve") {
            return;
        }
        println!(
            "  {}  {}",
            fmt_f(tret, 10, 6),
            fmt_f((N_VDotProd(&y, &y) / udata.N as f64).sqrt(), 10, 6)
        );
    }
    println!("   ----------------------");

    /* Print statistics */
    let mut stdout = std::io::stdout();
    println!("\nSplitting Stepper Statistics:");
    let flag = ARKodePrintAllStats(&mut arkode_mem, &mut stdout, SUN_OUTPUTFORMAT_TABLE);
    if check_flag(flag, "ARKodePrintAllStats") {
        return;
    }

    /* (borrow each partition integrator back out of the outer step
       memory for its statistics) */
    {
        let step_mem = arkode_mem
            .step_mem
            .as_mut()
            .unwrap()
            .downcast_mut::<ARKodeSplittingStepMem>()
            .unwrap();

        println!("\nAdvection Stepper Statistics:");
        let advection_mem = step_mem.steppers[0]
            .content
            .as_mut()
            .unwrap()
            .downcast_mut::<ARKodeMem>()
            .unwrap();
        let flag = ARKodePrintAllStats(advection_mem, &mut stdout, SUN_OUTPUTFORMAT_TABLE);
        if check_flag(flag, "ARKodePrintAllStats") {
            return;
        }

        println!("\nDiffusion Stepper Statistics:");
        let diffusion_mem = step_mem.steppers[1]
            .content
            .as_mut()
            .unwrap()
            .downcast_mut::<ARKodeMem>()
            .unwrap();
        let flag = ARKodePrintAllStats(diffusion_mem, &mut stdout, SUN_OUTPUTFORMAT_TABLE);
        if check_flag(flag, "ARKodePrintAllStats") {
            return;
        }

        println!("\nReaction Stepper Statistics:");
        let reaction_mem = step_mem.steppers[2]
            .content
            .as_mut()
            .unwrap()
            .downcast_mut::<ARKodeMem>()
            .unwrap();
        let flag = ARKodePrintAllStats(reaction_mem, &mut stdout, SUN_OUTPUTFORMAT_TABLE);
        if check_flag(flag, "ARKodePrintAllStats") {
            return;
        }
    }

    /* Clean up and return with successful completion */
    drop(y);
    let mut slot = Some(arkode_mem);
    ARKodeFree(&mut slot);
}
