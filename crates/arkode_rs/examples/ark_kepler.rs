/* -----------------------------------------------------------------
 * Translation of examples/arkode/C_serial/ark_kepler.c (+ the
 * ark_kepler.h helpers) from SUNDIALS 7.7.0.
 *
 * We consider the Kepler problem. We choose one body to be the
 * center of our coordinate system and then we use the coordinates
 * q = (q1, q2) to represent the position of the second body
 * relative to the first (center). This yields the ODE:
 *    dq/dt = [ p1 ]
 *            [ p2 ]
 *    dp/dt = [ -q1 / (q1^2 + q2^2)^(3/2) ]
 *          = [ -q2 / (q1^2 + q2^2)^(3/2) ]
 * with the initial conditions
 *    q(0) = [ 1 - e ],  p(0) = [        0          ]
 *           [   0   ]          [ sqrt((1+e)/(1-e)) ]
 * where e = 0.6 is the eccentricity.
 *
 * The Hamiltonian for the system,
 *    H(p,q) = 1/2 * (p1^2 + p2^2) - 1/sqrt(q1^2 + q2^2)
 * is conserved as well as the angular momentum,
 *    L(p,q) = q1*p2 - q2*p1.
 *
 * By default we solve the problem by letting y = [ q, p ]^T then
 * using a 4th order symplectic integrator via the SPRKStep
 * time-stepper of ARKODE with a fixed time-step size.
 *
 * The rootfinding feature of SPRKStep is used to count the number
 * of complete orbits (g(q) = q2).
 *
 * CLI arguments (same as C):
 *   --step-mode <fixed, adapt>
 *   --stepper <SPRK, ERK>
 *   --method <string>
 *   --use-compensated-sums
 *   --disable-tstop
 *   --dt <Real>
 *   --tf <Real>
 *   --nout <int>
 *   --count-orbits
 *   --check-order
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]

use arkode_rs::arkode::{ARKodeEvolve, ARKodeFree, ARKodeSStolerances};
use arkode_rs::arkode_arkstep::ARKStepCreate;
use arkode_rs::arkode_arkstep_io::ARKStepSetTableName;
use arkode_rs::arkode_io::{
    ARKodeGetRootInfo, ARKodePrintAllStats, ARKodeSetFixedStep, ARKodeSetMaxNumSteps,
    ARKodeSetStopTime, ARKodeSetUseCompensatedSums, ARKodeSetUserData,
};
use arkode_rs::arkode_root::ARKodeRootInit;
use arkode_rs::arkode_sprk::{ARKodeSPRKTable_Free, ARKodeSPRKTable_LoadByName};
use arkode_rs::arkode_sprkstep::SPRKStepCreate;
use arkode_rs::arkode_sprkstep_io::SPRKStepSetMethodName;
use arkode_rs::sundials_utils::{fmt_e, fmt_f, fmt_g};
use arkode_rs::*;
use std::io::Write;

const NUM_DT: usize = 8;

struct KeplerUserData {
    #[allow(dead_code)] /* C stores ecc in user data; the RHS fns don't read it */
    ecc: f64,
}

struct ProblemResult {
    sol: NVector,
    energy_error: f64,
}

#[derive(Clone)]
struct ProgramArgs {
    step_mode: i32,
    stepper: i32,
    num_output_times: i32,
    use_compsums: i32,
    use_tstop: i32,
    count_orbits: i32,
    check_order: i32,
    dt: f64,
    tf: f64,
    method_name: String,
}

fn PrintHelp() {
    eprintln!(
        "ark_kepler: an ARKODE example demonstrating the SPRKStep time-stepping module solving the Kepler problem"
    );
    eprintln!("  --step-mode <fixed, adapt>  should we use a fixed time-step or adaptive time-step (default fixed)");
    eprintln!("  --stepper <SPRK, ERK>       should we use SPRKStep or ARKStep with an ERK method (default SPRK)");
    eprintln!("  --method <string>           which method to use (default ARKODE_SPRK_MCLACHLAN_4_4)");
    eprintln!("  --use-compensated-sums      turns on compensated summation in ARKODE where applicable");
    eprintln!("  --disable-tstop             turns off tstop mode");
    eprintln!("  --dt <Real>                 the fixed-time step size to use if fixed time stepping is turned on (default 0.01)");
    eprintln!("  --tf <Real>                 the final time for the simulation (default 100)");
    eprintln!("  --nout <int>                the number of output times (default 100)");
    eprintln!("  --count-orbits              use rootfinding to count the number of completed orbits");
    eprintln!("  --check-order               compute the order of the method used and check if it is within range of the expected");
}

fn ParseArgs(argv: &[String], args: &mut ProgramArgs) -> i32 {
    args.step_mode = 0;
    args.stepper = 0;
    args.method_name = String::new();
    args.count_orbits = 0;
    args.use_compsums = 0;
    args.use_tstop = 1;
    args.dt = 1.0e-2;
    args.tf = 100.0;
    args.check_order = 0;
    args.num_output_times = 50;

    let mut argi = 0;
    while argi + 1 < argv.len() {
        argi += 1;
        let arg = &argv[argi];
        match arg.as_str() {
            "--step-mode" => {
                argi += 1;
                match argv[argi].as_str() {
                    "fixed" => args.step_mode = 0,
                    "adapt" => args.step_mode = 1,
                    _ => {
                        eprintln!("ERROR: --step-mode must be 'fixed' or 'adapt'");
                        return 1;
                    }
                }
            }
            "--stepper" => {
                argi += 1;
                match argv[argi].as_str() {
                    "SPRK" => args.stepper = 0,
                    "ERK" => args.stepper = 1,
                    _ => {
                        eprintln!("ERROR: --stepper must be 'SPRK' or 'ERK'");
                        return 1;
                    }
                }
            }
            "--method" => {
                argi += 1;
                args.method_name = argv[argi].clone();
            }
            "--dt" => {
                argi += 1;
                args.dt = argv[argi].parse().unwrap_or(0.0);
            }
            "--tf" => {
                argi += 1;
                args.tf = argv[argi].parse().unwrap_or(0.0);
            }
            "--nout" => {
                argi += 1;
                args.num_output_times = argv[argi].parse().unwrap_or(0);
            }
            "--count-orbits" => args.count_orbits = 1,
            "--disable-tstop" => args.use_tstop = 0,
            "--use-compensated-sums" => args.use_compsums = 1,
            "--check-order" => args.check_order = 1,
            "--help" => {
                PrintHelp();
                return 1;
            }
            _ => {
                eprintln!("ERROR: unrecognized argument {}", arg);
                PrintHelp();
                return 1;
            }
        }
    }

    if args.method_name.is_empty() {
        if args.stepper == 0 {
            args.method_name = "ARKODE_SPRK_MCLACHLAN_4_4".to_string();
        } else if args.stepper == 1 {
            args.method_name = "ARKODE_ZONNEVELD_5_3_4".to_string();
        }
    }

    0
}

fn PrintArgs(args: &ProgramArgs) {
    println!("Problem Arguments:");
    println!("  stepper:              {}", args.stepper);
    println!("  step mode:            {}", args.step_mode);
    println!("  use tstop:            {}", args.use_tstop);
    println!("  use compensated sums: {}", args.use_compsums);
    println!("  dt:                   {}", fmt_g(args.dt, 0, 6));
    println!("  Tf:                   {}", fmt_g(args.tf, 0, 6));
    println!("  nout:                 {}\n", args.num_output_times);
}

#[allow(clippy::too_many_arguments)] /* C helper signature kept */
fn ComputeConvergence(
    num_dt: usize,
    orders: &[f64],
    _expected_order: f64,
    a11: f64,
    a12: f64,
    a21: f64,
    a22: f64,
    b1: f64,
    b2: f64,
    ord_avg: &mut f64,
    ord_max: &mut f64,
    ord_est: &mut f64,
) -> i32 {
    /* Compute/print overall estimated convergence rate */
    *ord_avg = 0.0;
    *ord_max = 0.0;
    *ord_est = 0.0;
    for i in 1..num_dt {
        *ord_avg += orders[i - 1];
        *ord_max = SUNMAX(*ord_max, orders[i - 1]);
    }
    *ord_avg /= num_dt as f64 - 1.0;
    let det = a11 * a22 - a12 * a21;
    *ord_est = (a11 * b2 - a21 * b1) / det;
    let _ = (a12, a22);
    0
}

fn InitialConditions(y0vec: &mut NVector, ecc: f64) {
    let y0 = &mut y0vec.data;
    y0[0] = 1.0 - ecc;
    y0[1] = 0.0;
    y0[2] = 0.0;
    y0[3] = SUNRsqrt((1.0 + ecc) / (1.0 - ecc));
}

fn Hamiltonian(yvec: &NVector) -> f64 {
    let y = &yvec.data;
    let sqrt_q_t_q = SUNRsqrt(y[0] * y[0] + y[1] * y[1]);
    let p_t_p = y[2] * y[2] + y[3] * y[3];
    0.5 * p_t_p - 1.0 / sqrt_q_t_q
}

fn AngularMomentum(yvec: &NVector) -> f64 {
    let y = &yvec.data;
    let (q1, q2, p1, p2) = (y[0], y[1], y[2], y[3]);
    q1 * p2 - q2 * p1
}

fn dydt(t: f64, yvec: &NVector, ydotvec: &mut NVector, user_data: &mut UserData) -> i32 {
    let mut retval = 0;
    retval += force(t, yvec, ydotvec, user_data);
    retval += velocity(t, yvec, ydotvec, user_data);
    retval
}

fn velocity(_t: f64, yvec: &NVector, ydotvec: &mut NVector, _user_data: &mut UserData) -> i32 {
    let p1 = yvec.data[2];
    let p2 = yvec.data[3];
    ydotvec.data[0] = p1;
    ydotvec.data[1] = p2;
    0
}

fn force(_t: f64, yvec: &NVector, ydotvec: &mut NVector, _user_data: &mut UserData) -> i32 {
    let q1 = yvec.data[0];
    let q2 = yvec.data[1];
    let sqrt_q_t_q = SUNRsqrt(q1 * q1 + q2 * q2);
    ydotvec.data[2] = -q1 / SUNRpowerR(sqrt_q_t_q, 3.0);
    ydotvec.data[3] = -q2 / SUNRpowerR(sqrt_q_t_q, 3.0);
    0
}

fn rootfn(_t: f64, yvec: &NVector, gout: &mut [f64], _user_data: &mut UserData) -> i32 {
    gout[0] = yvec.data[1];
    0
}

/* N_VPrintFile_Serial: one "% .15e" value per line (SUN_FORMAT_E) */
fn n_v_print_file(x: &NVector, outfile: &mut std::fs::File) {
    for &xi in &x.data {
        let sp = if xi.is_sign_negative() { "" } else { " " };
        let _ = writeln!(outfile, "{}{}", sp, fmt_e(xi, 0, 15));
    }
}

fn SolveProblem(args: &ProgramArgs, result: &mut ProblemResult, ctx: &SUNContext) -> i32 {
    let count_orbits = args.count_orbits;
    let step_mode = args.step_mode;
    let stepper = args.stepper;
    let use_compsums = args.use_compsums;
    let num_output_times = args.num_output_times;
    let method_name = &args.method_name;
    let dt = args.dt;
    let tf = args.tf;

    /* Default problem parameters */
    let t0: f64 = 0.0;
    let dtout = (tf - t0) / (num_output_times as f64);
    let ecc: f64 = 0.6;

    println!("\n   Begin Kepler Problem\n");
    PrintArgs(args);

    /* Allocate our state vector; fill the initial conditions */
    let mut y = N_VNew_Serial(4, ctx);
    InitialConditions(&mut y, ecc);

    let udata = KeplerUserData { ecc };

    /* Create the integrator */
    let mut num_orbits: f64 = 0.0;
    let mut rootsfound = [0i32; 1];
    let mut arkode_mem = if stepper == 0 {
        let mut am = SPRKStepCreate(Some(force), Some(velocity), t0, &y, ctx)
            .expect("SPRKStepCreate");

        /* Optional: enable temporal root-finding */
        if count_orbits != 0 {
            let retval = ARKodeRootInit(&mut am, 1, Some(rootfn));
            assert!(retval >= 0, "ARKodeRootInit failed: {}", retval);
        }

        let retval = SPRKStepSetMethodName(&mut am, method_name);
        assert!(retval >= 0, "SPRKStepSetMethodName failed: {}", retval);

        let retval = ARKodeSetUseCompensatedSums(&mut am, use_compsums != 0);
        assert!(retval >= 0, "ARKodeSetUseCompensatedSums failed: {}", retval);

        if step_mode == 0 {
            let retval = ARKodeSetFixedStep(&mut am, dt);
            assert!(retval >= 0, "ARKodeSetFixedStep failed: {}", retval);

            let retval = ARKodeSetMaxNumSteps(&mut am, ((tf / dt).ceil() as i64) + 1);
            assert!(retval >= 0, "ARKodeSetMaxNumSteps failed: {}", retval);
        } else {
            eprintln!("ERROR: adaptive time-steps are not supported with SPRKStep");
            return 1;
        }

        let retval = ARKodeSetUserData(&mut am, Some(Box::new(udata)));
        assert!(retval >= 0, "ARKodeSetUserData failed: {}", retval);

        am
    } else {
        let mut am = ARKStepCreate(Some(dydt), None, t0, &y, ctx).expect("ARKStepCreate");

        let retval = ARKStepSetTableName(&mut am, "ARKODE_DIRK_NONE", method_name);
        assert!(retval >= 0, "ARKStepSetTableName failed: {}", retval);

        if count_orbits != 0 {
            let retval = ARKodeRootInit(&mut am, 1, Some(rootfn));
            assert!(retval >= 0, "ARKodeRootInit failed: {}", retval);
        }

        let retval = ARKodeSetUserData(&mut am, Some(Box::new(udata)));
        assert!(retval >= 0, "ARKodeSetUserData failed: {}", retval);

        let retval = ARKodeSetMaxNumSteps(&mut am, ((tf / dt).ceil() as i64) + 1);
        assert!(retval >= 0, "ARKodeSetMaxNumSteps failed: {}", retval);

        if step_mode == 0 {
            let retval = ARKodeSetFixedStep(&mut am, dt);
            assert!(retval >= 0, "ARKodeSetFixedStep failed: {}", retval);
        } else {
            let retval = ARKodeSStolerances(&mut am, dt, dt);
            assert!(retval >= 0, "ARKodeSStolerances failed: {}", retval);
        }

        am
    };

    /* Open output files */
    let mut conserved_fp = std::fs::File::create(format!(
        "ark_kepler_conserved_{}-dt-{}.txt",
        method_name,
        fmt_e(dt, 0, 2)
    ))
    .unwrap();
    let mut solution_fp = std::fs::File::create(format!(
        "ark_kepler_solution_{}-dt-{}.txt",
        method_name,
        fmt_e(dt, 0, 2)
    ))
    .unwrap();
    let mut times_fp = std::fs::File::create(format!(
        "ark_kepler_times_{}-dt-{}.txt",
        method_name,
        fmt_e(dt, 0, 2)
    ))
    .unwrap();

    /* Print out starting energy, momentum before integrating */
    let mut tret = t0;
    let mut tout = t0 + dtout;
    let h0 = Hamiltonian(&y);
    let l0 = AngularMomentum(&y);
    println!(
        "t = {}, H(p,q) = {}, L(p,q) = {}",
        fmt_f(tret, 0, 4),
        fmt_f(h0, 0, 16),
        fmt_f(l0, 0, 16)
    );
    let _ = writeln!(times_fp, "{}", fmt_f(tret, 0, 16));
    let _ = writeln!(conserved_fp, "{}, {}", fmt_f(h0, 0, 16), fmt_f(l0, 0, 16));
    n_v_print_file(&y, &mut solution_fp);

    /* Do integration */
    let mut iout = 0;
    while iout < num_output_times {
        /* Optional: if the stop time is not set, then its possible that
           the exact requested output time will not be hit (even with a
           fixed time-step due to roundoff error accumulation) and
           interpolation will be used to get the solution at the output
           time. */
        if args.use_tstop != 0 {
            ARKodeSetStopTime(&mut arkode_mem, tout);
        }
        let retval = ARKodeEvolve(&mut arkode_mem, tout, &mut y, &mut tret, ARK_NORMAL);

        if retval == ARK_ROOT_RETURN {
            num_orbits += 0.5;

            print!("ROOT RETURN:\t");
            ARKodeGetRootInfo(&mut arkode_mem, &mut rootsfound);
            println!(
                "  g[0] = {:3}, y[0] = {:>3}, y[1] = {:>3}, num. orbits is now {}",
                rootsfound[0],
                fmt_g(y.data[0], 3, 6),
                fmt_g(y.data[1], 3, 6),
                fmt_f(num_orbits, 0, 2)
            );
            println!(
                "t = {}, H(p,q)-H0 = {}, L(p,q)-L0 = {}",
                fmt_f(tret, 0, 4),
                fmt_e(Hamiltonian(&y) - h0, 0, 16),
                fmt_e(AngularMomentum(&y) - l0, 0, 16)
            );
        } else if retval >= 0 {
            /* Output current integration status */
            println!(
                "t = {}, H(p,q)-H0 = {}, L(p,q)-L0 = {}",
                fmt_f(tret, 0, 4),
                fmt_e(Hamiltonian(&y) - h0, 0, 16),
                fmt_e(AngularMomentum(&y) - l0, 0, 16)
            );
            let _ = writeln!(times_fp, "{}", fmt_f(tret, 0, 16));
            let _ = writeln!(
                conserved_fp,
                "{}, {}",
                fmt_f(Hamiltonian(&y), 0, 16),
                fmt_f(AngularMomentum(&y), 0, 16)
            );

            n_v_print_file(&y, &mut solution_fp);

            tout += dtout;
            tout = if tout > tf { tf } else { tout };
            iout += 1;
        } else {
            eprintln!("Solver failure, stopping integration");
            break;
        }
    }

    /* Copy results */
    result.sol.data.copy_from_slice(&y.data);
    result.energy_error = Hamiltonian(&y) - h0;

    drop(times_fp);
    drop(conserved_fp);
    drop(solution_fp);
    drop(y);
    let mut stdout = std::io::stdout();
    ARKodePrintAllStats(&mut arkode_mem, &mut stdout, SUN_OUTPUTFORMAT_TABLE);
    let mut slot = Some(arkode_mem);
    ARKodeFree(&mut slot);
    0
}

fn main() {
    /* Create the SUNDIALS context object for this simulation */
    let sunctx = SUNContext_Create();

    /* Parse the command line arguments */
    let argv: Vec<String> = std::env::args().collect();
    let mut args = ProgramArgs {
        step_mode: 0,
        stepper: 0,
        num_output_times: 0,
        use_compsums: 0,
        use_tstop: 0,
        count_orbits: 0,
        check_order: 0,
        dt: 0.0,
        tf: 0.0,
        method_name: String::new(),
    };
    if ParseArgs(&argv, &mut args) != 0 {
        std::process::exit(1);
    }

    /* Allocate space for result variables */
    let mut result = ProblemResult {
        sol: N_VNew_Serial(4, &sunctx),
        energy_error: 0.0,
    };

    if args.check_order == 0 {
        /* SolveProblem calls a stepper to evolve the problem to Tf */
        let retval = SolveProblem(&args, &mut result, &sunctx);
        if retval != 0 {
            std::process::exit(1);
        }
    } else {
        /* Compute the order of accuracy of the method by testing
           it with different step sizes. */
        let mut acc_orders = [0.0f64; NUM_DT];
        let mut con_orders = [0.0f64; NUM_DT];
        let mut acc_errors = [0.0f64; NUM_DT];
        let mut con_errors = [0.0f64; NUM_DT];
        let method = ARKodeSPRKTable_LoadByName(&args.method_name).unwrap();
        let expected_order = method.q;
        let mut ref_sol = N_VClone(&result.sol);
        let mut error = N_VClone(&result.sol);
        let (mut a11, mut a12, mut a21, mut a22) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let (mut b1, mut b2, mut b1e, mut b2e) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let (mut ord_max_acc, mut ord_max_conv, mut ord_avg, mut ord_est) =
            (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let refine: f64 = 0.5;
        let dt: f64 = if expected_order >= 3 { 1e-1 } else { 1e-3 };
        let mut dts = [0.0f64; NUM_DT];

        /* Create a reference solution using 8th order ERK with a small
           time step */
        let old_step_mode = args.step_mode;
        let old_stepper = args.stepper;
        let old_method_name = args.method_name.clone();
        args.dt = 1e-3;
        args.step_mode = 0;
        args.stepper = 1;
        args.method_name = "ARKODE_ARK548L2SAb_ERK_8_4_5".to_string();

        /* Free method, we just needed it to get its order */
        ARKodeSPRKTable_Free(method);

        /* SolveProblem calls a stepper to evolve the problem to Tf */
        let retval = SolveProblem(&args, &mut result, &sunctx);
        if retval != 0 {
            std::process::exit(1);
        }

        /* Store the reference solution */
        ref_sol.data.copy_from_slice(&result.sol.data);

        /* Restore the program args */
        args.step_mode = old_step_mode;
        args.stepper = old_stepper;
        args.method_name = old_method_name;

        for (i, d) in dts.iter_mut().enumerate() {
            *d = dt * SUNRpowerR(refine, i as f64);
        }

        /* Compute the error with various step sizes */
        for i in 0..NUM_DT {
            /* Set the dt to use for this solve */
            args.dt = dts[i];

            /* SolveProblem calls a stepper to evolve the problem to Tf */
            let retval = SolveProblem(&args, &mut result, &sunctx);
            if retval != 0 {
                std::process::exit(1);
            }

            println!();

            /* Compute the error */
            {
                let ProblemResult { sol, .. } = &result;
                N_VLinearSum(1.0, sol, -1.0, &ref_sol, &mut error);
            }
            acc_errors[i] =
                SUNRsqrt(N_VDotProd(&error, &error)) / (N_VGetLength(&error) as f64);
            con_errors[i] = SUNRabs(result.energy_error);

            a11 += 1.0;
            a12 += dts[i].ln();
            a21 += dts[i].ln();
            a22 += dts[i].ln() * dts[i].ln();
            b1 += acc_errors[i].ln();
            b2 += acc_errors[i].ln() * dts[i].ln();
            b1e += con_errors[i].ln();
            b2e += con_errors[i].ln() * dts[i].ln();

            if i >= 1 {
                acc_orders[i - 1] =
                    (acc_errors[i] / acc_errors[i - 1]).ln() / (dts[i] / dts[i - 1]).ln();
                con_orders[i - 1] =
                    (con_errors[i] / con_errors[i - 1]).ln() / (dts[i] / dts[i - 1]).ln();
            }
        }

        /* Compute the order of accuracy */
        ComputeConvergence(
            NUM_DT,
            &acc_orders,
            expected_order as f64,
            a11,
            a12,
            a21,
            a22,
            b1,
            b2,
            &mut ord_avg,
            &mut ord_max_acc,
            &mut ord_est,
        );
        println!(
            "Order of accuracy wrt solution:    expected = {}, max = {},  avg = {},  overall = {}",
            expected_order,
            fmt_f(ord_max_acc, 0, 4),
            fmt_f(ord_avg, 0, 4),
            fmt_f(ord_est, 0, 4)
        );

        /* Compute the order of accuracy with respect to conservation */
        ComputeConvergence(
            NUM_DT,
            &con_orders,
            expected_order as f64,
            a11,
            a12,
            a21,
            a22,
            b1e,
            b2e,
            &mut ord_avg,
            &mut ord_max_conv,
            &mut ord_est,
        );

        println!(
            "Order of accuracy wrt Hamiltonian: expected = {}, max = {},  avg = {},  overall = {}",
            expected_order,
            fmt_f(ord_max_conv, 0, 4),
            fmt_f(ord_avg, 0, 4),
            fmt_f(ord_est, 0, 4)
        );

        if ord_max_acc < (expected_order as f64 - 0.5) {
            println!(
                ">>> FAILURE: computed order of accuracy wrt solution is below expected ({})",
                expected_order
            );
            std::process::exit(1);
        }

        if ord_max_conv < (expected_order as f64 - 0.5) {
            println!(
                ">>> FAILURE: computed order of accuracy wrt Hamiltonian is below expected ({})",
                expected_order
            );
            std::process::exit(1);
        }

        drop(ref_sol);
        drop(error);
    }
}
