/* ----------------------------------------------------------------------------
 * Translation of examples/arkode/C_serial/ark_damped_harmonic_symplectic.c
 * (+ ark_damped_harmonic_symplectic.h) (SUNDIALS 7.7.0).
 *
 * In this example we consider the time-dependent damped harmonic oscillator
 *    q'(t) = p(t) exp(-F(t))
 *    p'(t) = -(F(t) * p + omega^2(t) * q)
 * With the initial conditions q(0) = 1, p(0) = 0.
 * The Hamiltonian for the system is
 *    H(p,q,t) = (p^2 * exp(-F(t)))/2 + (omega^2(t) * q^2 * exp(F(t)))/2
 * where omega(t) = cos(t/2), F(t) = 0.018*sin(t/pi).
 * We simulate the problem on t = [0, 30] using the symplectic methods in
 * SPRKStep.
 *
 * This is example 7.2 from:
 * Struckmeier, J., & Riedel, C. (2002). Canonical transformations and exact
 * invariants for time-dependent Hamiltonian systems. Annalen der Physik,
 * 11(1), 15-38.
 *
 * The example has the following command line arguments:
 *   --order <int>               the order of the method to use (default 4)
 *   --dt <Real>                 the fixed-time step size to use (default 0.01)
 *   --nout <int>                the number of output times (default 100)
 *   --disable-tstop             turns off tstop mode
 *   --use-compensated-sums      turns on compensated summation in ARKODE where
 *                               applicable
 * --------------------------------------------------------------------------*/

use arkode_rs::arkode::{ARKodeEvolve, ARKodeFree};
use arkode_rs::arkode_io::{
    ARKodePrintAllStats, ARKodeSetFixedStep, ARKodeSetMaxNumSteps, ARKodeSetOrder,
    ARKodeSetStopTime, ARKodeSetUseCompensatedSums,
};
use arkode_rs::arkode_sprkstep::SPRKStepCreate;
use arkode_rs::sundials_utils::fmt_f;
use arkode_rs::*;

#[allow(clippy::approx_constant)]
const PI: f64 = 3.14159265358979323846264338327950;

struct ProgramArgs {
    order: i32,
    num_output_times: i32,
    use_compsums: i32,
    use_tstop: i32,
    tf: f64,
    dt: f64,
}

fn print_help() {
    eprintln!(
        "ark_damped_harmonic_symplectic: an ARKODE example demonstrating the SPRKStep time-stepping module solving a time-dependent damped harmonic oscillator"
    );
    eprintln!("  --order <int>               the order of the method to use (default 4)");
    eprintln!("  --dt <Real>                 the fixed-time step size to use (default 0.01)");
    eprintln!("  --nout <int>                the number of output times (default 100)");
    eprintln!(
        "  --use-compensated-sums      turns on compensated summation in ARKODE where applicable"
    );
    eprintln!("  --disable-tstop             turns off tstop mode");
}

fn parse_args(argv: &[String], args: &mut ProgramArgs) -> bool {
    let mut argi = 1;
    while argi < argv.len() {
        match argv[argi].as_str() {
            "--order" => {
                argi += 1;
                args.order = argv[argi].parse().unwrap_or(0);
            }
            "--tf" => {
                argi += 1;
                args.tf = argv[argi].parse().unwrap_or(0.0);
            }
            "--dt" => {
                argi += 1;
                args.dt = argv[argi].parse().unwrap_or(0.0);
            }
            "--nout" => {
                argi += 1;
                args.num_output_times = argv[argi].parse().unwrap_or(0);
            }
            "--use-compensated-sums" => {
                args.use_compsums = 1;
            }
            "--disable-tstop" => {
                args.use_tstop = 0;
            }
            "--help" => {
                print_help();
                return true;
            }
            other => {
                eprintln!("ERROR: unrecognized argument {}", other);
                print_help();
                return true;
            }
        }
        argi += 1;
    }
    false
}

/* Check if a SUNDIALS function returned a negative value */
fn check_retval(retval: i32, funcname: &str) -> bool {
    if retval < 0 {
        eprintln!("\nERROR: {}() failed with retval = {}\n", funcname, retval);
        return true;
    }
    false
}

fn omega(t: f64) -> f64 {
    (t / 2.0).cos()
}

fn big_f(t: f64) -> f64 {
    0.018 * (t / PI).sin()
}

fn hamiltonian(yvec: &NVector, t: f64) -> f64 {
    let p = yvec.data[0];
    let q = yvec.data[1];

    (p * p * (-big_f(t)).exp()) / 2.0 + (omega(t) * omega(t) * q * q * big_f(t).exp()) / 2.0
}

fn qdot(t: f64, y: &NVector, ydot: &mut NVector, _user_data: &mut UserData) -> i32 {
    let p = y.data[0];

    ydot.data[1] = p * (-big_f(t)).exp();

    0
}

fn pdot(t: f64, y: &NVector, ydot: &mut NVector, _user_data: &mut UserData) -> i32 {
    let p = y.data[0];
    let q = y.data[1];

    ydot.data[0] = -(big_f(t) * p + omega(t) * omega(t) * q);

    0
}

fn main() {
    let t0: f64 = 0.0;

    /* Parse the command line arguments */
    let argv: Vec<String> = std::env::args().collect();
    let mut args = ProgramArgs {
        order: 4,
        num_output_times: 8,
        use_compsums: 0,
        use_tstop: 1,
        tf: 10.0 * PI,
        dt: 1e-3,
    };
    if parse_args(&argv, &mut args) {
        std::process::exit(1);
    }

    /* Default integrator options */
    let order = args.order;
    let use_compsums = args.use_compsums;
    let num_output_times = args.num_output_times;

    /* Default problem parameters */
    let tf = args.tf;
    let dt = args.dt;
    let d_tout = (tf - t0) / (num_output_times as f64);

    /* Create the SUNDIALS context object for this simulation */
    let ctx = SUNContext_Create();

    println!("\n   Begin time-dependent damped harmonic oscillator problem\n");

    /* Allocate our state vector */
    let mut y = N_VNew_Serial(2, &ctx);

    /* Fill the initial conditions */
    y.data[0] = 0.0; /* \dot{q} = p */
    y.data[1] = 1.0; /* \ddot{q} = \dot{p} */

    /* Create SPRKStep integrator */
    let mut arkode_mem =
        SPRKStepCreate(Some(qdot), Some(pdot), t0, &y, &ctx).expect("SPRKStepCreate");

    let retval = ARKodeSetOrder(&mut arkode_mem, order);
    if check_retval(retval, "ARKodeSetOrder") {
        return;
    }

    let retval = ARKodeSetUseCompensatedSums(&mut arkode_mem, use_compsums != 0);
    if check_retval(retval, "ARKodeSetUseCompensatedSums") {
        return;
    }

    let retval = ARKodeSetFixedStep(&mut arkode_mem, dt);
    if check_retval(retval, "ARKodeSetFixedStep") {
        return;
    }

    let retval = ARKodeSetMaxNumSteps(&mut arkode_mem, ((tf / dt).ceil() as i64) + 2);
    if check_retval(retval, "ARKodeSetMaxNumSteps") {
        return;
    }

    /* Print out starting Hamiltonian before integrating */
    let mut tret = t0;
    let mut tout = t0 + d_tout;
    /* Output current integration status */
    println!(
        "t = {}, q(t) = {}, H = {}",
        fmt_f(tret, 0, 6),
        fmt_f(y.data[1], 0, 6),
        fmt_f(hamiltonian(&y, tret), 0, 6)
    );

    /* Do integration */
    for _iout in 0..num_output_times {
        if args.use_tstop != 0 {
            ARKodeSetStopTime(&mut arkode_mem, tout);
        }
        let retval = ARKodeEvolve(&mut arkode_mem, tout, &mut y, &mut tret, ARK_NORMAL);

        /* Output current integration status */
        println!(
            "t = {}, q(t) = {}, H = {}",
            fmt_f(tret, 0, 6),
            fmt_f(y.data[1], 0, 6),
            fmt_f(hamiltonian(&y, tret), 0, 6)
        );

        /* Check if the solve was successful, if so, update the time and
           continue */
        if retval >= 0 {
            tout += d_tout;
            tout = if tout > tf { tf } else { tout };
        } else {
            eprintln!("Solver failure, stopping integration");
            break;
        }
    }

    println!();
    let mut stdout = std::io::stdout();
    let _ = ARKodePrintAllStats(&mut arkode_mem, &mut stdout, SUN_OUTPUTFORMAT_TABLE);
    drop(y);
    let mut slot = Some(arkode_mem);
    ARKodeFree(&mut slot);
}
