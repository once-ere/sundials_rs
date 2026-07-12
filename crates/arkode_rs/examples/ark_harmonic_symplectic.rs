/* ----------------------------------------------------------------------------
 * Translation of examples/arkode/C_serial/ark_harmonic_symplectic.c
 * (+ ark_harmonic_symplectic.h) (SUNDIALS 7.7.0).
 *
 * In this example we consider the simple harmonic oscillator
 *    x''(t) + omega^2*x(t) = 0.
 * We rewrite the second order ODE as the first order ODE model
 *    x'(t) = v(t)
 *    v'(t) = -omega^2*x(t).
 * With the initial conditions x(0) = x0 and v(0) = v0,
 * the analytical solution is
 *    x(t) = A*cos(t*omega + phi),
 *    v(t) = -A*omega*sin(t*omega + phi)
 * where A = sqrt(x0^2 + v0^2/omega) and tan(phi) = v0/(omega*x0).
 * The total energy (potential + kinetic) in this system is
 *    E = (v^2 + omega^2*x^2) / 2
 * E is conserved and is the system Hamiltonian.
 * We simulate the problem on t = [0, 2pi] using the symplectic methods
 * in SPRKStep. Symplectic methods will approximately conserve E.
 *
 * The example has the following command line arguments:
 *   --order <int>               the order of the method to use (default 4)
 *   --dt <Real>                 the fixed-time step size to use (default 0.01)
 *   --nout <int>                the number of output times (default 100)
 *   --use-compensated-sums      turns on compensated summation in ARKODE where
 *                               applicable
 *   --disable-tstop             turns off tstop mode
 * --------------------------------------------------------------------------*/
#![allow(non_snake_case)]

use arkode_rs::arkode::{ARKodeEvolve, ARKodeFree};
use arkode_rs::arkode_io::{
    ARKodePrintAllStats, ARKodeSetFixedStep, ARKodeSetMaxNumSteps, ARKodeSetOrder,
    ARKodeSetStopTime, ARKodeSetUseCompensatedSums, ARKodeSetUserData,
};
use arkode_rs::arkode_sprkstep::SPRKStepCreate;
use arkode_rs::sundials_utils::{fmt_e, fmt_f};
use arkode_rs::*;

#[allow(clippy::approx_constant)]
const PI: f64 = 3.14159265358979323846264338327950;

#[derive(Clone)]
struct HarmonicData {
    A: f64,
    phi: f64,
    omega: f64,
}

struct ProgramArgs {
    order: i32,
    num_output_times: i32,
    use_compsums: i32,
    use_tstop: i32,
    tf: f64,
    dt: f64,
}

fn print_help() {
    let defaults = ProgramArgs {
        order: 4,
        num_output_times: 8,
        use_compsums: 0,
        use_tstop: 1,
        dt: 1e-3,
        tf: 2.0 * PI,
    };
    eprintln!(
        "ark_harmonic_symplectic: an ARKODE example demonstrating the SPRKStep time-stepping module solving a simple harmonic oscillator"
    );
    eprintln!(
        "  --order <int>               the order of the method to use (default {})",
        defaults.order
    );
    eprintln!(
        "  --dt <Real>                 the fixed-time step size to use (default {})",
        fmt_e(defaults.dt, 0, 1)
    );
    eprintln!(
        "  --nout <int>                the number of output times (default {})",
        defaults.num_output_times
    );
    eprintln!(
        "  --use-compensated-sums      turns on compensated summation in ARKODE where applicable"
    );
    eprintln!("  --disable-tstop             turns off tstop mode");
    let _ = defaults.use_compsums;
    let _ = defaults.use_tstop;
    let _ = defaults.tf;
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

fn solution(t: f64, solvec: &mut NVector, udata: &HarmonicData) {
    /* compute solution */
    solvec.data[0] = udata.A * (udata.omega * t + udata.phi).cos();
    solvec.data[1] = -udata.A * udata.omega * (udata.omega * t + udata.phi).sin();
}

fn energy(yvec: &NVector, _dt: f64, udata: &HarmonicData) -> f64 {
    let x = yvec.data[0];
    let v = yvec.data[1];
    let omega2 = udata.omega * udata.omega;

    (v * v + omega2 * x * x) / 2.0
}

fn xdot(_t: f64, y: &NVector, ydot: &mut NVector, _user_data: &mut UserData) -> i32 {
    let v = y.data[1];

    ydot.data[0] = v;

    0
}

fn vdot(_t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let udata = user_data.as_ref().unwrap().downcast_ref::<HarmonicData>().unwrap();
    let x = y.data[0];
    let omega2 = udata.omega * udata.omega;

    ydot.data[1] = -omega2 * x;

    0
}

fn main() {
    let t0: f64 = 0.0;
    let big_a: f64 = 10.0;
    let phi: f64 = 0.0;
    let omega: f64 = 1.0;

    /* Parse the command line arguments */
    let argv: Vec<String> = std::env::args().collect();
    let mut args = ProgramArgs {
        order: 4,
        num_output_times: 8,
        use_compsums: 0,
        use_tstop: 1,
        dt: 1e-3,
        tf: 2.0 * PI,
    };
    if parse_args(&argv, &mut args) {
        std::process::exit(1);
    }

    /* Default integrator options and problem parameters */
    let order = args.order;
    let use_compsums = args.use_compsums;
    let num_output_times = args.num_output_times;
    let tf = args.tf;
    let dt = args.dt;
    let d_tout = (tf - t0) / (num_output_times as f64);

    /* Create the SUNDIALS context object for this simulation */
    let ctx = SUNContext_Create();

    println!("\n   Begin simple harmonic oscillator problem\n");

    /* Allocate and fill udata structure */
    let udata = HarmonicData {
        A: big_a,
        phi,
        omega,
    };

    /* Allocate our state vector [x, v]^T */
    let mut y = N_VNew_Serial(2, &ctx);
    let mut sol = N_VClone(&y);

    /* Fill the initial conditions (x0 then v0) */
    y.data[0] = big_a * phi.cos();
    y.data[1] = -big_a * omega * phi.sin();

    /* Create SPRKStep integrator */
    let mut arkode_mem =
        SPRKStepCreate(Some(xdot), Some(vdot), t0, &y, &ctx).expect("SPRKStepCreate");

    let retval = ARKodeSetOrder(&mut arkode_mem, order);
    if check_retval(retval, "ARKodeSetOrder") {
        return;
    }

    let retval = ARKodeSetUserData(&mut arkode_mem, Some(Box::new(udata.clone())));
    if check_retval(retval, "ARKodeSetUserData") {
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

    /* Print out starting energy, momentum before integrating */
    let mut tret = t0;
    let mut tout = t0 + d_tout;
    println!(
        "t = {}, x(t) = {}, E = {}, sol. err = {}",
        fmt_f(tret, 0, 6),
        fmt_f(y.data[0], 0, 6),
        fmt_f(energy(&y, dt, &udata), 0, 6),
        fmt_f(0.0, 0, 6)
    );

    /* Do integration */
    for _iout in 0..num_output_times {
        if args.use_tstop != 0 {
            ARKodeSetStopTime(&mut arkode_mem, tout);
        }
        let retval = ARKodeEvolve(&mut arkode_mem, tout, &mut y, &mut tret, ARK_NORMAL);

        /* Compute the analytical solution */
        solution(tret, &mut sol, &udata);

        /* Compute L2 error: the C N_VLinearSum(1, y, -1, solution, solution)
           aliases its output with solution -> in-place form (bitwise equal
           to C's VDiff kernel) */
        sol.linear_sum_with(-1.0, 1.0, &y);
        let err = N_VDotProd(&sol, &sol).sqrt();

        /* Output current integration status */
        println!(
            "t = {}, x(t) = {}, E = {}, sol. err = {}",
            fmt_f(tret, 0, 6),
            fmt_f(y.data[0], 0, 6),
            fmt_f(energy(&y, dt, &udata), 0, 6),
            fmt_e(err, 0, 16)
        );

        /* Check that solution error is within tolerance */
        if err > SUNMAX(dt / (10.0f64).powi(order - 2), 1000.0 * SUN_UNIT_ROUNDOFF) {
            eprintln!("FAILURE: solution error is too high");
            std::process::exit(1);
        }

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
    drop(y);
    drop(sol);
    let mut stdout = std::io::stdout();
    let _ = ARKodePrintAllStats(&mut arkode_mem, &mut stdout, SUN_OUTPUTFORMAT_TABLE);
    let mut slot = Some(arkode_mem);
    ARKodeFree(&mut slot);
}
