/* -----------------------------------------------------------------
 * Translation of examples/arkode/C_serial/ark_analytic_partitioned.c
 * (SUNDIALS 7.7.0).
 *
 * We consider the initial value problem
 *    y' + lambda*y = y^2, y(0) = 1
 * proposed in
 *
 * Estep, D., et al. "An a posteriori-a priori analysis of multiscale operator
 * splitting." SIAM Journal on Numerical Analysis 46.3 (2008): 1116-1146.
 *
 * The parameter lambda is positive, t is in [0, 1], and the exact solution is
 *
 *    y(t) = lambda*y / (y(0) - (y(0) - lambda)*exp(lambda*t))
 *
 * This program solves the problem with a splitting or forcing method which can
 * be specified with the command line syntax
 *
 * ./ark_analytic_partitioned <integrator> <coefficients>
 *    integrator: either 'splitting' or 'forcing'
 *    coefficients (splitting only): the SplittingStepCoefficients to load
 *
 * The linear term lambda*y and nonlinear term y^2 are treated as the two
 * partitions. The former is integrated using a time step of 5e-3, while the
 * later uses a time step of 1e-3. The overall splitting or forcing integrator
 * uses a time step of 1e-2. Once solved, the program prints the error and
 * statistics.
 *
 * Note on ownership: the C example retains its own linear_mem/
 * nonlinear_mem pointers alongside the SUNSteppers; in this port the
 * steppers own the wrapped integrators (and the outer integrator owns
 * the steppers), so the inner statistics are printed by borrowing the
 * integrators back out of the outer step memory.
 * -----------------------------------------------------------------*/

use arkode_rs::arkode::{ARKodeEvolve, ARKodeFree};
use arkode_rs::arkode_arkstep::ARKStepCreate;
use arkode_rs::arkode_erkstep::ERKStepCreate;
use arkode_rs::arkode_forcingstep::ForcingStepCreate;
use arkode_rs::arkode_forcingstep_impl::ARKodeForcingStepMem;
use arkode_rs::arkode_io::{ARKodePrintAllStats, ARKodeSetFixedStep, ARKodeSetUserData};
use arkode_rs::arkode_splittingstep::{SplittingStepCreate, SplittingStepSetCoefficients};
use arkode_rs::arkode_splittingstep_coefficients::{
    SplittingStepCoefficients_Destroy, SplittingStepCoefficients_LoadCoefficientsByName,
};
use arkode_rs::arkode_splittingstep_impl::ARKodeSplittingStepMem;
use arkode_rs::arkode_sunstepper::ARKodeCreateSUNStepper;
use arkode_rs::sundials_utils::fmt_g;
use arkode_rs::*;

#[derive(Clone, Copy)]
struct PartitionedUserData {
    lambda: f64,
}

/* RHS for f^1(t, y) = -lambda * y */
fn f_linear(_t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let lambda = user_data
        .as_mut()
        .unwrap()
        .downcast_mut::<PartitionedUserData>()
        .unwrap()
        .lambda;
    N_VScale(-lambda, y, ydot);
    0
}

/* RHS for f^2(t, y) = y^2 */
fn f_nonlinear(_t: f64, y: &NVector, ydot: &mut NVector, _user_data: &mut UserData) -> i32 {
    N_VProd(y, y, ydot);
    0
}

/* Compute the exact analytic solution */
fn exact_sol(y0: &NVector, tf: f64, user_data: &PartitionedUserData) -> NVector {
    let mut sol = N_VClone(y0);
    let y0_val = y0.data[0];
    let lambda = user_data.lambda;
    sol.data[0] = lambda * y0_val / (y0_val - (y0_val - lambda) * SUNRexp(lambda * tf));
    sol
}

fn main() {
    /* Parse arguments */
    let args: Vec<String> = std::env::args().collect();
    let integrator_name: &str = if args.len() > 1 { &args[1] } else { "splitting" };
    if integrator_name != "splitting" && integrator_name != "forcing" {
        eprintln!(
            "Invalid integrator: {}\nMust be 'splitting' or 'forcing'",
            integrator_name
        );
        std::process::exit(1);
    }
    let coefficients_name: Option<&String> = if args.len() > 2 { Some(&args[2]) } else { None };

    /* Problem parameters */
    let t0: f64 = 0.0; /* initial time */
    let tf: f64 = 1.0; /* final time */
    let dt: f64 = 0.01; /* outer time step */
    let dt_linear = dt / 5.0; /* linear integrator time step */
    let dt_nonlinear = dt / 10.0; /* nonlinear integrator time step */

    let user_data = PartitionedUserData { lambda: 2.0 };

    /* Create the SUNDIALS context object for this simulation */
    let ctx = SUNContext_Create();

    /* Initialize vector with initial condition */
    let mut y = N_VNew_Serial(1, &ctx);
    N_VConst(1.0, &mut y);

    let y_exact = exact_sol(&y, tf, &user_data);

    println!("\nAnalytical ODE test problem:");
    println!("   integrator = {} method", integrator_name);
    if let Some(name) = coefficients_name {
        println!("   coefficients = {}", name);
    }
    println!("   lambda     = {}", fmt_g(user_data.lambda, 0, 6));

    /* Create the integrator for the linear partition */
    let mut linear_mem = ERKStepCreate(f_linear, t0, &y, &ctx);

    let flag = ARKodeSetUserData(&mut linear_mem, Some(Box::new(user_data)));
    assert!(flag >= 0, "ARKodeSetUserData failed with flag = {}", flag);

    let flag = ARKodeSetFixedStep(&mut linear_mem, dt_linear);
    assert!(flag >= 0, "ARKodeSetFixedStep failed with flag = {}", flag);

    /* Create the integrator for the nonlinear partition */
    let mut nonlinear_mem =
        ARKStepCreate(Some(f_nonlinear), None, t0, &y, &ctx).expect("ARKStepCreate");

    let flag = ARKodeSetFixedStep(&mut nonlinear_mem, dt_nonlinear);
    assert!(flag >= 0, "ARKodeSetFixedStep failed with flag = {}", flag);

    /* Create SUNSteppers out of the integrators */
    let steppers = vec![
        ARKodeCreateSUNStepper(linear_mem),
        ARKodeCreateSUNStepper(nonlinear_mem),
    ];

    /* Create the outer integrator */
    let splitting = integrator_name == "splitting";
    let mut arkode_mem = if splitting {
        let mut arkode_mem =
            SplittingStepCreate(steppers, 2, t0, &y, &ctx).expect("SplittingStepCreate");

        if let Some(name) = coefficients_name {
            let mut coefficients = SplittingStepCoefficients_LoadCoefficientsByName(name);
            assert!(
                coefficients.is_some(),
                "SplittingStepCoefficients_LoadCoefficientsByName() failed - returned NULL pointer"
            );

            let flag = SplittingStepSetCoefficients(&mut arkode_mem, coefficients.as_deref().unwrap());
            assert!(flag >= 0, "SplittingStepSetCoefficients failed with flag = {}", flag);

            SplittingStepCoefficients_Destroy(&mut coefficients);
        }
        arkode_mem
    } else {
        let mut steppers_iter = steppers.into_iter();
        let stepper1 = steppers_iter.next().unwrap();
        let stepper2 = steppers_iter.next().unwrap();
        ForcingStepCreate(stepper1, stepper2, t0, &y, &ctx).expect("ForcingStepCreate")
    };

    let flag = ARKodeSetFixedStep(&mut arkode_mem, dt);
    assert!(flag >= 0, "ARKodeSetFixedStep failed with flag = {}", flag);

    /* Compute the numerical solution */
    let mut tret: f64 = 0.0;
    let flag = ARKodeEvolve(&mut arkode_mem, tf, &mut y, &mut tret, ARK_NORMAL);
    assert!(flag >= 0, "ARKodeEvolve failed with flag = {}", flag);

    /* Print the numerical error and statistics */
    let mut y_err = N_VClone(&y);
    N_VLinearSum(1.0, &y, -1.0, &y_exact, &mut y_err);
    println!("\nError: {}", fmt_g(N_VMaxNorm(&y_err), 0, 6));

    let mut stdout = std::io::stdout();
    println!("\nSplitting Stepper Statistics:");
    let flag = ARKodePrintAllStats(&mut arkode_mem, &mut stdout, SUN_OUTPUTFORMAT_TABLE);
    assert!(flag >= 0, "ARKodePrintAllStats failed with flag = {}", flag);

    /* borrow the inner integrators back out of the outer step memory
    (C keeps separate linear_mem/nonlinear_mem pointers) */
    if splitting {
        let step_mem = arkode_mem
            .step_mem
            .as_mut()
            .unwrap()
            .downcast_mut::<ARKodeSplittingStepMem>()
            .unwrap();

        println!("\nLinear Stepper Statistics:");
        let linear_mem = step_mem.steppers[0]
            .content
            .as_mut()
            .unwrap()
            .downcast_mut::<ARKodeMem>()
            .unwrap();
        let flag = ARKodePrintAllStats(linear_mem, &mut stdout, SUN_OUTPUTFORMAT_TABLE);
        assert!(flag >= 0, "ARKodePrintAllStats failed with flag = {}", flag);

        println!("\nNonlinear Stepper Statistics:");
        let nonlinear_mem = step_mem.steppers[1]
            .content
            .as_mut()
            .unwrap()
            .downcast_mut::<ARKodeMem>()
            .unwrap();
        let flag = ARKodePrintAllStats(nonlinear_mem, &mut stdout, SUN_OUTPUTFORMAT_TABLE);
        assert!(flag >= 0, "ARKodePrintAllStats failed with flag = {}", flag);
    } else {
        let step_mem = arkode_mem
            .step_mem
            .as_mut()
            .unwrap()
            .downcast_mut::<ARKodeForcingStepMem>()
            .unwrap();

        println!("\nLinear Stepper Statistics:");
        let linear_mem = step_mem.stepper[0]
            .content
            .as_mut()
            .unwrap()
            .downcast_mut::<ARKodeMem>()
            .unwrap();
        let flag = ARKodePrintAllStats(linear_mem, &mut stdout, SUN_OUTPUTFORMAT_TABLE);
        assert!(flag >= 0, "ARKodePrintAllStats failed with flag = {}", flag);

        println!("\nNonlinear Stepper Statistics:");
        let nonlinear_mem = step_mem.stepper[1]
            .content
            .as_mut()
            .unwrap()
            .downcast_mut::<ARKodeMem>()
            .unwrap();
        let flag = ARKodePrintAllStats(nonlinear_mem, &mut stdout, SUN_OUTPUTFORMAT_TABLE);
        assert!(flag >= 0, "ARKodePrintAllStats failed with flag = {}", flag);
    }

    /* Free memory (the outer integrator owns the steppers and inner
    integrators; dropping it frees everything) */
    drop(y);
    drop(y_exact);
    drop(y_err);
    let mut slot = Some(arkode_mem);
    ARKodeFree(&mut slot);
}
