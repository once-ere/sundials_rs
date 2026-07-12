/* ------------------------------------------------------------------
 * Translation of examples/arkode/C_serial/ark_reaction_diffusion_mri.c
 * (SUNDIALS 7.7.0).
 *
 * Example problem:
 *
 * The following test simulates a simple 1D reaction-diffusion
 * equation,
 *
 *   y_t = k * y_xx + y^2 * (1-y)
 *
 * for t in [0, 3], x in [0, L] with boundary conditions,
 *
 *   y_x(0,t) = y_x(L,t) = 0
 *
 * and initial condition,
 *
 *   y(x,0) = (1 + exp(lambda*(x-1))^(-1),
 *
 * with parameter k = 1e-4/ep, lambda = 0.5*sqrt(2*ep*1e4),
 * ep = 1e-2, and L = 5.
 *
 * The spatial derivatives are computed using second-order
 * centered differences, with the data distributed over N points
 * on a uniform spatial grid.
 *
 * This program solves the problem with the MRI stepper. Outputs are
 * printed at equal intervals of 0.1 and run statistics are printed
 * at the end.
 *
 * Note on ownership: the stepper owns the wrapped inner integrator
 * (and the outer integrator owns the stepper), so the final fast
 * statistics are printed by borrowing the inner integrator back out
 * of the outer step memory.  The shared C udata pointer becomes one
 * identical UserData clone per integrator (the problem data is
 * immutable during the run).
 * ----------------------------------------------------------------*/
#![allow(non_snake_case)]

use std::io::Write;

use arkode_rs::arkode::{ARKodeCreateMRIStepInnerStepper, ARKodeEvolve, ARKodeFree};
use arkode_rs::arkode_arkstep::ARKStepCreate;
use arkode_rs::arkode_arkstep_io::ARKStepSetTableNum;
use arkode_rs::arkode_butcher_erk::ARKODE_KNOTH_WOLKE_3_3;
use arkode_rs::arkode_io::{
    ARKodePrintAllStats, ARKodeSetFixedStep, ARKodeSetMaxNumSteps, ARKodeSetUserData,
};
use arkode_rs::arkode_mristep::MRIStepCreate;
use arkode_rs::arkode_mristep_impl::ARKodeMRIStepMem;
use arkode_rs::sundials_utils::{fmt_e, fmt_f, fmt_g};
use arkode_rs::*;

/* user data structure */
#[derive(Clone)]
struct UData {
    N: i64,  /* number of intervals   */
    dx: f64, /* mesh spacing          */
    k: f64,  /* diffusion coefficient */
    lam: f64,
}

/* ------------------------------
 * Functions called by the solver
 * ------------------------------*/

/* ff routine to compute the fast portion of the ODE RHS. */
fn ff(_t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let udata = user_data.as_ref().unwrap().downcast_ref::<UData>().unwrap();
    let N = udata.N; /* set variable shortcuts */

    /* iterate over domain, computing reaction term */
    for i in 0..N as usize {
        ydot.data[i] = y.data[i] * y.data[i] * (1.0 - y.data[i]);
    }

    /* Return with success */
    0
}

/* fs routine to compute the slow portion of the ODE RHS. */
fn fs(_t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let udata = user_data.as_ref().unwrap().downcast_ref::<UData>().unwrap();
    let N = udata.N; /* set variable shortcuts */
    let k = udata.k;
    let dx = udata.dx;

    /* iterate over domain, computing diffusion term */
    let c1 = k / dx / dx;
    let c2 = 2.0 * k / dx / dx;

    /* left boundary condition */
    ydot.data[0] = c2 * (y.data[1] - y.data[0]);

    /* interior points */
    for i in 1..(N - 1) as usize {
        ydot.data[i] = c1 * y.data[i - 1] - c2 * y.data[i] + c1 * y.data[i + 1];
    }

    /* right boundary condition */
    ydot.data[(N - 1) as usize] = c2 * (y.data[(N - 2) as usize] - y.data[(N - 1) as usize]);

    /* Return with success */
    0
}

/* -----------------------------------------
 * Private function to set initial condition
 * -----------------------------------------*/

fn SetInitialCondition(y: &mut NVector, udata: &UData) -> i32 {
    let N = udata.N; /* set variable shortcuts */
    let lam = udata.lam;
    let dx = udata.dx;

    /* set initial condition */
    for i in 0..N as usize {
        y.data[i] = 1.0 / (1.0 + (lam * (i as f64 * dx - 1.0)).exp());
    }

    /* Return with success */
    0
}

/* ------------------------------
 * Private helper functions
 * ------------------------------*/

/* Check if a SUNDIALS function returned a negative flag */
fn check_retval(retval: i32, funcname: &str) -> bool {
    if retval < 0 {
        eprintln!("\nSUNDIALS_ERROR: {}() failed with retval = {}\n", funcname, retval);
        return true;
    }
    false
}

/* Main Program */
fn main() {
    /* general problem parameters */
    let T0: f64 = 0.0; /* initial time */
    let Tf: f64 = 3.0; /* final time */
    let mut dTout: f64 = 0.1; /* time between outputs */
    let Nt: i32 = (Tf / dTout).ceil() as i32; /* number of output times */
    let hs: f64 = 0.001; /* slow step size */
    let hf: f64 = 0.00002; /* fast step size */

    let L: f64 = 5.0; /* domain length */
    let N: i64 = 1001; /* number of mesh points */
    let ep: f64 = 1e-2;

    /* Create the SUNDIALS context object for this simulation */
    let ctx = SUNContext_Create();

    /*
     * Initialization
     */

    /* allocate and fill user data structure */
    let udata = UData {
        N,
        dx: L / (1.0 * N as f64 - 1.0),
        k: 1e-4 / ep,
        lam: 0.5 * (2.0 * ep * 1e4).sqrt(),
    };

    /* Initial problem output */
    println!("\n1D reaction-diffusion PDE test problem:");
    println!("  N = {}", udata.N);
    println!("  diffusion coefficient:  k = {}", fmt_g(udata.k, 0, 6));

    /* Create and initialize serial vector for the solution */
    let mut y = N_VNew_Serial(N, &ctx);
    let retval = SetInitialCondition(&mut y, &udata);
    if check_retval(retval, "SetInitialCondition") {
        return;
    }

    /*
     * Create the fast integrator and set options
     */

    /* Initialize the fast integrator. Specify the explicit fast right-hand
       side function in y'=fe(t,y)+fi(t,y)+ff(t,y), the initial time T0, and
       the initial dependent variable vector y. */
    let mut inner_arkode_mem = ARKStepCreate(Some(ff), None, T0, &y, &ctx).expect("ARKStepCreate");

    /* Attach user data to fast integrator */
    let retval = ARKodeSetUserData(&mut inner_arkode_mem, Some(Box::new(udata.clone())));
    if check_retval(retval, "ARKodeSetUserData") {
        return;
    }

    /* Set the fast method */
    let retval = ARKStepSetTableNum(&mut inner_arkode_mem, -1, ARKODE_KNOTH_WOLKE_3_3);
    if check_retval(retval, "ARKStepSetTableNum") {
        return;
    }

    /* Set the fast step size */
    let retval = ARKodeSetFixedStep(&mut inner_arkode_mem, hf);
    if check_retval(retval, "ARKodeSetFixedStep") {
        return;
    }

    /* Create inner stepper */
    let inner_stepper =
        ARKodeCreateMRIStepInnerStepper(inner_arkode_mem).expect("ARKodeCreateMRIStepInnerStepper");

    /*
     * Create the slow integrator and set options
     */

    /* Initialize the slow integrator. Specify the explicit slow right-hand
       side function in y'=fe(t,y)+fi(t,y)+ff(t,y), the initial time T0, the
       initial dependent variable vector y, and the fast integrator. */
    let mut arkode_mem =
        MRIStepCreate(Some(fs), None, T0, &y, inner_stepper, &ctx).expect("MRIStepCreate");

    /* Pass udata to user functions */
    let retval = ARKodeSetUserData(&mut arkode_mem, Some(Box::new(udata.clone())));
    if check_retval(retval, "ARKodeSetUserData") {
        return;
    }

    /* Set the slow step size */
    let retval = ARKodeSetFixedStep(&mut arkode_mem, hs);
    if check_retval(retval, "ARKodeSetFixedStep") {
        return;
    }

    /* Increase max num steps */
    let retval = ARKodeSetMaxNumSteps(&mut arkode_mem, 10000);
    if check_retval(retval, "ARKodeSetMaxNumSteps") {
        return;
    }

    /*
     * Integrate ODE
     */

    /* output mesh to disk */
    let mut fid = std::fs::File::create("heat_mesh.txt").expect("fopen");
    for i in 0..N {
        let _ = writeln!(fid, "  {}", fmt_e(udata.dx * i as f64, 0, 16));
    }
    drop(fid);

    /* Open output stream for results */
    let mut ufid = std::fs::File::create("heat1D.txt").expect("fopen");

    /* output initial condition to disk */
    for i in 0..N as usize {
        let _ = write!(ufid, " {}", fmt_e(y.data[i], 0, 16));
    }
    let _ = writeln!(ufid);

    /* Main time-stepping loop: calls ARKodeEvolve to perform the integration,
       then prints results. Stops when the final time has been reached */
    let mut t = T0;
    dTout = (Tf - T0) / Nt as f64;
    let mut tout = T0 + dTout;
    println!("        t      ||u||_rms");
    println!("   -------------------------");
    println!(
        "  {}  {}",
        fmt_f(t, 10, 6),
        fmt_f((N_VDotProd(&y, &y) / N as f64).sqrt(), 10, 6)
    );
    for _iout in 0..Nt {
        /* call integrator */
        let retval = ARKodeEvolve(&mut arkode_mem, tout, &mut y, &mut t, ARK_NORMAL);
        if check_retval(retval, "ARKodeEvolve") {
            break;
        }

        /* print solution stats and output results to disk */
        println!(
            "  {}  {}",
            fmt_f(t, 10, 6),
            fmt_f((N_VDotProd(&y, &y) / N as f64).sqrt(), 10, 6)
        );
        for i in 0..N as usize {
            let _ = write!(ufid, " {}", fmt_e(y.data[i], 0, 16));
        }
        let _ = writeln!(ufid);

        /* successful solve: update output time */
        tout += dTout;
        tout = if tout > Tf { Tf } else { tout };
    }
    println!("   -------------------------");
    drop(ufid);

    /* Print final statistics to the screen (fast statistics via the
       inner integrator borrowed back out of the outer step memory) */
    let mut stdout = std::io::stdout();
    println!("\nFinal Slow Statistics:");
    let _ = ARKodePrintAllStats(&mut arkode_mem, &mut stdout, SUN_OUTPUTFORMAT_TABLE);
    println!("\nFinal Fast Statistics:");
    {
        let step_mem = arkode_mem
            .step_mem
            .as_mut()
            .unwrap()
            .downcast_mut::<ARKodeMRIStepMem>()
            .unwrap();
        let inner = step_mem
            .stepper
            .content
            .as_mut()
            .unwrap()
            .downcast_mut::<ARKodeMem>()
            .unwrap();
        let _ = ARKodePrintAllStats(inner, &mut stdout, SUN_OUTPUTFORMAT_TABLE);
    }

    /* Print final statistics to a file in CSV format */
    let mut fid =
        std::fs::File::create("ark_reaction_diffusion_mri_slow_stats.csv").expect("fopen");
    let _ = ARKodePrintAllStats(&mut arkode_mem, &mut fid, SUN_OUTPUTFORMAT_CSV);
    drop(fid);
    let mut fid =
        std::fs::File::create("ark_reaction_diffusion_mri_fast_stats.csv").expect("fopen");
    {
        let step_mem = arkode_mem
            .step_mem
            .as_mut()
            .unwrap()
            .downcast_mut::<ARKodeMRIStepMem>()
            .unwrap();
        let inner = step_mem
            .stepper
            .content
            .as_mut()
            .unwrap()
            .downcast_mut::<ARKodeMem>()
            .unwrap();
        let _ = ARKodePrintAllStats(inner, &mut fid, SUN_OUTPUTFORMAT_CSV);
    }
    drop(fid);

    /* Clean up and return */
    drop(y); /* Free y vector */
    let mut slot = Some(arkode_mem);
    ARKodeFree(&mut slot); /* Free integrator memory */
}
