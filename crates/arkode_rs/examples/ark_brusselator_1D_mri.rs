/* -----------------------------------------------------------------
 * Translation of examples/arkode/C_serial/ark_brusselator_1D_mri.c
 * (SUNDIALS 7.7.0).
 *
 * This program simulates a 1D advection-reaction problem. The
 * brusselator problem from chemical kinetics is used for the
 * reaction terms.
 *
 * This program uses the MRIStep module with an explicit slow
 * method and an implicit fast method. The explicit method uses a
 * fixed step size and the implicit method uses adaptive steps.
 * Implicit systems are solved using a Newton iteration with the
 * band linear solver, and a user-supplied Jacobian routine for the
 * fast RHS.
 *
 * 100 outputs are printed at equal intervals, and run statistics
 * are printed at the end.
 *
 * Note on ownership: the stepper owns the wrapped inner integrator
 * (see ark_brusselator_mri.rs); the fast-integrator statistics are
 * read by borrowing the inner integrator back out of the outer
 * step memory.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]

use arkode_rs::arkode::{
    ARKodeCreateMRIStepInnerStepper, ARKodeEvolve, ARKodeFree, ARKodeSStolerances,
};
use arkode_rs::arkode_arkstep::ARKStepCreate;
use arkode_rs::arkode_arkstep_io::ARKStepSetTableNum;
use arkode_rs::arkode_butcher_dirk::ARKODE_ARK324L2SA_DIRK_4_2_3;
use arkode_rs::arkode_io::{
    ARKodeGetNumErrTestFails, ARKodeGetNumRhsEvals, ARKodeGetNumStepAttempts, ARKodeGetNumSteps,
    ARKodeSetFixedStep, ARKodeSetUserData,
};
use arkode_rs::arkode_ls::{
    ARKodeGetNumJacEvals, ARKodeGetNumLinRhsEvals, ARKodeSetJacFn, ARKodeSetLinearSolver,
};
use arkode_rs::arkode_io::{ARKodeGetNumLinSolvSetups, ARKodeGetNumNonlinSolvConvFails,
    ARKodeGetNumNonlinSolvIters};
use arkode_rs::arkode_mristep::MRIStepCreate;
use arkode_rs::arkode_mristep_impl::ARKodeMRIStepMem;
use arkode_rs::sundials_utils::{fmt_e, fmt_f, fmt_g};
use arkode_rs::*;
use std::io::Write;

/* Define some constants */
const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;

/* accessor macro between (x,v) location and 1D NVector array */
fn IDX(x: usize, v: usize) -> usize {
    3 * x + v
}

/* user data structure */
#[derive(Clone, Copy)]
struct BrusselatorData {
    N: i64,  /* number of intervals      */
    dx: f64, /* mesh spacing             */
    a: f64,  /* constant forcing on u    */
    b: f64,  /* steady-state value of w  */
    c: f64,  /* advection coefficient    */
    ep: f64, /* stiffness parameter      */
}

/* ff routine to compute the fast portion of the ODE RHS. */
fn ff(_t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let udata = user_data
        .as_mut()
        .unwrap()
        .downcast_mut::<BrusselatorData>()
        .unwrap();
    let N = udata.N as usize;
    let a = udata.a;
    let b = udata.b;
    let ep = udata.ep;

    /* iterate over domain, computing reactions */
    for i in 0..N {
        /* set shortcuts */
        let u = y.data[IDX(i, 0)];
        let v = y.data[IDX(i, 1)];
        let w = y.data[IDX(i, 2)];

        /* u_t = a - (w+1)*u + v*u^2 */
        ydot.data[IDX(i, 0)] = a - (w + ONE) * u + v * u * u;

        /* v_t = w*u - v*u^2 */
        ydot.data[IDX(i, 1)] = w * u - v * u * u;

        /* w_t = (b-w)/ep - w*u */
        ydot.data[IDX(i, 2)] = (b - w) / ep - w * u;
    }

    /* return success */
    0
}

/* fs routine to compute the slow portion of the ODE RHS. */
fn fs(_t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let udata = user_data
        .as_mut()
        .unwrap()
        .downcast_mut::<BrusselatorData>()
        .unwrap();
    let N = udata.N as usize;
    let c = udata.c;
    let dx = udata.dx;

    /* iterate over domain, computing advection */
    let tmp = -c / dx;

    if c > ZERO {
        /*
         * right moving flow
         */

        /* left boundary Jacobian entries */
        ydot.data[IDX(0, 0)] = tmp * (y.data[IDX(0, 0)] - y.data[IDX(N - 1, 0)]);
        ydot.data[IDX(0, 1)] = tmp * (y.data[IDX(0, 1)] - y.data[IDX(N - 1, 1)]);
        ydot.data[IDX(0, 2)] = tmp * (y.data[IDX(0, 2)] - y.data[IDX(N - 1, 2)]);

        /* interior Jacobian entries */
        for i in 1..N {
            ydot.data[IDX(i, 0)] = tmp * (y.data[IDX(i, 0)] - y.data[IDX(i - 1, 0)]);
            ydot.data[IDX(i, 1)] = tmp * (y.data[IDX(i, 1)] - y.data[IDX(i - 1, 1)]);
            ydot.data[IDX(i, 2)] = tmp * (y.data[IDX(i, 2)] - y.data[IDX(i - 1, 2)]);
        }
    } else if c < ZERO {
        /*
         * left moving flow
         */

        /* interior Jacobian entries */
        for i in 0..(N - 1) {
            ydot.data[IDX(i, 0)] = tmp * (y.data[IDX(i + 1, 0)] - y.data[IDX(i, 0)]);
            ydot.data[IDX(i, 1)] = tmp * (y.data[IDX(i + 1, 1)] - y.data[IDX(i, 1)]);
            ydot.data[IDX(i, 2)] = tmp * (y.data[IDX(i + 1, 2)] - y.data[IDX(i, 2)]);
        }

        /* right boundary Jacobian entries */
        ydot.data[IDX(N - 1, 0)] = tmp * (y.data[IDX(N - 1, 0)] - y.data[IDX(0, 0)]);
        ydot.data[IDX(N - 1, 1)] = tmp * (y.data[IDX(N - 1, 1)] - y.data[IDX(0, 1)]);
        ydot.data[IDX(N - 1, 2)] = tmp * (y.data[IDX(N - 1, 2)] - y.data[IDX(0, 2)]);
    }

    /* return success */
    0
}

/* Jf routine to compute the Jacobian of the fast portion of the ODE RHS. */
fn Jf(
    _t: f64,
    y: &NVector,
    _fy: &NVector,
    Jac: &mut SUNMatrix,
    user_data: &mut UserData,
    _tmp1: &mut NVector,
    _tmp2: &mut NVector,
    _tmp3: &mut NVector,
) -> i32 {
    let udata = user_data
        .as_mut()
        .unwrap()
        .downcast_mut::<BrusselatorData>()
        .unwrap();
    let N = udata.N as usize;
    let ep = udata.ep;

    let jm = match Jac {
        SUNMatrix::Band(m) => m,
        _ => return 1,
    };

    /* iterate over nodes, filling in Jacobian entries */
    for i in 0..N {
        /* set nodal value shortcuts */
        let u = y.data[IDX(i, 0)];
        let v = y.data[IDX(i, 1)];
        let w = y.data[IDX(i, 2)];

        /* all vars wrt u */
        jm.set(IDX(i, 0) as i64, IDX(i, 0) as i64, TWO * u * v - (w + ONE));
        jm.set(IDX(i, 1) as i64, IDX(i, 0) as i64, w - TWO * u * v);
        jm.set(IDX(i, 2) as i64, IDX(i, 0) as i64, -w);

        /* all vars wrt v */
        jm.set(IDX(i, 0) as i64, IDX(i, 1) as i64, u * u);
        jm.set(IDX(i, 1) as i64, IDX(i, 1) as i64, -u * u);

        /* all vars wrt w */
        jm.set(IDX(i, 0) as i64, IDX(i, 2) as i64, -u);
        jm.set(IDX(i, 1) as i64, IDX(i, 2) as i64, u);
        jm.set(IDX(i, 2) as i64, IDX(i, 2) as i64, -ONE / ep - u);
    }

    /* return success */
    0
}

/* Set the initial condition */
fn SetIC(y: &mut NVector, udata: &BrusselatorData) -> i32 {
    let N = udata.N as usize;
    let a = udata.a;
    let b = udata.b;
    let dx = udata.dx;

    /* Set initial conditions into y */
    for i in 0..N {
        let x = i as f64 * dx;
        let p = 0.1 * SUNRexp(-(SUNSQR(x - 0.5)) / 0.1);
        y.data[IDX(i, 0)] = a + p;
        y.data[IDX(i, 1)] = b / a + p;
        y.data[IDX(i, 2)] = b + p;
    }

    /* return success */
    0
}

/* C " %.16e" (leading space printed separately) */
fn fmt_e16(x: f64) -> String {
    fmt_e(x, 0, 16)
}

/* Main Program */
fn main() {
    /* general problem parameters */
    let T0: f64 = 0.0; /* initial time                    */
    let Tf: f64 = 10.0; /* final time                      */
    let Nt: i32 = 100; /* total number of output times    */
    let Nvar: i64 = 3; /* number of solution fields       */
    let N: i64 = 200; /* spatial mesh size (N intervals) */
    let a: f64 = 1.0; /* problem parameters              */
    let b: f64 = 3.5;
    let c: f64 = 0.25;
    let ep: f64 = 1.0e-6; /* stiffness parameter */
    let reltol: f64 = 1.0e-6; /* tolerances          */
    let abstol: f64 = 1.0e-10;

    /* Create the SUNDIALS context object for this simulation */
    let ctx = SUNContext_Create();

    /* store the inputs in the UserData structure */
    let udata = BrusselatorData {
        N,
        a,
        b,
        c,
        ep,
        dx: 1.0 / N as f64, /* periodic BC, divide by N not N-1 */
    };

    /* set total allocated vector length */
    let NEQ: i64 = Nvar * udata.N;

    /* set the slow step size */
    let hs: f64 = 0.5 * (udata.dx / SUNRabs(c));

    /* Initial problem output */
    println!("\n1D Advection-Reaction example problem:");
    println!("    N = {},  NEQ = {}", udata.N, NEQ);
    println!(
        "    problem parameters:  a = {},  b = {},  ep = {}",
        fmt_g(udata.a, 0, 6),
        fmt_g(udata.b, 0, 6),
        fmt_g(udata.ep, 0, 6)
    );
    println!("    advection coefficient:  c = {}", fmt_g(udata.c, 0, 6));
    println!(
        "    reltol = {},  abstol = {}\n",
        fmt_e(reltol, 0, 1),
        fmt_e(abstol, 0, 1)
    );

    /* Create solution vector */
    let mut y = N_VNew_Serial(NEQ, &ctx); /* Create vector for solution */

    /* Set initial condition */
    let retval = SetIC(&mut y, &udata);
    assert!(retval >= 0, "SetIC failed");

    /* Create vector masks */
    let mut umask = N_VClone(&y);
    let mut vmask = N_VClone(&y);
    let mut wmask = N_VClone(&y);

    /* Set mask array values for each solution component */
    N_VConst(0.0, &mut umask);
    for i in 0..N as usize {
        umask.data[IDX(i, 0)] = 1.0;
    }

    N_VConst(0.0, &mut vmask);
    for i in 0..N as usize {
        vmask.data[IDX(i, 1)] = 1.0;
    }

    N_VConst(0.0, &mut wmask);
    for i in 0..N as usize {
        wmask.data[IDX(i, 2)] = 1.0;
    }

    /*
     * Create the fast integrator and set options
     */

    /* Initialize matrix and linear solver data structures */
    let A = SUNBandMatrix(NEQ, 4, 4, &ctx);
    let LS = SUNLinSol_Band(&y, &A, &ctx);

    /* Initialize the fast integrator. */
    let mut inner_arkode_mem = ARKStepCreate(None, Some(ff), T0, &y, &ctx).expect("ARKStepCreate");

    /* Attach user data to fast integrator */
    let retval = ARKodeSetUserData(
        &mut inner_arkode_mem,
        Some(Box::new(udata)),
    );
    assert!(retval >= 0, "ARKodeSetUserData failed with retval = {}", retval);

    /* Set the fast method */
    let retval = ARKStepSetTableNum(&mut inner_arkode_mem, ARKODE_ARK324L2SA_DIRK_4_2_3, -1);
    assert!(retval >= 0, "ARKStepSetTableNum failed with retval = {}", retval);

    /* Specify fast tolerances */
    let retval = ARKodeSStolerances(&mut inner_arkode_mem, reltol, abstol);
    assert!(retval >= 0, "ARKodeSStolerances failed with retval = {}", retval);

    /* Attach matrix and linear solver */
    let retval = ARKodeSetLinearSolver(&mut inner_arkode_mem, LS, Some(A));
    assert!(retval >= 0, "ARKodeSetLinearSolver failed with retval = {}", retval);

    /* Set the Jacobian routine */
    let retval = ARKodeSetJacFn(&mut inner_arkode_mem, Some(Jf));
    assert!(retval >= 0, "ARKodeSetJacFn failed with retval = {}", retval);

    /* Create inner stepper */
    let inner_stepper =
        ARKodeCreateMRIStepInnerStepper(inner_arkode_mem).expect("ARKodeCreateMRIStepInnerStepper");

    /*
     * Create the slow integrator and set options
     */

    /* Initialize the slow integrator. */
    let mut arkode_mem =
        MRIStepCreate(Some(fs), None, T0, &y, inner_stepper, &ctx).expect("MRIStepCreate");

    /* Pass udata to user functions */
    let retval = ARKodeSetUserData(&mut arkode_mem, Some(Box::new(udata)));
    assert!(retval >= 0, "ARKodeSetUserData failed with retval = {}", retval);

    /* Set the slow step size */
    let retval = ARKodeSetFixedStep(&mut arkode_mem, hs);
    assert!(retval >= 0, "ARKodeSetFixedStep failed with retval = {}", retval);

    /* output spatial mesh to disk (add extra point for periodic BC) */
    let mut FID = std::fs::File::create("mesh.txt").unwrap();
    for i in 0..=(N as usize) {
        let _ = writeln!(FID, "  {}", fmt_e16(udata.dx * i as f64));
    }
    drop(FID);

    /* Open output stream for results, access data arrays */
    let mut UFID = std::fs::File::create("u.txt").unwrap();
    let mut VFID = std::fs::File::create("v.txt").unwrap();
    let mut WFID = std::fs::File::create("w.txt").unwrap();

    /* output initial condition to disk (extra output for periodic BC) */
    for i in 0..N as usize {
        let _ = write!(UFID, " {}", fmt_e16(y.data[IDX(i, 0)]));
    }
    let _ = write!(UFID, " {}", fmt_e16(y.data[IDX(0, 0)]));
    let _ = writeln!(UFID);

    for i in 0..N as usize {
        let _ = write!(VFID, " {}", fmt_e16(y.data[IDX(i, 1)]));
    }
    let _ = write!(VFID, " {}", fmt_e16(y.data[IDX(0, 1)]));
    let _ = writeln!(VFID);

    for i in 0..N as usize {
        let _ = write!(WFID, " {}", fmt_e16(y.data[IDX(i, 2)]));
    }
    let _ = write!(WFID, " {}", fmt_e16(y.data[IDX(0, 2)]));
    let _ = writeln!(WFID);

    /* Main time-stepping loop */
    let mut t = T0;
    let dTout = (Tf - T0) / Nt as f64;
    let mut tout = T0 + dTout;
    println!("        t      ||u||_rms   ||v||_rms   ||w||_rms");
    println!("   ----------------------------------------------");
    for _iout in 0..Nt {
        /* call integrator */
        let retval = ARKodeEvolve(&mut arkode_mem, tout, &mut y, &mut t, ARK_NORMAL);
        if retval < 0 {
            eprintln!("\nSUNDIALS_ERROR: ARKodeEvolve() failed with retval = {}\n", retval);
            break;
        }

        /* access/print solution statistics */
        let mut u = N_VWL2Norm(&y, &umask);
        u = SUNRsqrt(u * u / N as f64);
        let mut v = N_VWL2Norm(&y, &vmask);
        v = SUNRsqrt(v * v / N as f64);
        let mut w = N_VWL2Norm(&y, &wmask);
        w = SUNRsqrt(w * w / N as f64);
        println!(
            "  {}  {}  {}  {}",
            fmt_f(t, 10, 6),
            fmt_f(u, 10, 6),
            fmt_f(v, 10, 6),
            fmt_f(w, 10, 6)
        );

        /* output results to disk (extr output for periodic BC) */
        for i in 0..N as usize {
            let _ = write!(UFID, " {}", fmt_e16(y.data[IDX(i, 0)]));
        }
        let _ = write!(UFID, " {}", fmt_e16(y.data[IDX(0, 0)]));
        let _ = writeln!(UFID);

        for i in 0..N as usize {
            let _ = write!(VFID, " {}", fmt_e16(y.data[IDX(i, 1)]));
        }
        let _ = write!(VFID, " {}", fmt_e16(y.data[IDX(0, 1)]));
        let _ = writeln!(VFID);

        for i in 0..N as usize {
            let _ = write!(WFID, " {}", fmt_e16(y.data[IDX(i, 2)]));
        }
        let _ = write!(WFID, " {}", fmt_e16(y.data[IDX(0, 2)]));
        let _ = writeln!(WFID);

        /* successful solve: update output time */
        tout += dTout;
        tout = if tout > Tf { Tf } else { tout };
    }
    println!("   ----------------------------------------------");
    drop(UFID);
    drop(VFID);
    drop(WFID);

    /* Get some slow integrator statistics */
    let mut nsts: i64 = 0;
    let mut nfse: i64 = 0;
    ARKodeGetNumSteps(&mut arkode_mem, &mut nsts);
    ARKodeGetNumRhsEvals(&mut arkode_mem, 0, &mut nfse);

    /* Get some fast integrator statistics (borrow the inner integrator
    back out of the outer step memory) */
    let mut nstf: i64 = 0;
    let mut nstf_a: i64 = 0;
    let mut nffi: i64 = 0;
    let mut nsetups: i64 = 0;
    let mut netf: i64 = 0;
    let mut nni: i64 = 0;
    let mut ncfn: i64 = 0;
    let mut nje: i64 = 0;
    let mut nfeLS: i64 = 0;
    {
        let step_mem = arkode_mem
            .step_mem
            .as_mut()
            .unwrap()
            .downcast_mut::<ARKodeMRIStepMem>()
            .unwrap();
        let inner_arkode_mem = step_mem
            .stepper
            .content
            .as_mut()
            .unwrap()
            .downcast_mut::<ARKodeMem>()
            .unwrap();
        ARKodeGetNumSteps(inner_arkode_mem, &mut nstf);
        ARKodeGetNumStepAttempts(inner_arkode_mem, &mut nstf_a);
        ARKodeGetNumRhsEvals(inner_arkode_mem, 1, &mut nffi);
        ARKodeGetNumLinSolvSetups(inner_arkode_mem, &mut nsetups);
        ARKodeGetNumErrTestFails(inner_arkode_mem, &mut netf);
        ARKodeGetNumNonlinSolvIters(inner_arkode_mem, &mut nni);
        ARKodeGetNumNonlinSolvConvFails(inner_arkode_mem, &mut ncfn);
        ARKodeGetNumJacEvals(inner_arkode_mem, &mut nje);
        ARKodeGetNumLinRhsEvals(inner_arkode_mem, &mut nfeLS);
    }

    /* Print some final statistics */
    println!("\nFinal Solver Statistics:");
    println!("   Slow Steps: nsts = {}", nsts);
    println!("   Fast Steps: nstf = {} (attempted = {})", nstf, nstf_a);
    println!("   Total RHS evals:  Fs = {},  Ff = {}", nfse, nffi);
    println!("   Total number of fast error test failures = {}", netf);
    println!("   Total linear solver setups = {}", nsetups);
    println!("   Total RHS evals for setting up the linear system = {}", nfeLS);
    println!("   Total number of Jacobian evaluations = {}", nje);
    println!("   Total number of Newton iterations = {}", nni);
    println!(
        "   Total number of nonlinear solver convergence failures = {}",
        ncfn
    );

    /* Clean up and return with successful completion */
    drop(y);
    drop(umask);
    drop(vmask);
    drop(wmask);
    let mut slot = Some(arkode_mem);
    ARKodeFree(&mut slot);
}
