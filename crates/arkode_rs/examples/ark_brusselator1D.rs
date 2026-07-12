/*---------------------------------------------------------------
 * Translation of examples/arkode/C_serial/ark_brusselator1D.c
 * (SUNDIALS 7.7.0).
 *
 * Example problem:
 *
 * The following test simulates a brusselator problem from chemical
 * kinetics.  This is n PDE system with 3 components, Y = [u,v,w],
 * satisfying the equations,
 *    u_t = du*u_xx + a - (w+1)*u + v*u^2
 *    v_t = dv*v_xx + w*u - v*u^2
 *    w_t = dw*w_xx + (b-w)/ep - w*u
 * for t in [0, 80], x in [0, 1], with initial conditions
 *    u(0,x) =  a  + 0.1*sin(pi*x)
 *    v(0,x) = b/a + 0.1*sin(pi*x)
 *    w(0,x) =  b  + 0.1*sin(pi*x),
 * and with stationary boundary conditions, i.e.
 *    u_t(t,0) = u_t(t,1) = 0,
 *    v_t(t,0) = v_t(t,1) = 0,
 *    w_t(t,0) = w_t(t,1) = 0.
 * Note: these can also be implemented as Dirichlet boundary
 * conditions with values identical to the initial conditions.
 *
 * The spatial derivatives are computed using second-order
 * centered differences, with the data distributed over N points
 * on a uniform spatial grid.
 *
 * This program solves the problem with the DIRK method, using a
 * Newton iteration with the SUNBAND band linear solver, and a
 * user-supplied Jacobian routine.
 *
 * 100 outputs are printed at equal intervals, and run statistics
 * are printed at the end.
 *---------------------------------------------------------------*/
#![allow(non_snake_case)]

use std::io::Write;

use arkode_rs::arkode::{ARKodeEvolve, ARKodeFree, ARKodeSStolerances};
use arkode_rs::arkode_arkstep::ARKStepCreate;
use arkode_rs::arkode_io::{
    ARKodeGetNumErrTestFails, ARKodeGetNumLinSolvSetups, ARKodeGetNumNonlinSolvConvFails,
    ARKodeGetNumNonlinSolvIters, ARKodeGetNumRhsEvals, ARKodeGetNumStepAttempts,
    ARKodeGetNumSteps, ARKodeSetAutonomous, ARKodeSetUserData,
};
use arkode_rs::arkode_ls::{
    ARKodeGetNumJacEvals, ARKodeGetNumLinRhsEvals, ARKodeSetJacFn, ARKodeSetLinearSolver,
};
use arkode_rs::sundials_utils::{fmt_e, fmt_f, fmt_g};
use arkode_rs::*;

/* accessor between (x,v) location and 1D NVector array */
fn IDX(x: usize, v: usize) -> usize {
    3 * x + v
}

/* user data structure */
struct UData {
    N: i64,  /* number of intervals     */
    dx: f64, /* mesh spacing            */
    a: f64,  /* constant forcing on u   */
    b: f64,  /* steady-state value of w */
    du: f64, /* diffusion coeff for u   */
    dv: f64, /* diffusion coeff for v   */
    dw: f64, /* diffusion coeff for w   */
    ep: f64, /* stiffness parameter     */
}

/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
fn f(_t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let udata = user_data.as_ref().unwrap().downcast_ref::<UData>().unwrap();
    let N = udata.N; /* set variable shortcuts */
    let a = udata.a;
    let b = udata.b;
    let ep = udata.ep;
    let du = udata.du;
    let dv = udata.dv;
    let dw = udata.dw;
    let dx = udata.dx;

    N_VConst(0.0, ydot); /* initialize ydot to zero */

    /* iterate over domain, computing all equations */
    let uconst = du / dx / dx;
    let vconst = dv / dx / dx;
    let wconst = dw / dx / dx;
    for i in 1..(N - 1) as usize {
        /* set shortcuts */
        let u = y.data[IDX(i, 0)];
        let ul = y.data[IDX(i - 1, 0)];
        let ur = y.data[IDX(i + 1, 0)];
        let v = y.data[IDX(i, 1)];
        let vl = y.data[IDX(i - 1, 1)];
        let vr = y.data[IDX(i + 1, 1)];
        let w = y.data[IDX(i, 2)];
        let wl = y.data[IDX(i - 1, 2)];
        let wr = y.data[IDX(i + 1, 2)];

        /* Fill in ODE RHS for u */
        ydot.data[IDX(i, 0)] = (ul - 2.0 * u + ur) * uconst + a - (w + 1.0) * u + v * u * u;

        /* Fill in ODE RHS for v */
        ydot.data[IDX(i, 1)] = (vl - 2.0 * v + vr) * vconst + w * u - v * u * u;

        /* Fill in ODE RHS for w */
        ydot.data[IDX(i, 2)] = (wl - 2.0 * w + wr) * wconst + (b - w) / ep - w * u;
    }

    /* enforce stationary boundaries */
    ydot.data[IDX(0, 0)] = 0.0;
    ydot.data[IDX(0, 1)] = 0.0;
    ydot.data[IDX(0, 2)] = 0.0;
    ydot.data[IDX((N - 1) as usize, 0)] = 0.0;
    ydot.data[IDX((N - 1) as usize, 1)] = 0.0;
    ydot.data[IDX((N - 1) as usize, 2)] = 0.0;

    0 /* Return with success */
}

/* Jacobian routine to compute J(t,y) = df/dy. */
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
    let udata = user_data.as_ref().unwrap().downcast_ref::<UData>().unwrap();
    SUNMatZero(j); /* Initialize Jacobian to zero */

    /* Fill in the Laplace matrix */
    laplace_matrix(1.0, j, udata);

    /* Add in the Jacobian of the reaction terms matrix */
    reaction_jac(1.0, y, j, udata);

    0 /* Return with success */
}

/*-------------------------------
 * Private helper functions
 *-------------------------------*/

/* Routine to compute the stiffness matrix from (L*y), scaled by the factor c.
   We add the result into Jac and do not erase what was already there */
fn laplace_matrix(c: f64, jac: &mut SUNMatrix, udata: &UData) -> i32 {
    let N = udata.N; /* set shortcuts */
    let dx = udata.dx;

    let jm = match jac {
        SUNMatrix::Band(m) => m,
        _ => return 1,
    };

    /* iterate over intervals, filling in Jacobian of (L*y) (the C
       SM_ELEMENT_B(J, r, c) += v accumulations) */
    for i in 1..(N - 1) as usize {
        let r0 = IDX(i, 0) as i64;
        let r1 = IDX(i, 1) as i64;
        let r2 = IDX(i, 2) as i64;
        let cl0 = IDX(i - 1, 0) as i64;
        let cl1 = IDX(i - 1, 1) as i64;
        let cl2 = IDX(i - 1, 2) as i64;
        let cr0 = IDX(i + 1, 0) as i64;
        let cr1 = IDX(i + 1, 1) as i64;
        let cr2 = IDX(i + 1, 2) as i64;

        jm.set(r0, cl0, jm.get(r0, cl0) + c * udata.du / dx / dx);
        jm.set(r1, cl1, jm.get(r1, cl1) + c * udata.dv / dx / dx);
        jm.set(r2, cl2, jm.get(r2, cl2) + c * udata.dw / dx / dx);
        jm.set(r0, r0, jm.get(r0, r0) + -c * 2.0 * udata.du / dx / dx);
        jm.set(r1, r1, jm.get(r1, r1) + -c * 2.0 * udata.dv / dx / dx);
        jm.set(r2, r2, jm.get(r2, r2) + -c * 2.0 * udata.dw / dx / dx);
        jm.set(r0, cr0, jm.get(r0, cr0) + c * udata.du / dx / dx);
        jm.set(r1, cr1, jm.get(r1, cr1) + c * udata.dv / dx / dx);
        jm.set(r2, cr2, jm.get(r2, cr2) + c * udata.dw / dx / dx);
    }

    0 /* Return with success */
}

/* Routine to compute the Jacobian matrix from R(y), scaled by the factor c.
   We add the result into Jac and do not erase what was already there */
fn reaction_jac(c: f64, y: &NVector, jac: &mut SUNMatrix, udata: &UData) -> i32 {
    let N = udata.N; /* set shortcuts */
    let ep = udata.ep;

    let jm = match jac {
        SUNMatrix::Band(m) => m,
        _ => return 1,
    };

    /* iterate over nodes, filling in Jacobian of reaction terms */
    for i in 1..(N - 1) as usize {
        let u = y.data[IDX(i, 0)]; /* set nodal value shortcuts */
        let v = y.data[IDX(i, 1)];
        let w = y.data[IDX(i, 2)];

        let r0 = IDX(i, 0) as i64;
        let r1 = IDX(i, 1) as i64;
        let r2 = IDX(i, 2) as i64;

        /* all vars wrt u */
        jm.set(r0, r0, jm.get(r0, r0) + c * (2.0 * u * v - (w + 1.0)));
        jm.set(r1, r0, jm.get(r1, r0) + c * (w - 2.0 * u * v));
        jm.set(r2, r0, jm.get(r2, r0) + c * (-w));

        /* all vars wrt v */
        jm.set(r0, r1, jm.get(r0, r1) + c * (u * u));
        jm.set(r1, r1, jm.get(r1, r1) + c * (-u * u));

        /* all vars wrt w */
        jm.set(r0, r2, jm.get(r0, r2) + c * (-u));
        jm.set(r1, r2, jm.get(r1, r2) + c * u);
        jm.set(r2, r2, jm.get(r2, r2) + c * (-1.0 / ep - u));
    }

    0 /* Return with success */
}

/* Check if a SUNDIALS function returned a negative flag */
fn check_flag(flag: i32, funcname: &str) -> bool {
    if flag < 0 {
        eprintln!("\nSUNDIALS_ERROR: {}() failed with flag = {}\n", funcname, flag);
        return true;
    }
    false
}

/* Main Program */
fn main() {
    /* general problem parameters */
    let T0: f64 = 0.0; /* initial time */
    let Tf: f64 = 10.0; /* final time */
    let Nt: i32 = 100; /* total number of output times */
    let Nvar: i64 = 3; /* number of solution fields */
    let N: i64 = 201; /* spatial mesh size */
    let a: f64 = 0.6; /* problem parameters */
    let b: f64 = 2.0;
    let du: f64 = 0.025;
    let dv: f64 = 0.025;
    let dw: f64 = 0.025;
    let ep: f64 = 1.0e-5; /* stiffness parameter */
    let reltol: f64 = 1.0e-6; /* tolerances */
    let abstol: f64 = 1.0e-10;

    /* Create the SUNDIALS context object for this simulation */
    let ctx = SUNContext_Create();

    /* allocate udata structure and store the inputs */
    let mut udata = UData {
        N,
        dx: 0.0,
        a,
        b,
        du,
        dv,
        dw,
        ep,
    };

    /* set total allocated vector length */
    let NEQ = Nvar * udata.N;

    /* Initial problem output */
    println!("\n1D Brusselator PDE test problem:");
    println!("    N = {},  NEQ = {}", udata.N, NEQ);
    println!(
        "    problem parameters:  a = {},  b = {},  ep = {}",
        fmt_g(udata.a, 0, 6),
        fmt_g(udata.b, 0, 6),
        fmt_g(udata.ep, 0, 6)
    );
    println!(
        "    diffusion coefficients:  du = {},  dv = {},  dw = {}",
        fmt_g(udata.du, 0, 6),
        fmt_g(udata.dv, 0, 6),
        fmt_g(udata.dw, 0, 6)
    );
    println!("    reltol = {},  abstol = {}\n", fmt_e(reltol, 0, 1), fmt_e(abstol, 0, 1));

    /* Initialize data structures */
    let mut y = N_VNew_Serial(NEQ, &ctx); /* Create serial vector for solution */
    let mut umask = N_VClone(&y);
    let mut vmask = N_VClone(&y);
    let mut wmask = N_VClone(&y);

    /* Set initial conditions into y */
    udata.dx = 1.0 / (N - 1) as f64; /* set spatial mesh spacing */
    let u_dx = udata.dx;
    let u_n = udata.N;

    let pi = 4.0 * (1.0f64).atan();
    for i in 0..N as usize {
        y.data[IDX(i, 0)] = a + 0.1 * (pi * i as f64 * udata.dx).sin(); /* u */
        y.data[IDX(i, 1)] = b / a + 0.1 * (pi * i as f64 * udata.dx).sin(); /* v */
        y.data[IDX(i, 2)] = b + 0.1 * (pi * i as f64 * udata.dx).sin(); /* w */
    }

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

    /* Call ARKStepCreate to initialize the ARK timestepper module and
       specify the right-hand side function in y'=f(t,y), the initial time
       T0, and the initial dependent variable vector y.  Note: since this
       problem is fully implicit, we set f_E to NULL and f_I to f. */
    let mut arkode_mem = ARKStepCreate(None, Some(f), T0, &y, &ctx).expect("ARKStepCreate");

    /* Set routines */
    let flag = ARKodeSetUserData(&mut arkode_mem, Some(Box::new(udata)));
    if check_flag(flag, "ARKodeSetUserData") {
        return;
    }
    let flag = ARKodeSStolerances(&mut arkode_mem, reltol, abstol); /* Specify tolerances */
    if check_flag(flag, "ARKodeSStolerances") {
        return;
    }

    /* Initialize band matrix data structure and solver -- A will be factored,
       so set smu to ml+mu */
    let a_mat = SUNBandMatrix(NEQ, 4, 4, &ctx);
    let ls = SUNLinSol_Band(&y, &a_mat, &ctx);

    /* Linear solver interface */
    let flag = ARKodeSetLinearSolver(&mut arkode_mem, ls, Some(a_mat)); /* Attach matrix and linear solver */
    if check_flag(flag, "ARKodeSetLinearSolver") {
        return;
    }

    let flag = ARKodeSetJacFn(&mut arkode_mem, Some(jac)); /* Set the Jacobian routine */
    if check_flag(flag, "ARKodeSetJacFn") {
        return;
    }

    let flag = ARKodeSetAutonomous(&mut arkode_mem, true);
    if check_flag(flag, "ARKodeSetAutonomous") {
        return;
    }

    /* output spatial mesh to disk */
    let mut fid = std::fs::File::create("bruss_mesh.txt").expect("fopen");
    for i in 0..N {
        let _ = writeln!(fid, "  {}", fmt_e(u_dx * i as f64, 0, 16));
    }
    drop(fid);

    /* Open output streams for results */
    let mut ufid = std::fs::File::create("bruss_u.txt").expect("fopen");
    let mut vfid = std::fs::File::create("bruss_v.txt").expect("fopen");
    let mut wfid = std::fs::File::create("bruss_w.txt").expect("fopen");

    /* output initial condition to disk */
    for i in 0..u_n as usize {
        let _ = write!(ufid, " {}", fmt_e(y.data[IDX(i, 0)], 0, 16));
    }
    for i in 0..u_n as usize {
        let _ = write!(vfid, " {}", fmt_e(y.data[IDX(i, 1)], 0, 16));
    }
    for i in 0..u_n as usize {
        let _ = write!(wfid, " {}", fmt_e(y.data[IDX(i, 2)], 0, 16));
    }
    let _ = writeln!(ufid);
    let _ = writeln!(vfid);
    let _ = writeln!(wfid);

    /* Main time-stepping loop: calls ARKodeEvolve to perform the integration,
       then prints results.  Stops when the final time has been reached */
    let mut t = T0;
    let dTout = (Tf - T0) / Nt as f64;
    let mut tout = T0 + dTout;
    println!("        t      ||u||_rms   ||v||_rms   ||w||_rms");
    println!("   ----------------------------------------------");
    for _iout in 0..Nt {
        let flag = ARKodeEvolve(&mut arkode_mem, tout, &mut y, &mut t, ARK_NORMAL); /* call integrator */
        if check_flag(flag, "ARKodeEvolve") {
            break;
        }
        let mut u = N_VWL2Norm(&y, &umask); /* access/print solution statistics */
        u = (u * u / N as f64).sqrt();
        let mut v = N_VWL2Norm(&y, &vmask);
        v = (v * v / N as f64).sqrt();
        let mut w = N_VWL2Norm(&y, &wmask);
        w = (w * w / N as f64).sqrt();
        println!(
            "  {}  {}  {}  {}",
            fmt_f(t, 10, 6),
            fmt_f(u, 10, 6),
            fmt_f(v, 10, 6),
            fmt_f(w, 10, 6)
        );
        if flag >= 0 {
            /* successful solve: update output time */
            tout += dTout;
            tout = if tout > Tf { Tf } else { tout };
        } else {
            /* unsuccessful solve: break */
            eprintln!("Solver failure, stopping integration");
            break;
        }

        /* output results to disk */
        for i in 0..u_n as usize {
            let _ = write!(ufid, " {}", fmt_e(y.data[IDX(i, 0)], 0, 16));
        }
        for i in 0..u_n as usize {
            let _ = write!(vfid, " {}", fmt_e(y.data[IDX(i, 1)], 0, 16));
        }
        for i in 0..u_n as usize {
            let _ = write!(wfid, " {}", fmt_e(y.data[IDX(i, 2)], 0, 16));
        }
        let _ = writeln!(ufid);
        let _ = writeln!(vfid);
        let _ = writeln!(wfid);
    }
    println!("   ----------------------------------------------");
    drop(ufid);
    drop(vfid);
    drop(wfid);

    /* Print some final statistics */
    let mut nst: i64 = 0;
    let mut nst_a: i64 = 0;
    let mut nfe: i64 = 0;
    let mut nfi: i64 = 0;
    let mut nsetups: i64 = 0;
    let mut nje: i64 = 0;
    let mut nfeLS: i64 = 0;
    let mut nni: i64 = 0;
    let mut ncfn: i64 = 0;
    let mut netf: i64 = 0;

    let flag = ARKodeGetNumSteps(&mut arkode_mem, &mut nst);
    check_flag(flag, "ARKodeGetNumSteps");
    let flag = ARKodeGetNumStepAttempts(&mut arkode_mem, &mut nst_a);
    check_flag(flag, "ARKodeGetNumStepAttempts");
    let flag = ARKodeGetNumRhsEvals(&mut arkode_mem, 0, &mut nfe);
    check_flag(flag, "ARKodeGetNumRhsEvals");
    let flag = ARKodeGetNumRhsEvals(&mut arkode_mem, 1, &mut nfi);
    check_flag(flag, "ARKodeGetNumRhsEvals");
    let flag = ARKodeGetNumLinSolvSetups(&mut arkode_mem, &mut nsetups);
    check_flag(flag, "ARKodeGetNumLinSolvSetups");
    let flag = ARKodeGetNumErrTestFails(&mut arkode_mem, &mut netf);
    check_flag(flag, "ARKodeGetNumErrTestFails");
    let flag = ARKodeGetNumNonlinSolvIters(&mut arkode_mem, &mut nni);
    check_flag(flag, "ARKodeGetNumNonlinSolvIters");
    let flag = ARKodeGetNumNonlinSolvConvFails(&mut arkode_mem, &mut ncfn);
    check_flag(flag, "ARKodeGetNumNonlinSolvConvFails");
    let flag = ARKodeGetNumJacEvals(&mut arkode_mem, &mut nje);
    check_flag(flag, "ARKodeGetNumJacEvals");
    let flag = ARKodeGetNumLinRhsEvals(&mut arkode_mem, &mut nfeLS);
    check_flag(flag, "ARKodeGetNumLinRhsEvals");

    println!("\nFinal Solver Statistics:");
    println!("   Internal solver steps = {} (attempted = {})", nst, nst_a);
    println!("   Total RHS evals:  Fe = {},  Fi = {}", nfe, nfi);
    println!("   Total linear solver setups = {}", nsetups);
    println!("   Total RHS evals for setting up the linear system = {}", nfeLS);
    println!("   Total number of Jacobian evaluations = {}", nje);
    println!("   Total number of Newton iterations = {}", nni);
    println!("   Total number of nonlinear solver convergence failures = {}", ncfn);
    println!("   Total number of error test failures = {}\n", netf);

    /* Clean up and return with successful completion */
    drop(y); /* Free vectors */
    drop(umask);
    drop(vmask);
    drop(wmask);
    let mut slot = Some(arkode_mem);
    ARKodeFree(&mut slot); /* Free integrator memory */
}
