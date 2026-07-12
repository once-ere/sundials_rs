/*---------------------------------------------------------------
 * Translation of examples/arkode/C_serial/ark_heat1D_adapt.c
 * (SUNDIALS 7.7.0).
 *
 * Example problem:
 *
 * The following test simulates a simple 1D heat equation,
 *    u_t = k*u_xx + f
 * for t in [0, 10], x in [0, 1], with initial conditions
 *    u(0,x) =  0
 * Dirichlet boundary conditions, i.e.
 *    u_t(t,0) = u_t(t,1) = 0,
 * and a heating term of the form
 *    f = 2*exp(-200*(x-0.25)*(x-0.25))
 *        - exp(-400*(x-0.7)*(x-0.7))
 *        + exp(-500*(x-0.4)*(x-0.4))
 *        - 2*exp(-600*(x-0.55)*(x-0.55));
 *
 * The spatial derivatives are computed using a three-point
 * centered stencil (second order for a uniform mesh).  The data
 * is initially uniformly distributed over N points in the interval
 * [0, 1], but as the simulation proceeds the mesh is adapted.
 *
 * This program solves the problem with a DIRK method, solved with
 * a Newton iteration and SUNLinSol_PCG linear solver, with a
 * user-supplied Jacobian-vector product routine.
 *---------------------------------------------------------------*/
#![allow(non_snake_case)]

use std::io::Write;

use arkode_rs::arkode::{ARKodeEvolve, ARKodeFree, ARKodeResize, ARKodeSStolerances};
use arkode_rs::arkode_arkstep::ARKStepCreate;
use arkode_rs::arkode_arkstep_io::ARKStepSetAdaptivityMethod;
use arkode_rs::arkode_io::{
    ARKodeGetCurrentStep, ARKodeGetLastStep, ARKodeGetNumNonlinSolvIters,
    ARKodeSetLinear, ARKodeSetMaxNumSteps, ARKodeSetPredictorMethod, ARKodeSetStopTime, ARKodeSetUserData,
};
use arkode_rs::arkode_ls::{
    ARKodeGetNumLinIters, ARKodeSetJacTimes, ARKodeSetLinearSolver,
};
use arkode_rs::sundials_utils::{fmt_e, fmt_g};
use arkode_rs::*;

/* constants */
const ZERO: f64 = 0.0;
const PT25: f64 = 0.25;
const PT4: f64 = 0.4;
const PT5: f64 = 0.5;
const PT55: f64 = 0.55;
const PT7: f64 = 0.7;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;
const TWOHUNDRED: f64 = 200.0;
const FOURHUNDRED: f64 = 400.0;
const FIVEHUNDRED: f64 = 500.0;
const SIXHUNDRED: f64 = 600.0;

/* user data structure */
struct HeatData {
    N: i64,          /* current number of intervals */
    x: Vec<f64>,     /* current mesh */
    k: f64,          /* diffusion coefficient */
    refine_tol: f64, /* adaptivity tolerance */
}

/*--------------------------------
 * Functions called by the solver
 *--------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
fn f(_t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let udata = user_data.as_ref().unwrap().downcast_ref::<HeatData>().unwrap();
    let N = udata.N; /* set variable shortcuts */
    let k = udata.k;
    let x = &udata.x;

    /* Initialize ydot to zero - also handles boundary conditions */
    N_VConst(ZERO, ydot);

    /* iterate over domain interior, computing all equations */
    for i in 1..(N - 1) as usize {
        /* interior */
        let dxL = x[i] - x[i - 1];
        let dxR = x[i + 1] - x[i];
        ydot.data[i] = y.data[i - 1] * k * TWO / (dxL * (dxL + dxR))
            - y.data[i] * k * TWO / (dxL * dxR)
            + y.data[i + 1] * k * TWO / (dxR * (dxL + dxR))
            + TWO * (-TWOHUNDRED * (x[i] - PT25) * (x[i] - PT25)).exp() /* source term */
            - (-FOURHUNDRED * (x[i] - PT7) * (x[i] - PT7)).exp()
            + (-FIVEHUNDRED * (x[i] - PT4) * (x[i] - PT4)).exp()
            - TWO * (-SIXHUNDRED * (x[i] - PT55) * (x[i] - PT55)).exp();
    }

    0 /* Return with success */
}

/* Jacobian routine to compute J(t,y) = df/dy. */
fn jac(
    v: &NVector,
    jv: &mut NVector,
    _t: f64,
    _y: &NVector,
    _fy: &NVector,
    user_data: &mut UserData,
    _tmp: &mut NVector,
) -> i32 {
    let udata = user_data.as_ref().unwrap().downcast_ref::<HeatData>().unwrap();
    let N = udata.N; /* variable shortcuts */
    let k = udata.k;
    let x = &udata.x;

    /* initialize Jv product to zero - also handles boundary conditions */
    N_VConst(ZERO, jv);

    /* iterate over domain, computing all Jacobian-vector products */
    for i in 1..(N - 1) as usize {
        let dxL = x[i] - x[i - 1];
        let dxR = x[i + 1] - x[i];
        jv.data[i] = v.data[i - 1] * k * TWO / (dxL * (dxL + dxR))
            - v.data[i] * k * TWO / (dxL * dxR)
            + v.data[i + 1] * k * TWO / (dxR * (dxL + dxR));
    }

    0 /* Return with success */
}

/*-------------------------------
 * Private helper functions
 *-------------------------------*/

/* Adapts the current mesh, using a simple adaptivity strategy of
   refining when an approximation of the scaled second-derivative is
   too large.  We only do this in one sweep, so no attempt is made to
   ensure the resulting mesh meets these same criteria after adaptivity:
      y [input] -- the current solution vector
      Nnew [output] -- the size of the new mesh
      udata [input] -- the current system information
   The return for this function is the new mesh (None on failure). */
fn adapt_mesh(y: &NVector, Nnew: &mut i64, udata: &HeatData) -> Option<Vec<f64>> {
    /* Access current solution and mesh arrays */
    let xold = &udata.x;
    let n = udata.N as usize;

    /* create marking array */
    let mut marks = vec![0i32; n - 1];

    /* perform marking:
        0 -> leave alone
        1 -> refine */
    for i in 1..n - 1 {
        /* approximate scaled second-derivative */
        let ydd = y.data[i - 1] - TWO * y.data[i] + y.data[i + 1];

        /* check for refinement */
        if SUNRabs(ydd) > udata.refine_tol {
            marks[i - 1] = 1;
            marks[i] = 1;
        }
    }

    /* allocate new mesh */
    let mut num_refine: i64 = 0;
    for &m in marks.iter() {
        if m == 1 {
            num_refine += 1;
        }
    }
    let N_new = udata.N + num_refine;
    *Nnew = N_new; /* Store new array length */
    let mut xnew = vec![0.0f64; N_new as usize];

    /* fill new mesh */
    xnew[0] = xold[0]; /* store endpoints */
    xnew[(N_new - 1) as usize] = xold[n - 1];
    let mut j = 1usize;
    /* iterate over old intervals */
    for i in 0..n - 1 {
        /* if mark is 0, reuse old interval */
        if marks[i] == 0 {
            xnew[j] = xold[i + 1];
            j += 1;
            continue;
        }

        /* if mark is 1, refine old interval */
        if marks[i] == 1 {
            xnew[j] = PT5 * (xold[i] + xold[i + 1]);
            j += 1;
            xnew[j] = xold[i + 1];
            j += 1;
            continue;
        }
    }

    /* verify that new mesh is legal */
    for i in 0..(N_new - 1) as usize {
        if xnew[i + 1] <= xnew[i] {
            eprintln!("adapt_mesh error: illegal mesh created");
            return None;
        }
    }

    Some(xnew) /* Return with success */
}

/* Projects one vector onto another:
      Nold [input] -- the size of the old mesh
      xold [input] -- the old mesh
      yold [input] -- the vector defined over the old mesh
      Nnew [input] -- the size of the new mesh
      xnew [input] -- the new mesh
      ynew [output] -- the vector defined over the new mesh
                       (allocated prior to calling project) */
fn project(
    Nold: i64,
    xold: &[f64],
    yold: &NVector,
    Nnew: i64,
    xnew: &[f64],
    ynew: &mut NVector,
) -> i32 {
    /* loop over new mesh, finding corresponding interval within old mesh,
       and perform piecewise linear interpolation from yold to ynew */
    let mut iv = 0usize;
    for i in 0..Nnew as usize {
        /* find old interval, start with previous value since sorted */
        for j in iv..(Nold - 1) as usize {
            if xnew[i] >= xold[j] && xnew[i] <= xold[j + 1] {
                iv = j;
                break;
            }
            iv = (Nold - 1) as usize; /* just in case it wasn't found above */
        }

        /* perform interpolation */
        ynew.data[i] = yold.data[iv] * (xnew[i] - xold[iv + 1]) / (xold[iv] - xold[iv + 1])
            + yold.data[iv + 1] * (xnew[i] - xold[iv]) / (xold[iv + 1] - xold[iv]);
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
    let Tf: f64 = 1.0; /* final time */
    let rtol: f64 = 1.0e-3; /* relative tolerance */
    let atol: f64 = 1.0e-10; /* absolute tolerance */
    let hscale: f64 = 1.0; /* time step change factor on resizes */
    let N: i64 = 21; /* initial spatial mesh size */
    let refine: f64 = 3.0e-3; /* adaptivity refinement tolerance */
    let k: f64 = 0.5; /* heat conductivity */
    let mut nni_tot: i64 = 0;
    let mut nli_tot: i64 = 0;
    let mut iout: i32 = 0;

    /* Create the SUNDIALS context object for this simulation */
    let ctx = SUNContext_Create();

    /* allocate and fill initial udata structure */
    let mut x = vec![0.0f64; N as usize];
    for (i, xi) in x.iter_mut().enumerate() {
        *xi = ONE * i as f64 / (N - 1) as f64;
    }
    let udata = HeatData {
        N,
        x,
        k,
        refine_tol: refine,
    };

    /* Initial problem output */
    println!("\n1D adaptive Heat PDE test problem:");
    println!("  diffusion coefficient:  k = {}", fmt_g(udata.k, 0, 6));
    println!("  initial N = {}", udata.N);

    /* Initialize data structures */
    let mut y = N_VNew_Serial(N, &ctx); /* Create initial serial vector for solution */
    N_VConst(ZERO, &mut y); /* Set initial conditions */

    /* output mesh to disk */
    let mut xfid = std::fs::File::create("heat_mesh.txt").expect("fopen");

    /* output initial mesh to disk */
    for i in 0..udata.N as usize {
        let _ = write!(xfid, " {}", fmt_e(udata.x[i], 0, 16));
    }
    let _ = writeln!(xfid);

    /* Open output stream for results */
    let mut ufid = std::fs::File::create("heat1D.txt").expect("fopen");

    /* output initial condition to disk */
    for i in 0..udata.N as usize {
        let _ = write!(ufid, " {}", fmt_e(y.data[i], 0, 16));
    }
    let _ = writeln!(ufid);

    /* Initialize the ARK timestepper */
    let mut arkode_mem = ARKStepCreate(None, Some(f), T0, &y, &ctx).expect("ARKStepCreate");

    /* Set routines */
    let flag = ARKodeSetUserData(&mut arkode_mem, Some(Box::new(udata)));
    if check_flag(flag, "ARKodeSetUserData") {
        return;
    }
    let flag = ARKodeSetMaxNumSteps(&mut arkode_mem, 10000); /* Increase max num steps */
    if check_flag(flag, "ARKodeSetMaxNumSteps") {
        return;
    }
    let flag = ARKodeSStolerances(&mut arkode_mem, rtol, atol); /* Specify tolerances */
    if check_flag(flag, "ARKodeSStolerances") {
        return;
    }
    let flag = ARKStepSetAdaptivityMethod(&mut arkode_mem, 2, 1, 0, None); /* Set adaptivity method */
    if check_flag(flag, "ARKodeSetAdaptivityMethod") {
        return;
    }
    let flag = ARKodeSetPredictorMethod(&mut arkode_mem, 0); /* Set predictor method */
    if check_flag(flag, "ARKodeSetPredictorMethod") {
        return;
    }

    /* Specify linearly implicit RHS, with time-dependent Jacobian */
    let flag = ARKodeSetLinear(&mut arkode_mem, 1);
    if check_flag(flag, "ARKodeSetLinear") {
        return;
    }

    /* Initialize PCG solver -- no preconditioning, with up to N iterations */
    let ls = SUNLinSol_PCG(&y, 0, N as i32, &ctx);

    /* Linear solver interface -- set user-supplied J*v routine (no 'jtsetup'
       required) */
    let flag = ARKodeSetLinearSolver(&mut arkode_mem, ls, None); /* Attach linear solver */
    if check_flag(flag, "ARKodeSetLinearSolver") {
        return;
    }
    let flag = ARKodeSetJacTimes(&mut arkode_mem, None, Some(jac)); /* Set the Jacobian routine */
    if check_flag(flag, "ARKodeSetJacTimes") {
        return;
    }

    /* Main time-stepping loop: calls ARKodeEvolve to perform the integration,
       then prints results.  Stops when the final time has been reached */
    let mut t = T0;
    let mut olddt: f64 = ZERO;
    let mut newdt: f64 = ZERO;
    let mut nni: i64 = 0;
    let mut nli: i64 = 0;
    let mut no_resize_data: UserData = None;
    println!(
        "  iout          dt_old                 dt_new               ||u||_rms       N   NNI  NLI"
    );
    println!(
        " ----------------------------------------------------------------------------------------"
    );
    {
        let cur_n = arkode_mem
            .user_data
            .as_ref()
            .unwrap()
            .downcast_ref::<HeatData>()
            .unwrap()
            .N;
        println!(
            " {:4}  {}  {}  {}  {}   {:2}  {:3}",
            iout,
            fmt_e(olddt, 19, 15),
            fmt_e(newdt, 19, 15),
            fmt_e((N_VDotProd(&y, &y) / cur_n as f64).sqrt(), 19, 15),
            cur_n,
            0,
            0
        );
    }
    while t < Tf {
        /* "set" routines */
        let flag = ARKodeSetStopTime(&mut arkode_mem, Tf);
        if check_flag(flag, "ARKodeSetStopTime") {
            return;
        }

        /* call integrator */
        let flag = ARKodeEvolve(&mut arkode_mem, Tf, &mut y, &mut t, ARK_ONE_STEP);
        if check_flag(flag, "ARKodeEvolve") {
            return;
        }

        /* "get" routines */
        let flag = ARKodeGetLastStep(&mut arkode_mem, &mut olddt);
        if check_flag(flag, "ARKodeGetLastStep") {
            return;
        }
        let flag = ARKodeGetCurrentStep(&mut arkode_mem, &mut newdt);
        if check_flag(flag, "ARKodeGetCurrentStep") {
            return;
        }
        let flag = ARKodeGetNumNonlinSolvIters(&mut arkode_mem, &mut nni);
        if check_flag(flag, "ARKodeGetNumNonlinSolvIters") {
            return;
        }
        let flag = ARKodeGetNumLinIters(&mut arkode_mem, &mut nli);
        if check_flag(flag, "ARKodeGetNumLinIters") {
            return;
        }

        /* print current solution stats */
        iout += 1;
        {
            let cur_n = arkode_mem
                .user_data
                .as_ref()
                .unwrap()
                .downcast_ref::<HeatData>()
                .unwrap()
                .N;
            println!(
                " {:4}  {}  {}  {}  {}   {:2}  {:3}",
                iout,
                fmt_e(olddt, 19, 15),
                fmt_e(newdt, 19, 15),
                fmt_e((N_VDotProd(&y, &y) / cur_n as f64).sqrt(), 19, 15),
                cur_n,
                nni,
                nli
            );
        }
        nni_tot += nni;
        nli_tot += nli;

        /* output results and current mesh to disk */
        {
            let udata = arkode_mem
                .user_data
                .as_ref()
                .unwrap()
                .downcast_ref::<HeatData>()
                .unwrap();
            for i in 0..udata.N as usize {
                let _ = write!(ufid, " {}", fmt_e(y.data[i], 0, 16));
            }
            let _ = writeln!(ufid);
            for i in 0..udata.N as usize {
                let _ = write!(xfid, " {}", fmt_e(udata.x[i], 0, 16));
            }
            let _ = writeln!(xfid);
        }

        /* adapt the spatial mesh; create N_Vector of new length; project
           solution onto new mesh; swap mesh/solution into place */
        let mut Nnew: i64 = 0;
        let mut y2;
        {
            let udata = arkode_mem
                .user_data
                .as_mut()
                .unwrap()
                .downcast_mut::<HeatData>()
                .unwrap();

            let xnew = match adapt_mesh(&y, &mut Nnew, udata) {
                Some(xn) => xn,
                None => {
                    eprintln!("\nSUNDIALS_ERROR: ark_adapt() failed - returned NULL pointer\n");
                    return;
                }
            };

            /* create N_Vector of new length */
            y2 = NVector::new(Nnew as usize);

            /* project solution onto new mesh */
            let flag = project(udata.N, &udata.x, &y, Nnew, &xnew, &mut y2);
            if check_flag(flag, "project") {
                return;
            }

            /* delete old vector and old mesh; swap x and xnew so that the new
               mesh is stored in the udata structure; store size of new mesh */
            udata.x = xnew;
            udata.N = Nnew;
        }

        /* swap y and y2 so that y holds new solution */
        y = y2;

        /* call ARKodeResize to notify integrator of change in mesh */
        let flag = ARKodeResize(&mut arkode_mem, &y, hscale, t, None, &mut no_resize_data);
        if check_flag(flag, "ARKodeResize") {
            return;
        }

        /* destroy and re-allocate linear solver memory; reattach to ARKODE
           interface (note: C passes the ORIGINAL N as the PCG maxl here) */
        let ls = SUNLinSol_PCG(&y, 0, N as i32, &ctx);
        let flag = ARKodeSetLinearSolver(&mut arkode_mem, ls, None);
        if check_flag(flag, "ARKodeSetLinearSolver") {
            return;
        }
        let flag = ARKodeSetJacTimes(&mut arkode_mem, None, Some(jac));
        if check_flag(flag, "ARKodeSetJacTimes") {
            return;
        }
    }
    println!(
        " ----------------------------------------------------------------------------------------"
    );

    /* print some final statistics */
    println!(" Final solver statistics:");
    println!("   Total number of time steps = {}", iout);
    println!("   Total nonlinear iterations = {}", nni_tot);
    println!("   Total linear iterations    = {}\n", nli_tot);

    /* Clean up and return with successful completion */
    drop(ufid);
    drop(xfid);
    drop(y); /* Free vectors */
    let mut slot = Some(arkode_mem);
    ARKodeFree(&mut slot); /* Free integrator memory */
}
