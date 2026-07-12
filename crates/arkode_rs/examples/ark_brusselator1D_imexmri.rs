/* -----------------------------------------------------------------------------
 * Translation of examples/arkode/C_serial/ark_brusselator1D_imexmri.c
 * (SUNDIALS 7.7.0).
 *
 * Example problem:
 *
 * The following test simulates a brusselator problem from chemical
 * kinetics.  This is n PDE system with 3 components, Y = [u,v,w],
 * satisfying the equations,
 *    u_t = du*u_xx - au*u_x +  a - (w+1)*u + v*u^2
 *    v_t = dv*v_xx - av*v_x +  w*u - v*u^2
 *    w_t = dw*w_xx - aw*w_x + (b-w)/ep - w*u
 * for t in [0, 10], x in [0, 1], with initial conditions
 *    u(0,x) =  a  + 0.1*sin(pi*x)
 *    v(0,x) = b/a + 0.1*sin(pi*x)
 *    w(0,x) =  b  + 0.1*sin(pi*x),
 * and with stationary boundary conditions.
 *
 * This program solves the problem with multiple solvers listed below.
 * We select method to used based on solve_type input:
 * 0. MIS with third order dirk inner
 * 1. 5th order dirk method for reference solution
 * 2. MRI-GARK34a with erk inner
 * 3. MRI-GARK34a with dirk inner
 * 4. IMEX-MRI3b with erk inner
 * 5. IMEX-MRI3b with dirk inner
 * 6. IMEX-MRI4 with erk inner
 * 7. IMEX-MRI4 with dirk inner
 *
 * We use Newton iteration with the SUNBAND linear solver and a user
 * supplied Jacobian routine for nonlinear solves.
 *
 * This program solves the problem with the MRI stepper. 10 outputs
 * are printed at equal intervals, and run statistics are printed at
 * the end.
 *
 * Note on ownership: the C example keeps its own inner_arkode_mem
 * pointer alongside the inner stepper; in this port the stepper owns
 * the wrapped inner integrator (and the outer integrator owns the
 * stepper), so the final fast statistics are read by borrowing the
 * inner integrator back out of the outer step memory.  The shared C
 * udata pointer becomes one identical UserData clone per integrator
 * (the problem data is immutable during the run).
 * ---------------------------------------------------------------------------*/
#![allow(non_snake_case)]

use std::io::Write;

use arkode_rs::arkode::{ARKodeCreateMRIStepInnerStepper, ARKodeEvolve, ARKodeFree};
use arkode_rs::arkode::ARKodeSStolerances;
use arkode_rs::arkode_arkstep::ARKStepCreate;
use arkode_rs::arkode_arkstep_io::{ARKStepSetTableNum, ARKStepSetTables};
use arkode_rs::arkode_butcher::ARKodeButcherTable_Alloc;
use arkode_rs::arkode_butcher_dirk::ARKODE_CASH_5_3_4;
use arkode_rs::arkode_io::{
    ARKodeGetNonlinSolvStats, ARKodeGetNumRhsEvals, ARKodeGetNumSteps, ARKodeSetFixedStep,
    ARKodeSetMaxNonlinIters, ARKodeSetMaxNumSteps, ARKodeSetOrder, ARKodeSetUserData,
};
use arkode_rs::arkode_ls::{ARKodeGetNumJacEvals, ARKodeSetJacFn, ARKodeSetLinearSolver};
use arkode_rs::arkode_mri_tables::{
    MRIStepCoupling_LoadTable, MRIStepCoupling_MIStoMRI, ARKODE_IMEX_MRI_GARK3b,
    ARKODE_IMEX_MRI_GARK4, ARKODE_MRI_GARK_ESDIRK34a,
};
use arkode_rs::arkode_mristep::MRIStepCreate;
use arkode_rs::arkode_mristep_impl::ARKodeMRIStepMem;
use arkode_rs::arkode_mristep_io::MRIStepSetCoupling;
use arkode_rs::sundials_utils::{fmt_e, fmt_f, fmt_g};
use arkode_rs::*;

/* Define some constants */
const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;

/* accessor between (x,v) location and 1D NVector array */
fn IDX(x: usize, v: usize) -> usize {
    3 * x + v
}

/* user data structure */
#[derive(Clone)]
struct UData {
    N: i64,  /* number of intervals     */
    dx: f64, /* mesh spacing            */
    a: f64,  /* constant forcing on u   */
    b: f64,  /* steady-state value of w */
    pi: f64, /* value of pi             */
    du: f64, /* diffusion coeff for u   */
    dv: f64, /* diffusion coeff for v   */
    dw: f64, /* diffusion coeff for w   */
    au: f64, /* advection coeff for u   */
    av: f64, /* advection coeff for v   */
    aw: f64, /* advection coeff for w   */
    ep: f64, /* stiffness parameter     */
}

/*------------------------------------
 * Functions called by the integrator
 *------------------------------------*/

/* ff routine to compute the fast portion of the ODE RHS. */
fn ff(_t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let udata = user_data.as_ref().unwrap().downcast_ref::<UData>().unwrap();
    let N = udata.N; /* set variable shortcuts */
    let a = udata.a;
    let b = udata.b;
    let ep = udata.ep;

    N_VConst(0.0, ydot); /* initialize ydot to zero */

    /* iterate over domain, computing all equations */
    for i in 1..(N - 1) as usize {
        /* set shortcuts */
        let u = y.data[IDX(i, 0)];
        let v = y.data[IDX(i, 1)];
        let w = y.data[IDX(i, 2)];

        /* Fill in ODE RHS for u */
        ydot.data[IDX(i, 0)] = a - (w + ONE) * u + v * u * u;

        /* Fill in ODE RHS for v */
        ydot.data[IDX(i, 1)] = w * u - v * u * u;

        /* Fill in ODE RHS for w */
        ydot.data[IDX(i, 2)] = (b - w) / ep - w * u;
    }

    /* enforce stationary boundaries */
    ydot.data[IDX(0, 0)] = 0.0;
    ydot.data[IDX(0, 1)] = 0.0;
    ydot.data[IDX(0, 2)] = 0.0;
    ydot.data[IDX((N - 1) as usize, 0)] = 0.0;
    ydot.data[IDX((N - 1) as usize, 1)] = 0.0;
    ydot.data[IDX((N - 1) as usize, 2)] = 0.0;

    /* Return with success */
    0
}

/* fse routine to compute the slow-explicit portion of the ODE RHS function. */
fn fse(_t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let udata = user_data.as_ref().unwrap().downcast_ref::<UData>().unwrap();
    let N = udata.N; /* set variable shortcuts */
    let au = udata.au;
    let av = udata.av;
    let aw = udata.aw;
    let dx = udata.dx;

    N_VConst(0.0, ydot); /* initialize ydot to zero */

    /* iterate over domain, computing all equations */
    let auconst = -au / 2.0 / dx;
    let avconst = -av / 2.0 / dx;
    let awconst = -aw / 2.0 / dx;
    for i in 1..(N - 1) as usize {
        /* set shortcuts */
        let ul = y.data[IDX(i - 1, 0)];
        let ur = y.data[IDX(i + 1, 0)];
        let vl = y.data[IDX(i - 1, 1)];
        let vr = y.data[IDX(i + 1, 1)];
        let wl = y.data[IDX(i - 1, 2)];
        let wr = y.data[IDX(i + 1, 2)];

        /* Fill in ODE RHS for u */
        ydot.data[IDX(i, 0)] = (ur - ul) * auconst;

        /* Fill in ODE RHS for v */
        ydot.data[IDX(i, 1)] = (vr - vl) * avconst;

        /* Fill in ODE RHS for w */
        ydot.data[IDX(i, 2)] = (wr - wl) * awconst;
    }

    /* enforce stationary boundaries */
    ydot.data[IDX(0, 0)] = 0.0;
    ydot.data[IDX(0, 1)] = 0.0;
    ydot.data[IDX(0, 2)] = 0.0;
    ydot.data[IDX((N - 1) as usize, 0)] = 0.0;
    ydot.data[IDX((N - 1) as usize, 1)] = 0.0;
    ydot.data[IDX((N - 1) as usize, 2)] = 0.0;

    /* Return with success */
    0
}

/* fsi routine to compute the slow-implicit portion of the ODE RHS. */
fn fsi(_t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let udata = user_data.as_ref().unwrap().downcast_ref::<UData>().unwrap();
    let N = udata.N; /* set variable shortcuts */
    let du = udata.du;
    let dv = udata.dv;
    let dw = udata.dw;
    let dx = udata.dx;

    N_VConst(0.0, ydot); /* initialize ydot to zero */

    /* iterate over domain, computing all equations */
    let duconst = du / dx / dx;
    let dvconst = dv / dx / dx;
    let dwconst = dw / dx / dx;
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
        ydot.data[IDX(i, 0)] = (ul - 2.0 * u + ur) * duconst;

        /* Fill in ODE RHS for v */
        ydot.data[IDX(i, 1)] = (vl - 2.0 * v + vr) * dvconst;

        /* Fill in ODE RHS for w */
        ydot.data[IDX(i, 2)] = (wl - 2.0 * w + wr) * dwconst;
    }

    /* enforce stationary boundaries */
    ydot.data[IDX(0, 0)] = 0.0;
    ydot.data[IDX(0, 1)] = 0.0;
    ydot.data[IDX(0, 2)] = 0.0;
    ydot.data[IDX((N - 1) as usize, 0)] = 0.0;
    ydot.data[IDX((N - 1) as usize, 1)] = 0.0;
    ydot.data[IDX((N - 1) as usize, 2)] = 0.0;

    /* Return with success */
    0
}

/* fs routine to compute the slow portion of the ODE RHS. */
fn fs(_t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let udata = user_data.as_ref().unwrap().downcast_ref::<UData>().unwrap();
    let N = udata.N; /* set variable shortcuts */
    let du = udata.du;
    let dv = udata.dv;
    let dw = udata.dw;
    let au = udata.au;
    let av = udata.av;
    let aw = udata.aw;
    let dx = udata.dx;

    N_VConst(0.0, ydot); /* initialize ydot to zero */

    /* iterate over domain, computing all equations */
    let duconst = du / dx / dx;
    let dvconst = dv / dx / dx;
    let dwconst = dw / dx / dx;
    let auconst = -au / TWO / dx;
    let avconst = -av / TWO / dx;
    let awconst = -aw / TWO / dx;
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
        ydot.data[IDX(i, 0)] = (ul - TWO * u + ur) * duconst + (ur - ul) * auconst;

        /* Fill in ODE RHS for v */
        ydot.data[IDX(i, 1)] = (vl - TWO * v + vr) * dvconst + (vr - vl) * avconst;

        /* Fill in ODE RHS for w */
        ydot.data[IDX(i, 2)] = (wl - TWO * w + wr) * dwconst + (wr - wl) * awconst;
    }

    /* enforce stationary boundaries */
    ydot.data[IDX(0, 0)] = 0.0;
    ydot.data[IDX(0, 1)] = 0.0;
    ydot.data[IDX(0, 2)] = 0.0;
    ydot.data[IDX((N - 1) as usize, 0)] = 0.0;
    ydot.data[IDX((N - 1) as usize, 1)] = 0.0;
    ydot.data[IDX((N - 1) as usize, 2)] = 0.0;

    /* Return with success */
    0
}

/* f routine to compute the full ODE RHS function. */
fn f(_t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let udata = user_data.as_ref().unwrap().downcast_ref::<UData>().unwrap();
    let N = udata.N; /* set variable shortcuts */
    let a = udata.a;
    let b = udata.b;
    let ep = udata.ep;
    let du = udata.du;
    let dv = udata.dv;
    let dw = udata.dw;
    let au = udata.au;
    let av = udata.av;
    let aw = udata.aw;
    let dx = udata.dx;

    N_VConst(0.0, ydot); /* initialize ydot to zero */

    /* iterate over domain, computing all equations */
    let duconst = du / dx / dx;
    let dvconst = dv / dx / dx;
    let dwconst = dw / dx / dx;
    let auconst = -au / 2.0 / dx;
    let avconst = -av / 2.0 / dx;
    let awconst = -aw / 2.0 / dx;
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
        ydot.data[IDX(i, 0)] =
            (ul - 2.0 * u + ur) * duconst + (ur - ul) * auconst + a - (w + 1.0) * u + v * u * u;

        /* Fill in ODE RHS for v */
        ydot.data[IDX(i, 1)] =
            (vl - 2.0 * v + vr) * dvconst + (vr - vl) * avconst + w * u - v * u * u;

        /* Fill in ODE RHS for w */
        ydot.data[IDX(i, 2)] =
            (wl - 2.0 * w + wr) * dwconst + (wr - wl) * awconst + (b - w) / ep - w * u;
    }

    /* enforce stationary boundaries */
    ydot.data[IDX(0, 0)] = 0.0;
    ydot.data[IDX(0, 1)] = 0.0;
    ydot.data[IDX(0, 2)] = 0.0;
    ydot.data[IDX((N - 1) as usize, 0)] = 0.0;
    ydot.data[IDX((N - 1) as usize, 1)] = 0.0;
    ydot.data[IDX((N - 1) as usize, 2)] = 0.0;

    /* Return with success */
    0
}

/* Placeholder function of zeroes */
fn f0(_t: f64, _y: &NVector, ydot: &mut NVector, _user_data: &mut UserData) -> i32 {
    N_VConst(ZERO, ydot);
    0
}

/* Jf routine to compute Jacobian of the fast portion of the ODE RHS */
#[allow(clippy::too_many_arguments)]
fn Jf(
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

    /* Add in the Jacobian of the reaction terms matrix */
    reaction_jac(1.0, y, j, udata);

    /* Return with success */
    0
}

/* Jsi routine to compute the Jacobian of the slow-implicit portion of the
   ODE RHS. */
#[allow(clippy::too_many_arguments)]
fn Jsi(
    _t: f64,
    _y: &NVector,
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

    /* Return with success */
    0
}

/* Js routine to compute the Jacobian of the slow portion of ODE RHS. */
#[allow(clippy::too_many_arguments)]
fn Js(
    _t: f64,
    _y: &NVector,
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

    /* Add Jacobian of the advection terms */
    advection_jac(1.0, j, udata);

    /* Return with success */
    0
}

/* Jac routine to compute the Jacobian of the full ODE RHS. */
#[allow(clippy::too_many_arguments)]
fn Jac(
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

    /* Add Jacobian of the advection terms */
    advection_jac(1.0, j, udata);

    /* Add in the Jacobian of the reaction terms matrix */
    reaction_jac(1.0, y, j, udata);

    /* Return with success */
    0
}

/*-------------------------------
 * Private helper functions
 *-------------------------------*/

/* Set the initial condition */
fn SetIC(y: &mut NVector, udata: &UData) -> i32 {
    let N = udata.N; /* set variable shortcuts */
    let a = udata.a;
    let b = udata.b;
    let dx = udata.dx;
    let pi = udata.pi;

    /* Set initial conditions into y */
    for i in 0..N as usize {
        y.data[IDX(i, 0)] = a + 0.1 * (pi * i as f64 * dx).sin(); /* u */
        y.data[IDX(i, 1)] = b / a + 0.1 * (pi * i as f64 * dx).sin(); /* v */
        y.data[IDX(i, 2)] = b + 0.1 * (pi * i as f64 * dx).sin(); /* w */
    }

    /* Return with success */
    0
}

/* Routine to compute the Jacobian matrix from fse(t,y), scaled by the factor c.
   We add the result into Jac and do not erase what was already there */
fn advection_jac(c: f64, jac: &mut SUNMatrix, udata: &UData) -> i32 {
    /* Set shortcuts */
    let N = udata.N;
    let dx = udata.dx;
    let au = udata.au;
    let av = udata.av;
    let aw = udata.aw;
    let auconst = -au / TWO / dx;
    let avconst = -av / TWO / dx;
    let awconst = -aw / TWO / dx;

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

        jm.set(r0, cl0, jm.get(r0, cl0) + -c * auconst);
        jm.set(r1, cl1, jm.get(r1, cl1) + -c * avconst);
        jm.set(r2, cl2, jm.get(r2, cl2) + -c * awconst);
        jm.set(r0, cr0, jm.get(r0, cr0) + c * auconst);
        jm.set(r1, cr1, jm.get(r1, cr1) + c * avconst);
        jm.set(r2, cr2, jm.get(r2, cr2) + c * awconst);
    }

    /* Return with success */
    0
}

/* Routine to compute the stiffness matrix from (L*y), scaled by the factor c.
   We add the result into Jac and do not erase what was already there */
fn laplace_matrix(c: f64, jac: &mut SUNMatrix, udata: &UData) -> i32 {
    let N = udata.N; /* set shortcuts */
    let dx = udata.dx;

    let jm = match jac {
        SUNMatrix::Band(m) => m,
        _ => return 1,
    };

    /* iterate over intervals, filling in Jacobian of (L*y) */
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

    /* Return with success */
    0
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
        jm.set(r0, r0, jm.get(r0, r0) + c * (TWO * u * v - (w + ONE)));
        jm.set(r1, r0, jm.get(r1, r0) + c * (w - TWO * u * v));
        jm.set(r2, r0, jm.get(r2, r0) + c * (-w));

        /* all vars wrt v */
        jm.set(r0, r1, jm.get(r0, r1) + c * (u * u));
        jm.set(r1, r1, jm.get(r1, r1) + c * (-u * u));

        /* all vars wrt w */
        jm.set(r0, r2, jm.get(r0, r2) + c * (-u));
        jm.set(r1, r2, jm.get(r1, r2) + c * u);
        jm.set(r2, r2, jm.get(r2, r2) + c * (-ONE / ep - u));
    }

    /* Return with success */
    0
}

/* Check if a SUNDIALS function returned a negative flag */
fn check_retval(retval: i32, funcname: &str) -> bool {
    if retval < 0 {
        eprintln!("\nSUNDIALS_ERROR: {}() failed with flag = {}\n", funcname, retval);
        return true;
    }
    false
}

/* Main Program */
fn main() {
    /* general problem parameters */
    let T0: f64 = ZERO; /* initial time                 */
    let Tf: f64 = 10.0; /* final time                   */
    let Nt: i32 = 10; /* total number of output times */
    let dTout: f64 = (Tf - T0) / Nt as f64; /* time between outputs */
    let Nvar: i64 = 3; /* number of solution fields    */
    let N: i64 = 101; /* spatial mesh size            */
    let m: f64 = 10.0; /* time-scale separation factor */
    let dx: f64 = ONE / (N - 1) as f64; /* set spatial mesh spacing     */
    let a: f64 = 0.6; /* problem parameters           */
    let b: f64 = 2.0;
    let pi: f64 = 4.0 * (ONE).atan();
    let du: f64 = 0.01;
    let dv: f64 = 0.01;
    let dw: f64 = 0.01;
    let au: f64 = -0.001;
    let av: f64 = -0.001;
    let aw: f64 = -0.001;
    let ep: f64 = 1.0e-2; /* stiffness parameter          */
    let mut reltol: f64 = 1.0e-12; /* tolerances                   */
    let mut abstol: f64 = 1.0e-14;

    /* Create the SUNDIALS context object for this simulation. */
    let ctx = SUNContext_Create();

    /*
     * Initialization
     */

    /* Retrieve the command-line options: solve_type h */
    let argv: Vec<String> = std::env::args().collect();
    if argv.len() < 2 {
        println!("ERROR: enter solve_type and hs ");
        std::process::exit(-1);
    }
    let solve_type: i64 = argv[1].parse().unwrap_or(0);
    let hs: f64 = argv[2].parse().unwrap_or(0.0);

    /* Check arguments for validity */
    /*   0 <= solve_type <= 7       */
    /*   h > 0                      */
    if !(0..=7).contains(&solve_type) {
        println!("ERROR: solve_type be an integer in [0,7] ");
        std::process::exit(-1);
    }
    let implicit_slow = solve_type > 1;
    let imex_slow = solve_type > 3;
    if hs <= ZERO {
        println!("ERROR: hs must be in positive");
        std::process::exit(-1);
    }
    let hf = hs / m;
    let NEQ = N * Nvar;

    /* Initial problem output */
    println!("\n1D Advection-Diffusion-Reaction (Brusselator) test problem:");
    println!("    time domain:  ({},{}]", fmt_g(T0, 0, 6), fmt_g(Tf, 0, 6));
    println!("    hs = {}", fmt_g(hs, 0, 6));
    println!("    hf = {}", fmt_g(hf, 0, 6));
    println!("    m  = {}", fmt_g(m, 0, 6));
    println!("    N  = {},  NEQ = {}", N, NEQ);
    println!("    dx = {}", fmt_g(dx, 0, 6));
    println!(
        "    problem parameters:  a = {},  b = {},  ep = {}",
        fmt_g(a, 0, 6),
        fmt_g(b, 0, 6),
        fmt_g(ep, 0, 6)
    );
    println!(
        "    diffusion coefficients:  du = {},  dv = {},  dw = {}",
        fmt_g(du, 0, 6),
        fmt_g(dv, 0, 6),
        fmt_g(dw, 0, 6)
    );
    println!(
        "    advection coefficients:  au = {},  av = {},  aw = {}",
        fmt_g(au, 0, 6),
        fmt_g(av, 0, 6),
        fmt_g(aw, 0, 6)
    );

    match solve_type {
        0 => {
            println!("    solver: exp-3/dirk-3 (MIS / ESDIRK-3-3)\n");
            println!("    reltol = {},  abstol = {}\n", fmt_e(reltol, 0, 2), fmt_e(abstol, 0, 2));
        }
        1 => {
            reltol = SUNMAX(hs * hs * hs * hs * hs, 1e-14);
            abstol = 1e-14;
            println!("    solver: none/dirk-5 (no slow, default 5th order dirk fast)\n");
            println!("    reltol = {},  abstol = {}\n", fmt_e(reltol, 0, 2), fmt_e(abstol, 0, 2));
        }
        2 => {
            println!(
                "    solver: dirk-3/exp-3 (MRI-GARK-ESDIRK34a / ERK-3-3) -- solve decoupled\n"
            );
            println!("    reltol = {},  abstol = {}\n", fmt_e(reltol, 0, 2), fmt_e(abstol, 0, 2));
        }
        3 => {
            println!(
                "    solver: dirk-3/dirk-3 (MRI-GARK-ESDIRK34a / ESDIRK-3-3) -- solve decoupled\n"
            );
            println!("    reltol = {},  abstol = {}\n", fmt_e(reltol, 0, 2), fmt_e(abstol, 0, 2));
        }
        4 => {
            println!("    solver: ars343/exp-3 (IMEX-MRI3b / ERK-3-3) -- solve decoupled\n");
            println!("    reltol = {},  abstol = {}\n", fmt_e(reltol, 0, 2), fmt_e(abstol, 0, 2));
        }
        5 => {
            println!("    solver: ars343/dirk-3 (IMEX-MRI3b / ESDIRK-3-3) -- solve decoupled\n");
            println!("    reltol = {},  abstol = {}\n", fmt_e(reltol, 0, 2), fmt_e(abstol, 0, 2));
        }
        6 => {
            println!("    solver: imexark4/exp-4 (IMEX-MRI4 / ERK-4-4) -- solve decoupled\n");
            println!("    reltol = {},  abstol = {}\n", fmt_e(reltol, 0, 2), fmt_e(abstol, 0, 2));
        }
        7 => {
            println!(
                "    solver: imexark4/dirk-4 (IMEX-MRI4 / CASH(5,3,4)-DIRK ) -- solve decoupled\n"
            );
            println!("    reltol = {},  abstol = {}\n", fmt_e(reltol, 0, 2), fmt_e(abstol, 0, 2));
        }
        _ => unreachable!(),
    }

    /* store the inputs in the UserData structure */
    let udata = UData {
        N,
        a,
        b,
        du,
        dv,
        dw,
        au,
        av,
        aw,
        ep,
        pi,
        dx,
    };

    /* Create solution vector and set initial condition */
    let mut y = N_VNew_Serial(NEQ, &ctx);
    let retval = SetIC(&mut y, &udata);
    if check_retval(retval, "SetIC") {
        return;
    }

    /* Create vector masks and set mask array values for each solution
       component */
    let mut umask = N_VClone(&y);
    let mut vmask = N_VClone(&y);
    let mut wmask = N_VClone(&y);

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

    /* Initialize the fast integrator. Specify the fast right-hand side
       function in y'=fs(t,y)+ff(t,y) = fse(t,y)+fsi(t,y)+ff(t,y), the initial
       time T0, and the initial dependent variable vector y. */
    let mut inner_arkode_mem = match solve_type {
        0 | 3 | 5 => {
            /* esdirk-3-3 fast solver */
            let mut inner = ARKStepCreate(None, Some(ff), T0, &y, &ctx).expect("ARKStepCreate");
            let mut B = ARKodeButcherTable_Alloc(3, false).expect("ARKodeButcherTable_Alloc");
            let beta = (3.0f64).sqrt() / 6.0 + 0.5;
            let gamma = (-ONE / 8.0) * ((3.0f64).sqrt() + ONE);
            B.A[1][0] = 4.0 * gamma + TWO * beta;
            B.A[1][1] = ONE - 4.0 * gamma - TWO * beta;
            B.A[2][0] = 0.5 - beta - gamma;
            B.A[2][1] = gamma;
            B.A[2][2] = beta;
            B.b[0] = ONE / 6.0;
            B.b[1] = ONE / 6.0;
            B.b[2] = TWO / 3.0;
            B.c[1] = ONE;
            B.c[2] = 0.5;
            B.q = 3;
            let retval = ARKStepSetTables(&mut inner, 3, 0, Some(&B), None);
            if check_retval(retval, "ARKStepSetTables") {
                return;
            }

            /* Initialize matrix and linear solver data structures */
            let Af = SUNBandMatrix(NEQ, 4, 4, &ctx);
            let LSf = SUNLinSol_Band(&y, &Af, &ctx);

            /* Specify fast tolerances */
            let retval = ARKodeSStolerances(&mut inner, reltol, abstol);
            if check_retval(retval, "ARKodeSStolerances") {
                return;
            }

            /* Attach matrix and linear solver */
            let retval = ARKodeSetLinearSolver(&mut inner, LSf, Some(Af));
            if check_retval(retval, "ARKodeSetLinearSolver") {
                return;
            }

            /* Set max number of nonlinear iters */
            let retval = ARKodeSetMaxNonlinIters(&mut inner, 10);
            if check_retval(retval, "ARKodeSetMaxNonlinIters") {
                return;
            }

            /* Set the Jacobian routine */
            let retval = ARKodeSetJacFn(&mut inner, Some(Jf));
            if check_retval(retval, "ARKodeSetJacFn") {
                return;
            }
            inner
        }
        1 => {
            /* dirk 5th order fast solver (full problem) */
            let mut inner = ARKStepCreate(None, Some(f), T0, &y, &ctx).expect("ARKStepCreate");

            /* Set method order to use */
            let retval = ARKodeSetOrder(&mut inner, 5);
            if check_retval(retval, "ARKodeSetOrder") {
                return;
            }

            /* Initialize matrix and linear solver data structures */
            let Af = SUNBandMatrix(NEQ, 4, 4, &ctx);
            let LSf = SUNLinSol_Band(&y, &Af, &ctx);

            /* Specify fast tolerances */
            let retval = ARKodeSStolerances(&mut inner, reltol, abstol);
            if check_retval(retval, "ARKodeSStolerances") {
                return;
            }

            /* Attach matrix and linear solver */
            let retval = ARKodeSetLinearSolver(&mut inner, LSf, Some(Af));
            if check_retval(retval, "ARKodeSetLinearSolver") {
                return;
            }

            /* Set the Jacobian routine */
            let retval = ARKodeSetJacFn(&mut inner, Some(Jac));
            if check_retval(retval, "ARKodeSetJacFn") {
                return;
            }
            inner
        }
        2 | 4 => {
            /* erk-3-3 fast solver */
            let mut inner = ARKStepCreate(Some(ff), None, T0, &y, &ctx).expect("ARKStepCreate");
            let mut B = ARKodeButcherTable_Alloc(3, true).expect("ARKodeButcherTable_Alloc");
            B.A[1][0] = 0.5;
            B.A[2][0] = -ONE;
            B.A[2][1] = TWO;
            B.b[0] = ONE / 6.0;
            B.b[1] = TWO / 3.0;
            B.b[2] = ONE / 6.0;
            B.d.as_mut().unwrap()[1] = ONE;
            B.c[1] = 0.5;
            B.c[2] = ONE;
            B.q = 3;
            B.p = 2;
            let retval = ARKStepSetTables(&mut inner, 3, 2, None, Some(&B));
            if check_retval(retval, "ARKStepSetTables") {
                return;
            }
            inner
        }
        6 => {
            /* erk-4-4 fast solver */
            let mut inner = ARKStepCreate(Some(ff), None, T0, &y, &ctx).expect("ARKStepCreate");
            let mut B = ARKodeButcherTable_Alloc(4, false).expect("ARKodeButcherTable_Alloc");
            B.A[1][0] = 0.5;
            B.A[2][1] = 0.5;
            B.A[3][2] = ONE;
            B.b[0] = ONE / 6.0;
            B.b[1] = ONE / 3.0;
            B.b[2] = ONE / 3.0;
            B.b[3] = ONE / 6.0;
            B.c[1] = 0.5;
            B.c[2] = 0.5;
            B.c[3] = ONE;
            B.q = 4;
            let retval = ARKStepSetTables(&mut inner, 4, 0, None, Some(&B));
            if check_retval(retval, "ARKStepSetTables") {
                return;
            }
            inner
        }
        7 => {
            /* Cash(5,3,4)-SDIRK fast solver */
            let mut inner = ARKStepCreate(None, Some(ff), T0, &y, &ctx).expect("ARKStepCreate");

            /* Set fast method */
            let retval = ARKStepSetTableNum(&mut inner, ARKODE_CASH_5_3_4, -1);
            if check_retval(retval, "ARKStepSetTableNum") {
                return;
            }

            /* Initialize matrix and linear solver data structures */
            let Af = SUNBandMatrix(NEQ, 4, 4, &ctx);
            let LSf = SUNLinSol_Band(&y, &Af, &ctx);

            /* Specify fast tolerances */
            let retval = ARKodeSStolerances(&mut inner, reltol, abstol);
            if check_retval(retval, "ARKodeSStolerances") {
                return;
            }

            /* Attach matrix and linear solver */
            let retval = ARKodeSetLinearSolver(&mut inner, LSf, Some(Af));
            if check_retval(retval, "ARKodeSetLinearSolver") {
                return;
            }

            /* Set max number of nonlinear iters */
            let retval = ARKodeSetMaxNonlinIters(&mut inner, 10);
            if check_retval(retval, "ARKodeSetMaxNonlinIters") {
                return;
            }

            /* Set the Jacobian routine */
            let retval = ARKodeSetJacFn(&mut inner, Some(Jf));
            if check_retval(retval, "ARKodeSetJacFn") {
                return;
            }
            inner
        }
        _ => unreachable!(),
    };

    /* Attach user data to fast integrator */
    let retval = ARKodeSetUserData(&mut inner_arkode_mem, Some(Box::new(udata.clone())));
    if check_retval(retval, "ARKodeSetUserData") {
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

    /* Initialize the slow integrator. Specify the slow right-hand side
       function in y'=fs(t,y)+ff(t,y) = fse(t,y)+fsi(t,y)+ff(t,y), the initial
       time T0, the initial dependent variable vector y, and the fast
       integrator. */
    let mut arkode_mem = match solve_type {
        0 => {
            /* use MIS outer integrator default for MRIStep */
            MRIStepCreate(Some(fs), None, T0, &y, inner_stepper, &ctx).expect("MRIStepCreate")
        }
        1 => {
            /* no slow dynamics (use ERK-2-2) */
            let mut outer =
                MRIStepCreate(Some(f0), None, T0, &y, inner_stepper, &ctx).expect("MRIStepCreate");
            let mut B = ARKodeButcherTable_Alloc(2, false).expect("ARKodeButcherTable_Alloc");
            B.A[1][0] = TWO / 3.0;
            B.b[0] = 0.25;
            B.b[1] = 0.75;
            B.c[1] = TWO / 3.0;
            B.q = 2;
            let C = MRIStepCoupling_MIStoMRI(&B, 2, 0).expect("MRIStepCoupling_MIStoMRI");
            let retval = MRIStepSetCoupling(&mut outer, &C);
            if check_retval(retval, "MRIStepSetCoupling") {
                return;
            }
            outer
        }
        2 | 3 => {
            /* MRI-GARK-ESDIRK34a, solve-decoupled slow solver */
            let mut outer =
                MRIStepCreate(None, Some(fs), T0, &y, inner_stepper, &ctx).expect("MRIStepCreate");

            let C =
                MRIStepCoupling_LoadTable(ARKODE_MRI_GARK_ESDIRK34a).expect("MRIStepCoupling_LoadTable");

            let retval = MRIStepSetCoupling(&mut outer, &C);
            if check_retval(retval, "MRIStepSetCoupling") {
                return;
            }

            /* Initialize matrix and linear solver data structures */
            let As = SUNBandMatrix(NEQ, 4, 4, &ctx);
            let LSs = SUNLinSol_Band(&y, &As, &ctx);

            /* Specify tolerances */
            let retval = ARKodeSStolerances(&mut outer, reltol, abstol);
            if check_retval(retval, "ARKodeSStolerances") {
                return;
            }

            /* Attach matrix and linear solver */
            let retval = ARKodeSetLinearSolver(&mut outer, LSs, Some(As));
            if check_retval(retval, "ARKodeSetLinearSolver") {
                return;
            }

            /* Set the Jacobian routine */
            let retval = ARKodeSetJacFn(&mut outer, Some(Js));
            if check_retval(retval, "ARKodeSetJacFn") {
                return;
            }
            outer
        }
        4 | 5 => {
            /* IMEX-MRI-GARK3b, solve-decoupled slow solver */
            let mut outer = MRIStepCreate(Some(fse), Some(fsi), T0, &y, inner_stepper, &ctx)
                .expect("MRIStepCreate");

            let C =
                MRIStepCoupling_LoadTable(ARKODE_IMEX_MRI_GARK3b).expect("MRIStepCoupling_LoadTable");

            let retval = MRIStepSetCoupling(&mut outer, &C);
            if check_retval(retval, "MRIStepSetCoupling") {
                return;
            }

            /* Initialize matrix and linear solver data structures */
            let As = SUNBandMatrix(NEQ, 4, 4, &ctx);
            let LSs = SUNLinSol_Band(&y, &As, &ctx);

            /* Specify tolerances */
            let retval = ARKodeSStolerances(&mut outer, reltol, abstol);
            if check_retval(retval, "ARKodeSStolerances") {
                return;
            }

            /* Attach matrix and linear solver */
            let retval = ARKodeSetLinearSolver(&mut outer, LSs, Some(As));
            if check_retval(retval, "ARKodeSetLinearSolver") {
                return;
            }

            /* Set the Jacobian routine */
            let retval = ARKodeSetJacFn(&mut outer, Some(Jsi));
            if check_retval(retval, "ARKodeSetJacFn") {
                return;
            }
            outer
        }
        6 | 7 => {
            /* IMEX-MRI-GARK4, solve-decoupled slow solver */
            let mut outer = MRIStepCreate(Some(fse), Some(fsi), T0, &y, inner_stepper, &ctx)
                .expect("MRIStepCreate");

            let C =
                MRIStepCoupling_LoadTable(ARKODE_IMEX_MRI_GARK4).expect("MRIStepCoupling_LoadTable");

            let retval = MRIStepSetCoupling(&mut outer, &C);
            if check_retval(retval, "MRIStepSetCoupling") {
                return;
            }

            /* Initialize matrix and linear solver data structures */
            let As = SUNBandMatrix(NEQ, 4, 4, &ctx);
            let LSs = SUNLinSol_Band(&y, &As, &ctx);

            /* Specify tolerances */
            let retval = ARKodeSStolerances(&mut outer, reltol, abstol);
            if check_retval(retval, "ARKodeSStolerances") {
                return;
            }

            /* Attach matrix and linear solver */
            let retval = ARKodeSetLinearSolver(&mut outer, LSs, Some(As));
            if check_retval(retval, "ARKodeSetLinearSolver") {
                return;
            }

            /* Set the Jacobian routine */
            let retval = ARKodeSetJacFn(&mut outer, Some(Jsi));
            if check_retval(retval, "ARKodeSetJacFn") {
                return;
            }
            outer
        }
        _ => unreachable!(),
    };

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

    /* Set maximum number of steps taken by solver */
    let retval = ARKodeSetMaxNumSteps(&mut arkode_mem, 1000000);
    if check_retval(retval, "ARKodeSetMaxNumSteps") {
        return;
    }

    /*
     * Integrate ODE
     */

    /* output spatial mesh to disk */
    let mut fid = std::fs::File::create("bruss1D_mesh.txt").expect("fopen");
    for i in 0..N {
        let _ = writeln!(fid, "  {}", fmt_e(udata.dx * i as f64, 0, 16));
    }
    drop(fid);

    /* Open output stream for results */
    let mut ufid =
        std::fs::File::create(format!("bruss1D_u_{}_{}.txt", argv[1], argv[2])).expect("fopen");
    let mut vfid =
        std::fs::File::create(format!("bruss1D_v_{}_{}.txt", argv[1], argv[2])).expect("fopen");
    let mut wfid =
        std::fs::File::create(format!("bruss1D_w_{}_{}.txt", argv[1], argv[2])).expect("fopen");

    /* output initial condition to disk */
    for i in 0..N as usize {
        let _ = write!(ufid, " {}", fmt_e(y.data[IDX(i, 0)], 0, 16));
    }
    for i in 0..N as usize {
        let _ = write!(vfid, " {}", fmt_e(y.data[IDX(i, 1)], 0, 16));
    }
    for i in 0..N as usize {
        let _ = write!(wfid, " {}", fmt_e(y.data[IDX(i, 2)], 0, 16));
    }
    let _ = writeln!(ufid);
    let _ = writeln!(vfid);
    let _ = writeln!(wfid);

    /* Main time-stepping loop: calls ARKodeEvolve to perform the integration,
       then prints results.  Stops when the final time has been reached */
    let mut t = T0;
    let mut tout = T0 + dTout;
    println!("        t      ||u||_rms   ||v||_rms   ||w||_rms");
    println!("   ----------------------------------------------");
    for _iout in 0..Nt {
        /* call integrator */
        let retval = ARKodeEvolve(&mut arkode_mem, tout, &mut y, &mut t, ARK_NORMAL);
        if check_retval(retval, "ARKodeEvolve") {
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

        /* successful solve: update output time */
        tout += dTout;
        tout = if tout > Tf { Tf } else { tout };

        /* output results to disk */
        for i in 0..N as usize {
            let _ = write!(ufid, " {}", fmt_e(y.data[IDX(i, 0)], 0, 16));
        }
        for i in 0..N as usize {
            let _ = write!(vfid, " {}", fmt_e(y.data[IDX(i, 1)], 0, 16));
        }
        for i in 0..N as usize {
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

    /*
     * Finalize
     */

    /* Get some slow integrator statistics */
    let mut nsts: i64 = 0;
    let mut nfse_c: i64 = 0;
    let mut nfsi_c: i64 = 0;
    ARKodeGetNumSteps(&mut arkode_mem, &mut nsts);
    ARKodeGetNumRhsEvals(&mut arkode_mem, 0, &mut nfse_c);
    ARKodeGetNumRhsEvals(&mut arkode_mem, 1, &mut nfsi_c);

    /* Get some fast integrator statistics (borrow the inner integrator
       back out of the outer step memory) */
    let mut nstf: i64 = 0;
    let mut nffe: i64 = 0;
    let mut nffi: i64 = 0;
    let mut nnif: i64 = 0;
    let mut nncf: i64 = 0;
    let mut njef: i64 = 0;
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
        ARKodeGetNumSteps(inner, &mut nstf);
        ARKodeGetNumRhsEvals(inner, 0, &mut nffe);
        ARKodeGetNumRhsEvals(inner, 1, &mut nffi);

        /* fast integrator implicit solver statistics (read here while the
           inner integrator is borrowed; printed below in C's order) */
        if solve_type == 0 || solve_type == 1 || solve_type == 3 || solve_type == 5
            || solve_type == 7
        {
            ARKodeGetNonlinSolvStats(inner, &mut nnif, &mut nncf);
            ARKodeGetNumJacEvals(inner, &mut njef);
        }
    }

    /* Print some final statistics */
    println!("\nFinal Solver Statistics:");
    println!("   Slow Steps: nsts = {}", nsts);
    println!("   Fast Steps: nstf = {}", nstf);
    if imex_slow {
        if solve_type == 0 || solve_type == 1 || solve_type == 3 || solve_type == 5
            || solve_type == 7
        {
            println!("   Total RHS evals:  Fse = {}, Fsi = {},  Ff = {}", nfse_c, nfsi_c, nffi);
        } else {
            println!("   Total RHS evals:  Fse = {}, Fsi = {},  Ff = {}", nfse_c, nfsi_c, nffe);
        }
    } else if implicit_slow {
        if solve_type == 0 || solve_type == 1 || solve_type == 3 || solve_type == 5
            || solve_type == 7
        {
            println!("   Total RHS evals:  Fs = {},  Ff = {}", nfsi_c, nffi);
        } else {
            println!("   Total RHS evals:  Fs = {},  Ff = {}", nfsi_c, nffe);
        }
    } else if solve_type == 0 || solve_type == 1 || solve_type == 3 || solve_type == 5
        || solve_type == 7
    {
        println!("   Total RHS evals:  Fs = {},  Ff = {}", nfse_c, nffi);
    } else {
        println!("   Total RHS evals:  Fs = {},  Ff = {}", nfse_c, nffe);
    }

    /* Get/print slow integrator decoupled implicit solver statistics */
    if solve_type > 1 {
        let mut nnis: i64 = 0;
        let mut nncs: i64 = 0;
        let mut njes: i64 = 0;
        ARKodeGetNonlinSolvStats(&mut arkode_mem, &mut nnis, &mut nncs);
        ARKodeGetNumJacEvals(&mut arkode_mem, &mut njes);
        println!("   Slow Newton iters = {}", nnis);
        println!("   Slow Newton conv fails = {}", nncs);
        println!("   Slow Jacobian evals = {}", njes);
    }

    /* Print fast integrator implicit solver statistics */
    if solve_type == 0 || solve_type == 1 || solve_type == 3 || solve_type == 5 || solve_type == 7
    {
        println!("   Fast Newton iters = {}", nnif);
        println!("   Fast Newton conv fails = {}", nncf);
        println!("   Fast Jacobian evals = {}", njef);
    }

    /* Clean up and return with successful completion */
    drop(y); /* Free vectors */
    drop(umask);
    drop(vmask);
    drop(wmask);
    let mut slot = Some(arkode_mem);
    ARKodeFree(&mut slot);
}
