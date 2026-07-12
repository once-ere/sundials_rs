/* -----------------------------------------------------------------
 * Translation of examples/cvodes/serial/cvsAdvDiff_ASAi_bnd.c
 * (SUNDIALS 7.7.0).
 *
 * Adjoint sensitivity example problem:
 *
 * The following is a simple example problem with a banded Jacobian.
 * The problem is the semi-discrete form of the advection-diffusion
 * equation in 2-D:
 *   du/dt = d^2 u / dx^2 + .5 du/dx + d^2 u / dy^2
 * on the rectangle 0 <= x <= 2, 0 <= y <= 1, and the time
 * interval 0 <= t <= 1. Homogeneous Dirichlet boundary conditions
 * are posed, and the initial condition is:
 *   u(x,y,t=0) = x(2-x)y(1-y)exp(5xy).
 * The PDE is discretized on a uniform MX+2 by MY+2 grid with
 * central differencing, leaving an ODE system of size NEQ = MX*MY.
 *
 * Additionally, CVODES integrates backwards in time the
 * semi-discrete form of the adjoint PDE:
 *   d(lambda)/dt = - d^2(lambda) / dx^2 + 0.5 d(lambda) / dx
 *                  - d^2(lambda) / dy^2 - 1.0
 * with homogeneous Dirichlet boundary conditions and final
 * conditions lambda(x,y,t=t_final) = 0.0, whose solution at t = 0
 * represents the sensitivity of
 *   G = int_0^t_final int_x int _y u(t,x,y) dx dy dt
 * with respect to the initial conditions of the original problem.
 *
 * Translation note: the C program shares one UserData pointer
 * between the forward and backward problems; here each side owns an
 * identical GridData copy (the struct is constant after setup).
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]

use cvodes_rs::cvodea::{
    CVodeAdjInit, CVodeB, CVodeCreateB, CVodeF, CVodeGetB, CVodeInitB, CVodeSStolerancesB,
};
use cvodes_rs::cvodea_io::CVodeSetUserDataB;
use cvodes_rs::cvodes::{CVodeCreate, CVodeFree, CVodeInit, CVodeSStolerances};
use cvodes_rs::cvodes_io::CVodeSetUserData;
use cvodes_rs::cvodes_ls::{
    CVodeSetJacFn, CVodeSetJacFnB, CVodeSetLinearSolver, CVodeSetLinearSolverB,
};
use cvodes_rs::sundials_utils::fmt_e;
use cvodes_rs::*;

/* Problem Constants */

const XMAX: f64 = 2.0; /* domain boundaries             */
const YMAX: f64 = 1.0;
const MX: i64 = 40; /* mesh dimensions               */
const MY: i64 = 20;
const NEQ: i64 = MX * MY; /* number of equations           */
const ATOL: f64 = 1.0e-5;
const RTOLB: f64 = 1.0e-6;
const T0: f64 = 0.0; /* initial time                  */
const TOUT: f64 = 1.0; /* final time                    */
const NSTEP: i64 = 50; /* check point saved every NSTEP */

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;

/* User-defined vector access macro IJth */

/* IJth(vdata,i,j) references the element in the vdata array for
   u at mesh point (i,j), where 1 <= i <= MX, 1 <= j <= MY.
   The variables are ordered by the y index j, then by the x index i. */

macro_rules! IJth {
    ($vdata:expr, $i:expr, $j:expr) => {
        $vdata[(($j - 1) + ($i - 1) * MY) as usize]
    };
}

/* Type : GridData (the C UserData)
   contains grid constants */

#[derive(Clone)]
struct GridData {
    dx: f64,
    dy: f64,
    hdcoef: f64,
    hacoef: f64,
    vdcoef: f64,
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY CVODES
 *--------------------------------------------------------------------
 */

/*
 * f routine. right-hand side of forward ODE.
 */
fn f(_t: f64, u: &NVector, udot: &mut NVector, user_data: &mut UserData) -> i32 {
    let udata = &u.data;

    /* Extract needed constants from data */

    let data = user_data.as_ref().unwrap().downcast_ref::<GridData>().unwrap();
    let hordc = data.hdcoef;
    let horac = data.hacoef;
    let verdc = data.vdcoef;

    let dudata = &mut udot.data;

    /* Loop over all grid points. */

    for j in 1..=MY {
        for i in 1..=MX {
            /* Extract u at x_i, y_j and four neighboring points */

            let uij = IJth!(udata, i, j);
            let udn = if j == 1 { ZERO } else { IJth!(udata, i, j - 1) };
            let uup = if j == MY { ZERO } else { IJth!(udata, i, j + 1) };
            let ult = if i == 1 { ZERO } else { IJth!(udata, i - 1, j) };
            let urt = if i == MX { ZERO } else { IJth!(udata, i + 1, j) };

            /* Set diffusion and advection terms and load into udot */

            let hdiff = hordc * (ult - TWO * uij + urt);
            let hadv = horac * (urt - ult);
            let vdiff = verdc * (uup - TWO * uij + udn);
            IJth!(dudata, i, j) = hdiff + hadv + vdiff;
        }
    }

    0
}

/*
 * Jac function. Jacobian of forward ODE.
 */
#[allow(clippy::too_many_arguments)]
fn jac(
    _t: f64,
    _u: &NVector,
    _fu: &NVector,
    j_mat: &mut SUNMatrix,
    user_data: &mut UserData,
    _tmp1: &mut NVector,
    _tmp2: &mut NVector,
    _tmp3: &mut NVector,
) -> i32 {
    /*
      The components of f = udot that depend on u(i,j) are
      f(i,j), f(i-1,j), f(i+1,j), f(i,j-1), f(i,j+1), with
        df(i,j)/du(i,j) = -2 (1/dx^2 + 1/dy^2)
        df(i-1,j)/du(i,j) = 1/dx^2 + .25/dx  (if i > 1)
        df(i+1,j)/du(i,j) = 1/dx^2 - .25/dx  (if i < MX)
        df(i,j-1)/du(i,j) = 1/dy^2           (if j > 1)
        df(i,j+1)/du(i,j) = 1/dy^2           (if j < MY)
    */

    let data = user_data.as_ref().unwrap().downcast_ref::<GridData>().unwrap();
    let hordc = data.hdcoef;
    let horac = data.hacoef;
    let verdc = data.vdcoef;

    let jm = match j_mat {
        SUNMatrix::Band(m) => m,
        _ => return -1,
    };

    for j in 1..=MY {
        for i in 1..=MX {
            let k = j - 1 + (i - 1) * MY;

            /* set the kth column of J */

            jm.set(k, k, -TWO * (verdc + hordc));
            if i != 1 {
                jm.set(k - MY, k, hordc + horac);
            }
            if i != MX {
                jm.set(k + MY, k, hordc - horac);
            }
            if j != 1 {
                jm.set(k - 1, k, verdc);
            }
            if j != MY {
                jm.set(k + 1, k, verdc);
            }
        }
    }

    0
}

/*
 * fB function. Right-hand side of backward ODE.
 */
fn fB(_tB: f64, _u: &NVector, uB: &NVector, uBdot: &mut NVector, user_dataB: &mut UserData) -> i32 {
    let uBdata = &uB.data;

    /* Extract needed constants from data */

    let data = user_dataB.as_ref().unwrap().downcast_ref::<GridData>().unwrap();
    let hordc = data.hdcoef;
    let horac = data.hacoef;
    let verdc = data.vdcoef;

    let duBdata = &mut uBdot.data;

    /* Loop over all grid points. */

    for j in 1..=MY {
        for i in 1..=MX {
            /* Extract u at x_i, y_j and four neighboring points */

            let uBij = IJth!(uBdata, i, j);
            let uBdn = if j == 1 { ZERO } else { IJth!(uBdata, i, j - 1) };
            let uBup = if j == MY { ZERO } else { IJth!(uBdata, i, j + 1) };
            let uBlt = if i == 1 { ZERO } else { IJth!(uBdata, i - 1, j) };
            let uBrt = if i == MX { ZERO } else { IJth!(uBdata, i + 1, j) };

            /* Set diffusion and advection terms and load into udot */

            let hdiffB = hordc * (-uBlt + TWO * uBij - uBrt);
            let hadvB = horac * (uBrt - uBlt);
            let vdiffB = verdc * (-uBup + TWO * uBij - uBdn);
            IJth!(duBdata, i, j) = hdiffB + hadvB + vdiffB - ONE;
        }
    }

    0
}

/*
 * JacB function. Jacobian of backward ODE
 */
#[allow(clippy::too_many_arguments)]
fn jacB(
    _tB: f64,
    _u: &NVector,
    _uB: &NVector,
    _fuB: &NVector,
    jB_mat: &mut SUNMatrix,
    user_dataB: &mut UserData,
    _tmp1B: &mut NVector,
    _tmp2B: &mut NVector,
    _tmp3B: &mut NVector,
) -> i32 {
    /* The Jacobian of the adjoint system is: JB = -J^T */

    let data = user_dataB.as_ref().unwrap().downcast_ref::<GridData>().unwrap();
    let hordc = data.hdcoef;
    let horac = data.hacoef;
    let verdc = data.vdcoef;

    let jm = match jB_mat {
        SUNMatrix::Band(m) => m,
        _ => return -1,
    };

    for j in 1..=MY {
        for i in 1..=MX {
            let k = j - 1 + (i - 1) * MY;

            /* set the kth column of J */

            jm.set(k, k, TWO * (verdc + hordc));
            if i != 1 {
                jm.set(k - MY, k, -hordc + horac);
            }
            if i != MX {
                jm.set(k + MY, k, -hordc - horac);
            }
            if j != 1 {
                jm.set(k - 1, k, -verdc);
            }
            if j != MY {
                jm.set(k + 1, k, -verdc);
            }
        }
    }

    0
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * Set initial conditions in u vector
 */
fn SetIC(u: &mut NVector, data: &GridData) {
    /* Extract needed constants from data */

    let dx = data.dx;
    let dy = data.dy;

    /* Set pointer to data array in vector u. */

    let udata = &mut u.data;

    /* Load initial profile into u vector */

    for j in 1..=MY {
        let y = j as f64 * dy;
        for i in 1..=MX {
            let x = i as f64 * dx;
            IJth!(udata, i, j) = x * (XMAX - x) * y * (YMAX - y) * (5.0 * x * y).exp();
        }
    }
}

/*
 * Print results after backward integration
 */
fn PrintOutput(uB: &NVector, data: &GridData) {
    let mut x = ZERO;
    let mut y = ZERO;

    let dx = data.dx;
    let dy = data.dy;

    let uBdata = &uB.data;

    let mut uBmax = ZERO;
    for j in 1..=MY {
        for i in 1..=MX {
            let uBij = IJth!(uBdata, i, j);
            if uBij.abs() > uBmax {
                uBmax = uBij;
                x = i as f64 * dx;
                y = j as f64 * dy;
            }
        }
    }

    println!("\nMaximum sensitivity");
    println!("  lambda max = {}", fmt_e(uBmax, 0, 6));
    println!("at");
    println!("  x = {}\n  y = {}", fmt_e(x, 0, 6), fmt_e(y, 0, 6));
}

/* Check if a SUNDIALS function returned a negative value */
fn check_retval(retval: i32, funcname: &str) -> bool {
    if retval < 0 {
        eprintln!("\nSUNDIALS_ERROR: {}() failed with retval = {}\n", funcname, retval);
        return true;
    }
    false
}

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */
fn main() {
    /* Allocate and initialize user data memory */

    let dx = XMAX / (MX + 1) as f64;
    let dy = YMAX / (MY + 1) as f64;
    let data = GridData {
        dx,
        dy,
        hdcoef: ONE / (dx * dx),
        hacoef: 1.5 / (TWO * dx),
        vdcoef: ONE / (dy * dy),
    };

    /* Set the tolerances for the forward integration */
    let reltol = ZERO;
    let abstol = ATOL;

    /* Create the SUNDIALS simulation context that all SUNDIALS objects require */
    let sunctx = SUNContext_Create();

    /* Allocate u vector */
    let mut u = N_VNew_Serial(NEQ, &sunctx);

    /* Initialize u vector */
    SetIC(&mut u, &data);

    /* Create and allocate CVODES memory for forward run */

    println!("\nCreate and allocate CVODES memory for forward runs");

    let mut cvode_mem = CVodeCreate(CV_BDF, &sunctx);

    let retval = CVodeSetUserData(&mut cvode_mem, Some(Box::new(data.clone())));
    if check_retval(retval, "CVodeSetUserData") {
        return;
    }

    let retval = CVodeInit(&mut cvode_mem, f, T0, &u);
    if check_retval(retval, "CVodeInit") {
        return;
    }

    let retval = CVodeSStolerances(&mut cvode_mem, reltol, abstol);
    if check_retval(retval, "CVodeSStolerances") {
        return;
    }

    /* Create banded SUNMatrix for the forward problem */
    let a_mat = SUNBandMatrix(NEQ, MY, MY, &sunctx);

    /* Create banded SUNLinearSolver for the forward problem */
    let ls = SUNLinSol_Band(&u, &a_mat, &sunctx);

    /* Attach the matrix and linear solver */
    let retval = CVodeSetLinearSolver(&mut cvode_mem, ls, Some(a_mat));
    if check_retval(retval, "CVodeSetLinearSolver") {
        return;
    }

    /* Set the user-supplied Jacobian routine for the forward problem */
    let retval = CVodeSetJacFn(&mut cvode_mem, Some(jac));
    if check_retval(retval, "CVodeSetJacFn") {
        return;
    }

    /* Allocate global memory */

    println!("\nAllocate global memory");

    let retval = CVodeAdjInit(&mut cvode_mem, NSTEP, CV_HERMITE);
    if check_retval(retval, "CVodeAdjInit") {
        return;
    }

    /* Perform forward run */
    println!("\nForward integration");
    let mut t = 0.0;
    let mut ncheck: i32 = 0;
    let retval = CVodeF(&mut cvode_mem, TOUT, &mut u, &mut t, CV_NORMAL, &mut ncheck);
    if check_retval(retval, "CVodeF") {
        return;
    }

    println!("\nncheck = {}", ncheck);

    /* Set the tolerances for the backward integration */
    let reltolB = RTOLB;
    let abstolB = ATOL;

    /* Allocate uB */
    let mut uB = N_VNew_Serial(NEQ, &sunctx);
    /* Initialize uB = 0 */
    N_VConst(ZERO, &mut uB);

    /* Create and allocate CVODES memory for backward run */

    println!("\nCreate and allocate CVODES memory for backward run");

    let mut indexB: i32 = 0;
    let retval = CVodeCreateB(&mut cvode_mem, CV_BDF, &mut indexB);
    if check_retval(retval, "CVodeCreateB") {
        return;
    }

    let retval = CVodeSetUserDataB(&mut cvode_mem, indexB, Some(Box::new(data.clone())));
    if check_retval(retval, "CVodeSetUserDataB") {
        return;
    }

    let retval = CVodeInitB(&mut cvode_mem, indexB, fB, TOUT, &uB);
    if check_retval(retval, "CVodeInitB") {
        return;
    }

    let retval = CVodeSStolerancesB(&mut cvode_mem, indexB, reltolB, abstolB);
    if check_retval(retval, "CVodeSStolerancesB") {
        return;
    }

    /* Create banded SUNMatrix for the backward problem */
    let aB_mat = SUNBandMatrix(NEQ, MY, MY, &sunctx);

    /* Create banded SUNLinearSolver for the backward problem */
    let lsB = SUNLinSol_Band(&uB, &aB_mat, &sunctx);

    /* Attach the matrix and linear solver */
    let retval = CVodeSetLinearSolverB(&mut cvode_mem, indexB, lsB, Some(aB_mat));
    if check_retval(retval, "CVodeSetLinearSolverB") {
        return;
    }

    /* Set the user-supplied Jacobian routine for the backward problem */
    let retval = CVodeSetJacFnB(&mut cvode_mem, indexB, Some(jacB));
    if check_retval(retval, "CVodeSetJacFnB") {
        return;
    }

    /* Perform backward integration */
    println!("\nBackward integration");
    let retval = CVodeB(&mut cvode_mem, T0, CV_NORMAL);
    if check_retval(retval, "CVodeB") {
        return;
    }

    let retval = CVodeGetB(&mut cvode_mem, indexB, &mut t, &mut uB);
    if check_retval(retval, "CVodeGetB") {
        return;
    }

    PrintOutput(&uB, &data);

    /* Free memory */
    drop(u);
    drop(uB);
    CVodeFree(cvode_mem);
}
