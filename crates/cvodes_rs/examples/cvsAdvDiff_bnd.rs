/* -----------------------------------------------------------------
 * Translated from examples/cvode/serial/cvsAdvDiff_bnd.c (CVODE 7.7.0)
 *
 * Example problem:
 *
 * The following is a simple example problem with a banded Jacobian,
 * with the program for its solution by CVODE.
 * The problem is the semi-discrete form of the advection-diffusion
 * equation in 2-D:
 *   du/dt = d^2 u / dx^2 + .5 du/dx + d^2 u / dy^2
 * on the rectangle 0 <= x <= 2, 0 <= y <= 1, and the time
 * interval 0 <= t <= 1. Homogeneous Dirichlet boundary conditions
 * are posed, and the initial condition is
 *   u(x,y,t=0) = x(2-x)y(1-y)exp(5xy).
 * The PDE is discretized on a uniform MX+2 by MY+2 grid with
 * central differencing, and with boundary values eliminated,
 * leaving an ODE system of size NEQ = MX*MY.
 * This program solves the problem with the BDF method, Newton
 * iteration with the SUNBAND linear solver, and a user-supplied
 * Jacobian routine.
 * It uses scalar relative and absolute tolerances.
 * Output is printed at t = .1, .2, ..., 1.
 * Run statistics (optional outputs) are printed at the end.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use cvodes_rs::sundials_utils::{fmt_e, fmt_f, fmt_g};
use cvodes_rs::*;

/* Problem Constants */

const XMAX: f64 = 2.0; /* domain boundaries         */
const YMAX: f64 = 1.0;
const MX: i64 = 10; /* mesh dimensions           */
const MY: i64 = 5;
const NEQ: i64 = MX * MY; /* number of equations       */
const ATOL: f64 = 1.0e-5; /* scalar absolute tolerance */
const T0: f64 = 0.0; /* initial time              */
const T1: f64 = 0.1; /* first output time         */
const DTOUT: f64 = 0.1; /* output time increment     */
const NOUT: i32 = 10; /* number of output times    */

const ZERO: f64 = 0.0;
const HALF: f64 = 0.5;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;
const FIVE: f64 = 5.0;

/* IJth references the element in the vdata array for
   u at mesh point (i,j), where 1 <= i <= MX, 1 <= j <= MY.
   The variables are ordered by the y index j, then by the x index i. */
macro_rules! IJth {
    ($vdata:expr, $i:expr, $j:expr) => {
        $vdata[(($j - 1) + ($i - 1) * MY) as usize]
    };
}

/* Type : UserDataStruct (contains grid constants) */
struct UserDataStruct {
    dx: f64,
    #[allow(dead_code)]
    dy: f64,
    hdcoef: f64,
    hacoef: f64,
    vdcoef: f64,
}

/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/* f routine. Compute f(t,u). */

fn f(_t: f64, u: &NVector, udot: &mut NVector, user_data: &mut UserData) -> i32 {
    let data = user_data
        .as_mut()
        .unwrap()
        .downcast_mut::<UserDataStruct>()
        .unwrap();

    let udata = &u.data;
    let dudata = &mut udot.data;

    /* Extract needed constants from data */

    let hordc = data.hdcoef;
    let horac = data.hacoef;
    let verdc = data.vdcoef;

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

/* Jacobian routine. Compute J(t,u). */

fn Jac(
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
     * The components of f = udot that depend on u(i,j) are
     * f(i,j), f(i-1,j), f(i+1,j), f(i,j-1), f(i,j+1), with
     *   df(i,j)/du(i,j) = -2 (1/dx^2 + 1/dy^2)
     *   df(i-1,j)/du(i,j) = 1/dx^2 + .25/dx  (if i > 1)
     *   df(i+1,j)/du(i,j) = 1/dx^2 - .25/dx  (if i < MX)
     *   df(i,j-1)/du(i,j) = 1/dy^2           (if j > 1)
     *   df(i,j+1)/du(i,j) = 1/dy^2           (if j < MY)
     */

    let data = user_data
        .as_mut()
        .unwrap()
        .downcast_mut::<UserDataStruct>()
        .unwrap();
    let hordc = data.hdcoef;
    let horac = data.hacoef;
    let verdc = data.vdcoef;

    let jm = match j_mat {
        SUNMatrix::Band(m) => m,
        _ => return -1,
    };

    /* set non-zero Jacobian entries */
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
 *-------------------------------
 * Private helper functions
 *-------------------------------
 */

/* Set initial conditions in u vector */

fn SetIC(u: &mut NVector, data: &UserDataStruct) {
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
            IJth!(udata, i, j) = x * (XMAX - x) * y * (YMAX - y) * (FIVE * x * y).exp();
        }
    }
}

/* Print first lines of output (problem description) */

fn PrintHeader(reltol: f64, abstol: f64, umax: f64) {
    println!("\n2-D Advection-Diffusion Equation");
    println!("Mesh dimensions = {} X {}", MX, MY);
    println!("Total system size = {}", NEQ);
    println!(
        "Tolerance parameters: reltol = {}   abstol = {}\n",
        fmt_g(reltol, 0, 6),
        fmt_g(abstol, 0, 6)
    );
    println!("At t = {}      max.norm(u) ={} ", fmt_g(T0, 0, 6), fmt_e(umax, 14, 6));
}

/* Print current value */

fn PrintOutput(t: f64, umax: f64, nst: i64) {
    println!(
        "At t = {}   max.norm(u) ={}   nst = {:4}",
        fmt_f(t, 4, 2),
        fmt_e(umax, 14, 6),
        nst
    );
}

/* Get and print some final statistics */

fn PrintFinalStats(cvode_mem: &mut CVodeMem) {
    let (mut nst, mut nfe, mut nsetups, mut netf, mut nni, mut ncfn) =
        (0i64, 0i64, 0i64, 0i64, 0i64, 0i64);
    let (mut nje, mut nfeLS) = (0i64, 0i64);

    let mut retval = CVodeGetNumSteps(cvode_mem, &mut nst);
    check_retval(retval, "CVodeGetNumSteps");
    retval = CVodeGetNumRhsEvals(cvode_mem, &mut nfe);
    check_retval(retval, "CVodeGetNumRhsEvals");
    retval = CVodeGetNumLinSolvSetups(cvode_mem, &mut nsetups);
    check_retval(retval, "CVodeGetNumLinSolvSetups");
    retval = CVodeGetNumErrTestFails(cvode_mem, &mut netf);
    check_retval(retval, "CVodeGetNumErrTestFails");
    retval = CVodeGetNumNonlinSolvIters(cvode_mem, &mut nni);
    check_retval(retval, "CVodeGetNumNonlinSolvIters");
    retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &mut ncfn);
    check_retval(retval, "CVodeGetNumNonlinSolvConvFails");

    retval = CVodeGetNumJacEvals(cvode_mem, &mut nje);
    check_retval(retval, "CVodeGetNumJacEvals");
    retval = CVodeGetNumLinRhsEvals(cvode_mem, &mut nfeLS);
    check_retval(retval, "CVodeGetNumLinRhsEvals");

    println!("\nFinal Statistics:");
    println!(
        "nst = {:<6} nfe  = {:<6} nsetups = {:<6} nfeLS = {:<6} nje = {}",
        nst, nfe, nsetups, nfeLS, nje
    );
    print!("nni = {:<6} ncfn = {:<6} netf = {}\n \n", nni, ncfn, netf);
}

/* Check function return value */
fn check_retval(retval: i32, funcname: &str) -> bool {
    if retval < 0 {
        eprintln!("\nSUNDIALS_ERROR: {}() failed with retval = {}\n", funcname, retval);
        return true;
    }
    false
}

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

fn main() {
    /* Create the SUNDIALS context */
    let sunctx = SUNContext_Create();

    /* Create a serial vector */

    let mut u = N_VNew_Serial(NEQ, &sunctx); /* Allocate u vector */

    let reltol = ZERO; /* Set the tolerances */
    let abstol = ATOL;

    /* Set grid coefficients in data */
    let dx = XMAX / (MX + 1) as f64;
    let dy = YMAX / (MY + 1) as f64;
    let data = UserDataStruct {
        dx,
        dy,
        hdcoef: ONE / (dx * dx),
        hacoef: HALF / (TWO * dx),
        vdcoef: ONE / (dy * dy),
    };

    SetIC(&mut u, &data); /* Initialize u vector */

    /* Call CVodeCreate to create the solver memory and specify the
     * Backward Differentiation Formula */

    let mut cvode_mem = CVodeCreate(CV_BDF, &sunctx);

    /* Call CVodeInit to initialize the integrator memory and specify the
     * user's right hand side function in u'=f(t,u), the initial time T0, and
     * the initial dependent variable vector u. */
    let mut retval = CVodeInit(&mut cvode_mem, f, T0, &u);
    if check_retval(retval, "CVodeInit") {
        std::process::exit(1);
    }

    /* Call CVodeSStolerances to specify the scalar relative tolerance
     * and scalar absolute tolerance */
    retval = CVodeSStolerances(&mut cvode_mem, reltol, abstol);
    if check_retval(retval, "CVodeSStolerances") {
        std::process::exit(1);
    }

    /* Set the pointer to user-defined data */
    retval = CVodeSetUserData(&mut cvode_mem, Some(Box::new(data)));
    if check_retval(retval, "CVodeSetUserData") {
        std::process::exit(1);
    }

    /* Create banded SUNMatrix for use in linear solves -- since this will be factored,
       set the storage bandwidth to be the sum of upper and lower bandwidths */
    let a = SUNBandMatrix(NEQ, MY, MY, &sunctx);

    /* Create banded SUNLinearSolver object for use by CVode */
    let ls = SUNLinSol_Band(&u, &a, &sunctx);

    /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
    retval = CVodeSetLinearSolver(&mut cvode_mem, ls, Some(a));
    if check_retval(retval, "CVodeSetLinearSolver") {
        std::process::exit(1);
    }

    /* Set the user-supplied Jacobian routine Jac */
    retval = CVodeSetJacFn(&mut cvode_mem, Some(Jac));
    if check_retval(retval, "CVodeSetJacFn") {
        std::process::exit(1);
    }

    /* In loop over output points: call CVode, print results, test for errors */

    let mut umax = N_VMaxNorm(&u);
    PrintHeader(reltol, abstol, umax);
    let mut t = 0.0;
    let mut tout = T1;
    for _iout in 1..=NOUT {
        retval = CVode(&mut cvode_mem, tout, &mut u, &mut t, CV_NORMAL);
        if check_retval(retval, "CVode") {
            break;
        }
        umax = N_VMaxNorm(&u);
        let mut nst = 0i64;
        retval = CVodeGetNumSteps(&mut cvode_mem, &mut nst);
        check_retval(retval, "CVodeGetNumSteps");
        PrintOutput(t, umax, nst);
        tout += DTOUT;
    }
    PrintFinalStats(&mut cvode_mem); /* Print some final statistics   */

    CVodeFree(cvode_mem); /* Free the integrator memory */
}
