/* -----------------------------------------------------------------
 * Translated from examples/idas/serial/idasHeat2D_bnd.c (IDAS 7.7.0)
 * Programmer(s): Allan Taylor, Alan Hindmarsh and Radu Serban @ LLNL
 *
 * Example problem for IDA: 2D heat equation, serial, banded.
 *
 * The DAE system is a spatial discretization of
 *          du/dt = d^2u/dx^2 + d^2u/dy^2
 * on the unit square, u = 0 on the boundary, initial condition
 * u = 16 x (1-x) y (1-y). Central differences on a uniform M x M grid
 * (M = 10); interior points give ODEs, boundary points give algebraic
 * equations u = 0. Solved with the band solver, half-bandwidths M,
 * default DQ Jacobian; IDACalcIC corrects boundary values; constraints
 * u >= 0. Output at t = 0, .01, .02, .04, ..., 10.24.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use idas_rs::sundials_utils::{fmt_e, fmt_f, fmt_g};
use idas_rs::*;

/* Problem Constants */
const NOUT: i32 = 11;
const MGRID: usize = 10;
const NEQ: usize = MGRID * MGRID;
const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;
const BVAL: f64 = 0.1;

/* Type: UserData */
struct HeatData {
    mm: usize,
    #[allow(dead_code)]
    dx: f64,
    coeff: f64,
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY IDA
 *--------------------------------------------------------------------
 */

/* heatres core: 5-point central differencing on interior points,
   res = u for boundary points (via the initial N_VScale). */
fn heatres_compute(mm: usize, coeff: f64, uu: &NVector, up: &NVector, resval: &mut NVector) {
    /* Initialize resval to uu, to take care of boundary equations. */
    N_VScale(ONE, uu, resval);

    /* Loop over interior points; set res = up - (central difference). */
    for j in 1..mm - 1 {
        let offset = mm * j;
        for i in 1..mm - 1 {
            let loc = offset + i;
            resval.data[loc] = up.data[loc]
                - coeff
                    * (uu.data[loc - 1]
                        + uu.data[loc + 1]
                        + uu.data[loc - mm]
                        + uu.data[loc + mm]
                        - 4.0 * uu.data[loc]);
        }
    }
}

fn heatres(_tres: f64, uu: &NVector, up: &NVector, resval: &mut NVector, user_data: &mut UserData) -> i32 {
    let data = user_data.as_ref().unwrap().downcast_ref::<HeatData>().unwrap();
    heatres_compute(data.mm, data.coeff, uu, up, resval);
    0
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/* SetInitialProfile: initialize u, up, and id vectors. */
fn SetInitialProfile(data: &HeatData, uu: &mut NVector, up: &mut NVector, id: &mut NVector, res: &mut NVector) {
    let mm = data.mm;
    let mm1 = mm - 1;

    /* Initialize id to 1's. */
    N_VConst(ONE, id);

    /* Initialize uu on all grid points. */
    for j in 0..mm {
        let yfact = data.dx * j as f64;
        let offset = mm * j;
        for i in 0..mm {
            let xfact = data.dx * i as f64;
            let loc = offset + i;
            uu.data[loc] = 16.0 * xfact * (ONE - xfact) * yfact * (ONE - yfact);
        }
    }

    /* Initialize up vector to 0. */
    N_VConst(ZERO, up);

    /* heatres sets res to negative of ODE RHS values at interior points. */
    heatres_compute(data.mm, data.coeff, uu, up, res);

    /* Copy -res into up to get correct interior initial up values. */
    N_VScale(-ONE, res, up);

    /* Finally, set values of u, up, and id at boundary points. */
    for j in 0..mm {
        let offset = mm * j;
        for i in 0..mm {
            let loc = offset + i;
            if j == 0 || j == mm1 || i == 0 || i == mm1 {
                uu.data[loc] = BVAL;
                up.data[loc] = ZERO;
                id.data[loc] = ZERO;
            }
        }
    }
}

/* Print first lines of output (problem description). */
fn PrintHeader(rtol: f64, atol: f64) {
    print!("\nidasHeat2D_bnd: Heat equation, serial example problem for IDA\n");
    print!("              Discretized heat equation on 2D unit square.\n");
    print!("              Zero boundary conditions,");
    print!(" polynomial initial conditions.\n");
    print!("              Mesh dimensions: {} x {}", MGRID, MGRID);
    print!("        Total system size: {}\n\n", NEQ);
    print!(
        "Tolerance parameters:  rtol = {}   atol = {}\n",
        fmt_g(rtol, 0, 6),
        fmt_g(atol, 0, 6)
    );
    print!("Constraints set to force all solution components >= 0. \n");
    print!("Linear solver: BAND, banded direct solver \n");
    print!("       difference quotient Jacobian, half-bandwidths = {} \n", MGRID);
    print!("IDACalcIC called with input boundary values = {} \n", fmt_g(BVAL, 0, 6));
    /* Print output table heading and initial line of table. */
    print!("\n   Output Summary (umax = max-norm of solution) \n\n");
    print!("  time       umax     k  nst  nni  nje   nre   nreLS    h      \n");
    print!(" .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . \n");
}

/* Print Output */
fn PrintOutput(mem: &mut IDAMem, t: f64, uu: &NVector) {
    let umax = N_VMaxNorm(uu);

    let mut kused = 0i32;
    let mut nst = 0i64;
    let mut nni = 0i64;
    let mut nre = 0i64;
    let mut hused = 0.0f64;
    let mut nje = 0i64;
    let mut nreLS = 0i64;

    IDAGetLastOrder(mem, &mut kused);
    IDAGetNumSteps(mem, &mut nst);
    IDAGetNumNonlinSolvIters(mem, &mut nni);
    IDAGetNumResEvals(mem, &mut nre);
    IDAGetLastStep(mem, &mut hused);
    IDAGetNumJacEvals(mem, &mut nje);
    IDAGetNumLinResEvals(mem, &mut nreLS);

    println!(
        " {} {}  {}  {:3}  {:3}  {:3}  {:4}  {:4}  {} ",
        fmt_f(t, 5, 2),
        fmt_e(umax, 13, 5),
        kused,
        nst,
        nni,
        nje,
        nre,
        nreLS,
        fmt_e(hused, 9, 2)
    );
}

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
    /* Create the SUNDIALS context object for this simulation */
    let sunctx = SUNContext_Create();

    /* Create vectors uu, up, res, constraints, id. */
    let mut uu = N_VNew_Serial(NEQ as i64, &sunctx);
    let mut up = N_VClone(&uu);
    let mut res = N_VClone(&uu);
    let mut constraints = N_VClone(&uu);
    let mut id = N_VClone(&uu);

    /* Create and load problem data block. */
    let data = HeatData {
        mm: MGRID,
        dx: ONE / (MGRID as f64 - ONE),
        coeff: ONE / ((ONE / (MGRID as f64 - ONE)) * (ONE / (MGRID as f64 - ONE))),
    };

    /* Initialize uu, up, id. */
    SetInitialProfile(&data, &mut uu, &mut up, &mut id, &mut res);

    /* Set constraints to all 1's for nonnegative solution values. */
    N_VConst(ONE, &mut constraints);

    /* Set remaining input parameters. */
    let t0 = ZERO;
    let t1 = 0.01;
    let rtol = ZERO;
    let atol = 1.0e-3;

    /* Call IDACreate and IDAMalloc to initialize solution */
    let mut mem = IDACreate(&sunctx);

    let mut retval = IDASetUserData(&mut mem, Some(Box::new(data)));
    if check_retval(retval, "IDASetUserData") {
        std::process::exit(1);
    }

    /* Set which components are algebraic or differential */
    retval = IDASetId(&mut mem, Some(&id));
    if check_retval(retval, "IDASetId") {
        std::process::exit(1);
    }

    retval = IDASetConstraints(&mut mem, Some(&constraints));
    if check_retval(retval, "IDASetConstraints") {
        std::process::exit(1);
    }
    drop(constraints);

    retval = IDAInit(&mut mem, heatres, t0, &uu, &up);
    if check_retval(retval, "IDAInit") {
        std::process::exit(1);
    }

    retval = IDASStolerances(&mut mem, rtol, atol);
    if check_retval(retval, "IDASStolerances") {
        std::process::exit(1);
    }

    /* Create banded SUNMatrix for use in linear solves */
    let mu = MGRID as i64;
    let ml = MGRID as i64;
    let a = SUNBandMatrix(NEQ as i64, mu, ml, &sunctx);

    /* Create banded SUNLinearSolver object */
    let ls = SUNLinSol_Band(&uu, &a, &sunctx);

    /* Attach the matrix and linear solver */
    retval = IDASetLinearSolver(&mut mem, ls, Some(a));
    if check_retval(retval, "IDASetLinearSolver") {
        std::process::exit(1);
    }

    /* Call IDACalcIC to correct the initial values. */
    retval = IDACalcIC(&mut mem, IDA_YA_YDP_INIT, t1);
    if check_retval(retval, "IDACalcIC") {
        std::process::exit(1);
    }

    /* Print output heading. */
    PrintHeader(rtol, atol);

    PrintOutput(&mut mem, t0, &uu);

    /* Loop over output times, call IDASolve, and print results. */
    let mut tret = 0.0;
    let mut tout = t1;
    for _iout in 1..=NOUT {
        retval = IDASolve(&mut mem, tout, &mut tret, &mut uu, &mut up, IDA_NORMAL);
        if check_retval(retval, "IDASolve") {
            std::process::exit(1);
        }

        PrintOutput(&mut mem, tret, &uu);
        tout *= TWO;
    }

    /* Print remaining counters. */
    let mut netf = 0i64;
    let mut ncfn = 0i64;
    IDAGetNumErrTestFails(&mut mem, &mut netf);
    IDAGetNumNonlinSolvConvFails(&mut mem, &mut ncfn);
    println!("\n netf = {},   ncfn = {} ", netf, ncfn);

    /* Free memory (RAII) */
    IDAFree(mem);
}
