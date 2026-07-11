/* -----------------------------------------------------------------
 * Translated from examples/idas/serial/idasAkzoNob_ASAi_dns.c
 * (IDAS 7.7.0)
 * Programmer(s): Radu Serban and Cosmin Petra @ LLNL
 *
 * Adjoint sensitivity example problem
 *
 * This IVP is a stiff system of 6 non-linear DAEs of index 1. The
 * problem originates from Akzo Nobel Central research in Arnhern,
 * The Netherlands, and describes a chemical process in which 2
 * species are mixed, while carbon dioxide is continuously added.
 * See http://pitagora.dm.uniba.it/~testset/report/chemakzo.pdf
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use idas_rs::sundials_utils::{fmt_e, fmt_f};
use idas_rs::*;

/* Problem Constants */
const NEQ: usize = 6;
const T0: f64 = 0.0;

const TF: f64 = 180.0; /* Final time. */

const RTOL: f64 = 1.0e-08;
const ATOL: f64 = 1.0e-10;
const RTOLB: f64 = 1.0e-06;
const ATOLB: f64 = 1.0e-08;
const RTOLQ: f64 = 1.0e-10;
const ATOLQ: f64 = 1.0e-12;

const ZERO: f64 = 0.0;
const HALF: f64 = 0.5;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;

const STEPS: i64 = 150;

#[derive(Clone)]
struct AkzoData {
    k1: f64,
    k2: f64,
    k3: f64,
    k4: f64,
    K: f64,
    klA: f64,
    Ks: f64,
    pCO2: f64,
    H: f64,
}

fn res(_t: f64, yy: &NVector, yd: &NVector, resval: &mut NVector, userdata: &mut UserData) -> i32 {
    let data = userdata.as_ref().unwrap().downcast_ref::<AkzoData>().unwrap();
    let k1 = data.k1;
    let k2 = data.k2;
    let k3 = data.k3;
    let k4 = data.k4;
    let K = data.K;
    let klA = data.klA;
    let Ks = data.Ks;
    let pCO2 = data.pCO2;
    let H = data.H;

    let y1 = yy.data[0];
    let y2 = yy.data[1];
    let y3 = yy.data[2];
    let y4 = yy.data[3];
    let y5 = yy.data[4];
    let y6 = yy.data[5];

    let yd1 = yd.data[0];
    let yd2 = yd.data[1];
    let yd3 = yd.data[2];
    let yd4 = yd.data[3];
    let yd5 = yd.data[4];

    let r1 = k1 * SUNRpowerI(y1, 4) * y2.sqrt();
    let r2 = k2 * y3 * y4;
    let r3 = k2 / K * y1 * y5;
    let r4 = k3 * y1 * y4 * y4;
    let r5 = k4 * y6 * y6 * y2.sqrt();
    let Fin = klA * (pCO2 / H - y2);

    resval.data[0] = yd1 + TWO * r1 - r2 + r3 + r4;
    resval.data[1] = yd2 + HALF * r1 + r4 + HALF * r5 - Fin;
    resval.data[2] = yd3 - r1 + r2 - r3;
    resval.data[3] = yd4 + r2 - r3 + TWO * r4;
    resval.data[4] = yd5 - r2 + r3 - r5;
    resval.data[5] = Ks * y1 * y4 - y6;

    0
}

/*
 * rhsQ routine. Computes quadrature(t,y).
 */
fn rhsQ(_t: f64, yy: &NVector, _yp: &NVector, qdot: &mut NVector, _user_data: &mut UserData) -> i32 {
    qdot.data[0] = yy.data[0];

    0
}

const QUARTER: f64 = 0.25;
const FOUR: f64 = 4.0;
const EIGHT: f64 = 8.0;

/*
 * resB routine. Residual for adjoint system.
 */
fn resB(_tt: f64, yy: &NVector, _yp: &NVector, yyB: &NVector, ypB: &NVector, rrB: &mut NVector,
        user_dataB: &mut UserData) -> i32 {
    let data = user_dataB.as_ref().unwrap().downcast_ref::<AkzoData>().unwrap();
    let k1 = data.k1;
    let k2 = data.k2;
    let k3 = data.k3;
    let k4 = data.k4;
    let K = data.K;
    let klA = data.klA;
    let Ks = data.Ks;

    let y1 = yy.data[0];
    let y2 = yy.data[1];
    let y3 = yy.data[2];
    let y4 = yy.data[3];
    let y5 = yy.data[4];
    let y6 = yy.data[5];

    let yB1 = yyB.data[0];
    let yB2 = yyB.data[1];
    let yB3 = yyB.data[2];
    let yB4 = yyB.data[3];
    let yB5 = yyB.data[4];
    let yB6 = yyB.data[5];

    let ypB1 = ypB.data[0];
    let ypB2 = ypB.data[1];
    let ypB3 = ypB.data[2];
    let ypB4 = ypB.data[3];
    let ypB5 = ypB.data[4];

    let y2tohalf = y2.sqrt();
    let y1to3 = y1 * y1 * y1;
    let k2overK = k2 / K;

    let tmp1 = k1 * y1to3 * y2tohalf;
    let tmp2 = k3 * y4 * y4;
    rrB.data[0] = 1.0 + ypB1 - (EIGHT * tmp1 + k2overK * y5 + tmp2) * yB1
        - (TWO * tmp1 + tmp2) * yB2
        + (FOUR * tmp1 + k2overK * y5) * yB3
        + k2overK * y5 * (yB4 - yB5)
        - TWO * tmp2 * yB4
        + Ks * y4 * yB6;

    let tmp1 = k1 * y1 * y1to3 * (y2tohalf / y2);
    let tmp2 = k4 * y6 * y6 * (y2tohalf / y2);
    rrB.data[1] = ypB2 - tmp1 * yB1 - (QUARTER * tmp1 + QUARTER * tmp2 + klA) * yB2
        + HALF * tmp1 * yB3
        + HALF * tmp2 * yB5;

    rrB.data[2] = ypB3 + k2 * y4 * (yB1 - yB3 - yB4 + yB5);

    let tmp1 = k3 * y1 * y4;
    let tmp2 = k2 * y3;
    rrB.data[3] = ypB4 + (tmp2 - TWO * tmp1) * yB1 - TWO * tmp1 * yB2 - tmp2 * yB3
        - (tmp2 + FOUR * tmp1) * yB4
        + tmp2 * yB5
        + Ks * y1 * yB6;

    rrB.data[4] = ypB5 - k2overK * y1 * (yB1 - yB3 - yB4 + yB5);

    rrB.data[5] = k4 * y6 * y2tohalf * (2.0 * yB5 - yB2) - yB6;

    0
}

/*
 * Print results after backward integration
 */
fn PrintOutput(_tfinal: f64, yB: &NVector, _ypB: &NVector) {
    print!(
        "dG/dy0: \t{}\n\t\t{}\n\t\t{}\n\t\t{}\n\t\t{}\n",
        fmt_e(yB.data[0], 12, 4),
        fmt_e(yB.data[1], 12, 4),
        fmt_e(yB.data[2], 12, 4),
        fmt_e(yB.data[3], 12, 4),
        fmt_e(yB.data[4], 12, 4)
    );
    print!("--------------------------------------------------------\n\n");
}

fn check_retval(retval: i32, funcname: &str) -> bool {
    if retval < 0 {
        eprintln!("\nSUNDIALS_ERROR: {}() failed with retval = {}\n", funcname, retval);
        return true;
    }
    false
}

/* Main program */
fn main() {
    /* Consistent IC for  y, y'. */
    let y01: f64 = 0.444;
    let y02: f64 = 0.00123;
    let y03: f64 = 0.0;
    let y04: f64 = 0.007;
    let y05: f64 = 0.0;

    print!("\nAdjoint Sensitivity Example for Akzo-Nobel Chemical Kinetics\n");
    print!("-------------------------------------------------------------\n");
    print!("Sensitivity of G = int_t0^tf (y1) dt with respect to IC.\n");
    print!("-------------------------------------------------------------\n\n");

    /* Create the SUNDIALS context object for this simulation */
    let sunctx = SUNContext_Create();

    /* Fill user's data with the appropriate values for coefficients. */
    let data = AkzoData {
        k1: 18.7,
        k2: 0.58,
        k3: 0.09,
        k4: 0.42,
        K: 34.4,
        klA: 3.3,
        Ks: 115.83,
        pCO2: 0.9,
        H: 737.0,
    };

    /* Allocate N-vectors. */
    let mut yy = N_VNew_Serial(NEQ as i64, &sunctx);
    let mut yp = N_VClone(&yy);

    /* Set IC */
    yy.data[0] = y01;
    yy.data[1] = y02;
    yy.data[2] = y03;
    yy.data[3] = y04;
    yy.data[4] = y05;
    yy.data[5] = data.Ks * y01 * y04;

    /* Get y' = - res(t0, y, 0) */
    N_VConst(ZERO, &mut yp);

    let mut ud: UserData = Some(Box::new(data.clone()));
    {
        let mut rr = N_VClone(&yy);
        res(T0, &yy, &yp, &mut rr, &mut ud);
        N_VScale(-ONE, &rr, &mut yp);
        /* (rr is dropped here — RAII) */
    }

    /* Create and initialize q0 for quadratures. */
    let mut q = N_VNew_Serial(1, &sunctx);
    q.data[0] = ZERO;

    /* Call IDACreate and IDAInit to initialize IDA memory */
    let mut mem = IDACreate(&sunctx);

    let retval = IDAInit(&mut mem, res, T0, &yy, &yp);
    if check_retval(retval, "IDAInit") {
        std::process::exit(1);
    }

    /* Set tolerances. */
    let retval = IDASStolerances(&mut mem, RTOL, ATOL);
    if check_retval(retval, "IDASStolerances") {
        std::process::exit(1);
    }

    /* Attach user data. */
    let retval = IDASetUserData(&mut mem, ud);
    if check_retval(retval, "IDASetUserData") {
        std::process::exit(1);
    }

    /* Create dense SUNMatrix for use in linear solves */
    let a = SUNDenseMatrix(NEQ as i64, NEQ as i64, &sunctx);

    /* Create dense SUNLinearSolver object */
    let ls = SUNLinSol_Dense(&yy, &a, &sunctx);

    /* Attach the matrix and linear solver */
    let retval = IDASetLinearSolver(&mut mem, ls, Some(a));
    if check_retval(retval, "IDASetLinearSolver") {
        std::process::exit(1);
    }

    /* Initialize QUADRATURE(S). */
    let retval = IDAQuadInit(&mut mem, rhsQ, &q);
    if check_retval(retval, "IDAQuadInit") {
        std::process::exit(1);
    }

    /* Set tolerances and error control for quadratures. */
    let retval = IDAQuadSStolerances(&mut mem, RTOLQ, ATOLQ);
    if check_retval(retval, "IDAQuadSStolerances") {
        std::process::exit(1);
    }

    let retval = IDASetQuadErrCon(&mut mem, true);
    if check_retval(retval, "IDASetQuadErrCon") {
        std::process::exit(1);
    }

    /* Prepare ADJOINT. */
    let retval = IDAAdjInit(&mut mem, STEPS, IDA_HERMITE);
    if check_retval(retval, "IDAAdjInit") {
        std::process::exit(1);
    }

    /* FORWARD run. */
    print!("Forward integration ... ");
    let mut time = 0.0;
    let mut ncheck = 0;
    let retval = IDASolveF(&mut mem, TF, &mut time, &mut yy, &mut yp, IDA_NORMAL, &mut ncheck);
    if check_retval(retval, "IDASolveF") {
        std::process::exit(1);
    }

    let mut nst = 0i64;
    let retval = IDAGetNumSteps(&mut mem, &mut nst);
    if check_retval(retval, "IDAGetNumSteps") {
        std::process::exit(1);
    }

    println!("done ( nst = {} )", nst);

    let retval = IDAGetQuad(&mem, &mut time, &mut q);
    if check_retval(retval, "IDAGetQuad") {
        std::process::exit(1);
    }

    println!("G:          {} ", fmt_f(q.data[0], 24, 16));
    print!("--------------------------------------------------------\n\n");

    /* BACKWARD run */

    /* Initialize yB */
    let mut yB = N_VClone(&yy);
    N_VConst(ZERO, &mut yB);

    let mut ypB = N_VClone(&yB);
    N_VConst(ZERO, &mut ypB);
    ypB.data[0] = -ONE;

    let mut indexB = 0;
    let retval = IDACreateB(&mut mem, &mut indexB);
    if check_retval(retval, "IDACreateB") {
        std::process::exit(1);
    }

    let retval = IDAInitB(&mut mem, indexB, resB, TF, &yB, &ypB);
    if check_retval(retval, "IDAInitB") {
        std::process::exit(1);
    }

    let retval = IDASStolerancesB(&mut mem, indexB, RTOLB, ATOLB);
    if check_retval(retval, "IDASStolerancesB") {
        std::process::exit(1);
    }

    let retval = IDASetUserDataB(&mut mem, indexB, Some(Box::new(data.clone())));
    if check_retval(retval, "IDASetUserDataB") {
        std::process::exit(1);
    }

    IDASetMaxNumStepsB(&mut mem, indexB, 1000);

    /* Create dense SUNMatrix for use in linear solves */
    let ab = SUNDenseMatrix(NEQ as i64, NEQ as i64, &sunctx);

    /* Create dense SUNLinearSolver object */
    let lsb = SUNLinSol_Dense(&yB, &ab, &sunctx);

    /* Attach the matrix and linear solver */
    let retval = IDASetLinearSolverB(&mut mem, indexB, lsb, Some(ab));
    if check_retval(retval, "IDASetLinearSolverB") {
        std::process::exit(1);
    }

    print!("Backward integration ... ");

    let retval = IDASolveB(&mut mem, T0, IDA_NORMAL);
    if check_retval(retval, "IDASolveB") {
        std::process::exit(1);
    }

    let mut nstB = 0i64;
    IDAGetNumSteps(IDAGetAdjIDABmem(&mut mem, indexB).unwrap(), &mut nstB);
    println!("done ( nst = {} )", nstB);

    let retval = IDAGetB(&mut mem, indexB, &mut time, &mut yB, &mut ypB);
    if check_retval(retval, "IDAGetB") {
        std::process::exit(1);
    }

    PrintOutput(time, &yB, &ypB);

    /* Free memory (RAII) */
    IDAFree(mem);
}
