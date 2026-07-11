/* -----------------------------------------------------------------
 * Translated from examples/idas/serial/idasAkzoNob_dns.c (IDAS 7.7.0)
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

use idas_rs::sundials_utils::{fmt_e, fmt_f, fmt_g};
use idas_rs::*;

/* Problem Constants */
const NEQ: usize = 6;
const T0: f64 = 0.0;
const T1: f64 = 1e-8; /* first time for output */

const TF: f64 = 180.0; /* Final time. */
const NF: i32 = 25; /* Total number of outputs. */

const RTOL: f64 = 1.0e-08;
const ATOL: f64 = 1.0e-10;
const RTOLQ: f64 = 1.0e-10;
const ATOLQ: f64 = 1.0e-12;

const ZERO: f64 = 0.0;
const HALF: f64 = 0.5;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;

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

    let r1 = k1 * SUNRpowerI(y1, 4) * SUNRsqrt(y2);
    let r2 = k2 * y3 * y4;
    let r3 = k2 / K * y1 * y5;
    let r4 = k3 * y1 * y4 * y4;
    let r5 = k4 * y6 * y6 * SUNRsqrt(y2);
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

fn PrintHeader(rtol: f64, avtol: f64) {
    print!("\nidasAkzoNob_dns: Akzo Nobel chemical kinetics DAE serial example problem for IDAS\n");
    print!("Linear solver: DENSE, Jacobian is computed by IDAS.\n");
    print!(
        "Tolerance parameters:  rtol = {}   atol = {}\n",
        fmt_g(rtol, 0, 6),
        fmt_g(avtol, 0, 6)
    );
    print!("---------------------------------------------------------------------------------\n");
    print!("   t        y1        y2       y3       y4       y5");
    print!("      y6    | nst  k      h\n");
    print!("---------------------------------------------------------------------------------\n");
}

fn PrintOutput(mem: &mut IDAMem, t: f64, y: &NVector) {
    let mut kused = 0i32;
    let mut nst = 0i64;
    let mut hused = 0.0f64;

    IDAGetLastOrder(mem, &mut kused);
    IDAGetNumSteps(mem, &mut nst);
    IDAGetLastStep(mem, &mut hused);

    println!(
        "{} {} {} {} {} {} {} | {:3}  {:1} {}",
        fmt_e(t, 8, 2),
        fmt_e(y.data[0], 8, 2),
        fmt_e(y.data[1], 8, 2),
        fmt_e(y.data[2], 8, 2),
        fmt_e(y.data[3], 8, 2),
        fmt_e(y.data[4], 8, 2),
        fmt_e(y.data[5], 8, 2),
        nst,
        kused,
        fmt_e(hused, 8, 2)
    );
}

fn PrintFinalStats(mem: &mut IDAMem) {
    let (mut nst, mut nni, mut nje, mut nre, mut nreLS, mut netf, mut ncfn) =
        (0i64, 0i64, 0i64, 0i64, 0i64, 0i64, 0i64);

    IDAGetNumSteps(mem, &mut nst);
    IDAGetNumResEvals(mem, &mut nre);
    IDAGetNumJacEvals(mem, &mut nje);
    IDAGetNumNonlinSolvIters(mem, &mut nni);
    IDAGetNumErrTestFails(mem, &mut netf);
    IDAGetNumNonlinSolvConvFails(mem, &mut ncfn);
    IDAGetNumLinResEvals(mem, &mut nreLS);

    print!("\nFinal Run Statistics: \n\n");
    println!("Number of steps                    = {}", nst);
    println!("Number of residual evaluations     = {}", nre + nreLS);
    println!("Number of Jacobian evaluations     = {}", nje);
    println!("Number of nonlinear iterations     = {}", nni);
    println!("Number of error test failures      = {}", netf);
    println!("Number of nonlinear conv. failures = {}", ncfn);
}

/* Main program */
fn main() {
    /* Consistent IC for  y, y'. */
    let y01: f64 = 0.444;
    let y02: f64 = 0.00123;
    let y03: f64 = 0.0;
    let y04: f64 = 0.007;
    let y05: f64 = 0.0;

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

    let mut ud: UserData = Some(Box::new(data));
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

    PrintHeader(RTOL, ATOL);
    /* Print initial states */
    PrintOutput(&mut mem, 0.0, &yy);

    let mut tout = T1;
    let mut nout = 0;
    let incr = SUNRpowerR(TF / T1, ONE / NF as f64);

    /* FORWARD run. */
    let mut time = 0.0;
    loop {
        let retval = IDASolve(&mut mem, tout, &mut time, &mut yy, &mut yp, IDA_NORMAL);
        if check_retval(retval, "IDASolve") {
            std::process::exit(1);
        }

        PrintOutput(&mut mem, time, &yy);

        nout += 1;
        tout *= incr;

        if nout > NF {
            break;
        }
    }

    let retval = IDAGetQuad(&mem, &mut time, &mut q);
    if check_retval(retval, "IDAGetQuad") {
        std::process::exit(1);
    }

    print!("\n--------------------------------------------------------\n");
    println!("G:          {} ", fmt_f(q.data[0], 24, 16));
    print!("--------------------------------------------------------\n\n");

    PrintFinalStats(&mut mem);

    /* Free memory (RAII) */
    IDAFree(mem);
}

fn check_retval(retval: i32, funcname: &str) -> bool {
    if retval < 0 {
        eprintln!("\nSUNDIALS_ERROR: {}() failed with retval = {}\n", funcname, retval);
        return true;
    }
    false
}
