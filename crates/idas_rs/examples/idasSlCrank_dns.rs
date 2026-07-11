/* -----------------------------------------------------------------
 * Translated from examples/idas/serial/idasSlCrank_dns.c (IDAS 7.7.0)
 * Programmer: Radu Serban and Cosmin Petra @ LLNL
 *
 * Simulation of a slider-crank mechanism modelled with 3 generalized
 * coordinates: crank angle, connecting bar angle, and slider location.
 * The equations of motion are formulated as stabilized index-2 DAEs
 * (Gear-Gupta-Leimkuhler). Dense linear solver, IDAS-computed Jacobian.
 *
 * IDAS also computes the average kinetic energy as the quadrature:
 *   G = int_t0^tend g(t,y,p) dt,
 * where
 *   g(t,y,p) = 0.5*J1*v1^2 + 0.5*J2*v3^2 + 0.5*m2*v2^2
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use idas_rs::sundials_utils::{fmt_e, fmt_f, fmt_g};
use idas_rs::*;

/* Problem Constants */
const NEQ: usize = 10;

const TBEGIN: f64 = 0.0;
const TEND: f64 = 10.0;

const NOUT: i32 = 25;

const RTOLF: f64 = 1.0e-06;
const ATOLF: f64 = 1.0e-07;

const RTOLQ: f64 = 1.0e-06;
const ATOLQ: f64 = 1.0e-08;

const ZERO: f64 = 0.00;
const HALF: f64 = 0.50;
const ONE: f64 = 1.00;
const TWO: f64 = 2.00;
const FOUR: f64 = 4.00;

struct SlData {
    a: f64,
    J1: f64,
    J2: f64,
    m1: f64,
    m2: f64,
    l0: f64,
    params: [f64; 2],
    F: f64,
}

fn force(yy: &NVector, data: &SlData) -> [f64; 3] {
    let a = data.a;
    let k = data.params[0];
    let c = data.params[1];
    let l0 = data.l0;
    let F = data.F;

    let q = yy.data[0];
    let x = yy.data[1];
    let p = yy.data[2];

    let qd = yy.data[3];
    let xd = yy.data[4];
    let pd = yy.data[5];

    let s1 = q.sin();
    let c1 = q.cos();
    let s2 = p.sin();
    let c2 = p.cos();
    let s21 = s2 * c1 - c2 * s1;
    let c21 = c2 * c1 + s2 * s1;

    let l2 = x * x - x * (c2 + a * c1) + (ONE + a * a) / FOUR + a * c21 / TWO;
    let l = l2.sqrt();
    let mut ld = TWO * x * xd - xd * (c2 + a * c1) + x * (s2 * pd + a * s1 * qd)
        - a * s21 * (pd - qd) / TWO;
    ld /= TWO * l;

    let f = k * (l - l0) + c * ld;
    let fl = f / l;

    let mut Q = [0.0f64; 3];
    Q[0] = -fl * a * (s21 / TWO + x * s1) / TWO;
    Q[1] = fl * (c2 / TWO - x + a * c1 / TWO) + F;
    Q[2] = -fl * (x * s2 - a * s21 / TWO) / TWO - F * s2;
    Q
}

fn setIC(yy: &mut NVector, yp: &mut NVector, data: &SlData) {
    N_VConst(ZERO, yy);
    N_VConst(ZERO, yp);

    let pi = FOUR * ONE.atan();

    let a = data.a;
    let J1 = data.J1;
    let m2 = data.m2;
    let J2 = data.J2;

    let q = pi / TWO;
    let p = (-a).asin();
    let x = p.cos();

    yy.data[0] = q;
    yy.data[1] = x;
    yy.data[2] = p;

    let Q = force(yy, data);

    yp.data[3] = Q[0] / J1;
    yp.data[4] = Q[1] / m2;
    yp.data[5] = Q[2] / J2;
}

fn ressc(_tres: f64, yy: &NVector, yp: &NVector, rr: &mut NVector, user_data: &mut UserData) -> i32 {
    let data = user_data.as_ref().unwrap().downcast_ref::<SlData>().unwrap();

    let a = data.a;
    let J1 = data.J1;
    let m2 = data.m2;
    let J2 = data.J2;

    let q = yy.data[0];
    let p = yy.data[2];

    let qd = yy.data[3];
    let xd = yy.data[4];
    let pd = yy.data[5];

    let lam1 = yy.data[6];
    let lam2 = yy.data[7];

    let mu1 = yy.data[8];
    let mu2 = yy.data[9];

    let s1 = q.sin();
    let c1 = q.cos();
    let s2 = p.sin();
    let c2 = p.cos();

    let Q = force(yy, data);

    let x = yy.data[1];

    rr.data[0] = yp.data[0] - qd + a * s1 * mu1 - a * c1 * mu2;
    rr.data[1] = yp.data[1] - xd + mu1;
    rr.data[2] = yp.data[2] - pd + s2 * mu1 - c2 * mu2;

    rr.data[3] = J1 * yp.data[3] - Q[0] + a * s1 * lam1 - a * c1 * lam2;
    rr.data[4] = m2 * yp.data[4] - Q[1] + lam1;
    rr.data[5] = J2 * yp.data[5] - Q[2] + s2 * lam1 - c2 * lam2;

    rr.data[6] = x - c2 - a * c1;
    rr.data[7] = -s2 - a * s1;

    rr.data[8] = a * s1 * qd + xd + s2 * pd;
    rr.data[9] = -a * c1 * qd - c2 * pd;

    0
}

fn rhsQ(_t: f64, yy: &NVector, _yp: &NVector, qdot: &mut NVector, user_data: &mut UserData) -> i32 {
    let data = user_data.as_ref().unwrap().downcast_ref::<SlData>().unwrap();
    let J1 = data.J1;
    let m2 = data.m2;
    let J2 = data.J2;

    let v1 = yy.data[3];
    let v2 = yy.data[4];
    let v3 = yy.data[5];

    qdot.data[0] = HALF * (J1 * v1 * v1 + m2 * v2 * v2 + J2 * v3 * v3);

    0
}

fn PrintHeader(rtol: f64, avtol: f64) {
    print!("\nidasSlCrank_dns: Slider-Crank DAE serial example problem for IDAS\n");
    print!("Linear solver: DENSE, Jacobian is computed by IDAS.\n");
    print!(
        "Tolerance parameters:  rtol = {}   atol = {}\n",
        fmt_g(rtol, 0, 6),
        fmt_g(avtol, 0, 6)
    );
    print!("-----------------------------------------------------------------------\n");
    print!("  t         y1          y2           y3");
    print!("      | nst  k      h\n");
    print!("-----------------------------------------------------------------------\n");
}

fn PrintOutput(mem: &mut IDAMem, t: f64, y: &NVector) {
    let mut kused = 0i32;
    let mut nst = 0i64;
    let mut hused = 0.0f64;

    IDAGetLastOrder(mem, &mut kused);
    IDAGetNumSteps(mem, &mut nst);
    IDAGetLastStep(mem, &mut hused);

    println!(
        "{} {} {} {} | {:3}  {:1} {}",
        fmt_f(t, 5, 2),
        fmt_e(y.data[0], 12, 4),
        fmt_e(y.data[1], 12, 4),
        fmt_e(y.data[2], 12, 4),
        nst,
        kused,
        fmt_e(hused, 12, 4)
    );
}

fn PrintFinalStats(mem: &mut IDAMem) {
    let (mut nst, mut nni, mut nnf, mut nje, mut nre, mut nreLS, mut netf, mut ncfn) =
        (0i64, 0i64, 0i64, 0i64, 0i64, 0i64, 0i64, 0i64);

    IDAGetNumSteps(mem, &mut nst);
    IDAGetNumResEvals(mem, &mut nre);
    IDAGetNumJacEvals(mem, &mut nje);
    IDAGetNumNonlinSolvIters(mem, &mut nni);
    IDAGetNumErrTestFails(mem, &mut netf);
    IDAGetNumNonlinSolvConvFails(mem, &mut nnf);
    IDAGetNumStepSolveFails(mem, &mut ncfn);
    IDAGetNumLinResEvals(mem, &mut nreLS);

    print!("\nFinal Run Statistics: \n\n");
    println!("Number of steps                    = {}", nst);
    println!("Number of residual evaluations     = {}", nre + nreLS);
    println!("Number of Jacobian evaluations     = {}", nje);
    println!("Number of nonlinear iterations     = {}", nni);
    println!("Number of error test failures      = {}", netf);
    println!("Number of nonlinear conv. failures = {}", nnf);
    println!("Number of step solver failures     = {}", ncfn);
}

fn main() {
    /* Create the SUNDIALS context object for this simulation */
    let sunctx = SUNContext_Create();

    let mut id = N_VNew_Serial(NEQ as i64, &sunctx);
    let mut yy = N_VClone(&id);
    let mut yp = N_VClone(&id);
    let mut q = N_VNew_Serial(1, &sunctx);

    /* User data */
    let data = SlData {
        a: 0.5,  /* half-length of crank */
        J1: 1.0, /* crank moment of inertia */
        m2: 1.0, /* mass of connecting rod */
        m1: 1.0,
        J2: 2.0,             /* moment of inertia of connecting rod */
        params: [1.0, 1.0],  /* spring constant; damper constant */
        l0: 1.0,             /* spring free length */
        F: 1.0,              /* external constant force */
    };
    let _ = data.m1;

    N_VConst(ONE, &mut id);
    id.data[9] = ZERO;
    id.data[8] = ZERO;
    id.data[7] = ZERO;
    id.data[6] = ZERO;

    /* Consistent IC*/
    setIC(&mut yy, &mut yp, &data);

    /* IDAS initialization */
    let mut mem = IDACreate(&sunctx);
    IDAInit(&mut mem, ressc, TBEGIN, &yy, &yp);
    IDASStolerances(&mut mem, RTOLF, ATOLF);
    IDASetUserData(&mut mem, Some(Box::new(data)));
    IDASetId(&mut mem, Some(&id));
    IDASetSuppressAlg(&mut mem, true);
    IDASetMaxNumSteps(&mut mem, 20000);

    /* Create dense SUNMatrix for use in linear solves */
    let a = SUNDenseMatrix(NEQ as i64, NEQ as i64, &sunctx);

    /* Create dense SUNLinearSolver object */
    let ls = SUNLinSol_Dense(&yy, &a, &sunctx);

    /* Attach the matrix and linear solver */
    let retval = IDASetLinearSolver(&mut mem, ls, Some(a));
    if check_retval(retval, "IDASetLinearSolver") {
        std::process::exit(1);
    }

    N_VConst(ZERO, &mut q);
    IDAQuadInit(&mut mem, rhsQ, &q);
    IDAQuadSStolerances(&mut mem, RTOLQ, ATOLQ);
    IDASetQuadErrCon(&mut mem, true);

    PrintHeader(RTOLF, ATOLF);

    /* Print initial states */
    PrintOutput(&mut mem, 0.0, &yy);

    /* Perform forward run */
    let mut tout = TEND / NOUT as f64;
    let mut tret = 0.0;

    loop {
        let retval = IDASolve(&mut mem, tout, &mut tret, &mut yy, &mut yp, IDA_NORMAL);
        if check_retval(retval, "IDASolve") {
            std::process::exit(1);
        }

        PrintOutput(&mut mem, tret, &yy);

        tout += TEND / NOUT as f64;

        if tret > TEND {
            break;
        }
    }

    PrintFinalStats(&mut mem);

    IDAGetQuad(&mem, &mut tret, &mut q);
    print!("--------------------------------------------\n");
    println!("  G = {}", fmt_f(q.data[0], 24, 16));
    print!("--------------------------------------------\n\n");

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
