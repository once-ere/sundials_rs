/* -----------------------------------------------------------------
 * Translated from examples/idas/serial/idasSlCrank_FSA_dns.c
 * (IDAS 7.7.0)
 * Programmer: Radu Serban and Cosmin Petra @ LLNL
 *
 * Simulation of a slider-crank mechanism modelled with 3 generalized
 * coordinates: crank angle, connecting bar angle, and slider location.
 * The mechanism moves under the action of a constant horizontal
 * force applied to the connecting rod and a spring-damper connecting
 * the crank and connecting rod.
 *
 * The equations of motion are formulated as a system of stabilized
 * index-2 DAEs (Gear-Gupta-Leimkuhler formulation).
 *
 * IDAS also computes sensitivities with respect to the problem
 * parameters k (spring constant) and c (damper constant) of the
 * kinetic energy:
 *   G = int_t0^tend g(t,y,p) dt,
 * where
 *   g(t,y,p) = 0.5*J1*v1^2 + 0.5*J2*v3^2 + 0.5*m2*v2^2
 *
 * (C attaches one UserData pointer and mutates params through it
 * between the finite-difference reruns; here the example keeps a
 * local copy and re-attaches a clone after each mutation — the
 * solver sees the same values at the same times.  Because the
 * sensitivity residual is the INTERNAL DQ (resS = NULL), the
 * parameters live in the FSAUserData wrapper so the DQ perturbations
 * of p reach ressc — see sundials_types.rs.)
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use idas_rs::sundials_utils::{fmt_e, fmt_f};
use idas_rs::*;

/* Problem Constants */

const NEQ: usize = 10;
const NP: usize = 2;

const TBEGIN: f64 = 0.0;
const TEND: f64 = 10.000;

const RTOLF: f64 = 1.0e-06;
const ATOLF: f64 = 1.0e-07;

const RTOLQ: f64 = 1.0e-06;
const ATOLQ: f64 = 1.0e-08;

const RTOLFD: f64 = 1.0e-06;
const ATOLFD: f64 = 1.0e-08;

const ZERO: f64 = 0.00;
const HALF: f64 = 0.50;
const ONE: f64 = 1.00;
const TWO: f64 = 2.00;
const FOUR: f64 = 4.00;

/* C UserData splits into FSAUserData { p: params, user: SlConsts }
   (internal-DQ convention) */
#[derive(Clone)]
struct SlConsts {
    a: f64,
    J1: f64,
    J2: f64,
    m1: f64,
    m2: f64,
    l0: f64,
    F: f64,
}

/* view of the C struct: constants + the parameter array */
struct SlData<'a> {
    c: &'a SlConsts,
    params: &'a [f64],
}

fn sl_view(user_data: &UserData) -> SlData<'_> {
    let f = user_data.as_ref().unwrap().downcast_ref::<FSAUserData>().unwrap();
    SlData { c: f.user.downcast_ref::<SlConsts>().unwrap(), params: &f.p }
}

fn force(yy: &NVector, data: &SlData) -> [f64; 3] {
    let a = data.c.a;
    let k = data.params[0];
    let c = data.params[1];
    let l0 = data.c.l0;
    let F = data.c.F;

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

    let a = data.c.a;
    let J1 = data.c.J1;
    let m2 = data.c.m2;
    let J2 = data.c.J2;

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
    let data = sl_view(user_data);

    let a = data.c.a;
    let J1 = data.c.J1;
    let m2 = data.c.m2;
    let J2 = data.c.J2;

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

    let Q = force(yy, &data);

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
    let data = sl_view(user_data);
    let J1 = data.c.J1;
    let m2 = data.c.m2;
    let J2 = data.c.J2;

    let v1 = yy.data[3];
    let v2 = yy.data[4];
    let v3 = yy.data[5];

    qdot.data[0] = HALF * (J1 * v1 * v1 + m2 * v2 * v2 + J2 * v3 * v3);

    0
}

#[allow(clippy::too_many_arguments)]
fn rhsQS(_Ns: i32, _t: f64, yy: &NVector, _yp: &NVector, yyS: &[NVector], _ypS: &[NVector],
         _rrQ: &NVector, rhsvalQS: &mut [NVector], user_data: &mut UserData,
         _yytmp: &mut NVector, _yptmp: &mut NVector, _tmpQS: &mut NVector) -> i32 {
    let data = sl_view(user_data);

    let J1 = data.c.J1;
    let m2 = data.c.m2;
    let J2 = data.c.J2;

    let v1 = yy.data[3];
    let v2 = yy.data[4];
    let v3 = yy.data[5];

    /* Sensitivities of v. */
    let s1 = yyS[0].data[3];
    let s2 = yyS[0].data[4];
    let s3 = yyS[0].data[5];

    rhsvalQS[0].data[0] = J1 * v1 * s1 + m2 * v2 * s2 + J2 * v3 * s3;

    let s1 = yyS[1].data[3];
    let s2 = yyS[1].data[4];
    let s3 = yyS[1].data[5];

    rhsvalQS[1].data[0] = J1 * v1 * s1 + m2 * v2 * s2 + J2 * v3 * s3;

    0
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

fn check_retval(retval: i32, funcname: &str) -> bool {
    if retval < 0 {
        eprintln!("\nSUNDIALS_ERROR: {}() failed with retval = {}\n", funcname, retval);
        return true;
    }
    false
}

/*
 *--------------------------------------------------------------------
 * Main Program
 *--------------------------------------------------------------------
 */
fn main() {
    /* Create the SUNDIALS context object for this simulation */
    let sunctx = SUNContext_Create();

    let mut id = N_VNew_Serial(NEQ as i64, &sunctx);
    let mut yy = N_VClone(&id);
    let mut yp = N_VClone(&id);
    let mut q = N_VNew_Serial(1, &sunctx);

    let mut yyS: Vec<NVector> = (0..NP).map(|_| N_VClone(&yy)).collect();
    let mut ypS: Vec<NVector> = (0..NP).map(|_| N_VClone(&yp)).collect();
    let mut qS: Vec<NVector> = (0..NP).map(|_| N_VClone(&q)).collect();

    let consts = SlConsts {
        a: 0.5,  /* half-length of crank */
        J1: 1.0, /* crank moment of inertia */
        m2: 1.0, /* mass of connecting rod */
        m1: 1.0,
        J2: 2.0, /* moment of inertia of connecting rod */
        l0: 1.0, /* spring free length */
        F: 1.0,  /* external constant force */
    };
    let _ = consts.m1;
    let mut params: [f64; 2] = [1.0, 1.0]; /* spring constant; damper constant */
    let wrap = |params: &[f64; 2]| -> UserData {
        Some(Box::new(FSAUserData { p: params.to_vec(), user: Box::new(consts.clone()) }))
    };

    N_VConst(ONE, &mut id);
    id.data[9] = ZERO;
    id.data[8] = ZERO;
    id.data[7] = ZERO;
    id.data[6] = ZERO;

    print!("\nSlider-Crank example for IDAS:\n");

    /* Consistent IC*/
    setIC(&mut yy, &mut yp, &SlData { c: &consts, params: &params });

    for is in 0..NP {
        N_VConst(ZERO, &mut yyS[is]);
        N_VConst(ZERO, &mut ypS[is]);
    }

    /* IDA initialization */
    let mut mem = IDACreate(&sunctx);
    IDAInit(&mut mem, ressc, TBEGIN, &yy, &yp);
    IDASStolerances(&mut mem, RTOLF, ATOLF);
    IDASetUserData(&mut mem, wrap(&params));
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

    IDASensInit(&mut mem, NP as i32, IDA_SIMULTANEOUS, None, &yyS, &ypS);
    let pbar = [params[0], params[1]];
    IDASetSensParams(&mut mem, Some(&params), Some(&pbar), None);
    IDASensEEtolerances(&mut mem);
    IDASetSensErrCon(&mut mem, true);

    N_VConst(ZERO, &mut q);
    IDAQuadInit(&mut mem, rhsQ, &q);
    IDAQuadSStolerances(&mut mem, RTOLQ, ATOLQ);
    IDASetQuadErrCon(&mut mem, true);

    for s in qS.iter_mut() {
        N_VConst(ZERO, s);
    }
    IDAQuadSensInit(&mut mem, Some(rhsQS), &qS);
    let atolS = [ATOLQ; NP];
    IDAQuadSensSStolerances(&mut mem, RTOLQ, &atolS);
    IDASetQuadSensErrCon(&mut mem, true);

    /* Perform forward run */
    print!("\nForward integration ... ");

    let mut tret = 0.0;
    let retval = IDASolve(&mut mem, TEND, &mut tret, &mut yy, &mut yp, IDA_NORMAL);
    if check_retval(retval, "IDASolve") {
        std::process::exit(1);
    }

    print!("done!\n");

    PrintFinalStats(&mut mem);

    IDAGetQuad(&mem, &mut tret, &mut q);
    print!("--------------------------------------------\n");
    println!("  G = {}", fmt_f(q.data[0], 24, 16));
    print!("--------------------------------------------\n\n");

    IDAGetQuadSens(&mem, &mut tret, &mut qS);
    print!("-------------F O R W A R D------------------\n");
    println!("   dG/dp:  {} {}", fmt_e(qS[0].data[0], 12, 4), fmt_e(qS[1].data[0], 12, 4));
    print!("--------------------------------------------\n\n");

    IDAFree(mem);

    /* Finite differences for dG/dp */
    let dp = 0.00001;
    params[0] = ONE;
    params[1] = ONE;

    let mut mem = IDACreate(&sunctx);

    setIC(&mut yy, &mut yp, &SlData { c: &consts, params: &params });
    IDAInit(&mut mem, ressc, TBEGIN, &yy, &yp);
    IDASStolerances(&mut mem, RTOLFD, ATOLFD);
    IDASetUserData(&mut mem, wrap(&params));
    IDASetId(&mut mem, Some(&id));
    IDASetSuppressAlg(&mut mem, true);

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

    IDASolve(&mut mem, TEND, &mut tret, &mut yy, &mut yp, IDA_NORMAL);

    IDAGetQuad(&mem, &mut tret, &mut q);
    let G = q.data[0];
    /*print!("  G  ={}\n", fmt_e(q.data[0], 12, 6));*/

    let mut Gm = [0.0f64; 2];
    let mut Gp = [0.0f64; 2];

    /******************************
     * BACKWARD for k
     ******************************/
    params[0] -= dp;
    setIC(&mut yy, &mut yp, &SlData { c: &consts, params: &params });
    IDASetUserData(&mut mem, wrap(&params));

    IDAReInit(&mut mem, TBEGIN, &yy, &yp);

    N_VConst(ZERO, &mut q);
    IDAQuadReInit(&mut mem, &q);

    IDASolve(&mut mem, TEND, &mut tret, &mut yy, &mut yp, IDA_NORMAL);
    IDAGetQuad(&mem, &mut tret, &mut q);
    Gm[0] = q.data[0];
    /*print!("Gm[0]={}\n", fmt_e(q.data[0], 12, 6));*/

    /****************************
     * FORWARD for k *
     ****************************/
    params[0] += TWO * dp;
    setIC(&mut yy, &mut yp, &SlData { c: &consts, params: &params });
    IDASetUserData(&mut mem, wrap(&params));
    IDAReInit(&mut mem, TBEGIN, &yy, &yp);

    N_VConst(ZERO, &mut q);
    IDAQuadReInit(&mut mem, &q);

    IDASolve(&mut mem, TEND, &mut tret, &mut yy, &mut yp, IDA_NORMAL);
    IDAGetQuad(&mem, &mut tret, &mut q);
    Gp[0] = q.data[0];
    /*print!("Gp[0]={}\n", fmt_e(q.data[0], 12, 6));*/

    /* Backward for c */
    params[0] = ONE;
    params[1] -= dp;
    setIC(&mut yy, &mut yp, &SlData { c: &consts, params: &params });
    IDASetUserData(&mut mem, wrap(&params));
    IDAReInit(&mut mem, TBEGIN, &yy, &yp);

    N_VConst(ZERO, &mut q);
    IDAQuadReInit(&mut mem, &q);

    IDASolve(&mut mem, TEND, &mut tret, &mut yy, &mut yp, IDA_NORMAL);
    IDAGetQuad(&mem, &mut tret, &mut q);
    Gm[1] = q.data[0];

    /* Forward for c */
    params[1] += TWO * dp;
    setIC(&mut yy, &mut yp, &SlData { c: &consts, params: &params });
    IDASetUserData(&mut mem, wrap(&params));
    IDAReInit(&mut mem, TBEGIN, &yy, &yp);

    N_VConst(ZERO, &mut q);
    IDAQuadReInit(&mut mem, &q);

    IDASolve(&mut mem, TEND, &mut tret, &mut yy, &mut yp, IDA_NORMAL);
    IDAGetQuad(&mem, &mut tret, &mut q);
    Gp[1] = q.data[0];

    IDAFree(mem);

    print!("\n\n   Checking using Finite Differences \n\n");

    print!("---------------BACKWARD------------------\n");
    println!("   dG/dp:  {} {}", fmt_e((G - Gm[0]) / dp, 12, 4), fmt_e((G - Gm[1]) / dp, 12, 4));
    print!("-----------------------------------------\n\n");

    print!("---------------FORWARD-------------------\n");
    println!("   dG/dp:  {} {}", fmt_e((Gp[0] - G) / dp, 12, 4), fmt_e((Gp[1] - G) / dp, 12, 4));
    print!("-----------------------------------------\n\n");

    print!("--------------CENTERED-------------------\n");
    println!(
        "   dG/dp:  {} {}",
        fmt_e((Gp[0] - Gm[0]) / (TWO * dp), 12, 4),
        fmt_e((Gp[1] - Gm[1]) / (TWO * dp), 12, 4)
    );
    print!("-----------------------------------------\n\n");
}
