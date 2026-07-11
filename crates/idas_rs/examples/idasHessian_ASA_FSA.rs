/* -----------------------------------------------------------------
 * Translated from examples/idas/serial/idasHessian_ASA_FSA.c
 * (IDAS 7.7.0)
 * Programmer(s): Radu Serban and Cosmin Petra @ LLNL
 *
 * Hessian using adjoint sensitivity example problem.
 *
 * This simple example problem for IDAS, due to Robertson,
 * is from chemical kinetics, and consists of the following three
 * equations:
 *
 *   y1' + p1 * y1 - p2 * y2 * y3             = 0
 *   y2' - p1 * y1 + p2 * y2 * y3 + p3 * y2^2 = 0
 *   y1 + y2 + y3 - 1                         = 0
 *
 *        [1]        [-p1]
 *   y(0)=[0]  y'(0)=[ p1]   p1 = 0.04   p2 = 1e4   p3 = 1e07
 *        [0]        [ 0 ]
 *
 *       80
 *      /
 *  G = | 0.5 * (y1^2 + y2^2 + y3^2) dt
 *      /
 *      0
 * Compute the gradient (using FSA and ASA) and Hessian (FSA over ASA)
 * of G with respect to parameters p1 and p2.
 *
 * Reference: D.B. Ozyurt and P.I. Barton, SISC 26(5) 1725-1743, 2005.
 *
 * Error handling was suppressed for code readability reasons.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use idas_rs::sundials_utils::fmt_e;
use idas_rs::*;

/* Problem Constants */
const NEQ: usize = 3; /* number of equations                  */
const NP: usize = 2; /* number of sensitivities              */

const T0: f64 = 0.0; /* Initial time. */
const TF: f64 = 80.0; /* Final time. */

/* Tolerances */
const RTOL: f64 = 1e-08; /* scalar relative tolerance            */
const ATOL: f64 = 1e-10; /* vector absolute tolerance components */
const RTOLA: f64 = 1e-08; /* for adjoint integration              */
const ATOLA: f64 = 1e-08; /* for adjoint integration              */

/* Parameters */
const P1: f64 = 0.04;
const P2: f64 = 1.0e4;
const P3: f64 = 3.0e7;

/* Predefined consts */
const HALF: f64 = 0.5;
const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;

/* User defined struct */
#[derive(Clone)]
struct RobData {
    p: [f64; 3],
}

/* residual for forward problem */
fn res(_tres: f64, yy: &NVector, yp: &NVector, rr: &mut NVector, user_data: &mut UserData) -> i32 {
    let y1 = yy.data[0];
    let y2 = yy.data[1];
    let y3 = yy.data[2];
    let yp1 = yp.data[0];
    let yp2 = yp.data[1];
    let rval = &mut rr.data;

    let data = user_data.as_ref().unwrap().downcast_ref::<RobData>().unwrap();
    let p1 = data.p[0];
    let p2 = data.p[1];
    let p3 = data.p[2];

    rval[0] = p1 * y1 - p2 * y2 * y3;
    rval[1] = -rval[0] + p3 * y2 * y2 + yp2;
    rval[0] += yp1;
    rval[2] = y1 + y2 + y3 - 1.0;

    0
}

#[allow(clippy::too_many_arguments)]
fn resS(_Ns: i32, _t: f64, yy: &NVector, _yp: &NVector, _resval: &NVector, yyS: &[NVector],
        ypS: &[NVector], resvalS: &mut [NVector], user_data: &mut UserData, _tmp1: &mut NVector,
        _tmp2: &mut NVector, _tmp3: &mut NVector) -> i32 {
    let data = user_data.as_ref().unwrap().downcast_ref::<RobData>().unwrap();
    let p1 = data.p[0];
    let p2 = data.p[1];
    let p3 = data.p[2];

    let y1 = yy.data[0];
    let y2 = yy.data[1];
    let y3 = yy.data[2];

    for is in 0..NP {
        let s1 = yyS[is].data[0];
        let s2 = yyS[is].data[1];
        let s3 = yyS[is].data[2];

        let sd1 = ypS[is].data[0];
        let sd2 = ypS[is].data[1];

        let mut rs1 = sd1 + p1 * s1 - p2 * y3 * s2 - p2 * y2 * s3;
        let mut rs2 = sd2 - p1 * s1 + p2 * y3 * s2 + p2 * y2 * s3 + TWO * p3 * y2 * s2;
        let rs3 = s1 + s2 + s3;

        match is {
            0 => {
                rs1 += y1;
                rs2 -= y1;
            }
            1 => {
                rs1 -= y2 * y3;
                rs2 += y2 * y3;
            }
            _ => {}
        }

        resvalS[is].data[0] = rs1;
        resvalS[is].data[1] = rs2;
        resvalS[is].data[2] = rs3;
    }

    0
}

fn rhsQ(_t: f64, yy: &NVector, _yp: &NVector, qdot: &mut NVector, _user_data: &mut UserData) -> i32 {
    let y1 = yy.data[0];
    let y2 = yy.data[1];
    let y3 = yy.data[2];
    qdot.data[0] = HALF * (y1 * y1 + y2 * y2 + y3 * y3);

    0
}

#[allow(clippy::too_many_arguments)]
fn rhsQS(_Ns: i32, _t: f64, yy: &NVector, _yp: &NVector, yyS: &[NVector], _ypS: &[NVector],
         _rrQ: &NVector, rhsvalQS: &mut [NVector], _user_data: &mut UserData,
         _yytmp: &mut NVector, _yptmp: &mut NVector, _tmpQS: &mut NVector) -> i32 {
    let y1 = yy.data[0];
    let y2 = yy.data[1];
    let y3 = yy.data[2];

    /* 1st sensitivity RHS */
    let s1 = yyS[0].data[0];
    let s2 = yyS[0].data[1];
    let s3 = yyS[0].data[2];
    rhsvalQS[0].data[0] = y1 * s1 + y2 * s2 + y3 * s3;

    /* 2nd sensitivity RHS */
    let s1 = yyS[1].data[0];
    let s2 = yyS[1].data[1];
    let s3 = yyS[1].data[2];
    rhsvalQS[1].data[0] = y1 * s1 + y2 * s2 + y3 * s3;

    0
}

/* Residuals for adjoint model. */
#[allow(clippy::too_many_arguments)]
fn resBS1(_tt: f64, yy: &NVector, _yp: &NVector, yyS: &[NVector], _ypS: &[NVector], yyB: &NVector,
          ypB: &NVector, rrBS: &mut NVector, user_dataB: &mut UserData) -> i32 {
    let data = user_dataB.as_ref().unwrap().downcast_ref::<RobData>().unwrap();

    /* The parameters. */
    let p1 = data.p[0];
    let p2 = data.p[1];
    let p3 = data.p[2];

    /* The y vector. */
    let y1 = yy.data[0];
    let y2 = yy.data[1];
    let y3 = yy.data[2];

    /* The lambda vector. */
    let l1 = yyB.data[0];
    let l2 = yyB.data[1];
    let l3 = yyB.data[2];
    /* The mu vector. */
    let m1 = yyB.data[3];
    let m2 = yyB.data[4];
    let m3 = yyB.data[5];

    /* The lambda dot vector. */
    let lp1 = ypB.data[0];
    let lp2 = ypB.data[1];
    /* The mu dot vector. */
    let mp1 = ypB.data[3];
    let mp2 = ypB.data[4];

    /* The sensitivity with respect to p1 */
    let s1 = yyS[0].data[0];
    let s2 = yyS[0].data[1];
    let s3 = yyS[0].data[2];

    /* Temporary variables */
    let l21 = l2 - l1;

    rrBS.data[0] = lp1 + p1 * l21 - l3 + y1;
    rrBS.data[1] = lp2 - p2 * y3 * l21 - TWO * p3 * y2 * l2 - l3 + y2;
    rrBS.data[2] = -p2 * y2 * l21 - l3 + y3;

    rrBS.data[3] = mp1 + p1 * (-m1 + m2) - m3 + l21 + s1;
    rrBS.data[4] = mp2 + p2 * y3 * m1 - (p2 * y3 + TWO * p3 * y2) * m2 - m3 + p2 * s3 * l1
        - (TWO * p3 * s2 + p2 * s3) * l2
        + s2;
    rrBS.data[5] = p2 * y2 * (m1 - m2) - m3 - p2 * s2 * l21 + s3;

    0
}

#[allow(clippy::too_many_arguments)]
fn rhsQBS1(_tt: f64, yy: &NVector, _yp: &NVector, yyS: &[NVector], _ypS: &[NVector],
           yyB: &NVector, _ypB: &NVector, rhsBQS: &mut NVector, _user_dataB: &mut UserData) -> i32 {
    /* The y vector */
    let y1 = yy.data[0];
    let y2 = yy.data[1];
    let y3 = yy.data[2];

    /* The lambda vector. */
    let l1 = yyB.data[0];
    let l2 = yyB.data[1];

    /* The mu vector. */
    let m1 = yyB.data[3];
    let m2 = yyB.data[4];

    /* The sensitivity with respect to p1 */
    let s1 = yyS[0].data[0];
    let s2 = yyS[0].data[1];
    let s3 = yyS[0].data[2];

    /* Temporary variables */
    let l21 = l2 - l1;

    rhsBQS.data[0] = -y1 * l21;
    rhsBQS.data[1] = y2 * y3 * l21;

    rhsBQS.data[2] = y1 * (m1 - m2) - s1 * l21;
    rhsBQS.data[3] = y2 * y3 * (m2 - m1) + (y3 * s2 + y2 * s3) * l21;

    0
}

#[allow(clippy::too_many_arguments)]
fn resBS2(_tt: f64, yy: &NVector, _yp: &NVector, yyS: &[NVector], _ypS: &[NVector], yyB: &NVector,
          ypB: &NVector, rrBS: &mut NVector, user_dataB: &mut UserData) -> i32 {
    let data = user_dataB.as_ref().unwrap().downcast_ref::<RobData>().unwrap();

    /* The parameters. */
    let p1 = data.p[0];
    let p2 = data.p[1];
    let p3 = data.p[2];

    /* The y vector. */
    let y1 = yy.data[0];
    let y2 = yy.data[1];
    let y3 = yy.data[2];

    /* The lambda vector. */
    let l1 = yyB.data[0];
    let l2 = yyB.data[1];
    let l3 = yyB.data[2];
    /* The mu vector. */
    let m1 = yyB.data[3];
    let m2 = yyB.data[4];
    let m3 = yyB.data[5];

    /* The lambda dot vector. */
    let lp1 = ypB.data[0];
    let lp2 = ypB.data[1];

    /* The mu dot vector. */
    let mp1 = ypB.data[3];
    let mp2 = ypB.data[4];

    /* The sensitivity with respect to p2 */
    let s1 = yyS[1].data[0];
    let s2 = yyS[1].data[1];
    let s3 = yyS[1].data[2];

    /* Temporary variables */
    let l21 = l2 - l1;

    rrBS.data[0] = lp1 + p1 * l21 - l3 + y1;
    rrBS.data[1] = lp2 - p2 * y3 * l21 - TWO * p3 * y2 * l2 - l3 + y2;
    rrBS.data[2] = -p2 * y2 * l21 - l3 + y3;

    rrBS.data[3] = mp1 + p1 * (-m1 + m2) - m3 + s1;
    rrBS.data[4] = mp2 + p2 * y3 * m1 - (p2 * y3 + TWO * p3 * y2) * m2 - m3
        + (y3 + p2 * s3) * l1
        - (y3 + TWO * p3 * s2 + p2 * s3) * l2
        + s2;
    rrBS.data[5] = p2 * y2 * (m1 - m2) - m3 - (y2 + p2 * s2) * l21 + s3;

    0
}

#[allow(clippy::too_many_arguments)]
fn rhsQBS2(_tt: f64, yy: &NVector, _yp: &NVector, yyS: &[NVector], _ypS: &[NVector],
           yyB: &NVector, _ypB: &NVector, rhsBQS: &mut NVector, _user_dataB: &mut UserData) -> i32 {
    /* The y vector */
    let y1 = yy.data[0];
    let y2 = yy.data[1];
    let y3 = yy.data[2];

    /* The lambda vector. */
    let l1 = yyB.data[0];
    let l2 = yyB.data[1];

    /* The mu vector. */
    let m1 = yyB.data[3];
    let m2 = yyB.data[4];

    /* The sensitivity with respect to p2 */
    let s1 = yyS[1].data[0];
    let s2 = yyS[1].data[1];
    let s3 = yyS[1].data[2];

    /* Temporary variables */
    let l21 = l2 - l1;

    rhsBQS.data[0] = -y1 * l21;
    rhsBQS.data[1] = y2 * y3 * l21;

    rhsBQS.data[2] = y1 * (m1 - m2) - s1 * l21;
    rhsBQS.data[3] = y2 * y3 * (m2 - m1) + (y3 * s2 + y2 * s3) * l21;

    0
}

fn check_retval(retval: i32, funcname: &str) -> bool {
    if retval < 0 {
        eprintln!("\nSUNDIALS_ERROR: {}() failed with retval = {}\n", funcname, retval);
        return true;
    }
    false
}

fn main() {
    /* Create the SUNDIALS context object for this simulation. */
    let sunctx = SUNContext_Create();

    /* Print problem description */
    print!("\nAdjoint Sensitivity Example for Chemical Kinetics\n");
    print!("---------------------------------------------------------\n");
    print!("DAE: dy1/dt + p1*y1 - p2*y2*y3 = 0\n");
    print!("     dy2/dt - p1*y1 + p2*y2*y3 + p3*(y2)^2 = 0\n");
    print!("               y1  +  y2  +  y3 = 0\n\n");
    print!("Find dG/dp and d^2G/dp^2, where p=[p1,p2] for\n");
    print!("     G = int_t0^tB0 g(t,p,y) dt\n");
    print!("     g(t,p,y) = y3\n\n\n");

    /* Allocate and initialize user data. */
    let mut data = RobData { p: [P1, P2, P3] };

    /* Consistent IC */
    let mut yy = N_VNew_Serial(NEQ as i64, &sunctx);
    let mut yp = N_VClone(&yy);
    yy.data[0] = ONE;
    yy.data[1] = ZERO;
    yy.data[2] = ZERO;
    yp.data[0] = -P1;
    yp.data[1] = P1;
    yp.data[2] = 0.0;

    let mut q = N_VNew_Serial(1, &sunctx);
    N_VConst(ZERO, &mut q);

    let mut yyS: Vec<NVector> = (0..NP).map(|_| N_VClone(&yy)).collect();
    let mut ypS: Vec<NVector> = (0..NP).map(|_| N_VClone(&yp)).collect();
    N_VConst(ZERO, &mut yyS[0]);
    N_VConst(ZERO, &mut yyS[1]);
    N_VConst(ZERO, &mut ypS[0]);
    N_VConst(ZERO, &mut ypS[1]);

    let mut qS: Vec<NVector> = (0..NP).map(|_| N_VClone(&q)).collect();
    for s in qS.iter_mut() {
        N_VConst(ZERO, s);
    }

    let mut ida_mem = IDACreate(&sunctx);

    let ti = T0;
    IDAInit(&mut ida_mem, res, ti, &yy, &yp);

    /* Forward problem's setup. */
    IDASStolerances(&mut ida_mem, RTOL, ATOL);

    /* Create dense SUNMatrix for use in linear solves */
    let a = SUNDenseMatrix(NEQ as i64, NEQ as i64, &sunctx);

    /* Create dense SUNLinearSolver object */
    let ls = SUNLinSol_Dense(&yy, &a, &sunctx);

    /* Attach the matrix and linear solver */
    let retval = IDASetLinearSolver(&mut ida_mem, ls, Some(a));
    if check_retval(retval, "IDASetLinearSolver") {
        std::process::exit(1);
    }

    IDASetUserData(&mut ida_mem, Some(Box::new(data.clone())));
    IDASetMaxNumSteps(&mut ida_mem, 1500);

    /* Quadrature's setup. */
    IDAQuadInit(&mut ida_mem, rhsQ, &q);
    IDAQuadSStolerances(&mut ida_mem, RTOL, ATOL);
    IDASetQuadErrCon(&mut ida_mem, true);

    /* Sensitivity's setup. */
    IDASensInit(&mut ida_mem, NP as i32, IDA_SIMULTANEOUS, Some(resS), &yyS, &ypS);
    IDASensEEtolerances(&mut ida_mem);
    IDASetSensErrCon(&mut ida_mem, true);

    /* Setup of quadrature's sensitivities */
    IDAQuadSensInit(&mut ida_mem, Some(rhsQS), &qS);
    IDAQuadSensEEtolerances(&mut ida_mem);
    IDASetQuadSensErrCon(&mut ida_mem, true);

    /* Initialize ASA. */
    IDAAdjInit(&mut ida_mem, 100, IDA_HERMITE);

    print!("---------------------------------------------------------\n");
    print!("Forward integration\n");
    print!("---------------------------------------------------------\n\n");

    let tf = TF;
    let mut time = 0.0;
    let mut nckp = 0;
    IDASolveF(&mut ida_mem, tf, &mut time, &mut yy, &mut yp, IDA_NORMAL, &mut nckp);

    IDAGetQuad(&ida_mem, &mut time, &mut q);
    let G = q.data[0];
    println!("     G:    {}", fmt_e(G, 12, 4));

    /* Sensitivities are needed for IC of backward problems. */
    IDAGetSensDky(&ida_mem, tf, 0, &mut yyS);
    IDAGetSensDky(&ida_mem, tf, 1, &mut ypS);

    IDAGetQuadSens(&ida_mem, &mut time, &mut qS);
    println!("   dG/dp:  {} {}", fmt_e(qS[0].data[0], 12, 4), fmt_e(qS[1].data[0], 12, 4));
    println!();

    /******************************
     * BACKWARD PROBLEM #1
     *******************************/

    /* Consistent IC. */
    let mut yyB1 = N_VNew_Serial(2 * NEQ as i64, &sunctx);
    let mut ypB1 = N_VClone(&yyB1);

    N_VConst(ZERO, &mut yyB1);
    yyB1.data[2] = yy.data[2];
    yyB1.data[5] = yyS[0].data[2];

    N_VConst(ZERO, &mut ypB1);
    ypB1.data[0] = yy.data[2] - yy.data[0];
    ypB1.data[1] = yy.data[2] - yy.data[1];
    ypB1.data[3] = yyS[0].data[2] - yyS[0].data[0];
    ypB1.data[4] = yyS[0].data[2] - yyS[0].data[1];

    let mut qB1 = N_VNew_Serial(2 * NP as i64, &sunctx);
    N_VConst(ZERO, &mut qB1);

    let mut indexB1 = 0;
    IDACreateB(&mut ida_mem, &mut indexB1);
    IDAInitBS(&mut ida_mem, indexB1, resBS1, tf, &yyB1, &ypB1);
    IDASStolerancesB(&mut ida_mem, indexB1, RTOLA, ATOLA);
    IDASetUserDataB(&mut ida_mem, indexB1, Some(Box::new(data.clone())));
    IDASetMaxNumStepsB(&mut ida_mem, indexB1, 5000);

    /* Create dense SUNMatrix for use in linear solves */
    let ab1 = SUNDenseMatrix(2 * NEQ as i64, 2 * NEQ as i64, &sunctx);

    /* Create dense SUNLinearSolver object */
    let lsb1 = SUNLinSol_Dense(&yyB1, &ab1, &sunctx);

    /* Attach the matrix and linear solver */
    let retval = IDASetLinearSolverB(&mut ida_mem, indexB1, lsb1, Some(ab1));
    if check_retval(retval, "IDASetLinearSolverB") {
        std::process::exit(1);
    }

    IDAQuadInitBS(&mut ida_mem, indexB1, rhsQBS1, &qB1);

    /******************************
     * BACKWARD PROBLEM #2
     *******************************/

    /* Consistent IC. */
    let mut yyB2 = N_VNew_Serial(2 * NEQ as i64, &sunctx);
    let mut ypB2 = N_VNew_Serial(2 * NEQ as i64, &sunctx);

    N_VConst(ZERO, &mut yyB2);
    yyB2.data[2] = yy.data[2];
    yyB2.data[5] = yyS[1].data[2];

    N_VConst(ZERO, &mut ypB2);
    ypB2.data[0] = yy.data[2] - yy.data[0];
    ypB2.data[1] = yy.data[2] - yy.data[1];
    ypB2.data[3] = yyS[1].data[2] - yyS[1].data[0];
    ypB2.data[4] = yyS[1].data[2] - yyS[1].data[1];

    let mut qB2 = N_VNew_Serial(2 * NP as i64, &sunctx);
    N_VConst(ZERO, &mut qB2);

    let mut indexB2 = 0;
    IDACreateB(&mut ida_mem, &mut indexB2);
    IDAInitBS(&mut ida_mem, indexB2, resBS2, tf, &yyB2, &ypB2);
    IDASStolerancesB(&mut ida_mem, indexB2, RTOLA, ATOLA);
    IDASetUserDataB(&mut ida_mem, indexB2, Some(Box::new(data.clone())));
    IDASetMaxNumStepsB(&mut ida_mem, indexB2, 2500);

    /* Create dense SUNMatrix for use in linear solves */
    let ab2 = SUNDenseMatrix(2 * NEQ as i64, 2 * NEQ as i64, &sunctx);

    /* Create dense SUNLinearSolver object */
    let lsb2 = SUNLinSol_Dense(&yyB2, &ab2, &sunctx);

    /* Attach the matrix and linear solver */
    let retval = IDASetLinearSolverB(&mut ida_mem, indexB2, lsb2, Some(ab2));
    if check_retval(retval, "IDASetLinearSolverB") {
        std::process::exit(1);
    }

    IDAQuadInitBS(&mut ida_mem, indexB2, rhsQBS2, &qB2);

    /* Integrate backward problems. */
    print!("---------------------------------------------------------\n");
    print!("Backward integration \n");
    print!("---------------------------------------------------------\n\n");

    IDASolveB(&mut ida_mem, ti, IDA_NORMAL);

    IDAGetB(&mut ida_mem, indexB1, &mut time, &mut yyB1, &mut ypB1);
    /*
       let mut nst = 0i64;
       IDAGetNumSteps(IDAGetAdjIDABmem(&mut ida_mem, indexB1).unwrap(), &mut nst);
       println!("at time={} \tpb 1 Num steps:{}", time, nst);
       IDAGetNumSteps(IDAGetAdjIDABmem(&mut ida_mem, indexB2).unwrap(), &mut nst);
       println!("at time={} \tpb 2 Num steps:{}\n", time, nst);
    */

    IDAGetQuadB(&mut ida_mem, indexB1, &mut time, &mut qB1);
    IDAGetQuadB(&mut ida_mem, indexB2, &mut time, &mut qB2);
    println!(
        "   dG/dp:  {} {}   (from backward pb. 1)",
        fmt_e(qB1.data[0], 12, 4),
        fmt_e(qB1.data[1], 12, 4)
    );
    println!(
        "   dG/dp:  {} {}   (from backward pb. 2)",
        fmt_e(qB2.data[0], 12, 4),
        fmt_e(qB2.data[1], 12, 4)
    );

    println!();
    print!("   H = d2G/dp2:\n");
    print!("        (1)            (2)\n");
    println!("  {}  {}", fmt_e(qB1.data[2], 12, 4), fmt_e(qB2.data[2], 12, 4));
    println!("  {}  {}", fmt_e(qB1.data[3], 12, 4), fmt_e(qB2.data[3], 12, 4));

    IDAFree(ida_mem);

    /*********************************
     * Use Finite Differences to verify
     **********************************/

    /* Perturbations are of different magnitudes as p1 and p2 are. */
    let dp1: f64 = 1.0e-3;
    let dp2: f64 = 2.5e+2;

    println!();
    print!("---------------------------------------------------------\n");
    println!(
        "Finite Differences ( dp1={} and dp2 = {} )",
        fmt_e(dp1, 6, 1),
        fmt_e(dp2, 6, 1)
    );
    print!("---------------------------------------------------------\n\n");

    let mut ida_mem = IDACreate(&sunctx);

    /********************
     * Forward FD for p1
     ********************/
    data.p[0] += dp1;

    yy.data[0] = ONE;
    yy.data[1] = ZERO;
    yy.data[2] = ZERO;
    yp.data[0] = -data.p[0];
    yp.data[1] = -yp.data[0];
    yp.data[2] = 0.0;
    N_VConst(ZERO, &mut q);
    let ti = T0;
    let tf = TF;

    IDAInit(&mut ida_mem, res, ti, &yy, &yp);

    let rtolFD: f64 = 1.0e-12;
    let atolFD: f64 = 1.0e-14;

    IDASStolerances(&mut ida_mem, rtolFD, atolFD);

    /* Create dense SUNMatrix for use in linear solves */
    let a = SUNDenseMatrix(NEQ as i64, NEQ as i64, &sunctx);

    /* Create dense SUNLinearSolver object */
    let ls = SUNLinSol_Dense(&yy, &a, &sunctx);

    /* Attach the matrix and linear solver */
    let retval = IDASetLinearSolver(&mut ida_mem, ls, Some(a));
    if check_retval(retval, "IDASetLinearSolver") {
        std::process::exit(1);
    }

    IDASetUserData(&mut ida_mem, Some(Box::new(data.clone())));
    IDASetMaxNumSteps(&mut ida_mem, 10000);

    IDAQuadInit(&mut ida_mem, rhsQ, &q);
    IDAQuadSStolerances(&mut ida_mem, rtolFD, atolFD);
    IDASetQuadErrCon(&mut ida_mem, true);

    IDASolve(&mut ida_mem, tf, &mut time, &mut yy, &mut yp, IDA_NORMAL);
    IDAGetQuad(&ida_mem, &mut time, &mut q);
    let mut Gp = q.data[0];

    /********************
     * Backward FD for p1
     ********************/
    data.p[0] -= 2.0 * dp1;
    IDASetUserData(&mut ida_mem, Some(Box::new(data.clone())));

    yy.data[0] = ONE;
    yy.data[1] = ZERO;
    yy.data[2] = ZERO;
    yp.data[0] = -data.p[0];
    yp.data[1] = -yp.data[0];
    yp.data[2] = 0.0;
    N_VConst(ZERO, &mut q);

    IDAReInit(&mut ida_mem, ti, &yy, &yp);
    IDAQuadReInit(&mut ida_mem, &q);

    IDASolve(&mut ida_mem, tf, &mut time, &mut yy, &mut yp, IDA_NORMAL);
    IDAGetQuad(&ida_mem, &mut time, &mut q);
    let mut Gm = q.data[0];

    /* Compute FD for p1. */
    let mut grdG_fwd = [0.0f64; 2];
    let mut grdG_bck = [0.0f64; 2];
    let mut grdG_cntr = [0.0f64; 2];
    grdG_fwd[0] = (Gp - G) / dp1;
    grdG_bck[0] = (G - Gm) / dp1;
    grdG_cntr[0] = (Gp - Gm) / (2.0 * dp1);
    let H11 = (Gp - 2.0 * G + Gm) / (dp1 * dp1);

    /********************
     * Forward FD for p2
     ********************/
    /*restore p1*/
    data.p[0] += dp1;
    data.p[1] += dp2;
    IDASetUserData(&mut ida_mem, Some(Box::new(data.clone())));

    yy.data[0] = ONE;
    yy.data[1] = ZERO;
    yy.data[2] = ZERO;
    yp.data[0] = -data.p[0];
    yp.data[1] = -yp.data[0];
    yp.data[2] = 0.0;
    N_VConst(ZERO, &mut q);

    IDAReInit(&mut ida_mem, ti, &yy, &yp);
    IDAQuadReInit(&mut ida_mem, &q);

    IDASolve(&mut ida_mem, tf, &mut time, &mut yy, &mut yp, IDA_NORMAL);
    IDAGetQuad(&ida_mem, &mut time, &mut q);
    Gp = q.data[0];

    /********************
     * Backward FD for p2
     ********************/
    data.p[1] -= 2.0 * dp2;
    IDASetUserData(&mut ida_mem, Some(Box::new(data.clone())));

    yy.data[0] = ONE;
    yy.data[1] = ZERO;
    yy.data[2] = ZERO;
    yp.data[0] = -data.p[0];
    yp.data[1] = -yp.data[0];
    yp.data[2] = 0.0;
    N_VConst(ZERO, &mut q);

    IDAReInit(&mut ida_mem, ti, &yy, &yp);
    IDAQuadReInit(&mut ida_mem, &q);

    IDASolve(&mut ida_mem, tf, &mut time, &mut yy, &mut yp, IDA_NORMAL);
    IDAGetQuad(&ida_mem, &mut time, &mut q);
    Gm = q.data[0];

    /* Compute FD for p2. */
    grdG_fwd[1] = (Gp - G) / dp2;
    grdG_bck[1] = (G - Gm) / dp2;
    grdG_cntr[1] = (Gp - Gm) / (2.0 * dp2);
    let H22 = (Gp - 2.0 * G + Gm) / (dp2 * dp2);

    println!();
    println!("   dG/dp:  {}  {}   (fwd FD)", fmt_e(grdG_fwd[0], 12, 4), fmt_e(grdG_fwd[1], 12, 4));
    println!("           {}  {}   (bck FD)", fmt_e(grdG_bck[0], 12, 4), fmt_e(grdG_bck[1], 12, 4));
    println!(
        "           {}  {}   (cntr FD)",
        fmt_e(grdG_cntr[0], 12, 4),
        fmt_e(grdG_cntr[1], 12, 4)
    );
    println!();
    println!("  H(1,1):  {}", fmt_e(H11, 12, 4));
    println!("  H(2,2):  {}", fmt_e(H22, 12, 4));

    /* Free memory (RAII) */
    IDAFree(ida_mem);
}
