/* -----------------------------------------------------------------
 * Translated from examples/idas/serial/idasRoberts_ASAi_dns.c
 * (IDAS 7.7.0)
 * Programmer(s): Radu Serban and Cosmin Petra @ LLNL
 *
 * Adjoint sensitivity example problem.
 *
 * This simple example problem for IDAS, due to Robertson,
 * is from chemical kinetics, and consists of the following three
 * equations:
 *
 *      dy1/dt + p1*y1 - p2*y2*y3            = 0
 *      dy2/dt - p1*y1 + p2*y2*y3 + p3*y2**2 = 0
 *                 y1  +  y2  +  y3  -  1    = 0
 *
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1, y2 = y3 = 0.The reaction rates are: p1=0.04,
 * p2=1e4, and p3=3e7
 *
 * It uses a scalar relative tolerance and a vector absolute
 * tolerance.
 *
 * IDAS can also compute sensitivities with respect to
 * the problem parameters p1, p2, and p3 of the following quantity:
 *   G = int_t0^t1 g(t,p,y) dt
 * where
 *   g(t,p,y) = y3
 *
 * The gradient dG/dp is obtained as:
 *   dG/dp = int_t0^t1 (g_p - lambda^T F_p ) dt -
 *           lambda^T*F_y'*y_p | _t0^t1
 *         = int_t0^t1 (lambda^T*F_p) dt
 * where lambda and are solutions of the adjoint system:
 *   d(lambda^T * F_y' )/dt -lambda^T F_y = -g_y
 *
 * During the backward integration, IDAS also evaluates G as
 *   G = - phi(t0)
 * where
 *   d(phi)/dt = g(t,y,p)
 *   phi(t1) = 0
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use idas_rs::sundials_utils::fmt_e;
use idas_rs::*;

/* Problem Constants */

const NEQ: usize = 3; /* number of equations                  */

const RTOL: f64 = 1e-06; /* scalar relative tolerance            */

const ATOL1: f64 = 1e-08; /* vector absolute tolerance components */
const ATOL2: f64 = 1e-12;
const ATOL3: f64 = 1e-08;

const ATOLA: f64 = 1e-08; /* absolute tolerance for adjoint vars. */
const ATOLQ: f64 = 1e-06; /* absolute tolerance for quadratures   */

const T0: f64 = 0.0; /* initial time                         */
const TOUT: f64 = 4e10; /* final time                           */

const TB1: f64 = 50.0; /* starting point for adjoint problem   */
const TB2: f64 = TOUT; /* starting point for adjoint problem   */

const T1B: f64 = 49.0; /* for IDACalcICB                       */

const STEPS: i64 = 100; /* number of steps between check points */

const NP: usize = 3; /* number of problem parameters         */

const ONE: f64 = 1.0;
const ZERO: f64 = 0.0;

/* Type : UserData */
#[derive(Clone)]
struct RobData {
    p: [f64; 3],
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY IDAS
 *--------------------------------------------------------------------
 */

/*
 * f routine. Compute f(t,y).
 */
fn res(_t: f64, yy: &NVector, yp: &NVector, resval: &mut NVector, user_data: &mut UserData) -> i32 {
    let y1 = yy.data[0];
    let y2 = yy.data[1];
    let y3 = yy.data[2];
    let yp1 = yp.data[0];
    let yp2 = yp.data[1];
    let rval = &mut resval.data;

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

/*
 * Jacobian routine. Compute J(t,y).
 */
#[allow(clippy::too_many_arguments)]
fn Jac(_t: f64, cj: f64, yy: &NVector, _yp: &NVector, _resvec: &NVector, J: &mut SUNMatrix,
       user_data: &mut UserData, _tmp1: &mut NVector, _tmp2: &mut NVector,
       _tmp3: &mut NVector) -> i32 {
    let y2 = yy.data[1];
    let y3 = yy.data[2];

    let data = user_data.as_ref().unwrap().downcast_ref::<RobData>().unwrap();
    let p1 = data.p[0];
    let p2 = data.p[1];
    let p3 = data.p[2];

    let m = match J {
        SUNMatrix::Dense(m) => m,
        _ => return -1,
    };

    /* IJth(J, i, j) is 1-based in the C example → set(i-1, j-1). */
    m.set(0, 0, p1 + cj);
    m.set(1, 0, -p1);
    m.set(2, 0, ONE);

    m.set(0, 1, -p2 * y3);
    m.set(1, 1, p2 * y3 + 2.0 * p3 * y2 + cj);
    m.set(2, 1, ONE);

    m.set(0, 2, -p2 * y2);
    m.set(1, 2, p2 * y2);
    m.set(2, 2, ONE);

    0
}

/*
 * rhsQ routine. Compute fQ(t,y).
 */
fn rhsQ(_t: f64, yy: &NVector, _yp: &NVector, qdot: &mut NVector, _user_data: &mut UserData) -> i32 {
    qdot.data[0] = yy.data[2];
    0
}

/*
 * EwtSet function. Computes the error weights at the current solution.
 */
fn ewt(y: &NVector, w: &mut NVector, _user_data: &mut UserData) -> i32 {
    let rtol = RTOL;
    let atol = [ATOL1, ATOL2, ATOL3];

    for i in 0..3 {
        let yy = y.data[i];
        let ww = rtol * SUNRabs(yy) + atol[i];
        if ww <= 0.0 {
            return -1;
        }
        w.data[i] = 1.0 / ww;
    }

    0
}

/*
 * resB routine.
 */
fn resB(_tt: f64, yy: &NVector, _yp: &NVector, yyB: &NVector, ypB: &NVector, rrB: &mut NVector,
        user_dataB: &mut UserData) -> i32 {
    let data = user_dataB.as_ref().unwrap().downcast_ref::<RobData>().unwrap();

    /* The p vector */
    let p1 = data.p[0];
    let p2 = data.p[1];
    let p3 = data.p[2];

    /* The y  vector */
    let y2 = yy.data[1];
    let y3 = yy.data[2];

    /* The lambda vector */
    let l1 = yyB.data[0];
    let l2 = yyB.data[1];
    let l3 = yyB.data[2];

    /* The lambda dot vector */
    let lp1 = ypB.data[0];
    let lp2 = ypB.data[1];

    /* Temporary variables */
    let l21 = l2 - l1;

    /* Load residual. */
    rrB.data[0] = lp1 + p1 * l21 - l3;
    rrB.data[1] = lp2 - p2 * y3 * l21 - 2.0 * p3 * y2 * l2 - l3;
    rrB.data[2] = -p2 * y2 * l21 - l3 + 1.0;

    0
}

/* Jacobian for backward problem. */
#[allow(clippy::too_many_arguments)]
fn JacB(_tt: f64, cj: f64, yy: &NVector, _yp: &NVector, _yyB: &NVector, _ypB: &NVector,
        _rrB: &NVector, JB: &mut SUNMatrix, user_data: &mut UserData, _tmp1B: &mut NVector,
        _tmp2B: &mut NVector, _tmp3B: &mut NVector) -> i32 {
    let y2 = yy.data[1];
    let y3 = yy.data[2];

    let data = user_data.as_ref().unwrap().downcast_ref::<RobData>().unwrap();
    let p1 = data.p[0];
    let p2 = data.p[1];
    let p3 = data.p[2];

    let m = match JB {
        SUNMatrix::Dense(m) => m,
        _ => return -1,
    };

    m.set(0, 0, -p1 + cj);
    m.set(0, 1, p1);
    m.set(0, 2, -ONE);

    m.set(1, 0, p2 * y3);
    m.set(1, 1, -(p2 * y3 + 2.0 * p3 * y2) + cj);
    m.set(1, 2, -ONE);

    m.set(2, 0, p2 * y2);
    m.set(2, 1, -p2 * y2);
    m.set(2, 2, -ONE);

    0
}

fn rhsQB(_tt: f64, yy: &NVector, _yp: &NVector, yyB: &NVector, _ypB: &NVector,
         rrQB: &mut NVector, _user_dataB: &mut UserData) -> i32 {
    /* The y vector */
    let y1 = yy.data[0];
    let y2 = yy.data[1];
    let y3 = yy.data[2];

    /* The lambda vector */
    let l1 = yyB.data[0];
    let l2 = yyB.data[1];

    /* Temporary variables */
    let l21 = l2 - l1;

    rrQB.data[0] = y1 * l21;
    rrQB.data[1] = -y3 * y2 * l21;
    rrQB.data[2] = -y2 * y2 * l2;

    0
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * Print results after backward integration
 */
fn PrintOutput(tfinal: f64, yB: &NVector, _ypB: &NVector, qB: &NVector) {
    print!("--------------------------------------------------------\n");
    println!("tB0:        {}", fmt_e(tfinal, 12, 4));
    println!(
        "dG/dp:      {} {} {}",
        fmt_e(-qB.data[0], 12, 4),
        fmt_e(-qB.data[1], 12, 4),
        fmt_e(-qB.data[2], 12, 4)
    );
    println!(
        "lambda(t0): {} {} {}",
        fmt_e(yB.data[0], 12, 4),
        fmt_e(yB.data[1], 12, 4),
        fmt_e(yB.data[2], 12, 4)
    );
    print!("--------------------------------------------------------\n");
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
    /* Print problem description */
    print!("\nAdjoint Sensitivity Example for Chemical Kinetics\n");
    print!("-------------------------------------------------\n\n");
    print!("DAE: dy1/dt + p1*y1 - p2*y2*y3 = 0\n");
    print!("     dy2/dt - p1*y1 + p2*y2*y3 + p3*(y2)^2 = 0\n");
    print!("               y1  +  y2  +  y3 = 0\n\n");
    print!("Find dG/dp for\n");
    print!("     G = int_t0^tB0 g(t,p,y) dt\n");
    print!("     g(t,p,y) = y3\n\n\n");

    /* Create the SUNDIALS context object for this simulation */
    let sunctx = SUNContext_Create();

    /* User data structure */
    let data = RobData { p: [0.04, 1.0e4, 3.0e7] };

    /* Initialize y */
    let mut yy = N_VNew_Serial(NEQ as i64, &sunctx);
    yy.data[0] = ONE;
    yy.data[1] = ZERO;
    yy.data[2] = ZERO;

    /* Initialize yprime */
    let mut yp = N_VClone(&yy);
    yp.data[0] = -0.04;
    yp.data[1] = 0.04;
    yp.data[2] = ZERO;

    /* Initialize q */
    let mut q = N_VNew_Serial(1, &sunctx);
    q.data[0] = ZERO;

    /* Set the scalar relative and absolute tolerances reltolQ and abstolQ */
    let reltolQ = RTOL;
    let abstolQ = ATOLQ;

    /* Create and allocate IDAS memory for forward run */
    print!("Create and allocate IDAS memory for forward runs\n");

    let mut ida_mem = IDACreate(&sunctx);

    let retval = IDAInit(&mut ida_mem, res, T0, &yy, &yp);
    if check_retval(retval, "IDAInit") {
        std::process::exit(1);
    }

    let retval = IDAWFtolerances(&mut ida_mem, ewt);
    if check_retval(retval, "IDAWFtolerances") {
        std::process::exit(1);
    }

    let retval = IDASetUserData(&mut ida_mem, Some(Box::new(data.clone())));
    if check_retval(retval, "IDASetUserData") {
        std::process::exit(1);
    }

    /* Create dense SUNMatrix for use in linear solves */
    let a = SUNDenseMatrix(NEQ as i64, NEQ as i64, &sunctx);

    /* Create dense SUNLinearSolver object */
    let ls = SUNLinSol_Dense(&yy, &a, &sunctx);

    /* Attach the matrix and linear solver */
    let retval = IDASetLinearSolver(&mut ida_mem, ls, Some(a));
    if check_retval(retval, "IDASetLinearSolver") {
        std::process::exit(1);
    }

    /* Set the user-supplied Jacobian routine */
    let retval = IDASetJacFn(&mut ida_mem, Some(Jac));
    if check_retval(retval, "IDASetJacFn") {
        std::process::exit(1);
    }

    /* Setup quadrature integration */
    let retval = IDAQuadInit(&mut ida_mem, rhsQ, &q);
    if check_retval(retval, "IDAQuadInit") {
        std::process::exit(1);
    }

    let retval = IDAQuadSStolerances(&mut ida_mem, reltolQ, abstolQ);
    if check_retval(retval, "IDAQuadSStolerances") {
        std::process::exit(1);
    }

    let retval = IDASetQuadErrCon(&mut ida_mem, true);
    if check_retval(retval, "IDASetQuadErrCon") {
        std::process::exit(1);
    }

    /* Call IDASetMaxNumSteps to set the maximum number of steps the
     * solver will take in an attempt to reach the next output time
     * during forward integration. */
    let retval = IDASetMaxNumSteps(&mut ida_mem, 2500);
    if check_retval(retval, "IDASetMaxNumSteps") {
        std::process::exit(1);
    }

    /* Allocate global memory */

    let steps = STEPS;
    let retval = IDAAdjInit(&mut ida_mem, steps, IDA_HERMITE);
    /*let retval = IDAAdjInit(&mut ida_mem, steps, IDA_POLYNOMIAL);*/
    if check_retval(retval, "IDAAdjInit") {
        std::process::exit(1);
    }

    /* Perform forward run */
    print!("Forward integration ...\n");

    /* Integrate till TB1 and get the solution (y, y') at that time. */
    let mut time = 0.0;
    let mut ncheck = 0;
    let retval = IDASolveF(&mut ida_mem, TB1, &mut time, &mut yy, &mut yp, IDA_NORMAL, &mut ncheck);
    if check_retval(retval, "IDASolveF") {
        std::process::exit(1);
    }

    let mut yyTB1 = N_VClone(&yy);
    let mut ypTB1 = N_VClone(&yp);
    /* Save the states at t=TB1. */
    N_VScale(ONE, &yy, &mut yyTB1);
    N_VScale(ONE, &yp, &mut ypTB1);

    /* Continue integrating till TOUT is reached. */
    let retval = IDASolveF(&mut ida_mem, TOUT, &mut time, &mut yy, &mut yp, IDA_NORMAL, &mut ncheck);
    if check_retval(retval, "IDASolveF") {
        std::process::exit(1);
    }

    let retval = IDAGetQuad(&ida_mem, &mut time, &mut q);
    if check_retval(retval, "IDAGetQuad") {
        std::process::exit(1);
    }

    print!("--------------------------------------------------------\n");
    println!("G:          {} ", fmt_e(q.data[0], 12, 4));
    print!("--------------------------------------------------------\n");

    /* Print final statistics to the screen */
    print!("\nFinal Statistics:\n");
    let mut stdout = std::io::stdout();
    IDAPrintAllStats(&mut ida_mem, &mut stdout, SUN_OUTPUTFORMAT_TABLE);

    /* Print final statistics to a file in CSV format */
    let mut fid = std::fs::File::create("idasRoberts_ASAi_dns_fwd_stats.csv").expect("create csv");
    IDAPrintAllStats(&mut ida_mem, &mut fid, SUN_OUTPUTFORMAT_CSV);
    drop(fid);

    /* Test check point linked list
       (uncomment next block to print check point information) */

    /*
    {
        println!("\nList of Check Points (ncheck = {})\n", ncheck);
        let mut ckpnt = vec![IDAadjCheckPointRec::default(); (ncheck + 1) as usize];
        IDAGetAdjCheckPointsInfo(&mut ida_mem, &mut ckpnt);
        for ck in ckpnt.iter() {
            println!("Address:       {:?}", ck.my_addr);
            println!("Next:          {:?}", ck.next_addr);
            println!("Time interval: {:e}  {:e}", ck.t0, ck.t1);
            println!("Step number:   {}", ck.nstep);
            println!("Order:         {}", ck.order);
            println!("Step size:     {:e}", ck.step);
            println!();
        }
    }
    */

    /* Create BACKWARD problem. */

    /* Allocate yB (i.e. lambda_0). */
    let mut yB = N_VClone(&yy);

    /* Consistently initialize yB. */
    yB.data[0] = ZERO;
    yB.data[1] = ZERO;
    yB.data[2] = ONE;

    /* Allocate ypB (i.e. lambda'_0). */
    let mut ypB = N_VClone(&yy);

    /* Consistently initialize ypB. */
    ypB.data[0] = ONE;
    ypB.data[1] = ONE;
    ypB.data[2] = ZERO;

    /* Set the scalar relative tolerance reltolB */
    let reltolB = RTOL;

    /* Set the scalar absolute tolerance abstolB */
    let abstolB = ATOLA;

    /* Set the scalar absolute tolerance abstolQB */
    let abstolQB = ATOLQ;

    /* Create and allocate IDAS memory for backward run */
    print!("\nCreate and allocate IDAS memory for backward run\n");

    let mut indexB = 0;
    let retval = IDACreateB(&mut ida_mem, &mut indexB);
    if check_retval(retval, "IDACreateB") {
        std::process::exit(1);
    }

    let retval = IDAInitB(&mut ida_mem, indexB, resB, TB2, &yB, &ypB);
    if check_retval(retval, "IDAInitB") {
        std::process::exit(1);
    }

    let retval = IDASStolerancesB(&mut ida_mem, indexB, reltolB, abstolB);
    if check_retval(retval, "IDASStolerancesB") {
        std::process::exit(1);
    }

    let retval = IDASetUserDataB(&mut ida_mem, indexB, Some(Box::new(data.clone())));
    if check_retval(retval, "IDASetUserDataB") {
        std::process::exit(1);
    }

    let retval = IDASetMaxNumStepsB(&mut ida_mem, indexB, 1000);
    if check_retval(retval, "IDASetMaxNumStepsB") {
        std::process::exit(1);
    }

    /* Create dense SUNMatrix for use in linear solves */
    let ab = SUNDenseMatrix(NEQ as i64, NEQ as i64, &sunctx);

    /* Create dense SUNLinearSolver object */
    let lsb = SUNLinSol_Dense(&yB, &ab, &sunctx);

    /* Attach the matrix and linear solver */
    let retval = IDASetLinearSolverB(&mut ida_mem, indexB, lsb, Some(ab));
    if check_retval(retval, "IDASetLinearSolverB") {
        std::process::exit(1);
    }

    /* Set the user-supplied Jacobian routine */
    let retval = IDASetJacFnB(&mut ida_mem, indexB, Some(JacB));
    if check_retval(retval, "IDASetJacFnB") {
        std::process::exit(1);
    }

    /* Quadrature for backward problem. */

    /* Initialize qB */
    let mut qB = N_VNew_Serial(NP as i64, &sunctx);
    qB.data[0] = ZERO;
    qB.data[1] = ZERO;
    qB.data[2] = ZERO;

    let retval = IDAQuadInitB(&mut ida_mem, indexB, rhsQB, &qB);
    if check_retval(retval, "IDAQuadInitB") {
        std::process::exit(1);
    }

    let retval = IDAQuadSStolerancesB(&mut ida_mem, indexB, reltolB, abstolQB);
    if check_retval(retval, "IDAQuadSStolerancesB") {
        std::process::exit(1);
    }

    /* Include quadratures in error control. */
    let retval = IDASetQuadErrConB(&mut ida_mem, indexB, true);
    if check_retval(retval, "IDASetQuadErrConB") {
        std::process::exit(1);
    }

    /* Backward Integration */
    print!("Backward integration ...\n");

    let retval = IDASolveB(&mut ida_mem, T0, IDA_NORMAL);
    if check_retval(retval, "IDASolveB") {
        std::process::exit(1);
    }

    let retval = IDAGetB(&mut ida_mem, indexB, &mut time, &mut yB, &mut ypB);
    if check_retval(retval, "IDAGetB") {
        std::process::exit(1);
    }

    let retval = IDAGetQuadB(&mut ida_mem, indexB, &mut time, &mut qB);
    if check_retval(retval, "IDAGetB") {
        std::process::exit(1);
    }

    PrintOutput(TB2, &yB, &ypB, &qB);

    /* Print final statistics to the screen */
    print!("\nFinal Statistics:\n");
    let mut stdout = std::io::stdout();
    IDAPrintAllStats(IDAGetAdjIDABmem(&mut ida_mem, indexB).unwrap(), &mut stdout,
                     SUN_OUTPUTFORMAT_TABLE);

    /* Print final statistics to a file in CSV format */
    let mut fid = std::fs::File::create("idasRoberts_ASAi_dns_bkw1_stats.csv").expect("create csv");
    IDAPrintAllStats(IDAGetAdjIDABmem(&mut ida_mem, indexB).unwrap(), &mut fid,
                     SUN_OUTPUTFORMAT_CSV);
    drop(fid);

    /* Reinitialize backward phase and start from a different time (TB1). */
    print!("\nRe-initialize IDAS memory for backward run\n");

    /* Both algebraic part from y and the entire y' are computed by IDACalcIC. */
    yB.data[0] = ZERO;
    yB.data[1] = ZERO;
    yB.data[2] = 0.50; /* not consistent */

    /* Rough guess for ypB. */
    ypB.data[0] = 0.80;
    ypB.data[1] = 0.75;
    ypB.data[2] = ZERO;

    /* Initialize qB */
    qB.data[0] = ZERO;
    qB.data[1] = ZERO;
    qB.data[2] = ZERO;

    let retval = IDAReInitB(&mut ida_mem, indexB, TB1, &yB, &ypB);
    if check_retval(retval, "IDAReInitB") {
        std::process::exit(1);
    }

    /* Also reinitialize quadratures. */
    let retval = IDAQuadReInitB(&mut ida_mem, indexB, &qB);
    if check_retval(retval, "IDAQuadReInitB") {
        std::process::exit(1);
    }

    /* Use IDACalcICB to compute consistent initial conditions
       for this backward problem. */

    let mut id = N_VClone(&yy);
    id.data[0] = 1.0;
    id.data[1] = 1.0;
    id.data[2] = 0.0;

    /* Specify which variables are differential (1) and which algebraic (0).*/
    let retval = IDASetIdB(&mut ida_mem, indexB, Some(&id));
    if check_retval(retval, "IDASetId") {
        std::process::exit(1);
    }

    let retval = IDACalcICB(&mut ida_mem, indexB, T1B, &yyTB1, &ypTB1);
    if check_retval(retval, "IDACalcICB") {
        std::process::exit(1);
    }

    /* Get the consistent IC found by IDAS. */
    let retval = IDAGetConsistentICB(&mut ida_mem, indexB, Some(&mut yB), Some(&mut ypB));
    if check_retval(retval, "IDAGetConsistentICB") {
        std::process::exit(1);
    }

    let retval = IDASolveB(&mut ida_mem, T0, IDA_NORMAL);
    if check_retval(retval, "IDASolveB") {
        std::process::exit(1);
    }

    let retval = IDAGetB(&mut ida_mem, indexB, &mut time, &mut yB, &mut ypB);
    if check_retval(retval, "IDAGetB") {
        std::process::exit(1);
    }

    let retval = IDAGetQuadB(&mut ida_mem, indexB, &mut time, &mut qB);
    if check_retval(retval, "IDAGetQuadB") {
        std::process::exit(1);
    }

    PrintOutput(TB1, &yB, &ypB, &qB);

    /* Print final statistics to the screen */
    print!("\nFinal Statistics:\n");
    let mut stdout = std::io::stdout();
    IDAPrintAllStats(IDAGetAdjIDABmem(&mut ida_mem, indexB).unwrap(), &mut stdout,
                     SUN_OUTPUTFORMAT_TABLE);

    /* Print final statistics to a file in CSV format */
    let mut fid = std::fs::File::create("idasRoberts_ASAi_dns_bkw1_stats.csv").expect("create csv");
    IDAPrintAllStats(IDAGetAdjIDABmem(&mut ida_mem, indexB).unwrap(), &mut fid,
                     SUN_OUTPUTFORMAT_CSV);
    drop(fid);

    /* Free memory (RAII) */
    IDAFree(ida_mem);
}
