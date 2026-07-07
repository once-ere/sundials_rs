/* -----------------------------------------------------------------
 * Translated from examples/ida/serial/idaRoberts_dns.c (IDA 7.7.0)
 * Programmer(s): Allan Taylor, Alan Hindmarsh and Radu Serban @ LLNL
 *
 * This simple example problem for IDA, due to Robertson, is from
 * chemical kinetics, and consists of the following three equations:
 *
 *      dy1/dt = -.04*y1 + 1.e4*y2*y3
 *      dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
 *         0   = y1 + y2 + y3 - 1
 *
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1, y2 = y3 = 0.
 *
 * While integrating the system, we also use the rootfinding feature
 * to find the points at which y1 = 1e-4 or at which y3 = 0.01.
 *
 * The problem is solved with IDA using the DENSE linear solver, with
 * a user-supplied Jacobian. Output is printed at t = .4, 4, 40, ...,
 * 4e10.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use ida_rs::sundials_utils::{fmt_e, fmt_g};
use ida_rs::*;

/* Problem Constants */
const NEQ: usize = 3;
const NOUT: i32 = 12;

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/*
 * Define the system residual function.
 */
fn resrob(_tres: f64, yy: &NVector, yp: &NVector, rr: &mut NVector, _ud: &mut UserData) -> i32 {
    let yval = &yy.data;
    let ypval = &yp.data;

    let mut rval0 = -0.04 * yval[0] + 1.0e4 * yval[1] * yval[2];
    let rval1 = -rval0 - 3.0e7 * yval[1] * yval[1] - ypval[1];
    rval0 -= ypval[0];
    let rval2 = yval[0] + yval[1] + yval[2] - ONE;

    rr.data[0] = rval0;
    rr.data[1] = rval1;
    rr.data[2] = rval2;

    0
}

/*
 * Root function routine. Compute functions g_i(t,y) for i = 0,1.
 */
fn grob(_t: f64, yy: &NVector, _yp: &NVector, gout: &mut [f64], _ud: &mut UserData) -> i32 {
    let yval = &yy.data;
    let y1 = yval[0];
    let y3 = yval[2];
    gout[0] = y1 - 0.0001;
    gout[1] = y3 - 0.01;

    0
}

/*
 * Define the Jacobian function.
 */
fn jacrob(
    _tt: f64,
    cj: f64,
    yy: &NVector,
    _yp: &NVector,
    _rr: &NVector,
    jj: &mut SUNMatrix,
    _ud: &mut UserData,
    _t1: &mut NVector,
    _t2: &mut NVector,
    _t3: &mut NVector,
) -> i32 {
    let yval = &yy.data;

    let m = match jj {
        SUNMatrix::Dense(m) => m,
        _ => return -1,
    };

    /* IJth(JJ, i, j) is 1-based in the C example → set(i-1, j-1). */
    m.set(0, 0, -0.04 - cj);
    m.set(1, 0, 0.04);
    m.set(2, 0, ONE);
    m.set(0, 1, 1.0e4 * yval[2]);
    m.set(1, 1, -1.0e4 * yval[2] - 6.0e7 * yval[1] - cj);
    m.set(2, 1, ONE);
    m.set(0, 2, 1.0e4 * yval[1]);
    m.set(1, 2, -1.0e4 * yval[1]);
    m.set(2, 2, ONE);

    0
}

/*
 *--------------------------------------------------------------------
 * Private functions
 *--------------------------------------------------------------------
 */

fn PrintHeader(rtol: f64, avtol: &NVector, y: &NVector) {
    let atval = &avtol.data;
    let yval = &y.data;

    print!("\nidaRoberts_dns: Robertson kinetics DAE serial example problem for IDA\n");
    print!("         Three equation chemical kinetics problem.\n\n");
    print!("Linear solver: DENSE, with user-supplied Jacobian.\n");
    print!(
        "Tolerance parameters:  rtol = {}   atol = {} {} {} \n",
        fmt_g(rtol, 0, 6),
        fmt_g(atval[0], 0, 6),
        fmt_g(atval[1], 0, 6),
        fmt_g(atval[2], 0, 6)
    );
    print!(
        "Initial conditions y0 = ({} {} {})\n",
        fmt_g(yval[0], 0, 6),
        fmt_g(yval[1], 0, 6),
        fmt_g(yval[2], 0, 6)
    );
    print!("Constraints and id not used.\n\n");
    print!("-----------------------------------------------------------------------\n");
    print!("  t             y1           y2           y3");
    print!("      | nst  k      h\n");
    print!("-----------------------------------------------------------------------\n");
}

fn PrintOutput(mem: &mut IDAMem, t: f64, y: &NVector) {
    let yval = &y.data;

    let mut kused = 0i32;
    let mut nst = 0i64;
    let mut hused = 0.0f64;

    IDAGetLastOrder(mem, &mut kused);
    IDAGetNumSteps(mem, &mut nst);
    IDAGetLastStep(mem, &mut hused);

    println!(
        "{} {} {} {} | {:3}  {:1} {}",
        fmt_e(t, 10, 4),
        fmt_e(yval[0], 12, 4),
        fmt_e(yval[1], 12, 4),
        fmt_e(yval[2], 12, 4),
        nst,
        kused,
        fmt_e(hused, 12, 4)
    );
}

fn PrintRootInfo(root_f1: i32, root_f2: i32) {
    println!("    rootsfound[] = {:3} {:3}", root_f1, root_f2);
}

fn check_retval(retval: i32, funcname: &str) -> bool {
    if retval < 0 {
        eprintln!("\nSUNDIALS_ERROR: {}() failed with retval = {}\n", funcname, retval);
        return true;
    }
    false
}

/* compare the solution at the final time 4e10s to a reference solution computed
   using a relative tolerance of 1e-8 and absolute tolerance of 1e-14 */
fn check_ans(y: &NVector, _t: f64, rtol: f64, atol: &NVector) -> i32 {
    let mut refv = N_VClone(y);
    let mut ewt = N_VClone(y);

    /* set the reference solution data */
    refv.data[0] = 5.2083474251394888e-08;
    refv.data[1] = 2.0833390772616859e-13;
    refv.data[2] = 9.9999994791631752e-01;

    /* compute the error weight vector, loosen atol */
    N_VAbs(&refv, &mut ewt);
    ewt.linear_sum_with(rtol, 10.0, atol);
    if N_VMin(&ewt) <= ZERO {
        eprintln!("\nSUNDIALS_ERROR: check_ans failed - ewt <= 0\n");
        return -1;
    }
    ewt.invert_inplace();

    /* compute the solution error: ref = y - ref */
    refv.linear_sum_with(-1.0, 1.0, y);
    let err = N_VWrmsNorm(&refv, &ewt);

    /* is the solution within the tolerances? */
    let passfail = if err < ONE { 0 } else { 1 };
    if passfail != 0 {
        print!("\nSUNDIALS_WARNING: check_ans error={}\n\n", fmt_g(err, 0, 6));
    }
    passfail
}

fn main() {
    /* Create SUNDIALS context */
    let sunctx = SUNContext_Create();

    /* Allocate N-vectors. */
    let mut yy = N_VNew_Serial(NEQ as i64, &sunctx);
    let mut yp = N_VClone(&yy);
    let mut avtol = N_VClone(&yy);

    /* Create and initialize y, y', and absolute tolerance vectors. */
    yy.data[0] = ONE;
    yy.data[1] = ZERO;
    yy.data[2] = ZERO;

    yp.data[0] = -0.04;
    yp.data[1] = 0.04;
    yp.data[2] = ZERO;

    let rtol = 1.0e-4;

    avtol.data[0] = 1.0e-8;
    avtol.data[1] = 1.0e-6;
    avtol.data[2] = 1.0e-6;

    /* Integration limits */
    let t0 = ZERO;
    let tout1 = 0.4;

    PrintHeader(rtol, &avtol, &yy);

    /* Call IDACreate and IDAInit to initialize IDA memory */
    let mut mem = IDACreate(&sunctx);
    let mut retval = IDAInit(&mut mem, resrob, t0, &yy, &yp);
    if check_retval(retval, "IDAInit") {
        std::process::exit(1);
    }
    /* Call IDASVtolerances to set tolerances */
    retval = IDASVtolerances(&mut mem, rtol, &avtol);
    if check_retval(retval, "IDASVtolerances") {
        std::process::exit(1);
    }

    /* Call IDARootInit to specify the root function grob with 2 components */
    retval = IDARootInit(&mut mem, 2, Some(grob));
    if check_retval(retval, "IDARootInit") {
        std::process::exit(1);
    }

    /* Create dense SUNMatrix for use in linear solves */
    let a = SUNDenseMatrix(NEQ as i64, NEQ as i64, &sunctx);

    /* Create dense SUNLinearSolver object */
    let ls = SUNLinSol_Dense(&yy, &a, &sunctx);

    /* Attach the matrix and linear solver */
    retval = IDASetLinearSolver(&mut mem, ls, Some(a));
    if check_retval(retval, "IDASetLinearSolver") {
        std::process::exit(1);
    }

    /* Set the user-supplied Jacobian routine */
    retval = IDASetJacFn(&mut mem, Some(jacrob));
    if check_retval(retval, "IDASetJacFn") {
        std::process::exit(1);
    }

    /* Create Newton SUNNonlinearSolver object. IDA uses a Newton
     * SUNNonlinearSolver by default, so it is unnecessary to create it
     * and attach it. It is done in this example code solely for
     * demonstration purposes. */
    let nls = SUNNonlinSol_Newton(&yy, &sunctx);

    /* Attach the nonlinear solver */
    retval = IDASetNonlinearSolver(&mut mem, nls);
    if check_retval(retval, "IDASetNonlinearSolver") {
        std::process::exit(1);
    }

    /* In loop, call IDASolve, print results, and test for error.
       Break out of loop when NOUT preset output times have been reached. */
    let mut iout = 0;
    let mut tout = tout1;
    let mut tret = 0.0;
    let mut rootsfound = [0i32; 2];
    loop {
        retval = IDASolve(&mut mem, tout, &mut tret, &mut yy, &mut yp, IDA_NORMAL);

        PrintOutput(&mut mem, tret, &yy);

        if check_retval(retval, "IDASolve") {
            std::process::exit(1);
        }

        if retval == IDA_ROOT_RETURN {
            let retvalr = IDAGetRootInfo(&mut mem, &mut rootsfound);
            check_retval(retvalr, "IDAGetRootInfo");
            PrintRootInfo(rootsfound[0], rootsfound[1]);
        }

        if retval == IDA_SUCCESS {
            iout += 1;
            tout *= 10.0;
        }

        if iout == NOUT {
            break;
        }
    }

    /* Print final statistics to the screen */
    println!("\nFinal Statistics:");
    let mut stdout = std::io::stdout();
    IDAPrintAllStats(&mut mem, &mut stdout, SUN_OUTPUTFORMAT_TABLE);

    /* Print final statistics to a file in CSV format */
    let mut fid = std::fs::File::create("idaRoberts_dns_stats.csv").expect("create csv");
    IDAPrintAllStats(&mut mem, &mut fid, SUN_OUTPUTFORMAT_CSV);
    drop(fid);

    /* check the solution error */
    let retval = check_ans(&yy, tret, rtol, &avtol);

    /* Free memory (RAII) */
    IDAFree(mem);

    std::process::exit(retval);
}
