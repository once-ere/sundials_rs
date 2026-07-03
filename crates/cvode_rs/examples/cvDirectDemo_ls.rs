/* -----------------------------------------------------------------
 * Translated from examples/cvode/serial/cvDirectDemo_ls.c (CVODE 7.7.0)
 *
 * Demonstration program for CVODE - direct linear solvers.
 * Two separate problems are solved using both the CV_ADAMS and CV_BDF
 * linear multistep methods in combination with the
 * SUNNONLINSOL_FIXEDPOINT and SUNNONLINSOL_NEWTON nonlinear solver
 * modules:
 *
 * Problem 1: Van der Pol oscillator
 *   xdotdot - 3*(1 - x^2)*xdot + x = 0, x(0) = 2, xdot(0) = 0.
 * This second-order ODE is converted to a first-order system by
 * defining y0 = x and y1 = xdot.
 * The NEWTON iteration cases use the following types of Jacobian
 * approximation: (1) dense, user-supplied, (2) dense, difference
 * quotient approximation, (3) diagonal approximation.
 *
 * Problem 2: ydot = A * y, where A is a banded lower triangular
 * matrix derived from 2-D advection PDE.
 * The NEWTON iteration cases use the following types of Jacobian
 * approximation: (1) band, user-supplied, (2) band, difference
 * quotient approximation, (3) diagonal approximation.
 *
 * For each problem, in the series of eight runs, CVodeInit is
 * called only once, for the first run, whereas CVodeReInit is
 * called for each of the remaining seven runs.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use cvode_rs::sundials_utils::{fmt_e, fmt_f, fmt_g};
use cvode_rs::*;

/* Shared Problem Constants */

const ATOL: f64 = 1.0e-6;
const RTOL: f64 = 0.0;

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;
const THIRTY: f64 = 30.0;

/* Problem #1 Constants */

const P1_NEQ: i64 = 2;
const P1_ETA: f64 = 3.0;
const P1_NOUT: i32 = 4;
const P1_T0: f64 = 0.0;
const P1_T1: f64 = 1.39283880203;
const P1_DTOUT: f64 = 2.214773875;
const P1_TOL_FACTOR: f64 = 1.0e4;

/* Problem #2 Constants */

const P2_MESHX: usize = 5;
const P2_MESHY: usize = 5;
const P2_NEQ: i64 = (P2_MESHX * P2_MESHY) as i64;
const P2_ALPH1: f64 = 1.0;
const P2_ALPH2: f64 = 1.0;
const P2_NOUT: i32 = 5;
const P2_ML: i64 = 5;
const P2_MU: i64 = 0;
const P2_T0: f64 = 0.0;
const P2_T1: f64 = 0.01;
const P2_TOUT_MULT: f64 = 10.0;
const P2_TOL_FACTOR: f64 = 1.0e3;

/* Linear Solver Options */

const FUNC: i32 = 0;
const DENSE_USER: i32 = 1;
const DENSE_DQ: i32 = 2;
const DIAG: i32 = 3;
const BAND_USER: i32 = 4;
const BAND_DQ: i32 = 5;

/* Implementation */

fn main() {
    let mut nerr = Problem1();
    nerr += Problem2();
    PrintErrInfo(nerr);
}

fn Problem1() -> i32 {
    let reltol = RTOL;
    let abstol = ATOL;
    let mut nerr: i32 = 0;

    /* Create the SUNDIALS context */
    let sunctx = SUNContext_Create();

    let mut y = N_VNew_Serial(P1_NEQ, &sunctx);
    PrintIntro1();

    let mut cvode_mem = CVodeCreate(CV_ADAMS, &sunctx);

    for miter in FUNC..=DIAG {
        let mut ero = ZERO;
        y.data[0] = TWO;
        y.data[1] = ZERO;

        let firstrun = miter == FUNC;
        if firstrun {
            /* initialize CVode */
            let retval = CVodeInit(&mut cvode_mem, f1, P1_T0, &y);
            if check_retval(retval, "CVodeInit") {
                std::process::exit(1);
            }

            /* set scalar tolerances */
            let retval = CVodeSStolerances(&mut cvode_mem, reltol, abstol);
            if check_retval(retval, "CVodeSStolerances") {
                std::process::exit(1);
            }
        } else {
            /* reinitialize CVode */
            let retval = CVodeReInit(&mut cvode_mem, P1_T0, &y);
            if check_retval(retval, "CVodeReInit") {
                std::process::exit(1);
            }
        }

        let retval = PrepareNextRun(&sunctx, &mut cvode_mem, CV_ADAMS, miter, &y, 0, 0);
        if check_retval(retval, "PrepareNextRun") {
            std::process::exit(1);
        }

        PrintHeader1();

        let mut t = ZERO;
        let mut qu: i32 = 0;
        let mut hu = ZERO;
        let mut tout = P1_T1;
        for iout in 1..=P1_NOUT {
            let retval = CVode(&mut cvode_mem, tout, &mut y, &mut t, CV_NORMAL);
            check_retval(retval, "CVode");
            let temp_retval = CVodeGetLastOrder(&mut cvode_mem, &mut qu);
            if check_retval(temp_retval, "CVodeGetLastOrder") {
                nerr += 1;
            }
            let temp_retval = CVodeGetLastStep(&mut cvode_mem, &mut hu);
            if check_retval(temp_retval, "CVodeGetLastStep") {
                nerr += 1;
            }
            PrintOutput1(t, y.data[0], y.data[1], qu, hu);
            if retval != CV_SUCCESS {
                nerr += 1;
                break;
            }
            if iout % 2 == 0 {
                let er = y.data[0].abs() / abstol;
                if er > ero {
                    ero = er;
                }
                if er > P1_TOL_FACTOR {
                    nerr += 1;
                    PrintErrOutput(P1_TOL_FACTOR);
                }
            }
            tout += P1_DTOUT;
        }

        PrintFinalStats(&mut cvode_mem, miter, ero);
    }

    CVodeFree(cvode_mem);

    cvode_mem = CVodeCreate(CV_BDF, &sunctx);

    for miter in FUNC..=DIAG {
        let mut ero = ZERO;
        y.data[0] = TWO;
        y.data[1] = ZERO;

        let firstrun = miter == FUNC;
        if firstrun {
            /* initialize CVode */
            let retval = CVodeInit(&mut cvode_mem, f1, P1_T0, &y);
            if check_retval(retval, "CVodeInit") {
                std::process::exit(1);
            }

            /* set scalar tolerances */
            let retval = CVodeSStolerances(&mut cvode_mem, reltol, abstol);
            if check_retval(retval, "CVodeSStolerances") {
                std::process::exit(1);
            }
        } else {
            /* reinitialize CVode */
            let retval = CVodeReInit(&mut cvode_mem, P1_T0, &y);
            if check_retval(retval, "CVodeReInit") {
                std::process::exit(1);
            }
        }

        let retval = PrepareNextRun(&sunctx, &mut cvode_mem, CV_BDF, miter, &y, 0, 0);
        if check_retval(retval, "PrepareNextRun") {
            std::process::exit(1);
        }

        PrintHeader1();

        let mut t = ZERO;
        let mut qu: i32 = 0;
        let mut hu = ZERO;
        let mut tout = P1_T1;
        for iout in 1..=P1_NOUT {
            let retval = CVode(&mut cvode_mem, tout, &mut y, &mut t, CV_NORMAL);
            check_retval(retval, "CVode");
            let temp_retval = CVodeGetLastOrder(&mut cvode_mem, &mut qu);
            if check_retval(temp_retval, "CVodeGetLastOrder") {
                nerr += 1;
            }
            let temp_retval = CVodeGetLastStep(&mut cvode_mem, &mut hu);
            if check_retval(temp_retval, "CVodeGetLastStep") {
                nerr += 1;
            }
            PrintOutput1(t, y.data[0], y.data[1], qu, hu);
            if retval != CV_SUCCESS {
                nerr += 1;
                break;
            }
            if iout % 2 == 0 {
                let er = y.data[0].abs() / abstol;
                if er > ero {
                    ero = er;
                }
                if er > P1_TOL_FACTOR {
                    nerr += 1;
                    PrintErrOutput(P1_TOL_FACTOR);
                }
            }
            tout += P1_DTOUT;
        }

        PrintFinalStats(&mut cvode_mem, miter, ero);
    }

    CVodeFree(cvode_mem);

    nerr
}

fn PrintIntro1() {
    println!("Demonstration program for CVODE package - direct linear solvers");
    print!("\n\n");
    println!("Problem 1: Van der Pol oscillator");
    println!(" xdotdot - 3*(1 - x^2)*xdot + x = 0, x(0) = 2, xdot(0) = 0");
    print!(
        " neq = {},  reltol = {},  abstol = {}",
        P1_NEQ,
        fmt_g(RTOL, 0, 2),
        fmt_g(ATOL, 0, 2)
    );
}

fn PrintHeader1() {
    println!("\n     t           x              xdot         qu     hu ");
}

fn PrintOutput1(t: f64, y0: f64, y1: f64, qu: i32, hu: f64) {
    println!(
        "{}    {}   {}   {:2}    {}",
        fmt_f(t, 10, 5),
        fmt_e(y0, 12, 5),
        fmt_e(y1, 12, 5),
        qu,
        fmt_e(hu, 6, 4)
    );
}

fn f1(_t: f64, y: &NVector, ydot: &mut NVector, _user_data: &mut UserData) -> i32 {
    let y0 = y.data[0];
    let y1 = y.data[1];

    ydot.data[0] = y1;
    ydot.data[1] = (ONE - y0 * y0) * P1_ETA * y1 - y0;

    0
}

fn Jac1(
    _tn: f64,
    y: &NVector,
    _fy: &NVector,
    j: &mut SUNMatrix,
    _user_data: &mut UserData,
    _tmp1: &mut NVector,
    _tmp2: &mut NVector,
    _tmp3: &mut NVector,
) -> i32 {
    let y0 = y.data[0];
    let y1 = y.data[1];

    let jm = match j {
        SUNMatrix::Dense(m) => m,
        _ => return -1,
    };
    jm.set(0, 1, ONE);
    jm.set(1, 0, -TWO * P1_ETA * y0 * y1 - ONE);
    jm.set(1, 1, P1_ETA * (ONE - y0 * y0));

    0
}

fn Problem2() -> i32 {
    let reltol = RTOL;
    let abstol = ATOL;
    let mut nerr: i32 = 0;

    /* Create SUNDIALS context */
    let sunctx = SUNContext_Create();

    let mut y = N_VNew_Serial(P2_NEQ, &sunctx);

    PrintIntro2();

    let mut cvode_mem = CVodeCreate(CV_ADAMS, &sunctx);

    for miter in FUNC..=BAND_DQ {
        if miter == DENSE_USER || miter == DENSE_DQ {
            continue;
        }
        let mut ero = ZERO;
        N_VConst(ZERO, &mut y);
        y.data[0] = ONE;

        let firstrun = miter == FUNC;
        if firstrun {
            /* initialize CVode */
            let retval = CVodeInit(&mut cvode_mem, f2, P2_T0, &y);
            if check_retval(retval, "CVodeInit") {
                std::process::exit(1);
            }

            /* set scalar tolerances */
            let retval = CVodeSStolerances(&mut cvode_mem, reltol, abstol);
            if check_retval(retval, "CVodeSStolerances") {
                std::process::exit(1);
            }
        } else {
            /* reinitialize CVode */
            let retval = CVodeReInit(&mut cvode_mem, P2_T0, &y);
            if check_retval(retval, "CVodeReInit") {
                std::process::exit(1);
            }
        }

        let retval = PrepareNextRun(&sunctx, &mut cvode_mem, CV_ADAMS, miter, &y, P2_MU, P2_ML);
        if check_retval(retval, "PrepareNextRun") {
            std::process::exit(1);
        }

        PrintHeader2();

        let mut t = ZERO;
        let mut qu: i32 = 0;
        let mut hu = ZERO;
        let mut tout = P2_T1;
        for _iout in 1..=P2_NOUT {
            let retval = CVode(&mut cvode_mem, tout, &mut y, &mut t, CV_NORMAL);
            check_retval(retval, "CVode");
            let erm = MaxError(&y, t);
            let temp_retval = CVodeGetLastOrder(&mut cvode_mem, &mut qu);
            if check_retval(temp_retval, "CVodeGetLastOrder") {
                nerr += 1;
            }
            let temp_retval = CVodeGetLastStep(&mut cvode_mem, &mut hu);
            if check_retval(temp_retval, "CVodeGetLastStep") {
                nerr += 1;
            }
            PrintOutput2(t, erm, qu, hu);
            if retval != CV_SUCCESS {
                nerr += 1;
                break;
            }
            let er = erm / abstol;
            if er > ero {
                ero = er;
            }
            if er > P2_TOL_FACTOR {
                nerr += 1;
                PrintErrOutput(P2_TOL_FACTOR);
            }
            tout *= P2_TOUT_MULT;
        }

        PrintFinalStats(&mut cvode_mem, miter, ero);
    }

    CVodeFree(cvode_mem);

    cvode_mem = CVodeCreate(CV_BDF, &sunctx);

    for miter in FUNC..=BAND_DQ {
        if miter == DENSE_USER || miter == DENSE_DQ {
            continue;
        }
        let mut ero = ZERO;
        N_VConst(ZERO, &mut y);
        y.data[0] = ONE;

        let firstrun = miter == FUNC;
        if firstrun {
            /* initialize CVode */
            let retval = CVodeInit(&mut cvode_mem, f2, P2_T0, &y);
            if check_retval(retval, "CVodeInit") {
                std::process::exit(1);
            }

            /* set scalar tolerances */
            let retval = CVodeSStolerances(&mut cvode_mem, reltol, abstol);
            if check_retval(retval, "CVodeSStolerances") {
                std::process::exit(1);
            }
        } else {
            /* reinitialize CVode */
            let retval = CVodeReInit(&mut cvode_mem, P2_T0, &y);
            if check_retval(retval, "CVodeReInit") {
                std::process::exit(1);
            }
        }

        let retval = PrepareNextRun(&sunctx, &mut cvode_mem, CV_BDF, miter, &y, P2_MU, P2_ML);
        if check_retval(retval, "PrepareNextRun") {
            std::process::exit(1);
        }

        PrintHeader2();

        let mut t = ZERO;
        let mut qu: i32 = 0;
        let mut hu = ZERO;
        let mut tout = P2_T1;
        for _iout in 1..=P2_NOUT {
            let retval = CVode(&mut cvode_mem, tout, &mut y, &mut t, CV_NORMAL);
            check_retval(retval, "CVode");
            let erm = MaxError(&y, t);
            let temp_retval = CVodeGetLastOrder(&mut cvode_mem, &mut qu);
            if check_retval(temp_retval, "CVodeGetLastOrder") {
                nerr += 1;
            }
            let temp_retval = CVodeGetLastStep(&mut cvode_mem, &mut hu);
            if check_retval(temp_retval, "CVodeGetLastStep") {
                nerr += 1;
            }
            PrintOutput2(t, erm, qu, hu);
            if retval != CV_SUCCESS {
                nerr += 1;
                break;
            }
            let er = erm / abstol;
            if er > ero {
                ero = er;
            }
            if er > P2_TOL_FACTOR {
                nerr += 1;
                PrintErrOutput(P2_TOL_FACTOR);
            }
            tout *= P2_TOUT_MULT;
        }

        PrintFinalStats(&mut cvode_mem, miter, ero);
    }

    CVodeFree(cvode_mem);

    nerr
}

fn PrintIntro2() {
    print!("\n\n-------------------------------------------------------------");
    print!("\n-------------------------------------------------------------");
    print!("\n\nProblem 2: ydot = A * y, where A is a banded lower\n");
    print!("triangular matrix derived from 2-D advection PDE\n\n");
    println!(" neq = {}, ml = {}, mu = {}", P2_NEQ, P2_ML, P2_MU);
    print!(
        " itol = {}, reltol = {}, abstol = {}",
        "CV_SS",
        fmt_g(RTOL, 0, 2),
        fmt_g(ATOL, 0, 2)
    );
}

fn PrintHeader2() {
    println!("\n      t        max.err      qu     hu ");
}

fn PrintOutput2(t: f64, erm: f64, qu: i32, hu: f64) {
    println!(
        "{}  {}   {:2}   {}",
        fmt_f(t, 10, 3),
        fmt_e(erm, 12, 4),
        qu,
        fmt_e(hu, 12, 4)
    );
}

fn f2(_t: f64, y: &NVector, ydot: &mut NVector, _user_data: &mut UserData) -> i32 {
    /*
       Excluding boundaries,

       ydot    = f    = -2 y    + alpha1 * y      + alpha2 * y
           i,j    i,j       i,j             i-1,j             i,j-1
    */

    for j in 0..P2_MESHY {
        for i in 0..P2_MESHX {
            let k = i + j * P2_MESHX;
            let mut d = -TWO * y.data[k];
            if i != 0 {
                d += P2_ALPH1 * y.data[k - 1];
            }
            if j != 0 {
                d += P2_ALPH2 * y.data[k - P2_MESHX];
            }
            ydot.data[k] = d;
        }
    }

    0
}

fn Jac2(
    _tn: f64,
    _y: &NVector,
    _fy: &NVector,
    j: &mut SUNMatrix,
    _user_data: &mut UserData,
    _tmp1: &mut NVector,
    _tmp2: &mut NVector,
    _tmp3: &mut NVector,
) -> i32 {
    let jb = match j {
        SUNMatrix::Band(m) => m,
        _ => return -1,
    };
    let s_mu = jb.s_mu as usize;

    /*
       The components of f(t,y) which depend on y_{i,j} are
       f_{i,j}, f_{i+1,j}, and f_{i,j+1}.
    */

    for jj in 0..P2_MESHY {
        for i in 0..P2_MESHX {
            let k = i + jj * P2_MESHX;
            let kthCol = jb.col_mut(k as i64);
            /* SM_COLUMN_ELEMENT_B(kthCol, i, j) = kthCol[(i)-(j)+s_mu] */
            kthCol[s_mu] = -TWO;
            if i != P2_MESHX - 1 {
                kthCol[1 + s_mu] = P2_ALPH1;
            }
            if jj != P2_MESHY - 1 {
                kthCol[P2_MESHX + s_mu] = P2_ALPH2;
            }
        }
    }

    0
}

fn MaxError(y: &NVector, t: f64) -> f64 {
    let mut ex = ZERO;
    let mut maxError = ZERO;
    let mut jfact_inv = ONE;

    if t == ZERO {
        return ZERO;
    }

    if t <= THIRTY {
        ex = (-TWO * t).exp();
    }

    for j in 0..P2_MESHY {
        let mut ifact_inv = ONE;
        for i in 0..P2_MESHX {
            let k = i + j * P2_MESHX;
            let yt = t.powf((i + j) as f64) * ex * ifact_inv * jfact_inv;
            let er = (y.data[k] - yt).abs();
            if er > maxError {
                maxError = er;
            }
            ifact_inv /= (i + 1) as f64;
        }
        jfact_inv /= (j + 1) as f64;
    }
    maxError
}

fn PrepareNextRun(
    sunctx: &SUNContext,
    cvode_mem: &mut CVodeMem,
    lmm: i32,
    miter: i32,
    y: &NVector,
    mu: i64,
    ml: i64,
) -> i32 {
    let mut retval;

    print!("\n\n-------------------------------------------------------------");

    print!("\n\nLinear Multistep Method : ");
    if lmm == CV_ADAMS {
        println!("ADAMS");
    } else {
        println!("BDF");
    }

    print!("Iteration               : ");
    if miter == FUNC {
        println!("FIXEDPOINT");

        /* create fixed point nonlinear solver object and attach to CVode */
        let nls = SUNNonlinSol_FixedPoint(y, 0, sunctx);
        retval = CVodeSetNonlinearSolver(cvode_mem, nls);
        if check_retval(retval, "CVodeSetNonlinearSolver") {
            return 1;
        }
    } else {
        println!("NEWTON");

        /* create Newton nonlinear solver object and attach to CVode */
        let nls = SUNNonlinSol_Newton(y, sunctx);
        retval = CVodeSetNonlinearSolver(cvode_mem, nls);
        if check_retval(retval, "CVodeSetNonlinearSolver") {
            return 1;
        }

        print!("Linear Solver           : ");

        match miter {
            DENSE_USER => {
                println!("Dense, User-Supplied Jacobian");

                /* Create dense SUNMatrix for use in linear solves */
                let a = SUNDenseMatrix(P1_NEQ, P1_NEQ, sunctx);

                /* Create dense SUNLinearSolver object for use by CVode */
                let ls = SUNLinSol_Dense(y, &a, sunctx);

                /* Attach the matrix and linear solver to CVode */
                retval = CVodeSetLinearSolver(cvode_mem, ls, Some(a));
                if check_retval(retval, "CVodeSetLinearSolver") {
                    return 1;
                }

                /* Set the user-supplied Jacobian routine Jac */
                retval = CVodeSetJacFn(cvode_mem, Some(Jac1));
                if check_retval(retval, "CVodeSetJacFn") {
                    return 1;
                }
            }

            DENSE_DQ => {
                println!("Dense, Difference Quotient Jacobian");

                /* Create dense SUNMatrix for use in linear solves */
                let a = SUNDenseMatrix(P1_NEQ, P1_NEQ, sunctx);

                /* Create dense SUNLinearSolver object for use by CVode */
                let ls = SUNLinSol_Dense(y, &a, sunctx);

                /* Attach the matrix and linear solver to CVode */
                retval = CVodeSetLinearSolver(cvode_mem, ls, Some(a));
                if check_retval(retval, "CVodeSetLinearSolver") {
                    return 1;
                }

                /* Use a difference quotient Jacobian */
                retval = CVodeSetJacFn(cvode_mem, None);
                if check_retval(retval, "CVodeSetJacFn") {
                    return 1;
                }
            }

            DIAG => {
                println!("Diagonal Jacobian");

                /* Call CVDiag to create/attach the CVODE-specific diagonal solver */
                retval = CVDiag(cvode_mem);
                if check_retval(retval, "CVDiag") {
                    return 1;
                }
            }

            BAND_USER => {
                println!("Band, User-Supplied Jacobian");

                /* Create band SUNMatrix for use in linear solves */
                let a = SUNBandMatrix(P2_NEQ, mu, ml, sunctx);

                /* Create banded SUNLinearSolver object for use by CVode */
                let ls = SUNLinSol_Band(y, &a, sunctx);

                /* Attach the matrix and linear solver to CVode */
                retval = CVodeSetLinearSolver(cvode_mem, ls, Some(a));
                if check_retval(retval, "CVodeSetLinearSolver") {
                    return 1;
                }

                /* Set the user-supplied Jacobian routine Jac */
                retval = CVodeSetJacFn(cvode_mem, Some(Jac2));
                if check_retval(retval, "CVodeSetJacFn") {
                    return 1;
                }
            }

            BAND_DQ => {
                println!("Band, Difference Quotient Jacobian");

                /* Create band SUNMatrix for use in linear solves */
                let a = SUNBandMatrix(P2_NEQ, mu, ml, sunctx);

                /* Create banded SUNLinearSolver object for use by CVode */
                let ls = SUNLinSol_Band(y, &a, sunctx);

                /* Attach the matrix and linear solver to CVode */
                retval = CVodeSetLinearSolver(cvode_mem, ls, Some(a));
                if check_retval(retval, "CVodeSetLinearSolver") {
                    return 1;
                }

                /* Use a difference quotient Jacobian */
                retval = CVodeSetJacFn(cvode_mem, None);
                if check_retval(retval, "CVodeSetJacFn") {
                    return 1;
                }
            }

            _ => {}
        }
    }

    retval
}

fn PrintErrOutput(tol_factor: f64) {
    print!(
        "\n\n Error exceeds {} * tolerance \n\n",
        fmt_g(tol_factor, 0, 6)
    );
}

fn PrintFinalStats(cvode_mem: &mut CVodeMem, miter: i32, ero: f64) {
    let (mut lenrw, mut leniw) = (0i64, 0i64);
    let (mut lenrwLS, mut leniwLS) = (0i64, 0i64);
    let (mut nst, mut nfe, mut nsetups, mut nni, mut ncfn, mut netf) =
        (0i64, 0i64, 0i64, 0i64, 0i64, 0i64);
    let mut nje = 0i64;
    let mut nfeLS = 0i64;

    let mut retval = CVodeGetWorkSpace(cvode_mem, &mut lenrw, &mut leniw);
    check_retval(retval, "CVodeGetWorkSpace");
    retval = CVodeGetNumSteps(cvode_mem, &mut nst);
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

    println!("\n Final statistics for this run:\n");
    println!(" CVode real workspace length              = {:4} ", lenrw);
    println!(" CVode integer workspace length           = {:4} ", leniw);
    println!(" Number of steps                          = {:4} ", nst);
    println!(" Number of f-s                            = {:4} ", nfe);
    println!(" Number of setups                         = {:4} ", nsetups);
    println!(" Number of nonlinear iterations           = {:4} ", nni);
    println!(" Number of nonlinear convergence failures = {:4} ", ncfn);
    println!(" Number of error test failures            = {:4} \n", netf);

    if miter != FUNC {
        if miter != DIAG {
            retval = CVodeGetNumJacEvals(cvode_mem, &mut nje);
            check_retval(retval, "CVodeGetNumJacEvals");
            retval = CVodeGetNumLinRhsEvals(cvode_mem, &mut nfeLS);
            check_retval(retval, "CVodeGetNumLinRhsEvals");
            retval = CVodeGetLinWorkSpace(cvode_mem, &mut lenrwLS, &mut leniwLS);
            check_retval(retval, "CVodeGetLinWorkSpace");
        } else {
            nje = nsetups;
            retval = CVDiagGetNumRhsEvals(cvode_mem, &mut nfeLS);
            check_retval(retval, "CVDiagGetNumRhsEvals");
            retval = CVDiagGetWorkSpace(cvode_mem, &mut lenrwLS, &mut leniwLS);
            check_retval(retval, "CVDiagGetWorkSpace");
        }
        println!(" Linear solver real workspace length      = {:4} ", lenrwLS);
        println!(" Linear solver integer workspace length   = {:4} ", leniwLS);
        println!(" Number of Jacobian evaluations           = {:4} ", nje);
        println!(" Number of f evals. in linear solver      = {:4} \n", nfeLS);
    }

    println!(" Error overrun = {} ", fmt_f(ero, 0, 3));
}

fn PrintErrInfo(nerr: i32) {
    print!("\n\n-------------------------------------------------------------");
    print!("\n-------------------------------------------------------------");
    println!("\n\n Number of errors encountered = {} ", nerr);
}

fn check_retval(retval: i32, funcname: &str) -> bool {
    if retval < 0 {
        eprintln!(
            "\nSUNDIALS_ERROR: {}() failed with retval = {}\n",
            funcname, retval
        );
        return true;
    }
    false
}
