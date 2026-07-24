/* -----------------------------------------------------------------
 * Translated from examples/idas/serial/idasRoberts_FSA_dns.c
 * (IDAS 7.7.0)
 * Programmer(s): Cosmin Petra and Radu Serban @ LLNL
 *
 * Example problem:
 *
 * This simple example problem for IDA, due to Robertson,
 * is from chemical kinetics, and consists of the following three
 * equations:
 *
 *      dy1/dt = -p1*y1 + p2*y2*y3
 *      dy2/dt = p1*y1 - p2*y2*y3 - p3*y2**2
 *         0   = y1 + y2 + y3 - 1
 *
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1, y2 = y3 = 0.The reaction rates are: p1=0.04,
 * p2=1e4, and p3=3e7
 *
 * Optionally, IDAS can compute sensitivities with respect to the
 * problem parameters p1, p2, and p3.
 * The sensitivity right hand side is given analytically through the
 * user routine fS (of type SensRhs1Fn).
 * Any of two sensitivity methods (SIMULTANEOUS and STAGGERED can be
 * used and sensitivities may be included in the error test or not
 *(error control set on SUNTRUE or SUNFALSE, respectively).
 *
 * Execution:
 *
 * If no sensitivities are desired:
 *    % idasRoberts_FSA_dns -nosensi
 * If sensitivities are to be computed:
 *    % idasRoberts_FSA_dns -sensi sensi_meth err_con
 * where sensi_meth is one of {sim, stg} and err_con is one of
 * {t, f}.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use idas_rs::sundials_utils::{fmt_e, fmt_g};
use idas_rs::*;

/* Problem Constants */

const NEQ: usize = 3; /* number of equations  */
const T0: f64 = 0.0; /* initial time */
const T1: f64 = 0.4; /* first output time */
const TMULT: f64 = 10.0; /* output time factor */
const NOUT: i32 = 12; /* number of output times */

const NS: usize = 3; /* number of sensitivities computed */

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/* Type : UserData */
struct RobData {
    p: [f64; 3], /* problem parameters */
    coef: f64,
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY IDAS
 *--------------------------------------------------------------------
 */

/*
 * Residual routine. Compute F(t,y,y',p).
 */
fn res(_t: f64, yy: &NVector, yp: &NVector, resval: &mut NVector, user_data: &mut UserData) -> i32 {
    let data = user_data.as_ref().unwrap().downcast_ref::<RobData>().unwrap();
    let p1 = data.p[0];
    let p2 = data.p[1];
    let p3 = data.p[2];

    let y1 = yy.data[0];
    let y2 = yy.data[1];
    let y3 = yy.data[2];

    let yp1 = yp.data[0];
    let yp2 = yp.data[1];

    resval.data[0] = yp1 + p1 * y1 - p2 * y2 * y3;
    resval.data[1] = yp2 - p1 * y1 + p2 * y2 * y3 + p3 * y2 * y2;
    resval.data[2] = y1 + y2 + y3 - ONE;

    0
}

/*
 * resS routine. Compute sensitivity r.h.s.
 */
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

    for is in 0..NS {
        let s1 = yyS[is].data[0];
        let s2 = yyS[is].data[1];
        let s3 = yyS[is].data[2];

        let sd1 = ypS[is].data[0];
        let sd2 = ypS[is].data[1];

        let mut rs1 = sd1 + p1 * s1 - p2 * y3 * s2 - p2 * y2 * s3;
        let mut rs2 = sd2 - p1 * s1 + p2 * y3 * s2 + p2 * y2 * s3 + 2.0 * p3 * y2 * s2;
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
            2 => {
                rs2 += y2 * y2;
            }
            _ => {}
        }

        resvalS[is].data[0] = rs1;
        resvalS[is].data[1] = rs2;
        resvalS[is].data[2] = rs3;
    }

    0
}

fn rhsQ(_t: f64, y: &NVector, _yp: &NVector, ypQ: &mut NVector, user_data: &mut UserData) -> i32 {
    let data = user_data.as_ref().unwrap().downcast_ref::<RobData>().unwrap();

    ypQ.data[0] = y.data[2];

    ypQ.data[1] =
        data.coef * (y.data[0] * y.data[0] + y.data[1] * y.data[1] + y.data[2] * y.data[2]);

    0
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * Process and verify arguments to idasfwddenx.
 */
fn ProcessArgs(args: &[String]) -> (bool, i32, bool) {
    /* C pre-initializes sensi to SUNFALSE; every path below assigns it or
    diverges (WrongArgs -> !), so Rust definite assignment replaces the
    initializer (and each path assigns exactly once, so no `mut`) */
    let sensi;
    let mut sensi_meth = -1;
    let mut err_con = false;

    if args.len() < 2 {
        WrongArgs(&args[0]);
    }

    if args[1] == "-nosensi" {
        sensi = false;
    } else if args[1] == "-sensi" {
        sensi = true;
    } else {
        WrongArgs(&args[0]);
    }

    if sensi {
        if args.len() != 4 {
            WrongArgs(&args[0]);
        }

        if args[2] == "sim" {
            sensi_meth = IDA_SIMULTANEOUS;
        } else if args[2] == "stg" {
            sensi_meth = IDA_STAGGERED;
        } else {
            WrongArgs(&args[0]);
        }

        if args[3] == "t" {
            err_con = true;
        } else if args[3] == "f" {
            err_con = false;
        } else {
            WrongArgs(&args[0]);
        }
    }

    (sensi, sensi_meth, err_con)
}

fn WrongArgs(name: &str) -> ! {
    println!("\nUsage: {} [-nosensi] [-sensi sensi_meth err_con]", name);
    println!("         sensi_meth = sim or stg");
    println!("         err_con    = t or f");

    std::process::exit(0);
}

fn PrintIC(y: &NVector, yp: &NVector) {
    let data = &y.data;
    print!("\n\nConsistent IC:\n");
    print!("\ty = ");
    print!("{} {} {} \n", fmt_e(data[0], 12, 4), fmt_e(data[1], 12, 4), fmt_e(data[2], 12, 4));

    let data = &yp.data;
    print!("\typ= ");
    print!("{} {} {} \n", fmt_e(data[0], 12, 4), fmt_e(data[1], 12, 4), fmt_e(data[2], 12, 4));
}

fn PrintSensIC(_y: &NVector, _yp: &NVector, yS: &[NVector], ypS: &[NVector]) {
    let sdata = &yS[0].data;
    print!("                  Sensitivity 1  ");

    print!("\n\ts1 = ");
    print!("{} {} {} \n", fmt_e(sdata[0], 12, 4), fmt_e(sdata[1], 12, 4), fmt_e(sdata[2], 12, 4));
    let sdata = &ypS[0].data;
    print!("\ts1'= ");
    print!("{} {} {} \n", fmt_e(sdata[0], 12, 4), fmt_e(sdata[1], 12, 4), fmt_e(sdata[2], 12, 4));

    print!("                  Sensitivity 2  ");
    let sdata = &yS[1].data;
    print!("\n\ts2 = ");
    print!("{} {} {} \n", fmt_e(sdata[0], 12, 4), fmt_e(sdata[1], 12, 4), fmt_e(sdata[2], 12, 4));
    let sdata = &ypS[1].data;
    print!("\ts2'= ");
    print!("{} {} {} \n", fmt_e(sdata[0], 12, 4), fmt_e(sdata[1], 12, 4), fmt_e(sdata[2], 12, 4));

    print!("                  Sensitivity 3  ");
    let sdata = &yS[2].data;
    print!("\n\ts3 = ");
    print!("{} {} {} \n", fmt_e(sdata[0], 12, 4), fmt_e(sdata[1], 12, 4), fmt_e(sdata[2], 12, 4));
    let sdata = &ypS[2].data;
    print!("\ts3'= ");
    print!("{} {} {} \n", fmt_e(sdata[0], 12, 4), fmt_e(sdata[1], 12, 4), fmt_e(sdata[2], 12, 4));
}

/*
 * Print current t, step count, order, stepsize, and solution.
 */
fn PrintOutput(ida_mem: &mut IDAMem, t: f64, u: &NVector) {
    let mut nst = 0i64;
    let mut qu = 0i32;
    let mut hu = 0.0f64;

    IDAGetNumSteps(ida_mem, &mut nst);
    IDAGetLastOrder(ida_mem, &mut qu);
    IDAGetLastStep(ida_mem, &mut hu);

    println!("{} {:2}  {} {:5}", fmt_e(t, 8, 3), qu, fmt_e(hu, 8, 3), nst);

    print!("                  Solution       ");

    print!(
        "{} {} {} \n",
        fmt_e(u.data[0], 12, 4),
        fmt_e(u.data[1], 12, 4),
        fmt_e(u.data[2], 12, 4)
    );
}

/*
 * Print sensitivities.
 */
fn PrintSensOutput(uS: &[NVector]) {
    let sdata = &uS[0].data;
    print!("                  Sensitivity 1  ");

    print!("{} {} {} \n", fmt_e(sdata[0], 12, 4), fmt_e(sdata[1], 12, 4), fmt_e(sdata[2], 12, 4));

    let sdata = &uS[1].data;
    print!("                  Sensitivity 2  ");

    print!("{} {} {} \n", fmt_e(sdata[0], 12, 4), fmt_e(sdata[1], 12, 4), fmt_e(sdata[2], 12, 4));

    let sdata = &uS[2].data;
    print!("                  Sensitivity 3  ");

    print!("{} {} {} \n", fmt_e(sdata[0], 12, 4), fmt_e(sdata[1], 12, 4), fmt_e(sdata[2], 12, 4));
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
    /* Process arguments */
    let args: Vec<String> = std::env::args().collect();
    let (sensi, sensi_meth, err_con) = ProcessArgs(&args);

    /* Create the SUNDIALS context object for this simulation */
    let sunctx = SUNContext_Create();

    /* User data structure */
    let data = RobData { p: [0.040, 1.0e4, 3.0e7], coef: 0.5 };
    let p = data.p;

    /* Initial conditions */
    let mut y = N_VNew_Serial(NEQ as i64, &sunctx);

    y.data[0] = ONE;
    y.data[1] = ZERO;
    y.data[2] = ZERO;

    let mut yp = N_VClone(&y);

    /* These initial conditions are NOT consistent. See IDACalcIC below. */
    yp.data[0] = 0.1;
    yp.data[1] = ZERO;
    yp.data[2] = ZERO;

    /* Create IDAS object */
    let mut ida_mem = IDACreate(&sunctx);

    /* Allocate space for IDAS */
    let retval = IDAInit(&mut ida_mem, res, T0, &y, &yp);
    if check_retval(retval, "IDAInit") {
        std::process::exit(1);
    }

    /* Specify scalar relative tol. and vector absolute tol. */
    let reltol = 1.0e-6;
    let mut abstol = N_VClone(&y);
    abstol.data[0] = 1.0e-8;
    abstol.data[1] = 1.0e-14;
    abstol.data[2] = 1.0e-6;
    let retval = IDASVtolerances(&mut ida_mem, reltol, &abstol);
    if check_retval(retval, "IDASVtolerances") {
        std::process::exit(1);
    }

    /* Set ID vector */
    let mut id = N_VClone(&y);
    id.data[0] = 1.0;
    id.data[1] = 1.0;
    id.data[2] = 0.0;
    let retval = IDASetId(&mut ida_mem, Some(&id));
    if check_retval(retval, "IDASetId") {
        std::process::exit(1);
    }

    /* Attach user data */
    let retval = IDASetUserData(&mut ida_mem, Some(Box::new(data)));
    if check_retval(retval, "IDASetUserData") {
        std::process::exit(1);
    }

    /* Create dense SUNMatrix for use in linear solves */
    let a = SUNDenseMatrix(NEQ as i64, NEQ as i64, &sunctx);

    /* Create dense SUNLinearSolver object */
    let ls = SUNLinSol_Dense(&y, &a, &sunctx);

    /* Attach the matrix and linear solver */
    let retval = IDASetLinearSolver(&mut ida_mem, ls, Some(a));
    if check_retval(retval, "IDASetLinearSolver") {
        std::process::exit(1);
    }

    print!("\n3-species chemical kinetics problem\n");

    /* Sensitivity-related settings */
    let mut yS: Vec<NVector> = Vec::new();
    let mut ypS: Vec<NVector> = Vec::new();
    if sensi {
        let pbar = [p[0], p[1], p[2]];

        yS = (0..NS).map(|_| N_VClone(&y)).collect();
        for s in yS.iter_mut() {
            N_VConst(ZERO, s);
        }

        ypS = (0..NS).map(|_| N_VClone(&y)).collect();
        for s in ypS.iter_mut() {
            N_VConst(ZERO, s);
        }

        /*
         * Only non-zero sensitivity I.C. are ypS[0]:
         * - Ith(ypS[0],1) = -ONE;
         * - Ith(ypS[0],2) =  ONE;
         *
         * They are not set. IDACalcIC also computes consistent IC for
         * sensitivities.
         */

        let retval = IDASensInit(&mut ida_mem, NS as i32, sensi_meth, Some(resS), &yS, &ypS);
        if check_retval(retval, "IDASensInit") {
            std::process::exit(1);
        }

        let retval = IDASensEEtolerances(&mut ida_mem);
        if check_retval(retval, "IDASensEEtolerances") {
            std::process::exit(1);
        }

        let retval = IDASetSensErrCon(&mut ida_mem, err_con);
        if check_retval(retval, "IDASetSensErrCon") {
            std::process::exit(1);
        }

        let retval = IDASetSensParams(&mut ida_mem, Some(&p), Some(&pbar), None);
        if check_retval(retval, "IDASetSensParams") {
            std::process::exit(1);
        }

        print!("Sensitivity: YES ");
        if sensi_meth == IDA_SIMULTANEOUS {
            print!("( SIMULTANEOUS +");
        } else {
            print!("( STAGGERED +");
        }
        if err_con {
            print!(" FULL ERROR CONTROL )");
        } else {
            print!(" PARTIAL ERROR CONTROL )");
        }
    } else {
        print!("Sensitivity: NO ");
    }

    /*----------------------------------------------------------
     *               Q U A D R A T U R E S
     * ---------------------------------------------------------*/
    let mut yQ = N_VNew_Serial(2, &sunctx);

    yQ.data[0] = 0.0;
    yQ.data[1] = 0.0;

    IDAQuadInit(&mut ida_mem, rhsQ, &yQ);

    let mut yQS: Vec<NVector> = Vec::new();
    if sensi {
        yQS = (0..NS).map(|_| N_VClone(&yQ)).collect();
        for s in yQS.iter_mut() {
            N_VConst(ZERO, s);
        }

        IDAQuadSensInit(&mut ida_mem, None, &yQS);
    }

    /* Call IDACalcIC to compute consistent initial conditions. If sensitivity
       is enabled, this function also try to find consistent IC for the
       sensitivities. */

    let retval = IDACalcIC(&mut ida_mem, IDA_YA_YDP_INIT, T1);
    if check_retval(retval, "IDACalcIC") {
        std::process::exit(1);
    }

    let retval = IDAGetConsistentIC(&mut ida_mem, Some(&mut y), Some(&mut yp));
    if check_retval(retval, "IDAGetConsistentIC") {
        std::process::exit(1);
    }

    PrintIC(&y, &yp);

    if sensi {
        IDAGetSensConsistentIC(&mut ida_mem, Some(&mut yS), Some(&mut ypS));
        PrintSensIC(&y, &yp, &yS, &ypS);
    }

    /* In loop over output points, call IDA, print results, test for error */

    print!("\n\n");
    print!("===========================================");
    print!("============================\n");
    print!("     T     Q       H      NST           y1");
    print!("           y2           y3    \n");
    print!("===========================================");
    print!("============================\n");

    let mut t = 0.0;
    let mut tout = T1;
    for _iout in 1..=NOUT {
        let retval = IDASolve(&mut ida_mem, tout, &mut t, &mut y, &mut yp, IDA_NORMAL);
        if check_retval(retval, "IDASolve") {
            break;
        }

        PrintOutput(&mut ida_mem, t, &y);

        if sensi {
            let retval = IDAGetSens(&ida_mem, &mut t, &mut yS);
            if check_retval(retval, "IDAGetSens") {
                break;
            }
            PrintSensOutput(&yS);
        }
        print!("-----------------------------------------");
        print!("------------------------------\n");

        tout *= TMULT;
    }

    print!("\nQuadrature:\n");
    IDAGetQuad(&ida_mem, &mut t, &mut yQ);
    println!("G:      {}", fmt_e(yQ.data[0], 10, 4));

    if sensi {
        IDAGetQuadSens(&ida_mem, &mut t, &mut yQS);
        println!("\nSensitivities at t={}:", fmt_g(t, 0, 6));
        println!("dG/dp1: {}", fmt_e(yQS[0].data[0], 11, 4));
        println!("dG/dp1: {}", fmt_e(yQS[1].data[0], 11, 4));
        println!("dG/dp1: {}", fmt_e(yQS[2].data[0], 11, 4));
    }

    /* Print final statistics to the screen */
    print!("\nFinal Statistics:\n");
    let mut stdout = std::io::stdout();
    IDAPrintAllStats(&mut ida_mem, &mut stdout, SUN_OUTPUTFORMAT_TABLE);

    /* Print final statistics to a file in CSV format */
    let mut fname = String::from("idasRoberts_FSA_dns_stats");
    if sensi {
        if sensi_meth == IDA_SIMULTANEOUS {
            fname.push_str("_-sensi_sim");
        } else {
            fname.push_str("_-sensi_stg");
        }
        if err_con {
            fname.push_str("_t");
        } else {
            fname.push_str("_f");
        }
    }
    fname.push_str(".csv");
    let mut fid = std::fs::File::create(&fname).expect("create csv");
    IDAPrintAllStats(&mut ida_mem, &mut fid, SUN_OUTPUTFORMAT_CSV);
    drop(fid);

    /* Free memory (RAII) */
    IDAFree(ida_mem);
}
