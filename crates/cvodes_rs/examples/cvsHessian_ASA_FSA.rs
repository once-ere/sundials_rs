/* -----------------------------------------------------------------
 * Translation of examples/cvodes/serial/cvsHessian_ASA_FSA.c
 * (SUNDIALS 7.7.0).
 *
 * Hessian through adjoint sensitivity example problem.
 *
 *        [ - p1 * y1^2 - y3 ]           [ 1 ]
 *   y' = [    - y2          ]    y(0) = [ 1 ]
 *        [ -p2^2 * y2 * y3  ]           [ 1 ]
 *
 *   p1 = 1.0
 *   p2 = 2.0
 *
 *           2
 *          /
 *   G(p) = |  0.5 * ( y1^2 + y2^2 + y3^2 ) dt
 *          /
 *          0
 *
 * Compute the gradient (ASA) and Hessian (FSA over ASA) of G(p).
 *
 * See D.B. Ozyurt and P.I. Barton, SISC 26(5) 1725-1743, 2005.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]

use cvodes_rs::cvodea::{
    CVodeAdjInit, CVodeB, CVodeCreateB, CVodeF, CVodeGetB, CVodeGetQuadB, CVodeInitBS,
    CVodeQuadInitBS, CVodeQuadSStolerancesB, CVodeSStolerancesB,
};
use cvodes_rs::cvodea_io::{CVodeGetAdjCVodeBmem, CVodeSetQuadErrConB, CVodeSetUserDataB};
use cvodes_rs::cvodes::{
    CVode, CVodeCreate, CVodeFree, CVodeGetQuad, CVodeGetQuadSens, CVodeGetSens, CVodeInit,
    CVodeQuadInit, CVodeQuadReInit, CVodeQuadSStolerances, CVodeQuadSensEEtolerances,
    CVodeQuadSensInit, CVodeReInit, CVodeSStolerances, CVodeSensEEtolerances, CVodeSensInit,
};
use cvodes_rs::cvodes_io::{
    CVodeGetIntegratorStats, CVodeGetNonlinSolvStats, CVodeGetQuadSensStats, CVodeGetQuadStats,
    CVodeGetSensStats, CVodeSetQuadErrCon, CVodeSetQuadSensErrCon, CVodeSetSensErrCon,
    CVodeSetUserData,
};
use cvodes_rs::cvodes_ls::{CVodeSetLinearSolver, CVodeSetLinearSolverB};
use cvodes_rs::sundials_utils::{fmt_e, fmt_g};
use cvodes_rs::*;

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

#[derive(Clone)]
struct HessData {
    p1: f64,
    p2: f64,
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY CVODES
 *--------------------------------------------------------------------
 */

fn f(_t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let data = user_data.as_ref().unwrap().downcast_ref::<HessData>().unwrap();
    let p1 = data.p1;
    let p2 = data.p2;

    let y1 = y.data[0];
    let y2 = y.data[1];
    let y3 = y.data[2];

    ydot.data[0] = -p1 * y1 * y1 - y3;
    ydot.data[1] = -y2;
    ydot.data[2] = -p2 * p2 * y2 * y3;

    0
}

fn fQ(_t: f64, y: &NVector, qdot: &mut NVector, _user_data: &mut UserData) -> i32 {
    let y1 = y.data[0];
    let y2 = y.data[1];
    let y3 = y.data[2];

    qdot.data[0] = 0.5 * (y1 * y1 + y2 * y2 + y3 * y3);

    0
}

#[allow(clippy::too_many_arguments)]
fn fS(
    _Ns: i32,
    _t: f64,
    y: &NVector,
    _ydot: &NVector,
    yS: &[NVector],
    ySdot: &mut [NVector],
    user_data: &mut UserData,
    _tmp1: &mut NVector,
    _tmp2: &mut NVector,
) -> i32 {
    let data = user_data.as_ref().unwrap().downcast_ref::<HessData>().unwrap();
    let p1 = data.p1;
    let p2 = data.p2;

    let y1 = y.data[0];
    let y2 = y.data[1];
    let y3 = y.data[2];

    /* 1st sensitivity RHS */

    let s1 = yS[0].data[0];
    let s2 = yS[0].data[1];
    let s3 = yS[0].data[2];

    let fys1 = -2.0 * p1 * y1 * s1 - s3;
    let fys2 = -s2;
    let fys3 = -p2 * p2 * y3 * s2 - p2 * p2 * y2 * s3;

    ySdot[0].data[0] = fys1 - y1 * y1;
    ySdot[0].data[1] = fys2;
    ySdot[0].data[2] = fys3;

    /* 2nd sensitivity RHS */

    let s1 = yS[1].data[0];
    let s2 = yS[1].data[1];
    let s3 = yS[1].data[2];

    let fys1 = -2.0 * p1 * y1 * s1 - s3;
    let fys2 = -s2;
    let fys3 = -p2 * p2 * y3 * s2 - p2 * p2 * y2 * s3;

    ySdot[1].data[0] = fys1;
    ySdot[1].data[1] = fys2;
    ySdot[1].data[2] = fys3 - 2.0 * p2 * y2 * y3;

    0
}

#[allow(clippy::too_many_arguments)]
fn fQS(
    _Ns: i32,
    _t: f64,
    y: &NVector,
    yS: &[NVector],
    _yQdot: &NVector,
    yQSdot: &mut [NVector],
    _user_data: &mut UserData,
    _tmp: &mut NVector,
    _tmpQ: &mut NVector,
) -> i32 {
    let y1 = y.data[0];
    let y2 = y.data[1];
    let y3 = y.data[2];

    /* 1st sensitivity RHS */

    let s1 = yS[0].data[0];
    let s2 = yS[0].data[1];
    let s3 = yS[0].data[2];

    yQSdot[0].data[0] = y1 * s1 + y2 * s2 + y3 * s3;

    /* 1st sensitivity RHS */

    let s1 = yS[1].data[0];
    let s2 = yS[1].data[1];
    let s3 = yS[1].data[2];

    yQSdot[1].data[0] = y1 * s1 + y2 * s2 + y3 * s3;

    0
}

fn fB1(
    _t: f64,
    y: &NVector,
    yS: &[NVector],
    yB: &NVector,
    yBdot: &mut NVector,
    user_dataB: &mut UserData,
) -> i32 {
    let data = user_dataB.as_ref().unwrap().downcast_ref::<HessData>().unwrap();
    let p1 = data.p1;
    let p2 = data.p2;

    let y1 = y.data[0]; /* solution */
    let y2 = y.data[1];
    let y3 = y.data[2];

    let s1 = yS[0].data[0]; /* sensitivity 1 */
    let s2 = yS[0].data[1];
    let s3 = yS[0].data[2];

    let l1 = yB.data[0]; /* lambda */
    let l2 = yB.data[1];
    let l3 = yB.data[2];

    let m1 = yB.data[3]; /* mu */
    let m2 = yB.data[4];
    let m3 = yB.data[5];

    yBdot.data[0] = 2.0 * p1 * y1 * l1 - y1;
    yBdot.data[1] = l2 + p2 * p2 * y3 * l3 - y2;
    yBdot.data[2] = l1 + p2 * p2 * y2 * l3 - y3;

    yBdot.data[3] = 2.0 * p1 * y1 * m1 + l1 * 2.0 * (y1 + p1 * s1) - s1;
    yBdot.data[4] = m2 + p2 * p2 * y3 * m3 + l3 * p2 * p2 * s3 - s2;
    yBdot.data[5] = m1 + p2 * p2 * y2 * m3 + l3 * p2 * p2 * s2 - s3;

    0
}

fn fQB1(
    _t: f64,
    y: &NVector,
    yS: &[NVector],
    yB: &NVector,
    qBdot: &mut NVector,
    user_dataB: &mut UserData,
) -> i32 {
    let data = user_dataB.as_ref().unwrap().downcast_ref::<HessData>().unwrap();

    let p2 = data.p2;

    let y1 = y.data[0]; /* solution */
    let y2 = y.data[1];
    let y3 = y.data[2];

    let s1 = yS[0].data[0]; /* sensitivity 1 */
    let s2 = yS[0].data[1];
    let s3 = yS[0].data[2];

    let l1 = yB.data[0]; /* lambda */
    let l3 = yB.data[2];

    let m1 = yB.data[3]; /* mu */
    let m3 = yB.data[5];

    qBdot.data[0] = -y1 * y1 * l1;
    qBdot.data[1] = -2.0 * p2 * y2 * y3 * l3;

    qBdot.data[2] = -y1 * y1 * m1 - l1 * 2.0 * y1 * s1;
    qBdot.data[3] = -2.0 * p2 * y2 * y3 * m3 - l3 * 2.0 * (p2 * y3 * s2 + p2 * y2 * s3);

    0
}

fn fB2(
    _t: f64,
    y: &NVector,
    yS: &[NVector],
    yB: &NVector,
    yBdot: &mut NVector,
    user_dataB: &mut UserData,
) -> i32 {
    let data = user_dataB.as_ref().unwrap().downcast_ref::<HessData>().unwrap();
    let p1 = data.p1;
    let p2 = data.p2;

    let y1 = y.data[0]; /* solution */
    let y2 = y.data[1];
    let y3 = y.data[2];

    let s1 = yS[1].data[0]; /* sensitivity 2 */
    let s2 = yS[1].data[1];
    let s3 = yS[1].data[2];

    let l1 = yB.data[0]; /* lambda */
    let l2 = yB.data[1];
    let l3 = yB.data[2];

    let m1 = yB.data[3]; /* mu */
    let m2 = yB.data[4];
    let m3 = yB.data[5];

    yBdot.data[0] = 2.0 * p1 * y1 * l1 - y1;
    yBdot.data[1] = l2 + p2 * p2 * y3 * l3 - y2;
    yBdot.data[2] = l1 + p2 * p2 * y2 * l3 - y3;

    yBdot.data[3] = 2.0 * p1 * y1 * m1 + l1 * 2.0 * p1 * s1 - s1;
    yBdot.data[4] = m2 + p2 * p2 * y3 * m3 + l3 * (2.0 * p2 * y3 + p2 * p2 * s3) - s2;
    yBdot.data[5] = m1 + p2 * p2 * y2 * m3 + l3 * (2.0 * p2 * y2 + p2 * p2 * s2) - s3;

    0
}

fn fQB2(
    _t: f64,
    y: &NVector,
    yS: &[NVector],
    yB: &NVector,
    qBdot: &mut NVector,
    user_dataB: &mut UserData,
) -> i32 {
    let data = user_dataB.as_ref().unwrap().downcast_ref::<HessData>().unwrap();

    let p2 = data.p2;

    let y1 = y.data[0]; /* solution */
    let y2 = y.data[1];
    let y3 = y.data[2];

    let s1 = yS[1].data[0]; /* sensitivity 2 */
    let s2 = yS[1].data[1];
    let s3 = yS[1].data[2];

    let l1 = yB.data[0]; /* lambda */
    let l3 = yB.data[2];

    let m1 = yB.data[3]; /* mu */
    let m3 = yB.data[5];

    qBdot.data[0] = -y1 * y1 * l1;
    qBdot.data[1] = -2.0 * p2 * y2 * y3 * l3;

    qBdot.data[2] = -y1 * y1 * m1 - l1 * 2.0 * y1 * s1;
    qBdot.data[3] =
        -2.0 * p2 * y2 * y3 * m3 - l3 * 2.0 * (p2 * y3 * s2 + p2 * y2 * s3 + y2 * y3);

    0
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

fn PrintFwdStats(cvode_mem: &mut CVodeMem) -> i32 {
    let (mut nst, mut nfe, mut nsetups, mut nni, mut ncfn, mut netf) = (0i64, 0, 0, 0, 0, 0);
    let (mut nfQe, mut netfQ) = (0i64, 0);
    let (mut nfSe, mut nfeS, mut nsetupsS, mut netfS) = (0i64, 0, 0, 0);
    let (mut nfQSe, mut netfQS) = (0i64, 0);

    let (mut qlast, mut qcur) = (0i32, 0);
    let (mut h0u, mut hlast, mut hcur, mut tcur) = (0.0, 0.0, 0.0, 0.0);

    let _ = CVodeGetIntegratorStats(
        cvode_mem, &mut nst, &mut nfe, &mut nsetups, &mut netf, &mut qlast, &mut qcur,
        &mut h0u, &mut hlast, &mut hcur, &mut tcur,
    );

    let _ = CVodeGetNonlinSolvStats(cvode_mem, &mut nni, &mut ncfn);

    let _ = CVodeGetQuadStats(cvode_mem, &mut nfQe, &mut netfQ);

    let _ = CVodeGetSensStats(cvode_mem, &mut nfSe, &mut nfeS, &mut netfS, &mut nsetupsS);

    let retval = CVodeGetQuadSensStats(cvode_mem, &mut nfQSe, &mut netfQS);

    println!(" Number steps: {:5}\n", nst);
    println!(" Function evaluations:");
    println!(
        "  f:        {:5}\n  fQ:       {:5}\n  fS:       {:5}\n  fQS:      {:5}",
        nfe, nfQe, nfSe, nfQSe
    );
    println!(" Error test failures:");
    println!(
        "  netf:     {:5}\n  netfQ:    {:5}\n  netfS:    {:5}\n  netfQS:   {:5}",
        netf, netfQ, netfS, netfQS
    );
    println!(" Linear solver setups:");
    println!("  nsetups:  {:5}\n  nsetupsS: {:5}", nsetups, nsetupsS);
    println!(" Nonlinear iterations:");
    println!("  nni:      {:5}", nni);
    println!(" Convergence failures:");
    println!("  ncfn:     {:5}", ncfn);

    println!();

    retval
}

fn PrintBckStats(cvode_mem: &mut CVodeMem, idx: i32) -> i32 {
    let (mut nst, mut nfe, mut nsetups, mut nni, mut ncfn, mut netf) = (0i64, 0, 0, 0, 0, 0);
    let (mut nfQe, mut netfQ) = (0i64, 0);

    let (mut qlast, mut qcur) = (0i32, 0);
    let (mut h0u, mut hlast, mut hcur, mut tcur) = (0.0, 0.0, 0.0, 0.0);

    let cvode_mem_bck = match CVodeGetAdjCVodeBmem(cvode_mem, idx) {
        Some(m) => m,
        None => return -1,
    };

    let _ = CVodeGetIntegratorStats(
        cvode_mem_bck, &mut nst, &mut nfe, &mut nsetups, &mut netf, &mut qlast, &mut qcur,
        &mut h0u, &mut hlast, &mut hcur, &mut tcur,
    );

    let _ = CVodeGetNonlinSolvStats(cvode_mem_bck, &mut nni, &mut ncfn);

    let retval = CVodeGetQuadStats(cvode_mem_bck, &mut nfQe, &mut netfQ);

    println!(" Number steps: {:5}\n", nst);
    println!(" Function evaluations:");
    println!("  f:        {:5}\n  fQ:       {:5}", nfe, nfQe);
    println!(" Error test failures:");
    println!("  netf:     {:5}\n  netfQ:    {:5}", netf, netfQ);
    println!(" Linear solver setups:");
    println!("  nsetups:  {:5}", nsetups);
    println!(" Nonlinear iterations:");
    println!("  nni:      {:5}", nni);
    println!(" Convergence failures:");
    println!("  ncfn:     {:5}", ncfn);

    println!();

    retval
}

/* Adjust the problem parameters held in the solver's user data
   (C mutates the shared UserData struct directly) */
fn adjust_params(cvode_mem: &mut CVodeMem, dp1: f64, dp2: f64) {
    if let Some(d) = cvode_mem.cv_user_data.as_mut() {
        if let Some(data) = d.downcast_mut::<HessData>() {
            data.p1 += dp1;
            data.p2 += dp2;
        }
    }
}

/* Check if a SUNDIALS function returned a negative value */
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
    /* User data structure */
    let data = HessData { p1: 1.0, p2: 2.0 };

    /* Problem size, integration interval, and tolerances */
    let neq: i64 = 3;
    let np: usize = 2;
    let np2: i64 = 2 * np as i64;

    let t0 = 0.0;
    let tf = 2.0;

    let reltol = 1.0e-8;

    let abstol = 1.0e-8;
    let abstolQ = 1.0e-8;

    let abstolB = 1.0e-8;
    let abstolQB = 1.0e-8;

    /* Create the SUNDIALS simulation context that all SUNDIALS objects require */
    let sunctx = SUNContext_Create();

    /* Initializations for forward problem */

    let mut y = N_VNew_Serial(neq, &sunctx);
    N_VConst(ONE, &mut y);

    let mut yQ = N_VNew_Serial(1, &sunctx);
    N_VConst(ZERO, &mut yQ);

    let mut yS: Vec<NVector> = (0..np).map(|_| N_VClone(&y)).collect();
    N_VConst(ZERO, &mut yS[0]);
    N_VConst(ZERO, &mut yS[1]);

    let mut yQS: Vec<NVector> = (0..np).map(|_| N_VClone(&yQ)).collect();
    N_VConst(ZERO, &mut yQS[0]);
    N_VConst(ZERO, &mut yQS[1]);

    /* Create and initialize forward problem */

    let mut cvode_mem = CVodeCreate(CV_BDF, &sunctx);

    let retval = CVodeInit(&mut cvode_mem, f, t0, &y);
    if check_retval(retval, "CVodeInit") {
        return;
    }

    let retval = CVodeSStolerances(&mut cvode_mem, reltol, abstol);
    if check_retval(retval, "CVodeSStolerances") {
        return;
    }

    let retval = CVodeSetUserData(&mut cvode_mem, Some(Box::new(data.clone())));
    if check_retval(retval, "CVodeSetUserData") {
        return;
    }

    /* Create a dense SUNMatrix and dense SUNLinearSolver */
    let a_mat = SUNDenseMatrix(neq, neq, &sunctx);
    let ls = SUNLinSol_Dense(&y, &a_mat, &sunctx);

    /* Attach the matrix and linear solver */
    let retval = CVodeSetLinearSolver(&mut cvode_mem, ls, Some(a_mat));
    if check_retval(retval, "CVodeSetLinearSolver") {
        return;
    }

    let retval = CVodeQuadInit(&mut cvode_mem, fQ, &yQ);
    if check_retval(retval, "CVodeQuadInit") {
        return;
    }

    let retval = CVodeQuadSStolerances(&mut cvode_mem, reltol, abstolQ);
    if check_retval(retval, "CVodeQuadSStolerances") {
        return;
    }

    let retval = CVodeSetQuadErrCon(&mut cvode_mem, SUNTRUE);
    if check_retval(retval, "CVodeSetQuadErrCon") {
        return;
    }

    let retval = CVodeSensInit(&mut cvode_mem, np as i32, CV_SIMULTANEOUS, Some(fS), &yS);
    if check_retval(retval, "CVodeSensInit") {
        return;
    }

    let retval = CVodeSensEEtolerances(&mut cvode_mem);
    if check_retval(retval, "CVodeSensEEtolerances") {
        return;
    }

    let retval = CVodeSetSensErrCon(&mut cvode_mem, SUNTRUE);
    if check_retval(retval, "CVodeSetSensErrCon") {
        return;
    }

    let retval = CVodeQuadSensInit(&mut cvode_mem, Some(fQS), &yQS);
    if check_retval(retval, "CVodeQuadSensInit") {
        return;
    }

    let retval = CVodeQuadSensEEtolerances(&mut cvode_mem);
    if check_retval(retval, "CVodeQuadSensEEtolerances") {
        return;
    }

    let retval = CVodeSetQuadSensErrCon(&mut cvode_mem, SUNTRUE);
    if check_retval(retval, "CVodeSetQuadSensErrCon") {
        return;
    }

    /* Initialize ASA */

    let steps = 100;
    let retval = CVodeAdjInit(&mut cvode_mem, steps, CV_POLYNOMIAL);
    if check_retval(retval, "CVodeAdjInit") {
        return;
    }

    /* Forward integration */

    println!("-------------------");
    println!("Forward integration");
    println!("-------------------\n");

    let mut time = 0.0;
    let mut ncheck: i32 = 0;
    let retval = CVodeF(&mut cvode_mem, tf, &mut y, &mut time, CV_NORMAL, &mut ncheck);
    if check_retval(retval, "CVodeF") {
        return;
    }

    let retval = CVodeGetQuad(&cvode_mem, &mut time, &mut yQ);
    if check_retval(retval, "CVodeGetQuad") {
        return;
    }

    let g_val = yQ.data[0];

    let retval = CVodeGetSens(&cvode_mem, &mut time, &mut yS);
    if check_retval(retval, "CVodeGetSens") {
        return;
    }

    let retval = CVodeGetQuadSens(&cvode_mem, &mut time, &mut yQS);
    if check_retval(retval, "CVodeGetQuadSens") {
        return;
    }

    println!("ncheck = {}", ncheck);
    println!();
    print!(
        "     y:    {} {} {}",
        fmt_e(y.data[0], 12, 4),
        fmt_e(y.data[1], 12, 4),
        fmt_e(y.data[2], 12, 4)
    );
    println!("     G:    {}", fmt_e(yQ.data[0], 12, 4));
    println!();
    println!(
        "     yS1:  {} {} {}",
        fmt_e(yS[0].data[0], 12, 4),
        fmt_e(yS[0].data[1], 12, 4),
        fmt_e(yS[0].data[2], 12, 4)
    );
    println!(
        "     yS2:  {} {} {}",
        fmt_e(yS[1].data[0], 12, 4),
        fmt_e(yS[1].data[1], 12, 4),
        fmt_e(yS[1].data[2], 12, 4)
    );
    println!();
    println!(
        "   dG/dp:  {} {}",
        fmt_e(yQS[0].data[0], 12, 4),
        fmt_e(yQS[1].data[0], 12, 4)
    );
    println!();

    println!("Final Statistics for forward pb.");
    println!("--------------------------------");
    let retval = PrintFwdStats(&mut cvode_mem);
    if check_retval(retval, "PrintFwdStats") {
        return;
    }

    /* Initializations for backward problems */

    let mut yB1 = N_VNew_Serial(2 * neq, &sunctx);
    N_VConst(ZERO, &mut yB1);

    let mut yQB1 = N_VNew_Serial(np2, &sunctx);
    N_VConst(ZERO, &mut yQB1);

    let mut yB2 = N_VNew_Serial(2 * neq, &sunctx);
    N_VConst(ZERO, &mut yB2);

    let mut yQB2 = N_VNew_Serial(np2, &sunctx);
    N_VConst(ZERO, &mut yQB2);

    /* Create and initialize backward problems (one for each column of the
       Hessian) */

    /* -------------------------
       First backward problem
       -------------------------*/

    let mut indexB1: i32 = 0;
    let retval = CVodeCreateB(&mut cvode_mem, CV_BDF, &mut indexB1);
    if check_retval(retval, "CVodeCreateB") {
        return;
    }

    let retval = CVodeInitBS(&mut cvode_mem, indexB1, fB1, tf, &yB1);
    if check_retval(retval, "CVodeInitBS") {
        return;
    }

    let retval = CVodeSStolerancesB(&mut cvode_mem, indexB1, reltol, abstolB);
    if check_retval(retval, "CVodeSStolerancesB") {
        return;
    }

    let retval = CVodeSetUserDataB(&mut cvode_mem, indexB1, Some(Box::new(data.clone())));
    if check_retval(retval, "CVodeSetUserDataB") {
        return;
    }

    let retval = CVodeQuadInitBS(&mut cvode_mem, indexB1, fQB1, &yQB1);
    if check_retval(retval, "CVodeQuadInitBS") {
        return;
    }

    let retval = CVodeQuadSStolerancesB(&mut cvode_mem, indexB1, reltol, abstolQB);
    if check_retval(retval, "CVodeQuadSStolerancesB") {
        return;
    }

    let retval = CVodeSetQuadErrConB(&mut cvode_mem, indexB1, SUNTRUE);
    if check_retval(retval, "CVodeSetQuadErrConB") {
        return;
    }

    /* Create a dense SUNMatrix and dense SUNLinearSolver */
    let aB1_mat = SUNDenseMatrix(2 * neq, 2 * neq, &sunctx);
    let lsB1 = SUNLinSol_Dense(&yB1, &aB1_mat, &sunctx);

    /* Attach the matrix and linear solver */
    let retval = CVodeSetLinearSolverB(&mut cvode_mem, indexB1, lsB1, Some(aB1_mat));
    if check_retval(retval, "CVodeSetLinearSolverB") {
        return;
    }

    /* -------------------------
       Second backward problem
       -------------------------*/

    let mut indexB2: i32 = 0;
    let retval = CVodeCreateB(&mut cvode_mem, CV_BDF, &mut indexB2);
    if check_retval(retval, "CVodeCreateB") {
        return;
    }

    let retval = CVodeInitBS(&mut cvode_mem, indexB2, fB2, tf, &yB2);
    if check_retval(retval, "CVodeInitBS") {
        return;
    }

    let retval = CVodeSStolerancesB(&mut cvode_mem, indexB2, reltol, abstolB);
    if check_retval(retval, "CVodeSStolerancesB") {
        return;
    }

    let retval = CVodeSetUserDataB(&mut cvode_mem, indexB2, Some(Box::new(data.clone())));
    if check_retval(retval, "CVodeSetUserDataB") {
        return;
    }

    let retval = CVodeQuadInitBS(&mut cvode_mem, indexB2, fQB2, &yQB2);
    if check_retval(retval, "CVodeQuadInitBS") {
        return;
    }

    let retval = CVodeQuadSStolerancesB(&mut cvode_mem, indexB2, reltol, abstolQB);
    if check_retval(retval, "CVodeQuadSStolerancesB") {
        return;
    }

    let retval = CVodeSetQuadErrConB(&mut cvode_mem, indexB2, SUNTRUE);
    if check_retval(retval, "CVodeSetQuadErrConB") {
        return;
    }

    /* Create a dense SUNMatrix and dense SUNLinearSolver */
    let aB2_mat = SUNDenseMatrix(2 * neq, 2 * neq, &sunctx);
    let lsB2 = SUNLinSol_Dense(&yB2, &aB2_mat, &sunctx);

    /* Attach the matrix and linear solver */
    let retval = CVodeSetLinearSolverB(&mut cvode_mem, indexB2, lsB2, Some(aB2_mat));
    if check_retval(retval, "CVodeSetLinearSolverB") {
        return;
    }

    /* Backward integration */

    println!("---------------------------------------------");
    println!("Backward integration ... (2 adjoint problems)");
    println!("---------------------------------------------\n");

    let retval = CVodeB(&mut cvode_mem, t0, CV_NORMAL);
    if check_retval(retval, "CVodeB") {
        return;
    }

    let retval = CVodeGetB(&mut cvode_mem, indexB1, &mut time, &mut yB1);
    if check_retval(retval, "CVodeGetB") {
        return;
    }

    let retval = CVodeGetQuadB(&mut cvode_mem, indexB1, &mut time, &mut yQB1);
    if check_retval(retval, "CVodeGetQuadB") {
        return;
    }

    let retval = CVodeGetB(&mut cvode_mem, indexB2, &mut time, &mut yB2);
    if check_retval(retval, "CVodeGetB") {
        return;
    }

    let retval = CVodeGetQuadB(&mut cvode_mem, indexB2, &mut time, &mut yQB2);
    if check_retval(retval, "CVodeGetQuadB") {
        return;
    }

    println!(
        "   dG/dp:  {} {}   (from backward pb. 1)",
        fmt_e(-yQB1.data[0], 12, 4),
        fmt_e(-yQB1.data[1], 12, 4)
    );
    println!(
        "           {} {}   (from backward pb. 2)",
        fmt_e(-yQB2.data[0], 12, 4),
        fmt_e(-yQB2.data[1], 12, 4)
    );
    println!();
    println!("   H = d2G/dp2:");
    println!("        (1)            (2)");
    println!(
        "  {}   {}",
        fmt_e(-yQB1.data[2], 12, 4),
        fmt_e(-yQB2.data[2], 12, 4)
    );
    println!(
        "  {}   {}",
        fmt_e(-yQB1.data[3], 12, 4),
        fmt_e(-yQB2.data[3], 12, 4)
    );
    println!();

    println!("Final Statistics for backward pb. 1");
    println!("-----------------------------------");
    let retval = PrintBckStats(&mut cvode_mem, indexB1);
    if check_retval(retval, "PrintBckStats") {
        return;
    }

    println!("Final Statistics for backward pb. 2");
    println!("-----------------------------------");
    let retval = PrintBckStats(&mut cvode_mem, indexB2);
    if check_retval(retval, "PrintBckStats") {
        return;
    }

    /* Free memory */

    CVodeFree(cvode_mem);

    /* Finite difference tests */

    let dp = 1.0e-2;

    println!("-----------------------");
    println!("Finite Difference tests");
    println!("-----------------------\n");

    println!("del_p = {}\n", fmt_g(dp, 0, 6));

    let mut cvode_mem = CVodeCreate(CV_BDF, &sunctx);

    N_VConst(ONE, &mut y);
    N_VConst(ZERO, &mut yQ);

    let retval = CVodeInit(&mut cvode_mem, f, t0, &y);
    if check_retval(retval, "CVodeInit") {
        return;
    }

    let retval = CVodeSStolerances(&mut cvode_mem, reltol, abstol);
    if check_retval(retval, "CVodeSStolerances") {
        return;
    }

    let retval = CVodeSetUserData(&mut cvode_mem, Some(Box::new(data.clone())));
    if check_retval(retval, "CVodeSetUserData") {
        return;
    }

    /* Create a dense SUNMatrix and dense SUNLinearSolver */
    let a_mat = SUNDenseMatrix(neq, neq, &sunctx);
    let ls = SUNLinSol_Dense(&y, &a_mat, &sunctx);

    /* Attach the matrix and linear solver */
    let retval = CVodeSetLinearSolver(&mut cvode_mem, ls, Some(a_mat));
    if check_retval(retval, "CVodeSetLinearSolver") {
        return;
    }

    let retval = CVodeQuadInit(&mut cvode_mem, fQ, &yQ);
    if check_retval(retval, "CVodeQuadInit") {
        return;
    }

    let retval = CVodeQuadSStolerances(&mut cvode_mem, reltol, abstolQ);
    if check_retval(retval, "CVodeQuadSStolerances") {
        return;
    }

    let retval = CVodeSetQuadErrCon(&mut cvode_mem, SUNTRUE);
    if check_retval(retval, "CVodeSetQuadErrCon") {
        return;
    }

    adjust_params(&mut cvode_mem, dp, 0.0); /* data->p1 += dp */

    let retval = CVode(&mut cvode_mem, tf, &mut y, &mut time, CV_NORMAL);
    if check_retval(retval, "CVode") {
        return;
    }

    let retval = CVodeGetQuad(&cvode_mem, &mut time, &mut yQ);
    if check_retval(retval, "CVodeGetQuad") {
        return;
    }

    let mut gp = yQ.data[0];

    print!(
        "p1+  y:   {} {} {}",
        fmt_e(y.data[0], 12, 4),
        fmt_e(y.data[1], 12, 4),
        fmt_e(y.data[2], 12, 4)
    );
    println!("     G:   {}", fmt_e(yQ.data[0], 12, 4));

    adjust_params(&mut cvode_mem, -2.0 * dp, 0.0); /* data->p1 -= 2.0*dp */

    N_VConst(ONE, &mut y);
    N_VConst(ZERO, &mut yQ);

    CVodeReInit(&mut cvode_mem, t0, &y);
    CVodeQuadReInit(&mut cvode_mem, &yQ);

    let retval = CVode(&mut cvode_mem, tf, &mut y, &mut time, CV_NORMAL);
    if check_retval(retval, "CVode") {
        return;
    }

    let retval = CVodeGetQuad(&cvode_mem, &mut time, &mut yQ);
    if check_retval(retval, "CVodeGetQuad") {
        return;
    }

    let mut gm = yQ.data[0];
    print!(
        "p1-  y:   {} {} {}",
        fmt_e(y.data[0], 12, 4),
        fmt_e(y.data[1], 12, 4),
        fmt_e(y.data[2], 12, 4)
    );
    println!("     G:   {}", fmt_e(yQ.data[0], 12, 4));

    adjust_params(&mut cvode_mem, dp, 0.0); /* data->p1 += dp */

    let mut grdG_fwd = [0.0; 2];
    let mut grdG_bck = [0.0; 2];
    let mut grdG_cntr = [0.0; 2];
    grdG_fwd[0] = (gp - g_val) / dp;
    grdG_bck[0] = (g_val - gm) / dp;
    grdG_cntr[0] = (gp - gm) / (2.0 * dp);
    let h11 = (gp - 2.0 * g_val + gm) / (dp * dp);

    adjust_params(&mut cvode_mem, 0.0, dp); /* data->p2 += dp */

    N_VConst(ONE, &mut y);
    N_VConst(ZERO, &mut yQ);

    CVodeReInit(&mut cvode_mem, t0, &y);
    CVodeQuadReInit(&mut cvode_mem, &yQ);

    let retval = CVode(&mut cvode_mem, tf, &mut y, &mut time, CV_NORMAL);
    if check_retval(retval, "CVode") {
        return;
    }

    let retval = CVodeGetQuad(&cvode_mem, &mut time, &mut yQ);
    if check_retval(retval, "CVodeGetQuad") {
        return;
    }

    gp = yQ.data[0];
    print!(
        "p2+  y:   {} {} {}",
        fmt_e(y.data[0], 12, 4),
        fmt_e(y.data[1], 12, 4),
        fmt_e(y.data[2], 12, 4)
    );
    println!("     G:   {}", fmt_e(yQ.data[0], 12, 4));

    adjust_params(&mut cvode_mem, 0.0, -2.0 * dp); /* data->p2 -= 2.0*dp */

    N_VConst(ONE, &mut y);
    N_VConst(ZERO, &mut yQ);

    CVodeReInit(&mut cvode_mem, t0, &y);
    CVodeQuadReInit(&mut cvode_mem, &yQ);

    let retval = CVode(&mut cvode_mem, tf, &mut y, &mut time, CV_NORMAL);
    if check_retval(retval, "CVode") {
        return;
    }

    let retval = CVodeGetQuad(&cvode_mem, &mut time, &mut yQ);
    if check_retval(retval, "CVodeGetQuad") {
        return;
    }

    gm = yQ.data[0];
    print!(
        "p2-  y:   {} {} {}",
        fmt_e(y.data[0], 12, 4),
        fmt_e(y.data[1], 12, 4),
        fmt_e(y.data[2], 12, 4)
    );
    println!("     G:   {}", fmt_e(yQ.data[0], 12, 4));

    adjust_params(&mut cvode_mem, 0.0, dp); /* data->p2 += dp */

    grdG_fwd[1] = (gp - g_val) / dp;
    grdG_bck[1] = (g_val - gm) / dp;
    grdG_cntr[1] = (gp - gm) / (2.0 * dp);
    let h22 = (gp - 2.0 * g_val + gm) / (dp * dp);

    println!();

    println!(
        "   dG/dp:  {} {}   (fwd FD)",
        fmt_e(grdG_fwd[0], 12, 4),
        fmt_e(grdG_fwd[1], 12, 4)
    );
    println!(
        "           {} {}   (bck FD)",
        fmt_e(grdG_bck[0], 12, 4),
        fmt_e(grdG_bck[1], 12, 4)
    );
    println!(
        "           {} {}   (cntr FD)",
        fmt_e(grdG_cntr[0], 12, 4),
        fmt_e(grdG_cntr[1], 12, 4)
    );
    println!();
    println!("  H(1,1):  {}", fmt_e(h11, 12, 4));
    println!("  H(2,2):  {}", fmt_e(h22, 12, 4));

    /* Free memory */

    CVodeFree(cvode_mem);

    drop(y);
    drop(yQ);
    drop(yS);
    drop(yQS);
    drop(yB1);
    drop(yQB1);
    drop(yB2);
    drop(yQB2);
}
