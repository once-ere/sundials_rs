/* -----------------------------------------------------------------
 * Translation of examples/cvodes/serial/cvsFoodWeb_ASAp_kry.c
 * (SUNDIALS 7.7.0).
 *
 * This program solves a stiff ODE system that arises from a system
 * of partial differential equations. The PDE system is a food web
 * population model, with predator-prey interaction and diffusion on
 * the unit square in two dimensions.
 *
 * The ODE system is solved by CVODES using Newton iteration and
 * the SUNLinSol_SPGMR linear solver (scaled preconditioned GMRES).
 * The preconditioner matrix used is the product of two matrices:
 * (1) a matrix, only defined implicitly, based on a fixed number
 * of Gauss-Seidel iterations using the diffusion terms only, and
 * (2) a block-diagonal matrix based on the partial derivatives of
 * the interaction terms f only, using block-grouping.
 *
 * Additionally, CVODES integrates backwards in time the
 * the semi-discrete form of the adjoint PDE:
 *   d(lambda)/dt = - D^T ( lambda_xx + lambda_yy )
 *                  - F_c^T lambda
 * with homogeneous Neumann boundary conditions and final conditions
 *   lambda(x,y,t=t_final) = - g_c^T(t_final)
 * whose solution at t = 0 represents the sensitivity of
 *   int_x int _y g(t_final,c) dx dy dt
 * with respect to the initial conditions of the original problem.
 * Here g(t,c) = c(ISPEC).
 *
 * Translation notes (same pinned pattern as cvsKrylovDemo_prec):
 * - The C code stores the cvode_mem pointer in the user data so the
 *   Precond routines can fetch current error weights with
 *   CVodeGetErrWeights. Rust callbacks cannot reach the integrator
 *   memory, so user efuns reproduce the CVodeSStolerances weights
 *   bit-for-bit and snapshot them into wdata.rewt / wdata.rewtB.
 *   For the backward problem the efun is installed directly on the
 *   backward CVodeMem (CVodeGetAdjCVodeBmem + CVodeWFtolerances,
 *   the exact composition CVodeSStolerancesB performs in C); its
 *   user data during CVodeB is the forward memory, through which it
 *   reaches the backward problem's own user data. The C code fetches
 *   both problems' weights into ONE wdata->rewt scratch at use time;
 *   the snapshots keep them in separate fields (rewt/rewtB) so each
 *   Precond reads exactly what its C fetch would return.
 * - The C Precond routines perturb the solution array in place and
 *   restore it; here the difference quotients read a mutable local
 *   copy.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use cvodes_rs::cvodea::{
    CVodeAdjInit, CVodeB, CVodeCreateB, CVodeF, CVodeGetB, CVodeInitB,
};
use cvodes_rs::cvodea_io::{CVodeGetAdjCVodeBmem, CVodeSetUserDataB};
use cvodes_rs::cvodes::{CVodeCreate, CVodeFree, CVodeInit, CVodeWFtolerances};
use cvodes_rs::cvodes_io::{CVodeSetMaxNumSteps, CVodeSetUserData};
use cvodes_rs::cvodes_ls::{
    CVodeSetLinearSolver, CVodeSetLinearSolverB, CVodeSetPreconditioner, CVodeSetPreconditionerB,
};
use cvodes_rs::sundials_dense::{
    SUNDlsMat_denseAddIdentity, SUNDlsMat_denseGETRF, SUNDlsMat_denseGETRS,
};
use cvodes_rs::sundials_utils::{fmt_e, fmt_f};
use cvodes_rs::*;
use std::cell::RefCell;
use std::rc::Rc;

/* The C program shares ONE WebData between the forward and backward
   problems (both CVodeSetUserData calls receive the same pointer), and
   that sharing is load-bearing: checkpoint-replay forward f/Precond
   calls interleave with fB/PrecondB writes to the same fsave and P
   blocks during the backward run. Both Rust user datas therefore hold
   the same Rc<RefCell<WebData>>. */
type WebRef = Rc<RefCell<WebData>>;

/* Constants */

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;

/* Problem Specification Constants */

const AA: f64 = ONE; /* AA = a */
const EE: f64 = 1.0e4; /* EE = e */
const GG: f64 = 0.5e-6; /* GG = g */
const BB: f64 = ONE; /* BB = b */
const DPREY: f64 = ONE;
const DPRED: f64 = 0.5;
const ALPHA: f64 = ONE;
const NP: usize = 3;
const NS: usize = 2 * NP;

/* Method Constants */

const MX: usize = 20;
const MY: usize = 20;
const MXNS: usize = MX * NS;
const AX: f64 = ONE;
const AY: f64 = ONE;
const DX: f64 = AX / (MX - 1) as f64;
const DY: f64 = AY / (MY - 1) as f64;
const MP: usize = NS;
const MQ: usize = MX * MY;
const MXMP: usize = MX * MP;
const NGX: usize = 2;
const NGY: usize = 2;
const NGRP: usize = NGX * NGY;
const ITMAX: i32 = 5;

/* CVodeInit Constants */

const NEQ: usize = NS * MX * MY;
const T0: f64 = ZERO;
const RTOL: f64 = 1.0e-5;
const ATOL: f64 = 1.0e-5;

/* Output Constants */

const TOUT: f64 = 10.0;

/* Adjoint calculation constants */
/* g = int_x int_y c(ISPEC) dy dx at t = Tfinal */

const NSTEPS: i64 = 80; /* check points every NSTEPS steps */
const ISPEC: usize = 6; /* species # in objective */

/* Note: The value for species i at mesh point (j,k) is stored in */
/* component number (i-1) + j*NS + k*NS*MX of an N_Vector,        */
/* where 1 <= i <= NS, 0 <= j < MX, 0 <= k < MY.                  */

/* Structure for user data */

struct WebData {
    P: Vec<DenseMatrix>, /* [NGRP] dense mp x mp blocks */
    pivot: Vec<[sunindextype; NS]>,
    ns: usize,
    mxns: usize,
    mp: usize,
    #[allow(dead_code)]
    mq: usize,
    mx: usize,
    my: usize,
    ngrp: usize,
    ngx: usize,
    ngy: usize,
    mxmp: usize,
    #[allow(dead_code)]
    jgx: [usize; NGX + 1],
    #[allow(dead_code)]
    jgy: [usize; NGY + 1],
    jigx: [usize; MX],
    jigy: [usize; MY],
    jxr: [usize; NGX],
    jyr: [usize; NGY],
    acoef: [[f64; NS]; NS],
    bcoef: [f64; NS],
    diff: [f64; NS],
    cox: [f64; NS],
    coy: [f64; NS],
    dx: f64,
    dy: f64,
    srur: f64,
    fsave: Vec<f64>,  /* [NEQ] */
    fBsave: Vec<f64>, /* [NEQ] */
    rewt: NVector,  /* NEQ (forward weights snapshot)          */
    vtemp: NVector, /* NEQ (shared GS temp, as the C wdata's)  */
    rewtB: NVector, /* NEQ (backward weights snapshot; C keeps
                       one rewt scratch and re-fetches at use) */
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY CVODES
 *--------------------------------------------------------------------
 */

/*
 * This routine computes the right-hand side of the ODE system and
 * returns it in cdot. The interaction rates are computed by calls to
 * WebRates, and these are saved in fsave for use in preconditioning.
 */
fn f(t: f64, c: &NVector, cdot: &mut NVector, user_data: &mut UserData) -> i32 {
    let cell = user_data
        .as_ref()
        .unwrap()
        .downcast_ref::<WebRef>()
        .unwrap()
        .clone();
    let wdata = &mut *cell.borrow_mut();
    let cdata = &c.data;

    let mxns = wdata.mxns;
    let ns = wdata.ns;
    let dx = wdata.dx;
    let dy = wdata.dy;
    {
        let cdotdata = &mut cdot.data;
        let WebData {
            fsave,
            cox,
            coy,
            acoef,
            bcoef,
            ..
        } = wdata;

        for jy in 0..MY {
            let y = jy as f64 * dy;
            let iyoff = mxns * jy;
            let idyu: i64 = if jy == MY - 1 {
                -(mxns as i64)
            } else {
                mxns as i64
            };
            let idyl: i64 = if jy == 0 { -(mxns as i64) } else { mxns as i64 };
            for jx in 0..MX {
                let x = jx as f64 * dx;
                let ic = iyoff + ns * jx;
                /* Get interaction rates at one point (x,y). */
                WebRates(x, y, t, &cdata[ic..], &mut fsave[ic..], acoef, bcoef, ns);
                let idxu: i64 = if jx == MX - 1 { -(ns as i64) } else { ns as i64 };
                let idxl: i64 = if jx == 0 { -(ns as i64) } else { ns as i64 };
                for i in 1..=ns {
                    let ici = ic + i - 1;
                    /* Do differencing in y. */
                    let dcyli = cdata[ici] - cdata[(ici as i64 - idyl) as usize];
                    let dcyui = cdata[(ici as i64 + idyu) as usize] - cdata[ici];
                    /* Do differencing in x. */
                    let dcxli = cdata[ici] - cdata[(ici as i64 - idxl) as usize];
                    let dcxui = cdata[(ici as i64 + idxu) as usize] - cdata[ici];
                    /* Collect terms and load cdot elements. */
                    cdotdata[ici] = coy[i - 1] * (dcyui - dcyli)
                        + cox[i - 1] * (dcxui - dcxli)
                        + fsave[ici];
                }
            }
        }
    }

    0
}

/*
 * This routine generates the block-diagonal part of the Jacobian
 * corresponding to the interaction rates, multiplies by -gamma,
 * adds the identity matrix, and calls SUNDlsMat_denseGETRF to do
 * the LU decomposition of each diagonal block. One block per group
 * is computed; the Jacobian elements are generated by difference
 * quotients using calls to the routine fblock.
 */
fn Precond(
    t: f64,
    c: &NVector,
    fc: &NVector,
    _jok: bool,
    jcurPtr: &mut bool,
    gamma: f64,
    user_data: &mut UserData,
) -> i32 {
    let cell = user_data
        .as_ref()
        .unwrap()
        .downcast_ref::<WebRef>()
        .unwrap()
        .clone();
    let wdata = &mut *cell.borrow_mut();

    /* C: CVodeGetErrWeights(wdata->cvode_mem, rewt); wdata.rewt is the
       snapshot the efun keeps identical to the integrator's ewt */
    let uround = SUN_UNIT_ROUNDOFF;

    let mp = wdata.mp;
    let srur = wdata.srur;
    let ngrp = wdata.ngrp;
    let ngx = wdata.ngx;
    let ngy = wdata.ngy;
    let mxmp = wdata.mxmp;
    let WebData {
        P,
        pivot,
        jxr,
        jyr,
        fsave,
        rewt,
        acoef,
        bcoef,
        ..
    } = wdata;

    /* Make mp calls to fblock to approximate each diagonal block of
       Jacobian. fsave contains the base value of the rate vector and r0
       is a minimum increment factor for the difference quotient. (The C
       code perturbs cdata in place and restores it; the local copy here
       sees the identical perturbed values.) */

    let mut cdata: Vec<f64> = c.data.clone();
    let mut f1 = [ZERO; NS];

    let fac = N_VWrmsNorm(fc, rewt);
    let rewtdata = &rewt.data;
    let mut r0 = 1000.0 * SUNRabs(gamma) * uround * NEQ as f64 * fac;
    if r0 == ZERO {
        r0 = ONE;
    }

    for igy in 0..ngy {
        let jy = jyr[igy];
        let if00 = jy * mxmp;
        for igx in 0..ngx {
            let jx = jxr[igx];
            let if0 = if00 + jx * mp;
            let ig = igx + igy * ngx;
            /* Generate ig-th diagonal block */
            for j in 0..mp {
                /* Generate the jth column as a difference quotient */
                let jj = if0 + j;
                let save = cdata[jj];
                let r = SUNMAX(srur * SUNRabs(save), r0 / rewtdata[jj]);
                cdata[jj] += r;
                let fac = -gamma / r;
                fblock(t, &cdata, jx, jy, &mut f1, acoef, bcoef);
                for i in 0..mp {
                    /* C: P[ig][j][i] (column j, row i) */
                    P[ig].set(i as i64, j as i64, (f1[i] - fsave[if0 + i]) * fac);
                }
                cdata[jj] = save;
            }
        }
    }

    /* Add identity matrix and do LU decompositions on blocks. */

    for ig in 0..ngrp {
        SUNDlsMat_denseAddIdentity(&mut P[ig]);
        let ier = SUNDlsMat_denseGETRF(&mut P[ig], &mut pivot[ig]);
        if ier != 0 {
            return 1;
        }
    }

    *jcurPtr = SUNTRUE;
    0
}

/*
 * This routine applies two inverse preconditioner matrices to the
 * vector r, using the interaction-only block-diagonal Jacobian with
 * block-grouping, denoted Jr, and Gauss-Seidel applied to the
 * diffusion contribution to the Jacobian, denoted Jd.
 */
#[allow(clippy::too_many_arguments)]
fn PSolve(
    _tn: f64,
    _c: &NVector,
    _fc: &NVector,
    r: &NVector,
    z: &mut NVector,
    gamma: f64,
    _delta: f64,
    _lr: i32,
    user_data: &mut UserData,
) -> i32 {
    let cell = user_data
        .as_ref()
        .unwrap()
        .downcast_ref::<WebRef>()
        .unwrap()
        .clone();
    let wdata = &mut *cell.borrow_mut();

    N_VScale(ONE, r, z);

    /* call GSIter for Gauss-Seidel iterations (the temporary vector
       wdata->vtemp is detached so z, x and wdata can be borrowed together) */

    let mut x = std::mem::take(&mut wdata.vtemp);
    GSIter(gamma, z, &mut x, wdata);
    wdata.vtemp = x;

    /* Do backsolves for inverse of block-diagonal preconditioner factor */

    let mx = wdata.mx;
    let my = wdata.my;
    let ngx = wdata.ngx;
    let mp = wdata.mp;

    let mut iv = 0usize;
    for jy in 0..my {
        let igy = wdata.jigy[jy];
        for jx in 0..mx {
            let igx = wdata.jigx[jx];
            let ig = igx + igy * ngx;
            SUNDlsMat_denseGETRS(&wdata.P[ig], &wdata.pivot[ig], &mut z.data[iv..iv + mp]);
            iv += mp;
        }
    }

    0
}

/*
 * This routine computes the right-hand side of the adjoint ODE system
 * and returns it in cBdot. The interaction rates are computed by calls
 * to WebRates, and these are saved in fsave for use in preconditioning.
 * The adjoint interaction rates are computed by calls to WebRatesB.
 */
fn fB(t: f64, c: &NVector, cB: &NVector, cBdot: &mut NVector, user_data: &mut UserData) -> i32 {
    let cell = user_data
        .as_ref()
        .unwrap()
        .downcast_ref::<WebRef>()
        .unwrap()
        .clone();
    let wdata = &mut *cell.borrow_mut();
    let cdata = &c.data;
    let cBdata = &cB.data;
    let cBdotdata = &mut cBdot.data;

    let mxns = wdata.mxns;
    let ns = wdata.ns;
    let dx = wdata.dx;
    let dy = wdata.dy;
    let WebData {
        fsave,
        fBsave,
        cox,
        coy,
        acoef,
        bcoef,
        ..
    } = wdata;

    for jy in 0..MY {
        let y = jy as f64 * dy;
        let iyoff = mxns * jy;
        let idyu: i64 = if jy == MY - 1 {
            -(mxns as i64)
        } else {
            mxns as i64
        };
        let idyl: i64 = if jy == 0 { -(mxns as i64) } else { mxns as i64 };
        for jx in 0..MX {
            let x = jx as f64 * dx;
            let ic = iyoff + ns * jx;
            /* Get interaction rates at one point (x,y). */
            WebRatesB(
                x,
                y,
                t,
                &cdata[ic..],
                &cBdata[ic..],
                &mut fsave[ic..],
                &mut fBsave[ic..],
                acoef,
                bcoef,
                ns,
            );
            let idxu: i64 = if jx == MX - 1 { -(ns as i64) } else { ns as i64 };
            let idxl: i64 = if jx == 0 { -(ns as i64) } else { ns as i64 };
            for i in 1..=ns {
                let ici = ic + i - 1;
                /* Do differencing in y. */
                let dcyli = cBdata[ici] - cBdata[(ici as i64 - idyl) as usize];
                let dcyui = cBdata[(ici as i64 + idyu) as usize] - cBdata[ici];
                /* Do differencing in x. */
                let dcxli = cBdata[ici] - cBdata[(ici as i64 - idxl) as usize];
                let dcxui = cBdata[(ici as i64 + idxu) as usize] - cBdata[ici];
                /* Collect terms and load cdot elements. */
                cBdotdata[ici] = -coy[i - 1] * (dcyui - dcyli)
                    - cox[i - 1] * (dcxui - dcxli)
                    - fBsave[ici];
            }
        }
    }

    0
}

/*
 * Preconditioner setup function for the backward problem
 */
#[allow(clippy::too_many_arguments)]
fn PrecondB(
    t: f64,
    c: &NVector,
    _cB: &NVector,
    fcB: &NVector,
    _jok: bool,
    jcurPtr: &mut bool,
    gamma: f64,
    user_data: &mut UserData,
) -> i32 {
    let cell = user_data
        .as_ref()
        .unwrap()
        .downcast_ref::<WebRef>()
        .unwrap()
        .clone();
    let wdata = &mut *cell.borrow_mut();

    /* C: CVodeGetErrWeights(CVodeGetAdjCVodeBmem(...), rewt); wdata.rewtB
       is the snapshot the backward efun keeps identical to that vector */
    let uround = SUN_UNIT_ROUNDOFF;

    let mp = wdata.mp;
    let srur = wdata.srur;
    let ngrp = wdata.ngrp;
    let ngx = wdata.ngx;
    let ngy = wdata.ngy;
    let mxmp = wdata.mxmp;
    let WebData {
        P,
        pivot,
        jxr,
        jyr,
        fsave,
        rewtB,
        acoef,
        bcoef,
        ..
    } = wdata;

    /* Make mp calls to fblock to approximate each diagonal block of
       Jacobian. fsave contains the base value of the rate vector and r0
       is a minimum increment factor for the difference quotient. */

    let mut cdata: Vec<f64> = c.data.clone();
    let mut f1 = [ZERO; NS];

    let fac = N_VWrmsNorm(fcB, rewtB);
    let rewtdata = &rewtB.data;
    let mut r0 = 1000.0 * SUNRabs(gamma) * uround * NEQ as f64 * fac;
    if r0 == ZERO {
        r0 = ONE;
    }

    for igy in 0..ngy {
        let jy = jyr[igy];
        let if00 = jy * mxmp;
        for igx in 0..ngx {
            let jx = jxr[igx];
            let if0 = if00 + jx * mp;
            let ig = igx + igy * ngx;
            /* Generate ig-th diagonal block */
            for j in 0..mp {
                /* Generate the jth column as a difference quotient */
                let jj = if0 + j;
                let save = cdata[jj];
                let r = SUNMAX(srur * SUNRabs(save), r0 / rewtdata[jj]);
                cdata[jj] += r;
                let fac = gamma / r;
                fblock(t, &cdata, jx, jy, &mut f1, acoef, bcoef);
                for i in 0..mp {
                    /* C: P[ig][i][j] (column i, row j — the transpose) */
                    P[ig].set(j as i64, i as i64, (f1[i] - fsave[if0 + i]) * fac);
                }
                cdata[jj] = save;
            }
        }
    }

    /* Add identity matrix and do LU decompositions on blocks. */

    for ig in 0..ngrp {
        SUNDlsMat_denseAddIdentity(&mut P[ig]);
        let ier = SUNDlsMat_denseGETRF(&mut P[ig], &mut pivot[ig]);
        if ier != 0 {
            return 1;
        }
    }

    *jcurPtr = SUNTRUE;
    0
}

/*
 * Preconditioner solve function for the backward problem
 */
#[allow(clippy::too_many_arguments)]
fn PSolveB(
    _tn: f64,
    _c: &NVector,
    _cB: &NVector,
    _fcB: &NVector,
    r: &NVector,
    z: &mut NVector,
    gamma: f64,
    _delta: f64,
    _lr: i32,
    user_data: &mut UserData,
) -> i32 {
    let cell = user_data
        .as_ref()
        .unwrap()
        .downcast_ref::<WebRef>()
        .unwrap()
        .clone();
    let wdata = &mut *cell.borrow_mut();

    N_VScale(ONE, r, z);

    /* call GSIter for Gauss-Seidel iterations
       (same routine but with gamma=-gamma) */

    let mut x = std::mem::take(&mut wdata.vtemp);
    GSIter(-gamma, z, &mut x, wdata);
    wdata.vtemp = x;

    /* Do backsolves for inverse of block-diagonal preconditioner factor */

    let mx = wdata.mx;
    let my = wdata.my;
    let ngx = wdata.ngx;
    let mp = wdata.mp;

    let mut iv = 0usize;
    for jy in 0..my {
        let igy = wdata.jigy[jy];
        for jx in 0..mx {
            let igx = wdata.jigx[jx];
            let ig = igx + igy * ngx;
            SUNDlsMat_denseGETRS(&wdata.P[ig], &wdata.pivot[ig], &mut z.data[iv..iv + mp]);
            iv += mp;
        }
    }

    0
}

/* Forward efun: set weights identically to the internal cvEwtSetSS for
   CVodeSStolerances(RTOL, ATOL) and keep the snapshot the C Precond
   would obtain from CVodeGetErrWeights. */
fn ewt(y: &NVector, w: &mut NVector, user_data: &mut UserData) -> i32 {
    let cell = user_data
        .as_ref()
        .unwrap()
        .downcast_ref::<WebRef>()
        .unwrap()
        .clone();
    let wdata = &mut *cell.borrow_mut();
    for i in 0..y.data.len() {
        w.data[i] = ONE / (RTOL * SUNRabs(y.data[i]) + ATOL);
    }
    wdata.rewt.data.copy_from_slice(&w.data);
    0
}

/* Backward efun (installed on the backward CVodeMem): during CVodeB its
   user data is the forward memory, through which the backward problem's
   own WebData is reached to keep the rewtB snapshot. */
fn ewtB(yB: &NVector, w: &mut NVector, user_data: &mut UserData) -> i32 {
    for i in 0..yB.data.len() {
        w.data[i] = ONE / (RTOL * SUNRabs(yB.data[i]) + ATOL);
    }
    if let Some(fwd) = user_data.as_mut().and_then(|d| d.downcast_mut::<CVodeMem>()) {
        if let Some(ca_mem) = fwd.cv_adj_mem.as_mut() {
            if let Some(cell) = ca_mem.cvB_mem[0]
                .cv_user_data
                .as_ref()
                .and_then(|d| d.downcast_ref::<WebRef>())
            {
                cell.borrow_mut().rewtB.data.copy_from_slice(&w.data);
            }
        }
    }
    0
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * Allocate space for user data structure
 */
fn AllocUserData(sunctx: &SUNContext) -> WebData {
    let mut P = Vec::with_capacity(NGRP);
    let mut pivot = Vec::with_capacity(NGRP);
    for _i in 0..NGRP {
        P.push(DenseMatrix::new(NS as i64, NS as i64));
        pivot.push([0 as sunindextype; NS]);
    }
    WebData {
        P,
        pivot,
        ns: 0,
        mxns: 0,
        mp: 0,
        mq: 0,
        mx: 0,
        my: 0,
        ngrp: 0,
        ngx: 0,
        ngy: 0,
        mxmp: 0,
        jgx: [0; NGX + 1],
        jgy: [0; NGY + 1],
        jigx: [0; MX],
        jigy: [0; MY],
        jxr: [0; NGX],
        jyr: [0; NGY],
        acoef: [[ZERO; NS]; NS],
        bcoef: [ZERO; NS],
        diff: [ZERO; NS],
        cox: [ZERO; NS],
        coy: [ZERO; NS],
        dx: ZERO,
        dy: ZERO,
        srur: ZERO,
        fsave: vec![ZERO; NEQ],
        fBsave: vec![ZERO; NEQ],
        rewt: N_VNew_Serial(NEQ as i64, sunctx),
        vtemp: N_VNew_Serial(NEQ as i64, sunctx),
        rewtB: N_VNew_Serial(NEQ as i64, sunctx),
    }
}

/*
 * Initialize user data structure
 */
fn InitUserData(wdata: &mut WebData) {
    wdata.ns = NS;
    let ns = wdata.ns;

    for j in 0..NS {
        for i in 0..NS {
            wdata.acoef[i][j] = ZERO;
        }
    }
    for j in 0..NP {
        for i in 0..NP {
            wdata.acoef[NP + i][j] = EE;
            wdata.acoef[i][NP + j] = -GG;
        }
        wdata.acoef[j][j] = -AA;
        wdata.acoef[NP + j][NP + j] = -AA;
        wdata.bcoef[j] = BB;
        wdata.bcoef[NP + j] = -BB;
        wdata.diff[j] = DPREY;
        wdata.diff[NP + j] = DPRED;
    }

    /* Set remaining problem parameters */

    wdata.mxns = MXNS;
    wdata.dx = DX;
    let dx = wdata.dx;
    wdata.dy = DY;
    let dy = wdata.dy;
    for i in 0..ns {
        wdata.cox[i] = wdata.diff[i] / (dx * dx);
        wdata.coy[i] = wdata.diff[i] / (dy * dy);
    }

    /* Set remaining method parameters */

    wdata.mp = MP;
    wdata.mq = MQ;
    wdata.mx = MX;
    wdata.my = MY;
    wdata.srur = SUN_UNIT_ROUNDOFF.sqrt();
    wdata.mxmp = MXMP;
    wdata.ngrp = NGRP;
    wdata.ngx = NGX;
    wdata.ngy = NGY;
    SetGroups(MX, NGX, &mut wdata.jgx, &mut wdata.jigx, &mut wdata.jxr);
    SetGroups(MY, NGY, &mut wdata.jgy, &mut wdata.jigy, &mut wdata.jyr);
}

/*
 * This routine sets arrays jg, jig, and jr describing
 * a uniform partition of (0,1,2,...,m-1) into ng groups.
 */
fn SetGroups(m: usize, ng: usize, jg: &mut [usize], jig: &mut [usize], jr: &mut [usize]) {
    let mper = m / ng; /* does integer division */
    for ig in 0..ng {
        jg[ig] = ig * mper;
    }
    jg[ng] = m;

    let ngm1 = ng - 1;
    let len1 = ngm1 * mper;
    for j in 0..len1 {
        jig[j] = j / mper;
    }
    for j in len1..m {
        jig[j] = ngm1;
    }

    for ig in 0..ngm1 {
        jr[ig] = ((2 * ig + 1) * mper - 1) / 2;
    }
    jr[ngm1] = (ngm1 * mper + m - 1) / 2;
}

/*
 * This routine computes and loads the vector of initial values.
 */
fn CInit(c: &mut NVector, wdata: &WebData) {
    let cdata = &mut c.data;
    let ns = wdata.ns;
    let mxns = wdata.mxns;
    let dx = wdata.dx;
    let dy = wdata.dy;

    let x_factor = 4.0 / (AX * AX);
    let y_factor = 4.0 / (AY * AY);
    for jy in 0..MY {
        let y = jy as f64 * dy;
        let argy = (y_factor * y * (AY - y)) * (y_factor * y * (AY - y));
        let iyoff = mxns * jy;
        for jx in 0..MX {
            let x = jx as f64 * dx;
            let argx = (x_factor * x * (AX - x)) * (x_factor * x * (AX - x));
            let ioff = iyoff + ns * jx;
            for i in 1..=ns {
                let ici = ioff + i - 1;
                cdata[ici] = 10.0 + i as f64 * argx * argy;
            }
        }
    }
}

/*
 * This function computes and loads the final values for the adjoint
 * variables
 */
fn CbInit(c: &mut NVector, _is: usize, wdata: &WebData) {
    let cdata = &mut c.data;
    let ns = wdata.ns;
    let mxns = wdata.mxns;

    let mut gu = [ZERO; NS];
    gu[ISPEC - 1] = ONE;

    for jy in 0..MY {
        let iyoff = mxns * jy;
        for jx in 0..MX {
            let ioff = iyoff + ns * jx;
            for i in 1..=ns {
                let ici = ioff + i - 1;
                cdata[ici] = gu[i - 1];
            }
        }
    }
}

/*
 * This routine computes the interaction rates for the species
 * c_1, ... ,c_ns (stored in c[0],...,c[ns-1]), at one spatial point
 * and at time t.
 */
#[allow(clippy::too_many_arguments)]
fn WebRates(
    x: f64,
    y: f64,
    _t: f64,
    c: &[f64],
    rate: &mut [f64],
    acoef: &[[f64; NS]; NS],
    bcoef: &[f64; NS],
    ns: usize,
) {
    for i in 0..ns {
        rate[i] = ZERO;
    }

    for j in 0..ns {
        for i in 0..ns {
            rate[i] += c[j] * acoef[i][j];
        }
    }

    let fac = ONE + ALPHA * x * y;
    for i in 0..ns {
        rate[i] = c[i] * (bcoef[i] * fac + rate[i]);
    }
}

/*
 * This routine computes the interaction rates for the backward problem
 */
#[allow(clippy::too_many_arguments)]
fn WebRatesB(
    x: f64,
    y: f64,
    _t: f64,
    c: &[f64],
    cB: &[f64],
    rate: &mut [f64],
    rateB: &mut [f64],
    acoef: &[[f64; NS]; NS],
    bcoef: &[f64; NS],
    ns: usize,
) {
    let fac = ONE + ALPHA * x * y;

    for i in 0..ns {
        rate[i] = bcoef[i] * fac;
    }

    for j in 0..ns {
        for i in 0..ns {
            rate[i] += acoef[i][j] * c[j];
        }
    }

    for i in 0..ns {
        rateB[i] = cB[i] * rate[i];
        rate[i] = c[i] * rate[i];
    }

    for j in 0..ns {
        for i in 0..ns {
            rateB[i] += acoef[j][i] * c[j] * cB[j];
        }
    }
}

/*
 * This routine computes one block of the interaction terms of the
 * system, namely block (jx,jy), for use in preconditioning.
 * Here jx and jy count from 0.
 */
fn fblock(
    t: f64,
    cdata: &[f64],
    jx: usize,
    jy: usize,
    cdotdata: &mut [f64],
    acoef: &[[f64; NS]; NS],
    bcoef: &[f64; NS],
) {
    let iblok = jx + jy * MX;
    let y = jy as f64 * DY;
    let x = jx as f64 * DX;
    let ic = NS * iblok;
    WebRates(x, y, t, &cdata[ic..], cdotdata, acoef, bcoef, NS);
}

/*
 * This routine performs ITMAX=5 Gauss-Seidel iterations to compute an
 * approximation to (P-inverse)*z, where P = I - gamma*Jd, and
 * Jd represents the diffusion contributions to the Jacobian.
 */
fn GSIter(gamma: f64, z: &mut NVector, x: &mut NVector, wdata: &WebData) {
    let ns = wdata.ns;
    let mx = wdata.mx;
    let my = wdata.my;
    let mxns = wdata.mxns;
    let cox = &wdata.cox;
    let coy = &wdata.coy;

    let mut beta = [ZERO; NS];
    let mut beta2 = [ZERO; NS];
    let mut cof1 = [ZERO; NS];
    let mut gam = [ZERO; NS];
    let mut gam2 = [ZERO; NS];

    /* Write matrix as P = D - L - U.
       Load local arrays beta, beta2, gam, gam2, and cof1. */

    for i in 0..ns {
        let temp = ONE / (ONE + TWO * gamma * (cox[i] + coy[i]));
        beta[i] = gamma * cox[i] * temp;
        beta2[i] = TWO * beta[i];
        gam[i] = gamma * coy[i] * temp;
        gam2[i] = TWO * gam[i];
        cof1[i] = temp;
    }

    /* Begin iteration loop.
       Load vector x with (D-inverse)*z for first iteration. */

    {
        let xd = &mut x.data;
        let zd = &z.data;
        for jy in 0..my {
            let iyoff = mxns * jy;
            for jx in 0..mx {
                let ic = iyoff + ns * jx;
                v_prod(xd, ic, &cof1, zd, ic, ns); /* x[ic+i] = cof1[i]z[ic+i] */
            }
        }
    }
    N_VConst(ZERO, z);

    /* Looping point for iterations. */

    for iter in 1..=ITMAX {
        /* Calculate (D-inverse)*U*x if not the first iteration. */

        if iter > 1 {
            let xd = &mut x.data;
            for jy in 0..my {
                let iyoff = mxns * jy;
                for jx in 0..mx {
                    /* order of loops matters */
                    let ic = iyoff + ns * jx;
                    let x_loc = if jx == 0 {
                        0
                    } else if jx == mx - 1 {
                        2
                    } else {
                        1
                    };
                    let y_loc = if jy == 0 {
                        0
                    } else if jy == my - 1 {
                        2
                    } else {
                        1
                    };
                    match 3 * y_loc + x_loc {
                        0 => {
                            /* jx == 0, jy == 0 */
                            /* x[ic+i] = beta2[i]x[ic+ns+i] + gam2[i]x[ic+mxns+i] */
                            v_sum_prods(xd, ic, &beta2, ic + ns, &gam2, ic + mxns, ns);
                        }
                        1 => {
                            /* 1 <= jx <= mx-2, jy == 0 */
                            /* x[ic+i] = beta[i]x[ic+ns+i] + gam2[i]x[ic+mxns+i] */
                            v_sum_prods(xd, ic, &beta, ic + ns, &gam2, ic + mxns, ns);
                        }
                        2 => {
                            /* jx == mx-1, jy == 0 */
                            /* x[ic+i] = gam2[i]x[ic+mxns+i] */
                            v_prod_self(xd, ic, &gam2, ic + mxns, ns);
                        }
                        3 => {
                            /* jx == 0, 1 <= jy <= my-2 */
                            /* x[ic+i] = beta2[i]x[ic+ns+i] + gam[i]x[ic+mxns+i] */
                            v_sum_prods(xd, ic, &beta2, ic + ns, &gam, ic + mxns, ns);
                        }
                        4 => {
                            /* 1 <= jx <= mx-2, 1 <= jy <= my-2 */
                            /* x[ic+i] = beta[i]x[ic+ns+i] + gam[i]x[ic+mxns+i] */
                            v_sum_prods(xd, ic, &beta, ic + ns, &gam, ic + mxns, ns);
                        }
                        5 => {
                            /* jx == mx-1, 1 <= jy <= my-2 */
                            /* x[ic+i] = gam[i]x[ic+mxns+i] */
                            v_prod_self(xd, ic, &gam, ic + mxns, ns);
                        }
                        6 => {
                            /* jx == 0, jy == my-1 */
                            /* x[ic+i] = beta2[i]x[ic+ns+i] */
                            v_prod_self(xd, ic, &beta2, ic + ns, ns);
                        }
                        7 => {
                            /* 1 <= jx <= mx-2, jy == my-1 */
                            /* x[ic+i] = beta[i]x[ic+ns+i] */
                            v_prod_self(xd, ic, &beta, ic + ns, ns);
                        }
                        8 => {
                            /* jx == mx-1, jy == my-1 */
                            /* x[ic+i] = 0.0 */
                            v_zero(xd, ic, ns);
                        }
                        _ => {}
                    }
                }
            }
        } /* end if (iter > 1) */

        /* Overwrite x with [(I - (D-inverse)*L)-inverse]*x. */

        {
            let xd = &mut x.data;
            for jy in 0..my {
                let iyoff = mxns * jy;
                for jx in 0..mx {
                    /* order of loops matters */
                    let ic = iyoff + ns * jx;
                    let x_loc = if jx == 0 {
                        0
                    } else if jx == mx - 1 {
                        2
                    } else {
                        1
                    };
                    let y_loc = if jy == 0 {
                        0
                    } else if jy == my - 1 {
                        2
                    } else {
                        1
                    };
                    match 3 * y_loc + x_loc {
                        0 => {
                            /* jx == 0, jy == 0 */
                        }
                        1 => {
                            /* 1 <= jx <= mx-2, jy == 0 */
                            /* x[ic+i] += beta[i]x[ic-ns+i] */
                            v_inc_by_prod(xd, ic, &beta, ic - ns, ns);
                        }
                        2 => {
                            /* jx == mx-1, jy == 0 */
                            /* x[ic+i] += beta2[i]x[ic-ns+i] */
                            v_inc_by_prod(xd, ic, &beta2, ic - ns, ns);
                        }
                        3 => {
                            /* jx == 0, 1 <= jy <= my-2 */
                            /* x[ic+i] += gam[i]x[ic-mxns+i] */
                            v_inc_by_prod(xd, ic, &gam, ic - mxns, ns);
                        }
                        4 => {
                            /* 1 <= jx <= mx-2, 1 <= jy <= my-2 */
                            /* x[ic+i] += beta[i]x[ic-ns+i] + gam[i]x[ic-mxns+i] */
                            v_inc_by_prod(xd, ic, &beta, ic - ns, ns);
                            v_inc_by_prod(xd, ic, &gam, ic - mxns, ns);
                        }
                        5 => {
                            /* jx == mx-1, 1 <= jy <= my-2 */
                            /* x[ic+i] += beta2[i]x[ic-ns+i] + gam[i]x[ic-mxns+i] */
                            v_inc_by_prod(xd, ic, &beta2, ic - ns, ns);
                            v_inc_by_prod(xd, ic, &gam, ic - mxns, ns);
                        }
                        6 => {
                            /* jx == 0, jy == my-1 */
                            /* x[ic+i] += gam2[i]x[ic-mxns+i] */
                            v_inc_by_prod(xd, ic, &gam2, ic - mxns, ns);
                        }
                        7 => {
                            /* 1 <= jx <= mx-2, jy == my-1 */
                            /* x[ic+i] += beta[i]x[ic-ns+i] + gam2[i]x[ic-mxns+i] */
                            v_inc_by_prod(xd, ic, &beta, ic - ns, ns);
                            v_inc_by_prod(xd, ic, &gam2, ic - mxns, ns);
                        }
                        8 => {
                            /* jx == mx-1, jy == my-1 */
                            /* x[ic+i] += beta2[i]x[ic-ns+i] + gam2[i]x[ic-mxns+i] */
                            v_inc_by_prod(xd, ic, &beta2, ic - ns, ns);
                            v_inc_by_prod(xd, ic, &gam2, ic - mxns, ns);
                        }
                        _ => {}
                    }
                }
            }
        }

        /* Add increment x to z : z <- z+x */

        z.linear_sum_with(ONE, ONE, x); /* N_VLinearSum(ONE, z, ONE, x, z) aliases z */
    }
}

/* Small Vector Kernels. In C u, q, w are pointers into the same array,
   so the kernels here take one slice plus offsets. */

/* u[i] += v[i]*w[i] with u = xd+uo, w = xd+wo */
fn v_inc_by_prod(xd: &mut [f64], uo: usize, v: &[f64], wo: usize, n: usize) {
    for i in 0..n {
        xd[uo + i] += v[i] * xd[wo + i];
    }
}

/* u[i] = p[i]*q[i] + v[i]*w[i] with u = xd+uo, q = xd+qo, w = xd+wo */
fn v_sum_prods(xd: &mut [f64], uo: usize, p: &[f64], qo: usize, v: &[f64], wo: usize, n: usize) {
    for i in 0..n {
        xd[uo + i] = p[i] * xd[qo + i] + v[i] * xd[wo + i];
    }
}

/* u[i] = v[i]*w[i] with u = xd+uo, w = xd+wo (v_prod, aliased form) */
fn v_prod_self(xd: &mut [f64], uo: usize, v: &[f64], wo: usize, n: usize) {
    for i in 0..n {
        xd[uo + i] = v[i] * xd[wo + i];
    }
}

/* u[i] = v[i]*w[i] with u = xd+uo, w = wd+wo (v_prod, two-array form) */
fn v_prod(xd: &mut [f64], uo: usize, v: &[f64], wd: &[f64], wo: usize, n: usize) {
    for i in 0..n {
        xd[uo + i] = v[i] * wd[wo + i];
    }
}

/* u[i] = 0 with u = xd+uo */
fn v_zero(xd: &mut [f64], uo: usize, n: usize) {
    for i in 0..n {
        xd[uo + i] = ZERO;
    }
}

/*
 * Print maximum sensitivity of G for each species
 */
fn PrintOutput(cB: &NVector, ns: usize, mxns: usize, wdata: &WebData) {
    let mut x = ZERO;
    let mut y = ZERO;

    let cdata = &cB.data;

    for i in 1..=ns {
        let mut cmax = ZERO;

        for jy in (0..MY).rev() {
            for jx in 0..MX {
                let cij = cdata[(i - 1) + jx * ns + jy * mxns];
                if SUNRabs(cij) > cmax {
                    cmax = cij;
                    x = jx as f64 * wdata.dx;
                    y = jy as f64 * wdata.dy;
                }
            }
        }

        println!("\nMaximum sensitivity with respect to I.C. of species {}", i);
        println!("  mu max = {}", fmt_e(cmax, 0, 6));
        println!("at");
        println!("  x = {}\n  y = {}", fmt_e(x, 0, 6), fmt_e(y, 0, 6));
    }
}

/*
 * Compute double space integral
 */
fn doubleIntgr(cdata: &[f64], i: usize, wdata: &WebData) -> f64 {
    let ns = wdata.ns;
    let mx = wdata.mx;
    let my = wdata.my;
    let mxns = wdata.mxns;
    let dx = wdata.dx;
    let dy = wdata.dy;

    let mut jy = 0usize;
    let mut intgr_x = cdata[(i - 1) + jy * mxns];
    for jx in 1..mx - 1 {
        intgr_x += TWO * cdata[(i - 1) + jx * ns + jy * mxns];
    }
    intgr_x += cdata[(i - 1) + (mx - 1) * ns + jy * mxns];
    intgr_x *= 0.5 * dx;

    let mut intgr_xy = intgr_x;

    for jy in 1..my - 1 {
        intgr_x = cdata[(i - 1) + jy * mxns];
        for jx in 1..mx - 1 {
            intgr_x += TWO * cdata[(i - 1) + jx * ns + jy * mxns];
        }
        intgr_x += cdata[(i - 1) + (mx - 1) * ns + jy * mxns];
        intgr_x *= 0.5 * dx;

        intgr_xy += TWO * intgr_x;
    }

    jy = my - 1;
    intgr_x = cdata[(i - 1) + jy * mxns];
    for jx in 1..mx - 1 {
        intgr_x += TWO * cdata[(i - 1) + jx * ns + jy * mxns];
    }
    intgr_x += cdata[(i - 1) + (mx - 1) * ns + jy * mxns];
    intgr_x *= 0.5 * dx;

    intgr_xy += intgr_x;

    intgr_xy *= 0.5 * dy;

    intgr_xy
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
    let abstol = ATOL;
    let reltol = RTOL;
    let reltolB = RTOL;
    let abstolB = ATOL;
    let _ = (reltolB, abstolB); /* used via the backward efun (ewtB) */

    /* Create the SUNDIALS simulation context that all SUNDIALS objects require */
    let sunctx = SUNContext_Create();

    /* Allocate and initialize user data (shared between the forward and
       backward problems, exactly as the C program shares one pointer) */

    let mut wdata = AllocUserData(&sunctx);
    InitUserData(&mut wdata);
    let wdata: WebRef = Rc::new(RefCell::new(wdata));

    /* Set-up forward problem */

    /* Initializations */
    let mut c = N_VNew_Serial(NEQ as i64, &sunctx);
    CInit(&mut c, &wdata.borrow());

    /* Call CVodeCreate/CVodeInit for forward run */
    println!("\nCreate and allocate CVODES memory for forward run");
    let mut cvode_mem = CVodeCreate(CV_BDF, &sunctx);
    let retval = CVodeSetUserData(&mut cvode_mem, Some(Box::new(wdata.clone())));
    if check_retval(retval, "CVodeSetUserData") {
        return;
    }
    let retval = CVodeInit(&mut cvode_mem, f, T0, &c);
    if check_retval(retval, "CVodeInit") {
        return;
    }
    /* C: CVodeSStolerances(cvode_mem, reltol, abstol); the user efun
       computes the identical weights and snapshots them for Precond */
    let _ = (reltol, abstol);
    let retval = CVodeWFtolerances(&mut cvode_mem, ewt);
    if check_retval(retval, "CVodeSStolerances") {
        return;
    }

    /* Create SUNLinSol_SPGMR linear solver for forward run */
    let ls = SUNLinSol_SPGMR(&c, SUN_PREC_LEFT, 0, &sunctx);

    /* Attach the linear solver */
    let retval = CVodeSetLinearSolver(&mut cvode_mem, ls, None);
    if check_retval(retval, "CVodeSetLinearSolver") {
        return;
    }

    /* Set the preconditioner solve and setup functions */
    let retval = CVodeSetPreconditioner(&mut cvode_mem, Some(Precond), Some(PSolve));
    if check_retval(retval, "CVodeSetPreconditioner") {
        return;
    }

    /* Call CVodeSetMaxNumSteps to set the maximum number of steps the
     * solver will take in an attempt to reach the next output time
     * during forward integration. */
    let retval = CVodeSetMaxNumSteps(&mut cvode_mem, 2500);
    if check_retval(retval, "CVodeSetMaxNumSteps") {
        return;
    }

    /* Set-up adjoint calculations */

    println!("\nAllocate global memory");
    let retval = CVodeAdjInit(&mut cvode_mem, NSTEPS, CV_HERMITE);
    if check_retval(retval, "CVodeAdjInit") {
        return;
    }

    /* Perform forward run */

    println!("\nForward integration");
    let mut t = 0.0;
    let mut ncheck: i32 = 0;
    let retval = CVodeF(&mut cvode_mem, TOUT, &mut c, &mut t, CV_NORMAL, &mut ncheck);
    if check_retval(retval, "CVodeF") {
        return;
    }

    println!("\nncheck = {}", ncheck);

    println!(
        "\n   g = int_x int_y c{}(Tfinal,x,y) dx dy = {} \n",
        ISPEC,
        fmt_f(doubleIntgr(&c.data, ISPEC, &wdata.borrow()), 0, 6)
    );

    /* Set-up backward problem */

    /* Allocate cB */
    let mut cB = N_VNew_Serial(NEQ as i64, &sunctx);
    /* Initialize cB */
    CbInit(&mut cB, ISPEC, &wdata.borrow());

    /* Create and allocate CVODES memory for backward run */
    println!("\nCreate and allocate CVODES memory for backward run");
    let mut indexB: i32 = 0;
    let retval = CVodeCreateB(&mut cvode_mem, CV_BDF, &mut indexB);
    if check_retval(retval, "CVodeCreateB") {
        return;
    }
    let retval = CVodeSetUserDataB(&mut cvode_mem, indexB, Some(Box::new(wdata.clone())));
    if check_retval(retval, "CVodeSetUserDataB") {
        return;
    }
    let retval = CVodeInitB(&mut cvode_mem, indexB, fB, TOUT, &cB);
    if check_retval(retval, "CVodeInitB") {
        return;
    }
    /* C: CVodeSStolerancesB(cvode_mem, indexB, reltolB, abstolB); the
       backward efun computes the identical weights and snapshots them
       for PrecondB (CVodeSStolerancesB is exactly CVodeSStolerances on
       the backward memory in C) */
    let retval = match CVodeGetAdjCVodeBmem(&mut cvode_mem, indexB) {
        Some(bmem) => CVodeWFtolerances(bmem, ewtB),
        None => -1,
    };
    if check_retval(retval, "CVodeSStolerancesB") {
        return;
    }

    /* Create SUNLinSol_SPGMR linear solver for backward run */
    let lsB = SUNLinSol_SPGMR(&cB, SUN_PREC_LEFT, 0, &sunctx);

    /* Attach the linear solver */
    let retval = CVodeSetLinearSolverB(&mut cvode_mem, indexB, lsB, None);
    if check_retval(retval, "CVodeSetLinearSolverB") {
        return;
    }

    /* Set the preconditioner solve and setup functions */
    let retval = CVodeSetPreconditionerB(&mut cvode_mem, indexB, Some(PrecondB), Some(PSolveB));
    if check_retval(retval, "CVodeSetPreconditionerB") {
        return;
    }

    /* Perform backward integration */

    println!("\nBackward integration");
    let retval = CVodeB(&mut cvode_mem, T0, CV_NORMAL);
    if check_retval(retval, "CVodeB") {
        return;
    }

    let retval = CVodeGetB(&mut cvode_mem, indexB, &mut t, &mut cB);
    if check_retval(retval, "CVodeGetB") {
        return;
    }
    PrintOutput(&cB, NS, MXNS, &wdata.borrow());

    /* Free all memory */
    CVodeFree(cvode_mem);

    drop(c);
    drop(cB);
}
