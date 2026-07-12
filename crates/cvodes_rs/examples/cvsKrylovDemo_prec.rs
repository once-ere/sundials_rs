/* -----------------------------------------------------------------
 * Translated from examples/cvode/serial/cvsKrylovDemo_prec.c (CVODE 7.7.0)
 *
 * Demonstration program for CVODES - Krylov linear solver.
 * ODE system from ns-species interaction PDE in 2 dimensions.
 *
 * This program solves a stiff ODE system that arises from a system
 * of partial differential equations. The PDE system is a food web
 * population model, with predator-prey interaction and diffusion on
 * the unit square in two dimensions. The dependent variable vector is:
 *
 *        1   2        ns
 *  c = (c , c , ..., c  )
 *
 * and the PDEs are as follows:
 *
 *    i               i      i
 *  dc /dt  =  d(i)*(c    + c   )  +  f (x,y,c)  (i=1,...,ns)
 *                    xx     yy        i
 *
 * where
 *
 *                 i          ns         j
 *  f (x,y,c)  =  c *(b(i) + sum a(i,j)*c )
 *   i                       j=1
 *
 * The number of species is ns = 2*np, with the first np being prey
 * and the last np being predators. The coefficients a(i,j), b(i),
 * d(i) are:
 *
 *  a(i,i) = -a  (all i)
 *  a(i,j) = -g  (i <= np, j > np)
 *  a(i,j) =  e  (i > np, j <= np)
 *  b(i) =  b*(1 + alpha*x*y)  (i <= np)
 *  b(i) = -b*(1 + alpha*x*y)  (i > np)
 *  d(i) = Dprey  (i <= np)
 *  d(i) = Dpred  (i > np)
 *
 * The spatial domain is the unit square. The final time is 10. The
 * boundary conditions are: normal derivative = 0. A polynomial in x
 * and y is used to set the initial conditions.
 *
 * The PDEs are discretized by central differencing on an MX by MY
 * mesh.
 *
 * The resulting ODE system is stiff.
 *
 * The ODE system is solved using Newton iteration and the
 * SUNLinSol_SPGMR linear solver (scaled preconditioned GMRES).
 *
 * The preconditioner matrix used is the product of two matrices:
 * (1) A matrix, only defined implicitly, based on a fixed number of
 * Gauss-Seidel iterations using the diffusion terms only. (2) A
 * block-diagonal matrix based on the partial derivatives of the
 * interaction terms f only, using block-grouping (computing only a
 * subset of the ns by ns blocks).
 *
 * Four different runs are made for this problem. The product
 * preconditioner is applied on the left and on the right. In each
 * case, both the modified and classical Gram-Schmidt options are
 * tested. In the series of runs, CVodeInit, SUNLinSol_SPGMR, and
 * CVSetLinearSolver are called only for the first run, whereas
 * CVodeReInit, SUNLinSol_SPGMRSetPrecType, and
 * SUNLinSol_SPGMRSetGSType are called for each of the remaining
 * three runs.
 *
 * A problem description, performance statistics at selected output
 * times, and final statistics are written to standard output. On
 * the first run, solution values are also printed at output times.
 * Error and warning messages are written to standard error, but
 * there should be no such messages.
 *
 * Translation notes:
 * - The C demo stores the cvode_mem pointer inside the user data so
 *   that Precond can fetch the current error weights with
 *   CVodeGetErrWeights. Rust callbacks cannot reach the integrator
 *   memory, so an equivalent user efun (CVodeWFtolerances) reproduces
 *   the CVodeSStolerances weights bit-for-bit
 *   (ewt_i = 1/(reltol*|y_i| + abstol), the exact cvEwtSetSS
 *   arithmetic) and keeps a snapshot in wdata.rewt, which is exactly
 *   what CVodeGetErrWeights would return at Precond time.
 * - The C Precond perturbs the solution array in place and restores
 *   it; here the difference quotients read a mutable local copy.
 * - SUNLinSol_SPGMRSetPrecType/SetGSType operate on the solver owned
 *   by the integrator (in C the LS handle remains usable after
 *   CVodeSetLinearSolver).
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use cvodes_rs::sundials_dense::{
    SUNDlsMat_denseAddIdentity, SUNDlsMat_denseGETRF, SUNDlsMat_denseGETRS,
};
use cvodes_rs::sundials_utils::{fmt_e, fmt_f, fmt_g};
use cvodes_rs::*;

/* Constants */

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

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

const MX: usize = 6;
const MY: usize = 6;
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

/* Spgmr/CVLS Constants */

const MAXL: i32 = 0; /* => use default = MIN(NEQ, 5)            */
const DELTA: f64 = ZERO; /* => use default = 0.05                   */

/* Output Constants */

const T1: f64 = 1.0e-8;
const TOUT_MULT: f64 = 10.0;
const DTOUT: f64 = ONE;
const NOUT: i32 = 18;

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
    fsave: [f64; NEQ],
    tmp: NVector,
    rewt: NVector,
}

/* Implementation */

fn main() {
    /* Create the SUNDIALS context */
    let sunctx = SUNContext_Create();

    /* Initializations */
    let mut c = N_VNew_Serial(NEQ as i64, &sunctx);
    let mut wdata = AllocUserData(&sunctx);
    InitUserData(&mut wdata);
    let ns = wdata.ns;
    let mxns = wdata.mxns;
    let dx = wdata.dx;
    let dy = wdata.dy;

    /* Print problem description */
    PrintIntro();

    /* In C cvode_mem is created on the first pass through the loop; here
       ownership requires creating it up front. */
    let mut cvode_mem = CVodeCreate(CV_BDF, &sunctx);
    let mut wdata_opt = Some(wdata);

    /* Loop over jpre and gstype (four cases) */
    for jpre in SUN_PREC_LEFT..=SUN_PREC_RIGHT {
        for gstype in SUN_MODIFIED_GS..=SUN_CLASSICAL_GS {
            /* Initialize c and print heading */
            CInit(&mut c, ns, mxns, dx, dy);
            PrintHeader(jpre, gstype);

            /* Call CVodeInit or CVodeReInit, then SUNLinSol_SPGMR to set up problem */

            let firstrun = (jpre == SUN_PREC_LEFT) && (gstype == SUN_MODIFIED_GS);
            if firstrun {
                let mut retval =
                    CVodeSetUserData(&mut cvode_mem, Some(Box::new(wdata_opt.take().unwrap())));
                if check_retval(retval, "CVodeSetUserData") {
                    std::process::exit(1);
                }

                retval = CVodeInit(&mut cvode_mem, f, T0, &c);
                if check_retval(retval, "CVodeInit") {
                    std::process::exit(1);
                }

                /* C: CVodeSStolerances(cvode_mem, reltol, abstol); the user efun
                   computes the identical weights and snapshots them for Precond */
                retval = CVodeWFtolerances(&mut cvode_mem, ewt);
                if check_retval(retval, "CVodeSStolerances") {
                    std::process::exit(1);
                }

                let LS = SUNLinSol_SPGMR(&c, jpre, MAXL, &sunctx);

                retval = CVodeSetLinearSolver(&mut cvode_mem, LS, None);
                if check_retval(retval, "CVodeSetLinearSolver") {
                    std::process::exit(1);
                }

                retval = SUNLinSol_SPGMRSetGSType(&mut cvode_mem, gstype);
                if check_retval(retval, "SUNLinSol_SPGMRSetGSType") {
                    std::process::exit(1);
                }

                retval = CVodeSetEpsLin(&mut cvode_mem, DELTA);
                if check_retval(retval, "CVodeSetEpsLin") {
                    std::process::exit(1);
                }

                retval = CVodeSetPreconditioner(&mut cvode_mem, Some(Precond), Some(PSolve));
                if check_retval(retval, "CVodeSetPreconditioner") {
                    std::process::exit(1);
                }
            } else {
                let mut retval = CVodeReInit(&mut cvode_mem, T0, &c);
                if check_retval(retval, "CVodeReInit") {
                    std::process::exit(1);
                }

                retval = SUNLinSol_SPGMRSetPrecType(&mut cvode_mem, jpre);
                check_retval(retval, "SUNLinSol_SPGMRSetPrecType");
                retval = SUNLinSol_SPGMRSetGSType(&mut cvode_mem, gstype);
                if check_retval(retval, "SUNLinSol_SPGMRSetGSType") {
                    std::process::exit(1);
                }
            }

            /* Print initial values */
            if firstrun {
                PrintAllSpecies(&c, ns, mxns, T0);
            }

            /* Loop over output points, call CVode, print sample solution values. */
            let mut t = ZERO;
            let mut tout = T1;
            for iout in 1..=NOUT {
                let retval = CVode(&mut cvode_mem, tout, &mut c, &mut t, CV_NORMAL);
                PrintOutput(&mut cvode_mem, t);
                if firstrun && (iout % 3 == 0) {
                    PrintAllSpecies(&c, ns, mxns, t);
                }
                if check_retval(retval, "CVode") {
                    break;
                }
                if tout > 0.9 {
                    tout += DTOUT;
                } else {
                    tout *= TOUT_MULT;
                }
            }

            /* Print final statistics, and loop for next case */
            PrintFinalStats(&mut cvode_mem);
        }
    }

    /* Free all memory */
    CVodeFree(cvode_mem);
}

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
        fsave: [ZERO; NEQ],
        rewt: N_VNew_Serial(NEQ as i64, sunctx),
        tmp: N_VNew_Serial(NEQ as i64, sunctx),
    }
}

fn InitUserData(wdata: &mut WebData) {
    wdata.ns = NS;
    let ns = wdata.ns;

    for j in 0..NS {
        for i in 0..NS {
            wdata.acoef[i][j] = 0.;
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
 This routine sets arrays jg, jig, and jr describing
 a uniform partition of (0,1,2,...,m-1) into ng groups.
 The arrays set are:
   jg    = length ng+1 array of group boundaries.
           Group ig has indices j = jg[ig],...,jg[ig+1]-1.
   jig   = length m array of group indices vs node index.
           Node index j is in group jig[j].
   jr    = length ng array of indices representing the groups.
           The index for group ig is j = jr[ig].
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

/* This routine computes and loads the vector of initial values.
   (ns, mxns, dx, dy are the wdata fields the C version reads.) */
fn CInit(c: &mut NVector, ns: usize, mxns: usize, dx: f64, dy: f64) {
    let cdata = &mut c.data;

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

fn PrintIntro() {
    print!("\n\nDemonstration program for CVODES - SPGMR linear solver\n\n");
    print!("Food web problem with ns species, ns = {}\n", NS);
    print!("Predator-prey interaction and diffusion on a 2-D square\n\n");
    print!(
        "Matrix parameters: a = {}   e = {}   g = {}\n",
        fmt_g(AA, 0, 2),
        fmt_g(EE, 0, 2),
        fmt_g(GG, 0, 2)
    );
    print!("b parameter = {}\n", fmt_g(BB, 0, 2));
    print!(
        "Diffusion coefficients: Dprey = {}   Dpred = {}\n",
        fmt_g(DPREY, 0, 2),
        fmt_g(DPRED, 0, 2)
    );
    print!("Rate parameter alpha = {}\n\n", fmt_g(ALPHA, 0, 2));
    print!("Mesh dimensions (mx,my) are {}, {}.  ", MX, MY);
    print!("Total system size is neq = {} \n\n", NEQ);
    print!(
        "Tolerances: reltol = {}, abstol = {} \n\n",
        fmt_g(RTOL, 0, 2),
        fmt_g(ATOL, 0, 2)
    );
    print!("Preconditioning uses a product of:\n");
    print!("  (1) Gauss-Seidel iterations with ");
    print!("itmax = {} iterations, and\n", ITMAX);
    print!("  (2) interaction-only block-diagonal matrix ");
    print!("with block-grouping\n");
    print!("  Number of diagonal block groups = ngrp = {}", NGRP);
    print!("  (ngx by ngy, ngx = {}, ngy = {})\n", NGX, NGY);
    print!("\n\n--------------------------------------------------------------");
    print!("--------------\n");
}

fn PrintHeader(jpre: i32, gstype: i32) {
    if jpre == SUN_PREC_LEFT {
        print!(
            "\n\nPreconditioner type is           jpre = {}\n",
            "SUN_PREC_LEFT"
        );
    } else {
        print!(
            "\n\nPreconditioner type is           jpre = {}\n",
            "SUN_PREC_RIGHT"
        );
    }

    if gstype == SUN_MODIFIED_GS {
        print!(
            "\nGram-Schmidt method type is    gstype = {}\n\n\n",
            "SUN_MODIFIED_GS"
        );
    } else {
        print!(
            "\nGram-Schmidt method type is    gstype = {}\n\n\n",
            "SUN_CLASSICAL_GS"
        );
    }
}

fn PrintAllSpecies(c: &NVector, ns: usize, mxns: usize, t: f64) {
    let cdata = &c.data;
    print!("c values at t = {}:\n\n", fmt_g(t, 0, 6));
    for i in 1..=ns {
        print!("Species {}\n", i);
        for jy in (0..MY).rev() {
            for jx in 0..MX {
                /* C "%-10.6g": left-justified, width 10 */
                print!("{:<10}", fmt_g(cdata[(i - 1) + jx * ns + jy * mxns], 0, 6));
            }
            print!("\n");
        }
        print!("\n");
    }
}

fn PrintOutput(cvode_mem: &mut CVodeMem, t: f64) {
    let mut nst = 0i64;
    let mut nfe = 0i64;
    let mut nni = 0i64;
    let mut qu = 0i32;
    let mut hu = ZERO;

    let mut retval = CVodeGetNumSteps(cvode_mem, &mut nst);
    check_retval(retval, "CVodeGetNumSteps");
    retval = CVodeGetNumRhsEvals(cvode_mem, &mut nfe);
    check_retval(retval, "CVodeGetNumRhsEvals");
    retval = CVodeGetNumNonlinSolvIters(cvode_mem, &mut nni);
    check_retval(retval, "CVodeGetNumNonlinSolvIters");
    retval = CVodeGetLastOrder(cvode_mem, &mut qu);
    check_retval(retval, "CVodeGetLastOrder");
    retval = CVodeGetLastStep(cvode_mem, &mut hu);
    check_retval(retval, "CVodeGetLastStep");

    print!(
        "t = {}  nst = {}  nfe = {}  nni = {}",
        fmt_e(t, 10, 2),
        nst,
        nfe,
        nni
    );
    print!("  qu = {}  hu = {}\n\n", qu, fmt_e(hu, 11, 2));
}

fn PrintFinalStats(cvode_mem: &mut CVodeMem) {
    let (mut lenrw, mut leniw) = (0i64, 0i64);
    let (mut lenrwLS, mut leniwLS) = (0i64, 0i64);
    let (mut nst, mut nfe, mut nsetups, mut nni, mut ncfn, mut netf) =
        (0i64, 0i64, 0i64, 0i64, 0i64, 0i64);
    let (mut nli, mut npe, mut nps, mut ncfl, mut nfeLS) = (0i64, 0i64, 0i64, 0i64, 0i64);

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

    retval = CVodeGetLinWorkSpace(cvode_mem, &mut lenrwLS, &mut leniwLS);
    check_retval(retval, "CVodeGetLinWorkSpace");
    retval = CVodeGetNumLinIters(cvode_mem, &mut nli);
    check_retval(retval, "CVodeGetNumLinIters");
    retval = CVodeGetNumPrecEvals(cvode_mem, &mut npe);
    check_retval(retval, "CVodeGetNumPrecEvals");
    retval = CVodeGetNumPrecSolves(cvode_mem, &mut nps);
    check_retval(retval, "CVodeGetNumPrecSolves");
    retval = CVodeGetNumLinConvFails(cvode_mem, &mut ncfl);
    check_retval(retval, "CVodeGetNumLinConvFails");
    retval = CVodeGetNumLinRhsEvals(cvode_mem, &mut nfeLS);
    check_retval(retval, "CVodeGetNumLinRhsEvals");

    print!("\n\n Final statistics for this run:\n\n");
    print!(" CVode real workspace length           = {:4} \n", lenrw);
    print!(" CVode integer workspace length        = {:4} \n", leniw);
    print!(" CVLS real workspace length            = {:4} \n", lenrwLS);
    print!(" CVLS integer workspace length         = {:4} \n", leniwLS);
    print!(" Number of steps                       = {:4} \n", nst);
    print!(" Number of f-s                         = {:4} \n", nfe);
    print!(" Number of f-s (SPGMR)                 = {:4} \n", nfeLS);
    print!(" Number of f-s (TOTAL)                 = {:4} \n", nfe + nfeLS);
    print!(" Number of setups                      = {:4} \n", nsetups);
    print!(" Number of nonlinear iterations        = {:4} \n", nni);
    print!(" Number of linear iterations           = {:4} \n", nli);
    print!(" Number of preconditioner evaluations  = {:4} \n", npe);
    print!(" Number of preconditioner solves       = {:4} \n", nps);
    print!(" Number of error test failures         = {:4} \n", netf);
    print!(" Number of nonlinear conv. failures    = {:4} \n", ncfn);
    print!(" Number of linear convergence failures = {:4} \n", ncfl);
    let avdim = if nni > 0 {
        nli as f64 / nni as f64
    } else {
        ZERO
    };
    print!(
        " Average Krylov subspace dimension     = {} \n",
        fmt_f(avdim, 0, 3)
    );
    print!("\n\n--------------------------------------------------------------");
    print!("--------------\n");
    print!("--------------------------------------------------------------");
    print!("--------------\n");
}

/*
 This routine computes the right-hand side of the ODE system and
 returns it in cdot. The interaction rates are computed by calls to WebRates,
 and these are saved in fsave for use in preconditioning.
*/
fn f(t: f64, c: &NVector, cdot: &mut NVector, user_data: &mut UserData) -> i32 {
    let wdata = user_data
        .as_mut()
        .unwrap()
        .downcast_mut::<WebData>()
        .unwrap();
    let cdata = &c.data;
    let cdotdata = &mut cdot.data;

    let mxns = wdata.mxns;
    let ns = wdata.ns;
    let dx = wdata.dx;
    let dy = wdata.dy;
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
                cdotdata[ici] =
                    coy[i - 1] * (dcyui - dcyli) + cox[i - 1] * (dcxui - dcxli) + fsave[ici];
            }
        }
    }

    0
}

/*
  This routine computes the interaction rates for the species
  c_1, ... ,c_ns (stored in c[0],...,c[ns-1]), at one spatial point
  and at time t.
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

/* Set weights identically to the internal cvEwtSetSS for
   CVodeSStolerances(RTOL, ATOL) and keep the snapshot the C Precond
   would obtain from CVodeGetErrWeights. */
fn ewt(y: &NVector, w: &mut NVector, user_data: &mut UserData) -> i32 {
    let wdata = user_data
        .as_mut()
        .unwrap()
        .downcast_mut::<WebData>()
        .unwrap();
    for i in 0..NEQ {
        w.data[i] = ONE / (RTOL * SUNRabs(y.data[i]) + ATOL);
    }
    wdata.rewt.data.copy_from_slice(&w.data);
    0
}

/* SUNLinSol_SPGMRSetPrecType on the SPGMR solver attached to cvode_mem */
fn SUNLinSol_SPGMRSetPrecType(cvode_mem: &mut CVodeMem, pretype: i32) -> i32 {
    match &mut cvode_mem.cv_lmem {
        LsModule::Ls(cvls_mem) => cvls_mem.LS.set_prec_type(pretype),
        _ => SUN_ERR_ARG_OUTOFRANGE,
    }
}

/* SUNLinSol_SPGMRSetGSType on the SPGMR solver attached to cvode_mem */
fn SUNLinSol_SPGMRSetGSType(cvode_mem: &mut CVodeMem, gstype: i32) -> i32 {
    match &mut cvode_mem.cv_lmem {
        LsModule::Ls(cvls_mem) => cvls_mem.LS.set_gs_type(gstype),
        _ => SUN_ERR_ARG_OUTOFRANGE,
    }
}

/*
 This routine generates the block-diagonal part of the Jacobian
 corresponding to the interaction rates, multiplies by -gamma, adds
 the identity matrix, and calls SUNDlsMat_denseGETRF to do the LU decomposition of
 each diagonal block. The computation of the diagonal blocks uses
 the preset block and grouping information. One block per group is
 computed. The Jacobian elements are generated by difference
 quotients using calls to the routine fblock.

 This routine can be regarded as a prototype for the general case
 of a block-diagonal preconditioner. The blocks are of size mp, and
 there are ngrp=ngx*ngy blocks computed in the block-grouping scheme.
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
    let wdata = user_data
        .as_mut()
        .unwrap()
        .downcast_mut::<WebData>()
        .unwrap();

    /* C: retval = CVodeGetErrWeights(wdata->cvode_mem, rewt); here wdata.rewt
       is the snapshot the efun keeps identical to the integrator's ewt */
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

    /* Make mp calls to fblock to approximate each diagonal block of Jacobian.
       Here, fsave contains the base value of the rate vector and
       r0 is a minimum increment factor for the difference quotient.
       (The C code perturbs cdata in place and restores it; the local copy
       here sees the identical perturbed values.) */

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
  This routine computes one block of the interaction terms of the
  system, namely block (jx,jy), for use in preconditioning.
  Here jx and jy count from 0.
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
  This routine applies two inverse preconditioner matrices
  to the vector r, using the interaction-only block-diagonal Jacobian
  with block-grouping, denoted Jr, and Gauss-Seidel applied to the
  diffusion contribution to the Jacobian, denoted Jd.
  It first calls GSIter for a Gauss-Seidel approximation to
  ((I - gamma*Jd)-inverse)*r, and stores the result in z.
  Then it computes ((I - gamma*Jr)-inverse)*z, using LU factors of the
  blocks in P, and pivot information in pivot, and returns the result in z.
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
    let wdata = user_data
        .as_mut()
        .unwrap()
        .downcast_mut::<WebData>()
        .unwrap();

    N_VScale(ONE, r, z);

    /* call GSIter for Gauss-Seidel iterations (the temporary vector
       wdata->tmp is detached so z, x and wdata can be borrowed together) */

    let mut x = std::mem::take(&mut wdata.tmp);
    GSIter(gamma, z, &mut x, wdata);
    wdata.tmp = x;

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
  This routine performs ITMAX=5 Gauss-Seidel iterations to compute an
  approximation to (P-inverse)*z, where P = I - gamma*Jd, and
  Jd represents the diffusion contributions to the Jacobian.
  The answer is stored in z on return, and x is a temporary vector.
  The dimensions below assume a global constant NS >= ns.
  Some inner loops of length ns are implemented with the small
  vector kernels v_sum_prods, v_prod, v_inc_by_prod.
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
        let temp = ONE / (ONE + 2.0 * gamma * (cox[i] + coy[i]));
        beta[i] = gamma * cox[i] * temp;
        beta2[i] = 2.0 * beta[i];
        gam[i] = gamma * coy[i] * temp;
        gam2[i] = 2.0 * gam[i];
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

/* Check function return value (retval < 0 means failure) */
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
