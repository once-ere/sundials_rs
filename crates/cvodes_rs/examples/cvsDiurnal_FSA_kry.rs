/* -----------------------------------------------------------------
 * Translation of examples/cvodes/serial/cvsDiurnal_FSA_kry.c
 * (SUNDIALS 7.7.0).
 *
 * Example problem:
 *
 * An ODE system is generated from the following 2-species diurnal
 * kinetics advection-diffusion PDE system in 2 space dimensions:
 *
 * dc(i)/dt = Kh*(d/dx)^2 c(i) + V*dc(i)/dx + (d/dz)(Kv(z)*dc(i)/dz)
 *                 + Ri(c1,c2,t)      for i = 1,2,   where
 *   R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 ,
 *   R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 ,
 *   Kv(z) = Kv0*exp(z/5) ,
 * Kh, V, Kv0, q1, q2, and c3 are constants, and q3(t) and q4(t)
 * vary diurnally. The problem is solved with CVODES, with the
 * BDF/GMRES method and the block-diagonal part of the Newton
 * matrix as a left preconditioner.
 *
 * Optionally, CVODES can compute sensitivities with respect to the
 * problem parameters q1 and q2 (internal difference-quotient
 * sensitivity right-hand side).
 *
 * Execution:
 *    % cvsDiurnal_FSA_kry -nosensi
 *    % cvsDiurnal_FSA_kry -sensi sensi_meth err_con
 * where sensi_meth is one of {sim, stg, stg1} and err_con is one
 * of {t, f}.
 *
 * Translation note: the internal DQ sensitivity RHS requires the
 * pinned FSAUserData wrapper (ARCHITECTURE.md §3.6); the C UserData
 * fields other than p[] live in the DiurnalData payload behind
 * FSAUserData.user.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use cvodes_rs::cvodes::{
    CVode, CVodeCreate, CVodeFree, CVodeGetSens, CVodeInit, CVodeSStolerances,
    CVodeSensEEtolerances, CVodeSensInit1,
};
use cvodes_rs::cvodes_io::{
    CVodeGetLastOrder, CVodeGetLastStep, CVodeGetNumErrTestFails,
    CVodeGetNumNonlinSolvConvFails, CVodeGetNumNonlinSolvIters, CVodeGetNumRhsEvals,
    CVodeGetNumRhsEvalsSens, CVodeGetNumLinSolvSetups, CVodeGetNumSteps,
    CVodeGetSensNumErrTestFails, CVodeGetSensNumLinSolvSetups,
    CVodeGetSensNumNonlinSolvConvFails, CVodeGetSensNumNonlinSolvIters,
    CVodeGetSensNumRhsEvals, CVodeSetMaxNumSteps, CVodeSetSensDQMethod, CVodeSetSensErrCon,
    CVodeSetSensParams, CVodeSetUserData,
};
use cvodes_rs::cvodes_ls::{
    CVodeGetNumLinConvFails, CVodeGetNumLinIters, CVodeGetNumPrecEvals, CVodeGetNumPrecSolves,
    CVodeSetLinearSolver, CVodeSetPreconditioner,
};
use cvodes_rs::sundials_dense::{
    SUNDlsMat_denseAddIdentity, SUNDlsMat_denseCopy, SUNDlsMat_denseGETRF, SUNDlsMat_denseGETRS,
    SUNDlsMat_denseScale,
};
use cvodes_rs::sundials_utils::fmt_e;
use cvodes_rs::*;

/* Problem Constants */

const NUM_SPECIES: i64 = 2; /* number of species */
const C1_SCALE: f64 = 1.0e6; /* coefficients in initial profiles */
const C2_SCALE: f64 = 1.0e12;

const T0: f64 = 0.0; /* initial time */
const NOUT: i32 = 12; /* number of output times */
const TWOHR: f64 = 7200.0; /* number of seconds in two hours  */
const HALFDAY: f64 = 4.32e4; /* number of seconds in a half day */
#[allow(clippy::approx_constant)]
const PI: f64 = 3.1415926535898; /* pi */

const XMIN: f64 = 0.0; /* grid boundaries in x  */
const XMAX: f64 = 20.0;
const ZMIN: f64 = 30.0; /* grid boundaries in z  */
const ZMAX: f64 = 50.0;
const XMID: f64 = 10.0; /* grid midpoints in x,z */
const ZMID: f64 = 40.0;

const MX: i64 = 15; /* MX = number of x mesh points */
const MZ: i64 = 15; /* MZ = number of z mesh points */
const NSMX: i64 = NUM_SPECIES * MX; /* NSMX = NUM_SPECIES*MX */
const MM: i64 = MX * MZ; /* MM = MX*MZ */

/* CVodeInit Constants */
const RTOL: f64 = 1.0e-5; /* scalar relative tolerance */
const FLOOR: f64 = 100.0; /* value of C1 or C2 at which tolerances */
/* change from relative to absolute      */
const ATOL: f64 = RTOL * FLOOR; /* scalar absolute tolerance */
const NEQ: i64 = NUM_SPECIES * MM; /* NEQ = number of equations */

/* Sensitivity Constants */
const NP: usize = 8;
const NS: usize = 2;

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/* User-defined vector and matrix accessor macro: IJKth */

macro_rules! IJKth {
    ($vdata:expr, $i:expr, $j:expr, $k:expr) => {
        $vdata[($i - 1 + ($j) * NUM_SPECIES + ($k) * NSMX) as usize]
    };
}

/* Type : DiurnalData (the C UserData minus p[], which lives in the
   FSAUserData wrapper)
   contains preconditioner blocks, pivot arrays, and problem constants */

struct DiurnalData {
    P: Vec<Vec<DenseMatrix>>,   /* [MX][MZ] 2x2 blocks */
    Jbd: Vec<Vec<DenseMatrix>>, /* [MX][MZ] 2x2 blocks */
    pivot: Vec<Vec<[i64; NUM_SPECIES as usize]>>,
    q4: f64,
    om: f64,
    dx: f64,
    dz: f64,
    hdco: f64,
    haco: f64,
    vdco: f64,
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY CVODES
 *--------------------------------------------------------------------
 */

/*
 * f routine. Compute f(t,y).
 */
fn f(t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let fsa = user_data
        .as_mut()
        .unwrap()
        .downcast_mut::<FSAUserData>()
        .unwrap();

    /* Load problem coefficients and parameters */
    let Q1 = fsa.p[0];
    let Q2 = fsa.p[1];
    let C3 = fsa.p[2];
    let A3 = fsa.p[3];
    let A4 = fsa.p[4];

    let data = fsa.user.downcast_mut::<DiurnalData>().unwrap();
    let ydata = &y.data;
    let dydata = &mut ydot.data;

    /* Set diurnal rate coefficients. */

    let s = (data.om * t).sin();
    let q3: f64;
    if s > ZERO {
        q3 = (-A3 / s).exp();
        data.q4 = (-A4 / s).exp();
    } else {
        q3 = ZERO;
        data.q4 = ZERO;
    }

    /* Make local copies of problem variables, for efficiency. */

    let q4coef = data.q4;
    let delz = data.dz;
    let verdco = data.vdco;
    let hordco = data.hdco;
    let horaco = data.haco;

    /* Loop over all grid points. */

    for jz in 0..MZ {
        /* Set vertical diffusion coefficients at jz +- 1/2 */

        let zdn = ZMIN + (jz as f64 - 0.5) * delz;
        let zup = zdn + delz;
        let czdn = verdco * (0.2 * zdn).exp();
        let czup = verdco * (0.2 * zup).exp();
        let idn: i64 = if jz == 0 { 1 } else { -1 };
        let iup: i64 = if jz == MZ - 1 { -1 } else { 1 };
        for jx in 0..MX {
            /* Extract c1 and c2, and set kinetic rate terms. */

            let c1 = IJKth!(ydata, 1, jx, jz);
            let c2 = IJKth!(ydata, 2, jx, jz);
            let qq1 = Q1 * c1 * C3;
            let qq2 = Q2 * c1 * c2;
            let qq3 = q3 * C3;
            let qq4 = q4coef * c2;
            let rkin1 = -qq1 - qq2 + 2.0 * qq3 + qq4;
            let rkin2 = qq1 - qq2 - qq4;

            /* Set vertical diffusion terms. */

            let c1dn = IJKth!(ydata, 1, jx, jz + idn);
            let c2dn = IJKth!(ydata, 2, jx, jz + idn);
            let c1up = IJKth!(ydata, 1, jx, jz + iup);
            let c2up = IJKth!(ydata, 2, jx, jz + iup);
            let vertd1 = czup * (c1up - c1) - czdn * (c1 - c1dn);
            let vertd2 = czup * (c2up - c2) - czdn * (c2 - c2dn);

            /* Set horizontal diffusion and advection terms. */

            let ileft: i64 = if jx == 0 { 1 } else { -1 };
            let iright: i64 = if jx == MX - 1 { -1 } else { 1 };
            let c1lt = IJKth!(ydata, 1, jx + ileft, jz);
            let c2lt = IJKth!(ydata, 2, jx + ileft, jz);
            let c1rt = IJKth!(ydata, 1, jx + iright, jz);
            let c2rt = IJKth!(ydata, 2, jx + iright, jz);
            let hord1 = hordco * (c1rt - 2.0 * c1 + c1lt);
            let hord2 = hordco * (c2rt - 2.0 * c2 + c2lt);
            let horad1 = horaco * (c1rt - c1lt);
            let horad2 = horaco * (c2rt - c2lt);

            /* Load all terms into ydot. */

            IJKth!(dydata, 1, jx, jz) = vertd1 + hord1 + horad1 + rkin1;
            IJKth!(dydata, 2, jx, jz) = vertd2 + hord2 + horad2 + rkin2;
        }
    }

    0
}

/*
 * Preconditioner setup routine. Generate and preprocess P.
 */
fn Precond(
    _tn: f64,
    y: &NVector,
    _fy: &NVector,
    jok: bool,
    jcurPtr: &mut bool,
    gamma: f64,
    user_data: &mut UserData,
) -> i32 {
    /* Make local copies of pointers in user_data, and of pointer to y's data */
    let fsa = user_data
        .as_mut()
        .unwrap()
        .downcast_mut::<FSAUserData>()
        .unwrap();

    /* Load problem coefficients and parameters */
    let Q1 = fsa.p[0];
    let Q2 = fsa.p[1];
    let C3 = fsa.p[2];

    let data = fsa.user.downcast_mut::<DiurnalData>().unwrap();
    let ydata = &y.data;

    if jok {
        /* jok = SUNTRUE: Copy Jbd to P */

        for jz in 0..MZ as usize {
            for jx in 0..MX as usize {
                SUNDlsMat_denseCopy(&data.Jbd[jx][jz], &mut data.P[jx][jz]);
            }
        }

        *jcurPtr = SUNFALSE;
    } else {
        /* jok = SUNFALSE: Generate Jbd from scratch and copy to P */

        /* Make local copies of problem variables, for efficiency. */

        let q4coef = data.q4;
        let delz = data.dz;
        let verdco = data.vdco;
        let hordco = data.hdco;

        /* Compute 2x2 diagonal Jacobian blocks (using q4 values
         computed on the last f call).  Load into P. */

        for jz in 0..MZ {
            let zdn = ZMIN + (jz as f64 - 0.5) * delz;
            let zup = zdn + delz;
            let czdn = verdco * (0.2 * zdn).exp();
            let czup = verdco * (0.2 * zup).exp();
            let diag = -(czdn + czup + 2.0 * hordco);
            for jx in 0..MX {
                let c1 = IJKth!(ydata, 1, jx, jz);
                let c2 = IJKth!(ydata, 2, jx, jz);
                let j = &mut data.Jbd[jx as usize][jz as usize];
                /* IJth(a,i,j) = a[j-1][i-1] (column-major, 1-based) */
                j.set(0, 0, (-Q1 * C3 - Q2 * c2) + diag);
                j.set(0, 1, -Q2 * c1 + q4coef);
                j.set(1, 0, Q1 * C3 - Q2 * c2);
                j.set(1, 1, (-Q2 * c1 - q4coef) + diag);
                SUNDlsMat_denseCopy(
                    &data.Jbd[jx as usize][jz as usize],
                    &mut data.P[jx as usize][jz as usize],
                );
            }
        }

        *jcurPtr = SUNTRUE;
    }

    /* Scale by -gamma */

    for jz in 0..MZ as usize {
        for jx in 0..MX as usize {
            SUNDlsMat_denseScale(-gamma, &mut data.P[jx][jz]);
        }
    }

    /* Add identity matrix and do LU decompositions on blocks in place. */

    for jx in 0..MX as usize {
        for jz in 0..MZ as usize {
            SUNDlsMat_denseAddIdentity(&mut data.P[jx][jz]);
            let retval = SUNDlsMat_denseGETRF(&mut data.P[jx][jz], &mut data.pivot[jx][jz]);
            if retval != 0 {
                return 1;
            }
        }
    }

    0
}

/*
 * Preconditioner solve routine
 */
#[allow(clippy::too_many_arguments)]
fn PSolve(
    _tn: f64,
    _y: &NVector,
    _fy: &NVector,
    r: &NVector,
    z: &mut NVector,
    _gamma: f64,
    _delta: f64,
    _lr: i32,
    user_data: &mut UserData,
) -> i32 {
    /* Extract the P and pivot arrays from user_data. */
    let fsa = user_data
        .as_mut()
        .unwrap()
        .downcast_mut::<FSAUserData>()
        .unwrap();
    let data = fsa.user.downcast_mut::<DiurnalData>().unwrap();

    N_VScale(ONE, r, z);

    let zdata = &mut z.data;

    /* Solve the block-diagonal system Px = r using LU factors stored
       in P and pivot data in pivot, and return the solution in z. */

    for jx in 0..MX {
        for jz in 0..MZ {
            let k = (jx * NUM_SPECIES + jz * NSMX) as usize;
            let v = &mut zdata[k..k + NUM_SPECIES as usize];
            SUNDlsMat_denseGETRS(
                &data.P[jx as usize][jz as usize],
                &data.pivot[jx as usize][jz as usize],
                v,
            );
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
 * Process and verify arguments to cvsfwdkryx.
 */
fn process_args(argv: &[String]) -> (bool, i32, bool) {
    let sensi;
    let mut sensi_meth: i32 = -1;
    let mut err_con = false;

    if argv.len() < 2 {
        wrong_args(&argv[0]);
    }

    if argv[1] == "-nosensi" {
        sensi = false;
    } else if argv[1] == "-sensi" {
        sensi = true;
    } else {
        wrong_args(&argv[0]);
    }

    if sensi {
        if argv.len() != 4 {
            wrong_args(&argv[0]);
        }

        if argv[2] == "sim" {
            sensi_meth = CV_SIMULTANEOUS;
        } else if argv[2] == "stg" {
            sensi_meth = CV_STAGGERED;
        } else if argv[2] == "stg1" {
            sensi_meth = CV_STAGGERED1;
        } else {
            wrong_args(&argv[0]);
        }

        if argv[3] == "t" {
            err_con = true;
        } else if argv[3] == "f" {
            err_con = false;
        } else {
            wrong_args(&argv[0]);
        }
    }

    (sensi, sensi_meth, err_con)
}

fn wrong_args(name: &str) -> ! {
    println!("\nUsage: {} [-nosensi] [-sensi sensi_meth err_con]", name);
    println!("         sensi_meth = sim, stg, or stg1");
    println!("         err_con    = t or f");
    std::process::exit(0);
}

/*
 * Allocate memory for data structure of type DiurnalData
 */
fn AllocUserData() -> DiurnalData {
    let mut P = Vec::with_capacity(MX as usize);
    let mut Jbd = Vec::with_capacity(MX as usize);
    let mut pivot = Vec::with_capacity(MX as usize);
    for _jx in 0..MX {
        let mut prow = Vec::with_capacity(MZ as usize);
        let mut jrow = Vec::with_capacity(MZ as usize);
        let mut pvrow = Vec::with_capacity(MZ as usize);
        for _jz in 0..MZ {
            prow.push(DenseMatrix::new(NUM_SPECIES, NUM_SPECIES));
            jrow.push(DenseMatrix::new(NUM_SPECIES, NUM_SPECIES));
            pvrow.push([0i64; NUM_SPECIES as usize]);
        }
        P.push(prow);
        Jbd.push(jrow);
        pivot.push(pvrow);
    }
    DiurnalData {
        P,
        Jbd,
        pivot,
        q4: ZERO,
        om: ZERO,
        dx: ZERO,
        dz: ZERO,
        hdco: ZERO,
        haco: ZERO,
        vdco: ZERO,
    }
}

/*
 * Load problem constants in data (returns the problem parameters p)
 */
fn InitUserData(data: &mut DiurnalData) -> Vec<f64> {
    /* Set problem parameters */
    let Q1 = 1.63e-16; /* Q1  coefficients q1, q2, c3             */
    let Q2 = 4.66e-16; /* Q2                                      */
    let C3 = 3.7e16; /* C3                                      */
    let A3 = 22.62; /* A3  coefficient in expression for q3(t) */
    let A4 = 7.601; /* A4  coefficient in expression for q4(t) */
    let KH = 4.0e-6; /* KH  horizontal diffusivity Kh           */
    let VEL = 0.001; /* VEL advection velocity V                */
    let KV0 = 1.0e-8; /* KV0 coefficient in Kv(z)                */

    data.om = PI / HALFDAY;
    data.dx = (XMAX - XMIN) / (MX - 1) as f64;
    data.dz = (ZMAX - ZMIN) / (MZ - 1) as f64;
    data.hdco = KH / (data.dx * data.dx);
    data.haco = VEL / (2.0 * data.dx);
    data.vdco = (ONE / (data.dz * data.dz)) * KV0;

    let mut p = vec![ZERO; NP];
    p[0] = Q1;
    p[1] = Q2;
    p[2] = C3;
    p[3] = A3;
    p[4] = A4;
    p[5] = KH;
    p[6] = VEL;
    p[7] = KV0;

    p
}

/*
 * Set initial conditions in y
 */
fn SetInitialProfiles(y: &mut NVector, dx: f64, dz: f64) {
    /* Set pointer to data array in vector y. */

    let ydata = &mut y.data;

    /* Load initial profiles of c1 and c2 into y vector */

    for jz in 0..MZ {
        let z = ZMIN + jz as f64 * dz;
        let mut cz = (0.1 * (z - ZMID)) * (0.1 * (z - ZMID));
        cz = ONE - cz + 0.5 * (cz * cz);
        for jx in 0..MX {
            let x = XMIN + jx as f64 * dx;
            let mut cx = (0.1 * (x - XMID)) * (0.1 * (x - XMID));
            cx = ONE - cx + 0.5 * (cx * cx);
            IJKth!(ydata, 1, jx, jz) = C1_SCALE * cx * cz;
            IJKth!(ydata, 2, jx, jz) = C2_SCALE * cx * cz;
        }
    }
}

/*
 * Print current t, step count, order, stepsize, and sampled c1,c2 values
 */
fn PrintOutput(cvode_mem: &mut CVodeMem, t: f64, y: &NVector) {
    let mut nst: i64 = 0;
    let mut qu: i32 = 0;
    let mut hu: f64 = 0.0;

    let ydata = &y.data;

    CVodeGetNumSteps(cvode_mem, &mut nst);
    CVodeGetLastOrder(cvode_mem, &mut qu);
    CVodeGetLastStep(cvode_mem, &mut hu);

    println!("{} {:2}  {} {:5}", fmt_e(t, 8, 3), qu, fmt_e(hu, 8, 3), nst);

    print!("                                Solution       ");
    println!(
        "{} {} ",
        fmt_e(IJKth!(ydata, 1, 0, 0), 12, 4),
        fmt_e(IJKth!(ydata, 1, MX - 1, MZ - 1), 12, 4)
    );
    print!("                                               ");
    println!(
        "{} {} ",
        fmt_e(IJKth!(ydata, 2, 0, 0), 12, 4),
        fmt_e(IJKth!(ydata, 2, MX - 1, MZ - 1), 12, 4)
    );
}

/*
 * Print sampled sensitivities
 */
fn PrintOutputS(uS: &[NVector]) {
    let sdata = &uS[0].data;

    print!("                                ");
    println!("----------------------------------------");
    print!("                                Sensitivity 1  ");
    println!(
        "{} {} ",
        fmt_e(IJKth!(sdata, 1, 0, 0), 12, 4),
        fmt_e(IJKth!(sdata, 1, MX - 1, MZ - 1), 12, 4)
    );
    print!("                                               ");
    println!(
        "{} {} ",
        fmt_e(IJKth!(sdata, 2, 0, 0), 12, 4),
        fmt_e(IJKth!(sdata, 2, MX - 1, MZ - 1), 12, 4)
    );

    let sdata = &uS[1].data;

    print!("                                ");
    println!("----------------------------------------");
    print!("                                Sensitivity 2  ");
    println!(
        "{} {} ",
        fmt_e(IJKth!(sdata, 1, 0, 0), 12, 4),
        fmt_e(IJKth!(sdata, 1, MX - 1, MZ - 1), 12, 4)
    );
    print!("                                               ");
    println!(
        "{} {} ",
        fmt_e(IJKth!(sdata, 2, 0, 0), 12, 4),
        fmt_e(IJKth!(sdata, 2, MX - 1, MZ - 1), 12, 4)
    );
}

/*
 * Print final statistics contained in iopt
 */
fn PrintFinalStats(cvode_mem: &mut CVodeMem, sensi: bool, err_con: bool, sensi_meth: i32) {
    let mut nst: i64 = 0;
    let mut nfe: i64 = 0;
    let mut nsetups: i64 = 0;
    let mut nni: i64 = 0;
    let mut ncfn: i64 = 0;
    let mut netf: i64 = 0;
    let mut nfSe: i64 = 0;
    let mut nfeS: i64 = 0;
    let mut nsetupsS: i64 = 0;
    let mut nniS: i64 = 0;
    let mut ncfnS: i64 = 0;
    let mut netfS: i64 = 0;
    let mut nli: i64 = 0;
    let mut ncfl: i64 = 0;
    let mut npe: i64 = 0;
    let mut nps: i64 = 0;

    CVodeGetNumSteps(cvode_mem, &mut nst);
    CVodeGetNumRhsEvals(cvode_mem, &mut nfe);
    CVodeGetNumLinSolvSetups(cvode_mem, &mut nsetups);
    CVodeGetNumErrTestFails(cvode_mem, &mut netf);
    CVodeGetNumNonlinSolvIters(cvode_mem, &mut nni);
    CVodeGetNumNonlinSolvConvFails(cvode_mem, &mut ncfn);

    if sensi {
        CVodeGetSensNumRhsEvals(cvode_mem, &mut nfSe);
        CVodeGetNumRhsEvalsSens(cvode_mem, &mut nfeS);
        CVodeGetSensNumLinSolvSetups(cvode_mem, &mut nsetupsS);
        if err_con {
            CVodeGetSensNumErrTestFails(cvode_mem, &mut netfS);
        } else {
            netfS = 0;
        }
        if (sensi_meth == CV_STAGGERED) || (sensi_meth == CV_STAGGERED1) {
            CVodeGetSensNumNonlinSolvIters(cvode_mem, &mut nniS);
            CVodeGetSensNumNonlinSolvConvFails(cvode_mem, &mut ncfnS);
        } else {
            nniS = 0;
            ncfnS = 0;
        }
    }

    CVodeGetNumLinIters(cvode_mem, &mut nli);
    CVodeGetNumLinConvFails(cvode_mem, &mut ncfl);
    CVodeGetNumPrecEvals(cvode_mem, &mut npe);
    CVodeGetNumPrecSolves(cvode_mem, &mut nps);

    println!("\nFinal Statistics\n");
    println!("nst     = {:5}\n", nst);
    println!("nfe     = {:5}", nfe);
    println!("netf    = {:5}    nsetups  = {:5}", netf, nsetups);
    println!("nni     = {:5}    ncfn     = {:5}", nni, ncfn);

    if sensi {
        println!();
        println!("nfSe    = {:5}    nfeS     = {:5}", nfSe, nfeS);
        println!("netfs   = {:5}    nsetupsS = {:5}", netfS, nsetupsS);
        println!("nniS    = {:5}    ncfnS    = {:5}", nniS, ncfnS);
    }

    println!();
    println!("nli     = {:5}    ncfl     = {:5}", nli, ncfl);
    println!("npe     = {:5}    nps      = {:5}", npe, nps);
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
    /* Process arguments */
    let argv: Vec<String> = std::env::args().collect();
    let (sensi, sensi_meth, err_con) = process_args(&argv);

    /* Problem parameters */
    let mut data = AllocUserData();
    let p = InitUserData(&mut data);

    /* Create the SUNDIALS simulation context that all SUNDIALS objects require */
    let sunctx = SUNContext_Create();

    /* Initial states */
    let mut y = N_VNew_Serial(NEQ, &sunctx);
    SetInitialProfiles(&mut y, data.dx, data.dz);

    /* Tolerances */
    let abstol = ATOL;
    let reltol = RTOL;

    /* Create CVODES object */
    let mut cvode_mem = CVodeCreate(CV_BDF, &sunctx);

    let retval = CVodeSetUserData(
        &mut cvode_mem,
        Some(Box::new(FSAUserData {
            p: p.clone(),
            user: Box::new(data),
        })),
    );
    if check_retval(retval, "CVodeSetUserData") {
        return;
    }

    let retval = CVodeSetMaxNumSteps(&mut cvode_mem, 2000);
    if check_retval(retval, "CVodeSetMaxNumSteps") {
        return;
    }

    /* Allocate CVODES memory */
    let retval = CVodeInit(&mut cvode_mem, f, T0, &y);
    if check_retval(retval, "CVodeInit") {
        return;
    }

    let retval = CVodeSStolerances(&mut cvode_mem, reltol, abstol);
    if check_retval(retval, "CVodeSStolerances") {
        return;
    }

    /* Create the SUNLinSol_SPGMR linear solver with left
       preconditioning and the default Krylov dimension */
    let ls = SUNLinSol_SPGMR(&y, SUN_PREC_LEFT, 0, &sunctx);

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

    println!("\n2-species diurnal advection-diffusion problem");

    /* Forward sensitivity analysis */
    let mut uS: Vec<NVector> = Vec::new();
    if sensi {
        let plist: Vec<i32> = (0..NS as i32).collect();

        let mut pbar = [ZERO; NS];
        for is in 0..NS {
            pbar[is] = p[plist[is] as usize];
        }

        uS = (0..NS).map(|_| N_VClone(&y)).collect();
        for us in uS.iter_mut() {
            N_VConst(ZERO, us);
        }

        let retval = CVodeSensInit1(&mut cvode_mem, NS as i32, sensi_meth, None, &uS);
        if check_retval(retval, "CVodeSensInit") {
            return;
        }

        let retval = CVodeSensEEtolerances(&mut cvode_mem);
        if check_retval(retval, "CVodeSensEEtolerances") {
            return;
        }

        let retval = CVodeSetSensErrCon(&mut cvode_mem, err_con);
        if check_retval(retval, "CVodeSetSensErrCon") {
            return;
        }

        let retval = CVodeSetSensDQMethod(&mut cvode_mem, CV_CENTERED, ZERO);
        if check_retval(retval, "CVodeSetSensDQMethod") {
            return;
        }

        let retval = CVodeSetSensParams(&mut cvode_mem, Some(&p), Some(&pbar), Some(&plist));
        if check_retval(retval, "CVodeSetSensParams") {
            return;
        }

        print!("Sensitivity: YES ");
        if sensi_meth == CV_SIMULTANEOUS {
            print!("( SIMULTANEOUS +");
        } else if sensi_meth == CV_STAGGERED {
            print!("( STAGGERED +");
        } else {
            print!("( STAGGERED1 +");
        }
        if err_con {
            print!(" FULL ERROR CONTROL )");
        } else {
            print!(" PARTIAL ERROR CONTROL )");
        }
    } else {
        print!("Sensitivity: NO ");
    }

    /* In loop over output points, call CVode, print results, test for error */

    println!("\n");
    println!("========================================================================");
    println!("     T     Q       H      NST                    Bottom left  Top right ");
    println!("========================================================================");

    let mut t = T0;
    let mut tout = TWOHR;
    for _iout in 1..=NOUT {
        let retval = CVode(&mut cvode_mem, tout, &mut y, &mut t, CV_NORMAL);
        if check_retval(retval, "CVode") {
            break;
        }
        PrintOutput(&mut cvode_mem, t, &y);
        if sensi {
            let retval = CVodeGetSens(&cvode_mem, &mut t, &mut uS);
            if check_retval(retval, "CVodeGetSens") {
                break;
            }
            PrintOutputS(&uS);
        }

        println!("------------------------------------------------------------------------");

        tout += TWOHR;
    }

    /* Print final statistics */
    PrintFinalStats(&mut cvode_mem, sensi, err_con, sensi_meth);

    /* Free memory */
    drop(y);
    drop(uS);
    CVodeFree(cvode_mem);
}
