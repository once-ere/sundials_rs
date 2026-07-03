/* -----------------------------------------------------------------
 * Translated from examples/cvode/serial/cvKrylovDemo_ls.c (CVODE 7.7.0)
 *
 * This example loops through the available iterative linear solvers:
 * SPGMR, SPFGMR, SPBCG and SPTFQMR.
 *
 * Example problem:
 *
 * An ODE system is generated from the following 2-species diurnal
 * kinetics advection-diffusion PDE system in 2 space dimensions:
 *
 * dc(i)/dt = Kh*(d/dx)^2 c(i) + V*dc(i)/dx + (d/dy)(Kv(y)*dc(i)/dy)
 *                 + Ri(c1,c2,t)      for i = 1,2,   where
 *   R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 ,
 *   R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 ,
 *   Kv(y) = Kv0*exp(y/5) ,
 * Kh, V, Kv0, q1, q2, and c3 are constants, and q3(t) and q4(t)
 * vary diurnally. The problem is posed on the square
 *   0 <= x <= 20,    30 <= y <= 50   (all in km),
 * with homogeneous Neumann boundary conditions, and for time t in
 *   0 <= t <= 86400 sec (1 day).
 * The PDE system is treated by central differences on a uniform
 * 10 x 10 mesh, with simple polynomial initial profiles.
 * The problem is solved with CVODE, with the BDF/GMRES, BDF/FGMRES
 * BDF/Bi-CGStab, and BDF/TFQMR methods (i.e. using the SUNLinSol_SPGMR,
 * SUNLinSol_SPFGMR, SUNLinSol_SPBCGS, and SUNLinSol_SPTFQMR linear solvers)
 * and the block-diagonal part of the Newton matrix as a left preconditioner.
 * A copy of the block-diagonal part of the Jacobian is saved and
 * conditionally reused within the Precond routine.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use cvode_rs::sundials_dense::{
    SUNDlsMat_denseAddIdentity, SUNDlsMat_denseCopy, SUNDlsMat_denseGETRF, SUNDlsMat_denseGETRS,
    SUNDlsMat_denseScale,
};
use cvode_rs::sundials_utils::fmt_e;
use cvode_rs::*;

/* Problem Constants */

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;

const NUM_SPECIES: i64 = 2; /* number of species         */
const KH: f64 = 4.0e-6; /* horizontal diffusivity Kh */
const VEL: f64 = 0.001; /* advection velocity V      */
const KV0: f64 = 1.0e-8; /* coefficient in Kv(y)      */
const Q1: f64 = 1.63e-16; /* coefficients q1, q2, c3   */
const Q2: f64 = 4.66e-16;
const C3: f64 = 3.7e16;
const A3: f64 = 22.62; /* coefficient in expression for q3(t) */
const A4: f64 = 7.601; /* coefficient in expression for q4(t) */
const C1_SCALE: f64 = 1.0e6; /* coefficients in initial profiles    */
const C2_SCALE: f64 = 1.0e12;

const T0: f64 = ZERO; /* initial time */
const NOUT: i32 = 12; /* number of output times */
const TWOHR: f64 = 7200.0; /* number of seconds in two hours  */
const HALFDAY: f64 = 4.32e4; /* number of seconds in a half day */
const PI: f64 = 3.1415926535898; /* pi */

const XMIN: f64 = ZERO; /* grid boundaries in x  */
const XMAX: f64 = 20.0;
const YMIN: f64 = 30.0; /* grid boundaries in y  */
const YMAX: f64 = 50.0;
const XMID: f64 = 10.0; /* grid midpoints in x,y */
const YMID: f64 = 40.0;

const MX: i64 = 10; /* MX = number of x mesh points */
const MY: i64 = 10; /* MY = number of y mesh points */
const NSMX: i64 = 20; /* NSMX = NUM_SPECIES*MX */
#[allow(dead_code)]
const MM: i64 = MX * MY; /* MM = MX*MY */

/* CVodeInit Constants */

const RTOL: f64 = 1.0e-5; /* scalar relative tolerance */
const FLOOR: f64 = 100.0; /* value of C1 or C2 at which tolerances */
/* change from relative to absolute      */
const ATOL: f64 = RTOL * FLOOR; /* scalar absolute tolerance */
const NEQ: i64 = NUM_SPECIES * MM; /* NEQ = number of equations */

/* Linear Solver Loop Constants */

const USE_SPGMR: i32 = 0;
const USE_SPFGMR: i32 = 1;
const USE_SPBCG: i32 = 2;
const USE_SPTFQMR: i32 = 3;

/* User-defined vector and matrix accessor macro: IJKth */

/* IJKth(vdata,i,j,k) references the element in the vdata array for
   species i at mesh point (j,k), where 1 <= i <= NUM_SPECIES,
   0 <= j <= MX-1, 0 <= k <= MY-1.
   For each mesh point (j,k), the elements for species i and i+1 are
   contiguous within vdata. */

macro_rules! IJKth {
    ($vdata:expr, $i:expr, $j:expr, $k:expr) => {
        $vdata[($i - 1 + ($j) * NUM_SPECIES + ($k) * NSMX) as usize]
    };
}

/* Type : UserDataStruct
   contains preconditioner blocks, pivot arrays, problem constants,
   and linsolver type */

struct UserDataStruct {
    P: Vec<Vec<DenseMatrix>>,   /* [MX][MY] 2x2 blocks */
    Jbd: Vec<Vec<DenseMatrix>>, /* [MX][MY] 2x2 blocks */
    pivot: Vec<Vec<[i64; NUM_SPECIES as usize]>>,
    q4: f64,
    om: f64,
    dx: f64,
    dy: f64,
    hdco: f64,
    haco: f64,
    vdco: f64,
    #[allow(dead_code)]
    linsolver: i32,
}

/*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 */

/* Allocate memory for data structure of type UserDataStruct */

fn AllocUserData() -> UserDataStruct {
    let mut P = Vec::with_capacity(MX as usize);
    let mut Jbd = Vec::with_capacity(MX as usize);
    let mut pivot = Vec::with_capacity(MX as usize);
    for _jx in 0..MX {
        let mut prow = Vec::with_capacity(MY as usize);
        let mut jrow = Vec::with_capacity(MY as usize);
        let mut pvrow = Vec::with_capacity(MY as usize);
        for _jy in 0..MY {
            prow.push(DenseMatrix::new(NUM_SPECIES, NUM_SPECIES));
            jrow.push(DenseMatrix::new(NUM_SPECIES, NUM_SPECIES));
            pvrow.push([0i64; NUM_SPECIES as usize]);
        }
        P.push(prow);
        Jbd.push(jrow);
        pivot.push(pvrow);
    }
    UserDataStruct {
        P,
        Jbd,
        pivot,
        q4: ZERO,
        om: ZERO,
        dx: ZERO,
        dy: ZERO,
        hdco: ZERO,
        haco: ZERO,
        vdco: ZERO,
        linsolver: 0,
    }
}

/* Load problem constants in data */

fn InitUserData(data: &mut UserDataStruct) {
    data.om = PI / HALFDAY;
    data.dx = (XMAX - XMIN) / (MX - 1) as f64;
    data.dy = (YMAX - YMIN) / (MY - 1) as f64;
    data.hdco = KH / (data.dx * data.dx);
    data.haco = VEL / (TWO * data.dx);
    data.vdco = (ONE / (data.dy * data.dy)) * KV0;
}

/* Set initial conditions in u */

fn SetInitialProfiles(u: &mut NVector, dx: f64, dy: f64) {
    /* Set pointer to data array in vector u. */

    let udata = &mut u.data;

    /* Load initial profiles of c1 and c2 into u vector */

    for jy in 0..MY {
        let y = YMIN + jy as f64 * dy;
        let mut cy = (0.1 * (y - YMID)) * (0.1 * (y - YMID));
        cy = ONE - cy + 0.5 * (cy * cy);
        for jx in 0..MX {
            let x = XMIN + jx as f64 * dx;
            let mut cx = (0.1 * (x - XMID)) * (0.1 * (x - XMID));
            cx = ONE - cx + 0.5 * (cx * cx);
            IJKth!(udata, 1, jx, jy) = C1_SCALE * cx * cy;
            IJKth!(udata, 2, jx, jy) = C2_SCALE * cx * cy;
        }
    }
}

/* Print current t, step count, order, stepsize, and sampled c1,c2 values */

fn PrintOutput(cvode_mem: &mut CVodeMem, u: &NVector, t: f64) {
    let mut nst = 0i64;
    let mut qu = 0i32;
    let mut hu = ZERO;
    let mxh = MX / 2 - 1;
    let myh = MY / 2 - 1;
    let mx1 = MX - 1;
    let my1 = MY - 1;

    let udata = &u.data;

    let mut retval = CVodeGetNumSteps(cvode_mem, &mut nst);
    check_retval(retval, "CVodeGetNumSteps");
    retval = CVodeGetLastOrder(cvode_mem, &mut qu);
    check_retval(retval, "CVodeGetLastOrder");
    retval = CVodeGetLastStep(cvode_mem, &mut hu);
    check_retval(retval, "CVodeGetLastStep");

    println!(
        "t = {}   no. steps = {}   order = {}   stepsize = {}",
        fmt_e(t, 0, 2),
        nst,
        qu,
        fmt_e(hu, 0, 2)
    );
    println!(
        "c1 (bot.left/middle/top rt.) = {}  {}  {}",
        fmt_e(IJKth!(udata, 1, 0, 0), 12, 3),
        fmt_e(IJKth!(udata, 1, mxh, myh), 12, 3),
        fmt_e(IJKth!(udata, 1, mx1, my1), 12, 3)
    );
    println!(
        "c2 (bot.left/middle/top rt.) = {}  {}  {}\n",
        fmt_e(IJKth!(udata, 2, 0, 0), 12, 3),
        fmt_e(IJKth!(udata, 2, mxh, myh), 12, 3),
        fmt_e(IJKth!(udata, 2, mx1, my1), 12, 3)
    );
}

/* Get and print final statistics */

fn PrintStats(cvode_mem: &mut CVodeMem, linsolver: i32, stats: i32) {
    let (mut lenrw, mut leniw) = (0i64, 0i64);
    let (mut lenrwLS, mut leniwLS) = (0i64, 0i64);
    let (mut nst, mut nfe, mut nsetups, mut nni, mut ncfn, mut netf) =
        (0i64, 0i64, 0i64, 0i64, 0i64, 0i64);
    let (mut nje, mut nli, mut npe, mut nps, mut ncfl, mut nfeLS, mut njts, mut njte) =
        (0i64, 0i64, 0i64, 0i64, 0i64, 0i64, 0i64, 0i64);

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
    CVodeGetLinSolveStats(
        cvode_mem, &mut nje, &mut nfeLS, &mut nli, &mut ncfl, &mut npe, &mut nps, &mut njts,
        &mut njte,
    );
    check_retval(retval, "CVodeGetLinWorkSpace");

    /* nje, njts and njte are retrieved but not printed (as in the C) */
    let _ = (nje, njts, njte);

    if stats != 0 {
        println!("\nFinal Statistics.. \n");
    } else {
        println!("\nIntermediate Statistics.. \n");
    }
    println!("lenrw   = {:5}     leniw   = {:5}", lenrw, leniw);
    println!("lenrwLS = {:5}     leniwLS = {:5}", lenrwLS, leniwLS);
    println!("nst     = {:5}", nst);
    println!("nfe     = {:5}     nfeLS   = {:5}", nfe, nfeLS);
    println!("nni     = {:5}     nli     = {:5}", nni, nli);
    println!("nsetups = {:5}     netf    = {:5}", nsetups, netf);
    println!("npe     = {:5}     nps     = {:5}", npe, nps);
    println!("ncfn    = {:5}     ncfl    = {:5}\n", ncfn, ncfl);

    if linsolver < 2 {
        println!("======================================================================\n");
    }
}

/* Check function return value */
fn check_retval(retval: i32, funcname: &str) -> bool {
    if retval < 0 {
        eprintln!("\nSUNDIALS_ERROR: {}() failed with retval = {}\n", funcname, retval);
        return true;
    }
    false
}

/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/* f routine. Compute RHS function f(t,u). */

fn f(t: f64, u: &NVector, udot: &mut NVector, user_data: &mut UserData) -> i32 {
    let data = user_data
        .as_mut()
        .unwrap()
        .downcast_mut::<UserDataStruct>()
        .unwrap();
    let udata = &u.data;
    let dudata = &mut udot.data;

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
    let dely = data.dy;
    let verdco = data.vdco;
    let hordco = data.hdco;
    let horaco = data.haco;

    /* Loop over all grid points. */

    for jy in 0..MY {
        /* Set vertical diffusion coefficients at jy +- 1/2 */

        let ydn = YMIN + (jy as f64 - 0.5) * dely;
        let yup = ydn + dely;
        let cydn = verdco * (0.2 * ydn).exp();
        let cyup = verdco * (0.2 * yup).exp();
        let idn: i64 = if jy == 0 { 1 } else { -1 };
        let iup: i64 = if jy == MY - 1 { -1 } else { 1 };
        for jx in 0..MX {
            /* Extract c1 and c2, and set kinetic rate terms. */

            let c1 = IJKth!(udata, 1, jx, jy);
            let c2 = IJKth!(udata, 2, jx, jy);
            let qq1 = Q1 * c1 * C3;
            let qq2 = Q2 * c1 * c2;
            let qq3 = q3 * C3;
            let qq4 = q4coef * c2;
            let rkin1 = -qq1 - qq2 + TWO * qq3 + qq4;
            let rkin2 = qq1 - qq2 - qq4;

            /* Set vertical diffusion terms. */

            let c1dn = IJKth!(udata, 1, jx, jy + idn);
            let c2dn = IJKth!(udata, 2, jx, jy + idn);
            let c1up = IJKth!(udata, 1, jx, jy + iup);
            let c2up = IJKth!(udata, 2, jx, jy + iup);
            let vertd1 = cyup * (c1up - c1) - cydn * (c1 - c1dn);
            let vertd2 = cyup * (c2up - c2) - cydn * (c2 - c2dn);

            /* Set horizontal diffusion and advection terms. */

            let ileft: i64 = if jx == 0 { 1 } else { -1 };
            let iright: i64 = if jx == MX - 1 { -1 } else { 1 };
            let c1lt = IJKth!(udata, 1, jx + ileft, jy);
            let c2lt = IJKth!(udata, 2, jx + ileft, jy);
            let c1rt = IJKth!(udata, 1, jx + iright, jy);
            let c2rt = IJKth!(udata, 2, jx + iright, jy);
            let hord1 = hordco * (c1rt - TWO * c1 + c1lt);
            let hord2 = hordco * (c2rt - TWO * c2 + c2lt);
            let horad1 = horaco * (c1rt - c1lt);
            let horad2 = horaco * (c2rt - c2lt);

            /* Load all terms into udot. */

            IJKth!(dudata, 1, jx, jy) = vertd1 + hord1 + horad1 + rkin1;
            IJKth!(dudata, 2, jx, jy) = vertd2 + hord2 + horad2 + rkin2;
        }
    }

    0
}

/* Preconditioner setup routine. Generate and preprocess P. */

fn Precond(
    _tn: f64,
    u: &NVector,
    _fu: &NVector,
    jok: bool,
    jcurPtr: &mut bool,
    gamma: f64,
    user_data: &mut UserData,
) -> i32 {
    /* Make local copies of pointers in user_data, and of pointer to u's data */

    let data = user_data
        .as_mut()
        .unwrap()
        .downcast_mut::<UserDataStruct>()
        .unwrap();
    let udata = &u.data;

    if jok {
        /* jok = SUNTRUE: Copy Jbd to P */

        for jy in 0..MY as usize {
            for jx in 0..MX as usize {
                SUNDlsMat_denseCopy(&data.Jbd[jx][jy], &mut data.P[jx][jy]);
            }
        }

        *jcurPtr = SUNFALSE;
    } else {
        /* jok = SUNFALSE: Generate Jbd from scratch and copy to P */

        /* Make local copies of problem variables, for efficiency. */

        let q4coef = data.q4;
        let dely = data.dy;
        let verdco = data.vdco;
        let hordco = data.hdco;

        /* Compute 2x2 diagonal Jacobian blocks (using q4 values
        computed on the last f call).  Load into P. */

        for jy in 0..MY {
            let ydn = YMIN + (jy as f64 - 0.5) * dely;
            let yup = ydn + dely;
            let cydn = verdco * (0.2 * ydn).exp();
            let cyup = verdco * (0.2 * yup).exp();
            let diag = -(cydn + cyup + TWO * hordco);
            for jx in 0..MX {
                let c1 = IJKth!(udata, 1, jx, jy);
                let c2 = IJKth!(udata, 2, jx, jy);
                let j = &mut data.Jbd[jx as usize][jy as usize];
                /* IJth(a,i,j) = a[j-1][i-1] (column-major, 1-based) */
                j.set(0, 0, (-Q1 * C3 - Q2 * c2) + diag);
                j.set(0, 1, -Q2 * c1 + q4coef);
                j.set(1, 0, Q1 * C3 - Q2 * c2);
                j.set(1, 1, (-Q2 * c1 - q4coef) + diag);
                SUNDlsMat_denseCopy(
                    &data.Jbd[jx as usize][jy as usize],
                    &mut data.P[jx as usize][jy as usize],
                );
            }
        }

        *jcurPtr = SUNTRUE;
    }

    /* Scale by -gamma */

    for jy in 0..MY as usize {
        for jx in 0..MX as usize {
            SUNDlsMat_denseScale(-gamma, &mut data.P[jx][jy]);
        }
    }

    /* Add identity matrix and do LU decompositions on blocks in place. */

    for jx in 0..MX as usize {
        for jy in 0..MY as usize {
            SUNDlsMat_denseAddIdentity(&mut data.P[jx][jy]);
            let retval = SUNDlsMat_denseGETRF(&mut data.P[jx][jy], &mut data.pivot[jx][jy]);
            if retval != 0 {
                return 1;
            }
        }
    }

    0
}

/* Preconditioner solve routine */

fn PSolve(
    _tn: f64,
    _u: &NVector,
    _fu: &NVector,
    r: &NVector,
    z: &mut NVector,
    _gamma: f64,
    _delta: f64,
    _lr: i32,
    user_data: &mut UserData,
) -> i32 {
    /* Extract the P and pivot arrays from user_data. */

    let data = user_data
        .as_mut()
        .unwrap()
        .downcast_mut::<UserDataStruct>()
        .unwrap();

    N_VScale(ONE, r, z);

    let zdata = &mut z.data;

    /* Solve the block-diagonal system Px = r using LU factors stored
    in P and pivot data in pivot, and return the solution in z. */

    for jx in 0..MX {
        for jy in 0..MY {
            let k = (jx * NUM_SPECIES + jy * NSMX) as usize;
            let v = &mut zdata[k..k + NUM_SPECIES as usize];
            SUNDlsMat_denseGETRS(
                &data.P[jx as usize][jy as usize],
                &data.pivot[jx as usize][jy as usize],
                v,
            );
        }
    }

    0
}

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

fn main() {
    let mut nrmfactor = 0; /* LS norm conversion factor flag */
    let mut monitor = 0; /* LS residual monitoring flag    */

    /* Retrieve the command-line options */
    let args: Vec<String> = std::env::args().collect();
    if args.len() > 1 {
        nrmfactor = args[1].parse::<i32>().unwrap_or(0);
    }
    if args.len() > 2 {
        monitor = args[2].parse::<i32>().unwrap_or(0);
    }

    /* Create the SUNDIALS context (the SUNLogger used for residual
     * monitoring in the C original is not ported) */
    let sunctx = SUNContext_Create();

    /* Allocate memory, and set problem data, initial values, tolerances */
    let mut u = N_VNew_Serial(NEQ, &sunctx);
    let mut data = AllocUserData();
    InitUserData(&mut data);
    SetInitialProfiles(&mut u, data.dx, data.dy);
    let abstol = ATOL;
    let reltol = RTOL;

    /* Call CVodeCreate to create the solver memory and specify the
     * Backward Differentiation Formula */
    let mut cvode_mem = CVodeCreate(CV_BDF, &sunctx);

    /* Set the pointer to user-defined data */
    let mut retval = CVodeSetUserData(&mut cvode_mem, Some(Box::new(data)));
    if check_retval(retval, "CVodeSetUserData") {
        std::process::exit(1);
    }

    /* Call CVodeInit to initialize the integrator memory and specify the
     * user's right hand side function in u'=f(t,u), the initial time T0, and
     * the initial dependent variable vector u. */
    retval = CVodeInit(&mut cvode_mem, f, T0, &u);
    if check_retval(retval, "CVodeInit") {
        std::process::exit(1);
    }

    /* Call CVodeSStolerances to specify the scalar relative tolerance
     * and scalar absolute tolerances */
    retval = CVodeSStolerances(&mut cvode_mem, reltol, abstol);
    if check_retval(retval, "CVodeSStolerances") {
        std::process::exit(1);
    }

    /* Create the SUNNonlinearSolver */
    let nls = SUNNonlinSol_Newton(&u, &sunctx);

    /* Call CVodeSetNonlinearSolver to attach the nonlinear solver to CVode */
    retval = CVodeSetNonlinearSolver(&mut cvode_mem, nls);
    if check_retval(retval, "CVodeSetNonlinearSolver") {
        std::process::exit(1);
    }

    let mut t = ZERO;

    /* START: Loop through SPGMR, SPFGMR, SPBCG and SPTFQMR linear solver modules */
    for linsolver in 0..4 {
        if linsolver != 0 {
            /* Re-initialize user data */
            let data = CVodeGetUserData(&mut cvode_mem)
                .as_mut()
                .unwrap()
                .downcast_mut::<UserDataStruct>()
                .unwrap();
            InitUserData(data);
            let (dx, dy) = (data.dx, data.dy);
            SetInitialProfiles(&mut u, dx, dy);

            /* Re-initialize CVode for the solution of the same problem, but
            using a different linear solver module */
            retval = CVodeReInit(&mut cvode_mem, T0, &u);
            if check_retval(retval, "CVodeReInit") {
                std::process::exit(1);
            }
        }

        /* The previous linear solver is freed automatically when a new
        linear solver module is attached (SUNLinSolFree in the C) */

        /* Set the linear solver type in user data */
        CVodeGetUserData(&mut cvode_mem)
            .as_mut()
            .unwrap()
            .downcast_mut::<UserDataStruct>()
            .unwrap()
            .linsolver = linsolver;

        match linsolver {
            /* (a) SPGMR */
            USE_SPGMR => {
                /* Print header */
                print!(" -------");
                print!(" \n| SPGMR |\n");
                print!(" -------\n");

                /* Call SUNLinSol_SPGMR to specify the linear solver SPGMR with
                left preconditioning and the default maximum Krylov dimension */
                let LS = SUNLinSol_SPGMR(&u, SUN_PREC_LEFT, 0, &sunctx);

                retval = CVodeSetLinearSolver(&mut cvode_mem, LS, None);
                if check_retval(retval, "CVodeSetLinearSolver") {
                    std::process::exit(1);
                }
            }

            /* (b) SPFGMR */
            USE_SPFGMR => {
                /* Print header */
                print!(" ---------");
                print!(" \n| SPFGMR |\n");
                print!(" ---------\n");

                /* Call SUNLinSol_SPFGMR to specify the linear solver SPFGMR with
                left preconditioning and the default maximum Krylov dimension */
                let LS = SUNLinSol_SPFGMR(&u, SUN_PREC_LEFT, 0, &sunctx);

                retval = CVodeSetLinearSolver(&mut cvode_mem, LS, None);
                if check_retval(retval, "CVodeSetLinearSolver") {
                    std::process::exit(1);
                }
            }

            /* (c) SPBCG */
            USE_SPBCG => {
                /* Print header */
                print!(" -------");
                print!(" \n| SPBCGS |\n");
                print!(" -------\n");

                /* Call SUNLinSol_SPBCGS to specify the linear solver SPBCGS with
                left preconditioning and the default maximum Krylov dimension */
                let LS = SUNLinSol_SPBCGS(&u, SUN_PREC_LEFT, 0, &sunctx);

                retval = CVodeSetLinearSolver(&mut cvode_mem, LS, None);
                if check_retval(retval, "CVodeSetLinearSolver") {
                    std::process::exit(1);
                }
            }

            /* (d) SPTFQMR */
            USE_SPTFQMR => {
                /* Print header */
                print!(" ---------");
                print!(" \n| SPTFQMR |\n");
                print!(" ---------\n");

                /* Call SUNLinSol_SPTFQMR to specify the linear solver SPTFQMR with
                left preconditioning and the default maximum Krylov dimension */
                let LS = SUNLinSol_SPTFQMR(&u, SUN_PREC_LEFT, 0, &sunctx);

                retval = CVodeSetLinearSolver(&mut cvode_mem, LS, None);
                if check_retval(retval, "CVodeSetLinearSolver") {
                    std::process::exit(1);
                }
            }

            _ => {}
        }

        /* Set preconditioner setup and solve routines Precond and PSolve,
        and the pointer to the user-defined block data */
        retval = CVodeSetPreconditioner(&mut cvode_mem, Some(Precond), Some(PSolve));
        if check_retval(retval, "CVodeSetPreconditioner") {
            std::process::exit(1);
        }

        /* Set the linear solver tolerance conversion factor */
        let nrmfac = match nrmfactor {
            1 => {
                /* use the square root of the vector length */
                (NEQ as f64).sqrt()
            }
            2 => {
                /* compute with dot product */
                -ONE
            }
            _ => {
                /* use the default */
                ZERO
            }
        };

        retval = CVodeSetLSNormFactor(&mut cvode_mem, nrmfac);
        if check_retval(retval, "CVodeSetLSNormFactor") {
            std::process::exit(1);
        }

        /* In loop over output points, call CVode, print results, and test for error */
        println!(" \n2-species diurnal advection-diffusion problem\n");
        let mut tout = TWOHR;
        for _iout in 1..=NOUT {
            retval = CVode(&mut cvode_mem, tout, &mut u, &mut t, CV_NORMAL);
            if monitor == 0 {
                PrintOutput(&mut cvode_mem, &u, t);
            }
            if check_retval(retval, "CVode") {
                break;
            }
            tout += TWOHR;
        }
        if monitor != 0 {
            PrintOutput(&mut cvode_mem, &u, t);
        }
        PrintStats(&mut cvode_mem, linsolver, 1);
    } /* END: Loop through SPGMR, SPFGMR, SPBCG and SPTFQMR linear solver modules */

    /* Free memory */
    CVodeFree(cvode_mem);
}
