/* -----------------------------------------------------------------
 * Translated from examples/ida/serial/idaFoodWeb_bnd.c (IDA 7.7.0)
 * Programmer(s): Allan Taylor, Alan Hindmarsh and Radu Serban @ LLNL
 *
 * Example program for IDA: Food web problem (serial, banded linear
 * solver, IDACalcIC for initial condition calculation).
 *
 * The mathematical problem is a DAE system arising from a
 * predator-prey food-web PDE model with diffusion on the unit square,
 * discretized by central differencing on an MX by MY mesh. np = 1,
 * ns = 2. Homogeneous Neumann boundary conditions.
 *
 * Output is printed at t = 0, .001, .01, .1, .4, .7, 1.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use ida_rs::sundials_utils::{fmt_e, fmt_g};
use ida_rs::*;

/* Problem Constants. */
const NPREY: usize = 1;
const NUM_SPECIES: usize = 2 * NPREY;

const PI: f64 = 3.1415926535898;
const FOURPI: f64 = 4.0 * PI;

const MX: usize = 20;
const MY: usize = 20;
const NSMX: usize = NUM_SPECIES * MX;
const NEQ: usize = NUM_SPECIES * MX * MY;
const AA: f64 = 1.0;
const EE: f64 = 10000.0;
const GG: f64 = 0.5e-6;
const BB: f64 = 1.0;
const DPREY: f64 = 1.0;
const DPRED: f64 = 0.05;
const ALPHA: f64 = 50.0;
const BETA: f64 = 1000.0;
const AX: f64 = 1.0;
const AY: f64 = 1.0;
const RTOL: f64 = 1.0e-5;
const ATOL: f64 = 1.0e-5;
const NOUT: i32 = 6;
const TMULT: f64 = 10.0;
const TADD: f64 = 0.3;
const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/* IJ_Vptr(vv,i,j): index into vv for species is=0, x-index i, y-index j. */
fn ij(i: usize, j: usize) -> usize {
    i * NUM_SPECIES + j * NSMX
}

/* Type: WebData.  Contains problem constants, etc. */
struct WebData {
    np: usize,
    dx: f64,
    dy: f64,
    acoef: [[f64; NUM_SPECIES]; NUM_SPECIES],
    cox: [f64; NUM_SPECIES],
    coy: [f64; NUM_SPECIES],
    bcoef: [f64; NUM_SPECIES],
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY IDA
 *--------------------------------------------------------------------
 */

/*
 * resweb: System residual function for predator-prey system.
 */
fn resweb(_tt: f64, cc: &NVector, cp: &NVector, res: &mut NVector, user_data: &mut UserData) -> i32 {
    let webdata = user_data.as_ref().unwrap().downcast_ref::<WebData>().unwrap();
    let np = webdata.np;

    /* Call Fweb to set res to vector of right-hand sides. */
    Fweb(cc, res, webdata);

    /* Loop over all grid points, setting residual values appropriately
       for differential or algebraic components. */
    for jy in 0..MY {
        let yloc = NSMX * jy;
        for jx in 0..MX {
            let loc = yloc + NUM_SPECIES * jx;
            for is in 0..NUM_SPECIES {
                if is < np {
                    res.data[loc + is] = cp.data[loc + is] - res.data[loc + is];
                } else {
                    res.data[loc + is] = -res.data[loc + is];
                }
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

/* InitUserData: Load problem constants in webdata. */
fn InitUserData(webdata: &mut WebData) {
    webdata.np = NPREY;
    webdata.dx = AX / (MX as f64 - 1.0);
    webdata.dy = AY / (MY as f64 - 1.0);

    let np = webdata.np;
    let dx2 = webdata.dx * webdata.dx;
    let dy2 = webdata.dy * webdata.dy;

    for i in 0..np {
        /*  Fill in the portion of acoef in the four quadrants, row by row. */
        for j in 0..np {
            webdata.acoef[i][np + j] = -GG;
            webdata.acoef[i + np][j] = EE;
            webdata.acoef[i][j] = ZERO;
            webdata.acoef[i + np][np + j] = ZERO;
        }

        /* Reset the diagonal elements of acoef to -AA. */
        webdata.acoef[i][i] = -AA;
        webdata.acoef[i + np][i + np] = -AA;

        /* Set coefficients for b and diffusion terms. */
        webdata.bcoef[i] = BB;
        webdata.bcoef[i + np] = -BB;
        webdata.cox[i] = DPREY / dx2;
        webdata.cox[i + np] = DPRED / dx2;
        webdata.coy[i] = DPREY / dy2;
        webdata.coy[i + np] = DPRED / dy2;
    }
}

/* SetInitialProfiles: Set initial conditions in cc, cp, and id. */
fn SetInitialProfiles(cc: &mut NVector, cp: &mut NVector, id: &mut NVector, webdata: &WebData) {
    let np = webdata.np;

    /* Loop over grid, load cc values and id values. */
    for jy in 0..MY {
        let yy = jy as f64 * webdata.dy;
        let yloc = NSMX * jy;
        for jx in 0..MX {
            let xx = jx as f64 * webdata.dx;
            let mut xyfactor = 16.0 * xx * (ONE - xx) * yy * (ONE - yy);
            xyfactor *= xyfactor;
            let loc = yloc + NUM_SPECIES * jx;

            for is in 0..NUM_SPECIES {
                if is < np {
                    cc.data[loc + is] = 10.0 + (is + 1) as f64 * xyfactor;
                    id.data[loc + is] = ONE;
                } else {
                    cc.data[loc + is] = 1.0e5;
                    id.data[loc + is] = ZERO;
                }
            }
        }
    }

    /* Set c' for the prey by calling the function Fweb. */
    Fweb(cc, cp, webdata);

    /* Set c' for predators to 0. */
    for jy in 0..MY {
        let yloc = NSMX * jy;
        for jx in 0..MX {
            let loc = yloc + NUM_SPECIES * jx;
            for is in np..NUM_SPECIES {
                cp.data[loc + is] = ZERO;
            }
        }
    }
}

/* Print first lines of output (problem description). */
fn PrintHeader(mu: i64, ml: i64, rtol: f64, atol: f64) {
    print!("\nidaFoodWeb_bnd: Predator-prey DAE serial example problem for IDA \n\n");
    print!("Number of species ns: {}", NUM_SPECIES);
    print!("     Mesh dimensions: {} x {}", MX, MY);
    print!("     System size: {}\n", NEQ);
    print!(
        "Tolerance parameters:  rtol = {}   atol = {}\n",
        fmt_g(rtol, 0, 6),
        fmt_g(atol, 0, 6)
    );
    print!("Linear solver: BAND,  Band parameters mu = {}, ml = {}\n", mu, ml);
    print!("CalcIC called to correct initial predator concentrations.\n\n");
    print!("-----------------------------------------------------------\n");
    print!("  t        bottom-left  top-right");
    print!("    | nst  k      h\n");
    print!("-----------------------------------------------------------\n\n");
}

/* PrintOutput: Print output values at output time t = tt. */
fn PrintOutput(mem: &mut IDAMem, c: &NVector, t: f64) {
    let mut kused = 0i32;
    let mut nst = 0i64;
    let mut hused = 0.0f64;

    IDAGetLastOrder(mem, &mut kused);
    IDAGetNumSteps(mem, &mut nst);
    IDAGetLastStep(mem, &mut hused);

    let bl = ij(0, 0);
    let tr = ij(MX - 1, MY - 1);

    println!(
        "{} {} {}   | {:3}  {:1} {}",
        fmt_e(t, 8, 2),
        fmt_e(c.data[bl], 12, 4),
        fmt_e(c.data[tr], 12, 4),
        nst,
        kused,
        fmt_e(hused, 12, 4)
    );
    for i in 1..NUM_SPECIES {
        println!("         {} {}   |", fmt_e(c.data[bl + i], 12, 4), fmt_e(c.data[tr + i], 12, 4));
    }

    println!();
}

/* PrintFinalStats: Print final run data. */
fn PrintFinalStats(mem: &mut IDAMem) {
    let (mut nst, mut nre, mut nreLS, mut nni, mut nnf, mut nje, mut netf, mut ncfn) =
        (0i64, 0i64, 0i64, 0i64, 0i64, 0i64, 0i64, 0i64);

    IDAGetNumSteps(mem, &mut nst);
    IDAGetNumNonlinSolvIters(mem, &mut nni);
    IDAGetNumResEvals(mem, &mut nre);
    IDAGetNumErrTestFails(mem, &mut netf);
    IDAGetNumNonlinSolvConvFails(mem, &mut nnf);
    IDAGetNumStepSolveFails(mem, &mut ncfn);
    IDAGetNumJacEvals(mem, &mut nje);
    IDAGetNumLinResEvals(mem, &mut nreLS);

    print!("-----------------------------------------------------------\n");
    print!("Final run statistics: \n\n");
    println!("Number of steps                    = {}", nst);
    println!("Number of residual evaluations     = {}", nre + nreLS);
    println!("Number of Jacobian evaluations     = {}", nje);
    println!("Number of nonlinear iterations     = {}", nni);
    println!("Number of error test failures      = {}", netf);
    println!("Number of nonlinear conv. failures = {}", nnf);
    println!("Number of step solver failures     = {}", ncfn);
}

/* Fweb: Rate function for the food-web problem. */
fn Fweb(cc: &NVector, crate_: &mut NVector, webdata: &WebData) {
    for jy in 0..MY {
        let yy = webdata.dy * jy as f64;
        let idyu: isize = if jy != MY - 1 { NSMX as isize } else { -(NSMX as isize) };
        let idyl: isize = if jy != 0 { NSMX as isize } else { -(NSMX as isize) };

        for jx in 0..MX {
            let xx = webdata.dx * jx as f64;
            let idxu: isize = if jx != MX - 1 { NUM_SPECIES as isize } else { -(NUM_SPECIES as isize) };
            let idxl: isize = if jx != 0 { NUM_SPECIES as isize } else { -(NUM_SPECIES as isize) };
            let loc = (jy * NSMX + jx * NUM_SPECIES) as isize;

            /* Get interaction vector at this grid point. */
            let mut ratesxy = [0.0f64; NUM_SPECIES];
            WebRates(xx, yy, cc, loc as usize, &mut ratesxy, webdata);

            /* Loop over species, do differencing, load crate segment. */
            for is in 0..NUM_SPECIES {
                let isz = is as isize;
                /* Differencing in y. */
                let dcyli = cc.data[(loc + isz) as usize] - cc.data[(loc - idyl + isz) as usize];
                let dcyui = cc.data[(loc + idyu + isz) as usize] - cc.data[(loc + isz) as usize];

                /* Differencing in x. */
                let dcxli = cc.data[(loc + isz) as usize] - cc.data[(loc - idxl + isz) as usize];
                let dcxui = cc.data[(loc + idxu + isz) as usize] - cc.data[(loc + isz) as usize];

                /* Compute the crate values at (xx,yy). */
                crate_.data[(loc + isz) as usize] = webdata.coy[is] * (dcyui - dcyli)
                    + webdata.cox[is] * (dcxui - dcxli)
                    + ratesxy[is];
            }
        }
    }
}

/* WebRates: Evaluate reaction rates at a given spatial point. */
fn WebRates(xx: f64, yy: f64, cc: &NVector, loc: usize, ratesxy: &mut [f64; NUM_SPECIES], webdata: &WebData) {
    for is in 0..NUM_SPECIES {
        ratesxy[is] = dotprod(cc, loc, &webdata.acoef[is]);
    }

    let fac = ONE + ALPHA * xx * yy + BETA * (FOURPI * xx).sin() * (FOURPI * yy).sin();

    for is in 0..NUM_SPECIES {
        ratesxy[is] = cc.data[loc + is] * (webdata.bcoef[is] * fac + ratesxy[is]);
    }
}

/* dotprod: dot product of cc[loc..loc+ns] with a matrix row. */
fn dotprod(cc: &NVector, loc: usize, arow: &[f64; NUM_SPECIES]) -> f64 {
    let mut temp = ZERO;
    for i in 0..NUM_SPECIES {
        temp += cc.data[loc + i] * arow[i];
    }
    temp
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
    /* Create the SUNDIALS context object for this simulation */
    let sunctx = SUNContext_Create();

    /* Allocate and initialize user data block webdata. */
    let mut webdata = WebData {
        np: 0,
        dx: 0.0,
        dy: 0.0,
        acoef: [[0.0; NUM_SPECIES]; NUM_SPECIES],
        cox: [0.0; NUM_SPECIES],
        coy: [0.0; NUM_SPECIES],
        bcoef: [0.0; NUM_SPECIES],
    };
    InitUserData(&mut webdata);

    /* Allocate N-vectors and initialize cc, cp, and id. */
    let mut cc = N_VNew_Serial(NEQ as i64, &sunctx);
    let mut cp = N_VClone(&cc);
    let mut id = N_VClone(&cc);

    SetInitialProfiles(&mut cc, &mut cp, &mut id, &webdata);

    /* Set remaining inputs to IDAMalloc. */
    let t0 = ZERO;
    let rtol = RTOL;
    let atol = ATOL;

    /* Call IDACreate and IDAMalloc to initialize IDA. */
    let mut mem = IDACreate(&sunctx);

    let mut retval = IDASetUserData(&mut mem, Some(Box::new(webdata)));
    if check_retval(retval, "IDASetUserData") {
        std::process::exit(1);
    }

    retval = IDASetId(&mut mem, Some(&id));
    if check_retval(retval, "IDASetId") {
        std::process::exit(1);
    }

    retval = IDAInit(&mut mem, resweb, t0, &cc, &cp);
    if check_retval(retval, "IDAInit") {
        std::process::exit(1);
    }

    retval = IDASStolerances(&mut mem, rtol, atol);
    if check_retval(retval, "IDASStolerances") {
        std::process::exit(1);
    }

    /* Create banded SUNMatrix for use in linear solves */
    let mu = NSMX as i64;
    let ml = NSMX as i64;
    let a = SUNBandMatrix(NEQ as i64, mu, ml, &sunctx);

    /* Create banded SUNLinearSolver object */
    let ls = SUNLinSol_Band(&cc, &a, &sunctx);

    /* Attach the matrix and linear solver */
    retval = IDASetLinearSolver(&mut mem, ls, Some(a));
    if check_retval(retval, "IDASetLinearSolver") {
        std::process::exit(1);
    }

    /* Call IDACalcIC (with default options) to correct the initial values. */
    let mut tout = 0.001;
    retval = IDACalcIC(&mut mem, IDA_YA_YDP_INIT, tout);
    if check_retval(retval, "IDACalcIC") {
        std::process::exit(1);
    }

    /* Print heading, basic parameters, and initial values. */
    PrintHeader(mu, ml, rtol, atol);
    PrintOutput(&mut mem, &cc, ZERO);

    /* Loop over iout, call IDASolve (normal mode), print selected output. */
    let mut tret = 0.0;
    for iout in 1..=NOUT {
        retval = IDASolve(&mut mem, tout, &mut tret, &mut cc, &mut cp, IDA_NORMAL);
        if check_retval(retval, "IDASolve") {
            std::process::exit(retval);
        }

        PrintOutput(&mut mem, &cc, tret);

        if iout < 3 {
            tout *= TMULT;
        } else {
            tout += TADD;
        }
    }

    /* Print final statistics and free memory. */
    PrintFinalStats(&mut mem);

    /* Free memory (RAII) */
    IDAFree(mem);

    std::process::exit(0);
}
