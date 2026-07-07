/* -----------------------------------------------------------------
 * Translated from examples/ida/serial/idaFoodWeb_kry.c (IDA 7.7.0)
 * Programmer(s): Ting Yan @ UMBC
 *
 * Example program for IDA: Food web problem, GMRES, user-supplied
 * preconditioner (serial).
 *
 * The DAE system (predator-prey food web, central-differenced on an
 * MX by MY mesh, np = 1, ns = 2, homogeneous Neumann BC) is solved by
 * IDA using the SPGMR linear solver with a user-supplied block-diagonal
 * preconditioner (one NUM_SPECIES x NUM_SPECIES dense block per grid
 * point) and IDACalcIC for initial condition calculation.
 *
 * Output is printed at t = 0, .001, .01, .1, .4, .7, 1.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use ida_rs::sundials_dense::{SUNDlsMat_denseGETRF, SUNDlsMat_denseGETRS};
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

/* index into the per-grid-point preconditioner block storage ([jx][jy]) */
fn pij(jx: usize, jy: usize) -> usize {
    jx * MY + jy
}

/* Type: WebData.  Contains problem constants, preconditioner blocks, and
   the reaction-rate scratch vector (rates from the last Fweb call, reused
   by Precond). */
struct WebData {
    np: usize,
    dx: f64,
    dy: f64,
    acoef: [[f64; NUM_SPECIES]; NUM_SPECIES],
    cox: [f64; NUM_SPECIES],
    coy: [f64; NUM_SPECIES],
    bcoef: [f64; NUM_SPECIES],
    rates: NVector,
    pp: Vec<DenseMatrix>, /* MX*MY dense NUM_SPECIES x NUM_SPECIES blocks */
    pivot: Vec<[sunindextype; NUM_SPECIES]>,
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY IDA
 *--------------------------------------------------------------------
 */

/* resweb: System residual function for predator-prey system. */
fn resweb(_tt: f64, cc: &NVector, cp: &NVector, res: &mut NVector, user_data: &mut UserData) -> i32 {
    let webdata = user_data.as_mut().unwrap().downcast_mut::<WebData>().unwrap();
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

/* Precond: setup the block-diagonal preconditioner (one dense block per
   grid point) via difference-quotient Jacobian of the reaction terms.
   In C the current error weights and step size are obtained inside the
   callback via IDAGetErrWeights / IDAGetCurrentStep; here the port hands
   them in as `ewt` and `hh` (see IDALsPrecSetupFn note). */
fn Precond(
    _tt: f64,
    cc: &NVector,
    cp: &NVector,
    _rr: &NVector,
    cj: f64,
    ewt: &NVector,
    hh: f64,
    user_data: &mut UserData,
) -> i32 {
    let webdata = user_data.as_mut().unwrap().downcast_mut::<WebData>().unwrap();
    let del_x = webdata.dx;
    let del_y = webdata.dy;

    let uround = SUN_UNIT_ROUNDOFF;
    let sqru = uround.sqrt();

    for jy in 0..MY {
        let yy = jy as f64 * del_y;

        for jx in 0..MX {
            let xx = jx as f64 * del_x;
            let loc = ij(jx, jy);
            let pidx = pij(jx, jy);

            /* unperturbed concentrations at this grid point */
            let cxy0 = [cc.data[loc], cc.data[loc + 1]];

            for js in 0..NUM_SPECIES {
                let inc = sqru
                    * SUNRabs(cxy0[js])
                        .max((hh * SUNRabs(cp.data[loc + js])).max(ONE / ewt.data[loc + js]));
                let fac = -ONE / inc;

                /* perturb species js and evaluate reaction rates */
                let mut cxy = cxy0;
                cxy[js] += inc;
                let mut perturb_rates = [0.0f64; NUM_SPECIES];
                WebRatesLocal(xx, yy, &cxy, &mut perturb_rates, webdata);

                for is in 0..NUM_SPECIES {
                    webdata.pp[pidx].set(
                        is as i64,
                        js as i64,
                        (perturb_rates[is] - webdata.rates.data[loc + is]) * fac,
                    );
                }

                if js < 1 {
                    let d = webdata.pp[pidx].get(js as i64, js as i64) + cj;
                    webdata.pp[pidx].set(js as i64, js as i64, d);
                }
            }

            let ret = SUNDlsMat_denseGETRF(&mut webdata.pp[pidx], &mut webdata.pivot[pidx]);
            if ret != 0 {
                return 1;
            }
        }
    }

    0
}

/* PSolve: apply the block-diagonal preconditioner. */
#[allow(clippy::too_many_arguments)]
fn PSolve(
    _tt: f64,
    _cc: &NVector,
    _cp: &NVector,
    _rr: &NVector,
    rvec: &NVector,
    zvec: &mut NVector,
    _cj: f64,
    _delta: f64,
    user_data: &mut UserData,
) -> i32 {
    let webdata = user_data.as_mut().unwrap().downcast_mut::<WebData>().unwrap();

    N_VScale(ONE, rvec, zvec);

    for jx in 0..MX {
        for jy in 0..MY {
            let loc = ij(jx, jy);
            let pidx = pij(jx, jy);
            SUNDlsMat_denseGETRS(
                &webdata.pp[pidx],
                &webdata.pivot[pidx],
                &mut zvec.data[loc..loc + NUM_SPECIES],
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

/* InitUserData: Load problem constants in webdata. */
fn InitUserData(webdata: &mut WebData) {
    webdata.np = NPREY;
    webdata.dx = AX / (MX as f64 - 1.0);
    webdata.dy = AY / (MY as f64 - 1.0);

    let np = webdata.np;
    let dx2 = webdata.dx * webdata.dx;
    let dy2 = webdata.dy * webdata.dy;

    for i in 0..np {
        for j in 0..np {
            webdata.acoef[i][np + j] = -GG;
            webdata.acoef[i + np][j] = EE;
            webdata.acoef[i][j] = ZERO;
            webdata.acoef[i + np][np + j] = ZERO;
        }

        webdata.acoef[i][i] = -AA;
        webdata.acoef[i + np][i + np] = -AA;

        webdata.bcoef[i] = BB;
        webdata.bcoef[i + np] = -BB;
        webdata.cox[i] = DPREY / dx2;
        webdata.cox[i + np] = DPRED / dx2;
        webdata.coy[i] = DPREY / dy2;
        webdata.coy[i + np] = DPRED / dy2;
    }
}

/* SetInitialProfiles: Set initial conditions in cc, cp, and id. */
fn SetInitialProfiles(cc: &mut NVector, cp: &mut NVector, id: &mut NVector, webdata: &mut WebData) {
    let np = webdata.np;

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
fn PrintHeader(maxl: i32, rtol: f64, atol: f64) {
    print!("\nidaFoodWeb_kry: Predator-prey DAE serial example problem for IDA \n\n");
    print!("Number of species ns: {}", NUM_SPECIES);
    print!("     Mesh dimensions: {} x {}", MX, MY);
    print!("     System size: {}\n", NEQ);
    print!(
        "Tolerance parameters:  rtol = {}   atol = {}\n",
        fmt_g(rtol, 0, 6),
        fmt_g(atol, 0, 6)
    );
    print!("Linear solver: SPGMR,  SPGMR parameters maxl = {}\n", maxl);
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
    let (mut nst, mut nre, mut sli, mut netf, mut nps, mut npevals, mut nrevalsLS) =
        (0i64, 0i64, 0i64, 0i64, 0i64, 0i64, 0i64);

    IDAGetNumSteps(mem, &mut nst);
    IDAGetNumLinIters(mem, &mut sli);
    IDAGetNumResEvals(mem, &mut nre);
    IDAGetNumErrTestFails(mem, &mut netf);
    IDAGetNumPrecSolves(mem, &mut nps);
    IDAGetNumPrecEvals(mem, &mut npevals);
    IDAGetNumLinResEvals(mem, &mut nrevalsLS);

    print!("-----------------------------------------------------------\n");
    print!("Final run statistics: \n\n");
    println!("Number of steps                       = {}", nst);
    println!("Number of residual evaluations        = {}", nre);
    println!("Number of Preconditioner evaluations  = {}", npevals);
    println!("Number of linear iterations           = {}", sli);
    println!("Number of error test failures         = {}", netf);
    println!("Number of precond solve fun called    = {}", nps);
}

/* Fweb: Rate function for the food-web problem. Writes the reaction rates
   into webdata.rates (reused by Precond) and the RHS into crate_. */
fn Fweb(cc: &NVector, crate_: &mut NVector, webdata: &mut WebData) {
    for jy in 0..MY {
        let yy = webdata.dy * jy as f64;
        let idyu: isize = if jy != MY - 1 { NSMX as isize } else { -(NSMX as isize) };
        let idyl: isize = if jy != 0 { NSMX as isize } else { -(NSMX as isize) };

        for jx in 0..MX {
            let xx = webdata.dx * jx as f64;
            let idxu: isize = if jx != MX - 1 { NUM_SPECIES as isize } else { -(NUM_SPECIES as isize) };
            let idxl: isize = if jx != 0 { NUM_SPECIES as isize } else { -(NUM_SPECIES as isize) };
            let loc = (jy * NSMX + jx * NUM_SPECIES) as isize;
            let locu = loc as usize;

            /* Get interaction vector at this grid point; store into rates. */
            let cxy = [cc.data[locu], cc.data[locu + 1]];
            let mut ratesxy = [0.0f64; NUM_SPECIES];
            WebRatesLocal(xx, yy, &cxy, &mut ratesxy, webdata);
            for is in 0..NUM_SPECIES {
                webdata.rates.data[locu + is] = ratesxy[is];
            }

            /* Loop over species, do differencing, load crate segment. */
            for is in 0..NUM_SPECIES {
                let isz = is as isize;
                let dcyli = cc.data[(loc + isz) as usize] - cc.data[(loc - idyl + isz) as usize];
                let dcyui = cc.data[(loc + idyu + isz) as usize] - cc.data[(loc + isz) as usize];

                let dcxli = cc.data[(loc + isz) as usize] - cc.data[(loc - idxl + isz) as usize];
                let dcxui = cc.data[(loc + idxu + isz) as usize] - cc.data[(loc + isz) as usize];

                crate_.data[(loc + isz) as usize] = webdata.coy[is] * (dcyui - dcyli)
                    + webdata.cox[is] * (dcxui - dcxli)
                    + ratesxy[is];
            }
        }
    }
}

/* WebRatesLocal: Evaluate reaction rates for one grid point given the
   NUM_SPECIES local concentrations. */
fn WebRatesLocal(xx: f64, yy: f64, cxy: &[f64; NUM_SPECIES], ratesxy: &mut [f64; NUM_SPECIES], webdata: &WebData) {
    for is in 0..NUM_SPECIES {
        ratesxy[is] = dotprod(cxy, &webdata.acoef[is]);
    }

    let fac = ONE + ALPHA * xx * yy + BETA * (FOURPI * xx).sin() * (FOURPI * yy).sin();

    for is in 0..NUM_SPECIES {
        ratesxy[is] = cxy[is] * (webdata.bcoef[is] * fac + ratesxy[is]);
    }
}

/* dotprod: dot product of two NUM_SPECIES arrays. */
fn dotprod(x1: &[f64; NUM_SPECIES], x2: &[f64; NUM_SPECIES]) -> f64 {
    let mut temp = ZERO;
    for i in 0..NUM_SPECIES {
        temp += x1[i] * x2[i];
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
    let mut pp = Vec::with_capacity(MX * MY);
    for _ in 0..(MX * MY) {
        pp.push(DenseMatrix::new(NUM_SPECIES as i64, NUM_SPECIES as i64));
    }
    let mut webdata = WebData {
        np: 0,
        dx: 0.0,
        dy: 0.0,
        acoef: [[0.0; NUM_SPECIES]; NUM_SPECIES],
        cox: [0.0; NUM_SPECIES],
        coy: [0.0; NUM_SPECIES],
        bcoef: [0.0; NUM_SPECIES],
        rates: N_VNew_Serial(NEQ as i64, &sunctx),
        pp,
        pivot: vec![[0; NUM_SPECIES]; MX * MY],
    };
    InitUserData(&mut webdata);

    /* Allocate N-vectors and initialize cc, cp, and id. */
    let mut cc = N_VNew_Serial(NEQ as i64, &sunctx);
    let mut cp = N_VClone(&cc);
    let mut id = N_VClone(&cc);

    SetInitialProfiles(&mut cc, &mut cp, &mut id, &mut webdata);

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

    /* Create the linear solver SUNLinSol_SPGMR with left preconditioning
       and maximum Krylov dimension maxl */
    let maxl = 16;
    let mut ls = SUNLinSol_SPGMR(&cc, SUN_PREC_LEFT, maxl, &sunctx);

    /* IDA recommends allowing up to 5 restarts (default is 0) */
    retval = match &mut ls {
        LinearSolver::Spgmr(s) => s.set_max_restarts(5),
        _ => -1,
    };
    if check_retval(retval, "SUNLinSol_SPGMRSetMaxRestarts") {
        std::process::exit(1);
    }

    /* Attach the linear solver */
    retval = IDASetLinearSolver(&mut mem, ls, None);
    if check_retval(retval, "IDASetLinearSolver") {
        std::process::exit(1);
    }

    /* Set the preconditioner solve and setup functions */
    retval = IDASetPreconditioner(&mut mem, Some(Precond), Some(PSolve));
    if check_retval(retval, "IDASetPreconditioner") {
        std::process::exit(1);
    }

    /* Call IDACalcIC (with default options) to correct the initial values. */
    let mut tout = 0.001;
    retval = IDACalcIC(&mut mem, IDA_YA_YDP_INIT, tout);
    if check_retval(retval, "IDACalcIC") {
        std::process::exit(1);
    }

    /* Print heading, basic parameters, and initial values. */
    PrintHeader(maxl, rtol, atol);
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
