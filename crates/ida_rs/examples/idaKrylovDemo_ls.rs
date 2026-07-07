/* -----------------------------------------------------------------
 * Translated from examples/ida/serial/idaKrylovDemo_ls.c (IDA 7.7.0)
 * Programmer(s): Allan Taylor, Alan Hindmarsh and Radu Serban @ LLNL
 *
 * Loops through the Krylov linear solvers SPGMR, SPBCGS and SPTFQMR
 * on the serial 2D heat-equation DAE with a diagonal preconditioner.
 * No IDACalcIC. Output at t = 0.01, .02, ..., 10.24.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use ida_rs::sundials_utils::{fmt_e, fmt_f, fmt_g};
use ida_rs::*;

/* Problem Constants */
const NOUT: i32 = 11;
const MGRID: usize = 10;
const NEQ: usize = MGRID * MGRID;
const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;
const FOUR: f64 = 4.0;

/* Linear Solver Loop Constants */
const USE_SPGMR: i32 = 0;
const USE_SPBCG: i32 = 1;
const USE_SPTFQMR: i32 = 2;

/* User data type */
struct HeatData {
    mm: usize,
    dx: f64,
    coeff: f64,
    pp: NVector, /* inverse-diagonal preconditioner vector */
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY IDA
 *--------------------------------------------------------------------
 */

fn resHeat_compute(mm: usize, coeff: f64, uu: &NVector, up: &NVector, rr: &mut NVector) {
    N_VScale(ONE, uu, rr);
    for j in 1..MGRID - 1 {
        let offset = mm * j;
        for i in 1..mm - 1 {
            let loc = offset + i;
            let dif1 = uu.data[loc - 1] + uu.data[loc + 1] - TWO * uu.data[loc];
            let dif2 = uu.data[loc - mm] + uu.data[loc + mm] - TWO * uu.data[loc];
            rr.data[loc] = up.data[loc] - coeff * (dif1 + dif2);
        }
    }
}

fn resHeat(_tt: f64, uu: &NVector, up: &NVector, rr: &mut NVector, user_data: &mut UserData) -> i32 {
    let data = user_data.as_ref().unwrap().downcast_ref::<HeatData>().unwrap();
    resHeat_compute(data.mm, data.coeff, uu, up, rr);
    0
}

/* PsetupHeat: diagonal preconditioner. Only cj and coeff are used. */
fn PsetupHeat(
    _tt: f64,
    _uu: &NVector,
    _up: &NVector,
    _rr: &NVector,
    c_j: f64,
    _ewt: &NVector,
    _hh: f64,
    prec_data: &mut UserData,
) -> i32 {
    let data = prec_data.as_mut().unwrap().downcast_mut::<HeatData>().unwrap();
    let mm = data.mm;

    N_VConst(ONE, &mut data.pp);
    let pelinv = ONE / (c_j + FOUR * data.coeff);
    for j in 1..mm - 1 {
        let offset = mm * j;
        for i in 1..mm - 1 {
            let loc = offset + i;
            data.pp.data[loc] = pelinv;
        }
    }
    0
}

/* PsolveHeat: z = pp .* r. */
#[allow(clippy::too_many_arguments)]
fn PsolveHeat(
    _tt: f64,
    _uu: &NVector,
    _up: &NVector,
    _rr: &NVector,
    rvec: &NVector,
    zvec: &mut NVector,
    _c_j: f64,
    _delta: f64,
    prec_data: &mut UserData,
) -> i32 {
    let data = prec_data.as_ref().unwrap().downcast_ref::<HeatData>().unwrap();
    N_VProd(&data.pp, rvec, zvec);
    0
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

fn SetInitialProfile(data: &HeatData, uu: &mut NVector, up: &mut NVector, res: &mut NVector) {
    let mm = data.mm;
    let mm1 = mm - 1;

    for j in 0..mm {
        let yfact = data.dx * j as f64;
        let offset = mm * j;
        for i in 0..mm {
            let xfact = data.dx * i as f64;
            let loc = offset + i;
            uu.data[loc] = 16.0 * xfact * (ONE - xfact) * yfact * (ONE - yfact);
        }
    }

    N_VConst(ZERO, up);
    resHeat_compute(data.mm, data.coeff, uu, up, res);
    N_VScale(-ONE, res, up);

    for j in 0..mm {
        let offset = mm * j;
        for i in 0..mm {
            let loc = offset + i;
            if j == 0 || j == mm1 || i == 0 || i == mm1 {
                up.data[loc] = ZERO;
            }
        }
    }
}

fn PrintHeader(rtol: f64, atol: f64, linsolver: i32) {
    print!("\nidaKrylovDemo_ls: Heat equation, serial example problem for IDA\n");
    print!("               Discretized heat equation on 2D unit square.\n");
    print!("               Zero boundary conditions,");
    print!(" polynomial initial conditions.\n");
    print!("         Mesh dimensions: {} x {}", MGRID, MGRID);
    print!("       Total system size: {}\n\n", NEQ);
    print!(
        "Tolerance parameters:  rtol = {}   atol = {}\n",
        fmt_g(rtol, 0, 6),
        fmt_g(atol, 0, 6)
    );
    print!("Constraints set to force all solution components >= 0. \n");

    match linsolver {
        USE_SPGMR => print!("Linear solver: SPGMR, preconditioner using diagonal elements. \n"),
        USE_SPBCG => print!("Linear solver: SPBCG, preconditioner using diagonal elements. \n"),
        USE_SPTFQMR => print!("Linear solver: SPTFQMR, preconditioner using diagonal elements. \n"),
        _ => {}
    }
}

fn PrintOutput(mem: &mut IDAMem, t: f64, uu: &NVector) {
    let umax = N_VMaxNorm(uu);

    let mut kused = 0i32;
    let (mut nst, mut nni, mut nje, mut nre, mut nreLS, mut nli, mut npe, mut nps) =
        (0i64, 0i64, 0i64, 0i64, 0i64, 0i64, 0i64, 0i64);
    let mut hused = 0.0f64;

    IDAGetLastOrder(mem, &mut kused);
    IDAGetNumSteps(mem, &mut nst);
    IDAGetNumNonlinSolvIters(mem, &mut nni);
    IDAGetNumResEvals(mem, &mut nre);
    IDAGetLastStep(mem, &mut hused);
    IDAGetNumJtimesEvals(mem, &mut nje);
    IDAGetNumLinIters(mem, &mut nli);
    IDAGetNumLinResEvals(mem, &mut nreLS);
    IDAGetNumPrecEvals(mem, &mut npe);
    IDAGetNumPrecSolves(mem, &mut nps);
    let _ = nli;

    println!(
        " {} {}  {}  {:3}  {:3}  {:3}  {:4}  {:4}  {}  {:3} {:3}",
        fmt_f(t, 5, 2),
        fmt_e(umax, 13, 5),
        kused,
        nst,
        nni,
        nje,
        nre,
        nreLS,
        fmt_e(hused, 9, 2),
        npe,
        nps
    );
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
    /* Retrieve the command-line options */
    let args: Vec<String> = std::env::args().collect();
    let nrmfactor: i32 = if args.len() > 1 { args[1].parse().unwrap_or(0) } else { 0 };

    /* Create the SUNDIALS context object for this simulation */
    let sunctx = SUNContext_Create();

    /* Allocate N-vectors and the user data structure. */
    let mut uu = N_VNew_Serial(NEQ as i64, &sunctx);
    let mut up = N_VClone(&uu);
    let mut res = N_VClone(&uu);
    let mut constraints = N_VClone(&uu);

    let data = HeatData {
        mm: MGRID,
        dx: ONE / (MGRID as f64 - ONE),
        coeff: ONE / ((ONE / (MGRID as f64 - ONE)) * (ONE / (MGRID as f64 - ONE))),
        pp: N_VClone(&uu),
    };

    /* Initialize uu, up. */
    SetInitialProfile(&data, &mut uu, &mut up, &mut res);

    /* Set constraints to all 1's for nonnegative solution values. */
    N_VConst(ONE, &mut constraints);

    /* Assign various parameters. */
    let t0 = ZERO;
    let t1 = 0.01;
    let rtol = ZERO;
    let atol = 1.0e-3;

    /* Call IDACreate and IDAMalloc to initialize solution */
    let mut mem = IDACreate(&sunctx);

    let mut retval = IDASetUserData(&mut mem, Some(Box::new(data)));
    if check_retval(retval, "IDASetUserData") {
        std::process::exit(1);
    }

    retval = IDASetConstraints(&mut mem, Some(&constraints));
    if check_retval(retval, "IDASetConstraints") {
        std::process::exit(1);
    }
    drop(constraints);

    retval = IDAInit(&mut mem, resHeat, t0, &uu, &up);
    if check_retval(retval, "IDAInit") {
        std::process::exit(1);
    }

    retval = IDASStolerances(&mut mem, rtol, atol);
    if check_retval(retval, "IDASStolerances") {
        std::process::exit(1);
    }

    /* START: Loop through SPGMR, SPBCG and SPTFQMR linear solver modules */
    for linsolver in 0..3 {
        if linsolver != 0 {
            /* Re-initialize uu, up (fetch HeatData scalars from user_data). */
            let (mm, dx, coeff) = {
                let d = mem.ida_user_data.as_ref().unwrap().downcast_ref::<HeatData>().unwrap();
                (d.mm, d.dx, d.coeff)
            };
            let tmp = HeatData { mm, dx, coeff, pp: N_VClone(&uu) };
            SetInitialProfile(&tmp, &mut uu, &mut up, &mut res);

            retval = IDAReInit(&mut mem, t0, &uu, &up);
            if check_retval(retval, "IDAReInit") {
                std::process::exit(1);
            }
        }

        /* Attach a new linear solver module */
        let ls = match linsolver {
            USE_SPGMR => {
                print!(" -------");
                print!(" \n| SPGMR |\n");
                print!(" -------\n");
                SUNLinSol_SPGMR(&uu, SUN_PREC_LEFT, 0, &sunctx)
            }
            USE_SPBCG => {
                print!(" -------");
                print!(" \n| SPBCGS |\n");
                print!(" -------\n");
                SUNLinSol_SPBCGS(&uu, SUN_PREC_LEFT, 0, &sunctx)
            }
            _ => {
                print!(" ---------");
                print!(" \n| SPTFQMR |\n");
                print!(" ---------\n");
                SUNLinSol_SPTFQMR(&uu, SUN_PREC_LEFT, 0, &sunctx)
            }
        };

        retval = IDASetLinearSolver(&mut mem, ls, None);
        if check_retval(retval, "IDASetLinearSolver") {
            std::process::exit(1);
        }

        /* Specify preconditioner */
        retval = IDASetPreconditioner(&mut mem, Some(PsetupHeat), Some(PsolveHeat));
        if check_retval(retval, "IDASetPreconditioner") {
            std::process::exit(1);
        }

        /* Set the linear solver tolerance conversion factor */
        let nrmfac = match nrmfactor {
            1 => (NEQ as f64).sqrt(),
            2 => -ONE,
            _ => ZERO,
        };
        retval = IDASetLSNormFactor(&mut mem, nrmfac);
        if check_retval(retval, "IDASetLSNormFactor") {
            std::process::exit(1);
        }

        /* Print output heading. */
        PrintHeader(rtol, atol, linsolver);

        print!("\n   Output Summary (umax = max-norm of solution) \n\n");
        print!("  time     umax       k  nst  nni  nje   nre   nreLS    h      npe nps\n");
        print!("----------------------------------------------------------------------\n");

        /* Loop over output times, call IDASolve, and print results. */
        let mut tret = 0.0;
        let mut tout = t1;
        for _iout in 1..=NOUT {
            retval = IDASolve(&mut mem, tout, &mut tret, &mut uu, &mut up, IDA_NORMAL);
            if check_retval(retval, "IDASolve") {
                std::process::exit(1);
            }
            PrintOutput(&mut mem, tret, &uu);
            tout *= TWO;
        }

        /* Print remaining counters. */
        let mut netf = 0i64;
        let mut ncfn = 0i64;
        let mut ncfl = 0i64;
        IDAGetNumErrTestFails(&mut mem, &mut netf);
        IDAGetNumNonlinSolvConvFails(&mut mem, &mut ncfn);
        IDAGetNumLinConvFails(&mut mem, &mut ncfl);
        println!("\nError test failures            = {}", netf);
        println!("Nonlinear convergence failures = {}", ncfn);
        println!("Linear convergence failures    = {}", ncfl);

        if linsolver < 2 {
            print!("\n======================================================================\n\n");
        }
    }

    /* Free memory (RAII) */
    IDAFree(mem);
}
