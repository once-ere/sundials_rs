/* -----------------------------------------------------------------
 * Translated from examples/ida/serial/idaHeat2D_kry.c (IDA 7.7.0)
 * Programmer(s): Allan Taylor, Alan Hindmarsh and Radu Serban @ LLNL
 *
 * Example problem for IDA: 2D heat equation, serial, Krylov (SPGMR)
 * with a diagonal preconditioner. Two cases are run (modified then
 * classical Gram-Schmidt). No IDACalcIC. Output at t = 0.01, .02, ...
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

/* resHeat core: 5-point central differencing; res = u at boundary. */
fn resHeat_compute(mm: usize, coeff: f64, uu: &NVector, up: &NVector, rr: &mut NVector) {
    /* Initialize rr to uu, to take care of boundary equations. */
    N_VScale(ONE, uu, rr);

    /* Loop over interior points; set res = up - (central difference). */
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

/* PsetupHeat: setup for diagonal preconditioner. Keeps only the diagonal
   of J = dF/du + cj*dF/du', stored as inverses in data.pp. Only cj and
   data.coeff are used (ewt/hh from the port's IDALsPrecSetupFn are unused
   here). */
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

    /* Initialize the entire vector to 1., then set the interior points. */
    N_VConst(ONE, &mut data.pp);

    /* Compute the inverse of the preconditioner diagonal elements. */
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

/* PsolveHeat: z = pp .* r (diagonal preconditioner apply). */
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

/* SetInitialProfile: initialize u and up vectors. */
fn SetInitialProfile(data: &HeatData, uu: &mut NVector, up: &mut NVector, res: &mut NVector) {
    let mm = data.mm;
    let mm1 = mm - 1;

    /* Initialize uu on all grid points. */
    for j in 0..mm {
        let yfact = data.dx * j as f64;
        let offset = mm * j;
        for i in 0..mm {
            let xfact = data.dx * i as f64;
            let loc = offset + i;
            uu.data[loc] = 16.0 * xfact * (ONE - xfact) * yfact * (ONE - yfact);
        }
    }

    /* Initialize up vector to 0. */
    N_VConst(ZERO, up);

    /* resHeat sets res to negative of ODE RHS values at interior points. */
    resHeat_compute(data.mm, data.coeff, uu, up, res);

    /* Copy -res into up to get correct interior initial up values. */
    N_VScale(-ONE, res, up);

    /* Set up at boundary points to zero. */
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

/* Print first lines of output (problem description). */
fn PrintHeader(rtol: f64, atol: f64) {
    print!("\nidaHeat2D_kry: Heat equation, serial example problem for IDA \n");
    print!("         Discretized heat equation on 2D unit square. \n");
    print!("         Zero boundary conditions,");
    print!(" polynomial initial conditions.\n");
    print!("         Mesh dimensions: {} x {}", MGRID, MGRID);
    print!("        Total system size: {}\n\n", NEQ);
    print!(
        "Tolerance parameters:  rtol = {}   atol = {}\n",
        fmt_g(rtol, 0, 6),
        fmt_g(atol, 0, 6)
    );
    print!("Constraints set to force all solution components >= 0. \n");
    print!("Linear solver: SPGMR, preconditioner using diagonal elements. \n");
}

/* PrintOutput: print max norm of solution and current solver statistics. */
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

    /* Create the linear solver SUNLinSol_SPGMR with left preconditioning
       and the default Krylov dimension */
    let mut ls = SUNLinSol_SPGMR(&uu, SUN_PREC_LEFT, 0, &sunctx);

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
    retval = IDASetPreconditioner(&mut mem, Some(PsetupHeat), Some(PsolveHeat));
    if check_retval(retval, "IDASetPreconditioner") {
        std::process::exit(1);
    }

    /* Print output heading. */
    PrintHeader(rtol, atol);

    /* CASE I */
    print!("\n\nCase 1: gsytpe = SUN_MODIFIED_GS\n");
    print!("\n   Output Summary (umax = max-norm of solution) \n\n");
    print!("  time     umax       k  nst  nni  nje   nre   nreLS    h      npe nps\n");
    print!("----------------------------------------------------------------------\n");

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

    let mut netf = 0i64;
    let mut ncfn = 0i64;
    let mut ncfl = 0i64;
    IDAGetNumErrTestFails(&mut mem, &mut netf);
    IDAGetNumNonlinSolvConvFails(&mut mem, &mut ncfn);
    IDAGetNumLinConvFails(&mut mem, &mut ncfl);
    println!("\nError test failures            = {}", netf);
    println!("Nonlinear convergence failures = {}", ncfn);
    println!("Linear convergence failures    = {}", ncfl);

    /* CASE II */
    SetInitialProfile_reinit(&mut mem, &mut uu, &mut up, &mut res);

    retval = IDAReInit(&mut mem, t0, &uu, &up);
    if check_retval(retval, "IDAReInit") {
        std::process::exit(1);
    }

    retval = match &mut mem.ida_lmem {
        LsModule::Ls(idals) => match &mut idals.LS {
            LinearSolver::Spgmr(s) => s.set_gs_type(SUN_CLASSICAL_GS),
            _ => -1,
        },
        _ => -1,
    };
    if check_retval(retval, "SUNLinSol_SPGMRSetGSType") {
        std::process::exit(1);
    }

    print!("\n\nCase 2: gstype = SUN_CLASSICAL_GS\n");
    print!("\n   Output Summary (umax = max-norm of solution) \n\n");
    print!("  time     umax       k  nst  nni  nje   nre   nreLS    h      npe nps\n");
    print!("----------------------------------------------------------------------\n");

    tout = t1;
    for _iout in 1..=NOUT {
        retval = IDASolve(&mut mem, tout, &mut tret, &mut uu, &mut up, IDA_NORMAL);
        if check_retval(retval, "IDASolve") {
            std::process::exit(1);
        }
        PrintOutput(&mut mem, tret, &uu);
        tout *= TWO;
    }

    IDAGetNumErrTestFails(&mut mem, &mut netf);
    IDAGetNumNonlinSolvConvFails(&mut mem, &mut ncfn);
    IDAGetNumLinConvFails(&mut mem, &mut ncfl);
    println!("\nError test failures            = {}", netf);
    println!("Nonlinear convergence failures = {}", ncfn);
    println!("Linear convergence failures    = {}", ncfl);

    /* Free memory (RAII) */
    IDAFree(mem);
}

/* Re-run SetInitialProfile using the HeatData stored in ida_user_data. */
fn SetInitialProfile_reinit(mem: &mut IDAMem, uu: &mut NVector, up: &mut NVector, res: &mut NVector) {
    let data = mem.ida_user_data.as_ref().unwrap().downcast_ref::<HeatData>().unwrap();
    /* borrow only the fields we need (copy scalars) to avoid aliasing mem */
    let (mm, dx, coeff) = (data.mm, data.dx, data.coeff);
    let tmp = HeatData { mm, dx, coeff, pp: N_VClone(uu) };
    SetInitialProfile(&tmp, uu, up, res);
}
