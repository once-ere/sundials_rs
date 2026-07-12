/* ----------------------------------------------------------------
 * Translation of examples/arkode/C_serial/ark_kpr_mri.c
 * (SUNDIALS 7.7.0).
 *
 * Multirate nonlinear Kvaerno-Prothero-Robinson ODE test problem:
 *
 *    [u]' = [ G  e ] [(-1+u^2-r)/(2u)] + [      r'(t)/(2u)        ]
 *    [v]    [ e -1 ] [(-2+v^2-s)/(2v)]   [ s'(t)/(2*sqrt(2+s(t))) ]
 *         = [ fs(t,u,v) ]
 *           [ ff(t,u,v) ]
 *
 * where r(t) = 0.5*cos(t),  s(t) = cos(w*t),  0 < t < 5.
 *
 * This problem has analytical solution given by
 *    u(t) = sqrt(1+r(t)),  v(t) = sqrt(2+s(t)).
 *
 * We use the parameters:
 *   e = 0.5 (fast/slow coupling strength) [default]
 *   G = -1e2 (stiffness at slow time scale) [default]
 *   w = 100  (time-scale separation factor) [default]
 *   hs = 0.01 (slow step size) [default]
 *
 * We select the MRI method to use based on additional inputs:
 *
 *   slow_type:
 *      0 - none (full problem at fast scale)
 *      1 - ARKODE_MIS_KW3
 *      2 - ARKODE_MRI_GARK_ERK45a
 *      3 - ARKODE_MERK21
 *      4 - ARKODE_MERK32
 *      5 - ARKODE_MERK43
 *      6 - ARKODE_MERK54
 *      7 - ARKODE_MRI_GARK_IRK21a
 *      8 - ARKODE_MRI_GARK_ESDIRK34a
 *      9 - ARKODE_IMEX_MRI_GARK3b
 *     10 - ARKODE_IMEX_MRI_GARK4
 *     11 - ARKODE_IMEX_MRI_SR21
 *     12 - ARKODE_IMEX_MRI_SR32
 *     13 - ARKODE_IMEX_MRI_SR43
 *
 *   fast_type:
 *      0 - none (full problem at slow scale)
 *      1 - esdirk-3-3 (manually entered non-embedded table)
 *      2 - ARKODE_HEUN_EULER_2_1_2
 *      3 - erk-3-3 (manually entered non-embedded table)
 *      4 - erk-4-4 (manually entered non-embeded table)
 *      5 - ARKODE_DORMAND_PRINCE_7_4_5
 *
 * The program should be run with arguments in the following order:
 *   $ ark_kpr_mri slow_type fast_type h G w e deduce_rhs
 * (trailing arguments may be omitted from end-to-beginning).
 *
 * This program solves the problem with the MRI stepper. Outputs are
 * printed at equal intervals of 0.1 and run statistics are printed
 * at the end.
 *
 * Note on ownership: the stepper owns the wrapped inner integrator
 * (and the outer integrator owns the stepper), so the final fast
 * statistics are read by borrowing the inner integrator back out of
 * the outer step memory.  The shared C rpar pointer becomes one
 * identical [G,w,e] copy per integrator (immutable during the run).
 * ----------------------------------------------------------------*/
#![allow(non_snake_case)]

use std::io::Write;

use arkode_rs::arkode::{ARKodeCreateMRIStepInnerStepper, ARKodeEvolve, ARKodeFree};
use arkode_rs::arkode::ARKodeSStolerances;
use arkode_rs::arkode_arkstep::ARKStepCreate;
use arkode_rs::arkode_arkstep_io::ARKStepSetTables;
use arkode_rs::arkode_butcher::ARKodeButcherTable_Alloc;
use arkode_rs::arkode_butcher_erk::{
    ARKodeButcherTable_LoadERK, ARKODE_DORMAND_PRINCE_7_4_5, ARKODE_HEUN_EULER_2_1_2,
};
use arkode_rs::arkode_cli::ARKodeSetOptions;
use arkode_rs::arkode_io::{
    ARKodeGetNonlinSolvStats, ARKodeGetNumRhsEvals, ARKodeGetNumSteps,
    ARKodeSetDeduceImplicitRhs, ARKodeSetFixedStep, ARKodeSetUserData,
};
use arkode_rs::arkode_ls::{ARKodeGetNumJacEvals, ARKodeSetJacFn, ARKodeSetLinearSolver};
use arkode_rs::arkode_mri_tables::{
    MRIStepCoupling_LoadTable, MRIStepCoupling_MIStoMRI, ARKODE_IMEX_MRI_GARK3b,
    ARKODE_IMEX_MRI_GARK4, ARKODE_IMEX_MRI_SR21, ARKODE_IMEX_MRI_SR32, ARKODE_IMEX_MRI_SR43,
    ARKODE_MERK21, ARKODE_MERK32, ARKODE_MERK43, ARKODE_MERK54, ARKODE_MIS_KW3,
    ARKODE_MRI_GARK_ERK45a, ARKODE_MRI_GARK_ESDIRK34a, ARKODE_MRI_GARK_IRK21a,
};
use arkode_rs::arkode_mristep::MRIStepCreate;
use arkode_rs::arkode_mristep_impl::ARKodeMRIStepMem;
use arkode_rs::arkode_mristep_io::MRIStepSetCoupling;
use arkode_rs::sundials_utils::{fmt_e, fmt_f, fmt_g};
use arkode_rs::*;

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;
const TWO: f64 = 2.0;

/* ------------------------------
 * Functions called by the solver
 * ------------------------------*/

fn get_rpar(user_data: &UserData) -> &[f64; 3] {
    user_data.as_ref().unwrap().downcast_ref::<[f64; 3]>().unwrap()
}

/* ff routine to compute the fast portion of the ODE RHS. */
fn ff(t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let rpar = *get_rpar(user_data);
    let e = rpar[2];
    let u = y.data[0];
    let v = y.data[1];

    /* fill in the RHS function:
       [0  0]*[(-1+u^2-r(t))/(2*u)] + [         0          ]
       [e -1] [(-2+v^2-s(t))/(2*v)]   [sdot(t)/(2*vtrue(t))] */
    let tmp1 = (-ONE + u * u - r(t, &rpar)) / (TWO * u);
    let tmp2 = (-TWO + v * v - s(t, &rpar)) / (TWO * v);
    ydot.data[0] = ZERO;
    ydot.data[1] = e * tmp1 - tmp2 + sdot(t, &rpar) / (TWO * vtrue(t, &rpar));

    /* Return with success */
    0
}

/* fs routine to compute the slow portion of the ODE RHS. */
fn fs(t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let rpar = *get_rpar(user_data);
    let G = rpar[0];
    let e = rpar[2];
    let u = y.data[0];
    let v = y.data[1];

    /* fill in the RHS function:
       [G e]*[(-1+u^2-r(t))/(2*u))] + [rdot(t)/(2*u)]
       [0 0] [(-2+v^2-s(t))/(2*v)]    [      0      ] */
    let tmp1 = (-ONE + u * u - r(t, &rpar)) / (TWO * u);
    let tmp2 = (-TWO + v * v - s(t, &rpar)) / (TWO * v);
    ydot.data[0] = G * tmp1 + e * tmp2 + rdot(t, &rpar) / (TWO * u);
    ydot.data[1] = ZERO;

    /* Return with success */
    0
}

/* fse routine to compute the slow portion of the ODE RHS. */
fn fse(t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let rpar = *get_rpar(user_data);
    let u = y.data[0];

    /* fill in the slow explicit RHS function:
       [rdot(t)/(2*u)]
       [      0      ] */
    ydot.data[0] = rdot(t, &rpar) / (TWO * u);
    ydot.data[1] = ZERO;

    /* Return with success */
    0
}

/* fsi routine to compute the slow portion of the ODE RHS. */
fn fsi(t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let rpar = *get_rpar(user_data);
    let G = rpar[0];
    let e = rpar[2];
    let u = y.data[0];
    let v = y.data[1];

    /* fill in the slow implicit RHS function:
       [G e]*[(-1+u^2-r(t))/(2*u))]
       [0 0] [(-2+v^2-s(t))/(2*v)]  */
    let tmp1 = (-ONE + u * u - r(t, &rpar)) / (TWO * u);
    let tmp2 = (-TWO + v * v - s(t, &rpar)) / (TWO * v);
    ydot.data[0] = G * tmp1 + e * tmp2;
    ydot.data[1] = ZERO;

    /* Return with success */
    0
}

fn fnrhs(t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32 {
    let rpar = *get_rpar(user_data);
    let G = rpar[0];
    let e = rpar[2];
    let u = y.data[0];
    let v = y.data[1];

    /* fill in the RHS function:
       [G e]*[(-1+u^2-r(t))/(2*u))] + [rdot(t)/(2*u)]
       [e -1] [(-2+v^2-s(t))/(2*v)]   [sdot(t)/(2*vtrue(t))] */
    let tmp1 = (-ONE + u * u - r(t, &rpar)) / (TWO * u);
    let tmp2 = (-TWO + v * v - s(t, &rpar)) / (TWO * v);
    ydot.data[0] = G * tmp1 + e * tmp2 + rdot(t, &rpar) / (TWO * u);
    ydot.data[1] = e * tmp1 - tmp2 + sdot(t, &rpar) / (TWO * vtrue(t, &rpar));

    /* Return with success */
    0
}

fn f0(_t: f64, _y: &NVector, ydot: &mut NVector, _user_data: &mut UserData) -> i32 {
    N_VConst(ZERO, ydot);
    0
}

#[allow(clippy::too_many_arguments)]
fn Js(
    t: f64,
    y: &NVector,
    _fy: &NVector,
    j: &mut SUNMatrix,
    user_data: &mut UserData,
    _tmp1: &mut NVector,
    _tmp2: &mut NVector,
    _tmp3: &mut NVector,
) -> i32 {
    let rpar = *get_rpar(user_data);
    let G = rpar[0];
    let e = rpar[2];
    let u = y.data[0];
    let v = y.data[1];

    /* fill in the Jacobian:
       [G/2 + (G*(1+r(t))-rdot(t))/(2*u^2)   e/2+e*(2+s(t))/(2*v^2)]
       [                 0                             0           ] */
    if let SUNMatrix::Dense(dm) = j {
        /* (column-major SM_ELEMENT_D order) */
        dm.data[0] = G / TWO + (G * (ONE + r(t, &rpar)) - rdot(t, &rpar)) / (2.0 * u * u);
        dm.data[2] = e / TWO + e * (TWO + s(t, &rpar)) / (TWO * v * v);
        dm.data[1] = ZERO;
        dm.data[3] = ZERO;
    }

    /* Return with success */
    0
}

#[allow(clippy::too_many_arguments)]
fn Jsi(
    t: f64,
    y: &NVector,
    _fy: &NVector,
    j: &mut SUNMatrix,
    user_data: &mut UserData,
    _tmp1: &mut NVector,
    _tmp2: &mut NVector,
    _tmp3: &mut NVector,
) -> i32 {
    let rpar = *get_rpar(user_data);
    let G = rpar[0];
    let e = rpar[2];
    let u = y.data[0];
    let v = y.data[1];

    /* fill in the Jacobian:
       [G/2 + (G*(1+r(t)))/(2*u^2)   e/2 + e*(2+s(t))/(2*v^2)]
       [                 0                       0           ] */
    if let SUNMatrix::Dense(dm) = j {
        dm.data[0] = G / TWO + (G * (ONE + r(t, &rpar))) / (2.0 * u * u);
        dm.data[2] = e / TWO + e * (TWO + s(t, &rpar)) / (TWO * v * v);
        dm.data[1] = ZERO;
        dm.data[3] = ZERO;
    }

    /* Return with success */
    0
}

#[allow(clippy::too_many_arguments)]
fn Jn(
    t: f64,
    y: &NVector,
    _fy: &NVector,
    j: &mut SUNMatrix,
    user_data: &mut UserData,
    _tmp1: &mut NVector,
    _tmp2: &mut NVector,
    _tmp3: &mut NVector,
) -> i32 {
    let rpar = *get_rpar(user_data);
    let G = rpar[0];
    let e = rpar[2];
    let u = y.data[0];
    let v = y.data[1];

    /* fill in the Jacobian:
       [G/2 + (G*(1+r(t))-rdot(t))/(2*u^2)     e/2 + e*(2+s(t))/(2*v^2)]
       [e/2+e*(1+r(t))/(2*u^2)                -1/2 - (2+s(t))/(2*v^2)  ] */
    if let SUNMatrix::Dense(dm) = j {
        dm.data[0] = G / TWO + (G * (ONE + r(t, &rpar)) - rdot(t, &rpar)) / (2.0 * u * u);
        dm.data[2] = e / TWO + e * (TWO + s(t, &rpar)) / (TWO * v * v);
        dm.data[1] = e / TWO + e * (ONE + r(t, &rpar)) / (TWO * u * u);
        dm.data[3] = -ONE / TWO - (TWO + s(t, &rpar)) / (TWO * v * v);
    }

    /* Return with success */
    0
}

#[allow(clippy::too_many_arguments)]
fn Jf(
    t: f64,
    y: &NVector,
    _fy: &NVector,
    j: &mut SUNMatrix,
    user_data: &mut UserData,
    _tmp1: &mut NVector,
    _tmp2: &mut NVector,
    _tmp3: &mut NVector,
) -> i32 {
    let rpar = *get_rpar(user_data);
    let e = rpar[2];
    let u = y.data[0];
    let v = y.data[1];

    /* fill in the Jacobian:
       [        0                           0        ]
       [e/2+e*(1+r(t))/(2*u^2)  -1/2-(2+s(t))/(2*v^2)] */
    if let SUNMatrix::Dense(dm) = j {
        dm.data[0] = ZERO;
        dm.data[2] = ZERO;
        dm.data[1] = e / TWO + e * (ONE + r(t, &rpar)) / (TWO * u * u);
        dm.data[3] = -ONE / TWO - (TWO + s(t, &rpar)) / (TWO * v * v);
    }

    /* Return with success */
    0
}

/* ------------------------------
 * Private helper functions
 * ------------------------------*/

fn r(t: f64, _rpar: &[f64; 3]) -> f64 {
    0.5 * t.cos()
}

fn s(t: f64, rpar: &[f64; 3]) -> f64 {
    (rpar[1] * t).cos()
}

fn rdot(t: f64, _rpar: &[f64; 3]) -> f64 {
    -0.5 * t.sin()
}

fn sdot(t: f64, rpar: &[f64; 3]) -> f64 {
    -rpar[1] * (rpar[1] * t).sin()
}

fn utrue(t: f64, rpar: &[f64; 3]) -> f64 {
    (ONE + r(t, rpar)).sqrt()
}

fn vtrue(t: f64, rpar: &[f64; 3]) -> f64 {
    (TWO + s(t, rpar)).sqrt()
}

fn Ytrue(t: f64, y: &mut NVector, rpar: &[f64; 3]) -> i32 {
    y.data[0] = utrue(t, rpar);
    y.data[1] = vtrue(t, rpar);
    0
}

/* Check if a SUNDIALS function returned a negative flag */
fn check_retval(retval: i32, funcname: &str) -> bool {
    if retval < 0 {
        eprintln!("\nSUNDIALS_ERROR: {}() failed with retval = {}\n", funcname, retval);
        return true;
    }
    false
}

/* Main Program */
fn main() {
    /* general problem parameters */
    let T0: f64 = 0.0; /* initial time */
    let Tf: f64 = 5.0; /* final time */
    let dTout: f64 = 0.1; /* time between outputs */
    let NEQ: i64 = 2; /* number of dependent vars. */
    let Nt: i32 = (Tf / dTout).ceil() as i32; /* number of output times */
    let mut hs: f64 = 0.01; /* slow step size */
    let mut e: f64 = 0.5; /* fast/slow coupling strength */
    let mut G: f64 = -100.0; /* stiffness at slow time scale */
    let mut w: f64 = 100.0; /* time-scale separation factor */
    let mut reltol: f64 = 0.01;
    let mut abstol: f64 = 1.0e-11;

    let mut implicit_slow = false;
    let mut imex_slow = false;
    let mut explicit_slow = false;
    let mut no_slow = false;
    let mut implicit_fast = false;
    let mut explicit_fast = false;
    let mut no_fast = false;
    let mut deduce_rhs = false;

    /*
     * Initialization
     */

    /* Retrieve the command-line options: slow_type fast_type h G w e deduce_rhs */
    let argv: Vec<String> = std::env::args().collect();
    if argv.len() < 3 {
        println!("ERROR: executable requires at least two arguments [slow_type fast_type]");
        println!("Usage:");
        print!("  ark_kpr_mri slow_type fast_type h G w e deduce_rhs");
        std::process::exit(-1);
    }
    let slow_type: i32 = argv[1].parse().unwrap_or(0);
    let fast_type: i32 = argv[2].parse().unwrap_or(0);
    if argv.len() > 3 {
        hs = argv[3].parse().unwrap_or(0.0);
    }
    if argv.len() > 4 {
        G = argv[4].parse().unwrap_or(0.0);
    }
    if argv.len() > 5 {
        w = argv[5].parse().unwrap_or(0.0);
    }
    if argv.len() > 6 {
        e = argv[6].parse().unwrap_or(0.0);
    }
    if argv.len() > 7 {
        deduce_rhs = argv[7].parse::<i32>().unwrap_or(0) != 0;
    }

    /* Check arguments for validity */
    if !(0..=13).contains(&slow_type) {
        println!("ERROR: slow_type be an integer in [0,13] ");
        std::process::exit(-1);
    }
    if !(0..=5).contains(&fast_type) {
        println!("ERROR: fast_type be an integer in [0,5] ");
        std::process::exit(-1);
    }
    if slow_type == 0 && fast_type == 0 {
        println!("ERROR: at least one of slow_type and fast_type must be nonzero");
        std::process::exit(-1);
    }
    if slow_type >= 9 && fast_type == 0 {
        println!("ERROR: example not configured for ImEx slow solver with no fast solver");
        std::process::exit(-1);
    }
    if G >= ZERO {
        println!("ERROR: G must be a negative real number");
        std::process::exit(-1);
    }
    if hs <= ZERO {
        println!("ERROR: hs must be in positive");
        std::process::exit(-1);
    }
    /* (C evaluates this with implicit_slow still SUNFALSE) */
    if hs > ONE / G.abs() && !implicit_slow {
        println!("ERROR: hs must be in (0, 1/|G|)");
        std::process::exit(-1);
    }
    if w < ONE {
        println!("ERROR: w must be >= 1.0");
        std::process::exit(-1);
    }
    let rpar: [f64; 3] = [G, w, e];
    let hf = hs / w;

    /* Initial problem output (and set implicit solver tolerances as needed) */
    println!("\nMultirate nonlinear Kvaerno-Prothero-Robinson test problem:");
    println!("    time domain:  ({},{}]", fmt_g(T0, 0, 6), fmt_g(Tf, 0, 6));
    println!("    hs = {}", fmt_g(hs, 0, 6));
    println!("    hf = {}", fmt_g(hf, 0, 6));
    println!("    G = {}", fmt_g(G, 0, 6));
    println!("    w = {}", fmt_g(w, 0, 6));
    println!("    e = {}", fmt_g(e, 0, 6));
    match slow_type {
        0 => {
            println!("    slow solver: none");
            no_slow = true;
        }
        1 => {
            println!("    slow solver: ARKODE_MIS_KW3");
            explicit_slow = true;
        }
        2 => {
            println!("    slow solver: ARKODE_MRI_GARK_ERK45a");
            explicit_slow = true;
        }
        3 => {
            println!("    slow solver: ARKODE_MERK21");
            explicit_slow = true;
        }
        4 => {
            println!("    slow solver: ARKODE_MERK32");
            explicit_slow = true;
        }
        5 => {
            println!("    slow solver: ARKODE_MERK43");
            explicit_slow = true;
        }
        6 => {
            println!("    slow solver: ARKODE_MERK54");
            explicit_slow = true;
        }
        7 => {
            println!("    slow solver: ARKODE_MRI_GARK_IRK21a");
            implicit_slow = true;
            reltol = SUNMAX(hs * hs, 1.0e-10);
            abstol = 1.0e-11;
            println!("      reltol = {},  abstol = {}", fmt_e(reltol, 0, 2), fmt_e(abstol, 0, 2));
        }
        8 => {
            println!("    slow solver: ARKODE_MRI_GARK_ESDIRK34a");
            implicit_slow = true;
            reltol = SUNMAX(hs * hs * hs, 1.0e-10);
            abstol = 1.0e-11;
            println!("      reltol = {},  abstol = {}", fmt_e(reltol, 0, 2), fmt_e(abstol, 0, 2));
        }
        9 => {
            println!("    slow solver: ARKODE_IMEX_MRI_GARK3b");
            imex_slow = true;
            reltol = SUNMAX(hs * hs * hs, 1.0e-10);
            abstol = 1.0e-11;
            println!("      reltol = {},  abstol = {}", fmt_e(reltol, 0, 2), fmt_e(abstol, 0, 2));
        }
        10 => {
            println!("    slow solver: ARKODE_IMEX_MRI_GARK4");
            imex_slow = true;
            reltol = SUNMAX(hs * hs * hs * hs, 1.0e-14);
            abstol = 1.0e-14;
            println!("      reltol = {},  abstol = {}", fmt_e(reltol, 0, 2), fmt_e(abstol, 0, 2));
        }
        11 => {
            println!("    slow solver: ARKODE_IMEX_MRI_SR21");
            imex_slow = true;
            reltol = SUNMAX(hs * hs, 1.0e-10);
            abstol = 1.0e-11;
            println!("      reltol = {},  abstol = {}", fmt_e(reltol, 0, 2), fmt_e(abstol, 0, 2));
        }
        12 => {
            println!("    slow solver: ARKODE_IMEX_MRI_SR32");
            imex_slow = true;
            reltol = SUNMAX(hs * hs * hs, 1.0e-10);
            abstol = 1.0e-11;
            println!("      reltol = {},  abstol = {}", fmt_e(reltol, 0, 2), fmt_e(abstol, 0, 2));
        }
        13 => {
            println!("    slow solver: ARKODE_IMEX_MRI_SR43");
            imex_slow = true;
            reltol = SUNMAX(hs * hs * hs * hs, 1.0e-14);
            abstol = 1.0e-14;
            println!("      reltol = {},  abstol = {}", fmt_e(reltol, 0, 2), fmt_e(abstol, 0, 2));
        }
        _ => unreachable!(),
    }
    match fast_type {
        0 => {
            println!("    fast solver: none");
            no_fast = true;
        }
        1 => {
            println!("    fast solver: esdirk-3-3");
            implicit_fast = true;
            reltol = SUNMAX(hs * hs * hs, 1.0e-10);
            abstol = 1.0e-11;
            println!("      reltol = {},  abstol = {}", fmt_e(reltol, 0, 2), fmt_e(abstol, 0, 2));
        }
        2 => {
            println!("    fast solver: ARKODE_HEUN_EULER_2_1_2");
            explicit_fast = true;
        }
        3 => {
            println!("    fast solver: erk-3-3");
            explicit_fast = true;
        }
        4 => {
            println!("    fast solver: erk-4-4");
            explicit_fast = true;
        }
        5 => {
            println!("    fast solver: ARKODE_DORMAND_PRINCE_7_4_5");
            explicit_fast = true;
        }
        _ => unreachable!(),
    }

    /* Create the SUNDIALS context object for this simulation */
    let ctx = SUNContext_Create();

    /* Create and initialize serial vector for the solution */
    let mut y = N_VNew_Serial(NEQ, &ctx);
    let retval = Ytrue(T0, &mut y, &rpar);
    if check_retval(retval, "Ytrue") {
        return;
    }

    /*
     * Create the fast integrator and set options
     */

    /* Initialize the fast integrator.  If the fast scale is implicit, set up
       matrix, linear solver, and Jacobian function */
    let mut inner_arkode_mem = if no_fast {
        ARKStepCreate(Some(f0), None, T0, &y, &ctx).expect("ARKStepCreate")
    } else if explicit_fast && !no_slow {
        ARKStepCreate(Some(ff), None, T0, &y, &ctx).expect("ARKStepCreate")
    } else if explicit_fast && no_slow {
        ARKStepCreate(Some(fnrhs), None, T0, &y, &ctx).expect("ARKStepCreate")
    } else if implicit_fast && no_slow {
        let mut inner = ARKStepCreate(None, Some(fnrhs), T0, &y, &ctx).expect("ARKStepCreate");
        let Af = SUNDenseMatrix(NEQ, NEQ, &ctx);
        let LSf = SUNLinSol_Dense(&y, &Af, &ctx);
        let retval = ARKodeSetLinearSolver(&mut inner, LSf, Some(Af));
        if check_retval(retval, "ARKodeSetLinearSolver") {
            return;
        }
        let retval = ARKodeSetJacFn(&mut inner, Some(Jn));
        if check_retval(retval, "ARKodeSetJacFn") {
            return;
        }
        inner
    } else {
        /* implicit_fast && !no_slow */
        let mut inner = ARKStepCreate(None, Some(ff), T0, &y, &ctx).expect("ARKStepCreate");
        let Af = SUNDenseMatrix(NEQ, NEQ, &ctx);
        let LSf = SUNLinSol_Dense(&y, &Af, &ctx);
        let retval = ARKodeSetLinearSolver(&mut inner, LSf, Some(Af));
        if check_retval(retval, "ARKodeSetLinearSolver") {
            return;
        }
        let retval = ARKodeSetJacFn(&mut inner, Some(Jf));
        if check_retval(retval, "ARKodeSetJacFn") {
            return;
        }
        inner
    };

    /* Set Butcher table for fast integrator */
    match fast_type {
        0 => {
            let mut B = ARKodeButcherTable_Alloc(3, true).expect("ARKodeButcherTable_Alloc");
            B.A[1][0] = 0.5;
            B.A[2][0] = -ONE;
            B.A[2][1] = TWO;
            B.b[0] = ONE / 6.0;
            B.b[1] = TWO / 3.0;
            B.b[2] = ONE / 6.0;
            B.d.as_mut().unwrap()[1] = ONE;
            B.c[1] = 0.5;
            B.c[2] = ONE;
            B.q = 3;
            B.p = 2;
            let retval = ARKStepSetTables(&mut inner_arkode_mem, 3, 2, None, Some(&B));
            if check_retval(retval, "ARKStepSetTables") {
                return;
            }
        }
        1 => {
            let mut B = ARKodeButcherTable_Alloc(3, false).expect("ARKodeButcherTable_Alloc");
            let beta = (3.0f64).sqrt() / 6.0 + 0.5;
            let gamma = (-ONE / 8.0) * ((3.0f64).sqrt() + ONE);
            B.A[1][0] = 4.0 * gamma + TWO * beta;
            B.A[1][1] = ONE - 4.0 * gamma - TWO * beta;
            B.A[2][0] = 0.5 - beta - gamma;
            B.A[2][1] = gamma;
            B.A[2][2] = beta;
            B.b[0] = ONE / 6.0;
            B.b[1] = ONE / 6.0;
            B.b[2] = TWO / 3.0;
            B.c[1] = ONE;
            B.c[2] = 0.5;
            B.q = 3;
            let retval = ARKStepSetTables(&mut inner_arkode_mem, 3, 0, Some(&B), None);
            if check_retval(retval, "ARKStepSetTables") {
                return;
            }
        }
        2 => {
            let B = ARKodeButcherTable_LoadERK(ARKODE_HEUN_EULER_2_1_2)
                .expect("ARKodeButcherTable_LoadERK");
            let retval = ARKStepSetTables(&mut inner_arkode_mem, 2, 1, None, Some(&B));
            if check_retval(retval, "ARKStepSetTables") {
                return;
            }
        }
        3 => {
            let mut B = ARKodeButcherTable_Alloc(3, true).expect("ARKodeButcherTable_Alloc");
            B.A[1][0] = 0.5;
            B.A[2][0] = -ONE;
            B.A[2][1] = TWO;
            B.b[0] = ONE / 6.0;
            B.b[1] = TWO / 3.0;
            B.b[2] = ONE / 6.0;
            B.d.as_mut().unwrap()[1] = ONE;
            B.c[1] = 0.5;
            B.c[2] = ONE;
            B.q = 3;
            B.p = 2;
            let retval = ARKStepSetTables(&mut inner_arkode_mem, 3, 2, None, Some(&B));
            if check_retval(retval, "ARKStepSetTables") {
                return;
            }
        }
        4 => {
            let mut B = ARKodeButcherTable_Alloc(4, false).expect("ARKodeButcherTable_Alloc");
            B.A[1][0] = 0.5;
            B.A[2][1] = 0.5;
            B.A[3][2] = ONE;
            B.b[0] = ONE / 6.0;
            B.b[1] = ONE / 3.0;
            B.b[2] = ONE / 3.0;
            B.b[3] = ONE / 6.0;
            B.c[1] = 0.5;
            B.c[2] = 0.5;
            B.c[3] = ONE;
            B.q = 4;
            let retval = ARKStepSetTables(&mut inner_arkode_mem, 4, 0, None, Some(&B));
            if check_retval(retval, "ARKStepSetTables") {
                return;
            }
        }
        5 => {
            let B = ARKodeButcherTable_LoadERK(ARKODE_DORMAND_PRINCE_7_4_5)
                .expect("ARKodeButcherTable_LoadERK");
            let retval = ARKStepSetTables(&mut inner_arkode_mem, 5, 4, None, Some(&B));
            if check_retval(retval, "ARKStepSetTables") {
                return;
            }
        }
        _ => unreachable!(),
    }

    /* Set the tolerances */
    let retval = ARKodeSStolerances(&mut inner_arkode_mem, reltol, abstol);
    if check_retval(retval, "ARKodeSStolerances") {
        return;
    }

    /* Set the user data pointer */
    let retval = ARKodeSetUserData(&mut inner_arkode_mem, Some(Box::new(rpar)));
    if check_retval(retval, "ARKodeSetUserData") {
        return;
    }

    /* Set the fast step size */
    let retval = ARKodeSetFixedStep(&mut inner_arkode_mem, hf);
    if check_retval(retval, "ARKodeSetFixedStep") {
        return;
    }

    /* Override any current settings with command-line options -- enforce
       the prefix "inner" */
    let retval = ARKodeSetOptions(&mut inner_arkode_mem, Some("inner"), Some(""), &argv);
    if check_retval(retval, "ARKodeSetOptions") {
        return;
    }

    /* Create inner stepper */
    let inner_stepper =
        ARKodeCreateMRIStepInnerStepper(inner_arkode_mem).expect("ARKodeCreateMRIStepInnerStepper");

    /*
     * Create the slow integrator and set options
     */

    /* Initialize the slow integrator.  If the slow scale contains an implicit
       component, set up matrix, linear solver, and Jacobian function. */
    let mut arkode_mem = if no_slow {
        MRIStepCreate(Some(f0), None, T0, &y, inner_stepper, &ctx).expect("MRIStepCreate")
    } else if explicit_slow && !no_fast {
        MRIStepCreate(Some(fs), None, T0, &y, inner_stepper, &ctx).expect("MRIStepCreate")
    } else if explicit_slow && no_fast {
        MRIStepCreate(Some(fnrhs), None, T0, &y, inner_stepper, &ctx).expect("MRIStepCreate")
    } else if implicit_slow && !no_fast {
        let mut outer =
            MRIStepCreate(None, Some(fs), T0, &y, inner_stepper, &ctx).expect("MRIStepCreate");
        let As = SUNDenseMatrix(NEQ, NEQ, &ctx);
        let LSs = SUNLinSol_Dense(&y, &As, &ctx);
        let retval = ARKodeSetLinearSolver(&mut outer, LSs, Some(As));
        if check_retval(retval, "ARKodeSetLinearSolver") {
            return;
        }
        let retval = ARKodeSetJacFn(&mut outer, Some(Js));
        if check_retval(retval, "ARKodeSetJacFn") {
            return;
        }
        outer
    } else if implicit_slow && no_fast {
        let mut outer =
            MRIStepCreate(None, Some(fnrhs), T0, &y, inner_stepper, &ctx).expect("MRIStepCreate");
        let As = SUNDenseMatrix(NEQ, NEQ, &ctx);
        let LSs = SUNLinSol_Dense(&y, &As, &ctx);
        let retval = ARKodeSetLinearSolver(&mut outer, LSs, Some(As));
        if check_retval(retval, "ARKodeSetLinearSolver") {
            return;
        }
        let retval = ARKodeSetJacFn(&mut outer, Some(Jn));
        if check_retval(retval, "ARKodeSetJacFn") {
            return;
        }
        outer
    } else {
        /* imex_slow */
        let mut outer = MRIStepCreate(Some(fse), Some(fsi), T0, &y, inner_stepper, &ctx)
            .expect("MRIStepCreate");
        let As = SUNDenseMatrix(NEQ, NEQ, &ctx);
        let LSs = SUNLinSol_Dense(&y, &As, &ctx);
        let retval = ARKodeSetLinearSolver(&mut outer, LSs, Some(As));
        if check_retval(retval, "ARKodeSetLinearSolver") {
            return;
        }
        let retval = ARKodeSetJacFn(&mut outer, Some(Jsi));
        if check_retval(retval, "ARKodeSetJacFn") {
            return;
        }
        outer
    };

    /* Set coupling table for slow integrator */
    let C = match slow_type {
        0 => {
            /* no slow dynamics (use ERK-2-2) */
            let mut B = ARKodeButcherTable_Alloc(2, false).expect("ARKodeButcherTable_Alloc");
            B.A[1][0] = TWO / 3.0;
            B.b[0] = 0.25;
            B.b[1] = 0.75;
            B.c[1] = TWO / 3.0;
            B.q = 2;
            MRIStepCoupling_MIStoMRI(&B, 2, 0).expect("MRIStepCoupling_MIStoMRI")
        }
        1 => MRIStepCoupling_LoadTable(ARKODE_MIS_KW3).expect("MRIStepCoupling_LoadTable"),
        2 => MRIStepCoupling_LoadTable(ARKODE_MRI_GARK_ERK45a).expect("MRIStepCoupling_LoadTable"),
        3 => MRIStepCoupling_LoadTable(ARKODE_MERK21).expect("MRIStepCoupling_LoadTable"),
        4 => MRIStepCoupling_LoadTable(ARKODE_MERK32).expect("MRIStepCoupling_LoadTable"),
        5 => MRIStepCoupling_LoadTable(ARKODE_MERK43).expect("MRIStepCoupling_LoadTable"),
        6 => MRIStepCoupling_LoadTable(ARKODE_MERK54).expect("MRIStepCoupling_LoadTable"),
        7 => MRIStepCoupling_LoadTable(ARKODE_MRI_GARK_IRK21a).expect("MRIStepCoupling_LoadTable"),
        8 => MRIStepCoupling_LoadTable(ARKODE_MRI_GARK_ESDIRK34a)
            .expect("MRIStepCoupling_LoadTable"),
        9 => MRIStepCoupling_LoadTable(ARKODE_IMEX_MRI_GARK3b).expect("MRIStepCoupling_LoadTable"),
        10 => MRIStepCoupling_LoadTable(ARKODE_IMEX_MRI_GARK4).expect("MRIStepCoupling_LoadTable"),
        11 => MRIStepCoupling_LoadTable(ARKODE_IMEX_MRI_SR21).expect("MRIStepCoupling_LoadTable"),
        12 => MRIStepCoupling_LoadTable(ARKODE_IMEX_MRI_SR32).expect("MRIStepCoupling_LoadTable"),
        13 => MRIStepCoupling_LoadTable(ARKODE_IMEX_MRI_SR43).expect("MRIStepCoupling_LoadTable"),
        _ => unreachable!(),
    };
    let retval = MRIStepSetCoupling(&mut arkode_mem, &C);
    if check_retval(retval, "MRIStepSetCoupling") {
        return;
    }
    drop(C); /* free coupling coefficients */

    /* Set the tolerances */
    let retval = ARKodeSStolerances(&mut arkode_mem, reltol, abstol);
    if check_retval(retval, "ARKodeSStolerances") {
        return;
    }

    /* Set the user data pointer */
    let retval = ARKodeSetUserData(&mut arkode_mem, Some(Box::new(rpar)));
    if check_retval(retval, "ARKodeSetUserData") {
        return;
    }

    let retval = ARKodeSetDeduceImplicitRhs(&mut arkode_mem, deduce_rhs);
    if check_retval(retval, "ARKodeSetDeduceImplicitRhs") {
        return;
    }

    /* Set the slow step size */
    let retval = ARKodeSetFixedStep(&mut arkode_mem, hs);
    if check_retval(retval, "ARKodeSetFixedStep") {
        return;
    }

    /* Override any current settings with command-line options -- enforce
       the prefix "outer" */
    let retval = ARKodeSetOptions(&mut arkode_mem, Some("outer"), Some(""), &argv);
    if check_retval(retval, "ARKodeSetOptions") {
        return;
    }

    /*
     * Integrate ODE
     */

    /* Open output stream for results, output comment line */
    let mut ufid = std::fs::File::create("ark_kpr_mri_solution.txt").expect("fopen");
    let _ = writeln!(ufid, "# t u v uerr verr");

    /* output initial condition to disk */
    let _ = writeln!(
        ufid,
        " {} {} {} {} {}",
        fmt_e(T0, 0, 16),
        fmt_e(y.data[0], 0, 16),
        fmt_e(y.data[1], 0, 16),
        fmt_e((y.data[0] - utrue(T0, &rpar)).abs(), 0, 16),
        fmt_e((y.data[1] - vtrue(T0, &rpar)).abs(), 0, 16)
    );

    /* Main time-stepping loop: calls ARKodeEvolve to perform the
       integration, then prints results. Stops when the final time
       has been reached */
    let mut t = T0;
    let mut tout = T0 + dTout;
    let mut uerr = ZERO;
    let mut verr = ZERO;
    let mut uerrtot = ZERO;
    let mut verrtot = ZERO;
    let mut errtot = ZERO;
    println!("        t           u           v       uerr      verr");
    println!("   ------------------------------------------------------");
    println!(
        "  {}  {}  {}  {}  {}",
        fmt_f(t, 10, 6),
        fmt_f(y.data[0], 10, 6),
        fmt_f(y.data[1], 10, 6),
        fmt_e(uerr, 0, 2),
        fmt_e(verr, 0, 2)
    );

    for _iout in 0..Nt {
        /* call integrator */
        let retval = ARKodeEvolve(&mut arkode_mem, tout, &mut y, &mut t, ARK_NORMAL);
        if check_retval(retval, "ARKodeEvolve") {
            break;
        }

        /* access/print solution and error */
        uerr = (y.data[0] - utrue(t, &rpar)).abs();
        verr = (y.data[1] - vtrue(t, &rpar)).abs();
        println!(
            "  {}  {}  {}  {}  {}",
            fmt_f(t, 10, 6),
            fmt_f(y.data[0], 10, 6),
            fmt_f(y.data[1], 10, 6),
            fmt_e(uerr, 0, 2),
            fmt_e(verr, 0, 2)
        );
        let _ = writeln!(
            ufid,
            " {} {} {} {} {}",
            fmt_e(t, 0, 16),
            fmt_e(y.data[0], 0, 16),
            fmt_e(y.data[1], 0, 16),
            fmt_e(uerr, 0, 16),
            fmt_e(verr, 0, 16)
        );
        uerrtot += uerr * uerr;
        verrtot += verr * verr;
        errtot += uerr * uerr + verr * verr;

        /* successful solve: update time */
        tout += dTout;
        tout = if tout > Tf { Tf } else { tout };
    }
    uerrtot = (uerrtot / Nt as f64).sqrt();
    verrtot = (verrtot / Nt as f64).sqrt();
    errtot = (errtot / Nt as f64 / 2.0).sqrt();
    println!("   ------------------------------------------------------");
    drop(ufid);

    /*
     * Finalize
     */

    /* Get some slow integrator statistics */
    let mut nsts: i64 = 0;
    let mut nfse_c: i64 = 0;
    let mut nfsi_c: i64 = 0;
    ARKodeGetNumSteps(&mut arkode_mem, &mut nsts);
    ARKodeGetNumRhsEvals(&mut arkode_mem, 0, &mut nfse_c);
    ARKodeGetNumRhsEvals(&mut arkode_mem, 1, &mut nfsi_c);

    /* Get some fast integrator statistics (borrow the inner integrator
       back out of the outer step memory) */
    let mut nstf: i64 = 0;
    let mut nff: i64 = 0;
    let mut nnif: i64 = 0;
    let mut nncf: i64 = 0;
    let mut njef: i64 = 0;
    {
        let step_mem = arkode_mem
            .step_mem
            .as_mut()
            .unwrap()
            .downcast_mut::<ARKodeMRIStepMem>()
            .unwrap();
        let inner = step_mem
            .stepper
            .content
            .as_mut()
            .unwrap()
            .downcast_mut::<ARKodeMem>()
            .unwrap();
        ARKodeGetNumSteps(inner, &mut nstf);
        ARKodeGetNumRhsEvals(inner, 0, &mut nff);
        if implicit_fast {
            ARKodeGetNonlinSolvStats(inner, &mut nnif, &mut nncf);
            ARKodeGetNumJacEvals(inner, &mut njef);
        }
    }

    /* Print some final statistics */
    println!("\nFinal Solver Statistics:");
    println!("   Steps: nsts = {}, nstf = {}", nsts, nstf);
    println!(
        "   u error = {}, v error = {}, total error = {}",
        fmt_e(uerrtot, 0, 3),
        fmt_e(verrtot, 0, 3),
        fmt_e(errtot, 0, 3)
    );
    if imex_slow {
        println!("   Total RHS evals:  Fse = {}, Fsi = {},  Ff = {}", nfse_c, nfsi_c, nff);
    } else if implicit_slow {
        println!("   Total RHS evals:  Fs = {},  Ff = {}", nfsi_c, nff);
    } else {
        println!("   Total RHS evals:  Fs = {},  Ff = {}", nfse_c, nff);
    }

    /* Get/print slow integrator decoupled implicit solver statistics */
    if implicit_slow || imex_slow {
        let mut nnis: i64 = 0;
        let mut nncs: i64 = 0;
        let mut njes: i64 = 0;
        ARKodeGetNonlinSolvStats(&mut arkode_mem, &mut nnis, &mut nncs);
        ARKodeGetNumJacEvals(&mut arkode_mem, &mut njes);
        println!("   Slow Newton iters = {}", nnis);
        println!("   Slow Newton conv fails = {}", nncs);
        println!("   Slow Jacobian evals = {}", njes);
    }

    /* Print fast integrator implicit solver statistics */
    if implicit_fast {
        println!("   Fast Newton iters = {}", nnif);
        println!("   Fast Newton conv fails = {}", nncf);
        println!("   Fast Jacobian evals = {}", njef);
    }

    /* Clean up and return */
    drop(y); /* Free y vector */
    let mut slot = Some(arkode_mem);
    ARKodeFree(&mut slot);
}
