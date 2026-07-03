/* -----------------------------------------------------------------
 * solar_system — new example for the cvode_rs pure-Rust port.
 *
 * Newtonian N-body integration of the outer solar system:
 * the Sun (with the inner planets lumped in) plus Jupiter, Saturn,
 * Uranus, Neptune and Pluto. This is the classical DETEST-style
 * "outer planet" problem: positions in AU, time in days, masses
 * relative to the Sun, G = 2.95912208286e-4 AU^3 / (day^2 M_sun).
 *
 *   dq_i/dt = v_i
 *   dv_i/dt = G * sum_{j != i} m_j (q_j - q_i) / |q_j - q_i|^3
 *
 * 6 bodies x 3 coordinates x (position + velocity) = 36 equations.
 * Integrated with CVODE Adams (non-stiff) + Newton iteration and the
 * dense linear solver with an internal difference-quotient Jacobian,
 * over 500 000 days (~1369 years, about 5.5 Pluto orbits).
 * Energy conservation is reported as the correctness check.
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]

use cvode_rs::sundials_utils::{fmt_e, fmt_f};
use cvode_rs::*;

const NBODY: usize = 6;
const NEQ: usize = 6 * NBODY; /* 36 */

const G: f64 = 2.95912208286e-4;

/* masses relative to the sun (sun includes the inner planets) */
const MASS: [f64; NBODY] = [
    1.00000597682,     /* Sun + inner planets */
    9.54786104043e-4,  /* Jupiter */
    2.85583733151e-4,  /* Saturn  */
    4.37273164546e-5,  /* Uranus  */
    5.17759138449e-5,  /* Neptune */
    1.0 / 1.3e8,       /* Pluto   */
];

/* initial positions (AU) and velocities (AU/day), J2000-era epoch
   (E. Hairer, S.P. Norsett, G. Wanner, "Solving ODEs I", Sect. I.1) */
const Q0: [[f64; 3]; NBODY] = [
    [0.0, 0.0, 0.0],
    [-3.5023653, -3.8169847, -1.5507963],
    [9.0755314, -3.0458353, -1.6483708],
    [8.3101420, -16.2901086, -7.2521278],
    [11.4707666, -25.7294829, -10.8169456],
    [-15.5387357, -25.2225594, -3.1902382],
];
const V0: [[f64; 3]; NBODY] = [
    [0.0, 0.0, 0.0],
    [0.00565429, -0.00412490, -0.00190589],
    [0.00168318, 0.00483525, 0.00192462],
    [0.00354178, 0.00137102, 0.00055029],
    [0.00288930, 0.00114527, 0.00039677],
    [0.00276725, -0.00170702, -0.00136504],
];

const T0: f64 = 0.0;
const TOUT_STEP: f64 = 50_000.0; /* print every 50 000 days */
const NOUT: i32 = 10; /* to t = 500 000 days */

const RTOL: f64 = 1.0e-10;
const ATOL: f64 = 1.0e-12;

/* state layout: y = [q1x q1y q1z v1x v1y v1z | q2x ... ] */
fn f(_t: f64, y: &NVector, ydot: &mut NVector, _user_data: &mut UserData) -> i32 {
    let d = &y.data;
    let out = &mut ydot.data;

    /* dq/dt = v */
    for i in 0..NBODY {
        for k in 0..3 {
            out[6 * i + k] = d[6 * i + 3 + k];
        }
    }

    /* dv/dt = gravitational acceleration */
    for i in 0..NBODY {
        let mut acc = [0.0f64; 3];
        for j in 0..NBODY {
            if i == j {
                continue;
            }
            let dx = d[6 * j] - d[6 * i];
            let dy = d[6 * j + 1] - d[6 * i + 1];
            let dz = d[6 * j + 2] - d[6 * i + 2];
            let r2 = dx * dx + dy * dy + dz * dz;
            let r3 = r2 * r2.sqrt();
            let c = G * MASS[j] / r3;
            acc[0] += c * dx;
            acc[1] += c * dy;
            acc[2] += c * dz;
        }
        for k in 0..3 {
            out[6 * i + 3 + k] = acc[k];
        }
    }

    0
}

/* total energy: kinetic + potential (conserved quantity) */
fn energy(y: &NVector) -> f64 {
    let d = &y.data;
    let mut e = 0.0;
    for i in 0..NBODY {
        let v2 = d[6 * i + 3] * d[6 * i + 3]
            + d[6 * i + 4] * d[6 * i + 4]
            + d[6 * i + 5] * d[6 * i + 5];
        e += 0.5 * MASS[i] * v2;
        for j in (i + 1)..NBODY {
            let dx = d[6 * j] - d[6 * i];
            let dy = d[6 * j + 1] - d[6 * i + 1];
            let dz = d[6 * j + 2] - d[6 * i + 2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            e -= G * MASS[i] * MASS[j] / r;
        }
    }
    e
}

fn main() {
    let sunctx = SUNContext_Create();

    /* initial state */
    let mut y = N_VNew_Serial(NEQ as i64, &sunctx);
    for i in 0..NBODY {
        for k in 0..3 {
            y.data[6 * i + k] = Q0[i][k];
            y.data[6 * i + 3 + k] = V0[i][k];
        }
    }

    let e0 = energy(&y);

    /* Adams method (non-stiff problem), Newton + dense DQ Jacobian */
    let mut cvode_mem = CVodeCreate(CV_ADAMS, &sunctx);

    let mut retval = CVodeInit(&mut cvode_mem, f, T0, &y);
    assert!(retval == CV_SUCCESS, "CVodeInit failed: {retval}");

    retval = CVodeSStolerances(&mut cvode_mem, RTOL, ATOL);
    assert!(retval == CV_SUCCESS, "CVodeSStolerances failed: {retval}");

    let a = SUNDenseMatrix(NEQ as i64, NEQ as i64, &sunctx);
    let ls = SUNLinSol_Dense(&y, &a, &sunctx);
    retval = CVodeSetLinearSolver(&mut cvode_mem, ls, Some(a));
    assert!(retval == CV_SUCCESS, "CVodeSetLinearSolver failed: {retval}");

    retval = CVodeSetMaxNumSteps(&mut cvode_mem, 500_000);
    assert!(retval == CV_SUCCESS);

    println!("\nOuter solar system N-body problem (6 bodies, 36 equations)");
    println!("Adams-Moulton + Newton + dense DQ Jacobian, rtol = {}, atol = {}\n", RTOL, ATOL);
    println!(
        "{:>10}  {:>14} {:>14} {:>14}  {:>16}",
        "t (days)", "Pluto x (AU)", "Pluto y (AU)", "Pluto z (AU)", "energy drift"
    );

    let mut t = T0;
    let mut tout = TOUT_STEP;
    for _ in 0..NOUT {
        retval = CVode(&mut cvode_mem, tout, &mut y, &mut t, CV_NORMAL);
        if retval < 0 {
            eprintln!("CVode failed with retval = {retval}");
            std::process::exit(1);
        }
        let e = energy(&y);
        let drift = ((e - e0) / e0).abs();
        println!(
            "{:>10}  {} {} {}  {:>16}",
            fmt_f(t, 10, 0),
            fmt_f(y.data[6 * 5], 14, 8),
            fmt_f(y.data[6 * 5 + 1], 14, 8),
            fmt_f(y.data[6 * 5 + 2], 14, 8),
            fmt_e(drift, 16, 6)
        );
        tout += TOUT_STEP;
    }

    /* final statistics */
    let mut nst = 0i64;
    let mut nfe = 0i64;
    let mut nni = 0i64;
    let mut netf = 0i64;
    CVodeGetNumSteps(&mut cvode_mem, &mut nst);
    CVodeGetNumRhsEvals(&mut cvode_mem, &mut nfe);
    CVodeGetNumNonlinSolvIters(&mut cvode_mem, &mut nni);
    CVodeGetNumErrTestFails(&mut cvode_mem, &mut netf);

    println!("\nFinal Statistics:");
    println!("  internal steps     = {nst}");
    println!("  rhs evaluations    = {nfe}");
    println!("  nonlinear iters    = {nni}");
    println!("  error test fails   = {netf}");

    let e = energy(&y);
    let drift = ((e - e0) / e0).abs();
    println!("  relative energy drift over 500000 days = {}", fmt_e(drift, 0, 6));

    /* correctness check: energy must be conserved to ~1e-6 relative */
    if drift < 1.0e-6 {
        println!("\nSUCCESS: energy conserved within tolerance");
    } else {
        println!("\nFAILURE: energy drift too large");
        std::process::exit(1);
    }

    CVodeFree(cvode_mem);
}
