/* -----------------------------------------------------------------
 * Translation of sundials-7.7.0/src/arkode/arkode_interp_impl.h.
 * Content structures for the two ARKODE temporal interpolation
 * implementations. The generic ARKInterp object (declared in
 * arkode_impl.h) is the ARKInterp enum in arkode_impl.rs; the
 * operations (arkInterpCreate_Hermite / _Lagrange etc.) live in
 * arkode_interp.rs (from arkode_interp.c). The HINT_ / LINT_
 * accessor macros become direct field access on the enum variants.
 * -----------------------------------------------------------------*/

use crate::nvector_serial::NVector;

/* Numeric constants */
pub const FOURTH: f64 = 0.25;
pub const THREE: f64 = 3.0;
pub const SIX: f64 = 6.0;
pub const TWELVE: f64 = 12.0;

/// struct _ARKInterpContent_Hermite
pub struct ARKInterpContent_Hermite {
    pub degree: i32,      /* maximum interpolant degree to use           */
    pub fold: NVector,    /* f(t,y) at beginning of last successful step */
    pub yold: NVector,    /* y at beginning of last successful step      */
    pub fa: NVector,      /* f(t,y) used in higher-order interpolation   */
    pub fb: NVector,      /* f(t,y) used in higher-order interpolation   */
    pub told: f64,        /* t at beginning of last successful step      */
    pub tnew: f64,        /* t at end of last successful step            */
    pub h: f64,           /* last successful step size                   */
}

/// struct _ARKInterpContent_Lagrange
pub struct ARKInterpContent_Lagrange {
    pub nmax: i32,         /* number of previous solutions to use      */
    pub nmaxalloc: i32,    /* vectors allocated for previous solutions */
    pub yhist: Vec<NVector>, /* previous solution vectors              */
    pub thist: Vec<f64>,   /* 't' values associated with yhist         */
    pub nhist: i32,        /* number of 'active' vectors in yhist      */
    pub tround: f64,       /* unit roundoff for 't' values             */
}
