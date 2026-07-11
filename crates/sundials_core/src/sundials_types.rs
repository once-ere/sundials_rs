/* -----------------------------------------------------------------
 * Translated from include/sundials/sundials_types.h (SUNDIALS 7.7.0)
 * Basic scalar types and floating-point constants.
 * -----------------------------------------------------------------*/

/// `sunrealtype` — SUNDIALS was configured with double precision.
pub type sunrealtype = f64;
/// `sunindextype` — 64-bit signed indices.
pub type sunindextype = i64;
/// `sunbooleantype`
pub type sunbooleantype = bool;
/// `suncountertype` — SUNDIALS_COUNTER_TYPE (long long, 64-bit counters).
pub type suncountertype = i64;

pub const SUNFALSE: bool = false;
pub const SUNTRUE: bool = true;

/// SUN_BIG_REAL — largest normalized double (BIG_REAL)
pub const SUN_BIG_REAL: f64 = f64::MAX;
/// SUN_SMALL_REAL — smallest normalized positive double
pub const SUN_SMALL_REAL: f64 = f64::MIN_POSITIVE;
/// SUN_UNIT_ROUNDOFF — DBL_EPSILON
pub const SUN_UNIT_ROUNDOFF: f64 = f64::EPSILON;

/// user_data passed through to all user callbacks (C `void*`).
pub type UserData = Option<Box<dyn std::any::Any>>;

/// FSA user-data convention for the internal difference-quotient
/// sensitivity residuals (CVODES cvSensRhs1InternalDQ / IDAS
/// IDASensRes1DQ and the QuadSens analogues).
///
/// In C, `CVodeSetSensParams`/`IDASetSensParams` store the user's `p`
/// POINTER, and the internal DQ routines perturb `p[which]` in place —
/// the perturbation reaches the user's RHS/residual because `p` aliases
/// the user's own parameter array inside their user data.  The Rust
/// ports keep `cv_p`/`ida_p` as owned copies (no aliasing into user
/// data is expressible), so user code that relies on the INTERNAL DQ
/// sensitivities (fS/resS = None) must store its parameter array in an
/// `FSAUserData` wrapper: the DQ routines downcast the user data to
/// this type and mirror each `p[which]` perturbation into `.p` so the
/// user callback observes it (and restore it afterwards, as C does).
/// Problem constants beyond the parameters live in `.user`.
///
/// User code with an analytic sensitivity residual does not need this
/// wrapper.
pub struct FSAUserData {
    /// the problem parameters (C: the array `p` handed to
    /// CVodeSetSensParams/IDASetSensParams)
    pub p: Vec<f64>,
    /// remaining user data (C: the rest of the user's structure)
    pub user: Box<dyn std::any::Any>,
}

/// SUNOutputFormat (sundials_types.h)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SUNOutputFormat {
    SUN_OUTPUTFORMAT_TABLE,
    SUN_OUTPUTFORMAT_CSV,
}
pub use SUNOutputFormat::{SUN_OUTPUTFORMAT_CSV, SUN_OUTPUTFORMAT_TABLE};
