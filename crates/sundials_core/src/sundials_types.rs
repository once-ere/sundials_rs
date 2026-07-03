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

/// SUNOutputFormat (sundials_types.h)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SUNOutputFormat {
    SUN_OUTPUTFORMAT_TABLE,
    SUN_OUTPUTFORMAT_CSV,
}
pub use SUNOutputFormat::{SUN_OUTPUTFORMAT_CSV, SUN_OUTPUTFORMAT_TABLE};
