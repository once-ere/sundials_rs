/* -----------------------------------------------------------------
 * Counterpart of src/sundials/sundials_utils.h (SUNDIALS 7.7.0):
 * small formatting helpers. In C these are snprintf wrappers
 * (sunsnprintf, SUN_FORMAT_G, ...). Here they reproduce C printf
 * float conversions exactly so translated examples match the
 * reference .out files byte-for-byte.
 * -----------------------------------------------------------------*/

/// C `%.*e` conversion (two-digit minimum exponent, as on macOS/Linux).
pub fn fmt_e(x: f64, width: usize, prec: usize) -> String {
    let s = fmt_e_core(x, prec);
    pad(s, width)
}

fn fmt_e_core(x: f64, prec: usize) -> String {
    if x.is_nan() {
        return "nan".to_string();
    }
    if x.is_infinite() {
        return if x < 0.0 { "-inf".into() } else { "inf".into() };
    }
    // Rust {:e} gives e.g. "1.5e-3" / "-2.25e10"; normalize to C form.
    let s = format!("{:.*e}", prec, x);
    let (mant, exp) = s.split_once('e').expect("exponent form");
    let expv: i32 = exp.parse().expect("exponent integer");
    let sign = if expv < 0 { '-' } else { '+' };
    format!("{}e{}{:02}", mant, sign, expv.abs())
}

/// C `%.*f` conversion.
pub fn fmt_f(x: f64, width: usize, prec: usize) -> String {
    pad(format!("{:.*}", prec, x), width)
}

/// C `%.*g` conversion (SUN_FORMAT_G uses %.17g; examples use %g).
pub fn fmt_g(x: f64, width: usize, prec: usize) -> String {
    pad(fmt_g_core(x, prec), width)
}

fn fmt_g_core(x: f64, prec: usize) -> String {
    if x.is_nan() {
        return "nan".to_string();
    }
    if x.is_infinite() {
        return if x < 0.0 { "-inf".into() } else { "inf".into() };
    }
    let p = if prec == 0 { 1 } else { prec };
    if x == 0.0 {
        return "0".to_string();
    }
    // decimal exponent X of the %e representation rounded to p-1 digits
    let e_repr = format!("{:.*e}", p - 1, x);
    let expv: i32 = e_repr
        .split_once('e')
        .map(|(_, e)| e.parse().unwrap())
        .unwrap();
    let out = if expv >= -4 && expv < p as i32 {
        // %f style with precision p-1-X, then strip trailing zeros
        let fprec = (p as i32 - 1 - expv).max(0) as usize;
        strip_zeros(format!("{:.*}", fprec, x))
    } else {
        // %e style with precision p-1, strip zeros in mantissa
        let (mant, exp) = e_repr.split_once('e').unwrap();
        let mant = strip_zeros(mant.to_string());
        let ev: i32 = exp.parse().unwrap();
        let sign = if ev < 0 { '-' } else { '+' };
        format!("{}e{}{:02}", mant, sign, ev.abs())
    };
    out
}

fn strip_zeros(mut s: String) -> String {
    if s.contains('.') {
        while s.ends_with('0') {
            s.pop();
        }
        if s.ends_with('.') {
            s.pop();
        }
    }
    s
}

fn pad(s: String, width: usize) -> String {
    if s.len() >= width {
        s
    } else {
        format!("{:>width$}", s, width = width)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn e_matches_c_printf() {
        assert_eq!(fmt_e(0.4, 12, 4), "  4.0000e-01");
        assert_eq!(fmt_e(-1.234567e-5, 0, 6), "-1.234567e-05");
        assert_eq!(fmt_e(2.6391e2, 0, 4), "2.6391e+02");
        assert_eq!(fmt_e(0.0, 0, 2), "0.00e+00");
        assert_eq!(fmt_e(1e100, 0, 3), "1.000e+100");
    }

    #[test]
    fn g_matches_c_printf() {
        assert_eq!(fmt_g(0.0, 0, 6), "0");
        assert_eq!(fmt_g(100000.0, 0, 6), "100000");
        assert_eq!(fmt_g(1000000.0, 0, 6), "1e+06");
        assert_eq!(fmt_g(0.0001234, 0, 6), "0.0001234");
        assert_eq!(fmt_g(1.5e-5, 0, 6), "1.5e-05");
        assert_eq!(fmt_g(0.4, 0, 6), "0.4");
        assert_eq!(fmt_g(123.456789, 0, 6), "123.457");
    }
}
