/* -----------------------------------------------------------------
 * Translated from src/sundials/sundials_version.c (SUNDIALS 7.7.0).
 * This file implements functions for getting SUNDIALS version
 * information.
 *
 * The SUNDIALS_VERSION* constants come from
 * include/sundials/sundials_config.in configured for release 7.7.0
 * (version label empty). C `char* version, int len` buffers become
 * `&mut String` plus the C `len` argument; the C NULL-pointer checks
 * (SUN_ERR_ARG_CORRUPT) vanish because references cannot be null,
 * while the length checks are kept verbatim (including the C
 * `(size_t)len` cast, so a negative `len` wraps to a huge value and
 * passes the check exactly as in C).
 * -----------------------------------------------------------------*/
use crate::sundials_errors::*;

/* sundials_config.h for release 7.7.0 */
pub const SUNDIALS_VERSION: &str = "7.7.0";
pub const SUNDIALS_VERSION_MAJOR: i32 = 7;
pub const SUNDIALS_VERSION_MINOR: i32 = 7;
pub const SUNDIALS_VERSION_PATCH: i32 = 0;
pub const SUNDIALS_VERSION_LABEL: &str = "";

/* note strlen does not include terminating null character hence the
   use of >= when checking len below and strncpy copies up to len
   characters including the terminating null character */

/// SUNDIALSGetVersion: fill string with SUNDIALS version information.
pub fn SUNDIALSGetVersion(version: &mut String, len: i32) -> SUNErrCode {
    if SUNDIALS_VERSION.len() >= len as usize {
        return SUN_ERR_ARG_OUTOFRANGE;
    }

    version.clear();
    version.push_str(SUNDIALS_VERSION);

    SUN_SUCCESS
}

/// SUNDIALSGetVersionNumber: fill integers with SUNDIALS major, minor,
/// and patch release numbers and fill a string with the release label.
pub fn SUNDIALSGetVersionNumber(
    major: &mut i32,
    minor: &mut i32,
    patch: &mut i32,
    label: &mut String,
    len: i32,
) -> SUNErrCode {
    if SUNDIALS_VERSION_LABEL.len() >= len as usize {
        return SUN_ERR_ARG_OUTOFRANGE;
    }

    *major = SUNDIALS_VERSION_MAJOR;
    *minor = SUNDIALS_VERSION_MINOR;
    *patch = SUNDIALS_VERSION_PATCH;
    label.clear();
    label.push_str(SUNDIALS_VERSION_LABEL);

    SUN_SUCCESS
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn get_version() {
        let mut version = String::new();
        assert_eq!(SUNDIALSGetVersion(&mut version, 25), SUN_SUCCESS);
        assert_eq!(version, "7.7.0");

        /* strlen("7.7.0") == 5, so len == 5 is too short (needs the
        terminating null in C), len == 6 is the minimum */
        let mut short = String::new();
        assert_eq!(SUNDIALSGetVersion(&mut short, 5), SUN_ERR_ARG_OUTOFRANGE);
        assert_eq!(SUNDIALSGetVersion(&mut short, 6), SUN_SUCCESS);
        assert_eq!(short, "7.7.0");
    }

    #[test]
    fn get_version_number() {
        let (mut major, mut minor, mut patch) = (-1, -1, -1);
        let mut label = String::from("junk");
        assert_eq!(
            SUNDIALSGetVersionNumber(&mut major, &mut minor, &mut patch, &mut label, 10),
            SUN_SUCCESS
        );
        assert_eq!((major, minor, patch), (7, 7, 0));
        assert_eq!(label, "");

        /* empty label still needs room for the terminating null */
        assert_eq!(
            SUNDIALSGetVersionNumber(&mut major, &mut minor, &mut patch, &mut label, 0),
            SUN_ERR_ARG_OUTOFRANGE
        );
    }
}
