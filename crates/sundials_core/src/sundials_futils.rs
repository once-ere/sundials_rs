/* -----------------------------------------------------------------
 * Translated from src/sundials/sundials_futils.c and
 * include/sundials/sundials_futils.h (SUNDIALS 7.7.0).
 * SUNDIALS FILE interface utility implementations.
 *
 * The C `FILE*` becomes the owned enum `SUNFileHandle`; the C NULL
 * file pointer is the `Null` variant (so the out-param convention
 * `FILE** fp_out` maps to `&mut SUNFileHandle`). A NULL `filename`
 * maps to `None`, in which case the incoming handle is kept, exactly
 * like the C code reads `fp = *fp_out`. Closing replaces a `File`
 * handle with `Null` (RAII closes the file); the C code leaves the
 * stale pointer behind, which safe Rust cannot express.
 *
 * `SUNFileHandle` implements `std::io::Write` so downstream loggers
 * can write to it; writing to `Null` discards the bytes (in C,
 * writing through a NULL `FILE*` is undefined behavior).
 * -----------------------------------------------------------------*/
use crate::sundials_errors::*;

/// C `FILE*` handle: stdout, stderr, an owned file, or NULL.
#[derive(Debug, Default)]
pub enum SUNFileHandle {
    Stdout,
    Stderr,
    File(std::fs::File),
    #[default]
    Null,
}

/// C `fopen(filename, mode)`: mode strings "r", "w", "a", "r+",
/// "w+", "a+" with an optional 'b' (binary is irrelevant here).
/// Returns None (C NULL) on unknown mode or open failure.
fn fopen(filename: &str, mode: &str) -> Option<std::fs::File> {
    let base = mode.chars().next()?;
    let plus = mode.contains('+');
    let mut opts = std::fs::OpenOptions::new();
    match base {
        'r' => {
            opts.read(true);
            if plus {
                opts.write(true);
            }
        }
        'w' => {
            opts.write(true).create(true).truncate(true);
            if plus {
                opts.read(true);
            }
        }
        'a' => {
            opts.append(true).create(true);
            if plus {
                opts.read(true);
            }
        }
        _ => return None,
    }
    opts.open(filename).ok()
}

/// SUNFileOpen: create a file handle with the given file name and mode.
pub fn SUNFileOpen(filename: Option<&str>, mode: &str, fp_out: &mut SUNFileHandle) -> SUNErrCode {
    let mut err = SUN_SUCCESS;

    if let Some(filename) = filename {
        if filename == "stdout" {
            *fp_out = SUNFileHandle::Stdout;
        } else if filename == "stderr" {
            *fp_out = SUNFileHandle::Stderr;
        } else {
            *fp_out = match fopen(filename, mode) {
                Some(f) => SUNFileHandle::File(f),
                None => SUNFileHandle::Null,
            };
        }
    }

    if matches!(fp_out, SUNFileHandle::Null) {
        err = SUN_ERR_FILE_OPEN;
    }

    err
}

/// SUNDIALSFileOpen
pub fn SUNDIALSFileOpen(
    filename: Option<&str>,
    mode: &str,
    fp_out: &mut SUNFileHandle,
) -> SUNErrCode {
    SUNFileOpen(filename, mode, fp_out)
}

/// SUNFileClose: close a file handle (stdout/stderr are not closed).
pub fn SUNFileClose(fp_ptr: &mut SUNFileHandle) -> SUNErrCode {
    if matches!(fp_ptr, SUNFileHandle::File(_)) {
        *fp_ptr = SUNFileHandle::Null; /* fclose(fp) via RAII */
    }
    SUN_SUCCESS
}

/// SUNDIALSFileClose
pub fn SUNDIALSFileClose(fp_ptr: &mut SUNFileHandle) -> SUNErrCode {
    SUNFileClose(fp_ptr)
}

impl std::io::Write for SUNFileHandle {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            SUNFileHandle::Stdout => std::io::stdout().write(buf),
            SUNFileHandle::Stderr => std::io::stderr().write(buf),
            SUNFileHandle::File(f) => f.write(buf),
            SUNFileHandle::Null => Ok(buf.len()),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            SUNFileHandle::Stdout => std::io::stdout().flush(),
            SUNFileHandle::Stderr => std::io::stderr().flush(),
            SUNFileHandle::File(f) => f.flush(),
            SUNFileHandle::Null => Ok(()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn open_write_close_roundtrip() {
        let path = std::env::temp_dir().join("sundials_rs_futils_test.txt");
        let path_str = path.to_str().unwrap();

        let mut fp = SUNFileHandle::Null;
        assert_eq!(SUNDIALSFileOpen(Some(path_str), "w", &mut fp), SUN_SUCCESS);
        assert!(matches!(fp, SUNFileHandle::File(_)));
        write!(fp, "hello sundials").unwrap();
        fp.flush().unwrap();
        assert_eq!(SUNDIALSFileClose(&mut fp), SUN_SUCCESS);
        assert!(matches!(fp, SUNFileHandle::Null));

        assert_eq!(std::fs::read_to_string(&path).unwrap(), "hello sundials");

        /* append mode adds to the end */
        let mut fp2 = SUNFileHandle::Null;
        assert_eq!(SUNFileOpen(Some(path_str), "a", &mut fp2), SUN_SUCCESS);
        write!(fp2, "!").unwrap();
        assert_eq!(SUNFileClose(&mut fp2), SUN_SUCCESS);
        assert_eq!(std::fs::read_to_string(&path).unwrap(), "hello sundials!");

        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn special_names_and_failures() {
        let mut fp = SUNFileHandle::Null;
        assert_eq!(SUNFileOpen(Some("stdout"), "w", &mut fp), SUN_SUCCESS);
        assert!(matches!(fp, SUNFileHandle::Stdout));
        /* closing stdout is a no-op */
        assert_eq!(SUNFileClose(&mut fp), SUN_SUCCESS);
        assert!(matches!(fp, SUNFileHandle::Stdout));

        let mut fp = SUNFileHandle::Null;
        assert_eq!(SUNFileOpen(Some("stderr"), "w", &mut fp), SUN_SUCCESS);
        assert!(matches!(fp, SUNFileHandle::Stderr));

        /* opening a nonexistent file for reading fails */
        let missing = std::env::temp_dir().join("sundials_rs_futils_missing.txt");
        std::fs::remove_file(&missing).ok();
        let mut fp = SUNFileHandle::Null;
        assert_eq!(
            SUNFileOpen(Some(missing.to_str().unwrap()), "r", &mut fp),
            SUN_ERR_FILE_OPEN
        );
        assert!(matches!(fp, SUNFileHandle::Null));

        /* NULL filename keeps the incoming handle */
        let mut fp = SUNFileHandle::Stderr;
        assert_eq!(SUNFileOpen(None, "w", &mut fp), SUN_SUCCESS);
        assert!(matches!(fp, SUNFileHandle::Stderr));
        let mut fp = SUNFileHandle::Null;
        assert_eq!(SUNFileOpen(None, "w", &mut fp), SUN_ERR_FILE_OPEN);

        /* writing to a Null handle discards the bytes */
        let mut fp = SUNFileHandle::Null;
        assert_eq!(fp.write(b"discarded").unwrap(), 9);
    }
}
