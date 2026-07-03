/* -----------------------------------------------------------------
 * Translated from src/sundials/sundials_cli.c and
 * src/sundials/sundials_cli.h (SUNDIALS 7.7.0).
 *
 * Command-line input utility routines: each helper scans a dispatch
 * table of (key, set-routine) pairs against argv[*argidx] (with the
 * solver-id prefix of length `offset` stripped), consumes the value
 * argument(s) that follow, and invokes the matching set routine.
 *
 * Rust adaptations (behavior otherwise follows the C line-for-line):
 *  - The C `void* mem` handle becomes a generic parameter `M`; the
 *    sun*SetFn typedefs are plain `fn` pointers taking `&mut M`
 *    (workspace callback convention). Note cvode_rs/src/cvode_cli.rs
 *    predates this module and carries its own CVodeMem-monomorphic
 *    copies of five of these helpers; it is left untouched.
 *  - `char* argv[]` becomes `&[String]`, `int* argidx` and
 *    `int* failedarg` become `&mut usize` (Rust slice indices).
 *    The C `int numpairs` parameter is dropped: the `testpairs`
 *    slice carries its own length.
 *  - atoi/atol are reproduced by sun_atoi/sun_atol (leading
 *    optionally-signed integer, 0 on failure); SUNStrToReal comes
 *    from sundials_math.rs (strtod equivalent, 0.0 on failure).
 *  - Where C indexes past the argv array (missing value arguments,
 *    or `argv[*argidx] + offset` past the end of the string — both
 *    undefined behavior in C, guarded by the callers), the Rust port
 *    substitutes "" so the parse yields 0 / 0.0 / no key match
 *    instead of faulting.
 * -----------------------------------------------------------------*/

use crate::sundials_errors::{SUNErrCode, SUN_SUCCESS};
use crate::sundials_math::SUNStrToReal;
use crate::sundials_types::{sunbooleantype, sunrealtype, SUNFALSE, SUNTRUE};

/* ----------------------------------------------------------------
 * local parsing helpers (C atoi/atol semantics)
 * ---------------------------------------------------------------- */

/// atol: parse a leading (optionally signed) integer, 0 on failure.
fn sun_atol(s: &str) -> i64 {
    let t = s.trim_start();
    let bytes = t.as_bytes();
    let mut i = 0;
    if i < bytes.len() && (bytes[i] == b'+' || bytes[i] == b'-') {
        i += 1;
    }
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        i += 1;
    }
    t[..i].parse().unwrap_or(0)
}

/// atoi
fn sun_atoi(s: &str) -> i32 {
    sun_atol(s) as i32
}

/// argv[idx] with C's out-of-bounds access replaced by "".
fn arg_at(argv: &[String], idx: usize) -> &str {
    argv.get(idx).map(String::as_str).unwrap_or("")
}

/// argv[*argidx] + offset with C's past-the-end pointer replaced by "".
fn arg_key(argv: &[String], idx: usize, offset: usize) -> &str {
    arg_at(argv, idx).get(offset..).unwrap_or("")
}

/*===============================================================
  Command-line input utility routines
  ===============================================================*/

/* utilities for integer "set" routines */
pub type sunIntSetFn<M> = fn(&mut M, i32) -> SUNErrCode;

pub struct sunKeyIntPair<M> {
    pub key: &'static str,
    pub set: sunIntSetFn<M>,
}

/// sunCheckAndSetIntArgs
pub fn sunCheckAndSetIntArgs<M>(
    mem: &mut M,
    argidx: &mut usize,
    argv: &[String],
    offset: usize,
    testpairs: &[sunKeyIntPair<M>],
    arg_used: &mut sunbooleantype,
    failedarg: &mut usize,
) -> SUNErrCode {
    for (j, pair) in testpairs.iter().enumerate() {
        *arg_used = SUNFALSE;
        if arg_key(argv, *argidx, offset) == pair.key {
            *argidx += 1;
            let iarg = sun_atoi(arg_at(argv, *argidx));
            let retval = (pair.set)(mem, iarg);
            if retval != SUN_SUCCESS {
                *failedarg = j;
                return retval;
            }
            *arg_used = SUNTRUE;
            return SUN_SUCCESS;
        }
    }
    SUN_SUCCESS
}

/* utilities for pair-of-integer "set" routines */
pub type sunTwoIntSetFn<M> = fn(&mut M, i32, i32) -> SUNErrCode;

pub struct sunKeyTwoIntPair<M> {
    pub key: &'static str,
    pub set: sunTwoIntSetFn<M>,
}

/// sunCheckAndSetTwoIntArgs
pub fn sunCheckAndSetTwoIntArgs<M>(
    mem: &mut M,
    argidx: &mut usize,
    argv: &[String],
    offset: usize,
    testpairs: &[sunKeyTwoIntPair<M>],
    arg_used: &mut sunbooleantype,
    failedarg: &mut usize,
) -> SUNErrCode {
    for (j, pair) in testpairs.iter().enumerate() {
        *arg_used = SUNFALSE;
        if arg_key(argv, *argidx, offset) == pair.key {
            *argidx += 1;
            let iarg1 = sun_atoi(arg_at(argv, *argidx));
            *argidx += 1;
            let iarg2 = sun_atoi(arg_at(argv, *argidx));
            let retval = (pair.set)(mem, iarg1, iarg2);
            if retval != SUN_SUCCESS {
                *failedarg = j;
                return retval;
            }
            *arg_used = SUNTRUE;
            return SUN_SUCCESS;
        }
    }
    SUN_SUCCESS
}

/* utilities for long int "set" routines */
pub type sunLongSetFn<M> = fn(&mut M, i64) -> SUNErrCode;

pub struct sunKeyLongPair<M> {
    pub key: &'static str,
    pub set: sunLongSetFn<M>,
}

/// sunCheckAndSetLongArgs
pub fn sunCheckAndSetLongArgs<M>(
    mem: &mut M,
    argidx: &mut usize,
    argv: &[String],
    offset: usize,
    testpairs: &[sunKeyLongPair<M>],
    arg_used: &mut sunbooleantype,
    failedarg: &mut usize,
) -> SUNErrCode {
    for (j, pair) in testpairs.iter().enumerate() {
        *arg_used = SUNFALSE;
        if arg_key(argv, *argidx, offset) == pair.key {
            *argidx += 1;
            let iarg = sun_atol(arg_at(argv, *argidx));
            let retval = (pair.set)(mem, iarg);
            if retval != SUN_SUCCESS {
                *failedarg = j;
                return retval;
            }
            *arg_used = SUNTRUE;
            return SUN_SUCCESS;
        }
    }
    SUN_SUCCESS
}

/* utilities for pair int/sunrealtype "set" routines */
pub type sunIntRealSetFn<M> = fn(&mut M, i32, sunrealtype) -> SUNErrCode;

pub struct sunKeyIntRealPair<M> {
    pub key: &'static str,
    pub set: sunIntRealSetFn<M>,
}

/// sunCheckAndSetIntRealArgs
pub fn sunCheckAndSetIntRealArgs<M>(
    mem: &mut M,
    argidx: &mut usize,
    argv: &[String],
    offset: usize,
    testpairs: &[sunKeyIntRealPair<M>],
    arg_used: &mut sunbooleantype,
    failedarg: &mut usize,
) -> SUNErrCode {
    for (j, pair) in testpairs.iter().enumerate() {
        *arg_used = SUNFALSE;
        if arg_key(argv, *argidx, offset) == pair.key {
            *argidx += 1;
            let iarg = sun_atoi(arg_at(argv, *argidx));
            *argidx += 1;
            let rarg = SUNStrToReal(arg_at(argv, *argidx));
            let retval = (pair.set)(mem, iarg, rarg);
            if retval != SUN_SUCCESS {
                *failedarg = j;
                return retval;
            }
            *arg_used = SUNTRUE;
            return SUN_SUCCESS;
        }
    }
    SUN_SUCCESS
}

/* utilities for triplet int/sunrealtype/sunrealtype "set" routines */
pub type sunIntRealRealSetFn<M> = fn(&mut M, i32, sunrealtype, sunrealtype) -> SUNErrCode;

pub struct sunKeyIntRealRealPair<M> {
    pub key: &'static str,
    pub set: sunIntRealRealSetFn<M>,
}

/// sunCheckAndSetIntRealRealArgs
pub fn sunCheckAndSetIntRealRealArgs<M>(
    mem: &mut M,
    argidx: &mut usize,
    argv: &[String],
    offset: usize,
    testpairs: &[sunKeyIntRealRealPair<M>],
    arg_used: &mut sunbooleantype,
    failedarg: &mut usize,
) -> SUNErrCode {
    for (j, pair) in testpairs.iter().enumerate() {
        *arg_used = SUNFALSE;
        if arg_key(argv, *argidx, offset) == pair.key {
            *argidx += 1;
            let iarg = sun_atoi(arg_at(argv, *argidx));
            *argidx += 1;
            let rarg1 = SUNStrToReal(arg_at(argv, *argidx));
            *argidx += 1;
            let rarg2 = SUNStrToReal(arg_at(argv, *argidx));
            let retval = (pair.set)(mem, iarg, rarg1, rarg2);
            if retval != SUN_SUCCESS {
                *failedarg = j;
                return retval;
            }
            *arg_used = SUNTRUE;
            return SUN_SUCCESS;
        }
    }
    SUN_SUCCESS
}

/* utilities for pair int/long int "set" routines */
pub type sunIntLongSetFn<M> = fn(&mut M, i32, i64) -> SUNErrCode;

pub struct sunKeyIntLongPair<M> {
    pub key: &'static str,
    pub set: sunIntLongSetFn<M>,
}

/// sunCheckAndSetIntLongArgs
pub fn sunCheckAndSetIntLongArgs<M>(
    mem: &mut M,
    argidx: &mut usize,
    argv: &[String],
    offset: usize,
    testpairs: &[sunKeyIntLongPair<M>],
    arg_used: &mut sunbooleantype,
    failedarg: &mut usize,
) -> SUNErrCode {
    for (j, pair) in testpairs.iter().enumerate() {
        *arg_used = SUNFALSE;
        if arg_key(argv, *argidx, offset) == pair.key {
            *argidx += 1;
            let iarg = sun_atoi(arg_at(argv, *argidx));
            *argidx += 1;
            let large = sun_atol(arg_at(argv, *argidx));
            let retval = (pair.set)(mem, iarg, large);
            if retval != SUN_SUCCESS {
                *failedarg = j;
                return retval;
            }
            *arg_used = SUNTRUE;
            return SUN_SUCCESS;
        }
    }
    SUN_SUCCESS
}

/* utilities for sunrealtype "set" routines */
pub type sunRealSetFn<M> = fn(&mut M, sunrealtype) -> SUNErrCode;

pub struct sunKeyRealPair<M> {
    pub key: &'static str,
    pub set: sunRealSetFn<M>,
}

/// sunCheckAndSetRealArgs
pub fn sunCheckAndSetRealArgs<M>(
    mem: &mut M,
    argidx: &mut usize,
    argv: &[String],
    offset: usize,
    testpairs: &[sunKeyRealPair<M>],
    arg_used: &mut sunbooleantype,
    failedarg: &mut usize,
) -> SUNErrCode {
    for (j, pair) in testpairs.iter().enumerate() {
        *arg_used = SUNFALSE;
        if arg_key(argv, *argidx, offset) == pair.key {
            *argidx += 1;
            let rarg = SUNStrToReal(arg_at(argv, *argidx));
            let retval = (pair.set)(mem, rarg);
            if retval != SUN_SUCCESS {
                *failedarg = j;
                return retval;
            }
            *arg_used = SUNTRUE;
            return SUN_SUCCESS;
        }
    }
    SUN_SUCCESS
}

/* utilities for pair-of-sunrealtype "set" routines */
pub type sunTwoRealSetFn<M> = fn(&mut M, sunrealtype, sunrealtype) -> SUNErrCode;

pub struct sunKeyTwoRealPair<M> {
    pub key: &'static str,
    pub set: sunTwoRealSetFn<M>,
}

/// sunCheckAndSetTwoRealArgs
pub fn sunCheckAndSetTwoRealArgs<M>(
    mem: &mut M,
    argidx: &mut usize,
    argv: &[String],
    offset: usize,
    testpairs: &[sunKeyTwoRealPair<M>],
    arg_used: &mut sunbooleantype,
    failedarg: &mut usize,
) -> SUNErrCode {
    for (j, pair) in testpairs.iter().enumerate() {
        *arg_used = SUNFALSE;
        if arg_key(argv, *argidx, offset) == pair.key {
            *argidx += 1;
            let rarg1 = SUNStrToReal(arg_at(argv, *argidx));
            *argidx += 1;
            let rarg2 = SUNStrToReal(arg_at(argv, *argidx));
            let retval = (pair.set)(mem, rarg1, rarg2);
            if retval != SUN_SUCCESS {
                *failedarg = j;
                return retval;
            }
            *arg_used = SUNTRUE;
            return SUN_SUCCESS;
        }
    }
    SUN_SUCCESS
}

/* utilities for char* "set" routines */
pub type sunCharSetFn<M> = fn(&mut M, &str) -> SUNErrCode;

pub struct sunKeyCharPair<M> {
    pub key: &'static str,
    pub set: sunCharSetFn<M>,
}

/// sunCheckAndSetCharArgs
pub fn sunCheckAndSetCharArgs<M>(
    mem: &mut M,
    argidx: &mut usize,
    argv: &[String],
    offset: usize,
    testpairs: &[sunKeyCharPair<M>],
    arg_used: &mut sunbooleantype,
    failedarg: &mut usize,
) -> SUNErrCode {
    for (j, pair) in testpairs.iter().enumerate() {
        *arg_used = SUNFALSE;
        if arg_key(argv, *argidx, offset) == pair.key {
            *argidx += 1;
            let retval = (pair.set)(mem, arg_at(argv, *argidx));
            if retval != SUN_SUCCESS {
                *failedarg = j;
                return retval;
            }
            *arg_used = SUNTRUE;
            return SUN_SUCCESS;
        }
    }
    SUN_SUCCESS
}

/* utilities for pair-of-char* "set" routines */
pub type sunTwoCharSetFn<M> = fn(&mut M, &str, &str) -> SUNErrCode;

pub struct sunKeyTwoCharPair<M> {
    pub key: &'static str,
    pub set: sunTwoCharSetFn<M>,
}

/// sunCheckAndSetTwoCharArgs — note the C original invokes the set
/// routine on argv[*argidx + 1] / argv[*argidx + 2] and only then
/// advances *argidx by 2; that order is preserved.
pub fn sunCheckAndSetTwoCharArgs<M>(
    mem: &mut M,
    argidx: &mut usize,
    argv: &[String],
    offset: usize,
    testpairs: &[sunKeyTwoCharPair<M>],
    arg_used: &mut sunbooleantype,
    failedarg: &mut usize,
) -> SUNErrCode {
    for (j, pair) in testpairs.iter().enumerate() {
        *arg_used = SUNFALSE;
        if arg_key(argv, *argidx, offset) == pair.key {
            let retval = (pair.set)(mem, arg_at(argv, *argidx + 1), arg_at(argv, *argidx + 2));
            *argidx += 2;
            if retval != SUN_SUCCESS {
                *failedarg = j;
                return retval;
            }
            *arg_used = SUNTRUE;
            return SUN_SUCCESS;
        }
    }
    SUN_SUCCESS
}

/* utilities for action "set" routines */
pub type sunActionSetFn<M> = fn(&mut M) -> SUNErrCode;

pub struct sunKeyActionPair<M> {
    pub key: &'static str,
    pub set: sunActionSetFn<M>,
}

/// sunCheckAndSetActionArgs
pub fn sunCheckAndSetActionArgs<M>(
    mem: &mut M,
    argidx: &mut usize,
    argv: &[String],
    offset: usize,
    testpairs: &[sunKeyActionPair<M>],
    arg_used: &mut sunbooleantype,
    failedarg: &mut usize,
) -> SUNErrCode {
    for (j, pair) in testpairs.iter().enumerate() {
        *arg_used = SUNFALSE;
        if arg_key(argv, *argidx, offset) == pair.key {
            let retval = (pair.set)(mem);
            if retval != SUN_SUCCESS {
                *failedarg = j;
                return retval;
            }
            *arg_used = SUNTRUE;
            return SUN_SUCCESS;
        }
    }
    SUN_SUCCESS
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Default)]
    struct TestMem {
        i: i32,
        two_i: (i32, i32),
        l: i64,
        r: f64,
        two_r: (f64, f64),
        int_real: (i32, f64),
        int_real_real: (i32, f64, f64),
        int_long: (i32, i64),
        c: String,
        two_c: (String, String),
        action_hits: i32,
    }

    fn argv(items: &[&str]) -> Vec<String> {
        items.iter().map(|s| s.to_string()).collect()
    }

    const OFFSET: usize = 5; /* strlen("test.") */

    #[test]
    fn int_args_accept_matching_key_and_consume_value() {
        let mut mem = TestMem::default();
        let pairs = [
            sunKeyIntPair::<TestMem> { key: "other", set: |_, _| SUN_SUCCESS },
            sunKeyIntPair::<TestMem> { key: "max_order", set: |m, v| { m.i = v; SUN_SUCCESS } },
        ];
        let args = argv(&["test.max_order", "5"]);
        let mut idx = 0usize;
        let mut used = SUNFALSE;
        let mut failed = 0usize;
        let ret = sunCheckAndSetIntArgs(&mut mem, &mut idx, &args, OFFSET, &pairs, &mut used, &mut failed);
        assert_eq!(ret, SUN_SUCCESS);
        assert!(used);
        assert_eq!(idx, 1); /* advanced onto the consumed value */
        assert_eq!(mem.i, 5);
    }

    #[test]
    fn unmatched_key_leaves_index_and_reports_unused() {
        let mut mem = TestMem::default();
        let pairs = [sunKeyIntPair::<TestMem> { key: "max_order", set: |m, v| { m.i = v; SUN_SUCCESS } }];
        let args = argv(&["test.unknown_key", "5"]);
        let mut idx = 0usize;
        let mut used = SUNTRUE;
        let mut failed = 0usize;
        let ret = sunCheckAndSetIntArgs(&mut mem, &mut idx, &args, OFFSET, &pairs, &mut used, &mut failed);
        assert_eq!(ret, SUN_SUCCESS);
        assert!(!used);
        assert_eq!(idx, 0);
        assert_eq!(mem.i, 0);
    }

    #[test]
    fn failing_set_routine_reports_failed_index() {
        let mut mem = TestMem::default();
        let pairs = [
            sunKeyIntPair::<TestMem> { key: "ok", set: |_, _| SUN_SUCCESS },
            sunKeyIntPair::<TestMem> { key: "bad", set: |_, _| -99 },
        ];
        let args = argv(&["test.bad", "7"]);
        let mut idx = 0usize;
        let mut used = SUNTRUE;
        let mut failed = 0usize;
        let ret = sunCheckAndSetIntArgs(&mut mem, &mut idx, &args, OFFSET, &pairs, &mut used, &mut failed);
        assert_eq!(ret, -99);
        assert_eq!(failed, 1);
        assert!(!used); /* C leaves arg_used SUNFALSE on failure */
    }

    #[test]
    fn non_numeric_value_parses_as_zero_like_atoi() {
        let mut mem = TestMem::default();
        mem.i = 42;
        let pairs = [sunKeyIntPair::<TestMem> { key: "n", set: |m, v| { m.i = v; SUN_SUCCESS } }];
        let args = argv(&["test.n", "not_a_number"]);
        let mut idx = 0usize;
        let mut used = SUNFALSE;
        let mut failed = 0usize;
        assert_eq!(
            sunCheckAndSetIntArgs(&mut mem, &mut idx, &args, OFFSET, &pairs, &mut used, &mut failed),
            SUN_SUCCESS
        );
        assert!(used);
        assert_eq!(mem.i, 0); /* atoi("not_a_number") == 0 */
    }

    #[test]
    fn two_int_and_int_long_and_long_args() {
        let mut mem = TestMem::default();
        let two_pairs = [sunKeyTwoIntPair::<TestMem> {
            key: "pair",
            set: |m, a, b| { m.two_i = (a, b); SUN_SUCCESS },
        }];
        let args = argv(&["test.pair", "3", "-4"]);
        let (mut idx, mut used, mut failed) = (0usize, SUNFALSE, 0usize);
        assert_eq!(
            sunCheckAndSetTwoIntArgs(&mut mem, &mut idx, &args, OFFSET, &two_pairs, &mut used, &mut failed),
            SUN_SUCCESS
        );
        assert!(used);
        assert_eq!(idx, 2);
        assert_eq!(mem.two_i, (3, -4));

        let long_pairs = [sunKeyLongPair::<TestMem> {
            key: "steps",
            set: |m, v| { m.l = v; SUN_SUCCESS },
        }];
        let args = argv(&["test.steps", "123456789012"]);
        let (mut idx, mut used, mut failed) = (0usize, SUNFALSE, 0usize);
        assert_eq!(
            sunCheckAndSetLongArgs(&mut mem, &mut idx, &args, OFFSET, &long_pairs, &mut used, &mut failed),
            SUN_SUCCESS
        );
        assert!(used);
        assert_eq!(mem.l, 123456789012);

        let il_pairs = [sunKeyIntLongPair::<TestMem> {
            key: "il",
            set: |m, a, b| { m.int_long = (a, b); SUN_SUCCESS },
        }];
        let args = argv(&["test.il", "2", "-30"]);
        let (mut idx, mut used, mut failed) = (0usize, SUNFALSE, 0usize);
        assert_eq!(
            sunCheckAndSetIntLongArgs(&mut mem, &mut idx, &args, OFFSET, &il_pairs, &mut used, &mut failed),
            SUN_SUCCESS
        );
        assert!(used);
        assert_eq!(idx, 2);
        assert_eq!(mem.int_long, (2, -30));
    }

    #[test]
    fn real_flavors_parse_with_strtod_semantics() {
        let mut mem = TestMem::default();
        let r_pairs = [sunKeyRealPair::<TestMem> { key: "tol", set: |m, v| { m.r = v; SUN_SUCCESS } }];
        let args = argv(&["test.tol", "1.5e-2"]);
        let (mut idx, mut used, mut failed) = (0usize, SUNFALSE, 0usize);
        assert_eq!(
            sunCheckAndSetRealArgs(&mut mem, &mut idx, &args, OFFSET, &r_pairs, &mut used, &mut failed),
            SUN_SUCCESS
        );
        assert!(used);
        assert_eq!(mem.r, 1.5e-2);

        let rr_pairs = [sunKeyTwoRealPair::<TestMem> {
            key: "tols",
            set: |m, a, b| { m.two_r = (a, b); SUN_SUCCESS },
        }];
        let args = argv(&["test.tols", "1e-4", "1e-9"]);
        let (mut idx, mut used, mut failed) = (0usize, SUNFALSE, 0usize);
        assert_eq!(
            sunCheckAndSetTwoRealArgs(&mut mem, &mut idx, &args, OFFSET, &rr_pairs, &mut used, &mut failed),
            SUN_SUCCESS
        );
        assert!(used);
        assert_eq!(idx, 2);
        assert_eq!(mem.two_r, (1e-4, 1e-9));

        let ir_pairs = [sunKeyIntRealPair::<TestMem> {
            key: "ir",
            set: |m, a, b| { m.int_real = (a, b); SUN_SUCCESS },
        }];
        let args = argv(&["test.ir", "3", "0.5"]);
        let (mut idx, mut used, mut failed) = (0usize, SUNFALSE, 0usize);
        assert_eq!(
            sunCheckAndSetIntRealArgs(&mut mem, &mut idx, &args, OFFSET, &ir_pairs, &mut used, &mut failed),
            SUN_SUCCESS
        );
        assert!(used);
        assert_eq!(mem.int_real, (3, 0.5));

        let irr_pairs = [sunKeyIntRealRealPair::<TestMem> {
            key: "irr",
            set: |m, a, b, c| { m.int_real_real = (a, b, c); SUN_SUCCESS },
        }];
        let args = argv(&["test.irr", "1", "2.5", "-3.5"]);
        let (mut idx, mut used, mut failed) = (0usize, SUNFALSE, 0usize);
        assert_eq!(
            sunCheckAndSetIntRealRealArgs(&mut mem, &mut idx, &args, OFFSET, &irr_pairs, &mut used, &mut failed),
            SUN_SUCCESS
        );
        assert!(used);
        assert_eq!(idx, 3);
        assert_eq!(mem.int_real_real, (1, 2.5, -3.5));
    }

    #[test]
    fn char_two_char_and_action_args() {
        let mut mem = TestMem::default();
        let c_pairs = [sunKeyCharPair::<TestMem> {
            key: "name",
            set: |m, s| { m.c = s.to_string(); SUN_SUCCESS },
        }];
        let args = argv(&["test.name", "hello"]);
        let (mut idx, mut used, mut failed) = (0usize, SUNFALSE, 0usize);
        assert_eq!(
            sunCheckAndSetCharArgs(&mut mem, &mut idx, &args, OFFSET, &c_pairs, &mut used, &mut failed),
            SUN_SUCCESS
        );
        assert!(used);
        assert_eq!(mem.c, "hello");

        let cc_pairs = [sunKeyTwoCharPair::<TestMem> {
            key: "files",
            set: |m, a, b| { m.two_c = (a.to_string(), b.to_string()); SUN_SUCCESS },
        }];
        let args = argv(&["test.files", "in.txt", "out.txt"]);
        let (mut idx, mut used, mut failed) = (0usize, SUNFALSE, 0usize);
        assert_eq!(
            sunCheckAndSetTwoCharArgs(&mut mem, &mut idx, &args, OFFSET, &cc_pairs, &mut used, &mut failed),
            SUN_SUCCESS
        );
        assert!(used);
        assert_eq!(idx, 2);
        assert_eq!(mem.two_c, ("in.txt".to_string(), "out.txt".to_string()));

        let a_pairs = [sunKeyActionPair::<TestMem> {
            key: "clear",
            set: |m| { m.action_hits += 1; SUN_SUCCESS },
        }];
        let args = argv(&["test.clear"]);
        let (mut idx, mut used, mut failed) = (0usize, SUNFALSE, 0usize);
        assert_eq!(
            sunCheckAndSetActionArgs(&mut mem, &mut idx, &args, OFFSET, &a_pairs, &mut used, &mut failed),
            SUN_SUCCESS
        );
        assert!(used);
        assert_eq!(idx, 0); /* action args consume no value */
        assert_eq!(mem.action_hits, 1);
    }

    #[test]
    fn missing_value_arguments_parse_as_empty() {
        /* C would index past argv (UB); the port substitutes "" -> 0 */
        let mut mem = TestMem::default();
        mem.i = 7;
        let pairs = [sunKeyIntPair::<TestMem> { key: "n", set: |m, v| { m.i = v; SUN_SUCCESS } }];
        let args = argv(&["test.n"]);
        let (mut idx, mut used, mut failed) = (0usize, SUNFALSE, 0usize);
        assert_eq!(
            sunCheckAndSetIntArgs(&mut mem, &mut idx, &args, OFFSET, &pairs, &mut used, &mut failed),
            SUN_SUCCESS
        );
        assert!(used);
        assert_eq!(mem.i, 0);
    }
}
