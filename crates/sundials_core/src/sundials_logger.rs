/* -----------------------------------------------------------------
 * Translated from src/sundials/sundials_logger.c and
 * src/sundials/sundials_logger_impl.h (SUNDIALS 7.7.0).
 *
 * SUNLogger: leveled (error/warning/info/debug) logging to per-level
 * output streams. Filenames are deduplicated through a SUNHashMap so
 * the same file used for multiple levels is only opened once.
 *
 * Adaptations for the serial pure-Rust build:
 *  - MPI is excluded: the `SUNComm comm` argument of SUNLogger_Create
 *    / SUNLogger_CreateFromEnv is dropped (same policy as the
 *    SUNContext port) and `sunLoggerIsOutputRank` always reports rank
 *    0 / SUNTRUE, exactly like the C `#else` (non-MPI) branch.
 *  - The C sources compile most bodies only when the cpp symbol
 *    SUNDIALS_LOGGING_LEVEL enables them. This port is always fully
 *    enabled (equivalent to a C build with SUNDIALS_LOGGING_LEVEL=4);
 *    the solver translations keep their log call sites compiled out,
 *    so this module is complete and correct standalone.
 *  - `SUNLogger_QueueMsg` takes a pre-formatted `&str` message instead
 *    of a printf format string + varargs (workspace-wide ABI
 *    adaptation), so `sunCreateLogMessage` loses its `va_list` and
 *    returns the composed `String` instead of filling a `char**`.
 *  - C `FILE*` streams become `Box<dyn std::io::Write>` sinks created
 *    from filename strings ("stdout" -> stdout, "stderr" -> stderr,
 *    anything else -> that file truncated, like fopen mode "w+").
 *    These live in the `filenames` hashmap, which owns and closes
 *    them (drop == fclose). Because safe Rust cannot alias the map's
 *    `FILE*` into the per-level fields the way C does, the fields
 *    hold a `SunLogStream` tag that is resolved against the map at
 *    write time — observable behavior is identical.
 *  - `SUNLogger_Set*File` accepts `Option<Box<dyn Write>>` in place of
 *    a nullable `FILE*`; the logger takes ownership of the sink.
 *  - The overridable operations keep their C names as plain `fn`
 *    pointer fields; `destroy` receives `&mut SUNLogger` (the handle
 *    itself is dropped by SUNLogger_Destroy afterwards).
 * -----------------------------------------------------------------*/

use std::io::Write;

use crate::sundials_errors::*;
use crate::sundials_hashmap::*;

/* sundials_logger_impl.h */
pub const SUNDIALS_LOGGING_ERROR: i32 = 1;
pub const SUNDIALS_LOGGING_WARNING: i32 = 2;
pub const SUNDIALS_LOGGING_INFO: i32 = 3;
pub const SUNDIALS_LOGGING_DEBUG: i32 = 4;

/// enum SUNLogLevel (include/sundials/sundials_logger.h)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum SUNLogLevel {
    SUN_LOGLEVEL_ALL = -1,
    SUN_LOGLEVEL_NONE = 0,
    SUN_LOGLEVEL_ERROR = 1,
    SUN_LOGLEVEL_WARNING = 2,
    SUN_LOGLEVEL_INFO = 3,
    SUN_LOGLEVEL_DEBUG = 4,
}
pub use SUNLogLevel::*;

/* default number of files that we allocate space for */
const SUN_DEFAULT_LOGFILE_HANDLES_: i64 = 8;

/// A per-level output stream slot (C: the nullable `FILE*` fields
/// `error_fp` / `warning_fp` / `info_fp` / `debug_fp`).
pub enum SunLogStream {
    /// C NULL — output disabled for this level.
    None,
    /// The process stderr stream (default for the error level).
    Stderr,
    /// The process stdout stream (default for the warning level).
    Stdout,
    /// A filename key into the logger's `filenames` hashmap.
    Key(String),
    /// A caller-supplied sink installed via `SUNLogger_Set*File`.
    Direct(Box<dyn Write>),
}

/// struct SUNLogger_ (sundials_logger_impl.h), minus the MPI `comm`.
pub struct SUNLogger {
    pub output_rank: i32,

    /* Output streams (C: FILE* fields) */
    debug_fp: SunLogStream,
    warning_fp: SunLogStream,
    info_fp: SunLogStream,
    error_fp: SunLogStream,

    /* Hashmap used to store filename, stream pairs */
    pub filenames: Option<SUNHashMap<Box<dyn Write>>>,

    /* Slic-style format string (kept for struct parity; unused, as in C) */
    pub format: Option<String>,

    /* Content for custom implementations */
    pub content: Option<Box<dyn std::any::Any>>,

    /* Overridable operations */
    pub queuemsg:
        Option<fn(&mut SUNLogger, SUNLogLevel, &str, &str, &str) -> SUNErrCode>,
    pub flush: Option<fn(&mut SUNLogger, SUNLogLevel) -> SUNErrCode>,
    pub destroy: Option<fn(&mut SUNLogger) -> SUNErrCode>,
}

/// sunCreateLogMessage — composes
/// `[<LEVEL>][rank <rank>][<scope>][<label>] <txt>\n`.
/// The C version formatted `txt` with `va_list` args into an allocated
/// `char**` out-parameter; here `txt` arrives pre-formatted and the
/// message is returned (the "message size too large" vsnprintf error
/// path drops out).
pub fn sunCreateLogMessage(
    lvl: SUNLogLevel,
    rank: i32,
    scope: &str,
    label: &str,
    txt: &str,
) -> String {
    let prefix = if lvl == SUN_LOGLEVEL_DEBUG {
        "DEBUG"
    } else if lvl == SUN_LOGLEVEL_WARNING {
        "WARNING"
    } else if lvl == SUN_LOGLEVEL_INFO {
        "INFO"
    } else if lvl == SUN_LOGLEVEL_ERROR {
        "ERROR"
    } else {
        /* C leaves the prefix NULL here (never reached by QueueMsg's
        level switch, which rejects ALL/NONE) */
        ""
    };

    format!("[{}][rank {}][{}][{}] {}\n", prefix, rank, scope, label, txt)
}

/// sunOpenLogFile — "stdout"/"stderr" map to the process streams, any
/// other name is opened with the effect of fopen(fname, "w+")
/// (create/truncate). Returns None if the file cannot be opened.
fn sunOpenLogFile(fname: &str) -> Option<Box<dyn Write>> {
    if fname.is_empty() {
        return None;
    }
    if fname == "stdout" {
        Some(Box::new(std::io::stdout()))
    } else if fname == "stderr" {
        Some(Box::new(std::io::stderr()))
    } else {
        match std::fs::File::create(fname) {
            Ok(f) => Some(Box::new(f)),
            Err(_) => None,
        }
    }
}

/* sunCloseLogFile: dropping the boxed sink closes it; the hashmap owns
its streams so the C helper (and sunLoggerFreeKeyValue) drop out. */

/// sunLoggerIsOutputRank — serial build: rank 0, always an output rank.
fn sunLoggerIsOutputRank(_logger: &SUNLogger, rank_ref: Option<&mut i32>) -> bool {
    if let Some(r) = rank_ref {
        *r = 0;
    }
    true
}

/// sunLoggerSetFilename — resolves `filename` to the stream slot value
/// the caller should install, opening the file and registering it in
/// the `filenames` map if it is not already there.
/// Ok(None) means "leave the current stream unchanged" (non-output
/// rank in C; unreachable in the serial build).
fn sunLoggerSetFilename(
    logger: &mut SUNLogger,
    filename: &str,
) -> Result<Option<SunLogStream>, SUNErrCode> {
    if !sunLoggerIsOutputRank(logger, None) {
        return Ok(None);
    }

    /* An empty or NULL filename disables output for this stream. */
    if filename.is_empty() {
        /* Don't close the file here, that is managed by the underlying hashmap */
        return Ok(Some(SunLogStream::None));
    }

    let map = match logger.filenames.as_mut() {
        Some(m) => m,
        /* C: SUNHashMap_GetValue on a NULL map yields SUNHASHMAP_ERROR */
        None => return Err(SUN_ERR_FILE_OPEN),
    };

    let mut fp: Option<&Box<dyn Write>> = None;
    let err = SUNHashMap_GetValue(map, filename, &mut fp);
    if err == SUNHASHMAP_ERROR {
        return Err(SUN_ERR_FILE_OPEN);
    } else if err == SUNHASHMAP_KEYNOTFOUND {
        let stream = match sunOpenLogFile(filename) {
            Some(s) => s,
            None => return Err(SUN_ERR_FILE_OPEN),
        };

        let err = SUNHashMap_Insert(map, filename, stream);
        if err != 0 {
            return Err(SUN_ERR_FILE_OPEN);
        }
    }

    Ok(Some(SunLogStream::Key(filename.to_string())))
}

/// sunLoggerSetFilePointer — install a caller-supplied sink (C FILE*).
fn sunLoggerSetFilePointer(
    logger: &SUNLogger,
    file_ptr: Option<Box<dyn Write>>,
) -> Option<SunLogStream> {
    if !sunLoggerIsOutputRank(logger, None) {
        return None;
    }
    Some(match file_ptr {
        Some(w) => SunLogStream::Direct(w),
        None => SunLogStream::None,
    })
}

/// Write `msg` to `stream` (C: `fprintf(fp, "%s", log_msg)`, whose
/// status is ignored). `Key` slots are resolved against `filenames`.
fn sunLoggerWriteStream(
    stream: &mut SunLogStream,
    filenames: &mut Option<SUNHashMap<Box<dyn Write>>>,
    msg: &str,
) {
    match stream {
        SunLogStream::None => {}
        SunLogStream::Stdout => {
            let _ = write!(std::io::stdout(), "{}", msg);
        }
        SunLogStream::Stderr => {
            let _ = write!(std::io::stderr(), "{}", msg);
        }
        SunLogStream::Key(key) => {
            if let Some(map) = filenames.as_mut() {
                if let Some(w) = map.get_mut(key) {
                    let _ = write!(w, "{}", msg);
                }
            }
        }
        SunLogStream::Direct(w) => {
            let _ = write!(w, "{}", msg);
        }
    }
}

/// Flush `stream` (C: `fflush(fp)`, status ignored).
fn sunLoggerFlushStream(
    stream: &mut SunLogStream,
    filenames: &mut Option<SUNHashMap<Box<dyn Write>>>,
) {
    match stream {
        SunLogStream::None => {}
        SunLogStream::Stdout => {
            let _ = std::io::stdout().flush();
        }
        SunLogStream::Stderr => {
            let _ = std::io::stderr().flush();
        }
        SunLogStream::Key(key) => {
            if let Some(map) = filenames.as_mut() {
                if let Some(w) = map.get_mut(key) {
                    let _ = w.flush();
                }
            }
        }
        SunLogStream::Direct(w) => {
            let _ = w.flush();
        }
    }
}

/// SUNLogger_Create (the `SUNComm comm` argument is dropped — serial
/// build only; the C non-MPI branch errors out for any non-NULL comm).
pub fn SUNLogger_Create(output_rank: i32, logger_ptr: &mut Option<SUNLogger>) -> SUNErrCode {
    *logger_ptr = None;

    let mut logger = SUNLogger {
        output_rank,
        content: None,

        /* use default routines */
        queuemsg: None,
        flush: None,
        destroy: None,

        /* set the output file handles */
        filenames: None,
        error_fp: SunLogStream::Stderr,
        warning_fp: SunLogStream::Stdout,
        debug_fp: SunLogStream::None,
        info_fp: SunLogStream::None,

        format: None,
    };

    if sunLoggerIsOutputRank(&logger, None) {
        /* We store the streams in a hash map so that we can ensure
        that we do not open a file twice if the same file is used
        for multiple output levels */
        SUNHashMap_New(SUN_DEFAULT_LOGFILE_HANDLES_, &mut logger.filenames);
    }

    *logger_ptr = Some(logger);
    SUN_SUCCESS
}

/// C atoi semantics for SUNLOGGER_OUTPUT_RANK: optional leading
/// whitespace and sign, then the longest run of digits; 0 otherwise.
fn sun_atoi(s: &str) -> i32 {
    let t = s.trim_start();
    let mut sign: i64 = 1;
    let mut digits = t;
    if let Some(rest) = t.strip_prefix('-') {
        sign = -1;
        digits = rest;
    } else if let Some(rest) = t.strip_prefix('+') {
        digits = rest;
    }
    let mut value: i64 = 0;
    for c in digits.chars() {
        match c.to_digit(10) {
            Some(d) => value = value.saturating_mul(10).saturating_add(d as i64),
            None => break,
        }
    }
    (sign * value) as i32
}

/// SUNLogger_CreateFromEnv (comm argument dropped, as in SUNLogger_Create).
pub fn SUNLogger_CreateFromEnv(logger_out: &mut Option<SUNLogger>) -> SUNErrCode {
    let mut err = SUN_SUCCESS;

    let output_rank_env = std::env::var("SUNLOGGER_OUTPUT_RANK").ok();
    let output_rank = match &output_rank_env {
        Some(s) => sun_atoi(s),
        None => 0,
    };
    let error_fname_env = std::env::var("SUNLOGGER_ERROR_FILENAME").ok();
    let warning_fname_env = std::env::var("SUNLOGGER_WARNING_FILENAME").ok();
    let info_fname_env = std::env::var("SUNLOGGER_INFO_FILENAME").ok();
    let debug_fname_env = std::env::var("SUNLOGGER_DEBUG_FILENAME").ok();

    let mut logger: Option<SUNLogger> = None;
    if SUNLogger_Create(output_rank, &mut logger) != SUN_SUCCESS {
        return SUN_ERR_CORRUPT;
    }

    if let Some(lg) = logger.as_mut() {
        /* C: do { ... } while (0) with early breaks */
        'set: {
            /* Only override the default logging if the env var is defined */
            if let Some(fname) = &error_fname_env {
                err = SUNLogger_SetErrorFilename(lg, fname);
                if err != SUN_SUCCESS {
                    break 'set;
                }
            }

            if let Some(fname) = &warning_fname_env {
                err = SUNLogger_SetWarningFilename(lg, fname);
                if err != SUN_SUCCESS {
                    break 'set;
                }
            }

            if let Some(fname) = &debug_fname_env {
                err = SUNLogger_SetDebugFilename(lg, fname);
                if err != SUN_SUCCESS {
                    break 'set;
                }
            }

            if let Some(fname) = &info_fname_env {
                err = SUNLogger_SetInfoFilename(lg, fname);
                if err != SUN_SUCCESS {
                    break 'set;
                }
            }
        }
    }

    if err != SUN_SUCCESS {
        SUNLogger_Destroy(&mut logger);
    } else {
        *logger_out = logger;
    }

    err
}

/// SUNLogger_SetErrorFilename
pub fn SUNLogger_SetErrorFilename(logger: &mut SUNLogger, error_filename: &str) -> SUNErrCode {
    match sunLoggerSetFilename(logger, error_filename) {
        Ok(Some(stream)) => {
            logger.error_fp = stream;
            SUN_SUCCESS
        }
        Ok(None) => SUN_SUCCESS,
        Err(e) => e,
    }
}

/// SUNLogger_SetErrorFile (C `FILE*` becomes an owned, optional sink)
pub fn SUNLogger_SetErrorFile(
    logger: &mut SUNLogger,
    error_fp: Option<Box<dyn Write>>,
) -> SUNErrCode {
    if let Some(stream) = sunLoggerSetFilePointer(logger, error_fp) {
        logger.error_fp = stream;
    }
    SUN_SUCCESS
}

/// SUNLogger_SetWarningFilename
pub fn SUNLogger_SetWarningFilename(
    logger: &mut SUNLogger,
    warning_filename: &str,
) -> SUNErrCode {
    match sunLoggerSetFilename(logger, warning_filename) {
        Ok(Some(stream)) => {
            logger.warning_fp = stream;
            SUN_SUCCESS
        }
        Ok(None) => SUN_SUCCESS,
        Err(e) => e,
    }
}

/// SUNLogger_SetWarningFile
pub fn SUNLogger_SetWarningFile(
    logger: &mut SUNLogger,
    warning_fp: Option<Box<dyn Write>>,
) -> SUNErrCode {
    if let Some(stream) = sunLoggerSetFilePointer(logger, warning_fp) {
        logger.warning_fp = stream;
    }
    SUN_SUCCESS
}

/// SUNLogger_SetInfoFilename
pub fn SUNLogger_SetInfoFilename(logger: &mut SUNLogger, info_filename: &str) -> SUNErrCode {
    match sunLoggerSetFilename(logger, info_filename) {
        Ok(Some(stream)) => {
            logger.info_fp = stream;
            SUN_SUCCESS
        }
        Ok(None) => SUN_SUCCESS,
        Err(e) => e,
    }
}

/// SUNLogger_SetInfoFile
pub fn SUNLogger_SetInfoFile(
    logger: &mut SUNLogger,
    info_fp: Option<Box<dyn Write>>,
) -> SUNErrCode {
    if let Some(stream) = sunLoggerSetFilePointer(logger, info_fp) {
        logger.info_fp = stream;
    }
    SUN_SUCCESS
}

/// SUNLogger_SetDebugFilename
pub fn SUNLogger_SetDebugFilename(logger: &mut SUNLogger, debug_filename: &str) -> SUNErrCode {
    match sunLoggerSetFilename(logger, debug_filename) {
        Ok(Some(stream)) => {
            logger.debug_fp = stream;
            SUN_SUCCESS
        }
        Ok(None) => SUN_SUCCESS,
        Err(e) => e,
    }
}

/// SUNLogger_SetDebugFile
pub fn SUNLogger_SetDebugFile(
    logger: &mut SUNLogger,
    debug_fp: Option<Box<dyn Write>>,
) -> SUNErrCode {
    if let Some(stream) = sunLoggerSetFilePointer(logger, debug_fp) {
        logger.debug_fp = stream;
    }
    SUN_SUCCESS
}

/// SUNLogger_QueueMsg — `msg_txt` arrives pre-formatted (the C
/// printf-style varargs are the documented Rust ABI adaptation).
pub fn SUNLogger_QueueMsg(
    logger: &mut SUNLogger,
    lvl: SUNLogLevel,
    scope: &str,
    label: &str,
    msg_txt: &str,
) -> SUNErrCode {
    let mut retval = SUN_SUCCESS;

    if let Some(queuemsg) = logger.queuemsg {
        retval = queuemsg(logger, lvl, scope, label, msg_txt);
    } else {
        /* Default implementation */
        let mut rank = 0;
        if sunLoggerIsOutputRank(logger, Some(&mut rank)) {
            let log_msg = sunCreateLogMessage(lvl, rank, scope, label, msg_txt);

            match lvl {
                SUN_LOGLEVEL_DEBUG => {
                    sunLoggerWriteStream(&mut logger.debug_fp, &mut logger.filenames, &log_msg);
                }
                SUN_LOGLEVEL_WARNING => {
                    sunLoggerWriteStream(&mut logger.warning_fp, &mut logger.filenames, &log_msg);
                }
                SUN_LOGLEVEL_INFO => {
                    sunLoggerWriteStream(&mut logger.info_fp, &mut logger.filenames, &log_msg);
                }
                SUN_LOGLEVEL_ERROR => {
                    sunLoggerWriteStream(&mut logger.error_fp, &mut logger.filenames, &log_msg);
                }
                _ => retval = SUN_ERR_UNREACHABLE,
            }
        }
    }

    retval
}

/// SUNLogger_Flush
pub fn SUNLogger_Flush(logger: &mut SUNLogger, lvl: SUNLogLevel) -> SUNErrCode {
    let mut retval = SUN_SUCCESS;

    if let Some(flush) = logger.flush {
        retval = flush(logger, lvl);
    } else {
        /* Default implementation */
        if sunLoggerIsOutputRank(logger, None) {
            match lvl {
                SUN_LOGLEVEL_DEBUG => {
                    sunLoggerFlushStream(&mut logger.debug_fp, &mut logger.filenames);
                }
                SUN_LOGLEVEL_WARNING => {
                    sunLoggerFlushStream(&mut logger.warning_fp, &mut logger.filenames);
                }
                SUN_LOGLEVEL_INFO => {
                    sunLoggerFlushStream(&mut logger.info_fp, &mut logger.filenames);
                }
                SUN_LOGLEVEL_ERROR => {
                    sunLoggerFlushStream(&mut logger.error_fp, &mut logger.filenames);
                }
                SUN_LOGLEVEL_ALL => {
                    sunLoggerFlushStream(&mut logger.debug_fp, &mut logger.filenames);
                    sunLoggerFlushStream(&mut logger.warning_fp, &mut logger.filenames);
                    sunLoggerFlushStream(&mut logger.info_fp, &mut logger.filenames);
                    sunLoggerFlushStream(&mut logger.error_fp, &mut logger.filenames);
                }
                _ => retval = SUN_ERR_UNREACHABLE,
            }
        }
    }

    retval
}

/// SUNLogger_GetOutputRank
pub fn SUNLogger_GetOutputRank(logger: &SUNLogger, output_rank: &mut i32) -> SUNErrCode {
    *output_rank = logger.output_rank;
    SUN_SUCCESS
}

/// SUNLogger_Destroy — the custom `destroy` callback (if any) runs
/// first; the handle itself is then dropped (Rust owns the memory the
/// C callback would have freed).
pub fn SUNLogger_Destroy(logger_ptr: &mut Option<SUNLogger>) -> SUNErrCode {
    let mut retval = SUN_SUCCESS;

    if let Some(logger) = logger_ptr.as_mut() {
        if let Some(destroy) = logger.destroy {
            retval = destroy(logger);
        } else {
            /* Default implementation */
            if sunLoggerIsOutputRank(logger, None) {
                SUNHashMap_Destroy(&mut logger.filenames);
            }
        }
    }
    *logger_ptr = None;

    retval
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    fn tmp_path(tag: &str) -> std::path::PathBuf {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .subsec_nanos();
        std::env::temp_dir().join(format!(
            "sundials_rs_logger_{}_{}_{}.log",
            tag,
            std::process::id(),
            nanos
        ))
    }

    fn read_to_string(path: &std::path::Path) -> String {
        let mut s = String::new();
        std::fs::File::open(path)
            .unwrap()
            .read_to_string(&mut s)
            .unwrap();
        s
    }

    #[test]
    fn test_each_level_writes_to_its_file() {
        let epath = tmp_path("err");
        let wpath = tmp_path("warn");
        let ipath = tmp_path("info");
        let dpath = tmp_path("dbg");

        let mut logger: Option<SUNLogger> = None;
        assert_eq!(SUNLogger_Create(0, &mut logger), SUN_SUCCESS);
        let lg = logger.as_mut().unwrap();

        let mut rank = -1;
        assert_eq!(SUNLogger_GetOutputRank(lg, &mut rank), SUN_SUCCESS);
        assert_eq!(rank, 0);

        assert_eq!(
            SUNLogger_SetErrorFilename(lg, epath.to_str().unwrap()),
            SUN_SUCCESS
        );
        assert_eq!(
            SUNLogger_SetWarningFilename(lg, wpath.to_str().unwrap()),
            SUN_SUCCESS
        );
        assert_eq!(
            SUNLogger_SetInfoFilename(lg, ipath.to_str().unwrap()),
            SUN_SUCCESS
        );
        assert_eq!(
            SUNLogger_SetDebugFilename(lg, dpath.to_str().unwrap()),
            SUN_SUCCESS
        );

        assert_eq!(
            SUNLogger_QueueMsg(lg, SUN_LOGLEVEL_ERROR, "scope", "label", "err msg"),
            SUN_SUCCESS
        );
        assert_eq!(
            SUNLogger_QueueMsg(lg, SUN_LOGLEVEL_WARNING, "scope", "label", "warn msg"),
            SUN_SUCCESS
        );
        assert_eq!(
            SUNLogger_QueueMsg(lg, SUN_LOGLEVEL_INFO, "scope", "label", "info msg"),
            SUN_SUCCESS
        );
        assert_eq!(
            SUNLogger_QueueMsg(lg, SUN_LOGLEVEL_DEBUG, "scope", "label", "dbg msg"),
            SUN_SUCCESS
        );

        /* queueing at NONE hits the C switch default */
        assert_eq!(
            SUNLogger_QueueMsg(lg, SUN_LOGLEVEL_NONE, "scope", "label", "x"),
            SUN_ERR_UNREACHABLE
        );

        assert_eq!(SUNLogger_Flush(lg, SUN_LOGLEVEL_ALL), SUN_SUCCESS);
        assert_eq!(SUNLogger_Destroy(&mut logger), SUN_SUCCESS);
        assert!(logger.is_none());

        assert_eq!(
            read_to_string(&epath),
            "[ERROR][rank 0][scope][label] err msg\n"
        );
        assert_eq!(
            read_to_string(&wpath),
            "[WARNING][rank 0][scope][label] warn msg\n"
        );
        assert_eq!(
            read_to_string(&ipath),
            "[INFO][rank 0][scope][label] info msg\n"
        );
        assert_eq!(
            read_to_string(&dpath),
            "[DEBUG][rank 0][scope][label] dbg msg\n"
        );

        for p in [&epath, &wpath, &ipath, &dpath] {
            let _ = std::fs::remove_file(p);
        }
    }

    #[test]
    fn test_shared_filename_is_opened_once() {
        /* the same file for two levels must be deduplicated through the
        hashmap: a second fopen("w+") would truncate the first line */
        let path = tmp_path("shared");
        let fname = path.to_str().unwrap();

        let mut logger: Option<SUNLogger> = None;
        assert_eq!(SUNLogger_Create(0, &mut logger), SUN_SUCCESS);
        let lg = logger.as_mut().unwrap();

        assert_eq!(SUNLogger_SetInfoFilename(lg, fname), SUN_SUCCESS);
        assert_eq!(
            SUNLogger_QueueMsg(lg, SUN_LOGLEVEL_INFO, "s", "l", "first"),
            SUN_SUCCESS
        );
        assert_eq!(SUNLogger_SetDebugFilename(lg, fname), SUN_SUCCESS);
        assert_eq!(
            SUNLogger_QueueMsg(lg, SUN_LOGLEVEL_DEBUG, "s", "l", "second"),
            SUN_SUCCESS
        );
        assert_eq!(SUNLogger_Flush(lg, SUN_LOGLEVEL_ALL), SUN_SUCCESS);
        assert_eq!(SUNLogger_Destroy(&mut logger), SUN_SUCCESS);

        assert_eq!(
            read_to_string(&path),
            "[INFO][rank 0][s][l] first\n[DEBUG][rank 0][s][l] second\n"
        );
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_empty_filename_disables_stream() {
        let path = tmp_path("disable");
        let fname = path.to_str().unwrap();

        let mut logger: Option<SUNLogger> = None;
        assert_eq!(SUNLogger_Create(0, &mut logger), SUN_SUCCESS);
        let lg = logger.as_mut().unwrap();

        assert_eq!(SUNLogger_SetInfoFilename(lg, fname), SUN_SUCCESS);
        assert_eq!(
            SUNLogger_QueueMsg(lg, SUN_LOGLEVEL_INFO, "s", "l", "kept"),
            SUN_SUCCESS
        );
        /* An empty filename disables output for this stream */
        assert_eq!(SUNLogger_SetInfoFilename(lg, ""), SUN_SUCCESS);
        assert_eq!(
            SUNLogger_QueueMsg(lg, SUN_LOGLEVEL_INFO, "s", "l", "dropped"),
            SUN_SUCCESS
        );
        assert_eq!(SUNLogger_Flush(lg, SUN_LOGLEVEL_INFO), SUN_SUCCESS);
        assert_eq!(SUNLogger_Destroy(&mut logger), SUN_SUCCESS);

        assert_eq!(read_to_string(&path), "[INFO][rank 0][s][l] kept\n");
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_set_file_direct_sink() {
        /* SUNLogger_Set*File installs a caller-supplied sink */
        let path = tmp_path("direct");
        let file = std::fs::File::create(&path).unwrap();

        let mut logger: Option<SUNLogger> = None;
        assert_eq!(SUNLogger_Create(0, &mut logger), SUN_SUCCESS);
        let lg = logger.as_mut().unwrap();

        assert_eq!(SUNLogger_SetErrorFile(lg, Some(Box::new(file))), SUN_SUCCESS);
        assert_eq!(
            SUNLogger_QueueMsg(lg, SUN_LOGLEVEL_ERROR, "s", "l", "direct"),
            SUN_SUCCESS
        );
        assert_eq!(SUNLogger_Flush(lg, SUN_LOGLEVEL_ERROR), SUN_SUCCESS);

        /* NULL file pointer disables the stream */
        assert_eq!(SUNLogger_SetErrorFile(lg, None), SUN_SUCCESS);
        assert_eq!(
            SUNLogger_QueueMsg(lg, SUN_LOGLEVEL_ERROR, "s", "l", "gone"),
            SUN_SUCCESS
        );
        assert_eq!(SUNLogger_Destroy(&mut logger), SUN_SUCCESS);

        assert_eq!(read_to_string(&path), "[ERROR][rank 0][s][l] direct\n");
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_custom_queuemsg_override() {
        fn custom(
            logger: &mut SUNLogger,
            lvl: SUNLogLevel,
            scope: &str,
            label: &str,
            msg: &str,
        ) -> SUNErrCode {
            logger.content = Some(Box::new(format!("{:?}|{}|{}|{}", lvl, scope, label, msg)));
            SUN_SUCCESS
        }

        let mut logger: Option<SUNLogger> = None;
        assert_eq!(SUNLogger_Create(0, &mut logger), SUN_SUCCESS);
        let lg = logger.as_mut().unwrap();
        lg.queuemsg = Some(custom);

        assert_eq!(
            SUNLogger_QueueMsg(lg, SUN_LOGLEVEL_INFO, "sc", "lb", "hello"),
            SUN_SUCCESS
        );
        let seen = lg
            .content
            .as_ref()
            .unwrap()
            .downcast_ref::<String>()
            .unwrap();
        assert_eq!(seen, "SUN_LOGLEVEL_INFO|sc|lb|hello");
        assert_eq!(SUNLogger_Destroy(&mut logger), SUN_SUCCESS);
    }

    #[test]
    fn test_create_from_env() {
        let path = tmp_path("env");
        std::env::set_var("SUNLOGGER_OUTPUT_RANK", "0");
        std::env::set_var("SUNLOGGER_INFO_FILENAME", path.to_str().unwrap());

        let mut logger: Option<SUNLogger> = None;
        assert_eq!(SUNLogger_CreateFromEnv(&mut logger), SUN_SUCCESS);
        std::env::remove_var("SUNLOGGER_OUTPUT_RANK");
        std::env::remove_var("SUNLOGGER_INFO_FILENAME");

        let lg = logger.as_mut().unwrap();
        assert_eq!(
            SUNLogger_QueueMsg(lg, SUN_LOGLEVEL_INFO, "s", "l", "from env"),
            SUN_SUCCESS
        );
        assert_eq!(SUNLogger_Flush(lg, SUN_LOGLEVEL_INFO), SUN_SUCCESS);
        assert_eq!(SUNLogger_Destroy(&mut logger), SUN_SUCCESS);

        assert_eq!(read_to_string(&path), "[INFO][rank 0][s][l] from env\n");
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_atoi() {
        assert_eq!(sun_atoi("42"), 42);
        assert_eq!(sun_atoi("  -7"), -7);
        assert_eq!(sun_atoi("+3abc"), 3);
        assert_eq!(sun_atoi("abc"), 0);
        assert_eq!(sun_atoi(""), 0);
    }
}
