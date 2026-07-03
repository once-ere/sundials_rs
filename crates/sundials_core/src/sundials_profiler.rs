/* -----------------------------------------------------------------
 * Translated from src/sundials/sundials_profiler.c together with
 * include/sundials/sundials_profiler.h and
 * src/sundials/sundials_profiler_impl.h (SUNDIALS 7.7.0).
 *
 * Serial build only: the MPI branches (sunCollectTimers,
 * sunTimerStructReduceMaxAndSum, MPI_Comm_dup/free) and the Caliper
 * backend are compiled out, exactly like a C build without
 * SUNDIALS_MPI_ENABLED / SUNDIALS_CALIPER_ENABLED. The
 * SUNDIALS_MARK_BEGIN/END macros expand to SUNProfiler_Begin/End
 * (the SUNDIALS_ENABLE_PROFILING branch of sundials_profiler.h).
 *
 * Rust adaptations (documented, everything else follows the C
 * line-for-line):
 *  - Timing uses std::time::Instant instead of
 *    clock_gettime(CLOCK_MONOTONIC)/QueryPerformanceCounter; Instant
 *    is the std monotonic clock. The elapsed-time bookkeeping
 *    (elapsed += sec_difference + nsec_difference * 1e-9; average and
 *    maximum mirror elapsed in serial; count) is kept exactly as in
 *    sunStopTiming. A timer whose tic was never set (C: tv_sec ==
 *    tv_nsec == 0) accumulates nothing instead of a time measured
 *    from the clock epoch.
 *  - The timer registry is a std::collections::HashMap<String,
 *    sunTimerStruct> instead of the C SUNHashMap, chosen locally to
 *    avoid coupling to sundials_hashmap.rs. The
 *    SUNPROFILER_MAX_ENTRIES capacity (default 2560) is enforced
 *    manually so a full map still yields SUN_ERR_PROFILER_MAPFULL,
 *    and the internal SUNHASHMAP_ERROR paths (SUN_ERR_PROFILER_MAPGET
 *    / SUN_ERR_PROFILER_MAPINSERT / SUN_ERR_PROFILER_MAPSORT) cannot
 *    occur.
 *  - `FILE* fp` becomes `&mut dyn std::io::Write` (workspace donor
 *    pattern, cvode_io.rs).
 *  - The SUN_ERR_PROFILER_* codes are not present in
 *    sundials_errors.rs (which renumbered the tail of the C enum);
 *    they are defined here with workspace-unique values in the C
 *    error range [-10000, -1000].
 *  - `typedef int SUNComm` / SUN_COMM_NULL (serial branch of
 *    sundials_types.h) are defined here; nothing else in the
 *    workspace declares them. SUNProfiler_Create returns -1 for a
 *    non-NULL comm exactly as the non-MPI C branch does.
 *  - C null-handle checks (`if (!p) return SUN_ERR_ARG_CORRUPT`)
 *    vanish where the argument is a Rust reference; Create/Free keep
 *    the C out-pointer shape via &mut Option<SUNProfiler>.
 *  - SUNProfiler_GetTimerResolution reports 1e-9 s (Instant counts
 *    nanoseconds); C queries clock_getres(CLOCK_MONOTONIC).
 *  - SUNDIALS_GIT_VERSION is "" as in a configured release tarball.
 * -----------------------------------------------------------------*/

use std::collections::HashMap;
use std::time::Instant;

use crate::sundials_errors::*;
use crate::sundials_utils::{fmt_f, fmt_g};

/// `SUNComm` (sundials_types.h, serial branch: `typedef int SUNComm`).
pub type SUNComm = i32;
/// `SUN_COMM_NULL` (sundials_types.h, serial branch).
pub const SUN_COMM_NULL: SUNComm = 0;

/* SUN_ERR_PROFILER_* (sundials_errors.h). In the C enum these are
   -9981..-9977; sundials_errors.rs renumbered the codes that follow
   the profiler block, so the C values would collide with
   SUN_ERR_SUNCTX_CORRUPT etc. Unique values in the same C error
   range are used instead (name lookup still reports "unknown
   error" from SUNGetErrMsg, which never handled these codes here). */
pub const SUN_ERR_PROFILER_MAPFULL: SUNErrCode = -9976;
pub const SUN_ERR_PROFILER_MAPGET: SUNErrCode = -9975;
pub const SUN_ERR_PROFILER_MAPINSERT: SUNErrCode = -9974;
pub const SUN_ERR_PROFILER_MAPKEYNOTFOUND: SUNErrCode = -9973;
pub const SUN_ERR_PROFILER_MAPSORT: SUNErrCode = -9972;

/// SUNDIALS_ROOT_TIMER (sundials_profiler_impl.h)
pub const SUNDIALS_ROOT_TIMER: &str = "From profiler epoch";

/// SUNDIALS_GIT_VERSION (sundials_config.h) — empty in release builds.
const SUNDIALS_GIT_VERSION: &str = "";

/*
  sunTimerStruct (sundials_profiler_impl.h).
  A private structure holding timing information.
 */
#[derive(Debug, Clone)]
pub struct sunTimerStruct {
    /// C: `sunTimespec *tic` (the matching `toc` is taken at stop
    /// time and never stored across calls in the Rust port).
    tic: Option<Instant>,
    pub average: f64,
    pub maximum: f64,
    pub elapsed: f64,
    /// C `long count`
    pub count: i64,
}

/// sunTimerStructNew
fn sunTimerStructNew() -> sunTimerStruct {
    sunTimerStruct { tic: None, elapsed: 0.0, average: 0.0, maximum: 0.0, count: 0 }
}

/// sunStartTiming
fn sunStartTiming(entry: &mut sunTimerStruct) {
    /* sunclock_gettime_monotonic(entry->tic) */
    entry.tic = Some(Instant::now());
}

/// sunStopTiming
fn sunStopTiming(entry: &mut sunTimerStruct) {
    if let Some(tic) = entry.tic {
        /* toc - tic split into whole seconds and nanoseconds, summed
           exactly as in C: s_difference + ns_difference * 1e-9 */
        let d = tic.elapsed();
        let s_difference = d.as_secs() as f64;
        let ns_difference = d.subsec_nanos() as f64;
        entry.elapsed += s_difference + ns_difference * 1e-9;
    }
    entry.average = entry.elapsed;
    entry.maximum = entry.elapsed;
}

/// sunResetTiming
fn sunResetTiming(entry: &mut sunTimerStruct) {
    entry.tic = None;
    entry.elapsed = 0.0;
    entry.average = 0.0;
    entry.maximum = 0.0;
    entry.count = 0;
}

/// struct SUNProfiler_ (sundials_profiler_impl.h); the C `SUNProfiler`
/// handle is a pointer to it.
pub struct SUNProfiler_ {
    pub comm: SUNComm,
    pub title: String,
    /// C: `SUNHashMap map` — see module header for the HashMap choice.
    map: HashMap<String, sunTimerStruct>,
    /// Capacity of the C hashmap (SUNPROFILER_MAX_ENTRIES).
    max_entries: usize,
    overhead: sunTimerStruct,
    pub sundials_time: f64,
}

pub type SUNProfiler = SUNProfiler_;

/// atoi (for the SUNPROFILER_MAX_ENTRIES environment variable):
/// parse a leading optionally-signed integer, 0 on failure.
fn sun_atoi(s: &str) -> i32 {
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

/// SUNProfiler_Create
pub fn SUNProfiler_Create(comm: SUNComm, title: &str, p: &mut Option<SUNProfiler>) -> SUNErrCode {
    /* C: *p = profiler = malloc(...); allocation failure paths
       (return SUN_SUCCESS with *p == NULL, SUN_ERR_MALLOC_FAIL) are
       unreachable in Rust. */
    let mut profiler = SUNProfiler_ {
        comm: SUN_COMM_NULL,
        title: String::new(),
        map: HashMap::new(),
        max_entries: 0,
        overhead: sunTimerStructNew(),
        sundials_time: 0.0,
    };

    sunStartTiming(&mut profiler.overhead);

    /* Check to see if max entries env variable was set, and use if it was. */
    let mut max_entries = 2560;
    if let Ok(max_entries_env) = std::env::var("SUNPROFILER_MAX_ENTRIES") {
        max_entries = sun_atoi(&max_entries_env);
    }
    if max_entries <= 0 {
        max_entries = 2560;
    }

    /* Create the hashmap used to store the timers */
    profiler.max_entries = max_entries as usize;

    /* Attach the comm (non-MPI branch: any non-NULL comm is an error). */
    if comm != SUN_COMM_NULL {
        *p = None;
        return -1;
    }
    profiler.comm = SUN_COMM_NULL;

    /* Copy the title of the profiler */
    profiler.title = title.to_string();

    /* Initialize the overall timer to 0. */
    profiler.sundials_time = 0.0;

    /* SUNDIALS_MARK_BEGIN(profiler, SUNDIALS_ROOT_TIMER) */
    let _ = SUNProfiler_Begin(&mut profiler, SUNDIALS_ROOT_TIMER);
    sunStopTiming(&mut profiler.overhead);

    *p = Some(profiler);
    SUN_SUCCESS
}

/// SUNProfiler_Free
pub fn SUNProfiler_Free(p: &mut Option<SUNProfiler>) -> SUNErrCode {
    if p.is_none() {
        return SUN_SUCCESS;
    }

    /* SUNDIALS_MARK_END(*p, SUNDIALS_ROOT_TIMER) */
    if let Some(profiler) = p.as_mut() {
        let _ = SUNProfiler_End(profiler, SUNDIALS_ROOT_TIMER);
    }

    /* map, overhead and title are freed by drop */
    *p = None;

    SUN_SUCCESS
}

/// SUNProfiler_Begin
pub fn SUNProfiler_Begin(p: &mut SUNProfiler, name: &str) -> SUNErrCode {
    sunStartTiming(&mut p.overhead);

    if !p.map.contains_key(name) {
        /* C: SUNHashMap_Insert fails with SUNHASHMAP_DUPLICATE when the
           map is at capacity -> SUN_ERR_PROFILER_MAPFULL. */
        if p.map.len() >= p.max_entries {
            sunStopTiming(&mut p.overhead);
            return SUN_ERR_PROFILER_MAPFULL;
        }
        p.map.insert(name.to_string(), sunTimerStructNew());
    }

    if let Some(timer) = p.map.get_mut(name) {
        timer.count += 1;
        sunStartTiming(timer);
    }

    sunStopTiming(&mut p.overhead);
    SUN_SUCCESS
}

/// SUNProfiler_End
pub fn SUNProfiler_End(p: &mut SUNProfiler, name: &str) -> SUNErrCode {
    sunStartTiming(&mut p.overhead);

    match p.map.get_mut(name) {
        None => {
            sunStopTiming(&mut p.overhead);
            SUN_ERR_PROFILER_MAPKEYNOTFOUND
        }
        Some(timer) => {
            sunStopTiming(timer);
            sunStopTiming(&mut p.overhead);
            SUN_SUCCESS
        }
    }
}

/// SUNProfiler_GetTimerResolution — std::time::Instant counts whole
/// nanoseconds, so report 1e-9 s (C queries clock_getres(CLOCK_MONOTONIC)).
pub fn SUNProfiler_GetTimerResolution(_p: &SUNProfiler, resolution: &mut f64) -> SUNErrCode {
    *resolution = 1e-9;
    SUN_SUCCESS
}

/// SUNProfiler_GetElapsedTime
pub fn SUNProfiler_GetElapsedTime(p: &SUNProfiler, name: &str, time: &mut f64) -> SUNErrCode {
    match p.map.get(name) {
        None => -1,
        Some(timer) => {
            *time = timer.elapsed;
            SUN_SUCCESS
        }
    }
}

/// SUNProfiler_Reset
pub fn SUNProfiler_Reset(p: &mut SUNProfiler) -> SUNErrCode {
    /* Reset the overhead timer */
    sunResetTiming(&mut p.overhead);
    sunStartTiming(&mut p.overhead);

    /* Reset all timers */
    for timer in p.map.values_mut() {
        sunResetTiming(timer);
    }

    /* Reset the overall timer. */
    p.sundials_time = 0.0;

    /* SUNDIALS_MARK_BEGIN(p, SUNDIALS_ROOT_TIMER) */
    let _ = SUNProfiler_Begin(p, SUNDIALS_ROOT_TIMER);
    sunStopTiming(&mut p.overhead);

    SUN_SUCCESS
}

/// SUNProfiler_Print
pub fn SUNProfiler_Print(p: &mut SUNProfiler, fp: &mut dyn std::io::Write) -> SUNErrCode {
    /* rank == 0 always (serial build) */

    sunStartTiming(&mut p.overhead);

    /* Get the total SUNDIALS time up to this point:
       SUNDIALS_MARK_END / SUNDIALS_MARK_BEGIN on the root timer */
    let _ = SUNProfiler_End(p, SUNDIALS_ROOT_TIMER);
    let _ = SUNProfiler_Begin(p, SUNDIALS_ROOT_TIMER);

    match p.map.get(SUNDIALS_ROOT_TIMER) {
        None => return SUN_ERR_PROFILER_MAPKEYNOTFOUND,
        Some(timer) => p.sundials_time = timer.elapsed,
    }

    {
        let mut resolution = 0.0;
        /* Sort the timers in descending order (SUNHashMap_Sort with
           sunCompareTimes) */
        let mut sorted: Vec<(String, sunTimerStruct)> =
            p.map.iter().map(|(k, v)| (k.clone(), v.clone())).collect();
        sorted.sort_by(sunCompareTimes);
        let _ = SUNProfiler_GetTimerResolution(p, &mut resolution);
        let _ = write!(fp, "\n{}\n", "=".repeat(112));
        let _ = write!(fp, "SUNDIALS GIT VERSION: {}\n", SUNDIALS_GIT_VERSION);
        let _ = write!(fp, "SUNDIALS PROFILER: {}\n", p.title);
        let _ = write!(fp, "TIMER RESOLUTION: {}s\n", fmt_g(resolution, 0, 6));
        let _ = write!(
            fp,
            "{:<40}\t % time (inclusive) \t max/rank \t average/rank \t count \n",
            "RESULTS:"
        );
        let _ = write!(fp, "{}\n", "=".repeat(112));

        /* Print all the other timers out */
        for (key, ts) in &sorted {
            sunPrintTimer(key, ts, fp, p);
        }
    }

    sunStopTiming(&mut p.overhead);

    {
        /* Print out the total time and the profiler overhead */
        let _ = write!(
            fp,
            "{:<40}\t {}% \t         {}s \t -- \t\t -- \n",
            "Est. profiler overhead",
            fmt_f(p.overhead.elapsed / p.sundials_time, 6, 2),
            fmt_f(p.overhead.elapsed, 0, 6)
        );

        /* End of output */
        let _ = write!(fp, "\n");
    }

    SUN_SUCCESS
}

/* Print out the: timer name, percentage of exec time (based on the max),
   max across ranks, average across ranks, and the timer counter. */
fn sunPrintTimer(key: &str, ts: &sunTimerStruct, fp: &mut dyn std::io::Write, p: &SUNProfiler) {
    let maximum = ts.maximum;
    let average = ts.average;
    let percent = if key != SUNDIALS_ROOT_TIMER {
        maximum / p.sundials_time * 100.0
    } else {
        100.0
    };
    let _ = write!(
        fp,
        "{:<40}\t {}% \t         {}s \t {}s \t {}\n",
        key,
        fmt_f(percent, 6, 2),
        fmt_f(maximum, 0, 6),
        fmt_f(average, 0, 6),
        ts.count
    );
}

/* Comparator for qsort that compares key-value pairs
   based on the maximum time in the sunTimerStruct.
   (The C NULL-slot cases cannot occur: every entry is populated.) */
fn sunCompareTimes(
    l: &(String, sunTimerStruct),
    r: &(String, sunTimerStruct),
) -> std::cmp::Ordering {
    let left_max = l.1.maximum;
    let right_max = r.1.maximum;

    if left_max < right_max {
        return std::cmp::Ordering::Greater;
    }
    if left_max > right_max {
        return std::cmp::Ordering::Less;
    }

    std::cmp::Ordering::Equal
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::thread::sleep;
    use std::time::Duration;

    #[test]
    fn create_rejects_non_null_comm() {
        let mut p: Option<SUNProfiler> = None;
        assert_eq!(SUNProfiler_Create(1, "bad", &mut p), -1);
        assert!(p.is_none());
    }

    #[test]
    fn begin_end_accumulates_positive_elapsed_time() {
        let mut opt: Option<SUNProfiler> = None;
        assert_eq!(SUNProfiler_Create(SUN_COMM_NULL, "test title", &mut opt), SUN_SUCCESS);
        let p = opt.as_mut().expect("profiler created");

        assert_eq!(SUNProfiler_Begin(p, "region A"), SUN_SUCCESS);
        sleep(Duration::from_millis(2));
        assert_eq!(SUNProfiler_End(p, "region A"), SUN_SUCCESS);

        let mut t = 0.0;
        assert_eq!(SUNProfiler_GetElapsedTime(p, "region A", &mut t), SUN_SUCCESS);
        assert!(t > 0.0, "elapsed time must be positive, got {}", t);

        /* second Begin/End accumulates into the same timer */
        assert_eq!(SUNProfiler_Begin(p, "region A"), SUN_SUCCESS);
        sleep(Duration::from_millis(2));
        assert_eq!(SUNProfiler_End(p, "region A"), SUN_SUCCESS);
        let mut t2 = 0.0;
        assert_eq!(SUNProfiler_GetElapsedTime(p, "region A", &mut t2), SUN_SUCCESS);
        assert!(t2 > t, "elapsed time must accumulate: {} vs {}", t2, t);
        assert_eq!(p.map.get("region A").map(|ts| ts.count), Some(2));

        /* ending an unknown timer is an error */
        assert_eq!(SUNProfiler_End(p, "no such region"), SUN_ERR_PROFILER_MAPKEYNOTFOUND);
        /* elapsed time of an unknown timer is an error (C returns -1) */
        let mut tu = 0.0;
        assert_eq!(SUNProfiler_GetElapsedTime(p, "no such region", &mut tu), -1);

        assert_eq!(SUNProfiler_Free(&mut opt), SUN_SUCCESS);
        assert!(opt.is_none());
        /* freeing twice is harmless */
        assert_eq!(SUNProfiler_Free(&mut opt), SUN_SUCCESS);
    }

    #[test]
    fn print_writes_region_names() {
        let mut opt: Option<SUNProfiler> = None;
        assert_eq!(SUNProfiler_Create(SUN_COMM_NULL, "my profiler", &mut opt), SUN_SUCCESS);
        let p = opt.as_mut().expect("profiler created");

        assert_eq!(SUNProfiler_Begin(p, "solve phase"), SUN_SUCCESS);
        sleep(Duration::from_millis(1));
        assert_eq!(SUNProfiler_End(p, "solve phase"), SUN_SUCCESS);

        let mut buf: Vec<u8> = Vec::new();
        assert_eq!(SUNProfiler_Print(p, &mut buf), SUN_SUCCESS);
        let out = String::from_utf8(buf).expect("utf8 output");

        assert!(out.contains("SUNDIALS PROFILER: my profiler"), "{}", out);
        assert!(out.contains("solve phase"), "{}", out);
        assert!(out.contains(SUNDIALS_ROOT_TIMER), "{}", out);
        assert!(out.contains("Est. profiler overhead"), "{}", out);
        assert!(out.contains("RESULTS:"), "{}", out);
        assert!(out.contains(&"=".repeat(112)), "{}", out);
    }

    #[test]
    fn reset_zeroes_timers() {
        let mut opt: Option<SUNProfiler> = None;
        assert_eq!(SUNProfiler_Create(SUN_COMM_NULL, "reset", &mut opt), SUN_SUCCESS);
        let p = opt.as_mut().expect("profiler created");

        assert_eq!(SUNProfiler_Begin(p, "r"), SUN_SUCCESS);
        sleep(Duration::from_millis(1));
        assert_eq!(SUNProfiler_End(p, "r"), SUN_SUCCESS);

        assert_eq!(SUNProfiler_Reset(p), SUN_SUCCESS);
        let mut t = 1.0;
        assert_eq!(SUNProfiler_GetElapsedTime(p, "r", &mut t), SUN_SUCCESS);
        assert_eq!(t, 0.0);
        assert_eq!(p.sundials_time, 0.0);
        /* root timer restarted by Reset (count back to 1) */
        assert_eq!(p.map.get(SUNDIALS_ROOT_TIMER).map(|ts| ts.count), Some(1));
    }

    #[test]
    fn map_full_returns_error() {
        let mut opt: Option<SUNProfiler> = None;
        assert_eq!(SUNProfiler_Create(SUN_COMM_NULL, "full", &mut opt), SUN_SUCCESS);
        let p = opt.as_mut().expect("profiler created");

        /* default capacity 2560; the root timer occupies one slot */
        for i in 0..(2560 - 1) {
            assert_eq!(SUNProfiler_Begin(p, &format!("t{}", i)), SUN_SUCCESS);
        }
        assert_eq!(SUNProfiler_Begin(p, "one too many"), SUN_ERR_PROFILER_MAPFULL);
        /* an existing timer can still be re-entered */
        assert_eq!(SUNProfiler_Begin(p, "t0"), SUN_SUCCESS);
    }
}
