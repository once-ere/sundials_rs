/* -----------------------------------------------------------------
 * Translated from src/sundials/sundials_context.c (SUNDIALS 7.7.0).
 * In the serial pure-Rust build the context carries no MPI
 * communicator; the profiler and logger backends are compiled out
 * (equivalent to a C build without SUNDIALS_BUILD_WITH_PROFILING /
 * SUNDIALS_LOGGING_LEVEL). The object is retained so the public API
 * signatures (`CVodeCreate(lmm, sunctx)`) match the C library.
 * -----------------------------------------------------------------*/
use crate::sundials_errors::*;

#[derive(Default, Debug, Clone)]
pub struct SUNContext {
    pub last_err: SUNErrCode,
}

/// SUNContext_Create (comm argument dropped: serial build only).
pub fn SUNContext_Create() -> SUNContext {
    SUNContext::default()
}

/// SUNContext_Free — RAII makes this a no-op; kept for API parity.
pub fn SUNContext_Free(_sunctx: &mut SUNContext) -> SUNErrCode {
    SUN_SUCCESS
}

/// SUNContext_GetLastError
pub fn SUNContext_GetLastError(sunctx: &mut SUNContext) -> SUNErrCode {
    let e = sunctx.last_err;
    sunctx.last_err = SUN_SUCCESS;
    e
}

/// SUNContext_PeekLastError
pub fn SUNContext_PeekLastError(sunctx: &SUNContext) -> SUNErrCode {
    sunctx.last_err
}
