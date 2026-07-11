/* -----------------------------------------------------------------
 * Translated from include/sundials/sundials_errors.h and
 * src/sundials/sundials_errors.c (SUNDIALS 7.7.0).
 * Only the error codes and name lookup are meaningful in a pure-Rust
 * serial build (the C error-handler stack manipulates void* handles).
 * -----------------------------------------------------------------*/

pub type SUNErrCode = i32;

pub const SUN_SUCCESS: SUNErrCode = 0;

/* sundials_errors.h: SUN_ERR_* codes live in [-10000, -1000] */
pub const SUN_ERR_MINIMUM: SUNErrCode = -10000;
pub const SUN_ERR_ARG_CORRUPT: SUNErrCode = -9999;
pub const SUN_ERR_ARG_INCOMPATIBLE: SUNErrCode = -9998;
pub const SUN_ERR_ARG_OUTOFRANGE: SUNErrCode = -9997;
pub const SUN_ERR_ARG_WRONGTYPE: SUNErrCode = -9996;
pub const SUN_ERR_ARG_DIMSMISMATCH: SUNErrCode = -9995;
pub const SUN_ERR_GENERIC: SUNErrCode = -9994;
pub const SUN_ERR_CORRUPT: SUNErrCode = -9993;
pub const SUN_ERR_OUTOFRANGE: SUNErrCode = -9992;
pub const SUN_ERR_FILE_OPEN: SUNErrCode = -9991;
pub const SUN_ERR_OP_FAIL: SUNErrCode = -9990;
pub const SUN_ERR_MEM_FAIL: SUNErrCode = -9989;
pub const SUN_ERR_MALLOC_FAIL: SUNErrCode = -9988;
pub const SUN_ERR_EXT_FAIL: SUNErrCode = -9987;
pub const SUN_ERR_DESTROY_FAIL: SUNErrCode = -9986;
pub const SUN_ERR_NOT_IMPLEMENTED: SUNErrCode = -9985;
pub const SUN_ERR_USER_FCN_FAIL: SUNErrCode = -9984;
pub const SUN_ERR_DATANODE_NODENOTFOUND: SUNErrCode = -9983;
pub const SUN_ERR_PROFILER_MAPFULL: SUNErrCode = -9982;
pub const SUN_ERR_PROFILER_MAPGET: SUNErrCode = -9981;
pub const SUN_ERR_PROFILER_MAPINSERT: SUNErrCode = -9980;
pub const SUN_ERR_PROFILER_MAPKEYNOTFOUND: SUNErrCode = -9979;
pub const SUN_ERR_PROFILER_MAPSORT: SUNErrCode = -9978;
pub const SUN_ERR_ADJOINT_STEPPERFAILED: SUNErrCode = -9977;
pub const SUN_ERR_ADJOINT_STEPPERINVALIDSTOP: SUNErrCode = -9976;
pub const SUN_ERR_CHECKPOINT_NOT_FOUND: SUNErrCode = -9975;
pub const SUN_ERR_CHECKPOINT_MISMATCH: SUNErrCode = -9974;
pub const SUN_ERR_SUNCTX_CORRUPT: SUNErrCode = -9973;
pub const SUN_ERR_MPI_FAIL: SUNErrCode = -9972;
pub const SUN_ERR_UNREACHABLE: SUNErrCode = -9971;
pub const SUN_ERR_UNKNOWN: SUNErrCode = -9970;
pub const SUN_ERR_MAXIMUM: SUNErrCode = -1000;

/// SUNGetErrMsg (sundials_errors.c)
pub fn SUNGetErrMsg(code: SUNErrCode) -> &'static str {
    match code {
        SUN_SUCCESS => "success",
        SUN_ERR_ARG_CORRUPT => "argument provided is NULL or corrupted",
        SUN_ERR_ARG_INCOMPATIBLE => "argument provided is not compatible",
        SUN_ERR_ARG_OUTOFRANGE => "argument is out of the valid range",
        SUN_ERR_ARG_WRONGTYPE => "argument provided is not the right type",
        SUN_ERR_ARG_DIMSMISMATCH => "argument dimensions do not agree",
        SUN_ERR_GENERIC => "an error occurred",
        SUN_ERR_CORRUPT => "value is NULL or corrupt",
        SUN_ERR_OUTOFRANGE => "Value is out of the expected range",
        SUN_ERR_FILE_OPEN => "Unable to open file",
        SUN_ERR_OP_FAIL => "an operation failed",
        SUN_ERR_MEM_FAIL => "a memory operation failed",
        SUN_ERR_MALLOC_FAIL => "malloc returned NULL",
        SUN_ERR_EXT_FAIL => "a failure occurred in an external library",
        SUN_ERR_DESTROY_FAIL => "a destroy function returned an error",
        SUN_ERR_NOT_IMPLEMENTED => "operation is not implemented: function pointer is NULL",
        SUN_ERR_USER_FCN_FAIL => "the user provided callback function failed",
        SUN_ERR_DATANODE_NODENOTFOUND => "the data node could not be found",
        SUN_ERR_PROFILER_MAPFULL => {
            "the number of profiler entries exceeded SUNPROFILER_MAX_ENTRIES"
        }
        SUN_ERR_PROFILER_MAPGET => "unknown error getting SUNProfiler timer",
        SUN_ERR_PROFILER_MAPINSERT => "unknown error inserting SUNProfiler timer",
        SUN_ERR_PROFILER_MAPKEYNOTFOUND => "timer was not found in SUNProfiler",
        SUN_ERR_PROFILER_MAPSORT => "error sorting SUNProfiler map",
        SUN_ERR_ADJOINT_STEPPERFAILED => {
            "SUNStepper stopped without successfully reaching the requested \
             output time when solving the adjoint system"
        }
        SUN_ERR_ADJOINT_STEPPERINVALIDSTOP => {
            "SUNStepper stopped with a flag not supported by the adjoint solver"
        }
        SUN_ERR_CHECKPOINT_NOT_FOUND => "the requested checkpoint was not found",
        SUN_ERR_CHECKPOINT_MISMATCH => {
            "the expected time for the checkpoint and the stored time do not match"
        }
        SUN_ERR_SUNCTX_CORRUPT => "SUNContext is NULL or corrupt",
        SUN_ERR_UNREACHABLE => {
            "reached code that should be unreachable: open an issue at: \
             https://github.com/LLNL/sundials"
        }
        SUN_ERR_UNKNOWN => "unknown error occurred: open an issue at: \
             https://github.com/LLNL/sundials",
        _ => "unknown error",
    }
}
