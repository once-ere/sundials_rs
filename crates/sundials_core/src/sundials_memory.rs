/* -----------------------------------------------------------------
 * Translated from src/sundials/sundials_memory.c and
 * include/sundials/sundials_memory.h (SUNDIALS 7.7.0).
 * SUNDIALS memory helper.
 *
 * The C generic SUNMemoryHelper (ops-table over void* content)
 * becomes an enum over the concrete helper implementations; only the
 * SYSTEM (host) helper is in scope for the pure-Rust serial build
 * (sundials_system_memory.rs). The ops-table plumbing functions
 * (SUNMemoryHelper_NewEmpty, SUNMemoryHelper_CopyOps) are subsumed
 * by the enum and have no translation.
 *
 * struct SUNMemory_ owns its data: the C `void* ptr` becomes
 * `Option<Vec<u8>>` (None = NULL pointer). Consequences (documented
 * ownership adaptations, function names / argument order / return
 * codes unchanged):
 *   - SUNMemoryHelper_Alias deep-copies the bytes instead of
 *     aliasing the pointer (safe Rust cannot share a mutable
 *     buffer); `own` is still SUNFALSE so Dealloc skips the
 *     allocation statistics exactly as in C.
 *   - SUNMemoryHelper_Wrap takes ownership of a caller-provided
 *     Vec<u8> instead of borrowing a raw pointer.
 *   - The C `void* queue` parameters (device execution queues, out
 *     of scope) become `Option<&mut dyn Any>`; the SYSTEM helper
 *     never dereferences them, so SUNMemoryHelper_SetDefaultQueue
 *     stores nothing.
 * SUNDIALS_ENABLE_PROFILING is compiled out (getSUNProfiler and the
 * SUNDIALS_MARK_FUNCTION_* macros are no-ops in this build).
 * -----------------------------------------------------------------*/
use crate::sundials_context::SUNContext;
use crate::sundials_errors::*;
use crate::sundials_system_memory::{self, SUNMemoryHelper_Content_Sys};
use crate::sundials_types::*;

/// enum SUNMemoryType (sundials_memory.h)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SUNMemoryType {
    /// pageable memory accessible on the host
    SUNMEMTYPE_HOST,
    /// page-locked memory accessible on the host
    SUNMEMTYPE_PINNED,
    /// memory accessible from the device
    SUNMEMTYPE_DEVICE,
    /// memory accessible from the host or device
    SUNMEMTYPE_UVM,
}
pub use SUNMemoryType::*;

/// C `void* queue` argument (device queues are out of scope; the host
/// paths never dereference it).
pub type SUNMemQueue<'a> = Option<&'a mut dyn std::any::Any>;

/// struct SUNMemory_ (sundials_memory.h): an abstraction of a
/// contiguous block of memory that tracks its type and ownership.
#[derive(Clone, Debug)]
pub struct SUNMemory {
    /// C `void* ptr`; None = NULL.
    pub ptr: Option<Vec<u8>>,
    /// C `SUNMemoryType type` (`type` is reserved in Rust).
    pub type_: SUNMemoryType,
    pub own: sunbooleantype,
    pub bytes: usize,
    pub stride: usize,
}

/// SUNMemoryNewEmpty: creates a new SUNMemory object with a NULL ptr.
/// (C leaves ptr/type/own uninitialized and sets only bytes/stride;
/// Rust must initialize every field: NULL ptr, HOST type, not owned.)
pub fn SUNMemoryNewEmpty(_sunctx: &SUNContext) -> SUNMemory {
    SUNMemory { ptr: None, type_: SUNMEMTYPE_HOST, own: SUNFALSE, bytes: 0, stride: 1 }
}

/// Generic SUNMemoryHelper: holds the implementation that can
/// allocate, deallocate, and copy SUNMemory.
#[derive(Clone, Debug)]
pub enum SUNMemoryHelper {
    /// SUNMemoryHelper_Sys (sunmemory_system.h) with its content.
    Sys(SUNMemoryHelper_Content_Sys),
}

/// SUNMemoryHelper_ImplementsRequiredOps: every enum variant provides
/// alloc, dealloc, and copy, so this is always SUNTRUE.
pub fn SUNMemoryHelper_ImplementsRequiredOps(_helper: &SUNMemoryHelper) -> sunbooleantype {
    SUNTRUE
}

/// SUNMemoryHelper_Alias: creates a new SUNMemory object which points
/// to the same data as another SUNMemory object; it does not own the
/// data, so Dealloc will not free it. Ownership adaptation: the bytes
/// are deep-copied (see module header). As in C, `bytes` stays 0 and
/// `stride` stays 1 (SUNMemoryNewEmpty defaults).
pub fn SUNMemoryHelper_Alias(_helper: &SUNMemoryHelper, mem: &SUNMemory) -> SUNMemory {
    let mut alias = SUNMemoryNewEmpty(&SUNContext::default());

    alias.ptr = mem.ptr.clone();
    alias.type_ = mem.type_;
    alias.own = SUNFALSE;

    alias
}

/// SUNMemoryHelper_Wrap: creates a new SUNMemory object wrapping user
/// provided data; it does not own the data (Dealloc will not count
/// it). The C mem_type validity check is subsumed by the exhaustive
/// Rust enum.
pub fn SUNMemoryHelper_Wrap(
    _helper: &SUNMemoryHelper,
    ptr: Vec<u8>,
    mem_type: SUNMemoryType,
) -> SUNMemory {
    let mut mem = SUNMemoryNewEmpty(&SUNContext::default());

    mem.ptr = Some(ptr);
    mem.own = SUNFALSE;
    mem.type_ = mem_type;

    mem
}

/// SUNMemoryHelper_GetAllocStats
pub fn SUNMemoryHelper_GetAllocStats(
    helper: &SUNMemoryHelper,
    mem_type: SUNMemoryType,
    num_allocations: &mut u64,
    num_deallocations: &mut u64,
    bytes_allocated: &mut usize,
    bytes_high_watermark: &mut usize,
) -> SUNErrCode {
    match helper {
        SUNMemoryHelper::Sys(_) => sundials_system_memory::SUNMemoryHelper_GetAllocStats_Sys(
            helper,
            mem_type,
            num_allocations,
            num_deallocations,
            bytes_allocated,
            bytes_high_watermark,
        ),
    }
}

/// SUNMemoryHelper_Alloc
pub fn SUNMemoryHelper_Alloc(
    helper: &mut SUNMemoryHelper,
    memptr: &mut Option<SUNMemory>,
    mem_size: usize,
    mem_type: SUNMemoryType,
    queue: SUNMemQueue<'_>,
) -> SUNErrCode {
    match helper {
        SUNMemoryHelper::Sys(_) => sundials_system_memory::SUNMemoryHelper_Alloc_Sys(
            helper, memptr, mem_size, mem_type, queue,
        ),
    }
}

/// SUNMemoryHelper_AllocStrided
pub fn SUNMemoryHelper_AllocStrided(
    helper: &mut SUNMemoryHelper,
    memptr: &mut Option<SUNMemory>,
    mem_size: usize,
    stride: usize,
    mem_type: SUNMemoryType,
    queue: SUNMemQueue<'_>,
) -> SUNErrCode {
    match helper {
        SUNMemoryHelper::Sys(_) => sundials_system_memory::SUNMemoryHelper_AllocStrided_Sys(
            helper, memptr, mem_size, stride, mem_type, queue,
        ),
    }
}

/// SUNMemoryHelper_Dealloc
pub fn SUNMemoryHelper_Dealloc(
    helper: &mut SUNMemoryHelper,
    mem: Option<SUNMemory>,
    queue: SUNMemQueue<'_>,
) -> SUNErrCode {
    if mem.is_none() {
        return SUN_SUCCESS;
    }
    match helper {
        SUNMemoryHelper::Sys(_) => {
            sundials_system_memory::SUNMemoryHelper_Dealloc_Sys(helper, mem, queue)
        }
    }
}

/// SUNMemoryHelper_Copy
pub fn SUNMemoryHelper_Copy(
    helper: &SUNMemoryHelper,
    dst: &mut SUNMemory,
    src: &SUNMemory,
    memory_size: usize,
    queue: SUNMemQueue<'_>,
) -> SUNErrCode {
    match helper {
        SUNMemoryHelper::Sys(_) => sundials_system_memory::SUNMemoryHelper_Copy_Sys(
            helper,
            dst,
            src,
            memory_size,
            queue,
        ),
    }
}

/// SUNMemoryHelper_CopyAsync: the SYSTEM helper provides no copyasync
/// op, so this falls back to SUNMemoryHelper_Copy (C default path).
pub fn SUNMemoryHelper_CopyAsync(
    helper: &SUNMemoryHelper,
    dst: &mut SUNMemory,
    src: &SUNMemory,
    memory_size: usize,
    queue: SUNMemQueue<'_>,
) -> SUNErrCode {
    match helper {
        SUNMemoryHelper::Sys(_) => SUNMemoryHelper_Copy(helper, dst, src, memory_size, queue),
    }
}

/// SUNMemoryHelper_Destroy: frees the SUNMemoryHelper (RAII).
pub fn SUNMemoryHelper_Destroy(helper: SUNMemoryHelper) -> SUNErrCode {
    match helper {
        SUNMemoryHelper::Sys(_) => sundials_system_memory::SUNMemoryHelper_Destroy_Sys(helper),
    }
}

/// SUNMemoryHelper_Clone
pub fn SUNMemoryHelper_Clone(helper: &SUNMemoryHelper) -> SUNMemoryHelper {
    match helper {
        SUNMemoryHelper::Sys(_) => sundials_system_memory::SUNMemoryHelper_Clone_Sys(helper),
    }
}

/// SUNMemoryHelper_SetDefaultQueue: sets the default queue to use for
/// memory helper operations. The SYSTEM helper never uses a queue, so
/// (matching the observable C behavior) nothing is stored.
pub fn SUNMemoryHelper_SetDefaultQueue(
    _helper: &mut SUNMemoryHelper,
    _queue: SUNMemQueue<'_>,
) -> SUNErrCode {
    SUN_SUCCESS
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sundials_context::SUNContext_Create;
    use crate::sundials_system_memory::SUNMemoryHelper_Sys;

    #[test]
    fn new_empty_and_wrap_and_alias() {
        let ctx = SUNContext_Create();
        let mem = SUNMemoryNewEmpty(&ctx);
        assert!(mem.ptr.is_none());
        assert_eq!(mem.bytes, 0);
        assert_eq!(mem.stride, 1);

        let helper = SUNMemoryHelper_Sys(&ctx);
        assert!(SUNMemoryHelper_ImplementsRequiredOps(&helper));

        let wrapped = SUNMemoryHelper_Wrap(&helper, vec![1u8, 2, 3], SUNMEMTYPE_HOST);
        assert_eq!(wrapped.ptr.as_deref(), Some(&[1u8, 2, 3][..]));
        assert!(!wrapped.own);

        let alias = SUNMemoryHelper_Alias(&helper, &wrapped);
        assert_eq!(alias.ptr.as_deref(), Some(&[1u8, 2, 3][..]));
        assert_eq!(alias.type_, SUNMEMTYPE_HOST);
        assert!(!alias.own);
        assert_eq!(alias.bytes, 0); /* as in C: NewEmpty default kept */
    }

    #[test]
    fn generic_alloc_copy_dealloc_and_stats() {
        let ctx = SUNContext_Create();
        let mut helper = SUNMemoryHelper_Sys(&ctx);

        let mut src = None;
        let mut dst = None;
        assert_eq!(
            SUNMemoryHelper_Alloc(&mut helper, &mut src, 8, SUNMEMTYPE_HOST, None),
            SUN_SUCCESS
        );
        assert_eq!(
            SUNMemoryHelper_AllocStrided(&mut helper, &mut dst, 8, 2, SUNMEMTYPE_HOST, None),
            SUN_SUCCESS
        );
        let mut src = src.unwrap();
        let mut dst = dst.unwrap();
        assert_eq!(dst.stride, 2);
        src.ptr.as_mut().unwrap().copy_from_slice(&[9, 8, 7, 6, 5, 4, 3, 2]);

        assert_eq!(SUNMemoryHelper_Copy(&helper, &mut dst, &src, 8, None), SUN_SUCCESS);
        assert_eq!(dst.ptr.as_deref(), src.ptr.as_deref());

        /* CopyAsync falls back to Copy */
        let mut dst2 = SUNMemoryHelper_Wrap(&helper, vec![0u8; 8], SUNMEMTYPE_HOST);
        assert_eq!(SUNMemoryHelper_CopyAsync(&helper, &mut dst2, &src, 8, None), SUN_SUCCESS);
        assert_eq!(dst2.ptr.as_deref(), src.ptr.as_deref());

        /* device types are rejected like the C SYSTEM helper */
        let mut dev = None;
        assert_eq!(
            SUNMemoryHelper_Alloc(&mut helper, &mut dev, 8, SUNMEMTYPE_DEVICE, None),
            SUN_ERR_ARG_INCOMPATIBLE
        );

        let (mut na, mut nd, mut ba, mut hw) = (0u64, 0u64, 0usize, 0usize);
        assert_eq!(
            SUNMemoryHelper_GetAllocStats(&helper, SUNMEMTYPE_HOST, &mut na, &mut nd,
                                          &mut ba, &mut hw),
            SUN_SUCCESS
        );
        assert_eq!((na, nd, ba, hw), (2, 0, 16, 16));

        assert_eq!(SUNMemoryHelper_Dealloc(&mut helper, Some(src), None), SUN_SUCCESS);
        assert_eq!(SUNMemoryHelper_Dealloc(&mut helper, None, None), SUN_SUCCESS);
        /* dst2 does not own its data: no stats update */
        assert_eq!(SUNMemoryHelper_Dealloc(&mut helper, Some(dst2), None), SUN_SUCCESS);

        assert_eq!(
            SUNMemoryHelper_GetAllocStats(&helper, SUNMEMTYPE_HOST, &mut na, &mut nd,
                                          &mut ba, &mut hw),
            SUN_SUCCESS
        );
        assert_eq!((na, nd, ba, hw), (2, 1, 8, 16));

        /* clone starts with fresh stats; destroy is RAII */
        let clone = SUNMemoryHelper_Clone(&helper);
        assert_eq!(
            SUNMemoryHelper_GetAllocStats(&clone, SUNMEMTYPE_HOST, &mut na, &mut nd,
                                          &mut ba, &mut hw),
            SUN_SUCCESS
        );
        assert_eq!((na, nd, ba, hw), (0, 0, 0, 0));
        assert_eq!(SUNMemoryHelper_Destroy(clone), SUN_SUCCESS);
        assert_eq!(SUNMemoryHelper_SetDefaultQueue(&mut helper, None), SUN_SUCCESS);
        assert_eq!(SUNMemoryHelper_Destroy(helper), SUN_SUCCESS);
    }
}
