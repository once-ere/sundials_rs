/* -----------------------------------------------------------------
 * Translated from src/sunmemory/system/sundials_system_memory.c and
 * include/sunmemory/sunmemory_system.h (SUNDIALS 7.7.0).
 * SUNDIALS memory helper implementation that uses the standard
 * system memory allocators.
 *
 * malloc/free become Vec allocation/RAII (see sundials_memory.rs for
 * the ownership adaptations). The C SUNAssert checks that produce
 * error returns are kept; the ones covering C undefined behavior
 * (NULL deref / buffer overrun in memcpy) become explicit checks
 * returning SUN_ERR_ARG_CORRUPT / SUN_ERR_ARG_OUTOFRANGE. malloc
 * leaves memory uninitialized; safe Rust zero-initializes instead,
 * and a failed allocation aborts rather than returning
 * SUN_ERR_MALLOC_FAIL.
 * -----------------------------------------------------------------*/
use crate::sundials_context::SUNContext;
use crate::sundials_errors::*;
use crate::sundials_memory::*;

/// struct SUNMemoryHelper_Content_Sys_
#[derive(Clone, Debug, Default)]
pub struct SUNMemoryHelper_Content_Sys {
    pub num_allocations: u64,
    pub num_deallocations: u64,
    pub bytes_allocated: usize,
    pub bytes_high_watermark: usize,
}

/// C macro SUNHELPER_CONTENT(h): access the SYSTEM helper content.
fn SUNHELPER_CONTENT(h: &mut SUNMemoryHelper) -> &mut SUNMemoryHelper_Content_Sys {
    let SUNMemoryHelper::Sys(content) = h;
    content
}

/// SUNMemoryHelper_Sys: create the SYSTEM memory helper. The C ops
/// wiring is subsumed by the SUNMemoryHelper enum.
pub fn SUNMemoryHelper_Sys(_sunctx: &SUNContext) -> SUNMemoryHelper {
    SUNMemoryHelper::Sys(SUNMemoryHelper_Content_Sys {
        num_allocations: 0,
        num_deallocations: 0,
        bytes_allocated: 0,
        bytes_high_watermark: 0,
    })
}

/// SUNMemoryHelper_Alloc_Sys
pub fn SUNMemoryHelper_Alloc_Sys(
    helper: &mut SUNMemoryHelper,
    memptr: &mut Option<SUNMemory>,
    mem_size: usize,
    mem_type: SUNMemoryType,
    _queue: SUNMemQueue<'_>,
) -> SUNErrCode {
    if !(mem_type == SUNMEMTYPE_HOST || mem_type == SUNMEMTYPE_UVM) {
        return SUN_ERR_ARG_INCOMPATIBLE;
    }

    let mut mem = SUNMemoryNewEmpty(&SUNContext::default());

    mem.ptr = None;
    mem.own = crate::sundials_types::SUNTRUE;
    mem.type_ = mem_type;
    mem.bytes = mem_size;

    if mem_type == SUNMEMTYPE_HOST || mem_type == SUNMEMTYPE_UVM {
        mem.ptr = Some(vec![0u8; mem_size]); /* C: malloc(mem_size) */
        let content = SUNHELPER_CONTENT(helper);
        content.bytes_allocated += mem_size;
        content.num_allocations += 1;
        content.bytes_high_watermark =
            content.bytes_allocated.max(content.bytes_high_watermark); /* SUNMAX */
    }

    *memptr = Some(mem);
    SUN_SUCCESS
}

/// SUNMemoryHelper_AllocStrided_Sys
pub fn SUNMemoryHelper_AllocStrided_Sys(
    helper: &mut SUNMemoryHelper,
    memptr: &mut Option<SUNMemory>,
    mem_size: usize,
    stride: usize,
    mem_type: SUNMemoryType,
    queue: SUNMemQueue<'_>,
) -> SUNErrCode {
    /* SUNCheckCall */
    let err = SUNMemoryHelper_Alloc_Sys(helper, memptr, mem_size, mem_type, queue);
    if err != SUN_SUCCESS {
        return err;
    }

    if let Some(mem) = memptr {
        mem.stride = stride;
    }

    SUN_SUCCESS
}

/// SUNMemoryHelper_Dealloc_Sys
pub fn SUNMemoryHelper_Dealloc_Sys(
    helper: &mut SUNMemoryHelper,
    mem: Option<SUNMemory>,
    _queue: SUNMemQueue<'_>,
) -> SUNErrCode {
    let Some(mut mem) = mem else { return SUN_SUCCESS };

    if !(mem.type_ == SUNMEMTYPE_HOST || mem.type_ == SUNMEMTYPE_UVM) {
        return SUN_ERR_ARG_INCOMPATIBLE;
    }

    if mem.ptr.is_some() && mem.own {
        if mem.type_ == SUNMEMTYPE_HOST || mem.type_ == SUNMEMTYPE_UVM {
            let content = SUNHELPER_CONTENT(helper);
            content.num_deallocations += 1;
            content.bytes_allocated -= mem.bytes;
            mem.ptr = None; /* free(mem->ptr) */
        }
    }

    /* free(mem): RAII drop */
    SUN_SUCCESS
}

/// SUNMemoryHelper_Copy_Sys: memcpy(dst->ptr, src->ptr, memory_size).
/// The NULL-deref / overrun cases (C undefined behavior) return
/// SUN_ERR_ARG_CORRUPT / SUN_ERR_ARG_OUTOFRANGE instead.
pub fn SUNMemoryHelper_Copy_Sys(
    _helper: &SUNMemoryHelper,
    dst: &mut SUNMemory,
    src: &SUNMemory,
    memory_size: usize,
    _queue: SUNMemQueue<'_>,
) -> SUNErrCode {
    if !(src.type_ == SUNMEMTYPE_HOST || src.type_ == SUNMEMTYPE_UVM) {
        return SUN_ERR_ARG_INCOMPATIBLE;
    }
    if !(dst.type_ == SUNMEMTYPE_HOST || dst.type_ == SUNMEMTYPE_UVM) {
        return SUN_ERR_ARG_INCOMPATIBLE;
    }
    match (dst.ptr.as_mut(), src.ptr.as_ref()) {
        (Some(d), Some(s)) => {
            if memory_size > d.len() || memory_size > s.len() {
                return SUN_ERR_ARG_OUTOFRANGE;
            }
            d[..memory_size].copy_from_slice(&s[..memory_size]);
        }
        _ => return SUN_ERR_ARG_CORRUPT,
    }
    SUN_SUCCESS
}

/// SUNMemoryHelper_GetAllocStats_Sys
pub fn SUNMemoryHelper_GetAllocStats_Sys(
    helper: &SUNMemoryHelper,
    mem_type: SUNMemoryType,
    num_allocations: &mut u64,
    num_deallocations: &mut u64,
    bytes_allocated: &mut usize,
    bytes_high_watermark: &mut usize,
) -> SUNErrCode {
    if !(mem_type == SUNMEMTYPE_HOST || mem_type == SUNMEMTYPE_UVM) {
        return SUN_ERR_ARG_INCOMPATIBLE;
    }
    let SUNMemoryHelper::Sys(content) = helper;
    *num_allocations = content.num_allocations;
    *num_deallocations = content.num_deallocations;
    *bytes_allocated = content.bytes_allocated;
    *bytes_high_watermark = content.bytes_high_watermark;
    SUN_SUCCESS
}

/// SUNMemoryHelper_Clone_Sys: a fresh SYSTEM helper (as in C, the
/// allocation statistics are not copied).
pub fn SUNMemoryHelper_Clone_Sys(_helper: &SUNMemoryHelper) -> SUNMemoryHelper {
    SUNMemoryHelper_Sys(&SUNContext::default())
}

/// SUNMemoryHelper_Destroy_Sys: frees content, ops, and helper (RAII).
pub fn SUNMemoryHelper_Destroy_Sys(helper: SUNMemoryHelper) -> SUNErrCode {
    drop(helper);
    SUN_SUCCESS
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sundials_context::SUNContext_Create;

    #[test]
    fn sys_alloc_copy_clone() {
        let ctx = SUNContext_Create();
        let mut helper = SUNMemoryHelper_Sys(&ctx);

        let mut a = None;
        assert_eq!(
            SUNMemoryHelper_Alloc_Sys(&mut helper, &mut a, 4, SUNMEMTYPE_HOST, None),
            SUN_SUCCESS
        );
        let mut a = a.unwrap();
        assert!(a.own);
        assert_eq!(a.bytes, 4);
        assert_eq!(a.type_, SUNMEMTYPE_HOST);
        a.ptr.as_mut().unwrap().copy_from_slice(&[10, 20, 30, 40]);

        let mut b = None;
        assert_eq!(
            SUNMemoryHelper_Alloc_Sys(&mut helper, &mut b, 4, SUNMEMTYPE_HOST, None),
            SUN_SUCCESS
        );
        let mut b = b.unwrap();
        assert_eq!(SUNMemoryHelper_Copy_Sys(&helper, &mut b, &a, 4, None), SUN_SUCCESS);
        assert_eq!(b.ptr.as_deref(), Some(&[10u8, 20, 30, 40][..]));

        /* overrun and NULL-src cases are rejected instead of C UB */
        assert_eq!(
            SUNMemoryHelper_Copy_Sys(&helper, &mut b, &a, 5, None),
            SUN_ERR_ARG_OUTOFRANGE
        );
        let empty = SUNMemoryNewEmpty(&ctx);
        assert_eq!(
            SUNMemoryHelper_Copy_Sys(&helper, &mut b, &empty, 1, None),
            SUN_ERR_ARG_CORRUPT
        );

        /* stats: 2 allocations of 4 bytes */
        let (mut na, mut nd, mut ba, mut hw) = (0u64, 0u64, 0usize, 0usize);
        assert_eq!(
            SUNMemoryHelper_GetAllocStats_Sys(&helper, SUNMEMTYPE_HOST, &mut na, &mut nd,
                                              &mut ba, &mut hw),
            SUN_SUCCESS
        );
        assert_eq!((na, nd, ba, hw), (2, 0, 8, 8));
        assert_eq!(
            SUNMemoryHelper_GetAllocStats_Sys(&helper, SUNMEMTYPE_DEVICE, &mut na, &mut nd,
                                              &mut ba, &mut hw),
            SUN_ERR_ARG_INCOMPATIBLE
        );

        assert_eq!(SUNMemoryHelper_Dealloc_Sys(&mut helper, Some(a), None), SUN_SUCCESS);
        assert_eq!(SUNMemoryHelper_Dealloc_Sys(&mut helper, Some(b), None), SUN_SUCCESS);
        assert_eq!(
            SUNMemoryHelper_GetAllocStats_Sys(&helper, SUNMEMTYPE_HOST, &mut na, &mut nd,
                                              &mut ba, &mut hw),
            SUN_SUCCESS
        );
        assert_eq!((na, nd, ba, hw), (2, 2, 0, 8));

        /* pinned/device allocations are incompatible with this helper */
        let mut p = None;
        assert_eq!(
            SUNMemoryHelper_Alloc_Sys(&mut helper, &mut p, 4, SUNMEMTYPE_PINNED, None),
            SUN_ERR_ARG_INCOMPATIBLE
        );
        assert!(p.is_none());

        let clone = SUNMemoryHelper_Clone_Sys(&helper);
        assert_eq!(SUNMemoryHelper_Destroy_Sys(clone), SUN_SUCCESS);
        assert_eq!(SUNMemoryHelper_Destroy_Sys(helper), SUN_SUCCESS);
    }
}
