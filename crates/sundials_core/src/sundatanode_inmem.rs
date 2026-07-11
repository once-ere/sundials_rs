/* -----------------------------------------------------------------
 * Translation of sundials-7.7.0/src/sundials/sundatanode/
 * sundatanode_inmem.c (+ sundatanode_inmem.h).
 *
 * In-memory SUNDataNode implementation. Modeling notes:
 *  - SUNStlVector_SUNDataNode (stl/sunstl_vector.h) becomes
 *    Vec<SUNDataNode>; SUNStlVector_Erase does NOT destroy the
 *    removed element (it shifts and drops the slot), so it maps to
 *    Vec::remove exactly, and the C free callbacks
 *    (sunDataNode_FreeKeyValue_InMem / sunDataNode_FreeValue_InMem)
 *    are subsumed by Drop.
 *  - The `parent` back-pointer is write-only in the entire C tree
 *    (set by AddChild/AddNamedChild, cleared by Remove*, never read),
 *    so it is not representable-and-not-needed: dropped.
 *  - `name` (C: aliases the caller's string) becomes an owned String.
 *  - `mem_helper` (C: aliases the checkpoint scheme's helper) is an
 *    owned clone per node; only observable difference is per-helper
 *    allocation statistics, which nothing in-tree reads.
 *  - Leaf data is a SUNMemory byte buffer laid out exactly as in C:
 *    data_ptr[0] = t (8 bytes), then the N_VBufPack'd vector.
 *  - GetChild on a node whose anon_children is NULL is undefined
 *    behavior in C (SUNStlVector_At(NULL, ..)); here it returns
 *    SUN_ERR_DATANODE_NODENOTFOUND.
 * -----------------------------------------------------------------*/

use crate::nvector_serial::{N_VBufPack, N_VBufSize, N_VBufUnpack, NVector};
use crate::sundials_context::SUNContext;
use crate::sundials_datanode::{
    sundataindex, SUNDataNode, SUNDataNodeContent, SUNDATANODE_LEAF, SUNDATANODE_LIST,
    SUNDATANODE_OBJECT,
};
use crate::sundials_errors::{
    SUNErrCode, SUN_ERR_ARG_INCOMPATIBLE, SUN_ERR_ARG_WRONGTYPE, SUN_ERR_DATANODE_NODENOTFOUND,
    SUN_ERR_OP_FAIL, SUN_SUCCESS,
};
use crate::sundials_hashmap::{
    SUNHashMap, SUNHashMap_New, SUNHashMap_Insert, SUNHashMap_Remove,
};
use crate::sundials_memory::{
    SUNMemory, SUNMemoryHelper, SUNMemoryHelper_Alloc, SUNMemoryHelper_AllocStrided,
    SUNMemoryHelper_Copy, SUNMemoryHelper_Dealloc, SUNMemoryHelper_Wrap, SUNMemoryType,
    SUNMEMTYPE_HOST,
};
use crate::sundials_types::{sunbooleantype, sunrealtype};

/// struct SUNDataNode_InMemContent_ (sundatanode_inmem.h)
pub struct SUNDataNode_InMemContent {
    /* (parent back-reference dropped — write-only in C, see header) */
    /* Properties for Leaf nodes (nodes that store data) */
    pub mem_helper: Option<SUNMemoryHelper>,
    pub leaf_data: Option<SUNMemory>,
    /* Properties for Object nodes (collections of named nodes) */
    pub name: Option<String>,
    pub named_children: Option<SUNHashMap<SUNDataNode>>,
    pub num_named_children: sundataindex,
    /* Properties for List nodes (collections of anonymous nodes) */
    pub anon_children: Option<Vec<SUNDataNode>>,
}

/// C GET_CONTENT macro.
fn get_content(node: &SUNDataNode) -> &SUNDataNode_InMemContent {
    let SUNDataNodeContent::InMem(content) = &node.content;
    content
}

fn get_content_mut(node: &mut SUNDataNode) -> &mut SUNDataNode_InMemContent {
    let SUNDataNodeContent::InMem(content) = &mut node.content;
    content
}

/// C sunDataNode_CreateCommon_InMem: CreateEmpty + install InMem ops
/// (subsumed by the enum) + zeroed content.
fn sunDataNode_CreateCommon_InMem(_sunctx: &SUNContext) -> SUNDataNode {
    SUNDataNode {
        dtype: SUNDATANODE_LEAF, /* C CreateEmpty sets dtype = 0 */
        content: SUNDataNodeContent::InMem(SUNDataNode_InMemContent {
            mem_helper: None,
            leaf_data: None,
            name: None,
            named_children: None,
            num_named_children: 0,
            anon_children: None,
        }),
    }
}

pub fn SUNDataNode_CreateList_InMem(
    init_size: sundataindex,
    sunctx: &SUNContext,
    node_out: &mut Option<SUNDataNode>,
) -> SUNErrCode {
    let mut node = sunDataNode_CreateCommon_InMem(sunctx);

    node.dtype = SUNDATANODE_LIST;
    get_content_mut(&mut node).anon_children =
        Some(Vec::with_capacity(init_size.max(0) as usize));

    *node_out = Some(node);
    SUN_SUCCESS
}

pub fn SUNDataNode_CreateObject_InMem(
    init_size: sundataindex,
    sunctx: &SUNContext,
    node_out: &mut Option<SUNDataNode>,
) -> SUNErrCode {
    let mut node = sunDataNode_CreateCommon_InMem(sunctx);

    node.dtype = SUNDATANODE_OBJECT;

    let mut map: Option<SUNHashMap<SUNDataNode>> = None;
    let err = SUNHashMap_New(init_size, &mut map);
    if err != SUN_SUCCESS {
        return err;
    }

    get_content_mut(&mut node).named_children = map;

    *node_out = Some(node);
    SUN_SUCCESS
}

pub fn SUNDataNode_CreateLeaf_InMem(
    mem_helper: &SUNMemoryHelper,
    sunctx: &SUNContext,
    node_out: &mut Option<SUNDataNode>,
) -> SUNErrCode {
    let mut node = sunDataNode_CreateCommon_InMem(sunctx);

    node.dtype = SUNDATANODE_LEAF;
    {
        let content = get_content_mut(&mut node);
        content.mem_helper = Some(mem_helper.clone());
        content.leaf_data = None;
    }

    *node_out = Some(node);
    SUN_SUCCESS
}

pub fn SUNDataNode_IsLeaf_InMem(node: &SUNDataNode, yes_or_no: &mut sunbooleantype) -> SUNErrCode {
    *yes_or_no = node.dtype == SUNDATANODE_LEAF;
    SUN_SUCCESS
}

pub fn SUNDataNode_IsList_InMem(node: &SUNDataNode, yes_or_no: &mut sunbooleantype) -> SUNErrCode {
    *yes_or_no = node.dtype == SUNDATANODE_LIST;
    SUN_SUCCESS
}

pub fn SUNDataNode_IsObject_InMem(
    node: &SUNDataNode,
    yes_or_no: &mut sunbooleantype,
) -> SUNErrCode {
    *yes_or_no = node.dtype == SUNDATANODE_OBJECT;
    SUN_SUCCESS
}

/* len() != 0 mirrors C's SUNStlVector_SUNDataNode_Size(..) != 0 */
#[allow(clippy::len_zero)]
pub fn SUNDataNode_HasChildren_InMem(
    node: &SUNDataNode,
    yes_or_no: &mut sunbooleantype,
) -> SUNErrCode {
    let content = get_content(node);
    *yes_or_no = (content.anon_children.is_some()
        && content.anon_children.as_ref().unwrap().len() != 0)
        || content.num_named_children != 0;
    SUN_SUCCESS
}

pub fn SUNDataNode_AddChild_InMem(node: &mut SUNDataNode, child_node: SUNDataNode) -> SUNErrCode {
    if node.dtype != SUNDATANODE_LIST {
        return SUN_ERR_ARG_WRONGTYPE;
    }
    get_content_mut(node)
        .anon_children
        .as_mut()
        .unwrap()
        .push(child_node);
    /* C: IMPL_MEMBER(child_node, parent) = self — parent is write-only */
    SUN_SUCCESS
}

pub fn SUNDataNode_AddNamedChild_InMem(
    node: &mut SUNDataNode,
    name: &str,
    mut child_node: SUNDataNode,
) -> SUNErrCode {
    if node.dtype != SUNDATANODE_OBJECT {
        return SUN_ERR_ARG_WRONGTYPE;
    }

    get_content_mut(&mut child_node).name = Some(name.to_string());
    let map = get_content_mut(node).named_children.as_mut().unwrap();
    if SUNHashMap_Insert(map, name, child_node) != 0 {
        return SUN_ERR_OP_FAIL;
    }

    get_content_mut(node).num_named_children += 1;

    SUN_SUCCESS
}

pub fn SUNDataNode_GetChild_InMem<'a>(
    node: &'a mut SUNDataNode,
    index: sundataindex,
    child_node: &mut Option<&'a mut SUNDataNode>,
) -> SUNErrCode {
    let mut has_children = false;
    let err = SUNDataNode_HasChildren_InMem(node, &mut has_children);
    if err != SUN_SUCCESS {
        return err;
    }

    if !has_children {
        return SUN_ERR_DATANODE_NODENOTFOUND;
    }

    match get_content_mut(node).anon_children.as_mut() {
        Some(children) if index >= 0 && (index as usize) < children.len() => {
            *child_node = Some(&mut children[index as usize]);
            SUN_SUCCESS
        }
        _ => SUN_ERR_DATANODE_NODENOTFOUND,
    }
}

pub fn SUNDataNode_GetNamedChild_InMem<'a>(
    node: &'a mut SUNDataNode,
    name: &str,
    child_node: &mut Option<&'a mut SUNDataNode>,
) -> SUNErrCode {
    *child_node = None;

    let mut has_children = false;
    let err = SUNDataNode_HasChildren_InMem(node, &mut has_children);
    if err != SUN_SUCCESS {
        return err;
    }

    if has_children {
        match get_content_mut(node)
            .named_children
            .as_mut()
            .and_then(|map| map.get_mut(name))
        {
            Some(child) => {
                *child_node = Some(child);
                SUN_SUCCESS
            }
            None => SUN_ERR_DATANODE_NODENOTFOUND,
        }
    } else {
        SUN_ERR_DATANODE_NODENOTFOUND
    }
}

pub fn SUNDataNode_RemoveChild_InMem(
    node: &mut SUNDataNode,
    index: sundataindex,
    child_node: &mut Option<SUNDataNode>,
) -> SUNErrCode {
    let mut has_children = false;
    let err = SUNDataNode_HasChildren_InMem(node, &mut has_children);
    if err != SUN_SUCCESS {
        return err;
    }

    if !has_children {
        *child_node = None;
        return SUN_SUCCESS;
    }

    match get_content_mut(node).anon_children.as_mut() {
        Some(children) if index >= 0 && (index as usize) < children.len() => {
            /* C: alias the element, clear its parent, then StlVector_Erase
            (which does not destroy the element) — i.e. move it out */
            *child_node = Some(children.remove(index as usize));
            SUN_SUCCESS
        }
        _ => SUN_ERR_DATANODE_NODENOTFOUND,
    }
}

pub fn SUNDataNode_RemoveNamedChild_InMem(
    node: &mut SUNDataNode,
    name: &str,
    child_node: &mut Option<SUNDataNode>,
) -> SUNErrCode {
    *child_node = None;

    let mut has_children = false;
    let err = SUNDataNode_HasChildren_InMem(node, &mut has_children);
    if err != SUN_SUCCESS {
        return err;
    }

    if has_children {
        let content = get_content_mut(node);
        if SUNHashMap_Remove(content.named_children.as_mut().unwrap(), name, child_node) != 0 {
            *child_node = None;
            return SUN_ERR_DATANODE_NODENOTFOUND;
        }
        content.num_named_children -= 1;
    }

    SUN_SUCCESS
}

pub fn SUNDataNode_GetData_InMem<'a>(
    node: &'a SUNDataNode,
    data: &mut Option<&'a [u8]>,
    data_stride: &mut usize,
    data_bytes: &mut usize,
) -> SUNErrCode {
    /* C dereferences leaf_data unconditionally (UB when NULL) */
    let leaf_data = match get_content(node).leaf_data.as_ref() {
        Some(m) => m,
        None => return SUN_ERR_ARG_WRONGTYPE,
    };

    *data_stride = leaf_data.stride;
    *data_bytes = leaf_data.bytes;
    *data = leaf_data.ptr.as_deref();

    SUN_SUCCESS
}

pub fn SUNDataNode_GetDataNvector_InMem(
    node: &mut SUNDataNode,
    v: &mut NVector,
    t: &mut sunrealtype,
) -> SUNErrCode {
    /* Use the default queue for the memory helper */
    let queue = None;

    let content = get_content_mut(node);
    let leaf_data = match content.leaf_data.as_ref() {
        Some(m) => m,
        None => return SUN_ERR_ARG_WRONGTYPE,
    };

    let leaf_mem_type = leaf_data.type_;

    let mut buffer_size: i64 = 0;
    let err = N_VBufSize(v, &mut buffer_size);
    if err != SUN_SUCCESS {
        return err;
    }
    if (buffer_size as usize) + std::mem::size_of::<sunrealtype>() != leaf_data.bytes {
        return SUN_ERR_ARG_INCOMPATIBLE;
    }

    if leaf_mem_type != SUNMEMTYPE_HOST {
        /* BufUnpack assumes the data is on the host. So if the leaf has it
        elsewhere, we need to move it to the host first. */
        let mut leaf_host_data: Option<SUNMemory> = None;
        let helper = content.mem_helper.as_mut().unwrap();
        let err = SUNMemoryHelper_Alloc(
            helper,
            &mut leaf_host_data,
            leaf_data.bytes,
            SUNMEMTYPE_HOST,
            queue,
        );
        if err != SUN_SUCCESS {
            return err;
        }
        let mut leaf_host_data = leaf_host_data.unwrap();

        let err = SUNMemoryHelper_Copy(
            helper,
            &mut leaf_host_data,
            leaf_data,
            buffer_size as usize,
            None,
        );
        if err != SUN_SUCCESS {
            return err;
        }

        {
            let data_ptr = leaf_host_data.ptr.as_deref().unwrap();
            let mut b = [0u8; 8];
            b.copy_from_slice(&data_ptr[0..8]);
            *t = f64::from_ne_bytes(b);
            let err = N_VBufUnpack(v, &data_ptr[8..]);
            if err != SUN_SUCCESS {
                return err;
            }
        }

        let err = SUNMemoryHelper_Dealloc(helper, Some(leaf_host_data), None);
        if err != SUN_SUCCESS {
            return err;
        }
    } else {
        let data_ptr = leaf_data.ptr.as_deref().unwrap();
        let mut b = [0u8; 8];
        b.copy_from_slice(&data_ptr[0..8]);
        *t = f64::from_ne_bytes(b);
        let err = N_VBufUnpack(v, &data_ptr[8..]);
        if err != SUN_SUCCESS {
            return err;
        }
    }

    SUN_SUCCESS
}

pub fn SUNDataNode_SetData_InMem(
    node: &mut SUNDataNode,
    src_mem_type: SUNMemoryType,
    node_mem_type: SUNMemoryType,
    data: Vec<u8>,
    data_stride: usize,
    data_bytes: usize,
) -> SUNErrCode {
    /* Use the default queue for the memory helper */
    let queue = None;

    if node.dtype != SUNDATANODE_LEAF {
        return SUN_ERR_ARG_WRONGTYPE;
    }

    let content = get_content_mut(node);
    let helper = content.mem_helper.as_mut().unwrap();

    let data_mem_src = SUNMemoryHelper_Wrap(helper, data, src_mem_type);

    let mut data_mem_dst: Option<SUNMemory> = None;
    let err = SUNMemoryHelper_AllocStrided(
        helper,
        &mut data_mem_dst,
        data_bytes,
        data_stride,
        node_mem_type,
        queue,
    );
    if err != SUN_SUCCESS {
        return err;
    }
    let mut dst = data_mem_dst.unwrap();

    let err = SUNMemoryHelper_Copy(helper, &mut dst, &data_mem_src, data_bytes, None);
    if err != SUN_SUCCESS {
        return err;
    }

    SUNMemoryHelper_Dealloc(helper, Some(data_mem_src), None);

    content.leaf_data = Some(dst);

    SUN_SUCCESS
}

pub fn SUNDataNode_SetDataNvector_InMem(
    node: &mut SUNDataNode,
    v: &NVector,
    t: sunrealtype,
) -> SUNErrCode {
    /* Use the default queue for the memory helper */
    let queue = None;

    let leaf_mem_type = SUNMEMTYPE_HOST;

    let mut buffer_size: i64 = 0;
    let err = N_VBufSize(v, &mut buffer_size);
    if err != SUN_SUCCESS {
        return err;
    }

    /* We allocate 1 extra sunrealtype for storing t */
    let content = get_content_mut(node);
    let mut leaf_data: Option<SUNMemory> = None;
    let err = SUNMemoryHelper_AllocStrided(
        content.mem_helper.as_mut().unwrap(),
        &mut leaf_data,
        (buffer_size as usize) + std::mem::size_of::<sunrealtype>(),
        std::mem::size_of::<sunrealtype>(),
        leaf_mem_type,
        queue,
    );
    if err != SUN_SUCCESS {
        return err;
    }
    let mut leaf_data = leaf_data.unwrap();

    /* BufPack will handle any necessary copies and fill data_ptr on the host */
    {
        let data_ptr = leaf_data.ptr.as_deref_mut().unwrap();
        data_ptr[0..8].copy_from_slice(&t.to_ne_bytes());
        let err = N_VBufPack(v, &mut data_ptr[8..]);
        if err != SUN_SUCCESS {
            return err;
        }
    }

    content.leaf_data = Some(leaf_data);

    SUN_SUCCESS
}

/// C SUNDataNode_Destroy_InMem: hashmap/list/leaf-specific teardown +
/// sunDataNode_DestroyCommon_InMem — all subsumed by recursive drop.
pub fn SUNDataNode_Destroy_InMem(node: &mut Option<SUNDataNode>) -> SUNErrCode {
    *node = None;
    SUN_SUCCESS
}

#[cfg(test)]
mod tests {
    use crate::nvector_serial::{N_VConst, N_VNew_Serial};
    use crate::sundials_context::SUNContext_Create;
    use crate::sundials_datanode::*;
    use crate::sundials_errors::*;
    use crate::sundials_system_memory::SUNMemoryHelper_Sys;
    use crate::sundials_types::SUNDATAIOMODE_INMEM;

    /// Mirrors the fixed checkpoint scheme's usage: object root holding
    /// named list nodes, each holding leaf snapshots of (t, y).
    #[test]
    fn tree_insert_load_remove_roundtrip() {
        let ctx = SUNContext_Create();
        let helper = SUNMemoryHelper_Sys(&ctx);

        let mut root = None;
        assert_eq!(
            SUNDataNode_CreateObject(SUNDATAIOMODE_INMEM, 16, &ctx, &mut root),
            SUN_SUCCESS
        );
        let mut root = root.unwrap();
        assert_eq!(root.dtype, SUNDATANODE_OBJECT);

        let mut yes = true;
        assert_eq!(SUNDataNode_HasChildren(&root, &mut yes), SUN_SUCCESS);
        assert!(!yes);
        assert_eq!(SUNDataNode_IsLeaf(&root, &mut yes), SUN_SUCCESS);
        assert!(!yes);

        /* step 0: list node with two stage snapshots */
        let mut step_node = None;
        assert_eq!(
            SUNDataNode_CreateList(SUNDATAIOMODE_INMEM, 0, &ctx, &mut step_node),
            SUN_SUCCESS
        );
        let mut step_node = step_node.unwrap();
        assert_eq!(SUNDataNode_IsList(&step_node, &mut yes), SUN_SUCCESS);
        assert!(yes);

        let mut y = N_VNew_Serial(4, &ctx);
        for stage in 0..2 {
            let mut leaf = None;
            assert_eq!(
                SUNDataNode_CreateLeaf(SUNDATAIOMODE_INMEM, &helper, &ctx, &mut leaf),
                SUN_SUCCESS
            );
            let mut leaf = leaf.unwrap();
            N_VConst(1.5 + stage as f64, &mut y);
            assert_eq!(
                SUNDataNode_SetDataNvector(&mut leaf, &y, 0.25 * stage as f64),
                SUN_SUCCESS
            );
            assert_eq!(SUNDataNode_AddChild(&mut step_node, leaf), SUN_SUCCESS);
        }
        assert_eq!(SUNDataNode_AddNamedChild(&mut root, "0", step_node), SUN_SUCCESS);

        assert_eq!(SUNDataNode_HasChildren(&root, &mut yes), SUN_SUCCESS);
        assert!(yes);

        /* load stage 1 back through the aliased accessors */
        {
            let mut step = None;
            assert_eq!(SUNDataNode_GetNamedChild(&mut root, "0", &mut step), SUN_SUCCESS);
            let step = step.unwrap();
            let mut leaf = None;
            assert_eq!(SUNDataNode_GetChild(step, 1, &mut leaf), SUN_SUCCESS);
            let leaf = leaf.unwrap();
            let mut t = -1.0;
            let mut yout = N_VNew_Serial(4, &ctx);
            assert_eq!(SUNDataNode_GetDataNvector(leaf, &mut yout, &mut t), SUN_SUCCESS);
            assert_eq!(t, 0.25);
            assert_eq!(yout.data, vec![2.5; 4]);

            /* raw view: 8 bytes of t + 32 bytes of packed vector, stride 8 */
            let (mut data, mut stride, mut bytes) = (None, 0usize, 0usize);
            assert_eq!(
                SUNDataNode_GetData(leaf, &mut data, &mut stride, &mut bytes),
                SUN_SUCCESS
            );
            assert_eq!((stride, bytes), (8, 40));
            assert_eq!(data.unwrap().len(), 40);
        }

        /* removal semantics: stage 0 comes out owned, list shrinks */
        {
            let mut step = None;
            assert_eq!(SUNDataNode_GetNamedChild(&mut root, "0", &mut step), SUN_SUCCESS);
            let step = step.unwrap();
            let mut removed = None;
            assert_eq!(SUNDataNode_RemoveChild(step, 0, &mut removed), SUN_SUCCESS);
            let mut removed = removed;
            let mut t = -1.0;
            let mut yout = N_VNew_Serial(4, &ctx);
            assert_eq!(
                SUNDataNode_GetDataNvector(removed.as_mut().unwrap(), &mut yout, &mut t),
                SUN_SUCCESS
            );
            assert_eq!(t, 0.0);
            assert_eq!(yout.data, vec![1.5; 4]);
            assert_eq!(SUNDataNode_Destroy(&mut removed), SUN_SUCCESS);
            assert!(removed.is_none());

            /* now index 1 no longer exists */
            let mut gone = None;
            assert_eq!(
                SUNDataNode_RemoveChild(step, 1, &mut gone),
                SUN_ERR_DATANODE_NODENOTFOUND
            );
        }

        /* remove the whole step node by name */
        let mut step_owned = None;
        assert_eq!(
            SUNDataNode_RemoveNamedChild(&mut root, "0", &mut step_owned),
            SUN_SUCCESS
        );
        assert!(step_owned.is_some());
        assert_eq!(SUNDataNode_HasChildren(&root, &mut yes), SUN_SUCCESS);
        assert!(!yes);
        assert_eq!(
            SUNDataNode_GetNamedChild(&mut root, "0", &mut None),
            SUN_ERR_DATANODE_NODENOTFOUND
        );
    }

    #[test]
    fn wrong_type_and_missing_children_errors() {
        let ctx = SUNContext_Create();
        let helper = SUNMemoryHelper_Sys(&ctx);

        let mut leaf = None;
        assert_eq!(
            SUNDataNode_CreateLeaf(SUNDATAIOMODE_INMEM, &helper, &ctx, &mut leaf),
            SUN_SUCCESS
        );
        let mut leaf = leaf.unwrap();

        let mut list = None;
        assert_eq!(
            SUNDataNode_CreateList(SUNDATAIOMODE_INMEM, 2, &ctx, &mut list),
            SUN_SUCCESS
        );
        let mut list = list.unwrap();

        /* AddChild on a non-list / AddNamedChild on a non-object */
        let mut other = None;
        assert_eq!(
            SUNDataNode_CreateList(SUNDATAIOMODE_INMEM, 0, &ctx, &mut other),
            SUN_SUCCESS
        );
        assert_eq!(
            SUNDataNode_AddChild(&mut leaf, other.take().unwrap()),
            SUN_ERR_ARG_WRONGTYPE
        );
        assert_eq!(
            SUNDataNode_CreateList(SUNDATAIOMODE_INMEM, 0, &ctx, &mut other),
            SUN_SUCCESS
        );
        assert_eq!(
            SUNDataNode_AddNamedChild(&mut list, "x", other.take().unwrap()),
            SUN_ERR_ARG_WRONGTYPE
        );

        /* accessors on childless nodes */
        assert_eq!(
            SUNDataNode_GetChild(&mut list, 0, &mut None),
            SUN_ERR_DATANODE_NODENOTFOUND
        );
        let mut removed = Some(None).unwrap();
        assert_eq!(SUNDataNode_RemoveChild(&mut list, 0, &mut removed), SUN_SUCCESS);
        assert!(removed.is_none());
        assert_eq!(
            SUNDataNode_RemoveNamedChild(&mut list, "x", &mut removed),
            SUN_SUCCESS
        );
        assert!(removed.is_none());

        /* SetData on a non-leaf */
        assert_eq!(
            SUNDataNode_SetData(
                &mut list,
                crate::sundials_memory::SUNMEMTYPE_HOST,
                crate::sundials_memory::SUNMEMTYPE_HOST,
                vec![0u8; 8],
                8,
                8
            ),
            SUN_ERR_ARG_WRONGTYPE
        );

        /* GetDataNvector size mismatch */
        let y3 = N_VNew_Serial(3, &ctx);
        assert_eq!(SUNDataNode_SetDataNvector(&mut leaf, &y3, 1.0), SUN_SUCCESS);
        let mut y4 = N_VNew_Serial(4, &ctx);
        let mut t = 0.0;
        assert_eq!(
            SUNDataNode_GetDataNvector(&mut leaf, &mut y4, &mut t),
            SUN_ERR_ARG_INCOMPATIBLE
        );
        let mut y3b = N_VNew_Serial(3, &ctx);
        assert_eq!(SUNDataNode_GetDataNvector(&mut leaf, &mut y3b, &mut t), SUN_SUCCESS);
        assert_eq!(t, 1.0);
    }

    /// Raw SetData/GetData round-trip (no in-tree C caller, but part of
    /// the public API).
    #[test]
    fn raw_set_get_data() {
        let ctx = SUNContext_Create();
        let helper = SUNMemoryHelper_Sys(&ctx);

        let mut leaf = None;
        assert_eq!(
            SUNDataNode_CreateLeaf(SUNDATAIOMODE_INMEM, &helper, &ctx, &mut leaf),
            SUN_SUCCESS
        );
        let mut leaf = leaf.unwrap();

        let payload: Vec<u8> = (0u8..16).collect();
        assert_eq!(
            SUNDataNode_SetData(
                &mut leaf,
                crate::sundials_memory::SUNMEMTYPE_HOST,
                crate::sundials_memory::SUNMEMTYPE_HOST,
                payload.clone(),
                1,
                16
            ),
            SUN_SUCCESS
        );

        let (mut data, mut stride, mut bytes) = (None, 0usize, 0usize);
        assert_eq!(
            SUNDataNode_GetData(&leaf, &mut data, &mut stride, &mut bytes),
            SUN_SUCCESS
        );
        assert_eq!((stride, bytes), (1, 16));
        assert_eq!(data.unwrap(), &payload[..]);
    }
}
