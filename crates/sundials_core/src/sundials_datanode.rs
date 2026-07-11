/* -----------------------------------------------------------------
 * Translation of sundials-7.7.0/src/sundials/sundials_datanode.c
 * (+ src/sundials/sundials_datanode.h).
 *
 * The SUNDataNode class is a hierarchical object which could be used
 * to build something like a JSON tree. Nodes can be lists, objects,
 * or leaves. Object nodes hold named children, list nodes hold
 * anonymous children, leaf nodes hold values. It is primarily the
 * backbone for checkpointing states in adjoint sensitivity analysis.
 *
 * Modeling (see ARCHITECTURE.md):
 *  - The ops table becomes enum dispatch over SUNDataNodeContent
 *    (only the InMem implementation exists in C 7.7.0), matching the
 *    other core base classes. SUNDataNode_CreateEmpty is C object-
 *    model plumbing subsumed by the constructors.
 *  - The C tree hands out ALIASED child pointers while the tree keeps
 *    ownership. In Rust the tree owns its children; GetChild /
 *    GetNamedChild lend `&mut SUNDataNode` borrows (C's accessor
 *    gives full access through the alias), and RemoveChild /
 *    RemoveNamedChild move the child out to the caller (C: the tree
 *    forgets the pointer, caller must destroy it) — exactly the
 *    SUNStlVector_Erase / SUNHashMap_Remove semantics.
 *  - SUNDataNode_Destroy is ownership drop (children, hashmaps and
 *    leaf buffers drop recursively; C walks the tree with destroy
 *    callbacks).
 *  - The constructors' io_mode switch default (SUN_ERR_ARG_OUTOFRANGE)
 *    is unreachable: the Rust SUNDataIOMode enum has only INMEM.
 * -----------------------------------------------------------------*/

use crate::nvector_serial::NVector;
use crate::sundials_context::SUNContext;
use crate::sundials_errors::{SUNErrCode, SUN_SUCCESS};
use crate::sundials_memory::{SUNMemoryHelper, SUNMemoryType};
use crate::sundials_types::{sunbooleantype, sunrealtype, SUNDataIOMode, SUNDATAIOMODE_INMEM};
use crate::sundatanode_inmem::{self, SUNDataNode_InMemContent};

/// C `typedef int64_t sundataindex`
pub type sundataindex = i64;

/// enum SUNDataNodeType
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SUNDataNodeType {
    SUNDATANODE_LEAF,
    SUNDATANODE_LIST,
    SUNDATANODE_OBJECT,
}
pub use SUNDataNodeType::{SUNDATANODE_LEAF, SUNDATANODE_LIST, SUNDATANODE_OBJECT};

/// The implementation content (C `void* content` + ops table): enum
/// dispatch, one variant per implementation.
pub enum SUNDataNodeContent {
    InMem(SUNDataNode_InMemContent),
}

/// struct SUNDataNode_ (sundials_datanode.h)
pub struct SUNDataNode {
    pub dtype: SUNDataNodeType,
    pub content: SUNDataNodeContent,
}

/* -----------------------------------------------------------------
 * sundials_datanode.c
 * ----------------------------------------------------------------- */

pub fn SUNDataNode_CreateLeaf(
    io_mode: SUNDataIOMode,
    mem_helper: &SUNMemoryHelper,
    sunctx: &SUNContext,
    node_out: &mut Option<SUNDataNode>,
) -> SUNErrCode {
    match io_mode {
        SUNDATAIOMODE_INMEM => {
            sundatanode_inmem::SUNDataNode_CreateLeaf_InMem(mem_helper, sunctx, node_out)
        }
    }
}

pub fn SUNDataNode_CreateList(
    io_mode: SUNDataIOMode,
    num_elements: sundataindex,
    sunctx: &SUNContext,
    node_out: &mut Option<SUNDataNode>,
) -> SUNErrCode {
    match io_mode {
        SUNDATAIOMODE_INMEM => {
            sundatanode_inmem::SUNDataNode_CreateList_InMem(num_elements, sunctx, node_out)
        }
    }
}

pub fn SUNDataNode_CreateObject(
    io_mode: SUNDataIOMode,
    num_elements: sundataindex,
    sunctx: &SUNContext,
    node_out: &mut Option<SUNDataNode>,
) -> SUNErrCode {
    match io_mode {
        SUNDATAIOMODE_INMEM => {
            sundatanode_inmem::SUNDataNode_CreateObject_InMem(num_elements, sunctx, node_out)
        }
    }
}

pub fn SUNDataNode_IsLeaf(node: &SUNDataNode, yes_or_no: &mut sunbooleantype) -> SUNErrCode {
    match &node.content {
        SUNDataNodeContent::InMem(_) => {
            sundatanode_inmem::SUNDataNode_IsLeaf_InMem(node, yes_or_no)
        }
    }
}

pub fn SUNDataNode_IsList(node: &SUNDataNode, yes_or_no: &mut sunbooleantype) -> SUNErrCode {
    match &node.content {
        SUNDataNodeContent::InMem(_) => {
            sundatanode_inmem::SUNDataNode_IsList_InMem(node, yes_or_no)
        }
    }
}

pub fn SUNDataNode_HasChildren(node: &SUNDataNode, yes_or_no: &mut sunbooleantype) -> SUNErrCode {
    match &node.content {
        SUNDataNodeContent::InMem(_) => {
            sundatanode_inmem::SUNDataNode_HasChildren_InMem(node, yes_or_no)
        }
    }
}

pub fn SUNDataNode_AddChild(node: &mut SUNDataNode, child_node: SUNDataNode) -> SUNErrCode {
    match &node.content {
        SUNDataNodeContent::InMem(_) => {
            sundatanode_inmem::SUNDataNode_AddChild_InMem(node, child_node)
        }
    }
}

pub fn SUNDataNode_AddNamedChild(
    node: &mut SUNDataNode,
    name: &str,
    child_node: SUNDataNode,
) -> SUNErrCode {
    match &node.content {
        SUNDataNodeContent::InMem(_) => {
            sundatanode_inmem::SUNDataNode_AddNamedChild_InMem(node, name, child_node)
        }
    }
}

pub fn SUNDataNode_GetChild<'a>(
    node: &'a mut SUNDataNode,
    index: sundataindex,
    child_node: &mut Option<&'a mut SUNDataNode>,
) -> SUNErrCode {
    match &node.content {
        SUNDataNodeContent::InMem(_) => {
            sundatanode_inmem::SUNDataNode_GetChild_InMem(node, index, child_node)
        }
    }
}

pub fn SUNDataNode_GetNamedChild<'a>(
    node: &'a mut SUNDataNode,
    name: &str,
    child_node: &mut Option<&'a mut SUNDataNode>,
) -> SUNErrCode {
    match &node.content {
        SUNDataNodeContent::InMem(_) => {
            sundatanode_inmem::SUNDataNode_GetNamedChild_InMem(node, name, child_node)
        }
    }
}

pub fn SUNDataNode_RemoveChild(
    node: &mut SUNDataNode,
    index: sundataindex,
    child_node: &mut Option<SUNDataNode>,
) -> SUNErrCode {
    match &node.content {
        SUNDataNodeContent::InMem(_) => {
            sundatanode_inmem::SUNDataNode_RemoveChild_InMem(node, index, child_node)
        }
    }
}

pub fn SUNDataNode_RemoveNamedChild(
    node: &mut SUNDataNode,
    name: &str,
    child_node: &mut Option<SUNDataNode>,
) -> SUNErrCode {
    match &node.content {
        SUNDataNodeContent::InMem(_) => {
            sundatanode_inmem::SUNDataNode_RemoveNamedChild_InMem(node, name, child_node)
        }
    }
}

pub fn SUNDataNode_GetData<'a>(
    node: &'a SUNDataNode,
    data: &mut Option<&'a [u8]>,
    data_stride: &mut usize,
    data_bytes: &mut usize,
) -> SUNErrCode {
    match &node.content {
        SUNDataNodeContent::InMem(_) => {
            sundatanode_inmem::SUNDataNode_GetData_InMem(node, data, data_stride, data_bytes)
        }
    }
}

/// C signature: (const SUNDataNode, N_Vector v, sunrealtype* t). The
/// node is `&mut` here because the non-host staging path allocates
/// through the node's memory helper (alloc-stats update).
pub fn SUNDataNode_GetDataNvector(
    node: &mut SUNDataNode,
    v: &mut NVector,
    t: &mut sunrealtype,
) -> SUNErrCode {
    match &node.content {
        SUNDataNodeContent::InMem(_) => {
            sundatanode_inmem::SUNDataNode_GetDataNvector_InMem(node, v, t)
        }
    }
}

pub fn SUNDataNode_SetData(
    node: &mut SUNDataNode,
    src_mem_type: SUNMemoryType,
    node_mem_type: SUNMemoryType,
    data: Vec<u8>,
    data_stride: usize,
    data_bytes: usize,
) -> SUNErrCode {
    match &node.content {
        SUNDataNodeContent::InMem(_) => sundatanode_inmem::SUNDataNode_SetData_InMem(
            node,
            src_mem_type,
            node_mem_type,
            data,
            data_stride,
            data_bytes,
        ),
    }
}

pub fn SUNDataNode_SetDataNvector(node: &mut SUNDataNode, v: &NVector, t: sunrealtype) -> SUNErrCode {
    match &node.content {
        SUNDataNodeContent::InMem(_) => {
            sundatanode_inmem::SUNDataNode_SetDataNvector_InMem(node, v, t)
        }
    }
}

pub fn SUNDataNode_Destroy(node: &mut Option<SUNDataNode>) -> SUNErrCode {
    if let Some(n) = node.as_ref() {
        match &n.content {
            SUNDataNodeContent::InMem(_) => {
                return sundatanode_inmem::SUNDataNode_Destroy_InMem(node);
            }
        }
    }
    /* no destroy op: free base object only (drop) */
    *node = None;
    SUN_SUCCESS
}
