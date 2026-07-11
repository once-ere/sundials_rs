/* -----------------------------------------------------------------
 * Translation of sundials-7.7.0/src/sunadjointcheckpointscheme/fixed/
 * sunadjointcheckpointscheme_fixed.c
 * (+ include/sunadjointcheckpointscheme/
 *  sunadjointcheckpointscheme_fixed.h).
 *
 * Fixed-interval checkpointing scheme over a SUNDataNode tree:
 * an object root node holds one named list node per step (key =
 * decimal step number), each list holds the step/stage snapshots as
 * leaves.
 *
 * Adaptation: C caches raw ALIASES into the tree
 * (current_insert_step_node / current_load_step_node) purely to skip
 * the hashmap lookup; the alias is always root[key(step_num)], where
 * step_num is recorded in step_num_of_current_insert/_load. The Rust
 * port keeps those step-number fields (same branch structure as C)
 * and re-derives the node by key lookup instead of holding a borrow.
 * One divergence where C is undefined behavior: after the !keep load
 * path removes a fully-consumed step node, C leaves
 * step_num_of_current_load pointing at the freed node (re-requesting
 * that step would use-after-free); here the key lookup simply fails
 * and LoadVector returns SUN_ERR_CHECKPOINT_NOT_FOUND.
 *
 * sunSignedToString(step_num) is format!("{}"). The SUNLogExtraDebug
 * logging (debug builds only) is compiled out in this port.
 * -----------------------------------------------------------------*/

use crate::nvector_serial::NVector;
use crate::sundials_adjointcheckpointscheme::{
    SUNAdjointCheckpointScheme, SUNAdjointCheckpointScheme_NewEmpty,
};
use crate::sundials_context::SUNContext;
use crate::sundials_datanode::{
    SUNDataNode, SUNDataNode_AddChild, SUNDataNode_AddNamedChild, SUNDataNode_CreateLeaf,
    SUNDataNode_CreateList, SUNDataNode_CreateObject, SUNDataNode_Destroy,
    SUNDataNode_GetChild, SUNDataNode_GetDataNvector, SUNDataNode_GetNamedChild,
    SUNDataNode_HasChildren, SUNDataNode_RemoveChild, SUNDataNode_RemoveNamedChild,
    SUNDataNode_SetDataNvector,
};
use crate::sundials_errors::{
    SUNErrCode, SUN_ERR_CHECKPOINT_NOT_FOUND, SUN_ERR_DATANODE_NODENOTFOUND, SUN_SUCCESS,
};
use crate::sundials_memory::SUNMemoryHelper;
use crate::sundials_types::{
    suncountertype, sunbooleantype, sunrealtype, SUNDataIOMode, SUNFALSE, SUNTRUE,
};

/// struct SUNAdjointCheckpointScheme_Fixed_Content_
pub struct SUNAdjointCheckpointScheme_Fixed_Content {
    pub backup_interval: suncountertype,
    pub interval: suncountertype,
    pub step_num_of_current_insert: suncountertype,
    pub step_num_of_current_load: suncountertype,
    pub mem_helper: SUNMemoryHelper,
    pub root_node: Option<SUNDataNode>,
    /* current_insert_step_node / current_load_step_node: re-derived
    from the step_num fields by key lookup (see module header) */
    pub io_mode: SUNDataIOMode,
    pub keep: sunbooleantype,
}

/// C GET_CONTENT / IMPL_MEMBER macros (unchecked cast in C).
fn get_content_mut(
    scheme: &mut SUNAdjointCheckpointScheme,
) -> &mut SUNAdjointCheckpointScheme_Fixed_Content {
    scheme
        .content
        .as_mut()
        .unwrap()
        .downcast_mut::<SUNAdjointCheckpointScheme_Fixed_Content>()
        .unwrap()
}

pub fn SUNAdjointCheckpointScheme_Create_Fixed(
    io_mode: SUNDataIOMode,
    mem_helper: &SUNMemoryHelper,
    interval: suncountertype,
    estimate: suncountertype,
    keep: sunbooleantype,
    sunctx: &SUNContext,
    check_scheme_ptr: &mut Option<SUNAdjointCheckpointScheme>,
) -> SUNErrCode {
    let mut check_scheme = SUNAdjointCheckpointScheme_NewEmpty(sunctx);

    check_scheme.ops.needssaving = Some(SUNAdjointCheckpointScheme_NeedsSaving_Fixed);
    check_scheme.ops.insertvector = Some(SUNAdjointCheckpointScheme_InsertVector_Fixed);
    check_scheme.ops.loadvector = Some(SUNAdjointCheckpointScheme_LoadVector_Fixed);
    check_scheme.ops.enableDense = Some(SUNAdjointCheckpointScheme_EnableDense_Fixed);
    check_scheme.ops.destroy = Some(SUNAdjointCheckpointScheme_Destroy_Fixed);

    let mut content = SUNAdjointCheckpointScheme_Fixed_Content {
        backup_interval: interval,
        interval,
        step_num_of_current_insert: -2,
        step_num_of_current_load: -2,
        mem_helper: mem_helper.clone(),
        root_node: None,
        io_mode,
        keep,
    };

    let err = SUNDataNode_CreateObject(io_mode, estimate, sunctx, &mut content.root_node);
    if err != SUN_SUCCESS {
        return err;
    }

    check_scheme.content = Some(Box::new(content));
    *check_scheme_ptr = Some(check_scheme);

    SUN_SUCCESS
}

pub fn SUNAdjointCheckpointScheme_NeedsSaving_Fixed(
    scheme: &mut SUNAdjointCheckpointScheme,
    step_num: suncountertype,
    _stage_num: suncountertype,
    _t: sunrealtype,
    yes_or_no: &mut sunbooleantype,
) -> SUNErrCode {
    let content = get_content_mut(scheme);
    if step_num % content.interval == 0 {
        *yes_or_no = SUNTRUE;
    } else {
        *yes_or_no = SUNFALSE;
    }
    SUN_SUCCESS
}

pub fn SUNAdjointCheckpointScheme_InsertVector_Fixed(
    scheme: &mut SUNAdjointCheckpointScheme,
    step_num: suncountertype,
    _stage_num: suncountertype,
    t: sunrealtype,
    y: &NVector,
) -> SUNErrCode {
    let sunctx = SUNContext::default();
    let content = get_content_mut(scheme);

    /* If this is the first state for a step, then we need to create a
    list node first to store the step and all stage solutions in.
    (C keeps a pointer to the list node for fast stage access; the
    port re-derives it from the key.) */
    let key = format!("{}", step_num);
    if step_num != content.step_num_of_current_insert {
        let mut step_data_node = None;
        let err = SUNDataNode_CreateList(content.io_mode, 0, &sunctx, &mut step_data_node);
        if err != SUN_SUCCESS {
            return err;
        }
        content.step_num_of_current_insert = step_num;

        /* Store the step node in the root node object. */
        let err = SUNDataNode_AddNamedChild(
            content.root_node.as_mut().unwrap(),
            &key,
            step_data_node.unwrap(),
        );
        if err != SUN_SUCCESS {
            return err;
        }
    }

    /* Add the state data as a leaf node in the step node's list of children. */
    let mut solution_node = None;
    let err = SUNDataNode_CreateLeaf(content.io_mode, &content.mem_helper, &sunctx, &mut solution_node);
    if err != SUN_SUCCESS {
        return err;
    }
    let mut solution_node = solution_node.unwrap();
    let err = SUNDataNode_SetDataNvector(&mut solution_node, y, t);
    if err != SUN_SUCCESS {
        return err;
    }

    let mut step_data_node = None;
    let err = SUNDataNode_GetNamedChild(content.root_node.as_mut().unwrap(), &key, &mut step_data_node);
    if err != SUN_SUCCESS {
        return err;
    }
    let err = SUNDataNode_AddChild(step_data_node.unwrap(), solution_node);
    if err != SUN_SUCCESS {
        return err;
    }

    SUN_SUCCESS
}

pub fn SUNAdjointCheckpointScheme_LoadVector_Fixed(
    scheme: &mut SUNAdjointCheckpointScheme,
    step_num: suncountertype,
    stage_num: suncountertype,
    peek: sunbooleantype,
    yout: &mut NVector,
    tout: &mut sunrealtype,
) -> SUNErrCode {
    let content = get_content_mut(scheme);
    let key = format!("{}", step_num);

    /* If we are trying to load the step solution, we need to load the
    list which holds the step and stage solutions. (C caches the node
    pointer; the port re-checks by key — same branch structure.) */
    let step_found = {
        let mut step_data_node = None;
        let err = SUNDataNode_GetNamedChild(
            content.root_node.as_mut().unwrap(),
            &key,
            &mut step_data_node,
        );
        if err == SUN_SUCCESS {
            if step_num != content.step_num_of_current_load {
                content.step_num_of_current_load = step_num;
            }
            true
        } else if err == SUN_ERR_DATANODE_NODENOTFOUND {
            false
        } else {
            return err;
        }
    };

    if !step_found {
        return SUN_ERR_CHECKPOINT_NOT_FOUND;
    }

    if content.keep || peek {
        /* aliased access: read the stage in place */
        let mut step_data_node = None;
        let err = SUNDataNode_GetNamedChild(
            content.root_node.as_mut().unwrap(),
            &key,
            &mut step_data_node,
        );
        if err != SUN_SUCCESS {
            return err;
        }
        let step_data_node = step_data_node.unwrap();

        let mut solution_node = None;
        let err = SUNDataNode_GetChild(step_data_node, stage_num, &mut solution_node);
        if err != SUN_SUCCESS && err != SUN_ERR_DATANODE_NODENOTFOUND {
            return err;
        }

        match solution_node {
            None => SUN_ERR_CHECKPOINT_NOT_FOUND,
            Some(solution_node) => {
                let err = SUNDataNode_GetDataNvector(solution_node, yout, tout);
                if err != SUN_SUCCESS {
                    return err;
                }
                SUN_SUCCESS
            }
        }
    } else {
        /* consume the stage: remove it (ownership moves out) */
        let mut solution_node: Option<SUNDataNode> = None;
        let step_now_empty;
        {
            let mut step = None;
            let err =
                SUNDataNode_GetNamedChild(content.root_node.as_mut().unwrap(), &key, &mut step);
            if err != SUN_SUCCESS {
                return err;
            }
            let step = step.unwrap();

            let mut has_children = SUNFALSE;
            let err = SUNDataNode_HasChildren(step, &mut has_children);
            if err != SUN_SUCCESS {
                return err;
            }

            if has_children {
                let err = SUNDataNode_RemoveChild(step, stage_num, &mut solution_node);
                if err != SUN_SUCCESS && err != SUN_ERR_DATANODE_NODENOTFOUND {
                    return err;
                }
            }

            /* If we just removed the last stage (so has_children==false),
            then we should remove the step too. */
            let err = SUNDataNode_HasChildren(step, &mut has_children);
            if err != SUN_SUCCESS {
                return err;
            }
            step_now_empty = !has_children;
        }

        if step_now_empty {
            let mut step_data_node = None;
            let err = SUNDataNode_RemoveNamedChild(
                content.root_node.as_mut().unwrap(),
                &key,
                &mut step_data_node,
            );
            if err != SUN_SUCCESS {
                return err;
            }
            let err = SUNDataNode_Destroy(&mut step_data_node);
            if err != SUN_SUCCESS {
                return err;
            }
        }

        match solution_node.as_mut() {
            None => SUN_ERR_CHECKPOINT_NOT_FOUND,
            Some(sol) => {
                let err = SUNDataNode_GetDataNvector(sol, yout, tout);
                if err != SUN_SUCCESS {
                    return err;
                }
                /* Cleanup the checkpoint memory */
                let err = SUNDataNode_Destroy(&mut solution_node);
                if err != SUN_SUCCESS {
                    return err;
                }
                SUN_SUCCESS
            }
        }
    }
}

/// C SUNAdjointCheckpointScheme_Destroy_Fixed: destroys the root node
/// and frees content/ops/self — all ownership drop here.
pub fn SUNAdjointCheckpointScheme_Destroy_Fixed(
    self_ptr: &mut Option<SUNAdjointCheckpointScheme>,
) -> SUNErrCode {
    *self_ptr = None;
    SUN_SUCCESS
}

pub fn SUNAdjointCheckpointScheme_EnableDense_Fixed(
    check_scheme: &mut SUNAdjointCheckpointScheme,
    on_or_off: sunbooleantype,
) -> SUNErrCode {
    let content = get_content_mut(check_scheme);
    if on_or_off {
        content.backup_interval = content.interval;
        content.interval = 1;
    } else {
        content.interval = content.backup_interval;
    }
    SUN_SUCCESS
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nvector_serial::{N_VConst, N_VNew_Serial};
    use crate::sundials_adjointcheckpointscheme::*;
    use crate::sundials_context::SUNContext_Create;
    use crate::sundials_errors::SUN_ERR_NOT_IMPLEMENTED;
    use crate::sundials_system_memory::SUNMemoryHelper_Sys;
    use crate::sundials_types::SUNDATAIOMODE_INMEM;

    fn make_scheme(interval: suncountertype, keep: bool) -> SUNAdjointCheckpointScheme {
        let ctx = SUNContext_Create();
        let helper = SUNMemoryHelper_Sys(&ctx);
        let mut scheme = None;
        assert_eq!(
            SUNAdjointCheckpointScheme_Create_Fixed(
                SUNDATAIOMODE_INMEM,
                &helper,
                interval,
                8,
                keep,
                &ctx,
                &mut scheme
            ),
            SUN_SUCCESS
        );
        scheme.unwrap()
    }

    #[test]
    fn empty_scheme_ops_not_implemented() {
        let ctx = SUNContext_Create();
        let mut scheme = SUNAdjointCheckpointScheme_NewEmpty(&ctx);
        let mut yes = false;
        assert_eq!(
            SUNAdjointCheckpointScheme_NeedsSaving(&mut scheme, 0, 0, 0.0, &mut yes),
            SUN_ERR_NOT_IMPLEMENTED
        );
        let mut slot = Some(scheme);
        assert_eq!(SUNAdjointCheckpointScheme_Destroy(&mut slot), SUN_SUCCESS);
        assert!(slot.is_none());
    }

    #[test]
    fn needs_saving_interval_and_enable_dense() {
        let mut scheme = make_scheme(3, true);
        let mut yes = false;
        for (step, expect) in [(0, true), (1, false), (2, false), (3, true), (4, false)] {
            assert_eq!(
                SUNAdjointCheckpointScheme_NeedsSaving(&mut scheme, step, 0, 0.0, &mut yes),
                SUN_SUCCESS
            );
            assert_eq!(yes, expect, "step {}", step);
        }
        /* dense mode: every step saves; disabling restores interval */
        assert_eq!(SUNAdjointCheckpointScheme_EnableDense(&mut scheme, true), SUN_SUCCESS);
        assert_eq!(
            SUNAdjointCheckpointScheme_NeedsSaving(&mut scheme, 4, 0, 0.0, &mut yes),
            SUN_SUCCESS
        );
        assert!(yes);
        assert_eq!(SUNAdjointCheckpointScheme_EnableDense(&mut scheme, false), SUN_SUCCESS);
        assert_eq!(
            SUNAdjointCheckpointScheme_NeedsSaving(&mut scheme, 4, 0, 0.0, &mut yes),
            SUN_SUCCESS
        );
        assert!(!yes);
    }

    #[test]
    fn insert_load_keep_and_peek() {
        let ctx = SUNContext_Create();
        let mut scheme = make_scheme(1, true);
        let mut y = N_VNew_Serial(3, &ctx);

        /* two steps, two stages each */
        for step in 0..2i64 {
            for stage in 0..2i64 {
                N_VConst((10 * step + stage) as f64, &mut y);
                assert_eq!(
                    SUNAdjointCheckpointScheme_InsertVector(
                        &mut scheme,
                        step,
                        stage,
                        0.1 * (2 * step + stage) as f64,
                        &y
                    ),
                    SUN_SUCCESS
                );
            }
        }

        /* keep=true: repeated loads all succeed */
        let mut yout = N_VNew_Serial(3, &ctx);
        let mut tout = -1.0;
        for _ in 0..2 {
            assert_eq!(
                SUNAdjointCheckpointScheme_LoadVector(&mut scheme, 1, 1, false, &mut yout, &mut tout),
                SUN_SUCCESS
            );
            assert_eq!(tout, 0.1 * 3.0);
            assert_eq!(yout.data, vec![11.0; 3]);
        }
        assert_eq!(
            SUNAdjointCheckpointScheme_LoadVector(&mut scheme, 0, 0, false, &mut yout, &mut tout),
            SUN_SUCCESS
        );
        assert_eq!(tout, 0.0);
        assert_eq!(yout.data, vec![0.0; 3]);

        /* missing step / stage */
        assert_eq!(
            SUNAdjointCheckpointScheme_LoadVector(&mut scheme, 7, 0, false, &mut yout, &mut tout),
            SUN_ERR_CHECKPOINT_NOT_FOUND
        );
        assert_eq!(
            SUNAdjointCheckpointScheme_LoadVector(&mut scheme, 1, 5, false, &mut yout, &mut tout),
            SUN_ERR_CHECKPOINT_NOT_FOUND
        );
    }

    #[test]
    fn load_consumes_when_not_keeping() {
        let ctx = SUNContext_Create();
        let mut scheme = make_scheme(1, false);
        let mut y = N_VNew_Serial(2, &ctx);

        N_VConst(5.0, &mut y);
        assert_eq!(
            SUNAdjointCheckpointScheme_InsertVector(&mut scheme, 0, 0, 1.0, &y),
            SUN_SUCCESS
        );
        N_VConst(6.0, &mut y);
        assert_eq!(
            SUNAdjointCheckpointScheme_InsertVector(&mut scheme, 0, 1, 2.0, &y),
            SUN_SUCCESS
        );

        let mut yout = N_VNew_Serial(2, &ctx);
        let mut tout = 0.0;

        /* peek does not consume */
        assert_eq!(
            SUNAdjointCheckpointScheme_LoadVector(&mut scheme, 0, 0, true, &mut yout, &mut tout),
            SUN_SUCCESS
        );
        assert_eq!((tout, yout.data[0]), (1.0, 5.0));
        assert_eq!(
            SUNAdjointCheckpointScheme_LoadVector(&mut scheme, 0, 0, false, &mut yout, &mut tout),
            SUN_SUCCESS
        );
        assert_eq!((tout, yout.data[0]), (1.0, 5.0));

        /* stage 0 is gone now; after consuming the remaining stage the
        whole step node is removed */
        assert_eq!(
            SUNAdjointCheckpointScheme_LoadVector(&mut scheme, 0, 0, false, &mut yout, &mut tout),
            SUN_SUCCESS
        );
        assert_eq!((tout, yout.data[0]), (2.0, 6.0));
        assert_eq!(
            SUNAdjointCheckpointScheme_LoadVector(&mut scheme, 0, 0, false, &mut yout, &mut tout),
            SUN_ERR_CHECKPOINT_NOT_FOUND
        );
    }
}
