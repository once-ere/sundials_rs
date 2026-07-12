/* -----------------------------------------------------------------
 * Translated from src/arkode/arkode_user_controller.c (+ .h)
 * (ARKODE 7.7.0).  ARKUserControl SUNAdaptController module: wraps a
 * user-supplied (deprecated) ARKAdaptFn as a temporal controller.
 *
 * Storage adaptation (pinned, mirrors arkode_mristep_controller.rs):
 * C builds a SUNAdaptController whose content carries an ark_mem
 * back-pointer and stores it in hadapt_mem->hcontroller.  Safe Rust
 * cannot store the back-pointer, so the wrapper becomes the
 * ARKUserControlContent box held in the Rust-only
 * hadapt_mem.usercontrol slot ("hcontroller is an owned ARKUserControl
 * wrapper"); arkAdapt / arkCompleteStep dispatch its estimatestep and
 * updateh ops with ark_mem in scope, and the reset/write/space ops
 * are forwarded from the hcontroller call sites.
 * -----------------------------------------------------------------*/
use crate::arkode_impl::{ARKAdaptFn, ARKodeMem};
use crate::sundials_adaptcontroller::{
    SUNAdaptController_Type, SUN_ADAPTCONTROLLER_H,
};
use crate::sundials_errors::{SUNErrCode, SUN_ERR_USER_FCN_FAIL, SUN_SUCCESS};
use crate::sundials_types::UserData;
use crate::sundials_utils::fmt_g;

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/* ARKUserControlContent (arkode_user_controller.h); the C ark_mem
   back-pointer is replaced by an &mut ARKodeMem argument at the two
   call sites that need it. */
pub struct ARKUserControlContent {
    pub hp: f64,  /* h from previous step            */
    pub hpp: f64, /* h from 2 steps ago              */
    pub ep: f64,  /* error from previous step        */
    pub epp: f64, /* error from 2 steps ago          */
    pub hadapt: ARKAdaptFn, /* user-provided adaptivity fn */
    pub hadapt_data: UserData, /* user-provided data pointer */
}

/* -----------------------------------------------------------------
 * Function to create a new ARKUserControl controller
 * ----------------------------------------------------------------- */
pub fn ARKUserControl(hadapt: ARKAdaptFn, hadapt_data: UserData) -> Box<ARKUserControlContent> {
    let mut content = Box::new(ARKUserControlContent {
        hp: ZERO,
        hpp: ZERO,
        ep: ONE,
        epp: ONE,
        hadapt,
        hadapt_data,
    });

    /* Fill content with default/reset values */
    let _ = SUNAdaptController_Reset_ARKUserControl(&mut content);

    content
}

/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

pub fn SUNAdaptController_GetType_ARKUserControl(
    _c: &ARKUserControlContent,
) -> SUNAdaptController_Type {
    SUN_ADAPTCONTROLLER_H
}

/* C SUNAdaptController_EstimateStep_ARKUserControl; q/p are the
   hadapt_mem method/embedding orders (read at the call site since
   hadapt_mem holds this content). */
pub fn SUNAdaptController_EstimateStep_ARKUserControl(
    ark_mem: &mut ARKodeMem,
    c: &mut ARKUserControlContent,
    q: i32,
    p: i32,
    h: f64,
    dsm: f64,
    hnew: &mut f64,
) -> SUNErrCode {
    /* call user-provided function to compute new step */
    let ttmp = if dsm <= ONE {
        ark_mem.tn + ark_mem.h
    } else {
        ark_mem.tn
    };
    let hadapt = c.hadapt;
    let retval = {
        let ycur = std::mem::replace(&mut ark_mem.ycur, crate::nvector_serial::NVector::new(0));
        let ret = hadapt(&ycur, ttmp, h, c.hp, c.hpp, dsm, c.ep, c.epp, q, p, hnew, &mut c.hadapt_data);
        ark_mem.ycur = ycur;
        ret
    };
    if retval != SUN_SUCCESS {
        return SUN_ERR_USER_FCN_FAIL;
    }
    SUN_SUCCESS
}

pub fn SUNAdaptController_Reset_ARKUserControl(c: &mut ARKUserControlContent) -> SUNErrCode {
    c.ep = 1.0;
    c.epp = 1.0;
    c.hp = 0.0;
    c.hpp = 0.0;
    SUN_SUCCESS
}

pub fn SUNAdaptController_Write_ARKUserControl(
    c: &ARKUserControlContent,
    fptr: &mut dyn std::io::Write,
) -> SUNErrCode {
    let _ = writeln!(fptr, "ARKUserControl module:");
    let _ = writeln!(fptr, "  hp = {}", fmt_g(c.hp, 0, 15));
    let _ = writeln!(fptr, "  hpp = {}", fmt_g(c.hpp, 0, 15));
    let _ = writeln!(fptr, "  ep = {}", fmt_g(c.ep, 0, 15));
    let _ = writeln!(fptr, "  epp = {}", fmt_g(c.epp, 0, 15));
    let _ = writeln!(fptr, "  hadapt_data = {:p}", &c.hadapt_data);
    SUN_SUCCESS
}

pub fn SUNAdaptController_UpdateH_ARKUserControl(
    c: &mut ARKUserControlContent,
    h: f64,
    dsm: f64,
) -> SUNErrCode {
    c.hpp = c.hp;
    c.hp = h;
    c.epp = c.ep;
    c.ep = dsm;
    SUN_SUCCESS
}

pub fn SUNAdaptController_Space_ARKUserControl(
    _c: &ARKUserControlContent,
    lenrw: &mut i64,
    leniw: &mut i64,
) -> SUNErrCode {
    *lenrw = 4;
    *leniw = 2;
    SUN_SUCCESS
}
