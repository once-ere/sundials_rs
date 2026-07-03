/* -----------------------------------------------------------------
 * Translation of sundials-7.7.0/src/sundials/sundials_adaptcontroller.c
 * (+ include/sundials/sundials_adaptcontroller.h).
 *
 * The C generic SUNAdaptController is a base "class" holding an ops
 * table; per the workspace architecture, the ops table becomes enum
 * dispatch. The C-object-model plumbing NewEmpty/DestroyEmpty has no
 * Rust counterpart (enum variants are constructed directly by the
 * implementation constructors in sunadaptcontroller_*.rs); Destroy is
 * ownership drop. All other base-class entry points are ported 1:1.
 * -----------------------------------------------------------------*/

use crate::sundials_errors::{SUNErrCode, SUN_SUCCESS, SUN_ERR_ARG_INCOMPATIBLE, SUN_ERR_ARG_OUTOFRANGE};
use crate::sundials_math::SUNStrToReal;
use crate::sunadaptcontroller_imexgus::SUNAdaptControllerContent_ImExGus;
use crate::sunadaptcontroller_mrihtol::SUNAdaptControllerContent_MRIHTol;
use crate::sunadaptcontroller_soderlind::SUNAdaptControllerContent_Soderlind;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SUNAdaptController_Type {
    SUN_ADAPTCONTROLLER_NONE,
    SUN_ADAPTCONTROLLER_H,
    SUN_ADAPTCONTROLLER_MRI_H_TOL,
}
pub use SUNAdaptController_Type::*;

/// The generic controller: C's ops-table polymorphism as enum dispatch.
pub enum SUNAdaptController {
    Soderlind(SUNAdaptControllerContent_Soderlind),
    ImExGus(SUNAdaptControllerContent_ImExGus),
    MRIHTol(SUNAdaptControllerContent_MRIHTol),
}

/* -----------------------------------------------------------------
 * Required functions
 * ----------------------------------------------------------------- */

pub fn SUNAdaptController_GetType(c: &SUNAdaptController) -> SUNAdaptController_Type {
    match c {
        SUNAdaptController::Soderlind(_) => {
            crate::sunadaptcontroller_soderlind::SUNAdaptController_GetType_Soderlind()
        }
        SUNAdaptController::ImExGus(_) => {
            crate::sunadaptcontroller_imexgus::SUNAdaptController_GetType_ImExGus()
        }
        SUNAdaptController::MRIHTol(_) => {
            crate::sunadaptcontroller_mrihtol::SUNAdaptController_GetType_MRIHTol()
        }
    }
}

/* -----------------------------------------------------------------
 * internal utility routines
 * ----------------------------------------------------------------- */

/// C sunadctrlSetFromCommandLine: process base-class options. `args`
/// corresponds to argv (args[0] is the program name and is skipped).
fn sunadctrlSetFromCommandLine(
    c: &mut SUNAdaptController,
    cid: Option<&str>,
    args: &[String],
) -> SUNErrCode {
    let default_id = "sunadaptcontroller";
    let prefix = match cid {
        Some(id) if !id.is_empty() => format!("{}.", id),
        _ => format!("{}.", default_id),
    };
    let offset = prefix.len();

    let mut idx = 1;
    while idx < args.len() {
        let arg = &args[idx];
        if !arg.starts_with(&prefix) {
            idx += 1;
            continue;
        }
        let key = &arg[offset..];

        /* control over SetDefaults function */
        if key == "defaults" {
            let retval = SUNAdaptController_SetDefaults(c);
            if retval != SUN_SUCCESS {
                return retval;
            }
            idx += 1;
            continue;
        }

        /* control over SetErrorBias function */
        if key == "error_bias" {
            idx += 1;
            let rarg = SUNStrToReal(&args[idx]);
            let retval = SUNAdaptController_SetErrorBias(c, rarg);
            if retval != SUN_SUCCESS {
                return retval;
            }
            idx += 1;
            continue;
        }

        /* Note: SUNAdaptController_Write is processed in the implementations,
        not here (it must run after all options are set). */
        idx += 1;
    }
    SUN_SUCCESS
}

/* -----------------------------------------------------------------
 * Optional functions
 * ----------------------------------------------------------------- */

/// C SUNAdaptController_Destroy: ownership drop in Rust.
pub fn SUNAdaptController_Destroy(c: SUNAdaptController) -> SUNErrCode {
    drop(c);
    SUN_SUCCESS
}

pub fn SUNAdaptController_EstimateStep(
    c: &mut SUNAdaptController,
    h: f64,
    p: i32,
    dsm: f64,
    hnew: &mut f64,
) -> SUNErrCode {
    if !h.is_finite() || p < 0 || dsm < 0.0 {
        return SUN_ERR_ARG_OUTOFRANGE;
    }
    *hnew = h; /* initialize output with identity */
    match c {
        SUNAdaptController::Soderlind(_) => {
            crate::sunadaptcontroller_soderlind::SUNAdaptController_EstimateStep_Soderlind(
                c, h, p, dsm, hnew,
            )
        }
        SUNAdaptController::ImExGus(_) => {
            crate::sunadaptcontroller_imexgus::SUNAdaptController_EstimateStep_ImExGus(
                c, h, p, dsm, hnew,
            )
        }
        SUNAdaptController::MRIHTol(_) => SUN_SUCCESS, /* op not provided */
    }
}

#[allow(clippy::too_many_arguments)]
pub fn SUNAdaptController_EstimateStepTol(
    c: &mut SUNAdaptController,
    big_h: f64,
    tolfac: f64,
    big_p: i32,
    big_dsm: f64,
    dsm: f64,
    hnew: &mut f64,
    tolfacnew: &mut f64,
) -> SUNErrCode {
    *hnew = big_h; /* initialize outputs with identity */
    *tolfacnew = tolfac;
    match c {
        SUNAdaptController::MRIHTol(_) => {
            crate::sunadaptcontroller_mrihtol::SUNAdaptController_EstimateStepTol_MRIHTol(
                c, big_h, tolfac, big_p, big_dsm, dsm, hnew, tolfacnew,
            )
        }
        _ => SUN_SUCCESS, /* op not provided */
    }
}

pub fn SUNAdaptController_Reset(c: &mut SUNAdaptController) -> SUNErrCode {
    match c {
        SUNAdaptController::Soderlind(_) => {
            crate::sunadaptcontroller_soderlind::SUNAdaptController_Reset_Soderlind(c)
        }
        SUNAdaptController::ImExGus(_) => {
            crate::sunadaptcontroller_imexgus::SUNAdaptController_Reset_ImExGus(c)
        }
        SUNAdaptController::MRIHTol(_) => {
            crate::sunadaptcontroller_mrihtol::SUNAdaptController_Reset_MRIHTol(c)
        }
    }
}

/// C SUNAdaptController_SetOptions(C, Cid, file_name, argc, argv).
pub fn SUNAdaptController_SetOptions(
    c: &mut SUNAdaptController,
    cid: Option<&str>,
    file_name: Option<&str>,
    args: &[String],
) -> SUNErrCode {
    /* File-based option control is currently unimplemented */
    if let Some(f) = file_name {
        if !f.is_empty() {
            return SUN_ERR_ARG_INCOMPATIBLE;
        }
    }

    /* First, process all base-class options */
    if !args.is_empty() {
        let retval = sunadctrlSetFromCommandLine(c, cid, args);
        if retval != SUN_SUCCESS {
            return retval;
        }
    }

    /* Second, ask the implementation to process any remaining options */
    match c {
        SUNAdaptController::Soderlind(_) => {
            crate::sunadaptcontroller_soderlind::SUNAdaptController_SetOptions_Soderlind(
                c, cid, file_name, args,
            )
        }
        SUNAdaptController::ImExGus(_) => {
            crate::sunadaptcontroller_imexgus::SUNAdaptController_SetOptions_ImExGus(
                c, cid, file_name, args,
            )
        }
        SUNAdaptController::MRIHTol(_) => {
            crate::sunadaptcontroller_mrihtol::SUNAdaptController_SetOptions_MRIHTol(
                c, cid, file_name, args,
            )
        }
    }
}

pub fn SUNAdaptController_SetDefaults(c: &mut SUNAdaptController) -> SUNErrCode {
    match c {
        SUNAdaptController::Soderlind(_) => {
            crate::sunadaptcontroller_soderlind::SUNAdaptController_SetDefaults_Soderlind(c)
        }
        SUNAdaptController::ImExGus(_) => {
            crate::sunadaptcontroller_imexgus::SUNAdaptController_SetDefaults_ImExGus(c)
        }
        SUNAdaptController::MRIHTol(_) => {
            crate::sunadaptcontroller_mrihtol::SUNAdaptController_SetDefaults_MRIHTol(c)
        }
    }
}

pub fn SUNAdaptController_Write(
    c: &SUNAdaptController,
    fptr: &mut dyn std::io::Write,
) -> SUNErrCode {
    match c {
        SUNAdaptController::Soderlind(_) => {
            crate::sunadaptcontroller_soderlind::SUNAdaptController_Write_Soderlind(c, fptr)
        }
        SUNAdaptController::ImExGus(_) => {
            crate::sunadaptcontroller_imexgus::SUNAdaptController_Write_ImExGus(c, fptr)
        }
        SUNAdaptController::MRIHTol(_) => {
            crate::sunadaptcontroller_mrihtol::SUNAdaptController_Write_MRIHTol(c, fptr)
        }
    }
}

pub fn SUNAdaptController_SetErrorBias(c: &mut SUNAdaptController, bias: f64) -> SUNErrCode {
    match c {
        SUNAdaptController::Soderlind(_) => {
            crate::sunadaptcontroller_soderlind::SUNAdaptController_SetErrorBias_Soderlind(c, bias)
        }
        SUNAdaptController::ImExGus(_) => {
            crate::sunadaptcontroller_imexgus::SUNAdaptController_SetErrorBias_ImExGus(c, bias)
        }
        SUNAdaptController::MRIHTol(_) => {
            crate::sunadaptcontroller_mrihtol::SUNAdaptController_SetErrorBias_MRIHTol(c, bias)
        }
    }
}

pub fn SUNAdaptController_UpdateH(c: &mut SUNAdaptController, h: f64, dsm: f64) -> SUNErrCode {
    if !h.is_finite() || dsm < 0.0 {
        return SUN_ERR_ARG_OUTOFRANGE;
    }
    match c {
        SUNAdaptController::Soderlind(_) => {
            crate::sunadaptcontroller_soderlind::SUNAdaptController_UpdateH_Soderlind(c, h, dsm)
        }
        SUNAdaptController::ImExGus(_) => {
            crate::sunadaptcontroller_imexgus::SUNAdaptController_UpdateH_ImExGus(c, h, dsm)
        }
        SUNAdaptController::MRIHTol(_) => SUN_SUCCESS, /* op not provided */
    }
}

pub fn SUNAdaptController_UpdateMRIHTol(
    c: &mut SUNAdaptController,
    big_h: f64,
    tolfac: f64,
    big_dsm: f64,
    dsm: f64,
) -> SUNErrCode {
    match c {
        SUNAdaptController::MRIHTol(_) => {
            crate::sunadaptcontroller_mrihtol::SUNAdaptController_UpdateMRIHTol_MRIHTol(
                c, big_h, tolfac, big_dsm, dsm,
            )
        }
        _ => SUN_SUCCESS, /* op not provided */
    }
}

pub fn SUNAdaptController_Space(
    c: &SUNAdaptController,
    lenrw: &mut i64,
    leniw: &mut i64,
) -> SUNErrCode {
    *lenrw = 0; /* initialize outputs with identity */
    *leniw = 0;
    match c {
        SUNAdaptController::Soderlind(_) => {
            crate::sunadaptcontroller_soderlind::SUNAdaptController_Space_Soderlind(lenrw, leniw)
        }
        SUNAdaptController::ImExGus(_) => {
            crate::sunadaptcontroller_imexgus::SUNAdaptController_Space_ImExGus(lenrw, leniw)
        }
        SUNAdaptController::MRIHTol(_) => {
            crate::sunadaptcontroller_mrihtol::SUNAdaptController_Space_MRIHTol(c, lenrw, leniw)
        }
    }
}
