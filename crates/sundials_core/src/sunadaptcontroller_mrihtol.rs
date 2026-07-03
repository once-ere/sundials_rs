/* -----------------------------------------------------------------
 * Translation of
 * sundials-7.7.0/src/sunadaptcontroller/mrihtol/sunadaptcontroller_mrihtol.c
 *
 * SUNAdaptController_MRIHTol: multirate H-Tol controller wrapping a
 * slow step controller (HControl) and a fast tolerance controller
 * (TolControl). In C these are SUNAdaptController pointers stored in
 * the content; here the content owns them as Boxes. Implementation
 * ops called on a controller of another enum variant return
 * SUN_ERR_ARG_INCOMPATIBLE (mismatched content cast is UB in C).
 * -----------------------------------------------------------------*/

use crate::sundials_adaptcontroller::{
    SUNAdaptController, SUNAdaptController_EstimateStep, SUNAdaptController_GetType,
    SUNAdaptController_Reset, SUNAdaptController_SetDefaults, SUNAdaptController_SetErrorBias,
    SUNAdaptController_Space, SUNAdaptController_Type, SUNAdaptController_UpdateH,
    SUNAdaptController_Write, SUN_ADAPTCONTROLLER_H, SUN_ADAPTCONTROLLER_MRI_H_TOL,
};
use crate::sundials_errors::{
    SUNErrCode, SUN_ERR_ARG_CORRUPT, SUN_ERR_ARG_INCOMPATIBLE, SUN_ERR_ARG_OUTOFRANGE, SUN_SUCCESS,
};
use crate::sundials_math::SUNStrToReal;
use crate::sundials_utils::fmt_g;

/* ------------------
 * Default parameters
 * ------------------ */

/// maximum relative change for inner tolerance factor
const INNER_MAX_RELCH: f64 = 20.0;
/// minimum tolerance factor for inner solver
const INNER_MIN_TOLFAC: f64 = 1.0e-5;
/// maximum tolerance factor for inner solver
const INNER_MAX_TOLFAC: f64 = 1.0;

/// C struct SUNAdaptControllerContent_MRIHTol_
pub struct SUNAdaptControllerContent_MRIHTol {
    pub HControl: Box<SUNAdaptController>,
    pub TolControl: Box<SUNAdaptController>,
    pub inner_max_relch: f64,
    pub inner_min_tolfac: f64,
    pub inner_max_tolfac: f64,
}

fn content(c: &SUNAdaptController) -> Result<&SUNAdaptControllerContent_MRIHTol, SUNErrCode> {
    match c {
        SUNAdaptController::MRIHTol(s) => Ok(s),
        _ => Err(SUN_ERR_ARG_INCOMPATIBLE),
    }
}

fn content_mut(
    c: &mut SUNAdaptController,
) -> Result<&mut SUNAdaptControllerContent_MRIHTol, SUNErrCode> {
    match c {
        SUNAdaptController::MRIHTol(s) => Ok(s),
        _ => Err(SUN_ERR_ARG_INCOMPATIBLE),
    }
}

/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/// Function to create a new MRIHTol controller. Takes ownership of the
/// slow (HControl) and fast (TolControl) sub-controllers; both must be
/// of type SUN_ADAPTCONTROLLER_H, else None (C returns NULL).
pub fn SUNAdaptController_MRIHTol(
    HControl: SUNAdaptController,
    TolControl: SUNAdaptController,
) -> Option<SUNAdaptController> {
    /* Verify that input controllers have the appropriate type */
    if SUNAdaptController_GetType(&HControl) != SUN_ADAPTCONTROLLER_H {
        return None;
    }
    if SUNAdaptController_GetType(&TolControl) != SUN_ADAPTCONTROLLER_H {
        return None;
    }

    Some(SUNAdaptController::MRIHTol(SUNAdaptControllerContent_MRIHTol {
        HControl: Box::new(HControl),
        TolControl: Box::new(TolControl),
        /* Set parameters to default values */
        inner_max_relch: INNER_MAX_RELCH,
        inner_min_tolfac: INNER_MIN_TOLFAC,
        inner_max_tolfac: INNER_MAX_TOLFAC,
    }))
}

/// Function to control set routines via the command line or file
pub fn SUNAdaptController_SetOptions_MRIHTol(
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
    if !args.is_empty() {
        return setFromCommandLine_MRIHTol(c, cid, args);
    }
    SUN_SUCCESS
}

/// Function to control MRIHTol parameters from the command line
fn setFromCommandLine_MRIHTol(
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

    let mut write_parameters = false;
    let mut idx = 1;
    while idx < args.len() {
        let arg = &args[idx];
        if !arg.starts_with(&prefix) {
            idx += 1;
            continue;
        }
        let key = &arg[offset..];

        /* control over SetParams function */
        if key == "params_mrihtol" {
            let r1 = SUNStrToReal(&args[idx + 1]);
            let r2 = SUNStrToReal(&args[idx + 2]);
            let r3 = SUNStrToReal(&args[idx + 3]);
            idx += 3;
            let retval = SUNAdaptController_SetParams_MRIHTol(c, r1, r2, r3);
            if retval != SUN_SUCCESS {
                return retval;
            }
            idx += 1;
            continue;
        }

        /* check whether it was requested that all parameters be printed */
        if key == "write_parameters" {
            write_parameters = true;
        }
        idx += 1;
    }

    /* Call SUNAdaptController_Write (if requested) now that all
    command-line options have been set. */
    if write_parameters {
        let retval = SUNAdaptController_Write(c, &mut std::io::stdout());
        if retval != SUN_SUCCESS {
            return retval;
        }
    }
    SUN_SUCCESS
}

/// Function to set MRIHTol parameters
pub fn SUNAdaptController_SetParams_MRIHTol(
    c: &mut SUNAdaptController,
    inner_max_relch: f64,
    inner_min_tolfac: f64,
    inner_max_tolfac: f64,
) -> SUNErrCode {
    let s = match content_mut(c) {
        Ok(s) => s,
        Err(e) => return e,
    };
    if inner_max_tolfac <= inner_min_tolfac {
        return SUN_ERR_ARG_OUTOFRANGE;
    }
    if inner_max_relch < 1.0 {
        s.inner_max_relch = INNER_MAX_RELCH;
    } else {
        s.inner_max_relch = inner_max_relch;
    }
    if inner_min_tolfac <= 0.0 {
        s.inner_min_tolfac = INNER_MIN_TOLFAC;
    } else {
        s.inner_min_tolfac = inner_min_tolfac;
    }
    if inner_max_tolfac <= 0.0 {
        s.inner_max_tolfac = INNER_MAX_TOLFAC;
    } else {
        s.inner_max_tolfac = inner_max_tolfac;
    }
    SUN_SUCCESS
}

/// Function to get the slow sub-controller (C returns a borrowed pointer).
pub fn SUNAdaptController_GetSlowController_MRIHTol(
    c: &mut SUNAdaptController,
) -> Result<&mut SUNAdaptController, SUNErrCode> {
    let s = content_mut(c)?;
    Ok(&mut s.HControl)
}

/// Function to get the fast sub-controller (C returns a borrowed pointer).
pub fn SUNAdaptController_GetFastController_MRIHTol(
    c: &mut SUNAdaptController,
) -> Result<&mut SUNAdaptController, SUNErrCode> {
    let s = content_mut(c)?;
    Ok(&mut s.TolControl)
}

/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

pub fn SUNAdaptController_GetType_MRIHTol() -> SUNAdaptController_Type {
    SUN_ADAPTCONTROLLER_MRI_H_TOL
}

#[allow(clippy::too_many_arguments)]
pub fn SUNAdaptController_EstimateStepTol_MRIHTol(
    c: &mut SUNAdaptController,
    big_h: f64,
    tolfac: f64,
    big_p: i32,
    big_dsm: f64,
    dsm: f64,
    hnew: &mut f64,
    tolfacnew: &mut f64,
) -> SUNErrCode {
    let s = match content_mut(c) {
        Ok(s) => s,
        Err(e) => return e,
    };
    let mut tolfacest = 0.0;

    /* Call slow time scale sub-controller to fill Hnew -- note that all
    heuristics bounds on Hnew will be enforced by the time integrator itself */
    let retval = SUNAdaptController_EstimateStep(&mut s.HControl, big_h, big_p, big_dsm, hnew);
    if retval != SUN_SUCCESS {
        return retval;
    }

    /* Call fast time scale sub-controller with order=1: no matter the
    integrator order, we expect its error to be proportional to the
    tolerance factor */
    let retval = SUNAdaptController_EstimateStep(&mut s.TolControl, tolfac, 0, dsm, &mut tolfacest);
    if retval != SUN_SUCCESS {
        return retval;
    }

    /* Enforce bounds on estimated tolerance factor */
    /*     keep relative change within bounds */
    tolfacest = f64::max(tolfacest, tolfac / s.inner_max_relch);
    tolfacest = f64::min(tolfacest, tolfac * s.inner_max_relch);
    /*     enforce absolute min/max bounds */
    tolfacest = f64::max(tolfacest, s.inner_min_tolfac);
    tolfacest = f64::min(tolfacest, s.inner_max_tolfac);

    /* Set result and return */
    *tolfacnew = tolfacest;
    SUN_SUCCESS
}

pub fn SUNAdaptController_Reset_MRIHTol(c: &mut SUNAdaptController) -> SUNErrCode {
    let s = match content_mut(c) {
        Ok(s) => s,
        Err(e) => return e,
    };
    let retval = SUNAdaptController_Reset(&mut s.HControl);
    if retval != SUN_SUCCESS {
        return retval;
    }
    SUNAdaptController_Reset(&mut s.TolControl)
}

pub fn SUNAdaptController_SetDefaults_MRIHTol(c: &mut SUNAdaptController) -> SUNErrCode {
    let s = match content_mut(c) {
        Ok(s) => s,
        Err(e) => return e,
    };
    let retval = SUNAdaptController_SetDefaults(&mut s.HControl);
    if retval != SUN_SUCCESS {
        return retval;
    }
    let retval = SUNAdaptController_SetDefaults(&mut s.TolControl);
    if retval != SUN_SUCCESS {
        return retval;
    }
    s.inner_max_relch = INNER_MAX_RELCH;
    s.inner_min_tolfac = INNER_MIN_TOLFAC;
    s.inner_max_tolfac = INNER_MAX_TOLFAC;
    SUN_SUCCESS
}

pub fn SUNAdaptController_Write_MRIHTol(
    c: &SUNAdaptController,
    fptr: &mut dyn std::io::Write,
) -> SUNErrCode {
    let s = match content(c) {
        Ok(s) => s,
        Err(e) => return e,
    };
    let r = (|| -> std::io::Result<()> {
        writeln!(fptr, "Multirate H-Tol SUNAdaptController module:")?;
        writeln!(fptr, "  inner_max_relch  = {}", fmt_g(s.inner_max_relch, 0, 15))?;
        writeln!(fptr, "  inner_min_tolfac = {}", fmt_g(s.inner_min_tolfac, 0, 15))?;
        writeln!(fptr, "  inner_max_tolfac = {}", fmt_g(s.inner_max_tolfac, 0, 15))?;
        Ok(())
    })();
    if r.is_err() {
        return SUN_ERR_ARG_CORRUPT;
    }
    let r = writeln!(fptr, "\nSlow step controller:");
    if r.is_err() {
        return SUN_ERR_ARG_CORRUPT;
    }
    let retval = SUNAdaptController_Write(&s.HControl, fptr);
    if retval != SUN_SUCCESS {
        return retval;
    }
    let r = writeln!(fptr, "\nFast tolerance controller:");
    if r.is_err() {
        return SUN_ERR_ARG_CORRUPT;
    }
    SUNAdaptController_Write(&s.TolControl, fptr)
}

pub fn SUNAdaptController_SetErrorBias_MRIHTol(
    c: &mut SUNAdaptController,
    bias: f64,
) -> SUNErrCode {
    let s = match content_mut(c) {
        Ok(s) => s,
        Err(e) => return e,
    };
    let retval = SUNAdaptController_SetErrorBias(&mut s.HControl, bias);
    if retval != SUN_SUCCESS {
        return retval;
    }
    SUNAdaptController_SetErrorBias(&mut s.TolControl, bias)
}

pub fn SUNAdaptController_UpdateMRIHTol_MRIHTol(
    c: &mut SUNAdaptController,
    big_h: f64,
    tolfac: f64,
    big_dsm: f64,
    dsm: f64,
) -> SUNErrCode {
    let s = match content_mut(c) {
        Ok(s) => s,
        Err(e) => return e,
    };
    let retval = SUNAdaptController_UpdateH(&mut s.HControl, big_h, big_dsm);
    if retval != SUN_SUCCESS {
        return retval;
    }
    SUNAdaptController_UpdateH(&mut s.TolControl, tolfac, dsm)
}

pub fn SUNAdaptController_Space_MRIHTol(
    c: &SUNAdaptController,
    lenrw: &mut i64,
    leniw: &mut i64,
) -> SUNErrCode {
    let s = match content(c) {
        Ok(s) => s,
        Err(e) => return e,
    };
    let (mut lrw, mut liw) = (0i64, 0i64);
    let retval = SUNAdaptController_Space(&s.HControl, lenrw, leniw);
    if retval != SUN_SUCCESS {
        return retval;
    }
    let retval = SUNAdaptController_Space(&s.TolControl, &mut lrw, &mut liw);
    if retval != SUN_SUCCESS {
        return retval;
    }
    *lenrw += lrw;
    *leniw += liw;
    SUN_SUCCESS
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sundials_adaptcontroller::*;
    use crate::sunadaptcontroller_soderlind::{SUNAdaptController_I, SUNAdaptController_PI};

    fn make() -> SUNAdaptController {
        SUNAdaptController_MRIHTol(SUNAdaptController_PI(), SUNAdaptController_I()).unwrap()
    }

    #[test]
    fn create_and_type() {
        let c = make();
        assert_eq!(SUNAdaptController_GetType(&c), SUN_ADAPTCONTROLLER_MRI_H_TOL);
        /* wrapping a non-H controller is rejected */
        assert!(SUNAdaptController_MRIHTol(make(), SUNAdaptController_I()).is_none());
    }

    #[test]
    fn estimate_step_tol_bounds() {
        let mut c = make();
        let (mut hnew, mut tolfacnew) = (0.0, 0.0);
        assert_eq!(
            SUNAdaptController_EstimateStepTol(&mut c, 1.0, 1.0, 3, 0.5, 0.5, &mut hnew, &mut tolfacnew),
            SUN_SUCCESS
        );
        assert!(hnew > 0.0);
        /* fast controller called with p=0 (ord=1): raw I estimate = 1.0 * 0.5^-1 = 2.0,
        clipped to inner_max_tolfac = 1.0 */
        assert_eq!(tolfacnew, 1.0);
        /* params: min>=max rejected */
        assert_eq!(
            SUNAdaptController_SetParams_MRIHTol(&mut c, 10.0, 2.0, 1.0),
            SUN_ERR_ARG_OUTOFRANGE
        );
        assert_eq!(SUNAdaptController_SetParams_MRIHTol(&mut c, 0.5, -1.0, 2.0), SUN_SUCCESS);
        match &c {
            SUNAdaptController::MRIHTol(s) => {
                assert_eq!(s.inner_max_relch, 20.0); /* <1 → default */
                assert_eq!(s.inner_min_tolfac, 1.0e-5); /* <=0 → default */
                assert_eq!(s.inner_max_tolfac, 2.0);
            }
            _ => panic!(),
        }
    }

    #[test]
    fn update_and_space() {
        let mut c = make();
        assert_eq!(SUNAdaptController_UpdateMRIHTol(&mut c, 0.1, 0.5, 0.4, 0.3), SUN_SUCCESS);
        let (mut lrw, mut liw) = (0i64, 0i64);
        assert_eq!(SUNAdaptController_Space(&c, &mut lrw, &mut liw), SUN_SUCCESS);
        assert_eq!((lrw, liw), (20, 4)); /* two Soderlind sub-controllers */
        assert_eq!(SUNAdaptController_Reset(&mut c), SUN_SUCCESS);
    }

    #[test]
    fn write_output() {
        let c = make();
        let mut buf: Vec<u8> = Vec::new();
        assert_eq!(SUNAdaptController_Write(&c, &mut buf), SUN_SUCCESS);
        let s = String::from_utf8(buf).unwrap();
        assert!(s.starts_with("Multirate H-Tol SUNAdaptController module:\n"));
        assert!(s.contains("\nSlow step controller:\nSoderlind SUNAdaptController module:\n"));
        assert!(s.contains("\nFast tolerance controller:\nSoderlind SUNAdaptController module:\n"));
    }
}
