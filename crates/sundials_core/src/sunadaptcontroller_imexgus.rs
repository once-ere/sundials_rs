/* -----------------------------------------------------------------
 * Translation of
 * sundials-7.7.0/src/sunadaptcontroller/imexgus/sunadaptcontroller_imexgus.c
 *
 * SUNAdaptController_ImExGus: the ImEx Gustafsson controller.
 * Implementation ops called on a controller of another enum variant
 * return SUN_ERR_ARG_INCOMPATIBLE (mismatched content cast is UB in C).
 * -----------------------------------------------------------------*/

use crate::sundials_adaptcontroller::{
    SUNAdaptController, SUNAdaptController_Type, SUNAdaptController_Write,
    SUN_ADAPTCONTROLLER_H,
};
use crate::sundials_errors::{SUNErrCode, SUN_ERR_ARG_CORRUPT, SUN_ERR_ARG_INCOMPATIBLE, SUN_SUCCESS};
use crate::sundials_math::{SUNRpowerR, SUNStrToReal};
use crate::sundials_utils::fmt_g;

/* ------------------
 * Default parameters
 * ------------------ */

const DEFAULT_K1E: f64 = 0.367;
const DEFAULT_K2E: f64 = 0.268;
const DEFAULT_K1I: f64 = 0.95;
const DEFAULT_K2I: f64 = 0.95;
const DEFAULT_BIAS: f64 = 1.0;

/// C struct SUNAdaptControllerContent_ImExGus_
#[derive(Clone, Debug, Default)]
pub struct SUNAdaptControllerContent_ImExGus {
    pub k1e: f64,
    pub k2e: f64,
    pub k1i: f64,
    pub k2i: f64,
    pub bias: f64,
    pub ep: f64,
    pub hp: f64,
    pub firststep: bool,
}

fn content(c: &SUNAdaptController) -> Result<&SUNAdaptControllerContent_ImExGus, SUNErrCode> {
    match c {
        SUNAdaptController::ImExGus(s) => Ok(s),
        _ => Err(SUN_ERR_ARG_INCOMPATIBLE),
    }
}

fn content_mut(
    c: &mut SUNAdaptController,
) -> Result<&mut SUNAdaptControllerContent_ImExGus, SUNErrCode> {
    match c {
        SUNAdaptController::ImExGus(s) => Ok(s),
        _ => Err(SUN_ERR_ARG_INCOMPATIBLE),
    }
}

/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/// Function to create a new ImExGus controller
pub fn SUNAdaptController_ImExGus() -> SUNAdaptController {
    let mut c = SUNAdaptController::ImExGus(SUNAdaptControllerContent_ImExGus::default());
    /* Fill content with default/reset values */
    let _ = SUNAdaptController_SetDefaults_ImExGus(&mut c);
    let _ = SUNAdaptController_Reset_ImExGus(&mut c);
    c
}

/// Function to control set routines via the command line or file
pub fn SUNAdaptController_SetOptions_ImExGus(
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
        return setFromCommandLine_ImExGus(c, cid, args);
    }
    SUN_SUCCESS
}

/// Function to control ImExGus parameters from the command line
fn setFromCommandLine_ImExGus(
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
        if key == "params_imexgus" {
            let r1 = SUNStrToReal(&args[idx + 1]);
            let r2 = SUNStrToReal(&args[idx + 2]);
            let r3 = SUNStrToReal(&args[idx + 3]);
            let r4 = SUNStrToReal(&args[idx + 4]);
            idx += 4;
            let retval = SUNAdaptController_SetParams_ImExGus(c, r1, r2, r3, r4);
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

/// Function to set ImExGus parameters
pub fn SUNAdaptController_SetParams_ImExGus(
    c: &mut SUNAdaptController,
    k1e: f64,
    k2e: f64,
    k1i: f64,
    k2i: f64,
) -> SUNErrCode {
    let s = match content_mut(c) {
        Ok(s) => s,
        Err(e) => return e,
    };
    s.k1e = k1e;
    s.k2e = k2e;
    s.k1i = k1i;
    s.k2i = k2i;
    SUN_SUCCESS
}

/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

pub fn SUNAdaptController_GetType_ImExGus() -> SUNAdaptController_Type {
    SUN_ADAPTCONTROLLER_H
}

pub fn SUNAdaptController_EstimateStep_ImExGus(
    c: &mut SUNAdaptController,
    h: f64,
    p: i32,
    dsm: f64,
    hnew: &mut f64,
) -> SUNErrCode {
    let s = match content(c) {
        Ok(s) => s,
        Err(e) => return e,
    };

    /* order parameter to use */
    let ord = p + 1;

    /* compute estimated time step size, modifying the first step formula */
    if s.firststep {
        /* set usable time-step adaptivity parameters -- first step */
        let k = -1.0 / ord as f64;
        let e = s.bias * dsm;
        *hnew = h * SUNRpowerR(e, k);
    } else {
        /* set usable time-step adaptivity parameters -- subsequent steps */
        let k1e = -s.k1e / ord as f64;
        let k2e = -s.k2e / ord as f64;
        let k1i = -s.k1i / ord as f64;
        let k2i = -s.k2i / ord as f64;
        let e1 = s.bias * dsm;
        let e2 = e1 / s.ep;
        let hrat = h / s.hp;
        *hnew = h * f64::min(
            hrat * SUNRpowerR(e1, k1i) * SUNRpowerR(e2, k2i),
            SUNRpowerR(e1, k1e) * SUNRpowerR(e2, k2e),
        );
    }

    SUN_SUCCESS
}

pub fn SUNAdaptController_Reset_ImExGus(c: &mut SUNAdaptController) -> SUNErrCode {
    let s = match content_mut(c) {
        Ok(s) => s,
        Err(e) => return e,
    };
    s.ep = 1.0;
    s.firststep = true;
    SUN_SUCCESS
}

pub fn SUNAdaptController_SetDefaults_ImExGus(c: &mut SUNAdaptController) -> SUNErrCode {
    {
        let s = match content_mut(c) {
            Ok(s) => s,
            Err(e) => return e,
        };
        s.bias = DEFAULT_BIAS;
    }
    SUNAdaptController_SetParams_ImExGus(c, DEFAULT_K1E, DEFAULT_K2E, DEFAULT_K1I, DEFAULT_K2I)
}

pub fn SUNAdaptController_Write_ImExGus(
    c: &SUNAdaptController,
    fptr: &mut dyn std::io::Write,
) -> SUNErrCode {
    let s = match content(c) {
        Ok(s) => s,
        Err(e) => return e,
    };
    let r = (|| -> std::io::Result<()> {
        writeln!(fptr, "ImEx Gustafsson SUNAdaptController module:")?;
        writeln!(fptr, "  k1e = {}", fmt_g(s.k1e, 0, 15))?;
        writeln!(fptr, "  k2e = {}", fmt_g(s.k2e, 0, 15))?;
        writeln!(fptr, "  k1i = {}", fmt_g(s.k1i, 0, 15))?;
        writeln!(fptr, "  k2i = {}", fmt_g(s.k2i, 0, 15))?;
        writeln!(fptr, "  bias factor = {}", fmt_g(s.bias, 0, 15))?;
        writeln!(fptr, "  previous error = {}", fmt_g(s.ep, 0, 15))?;
        writeln!(fptr, "  previous step = {}", fmt_g(s.hp, 0, 15))?;
        Ok(())
    })();
    if r.is_err() {
        return SUN_ERR_ARG_CORRUPT;
    }
    SUN_SUCCESS
}

pub fn SUNAdaptController_SetErrorBias_ImExGus(
    c: &mut SUNAdaptController,
    bias: f64,
) -> SUNErrCode {
    let s = match content_mut(c) {
        Ok(s) => s,
        Err(e) => return e,
    };
    /* set allowed value, otherwise set default */
    if bias <= 0.0 {
        s.bias = DEFAULT_BIAS;
    } else {
        s.bias = bias;
    }
    SUN_SUCCESS
}

pub fn SUNAdaptController_UpdateH_ImExGus(
    c: &mut SUNAdaptController,
    h: f64,
    dsm: f64,
) -> SUNErrCode {
    let s = match content_mut(c) {
        Ok(s) => s,
        Err(e) => return e,
    };
    s.ep = s.bias * dsm;
    s.hp = h;
    s.firststep = false;
    SUN_SUCCESS
}

pub fn SUNAdaptController_Space_ImExGus(lenrw: &mut i64, leniw: &mut i64) -> SUNErrCode {
    *lenrw = 7;
    *leniw = 1;
    SUN_SUCCESS
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sundials_adaptcontroller::*;

    #[test]
    fn defaults_and_first_step() {
        let mut c = SUNAdaptController_ImExGus();
        assert_eq!(SUNAdaptController_GetType(&c), SUN_ADAPTCONTROLLER_H);
        let mut hnew = 0.0;
        /* first step: I-controller formula */
        assert_eq!(SUNAdaptController_EstimateStep(&mut c, 0.1, 2, 0.5, &mut hnew), SUN_SUCCESS);
        assert_eq!(hnew, 0.1 * SUNRpowerR(0.5, -1.0 / 3.0));
        /* after UpdateH the two-term Gustafsson formula engages */
        assert_eq!(SUNAdaptController_UpdateH(&mut c, 0.1, 0.5), SUN_SUCCESS);
        let mut hnew2 = 0.0;
        assert_eq!(SUNAdaptController_EstimateStep(&mut c, 0.1, 2, 0.5, &mut hnew2), SUN_SUCCESS);
        assert_ne!(hnew2, hnew);
        /* space counts from the C source */
        let (mut lrw, mut liw) = (0i64, 0i64);
        assert_eq!(SUNAdaptController_Space(&c, &mut lrw, &mut liw), SUN_SUCCESS);
        assert_eq!((lrw, liw), (7, 1));
    }

    #[test]
    fn write_output() {
        let c = SUNAdaptController_ImExGus();
        let mut buf: Vec<u8> = Vec::new();
        assert_eq!(SUNAdaptController_Write(&c, &mut buf), SUN_SUCCESS);
        let s = String::from_utf8(buf).unwrap();
        assert!(s.starts_with("ImEx Gustafsson SUNAdaptController module:\n"));
        assert!(s.contains("  k1e = 0.367\n"));
    }
}
