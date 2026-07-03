/* -----------------------------------------------------------------
 * Translation of
 * sundials-7.7.0/src/sunadaptcontroller/soderlind/sunadaptcontroller_soderlind.c
 *
 * SUNAdaptController_Soderlind (a.k.a. H_{0}321) and its derived
 * parameterizations (PID, PI, I, ExpGus, ImpGus, H0211, H0321, H211,
 * H312). The C content struct is the Soderlind variant payload of the
 * SUNAdaptController enum; implementation ops called on a controller
 * of another variant return SUN_ERR_ARG_INCOMPATIBLE (in C this would
 * be undefined behavior through a mismatched content cast).
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

const DEFAULT_K1: f64 = 1.25; /* H_{0}321 parameters */
const DEFAULT_K2: f64 = 0.5;
const DEFAULT_K3: f64 = -0.75;
const DEFAULT_K4: f64 = 0.25;
const DEFAULT_K5: f64 = 0.75;
const DEFAULT_PID_K1: f64 = 0.58; /* PID parameters */
const DEFAULT_PID_K2: f64 = -0.21;
const DEFAULT_PID_K3: f64 = 0.1;
const DEFAULT_PI_K1: f64 = 0.8; /* PI parameters */
const DEFAULT_PI_K2: f64 = -0.31;
const DEFAULT_I_K1: f64 = 1.0; /* I parameters */
const DEFAULT_EXPGUS_K1: f64 = 0.367; /* Explicit Gustafsson parameters */
const DEFAULT_EXPGUS_K2: f64 = 0.268;
const DEFAULT_IMPGUS_K1: f64 = 0.98; /* Implicit Gustafsson parameters */
const DEFAULT_IMPGUS_K2: f64 = 0.95;
const DEFAULT_BIAS: f64 = 1.0;

/// C struct SUNAdaptControllerContent_Soderlind_
#[derive(Clone, Debug, Default)]
pub struct SUNAdaptControllerContent_Soderlind {
    pub k1: f64,
    pub k2: f64,
    pub k3: f64,
    pub k4: f64,
    pub k5: f64,
    pub bias: f64,
    pub ep: f64,
    pub epp: f64,
    pub hp: f64,
    pub hpp: f64,
    pub firststeps: i32,
    pub historysize: i32,
}

fn content(c: &SUNAdaptController) -> Result<&SUNAdaptControllerContent_Soderlind, SUNErrCode> {
    match c {
        SUNAdaptController::Soderlind(s) => Ok(s),
        _ => Err(SUN_ERR_ARG_INCOMPATIBLE),
    }
}

fn content_mut(
    c: &mut SUNAdaptController,
) -> Result<&mut SUNAdaptControllerContent_Soderlind, SUNErrCode> {
    match c {
        SUNAdaptController::Soderlind(s) => Ok(s),
        _ => Err(SUN_ERR_ARG_INCOMPATIBLE),
    }
}

/* -----------------------------------------------------------------
 * exported functions
 * ----------------------------------------------------------------- */

/// Function to create a new Soderlind controller (a.k.a., H_{0}321)
pub fn SUNAdaptController_Soderlind() -> SUNAdaptController {
    let mut c = SUNAdaptController::Soderlind(SUNAdaptControllerContent_Soderlind::default());
    /* Fill content with default/reset values */
    let _ = SUNAdaptController_SetDefaults_Soderlind(&mut c);
    let _ = SUNAdaptController_Reset_Soderlind(&mut c);
    c
}

/// Function to control set routines via the command line or file
pub fn SUNAdaptController_SetOptions_Soderlind(
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
        return setFromCommandLine_Soderlind(c, cid, args);
    }
    SUN_SUCCESS
}

/// Function to control Soderlind parameters from the command line
fn setFromCommandLine_Soderlind(
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

        match key {
            "params_soderlind" => {
                let r1 = SUNStrToReal(&args[idx + 1]);
                let r2 = SUNStrToReal(&args[idx + 2]);
                let r3 = SUNStrToReal(&args[idx + 3]);
                let r4 = SUNStrToReal(&args[idx + 4]);
                let r5 = SUNStrToReal(&args[idx + 5]);
                idx += 5;
                let retval = SUNAdaptController_SetParams_Soderlind(c, r1, r2, r3, r4, r5);
                if retval != SUN_SUCCESS {
                    return retval;
                }
            }
            "params_pid" => {
                let r1 = SUNStrToReal(&args[idx + 1]);
                let r2 = SUNStrToReal(&args[idx + 2]);
                let r3 = SUNStrToReal(&args[idx + 3]);
                idx += 3;
                let retval = SUNAdaptController_SetParams_PID(c, r1, r2, r3);
                if retval != SUN_SUCCESS {
                    return retval;
                }
            }
            "params_pi" => {
                let r1 = SUNStrToReal(&args[idx + 1]);
                let r2 = SUNStrToReal(&args[idx + 2]);
                idx += 2;
                let retval = SUNAdaptController_SetParams_PI(c, r1, r2);
                if retval != SUN_SUCCESS {
                    return retval;
                }
            }
            "params_i" => {
                let r1 = SUNStrToReal(&args[idx + 1]);
                idx += 1;
                let retval = SUNAdaptController_SetParams_I(c, r1);
                if retval != SUN_SUCCESS {
                    return retval;
                }
            }
            "params_expgus" => {
                let r1 = SUNStrToReal(&args[idx + 1]);
                let r2 = SUNStrToReal(&args[idx + 2]);
                idx += 2;
                let retval = SUNAdaptController_SetParams_ExpGus(c, r1, r2);
                if retval != SUN_SUCCESS {
                    return retval;
                }
            }
            "params_impgus" => {
                let r1 = SUNStrToReal(&args[idx + 1]);
                let r2 = SUNStrToReal(&args[idx + 2]);
                idx += 2;
                let retval = SUNAdaptController_SetParams_ImpGus(c, r1, r2);
                if retval != SUN_SUCCESS {
                    return retval;
                }
            }
            "write_parameters" => {
                write_parameters = true;
            }
            _ => {}
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

/// Function to set Soderlind parameters
pub fn SUNAdaptController_SetParams_Soderlind(
    c: &mut SUNAdaptController,
    k1: f64,
    k2: f64,
    k3: f64,
    k4: f64,
    k5: f64,
) -> SUNErrCode {
    let s = match content_mut(c) {
        Ok(s) => s,
        Err(e) => return e,
    };
    s.k1 = k1;
    s.k2 = k2;
    s.k3 = k3;
    s.k4 = k4;
    s.k5 = k5;

    if k5 != 0.0 || k3 != 0.0 {
        s.historysize = 2;
    } else if k4 != 0.0 || k2 != 0.0 {
        s.historysize = 1;
    } else {
        s.historysize = 0;
    }
    SUN_SUCCESS
}

/// Function to create a PID controller (subset of Soderlind)
pub fn SUNAdaptController_PID() -> SUNAdaptController {
    let mut c = SUNAdaptController_Soderlind();
    let _ = SUNAdaptController_SetParams_PID(&mut c, DEFAULT_PID_K1, DEFAULT_PID_K2, DEFAULT_PID_K3);
    c
}

/// Function to set PID parameters
pub fn SUNAdaptController_SetParams_PID(
    c: &mut SUNAdaptController,
    k1: f64,
    k2: f64,
    k3: f64,
) -> SUNErrCode {
    SUNAdaptController_SetParams_Soderlind(c, k1, k2, k3, 0.0, 0.0)
}

/// Function to create a PI controller (subset of Soderlind)
pub fn SUNAdaptController_PI() -> SUNAdaptController {
    let mut c = SUNAdaptController_Soderlind();
    let _ = SUNAdaptController_SetParams_PI(&mut c, DEFAULT_PI_K1, DEFAULT_PI_K2);
    c
}

/// Function to set PI parameters
pub fn SUNAdaptController_SetParams_PI(c: &mut SUNAdaptController, k1: f64, k2: f64) -> SUNErrCode {
    SUNAdaptController_SetParams_Soderlind(c, k1, k2, 0.0, 0.0, 0.0)
}

/// Function to create an I controller (subset of Soderlind)
pub fn SUNAdaptController_I() -> SUNAdaptController {
    let mut c = SUNAdaptController_Soderlind();
    let _ = SUNAdaptController_SetParams_I(&mut c, DEFAULT_I_K1);
    c
}

/// Function to set I parameters
pub fn SUNAdaptController_SetParams_I(c: &mut SUNAdaptController, k1: f64) -> SUNErrCode {
    SUNAdaptController_SetParams_Soderlind(c, k1, 0.0, 0.0, 0.0, 0.0)
}

/// Function to create an explicit Gustafsson controller
pub fn SUNAdaptController_ExpGus() -> SUNAdaptController {
    let mut c = SUNAdaptController_Soderlind();
    let _ = SUNAdaptController_SetParams_ExpGus(&mut c, DEFAULT_EXPGUS_K1, DEFAULT_EXPGUS_K2);
    c
}

/// Function to set explicit Gustafsson parameters
pub fn SUNAdaptController_SetParams_ExpGus(
    c: &mut SUNAdaptController,
    k1: f64,
    k2: f64,
) -> SUNErrCode {
    SUNAdaptController_SetParams_Soderlind(c, k1 + k2, -k2, 0.0, 0.0, 0.0)
}

/// Function to create an implicit Gustafsson controller
pub fn SUNAdaptController_ImpGus() -> SUNAdaptController {
    let mut c = SUNAdaptController_Soderlind();
    let _ = SUNAdaptController_SetParams_ImpGus(&mut c, DEFAULT_IMPGUS_K1, DEFAULT_IMPGUS_K2);
    c
}

/// Function to set implicit Gustafsson parameters
pub fn SUNAdaptController_SetParams_ImpGus(
    c: &mut SUNAdaptController,
    k1: f64,
    k2: f64,
) -> SUNErrCode {
    SUNAdaptController_SetParams_Soderlind(c, k1 + k2, -k2, 0.0, 1.0, 0.0)
}

/// Function to create an H_{0}211 controller (subset of Soderlind)
pub fn SUNAdaptController_H0211() -> SUNAdaptController {
    let mut c = SUNAdaptController_Soderlind();
    let _ = SUNAdaptController_SetParams_Soderlind(&mut c, 0.5, 0.5, 0.0, -0.5, 0.0);
    c
}

/// Function to create an H_{0}321 controller (subset of Soderlind)
pub fn SUNAdaptController_H0321() -> SUNAdaptController {
    let mut c = SUNAdaptController_Soderlind();
    let _ = SUNAdaptController_SetParams_Soderlind(&mut c, 1.25, 0.5, -0.75, 0.25, 0.75);
    c
}

/// Function to create an H211 controller (subset of Soderlind)
pub fn SUNAdaptController_H211() -> SUNAdaptController {
    let mut c = SUNAdaptController_Soderlind();
    let _ = SUNAdaptController_SetParams_Soderlind(&mut c, 0.25, 0.25, 0.0, -0.25, 0.0);
    c
}

/// Function to create an H312 controller (subset of Soderlind)
pub fn SUNAdaptController_H312() -> SUNAdaptController {
    let mut c = SUNAdaptController_Soderlind();
    let _ = SUNAdaptController_SetParams_Soderlind(
        &mut c,
        1.0 / 8.0,
        0.25,
        1.0 / 8.0,
        -3.0 / 8.0,
        -1.0 / 8.0,
    );
    c
}

/* -----------------------------------------------------------------
 * implementation of controller operations
 * ----------------------------------------------------------------- */

pub fn SUNAdaptController_GetType_Soderlind() -> SUNAdaptController_Type {
    SUN_ADAPTCONTROLLER_H
}

pub fn SUNAdaptController_EstimateStep_Soderlind(
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
    let e1 = s.bias * dsm;

    /* Handle the case of insufficient history */
    if s.firststeps < s.historysize {
        /* Fall back onto an I controller */
        *hnew = h * SUNRpowerR(e1, -1.0 / ord as f64);
        return SUN_SUCCESS;
    }

    let k1 = -s.k1 / ord as f64;
    *hnew = h * SUNRpowerR(e1, k1);

    /* This branching is not ideal, but it's more efficient than computing
     * extra math operations with degenerate k values. */
    if s.historysize > 0 {
        let k2 = -s.k2 / ord as f64;
        let hrat1 = h / s.hp;
        *hnew *= SUNRpowerR(s.ep, k2) * SUNRpowerR(hrat1, s.k4);

        if s.historysize > 1 {
            let k3 = -s.k3 / ord as f64;
            let hrat2 = s.hp / s.hpp;
            *hnew *= SUNRpowerR(s.epp, k3) * SUNRpowerR(hrat2, s.k5);
        }
    }

    SUN_SUCCESS
}

pub fn SUNAdaptController_Reset_Soderlind(c: &mut SUNAdaptController) -> SUNErrCode {
    let s = match content_mut(c) {
        Ok(s) => s,
        Err(e) => return e,
    };
    s.ep = 1.0;
    s.epp = 1.0;
    s.hp = 1.0;
    s.hpp = 1.0;
    s.firststeps = 0;
    SUN_SUCCESS
}

pub fn SUNAdaptController_SetDefaults_Soderlind(c: &mut SUNAdaptController) -> SUNErrCode {
    {
        let s = match content_mut(c) {
            Ok(s) => s,
            Err(e) => return e,
        };
        s.bias = DEFAULT_BIAS;
    }
    SUNAdaptController_SetParams_Soderlind(c, DEFAULT_K1, DEFAULT_K2, DEFAULT_K3, DEFAULT_K4,
                                           DEFAULT_K5)
}

pub fn SUNAdaptController_Write_Soderlind(
    c: &SUNAdaptController,
    fptr: &mut dyn std::io::Write,
) -> SUNErrCode {
    let s = match content(c) {
        Ok(s) => s,
        Err(e) => return e,
    };
    let r = (|| -> std::io::Result<()> {
        writeln!(fptr, "Soderlind SUNAdaptController module:")?;
        writeln!(fptr, "  k1 = {}", fmt_g(s.k1, 0, 15))?;
        writeln!(fptr, "  k2 = {}", fmt_g(s.k2, 0, 15))?;
        writeln!(fptr, "  k3 = {}", fmt_g(s.k3, 0, 15))?;
        writeln!(fptr, "  k4 = {}", fmt_g(s.k4, 0, 15))?;
        writeln!(fptr, "  k5 = {}", fmt_g(s.k5, 0, 15))?;
        writeln!(fptr, "  bias factor = {}", fmt_g(s.bias, 0, 15))?;
        writeln!(fptr, "  previous error = {}", fmt_g(s.ep, 0, 15))?;
        writeln!(fptr, "  previous-previous error = {}", fmt_g(s.epp, 0, 15))?;
        writeln!(fptr, "  previous step = {}", fmt_g(s.hp, 0, 15))?;
        writeln!(fptr, "  previous-previous step = {}", fmt_g(s.hpp, 0, 15))?;
        writeln!(fptr, "  firststeps = {}", s.firststeps)?;
        writeln!(fptr, "  historysize = {}", s.historysize)?;
        Ok(())
    })();
    if r.is_err() {
        return SUN_ERR_ARG_CORRUPT;
    }
    SUN_SUCCESS
}

pub fn SUNAdaptController_SetErrorBias_Soderlind(
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

pub fn SUNAdaptController_UpdateH_Soderlind(
    c: &mut SUNAdaptController,
    h: f64,
    dsm: f64,
) -> SUNErrCode {
    let s = match content_mut(c) {
        Ok(s) => s,
        Err(e) => return e,
    };
    s.epp = s.ep;
    s.ep = s.bias * dsm;
    s.hpp = s.hp;
    s.hp = h;
    if s.firststeps < s.historysize {
        s.firststeps += 1;
    }
    SUN_SUCCESS
}

pub fn SUNAdaptController_Space_Soderlind(lenrw: &mut i64, leniw: &mut i64) -> SUNErrCode {
    *lenrw = 10;
    *leniw = 2;
    SUN_SUCCESS
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sundials_adaptcontroller::*;

    #[test]
    fn defaults_and_reset() {
        let c = SUNAdaptController_Soderlind();
        match &c {
            SUNAdaptController::Soderlind(s) => {
                assert_eq!(s.k1, 1.25);
                assert_eq!(s.k2, 0.5);
                assert_eq!(s.k3, -0.75);
                assert_eq!(s.k4, 0.25);
                assert_eq!(s.k5, 0.75);
                assert_eq!(s.bias, 1.0);
                assert_eq!(s.ep, 1.0);
                assert_eq!(s.historysize, 2);
                assert_eq!(s.firststeps, 0);
            }
            _ => panic!(),
        }
        assert_eq!(SUNAdaptController_GetType(&c), SUN_ADAPTCONTROLLER_H);
    }

    #[test]
    fn i_controller_estimate() {
        /* An I controller has no history: hnew = h * e^(-1/(p+1)) */
        let mut c = SUNAdaptController_I();
        let mut hnew = 0.0;
        let h = 0.1;
        let p = 2;
        let dsm = 0.5;
        assert_eq!(SUNAdaptController_EstimateStep(&mut c, h, p, dsm, &mut hnew), SUN_SUCCESS);
        let expect = h * SUNRpowerR(dsm, -1.0 / 3.0);
        assert_eq!(hnew, expect);
    }

    #[test]
    fn history_progression() {
        let mut c = SUNAdaptController_Soderlind();
        let mut hnew = 0.0;
        /* first two calls fall back to I controller until history fills */
        assert_eq!(SUNAdaptController_EstimateStep(&mut c, 0.1, 2, 0.5, &mut hnew), SUN_SUCCESS);
        let i_est = 0.1 * SUNRpowerR(0.5, -1.0 / 3.0);
        assert_eq!(hnew, i_est);
        assert_eq!(SUNAdaptController_UpdateH(&mut c, 0.09, 0.3), SUN_SUCCESS);
        assert_eq!(SUNAdaptController_EstimateStep(&mut c, 0.1, 2, 0.5, &mut hnew), SUN_SUCCESS);
        assert_eq!(hnew, i_est);
        assert_eq!(SUNAdaptController_UpdateH(&mut c, 0.11, 0.7), SUN_SUCCESS);
        /* history now full (asymmetric, so the H0321 formula cannot coincide
        with the I fallback): the full formula engages */
        assert_eq!(SUNAdaptController_EstimateStep(&mut c, 0.1, 2, 0.5, &mut hnew), SUN_SUCCESS);
        assert_ne!(hnew, i_est);
        /* cross-check against the explicit H0321 expression */
        let ord = 3.0;
        let expect = 0.1
            * SUNRpowerR(0.5, -1.25 / ord)
            * (SUNRpowerR(0.7, -0.5 / ord) * SUNRpowerR(0.1 / 0.11, 0.25))
            * (SUNRpowerR(0.3, 0.75 / ord) * SUNRpowerR(0.11 / 0.09, 0.75));
        assert!((hnew - expect).abs() <= 1e-15 * expect.abs());
        /* reset restores the fallback */
        assert_eq!(SUNAdaptController_Reset(&mut c), SUN_SUCCESS);
        assert_eq!(SUNAdaptController_EstimateStep(&mut c, 0.1, 2, 0.5, &mut hnew), SUN_SUCCESS);
        assert_eq!(hnew, i_est);
    }

    #[test]
    fn expgus_impgus_param_mapping() {
        let c = SUNAdaptController_ExpGus();
        match &c {
            SUNAdaptController::Soderlind(s) => {
                assert_eq!(s.k1, 0.367 + 0.268);
                assert_eq!(s.k2, -0.268);
                assert_eq!(s.historysize, 1);
            }
            _ => panic!(),
        }
        let c = SUNAdaptController_ImpGus();
        match &c {
            SUNAdaptController::Soderlind(s) => {
                assert_eq!(s.k1, 0.98 + 0.95);
                assert_eq!(s.k2, -0.95);
                assert_eq!(s.k4, 1.0);
                assert_eq!(s.historysize, 1);
            }
            _ => panic!(),
        }
    }

    #[test]
    fn error_bias_and_space() {
        let mut c = SUNAdaptController_PID();
        assert_eq!(SUNAdaptController_SetErrorBias(&mut c, 1.5), SUN_SUCCESS);
        match &c {
            SUNAdaptController::Soderlind(s) => assert_eq!(s.bias, 1.5),
            _ => panic!(),
        }
        assert_eq!(SUNAdaptController_SetErrorBias(&mut c, -1.0), SUN_SUCCESS);
        match &c {
            SUNAdaptController::Soderlind(s) => assert_eq!(s.bias, 1.0),
            _ => panic!(),
        }
        let (mut lrw, mut liw) = (0i64, 0i64);
        assert_eq!(SUNAdaptController_Space(&c, &mut lrw, &mut liw), SUN_SUCCESS);
        assert_eq!((lrw, liw), (10, 2));
    }

    #[test]
    fn write_output() {
        let c = SUNAdaptController_Soderlind();
        let mut buf: Vec<u8> = Vec::new();
        assert_eq!(SUNAdaptController_Write(&c, &mut buf), SUN_SUCCESS);
        let s = String::from_utf8(buf).unwrap();
        assert!(s.starts_with("Soderlind SUNAdaptController module:\n"));
        assert!(s.contains("  k1 = 1.25\n"));
        assert!(s.contains("  historysize = 2\n"));
    }
}
