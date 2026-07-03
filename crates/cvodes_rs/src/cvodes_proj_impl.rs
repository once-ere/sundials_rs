/* -----------------------------------------------------------------
 * Translated from src/cvodes/cvodes_proj_impl.h (CVODES 7.7.0).
 * Data structure for projections in CVODE.
 * -----------------------------------------------------------------*/
use crate::nvector_serial::NVector;
use crate::sundials_types::UserData;

/* Default Projection Constants */
pub const PROJ_MAX_FAILS: i32 = 10; /* max projection failures in one step attempt */
pub const PROJ_EPS: f64 = 0.1; /* projection solve tolerance */
pub const PROJ_FAIL_ETA: f64 = 0.25; /* max step size decrease on projection failure */

/* CVProjFn (cvodes_proj.h): user projection function */
pub type CVProjFn = fn(
    t: f64,
    ycur: &NVector,
    corr: &mut NVector,
    eps_proj: f64,
    err: Option<&mut NVector>,
    user_data: &mut UserData,
) -> i32;

/* -----------------------------------------------------------------
   Types : struct CVodeProjMemRec, CVodeProjMem
   -----------------------------------------------------------------*/
pub struct CVodeProjMem {
    pub internal_proj: bool, /* use the internal projection algorithm? */
    pub err_proj: bool,      /* is error projection enabled?           */
    pub first_proj: bool,    /* is this the first time we project?     */

    pub freq: i64,    /* projection frequency           */
    pub nstlprj: i64, /* step number of last projection */

    pub max_fails: i32, /* maximum number of projection failures */

    pub pfun: Option<CVProjFn>, /* function to perform projection */

    pub eps_proj: f64,  /* projection solve tolerance               */
    pub eta_pfail: f64, /* projection failure step reduction factor */

    pub nproj: i64,   /* number of projections performed */
    pub npfails: i64, /* number of projection failures   */
}
