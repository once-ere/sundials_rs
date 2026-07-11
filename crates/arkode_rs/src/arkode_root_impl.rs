/* -----------------------------------------------------------------
 * Translation of sundials-7.7.0/src/arkode/arkode_root_impl.h.
 * Implementation header for ARKODE's root-finding; the routines
 * (arkRootCheck1..3, arkRootfind, ...) live in arkode_root.rs (from
 * arkode_root.c).
 * -----------------------------------------------------------------*/

use crate::arkode_impl::ARKRootFn;
use crate::sundials_types::UserData;

pub const ARK_ROOT_LRW: i64 = 5;
pub const ARK_ROOT_LIW: i64 = 12;

/* Numeric constants */
pub const HUND: f64 = 100.0;

/// struct ARKodeRootMemRec (arkode_root_impl.h)
pub struct ARKodeRootMem {
    pub gfun: Option<ARKRootFn>, /* function g for roots sought                  */
    pub nrtfn: i32,              /* number of components of g                    */
    pub iroots: Vec<i32>,        /* array for root information                   */
    pub rootdir: Vec<i32>,       /* array specifying direction of zero-crossing  */
    pub tlo: f64,                /* nearest endpoint of interval in root search  */
    pub thi: f64,                /* farthest endpoint of interval in root search */
    pub trout: f64,              /* t value returned by rootfinding routine      */
    pub glo: Vec<f64>,           /* saved array of g values at t = tlo           */
    pub ghi: Vec<f64>,           /* saved array of g values at t = thi           */
    pub grout: Vec<f64>,         /* array of g values at t = trout               */
    pub ttol: f64,               /* tolerance on root location                   */
    pub irfnd: i32,              /* flag showing whether last step had a root    */
    pub nge: i64,                /* counter for g evaluations                    */
    pub gactive: Vec<bool>,      /* array with active/inactive event functions   */
    pub mxgnull: i32,            /* num. warning messages about possible g==0    */
    pub root_data: UserData,     /* pointer to user_data                         */
}
