/* -----------------------------------------------------------------
 * Translation of sundials-7.7.0/src/arkode/arkode_adapt_impl.h.
 * Implementation header for ARKODE's time step adaptivity utilities;
 * the routines (arkAdaptInit / arkPrintAdaptMem / arkAdapt) live in
 * arkode_adapt.rs (from arkode_adapt.c).
 * -----------------------------------------------------------------*/

use crate::arkode_impl::ARKExpStabFn;
use crate::sundials_adaptcontroller::SUNAdaptController;
use crate::sundials_types::UserData;

/* size constants for the adaptivity memory structure */
pub const ARK_ADAPT_LRW: i64 = 10;
pub const ARK_ADAPT_LIW: i64 = 7; /* includes function/data pointers */

/* Time step controller default values */
pub const CFLFAC: f64 = 0.5;
pub const SAFETY: f64 = 0.9; /* CVODE uses 1.0  */
pub const GROWTH: f64 = 20.0; /* CVODE uses 10.0 */
pub const HFIXED_LB: f64 = 1.0; /* CVODE uses 1.0  */
pub const HFIXED_UB: f64 = 1.0; /* CVODE uses 1.5  */

/* maximum step size change on first step */
pub const ETAMX1: f64 = 10000.0;
/* step size reduction factor on multiple error test failures (multiple implies >= SMALL_NEF) */
pub const ETAMXF: f64 = 0.3;
/* smallest allowable step size reduction factor on an error test failure */
pub const ETAMIN: f64 = 0.1;
/* step size reduction factor on nonlinear convergence failure */
pub const ETACF: f64 = 0.25;
/* if an error failure occurs and SMALL_NEF <= nef, then reset eta = MIN(eta, ETAMXF) */
pub const SMALL_NEF: i32 = 2;
/* order to use for controller: 0=embedding, 1=method, otherwise
   min(method,embedding). DEPRECATED, REMOVE AT SAME TIME AS
   ARKStepSetAdaptivityMethod */
pub const PQ: i32 = 0;
/* adjustment to apply within controller to method order of accuracy */
pub const ADJUST: i32 = 0;

/// struct ARKodeHAdaptMemRec (arkode_adapt_impl.h)
pub struct ARKodeHAdaptMem {
    pub etamax: f64,    /* eta <= etamax                              */
    pub etamx1: f64,    /* max step size change on first step         */
    pub etamxf: f64,    /* h reduction factor on multiple error fails */
    pub etamin: f64,    /* eta >= etamin on error test fail           */
    pub small_nef: i32, /* bound to determine 'multiple' above        */
    pub etacf: f64,     /* h reduction factor on nonlinear conv fail  */
    pub cfl: f64,       /* cfl safety factor                          */
    pub safety: f64,    /* accuracy safety factor on h                */
    pub growth: f64,    /* maximum step growth safety factor          */
    pub lbound: f64,    /* eta lower bound to leave h unchanged       */
    pub ubound: f64,    /* eta upper bound to leave h unchanged       */
    pub p: i32,         /* embedding order                            */
    pub q: i32,         /* method order                               */
    pub pq: i32,        /* decision flag for controller order         */
    pub adjust: i32,    /* controller order adjustment factor         */

    pub hcontroller: Option<SUNAdaptController>, /* temporal error controller  */
    pub owncontroller: bool,                     /* hcontroller ownership flag */
    /* Rust-only slot: "hcontroller is an owned ARKUserControl wrapper"
       (arkode_user_controller.rs; C stores the wrapper in hcontroller
       with an ark_mem back-pointer that safe Rust cannot express).
       Invariant: usercontrol.is_some() => hcontroller.is_none(). */
    pub usercontrol: Option<Box<crate::arkode_user_controller::ARKUserControlContent>>,

    pub expstab: Option<ARKExpStabFn>, /* step stability function        */
    pub estab_data: UserData,          /* user pointer passed to expstab */

    pub nst_acc: i64, /* num accuracy-limited internal steps  */
    pub nst_exp: i64, /* num stability-limited internal steps */
}
