/* -----------------------------------------------------------------
 * Translation of sundials-7.7.0/src/arkode/arkode_root.c.
 *
 * Adaptations per the crate conventions (ARCHITECTURE Addendum C.1):
 *  - `void* arkode_mem` becomes &mut ARKodeMem (the C NULL check has
 *    no translation).
 *  - The check/find routines take root_mem out of ark_mem for the
 *    call and put it back (public wrappers + _inner workers), so the
 *    root arrays and the ark_mem vectors borrow disjointly.
 *  - C caches root_mem->root_data = ark_mem->user_data (a pointer
 *    alias that ARKodeSetUserData keeps in sync); the Rust port
 *    passes &mut ark_mem.user_data to gfun directly and leaves the
 *    root_data field unused.
 *  - ARKodeGetDky writes into ark_mem-owned vectors (tempv4/ycur);
 *    those call sites take the destination vector out around the
 *    call.
 * -----------------------------------------------------------------*/

use crate::arkode::{arkAllocVec, ARKodeGetDky};
use crate::arkode_impl::{
    ark_step_fullrhs_yn_fn, arkProcessError, ARKRootFn, ARKodeMem, ARK_FULLRHS_START,
    ARK_ILL_INPUT, ARK_NORMAL, ARK_ONE_STEP, ARK_RHSFUNC_FAIL, ARK_RTFUNC_FAIL,
    ARK_SUCCESS, CLOSERT, FIVE, HALF, MSG_ARK_MISSING_FULLRHS, MSG_ARK_NULL_G, ONE, RTFOUND,
    TENTH, TWO, ZERO,
};
use crate::arkode_root_impl::{ARKodeRootMem, ARK_ROOT_LIW, ARK_ROOT_LRW, HUND};
use crate::nvector_serial::{N_VLinearSum, N_VScale, NVector};
use crate::sundials_math::{SUNMAX, SUNRabs, SUNRdifferentsign};
use crate::sundials_utils::fmt_g;

/*===============================================================
  Exported functions
  ===============================================================*/

/*---------------------------------------------------------------
  ARKodeRootInit:

  ARKodeRootInit initializes a rootfinding problem to be solved
  during the integration of the ODE system.  It loads the root
  function pointer and the number of root functions, notifies
  ARKODE that the "fullrhs" function is required, and allocates
  workspace memory.  The return value is ARK_SUCCESS = 0 if no
  errors occurred, or a negative value otherwise.
  ---------------------------------------------------------------*/
/* C compares the incoming g against the stored gfun pointer; Rust fn
   addresses are not guaranteed unique, but a spurious mismatch only
   re-assigns the same function (harmless) */
#[allow(unpredictable_function_pointer_comparisons)]
pub fn ARKodeRootInit(ark_mem: &mut ARKodeMem, nrtfn: i32, g: Option<ARKRootFn>) -> i32 {
    let nrt = if nrtfn < 0 { 0 } else { nrtfn };

    /* Ensure that stepper provides fullrhs function */
    if nrt > 0 {
        if ark_mem.step_fullrhs.is_none() {
            arkProcessError(
                Some(ark_mem),
                ARK_ILL_INPUT,
                line!(),
                "ARKodeRootInit",
                file!(),
                MSG_ARK_MISSING_FULLRHS,
            );
            return ARK_ILL_INPUT;
        }

        let yn_len = ark_mem.yn.data.len();
        let mut fn_ = std::mem::replace(&mut ark_mem.fn_, NVector::new(0));
        arkAllocVec(ark_mem, yn_len, &mut fn_);
        ark_mem.fn_ = fn_;
    }

    /* If unallocated, allocate rootfinding structure, set defaults, update
    space */
    if ark_mem.root_mem.is_none() {
        ark_mem.root_mem = Some(ARKodeRootMem {
            glo: Vec::new(),
            ghi: Vec::new(),
            grout: Vec::new(),
            iroots: Vec::new(),
            rootdir: Vec::new(),
            gfun: None,
            nrtfn: 0,
            irfnd: 0,
            gactive: Vec::new(),
            mxgnull: 1,
            /* C: root_data = ark_mem->user_data (alias; unused here) */
            root_data: None,
            tlo: 0.0,
            thi: 0.0,
            trout: 0.0,
            ttol: 0.0,
            nge: 0,
        });

        ark_mem.lrw += ARK_ROOT_LRW;
        ark_mem.liw += ARK_ROOT_LIW;
    }

    let rootmem = ark_mem.root_mem.as_mut().unwrap();

    /* If rerunning ARKodeRootInit() with a different number of root
    functions (changing number of gfun components), then free
    currently held memory resources */
    if (nrt != rootmem.nrtfn) && (rootmem.nrtfn > 0) {
        rootmem.glo = Vec::new();
        rootmem.ghi = Vec::new();
        rootmem.grout = Vec::new();
        rootmem.iroots = Vec::new();
        rootmem.rootdir = Vec::new();
        rootmem.gactive = Vec::new();

        let old = rootmem.nrtfn as i64;
        ark_mem.lrw -= 3 * old;
        ark_mem.liw -= 3 * old;
    }

    let rootmem = ark_mem.root_mem.as_mut().unwrap();

    /* If ARKodeRootInit() was called with nrtfn == 0, then set
    nrtfn to zero and gfun to NULL before returning */
    if nrt == 0 {
        rootmem.nrtfn = nrt;
        rootmem.gfun = None;
        return ARK_SUCCESS;
    }

    /* If rerunning ARKodeRootInit() with the same number of root
    functions (not changing number of gfun components), then
    check if the root function argument has changed */
    /* If g != NULL then return as currently reserved memory
    resources will suffice */
    if nrt == rootmem.nrtfn {
        if g != rootmem.gfun {
            if g.is_none() {
                rootmem.glo = Vec::new();
                rootmem.ghi = Vec::new();
                rootmem.grout = Vec::new();
                rootmem.iroots = Vec::new();
                rootmem.rootdir = Vec::new();
                rootmem.gactive = Vec::new();

                ark_mem.lrw -= 3 * nrt as i64;
                ark_mem.liw -= 3 * nrt as i64;

                arkProcessError(
                    Some(ark_mem),
                    ARK_ILL_INPUT,
                    line!(),
                    "ARKodeRootInit",
                    file!(),
                    MSG_ARK_NULL_G,
                );
                return ARK_ILL_INPUT;
            } else {
                rootmem.gfun = g;
                return ARK_SUCCESS;
            }
        } else {
            return ARK_SUCCESS;
        }
    }

    /* Set variable values in ARKODE memory block */
    rootmem.nrtfn = nrt;
    if g.is_none() {
        arkProcessError(
            Some(ark_mem),
            ARK_ILL_INPUT,
            line!(),
            "ARKodeRootInit",
            file!(),
            MSG_ARK_NULL_G,
        );
        return ARK_ILL_INPUT;
    } else {
        rootmem.gfun = g;
    }

    /* Allocate necessary memory and return (allocation cannot fail) */
    rootmem.glo = vec![0.0; nrt as usize];
    rootmem.ghi = vec![0.0; nrt as usize];
    rootmem.grout = vec![0.0; nrt as usize];
    rootmem.iroots = vec![0; nrt as usize];
    rootmem.rootdir = vec![0; nrt as usize];
    rootmem.gactive = Vec::new();

    /* Set default values for rootdir (both directions) */
    for i in 0..nrt as usize {
        rootmem.rootdir[i] = 0;
    }

    /* Set default values for gactive (all active) */
    rootmem.gactive = vec![true; nrt as usize];

    ark_mem.lrw += 3 * nrt as i64;
    ark_mem.liw += 3 * nrt as i64;

    ARK_SUCCESS
}

/*===============================================================
  Private functions
  ===============================================================*/

/*---------------------------------------------------------------
  arkRootFree

  This routine frees all memory associated with ARKODE's
  rootfinding module.
  ---------------------------------------------------------------*/
pub fn arkRootFree(ark_mem: &mut ARKodeMem) -> i32 {
    if let Some(rootmem) = ark_mem.root_mem.take() {
        if rootmem.nrtfn > 0 {
            ark_mem.lrw -= 3 * rootmem.nrtfn as i64;
            ark_mem.liw -= 3 * rootmem.nrtfn as i64;
        }
        drop(rootmem);
        ark_mem.lrw -= ARK_ROOT_LRW;
        ark_mem.liw -= ARK_ROOT_LIW;
    }
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkPrintRootMem

  This routine outputs the root-finding memory structure to a
  specified file pointer.
  ---------------------------------------------------------------*/
pub fn arkPrintRootMem(ark_mem: &ARKodeMem, outfile: &mut dyn std::io::Write) -> i32 {
    if let Some(rootmem) = ark_mem.root_mem.as_ref() {
        let _ = write!(outfile, "ark_nrtfn = {}\n", rootmem.nrtfn);
        let _ = write!(outfile, "ark_nge = {}\n", rootmem.nge);
        if !rootmem.iroots.is_empty() {
            for i in 0..rootmem.nrtfn as usize {
                let _ = write!(outfile, "ark_iroots[{}] = {}\n", i, rootmem.iroots[i]);
            }
        }
        if !rootmem.rootdir.is_empty() {
            for i in 0..rootmem.nrtfn as usize {
                let _ = write!(outfile, "ark_rootdir[{}] = {}\n", i, rootmem.rootdir[i]);
            }
        }
        let _ = write!(outfile, "ark_irfnd = {}\n", rootmem.irfnd);
        let _ = write!(outfile, "ark_mxgnull = {}\n", rootmem.mxgnull);
        if !rootmem.gactive.is_empty() {
            for i in 0..rootmem.nrtfn as usize {
                let _ = write!(outfile, "ark_gactive[{}] = {}\n", i, rootmem.gactive[i] as i32);
            }
        }
        let _ = write!(outfile, "ark_tlo = {}\n", fmt_g(rootmem.tlo, 0, 15));
        let _ = write!(outfile, "ark_thi = {}\n", fmt_g(rootmem.thi, 0, 15));
        let _ = write!(outfile, "ark_trout = {}\n", fmt_g(rootmem.trout, 0, 15));
        if !rootmem.glo.is_empty() {
            for i in 0..rootmem.nrtfn as usize {
                let _ = write!(outfile, "ark_glo[{}] = {}\n", i, fmt_g(rootmem.glo[i], 0, 15));
            }
        }
        if !rootmem.ghi.is_empty() {
            for i in 0..rootmem.nrtfn as usize {
                let _ = write!(outfile, "ark_ghi[{}] = {}\n", i, fmt_g(rootmem.ghi[i], 0, 15));
            }
        }
        if !rootmem.grout.is_empty() {
            for i in 0..rootmem.nrtfn as usize {
                let _ = write!(
                    outfile,
                    "ark_grout[{}] = {}\n",
                    i,
                    fmt_g(rootmem.grout[i], 0, 15)
                );
            }
        }
        let _ = write!(outfile, "ark_ttol = {}\n", fmt_g(rootmem.ttol, 0, 15));
    }
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkRootCheck1

  This routine completes the initialization of rootfinding memory
  information, and checks whether g has a zero both at and very near
  the initial point of the IVP.

  This routine returns an int equal to:
    ARK_RTFUNC_FAIL < 0  if the g function failed, or
    ARK_SUCCESS     = 0  otherwise.
  ---------------------------------------------------------------*/
pub fn arkRootCheck1(ark_mem: &mut ARKodeMem) -> i32 {
    let mut rootmem = ark_mem.root_mem.take().unwrap();
    let ret = arkRootCheck1_inner(ark_mem, &mut rootmem);
    ark_mem.root_mem = Some(rootmem);
    ret
}

fn arkRootCheck1_inner(ark_mem: &mut ARKodeMem, rootmem: &mut ARKodeRootMem) -> i32 {
    for i in 0..rootmem.nrtfn as usize {
        rootmem.iroots[i] = 0;
    }
    rootmem.tlo = ark_mem.tcur;
    rootmem.ttol = (SUNRabs(ark_mem.tcur) + SUNRabs(ark_mem.h)) * ark_mem.uround * HUND;

    /* Evaluate g at initial t and check for zero values. */
    let gfun = rootmem.gfun.unwrap();
    let retval = gfun(
        rootmem.tlo,
        &ark_mem.yn,
        &mut rootmem.glo,
        &mut ark_mem.user_data,
    );
    rootmem.nge = 1;
    if retval != 0 {
        arkProcessError(
            Some(ark_mem),
            ARK_RTFUNC_FAIL,
            line!(),
            "arkRootCheck1",
            file!(),
            &format!(
                "At t = {}, the rootfinding routine failed in an unrecoverable manner.",
                fmt_g(ark_mem.tcur, 0, 15)
            ),
        );
        return ARK_RTFUNC_FAIL;
    }

    let mut zroot = false;
    for i in 0..rootmem.nrtfn as usize {
        if SUNRabs(rootmem.glo[i]) == ZERO {
            zroot = true;
            rootmem.gactive[i] = false;
        }
    }
    if !zroot {
        return ARK_SUCCESS;
    }

    /* call full RHS if needed */
    if !ark_mem.fn_is_current {
        let retval = ark_step_fullrhs_yn_fn(ark_mem, ark_mem.tn, ARK_FULLRHS_START);
        if retval != 0 {
            arkProcessError(
                Some(ark_mem),
                ARK_RHSFUNC_FAIL,
                line!(),
                "arkRootCheck1",
                file!(),
                &format!(
                    "At t = {}, the right-hand side routine failed in an unrecoverable manner.",
                    fmt_g(ark_mem.tcur, 0, 15)
                ),
            );
            return ARK_RHSFUNC_FAIL;
        }
        ark_mem.fn_is_current = true;
    }

    /* Some g_i is zero at t0; look at g at t0+(small increment). */
    let hratio = SUNMAX(rootmem.ttol / SUNRabs(ark_mem.h), TENTH);
    let smallh = hratio * ark_mem.h;
    let tplus = rootmem.tlo + smallh;
    N_VLinearSum(ONE, &ark_mem.yn, smallh, &ark_mem.fn_, &mut ark_mem.tempv4);
    let retval = gfun(tplus, &ark_mem.tempv4, &mut rootmem.ghi, &mut ark_mem.user_data);
    rootmem.nge += 1;
    if retval != 0 {
        arkProcessError(
            Some(ark_mem),
            ARK_RTFUNC_FAIL,
            line!(),
            "arkRootCheck1",
            file!(),
            &format!(
                "At t = {}, the rootfinding routine failed in an unrecoverable manner.",
                fmt_g(ark_mem.tcur, 0, 15)
            ),
        );
        return ARK_RTFUNC_FAIL;
    }

    /* We check now only the components of g which were exactly 0.0 at t0
     * to see if we can 'activate' them. */
    for i in 0..rootmem.nrtfn as usize {
        if !rootmem.gactive[i] && SUNRabs(rootmem.ghi[i]) != ZERO {
            rootmem.gactive[i] = true;
            rootmem.glo[i] = rootmem.ghi[i];
        }
    }
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkRootCheck2

  This routine checks for exact zeros of g at the last root found,
  if the last return was a root.  It then checks for a close pair of
  zeros (an error condition), and for a new root at a nearby point.
  The array glo = g(tlo) at the left endpoint of the search interval
  is adjusted if necessary to assure that all g_i are nonzero
  there, before returning to do a root search in the interval.

  On entry, tlo = tretlast is the last value of tret returned by
  ARKODE.  This may be the previous tn, the previous tout value, or
  the last root location.

  This routine returns an int equal to:
    ARK_RTFUNC_FAIL < 0 if the g function failed, or
    CLOSERT         = 3 if a close pair of zeros was found, or
    RTFOUND         = 1 if a new zero of g was found near tlo, or
    ARK_SUCCESS     = 0 otherwise.
  ---------------------------------------------------------------*/
pub fn arkRootCheck2(ark_mem: &mut ARKodeMem) -> i32 {
    let mut rootmem = ark_mem.root_mem.take().unwrap();
    let ret = arkRootCheck2_inner(ark_mem, &mut rootmem);
    ark_mem.root_mem = Some(rootmem);
    ret
}

fn arkRootCheck2_inner(ark_mem: &mut ARKodeMem, rootmem: &mut ARKodeRootMem) -> i32 {
    /* return if no roots in previous step */
    if rootmem.irfnd == 0 {
        return ARK_SUCCESS;
    }

    /* Set tempv4 = y(tlo) */
    let mut tempv4 = std::mem::replace(&mut ark_mem.tempv4, NVector::new(0));
    let _ = ARKodeGetDky(ark_mem, rootmem.tlo, 0, &mut tempv4);
    ark_mem.tempv4 = tempv4;

    /* Evaluate root-finding function: glo = g(tlo, y(tlo)) */
    let gfun = rootmem.gfun.unwrap();
    let retval = gfun(
        rootmem.tlo,
        &ark_mem.tempv4,
        &mut rootmem.glo,
        &mut ark_mem.user_data,
    );
    rootmem.nge += 1;
    if retval != 0 {
        return ARK_RTFUNC_FAIL;
    }

    /* reset root-finding flags (overall, and for specific eqns) */
    let mut zroot = false;
    for i in 0..rootmem.nrtfn as usize {
        rootmem.iroots[i] = 0;
    }

    /* for all active roots, check if glo_i == 0 to mark roots found */
    for i in 0..rootmem.nrtfn as usize {
        if !rootmem.gactive[i] {
            continue;
        }
        if SUNRabs(rootmem.glo[i]) == ZERO {
            zroot = true;
            rootmem.iroots[i] = 1;
        }
    }
    if !zroot {
        return ARK_SUCCESS; /* return if no roots */
    }

    /* One or more g_i has a zero at tlo.  Check g at tlo+smallh. */
    /*     set time tolerance */
    rootmem.ttol = (SUNRabs(ark_mem.tcur) + SUNRabs(ark_mem.h)) * ark_mem.uround * HUND;
    /*     set tplus = tlo + smallh */
    let smallh = if ark_mem.h > ZERO {
        rootmem.ttol
    } else {
        -rootmem.ttol
    };
    let tplus = rootmem.tlo + smallh;
    /*     update ark_ycur with small explicit Euler step (if tplus is past
    tn) */
    if (tplus - ark_mem.tcur) * ark_mem.h >= ZERO {
        /* hratio = smallh/ark_mem->h; */
        N_VLinearSum(ONE, &ark_mem.tempv4, smallh, &ark_mem.fn_, &mut ark_mem.ycur);
    } else {
        /*   set ark_ycur = y(tplus) via interpolation */
        let mut ycur = std::mem::replace(&mut ark_mem.ycur, NVector::new(0));
        let _ = ARKodeGetDky(ark_mem, tplus, 0, &mut ycur);
        ark_mem.ycur = ycur;
    }
    /*     set ghi = g(tplus,y(tplus)) */
    let retval = gfun(tplus, &ark_mem.ycur, &mut rootmem.ghi, &mut ark_mem.user_data);
    rootmem.nge += 1;
    if retval != 0 {
        return ARK_RTFUNC_FAIL;
    }

    /* Check for close roots (error return), for a new zero at tlo+smallh,
    and for a g_i that changed from zero to nonzero. */
    let mut zroot = false;
    for i in 0..rootmem.nrtfn as usize {
        if !rootmem.gactive[i] {
            continue;
        }
        if SUNRabs(rootmem.ghi[i]) == ZERO {
            if rootmem.iroots[i] == 1 {
                return CLOSERT;
            }
            zroot = true;
            rootmem.iroots[i] = 1;
        } else {
            if rootmem.iroots[i] == 1 {
                rootmem.glo[i] = rootmem.ghi[i];
            }
        }
    }
    if zroot {
        return RTFOUND;
    }
    ARK_SUCCESS
}

/*---------------------------------------------------------------
  arkRootCheck3

  This routine interfaces to arkRootfind to look for a root of g
  between tlo and either tn or tout, whichever comes first.
  Only roots beyond tlo in the direction of integration are sought.

  This routine returns an int equal to:
    ARK_RTFUNC_FAIL < 0 if the g function failed, or
    RTFOUND         = 1 if a root of g was found, or
    ARK_SUCCESS     = 0 otherwise.
  ---------------------------------------------------------------*/
pub fn arkRootCheck3(ark_mem: &mut ARKodeMem, tout: f64, itask: i32) -> i32 {
    let mut rootmem = ark_mem.root_mem.take().unwrap();
    let ret = arkRootCheck3_inner(ark_mem, &mut rootmem, tout, itask);
    ark_mem.root_mem = Some(rootmem);
    ret
}

fn arkRootCheck3_inner(
    ark_mem: &mut ARKodeMem,
    rootmem: &mut ARKodeRootMem,
    tout: f64,
    itask: i32,
) -> i32 {
    /* Set thi = tn or tout, whichever comes first; set y = y(thi). */
    if itask == ARK_ONE_STEP {
        rootmem.thi = ark_mem.tcur;
        N_VScale(ONE, &ark_mem.yn, &mut ark_mem.tempv4);
    }
    if itask == ARK_NORMAL {
        if (tout - ark_mem.tcur) * ark_mem.h >= ZERO {
            rootmem.thi = ark_mem.tcur;
            N_VScale(ONE, &ark_mem.yn, &mut ark_mem.tempv4);
        } else {
            rootmem.thi = tout;
            let mut tempv4 = std::mem::replace(&mut ark_mem.tempv4, NVector::new(0));
            let _ = ARKodeGetDky(ark_mem, rootmem.thi, 0, &mut tempv4);
            ark_mem.tempv4 = tempv4;
        }
    }

    /* Set rootmem->ghi = g(thi) and call arkRootfind to search (tlo,thi) for
    roots. */
    let gfun = rootmem.gfun.unwrap();
    let retval = gfun(
        rootmem.thi,
        &ark_mem.tempv4,
        &mut rootmem.ghi,
        &mut ark_mem.user_data,
    );
    rootmem.nge += 1;
    if retval != 0 {
        return ARK_RTFUNC_FAIL;
    }

    rootmem.ttol = (SUNRabs(ark_mem.tcur) + SUNRabs(ark_mem.h)) * ark_mem.uround * HUND;
    let ier = arkRootfind_inner(ark_mem, rootmem);
    if ier == ARK_RTFUNC_FAIL {
        return ARK_RTFUNC_FAIL;
    }
    for i in 0..rootmem.nrtfn as usize {
        if !rootmem.gactive[i] && rootmem.grout[i] != ZERO {
            rootmem.gactive[i] = true;
        }
    }
    rootmem.tlo = rootmem.trout;
    for i in 0..rootmem.nrtfn as usize {
        rootmem.glo[i] = rootmem.grout[i];
    }

    /* If no root found, return ARK_SUCCESS. */
    if ier == ARK_SUCCESS {
        return ARK_SUCCESS;
    }

    /* If a root was found, interpolate to get y(trout) and return.  */
    let mut ycur = std::mem::replace(&mut ark_mem.ycur, NVector::new(0));
    let _ = ARKodeGetDky(ark_mem, rootmem.trout, 0, &mut ycur);
    ark_mem.ycur = ycur;
    RTFOUND
}

/*---------------------------------------------------------------
  arkRootfind

  This routine solves for a root of g(t) between tlo and thi, if
  one exists.  Only roots of odd multiplicity (i.e. with a change
  of sign in one of the g_i), or exact zeros, are found.
  Here the sign of tlo - thi is arbitrary, but if multiple roots
  are found, the one closest to tlo is returned.

  The method used is the Illinois algorithm, a modified secant
  method (see the C source for the full parameter description).

  This routine returns an int equal to:
    ARK_RTFUNC_FAIL < 0 if the g function failed, or
    RTFOUND         = 1 if a root of g was found, or
    ARK_SUCCESS     = 0 otherwise.
  ---------------------------------------------------------------*/
pub fn arkRootfind(ark_mem: &mut ARKodeMem) -> i32 {
    let mut rootmem = ark_mem.root_mem.take().unwrap();
    let ret = arkRootfind_inner(ark_mem, &mut rootmem);
    ark_mem.root_mem = Some(rootmem);
    ret
}

fn arkRootfind_inner(ark_mem: &mut ARKodeMem, rootmem: &mut ARKodeRootMem) -> i32 {
    let mut imax = 0usize;

    /* First check for change in sign in ghi or for a zero in ghi. */
    let mut maxfrac = ZERO;
    let mut zroot = false;
    let mut sgnchg = false;
    for i in 0..rootmem.nrtfn as usize {
        if !rootmem.gactive[i] {
            continue;
        }
        if SUNRabs(rootmem.ghi[i]) == ZERO {
            if rootmem.rootdir[i] as f64 * rootmem.glo[i] <= ZERO {
                zroot = true;
            }
        } else {
            if SUNRdifferentsign(rootmem.glo[i], rootmem.ghi[i])
                && (rootmem.rootdir[i] as f64 * rootmem.glo[i] <= ZERO)
            {
                let gfrac = SUNRabs(rootmem.ghi[i] / (rootmem.ghi[i] - rootmem.glo[i]));
                if gfrac > maxfrac {
                    sgnchg = true;
                    maxfrac = gfrac;
                    imax = i;
                }
            }
        }
    }

    /* If no sign change was found, reset trout and grout.  Then return
    ARK_SUCCESS if no zero was found, or set iroots and return RTFOUND.  */
    if !sgnchg {
        rootmem.trout = rootmem.thi;
        for i in 0..rootmem.nrtfn as usize {
            rootmem.grout[i] = rootmem.ghi[i];
        }
        if !zroot {
            return ARK_SUCCESS;
        }
        for i in 0..rootmem.nrtfn as usize {
            rootmem.iroots[i] = 0;
            if !rootmem.gactive[i] {
                continue;
            }
            if SUNRabs(rootmem.ghi[i]) == ZERO {
                rootmem.iroots[i] = if rootmem.glo[i] > 0.0 { -1 } else { 1 };
            }
        }
        return RTFOUND;
    }

    /* Initialize alpha to avoid compiler warning */
    let mut alpha = ONE;

    /* A sign change was found.  Loop to locate nearest root. */
    let mut side = 0;
    let mut sideprev = -1;
    loop {
        /* Looping point */

        /* If interval size is already less than tolerance ttol, break. */
        if SUNRabs(rootmem.thi - rootmem.tlo) <= rootmem.ttol {
            break;
        }

        /* Set weight alpha.
        On the first two passes, set alpha = 1.  Thereafter, reset alpha
        according to the side (low vs high) of the subinterval in which
        the sign change was found in the previous two passes.
        If the sides were opposite, set alpha = 1.
        If the sides were the same, then double alpha (if high side),
        or halve alpha (if low side).
        The next guess tmid is the secant method value if alpha = 1, but
        is closer to tlo if alpha < 1, and closer to thi if alpha > 1.    */
        if sideprev == side {
            alpha = if side == 2 { alpha * TWO } else { alpha * HALF };
        } else {
            alpha = ONE;
        }

        /* Set next root approximation tmid and get g(tmid).
        If tmid is too close to tlo or thi, adjust it inward,
        by a fractional distance that is between 0.1 and 0.5.  */
        let mut tmid = rootmem.thi
            - (rootmem.thi - rootmem.tlo) * rootmem.ghi[imax]
                / (rootmem.ghi[imax] - alpha * rootmem.glo[imax]);
        if SUNRabs(tmid - rootmem.tlo) < HALF * rootmem.ttol {
            let fracint = SUNRabs(rootmem.thi - rootmem.tlo) / rootmem.ttol;
            let fracsub = if fracint > FIVE { TENTH } else { HALF / fracint };
            tmid = rootmem.tlo + fracsub * (rootmem.thi - rootmem.tlo);
        }
        if SUNRabs(rootmem.thi - tmid) < HALF * rootmem.ttol {
            let fracint = SUNRabs(rootmem.thi - rootmem.tlo) / rootmem.ttol;
            let fracsub = if fracint > FIVE { TENTH } else { HALF / fracint };
            tmid = rootmem.thi - fracsub * (rootmem.thi - rootmem.tlo);
        }

        let mut tempv4 = std::mem::replace(&mut ark_mem.tempv4, NVector::new(0));
        let _ = ARKodeGetDky(ark_mem, tmid, 0, &mut tempv4);
        ark_mem.tempv4 = tempv4;
        let gfun = rootmem.gfun.unwrap();
        let retval = gfun(
            tmid,
            &ark_mem.tempv4,
            &mut rootmem.grout,
            &mut ark_mem.user_data,
        );
        rootmem.nge += 1;
        if retval != 0 {
            return ARK_RTFUNC_FAIL;
        }

        /* Check to see in which subinterval g changes sign, and reset imax.
        Set side = 1 if sign change is on low side, or 2 if on high side.  */
        maxfrac = ZERO;
        zroot = false;
        sgnchg = false;
        sideprev = side;
        for i in 0..rootmem.nrtfn as usize {
            if !rootmem.gactive[i] {
                continue;
            }
            if SUNRabs(rootmem.grout[i]) == ZERO {
                if rootmem.rootdir[i] as f64 * rootmem.glo[i] <= ZERO {
                    zroot = true;
                }
            } else {
                if SUNRdifferentsign(rootmem.glo[i], rootmem.grout[i])
                    && (rootmem.rootdir[i] as f64 * rootmem.glo[i] <= ZERO)
                {
                    let gfrac = SUNRabs(rootmem.grout[i] / (rootmem.grout[i] - rootmem.glo[i]));
                    if gfrac > maxfrac {
                        sgnchg = true;
                        maxfrac = gfrac;
                        imax = i;
                    }
                }
            }
        }
        if sgnchg {
            /* Sign change found in (tlo,tmid); replace thi with tmid. */
            rootmem.thi = tmid;
            for i in 0..rootmem.nrtfn as usize {
                rootmem.ghi[i] = rootmem.grout[i];
            }
            side = 1;
            /* Stop at root thi if converged; otherwise loop. */
            if SUNRabs(rootmem.thi - rootmem.tlo) <= rootmem.ttol {
                break;
            }
            continue; /* Return to looping point. */
        }

        if zroot {
            /* No sign change in (tlo,tmid), but g = 0 at tmid; return root
            tmid. */
            rootmem.thi = tmid;
            for i in 0..rootmem.nrtfn as usize {
                rootmem.ghi[i] = rootmem.grout[i];
            }
            break;
        }

        /* No sign change in (tlo,tmid), and no zero at tmid.
        Sign change must be in (tmid,thi).  Replace tlo with tmid. */
        rootmem.tlo = tmid;
        for i in 0..rootmem.nrtfn as usize {
            rootmem.glo[i] = rootmem.grout[i];
        }
        side = 2;
        /* Stop at root thi if converged; otherwise loop back. */
        if SUNRabs(rootmem.thi - rootmem.tlo) <= rootmem.ttol {
            break;
        }
    } /* End of root-search loop */

    /* Reset trout and grout, set iroots, and return RTFOUND. */
    rootmem.trout = rootmem.thi;
    for i in 0..rootmem.nrtfn as usize {
        rootmem.grout[i] = rootmem.ghi[i];
        rootmem.iroots[i] = 0;
        if !rootmem.gactive[i] {
            continue;
        }
        if (SUNRabs(rootmem.ghi[i]) == ZERO)
            && (rootmem.rootdir[i] as f64 * rootmem.glo[i] <= ZERO)
        {
            rootmem.iroots[i] = if rootmem.glo[i] > 0.0 { -1 } else { 1 };
        }
        if SUNRdifferentsign(rootmem.glo[i], rootmem.ghi[i])
            && (rootmem.rootdir[i] as f64 * rootmem.glo[i] <= ZERO)
        {
            rootmem.iroots[i] = if rootmem.glo[i] > 0.0 { -1 } else { 1 };
        }
    }
    RTFOUND
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arkode_interp::{arkInterpCreate_Hermite, arkInterpInit, arkInterpUpdate};
    use crate::sundials_types::SUN_UNIT_ROUNDOFF;

    /* full RHS for y' = 3t^2 (y = t^3), independent of y */
    fn cubic_rhs(_a: &mut ARKodeMem, t: f64, _y: &NVector, f: &mut NVector, _m: i32) -> i32 {
        f.data[0] = 3.0 * t * t;
        0
    }

    /* g(t,y) = y - 3.375: root at t = 1.5 on y = t^3, increasing */
    fn g_shift(_t: f64, y: &NVector, gout: &mut [f64], _ud: &mut crate::sundials_types::UserData) -> i32 {
        gout[0] = y.data[0] - 3.375;
        0
    }

    /* g(t,y) = y - 1: exactly zero at the initial point (t=1, y=1) */
    fn g_zero_at_t0(_t: f64, y: &NVector, gout: &mut [f64], _ud: &mut crate::sundials_types::UserData) -> i32 {
        gout[0] = y.data[0] - 1.0;
        0
    }

    /* ark_mem holding y = t^3 state at t = 1 with a cubic Hermite
    interpolant, mirroring the arkode_interp tests */
    fn mem_at_t1() -> ARKodeMem {
        let mut ark_mem = ARKodeMem::default();
        ark_mem.uround = SUN_UNIT_ROUNDOFF;
        ark_mem.lrw1 = 1;
        ark_mem.liw1 = 2;
        ark_mem.step_fullrhs = Some(cubic_rhs);
        ark_mem.tn = 1.0;
        ark_mem.tcur = 1.0;
        ark_mem.h = 1.0;
        ark_mem.yn = NVector::new(1);
        ark_mem.yn.data[0] = 1.0;
        ark_mem.fn_ = NVector::new(1);
        ark_mem.tempv4 = NVector::new(1);
        ark_mem.ycur = NVector::new(1);
        ark_mem.interp = arkInterpCreate_Hermite(&mut ark_mem, 3);
        assert_eq!(arkInterpInit(&mut ark_mem, 1.0), ARK_SUCCESS);
        ark_mem
    }

    /* advance the manufactured step from t=1 to t=2 (y: 1 -> 8) */
    fn take_step(ark_mem: &mut ARKodeMem) {
        assert_eq!(arkInterpUpdate(ark_mem, 2.0), ARK_SUCCESS);
        ark_mem.tn = 2.0;
        ark_mem.tcur = 2.0;
        ark_mem.hold = 1.0;
        ark_mem.yn.data[0] = 8.0;
        ark_mem.fn_is_current = false;
    }

    #[test]
    fn root_init_alloc_and_rerun_paths() {
        let mut ark_mem = mem_at_t1();
        assert_eq!(ARKodeRootInit(&mut ark_mem, 1, Some(g_shift)), ARK_SUCCESS);
        let rm = ark_mem.root_mem.as_ref().unwrap();
        assert_eq!(rm.nrtfn, 1);
        assert_eq!(rm.rootdir, vec![0]);
        assert_eq!(rm.gactive, vec![true]);

        /* same count, different g: just swaps the function */
        assert_eq!(ARKodeRootInit(&mut ark_mem, 1, Some(g_zero_at_t0)), ARK_SUCCESS);
        #[allow(unpredictable_function_pointer_comparisons)]
        {
            assert_eq!(
                ark_mem.root_mem.as_ref().unwrap().gfun,
                Some(g_zero_at_t0 as ARKRootFn)
            );
        }

        /* same count, NULL g: error */
        assert_eq!(ARKodeRootInit(&mut ark_mem, 1, None), ARK_ILL_INPUT);

        /* nrtfn = 0 disables rootfinding */
        assert_eq!(ARKodeRootInit(&mut ark_mem, 0, None), ARK_SUCCESS);
        assert_eq!(ark_mem.root_mem.as_ref().unwrap().nrtfn, 0);

        /* free restores the workspace counters */
        assert_eq!(arkRootFree(&mut ark_mem), ARK_SUCCESS);
        assert!(ark_mem.root_mem.is_none());
    }

    #[test]
    fn rootfind_locates_interior_root() {
        let mut ark_mem = mem_at_t1();
        assert_eq!(ARKodeRootInit(&mut ark_mem, 1, Some(g_shift)), ARK_SUCCESS);
        assert_eq!(arkRootCheck1(&mut ark_mem), ARK_SUCCESS);
        {
            let rm = ark_mem.root_mem.as_ref().unwrap();
            assert_eq!(rm.nge, 1);
            assert!((rm.glo[0] + 2.375).abs() < 1e-15);
        }

        take_step(&mut ark_mem);
        assert_eq!(arkRootCheck3(&mut ark_mem, 0.0, ARK_ONE_STEP), RTFOUND);
        let rm = ark_mem.root_mem.as_ref().unwrap();
        assert!(
            (rm.trout - 1.5).abs() <= rm.ttol,
            "trout = {}, ttol = {}",
            rm.trout,
            rm.ttol
        );
        assert_eq!(rm.iroots[0], 1); /* increasing crossing */
        assert!((ark_mem.ycur.data[0] - 3.375).abs() < 1e-6);
    }

    #[test]
    fn rootdir_filters_direction() {
        let mut ark_mem = mem_at_t1();
        assert_eq!(ARKodeRootInit(&mut ark_mem, 1, Some(g_shift)), ARK_SUCCESS);
        ark_mem.root_mem.as_mut().unwrap().rootdir[0] = -1; /* only decreasing */
        assert_eq!(arkRootCheck1(&mut ark_mem), ARK_SUCCESS);
        take_step(&mut ark_mem);
        assert_eq!(arkRootCheck3(&mut ark_mem, 0.0, ARK_ONE_STEP), ARK_SUCCESS);
        assert_eq!(ark_mem.root_mem.as_ref().unwrap().iroots[0], 0);
    }

    #[test]
    fn check1_reactivates_zero_at_t0() {
        let mut ark_mem = mem_at_t1();
        assert_eq!(ARKodeRootInit(&mut ark_mem, 1, Some(g_zero_at_t0)), ARK_SUCCESS);
        assert_eq!(arkRootCheck1(&mut ark_mem), ARK_SUCCESS);
        let rm = ark_mem.root_mem.as_ref().unwrap();
        /* g was exactly zero at t0 but moved away at t0+smallh, so the
        component is re-activated with glo = ghi */
        assert_eq!(rm.nge, 2);
        assert!(rm.gactive[0]);
        assert!(rm.glo[0] > 0.0);
    }
}
