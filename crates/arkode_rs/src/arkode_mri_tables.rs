/* -----------------------------------------------------------------
 * Translation of src/arkode/arkode_mri_tables.c (+ the
 * arkode_mri_tables.def X-macro table, arkode_mri_tables_impl.h and
 * the coupling-table parts of include/arkode/arkode_mristep.h),
 * SUNDIALS 7.7.0.
 *
 * This is the implementation file for ARKODE's MRIStepCoupling
 * tables.
 *
 * C's MRIStepCoupling is a pointer to a Mem struct whose W/G
 * tensors are ragged calloc'd arrays indexed [k][i][j] (matrix,
 * row 0..=stages, column 0..stages); the Rust port models them as
 * nested Vecs with the same indexing, where an empty outer Vec
 * plays the role of a NULL W/G pointer.  Functions that return
 * NULL on failure return Option<MRIStepCoupling> (a boxed Mem).
 * The C field `type` is renamed `type_`.
 * -----------------------------------------------------------------*/

use crate::arkode_butcher::{ARKodeButcherTable, ARKodeButcherTable_Alloc};
use crate::arkode_butcher_erk::{
    ARKodeButcherTable_LoadERK, ARKODE_EXPLICIT_MIDPOINT_EULER_2_1_2, ARKODE_FORWARD_EULER_1_1,
    ARKODE_HEUN_EULER_2_1_2, ARKODE_KNOTH_WOLKE_3_3, ARKODE_RALSTON_EULER_2_1_2,
};
use crate::arkode_impl::{arkProcessError, ARK_ILL_INPUT, ARK_INVALID_TABLE, ARK_SUCCESS, ZERO};
use crate::sundials_math::SUNRabs;
use crate::sundials_types::SUN_UNIT_ROUNDOFF;
use crate::sundials_utils::fmt_e;

const ONE: f64 = 1.0;
const TWO: f64 = 2.0;

/* enum MRISTEP_METHOD_TYPE (arkode_mristep.h) */
pub type MRISTEP_METHOD_TYPE = i32;
pub const MRISTEP_EXPLICIT: i32 = 0;
pub const MRISTEP_IMPLICIT: i32 = 1;
pub const MRISTEP_IMEX: i32 = 2;
pub const MRISTEP_MERK: i32 = 3;
pub const MRISTEP_SR: i32 = 4;

/* enum ARKODE_MRITableID (arkode_mristep.h) */
pub type ARKODE_MRITableID = i32;
pub const ARKODE_MRI_NONE: i32 = -1; /* ensure enum is signed int */
pub const ARKODE_MIS_KW3: i32 = 200;
pub const ARKODE_MIN_MRI_NUM: i32 = 200;
pub const ARKODE_MRI_GARK_ERK33a: i32 = 201;
pub const ARKODE_MRI_GARK_ERK45a: i32 = 202;
pub const ARKODE_MRI_GARK_IRK21a: i32 = 203;
pub const ARKODE_MRI_GARK_ESDIRK34a: i32 = 204;
pub const ARKODE_MRI_GARK_ESDIRK46a: i32 = 205;
pub const ARKODE_IMEX_MRI_GARK3a: i32 = 206;
pub const ARKODE_IMEX_MRI_GARK3b: i32 = 207;
pub const ARKODE_IMEX_MRI_GARK4: i32 = 208;
pub const ARKODE_MRI_GARK_FORWARD_EULER: i32 = 209;
pub const ARKODE_MRI_GARK_RALSTON2: i32 = 210;
pub const ARKODE_MRI_GARK_ERK22a: i32 = 211;
pub const ARKODE_MRI_GARK_ERK22b: i32 = 212;
pub const ARKODE_MRI_GARK_RALSTON3: i32 = 213;
pub const ARKODE_MRI_GARK_BACKWARD_EULER: i32 = 214;
pub const ARKODE_MRI_GARK_IMPLICIT_MIDPOINT: i32 = 215;
pub const ARKODE_IMEX_MRI_GARK_EULER: i32 = 216;
pub const ARKODE_IMEX_MRI_GARK_TRAPEZOIDAL: i32 = 217;
pub const ARKODE_IMEX_MRI_GARK_MIDPOINT: i32 = 218;
pub const ARKODE_MERK21: i32 = 219;
pub const ARKODE_MERK32: i32 = 220;
pub const ARKODE_MERK43: i32 = 221;
pub const ARKODE_MERK54: i32 = 222;
pub const ARKODE_IMEX_MRI_SR21: i32 = 223;
pub const ARKODE_IMEX_MRI_SR32: i32 = 224;
pub const ARKODE_IMEX_MRI_SR43: i32 = 225;
pub const ARKODE_MAX_MRI_NUM: i32 = ARKODE_IMEX_MRI_SR43;

/* Default MRI coupling tables for each order and type (arkode_mristep.h) */
pub const MRISTEP_DEFAULT_EXPL_1: i32 = ARKODE_MRI_GARK_FORWARD_EULER;
pub const MRISTEP_DEFAULT_EXPL_2: i32 = ARKODE_MRI_GARK_ERK22b;
pub const MRISTEP_DEFAULT_EXPL_3: i32 = ARKODE_MIS_KW3;
pub const MRISTEP_DEFAULT_EXPL_4: i32 = ARKODE_MRI_GARK_ERK45a;

pub const MRISTEP_DEFAULT_EXPL_2_AD: i32 = ARKODE_MRI_GARK_ERK22b;
pub const MRISTEP_DEFAULT_EXPL_3_AD: i32 = ARKODE_MRI_GARK_ERK33a;
pub const MRISTEP_DEFAULT_EXPL_4_AD: i32 = ARKODE_MRI_GARK_ERK45a;
pub const MRISTEP_DEFAULT_EXPL_5_AD: i32 = ARKODE_MERK54;

pub const MRISTEP_DEFAULT_IMPL_SD_1: i32 = ARKODE_MRI_GARK_BACKWARD_EULER;
pub const MRISTEP_DEFAULT_IMPL_SD_2: i32 = ARKODE_MRI_GARK_IRK21a;
pub const MRISTEP_DEFAULT_IMPL_SD_3: i32 = ARKODE_MRI_GARK_ESDIRK34a;
pub const MRISTEP_DEFAULT_IMPL_SD_4: i32 = ARKODE_MRI_GARK_ESDIRK46a;

pub const MRISTEP_DEFAULT_IMEX_SD_1: i32 = ARKODE_IMEX_MRI_GARK_EULER;
pub const MRISTEP_DEFAULT_IMEX_SD_2: i32 = ARKODE_IMEX_MRI_GARK_TRAPEZOIDAL;
pub const MRISTEP_DEFAULT_IMEX_SD_3: i32 = ARKODE_IMEX_MRI_GARK3b;
pub const MRISTEP_DEFAULT_IMEX_SD_4: i32 = ARKODE_IMEX_MRI_GARK4;

pub const MRISTEP_DEFAULT_IMEX_SD_2_AD: i32 = ARKODE_IMEX_MRI_SR21;
pub const MRISTEP_DEFAULT_IMEX_SD_3_AD: i32 = ARKODE_IMEX_MRI_SR32;
pub const MRISTEP_DEFAULT_IMEX_SD_4_AD: i32 = ARKODE_IMEX_MRI_SR43;

/// struct MRIStepCouplingMem
pub struct MRIStepCouplingMem {
    /// flag to encode the MRI method type (C field `type`)
    pub type_: MRISTEP_METHOD_TYPE,
    /// number of MRI coupling matrices
    pub nmat: i32,
    /// size of coupling matrices ((stages+1) * stages)
    pub stages: i32,
    /// method order of accuracy
    pub q: i32,
    /// embedding order of accuracy
    pub p: i32,
    /// stage abscissae
    pub c: Vec<f64>,
    /// explicit coupling matrices \[nmat\]\[stages+1\]\[stages\] (empty = NULL)
    pub W: Vec<Vec<Vec<f64>>>,
    /// implicit coupling matrices \[nmat\]\[stages+1\]\[stages\] (empty = NULL)
    pub G: Vec<Vec<Vec<f64>>>,
    /// number of stage groups (MERK-specific)
    pub ngroup: i32,
    /// stages to integrate together (MERK-specific; empty = NULL)
    pub group: Vec<Vec<i32>>,
}

/// C `MRIStepCoupling` (pointer type); NULL <-> None at use sites.
pub type MRIStepCoupling = Box<MRIStepCouplingMem>;

/*---------------------------------------------------------------
  Routine to allocate an empty MRIStepCoupling structure
  ---------------------------------------------------------------*/
pub fn MRIStepCoupling_Alloc(
    nmat: i32,
    stages: i32,
    type_: MRISTEP_METHOD_TYPE,
) -> Option<MRIStepCoupling> {
    /* Check for legal input values */
    if nmat < 1 || stages < 1 {
        return None;
    }

    /* Determine general storage format */
    let hasOmegas = type_ == MRISTEP_EXPLICIT
        || type_ == MRISTEP_IMEX
        || type_ == MRISTEP_MERK
        || type_ == MRISTEP_SR;
    let hasGammas = type_ == MRISTEP_IMPLICIT || type_ == MRISTEP_IMEX || type_ == MRISTEP_SR;

    /* Allocate abscissae and coupling coefficients (calloc'd in C, so
    only non-zero coefficients need to be set) */
    let mut MRIC = Box::new(MRIStepCouplingMem {
        type_,
        nmat,
        stages,
        q: 0,
        p: 0,
        c: vec![0.0; stages as usize],
        W: Vec::new(),
        G: Vec::new(),
        ngroup: 0,
        group: Vec::new(),
    });

    if hasOmegas {
        MRIC.W = vec![vec![vec![0.0; stages as usize]; (stages + 1) as usize]; nmat as usize];
    }

    if hasGammas {
        MRIC.G = vec![vec![vec![0.0; stages as usize]; (stages + 1) as usize]; nmat as usize];
    }

    /* for MERK methods, allocate maximum possible number/sizes of stage groups */
    if type_ == MRISTEP_MERK {
        MRIC.ngroup = stages;
        MRIC.group = vec![vec![-1; stages as usize]; stages as usize];
    }

    Some(MRIC)
}

/*---------------------------------------------------------------
  Routine to allocate and fill an explicit, implicit, or ImEx
  MRIGARK MRIStepCoupling structure (W/G passed flattened in C's
  1D layouts).
  ---------------------------------------------------------------*/
pub fn MRIStepCoupling_Create(
    nmat: i32,
    stages: i32,
    q: i32,
    p: i32,
    W: Option<&[f64]>,
    G: Option<&[f64]>,
    c: &[f64],
) -> Option<MRIStepCoupling> {
    /* Check for legal inputs (C also rejects c == NULL) */
    if nmat < 1 || stages < 1 {
        return None;
    }

    /* Check for method coefficients and set method type */
    let type_ = match (W.is_some(), G.is_some()) {
        (true, true) => MRISTEP_IMEX,
        (true, false) => MRISTEP_EXPLICIT,
        (false, true) => MRISTEP_IMPLICIT,
        (false, false) => return None,
    };

    /* Allocate MRIStepCoupling structure */
    let mut MRIC = MRIStepCoupling_Alloc(nmat, stages, type_)?;

    /* Method and embedding order */
    MRIC.q = q;
    MRIC.p = p;

    /* Abscissae */
    MRIC.c[..stages as usize].copy_from_slice(&c[..stages as usize]);

    /* Coupling coefficients stored as 1D arrays, based on whether they
       include embedding coefficients */
    let s = stages as usize;
    if p == 0 {
        /* non-embedded method:  coupling coefficient 1D arrays have
           length nmat * stages * stages, with each stages * stages
           matrix stored in C (row-major) order */
        for k in 0..nmat as usize {
            for i in 0..s {
                for j in 0..s {
                    if let Some(w) = W {
                        MRIC.W[k][i][j] = w[s * (s * k + i) + j];
                    }
                    if let Some(g) = G {
                        MRIC.G[k][i][j] = g[s * (s * k + i) + j];
                    }
                }
            }
        }
    } else {
        /* embedded method:  coupling coefficient 1D arrays have
           length nmat * (stages+1) * stages, with each (stages+1) * stages
           matrix stored in C (row-major) order */
        for k in 0..nmat as usize {
            for i in 0..=s {
                for j in 0..s {
                    if let Some(w) = W {
                        MRIC.W[k][i][j] = w[(s + 1) * (s * k + i) + j];
                    }
                    if let Some(g) = G {
                        MRIC.G[k][i][j] = g[(s + 1) * (s * k + i) + j];
                    }
                }
            }
        }
    }
    Some(MRIC)
}

/*---------------------------------------------------------------
  Construct the MRIGARK coupling matrix for an MIS method based
  on a given "slow" Butcher table.
  ---------------------------------------------------------------*/
pub fn MRIStepCoupling_MIStoMRI(B: &ARKodeButcherTable, q: i32, p: i32) -> Option<MRIStepCoupling> {
    let tol: f64 = 100.0 * SUN_UNIT_ROUNDOFF;

    /* If p>0, check that input table includes embedding coefficients */
    if p > 0 && B.d.is_none() {
        return None;
    }

    /* -----------------------------------
     * Check that the input table is valid
     * ----------------------------------- */

    let bs = B.stages as usize;

    /* First stage is just old solution */
    let mut Asum = SUNRabs(B.c[0]);
    for j in 0..bs {
        Asum += SUNRabs(B.A[0][j]);
    }
    if Asum > tol {
        return None;
    }

    /* Last stage exceeds 1 */
    if B.c[bs - 1] > ONE + tol {
        return None;
    }

    /* All stages are sorted */
    for j in 1..bs {
        if (B.c[j] - B.c[j - 1]) < -tol {
            return None;
        }
    }

    /* Each stage at most diagonally implicit */
    Asum = ZERO;
    for i in 0..bs {
        for j in (i + 1)..bs {
            Asum += SUNRabs(B.A[i][j]);
        }
    }
    if Asum > tol {
        return None;
    }

    /* -----------------------------------------
     * determine whether the table needs padding
     * ----------------------------------------- */

    let mut padding = false;

    /* Pad if last stage does not equal 1 */
    if SUNRabs(B.c[bs - 1] - ONE) > tol {
        padding = true;
    }

    /* Pad if last row of A does not equal b */
    for j in 0..bs {
        if SUNRabs(B.A[bs - 1][j] - B.b[j]) > tol {
            padding = true;
        }
    }

    /* If final stage is implicit and the method contains an embedding,
       we require padding since d != b */
    if p > 0 && SUNRabs(B.A[bs - 1][bs - 1]) > tol {
        padding = true;
    }
    let stages = if padding { B.stages + 1 } else { B.stages };

    /* -------------------------
     * determine the method type
     * ------------------------- */

    /* Check if the table is strictly lower triangular (explicit) */
    let mut type_ = MRISTEP_EXPLICIT;

    for i in 0..bs {
        for j in i..bs {
            if SUNRabs(B.A[i][j]) > tol {
                type_ = MRISTEP_IMPLICIT;
            }
        }
    }

    /* ----------------------------
     * construct coupling structure
     * ---------------------------- */

    let mut MRIC = MRIStepCoupling_Alloc(1, stages, type_)?;

    /* Copy method/embedding orders */
    MRIC.q = q;
    MRIC.p = p;

    /* Copy abscissae, padding if needed */
    for i in 0..bs {
        MRIC.c[i] = B.c[i];
    }

    if padding {
        MRIC.c[stages as usize - 1] = ONE;
    }

    /* Construct the coupling table */
    let C = if type_ == MRISTEP_EXPLICIT {
        &mut MRIC.W
    } else {
        &mut MRIC.G
    };

    /* First row is identically zero */
    for i in 0..stages as usize {
        for j in 0..stages as usize {
            C[0][i][j] = ZERO;
        }
    }

    /* Remaining rows = A(2:end,:) - A(1:end-1,:) */
    for i in 1..bs {
        for j in 0..bs {
            C[0][i][j] = B.A[i][j] - B.A[i - 1][j];
        }
    }

    /* Padded row = b(:) - A(end,:) */
    if padding {
        for j in 0..bs {
            C[0][stages as usize - 1][j] = B.b[j] - B.A[bs - 1][j];
        }
    }

    /* Embedded row = d(:) - A(end,:) */
    if p > 0 {
        let d = B.d.as_ref().unwrap();
        for j in 0..bs {
            C[0][stages as usize][j] = d[j] - B.A[bs - 1][j];
        }
    }

    Some(MRIC)
}

/*---------------------------------------------------------------
  Routine to copy a MRIStepCoupling structure
  ---------------------------------------------------------------*/
pub fn MRIStepCoupling_Copy(MRIC: &MRIStepCouplingMem) -> Option<MRIStepCoupling> {
    /* Check for stage times */
    if MRIC.c.is_empty() {
        return None;
    }

    /* Allocate coupling structure */
    let mut MRICcopy = MRIStepCoupling_Alloc(MRIC.nmat, MRIC.stages, MRIC.type_)?;

    /* Copy method and embedding orders */
    MRICcopy.q = MRIC.q;
    MRICcopy.p = MRIC.p;

    /* Copy abscissae */
    MRICcopy.c.copy_from_slice(&MRIC.c);

    /* Copy explicit coupling matrices W */
    if !MRIC.W.is_empty() {
        MRICcopy.W.clone_from(&MRIC.W);
    }

    /* Copy implicit coupling matrices G */
    if !MRIC.G.is_empty() {
        MRICcopy.G.clone_from(&MRIC.G);
    }

    /* Copy MERK stage groups */
    if !MRIC.group.is_empty() {
        MRICcopy.ngroup = MRIC.ngroup;
        MRICcopy.group.clone_from(&MRIC.group);
    }

    Some(MRICcopy)
}

/*---------------------------------------------------------------
  Routine to query the MRIStepCoupling structure workspace size
  ---------------------------------------------------------------*/
pub fn MRIStepCoupling_Space(MRIC: Option<&MRIStepCouplingMem>, liw: &mut i64, lrw: &mut i64) {
    /* initialize outputs and return if MRIC is not allocated */
    *liw = 0;
    *lrw = 0;
    let MRIC = match MRIC {
        Some(m) => m,
        None => return,
    };

    /* fill outputs based on MRIC */
    *liw = 5;
    if !MRIC.c.is_empty() {
        *lrw += MRIC.stages as i64;
    }
    if !MRIC.W.is_empty() {
        *lrw += (MRIC.nmat * (MRIC.stages + 1) * MRIC.stages) as i64;
    }
    if !MRIC.G.is_empty() {
        *lrw += (MRIC.nmat * (MRIC.stages + 1) * MRIC.stages) as i64;
    }
    if !MRIC.group.is_empty() {
        *liw += (MRIC.stages * MRIC.stages) as i64;
    }
}

/*---------------------------------------------------------------
  Routine to free a MRIStepCoupling structure
  ---------------------------------------------------------------*/
pub fn MRIStepCoupling_Free(MRIC: &mut Option<MRIStepCoupling>) {
    *MRIC = None;
}

/* SUN_FORMAT_E == "% .15e" for double precision */
fn sun_format_e(x: f64) -> String {
    let e = fmt_e(x, 0, 15);
    if e.starts_with('-') {
        e
    } else {
        format!(" {}", e)
    }
}

/*---------------------------------------------------------------
  Routine to print a MRIStepCoupling structure
  ---------------------------------------------------------------*/
pub fn MRIStepCoupling_Write(MRIC: &MRIStepCouplingMem, outfile: &mut dyn std::io::Write) {
    /* check for valid coupling structure */
    if MRIC.W.is_empty() && MRIC.G.is_empty() {
        return;
    }
    if MRIC.c.is_empty() {
        return;
    }

    match MRIC.type_ {
        MRISTEP_EXPLICIT => {
            let _ = writeln!(outfile, "  type = explicit MRI");
        }
        MRISTEP_IMPLICIT => {
            let _ = writeln!(outfile, "  type = implicit MRI");
        }
        MRISTEP_IMEX => {
            let _ = writeln!(outfile, "  type = ImEx MRI");
        }
        MRISTEP_MERK => {
            let _ = writeln!(outfile, "  type = MERK");
        }
        MRISTEP_SR => {
            let _ = writeln!(outfile, "  type = MRISR");
        }
        _ => {
            let _ = writeln!(outfile, "  type = unknown");
        }
    }
    let _ = writeln!(outfile, "  nmat = {}", MRIC.nmat);
    let _ = writeln!(outfile, "  stages = {}", MRIC.stages);
    let _ = writeln!(outfile, "  method order (q) = {}", MRIC.q);
    let _ = writeln!(outfile, "  embedding order (p) = {}", MRIC.p);
    let _ = write!(outfile, "  c = ");
    for i in 0..MRIC.stages as usize {
        let _ = write!(outfile, "{}  ", sun_format_e(MRIC.c[i]));
    }
    let _ = writeln!(outfile);

    if !MRIC.W.is_empty() {
        for k in 0..MRIC.nmat as usize {
            let _ = writeln!(outfile, "  W[{}] = ", k);
            for i in 0..=MRIC.stages as usize {
                let _ = write!(outfile, "      ");
                for j in 0..MRIC.stages as usize {
                    let _ = write!(outfile, "{}  ", sun_format_e(MRIC.W[k][i][j]));
                }
                let _ = writeln!(outfile);
            }
            let _ = writeln!(outfile);
        }
    }

    if !MRIC.G.is_empty() {
        for k in 0..MRIC.nmat as usize {
            let _ = writeln!(outfile, "  G[{}] = ", k);
            for i in 0..=MRIC.stages as usize {
                let _ = write!(outfile, "      ");
                for j in 0..MRIC.stages as usize {
                    let _ = write!(outfile, "{}  ", sun_format_e(MRIC.G[k][i][j]));
                }
                let _ = writeln!(outfile);
            }
            let _ = writeln!(outfile);
        }
    }

    if !MRIC.group.is_empty() {
        let _ = writeln!(outfile, "  ngroup = {}", MRIC.ngroup);
        for i in 0..MRIC.ngroup as usize {
            let _ = write!(outfile, "  group[{}] = ", i);
            for j in 0..MRIC.stages as usize {
                if MRIC.group[i][j] >= 0 {
                    let _ = write!(outfile, "{} ", MRIC.group[i][j]);
                }
            }
            let _ = writeln!(outfile);
        }
    }
}

/* ===========================================================================
 * Private Functions
 * ===========================================================================*/

/* ---------------------------------------------------------------------------
 * Stage type identifier: returns one of the constants
 *
 * MRISTAGE_ERK_FAST    -- standard MIS-like stage
 * MRISTAGE_ERK_NOFAST  -- standard ERK stage
 * MRISTAGE_DIRK_NOFAST -- standard DIRK stage
 * MRISTAGE_DIRK_FAST   -- coupled DIRK + MIS-like stage
 * MRISTAGE_STIFF_ACC   -- "extra" stiffly-accurate stage
 *
 * for each nontrivial stage, or embedding stage, in an MRI-like method.
 * Otherwise (i.e., stage is not in [1,MRIC->stages]), returns
 * ARK_INVALID_TABLE (<0).
 * ---------------------------------------------------------------------------*/
pub fn mriStepCoupling_GetStageType(MRIC: &MRIStepCouplingMem, is: i32) -> i32 {
    use crate::arkode_mristep_impl::{
        MRISTAGE_DIRK_FAST, MRISTAGE_DIRK_NOFAST, MRISTAGE_ERK_FAST, MRISTAGE_ERK_NOFAST,
        MRISTAGE_FIRST, MRISTAGE_STIFF_ACC,
    };
    let tol: f64 = 100.0 * SUN_UNIT_ROUNDOFF;

    if is < 0 || is > MRIC.stages {
        return ARK_INVALID_TABLE;
    }

    if is == 0 {
        return MRISTAGE_FIRST;
    }

    /* report MRISTAGE_ERK_FAST for MERK and MRI-SR methods */
    if MRIC.type_ == MRISTEP_SR || MRIC.type_ == MRISTEP_MERK {
        return MRISTAGE_ERK_FAST;
    }

    let isu = is as usize;
    let mut Gdiag = false;
    let mut Grow = false;
    let mut Wrow = false;
    let cdiff;

    /* separately handle an embedding "stage" from normal stages */
    if is < MRIC.stages {
        /* normal */
        if !MRIC.G.is_empty() {
            for i in 0..MRIC.nmat as usize {
                Gdiag = Gdiag || (SUNRabs(MRIC.G[i][isu][isu]) > tol);
                for j in 0..MRIC.stages as usize {
                    Grow = Grow || (SUNRabs(MRIC.G[i][isu][j]) > tol);
                }
            }
        }
        if !MRIC.W.is_empty() {
            for i in 0..MRIC.nmat as usize {
                for j in 0..MRIC.stages as usize {
                    Wrow = Wrow || (SUNRabs(MRIC.W[i][isu][j]) > tol);
                }
            }
        }

        /* abscissae difference */
        cdiff = SUNRabs(MRIC.c[isu] - MRIC.c[isu - 1]) > tol;
    } else {
        /* embedding */
        if !MRIC.G.is_empty() {
            for i in 0..MRIC.nmat as usize {
                Gdiag = Gdiag || (SUNRabs(MRIC.G[i][isu][isu - 1]) > tol);
                for j in 0..MRIC.stages as usize {
                    Grow = Grow || (SUNRabs(MRIC.G[i][isu][j]) > tol);
                }
            }
        }
        if !MRIC.W.is_empty() {
            for i in 0..MRIC.nmat as usize {
                for j in 0..MRIC.stages as usize {
                    Wrow = Wrow || (SUNRabs(MRIC.W[i][isu][j]) > tol);
                }
            }
        }
        cdiff = SUNRabs(MRIC.c[isu - 1] - MRIC.c[isu - 2]) > tol;
    }

    /* make determination */
    if !(Gdiag || Grow || Wrow || cdiff) && is > 0 {
        /* stiffly-accurate stage */
        return MRISTAGE_STIFF_ACC;
    }
    if Gdiag {
        /* DIRK */
        if cdiff {
            /* Fast */
            MRISTAGE_DIRK_FAST
        } else {
            MRISTAGE_DIRK_NOFAST
        }
    } else {
        /* ERK */
        if cdiff {
            /* Fast */
            MRISTAGE_ERK_FAST
        } else {
            MRISTAGE_ERK_NOFAST
        }
    }
}

/* ---------------------------------------------------------------------------
 * Computes the stage RHS vector storage maps. With repeated abscissae the
 * first stage of the pair generally corresponds to a column of zeros and so
 * does not need to be computed and stored. The stage_map indicates if the RHS
 * needs to be computed and where to store it i.e., stage_map[i] > -1.
 *
 * Note: for MERK and MRI-SR methods, this should be an "identity" map, and all
 * stage vectors should be allocated.
 * ---------------------------------------------------------------------------*/
pub fn mriStepCoupling_GetStageMap(
    MRIC: &MRIStepCouplingMem,
    stage_map: &mut [i32],
    nstages_active: &mut i32,
) -> i32 {
    let tol: f64 = 100.0 * SUN_UNIT_ROUNDOFF;

    /* Check for valid inputs */
    if MRIC.W.is_empty() && MRIC.G.is_empty() {
        return ARK_ILL_INPUT;
    }

    /* MERK and MRI-SR have "identity" storage map */
    if MRIC.type_ == MRISTEP_MERK || MRIC.type_ == MRISTEP_SR {
        /* Number of stage RHS vectors active */
        *nstages_active = MRIC.stages;

        /* Create an identity map (all columns are non-zero) */
        for (j, sm) in stage_map.iter_mut().enumerate().take(MRIC.stages as usize) {
            *sm = j as i32;
        }
        return ARK_SUCCESS;
    }

    /* Compute storage map for MRI-GARK methods */

    /* Number of stage RHS vectors active */
    *nstages_active = 0;

    /* Initial storage index */
    let mut idx: i32 = 0;

    /* Check if a stage corresponds to a column of zeros for all coupling
     * matrices by computing the column sums */
    for j in 0..MRIC.stages as usize {
        let mut Wsum = ZERO;
        let mut Gsum = ZERO;

        if !MRIC.W.is_empty() {
            for k in 0..MRIC.nmat as usize {
                for i in 0..=MRIC.stages as usize {
                    Wsum += SUNRabs(MRIC.W[k][i][j]);
                }
            }
        }

        if !MRIC.G.is_empty() {
            for k in 0..MRIC.nmat as usize {
                for i in 0..=MRIC.stages as usize {
                    Gsum += SUNRabs(MRIC.G[k][i][j]);
                }
            }
        }

        if Wsum > tol || Gsum > tol {
            stage_map[j] = idx;
            idx += 1;
        } else {
            stage_map[j] = -1;
        }
    }

    /* Check and set number of stage RHS vectors active */
    if idx < 1 {
        return ARK_ILL_INPUT;
    }

    *nstages_active = idx;

    ARK_SUCCESS
}

/* ===========================================================================
 * Coupling tables (arkode_mri_tables.def)
 * ===========================================================================*/

fn mri_table_from_erk(id: crate::arkode_butcher_erk::ARKODE_ERKTableID) -> Option<MRIStepCoupling> {
    let B = ARKodeButcherTable_LoadERK(id)?;
    MRIStepCoupling_MIStoMRI(&B, B.q, B.p)
}

/* ARKODE_MRI_GARK_IRK21a: A. Sandu, SINUM 57:2300-2327, 2019 */
fn mri_table_irk21a() -> Option<MRIStepCoupling> {
    let mut B = ARKodeButcherTable_Alloc(3, true)?;

    B.q = 2;
    B.p = 1;

    B.c[1] = ONE;
    B.c[2] = ONE;

    B.A[1][0] = ONE;
    B.A[2][0] = 0.5;
    B.A[2][2] = 0.5;

    B.b[0] = 0.5;
    B.b[2] = 0.5;

    B.d.as_mut().unwrap()[2] = 1.0;

    MRIStepCoupling_MIStoMRI(&B, B.q, B.p)
}

/* ARKODE_MRI_GARK_ERK33a: A. Sandu, SINUM 57:2300-2327, 2019 */
fn mri_table_erk33a() -> Option<MRIStepCoupling> {
    let mut C = MRIStepCoupling_Alloc(2, 4, MRISTEP_EXPLICIT)?;

    C.q = 3;
    C.p = 2;

    C.c[1] = ONE / 3.0;
    C.c[2] = TWO / 3.0;
    C.c[3] = ONE;

    C.W[0][1][0] = ONE / 3.0;
    C.W[0][2][0] = -ONE / 3.0;
    C.W[0][2][1] = TWO / 3.0;
    C.W[0][3][1] = -TWO / 3.0;
    C.W[0][3][2] = ONE;
    C.W[0][4][0] = ONE / 12.0;
    C.W[0][4][1] = -ONE / 3.0;
    C.W[0][4][2] = 7.0 / 12.0;

    C.W[1][3][0] = ONE / TWO;
    C.W[1][3][2] = -ONE / TWO;
    Some(C)
}

/* ARKODE_MRI_GARK_RALSTON3: Roberts et al., SISC 44:A1405 - A1427, 2022 */
fn mri_table_ralston3() -> Option<MRIStepCoupling> {
    let mut C = MRIStepCoupling_Alloc(2, 4, MRISTEP_EXPLICIT)?;

    C.q = 3;
    C.p = 2;

    C.c[1] = ONE / TWO;
    C.c[2] = 3.0 / 4.0;
    C.c[3] = ONE;

    C.W[0][1][0] = ONE / TWO;
    C.W[0][2][0] = -11.0 / 4.0;
    C.W[0][2][1] = 3.0;
    C.W[0][3][0] = 47.0 / 36.0;
    C.W[0][3][1] = -ONE / 6.0;
    C.W[0][3][2] = -8.0 / 9.0;
    C.W[0][4][0] = ONE / 40.0;
    C.W[0][4][1] = 7.0 / 40.0;
    C.W[0][4][2] = ONE / 20.0;

    C.W[1][2][0] = 9.0 / TWO;
    C.W[1][2][1] = -9.0 / TWO;
    C.W[1][3][0] = -13.0 / 6.0;
    C.W[1][3][1] = -ONE / TWO;
    C.W[1][3][2] = 8.0 / 3.0;
    Some(C)
}

/* ARKODE_MRI_GARK_ERK45a: A. Sandu, SINUM 57:2300-2327, 2019
   (embedding coefficients CORRECTED in A. Sandu, arxiv:1808.02759, 2018) */
fn mri_table_erk45a() -> Option<MRIStepCoupling> {
    let mut C = MRIStepCoupling_Alloc(2, 6, MRISTEP_EXPLICIT)?;

    C.q = 4;
    C.p = 3;

    C.c[1] = 0.2;
    C.c[2] = 0.4;
    C.c[3] = 0.6;
    C.c[4] = 0.8;
    C.c[5] = ONE;

    C.W[0][1][0] = 0.2;
    C.W[0][2][0] = -53.0 / 16.0;
    C.W[0][2][1] = 281.0 / 80.0;
    C.W[0][3][0] = -36562993.0 / 71394880.0;
    C.W[0][3][1] = 34903117.0 / 17848720.0;
    C.W[0][3][2] = -88770499.0 / 71394880.0;
    C.W[0][4][0] = -7631593.0 / 71394880.0;
    C.W[0][4][1] = -166232021.0 / 35697440.0;
    C.W[0][4][2] = 6068517.0 / 1519040.0;
    C.W[0][4][3] = 8644289.0 / 8924360.0;
    C.W[0][5][0] = 277061.0 / 303808.0;
    C.W[0][5][1] = -209323.0 / 1139280.0;
    C.W[0][5][2] = -1360217.0 / 1139280.0;
    C.W[0][5][3] = -148789.0 / 56964.0;
    C.W[0][5][4] = 147889.0 / 45120.0;
    C.W[0][6][0] = -88227.0 / 47470.0;
    C.W[0][6][1] = 756870829.0 / 340217490.0;
    C.W[0][6][2] = -713704111.0 / 1360869960.0;
    C.W[0][6][3] = -31967827.0 / 340217490.0;
    C.W[0][6][4] = 129673.0 / 286680.0;

    C.W[1][2][0] = 503.0 / 80.0;
    C.W[1][2][1] = -503.0 / 80.0;
    C.W[1][3][0] = -1365537.0 / 35697440.0;
    C.W[1][3][1] = 4963773.0 / 7139488.0;
    C.W[1][3][2] = -1465833.0 / 2231090.0;
    C.W[1][4][0] = 66974357.0 / 35697440.0;
    C.W[1][4][1] = 21445367.0 / 7139488.0;
    C.W[1][4][2] = -3.0;
    C.W[1][4][3] = -8388609.0 / 4462180.0;
    C.W[1][5][0] = -18227.0 / 7520.0;
    C.W[1][5][1] = TWO;
    C.W[1][5][2] = ONE;
    C.W[1][5][3] = 5.0;
    C.W[1][5][4] = -41933.0 / 7520.0;
    C.W[1][6][0] = 6213.0 / 1880.0;
    C.W[1][6][1] = -6213.0 / 1880.0;
    Some(C)
}

/* ARKODE_MRI_GARK_BACKWARD_EULER */
fn mri_table_backward_euler() -> Option<MRIStepCoupling> {
    let mut C = MRIStepCoupling_Alloc(1, 3, MRISTEP_IMPLICIT)?;

    C.q = 1;
    C.p = 0;

    C.c[1] = ONE;
    C.c[2] = ONE;

    C.G[0][1][0] = ONE;
    C.G[0][2][0] = -ONE;
    C.G[0][2][2] = ONE;
    Some(C)
}

/* ARKODE_MRI_GARK_IMPLICIT_MIDPOINT */
fn mri_table_implicit_midpoint() -> Option<MRIStepCoupling> {
    let mut C = MRIStepCoupling_Alloc(1, 4, MRISTEP_IMPLICIT)?;

    C.q = 2;
    C.p = 0;

    C.c[1] = ONE / TWO;
    C.c[2] = ONE / TWO;
    C.c[3] = ONE;

    C.G[0][1][0] = ONE / TWO;
    C.G[0][2][0] = -ONE / TWO;
    C.G[0][2][2] = ONE / TWO;
    C.G[0][3][2] = ONE / TWO;
    Some(C)
}

/* ARKODE_MRI_GARK_ESDIRK34a: A. Sandu, SINUM 57:2300-2327, 2019 */
fn mri_table_esdirk34a() -> Option<MRIStepCoupling> {
    let mut C = MRIStepCoupling_Alloc(1, 8, MRISTEP_IMPLICIT)?;
    let beta: f64 = 0.4358665215084589994160194511935568425;

    C.q = 3;
    C.p = 2;

    C.c[1] = ONE / 3.0;
    C.c[2] = ONE / 3.0;
    C.c[3] = TWO / 3.0;
    C.c[4] = TWO / 3.0;
    C.c[5] = ONE;
    C.c[6] = ONE;
    C.c[7] = ONE;

    C.G[0][1][0] = ONE / 3.0;
    C.G[0][2][0] = -beta;
    C.G[0][2][2] = beta;
    C.G[0][3][0] = -0.3045790611944504970424837655380884888;
    C.G[0][3][2] = 0.6379123945277838303758170988714218222;
    C.G[0][4][0] = 0.2116913105640266601676536489364004869;
    C.G[0][4][2] = -0.6475578320724856595836731001299573294;
    C.G[0][4][4] = beta;
    C.G[0][5][0] = 0.4454209388055495029575162344619115112;
    C.G[0][5][2] = 0.8813784805616198280398949036456491923;
    C.G[0][5][4] = -0.9934660860338359976640778047742273701;
    C.G[0][6][0] = -beta;
    C.G[0][6][6] = beta;
    C.G[0][8][0] = 0.2453831999117524372455680781104585876241;
    C.G[0][8][2] = 0.4204215033044044563073464989473988121422;
    C.G[0][8][4] = -1.576992606344066224351397232226173387157;
    C.G[0][8][6] = 0.9111879031279093307984826551683159873903;
    Some(C)
}

/* ARKODE_MRI_GARK_ESDIRK46a: A. Sandu, SINUM 57:2300-2327, 2019 */
fn mri_table_esdirk46a() -> Option<MRIStepCoupling> {
    let mut C = MRIStepCoupling_Alloc(2, 12, MRISTEP_IMPLICIT)?;

    C.q = 4;
    C.p = 3;

    C.c[1] = ONE / 5.0;
    C.c[2] = ONE / 5.0;
    C.c[3] = TWO / 5.0;
    C.c[4] = TWO / 5.0;
    C.c[5] = 3.0 / 5.0;
    C.c[6] = 3.0 / 5.0;
    C.c[7] = 4.0 / 5.0;
    C.c[8] = 4.0 / 5.0;
    C.c[9] = ONE;
    C.c[10] = ONE;
    C.c[11] = ONE;

    C.G[0][1][0] = ONE / 5.0;
    C.G[0][2][0] = -ONE / 4.0;
    C.G[0][2][2] = ONE / 4.0;
    C.G[0][3][0] = 1771023115159.0 / 1929363690800.0;
    C.G[0][3][2] = -1385150376999.0 / 1929363690800.0;
    C.G[0][4][0] = 914009.0 / 345800.0;
    C.G[0][4][2] = -1000459.0 / 345800.0;
    C.G[0][4][4] = ONE / 4.0;
    C.G[0][5][0] = 18386293581909.0 / 36657910125200.0;
    C.G[0][5][2] = 5506531089.0 / 80566835440.0;
    C.G[0][5][4] = -178423463189.0 / 482340922700.0;
    C.G[0][6][0] = 36036097.0 / 8299200.0;
    C.G[0][6][2] = 4621.0 / 118560.0;
    C.G[0][6][4] = -38434367.0 / 8299200.0;
    C.G[0][6][6] = ONE / 4.0;
    C.G[0][7][0] = -247809665162987.0 / 146631640500800.0;
    C.G[0][7][2] = 10604946373579.0 / 14663164050080.0;
    C.G[0][7][4] = 10838126175385.0 / 5865265620032.0;
    C.G[0][7][6] = -24966656214317.0 / 36657910125200.0;
    C.G[0][8][0] = 38519701.0 / 11618880.0;
    C.G[0][8][2] = 10517363.0 / 9682400.0;
    C.G[0][8][4] = -23284701.0 / 19364800.0;
    C.G[0][8][6] = -10018609.0 / 2904720.0;
    C.G[0][8][8] = ONE / 4.0;
    C.G[0][9][0] = -52907807977903.0 / 33838070884800.0;
    C.G[0][9][2] = 74846944529257.0 / 73315820250400.0;
    C.G[0][9][4] = 365022522318171.0 / 146631640500800.0;
    C.G[0][9][6] = -20513210406809.0 / 109973730375600.0;
    C.G[0][9][8] = -2918009798.0 / 1870301537.0;
    C.G[0][10][0] = 19.0 / 100.0;
    C.G[0][10][2] = -73.0 / 300.0;
    C.G[0][10][4] = 127.0 / 300.0;
    C.G[0][10][6] = 127.0 / 300.0;
    C.G[0][10][8] = -313.0 / 300.0;
    C.G[0][10][10] = ONE / 4.0;
    C.G[0][12][0] = -ONE / 4.0;
    C.G[0][12][2] = 5595.0 / 8804.0;
    C.G[0][12][4] = -2445.0 / 8804.0;
    C.G[0][12][6] = -4225.0 / 8804.0;
    C.G[0][12][8] = 2205.0 / 4402.0;
    C.G[0][12][10] = -567.0 / 4402.0;

    C.G[1][3][0] = -1674554930619.0 / 964681845400.0;
    C.G[1][3][2] = 1674554930619.0 / 964681845400.0;
    C.G[1][4][0] = -1007739.0 / 172900.0;
    C.G[1][4][2] = 1007739.0 / 172900.0;
    C.G[1][5][0] = -8450070574289.0 / 18328955062600.0;
    C.G[1][5][2] = -39429409169.0 / 40283417720.0;
    C.G[1][5][4] = 173621393067.0 / 120585230675.0;
    C.G[1][6][0] = -122894383.0 / 16598400.0;
    C.G[1][6][2] = 14501.0 / 237120.0;
    C.G[1][6][4] = 121879313.0 / 16598400.0;
    C.G[1][7][0] = 32410002731287.0 / 15434909526400.0;
    C.G[1][7][2] = -46499276605921.0 / 29326328100160.0;
    C.G[1][7][4] = -34914135774643.0 / 11730531240064.0;
    C.G[1][7][6] = 45128506783177.0 / 18328955062600.0;
    C.G[1][8][0] = -128357303.0 / 23237760.0;
    C.G[1][8][2] = -35433927.0 / 19364800.0;
    C.G[1][8][4] = 71038479.0 / 38729600.0;
    C.G[1][8][6] = 8015933.0 / 1452360.0;
    C.G[1][9][0] = 136721604296777.0 / 67676141769600.0;
    C.G[1][9][2] = -349632444539303.0 / 146631640500800.0;
    C.G[1][9][4] = -1292744859249609.0 / 293263281001600.0;
    C.G[1][9][6] = 8356250416309.0 / 54986865187800.0;
    C.G[1][9][8] = 17282943803.0 / 3740603074.0;
    C.G[1][10][0] = 3.0 / 25.0;
    C.G[1][10][2] = -29.0 / 300.0;
    C.G[1][10][4] = 71.0 / 300.0;
    C.G[1][10][6] = 71.0 / 300.0;
    C.G[1][10][8] = -149.0 / 300.0;
    Some(C)
}

/* ARKODE_IMEX_MRI_GARK_EULER */
fn mri_table_imex_euler() -> Option<MRIStepCoupling> {
    let mut C = MRIStepCoupling_Alloc(1, 3, MRISTEP_IMEX)?;

    C.q = 1;
    C.p = 0;

    C.c[1] = ONE;
    C.c[2] = ONE;

    C.W[0][1][0] = ONE;

    C.G[0][1][0] = ONE;
    C.G[0][2][0] = -ONE;
    C.G[0][2][2] = ONE;
    Some(C)
}

/* ARKODE_IMEX_MRI_GARK_TRAPEZOIDAL */
fn mri_table_imex_trapezoidal() -> Option<MRIStepCoupling> {
    let mut C = MRIStepCoupling_Alloc(1, 4, MRISTEP_IMEX)?;

    C.q = 2;
    C.p = 0;

    C.c[1] = ONE;
    C.c[2] = ONE;
    C.c[3] = ONE;

    C.W[0][1][0] = ONE;
    C.W[0][3][0] = -ONE / TWO;
    C.W[0][3][2] = ONE / TWO;

    C.G[0][1][0] = ONE;
    C.G[0][2][0] = -ONE / TWO;
    C.G[0][2][2] = ONE / TWO;
    Some(C)
}

/* ARKODE_IMEX_MRI_GARK_MIDPOINT */
fn mri_table_imex_midpoint() -> Option<MRIStepCoupling> {
    let mut C = MRIStepCoupling_Alloc(1, 4, MRISTEP_IMEX)?;

    C.q = 2;
    C.p = 0;

    C.c[1] = ONE / TWO;
    C.c[2] = ONE / TWO;
    C.c[3] = ONE;

    C.W[0][1][0] = ONE / TWO;
    C.W[0][3][0] = -ONE / TWO;
    C.W[0][3][2] = ONE;

    C.G[0][1][0] = ONE / TWO;
    C.G[0][2][0] = -ONE / TWO;
    C.G[0][2][2] = ONE / TWO;
    C.G[0][3][2] = ONE / TWO;
    Some(C)
}

/* ARKODE_IMEX_MRI_GARK3a: R. Chinomona & D. Reynolds SINUM 43(5):A3082-A3113, 2021 */
fn mri_table_imex_gark3a() -> Option<MRIStepCoupling> {
    let mut C = MRIStepCoupling_Alloc(1, 8, MRISTEP_IMEX)?;
    let beta: f64 = 0.4358665215084589994160194511935568425;

    C.q = 3;
    C.p = 0;

    C.c[1] = beta;
    C.c[2] = beta;
    C.c[3] = 0.7179332607542294997080097255967784213;
    C.c[4] = 0.7179332607542294997080097255967784213;
    C.c[5] = ONE;
    C.c[6] = ONE;
    C.c[7] = ONE;

    C.W[0][1][0] = beta;
    C.W[0][3][0] = -0.5688715801234400928465032925317932021;
    C.W[0][3][2] = 0.8509383193692105931384935669350147809;
    C.W[0][4][0] = 0.454283944643608855878770886900124654;
    C.W[0][4][2] = -0.454283944643608855878770886900124654;
    C.W[0][5][0] = -0.4271371821005074011706645050390732474;
    C.W[0][5][2] = 0.1562747733103380821014660497037023496;
    C.W[0][5][4] = 0.5529291480359398193611887297385924765;
    C.W[0][7][0] = 0.105858296071879638722377459477184953;
    C.W[0][7][2] = 0.655567501140070250975288954324730635;
    C.W[0][7][4] = -1.197292318720408889113685864995472431;
    C.W[0][7][6] = beta;

    C.G[0][1][0] = beta;
    C.G[0][2][0] = -beta;
    C.G[0][2][2] = beta;
    C.G[0][3][0] = -0.4103336962288525014599513720161078937;
    C.G[0][3][2] = 0.6924004354746230017519416464193294724;
    C.G[0][4][0] = 0.4103336962288525014599513720161078937;
    C.G[0][4][2] = -0.8462002177373115008759708232096647362;
    C.G[0][4][4] = beta;
    C.G[0][5][0] = beta;
    C.G[0][5][2] = 0.9264299099302395700444874096601015328;
    C.G[0][5][4] = -1.080229692192928069168516586450436797;
    C.G[0][6][0] = -beta;
    C.G[0][6][6] = beta;
    Some(C)
}

/* ARKODE_IMEX_MRI_GARK3b: R. Chinomona & D. Reynolds SINUM 43(5):A3082-A3113, 2021 */
fn mri_table_imex_gark3b() -> Option<MRIStepCoupling> {
    let mut C = MRIStepCoupling_Alloc(1, 8, MRISTEP_IMEX)?;
    let beta: f64 = 0.4358665215084589994160194511935568425;

    C.q = 3;
    C.p = 0;

    C.c[1] = beta;
    C.c[2] = beta;
    C.c[3] = 0.7179332607542294997080097255967784213;
    C.c[4] = 0.7179332607542294997080097255967784213;
    C.c[5] = ONE;
    C.c[6] = ONE;
    C.c[7] = ONE;

    C.W[0][1][0] = beta;
    C.W[0][3][0] = -0.1750145285570467590610670000018749059;
    C.W[0][3][2] = 0.4570812678028172593530572744050964846;
    C.W[0][4][0] = 0.06042689307721552209333459437020635774;
    C.W[0][4][2] = -0.06042689307721552209333459437020635774;
    C.W[0][5][0] = 0.1195213959425454440038786034027936869;
    C.W[0][5][2] = -1.84372522668966191789853395029629765;
    C.W[0][5][4] = 2.006270569992886974186645621296725542;
    C.W[0][6][0] = -0.5466585780430528451745431084418669343;
    C.W[0][6][2] = 2.0;
    C.W[0][6][4] = -1.453341421956947154825456891558133066;
    C.W[0][7][0] = 0.105858296071879638722377459477184953;
    C.W[0][7][2] = 0.655567501140070250975288954324730635;
    C.W[0][7][4] = -1.197292318720408889113685864995472431;
    C.W[0][7][6] = beta;

    C.G[0][1][0] = beta;
    C.G[0][2][0] = -beta;
    C.G[0][2][2] = beta;
    C.G[0][3][0] = 0.0414273753564414837153799230278275639;
    C.G[0][3][2] = 0.2406393638893290165766103513753940148;
    C.G[0][4][0] = -0.0414273753564414837153799230278275639;
    C.G[0][4][2] = -0.3944391461520175157006395281657292786;
    C.G[0][4][4] = beta;
    C.G[0][5][0] = 0.1123373143006047802633543416889605123;
    C.G[0][5][2] = 1.051807513648115027700693049638099167;
    C.G[0][5][4] = -0.8820780887029493076720571169238381009;
    C.G[0][6][0] = -0.1123373143006047802633543416889605123;
    C.G[0][6][2] = -0.1253776037178754576562056399779976346;
    C.G[0][6][4] = -0.1981516034899787614964594695265986957;
    C.G[0][6][6] = beta;
    Some(C)
}

/* ARKODE_IMEX_MRI_GARK4: R. Chinomona & D. Reynolds SINUM 43(5):A3082-A3113, 2021 */
fn mri_table_imex_gark4() -> Option<MRIStepCoupling> {
    let mut C = MRIStepCoupling_Alloc(2, 12, MRISTEP_IMEX)?;

    C.q = 4;
    C.p = 0;

    C.c[1] = 0.5;
    C.c[2] = 0.5;
    C.c[3] = 0.625;
    C.c[4] = 0.625;
    C.c[5] = 0.75;
    C.c[6] = 0.75;
    C.c[7] = 0.875;
    C.c[8] = 0.875;
    C.c[9] = ONE;
    C.c[10] = ONE;
    C.c[11] = ONE;

    C.W[0][1][0] = 0.5;
    C.W[0][3][0] = -1.91716534363662868878172216064946905;
    C.W[0][3][2] = 2.04216534363662868878172216064946905;
    C.W[0][4][0] = -0.4047510318011059426979159070469904691;
    C.W[0][4][2] = 0.4047510318011059426979159070469904691;
    C.W[0][5][0] = 11.45146602249221636665698028602631728;
    C.W[0][5][2] = -30.21075747526504271440647815573950607;
    C.W[0][5][4] = 18.88429145277282634774949786971318879;
    C.W[0][6][0] = -0.7090335647602614506847116729463301439;
    C.W[0][6][2] = 1.03030720858751876652616190884004718;
    C.W[0][6][4] = -0.3212736438272573158414502358937170357;
    C.W[0][7][0] = -29.99548716455828439840910684944199275;
    C.W[0][7][2] = 37.60598277499180180536489685624385701;
    C.W[0][7][4] = 0.3212736438272573158414502358937170357;
    C.W[0][7][6] = -7.806769254260774722797240242695581295;
    C.W[0][8][0] = 3.104665054272962116338769391849124223;
    C.W[0][8][2] = -2.430325019757162297132065927415566359;
    C.W[0][8][4] = -1.905479301151524635219201659483842131;
    C.W[0][8][6] = 1.231139266635724816012498195050284266;
    C.W[0][9][0] = -2.424429547752047869875875914355514008;
    C.W[0][9][2] = 2.430325019757162297132065927415566359;
    C.W[0][9][4] = 1.905479301151524635219201659483842131;
    C.W[0][9][6] = -1.231139266635724816012498195050284266;
    C.W[0][9][8] = -0.555235506520914246462893477493610215;
    C.W[0][10][0] = -0.01044135044479748590294518945165354204;
    C.W[0][10][2] = 0.07260303614655074505152104505488141613;
    C.W[0][10][4] = -0.1288275951677260952239454098576424313;
    C.W[0][10][6] = 0.1129355350093823566139440107122154084;
    C.W[0][10][8] = -0.04626962554340952053857445645780085125;
    C.W[0][11][0] = -0.8108522787762101328175789228607932098;
    C.W[0][11][2] = 0.2560073199220492435001562192140882299;
    C.W[0][11][4] = 0.8068294072697527893665866422787819475;
    C.W[0][11][6] = -0.4557148228721823795105894821742761164;
    C.W[0][11][8] = -0.04626962554340952053857445645780085125;
    C.W[0][11][10] = 0.25;

    C.W[1][3][0] = 4.084330687273257377563444321298938099;
    C.W[1][3][2] = -4.084330687273257377563444321298938099;
    C.W[1][5][0] = -21.84342998138222084791812875795865363;
    C.W[1][5][2] = 59.61201288692787354341712449738503121;
    C.W[1][5][4] = -37.76858290554565269549899573942637758;
    C.W[1][7][0] = 61.65904145863709169818763704477664579;
    C.W[1][7][2] = -77.27257996715864114378211753016780838;
    C.W[1][7][6] = 15.61353850852154944559448048539116259;
    C.W[1][9][0] = -1.11047101304182849292578695498722043;
    C.W[1][9][8] = 1.11047101304182849292578695498722043;

    C.G[0][1][0] = 0.5;
    C.G[0][2][0] = -0.25;
    C.G[0][2][2] = 0.25;
    C.G[0][3][0] = -3.977281248108488183067033851462278892;
    C.G[0][3][2] = 4.102281248108488183067033851462278892;
    C.G[0][4][0] = -0.06905388741401691232724147084809374064;
    C.G[0][4][2] = -0.1809461125859830876727585291519062594;
    C.G[0][4][4] = 0.25;
    C.G[0][5][0] = -1.761767663757920528863378964822412405;
    C.G[0][5][2] = 2.694524698377298610155338150791461384;
    C.G[0][5][4] = -0.8077570346193780812919591859690489783;
    C.G[0][6][0] = 0.5558721791553969487305081009588084962;
    C.G[0][6][2] = -0.6799140501579995013958501527883486949;
    C.G[0][6][4] = -0.1259581289973974473346579481704598013;
    C.G[0][6][6] = 0.25;
    C.G[0][7][0] = -5.840176028724955954446426657541065113;
    C.G[0][7][2] = 8.174456684291915089191270805710716374;
    C.G[0][7][4] = 0.1259581289973974473346579481704598013;
    C.G[0][7][6] = -2.335238784564356582079502096340111063;
    C.G[0][8][0] = -1.906792645167811808094759305036052304;
    C.G[0][8][2] = -1.547057811385123933632984579249388443;
    C.G[0][8][4] = 4.129888013149350305954491738020313225;
    C.G[0][8][6] = -0.9260375565964145642267478537348724775;
    C.G[0][8][8] = 0.25;
    C.G[0][9][0] = 3.337028151688726054557652782529662519;
    C.G[0][9][2] = 1.547057811385123933632984579249388443;
    C.G[0][9][4] = -4.129888013149350305954491738020313225;
    C.G[0][9][6] = 0.9260375565964145642267478537348724775;
    C.G[0][9][8] = -1.555235506520914246462893477493610215;
    C.G[0][10][0] = -0.8212936292210076187205241123124467518;
    C.G[0][10][2] = 0.328610356068599988551677264268969646;
    C.G[0][10][4] = 0.6780018121020266941426412324211395162;
    C.G[0][10][6] = -0.3427792878628000228966454714620607079;
    C.G[0][10][8] = -0.0925392510868190410771489129156017025;
    C.G[0][10][10] = 0.25;

    C.G[1][3][0] = 8.704562496216976366134067702924557783;
    C.G[1][3][2] = -8.704562496216976366134067702924557783;
    C.G[1][5][0] = 3.911643102343874882381240871341012292;
    C.G[1][5][2] = -5.027157171582631044965159243279110249;
    C.G[1][5][4] = 1.115514069238756162583918371938097957;
    C.G[1][7][0] = 10.81860769913911801143183711316451323;
    C.G[1][7][2] = -14.98908526826783117559084130584473536;
    C.G[1][7][6] = 4.170477569128713164159004192680222125;
    C.G[1][9][0] = -2.61047101304182849292578695498722043;
    C.G[1][9][8] = 2.61047101304182849292578695498722043;
    Some(C)
}

/* ARKODE_IMEX_MRI_SR21: A.C. Fish, D.R. Reynolds, S.B. Roberts, JCAM 438:115534, 2024 */
fn mri_table_imex_sr21() -> Option<MRIStepCoupling> {
    let mut C = MRIStepCoupling_Alloc(1, 4, MRISTEP_SR)?;

    C.q = 2;
    C.p = 1;

    C.c[1] = 3.0 / 5.0;
    C.c[2] = 4.0 / 15.0;
    C.c[3] = ONE;

    C.W[0][1][0] = 3.0 / 5.0;
    C.W[0][2][0] = 14.0 / 165.0;
    C.W[0][2][1] = 2.0 / 11.0;
    C.W[0][3][0] = -13.0 / 54.0;
    C.W[0][3][1] = 137.0 / 270.0;
    C.W[0][3][2] = 11.0 / 15.0;
    C.W[0][4][0] = -0.25;
    C.W[0][4][1] = 0.5;
    C.W[0][4][2] = 0.75;

    C.G[0][1][0] = -11.0 / 23.0;
    C.G[0][1][1] = 11.0 / 23.0;
    C.G[0][2][0] = -6692.0 / 52371.0;
    C.G[0][2][1] = -18355.0 / 52371.0;
    C.G[0][2][2] = 11.0 / 23.0;
    C.G[0][3][0] = 11621.0 / 90666.0;
    C.G[0][3][1] = -215249.0 / 226665.0;
    C.G[0][3][2] = 17287.0 / 50370.0;
    C.G[0][3][3] = 11.0 / 23.0;
    C.G[0][4][0] = -31.0 / 12.0;
    C.G[0][4][1] = -ONE / 6.0;
    C.G[0][4][2] = 11.0 / 4.0;
    Some(C)
}

/* ARKODE_IMEX_MRI_SR32: A.C. Fish, D.R. Reynolds, S.B. Roberts, JCAM 438:115534, 2024 */
fn mri_table_imex_sr32() -> Option<MRIStepCoupling> {
    let mut C = MRIStepCoupling_Alloc(2, 5, MRISTEP_SR)?;

    C.q = 3;
    C.p = 2;

    C.c[1] = 23.0 / 34.0;
    C.c[2] = 4.0 / 5.0;
    C.c[3] = 17.0 / 15.0;
    C.c[4] = ONE;

    C.W[0][1][0] = 23.0 / 34.0;
    C.W[0][2][0] = 71.0 / 70.0;
    C.W[0][2][1] = -3.0 / 14.0;
    C.W[0][3][0] = 124.0 / 1155.0;
    C.W[0][3][1] = 4.0 / 7.0;
    C.W[0][3][2] = 5.0 / 11.0;
    C.W[0][4][0] = 162181.0 / 187680.0;
    C.W[0][4][1] = 119.0 / 1380.0;
    C.W[0][4][2] = 11.0 / 32.0;
    C.W[0][4][3] = -5.0 / 17.0;
    C.W[0][5][0] = 76355.0 / 74834.0;
    C.W[0][5][1] = -46.0 / 31.0;
    C.W[0][5][2] = 67.0 / 34.0;
    C.W[0][5][3] = -36.0 / 71.0;

    C.W[1][2][0] = -14453.0 / 63825.0;
    C.W[1][2][1] = 14453.0 / 63825.0;
    C.W[1][3][0] = -2101267877.0 / 1206582300.0;
    C.W[1][3][1] = 2476735438.0 / 301645575.0;
    C.W[1][3][2] = -13575085.0 / 2098404.0;
    C.W[1][4][0] = -762580446799.0 / 588660102960.0;
    C.W[1][4][1] = 11083240219.0 / 4328383110.0;
    C.W[1][4][2] = -211274129.0 / 100368304.0;
    C.W[1][4][3] = 89562055.0 / 106641323.0;
    C.W[1][5][0] = -3732974.0 / 2278035.0;
    C.W[1][5][1] = 13857574.0 / 2278035.0;
    C.W[1][5][2] = -52.0 / 9.0;
    C.W[1][5][3] = 4.0 / 3.0;

    C.G[0][1][0] = -4.0 / 7.0;
    C.G[0][1][1] = 4.0 / 7.0;
    C.G[0][2][0] = -2707004.0 / 3127425.0;
    C.G[0][2][1] = 919904.0 / 3127425.0;
    C.G[0][2][2] = 4.0 / 7.0;
    C.G[0][3][0] = 852879271.0 / 703839675.0;
    C.G[0][3][1] = -1575000496.0 / 703839675.0;
    C.G[0][3][2] = 5.0 / 11.0;
    C.G[0][3][3] = 4.0 / 7.0;
    C.G[0][4][0] = 43136869.0 / 2019912118.0;
    C.G[0][4][1] = -73810600.0 / 1009956059.0;
    C.G[0][4][2] = -17653551.0 / 87822266.0;
    C.G[0][4][3] = -13993902.0 / 43911133.0;
    C.G[0][4][4] = 4.0 / 7.0;
    C.G[0][5][0] = -179.0 / 4140.0;
    C.G[0][5][1] = 799.0 / 14490.0;
    C.G[0][5][2] = ONE / 14.0;
    C.G[0][5][3] = -ONE / 12.0;
    Some(C)
}

/* ARKODE_IMEX_MRI_SR43: A.C. Fish, D.R. Reynolds, S.B. Roberts, arXiv:2301.00865, 2023 */
fn mri_table_imex_sr43() -> Option<MRIStepCoupling> {
    let mut C = MRIStepCoupling_Alloc(2, 7, MRISTEP_SR)?;

    C.q = 4;
    C.p = 3;

    C.c[1] = ONE / 4.0;
    C.c[2] = 3.0 / 4.0;
    C.c[3] = 11.0 / 20.0;
    C.c[4] = ONE / 2.0;
    C.c[5] = ONE;
    C.c[6] = ONE;

    C.W[0][1][0] = ONE / 4.0;
    C.W[0][2][0] = 9.0 / 8.0;
    C.W[0][2][1] = -3.0 / 8.0;
    C.W[0][3][0] = 187.0 / 2340.0;
    C.W[0][3][1] = 7.0 / 9.0;
    C.W[0][3][2] = -4.0 / 13.0;
    C.W[0][4][0] = 64.0 / 165.0;
    C.W[0][4][1] = ONE / 6.0;
    C.W[0][4][2] = -3.0 / 5.0;
    C.W[0][4][3] = 6.0 / 11.0;
    C.W[0][5][0] = 1816283.0 / 549120.0;
    C.W[0][5][1] = -2.0 / 9.0;
    C.W[0][5][2] = -4.0 / 11.0;
    C.W[0][5][3] = -ONE / 6.0;
    C.W[0][5][4] = -2561809.0 / 1647360.0;
    C.W[0][6][1] = 7.0 / 11.0;
    C.W[0][6][2] = -2203.0 / 264.0;
    C.W[0][6][3] = 10825.0 / 792.0;
    C.W[0][6][4] = -85.0 / 12.0;
    C.W[0][6][5] = 841.0 / 396.0;
    C.W[0][7][0] = ONE / 400.0;
    C.W[0][7][1] = 49.0 / 12.0;
    C.W[0][7][2] = 43.0 / 6.0;
    C.W[0][7][3] = -7.0 / 10.0;
    C.W[0][7][4] = -85.0 / 12.0;
    C.W[0][7][5] = -2963.0 / 1200.0;

    C.W[1][2][0] = -11.0 / 4.0;
    C.W[1][2][1] = 11.0 / 4.0;
    C.W[1][3][0] = -1228.0 / 2925.0;
    C.W[1][3][1] = -92.0 / 225.0;
    C.W[1][3][2] = 808.0 / 975.0;
    C.W[1][4][0] = -2572.0 / 2805.0;
    C.W[1][4][1] = 167.0 / 255.0;
    C.W[1][4][2] = 199.0 / 136.0;
    C.W[1][4][3] = -1797.0 / 1496.0;
    C.W[1][5][0] = -1816283.0 / 274560.0;
    C.W[1][5][1] = 253.0 / 36.0;
    C.W[1][5][2] = -23.0 / 44.0;
    C.W[1][5][3] = 76.0 / 3.0;
    C.W[1][5][4] = -20775791.0 / 823680.0;
    C.W[1][6][1] = 107.0 / 132.0;
    C.W[1][6][2] = 1289.0 / 88.0;
    C.W[1][6][3] = -9275.0 / 792.0;
    C.W[1][6][5] = -371.0 / 99.0;
    C.W[1][7][0] = -ONE / 200.0;
    C.W[1][7][1] = -137.0 / 24.0;
    C.W[1][7][2] = -235.0 / 16.0;
    C.W[1][7][3] = 1237.0 / 80.0;
    C.W[1][7][5] = 2963.0 / 600.0;

    C.G[0][1][0] = -ONE / 4.0;
    C.G[0][1][1] = ONE / 4.0;
    C.G[0][2][0] = ONE / 4.0;
    C.G[0][2][1] = -ONE / 2.0;
    C.G[0][2][2] = ONE / 4.0;
    C.G[0][3][0] = 13.0 / 100.0;
    C.G[0][3][1] = -7.0 / 30.0;
    C.G[0][3][2] = -11.0 / 75.0;
    C.G[0][3][3] = ONE / 4.0;
    C.G[0][4][0] = 6.0 / 85.0;
    C.G[0][4][1] = -301.0 / 1360.0;
    C.G[0][4][2] = -99.0 / 544.0;
    C.G[0][4][3] = 45.0 / 544.0;
    C.G[0][4][4] = ONE / 4.0;
    C.G[0][5][1] = -9.0 / 4.0;
    C.G[0][5][2] = -19.0 / 48.0;
    C.G[0][5][3] = -75.0 / 16.0;
    C.G[0][5][4] = 85.0 / 12.0;
    C.G[0][5][5] = ONE / 4.0;
    Some(C)
}

/* ARKODE_MERK21: A.C. Fish, D.R. Reynolds, S.B. Roberts, JCAM 438:115534, 2024 (embedding unpublished) */
fn mri_table_merk21() -> Option<MRIStepCoupling> {
    let mut C = MRIStepCoupling_Alloc(2, 3, MRISTEP_MERK)?;
    let c2: f64 = 0.5;

    C.q = 2;
    C.p = 1;
    C.ngroup = 2;
    C.group[0][0] = 1;
    C.group[0][1] = 3;
    C.group[1][0] = 2;

    C.c[1] = c2;
    C.c[2] = ONE;

    C.W[0][1][0] = ONE;
    C.W[0][2][0] = ONE;
    C.W[0][3][0] = ONE;

    C.W[1][2][0] = -ONE / c2;
    C.W[1][2][1] = ONE / c2;
    Some(C)
}

/* ARKODE_MERK32: A.C. Fish, D.R. Reynolds, S.B. Roberts, JCAM 438:115534, 2024 (embedding unpublished) */
fn mri_table_merk32() -> Option<MRIStepCoupling> {
    let mut C = MRIStepCoupling_Alloc(2, 4, MRISTEP_MERK)?;
    let c2: f64 = 0.5;

    C.q = 3;
    C.p = 2;
    C.ngroup = 3;
    C.group[0][0] = 1;
    C.group[1][0] = 2;
    C.group[1][1] = 4;
    C.group[2][0] = 3;

    C.c[1] = c2;
    C.c[2] = 2.0 / 3.0;
    C.c[3] = ONE;

    C.W[0][1][0] = ONE;
    C.W[0][2][0] = ONE;
    C.W[0][3][0] = ONE;
    C.W[0][4][0] = ONE;

    C.W[1][2][0] = -ONE / c2;
    C.W[1][2][1] = ONE / c2;
    C.W[1][3][0] = -1.5;
    C.W[1][3][2] = 1.5;
    C.W[1][4][0] = -ONE / c2;
    C.W[1][4][1] = ONE / c2;
    Some(C)
}

/* ARKODE_MERK43: A.C. Fish, D.R. Reynolds, S.B. Roberts, JCAM 438:115534, 2024 (embedding unpublished) */
fn mri_table_merk43() -> Option<MRIStepCoupling> {
    let mut C = MRIStepCoupling_Alloc(3, 7, MRISTEP_MERK)?;
    let c2: f64 = 0.5;
    let c3: f64 = 0.5;
    let c4: f64 = ONE / 3.0;
    let c5: f64 = 5.0 / 6.0;
    let c6: f64 = ONE / 3.0;

    C.q = 4;
    C.p = 3;
    C.ngroup = 4;
    C.group[0][0] = 1;
    C.group[1][0] = 3;
    C.group[1][1] = 2;
    C.group[2][0] = 5;
    C.group[2][1] = 4;
    C.group[2][2] = 7;
    C.group[3][0] = 6;

    C.c[1] = c2;
    C.c[2] = c3;
    C.c[3] = c4;
    C.c[4] = c5;
    C.c[5] = c6;
    C.c[6] = ONE;

    C.W[0][1][0] = ONE;
    C.W[0][2][0] = ONE;
    C.W[0][3][0] = ONE;
    C.W[0][4][0] = ONE;
    C.W[0][5][0] = ONE;
    C.W[0][6][0] = ONE;
    C.W[0][7][0] = ONE;

    C.W[1][2][0] = -ONE / c2;
    C.W[1][2][1] = ONE / c2;
    C.W[1][3][0] = -ONE / c2;
    C.W[1][3][1] = ONE / c2;
    C.W[1][4][0] = c4 / c3 / (c3 - c4) - c3 / c4 / (c3 - c4);
    C.W[1][4][2] = -c4 / c3 / (c3 - c4);
    C.W[1][4][3] = c3 / c4 / (c3 - c4);
    C.W[1][5][0] = c4 / c3 / (c3 - c4) - c3 / c4 / (c3 - c4);
    C.W[1][5][2] = -c4 / c3 / (c3 - c4);
    C.W[1][5][3] = c3 / c4 / (c3 - c4);
    C.W[1][6][0] = c6 / c5 / (c5 - c6) - c5 / c6 / (c5 - c6);
    C.W[1][6][4] = -c6 / c5 / (c5 - c6);
    C.W[1][6][5] = c5 / c6 / (c5 - c6);
    C.W[1][7][0] = c4 / c3 / (c3 - c4) - c3 / c4 / (c3 - c4);
    C.W[1][7][2] = -c4 / c3 / (c3 - c4);
    C.W[1][7][3] = c3 / c4 / (c3 - c4);

    C.W[2][4][0] = ONE / c4 / (c3 - c4) - ONE / c3 / (c3 - c4);
    C.W[2][4][2] = ONE / c3 / (c3 - c4);
    C.W[2][4][3] = -ONE / c4 / (c3 - c4);
    C.W[2][5][0] = ONE / c4 / (c3 - c4) - ONE / c3 / (c3 - c4);
    C.W[2][5][2] = ONE / c3 / (c3 - c4);
    C.W[2][5][3] = -ONE / c4 / (c3 - c4);
    C.W[2][6][0] = ONE / c6 / (c5 - c6) - ONE / c5 / (c5 - c6);
    C.W[2][6][4] = ONE / c5 / (c5 - c6);
    C.W[2][6][5] = -ONE / c6 / (c5 - c6);
    C.W[2][7][0] = ONE / c4 / (c3 - c4) - ONE / c3 / (c3 - c4);
    C.W[2][7][2] = ONE / c3 / (c3 - c4);
    C.W[2][7][3] = -ONE / c4 / (c3 - c4);

    Some(C)
}

/* ARKODE_MERK54: A.C. Fish, D.R. Reynolds, S.B. Roberts, JCAM 438:115534, 2024 (embedding unpublished) */
fn mri_table_merk54() -> Option<MRIStepCoupling> {
    let mut C = MRIStepCoupling_Alloc(4, 11, MRISTEP_MERK)?;
    let c2: f64 = 0.5;
    let c3: f64 = 0.5;
    let c4: f64 = ONE / 3.0;
    let c5: f64 = 0.5;
    let c6: f64 = ONE / 3.0;
    let c7: f64 = 0.25;
    let c8: f64 = 0.7;
    let c9: f64 = 0.5;
    let c10: f64 = 2.0 / 3.0;
    let a2: f64 = ONE / c2;
    let a3: f64 = c4 / c3 / (c4 - c3);
    let a4: f64 = c3 / c4 / (c3 - c4);
    let a5: f64 = c6 * c7 / c5 / (c5 - c6) / (c5 - c7);
    let a6: f64 = c5 * c7 / c6 / (c6 - c5) / (c6 - c7);
    let a7: f64 = c5 * c6 / c7 / (c7 - c5) / (c7 - c6);
    let a8: f64 = c9 * c10 / c8 / (c8 - c9) / (c8 - c10);
    let a9: f64 = c8 * c10 / c9 / (c9 - c8) / (c9 - c10);
    let a10: f64 = c8 * c9 / c10 / (c10 - c8) / (c10 - c9);
    let b3: f64 = ONE / c3 / (c3 - c4);
    let b4: f64 = ONE / c4 / (c3 - c4);
    let b5: f64 = (c6 + c7) / c5 / (c5 - c6) / (c5 - c7);
    let b6: f64 = (c5 + c7) / c6 / (c6 - c5) / (c6 - c7);
    let b7: f64 = (c5 + c6) / c7 / (c7 - c5) / (c7 - c6);
    let b8: f64 = (c9 + c10) / c8 / (c8 - c9) / (c8 - c10);
    let b9: f64 = (c8 + c10) / c9 / (c9 - c8) / (c9 - c10);
    let b10: f64 = (c8 + c9) / c10 / (c10 - c8) / (c10 - c9);
    let g5: f64 = ONE / c5 / (c5 - c6) / (c5 - c7);
    let g6: f64 = ONE / c6 / (c6 - c5) / (c6 - c7);
    let g7: f64 = ONE / c7 / (c7 - c5) / (c7 - c6);
    let g8: f64 = ONE / c8 / (c8 - c9) / (c8 - c10);
    let g9: f64 = ONE / c9 / (c9 - c8) / (c9 - c10);
    let g10: f64 = ONE / c10 / (c10 - c8) / (c10 - c9);

    C.q = 5;
    C.p = 4;
    C.ngroup = 5;
    C.group[0][0] = 1;
    C.group[1][0] = 3;
    C.group[1][1] = 2;
    C.group[2][0] = 6;
    C.group[2][1] = 5;
    C.group[2][2] = 4;
    C.group[3][0] = 8;
    C.group[3][1] = 9;
    C.group[3][2] = 7;
    C.group[3][3] = 11;
    C.group[4][0] = 10;

    C.c[1] = c2;
    C.c[2] = c3;
    C.c[3] = c4;
    C.c[4] = c5;
    C.c[5] = c6;
    C.c[6] = c7;
    C.c[7] = c8;
    C.c[8] = c9;
    C.c[9] = c10;
    C.c[10] = ONE;

    for i in 1..=11usize {
        C.W[0][i][0] = ONE;
    }

    C.W[1][2][0] = -a2;
    C.W[1][2][1] = a2;
    C.W[1][3][0] = -a2;
    C.W[1][3][1] = a2;
    for i in [4usize, 5, 6] {
        C.W[1][i][0] = -(a3 + a4);
        C.W[1][i][2] = a3;
        C.W[1][i][3] = a4;
    }
    for i in [7usize, 8, 9, 11] {
        C.W[1][i][0] = -(a5 + a6 + a7);
        C.W[1][i][4] = a5;
        C.W[1][i][5] = a6;
        C.W[1][i][6] = a7;
    }
    C.W[1][10][0] = -(a8 + a9 + a10);
    C.W[1][10][7] = a8;
    C.W[1][10][8] = a9;
    C.W[1][10][9] = a10;

    for i in [4usize, 5, 6] {
        C.W[2][i][0] = b4 - b3;
        C.W[2][i][2] = b3;
        C.W[2][i][3] = -b4;
    }
    for i in [7usize, 8, 9, 11] {
        C.W[2][i][0] = b5 + b6 + b7;
        C.W[2][i][4] = -b5;
        C.W[2][i][5] = -b6;
        C.W[2][i][6] = -b7;
    }
    C.W[2][10][0] = b8 + b9 + b10;
    C.W[2][10][7] = -b8;
    C.W[2][10][8] = -b9;
    C.W[2][10][9] = -b10;

    for i in [7usize, 8, 9, 11] {
        C.W[3][i][0] = -(g5 + g6 + g7);
        C.W[3][i][4] = g5;
        C.W[3][i][5] = g6;
        C.W[3][i][6] = g7;
    }
    C.W[3][10][0] = -(g8 + g9 + g10);
    C.W[3][10][7] = g8;
    C.W[3][10][8] = g9;
    C.W[3][10][9] = g10;

    Some(C)
}

/*---------------------------------------------------------------
  Returns MRIStepCoupling table structure for pre-set MRI methods.

  Input:  imeth -- integer key for the desired method
  ---------------------------------------------------------------*/
pub fn MRIStepCoupling_LoadTable(method: ARKODE_MRITableID) -> Option<MRIStepCoupling> {
    match method {
        ARKODE_MRI_NONE => None,
        ARKODE_MRI_GARK_FORWARD_EULER => mri_table_from_erk(ARKODE_FORWARD_EULER_1_1),
        ARKODE_MRI_GARK_RALSTON2 => mri_table_from_erk(ARKODE_RALSTON_EULER_2_1_2),
        ARKODE_MIS_KW3 => mri_table_from_erk(ARKODE_KNOTH_WOLKE_3_3),
        ARKODE_MRI_GARK_ERK22a => mri_table_from_erk(ARKODE_EXPLICIT_MIDPOINT_EULER_2_1_2),
        ARKODE_MRI_GARK_ERK22b => mri_table_from_erk(ARKODE_HEUN_EULER_2_1_2),
        ARKODE_MRI_GARK_ERK33a => mri_table_erk33a(),
        ARKODE_MRI_GARK_RALSTON3 => mri_table_ralston3(),
        ARKODE_MRI_GARK_ERK45a => mri_table_erk45a(),
        ARKODE_MRI_GARK_BACKWARD_EULER => mri_table_backward_euler(),
        ARKODE_MRI_GARK_IRK21a => mri_table_irk21a(),
        ARKODE_MRI_GARK_IMPLICIT_MIDPOINT => mri_table_implicit_midpoint(),
        ARKODE_MRI_GARK_ESDIRK34a => mri_table_esdirk34a(),
        ARKODE_MRI_GARK_ESDIRK46a => mri_table_esdirk46a(),
        ARKODE_IMEX_MRI_GARK_EULER => mri_table_imex_euler(),
        ARKODE_IMEX_MRI_GARK_TRAPEZOIDAL => mri_table_imex_trapezoidal(),
        ARKODE_IMEX_MRI_GARK_MIDPOINT => mri_table_imex_midpoint(),
        ARKODE_IMEX_MRI_GARK3a => mri_table_imex_gark3a(),
        ARKODE_IMEX_MRI_GARK3b => mri_table_imex_gark3b(),
        ARKODE_IMEX_MRI_GARK4 => mri_table_imex_gark4(),
        ARKODE_IMEX_MRI_SR21 => mri_table_imex_sr21(),
        ARKODE_IMEX_MRI_SR32 => mri_table_imex_sr32(),
        ARKODE_IMEX_MRI_SR43 => mri_table_imex_sr43(),
        ARKODE_MERK21 => mri_table_merk21(),
        ARKODE_MERK32 => mri_table_merk32(),
        ARKODE_MERK43 => mri_table_merk43(),
        ARKODE_MERK54 => mri_table_merk54(),
        _ => {
            arkProcessError(
                None,
                ARK_ILL_INPUT,
                line!(),
                "MRIStepCoupling_LoadTable",
                file!(),
                "Unknown coupling table",
            );
            None
        }
    }
}

/*---------------------------------------------------------------
  Returns MRIStepCoupling table structure for pre-set MRI methods.

  Input:  method -- string key for the desired method
  ---------------------------------------------------------------*/
pub fn MRIStepCoupling_LoadTableByName(method: &str) -> Option<MRIStepCoupling> {
    match method {
        "ARKODE_MRI_NONE" => None,
        "ARKODE_MRI_GARK_FORWARD_EULER" => mri_table_from_erk(ARKODE_FORWARD_EULER_1_1),
        "ARKODE_MRI_GARK_RALSTON2" => mri_table_from_erk(ARKODE_RALSTON_EULER_2_1_2),
        "ARKODE_MIS_KW3" => mri_table_from_erk(ARKODE_KNOTH_WOLKE_3_3),
        "ARKODE_MRI_GARK_ERK22a" => mri_table_from_erk(ARKODE_EXPLICIT_MIDPOINT_EULER_2_1_2),
        "ARKODE_MRI_GARK_ERK22b" => mri_table_from_erk(ARKODE_HEUN_EULER_2_1_2),
        "ARKODE_MRI_GARK_ERK33a" => mri_table_erk33a(),
        "ARKODE_MRI_GARK_RALSTON3" => mri_table_ralston3(),
        "ARKODE_MRI_GARK_ERK45a" => mri_table_erk45a(),
        "ARKODE_MRI_GARK_BACKWARD_EULER" => mri_table_backward_euler(),
        "ARKODE_MRI_GARK_IRK21a" => mri_table_irk21a(),
        "ARKODE_MRI_GARK_IMPLICIT_MIDPOINT" => mri_table_implicit_midpoint(),
        "ARKODE_MRI_GARK_ESDIRK34a" => mri_table_esdirk34a(),
        "ARKODE_MRI_GARK_ESDIRK46a" => mri_table_esdirk46a(),
        "ARKODE_IMEX_MRI_GARK_EULER" => mri_table_imex_euler(),
        "ARKODE_IMEX_MRI_GARK_TRAPEZOIDAL" => mri_table_imex_trapezoidal(),
        "ARKODE_IMEX_MRI_GARK_MIDPOINT" => mri_table_imex_midpoint(),
        "ARKODE_IMEX_MRI_GARK3a" => mri_table_imex_gark3a(),
        "ARKODE_IMEX_MRI_GARK3b" => mri_table_imex_gark3b(),
        "ARKODE_IMEX_MRI_GARK4" => mri_table_imex_gark4(),
        "ARKODE_IMEX_MRI_SR21" => mri_table_imex_sr21(),
        "ARKODE_IMEX_MRI_SR32" => mri_table_imex_sr32(),
        "ARKODE_IMEX_MRI_SR43" => mri_table_imex_sr43(),
        "ARKODE_MERK21" => mri_table_merk21(),
        "ARKODE_MERK32" => mri_table_merk32(),
        "ARKODE_MERK43" => mri_table_merk43(),
        "ARKODE_MERK54" => mri_table_merk54(),
        _ => {
            arkProcessError(
                None,
                ARK_ILL_INPUT,
                line!(),
                "MRIStepCoupling_LoadTableByName",
                file!(),
                "Unknown coupling table",
            );
            None
        }
    }
}
