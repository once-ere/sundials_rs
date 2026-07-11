/* -----------------------------------------------------------------
 * Translation of sundials-7.7.0/src/arkode/arkode_butcher.c
 * (+ include/arkode/arkode_butcher.h).
 *
 * PART I: the ARKodeButcherTable object (Alloc/Create/Copy/Space/
 * Free/Write/IsStifflyAccurate). PART II: CheckOrder / CheckARKOrder
 * with the order-condition helper family (rowsum, order1..order6s,
 * __ButcherSimplifyingAssumptions).
 *
 * Modeling: C's ARKodeButcherTable is a heap pointer that may be
 * NULL — Rust functions return Option<ARKodeButcherTable>.
 * `sunrealtype** A` → Vec<Vec<f64>> (row-major, stages×stages);
 * `d` (embedding) is Some only for embedded tables, exactly C's
 * NULL test. ARKodeButcherTable_Free = drop.
 * SUN_FORMAT_E is "% .15e": fmt_e with a leading space for
 * non-negative values (established pattern).
 *
 * PART II notes:
 *  - C's method-order and embedding-order blocks in CheckOrder (and
 *    in CheckARKOrder) are literal duplicates differing only in the
 *    printed label ("method"/"embedding") and the weight vector
 *    (b or d); they are factored into arkode_butcher_check_order /
 *    arkode_butcher_check_ark_order helpers that reproduce the C
 *    message text byte-for-byte (condition sequence, short-circuit
 *    structure and alltrue accumulation preserved exactly).
 *  - CheckARKOrder quirk preserved: C sets d[1] = B1->d (not B2->d),
 *    so the embedding check runs only if B1 is embedded and tests
 *    B1's d twice (arkode_butcher.c:1168).
 *  - The mv/vv/vp/dot utilities keep C's s<1 failure returns; the
 *    C NULL-pointer failure branches are unrepresentable.
 * -----------------------------------------------------------------*/

use crate::sundials_math::{SUNRabs, SUNRpowerI, SUNRsqrt};
use crate::sundials_types::SUN_UNIT_ROUNDOFF;
use crate::sundials_utils::fmt_e;

/* tolerance for checking order conditions */
#[allow(non_snake_case)]
fn TOL() -> f64 {
    SUNRsqrt(SUN_UNIT_ROUNDOFF)
}

/// struct ARKodeButcherTableMem (arkode_butcher.h)
pub struct ARKodeButcherTable {
    pub q: i32,           /* method order of accuracy       */
    pub p: i32,           /* embedding order of accuracy    */
    pub stages: i32,      /* number of stages               */
    pub A: Vec<Vec<f64>>, /* Butcher table coefficients     */
    pub c: Vec<f64>,      /* canopy node coefficients       */
    pub b: Vec<f64>,      /* root node coefficients         */
    pub d: Option<Vec<f64>>, /* embedding coefficients      */
}

/* C SUN_FORMAT_E = "% .15e" */
fn fmt_e_spaced(x: f64) -> String {
    let e = fmt_e(x, 0, 15);
    if e.starts_with('-') {
        e
    } else {
        format!(" {}", e)
    }
}

/*---------------------------------------------------------------
  Routine to allocate an empty Butcher table structure
  ---------------------------------------------------------------*/
pub fn ARKodeButcherTable_Alloc(stages: i32, embedded: bool) -> Option<ARKodeButcherTable> {
    /* Check for legal 'stages' value */
    if stages < 1 {
        return None;
    }

    let s = stages as usize;
    Some(ARKodeButcherTable {
        stages,
        A: vec![vec![0.0; s]; s],
        b: vec![0.0; s],
        c: vec![0.0; s],
        d: if embedded { Some(vec![0.0; s]) } else { None },
        /* initialize order parameters */
        q: 0,
        p: 0,
    })
}

/*---------------------------------------------------------------
  Routine to allocate and fill a Butcher table structure
  ---------------------------------------------------------------*/
pub fn ARKodeButcherTable_Create(
    s: i32,
    q: i32,
    p: i32,
    c: &[f64],
    A: &[f64],
    b: &[f64],
    d: Option<&[f64]>,
) -> Option<ARKodeButcherTable> {
    /* Check for legal number of stages */
    if s < 1 {
        return None;
    }

    /* Does the table have an embedding? */
    let embedded = d.is_some();

    /* Allocate Butcher table structure */
    let mut B = ARKodeButcherTable_Alloc(s, embedded)?;

    /* set the relevant parameters */
    B.stages = s;
    B.q = q;
    B.p = p;

    let s = s as usize;
    for i in 0..s {
        B.c[i] = c[i];
        B.b[i] = b[i];
        for j in 0..s {
            B.A[i][j] = A[i * s + j];
        }
    }

    if let Some(d) = d {
        let bd = B.d.as_mut().unwrap();
        for i in 0..s {
            bd[i] = d[i];
        }
    }

    Some(B)
}

/*---------------------------------------------------------------
  Routine to copy a Butcher table structure
  ---------------------------------------------------------------*/
pub fn ARKodeButcherTable_Copy(B: &ARKodeButcherTable) -> Option<ARKodeButcherTable> {
    /* Get the number of stages */
    let s = B.stages;

    /* Does the table have an embedding? */
    let embedded = B.d.is_some();

    /* Allocate Butcher table structure */
    let mut Bcopy = ARKodeButcherTable_Alloc(s, embedded)?;

    /* set the relevant parameters */
    Bcopy.stages = B.stages;
    Bcopy.q = B.q;
    Bcopy.p = B.p;

    /* Copy Butcher table */
    let s = s as usize;
    for i in 0..s {
        Bcopy.c[i] = B.c[i];
        Bcopy.b[i] = B.b[i];
        for j in 0..s {
            Bcopy.A[i][j] = B.A[i][j];
        }
    }

    if embedded {
        let bd = Bcopy.d.as_mut().unwrap();
        let sd = B.d.as_ref().unwrap();
        for i in 0..s {
            bd[i] = sd[i];
        }
    }

    Some(Bcopy)
}

/*---------------------------------------------------------------
  Routine to query the Butcher table structure workspace size
  ---------------------------------------------------------------*/
pub fn ARKodeButcherTable_Space(B: &ARKodeButcherTable, liw: &mut i64, lrw: &mut i64) {
    /* fill outputs based on B */
    *liw = 3;
    if B.d.is_some() {
        *lrw = (B.stages as i64) * (B.stages as i64 + 3);
    } else {
        *lrw = (B.stages as i64) * (B.stages as i64 + 2);
    }
}

/*---------------------------------------------------------------
  Routine to free a Butcher table structure (ownership drop)
  ---------------------------------------------------------------*/
pub fn ARKodeButcherTable_Free(B: Option<ARKodeButcherTable>) {
    drop(B);
}

/*---------------------------------------------------------------
  Routine to print a Butcher table structure
  ---------------------------------------------------------------*/
pub fn ARKodeButcherTable_Write(B: &ARKodeButcherTable, outfile: &mut dyn std::io::Write) {
    let stages = B.stages as usize;

    let _ = write!(outfile, "  A = \n");
    for i in 0..stages {
        let _ = write!(outfile, "      ");
        for j in 0..stages {
            let _ = write!(outfile, "{}  ", fmt_e_spaced(B.A[i][j]));
        }
        let _ = write!(outfile, "\n");
    }

    let _ = write!(outfile, "  c = ");
    for i in 0..stages {
        let _ = write!(outfile, "{}  ", fmt_e_spaced(B.c[i]));
    }
    let _ = write!(outfile, "\n");

    let _ = write!(outfile, "  b = ");
    for i in 0..stages {
        let _ = write!(outfile, "{}  ", fmt_e_spaced(B.b[i]));
    }
    let _ = write!(outfile, "\n");

    if let Some(d) = B.d.as_ref() {
        let _ = write!(outfile, "  d = ");
        for i in 0..stages {
            let _ = write!(outfile, "{}  ", fmt_e_spaced(d[i]));
        }
        let _ = write!(outfile, "\n");
    }
}

pub fn ARKodeButcherTable_IsStifflyAccurate(B: &ARKodeButcherTable) -> bool {
    let stages = B.stages as usize;
    for i in 0..stages {
        if (B.b[i] - B.A[stages - 1][i]).abs() > 100.0 * SUN_UNIT_ROUNDOFF {
            return false;
        }
    }
    true
}

/* -----------------------------------------------------------------
 * PART II — order-of-accuracy analysis
 * ----------------------------------------------------------------- */

/// C `if (outfile) fprintf(outfile, ...)`.
fn outmsg(outfile: &mut Option<&mut dyn std::io::Write>, s: &str) {
    if let Some(f) = outfile.as_mut() {
        let _ = f.write_all(s.as_bytes());
    }
}

/*---------------------------------------------------------------
  Routine to determine the analytical order of accuracy for a
  specified Butcher table.  We check the analytical [necessary]
  order conditions up through order 6.  After that, we revert to
  the [sufficient] Butcher simplifying assumptions.

  Return values:
     0 (success): internal {q,p} values match analytical order
     1 (warning): internal {q,p} values are lower than analytical
        order, or method achieves maximum order possible with this
        routine and internal {q,p} are higher.
    -1 (failure): internal p and q values are higher than analytical
         order
    -2 (failure): NULL-valued B (or critical contents)

  Note: for embedded methods, if the return flags for p and q would
  differ, failure takes precedence over warning, which takes
  precedence over success.
  ---------------------------------------------------------------*/
pub fn ARKodeButcherTable_CheckOrder(
    B: &ARKodeButcherTable,
    q: &mut i32,
    p: &mut i32,
    outfile: Option<&mut dyn std::io::Write>,
) -> i32 {
    let mut outfile = outfile;
    *q = 0;
    *p = 0;

    /* verify Butcher table contents (C NULL checks unrepresentable) */
    if B.stages < 1 {
        return -2;
    }

    /* set shortcuts for Butcher table components */
    let (A, b, c, d, s) = (&B.A, &B.b, &B.c, B.d.as_ref(), B.stages);

    /* check method order */
    outmsg(&mut outfile, "ARKodeButcherTable_CheckOrder:\n");
    *q = arkode_butcher_check_all_conditions(A, b, c, s, &mut outfile, "method");

    /* check embedding order */
    if let Some(d) = d {
        outmsg(&mut outfile, "\n");
        *p = arkode_butcher_check_all_conditions(A, d, c, s, &mut outfile, "embedding");
    }

    /* compare results against stored values and return */

    /*    check failure modes first */
    if (*q < B.q) && (*q < 6) {
        return -1;
    }
    if d.is_some() && (*p < B.p) && (*p < 6) {
        return -1;
    }

    /*    check warning modes */
    if *q > B.q {
        return 1;
    }
    if d.is_some() && (*p > B.p) {
        return 1;
    }
    if (*q < B.q) && (*q >= 6) {
        return 1;
    }
    if d.is_some() && (*p < B.p) && (*p >= 6) {
        return 1;
    }

    /*    return success */
    0
}

/// The shared method/embedding condition sequence of C's
/// ARKodeButcherTable_CheckOrder (literal duplicates in C differing
/// only in the label and weight vector). Returns the measured order,
/// including the simplifying-assumptions extension past order 6.
fn arkode_butcher_check_all_conditions(
    A: &[Vec<f64>],
    b: &[f64],
    c: &[f64],
    s: i32,
    outfile: &mut Option<&mut dyn std::io::Write>,
    label: &str,
) -> i32 {
    let mut q: i32;

    /*    row sum condition */
    if arkode_butcher_rowsum(A, c, s) {
        q = 0;
    } else {
        q = -1;
        outmsg(outfile, &format!("  {} fails row sum condition\n", label));
    }
    /*    order 1 condition */
    if q == 0 {
        if arkode_butcher_order1(b, s) {
            q = 1;
        } else {
            outmsg(outfile, &format!("  {} fails order 1 condition\n", label));
        }
    }
    /*    order 2 condition */
    if q == 1 {
        if arkode_butcher_order2(b, c, s) {
            q = 2;
        } else {
            outmsg(outfile, &format!("  {} fails order 2 condition\n", label));
        }
    }
    /*    order 3 conditions */
    if q == 2 {
        let mut alltrue = true;
        if !arkode_butcher_order3a(b, c, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 3 condition A\n", label));
        }
        if !arkode_butcher_order3b(b, A, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 3 condition B\n", label));
        }
        if alltrue {
            q = 3;
        }
    }
    /*    order 4 conditions */
    if q == 3 {
        let mut alltrue = true;
        if !arkode_butcher_order4a(b, c, c, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 4 condition A\n", label));
        }
        if !arkode_butcher_order4b(b, c, A, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 4 condition B\n", label));
        }
        if !arkode_butcher_order4c(b, A, c, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 4 condition C\n", label));
        }
        if !arkode_butcher_order4d(b, A, A, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 4 condition D\n", label));
        }
        if alltrue {
            q = 4;
        }
    }
    /*    order 5 conditions */
    if q == 4 {
        let mut alltrue = true;
        if !arkode_butcher_order5a(b, c, c, c, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 5 condition A\n", label));
        }
        if !arkode_butcher_order5b(b, c, c, A, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 5 condition B\n", label));
        }
        if !arkode_butcher_order5c(b, A, c, A, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 5 condition C\n", label));
        }
        if !arkode_butcher_order5d(b, c, A, c, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 5 condition D\n", label));
        }
        if !arkode_butcher_order5e(b, A, c, c, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 5 condition E\n", label));
        }
        if !arkode_butcher_order5f(b, c, A, A, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 5 condition F\n", label));
        }
        if !arkode_butcher_order5g(b, A, c, A, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 5 condition G\n", label));
        }
        if !arkode_butcher_order5h(b, A, A, c, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 5 condition H\n", label));
        }
        if !arkode_butcher_order5i(b, A, A, A, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 5 condition I\n", label));
        }
        if alltrue {
            q = 5;
        }
    }
    /*    order 6 conditions */
    if q == 5 {
        let mut alltrue = true;
        if !arkode_butcher_order6a(b, c, c, c, c, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 6 condition A\n", label));
        }
        if !arkode_butcher_order6b(b, c, c, c, A, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 6 condition B\n", label));
        }
        if !arkode_butcher_order6c(b, c, A, c, A, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 6 condition C\n", label));
        }
        if !arkode_butcher_order6d(b, c, c, A, c, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 6 condition D\n", label));
        }
        if !arkode_butcher_order6e(b, c, c, A, A, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 6 condition E\n", label));
        }
        if !arkode_butcher_order6f(b, A, A, c, A, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 6 condition F\n", label));
        }
        if !arkode_butcher_order6g(b, c, A, c, c, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 6 condition G\n", label));
        }
        if !arkode_butcher_order6h(b, c, A, c, A, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 6 condition H\n", label));
        }
        if !arkode_butcher_order6i(b, c, A, A, c, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 6 condition I\n", label));
        }
        if !arkode_butcher_order6j(b, c, A, A, A, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 6 condition J\n", label));
        }
        if !arkode_butcher_order6k(b, A, c, c, c, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 6 condition K\n", label));
        }
        if !arkode_butcher_order6l(b, A, c, c, A, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 6 condition L\n", label));
        }
        if !arkode_butcher_order6m(b, A, A, c, A, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 6 condition M\n", label));
        }
        if !arkode_butcher_order6n(b, A, c, A, c, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 6 condition N\n", label));
        }
        if !arkode_butcher_order6o(b, A, c, A, A, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 6 condition O\n", label));
        }
        if !arkode_butcher_order6p(b, A, A, c, c, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 6 condition P\n", label));
        }
        if !arkode_butcher_order6q(b, A, A, c, A, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 6 condition Q\n", label));
        }
        if !arkode_butcher_order6r(b, A, A, A, c, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 6 condition R\n", label));
        }
        if !arkode_butcher_order6s(b, A, A, A, A, c, s) {
            alltrue = false;
            outmsg(outfile, &format!("  {} fails order 6 condition S\n", label));
        }
        if alltrue {
            q = 6;
        }
    }
    /*    higher order conditions (via simplifying assumptions) */
    if q == 6 {
        outmsg(
            outfile,
            &format!("  {} order >= 6; reverting to simplifying assumptions\n", label),
        );
        let q_sa = __ButcherSimplifyingAssumptions(A, b, c, s);
        q = q.max(q_sa);
        outmsg(outfile, &format!("  {} order = {}\n", label, q));
    }

    q
}

/*---------------------------------------------------------------
  Routine to determine the analytical order of accuracy for a
  specified pair of Butcher tables in an ARK pair.  We check the
  analytical order conditions up through order 6.

  Return values:
     0 (success): completed checks
     1 (warning): internal {q,p} values are lower than analytical
        order, or method achieves maximum order possible with this
        routine and internal {q,p} are higher.
    -1 (failure): NULL-valued B1, B2 (or critical contents)
  ---------------------------------------------------------------*/
pub fn ARKodeButcherTable_CheckARKOrder(
    B1: &ARKodeButcherTable,
    B2: &ARKodeButcherTable,
    q: &mut i32,
    p: &mut i32,
    outfile: Option<&mut dyn std::io::Write>,
) -> i32 {
    let mut outfile = outfile;
    *q = 0;
    *p = 0;

    /* verify Butcher table contents (C NULL checks unrepresentable) */
    if B1.stages < 1 {
        return -1;
    }
    if B2.stages < 1 {
        return -1;
    }
    if B1.stages != B2.stages {
        return -1;
    }

    /* set shortcuts for Butcher table components */
    let A = [&B1.A, &B2.A];
    let b: [&[f64]; 2] = [&B1.b, &B2.b];
    let c: [&[f64]; 2] = [&B1.c, &B2.c];
    /* C quirk preserved: d[1] = B1->d, not B2->d (arkode_butcher.c:1168) */
    let d: [Option<&Vec<f64>>; 2] = [B1.d.as_ref(), B1.d.as_ref()];
    let s = B1.stages;

    /* check method order */
    outmsg(&mut outfile, "ARKodeButcherTable_CheckARKOrder:\n");
    *q = arkode_butcher_check_all_ark_conditions(A, b, c, s, &mut outfile, "method");

    /* check embedding order */
    let embedded = d[0].is_some() && d[1].is_some();
    if embedded {
        outmsg(&mut outfile, "\n");
        let w: [&[f64]; 2] = [d[0].unwrap(), d[1].unwrap()];
        *p = arkode_butcher_check_all_ark_conditions(A, w, c, s, &mut outfile, "embedding");
    }

    /* compare results against stored values and return */

    /*    check warning modes */
    if *q > B1.q {
        return 1;
    }
    if *q > B2.q {
        return 1;
    }
    if embedded {
        if *p > B1.p {
            return 1;
        }
        if *p > B2.p {
            return 1;
        }
    }
    if (*q < B1.q) && (*q == 6) {
        return 1;
    }
    if (*q < B2.q) && (*q == 6) {
        return 1;
    }
    if embedded {
        if (*p < B1.p) && (*p == 6) {
            return 1;
        }
        if (*p < B2.p) && (*p == 6) {
            return 1;
        }
    }

    /*    return success */
    0
}

/// The shared method/embedding condition sequence of C's
/// ARKodeButcherTable_CheckARKOrder: every condition is checked for
/// all combinations of the two tables' coefficients. `alltrue`
/// accumulates across the condition families of each order block
/// exactly as in C (a failed family keeps later families' messages
/// printing). Returns the measured order (no simplifying-assumptions
/// step in the ARK variant).
fn arkode_butcher_check_all_ark_conditions(
    A: [&Vec<Vec<f64>>; 2],
    w: [&[f64]; 2],
    c: [&[f64]; 2],
    s: i32,
    outfile: &mut Option<&mut dyn std::io::Write>,
    label: &str,
) -> i32 {
    let mut q: i32;

    /*    row sum conditions */
    if arkode_butcher_rowsum(A[0], c[0], s) && arkode_butcher_rowsum(A[1], c[1], s) {
        q = 0;
    } else {
        q = -1;
        outmsg(outfile, &format!("  {} fails row sum conditions\n", label));
    }
    /*    order 1 conditions */
    if q == 0 {
        if arkode_butcher_order1(w[0], s) && arkode_butcher_order1(w[1], s) {
            q = 1;
        } else {
            outmsg(outfile, &format!("  {} fails order 1 conditions\n", label));
        }
    }
    /*    order 2 conditions */
    if q == 1 {
        let mut alltrue = true;
        for i in 0..2 {
            for j in 0..2 {
                alltrue = alltrue && arkode_butcher_order2(w[i], c[j], s);
            }
        }
        if alltrue {
            q = 2;
        } else {
            outmsg(outfile, &format!("  {} fails order 2 conditions\n", label));
        }
    }
    /*    order 3 conditions */
    if q == 2 {
        let mut alltrue = true;
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    alltrue = alltrue && arkode_butcher_order3a(w[i], c[j], c[k], s);
                }
            }
        }
        if !alltrue {
            outmsg(outfile, &format!("  {} fails order 3 conditions A\n", label));
        }
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    alltrue = alltrue && arkode_butcher_order3b(w[i], A[j], c[k], s);
                }
            }
        }
        if !alltrue {
            outmsg(outfile, &format!("  {} fails order 3 conditions B\n", label));
        }
        if alltrue {
            q = 3;
        }
    }
    /*    order 4 conditions */
    if q == 3 {
        let mut alltrue = true;
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    for l in 0..2 {
                        alltrue = alltrue && arkode_butcher_order4a(w[i], c[j], c[k], c[l], s);
                    }
                }
            }
        }
        if !alltrue {
            outmsg(outfile, &format!("  {} fails order 4 conditions A\n", label));
        }
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    for l in 0..2 {
                        alltrue = alltrue && arkode_butcher_order4b(w[i], c[j], A[k], c[l], s);
                    }
                }
            }
        }
        if !alltrue {
            outmsg(outfile, &format!("  {} fails order 4 conditions B\n", label));
        }
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    for l in 0..2 {
                        alltrue = alltrue && arkode_butcher_order4c(w[i], A[j], c[k], c[l], s);
                    }
                }
            }
        }
        if !alltrue {
            outmsg(outfile, &format!("  {} fails order 4 conditions C\n", label));
        }
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    for l in 0..2 {
                        alltrue = alltrue && arkode_butcher_order4d(w[i], A[j], A[k], c[l], s);
                    }
                }
            }
        }
        if !alltrue {
            outmsg(outfile, &format!("  {} fails order 4 conditions D\n", label));
        }
        if alltrue {
            q = 4;
        }
    }
    /*    order 5 conditions */
    if q == 4 {
        let mut alltrue = true;
        macro_rules! ark5 {
            ($f:ident, $s1:ident, $s2:ident, $s3:ident, $letter:expr) => {
                for i in 0..2 {
                    for j in 0..2 {
                        for k in 0..2 {
                            for l in 0..2 {
                                for m in 0..2 {
                                    alltrue = alltrue
                                        && $f(w[i], $s1[j], $s2[k], $s3[l], c[m], s);
                                }
                            }
                        }
                    }
                }
                if !alltrue {
                    outmsg(
                        outfile,
                        &format!("  {} fails order 5 conditions {}\n", label, $letter),
                    );
                }
            };
        }
        ark5!(arkode_butcher_order5a, c, c, c, "A");
        ark5!(arkode_butcher_order5b, c, c, A, "B");
        ark5!(arkode_butcher_order5c, A, c, A, "C");
        ark5!(arkode_butcher_order5d, c, A, c, "D");
        ark5!(arkode_butcher_order5e, A, c, c, "E");
        ark5!(arkode_butcher_order5f, c, A, A, "F");
        ark5!(arkode_butcher_order5g, A, c, A, "G");
        ark5!(arkode_butcher_order5h, A, A, c, "H");
        ark5!(arkode_butcher_order5i, A, A, A, "I");
        if alltrue {
            q = 5;
        }
    }
    /*    order 6 conditions */
    if q == 5 {
        let mut alltrue = true;
        macro_rules! ark6 {
            ($f:ident, $s1:ident, $s2:ident, $s3:ident, $s4:ident, $letter:expr) => {
                for i in 0..2 {
                    for j in 0..2 {
                        for k in 0..2 {
                            for l in 0..2 {
                                for m in 0..2 {
                                    for n in 0..2 {
                                        alltrue = alltrue
                                            && $f(w[i], $s1[j], $s2[k], $s3[l], $s4[m], c[n], s);
                                    }
                                }
                            }
                        }
                    }
                }
                if !alltrue {
                    outmsg(
                        outfile,
                        &format!("  {} fails order 6 conditions {}\n", label, $letter),
                    );
                }
            };
        }
        ark6!(arkode_butcher_order6a, c, c, c, c, "A");
        ark6!(arkode_butcher_order6b, c, c, c, A, "B");
        ark6!(arkode_butcher_order6c, c, A, c, A, "C");
        ark6!(arkode_butcher_order6d, c, c, A, c, "D");
        ark6!(arkode_butcher_order6e, c, c, A, A, "E");
        ark6!(arkode_butcher_order6f, A, A, c, A, "F");
        ark6!(arkode_butcher_order6g, c, A, c, c, "G");
        ark6!(arkode_butcher_order6h, c, A, c, A, "H");
        ark6!(arkode_butcher_order6i, c, A, A, c, "I");
        ark6!(arkode_butcher_order6j, c, A, A, A, "J");
        ark6!(arkode_butcher_order6k, A, c, c, c, "K");
        ark6!(arkode_butcher_order6l, A, c, c, A, "L");
        ark6!(arkode_butcher_order6m, A, A, c, A, "M");
        ark6!(arkode_butcher_order6n, A, c, A, c, "N");
        ark6!(arkode_butcher_order6o, A, c, A, A, "O");
        ark6!(arkode_butcher_order6p, A, A, c, c, "P");
        ark6!(arkode_butcher_order6q, A, A, c, A, "Q");
        ark6!(arkode_butcher_order6r, A, A, A, c, "R");
        ark6!(arkode_butcher_order6s, A, A, A, A, "S");
        if alltrue {
            q = 6;
        }
    }

    q
}

/*---------------------------------------------------------------
  Private utility routines for checking method order
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  Utility routine to compute small dense matrix-vector product
       b = A*x
  Here A is (s x s), x and b are (s x 1).  Returns 0 on success,
  nonzero on failure.
  ---------------------------------------------------------------*/
fn arkode_butcher_mv(A: &[Vec<f64>], x: &[f64], s: i32, b: &mut [f64]) -> i32 {
    if s < 1 {
        return 1;
    }
    let s = s as usize;
    for i in 0..s {
        b[i] = 0.0;
    }
    for i in 0..s {
        for j in 0..s {
            b[i] += A[i][j] * x[j];
        }
    }
    0
}

/*---------------------------------------------------------------
  Utility routine to compute small vector .* vector product
       z = x.*y   [Matlab notation]
  ---------------------------------------------------------------*/
fn arkode_butcher_vv(x: &[f64], y: &[f64], s: i32, z: &mut [f64]) -> i32 {
    if s < 1 {
        return 1;
    }
    for i in 0..s as usize {
        z[i] = x[i] * y[i];
    }
    0
}

/*---------------------------------------------------------------
  Utility routine to compute small vector .^ int
       z = x.^l   [Matlab notation]
  ---------------------------------------------------------------*/
fn arkode_butcher_vp(x: &[f64], l: i32, s: i32, z: &mut [f64]) -> i32 {
    if s < 1 {
        return 1;
    }
    for i in 0..s as usize {
        z[i] = SUNRpowerI(x[i], l);
    }
    0
}

/*---------------------------------------------------------------
  Utility routine to compute small vector dot product:
       d = dot(x,y)
  ---------------------------------------------------------------*/
fn arkode_butcher_dot(x: &[f64], y: &[f64], s: i32, d: &mut f64) -> i32 {
    if s < 1 {
        return 1;
    }
    *d = 0.0;
    for i in 0..s as usize {
        *d += x[i] * y[i];
    }
    0
}

/*---------------------------------------------------------------
  Utility routines to check specific order conditions.  Each
  returns SUNTRUE on success, SUNFALSE on failure.
     Order 0:  arkode_butcher_rowsum
     Order 1:  arkode_butcher_order1
     Order 2:  arkode_butcher_order2
     Order 3:  arkode_butcher_order3a and arkode_butcher_order3b
     Order 4:  arkode_butcher_order4a through arkode_butcher_order4d
     Order 5:  arkode_butcher_order5a through arkode_butcher_order5i
     Order 6:  arkode_butcher_order6a through arkode_butcher_order6s
  (comparisons are written !(err > TOL) so NaN behaves as in C)
  ---------------------------------------------------------------*/

/* c(i) = sum(A(i,:)) */
fn arkode_butcher_rowsum(A: &[Vec<f64>], c: &[f64], s: i32) -> bool {
    for i in 0..s.max(0) as usize {
        let mut rsum = 0.0;
        for j in 0..s as usize {
            rsum += A[i][j];
        }
        if SUNRabs(rsum - c[i]) > TOL() {
            return false;
        }
    }
    true
}

/* b'*e = 1 */
fn arkode_butcher_order1(b: &[f64], s: i32) -> bool {
    let mut err = 1.0;
    for i in 0..s.max(0) as usize {
        err -= b[i];
    }
    !(SUNRabs(err) > TOL())
}

/* b'*c = 1/2 */
fn arkode_butcher_order2(b: &[f64], c: &[f64], s: i32) -> bool {
    let mut bc = 0.0;
    if arkode_butcher_dot(b, c, s, &mut bc) != 0 {
        return false;
    }
    !(SUNRabs(bc - 0.5) > TOL())
}

/* b'*(c1.*c2) = 1/3 */
fn arkode_butcher_order3a(b: &[f64], c1: &[f64], c2: &[f64], s: i32) -> bool {
    let n = s.max(0) as usize;
    let mut tmp = vec![0.0; n];
    let mut bcc = 0.0;
    if arkode_butcher_vv(c1, c2, s, &mut tmp) != 0 {
        return false;
    }
    if arkode_butcher_dot(b, &tmp, s, &mut bcc) != 0 {
        return false;
    }
    !(SUNRabs(bcc - 1.0 / 3.0) > TOL())
}

/* b'*(A*c) = 1/6 */
fn arkode_butcher_order3b(b: &[f64], A: &[Vec<f64>], c: &[f64], s: i32) -> bool {
    let n = s.max(0) as usize;
    let mut tmp = vec![0.0; n];
    let mut bAc = 0.0;
    if arkode_butcher_mv(A, c, s, &mut tmp) != 0 {
        return false;
    }
    if arkode_butcher_dot(b, &tmp, s, &mut bAc) != 0 {
        return false;
    }
    !(SUNRabs(bAc - 1.0 / 6.0) > TOL())
}

/* b'*(c1.*c2.*c3) = 1/4 */
fn arkode_butcher_order4a(b: &[f64], c1: &[f64], c2: &[f64], c3: &[f64], s: i32) -> bool {
    let n = s.max(0) as usize;
    let mut tmp1 = vec![0.0; n];
    let mut tmp2 = vec![0.0; n];
    let mut bccc = 0.0;
    if arkode_butcher_vv(c1, c2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_vv(c3, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_dot(b, &tmp2, s, &mut bccc) != 0 {
        return false;
    }
    !(SUNRabs(bccc - 0.25) > TOL())
}

/* (b.*c1)'*(A*c2) = 1/8 */
fn arkode_butcher_order4b(b: &[f64], c1: &[f64], A: &[Vec<f64>], c2: &[f64], s: i32) -> bool {
    let n = s.max(0) as usize;
    let mut tmp1 = vec![0.0; n];
    let mut tmp2 = vec![0.0; n];
    let mut bcAc = 0.0;
    if arkode_butcher_vv(b, c1, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_mv(A, c2, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_dot(&tmp1, &tmp2, s, &mut bcAc) != 0 {
        return false;
    }
    !(SUNRabs(bcAc - 0.125) > TOL())
}

/* b'*A*(c1.*c2) = 1/12 */
fn arkode_butcher_order4c(b: &[f64], A: &[Vec<f64>], c1: &[f64], c2: &[f64], s: i32) -> bool {
    let n = s.max(0) as usize;
    let mut tmp1 = vec![0.0; n];
    let mut tmp2 = vec![0.0; n];
    let mut bAcc = 0.0;
    if arkode_butcher_vv(c1, c2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_mv(A, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_dot(b, &tmp2, s, &mut bAcc) != 0 {
        return false;
    }
    !(SUNRabs(bAcc - 1.0 / 12.0) > TOL())
}

/* b'*A1*A2*c = 1/24 */
fn arkode_butcher_order4d(
    b: &[f64],
    A1: &[Vec<f64>],
    A2: &[Vec<f64>],
    c: &[f64],
    s: i32,
) -> bool {
    let n = s.max(0) as usize;
    let mut tmp1 = vec![0.0; n];
    let mut tmp2 = vec![0.0; n];
    let mut bAAc = 0.0;
    if arkode_butcher_mv(A2, c, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_mv(A1, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_dot(b, &tmp2, s, &mut bAAc) != 0 {
        return false;
    }
    !(SUNRabs(bAAc - 1.0 / 24.0) > TOL())
}

/* b'*(c1.*c2.*c3.*c4) = 1/5 */
fn arkode_butcher_order5a(
    b: &[f64],
    c1: &[f64],
    c2: &[f64],
    c3: &[f64],
    c4: &[f64],
    s: i32,
) -> bool {
    let n = s.max(0) as usize;
    let mut tmp1 = vec![0.0; n];
    let mut tmp2 = vec![0.0; n];
    let mut bcccc = 0.0;
    if arkode_butcher_vv(c1, c2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_vv(c3, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_vv(c4, &tmp2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_dot(b, &tmp1, s, &mut bcccc) != 0 {
        return false;
    }
    !(SUNRabs(bcccc - 0.2) > TOL())
}

/* (b.*c1.*c2)'*(A*c3) = 1/10 */
fn arkode_butcher_order5b(
    b: &[f64],
    c1: &[f64],
    c2: &[f64],
    A: &[Vec<f64>],
    c3: &[f64],
    s: i32,
) -> bool {
    let n = s.max(0) as usize;
    let mut tmp1 = vec![0.0; n];
    let mut tmp2 = vec![0.0; n];
    let mut bccAc = 0.0;
    if arkode_butcher_vv(c1, c2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_vv(b, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_mv(A, c3, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_dot(&tmp1, &tmp2, s, &mut bccAc) != 0 {
        return false;
    }
    !(SUNRabs(bccAc - 0.1) > TOL())
}

/* b'*((A1*c1).*(A2*c2)) = 1/20 */
fn arkode_butcher_order5c(
    b: &[f64],
    A1: &[Vec<f64>],
    c1: &[f64],
    A2: &[Vec<f64>],
    c2: &[f64],
    s: i32,
) -> bool {
    let n = s.max(0) as usize;
    let mut tmp1 = vec![0.0; n];
    let mut tmp2 = vec![0.0; n];
    let mut tmp3 = vec![0.0; n];
    let mut bAcAc = 0.0;
    if arkode_butcher_mv(A1, c1, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_mv(A2, c2, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_vv(&tmp1, &tmp2, s, &mut tmp3) != 0 {
        return false;
    }
    if arkode_butcher_dot(b, &tmp3, s, &mut bAcAc) != 0 {
        return false;
    }
    !(SUNRabs(bAcAc - 0.05) > TOL())
}

/* (b.*c1)'*A*(c2.*c3) = 1/15 */
fn arkode_butcher_order5d(
    b: &[f64],
    c1: &[f64],
    A: &[Vec<f64>],
    c2: &[f64],
    c3: &[f64],
    s: i32,
) -> bool {
    let n = s.max(0) as usize;
    let mut tmp1 = vec![0.0; n];
    let mut tmp2 = vec![0.0; n];
    let mut bcAcc = 0.0;
    if arkode_butcher_vv(c2, c3, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_mv(A, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_vv(b, c1, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_dot(&tmp1, &tmp2, s, &mut bcAcc) != 0 {
        return false;
    }
    !(SUNRabs(bcAcc - 1.0 / 15.0) > TOL())
}

/* b'*A*(c1.*c2.*c3) = 1/20 */
fn arkode_butcher_order5e(
    b: &[f64],
    A: &[Vec<f64>],
    c1: &[f64],
    c2: &[f64],
    c3: &[f64],
    s: i32,
) -> bool {
    let n = s.max(0) as usize;
    let mut tmp1 = vec![0.0; n];
    let mut tmp2 = vec![0.0; n];
    let mut bAccc = 0.0;
    if arkode_butcher_vv(c1, c2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_vv(c3, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_mv(A, &tmp2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_dot(b, &tmp1, s, &mut bAccc) != 0 {
        return false;
    }
    !(SUNRabs(bAccc - 0.05) > TOL())
}

/* (b.*c1)'*A1*A2*c2 = 1/30 */
fn arkode_butcher_order5f(
    b: &[f64],
    c1: &[f64],
    A1: &[Vec<f64>],
    A2: &[Vec<f64>],
    c2: &[f64],
    s: i32,
) -> bool {
    let n = s.max(0) as usize;
    let mut tmp1 = vec![0.0; n];
    let mut tmp2 = vec![0.0; n];
    let mut bcAAc = 0.0;
    if arkode_butcher_mv(A2, c2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_mv(A1, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_vv(b, c1, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_dot(&tmp1, &tmp2, s, &mut bcAAc) != 0 {
        return false;
    }
    !(SUNRabs(bcAAc - 1.0 / 30.0) > TOL())
}

/* b'*A1*(c1.*(A2*c2)) = 1/40 */
fn arkode_butcher_order5g(
    b: &[f64],
    A1: &[Vec<f64>],
    c1: &[f64],
    A2: &[Vec<f64>],
    c2: &[f64],
    s: i32,
) -> bool {
    let n = s.max(0) as usize;
    let mut tmp1 = vec![0.0; n];
    let mut tmp2 = vec![0.0; n];
    let mut bAcAc = 0.0;
    if arkode_butcher_mv(A2, c2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_vv(c1, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_mv(A1, &tmp2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_dot(b, &tmp1, s, &mut bAcAc) != 0 {
        return false;
    }
    !(SUNRabs(bAcAc - 1.0 / 40.0) > TOL())
}

/* b'*A1*A2*(c1.*c2) = 1/60 */
fn arkode_butcher_order5h(
    b: &[f64],
    A1: &[Vec<f64>],
    A2: &[Vec<f64>],
    c1: &[f64],
    c2: &[f64],
    s: i32,
) -> bool {
    let n = s.max(0) as usize;
    let mut tmp1 = vec![0.0; n];
    let mut tmp2 = vec![0.0; n];
    let mut bAAcc = 0.0;
    if arkode_butcher_vv(c1, c2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_mv(A2, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_mv(A1, &tmp2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_dot(b, &tmp1, s, &mut bAAcc) != 0 {
        return false;
    }
    !(SUNRabs(bAAcc - 1.0 / 60.0) > TOL())
}

/* b'*A1*A2*A3*c = 1/120 */
fn arkode_butcher_order5i(
    b: &[f64],
    A1: &[Vec<f64>],
    A2: &[Vec<f64>],
    A3: &[Vec<f64>],
    c: &[f64],
    s: i32,
) -> bool {
    let n = s.max(0) as usize;
    let mut tmp1 = vec![0.0; n];
    let mut tmp2 = vec![0.0; n];
    let mut bAAAc = 0.0;
    if arkode_butcher_mv(A3, c, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_mv(A2, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_mv(A1, &tmp2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_dot(b, &tmp1, s, &mut bAAAc) != 0 {
        return false;
    }
    !(SUNRabs(bAAAc - 1.0 / 120.0) > TOL())
}

/* The 19 order-6 conditions share one shape: a chain of vv/mv steps
   into two (or three) temporaries, a final dot, and a target. Each
   function below transcribes its C chain step-for-step. */

/* b'*(c1.*c2.*c3.*c4.*c5) = 1/6 */
#[allow(clippy::too_many_arguments)]
fn arkode_butcher_order6a(
    b: &[f64],
    c1: &[f64],
    c2: &[f64],
    c3: &[f64],
    c4: &[f64],
    c5: &[f64],
    s: i32,
) -> bool {
    let n = s.max(0) as usize;
    let (mut tmp1, mut tmp2) = (vec![0.0; n], vec![0.0; n]);
    let mut bccccc = 0.0;
    if arkode_butcher_vv(c1, c2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_vv(c3, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_vv(c4, &tmp2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_vv(c5, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_dot(b, &tmp2, s, &mut bccccc) != 0 {
        return false;
    }
    !(SUNRabs(bccccc - 1.0 / 6.0) > TOL())
}

/* (b.*c1.*c2.*c3)'*(A*c4) = 1/12 */
#[allow(clippy::too_many_arguments)]
fn arkode_butcher_order6b(
    b: &[f64],
    c1: &[f64],
    c2: &[f64],
    c3: &[f64],
    A: &[Vec<f64>],
    c4: &[f64],
    s: i32,
) -> bool {
    let n = s.max(0) as usize;
    let (mut tmp1, mut tmp2) = (vec![0.0; n], vec![0.0; n]);
    let mut bcccAc = 0.0;
    if arkode_butcher_vv(b, c1, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_vv(c2, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_vv(c3, &tmp2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_mv(A, c4, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_dot(&tmp1, &tmp2, s, &mut bcccAc) != 0 {
        return false;
    }
    !(SUNRabs(bcccAc - 1.0 / 12.0) > TOL())
}

/* b'*(c1.*(A1*c2).*(A2*c3)) = 1/24 */
#[allow(clippy::too_many_arguments)]
fn arkode_butcher_order6c(
    b: &[f64],
    c1: &[f64],
    A1: &[Vec<f64>],
    c2: &[f64],
    A2: &[Vec<f64>],
    c3: &[f64],
    s: i32,
) -> bool {
    let n = s.max(0) as usize;
    let (mut tmp1, mut tmp2, mut tmp3) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    let mut bcAc2 = 0.0;
    if arkode_butcher_mv(A2, c3, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_mv(A1, c2, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_vv(&tmp1, &tmp2, s, &mut tmp3) != 0 {
        return false;
    }
    if arkode_butcher_vv(c1, &tmp3, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_dot(b, &tmp1, s, &mut bcAc2) != 0 {
        return false;
    }
    !(SUNRabs(bcAc2 - 1.0 / 24.0) > TOL())
}

/* (b.*c1.*c2)'*A*(c3.*c4) = 1/18 */
#[allow(clippy::too_many_arguments)]
fn arkode_butcher_order6d(
    b: &[f64],
    c1: &[f64],
    c2: &[f64],
    A: &[Vec<f64>],
    c3: &[f64],
    c4: &[f64],
    s: i32,
) -> bool {
    let n = s.max(0) as usize;
    let (mut tmp1, mut tmp2, mut tmp3) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    let mut bccAcc = 0.0;
    if arkode_butcher_vv(c3, c4, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_mv(A, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_vv(c1, c2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_vv(b, &tmp1, s, &mut tmp3) != 0 {
        return false;
    }
    if arkode_butcher_dot(&tmp2, &tmp3, s, &mut bccAcc) != 0 {
        return false;
    }
    !(SUNRabs(bccAcc - 1.0 / 18.0) > TOL())
}

/* (b.*(c1.*c2))'*A1*A2*c3 = 1/36 */
#[allow(clippy::too_many_arguments)]
fn arkode_butcher_order6e(
    b: &[f64],
    c1: &[f64],
    c2: &[f64],
    A1: &[Vec<f64>],
    A2: &[Vec<f64>],
    c3: &[f64],
    s: i32,
) -> bool {
    let n = s.max(0) as usize;
    let (mut tmp1, mut tmp2, mut tmp3) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    let mut bccAAc = 0.0;
    if arkode_butcher_vv(c1, c2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_vv(b, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_mv(A2, c3, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_mv(A1, &tmp1, s, &mut tmp3) != 0 {
        return false;
    }
    if arkode_butcher_dot(&tmp2, &tmp3, s, &mut bccAAc) != 0 {
        return false;
    }
    !(SUNRabs(bccAAc - 1.0 / 36.0) > TOL())
}

/* b'*((A1*A2*c1).*(A3*c2)) = 1/72 */
#[allow(clippy::too_many_arguments)]
fn arkode_butcher_order6f(
    b: &[f64],
    A1: &[Vec<f64>],
    A2: &[Vec<f64>],
    c1: &[f64],
    A3: &[Vec<f64>],
    c2: &[f64],
    s: i32,
) -> bool {
    let n = s.max(0) as usize;
    let (mut tmp1, mut tmp2, mut tmp3) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    let mut bAAcAc = 0.0;
    if arkode_butcher_mv(A2, c1, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_mv(A1, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_mv(A3, c2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_vv(&tmp1, &tmp2, s, &mut tmp3) != 0 {
        return false;
    }
    if arkode_butcher_dot(b, &tmp3, s, &mut bAAcAc) != 0 {
        return false;
    }
    !(SUNRabs(bAAcAc - 1.0 / 72.0) > TOL())
}

/* b'*(c1.*(A*(c2.*c3.*c4))) = 1/24 */
#[allow(clippy::too_many_arguments)]
fn arkode_butcher_order6g(
    b: &[f64],
    c1: &[f64],
    A: &[Vec<f64>],
    c2: &[f64],
    c3: &[f64],
    c4: &[f64],
    s: i32,
) -> bool {
    let n = s.max(0) as usize;
    let (mut tmp1, mut tmp2) = (vec![0.0; n], vec![0.0; n]);
    let mut bcAccc = 0.0;
    if arkode_butcher_vv(c2, c3, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_vv(c4, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_mv(A, &tmp2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_vv(c1, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_dot(b, &tmp2, s, &mut bcAccc) != 0 {
        return false;
    }
    !(SUNRabs(bcAccc - 1.0 / 24.0) > TOL())
}

/* b'*(c1.*(A1*(c2.*(A2*c3)))) = 1/48 */
#[allow(clippy::too_many_arguments)]
fn arkode_butcher_order6h(
    b: &[f64],
    c1: &[f64],
    A1: &[Vec<f64>],
    c2: &[f64],
    A2: &[Vec<f64>],
    c3: &[f64],
    s: i32,
) -> bool {
    let n = s.max(0) as usize;
    let (mut tmp1, mut tmp2) = (vec![0.0; n], vec![0.0; n]);
    let mut bcAcAc = 0.0;
    if arkode_butcher_mv(A2, c3, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_vv(c2, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_mv(A1, &tmp2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_vv(c1, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_dot(b, &tmp2, s, &mut bcAcAc) != 0 {
        return false;
    }
    !(SUNRabs(bcAcAc - 1.0 / 48.0) > TOL())
}

/* b'*(c1.*(A1*A2*(c2.*c3))) = 1/72 */
#[allow(clippy::too_many_arguments)]
fn arkode_butcher_order6i(
    b: &[f64],
    c1: &[f64],
    A1: &[Vec<f64>],
    A2: &[Vec<f64>],
    c2: &[f64],
    c3: &[f64],
    s: i32,
) -> bool {
    let n = s.max(0) as usize;
    let (mut tmp1, mut tmp2) = (vec![0.0; n], vec![0.0; n]);
    let mut bcAAcc = 0.0;
    if arkode_butcher_vv(c2, c3, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_mv(A2, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_mv(A1, &tmp2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_vv(c1, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_dot(b, &tmp2, s, &mut bcAAcc) != 0 {
        return false;
    }
    !(SUNRabs(bcAAcc - 1.0 / 72.0) > TOL())
}

/* b'*(c1.*(A1*A2*A3*c2)) = 1/144 */
#[allow(clippy::too_many_arguments)]
fn arkode_butcher_order6j(
    b: &[f64],
    c1: &[f64],
    A1: &[Vec<f64>],
    A2: &[Vec<f64>],
    A3: &[Vec<f64>],
    c2: &[f64],
    s: i32,
) -> bool {
    let n = s.max(0) as usize;
    let (mut tmp1, mut tmp2) = (vec![0.0; n], vec![0.0; n]);
    let mut bcAAAc = 0.0;
    if arkode_butcher_mv(A3, c2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_mv(A2, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_mv(A1, &tmp2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_vv(c1, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_dot(b, &tmp2, s, &mut bcAAAc) != 0 {
        return false;
    }
    !(SUNRabs(bcAAAc - 1.0 / 144.0) > TOL())
}

/* b'*A*(c1.*c2.*c3.*c4) = 1/30 */
#[allow(clippy::too_many_arguments)]
fn arkode_butcher_order6k(
    b: &[f64],
    A: &[Vec<f64>],
    c1: &[f64],
    c2: &[f64],
    c3: &[f64],
    c4: &[f64],
    s: i32,
) -> bool {
    let n = s.max(0) as usize;
    let (mut tmp1, mut tmp2) = (vec![0.0; n], vec![0.0; n]);
    let mut bAcccc = 0.0;
    if arkode_butcher_vv(c1, c2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_vv(c3, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_vv(c4, &tmp2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_mv(A, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_dot(b, &tmp2, s, &mut bAcccc) != 0 {
        return false;
    }
    !(SUNRabs(bAcccc - 1.0 / 30.0) > TOL())
}

/* b'*A1*(c1.*c2.*(A2*c3)) = 1/60 */
#[allow(clippy::too_many_arguments)]
fn arkode_butcher_order6l(
    b: &[f64],
    A1: &[Vec<f64>],
    c1: &[f64],
    c2: &[f64],
    A2: &[Vec<f64>],
    c3: &[f64],
    s: i32,
) -> bool {
    let n = s.max(0) as usize;
    let (mut tmp1, mut tmp2) = (vec![0.0; n], vec![0.0; n]);
    let mut bAccAc = 0.0;
    if arkode_butcher_mv(A2, c3, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_vv(c2, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_vv(c1, &tmp2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_mv(A1, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_dot(b, &tmp2, s, &mut bAccAc) != 0 {
        return false;
    }
    !(SUNRabs(bAccAc - 1.0 / 60.0) > TOL())
}

/* b'*A1*((A2*c1).*(A3*c2)) = 1/120 */
#[allow(clippy::too_many_arguments)]
fn arkode_butcher_order6m(
    b: &[f64],
    A1: &[Vec<f64>],
    A2: &[Vec<f64>],
    c1: &[f64],
    A3: &[Vec<f64>],
    c2: &[f64],
    s: i32,
) -> bool {
    let n = s.max(0) as usize;
    let (mut tmp1, mut tmp2, mut tmp3) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    let mut bAAcAc = 0.0;
    if arkode_butcher_mv(A3, c2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_mv(A2, c1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_vv(&tmp1, &tmp2, s, &mut tmp3) != 0 {
        return false;
    }
    if arkode_butcher_mv(A1, &tmp3, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_dot(b, &tmp1, s, &mut bAAcAc) != 0 {
        return false;
    }
    !(SUNRabs(bAAcAc - 1.0 / 120.0) > TOL())
}

/* b'*A1*(c1.*(A2*(c2.*c3))) = 1/90 */
#[allow(clippy::too_many_arguments)]
fn arkode_butcher_order6n(
    b: &[f64],
    A1: &[Vec<f64>],
    c1: &[f64],
    A2: &[Vec<f64>],
    c2: &[f64],
    c3: &[f64],
    s: i32,
) -> bool {
    let n = s.max(0) as usize;
    let (mut tmp1, mut tmp2) = (vec![0.0; n], vec![0.0; n]);
    let mut bAcAcc = 0.0;
    if arkode_butcher_vv(c2, c3, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_mv(A2, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_vv(c1, &tmp2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_mv(A1, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_dot(b, &tmp2, s, &mut bAcAcc) != 0 {
        return false;
    }
    !(SUNRabs(bAcAcc - 1.0 / 90.0) > TOL())
}

/* b'*A1*(c1.*(A2*A3*c2)) = 1/180 */
#[allow(clippy::too_many_arguments)]
fn arkode_butcher_order6o(
    b: &[f64],
    A1: &[Vec<f64>],
    c1: &[f64],
    A2: &[Vec<f64>],
    A3: &[Vec<f64>],
    c2: &[f64],
    s: i32,
) -> bool {
    let n = s.max(0) as usize;
    let (mut tmp1, mut tmp2) = (vec![0.0; n], vec![0.0; n]);
    let mut bAcAAc = 0.0;
    if arkode_butcher_mv(A3, c2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_mv(A2, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_vv(c1, &tmp2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_mv(A1, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_dot(b, &tmp2, s, &mut bAcAAc) != 0 {
        return false;
    }
    !(SUNRabs(bAcAAc - 1.0 / 180.0) > TOL())
}

/* b'*A1*A2*(c1.*c2.*c3) = 1/120 */
#[allow(clippy::too_many_arguments)]
fn arkode_butcher_order6p(
    b: &[f64],
    A1: &[Vec<f64>],
    A2: &[Vec<f64>],
    c1: &[f64],
    c2: &[f64],
    c3: &[f64],
    s: i32,
) -> bool {
    let n = s.max(0) as usize;
    let (mut tmp1, mut tmp2) = (vec![0.0; n], vec![0.0; n]);
    let mut bAAccc = 0.0;
    if arkode_butcher_vv(c1, c2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_vv(c3, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_mv(A2, &tmp2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_mv(A1, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_dot(b, &tmp2, s, &mut bAAccc) != 0 {
        return false;
    }
    !(SUNRabs(bAAccc - 1.0 / 120.0) > TOL())
}

/* b'*A1*A2*(c1.*(A3*c2)) = 1/240 */
#[allow(clippy::too_many_arguments)]
fn arkode_butcher_order6q(
    b: &[f64],
    A1: &[Vec<f64>],
    A2: &[Vec<f64>],
    c1: &[f64],
    A3: &[Vec<f64>],
    c2: &[f64],
    s: i32,
) -> bool {
    let n = s.max(0) as usize;
    let (mut tmp1, mut tmp2) = (vec![0.0; n], vec![0.0; n]);
    let mut bAAcAc = 0.0;
    if arkode_butcher_mv(A3, c2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_vv(c1, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_mv(A2, &tmp2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_mv(A1, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_dot(b, &tmp2, s, &mut bAAcAc) != 0 {
        return false;
    }
    !(SUNRabs(bAAcAc - 1.0 / 240.0) > TOL())
}

/* b'*A1*A2*A3*(c1.*c2) = 1/360 */
#[allow(clippy::too_many_arguments)]
fn arkode_butcher_order6r(
    b: &[f64],
    A1: &[Vec<f64>],
    A2: &[Vec<f64>],
    A3: &[Vec<f64>],
    c1: &[f64],
    c2: &[f64],
    s: i32,
) -> bool {
    let n = s.max(0) as usize;
    let (mut tmp1, mut tmp2) = (vec![0.0; n], vec![0.0; n]);
    let mut bAAAcc = 0.0;
    if arkode_butcher_vv(c1, c2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_mv(A3, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_mv(A2, &tmp2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_mv(A1, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_dot(b, &tmp2, s, &mut bAAAcc) != 0 {
        return false;
    }
    !(SUNRabs(bAAAcc - 1.0 / 360.0) > TOL())
}

/* b'*A1*A2*A3*A4*c = 1/720 */
#[allow(clippy::too_many_arguments)]
fn arkode_butcher_order6s(
    b: &[f64],
    A1: &[Vec<f64>],
    A2: &[Vec<f64>],
    A3: &[Vec<f64>],
    A4: &[Vec<f64>],
    c: &[f64],
    s: i32,
) -> bool {
    let n = s.max(0) as usize;
    let (mut tmp1, mut tmp2) = (vec![0.0; n], vec![0.0; n]);
    let mut bAAAAc = 0.0;
    if arkode_butcher_mv(A4, c, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_mv(A3, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_mv(A2, &tmp2, s, &mut tmp1) != 0 {
        return false;
    }
    if arkode_butcher_mv(A1, &tmp1, s, &mut tmp2) != 0 {
        return false;
    }
    if arkode_butcher_dot(b, &tmp2, s, &mut bAAAAc) != 0 {
        return false;
    }
    !(SUNRabs(bAAAAc - 1.0 / 720.0) > TOL())
}

/*---------------------------------------------------------------
  Utility routine to check Butcher's simplifying assumptions.
  Returns the maximum predicted order.
  ---------------------------------------------------------------*/
#[allow(non_snake_case)]
fn __ButcherSimplifyingAssumptions(A: &[Vec<f64>], b: &[f64], c: &[f64], s: i32) -> i32 {
    let n = s.max(0) as usize;
    let mut tmp = vec![0.0; n];

    /* B(P) */
    let mut P = 0;
    for i in 1..1000 {
        if arkode_butcher_vp(c, i - 1, s, &mut tmp) != 0 {
            return 0;
        }
        let mut LHS = 0.0;
        if arkode_butcher_dot(b, &tmp, s, &mut LHS) != 0 {
            return 0;
        }
        let RHS = 1.0 / (i as f64);
        if SUNRabs(RHS - LHS) > TOL() {
            break;
        }
        P += 1;
    }

    /* C(Q) */
    let mut Q = 0;
    for k in 1..1000 {
        let mut alltrue = true;
        for i in 0..n {
            if arkode_butcher_vp(c, k - 1, s, &mut tmp) != 0 {
                return 0;
            }
            let mut LHS = 0.0;
            if arkode_butcher_dot(&A[i], &tmp, s, &mut LHS) != 0 {
                return 0;
            }
            let RHS = SUNRpowerI(c[i], k) / (k as f64);
            if SUNRabs(RHS - LHS) > TOL() {
                alltrue = false;
                break;
            }
        }
        if alltrue {
            Q += 1;
        } else {
            break;
        }
    }

    /* D(R) */
    let mut R = 0;
    for k in 1..1000 {
        let mut alltrue = true;
        for j in 0..n {
            let mut LHS = 0.0;
            for i in 0..n {
                LHS += A[i][j] * b[i] * SUNRpowerI(c[i], k - 1);
            }
            let RHS = b[j] / (k as f64) * (1.0 - SUNRpowerI(c[j], k));
            if SUNRabs(RHS - LHS) > TOL() {
                alltrue = false;
                break;
            }
        }
        if alltrue {
            R += 1;
        } else {
            break;
        }
    }

    /* determine q, clean up and return */
    let mut q = 0;
    for _i in 1..=P {
        if (q > Q + R + 1) || (q > 2 * Q + 2) {
            break;
        }
        q += 1;
    }
    q
}

#[cfg(test)]
mod tests {
    use super::*;

    /* classic RK4 with a trivial embedding vector */
    fn rk4() -> ARKodeButcherTable {
        let c = [0.0, 0.5, 0.5, 1.0];
        #[rustfmt::skip]
        let a = [
            0.0, 0.0, 0.0, 0.0,
            0.5, 0.0, 0.0, 0.0,
            0.0, 0.5, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0,
        ];
        let b = [1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0];
        ARKodeButcherTable_Create(4, 4, 0, &c, &a, &b, None).unwrap()
    }

    #[test]
    fn alloc_create_copy_space() {
        assert!(ARKodeButcherTable_Alloc(0, false).is_none());
        let B = rk4();
        assert_eq!((B.q, B.p, B.stages), (4, 0, 4));
        assert_eq!(B.A[1][0], 0.5);
        assert!(B.d.is_none());

        let Bc = ARKodeButcherTable_Copy(&B).unwrap();
        assert_eq!(Bc.b, B.b);
        assert_eq!(Bc.A, B.A);

        let (mut liw, mut lrw) = (0i64, 0i64);
        ARKodeButcherTable_Space(&B, &mut liw, &mut lrw);
        assert_eq!((liw, lrw), (3, 4 * 6));

        /* embedded table gets stages*(stages+3) */
        let d = [0.25; 4];
        let Be = ARKodeButcherTable_Create(4, 4, 3, &B.c, &[0.0; 16], &B.b, Some(&d)).unwrap();
        ARKodeButcherTable_Space(&Be, &mut liw, &mut lrw);
        assert_eq!((liw, lrw), (3, 4 * 7));
        ARKodeButcherTable_Free(Some(Be));
    }

    #[test]
    fn stiffly_accurate() {
        /* backward Euler: A = [1], b = [1] — stiffly accurate */
        let B = ARKodeButcherTable_Create(1, 1, 0, &[1.0], &[1.0], &[1.0], None).unwrap();
        assert!(ARKodeButcherTable_IsStifflyAccurate(&B));
        /* RK4 is not */
        assert!(!ARKodeButcherTable_IsStifflyAccurate(&rk4()));
    }

    #[test]
    fn write_format() {
        let B = ARKodeButcherTable_Create(1, 1, 0, &[0.5], &[0.25], &[1.0], None).unwrap();
        let mut out = Vec::new();
        ARKodeButcherTable_Write(&B, &mut out);
        let s = String::from_utf8(out).unwrap();
        assert_eq!(
            s,
            "  A = \n       2.500000000000000e-01  \n  c =  5.000000000000000e-01  \n  b =  1.000000000000000e+00  \n"
        );
    }

    /* ------------------------- PART II ------------------------- */

    /* Heun-Euler embedded pair: q=2 method, p=1 embedding */
    fn heun_euler() -> ARKodeButcherTable {
        let c = [0.0, 1.0];
        let a = [0.0, 0.0, 1.0, 0.0];
        let b = [0.5, 0.5];
        let d = [1.0, 0.0];
        ARKodeButcherTable_Create(2, 2, 1, &c, &a, &b, Some(&d)).unwrap()
    }

    #[test]
    fn check_order_rk4() {
        let B = rk4();
        let (mut q, mut p) = (0, 0);
        assert_eq!(ARKodeButcherTable_CheckOrder(&B, &mut q, &mut p, None), 0);
        assert_eq!((q, p), (4, 0));
    }

    #[test]
    fn check_order_heun_euler_embedded() {
        let B = heun_euler();
        let (mut q, mut p) = (0, 0);
        assert_eq!(ARKodeButcherTable_CheckOrder(&B, &mut q, &mut p, None), 0);
        assert_eq!((q, p), (2, 1));
    }

    #[test]
    fn check_order_flags() {
        /* declared order higher than analytical -> failure (-1) */
        let mut B = heun_euler();
        B.q = 3;
        let (mut q, mut p) = (0, 0);
        assert_eq!(ARKodeButcherTable_CheckOrder(&B, &mut q, &mut p, None), -1);
        assert_eq!(q, 2);

        /* declared order lower than analytical -> warning (1) */
        B.q = 1;
        assert_eq!(ARKodeButcherTable_CheckOrder(&B, &mut q, &mut p, None), 1);
        assert_eq!(q, 2);

        /* embedding declared too high -> failure */
        B.q = 2;
        B.p = 2;
        assert_eq!(ARKodeButcherTable_CheckOrder(&B, &mut q, &mut p, None), -1);
        assert_eq!(p, 1);
    }

    #[test]
    fn check_order_output_text() {
        /* explicit Euler: passes order 1, fails order 2 */
        let B = ARKodeButcherTable_Create(1, 1, 0, &[0.0], &[0.0], &[1.0], None).unwrap();
        let (mut q, mut p) = (0, 0);
        let mut out: Vec<u8> = Vec::new();
        assert_eq!(
            ARKodeButcherTable_CheckOrder(&B, &mut q, &mut p, Some(&mut out)),
            0
        );
        assert_eq!((q, p), (1, 0));
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "ARKodeButcherTable_CheckOrder:\n  method fails order 2 condition\n"
        );
    }

    #[test]
    fn check_ark_order_identical_pair() {
        /* B1 = B2 = RK4: all cross conditions collapse to the plain ones */
        let B1 = rk4();
        let B2 = rk4();
        let (mut q, mut p) = (0, 0);
        assert_eq!(
            ARKodeButcherTable_CheckARKOrder(&B1, &B2, &mut q, &mut p, None),
            0
        );
        assert_eq!((q, p), (4, 0));

        /* stage mismatch is rejected */
        let B3 = heun_euler();
        assert_eq!(
            ARKodeButcherTable_CheckARKOrder(&B1, &B3, &mut q, &mut p, None),
            -1
        );
    }

    #[test]
    fn check_ark_order_embedding_and_quirk() {
        /* B1 = B2 = Heun-Euler: embedding path measures p=1 */
        let B1 = heun_euler();
        let B2 = heun_euler();
        let (mut q, mut p) = (0, 0);
        assert_eq!(
            ARKodeButcherTable_CheckARKOrder(&B1, &B2, &mut q, &mut p, None),
            0
        );
        assert_eq!((q, p), (2, 1));

        /* C quirk (d[1] = B1->d): embedding is gated on B1 alone, so a
        non-embedded B1 skips the embedding check even though B2 has d */
        let mut B1n = heun_euler();
        B1n.d = None;
        B1n.p = 0;
        assert_eq!(
            ARKodeButcherTable_CheckARKOrder(&B1n, &B2, &mut q, &mut p, None),
            0
        );
        assert_eq!((q, p), (2, 0));
    }

    #[test]
    fn simplifying_assumptions_gauss_legendre_2() {
        /* 2-stage Gauss-Legendre collocation: order 4 with B(P)/C(Q)/D(R)
        satisfied to P=4, Q=2, R=2 -> the SA routine itself reports 4 */
        let r3 = 3.0_f64.sqrt();
        let a = vec![
            vec![0.25, 0.25 - r3 / 6.0],
            vec![0.25 + r3 / 6.0, 0.25],
        ];
        let b = [0.5, 0.5];
        let c = [0.5 - r3 / 6.0, 0.5 + r3 / 6.0];
        assert_eq!(__ButcherSimplifyingAssumptions(&a, &b, &c, 2), 4);
    }
}
