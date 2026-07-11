/* -----------------------------------------------------------------
 * Translation of sundials-7.7.0/src/arkode/arkode_butcher.c
 * (+ include/arkode/arkode_butcher.h).
 *
 * PART I: the ARKodeButcherTable object (Alloc/Create/Copy/Space/
 * Free/Write/IsStifflyAccurate). PART II (the CheckOrder /
 * CheckARKOrder analysis with its order-condition helper family)
 * follows in this module.
 *
 * Modeling: C's ARKodeButcherTable is a heap pointer that may be
 * NULL — Rust functions return Option<ARKodeButcherTable>.
 * `sunrealtype** A` → Vec<Vec<f64>> (row-major, stages×stages);
 * `d` (embedding) is Some only for embedded tables, exactly C's
 * NULL test. ARKodeButcherTable_Free = drop.
 * SUN_FORMAT_E is "% .15e": fmt_e with a leading space for
 * non-negative values (established pattern).
 * -----------------------------------------------------------------*/

use crate::sundials_types::SUN_UNIT_ROUNDOFF;
use crate::sundials_utils::fmt_e;

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
}
