/* -----------------------------------------------------------------
 * Translation of src/arkode/arkode_splittingstep_coefficients.c
 * (+ the arkode_splittingstep_coefficients.def X-macro table and the
 *  coefficient parts of include/arkode/arkode_splittingstep.h),
 * SUNDIALS 7.7.0.
 *
 * This is the implementation file for splitting coefficients.
 *
 * C's SplittingStepCoefficients is a pointer to a Mem struct whose
 * beta tensor is one contiguous calloc with two layers of row
 * pointers so it can be indexed beta[i][j][k]; the Rust port models
 * it as nested Vecs with the same indexing (sequential method i,
 * stage row j in 0..=stages, partition k).  Functions that return
 * NULL on failure return Option<SplittingStepCoefficients> (a boxed
 * Mem).  C NULL-pointer argument checks with no Rust equivalent
 * (alpha/beta arrays) are noted at their sites.
 * -----------------------------------------------------------------*/

use crate::arkode_impl::{arkProcessError, ARK_ILL_INPUT};
use crate::sundials_math::{SUNIpowerI, SUNRpowerR, SUNRsqrt};
use crate::sundials_utils::fmt_e;

/* enum ARKODE_SplittingCoefficientsID (arkode_splittingstep.h).
Splitting names use the convention
ARKODE_SPLITTING_<name>_<stages>_<order>_<partitions> */
pub type ARKODE_SplittingCoefficientsID = i32;
pub const ARKODE_SPLITTING_NONE: i32 = -1; /* ensure enum is signed int */
pub const ARKODE_SPLITTING_LIE_TROTTER_1_1_2: i32 = 0;
pub const ARKODE_MIN_SPLITTING_NUM: i32 = 0;
pub const ARKODE_SPLITTING_STRANG_2_2_2: i32 = 1;
pub const ARKODE_SPLITTING_BEST_2_2_2: i32 = 2;
pub const ARKODE_SPLITTING_SUZUKI_3_3_2: i32 = 3;
pub const ARKODE_SPLITTING_RUTH_3_3_2: i32 = 4;
pub const ARKODE_SPLITTING_YOSHIDA_4_4_2: i32 = 5;
pub const ARKODE_SPLITTING_YOSHIDA_8_6_2: i32 = 6;
pub const ARKODE_MAX_SPLITTING_NUM: i32 = ARKODE_SPLITTING_YOSHIDA_8_6_2;

/// struct SplittingStepCoefficientsMem
pub struct SplittingStepCoefficientsMem {
    /// weights for sum over sequential splitting methods
    pub alpha: Vec<f64>,
    /// subintegration nodes, indexed by the sequential method, stage, and partition
    pub beta: Vec<Vec<Vec<f64>>>,
    /// number of sequential splitting methods
    pub sequential_methods: i32,
    /// number of stages within each sequential splitting method
    pub stages: i32,
    /// number of RHS partitions
    pub partitions: i32,
    /// order of convergence
    pub order: i32,
}

/// C `SplittingStepCoefficients` (pointer type); NULL <-> None at use
/// sites.
pub type SplittingStepCoefficients = Box<SplittingStepCoefficientsMem>;

/*---------------------------------------------------------------
  Routine to allocate splitting coefficients with zero values for
  alpha and beta
  ---------------------------------------------------------------*/
pub fn SplittingStepCoefficients_Alloc(
    sequential_methods: i32,
    stages: i32,
    partitions: i32,
) -> Option<SplittingStepCoefficients> {
    if sequential_methods < 1 || stages < 1 || partitions < 1 {
        return None;
    }

    /* beta is calloc'd in C so only non-zero coefficients need to be set */
    Some(Box::new(SplittingStepCoefficientsMem {
        alpha: vec![0.0; sequential_methods as usize],
        beta: vec![
            vec![vec![0.0; partitions as usize]; (stages + 1) as usize];
            sequential_methods as usize
        ],
        sequential_methods,
        stages,
        partitions,
        order: 0,
    }))
}

/*---------------------------------------------------------------
  Routine to create splitting coefficients which performs a copy
  of the alpha and beta parameters (beta is passed flattened in
  C's contiguous beta[i][j][k] layout)
  ---------------------------------------------------------------*/
pub fn SplittingStepCoefficients_Create(
    sequential_methods: i32,
    stages: i32,
    partitions: i32,
    order: i32,
    alpha: &[f64],
    beta: &[f64],
) -> Option<SplittingStepCoefficients> {
    /* C also rejects alpha == NULL || beta == NULL (inexpressible here) */
    if order < 1 {
        return None;
    }

    let mut coefficients = SplittingStepCoefficients_Alloc(sequential_methods, stages, partitions)?;

    coefficients.order = order;
    coefficients.alpha[..sequential_methods as usize]
        .copy_from_slice(&alpha[..sequential_methods as usize]);
    for i in 0..sequential_methods as usize {
        for j in 0..=stages as usize {
            for k in 0..partitions as usize {
                coefficients.beta[i][j][k] =
                    beta[(i * (stages as usize + 1) + j) * partitions as usize + k];
            }
        }
    }

    Some(coefficients)
}

/*---------------------------------------------------------------
  Routine to free splitting coefficients
  ---------------------------------------------------------------*/
pub fn SplittingStepCoefficients_Destroy(coefficients: &mut Option<SplittingStepCoefficients>) {
    *coefficients = None;
}

/*---------------------------------------------------------------
  Routine to create a copy of splitting coefficients
  ---------------------------------------------------------------*/
pub fn SplittingStepCoefficients_Copy(
    coefficients: &SplittingStepCoefficientsMem,
) -> Option<SplittingStepCoefficients> {
    let mut coefficientsCopy = SplittingStepCoefficients_Alloc(
        coefficients.sequential_methods,
        coefficients.stages,
        coefficients.partitions,
    )?;

    coefficientsCopy.order = coefficients.order;
    coefficientsCopy.alpha.copy_from_slice(&coefficients.alpha);
    coefficientsCopy.beta.clone_from(&coefficients.beta);

    Some(coefficientsCopy)
}

/*---------------------------------------------------------------
  Routine to load coefficients from an ID
  (X-macro over arkode_splittingstep_coefficients.def)
  ---------------------------------------------------------------*/
pub fn SplittingStepCoefficients_LoadCoefficients(
    method: ARKODE_SplittingCoefficientsID,
) -> Option<SplittingStepCoefficients> {
    match method {
        ARKODE_SPLITTING_NONE => None,
        ARKODE_SPLITTING_LIE_TROTTER_1_1_2 => SplittingStepCoefficients_LieTrotter(2),
        ARKODE_SPLITTING_STRANG_2_2_2 => SplittingStepCoefficients_Strang(2),
        ARKODE_SPLITTING_BEST_2_2_2 => splittingStepCoefficients_Best_2_2_2(),
        ARKODE_SPLITTING_SUZUKI_3_3_2 => SplittingStepCoefficients_ThirdOrderSuzuki(2),
        ARKODE_SPLITTING_RUTH_3_3_2 => splittingStepCoefficients_Ruth_3_3_2(),
        ARKODE_SPLITTING_YOSHIDA_4_4_2 => SplittingStepCoefficients_TripleJump(2, 4),
        ARKODE_SPLITTING_YOSHIDA_8_6_2 => SplittingStepCoefficients_TripleJump(2, 6),
        _ => {
            arkProcessError(
                None,
                ARK_ILL_INPUT,
                line!(),
                "SplittingStepCoefficients_LoadCoefficients",
                file!(),
                "Unknown splitting coefficients",
            );
            None
        }
    }
}

/*---------------------------------------------------------------
  Routine to load coefficients using a string representation of
  an enum entry in ARKODE_SplittingCoefficientsID
  ---------------------------------------------------------------*/
pub fn SplittingStepCoefficients_LoadCoefficientsByName(
    method: &str,
) -> Option<SplittingStepCoefficients> {
    match method {
        "ARKODE_SPLITTING_NONE" => None,
        "ARKODE_SPLITTING_LIE_TROTTER_1_1_2" => SplittingStepCoefficients_LieTrotter(2),
        "ARKODE_SPLITTING_STRANG_2_2_2" => SplittingStepCoefficients_Strang(2),
        "ARKODE_SPLITTING_BEST_2_2_2" => splittingStepCoefficients_Best_2_2_2(),
        "ARKODE_SPLITTING_SUZUKI_3_3_2" => SplittingStepCoefficients_ThirdOrderSuzuki(2),
        "ARKODE_SPLITTING_RUTH_3_3_2" => splittingStepCoefficients_Ruth_3_3_2(),
        "ARKODE_SPLITTING_YOSHIDA_4_4_2" => SplittingStepCoefficients_TripleJump(2, 4),
        "ARKODE_SPLITTING_YOSHIDA_8_6_2" => SplittingStepCoefficients_TripleJump(2, 6),
        _ => {
            arkProcessError(
                None,
                ARK_ILL_INPUT,
                line!(),
                "SplittingStepCoefficients_LoadCoefficientsByName",
                file!(),
                "Unknown splitting coefficients",
            );
            None
        }
    }
}

/*---------------------------------------------------------------
  Routine to convert a coefficient enum value to its string
  representation
  ---------------------------------------------------------------*/
pub fn SplittingStepCoefficients_IDToName(
    id: ARKODE_SplittingCoefficientsID,
) -> Option<&'static str> {
    match id {
        ARKODE_SPLITTING_NONE => Some("ARKODE_SPLITTING_NONE"),
        ARKODE_SPLITTING_LIE_TROTTER_1_1_2 => Some("ARKODE_SPLITTING_LIE_TROTTER_1_1_2"),
        ARKODE_SPLITTING_STRANG_2_2_2 => Some("ARKODE_SPLITTING_STRANG_2_2_2"),
        ARKODE_SPLITTING_BEST_2_2_2 => Some("ARKODE_SPLITTING_BEST_2_2_2"),
        ARKODE_SPLITTING_SUZUKI_3_3_2 => Some("ARKODE_SPLITTING_SUZUKI_3_3_2"),
        ARKODE_SPLITTING_RUTH_3_3_2 => Some("ARKODE_SPLITTING_RUTH_3_3_2"),
        ARKODE_SPLITTING_YOSHIDA_4_4_2 => Some("ARKODE_SPLITTING_YOSHIDA_4_4_2"),
        ARKODE_SPLITTING_YOSHIDA_8_6_2 => Some("ARKODE_SPLITTING_YOSHIDA_8_6_2"),
        _ => {
            arkProcessError(
                None,
                ARK_ILL_INPUT,
                line!(),
                "SplittingStepCoefficients_IDToName",
                file!(),
                "Unknown splitting coefficients",
            );
            None
        }
    }
}

/* .def inline table entry ARKODE_SPLITTING_BEST_2_2_2 */
fn splittingStepCoefficients_Best_2_2_2() -> Option<SplittingStepCoefficients> {
    let mut coefficients = SplittingStepCoefficients_Alloc(1, 2, 2)?;
    coefficients.order = 2;
    coefficients.alpha[0] = 1.0;
    coefficients.beta[0][1][0] = 1.0 - SUNRsqrt(0.5);
    coefficients.beta[0][1][1] = SUNRsqrt(0.5);
    coefficients.beta[0][2][0] = 1.0;
    coefficients.beta[0][2][1] = 1.0;
    Some(coefficients)
}

/* .def inline table entry ARKODE_SPLITTING_RUTH_3_3_2 */
fn splittingStepCoefficients_Ruth_3_3_2() -> Option<SplittingStepCoefficients> {
    let mut coefficients = SplittingStepCoefficients_Alloc(1, 3, 2)?;
    coefficients.order = 3;
    coefficients.alpha[0] = 1.0;
    coefficients.beta[0][1][0] = 1.0;
    coefficients.beta[0][1][1] = -1.0 / 24.0;
    coefficients.beta[0][2][0] = 1.0 / 3.0;
    coefficients.beta[0][2][1] = 17.0 / 24.0;
    coefficients.beta[0][3][0] = 1.0;
    coefficients.beta[0][3][1] = 1.0;
    Some(coefficients)
}

/*---------------------------------------------------------------
  Routine to construct the standard Lie-Trotter splitting
  ---------------------------------------------------------------*/
pub fn SplittingStepCoefficients_LieTrotter(partitions: i32) -> Option<SplittingStepCoefficients> {
    let mut coefficients = SplittingStepCoefficients_Alloc(1, 1, partitions)?;

    coefficients.order = 1;
    coefficients.alpha[0] = 1.0;
    for i in 0..partitions as usize {
        coefficients.beta[0][1][i] = 1.0;
    }

    Some(coefficients)
}

/*---------------------------------------------------------------
  Routine to construct the standard Stang splitting
  ---------------------------------------------------------------*/
pub fn SplittingStepCoefficients_Strang(partitions: i32) -> Option<SplittingStepCoefficients> {
    SplittingStepCoefficients_TripleJump(partitions, 2)
}

/*---------------------------------------------------------------
  Routine to construct a parallel splitting method
  Phi_1(h) + Phi_2(h) + ... + Phi_p(h) - (p - 1) * y_n
  where Phi_i is the flow of partition i and p = partitions.
  ---------------------------------------------------------------*/
pub fn SplittingStepCoefficients_Parallel(partitions: i32) -> Option<SplittingStepCoefficients> {
    let mut coefficients = SplittingStepCoefficients_Alloc(partitions + 1, 1, partitions)?;

    coefficients.order = 1;
    for i in 0..partitions as usize {
        coefficients.alpha[i] = 1.0;
        coefficients.beta[i][1][i] = 1.0;
    }

    coefficients.alpha[partitions as usize] = (1 - partitions) as f64;

    Some(coefficients)
}

/*---------------------------------------------------------------
  Routine to construct a symmetric parallel splitting which is
  the average of the Lie-Trotter method and its adjoint
  ---------------------------------------------------------------*/
pub fn SplittingStepCoefficients_SymmetricParallel(
    partitions: i32,
) -> Option<SplittingStepCoefficients> {
    let mut coefficients = SplittingStepCoefficients_Alloc(2, partitions, partitions)?;

    coefficients.order = 2;
    coefficients.alpha[0] = 0.5;
    coefficients.alpha[1] = 0.5;

    for i in 0..partitions as usize {
        coefficients.beta[0][partitions as usize][i] = 1.0;
        for j in (partitions as usize - i - 1)..partitions as usize {
            coefficients.beta[1][i + 1][j] = 1.0;
        }
    }

    Some(coefficients)
}

/*---------------------------------------------------------------
  Routine to construct a 3rd order method of Suzuki of the form
  L(p1 h) * L*(p2 h) * L(p3 h) * L*(p4 h) * L(p5 h)
  where L is a Lie-Trotter splitting and L* is its adjoint.
  Composition is denoted by *.
  ---------------------------------------------------------------*/
pub fn SplittingStepCoefficients_ThirdOrderSuzuki(
    partitions: i32,
) -> Option<SplittingStepCoefficients> {
    let mut coefficients = SplittingStepCoefficients_Alloc(1, 2 * partitions - 1, partitions)?;

    coefficients.order = 3;
    coefficients.alpha[0] = 1.0;

    for i in 1..partitions {
        for j in 0..partitions {
            // Constants from https://doi.org/10.1143/JPSJ.61.3015 pg. 3019
            let p1: f64 = 0.2683300957817599249569552299254991394812;
            let p2: f64 = 0.6513314272356399320939424082278836500821;

            coefficients.beta[0][i as usize][j as usize] =
                if i + j < partitions { p1 } else { p1 + p2 };
            coefficients.beta[0][(partitions + i - 1) as usize][j as usize] =
                1.0 - (if i + j < partitions { p1 + p2 } else { p1 });
        }
    }

    for i in 0..partitions as usize {
        coefficients.beta[0][(2 * partitions - 1) as usize][i] = 1.0;
    }

    Some(coefficients)
}

/*---------------------------------------------------------------
  Routine to construct a composition method of the form
  S(gamma_0 h)^c * S(gamma_1 h) * S(gamma_0)^c
  where S is a lower order splitting (with Stang as the base case),
  * and ^ denote composition, and c = composition_stages. This
  covers both the triple jump (c=1) and Suzuki fractal (c=2).

  C walks a row pointer through the beta[0] matrix and returns the
  advanced pointer; the Rust port passes the row slice and returns
  the row offset the next call starts at (relative to that slice).
  ---------------------------------------------------------------*/
fn splittingStepCoefficients_ComposeStrangHelper(
    partitions: i32,
    order: i32,
    composition_stages: i32,
    start: f64,
    end: f64,
    beta: &mut [Vec<f64>],
) -> usize {
    let diff = end - start;
    if order == 2 {
        /* The base case is an order 2 Strang splitting */
        let mid = start + diff / 2.0;
        for j in 1..=partitions as usize {
            for k in 0..partitions as usize {
                beta[j][k] = if (k + j) < partitions as usize { mid } else { end };
            }
        }

        return (partitions - 1) as usize;
    }

    let mut beta_cur: usize = 0;
    let mut start_cur = start;
    /* This is essentially the gamma coefficient from Geometric Numerical
     * Integration (https://doi.org/10.1007/3-540-30666-8) pg 44-45 scaled by the
     * current interval */
    let gamma = diff
        / ((composition_stages - 1) as f64
            - SUNRpowerR((composition_stages - 1) as f64, 1.0 / (order - 1) as f64));
    for i in 1..=composition_stages {
        /* To avoid roundoff issues, this ensures end_cur=1 for the last value of i*/
        let end_cur = if 2 * i < composition_stages {
            start + i as f64 * gamma
        } else {
            end + (i - composition_stages) as f64 * gamma
        };
        /* Recursively generate coefficients and shift beta_cur */
        beta_cur += splittingStepCoefficients_ComposeStrangHelper(
            partitions,
            order - 2,
            composition_stages,
            start_cur,
            end_cur,
            &mut beta[beta_cur..],
        );
        start_cur = end_cur;
    }

    beta_cur
}

/*---------------------------------------------------------------
  Routine which does validation and setup before calling
  SplittingStepCoefficients_ComposeStrangHelper to fill in the
  beta coefficients
  ---------------------------------------------------------------*/
fn splittingStepCoefficients_ComposeStrang(
    partitions: i32,
    order: i32,
    composition_stages: i32,
) -> Option<SplittingStepCoefficients> {
    if order < 2 || order % 2 != 0 {
        // Only even orders allowed
        return None;
    }

    let stages = 1 + (partitions - 1) * SUNIpowerI(composition_stages, order / 2 - 1);
    let mut coefficients = SplittingStepCoefficients_Alloc(1, stages, partitions)?;

    coefficients.order = order;
    coefficients.alpha[0] = 1.0;

    splittingStepCoefficients_ComposeStrangHelper(
        partitions,
        order,
        composition_stages,
        0.0,
        1.0,
        &mut coefficients.beta[0],
    );

    Some(coefficients)
}

pub fn SplittingStepCoefficients_TripleJump(
    partitions: i32,
    order: i32,
) -> Option<SplittingStepCoefficients> {
    splittingStepCoefficients_ComposeStrang(partitions, order, 3)
}

pub fn SplittingStepCoefficients_SuzukiFractal(
    partitions: i32,
    order: i32,
) -> Option<SplittingStepCoefficients> {
    splittingStepCoefficients_ComposeStrang(partitions, order, 5)
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
  Routine to print a splitting coefficient structure
  ---------------------------------------------------------------*/
pub fn SplittingStepCoefficients_Write(
    coefficients: &SplittingStepCoefficientsMem,
    outfile: &mut dyn std::io::Write,
) {
    /* C also guards outfile/coefficients/alpha/beta against NULL */

    let _ = writeln!(
        outfile,
        "  sequential methods = {}",
        coefficients.sequential_methods
    );
    let _ = writeln!(outfile, "  stages = {}", coefficients.stages);
    let _ = writeln!(outfile, "  partitions = {}", coefficients.partitions);
    let _ = writeln!(outfile, "  order = {}", coefficients.order);
    let _ = write!(outfile, "  alpha = ");
    for i in 0..coefficients.sequential_methods as usize {
        let _ = write!(outfile, "{}  ", sun_format_e(coefficients.alpha[i]));
    }
    let _ = writeln!(outfile);

    for i in 0..coefficients.sequential_methods as usize {
        let _ = writeln!(outfile, "  beta[{}] = ", i);
        for j in 0..=coefficients.stages as usize {
            let _ = write!(outfile, "      ");
            for k in 0..coefficients.partitions as usize {
                let _ = write!(outfile, "{}  ", sun_format_e(coefficients.beta[i][j][k]));
            }
            let _ = writeln!(outfile);
        }
        let _ = writeln!(outfile);
    }
}
