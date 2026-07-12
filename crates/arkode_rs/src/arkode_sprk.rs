/* -----------------------------------------------------------------
 * Translation of src/arkode/arkode_sprk.c and the type definitions
 * of include/arkode/arkode_sprk.h (ARKODE 7.7.0).
 * Symplectic partitioned Runge-Kutta coefficient tables.
 * -----------------------------------------------------------------*/

use crate::arkode_butcher::{ARKodeButcherTable, ARKodeButcherTable_Alloc};
use crate::sundials_math::{SUNRpowerR, SUNRsqrt};

/* ARKODE_SPRKMethodID (arkode_sprk.h) */
pub type ARKODE_SPRKMethodID = i32;
pub const ARKODE_SPRK_NONE: ARKODE_SPRKMethodID = -1;
pub const ARKODE_SPRK_EULER_1_1: ARKODE_SPRKMethodID = 0;
pub const ARKODE_MIN_SPRK_NUM: ARKODE_SPRKMethodID = 0;
pub const ARKODE_SPRK_LEAPFROG_2_2: ARKODE_SPRKMethodID = 1;
pub const ARKODE_SPRK_PSEUDO_LEAPFROG_2_2: ARKODE_SPRKMethodID = 2;
pub const ARKODE_SPRK_RUTH_3_3: ARKODE_SPRKMethodID = 3;
pub const ARKODE_SPRK_MCLACHLAN_2_2: ARKODE_SPRKMethodID = 4;
pub const ARKODE_SPRK_MCLACHLAN_3_3: ARKODE_SPRKMethodID = 5;
pub const ARKODE_SPRK_CANDY_ROZMUS_4_4: ARKODE_SPRKMethodID = 6;
pub const ARKODE_SPRK_MCLACHLAN_4_4: ARKODE_SPRKMethodID = 7;
pub const ARKODE_SPRK_MCLACHLAN_5_6: ARKODE_SPRKMethodID = 8;
pub const ARKODE_SPRK_YOSHIDA_6_8: ARKODE_SPRKMethodID = 9;
pub const ARKODE_SPRK_SUZUKI_UMENO_8_16: ARKODE_SPRKMethodID = 10;
pub const ARKODE_SPRK_SOFRONIOU_10_36: ARKODE_SPRKMethodID = 11;
pub const ARKODE_MAX_SPRK_NUM: ARKODE_SPRKMethodID = ARKODE_SPRK_SOFRONIOU_10_36;

/// struct ARKodeSPRKTableMem (arkode_sprk.h); None = C NULL handle.
#[derive(Clone)]
pub struct ARKodeSPRKTable {
    /* method order of accuracy */
    pub q: i32,
    /* number of stages */
    pub stages: i32,
    /* the a_i coefficients generate the explicit Butcher table */
    pub a: Vec<f64>,
    /* the ahat_i coefficients generate the diagonally-implicit Butcher table */
    pub ahat: Vec<f64>,
}

fn arkodeSymplecticEuler() -> ARKodeSPRKTable {
    let mut sprk_table = ARKodeSPRKTable_Alloc(1);
    sprk_table.q = 1;
    sprk_table.stages = 1;
    sprk_table.a[0] = 1.0;
    sprk_table.ahat[0] = 1.0;
    sprk_table
}

/*
  The following methods are from:

  J Candy, W Rozmus, A symplectic integration algorithm for separable
  Hamiltonian functions, Journal of Computational Physics, Volume 92,
  Issue 1, 1991, Pages 230-256, ISSN 0021-9991,
  https://doi.org/10.1016/0021-9991(91)90299-Z.
 */

fn arkodeSymplecticLeapfrog2() -> ARKodeSPRKTable {
    let mut sprk_table = ARKodeSPRKTable_Alloc(2);
    sprk_table.q = 2;
    sprk_table.stages = 2;
    sprk_table.a[0] = 0.5;
    sprk_table.a[1] = 0.5;
    sprk_table.ahat[0] = 0.0;
    sprk_table.ahat[1] = 1.0;
    sprk_table
}

fn arkodeSymplecticPseudoLeapfrog2() -> ARKodeSPRKTable {
    let mut sprk_table = ARKodeSPRKTable_Alloc(2);
    sprk_table.q = 2;
    sprk_table.stages = 2;
    sprk_table.a[0] = 1.0;
    sprk_table.a[1] = 0.0;
    sprk_table.ahat[0] = 0.5;
    sprk_table.ahat[1] = 0.5;
    sprk_table
}

fn arkodeSymplecticCandyRozmus4() -> ARKodeSPRKTable {
    let mut sprk_table = ARKodeSPRKTable_Alloc(4);
    sprk_table.q = 4;
    sprk_table.stages = 4;
    sprk_table.a[0] =
        (2.0 + SUNRpowerR(2.0, 1.0 / 3.0) + SUNRpowerR(2.0, -1.0 / 3.0)) / 6.0;
    sprk_table.a[1] =
        (1.0 - SUNRpowerR(2.0, 1.0 / 3.0) - SUNRpowerR(2.0, -1.0 / 3.0)) / 6.0;
    sprk_table.a[2] = sprk_table.a[1];
    sprk_table.a[3] = sprk_table.a[0];
    sprk_table.ahat[0] = 0.0;
    sprk_table.ahat[1] = 1.0 / (2.0 - SUNRpowerR(2.0, 1.0 / 3.0));
    sprk_table.ahat[2] = 1.0 / (1.0 - SUNRpowerR(2.0, 2.0 / 3.0));
    sprk_table.ahat[3] = sprk_table.ahat[1];
    sprk_table
}

/*
  The following methods are from:

  Ruth, R. D. (1983). A CANONICAL INTEGRATION TECHNIQUE.
  IEEE Transactions on Nuclear Science, 30(4).
  https://accelconf.web.cern.ch/p83/PDF/PAC1983_2669.PDF
 */

fn arkodeSymplecticRuth3() -> ARKodeSPRKTable {
    let mut sprk_table = ARKodeSPRKTable_Alloc(3);
    sprk_table.q = 3;
    sprk_table.stages = 3;
    sprk_table.a[0] = 2.0 / 3.0;
    sprk_table.a[1] = -2.0 / 3.0;
    sprk_table.a[2] = 1.0;
    sprk_table.ahat[0] = 7.0 / 24.0;
    sprk_table.ahat[1] = 3.0 / 4.0;
    sprk_table.ahat[2] = -1.0 / 24.0;
    sprk_table
}

/*
  The following methods are from:

  McLachlan, R.I., Atela, P.: The accuracy of symplectic integrators.
  Nonlinearity. 5, 541-562 (1992).
  https://doi.org/10.1088/0951-7715/5/2/011
 */

fn arkodeSymplecticMcLachlan2() -> ARKodeSPRKTable {
    let mut sprk_table = ARKodeSPRKTable_Alloc(2);
    sprk_table.q = 2;
    sprk_table.stages = 2;
    sprk_table.a[1] = 1.0 - (1.0 / 2.0) * SUNRsqrt(2.0);
    sprk_table.a[0] = 1.0 - sprk_table.a[1];
    sprk_table.ahat[1] = 1.0 / (2.0 * (1.0 - sprk_table.a[1]));
    sprk_table.ahat[0] = 1.0 - sprk_table.ahat[1];
    sprk_table
}

fn arkodeSymplecticMcLachlan3() -> ARKodeSPRKTable {
    let mut sprk_table = ARKodeSPRKTable_Alloc(3);

    sprk_table.q = 3;
    sprk_table.stages = 3;

    let z = -SUNRpowerR(
        (2.0 / 27.0) - 1.0 / (9.0 * SUNRsqrt(3.0)),
        1.0 / 3.0,
    );
    let w = -2.0 / 3.0 + 1.0 / (9.0 * z) + z;
    let y = (1.0 + w * w) / 4.0;
    sprk_table.a[0] =
        SUNRsqrt(1.0 / (9.0 * y) - w / 2.0 + SUNRsqrt(y)) - 1.0 / (3.0 * SUNRsqrt(y));
    sprk_table.a[1] = 0.25 / sprk_table.a[0] - sprk_table.a[0] / 2.0;
    sprk_table.a[2] = 1.0 - sprk_table.a[0] - sprk_table.a[1];
    sprk_table.ahat[0] = sprk_table.a[2];
    sprk_table.ahat[1] = sprk_table.a[1];
    sprk_table.ahat[2] = sprk_table.a[0];
    sprk_table
}

fn arkodeSymplecticMcLachlan4() -> ARKodeSPRKTable {
    let mut sprk_table = ARKodeSPRKTable_Alloc(4);
    sprk_table.q = 4;
    sprk_table.stages = 4;
    sprk_table.a[0] = 0.515352837431122936;
    sprk_table.a[1] = -0.085782019412973646;
    sprk_table.a[2] = 0.441583023616466524;
    sprk_table.a[3] = 0.128846158365384185;
    sprk_table.ahat[0] = 0.134496199277431089;
    sprk_table.ahat[1] = -0.224819803079420806;
    sprk_table.ahat[2] = 0.756320000515668291;
    sprk_table.ahat[3] = 0.33400360328632142;
    sprk_table
}

fn arkodeSymplecticMcLachlan5() -> ARKodeSPRKTable {
    let mut sprk_table = ARKodeSPRKTable_Alloc(6);
    sprk_table.q = 5;
    sprk_table.stages = 6;
    sprk_table.a[0] = 0.339839625839110000;
    sprk_table.a[1] = -0.088601336903027329;
    sprk_table.a[2] = 0.5858564768259621188;
    sprk_table.a[3] = -0.603039356536491888;
    sprk_table.a[4] = 0.3235807965546976394;
    sprk_table.a[5] = 0.4423637942197494587;
    sprk_table.ahat[0] = 0.1193900292875672758;
    sprk_table.ahat[1] = 0.6989273703824752308;
    sprk_table.ahat[2] = -0.1713123582716007754;
    sprk_table.ahat[3] = 0.4012695022513534480;
    sprk_table.ahat[4] = 0.0107050818482359840;
    sprk_table.ahat[5] = -0.0589796254980311632;
    sprk_table
}

/*
  The following methods are from:

  Yoshida, H.: Construction of higher order symplectic integrators.
  Phys Lett A. 150, 262-268 (1990).
  https://doi.org/10.1016/0375-9601(90)90092-3
 */

fn arkodeSymplecticYoshida6() -> ARKodeSPRKTable {
    let mut sprk_table = ARKodeSPRKTable_Alloc(8);
    sprk_table.q = 6;
    sprk_table.stages = 8;
    sprk_table.a[0] = 0.7845136104775572638194976338663498757768;
    sprk_table.a[1] = 0.2355732133593581336847931829785346016865;
    sprk_table.a[2] = -1.177679984178871006946415680964315734639;
    sprk_table.a[3] = 1.315186320683911218884249728238862514352;
    sprk_table.a[4] = sprk_table.a[2];
    sprk_table.a[5] = sprk_table.a[1];
    sprk_table.a[6] = sprk_table.a[0];
    sprk_table.a[7] = 0.0;
    sprk_table.ahat[0] = sprk_table.a[0] / 2.0;
    sprk_table.ahat[1] = (sprk_table.a[0] + sprk_table.a[1]) / 2.0;
    sprk_table.ahat[2] = (sprk_table.a[1] + sprk_table.a[2]) / 2.0;
    sprk_table.ahat[3] = (sprk_table.a[2] + sprk_table.a[3]) / 2.0;
    sprk_table.ahat[4] = sprk_table.ahat[3];
    sprk_table.ahat[5] = sprk_table.ahat[2];
    sprk_table.ahat[6] = sprk_table.ahat[1];
    sprk_table.ahat[7] = sprk_table.ahat[0];
    sprk_table
}

/*
  The following methods are from:

  (Original) Suzuki, M., & Umeno, K. (1993). Higher-order decomposition
  theory of exponential operators and its applications to QMC and
  nonlinear dynamics. Computer simulation studies in condensed-matter
  physics VI, 74-86. https://doi.org/10.1007/978-3-642-78448-4_7

  McLachlan, R.I.: On the Numerical Integration of Ordinary Differential
  Equations by Symmetric Composition Methods. Siam J Sci Comput. 16,
  151-168 (1995). https://doi.org/10.1137/0916010
 */

fn arkodeSymplecticSuzukiUmeno816() -> ARKodeSPRKTable {
    let mut sprk_table = ARKodeSPRKTable_Alloc(16);
    sprk_table.q = 8;
    sprk_table.stages = 16;
    sprk_table.a[0] = 0.7416703643506129534482278017838063156035;
    sprk_table.a[1] = -0.4091008258000315939973000958935634173099;
    sprk_table.a[2] = 0.1907547102962383799538762564503716627355;
    sprk_table.a[3] = -0.5738624711160822666563877266355357421595;
    sprk_table.a[4] = 0.2990641813036559238444635406886029882258;
    sprk_table.a[5] = 0.3346249182452981837849579798821822886337;
    sprk_table.a[6] = 0.3152930923967665966320566638110024309941;
    sprk_table.a[7] = -0.7968879393529163540197888401737330534463;
    sprk_table.a[8] = sprk_table.a[6];
    sprk_table.a[9] = sprk_table.a[5];
    sprk_table.a[10] = sprk_table.a[4];
    sprk_table.a[11] = sprk_table.a[3];
    sprk_table.a[12] = sprk_table.a[2];
    sprk_table.a[13] = sprk_table.a[1];
    sprk_table.a[14] = sprk_table.a[0];
    sprk_table.a[15] = 0.0;
    sprk_table.ahat[0] = sprk_table.a[0] / 2.0;
    sprk_table.ahat[1] = (sprk_table.a[0] + sprk_table.a[1]) / 2.0;
    sprk_table.ahat[2] = (sprk_table.a[1] + sprk_table.a[2]) / 2.0;
    sprk_table.ahat[3] = (sprk_table.a[2] + sprk_table.a[3]) / 2.0;
    sprk_table.ahat[4] = (sprk_table.a[3] + sprk_table.a[4]) / 2.0;
    sprk_table.ahat[5] = (sprk_table.a[4] + sprk_table.a[5]) / 2.0;
    sprk_table.ahat[6] = (sprk_table.a[5] + sprk_table.a[6]) / 2.0;
    sprk_table.ahat[7] = (sprk_table.a[6] + sprk_table.a[7]) / 2.0;
    sprk_table.ahat[8] = sprk_table.ahat[7];
    sprk_table.ahat[9] = sprk_table.ahat[6];
    sprk_table.ahat[10] = sprk_table.ahat[5];
    sprk_table.ahat[11] = sprk_table.ahat[4];
    sprk_table.ahat[12] = sprk_table.ahat[3];
    sprk_table.ahat[13] = sprk_table.ahat[2];
    sprk_table.ahat[14] = sprk_table.ahat[1];
    sprk_table.ahat[15] = sprk_table.ahat[0];
    sprk_table
}

/*
  The following methods are from:

  Sofroniou, M., Spaletta, G.: Derivation of symmetric composition
  constants for symmetric integrators. Optim Methods Softw. 20, 597-613
  (2005). https://doi.org/10.1080/10556780500140664
 */

fn arkodeSymplecticSofroniou10() -> ARKodeSPRKTable {
    let mut sprk_table = ARKodeSPRKTable_Alloc(36);
    sprk_table.q = 10;
    sprk_table.stages = 36;

    sprk_table.a[0] = 0.078795722521686419263907679337684;
    sprk_table.a[1] = 0.31309610341510852776481247192647;
    sprk_table.a[2] = 0.027918383235078066109520273275299;
    sprk_table.a[3] = -0.22959284159390709415121339679655;
    sprk_table.a[4] = 0.13096206107716486317465685927961;
    sprk_table.a[5] = -0.26973340565451071434460973222411;
    sprk_table.a[6] = 0.074973343155891435666137105641410;
    sprk_table.a[7] = 0.11199342399981020488957508073640;
    sprk_table.a[8] = 0.36613344954622675119314812353150;
    sprk_table.a[9] = -0.39910563013603589787862981058340;
    sprk_table.a[10] = 0.10308739852747107731580277001372;
    sprk_table.a[11] = 0.41143087395589023782070411897608;
    sprk_table.a[12] = -0.0048663605831352617621956593099771;
    sprk_table.a[13] = -0.39203335370863990644808193642610;
    sprk_table.a[14] = 0.051942502962449647037182904015976;
    sprk_table.a[15] = 0.050665090759924496335874344156866;
    sprk_table.a[16] = 0.049674370639729879054568800279461;
    sprk_table.a[17] = 0.049317735759594537917680008339338;
    for i in 18..35 {
        sprk_table.a[i] = sprk_table.a[34 - i];
    }
    sprk_table.a[35] = 0.0;
    sprk_table.ahat[0] = sprk_table.a[0] / 2.0;
    for i in 1..=17 {
        sprk_table.ahat[i] = (sprk_table.a[i - 1] + sprk_table.a[i]) / 2.0;
    }
    for i in 18..=35 {
        sprk_table.ahat[i] = sprk_table.ahat[35 - i];
    }

    sprk_table
}

pub fn ARKodeSPRKTable_Create(
    s: i32,
    q: i32,
    a: &[f64],
    ahat: &[f64],
) -> Option<ARKodeSPRKTable> {
    if s < 1 {
        return None;
    }

    let mut sprk_table = ARKodeSPRKTable_Alloc(s);

    sprk_table.stages = s;
    sprk_table.q = q;

    for i in 0..s as usize {
        sprk_table.a[i] = a[i];
        sprk_table.ahat[i] = ahat[i];
    }

    Some(sprk_table)
}

pub fn ARKodeSPRKTable_Alloc(stages: i32) -> ARKodeSPRKTable {
    ARKodeSPRKTable {
        q: 0,
        stages,
        a: vec![0.0; stages as usize],
        ahat: vec![0.0; stages as usize],
    }
}

pub fn ARKodeSPRKTable_Load(id: ARKODE_SPRKMethodID) -> Option<ARKodeSPRKTable> {
    match id {
        ARKODE_SPRK_EULER_1_1 => Some(arkodeSymplecticEuler()),
        ARKODE_SPRK_LEAPFROG_2_2 => Some(arkodeSymplecticLeapfrog2()),
        ARKODE_SPRK_PSEUDO_LEAPFROG_2_2 => Some(arkodeSymplecticPseudoLeapfrog2()),
        ARKODE_SPRK_RUTH_3_3 => Some(arkodeSymplecticRuth3()),
        ARKODE_SPRK_MCLACHLAN_2_2 => Some(arkodeSymplecticMcLachlan2()),
        ARKODE_SPRK_MCLACHLAN_3_3 => Some(arkodeSymplecticMcLachlan3()),
        ARKODE_SPRK_MCLACHLAN_4_4 => Some(arkodeSymplecticMcLachlan4()),
        ARKODE_SPRK_CANDY_ROZMUS_4_4 => Some(arkodeSymplecticCandyRozmus4()),
        ARKODE_SPRK_MCLACHLAN_5_6 => Some(arkodeSymplecticMcLachlan5()),
        ARKODE_SPRK_YOSHIDA_6_8 => Some(arkodeSymplecticYoshida6()),
        ARKODE_SPRK_SUZUKI_UMENO_8_16 => Some(arkodeSymplecticSuzukiUmeno816()),
        ARKODE_SPRK_SOFRONIOU_10_36 => Some(arkodeSymplecticSofroniou10()),
        _ => None,
    }
}

pub fn ARKodeSPRKTable_LoadByName(method: &str) -> Option<ARKodeSPRKTable> {
    match method {
        "ARKODE_SPRK_EULER_1_1" => Some(arkodeSymplecticEuler()),
        "ARKODE_SPRK_LEAPFROG_2_2" => Some(arkodeSymplecticLeapfrog2()),
        "ARKODE_SPRK_PSEUDO_LEAPFROG_2_2" => Some(arkodeSymplecticPseudoLeapfrog2()),
        "ARKODE_SPRK_RUTH_3_3" => Some(arkodeSymplecticRuth3()),
        "ARKODE_SPRK_MCLACHLAN_2_2" => Some(arkodeSymplecticMcLachlan2()),
        "ARKODE_SPRK_MCLACHLAN_3_3" => Some(arkodeSymplecticMcLachlan3()),
        "ARKODE_SPRK_MCLACHLAN_4_4" => Some(arkodeSymplecticMcLachlan4()),
        "ARKODE_SPRK_CANDY_ROZMUS_4_4" => Some(arkodeSymplecticCandyRozmus4()),
        "ARKODE_SPRK_MCLACHLAN_5_6" => Some(arkodeSymplecticMcLachlan5()),
        "ARKODE_SPRK_YOSHIDA_6_8" => Some(arkodeSymplecticYoshida6()),
        "ARKODE_SPRK_SUZUKI_UMENO_8_16" => Some(arkodeSymplecticSuzukiUmeno816()),
        "ARKODE_SPRK_SOFRONIOU_10_36" => Some(arkodeSymplecticSofroniou10()),
        _ => None,
    }
}

pub fn ARKodeSPRKTable_Copy(that_sprk_table: &ARKodeSPRKTable) -> ARKodeSPRKTable {
    let mut sprk_table = ARKodeSPRKTable_Alloc(that_sprk_table.stages);

    sprk_table.q = that_sprk_table.q;

    for i in 0..sprk_table.stages as usize {
        sprk_table.ahat[i] = that_sprk_table.ahat[i];
        sprk_table.a[i] = that_sprk_table.a[i];
    }

    sprk_table
}

pub fn ARKodeSPRKTable_Space(sprk_table: &ARKodeSPRKTable, liw: &mut i64, lrw: &mut i64) {
    *liw = 2;
    *lrw = sprk_table.stages as i64 * 2;
}

pub fn ARKodeSPRKTable_Write(sprk_table: &ARKodeSPRKTable, outfile: &mut dyn std::io::Write) {
    let mut a: Option<ARKodeButcherTable> = None;
    let mut b: Option<ARKodeButcherTable> = None;

    ARKodeSPRKTable_ToButcher(sprk_table, &mut a, &mut b);

    crate::arkode_butcher::ARKodeButcherTable_Write(a.as_ref().unwrap(), outfile);
    crate::arkode_butcher::ARKodeButcherTable_Write(b.as_ref().unwrap(), outfile);
}

pub fn ARKodeSPRKTable_ToButcher(
    sprk_table: &ARKodeSPRKTable,
    a_ptr: &mut Option<ARKodeButcherTable>,
    b_ptr: &mut Option<ARKodeButcherTable>,
) -> i32 {
    let stages = sprk_table.stages as usize;
    let mut a = ARKodeButcherTable_Alloc(sprk_table.stages, false).unwrap();
    let mut b = ARKodeButcherTable_Alloc(sprk_table.stages, false).unwrap();

    /* NOTE: the C routine's outer `for (i = 0; ...)` loop reuses `i` in
       its nested loops, so the outer body executes exactly once (with
       i == 0) and only b->b[0] / b->A[0][0] of the DIRK weights are
       populated; this quirk is replicated as straight-line code. */

    /* DIRK table (outer iteration i == 0 only) */
    b.b[0] = sprk_table.ahat[0];
    b.A[0][0] = sprk_table.ahat[0];

    /* Time weights: C_j = sum_{i=0}^{j} b_i */
    for j in 0..stages {
        for i in 0..=j {
            b.c[j] += sprk_table.ahat[i];
        }
    }

    /* Explicit table */
    for i in 0..stages {
        a.b[i] = sprk_table.a[i];
        for j in 0..i {
            a.A[i][j] = sprk_table.a[j];
        }
    }

    /* Time weights: c_j = sum_{i=0}^{j-1} a_i */
    for j in 0..stages {
        for i in 0..j {
            a.c[j] += sprk_table.a[i];
        }
    }

    /* Set method order */
    a.q = sprk_table.q;
    b.q = sprk_table.q;

    /* No embedding, so set embedding order to 0 */
    a.p = 0;
    b.p = 0;

    *a_ptr = Some(a);
    *b_ptr = Some(b);

    ARK_SUCCESS
}

use crate::arkode_impl::ARK_SUCCESS;

pub fn ARKodeSPRKTable_Free(sprk_table: ARKodeSPRKTable) {
    drop(sprk_table);
}
