/* -----------------------------------------------------------------
 * Translated from src/sundials/sundials_direct.c and
 * include/sundials/sundials_direct.h (SUNDIALS 7.7.0): the legacy
 * SUNDlsMat generic dense/band matrix type with its constructors,
 * destructors and helpers.
 *
 * The C `struct _DlsMat` (type tag + M, N, ldim, mu, ml, s_mu, data,
 * ldata, cols) becomes the `SUNDlsMat` enum over the existing
 * `DenseMatrix` / `BandMatrix` structs (sunmatrix_dense.rs /
 * sunmatrix_band.rs), which already carry the same column-major
 * storage: dense element (i,j) at data[j*m + i] (SUNDLS_DENSE_ELEM),
 * band element (i,j) at data[j*ldim + s_mu + i - j]
 * (SUNDLS_BAND_ELEM). The header's accessor macros are the
 * `get`/`set` methods below.
 *
 * The lowercase `sunrealtype**` variants (SUNDlsMat_newDenseMat,
 * SUNDlsMat_newBandMat) return plain `Vec<Vec<f64>>` column arrays,
 * mirroring the C column-pointer-array access a[j][i].
 *
 * C returns NULL on non-positive dimensions; here that is `None`
 * (allocation failure itself aborts in safe Rust). C mallocs
 * uninitialized storage; here all allocations are zero-filled.
 *
 * The SUNDlsMat LU wrappers of sundials_dense.c / sundials_band.c
 * (SUNDlsMat_DenseGETRF, SUNDlsMat_BandGBTRF, ...) are NOT here: they
 * collapsed into the lowercase kernels of sundials_dense.rs and
 * sundials_band.rs, which operate on DenseMatrix/BandMatrix directly.
 * -----------------------------------------------------------------*/
use crate::sundials_types::sunindextype;
use crate::sundials_utils::fmt_e;
use crate::sunmatrix_band::BandMatrix;
use crate::sunmatrix_dense::DenseMatrix;

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/// SUNDIALS_DENSE: dense matrix type tag.
pub const SUNDIALS_DENSE: i32 = 1;
/// SUNDIALS_BAND: banded matrix type tag.
pub const SUNDIALS_BAND: i32 = 2;

/// The legacy generic direct-solver matrix (C `SUNDlsMat`, a tagged
/// struct); the `type` field becomes the enum discriminant.
#[derive(Clone, Debug)]
pub enum SUNDlsMat {
    Dense(DenseMatrix),
    Band(BandMatrix),
}

impl SUNDlsMat {
    /// The C `A->type` field (SUNDIALS_DENSE or SUNDIALS_BAND).
    #[inline]
    pub fn mat_type(&self) -> i32 {
        match self {
            SUNDlsMat::Dense(_) => SUNDIALS_DENSE,
            SUNDlsMat::Band(_) => SUNDIALS_BAND,
        }
    }

    /// SUNDLS_DENSE_ELEM(A,i,j) / SUNDLS_BAND_ELEM(A,i,j).
    /// For band matrices (i,j) must satisfy j-mu <= i <= j+ml
    /// (up to j-s_mu during factorizations).
    #[inline]
    pub fn get(&self, i: sunindextype, j: sunindextype) -> f64 {
        match self {
            SUNDlsMat::Dense(a) => a.get(i, j),
            SUNDlsMat::Band(a) => a.get(i, j),
        }
    }

    /// SUNDLS_DENSE_ELEM(A,i,j) = v / SUNDLS_BAND_ELEM(A,i,j) = v.
    #[inline]
    pub fn set(&mut self, i: sunindextype, j: sunindextype, v: f64) {
        match self {
            SUNDlsMat::Dense(a) => a.set(i, j, v),
            SUNDlsMat::Band(a) => a.set(i, j, v),
        }
    }
}

/// SUNDlsMat_NewDenseMat: allocate an M-by-N dense SUNDlsMat
/// (C returns NULL for non-positive dimensions -> None).
pub fn SUNDlsMat_NewDenseMat(m: sunindextype, n: sunindextype) -> Option<SUNDlsMat> {
    if m <= 0 || n <= 0 {
        return None;
    }
    Some(SUNDlsMat::Dense(DenseMatrix::new(m, n)))
}

/// SUNDlsMat_newDenseMat: allocate an m-by-n column array
/// (`sunrealtype**` in C): n columns, each of length m.
pub fn SUNDlsMat_newDenseMat(m: sunindextype, n: sunindextype) -> Option<Vec<Vec<f64>>> {
    if n <= 0 || m <= 0 {
        return None;
    }
    Some(vec![vec![ZERO; m as usize]; n as usize])
}

/// SUNDlsMat_NewBandMat: allocate an N-by-N band SUNDlsMat with upper
/// bandwidth mu, lower bandwidth ml and storage upper bandwidth smu
/// (pass smu = mu if A will not be factored, smu = min(N-1, mu+ml)
/// if it will). Column j has leading dimension smu+ml+1 and its
/// diagonal element sits at offset s_mu (element (i,j) at offset
/// s_mu + i - j).
pub fn SUNDlsMat_NewBandMat(
    n: sunindextype,
    mu: sunindextype,
    ml: sunindextype,
    smu: sunindextype,
) -> Option<SUNDlsMat> {
    if n <= 0 {
        return None;
    }
    Some(SUNDlsMat::Band(BandMatrix::new(n, mu, ml, smu)))
}

/// SUNDlsMat_newBandMat: allocate the raw band column array
/// (`sunrealtype**` in C): n columns, each of length smu+ml+1.
pub fn SUNDlsMat_newBandMat(
    n: sunindextype,
    smu: sunindextype,
    ml: sunindextype,
) -> Option<Vec<Vec<f64>>> {
    if n <= 0 {
        return None;
    }
    let col_size = smu + ml + 1;
    Some(vec![vec![ZERO; col_size as usize]; n as usize])
}

/// SUNDlsMat_DestroyMat: free a SUNDlsMat (ownership drop in Rust).
pub fn SUNDlsMat_DestroyMat(_a: SUNDlsMat) {}

/// SUNDlsMat_destroyMat: free a column array (ownership drop).
pub fn SUNDlsMat_destroyMat(_a: Vec<Vec<f64>>) {}

/// SUNDlsMat_NewIntArray: allocate an array of N ints
/// (None for N <= 0, as C returns NULL).
pub fn SUNDlsMat_NewIntArray(n: i32) -> Option<Vec<i32>> {
    if n <= 0 {
        return None;
    }
    Some(vec![0; n as usize])
}

/// SUNDlsMat_newIntArray: allocate an array of n ints.
pub fn SUNDlsMat_newIntArray(n: i32) -> Option<Vec<i32>> {
    if n <= 0 {
        return None;
    }
    Some(vec![0; n as usize])
}

/// SUNDlsMat_NewIndexArray: allocate an array of N sunindextype.
pub fn SUNDlsMat_NewIndexArray(n: sunindextype) -> Option<Vec<sunindextype>> {
    if n <= 0 {
        return None;
    }
    Some(vec![0; n as usize])
}

/// SUNDlsMat_newIndexArray: allocate an array of n sunindextype.
pub fn SUNDlsMat_newIndexArray(n: sunindextype) -> Option<Vec<sunindextype>> {
    if n <= 0 {
        return None;
    }
    Some(vec![0; n as usize])
}

/// SUNDlsMat_NewRealArray: allocate an array of N sunrealtype.
pub fn SUNDlsMat_NewRealArray(n: sunindextype) -> Option<Vec<f64>> {
    if n <= 0 {
        return None;
    }
    Some(vec![ZERO; n as usize])
}

/// SUNDlsMat_newRealArray: allocate an array of m sunrealtype.
pub fn SUNDlsMat_newRealArray(m: sunindextype) -> Option<Vec<f64>> {
    if m <= 0 {
        return None;
    }
    Some(vec![ZERO; m as usize])
}

/// SUNDlsMat_DestroyArray: free an array allocated by
/// SUNDlsMat_New{Int,Index,Real}Array (ownership drop).
pub fn SUNDlsMat_DestroyArray<T>(_v: Vec<T>) {}

/// SUNDlsMat_destroyArray: free an array (ownership drop).
pub fn SUNDlsMat_destroyArray<T>(_v: Vec<T>) {}

/// SUNDlsMat_AddIdentity: A_ii += 1 on the leading diagonal
/// (i = 0..N-1 dense, i = 0..M-1 band; M == N for band matrices).
pub fn SUNDlsMat_AddIdentity(a: &mut SUNDlsMat) {
    match a {
        SUNDlsMat::Dense(a) => {
            /* case SUNDIALS_DENSE: A->cols[i][i] += ONE */
            let m = a.m as usize;
            for i in 0..a.n as usize {
                a.data[i * m + i] += ONE;
            }
        }
        SUNDlsMat::Band(a) => {
            /* case SUNDIALS_BAND: A->cols[i][A->s_mu] += ONE */
            let ldim = a.ldim as usize;
            let smu = a.s_mu as usize;
            for i in 0..a.n as usize {
                a.data[i * ldim + smu] += ONE;
            }
        }
    }
}

/// SUNDlsMat_SetToZero: zero all elements of A. For band matrices
/// only the declared band (mu+ml+1 entries per column, starting at
/// offset s_mu - mu) is cleared, exactly as in the C source — the
/// factorization fill-in rows above the band are left untouched.
pub fn SUNDlsMat_SetToZero(a: &mut SUNDlsMat) {
    match a {
        SUNDlsMat::Dense(a) => {
            for j in 0..a.n as usize {
                let cj = j * a.m as usize;
                for i in 0..a.m as usize {
                    a.data[cj + i] = ZERO;
                }
            }
        }
        SUNDlsMat::Band(a) => {
            let col_size = (a.mu + a.ml + 1) as usize;
            for j in 0..a.n as usize {
                let off = j * a.ldim as usize + (a.s_mu - a.mu) as usize;
                for i in 0..col_size {
                    a.data[off + i] = ZERO;
                }
            }
        }
    }
}

/// C `SUN_FORMAT_E` for double: `"% .15e"` (space flag, DBL_DIG = 15).
fn fmt_sun_e(x: f64) -> String {
    let s = fmt_e(x, 0, 15);
    if s.starts_with('-') {
        s
    } else {
        format!(" {}", s)
    }
}

/// SUNDlsMat_PrintMat: print the M-by-N (dense or band) matrix A to
/// `outfile` as it would normally appear on paper (debugging aid; the
/// C FILE* becomes any `std::io::Write`, write errors are ignored as
/// fprintf's return value is in C).
pub fn SUNDlsMat_PrintMat(a: &SUNDlsMat, outfile: &mut dyn std::io::Write) {
    match a {
        SUNDlsMat::Dense(a) => {
            let _ = writeln!(outfile);
            for i in 0..a.m {
                for j in 0..a.n {
                    let _ = write!(outfile, "{}  ", fmt_sun_e(a.get(i, j)));
                }
                let _ = writeln!(outfile);
            }
        }
        SUNDlsMat::Band(a) => {
            let _ = writeln!(outfile);
            for i in 0..a.n {
                let start = sunindextype::max(0, i - a.ml);
                let finish = sunindextype::min(a.n - 1, i + a.mu);
                for _ in 0..start {
                    /* fprintf(outfile, "%12s  ", "") */
                    let _ = write!(outfile, "{:12}  ", "");
                }
                for j in start..=finish {
                    /* a[j][i - j + A->s_mu] */
                    let _ = write!(outfile, "{}  ", fmt_sun_e(a.get(i, j)));
                }
                let _ = writeln!(outfile);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_dense_mat_shape_and_zero_fill() {
        let a = SUNDlsMat_NewDenseMat(3, 2).expect("valid dims");
        assert_eq!(a.mat_type(), SUNDIALS_DENSE);
        match &a {
            SUNDlsMat::Dense(d) => {
                assert_eq!(d.m, 3);
                assert_eq!(d.n, 2);
                assert_eq!(d.data.len(), 6);
                assert!(d.data.iter().all(|&v| v == 0.0));
            }
            SUNDlsMat::Band(_) => panic!("expected dense"),
        }
        /* C returns NULL on non-positive dimensions */
        assert!(SUNDlsMat_NewDenseMat(0, 2).is_none());
        assert!(SUNDlsMat_NewDenseMat(3, -1).is_none());
        SUNDlsMat_DestroyMat(a);
    }

    #[test]
    fn new_band_mat_shape_and_zero_fill() {
        let a = SUNDlsMat_NewBandMat(5, 1, 2, 3).expect("valid dims");
        assert_eq!(a.mat_type(), SUNDIALS_BAND);
        match &a {
            SUNDlsMat::Band(b) => {
                assert_eq!(b.n, 5);
                assert_eq!(b.mu, 1);
                assert_eq!(b.ml, 2);
                assert_eq!(b.s_mu, 3);
                assert_eq!(b.ldim, 6); /* colSize = smu + ml + 1 */
                assert_eq!(b.data.len(), 30); /* ldata = N * colSize */
                assert!(b.data.iter().all(|&v| v == 0.0));
            }
            SUNDlsMat::Dense(_) => panic!("expected band"),
        }
        assert!(SUNDlsMat_NewBandMat(0, 1, 1, 2).is_none());
        SUNDlsMat_DestroyMat(a);
    }

    #[test]
    fn lowercase_column_array_constructors() {
        let a = SUNDlsMat_newDenseMat(3, 2).expect("valid dims");
        assert_eq!(a.len(), 2); /* n column pointers */
        assert!(a.iter().all(|col| col.len() == 3)); /* each of length m */
        assert!(SUNDlsMat_newDenseMat(0, 2).is_none());
        assert!(SUNDlsMat_newDenseMat(3, 0).is_none());
        SUNDlsMat_destroyMat(a);

        let b = SUNDlsMat_newBandMat(4, 2, 1).expect("valid dims");
        assert_eq!(b.len(), 4); /* n columns */
        assert!(b.iter().all(|col| col.len() == 4)); /* smu + ml + 1 */
        assert!(SUNDlsMat_newBandMat(-1, 2, 1).is_none());
        SUNDlsMat_destroyMat(b);
    }

    #[test]
    fn array_allocators() {
        assert_eq!(SUNDlsMat_NewIntArray(3).expect("n > 0"), vec![0i32; 3]);
        assert!(SUNDlsMat_NewIntArray(0).is_none());
        assert_eq!(SUNDlsMat_newIntArray(2).expect("n > 0"), vec![0i32; 2]);
        assert!(SUNDlsMat_newIntArray(-4).is_none());

        assert_eq!(SUNDlsMat_NewIndexArray(4).expect("n > 0"), vec![0i64; 4]);
        assert!(SUNDlsMat_NewIndexArray(0).is_none());
        assert_eq!(SUNDlsMat_newIndexArray(1).expect("n > 0"), vec![0i64; 1]);
        assert!(SUNDlsMat_newIndexArray(0).is_none());

        assert_eq!(SUNDlsMat_NewRealArray(2).expect("n > 0"), vec![0.0; 2]);
        assert!(SUNDlsMat_NewRealArray(-1).is_none());
        assert_eq!(SUNDlsMat_newRealArray(5).expect("m > 0"), vec![0.0; 5]);
        assert!(SUNDlsMat_newRealArray(0).is_none());

        SUNDlsMat_DestroyArray(SUNDlsMat_NewIntArray(1).expect("n > 0"));
        SUNDlsMat_destroyArray(SUNDlsMat_newRealArray(1).expect("m > 0"));
    }

    #[test]
    fn dense_elem_roundtrip_and_add_identity() {
        let mut a = SUNDlsMat_NewDenseMat(3, 3).expect("valid dims");
        a.set(0, 0, 2.0);
        a.set(2, 1, -7.5); /* off-diagonal */
        a.set(1, 2, 4.25);
        assert_eq!(a.get(0, 0), 2.0);
        assert_eq!(a.get(2, 1), -7.5);
        assert_eq!(a.get(1, 2), 4.25);

        SUNDlsMat_AddIdentity(&mut a);
        assert_eq!(a.get(0, 0), 3.0); /* 2 + 1 */
        assert_eq!(a.get(1, 1), 1.0); /* 0 + 1 */
        assert_eq!(a.get(2, 2), 1.0);
        assert_eq!(a.get(2, 1), -7.5); /* off-diagonal untouched */
    }

    #[test]
    fn band_elem_roundtrip_offset_edges() {
        /* n = 5, mu = 1, ml = 2, smu = min(n-1, mu+ml) = 3 */
        let mut a = SUNDlsMat_NewBandMat(5, 1, 2, 3).expect("valid dims");

        /* top of band in column j: i = j - mu (storage offset s_mu - mu) */
        a.set(1, 2, 12.0);
        /* diagonal (storage offset s_mu) */
        a.set(2, 2, 22.0);
        /* bottom of band: i = j + ml (storage offset s_mu + ml) */
        a.set(4, 2, 42.0);
        /* corner columns */
        a.set(0, 1, 1.5); /* (j-mu, j) at j = 1 */
        a.set(4, 4, 9.0); /* last diagonal */

        assert_eq!(a.get(1, 2), 12.0);
        assert_eq!(a.get(2, 2), 22.0);
        assert_eq!(a.get(4, 2), 42.0);
        assert_eq!(a.get(0, 1), 1.5);
        assert_eq!(a.get(4, 4), 9.0);

        /* verify raw storage location: element (i,j) at
         * data[j*ldim + s_mu + i - j] */
        match &a {
            SUNDlsMat::Band(b) => {
                assert_eq!(b.data[(2 * b.ldim + b.s_mu + 1 - 2) as usize], 12.0);
                assert_eq!(b.data[(2 * b.ldim + b.s_mu) as usize], 22.0);
                assert_eq!(b.data[(2 * b.ldim + b.s_mu + 4 - 2) as usize], 42.0);
            }
            SUNDlsMat::Dense(_) => panic!("expected band"),
        }

        SUNDlsMat_AddIdentity(&mut a);
        assert_eq!(a.get(2, 2), 23.0);
        assert_eq!(a.get(0, 0), 1.0);
        assert_eq!(a.get(4, 4), 10.0);
        assert_eq!(a.get(1, 2), 12.0); /* off-diagonal untouched */
    }

    #[test]
    fn set_to_zero_dense_and_band() {
        let mut d = SUNDlsMat_NewDenseMat(2, 2).expect("valid dims");
        d.set(0, 0, 1.0);
        d.set(1, 1, 2.0);
        SUNDlsMat_SetToZero(&mut d);
        match &d {
            SUNDlsMat::Dense(m) => assert!(m.data.iter().all(|&v| v == 0.0)),
            SUNDlsMat::Band(_) => panic!("expected dense"),
        }

        /* band: only the declared band (s_mu - mu .. s_mu + ml) is
         * cleared; the fill-in rows above it stay untouched (C loops
         * over colSize = mu + ml + 1 from cols[j] + s_mu - mu). */
        let mut b = SUNDlsMat_NewBandMat(4, 1, 1, 2).expect("valid dims");
        b.set(0, 1, 5.0); /* in band */
        b.set(2, 2, 6.0); /* diagonal */
        b.set(3, 2, 7.0); /* in band */
        match &mut b {
            SUNDlsMat::Band(m) => m.data[2 * m.ldim as usize] = 99.0, /* fill-in row of column 2 */
            SUNDlsMat::Dense(_) => panic!("expected band"),
        }
        SUNDlsMat_SetToZero(&mut b);
        assert_eq!(b.get(0, 1), 0.0);
        assert_eq!(b.get(2, 2), 0.0);
        assert_eq!(b.get(3, 2), 0.0);
        match &b {
            SUNDlsMat::Band(m) => {
                assert_eq!(m.data[2 * m.ldim as usize], 99.0); /* fill-in row untouched */
            }
            SUNDlsMat::Dense(_) => panic!("expected band"),
        }
    }

    #[test]
    fn print_mat_formats() {
        let mut d = SUNDlsMat_NewDenseMat(2, 2).expect("valid dims");
        d.set(0, 0, 1.0);
        d.set(0, 1, -2.5);
        d.set(1, 0, 0.0);
        d.set(1, 1, 4.0);
        let mut out = Vec::new();
        SUNDlsMat_PrintMat(&d, &mut out);
        let s = String::from_utf8(out).expect("utf8");
        assert_eq!(
            s,
            "\n 1.000000000000000e+00  -2.500000000000000e+00  \n\
             \u{20}0.000000000000000e+00   4.000000000000000e+00  \n"
        );

        /* band: rows past the first get a 12-blank (+2) leading pad
         * per skipped column */
        let mut b = SUNDlsMat_NewBandMat(3, 1, 1, 2).expect("valid dims");
        for j in 0..3i64 {
            b.set(j, j, 2.0);
            if j > 0 {
                b.set(j, j - 1, -1.0);
                b.set(j - 1, j, 1.0);
            }
        }
        let mut out = Vec::new();
        SUNDlsMat_PrintMat(&b, &mut out);
        let s = String::from_utf8(out).expect("utf8");
        let lines: Vec<&str> = s.split('\n').collect();
        assert_eq!(lines[0], ""); /* leading blank line */
        /* row 0: columns 0..=1 */
        assert_eq!(
            lines[1],
            " 2.000000000000000e+00   1.000000000000000e+00  "
        );
        /* row 2: column 0 skipped -> "%12s  " pad (14 blanks), then
         * columns 1..=2 */
        assert_eq!(
            lines[3],
            "              -1.000000000000000e+00   2.000000000000000e+00  "
        );
    }
}
