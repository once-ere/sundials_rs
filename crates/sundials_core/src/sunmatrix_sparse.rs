/* -----------------------------------------------------------------
 * Translated from src/sunmatrix/sparse/sunmatrix_sparse.c
 * (SUNDIALS/CVODE 7.7.0).
 *
 * The C SUNMatrixContent_Sparse (M, N, NNZ, NP, data, sparsetype,
 * indexvals, indexptrs + the rowvals/colptrs alias pointers) becomes
 * `SparseMatrix`; NP is derived (`np()`): N for CSC, M for CSR.
 * The ops-table entries (SUNMatZero_Sparse, ...) are free functions,
 * dispatched from the `SUNMatrix` enum in sundials_matrix.rs.
 * C `realloc` growth becomes `Vec::resize` (which preserves existing
 * entries exactly like realloc, zero-filling new ones).
 * -----------------------------------------------------------------*/
use crate::nvector_serial::NVector;
use crate::sundials_context::SUNContext;
use crate::sundials_errors::{SUNErrCode, SUN_ERR_ARG_DIMSMISMATCH, SUN_SUCCESS};
use crate::sundials_math::SUNRabs;
use crate::sundials_matrix::SUNMatrix;
use crate::sunmatrix_band::BandMatrix;
use crate::sunmatrix_dense::DenseMatrix;

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/* SUN_CSC_MAT / SUN_CSR_MAT (sunmatrix_sparse.h) */
pub const CSC_MAT: i32 = 0;
pub const CSR_MAT: i32 = 1;
pub const SUN_CSC_MAT: i32 = CSC_MAT;
pub const SUN_CSR_MAT: i32 = CSR_MAT;

/// Compressed-sparse matrix (CSC or CSR).
#[derive(Clone, Debug)]
pub struct SparseMatrix {
    pub m: i64,          /* rows (SM_ROWS_S)                 */
    pub n: i64,          /* columns (SM_COLUMNS_S)           */
    pub nnz: i64,        /* allocated nonzeros (SM_NNZ_S)    */
    pub sparsetype: i32, /* CSC_MAT or CSR_MAT               */
    pub indexvals: Vec<i64>, /* row (CSC) / column (CSR) indices, len nnz */
    pub indexptrs: Vec<i64>, /* start of each column (CSC) / row (CSR), len np+1 */
    pub data: Vec<f64>,      /* nonzero values, len nnz       */
}

impl SparseMatrix {
    /// Zero-filled sparse matrix with storage for `nnz` nonzeros.
    pub fn new(m: i64, n: i64, nnz: i64, sparsetype: i32) -> Self {
        let np = if sparsetype == CSC_MAT { n } else { m };
        SparseMatrix {
            m,
            n,
            nnz,
            sparsetype,
            indexvals: vec![0; nnz as usize],
            indexptrs: vec![0; (np + 1) as usize],
            data: vec![ZERO; nnz as usize],
        }
    }

    /// SM_NP_S: number of index pointers (N for CSC, M for CSR).
    #[inline]
    pub fn np(&self) -> i64 {
        if self.sparsetype == CSC_MAT {
            self.n
        } else {
            self.m
        }
    }
}

/// SUNSparseMatrix: create a new sparse matrix wrapped in the generic
/// `SUNMatrix` enum. Panics on illegal input (C returns NULL under
/// SUN_ERR_ARG_OUTOFRANGE).
pub fn SUNSparseMatrix(
    m: i64,
    n: i64,
    nnz: i64,
    sparsetype: i32,
    _sunctx: &SUNContext,
) -> SUNMatrix {
    assert!(m > 0 && n > 0, "SUNSparseMatrix: SUN_ERR_ARG_OUTOFRANGE");
    assert!(nnz >= 0, "SUNSparseMatrix: SUN_ERR_ARG_OUTOFRANGE");
    assert!(
        sparsetype == CSC_MAT || sparsetype == CSR_MAT,
        "SUNSparseMatrix: SUN_ERR_ARG_OUTOFRANGE"
    );
    SUNMatrix::Sparse(SparseMatrix::new(m, n, nnz, sparsetype))
}

/// SUNSparseFromDenseMatrix: create a new sparse matrix from an
/// existing dense matrix by copying all values with magnitude
/// strictly greater than `droptol`.
pub fn SUNSparseFromDenseMatrix(ad: &DenseMatrix, droptol: f64, sparsetype: i32) -> SUNMatrix {
    assert!(
        sparsetype == CSC_MAT || sparsetype == CSR_MAT,
        "SUNSparseFromDenseMatrix: SUN_ERR_ARG_OUTOFRANGE"
    );
    assert!(
        droptol >= ZERO,
        "SUNSparseFromDenseMatrix: SUN_ERR_ARG_OUTOFRANGE"
    );

    /* set size of new matrix */
    let m = ad.m;
    let n = ad.n;

    /* determine total number of nonzeros */
    let mut nnz: i64 = 0;
    for j in 0..n {
        for i in 0..m {
            nnz += (SUNRabs(ad.get(i, j)) > droptol) as i64;
        }
    }

    /* allocate sparse matrix */
    let mut asp = SparseMatrix::new(m, n, nnz, sparsetype);

    /* copy nonzeros from Ad into As, based on CSR/CSC type */
    let mut nnz: i64 = 0;
    if sparsetype == CSC_MAT {
        for j in 0..n {
            asp.indexptrs[j as usize] = nnz;
            for i in 0..m {
                if SUNRabs(ad.get(i, j)) > droptol {
                    asp.indexvals[nnz as usize] = i;
                    asp.data[nnz as usize] = ad.get(i, j);
                    nnz += 1;
                }
            }
        }
        asp.indexptrs[n as usize] = nnz;
    } else {
        /* CSR_MAT */
        for i in 0..m {
            asp.indexptrs[i as usize] = nnz;
            for j in 0..n {
                if SUNRabs(ad.get(i, j)) > droptol {
                    asp.indexvals[nnz as usize] = j;
                    asp.data[nnz as usize] = ad.get(i, j);
                    nnz += 1;
                }
            }
        }
        asp.indexptrs[m as usize] = nnz;
    }

    SUNMatrix::Sparse(asp)
}

/// SUNSparseFromBandMatrix: create a new sparse matrix from an
/// existing band matrix by copying all values with magnitude strictly
/// greater than `droptol`.
pub fn SUNSparseFromBandMatrix(ab: &BandMatrix, droptol: f64, sparsetype: i32) -> SUNMatrix {
    assert!(
        sparsetype == CSC_MAT || sparsetype == CSR_MAT,
        "SUNSparseFromBandMatrix: SUN_ERR_ARG_OUTOFRANGE"
    );
    assert!(
        droptol >= ZERO,
        "SUNSparseFromBandMatrix: SUN_ERR_ARG_OUTOFRANGE"
    );

    /* set size of new matrix */
    let m = ab.n;
    let n = ab.n;

    /* determine total number of nonzeros */
    let mut nnz: i64 = 0;
    for j in 0..n {
        let is = i64::max(0, j - ab.mu);
        let ie = i64::min(m - 1, j + ab.ml);
        for i in is..=ie {
            nnz += (SUNRabs(ab.get(i, j)) > droptol) as i64;
        }
    }

    /* allocate sparse matrix */
    let mut asp = SparseMatrix::new(m, n, nnz, sparsetype);

    /* copy nonzeros from Ab into As, based on CSR/CSC type */
    let mut nnz: i64 = 0;
    if sparsetype == CSC_MAT {
        for j in 0..n {
            asp.indexptrs[j as usize] = nnz;
            let is = i64::max(0, j - ab.mu);
            let ie = i64::min(m - 1, j + ab.ml);
            for i in is..=ie {
                if SUNRabs(ab.get(i, j)) > droptol {
                    asp.indexvals[nnz as usize] = i;
                    asp.data[nnz as usize] = ab.get(i, j);
                    nnz += 1;
                }
            }
        }
        asp.indexptrs[n as usize] = nnz;
    } else {
        /* CSR_MAT */
        for i in 0..m {
            asp.indexptrs[i as usize] = nnz;
            let js = i64::max(0, i - ab.ml);
            let je = i64::min(n - 1, i + ab.mu);
            for j in js..=je {
                if SUNRabs(ab.get(i, j)) > droptol {
                    asp.indexvals[nnz as usize] = j;
                    asp.data[nnz as usize] = ab.get(i, j);
                    nnz += 1;
                }
            }
        }
        asp.indexptrs[m as usize] = nnz;
    }

    SUNMatrix::Sparse(asp)
}

/// SUNSparseMatrix_ToCSR: create a new CSR matrix from a CSC matrix.
pub fn SUNSparseMatrix_ToCSR(a: &SparseMatrix) -> SparseMatrix {
    assert!(
        a.sparsetype == CSC_MAT,
        "SUNSparseMatrix_ToCSR: SUN_ERR_ARG_OUTOFRANGE"
    );
    let mut b = SparseMatrix::new(a.m, a.n, a.nnz, CSR_MAT);
    format_convert(a, &mut b);
    b
}

/// SUNSparseMatrix_ToCSC: create a new CSC matrix from a CSR matrix.
pub fn SUNSparseMatrix_ToCSC(a: &SparseMatrix) -> SparseMatrix {
    assert!(
        a.sparsetype == CSR_MAT,
        "SUNSparseMatrix_ToCSC: SUN_ERR_ARG_OUTOFRANGE"
    );
    let mut b = SparseMatrix::new(a.m, a.n, a.nnz, CSC_MAT);
    format_convert(a, &mut b);
    b
}

/// SUNSparseMatrix_Realloc: shrink (or grow) the storage arrays so the
/// matrix holds exactly indexptrs[NP] nonzeros.
pub fn SUNSparseMatrix_Realloc(a: &mut SparseMatrix) -> SUNErrCode {
    let nzmax = a.indexptrs[a.np() as usize];
    a.indexvals.resize(nzmax as usize, 0);
    a.data.resize(nzmax as usize, ZERO);
    a.nnz = nzmax;
    SUN_SUCCESS
}

/// SUNSparseMatrix_Reallocate: resize the storage arrays to hold a
/// specified number of nonzeros.
pub fn SUNSparseMatrix_Reallocate(a: &mut SparseMatrix, nnz: i64) -> SUNErrCode {
    a.indexvals.resize(nnz as usize, 0);
    a.data.resize(nnz as usize, ZERO);
    a.nnz = nnz;
    SUN_SUCCESS
}

/// SUNMatClone_Sparse: same shape / storage, zeroed.
pub fn SUNMatClone_Sparse(a: &SparseMatrix) -> SparseMatrix {
    SparseMatrix::new(a.m, a.n, a.nnz, a.sparsetype)
}

/// SUNMatZero_Sparse.
pub fn SUNMatZero_Sparse(a: &mut SparseMatrix) -> SUNErrCode {
    for i in 0..a.nnz as usize {
        a.data[i] = ZERO;
        a.indexvals[i] = 0;
    }
    let np = a.np() as usize;
    for i in 0..np {
        a.indexptrs[i] = 0;
    }
    a.indexptrs[np] = 0;
    SUN_SUCCESS
}

/// SUNMatCopy_Sparse: B = A (growing B's storage if needed).
pub fn SUNMatCopy_Sparse(a: &SparseMatrix, b: &mut SparseMatrix) -> SUNErrCode {
    /* both matrices must have the same shape and sparsity type */
    if a.m != b.m || a.n != b.n || a.sparsetype != b.sparsetype {
        return SUN_ERR_ARG_DIMSMISMATCH;
    }

    /* Perform operation */
    let a_nz = a.indexptrs[a.np() as usize];

    /* ensure that B is allocated with at least as
    much memory as we have nonzeros in A */
    if b.nnz < a_nz {
        b.indexvals.resize(a_nz as usize, 0);
        b.data.resize(a_nz as usize, ZERO);
        b.nnz = a_nz;
    }

    /* zero out B so that copy works correctly */
    SUNMatZero_Sparse(b);

    /* copy the data and row indices over */
    for i in 0..a_nz as usize {
        b.data[i] = a.data[i];
        b.indexvals[i] = a.indexvals[i];
    }

    /* copy the column pointers over */
    let np = a.np() as usize;
    for i in 0..np {
        b.indexptrs[i] = a.indexptrs[i];
    }
    b.indexptrs[np] = a_nz;

    SUN_SUCCESS
}

/// SUNMatScaleAddI_Sparse: A = c*A + I, inserting new nonzeros for any
/// missing diagonal entries (reallocating if required).
pub fn SUNMatScaleAddI_Sparse(c: f64, a: &mut SparseMatrix) -> SUNErrCode {
    let big_n: i64 = if a.sparsetype == CSC_MAT { a.n } else { a.m };
    let big_m: i64 = if a.sparsetype == CSC_MAT { a.m } else { a.n };

    let mut newvals: i64 = 0;
    for j in 0..big_n {
        /* scan column (row if CSR) of A, searching for diagonal value */
        let mut found = false;
        for i in a.indexptrs[j as usize]..a.indexptrs[(j + 1) as usize] {
            let iu = i as usize;
            if a.indexvals[iu] == j {
                found = true;
                a.data[iu] = ONE + c * a.data[iu];
            } else {
                a.data[iu] *= c;
            }
        }
        /* If no diagonal element found and the current column (row) can
         * actually contain a diagonal element, increment the counter */
        if !found && j < big_m {
            newvals += 1;
        }
    }

    /* At this point, A has the correctly updated values except for any
     * new diagonal elements that need to be added (of which there are
     * newvals). Now, we allocate additional storage if needed */
    let new_nnz = a.indexptrs[big_n as usize] + newvals;
    if new_nnz > a.nnz {
        SUNSparseMatrix_Reallocate(a, new_nnz);
    }

    let mut j = big_n - 1;
    while newvals > 0 {
        let mut found = false;
        let mut i = a.indexptrs[(j + 1) as usize] - 1;
        while i >= a.indexptrs[j as usize] {
            if a.indexvals[i as usize] == j {
                found = true;
            }

            /* Shift elements to make room for diagonal elements */
            a.indexvals[(i + newvals) as usize] = a.indexvals[i as usize];
            a.data[(i + newvals) as usize] = a.data[i as usize];
            i -= 1;
        }

        a.indexptrs[(j + 1) as usize] += newvals;
        if !found && j < big_m {
            /* This column (row) needs a diagonal element added */
            newvals -= 1;
            a.indexvals[(a.indexptrs[j as usize] + newvals) as usize] = j;
            a.data[(a.indexptrs[j as usize] + newvals) as usize] = ONE;
        }
        j -= 1;
    }

    SUN_SUCCESS
}

/// SUNMatScaleAdd_Sparse: A = c*A + B, expanding A's sparsity pattern
/// (and reallocating) as necessary.
pub fn SUNMatScaleAdd_Sparse(c: f64, a: &mut SparseMatrix, b: &SparseMatrix) -> SUNErrCode {
    if a.m != b.m || a.n != b.n || a.sparsetype != b.sparsetype {
        return SUN_ERR_ARG_DIMSMISMATCH;
    }

    /* store shortcuts to matrix dimensions (M is inner, N is outer) */
    let (big_m, big_n) = if a.sparsetype == CSC_MAT {
        (a.m, a.n)
    } else {
        (a.n, a.m)
    };

    /* create work arrays for row indices and nonzero column values */
    let mut w: Vec<i64> = vec![0; big_m as usize];
    let mut x: Vec<f64> = vec![ZERO; big_m as usize];

    /* determine if A already contains the sparsity pattern of B */
    let mut newvals: i64 = 0;
    for j in 0..big_n as usize {
        /* clear work array */
        for wi in w.iter_mut() {
            *wi = 0;
        }

        /* scan column of A, incrementing w by one */
        for i in a.indexptrs[j]..a.indexptrs[j + 1] {
            w[a.indexvals[i as usize] as usize] += 1;
        }

        /* scan column of B, decrementing w by one */
        for i in b.indexptrs[j]..b.indexptrs[j + 1] {
            w[b.indexvals[i as usize] as usize] -= 1;
        }

        /* if any entry of w is negative, A doesn't contain B's sparsity */
        for wi in w.iter() {
            if *wi < 0 {
                newvals += 1;
            }
        }
    }

    /* If extra nonzeros required, check whether A has sufficient
    storage space for new nonzero entries */
    let newmat = newvals > (a.nnz - a.indexptrs[big_n as usize]);

    /* perform operation based on existing/necessary structure */

    if newvals == 0 {
        /* case 1: A already contains sparsity pattern of B */
        for j in 0..big_n as usize {
            /* clear work array */
            for xi in x.iter_mut() {
                *xi = ZERO;
            }

            /* scan column of B, updating work array */
            for i in b.indexptrs[j]..b.indexptrs[j + 1] {
                x[b.indexvals[i as usize] as usize] = b.data[i as usize];
            }

            /* scan column of A, updating array entries appropriately */
            for i in a.indexptrs[j]..a.indexptrs[j + 1] {
                let iu = i as usize;
                a.data[iu] = c * a.data[iu] + x[a.indexvals[iu] as usize];
            }
        }
    } else if !newmat {
        /* case 2: A has sufficient storage, but does not already
         * contain B's sparsity */

        /* determine storage location where last column (row) should end */
        let mut nz = a.indexptrs[big_n as usize] + newvals;

        /* store pointer past last column (row) from original A,
        and store updated value in revised A */
        let mut cend = a.indexptrs[big_n as usize];
        a.indexptrs[big_n as usize] = nz;

        /* iterate through columns (rows) backwards */
        let mut j = big_n - 1;
        while j >= 0 {
            let ju = j as usize;

            /* clear out temporary arrays for this column (row) */
            for i in 0..big_m as usize {
                w[i] = 0;
                x[i] = ZERO;
            }

            /* iterate down column (row) of A, collecting nonzeros */
            let mut p = a.indexptrs[ju];
            while p < cend {
                let pu = p as usize;
                w[a.indexvals[pu] as usize] += 1;
                x[a.indexvals[pu] as usize] = c * a.data[pu];
                p += 1;
            }

            /* iterate down column of B, collecting nonzeros */
            for p in b.indexptrs[ju]..b.indexptrs[ju + 1] {
                let pu = p as usize;
                w[b.indexvals[pu] as usize] += 1;
                x[b.indexvals[pu] as usize] += b.data[pu];
            }

            /* fill entries of A with this column's (row's) data */
            let mut i = big_m - 1;
            while i >= 0 {
                if w[i as usize] > 0 {
                    nz -= 1;
                    a.indexvals[nz as usize] = i;
                    a.data[nz as usize] = x[i as usize];
                }
                i -= 1;
            }

            /* store ptr past this col (row) from orig A, update new A */
            cend = a.indexptrs[ju];
            a.indexptrs[ju] = nz;

            j -= 1;
        }
    } else {
        /* case 3: A must be reallocated with sufficient storage */

        /* create new matrix for sum */
        let mut cm = SparseMatrix::new(
            a.m,
            a.n,
            a.indexptrs[big_n as usize] + newvals,
            a.sparsetype,
        );

        /* initialize total nonzero count */
        let mut nz: i64 = 0;

        /* iterate through columns (rows) */
        for j in 0..big_n as usize {
            /* set current column (row) pointer to current # nonzeros */
            cm.indexptrs[j] = nz;

            /* clear out temporary arrays for this column (row) */
            for i in 0..big_m as usize {
                w[i] = 0;
                x[i] = ZERO;
            }

            /* iterate down column of A, collecting nonzeros */
            for p in a.indexptrs[j]..a.indexptrs[j + 1] {
                let pu = p as usize;
                w[a.indexvals[pu] as usize] += 1;
                x[a.indexvals[pu] as usize] = c * a.data[pu];
            }

            /* iterate down column of B, collecting nonzeros */
            for p in b.indexptrs[j]..b.indexptrs[j + 1] {
                let pu = p as usize;
                w[b.indexvals[pu] as usize] += 1;
                x[b.indexvals[pu] as usize] += b.data[pu];
            }

            /* fill entries of C with this column's data */
            for i in 0..big_m {
                if w[i as usize] > 0 {
                    cm.indexvals[nz as usize] = i;
                    cm.data[nz as usize] = x[i as usize];
                    nz += 1;
                }
            }
        }

        /* indicate end of data */
        cm.indexptrs[big_n as usize] = nz;

        /* update A's structure with C's values */
        a.nnz = cm.nnz;
        a.data = cm.data;
        a.indexvals = cm.indexvals;
        a.indexptrs = cm.indexptrs;
    }

    SUN_SUCCESS
}

/// SUNMatMatvec_Sparse: y = A*x.
pub fn SUNMatMatvec_Sparse(a: &SparseMatrix, x: &NVector, y: &mut NVector) -> SUNErrCode {
    if x.len() as i64 != a.n || y.len() as i64 != a.m {
        return SUN_ERR_ARG_DIMSMISMATCH;
    }

    if a.sparsetype == CSC_MAT {
        Matvec_SparseCSC(a, x, y)
    } else {
        Matvec_SparseCSR(a, x, y)
    }
}

/// SUNMatHermitianTransposeVec_Sparse: y = A^T x.
pub fn SUNMatHermitianTransposeVec_Sparse(
    a: &SparseMatrix,
    x: &NVector,
    y: &mut NVector,
) -> SUNErrCode {
    if x.len() as i64 != a.m || y.len() as i64 != a.n {
        return SUN_ERR_ARG_DIMSMISMATCH;
    }

    if a.sparsetype == CSC_MAT {
        MatTransposeVec_SparseCSC(a, x, y)
    } else {
        MatTransposeVec_SparseCSR(a, x, y)
    }
}

/// SUNMatSpace_Sparse: (lenrw, leniw).
pub fn SUNMatSpace_Sparse(a: &SparseMatrix) -> (i64, i64) {
    (a.nnz, 10 + a.np() + a.nnz)
}

/* -----------------------------------------------------------------
 * private functions
 * -----------------------------------------------------------------*/

/// Computes y = A*x for a CSC matrix.
fn Matvec_SparseCSC(a: &SparseMatrix, x: &NVector, y: &mut NVector) -> SUNErrCode {
    let xd = &x.data;
    let yd = &mut y.data;

    /* initialize result */
    for yi in yd.iter_mut().take(a.m as usize) {
        *yi = ZERO;
    }

    /* iterate through matrix columns */
    for j in 0..a.n as usize {
        /* iterate down column of A, performing product */
        for i in a.indexptrs[j]..a.indexptrs[j + 1] {
            let iu = i as usize;
            yd[a.indexvals[iu] as usize] += a.data[iu] * xd[j];
        }
    }
    SUN_SUCCESS
}

/// Computes y = A^T*x for a CSC matrix.
fn MatTransposeVec_SparseCSC(a: &SparseMatrix, x: &NVector, y: &mut NVector) -> SUNErrCode {
    let xd = &x.data;
    let yd = &mut y.data;

    /* initialize result vector */
    for yi in yd.iter_mut().take(a.n as usize) {
        *yi = ZERO;
    }

    /* iterate through matrix columns (rows of the transpose) */
    for j in 0..a.n as usize {
        for i in a.indexptrs[j]..a.indexptrs[j + 1] {
            let iu = i as usize;
            yd[j] += a.data[iu] * xd[a.indexvals[iu] as usize];
        }
    }
    SUN_SUCCESS
}

/// Computes y = A*x for a CSR matrix.
fn Matvec_SparseCSR(a: &SparseMatrix, x: &NVector, y: &mut NVector) -> SUNErrCode {
    let xd = &x.data;
    let yd = &mut y.data;

    /* initialize result */
    for yi in yd.iter_mut().take(a.m as usize) {
        *yi = ZERO;
    }

    /* iterate through matrix rows */
    for i in 0..a.m as usize {
        for j in a.indexptrs[i]..a.indexptrs[i + 1] {
            let ju = j as usize;
            yd[i] += a.data[ju] * xd[a.indexvals[ju] as usize];
        }
    }
    SUN_SUCCESS
}

/// Computes y = A^T*x for a CSR matrix.
fn MatTransposeVec_SparseCSR(a: &SparseMatrix, x: &NVector, y: &mut NVector) -> SUNErrCode {
    let xd = &x.data;
    let yd = &mut y.data;

    /* initialize result vector */
    for yi in yd.iter_mut().take(a.n as usize) {
        *yi = ZERO;
    }

    /* iterate over rows of A (columns of the transpose) */
    for i in 0..a.m as usize {
        for j in a.indexptrs[i]..a.indexptrs[i + 1] {
            let ju = j as usize;
            yd[a.indexvals[ju] as usize] += a.data[ju] * xd[i];
        }
    }
    SUN_SUCCESS
}

/// format_convert (private in C): copies A into B, where B is in the
/// opposite storage format of A.
fn format_convert(a: &SparseMatrix, b: &mut SparseMatrix) -> SUNErrCode {
    if a.sparsetype == b.sparsetype {
        return SUNMatCopy_Sparse(a, b);
    }

    let n_row = if a.sparsetype == CSR_MAT { a.m } else { a.n };
    let n_col = if a.sparsetype == CSR_MAT { a.n } else { a.m };

    let nnz = a.indexptrs[n_row as usize];

    SUNMatZero_Sparse(b);

    /* compute number of non-zero entries per column (if CSR) or per
     * row (if CSC) of A */
    for n in 0..nnz as usize {
        b.indexptrs[a.indexvals[n] as usize] += 1;
    }

    /* cumulative sum the nnz per column to get Bp[] */
    let mut csum: i64 = 0;
    for col in 0..n_col as usize {
        let temp = b.indexptrs[col];
        b.indexptrs[col] = csum;
        csum += temp;
    }
    b.indexptrs[n_col as usize] = nnz;

    for row in 0..n_row {
        for jj in a.indexptrs[row as usize]..a.indexptrs[(row + 1) as usize] {
            let jju = jj as usize;
            let col = a.indexvals[jju] as usize;
            let dest = b.indexptrs[col] as usize;

            b.indexvals[dest] = row;
            b.data[dest] = a.data[jju];

            b.indexptrs[col] += 1;
        }
    }

    let mut last: i64 = 0;
    for col in 0..=n_col as usize {
        let temp = b.indexptrs[col];
        b.indexptrs[col] = last;
        last = temp;
    }

    SUN_SUCCESS
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_dense_csc_layout_and_matvec() {
        let mut ad = DenseMatrix::new(3, 3);
        /* A = [1 0 2; 0 3 0; 4 0 5] */
        ad.set(0, 0, 1.0);
        ad.set(2, 0, 4.0);
        ad.set(1, 1, 3.0);
        ad.set(0, 2, 2.0);
        ad.set(2, 2, 5.0);

        let asp = match SUNSparseFromDenseMatrix(&ad, 0.0, CSC_MAT) {
            SUNMatrix::Sparse(s) => s,
            _ => unreachable!(),
        };
        assert_eq!(asp.nnz, 5);
        assert_eq!(asp.indexptrs, vec![0, 2, 3, 5]);
        assert_eq!(asp.indexvals, vec![0, 2, 1, 0, 2]);
        assert_eq!(asp.data, vec![1.0, 4.0, 3.0, 2.0, 5.0]);

        let x = NVector::from_slice(&[1.0, 2.0, 3.0]);
        let mut y = NVector::new(3);
        assert_eq!(SUNMatMatvec_Sparse(&asp, &x, &mut y), SUN_SUCCESS);
        for (yi, ei) in y.data.iter().zip([7.0, 6.0, 19.0].iter()) {
            assert!((yi - ei).abs() < 1e-15, "got {yi}, want {ei}");
        }

        /* CSR round-trip through format conversion */
        let acsr = match SUNSparseFromDenseMatrix(&ad, 0.0, CSR_MAT) {
            SUNMatrix::Sparse(s) => s,
            _ => unreachable!(),
        };
        let back = SUNSparseMatrix_ToCSC(&acsr);
        assert_eq!(back.indexptrs, asp.indexptrs);
        assert_eq!(back.indexvals, asp.indexvals);
        assert_eq!(back.data, asp.data);
    }

    #[test]
    fn scale_addi_inserts_missing_diagonal() {
        /* A = [0 1; 1 0] (CSC, no stored diagonal); 2*A + I = [1 2; 2 1] */
        let mut a = SparseMatrix::new(2, 2, 2, CSC_MAT);
        a.indexptrs = vec![0, 1, 2];
        a.indexvals = vec![1, 0];
        a.data = vec![1.0, 1.0];

        assert_eq!(SUNMatScaleAddI_Sparse(2.0, &mut a), SUN_SUCCESS);
        assert_eq!(a.nnz, 4);

        let mut y = NVector::new(2);
        let e0 = NVector::from_slice(&[1.0, 0.0]);
        assert_eq!(SUNMatMatvec_Sparse(&a, &e0, &mut y), SUN_SUCCESS);
        assert!((y.ith(0) - 1.0).abs() < 1e-15);
        assert!((y.ith(1) - 2.0).abs() < 1e-15);
        let e1 = NVector::from_slice(&[0.0, 1.0]);
        assert_eq!(SUNMatMatvec_Sparse(&a, &e1, &mut y), SUN_SUCCESS);
        assert!((y.ith(0) - 2.0).abs() < 1e-15);
        assert!((y.ith(1) - 1.0).abs() < 1e-15);
    }

    #[test]
    fn scale_add_expands_sparsity() {
        /* A diagonal, B off-diagonal: A = c*A + B needs new pattern */
        let mut a = SparseMatrix::new(2, 2, 2, CSC_MAT);
        a.indexptrs = vec![0, 1, 2];
        a.indexvals = vec![0, 1];
        a.data = vec![1.0, 2.0];

        let mut b = SparseMatrix::new(2, 2, 2, CSC_MAT);
        b.indexptrs = vec![0, 1, 2];
        b.indexvals = vec![1, 0];
        b.data = vec![3.0, 4.0];

        assert_eq!(SUNMatScaleAdd_Sparse(10.0, &mut a, &b), SUN_SUCCESS);
        /* expect A = [10 4; 3 20] */
        let mut y = NVector::new(2);
        let e0 = NVector::from_slice(&[1.0, 0.0]);
        SUNMatMatvec_Sparse(&a, &e0, &mut y);
        assert!((y.ith(0) - 10.0).abs() < 1e-15);
        assert!((y.ith(1) - 3.0).abs() < 1e-15);
        let e1 = NVector::from_slice(&[0.0, 1.0]);
        SUNMatMatvec_Sparse(&a, &e1, &mut y);
        assert!((y.ith(0) - 4.0).abs() < 1e-15);
        assert!((y.ith(1) - 20.0).abs() < 1e-15);

        /* case 2: same operation, but A already has spare storage so
         * the in-place backward insertion path runs */
        let mut a2 = SparseMatrix::new(2, 2, 4, CSC_MAT);
        a2.indexptrs = vec![0, 1, 2];
        a2.indexvals = vec![0, 1, 0, 0];
        a2.data = vec![1.0, 2.0, 0.0, 0.0];
        assert_eq!(SUNMatScaleAdd_Sparse(10.0, &mut a2, &b), SUN_SUCCESS);
        SUNMatMatvec_Sparse(&a2, &e0, &mut y);
        assert!((y.ith(0) - 10.0).abs() < 1e-15);
        assert!((y.ith(1) - 3.0).abs() < 1e-15);
        SUNMatMatvec_Sparse(&a2, &e1, &mut y);
        assert!((y.ith(0) - 4.0).abs() < 1e-15);
        assert!((y.ith(1) - 20.0).abs() < 1e-15);
    }

    #[test]
    fn from_band_matches_band_entries() {
        let mut ab = BandMatrix::new(4, 1, 1, 2);
        for j in 0..4i64 {
            ab.set(j, j, 4.0 + j as f64);
            if j > 0 {
                ab.set(j, j - 1, -1.0);
                ab.set(j - 1, j, 1.0);
            }
        }
        let asp = match SUNSparseFromBandMatrix(&ab, 0.0, CSC_MAT) {
            SUNMatrix::Sparse(s) => s,
            _ => unreachable!(),
        };
        assert_eq!(asp.nnz, 4 + 3 + 3);
        /* column 1 holds (0,1)=1, (1,1)=5, (2,1)=-1 */
        let start = asp.indexptrs[1] as usize;
        assert_eq!(&asp.indexvals[start..start + 3], &[0, 1, 2]);
        assert_eq!(&asp.data[start..start + 3], &[1.0, 5.0, -1.0]);
    }
}
