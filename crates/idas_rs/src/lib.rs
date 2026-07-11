/* -----------------------------------------------------------------
 * idas_rs — pure-Rust, memory-safe translation of IDAS v7.7.0
 * (SUNDIALS v7.7.0). Original C sources:
 *   Lawrence Livermore National Security, Southern Methodist
 *   University, University of Maryland Baltimore County and the
 *   SUNDIALS contributors. SPDX-License-Identifier: BSD-3-Clause
 *
 * File map: each module keeps the base name of the C file it was
 * translated from (see ARCHITECTURE.md at the workspace root).
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(non_upper_case_globals)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::needless_range_loop)]
/* clippy 1.94 stylistic lints on deliberately C-faithful constructs
   (same rationale as the sundials_core crate-level allows):
   - unnecessary_unwrap: `if x.is_some()` + `unwrap()` mirrors the C
     `if (ptr)` guard followed by direct use, keeping bodies line-for-line
   - needless_borrow / explicit_auto_deref: the `&mut *guard` RefCell
     reborrow in the iterative-solve closures (donor pattern)
   - ptr_arg: `&mut Vec<NVector>` where an empty Vec plays the C NULL
     N_Vector* and the callee may resize (idaLsGetY / idaa interpolation)
   - field_reassign_with_default: Default::default() + field assignments
     mirrors the C memset(0) + explicit-assignment initialization style */
#![allow(clippy::unnecessary_unwrap)]
#![allow(clippy::needless_borrow)]
#![allow(clippy::explicit_auto_deref)]
#![allow(clippy::ptr_arg)]
#![allow(clippy::field_reassign_with_default)]
#![forbid(unsafe_code)]

// Shared SUNDIALS core (re-exported so donor `crate::<mod>` paths resolve)
pub use sundials_core::sundials_types;
pub use sundials_core::sundials_math;
pub use sundials_core::sundials_errors;
pub use sundials_core::sundials_context;
pub use sundials_core::sundials_utils;
pub use sundials_core::nvector_serial;
pub use sundials_core::sundials_dense;
pub use sundials_core::sundials_band;
pub use sundials_core::sunmatrix_dense;
pub use sundials_core::sunmatrix_band;
pub use sundials_core::sunmatrix_sparse;
pub use sundials_core::sundials_matrix;
pub use sundials_core::sundials_iterative;
pub use sundials_core::sundials_linearsolver;
pub use sundials_core::sunlinsol_dense;
pub use sundials_core::sunlinsol_band;
pub use sundials_core::sunlinsol_spgmr;
pub use sundials_core::sunlinsol_spfgmr;
pub use sundials_core::sunlinsol_spbcgs;
pub use sundials_core::sunlinsol_sptfqmr;
pub use sundials_core::sunlinsol_pcg;
pub use sundials_core::sundials_nonlinearsolver;
pub use sundials_core::sunnonlinsol_newton;
pub use sundials_core::sunnonlinsol_fixedpoint;

// IDAS proper (modules land phase by phase; see ../../PROGRESS.md)
pub mod idas_impl;
pub mod idas_ls;
pub mod idas_ls_impl;

// Flat prelude so examples can `use idas_rs::*;` like a C `#include`.
pub use crate::idas_impl::*;
pub use crate::idas_ls::*;
pub use crate::idas_ls_impl::*;
pub use crate::nvector_serial::*;
pub use crate::sundials_context::*;
pub use crate::sundials_errors::*;
pub use crate::sundials_linearsolver::*;
pub use crate::sundials_math::*;
pub use crate::sundials_matrix::*;
pub use crate::sundials_nonlinearsolver::*;
pub use crate::sundials_types::*;
pub use crate::sunlinsol_band::*;
pub use crate::sunlinsol_dense::*;
pub use crate::sunlinsol_pcg::*;
pub use crate::sunlinsol_spbcgs::*;
pub use crate::sunlinsol_spfgmr::*;
pub use crate::sunlinsol_spgmr::*;
pub use crate::sunlinsol_sptfqmr::*;
pub use crate::sunmatrix_band::*;
pub use crate::sunmatrix_dense::*;
pub use crate::sunmatrix_sparse::*;
