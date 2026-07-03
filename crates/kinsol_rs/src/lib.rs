/* -----------------------------------------------------------------
 * kinsol_rs — pure-Rust, memory-safe translation of KINSOL v7.7.0
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

// KINSOL proper (modules land phase by phase; see ../../PROGRESS.md)
pub mod kinsol_impl;
pub mod kinsol_ls_impl;
pub mod kinsol;
pub mod kinsol_aa;
pub mod kinsol_orth;
pub mod kinsol_io;
pub mod kinsol_ls;
pub mod kinsol_bbdpre_impl;
pub mod kinsol_bbdpre;
pub mod kinsol_cli;

// Flat prelude so examples can `use kinsol_rs::*;` like a C `#include`.
pub use crate::kinsol_impl::*;
pub use crate::kinsol_ls_impl::*;
pub use crate::kinsol::*;
pub use crate::kinsol_aa::*;
pub use crate::kinsol_orth::*;
pub use crate::kinsol_io::*;
pub use crate::kinsol_ls::*;
pub use crate::kinsol_bbdpre_impl::*;
pub use crate::kinsol_bbdpre::*;
pub use crate::kinsol_cli::*;
pub use crate::nvector_serial::*;
pub use crate::sundials_context::*;
pub use crate::sundials_errors::*;
pub use crate::sundials_linearsolver::*;
pub use crate::sundials_math::*;
pub use crate::sundials_matrix::*;
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
