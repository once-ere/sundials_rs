/* -----------------------------------------------------------------
 * cvodes_rs — pure-Rust, memory-safe translation of CVODES v7.7.0
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

// CVODES proper (modules land phase by phase; see ../../PROGRESS.md)
pub mod cvodes_impl;
pub mod cvodes_ls_impl;
pub mod cvodes_diag_impl;
pub mod cvodes_proj_impl;
pub mod cvodes_bandpre_impl;
pub mod cvodes_bbdpre_impl;
pub mod cvodes;
pub mod cvodes_io;
pub mod cvodes_ls;
pub mod cvodes_nls;
pub mod cvodes_nls_sim;
pub mod cvodes_nls_stg;
pub mod cvodes_nls_stg1;
pub mod cvodes_diag;
pub mod cvodes_proj;
pub mod cvodes_bandpre;
pub mod cvodes_bbdpre;
pub mod cvodes_resize;
pub mod cvodes_cli;
pub mod cvodea;
pub mod cvodea_io;

// Flat prelude so examples can `use cvodes_rs::*;` like a C `#include`.
pub use crate::cvodes_impl::*;
pub use crate::cvodes_ls_impl::*;
pub use crate::cvodes_diag_impl::*;
pub use crate::cvodes_proj_impl::*;
pub use crate::cvodes_bandpre_impl::*;
pub use crate::cvodes_bbdpre_impl::*;
pub use crate::cvodes::*;
pub use crate::cvodes_io::*;
pub use crate::cvodes_ls::*;
pub use crate::cvodes_nls::*;
pub use crate::cvodes_nls_sim::*;
pub use crate::cvodes_nls_stg::*;
pub use crate::cvodes_nls_stg1::*;
pub use crate::cvodes_diag::*;
pub use crate::cvodes_proj::*;
pub use crate::cvodes_bandpre::*;
pub use crate::cvodes_bbdpre::*;
pub use crate::cvodes_resize::*;
pub use crate::cvodes_cli::*;
pub use crate::cvodea::*;
pub use crate::cvodea_io::*;
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
pub use crate::sunnonlinsol_fixedpoint::*;
pub use crate::sunnonlinsol_newton::*;
