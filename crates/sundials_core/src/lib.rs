/* -----------------------------------------------------------------
 * sundials_core — pure-Rust, memory-safe translation of the shared
 * SUNDIALS v7.7.0 core. Original C sources:
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

// SUNDIALS support layer
pub mod sundials_types;
pub mod sundials_math;
pub mod sundials_errors;
pub mod sundials_context;
pub mod sundials_utils;
pub mod nvector_serial;
pub mod sundials_dense;
pub mod sundials_band;
pub mod sunmatrix_dense;
pub mod sunmatrix_band;
pub mod sunmatrix_sparse;
pub mod sundials_matrix;
pub mod sundials_iterative;
pub mod sundials_linearsolver;
pub mod sunlinsol_dense;
pub mod sunlinsol_band;
pub mod sunlinsol_spgmr;
pub mod sunlinsol_spfgmr;
pub mod sunlinsol_spbcgs;
pub mod sunlinsol_sptfqmr;
pub mod sunlinsol_pcg;
pub mod sundials_nonlinearsolver;
pub mod sunnonlinsol_newton;
pub mod sunnonlinsol_fixedpoint;

// Phase 1 — remaining shared core (stubs being filled in)
pub mod sundials_version;
pub mod sundials_futils;
pub mod sundials_memory;
pub mod sundials_system_memory;
pub mod sundials_direct;
pub mod sundials_hashmap;
pub mod sundials_logger;
pub mod sundials_profiler;
pub mod sundials_cli;
pub mod sundials_adaptcontroller;
pub mod sunadaptcontroller_soderlind;
pub mod sunadaptcontroller_imexgus;
pub mod sunadaptcontroller_mrihtol;
pub mod sundials_domeigestimator;
pub mod sundomeigest_power;
pub mod sundomeigest_arnoldi;
pub mod sundials_nvector_senswrapper;
