/* -----------------------------------------------------------------
 * arkode_rs — pure-Rust, memory-safe translation of ARKODE v7.7.0
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
/* C-faithful constructs (line-for-line translation rule, CLAUDE.md):
   fprintf strings kept byte-identical to C, manual element loops
   mirroring C loops. Same rationale as the sundials_core allows. */
#![allow(clippy::write_with_newline)]
#![allow(clippy::manual_memcpy)]
/* C `(x > TOL) ? FALSE : TRUE` becomes !(x > TOL) so NaN keeps the C
   comparison semantics (same allow as sundials_core). */
#![allow(clippy::neg_cmp_op_on_partial_ord)]
/* C `return B;` inside the .def X-macro arms stays a return statement */
#![allow(clippy::needless_return)]
/* Butcher-table coefficients keep the C source's full decimal text */
#![allow(clippy::excessive_precision)]
/* C-shaped `x < lo || x > hi` bounds checks kept verbatim */
#![allow(clippy::manual_range_contains)]
/* C's nested if structure is kept as written */
#![allow(clippy::collapsible_if)]
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
pub use sundials_core::sundials_adaptcontroller;
pub use sundials_core::sunadaptcontroller_soderlind;
pub use sundials_core::sunadaptcontroller_imexgus;
pub use sundials_core::sunadaptcontroller_mrihtol;
pub use sundials_core::sundials_domeigestimator;
pub use sundials_core::sundomeigest_power;
pub use sundials_core::sundomeigest_arnoldi;
pub use sundials_core::sundials_stepper;
pub use sundials_core::sundials_datanode;
pub use sundials_core::sundatanode_inmem;
pub use sundials_core::sundials_adjointcheckpointscheme;
pub use sundials_core::sunadjointcheckpointscheme_fixed;
pub use sundials_core::sundials_adjointstepper;
pub use sundials_core::sundials_memory;
pub use sundials_core::sundials_system_memory;
pub use sundials_core::sundials_cli;

// ARKODE proper (modules land phase by phase; see ../../PROGRESS.md)
pub mod arkode_impl;
pub mod arkode_adapt_impl;
pub mod arkode_root_impl;
pub mod arkode_relaxation_impl;
pub mod arkode_interp_impl;
pub mod arkode_butcher;
pub mod arkode_butcher_erk;
pub mod arkode_butcher_dirk;
pub mod arkode;
pub mod arkode_interp;
pub mod arkode_adapt;
pub mod arkode_root;
pub mod arkode_io;
pub mod arkode_ls_impl;
pub mod arkode_ls;
pub mod arkode_arkstep_impl;
pub mod arkode_arkstep;
pub mod arkode_arkstep_nls;
pub mod arkode_arkstep_io;
pub mod arkode_sprk;
pub mod arkode_sprkstep_impl;
pub mod arkode_sprkstep;
pub mod arkode_sprkstep_io;
pub mod arkode_erkstep_impl;
pub mod arkode_erkstep;
pub mod arkode_erkstep_io;

// Flat prelude so examples can `use arkode_rs::*;` like a C `#include`.
pub use crate::arkode_impl::*;
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
