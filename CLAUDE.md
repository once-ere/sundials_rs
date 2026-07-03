# CLAUDE.md — sundials_rs workspace

Pure-Rust, memory-safe translation of SUNDIALS v7.7.0. No `unsafe`, no FFI,
no external crate dependencies, no warnings. Reference C sources live in
`../sundials-7.7.0/` (and `../cvode-7.7.0/`) — **never modify those trees**.
The original verified CVODE port `../cvode_rs/` is a read-only donor.

## Commands

- Build everything: `cargo build --workspace 2>&1 | tee /tmp/build.log`
- Tests: `cargo test --workspace`
- Verify examples: `tools/verify_examples.sh [crate|all]`, then Read
  `logs/summary.txt`; per-example diffs land in `logs/<name>.diff`.
- Run one example: `cargo run -p <crate> --example <name>`

## Layout

- Workspace crates in `crates/`: `sundials_core` (shared core) + one crate
  per solver (`cvode_rs`, `cvodes_rs`, `kinsol_rs`, `ida_rs`, `idas_rs`,
  `arkode_rs`). Module = original C file base name + `.rs`.
- Solver crates re-export every shared module from `sundials_core` at the
  crate root (`pub use sundials_core::nvector_serial;` …) so module bodies
  can use `crate::<mod>` paths uniformly; each also glob-re-exports a flat
  prelude so examples can `use <crate>::*;` like a C `#include`.
- `ARCHITECTURE.md` — pinned cross-module contracts; read before touching
  cross-module types. `PROGRESS.md` — per-file port status. `VERIFICATION.md`
  — per-example verification matrix. Keep both current.

## Hard rules

1. **Fidelity first**: control flow, constants, tolerances, heuristics must
   match the C source line-for-line in behavior. Preserve arithmetic order
   (floating point is not associative).
2. **Zero `unsafe`**, zero external dependencies, zero warnings. Crate roots
   allow `non_snake_case`, `non_camel_case_types`, `non_upper_case_globals`.
3. Public API keeps exact C names and return-flag conventions
   (`CV_SUCCESS = 0`, `IDA_SUCCESS = 0`, negatives fatal, positive
   recoverable).
4. `user_data` is `Option<Box<dyn Any>>`; callbacks are plain `fn` pointers
   (ARCHITECTURE.md §3.6). Do not change signatures without updating every
   example.
5. Aliasing-sensitive vector ops: in-place method family
   (`linear_sum_with`, …) when the C call aliases operands; free functions
   (`N_VLinearSum`, …) otherwise. Where C aliases user buffers with
   internal state (CVODE's `cv_y`/`yout`), copy back at **every** return
   path.
6. Example output must match the C reference `.out` byte-for-byte via
   `sundials_utils::fmt_e/fmt_f/fmt_g` (never `{:e}`); documented
   exceptions only (see VERIFICATION.md legend).
7. Excluded backends (GPU/MPI/KLU/SuperLU/LAPACK/Fortran/XBraid) stay
   excluded — do not add FFI to "complete" them.
8. When a symbol is missing, its definition is under
   `../sundials-7.7.0/src/` or `include/` — port it into `sundials_core`,
   don't invent it.

## Workflow

- After EVERY build/test/run command: `2>&1 | tee <log>`, then Read the log
  before editing. Never re-run a command that returned no visible output.
  ≤2 attempts per failing command, then switch strategy.
- Commit after every ported file or small group; tag phase gates
  (`phase0-cvode-green`, `phase2-cvodes-green`, …). This repo's git history
  is the undo mechanism.
- Phase gate = workspace builds warning-free + tests green + harness
  summary clean for all crates ported so far (cvode regression included).
- Resume after compaction from this file + PROGRESS.md + `git log` — do not
  re-explore.
