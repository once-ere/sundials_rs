# cvode_rs — Architectural Blueprint

A pure-Rust, memory-safe translation of **CVODE v7.7.0** (from SUNDIALS v7.7.0),
the adaptive-order, adaptive-step multistep solver for stiff (BDF) and
non-stiff (Adams) ordinary differential equation initial value problems:

```
dy/dt = f(t, y),   y(t0) = y0,   y ∈ R^N
```

## 1. Source → Target file map

Every translated file keeps its original base name with `.rs` appended
(extension dropped), and lives in `cvode_rs/src/`.

### CVODE proper (from `cvode-7.7.0/src/cvode/`)

| C source                | Rust module            | Contents |
|-------------------------|------------------------|----------|
| `cvode_impl.h`          | `cvode_impl.rs`        | `CVodeMem` struct, internal constants (Q_MAX, L_MAX, ETA_*, BIAS*, MXNCF…), control flags |
| `cvode.c`               | `cvode.rs`             | `CVodeCreate/Init/ReInit`, tolerances, `CVode` main loop, `cvStep`, `cvAdjustParams`, `cvSet/cvSetBDF/cvSetAdams`, `cvPredict`, `cvDoErrorTest`, `cvCompleteStep`, `cvPrepareNextStep`, `cvBDFStab` (STALD), `cvRcheck1/2/3`, `cvRootfind`, `CVodeGetDky`, `cvHin`, `cvYddNorm`, `cvEwtSet*`, `cvHandleFailure` |
| `cvode_io.c`            | `cvode_io.rs`          | All `CVodeSet*` / `CVodeGet*` optional input/output functions, `CVodePrintAllStats`, `CVodeGetReturnFlagName` |
| `cvode_ls.h/_impl.h`    | `cvode_ls_impl.rs`     | `CVLsMem` struct, CVLS return codes |
| `cvode_ls.c`            | `cvode_ls.rs`          | `CVodeSetLinearSolver`, Jacobian/preconditioner setters, DQ Jacobian (dense & band), DQ Jv product, `cvLsInitialize/Setup/Solve/Free`, all `CVodeGet*` LS stats |
| `cvode_nls.c`           | `cvode_nls.rs`         | Nonlinear-system residual/fixed-point functions, Newton & fixed-point iteration drivers wired to `CVodeMem` (`cvNlsInit`, `cvNlsSolve` = SUNNonlinSolSolve specialisation) |
| `cvode_diag.c/_impl.h`  | `cvode_diag.rs`, `cvode_diag_impl.rs` | `CVDiag` diagonal-Jacobian linear solver module |
| `cvode_bandpre.c/_impl.h` | `cvode_bandpre.rs`, `cvode_bandpre_impl.rs` | Banded difference-quotient preconditioner for Krylov methods |
| `cvode_bbdpre.c/_impl.h`  | `cvode_bbdpre.rs`, `cvode_bbdpre_impl.rs`   | Band-block-diagonal preconditioner (serial reduction) |
| `cvode_proj.c/_impl.h`  | `cvode_proj.rs`, `cvode_proj_impl.rs` | Projection onto constraint manifolds (`CVodeSetProjFn`, `cvDoProjection`) |
| `cvode_resize.c`        | `cvode_resize.rs`      | `CVodeResizeHistory` |
| `cvode_fused_stubs.c`   | `cvode_fused_stubs.rs` | Fused-kernel helpers (portable scalar versions) |
| `cvode_cli.c`           | `cvode_cli.rs`         | `CVodeSetOptions` command-line/option-string processing |

### SUNDIALS support code (from `cvode-7.7.0/src/{sundials,nvector,sunmatrix,sunlinsol,sunnonlinsol}` and `include/`)

| C source                          | Rust module              |
|-----------------------------------|--------------------------|
| `sundials_types.h`                | `sundials_types.rs`      |
| `sundials_math.c/h`               | `sundials_math.rs`       |
| `sundials_errors.c/h`             | `sundials_errors.rs`     |
| `sundials_context.c/h`            | `sundials_context.rs`    |
| `nvector_serial.c/h`              | `nvector_serial.rs`      |
| `sundials_dense.c` (`sundials_direct.h`) | `sundials_dense.rs` (dense LU: `denseGETRF/GETRS`, `densePOTRF/POTRS`, `denseGEQRF/ORMQR`) |
| `sundials_band.c`                 | `sundials_band.rs` (band LU: `bandGBTRF/GBTRS`) |
| `sunmatrix_dense.c`               | `sunmatrix_dense.rs`     |
| `sunmatrix_band.c`                | `sunmatrix_band.rs`      |
| `sunmatrix_sparse.c`              | `sunmatrix_sparse.rs`    |
| `sundials_matrix.c`               | `sundials_matrix.rs` (the `SUNMatrix` enum + generic ops) |
| `sundials_iterative.c`            | `sundials_iterative.rs` (modified/classical Gram-Schmidt, `QRfact`, `QRsol`) |
| `sunlinsol_dense.c`               | `sunlinsol_dense.rs`     |
| `sunlinsol_band.c`                | `sunlinsol_band.rs`      |
| `sunlinsol_spgmr.c`               | `sunlinsol_spgmr.rs`     |
| `sunlinsol_spfgmr.c`              | `sunlinsol_spfgmr.rs`    |
| `sunlinsol_spbcgs.c`              | `sunlinsol_spbcgs.rs`    |
| `sunlinsol_sptfqmr.c`             | `sunlinsol_sptfqmr.rs`   |
| `sunlinsol_pcg.c`                 | `sunlinsol_pcg.rs`       |
| `sundials_linearsolver.c`         | `sundials_linearsolver.rs` (the `LinearSolver` enum + generic dispatch + `SUNLS_*` flags) |
| `sundials_nonlinearsolver.c`, `sunnonlinsol_newton.c`, `sunnonlinsol_fixedpoint.c` | `sundials_nonlinearsolver.rs`, `sunnonlinsol_newton.rs`, `sunnonlinsol_fixedpoint.rs` |

### Deliberately out of scope (documented, not silent)

GPU (CUDA/HIP/SYCL/RAJA/Kokkos/Ginkgo/oneMKL/magma), MPI/parallel/hypre/PETSc/
SuperLU/KLU-backed modules, OpenMP vectors, the Fortran (`fmod_*`) bindings,
XBraid, profiler/logger backends, and the C unit-test harnesses cannot be part
of a *pure Rust, dependency-free, serial* library; they bind foreign runtimes.
Their serial equivalents above preserve 100% of CVODE's mathematical
functionality. Examples requiring KLU/SuperLU (`cvRoberts_klu`,
`cvRoberts_block_klu`, `cvRoberts_sps`) are therefore not portable as-is;
every other serial example is translated.

## 2. Type system mapping

| C                        | Rust |
|--------------------------|------|
| `sunrealtype` (double)   | `f64` (`pub type Sunrealtype = f64`) |
| `sunindextype`           | `i64` |
| `sunbooleantype`         | `bool` |
| `long int` counters      | `i64` |
| `int` return flags       | `i32` |
| `N_Vector` (ptr + ops table) | `NVector` struct, concrete serial implementation, `Clone` |
| `SUNMatrix` (ptr + ops table) | `enum SUNMatrix { Dense, Band, Sparse }` |
| `SUNLinearSolver`        | `enum LinearSolver { Dense, Band, Spgmr, Spfgmr, Spbcgs, Sptfqmr, Pcg, Custom(Box<dyn CustomLinSol>) }` |
| `SUNNonlinearSolver`     | `enum NonlinearSolver { Newton(..), FixedPoint(..) }` |
| `void* cvode_mem`        | `&mut CVodeMem` (created by `CVodeCreate`, owned by caller) |
| `void* user_data`        | `UserData = Option<Box<dyn Any>>`, examples downcast with `downcast_mut::<T>()` |
| function pointers        | plain `fn` pointers (no captures needed — state travels in `user_data`) |
| `SUNContext`             | `SUNContext` unit-like struct (serial; no profiler/logger backends) |

**Polymorphism**: the C ops-table (vtable-in-struct) pattern becomes enum
dispatch. This is memory-safe, faster (no indirect calls through `void*`),
and keeps exhaustiveness checked by the compiler.

**Aliasing**: C freely calls `N_VLinearSum(a, x, b, z, z)` (output aliases an
input). Rust forbids `&x, &mut x`. The vector module therefore provides both
free functions for distinct operands *and* in-place methods
(`linear_sum_with`, `scale_inplace`, …). Each call site in the C code is
statically one or the other; the translation picks the matching form.

**`cv_y` aliasing**: in C, `cv_mem->cv_y` aliases the user's `yout` during
`CVode()`. The Rust port keeps `cv_y` owned inside `CVodeMem` and copies to
`yout` on return — identical observable behavior, no aliasing.

## 3. Core data structure contracts (pinned interfaces)

### 3.1 `NVector` (nvector_serial.rs)

```rust
pub struct NVector { pub data: Vec<f64> }
impl NVector {
    pub fn new(n: usize) -> Self;              // zero-filled, N_VNew_Serial
    pub fn from_slice(s: &[f64]) -> Self;      // N_VMake_Serial
    pub fn len(&self) -> usize;                // N_VGetLength
    pub fn ith(&self, i: usize) -> f64;        // Ith()/NV_Ith_S 0-based
    pub fn set_ith(&mut self, i: usize, v: f64);
}
```
Free functions (distinct operands), 1:1 with N_V ops on serial vectors:
`N_VLinearSum(a,x,b,y,z)`, `N_VConst(c,z)`, `N_VProd(x,y,z)`, `N_VDiv(x,y,z)`,
`N_VScale(c,x,z)`, `N_VAbs(x,z)`, `N_VInv(x,z)`, `N_VAddConst(x,b,z)`,
`N_VDotProd(x,y)->f64`, `N_VMaxNorm(x)->f64`, `N_VWrmsNorm(x,w)->f64`,
`N_VWrmsNormMask(x,w,id)->f64`, `N_VMin(x)->f64`, `N_VWL2Norm(x,w)->f64`,
`N_VL1Norm(x)->f64`, `N_VCompare(c,x,z)`, `N_VInvTest(x,z)->bool`,
`N_VConstrMask(c,x,m)->bool`, `N_VMinQuotient(num,denom)->f64`.
In-place methods for aliased call sites:
`z.scale_inplace(c)`, `z.linear_sum_with(a, b, &y)  // z = a*z + b*y`,
`z.add_const_inplace(b)`, `z.prod_with(&x) // z *= x elementwise`,
`z.invert_inplace()`, `z.abs_inplace()`.
Direct `.data` slice access is allowed inside library modules when an
aliasing pattern fits neither form (e.g. Nordsieck loops via `split_at_mut`).

### 3.2 Matrices

```rust
pub struct DenseMatrix { pub m: i64, pub n: i64, pub data: Vec<f64> } // column-major
pub struct BandMatrix  { pub n: i64, pub mu: i64, pub ml: i64, pub s_mu: i64,
                         pub ldim: i64, pub data: Vec<f64> } // col-major, ldim = s_mu+ml+1
pub struct SparseMatrix{ pub m: i64, pub n: i64, pub nnz: i64, pub sparsetype: i32, // CSC=0, CSR=1
                         pub indexvals: Vec<i64>, pub indexptrs: Vec<i64>, pub data: Vec<f64> }
pub enum SUNMatrix { Dense(DenseMatrix), Band(BandMatrix), Sparse(SparseMatrix) }
```
Generic ops on `SUNMatrix` (sundials_matrix.rs): `zero()`, `copy_to(&mut B)`,
`scale_addi(c)  // A = c*A + I`, `scale_add(c, &B)`, `matvec(&x, &mut y)`,
`clone_empty()`. Dense/Band element access: `a.get(i,j)`, `a.set(i,j,v)`,
`DenseMatrix::col_mut(j) -> &mut [f64]`; Band uses the same indexing rules as
`SM_ELEMENT_B` (column `j` storage offset `s_mu + i - j`).

LU kernels (sundials_dense.rs / sundials_band.rs), names preserved:
`SUNDlsMat_denseGETRF(a, m, n, p) -> i64` (0 success, k>0 = zero pivot at k),
`SUNDlsMat_denseGETRS(a, n, p, b)`, `SUNDlsMat_bandGBTRF(a, p) -> i64`,
`SUNDlsMat_bandGBTRS(a, p, b)` operating on the structs above.

### 3.3 Linear solvers

```rust
pub const SUN_PREC_NONE: i32 = 0;  pub const SUN_PREC_LEFT: i32 = 1;
pub const SUN_PREC_RIGHT: i32 = 2; pub const SUN_PREC_BOTH: i32 = 3;
pub const SUN_MODIFIED_GS: i32 = 1; pub const SUN_CLASSICAL_GS: i32 = 2;

pub type ATimesFn<'a> = dyn FnMut(&NVector, &mut NVector) -> i32 + 'a;
/// psolve(r, z, tol, lr): solve P z = r ; lr=1 left, lr=2 right
pub type PSolveFn<'a> = dyn FnMut(&NVector, &mut NVector, f64, i32) -> i32 + 'a;

pub enum LinearSolver { Dense(DenseLS), Band(BandLS), Spgmr(SpgmrLS), Spfgmr(SpfgmrLS),
                        Spbcgs(SpbcgsLS), Sptfqmr(SptfqmrLS), Pcg(PcgLS),
                        Custom(Box<dyn CustomLinSol>) }

impl LinearSolver {
    pub fn ls_type(&self) -> LinearSolverType;   // Direct | Iterative | MatrixIterative | MatrixEmbedded
    pub fn initialize(&mut self) -> i32;
    pub fn setup(&mut self, a: Option<&mut SUNMatrix>) -> i32;
    pub fn solve(&mut self, a: Option<&mut SUNMatrix>, x: &mut NVector, b: &NVector,
                 tol: f64, atimes: &mut ATimesFn, psolve: Option<&mut PSolveFn>,
                 s1: Option<&NVector>, s2: Option<&NVector>) -> i32;
    pub fn set_prec_type(&mut self, pretype: i32) -> i32;
    pub fn set_gs_type(&mut self, gstype: i32) -> i32;
    pub fn set_maxl(&mut self, maxl: i32) -> i32;
    pub fn set_zero_guess(&mut self, onoff: bool);
    pub fn num_iters(&self) -> i32;
    pub fn res_norm(&self) -> f64;
    pub fn last_flag(&self) -> i64;
}
```
Constructors keep C names: `SUNLinSol_Dense(&y, &A)`, `SUNLinSol_Band(&y, &A)`,
`SUNLinSol_SPGMR(&y, pretype, maxl)`, `SUNLinSol_SPFGMR`, `SUNLinSol_SPBCGS`,
`SUNLinSol_SPTFQMR`, `SUNLinSol_PCG` → return `LinearSolver`.
`SUNLS_*` status codes preserved in `sundials_linearsolver.rs`
(e.g. `SUNLS_CONV_FAIL=804 → 4`… use the v7 numbering from
`sundials_errors.h` / `sundials_linearsolver.h`).

`CustomLinSol` trait (for MATRIX_EMBEDDED user solvers like cvAnalytic_mels):
```rust
pub trait CustomLinSol {
    fn ls_type(&self) -> LinearSolverType; // MatrixEmbedded
    fn solve(&mut self, x: &mut NVector, b: &NVector, tol: f64,
             t: f64, gamma: f64, user_data: &mut UserData) -> i32;
    fn last_flag(&self) -> i64 { 0 }
}
```

### 3.4 Nonlinear solvers

```rust
pub enum NonlinearSolver {
    Newton(NewtonSolver),        // SUNNonlinSol_Newton
    FixedPoint(FixedPointSolver) // SUNNonlinSol_FixedPoint(m) with Anderson acceleration
}
```
The iteration loops live in `cvode_nls.rs` as functions over `&mut CVodeMem`
(the C design routes SUNNonlinearSolver callbacks straight back into
`CVodeMem`; collapsing the indirection is behaviour-preserving and removes a
whole class of `void*` casts). Counters (`niters`, `nconvfails`) and options
(`maxcors`, Anderson depth `m`) live in the solver structs.

### 3.5 `CVodeMem` (cvode_impl.rs)

Direct field-for-field translation of `struct CVodeMemRec` (~120 fields),
with these ownership adaptations:
- `cv_zn: [NVector; L_MAX]`, all work vectors owned `NVector`s.
- `cv_lmem: Option<CVLsMem>` / `cv_diag_mem: Option<CVDiagMem>` — the
  `linit/lsetup/lsolve/lfree` function-pointer table becomes an enum
  `LsModule { Ls, Diag, None }` dispatched in `cvode.rs`; modules are `take()`n
  out of `CVodeMem` for the duration of a call so callbacks may borrow
  `&mut CVodeMem` (replaces the C aliasing of `cv_mem` inside `lmem`).
- Rootfinding arrays are `Vec`s.
- `proj_mem: Option<CVodeProjMem>`.
- `CVodeCreate(lmm, &ctx) -> Box<CVodeMem>`; `CVodeFree` = `drop`.

`CVLsMem` (cvode_ls_impl.rs): owns `LinearSolver`, optional `SUNMatrix` `A` and
`savedJ`, `ytemp`, `x`, Jacobian fn pointers + DQ flags, preconditioner fn
pointers, counters (`nje`, `nfeDQ`, `nstlj`, `npe`, `nli`, `nps`, `ncfl`,
`njtsetup`, `njtimes`), `msbj`, `dgmax_jbad`, `jacfn`, `jtimesDQ`, `fusefactor`
etc. — field-for-field from `CVLsMemRec`.

### 3.6 User-supplied function ABI (public, used by examples)

```rust
pub type UserData    = Option<Box<dyn std::any::Any>>;
pub type CVRhsFn     = fn(t: f64, y: &NVector, ydot: &mut NVector, user_data: &mut UserData) -> i32;
pub type CVRootFn    = fn(t: f64, y: &NVector, gout: &mut [f64], user_data: &mut UserData) -> i32;
pub type CVEwtFn     = fn(y: &NVector, ewt: &mut NVector, user_data: &mut UserData) -> i32;
pub type CVLsJacFn   = fn(t: f64, y: &NVector, fy: &NVector, jac: &mut SUNMatrix,
                          user_data: &mut UserData, tmp1: &mut NVector, tmp2: &mut NVector,
                          tmp3: &mut NVector) -> i32;
pub type CVLsPrecSetupFn = fn(t: f64, y: &NVector, fy: &NVector, jok: bool, jcur: &mut bool,
                              gamma: f64, user_data: &mut UserData) -> i32;
pub type CVLsPrecSolveFn = fn(t: f64, y: &NVector, fy: &NVector, r: &NVector, z: &mut NVector,
                              gamma: f64, delta: f64, lr: i32, user_data: &mut UserData) -> i32;
pub type CVLsJacTimesSetupFn = fn(t: f64, y: &NVector, fy: &NVector, user_data: &mut UserData) -> i32;
pub type CVLsJacTimesVecFn = fn(v: &NVector, jv: &mut NVector, t: f64, y: &NVector, fy: &NVector,
                                user_data: &mut UserData, tmp: &mut NVector) -> i32;
pub type CVProjFn    = fn(t: f64, ycur: &NVector, corr: &mut NVector, epsProj: f64,
                          err: Option<&mut NVector>, user_data: &mut UserData) -> i32;
pub type CVLocalFn   = fn(nlocal: i64, t: f64, y: &NVector, g: &mut NVector,
                          user_data: &mut UserData) -> i32;  // BBD
```
Return convention preserved: `0` success, `>0` recoverable, `<0` fatal.

**FSA internal-DQ user-data convention (pinned 2026-07-11, sundials_types.rs
`FSAUserData`).** In C, `CVodeSetSensParams`/`IDASetSensParams` store the
user's `p` POINTER and the internal difference-quotient sensitivity residuals
(cvSensRhs1InternalDQ / IDASensRes1DQ / the QuadSens DQs) perturb `p[which]`
in place — the perturbation reaches the user RHS/residual because `p` aliases
the user's own parameter array. The Rust ports keep `cv_p`/`ida_p` as owned
copies, so user code that relies on the INTERNAL DQ sensitivities (fS/resS =
None) must wrap its data as `FSAUserData { p, user }`; the DQ routines mirror
each `p[which]` perturbation into `.p` through the user-data downcast (helper
`ida_dq_set_p`; cvodes must adopt the same helper before its FSA examples —
its DQ currently perturbs only the dead copy). `IDASolve` guards this: internal
DQ + non-FSAUserData user data is rejected with IDA_ILL_INPUT instead of
silently producing zero sensitivities. Analytic-resS users are unaffected.

All public API functions keep their exact C names (`CVodeCreate`, `CVodeInit`,
`CVodeSStolerances`, `CVodeSetLinearSolver`, `CVode`, `CVodeGetDky`, …) with
`cvode_mem: &mut CVodeMem` as first argument and `i32` flag returns;
out-pointers become `&mut T`. Crate root sets
`#![allow(non_snake_case)]` / `non_camel_case_types` / `non_upper_case_globals`
so original identifiers survive verbatim.

### 3.7 C-style output formatting

Examples must reproduce reference `.out` files. `sundials_utils.rs` provides
`printf`-compatible float formatting: `fmt_e(x, width, prec)` (`%W.Pe`),
`fmt_g(x, width, prec)` (`%W.Pg`), `fmt_f(x, width, prec)` — matching C's
two/three-digit exponent and rounding rules on macOS (two-digit exponent).

## 4. Control-flow architecture of the integrator (faithful to C)

```
CVode(mem, tout, yout, &tret, itask)
 ├─ initial-step: cvInitialSetup → f(t0,y0) → cvHin (initial h) → cvRcheck1
 ├─ loop:
 │   ├─ cvEwtSet (error weights), overflow/‖h‖ guards
 │   ├─ cvStep
 │   │   ├─ cvAdjustParams (order change pending? → cvAdjustOrder, cvRescale)
 │   │   ├─ cvPredict  (advance Nordsieck array zn by Pascal-triangle shifts)
 │   │   ├─ cvSet      (l coefficients, tq test quantities: cvSetBDF/cvSetAdams)
 │   │   ├─ cvNls      (Newton via lsetup/lsolve, or functional iteration)
 │   │   ├─ cvCheckConstraints
 │   │   ├─ cvDoProjection (optional)
 │   │   ├─ cvDoErrorTest (dsm = ‖acor‖·tq[2]; fail → cvRestore, shrink h, retry)
 │   │   ├─ cvCompleteStep (update zn, tau, counters, qwait)
 │   │   └─ cvPrepareNextStep (etaqm1/etaq/etaqp1 → choose q', h'; cvBDFStab)
 │   ├─ cvRcheck3 (rootfinding on each step)
 │   └─ itask tests (CV_NORMAL: interpolate with CVodeGetDky at tout)
 └─ failure handling: cvHandleFailure maps nflag → CV_* return
```

The linear-solver interface (`cvode_ls.rs`) supplies `lsetup`/`lsolve` for
matrix-based (dense/band/sparse via LU) and matrix-free (SPGMR/SPFGMR/SPBCGS/
SPTFQMR/PCG with optional user preconditioner, DQ Jacobian-vector products)
methods; `cvode_diag.rs` supplies the cheap diagonal approximation;
`cvode_bandpre`/`cvode_bbdpre` provide reusable preconditioners.

## 5. Memory-safety & optimization notes

- Zero `unsafe` in the entire crate; no FFI; no external dependencies.
- All C `malloc/free` pairs become RAII (`Vec`, `Box`, `Option`); leak class
  of bugs is structurally impossible; double-free/use-after-free eliminated.
- Enum dispatch replaces function-pointer tables (devirtualised, inlinable).
- Elementwise vector kernels compile to auto-vectorised slice loops
  (iterator chains over `&[f64]` — bounds checks hoisted).
- Counters are `i64` (matching `long int` on LP64) — no overflow in practice;
  arithmetic on `f64` is bit-for-bit the same operations as the C code.

## 6. Examples (from `examples/cvode/serial/`)

Each translated to `examples/<name>.rs` and registered as `[[example]]` in
`Cargo.toml`: cvRoberts_dns, cvRoberts_dnsL, cvRoberts_dns_uw,
cvRoberts_dns_constraints, cvRoberts_dns_negsol, cvAdvDiff_bnd, cvAdvDiff_bndL,
cvDiurnal_kry, cvDiurnal_kry_bp, cvDirectDemo_ls, cvKrylovDemo_ls,
cvKrylovDemo_prec, cvDisc_dns, cvAnalytic_mels, cvParticle_dns, cvPendulum_dns,
cvRocket_dns — plus a new `solar_system` example (outer-planet N-body problem
integrated with BDF + dense Newton).
(`cvRoberts_dnsL`→dense LU, `cvAdvDiff_bndL`→band LU: the `…L` LAPACK variants
run through the native Rust LU, preserving their numerics.)

---

# Workspace addendum (sundials_rs)

This blueprint was inherited from the verified cvode_rs donor and governs
the whole workspace. The shared layer (§§1–3) now lives in the
`sundials_core` crate; solver crates re-export it so every `crate::<mod>`
path shown above still resolves. Per-solver addenda (cvodes sensitivities,
ida/idas, kinsol, arkode steppers, adaptcontrollers, domeigest, adjoint
stack) are appended below as each phase lands.

## Addendum A — shared-core completion (Phase 1)

- `SUNAdaptController` (sundials_adaptcontroller.rs): ops table → enum
  { Soderlind, ImExGus, MRIHTol }; content structs live in the
  sunadaptcontroller_*.rs impl modules; wrong-variant impl calls return
  SUN_ERR_ARG_INCOMPATIBLE; MRIHTol owns its two sub-controllers as Boxes.
  Constructors take no SUNContext (donor convention). NewEmpty/DestroyEmpty
  C-object plumbing has no Rust counterpart; Destroy = drop.
- `SUNDomEigEstimator` (sundials_domeigestimator.rs): same pattern,
  enum { Power, Arnoldi }.
- SUNLogger / SUNProfiler are complete standalone modules; SUNContext stays
  a unit-like struct (donor contract) because the default SUNDIALS build
  compiles solver logging/profiling macros out — solvers do not call them.
- SUNHashMap is generic over V; logger keys streams by filename, profiler
  uses std HashMap internally (documented adaptation).
- FILE* → &mut dyn std::io::Write everywhere; SUNFileHandle enum
  (sundials_futils.rs) models stdout/stderr/file/NULL.
- N_VectorSensWrapper (sundials_nvector_senswrapper.rs): struct owning
  Vec<NVector>; C's borrow-wrapping collapses to ownership.

## Addendum B — cvodes design contract (Phase 2)

- cvodes_impl.rs: CVodeMem is the donor struct + quad/sens/quadsens/adjoint
  extensions; `N_Vector*` → Vec<NVector>; znS-style arrays → [Vec<NVector>; L_MAX].
- Four correctors: NLS, NLSsim, NLSstg, NLSstg1 each Option<NonlinearSolver>.
  The C senswrapper aliases zn0Sim/ycorSim/ewtSim etc. are NOT stored:
  cvodes_nls{,_sim,_stg,_stg1}.rs operate directly on CVodeMem fields
  (same collapse the donor used for cvode_nls.rs).
- Adjoint: linked lists → Vecs (ck_mem: Vec<CVckpntMem>, cvB_mem:
  Vec<CVodeBMem>, dt_mem: Vec<CVdtpntMem>); current-pointer fields become
  indices; dt content void* → enum DtpntContent { Hermite, Polynomial };
  IM fn pointers → dispatch on ca_IMtype; CVodeBMem owns a nested
  Box<CVodeMem> for the backward problem.
