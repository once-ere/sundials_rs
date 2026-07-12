# sundials_rs — The Complete Guide

*A pure-Rust translation of SUNDIALS v7.7.0 — six solvers, one shared core —
explained from the ground up.*

This guide assumes you know basic algebra and a little calculus (what a
derivative is), and have seen at least one program before. By the end you
should be able to explain to your class — and your teacher — what this
library does, why it works, how the pieces fit together, and how to write a
program that uses any of the six solvers. For a much deeper dive into one
solver, `../cvode_rs/cvode.md` walks through CVODE alone at the same level
of detail; everything there applies here unchanged.

---

## Part 1 — What problems does SUNDIALS solve?

SUNDIALS ("SUite of Nonlinear and DIfferential/ALgebraic equation Solvers",
from Lawrence Livermore National Laboratory) is six related solvers:

| Crate | Solves | One-line description |
|---|---|---|
| `cvode_rs` | ODEs | adaptive multistep integrator (Adams + BDF) |
| `cvodes_rs` | ODEs + sensitivities | CVODE plus forward & adjoint sensitivity analysis |
| `arkode_rs` | ODEs | Runge-Kutta family: explicit, implicit, mixed, multirate, symplectic |
| `kinsol_rs` | nonlinear systems | solve `F(u) = 0` (no time at all) |
| `ida_rs` | DAEs | implicit systems `F(t, y, y′) = 0` |
| `idas_rs` | DAEs + sensitivities | IDA plus forward & adjoint sensitivity analysis |

All six sit on one shared crate, `sundials_core`, holding the vectors,
matrices, linear solvers, nonlinear solvers and utilities they have in
common — exactly mirroring how the C library shares `src/sundials/`,
`src/nvector/`, `src/sunmatrix/`, `src/sunlinsol/`, `src/sunnonlinsol/`.

### 1.1 Ordinary differential equations (ODEs)

Many things in science change continuously: chemical concentrations, planet
positions, currents, populations. A **differential equation** describes *how
fast* something changes instead of *what its value is*:

```
dy/dt = f(t, y),        y(t0) = y0
```

Read it as: "the rate of change of `y` at time `t` equals the formula
`f(t, y)`, and at the starting time `t0` the value is `y0`." This pair — a
rate equation plus a starting value — is an **initial value problem (IVP)**.
`y` is usually a whole list (a *vector*) of `N` numbers.

For most real problems there is *no formula* for the answer, so a computer
marches forward in small time steps, using `f` to predict where `y` goes
next. CVODE and ARKODE do this marching, in different ways.

### 1.2 Stiffness — the reason there are several integrators

A problem is **stiff** when it contains both very fast and very slow
behavior at once (a reaction that finishes in microseconds inside a
simulation that runs for days). Simple ("explicit") methods are forced to
take microsecond-sized steps *forever* or they blow up. **Implicit** methods
stay stable with big steps, but each step requires solving an equation
system — which is why the suite carries its own linear and nonlinear
solvers.

- **CVODE/CVODES**: Adams multistep methods (orders 1–12) for non-stiff
  problems, BDF multistep methods (orders 1–5) for stiff ones.
- **ARKODE**: one-step Runge-Kutta methods — explicit (ERK), diagonally
  implicit (DIRK), *additive* mixes of both (ARK/IMEX, splitting the RHS
  into a stiff implicit part and a cheap explicit part), low-storage
  super-time-stepping (LSRK), symplectic partitioned RK for mechanics
  (SPRK), multirate methods that run a fast integrator inside each slow
  stage (MRI), and operator-splitting/forcing composition of arbitrary
  steppers.

### 1.3 Differential-algebraic equations (DAEs)

Some models cannot be written as `dy/dt = formula`. A pendulum constrained
to a rod, or an electric circuit obeying Kirchhoff's laws, mixes
*differential* equations with pure *algebraic* conditions ("the length is
always 1"). The general form is

```
F(t, y, y′) = 0,        y(t0) = y0,  y′(t0) = y0′
```

IDA/IDAS integrate these directly with BDF methods. A DAE adds a twist: the
initial values must be *consistent* (they must already satisfy the algebraic
conditions), so IDA ships a starter, `IDACalcIC`, that computes consistent
initial conditions for you.

### 1.4 Plain nonlinear systems

Sometimes there is no time at all — you just need the `u` where
`F(u) = 0` (a steady state, an equilibrium). KINSOL solves exactly this,
with Newton's method (with line search), Picard iteration, or fixed-point
iteration accelerated by Anderson acceleration.

### 1.5 Sensitivity analysis — the "S" in CVODES/IDAS

Real models have parameters `p` (rate constants, lengths, temperatures).
**Sensitivity analysis** asks: *how much does the solution change if a
parameter changes?* Two strategies:

- **Forward sensitivity analysis (FSA)**: integrate the derivative vectors
  `s_i(t) = ∂y(t)/∂p_i` *alongside* `y` itself. Cost grows with the number
  of parameters; great for a few parameters and many outputs.
- **Adjoint sensitivity analysis (ASA)**: integrate a single companion
  system *backwards* in time from the final state to compute the gradient
  of one scalar functional `G` with respect to *all* parameters (or initial
  conditions) at once. Cost is roughly independent of the number of
  parameters; great for many parameters and one output (the workhorse of
  optimization and machine learning — this is "backpropagation through an
  ODE"). Because the backward pass needs the forward solution at every
  time, the forward pass saves **checkpoints** and replays/interpolates
  between them (cubic Hermite or variable-degree polynomial interpolation).

CVODES layers both on CVODE; IDAS layers both on IDA. FSA offers three
correctors — `CV_SIMULTANEOUS` (solve states+sensitivities as one big
nonlinear system), `CV_STAGGERED` (states first, then all sensitivities),
`CV_STAGGERED1` (states, then each sensitivity one at a time) — and either
user-supplied sensitivity right-hand sides or internal difference-quotient
(DQ) approximations. Second-order information (Hessians) comes from running
ASA *over* FSA (`CVodeInitBS`/`CVodeQuadInitBS`, see the
`cvsHessian_ASA_FSA` example).

---

## Part 2 — How each solver works, in five paragraphs each

### 2.1 CVODE / CVODES — variable-order multistep

CVODE stores the recent history of the solution as a **Nordsieck array**
`zn` — the Taylor-like column `[y, h·y′, h²·y″/2, …]` — and each step (a)
*predicts* by shifting that array forward, (b) *corrects* by solving the
implicit method equation with the chosen nonlinear solver, (c) runs a
weighted-error **test**, and (d) picks the next step size *and order* by
comparing the efficiency of order q−1, q, q+1. Everything else — Jacobian
reuse (`jok/jcur`), the γ-change heuristics, the `etamax` ladders, stability
limit detection (STALD), rootfinding on user functions `g(t,y)`, inequality
constraints, projection onto invariants — hangs off this skeleton.
CVODES adds quadrature integration (extra variables `dq/dt = fQ(t,y)` that
do not feed back into `y`), FSA, and the checkpoint-based ASA described
above. The full story, structure by structure, is in `cvode_rs/cvode.md`.

### 2.2 ARKODE — one-step Runge-Kutta, seven steppers

ARKODE is a *framework*: a shared driver (`arkode.rs` — step loop, error
test, adaptive controller, rootfinding, relaxation) into which pluggable
**steppers** install function pointers. The steppers:

- **ARKStep** — additive RK: `y′ = fE(t,y) + fI(t,y)`, explicit table for
  `fE`, diagonally-implicit table for `fI`, Newton or fixed-point corrector,
  ARKLS linear-solver interface, predictors, IMEX pairs up to order 5.
- **ERKStep** — pure explicit RK, lighter bookkeeping than ARKStep.
- **SPRKStep** — symplectic partitioned RK for separable Hamiltonian
  systems (`ark_harmonic_symplectic`), with optional compensated summation.
- **LSRKStep** — low-storage super-time-stepping (RKC/RKL) whose stage
  count adapts to the dominant eigenvalue, estimated by a user function or
  a `SUNDomEigEstimator` (power iteration / Arnoldi).
- **MRIStep** — multirate infinitesimal methods: each slow stage advances a
  *fast* inner integrator (any `SUNStepper`, usually a wrapped ARKStep)
  across the stage interval, with MIS/MRI-GARK/IMEX-MRI-SR coupling tables
  and an optional MRIHTol controller that adapts the fast tolerance.
- **SplittingStep / ForcingStep** — operator splitting (Lie-Trotter,
  Strang, higher-order compositions) and forcing composition of two or
  more `SUNStepper`s.

Adaptivity is delegated to a `SUNAdaptController` (Soderlind PID/PI/I
family, ImEx-Gustafsson, MRI-HTol wrapper), and *relaxation* (`arkRelax`)
can enforce conservation of an entropy/energy functional each step by
solving a scalar Newton/Brent problem for a relaxation factor.

### 2.3 KINSOL — Newton, Picard, fixed point

KINSOL solves `F(u) = 0`. Its modified/exact Newton iteration solves
`J(u) δ = −F(u)`, optionally with a line search that backtracks along δ
until the Goldstein-Armijo conditions hold; the linear solve goes through
the same SUNLinSol interface (dense, band, or preconditioned Krylov via
KINLS). Picard iteration splits `F = L·u − N(u)` with constant `L`, and
fixed-point iteration solves `u = G(u)`; both accept **Anderson
acceleration** (least-squares recombination of the last `m` iterates, with
QR downdating, damping, and optional delayed start) to speed convergence.

### 2.4 IDA / IDAS — BDF for implicit systems

IDA predicts with the polynomial extrapolating the last steps, then corrects
by solving `F(t, y, y′) = 0` with Newton, where `y` and `y′` are tied
together by the BDF relation `y′ = (α/h)·y + β` — so the iteration matrix is
`J = ∂F/∂y + cj·∂F/∂y′` with scalar `cj` supplied to the user's Jacobian.
`IDACalcIC` computes consistent initial conditions in two modes: `IDA_YA_YDP_INIT`
(given the differential components, find the algebraic components and
`y′`) and `IDA_Y_INIT` (given `y′`, find `y`). IDAS adds quadratures, FSA
(with the same simultaneous/staggered options) and checkpointed ASA exactly
parallel to CVODES.

---

## Part 3 — The shared core (`sundials_core`)

Every type below keeps the exact C name and semantics; the C
function-pointer "ops" tables become Rust enums or trait-free structs with
plain function fields (see Part 6 for why).

### 3.1 `NVector` — the vector

```rust
pub struct NVector { pub data: Vec<f64>, /* + context handle */ }
```

The serial N_Vector is a plain `Vec<f64>` you may index directly
(`y.data[i]` replaces C's `NV_Ith_S(y,i)`). The full kernel set is ported:
`N_VLinearSum`, `N_VConst`, `N_VProd`, `N_VDiv`, `N_VScale`, `N_VAbs`,
`N_VInv`, `N_VAddConst`, `N_VDotProd`, `N_VMaxNorm`, `N_VWrmsNorm` (the
error-test norm: RMS of `v[i]·w[i]`), `N_VWrmsNormMask`, `N_VMin`,
`N_VWL2Norm`, `N_VL1Norm`, `N_VCompare`, `N_VInvTest`, `N_VConstrMask`,
`N_VMinQuotient`, `N_VClone`, `N_VPrint`, and the fused/array operations.
**Aliasing rule**: where a C call writes its result over one of its own
inputs (`N_VLinearSum(a,x,b,z,z)`), Rust cannot take `&x` and `&mut z` on
the same vector, so in-place methods (`z.linear_sum_with(a,b,x)`,
`z.scale_inplace(c)`, …) provide the same bitwise arithmetic.

### 3.2 `SUNMatrix` — dense, band, sparse

```rust
pub enum SUNMatrix { Dense(DenseMatrix), Band(BandMatrix), Sparse(SparseMatrix) }
```

`DenseMatrix` stores column-major `data: Vec<f64>` (`A(i,j)` =
`data[j*M + i]`, matching `SM_ELEMENT_D`), `BandMatrix` the LAPACK-style
band storage with `set(i, j, v)` accessors, `SparseMatrix` CSR/CSC. The
small-dense kernels used inside preconditioners
(`SUNDlsMat_denseCopy/Scale/AddIdentity/GETRF/GETRS`) are ported too.

### 3.3 `LinearSolver` — direct and Krylov

```rust
pub enum LinearSolver { Dense(..), Band(..), SPGMR(..), SPFGMR(..),
                        SPBCGS(..), SPTFQMR(..), PCG(..), Custom(Box<dyn ..>) }
```

Constructors keep their C names (`SUNLinSol_Dense(&y, &A, &ctx)`,
`SUNLinSol_SPGMR(&y, pretype, maxl, &ctx)`, …). The Krylov solvers support
left/right preconditioning (`SUN_PREC_LEFT/RIGHT/BOTH/NONE`), modified or
classical Gram-Schmidt, and scaling vectors, byte-faithfully to the C
`SUNLinSolSolve` iterations. `Custom` lets an example supply a
matrix-embedded solver (see `cvsAnalytic_mels`).

### 3.4 `NonlinearSolver` — Newton and fixed point

`SUNNonlinSol_Newton(&y, &ctx)` and `SUNNonlinSol_FixedPoint(&y, m, &ctx)`
(plus the `*Sens` variants sized `(Ns+1)·N` for staggered/simultaneous
sensitivity correctors) implement the C `SUNNonlinearSolver` API: `Setup`,
`Solve`, `SetSysFn`, `SetLSetupFn`, `SetLSolveFn`, `SetConvTestFn`, counters.

### 3.5 Everything else in the core

- `SUNContext_Create()` — the context object every allocation takes.
- `SUNAdaptController_*` — Soderlind family (PID/PI/I/ExpGus/ImpGus/H0211…),
  ImExGus, and MRIHTol; used by ARKODE.
- `SUNStepper` — the "any integrator as a black box" interface used by
  MRIStep/SplittingStep/ForcingStep, with `evolve/one_step/full_rhs/reset`
  function slots and a `content: Option<Box<dyn Any>>` payload.
- `SUNAdjointStepper` / `SUNAdjointCheckpointScheme_Fixed`,
  `SUNDataNode`/`SUNMemoryHelper` (in-memory checkpoint storage).
- `SUNDomEigEstimator_Power` / `_Arnoldi` for LSRK.
- `sundials_cli` — the string-keyed option engine behind every
  `<Solver>SetOptions(mem, id, file, args)` ("cvodes.max_order 3" on the
  command line reaches `CVodeSetMaxOrd`).
- **`sundials_utils::fmt_e / fmt_f / fmt_g`** — printf-faithful float
  formatting (`%12.4e`, `%10.6f`, `%g`). Every example prints through these
  (never Rust's `{:e}`) — this is what makes byte-identical output
  possible.
- `SUNRabs/SUNRsqrt/SUNMAX/...` math helpers, `SUN_UNIT_ROUNDOFF`, error
  codes and `cvProcessError`-style reporting per solver.

---

## Part 4 — The user-facing API, solver by solver

Return-flag convention everywhere: `*_SUCCESS == 0`, negative values are
fatal errors (`CV_ILL_INPUT = -22`, `CV_MEM_NULL`, …), certain positive
values are informative (`CV_ROOT_RETURN = 2`, `CV_TSTOP_RETURN = 1`,
`IDA_ROOT_RETURN`, `ARK_ROOT_RETURN`, …). Function names, argument order
and flag values match the C headers exactly; only the *types* are Rust-ified
(`&mut CVodeMem` instead of `void*`, `&NVector` instead of `N_Vector`,
`Option<fn>` instead of nullable function pointers).

### 4.1 CVODE (and the shared CVODES base)

```rust
let mut mem = CVodeCreate(CV_BDF, &sunctx);          // or CV_ADAMS -> Box<CVodeMem>
CVodeInit(&mut mem, f, t0, &y0);                     // f: fn(t,&y,&mut ydot,&mut UserData)->i32
CVodeSStolerances(&mut mem, reltol, abstol);         // or SVtolerances / WFtolerances(efun)
CVodeSetUserData(&mut mem, Some(Box::new(data)));    // data: any 'static type
CVodeSetLinearSolver(&mut mem, LS, Some(A));         // attach CVLS
CVodeSetJacFn(&mut mem, Some(jac));                  // else internal DQ Jacobian
let flag = CVode(&mut mem, tout, &mut y, &mut t, CV_NORMAL);  // or CV_ONE_STEP
CVodeFree(mem);                                      // consumes the Box
```

Families (every C function exists under its exact name):
- *Lifecycle*: `CVodeCreate/Init/ReInit/Free`, `CVodeRootInit`.
- *Tolerances*: `CVodeSStolerances/SVtolerances/WFtolerances`.
- *Optional inputs* (`CVodeSet…`): user data, max order, max steps, initial/
  min/max step, stop time (+`Clear`), max error-test/conv fails, nonlinear
  coefficients, constraints, eta tuning knobs, `SetMonitorFn/Frequency`,
  projection (`CVodeSetProjFn/Frequency/ErrEst/…`), STALD, Jacobian/precond
  (`CVodeSetJacFn`, `CVodeSetPreconditioner`, `CVodeSetJacTimes`,
  `CVodeSetLinSysFn`, `CVodeSetEpsLin`, `CVodeSetLSNormFactor`, …).
- *Solve/extract*: `CVode`, `CVodeGetDky` (interpolate derivative `k` at any
  recent `t`), `CVodeComputeState`.
- *Statistics* (`CVodeGet…`): steps, RHS evals, linear-solver setups, error
  test fails, Newton iterations/failures, rootfinding info, `LastOrder/
  CurrentOrder`, `LastStep/CurrentStep/ActualInitStep`, workspace sizes,
  `CVodePrintAllStats(mem, writer, SUN_OUTPUTFORMAT_TABLE|CSV)`.
- *Preconditioner modules*: `CVBandPrecInit`, `CVBBDPrecInit/ReInit` and
  `CVDiag` (diagonal approximate Jacobian solver).
- *Resize & CLI*: `CVodeResizeHistory`, `CVodeSetOptions`.

### 4.2 CVODES — everything above **plus**

Quadratures — `CVodeQuadInit(mem, fQ, &yQ0)`, `CVodeQuadReInit`,
`CVodeQuadSStolerances/SVtolerances`, `CVodeSetQuadErrCon`,
`CVodeGetQuad/GetQuadDky`, quad statistics.

Forward sensitivities —
```rust
CVodeSensInit(&mut mem, Ns, ism, Some(fS), &yS0);    // all-at-once fS
CVodeSensInit1(&mut mem, Ns, ism, Some(fS1), &yS0);  // one-at-a-time fS1 (or None => internal DQ)
CVodeSensEEtolerances(&mut mem);                     // or SensSS/SensSV
CVodeSetSensErrCon(&mut mem, true);
CVodeSetSensParams(&mut mem, Some(&p), Some(&pbar), Some(&plist));
CVodeSetSensDQMethod(&mut mem, CV_CENTERED, 0.0);
CVodeGetSens(&mem, &mut t, &mut yS);                 // + GetSens1/GetSensDky
CVodeSensToggleOff / CVodeSensFree / CVodeSensReInit
```
plus quadrature sensitivities (`CVodeQuadSensInit`, `…EEtolerances`,
`CVodeSetQuadSensErrCon`, `CVodeGetQuadSens`) and the full sensitivity
statistics family. **Internal-DQ convention (pinned)**: when `fS` is
`None`, the user data attached to the solver *must* be the
`FSAUserData { p: Vec<f64>, user: Box<dyn Any> }` wrapper from
`sundials_core`, because the DQ machinery perturbs `p` through it
(ARCHITECTURE.md §3.6).

Adjoint (ASA) —
```rust
CVodeAdjInit(&mut mem, steps, CV_HERMITE /* or CV_POLYNOMIAL */);
CVodeF(&mut mem, tout, &mut y, &mut t, CV_NORMAL, &mut ncheck);  // forward + checkpoints
CVodeCreateB(&mut mem, CV_BDF, &mut which);
CVodeInitB(&mut mem, which, fB, tB0, &yB0);           // fB(t,&y,&yB,&mut yBdot,&mut UserData)
CVodeSStolerancesB / CVodeSetUserDataB / CVodeSetLinearSolverB / CVodeSetJacFnB
CVodeQuadInitB(&mut mem, which, fQB, &qB0);           // backward quadratures (dG/dp)
CVodeB(&mut mem, tBout, CV_NORMAL);                   // backward sweep (all B problems)
CVodeGetB / CVodeGetQuadB / CVodeGetAdjY / CVodeReInitB / CVodeQuadReInitB
CVodeGetAdjCVodeBmem(&mut mem, which) -> Option<&mut CVodeMem>  // stats of a B problem
```
Sensitivity-aware backward problems use `CVodeInitBS`/`CVodeQuadInitBS`
(their callbacks also receive `yS`), enabling Hessian computation.

### 4.3 KINSOL

```rust
let mut kmem = KINCreate(&sunctx);
KINInit(&mut kmem, func, &tmpl);                     // func: fn(&u,&mut fval,&mut UserData)->i32
KINSetUserData / KINSetFuncNormTol / KINSetScaledStepTol / KINSetConstraints
KINSetLinearSolver(&mut kmem, LS, Some(A)); KINSetJacFn(...);
KINSol(&mut kmem, &mut u, KIN_LINESEARCH, &u_scale, &f_scale);
// strategies: KIN_NONE (plain Newton), KIN_LINESEARCH, KIN_PICARD, KIN_FP
```
Fixed-point/Picard knobs: `KINSetMAA(m)` (Anderson depth), `KINSetDamping`,
`KINSetDampingAA`, `KINSetDelayAA`, `KINSetOrthAA` (MGS/ICWY/CGS2/DCGS2),
`KINSetNumMaxIters`, `KINSetMaxSetupCalls/MaxSubSetupCalls`, `KINSetEtaForm`
(`KIN_ETACHOICE1/2/CONSTANT`) …; statistics via `KINGet…` incl.
`KINGetFuncNorm`, `KINGetNumNonlinSolvIters`, `KINPrintAllStats`. BBD
preconditioner: `KINBBDPrecInit`.

### 4.4 IDA

```rust
let mut imem = IDACreate(&sunctx);
IDAInit(&mut imem, res, t0, &y0, &yp0);              // res: fn(t,&y,&yp,&mut rr,&mut UserData)->i32
IDASStolerances / IDASVtolerances / IDAWFtolerances
IDASetUserData / IDASetId (differential-vs-algebraic mask) / IDASetSuppressAlg
IDASetLinearSolver(&mut imem, LS, Some(A)); IDASetJacFn(...);   // jac gets cj
IDACalcIC(&mut imem, IDA_YA_YDP_INIT, tout1);        // consistent ICs
IDASolve(&mut imem, tout, &mut t, &mut y, &mut yp, IDA_NORMAL);
```
plus rootfinding (`IDARootInit`), constraints, `IDAGetDky`, the full
`IDASet…/IDAGet…` families, `IDAPrintAllStats`, and BBD preconditioning.

### 4.5 IDAS — IDA **plus**

`IDAQuadInit/ReInit/SStolerances/SetQuadErrCon/GetQuad`;
`IDASensInit` (simultaneous/staggered), `IDASensEEtolerances`,
`IDASetSensParams` (same pinned `FSAUserData` convention for internal DQ),
`IDAGetSens`, sensitivity statistics; `IDAQuadSensInit` family; and the
adjoint module `IDAAdjInit/IDASolveF/IDACreateB/IDAInitB(S)/IDACalcICB/
IDASolveB/IDAGetB/IDAGetQuadB/IDAGetAdjY` — same architecture as CVODES.

### 4.6 ARKODE

One driver API over every stepper — create with the stepper of your choice:
```rust
let mut am = ARKStepCreate(Some(fe), Some(fi), t0, &y0, &ctx)?;  // ERK/DIRK/ARK
let mut em = ERKStepCreate(f, t0, &y0, &ctx)?;
let mut sm = SPRKStepCreate(f1, f2, t0, &y0, &ctx)?;
let mut lm = LSRKStepCreateSTS(f, t0, &y0, &ctx)?;               // or SSP
let mut mm = MRIStepCreate(Some(fs), None, t0, &y0, inner_stepper, &ctx)?;
let mut pm = SplittingStepCreate(steppers, t0, &y0, &ctx)?;      // / ForcingStepCreate
```
then drive with the shared calls: `ARKodeSStolerances`, `ARKodeSetUserData`,
`ARKodeSetLinearSolver`/`ARKodeSetJacFn` (implicit steppers),
`ARKodeSetFixedStep` or the adaptive knobs (`ARKodeSetAdaptController`,
`ARKodeSetAdaptivityFn`, `ARKodeSetInitStep/MaxStep/MinStep`, safety/bias/
growth factors), `ARKodeEvolve(&mut mem, tout, &mut y, &mut t, ARK_NORMAL)`,
`ARKodeGetDky`, `ARKodeRootInit`, statistics (`ARKodeGetNumSteps`,
`ARKodeGetNumRhsEvals(mem, idx, &mut n)`, `ARKodePrintAllStats`, …),
`ARKodeResize`, `ARKodeReset`, `ARKodeFree(&mut slot)`.

Stepper-specific extras keep their C names: Butcher-table selection
(`ARKStepSetTableNum/TableName`, `ERKStepSetTable`, `MRIStepSetCoupling` and
the `ARKODE_MRI_GARK_*`/`ARKODE_IMEX_MRI_SR*` tables, `ARKodeSetOrder`),
IMEX splitting (`ARKStepSetExplicit/Implicit/ImEx`), `ARKStepSetLinear`,
predictor choice, `MRIStepGetNumInnerStepperFails`, SPRK method selection
(`ARKODE_SPRK_MCLACHLAN_4_4`, compensated sums), LSRK dominant-eigenvalue
hookup (`LSRKStepSetDomEigFn` or a `SUNDomEigEstimator`), relaxation
(`ARKodeSetRelaxFn(rfn, rjac)` + `ARKodeSetRelax*` tuning family), the
`ARKBandPrecInit`/`ARKBBDPrecInit` preconditioner modules, and
inner-stepper wrapping via `ARKodeCreateMRIStepInnerStepper` /
`ARKodeCreateSUNStepper`.

---

## Part 5 — A worked example per solver (sketches)

Complete, runnable versions live in each crate's `examples/` directory —
one `.rs` per upstream C example, byte-for-byte output-compatible. Start
with these:

| To learn… | Read |
|---|---|
| basic stiff ODE + dense Newton | `cvode_rs/examples/cvRoberts_dns.rs` |
| Krylov + user preconditioner | `cvodes_rs/examples/cvsDiurnal_kry.rs` |
| forward sensitivities (user fS) | `cvodes_rs/examples/cvsRoberts_FSA_dns.rs` |
| forward sensitivities (internal DQ + FSAUserData) | `cvodes_rs/examples/cvsAdvDiff_FSA_non.rs` |
| adjoint gradient dG/dp | `cvodes_rs/examples/cvsRoberts_ASAi_dns.rs` |
| Hessian via ASA-over-FSA | `cvodes_rs/examples/cvsHessian_ASA_FSA.rs` |
| nonlinear system + line search | `kinsol_rs/examples/kinRoboKin_dns.rs` |
| Anderson-accelerated fixed point | `kinsol_rs/examples/kinAnalytic_fp.rs` |
| DAE + consistent ICs | `ida_rs/examples/idaRoberts_dns.rs` |
| DAE adjoint | `idas_rs/examples/idasRoberts_ASAi_dns.rs` |
| IMEX ARK | `arkode_rs/examples/ark_brusselator1D_imexmri.rs` |
| multirate MRI | `arkode_rs/examples/ark_kpr_mri.rs` |
| symplectic mechanics | `arkode_rs/examples/ark_harmonic_symplectic.rs` |
| conservation via relaxation | `arkode_rs/examples/ark_conserved_exp_entropy_ark.rs` |

The universal program shape (CVODE flavor):

```rust
use cvode_rs::*;
struct MyData { k: f64 }

fn f(_t: f64, y: &NVector, ydot: &mut NVector, ud: &mut UserData) -> i32 {
    let d = ud.as_ref().unwrap().downcast_ref::<MyData>().unwrap();
    ydot.data[0] = -d.k * y.data[0];
    0
}

fn main() {
    let ctx = SUNContext_Create();
    let mut y = N_VNew_Serial(1, &ctx);
    y.data[0] = 1.0;
    let mut mem = CVodeCreate(CV_BDF, &ctx);
    CVodeInit(&mut mem, f, 0.0, &y);
    CVodeSStolerances(&mut mem, 1e-6, 1e-10);
    CVodeSetUserData(&mut mem, Some(Box::new(MyData { k: 2.0 })));
    let a = SUNDenseMatrix(1, 1, &ctx);
    let ls = SUNLinSol_Dense(&y, &a, &ctx);
    CVodeSetLinearSolver(&mut mem, ls, Some(a));
    let mut t = 0.0;
    CVode(&mut mem, 1.0, &mut y, &mut t, CV_NORMAL);
    println!("y(1) = {}", sundials_utils::fmt_e(y.data[0], 12, 4));
    CVodeFree(mem);
}
```

---

## Part 6 — How the translation stays safe *and* faithful

- **Zero `unsafe`, zero FFI, zero external crates, zero warnings** —
  `#![forbid(unsafe_code)]` at every crate root.
- **One C file → one Rust module** with the same base name; every public C
  function keeps its exact name, argument order and return flags.
- **Function-pointer tables become enums** (`LinearSolver`, `SUNMatrix`,
  `LsModule`, `PrecModule` for the band/BBD preconditioners) dispatched by
  `match` — same call graph, no vtables to corrupt.
- **`user_data` is `Option<Box<dyn Any>>`** and callbacks are plain `fn`
  pointers, so the library stays object-safe and Send-free without
  lifetimes leaking into the public API.
- **Borrow-checker patterns**: `std::mem::take` detaches a sub-structure
  (a temp vector, a step memory, an inner integrator) so it and its parent
  can be used simultaneously, and it is put back before returning — the
  Rust equivalent of C's aliased pointers, without aliasing UB.
- **Floating-point fidelity**: arithmetic order preserved statement by
  statement; C ternaries on comparisons keep their NaN semantics; aliased
  vector kernels use dedicated in-place methods bitwise-equal to the C
  loops; all output goes through `fmt_e/fmt_f/fmt_g`.
- **When C sharing is load-bearing, Rust shares too**: e.g. the food-web
  adjoint examples put one `WebData` behind `Rc<RefCell<…>>` for both the
  forward and backward problems because the C program's single shared
  struct is written by interleaved forward-replay and backward callbacks —
  separate copies provably change the trajectory.
- **Excluded by design** (hard exclusions, documented per ledger): GPU/MPI
  backends, KLU/SuperLU/LAPACK linear solvers, Fortran interfaces, XBraid,
  and the ARKODE adjoint stepper halves that require the ManyVector
  composite vector (ASA is covered end-to-end through CVODES/IDAS instead).

## Part 7 — How correctness is proven

`tools/verify_examples.sh [crate|all]` builds every example in release
mode, runs it (including argument-encoded reference variants like
`…_-sensi_stg1_t.out`), and diffs stdout byte-for-byte against the upstream
`.out` reference; `logs/summary.txt` gets one line per run. Where the
*shipped* reference is stale (built from older sources, a foreign libm, or
different compiler flags), the same example source is compiled from the
untouched `sundials-7.7.0` C tree on this machine with
`-O3 -ffp-contract=off`, the Rust output must match **that** byte-for-byte,
and the local reference is committed under `localref/`. Every row of
`VERIFICATION.md` is `identical`, `LOCAL-C`, or a one-line-justified
documented exception; `PROGRESS.md` tracks every C file to its Rust module.
Current full-workspace sweep: 150+ verified example runs, 94 IDENTICAL to
shipped references, 59 byte-identical to local C builds, 4 documented
exceptions, plus 60+ library unit tests.

## Part 8 — Cheat sheet: which solver do I need?

1. "I have `dy/dt = f(t,y)`" → **CVODE** (or **ARKODE** if you want RK,
   IMEX splitting, multirate, or symplectic structure).
2. "…and it's stiff / has fast transients" → CVODE with `CV_BDF` + Newton +
   a linear solver, or ARKODE ARKStep implicit/IMEX.
3. "…and I need d(solution)/d(parameters)" → **CVODES** FSA (few
   parameters) or ASA (gradient of one functional, many parameters).
4. "My model is `F(t,y,y′) = 0` (constraints, circuits, mechanics)" →
   **IDA**; with sensitivities → **IDAS**.
5. "No time, just solve `F(u)=0`" → **KINSOL**.
6. Then: pick tolerances (start `reltol 1e-6`, `abstol` ~ six orders below
   the size of each variable), attach a linear solver if implicit, call the
   solve routine in a loop over your output times, and read the statistics
   at the end (`…PrintAllStats`) to see what the integrator actually did.
