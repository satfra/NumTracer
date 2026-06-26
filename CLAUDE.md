# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

NumTracer is a standalone, dependency-free C++20 engine for the tensor traces of fRG
(functional Renormalization Group) loop integrands — Lorentz, Dirac, and SU(N) (colour
and flavour). It **generates flat, FORM-style C++ kernels** without a symbolic-algebra
runtime and without a FORM dependency: each diagram is contracted *numerically* over a
compact loop frame — its Dirac trace as 4×4 chiral matrix products, its Lorentz network by
bounded index elimination, its colour factor folded to a number — so each diagram collapses
to a small polynomial that lowers to straight-line real arithmetic. It is a *general*
engine; the physics lives in the network you hand it, not in the engine.

The actual C++ project lives in `numtracer/` (header-only library under
`numtracer/include/numtracer/`). `numtracer/README.md` is the authoritative reference for
the layout and validation status, and `numtracer/documentation/` is a Sphinx+Doxygen site —
read those before deep work. `Yang-Mills.nb` / `Yang-Mills.m` at the repo root are
Wolfram-side scratch, not part of the build.

> History: an earlier version contracted networks *inside the consumer's compiler* via a
> symbolic expression-tensor engine (`et::ETensor`) and an invariant-basis `"Generate"`
> codegen backend. Both are **removed**. The numeric matrix-product backend is the sole
> generation path. There is no `et/`, `inv/`, or `lorentz/` directory, and the Mathematica
> `"ET"` / `"Generate"` backends now error.

## Build & test

The CMake project root is `numtracer/`, **not** the repo root:

```bash
cd numtracer
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j4
ctest --test-dir build --output-on-failure
```

Builds use `-O3 -march=native -Wall -Wextra`. Tests build only top-level
(`NUMTRACER_BUILD_TESTS`, default ON when top-level). Other options: `NUMTRACER_GPU_TESTS`
and `NUMTRACER_KOKKOS_TESTS` (CUDA/Kokkos loop-integral tests, off), `NUMTRACER_BUILD_ZA4_147`
(the ~15 MB four-gluon 1/4/7 kernel, off), `NUMTRACER_INSTALL_MATHEMATICA` (default ON),
`NumTracer_DOCUMENTATION` (off). Each test target is a standalone executable that prints
`ALL TESTS PASSED` / `TESTS FAILED` and exits non-zero on failure:

```bash
ctest --test-dir build -R numeric_contract --output-on-failure   # the numeric engine gate
./build/test_numeric_contract                                    # or run the binary directly
./build/sun_contract 3                                           # colour seam: prints f^acd f^bcd = N
```

ctest names: `numeric_contract col_numeric projector_identity finiteT_quark`
(engine unit tests; `projector_identity` covers the transverse/longitudinal/electric/magnetic
projector algebra, `finiteT_quark` the heat-bath-frame Dirac trace × electric/magnetic gluon line)
and `flow_zq_num flow_za_num flow_za3_num flow_za3_147_num flow_za3_147_form flow_za4_num
flow_zaqbq1 flow_pion_quark flow_glu_quark flow_flavour_split flow_flavour_ingroup
flow_gluon_condensate color_su3 proj_num ftproj_num` (each generated kernel vs a FORM, dense, or
equivalence oracle; `ftproj_num` is the finite-T electric/magnetic + `ntSPS`/`ntVec[q,0]` path).
`flow_flavour_split` dresses the u/d-quark loops with independent runtime dressings as one diagram
per flavour (no SU(Nf) group). `flow_flavour_ingroup` / `flow_gluon_condensate` dress components
*within* a group via a **group-diagonal dressing** (`ntSUNDiagFund`/`ntSUNDiagAdj`): a δ tagged with
a per-component dressing array, folded by `sun_value_dressed` (network/sun_net.hpp) to a polynomial
`Σ_a c_a Z_a` over runtime per-component dressing leaves `Z[a](scale)` — the Dirac trace is computed
ONCE, not one diagram per component. `flow_flavour_ingroup` reproduces `flow_flavour_split` exactly
(u/d inside the SU(2) trace); `flow_gluon_condensate` is the adjoint case (a condensed gluon
`Z_A^a`, 8 components), validated by the collapse property (all `Z_a` equal == the colour-blind
kernel). Existing flows carry no diagonal-dress factor, so their kernels regenerate byte-identical.

Install (headers + CMake package config; default prefix `~/.local/share/NumTracer`):

```bash
cmake --install build
```

Consumers use `find_package(NumTracer REQUIRED HINTS ~/.local/share/NumTracer)` and link
`NumTracer::NumTracer` (plus `NumTracer::constexpr_budget` for raised GCC/Clang constexpr
limits if long lowering chains hit them). If `wolframscript` is found, the install also
places the Mathematica package in `$UserBaseDirectory/Applications/NumTracer` so
`Needs["NumTracer`"]` resolves. **`cmake --install` does not prune headers removed from the
source**, so after a reorg a stale prefix can shadow the in-tree tree via `find_package` —
reinstall (and delete orphaned dirs) so the installed tree matches `include/numtracer/`.

The repo-root `Tutorials/` folder is a standalone CMake project of the documentation's
tutorial programs (the runnable, hand-coded numeric path). It consumes the installed package,
falling back to the in-tree library:

```bash
cmake -S Tutorials -B Tutorials/build && cmake --build Tutorials/build
ctest --test-dir Tutorials/build        # runs every tutorial
```

## Architecture: the numeric contraction engine

The library is one pipeline. Understanding these stages is enough to navigate everything.

**1. Front-end** (`mathematica/`). `NumTrace[net, …]` analyses a tensor network written in
DSL heads (`ntMetric`, `ntVec`, `ntTransProj`, `ntSUNf`, `ntSUNDeltaAdj`, `ntGamma`, `ntDress`, …)
into an `NTKernel`: a list of diagrams, each a scalar coefficient times a contraction, plus
the loop frame and env-id layout. `FromFunKit` imports a FunKit flow into the same DSL.
Vertex structure-sums are kept eager (no product-of-sums monomial blow-up); a colour↔Lorentz
bridging sum is distributed so the sectors never fuse.

**2. Numeric contraction** (`numeric/`, `network/`, namespaces `numtracer::numeric` /
`numtracer::network`). `MakeNTKernel[ntk, genFile, kernelFile, tracesFile, "Backend" ->
"Numeric"]` emits a small C++ generator and runs it. Per diagram it contracts:
   - the **Dirac trace** as 4×4 chiral matrix products (`numeric/mpoly.hpp`: `Mat4`,
     `gammaC`/`slashC`/`matmul`/`mtrace`; incl. γ5, no `(2n−1)!!` Wick blow-up), tokens
     `network/dirac.hpp` (`DiracNet` of `dgamma`/`dslash`/`dg5`/`dcomm*`);
   - the **Lorentz network** of metrics/vectors/projectors/Levi-Civita by bounded index
     elimination (`numeric/numeric_contract.hpp`: `numeric_value` / `numeric_value_netval`,
     builders `nmet`/`nvec`/`nproj`/`neps`; no `2^np` projector-mask blow-up);
   - the **colour factor** to a number (`network/sun_net.hpp`: `sun_value_cx`).
   Each diagram becomes one `MPoly` (`numeric/mpoly.hpp`) — a polynomial in the frame's
   scalar symbols carrying surviving `1/k²` atoms. `numeric/numeric_driver.hpp` drives all
   diagrams.

**3. Lowering** (`expr/`, `codegen/`). The complex polynomial is lowered to a flat
straight-line real kernel: `expr/cse.hpp` (value-level hash-consed DAG), `expr/real_cse.hpp`
(split re/im, hash-indexed CSE, DCE), then `codegen/lower.hpp` (Horner) + `codegen/gen.hpp`
(shared `f[]` symbol env across diagrams). The emitted kernel is `trN(const double* f)` trace
functions + a `fill()` + the per-diagram assembly — the shape of a FORM reference kernel,
self-contained against `codegen/runtime.hpp` (`numtracer::complex`, `compute::{powr,pow,
sqrt,fma}`).

**Supporting pieces.** `dense/dtensor.hpp` (`dense::DTensor`) is a numeric tensor over
labelled `Ax<Id, Dim>` axes that contracts entry-for-entry with no pruning — the **oracle**
generated kernels are validated against, and the `"Dense"` codegen target. `dirac/` holds the
typed-out Euclidean γ tables (`dirac_data.hpp`) and a chiral 2×2-block dense trace
(`dense_trace.hpp`, `chiral_gamma_trace<K>`). `sun/` holds the typed-out SU(2)/SU(3) tables
and the `SUN<N>` runtime oracle. `core/` is `cx.hpp` (the `constexpr` complex `Cx`, an NTTP),
`axplan.hpp` (the `Ax<Id, Dim>` axis language and contraction planner), and `config.hpp`. `cuda/`
is a two-phase quadrature integrator over generated `__host__ __device__` kernels.

Namespaces: each library dir maps to its namespace — `numtracer::{core, expr, dense, dirac,
sun, network, numeric, cuda}` (so the γ tables / `chiral_gamma_trace` live in `numtracer::dirac`
and the `SUNBuilder`/`SUNData` oracle in `numtracer::sun`). `Cx` (a `constexpr` NTTP) stays at
top level `numtracer::Cx`. Two deliberate exceptions in `codegen/`: `runtime.hpp` defines
`numtracer::complex` / `numtracer::compute` — the **default support namespace** every generated
kernel imports (`using namespace numtracer::compute;`), so it must stay in `numtracer`; and the
build-time generator helpers `gen.hpp` / `lower.hpp` live in `numtracer::network` (they operate
on `network` types and are assembled into the `numtracer::network` generator program).

## Hand-coding & validation

The numeric path is fully callable from standalone C++ — `tests/test_numeric_contract.cpp`
is the canonical worked example (build a `DiracNet` + `NNet`, contract with `numeric_value`,
evaluate the `MPoly` at a point, compare to the dense `chiral_gamma_trace<K>` oracle). The
repo-root `Tutorials/` programs build the quark self-energy (FlowZq) numerator this way and
self-check. Generated flows are validated against **copied FormTracer oracles**
(`tests/refshim/*_kernel.hh`) and the dense backend by the `compare_*_num` drivers over
~200k random points; `tests/refshim/shim.hpp` provides standalone dressing/regulator stubs so
the copied kernels compile. The committed kernels under `tests/gen/*.hh` keep the C++ build
dependency-free; regenerate only when a flow or the codegen changes (needs wolframscript +
FunKit, and FORM for the YM/QCD tensor bases).

## Conventions

- **Header-only.** Generation runs as a build-time C++ program (the trace is plain numeric
  C++, not a constant-evaluation in the consumer's compiler), printing a flat kernel; the
  consumer only compiles that. Keep runtime work in the lowered kernel and the dense path.
- **Typed-out tables.** Fixed group/dimension constants (SU(2)/SU(3) tables, gamma matrices)
  are written out as `constexpr`/typed literals in `*_data.hpp` headers and cross-checked
  against a runtime oracle (`SUN<N>`, `chiral_gamma_trace`, `dense::DTensor`), rather than
  computed at runtime in the hot path. `Nf`/`Nc` are compile-time integers.
- Prefer templated/compile-time, stack memory, and sparsity in the hot path.
- Documentation is Doxygen-commented (`@file`, `@brief`, `@tparam`); the Sphinx+Doxygen site
  lives in `numtracer/documentation/` (off by default — see README for the toolchain).
