# NumTracer

NumTracer is a C++20 engine that builds and contracts the tensor networks that appear in
quantum-field-theory loop integrands — Lorentz, Dirac, and SU(N) (colour and flavour) — and
**generates flat, FORM-style C++ kernels** from them. It does the tensor algebra without a
symbolic-algebra runtime and without a FORM dependency: each diagram is contracted numerically over
a compact loop frame and lowered to straight-line real arithmetic.

It is a *general* tensor-tracing engine — the physics lives in the network you hand it, not in the
engine. The reference fixtures here are functional-Renormalization-Group (fRG) flows for Yang–Mills
and QCD, but nothing in the contraction or codegen is specific to them.

## Architecture

A kernel is produced by a four-stage pipeline, build-time generation followed by a dependency-free
runtime kernel:

1. **Front-end** (`mathematica/`). `NumTrace[net, …]` parses a tensor network written in DSL heads
   (`ntMetric`, `ntVec`, `ntTransProj`, `ntSUNf`, …) into an `NTKernel`: a list of diagrams (a scalar
   coefficient × a contraction) plus the loop frame. `FromFunKit` imports a FunKit flow into the same DSL.
2. **Numeric contraction** (`MakeNTKernel`). Emits, compiles, and runs a small C++ generator that
   contracts each diagram numerically — the Dirac trace as 4×4 chiral matrix products, the Lorentz
   network by bounded index elimination, the colour/flavour factor folded to a number. Each diagram
   collapses to one small polynomial in the frame's scalar symbols.
3. **Lowering** (`codegen/`). Greedy Horner factoring plus real value-numbering (CSE) lower that
   polynomial to a flat straight-line kernel — `trN(const double* f)` trace functions, a `fill()`, and
   the per-diagram assembly.
4. **Consumer.** The committed kernel is plain C++ that includes only `codegen/runtime.hpp`, so the
   runtime build is dependency-free.

The same contraction primitives are exposed directly in C++ for hand-built traces and as the numeric
oracle the backend is validated against.

### Module layout

Headers live under `include/numtracer/`; include paths use the full prefix, e.g.
`#include "numtracer/numeric/numeric_contract.hpp"`.

| path       | role                                                                                                                                                                                                                                                 |
| ---------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `core/`    | foundations: the `constexpr` complex `Cx` (usable as a template parameter) and its `Lit<C>` carrier, the dense complex matrix `Mat<N>`, the per-axis contraction planner, and build tunables                                                         |
| `dirac/`   | typed-out Euclidean γ matrices (Weyl basis) — the `constexpr` tables every Dirac contraction reads                                                                                                                                                   |
| `sun/`     | typed-out SU(2)/SU(3) colour/flavour tables plus the `SUNBuilder<N>` runtime oracle they are checked against                                                                                                                                         |
| `network/` | the Lorentz network value (`NetVal`) and its builders, plus the generator-side numeric colour and Dirac contractions                                                                                                                                 |
| `numeric/` | the **numeric matrix-product backend**: fold a diagram's Dirac trace by 4×4 spinor products and contract the Lorentz network to one polynomial (`MPoly`), then hand it to the lowering driver                                                        |
| `dense/`   | a brute-force dense tensor (`DTensor`) that contracts labelled Lorentz⊗Dirac⊗SU(N) axes entry-for-entry — the naive baseline the generated kernels are validated against                                                                             |
| `codegen/` | build-time emission: Horner-factor a polynomial (`lower.hpp`) through the real value-numbering CSE builder (`real_cse.hpp`) into a straight-line program and print the kernel header (`gen.hpp`); `runtime.hpp` is the minimal consumer-side support |
| `cuda/`    | a two-phase CUDA quadrature integrator over a grid of external momenta, consuming `__host__ __device__` generated kernels                                                                                                                            |

## Build & test

The CMake project root is `numtracer/`, **not** the repo root. It builds a small static library by
default (header-only on opt-in via `-DNUMTRACER_HEADER_ONLY=ON`).

```bash
cmake -S numtracer -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Builds use `-O3 -march=native -Wall -Wextra`. Tests and benchmarks build only when NumTracer is the
top-level project (`-DNUMTRACER_BUILD_TESTS=OFF` to skip). Each generated kernel under `tests/gen/` is
gated against a FORM or equivalence oracle over random points; regenerate them only when a flow or the
codegen changes. Opt-in `-DNUMTRACER_SANITIZE="ADDRESS;UNDEFINED"` runs the suite under ASan/UBSan.

GPU integration tests (CUDA + GSL) are off by default — see `tests/gpu/README.md`.

## Install & use from other projects

```bash
cmake --install build        # default prefix: ~/.local/share/NumTracer
```

```cmake
find_package(NumTracer REQUIRED HINTS ~/.local/share/NumTracer)
target_link_libraries(my_target PRIVATE NumTracer::NumTracer)
```

If `wolframscript` is found at configure time, the Mathematica front-end is also installed so
`Needs["NumTracer`"]` resolves from anywhere (disable with `-DNUMTRACER_INSTALL_MATHEMATICA=OFF`).

## Generating kernels

The codegen runs under `wolframscript`. The two entry points:

- `NumTrace[net, "Frame" -> frame, "Args" -> {…}]` — analyse a DSL network into an `NTKernel`.
- `MakeNTKernel[ntk, genFile, kernelFile, tracesFile, …]` — emit, compile, and run the numeric
  generator, producing the committed kernel.

Generated kernels are self-contained by default (they include `codegen/runtime.hpp`, are wrapped in a
neutral namespace, and take generic dressing parameters). A consumer with its own support API redirects
the codegen via the `"RuntimeInclude"` / `"SupportNamespace"` / `"DressingType"` options.

## Performance

NumTracer trades a small amount of runtime for a large amount of generation time. On the quark–gluon
vertex flows measured (see `numtracer/PERFORMANCE.md`):

- **Generation** is roughly **80–175× faster** than the FORM path (~18 s per structure vs ~25–56 min).
- **Runtime** kernels run at roughly **1.0–1.6× FORM** per evaluation. The lowered kernel is
  straight-line real arithmetic with no integer division, so it is GPU-friendly.

## Documentation

A Sphinx + Doxygen site (narrative guide + C++ API) lives in `numtracer/documentation/`; build it with
`documentation/build.sh` (needs `doxygen` plus a Python 3.9–3.12 environment — see that directory).
