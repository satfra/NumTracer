# Installation

NumTracer is a header-only C++20 library with a CMake package: build once,
`cmake --install`, and consume it from any project on the system via
`find_package`. The same install also places the Mathematica front-end on
Wolfram's application path.

## Requirements

* A C++20 compiler (GCC or Clang). The larger Dirac-trace chains push the
  compiler's constexpr evaluation budget; the exported
  `NumTracer::constexpr_budget` target (see below) raises the limits.
* CMake ≥ 3.20.
* Optional, for the Mathematica codegen: a Wolfram installation
  (`wolframscript`) and [FunKit](https://github.com/satfra/FunKit).

## Install

From the repository root:

```bash
cmake -S numtracer -B numtracer/build -DCMAKE_BUILD_TYPE=Release
cmake --build numtracer/build -j4
cmake --install numtracer/build
```

By default this installs to `~/.local/share/NumTracer/` — no `sudo`, nothing
system-wide (pass `-DCMAKE_INSTALL_PREFIX=...` for another location). The
install contains:

* `include/numtracer/...` — the headers;
* `lib/cmake/NumTracer/` — the CMake package config (`find_package` support).

If `wolframscript` is found at configure time, the Mathematica package
(`NumTracer.m` and its modules) is additionally installed to
`$UserBaseDirectory/Applications/NumTracer` — the same location FunKit uses —
so that `Needs["NumTracer`"]` works from any notebook or script. Disable this
with `-DNUMTRACER_INSTALL_MATHEMATICA=OFF`; if no Wolfram installation is
found, the step is skipped with a warning and the C++ install is unaffected.

## Use from your own project

```cmake
find_package(NumTracer REQUIRED HINTS ~/.local/share/NumTracer)

add_executable(my_flow my_flow.cpp)
target_link_libraries(my_flow PRIVATE NumTracer::NumTracer)
```

`NumTracer::NumTracer` is an interface target carrying the include path and the
C++20 requirement. If your contraction chains hit GCC's default constexpr
budget ("constexpr evaluation exceeded ... operations"), additionally link

```cmake
target_link_libraries(my_flow PRIVATE NumTracer::constexpr_budget)
```

which raises `-fconstexpr-ops-limit` / `-fconstexpr-loop-limit` /
`-fconstexpr-depth` (GCC) or `-fconstexpr-steps` / `-fconstexpr-depth` (Clang).

The [tutorial programs](../tutorials/index.md) in the repository's `Tutorials/`
folder are a standalone CMake project that consumes NumTracer the way any external
project would, and double as a consumer template. Kernels emitted by the codegen are
self-contained by default — they include only `numtracer/codegen/runtime.hpp` (which
supplies `numtracer::complex` and `numtracer::compute::{powr,pow,sqrt,fma}`) — so a
generated kernel drops into a consumer build with no extra dependency. A consumer that
provides its own support API redirects the codegen via the `"RuntimeInclude"` /
`"SupportNamespace"` options (see [codegen](../internals/codegen.md)).

## Build and test (in-tree)

```bash
cmake -S numtracer -B numtracer/build
cmake --build numtracer/build -j4
ctest --test-dir numtracer/build
```

This builds the unit tests (the numeric contraction engine, the Dirac trace, the
colour fold, the lowering passes), the flow-validation harnesses (each generated
kernel compared against a FORM or dense numeric oracle over random points), and the
SU(N) colour-contraction helper, and runs them all through CTest. Tests are built only
when NumTracer is the top-level project (override with `-DNUMTRACER_BUILD_TESTS=ON/OFF`).
To run a subset by name:

```bash
ctest --test-dir numtracer/build -R numeric --output-on-failure
```

### Build options

| option | default | effect |
|---|---|---|
| `NUMTRACER_BUILD_TESTS` | on top-level | the unit tests, flow-validation harnesses, and benchmarks |
| `NUMTRACER_INSTALL_MATHEMATICA` | `ON` | install the Wolfram front-end on Wolfram's application path |
| `NUMTRACER_SANITIZE` | `OFF` | opt-in sanitizer list for the C++ build (e.g. `-DNUMTRACER_SANITIZE="ADDRESS;UNDEFINED"`) |
| `NUMTRACER_GPU_TESTS` | `OFF` | CUDA loop-integral integration tests (needs CUDA + GSL) |
| `NUMTRACER_KOKKOS_TESTS` | `OFF` | Kokkos twins of the GPU integration tests |
| `NUMTRACER_BUILD_ZA4_147` | `OFF` | the ~15 MB four-gluon 1/4/7 kernel test (needs a locally generated kernel) |
| `NUMTRACER_BUILD_ZA3_147_DENSE` | `OFF` | the dense ZA3 1/4/7 reproducer + dense timing harnesses |
| `NumTracer_DOCUMENTATION` | `OFF` | build this documentation site (see below) |

## Build the documentation

The documentation needs a system `doxygen` plus a Python 3.9–3.12 environment
(doxysphinx does not support 3.13+):

```bash
python3.12 -m venv .venv && source .venv/bin/activate
pip install -r documentation/requirements.txt
```

Then either run the script directly:

```bash
documentation/build.sh
```

or drive it through CMake:

```bash
cmake -B build -DNumTracer_DOCUMENTATION=ON
cmake --build build --target documentation
```

The site is written to `documentation/_build/html/index.html`.
