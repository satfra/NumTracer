# Installation

NumTracer builds as a small static C++20 library with a CMake package. It also installs the
Mathematica front-end on Wolfram's application path.

## Requirements

* A C++20 compiler (GCC or Clang).
* CMake ≥ 3.20.
* Optional, only for the Mathematica code generator: a Wolfram installation (`wolframscript`) and
  [FunKit](https://github.com/satfra/FunKit).

```{admonition} Two usage paths — what each needs
:class: tip
- **C++ API** (trace and contract networks, evaluate them, lower a polynomial to straight-line C++):
  needs *nothing but a C++20 compiler*.
- **Mathematica code generator** (`MakeNTKernel`, which writes a complete kernel with its `fill()`,
  signature and dressings): needs a Wolfram kernel **and FunKit** — `MakeNTKernel` uses FunKit's C++
  emitter, even when the network itself is hand-built. FunKit's flow importer `FromFunKit` is
  optional on top of that.

See also [Bring your own network](bring-your-own-network.md).
```

## Install

Get the repository and build + install:

```bash
git clone https://github.com/satfra/NumTracer.git
cmake -S NumTracer/numtracer -B NumTracer/build -DCMAKE_BUILD_TYPE=Release
cmake --build NumTracer/build -j4
cmake --install NumTracer/build
```

By default this installs to `~/.local/share/NumTracer/` - (pass `-DCMAKE_INSTALL_PREFIX=...` for another location). The install contains:

* `include/numtracer/...` — the headers;
* `lib/cmake/NumTracer/` — the CMake package config (`find_package` support).

If a `wolframscript`, `wolfram` or `math` executable is found at configure time, the Mathematica package is additionally installed to `$UserBaseDirectory/Applications/NumTracer`, so that `Get["NumTracer`"]` works from any notebook or script. Disable this with `-DNUMTRACER_INSTALL_MATHEMATICA=OFF`; if no Wolfram installation is found, the step is skipped with a warning and the C++ install is unaffected.

## Use from your own project

```cmake
find_package(NumTracer REQUIRED HINTS ~/.local/share/NumTracer)

add_executable(my_flow my_flow.cpp)
target_link_libraries(my_flow PRIVATE NumTracer::NumTracer)
```

For convenience the whole public API is available through one umbrella header —
`#include <numtracer.hpp>` — instead of the individual subsystem headers.

The [tutorial programs](../tutorials/index.md) in the repository's `Tutorials/` folder are a standalone CMake project that consumes NumTracer the way any external project would.
Kernels emitted by the codegen are self-contained by default — they include only two NumTracer
headers, `numtracer/codegen/runtime.hpp` (which supplies `numtracer::complex` and
`numtracer::compute::{powr,pow,sqrt,fma}`) and `numtracer/sun/sun_data.hpp` — so a generated kernel
drops into a consumer build with no other dependency.
A consumer that provides its own support API redirects the codegen via the `"RuntimeInclude"` / `"SupportNamespace"` options (see [codegen](../internals/codegen.md)).

## Build and test (in-tree)

```bash
ctest --test-dir NumTracer/build
```

This runs the unit tests through CTest. The tests are built whenever `NUMTRACER_BUILD_TESTS` is on,
which is the default when NumTracer is the top-level project, and `cmake --install` runs them first
unless you pass `-DNUMTRACER_TEST_ON_INSTALL=OFF`.
To run a subset by name:

```bash
ctest --test-dir NumTracer/build -R numeric --output-on-failure
```

### Build options

| option | default | effect |
|---|---|---|
| `NUMTRACER_BUILD_TESTS` | on top-level | the unit tests, flow-validation harnesses, and benchmarks |
| `NUMTRACER_INSTALL_MATHEMATICA` | `ON` | install the Wolfram front-end on Wolfram's application path |
| `NUMTRACER_TEST_ON_INSTALL` | `ON` | run the test suite as part of `cmake --install` |
| `NUMTRACER_TEST_CODEGEN` | `OFF` | register the Wolfram regeneration gate (`ctest -L codegen`, run with `NT_RUN_CODEGEN=1`) |
| `NUMTRACER_SANITIZE` | `OFF` | opt-in sanitizer list for the C++ build (e.g. `-DNUMTRACER_SANITIZE="ADDRESS;UNDEFINED"`) |
| `NUMTRACER_GPU_TESTS` | `OFF` | CUDA loop-integral integration tests (needs CUDA + GSL) |
| `NUMTRACER_KOKKOS_TESTS` | `OFF` | Kokkos twins of the GPU integration tests |
| `NUMTRACER_BUILD_ZA4_147` | `OFF` | the ~15 MB four-gluon 1/4/7 kernel test (needs a locally generated kernel) |
| `NumTracer_DOCUMENTATION` | `OFF` | build this documentation site (see below) |

## Build the documentation

The documentation needs a system `doxygen` plus a Python 3.9–3.12 environment (doxysphinx does not support 3.13+):

```bash
python3.12 -m venv .venv && source .venv/bin/activate # replace python3.12 by your executable
pip install -r documentation/requirements.txt
```

Then either run the script directly:

```bash
bash NumTracer/numtracer/documentation/build.sh
```

or drive it through CMake:

```bash
cmake -B NumTracer/build -DNumTracer_DOCUMENTATION=ON
cmake --build NumTracer/build --target documentation
```

The site is then written either to `NumTracer/build/documentation/html/index.html` or `NumTracer/numtracer/documentation/_build/html/index.html`.
