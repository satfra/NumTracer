# NumTracer

NumTracer is a C++20 engine that contracts the tensor networks of quantum-field-theory loop
integrands — Lorentz, Dirac, and SU(N) structure — and **generates flat, straight-line C++ kernels**
from them. It needs no symbolic-algebra system at run time: each diagram is contracted numerically
over a kinematic frame, and the resulting polynomial is lowered to plain arithmetic.

It is a *general* engine: the physics lives in the network you hand it. The reference fixtures here
are functional-Renormalization-Group (fRG) flows for Yang–Mills and QCD, but nothing in the
contraction or the code generation is specific to them.

## A first trace

$\mathrm{tr}[\slashed p\,\gamma^\mu \slashed q\,\gamma^\nu]\,P^T_{\mu\nu}(l)$ with $q = l - p$, in a
one-angle frame:

```cpp
#include <numtracer.hpp>
namespace nt = numtracer;

nt::Frame F;
auto P = F.symbol("p"), L = F.symbol("l");
auto [C, S] = F.angle("theta");                        // sin is derived from cos
nt::Momentum p = F.momentum(P, 0, 0, 0);
nt::Momentum l = F.momentum(L * C, L * S, 0, 0);
auto [mu, nu] = F.indices<2>();

nt::Poly T = F.trace({nt::slash(p), nt::gamma(mu), nt::slash(l - p), nt::gamma(nu)},
                     nt::projT(mu, nu, l));         // a polynomial in p, l, cos
double v = F.eval(T, F.at(1.3, 0.86, 0.58)).re;      // = 4p(-3 cos l + p + 2 cos^2 p)
```

SU(N) factors fold to exact numbers through a group object:

```cpp
nt::SUN su3(3);
auto [a] = su3.adjoint<1>();
auto [i, j] = su3.fundamental<2>();
nt::Cx CFN = su3.value(su3.T(a, i, j) * su3.T(a, j, i));   // tr(T^a T^a) = 4
```

## Two ways to use it

- **C++ API** (needs only a C++20 compiler): build networks as above, contract and evaluate them,
  lower a polynomial to a straight-line C++ function (`nt::to_genprog`, `nt::emit_cpp`).
- **Mathematica code generator** (needs Wolfram and [FunKit](https://github.com/satfra/FunKit)):
  write the network in a small DSL (`ntVec`, `ntTransProj`, `ntGamma`, `ntSUNT`, …) or import a
  FunKit flow, and `MakeNTKernel` writes a complete kernel — every diagram, a `fill()` for the frame
  symbols, dressings, and the integrator-facing signature.

A generated kernel includes only two small NumTracer headers (`codegen/runtime.hpp`,
`sun/sun_data.hpp`), so the consumer build has no other dependency.

## Build & test

The CMake project root is `numtracer/`, **not** the repository root. It builds a small static
library; a header-only variant is the CMake target `NumTracer::NumTracer_headeronly`.

```bash
cmake -S numtracer -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j4
ctest --test-dir build --output-on-failure
```

Tests and benchmarks build only when NumTracer is the top-level project
(`-DNUMTRACER_BUILD_TESTS=OFF` to skip). Each generated kernel under `numtracer/tests/gen/` is
gated against a FORM or equivalence oracle over random points. GPU integration tests (CUDA + GSL)
are off by default — see `numtracer/tests/gpu/README.md`.

## Install & use from other projects

```bash
cmake --install build        # default prefix: ~/.local/share/NumTracer
```

```cmake
find_package(NumTracer REQUIRED HINTS ~/.local/share/NumTracer)
target_link_libraries(my_target PRIVATE NumTracer::NumTracer)
```

If a Wolfram kernel is found at configure time, the Mathematica front-end is also installed so
`Needs["NumTracer`"]` resolves from anywhere (disable with `-DNUMTRACER_INSTALL_MATHEMATICA=OFF`).

## Documentation

A Sphinx + Doxygen site (getting started, 22 tutorials, internals, C++ reference) lives in
`numtracer/documentation/`; build it with `documentation/build.sh`. Coming from FORM? Start with
*Getting started → Coming from FORM*. The tutorial programs are a standalone CMake project in
`Tutorials/` (`cmake -S Tutorials -B Tutorials/build && ctest --test-dir Tutorials/build`).

## Layout

| path | contents |
|---|---|
| `numtracer/include/numtracer/` | the library headers; `#include <numtracer.hpp>` pulls in the whole API |
| `numtracer/mathematica/` | the Mathematica front-end (`NumTrace`, `MakeNTKernel`, `FromFunKit`) |
| `numtracer/tests/` | unit tests, generated-kernel gates, and their fixtures |
| `numtracer/documentation/` | the documentation site |
| `Tutorials/` | the tutorial programs the documentation walks through |
