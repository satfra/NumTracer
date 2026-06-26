# NumTracer

**NumTracer** is a header-only C++20 engine that builds and contracts the tensor networks
that appear in quantum-field-theory loop integrands — Lorentz, Dirac, and SU(N) (colour and
flavour) — and **generates flat, FORM-style C++ kernels** from them, without a symbolic-algebra
runtime and without a FORM dependency.

NumTracer generates kernels by **numeric tracing**. A diagram is contracted *numerically* over a
compact loop frame — its Dirac trace as 4×4 chiral matrix products, its Lorentz network by
bounded index elimination, its colour factor folded to a number — so each diagram collapses
to a small polynomial in the frame's scalar symbols, which is then lowered to straight-line
real arithmetic. The result is the same shape a FORM reference kernel has, generated in
seconds.

It is a *general* engine: the physics lives in the network you hand it, not in the engine.

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} Getting started
:link: getting_started/index
:link-type: doc
What NumTracer does, the mental model behind a traced network, and how to build, install,
and consume the library.
:::

:::{grid-item-card} Tutorials
:link: tutorials/index
:link-type: doc
A worked flow — the quark self-energy — from a tensor network in the DSL to a generated,
validated kernel.
:::

:::{grid-item-card} Under the Hood
:link: internals/index
:link-type: doc
A developer's deep dive: the numeric contraction engine, the front-end and codegen, the
CSE / real-lowering pipeline, and the typed-out sector data — with code walk-throughs.
:::

:::{grid-item-card} C++ API
:link: doxygen/NumTracer/html/index
:link-type: doc
The full header reference, rendered by Doxygen.
:::

::::

```{toctree}
:hidden:
:maxdepth: 2

Getting Started <getting_started/index>
Tutorials <tutorials/index>
Under the Hood <internals/index>
C++ API <doxygen/NumTracer/html/index>
```
