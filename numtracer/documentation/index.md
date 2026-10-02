# NumTracer

**NumTracer** is a C++20 engine that builds and contracts tensor networks —
metrics, momentum vectors, projectors, Dirac gamma matrices, and SU($N$) factors —
and **generates flat, straight-line C++ kernels** from them, with no symbolic-algebra runtime and
no external dependency.

It is a *general* engine: the physics lives in the network you hand it, not in the engine. The
worked examples in this guide are quantum-field-theory loop integrands (functional-Renormalization-
Group flows for Yang–Mills and QCD), but nothing in the contraction or the codegen is specific to
them — see [Bring your own network](getting_started/bring-your-own-network.md).

NumTracer generates kernels by **numeric tracing**. Every momentum is given explicit components
in a kinematic *frame* (some of them symbols: magnitudes, angles), and each diagram is contracted
over those components — the Dirac trace as $4\times4$ matrix products, the Lorentz indices summed,
the SU($N$) factor folded to an exact number. What comes out is a small polynomial in the frame's
symbols, which is lowered to straight-line real arithmetic in seconds.

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} Getting started
:link: getting_started/index
:link-type: doc
What NumTracer does, the mental model behind a traced network, a FORM-to-NumTracer dictionary,
and how to build, install, and consume the library.
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
