# Under the Hood

This section is the developer's guide to how NumTracer actually works. The
[overview](../getting_started/overview.md) sketches the idea in a page; here we go all the way
down to the code.

```{toctree}
:maxdepth: 1

codegen
numeric-engine
expression-algebra
cse-and-lowering
sectors
worked-example
```

## What problem is the engine solving?

NumTracer computes the *scalar kernels* that appear in quantum-field-theory loop integrands. A
loop integrand is a product of tensors — gamma matrices, momentum vectors, metrics, projectors,
colour generators — contracted together down to a single number that still depends on a few
runtime quantities (momentum magnitudes, angles, propagator dressings). That number is fed to a
numerical integrator and evaluated at hundreds of thousands of grid points.

The traditional approach (FORM, FormTracer) does the tensor algebra symbolically *ahead of
time* and emits a flat polynomial in the scalar products, which a C++ compiler turns into a
fast kernel. NumTracer produces the same kind of kernel, but does the contraction in plain
C++ as a build-time step — no symbolic-algebra runtime, no FORM trace.

## The pipeline

```text
   DSL network ──NumTrace──▶ diagrams (coeff × contraction) + frame + env layout
        │
        ▼  MakeNTKernel  "Backend" -> "Numeric"
   numeric contraction (per diagram):
        Dirac trace  ─ 4×4 chiral matrix products ─▶ polynomial
        Lorentz net  ─ bounded index elimination  ─▶ polynomial
        colour       ─ folded to a number
        │
        ▼  one small polynomial (MPoly) per diagram, in the frame's scalar symbols
   lowering: CSE + Horner  ─▶  flat straight-line kernel  (trN(f) + fill + assembly)
```

The reading order of this section follows that pipeline:

| Page | Layer |
|---|---|
| [Front-end & codegen](codegen.md) | the Mathematica DSL, `NumTrace` / `MakeNTKernel`, eager summation, the FunKit adapter |
| [The numeric contraction engine](numeric-engine.md) | matrix-product Dirac trace, Lorentz network reduction, colour fold |
| [The expression algebra](expression-algebra.md) | the compile-time expression types and structural pruning the lowering is built on |
| [CSE and real lowering](cse-and-lowering.md) | turn a polynomial into a fast, straight-line real kernel |
| [Sector data](sectors.md) | the typed-out gamma and SU(N) tables, the chiral dense trace, the numeric oracle |
| [Worked example](worked-example.md) | the quark self-energy, end to end |

## Why it runs at *build* time, not in the consumer's compiler

An earlier version of NumTracer contracted networks *inside the consumer's C++ compiler*, with
every tensor entry encoded as a type. That is correct but expensive at scale: a trace of
several transverse projectors expands explosively in the frame-component basis, and the
compiler never reclaims the intermediate memory it allocates during constant evaluation, so
peak RAM grows with *total* allocations rather than the live set.

The numeric path runs the identical contraction as a *generator program* instead — numerically,
over a fixed loop frame. The generator runs on the native CPU in seconds and tens of megabytes,
prints a small flat kernel, and the consumer only ever compiles that. This is exactly how FORM
kernels are produced; the difference is that NumTracer does the trace in C++ rather than in a
symbolic-algebra system. The expression algebra and CSE machinery survive as the **lowering
stage** that turns the generator's polynomial into straight-line real arithmetic.

## Terminology

Frame
: A choice of reference components for every momentum in a diagram (e.g. a symmetric point with
  one loop angle). It fixes how each scalar product is written in terms of the kernel's runtime
  arguments.

Scalar symbol
: A scalar product (`l·p`, `l²`, …) kept as a named quantity. The generated kernel computes
  each once per call and the lowered arithmetic is a polynomial in them — this is the
  *invariant basis*, the same one FORM uses.

Diagram coefficient
: The scalar (dressings, regulators, propagator denominators) multiplying a diagram's
  contraction. NumTracer owns the contraction; the coefficient is ordinary C++, emitted by
  FunKit's COEN.

Structural pruning
: Dropping `Zero`/`One` sub-branches from an expression as it is built, so a product touching a
  zero vanishes. It is what keeps the lowering of a sparse trace small. See
  [the expression algebra](expression-algebra.md).
