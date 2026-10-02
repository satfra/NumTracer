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

A symbolic tensor-algebra system does this algebra *ahead of time* and emits a flat polynomial in
the scalar products, which a C++ compiler turns into a fast kernel. NumTracer produces the same
kind of kernel, but does the contraction in plain C++ as a build-time step — no symbolic-algebra
runtime.

```{admonition} Relationship to symbolic tracers
:class: note
This is the one place we compare NumTracer to a symbolic tensor-algebra system such as FORM (used,
via FormTracer, only as a validation oracle for the test suite). Such a tool does the tensor
algebra symbolically ahead of time and emits a flat polynomial in the scalar products; NumTracer
produces the *same kind* of flat scalar kernel, but does the contraction numerically in C++ over a
fixed frame. In practice it generates a kernel **~80–175× faster**, and the generated kernel is
**competitive with or faster than** the symbolic one: on the quark–gluon vertex `ZAqbq{1,4,7}_147`
it runs at **0.96× / 0.99× / 0.62×** the FORM kernel's time (`tests/refshim/bench_aqbq147.cpp`),
and on the pure-gauge `ZA3_147` at 1.01×. See [PERFORMANCE.md](../../PERFORMANCE.md) for the
per-flow table. Nothing downstream depends on that tool: the emitted kernel is self-contained C++.

A fixed-frame contraction was long assumed to be structurally unable to do the partial-fractioning
(integration-by-parts-like) step a symbolic tracer performs on scalar products before the frame is
substituted, and that was believed to leave an irreducible residual. It does not: the cancellation
can be done *in* the frame, by exact polynomial division of a shifted-line propagator denominator
into the numerator that contains it (`divThroughPolyAtoms`, see
[the numeric engine](numeric-engine.md)). That closed the residual, and improved
accuracy at the same time — the division and the terms that cancel against it both disappear.
```

## The pipeline

```text
   DSL network ──NumTrace──▶ diagrams (coeff × contraction) + frame + env layout
        │
        ▼  MakeNTKernel
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
| [Cx: the compile-time complex type](expression-algebra.md) | the `Cx` NTTP complex type (tensor entries, coefficients) and the `Lit` constant carrier |
| [CSE and Horner lowering](cse-and-lowering.md) | turn a polynomial into a fast, straight-line real kernel (Horner + real value-numbering) |
| [Sector data](sectors.md) | the typed-out gamma and SU(N) tables and the chiral dense trace |
| [Worked example](worked-example.md) | the quark self-energy, end to end |

## Why it runs at *build* time, not in the consumer's compiler

Contracting networks *inside the consumer's C++ compiler* — every tensor entry encoded as a
type — is correct but expensive at scale: a trace of several transverse projectors expands
explosively in the frame-component basis, and the compiler never reclaims the intermediate memory
it allocates during constant evaluation, so peak RAM grows with *total* allocations rather than
the live set.

The numeric path runs the identical contraction as a *generator program* instead — numerically,
over a fixed loop frame. The generator runs on the native CPU in seconds and tens of megabytes,
prints a small flat kernel, and the consumer only ever compiles that. The generator's polynomial is
turned into straight-line real arithmetic by the **lowering stage** — greedy Horner factoring plus
real value-numbering (CSE).

## Terminology

Frame
: A choice of reference components for every momentum in a diagram (e.g. a symmetric point with
  one loop angle). It fixes how each scalar product is written in terms of the kernel's runtime
  arguments.

Scalar symbol
: A scalar product (`l·p`, `l²`, …) kept as a named quantity. The generated kernel computes
  each once per call and the lowered arithmetic is a polynomial in them — the *invariant basis* of
  scalar products.

Diagram coefficient
: The scalar (dressings, regulators, propagator denominators) multiplying a diagram's
  contraction. NumTracer owns the contraction; the coefficient is ordinary C++, emitted by
  FunKit's COEN.

Numeric pruning
: A zero (or round-off-tiny) term never materialises in the contracted polynomial — a product
  over largely-sparse gamma matrices collapses to the few surviving monomials instead of $4^k$ of
  them. It is what keeps the lowering of a sparse trace small. See
  [the numeric engine](numeric-engine.md).
