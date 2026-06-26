# Overview

NumTracer assembles and evaluates the tensor networks that appear in quantum-field-theory
loop integrands, and **generates the scalar kernels** that a numerical integrator then runs.
It exists to reproduce the kernels that tools like FORM emit — same answer, comparable speed —
but it does the tensor algebra in plain compile-time C++, with no symbolic-algebra runtime and
no FORM dependency.

## What it computes

A loop integrand is a product of tensors — gamma matrices, momentum vectors, metrics,
projectors, colour generators — contracted together (summed over shared indices) down to a
single number that still depends on a few runtime quantities: momentum magnitudes, angles, and
propagator dressings. That number is fed to an integrator and evaluated at hundreds of
thousands of grid points.

NumTracer is a *general* tensor-tracing engine. The physics is entirely in the network you
hand it; the engine only contracts indices and folds the result. The worked flows here are
functional-Renormalization-Group (fRG) flows for Yang–Mills and QCD, but nothing in the
contraction or the codegen is specific to them. The numeric path handles Lorentz, Dirac, and
SU($N$) colour *and* flavour structure, dresses individual flavours or colour/flavour
components differently (per-flavour split and group-diagonal dressings), and factors a diagram
that disconnects into several closed traces into a product — see [Key concepts](concepts.md)
and the [dressed-flavour tutorial](../tutorials/dressed-flavour.md).

## The numeric path, end to end

NumTracer generates kernels by **numeric tracing** — a build-time generator driven from a small
Mathematica front-end:

1. **Front-end.** `NumTrace` analyses a tensor network written in DSL heads (`ntMetric`,
   `ntVec`, `ntTransProj`, `ntSUNf`, `ntSUNDeltaAdj`, `ntDress`, …) into a list of diagrams — each a
   scalar coefficient times a contraction — together with the loop frame that fixes every
   momentum's components and the env-id layout. `FromFunKit` imports a FunKit flow into the
   same DSL. Vertex structure-sums are kept eager, so no product-of-sums monomial blow-up ever
   forms; the same holds for dressed propagator numerators, which are collected into a single
   trace rather than distributed into `2^D` diagrams (`"DressingCollection"`, on by default).

2. **Numeric contraction.** `MakeNTKernel` emits a small C++ generator and runs it. The
   generator contracts each diagram numerically over the frame: the Dirac trace as 4×4 chiral
   matrix products (including γ5; no `(2n−1)!!` Wick blow-up), the Lorentz network by bounded
   index elimination (no `2^np` projector-mask blow-up), and the colour factor folded to a
   number. Each diagram becomes a small polynomial in the frame's scalar symbols.

3. **Lowering.** Common-subexpression elimination and Horner factoring lower that polynomial
   to a flat straight-line kernel — trace functions over the frame symbols, a `fill` that
   computes the symbols once per call, and the per-diagram assembly. The committed kernel is
   plain C++ that depends only on a tiny runtime-support header.

The same contraction primitives are exposed directly in C++, both for hand-built traces and
as the numeric oracle the generated kernels are validated against.

## Why generate, rather than evaluate symbolically

Doing the full tensor contraction inside the consumer's compiler (as an earlier symbolic
expression-tensor engine did) is correct but expensive: a trace of several transverse
projectors expands explosively in the frame-component basis, and the compiler never reclaims
intermediate memory. Running the *same* contraction as a build-time step — numerically, over a
fixed frame — sidesteps both: the generator runs in seconds and tens of megabytes, and the
consumer only ever compiles the small, flat result. This is exactly how FORM kernels are
produced; NumTracer just does the trace in C++ instead.

Next, [Key concepts](concepts.md) builds the mental model the
[tutorials](../tutorials/index.md) rely on. See the [C++ API](../doxygen/NumTracer/html/index)
for the full reference.
