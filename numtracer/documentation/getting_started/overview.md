# Overview

NumTracer assembles and evaluates the tensor networks that appear in quantum-field-theory
loop integrands, and **generates the scalar kernels** that a numerical integrator then runs.
It does the tensor algebra in plain compile-time C++ — no symbolic-algebra runtime, no external
tools — and emits a flat, straight-line real-arithmetic kernel that a consumer compiles directly.

## What it computes

A loop integrand is a product of tensors — gamma matrices, momentum vectors, metrics,
projectors, SU($N$) generators — contracted together (summed over shared indices) down to a
single number that still depends on a few runtime quantities: momentum magnitudes, angles, and
propagator dressings. That number is fed to an integrator and evaluated at hundreds of
thousands of grid points.

NumTracer is a *general* tensor-tracing engine. The physics is entirely in the network you
hand it; the engine only contracts indices and folds the result. The numeric path handles
Lorentz and Dirac structure and any number of SU($N$) groups (colour, flavour, …), dresses
individual SU($N$) components differently (per-component split and group-diagonal dressings), and factors
a diagram that disconnects into several closed traces into a product — see
[Key concepts](concepts.md) and the [step-17](../tutorials/step-17.md).

```{admonition} Is this for me?
:class: tip
If you have a tensor network built from metrics, vectors, projectors, gamma matrices, and/or
SU($N$) factors — in a **4D Euclidean** setting — NumTracer will contract it to a scalar kernel,
whatever theory it comes from. The worked flows here are fRG flows for Yang–Mills and QCD, but that
is just the authors' application. Before assuming it is (or isn't) for you, read
[Scope & conventions](scope-and-conventions.md) for the exact boundary, then
[Bring your own network](bring-your-own-network.md) for a domain-neutral hello-world and the
dictionary that maps your objects onto the engine's.
```

## What goes in, what comes out

There are two ways in, and both end in the same contraction engine.

**From C++** (no other tools needed). You describe a network directly — a Dirac chain such as
$\mathrm{tr}(\slashed p\,\gamma^\mu \slashed q\,\gamma^\nu)$, a Lorentz network such as a
projector $P_T(l)_{\mu\nu}$, an SU($N$) factor such as $T^a_{ij}T^a_{ji}$ — on a *frame* that says
what the momenta are. NumTracer contracts it and hands back a polynomial in the frame's symbols,
which you can evaluate, or lower to a straight-line C++ function. Tutorials
[1–5](../tutorials/step-01.md) do exactly this.

**From Mathematica** (needs Wolfram and FunKit). You write the network in a small DSL (or import a
FunKit flow), and `MakeNTKernel` writes a complete, self-contained C++ kernel: one function per
trace, a `fill()` that computes the frame symbols once per call, and the assembly of all diagrams.
[Tutorial 6](../tutorials/step-06.md) is the first one.

In both cases the result does not depend on NumTracer at run time: a generated kernel is plain C++.

## Why generate, rather than evaluate symbolically

Doing the full tensor contraction inside the consumer's compiler is correct but expensive: a
trace of several transverse projectors expands explosively in the frame-component basis, and the
compiler never reclaims intermediate memory. Running the *same* contraction as a build-time step —
numerically, over a fixed frame — sidesteps both: the generator runs in seconds and tens of
megabytes, and the consumer only ever compiles the small, flat result. How the engine avoids the
combinatorial blow-ups is described in [Under the Hood](../internals/index.md).

## Where to go next

- Coming from FORM? [Coming from FORM](coming-from-form.md) maps the vocabulary in one table.
- [Key concepts](concepts.md) builds the mental model the [tutorials](../tutorials/index.md) rely on.
- [Scope & conventions](scope-and-conventions.md) says exactly what is supported (4D, Euclidean, SU($N$)).
- The [C++ API](../doxygen/NumTracer/html/index) is the full reference.
