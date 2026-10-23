# Tutorials

The tutorials are a **numbered series**. Each step is a complete, runnable, self-checking program
plus the narrative that explains it, and each builds on concepts introduced earlier. Step 1 assumes
nothing but a C++ compiler; by step 22 you are generating GPU kernels for four-point vertex flows
and validating them against independent oracles.

Every step follows the same four-part shape:

1. **Introduction** — the problem, the mathematics, and *why* the API is shaped the way it is.
2. **The commented program** — the real source, walked through chunk by chunk.
3. **Results** — the program's actual output, and what to notice in it.
4. **Possibilities for extensions** — graduated exercises, including at least one deliberate
   breakage so you learn to read the engine's error messages.

```{admonition} New here? Read these two pages first
:class: tip
[**Key concepts**](../getting_started/concepts.md) defines how a network is described — indices
contract when they share a label, the frame fixes each momentum's components, scalar products
become the kernel's symbols. [**Scope & conventions**](../getting_started/scope-and-conventions.md)
states the hard assumptions (4 dimensions, **Euclidean** metric, chiral-basis gammas).
```

```{admonition} Not here for QCD?
:class: note
The worked examples are quantum-field-theory loop integrands, because that is what the authors
needed traced. The engine itself knows no physics — it contracts indices. Steps 1–9 use almost no
physics vocabulary; [**Bring your own network**](../getting_started/bring-your-own-network.md) is
the domain-neutral on-ramp, and the [**Glossary**](../getting_started/glossary.md) marks which
terms are *engine* concepts and which are merely the example's vocabulary.
```

## Connections between the tutorial programs

Each arrow means "introduces a concept the later step relies on". Nothing forces you to read in
order — pick a target and walk back along its arrows.

```{graphviz}
:alt: How the tutorial steps build on one another

digraph tutorials {
  bgcolor="transparent";
  rankdir=TB;
  node [shape=box, style="rounded,filled", fontname="Helvetica", fontsize=10,
        target="_top", penwidth=0.8, height=0.34];
  edge [color="#888888", arrowsize=0.7, penwidth=0.9];

  subgraph cluster_1 {
    label="I · the engine, by hand"; fontname="Helvetica"; fontsize=10;
    color="#cccccc"; style=dashed;
    node [fillcolor="#dbeafe", color="#60a5fa"];
    s01 [label="1 · hello network",     URL="../tutorials/step-01.html"];
    s02 [label="2 · colour factors",    URL="../tutorials/step-02.html"];
    s03 [label="3 · Dirac traces",      URL="../tutorials/step-03.html"];
    s04 [label="4 · Lorentz networks",  URL="../tutorials/step-04.html"];
    s05 [label="5 · a full diagram",    URL="../tutorials/step-05.html"];
  }

  subgraph cluster_2 {
    label="II · the front-end"; fontname="Helvetica"; fontsize=10;
    color="#cccccc"; style=dashed;
    node [fillcolor="#dcfce7", color="#4ade80"];
    s06 [label="6 · first generated kernel", URL="../tutorials/step-06.html"];
    s07 [label="7 · frames & kinematics",    URL="../tutorials/step-07.html"];
    s08 [label="8 · the DSL, head by head",  URL="../tutorials/step-08.html"];
    s09 [label="9 · codegen options",        URL="../tutorials/step-09.html"];
  }

  subgraph cluster_3 {
    label="III · physics workflows"; fontname="Helvetica"; fontsize=10;
    color="#cccccc"; style=dashed;
    node [fillcolor="#fef3c7", color="#fbbf24"];
    s10 [label="10 · the Wetterich equation", URL="../tutorials/step-10.html"];
    s11 [label="11 · projectors & bases",     URL="../tutorials/step-11.html"];
    s12 [label="12 · a propagator flow",      URL="../tutorials/step-12.html"];
    s13 [label="13 · counterterms",           URL="../tutorials/step-13.html"];
    s14 [label="14 · 3- and 4-point vertices",URL="../tutorials/step-14.html"];
    s15 [label="15 · wiring into DiFfRG",     URL="../tutorials/step-15.html"];
  }

  subgraph cluster_4 {
    label="IV · advanced sectors"; fontname="Helvetica"; fontsize=10;
    color="#cccccc"; style=dashed;
    node [fillcolor="#fae8ff", color="#e879f9"];
    s16 [label="16 · finite temperature",   URL="../tutorials/step-16.html"];
    s17 [label="17 · flavour dressings",    URL="../tutorials/step-17.html"];
    s18 [label="18 · dressed numerators",   URL="../tutorials/step-18.html"];
    s19 [label="19 · four-fermi & Fierz",   URL="../tutorials/step-19.html"];
  }

  subgraph cluster_5 {
    label="V · scale and trust"; fontname="Helvetica"; fontsize=10;
    color="#cccccc"; style=dashed;
    node [fillcolor="#fee2e2", color="#f87171"];
    s20 [label="20 · making it generate", URL="../tutorials/step-20.html"];
    s21 [label="21 · GPU & complex",      URL="../tutorials/step-21.html"];
    s22 [label="22 · validate & debug",   URL="../tutorials/step-22.html"];
  }

  s01 -> s02; s01 -> s03; s01 -> s04;
  s02 -> s05; s03 -> s05; s04 -> s05;
  s02 -> s17; s03 -> s18; s04 -> s16; s05 -> s18;
  s05 -> s06;
  s06 -> s07; s06 -> s08; s06 -> s09;
  s07 -> s10; s07 -> s16;
  s08 -> s11; s08 -> s19;
  s09 -> s20; s09 -> s21;
  s10 -> s11; s10 -> s12; s11 -> s12; s11 -> s14;
  s12 -> s13; s13 -> s14;
  s14 -> s15; s14 -> s19;
  s15 -> s20; s15 -> s22; s20 -> s22;
}
```

## The steps

| # | Step | What it teaches | Tags |
|---|---|---|---|
| 1 | [Hello, tensor network](step-01.md) | Indices contract when labels match; a closed network is a number. | `network` `projector` |
| 2 | [Colour factors fold to a number](step-02.md) | SU($N$) contraction without materialising the tensor: $C_F = 4/3$. | `sun` `colour` |
| 3 | [Dirac traces: matrices, then tokens](step-03.md) | A γ-trace is a 4×4 matrix product — twice: by hand, then via the token API. | `dirac` `trace` |
| 4 | [Lorentz networks and inverse atoms](step-04.md) | Projectors, and how $1/l^2$ survives contraction as an *atom*. | `lorentz` `projector` |
| 5 | [A full diagram, and what lowering does to it](step-05.md) | Assemble a whole numerator, validate it, then watch it become straight-line C++. | `capstone` `codegen` `cse` |
| 6 | [Your first generated kernel](step-06.md) | `NumTrace` + `MakeNTKernel`: the front-end writes steps 1–5 for you. | `mathematica` `codegen` |
| 7 | [Frames and kinematics](step-07.md) | What a frame *is*; the five builders; when the symmetric point is degenerate. | `frame` `kinematics` |
| 8 | [The DSL, head by head](step-08.md) | The complete network vocabulary, each head with a tiny net and its closed form. | `dsl` `reference` |
| 9 | [Codegen options: what each one changes](step-09.md) | The `MakeNTKernel` option surface, shown by diffing the emitted C++. | `codegen` `options` |
| 10 | [From the Wetterich equation to a network](step-10.md) | The FunKit path: derivatives, routing, dressing rules, `FromFunKit`. | `funkit` `frg` |
| 11 | [Projectors and tensor bases](step-11.md) | `TBGetProjector` vs `TBGetVertex`, index order, signs — and how to assert them. | `tensorbases` `projector` |
| 12 | [A propagator flow, end to end](step-12.md) | The ghost and gluon two-point flows; symmetries and when reduction is legal. | `frg` `flow` |
| 13 | [Counterterms and the constant slot](step-13.md) | The RG counterterm that pins a dressing, and the Newton loop that consumes it. | `frg` `diffrg` |
| 14 | [Three- and four-point vertices](step-14.md) | Multi-angle frames, `AngleDefs`, basis restriction, one derivative many projectors. | `vertex` `frame` |
| 15 | [Wiring a flow into DiFfRG](step-15.md) | `MakeNTKernelDiFfRG`, the `flows/` tree, the parameter-order contract, `model.hh`. | `diffrg` `build` |
| 16 | [Finite temperature](step-16.md) | The heat-bath split: electric/magnetic projectors, spatial products, Matsubara. | `finite-T` `projector` |
| 17 | [Per-flavour and per-component dressings](step-17.md) | Dressing u and d differently *inside* one SU($N$) trace. | `sun` `dressing` |
| 18 | [Dressed numerators and vertex collection](step-18.md) | One collected trace instead of $2^D$ diagrams. | `dressing` `dressed-poly` |
| 19 | [Four-fermi, Fierz bases, fixed-component γ](step-19.md) | Two independent spinor lines, contact structures, γ⁰ as a unit-vector slash. | `four-fermi` `fierz` |
| 20 | [Making a hard flow generate](step-20.md) | Phase A/B, thread caps, RAM, chunking, dedup, offline generation. | `performance` `build` |
| 21 | [GPU kernels and complex traces](step-21.md) | Device decorators, the two-phase integrator, the real-probe verdict. | `gpu` `cuda` `complex` |
| 22 | [Validating and debugging a kernel](step-22.md) | Oracles, the codegen blind spot in `ctest`, and how to read a failure. | `testing` `debugging` |

## By topic

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} I want to understand the engine
Steps [1](step-01.md), [3](step-03.md), [4](step-04.md), [5](step-05.md) — then
[Under the Hood](../internals/index.md).
:::

:::{grid-item-card} I want to generate a kernel
Steps [6](step-06.md), [7](step-07.md), [8](step-08.md), [9](step-09.md).
:::

:::{grid-item-card} I want to add a flow to my fRG project
Steps [10](step-10.md), [11](step-11.md), [12](step-12.md), [13](step-13.md),
[14](step-14.md), [15](step-15.md).
:::

:::{grid-item-card} My generation is too slow, or it OOMs
Step [20](step-20.md) first; then [21](step-21.md) if you are targeting a GPU.
:::

:::{grid-item-card} My kernel gives the wrong number
Step [22](step-22.md), then [11](step-11.md) — most "engine bugs" turn out to be
projector conventions.
:::

:::{grid-item-card} I need a specific sector
Colour/flavour [2](step-02.md), [17](step-17.md) · Dirac [3](step-03.md),
[18](step-18.md), [19](step-19.md) · finite $T$ [16](step-16.md).
:::

::::

## What each step needs

The series is deliberately layered so that the early steps have **no dependencies beyond a C++
compiler**. Every step states its tier at the top.

| Tier | Requirements | Steps |
|---|---|---|
| **A** | A C++20 compiler and CMake. Nothing else. | 1–5, 16, 17, 18 |
| **B** | …plus a Wolfram kernel and [FunKit](https://github.com/satfra/FunKit). The scripts build their networks *by hand*, so no TensorBases, no FORM, no DiFfRG. | 6, 7, 8, 9, 19, 21 |
| **C** | …plus TensorBases and DiFfRG — the rest of the fRG toolchain. | 10–15, 20, 22 |

```{admonition} Why tier B needs FunKit
:class: note
NumTracer owns the *tensor* half of code generation. The *scalar* half — expression formatting,
CSE of the dressing coefficients, and all the function/class/header boilerplate — is delegated to
FunKit's mature COEN emitter rather than reimplemented. So `NumTrace` will analyse a network
without FunKit, but `MakeNTKernel` cannot write a kernel without it. **FORM is never needed** for
your own kernels; it enters only when regenerating the project's committed reference flows, where
FunKit uses it to compute tensor *bases* — never to do NumTracer's trace.
```

## Building and running the programs

The programs live in the repository's `Tutorials/` folder as a **standalone CMake project** that
consumes NumTracer the way any external user would (`find_package(NumTracer)`, with a fallback to
the in-tree library so it also works straight from a checkout):

```bash
cd Tutorials
cmake -S . -B build            # configure once
cmake --build build -j4        # build every tier-A tutorial
ctest --test-dir build         # run them all as a self-checking test suite
```

Tier-B steps additionally generate their kernels through Wolfram. They are wired into the same
project but only when a Wolfram kernel is on `PATH`:

```bash
cmake --build build --target tutorials_codegen   # run every codegen tutorial
ctest --test-dir build                           # now includes the tier-B checks
```

Each page also gives the one-target command (`cmake --build build --target <name>`) and the
program's real output.

```{admonition} The code on these pages is the code that runs
:class: important
Every listing is pulled directly out of the file in `Tutorials/` at documentation-build time — the
pages cannot drift from the programs `ctest` executes. If a listing looks wrong, the program is
wrong, and `ctest` will say so.
```

For the design behind the engine, see [Under the Hood](../internals/index.md).

```{toctree}
:hidden:
:maxdepth: 1

step-01
step-02
step-03
step-04
step-05
step-06
step-07
step-08
step-09
step-10
step-11
step-12
step-13
step-14
step-15
step-16
step-17
step-18
step-19
step-20
step-21
step-22
```
