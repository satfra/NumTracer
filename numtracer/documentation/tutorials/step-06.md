# step-06: Your first generated kernel

*Builds on: [step-05](step-05.md) · Built on by: [step-07](step-07.md), [step-08](step-08.md),
[step-09](step-09.md) · Tags: `mathematica`, `codegen` · **Tier B** (a Wolfram kernel and FunKit,
in addition to a C++ compiler)*

```{admonition} What this step needs
:class: important
A Wolfram kernel (`wolframscript`, `wolfram` or `math`) **and** [FunKit](https://github.com/satfra/FunKit).
NumTracer owns the *tensor* half of code generation and delegates the *scalar* half — expression
formatting, CSE of the dressing coefficients, function/class/header boilerplate — to FunKit's COEN
emitter, so `MakeNTKernel` will not run without it.

You do **not** need TensorBases, FORM, or DiFfRG. This script builds its network by hand.
```

## Introduction

Steps 01–05 did everything a kernel needs, by hand: choose a frame, describe the Dirac chain,
describe the Lorentz network, contract, lower, emit. Every one of those actions was mechanical, and
every one of them is *derivable from the flow equation you started with*. So the front-end does
them for you.

This step generates a kernel for exactly the network [step-04](step-04.md) contracted — deliberately,
so that you can compare the emitted arithmetic against a result you already understand line by line.

### What "generating" actually means

`MakeNTKernel` does something slightly unusual, and it is worth being clear about because it
explains most of the option surface later:

1. It writes a small **C++ generator program** to disk.
2. It **compiles** that program (with `clang++` or `g++`, in parallel translation units).
3. It **runs** it. *This* is where the contraction of steps 01–05 happens — in compiled C++, not in
   Mathematica.
4. The generator writes the two committed headers: the lowered `trN()` traces, and the kernel class
   that fills the frame symbols and calls them.

So the Wolfram layer is an *analyser and emitter*, not a computer algebra engine doing the trace.
That split is the whole performance story: symbolic tensor algebra in Mathematica would be
hopelessly slow on a four-point vertex, while the same contraction as `constexpr`-ish C++ over a
compact frame takes seconds ([step-20](step-20.md)).

```{admonition} Why not do the contraction in the consumer's compiler?
:class: note
It was tried. Doing the contraction as C++ template metaprogramming in the consumer's build means
the compiler never reclaims the intermediate memory it allocates during constant evaluation, so
peak RAM grows with *total* allocations rather than the live set. Generating at build time and
shipping flat arithmetic sidesteps it entirely — see
[Under the Hood](../internals/index.md).
```

## The commented program

`Tutorials/step-06-first-kernel/first_kernel.wls`.

### Loading the front-end

```{literalinclude} ../../../Tutorials/step-06-first-kernel/first_kernel.wls
:language: mathematica
:start-after: "@snip begin: load"
:end-before: "@snip end: load"
```

`Needs["NumTracer`"]` resolves anywhere once the library is installed — the CMake install puts the
paclet on Wolfram's application path (disable with `-DNUMTRACER_INSTALL_MATHEMATICA=OFF`). The
`NUMTRACER_DIR` fallback is what lets the same script run from a bare checkout, and it is what the
tutorial CMake project sets.

`Global`Nc = 3` looks gratuitous for a network with no colour in it, and in a sense it is — but
every SU($N$) head carries its rank as a leading argument, `NumTrace` validates all of them up
front, and the default rank is read from `Nc`. Setting it costs nothing and prevents a confusing
abort the first time you add a `ntSUNT`.

### The frame

```{literalinclude} ../../../Tutorials/step-06-first-kernel/first_kernel.wls
:language: mathematica
:start-after: "@snip begin: frame"
:end-before: "@snip end: frame"
```

`propFrame[p, l1, cos1, qp, ql]` is the same one-angle kinematics [step-04](step-04.md) wrote out by
hand: `qp` along axis 0 with magnitude `p`, `ql` at angle `cos1` with magnitude `l1`. The difference
is that you no longer write components — you name the kinematics and the frame supplies them.

The first line binds the names to symbols in the `Global`` context explicitly. In a plain script
that is what `p` would mean anyway; spelling it out keeps the script correct when it is loaded from
inside another context (a package, a notebook with its own context). It matters because the code
generator recognises the one-angle propagator frame by these symbols — the external magnitude
`Global`p` and the loop magnitude `Global`l1` — and takes a faster path for it.

`"Args"` lists the same three scalars in a *different* order from `propFrame`, and on purpose:
`propFrame`'s arguments describe the frame, while `"Args"` is the signature of the generated kernel,
which follows the integrator's convention — integration variables (`l1`, `cos1`) first, external
momenta (`p`) after.

**The two kinds of name matter.** `qp` and `ql` are *momentum tags*: they appear inside the network
and the frame resolves them to components. `p`, `l1`, `cos1` are *scalars*: they appear in the
generated C++ signature, in the order given by `"Args"`. Mixing the two up — putting a momentum tag
in `"Args"`, or a scalar in the network — is the classic first-script error.
[step-07](step-07.md) is entirely about this layer.

### The network

```{literalinclude} ../../../Tutorials/step-06-first-kernel/first_kernel.wls
:language: mathematica
:start-after: "@snip begin: net"
:end-before: "@snip end: net"
```

Compare this one line with the C++ of [step-04](step-04.md). It is the same network: two `ntVec`
factors on labels `mu` and `nu`, an `ntTransProj` between them, and repeated labels summed. The DSL
is deliberately close to the FORM/FormTracer vocabulary so that a flow written for either reads the
same.

`ntSP[qp, qp]` is a **scalar coefficient**, not a tensor factor. Scalars multiply a diagram; they
never contract. The frame resolves `ntSP[qp, qp]` to $p^2$, and dividing by it is what makes the
answer the pure angular factor. The full head vocabulary is [step-08](step-08.md).

### The two entry points

```{literalinclude} ../../../Tutorials/step-06-first-kernel/first_kernel.wls
:language: mathematica
:start-after: "@snip begin: emit"
:end-before: "@snip end: emit"
```

`NumTrace` **analyses** — it produces an `NTKernel`, which is a list of diagrams (each a scalar
coefficient times a contraction), the frame, and the layout of the runtime array the kernel fills. Nothing has been contracted at
this point. Printing `Length[ntk[[1]]["Diagrams"]]` is the cheapest possible sanity check and worth
doing in every script: a network you expected to have 40 diagrams that reports 0 has been pruned by
a rule you did not intend.

`MakeNTKernel` **emits**. Three output files, always — the signature enforces it. If you pass two,
you get `MakeNTKernel::nfiles` and an abort.

## Results

```bash
cmake --build build --target first_kernel && ./build/first_kernel
```

```{literalinclude} ../../../Tutorials/step-06-first-kernel/first_kernel.expected.txt
:language: text
```

### Reading what was generated

This is the part worth lingering on. The traces header, as generated (a snapshot committed next to
the script; the tutorial test checks that a fresh generation still matches it):

```{literalinclude} ../../../Tutorials/step-06-first-kernel/snapshot/first_kernels.hh
:language: cpp
```

`tr0` computes $p^2(1 - \cos^2\theta)$ — **exactly the two-monomial polynomial
[step-04](step-04.md) produced**, lowered by the same Horner + CSE passes [step-05](step-05.md)
printed. Nothing new happened in the engine; the front-end merely wrote the driver.

Notice `nenv = 2`. The frame has three scalar arguments, but the contraction only ever needed two of
them — `l1` cancelled out of this particular network, so no env slot was allocated for it. The
`fill` signature still takes all three, because the signature is the *flow's* contract, not a
minimal set. (This is why generated headers are compiled with `-Wno-unused-parameter`; an unused
argument is normal, not a defect.)

And the kernel class (abridged to its two entry points; the full file is
`Tutorials/step-06-first-kernel/snapshot/first_kernel.hh`):

```{literalinclude} ../../../Tutorials/step-06-first-kernel/snapshot/first_kernel.hh
:language: cpp
:start-at: "static inline auto kernel"
```

Three things to take from this shape:

* **The signature is the `"Args"` list, in order.** That ordering is a contract with whatever calls
  the kernel, and [step-15](step-15.md) shows what happens in a real project when the two drift.
* **The $1/p^2$ became `powr<-2>(p)`, outside the trace.** The scalar coefficient is emitted as
  ordinary arithmetic multiplying the trace, exactly the
  $\sum_{\text{diagrams}} \text{coeff} \times \text{trace}$ shape every flow has. Scalars never
  enter the tensor algebra.
* **`constant()` exists and returns 0.** It is the loop-independent term a flow adds flat to the
  integral. This network has none; [step-13](step-13.md) is the case where it is the whole point.

### Self-contained by default

The header includes only two small NumTracer headers — `numtracer/codegen/runtime.hpp` (`powr`,
`pow`, `sqrt`, `fma`, a small complex type) and `numtracer/sun/sun_data.hpp` (the SU(2)/SU(3)
tables) — and is wrapped in the neutral namespace `numtracer_kernels`. A consumer with
its own support API redirects all of that with the `"RuntimeInclude"`, `"SupportNamespace"` and
`"KernelNamespace"` options — which is how the in-repo test kernels emit against the test shim, and
how [step-15](step-15.md) emits against DiFfRG.

## Possibilities for extensions

1. **Make `l1` matter.** Drop the `1/ntSP[qp, qp]` coefficient and regenerate. Predict what `nenv`
   becomes and what `tr0` returns before you look. Then replace the projector's momentum with a
   combination, e.g. contract against $P(ql - qp)$, and watch `nenv` grow to three.

2. **Add a second diagram.** A network that is a `Plus` at top level is a *sum of diagrams*:
   ```mathematica
   net = (1/ntSP[qp,qp]) ntVec[qp,mu] ntTransProj[ql,mu,nu] ntVec[qp,nu]
       + (1/ntSP[qp,qp]) ntVec[qp,mu] ntLongProj[ql,mu,nu] ntVec[qp,nu];
   ```
   The answer must be identically 1. Check both the value and how many `trN` functions were
   emitted — and then think about why the generator did not simply fold the two into one.

3. **Add a dressing.** Pass `"Dressings" -> {ZA}` and multiply the network by `ZA[l1]`. Look at how
   the kernel signature changes, and at where `ZA` appears in the body. Dressings are opaque
   callables the consumer supplies; the engine never looks inside one.

4. **Break it deliberately — a scalar in the network.** Write `ntVec[p, mu]` instead of
   `ntVec[qp, mu]`, i.e. put the *scalar* `p` where a momentum tag belongs. The frame has no entry
   for `p`, so it cannot be resolved. Read the error; then do the reverse and put `qp` in `"Args"`.
   These two mistakes account for most first-script failures.

5. **Break it deliberately — an unbalanced label.** Change the second `ntVec`'s label from `nu` to
   `rho`. The network no longer closes. Compare the message with the C++-side failure from
   [step-04](step-04.md) extension 5.

6. **Look at the generator.** Open `gen_first_kernel.cpp` in the build directory. It is a readable
   C++ program that sets up a `Frame`, builds the very `LorentzNet`/`DiracChain` objects you wrote
   by hand in Group I, traces them with `frame.trace`, and prints the header. Seeing that the front-end's output is *just the code you
   already know how to write* removes most of the mystery from the rest of this group.

## The plain program

```{literalinclude} ../../../Tutorials/step-06-first-kernel/first_kernel.wls
:language: mathematica
```

```{literalinclude} ../../../Tutorials/step-06-first-kernel/check_first_kernel.cpp
:language: cpp
```

Next, [step-07](step-07.md): frames in earnest — the frame builders, multi-angle kinematics, and the
one that is silently degenerate at the symmetric point.
