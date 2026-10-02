# step-05: A full diagram, and what lowering does to it

*Builds on: [step-02](step-02.md), [step-03](step-03.md), [step-04](step-04.md) · Built on by:
[step-06](step-06.md) · Tags: `capstone`, `codegen`, `cse` · **Tier A** (a C++20 compiler, nothing
else)*

## Introduction

This is the capstone of Group I. Everything so far has been one sector at a time; here they meet,
and then the result is pushed one stage further — through the **lowering** pass that turns a
polynomial into straight-line C++. By the end of this page you will have seen, in one program, the
entire path a generated kernel takes.

### The diagram

The one-loop quark self-energy: a quark of momentum $p$ emits a gluon, propagates with momentum
$q$, and reabsorbs it. Projected onto the wave-function structure, the Dirac numerator is

$$
T_{\text{num}} \;=\; \mathrm{tr}\!\big[\slashed p\,\gamma^\mu\,\slashed q\,\gamma^\nu\big]\,
P_{\mu\nu}(l), \qquad q = l - p .
$$

Read it as the composition of the two previous steps. The trace is a Dirac chain
([step-03](step-03.md)) with two open legs $\mu,\nu$ — the two ends of the gluon line. The
projector is a Lorentz network ([step-04](step-04.md)) that ties those legs off, carrying the
gluon's transversality and its $1/l^2$. The colour factor $C_F$ ([step-02](step-02.md)) multiplies
the whole thing and is a separate number; it shares no index with the rest.

That factorisation — Dirac chain × Lorentz net × colour number × scalar dressings — is not specific
to this diagram. It is the shape *every* one-loop flow has, and the shape the code generator
assumes.

### The closed form

For $q = l - p$ in the one-angle frame, with $p = |p|$, $l = |l|$ and $c = \cos\theta$:

$$
T_{\text{num}} \;=\; 4\,p\,\big(-3\,c\,l + p + 2\,c^2 p\big).
$$

Having an independent closed form is the point of the exercise: it pins down the engine's algebra,
including the sign conventions, the $\mathrm{tr}\,\mathbb{1} = 4$, and the Euclidean metric, all at
once.

### Lowering

A contracted diagram is a polynomial (`nt::Poly`): a sum of monomials over the frame's symbols and
inverse atoms.
Evaluating it as written would mean recomputing shared subexpressions over and over, for every one
of the hundreds of thousands of quadrature points the integrator visits. So the generator does two
build-time passes:

1. **Horner factoring** (`codegen/lower.hpp`) — rewrite the monomial set as a nest of
   `pivot * (…) + (…)` so that shared factors are computed once.
2. **Real value-numbering / CSE** (`codegen/real_cse.hpp`) — accumulate the factored arithmetic
   into a real SSA form that deduplicates every repeated subexpression and folds away trivial
   operations.

The output is a flat `double f(const double*)` with no branches, no divisions, and no temporaries
beyond named scalars. This program runs both passes, prints the result so you can read the kernel a
real flow would ship — and *runs* it, so you can see that it computes the right number.

## The commented program

`Tutorials/step-05-full-diagram/full_diagram.cpp` computes $T_{\text{num}}$ three ways: the
contraction, the closed form, and the lowered program run in-process.

### The frame

```{literalinclude} ../../../Tutorials/step-05-full-diagram/full_diagram.cpp
:language: cpp
:start-after: "@snip begin: frame"
:end-before: "@snip end: frame"
```

The one-angle frame of [step-04](step-04.md), with symbols $p$, $l$, $\cos\theta$. The internal
momentum is not declared: `q = l - p` is a combination of the two that are.

### The two networks, contracted together

```{literalinclude} ../../../Tutorials/step-05-full-diagram/full_diagram.cpp
:language: cpp
:start-after: "@snip begin: tokens"
:end-before: "@snip end: tokens"
```

The chain and the Lorentz network meet in one call. The engine first closes the Dirac chain into
$4\times4$ products and traces it, leaving a tensor with free indices $\mu$ and $\nu$; then it
contracts that tensor against the projector, summing $\mu$ and $\nu$ away; then it returns the
surviving scalar. The $2^{n_p}$ projector expansion and the $(2n-1)!!$ Wick expansion are both
avoided — neither intermediate is ever formed.

### The closed form

```{literalinclude} ../../../Tutorials/step-05-full-diagram/full_diagram.cpp
:language: cpp
:start-after: "@snip begin: closed"
:end-before: "@snip end: closed"
```

### Lowering to a kernel, and running it

```{literalinclude} ../../../Tutorials/step-05-full-diagram/full_diagram.cpp
:language: cpp
:start-after: "@snip begin: lower"
:end-before: "@snip end: lower"
```

`to_genprog` runs both lowering passes into a `GenProg` — the instruction list of a straight-line
program — over a shared environment `g`. The program reads its inputs from an array `f[]`; `g`
records which slot holds which symbol (or which $1/k^2$). `F.emit_fill` prints the function that
fills `f[]` from the frame's symbols, and `emit_cpp` prints the program itself. Together they are a
complete, self-contained C++ function of $(p, l, \cos\theta)$.

`interpret` executes the same instruction list on the values `F.fill_values` computes at our point.
That is the third check, and the one that tests the lowering.

## Results

```bash
cmake --build build --target full_diagram && ./build/full_diagram
```

```{literalinclude} ../../../Tutorials/step-05-full-diagram/full_diagram.expected.txt
:language: text
```

### Reading the emitted kernel

This is the artefact the whole library exists to produce, so it repays a careful look.

**`fill` is the calling convention.** `f[0..2]` are the three frame symbols. A real kernel's
`fill()` computes every slot from the integrator's arguments once per call; here that is just a
copy. There is no slot for $1/l^2$: the frame knows $l^2$ is the single monomial $l^2$, so every
$1/l^2$ the projector produced cancelled against an $l^2$ in the numerator during contraction. Had
an atom survived, it would have its own slot, filled as `1.0/(…)` — one division per call, and none
in the trace function.

**Every operation appears once.** `s2 = s0*s1` is $p\cos\theta$, computed once; the nesting
`s10 = fma(s0, s6, s9)` with `s6` itself an `fma` is Horner's scheme. The polynomial had a handful
of monomials; the emitted form shares every common factor between them.

**Everything is `const double`.** No arrays beyond `f`, no loops, no branches, no function calls
other than `fma`. This is what "straight-line" means, and it is why the same emitted code compiles
unchanged for a GPU ([step-21](step-21.md)).

**The constants are folded.** Coefficients appear as literals; the overall $-12$ was pulled out by
the factorisation. Nothing is computed at run time that could be computed at build time.

### And the value is right, three ways

The contraction, the closed form, and the *executed* lowered program agree. The third check is the
guarantee that lowering is *value-preserving*: a Horner or CSE bug would show up there and nowhere
else, since the first two paths never touch `to_genprog`. The comparison uses a tolerance rather
than bitwise equality because the emitted program fuses multiply-adds (`fma`), which round
differently from the separate operations of the polynomial evaluation.

## Possibilities for extensions

1. **Add the colour factor.** Multiply by the $C_F = 4/3$ from [step-02](step-02.md) and confirm
   the emitted kernel changes by exactly one constant. In a real flow the colour factor is folded
   into the diagram's scalar coefficient, never into the trace — see if you can see why from the
   emitted code.

2. **The second diagram.** The real quark self-energy has *two* diagrams: the regulator insertion
   sits on the gluon line in one and on the quark line in the other, which amounts to $q = l - p$
   versus $q = l + p$. Redo the contraction with the sign flipped and compare the polynomials. They
   should share most of their structure — which is exactly the redundancy the generator's
   cross-diagram deduplication exploits ([step-20](step-20.md)).

3. **Watch CSE earn its keep.** Replace the projector by a plain metric and re-emit. Then put it back but add a second projector. Count the emitted lines each time. The
   growth is sublinear in the monomial count, and that gap is the whole value of the pass.

4. **Change the frame and re-emit.** Give $l$ a third nonzero component. The value at a fixed point
   must not change; the emitted kernel will grow. This is [step-04](step-04.md) extension 6 seen
   from the other end — you are now looking at what the extra components *cost*.

5. **Break it deliberately — evaluate somewhere singular.** Evaluate $p\cdot P(l)\cdot q$ with a
   second external momentum (so that the $1/l^2$ survives) at $l = 0$. `F.eval` refuses: the
   projector denominator vanishes there. (Before the frame owned the denominators, a hand-supplied
   atom value was an unchecked contract; a wrong one gave a wrong number without complaint.)

6. **Read the lowering source.** With the emitted output in front of you,
   [CSE and Horner lowering](../internals/cse-and-lowering.md) will make considerably more sense
   than it would have before.

## The plain program

```{literalinclude} ../../../Tutorials/step-05-full-diagram/full_diagram.cpp
:language: cpp
```

You have now done by hand everything the code generator does. [step-06](step-06.md) shows the
front-end doing all of it from a one-line description — and the kernel it emits will look extremely
familiar.
