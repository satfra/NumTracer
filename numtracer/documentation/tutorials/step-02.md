# step-02: Colour factors fold to a number

*Builds on: [step-01](step-01.md) · Built on by: [step-05](step-05.md), [step-17](step-17.md) ·
Tags: `sun`, `colour` · **Tier A** (a C++20 compiler, nothing else)*

## Introduction

[step-01](step-01.md) contracted Lorentz indices. This step does the same thing in a different
sector — SU($N$), here as QCD colour — and the point is that the rule does not change: indices
contract by label. What is new is that SU($N$) has two kinds of index, adjoint (extent $N^2-1$) and
fundamental (extent $N$), and the labels carry that kind in their type.

What *is* new is how the contraction is carried out, and it is worth understanding because it is
the reason the library exists.

### The physics, in one line

The quark self-energy exchanges one gluon. The quark line picks up a generator $T^a$ where the
gluon is emitted and another where it is absorbed, with the gluon's colour index $a$ shared. The
colour structure is therefore

$$
(T^a T^a)_{ik} \;=\; C_F\,\delta_{ik},
\qquad
C_F \;=\; \frac{N^2-1}{2N} \;=\; \frac{4}{3} \ \ \text{for } N = 3 .
$$

$C_F$ is the *Casimir* — a pure number that multiplies the whole diagram. Every one-gluon-exchange
diagram in QCD carries it. We will compute it by closing the quark line into a loop, which turns
the free indices $i,k$ into a trace:

$$
\mathrm{tr}(T^aT^a) \;=\; C_F\,\mathrm{tr}\,\delta \;=\; C_F\,N \;=\; \frac{N^2-1}{2} \;=\; 4 .
$$

We also do a purely gluonic one, $f^{abc}f^{abc} = N(N^2-1) = 24$, which has no free indices at
all.

### Why not just build the tensor?

You could form $T^a$ as three $3\times3$ complex matrices, multiply, and trace. That is 8 matrices
of 9 entries — trivial. Now consider a four-gluon vertex flow, whose colour structure is a network
of half a dozen $f^{abc}$'s over adjoint indices running $1..8$. Materialising the intermediate
tensor means $8^6 = 262\,144$ complex entries, almost all of which are zero, and the answer is a
single number. The dense object is built only to be summed straight back down.

NumTracer instead **folds** the network: it walks the factors, contracting one index at a time
against the typed-out sparse tables in `sun/sun_data.hpp`, and never allocates an intermediate
larger than it needs. The generators are sparse ($\lambda_1$ has 2 nonzero entries out of 9), and
$f^{abc}$ is totally antisymmetric with only 9 independent nonzero components out of 512 — the
fold exploits both.

```{admonition} Two folds, two return types
:class: note
`su3.value(net)` returns a single `Cx` — the case for every ordinary colour factor. A network mixing
several groups (colour SU(3) ⊗ flavour SU(2)) is folded by the free function `nt::sun_value(net)`;
the value factorises into the per-group values. A second entry point, `nt::sun_value_dressed(net)`,
returns a *polynomial* over per-component dressings; that is [step-17](step-17.md), and you do not
need it yet.
```

## The commented program

`Tutorials/step-02-colour-factors/colour_factors.cpp`.

### The group and its labels

```{literalinclude} ../../../Tutorials/step-02-colour-factors/colour_factors.cpp
:language: cpp
:start-after: "@snip begin: labels"
:end-before: "@snip end: labels"
```

`nt::SUN su3(3)` is the group. It hands out labels of two types: `su3.adjoint<3>()` gives three
`AdjIndex` (extent $N^2-1 = 8$), `su3.fundamental<2>()` two `FundIndex` (extent $N = 3$). Every
factor builder takes the label types its slots need — `T(AdjIndex, FundIndex, FundIndex)` — so an
adjoint label in a fundamental slot is a compile error, not a wrong number.

Each `SUN` object is its own group. A QCD flow with two light flavours has an SU(3) colour group
*and* an SU(2) isospin group in the same network; make one object for each. Their labels never
contract with each other, and using one group's label in the other's factor throws.

### The Casimir

```{literalinclude} ../../../Tutorials/step-02-colour-factors/colour_factors.cpp
:language: cpp
:start-after: "@snip begin: cf"
:end-before: "@snip end: cf"
```

`*` multiplies factors into a network, `su3.value` folds it. Read the index pattern:
`su3.T(a, i, j)` then `su3.T(a, j, i)`. The adjoint label `a` is shared, so the gluon index is
summed — that is the exchange. The fundamental labels go $i \to j$ then $j \to i$, closing the
quark line into a loop. A closed loop with no free indices is a number.

### A purely gluonic one

```{literalinclude} ../../../Tutorials/step-02-colour-factors/colour_factors.cpp
:language: cpp
:start-after: "@snip begin: ff"
:end-before: "@snip end: ff"
```

Two structure constants sharing all three adjoint labels. There is nothing left free, so this is
already a scalar with no closing delta needed. This is the colour factor of the ghost loop in the
gluon self-energy, and of the gluon loop's leading piece.

## Results

```bash
cmake --build build --target colour_factors && ./build/colour_factors
```

```{literalinclude} ../../../Tutorials/step-02-colour-factors/colour_factors.expected.txt
:language: text
```

All three are *exact* — not "agrees to $10^{-12}$", but exactly the integers, because the fold sums
rationals stored as `Cx` and nothing here introduces a rounding error. That is a useful property to
know about: when a colour factor comes back as `1.3333333333333335` you have a bug in the network,
not a precision problem.

Note also what the second line is doing. $C_F$ is *defined* by $T^aT^a = C_F\,\delta$, so extracting
it from the trace requires dividing by $\mathrm{tr}\,\delta = N$. The engine gave us the trace; the
interpretation is ours. This division-by-$N$ step is a classic source of factor-of-3 errors in real
flows, and it is why [step-11](step-11.md) insists on asserting projector normalisations rather
than trusting them.

## Possibilities for extensions

1. **The adjoint Casimir.** $f^{acd}f^{bcd} = N\,\delta^{ab} = C_A\,\delta^{ab}$. Build it with
   `a`, `b` free and `c`, `d` shared, close it with `su3.delta(b, a)`, and check you get
   $N(N^2-1) = 24$ — the same number as `ff` above, which is not a coincidence. Then explain why.

2. **Change the group.** Make an `nt::SUN su2(2)` and recompute both quantities. Predict first:
   $C_F = 3/4$, $\mathrm{tr}(T^aT^a) = 3/2$, $f^{abc}f^{abc} = 6$. SU(2) is the isospin group
   [step-17](step-17.md) uses.

3. **The fundamental identity.** $T^a_{ij}T^a_{kl} = \tfrac12(\delta_{il}\delta_{kj} -
   \tfrac1N \delta_{ij}\delta_{kl})$ is the Fierz identity that underlies every four-quark
   calculation. Verify it by contracting both sides against a fixed pair of deltas and comparing.
   This is a miniature of what [step-19](step-19.md) does with a whole basis.

4. **Break it deliberately — mix the sectors.** Give a fundamental slot an adjoint label,
   `su3.T(a, a, j)`. It does not compile: `AdjIndex` is not a `FundIndex`. Now make a second group
   `nt::SUN su2(2)` with its own `auto [x] = su2.adjoint<1>();` and write `su3.T(x, i, j)`. That
   compiles — both are adjoint labels — and throws at run time, because `x` belongs to another
   group.

5. **Feel the blow-up you are avoiding.** Contract a chain of six $f$'s (a ring:
   $f^{a b c}f^{c d e}f^{e f g}\dots$ closing back on $a$) and time it. Then estimate how many
   entries the dense intermediate would have had.

## The plain program

```{literalinclude} ../../../Tutorials/step-02-colour-factors/colour_factors.cpp
:language: cpp
```

Next, [step-03](step-03.md): the Dirac sector, where the same "fold, don't materialise" idea meets
a genuinely explosive alternative.
