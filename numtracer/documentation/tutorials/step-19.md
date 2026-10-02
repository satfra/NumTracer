# step-19: Four-fermi, Fierz bases, fixed-component gammas

*Builds on: [step-08](step-08.md), [step-14](step-14.md) · Built on by: — ·
Tags: `four-fermi`, `fierz` · **Tier B** (a Wolfram kernel and FunKit)*

## Introduction

Four-quark operators are where several things the engine does quietly all become load-bearing at
once. A four-fermi diagram has **four external legs**, **two independent spinor lines**, structures
with **no momentum dependence at all**, and — because the standard Fierz-complete bases are written
in a finite-$T$ 3+1 split — Dirac structures carrying $\gamma^0$.

That last one is the interesting story, and it is the runnable part of this step.

### $\gamma^0$ is not a gamma with a label named 0

`ntGamma[mu, d1, d2]` takes a Lorentz *label* — a name to be contracted against something else.
`ntGamma[0, d1, d2]` means something categorically different: a gamma matrix with its Lorentz index
**pinned to a concrete component**, the temporal one. It is not waiting to be contracted; it is
already specific.

The DSL originally had no notion of that. The label machinery read the literal `0` as an ordinary
contraction label, found it occurring in many factors at once, and aborted:

```text
NumTrace::privclash: label(s) 0 are private dummies of one factor but also occur outside it
```

**That guard was behaving correctly.** The label bookkeeping was right; the *feature* was missing,
and the standard four-quark bases were simply un-traceable.

### The fix, and why it needed no C++

The resolution is a rewrite inside `NumTrace` (`expandFixedComponents` in `DSL.m`): a fixed
component is a contraction with the **constant unit basis vector** $e_i$,

$$
\gamma^0 \;=\; \gamma^\mu\,(e_0)_\mu ,
$$

so $\gamma^0$ becomes an ordinary **slash** against $e_0$, and $e_0$'s components are injected into
the frame. No new token kind, no new contraction rule, **zero C++ change** — the engine already knew
how to slash a momentum, and $e_0$ is just a momentum whose components happen to be $(1,0,0,0)$.

The head `ntUnitVec[i]` names it. You never write it: the engine introduces it.

The same rewrite applies to a fixed component on a **metric** (`ntMetric[0, nu]`) and on a
$\sigma^{\mu\nu}$ **leg** — for sigma it must convert a "free" leg into a "slash" leg, since a fixed
component is not free.

## The commented program

`Tutorials/step-19-four-fermi/fixed_components.wls`. The frame is `propFrameFT`, and that choice is
not incidental:

> The frame here is `propFrameFT`, whose **TEMPORAL components are nonzero** — at a $T=0$ vacuum
> frame every $q_0$ vanishes and half the oracles below would collapse to 0, i.e. the gate would
> pass vacuously no matter what the tracer did.

### The reference frame

```{literalinclude} ../../../Tutorials/step-19-four-fermi/fixed_components.wls
:language: mathematica
:start-after: "@snip begin: refframe"
:end-before: "@snip end: refframe"
```

This one line is the most valuable idea on the page. The *same* traces can be written **without the
feature at all**, by adding an explicit unit vector as an ordinary user momentum: $\gamma^0$ becomes
$\gamma^\mu (e_0)_\mu$, which is a contraction the engine has always been able to do.

Grading the two spellings against **each other** is convention-independent. An oracle I wrote by
hand can share a sign error with the code; two independent code paths through the engine cannot
agree by accident.

### A. The rewrite happens, and reaches the env

```{literalinclude} ../../../Tutorials/step-19-four-fermi/fixed_components.wls
:language: mathematica
:start-after: "@snip begin: rewrite"
:end-before: "@snip end: rewrite"
```

Checking that `ntUnitVec[0]` is in the env is checking the *mechanism*, not the answer. Without an
env entry the slash emission has nothing to reference — so this catches a rewrite that fired
half-way, which would otherwise show up much later as a confusing emission failure.

### B and C. Component indexing

```{literalinclude} ../../../Tutorials/step-19-four-fermi/fixed_components.wls
:language: mathematica
:start-after: "@snip begin: components"
:end-before: "@snip end: components"
```

Both $\gamma^0$ and $\gamma^1$ are emitted, and that is deliberate. In `propFrameFT` the temporal
component `p0` and the spatial magnitude `p` are **independent kernel arguments**, so an off-by-one
in the unit vector — $e_i$ placed at slot $i+1$, or 1-based components — *swaps the two answers*.
One check alone would not see it.

Assertion C, $\mathrm{tr}(\gamma^0\gamma^0) = 4$, pins something subtler: two fixed components on
one chain must get two **distinct fresh dummies** but the **same env momentum**. Fuse the dummies
and the pair becomes a self-contraction; allocate $e_0$ twice and the env is wrong.

### D and E. The cross-check

```{literalinclude} ../../../Tutorials/step-19-four-fermi/fixed_components.wls
:language: mathematica
:start-after: "@snip begin: crosscheck"
:end-before: "@snip end: crosscheck"
```

### F. Not gamma-specific

```{literalinclude} ../../../Tutorials/step-19-four-fermi/fixed_components.wls
:language: mathematica
:start-after: "@snip begin: metric"
:end-before: "@snip end: metric"
```

## Results

```bash
cmake --build build --target fixed_components && ./build/fixed_components
```

```text
step-19: fixed component rewritten, ntUnitVec[0] present in the env.

fixed-component gammas at (p0, p, l0, l1, cos1) = (0.31, 2.7, 0.83, 1.9, 0.4)
  tr(g^0 q1/)  =         1.24   expected         1.24   4 q1_0 = 4 p0              ok
  tr(g^1 q1/)  =         10.8   expected         10.8   4 q1_1 = 4 p               ok
  tr(g^0 g^0)  =            4   expected            4   4 (e_0.e_0) = 4            ok
  tr(g^0 q1/ g^0 ql/) =      -7.1788   expected      -7.1788   4[2 q1_0 ql_0 - q1.ql]     ok
  g^{0n} q1_n ... =       1.0292   expected       1.0292   4 ql_0 p0                  ok
the convention-independent cross-check
  fixed-component  =      -7.1788
  explicit e_0     =      -7.1788   identical
ALL TESTS PASSED
```

Note the deliberately separated scales: $p_0 = 0.31$ against $p = 2.7$. With $p_0 \approx p$ an
index swap could hide inside the tolerance; a factor of nine cannot.

The last block is the one that would survive my being wrong about everything else.

## The rest of the four-fermi story

The runnable part above is the *enabling feature*. The full four-quark path adds three more things,
all exercised by `numtracer/tests/gen/gen_fierz_ortho_numeric.wls` (tier C — it needs the
`FierzCompleteNf2Nc3NJL` basis from TensorBases):

**Two independent spinor lines.** A four-quark diagram has two quark lines that are *not* joined
into one trace. The `loop_sep` token marks the boundary: the contraction traces each loop
separately and multiplies the resulting Lorentz tensors, contracting their shared gluon legs through
the Lorentz net. It is emphatically **not** a Wick pairing across the loops — the spinor indices of
one line never meet the other's.

**A two-term inner product with a minus sign.** This basis's inner product is

```text
2 T1[1,2,3,4] T2[2,1,4,3]  -  2 T1[1,2,3,4] T2[2,3,4,1]
```

so the contraction is itself a `Plus` — exactly the shape of the historical multi-term-projector
bug, where an eager Dirac `Plus` of mixed $\gamma$-parity traced to an identically zero kernel. That
is why the gate exists.

**Momentum-independent contact structures.** Four-fermi operators are contact terms, so the trace
must come out **constant**. A kernel that depends on the kinematics is wrong by inspection — an
unusually cheap correctness signal.

### The symmetry trap

[step-14](step-14.md) said the field-list form of `FMakeSymmetryList` generates the full permutation
group of the listed legs, and that you should reach for it only when you have checked the projector
is invariant. Four-fermi is the case where the naive reading fails:

```{admonition} From QCD_Nf2/with_mesons/QCD.wl
:class: warning
> `FierzCompleteNf2Nc34D` is momentum-independent, so covariance is purely algebraic here: ONLY the
> combined exchange `qb1<->qb2` TOGETHER WITH `q1<->q2` is respected (**78/78 diagrams**). A single
> `qb<->qb` (or `q<->q`) swap is a **Fierz rearrangement into OTHER basis elements** and is
> respected by just **14/78** — so the auto field-list form, which would build the full fermionic
> group including those single swaps, must NOT be used here.

```mathematica
sym = FMakeSymmetryList[FSymmetry[Symmetric, {i1, i3}, {i2, i4}]];
```
```

Read that carefully. The four legs *are* two identical quarks and two identical antiquarks, so "the
full permutation group of the same-field legs" looks obviously right. It is wrong, because a single
$\bar q \leftrightarrow \bar q$ swap is not a symmetry of a basis *element* — it is a Fierz
transformation that maps it onto a combination of *other* elements. The reduction would then be
valid only for the 14 diagrams where the difference happens to be annihilated.

This is the concrete case [step-12](step-12.md)'s condition ($P\cdot\sigma(X) = P\cdot X$ for all
$X$) was warning about, and 78 versus 14 is what it costs to get wrong.

### The dressing-rule ordering trap

Also four-fermi, also from the production tree, and also silent:

```{admonition} From QCD_Nf2/with_mesons/QCD.wl
:class: warning
> The field list is FunKit's canonical Grassmann ordering `{qb,qb,q,q}` (what `FOrderFields` leaves
> behind), **NOT** the `{qb,q,qb,q}` spelling of the truncation entry: written the latter way the
> pattern **never fires** and the raw `dressing[...]` head survives into the integrand, where it is
> CForm'd into the kernel as an undefined `dressing(GammaN, List(...))`. The leading MINUS is the
> sign of the single `q<->qb` transposition that relates the two orderings.
```

Two lessons: write four-fermi dressing rules in the canonical ordering, and remember that the
transposition relating the two spellings carries a sign. [step-22](step-22.md) covers recognising
the `dressing(...)` symptom.

## Possibilities for extensions

1. **Break the indexing.** In a copy of the script, replace `ntGamma[0, …]` with `ntGamma[1, …]`
   in assertion B's first net and confirm the checker catches the swap. Then set $p_0 = p$ and watch
   the same error become invisible — the argument for separated scales.

2. **A fixed component on a sigma leg.** Build
   `ntSigma[{"free", 0}, {"slash", {{1, q1}}}, d1, d2]` and confirm via
   `FreeQ[ntk[[1]]["Diagrams"], {"free", _Integer}]` that the fixed component was rewritten to a
   *slash* leg rather than left as a free one. A fixed component is not a free leg.

3. **All four components.** Emit $\mathrm{tr}(\gamma^i \slashed q_1)$ for $i = 0,1,2,3$ and check
   each against $4 q_{1,i}$. In `propFrameFT`, $q_1$ has zero components 2 and 3, so two of these
   must be *exactly* zero — a good test of whether your net degenerates gracefully.

4. **Do it without the feature.** Rewrite one of your own $\gamma^0$-carrying nets entirely in the
   explicit-$e_0$ style. It works, it is more verbose, and it makes concrete what the rewrite is
   doing on your behalf.

5. **Read the Fierz gate.** With this page's small version understood,
   `tests/gen/gen_fierz_ortho_numeric.wls` becomes readable: the weighted-prime trick from
   [step-11](step-11.md), applied to a basis whose structures carry $\gamma^0$ and whose diagrams
   have two spinor lines.

## The plain program

```{literalinclude} ../../../Tutorials/step-19-four-fermi/fixed_components.wls
:language: mathematica
```

```{literalinclude} ../../../Tutorials/step-19-four-fermi/check_fixed_components.cpp
:language: cpp
```

Next, [step-20](step-20.md): what to do when a flow like the ones described here will not generate
at all.
