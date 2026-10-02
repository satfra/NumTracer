# step-04: Lorentz networks and inverse atoms

*Builds on: [step-01](step-01.md) · Built on by: [step-05](step-05.md), [step-16](step-16.md) ·
Tags: `lorentz`, `projector` · **Tier A** (a C++20 compiler, nothing else)*

## Introduction

[step-01](step-01.md) contracted a projector with numeric components. This step does it with the
components **symbolic**, which changes the character of the result completely: instead of a number
you get a polynomial, and that polynomial is what a generated kernel evaluates. Two ideas appear
for the first time — the *frame*, and the *inverse atom*.

### The frame

A loop integrand depends on the external momenta and the loop momentum, but not on all of their
components independently: rotational invariance means only magnitudes and relative angles matter.
A **frame** is the choice of concrete components that encodes this. For a one-loop propagator
diagram — one external momentum $p$, one loop momentum $l$ — the natural choice is

$$
p = (p, 0, 0, 0), \qquad l = (l\cos\theta,\ l\sin\theta,\ 0, 0),
$$

i.e. put $p$ along an axis and let $l$ live in the plane it spans with that axis. Two components of
$l$ are identically zero and can be dropped; the integrand is a function of exactly three scalars,
$p$, $l$ and $\cos\theta$. ($\sin\theta$ is not a fourth one: it is $\sqrt{1-\cos^2\theta}$.)

This is worth stating carefully because it is the single biggest lever on kernel size. The engine
contracts *over the frame's components*. Choosing a frame with two zero components means those
entire index-sum branches vanish before any polynomial is built — the saving is not a
simplification afterwards, it is work never done. [step-07](step-07.md) is entirely about frames
and the five builders the front-end provides; here we hand-write one.

### The inverse atom

A transverse projector

$$
P(l)_{\mu\nu} = \delta_{\mu\nu} - \frac{l_\mu l_\nu}{l^2}
$$

carries a denominator. In a symbolic frame $l^2$ is a polynomial in the frame's symbols — and
dividing one polynomial by another does not in general give a polynomial. So the engine cannot
simply "do the division"; it must carry $1/l^2$ as an opaque quantity.

It does this by giving each distinct denominator an **atom id**, and letting monomials carry
powers of atoms alongside powers of symbols. A monomial is then

$$
c \cdot \prod_i x_i^{e_i} \cdot \prod_j \mathrm{atom}_j^{a_j},
$$

and multiplication just adds exponents. Crucially the engine is *also* told what each atom is the
reciprocal of, which lets it cancel: when a monomial acquires both an $l^2$ from a numerator and an
$\mathrm{atom}_0$, the two annihilate instead of both being carried. Without that cancellation the
polynomials would grow without bound through a long chain of propagators.

```{admonition} Who numbers the atoms
:class: note
The frame does. Each projector's momentum is known, so the first time the frame meets a projector on
$l$ (or on $-l$, which has the same $l^2$) it registers one atom with denominator $l^2$, computed
from $l$'s components; every later projector on the same momentum reuses it. That sharing is what
lets the denominators cancel against the same numerators. Generated kernels do the same bookkeeping
in the front-end, one atom per distinct inverse propagator across the whole diagram.
```

## The commented program

`Tutorials/step-04-lorentz-networks/lorentz_networks.cpp` computes $p\cdot P(l)\cdot p$ in the
one-angle frame.

### Writing the frame down

```{literalinclude} ../../../Tutorials/step-04-lorentz-networks/lorentz_networks.cpp
:language: cpp
:start-after: "@snip begin: frame"
:end-before: "@snip end: frame"
```

Three symbols: the magnitudes `P`, `L` and the angle cosine `C`. `F.angle("cos")` returns the
cosine *and* the sine as symbols, but only the cosine is an input: the frame derives
$\sin\theta = \sqrt{1-\cos^2\theta}$ at evaluation time, and it uses $\cos^2+\sin^2 = 1$ during
contraction. That is what makes $l^2 = l^2\cos^2\theta + l^2\sin^2\theta$ collapse to the single
monomial $l^2$ — and a denominator that is a single monomial can be cancelled exactly.

Components given as `0` are *structural* zeros: they produce polynomials with no monomials, so any
product they enter is dropped immediately rather than carried as a term with coefficient zero.

Compare this with [step-03](step-03.md), which used eight independent symbols for two momenta. Both
are legitimate; they answer different questions. Eight symbols verify an identity *in general*;
three symbols compute the thing a kernel actually needs. Real generation always uses the frame.

### The network, and the contraction

```{literalinclude} ../../../Tutorials/step-04-lorentz-networks/lorentz_networks.cpp
:language: cpp
:start-after: "@snip begin: net"
:end-before: "@snip end: net"
```

`nt::projT(mu, nu, l)` is the projector on the loop momentum $l$, legs `mu` and `nu`. Both labels
appear twice across the three factors, so both are summed and the network closes. `F.contract` is
the pure-Lorentz form of `F.trace` ([step-03](step-03.md)); a whole diagram, Dirac and Lorentz
together, is [step-05](step-05.md).

## Results

```bash
cmake --build build --target lorentz_networks && ./build/lorentz_networks
```

```{literalinclude} ../../../Tutorials/step-04-lorentz-networks/lorentz_networks.expected.txt
:language: text
```

**Two monomials.** That is the headline. The contraction ran over $\mu,\nu \in \{0,1,2,3\}$ —
sixteen index combinations — and what came back is

$$
p^2 \;-\; p^2\cos^2\theta ,
$$

two terms and no atom left. The frame did the work twice. Because $p$ has only component 0, every
term involving the other components of $p$ was structurally absent. And the contraction produced
$p^2 l^2\cos^2\theta \cdot (1/l^2)$, in which the $l^2$ cancelled against the atom because the frame
knows that atom's denominator is exactly $l^2$. This compactness is what makes generated kernels
small.

**The physics.** $p\cdot P(l)\cdot p / p^2 = 1 - \cos^2\theta$, where $\theta$ is the angle between
$p$ and $l$. This is the angular weight of a transverse gluon exchange — the factor that appears in
the ghost and gluon loops of every Yang–Mills propagator flow. In a real kernel it is multiplied by
dressings and a regulator and handed to a quadrature over $\cos\theta$.

**Where the division happens.** An atom that survives contraction is evaluated by `F.eval` as
$1/l^2$ at the point — one division per atom. A generated kernel does the same: it computes each
`1/l²` once per call and then never divides again, which is why the emitted arithmetic is
division-free and GPU-friendly. Change the network to $p\cdot P(l)\cdot q$ with a second external
momentum and you will see an atom survive.

## Possibilities for extensions

1. **Watch the cancellation happen.** Contract $l\cdot P(l)\cdot l$, which is analytically zero.
   Check that `poly.size() == 0` — the terms cancel *during* contraction, not numerically at
   evaluation.

2. **Two projectors.** Add a second `projT` on the same momentum sharing a middle label, i.e.
   $p\cdot P(l)P(l)\cdot p$. Idempotence says the answer is unchanged. Check both the value *and*
   the monomial count, and `F.denominators().size()`: both projectors share one atom.

3. **The longitudinal complement.** Compute $p\cdot P^L(l)\cdot p$ with `nt::projL` and confirm it
   equals $p^2\cos^2\theta$, and that the two add to $p^2$.

4. **Change the frame.** Put $p$ along axis 1 instead of axis 0 and re-run. The value must not
   change (it is a scalar), but look at the polynomial — the *monomials* may differ. Then try
   declaring $l$'s two components as independent symbols instead of $l\cos\theta$, $l\sin\theta$:
   the atom no longer cancels, and the polynomial grows. This is the frame lever from the
   introduction, made visible.

5. **Break it deliberately — a free index.** Drop the second `vec`, leaving `nu` open. The network
   no longer closes, and `F.contract` throws:

   ```text
   numtracer: Lorentz index id 1 is OPEN (occurs once) — the network does not close to a scalar.
   ```

   Read it and remember what it looks like; an unclosed index in a hand-built net is one of the two
   mistakes ([step-01](step-01.md) extension 3 is the other) that account for most first-day
   confusion. The guard matters more than it looks: the contraction core sums *every* index it is
   given over $0..3$, so before it existed this net returned $\sum_\nu (p\cdot P)_\nu$ — a
   finite, plausible, meaningless number, with no diagnostic at all.

6. **A frame with no zero components.** Use four symbols for $l$ and four for $p$ and redo the
   contraction, then compare the monomial count to the three-symbol version. The ratio is the price
   of not choosing a frame — and for a four-point vertex it is the difference between a kernel that
   generates and one that does not ([step-20](step-20.md)).

## The plain program

```{literalinclude} ../../../Tutorials/step-04-lorentz-networks/lorentz_networks.cpp
:language: cpp
```

Next, [step-05](step-05.md): Dirac and Lorentz together — a whole diagram — and then the lowering
pass that turns the resulting polynomial into C++.
