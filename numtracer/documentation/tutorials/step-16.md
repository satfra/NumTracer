# step-16: Finite temperature

*Builds on: [step-04](step-04.md), [step-07](step-07.md) · Built on by: — ·
Tags: `finite-T`, `projector` · **Tier A** (a C++20 compiler, nothing else)*

## Introduction

At $T > 0$ the Euclidean $O(4)$ rotation symmetry breaks to spatial $O(3)$. The heat bath singles
out a rest frame $u^\mu = (1,0,0,0)$, and the **zeroth** — Euclidean-time, Matsubara — component of
every momentum becomes physically distinct from the spatial three.

The remarkable thing, from the engine's point of view, is how little this changes.

* **The metric is still the Euclidean $\delta_{\mu\nu}$.** In the Matsubara formulation nothing
  about the contraction changes. There is no new signature and no new index type.
* **The Dirac algebra is untouched.** The chain in this program is character-for-character the
  vacuum chain of [step-05](step-05.md).
* **NumTracer already contracts on explicit frame components**, so a momentum with an independent
  temporal component costs exactly nothing extra. You just choose a frame whose temporal slot is an
  independent symbol.

What *does* change is the **projectors**, and that is genuinely new physics.

```{admonition} Convention
:class: note
Component index `0` is the Euclidean-time / Matsubara direction (heat bath $u = (1,0,0,0)$);
`1, 2, 3` are spatial. This is the contract between a finite-$T$ frame and the electric/magnetic
projectors, and it is not negotiable — the projectors have it baked in.
```

### The thermal split

The 4d transverse projector is unchanged and remains valid at finite $T$:

$$
P^T_{\mu\nu}(k) = \delta_{\mu\nu} - \frac{k_\mu k_\nu}{k^2}.
$$

But it is no longer *irreducible*: with a preferred rest frame available, it splits into a
**magnetic** (spatial-transverse) part and an **electric** (time-like-transverse) part,

$$
P^M_{ij}(k) = \delta_{ij} - \frac{k_i k_j}{|\vec k|^2}\ (i,j \in \{1,2,3\}),\qquad
P^M_{0\nu} = P^M_{\mu 0} = 0,
$$
$$
P^E = P^T - P^M, \qquad\text{so}\qquad \operatorname{tr}P^M = 2,\quad \operatorname{tr}P^E = 1 .
$$

Physically these are the two gluon polarisations that no longer have to share a dressing: the
chromo-magnetic and chromo-electric sectors screen differently, and a finite-$T$ Yang–Mills flow
solves for $Z_A^E$ and $Z_A^M$ separately.

### The trap: two different denominators

$P^M$ divides by the **spatial** $|\vec k|^2 = k^2 - k_0^2$, not by $k^2$. $P^E$ needs *both*. So a
finite-$T$ network carries **two inverse atoms** where a vacuum one carried one, and confusing them
is the characteristic finite-$T$ bug — it produces a kernel that is correct in the $k_0 \to 0$ limit
and wrong everywhere else.

### The new heads

| Head | Meaning |
|---|---|
| `ntElectricProj[q, mu, nu]` | $P^E$ — needs both the full and the spatial atom |
| `ntMagneticProj[q, mu, nu]` | $P^M$ — needs only the spatial atom |
| `ntSPS[a, b]` | the **spatial** scalar product $\vec a\cdot\vec b = a\cdot b - a_0 b_0$ |
| `ntSpatialVec[q]` | the **spatial vector** $\bar q = \{0, q_1, q_2, q_3\}$ — a momentum, not a scalar |
| `ntVec[q, 0]` | an **integer** second argument: the *scalar* component $q_0$ (e.g. $\pi T$) |
| `propFrameFT[p0, p, l0, l1, cos1, …]` | a frame with independent temporal slots |

In C++ these are `nt::projE(mu, nu, l)` and `nt::projM(mu, nu, l)`; the frame registers both
denominators ($l^2$ and $|\vec l|^2$) itself.

### The spatial *vector* — and the spatial pslash

FormTracer's finite-$T$ vocabulary has two spatial objects, and they are not interchangeable:
`sps[a, b]` (a scalar) and `vecs[q, mu]` (a vector). `FromFunKit` maps them to `ntSPS` and to
`ntVec[ntSpatialVec[q], mu]` respectively.

`ntSpatialVec[q]` is **engine-introduced** — you do not write it by hand any more than you write
`ntUnitVec`. It is an ordinary momentum whose frame components are the parent's with slot 0 zeroed,
and that single fact is the whole implementation: `NumTrace` pushes it through sums (it is linear,
so `vecs[p-l, mu]` becomes $\bar p - \bar l$ rather than a third momentum) and injects the components
into the frame next to the unit basis vectors.

The payoff is the **spatial pslash**, which is how a finite-$T$ quark line is written:

```wl
(* vecs[p, mu] gamma[mu, d1, d2]  —  the FormTracer form *)
ntGamma[mu, d1, d2] ntVec[ntSpatialVec[pp], mu]
```

This needs no new engine token: a $\gamma$ whose Lorentz index carries an `ntVec` **is** a slash, so
it is emitted as an ordinary `slash` against a different momentum. FormTracer's one-argument
shorthand `gamma[..., vecs[p], ...]` is expanded to the two-argument form by `FromFunKit` before
anything else looks at it.

```{admonition} `sps` is a coefficient, `vecs` is a leg
:class: important
`ntSPS[a, b]` never reaches the engine — `scalarQ` classifies it as a scalar and the frame folds it
into the diagram coefficient. `ntSpatialVec` does reach the engine, as a momentum with its own env
`Base`. The codegen gate `gen_spatialvec_numeric.wls` pins them against each other: one of its nets
is the identity $\mathrm{tr}(\bar{p}\!\!\!/\,\bar{l}\!\!\!/) = 4\,\texttt{sps}(p,l)$, i.e. the two
implementations of "spatial" required to agree.
```

```{admonition} `ntVec[q, mu]` versus `ntVec[q, 0]`
:class: important
The second argument decides which object you get. A **symbolic** label (`ntVec[q, mu]`) is the usual
tensor leg, to be contracted. A **literal integer** (`ntVec[q, 0]`) is the *scalar component* $q_0$,
resolved by the frame like `ntSP` — a coefficient, not a factor. This is the same
tensor-versus-scalar distinction as [step-08](step-08.md), controlled by the argument's type.
```

## The commented program

`Tutorials/step-16-finite-temperature/finite_temperature.cpp`.

### Two atoms

```{literalinclude} ../../../Tutorials/step-16-finite-temperature/finite_temperature.cpp
:language: cpp
:start-after: "@snip begin: atoms"
:end-before: "@snip end: atoms"
```

There is nothing to write here for the denominators: a projector knows its momentum, the frame
knows the momentum's components, and the projector's *kind* says which denominator it needs —
$l^2 = \sum_{\mu=0}^3 l_\mu^2$ for $P^T$ and $P^E$, $|\vec l|^2 = \sum_{\mu=1}^3 l_\mu^2$ for $P^E$ and
$P^M$. Component 0 is the heat-bath direction. After the first contraction,
`F.denominators()` holds both.

Note also that *all eight* components are independent symbols. In a vacuum frame you would exploit
$O(4)$ to zero several of them ([step-04](step-04.md)); here the temporal ones are physical and must
survive.

### Checking the split through the engine

```{literalinclude} ../../../Tutorials/step-16-finite-temperature/finite_temperature.cpp
:language: cpp
:start-after: "@snip begin: traces"
:end-before: "@snip end: traces"
```

This is the check worth stealing. Contracting a projector against a metric on the *same two legs*
is its trace, $\delta_{\mu\nu}P_{\mu\nu} = \operatorname{tr}P$. Those traces must come out as the
pure numbers 1, 2 and 3 — **with no kinematic dependence at all**, at every point, for any $l_0$.

It is a strong check precisely because it is so rigid. Any confusion between the two atoms, any
leakage of a temporal row into $P^M$, any sign error, breaks the constancy immediately — and it
tests the *engine's* contraction rather than agreement with a hand-written oracle that could share
the same mistake.

### The self-energy

```{literalinclude} ../../../Tutorials/step-16-finite-temperature/finite_temperature.cpp
:language: cpp
:start-after: "@snip begin: chain"
:end-before: "@snip end: chain"
```

The Dirac chain is identical to [step-05](step-05.md)'s. The gluon line is now a **two-term**
Lorentz network — one term per projector, each with its own dressing coefficient — written exactly
as $G_{\mu\nu}(l) = Z_A^E P^E + Z_A^M P^M$: `ZAE * projE(…) + ZAM * projM(…)`.

## Results

```bash
cmake --build build --target finite_temperature && ./build/finite_temperature
```

```text
A. the thermal projector split, contracted through the engine
   tr P_E = 1, tr P_M = 2, tr P_T = 3 : worst |error| over 5000 points = 0.000e+00
B. finite-T quark self-energy (p_0 = πT, independent loop l_0)
   engine vs closed-form trace identity : worst relative = 1.339e-12
ALL TESTS PASSED
```

(The error estimates are rounding-level; their last digits depend on the compiler.)

The traces are constant to machine epsilon across 5000 random kinematic points with random $T$,
random $\vec p$ and a random independent $l_0$. That is the projector algebra confirmed, not
assumed.

The self-energy agrees with

$$
N = 4\big[\,2\,(p\cdot G\cdot q) - (p\cdot q)\operatorname{tr}G\,\big],
\qquad \operatorname{tr}G = Z_A^E + 2Z_A^M,
$$

to about $10^{-12}$ relative — the loss relative to the trace check is ordinary cancellation in the
oracle's explicit $4\times4$ sums, not engine error.

```{admonition} Why the external momentum sits at πT
:class: note
$p_0 = \pi T$ is the lowest fermionic Matsubara frequency, and it is chosen here for a *testing*
reason as much as a physical one: at a $T = 0$ vacuum frame every $q_0$ vanishes, the electric and
magnetic projectors degenerate, and **both checks would pass vacuously** no matter what the engine
did. A finite-$T$ gate evaluated at zero temporal components tests nothing.
```

## From the Mathematica front-end

The same objects exist in the DSL, so a generated finite-$T$ kernel needs only a finite-$T$ frame
and the new heads:

```mathematica
frame = propFrameFT[p0, p, l0, l1, cos1, pp, ll];   (* component 0 = Matsubara/temporal *)

net = ntVec[pp, i1] ntElectricProj[ll, i1, i2] ntVec[pp, i2]
    + ntVec[pp, j1] ntMagneticProj[ll, j1, j2] ntVec[pp, j2]
    + ntSPS[pp, ll] ntVec[ll, 0] *
        ntVec[pp, k1] ntTransProj[ll, k1, k2] ntVec[pp, k2];

ntk = NumTrace[net, "Frame" -> frame, "Args" -> {p0, p, l0, l1, cos1}];
```

(Remember the parentheses if you write that sum across lines — [step-09](step-09.md).)

`propFrameFT[p0, p, l0, l1, cos1, pp, ll]` sets

```text
pp = {p0, p, 0, 0}
ll = {l0, l1 cos1, l1 sqrt(1-cos1^2), 0}
```

with `p0` and `l0` independent of the magnitudes.

```{admonition} A finite-T frame leaves the unit-loop fast path
:class: note
Because `ll`'s temporal component is **not** proportional to the loop magnitude `l1`, the compact
unit-loop optimisation of [step-07](step-07.md) does not apply, and the backend automatically takes
the general frame path (`polyFrameSpec`), minting extra sine/cosine/radical symbols.

Nothing is wrong when this happens — but the kernel is bigger, and if you are wondering why a
finite-$T$ flow generates more slowly than its vacuum counterpart of the same order, this is why.
The projectors also each allocate their spatial $1/|\vec q|^2$ env slot automatically.
```

The generated finite-$T$ kernel is gated by the `ftproj_num` test
(`numtracer/tests/gen/gen_ftproj_numeric.wls`).

## Possibilities for extensions

1. **Verify the split as a network identity.** Build $P^E + P^M - P^T$ as a three-term network and
   contract it against two arbitrary vectors. The result must be identically zero — check
   `poly.size() == 0`, not just that it evaluates small, exactly as in
   [step-04](step-04.md) extension 1.

2. **See why the two denominators matter.** Replace $P^M$ by $P^T - P^E$ written with the *full*
   $l^2$ throughout (build the tensor yourself from `vec` and `metric` with a $1/l^2$ you supply as a
   symbol). Watch `tr P_M` stop being 2 away from $l_0 = 0$ and come back at $l_0 \to 0$ — the
   characteristic finite-$T$ bug is invisible in the vacuum limit.

3. **Check magnetic transversality.** $P^M$ must annihilate the *spatial* part of $l$ and must have
   vanishing temporal rows. Contract it with `vec` on $l$ and confirm; then contract it with a
   purely temporal vector and confirm you get zero.

4. **Dress them differently.** Set `ZAE` and `ZAM` to very different values and watch the
   self-energy change. In a real finite-$T$ flow these are separate interpolators solved for
   separately — the whole reason the split exists.

5. **Generate it.** Port the network to the DSL with `propFrameFT` and generate a kernel
   ([step-06](step-06.md)), then compare `nenv` against the vacuum version of the same network. The
   difference is what leaving the unit-loop path costs.

## The plain program

```{literalinclude} ../../../Tutorials/step-16-finite-temperature/finite_temperature.cpp
:language: cpp
```

```{admonition} Physics dressings
:class: tip
A full finite-$T$ flow additionally needs finite-$T$ dressings and regulators — functions of $l_0$
and $|\vec l|$ separately — supplied in the per-theory setup. Those are a modelling choice and do
not affect the tensor algebra shown here.
```

Next, [step-17](step-17.md): dressing individual colour or flavour components differently, inside a
single trace.
