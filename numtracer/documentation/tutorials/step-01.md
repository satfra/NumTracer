# step-01: Hello, tensor network

*Builds on: nothing · Built on by: [step-02](step-02.md), [step-03](step-03.md),
[step-04](step-04.md) · Tags: `network`, `projector` · **Tier A** (a C++20 compiler, nothing else)*

## Introduction

NumTracer exists to answer one question over and over: *given a product of tensors with shared
indices, what single number does it collapse to?* Everything else in the library — the code
generator, the Dirac algebra, the SU($N$) tables, the CSE pass — is scaffolding around that one
operation. So the first program does exactly it, with no physics in sight.

### The one rule

A **network** is a product of tensor factors, each carrying *index labels*. Two index slots are
summed over exactly when they carry **the same label** — the Einstein convention, made literal. A
label that appears once stays free; a network with no free labels left is a scalar.

That rule is the entire mental model, and it is worth being precise about what it does *not* say.
It says nothing about "upper" and "lower" indices — NumTracer is Euclidean, so raising and lowering
is the identity and there is no distinction to track. And labels of different kinds never meet: a
Lorentz label, an adjoint SU($N$) label and a fundamental one are different C++ types, so the
compiler keeps the sectors apart (see [Scope & conventions](../getting_started/scope-and-conventions.md)).

### The two networks

We contract two networks and check each against its closed form. First, the plainest possible
contraction — two vectors tied through a metric:

$$
a\cdot b \;=\; a_\mu\,\delta^{\mu\nu}\,b_\nu .
$$

Written as a network that is three factors — $a$ on label $\mu$, the metric on $\mu$ and $\nu$,
$b$ on $\nu$ — and both labels are used twice, so both are summed away.

Second, something with structure: a **transverse projector**,

$$
P(k)_{\mu\nu} \;=\; \delta_{\mu\nu} - \frac{k_\mu k_\nu}{k^2},
\qquad\text{so}\qquad
a\cdot P(k)\cdot a \;=\; a^2 - \frac{(a\cdot k)^2}{k^2}.
$$

$P(k)$ is the object that projects out the component of a vector along $k$; in a gauge theory it is
what sits on an internal gluon line. It matters here for a structural reason: it carries a
**denominator**, $1/k^2$. Contraction cannot evaluate that — $k^2$ is generally a runtime quantity —
so the engine must carry it through symbolically. The way it does so (an *inverse atom*) is the
subject of [step-04](step-04.md); this program just meets the idea.

```{admonition} Why not simply build the tensors and sum?
:class: note
You could materialise $P(k)$ as sixteen numbers and loop. That works here and does not work later:
a Dirac trace of $n$ gammas has $(2n-1)!!$ index pairings, and a Lorentz trace of $n_p$ transverse
projectors has $2^{n_p}$ terms if you expand the definition. The engine's contraction never forms
either intermediate. Step 1 uses the small case to introduce the vocabulary; steps 2–4 show the
same vocabulary doing the thing you could not do by hand.
```

## The commented program

The whole program is `Tutorials/step-01-hello-network/hello_network.cpp`. We take it in three
pieces.

### The frame, the momenta and the labels

```{literalinclude} ../../../Tutorials/step-01-hello-network/hello_network.cpp
:language: cpp
:start-after: "@snip begin: frame"
:end-before: "@snip end: frame"
```

A `Frame` holds everything the contraction needs to know besides the network itself: which momenta
exist and what their four components are. Here the components are plain numbers. In a real kernel
they are *symbolic* — functions of the loop magnitude and angles that the integrator supplies at
runtime — and the frame declares those symbols too ([step-03](step-03.md) onwards).
`F.momentum(...)` returns a `Momentum`, the handle you use in the network.

`F.indices<2>()` hands out two fresh Lorentz index labels. Each call returns labels that differ
from every label handed out before, so two factors contract only where you deliberately reuse one.
A label is its own type, `LorentzIndex`: an integer or a momentum cannot be passed where an index is
expected.

### Network one: a metric

```{literalinclude} ../../../Tutorials/step-01-hello-network/hello_network.cpp
:language: cpp
:start-after: "@snip begin: dot"
:end-before: "@snip end: dot"
```

`nt::vec(mu, va)` is the vector $a_\mu$, `nt::metric(mu, nu)` is $\delta_{\mu\nu}$, and `*` is
the tensor product. The shared labels are summed when `F.contract` collapses the network. The
result is an `nt::Poly`, a polynomial in the frame's symbols — here a constant, since there are no
symbols. `F.eval(poly, point)` evaluates it; `F.at(...)` builds the point from the symbol values,
none in this case.

A network can also be a sum: `+` adds two networks, and a number times a network scales it.

### Network two: a projector

```{literalinclude} ../../../Tutorials/step-01-hello-network/hello_network.cpp
:language: cpp
:start-after: "@snip begin: proj"
:end-before: "@snip end: proj"
```

`nt::projT(mu, nu, vk)` is $P(k)_{\mu\nu}$. Its $1/k^2$ cannot be multiplied out — $k^2$ is in
general a runtime quantity — so it rides through the contraction as a separate factor, an **inverse
atom**. The frame registers that atom the first time it meets the projector: it knows $k$'s
components, so it knows the denominator $k^2$, and `F.eval` computes $1/k^2$ at the evaluation
point. Knowing the denominator also lets the engine cancel a $k^2$ against a $1/k^2$ inside a
monomial instead of carrying both. [step-04](step-04.md) looks at atoms more closely.

## Results

```bash
cmake --build build --target hello_network && ./build/hello_network
```

```{literalinclude} ../../../Tutorials/step-01-hello-network/hello_network.expected.txt
:language: text
```

Both agree with the closed form to machine precision, which is the point: the engine was told
nothing about what a metric or a projector *means*. It was told that $\delta_{\mu\nu}$ has entries
$1$ on the diagonal and that $P$ has entries $\delta_{\mu\nu} - k_\mu k_\nu \cdot (1/k^2)$,
and it summed the shared labels. The identity $a\cdot P(k)\cdot a = a^2 - (a\cdot k)^2/k^2$ is not
knowledge the engine has — it is what falls out.

Notice also what the second result *is*. With $k$ fixed and $a$ varying, $a\cdot P(k)\cdot a / a^2$
is $1 - \cos^2\theta$, the angular weight that appears in every one-loop self-energy. You have
already computed a piece of real physics; [step-04](step-04.md) does it with the angle symbolic.

## Possibilities for extensions

1. **Check that the projector projects.** $P(k)$ should annihilate $k$: contract
   `vec(mu, vk) * projT(mu, nu, vk) * vec(nu, vk)` and confirm you get 0 to machine precision. Then
   check idempotence, $P\cdot P = P$, by tying two projectors through a third label
   (`auto [rho] = F.indices<1>();`) and comparing to a single one.

2. **A momentum that is a combination.** Replace the second `va` in network two by `va - vb`.
   Predict the closed form first, then check it. Momenta add and scale like vectors; this is how
   every propagator momentum $l - p$ is written.

3. **Reuse a label within one factor.** Change the metric to `metric(mu, mu)` and drop the two
   vectors. You get $\mathrm{tr}\,\delta = 4$: a label used twice within one factor is a
   self-contraction, perfectly legal. The engine cannot tell an intended contraction from a typo —
   labels are the whole contract, which is why each label should come from `F.indices`.

4. **Break it deliberately — leave an index open.** Contract `vec(mu, va) * metric(mu, nu)` (drop
   the second vector). Summing a free index over 0–3 would return a meaningless number, so the
   engine refuses:

   ```text
   numtracer: Lorentz index id 1 is OPEN (occurs once) — the network does not close to a scalar. ...
   ```

   Now pass a wrong number of values to `F.at` (say `F.at(1.0)`): the frame has no symbols, and
   says so.

5. **Add a longitudinal projector.** `nt::projL(mu, nu, vk)` is $k_\mu k_\nu / k^2$. Verify
   $P^T + P^L = \delta$ by contracting `projT(...) + projL(...)` against two vectors and comparing
   with the metric.

## The plain program

```{literalinclude} ../../../Tutorials/step-01-hello-network/hello_network.cpp
:language: cpp
```

Next, [step-02](step-02.md): the same idea in the SU($N$) sector, where the network folds to a
number without the tensor ever existing.
