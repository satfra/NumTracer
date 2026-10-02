# step-18: Dressed numerators and vertex collection

*Builds on: [step-03](step-03.md), [step-05](step-05.md) · Built on by: — ·
Tags: `dressing`, `dpoly` · **Tier A** (a C++20 compiler, nothing else)*

## Introduction

[step-10](step-10.md)'s generated kernel had four mysterious quantities computed before `fill`:

```cpp
const double dr_0 = Mq(l1);
const double dr_1 = -sqrt(powr<-1>(powr<2>(l1))) * RF(powr<2>(k), powr<2>(l1)) * Zq(k) - Zq(l1);
const double dr_2 = Mq(sqrt(powr<2>(l1) - 2.*cos1*l1*p + powr<2>(p)));
const double dr_3 = …;
```

Dressings, being passed **into the trace** rather than multiplying it. This step is why.

### The problem: a numerator that is a sum

A dressed quark propagator numerator is not a single Dirac structure. It is a sum of them, with
*runtime* coefficients:

$$
S(p) \;=\; M_q\,\mathbb{1} \;+\; Z(p)\,\slashed p .
$$

Nothing is wrong with that until you put several in one diagram. A trace containing $D$ such
numerators, if you distribute the sums, is $2^D$ separate traces — one per choice of structure in
each numerator. A four-point flow with four dressed propagators is 16 traces where the tensor
structure only ever had one shape; `with_mesons`'s dense flows would be far worse.

And it is *wasteful in a specific way*: the $2^D$ traces are not independent computations. They are
the same Dirac chain with different tokens spliced in, and the expensive part — the $4\times4$
matrix products and the Lorentz contraction — is nearly identical across all of them.

### The mechanism: keep it eager

NumTracer does not distribute. It keeps each dressed numerator **eager** as a *slot* — the list of
its structure options — contracts the chain **once**, and collects the result into a `DPoly`:

> a polynomial whose **variables** are the dressing calls and whose **coefficients** are the
> kinematic polynomials (`nt::Poly`) the engine already computes.

The dressings ride along as opaque atom-ids and never enter the trace arithmetic. So the Dirac and
Lorentz work is done a single time regardless of how many structures each numerator carries, and the
whole diagram lowers to **one** trace function whose dressing factors the shared CSE collects.

This is the Dirac-side twin of [step-17](step-17.md): there a group-diagonal $\delta$ folds to a
`SUNPoly`, here a dressed numerator folds to a `DPoly`. Same idea, different sector.

```{admonition} It is on by default
:class: note
`NumTrace` takes `"DressingCollection" -> True` by default, so any flow with dressed propagator
numerators already benefits. The `dr_N` variables in [step-10](step-10.md)'s kernel are exactly this
mechanism's footprint.
```

### The data model

| Type | Meaning |
|---|---|
| `DSlotOpt{coeff, dress, toks, netFacs}` | one structure option: a coefficient × a product of dressing atoms × a Dirac structure |
| `DSlot = vector<DSlotOpt>` | one dressed numerator = the sum of its options |
| `DChainTok` | one token of a dressed chain: either a **fixed** factor (`dtfix`) or a **slot reference** (`dtslot i`) |
| `DPoly` | the collected result: dressing monomials → kinematic `Poly` coefficients |

The `toks` / `netFacs` pair is more general than this step needs, and the generality is the point:

* `toks` is a **Dirac-token chain** spliced in place of the slot. Empty = the identity $\mathbb{1}$;
  `{nt::slash(p)}` = a slash; `{nt::gamma(mu)}` = an **open-legged** $\gamma^\mu$.
* `netFacs` is a set of **Lorentz-net factors** carrying any remaining open legs (a vector $p^\mu$,
  a metric $g^{\mu\nu}$, …).

The union of their free Lorentz ids is the slot's open-leg set, closed by the surrounding net — and
it is the *same set for every option*, so the net contracts a fixed leg set regardless of which
structure is chosen. That one form covers every case at any leg count: a propagator numerator
($k=0$), a single-gluon vertex $\gamma^\mu$ or $\slashed p_1\gamma^\mu$ or
$\sigma^{\mu\nu}\slashed p_\nu$ ($k=1$), a two-gluon vertex ($k=2$), and so on. Which is why the
same machinery does **vertex collection**, below.

## The commented program

`Tutorials/step-18-dressed-numerators/dressed_numerators.cpp` takes a quark bubble

$$
\mathrm{tr}\big[\gamma^\mu\,S(p)\,\gamma^\nu\,S(q)\big]\,\delta_{\mu\nu},
\qquad S = M_q\mathbb{1} + Z\,\slashed{\cdot}
$$

with two dressed numerators. Distributed that is $2\times2 = 4$ traces; collected it is one `DPoly`.

### The slots

```{literalinclude} ../../../Tutorials/step-18-dressed-numerators/dressed_numerators.cpp
:language: cpp
:start-after: "@snip begin: slots"
:end-before: "@snip end: slots"
```

The identity structure is the **empty token list** — a slot option that contributes nothing to the
chain but its dressing. That is worth pausing on: `𝟙` is not a token you insert, it is the absence
of one.

Note that atom 0 is **shared** between the two numerators: both mass terms are the same runtime
$M_q$. Atoms 1 and 2 are the per-leg $Z$'s. Sharing an id is how the collection knows two terms
carry the same dressing and can be merged.

### Collecting

```{literalinclude} ../../../Tutorials/step-18-dressed-numerators/dressed_numerators.cpp
:language: cpp
:start-after: "@snip begin: collect"
:end-before: "@snip end: collect"
```

A dressed chain is a list of `DChainTok`: `dtfix` for an ordinary factor, `dtslot i` for a reference
into the slot list. `F.trace(chain, slots, net)` contracts it once and returns the `DPoly`.

The reference computation in the program does the opposite — it enumerates the $2\times2$ choices,
builds each concrete undressed chain, traces it the *ordinary* way with `F.trace`, and weights
by the product of that choice's dressings. That is the $2^D$ path, written out so the two can be
compared.

## Results

```bash
cmake --build build --target dressed_numerators && ./build/dressed_numerators
```

```text
Dressed quark bubble  tr( γ^μ S(p) γ^ν S(q) ) δ_{μν},  S = Mq·𝟙 + Z·slash
  DPoly dressing monomials (one trace collected) : 2   (distributed would be 4 traces)
  collected == distributed over 5000 random points : worst |Δ| = 0.00e+00  (0 bad)
ALL TESTS PASSED
```

**Two monomials, not four.** The four structure choices produced only two distinct dressing
monomials, because $\mathrm{tr}[\gamma^\mu\,\slashed p\,\gamma^\nu\,\mathbb{1}]\delta_{\mu\nu}$ and
its partner vanish — a trace of three gammas is zero. The collection discovered that during the
contraction and never formed the terms. A $2^D$ distribution would have built all four and thrown
two away.

**The agreement is exactly zero**, not $10^{-12}$. Collected and distributed are not two algorithms
that happen to agree numerically; they are the same arithmetic in a different order, and the `DPoly`
coefficients are literally the polynomials the distributed path would have produced.

```{admonition} Sizing the dressing values
:class: warning
`F.eval(dp, pt, drVal)` indexes `drVal` **by dressing id**, so it must be sized to `max-id + 1` —
not to "the number of dressings I used". The program uses ids {0, 1, 2} and therefore a 3-element
array; a 2-element one is refused with a message naming the missing id.
```

## Vertex collection: the same trick, one level up

A dressed *vertex* has the same shape as a dressed propagator numerator — a sum of structures with
runtime coefficients — but with **open legs**. The quark–gluon vertex of [step-14](step-14.md) has

$$
\Gamma^\mu \;=\; \lambda_1\,\gamma^\mu \;+\; \lambda_4\,\slashed p_1\gamma^\mu
              \;+\; \lambda_7\,\sigma^{\mu\nu}p_\nu + \dots
$$

which is exactly `toks = {gamma(mu)}`, `{slash(p1), gamma(mu)}`, `{comm(mu, p)}` — the
general `DSlotOpt` form. So the same collection applies, and the front-end calls it `ntDiracSlot`.

It is **opt-in**, via `NT_VERTEX_COLLECT=1` or `$ntVertexCollect = True`, because the trade is not
always favourable:

* Measured **36.8× diagram reduction** on the quark–gluon vertex flow (`probe_za3_147_count.wls`).
* But on `za3_147` the collected generation is 8.9 s versus 13.3 s distributed — a win — while on
  flows with few structures per slot the collection overhead is not repaid.

A round-trip assertion (`tests/gen/test_diracslot_roundtrip.wls`) checks that collected and
distributed agree, which is the same guarantee this program checks at the C++ level.

```{admonition} The σ^{μν} fold
:class: note
FunKit has no $\sigma$ primitive: the struct-7 quark–gluon tensor arrives as a bare two-term
antisymmetric $\gamma$-pair, $s(\gamma^X\gamma^Y - \gamma^Y\gamma^X)$. Left alone that `Plus`
distributes into **two** traces. `foldDiracSigma` recognises it and keeps it as ONE `ntSigma` token,
so the antisymmetric pair is never split. Disable with `NT_NO_SIGMA_FOLD` if you need to compare.
```

## Possibilities for extensions

1. **Add a third structure.** Give `S(p)` a third option — say $Z_2(p)\,p^2\,\mathbb{1}$ with its
   own dressing id. Distributed that is $3\times2 = 6$ traces; check how many `DPoly` monomials
   result and whether the reference loop still agrees.

2. **Why only two monomials?** Count the gammas of each structure choice: the two mixed
   mass–slash choices carry an odd number and vanish identically, which is physics, not a dropped
   term. Insert one fixed `nt::slash(p)` into the chain and predict which two monomials survive
   now; then check.

3. **Share an id across three terms.** Give all options the same dressing id and confirm the
   `DPoly` collapses to a single monomial whose coefficient is the sum of the individual traces.

4. **An open-legged slot.** Build a slot whose options are `{gamma(mu)}` and
   `{slash(p), gamma(mu)}` — a miniature dressed vertex — and close `mu` with a metric. This is
   vertex collection by hand, and it demonstrates that the `toks`/`netFacs` generality is not
   speculative.

5. **Break it deliberately — mis-size `drVal`.** Pass a 2-element array with ids {0,1,2} in play.
   `F.eval` refuses: "the polynomial carries dressing 2 but only 2 dressing values were given".

6. **Measure the collection.** Time the collected path against the distributed reference as $D$
   grows (2, 3, 4 dressed numerators). The distributed cost doubles each time; the collected cost
   should barely move.

## The plain program

```{literalinclude} ../../../Tutorials/step-18-dressed-numerators/dressed_numerators.cpp
:language: cpp
```

Next, [step-19](step-19.md): four-fermi operators — two independent spinor lines, contact
structures with no momentum dependence, and $\gamma^0$ as a slash against a unit vector.
