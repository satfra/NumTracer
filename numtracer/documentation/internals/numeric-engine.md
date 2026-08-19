# The numeric contraction engine

> Headers: `numeric/{numeric_contract,mpoly}.hpp`,
> `network/{network,dirac,sun_net}.hpp` · Namespaces: `numtracer::numeric`, `numtracer::network`

This is the heart of how NumTracer generates a kernel. Given one diagram — a closed Dirac chain, a
pure-Lorentz network of projectors and metrics, and a colour network — it contracts everything
**numerically over a fixed loop frame** and returns one small polynomial. It is run at build
time by the [codegen](codegen.md) generator, not in the consumer's compiler.

The contraction avoids the two blow-ups that make this hard:

* the **`(2n−1)!!` Wick blow-up** of a Dirac trace (the number of index pairings), and
* the **`2^np` projector-mask blow-up** of a Lorentz trace of many transverse projectors.

Neither intermediate is ever formed. The result is already a small polynomial.

## The polynomial it works in

`numeric/mpoly.hpp` defines `MPoly`, the arithmetic currency of the whole engine: a multivariate
polynomial over the frame's scalar symbols that also tracks surviving inverse propagator
denominators `1/k²`. Everything below produces and combines `MPoly`s through addition,
multiplication, monomial-level cancellation of `k²·(1/k²)` (`divThroughMonomialAtoms`), and
collection of equal-momentum propagator groups (`reduce_units`).

Every `MPoly` (and `DPoly`) is minted through a `numtracer::numeric::LorentzEnv` (`numeric/env.hpp`),
a small factory that binds the symbol-space size `nsym` **once**: you write `nm::LorentzEnv env(nsym)`
and then call `env.mono(...)`, `env.atom(...)`, `env.var(i)`, `env.zero()`, `env.numeric_value(...)`.
Because the env carries `nsym`, the raw size-taking factories never appear in caller code, and two
polynomials built from the same env are guaranteed to share a symbol space. The illustrative
constructors below are shown as env methods for that reason.

### Data model

An `MPoly` is a sorted list of `(monomial, complex coefficient)` pairs (`MPoly::t`) — nothing
more exotic. A **monomial** (`Mono`) has two independent parts that live in two different
namespaces:

* **An exponent vector** (`Mono::e`), length `nsym`. After γ-matrices, metrics and projectors
  have been contracted numerically, the only symbolic quantities left are a handful of scalar
  *user symbols* — for a typical flow `nsym` is single digits: `l1` (loop momentum), `cos1`,
  angular components, `p` (external momentum). The exponent vector is the list of powers of each
  symbol in one product term: with symbol order `[l1, cos1, …, p]`, the plain product
  `l1²·cos1¹·p³` is stored as `e = [2,1,0,0,3]`. Position = which symbol, value = its power. (This
  is *not* an exponent that is itself an expression — just a tuple of integer powers. The
  generated component table writes these directly, e.g. `env.mono({1,1,0,0,0},…)` on a
  `LorentzEnv env(5)` (five symbols) is `l1·cos1`.)

  It *behaves* as that vector — `e[k]` reads and writes a power — but it is **stored** as a packed
  128-bit key: 5 bits per symbol, 12 symbols per 64-bit word, up to 24 symbols. That takes the
  stored term `pair<Mono,Cx>` from 128 B to 72 B and turns the monomial comparison that dominates
  `MPoly`'s sort into two integer compares instead of an `nsym`-long walk. The packing is
  **big-endian within each word** precisely so that comparing the two words *is* the old
  element-wise lexicographic order — symbol 0 dominates, then symbol 1, … — which is what keeps the
  emitted kernels byte-identical across the change. A symbol index past 24, a power above 31, or a
  negative power transparently spills to a heap exponent list, so `nsym` and degree are bounded only
  by `MonoExpT` (`int16_t`); no committed flow comes close (the largest `nsym` is 6), so in practice
  every monomial stays inline.

* **A sorted multiset of atom ids** (`Mono::atoms`). An **atom id** is an integer that *names a
  projector's denominator `k²`* — it is neither a symbol index nor the value `1/k²`. Each
  transverse/longitudinal/electric/magnetic projector in the network carries an id (`Elem::inv`,
  and `Elem::invS` for the spatial `1/|k⃗|²`); `collect_atom_denoms` interns the actual
  denominator polynomial once into a separate table, `atomDen[aid] = k² = Σ_μ comp[μ]²` (itself
  an `MPoly` in the symbols). A term multiplied by `1/D₃·1/D₃·1/D₇` carries `atoms = {3,3,7}`
  (sorted, with multiplicity ⇒ `1/D₃²` is `{3,3}`). `env.atom(aid)` builds a bare `1/D`
  as a monomial with empty exponents and `atoms = {aid}`.

Two monomials are "the same" (and so combine) iff their exponent vectors *and* their id multisets
match. Multiplication adds exponent vectors and merge-sorts the two id lists (`operator*`); both
are cheap `int16_t` operations.

**Why store an id rather than the denominator itself?** Three reasons. (a) At contraction time
`k²` is a *polynomial* in the symbols, not a number, and a reciprocal of a whole polynomial is
neither a symbol nor an integer power — so it cannot fold into the exponent vector and needs its
own slot in the monomial key. (b) The id is a cheap, exact handle: interning each denominator
once into `atomDen` keeps the monomial key small and makes monomial comparison an integer compare
instead of a polynomial compare (this matters — sorting/combining runs over millions of terms).
(c) The id is exactly what cancellation needs: `divThroughMonomialAtoms` looks up `atomDen[aid]`,
and if a bare loop's `k²` has collapsed to a single monomial (e.g. `l1²`, via the unit-vector
identity in `reduce_units`) whose powers the numerator dominates, it absorbs the factor —
subtract exponents, divide the coefficient, drop that id. A denominator that is a genuine
polynomial (a shifted line `k=l−q`) never matches that test, and is handled by the second pass,
`divThroughPolyAtoms`, which trial-divides it into the numerator and drops the id on an exact
division; only an atom surviving *both* passes is lowered to a runtime `inv` env slot. The id is
the link between the on-monomial factor and the denominator table that lets the cancellation
passes decide whether the factor dies or lives.

### Two dressing layers reuse the same trick

`MPoly` interns `1/k²` factors as ids and merges them on multiply. Two *dressing* layers reuse
that exact pattern — a sorted multiset of runtime-call ids carried on each term — for structures
that must stay symbolic to the end of the trace rather than cancel:

* **`DPoly`** (`numeric/dpoly.hpp`) — a **dressed Dirac numerator** like `Mq·δ + Z(p)·γ·p` is a
  *sum* of structures with runtime coefficients. Rather than distribute the diagram into `2^D`
  traces, the front-end keeps it eager and the engine collects one `DPoly`: a map from a
  **dressing-atom multiset** (`DMono`) to the kinematic `MPoly` it multiplies. The coefficient
  *is* an `MPoly`, so `DPoly` reuses `MPoly::operator*`/`+` verbatim and undressed flows are
  byte-identical (their `DPoly` is a single empty-dressing term). Walked through in the
  [step-18](../tutorials/step-18.md).

  **Where the `DPoly` is assembled matters.** The generator's *trace table* is plain `MPoly` even
  for dressed flows: at codegen time each structure×dressing combination is stripped to its bare
  structure, and the dressing is carried alongside as a per-sub-term scalar plus a `DMono`. So
  combinations that share a concrete Dirac structure and differ only in dressing collapse to **one**
  contraction — a 6.2× reduction in distinct traces on the dense quark–gluon flows, and what makes
  full-basis `ZAAqbq1` generate at all. The `DPoly` is then built in phase B, where
  `fold_net_dressed` (`numeric/trace_fold.hpp`) routes each sub-term's scaled `MPoly` into its
  dressing channel. The alternative — collecting the `DPoly` during contraction
  (`numeric_value_dressed_netval`) — is still the reference implementation and what
  `tests/test_dpoly.cpp` grades against, but it re-contracts once per dressing channel.
* **`SUNPoly`** (`network/sun_net.hpp`) — the colour/flavour analogue described in
  [step 3](#step-3-the-colour-fold): a group-diagonal `δ` folds to `Σ_a c_a Z_a` over named
  dressing ids instead of one flavour-blind number.

The decisive contrast with `MPoly`'s atoms: dressing ids **never cancel** — they are opaque
runtime values that ride untouched to the lowering (a `dress` env leaf, evaluated once like an
`inv` leaf). That is why they are a separate layer and not more entries in `Mono::atoms`: mixing
them would tax the hot undressed path and blur the "this factor can cancel" invariant.

## Step 1 — the Dirac trace as matrix products

A closed Dirac chain is traced by **multiplying 4×4 matrices**, not by enumerating index
pairings. `mpoly.hpp` builds each gamma `gammaC(mu)` and each slashed propagator
`slashC(components)` as a 4×4 matrix whose *entries are `MPoly`s* (numeric gamma data, symbolic
momenta), multiplies the chain with `matmul`, and reads off `mtrace`. Because the gamma matrices
are Hermitian and chiral (block-antidiagonal in the Weyl basis), the products stay sparse and
γ5 is free; a chain with an odd number of gammas traces to zero structurally.

When the chain has **free gamma legs** — open Lorentz indices $\gamma^\mu$ that must contract
into the Lorentz network — the trace is not a scalar but a Lorentz tensor; those open legs are
carried through and contracted into the Lorentz network in the next step.

## Step 2 — the Lorentz network reduction

The pure-Lorentz half of a diagram — metrics $\delta_{\mu\nu}$, vectors $a^\mu$, transverse
projectors $P_{\mu\nu}(l)$, Levi-Civita $\varepsilon$ — is assembled with the
`network::` builders (`met`, `vec`, `proj`, `contract`, `add`, `scale`) into a `NetVal`, and
joined with the Dirac tensor from step 1. `numeric/numeric_contract.hpp` contracts it by
**bounded variable elimination** (`eliminate` / `contract_factors`): it sums one shared Lorentz
index at a time over its four values, each element evaluated through the frame's component table,
so the work is bounded by the network's treewidth rather than by $2^{np}$. The result is a single
`MPoly` (`numeric_value` / `numeric_value_netval`), carrying the surviving `1/k²` atoms.

Before elimination, same-momentum projectors that share a dummy index are **fused** by
projector algebra (`fuse_projectors`): idempotency $P\cdot P\to P$ (with $\operatorname{tr}P^T=3$,
$P^L=1$, $P^E=1$, $P^M=2$) and the orthogonalities ($P^L$ ⟂ {$P^T,P^E,P^M$}, $P^E$ ⟂ $P^M$).
This collapses a chain of transverse/longitudinal/electric/magnetic projectors *before* it can
expand in the component basis — e.g. the regulator-dot fold $P^T\!\cdot\partial_t R\cdot P^T$ —
keeping the heaviest gluon traces small.

## Step 3 — the colour fold

A diagram's colour network — structure constants `f^{abc}`, adjoint and fundamental deltas,
generators `T^a` — is just a number once fully contracted. `network/sun_net.hpp`
(`sun_value` / `sun_value_cx`) contracts it numerically at build time over the typed-out
SU(N) tables, so the kernel never carries a colour tensor. Folding colour to a number also lets
the generator group diagrams that share a (factored) dressing coefficient, so identical
monomials from different diagrams merge (cross-diagram collection).

When a network carries a **group-diagonal dressing** (`ntSUNDiagFund` / `ntSUNDiagAdj`), the
colour fold keeps it as a polynomial instead of a single number: `sun_value_dressed` folds it to
$\sum_a c_a Z_a(\text{scale})$ over the named runtime scalar dressings the `spec` selects (the
$\delta$'s `comp2dr` map assigns each surviving component a dressing, and drops the rest). The Dirac
trace is still contracted once; only the colour weight is left component-resolved. (`sun_value_cx`
itself is unchanged, so colour-blind flows are byte-identical.)

## Step 4 — out to the lowering

Each diagram's `MPoly` is handed to `to_genprog` (`numeric_contract.hpp`), which lowers it into
one shared fundamental-symbol environment (`network`'s `GlobalEnv`) via the
[CSE + Horner lowering](cse-and-lowering.md), returning a `GenProg`. The emitted generator drives
that over every diagram through the two parallel phases of `numeric/trace_fold.hpp`
(`contract_traces` then `fold_groups_streaming`), and prints the `trN(const double* f)` functions
plus the `fill` that computes the frame symbols once per call. `mpoly_to_cpp`
(`numeric/numeric_driver.hpp`) renders a momentum component or a propagator denominator as a C++
expression for the `fill` body.

## Why the kernel is compact and fast

By contracting over a fixed frame the engine produces a polynomial in the frame's scalar symbols
directly, with the projector and Wick combinatorics resolved numerically rather than expanded — the
same invariant scalar-product basis a symbolic tracer arrives at, without ever forming the
intermediate blow-ups.

Three things then keep the emitted polynomial small. Two are cancellations, and they are where the
compactness actually comes from:

- **Unit-vector constraints** (`reduce_units`, `numeric/mpoly.hpp`). The loop momentum is written as
  magnitude × unit direction, so `ΣUμ² = 1` holds exactly; rewriting `U_last² → 1 − Σ_{μ<last} Uμ²`
  bounds the last component to power ≤ 1 and collapses the bare-loop denominator `Σ(l·Uμ)² = l²·ΣUμ²`
  to the monomial `l²`, which then cancels.
- **Propagator-denominator cancellation** (`divThroughMonomialAtoms` and `divThroughPolyAtoms`,
  same header). A surviving `1/k²` rides in the monomial key as an *atom*. When its denominator is a
  single monomial it cancels term-by-term. When it is a genuine polynomial — a shifted line
  `k = l − q`, where `k² = l² − 2l·q + q²` — the terms are grouped by their atom multiset and the
  denominator is **trial-divided into the numerator**; an exact division (vanishing remainder) drops
  the atom entirely. Dirac-trace numerators routinely do contain a factor of the very `k²` that sits
  underneath them, so this fires often and removes both a division and every term that cancelled
  against it.
- **Cross-trace CSE** (`"CrossTraceCSE"`, see [CSE and lowering](cse-and-lowering.md)), which shares
  subexpressions across all of a flow's traces rather than lowering each independently.

The second of these is the frame-space analogue of the partial-fractioning (integration-by-parts-like)
step a symbolic tracer does on scalar products *before* substituting a frame. It was previously
assumed a fixed-frame contraction could not do it at all, and that the difference left the generated
kernel slower at runtime in exchange for generating far faster. That is no longer the case: the
cancellation is available in the frame, and with it the numeric kernels are competitive with or faster
than the symbolic ones — see the
[relationship note](../getting_started/overview.md#why-generate-rather-than-evaluate-symbolically)
and `PERFORMANCE.md`. The win is concentrated in quark and ghost loops, where shifted-line
propagators meet matching numerators; pure-gauge flows have little to cancel and are unchanged.
