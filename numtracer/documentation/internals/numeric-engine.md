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

`numeric/mpoly.hpp` defines `MPoly`, a multivariate polynomial whose variables are the frame's
*scalar symbols* (momentum components / scalar products) and which carries surviving inverse
atoms `1/k²` (propagator denominators) on each monomial. It is the arithmetic currency of the
whole engine: addition, multiplication, monomial-level cancellation of `k²·(1/k²)`
(`divThroughMonomialAtoms`), and collection of equal-momentum propagator groups
(`reduce_units`). Everything below produces and combines `MPoly`s.

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
monomials from different diagrams merge — the cross-diagram collection FORM does.

When a network carries a **group-diagonal dressing** (`ntSUNDiagFund` / `ntSUNDiagAdj`), the
colour fold keeps it as a polynomial instead of a single number: `sun_value_dressed` folds it to
$\sum_a c_a Z_a$ over the runtime per-component dressing leaves $Z[a](\text{scale})$. The Dirac
trace is still contracted once; only the colour weight is left component-resolved. (`sun_value_cx`
itself is unchanged, so colour-blind flows are byte-identical.)

## Step 4 — out to the lowering

Each diagram's `MPoly` is handed to `to_genprog` (`numeric_contract.hpp`), which lowers it into
one shared fundamental-symbol environment (`network`'s `GlobalEnv`) via the
[CSE + Horner lowering](cse-and-lowering.md), returning a `GenProg`. The driver
`numeric/numeric_driver.hpp` (`run_numeric`) does this for every diagram and the generator prints
the `trN(const double* f)` functions plus the `fill` that computes the frame symbols once per
call. `mpoly_to_cpp` renders a momentum component or a propagator denominator as a C++ expression
for the `fill` body.

## Why this reaches FORM

FORM keeps scalar products symbolic and emits a flat polynomial in them. The numeric engine
reaches the *same* basis: by contracting over a fixed frame it produces a polynomial in the
frame's scalar symbols directly, with the projector and Wick combinatorics resolved numerically
rather than expanded. The residual runtime gap to FORM is the one optimisation FORM does that the
fixed-frame contraction does not — integration-by-parts on the symbolic scalar products *before*
the frame is substituted — which on the heaviest vertices leaves NumTracer ~1.2–1.6× FORM at
runtime, in exchange for generating the kernel ~80–175× faster (see the README).
