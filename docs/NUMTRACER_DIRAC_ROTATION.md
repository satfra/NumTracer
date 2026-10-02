# Dirac fold: cyclic rotation of the trace (design, 2026-10-02)

Status: option A LANDED ("perf: walk each Dirac trace from its cheapest cyclic start"). It changes the emitted bits (rounding only); it was graded on
the oracles, and the committed kernels get one regeneration together with D4.

Measured result (instructions vs the previous commit):
- za4_147 ×0.375 (phase A 46 → 22 s)
- aqbq147 ×0.61
- za3_147 ×0.81
- ZAAqbq1Small(Ref) ×0.96–0.97

Kernels: maximum relative deviation 7e-13 (za4_147); ctest 61/61 on 32 regenerated flows; za4_147
peak RSS +16%.

## Where the time is (HEAD 0658d59)

| flow | phase A share in `numeric_dirac` |
|---|---|
| za4_147 | 85% |
| ZAAqbq1SmallRef | 25% |

The cost inside `numeric_dirac` is the polynomial arithmetic of the 2×2 Weyl-block products.
After the four bit-identical MPoly levers, no structural waste is left in that arithmetic itself.

## The DFS and its cost

`numeric_dirac(chain)` computes the trace of one closed spinor loop for all 4^f assignments of its
f free γ legs. It walks the chain in token order with a DFS: free legs branch 4 ways, and the
running block product is shared along each prefix.

A token that sits after k free legs is therefore multiplied 4^k times. In Weyl blocks the costs differ:
- γ^μ (and C) is a signed permutation, i.e. a cheap scaled copy.
- A slash `p̸` is a full polynomial multiply.

Cost model: C = Σ over slash/C/Comm tokens of 4^k(token).

The chain is a trace, so any cyclic rotation gives the same value: tr(A·B) = tr(B·A). This holds for
the full 4×4 product, including γ5, C, transposed tokens and Comm. Rotating the chain so that the
slashes come before the free legs moves the expensive multiplies to small k.

Example, with quark-box-like alternation and S a propagator slash:

    γ¹ S₁ γ² S₂ γ³ S₃ γ⁴ S₄   →  C = 4 + 16 + 64 + 256 = 340
    S₄ γ¹ S₁ γ² S₂ γ³ S₃ γ⁴   →  C = 1 + 4 + 16 + 64   = 85

## Measured on the real chains

The probe computed C for the current order and for the best rotation over all chains of a flow. It
was temporary and has been removed.

| flow | free legs f | C(best) / C(current) |
|---|---|---|
| za4_147 | 4 (65%), 5 (35%) | 0.11 |
| aqbq147_1 | 3–6 | 0.02 |
| ZAAqbq1SmallRef | 4, 5 | 0.20 |
| ZAAqbq1Small | 4, 5 | 0.27 |
| za3_147 | 3 | 0.41 |

The model ignores two things, so the real gain will be smaller than 1/ratio:
- The running product grows along the chain.
- The γ copies of that product also scale with 4^k.

Rough expectation for za4_147 phase A: 2–4×.

## Proposed change (option A, ~40 lines in `numeric_dirac`)

1. Pick the start token that minimises C (O(L²) over L ≤ ~20 tokens; negligible).
2. Walk the rotated chain. Free legs are then met in rotated order, so each leaf's flat index must
   still be built in the ORIGINAL `freeLegs` order: give each leg the digit weight
   4^(f−1−originalPosition). `F.ids` is unchanged.
3. Nothing else changes: parity, `started`, γ5, C, Comm and transposed tokens are all
   rotation-invariant.

An optional later refinement is option B, the leaf fusion: tr(M·γ^μ) at the last free leg is a
signed sum of two entries per block, so the final block multiply per leaf is not needed. It is
small; measure it after A.

Not proposed now: moving the spinor indices into the Lorentz elimination network. That is a large
redesign.

## Why the bits change, and how to grade it

A different multiplication association rounds the coefficients differently, at the 1e-16 relative
level. The monomial set is mathematically the same. The noise prune (1e-9 relative) absorbs any
inexact cancellation.

Grading:
1. Oracle tests: `ctest -LE codegen` plus `NT_RUN_CODEGEN=1 ctest -L codegen` (FORM and dense
   oracles), with tolerances unchanged.
2. On the 7 bench flows, compare the base and new emitted kernels numerically at random points. The
   expected relative difference is ~1e-14 or smaller.
3. Performance: phase-A instructions at W=1 (W=32 for za4_147), as before.
   - Keep the change if za4_147 instructions are ≤ 0.6× and every oracle passes.

## Regeneration

The committed `tests/gen` kernels drift once, so they need a dedicated regen commit. D4
(`-ffp-contract=off` on the library) has the same consequence. Doing both before a single regen
gives one attributable drift commit instead of two.
