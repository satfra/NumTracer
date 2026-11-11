# Trace tables — Stage 0 probe (2026-11-11)

Question: can NumTracer precompute the dressing-free traces once per run, instead of evaluating
them inline at every quadrature node of every RHS evaluation?

Test case: `qcd-codes/vacuum/YangMills/full_AAcbc_sigma_split`, float kernels, RTX 4070 Laptop.
Scripts: `numtracer/gpubench/tools/tracetab/` (they read the emitted `*_kernels.hh`, no regeneration).

## Result

All four gates pass, **but only in the homogeneous form below**. The originally planned
"polynomial in k per (external × loop) node" table would need ~10 GB and does not fit the 4070.

| gate | result |
|---|---|
| G1 trace share | **35%** of GPU time per RHS (562 → 364 ms with traces stubbed); dominant flow ZA4tadpoleP1 41%. Ceiling 1.54x. |
| G2 predicted gain | ~1.4–1.5x on the 4070: reconstruction costs ~2 flop per table entry, 1–2% of the inline flops. |
| G3 memory | ≈0.5 GB in float for all running flows (homogeneous form). |
| G4 accuracy | float32 table matches inline **fp64** to within an order of magnitude, and is 10²–10⁸x better than today's float32 traces. |

## 1. Trace share (ablation, nsys, main kernel per flow, median ms per RHS)

Variant: every `trN(fenv)` → constant and `fill` dropped (`gpubench/tools/attrprobe.py::stub_nont`).
The baseline ran 206 RHS evaluations, the stubbed variant 2×6 (it diverges; per-launch timings stay
valid). The stubbed figures are the minima of two runs, interleaved with the baseline.

| flow | baseline | traces stubbed | trace share |
|---|--:|--:|--:|
| ZA4tadpoleP1 | 137.0 | 80.3 | 41% |
| hA | 107.8 | 68.2 | 37% |
| ZA3 | 77.4 | 34.2 | 56% |
| ZAAcbcT | 64.7 | 61.3 | 5% |
| Zcbcbcc | 60.4 | 32.6 | 46% |
| hc | 59.5 | 45.4 | 24% |
| ZAcbc | 29.6 | 21.0 | 29% |
| lamA2 / lamA2Sig / lamcczc | 14.7 / 15.1 / 7.1 | 9.8 / 13.7 / 6.8 | 33 / 9 / 5% |
| **per-RHS sum** | **562** | **364** | **35%** |

The run is fully GPU-bound (GPU kernel time equals residual wall time). ZA4 is 0.1% in this
variant. The agent's "traces are 100% of the cost" does not hold here. The remaining 65% is
interpolator lookups, regulators and the integrand combination.

## 2. The homogeneous form

- After `l1 = k·q̂`, every trace is a polynomial in k and the `inv` slots, with coefficients fixed
  per node (`analyse.py`).
- Every trace in every flow is **homogeneous in (l1, S0) jointly**, of degree d ∈ {0, 2, 4}
  (checked numerically, 0 exceptions). S0 is the only dimensionful external coordinate; S1 and
  SPhi are shape variables.
- Clearing the denominators Q_j = 1/inv_j (mass dimension w_j ∈ {2, 0}) gives

      T = N(k, S0) · ∏ inv_j^{m_j},   N = Σ_n c_n k^n S0^{D−n},   D = d + Σ w_j m_j

  The c_n depend only on (q̂, angles, S1, SPhi), **not on S0 and not on k**.
- Table axes are therefore loop nodes × non-scale external nodes: 864 × 60 = 51 840 for the vertex
  flows, and only the 5 184 loop nodes for ZA4. That is 48x fewer nodes than the naive table.
- The agent's "fixed power of k" is the special case of this with k-relative external grids.
- Reconstruction is overflow-free. With s = max(k, S0) and r = min/max ∈ [0, 1]:

      T = s^d · ∏ (s^{w_j} inv_j)^{m_j} · Horner_r(c)

  (Horner in r, with the coefficient order reversed when k > S0.)

Table entries per node (`E`, sum over the flow's traces) and inline cost:

| flow | traces | inline flop/node | E | float table |
|---|--:|--:|--:|--:|
| ZA4tadpoleP1 | 58 | 103 115 | 584 | 121 MB |
| hA | 19 | 10 942 | 143 | 30 MB |
| ZA3 | 12 | 18 998 | 108 | 22 MB |
| ZAAcbcT | 21 | 6 574 | 134 | 28 MB |
| Zcbcbcc | 12 | 2 849 | 66 | 14 MB |
| hc / ZAcbc | 9 / 6 | 788 / 908 | 28 / 28 | 6 / 6 MB |
| lamA2 / lamA2Sig / lamcczc | 58 / 105 / 15 | 103 068 / 53 709 / 1 845 | 584 / 779 / 60 | ≤ 121 / 161 / 12 MB |
| ZA4 | 5 | 2 760 | 33 | 0.7 MB |

E keeps the inv slots symbolic and counts (k-power, inv-monomial) pairs. In the univariate
cleared form above, E is D+1 per trace. The cleared form is the smaller of the two, except for
1–3 entry traces.

## 3. Accuracy (`validate_h.py --exact`)

- **Reference:** the emitted trace code evaluated in exact rational arithmetic at the double slot
  values.
- **Metric:** max over loop nodes of |T − T_ref| / max |T_ref|, then the max over shape, S0 ∈
  [0.005, 250] and k ∈ [0.034, 1000].

| flow | inline float32 (production) | table float32 | table fp64 | inline fp64 |
|---|--:|--:|--:|--:|
| ZA3 | 1.5e-2 | 8.2e-7 | 1.4e-10 | 3.0e-11 |
| ZA4tadpoleP1 | 1.1e-1 | 5.0e-5 | 5.3e-10 | 3.3e-10 |
| hA | 7.2e-3 | 8.9e-7 | 2.5e-11 | 9.2e-12 |
| Zcbcbcc | **2.4e+4** | 1.9e-4 | 1.9e-4 | 5.0e-5 |
| hc | 2.1e-1 | 1.3e-5 | 6.4e-10 | 1.1e-10 |
| ZAcbc | 7.4e-1 | 4.8e-7 | 1.9e-9 | 2.0e-10 |
| ZAAcbcT | 2.8e-1 | 4.2e-6 | 5.1e-10 | 1.4e-10 |
| ZA4 | 1.3e-6 | 2.4e-6 | 3.9e-13 | 5.0e-16 |
| lamA2 | 1.1e-1 | 5.0e-5 | 5.3e-10 | 3.3e-10 |
| lamA2Sig | 1.4e-4 | 9.0e-6 | 4.6e-13 | 2.4e-13 |
| lamcczc | **7.3e+4** | 2.0e-4 | 2.0e-4 | 9.1e-5 |

**Side finding: the production float32 traces are badly conditioned away from S0 ≈ k**
(`where32.py`, error per S0/k decade):

| S0/k | 1e-6 | 1e-4 | 1e-3 | 1e-2 | 1 | 1e2 | 1e3 |
|---|--:|--:|--:|--:|--:|--:|--:|
| Zcbcbcc | 2e5 | 2e2 | 5 | 3e-2 | 2e-4 | 2e-2 | 2 |
| ZAcbc | 2e-1 | 1e-2 | 1e-3 | 1e-4 | 9e-6 | 7e-4 | 6e-3 |
| ZA4tadpoleP1 | 4e-6 | 4e-6 | 5e-6 | 4e-6 | 5e-5 | 2e-3 | 1e-1 |

- The cause is cancellation in the monomial form: terms of order (k/S0)^n cancel down to the result.
- S0/k ≤ 1e-3 occurs for the low-S0 grid points during most of the flow (k from 1000 down), so
  Zcbcbcc there is noise in float32.
- The table form evaluates in r = min/max ≤ 1 and has no such cancellation.

**Caveat:** near (p−l)² → 0 the factor ∏(s² inv)^m is large while the cleared numerator is small.
One random float32 sample overflowed to inf (ZA3), and none did with the exact reference set. The
emitter has to multiply the factors in an order that cannot overflow, or keep high inv powers as
runtime factors.

## 4. Consequences for Stage 1

- **NumTracer** must know each slot's mass dimension (l1, S0/p: 1; angles and shape: 0; inv: −w).
  From that it can find the scale variable, check homogeneity, and emit c_n as functions of the
  loop and shape slots.
  - A flow with more than one dimensionful external coordinate needs a ratio axis, or falls back
    to inline.
  - Dressed traces (`"DressingCollection" -> True`) stay inline.
- **DiFfRG** must:
  - hand the kernel a (loop node, shape node) table index, since today it passes values only;
  - keep the table across RHS evaluations (it depends only on the grids and x_extent);
  - map threads so that the S0 index varies fastest within a warp, so the 32 lanes broadcast one
    coefficient load.
- **Verified:** q̂ = q/k is fixed for the whole run. The radial box is [0, √x_extent·k], x_extent is
  set once at construction, and `set_k` only rescales the box (`vacuum/integrator_p2_4d_2ang.hh:132-149`).

---

# Stage 1 — prototype (2026-11-11)

All 10 `map()`-based 4D-2-angle flows of `full_AAcbc_sigma_split` were rewritten to read their traces from a per-run
table. It is a post-processing prototype: there is no NumTracer change, and the DiFfRG changes are a header overlay.

| piece | where |
|---|---|
| kernel rewrite (coefficients, reconstruction, `tab()`/`kernel_tab()`/`post_tab()`) | `gpubench/tools/tracetab/emit_tables.py` |
| host-side accuracy test of a rewritten flow | `gpubench/tools/tracetab/host_test.py` |
| DiFfRG: `kernel_has_tab`, `TabTable`, `prepare_tab()`, map() branch, Transform forwarding | `gpubench/tools/tracetab/diffrg_tab.patch` (against the installed headers) |

Runtime switches: `DIFFRG_NO_TAB=1` gives the untouched inline path in the same binary; `DIFFRG_TAB_LAYOUT=loop` selects the
alternative table layout.

## Speed (4070, nsys, main kernel per flow, median ms per RHS)

| flow | inline | table | speedup | ablation ceiling |
|---|--:|--:|--:|--:|
| ZA4tadpoleP1 | 137.0 | 94.8 | 1.44x | 1.71x |
| hA | 107.8 | 71.3 | 1.51x | 1.58x |
| ZA3 | 77.4 | 36.0 | 2.15x | 2.26x |
| ZAAcbcT | 64.7 | 64.0 | 1.01x | 1.05x |
| Zcbcbcc | 60.4 | 40.8 | 1.48x | 1.86x |
| hc | 59.5 | 45.6 | 1.31x | 1.31x |
| ZAcbc | 29.6 | 22.2 | 1.33x | 1.41x |
| **per-RHS sum** | **562** | **400** | **1.40x** | 1.54x |

- **Layout matters.** On CUDA, Kokkos iterates the external index fastest within a warp. A table stored loop-node-fastest
  gave no gain on ZA4tadpoleP1 (1.02x) and 1.25x overall, because of uncoalesced loads. Storing it shape-fastest gives
  the numbers above.
- **One-off host fill:** 49 s per run (ZA4tadpoleP1 37 s). It uses complex long double: 51 840 nodes × (max k-power + 3)
  samples × every trace.
- **Tables:** 73 MB for ZA4tadpoleP1 (354 entries per node) and 128 MB for all 7 running flows. Storing only the live
  k-power range shrinks them by a further ~40% compared with section 2.

## Accuracy

`host_test.py`, max over cells of max_nodes |T − T_ref| / max |T_ref|. The reference is the emitted code evaluated in long double.

| flow | all cells: inline f32 | all cells: table | bulk 1e-2 ≤ S0/k ≤ 1e2: inline f32 | bulk: table |
|---|--:|--:|--:|--:|
| ZA3 | 2.5e-2 | 4.3e-4 | 6.5e-6 | 2.0e-6 |
| ZA4tadpoleP1 | 6.2e-2 | 6.8e-5 | 1.1e-5 | 2.5e-5 |
| Zcbcbcc | 2.8e+4 | 8.5e-2 | 3.7e-3 | 4.1e-5 |
| hA | 3.5e-2 | 1.7e-2 | 4.2e-4 | 2.0e-4 |

**Physics A/B** (same binary, t = 0…1.25, k = 1000 → 287):
- ZA, Zc, ZA4, ZAcbc, ZA4tadpole, ZAAcbcT and m2A agree to ≤ 1e-6; hA and hc to 1e-4.
- ZA3 differs by up to 3e-2 and Zcbcbcc by O(1), both only at the lowest S0 nodes. That is where the inline float32
  traces are wrong (section 3). Above S0 ≈ 1 the differences are ≤ 1e-5.
- This does not prove the table is right at small S0. A double-literal re-emission of those flows would decide it.

## What the prototype exposed for the production design

1. **The coefficients must come from NumTracer's exact polynomials, not from the emitted code.**
   - The float32 emission rounds its constants to 8 digits.
   - Heavily cancelling traces (Zcbcbcc: about 1e11 internal cancellation) then carry residues of ~5e-4 of their largest
     coefficient, in exactly the k-powers where the true coefficient is zero.
   - The prototype drops coefficients outside the live range, which removes the worst of this but cannot repair the
     in-range ones.
   - In the generator, the coefficient of k^n S^(D−n) is an MPoly over the shape and angle slots. Producing it needs
     Q_j as a polynomial in the slots, plus the slot mass dimensions, both passed in from the Wolfram side. That also
     makes the fill a cheap polynomial evaluation instead of a complex DFT.
2. **The live k-power range [n0, n1] is essential.** Without it, rounding noise in identically vanishing coefficients is
   amplified by (k/S0)^n: 1e11 relative errors were seen.
3. **DiFfRG:**
   - The table is keyed on grids and coordinates, never on k or the interpolator epoch.
   - It is stored shape-fastest, matching Kokkos' GPU iteration order.
   - It plugs into the hoisted path as `post_tab()`.
   - Still to add: the other integrators (1/3 angles, `map_points`, finite T).
4. **Still open:**
   - dressed traces (the large QCD kernels);
   - the balance switch, which is a flop/byte number in `emit_tables.py` (`--balance`) and has not been exercised,
     since every trace here pays off;
   - measurements on A100/H100.

## Reference with double-precision constants (2026-11-11)

The ZA3 and Zcbcbcc generators (`gen/gen_*_num.cpp`) were rebuilt to emit double precision (17-digit literals, double
`fill`). The traces and slot layout are identical to the float32 emission. Reference: that code evaluated in long double.
Metric: max over traces and nodes of |T − T_ref| / max |T_ref|, per S0/k decade.

| Zcbcbcc, S0/k | 1e-6 | 1e-5 | 1e-4 | 1e-3 | 1e-2 | 1 | 1e3 |
|---|--:|--:|--:|--:|--:|--:|--:|
| inline float32 (production) | 7.9e+4 | 1.5e+2 | 3.2 | 5.2e-2 | 9.4e-5 | 8.5e-4 | 1.9e-2 |
| inline fp64 | 1.1e-4 | 3.3e-7 | 4.6e-9 | 1.1e-10 | 2.3e-13 | 1.9e-12 | 3.7e-11 |
| table, coefficients from the float32 code | 2.0e-6 | 2.0e-6 | 2.0e-6 | 2.0e-6 | 2.0e-6 | 5.7e-5 | 9.7e-7 |
| table, coefficients from the double code | 3.8e-7 | 2.0e-7 | 4.1e-7 | 3.7e-7 | 3.0e-7 | 1.5e-6 | 9.7e-7 |

ZA3 shows the same pattern but milder: production float32 reaches 9.3e-3 at S0/k = 1e-6, and both tables stay ≤ 1.3e-6.

- **The physics A/B deviations at the lowest S0 nodes are errors of the production float32 traces.** The table run is the
  correct one there.
- With exact constants the table holds float32 precision at every S0/k.
- The production float32 Zcbcbcc kernel is O(1)-wrong for S0/k ≲ 1e-4. At k = 1000 that is S0 ≲ 0.1, about the lowest
  10 of the 48 S0 nodes.
- Whether this moves anything in the IR (Zcbcbcc feeds Zc through the ghost tadpole) is not yet measured. Up to
  t = 1.25, Zc agrees to 1e-6.

## Dressed traces: the large QCD flows (2026-11-11)

`gpubench/tools/tracetab/dressed_count.py` splits each trace into Σ_α D^α · T_α, where D are the dressing slots
(`dr_*`). Each T_α is purely kinematic and tabulable. Entries per node are counted per T_α after clearing denominators.

| flow | traces | inline flop/node | (trace, α) pieces | distinct α | entries/node | flop per entry |
|---|--:|--:|--:|--:|--:|--:|
| Nf2 AAqbq18 ZAAqbq1 | 460 | 737 079 | 6 344 | 161 | 24 661 | 30 |
| Nf2 AAqbq18 ZAAqbq11 | 517 | 858 403 | 6 896 | 161 | 26 351 | 33 |
| Nf2 AAqbq18 ZA4 | 86 | 407 206 | 6 567 | 163 | 24 007 | 17 |
| Nf2 fullAqbq_FP32 ZA4tadpole | 221 | 1 026 267 | 1 877 | 77 | 5 881 | 175 |
| finite-T 3Dquark ZAqbq4 | 351 | 380 576 | 2 100 | 43 | 11 291 | 34 |
| finite-T 3Dquark hPhiL | 1 207 | 384 064 | 6 329 | 56 | 29 824 | 13 |
| finite-T Tmu hPiL | 1 248 | 527 968 | 6 541 | 82 | 30 800 | 17 |

(YangMills, for comparison: 175–550 flop per entry.)

- **Most of the dressed giants are expensive because of their many dressing combinations, not their kinematics.** Each
  T_α is a small polynomial (~4 entries).
- A table then reads one coefficient for every 13–34 flops it saves. In fp64 that is 2–4 flop/byte, below the A100's
  ~6. Per-α tables would therefore be roughly neutral there, and a modest gain only where fp64 is slow (consumer GPUs).
- ZA4tadpole (175 flop per entry) is the exception and looks like the YangMills flows.
- Not measured: whether merging the outer dressing factors of the kernel body (Σ_t G_t · Σ_α D^α T_tα → Σ_β D^β K_β)
  reduces the number of pieces.
- The finite-T flows also need a Matsubara axis and a second scale T (bivariate tables), which this count ignores.
- The balance switch decides this per trace: on these flows it would leave most traces inline.
