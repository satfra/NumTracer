# NumTracer — performance summary

> **GPU (2026-08-09).** A full investigation of GPU kernel performance for the two production QCD
> codes lives in `../NUMTRACER_GPU_INVESTIGATION.md`, with the harness and raw data in
> `gpubench/`. Headline: NumTracer's emitted arithmetic is only 1.5-17 % of the fp64 in a
> YangMills kernel; the one large win (1.12-1.24x) was sharing the interpolator coordinate
> transform in DiFfRG, not a codegen change.

Generated numeric kernels vs their FormTracer (FORM) reference, on this machine
(`-O3 -march=native`, `bench_aqbq147` over 200k random points). `num/FORM` is the ratio of
kernel evaluation time; correctness is the max relative error against FORM.

```{note}
The kernel-runtime numbers below were last measured **2026-07-19** and still stand — nothing since
has changed how a kernel evaluates on the host. What has changed is **generation** cost and the GPU
emission default; see "Generation (2026-07-22/23)" at the end of this file. Do not read the
2026-07-19 date as covering those.
```

## Current headline (2026-07-19): the numeric backend is FASTER than FORM

`bench_aqbq147`, 200k points, `taskset`-pinned:

| struct  | FORM ns | numeric ns | num/FORM | was  | rel-err vs FORM |
|---------|---------|------------|----------|------|-----------------|
| ZAqbq1  | 1237.9  | **1184.0** | **0.96x** | 2.90x | 3.04e-09 |
| ZAqbq4  | 1130.0  | **1117.6** | **0.99x** | 3.02x | 9.74e-11 |
| ZAqbq7  | 2040.0  | **1262.8** | **0.62x** | 1.84x | 2.78e-09 |

Three changes landed to get there, in order of size:

1. **Multi-term denominator cancellation** (`divThroughPolyAtoms`, `numeric/mpoly.hpp`) — the big one.
   `divThroughMonomialAtoms` only ever cancelled an inverse atom whose denominator was a single
   *monomial*, so a shifted line `k = l − q` (with `k²` multi-term) kept its `1/k²` atom even when the
   numerator carried a factor of that same `k²` — which Dirac-trace numerators routinely do. Terms are
   grouped by their atom multiset and each multi-term denominator is trial-divided in; an exact
   division (vanishing remainder, same relative tolerance as the noise prune) drops the atom. On
   ZAqbq1_147 Mq-in: 384 of 1,168 trial divisions are exact, monomials 13,269 → 4,832, SSA
   33,775 → 7,649, runtime 3021 → 1102 ns. **Accuracy improves** (one division and a pile of
   cancelling terms are gone). Generation gets *faster*. `NT_GEN_NO_POLYDIV=1` disables.
2. **Coefficient snapping** (`snap_coeff`, `numeric/numeric_contract.hpp`) — 2.90x → 2.64x, free.

This supersedes the 1.35/1.37/1.05 row below, which had been unchanged since the initial commit and
predates `tests/gen/gen_qcd_aqbq147_numeric.wls` — it never measured this flow.

### Emitted-size effect per flow (`kernels.hh`, before → after the denominator cancellation)

The win is concentrated where shifted-line propagators meet matching numerators — quark and ghost
loops. Pure-gauge and small unit-test flows are untouched.

| flow | before | after | |
|---|---|---|---|
| ZAqbq1_147 / ZAqbq4_147 / ZAqbq7_147 | 1.92M / 1.81M / 2.01M | 431K / 389K / 451K | **0.22x** |
| ZAqbq1, ZAcbc | 74,902 / 50,931 | 20,482 / 13,944 | **0.27x** |
| ZAAqbq1Small (+Ref) | 88,213 / 89,341 | 42,660 / 43,000 | 0.48x |
| Zq | 4,375 | 2,693 | 0.62x |
| ZA4, YM_ZA4 | 179,714 / 163,486 | 137,266 / 125,572 | 0.76x |
| ZA4_147 | 5,122,566 | 4,261,663 | 0.83x |
| ZA3, ZA3_147 | 21,249 / 232,498 | 18,501 / 211,930 | 0.87x / 0.91x |
| ZA, EtaPiL, Discdirac, Glu_quark, Pion_quark, unit flows | — | unchanged | 1.00x |

`ZA3_147` runtime is unchanged at 1.01x num/FORM (`bench_147`) — a pure-gauge 3-point has little to
cancel. `ctest` is 36/36 including `codegen_regen` after regenerating every flow above.

## Symbolic dressing collection (default-on, no gate)

Keeping a dressed propagator-numerator sum eager as one collected `DPoly` trace (instead of
distributing into `2^D` diagrams) is now a **pure win** — faster than the previous GlobalCollect
path on both runtime and generation, so it stays default with no gate.

**(historical — see the current headline above; these ratios were measured on a flow definition that
predates `tests/gen/gen_qcd_aqbq147_numeric.wls` and do not reproduce today.)**

Fundamental-colour quark–gluon vertex `ZAqbq{1,4,7}_147` (struct 1 / 4 / 7):

| metric                | GlobalCollect (old) | collection, unfixed | **collection (current)** |
|-----------------------|---------------------|---------------------|--------------------------|
| runtime, num/FORM     | 1.42 / 1.45 / 1.10  | 1.74 / 1.72 / 1.35  | **1.35 / 1.37 / 1.05**   |
| generation (struct 1) | 27 s                | 45 s                | **27 s**                 |
| rel-err vs FORM       | —                   | —                   | 2.5e-7 / 5.1e-9 / 3.5e-8 |

What made it a win (both fixes in `mathematica/Codegen.m`, source-only — non-dressed kernels
regenerate byte-identical):

- **Runtime:** the imaginary-vanishing probe (`numericImagProbeRealQ`) is now dressed-aware. It
  previously failed to compile against a dressed `fill()` (missing the `dr_<id>` dressing-atom
  args), so every dressed kernel was conservatively kept complex and computed-then-discarded its
  imaginary half. It now stubs and passes the dressing atoms, confirms `Im≈0`, and emits the
  lossless **RePart** double kernel (re/im split; the dead imaginary halves are DCE'd by the
  compiler) — exactly what the distributed path already got.
- **Generation:** the dressed generator's main TU compiles at `-O1` (env `NT_GEN_MAIN_OPT`) — it
  runs once (~1 s) so `-O2`'s extra cost on the `DPoly`/`dch`/`dsl` templates was wasted (28.7 s
  vs 17.5 s, identical run). The `lc/dc` net-CSE also now runs for the dressed path.

`to_genprog(DPoly)` prunes round-off **per dressing channel** (each kinematic poly against its own
max), not against one global tolerance — a channel reweighted small at generation time can dominate
at runtime once a regulator suppresses the others, so a global tolerance could silently delete its
genuine terms.

## Other flows

All committed numeric kernels match their FORM / equivalence oracle — the full ctest
suite passes, including the per-flavour / per-component / disconnected flows
(`flow_flavour_split`, `flow_flavour_ingroup`, `flow_gluon_condensate`, `discdirac_num`,
`flow_sigl_dirac`). For indicative ratios of the larger gluon kernels (quark box `ZA4` ≈ 1.15×,
three-gluon `ZA3` 1/4/7 within a small factor of FORM), see the per-flow `compare_*` drivers and
`tests/refshim/`. The generation pipeline (the disposable build-time generator program) is
described in the repo's top-level `README.md` (the *Architecture* section).

## Generation (2026-07-22/23) — cost of *producing* a kernel

Distinct from everything above, which is how a kernel *runs*. These changes moved the generator, not
the emitted arithmetic; the host kernels are byte-identical across all of them.

| change | effect |
|---|---|
| Packed-exponent `MonoExp` + `NarrowAlloc` (`numeric/mpoly.hpp`) | stored polynomial term `pair<Mono,Cx>` **128 B → 72 B** (−43.75%), ~30% faster on a contraction microbench. Monomial compare becomes two integer compares. |
| Blocked `MPoly::operator*` above a scratch threshold | bounds the transient peak of a "collapsing" multiply (where the result is far smaller than the `n·m` scratch) instead of materialising all of it. |
| `MPolyFactory::scaled` | scaling a trace by a constant no longer pays `operator*`'s scratch and its dead `std::sort`: 3.7–4.1× on that operation, bit-identical. |
| **Lever (b)** — dedup dressed traces on structure | the trace table loses its dressing dimension: combinations sharing a concrete Dirac structure collapse to one contraction (**6.2×** fewer on the dense flows). Full-basis `ZAAqbq1` went from *not generating at all* (OOM at 20 GB after 8.66 CPU-hours, still inside phase A) to **2.4 GB at W=8**. |
| Vertex collection (`ntDiracSlot`, opt-in) | `za3_147` generates in **8.9 s** collected vs **13.3 s** distributed. Off by default — see the codegen internals doc for why. |

The generator-side dials that trade RAM against time (`NT_GEN_MAXW`, `NT_GEN_MAXW_B`,
`NT_GEN_GROUP_WINDOW`, `NT_GEN_MEMO_MAX`) are tabulated in
[documentation/internals/codegen.md](documentation/internals/codegen.md).

**GPU emission default changed.** Device trace functions above ~500 SSA instructions are now emitted
out-of-line, which is both faster and cheaper to compile at that size; below it they stay inlined.
The measurements and the two overrides are in [tests/gpu/README.md](tests/gpu/README.md).

## Reproducing

```bash
# regenerate a flow (needs wolframscript + FunKit, and FORM for the YM/QCD tensor bases):
wolfram -script tests/gen/gen_qcd_aqbq147_numeric.wls
# runtime + rel-err vs FORM (not a ctest — compile by hand):
g++ -O3 -march=native -std=c++20 -Iinclude -Itests/gen -Itests/refshim \
    tests/refshim/bench_aqbq147.cpp -o /tmp/b && /tmp/b
# the full numeric gate:
ctest --test-dir build --output-on-failure
```
