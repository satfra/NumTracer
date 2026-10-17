# GPU kernel performance: measured findings

> **Round 2 (2026-08-09/10) supersedes parts of this document** — see
> `../../NUMTRACER_GPU_ROUND2.md`. In particular: the §3c "Nf2 ZA4 = 65% NumTracer" figure
> measures 31% in SASS and the Nf2 IBP ceiling built on it is void; the first real Nf2 solve
> profile puts ZA4 at 39% and the ZAAqbq pair at 44% of GPU time; and the transform-side levers
> (position cache, k-only hoisting, fastInterpLookups) have landed. Raw round-2 data:
> `results/baseline_2026-08-09/NOTES.md`, `results/attr_2026-08-10.csv`.


Hardware for runtime numbers: RTX 4070 Laptop (sm_89), CUDA 12.9. Static/SASS numbers additionally
for sm_80 (A100) and sm_90 (H100) — `ptxas`/`cuobjdump` cross-target without the hardware.
Target deployment is A100/H100, so every conclusion is stated per architecture.

Baseline: `qcd-codes` @ `dd46312`, DiFfRG install built from `DiFfRG_KT` @ `8eb41dc1a`.

---

## 1. The headline: NumTracer's arithmetic is a small minority of the kernel

Measured on the **real production kernels** (full TU: wrapper + interpolators + regulators + traces),
against a probe containing only the NumTracer emission (`fill` + all `trN`):

| flow | NumTracer fp64 | full kernel fp64 | **NumTracer share** |
|---|--:|--:|--:|
| ym ZA    |    50 |  3,359 | **1.5 %** |
| ym Zc    |    46 |  2,981 | **1.5 %** |
| ym ZA3   |   265 |  5,137 | **5.2 %** |
| ym ZAcbc |   263 |  6,727 | **3.9 %** |
| ym ZA4   | 1,811 | 10,641 | **17 %**  |

So for four of five YangMills flows, **generating better kernels from NumTracer cannot buy more than
a few percent** — 83–98.5 % of the fp64 work is the DiFfRG layer around it. Only ZA4 (and the large
Nf2 vertex flows) has enough NumTracer content for emission-side work to matter.

### Where the other 83–98.5 % goes

Subtracting the NumTracer part and dividing by the interpolator count gives **207, 210, 212, 209,
252 fp64 ops per interpolator lookup** across the five flows — consistent to ~2 %. The cause is
`LogarithmicCoordinates1D::backward` (`coordinates.hh:487`):

```cpp
return log1p((y - c) / b) * gem1 / a;
```

A fp64 `log1p` is ~200 SASS instructions on these parts. The spline evaluation itself is ~10.
**Every dressing lookup pays a full fp64 log1p, and nothing shares it.**

---

## 2. The lever that worked: share the coordinate transform  ✅ LANDED

A generated kernel evaluates many dressings at the *same* momentum. nvcc **cannot** CSE the shared
`log1p`, because each interpolator owns its own `coordinates` members and the compiler cannot prove
they are equal across objects. Confirmed by the numbers: ym ZA has 13 lookups over 8 distinct
arguments and still pays 13 transforms.

**Change** (additive, in `physics/interpolation/{spline,linear}_interpolator_1d.hh`): split
`operator()` into

- `index(x)` — the coordinate transform only, depending **solely** on the coordinate system;
- `at(idx)` — clamping, stencil resolution and evaluation;
- `operator()(x) { return at(index(x)); }` — unchanged behaviour for every existing caller.

Clamping deliberately lives in `at()`, not `index()`: clamping is size-dependent, and folding it
into `index()` would silently break sharing between interpolators of different length.

The consumer then hoists one `index()` per distinct argument. `tools/interpshare.py` performs that
rewrite on an emitted `kernel.hh`.

**Superseded (2026-10):** the regex pass described below is replaced by `ntShareGridLookups` (CodegenKernel.m), which also covers the dressing atoms, shares per dressing type instead of all-or-nothing, and takes `"SharedGrids"` for interpolators on one grid.

**Now emitted natively.** `Codegen.m` performs the rewrite (`ntShareInterpIndices`), gated on a new
`MakeNTKernel` option `"ShareInterpolatorIndex"` (default **False**) that `DiFfRG_compat.m` turns on
after checking the dressing types. The layering matters: `index()`/`at()` is a *DiFfRG* property, not
a NumTracer one, so the generic emitter stays backend-agnostic — it takes the callee names from the
`dress` list the caller already supplies rather than pattern-matching DiFfRG's type spelling. Only
the 1-D interpolators are split so far, so a flow reading any dressing off a 2-D/3-D or stack
interpolator opts out entirely (the pass is all-or-nothing per kernel, since the emitted index is an
untyped `const auto`). `NT_NO_INTERP_SHARE=1` is the A/B control.

Verified: the Wolfram pass reproduces the validated Python rewrite on all five YangMills flows,
identical up to `_ixN` numbering and a trailing newline; building from the Wolfram output gives the
same checksums and the same timings.

### Measured (sm_89, production grid orders, median of 11, all checksums **bit-identical**)

| flow | baseline ms | shared ms | speedup | fp64 reduction |
|---|--:|--:|--:|--:|
| ym ZA    |  0.3182 |  0.2683 | **1.19x** | 23.2 % |
| ym Zc    |  0.2933 |  0.2553 | **1.15x** | 21.0 % |
| ym ZA3   |  2.8127 |  2.2238 | **1.26x** | 24.0 % |
| ym ZAcbc |  3.8813 |  3.1855 | **1.22x** | 20.6 % |
| ym ZA4   | 31.6012 | 28.0500 | **1.13x** | 14.5 % |

Runtime tracks the fp64 reduction almost exactly, and the saving per eliminated `log1p` is
**~155 fp64 ops**, consistent to 2 % across all five flows (156/157/154/154/154).

**End to end (L5).** Full YangMills ODE solve, alternating run order to cancel the laptop's thermal
drift (the very first run of a session is ~25 % fast on a cold GPU — an ordering bias that initially
made the *faster* build look slower):

| variant | round 2 | round 3 |
|---|--:|--:|
| baseline | 36.02 s | 36.02 s |
| index sharing | 31.97 s | 32.00 s |

**1.13x on the whole solve, physics bit-identical**: `m2A = 6037.240528758876`,
`Zc(0) = 0.9316141347304615` in both. (Run against `HEAD`'s `model.hh`; see §6.)

**Nf2 — now measured, not modelled.** All 11 flows built and run (the earlier claim that Nf2 could
not build was wrong: `JSONValue` is a deprecated alias for `ConfigTree`). All 11 checksums
bit-identical. fp64 SASS count is the grade here rather than wall clock — a sustained Nf2 sweep
throttles this laptop hard enough to *invert* a wall-clock comparison between rounds:

| flow | base fp64 | shared fp64 | reduction |
|---|--:|--:|--:|
| nf2 ZA      |  4,755 |  3,101 | **34.8 %** |
| nf2 ZA3     |  9,719 |  7,161 | **26.3 %** |
| nf2 ZAqbq1  | 15,151 | 11,941 | 21.2 % |
| nf2 ZAqbq7  | 14,927 | 11,765 | 21.2 % |
| nf2 ZAqbq4  | 17,273 | 13,689 | 20.7 % |
| nf2 Zq      |  3,571 |  2,831 | 20.7 % |
| nf2 Zc      |  2,975 |  2,371 | 20.3 % |
| nf2 ZAcbc   |  6,721 |  5,363 | 20.2 % |
| nf2 ZA4     | 31,221 | 26,969 | 13.6 % |
| nf2 ZAAqbq1 | 15,343 | 13,357 | 12.9 % |
| nf2 ZAAqbq2 | 22,921 | 20,663 | 9.9 % |
| **TOTAL**   | **144,577** | **119,211** | **17.5 %** |

Wall clock agrees where it can be trusted: bracketing both run orders gives ~1.25-1.30x on ZA3 and
~1.08-1.11x on ZA4/ZAAqbq2. Full data in `results/l2_nf2_interpshare.csv`.

> **Reporting caveat.** `interpshare.py`'s printed "log1p saved" percentage counts lookups across
> *all* branches of the `#if <VERDICT>` alternatives that complex flows emit, while only one branch
> is ever compiled. That inflates the count roughly 3x for those flows. The fp64 reduction column
> above is measured on the compiled object and is the number to trust.

**Also covered per request:** `LinearInterpolator1D` got the same split. `FocusedLogCoordinates1D`
(present in the `DiFfRG_KT` source, not yet in the install) benefits automatically — its `backward`
is `log(y)` plus `s_of_u(u)`, i.e. at least as expensive as `log1p`; the split is in the
interpolator, not the coordinate class, so any coordinate system inherits it.
**Not yet done:** `SplineInterpolator1DStack` (2D: stack index + coordinate) and
`LinearInterpolator{2D,3D}` — same idea, but the hoisted quantity is an index tuple.

---

## 2b. `LaunchBounds` on the range policies  ✅ LANDED

DiFfRG's policies carried no launch bounds at all, so ptxas could not see the block size and pinned
registers accordingly. Adding `Kokkos::LaunchBounds<128>` to `KokkosNDRangeHelper`
(`common/kokkos.hh`); `-DDIFFRG_LAUNCH_BOUNDS=0` restores the old behaviour and is the A/B control.

Measured as LB=128 vs LB=0 on identical kernels, alternating, three rounds. Absolute times drift
upward across rounds as the part throttles, so only the **ratio** is meaningful — it is stable:

| flow | r1 | r2 | r3 |
|---|--:|--:|--:|
| ym ZA    | 1.04x | 1.01x | 1.04x |
| ym Zc    | 1.04x | 1.04x | 1.09x |
| ym ZA3   | 1.04x | 1.05x | 1.06x |
| ym ZAcbc | **1.11x** | **1.12x** | **1.13x** |
| ym ZA4   | 1.03x | 1.04x | 1.06x |

Checksums unchanged. Combined with index sharing the flow set runs 1.16-1.33x faster than the
original configuration.

### Rejected: an explicit MDRange tile

`map()` passes no tile, so Kokkos uses its own heuristic. Supplying a tile matched to the launch
bounds *appeared* to give 1.3-1.5x on the 3-angle flows. **It is inert** — that was GPU clock ramp
(see below). Reverted, and deliberately so: Kokkos hard-aborts when the tile product exceeds
LaunchBounds, so an explicit tile converts any future launch-bounds change into a runtime abort for
no measured gain.

### Measurement hazard: GPU clock ramp (this bit three times)

The 4070 Laptop idles at ~210 MHz and needs order a **second** of sustained load to reach its
2445 MHz boost. Whichever configuration is measured first therefore comes out 20-40 % slow. This
manufactured two speedups that did not exist (the tile "1.3-1.5x"), and once made the *faster*
end-to-end build look 25 % *slower*.

`bench.cc` now burns a fixed 2.5 s on ZA4 before timing (`NTB_BURN_MS`). It must be a fixed
wall-clock burn: an "iterate until it stops improving" heuristic is satisfied by a single noisy
iteration and exits while the clock is still ramping — that is exactly how the bogus tile result
survived its first check. Protocol: alternate variants, run >= 2 rounds, discard round 1, and
compare ratios rather than absolute times.

Three independent lines back the interpolator result, which is why it survived this: wall clock,
the fp64 SASS count, and the end-to-end solve. Single-shot wall clock alone establishes nothing
here.

## 3. Measured and killed

### `fill()` CSE — inert at runtime
`fill()` genuinely never goes through the hash-consed SSA the trace bodies use: nf2 ZA4 issues 51
transcendental calls for **4** distinct values, ZAAqbq1/2 issue **103 for 4**, and slots already
stored in `f[4..6]` are re-substituted verbatim into `f[1,2,3,8,9]`. It looks like an obvious win.
It is not:

| flow | arch | base fp64 | fill-CSE fp64 | MUFU base → CSE |
|---|---|--:|--:|---|
| nf2 ZAAqbq2 | sm_80 | 1,680 | 1,685 | 14 → 14 |
| nf2 ZAAqbq1 | sm_80 |   771 |   772 | 12 → 12 |

nvcc's GVN already collapses them. Verified bit-identical over 200k–320k env slots per flow.
It stays a **code-quality** issue (source size, compile time), not a performance lever.
Kept as `tools/fillcse.py` so the negative result is reproducible.

### Fusing the two-phase reduction — not worth it
DiFfRG's `map()` is a `parallel_for` that materialises every sample to global memory, then a
64-team `TeamPolicy` reduction. Measured share of GPU time in phase 2:

| flow | phase 1 | phase 2 | phase-2 share |
|---|--:|--:|--:|
| ym ZA4 | 47.83 ms | 162 µs | **0.3 %** |
| ym ZA3 |  2.75 ms |  26 µs | 0.9 % |
| ym ZA  |   272 µs |  12 µs | 4.1 % |

I expected this to dominate the small η-loop flows. It does not — even for ZA it is 4 %. The
unfused reduction and the 64-team launch are **not** the problem; phase 1 is ~96–99.7 % throughout.

### `CrossTraceCSE` — inert; the design note was right

Implemented, off by default, and never plumbed through `MakeNTKernelDiFfRG` at all. Now plumbed
(`"CrossTraceCSE"`, default False) and measured by regenerating the whole YangMills flow set twice
from `YangMills.wl` — same engine, same everything, only the flag differing. It emits correctly:
per-trace `trN` functions are replaced by one fused `trace_all` (ZAcbc 310 -> 232 lines, 6 trace
functions -> 0).

| flow | base fp64 | CrossTraceCSE fp64 | delta |
|---|--:|--:|--:|
| ym ZA    | 2,601 | 2,605 | +0.2 % |
| ym Zc    | 2,371 | 2,371 |  0.0 % |
| ym ZA3   | 3,935 | 3,935 |  0.0 % |
| ym ZAcbc | 5,377 | 5,357 | -0.4 % |
| ym ZA4   | 9,182 | 9,147 | -0.4 % |

Runtime is inside the noise band (ZA4 27.47 -> 27.33 ms, ZA slightly worse). **Not worth enabling.**
`Codegen.m:3347` already said so — "GlobalCollect ... crossCSE is subsumed by it" — and since
GlobalCollect is on by default, there is nothing left for CrossTraceCSE to find. Nothing in the code
enforces that claim, so it was worth checking rather than assuming; it checks out.

Note it is **not** bit-identical: checksums move in the last 1-2 ULP (ZA `...740e+09` vs
`...741e+09`, ~1e-15 relative). That is the expected consequence of reassociating the arithmetic
into one CSE program, not a defect — but it does mean the flag is not a free swap for an existing
tuned run.

### Regeneration cross-check (closes the emitter's verification gap)

The same regeneration confirms the index-sharing emission end to end: freshly generated kernels
carry `index()`/`at()` with exactly the hoist counts the Python tool predicts (ZA 5 hoists / 10
lookups, ZA3 8/16, ZA4 10/20), and their checksums are **identical** to the post-hoc rewritten
variant. Emitter and prototype agree exactly.

### `-maxrregcount` on the big kernel — inert
On sm_90, nf2 ZA4 ignores `-maxrregcount` entirely (255 and 168 give identical output); only
`__launch_bounds__` moves the allocation, and raising it to 256 makes spilling *worse*.

---

## 3b. The eta loop is a non-problem — and ZA4 is 87 % of the solve

I planned CUDA-graph batching of the eta self-consistency loop on the premise that it issues up to
50 (YangMills) / 80 (Nf2) sequential map() pairs of tiny flows per RHS. **That premise was wrong.**
`eta_iter_max` is a cap that is never approached. nsys over a full YangMills solve, using instance
counts to identify the flows (vertex flows run once per RHS, eta-loop flows once per iteration):

| flow | GPU s | share | launches |
|---|--:|--:|--:|
| ZA4    | 13.611 | **86.9 %** | 494 |
| ZAcbc  |  0.953 | 6.1 % | 494 |
| ZA3    |  0.697 | 4.4 % | 494 |
| ZA     |  0.150 | 1.0 % | 938 |
| Zc     |  0.145 | 0.9 % | 938 |
| phase-2 reduces | 0.110 | 0.7 % | |

938 / 494 = **1.90 eta iterations per RHS**, against a cap of 50. The eta-loop flows are **1.9 %** of
GPU time combined, so CUDA-graph batching would target 1.9 % — not worth building. Dropped.

The useful part of this measurement is the other column: **ZA4 alone is 86.9 % of the solve.** For
YangMills, optimising ZA4 *is* the optimisation problem; everything else is rounding.

## 3c. Rung 6 — the op-count ceiling (the IBP question)

The gauge-flow IBP project is expensive, so the question is what it could possibly be worth. Upper
bound = (NumTracer's share of each kernel's fp64) x (that flow's share of the runtime), i.e. what
you would gain by removing **all** NumTracer arithmetic — IBP would achieve a fraction of it.

| | ceiling |
|---|--:|
| **YangMills** (weighted by the nsys solve profile) | **15.3 %** |
| **QCD_Nf2** (weighted by measured per-flow map() time) | **35.3 %** |

The two are dominated by their ZA4, but for opposite reasons: YangMills ZA4 is 87 % of the solve yet
only 17 % NumTracer, while Nf2 ZA4 is 65 % NumTracer (75 traces against YangMills' 5) and 41 % of its
flow set — contributing 26.6 of the 35.3 points on its own.

**Verdict: IBP is not worth pursuing for YangMills** (15 % ceiling, and IBP delivers only part of it).
For Nf2 there is real headroom, and it is concentrated almost entirely in **ZA4 plus the ZAqbq\*
family** — those are the only flows where emission-side algebra can still move the needle.

## 4. Open, sized, not yet landed

### 4a. The `__noinline__` size gate is dead — and its default is wrong for A100/H100
`grep -c noinline` is **0** across all 16 production flows: these kernels predate the
`ntDeviceEnvPrefix` fix in `Codegen.m`, so `NT_GEN_DEVICE` never reached the generator.
Reproducing the gate on nf2 ZA4 (`tools/noinline.py`):

| threshold | sm_89 spill | sm_90 spill | sm_90 regs | sm_90 occ | fp64 |
|---|--:|--:|--:|--:|--:|
| gate dead (today) | 1,600 B | 11,636 B | 168 | 18 % | 8,589 |
| **min=500 (current default)** | **0** | 1,488 B | 255 | **12 %** | 8,535 |
| min=300 | 0 | **0** | 168 | 18 % | 8,810 (+2.6 %) |
| min=200 | 0 | **0** | 126 | **25 %** | 9,019 (+5 %) |
| min=100 | 0 | 0 |  94 | **32 %** | 9,246 (+7.6 %) |

Two conclusions:
1. On **sm_89** the 500 default is correctly tuned — it eliminates spill outright.
2. On **sm_90** it does not: it leaves 1.5 KB of spill and drops occupancy *below* the ungated
   baseline. Enabling the gate at its current default would **hurt** H100. The threshold has to
   become architecture-dependent (300 or 200 look right for sm_80/sm_90), and that needs runtime
   confirmation on the real part before committing.

### 4a-bis. Register control: `-maxrregcount` is inert, launch bounds are the knob, and neither fixes sm_90

Swept `-maxrregcount` ∈ {none, 168, 128, 96, 64} on nf2 ZA4 and ZAqbq4, sm_80 and sm_90:
**every value produces byte-identical output** — same registers, stack, spill and fp64 count. That
is expected once stated: CUDA gives `__launch_bounds__` precedence over `-maxrregcount`, and since
LaunchBounds is now always present, the launch-bounds *value* is the only register knob. This
closes the "re-run the maxrregcount sweep on the datacenter parts" question — there is nothing there.

Sweeping the launch-bounds value instead (nf2 ZA4):

| LB | sm_80 regs / spill | sm_90 regs / spill |
|---|---|---|
| 64  | 254 / 1,824 | 168 / 11,636 |
| 128 | 254 / 1,824 | 168 / 11,636 |
| 256 | 254 / 1,824 | 128 / 14,504 |
| 512 | 64 / 9,480  | 64 / 22,160 |
| 1024| 64 / 9,424  | 64 / 22,160 |

LB <= 128 is optimal on both; anything larger is strictly worse. And crucially: **no launch-bounds
value avoids the sm_90 spill at all** — the floor is 11.6 KB. The only thing that removed it was the
per-function `__noinline__` gate at threshold ~300 (0 spill, +2.6 % ops; §4a). So on H100 the gate
is not an alternative to register tuning, it is the *only* available lever for this kernel.

### 4b. Occupancy is 12–18 % on the datacenter parts
The repo's standing conclusion ("registers/occupancy are not the bottleneck") was established at
1:64 fp64, where slow arithmetic hides everything. It does not transfer. Cross-arch census:

| flow | sm_80 regs / occ | sm_89 regs / occ | sm_90 regs / occ |
|---|---|---|---|
| nf2 ZA4    | 254 / 12 % | 254 / 16 % | 168 / 18 % (11.6 KB spill) |
| nf2 ZAqbq4 | 252 / 12 % | 252 / 16 % | 254 / 12 % |
| nf2 ZAqbq1 | 252 / 12 % | 252 / 16 % | 252 / 12 % |
| ym ZA4     | 128 / 25 % |  96 / 43 % | 128 / 25 % |

sm_89 tolerates a given register count better than sm_80/sm_90 (1536 vs 2048 max threads/SM), so
the datacenter parts sit *lower* on occupancy for the same kernel.

### 4c. Not investigated
`--use_fast_math` is dead in the build (`DiFfRGTargets.cmake:63` gates it on
`COMPILE_LANGUAGE:CUDA` while every TU compiles as CXX through `nvcc_wrapper`). Left alone
deliberately: it affects fp32 intrinsics and denormals, essentially nothing for fp64.

`CrossTraceCSE` (off by default) — not yet measured on these flows. Note the probe already shows
nvcc doing substantial cross-trace CSE by itself: nf2 ZA4's 20,318 SSA ops compile to 8,589 fp64
SASS instructions (2.4x), so the headroom is smaller than the source suggests.

`powr<N>` with N<0 (`Codegen.m:2537` has no negative branch and would silently return 1.0):
`census.py` checks every production `fill()` — **zero** negative exponents. Latent, not live.

---

## 5. What this means for the original question

> *Can we generate kernels from NumTracer better in order to improve performance?*

Mostly **no, and that is the useful answer**: NumTracer's emission is 1.5–17 % of the fp64 work in
these kernels, and the two emission-side defects that looked like blunders (`fill()` without CSE,
the dead `__noinline__` gate) are respectively inert and a 2.6 %-op tradeoff. The one large,
verified win — 1.13–1.26x measured, up to 26 % fewer fp64 ops — was in DiFfRG's interpolation
layer, not in NumTracer at all.

The emission-side work that *is* still worth doing is confined to the large vertex flows
(ZA4, ZAqbq*, ZAAqbq*) and to the architecture-dependent `__noinline__` threshold.

Landed total on the YangMills flow set: **1.16-1.33x per flow, 1.13x on the full solve**, all with
bit-identical results — none of it from changing how NumTracer lowers algebra.

---

## 6. Environment note (blocks a clean L5 on the working tree)

`vacuum/YangMills/SP/model.hh:67` uses `FocusedLogCoordinates1D`, which exists **only in the
`DiFfRG_KT` source** and not in the installed DiFfRG — and that line is an *uncommitted local edit*
(`HEAD` still uses `LogarithmicCoordinates1D`). So the working tree does not currently build against
the install, and the end-to-end numbers above were taken against `HEAD`'s `model.hh`. Reinstalling
DiFfRG from `DiFfRG_KT` would resolve it.

This makes the index/at split *more* valuable, not less: `FocusedLogCoordinates1D::backward` is
`log(y)` plus `s_of_u(u)`, i.e. at least as expensive as the `log1p` measured here. The split lives
in the interpolator rather than the coordinate class, so it applies unchanged.

## 7. Still not done

- `CrossTraceCSE` still unmeasured on these flows (needs regeneration).
- The `__noinline__` threshold retune needs runtime confirmation on a real A100/H100.
- `SplineInterpolator1DStack` and `LinearInterpolator{2D,3D}` not split (multi-index).
- MDRange tile-size sweep, flat `RangePolicy` variant, CUDA-graph batching of the eta loop.
- Nf2 not wired into `bench.cc`; its interpolator-sharing win is static-only so far.
- The op-count report that justifies or kills the gauge-flow IBP project.
