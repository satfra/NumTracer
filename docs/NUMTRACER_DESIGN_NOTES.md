# NumTracer design notes

Measurements and history moved out of the header comments in `numtracer/include/numtracer/`. The
headers keep the invariant each item protects; this file keeps the numbers behind it, so a measured
dead end is not retried. Sections follow the header they came from.

## core/config.hpp — `NT_THROW` and `-fno-exceptions`

The net-builder generator TUs are compiled with `-fno-exceptions` (mathematica/CodegenBuild.m).
Emitting exception-cleanup landing pads for their tens of thousands of destructible temporaries
dominated their `-O0` compile: turning exceptions off cut a representative unit from 15.3 s to 1.3 s.

## codegen/lower.hpp — Horner lowering

- `choose_pivot` tally: the dense id-indexed count is O(Σ|vp| + V). The previous per-factor linear
  scan was O(Σ|vp|·V); at V≈40 distinct ids it was the dominant term of the whole lowering
  (O(N·|vp|·V²) per trace). The first-appearance tie-break reproduces that scan exactly.
- `partition_pivot`: dividing the pivot out in place removed O(N·V) heap allocations per trace.
- Common-factor extraction at the Horner combine (`A*X + A*Y*Z -> A*(X+Y*Z)`) was implemented and
  rejected. At a combine the two addends are `pivot^n * W` and `ratio * Wo`; their only factors are
  the pivot chain (absent from the other by construction) and the branch sums, so a shared factor
  needs `W == Wo`. Measured 1.000x / 0.996x / 0.999x ops on ZAqbq{1,4,7}_147, no runtime gain on CPU
  or GPU. The structure exists elsewhere (27% of add nodes on the production with_mesons/ZA4 kernel,
  4-9% at this combine), so it would need a post-pass over the finished instruction stream.
- Rounding the normalised Horner ratio would merge coefficient noise (`-31.999999999999002` -> `-32`)
  and buys ~10% more sharing. The noise is ~1e-14 relative, so any rounding coarse enough to merge it
  perturbs by ~1e-14, and unlike `snap_coeff` (once per monomial, at a leaf) it lands at every
  interior node, where the ~1e6 cancellation amplification acts. ZA3_147 against
  compare_za3_147_num's 1e-9 gate: unrounded 5.6e-11, rounded to 13 significant digits 2.7e-08
  (failing). No intermediate setting helps; recovering the 10% needs exact rational coefficients
  upstream.

## codegen/gen.hpp — Horner ordering sweep, normalisation guard, emission

- `make_orderings` builds only the first `maxOrders` variants: each is a full deep copy of the
  monomial set, and `best_into` sweeps 3 (500-2000 monomials) or 1 (>2000), so building all 8 was
  waste. On the dense 1/4/7 traces (tens of thousands of monomials) the 8-way sweep dominated
  generation while changing the op count only ~1%: the op count tracks the monomial count, not the
  pivot order.
- `kNormGuardMax = 2000` (cost both the normalised and the plain Horner lowering only below it):
  on the small fixtures nearly every trace grows ~5-10% under normalisation (73 of 75 on the
  vendored nf2 ZA4); on the big flows nearly every trace shrinks (138 of 143 on with_mesons/ZA4).
  Costing both arms is why no regenerated kernel came out larger than its baseline.
- `EmitPlan` fma folding and single-use-constant inlining used to sit behind
  `NT_GEN_NO_FMA` / `NT_GEN_NO_CONST_INLINE`. Both measured as wins; the hatches were referenced
  nowhere and spelled `getenv(x) != nullptr`, so `NT_GEN_NO_FMA=0` turned fma folding OFF. To A/B
  them again, delete the branch rather than reintroduce a variable.

### The device `__noinline__` gate (`edetail::eff_decor`)

- Original sweep (sm_89, RTX 4070, GPU runtime, 12 flows): out-of-lining above ~500 SSA
  instructions per function looked runtime-faster (ZA4 655 instr/fn "0.73x", ZAAqbq1 772 "0.70x")
  and ~3x cheaper to compile; below it inlining won (ZAqbq1_147's 108 functions of ~125
  instructions ran 1.66x faster inlined). Spills reached 26 KB/thread on the dense 4-point flows.
- The ZA4 figure was later re-measured at −4.8%, wrong by ~5x, so the sweep is not a sound basis for
  500. A static SASS sweep on nf2 ZA4 (sm_90): ungated leaves 11,636 B of spill at 18% occupancy;
  min=300 leaves none; min=200 reaches 25% and min=100 32% occupancy, at +2.6% / +5% / +7.6% ops.
  That argues for 200-300 on datacenter parts. Not changed because the one flow claimed to lose
  from out-of-lining (ZAqbq1_147) sits below 500 and that claim was never reproduced.
- The gate was dead twice. First, it sniffed the decorator for `__device__`, but `ntKokkosDecor`
  (CodegenKernel.m) rewrites `__host__ __device__ inline` to `KOKKOS_INLINE_FUNCTION` before emission
  and DiFfRG_compat.m passes the Kokkos spelling directly (2026-08-08: QCD_Nf2/no_mesons ZA4 had 75
  trace functions, 7 over 500 lines, `tr0` at 2281, zero `noinline`). Sniffing `KOKKOS_` instead is
  wrong too: those macros expand to plain `inline` on host-only Kokkos builds. Hence `NT_GEN_DEVICE`.
  Second, the offline path (`cmake -P` NumTracerNumtraceRun.cmake) inherits nothing of the Wolfram
  kernel's environment; closed 2026-08-11 by the `"device"` field in numtrace.json (a manifest
  without it reads as false).
- nvcc does not force-inline everything: on SP_EM ZA4 (55 traces, 435k SSA) it out-of-lines 26 on
  its own; the gate takes that to 51.

## numeric/trace_fold.hpp — phase A / phase B

- Trace redundancy on dense flows: 30,807 contractions for 6,041 distinct traces (5.1x), and
  246,456 for 32,784 (7.5x). Hence the distinct-trace table.
- Per-net scheduling could not use the machine: sub-terms per net are skewed (max 2880 vs median 27),
  so the single biggest net exceeded the ideal per-thread load and pinned utilisation at ~33%
  regardless of core count. Hence the flat work list.
- Tree fold vs left fold: worst relative deviation of the emitted kernel ~7e-15 (last ulp).
- Phase-B granularity: making the GROUP the work item needs no per-net storage, but serialises every
  net of a group onto one thread. ZA4 (54 nets in 6 groups, one of 27 nets): generator run
  0.78 s -> 3.8 s. Parallelise over nets.
- Streaming phase B: the batch `fold_nets` + per-group sum kept both sets live and never released,
  20+ GB on the dense 4-point flows (488 nets x ~41 MB) against a 20 MB trace table.
- `check_group_partition` used to run only under `NT_GEN_PROFILE` and only warned, so in practice it
  never ran. Before it was made fatal, a full regeneration of all 29 DEFAULT_FLOWS reported zero
  violations.

## numeric/mpoly.hpp — polynomial storage and the reductions

- Inline `Mono` buffers: with `std::vector` members every product term in `operator*` was a heap
  allocation; profiling the generator run showed malloc/free at ~17% of total time.
- Monomial size chain: the stored term `pair<Mono,Cx>` went 128 B -> 112 B (`NarrowAlloc` on the
  atom list, 24 B -> 16 B header) -> 72 B (packed 128-bit `MonoExp`), −43.75% overall; `Mono` is 56 B.
  The packed order reproduces the previous element-wise lexicographic order exactly, so the change
  was byte-identical on every kernel.
- `MPolyScratch` inline storage: `from_scratch` is called tens of millions of times with a mean of
  ~4 terms.
- `MPoly::terms` inline storage (capacity 4): ~2% faster, but sizeof(Mat4) ~0.5 KB -> ~8.7 KB and peak
  RSS +~20%. Rejected.
- Rvalue `operator+`: the const& overload copied the surviving side on 70-85% of all calls.
- Rvalue `divThroughMonomialAtoms` / `reduce_units`: the pass-through fires on ~95% / ~88% of calls.
- An unconditional `reserve()` in the atom merge of `operator*` cost ~7% of the run.
- `divThroughPolyAtoms` (exact multi-term atom division), ZAqbq1_147 Mq-in (real part): monomials
  13,269 -> 4,832, atom factors 1,168 -> 832, fused SSA 33,775 -> 7,649 (0.23x); 384 of 1,168 trial
  divisions are exact. The lead pre-filter rejects >= 2/3 of trials (100% on some flows).

## numeric/numeric_contract.hpp — Dirac fold, elimination, lowering

- Dirac DFS prefix sharing: re-folding each of the 4^f free-leg assignments from scratch (the former
  `NT_DIRAC_FLAT=1` A/B control, removed) took 3.64 B mul2 calls against the DFS's 522 M on a
  production flow. Entries are bit-identical either way.
- Lazy transposed-gamma blocks: building ~40 MPoly by value on every `numeric_dirac` call cost +1.4%
  instructions on an all-untransposed fold.
- `eliminate` seeded from the first factor instead of `constant(1) * e`: with 2-4 factors per group
  the identity multiply was 1/4 to 1/2 of all elimination multiplies.
- Atom-denominator table aliasing: the per-trace reduction used to deep-copy the whole table plus one
  copy per entry on 10^5 traces x combinations, although the generated driver had already reduced it.
- The rank-1 electric-projector split used to sit behind `NT_NO_RANK1_PROJE=1`, which nothing set.

### `kCoeffSnapDigits` (coefficient snapping before lowering)

ZAqbq1_147 (Mq in), 200k points, against the FormTracer oracle:

| digits | multiplies | distinct literals | ns/eval | NT/FORM | rel. error vs FORM |
|-------:|-----------:|------------------:|--------:|--------:|-------------------:|
| off    | 29,395     | 2931              | 3508    | 2.90x   | 4.26e-09           |
| 15     | 27,049     | 740               | 3371    | 2.75x   | 2.98e-09           |
| 14     | 26,096     | 378               | 3283    | 2.64x   | 6.43e-09 (default) |
| 12     | 25,988     | 338               | 3220    | 2.63x   | 4.55e-06           |

14 captures essentially the whole speed win at no accuracy cost; 12 buys a further 2% for three
orders of magnitude of accuracy.
