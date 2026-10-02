# step-20: Making a hard flow generate

*Builds on: [step-09](step-09.md), [step-15](step-15.md) · Built on by: [step-22](step-22.md) ·
Tags: `performance`, `build` · **Tier C** (the full fRG toolchain)*

```{admonition} What this step needs
:class: important
Nothing new to install — but the flows worth applying it to are the production ones, so you need a
project like `QCD_Nf2`. The numbers quoted are measured; re-measure before quoting them yourself
([step-22](step-22.md) on stale numbers).
```

## Introduction

[step-14](step-14.md) ended with a table spanning two orders of magnitude in emitted size, from a
1.7 kB ghost propagator to a 467 kB four-gluon vertex — and in `with_mesons`, an 81 MB
`lambda1L3D` with 989 trace functions. Somewhere on that curve, generation stops being instant and
starts being a thing you have to manage.

This step is the diagnostic workflow: **which symptom points at which lever**. It is deliberately
organised that way rather than as a list of knobs, because the knobs are already listed — the
authoritative table is in
[Front-end & codegen](../internals/codegen.md), and it moves faster than this page can.

### The shape of the problem

Generation is a two-phase C++ program that the Wolfram front-end writes, compiles and runs:

| Phase | What it does | Scales with | Parallel over |
|---|---|---|---|
| **A** | contracts each *distinct* trace | number of distinct traces × per-trace size | traces |
| **B** | folds each net, then **lowers** it to straight-line SSA | number of nets × polynomial size | nets |

Both are parallel, and their resource profiles are opposite:

* **Phase A is memory-bound.** Peak RSS is roughly `W ×` the per-contraction working set, so the
  worker count *is* the RAM ceiling.
* **Phase B is compute-bound** — the CSE and Horner lowering is real work — but each worker also
  holds one full contraction, so on dense flows it is a RAM dial too.

That asymmetry is the single most important thing to internalise, because it makes the obvious move
wrong.

## Symptom → lever

### "It OOMs"

Cap the phase workers. From Wolfram:

```mathematica
SetNumTracerThreads[nA, nB]   (* or SetNumTracerThreads[n] for both *)
GetNumTracerThreads[]         (* returns {nA, nB}; Automatic = hardware concurrency *)
```

These export `NT_GEN_MAXW` / `NT_GEN_MAXW_B`, which the generator reads at run time. They only ever
*lower* the count below hardware concurrency.

```{admonition} Do NOT set a low global cap at the top of a multi-flow script
:class: warning
Phase B is *both* the fold **and** the lowering. Throttling it starves a compute-bound flow:

> the pure-gauge `ZA4` / `ZA4`-full-basis: its phase-B + lower is **~600 s at `nB = 2` vs tens of
> seconds unthrottled**.

Throttle **only around the RAM-heavy flows** — in practice the four-point two-quark-two-gluon
vertices — and leave `nB` high or unset everywhere else. A `SetNumTracerThreads[2]` at the top of a
`.wl` will make every other flow in the file an order of magnitude slower for no benefit.
```

Under `"Offline" -> True` there is no `Run[]` to inherit environment variables through, so the
values in force at emit time are **recorded per flow** into `numtrace.json` (`"maxw"` / `"maxw_b"`,
0 = unset) and re-applied by the build, where they lower — never raise — the build's own `-jN`. You
saw `"maxw": 8` in [step-15](step-15.md)'s manifest.

Other RAM levers, when the thread cap is not enough:

* `NT_GEN_GROUP_WINDOW` — how many net polynomials may be in flight at once (the streaming-fold
  bound).
* `NT_GEN_MEMO_MAX` — cap how many contracted traces stay resident; the rest are recomputed in
  phase B. This is the dial the memory-bound dressed flows use.

### "The generator takes forever to compile"

Generation compiles a C++ program before it runs one, and on a big flow that program is large. Two
levers:

* **`NT_GEN_MAIN_OPT`** — the optimisation level on the *main* translation unit. Measured:
  **`-O2` is over 600 s where `-O1` is ~6 s**, for no runtime benefit; `-O1` is the default recorded
  in `numtrace.json`, and `-O0` is catastrophically slow on dressed flows for a different reason
  (the generator then *runs* slowly). This is one of the few settings where the middle option is
  the right one.
* **Chunking.** Big net-builder literal tables are split into top-level functions of bounded size
  (`ntChunkDef`, `$ntDefChunk`), and the chunks are bin-packed into several
  translation units (`gen_*_u0.cpp … _uN.cpp`) so no single emitted function blows the compiler.
  This is automatic; you only meet it when a flow is large enough to produce eight units, as
  [step-15](step-15.md)'s `ZA` manifest does.

### "Generation runs forever"

Reduce the number of traces. This is a *physics* lever and it beats every resource knob:

* **Restrict the basis** (`TBRestrictBasis`) or **restrict the trace** (`dropAqbq47`,
  `dropFourQuark` — [step-14](step-14.md)). Measured on the four-quark internal quark–gluon
  tensors: **133 / 439 / 1349 diagrams and 13 / 38 / 87 s** for `{1}` / `{1,7}` / `{1,4,7}`.
* **Vertex collection** (`NT_VERTEX_COLLECT=1`, [step-18](step-18.md)) — a measured **36.8×**
  diagram reduction on the quark–gluon vertex, and 8.9 s vs 13.3 s on `za3_147`. Opt-in because on
  high-multiplicity flows it OOMs instead.
* **Sub-term dedup** contracts each *distinct* trace once into a shared table, and is on
  unconditionally. It matters more than it sounds: measured 5.1× and 7.5× redundancy on dense flows,
  and it also flattens the parallel work items — per-net scheduling had pinned utilisation at ~33 %.
  `NT_GEN_NO_DEDUP=1` is the A/B control, not a tuning knob.

```{admonition} Some flows are intractable, and that is information
:class: note
The full-basis `ZAAqbq1` two-quark-two-gluon vertex once failed to complete phase A at all —
measured at **> 8.66 CPU-hours and > 20 GB**, with the wall time explained by 51 198 distinct dense
trace contractions. No RAM knob fixes that; the fix was to reduce the *trace count*. If your flow is
in that regime, stop turning dials and go back to the truncation.
```

### "The emitted kernel is enormous"

Usually it is not your fault and there is nothing to do — a four-point vertex is big. But two things
are worth knowing:

* **Denominator cancellation** (`divThroughPolyAtoms`) is on by default and is the single largest
  size lever, concentrated where shifted-line propagators meet matching numerators — quark and ghost
  loops. Measured on `kernels.hh`: `ZAqbq{1,4,7}_147` went **1.92M/1.81M/2.01M → 431K/389K/451K**
  (0.22×), `ZAqbq1`/`ZAcbc` 0.27×, `ZA4` 0.76×. Accuracy *improved* and generation got faster.
  `NT_GEN_NO_POLYDIV=1` disables it — as a control, never as a setting.

## Offline generation

For anything expensive, do not generate inside the Wolfram script. `MakeNTKernelDiFfRG` defaults to
`"Offline" -> True`, which writes the generator sources plus a `numtrace.json` with
`"generated": 0` and stops. The consumer's CMake then does the work:

```cmake
numtracer_add_numtrace(QCDFlows ${CMAKE_CURRENT_SOURCE_DIR})
```

Three reasons this is the right default:

1. **It parallelises correctly.** The module chains flows **one at a time**, so each gets the whole
   `-jN`. Generation is memory-bound; running eight flows concurrently is how you OOM.
2. **It is incremental.** Only flows with `"generated": 0` are rebuilt.
3. **It puts the cost where the resources are.** Deriving a flow is a laptop job; generating a
   four-point vertex is not.

## Measuring, rather than guessing

`NT_GEN_PROFILE` gives per-phase timing and RSS from the generator (`=2` adds a per-wave RSS trace),
and `NT_GEN_VERBOSE` turns on the Wolfram-side `[prof]`/`[cse]`/`[time]` lines.

```{admonition} The profiling trap
:class: warning
`NT_GEN_PROFILE` output is **lost** through `wolframscript`'s `Run[]` — the generator is a child
process whose stdout does not come back. You must run the compiled `gen_<namespace>` binary
**directly**:

```bash
NT_GEN_PROFILE=2 ./gen_za_qcd_num
```

And note that the phase-A summary line prints only on *completion*, so a run you kill mid-phase-A
looks like it produced no profiling at all — which is exactly the run you most wanted to profile.
Under `"Offline" -> True` the binary is right there in the build tree, which is another reason to
prefer it.
```

A worked ordering for a flow that will not generate:

1. **Count the diagrams first** — `Print[Length[ntk[[1]]["Diagrams"]]]`. Free, and it tells you
   whether you have a resource problem or a truncation problem.
2. **Profile the binary directly** to find out whether phase A or phase B dominates.
3. **If A dominates**: reduce the trace count (basis restriction, vertex collection). Cap `nA` only
   to survive, not to go faster.
4. **If B dominates**: leave `nB` alone. Look at `NT_GEN_GROUP_WINDOW` for RAM and at the kernel
   size for the real cause.
5. **If the compile dominates**: check `main_opt` is `-O1` and that the chunking produced several
   units.

## Things that were measured and do NOT work

Worth stating explicitly, because they are the natural next ideas:

* **Parallelising the Wolfram front-end over diagram terms.** Measured **19× slower** — the
  cross-process marshalling dominates. The front-end phase is ~98 % FunKit's `MakeCppHeader`, so
  there is nothing there for NumTracer to parallelise anyway. Process-level splitting is the only
  real lever.
* **Memoising the per-diagram Dirac/colour compilation.** A measured loss.
* **Turning `nB` down globally** to save RAM — see the 600 s figure above.

## Possibilities for extensions

1. **Find the knee.** Generate `ZAcbc`, `ZA3`, `ZAqbq1`, `ZA4` in turn, timing each and recording
   peak RSS. Plot against the emitted size. The curve tells you where on it your own flow sits.

2. **Profile one flow properly.** Build a generator offline, run it directly with
   `NT_GEN_PROFILE=2`, and identify the phase split. Then change one physics lever (a basis
   restriction) and re-measure. You now have a cost model for your truncation.

3. **Reproduce a dead end.** Set `NT_GEN_MAXW_B=2` on a pure-gauge `ZA4` and confirm the slowdown
   is as dramatic as documented. Doing this once is what stops you reaching for a global cap later.

4. **Check the size lever.** Regenerate a quark-loop flow with `NT_GEN_NO_POLYDIV=1` and compare
   `kernels.hh` sizes against the 0.22× figure. Confirm the accuracy claim too — the cancelled form
   should be *more* accurate, not less.

5. **Watch the offline path work.** Set `"generated": 0` in one `numtrace.json`, rebuild, and
   observe that exactly one flow regenerates and the flag flips back.

Next, [step-21](step-21.md): the same kernels on a GPU, and what the generator does about traces
that are complex-valued.
