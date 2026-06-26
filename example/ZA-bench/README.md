# ZA-bench — NumTracer vs FormTracer on the full-QCD (AAqbq) gluon self-energy

This isolates the **ZA** flow of the full QCD `147 + AAqbq` truncation — the flow whose
NumTracer *"Generate"* output blew up — and benchmarks it against the FormTracer analogue on
**generation time, emitted size, compile time, binary size, and runtime ns/eval**, plus a
**faithfulness gate** that the global-collection rework is output-preserving.

## Why this exists

| | NumTracer (per-diagram "Generate") | FormTracer |
|---|---|---|
| emitted | `kernel.hh` 375 KB + `kernels.hh` 212 KB | `kernel.hh` 19 KB |
| structure | **273** trace functions, **6 952** temporaries | **1** collected kernel, **37** temporaries |
| runtime | 589 ns/eval (3.26× FORM) | 181 ns/eval |

**Root cause** (see the repo plan / `numtracer/README.md`): the "Generate" backend strips the
dressing functions out *before* tracing and emits one `trN(fenv)` per (diagram × colour channel);
diagrams fuse only when their dressing coefficient is *identical* and quark-loop fundamental-colour
diagrams are forced singleton → 273 singletons, and the C++-level `Σ` combiner can't collect across
the opaque trace calls. FORM instead keeps dressings symbolic and collects the **whole** integrand
into one factored expression. The rework makes NumTracer collect the same way (dressings become
inert symbols in one global invariant-basis polynomial).

## Layout

- `bench/za_bench_shim.hpp` — standalone shim (no DiFfRG/Kokkos/CUDA): math names, interpolator
  stubs, a fixed regulator, and the 15 QCD dressings used identically by both paths.
- `bench/shiminc/DiFfRG/physics/{interpolation,physics}.hh` — stubs so the DiFfRG-targeted
  generated kernels `#include "DiFfRG/physics/..."` resolve to the shim.
- `bench/compare_za_full.cpp` — driver, compiled once per NumTracer variant via
  `-DNT_HDR -DNT_CLASS -DNT_NAME`; reports NT-vs-Form agreement + ns/eval and dumps NT values
  on a deterministic grid for the faithfulness gate.
- `generated/{form,nt,nt_before}/` — generated kernels (populated by the scripts below).
- `gen_form_za.sh`, `gen_nt_za.sh`, `run_bench.sh` — live regeneration + the full report.

## Run

**Self-contained (uses whatever is in `generated/`)** — no toolchain needed:
```bash
cmake -S . -B build && cmake --build build
ctest --test-dir build --output-on-failure
```

**Full live regeneration + report** (needs `wolfram` + FunKit + DiFfRG and FORM/`tform`):
```bash
./run_bench.sh                 # regenerate both, compile, run, size/time table, gate
SKIP_GEN_FORM=1 ./run_bench.sh # reuse the cached FormTracer kernel
SKIP_GEN_NT=1 SKIP_GEN_FORM=1 ./run_bench.sh   # benchmark already-generated kernels only
```

## Correctness model

- **Primary gate** (`run_bench.sh`): the reworked `nt` kernel must reproduce the pre-rework
  `nt_before` kernel to **< 1e-10** — proving global collection is an output-preserving algebraic
  regrouping (the "no regression" guarantee).
- **Informational**: NT-vs-FormTracer agreement. The `QCD_NumTracer` derivation is a
  work-in-progress and currently differs from the `QCD_vacuum_147` FORM derivation by ~20%
  (ratio 0.8–1.3 pointwise). That gap is **pre-existing** (a derivation/convention difference
  between two independent projects), not introduced by the collection rework.
