# GPU integration tests (CUDA)

Full 4D loop integrals of the generated ZA3/ZA4 invariant-basis kernels over a grid of 64
log-spaced external momenta, run on the GPU as DiFfRG-style two-phase map/reduce
(`include/numtracer/cuda/integrator.cuh`):

1. **phase 1** — ONE grid-stride launch evaluates `prefactor · l1³ · Πwᵢ · Re kernel(...)` at
   every (grid point, quadrature node) pair into a device buffer (one contiguous segment per
   grid point);
2. **phase 2** — a PURE reduction: `cub::DeviceSegmentedReduce::Sum`, one segment per grid
   point, plus the loop-independent `constant()` per point.

Both phases are timed with `cudaEvent`s, the whole `map()` with `std::chrono`; each test
prints per-phase times and evals/s.

## Quadrature (GSL, host-side — `quadrature.hpp`)

Mirrors DiFfRG's `Integrator_p2_4D_{2,3}ang`:

| axis | rule | interval |
|------|------|----------|
| `l1` (radial) | Gauss-Legendre | `[0, √x_extent·k]`, `x_extent = 2` |
| `cos1` | Gauss-Chebyshev-2 | `[-1, 1]` — the `√(1-cos1²)` jacobian is **in the weights**, not the integrand |
| `cos2` | Gauss-Legendre | `[-1, 1]` |
| `phi` (ZA4 only) | Gauss-Legendre | `[0, 2π]` |

Prefactor: ZA3 `2π/(2π)⁴` (trivial phi), ZA4 `1/(2π)⁴`.

## Kernels

The tests consume the committed `__host__ __device__` kernel variants
`tests/gen/ZA{3,4}_gpu_kernel{,s}.hh` — generated with the codegen option
`"Decorator" -> "static __host__ __device__ inline"` (distinct class/namespace
`ZA3_gpu_kernel`/`za3_gpu` so they link next to the host kernels without ODR collisions).
Regenerate with:

```bash
wolfram -script numtracer/tests/gen/gen_qcd_za3_gpu.wls
wolfram -script numtracer/tests/gen/gen_qcd_za4_gpu.wls
```

The committed host generators can also re-emit the straight-line header of either variant
directly: `./build/gen_za4_inv -n za4_gpu -d 'static __host__ __device__ inline'`.

`shim.hpp` here is the host/device twin of `tests/refshim/shim.hpp` (identical formulas;
`NT_DH` decoration, id-dispatch `Fn` instead of host function pointers). The GPU targets put
`tests/gpu` on the include path instead of `tests/refshim`, so BOTH the `.cu` TU and the CPU
reference TU pick it up.

## Validation

`ref_za3.cpp` / `ref_za4.cpp` compute the same quadrature sums on the CPU with the committed
HOST kernels (independently compiled path, long-double accumulation); the GPU result must
match to < 1e-8 relative on every grid point (mechanics bugs are O(1); legitimate
FMA-contraction/summation-order differences sit around 1e-13). Validation and timing both run
at the production orders 32×8ⁿ (radial 32, each angle 8 — the typical DiFfRG production
sizes).

## Kokkos twins

`integrate_za{3,4}_kokkos.cpp` + `integrator_kokkos.hpp` are the same tests expressed in
Kokkos (CUDA backend) instead of raw CUDA + cub — same generated kernels, quadrature, orders
and buffer layout; phase 2 is a per-grid-point team reduction (DiFfRG `map()` style). They
build against the DiFfRG-bundled Kokkos (`NT_KOKKOS_ROOT`, default
`~/.local/share/DiFfRG/bundled`) via `kokkos_compilation(TARGET ...)` →
`nvcc_wrapper(g++-14)`. Enable with `-DNUMTRACER_KOKKOS_TESTS=ON`; ctest labels `gpu;kokkos`.

Measured comparison (RTX 4070 Laptop, orders 32×8ⁿ, 64 p-points, best of 5, back-to-back):

| | raw CUDA | Kokkos | |
|---|---|---|---|
| ZA3 phase 1 | 1.3–2.0 ms (65–100 Me/s) | 1.4–2.1 ms (63–96 Me/s) | parity (scatter = GPU clock ramp on the 2 ms workload) |
| ZA4 phase 1 | 30.1 ms (34.8 Me/s) | 30.0 ms (35.0 Me/s) | parity |
| phase 2 (reduce) | 7 / 13 µs (cub) | 109 / 557 µs (team) | ~15–40× slower, still negligible |
| phase-1 kernel regs (ZA3/ZA4) | 224 / 254 | 150 / 254 | all spill-free, 40 B stack |
| ZA4 TU compile (nvcc) | ~8 min | 2:47 | host-side g++ pass dominates the raw TU (527 s; ptxas only 2.5 s) |

**Launch bounds matter**: without `Kokkos::LaunchBounds<128>` on the phase-1 RangePolicy the
Kokkos ZA4 ran ~14% slower (37.6 ms) — and the raw CUDA kernel WITHOUT `__launch_bounds__(128)`
showed the same deficit (39.0 ms), so the gap was codegen pinning, not framework overhead.
With matched launch bounds the two implementations are equivalent in runtime.

Kokkos quirks encountered (already worked around in the sources): extended `__device__`
lambdas cannot first-capture a variable inside `if constexpr` (arity dispatch hoisted into
`eval_at`), and the integrand is passed as a trivially-copyable functor instead of a
parameter pack (packs in extended lambdas are fragile under nvcc).

## vs the FormTracer/DiFfRG production kernels (YangMills/SP example, sm_89)

Like-for-like phase-1 (pure integrand-to-memory) kernels — DiFfRG `QuadratureIntegrator::map`
with the FormTracer-generated YM kernels vs `nt_phase1` with the NumTracer QCD kernels:

| phase-1 kernel | registers | stack frame | local-mem instr (STL+LDL) | fp64 arith instr |
|---|---|---|---|---|
| NumTracer ZA3 (QCD) | 224 | 40 B | 5 | 4814 |
| FormTracer ZA3 (YM) | 255 (cap) | 1928 B | 717 (10.6%) | 2198 |
| NumTracer ZA4 (QCD) | 254 | 40 B | 5 | 9032 |
| FormTracer ZA4 (YM) | 255 (cap) | 2712 B | 963 (8.7%) | 4244 |

The FormTracer `map()` phase-1 numbers are essentially identical to the fused `get()` kernels
(1960/2704 B stack), and DiFfRG's phase-2 team reduction is clean (52–54 regs, 0 stack) — so
the 255-reg cap + ~2 KB/thread spill-via-stack churn is the monolithic FormTracer integrand
expression itself, not the integrator machinery. NumTracer's fill/trN-structured kernels stay
below the cap with no local-memory traffic while carrying the additional QCD quark-loop
diagrams. Both designs end at the same 16.7 % register-limited occupancy, which the
spill-experiment above shows is not the bottleneck (fp64 pipe is).

## Large traces: does the flat `const double sN = …` emission spill, and is that a problem?

The ZA3/ZA4 kernels above stay under the register cap, but four-gluon-with-quark kernels are
far bigger (`ZAqbq1_147` 23 k lines / 248 traces; `ZA4_147` 147 k lines / 977 traces). Probing
the trace arithmetic in isolation — a bare `__global__` that calls `fill()` then accumulates
every `trN` with a distinct runtime weight (so cross-trace *cancellation* can't happen but
legitimate cross-trace *CSE* still can, exactly like `nt_phase1`) — on the RTX 4070 / nvcc 12.9:

| kernel | lines | traces | regs | stack | spill st/ld | local-mem instr | fp64 |
|---|--:|--:|--:|--:|--:|--:|--:|
| ZA4_num | 5 059 | 17 | 56 | 40 B | 0 | 5 | 2 065 |
| ZA3_147 | 7 877 | 165 | 252 | 0 | 0 | 0 | 2 761 |
| ZAqbq1_147 | 23 315 | 248 | 254 | 2.4 KB | 2.7 / 2.9 KB | 699 | 9 506 |
| ZA4_147 | 147 200 | 977 | 255 | 17 KB | 38 / 39 KB | 9 646 | 56 842 |

So past ~20 k lines the all-inlined kernel **does** exceed the 255-register file and spills —
heavily at four-gluon scale. **The spill is benign, though**: emitting each `trN` as a
`__noinline__ __device__` function eliminates it completely (global kernel → 56 regs, every
trace ≤ 56 regs, 0 spill) yet is *slower* at every scale —

| best-of phase-1 | inlined (spills) | `__noinline__` (0 spill) |
|---|--:|--:|
| ZAqbq1_147 | **64.6 ms** | 77.2 ms (+19 %) |
| ZA4_147 | **202 ms** | 217 ms (+7 %) |

— because the regime is fp64-pipe-bound, not occupancy/memory-bound (three confirmations: a
`-maxrregcount` sweep 255→128→96→64 on the spilling kernel only gets slower; out-of-lining does
*more* total work — it sheds cross-trace CSE, fp64 9 506 → 11 269, and pushes `f[]` into local
memory; and the spill latency hides under the arithmetic pipe). On the CPU side the point is
even cleaner: `tr0` of ZA4_num (1186 fp ops) compiles to a **0-byte stack frame, 0 spills** —
Horner form has bounded register pressure by construction, so each flat-SSA trace is
register-clean regardless of length.

**Verdict: the flat all-inlined emission is the runtime-optimal one; don't change the default.**
The only lever that would speed these kernels up is fewer arithmetic ops (the IBP/collection
problem), which is orthogonal to emission style; a register-pressure-aware emission order is
not worth pursuing because register pressure is not the binding constraint.

The one cost of all-inlining is *compile* time: the 147 k-line ZA4_147 TU takes ~2.5 min /
2.7 GB through nvcc (host pass dominates). For users who hit nvcc time/RAM limits on a monster
kernel, the codegen exposes an opt-in escape hatch — set `NT_GEN_NOINLINE_TRACES=1` when running
the generator and the per-diagram `trN` are emitted `__attribute__((noinline))` (`fill`/`powr`
stay inline). That roughly halves the ZA4_147 compile time and RAM, at the ~7–19 % runtime cost
above. It is **not** the default and leaves the emitted bytes unchanged when unset.

## Build & run

```bash
cmake -S numtracer -B numtracer/build-gpu -DNUMTRACER_GPU_TESTS=ON -DNUMTRACER_KOKKOS_TESTS=ON \
      -DCMAKE_CUDA_HOST_COMPILER=g++-14    # nvcc 12.9 rejects the system gcc 16
cmake --build numtracer/build-gpu -j
ctest --test-dir numtracer/build-gpu -L gpu --output-on-failure     # all four
ctest --test-dir numtracer/build-gpu -L kokkos --output-on-failure  # Kokkos twins only
```

clang CUDA works as an alternative: `-DCMAKE_CUDA_COMPILER=clang++` (the CMake block sets
`-fconstexpr-steps` etc. for it).

## nvcc pitfall (why the codegen avoids `decltype`)

nvcc's frontend loses `decltype(var)` when `var`'s initializer has a DEDUCED return type
(`et::contract_all`) inside a class-template member — the gcc host pass then sees an
undeclared `__T0`. The kernel wrappers therefore evaluate compile-time colour factors via
argument deduction (`numtracer::et::escalar(_colnetN, nullptr)`) instead of
`expr::evaluate<escalar_t<decltype(_colnetN)>>` — same constant, nvcc-safe.
