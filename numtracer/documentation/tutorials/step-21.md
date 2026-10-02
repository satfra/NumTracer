# step-21: GPU kernels and complex traces

*Builds on: [step-09](step-09.md) · Built on by: — · Tags: `gpu`, `cuda`, `complex` ·
**Tier B** (a Wolfram kernel and FunKit; CUDA only for the integration tests)*

## Introduction

Two topics that seem unrelated and are not: what it takes to run a generated kernel on a GPU, and
what the generator does about traces that come out complex. They belong together because both are
consequences of the same design choice — **emit flat, branch-free, division-free real arithmetic** —
and because the complex-trace machinery is what makes a `double`-returning device kernel possible at
all.

## Part 1: the GPU

### The change is one option

```mathematica
MakeNTKernel[ntk, …, "Decorator" -> "static __host__ __device__ inline"]
```

That is it. The emitted arithmetic is unchanged — there is nothing in
`const double s0 = …; const double s1 = …; return sN;` that a device compiler can object to. No
branches, no loops, no divisions, no function calls, no allocation. The straight-line form
[step-05](step-05.md) printed *is* the GPU form.

The committed device variants live in `tests/gen/ZA{3,4}_gpu_kernel{,s}.hh` and are generated with a
distinct class and namespace (`ZA3_gpu_kernel` / `za3_gpu`) so they link alongside the host kernels
without ODR collisions. `"DeviceTarget"` and the Kokkos decorator spelling are handled by the
emitter; you do not normally set them by hand.

### The integrator

`tests/gpu/integrator.cuh` runs a full 4D loop integral over a grid of external momenta as a
DiFfRG-style **two-phase map/reduce**:

1. **Phase 1** — one grid-stride launch evaluates
   $\text{prefactor}\cdot l_1^3\cdot\prod_i w_i\cdot\operatorname{Re}\text{kernel}(\dots)$ at every
   (grid point, quadrature node) pair into a device buffer, one contiguous segment per grid point.
2. **Phase 2** — a pure reduction: `cub::DeviceSegmentedReduce::Sum`, one segment per grid point,
   plus the loop-independent `constant()` per point ([step-13](step-13.md)).

The quadrature mirrors DiFfRG's:

| axis | rule | interval |
|---|---|---|
| `l1` (radial) | Gauss–Legendre | $[0, \sqrt{x_{\text{extent}}}\,k]$ |
| `cos1` | Gauss–**Chebyshev-2** | $[-1,1]$ |
| `cos2` | Gauss–Legendre | $[-1,1]$ |
| `phi` (4-point only) | Gauss–Legendre | $[0, 2\pi]$ |

```{admonition} The Chebyshev row is not a detail
:class: warning
The $\sqrt{1-\cos_1^2}$ Jacobian is **in the weights**, not in the integrand. Using Gauss–Legendre
on that axis silently integrates against the wrong measure. This exact mistake — a $d=3$ measure
used for a $d=4$ integral — was the root cause of a ~2× gauge over-screening in a production
propagator integrator. A kernel can be perfectly correct and still give wrong physics if the
quadrature it is handed disagrees with the frame it was generated for.
```

### What is actually worth optimising

```{admonition} Measure before you optimise the emission
:class: important
On the two production QCD codes, **NumTracer's emitted arithmetic is only 1.5–17 % of the fp64 in a
Yang–Mills kernel.** The rest is interpolator lookups. The one large measured win (1.12–1.24×) was
sharing the interpolator coordinate transform *in DiFfRG*, not a codegen change.

If your GPU kernel is slow, the emission is probably not why.
```

Two more measured cautions:

* **Occupancy and spill are architecture-specific.** Results tuned on sm_89 (fp64 1:64) **do not
  transfer** to A100/H100 (fp64 1:2). There the production flows sit at 12–18 % occupancy and a
  four-point flow spills 11.6 kB/thread on sm_90 versus 1.6 kB on sm_89. The per-function
  `__noinline__` size gate helps, at a threshold different from the one tuned on sm_89, and no
  `launch_bounds` or `-maxrregcount` setting substitutes for it.
* **GPU clock ramp will lie to you.** An idle GPU clocks down hard, so whichever variant runs
  *first* looks 20–40 % slow. That artefact has faked a "1.3–1.5×" tile speedup and inverted an
  end-to-end result. Always: fixed-duration burn-in, alternate the variants, and compare ratios —
  never a single A-then-B run.

## Part 2: complex traces

### Why a trace is complex

The Euclidean gammas are Hermitian in the chiral basis and carry explicit $\pm i$ entries; SU($N$)
structure constants enter with factors of $i$; projectors can too. So the natural result of
contracting a diagram is a **complex** polynomial, and `Cx` is the coefficient type the engine
carries throughout ([step-03](step-03.md)).

The *physical* integrand, though, is real. For a flow that is manifestly real, the imaginary parts
cancel between diagrams — but only after the sum, and the generator has to decide what to emit
before it can know that.

### Three bodies and a probe

The generator emits **three** candidate kernel bodies:

| body | what it computes | valid when |
|---|---|---|
| the untouched complex form | $\sum c_i\,\mathrm{tr}_i$ | always |
| `ntPureIntegrand` | drops the imaginary coefficients | $\sum \operatorname{Im}(c_i)\,\mathrm{tr}_i \approx 0$ |
| `ntRePartIntegrand` | $\operatorname{Re}\big(\sum c_i \mathrm{tr}_i\big) = \sum[\operatorname{Re}(c_i)\,\mathrm{ntRe}_i - \operatorname{Im}(c_i)\,\mathrm{ntIm}_i]$ | always, but cheaper than complex |

It then emits a small **probe** program, runs it against the just-generated traces at random points
(4000 by default, tolerance $10^{-9}$), decides which body is valid, and writes the verdict to
`numtrace_verdict.hh`. The kernel selects with `#if <NS>_VERDICT`.

So the decision is made **by measurement on the actual generated traces**, not by a symbolic
argument that could be wrong. `"RealProbe" -> True` is the default and should stay on.

```{admonition} `PruneRealTraces` is hazardous with the probe
:class: warning
`"PruneRealTraces"` (default `False`) drops traces the generator believes are purely real. Combined
with the probe it has produced an O(1)-wrong kernel: the probe then runs on the **pruned** trace
set, sees no imaginary part left, and flips the verdict to "pure" — validating exactly the
assumption the pruning made. A self-confirming check is worse than no check.

```

```{admonition} `ntRePartIntegrand` is not the naive rewrite
:class: note
Taking the real part looks like it should be "replace each coefficient by its real part". It is not:
$\operatorname{Re}(c\cdot t) = \operatorname{Re}(c)\operatorname{Re}(t) -
\operatorname{Im}(c)\operatorname{Im}(t)$, and the naive rewrite silently drops the second term. The
generator builds it via an $i$-substitution for that reason. If you ever hand-edit a complex kernel,
this is the trap.
```

You met the footprint already: [step-10](step-10.md)'s kernel contained
`const auto _interp1 = ntRe(tr0(fenv));`.

### A measured curiosity

On the flows in this repository, **every** complex trace turned out to be purely imaginary — 218 of
218. The real halves are fully Horner-lowered and then multiplied by zero. Skipping them is
tempting and is *not* byte-identical (it changes the CSE's `best_into` ordering), which is why it
has not simply been switched on.

## Possibilities for extensions

1. **Emit a device kernel.** Take your [step-06](step-06.md) network and re-emit with
   `"Decorator" -> "static __host__ __device__ inline"`. Diff the two headers: the decorator changes
   and nothing else does. That is the whole claim of this page, verifiable in one `diff`.

2. **Read the verdict header.** Generate a complex flow and inspect `numtrace_verdict.hh` and the
   `#if` in `kernel.hh`. Then set `"RealProbe" -> False` and see what is emitted instead.

3. **Reproduce the clock-ramp artefact.** Time variant A then variant B on an idle GPU; then B then
   A. If your conclusion flips, you have reproduced the artefact — and you now know why the
   fixed-duration burn-in is not optional.

4. **Check the quadrature measure.** Integrate a function with a known closed form over `cos1` with
   Chebyshev-2 weights and again with Legendre. The ratio is the error the warning above describes.

5. **Find the fp64 fraction.** For one of your kernels, count the emitted arithmetic operations and
   compare against the interpolator lookups per call. If emission is 2 % of the work, optimising it
   is 2 % of the opportunity.

6. **Break the probe.** Set `"PruneRealTraces" -> True` on a genuinely complex flow and compare the
   result against the unpruned kernel. This is a controlled reproduction of a real O(1) error, and
   it is worth doing once.

Next, [step-22](step-22.md): how to know that any of this is right.
