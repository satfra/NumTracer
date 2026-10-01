# step-09: Codegen options — what each one changes

*Builds on: [step-06](step-06.md) · Built on by: [step-20](step-20.md), [step-21](step-21.md) ·
Tags: `codegen`, `options` · **Tier B** (a Wolfram kernel and FunKit)*

## Introduction

`MakeNTKernel` has 36 options. Listing them would be a reference page; this is a tutorial, so
instead we emit **the same network four ways** and look at what actually changed in the generated
C++.

The organising fact is this: **almost every option is semantics-preserving.** It changes how the
kernel is spelled — how diagrams are grouped, how subexpressions are shared, what namespace and
runtime support it targets — without changing the number it returns. So the natural test is that
all variants agree bit for bit, and that is what the checker asserts. The differences you are meant
to study are in the emitted source.

### The full option list

From `Codegen.m`:

```
"ComputeType","Name","Namespace","Dressings","ScalarParams","ADParams","ParameterOrder",
"Decorator","DeviceTarget","IncludeDir","RunGenerator","FullParallel","AngleDefs","CrossTraceCSE",
"Components","SymbolDefs","RuntimeInclude","ExtraIncludes","KernelNamespace","SupportNamespace",
"DressingType","ShareInterpolatorIndex","HoistLoopConstLookups","RegulatorTemplate",
"RegulatorAlias","RealProbe","PruneRealTraces","ComplexRuntimeProjection","ComplexEndProjection",
"RealOutput","Constant","Offline","CoordinateArgs","MatsubaraVar","DecayingRegulators",
"MatsubaraFiniteExtent"
```

They fall into seven groups:

| Group | Options | Covered in |
|---|---|---|
| **Identity** — what the kernel is called and where it lives | `Name`, `Namespace`, `KernelNamespace`, `SupportNamespace`, `Decorator`, `DeviceTarget` | here, and [step-21](step-21.md) for the device decorator |
| **Interface** — what it takes and returns | `Dressings`, `DressingType`, `ScalarParams`, `ADParams`, `Constant`, `CoordinateArgs`, `AngleDefs` | here, [step-13](step-13.md), [step-14](step-14.md) |
| **Target** — what it compiles against | `RuntimeInclude`, `ExtraIncludes`, `IncludeDir`, `RegulatorTemplate`, `RegulatorAlias` | [step-15](step-15.md) |
| **Emission strategy** — how it is spelled | `CrossTraceCSE`, `FullParallel`, `Components`, `SymbolDefs`, `ShareInterpolatorIndex`, `HoistLoopConstLookups` | here, [step-20](step-20.md), [step-21](step-21.md) |
| **Build orchestration & correctness probes** | `Offline`, `RunGenerator`, `RealProbe`, `PruneRealTraces` | [step-20](step-20.md), [step-21](step-21.md) |
| **Complex flows** — how a complex integrand becomes a real kernel | `ComplexRuntimeProjection`, `ComplexEndProjection`, `RealOutput` | below |
| **Finite temperature & precision** | `MatsubaraVar`, `DecayingRegulators`, `MatsubaraFiniteExtent`, `ComputeType` | below |

The options not covered by a later step:

| Option | Default | What it does |
|---|---|---|
| `ParameterOrder` | `Automatic` | Order of the runtime parameters in the kernel signature. `Automatic` keeps scalars first, then dressings; `MakeNTKernelDiFfRG` passes DiFfRG's `Parameters` order so the signature matches the integrator's forwarded tuple. |
| `ComplexRuntimeProjection` | `False` | For a coefficient whose `I` sits where the symbolic real/imaginary split cannot reach (a finite-density denominator `l0 + I muq`): `False` aborts (`MakeNTKernel::cplxnest`); `True` keeps it factored and takes `ntRe`/`ntIm` of it at kernel runtime, with complex arithmetic in the kernel. |
| `ComplexEndProjection` | `False` | Skip the symbolic Pure/RePart projections and the imaginary-part probe; emit one body `ntRe[integrand]`, i.e. the pointwise real part of the full complex integrand. Requires `RealOutput -> True`. |
| `RealOutput` | `False` | The consumer takes a real value, so the complex kernel body is not lowered at all (it is the expensive one). If the probe then finds a surviving imaginary part, the kernel returns the RePart body — a truncation, flagged by a `#warning` in the emitted header. |
| `MatsubaraVar` | `None` | Name of the Matsubara-frequency variable of a finite-T flow. When set, the generator proves whether the kernel is even in it and, if so, emits DiFfRG's `matsubara_even` trait. |
| `DecayingRegulators` | `Automatic` | Regulator functions that decay super-polynomially in their argument (`Automatic`: DiFfRG's six regulator wrappers). A summand multiplied by one has finite frequency extent, which decides the `matsubara_finite_extent` trait. List nothing for a regulator that decays only algebraically. |
| `MatsubaraFiniteExtent` | `Automatic` | `Automatic` derives the `matsubara_finite_extent` trait from the algebra; `True`/`False` force it. Forcing `True` on a summand that does not die above the regulator's support silently truncates the Matsubara sum. |
| `ComputeType` | `"double"` | Precision of the EMITTED kernel (`"double"`, `"float"`, or a complex type of either); derivation and the generator stay in double. |

```{admonition} This list moves
:class: note
The option set is not frozen — `GlobalCollect` and `NumericContract` were options and are now
unconditional behaviour, and `ShareInterpolatorIndex` / `HoistLoopConstLookups` are recent
additions. If this table and `Options[MakeNTKernel]` in `Codegen.m` disagree, believe `Codegen.m`.
The *listings* on this page cannot go stale — they are pulled from the script that runs — but prose
about a specific option can.
```

## The commented program

`Tutorials/step-09-codegen-options/options.wls`.

### The network

```{literalinclude} ../../../Tutorials/step-09-codegen-options/options.wls
:language: mathematica
:start-after: "@snip begin: net"
:end-before: "@snip end: net"
```

````{admonition} A real `.wls` trap, and the assertion that catches it
:class: warning
Note the parentheses around the sum. In a `.wls` script **a line that is already a complete
expression ends there**, so

```mathematica
net = termA
    + termB;
```

assigns only `termA` and leaves `+termB` as a separate, discarded expression. The kernel that comes
out is perfectly valid — it compiles, it runs, it returns finite numbers — and it computes the
wrong network. Nothing warns you.

Both this script and [step-07](step-07.md) hit exactly this while being written. The fix is the
parentheses; the *guard* is asserting the expected diagram count:

```mathematica
If[nd =!= 3, Print["FAIL: expected 3 diagrams, got ", nd]; Exit[1]];
```

Assert the number you expect, not merely that it is nonzero — a truncated sum still has diagrams.
````

### The four variants

```{literalinclude} ../../../Tutorials/step-09-codegen-options/options.wls
:language: mathematica
:start-after: "@snip begin: variants"
:end-before: "@snip end: variants"
```

## Results

```bash
cmake --build build --target options && ./build/options
```

```text
one network, four emissions, at (p, l1, cos1) = (1.7, 0.9, 0.35)
  default                      =       5.608045   (closed form 5.608045)
  CrossTraceCSE -> True        =       5.608045   == default
  Constant -> myZ[p]           =       5.608045   == default
  renamed namespaces           =       5.608045   == default
the constant() entry point
  default constant(p)          =              0   (expected 0)
  Constant -> myZ[p]           =          2.445   (expected myZ(p) = 2.445)
ALL TESTS PASSED
```

All four agree, as they must. Now the interesting part.

### What the default emits

```cpp
return fma(_interp1, tr1(fenv),
       fma(tr0(fenv), cos1 * l1 * p, 0.));
```

Two trace functions. `tr1` is the collected $P^T + P^L$ pair — the two diagrams sharing the
dressing coefficient `myZ(l1)`, folded into one trace so that `_interp1` multiplies them once. `tr0`
is the third diagram, multiplied by its own scalar coefficient `cos1*l1*p`. This is the
$\sum_{\text{diagrams}}\text{coeff}\times\text{trace}$ shape, made concrete.

### `CrossTraceCSE -> True` — a genuinely different shape

This one is visible immediately. Instead of separate `trN` functions, the traces header now
contains a single fused program:

```cpp
static inline constexpr int nenv = 1;          // was 2
static inline void fill(double *f, double l1, double cos1, double p) {
  f[0] = p;
}
using trace_all_t = double;
static inline void trace_all(const double *f, double *t) {
  const double s1 = f[0];
  const double s2 = s1*s1;
  t[0] = s2;
  t[1] = s2;
}
```

and the kernel calls it once:

```cpp
trace_all_t tarr[2];
trace_all(fenv, tarr);
return fma(_interp1, tarr[1], fma(tarr[0], cos1 * l1 * p, 0.));
```

Two things happened. The CSE now runs **across** traces, so it noticed that both traces are the
same value $p^2$ and computed `s1*s1` once instead of twice. And `nenv` dropped from 2 to 1 — the
fused program needed one fewer env slot.

That looks like a pure win, and on this toy it is. It is nevertheless **off by default**, because
on a real flow the fused program is one enormous function: register pressure goes up, the compiler's
optimiser gets slower (sometimes dramatically), and GPU occupancy can fall. Measured on the flows
in this repository it is roughly neutral-to-slightly-negative (~0.87×). Treat it as something to
measure per flow, not to switch on globally — [step-20](step-20.md).

### `Constant -> expr` — the one that changes behaviour

```cpp
static inline auto constant(const double& p, const auto& myZ)
{
  return myZ(p);
}
```

`kernel()` is unchanged; only `constant()` differs. This is a **loop-independent** term: it does not
enter the integrand, it is flat-added by the consumer's integrator after quadrature. That is
exactly what an RG counterterm needs, and [step-13](step-13.md) is the story of one.

Note that `constant()` gets its own parameter list — the coordinate arguments and the dressings,
but not the integration variables. There is nothing to integrate over.

### Renamed namespaces — the retargeting hook

```cpp
namespace step09_kernels { class Opt_ns_kernel { … }; }
using step09_kernels::Opt_ns_kernel;
```

A generated kernel is **self-contained by default**: namespace `numtracer_kernels`, including only
`numtracer/codegen/runtime.hpp`. `KernelNamespace`, `SupportNamespace`, `RuntimeInclude` and
`ExtraIncludes` redirect all of that at a host framework. That is how the in-repo test kernels emit
against the test shim, and how [step-15](step-15.md) emits into DiFfRG's namespace against DiFfRG's
support API.

## Options you will meet later

* **`Decorator`** — the C++ decoration on every emitted function. `"static inline"` by default;
  `"static __host__ __device__ inline"` makes the same kernel compile for CUDA
  ([step-21](step-21.md)). Nothing else changes, which is the point of emitting branch-free,
  division-free arithmetic.
* **`AngleDefs`** — the map from the scalar half's symbolic loop–external cosines to the frame's
  angle variables ([step-07](step-07.md), used for real in [step-14](step-14.md)).
* **`Offline`** — do not compile and run the generator now; write a manifest and let the consumer's
  CMake build do it. Essential for expensive flows ([step-20](step-20.md)).
* **`RealProbe` / `PruneRealTraces`** — how a complex-valued flow decides which of three emitted
  kernel bodies is valid ([step-21](step-21.md)). `PruneRealTraces` is documented as hazardous in
  combination with the probe; do not reach for it casually.
* **`DressingType` / `ScalarParams` / `ADParams`** — the C++ types of the kernel's parameters, and
  which of them participate in automatic differentiation ([step-15](step-15.md)).

## Possibilities for extensions

1. **Reproduce the truncation bug on purpose.** Remove the parentheses around `net` and re-run. The
   assertion fires with "expected 3 diagrams, got 1". Then comment the assertion out, rebuild, and
   watch the checker fail against the closed form instead. Compare how much easier the first
   failure was to diagnose than the second.

2. **See the dressing-coefficient grouping do more work.** Add several more diagrams sharing the
   same dressing product, and count the `_interp` hoists in the emission before and after.

3. **Measure `CrossTraceCSE` honestly.** Time both variants in a tight loop over $10^7$ points. On
   this toy the difference is noise; the exercise is to build the harness, because
   [step-20](step-20.md) asks you to do it for real.

4. **A `Decorator` that does not compile.** Set `"Decorator" -> "constexpr"` and see what the
   emitter does with it. The decorator is pasted through verbatim — it is a string, not a validated
   enum — which is the source of its flexibility and of some confusing errors.

5. **Retarget for real.** Write a tiny `my_runtime.hpp` providing `powr`, `pow`, `sqrt` and `fma`
   in namespace `mine`, then emit with `"RuntimeInclude" -> "my_runtime.hpp"` and
   `"SupportNamespace" -> "mine"`. If it compiles and agrees, you have done in miniature what
   [step-15](step-15.md) does for DiFfRG.

## The plain program

```{literalinclude} ../../../Tutorials/step-09-codegen-options/options.wls
:language: mathematica
```

```{literalinclude} ../../../Tutorials/step-09-codegen-options/check_options.cpp
:language: cpp
```

That completes the front-end group. [step-10](step-10.md) begins the physics workflows: where the
network comes from in the first place.
