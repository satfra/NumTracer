# Generating kernels from Mathematica

The previous tutorials walked through each step of the running example by hand. The Mathematica
front-end **writes all of it for you** from a tensor network in a small DSL: it analyses the
network, fixes the frame, contracts each diagram numerically, lowers the result, and emits a
flat C++ kernel.

## Why generate?

Everything the hand path required — choosing a frame and env layout, tracing the Dirac chain,
contracting the projectors, lowering, transcribing the dressing scaffolding — is mechanical and
*derivable* from the flow equation, which in a DiFfRG/FunKit workflow already exists as a symbolic
object. The codegen closes that gap:

$$
\text{kernel} \;=\; \sum_{\text{diagrams}} \text{coeff}(\text{dressings}\dots)
\times \text{trace}(\text{momenta}),
$$

`NumTrace` analyses a DSL network into that shape, and `MakeNTKernel` contracts it and emits it as
a committed C++ header.

```{note}
**Requirements.** The codegen runs under `wolframscript` and reuses FunKit's COEN emitter for the
scalar half, so it needs a Wolfram kernel and FunKit (and, for the YM/QCD flows, a FORM executable
— FunKit uses FORM to compute the tensor *bases*, not to do NumTracer's trace). The generated
`tests/gen/*.hh` are committed, so the ordinary C++ build and `ctest` never need any of this —
regenerate only when a flow or the codegen changes.
```

## A minimal driver

The smallest end-to-end example generates a kernel for the warm-up network
$p\cdot P(l)\cdot p / p^2 = 1 - \cos^2\theta$:

```mathematica
Needs["FunKit`"];                                    (* COEN C++ emitter *)
Needs["NumTracer`"];                                 (* loads DSL.m / Frames.m / Codegen.m *)

SetNc[3];                                            (* NumTrace needs Nc; this net is pure-Lorentz *)
{p, l1, cos1} = {Global`p, Global`l1, Global`cos1}; (* runtime scalar kernel args *)
{qp, ql}      = {Global`qp, Global`ql};             (* momentum tags: external, loop *)

frame = NumTracer`propFrame[p, l1, cos1, qp, ql];          (* symmetric-point, 1 angle *)
net   = (1/ntSP[qp, qp]) ntVec[qp, mu] ntTransProj[ql, mu, nu] ntVec[qp, nu];
ntk   = NumTrace[net, "Frame" -> frame, "Args" -> {l1, cos1, p}];

MakeNTKernel[ntk, "selftest_kernel.hh", "Name" -> "selftest_kernel", "Dressings" -> {}];
```

Reading the network `(1/ntSP[qp, qp]) ntVec[qp, mu] ntTransProj[ql, mu, nu] ntVec[qp, nu]`:
repeated index labels (`mu`, `nu`) are summed (Einstein), `ntSP[qp, qp]` is the scalar product
$qp\cdot qp = p^2$, and the top-level sum would be the sum of diagrams. The full head vocabulary
and the entry points are on the [codegen page](../internals/codegen.md).

- **`NumTrace[net, "Frame" -> …, "Args" -> …]`** analyses the network into an `NTKernel`: the
  list of diagrams, the env-id layout, and the loop frame.
- **`MakeNTKernel[ntk, file, …]`** serialises it to a C++ header.

## The generated kernel

The emitted header is a flat, templated kernel — the same shape as a FORM reference kernel: a
namespace, a class templated on the regulator, and a `kernel(args…)` that fills the frame symbols
once and returns the value. The kernel arguments are exactly the `"Args"` list, in order; a flow
with dressings takes them as extra arguments after the scalars. By default the header is
self-contained (it includes only `numtracer/codegen/runtime.hpp`); the in-repo flows emit against
the test shim so they can be compared to the FormTracer oracles.

The committed numeric kernels live under `numtracer/tests/gen/` — for example
`Zq_num_kernel.hh`, the quark self-energy generated below. Inspecting one shows the
`fill()` + `trN(const double* f)` structure of a [generated kernel](../internals/codegen.md).

## Generating the quark self-energy

The full running example is not typed in by hand: FunKit derives the quark two-point flow from
the Wetterich equation, and `FromFunKit` rewrites the traced loop — gammas and all — into the DSL:

```mathematica
Needs["FunKit`"]; Get["DiFfRG`"];
Get["tests/qcd_setup.m"]; Needs["NumTracer`"];

frame = NumTracer`propFrame[p, l1, cos1, p1, l1];   (* 1 angle: p1 along axis 0, l1 at cos1 *)
args  = {l1, cos1, p, k};
dress = {ZA, Zq, Mq, ZAqbq1, dtZA, dtZq};

fRGqbq = FTakeDerivatives[WetterichEquation, {qb[i1], q[i2]}] // FTruncate // FRoute;
trace  = FTerm[-TBGetProjector["qbq", 1, …]] ** (fRGqbq["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]);
net    = FromFunKit[PropParam[trace // dressingRules]];   (* FunKit flow -> NumTracer DSL *)

ntk = NumTrace[net, "Frame" -> frame, "Args" -> args];
MakeNTKernel[ntk, "gen_zq.cpp", "Zq_num_kernel.hh", "Zq_num_nets.hh",
             "Name" -> "Zq_num_kernel", "Dressings" -> dress];          (* numeric (default, 3-file) *)
```

`"Dressings" -> dress` declares the interpolator arguments; the emitted kernel takes them after
the `"Args"` scalars, exactly like a FormTracer kernel.

## The generation backend

`MakeNTKernel` uses the **`"Numeric"`** matrix-product engine (the three-file form above, the
default): it emits a small generator, compiles and runs it, and writes the committed straight-line
kernel plus its trace functions. This is the canonical — and only — generation path.

```{note}
The earlier symbolic backends are retired: `"Backend" -> "ET"` (the expression-tensor path) and
`"Backend" -> "Generate"` (the symbolic invariant-basis generator) now error and point at
`"Numeric"`.
```

## Verify

```bash
# regenerate the committed Zq kernel (needs wolframscript + FunKit + FORM):
cd numtracer && wolfram -script tests/gen/gen_zq_numeric.wls
# validate the committed kernels (no Wolfram needed):
cmake --build build -j4 && ctest -R "numeric|flow" --output-on-failure
```

Each generated kernel is checked against a FORM or dense numeric oracle over random points; the
two agree to many significant figures, which pins down the engine's algebra.
