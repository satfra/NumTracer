# Front-end & codegen

The Wolfram layer (`mathematica/`) drives NumTracer from a symbolic front end: it takes a traced
tensor network, **contracts it numerically at build time**, and emits a flat C++ kernel. It is
the automation of the whole path — analyse the network, fix the frame, fold colour, trace the
Dirac chain, reduce the Lorentz network, lower to straight-line arithmetic — with the
scalar (dressing) half delegated to FunKit's mature COEN emitter.

A generated kernel has the shape every flow shares,

$$\text{kernel} = \sum_{\text{diagrams}} \text{coeff}(\text{dressings}\dots)\times\text{trace}(\text{momenta}),$$

where the *trace* is the contraction of projectors, metrics, momentum vectors, gammas and
structure constants, and the *scalar coefficient* (dressings, regulator, propagator
denominators) is flat C++. NumTracer owns only the tensor half; everything scalar (CSE,
`powr<n>`, function/class/header boilerplate, formatting) goes through FunKit's COEN
(`CppForm` / `MakeCppFunction` / `MakeCppClass` / `MakeCppHeader`).

## The two entry points

`NumTrace[net, "Frame" -> frame, "Args" -> {…}]` analyses a tensor network written in the DSL
heads below into an `NTKernel` — a list of diagrams, each a scalar coefficient times one or more
contraction components, plus the env-id layout and the loop frame that fixes each momentum's
four components.

`MakeNTKernel[ntk, genFile, kernelFile, tracesFile, "Backend" -> "Numeric", …]` serialises that
into a kernel. The numeric backend emits a small C++ **generator** to `genFile`, then *compiles
and runs* it (`RunProcess`) to produce the committed straight-line traces header `tracesFile`
(a set of `trN(const double* f)` functions) and the kernel `kernelFile`, which fills the few
frame symbols and calls `trN(f)`. The colour and dressing halves are emitted by FunKit COEN.

The one-file form `MakeNTKernel[ntk, file, …]` defaults to `"Backend" -> "Dense"`: the
brute-force `numtracer::dense::DTensor` mirror, read via `.scalar_value().re`. It carries no
optimisation and is the **numeric oracle** the generated kernel is cross-checked against — a
driver can emit both from one `ntk`.

```{note}
The earlier backends are retired and now error: `"Backend" -> "ET"` (the symbolic
expression-tensor path) and `"Backend" -> "Generate"` (the symbolic invariant-basis generator)
both point at `"Numeric"`, the single canonical generation path.
```

## What the generator does

`MakeNTKernel "Numeric"` builds, per diagram, a Dirac chain and a pure-Lorentz network from the
DSL heads, then calls the [numeric contraction engine](numeric-engine.md): the Dirac trace folds
by 4×4 chiral matrix products, the Lorentz network reduces by bounded index elimination, and the
colour factor folds to a number. Each diagram becomes one small polynomial (`MPoly`) in the
frame's scalar symbols, which is lowered (CSE + Horner) into the shared `f[]` symbol layout and
printed.

The kinematic data the generator needs — each momentum's components in the frame, the propagator
denominators, and the symbol layout — is **derived from the frame**, not hand-written. Because
the contraction is numeric over a fixed frame, there is no symbolic blow-up to manage: the
generator is serial, runs in seconds, and is negligible in RAM. To keep the generator itself
cheap to compile, its per-network builders are split across several translation units compiled
`-O0` in parallel.

## The DSL

Open-index tensor heads, native to NumTracer but mirroring the FunKit/FormTracer vocabulary:

| head | meaning |
|---|---|
| `ntMetric[mu, nu]` | Lorentz metric $\delta_{\mu\nu}$ |
| `ntVec[q, mu]` | momentum component $q_\mu$ |
| `ntTransProj[q, mu, nu]` | transverse projector $P_{\mu\nu}(q)$ |
| `ntLongProj[q, mu, nu]` | longitudinal projector |
| `ntGamma[mu, …]`, slashed momenta | Dirac gamma matrices and $\slashed q = q^\mu\gamma_\mu$ |
| `ntSUNf[N, a, b, c]` | SU(N) structure constant $f^{abc}$ (rank $N$) |
| `ntSUNT[N, a, i, j]` | SU(N) fundamental generator $(T^a)_{ij}$ |
| `ntSUNDeltaAdj[N, a, b]` / `ntSUNDeltaFund[N, i, j]` | SU(N) adjoint / fundamental $\delta$ |
| `ntSUNDiagFund[N, i, j, Z, s]` / `ntSUNDiagAdj[N, a, b, Z, s]` | group-diagonal $\delta$ carrying a per-component dressing $Z[a](s)$ |
| `ntSP[q1, q2]` / `ntDress[h, args]` | scalar coefficients (dot product, opaque dressing) |

The SU(N) heads are one $N$-parameterized family (rank $N$ first), so colour SU($N_c$) and
flavour SU($N_f$) coexist in one network. The group-diagonal heads `ntSUNDiagFund` /
`ntSUNDiagAdj` fold (via `sun_value_dressed`) to $\sum_a c_a Z_a(s)$ over a runtime
per-component dressing array — dressing each colour/flavour component differently *without*
splitting into one diagram per component; the Dirac trace is still computed once. See the
[dressed-flavour tutorial](../tutorials/dressed-flavour.md).

A contraction is a product of heads with repeated index labels summed (Einstein); the top-level
sum is the linear sum of diagrams.

## Eager summation

The reason the generated kernels are competitive is that `NumTrace` does **not** expand a
diagram's product-of-sums into monomials. Each vertex is a sum of tensor structures; `NumTrace`
keeps it as a single summed object and the diagram as a single contraction. Using
$(A_1+A_2)\cdot B = A_1 B + A_2 B$, the cost is bounded by the contraction path, independent of
structures-per-vertex and vertices — the $t^V$ monomial blow-up never forms.

One subtlety is essential: a sum that **bridges colour and Lorentz** — the four-gluon vertex,
whose colour structure-sum is correlated with its Lorentz structure-sum — is **distributed**
(only those few colour structures), so colour still folds to a number and the Lorentz network
stays small; single-sector vertex sums (the three-gluon Lorentz vertex) stay eager.

## Dressed numerators

A dressed propagator numerator — e.g. $M_q\,\delta + Z(p)\,\slashed p$ with the *same* open
spinor indices — is another sum that must not be distributed: with $D$ such factors it would
explode into $2^D$ diagrams. The `"DressingCollection"` option (on `NumTrace` and `FromFunKit`,
default `True`) keeps it eager as a single collected trace (an `ntDressedNum` token carrying the
dressing coefficients symbolically), so the Dirac trace is computed once and the dressing leaves
appear as runtime atoms. This is what lets a flow dress individual flavours or colour/flavour
components differently in one trace. Flows with no dressed-numerator sum regenerate
byte-identical with the option on or off.

## Disconnected diagrams

A single diagram can disconnect into several index-disjoint closed pieces (a quark loop *and* a
flavour trace, or two independent traces). Their scalar values **multiply**. When a diagram has
two or more non-constant components, the generator emits one as the additive base and each other
as its own fused trace group, and the assembly forms the product $\prod_k \text{trace}_k$ — each
$P$ computed once, no trace-polynomial blow-up. This mirrors the dense path's
$\text{coeff}\times\prod \text{toks}$. (`discdirac_num`, `flow_zacbc_num` validate it; flows
without a disconnected diagram are unaffected.)

## The FunKit adapter

`FromFunKit[flow]` (`FunKitAdapter.m`) rewrites a FunKit trace — after `// dressingRules //
PropParam` — into the DSL: a pure rewrite over the closed FunKit head vocabulary
(`transProj → ntTransProj`, `FCol → ntSUNf[Nc, …]`, `deltaLorentz → ntMetric`, …; the SU(N)
group tokens get their rank N injected — colour from `Nc`, the isospin group from the
`"FlavourGroup"` option), so the front end
consumes a flow with no dependence on the TensorBases internals. The drivers in `tests/gen/`
load the field-space setup, run `FromFunKit → NumTrace → MakeNTKernel`, and write the committed
`tests/gen/*.hh` kernels.

```{note}
Generation runs under `wolframscript` and reuses FunKit's COEN, so it needs a Wolfram kernel and
FunKit (and, for the YM/QCD flows, a FORM executable — FunKit uses FORM to compute the tensor
*bases*, not to do NumTracer's trace). The generated kernels are committed under `tests/gen/`, so
the ordinary C++ build and tests stay dependency-free; regenerate only when the flows or the
codegen change.
```

## Self-contained kernels

Generated kernels are self-contained by default: they `#include
"numtracer/codegen/runtime.hpp"` (which supplies `numtracer::complex` and
`numtracer::compute::{powr,pow,sqrt,fma}`), are wrapped in a neutral namespace, and take generic
dressing parameters — no consumer dependency. A consumer that provides its own support API points
the codegen at it via the `"RuntimeInclude"` / `"SupportNamespace"` / `"KernelNamespace"` /
`"DressingType"` / `"ExtraIncludes"` options; the in-repo tests use these to emit against the
test shim. A `"Decorator"` option annotates the kernel for GPU use (`__host__ __device__`), and
`"AngleDefs"` hoists shared kinematic angle definitions into the kernel body so a common angular
`sqrt` is computed once.

A few further `MakeNTKernel` options shape the emitted kernel: `"GlobalCollect"` (default
`True`) folds fundamental colour numerically so diagrams sharing a dressing coefficient merge
into one trace; `"ScalarParams"` threads extra loop-independent doubles into the signature;
`"Constant"` flat-adds a loop-independent term (DiFfRG's `constant(...)`); and `"Components"`
restricts which external-projection components are emitted.
