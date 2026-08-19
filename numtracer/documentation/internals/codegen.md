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

`MakeNTKernel[ntk, genFile, kernelFile, tracesFile, …]` serialises that into a kernel. It emits a
small C++ **generator** to `genFile`, then *compiles and runs* it (`RunProcess`) to produce the
committed straight-line traces header `tracesFile` (a set of `trN(const double* f)` functions) and
the kernel `kernelFile`, which fills the few frame symbols and calls `trN(f)`. The colour and
dressing halves are emitted by FunKit COEN.

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
| `ntSUNDiagFund[N, i, j, spec]` / `ntSUNDiagAdj[N, a, b, spec]` | group-diagonal $\delta$ dressing selected components via a rules-list `spec` |
| `ntSP[q1, q2]` / `ntDress[h, args]` | scalar coefficients (dot product, opaque dressing) |

The SU(N) heads are one $N$-parameterized family (rank $N$ first), so colour SU($N_c$) and
flavour SU($N_f$) coexist in one network. The group-diagonal heads `ntSUNDiagFund` /
`ntSUNDiagAdj` take a `spec` rules list `{c -> name, …, Default -> defName}` (1-based component
indices → distinctly-named scalar dressings; unnamed components collapse to `Default` or drop)
and fold (via `sun_value_dressed`) to $\sum_a c_a Z_a(s)$ over those named runtime scalar
dressings — dressing selected colour/flavour components differently *without* splitting into one
diagram per component; the Dirac trace is still computed once. See the
[step-17](../tutorials/step-17.md).

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

### Collected vertices (open Lorentz legs)

The same argument applies to an internal **vertex** sum — a quark–gluon vertex written over a
tensor basis, $\sum_a c_a T_a$ — except that its structures carry *open gluon legs* that the
surrounding net must close. The `ntDiracSlot` token generalises `ntDressedNum` to exactly that: a
coefficient-weighted sum of structures sharing both the spinor in/out pair **and** the same set of
$k \ge 1$ open Lorentz legs. ($k=0$ is the propagator numerator above; a two-gluon vertex is
$k=2$, and so on.) Each option is split into a Dirac-token chain plus any Lorentz-net factors —
see the [step-18](../tutorials/step-18.md) for the data model.

**This path is OFF by default and must be enabled per flow.** Collection replaces $3^n$ diagrams
with one, which is a large win where the combination count stays bounded — `za3_147` generates in
8.9 s collected against 13.3 s distributed — but the codegen expansion materialises the
structure×dressing Cartesian product, and on a high-multiplicity flow (full-basis `ZAAqbq`: order
$10^6$ combinations per net) that OOMs the Wolfram kernel. Such flows *must* distribute. Enable it
with `NT_VERTEX_COLLECT=1`, or by setting `NumTracer`Private`$ntVertexCollect = True` in the
generator script, as `tests/gen/gen_qcd_za3_147_numeric.wls` does. The propagator-numerator
collection ($k=0$, at most two options apiece) is bounded and stays on unconditionally.

## Disconnected diagrams

A single diagram can disconnect into several index-disjoint closed pieces (a quark loop *and* a
flavour trace, or two independent traces). Their scalar values **multiply**. When a diagram has
two or more non-constant components, the generator emits one as the additive base and each other
as its own fused trace group, and the assembly forms the product $\prod_k \text{trace}_k$ — each
$P$ computed once, no trace-polynomial blow-up, matching the factored form
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
Generation runs under `wolframscript` (a Wolfram kernel) and reuses FunKit's COEN for the scalar
half. The *derivation* of the flow is separate and optional: `FromFunKit` imports one, but a
hand-built DSL network needs no FunKit derivation. **FORM is not part of the trace or the general
path** — it enters only when regenerating the committed YM/QCD reference-test flows, where FunKit
uses it to compute the tensor *bases*. The generated kernels are committed under `tests/gen/`, so
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

### Regulators are the consumer's

The emitted kernel is a **plain class** — no regulator template parameter, no regulator definitions.
NumTracer emits the kernel and nothing else. A flow whose dressing rules mention the regulators (the
inverse propagators typically do) emits **unqualified** calls to `RB`, `RF`, `RBdot`, `RFdot`,
`dq2RB`, `dq2RF`, exactly as it emits unqualified calls to the dressing parameters and to `powr`.
Supplying them is the consumer's job: put them in a header and pull it in with

```wolfram
"ExtraIncludes" -> {"my_regulators.hpp"}
```

(or define them in the `"RuntimeInclude"` header). Definitions at global scope are found from a
kernel emitted into any `"KernelNamespace"`, since unqualified lookup runs class → namespace →
global. The in-repo tests do exactly this — see `tests/refshim/nt_regulators.hpp`, which
`tests/refshim/shim.hpp` includes so every shim-emitted kernel gets it for free.

`"RegulatorTemplate" -> True` restores the fRG/DiFfRG shape instead: `template<typename REG> class
…` plus private wrappers forwarding to `REG::RB` and friends. `MakeNTKernelDiFfRG` sets it (with
`"RegulatorAlias" -> True`, which emits `using Regulator = REG;` and implies it), because DiFfRG's
scaffold forward-declares the kernel as a template and instantiates it as `KERNEL<Regulator>`.

Folding fundamental colour numerically — so that diagrams sharing a dressing coefficient merge into
one trace — is unconditional; there is no option for it.

A few further `MakeNTKernel` options shape the emitted kernel: `"ScalarParams"` threads extra
loop-independent doubles into the signature; `"Constant"` flat-adds a loop-independent term
(DiFfRG's `constant(...)`); and `"Components"` restricts which external-projection components are
emitted.

## Generation knobs (environment variables)

Generation is a two-phase C++ program the Wolfram front end compiles and runs (phase A contracts
each distinct trace; phase B folds each net and lowers it to straight-line SSA). Both phases, and
a few front-end steps, are steered by environment variables rather than options, because they are
resource dials rather than semantics — set them around one flow, not globally.

**How a boolean is read.** ONE rule, everywhere, on both sides of the pipeline: a flag is ON when
the variable is **set, non-empty, and not the single character `0`**. `FOO=1`, `FOO=true`, `FOO=on`
and `FOO=yes` all enable; `FOO=0`, `FOO=` and unset all disable. C++ reads it through
`env_flag` (`core/envvar.hpp`), Wolfram through `ntEnvFlag` (`mathematica/DSL.m`). This is worth
stating because three other spellings used to coexist, under two of which `FOO=0` turned the flag
**on** — which is how `NT_NO_LABEL_CHECK=0` once disabled a correctness guard. Numeric knobs are
different: for them `0` is a value, so *empty or absent* is what means "unset".

**When to set them.** The Wolfram-side flags are read lazily, so `SetEnvironment[...]` works from a
`.wls` that has already loaded the package. (They used to latch at load time, which made the
recipe below silently ineffective.)

| variable | affects | default | effect |
|---|---|---|---|
| `NT_GEN_MAXW` | phase A | hardware concurrency | Cap phase-A workers. Peak RSS is ~`W` × the per-contraction working set, so this is the hard RAM ceiling for the contraction. Only ever *lowers* the count. |
| `NT_GEN_MAXW_B` | phase B | `NT_GEN_MAXW` | Cap phase-B workers separately. Each holds one full contraction, so it is the RAM dial for dense quark flows — but the lowering is compute-heavy and parallel, so a low value badly slows big-kernel flows (`ZA4_147`: ~600 s at `nB=2` vs tens of seconds unthrottled). Throttle around the RAM-heavy flows only. |
| | | | *Both are also settable from Wolfram via `SetNumTracerThreads[nA, nB]` / `GetNumTracerThreads[]`, which export them for the generator's child process.* |
| `NT_GEN_GROUP_WINDOW` | phase B | `max(64, 8·W)` | How many net polynomials may be in flight at once — the streaming-fold RAM bound. A value ≥ the net count reproduces the old all-resident schedule. |
| `NT_GEN_MEMO_MAX` | phase A/B | all distinct traces | Cap how many contracted traces stay resident; the rest are recomputed in phase B. The memory-bound dressed flows (`ZAAqbq`) dial this back. |
| `NT_GEN_NO_DEDUP` | codegen | off | Disable the global sub-term dedup. The A/B control for the dedup speedup — `gen_zaaqbq1_small_numeric.wls` emits a dedup-off kernel with it as a graded reference. |
| `NT_GEN_NO_POLYDIV` | phase A | off | Disable multi-term denominator cancellation (`divThroughPolyAtoms`). Costs a large amount of kernel size on quark/ghost loops; the control for that measurement. |
| `NT_GEN_CC_CHUNK` | codegen | full fusion | Chunk size for fused-trace (`CrossTraceCSE`) lowering. Measured: full fusion wins. |
| `NT_VERTEX_COLLECT` | front end | off | Enable the open-leg vertex collection (`ntDiracSlot`, above). Off because high-multiplicity flows OOM; on where it wins. |
| `NT_NO_LABEL_CHECK` | front end | check on | Skip the per-diagram label census. ~14% of front-end time, but it is the guard that catches a label occurring more than twice — which otherwise becomes a silently wrong contraction. |
| `NT_GEN_NOINLINE_MIN` | emission | `500` | Per-function instruction threshold above which a device trace function is emitted out-of-line. See [tests/gpu/README.md](../../tests/gpu/README.md). |
| `NT_GEN_NOINLINE_TRACES` | emission | off | Force out-of-line for every trace function, host and device — the nvcc compile-cost lever. |
| `NT_GEN_PROFILE` | both | off | Per-phase timing/RSS diagnostics from the generator binary. `=2` adds the per-wave RSS trace. Note it must be run against the compiled `gen_<name>` binary directly — the output is lost through `wolframscript`'s `Run[]`. |
| `NT_GEN_VERBOSE` | front end | off | Enable the `[prof]`/`[cse]`/`[time]` Wolfram-side diagnostics (`ntLog`). `tests/gen/regen_check.sh` sets it because its density guard greps the `[cse]` line. |
| `NT_NO_SIGN_CANON` | front end | canon on | Stop ±q sharing one env base/inverse slot. A/B control; `gen_lambda3d_small_numeric.wls` builds its graded control kernel with it. |
| `NT_NO_UNIT_GROUPS` | front end | on where the frame qualifies | Force the general `polyFrameSpec` path instead of the unit-loop rewrite (measured 70× on `lambda3d`). Same fixture, same purpose. |
| `NT_NO_SIGMA_FOLD` | front end | fold on | Distribute a recognised γ-commutator into two Dirac traces instead of folding it to one `ntSigma` token. The baseline for measuring the sub-term reduction. |
| `NT_NO_INTERP_SHARE` | emission | share on | Disable interpolator index sharing (`"ShareInterpolatorIndex"`). A/B control; results are bit-identical either way. |
| `NT_GEN_NO_KHOIST` | emission | hoist on | Disable launch-constant dressing-lookup hoisting (`"HoistLoopConstLookups"`). A/B control — **not** bit-identical (host libm vs device libdevice, last ulp). |
| `NT_GEN_HORNER_ORDERS` | lowering | 8, auto-reduced by size | Force N trial Horner orderings instead of the size-scaled default (3 above 500 monomials, 1 above 2000). |
| `NT_GEN_SNAP_DIGITS` | lowering | 14 | Decimal digits emitted coefficients are snapped to. `0` disables snapping; outside `[0,17]` falls back to the default. |
| `NT_GEN_POLYSTATS` | lowering | off | `=1` monomial/SSA counts to stderr; `=2` dumps every monomial key instead (pipe through `sort -u` to count distinct monomials across traces). |
| `NT_GEN_DEVICE` | emission | sniffed from the decorator | Declare the emitted kernel DEVICE code, which is what enables the size-gated `__noinline__`. Set by `Codegen.m` online and by the manifest offline; the decorator sniff is a back-compat fallback that no production flow reaches. |

### Deployment configuration

Not resource dials — these say *where things are* and *how generation is driven*. They were
undocumented for a long time, which is the only reason they are listed separately.

| variable | default | effect |
|---|---|---|
| `NT_GEN_CXX` | a `CCompilerDriver`-known clang++, else g++ | Compiler for the generator TUs. The only override. |
| `NT_GEN_LIB` | searched | Full path to a specific `libNumTracer.a`. |
| `NT_ALLOW_STALE_LIB` | guard on | Proceed even when that archive is OLDER than the headers it will compile against. The guard exists because the mismatch is an ODR/ABI one that silently produces wrong traces while leaving `kernel.hh` byte-identical; overriding it is almost never right. |
| `NUMTRACER_INCLUDE_DIR` | in-tree, then the installed record, then `~/.local/share/NumTracer/include` | Header root. The only name here without the `NT_` prefix. |
| `NT_OFFLINE` | the `"Offline"` option | Tri-state override: set forces offline; `0`/`false`/`no`/`off` forces **online**, overriding `"Offline" -> True`. |
| `NT_GEN_MAIN_OPT` | `-O1` | Optimisation level for the generator's main TU. `-O0` is a large win on small flows and a large loss on dense ones. |
| `NT_GEN_JOBS` | `max(2, min(24, cores, freeGB/1.5))` | Pin the parallel-compile job count, bypassing both the core and the memory bound. |
| `NT_FUNKIT_BACKEND` | FunKit's own choice | `Mathematica` or `Cpp`; pins FunKit's backend for bisecting against frozen oracles (`tests/backend_pin.m`). |
| `NT_KEEP_TRACE_CACHE` | cache wiped | Reuse the FORM trace cache instead of deleting it (the `*_form.wls` oracles; a cold run is ~100 min). |
| `NT_TEST_NDIAG` | per fixture | Cap the diagram count so a fixture stays affordable in ctest; `0` keeps all. |
| `NT_DENSE_NP` | per bench | Point count for the dense benchmarks. |
| `NT_RUN_CODEGEN` | unset ⇒ skip | Opt in to the `codegen_regen` ctest gate (it exits 77 = Skipped without it). |
| `NT_BUILD_DIR` | `numtracer/build` | Build directory `regen_check.sh` rebuilds into; set by CTest. |

### Removed

These were escape hatches with no reference anywhere — no test, no fixture, no CMake, no doc — each
selecting a path that was measured *bit-identical but slower*, or simply worse. They and the code
they selected are gone; re-running such an A/B now means deleting the branch, not setting a variable.
`NT_DIRAC_FLAT`, `NT_RANK1_PROJM`, `NT_NO_RANK1_PROJE`, `NT_GEN_NO_FMA`, `NT_GEN_NO_CONST_INLINE`,
`NT_STATS_CSV`, `NT_NO_SLOT_COLLECT_NUMERIC`, `NT_GEN_NO_TABLE_DEDUP`, `NT_GEN_NO_PCH`,
`NT_GEN_NO_COMPILE_CACHE`, `NT_GEN_NO_TCMALLOC`, `NT_GEN_JOB_MEM`, `NT_TEST_SKIP_CONTROL`.
