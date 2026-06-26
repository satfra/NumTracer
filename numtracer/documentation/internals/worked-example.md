# Worked example: the quark self-energy

This page follows one real calculation through the whole pipeline: the one-loop quark
self-energy contribution to the flow of the quark propagator in QCD (call it **FlowZq**). It is
the example the [tutorials](../tutorials/index.md) build piece by piece; here we look at it as a
single end-to-end trace.

## What is being computed

Deriving the flow of the quark two-point function produces **two diagrams** — the same
quark–gluon triangle, once with the regulator insertion on the **gluon** line and once on the
**quark** line. After the colour, Dirac and Lorentz algebra is done, the integrand at one
phase-space point factorises into

$$
\text{integrand} \;=\; \underbrace{C_F}_{\text{colour}} \times
\Big(\,\text{term}_1\big|_{q=l-p} + \text{term}_2\big|_{q=l+p}\Big),
$$

where $C_F = (N^2-1)/2N = 4/3$ for SU(3) is the common colour factor and $q = l \mp p$ is the
internal quark momentum in each diagram. Each term carries the Dirac-trace numerator

$$
T_{\text{num}} = \mathrm{tr}\!\big[\slashed p\,\gamma^\mu\,\slashed q\,\gamma^\nu\big]\,
P_{\mu\nu}(l),
\qquad P_{\mu\nu}(l) = \delta_{\mu\nu} - \frac{l_\mu l_\nu}{l^2},
$$

and a scalar structure: the dressings ($Z_A$, $Z_q$, the quark–gluon vertex), the bosonic and
fermionic regulators, and the dressed quark denominator with a live mass term.

**The colour factor and the Dirac trace are produced by NumTracer**, not hardcoded. The
remaining scalar structure is ordinary C++.

## How the pipeline produces it

The flow is handed to the front end as a tensor network (typically imported from a FunKit
derivation with `FromFunKit`). `NumTrace` analyses it into the two diagrams, attaches the
one-angle loop frame, and lays out the env. Each diagram is then contracted by the
[numeric engine](numeric-engine.md):

* **Colour.** $T^a_{ij}T^a_{jk} = C_F\,\delta_{ik}$ has no momentum content, so the colour
  network folds to the single number $C_F = 4/3$ at build time (`network/sun_net.hpp`). Both
  diagrams share it, so it factors out of the sum.

* **Dirac trace.** $\mathrm{tr}[\slashed p\,\gamma^\mu\slashed q\,\gamma^\nu]$ has two free
  Lorentz legs ($\mu, \nu$). The engine traces the closed spinor chain by 4×4 chiral matrix
  products and reconstructs it as a compact Lorentz tensor — no $(2n-1)!!$ Wick expansion.

* **Lorentz.** That Dirac tensor is contracted against the transverse projector $P_{\mu\nu}(l)$
  by bounded index elimination, leaving one polynomial in the frame's scalar symbols
  ($l\!\cdot\!p$, $l^2$, …) with the $1/l^2$ propagator atom tracked per monomial.

The two diagram terms are kept as an **eager sum**, not multiplied out — exactly how the
reference kernel is structured, which keeps the polynomial small. The same compile-time trace
serves both diagrams; $q = l \mp p$ is just a different frame fill.

## From polynomial to kernel

Each diagram's polynomial is lowered by [CSE and real lowering](cse-and-lowering.md) into the
shared symbol environment, and the codegen prints a flat kernel: the trace functions over the
frame symbols, a `fill` that evaluates the symbols from the runtime arguments
$(|l|, \cos\theta, |p|, …)$ once per call, and the assembly
$\sum_{\text{diagrams}} \text{coeff}\times\text{trace}$ with the scalar coefficients (dressings,
regulators, denominators) emitted by FunKit's COEN. The committed kernel is plain C++ depending
only on `codegen/runtime.hpp`.

## How it is validated

The same network is also emitted through the `"Dense"` backend — the brute-force
`dense::DTensor` numeric mirror, with no optimisation — and the two kernels are evaluated against
a copied FormTracer oracle over a large set of random phase-space points; they agree to many
significant figures. Because the scalar scaffolding is shared, any disagreement could only come
from the parts NumTracer computes independently — the colour contraction and the Dirac trace — so
the match pins down the engine's algebra. The colour fold is additionally checked against
$C_F = 4/3$, and the bare trace against its analytic scalar-product polynomial.

## The whole pipeline, in one place

This example touches every layer:

1. **[Front-end & codegen](codegen.md)** — `NumTrace` analyses the network; `MakeNTKernel`
   contracts and emits the kernel.
2. **[The numeric engine](numeric-engine.md)** — the colour fold, the matrix-product Dirac
   trace, and the Lorentz reduction.
3. **[Sector data](sectors.md)** — the typed-out gamma and SU(N) tables the contraction reads.
4. **[The expression algebra](expression-algebra.md)** and **[CSE + lowering](cse-and-lowering.md)**
   — turn each diagram's polynomial into straight-line real code.

That is the complete picture: describe the network, contract and fold it numerically at build
time, and lower the result to straight-line real arithmetic — reproducing a FORM-generated kernel
without a symbolic-algebra runtime.

For per-symbol API details, see the [C++ API reference](../doxygen/NumTracer/html/index).
