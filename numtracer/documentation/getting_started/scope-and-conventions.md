# Scope & conventions

NumTracer is a *general* tensor-tracing engine, but "general" has an exact meaning: it contracts
networks over a fixed set of **sectors**, in a fixed set of **conventions**. This page states them
up front so you can tell, before writing anything, whether your problem fits and how your objects
map onto the engine's. Everything here is a hard convention baked into the typed-out tables (see
[Sector data](../internals/sectors.md)), not a runtime option.

## What the engine contracts

| Sector | Objects | Index extent |
|---|---|---|
| **Lorentz** (spacetime) | metric $\delta_{\mu\nu}$, vectors $q_\mu$, transverse/longitudinal/electric/magnetic projectors, Levi-Civita $\varepsilon$ | **4** ($\mu = 0..3$) |
| **Dirac** (spin-½) | gamma matrices $\gamma^\mu$, slashed vectors $\slashed q$, $\gamma_5$, commutators $\sigma^{\mu\nu}$ | **4** (spinor) |
| **SU($N$)** (any group: colour, flavour, …) | generators $T^a$, structure constants $f^{abc}$, fundamental/adjoint $\delta$ | $N$ (fund.) / $N^2-1$ (adj.) |

A network may mix all three at once; the sectors never interfere because contraction is by
**index label**, not by extent (see [Key concepts](concepts.md)). Scalar coefficients — dot
products, opaque runtime factors — ride along and are emitted as ordinary arithmetic.

## The conventions, exactly

- **Dimension is 4.** Every Lorentz index runs over four components, addressed `0..3`. There is no
  general-$d$ or arbitrary-dimension mode.
- **The metric is the Euclidean $\delta_{\mu\nu}$ — not Minkowski.** There is **no** $-+++$
  signature and no metric factors of $-1$: a scalar product is a plain sum,
  $a\cdot b = \sum_{\mu} a_\mu b_\mu$, and raising/lowering an index is the identity. If your theory
  is written in a Minkowski signature you must translate it to Euclidean first (the engine will not
  do it for you).
- **Dirac gammas are $4\times4$, Hermitian, in the chiral (Weyl) basis**, satisfying
  $\{\gamma^\mu,\gamma^\nu\} = 2\,\delta^{\mu\nu}I$. $\gamma_5=\mathrm{diag}(1,1,-1,-1)$ is
  block-diagonal (free to include), and a trace of an odd number of gammas vanishes structurally.
  The exact matrices are in [Sector data](../internals/sectors.md#dirac-gamma-matrices).
- **Groups are SU($N$) only.** Colour and flavour are both handled by one $N$-parameterised family
  of SU($N$) heads; $N$ is a **compile-time integer** (`SetNc` / `SetNf` in the DSL; a template
  parameter in C++), not a runtime argument. Other Lie groups are not built in.
- **A finished network is a real scalar.** A network with no free indices folds to a scalar
  polynomial in the frame's scalar symbols. A trace over real vectors is real, so the generated
  kernel carries only the **real part**; the imaginary half is dropped when it vanishes (the
  generic case) and kept only if it genuinely survives.
- **Closure is enforced, in every sector.** There is no tensor-valued result: an index that occurs
  only once is a caller error, not a free index, and both contraction cores refuse it. Reach for a
  tensor and you get an exception naming the offending index, never a number. To extract tensor
  information, either **project** — contract with a basis element and emit one scalar per
  $\langle P_i, T_j\rangle$ pair — or **pin the index to a fixed component** (`ntGamma[0, d1, d2]`
  is $\gamma^0$; see [step-19](../tutorials/step-19.md)) and emit one closed kernel per component.

## What it is — and isn't

NumTracer is *general index contraction over the sectors above*: hand it any network built from
those objects, in those conventions, and it folds it to a scalar kernel — the physics is entirely
in the network you write, not in the engine. It is **not** a general-metric or arbitrary-dimension
tensor system, a symbolic CAS, or a solver for anything beyond contracting the network you give it.

Next: [Bring your own network](bring-your-own-network.md) maps a concrete network onto the engine
from both front-ends and runs it end to end.
