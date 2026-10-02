# Key concepts: how a trace is described

This page builds the mental model the rest of the guide relies on: what a tensor network *is*
in NumTracer, how indices decide what gets summed, where the runtime numbers come from, and
what the generated kernel ends up looking like. It is short on purpose. For the exact sectors and
conventions it assumes — in particular that spacetime is 4D and the metric is the **Euclidean**
$\delta_{\mu\nu}$, not Minkowski — see [Scope & conventions](scope-and-conventions.md).

## A network is tensors with shared indices

A tensor is a list of *indices* plus one **entry** for every combination of index values. You
describe a network as a product of tensor *heads* with index labels; repeated labels are summed
(Einstein convention), and a label that appears once stays free. In the front-end DSL:

```mathematica
(* p · P(l) · p  —  the loop momentum's transverse projector sandwiched between two p's *)
ntVec[qp, mu] ntTransProj[ql, mu, nu] ntVec[qp, nu]
```

Here `mu` and `nu` each appear twice, so both are summed away; the network closes to a scalar.
The DSL heads mirror the FORM/FormTracer vocabulary:

| head | meaning |
|---|---|
| `ntMetric[mu, nu]` | Lorentz metric $\delta_{\mu\nu}$ |
| `ntVec[q, mu]` | momentum component $q_\mu$ |
| `ntTransProj[q, mu, nu]` | transverse projector $P^T_{\mu\nu}(q) = \delta_{\mu\nu} - q_\mu q_\nu / q^2$ |
| `ntLongProj[q, mu, nu]` | longitudinal projector $P^L_{\mu\nu}(q) = q_\mu q_\nu / q^2$ |
| `ntSUNf[N, a, b, c]` | SU(N) structure constant $f^{abc}$ (rank $N$) |
| `ntSUNT[N, a, i, j]` | SU(N) fundamental generator $(T^a)_{ij}$ |
| `ntSUNDeltaAdj[N, a, b]` | SU(N) adjoint $\delta^{ab}$ |
| `ntSUNDeltaFund[N, i, j]` | SU(N) fundamental $\delta_{ij}$ |
| `ntSUNDiagFund[N, i, j, spec]` | a fundamental $\delta_{ij}$ whose components carry different dressings: `spec = {1 -> Zu, 2 -> Zd}` gives $\mathrm{diag}(Z_u, Z_d)$ ([step-17](../tutorials/step-17.md)) |
| `ntSUNDiagAdj[N, a, b, spec]` | the same for an adjoint $\delta^{ab}$ |
| `ntGamma[mu, d1, d2]` | Dirac $\gamma^\mu$ with spinor indices $d_1, d_2$ |
| `ntVec[q, mu] ntGamma[mu, d1, d2]` | a slashed momentum $\slashed q$ (there is no separate slash head) |
| `ntSP[q1, q2]`, `ntDress[h, …]` | scalar coefficients (dot product, opaque dressing) |

The SU(N) heads form **one $N$-parameterized family**: the rank $N$ is always the first
argument, so a single network can carry colour SU($N_c$) and flavour SU($N_f$) at once —
they stay distinct because their indices carry different labels (see below), not because the
heads differ. Finite-temperature work adds the electric/magnetic projectors `ntElectricProj` /
`ntMagneticProj`, the spatial scalar product `ntSPS`, and the spatial *vector* `ntSpatialVec`
(FormTracer's `vecs`, which is what a spatial pslash is built from); see the
[step-16](../tutorials/step-16.md).

## Indices contract by label, not by extent

Every index carries a *label* and an *extent* (how many values it runs over):

| sector | index | extent |
|---|---|---|
| Lorentz (spacetime) | $\mu$ | 4 |
| Dirac | spinor | 4 |
| SU($N$) | adjoint | $N^2-1$ |
| SU($N$) | fundamental | $N$ |

The single rule of the whole engine is: **two indices are summed together exactly when they
carry the same label**, whether they sit on one tensor or on two different ones. Labels of
different sectors never meet, so in the DSL one expression can carry Lorentz, Dirac and SU($N$)
structure at once.

The value of such a network *factorises*: the SU($N$) part shares no index with the rest, so its
value is a number that multiplies the Dirac ⊗ Lorentz contraction. The C++ engine computes the two
pieces separately — `Frame::trace` for Dirac ⊗ Lorentz, `SUN::value` for each SU($N$) group — and
the diagram's value is their product.

In C++ the labels are typed, and they are handed out rather than chosen, so they are distinct by
construction:

```cpp
nt::Frame F;                          // the kinematic frame (next section)
auto [mu, nu] = F.indices<2>();       // nt::LorentzIndex
nt::SUN su3(3);                       // one object per SU(N) group
auto [a, b]   = su3.adjoint<2>();     // nt::AdjIndex
auto [i, j]   = su3.fundamental<2>(); // nt::FundIndex
// nt::gamma(mu), nt::vec(nu, p), su3.T(a, i, j), …
```

A label of the wrong kind in a slot does not compile (`su3.T(a, a, j)`), and a label of another
SU($N$) group throws.

## Runtime numbers come from the frame

A tensor head holds no numbers of its own. To get numbers in, you place each momentum in a
chosen reference **frame**, which fixes its components. For a one-angle loop integrand the
natural choice is:

```text
        axis 1
          ^
          |        l = (l·cosθ, l·sinθ, 0, 0)   -> components {0,1} present
       l /
        /  θ
   -----+--------> axis 0
        p = (|p|, 0, 0, 0)                      -> only component 0 present
```

The frame is the bridge between the symbolic network and the runtime kernel: it determines how
each scalar product (`l·p`, `l²`, …) is written in terms of the kernel's actual arguments
($|l|$, $\cos\theta$, $|p|$, …). The generated kernel evaluates those scalar symbols once per
call from the runtime arguments, then runs the lowered arithmetic.

In C++ the same frame is an `nt::Frame` — the symbols, the momenta built from them, and (filled in
automatically) the denominators $1/k^2$ of every projector:

```cpp
nt::Frame F;
auto P = F.symbol("p"), L = F.symbol("l");
auto [C, S] = F.angle("theta");             // sin θ is derived from cos θ
nt::Momentum p = F.momentum(P, 0, 0, 0);
nt::Momentum l = F.momentum(L * C, L * S, 0, 0);
```

In the Mathematica front-end, frame builders such as `propFrame` do this for you
([step-07](../tutorials/step-07.md)).

## What you get back

A network with no free indices is a single number — a scalar. A network with no runtime
content at all (a pure SU($N$) factor, a pure gamma trace) folds to a compile-time constant. A
network carrying momenta and dressings becomes a small polynomial in the frame's scalar
symbols, which the codegen lowers to a flat, straight-line kernel: trace functions over a few
symbols, plus the per-diagram assembly.

Next, the [tutorials](../tutorials/index.md) build a complete physical integrand — the quark
self-energy — from a network in the DSL all the way to a generated, validated kernel.
