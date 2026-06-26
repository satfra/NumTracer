# Key concepts: how a trace is described

This page builds the mental model the rest of the guide relies on: what a tensor network *is*
in NumTracer, how indices decide what gets summed, where the runtime numbers come from, and
what the generated kernel ends up looking like. It is short on purpose.

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
| `ntSUNDiagFund[N, i, j, Z, s]` / `ntSUNDiagAdj[N, a, b, Z, s]` | group-diagonal $\delta$ carrying a per-component dressing $Z[a](s)$ (see [dressed-flavour](../tutorials/dressed-flavour.md)) |
| `ntGamma`, slashed momenta | Dirac gamma matrices and $\slashed p = p^\mu\gamma_\mu$ |
| `ntSP[q1, q2]`, `ntDress[h, …]` | scalar coefficients (dot product, opaque dressing) |

The SU(N) heads form **one $N$-parameterized family**: the rank $N$ is always the first
argument, so a single network can carry colour SU($N_c$) and flavour SU($N_f$) at once —
they stay distinct because their indices carry different labels (see below), not because the
heads differ. Finite-temperature work adds the electric/magnetic projectors `ntElectricProj` /
`ntMagneticProj` and the spatial scalar product `ntSPS`; see the
[finite-temperature tutorial](../tutorials/finite-temperature.md).

## Indices contract by label, not by extent

Every index carries an integer *label* and an *extent* (how many values it runs over):

| sector | index | extent |
|---|---|---|
| Lorentz (spacetime) | $\mu$ | 4 |
| Dirac | spinor | 4 |
| SU($N$) colour | adjoint | $N^2-1$ |
| SU($N$) colour | fundamental | $N$ |

The single rule of the whole engine is: **two indices are summed together exactly when they
carry the same label**, whether they sit on one tensor or on two different ones. Because
matching is purely by label and the extent travels with each index, **different sectors never
interfere** — a Lorentz index and a colour index can sit in the same network and will only
contract if you deliberately give them the same label. So one network can carry
Lorentz ⊗ Dirac ⊗ SU(N) at once, and you control the summation entirely by which labels you
reuse.

In the C++ engine those labels are plain `int`s. Rather than scatter raw numbers through the
builders, the [tutorials](../tutorials/index.md) name them with a single unscoped `enum`, so the
contraction reads in index names while the compiler hands out a distinct value to each:

```cpp
enum {
  mu, nu,    // Lorentz
  A, B,      // fundamental colour
  a, b       // adjoint colour
};
// ... net::dgamma(mu), nm::nvec(nu, {{1.0, 0}}), SUN::T(3, a, A, B), ...
```

One enum across all sectors makes every label unique by construction — exactly the property the
contraction rule needs, with no hand-picked offsets to keep the sectors from colliding.

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

## What you get back

A network with no free indices is a single number — a scalar. A network with no runtime
content at all (a pure colour factor, a pure gamma trace) folds to a compile-time constant. A
network carrying momenta and dressings becomes a small polynomial in the frame's scalar
symbols, which the codegen lowers to a flat, straight-line kernel. The end product has the same
shape as a FORM reference kernel: trace functions over a few symbols, plus the per-diagram
assembly.

Next, the [tutorials](../tutorials/index.md) build a complete physical integrand — the quark
self-energy — from a network in the DSL all the way to a generated, validated kernel.
