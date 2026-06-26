# Finite temperature

At finite temperature the Euclidean $O(4)$ rotation symmetry breaks to spatial $O(3)$: the
heat bath singles out a rest frame $u^\mu = (1,0,0,0)$, and the **zeroth (Euclidean-time /
Matsubara) component** of every momentum becomes physically distinct from the spatial three. In
the Matsubara formulation the metric stays the Euclidean $\delta_{\mu\nu}$ — what changes is that a
momentum now carries an independent temporal component (a Matsubara frequency) alongside its spatial
3-vector.

NumTracer already contracts on a momentum's explicit frame components, so a distinct zeroth component
costs nothing extra — you just choose a frame whose temporal slot is an independent symbol. This page
continues the running [quark self-energy](full-diagram.md) example into finite $T$, introducing the
four pieces the broken symmetry needs:

* the **electric** $P_E$ and **magnetic** $P_M$ projectors the gluon transverse space splits into;
* the **spatial scalar product** `ntSPS` ($\vec a\cdot\vec b$);
* the **temporal component** `ntVec[q, 0]` (a scalar, e.g. the Matsubara frequency $\pi T$);
* a **finite-$T$ frame** with independent temporal components.

```{admonition} Convention
:class: note
Component index `0` is the Euclidean-time / Matsubara direction (heat bath $u=(1,0,0,0)$); `1,2,3`
are spatial. This is the contract between the finite-$T$ frame and the electric/magnetic projectors.
```

## The projectors

The 4d transverse projector is unchanged and remains valid at finite $T$:

$$ P^T_{\mu\nu}(k) = \delta_{\mu\nu} - \frac{k_\mu k_\nu}{k^2}. $$

It splits into a **magnetic** (spatial-transverse) and an **electric** (time-like-transverse) part,

$$
P^M_{ij}(k) = \delta_{ij} - \frac{k_i k_j}{|\vec k|^2}\ (i,j\in\{1,2,3\}),\quad P^M_{0\nu}=P^M_{\mu0}=0,
\qquad
P^E = P^T - P^M,
$$

so that $P^E + P^M = P^T$, with $\operatorname{tr}P^M = 2$ and $\operatorname{tr}P^E = 1$. The
magnetic projector uses the **spatial** denominator $|\vec k|^2 = k^2 - k_0^2$; the electric one uses
both $1/k^2$ and $1/|\vec k|^2$. In the numeric engine these are the factors `nprojE` and `nprojM`,
each carrying the loop momentum and the relevant inverse-atom ids (see
[`numeric_contract.hpp`](../internals/index.md)); in the dense oracle they are
`electric_projector` / `magnetic_projector`.

## The example: a finite-$T$ quark self-energy

In the heat-bath frame the gluon line splits into electric and magnetic parts, so the wave-function
numerator becomes a Dirac trace contracted with $G_{\mu\nu}(l) = Z_E\,P^E_{\mu\nu}(l) +
Z_M\,P^M_{\mu\nu}(l)$:

$$
N(p,q,l) = \mathrm{tr}\!\big[\slashed p\,\gamma^\mu\,\slashed q\,\gamma^\nu\big]\,G_{\mu\nu}(l),
\qquad q = l-p,
$$

with the external quark at the lowest fermionic Matsubara frequency $p_0 = \pi T$. The runnable
program is [`tests/test_finiteT_quark.cpp`](https://github.com/) — its core builds the trace through
the numeric engine and checks it against the closed-form trace identity:

```cpp
// 8 symbols: external p = vars 0..3 (vid 0), loop l = vars 4..7 (vid 1).
// q = l − p is a linear combination; atoms: 0 = l², 1 = |l⃗|² (component 0 dropped).
const network::DiracNet chain = {network::dslash({{1.0, pVid}}), network::dgamma(100),
                                 network::dslash({{1.0, lVid}, {-1.0, pVid}}), network::dgamma(101)};
// finite-T gluon line G_{μν}(l) = ZAE·P_E(l) + ZAM·P_M(l)
const nm::NNet gluon = {nm::NTerm{Cx{ZAE, 0}, {nm::nprojE(100, 101, {{1.0, lVid}}, 0, 1)}},
                        nm::NTerm{Cx{ZAM, 0}, {nm::nprojM(100, 101, {{1.0, lVid}}, 1)}}};
const nm::MPoly N = nm::numeric_value(nsym, chain, gluon, comp, atomDen);
```

| piece | meaning |
|---|---|
| `dslash({{1.0,1},{-1.0,0}})` | the slashed internal momentum $\slashed q$ with $q = l - p$ (a linear combination) |
| `nprojE(100,101,{{1.0,1}},0,1)` | $P_E(l)$: legs `100`/`101`, momentum $l$ (vid 1), full atom 0 ($l^2$), spatial atom 1 ($|\vec l|^2$) |
| `nprojM(100,101,{{1.0,1}},1)` | $P_M(l)$: only the spatial atom 1 |
| `atomDen = {l², |l⃗|²}` | atom 0 is $l^2$, atom 1 is $|\vec l|^2 = \sum_{\mu=1}^{3} l_\mu^2$ |

Evaluated over random kinematics with $p_0 = \pi T$ (and an independent loop energy $l_0$), the numeric
engine agrees with the closed form $N = 4\big[\,2\,(p\cdot G\cdot q) - (p\cdot q)\operatorname{tr}G\,\big]$
to machine precision:

```text
== finite-T quark self-energy (5000 random points, vec[p,0]=πT) ==
  numeric engine vs closed-form trace identity:  worst relative = 1.742e-11

ALL TESTS PASSED
```

## From the Mathematica front-end

The same objects exist in the DSL, so a generated finite-$T$ kernel needs only a finite-$T$ frame and
the new heads:

```mathematica
(* frames are now public — no NumTracer`Private` qualifier needed *)
frame = propFrameFT[p0, p, l0, l1, cos1, pp, ll];   (* component 0 = Matsubara/temporal *)

net = ntVec[pp, i1] ntElectricProj[ll, i1, i2] ntVec[pp, i2]            (* pp·P_E(ll)·pp *)
    + ntVec[pp, j1] ntMagneticProj[ll, j1, j2] ntVec[pp, j2]            (* pp·P_M(ll)·pp *)
    + ntSPS[pp, ll] ntVec[ll, 0] *                                       (* spatial dot × temporal comp *)
        ntVec[pp, k1] ntTransProj[ll, k1, k2] ntVec[pp, k2];

ntk = NumTrace[net, "Frame" -> frame, "Args" -> {p0, p, l0, l1, cos1}];
```

* `propFrameFT[p0, p, l0, l1, cos1, pp, ll]` sets `pp = {p0, p, 0, 0}` (with `p0` the Matsubara
  frequency $\pi T$) and `ll = {l0, l1 cos1, l1 sqrt(1-cos1^2), 0}` — independent temporal slots.
  Because `ll`'s temporal component is not proportional to the loop magnitude, the numeric backend
  automatically takes the general frame path (the compact unit-loop optimisation applies only to the
  $O(4)$-symmetric vacuum frames).
* `ntElectricProj` / `ntMagneticProj` are the new projector heads (lowered to `eproj`/`mproj` in the
  numeric backend and `electric_projector`/`magnetic_projector` in the dense one, each allocating the
  spatial $1/|\vec q|^2$ env slot automatically).
* `ntSPS[a, b]` is the spatial scalar product $\vec a\cdot\vec b$, resolved by the frame as the
  components-`1..3` dot — a scalar coefficient like `ntSP`.
* `ntVec[q, 0]` — an **integer** second argument — is the scalar component $q_0$ (here the temporal /
  Matsubara component), resolved by the frame; a *symbolic* label (`ntVec[q, mu]`) is still the usual
  tensor leg. A FunKit flow's `sps[...]` and `vec[q,0]` import straight through `FromFunKit`.

The generated kernel is validated numeric-vs-dense by the `ftproj_num` test
(`tests/gen/gen_ftproj_numeric.wls`).

```{admonition} Physics dressings
:class: tip
A full finite-$T$ flow additionally needs finite-$T$ dressings/regulators (functions of $l_0$ and
$|\vec l|$ separately) in the per-theory setup; those are a modelling choice and do not affect the
tensor algebra shown here.
```
