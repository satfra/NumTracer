# Per-flavour and per-component dressings

The flows so far dressed every quark the same way — one propagator dressing $Z_q$ for the
whole colour/flavour multiplet. Physics often needs more: the **u** and **d** quarks dressed
differently (broken isospin), or a colour component dressed on its own (a gluon condensate).
NumTracer supports this two ways, and the choice matters for how big the kernel gets.

| approach | how | when |
|---|---|---|
| **per-flavour split** | one diagram (and one dressing) per flavour | a handful of fields, dressed independently; isospin is *not* a symmetry |
| **group-diagonal dressing** | one trace, a $\delta$ tagged with a per-component dressing array | dress components *within* an SU($N$) group without splitting the trace |

The dividing line: folding a closed flavour loop as an SU($N_f$) **group** collapses it to the
single flavour-blind number $N_f$ — which can carry only *one* dressing. So to give components
different dressings you either avoid the group (split into separate diagrams) or keep the group
but make its $\delta$ carry a per-component dressing array.

## Per-flavour split

Here u and d are two **separate fields**, each its own diagram with its own dressing. There is
no SU($N_f$) flavour group at all — only colour SU(3) is folded. The generator
`tests/gen/gen_flavour_split.wls` builds a scalar (sigma) two-point self-energy from a quark
loop, summed over the two flavours, with hand-rolled Yukawa rules:

```mathematica
(* one scalar sigma + TWO quark flavours, each a field of its own *)
fields = <|"Commuting" -> {sig[p]},
  "Grassmann" -> {{qub[p, {di, ci}], qu[p, {di, ci}]},
                  {qdb[p, {di, ci}], qd[p, {di, ci}]}}|>;

(* the only difference between u and d is the dressing name on the propagator / Rdot insertion *)
dressed = (expr /. myRules) //. {
  dressing[GammaN, {qub, qu}, 1, {m1_, m2_}] :> Gu[Sqrt[sp[m2, m2]]],
  dressing[GammaN, {qdb, qd}, 1, {m1_, m2_}] :> Gd[Sqrt[sp[m2, m2]]],
  dressing[Rdot,   {qub, qu}, 1, {m1_, m2_}] :> GuDot[Sqrt[sp[m2, m2]]],
  dressing[Rdot,   {qdb, qd}, 1, {m1_, m2_}] :> GdDot[Sqrt[sp[m2, m2]]]};

SetNc[3];                          (* colour SU(3); no flavour group — no flavour heads in the net *)
net = FromFunKit[dressed];
ntk = NumTrace[net, "Frame" -> frame, "Args" -> {l1, cos1, p}];
MakeNTKernel[ntk, "gen_flavour_split_num.cpp",
  "Flavour_split_num_kernel.hh", "Flavour_split_num_kernels.hh",
  "Dressings" -> {Gu, Gd, GuDot, GdDot}];        (* four INDEPENDENT runtime dressings *)
```

The kernel returns `(u-loop with Gu) + (d-loop with Gd)` and exposes `Gu, Gd, GuDot, GdDot`
as independent runtime parameters. `compare_flavour_split.cpp` (ctest `flow_flavour_split`)
validates it against the dense oracle.

## Group-diagonal dressing

The same physics, but built as a **single** quark loop carrying an SU(2) flavour index whose
propagators carry a *flavour-diagonal* dressing $\mathrm{diag}(Z[u], Z[d])$ instead of a
flavour-blind $\delta$. The new heads are

```mathematica
ntSUNDiagFund[N, i, j, Z, s]   (* fundamental δ_ij tagged with per-component dressing Z[i](s) *)
ntSUNDiagAdj [N, a, b, Z, s]   (* adjoint     δ^ab tagged with per-component dressing Z[a](s) *)
```

The SU($N$) trace then folds — via `sun_value_dressed` (`network/sun_net.hpp`) — to
$\sum_a Z_a(s)$ over the runtime per-component dressing array, **and the Dirac/colour trace is
computed once**, not once per component. From `tests/gen/gen_flavour_ingroup.wls`:

```mathematica
fields = <|"Commuting" -> {sig[p]},
  "Grassmann" -> {{qb[p, {di, ci, fi}], q[p, {di, ci, fi}]}}|>;   (* ONE flavour doublet *)

diagFlav[fb_, f_, name_, m_] := NumTracer`ntSUNDiagFund[2, fb, f, name, Sqrt[sp[m, m]]];
myRules = {
  (* Yukawa vertex: flavour-BLIND δ (deltaFundFlav -> ntSUNDeltaFund) *)
  GammaN[{sig, qb, q}, …] :> I gamma5[dqb, dq] deltaFundCol[cqb, cq] ntSUNDeltaFund[2, fqb, fq],
  (* propagator: flavour-DIAGONAL dressing diag(G[u], G[d]) in place of the blind δ × dressing *)
  Propagator[{q, qb} | {qb, q}, …] :> deltaFundCol[c1, c2] diagFlav[f1, f2, G, m2] (-I) gamma[mu, d1, d2] vec[m2, mu],
  Rdot[{qb, q}, …]               :> deltaFundCol[c1, c2] diagFlav[f1, f2, Gdot, m2] (-I) gamma[mu, d1, d2] vec[m2, mu]};

MakeNTKernel[ntk, …, "Dressings" -> {G, Gdot}];   (* G, Gdot are ARRAYS indexed by component *)
```

This is the in-group counterpart of the split flow and **reproduces it exactly**:
`compare_flavour_ingroup.cpp` (ctest `flow_flavour_ingroup`) checks the two kernels agree
pointwise with `G[0]=Gu, G[1]=Gd, Gdot[0]=GuDot, Gdot[1]=GdDot`.

### The adjoint case: a gluon condensate

`ntSUNDiagAdj` does the same for an adjoint index — the natural home of a condensed gluon
dressing $Z_A^a$, with $N^2-1 = 8$ components for SU(3). `tests/gen/gen_gluon_condensate.wls`
emits two kernels: `_diag`, whose adjoint propagators carry `ntSUNDiagAdj[3, a, b, ZA, scale]`,
and `_plain`, whose propagators carry the colour-blind `ntSUNDeltaAdj[3, a, b]` times a scalar.
The test (`flow_gluon_condensate`) validates the **collapse property**: feeding the `_diag`
kernel an 8-component array whose components are all equal reproduces the `_plain` kernel
exactly ($\sum_a Z^2 = 8\,Z^2$) — i.e. the component-resolved path reduces to the colour-blind
one when there is nothing to resolve.

```{note}
Group-diagonal dressings leave `sun_value_cx` (the plain colour fold) untouched, so any flow
that does *not* use them regenerates byte-identical.
```

## Verify

```bash
cd numtracer
# regenerate (needs wolframscript + FunKit + FORM):
wolfram -script tests/gen/gen_flavour_split.wls
wolfram -script tests/gen/gen_flavour_ingroup.wls
wolfram -script tests/gen/gen_gluon_condensate.wls
# validate the committed kernels (no Wolfram needed):
cmake --build build -j4
ctest -R "flow_flavour_split|flow_flavour_ingroup|flow_gluon_condensate" --output-on-failure
```

For the design — how the colour fold keeps a group-diagonal $\delta$ as a polynomial and how a
diagram that disconnects into several traces is emitted as a product — see
[the numeric engine](../internals/numeric-engine.md) and [front-end & codegen](../internals/codegen.md).
