#!/usr/bin/env wolframscript
(* DEMONSTRATION — dressing u and d quarks differently WITHIN the SU(2) flavour group (no per-flavour
   diagram split).  Same physics as gen_flavour_split.wls (a scalar sigma two-point from a quark
   loop), but built as a SINGLE quark loop carrying an SU(2) flavour index, whose propagators carry a
   FLAVOUR-DIAGONAL dressing diag(G[0],G[1]) instead of a flavour-blind δ.  The SU(2) flavour trace
   then folds (via sun_value_dressed) to Σ_{f∈{u,d}} (per-flavour dressing product) — the
   per-component dressings G[f]/Gdot[f] flow through the group as a runtime colour-sum token, and the
   Dirac/colour trace is computed ONCE (not once per flavour, not as two fields).

   This is the in-group counterpart of the separate-diagram flow_flavour_split, and MUST reproduce it
   exactly: compare_flavour_ingroup.cpp checks the two kernels agree pointwise with G[0]=gu, G[1]=gd,
   Gdot[0]=gudot, Gdot[1]=gddot.

   Run from numtracer/:   wolfram -script tests/gen/gen_flavour_ingroup.wls *)

If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];

$ntShimOpts = {"RuntimeInclude" -> "shim.hpp", "ExtraIncludes" -> {},
  "KernelNamespace" -> "DiFfRG", "SupportNamespace" -> "DiFfRG",
  "DressingType" -> "SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>"};
SetOptions[MakeNTKernel, Sequence @@ $ntShimOpts];

(* one scalar sigma + one quark flavour-DOUBLET (Dirac di, colour ci, flavour fi). *)
fields = <|"Commuting" -> {sig[p]},
  "Grassmann" -> {{qb[p, {di, ci, fi}], q[p, {di, ci, fi}]}}|>;
truncation = <|
  GammaN -> {{sig, sig}, {sig, qb, q}, {qb, q}},
  Propagator -> {{sig, sig}, {qb, q}},
  R -> {{qb, q}}, Rdot -> {{qb, q}}, Field -> {{}}|>;
FSetGlobalSetup[<|"FieldSpace" -> fields, "Truncation" -> truncation,
  "FeynmanRules" -> <|GammaN -> {}, Propagator -> {}, R -> {}, Rdot -> {}|>|>];

fRG = FTakeDerivatives[WetterichEquation, {sig[i1], sig[i2]}] // FTruncate // FRoute;
expr = fRG["1-Loop"]["Expression"];

(* hand-rolled rules. The Yukawa vertex carries a FLAVOUR-blind δ (deltaFundFlav, a plain rank-2
   fundamental Kronecker emitted directly as the NumTracer head). Each quark propagator / Rdot
   insertion carries a FLAVOUR-DIAGONAL dressing ntSUNDiagFund[2, f1, f2, name, scale] — the per-
   flavour propagator dressing diag(name[u], name[d]) — in place of the flavour-blind δ × scalar
   dressing of flow_flavour_split. Colour (deltaFundCol) and Dirac (gamma5/gamma) are unchanged. *)
dFlav[fb_, f_] := NumTracer`ntSUNDeltaFund[2, fb, f];
diagFlav[fb_, f_, name_, m_] := NumTracer`ntSUNDiagFund[2, fb, f, name, Sqrt[sp[m, m]]];
myRules = {
  GammaN[{sig, qb, q}, {{_}, {_, {dqb_, cqb_, fqb_}}, {_, {dq_, cq_, fq_}}}] :>
    I gamma5[dqb, dq] deltaFundCol[cqb, cq] dFlav[fqb, fq],
  Propagator[{q, qb} | {qb, q}, {{m1_, {d1_, c1_, f1_}}, {m2_, {d2_, c2_, f2_}}}] :>
    With[{mu = Unique["lm"]},
      deltaFundCol[c1, c2] diagFlav[f1, f2, Global`G, m2] (-I) gamma[mu, d1, d2] vec[m2, mu]],
  Rdot[{qb, q}, {{m1_, {d1_, c1_, f1_}}, {m2_, {d2_, c2_, f2_}}}] :>
    With[{mu = Unique["lm"]},
      deltaFundCol[c1, c2] diagFlav[f1, f2, Global`Gdot, m2] (-I) gamma[mu, d1, d2] vec[m2, mu]]};

dressed = UseLorentzLinearity[expr /. myRules] //. {lf1 -> l1};

SetNc[3];
net = FromFunKit[dressed];   (* colour SU(3); flavour heads carry literal rank 2, passed through *)

frame = NumTracer`propFrame[Global`p, Global`l1, Global`cos1, Global`p1, Global`l1];
ntk   = NumTrace[net, "Frame" -> frame, "Args" -> {Global`l1, Global`cos1, Global`p}];
dress = {Global`G, Global`Gdot};

MakeNTKernel[ntk, FileNameJoin[{$gen, "gen_flavour_ingroup_num.cpp"}],
  FileNameJoin[{$gen, "Flavour_ingroup_num_kernel.hh"}], FileNameJoin[{$gen, "Flavour_ingroup_num_kernels.hh"}],
  "Backend" -> "Numeric", "Name" -> "Flavour_ingroup_num_kernel", "Namespace" -> "flavour_ingroup_num",
  "Dressings" -> dress];
Print["flavour-ingroup (u vs d via SU(2)-diagonal dressing) numeric kernel generated."];
