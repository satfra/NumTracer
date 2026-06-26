#!/usr/bin/env wolframscript
(* DEMONSTRATION — dressing the u-quark and d-quark DIFFERENTLY (broken flavour symmetry).

   A scalar (sigma) two-point self-energy from a quark loop, summed over TWO explicit quark
   flavours u and d.  Each flavour is its OWN field / diagram and carries its OWN propagator
   dressing (Gu vs Gd) and regulator-dot insertion (GuDot vs GdDot).  There is NO SU(Nf) flavour
   GROUP here: folding flavour as a group collapses a closed flavour loop to the single
   flavour-blind number Nf (sun_net.hpp), which can carry only ONE quark dressing.  When u and d
   are dressed differently SU(Nf) is no longer a symmetry, so the physically-correct
   representation is one diagram per flavour — exactly what this flow builds.  Only colour SU(3)
   is folded (tr(1)=Nc per loop).

   The generated kernel therefore exposes Gu, Gd, GuDot, GdDot as INDEPENDENT runtime dressing
   parameters and returns (u-loop with Gu) + (d-loop with Gd).

   Method mirrors gen_pion_quark.wls: hand-rolled Yukawa rules, a massless q-slash propagator
   numerator + an R_F q-slash (Rdot) insertion so the gamma count per loop is even (4) and the
   trace does not vanish.  Differences: scalar (isoscalar) external instead of the pion triplet,
   two explicit flavour fields instead of one SU(2) doublet, and the per-flavour dressings are
   kept SYMBOLIC (-> runtime kernel params) instead of baked.

   Run from numtracer/:   wolfram -script tests/gen/gen_flavour_split.wls *)

If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];

(* shim opt-in (same as the QCD test setups): emit dressings as DiFfRG::Fn runtime params so the
   compare driver can pass DISTINCT u/d dressing functions at run time. *)
$ntShimOpts = {"RuntimeInclude" -> "shim.hpp", "ExtraIncludes" -> {},
  "KernelNamespace" -> "DiFfRG", "SupportNamespace" -> "DiFfRG",
  "DressingType" -> "SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>"};
SetOptions[MakeNTKernel, Sequence @@ $ntShimOpts];

(* one scalar sigma + two quark flavours (Dirac index di, fundamental-colour index ci). The index
   NAMES are positional slot labels only; the hand-rolled rules below assign their meaning. *)
fields = <|"Commuting" -> {sig[p]},
  "Grassmann" -> {{qub[p, {di, ci}], qu[p, {di, ci}]}, {qdb[p, {di, ci}], qd[p, {di, ci}]}}|>;
truncation = <|
  GammaN -> {{sig, sig}, {sig, qub, qu}, {sig, qdb, qd}, {qub, qu}, {qdb, qd}},
  Propagator -> {{sig, sig}, {qub, qu}, {qdb, qd}},
  R -> {{qub, qu}, {qdb, qd}}, Rdot -> {{qub, qu}, {qdb, qd}}, Field -> {{}}|>;
FSetGlobalSetup[<|"FieldSpace" -> fields, "Truncation" -> truncation,
  "FeynmanRules" -> <|GammaN -> {}, Propagator -> {}, R -> {}, Rdot -> {}|>|>];

fRG = FTakeDerivatives[WetterichEquation, {sig[i1], sig[i2]}] // FTruncate // FRoute;
expr = fRG["1-Loop"]["Expression"];

(* hand-rolled diagrammatic rules, one set per flavour. The ONLY difference between u and d is the
   dressing name attached to the propagator / Rdot insertion.  Yukawa: i gamma5 (Dirac) x colour
   delta.  Quark line: massless q-slash numerator x colour delta x the per-flavour dressing. *)
(* the scalar sigma leg carries a momentum but NO index list, hence the leg pattern {_}. *)
yuk[fb_, f_] := GammaN[{sig, fb, f}, {{_}, {_, {dqb_, cqb_}}, {_, {dq_, cq_}}}] :>
  I gamma5[dqb, dq] deltaFundCol[cqb, cq];
prop[fb_, f_] := Propagator[{f, fb} | {fb, f}, {{m1_, {d1_, c1_}}, {m2_, {d2_, c2_}}}] :>
  With[{mu = Unique["lm"]},
    deltaFundCol[c1, c2] (-I) gamma[mu, d1, d2] vec[m2, mu] * dressing[GammaN, {fb, f}, 1, {m1, m2}]];
rdot[fb_, f_] := Rdot[{fb, f}, {{m1_, {d1_, c1_}}, {m2_, {d2_, c2_}}}] :>
  With[{mu = Unique["lm"]},
    deltaFundCol[c1, c2] (-I) gamma[mu, d1, d2] vec[m2, mu] * dressing[Rdot, {fb, f}, 1, {m1, m2}]];
myRules = {yuk[qub, qu], yuk[qdb, qd], prop[qub, qu], prop[qdb, qd], rdot[qub, qu], rdot[qdb, qd]};

(* per-flavour SYMBOLIC dressings -> become the runtime kernel params Gu/Gd/GuDot/GdDot. *)
dressed = (expr /. myRules) //. {
  dressing[GammaN, {qub, qu}, 1, {m1_, m2_}] :> Gu[Sqrt[sp[m2, m2]]],
  dressing[GammaN, {qdb, qd}, 1, {m1_, m2_}] :> Gd[Sqrt[sp[m2, m2]]],
  dressing[Rdot,   {qub, qu}, 1, {m1_, m2_}] :> GuDot[Sqrt[sp[m2, m2]]],
  dressing[Rdot,   {qdb, qd}, 1, {m1_, m2_}] :> GdDot[Sqrt[sp[m2, m2]]]};
dressed = UseLorentzLinearity[dressed] //. {lf1 -> l1};

SetNc[3];                       (* colour SU(3); NO flavour group (no flavour heads in the net) *)
net = FromFunKit[dressed];

frame = NumTracer`propFrame[Global`p, Global`l1, Global`cos1, Global`p1, Global`l1];
ntk   = NumTrace[net, "Frame" -> frame, "Args" -> {Global`l1, Global`cos1, Global`p}];
dress = {Gu, Gd, GuDot, GdDot};

MakeNTKernel[ntk, FileNameJoin[{$gen, "gen_flavour_split_num.cpp"}],
  FileNameJoin[{$gen, "Flavour_split_num_kernel.hh"}], FileNameJoin[{$gen, "Flavour_split_num_kernels.hh"}],
  "Backend" -> "Numeric", "Name" -> "Flavour_split_num_kernel", "Namespace" -> "flavour_split_num",
  "Dressings" -> dress];
MakeNTKernel[ntk, FileNameJoin[{$gen, "Flavour_split_num_dense_kernel.hh"}],
  "Backend" -> "Dense", "Name" -> "Flavour_split_num_dense_kernel", "Dressings" -> dress];
Print["flavour-split (u vs d distinct dressings) kernels generated."];
