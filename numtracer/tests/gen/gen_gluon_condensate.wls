#!/usr/bin/env wolframscript
(* DEMONSTRATION — a per-component ADJOINT dressing (the gluon-condensation use case): a loop whose
   internal line carries an SU(3)-ADJOINT index dressed component-by-component, diag(Z[0],..,Z[7]),
   instead of the colour-blind δ^{ab}. The adjoint trace folds (via sun_value_dressed) to
   Σ_{a=0}^{7} (per-component dressing product) — a runtime adjoint-sum token — with the Dirac trace
   computed ONCE, NOT 8 separate diagrams. This is exactly how a condensed gluon Z_A^a would enter.

   Built as a scalar (sigma) two-point from a loop carrying a Dirac index AND an SU(3)-adjoint index
   (the adjoint slot stands in for the gluon colour), so the adjoint structure exercises the
   N^2-1 = 8 -dimensional diagonal-dressing codegen path end to end. Two kernels are emitted:
     * _diag : the adjoint propagators carry ntSUNDiagAdj[3, a, b, ZA, scale]  (per-component);
     * _plain: the adjoint propagators carry the colour-blind ntSUNDeltaAdj[3, a, b] × ZAflat[scale].
   compare_gluon_condensate.cpp validates the COLLAPSE: feeding the _diag kernel an array whose 8
   components are all the same function ZAflat reproduces the _plain kernel exactly (Σ_a Z^2 = 8 Z^2).

   Run from numtracer/:   wolfram -script tests/gen/gen_gluon_condensate.wls *)

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

(* sigma + a loop field with a Dirac index di, a colour index ci, and an SU(3) ADJOINT index ai. *)
fields = <|"Commuting" -> {sig[p]},
  "Grassmann" -> {{qb[p, {di, ci, ai}], q[p, {di, ci, ai}]}}|>;
truncation = <|
  GammaN -> {{sig, sig}, {sig, qb, q}, {qb, q}},
  Propagator -> {{sig, sig}, {qb, q}},
  R -> {{qb, q}}, Rdot -> {{qb, q}}, Field -> {{}}|>;
FSetGlobalSetup[<|"FieldSpace" -> fields, "Truncation" -> truncation,
  "FeynmanRules" -> <|GammaN -> {}, Propagator -> {}, R -> {}, Rdot -> {}|>|>];

fRG = FTakeDerivatives[WetterichEquation, {sig[i1], sig[i2]}] // FTruncate // FRoute;
expr = fRG["1-Loop"]["Expression"];

(* the adjoint structure: vertex carries the colour-blind adjoint δ^{ab}; each propagator/Rdot
   carries either a per-component adjoint dressing (diag) or the colour-blind δ × a scalar (plain). *)
aAdj[ab_, a_] := NumTracer`ntSUNDeltaAdj[3, ab, a];
diagAdj[a1_, a2_, name_, m_] := NumTracer`ntSUNDiagAdj[3, a1, a2, name, Sqrt[sp[m, m]]];
vtx = GammaN[{sig, qb, q}, {{_}, {_, {dqb_, cqb_, aqb_}}, {_, {dq_, cq_, aq_}}}] :>
        I gamma5[dqb, dq] deltaFundCol[cqb, cq] aAdj[aqb, aq];
prop[adjStruct_] := Propagator[{q, qb} | {qb, q}, {{m1_, {d1_, c1_, a1_}}, {m2_, {d2_, c2_, a2_}}}] :>
  With[{mu = Unique["lm"]}, deltaFundCol[c1, c2] adjStruct[a1, a2, m2] (-I) gamma[mu, d1, d2] vec[m2, mu]];
rdot[adjStruct_] := Rdot[{qb, q}, {{m1_, {d1_, c1_, a1_}}, {m2_, {d2_, c2_, a2_}}}] :>
  With[{mu = Unique["lm"]}, deltaFundCol[c1, c2] adjStruct[a1, a2, m2] (-I) gamma[mu, d1, d2] vec[m2, mu]];

SetNc[3];
mkNet[propAdj_, rdotAdj_] := Module[{d},
  d = UseLorentzLinearity[expr /. {vtx, prop[propAdj], rdot[rdotAdj]}] //. {lf1 -> l1};
  FromFunKit[d]];

(* per-component diagonal dressings ZA[a], ZAdot[a] (the condensed-gluon dressings). *)
netDiag = mkNet[diagAdj[#1, #2, Global`ZA, #3] &, diagAdj[#1, #2, Global`ZAdot, #3] &];
(* colour-blind reference: δ^{ab} × a scalar dressing ZAflat / ZAflatDot. *)
netPlain = mkNet[(aAdj[#1, #2] Global`ZAflat[Sqrt[sp[#3, #3]]]) &,
                 (aAdj[#1, #2] Global`ZAflatDot[Sqrt[sp[#3, #3]]]) &];

frame = NumTracer`propFrame[Global`p, Global`l1, Global`cos1, Global`p1, Global`l1];
args  = {Global`l1, Global`cos1, Global`p};

ntkDiag  = NumTrace[netDiag,  "Frame" -> frame, "Args" -> args];
ntkPlain = NumTrace[netPlain, "Frame" -> frame, "Args" -> args];

MakeNTKernel[ntkDiag, FileNameJoin[{$gen, "gen_gluon_condensate_diag.cpp"}],
  FileNameJoin[{$gen, "Gluon_condensate_diag_kernel.hh"}], FileNameJoin[{$gen, "Gluon_condensate_diag_kernels.hh"}],
  "Backend" -> "Numeric", "Name" -> "Gluon_condensate_diag_kernel", "Namespace" -> "gluon_condensate_diag",
  "Dressings" -> {Global`ZA, Global`ZAdot}];
MakeNTKernel[ntkPlain, FileNameJoin[{$gen, "gen_gluon_condensate_plain.cpp"}],
  FileNameJoin[{$gen, "Gluon_condensate_plain_kernel.hh"}], FileNameJoin[{$gen, "Gluon_condensate_plain_kernels.hh"}],
  "Backend" -> "Numeric", "Name" -> "Gluon_condensate_plain_kernel", "Namespace" -> "gluon_condensate_plain",
  "Dressings" -> {Global`ZAflat, Global`ZAflatDot}];
Print["gluon-condensate (per-component ADJOINT dressing) diag + plain kernels generated."];
