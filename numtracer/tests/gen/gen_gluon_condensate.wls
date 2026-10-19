#!/usr/bin/env wolframscript
(* DEMONSTRATION — a per-component ADJOINT dressing (the gluon-condensation use case): a loop whose
   internal line carries an SU(3)-ADJOINT index dressed only along the CARTAN directions (λ3, λ8),
   with the other 6 colours dropped — instead of the colour-blind δ^{ab}. The adjoint trace folds
   (via sun_value_dressed) to Σ over the surviving components of the per-component dressing product —
   a runtime adjoint-sum token — with the Dirac trace computed ONCE, NOT one diagram per component.
   This is exactly how a gluon condensed in the Cartan subalgebra (Z_A^{a=3}, Z_A^{a=8}) enters.

   Built as a scalar (sigma) two-point from a loop carrying a Dirac index AND an SU(3)-adjoint index
   (the adjoint slot stands in for the gluon colour). Three kernels are emitted, exercising the named
   / Default-collapse / drop paths of the new spec-based diagonal dressing:
     * _plain : the adjoint propagators carry the colour-blind ntSUNDeltaAdj[3, a, b] × ZAflat[scale].
     * _diag  : ntSUNDiagAdj[3, a, b, {7 -> ZA3[scale], 8 -> ZA8[scale], Default -> ZAoff[scale]}] — λ3 & λ8
                named (distinct scalar dressings), the other 6 collapsed into a shared ZAoff.
     * _cartan: ntSUNDiagAdj[3, a, b, {7 -> ZA3[scale], 8 -> ZA8[scale]}] — only λ3 & λ8 survive, the other 6
                DROP (no dead terms). This is the headline use case.
   A spec index is 1-based in NumTracer's generator order (off-diagonal first, Cartan last), so λ3 and
   λ8 are components 7 and 8, not 3 and 8 (component 3 is λ4).
   compare_gluon_condensate.cpp validates: (1) COLLAPSE — feeding _diag one function for all of
   ZA3/ZA8/ZAoff reproduces _plain (Σ over 8 = 8 Z^2); (2) DROP — _cartan equals _diag with ZAoff≡0;
   (3) per-component independence.

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
diagAdj[a1_, a2_, spec_, m_] := NumTracer`ntSUNDiagAdj[3, a1, a2, spec[Sqrt[sp[m, m]]]];
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

(* Cartan-dressed (λ3 = 7 & λ8 = 8 named, other 6 via a shared Default ZAoff) — covers all 8 components. *)
specDiag[nm3_, nm8_, off_]  := Function[s, {7 -> nm3[s], 8 -> nm8[s], Default -> off[s]}];
(* Cartan-only: name λ3 & λ8, DROP the other 6 (no Default). *)
specCartan[nm3_, nm8_]      := Function[s, {7 -> nm3[s], 8 -> nm8[s]}];
netDiag = mkNet[
  diagAdj[#1, #2, specDiag[Global`ZA3, Global`ZA8, Global`ZAoff], #3] &,
  diagAdj[#1, #2, specDiag[Global`ZA3dot, Global`ZA8dot, Global`ZAoffdot], #3] &];
netCartan = mkNet[
  diagAdj[#1, #2, specCartan[Global`ZA3, Global`ZA8], #3] &,
  diagAdj[#1, #2, specCartan[Global`ZA3dot, Global`ZA8dot], #3] &];
(* colour-blind reference: δ^{ab} × a scalar dressing ZAflat / ZAflatDot. *)
netPlain = mkNet[(aAdj[#1, #2] Global`ZAflat[Sqrt[sp[#3, #3]]]) &,
                 (aAdj[#1, #2] Global`ZAflatDot[Sqrt[sp[#3, #3]]]) &];

frame = NumTracer`propFrame[Global`p, Global`l1, Global`cos1, Global`p1, Global`l1];
args  = {Global`l1, Global`cos1, Global`p};

ntkDiag   = NumTrace[netDiag,   "Frame" -> frame, "Args" -> args];
ntkCartan = NumTrace[netCartan, "Frame" -> frame, "Args" -> args];
ntkPlain  = NumTrace[netPlain,  "Frame" -> frame, "Args" -> args];

MakeNTKernel[ntkDiag, FileNameJoin[{$gen, "gen_gluon_condensate_diag.cpp"}],
  FileNameJoin[{$gen, "Gluon_condensate_diag_kernel.hh"}], FileNameJoin[{$gen, "Gluon_condensate_diag_kernels.hh"}],
  "Name" -> "Gluon_condensate_diag_kernel", "Namespace" -> "gluon_condensate_diag",
  "Dressings" -> {Global`ZA3, Global`ZA8, Global`ZAoff, Global`ZA3dot, Global`ZA8dot, Global`ZAoffdot}];
MakeNTKernel[ntkCartan, FileNameJoin[{$gen, "gen_gluon_condensate_cartan.cpp"}],
  FileNameJoin[{$gen, "Gluon_condensate_cartan_kernel.hh"}], FileNameJoin[{$gen, "Gluon_condensate_cartan_kernels.hh"}],
  "Name" -> "Gluon_condensate_cartan_kernel", "Namespace" -> "gluon_condensate_cartan",
  "Dressings" -> {Global`ZA3, Global`ZA8, Global`ZA3dot, Global`ZA8dot}];
MakeNTKernel[ntkPlain, FileNameJoin[{$gen, "gen_gluon_condensate_plain.cpp"}],
  FileNameJoin[{$gen, "Gluon_condensate_plain_kernel.hh"}], FileNameJoin[{$gen, "Gluon_condensate_plain_kernels.hh"}],
  "Name" -> "Gluon_condensate_plain_kernel", "Namespace" -> "gluon_condensate_plain",
  "Dressings" -> {Global`ZAflat, Global`ZAflatDot}];
Print["gluon-condensate (Cartan-dressed ADJOINT) plain + diag + cartan kernels generated."];
