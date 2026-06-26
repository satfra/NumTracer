#!/usr/bin/env wolframscript
(* Phase-2 example 3 — FUNDAMENTAL COLOUR via the size-guard NUMERIC path.
   The gluon self-energy from a quark loop (gluon -> quark loop -> gluon) with the quark-gluon
   vertex  i gamma^mu T^a.  Tracing it gives a Dirac trace, a trivial flavour loop (delta -> Nf),
   and a FUNDAMENTAL-colour trace tr(T^{a} T^{b}) = (1/2) delta^{ab}.

   The point of this example: the inv (build-time) backend is generated with the et-vs-numeric size
   guard forced to ZERO, so the fundamental colour is contracted by the NUMERIC SUNNet path
   (sun_value_cx in et/inv/col.hpp) instead of an et `constexpr`. The compare test then checks
   inv (numeric colour) == ET (et colour) == Dense — landing the numeric fundamental-colour path on a
   real flow. (The ET/Dense backends ignore the threshold and always fold constants via the et engine.)

   Method: derive with FunKit, hand-roll the diagrammatic rules. Massless quark q-slash propagator;
   R_F ~ q-slash so the Rdot insertion carries the even 4th gamma.

   Run from numtracer/:   wolfram -script tests/gen/gen_glu_quark.wls *)

If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];

fields = <|"Commuting" -> {A[p, {v, c}]}, "Grassmann" -> {{qb[p, {d, A, F}], q[p, {d, A, F}]}}|>;
truncation = <|GammaN -> {{A, A}, {A, qb, q}, {qb, q}}, Propagator -> {{A, A}, {qb, q}},
  R -> {{qb, q}}, Rdot -> {{qb, q}}, Field -> {{}}|>;
FSetGlobalSetup[<|"FieldSpace" -> fields, "Truncation" -> truncation,
  "FeynmanRules" -> <|GammaN -> {}, Propagator -> {}, R -> {}, Rdot -> {}|>|>];

fRG = FTakeDerivatives[WetterichEquation, {A[i1], A[i2]}] // FTruncate // FRoute;
expr = fRG["1-Loop"]["Expression"];

(* hand-rolled rules. QG vertex: i gamma^mu (Dirac+Lorentz) x T^a (fundamental colour) x delta (trivial
   flavour -> Nf, no tau). Quark prop / Rdot: massless q-slash x delta_color x delta_flav. *)
myRules = {
  GammaN[{A, qb, q}, {{_, {lor_, adjc_}}, {_, {dqb_, cqb_, fqb_}}, {_, {dq_, cq_, fq_}}}] :>
    I gamma[lor, dqb, dq] TCol[adjc, cqb, cq] deltaFundFlav[fqb, fq],
  Propagator[{q, qb} | {qb, q}, {{m1_, {d1_, c1_, f1_}}, {m2_, {d2_, c2_, f2_}}}] :>
    With[{mu = Unique["lm"]},
      deltaFundCol[c1, c2] deltaFundFlav[f1, f2] *
      (-I) gamma[mu, d1, d2] vec[m2, mu] * dressing[GammaN, {qb, q}, 1, {m1, m2}]],
  Rdot[{qb, q}, {{m1_, {d1_, c1_, f1_}}, {m2_, {d2_, c2_, f2_}}}] :>
    With[{mu = Unique["lm"]},
      deltaFundCol[c1, c2] deltaFundFlav[f1, f2] *
      (-I) gamma[mu, d1, d2] vec[m2, mu] * dressing[Rdot, {qb, q}, 1, {m1, m2}]]};

(* project the external gluon (i1,i2): transverse projector / (d-1) on Lorentz, adjoint delta /
   (Nc^2-1) on colour;  1/((4-1)(3^2-1)) = 1/24. *)
trace = FTerm[1/24, ntTransProjTok[p1, v1, v2] deltaAdjCol[c1, c2]] ** (expr /. myRules) /.
  ntTransProjTok -> transProj;

dressed = trace //. {
  dressing[GammaN, {qb, q}, 1, {m1_, m2_}] :> 1/(sp[m2, m2] + RF[k^2, sp[m2, m2]]),
  dressing[Rdot, {qb, q}, 1, {m1_, m2_}] :> RFdot[k^2, sp[m2, m2]]};
dressed = UseLorentzLinearity[dressed] //. {lf1 -> l1};

(* SetNc/SetNf before FromFunKit: the colour rank Nc is baked into the SU(N) heads at import,
   and the closed quark flavour loop folds to Nf = 2 (a scalar; no isospin structure here). *)
SetNc[3];
SetNf[2];
net = FromFunKit[dressed];

frame = NumTracer`propFrame[Global`p, Global`l1, Global`cos1, Global`p1, Global`l1];
ntk = NumTrace[net, "Frame" -> frame, "Args" -> {Global`l1, Global`cos1, Global`p, Global`k}];

MakeNTKernel[ntk, FileNameJoin[{$gen, "Glu_quark_dense_kernel.hh"}],
  "Backend" -> "Dense", "Name" -> "Glu_quark_dense", "Namespace" -> "glu_quark_dense"];
(* force the numeric backend's fundamental colour through the NUMERIC SUNNet path *)
NumTracer`Private`$NumTracerFundETMaxEntries = 0;
MakeNTKernel[ntk, FileNameJoin[{$gen, "gen_glu_quark_num.cpp"}],
  FileNameJoin[{$gen, "Glu_quark_num_kernel.hh"}], FileNameJoin[{$gen, "Glu_quark_num_kernels.hh"}],
  "Backend" -> "Numeric", "Name" -> "Glu_quark_num_kernel", "Namespace" -> "glu_quark_num"];
Print["gluon->quark-loop->gluon kernels generated (Dense oracle + numeric)."];
