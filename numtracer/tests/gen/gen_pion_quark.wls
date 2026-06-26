#!/usr/bin/env wolframscript
(* Phase-2 example 2 — the FUNDAMENTAL-flavour + Dirac + fundamental-colour sectors at once.
   A pion (isospin-triplet) self-energy from a QUARK LOOP via the Yukawa coupling
       h * qbar * i gamma_5 * (tau^a pi^a) * q .
   The one-loop graph is a closed quark line with two Yukawa vertices and a regulator insertion;
   tracing it produces, in one shot:
     - a Dirac trace tr[gamma_5 (q-slash) gamma_5 (q-slash) ...]      (the inv backend's gamma-trace theorem),
     - a FUNDAMENTAL-flavour trace tr(tau^{i1} tau^{i2}) = (1/2) delta^{i1 i2}  (SU(2) generators, the 2nd group),
     - a FUNDAMENTAL-colour trace tr(1) = Nc.
   The quark-propagator flavour delta sits INSIDE the tau-trace, so it is emitted as deltaFlavFundGen
   (-> ntSUNDeltaFund[Nf, ...], kept in the engine) rather than the trivial deltaFundFlav (-> Nf).

   Method: derive with FunKit, HAND-ROLL the diagrammatic rules. Massless quark with a q-slash
   propagator numerator and a fermion regulator R_F ~ q-slash (so the Rdot insertion carries the
   even 4th gamma — a massless 3-propagator loop alone would have an odd gamma count and vanish).

   Run from numtracer/:   wolfram -script tests/gen/gen_pion_quark.wls *)

If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];

fields = <|"Commuting" -> {Pi[p, {i}]}, "Grassmann" -> {{qb[p, {d, A, F}], q[p, {d, A, F}]}}|>;
truncation = <|GammaN -> {{Pi, Pi}, {Pi, qb, q}, {qb, q}}, Propagator -> {{Pi, Pi}, {qb, q}},
  R -> {{qb, q}}, Rdot -> {{qb, q}}, Field -> {{}}|>;
FSetGlobalSetup[<|"FieldSpace" -> fields, "Truncation" -> truncation,
  "FeynmanRules" -> <|GammaN -> {}, Propagator -> {}, R -> {}, Rdot -> {}|>|>];

fRG = FTakeDerivatives[WetterichEquation, {Pi[i1], Pi[i2]}] // FTruncate // FRoute;
expr = fRG["1-Loop"]["Expression"];

(* ---- hand-rolled rules. Yukawa: i gamma5 (Dirac) x tau^a (fundamental flavour) x delta (colour).
        Quark prop / Rdot: massless q-slash numerator x delta_color x delta_flav(in-trace). ---- *)
myRules = {
  GammaN[{Pi, qb, q}, {{_, {iso_}}, {_, {dqb_, cqb_, fqb_}}, {_, {dq_, cq_, fq_}}}] :>
    I gamma5[dqb, dq] tauFlav[iso, fqb, fq] deltaFundCol[cqb, cq],
  Propagator[{q, qb} | {qb, q}, {{m1_, {d1_, c1_, f1_}}, {m2_, {d2_, c2_, f2_}}}] :>
    With[{mu = Unique["lm"]},
      deltaFundCol[c1, c2] deltaFlavFundGen[f1, f2] *
      (-I) gamma[mu, d1, d2] vec[m2, mu] * dressing[GammaN, {qb, q}, 1, {m1, m2}]],
  Rdot[{qb, q}, {{m1_, {d1_, c1_, f1_}}, {m2_, {d2_, c2_, f2_}}}] :>
    With[{mu = Unique["lm"]},
      deltaFundCol[c1, c2] deltaFlavFundGen[f1, f2] *
      (-I) gamma[mu, d1, d2] vec[m2, mu] * dressing[Rdot, {qb, q}, 1, {m1, m2}]]};

(* project the external isospin (i1,i2) onto the scalar self-energy *)
trace = FTerm[1/3, deltaAdjFlav[Global`i1, Global`i2]] ** (expr /. myRules);

(* baked closed-form fermion dressings: G(q^2) = 1/(q^2 + RF(k^2,q^2)) ; regulator-dot RFdot. *)
dressed = trace //. {
  dressing[GammaN, {qb, q}, 1, {m1_, m2_}] :> 1/(sp[m2, m2] + RF[k^2, sp[m2, m2]]),
  dressing[Rdot, {qb, q}, 1, {m1_, m2_}] :> RFdot[k^2, sp[m2, m2]]};
dressed = UseLorentzLinearity[dressed] //. {lf1 -> l1};

(* SetNc[3] before FromFunKit: the colour rank Nc is baked into the SU(N) heads at import time,
   and the isospin rank is passed via FromFunKit's "FlavourGroup" option. *)
SetNc[3];
net = FromFunKit[dressed, "FlavourGroup" -> 2];

frame = NumTracer`propFrame[Global`p, Global`l1, Global`cos1, Global`p1, Global`l1];
ntk = NumTrace[net, "Frame" -> frame, "Args" -> {Global`l1, Global`cos1, Global`p, Global`k}];

MakeNTKernel[ntk, FileNameJoin[{$gen, "Pion_quark_dense_kernel.hh"}],
  "Backend" -> "Dense", "Name" -> "Pion_quark_dense", "Namespace" -> "pion_quark_dense"];
MakeNTKernel[ntk, FileNameJoin[{$gen, "gen_pion_quark_num.cpp"}],
  FileNameJoin[{$gen, "Pion_quark_num_kernel.hh"}], FileNameJoin[{$gen, "Pion_quark_num_kernels.hh"}],
  "Backend" -> "Numeric", "Name" -> "Pion_quark_num_kernel", "Namespace" -> "pion_quark_num"];
Print["pion->quark-loop->pion Yukawa kernels generated (Dense oracle + numeric)."];
