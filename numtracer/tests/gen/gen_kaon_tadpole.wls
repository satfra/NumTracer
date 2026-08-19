#!/usr/bin/env wolframscript
(* REGRESSION FIXTURE — a FUNDAMENTAL flavour index with NO flavour generators anywhere.
   This is the third case: every other flavour flow in tests/gen either carries generators
   (gen_pion_quark: tauFlav + deltaFlavFundGen), spells the engine head itself
   (gen_flavour_ingroup: ntSUNDeltaFund built by hand), is flavour-blind (gen_flavour_split), or
   has a flavour line that CLOSES cheaply (gen_glu_quark's quark loop -> Nf). None of them
   exercises a fundamental flavour delta that the blind contraction cannot finish.

   THE BUG IT PINS. A flavour label rides along a line as a Kronecker delta, and NumTracer has two
   ways to handle it: the cheap rewrite rules (stitch a chain, close a loop to Nf) and the SU(N)
   engine. It used to pick the engine only when it also saw flavour GENERATORS. A kaon-like field
   (K_i ~ sbar g5 q_i) carries a fundamental flavour index but couples through no generator, so it
   got the cheap path — and in the TADPOLE, where the quartic contributes
   delta_ac delta_bd + delta_ad delta_bc against the loop's own delta, the deltas form a WEB rather
   than a chain. No local rewrite closes a web, so the residue was CForm'd into the kernel as
   `NumTracer_Private_flavDelta(F1, F2)`, an undeclared identifier that GCC and Clang both parse.

   THE ORACLE. TFlav[0,f1,f2] = deltaFundFlav[f1,f2]/Sqrt[2 Nf], so the two spellings below are the
   SAME object: `fund` writes the delta plainly, `tflav` writes it as Sqrt[2 Nf] TFlav[0,..], which
   always worked because it trips the generator gate. compare_kaon_tadpole.cpp requires the two
   kernels to agree pointwise. That is an exact self-contained equivalence — no external truth, no
   dense backend, no FORM oracle.

   The TADPOLE is the whole point: the self-energy topology has the boson's flavour line EXTERNAL
   (it closes trivially), and only the quartic's internal loop produces the web.

   Run from numtracer/:   wolfram -script tests/gen/gen_kaon_tadpole.wls *)

If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];

(* shim opt-in (as gen_flavour_split / gen_flavour_ingroup): dressings become DiFfRG::Fn runtime
   params so the compare driver can feed both kernels the same functions. *)
$ntShimOpts = {"RuntimeInclude" -> "shim.hpp", "ExtraIncludes" -> {},
  "KernelNamespace" -> "DiFfRG", "SupportNamespace" -> "DiFfRG",
  "DressingType" -> "SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>"};
SetOptions[MakeNTKernel, Sequence @@ $ntShimOpts];

SetNc[3]; SetNf[2];   (* the light doublet the kaon's index runs over *)

(* One complex boson pair carrying ONE fundamental flavour index. No colour, no Dirac, no
   TensorBases basis: the flavour sector is the entire content of this flow. *)
fields     = <|"Commuting" -> {{Kb[p, {F}], K[p, {F}]}}|>;
truncation = <|GammaN -> {{Kb, K}, {Kb, K, Kb, K}}, Propagator -> {{Kb, K}},
               R -> {{Kb, K}}, Rdot -> {{Kb, K}}, Field -> {{}}|>;
FSetGlobalSetup[<|"FieldSpace" -> fields, "Truncation" -> truncation,
  "FeynmanRules" -> <|GammaN -> {}, Propagator -> {}, R -> {}, Rdot -> {}|>|>];

fRG  = FTakeDerivatives[WetterichEquation, {Kb[i1], K[i2]}] // FTruncate // FRoute;
expr = fRG["1-Loop"]["Expression"];

(* Hand-rolled rules, parameterised by HOW the flavour delta is spelled. Everything else — momenta,
   dressings, the quartic's two pairings — is identical between the two variants, so any difference
   in the emitted kernels is a difference in the flavour routing and nothing else. *)
rulesFor[d_] := {
  Propagator[{K, Kb} | {Kb, K}, {{q1_, {a_}}, {q2_, {b_}}}] :>
    d[a, b] * dressing[GammaN, {Kb, K}, 1, {q1, q2}],
  Rdot[{Kb, K}, {{q1_, {a_}}, {q2_, {b_}}}] :>
    d[a, b] * dressing[Rdot, {Kb, K}, 1, {q1, q2}],
  (* FRoute can offer a 2-leg pattern for the 4-point vertex; only the 4-leg one is meant. *)
  GammaN[{Kb, Kb, K, K}, {{q1_, {a_}}, {q2_, {b_}}}] :> 0,
  (* The quartic from V(rho), rho = Kbar.K: the two pairings of the Kb indices with the K indices.
     Kb outranks K, so the canonical leg order is {Kb,Kb,K,K}. THIS is the delta web. *)
  GammaN[{Kb, Kb, K, K}, {{q1_, {a_}}, {q2_, {b_}}, {q3_, {c_}}, {q4_, {d4_}}}] :>
    d[a, c] d[b, d4] + d[a, d4] d[b, c]};

(* External projector: close the two external flavour legs and normalise by Nf. *)
buildNet[d_] := Module[{tr},
  tr = FTerm[d[F1, F2]/Global`Nf, (expr /. rulesFor[d])] //. p2 -> -p1;
  tr = (Times @@ # &) /@ List @@ tr;
  tr = FEx @@ FTerm /@ tr;
  tr = tr //. {
    dressing[GammaN, {Kb, K}, 1, {m1_, m2_}] :> GK[Sqrt[sp[m2, m2]]],
    dressing[Rdot,   {Kb, K}, 1, {m1_, m2_}] :> GKdot[Sqrt[sp[m2, m2]]]};
  FromFunKit[UseLorentzLinearity[tr] //. {lf1 -> l1}]];

frame = NumTracer`propFrame[Global`p, Global`l1, Global`cos1, Global`p1, Global`l1];
emit[tag_, d_] := Module[{ntk},
  ntk = NumTrace[buildNet[d], "Frame" -> frame,
                 "Args" -> {Global`l1, Global`cos1, Global`p}];
  MakeNTKernel[ntk, FileNameJoin[{$gen, "gen_kaon_tadpole_" <> tag <> ".cpp"}],
    FileNameJoin[{$gen, "Kaon_tadpole_" <> tag <> "_kernel.hh"}],
    FileNameJoin[{$gen, "Kaon_tadpole_" <> tag <> "_kernels.hh"}],
    "Name" -> "Kaon_tadpole_" <> tag <> "_kernel", "Namespace" -> "kaon_tadpole_" <> tag,
    "Dressings" -> {GK, GKdot}]];

(* (1) the plain fundamental delta — the spelling that used to leak. *)
emit["fund", Global`deltaFundFlav];
(* (2) the SAME object through the generator route, which always worked: the oracle. *)
emit["tflav", (Sqrt[2 Global`Nf] Global`TFlav[0, #1, #2] &)];

(* Belt and braces: the emitted headers must not mention the blind head under any spelling. The
   textual $ntCppLeakPatterns guard in Codegen.m would already have aborted, but this fixture exists
   precisely for this class, so it says so itself rather than trusting a guard elsewhere. *)
Do[With[{txt = Import[FileNameJoin[{$gen, "Kaon_tadpole_" <> t <> "_kernel.hh"}], "Text"]},
    If[StringContainsQ[txt, "flavDelta"] || StringContainsQ[txt, "_Private_"],
      Print["[kaon_tadpole] FAIL: an uncontracted flavour delta reached ", t, "'s kernel."];
      Exit[1]]],
  {t, {"fund", "tflav"}}];

Print["kaon-tadpole kernels generated (fund + tflav oracle); no flavour delta survived."];
