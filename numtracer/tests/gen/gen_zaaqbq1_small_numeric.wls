#!/usr/bin/env wolframscript
(* QCD generator: the DENSE test flow — the two-gluon-two-quark vertex (ZAAqbq1), cut down.
   Run from numtracer/:   wolfram -script tests/gen/gen_zaaqbq1_small_numeric.wls

   THIS IS THE ONLY FLOW IN THE SUITE THAT EXERCISES THE SUB-TERM DEDUP. The generator contracts one
   Dirac trace per (net, sub-term), and the global dedup in CodegenGenerator.m exists because the SAME trace
   recurs across nets. Every other flow in tests/gen has exactly ONE sub-term per net, so the merge,
   the trace cache and the balanced tree fold are no-ops there — a green suite says nothing about
   them. This flow has hundreds of sub-terms per net and real cross-net redundancy, so it does.

   Why it is dense, and why it is still cheap: see the header of tests/qcd_aaqbq_setup.m. In short,
   the density is the Tmunurho = {[g_mu,g_nu], g_rho} structure of transAAqbqMinimal's element 2,
   which survives the sigma-fold and the dressed-sum split, and which enters the trace twice (as the
   in-loop vertex and inside the projector's dual). The setup restricts the basis to that element,
   which keeps the density while making the Gram 1x1; the real ZAAqbq1 in QCD.nb carries both basis
   elements, four external legs' worth of angles, and ~30 min of generation.

   The guard rails this flow must keep clearing (assert them with tools/redundancy.py after any
   change here) — if it drifts back to ~1 sub-term per net it is silently testing nothing:
       redundancy > 2x        sub-terms/net (max) > 100                                            *)

If[Quiet[Needs["FunKit`"]] === $Failed, (* installed (Applications dir) or the local checkout *)
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];

FSetRegisterSize[64];

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "tests", "qcd_aaqbq_setup.m"}]];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
SetOptions[MakeNTKernel, Sequence @@ $ntShimOpts];  (* tests: emit against the shim + DiFfRG namespace *)

(* 4-point symmetric frame: three angles (cos1, cos2, phi) — the AAqbq vertex has four legs. *)
frame = NumTracer`sp4Frame[Global`p, Global`l1, Global`cos1, Global`cos2, Global`phi,
          Global`p1, Global`p2, Global`p3, Global`p4, Global`l1];
args  = {Global`l1, Global`cos1, Global`cos2, Global`phi, Global`p, Global`k};
dress = {ZA3, ZAcbc, ZA4, ZAqbq1, ZAAqbq1, dtZc, Zc, dtZA, ZA, dtZq, Zq, Mq};

(* Symmetric-point angle relations for four external legs (frame geometry), emitted once in the
   kernel body via "AngleDefs" rather than inlined per occurrence. *)
spAngles4 = {
  Global`cosl1p1 -> Sqrt[1 - Global`cos1^2] Global`cos2,
  Global`cosl1p2 -> -(1/3) Sqrt[1 - Global`cos1^2] (Global`cos2 -
     2 Sqrt[2 - 2 Global`cos2^2] Cos[Global`phi]),
  Global`cosl1p3 -> -(1/3) Sqrt[1 - Global`cos1^2] (Global`cos2 +
     Sqrt[2 - 2 Global`cos2^2] (Cos[Global`phi] - Sqrt[3] Sin[Global`phi])),
  Global`cosl1p4 -> -(1/3) Sqrt[1 - Global`cos1^2] (Global`cos2 +
     Sqrt[2 - 2 Global`cos2^2] (Cos[Global`phi] + Sqrt[3] Sin[Global`phi]))};

(* Full permutation group of the identical external legs; the ClassTrans/transAAqbq
   projectors are covariant under it. Symmetries must be stated by hand since FunKit
   0bacdd1 -- without this the derivative yields every leg ordering as its own diagram
   (measured on za4_147: 6183 diagrams / 164 s, vs 294 / 7.6 s with it). *)
symLegs = FMakeSymmetryList[{A[i1], A[i2], qb[i3], q[i4]}];
fRGAAqbq = FTakeDerivatives[WetterichEquation, {A[i1], A[i2], qb[i3], q[i4]}, "Symmetries" -> symLegs] // FTruncate // FRoute;
proj = FTerm[TBGetProjector["transAAqbqR", 1,
   {i1, i2, i4, i3} /. fRGAAqbq["1-Loop"]["ExternalIndices"]]];
trace = FTerm[proj, fRGAAqbq["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]];

flow = SPParam[dressingRules[trace]];
net  = FromFunKit[flow];
Print["ZAAqbq1(small) net: top terms = ", If[Head[net] === Plus, Length[net], 1]];

(* CUT-DOWN DIAL. NumTrace splits the net's top-level Plus into diagrams, and sub-term density is a
   property of EACH diagram (the Tmunurho x Tmunurho product inside it), not of how many there are.
   So dropping diagrams scales the cost down linearly while leaving sub-terms/net — the thing under
   test — untouched. Sorted first, so the choice is deterministic and not an accident of term order.
   The full 43 diagrams generate in ~5 min and their dedup-OFF control run would be ~24x that; this
   keeps both affordable in ctest. NT_TEST_NDIAG overrides (0 = keep all) for calibration. *)
nKeep = With[{e = Environment["NT_TEST_NDIAG"]},
  If[StringQ[e] && IntegerQ[Quiet@ToExpression[e]], ToExpression[e], 6]];
If[nKeep > 0 && Head[net] === Plus && Length[net] > nKeep,
  net = Plus @@ Take[Sort[List @@ net], nKeep];
  Print["ZAAqbq1(small) cut down to ", nKeep, " diagrams (NT_TEST_NDIAG)"]];

{tA, ntk} = AbsoluteTiming[ NumTrace[net, "Frame" -> frame, "Args" -> args] ];
tG = First @ AbsoluteTiming[ MakeNTKernel[ntk,
  FileNameJoin[{$gen, "gen_zaaqbq1_small_num.cpp"}],
  FileNameJoin[{$gen, "ZAAqbq1Small_num_kernel.hh"}],
  FileNameJoin[{$gen, "ZAAqbq1Small_num_kernels.hh"}],
  "Name" -> "ZAAqbq1Small_num_kernel", "Namespace" -> "zaaqbq1_small_num", "Dressings" -> dress,
  "AngleDefs" -> spAngles4] ];
Print["[time] NumTracer generation TOTAL (analyse+codegen): ", tA + tG, " s"];

(* ---- THE REFERENCE KERNEL: the same flow with the dedup switched OFF -------------------------
   NT_GEN_NO_DEDUP makes the generator contract every (net, sub-term) occurrence separately, merging
   nothing and caching nothing — the computation as it was before the dedup existed. Emitting it here
   as a second kernel gives compare_zaaqbq1_small.cpp a control to grade against, and it is the only
   honest way to grade: a BYTE diff of the two kernels is meaningless, because the dedup changes the
   order traces are lowered and hence GlobalEnv's first-seen interning, which renumbers every sN even
   where the values are identical. Both kernels come from the SAME FunKit derivation and the same
   loop-momentum routing, so unlike the NumTracer-vs-FORM comparisons this one is legitimately
   POINTWISE-exact (no FRoute reparametrisation between them — see FUNKIT_ROUTING_ISSUE.md).

   It is the slow half of this script (no dedup = 9x the contractions here), which is precisely the
   speedup being tested. *)
SetEnvironment["NT_GEN_NO_DEDUP" -> "1"];
tR = First @ AbsoluteTiming[ MakeNTKernel[ntk,
  FileNameJoin[{$gen, "gen_zaaqbq1_small_ref.cpp"}],
  FileNameJoin[{$gen, "ZAAqbq1SmallRef_num_kernel.hh"}],
  FileNameJoin[{$gen, "ZAAqbq1SmallRef_num_kernels.hh"}],
  "Name" -> "ZAAqbq1SmallRef_num_kernel", "Namespace" -> "zaaqbq1_small_ref", "Dressings" -> dress,
  "AngleDefs" -> spAngles4] ];
SetEnvironment["NT_GEN_NO_DEDUP" -> None];
Print["[time] dedup-OFF reference kernel: ", tR, " s  (dedup-ON was ", tG, " s)"];
Print["ZAAqbq1(small) numeric kernel + dedup-off reference generated."];
