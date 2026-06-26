#!/usr/bin/env wolframscript
(* QCD generator: emit the NUMERIC matrix-product kernel for the quark-gluon vertex (ZAqbq1) flow.
   The Dirac trace is taken by the gamma-trace THEOREM (NumTracer, not FORM) -> a Lorentz network
   contracted component-numerically over the loop frame -> committed straight-line polynomials.
   Run from numtracer/:   wolfram -script tests/gen/gen_zaqbq1_numeric.wls                       *)

If[Quiet[Needs["FunKit`"]] === $Failed, (* installed (Applications dir) or the local checkout *)
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];

(* Match the FormTracer oracle's COEN budget (QCD_vacuum_base/QCD.m: FSetRegisterSize[64]):
   modest per-sub-kernel register pressure (GPU-friendly), NOT a mass-hoist. The repeated
   sqrt(3-3 cos1^2) that used to flood the integrand is removed at the source — we keep the
   loop/external angles cosl1p1/2/3 SYMBOLIC (like the oracle's SPParam) and define them once
   in the kernel body, instead of expanding them via spAngles into 86 inline sqrt. *)
FSetRegisterSize[64];

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "tests", "qcd_setup.m"}]];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
SetOptions[MakeNTKernel, Sequence @@ $ntShimOpts];  (* tests: emit against the shim + DiFfRG namespace *)

frame = NumTracer`sp3Frame[Global`p, Global`l1, Global`cos1, Global`cos2,
          Global`p1, Global`p2, Global`p3, Global`l1];
args  = {Global`l1, Global`cos1, Global`cos2, Global`p, Global`k};
dress = {ZA3, ZAcbc, ZA4, ZAqbq1, dtZc, Zc, dtZA, ZA, dtZq, Zq, Mq};

fRGAqbq = FTakeDerivatives[WetterichEquation, {A[i1], qb[i2], q[i3]}] // FTruncate // FRoute;
proj = FTerm[TBGetProjector["AqbqDirect1", 1, {i1, i3, i2} /. fRGAqbq["1-Loop"]["ExternalIndices"]]];
trace = proj ** (fRGAqbq["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]);

(* Symmetric-point angle relations (frame geometry). The dressing keeps cosl1p1/2/3 SYMBOLIC
   (SPParam in qcd_setup.m maps cos[l1,pj] :> cosl1pj); these defs are emitted ONCE in the
   kernel body via "AngleDefs" so the shared sqrt(3-3 cos1^2) is computed once, not inlined
   per occurrence — exactly the FormTracer oracle's KernelBody/SP3Defs approach. *)
spAngles = {
  Global`cosl1p1 -> Global`cos1,
  Global`cosl1p2 -> (-Global`cos1 + Sqrt[3 - 3 Global`cos1^2] Global`cos2)/2,
  Global`cosl1p3 -> (-Global`cos1 - Sqrt[3 - 3 Global`cos1^2] Global`cos2)/2};
flow = SPParam[trace // dressingRules];
net  = FromFunKit[flow];
Print["ZAqbq1 net: top terms = ", If[Head[net] === Plus, Length[net], 1]];

{tA, ntk} = AbsoluteTiming[ NumTrace[net, "Frame" -> frame, "Args" -> args] ];
tG = First @ AbsoluteTiming[ MakeNTKernel[ntk,
  FileNameJoin[{$gen, "gen_zaqbq1_num.cpp"}],
  FileNameJoin[{$gen, "ZAqbq1_num_kernel.hh"}],
  FileNameJoin[{$gen, "ZAqbq1_num_kernels.hh"}],
  "Backend" -> "Numeric", "Name" -> "ZAqbq1_num_kernel", "Namespace" -> "zaqbq1_num", "Dressings" -> dress,
  "AngleDefs" -> spAngles] ];
Print["[time] NumTracer generation TOTAL (analyse+codegen): ", tA + tG, " s"];
Print["ZAqbq1 numeric kernel generated."];
