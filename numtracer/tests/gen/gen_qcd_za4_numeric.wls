#!/usr/bin/env wolframscript
(* QCD ZA4 (four-gluon vertex, quark box) via the NUMERIC matrix-product backend — rank-4 stress
   test (4 free legs, sp4Frame with phi). Emits the Numeric kernel; validate vs the FORM oracle. *)
If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"]; Get["DiFfRG`"]; FSetRegisterSize[64];
$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "tests", "qcd_setup.m"}]];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
SetOptions[MakeNTKernel, Sequence @@ $ntShimOpts];

frame = NumTracer`sp4Frame[Global`p, Global`l1, Global`cos1, Global`cos2, Global`phi,
          Global`p1, Global`p2, Global`p3, Global`p4, Global`l1];
args  = {Global`l1, Global`cos1, Global`cos2, Global`phi, Global`p, Global`k};
dress = {ZA3, ZAcbc, ZA4, ZAqbq1, dtZc, Zc, dtZA, ZA, dtZq, Zq, Mq};
spAngles4 = {
  Global`cosl1p1 -> Sqrt[1 - Global`cos1^2] Global`cos2,
  Global`cosl1p2 -> -(1/3) Sqrt[1 - Global`cos1^2] (Global`cos2 - 2 Sqrt[2 - 2 Global`cos2^2] Cos[Global`phi]),
  Global`cosl1p3 -> -(1/3) Sqrt[1 - Global`cos1^2] (Global`cos2 + Sqrt[2 - 2 Global`cos2^2] (Cos[Global`phi] - Sqrt[3] Sin[Global`phi])),
  Global`cosl1p4 -> -(1/3) Sqrt[1 - Global`cos1^2] (Global`cos2 + Sqrt[2 - 2 Global`cos2^2] (Cos[Global`phi] + Sqrt[3] Sin[Global`phi]))};

fRG  = FTakeDerivatives[WetterichEquation, {A[i1], A[i2], A[i3], A[i4]}] // FTruncate // FPlot // FRoute // FPrint;
proj = FTerm[TBGetProjector["AAAAClassTrans", 1, {i1, i2, i3, i4} /. fRG["1-Loop"]["ExternalIndices"]]];
traceA4 = proj ** (fRG["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]);
net = FromFunKit[ SPParam[traceA4 // dressingRules] ];
Print["ZA4 net top terms = ", If[Head[net] === Plus, Length[net], 1]];

{tA, ntk} = AbsoluteTiming[ NumTrace[net, "Frame" -> frame, "Args" -> args] ];
Print["[time] NumTrace (analyse): ", tA, " s   diagrams = ", Length[ntk[[1]]["Diagrams"]]];

tN = First @ AbsoluteTiming[ MakeNTKernel[ntk,
  FileNameJoin[{$gen, "gen_qcd_za4_num.cpp"}], FileNameJoin[{$gen, "ZA4_num_kernel.hh"}], FileNameJoin[{$gen, "ZA4_num_kernels.hh"}],
  "Backend" -> "Numeric", "Name" -> "ZA4_num_kernel", "Namespace" -> "za4_num", "Dressings" -> dress, "AngleDefs" -> spAngles4] ];
Print["[time] Numeric backend: ", tN, " s"];
Print["ZA4 numeric kernel generated."];
