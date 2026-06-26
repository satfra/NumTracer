#!/usr/bin/env wolframscript
(* QCD ZA3 (three-gluon vertex, with quark triangle) via the NUMERIC matrix-product backend (task
   #22) — the headline target. Emits the Numeric kernel + a Dense oracle. Run from numtracer/:
     wolfram -script tests/gen/gen_qcd_za3_numeric.wls                                            *)

If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];
FSetRegisterSize[64];

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "tests", "qcd_setup.m"}]];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
SetOptions[MakeNTKernel, Sequence @@ $ntShimOpts];

frame = NumTracer`sp3Frame[Global`p, Global`l1, Global`cos1, Global`cos2,
          Global`p1, Global`p2, Global`p3, Global`l1];
args  = {Global`l1, Global`cos1, Global`cos2, Global`p, Global`k};
dress = {ZA3, ZAcbc, ZA4, ZAqbq1, dtZc, Zc, dtZA, ZA, dtZq, Zq, Mq};

fRGA3 = FTakeDerivatives[WetterichEquation, {A[i1], A[i2], A[i3]}] // FTruncate // FPlot // FRoute // FPrint;
projA3 = FTerm[TBGetProjector["AAAClassTrans", 1, {i1, i2, i3} /. fRGA3["1-Loop"]["ExternalIndices"]]];
traceA3 = projA3 ** (fRGA3["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]);
spAngles = {
  Global`cosl1p1 -> Global`cos1,
  Global`cosl1p2 -> (-Global`cos1 + Sqrt[3 - 3 Global`cos1^2] Global`cos2)/2,
  Global`cosl1p3 -> (-Global`cos1 - Sqrt[3 - 3 Global`cos1^2] Global`cos2)/2};
flow = SPParam[traceA3 // dressingRules];
net  = FromFunKit[flow];
Print["ZA3 net top terms = ", If[Head[net] === Plus, Length[net], 1]];

{tA, ntk} = AbsoluteTiming[ NumTrace[net, "Frame" -> frame, "Args" -> args] ];
Print["[time] NumTrace (analyse): ", tA, " s   diagrams = ", Length[ntk[[1]]["Diagrams"]]];

tN = First @ AbsoluteTiming[ MakeNTKernel[ntk,
  FileNameJoin[{$gen, "gen_qcd_za3_num.cpp"}],
  FileNameJoin[{$gen, "ZA3_num_kernel.hh"}],
  FileNameJoin[{$gen, "ZA3_num_kernels.hh"}],
  "Backend" -> "Numeric", "Name" -> "ZA3_num_kernel", "Namespace" -> "za3_num",
  "Dressings" -> dress, "AngleDefs" -> spAngles] ];
Print["[time] Numeric backend (contract+compile+run): ", tN, " s"];

tD = First @ AbsoluteTiming[ MakeNTKernel[ntk, FileNameJoin[{$gen, "ZA3_num_dense_kernel.hh"}],
  "Name" -> "ZA3_num_dense_kernel", "Dressings" -> dress, "Backend" -> "Dense", "AngleDefs" -> spAngles] ];
Print["[time] Dense oracle: ", tD, " s"];
Print["ZA3 numeric kernels generated."];
