#!/usr/bin/env wolframscript
(* Pure-YM ZAcbc (ghost-gluon vertex) via the NUMERIC matrix-product backend, plus the Dense oracle.
   This vertex has a diagram with a DISCONNECTED Lorentz network (two index-disjoint trace factors that
   MULTIPLY) — the regression guard for the "components summed instead of multiplied" codegen bug.
   Run from numtracer/:  wolfram -script tests/gen/gen_ym_zacbc_numeric.wls                          *)

If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];
FSetRegisterSize[64];

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "tests", "ym_setup.m"}]];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
SetOptions[MakeNTKernel, Sequence @@ $ntShimOpts];

frame = NumTracer`sp3Frame[Global`p, Global`l1, Global`cos1, Global`cos2,
          Global`p1, Global`p2, Global`p3, Global`l1];
args  = {Global`l1, Global`cos1, Global`cos2, Global`p, Global`k};
dress = {ZA3, ZAcbc, ZA4, dtZc, Zc, dtZA, ZA};
spAngles = {
  Global`cosl1p1 -> Global`cos1,
  Global`cosl1p2 -> (-Global`cos1 + Sqrt[3 - 3 Global`cos1^2] Global`cos2)/2,
  Global`cosl1p3 -> (-Global`cos1 - Sqrt[3 - 3 Global`cos1^2] Global`cos2)/2};

fRG   = FTakeDerivatives[WetterichEquation, {A[i1], cb[i2], c[i3]}] // FTruncate // FPlot // FRoute // FPrint;
proj  = FTerm[TBGetProjector["Acbc", 1, {i1, i2, i3} /. fRG["1-Loop"]["ExternalIndices"]]];
trace = proj ** (fRG["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]);
flow  = SPParam[trace // dressingRules];
net   = FromFunKit[flow];
Print["ZAcbc net top terms = ", If[Head[net] === Plus, Length[net], 1]];

{tA, ntk} = AbsoluteTiming[ NumTrace[net, "Frame" -> frame, "Args" -> args] ];
Print["[time] NumTrace (analyse): ", tA, " s   diagrams = ", Length[ntk[[1]]["Diagrams"]]];

tN = First @ AbsoluteTiming[ MakeNTKernel[ntk,
  FileNameJoin[{$gen, "gen_ym_zacbc_num.cpp"}],
  FileNameJoin[{$gen, "ZAcbc_num_kernel.hh"}],
  FileNameJoin[{$gen, "ZAcbc_num_kernels.hh"}],
  "Backend" -> "Numeric", "Name" -> "ZAcbc_num_kernel", "Namespace" -> "zacbc_num",
  "Dressings" -> dress, "AngleDefs" -> spAngles] ];
Print["[time] Numeric backend (contract+compile+run): ", tN, " s"];

tD = First @ AbsoluteTiming[ MakeNTKernel[ntk, FileNameJoin[{$gen, "ZAcbc_num_dense_kernel.hh"}],
  "Name" -> "ZAcbc_num_dense_kernel", "Dressings" -> dress, "Backend" -> "Dense", "AngleDefs" -> spAngles] ];
Print["[time] Dense oracle: ", tD, " s"];
Print["ZAcbc numeric kernels generated."];
