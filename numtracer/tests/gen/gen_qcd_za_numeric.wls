#!/usr/bin/env wolframscript
(* QCD ZA (gluon self-energy, with quark loop) via the NUMERIC matrix-product backend (task #22).
   Emits the Numeric kernel + a Dense oracle for the comparison test. Run from numtracer/:
     wolfram -script tests/gen/gen_qcd_za_numeric.wls                                             *)

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

frame = NumTracer`propFrame[Global`p, Global`l1, Global`cos1, Global`p1, Global`l1];
args  = {Global`l1, Global`cos1, Global`p, Global`k};
dress = {ZA3, ZAcbc, ZA4, ZAqbq1, dtZc, Zc, dtZA, ZA, dtZq, Zq, Mq};

fRGAA = FTakeDerivatives[WetterichEquation, {A[i1], A[i2]}] // FTruncate // FPlot // FRoute // FPrint;
trace = FTerm[1/sp[p1, p1], TBGetProjector["AA", 1, {i1, i2} /. fRGAA["1-Loop"]["ExternalIndices"]]] **
        (fRGAA["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]) //. p2 -> -p1;
flow = PropParam[trace // dressingRules];
net  = FromFunKit[flow];
Print["ZA net top terms = ", If[Head[net] === Plus, Length[net], 1]];

ntk = NumTrace[net, "Frame" -> frame, "Args" -> args];

(* NUMERIC backend (frame fallback -> Components polynomialised; sin1 = sqrt(1-cos1^2)) *)
tN = First @ AbsoluteTiming[ MakeNTKernel[ntk,
  FileNameJoin[{$gen, "gen_qcd_za_num.cpp"}],
  FileNameJoin[{$gen, "ZA_num_kernel.hh"}],
  FileNameJoin[{$gen, "ZA_num_kernels.hh"}],
  "Backend" -> "Numeric", "Name" -> "ZA_num_kernel", "Namespace" -> "za_num", "Dressings" -> dress] ];
Print["[time] Numeric backend (contract+compile+run): ", tN, " s"];

(* Dense oracle for the comparison *)
tD = First @ AbsoluteTiming[ MakeNTKernel[ntk, FileNameJoin[{$gen, "ZA_num_dense_kernel.hh"}],
  "Name" -> "ZA_num_dense_kernel", "Dressings" -> dress, "Backend" -> "Dense"] ];
Print["[time] Dense oracle: ", tD, " s"];
Print["ZA numeric kernels generated."];
