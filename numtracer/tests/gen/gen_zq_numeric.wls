#!/usr/bin/env wolframscript
(* Zq (quark wave-function self-energy, fermion flow) via the NUMERIC matrix-product backend. *)
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

frame = NumTracer`propFrame[Global`p, Global`l1, Global`cos1, Global`p1, Global`l1];
args  = {Global`l1, Global`cos1, Global`p, Global`k};
dress = {ZA, Zq, Mq, ZAqbq1, dtZA, dtZq};
fRGqbq = FTakeDerivatives[WetterichEquation, {qb[i1], q[i2]}] // FTruncate // FRoute;
trace = FTerm[-TBGetProjector["qbq", 1, {i2, i1} /. fRGqbq["1-Loop"]["ExternalIndices"]]] **
        (fRGqbq["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]);
flow = PropParam[trace // dressingRules];
net  = FromFunKit[flow];
ntk  = NumTrace[net, "Frame" -> frame, "Args" -> args];

tN = First @ AbsoluteTiming[ MakeNTKernel[ntk,
  FileNameJoin[{$gen, "gen_zq_num.cpp"}], FileNameJoin[{$gen, "Zq_num_kernel.hh"}], FileNameJoin[{$gen, "Zq_num_kernels.hh"}],
  "Backend" -> "Numeric", "Name" -> "Zq_num_kernel", "Namespace" -> "zq_num", "Dressings" -> dress] ];
Print["[time] Numeric backend: ", tN, " s"];
tD = First @ AbsoluteTiming[ MakeNTKernel[ntk, FileNameJoin[{$gen, "Zq_num_dense_kernel.hh"}],
  "Name" -> "Zq_num_dense_kernel", "Dressings" -> dress, "Backend" -> "Dense"] ];
Print["[time] Dense oracle: ", tD, " s"];
Print["Zq numeric kernels generated."];
