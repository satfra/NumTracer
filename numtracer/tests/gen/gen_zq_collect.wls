#!/usr/bin/env wolframscript
(* Zq via the NUMERIC backend with SYMBOLIC DRESSING COLLECTION (DressingCollection -> True): the
   quark propagator numerator Mq·δ + Zq·γ·l is kept EAGER (one ntDressedNum slot) and folds to ONE
   DPoly trace instead of distributing into separate diagrams. Validation harness for the dressed
   path — the emitted kernel must match the FORM reference exactly like the distributed Zq does. *)
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
(* DressingCollection -> True on BOTH FromFunKit (its expandBridges runs first) and NumTrace *)
net  = FromFunKit[flow, "DressingCollection" -> True];
ntk  = NumTrace[net, "Frame" -> frame, "Args" -> args, "DressingCollection" -> True];

(* distinct namespace zq_collect so the collected kernel can sit beside the committed (distributed)
   zq_num kernel in one compare TU (compare_zq_collect.cpp asserts they agree). *)
tN = First @ AbsoluteTiming[ MakeNTKernel[ntk,
  FileNameJoin[{$gen, "gen_zq_collect.cpp"}], FileNameJoin[{$gen, "Zq_collect_kernel.hh"}], FileNameJoin[{$gen, "Zq_collect_kernels.hh"}],
  "Backend" -> "Numeric", "Name" -> "Zq_collect_kernel", "Namespace" -> "zq_collect", "Dressings" -> dress] ];
Print["[time] Numeric (collected) backend: ", tN, " s"];
Print["Zq collected kernel generated."];
