#!/usr/bin/env wolframscript
(* Stage-3 generator: emit the DIRECT ET kernel for the quark wave-function (Zq) flow.
   The Dirac trace is NOT taken in FORM — the un-traced fermion loop is rewritten into the
   NumTracer DSL (gamma -> ntGamma -> Dirac::gamma_axis, etc.) and the trace is closed by
   the C++ et engine. Run from numtracer/:   wolfram -script tests/gen/gen_zq.wls           *)

If[Quiet[Needs["FunKit`"]] === $Failed, (* installed (Applications dir) or the local checkout *)
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
SetDirectory[$repo];                                  (* TBCache/ (incl. qbq) resolves here *)
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "tests", "qcd_setup.m"}]];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
SetOptions[MakeNTKernel, Sequence @@ $ntShimOpts];  (* tests: emit against the shim + DiFfRG namespace *)

(* symmetric-point frame, 1 angle: p1 along axis 0, loop l1 at angle cos1 *)
frame = NumTracer`propFrame[Global`p, Global`l1, Global`cos1, Global`p1, Global`l1];
args  = {Global`l1, Global`cos1, Global`p, Global`k};
dress = {ZA, Zq, Mq, ZAqbq1, dtZA, dtZq};

fRGqbq = FTakeDerivatives[WetterichEquation, {qb[i1], q[i2]}] // FTruncate // FRoute;
trace = FTerm[-TBGetProjector["qbq", 1, {i2, i1} /. fRGqbq["1-Loop"]["ExternalIndices"]]] **
        (fRGqbq["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]);
flow = PropParam[trace // dressingRules];
net  = FromFunKit[flow];

ntk = NumTrace[net, "Frame" -> frame, "Args" -> args];
tE = First @ AbsoluteTiming[ MakeNTKernel[ntk, FileNameJoin[{$gen, "Zq_gen_kernel.hh"}], "Name" -> "Zq_gen_kernel", "Dressings" -> dress] ];
Print["[time] MakeNTKernel (ET): ", tE, " s"];
Print["Zq ET kernel generated."];
(* dense brute-force baseline: same network, numeric DTensor backend (numtracer::dense::). *)
tD = First @ AbsoluteTiming[ MakeNTKernel[ntk, FileNameJoin[{$gen, "Zq_dense_kernel.hh"}], "Name" -> "Zq_dense_kernel", "Dressings" -> dress,
  "Backend" -> "Dense"] ];
Print["[time] MakeNTKernel (Dense): ", tD, " s"];
Print["Zq dense kernel generated."];
