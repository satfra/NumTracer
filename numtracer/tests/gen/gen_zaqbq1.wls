#!/usr/bin/env wolframscript
(* QCD generator: emit the ET + DENSE kernels for the quark-gluon vertex (ZAqbq1) flow —
   a 2-loop-angle Dirac vertex. The Dirac trace is taken by the C++ et engine (not FORM):
   the un-traced loop is rewritten into the NumTracer DSL and closed downstream.
   Symmetric point: same SP3 frame as the 3-gluon vertex (any 3-point vertex shares it).
   Run from numtracer/:   wolfram -script tests/gen/gen_zaqbq1.wls                          *)

If[Quiet[Needs["FunKit`"]] === $Failed, (* installed (Applications dir) or the local checkout *)
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
SetDirectory[$repo];                                  (* TBCache/ (incl. AqbqDirect) resolves here *)
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "tests", "qcd_setup.m"}]];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
SetOptions[MakeNTKernel, Sequence @@ $ntShimOpts];  (* tests: emit against the shim + DiFfRG namespace *)

(* SP3 frame: p1,p2,p3 at 120deg in the 0-1 plane, loop at (cos1,cos2). *)
frame = NumTracer`sp3Frame[Global`p, Global`l1, Global`cos1, Global`cos2,
          Global`p1, Global`p2, Global`p3, Global`l1];
args  = {Global`l1, Global`cos1, Global`cos2, Global`p, Global`k};
dress = {ZA3, ZAcbc, ZA4, ZAqbq1, dtZc, Zc, dtZA, ZA, dtZq, Zq, Mq};

(* quark-gluon vertex flow; project onto the classical (AqbqDirect1) tensor structure.
   Field order {A,qb,q} -> projector index order {i1,i3,i2} (gluon, quark, antiquark). *)
fRGAqbq = FTakeDerivatives[WetterichEquation, {A[i1], qb[i2], q[i3]}] // FTruncate // FRoute;
proj = FTerm[TBGetProjector["AqbqDirect1", 1, {i1, i3, i2} /. fRGAqbq["1-Loop"]["ExternalIndices"]]];
trace = proj ** (fRGAqbq["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]);

(* same SP loop-external cosines -> (cos1,cos2) map as the 3-gluon vertex *)
spAngles = {
  Global`cosl1p1 -> Global`cos1,
  Global`cosl1p2 -> (-Global`cos1 + Sqrt[3 - 3 Global`cos1^2] Global`cos2)/2,
  Global`cosl1p3 -> (-Global`cos1 - Sqrt[3 - 3 Global`cos1^2] Global`cos2)/2};
flow = SPParam[trace // dressingRules] /. spAngles;
net  = FromFunKit[flow];
Print["ZAqbq1 net: top terms = ", If[Head[net] === Plus, Length[net], 1]];

ntk = NumTrace[net, "Frame" -> frame, "Args" -> args];
(* NOTE: the symbolic ET kernel for this 2-angle Dirac VERTEX is too large to compile (OOM) —
   the et backend is infeasible here. Only the dense (numeric) and invariant-basis paths are
   generated for ZAqbq1; et remains the feasible path only for the smaller flows (e.g. Zq). *)
tD = First @ AbsoluteTiming[ MakeNTKernel[ntk, FileNameJoin[{$gen, "ZAqbq1_dense_kernel.hh"}], "Name" -> "ZAqbq1_dense_kernel", "Dressings" -> dress,
  "Backend" -> "Dense"] ];
Print["[time] MakeNTKernel (Dense): ", tD, " s"];
Print["ZAqbq1 dense kernel generated."];
