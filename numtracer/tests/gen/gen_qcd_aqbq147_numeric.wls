#!/usr/bin/env wolframscript
(* QCD quark-gluon vertex flow (A qb q 3-point) with the FULL {1,4,7} vertex tensor basis, projected
   ONTO EACH external tensor structure (1, 4, 7) in turn, via the NUMERIC matrix-product backend.
   Sibling of gen_zaqbq1_inv.wls (struct-1 only): loads qcd147_setup.m so the LOOP vertices use
   AqbqDirect147, and emits one kernel per projection
     ZAqbq1_147_num_kernel / ZAqbq4_147_num_kernel / ZAqbq7_147_num_kernel.
   2 loop angles (sp3Frame, cos1/cos2). Run from numtracer/:
     wolfram -script tests/gen/gen_qcd_aqbq147_numeric.wls                                          *)

If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];
FSetRegisterSize[64];

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "tests", "qcd147_setup.m"}]];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
SetOptions[MakeNTKernel, Sequence @@ $ntShimOpts];
$RecursionLimit = 1000000;  (* the dense 1/4/7 codegen tree exceeds the default 1024 *)

frame = NumTracer`sp3Frame[Global`p, Global`l1, Global`cos1, Global`cos2,
          Global`p1, Global`p2, Global`p3, Global`l1];
args  = {Global`l1, Global`cos1, Global`cos2, Global`p, Global`k};
dress = {ZA3, ZAcbc, ZA4, ZAqbq1, ZAqbq4, ZAqbq7, dtZc, Zc, dtZA, ZA, dtZq, Zq, Mq};
spAngles = {
  Global`cosl1p1 -> Global`cos1,
  Global`cosl1p2 -> (-Global`cos1 + Sqrt[3 - 3 Global`cos1^2] Global`cos2)/2,
  Global`cosl1p3 -> (-Global`cos1 - Sqrt[3 - 3 Global`cos1^2] Global`cos2)/2};

(* shared flow: the A,qb,q 3-point 1-loop, before projection *)
fRGAqbq = FTakeDerivatives[WetterichEquation, {A[i1], qb[i2], q[i3]}] // FTruncate // FRoute;
ext  = fRGAqbq["1-Loop"]["ExternalIndices"];
expr = fRGAqbq["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[];

(* basis position k -> physical structure label: AqbqDirect147 = TBRestrictBasis[...,{1,4,7}] *)
projIdx = {1, 2, 3};
labels  = {1, 4, 7};
Do[
  Module[{k = projIdx[[s]], lab = labels[[s]], proj, trace, flow, net, ntk, tA, tN},
    proj  = FTerm[TBGetProjector["AqbqDirect147", k, {i1, i3, i2} /. ext]];
    trace = proj ** expr;
    flow  = SPParam[trace // dressingRules];
    net   = FromFunKit[flow];
    Print["ZAqbq", lab, " (147) net top terms = ", If[Head[net] === Plus, Length[net], 1]];
    {tA, ntk} = AbsoluteTiming[ NumTrace[net, "Frame" -> frame, "Args" -> args] ];
    Print["[time] NumTrace ZAqbq", lab, ": ", tA, " s   diagrams = ", Length[ntk[[1]]["Diagrams"]]];
    tN = First @ AbsoluteTiming[ MakeNTKernel[ntk,
      FileNameJoin[{$gen, "gen_qcd_aqbq147_" <> ToString[lab] <> "_num.cpp"}],
      FileNameJoin[{$gen, "ZAqbq" <> ToString[lab] <> "_147_num_kernel.hh"}],
      FileNameJoin[{$gen, "ZAqbq" <> ToString[lab] <> "_147_num_kernels.hh"}],
      "Backend" -> "Numeric", "Name" -> "ZAqbq" <> ToString[lab] <> "_147_num_kernel",
      "Namespace" -> "zaqbq" <> ToString[lab] <> "_147_num", "Dressings" -> dress, "AngleDefs" -> spAngles] ];
    (* GlobalCollect is ON by default (MakeNTKernel): groups the 306 diagrams by dressing coeff
       (-> 216 traces), COEN collects the factored coefficient across diagrams like FORM -> ~30%
       runtime cut (struct 1/4/7: 2.13/2.31/1.71x -> 1.53/1.58/1.22x FORM), generation cost unchanged. *)
    Print["[time] Numeric backend ZAqbq", lab, ": ", tN, " s"]];
  , {s, 1, 3}];
Print["DONE: ZAqbq{1,4,7} 147 numeric kernels generated."];
