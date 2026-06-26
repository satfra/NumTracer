#!/usr/bin/env wolframscript
(* QCD generator: the CUDA (__host__ __device__) variant of the ZA3 invariant-basis kernel.
   Identical physics/pipeline to gen_qcd_za3.wls; only the emission differs:
     * "Decorator" -> "static __host__ __device__ inline" so every emitted function
       (kernel/constant, regulator wrappers, fill/trN/powr) is device-callable,
     * distinct class/namespace (ZA3_gpu_kernel / za3_gpu) so the GPU test binary can link
       the host AND device kernels together without ODR collisions,
     * the generator program goes to $TemporaryDirectory: its content differs from the
       committed gen_qcd_za3_inv.cpp only in the baked namespace default, so committing a
       near-duplicate ~100k-line source per variant would be pure churn. (The committed host
       generator can re-emit this header too: gen_za3_inv -n za3_gpu -d '<decorator>'.)
   Run from numtracer/:   wolfram -script tests/gen/gen_qcd_za3_gpu.wls                        *)

If[Quiet[Needs["FunKit`"]] === $Failed, (* installed (Applications dir) or the local checkout *)
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
SetOptions[MakeNTKernel, Sequence @@ $ntShimOpts];  (* tests: emit against the shim + DiFfRG namespace *)

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
Print["ZA3 net: top terms = ", If[Head[net] === Plus, Length[net], 1]];
Print["Nf baked as integer = ", Global`Nf];

{tA, ntk} = AbsoluteTiming[ NumTrace[net, "Frame" -> frame, "Args" -> args] ];
Print["[time] NumTrace (analyse): ", tA, " s   diagrams = ", Length[ntk[[1]]["Diagrams"]]];

tG = First @ AbsoluteTiming[ MakeNTKernel[ntk,
  FileNameJoin[{$TemporaryDirectory, "gen_qcd_za3_gpu_num.cpp"}],
  FileNameJoin[{$gen, "ZA3_gpu_kernel.hh"}],
  FileNameJoin[{$gen, "ZA3_gpu_kernels.hh"}],
  "Backend" -> "Numeric", "Name" -> "ZA3_gpu_kernel", "Namespace" -> "za3_gpu", "Dressings" -> dress,
  "AngleDefs" -> spAngles, "Decorator" -> "static __host__ __device__ inline"] ];
Print["[time] MakeNTKernel (Generate) (reduce+lower+compile+run): ", tG, " s"];
Print["ZA3 GPU (__host__ __device__) invariant-basis kernel generated."];
