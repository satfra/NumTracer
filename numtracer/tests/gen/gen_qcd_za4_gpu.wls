#!/usr/bin/env wolframscript
(* QCD generator: the CUDA (__host__ __device__) variant of the ZA4 invariant-basis kernel.
   Identical physics/pipeline to gen_qcd_za4.wls; only the emission differs:
     * "Decorator" -> "static __host__ __device__ inline" so every emitted function
       (kernel/constant, regulator wrappers, fill/trN/powr) is device-callable,
     * distinct class/namespace (ZA4_gpu_kernel / za4_gpu) so the GPU test binary can link
       the host AND device kernels together without ODR collisions,
     * the generator program goes to $TemporaryDirectory: its content differs from the
       committed gen_qcd_za4_inv.cpp only in the baked namespace default, so committing a
       near-duplicate ~650K source set per variant would be pure churn. (The committed host
       generator can re-emit this header too: gen_za4_inv -n za4_gpu -d '<decorator>'.)
   Run from numtracer/:   wolfram -script tests/gen/gen_qcd_za4_gpu.wls                        *)

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

frame = NumTracer`sp4Frame[Global`p, Global`l1, Global`cos1, Global`cos2, Global`phi,
          Global`p1, Global`p2, Global`p3, Global`p4, Global`l1];
args  = {Global`l1, Global`cos1, Global`cos2, Global`phi, Global`p, Global`k};
dress = {ZA3, ZAcbc, ZA4, ZAqbq1, dtZc, Zc, dtZA, ZA, dtZq, Zq, Mq};

(* SP4 loop-external cosines (DeclareSymmetricPoints4DP4) — the SAME formulas sp4Frame uses. *)
spAngles4 = {
  Global`cosl1p1 -> Sqrt[1 - Global`cos1^2] Global`cos2,
  Global`cosl1p2 -> -(1/3) Sqrt[1 - Global`cos1^2] (Global`cos2 - 2 Sqrt[2 - 2 Global`cos2^2] Cos[Global`phi]),
  Global`cosl1p3 -> -(1/3) Sqrt[1 - Global`cos1^2] (Global`cos2 + Sqrt[2 - 2 Global`cos2^2] (Cos[Global`phi] - Sqrt[3] Sin[Global`phi])),
  Global`cosl1p4 -> -(1/3) Sqrt[1 - Global`cos1^2] (Global`cos2 + Sqrt[2 - 2 Global`cos2^2] (Cos[Global`phi] + Sqrt[3] Sin[Global`phi]))};

{tTrace, traceA4} = AbsoluteTiming[
  Module[{fRG, proj},
    fRG  = FTakeDerivatives[WetterichEquation, {A[i1], A[i2], A[i3], A[i4]}] // FTruncate // FPlot // FRoute // FPrint;
    proj = FTerm[TBGetProjector["AAAAClassTrans", 1, {i1, i2, i3, i4} /. fRG["1-Loop"]["ExternalIndices"]]];
    proj ** (fRG["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[])]];
Print["[time] A4 trace build (shared): ", tTrace, " s"];

{tNet, net} = AbsoluteTiming[ FromFunKit[ SPParam[traceA4 // dressingRules] ] ];
Print["[time] FromFunKit+SPParam: ", tNet, " s"];
Print["ZA4 net: top terms = ", If[Head[net] === Plus, Length[net], 1]];
Print["Nf baked as integer = ", Global`Nf];

{tAnalyse, ntk} = AbsoluteTiming[ NumTrace[net, "Frame" -> frame, "Args" -> args] ];
Print["[time] NumTrace (analyse): ", tAnalyse, " s   diagrams = ", Length[ntk[[1]]["Diagrams"]]];

tGen = First @ AbsoluteTiming[ MakeNTKernel[ntk,
  FileNameJoin[{$TemporaryDirectory, "gen_qcd_za4_gpu_num.cpp"}],
  FileNameJoin[{$gen, "ZA4_gpu_kernel.hh"}],
  FileNameJoin[{$gen, "ZA4_gpu_kernels.hh"}],
  "Backend" -> "Numeric", "Name" -> "ZA4_gpu_kernel", "Namespace" -> "za4_gpu", "Dressings" -> dress,
  "AngleDefs" -> spAngles4, "Decorator" -> "static __host__ __device__ inline"] ];
Print["[time] MakeNTKernel (Numeric) (contract+lower+compile+run generator): ", tGen, " s"];
Print["ZA4 GPU (__host__ __device__) invariant-basis kernel generated."];
