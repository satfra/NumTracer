#!/usr/bin/env wolframscript
(* FormTracer (DiFfRG FORM-codegen) clean-cache TIMING for the QCD quark wave-function (Zq),
   replicating QCD_vacuum_base/QCD.m's Zq block (1 angle). Fresh empty TraceCache so FormTrace
   recomputes from scratch; TBCache stays warm. Times FORM trace + MakeKernel codegen.
   Run from numtracer/:   wolfram -script tests/gen/gen_qcd_zq_form.wls                          *)

If[Quiet[Needs["FunKit`"]] === $Failed, (* installed (Applications dir) or the local checkout *)
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];
$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "tests", "gen", "gen_qcd_form_setup.m"}]];

(* clean cold run: empty this kernel's trace cache so FormTrace reruns tform (the user's
   "EMPTY tracecache" requirement); TBCache is shared setup and left warm. *)
$cacheDir = FileNameJoin[{$repo, "numtracer", "tests", "gen", "TraceCache_zqform"}];
If[DirectoryQ[$cacheDir], DeleteDirectory[$cacheDir, DeleteContents -> True]];
FSetCacheDirectory[$cacheDir];

{tTrace, Flowqbq} = AbsoluteTiming[
  Module[{fRG, tr},
    fRG = FTakeDerivatives[WetterichEquation, {qb[i1], q[i2]}] // FTruncate // FSimplify // FPlot // FRoute // FPrint;
    tr  = FTerm[-TBGetProjector["qbq", 1, {i2, i1} /. fRG["1-Loop"]["ExternalIndices"]]] **
            (fRG["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]);
    (FormTrace["Zq", tr] // dressingRules) // PropParam // Simplify]];
Print["[time] FORM trace (FTakeDerivatives..FormTrace..PropParam): ", tTrace, " s"];

tKernel = First @ AbsoluteTiming[
  MakeKernel[Flowqbq,
    "Name" -> "Zq", "Integrator" -> "Integrator_p2_1ang", "d" -> 4, "AD" -> False,
    "ctype" -> "double", "Device" -> "GPU", "Type" -> "double",
    "Parameters" -> kernelParameterList,
    "IntegrationVariables" -> {"l1", "cos1"},
    "Coordinates" -> {"LogarithmicCoordinates1D<double>"}, "CoordinateArguments" -> {"p"}]];
Print["[time] FormTracer MakeKernel (FORM codegen + CSE): ", tKernel, " s"];
Print["[time] FormTracer generation TOTAL: ", tTrace + tKernel, " s"];
Print["DONE Zq FORM kernel."];
