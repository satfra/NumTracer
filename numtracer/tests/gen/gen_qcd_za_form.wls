#!/usr/bin/env wolframscript
(* FormTracer (DiFfRG FORM-codegen) clean-cache TIMING for the QCD gluon self-energy (ZA),
   replicating QCD_vacuum_base/QCD.m's ZA block. A fresh empty TraceCache dir forces FormTrace to
   actually run tform (the "full, real generation time"); the shared TBCache stays warm. Times the
   FORM trace and the C++ codegen (MakeKernel) so the FormTracer cost sits beside NumTracer's.
   Run from numtracer/:   wolfram -script tests/gen/gen_qcd_za_form.wls                          *)

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
$cacheDir = FileNameJoin[{$repo, "numtracer", "tests", "gen", "TraceCache_zaform"}];
If[DirectoryQ[$cacheDir], DeleteDirectory[$cacheDir, DeleteContents -> True]];
FSetCacheDirectory[$cacheDir];

{tTrace, FlowAA} = AbsoluteTiming[
  Module[{fRG, tr},
    fRG = FTakeDerivatives[WetterichEquation, {A[i1], A[i2]}] // FTruncate // FPlot // FRoute // FPrint;
    tr  = FTerm[TBGetProjector["AA", 1, {i1, i2} /. fRG["1-Loop"]["ExternalIndices"]]] **
            (fRG["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]);
    FormTrace["ZA", tr] // dressingRules // PropParam // Simplify]];
Print["[time] FORM trace (FTakeDerivatives..FormTrace..PropParam): ", tTrace, " s"];

tKernel = First @ AbsoluteTiming[
  MakeKernel[FlowAA/p^2,
    "Name" -> "ZA", "Integrator" -> "Integrator_p2_1ang", "d" -> 4, "AD" -> False,
    "ctype" -> "double", "Device" -> "GPU", "Type" -> "double",
    "Parameters" -> kernelParameterList,
    "IntegrationVariables" -> {"l1", "cos1"},
    "Coordinates" -> {"LogarithmicCoordinates1D<double>"}, "CoordinateArguments" -> {"p"}]];
Print["[time] FormTracer MakeKernel (FORM codegen + CSE): ", tKernel, " s"];
Print["[time] FormTracer generation TOTAL: ", tTrace + tKernel, " s"];
Print["DONE ZA FORM kernel."];
