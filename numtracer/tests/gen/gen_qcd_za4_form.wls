#!/usr/bin/env wolframscript
(* FormTracer (DiFfRG FORM-codegen) clean-cache TIMING for the QCD four-gluon vertex (ZA4),
   replicating QCD_vacuum_base/QCD.m's ZA4 block (3 angles, SP4 frame, ParallelMap[Simplify]).
   Fresh empty TraceCache so FormTrace recomputes from scratch; TBCache stays warm. Times FORM
   trace + MakeKernel codegen.
   Run from numtracer/:   wolfram -script tests/gen/gen_qcd_za4_form.wls                         *)

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
$cacheDir = FileNameJoin[{$repo, "numtracer", "tests", "gen", "TraceCache_za4form"}];
If[DirectoryQ[$cacheDir], DeleteDirectory[$cacheDir, DeleteContents -> True]];
FSetCacheDirectory[$cacheDir];

{tTrace, FlowA4} = AbsoluteTiming[
  Module[{fRG, proj, tr, fl},
    fRG  = FTakeDerivatives[WetterichEquation, {A[i1], A[i2], A[i3], A[i4]}] // FTruncate // FPlot // FRoute // FPrint;
    proj = TBGetProjector["AAAAClassTrans", 1, {i1, i2, i3, i4} /. fRG["1-Loop"]["ExternalIndices"]] //
             TBProjectToSymmetricPoint[#, l1, p, p1, p2, p3, p4] & // Simplify;
    tr   = FTerm[proj] ** (fRG["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]);
    fl   = (FormTrace["ZA4", tr, {}, SP4FormRule] // dressingRules //
              TBProjectToSymmetricPoint[#, l1, p, p1, p2, p3, p4] &) // SPParam;
    ParallelMap[Simplify, fl]]];
Print["[time] FORM trace (FTakeDerivatives..FormTrace..SPParam..Simplify): ", tTrace, " s"];

tKernel = First @ AbsoluteTiming[
  MakeKernel[FlowA4,
    "Name" -> "ZA4", "Integrator" -> "Integrator_p2_4D_3ang", "d" -> 4, "AD" -> False,
    "ctype" -> "double", "Device" -> "GPU", "Type" -> "double",
    "Parameters" -> kernelParameterList, "KernelBody" -> SP4Defs,
    "IntegrationVariables" -> {"l1", "cos1", "cos2", "phi"},
    "Coordinates" -> {"LogarithmicCoordinates1D<double>"}, "CoordinateArguments" -> {"p"}]];
Print["[time] FormTracer MakeKernel (FORM codegen + CSE): ", tKernel, " s"];
Print["[time] FormTracer generation TOTAL: ", tTrace + tKernel, " s"];
Print["DONE ZA4 FORM kernel."];
