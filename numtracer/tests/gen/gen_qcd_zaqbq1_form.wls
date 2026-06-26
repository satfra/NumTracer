#!/usr/bin/env wolframscript
(* FormTracer (DiFfRG FORM-codegen) clean-cache TIMING for the QCD quark-gluon vertex (ZAqbq1),
   replicating QCD_vacuum_base/QCD.m's ZAqbq1 block (2 angles, SP3 frame). Fresh empty TraceCache
   so FormTrace recomputes from scratch; TBCache stays warm. Times FORM trace + MakeKernel codegen.
   Run from numtracer/:   wolfram -script tests/gen/gen_qcd_zaqbq1_form.wls                       *)

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
$cacheDir = FileNameJoin[{$repo, "numtracer", "tests", "gen", "TraceCache_zaqbq1form"}];
If[DirectoryQ[$cacheDir], DeleteDirectory[$cacheDir, DeleteContents -> True]];
FSetCacheDirectory[$cacheDir];

{tTrace, FlowAqbq1} = AbsoluteTiming[
  Module[{fRG, proj, tr},
    fRG  = FTakeDerivatives[WetterichEquation, {A[i1], qb[i2], q[i3]}] // FTruncate // FPlot // FRoute // FPrint;
    proj = TBGetProjector["AqbqDirect1", 1, {i1, i3, i2} /. fRG["1-Loop"]["ExternalIndices"]];
    tr   = FTerm[proj] ** (fRG["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]) //
             TBProjectToSymmetricPoint[#, l1, p, p1, p2, p3] & // Simplify;
    FormTrace["ZAqbq1", tr, {}, SP3FormRule] // dressingRules //
      TBProjectToSymmetricPoint[#, l1, p, p1, p2, p3] & // SPParam // Simplify]];
Print["[time] FORM trace (FTakeDerivatives..FormTrace..SPParam): ", tTrace, " s"];

tKernel = First @ AbsoluteTiming[
  MakeKernel[FlowAqbq1,
    "Name" -> "ZAqbq1", "Integrator" -> "Integrator_p2_4D_2ang", "d" -> 4, "AD" -> False,
    "ctype" -> "double", "Device" -> "GPU", "Type" -> "double",
    "Parameters" -> kernelParameterList, "KernelBody" -> SP3Defs,
    "IntegrationVariables" -> {"l1", "cos1", "cos2"},
    "Coordinates" -> {"LogarithmicCoordinates1D<double>"}, "CoordinateArguments" -> {"p"}]];
Print["[time] FormTracer MakeKernel (FORM codegen + CSE): ", tKernel, " s"];
Print["[time] FormTracer generation TOTAL: ", tTrace + tKernel, " s"];
Print["DONE ZAqbq1 FORM kernel."];
