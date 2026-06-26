#!/usr/bin/env wolframscript
(* FormTracer (DiFfRG FORM-codegen) oracle for the QCD three-gluon vertex (ZA3) with the FULL
   quark-gluon vertex {1,4,7} tensor basis — the physical cross-check for the dense `1/4/7` numeric
   kernel. Emits ZA3_147_kernel.hh into refshim/. Run from numtracer/:
     wolfram -script tests/gen/gen_qcd_za3_147_form.wls                                            *)

If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];
$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "tests", "gen", "gen_qcd_za3_147_form_setup.m"}]];

$cacheDir = FileNameJoin[{$repo, "numtracer", "tests", "gen", "TraceCache_za3_147_form"}];
If[DirectoryQ[$cacheDir], DeleteDirectory[$cacheDir, DeleteContents -> True]];
FSetCacheDirectory[$cacheDir];

{tTrace, FlowA3} = AbsoluteTiming[
  Module[{fRG, proj, tr},
    fRG  = FTakeDerivatives[WetterichEquation, {A[i1], A[i2], A[i3]}] // FTruncate // FPlot // FRoute // FPrint;
    proj = TBGetProjector["AAAClassTrans", 1, {i1, i2, i3} /. fRG["1-Loop"]["ExternalIndices"]] //
             TBProjectToSymmetricPoint[#, l1, p, p1, p2, p3] & // Simplify;
    tr   = FTerm[proj] ** (fRG["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]);
    (FormTrace["ZA3_147", tr, {}, SP3FormRule] // dressingRules //
       TBProjectToSymmetricPoint[#, l1, p, p1, p2, p3] &) // SPParam // Simplify]];
Print["[time] FORM trace (147): ", tTrace, " s"];

tKernel = First @ AbsoluteTiming[
  MakeKernel[FlowA3,
    "Name" -> "ZA3_147", "Integrator" -> "Integrator_p2_4D_2ang", "d" -> 4, "AD" -> False,
    "ctype" -> "double", "Device" -> "GPU", "Type" -> "double",
    "Parameters" -> kernelParameterList, "KernelBody" -> SP3Defs,
    "IntegrationVariables" -> {"l1", "cos1", "cos2"},
    "Coordinates" -> {"LogarithmicCoordinates1D<double>"}, "CoordinateArguments" -> {"p"}]];
Print["[time] FormTracer MakeKernel (147): ", tKernel, " s"];
(* place the generated kernel header where the compare test includes it. *)
$out = FileNameJoin[{$repo, "numtracer", "tests", "refshim", "ZA3_147_kernel.hh"}];
If[FileExistsQ["ZA3_147.hh"], CopyFile["ZA3_147.hh", $out, OverwriteTarget -> True];
   Print["wrote oracle: ", $out],
   Print["WARN: ZA3_147.hh not found in ", Directory[]]];
Print["DONE ZA3 147 FORM kernel."];
