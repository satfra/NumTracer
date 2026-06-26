#!/usr/bin/env wolframscript
(* FormTracer (DiFfRG FORM-codegen) ORACLES for the QCD quark-gluon vertex flow with the FULL
   {1,4,7} vertex tensor basis, projected onto each external structure (1, 4, 7). The physical
   cross-check for the ZAqbq{1,4,7}_147 numeric kernels. Emits ZAqbq{1,4,7}_147_kernel.hh into
   refshim/. Reuses gen_qcd_za3_147_form_setup.m (loads qcd147_setup.m -> AqbqDirect147, carries
   ZAqbq1/4/7 dressings + SP3 helpers). Run from numtracer/:
     wolfram -script tests/gen/gen_qcd_aqbq147_form.wls                                            *)

If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];
$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "tests", "gen", "gen_qcd_za3_147_form_setup.m"}]];

(* shared flow derivation *)
fRG  = FTakeDerivatives[WetterichEquation, {A[i1], qb[i2], q[i3]}] // FTruncate // FPlot // FRoute // FPrint;
ext  = fRG["1-Loop"]["ExternalIndices"];
expr = fRG["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[];

projIdx = {1, 2, 3};
labels  = {1, 4, 7};
Do[
  Module[{k = projIdx[[s]], lab = labels[[s]], cacheDir, proj, tr, Flow, tT, tK, out},
    cacheDir = FileNameJoin[{$repo, "numtracer", "tests", "gen", "TraceCache_aqbq" <> ToString[lab] <> "_147form"}];
    If[DirectoryQ[cacheDir], DeleteDirectory[cacheDir, DeleteContents -> True]];
    FSetCacheDirectory[cacheDir];
    {tT, Flow} = AbsoluteTiming[
      proj = TBGetProjector["AqbqDirect147", k, {i1, i3, i2} /. ext] //
               TBProjectToSymmetricPoint[#, l1, p, p1, p2, p3] & // Simplify;
      tr   = FTerm[proj] ** expr;
      (FormTrace["ZAqbq" <> ToString[lab] <> "_147", tr, {}, SP3FormRule] // dressingRules //
         TBProjectToSymmetricPoint[#, l1, p, p1, p2, p3] &) // SPParam // Simplify];
    Print["[time] FORM trace ZAqbq", lab, "_147: ", tT, " s"];
    tK = First @ AbsoluteTiming[
      MakeKernel[Flow,
        "Name" -> "ZAqbq" <> ToString[lab] <> "_147", "Integrator" -> "Integrator_p2_4D_2ang",
        "d" -> 4, "AD" -> False, "ctype" -> "double", "Device" -> "GPU", "Type" -> "double",
        "Parameters" -> kernelParameterList, "KernelBody" -> SP3Defs,
        "IntegrationVariables" -> {"l1", "cos1", "cos2"},
        "Coordinates" -> {"LogarithmicCoordinates1D<double>"}, "CoordinateArguments" -> {"p"}]];
    Print["[time] FormTracer MakeKernel ZAqbq", lab, "_147: ", tK, " s"];
    out = FileNameJoin[{$repo, "numtracer", "tests", "refshim", "ZAqbq" <> ToString[lab] <> "_147_kernel.hh"}];
    If[FileExistsQ["ZAqbq" <> ToString[lab] <> "_147.hh"],
       CopyFile["ZAqbq" <> ToString[lab] <> "_147.hh", out, OverwriteTarget -> True];
       Print["wrote oracle: ", out],
       Print["WARN: ZAqbq", lab, "_147.hh not found in ", Directory[]]]];
  , {s, 1, 3}];
Print["DONE: ZAqbq{1,4,7} 147 FORM oracles generated."];
