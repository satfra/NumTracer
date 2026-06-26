(* ::Package:: *)

(* DiFfRG_compat.m — DiFfRG integration layer for NumTracer.

   Two public functions hide all the DiFfRG-specific wiring so a flow notebook reads as cleanly
   as the FormTracer reference (one trace line + one emit line + one UpdateFlows per flow):

     MakeNTKernelDiFfRG[ntk, "Name"->..., "Integrator"->..., "Parameters"->..., ...]
        = DiFfRG MakeKernel[0.] scaffold  +  NumTracer numeric MakeNTKernel  (with the DiFfRG
          emission constants baked in and the Regulator alias emitted natively).
     UpdateNTFlows[name]
        = DiFfRG UpdateFlows[name]  +  an idempotent flows/CMakeLists.txt patch (find_package +
          NumTracer link libs + UNITY_BUILD OFF), bundled so the patch can never be left un-applied.

   Loaded by NumTracer.m inside Begin["`Private`"]; the public symbols are declared there. DiFfRG
   symbols are referenced by explicit context so they bind to DiFfRG's whether it is loaded before
   or after NumTracer. *)


(* ---- DiFfRG flow directory (public delayed symbol, trailing slash) -------------------------- *)
ntFlowDir[dir_String] := dir;
ntFlowDir[Automatic]  := Module[{d = DiFfRG`CodeTools`Directory`flowDir},
  If[StringQ[d], d,
    (* fall back to DiFfRG's own default computation if the symbol is unset *)
    FileNameJoin[{If[$Notebooks, NotebookDirectory[], Directory[]], "flows"}]]];


(* ============================================================================================
   MakeNTKernelDiFfRG — scaffold + numeric kernel emission for one flow.
   ============================================================================================ *)

Options[MakeNTKernelDiFfRG] = {
  "Name"                 -> Automatic,   (* REQUIRED flow name, e.g. "ZA" -> dir flows/ZA, class ZA_kernel *)
  "Integrator"           -> Automatic,   (* REQUIRED DiFfRG integrator template, e.g. "Integrator_p2_1ang" *)
  "IntegrationVariables" -> Automatic,   (* REQUIRED, e.g. {"l1","cos1"} *)
  "Parameters"           -> Automatic,   (* REQUIRED kernelParameterList; also the source of Dressings/DressingType *)
  "Namespace"            -> Automatic,   (* C++ gen namespace tag; Automatic -> ToLowerCase[Name] *)
  "AngleDefs"            -> {},          (* NumTracer symbolic angle map (spAngles3/4) *)
  "Constant"             -> 0.,          (* loop-independent flat-added term -> constant(p,k,dressings) *)
  (* DiFfRG emission constants — baked, overridable *)
  "Coordinates"          -> {"LogarithmicCoordinates1D<double>"},
  "CoordinateArguments"  -> {"p"},
  "Device"               -> "GPU",
  "d"                    -> 4,
  "AD"                   -> False,
  "ctype"                -> "double",
  "Type"                 -> "double",
  "Decorator"            -> Automatic,   (* Automatic -> derived from Device *)
  "FlowDirectory"        -> Automatic,   (* Automatic -> DiFfRG`CodeTools`flowDir *)
  "GenDirectory"         -> Automatic    (* Automatic -> a "gen" sibling of the flow directory *)
};

MakeNTKernelDiFfRG::noname   = "\"Name\" is required (the flow name, e.g. \"ZA\").";
MakeNTKernelDiFfRG::noparams = "\"Parameters\" is required (the DiFfRG kernelParameterList).";
MakeNTKernelDiFfRG::nointeg  = "\"Integrator\" and \"IntegrationVariables\" are required.";
MakeNTKernelDiFfRG::mixtype  = "Parameters declare more than one interpolator type `1`; using the first.";

(* Positional second argument: the constant, mirroring DiFfRG's MakeKernel[kernelExpr, constExpr, ...].
   `constExpr` is a plain Mathematica expression in p/k and the dressing names (e.g. ZA[p]), NOT an
   NTKernel — the integrand still comes from `ntk`. Pass the constant EITHER positionally OR via the
   "Constant" option; if both are given the positional argument wins (it is spliced ahead of opts). *)
MakeNTKernelDiFfRG[ntk_NTKernel, constExpr_ /; Head[constExpr] =!= Rule && Head[constExpr] =!= NTKernel,
    opts : OptionsPattern[]] :=
  MakeNTKernelDiFfRG[ntk, "Constant" -> constExpr, opts];

MakeNTKernelDiFfRG[ntk_NTKernel, opts : OptionsPattern[]] := Module[
  {name, nsTag, params, dress, dressTys, dressTy, scalarParams, device, decor, body,
   flowDir, genDir, kernelDir, genFile, kernelFile, tracesFile},

  name = OptionValue["Name"];
  If[name === Automatic, Message[MakeNTKernelDiFfRG::noname]; Abort[]];
  params = OptionValue["Parameters"];
  If[params === Automatic, Message[MakeNTKernelDiFfRG::noparams]; Abort[]];
  If[OptionValue["Integrator"] === Automatic || OptionValue["IntegrationVariables"] === Automatic,
    Message[MakeNTKernelDiFfRG::nointeg]; Abort[]];

  nsTag  = OptionValue["Namespace"] /. Automatic :> ToLowerCase[name];
  device = OptionValue["Device"];
  decor  = OptionValue["Decorator"] /. Automatic :>
             If[device === "GPU", "static __host__ __device__ inline", "static inline"];
  body   = 0.;   (* placeholder scaffold body — NumTracer overwrites kernel.hh below *)

  (* auto-derive the dressing names + their common interpolator type from the non-"double" params *)
  dress    = Cases[params, a_?AssociationQ /; a["Type"] =!= "double" :> a["Name"]];
  dressTys = DeleteDuplicates[Cases[params, a_?AssociationQ /; a["Type"] =!= "double" :> a["Type"]]];
  dressTy  = Switch[Length[dressTys], 0, Automatic, 1, First[dressTys],
               _, Message[MakeNTKernelDiFfRG::mixtype, dressTys]; First[dressTys]];

  (* scalar "double" params (etaPiL, d1V, rhoL, ...) beyond k/p: these are forwarded by DiFfRG's
     integrator between k and the interpolators, so NumTracer must declare them in the kernel /
     constant signature (k and p are already supplied through the NumTrace Args). *)
  scalarParams = Cases[params, a_?AssociationQ /; a["Type"] === "double" &&
     ! MemberQ[{"k", "p"}, If[StringQ[a["Name"]], a["Name"], SymbolName[a["Name"]]]] :> a["Name"]];

  (* resolve paths *)
  flowDir   = ntFlowDir[OptionValue["FlowDirectory"]];
  genDir    = OptionValue["GenDirectory"] /. Automatic :> FileNameJoin[{ParentDirectory[flowDir], "gen"}];
  If[! DirectoryQ[genDir], CreateDirectory[genDir, CreateIntermediateDirectories -> True]];
  kernelDir  = FileNameJoin[{flowDir, name}];
  genFile    = FileNameJoin[{genDir, "gen_" <> nsTag <> "_num.cpp"}];
  kernelFile = FileNameJoin[{kernelDir, "kernel.hh"}];
  tracesFile = FileNameJoin[{kernelDir, "kernels.hh"}];

  (* (1) DiFfRG scaffold FIRST: lays down flows/<name>/ incl. the integrator TUs + a placeholder kernel.hh *)
  DiFfRG`CodeTools`MakeKernel`MakeKernel[body,
    "Name" -> name, "Integrator" -> OptionValue["Integrator"],
    "d" -> OptionValue["d"], "AD" -> OptionValue["AD"], "ctype" -> OptionValue["ctype"],
    "Device" -> device, "Type" -> OptionValue["Type"],
    "Parameters" -> params, "IntegrationVariables" -> OptionValue["IntegrationVariables"],
    "Coordinates" -> OptionValue["Coordinates"],
    "CoordinateArguments" -> OptionValue["CoordinateArguments"]];

  (* (2) NumTracer overwrites kernel.hh + writes kernels.hh with the real, numerically-traced kernel *)
  MakeNTKernel[ntk, genFile, kernelFile, tracesFile,
    "Backend" -> "Numeric", "Name" -> name <> "_kernel", "Namespace" -> nsTag,
    "AngleDefs" -> OptionValue["AngleDefs"], "Decorator" -> decor,
    "Dressings" -> dress, "DressingType" -> dressTy, "ScalarParams" -> scalarParams,
    "Constant" -> OptionValue["Constant"],
    "RuntimeInclude" -> None,
    "ExtraIncludes" -> {"DiFfRG/physics/interpolation.hh", "DiFfRG/physics/physics.hh"},
    "KernelNamespace" -> "DiFfRG", "SupportNamespace" -> "DiFfRG",
    "RegulatorAlias" -> True];

  kernelFile];


(* ============================================================================================
   UpdateNTFlows — DiFfRG UpdateFlows + idempotent NumTracer CMake patch (atomic).
   ============================================================================================ *)

Options[UpdateNTFlows] = {
  "FlowDirectory"  -> Automatic,
  "NumTracerHints" -> "~/.local/share/NumTracer",
  "UnityBuild"     -> False
};

UpdateNTFlows::nocmake   = "Expected the flows CMakeLists at `1` (generate a flow first).";
UpdateNTFlows::patchfail = "Patched `1` but it does not reference NumTracer — the DiFfRG CMake \
template changed; the find_package / link strings in DiFfRG_compat.m need updating.";

UpdateNTFlows[name_String, opts : OptionsPattern[]] := Module[
  {flowDir, f, txt},
  flowDir = ntFlowDir[OptionValue["FlowDirectory"]];

  (* (1) DiFfRG aggregation — regenerates flows/CMakeLists.txt from its template (wipes any prior patch) *)
  DiFfRG`CodeTools`UpdateFlows[name];

  f = FileNameJoin[{flowDir, "CMakeLists.txt"}];
  If[! FileExistsQ[f], Message[UpdateNTFlows::nocmake, f]; Abort[]];
  txt = Import[f, "Text"];

  (* (2) idempotent patch: find_package(NumTracer) + link NumTracer::* + UNITY_BUILD OFF *)
  If[! StringContainsQ[txt, "find_package(NumTracer"],
    txt = StringReplace[txt,
      "add_library(" <> name <> " STATIC ${" <> name <> "_SOURCES})" ->
      "find_package(NumTracer REQUIRED HINTS " <> OptionValue["NumTracerHints"] <> ")\n\n" <>
      "add_library(" <> name <> " STATIC ${" <> name <> "_SOURCES})"];
    txt = StringReplace[txt,
      "target_link_libraries(" <> name <> " DiFfRG::DiFfRG " <> name <> "_nowarn)" ->
      "target_link_libraries(" <> name <> " DiFfRG::DiFfRG " <> name <>
        "_nowarn NumTracer::NumTracer NumTracer::constexpr_budget)"]];
  If[! TrueQ[OptionValue["UnityBuild"]],
    txt = StringReplace[txt, "UNITY_BUILD ON" -> "UNITY_BUILD OFF"]];
  Export[f, txt, "Text"];

  (* (3) loud failure on template drift — StringReplace no-ops silently on a mismatch *)
  If[! StringContainsQ[txt, "NumTracer::NumTracer"],
    Message[UpdateNTFlows::patchfail, f]; Abort[]];

  f];
