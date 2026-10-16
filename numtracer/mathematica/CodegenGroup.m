(* ::Package:: *)

(* ============================================================================================
   Kernel groups: several flows that are projections of the same loop, lowered into ONE kernel that
   returns all of them as a DiFfRG::MultiValue. Each member is generated as an ordinary flow (nets,
   traces, generator, probe, manifest -- see mkGenerateKernel), and its "Front" is lowered here
   together with the others, so that everything the members share is computed once: the angles,
   every dressing lookup and coordinate transform (ntShareGridLookups over all members), every
   dressing atom, and the common subexpressions of the integrands (one COEN function over the list
   of integrands).

   What stays per member: its trace functions (its own namespace and fenv buffer) and its constant.
   ============================================================================================ *)

MakeNTKernelGroupDiFfRG::differ = "The members of kernel group \"`1`\" differ in `2`; group only flows that are projections of the same loop with the same parameters, frame and Matsubara traits.";

(* The fenv buffer of member j; with "FenvReuse" all members share one, filled and traced in turn. *)
ntGroupFenv[j_Integer, reuse_] := If[TrueQ[reuse], "fenv", "fenv" <> ToString[j - 1]];

(* The dressing atoms of all members, each distinct one declared once. Returns the declarations and,
   per member, the fill() arguments in that member's atom order. *)
ntGroupDrAtoms[fronts_List, render_] :=
  Module[{atoms, union, pos, names},
    atoms = Map[If[#["HasDr"], #["DrAtoms"], {}]&, fronts];
    union = DeleteDuplicates[Join @@ atoms];
    pos = AssociationThread[union -> Range[Length[union]]];
    names = Map[
      Flatten[Function[a, With[{nm = "dr_" <> ToString[pos[a] - 1]},
        If[ntDrComplexQ[a], {nm, nm <> "_im"}, {nm}]]] /@ #]&, atoms];
    <|"Decls" -> ntDrDecls[union, render], "FillArgs" -> names|>];

(* The `#if` condition that every member's probe verdict is `v`. *)
ntGroupVerdict[macros_List, v_String] := StringRiffle[("(" <> # <> " == " <> v <> ")")& /@ macros, " && "];

ntGroupKernel[fronts_List, className_String, consts_List, multiReturn_String, shareDress_List, sharedGrids_List,
              fenvReuse_, groupDir_String, kernelFile_String] :=
  Block[{$ntComplexRuntimeProjection = fronts[[1, "Options"]]["ComplexRuntimeProjection"]},
  Module[{o = fronts[[1, "Options"]], n = Length[fronts], sig = fronts[[1, "Signature"]], groupName, macros,
          complexQ, evens, parts, lookups, dr, ints, traceDecls, preamble, body, bodies, kernelFn, constFn,
          traits, tracesIncludes, verdictIncludes, header},
    groupName = StringReplace[className, "_kernel" ~~ EndOfString -> ""];
    (* -- the members must agree in everything but their projection -- *)
    Scan[
      If[!SameQ @@ #[[2]], Message[MakeNTKernelGroupDiFfRG::differ, groupName, #[[1]]]; Abort[]]&,
      {{"their kernel parameters", Map[#["Signature"]["KernelParams"]&, fronts]},
       {"their constant() parameters", Map[#["Signature"]["ConstParams"]&, fronts]},
       {"their frame (fill arguments)", fronts[[All, "FillArgs"]]},
       {"being real or complex", fronts[[All, "ComplexQ"]]}}];
    complexQ = fronts[[1, "ComplexQ"]];
    parts = fronts[[All, "Partition"]];
    If[AnyTrue[parts, #["Split"]&] || !SameQ @@ Map[#["FiniteExtent"]&, parts],
      Message[MakeNTKernelGroupDiFfRG::differ, groupName, "their Matsubara finite-extent classification"]; Abort[]];
    evens = Map[ntMatsubaraEvenQ[#["Integrand"], #["MatsubaraSym"], #["MVarIdx"], #["SymDefs"]]&, fronts];
    If[!SameQ @@ evens,
      Message[MakeNTKernelGroupDiFfRG::differ, groupName, "being even in the Matsubara frequency"]; Abort[]];
    macros = Map[#["Options"]["VerdictMacro"]&, fronts];

    (* -- the shared lookups, over the dressing atoms and integrands of all members -- *)
    lookups =
      If[shareDress =!= {},
        ntShareGridLookups[Join @@ Map[If[#["HasDr"], #["DrAtoms"], {}]&, fronts], fronts[[All, "Integrand"]],
          shareDress, sharedGrids],
        <|"Decls" -> {}, "Rules" -> {}|>];
    dr = ntGroupDrAtoms[fronts, cppFlat[# /. lookups["Rules"]]&];

    (* -- each member's integrand reads its traces from its own buffer, or, reusing one buffer, from
          scalars taken right after its fill -- *)
    traceDecls = ConstantArray[{}, n];
    ints = Table[
      With[{fenv = ntGroupFenv[j, fenvReuse]},
        If[TrueQ[fenvReuse],
          Module[{refs = DeleteDuplicates @ Cases[fronts[[j, "Integrand"]], s_String /; StringContainsQ[s, "::tr"], {0, Infinity}],
                  syms},
            syms = Table["ntt" <> ToString[j - 1] <> "_" <> ToString[i - 1], {i, Length[refs]}];
            traceDecls[[j]] = MapThread["[[maybe_unused]] const auto " <> #2 <> " = " <> #1 <> ";"&, {refs, syms}];
            fronts[[j, "Integrand"]] /. Dispatch[Thread[refs -> syms]]],
          fronts[[j, "Integrand"]] /. s_String :> StringReplace[s, "(fenv)" -> "(" <> fenv <> ")"]]],
      {j, n}];

    preamble = StringRiffle[Join[
        ntSupportUsings[o["SupportNamespace"]],
        ("const " <> $ntRealT <> " " <> SymbolName[First[#]] <> " = " <> cppFlat[Last[#]] <> ";")& /@ o["AngleDefs"],
        lookups["Decls"],
        If[TrueQ[fenvReuse],
          {$ntRealT <> " fenv[" <> Fold["((" <> #1 <> ") > (" <> #2 <> ") ? (" <> #1 <> ") : (" <> #2 <> "))"&, "1",
             Map[#["Options"]["NsHome"] <> "::nenv"&, fronts]] <> "];"},
          Table[$ntRealT <> " " <> ntGroupFenv[j, False] <> "[(" <> fronts[[j, "Options"]]["NsHome"] <> "::nenv) > 0 ? (" <>
            fronts[[j, "Options"]]["NsHome"] <> "::nenv) : 1];", {j, n}]],
        dr["Decls"],
        Flatten @ Table[{
            fronts[[j, "Options"]]["NsHome"] <> "::fill(" <> ntGroupFenv[j, fenvReuse] <> ", " <>
              StringRiffle[Join[SymbolName /@ fronts[[j, "FillArgs"]], dr["FillArgs"][[j]]], ", "] <> ");",
            traceDecls[[j]]}, {j, n}]],
      "\n"];

    (* -- the kernel: one body per probe verdict, each member taking the same branch -- *)
    body = Function[{label, proj},
      FunKit`MakeCppFunction[Map[proj[#] /. lookups["Rules"]&, ints], "Name" -> "kernel", "Prefix" -> o["Decorator"],
        "Return" -> "auto", "CodeParser" -> "Cpp", "Parameters" -> sig["KernelParams"], "Body" -> preamble,
        "MultiReturn" -> multiReturn]];
    kernelFn =
      Which[
        !complexQ, body["real", Identity],
        o["ComplexEndProjection"], body["EndRe", Global`ntRe[#]&],
        o["RealOutput"],
          StringRiffle[{
            "#if " <> ntGroupVerdict[macros, "2"],
            body["Pure", ntPureIntegrand],
            "#elif " <> StringRiffle[("(" <> # <> " != 2)")& /@ macros, " && "],
            "#  if " <> StringRiffle[("(" <> # <> " == 0)")& /@ macros, " || "],
            "#    warning \"NumTracer: a member of kernel group '" <> groupName <> "' probed GENUINELY COMPLEX (verdict 0) but was generated with RealOutput -> True; its imaginary part is discarded.\"",
            "#  endif",
            body["RePart", ntRePartIntegrand],
            "#else",
            "#error \"NumTracer: the members of kernel group '" <> groupName <> "' differ in their probe verdict; generate them as separate flows.\"",
            "#endif"}, "\n"],
        True,
          StringRiffle[{
            "#if " <> ntGroupVerdict[macros, "2"] <> "   // Pure: the Complex -> Re projection is exact",
            body["Pure", ntPureIntegrand],
            "#elif " <> ntGroupVerdict[macros, "1"] <> "   // RePart",
            body["RePart", ntRePartIntegrand],
            "#elif " <> ntGroupVerdict[macros, "0"] <> "   // genuinely complex",
            body["Complex", Identity],
            "#else",
            "#error \"NumTracer: the members of kernel group '" <> groupName <> "' differ in their probe verdict; generate them as separate flows.\"",
            "#endif"}, "\n"]];
    constFn = FunKit`MakeCppFunction[consts, "Name" -> "constant", "Prefix" -> o["Decorator"], "Return" -> "auto",
      "CodeParser" -> "Cpp", "Parameters" -> sig["ConstParams"], "MultiReturn" -> multiReturn,
      "Body" -> If[AllTrue[consts, MatchQ[#, 0 | 0.]&], "", StringRiffle[ntSupportUsings[o["SupportNamespace"]], "\n"]]];

    traits = Join[
      If[First[evens],
        {"static constexpr bool matsubara_even = " <>
          StringRiffle[Map[If[#["MVarIdx"] >= 0, o["KernelNamespace"] <> "::" <> #["Options"]["Namespace"] <> "::matsubara_even", "true"]&, fronts], " && "] <> ";"},
        {}],
      If[parts[[1, "FiniteExtent"]], {"static constexpr bool matsubara_finite_extent = true;"}, {}]];

    tracesIncludes = Map[ntRelativePath[groupDir, #["TracesFile"]]&, fronts];
    verdictIncludes =
      If[complexQ && !o["ComplexEndProjection"],
        Map[FileNameJoin[{FileNameDrop[ntRelativePath[groupDir, #["TracesFile"]], -1], ntVerdictFile}]&, fronts],
        {}];
    header =
      ntApplyTraceComplexOverride[
        FunKit`MakeCppHeader[
          "Includes" -> Join[o["ExtraIncludes"], {"DiFfRG/common/multi_value.hh"},
            If[o["RuntimeInclude"] === None || o["RuntimeInclude"] === "", {}, {o["RuntimeInclude"]}],
            {"numtracer/sun/sun_data.hpp"}, tracesIncludes, verdictIncludes],
          "Body" -> {"namespace " <> o["KernelNamespace"] <> "\n{",
            ntKernelClass[className, Join[traits, {kernelFn, constFn}], o["Decorator"], o["RegulatorTemplate"],
              o["RegulatorAlias"], If[complexQ, {ntReImAccessors[o["Decorator"]]}, {}]],
            "}", "using " <> o["KernelNamespace"] <> "::" <> className <> ";"}],
        First[tracesIncludes], o["KernelNamespace"], o["SupportNamespace"], complexQ];
    ntExportCpp[kernelFile, header];
    Print["wrote group kernel: ", kernelFile, " (", n, " outputs)"];
    kernelFile]];
