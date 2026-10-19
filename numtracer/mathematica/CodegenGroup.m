(* ::Package:: *)

(* ============================================================================================
   Kernel groups: several flows that are projections of the same loop, lowered into ONE kernel that
   returns all of them as a DiFfRG::MultiValue. Each member is generated as an ordinary flow (nets,
   traces, generator, probe, manifest -- see mkGenerateKernel), and its "Front" is lowered here
   together with the others, so that everything the members share is computed once: the angles,
   every dressing lookup and coordinate transform (ntShareGridLookups over all members), every
   dressing atom, and the common subexpressions of the integrands (one COEN function over the list
   of integrands).

   What stays per member: its trace functions (its own namespace) and its constant. The members are
   traced in turn in one fenv buffer, each one's trace values copied to scalars right after its fill;
   a buffer per member, all live at once, costs more stack than the sharing saves.
   ============================================================================================ *)

MakeNTKernelGroupDiFfRG::differ = "The members of kernel group \"`1`\" differ in `2`; group only flows that are projections of the same loop with the same parameters, frame and Matsubara traits.";

MakeNTKernelGroupDiFfRG::split = "A member of kernel group \"`1`\" splits into a finite-extent and a tail part (matsubara_split); kernel groups do not support that, generate it as a separate flow.";

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
              groupDir_String, kernelFile_String] :=
  Block[{$ntComplexRuntimeProjection = fronts[[1, "Options"]]["ComplexRuntimeProjection"]},
  Module[{o = fronts[[1, "Options"]], n = Length[fronts], sig = fronts[[1, "Signature"]], groupName, macros,
          complexQ, evens, parts, lookups, dr, ints, traceDecls, nenv, preamble, body, kernelFn, hoistTrait, constFn,
          traits, tracesIncludes, verdictIncludes, header},
    groupName = StringReplace[className, "_kernel" ~~ EndOfString -> ""];
    (* -- the members must agree in everything but their projection -- *)
    Scan[
      If[!SameQ @@ #[[2]], Message[MakeNTKernelGroupDiFfRG::differ, groupName, #[[1]]]; Abort[]]&,
      {{"their kernel parameters", Map[#["Signature"]["KernelParams"]&, fronts]},
       {"their constant() parameters", Map[#["Signature"]["ConstParams"]&, fronts]},
       {"their frame (fill arguments)", fronts[[All, "FillArgs"]]},
(* Real and complex members can share the one end-projected body; the probe-verdict bodies cannot. *)
       {"being real or complex", If[o["ComplexEndProjection"], {}, fronts[[All, "ComplexQ"]]]}}];
    complexQ = AnyTrue[fronts[[All, "ComplexQ"]], TrueQ];
    parts = fronts[[All, "Partition"]];
    If[AnyTrue[parts, #["Split"]&], Message[MakeNTKernelGroupDiFfRG::split, groupName]; Abort[]];
    If[!SameQ @@ Map[#["FiniteExtent"]&, parts],
      Message[MakeNTKernelGroupDiFfRG::differ, groupName, "their Matsubara finite-extent classification"]; Abort[]];
    evens = Map[ntMatsubaraEvenQ[#["Integrand"], #["MatsubaraSym"], #["MVarIdx"], #["SymDefs"]]&, fronts];
    If[!SameQ @@ evens,
      Message[MakeNTKernelGroupDiFfRG::differ, groupName, "being even in the Matsubara frequency"]; Abort[]];
    macros = Map[#["Options"]["VerdictMacro"]&, fronts];

    (* -- the shared lookups, over the dressing atoms and integrands of all members -- *)
    lookups =
      If[shareDress =!= {},
        ntShareGridLookups[DeleteDuplicates[Join @@ Map[If[#["HasDr"], #["DrAtoms"], {}]&, fronts]], fronts[[All, "Integrand"]],
          shareDress, sharedGrids],
        <|"Decls" -> {}, "Rules" -> {}|>];
    dr = ntGroupDrAtoms[fronts, cppFlat[# /. lookups["Rules"]]&];

    (* -- each member's integrand reads its trace values from scalars taken right after its fill -- *)
    traceDecls = ConstantArray[{}, n];
    ints = Table[
      Module[{refs = DeleteDuplicates @ Cases[fronts[[j, "Integrand"]],
                s_String /; StringStartsQ[s, fronts[[j, "Options"]]["NsHome"] <> "::tr"], {0, Infinity}],
              syms},
        syms = Table["ntt" <> ToString[j - 1] <> "_" <> ToString[i - 1], {i, Length[refs]}];
        traceDecls[[j]] = MapThread["[[maybe_unused]] const auto " <> #2 <> " = " <> #1 <> ";"&, {refs, syms}];
        fronts[[j, "Integrand"]] /. Dispatch[Thread[refs -> syms]]],
      {j, n}];

    (* the buffer is as large as the largest member's: ntnenv<n-1> = max(1, nenv_1, ..., nenv_n) *)
    nenv = Table["constexpr auto ntnenv" <> ToString[j] <> " = " <>
        If[j == 0, "1", "ntnenv" <> ToString[j - 1] <> " > " <> fronts[[j, "Options"]]["NsHome"] <> "::nenv ? ntnenv" <>
          ToString[j - 1] <> " : " <> fronts[[j, "Options"]]["NsHome"] <> "::nenv"] <> ";", {j, 0, n}];
    preamble = StringRiffle[Join[
        ntSupportUsings[o["SupportNamespace"]],
        ("const " <> $ntRealT <> " " <> SymbolName[First[#]] <> " = " <> cppFlat[Last[#]] <> ";")& /@ o["AngleDefs"],
        lookups["Decls"],
        nenv,
        {$ntRealT <> " fenv[ntnenv" <> ToString[n] <> "];"},
        dr["Decls"],
        Flatten @ Table[{
            fronts[[j, "Options"]]["NsHome"] <> "::fill(fenv, " <>
              StringRiffle[Join[SymbolName /@ fronts[[j, "FillArgs"]], dr["FillArgs"][[j]]], ", "] <> ");",
            traceDecls[[j]]}, {j, n}]],
      "\n"];

    (* -- the kernel: one body per probe verdict, each member taking the same branch -- *)
    body = Function[{proj, transform},
      FunKit`MakeCppFunction[Map[proj[#] /. lookups["Rules"]&, ints], "Name" -> "kernel", "Prefix" -> o["Decorator"],
        "Return" -> "auto", "CodeParser" -> "Cpp", "Parameters" -> sig["KernelParams"], "Body" -> preamble,
        "MultiReturn" -> multiReturn, "ReturnTransform" -> transform]];
    kernelFn =
      Which[
        !complexQ, body[Identity, Identity],
(* ntRe goes on the returned values, not around the integrands: COEN finds no common subexpressions
   in a list of outputs that are each one opaque ntRe[...] call, and then emits them unshared, which
   for the QuarkMeson meson group exhausted 20 GB. *)
        o["ComplexEndProjection"], body[Identity, Global`ntRe[#]&],
        o["RealOutput"],
(* RePart is exact for a Pure (verdict 2) member too, so any mix of verdicts takes it; only verdict 0
   (genuinely complex) loses its imaginary part *)
          StringRiffle[{
            "#if " <> ntGroupVerdict[macros, "2"],
            body[ntPureIntegrand, Identity],
            "#else",
            "#  if " <> StringRiffle[("(" <> # <> " == 0)")& /@ macros, " || "],
            "#    warning \"NumTracer: a member of kernel group '" <> groupName <> "' probed GENUINELY COMPLEX (verdict 0) but was generated with RealOutput -> True; its imaginary part is discarded.\"",
            "#  endif",
            body[ntRePartIntegrand, Identity],
            "#endif"}, "\n"],
        True,
          StringRiffle[{
            "#if " <> ntGroupVerdict[macros, "2"] <> "   // Pure: the Complex -> Re projection is exact",
            body[ntPureIntegrand, Identity],
            "#elif " <> StringRiffle[("(" <> # <> " != 0)")& /@ macros, " && "] <> "   // RePart, exact for every real member",
            body[ntRePartIntegrand, Identity],
            "#elif " <> ntGroupVerdict[macros, "0"] <> "   // genuinely complex",
            body[Identity, Identity],
            "#else",
            "#error \"NumTracer: kernel group '" <> groupName <> "' mixes genuinely complex (verdict 0) and real members; generate them as separate flows.\"",
            "#endif"}, "\n"]];
    {kernelFn, hoistTrait} = ntHoistedKernelFn[kernelFn, o, parts[[1, "Split"]], fronts[[1, "MatsubaraSym"]]];
    constFn = FunKit`MakeCppFunction[consts, "Name" -> "constant", "Prefix" -> o["Decorator"], "Return" -> "auto",
      "CodeParser" -> "Cpp", "Parameters" -> sig["ConstParams"], "MultiReturn" -> multiReturn,
      "Body" -> If[AllTrue[consts, MatchQ[#, 0 | 0.]&], "", StringRiffle[ntSupportUsings[o["SupportNamespace"]], "\n"]]];

    traits = Join[
      If[First[evens],
        {"static constexpr bool matsubara_even = " <>
          StringRiffle[Map[If[#["MVarIdx"] >= 0, #["Options"]["NsHome"] <> "::matsubara_even", "true"]&, fronts], " && "] <> ";"},
        {}],
      If[parts[[1, "FiniteExtent"]], {"static constexpr bool matsubara_finite_extent = true;"}, {}],
      hoistTrait];

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
