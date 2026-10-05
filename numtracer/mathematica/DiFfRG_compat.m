(* ::Package:: *)

(* DiFfRG_compat.m — DiFfRG integration layer for NumTracer.

   Two public functions hide all the DiFfRG-specific wiring so a flow notebook reads as cleanly
   as the FormTracer reference (one trace line + one emit line + one UpdateFlows per flow):

     MakeNTKernelDiFfRG[ntk, "Name"->..., "Integrator"->..., "Parameters"->..., ...]
        = DiFfRG MakeKernel[0.] scaffold  +  NumTracer numeric MakeNTKernel  (with the DiFfRG
          emission constants baked in and the DiFfRG kernel shape — `template<typename REG>` plus
          the `using Regulator = REG;` alias — emitted natively; the general NumTracer emission is a
          plain class with consumer-supplied regulators).
     UpdateNTFlows[name]
        = DiFfRG UpdateFlows[name]  +  an idempotent flows/CMakeLists.txt patch (find_package +
          NumTracer link libs + UNITY_BUILD OFF), bundled so the patch can never be left un-applied.

   Loaded by NumTracer.m inside Begin["`Private`"]; the public symbols are declared there. DiFfRG
   symbols are referenced by explicit context so they bind to DiFfRG's whether it is loaded before
   or after NumTracer. *)

(* ---- DiFfRG chatter capture ------------------------------------------------------------------
   DiFfRG Prints a line per emitted file plus a "Please run UpdateFlows[]" nudge, which is WRONG here:
   bare UpdateFlows drops the NumTracer CMake patch (use UpdateNTFlows). So capture the chatter, print
   one NumTracer line instead, and keep the raw lines behind $NumTracerVerbose.
   Internal`InheritedBlock restores Print even if DiFfRG aborts. *)

SetAttributes[ntCapturePrint, HoldFirst];

ntCapturePrint[expr_] :=
  Module[{lines = {}, res},
    res =
      Internal`InheritedBlock[
        {Print},
        Unprotect[Print];
        Print[args___] := AppendTo[lines, StringJoin[ToString /@ {args}]];
        expr
      ];
    {res, lines}
  ];

(* Does flows/CMakeLists.txt still need an UpdateNTFlows pass for this flow? True when the file is
   missing, has lost the NumTracer patch, does not mention this flow, or DiFfRG wrote files. *)

ntCMakeStaleQ[flowDir_, name_String, wrote_] :=
  Module[{f = FileNameJoin[{flowDir, "CMakeLists.txt"}], txt},
    If[!FileExistsQ[f],
      Return[True]
    ];
    txt = Quiet @ Check[Import[f, "Text"], ""];
    !StringContainsQ[txt, "find_package(NumTracer"] || !StringContainsQ[txt, name] || wrote =!= {}
  ];

(* One NumTracer-native line in place of DiFfRG's per-file chatter. *)

ntReportDiFfRG[name_String, flowDir_, lines_List] :=
  Module[{unchanged, wrote},
    unchanged = Select[lines, StringContainsQ[#, "unchanged"]&];
    (* everything that is neither an "unchanged" note nor the UpdateFlows nudge = a file DiFfRG wrote *)
    wrote = Select[lines, !StringContainsQ[#, "unchanged"] && !StringContainsQ[#, "UpdateFlows"]&];
    ntLog["[DiFfRG] ", #]& /@ lines;(* raw chatter, verbose only *)
    Print["[NumTracer] " <> name <> ": DiFfRG scaffold — " <> ToString[Length[wrote]] <> " written, " <> ToString[Length[unchanged]] <> " unchanged"];
    If[ntCMakeStaleQ[flowDir, name, wrote],
      Print["[NumTracer] " <> name <> ": flows/CMakeLists.txt needs refreshing — run UpdateNTFlows[\"" <> name <> "\"]  (NOT UpdateFlows[], which would drop the NumTracer CMake patch)"]
    ];
    wrote
  ];

(* Resolve the "MatsubaraVar" option to a symbol name (or None). Automatic follows DiFfRG's own
   convention — the Matsubara frequency is the LAST integration variable, which is what MakeKernel
   reads for its "MatsubaraEven" option — so a flow declaring
   "IntegrationVariables" -> {"l1", "cos1", "phi", "f0"} needs no extra option to get the check.
   None opts out; an empty variable list yields None.
   Automatic applies only to finite-T ("_fT") integrators: a vacuum flow's last variable is an angle,
   and a kernel even in cos1 would otherwise be stamped `matsubara_even`. *)
ntMatsubaraVar[None, _, _] := None;
ntMatsubaraVar[Automatic, _, {}] := None;
ntMatsubaraVar[Automatic, integrator_, vars_List] :=
  If[!StringQ[integrator] || !StringContainsQ[integrator, "_fT"],
    None,
    With[{v = Last[vars]},
      Which[
        StringQ[v], v,
        AssociationQ[v] && KeyExistsQ[v, "Name"], v["Name"],
        True, ToString[v]]]];
ntMatsubaraVar[v_, _, _] := v;

(* ---- DiFfRG flow directory: an explicit path, or DiFfRG's session flowDir ---- *)
ntFlowDir[dir_String] :=
  dir;

ntFlowDir[Automatic] :=
  Module[{d = DiFfRG`CodeTools`Directory`flowDir},
    If[StringQ[d],
      d,
      (* fall back to DiFfRG's own default computation if the symbol is unset *)
      FileNameJoin[
        {
          If[$Notebooks,
            NotebookDirectory[],
            Directory[]
          ],
          "flows"
        }
      ]
    ]
  ];

(* ============================================================================================
   MakeNTKernelDiFfRG — scaffold + numeric kernel emission for one flow.
   ============================================================================================ *)

Options[MakeNTKernelDiFfRG] =
  {
    "Name" -> Automatic, (* REQUIRED flow name, e.g. "ZA" -> dir flows/ZA, class ZA_kernel *)
    "Integrator" -> Automatic, (* REQUIRED DiFfRG integrator template, e.g. "Integrator_p2_1ang" *)
(* Name of the Matsubara-frequency symbol, for a finite-T flow. NumTracer PROVES whether the kernel
   is even in it and, if so, emits DiFfRG's `matsubara_even` trait (one evaluation per mode instead
   of two). Automatic = Last["IntegrationVariables"] (see ntMatsubaraVar); None skips the check; a
   name matching no frame symbol is reported loudly. MakeKernel's own "MatsubaraEven" cannot be used:
   it sees the placeholder `body = 0.`, so it would stamp the trait on every kernel. *)
    "MatsubaraVar" -> Automatic,
(* Passed straight through to MakeNTKernel; see the options there. Together they decide whether the
   `matsubara_finite_extent` trait is emitted, i.e. whether DiFfRG replaces the Gaussian Matsubara rule
   with the exact sum over the modes inside the regulator's support. *)
    "DecayingRegulators" -> Automatic,
    "MatsubaraFiniteExtent" -> Automatic,
(* MakeNTKernel's "ShareInterpolatorIndex" / "HoistLoopConstLookups". Automatic = on exactly when every
   dressing is a DiFfRG interpolator with a split index() (see shareInterpIdx below); False turns the
   pass off; True forces it on. *)
    "ShareInterpolatorIndex" -> Automatic,
    "HoistLoopConstLookups" -> Automatic,
    "IntegrationVariables" -> Automatic, (* REQUIRED, e.g. {"l1","cos1"} *)
    "Parameters" -> Automatic, (* REQUIRED kernelParameterList; also the source of Dressings/DressingType *)
    "ParameterOrder" -> Automatic, (* positional runtime ABI; Automatic follows Parameters order *)
    "Namespace" -> Automatic, (* C++ gen namespace tag; Automatic -> ToLowerCase[Name] *)
    "SymbolDefs" -> <||>, (* derived frame symbols -> closed form, forwarded to MakeNTKernel's "SymbolDefs" (fill) *)
    "AngleDefs" -> {}, (* NumTracer symbolic angle map (spAngles3/4) *)
    "Constant" -> 0., (* loop-independent flat-added term -> constant(p,k,dressings) *)
    (* DiFfRG emission constants — baked, overridable *)
    "Coordinates" -> {"LogarithmicCoordinates1D<double>"},
    "CoordinateArguments" -> {"p"},
(* Regulator for the emitted integrator class. Forwarded verbatim to DiFfRG's MakeKernel, which
   turns the pair into

     <optDef>
     using Regulator = <Regulator><<optName>>;

   inside flows/<Name>/<Name>.hh via DiFfRG`CodeTools`Regulator`getRegulator. "RegulatorOpts" is
   {optName, optDef}: leave it at {"", ""} to instantiate the regulator's own default option struct
   (`Regulator<>`), or pass e.g.

     "Regulator" -> "DiFfRG::RationalExpRegulator",
     "RegulatorOpts" -> {"REGOPTS", "struct REGOPTS { static constexpr int order = 8; ... };"}

   to bake a custom one. *)
    "Regulator" -> "DiFfRG::PolynomialExpRegulator",
    "RegulatorOpts" -> {"", ""},
    "Device" -> "GPU",
    "d" -> 4,
    "AD" -> False,
    (* The value type the kernel is evaluated and integrated in -- "double", "float", or a complex type
       of either. "float" makes the whole runtime path single precision (NumTracer's traces and kernel
       body, DiFfRG's integrator); map()/get() still write double results. *)
    "ComputeType" -> "double",
    "ctype" -> Automatic, (* deprecated spelling of "ComputeType" *)
    "Type" -> "double",
    "Decorator" -> Automatic, (* Automatic -> derived from Device *)
    "FlowDirectory" -> Automatic, (* Automatic -> DiFfRG`CodeTools`flowDir *)
    "GenDirectory" -> Automatic, (* Automatic -> a "gen" sibling of the flow directory *)
    (* Offline by default: emitting a flow writes the generator/probe sources and a numtrace.json
       switch set to 0, and the `numtrace` CMake target (wired in by UpdateNTFlows) compiles and runs
       them when the flows library is built — keeping the notebook out of the C++ build and giving the
       generation `make -j` parallelism across flows. "Offline" -> False (or NT_OFFLINE=0) generates
       inline. *)
    "Offline" -> True,
    (* Forwarded to MakeNTKernel, DELIBERATELY False even though DiFfRG integrators are real-valued
       (a complex verdict-0 kernel fails to compile there): True silently turns verdict 0 into
       Re[Integral[flow]], and truncating the flow equation is the user's call, not the scaffold's. *)
    "RealOutput" -> False,
    (* Forwarded to MakeNTKernel. Finite-density denominators such as l0 + I muq
       are not polynomial in the imaginary unit stand-in used by the symbolic
       real/imaginary split. Keep them factored and project with ntRe/ntIm in
       generated C++ when explicitly requested. *)
    "ComplexRuntimeProjection" -> False,
    (* Forwarded to MakeNTKernel. For real-valued consumers, keep the full
       assembled integrand complex and return its endpoint real part, avoiding
       the symbolic Pure/RePart projection and probe. *)
    "ComplexEndProjection" -> False
  };

MakeNTKernelDiFfRG::ctypedeprecated = "The option \"ctype\" is deprecated; use \"ComputeType\" -> `1` instead.";

MakeNTKernelDiFfRG::noname = "\"Name\" is required (the flow name, e.g. \"ZA\").";

MakeNTKernelDiFfRG::noparams = "\"Parameters\" is required (the DiFfRG kernelParameterList).";

MakeNTKernelDiFfRG::nointeg = "\"Integrator\" and \"IntegrationVariables\" are required.";

MakeNTKernelDiFfRG::mixtype = "Parameters declare more than one interpolator type `1`; emitting the dressing parameters as `const auto&`.";

(* ---- k-only lookup hoisting: patch the DiFfRG-generated wrappers ---------------------------
   When MakeNTKernel hoisted M loop-constant dressing lookups ("HoistCount" of ntMakeNTKernel, see
   "HoistLoopConstLookups" in CodegenKernel.m), the kernel signature carries M trailing `const double&
   nthk<i>` parameters and the kernel class a static host evaluator ntHoisted(k, scalars...,
   dressings...). The integrator forwards its args... verbatim into the kernel call, so the ONLY
   remaining wiring is the wrapper TUs DiFfRG's MakeKernel scaffold already wrote:

     CT_map_*.cc:  return integrator.map(dest, coordinates, <tail>);
       ->          const auto _nth = <name>_kernel<Regulator>::ntHoisted(<tail>);
                   return integrator.map(dest, coordinates, <tail>, _nth[0], ..., _nth[M-1]);
     CT_get.cc:    integrator.get(dest, <coordinates...>, <tail>);   (same treatment)

   <tail> is exactly (k, scalars..., dressings...) — ntHoisted's parameter list by construction. The
   coordinates are the flow's "CoordinateArguments" (none, `p`, `p0, p`, ...), so the get form is
   anchored on `k`. A wrapper left unpatched would only fail when the MODEL compiles, inside a DiFfRG
   integrator header ("too few arguments"), so an unrecognised wrapper aborts here instead.
   `Regulator` is the public alias the scaffold puts in the integrator class; the wrapper body
   resolves <name>_kernel through the class's enclosing namespace (DiFfRG). The patch is
   idempotent by the ntHoisted guard (a re-run of the scaffold rewrites the files fresh anyway). *)

ntPatchHoistWrappers::unpatched = "k-only hoist: the wrapper(s) `1` hold no integrator call of the form `integrator.map(dest, coordinates, k, ...)` / `integrator.get(dest, <coordinates...>, k, ...)`, so they cannot pass the `2` hoisted lookup(s) the kernel now expects. The model would fail to compile inside the DiFfRG integrator. Rerun with \"HoistLoopConstLookups\" -> False, or adapt the patterns in ntPatchHoistWrappers.";

ntPatchHoistWrappers[kernelDir_String, name_String, m_Integer] :=
  Module[{kcls = name <> "_kernel<Regulator>", idxs, files, txt, patched, nPatched = 0, unpatched = {}, hdr, hdrTxt, hdrPatched},
    idxs = StringRiffle[("_nth[" <> ToString[#] <> "]")& /@ Range[0, m - 1], ", "];
(* Make the <name>.hh tuple forwarders take the dressings by reference (`const auto&...t`), so
   ntHoisted gets no interpolator copies. A pure cost fix (DiFfRG interpolators are shallow-copyable);
   a no-op on scaffolds that already emit the reference form. *)
    hdr = FileNameJoin[{kernelDir, name <> ".hh"}];
    If[FileExistsQ[hdr],
      hdrTxt = Import[hdr, "Text"];
      hdrPatched = StringReplace[hdrTxt, "(const auto...t)" -> "(const auto&...t)"];
      Which[
        StringContainsQ[hdrTxt, "(const auto&...t)"],
          Null,(* already by-reference *)
        hdrPatched =!= hdrTxt,
          Export[hdr, hdrPatched, "Text"],
        True,
          Print["[NumTracer] ", name, ": note — tuple-forwarder pattern not found in ",
            hdr, "; the hoisted-lookup host evaluation may receive interpolator copies (a cost, ",
            "not an error)."]]];
    files = FileNames[{"CT_map_*.cc", "CT_get.cc"}, FileNameJoin[{kernelDir, "src"}]];
    Do[
      txt = Import[f, "Text"];
      If[StringContainsQ[txt, "ntHoisted"], Continue[]];
      patched = StringReplace[txt, {
        RegularExpression["return (integrator\\w*)\\.map\\(dest, coordinates,\\s*(.*?)\\);"] :>
          "const auto _nth = " <> kcls <> "::ntHoisted($2);\n  return $1.map(dest, coordinates, $2, " <> idxs <> ");",
        (* $2 = the coordinates, each with its leading whitespace (the scaffold puts two spaces
           before k); re-emitted single-spaced, so `get(dest, p, k, ...)` is spelled as before *)
        RegularExpression["(integrator\\w*)\\.get\\(dest,((?:\\s*\\w+,)*?)\\s*k,\\s*(.*?)\\);"] :>
          "const auto _nth = " <> kcls <> "::ntHoisted(k, $3);\n  $1.get(dest,$2 k, $3, " <> idxs <> ");"}];
      If[patched =!= txt,
        Export[f, patched, "Text"];
        nPatched++,
        AppendTo[unpatched, FileNameTake[f]]],
      {f, files}];
    If[unpatched =!= {},
      Message[ntPatchHoistWrappers::unpatched, StringRiffle[unpatched, ", "], m];
      Abort[]];
    Print["[NumTracer] ", name, ": k-only hoist — ", m, " lookup(s) hoisted; ", nPatched,
      " wrapper TU(s) patched to pass host-evaluated values."]];

(* A plain real scalar parameter (k, T, etaQ, ...), as opposed to a dressing's interpolator type. *)
ntScalarTypeQ[t_] := MemberQ[{"double", "float"}, t];

(* Positional second argument: the constant, mirroring DiFfRG's MakeKernel[kernelExpr, constExpr, ...].
   `constExpr` is a plain expression in p/k and the dressing names (e.g. ZA[p]), NOT an NTKernel.
   If it is also given as "Constant", the positional argument wins (it is spliced ahead of opts). *)

MakeNTKernelDiFfRG[ntk_NTKernel, constExpr_ /; Head[constExpr] =!= Rule && Head[constExpr] =!= NTKernel, opts : OptionsPattern[]] :=
  MakeNTKernelDiFfRG[ntk, "Constant" -> constExpr, opts];

MakeNTKernelDiFfRG[ntk_NTKernel, opts : OptionsPattern[]] :=
  Module[{computeType, name, nsTag, params, parameterOrder, dress, dressTys, dressTy, shareInterpIdx, shareIdxOpt, hoistOpt, hoistCount, diffrgOpts, ntOpts, scalarParams, adParams, device, decor, body, flowDir, genDir, kernelDir, genFile, kernelFile, tracesFile},
    computeType = OptionValue["ComputeType"];
    If[OptionValue["ctype"] =!= Automatic,
      Message[MakeNTKernelDiFfRG::ctypedeprecated, ToString[OptionValue["ctype"], InputForm]];
      computeType = OptionValue["ctype"]];
    name = OptionValue["Name"];
    If[name === Automatic,
      Message[MakeNTKernelDiFfRG::noname];
      Abort[]
    ];
    params = OptionValue["Parameters"];
    If[params === Automatic,
      Message[MakeNTKernelDiFfRG::noparams];
      Abort[]
    ];
    parameterOrder = OptionValue["ParameterOrder"] /. Automatic :> Cases[params, a_?AssociationQ :> a];
    If[OptionValue["Integrator"] === Automatic || OptionValue["IntegrationVariables"] === Automatic,
      Message[MakeNTKernelDiFfRG::nointeg];
      Abort[]
    ];
    nsTag = OptionValue["Namespace"] /. Automatic :> ToLowerCase[name];
    device = OptionValue["Device"];
    (* Kokkos, not raw CUDA: DiFfRG links Kokkos unconditionally and decorates its own kernels this
       way, so the emitted code works on every backend it supports rather than only nvcc. NumTracer's
       ntKokkosDecor normalises any raw __host__ __device__ spelling to the Kokkos macros regardless. *)
    decor =
      OptionValue["Decorator"] /.
        Automatic :>
          If[device === "GPU",
            "static KOKKOS_INLINE_FUNCTION",
            "static inline"
          ];
    body = 0.; (* placeholder scaffold body — NumTracer overwrites kernel.hh below *)
    (* auto-derive the dressing names + their common interpolator type from the non-scalar params *)
    dress = Cases[params, a_?AssociationQ /; !ntScalarTypeQ[a["Type"]] :> a["Name"]];
    dressTys = DeleteDuplicates[Cases[params, a_?AssociationQ /; !ntScalarTypeQ[a["Type"]] :> a["Type"]]];
(* Interpolator index sharing (NumTracer's "ShareInterpolatorIndex"): dressings evaluated at the SAME
   momentum share one coordinate transform (a fp64 log1p, the bulk of a lookup's cost; 1.13-1.26x on
   the YangMills flows). It needs DiFfRG's split index()/at() interface, which all five interpolators
   below provide; the pass is all-or-nothing per kernel, so any other dressing type opts it out. *)
    shareInterpIdx =
      dressTys =!= {} &&
        AllTrue[dressTys,
          StringQ[#] &&
            StringMatchQ[#,
              ("SplineInterpolator1D" | "SplineInterpolator1DStack" | "LinearInterpolator1D" |
                "LinearInterpolator2D" | "LinearInterpolator3D") ~~ "<" ~~ ___] &];
    shareIdxOpt = Replace[OptionValue["ShareInterpolatorIndex"], Automatic -> shareInterpIdx];
    hoistOpt = Replace[OptionValue["HoistLoopConstLookups"], Automatic -> shareInterpIdx];
    dressTy =
      Switch[Length[dressTys],
        0,
          Automatic,
        1,
          First[dressTys],
        _,
(* Mixed interpolator types (e.g. 1-D momentum grids plus a 3-D vertex grid) are passed through PER
   PARAMETER. `const auto&` would also bind them, but makes the kernel an abbreviated function
   template inside a class template, which crashes nvcc's front end (internal error in
   check_name_hiding_by_template_parameters). *)
          Association[Cases[params, a_?AssociationQ /; !ntScalarTypeQ[a["Type"]] :>
            (ntParamName[a["Name"]] -> a["Type"])]]
      ];
(* scalar "double" params (etaPiL, d1V, rhoL, ...) beyond k/p: these are forwarded by DiFfRG's
   integrator between k and the interpolators, so NumTracer must declare them in the kernel /
   constant signature (k and p are already supplied through the NumTrace Args). *)
    scalarParams =
      Cases[
        params,
        a_?AssociationQ /;
            ntScalarTypeQ[a["Type"]] && !MemberQ[{"k", "p"}, ntParamName[a["Name"]]] :> a["Name"]
      ];
(* AD-flagged scalar params (e.g. d1V, d2V for the FE-potential flows): the DiFfRG integrator has an
   autodiff twin (integrator_AD) that forwards these as autodiff::real, so the kernel must declare them
   as `const auto&`, not `const double&`. Thread their names to MakeNTKernel so it types them auto. *)
    adParams =
      Cases[
        params,
        a_?AssociationQ /;
            ntScalarTypeQ[a["Type"]] && TrueQ[a["AD"]] && !MemberQ[{"k", "p"}, ntParamName[a["Name"]]] :> a["Name"]
      ];
    (* resolve paths *)
    flowDir = ntFlowDir[OptionValue["FlowDirectory"]];
    (* create flows/ up front so the DiFfRG scaffold and the kernel/manifest writes have a target *)
    If[!DirectoryQ[flowDir],
      CreateDirectory[flowDir, CreateIntermediateDirectories -> True]
    ];
    (* gen/ as a sibling of flows/, by string surgery: ParentDirectory returns UNEVALUATED for a
       non-existent dir. FileNameDrop tolerates a trailing slash and preserves the filesystem root. *)
    genDir = OptionValue["GenDirectory"] /. Automatic :> FileNameJoin[{FileNameDrop[flowDir, -1], "gen"}];
    If[!DirectoryQ[genDir],
      CreateDirectory[genDir, CreateIntermediateDirectories -> True]
    ];
    kernelDir = FileNameJoin[{flowDir, name}];
    If[!DirectoryQ[kernelDir],
      CreateDirectory[kernelDir, CreateIntermediateDirectories -> True]
    ];
    genFile = FileNameJoin[{genDir, "gen_" <> nsTag <> "_num.cpp"}];
    kernelFile = FileNameJoin[{kernelDir, "kernel.hh"}];
    tracesFile = FileNameJoin[{kernelDir, "kernels.hh"}];
(* (1) DiFfRG scaffold FIRST: lays down flows/<name>/ incl. the integrator TUs + a placeholder kernel.hh.
   Its per-file Prints are captured and replaced by one NumTracer line (see ntReportDiFfRG).

   flowDir is BLOCKED across the call: MakeKernel has no output-directory option and writes to
   DiFfRG`CodeTools`Directory`flowDir, so otherwise "FlowDirectory" would redirect only step (2)
   while the `body = 0.` placeholder overwrote the same-named flow in the session's flow tree. *)
    diffrgOpts = {
      "Name" -> name,
      "Integrator" -> OptionValue["Integrator"],
      "d" -> OptionValue["d"],
      "AD" -> OptionValue["AD"],
      "ComputeType" -> computeType,
      "Device" -> device,
      "Type" -> OptionValue["Type"],
      "Parameters" -> params,
      "IntegrationVariables" -> OptionValue["IntegrationVariables"],
      "Coordinates" -> OptionValue["Coordinates"],
      "CoordinateArguments" -> OptionValue["CoordinateArguments"],
      "Regulator" -> OptionValue["Regulator"],
      "RegulatorOpts" -> OptionValue["RegulatorOpts"]};
    Internal`InheritedBlock[{DiFfRG`CodeTools`Directory`flowDir},
      DiFfRG`CodeTools`Directory`flowDir = flowDir;
      ntReportDiFfRG[name, flowDir, Last @ ntCapturePrint[DiFfRG`CodeTools`MakeKernel`MakeKernel[body, Sequence @@ diffrgOpts]]]];
(* (2) NumTracer overwrites kernel.hh + writes kernels.hh with the real, numerically-traced kernel.

   Guarded: if step (2) does not complete, step (1)'s placeholder kernel.hh (body `0.`) would stay
   on disk as valid C++ contributing a silent ZERO. So on any non-completion it is replaced by a
   `#error` naming the flow. The manifest is not written either (MakeNTKernel writes it last). *)
    hoistCount = 0;
    ntOpts = {
      "ComputeType" -> computeType,
      "Name" -> name <> "_kernel",
      "Namespace" -> nsTag,
      "AngleDefs" -> OptionValue["AngleDefs"],
      "SymbolDefs" -> OptionValue["SymbolDefs"],
      "Decorator" -> decor,
      "DeviceTarget" -> (device === "GPU"),
      "Dressings" -> dress,
      "DressingType" -> dressTy,
      "ShareInterpolatorIndex" -> shareIdxOpt,
      "HoistLoopConstLookups" -> hoistOpt,
      "RealOutput" -> OptionValue["RealOutput"],
      "ComplexRuntimeProjection" -> OptionValue["ComplexRuntimeProjection"],
      "ComplexEndProjection" -> OptionValue["ComplexEndProjection"],
      "ScalarParams" -> scalarParams,
      "ADParams" -> adParams,
      "ParameterOrder" -> parameterOrder,
      "Constant" -> OptionValue["Constant"],
      "Offline" -> OptionValue["Offline"],
      "CoordinateArgs" -> OptionValue["CoordinateArguments"],
      "MatsubaraVar" -> ntMatsubaraVar[OptionValue["MatsubaraVar"], OptionValue["Integrator"], OptionValue["IntegrationVariables"]],
      "DecayingRegulators" -> OptionValue["DecayingRegulators"],
      "MatsubaraFiniteExtent" -> OptionValue["MatsubaraFiniteExtent"],
      "RuntimeInclude" -> None,
      "ExtraIncludes" -> {"DiFfRG/physics/interpolation.hh", "DiFfRG/physics/physics.hh"},
      "KernelNamespace" -> "DiFfRG",
      "SupportNamespace" -> "DiFfRG",
      "RegulatorTemplate" -> True,
      "RegulatorAlias" -> True};
    If[TrueQ @ CheckAbort[
         hoistCount = ntMakeNTKernel[ntk, genFile, kernelFile, tracesFile, Sequence @@ ntOpts]["HoistCount"];
         True,
         False],
      Null,
      Export[kernelFile, "#error NumTracer generation for flow \"" <> name <> "\" did not complete; this kernel.hh is the DiFfRG placeholder (body 0.), not a traced kernel. Re-run the generation and fix the reported error.\n", "Text"];
      Print["[NumTracer] ", name, ": generation ABORTED — ", kernelFile, " poisoned with #error so the build cannot silently use the zero placeholder."];
      Abort[]];
(* (3) wire the hoisted k-only lookups through the scaffold's wrapper TUs (no-op when none). *)
    If[hoistCount > 0,
      ntPatchHoistWrappers[kernelDir, name, hoistCount]];
    kernelFile
  ];

(* ============================================================================================
   NTFlowGate — per-flow selection, so one derivation file can be sliced across processes.

   The Wolfram layer is serial, so the wall-clock lever for a big derivation file is to run it several
   times over disjoint flow subsets. Aggregation is filesystem-based (UpdateFlows rebuilds
   flows/CMakeLists.txt from the flow directory), so a single UpdateNTFlows afterwards sees them all.
   The gate must wrap the WHOLE flow, FunKit derivative and NumTrace included (the front end is a
   large share of a flow's cost), hence HoldRest.
   ============================================================================================ *)

(* The selection, re-read on every call (delayed: an immediate `=` would latch at Get[] time and
   ignore a later SetEnvironment). All = no selection = emit everything.
   With, not Module, so the gate does not bump $ModuleNumber (which would rename every later Module
   local relative to an ungated run). *)
$ntFlowSelection :=
  With[{s = Environment["NT_FLOWS"]},
    If[! StringQ[s] || StringTrim[s] === "",
      All,
      DeleteCases[StringTrim /@ StringSplit[s, ","], ""]]];

(* Names actually seen by a gate this run — the other half of the typo check. A selection entry that
   matches no gate produces no output and no error, which is indistinguishable from a run that had
   nothing to do; NTFlowSelectionReport turns that into a message. *)
$ntFlowGateSeen = {};

NTFlowSelectedQ[name_String] :=
  With[{sel = $ntFlowSelection},
    sel === All || MemberQ[sel, name]];

(* Is this run a SLICE (a selection is in force) rather than the whole file? A derivation file uses
   it to decide whether to run the final UpdateNTFlows itself: a slice must not, because the flows
   its siblings are still writing are not on disk yet, and CMakeLists.txt is built from what IS. *)
NTFlowSliceQ[] := $ntFlowSelection =!= All;

SetAttributes[NTFlowGate, HoldRest];

NTFlowGate[name_String, body_] := NTFlowGate[{name}, body];

(* List form, for a block whose FRONT END is shared by several flows (etaPiL/zPiL take one FunKit
   derivative and emit two kernels from it). Gate the block on the whole set and each emit on its own
   name: selecting one of the pair then pays the shared derivative once and emits only what was asked
   for, instead of forcing an all-or-nothing split of a block that is deliberately fused. *)
NTFlowGate[names : {__String}, body_] :=
  ($ntFlowGateSeen = Union[$ntFlowGateSeen, names];
   If[AnyTrue[names, NTFlowSelectedQ],
     body,
     Print["[NumTracer] " <> StringRiffle[names, "/"] <> ": skipped (NT_FLOWS)"];
     Null]);

NTFlowSelectionReport::unmatched = "NT_FLOWS asked for `1`, which no NTFlowGate in this run matched — check the spelling against the gate names in the derivation file. Nothing was emitted for it.";

NTFlowSelectionReport[] :=
  With[{sel = $ntFlowSelection},
    If[sel === All,
      {},
      With[{missing = Complement[sel, $ntFlowGateSeen]},
        If[missing =!= {}, Message[NTFlowSelectionReport::unmatched, missing]];
        missing]]];

(* ============================================================================================
   UpdateNTFlows — DiFfRG UpdateFlows + idempotent NumTracer CMake patch (atomic).
   ============================================================================================ *)

Options[UpdateNTFlows] = {"FlowDirectory" -> Automatic, "NumTracerHints" -> "~/.local/share/NumTracer", "UnityBuild" -> False};

UpdateNTFlows::nocmake = "Expected the flows CMakeLists at `1` (generate a flow first).";

UpdateNTFlows::patchfail = "Patched `1` but it does not reference NumTracer / numtracer_add_numtrace — the DiFfRG CMake template changed; the find_package / link / numtrace strings in DiFfRG_compat.m need updating.";

UpdateNTFlows[name_String, opts : OptionsPattern[]] :=
  Module[{flowDir, f, txt},
    flowDir = ntFlowDir[OptionValue["FlowDirectory"]];
(* (1) DiFfRG aggregation — regenerates flows/CMakeLists.txt from its template (wipes any prior patch).
   Chatter captured; step (4) prints the NumTracer-native summary once the patch is back on. *)
    With[{chatter = Last @ ntCapturePrint[DiFfRG`CodeTools`UpdateFlows[name]]},
      ntLog["[DiFfRG] ", #]& /@ chatter];
    f = FileNameJoin[{flowDir, "CMakeLists.txt"}];
    If[!FileExistsQ[f],
      Message[UpdateNTFlows::nocmake, f];
      Abort[]
    ];
    txt = Import[f, "Text"];
    (* (2) idempotent patch: find_package(NumTracer) + link NumTracer::* + UNITY_BUILD OFF *)
    If[!StringContainsQ[txt, "find_package(NumTracer"],
      txt = StringReplace[txt, "add_library(" <> name <> " STATIC ${" <> name <> "_SOURCES})" -> "find_package(NumTracer REQUIRED HINTS " <> OptionValue["NumTracerHints"] <> ")\n\n" <> "add_library(" <> name <> " STATIC ${" <> name <> "_SOURCES})"];
      txt = StringReplace[txt, "target_link_libraries(" <> name <> " DiFfRG::DiFfRG " <> name <> "_nowarn)" -> "target_link_libraries(" <> name <> " DiFfRG::DiFfRG " <> name <> "_nowarn NumTracer::NumTracer)"]
    ];
    If[!TrueQ[OptionValue["UnityBuild"]],
      txt = StringReplace[txt, "UNITY_BUILD ON" -> "UNITY_BUILD OFF"]
    ];
    (* (2b) the `numtrace` target: reads every flows/<Name>/numtrace.json, builds and runs the
       generator (and probe) of each flow whose switch is still 0, and makes the flows library depend
       on the lot. A no-op once every switch is 1. numtracer_add_numtrace comes from
       NumTracerNumtrace.cmake, pulled in by the find_package(NumTracer) added above. *)
    If[!StringContainsQ[txt, "numtracer_add_numtrace"],
      txt = StringReplace[txt, "target_link_libraries(" <> name <> " DiFfRG::DiFfRG " <> name <> "_nowarn NumTracer::NumTracer)" -> "target_link_libraries(" <> name <> " DiFfRG::DiFfRG " <> name <> "_nowarn NumTracer::NumTracer)\n\n" <> "numtracer_add_numtrace(" <> name <> " ${CMAKE_CURRENT_SOURCE_DIR})"]
    ];
    Export[f, txt, "Text"];
    (* (3) loud failure on template drift — StringReplace no-ops silently on a mismatch *)
    If[!StringContainsQ[txt, "NumTracer::NumTracer"] || !StringContainsQ[txt, "numtracer_add_numtrace"],
      Message[UpdateNTFlows::patchfail, f];
      Abort[]
    ];
    (* (4) NumTracer-native confirmation, in place of DiFfRG's captured chatter *)
    Print[
      "[NumTracer] " <> name <> ": flows/CMakeLists.txt regenerated + NumTracer patch applied " <> "(find_package + NumTracer::* link" <>
        If[!TrueQ[OptionValue["UnityBuild"]],
          " + UNITY_BUILD OFF",
          ""
        ] <> ")"
    ];
    f
  ];
