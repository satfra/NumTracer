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
   DiFfRG's generators Print a line per emitted file plus a "Please run UpdateFlows[]" nudge. Under
   NumTracer that nudge is actively WRONG: bare UpdateFlows regenerates flows/CMakeLists.txt from
   DiFfRG's template and drops the NumTracer patch (find_package + link libs + UNITY_BUILD OFF) —
   UpdateNTFlows exists to do both atomically. So capture the chatter, print one NumTracer-native line
   instead, and keep the raw lines behind $NumTracerVerbose for debugging.
   Internal`InheritedBlock scopes the Print override, so it is restored even if DiFfRG aborts. *)

SetAttributes[ntCapturePrint, HoldFirst];

ntCapturePrint[expr_] :=
  Module[{lines = {}, res},
    res =
      Internal`InheritedBlock[
        {Print}
        ,
        Unprotect[Print];
        Print[args___] := AppendTo[lines, StringJoin[ToString /@ {args}]];
        expr
      ];
    {res, lines}
  ];

(* Does flows/CMakeLists.txt still need an UpdateNTFlows pass for this flow? DiFfRG nudges every time;
   we only nudge when it is actually true — the file is missing, has lost the NumTracer patch, or does
   not yet mention this flow (i.e. a newly added one). *)

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

(* ---- DiFfRG flow directory (public delayed symbol, trailing slash) -------------------------- *)

(* Resolve the "MatsubaraVar" option to a symbol name (or None). Automatic follows DiFfRG's own
   convention — the Matsubara frequency is the LAST integration variable, which is what MakeKernel
   reads for its "MatsubaraEven" option — so a flow declaring
   "IntegrationVariables" -> {"l1", "cos1", "phi", "f0"} needs no extra option to get the check.
   A flow whose frame names the frequency something else passes it explicitly; None opts out.
   An empty variable list yields None rather than a Last[] error, so vacuum flows are unaffected.

   Automatic is gated on the integrator actually being a finite-T one. Without that gate a VACUUM
   flow's last integration variable is an angle ("cos1"), which is a perfectly good frame symbol —
   so the check would run, and a kernel that happens to be even in cos1 would be stamped
   `matsubara_even`. Inert today (no vacuum integrator reads the trait) but exactly the kind of
   latent mis-binding that surfaces later as wrong physics. *)
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

ntFlowDir[dir_String] :=
  dir;

ntFlowDir[Automatic] :=
  Module[{d = DiFfRG`CodeTools`Directory`flowDir},
    If[StringQ[d],
      d
      ,
      (* fall back to DiFfRG's own default computation if the symbol is unset *)
      FileNameJoin[
        {
          If[$Notebooks,
            NotebookDirectory[]
            ,
            Directory[]
          ]
          ,
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
    "Name" -> Automatic
    ,(* REQUIRED flow name, e.g. "ZA" -> dir flows/ZA, class ZA_kernel *)
    "Integrator" -> Automatic
    ,(* REQUIRED DiFfRG integrator template, e.g. "Integrator_p2_1ang" *)
(* Name of the Matsubara-frequency symbol, for a finite-T flow. It asks NumTracer to PROVE whether
   the kernel is even in that frequency and, if so, emit DiFfRG's `matsubara_even` trait: the
   integrator then evaluates the kernel once per Matsubara mode instead of twice, halving both that
   work and the number of inlined kernel copies per launch.

   Automatic = Last["IntegrationVariables"], the same convention DiFfRG's own MakeKernel uses for
   its "MatsubaraEven" option. Give an explicit symbol name to override, or None to skip the check
   entirely. If the name matches no symbol in the flow's frame, NumTracer says so loudly and emits
   no trait rather than quietly skipping the optimisation.

   This cannot be delegated to MakeKernel's "MatsubaraEven" option: on the numeric path MakeKernel
   is handed the placeholder `body = 0.` (NumTracer overwrites kernel.hh afterwards), so its
   symbolic check is trivially satisfied for EVERY flow and would stamp the trait on kernels that
   are not even -- silently dropping the odd half of the Matsubara sum. The verdict has to come
   from the polynomials, which is where NumTracer proves it. *)
    "MatsubaraVar" -> Automatic,
(* Passed straight through to MakeNTKernel; see the options there. Together they decide whether the
   `matsubara_finite_extent` trait is emitted, i.e. whether DiFfRG replaces the Gaussian Matsubara rule
   with the exact sum over the modes inside the regulator's support. *)
    "DecayingRegulators" -> Automatic,
    "MatsubaraFiniteExtent" -> Automatic,
    "IntegrationVariables" -> Automatic
    ,(* REQUIRED, e.g. {"l1","cos1"} *)
    "Parameters" -> Automatic
    ,(* REQUIRED kernelParameterList; also the source of Dressings/DressingType *)
    "ParameterOrder" -> Automatic
    ,(* positional runtime ABI; Automatic follows Parameters order *)
    "Namespace" -> Automatic
    ,(* C++ gen namespace tag; Automatic -> ToLowerCase[Name] *)
    "AngleDefs" -> {}
    ,(* NumTracer symbolic angle map (spAngles3/4) *)
    "Constant" -> 0.
    ,(* loop-independent flat-added term -> constant(p,k,dressings) *)
    (* DiFfRG emission constants — baked, overridable *)
    "Coordinates" -> {"LogarithmicCoordinates1D<double>"}
    ,
    "CoordinateArguments" -> {"p"}
    ,
(* Regulator for the emitted integrator class. Forwarded verbatim to DiFfRG's MakeKernel, which
   turns the pair into

     <optDef>
     using Regulator = <Regulator><<optName>>;

   inside flows/<Name>/<Name>.hh via DiFfRG`CodeTools`Regulator`getRegulator. "RegulatorOpts" is
   {optName, optDef}: leave it at {"", ""} to instantiate the regulator's own default option struct
   (`Regulator<>`), or pass e.g.

     "Regulator" -> "DiFfRG::RationalExpRegulator",
     "RegulatorOpts" -> {"REGOPTS", "struct REGOPTS { static constexpr int order = 8; ... };"}

   to bake a custom one. Without this passthrough every NumTracer-emitted kernel silently inherited
   MakeKernel's default, so an app could not choose its own regulator at all. *)
    "Regulator" -> "DiFfRG::PolynomialExpRegulator"
    ,
    "RegulatorOpts" -> {"", ""}
    ,
    "Device" -> "GPU"
    ,
    "d" -> 4
    ,
    "AD" -> False
    ,
    "ctype" -> "double"
    ,
    "Type" -> "double"
    ,
(* Fuse every trace into one shared CSE program (`trace_all`) instead of emitting one function per
   trace. Off by default, as in MakeNTKernel. Note the design comment in Codegen.m: GlobalCollect
   (on by default) is said to subsume it -- but nothing in the code enforces that, so the two are
   independently settable and the combination is worth measuring rather than assuming. *)
    "CrossTraceCSE" -> False,
    "Decorator" -> Automatic
    ,(* Automatic -> derived from Device *)
    "FlowDirectory" -> Automatic
    ,(* Automatic -> DiFfRG`CodeTools`flowDir *)
    "GenDirectory" -> Automatic
    ,(* Automatic -> a "gen" sibling of the flow directory *)
    (* Offline by default: emitting a flow writes the generator/probe sources and a numtrace.json
       switch set to 0, and the `numtrace` CMake target (wired in by UpdateNTFlows) compiles and runs
       them when the flows library is built — keeping the notebook out of the C++ build and giving the
       generation `make -j` parallelism across every flow at once. Pass "Offline" -> False (or set
       NT_OFFLINE=0) to generate inline as before. *)
    "Offline" -> True,
    (* Forwarded to MakeNTKernel, and DELIBERATELY still False here. A DiFfRG flow consumer really is
       real-valued — the integrator instantiates KERNEL::kernel(...) into a `double` accumulator, so a
       verdict-0 flow does not merely waste a lowered body, it fails to compile at all (a concept wall
       in QuadratureIntegrator, naming the integrator rather than the flow). It is therefore tempting
       to default this to True here. Do not: with the complex body gone, verdict 0 silently becomes
       Re[Integral[flow]], and choosing to truncate the flow equation belongs to whoever is doing the
       physics, not to the scaffolding. Set it explicitly when that is what you mean. *)
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

MakeNTKernelDiFfRG::noname = "\"Name\" is required (the flow name, e.g. \"ZA\").";

MakeNTKernelDiFfRG::noparams = "\"Parameters\" is required (the DiFfRG kernelParameterList).";

MakeNTKernelDiFfRG::nointeg = "\"Integrator\" and \"IntegrationVariables\" are required.";

MakeNTKernelDiFfRG::mixtype = "Parameters declare more than one interpolator type `1`; emitting the dressing parameters as `const auto&`.";

(* ---- k-only lookup hoisting: patch the DiFfRG-generated wrappers ---------------------------
   When MakeNTKernel hoisted M loop-constant dressing lookups ($ntLastHoistCount, see
   "HoistLoopConstLookups" in Codegen.m), the kernel signature carries M trailing `const double&
   nthk<i>` parameters and the kernel class a static host evaluator ntHoisted(k, scalars...,
   dressings...). The integrator forwards its args... verbatim into the kernel call, so the ONLY
   remaining wiring is the wrapper TUs DiFfRG's MakeKernel scaffold already wrote:

     CT_map_*.cc:  return integrator.map(dest, coordinates, <tail>);
       ->          const auto _nth = <name>_kernel<Regulator>::ntHoisted(<tail>);
                   return integrator.map(dest, coordinates, <tail>, _nth[0], ..., _nth[M-1]);
     CT_get.cc:    integrator.get(dest, p, <tail>);              (same treatment)

   <tail> is exactly (k, scalars..., dressings...) — ntHoisted's parameter list by construction.
   `Regulator` is the public alias the scaffold puts in the integrator class; the wrapper body
   resolves <name>_kernel through the class's enclosing namespace (DiFfRG). The patch is
   idempotent by the ntHoisted guard (a re-run of the scaffold rewrites the files fresh anyway). *)

ntPatchHoistWrappers[kernelDir_String, name_String, m_Integer] :=
  Module[{kcls = name <> "_kernel<Regulator>", idxs, files, txt, patched, nPatched = 0, hdr, hdrTxt, hdrPatched},
    idxs = StringRiffle[("_nth[" <> ToString[#] <> "]")& /@ Range[0, m - 1], ", "];
(* Older MakeKernel scaffolds unpacked the dressing tuple BY VALUE in the <name>.hh forwarders
   (`const auto...t`), handing ntHoisted a full copy of every interpolator per call. Current
   DiFfRG interpolators are shallow-copyable (Kokkos views, host AND device buffers survive a
   copy), so this is no longer a correctness matter — but the copies are pure cost, and current
   scaffolds already emit `const auto&...t`, in which case this rewrite is a no-op. *)
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
        RegularExpression["(integrator\\w*)\\.get\\(dest, p,\\s*(.*?)\\);"] :>
          "const auto _nth = " <> kcls <> "::ntHoisted($2);\n  $1.get(dest, p, $2, " <> idxs <> ");"}];
      If[patched =!= txt,
        Export[f, patched, "Text"];
        nPatched++],
      {f, files}];
    Print["[NumTracer] ", name, ": k-only hoist — ", m, " lookup(s) hoisted; ", nPatched,
      " wrapper TU(s) patched to pass host-evaluated values."]];

(* Positional second argument: the constant, mirroring DiFfRG's MakeKernel[kernelExpr, constExpr, ...].
   `constExpr` is a plain Mathematica expression in p/k and the dressing names (e.g. ZA[p]), NOT an
   NTKernel — the integrand still comes from `ntk`. Pass the constant EITHER positionally OR via the
   "Constant" option; if both are given the positional argument wins (it is spliced ahead of opts). *)

MakeNTKernelDiFfRG[ntk_NTKernel, constExpr_ /; Head[constExpr] =!= Rule && Head[constExpr] =!= NTKernel, opts : OptionsPattern[]] :=
  MakeNTKernelDiFfRG[ntk, "Constant" -> constExpr, opts];

MakeNTKernelDiFfRG[ntk_NTKernel, opts : OptionsPattern[]] :=
  Module[{name, nsTag, params, parameterOrder, dress, dressTys, dressTy, shareInterpIdx, scalarParams, adParams, device, decor, body, flowDir, genDir, kernelDir, genFile, kernelFile, tracesFile},
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
            "static KOKKOS_INLINE_FUNCTION"
            ,
            "static inline"
          ];
    body = 0.; (* placeholder scaffold body — NumTracer overwrites kernel.hh below *)
    (* auto-derive the dressing names + their common interpolator type from the non-"double" params *)
    dress = Cases[params, a_?AssociationQ /; a["Type"] =!= "double" :> a["Name"]];
    dressTys = DeleteDuplicates[Cases[params, a_?AssociationQ /; a["Type"] =!= "double" :> a["Type"]]];
(* Interpolator index sharing (NumTracer's "ShareInterpolatorIndex"): a flow evaluates several
   dressings at the SAME momentum, and each lookup otherwise repeats the coordinate transform --
   a fp64 log1p for a logarithmic axis, ~200 of the ~210 fp64 instructions a lookup costs. The
   compiler cannot share it, because every interpolator owns its own `coordinates` members. Paying
   it once measured 1.13-1.26x on the YangMills flow set, results bit-identical
   (numtracer/gpubench/FINDINGS.md).

   This is a DiFfRG-side capability, which is why NumTracer keeps the pass opt-in and it is enabled
   HERE. All five DiFfRG interpolators are split: the 1-D ones return a scalar index, the 2-D/3-D
   and stack ones return a device::array of indices, so the emitted `const auto _ixN = h.index(...)`
   is uniform and the rewrite needs no per-arity special case. The pass is still all-or-nothing per
   kernel, so any dressing handle NOT in this list opts the whole kernel out. *)
    shareInterpIdx =
      dressTys =!= {} &&
        AllTrue[dressTys,
          StringQ[#] &&
            StringMatchQ[#,
              ("SplineInterpolator1D" | "SplineInterpolator1DStack" | "LinearInterpolator1D" |
                "LinearInterpolator2D" | "LinearInterpolator3D") ~~ "<" ~~ ___] &];
    dressTy =
      Switch[Length[dressTys],
        0,
          Automatic
        ,
        1,
          First[dressTys]
        ,
        _,
(* Mixed interpolator types are legitimate: a flow may read most dressings off 1-D momentum grids
   and one or two off a 3-D vertex grid (S0,S1,SPhi). Taking the first would declare the 3-D grid as
   a 1-D spline and the kernel would fail to compile at the first three-argument call.

   Pass the types through PER PARAMETER rather than collapsing them. `const auto&` would also bind
   all of them, but it turns every kernel into an abbreviated function template, and nvcc's front end
   cannot handle that inside a class template -- it dies with
       internal error: assertion failed at: "symbol_ref.c", line 1629
       in check_name_hiding_by_template_parameters
   (reproducible in ~20 lines: a class template whose static member takes `const auto&` params).
   The parameter list already carries each declared type, so there is nothing to infer. *)
          Association[Cases[params, a_?AssociationQ /; a["Type"] =!= "double" :>
            (If[StringQ[a["Name"]], a["Name"], SymbolName[a["Name"]]] -> a["Type"])]]
      ];
(* scalar "double" params (etaPiL, d1V, rhoL, ...) beyond k/p: these are forwarded by DiFfRG's
   integrator between k and the interpolators, so NumTracer must declare them in the kernel /
   constant signature (k and p are already supplied through the NumTrace Args). *)
    scalarParams =
      Cases[
        params
        ,
        a_?AssociationQ /;
            a["Type"] === "double" &&
              !MemberQ[
                {"k", "p"}
                ,
                If[StringQ[a["Name"]],
                  a["Name"]
                  ,
                  SymbolName[a["Name"]]
                ]
              ] :> a["Name"]
      ];
(* AD-flagged scalar params (e.g. d1V, d2V for the FE-potential flows): the DiFfRG integrator has an
   autodiff twin (integrator_AD) that forwards these as autodiff::real, so the kernel must declare them
   as `const auto&`, not `const double&`. Thread their names to MakeNTKernel so it types them auto. *)
    adParams =
      Cases[
        params
        ,
        a_?AssociationQ /;
            a["Type"] === "double" && TrueQ[a["AD"]] &&
              !MemberQ[
                {"k", "p"}
                ,
                If[StringQ[a["Name"]],
                  a["Name"]
                  ,
                  SymbolName[a["Name"]]
                ]
              ] :> a["Name"]
      ];
    (* resolve paths *)
    flowDir = ntFlowDir[OptionValue["FlowDirectory"]];
    (* create the flows/ directory up front if it does not exist yet. This must precede the gen dir
       computation: WL14.3's ParentDirectory[dir] returns UNEVALUATED for a non-existent dir (rather
       than doing pure-string path math), which used to poison genDir — and every genFile/tracesFile
       derived from it — with a held expression, so the whole emission failed on a fresh checkout. We
       sidestep ParentDirectory entirely (pure-string parent below), but still ensure flows/ exists so
       the DiFfRG scaffold and the kernel/manifest writes below have somewhere to land. *)
    If[!DirectoryQ[flowDir],
      CreateDirectory[flowDir, CreateIntermediateDirectories -> True]
    ];
    (* gen/ as a sibling of flows/, computed by string surgery (no ParentDirectory — see above):
       FileNameDrop[..,-1] drops the last path segment, tolerating a trailing slash and preserving the
       filesystem root (a plain split+DeleteCases[""] would strip the leading root marker). *)
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

   flowDir is BLOCKED across the call. DiFfRG's MakeKernel has no output-directory option — it writes to
   FileNameJoin[DiFfRG`CodeTools`Directory`flowDir, name] — so without this our "FlowDirectory" option
   would redirect only NumTracer's own writes (step 2) while the scaffold silently kept landing in the
   session-global flow directory. That is a live footgun, not a hypothetical: emitting one flow into a
   scratch directory then overwrote the REAL committed kernel.hh of the same-named flow in the default
   tree with the `body = 0.` placeholder — valid C++ that contributes zero to the flow, with nothing
   saying so. Scoping it makes "FlowDirectory" mean what it says for both halves of the emission. *)
    Internal`InheritedBlock[{DiFfRG`CodeTools`Directory`flowDir},
      DiFfRG`CodeTools`Directory`flowDir = flowDir;
      ntReportDiFfRG[name, flowDir, Last @ ntCapturePrint[DiFfRG`CodeTools`MakeKernel`MakeKernel[body, "Name" -> name, "Integrator" -> OptionValue["Integrator"], "d" -> OptionValue["d"], "AD" -> OptionValue["AD"], "ctype" -> OptionValue["ctype"], "Device" -> device, "Type" -> OptionValue["Type"], "Parameters" -> params, "IntegrationVariables" -> OptionValue["IntegrationVariables"], "Coordinates" -> OptionValue["Coordinates"], "CoordinateArguments" -> OptionValue["CoordinateArguments"], "Regulator" -> OptionValue["Regulator"], "RegulatorOpts" -> OptionValue["RegulatorOpts"]]]]];
(* (2) NumTracer overwrites kernel.hh + writes kernels.hh with the real, numerically-traced kernel.

   Guarded, because step (1) has already written a PLACEHOLDER kernel.hh whose body is literally
   `0.`. If step (2) does not complete — a leak detected at the emission chokepoint, a singular Gram,
   an interrupt, an OOM — that placeholder is what stays on disk, and it is valid C++: the flow
   builds cleanly and contributes ZERO to the RG flow, with nothing anywhere saying so. A silent zero
   is the worst failure this pipeline can produce, so on any non-completion we replace the file with
   a `#error`, which turns it into a loud compile failure naming the flow. The manifest is not
   written in that case either (MakeNTKernel writes it last), so `make numtrace` also still owes the
   kernels. *)
    $ntLastHoistCount = 0;
    If[TrueQ @ CheckAbort[
         MakeNTKernel[ntk, genFile, kernelFile, tracesFile, "Name" -> name <> "_kernel", "Namespace" -> nsTag, "AngleDefs" -> OptionValue["AngleDefs"], "Decorator" -> decor, "DeviceTarget" -> (device === "GPU"), "Dressings" -> dress, "DressingType" -> dressTy, "ShareInterpolatorIndex" -> shareInterpIdx, "HoistLoopConstLookups" -> shareInterpIdx, "CrossTraceCSE" -> OptionValue["CrossTraceCSE"], "RealOutput" -> OptionValue["RealOutput"], "ComplexRuntimeProjection" -> OptionValue["ComplexRuntimeProjection"], "ComplexEndProjection" -> OptionValue["ComplexEndProjection"], "ScalarParams" -> scalarParams, "ADParams" -> adParams, "ParameterOrder" -> parameterOrder, "Constant" -> OptionValue["Constant"], "Offline" -> OptionValue["Offline"], "CoordinateArgs" -> OptionValue["CoordinateArguments"], "MatsubaraVar" -> ntMatsubaraVar[OptionValue["MatsubaraVar"], OptionValue["Integrator"], OptionValue["IntegrationVariables"]], "DecayingRegulators" -> OptionValue["DecayingRegulators"], "MatsubaraFiniteExtent" -> OptionValue["MatsubaraFiniteExtent"], "RuntimeInclude" -> None, "ExtraIncludes" -> {"DiFfRG/physics/interpolation.hh", "DiFfRG/physics/physics.hh"}, "KernelNamespace" -> "DiFfRG", "SupportNamespace" -> "DiFfRG", "RegulatorTemplate" -> True, "RegulatorAlias" -> True];
         True,
         False],
      Null,
      Export[kernelFile, "#error NumTracer generation for flow \"" <> name <> "\" did not complete; this kernel.hh is the DiFfRG placeholder (body 0.), not a traced kernel. Re-run the generation and fix the reported error.\n", "Text"];
      Print["[NumTracer] ", name, ": generation ABORTED — ", kernelFile, " poisoned with #error so the build cannot silently use the zero placeholder."];
      Abort[]];
(* (3) wire the hoisted k-only lookups through the scaffold's wrapper TUs (no-op when none).
   $ntLastHoistCount is a Private`-context handoff (the NumTracer` context is Protected). *)
    If[IntegerQ[$ntLastHoistCount] && $ntLastHoistCount > 0,
      ntPatchHoistWrappers[kernelDir, name, $ntLastHoistCount]];
    kernelFile
  ];

(* ============================================================================================
   NTFlowGate — per-flow selection, so one derivation file can be sliced across processes.

   A derivation file emits every flow it defines, top to bottom, in one kernel. That is fine until
   the file is big: with_mesons/QCD.wl spends ~63% of a ~30 min regeneration on nine four-fermion
   channels that differ only in a projector index. Nothing in the Wolfram layer is parallel (the
   generator's own phases and its compile are, but they run later, under the `numtrace` target), so
   the only lever on wall clock is to run the FILE several times over disjoint flow subsets.

   That works because aggregation is filesystem-based: DiFfRG's UpdateFlows rebuilds
   flows/CMakeLists.txt from `FileNames["*", flowDir, 1]`, not from an in-process registry. So the
   workers need share nothing, and a single UpdateNTFlows afterwards — in any process — sees them all.

   The gate must wrap the WHOLE flow, FunKit derivative and NumTrace included. Gating only the
   MakeNTKernelDiFfRG call would skip the emission and still pay the front end, and that is not a
   rounding error: measured on the four-fermion channels it is 35-43 s (of which NumTrace alone is
   18-21 s), against a whole-flow cost of 83-210 s. Hence HoldRest.
   ============================================================================================ *)

(* The selection, re-read on every call (delayed, like every other env-derived flag here — with an
   immediate `=` the value would latch at Get[] time and SetEnvironment from a driver .wls would
   silently do nothing). All = no selection = emit everything, so an ungated file is unaffected. *)
(* With, not Module: this runs once per gate per flow, and a Module would mint a fresh `s$nnn` each
   time, bumping $ModuleNumber and renaming every Module local below the gates relative to an ungated
   run. Nothing observable was traced to that -- the gate is byte-inert either way (verified: seven
   flows regenerate identically to a serial run) -- but a scoping construct on a hot path should not
   churn a global counter for a binding that is never captured. Matches ntEnvFlag in DSL.m. *)
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
          " + UNITY_BUILD OFF"
          ,
          ""
        ] <> ")"
    ];
    f
  ];
