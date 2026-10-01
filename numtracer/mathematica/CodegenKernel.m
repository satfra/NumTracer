(* CodegenKernel.m — the kernel-level driver: mkGenerateKernel (one flow, from an analysed NTKernel
   to generator sources, traces header and kernel header), the kernel-class/header emission helpers,
   and the public MakeNTKernel entry point. Loaded by NumTracer.m via ntLoadPart, in the
   NumTracer`Private` context. *)

(* ---- decorator normalisation: raw CUDA qualifiers -> Kokkos spelling -----------------------
   The decorator is emitted onto every generated function (kernel/constant, regulator wrappers, and
   via the generator's `-d` flag fill/trN/powr). The Kokkos macros expand correctly for every
   backend, including plain `inline` on a host-only build; `__host__ __device__` would hard-code
   CUDA. Longest match first: the inline/__forceinline__ variants must be rewritten before the bare
   rule claims them. *)
ntKokkosDecor[decor_String] := StringReplace[decor, {
  "__host__ __device__ inline"          -> "KOKKOS_INLINE_FUNCTION",
  "__host__ __device__ __forceinline__" -> "KOKKOS_FORCEINLINE_FUNCTION",
  "__host__ __device__"                 -> "KOKKOS_FUNCTION"}];
ntKokkosDecor[d_] := d;

(* The Kokkos macros are not in scope in the standalone build-time programs (e.g. the probe compiles
   the traces header with a bare g++), so this preamble neutralises them there. *)
ntKokkosStubDefs = "#define KOKKOS_INLINE_FUNCTION inline\n#define KOKKOS_FORCEINLINE_FUNCTION \
inline\n#define KOKKOS_FUNCTION\n#define __host__\n#define __device__\n";

(* regulator wrappers, prefixed with the user-chosen decorator (default plain "static inline";
   pass e.g. "Decorator" -> "static KOKKOS_INLINE_FUNCTION" for device-callable kernels).
   Emitted ONLY under "RegulatorTemplate" -> True (the fRG/DiFfRG shape); see ntKernelClass. *)
privDefs[decor_] :=
  StringRiffle[(decor <> " auto " <> # <> "(const auto &k2, const auto &p2) { return REG::" <> # <> "(k2, p2); }")& /@ {"RB", "RF", "RBdot", "RFdot", "dq2RB", "dq2RF"}, "\n"];

(* Real/imag accessors for the kernel class and the probe. They may be applied to an autodiff-typed
   integrand ("ComplexEndProjection"), so they are type-PRESERVING: a hard `-> double` would drop the
   derivative. The fallback is an UNQUALIFIED real(z)/imag(z) so a DiFfRG consumer's DiFfRG::real/imag
   is found by ordinary lookup. Trailing return types because kernel() calls these before their
   in-class definition, where g++ rejects a deduced `auto`. *)
ntReImAccessors[decor_] :=
  StringRiffle[{
    decor <> " " <> $ntRealT <> " ntRe(" <> $ntRealT <> " x) { return x; }",
    decor <> " " <> $ntRealT <> " ntIm(" <> $ntRealT <> ") { return " <> ntZeroLit[] <> "; }",
    "template <class T> " <> decor <> " auto ntRe(const T &z) -> decltype(z.real()) { return z.real(); }",
    "template <class T> " <> decor <> " auto ntRe(const T &z) -> decltype(real(z)) requires (!requires { z.real(); }) { return real(z); }",
    "template <class T> " <> decor <> " auto ntIm(const T &z) -> decltype(z.imag()) { return z.imag(); }",
    "template <class T> " <> decor <> " auto ntIm(const T &z) -> decltype(imag(z)) requires (!requires { z.imag(); }) { return imag(z); }"}, "\n"];

(* ---- interpolator index sharing --------------------------------------------------------------
   COEN hoists every dressing lookup into `const auto _interpN = <Interp>(<arg>);`, and several
   dressings are often evaluated at the SAME argument. The coordinate transform inside a lookup is
   expensive (an fp64 log1p) and the compiler cannot share it across interpolators, so each such
   group is rewritten to pay it once (1.13-1.26x end to end on YangMills, bit-identical):

     const auto _interp3 = ZA3(<arg>);        ->   const auto _ix0 = ZA3.index(<arg>);
     const auto _interp7 = ZAcbc(<arg>);           const auto _interp3 = ZA3.at(_ix0);
                                                   const auto _interp7 = ZAcbc.at(_ix0);

   Legal only because index() depends solely on the coordinate system (clamping lives in at()) and
   all interpolators of one flow share the consumer's coordinate object — a DiFfRG property, hence
   opt-in via "ShareInterpolatorIndex". Applied per kernel BODY, so a hoist never crosses an #if
   branch. Callees are restricted to `dress`: regulator calls and ntRe(...) have the same
   `const auto _interpN = f(...)` shape and must NOT be rewritten. *)
ntInterpLine[ln_String] :=
  Module[{m = StringCases[ln,
      RegularExpression["^(\\s*)const auto (_interp\\d+) = (\\w+)\\((.*)\\);\\s*$"] :> {"$1", "$2", "$3", "$4"}]},
    If[m === {}, None, First[m]]];

ntShareInterpIndices[text_String, splines_List] :=
  Module[{lines, parsed, counts, shared, idxName, emitted, out, nrw = 0},
    If[splines === {}, Return[text]];
    lines = StringSplit[text, "\n", All];
    parsed = ntInterpLine /@ lines;
(* count arguments over interpolator lookups only *)
    counts = Counts[Cases[parsed, p_List /; MemberQ[splines, p[[3]]] :> StringTrim[p[[4]]]]];
    shared = Select[Keys[counts], counts[#] >= 2 &];
    If[shared === {}, Return[text]];
    idxName = AssociationThread[Sort[shared] -> ("_ix" <> ToString[#] & /@ Range[0, Length[shared] - 1])];
    emitted = <||>;
    out = Reap[
        MapThread[
          Function[{ln, p},
            If[p === None || !MemberQ[splines, p[[3]]] || !KeyExistsQ[idxName, StringTrim[p[[4]]]],
              Sow[ln],
              Module[{ind = p[[1]], nm = p[[2]], cal = p[[3]], arg = StringTrim[p[[4]]], ix},
                ix = idxName[arg];
(* the first user of an argument owns the transform; every later one reuses the index *)
                If[!KeyExistsQ[emitted, arg],
                  emitted[arg] = True;
                  Sow[ind <> "const auto " <> ix <> " = " <> cal <> ".index(" <> arg <> ");"]];
                nrw++;
                Sow[ind <> "const auto " <> nm <> " = " <> cal <> ".at(" <> ix <> ");"]]]],
          {lines, parsed}]][[2]];
    ntLog["[interp] shared ", Length[shared], " of ", Length[counts],
      " distinct interpolator arguments (", nrw, " lookups rewritten)"];
    StringRiffle[If[out === {}, lines, First[out]], "\n"]];

(* FINITE MATSUBARA EXTENT. Each term of a finite-T kernel carries one dR/dt insertion. Where the
   insertion's ARGUMENT is coercive in the Matsubara symbol, the term dies super-polynomially in it,
   so `T Sum_n f(w_n)` is a finite sum DiFfRG can enumerate exactly instead of using a Gaussian rule.
   Whether the frequency is confined is read off the ARGUMENT, not the function name (one `RB` may
   regulate a 4D species and a 3D one, `RB[k^2, l1^2]`). The decay SHAPE cannot be read off the
   argument (Callan-Symanzik R = k^2 confines nothing), so it is declared per head via
   "DecayingRegulators". *)
ntCoerciveInQ[a_, ms_Symbol] :=
  With[{p = Expand[a]},
    PolynomialQ[p, ms] && Exponent[p, ms] === 2 && TrueQ[Positive[N[Coefficient[p, ms, 2]]]]];

(* True if `e` vanishes for large |ms|. Structural and conservative in the safe direction: a false
   "unbounded" costs an optimisation, a false "finite extent" truncates the sum.
     Plus  -- every summand must vanish (one surviving term is enough to ruin it),
     Times -- ONE vanishing factor suffices, because the declared decay is super-polynomial and the
              other factors (traces, propagators) grow at most polynomially,
     Power -- a positive integer power of something vanishing still vanishes. *)
ntFiniteExtentQ[e_, ms_Symbol, hs_List] :=
  Which[
    Head[e] === Plus, AllTrue[List @@ e, ntFiniteExtentQ[#, ms, hs] &],
    Head[e] === Times, AnyTrue[List @@ e, ntFiniteExtentQ[#, ms, hs] &],
    Head[e] === Power && IntegerQ[e[[2]]] && Positive[e[[2]]], ntFiniteExtentQ[e[[1]], ms, hs],
    (* Matched by NAME, not by symbol identity: the caller's RBdot is Global`RBdot while this file
       reads it from a package context, and MatsubaraVar already resolves its symbol the same way. *)
    Length[e] >= 1 && Head[Head[e]] === Symbol && MemberQ[hs, SymbolName[Head[e]]],
      ntCoerciveInQ[Last[e], ms],
    True, False];

(* ---- the kernel CLASS. Two shapes:
     regTemplate = False (default): a PLAIN class. Regulator calls (RB/RF/RBdot/...) in dressing
       rules are emitted UNQUALIFIED and supplied by the consumer (e.g. via "ExtraIncludes"); lookup
       runs class -> kernel namespace -> global, so global free functions are always found.
     regTemplate = True (fRG/DiFfRG shape): `template<typename REG> class ...` plus private wrappers
       forwarding to REG::, since DiFfRG instantiates the kernel as KERNEL<Regulator>.
   regAlias adds `using Regulator = REG;` and only makes sense with the template, so it implies it. *)
ntKernelClass[name_, members_List, decor_, regTemplate_, regAlias_, extraPriv_List] :=
  FunKit`MakeCppClass[
    Sequence @@
      If[TrueQ[regTemplate],
        {"TemplateTypes" -> {"REG"}},
        {}],
    "Name" -> name,
    "MembersPublic" ->
      If[TrueQ[regAlias],
        Prepend[members, "using Regulator = REG;"],
        members],
    "MembersPrivate" ->
      If[TrueQ[regTemplate],
        Join[{privDefs[decor]}, extraPriv],
        extraPriv]];

(* ---- emit helpers: support `using`s, namespace wrapping, runtime include -------------------
   A generated kernel pulls its math helpers (complex, powr, pow, sqrt, fma) from a configurable
   SUPPORT namespace and is wrapped in a configurable HOST namespace. The defaults (numtracer /
   numtracer/codegen/runtime.hpp) make the emitted code self-contained against NumTracer's own
   headers; a consumer that already provides equivalents (e.g. DiFfRG) points the codegen at
   them via "SupportNamespace" / "KernelNamespace" / "RuntimeInclude". *)
ntSupportUsings[sns_] :=
  DeleteDuplicates[{"using namespace " <> sns <> ";", "using namespace " <> sns <> "::compute;", "using namespace numtracer;"}];

(* the loop-independent `constant(p, k, dressings...)` function. DiFfRG flat-adds its return to the
   integral (constExpr second arg of MakeKernel). A 0 constant gets no namespace usings (keeps the
   emitted bytes unchanged); a nonzero expr may call compute helpers, so it gets the kernel body's. *)
ntConstFn[constExpr_, decor_, constParams_, sns_] := FunKit`MakeCppFunction[
    constExpr,
    "Name" -> "constant",
    "Prefix" -> decor,
    "Return" -> "auto",
    "CodeParser" -> "Cpp",
    "Parameters" -> constParams,
    "Body" ->
      If[MatchQ[constExpr, 0 | 0.],
        "",
        StringRiffle[ntSupportUsings[sns], "\n"]]];

ntWrapBody[kns_, classStr_, name_] := If[kns === None || kns === "",
    {classStr},
    {"namespace " <> kns <> "\n{", classStr, "}", "using " <> kns <> "::" <> name <> ";"}];

ntRuntimeIncludes[runInc_] := If[runInc === None || runInc === "",
    {},
    {runInc}];

ntApplyTraceComplexOverride[header_String, hdrInc_String, kns_String, sns_String, complexQ_] :=
  If[TrueQ[complexQ] && kns === "DiFfRG" && sns === "DiFfRG",
    StringReplace[header,
      "#include \"" <> hdrInc <> "\"" ->
        "#ifndef NT_TRACE_COMPLEX\n#define NT_TRACE_COMPLEX DiFfRG::complex<" <> $ntRealT <> ">\n#endif\n#include \"" <> hdrInc <> "\"",
      1],
    header];

(* the dressing-parameter type: Automatic -> `const auto&` (fully generic, self-contained);
   else the given concrete type string (e.g. a consumer's interpolator type). *)
ntDressType[dressTy_] := If[dressTy === Automatic,
    "auto",
    dressTy];

(* ---- mkGenerateKernel options (MakeNTKernel forwards its own; see Options[MakeNTKernel]). By default
   the emitted code is self-contained against NumTracer's headers; the namespace/include options
   point it at a consumer's support API instead. *)
Options[mkGenerateKernel] =
  {
    "Name" -> "nt_inv_kernel",
    "Namespace" -> Automatic,
    "Dressings" -> {},
    "ScalarParams" -> {},
    "ADParams" -> {},
    "ParameterOrder" -> Automatic,
    "IncludeDir" -> Automatic,
(* False: emit the generator sources only; the committed traces header is left untouched. *)
    "RunGenerator" -> True,
(* {sym -> expr}: kinematic angle symbols the dressing keeps SYMBOLIC, emitted once as
   `const double sym = ...;` so a shared sub-expression (a sqrt) is computed once. *)
    "AngleDefs" -> {},
(* True: lower all traces through ONE shared CSE program, trace_all(f, t[]), sharing subexpressions
   across traces. t[] is typed trace_all_t (complex iff some trace is), so nothing is truncated.
   Merging many small functions into one large block can spill and cost runtime. *)
    "CrossTraceCSE" -> False,
(* <|mom -> {e0,e1,e2,e3}|>: explicit momentum components, taken verbatim (polynomialised);
   Automatic picks the most compact frame parametrisation. *)
    "Components" -> Automatic,
(* <|sym -> expr|>: the C++ fill for DERIVED symbols (e.g. sin1 -> Sqrt[1-cos1^2]); plain free
   symbols are kernel arguments. *)
    "SymbolDefs" -> <||>,
(* Prefix on EVERY emitted function (kernel/constant, regulator wrappers, and via the generator's
   `-d` flag fill/trN/powr). Normalised by ntKokkosDecor. *)
    "Decorator" -> "static inline",
(* Does the kernel target DEVICE code? Drives gen.hpp's size-gated `__noinline__`, a device-only
   lever. Stated, not sniffed from the decorator: ntKokkosDecor rewrites `__device__` away, and the
   KOKKOS_ macros expand to plain `inline` on a host-only build. Automatic = host unless the raw CUDA
   spelling is passed; MakeNTKernelDiFfRG passes its own Device option. *)
    "DeviceTarget" -> Automatic,
(* Matsubara-frequency symbol of a finite-T flow. Asks the generator to PROVE the kernel even in it
   and, if so, emit DiFfRG's `matsubara_even` trait (one evaluation per mode instead of two).
   None = not a finite-T flow. *)
    "MatsubaraVar" -> None,
(* Regulator heads that decay super-polynomially in their argument; see ntFiniteExtentQ. Automatic =
   all six DiFfRG wrappers, not just the `dot` pair: a dressed insertion dt(Z R) = Zdot R + Z Rdot
   carries RB itself in half its terms. Only a factor in a PRODUCT counts (RB in a denominator makes
   the term grow). A regulator that decays only algebraically, or not at all, must NOT be listed. *)
    "DecayingRegulators" -> Automatic,
(* Automatic = derive the `matsubara_finite_extent` trait from the algebra. True/False force it;
   forcing True on a summand that does not die above the regulator's support silently truncates the
   Matsubara sum. *)
    "MatsubaraFiniteExtent" -> Automatic,
(* Support header #included first, providing `complex` and compute::{powr,pow,sqrt,fma}; None omits it. *)
    "RuntimeInclude" -> "numtracer/codegen/runtime.hpp",
(* Extra #includes ahead of everything, e.g. a header supplying the regulator functions. *)
    "ExtraIncludes" -> {},
(* Namespace wrapping the kernel class and the trace functions; None emits at the includer's scope. *)
    "KernelNamespace" -> "numtracer_kernels",
(* Where `complex`/`compute` are looked up via `using`. *)
    "SupportNamespace" -> "numtracer",
(* Dressing-parameter type: Automatic = `const auto&`, or a concrete type string. *)
    "DressingType" -> Automatic,
(* Opt-in: share one coordinate transform among lookups at the same argument (see
   ntShareInterpIndices). Needs index()/at() on the interpolator, a DiFfRG property. *)
    "ShareInterpolatorIndex" -> False,
(* Opt-in: turn launch-constant dressing lookups into host-evaluated kernel parameters. *)
    "HoistLoopConstLookups" -> False,
(* Kernel class shape; see ntKernelClass. RegulatorAlias implies RegulatorTemplate. *)
    "RegulatorTemplate" -> False,
    "RegulatorAlias" -> False,
(* When some coefficient carries an `i`, compile+run a probe on the generated traces to test whether
   Im(integrand) actually vanishes, and select a real kernel body if it does. *)
    "RealProbe" -> True,
(* Emit a `double` trace for groups whose dressing coefficient is real (only Re(trace) is consumed).
   Must be sequenced after the probe, which needs the unpruned traces; see ntkPruneSpec. *)
    "PruneRealTraces" -> False,
(* See $ntComplexRuntimeProjection. *)
    "ComplexRuntimeProjection" -> False,
(* See the note above MakeNTKernel::endproj. *)
    "ComplexEndProjection" -> False,
    "RealOutput" -> False,
(* The loop-INDEPENDENT piece flat-added to the integral, emitted as the body of
   `constant(p, k, dressings...)` like DiFfRG MakeKernel's constExpr. A plain Mathematica expression
   (ZA[p] -> ZA(p)), not an NTKernel. *)
    "Constant" -> 0.,
    (* True: emit the generator + probe sources and a numtrace.json switch set to 0, but compile and run
       nothing; the `numtrace` CMake target does that as a build step. NT_OFFLINE overrides. *)
    "Offline" -> False,
    (* coordinate argument NAMES of this flow's grid; see constArgQ *)
    "CoordinateArgs" -> Automatic
  };

(* The dressing-collection path needs no option: mkGenerateKernel detects the `ntDressedNum` tokens
   NumTrace emits under "DressingCollection" -> True and routes to the DPoly generator branch. *)

(* "ComplexRuntimeProjection": what to do when a coefficient carries a Complex below a head the
   symbolic real/imag split cannot traverse — in practice a finite-density denominator l0 + I muq,
   where the split is silently wrong (see ntIiSafeQ). OFF: refuse (MakeNTKernel::cplxnest). ON: keep
   the coefficient intact and project with ntRe/ntIm at C++ runtime (complex arithmetic in the
   kernel, needs the type-generic powr). Every summand, including the Pure body, then takes the exact
   per-summand path. A package global because ntProjectIntegrand reads it several call layers down. *)
$ntComplexRuntimeProjection = False;

(* "ComplexEndProjection": skip the symbolic Pure/RePart projections and the probe; emit one body
   ntRe[integrand], keeping finite-density denominators complex until runtime. This is the pointwise
   real part, not a proof that Im cancels. Requires RealOutput -> True. No fixture exercises it.
   "RealOutput": the consumer takes a REAL value, so the complex body (the expensive one to lower) is
   not emitted. A genuinely complex verdict then falls through to RePart, i.e. Re[flow]: a TRUNCATION
   of the flow equation, flagged by a #warning in the header. Not the default, and not set by
   MakeNTKernelDiFfRG: that truncation must be the consumer's explicit choice. *)

MakeNTKernel::endproj = "ComplexEndProjection requires RealOutput -> True; otherwise the consumer would receive a complex kernel return value instead of the requested endpoint real projection.";

MakeNTKernel::nonnumeric = "the integrand is not numeric: it contains `1` Indeterminate/Infinity value(s), which would be emitted as bare Mathematica symbols in the generated C++ (e.g. `return Indeterminate;`). This almost always means a SINGULAR Gram at the chosen kinematics: the basis's structures are linearly dependent there, so the inverse metric, and every dual projector built from it, carries 0/0. Check Det[TBGetMetric[basis]] under the frame's kinematics (e.g. the symmetric point), and project with a restricted sub-basis whose Gram is non-degenerate. Offending value(s): `2`";

mkGenerateKernel::genfail = "Generator compile/run failed: `1`";
mkGenerateKernel::pruneoff = "PruneRealTraces ignored: the RealProbe cannot run in this mode (Offline or RunGenerator->False), so pruned traces could not be probe-validated. Emitting all-complex traces; set RealProbe->False to assert the flow is safe to prune without a probe.";

mkGenerateKernel::emptynets = "Flow `1` produced no generator nets (nets=`2`, groups=`3`) — nothing to emit. Aborting instead of writing a placeholder kernel. (Either NumTrace returned no usable diagrams, or every diagram was dropped during the net build.)";

mkGenerateKernel::scalarleak = "Diagram `1`: a non-numeric factor `2` reached the generator scalar coefficient (a Lorentz tensor that was not resolved by the net builder, e.g. an un-anchored metric contraction). It would be emitted as undeclared C++. Aborting; fix the net build (compileLorentz) so the contraction folds numerically.";

(* ---- runtime-parameter typing ----------------------------------------------------------------
   A parameter name is spelled as a Symbol by a hand-written flow file and as a String by
   DiFfRG_compat (which reads it out of a DiFfRG parameter Association), and option names may be
   either too. Every list that gets MemberQ'd against another must therefore be normalised the same
   way first; ntParamName is that one normalisation. *)
ntParamName[nm_String] := nm;
ntParamName[nm_Symbol] := SymbolName[nm];
ntParamName[nm_] := ToString[nm];

(* The C++ type of ONE runtime parameter (scalar or dressing), as a string for mkParam. Every input
   is an explicit argument so the decision is unit-testable (tests/test_ad_param_typing.wls).
   The AD test comes FIRST: a caller may pass ADParams without listing them in ScalarParams, and an
   AD scalar's own entry says "Type" -> "double" (the AD flag is a separate key), so neither the
   scalarNames gate nor the entry's type may decide it. *)
ntRuntimeParamType[entry_, adNames_List, scalarNames_List, dressTy_] :=
  With[{nm = ntParamName[If[AssociationQ[entry], entry["Name"], entry]]},
    Which[
(* `auto` is not cosmetic: it makes the emitted function an abbreviated template, which is the only
   reason kernel()/constant() can bind autodiff::real from DiFfRG's integrator_AD twin. *)
      MemberQ[adNames, nm], "auto",
      MemberQ[scalarNames, nm], $ntRealT,
      AssociationQ[entry], entry["Type"],
      True, dressTy[nm]]];

(* POST-CONDITION on a finished mkParam list: every ADParams name in the signature is typed auto.
   kernel(), constant() and ntHoisted share runtimeParams, so one check covers all three. Nothing
   else sees a mistyped AD scalar: generation succeeds and only the consumer's AD twin fails to
   compile. Returns the list unchanged so it can be used inline. *)
ntAssertADTyped[params_List, adNames_List] :=
  With[{missed = Select[params, MemberQ[adNames, #["Name"]] && #["Type"] =!= "auto" &]},
    If[missed =!= {},
      Message[MakeNTKernel::adtype, Length[missed],
        StringRiffle[(#["Name"] <> " got " <> ToString[#["Type"]])& /@ missed, ", "]];
      Abort[]];
    params];

(* ---- STAGE: which trace groups may drop their imaginary half ("PruneRealTraces") --------------
   A group whose dressing coefficient is REAL has only Re(trace) consumed (the consumer takes Re of
   the whole kernel; a real coefficient cannot move Im(trace) into the real part). Flag it and the
   generator emits a `double` trace and never computes the dead imaginary half. The default is the
   empty list => all-complex traces, which is also what the probe needs in order to verify the
   cancellation it is checking for.

   diagData is an ARGUMENT, never read from an enclosing scope: an unassigned diagData makes
   `FreeQ[..., Complex]` vacuously True, silently marking every group real. The guards below make
   that failure loud.

   ORDERING: the probe checks Im(integrand)~0 on the GENERATED traces, so it must see the UNPRUNED
   set (pruned traces can misclassify the flow — an O(1) kernel error). When the probe will run, pass
   1 generates unpruned and the prune is applied in a second pass after the verdict. Without a probe:
   a syntactically real flow prunes directly, RealProbe->False prunes on the caller's assertion, and
   offline/no-generator drops the request with a warning. *)
mkGenerateKernel::prunedata = "PruneRealTraces: the diagram-coefficient table has `1` entries but the trace grouping references diagram index `2`. The table and the grouping have gone out of step, so the per-group real/complex verdict below would be read off the wrong diagram — or off nothing at all. This is the shape of the bug this guard exists for (an unassigned diagData read as vacuously real).";

mkGenerateKernel::pruneall = "PruneRealTraces: the flow is COMPLEX, yet every one of its `1` trace groups was judged real and would have its imaginary half dropped. That is exactly what the historic bug produced from an unassigned coefficient table, and it is indistinguishable by inspection from a legitimately all-real grouping. If it is legitimate, the flow should not have tripped the complex test at all. Refusing to emit; re-run without \"PruneRealTraces\" -> True to get the full complex kernel.";

ntkPruneSpec[diagData_, groups_, complexQ_, offline_, pruneRequested_, realProbe_, runGenerator_] :=
  Module[{pruneG, realOnlyG, probeWillRun},
    pruneG =
      If[pruneRequested,
        Module[{maxDiag = Max[Append[#[[1]]& /@ groups, -1]]},
(* Index guard: an out-of-range Part would stay unevaluated, and FreeQ would call it real. *)
          If[maxDiag >= Length[diagData],
            Message[mkGenerateKernel::prunedata, Length[diagData], maxDiag]; Abort[]];
          (FreeQ[diagData[[#[[1]] + 1]], Complex])& /@ groups],
        {}];
(* "Every group is real" on a flow the complex test flagged is the signature of a lost coefficient
   table, and no oracle would catch it (only the DROPPED imaginary parts are wrong). Refuse it. *)
    If[pruneRequested && TrueQ[complexQ] && pruneG =!= {} && AllTrue[pruneG, TrueQ],
      Message[mkGenerateKernel::pruneall, Length[pruneG]]; Abort[]];
    probeWillRun = TrueQ[complexQ] && realProbe && runGenerator && !TrueQ[offline];
    realOnlyG =
      Which[
        !pruneRequested,   {},
        !TrueQ[complexQ],  pruneG,
        probeWillRun,      {},
        realProbe,         Message[mkGenerateKernel::pruneoff]; {},
        True,              pruneG];
    ntStageResult["ntkPruneSpec", {"pruneG", "realOnlyG", "probeWillRun"},
      <|"pruneG" -> pruneG, "realOnlyG" -> realOnlyG, "probeWillRun" -> probeWillRun|>]];

(* ---- mkGenerateKernel: STAGE MAP --------------------------------------------------------------
   The whole generation pipeline for ONE flow, from an analysed NTKernel to three files on disk: the
   build-time generator sources (genFile + its units + decl), the committed straight-line traces
   header (headerFile, printed by RUNNING that generator), and the kernel header the consumer
   includes (kernelFile).

   The stages, in order:

      1. OPTIONS. Normalise every option and every parameter name exactly once (Symbol and String
         spellings collapse here, and nowhere else).
      2. FRAME SPEC. Probe which frame parametrisation the flow qualifies for (unit-loop / mixed
         unit-loop / general polynomial) and build the component table over it.
      3. RESET. Clear the per-generation memo caches and stamp $ctCtx.
      4. NET BUILD. Per diagram: split into colour groups, compile the Dirac chain and the Lorentz
         remainder, and accumulate the nets, the colour tokens and the per-diagram coefficient data.
      5. GROUPING. Partition the nets into trace groups — additive groups (summed into one trace)
         and factor groups (kept separate because a dressed factor entry carries its own token).
      6. INTEGRAND. Assemble the symbolic integrand over the groups, hoist launch-constant dressing
         lookups, and decide the prune spec (ntkPruneSpec, above).
      7. KERNEL LOWERING. Lower the integrand through FunKit into the kernel class/header text.
      8. EMIT + BUILD (genPass). Emit the generator (emitNumericGenerator), write it, then either
         compile and RUN it here (online) or leave it to the `numtrace` CMake target (offline).
      9. PROBE + RE-PRUNE. On a complex flow, probe whether the imaginary part cancels, write the
         verdict header, and — if PruneRealTraces was requested — re-run stage 8 with the prune
         applied. genPass is a closure precisely so it can be re-run.
     10. WRITE + MANIFEST. Write the kernel header (write-if-changed) and record the numtrace.json
         manifest for the build.

   Several Module locals are DECLARED here but ASSIGNED inside inner scopes (e.g. diagData,
   lorFacOf, pGroupOf, nAdd). Re-declaring one of them in an inner Module shadows it and silently
   leaves the outer one unassigned. Extracted stages return an Association checked by ntStageResult.
   The Block scopes the package globals assigned below to this one generation. *)
mkGenerateKernel[NTKernel[k_], genFile_, kernelFile_, headerFile_, OptionsPattern[]] :=
  Block[{$RecursionLimit = $RecursionLimit, $ctCtx = $ctCtx, $ntDressResolve = $ntDressResolve,
         $ntCanonIdsSrc = $ntCanonIdsSrc, $ntCanonRules = $ntCanonRules,
         $ntComplexRuntimeProjection = $ntComplexRuntimeProjection},
  Module[{name, ns, dress, scalarParams, adParams, parameterOrder, adNames, scalarParamNames, args, sigArgs, frame, env, nonzeroCompMask, ncomp, fillArgs, fillArgSig, constArgQ, invNets, invRest, g, colourNets, preamble, integrand, kernelParams, runtimeParams, constParams, mkParam, kernelFn, constFn, classStr, header, hdrInc, incDir, genPre, genUnits, genDecl, genMain, declFile, pchFile, unitFiles, genSrc, bin, complexQ, angleDefs, angleDecls, crossCSE, traceRef, nGrp, decor, tarrDecl, kns, sns, runInc, extraInc, interpTy, nsHome, regTemplate, regAlias, offline, realOut, endProject, verdictMacro, probeFile = None, mainOptForManifest, symDefs = <||>, realOnlyG = {}, pruneG = {}, probeVerdict = None, genPass,
    hoistCalls = {}, hoistSyms = {}, hoistFnStr = "", mVarIdx = -1, mSym = None, mEvenBody = False, mFiniteExtentBody = False, mSplit = False, feExpr = 0, tailExpr = 0, splitFns = {}, mkKernelFnNamed, timedBodyNamed, bodyFor, dressedIdx = {}, diagTokExpr = {}, factorNets = {}, lorFacOf = {}, pGroupOf = <||>, nAdd = 0, factorCompOf = <||>,
(* diagData is assigned by the net-build Module below but read after it closes (ntkPruneSpec), so
   it must be declared HERE; declared in the inner Module it would be unassigned at that read. *)
    diagData = {}},
    Needs["FunKit`"];
(* The integrand Sum and COEN's lowering recurse ~linearly in the number of trace groups (>1000 on
   large flows). Hitting $RecursionLimit does NOT abort: it returns a held expression and the kernel
   is silently skipped. Raise the limit; a real runaway still hits the ceiling. *)
    $RecursionLimit = Max[$RecursionLimit, 1048576];
    name = OptionValue["Name"];
    dress = OptionValue["Dressings"];
(* loop-independent scalar doubles threaded into the signature *)
    scalarParams = OptionValue["ScalarParams"];
    parameterOrder = OptionValue["ParameterOrder"];
(* AD-flagged scalars (d1V, d2V for FE-potential flows) must be `const auto&` so the kernel also
   accepts autodiff::real from the integrator_AD twin; everything else stays `const double&`. *)
    adParams = OptionValue["ADParams"];
(* Both name lists are normalised ONCE, here, as Module locals; nothing below rebinds them. *)
    adNames = ntParamName /@ adParams;
    scalarParamNames = ntParamName /@ scalarParams;
    angleDefs = OptionValue["AngleDefs"];
    crossCSE = OptionValue["CrossTraceCSE"];
    (* normalise raw CUDA qualifiers to the Kokkos macros — see ntKokkosDecor. *)
    decor = ntKokkosDecor[OptionValue["Decorator"]];
    (* Offline: emit sources + the per-flow numtrace.json switch and let the `numtrace` build target do
       the compiling/running. NT_OFFLINE overrides the option ("0"/"false" force online). *)
    offline =
      With[{e = Environment["NT_OFFLINE"]},
        If[StringQ[e] && StringTrim[e] =!= "",
          !MemberQ[{"0", "false", "no", "off"}, ToLowerCase[StringTrim[e]]],
          TrueQ[OptionValue["Offline"]]]];
    (* main-TU -O level, recorded in the manifest so the offline build matches the online generator.
       Must mirror `mainOpt` in the RunGenerator block. *)
    mainOptForManifest = With[{e = Environment["NT_GEN_MAIN_OPT"]}, If[StringQ[e] && e =!= "", e, "-O1"]];
    kns = OptionValue["KernelNamespace"];
    sns = OptionValue["SupportNamespace"];
    runInc = OptionValue["RuntimeInclude"];
    extraInc = OptionValue["ExtraIncludes"];
(* the fRG/DiFfRG kernel shape (template<typename REG> + the private REG:: wrappers) is opt-in; the
   Regulator alias is meaningless without it, so it implies it. *)
    regAlias = TrueQ[OptionValue["RegulatorAlias"]];
    regTemplate = TrueQ[OptionValue["RegulatorTemplate"]] || regAlias;
    interpTy = ntDressType[OptionValue["DressingType"]];
(* Which of `args` the loop-independent constant() receives. Automatic = p and k (a 1-D grid); a
   caller with another grid passes its coordinate names via "CoordinateArgs". k is always included. *)
    constArgQ =
      With[{ca = OptionValue["CoordinateArgs"]},
        If[ca === Automatic,
          Function[a, a === Global`p || a === Global`k],
          Function[a, a === Global`k || MemberQ[ca, If[StringQ[a], a, ToString[a]]]]]];
    ns = OptionValue["Namespace"] /. Automatic -> ToLowerCase[name];
(* where the generated trace fns / nenv / fill live *)
    nsHome = kns <> "::" <> ns;
(* the #if macro selecting one of the 3 complex-kernel bodies *)
    verdictMacro = ntVerdictMacro[ns];
    args = k["Args"];
    frame = k["Frame"];
    env = k["Env"];
    nonzeroCompMask = Association @ KeyValueMap[#1 -> frameMask[resolveComponents[#1, frame]]&, env];
(* scalars the fill needs *)
    fillArgs = Select[args, # =!= Global`k&];
(* FRAME SPEC + component table, timed as one block: each frame-spec test runs a Simplify sweep over
   the whole frame, and a general frame pays all of them before numericComponents starts. *)
    With[{ntT = First @ AbsoluteTiming[
    ncomp =
      Module[{uc = OptionValue["Components"], ud = OptionValue["SymbolDefs"], pf, ad, ug},
(* WHICH PARAMETRISATION THE FLOW GETS, most compact first. All three produce the same {components,
   symbol defs, unit groups} triple; they differ only in how the LOOP momentum is written, and that
   decides how big every polynomial downstream is.

     polyFrameSpec[uc]        explicit user "Components". Taken verbatim (polynomialised for any
                              Sqrt/trig). No unit groups — the caller owns the parametrisation.
     unitLoopFrameSpec        the compact vacuum symmetric-point case: BOTH the externals and the
                              loop ride opaque unit-direction symbols with a ΣU²=1 constraint, so a
                              denominator collapses to `l0² + l1²` instead of an angle polynomial.
     unitLoopMixedFrameSpec   the loop rides unit directions but the externals are NUMERIC vectors
                              (a general external configuration, or finite T where the heat bath
                              picks out slot 0). Far more compact than the fallback.
     polyFrameSpec[frame]     the general fallback: every component a polynomial in the frame's
                              scalars, no unit constraint. Always correct, never compact.

   The ORDER is the specification: each test is strictly narrower than the next, so the first that
   qualifies is the most compact one available. NT_NO_UNIT_GROUPS disables both unit-loop branches
   (tests/gen/gen_lambda3d_small_numeric.wls builds its control that way). *)
        {pf, ad, ug} =
          Which[
            uc =!= Automatic,
              Append[polyFrameSpec[uc], {}],
            unitLoopOkQ[frame, Global`p, Global`l1],
              unitLoopFrameSpec[frame, Global`p, Global`l1],
            unitLoopMixedOkQ[frame, Global`l1],
              unitLoopMixedFrameSpec[frame, Global`l1],
            True,
              Append[polyFrameSpec[frame], {}]
          ];
        symDefs = Join[ad, ud];
        numericComponents[env, pf, symDefs, ug]];]},
      ntLog["[prof] numericComponents + frame spec: ", ntT, " s"]];
(* The Matsubara frequency, in two forms. `mSym` is the SYMBOL, used by the Mathematica-side
   evenness test and finite-extent partition. `mVarIdx` is its MPoly variable index, used only by the
   generator to prove evenness of the TRACES. -1 is not an error: a purely SCALAR integrand has no
   momentum components, yet still depends on the frequency through its coefficient (denominators,
   regulator arguments). The lookup therefore spans the frame's symbols AND the fill arguments. *)
    mSym =
      With[{mv = OptionValue["MatsubaraVar"]},
        If[mv === None || mv === Automatic,
          None,
          With[{cands = Select[Join[ncomp["usyms"], fillArgs], SymbolName[#] === ToString[mv]&]},
            If[cands === {}, $Failed, First[cands]]]]];
(* A name that matches nothing is reported loudly, not turned into a quiet None: otherwise the
   caller gets a valid kernel without the requested optimisation and no way to tell. *)
    If[mSym === $Failed,
      Print["[NumTracer] WARNING: \"MatsubaraVar\" -> ", OptionValue["MatsubaraVar"],
        " names neither a frame symbol nor a kernel fill argument, so no evenness check was run ",
        "and no Matsubara trait will be emitted. The frame's symbols are: ", ncomp["usyms"],
        ", the fill arguments are: ", fillArgs,
        ". (The integration-variable name DiFfRG uses, e.g. \"f\", is often NOT the frame symbol, ",
        "e.g. \"f0\" — this option wants the frame symbol.)"];
      mSym = None];
    mVarIdx =
      If[mSym === None,
        -1,
        With[{pos = Position[ncomp["usyms"], mSym, {1}]}, If[pos === {}, -1, pos[[1, 1]] - 1]]];
    If[mSym =!= None,
      ntLog["[matsubara] frequency symbol ", mSym, " = ",
        If[mVarIdx >= 0,
          "MPoly var " <> ToString[mVarIdx] <> " — trace evenness will be proven at generation time",
          "not a momentum-component variable (scalar integrand) — the traces carry no MPoly " <>
            "variable, so the traits are decided entirely here"]]];
(* [[maybe_unused]]: a frame may not reference every fill() argument (e.g. an angle or dressing atom
   that only some diagrams use), so mark each parameter to keep the emitted kernel -Wunused-clean. *)
    fillArgSig = StringRiffle[("[[maybe_unused]] " <> $ntRealT <> " " <> SymbolName[#])& /@ fillArgs, ", "];
(* A syntactic `i` anywhere (e.g. projector i x imaginary non-abelian colour f^abc T^b T^c =
   (iN/2) T^a) makes the flow complexQ; the colour constant stays COMPLEX and the probe decides
   whether the assembled integrand is actually real. *)
    complexQ = !FreeQ[k["Diagrams"], Complex];
    preamble = {};
(* net -> factor-id *)
    factorCompOf = <||>;
(* Every net-builder memo depends on generation-fixed state ($ntDressResolve, env, nonzeroCompMask,
   frame), so all of them are cleared here, at the start of each generation. *)
    $ctCache = <||>;
    $dsCache = <||>;
    $odCache = <||>;
    $dslCache = <||>;
(* the generation-fixed half of every net-builder memo key, hashed ONCE here instead of on each
   (recursive) call. Covers everything those builders read besides the expression and its ids. *)
    $ctCtx = Hash[{env, nonzeroCompMask, frame}];
(* clear the per-component diagonal-dressing and scalar-dressing (ntDressedNum) registries *)
    resetDiagDr[];
    resetDr[];
(* frame resolver for dressed-numerator option coefficients (compileDirac → dressedSlotStr): the same
   ntSP/ntSPS/ntVec[q,i] → component substitution used for diag["Coeff"] below. *)
    $ntDressResolve =
      Function[s,
        s /. {ntSP[x_, y_] :> resolveComponents[x, frame] . resolveComponents[y, frame], ntSPS[x_, y_] :> Rest[resolveComponents[x, frame]] . Rest[resolveComponents[y, frame]], ntVec[q_, ii_Integer] :> resolveComponents[q, frame][[ii + 1]]}
      ];
(* NET BUILD: walk the diagrams; each is a Lorentz/Dirac trace x a colour factor x a dressing
   coefficient, recorded per net as {colNet, cores, coeff, rest, lorFac, factorId}.
   The accumulators are BAGS, not lists: Append copies, so N appends would be O(N^2) in the net
   count (often ~30x the diagram count). Only the running net index `nNetAcc` is read during the
   loop; the lists are materialised once afterwards.
   NB no `diagData` in this local list: it is declared in the OUTER Module because ntkPruneSpec
   reads it after this Module closes; re-declaring it here would shadow it. *)
    Module[{bInvNets = Internal`Bag[], bInvRest = Internal`Bag[], bColourNets = Internal`Bag[], bDiagData = Internal`Bag[], bLorFacOf = Internal`Bag[], bFactorNets = Internal`Bag[], nNetAcc = 0},
(* The net build itself is bound here, not passed as an ntLog argument: it is the work, not a
   diagnostic. See ntExportCpp for why load-bearing work must stay outside ntLog. *)
      $ntProf = <||>;
      With[{ntT =
        First @
          AbsoluteTiming[
            Block[{$ntProfOn = TrueQ[$NumTracerVerbose]}, MapIndexed[
              Function[{diag, di},
                Module[{coeff, colBr, constAcc = {}, d = di[[1]] - 1, pureLorAcc = {}, diracComps = {}},
(* cache this diagram's canonicalisation rules for ntCanonIds (see there) — one Dispatch per
   diagram instead of one Normal[KeyDrop[...]] per compileLorentz/diracSlotStr call. *)
                  $ntCanonIdsSrc = diag["Ids"];
                  $ntCanonRules = Dispatch[Normal[KeyDrop[diag["Ids"], Keys[env]]]];
                  coeff = diag["Coeff"] /. {ntSP[x_, y_] :> resolveComponents[x, frame] . resolveComponents[y, frame], ntSPS[x_, y_] :> Rest[resolveComponents[x, frame]] . Rest[resolveComponents[y, frame]], ntVec[q_, i_Integer] :> resolveComponents[q, frame][[i + 1]]};
                  MapIndexed[
                    Function[{comp, ci},
                        If[comp["Constant"],
(* Constant SU(N) component (colour and/or flavour; each head carries its own rank N). The fold
   is COMPLEX (sun_value_cx), so an imaginary non-abelian colour survives into the trace. A
   diagram may carry SEVERAL constant components (e.g. a colour AND a flavour trace): ACCUMULATE
   all factors and compile their product once after the loop (colBr), since a component may be a
   PLUS (e.g. the Fierz flavour structure δδ - 4·T·T) with no single-net representation. *)
                          constAcc = Join[constAcc, comp["Factors"]],
(* Non-constant component. The DISCONNECTED components of ONE diagram MULTIPLY — they are NOT
   separate summed diagrams. Collect them for the post-loop assembly. Route by structure:
     - any colour (entangled in a Plus, or a top-level T^a × …) or a gamma chain: collect the
       component's splitColourGroups entries (dirac_value net per colour branch / chunked
       Lorentz net) as ONE Dirac/colour component in diracComps;
     - pure Lorentz (no colour, no gamma): accumulate factors — ALL pure-Lorentz components
       fold into ONE product net (disjoint ids make the C++ contract_factors multiply them). *)
                          If[colourEntangledQ[comp["Factors"]] || !FreeQ[comp["Factors"], _ntGamma | _ntGamma5 | _ntC | _ntDeltaDirac | _ntDressedNum | _ntDiracSlot],
                            AppendTo[diracComps, ntProfTimed["splitColourGroups", splitColourGroups[comp["Factors"], diag["Ids"], env, nonzeroCompMask]]],
                            pureLorAcc = Join[pureLorAcc, comp["Factors"]]]]],
                    diag["Components"]];
(* The diagram's CONSTANT colour/flavour part, as a list of {netString, scalar} branches — one
   branch unless a constant component was a sum. Colour folds to a scalar (sun_value_cx -> Cx)
   and the generator already sums colour by emitting several nets into one group, so a summed
   colour component costs one extra net record per branch and needs no C++ support. *)
                  colBr =
                    If[constAcc === {},
                      {{"SUNNet{}", 1}},
                      ntProfTimed["compileColourSum", compileColourSum[Times @@ constAcc, diag["Ids"]]]];
(* ---- assemble the diagram's nets from its non-constant components -----------------------------
   A diagram with K disconnected non-constant components is a PRODUCT of K independent closed
   scalars (each a Dirac/colour trace or a pure-Lorentz scalar): `coeff * Times @@ toks`.
   ALL pure-Lorentz components fold into ONE product factor; each
   Dirac/colour component is its own factor. Two regimes:
     - <= 1 non-constant factor: the additive path. The single Dirac component's entries
       (colour folded, GlobalCollect-fusible) OR the single combined pure-Lorentz product net is
       appended with diagData = coeff*scal.
     - >= 2 factors: FACTORED product. One Dirac component is the additive BASE (traceRef[gi],
       diagData = coeff*scal); every other component is emitted as its own fused trace GROUP
       (its per-entry scalar folded into the net, its colour folded by the group sum) and tagged
       with a factor id. The base carries the list of factor ids; the assembly multiplies in
       Π traceRef[factor groups] — P computed ONCE per component, no trace-polynomial blow-up.
   scalarleakCheck: each entry's restNet scalars (e[[4]] {restStr,scal}) become the generator's
   `dsc[]` numeric constants — a symbolic Lorentz tensor the net builder failed to fold would be
   CForm'd into undeclared C++. Catch it loudly, with the offender. *)
                  Module[{nDir = Length[diracComps], hasLor = pureLorAcc =!= {}, factorComps, factorIds = {}, scalarleakCheck, appendRec},
                    scalarleakCheck =
                      Function[es,
                        Do[
                          Module[{badS = FirstCase[ee[[4]], {_, s_} /; !NumericQ[s] :> s, Missing[]]},
                            If[!MissingQ[badS],
                              Message[mkGenerateKernel::scalarleak, d, badS];
                              Abort[]]],
                          {ee, es}]];
                    (* rec = {colNet, cores, dData, restList, lorFac(None|{ids..}), factorId(None|id)}; returns net idx *)
                    appendRec =
                      Function[rec,
                        Module[{ni = nNetAcc},
                          Internal`StuffBag[bInvNets, rec[[2]]];
                          Internal`StuffBag[bInvRest, rec[[4]]];
                          Internal`StuffBag[bColourNets, rec[[1]]];
                          Internal`StuffBag[bDiagData, rec[[3]]];
                          Internal`StuffBag[bLorFacOf, rec[[5]]];
                          nNetAcc = ni + 1;
                          If[rec[[6]] =!= None,
                            Internal`StuffBag[bFactorNets, ni];
                            factorCompOf[ni] = rec[[6]]];
                          ni]];
                    If[nDir + Boole[hasLor] <= 1,
                      (* ---- single non-constant factor (or none): additive path ---- *)
                      Module[{
                        baseEntries =
                          Which[
                            nDir == 1,
                              diracComps[[1]],
                            hasLor,
                              ({"SUNNet{}", {#[[1]]}, #[[2]], {{"", 1}}}& /@ ntProfTimed["chunkLorentz", chunkLorentz[Times @@ pureLorAcc, diag["Ids"], env, nonzeroCompMask]]),
                            True,
                              {}]},
                        scalarleakCheck[baseEntries];
                        Do[
                          appendRec[{mergeColNet[cb[[1]], e[[1]]], e[[2]], coeff cb[[2]] e[[3]], e[[4]], None, None}],
                          {cb, colBr},
                          {
                            e,
                            If[baseEntries === {},
                              {{"SUNNet{}", {"konst(1.0)"}, 1, {{"", 1}}}},
                              baseEntries]}]],
(* ---- >= 2 factors: factored product of disconnected components ----
   EVERY non-constant component becomes its OWN fused trace group: its colour-branch entries
   are summed WITHIN the group (GlobalCollect-style, colour folded by the group sum, the entry
   scalar folded into the net) so the group trace IS the component scalar. The diagram then
   contributes ONE additive ANCHOR term = coeff * colv(col) * Π(component group traces).
   Each component stays exactly ONE trace; splitting its entries into singletons would defeat
   colour-channel fusion and explode the trace count. *)
                      (
                        factorComps =
                          If[hasLor,
                            Append[diracComps, {{"SUNNet{}", {#[[1]]}, #[[2]], {{"", 1}}}}&[compileLorentz[Times @@ pureLorAcc, diag["Ids"], env, nonzeroCompMask]]],
                            diracComps];
                        Do[
                          Module[{compEntries = factorComps[[ci]], fid = nNetAcc},
                            scalarleakCheck[compEntries];
                            (* here the entry scalar e[[3]] is folded into the rest scalars below, so it
                               must be numeric too *)
                            Do[If[! NumericQ[ee[[3]]], Message[mkGenerateKernel::scalarleak, d, ee[[3]]]; Abort[]],
                              {ee, compEntries}];
                            AppendTo[factorIds, fid];
                            Do[appendRec[{e[[1]], e[[2]], 1, ({#[[1]], #[[2]] e[[3]]}&) /@ e[[4]], None, fid}], {e, compEntries}]
                          ],
                          {ci, 1, Length[factorComps]}];
(* the anchor: a trivial unit net carrying the diagram coeff, the constant colour branch
   (folded once, via the group sum colv(net)*1), and the list of component factor ids.
   ONE anchor per colour branch — each is an independent additive term
   coeff·scal_b·colv(net_b)·Π(factor traces), and they all reference the SAME factorIds,
   so the expensive component traces are computed once and shared. *)
                        Do[appendRec[{cb[[1]], {"konst(1.0)"}, coeff cb[[2]], {{"", 1}}, factorIds, None}], {cb, colBr}]
                      )]]]],
              k["Diagrams"]]]]},
        ntLog[
          "[prof] per-diagram net-build (", Length[k["Diagrams"]], " diagrams): ", ntT, " s"]];
      ntProfReport["[prof]   net-build part "];
      (* materialise the bags once — everything downstream indexes these as plain lists *)
      invNets = Internal`BagPart[bInvNets, All];
      invRest = Internal`BagPart[bInvRest, All];
      colourNets = Internal`BagPart[bColourNets, All];
      diagData = Internal`BagPart[bDiagData, All];
      lorFacOf = Internal`BagPart[bLorFacOf, All];
      factorNets = Internal`BagPart[bFactorNets, All];
(* memo sizes should track the number of DISTINCT Dirac/Lorentz structures, not the call count;
   if they grow with the call count, a memo has stopped hitting. *)
      ntLog["[prof]   memo sizes: compileLorentz ", Length[$ctCache], " | orderDiracLoops ", Length[$odCache], " | dressedSlotStr ", Length[$dsCache], " | diracSlotStr ", Length[$dslCache], " (nets ", nNetAcc, ")"];
(* ---- per-component diagonal dressings (ntSUNDiag{Fund,Adj}) ----------------------------------
   A diagram whose colour net carries a diag factor folds (via the validated C++ engine,
   sun_value_dressed, run through the build-time seam) to a SUNPoly Σ_t coeff_t Π D^{dr}, where
   each dr names a distinct SCALAR runtime dressing (the surviving/named components; dropped ones
   vanish from the sum). Since the dressings are runtime, the colour can't be a compile-time
   constant baked into the trace: instead we (a) replace the diagram's colour net by the IDENTITY
   (so the generator folds colv=1 and the trace stays colour-free), and (b) build a runtime token
   `Σ_t coeff_t Π name(scale)` — ordinary scalar-dressing tokens — multiplied into the integrand.
   The Dirac/Lorentz trace is still computed ONCE, not a diagram per component. *)
      diagTokExpr = Table[1, {Length[colourNets]}];
      dressedIdx = Select[Range[Length[colourNets]], StringContainsQ[colourNets[[#]], ".diag"]&];
      If[dressedIdx =!= {},
        Module[{polys, resolveScale, incDir = OptionValue["IncludeDir"] /. Automatic :> resolveIncludeDir[]},
          resolveScale[s_] := s /. {ntSP[x_, y_] :> resolveComponents[x, frame] . resolveComponents[y, frame], ntSPS[x_, y_] :> Rest[resolveComponents[x, frame]] . Rest[resolveComponents[y, frame]], ntVec[q_, ii_Integer] :> resolveComponents[q, frame][[ii + 1]]};
          polys = diagColPolys[colourNets[[dressedIdx]], incDir];
          MapThread[
            Function[{d, p},
              diagTokExpr[[d]] =
                Total[
                  Function[term,
                      (term[[1]] + I term[[2]]) *
                        (
                          Times @@
                            (
                              Function[dr,
                                  resolveScale[$diagDrTable[dr]["Expr"]]
                                ] /@ term[[3]]))
                    ] /@ p];
              colourNets[[d]] = "SUNNet{}"],
            {dressedIdx, polys}]];
        ntLog["[prof] diagonal-dressed diagrams: ", Length[dressedIdx], " (per-component colour-sum folded via sun_value_dressed seam)"]
      ];
(* trace reference: with CrossTraceCSE the kernel fills a `tarr[]` once via trace_all() and reads
   tarr[i]; otherwise it calls the independent tr_i(fenv). *)
      traceRef =
        If[crossCSE,
          "tarr[" <> ToString[#] <> "]",
          nsHome <> "::tr" <> ToString[#] <> "(fenv)"
        ]&;
(* Group diagrams (0-based) into traces. Colour is folded numerically into the generator polynomial,
   so diagrams FUSE by identical dressing coeff and the kernel evaluates ~one polynomial per Feynman
   graph. Exceptions stay singletons: diag-dressed diagrams (diagTokExpr =!= 1) carry a per-diagram
   RUNTIME colour-sum token (they cannot fuse by dressing coefficient alone, the token differs).
   lorFac diagrams (lorFacOf =!= None — a Dirac trace times a disconnected pure-Lorentz scalar)
   likewise stay singletons: each carries a per-diagram multiplicative trace, so it must not fuse
   with another diagram's entries.
   The FACTOR nets (P, indices in factorNets) are EXCLUDED from the additive groups and appended as
   their own singleton trace groups at the tail of g — generated as traces but referenced only
   multiplicatively via lorFac, never summed into the integrand. nAdd marks the additive/factor
   boundary; pGroupOf maps a factor COMPONENT id to the list of {colour token, (0-based) trace-group
   ordinal} pairs it was split into — one per distinct token, see the GatherBy below. *)
      Module[{adj, fund, additivePos, factorPos = (# + 1)& /@ factorNets, gAdd, gFactor},
        additivePos = Complement[Range[Length[diagData]], factorPos];
        adj = Select[additivePos, diagTokExpr[[#]] === 1 && lorFacOf[[#]] === None&];
        fund = Select[additivePos, diagTokExpr[[#]] =!= 1 || lorFacOf[[#]] =!= None&];
        gAdd = Join[(# - 1)& /@ GatherBy[adj, diagData[[#]]&], List /@ (fund - 1)];
(* Each disconnected factor COMPONENT fuses into trace groups gathered by {component, COLOUR
   TOKEN}. A diag-dressed entry carries a runtime colour-sum token, which the assembly can only
   apply to a whole group; gathering by component alone would silently drop it. Per component:
       scalar = Σ_t token_t · traceRef[subgroup_t]     ( = Σ_entries token_e · trace_e )
   GatherBy is stable, so with a single token (all undressed flows) this yields one group per
   component, in order. factorNets are 0-based net indices. *)
        gFactor = GatherBy[factorNets, {factorCompOf[#], diagTokExpr[[# + 1]]}&];
        g = Join[gAdd, gFactor];
        nAdd = Length[gAdd];
(* One component maps to a LIST of {token, group ordinal} pairs, one per distinct token. *)
        pGroupOf = Merge[
            MapIndexed[
                (factorCompOf[#1[[1]]] -> {diagTokExpr[[#1[[1]] + 1]], nAdd + #2[[1]] - 1})&,
                gFactor],
            Identity]];
(* The dressing coefficient stays FACTORED in `diagData` (COEN CSEs it), so each group is one
   collected kinematic trace × its dressing — not a flat polynomial.
   Sum the ADDITIVE groups only (1..nAdd); factor groups (nAdd+1..) are referenced multiplicatively
   via lorFac, each computed ONCE as a separate trace. *)
      integrand =
        Sum[
          With[{rep = g[[gi, 1]]},
            diagData[[rep + 1]] * diagTokExpr[[rep + 1]] *
              If[lorFacOf[[rep + 1]] === None,
                1,
(* Each factor component contributes Σ_t token_t · traceRef[subgroup_t]. *)
                Times @@ (
                  Function[cid, Total[(First[#] traceRef[Last[#]])& /@ pGroupOf[cid]]] /@
                    lorFacOf[[rep + 1]])
              ] * traceRef[gi - 1]],
          {gi, nAdd}]];
(* ---- k-only dressing-lookup hoisting ("HoistLoopConstLookups") -----------------------------
   A lookup whose argument contains NO integration variable and NO grid coordinate is a LAUNCH
   CONSTANT: Zc[k] or ZA[(1+k^6)^(1/6)] is the same number for every thread of a map() launch,
   yet each thread pays the full coordinate transform (a fp64 log1p/log+asinh) plus the spline
   evaluation for it. Replace each DISTINCT such call with a scalar kernel parameter nthk<i>,
   evaluated once on the host by the generated static helper ntHoisted() (below) that the
   DiFfRG-side wrapper calls before launching.
   Applied to the integrand BEFORE the complex-branch split, so all #if branches see the same
   substitution and the kernel signature is branch-independent.
   NOT bit-identical (host libm vs device libdevice differ in the last ulp). Opt-in;
   MakeNTKernelDiFfRG enables it after checking the dressing types are DiFfRG interpolators. *)
    hoistCalls = {};
    hoistSyms = {};
    If[TrueQ[OptionValue["HoistLoopConstLookups"]] && dress =!= {},
      Module[{loopSyms = DeleteCases[args, Global`k], dressPat},
(* "Dressings" may arrive as strings (DiFfRG_compat derives them from the parameter list); in the
   integrand the calls carry SYMBOL heads, so normalise before matching. *)
        dressPat = Alternatives @@ (If[StringQ[#], Symbol["Global`" <> #], #]& /@ dress);
        hoistCalls =
          DeleteDuplicates @
            With[{loopPat = Alternatives @@ loopSyms},
              Cases[integrand, (d : dressPat)[a_] /; FreeQ[a, loopPat], {0, Infinity}]];
        If[hoistCalls =!= {},
          hoistSyms = Table[Symbol["Global`nthk" <> ToString[i - 1]], {i, Length[hoistCalls]}];
          integrand = integrand /. Thread[hoistCalls -> hoistSyms];
          ntLog["[khoist] hoisted ", Length[hoistCalls],
            " loop-constant dressing lookup(s) to host-evaluated kernel parameters"]]]];
(* $ntComplexRuntimeProjection is a package global read several call layers down
   (ntPureIntegrand/ntRePartIntegrand -> ntProjectIntegrand). Assigned UNCONDITIONALLY so a flow
   never inherits the previous flow's setting. *)
    endProject = TrueQ[OptionValue["ComplexEndProjection"]];
    If[endProject && !TrueQ[OptionValue["RealOutput"]],
      Message[MakeNTKernel::endproj];
      Abort[]];
    $ntComplexRuntimeProjection = TrueQ[OptionValue["ComplexRuntimeProjection"]] || endProject;
(* diagData is passed IN, not read from an enclosing scope — see ntkPruneSpec. *)
    With[{spec = ntkPruneSpec[diagData, g, complexQ, offline,
                   TrueQ[OptionValue["PruneRealTraces"]], TrueQ[OptionValue["RealProbe"]],
                   TrueQ[OptionValue["RunGenerator"]]]},
      pruneG       = spec["pruneG"];
      realOnlyG    = spec["realOnlyG"]];
    nGrp = Length[g];
(* The tarr declaration+fill, used by BOTH the kernel's coreBlock and the RealProbe TU (which
   evaluates the same integrand, so it needs the same tokens in scope). `trace_all_t` is emitted by
   emit_cpp_fused from the ACTUAL lowered roots (complex iff some trace is complex), so the array
   type can never disagree with what trace_all stores — and ntIm(double)=0.0 is then correct rather
   than lossy, because the type is double only when every trace really is real. *)
    tarrDecl = nsHome <> "::trace_all_t tarr[" <> ToString[nGrp] <> "]; " <> nsHome <> "::trace_all(fenv, tarr);";
(* surface the post-net-build shape, and abort on empty nets / empty grouping rather than emit a
   placeholder kernel. *)
    ntLog["[prof] post-net-build: nets=", Length[invNets], " groups(nGrp)=", nGrp, " complexQ=", complexQ];
    If[Length[invNets] === 0 || nGrp === 0,
      Message[mkGenerateKernel::emptynets, name, Length[invNets], nGrp];
      Abort[]];
    (* kinematic angle defs (kept symbolic in the dressing): emit once as named temporaries. *)
    angleDecls = ("const " <> $ntRealT <> " " <> SymbolName[First[#]] <> " = " <> cppFlat[Last[#]] <> ";")& /@ angleDefs;
(* NB: deliberately NO `using std::complex;` — unqualified complex<double> resolves to the
   support namespace's `complex` alias, so a device support header can substitute a device-safe
   complex (nvcc silently miscompiles std::complex arithmetic to 0 in device code).
   The fenv setup block: declare fenv, (dressed only) compute each dressing atom into dr_<id>,
   fill, and (CrossTraceCSE) precompute the traces. *)
    With[{hasDr = !FreeQ[invNets, _ntDressedCore]},
      Module[{coreBlock},
        coreBlock =
          {
            $ntRealT <> " fenv[(" <> nsHome <> "::nenv) > 0 ? (" <> nsHome <> "::nenv) : 1];",
            Sequence @@
              If[hasDr,
                KeyValueMap[
                  Function[{id, atom},
                    "const " <> $ntRealT <> " dr_" <> ToString[id] <> " = " <> cppFlat[atom] <> ";"],
                  $drTable],
                {}],
            With[{
              fillCallArgs =
                If[hasDr,
                  Join[SymbolName /@ fillArgs, ("dr_" <> ToString[#])& /@ Sort[Keys[$drTable]]],
                  SymbolName /@ fillArgs]},
              nsHome <> "::fill(fenv, " <> StringRiffle[fillCallArgs, ", "] <> ");"],
            If[crossCSE,
              tarrDecl,
              Nothing]};
(* DRESSED: the dr_<id> dressing expressions can reference the derived kinematic angles, so the
   angle (and colour) decls must precede the fenv block. NON-dressed: fenv first, then the
   angle/colour decls (keeps those kernels byte-identical). *)
        preamble =
          StringRiffle[
            If[hasDr,
              Join[ntSupportUsings[sns], angleDecls, coreBlock, preamble],
              Join[ntSupportUsings[sns], coreBlock, angleDecls, preamble]],
            "\n"]]];
    mkParam[nm_, ty_] := <|
        "Name" ->
          If[StringQ[nm],
            nm,
            SymbolName[nm]],
        "Type" -> ty,
        "Const" -> True,
        "Reference" -> True
      |>;
(* every dressing — including the named per-component diagonal dressings (ntSUNDiag{Fund,Adj}) —
   is an ordinary scalar interpolator kernel parameter. *)
    With[{
(* interpTy is normally one type string shared by every dressing. It may instead be an Association
   name -> type, for a flow that mixes 1-D momentum-grid interpolators with a 3-D vertex grid; a
   name not in the map (e.g. a NumTracer-internal ntSUNDiag dressing) falls back to the first
   declared type. *)
      dressTy =
        Function[nm,
          If[AssociationQ[interpTy],
            Lookup[interpTy, If[StringQ[nm], nm, ToString[nm]], First[Values[interpTy]]],
            interpTy]]},
(* the hoisted k-only lookup values ride at the END of the parameter list, so the DiFfRG wrapper
   can append them after the dressings without disturbing any existing argument position. The
   loop-independent constant() is called with the same argument tail (tuple_cat(pos, m_args)), so
   it must accept them too — unused there. *)
(* A scalar that is BOTH a runtime parameter and a frame coordinate (e.g. the temperature T at
   finite T, with an external leg pinned to vec[p,0] = pi T) must be declared once. Drop it from the
   ARGS side only: constParams and ntHoisted are also built from scalarParams, and fillArgs keeps
   the full args since the frame needs the symbol.
   Default order is scalars then dressings. A backend with a positional ABI supplies ParameterOrder;
   MakeNTKernelDiFfRG passes DiFfRG's Parameters order to match the integrator's forwarded tuple. *)
      sigArgs = DeleteCases[args, a_ /; MemberQ[scalarParamNames, ntParamName[a]]];
      runtimeParams =
        With[{runtimeNames = Join[scalarParams, dress]},
          With[{orderedEntries =
              If[parameterOrder === Automatic,
                runtimeNames,
                Join[
                  Select[
                    parameterOrder,
                    MemberQ[
                      ToString /@ runtimeNames,
                      ToString[If[AssociationQ[#], #["Name"], #]]
                    ] &
                  ],
                  Select[
                    runtimeNames,
                    !MemberQ[
                      ToString /@ (If[AssociationQ[#], #["Name"], #] & /@ parameterOrder),
                      ToString[#]
                    ] &
                  ]
                ]
              ]},
            Map[
              Function[entry,
                mkParam[
                  If[AssociationQ[entry], entry["Name"], entry],
                  ntRuntimeParamType[entry, adNames, scalarParamNames, dressTy]
                ]
              ],
              orderedEntries
            ]
          ]
        ];
      ntAssertADTyped[runtimeParams, adNames];
      kernelParams = Join[mkParam[#, $ntRealT]& /@ sigArgs, runtimeParams, mkParam[#, $ntRealT]& /@ hoistSyms];
(* The loop-independent `constant` is called by DiFfRG as constant(pos..., k, scalars..., dressings...),
   where pos is the FULL coordinate tuple of the flow's grid (quadrature_integrator.hh builds
   full_args = tuple_cat(coordinates.forward(idx), m_args)), so every grid coordinate must be a
   parameter (see constArgQ). `args` lists coordinates before k, so filtering preserves the order
   DiFfRG passes them in. *)
      constParams = Join[mkParam[#, $ntRealT]& /@ Select[args, constArgQ], runtimeParams, mkParam[#, $ntRealT]& /@ hoistSyms];
(* the host-side evaluator for the hoisted k-only lookups. The DiFfRG wrapper (patched by
   DiFfRG_compat.m) calls it once per map()/get() invocation and appends its results to the
   integrator call, in hoistSyms order. The lookups are plain `h(x)` calls: a DiFfRG interpolator's
   operator() picks the host or device buffer itself, so this un-decorated (host) function reads
   the host mirror. Same lowering as in-kernel, so results differ only in last-ulp rounding. *)
      hoistFnStr =
        If[hoistCalls === {},
          None,
          Module[{hkParams, vals},
            hkParams = Join[
              mkParam[#, $ntRealT]& /@ Select[args, # === Global`k&],
              runtimeParams];
            vals = (SymbolName[Head[#]] <> "(" <> cppFlat[#[[1]]] <> ")")& /@ hoistCalls;
            "static device::array<" <> $ntRealT <> ", " <> ToString[Length[hoistCalls]] <> "> ntHoisted(" <>
              StringRiffle[FunKit`MakeParameterString /@ hkParams, ", "] <> ")\n{\n  " <>
              StringRiffle[ntSupportUsings[sns], "\n  "] <> "\n  return {{" <>
              StringRiffle[vals, ",\n    "] <> "}};\n}"]];
(* dressed kernels: fill() takes one `double dr_<id>` per dressing atom — the kernel body computes
   the atom's value (regulators / interpolators in scope there) and passes it. Matches fm.dress. *)
      If[!FreeQ[invNets, _ntDressedCore],
        fillArgSig = fillArgSig <> StringJoin[(", [[maybe_unused]] " <> $ntRealT <> " dr_" <> ToString[#])& /@ Sort[Keys[$drTable]]]
      ]];
(* LOUD GUARD: the integrand must be numeric-valued before it is lowered to C++. A DEGENERATE input
   — most often a basis whose Gram is singular at the chosen kinematics, so its inverse metric (and
   hence every dual projector) carries 0/0 — leaves Indeterminate / ComplexInfinity / DirectedInfinity
   in the coefficients. FunKit's lowering prints those as bare identifiers (e.g.
   `return Indeterminate;`), and a differently-named leak could compile into a silently wrong
   kernel. Refuse rather than emit. *)
    With[{bad = Cases[integrand, Indeterminate | _DirectedInfinity | ComplexInfinity, {0, Infinity}]},
      If[bad =!= {},
        Message[MakeNTKernel::nonnumeric, Length[bad], Short[DeleteDuplicates[bad], 4]];
        Abort[]]];
(* FINITE-EXTENT PARTITION. The kernel is a flat sum of per-diagram terms, so classification is a
   Select, not a rewrite. Three outcomes:

     every term finite extent -> `matsubara_finite_extent` on the whole kernel (the pure-gauge case,
                                 and every flow while the quark is still 4D-regulated),
     no term finite extent    -> nothing; the Gaussian rule as before,
     MIXED                    -> emit BOTH halves as separate C++ entry points and let the
                                 integrator run each on the rule that suits it.

   The mixed case pays because term costs are very unequal: the expensive traces get the few exact
   modes, and only the unbounded terms run the full Gaussian rule. COEN's CSE is per-function, so
   each half computes only the traces and lookups it uses. *)
    With[{ms0 = mSym,
          hs0 = (If[Head[#] === Symbol, SymbolName[#], ToString[#]] & ) /@
                  Flatten[{Replace[OptionValue["DecayingRegulators"],
                    Automatic -> {"RB", "RF", "RBdot", "RFdot", "dq2RB", "dq2RF"}]}],
          forced = OptionValue["MatsubaraFiniteExtent"]},
      Which[
        mSym === None, Null,
        forced =!= Automatic,
          mFiniteExtentBody = TrueQ[forced],
        True,
          Module[{terms, feT, tlT},
            terms = If[Head[integrand] === Plus, List @@ integrand, {integrand}];
            {feT, tlT} = Lookup[GroupBy[terms, TrueQ[ntFiniteExtentQ[#, ms0, hs0]] &], {True, False}, {}];
            mFiniteExtentBody = (tlT === {}) && (feT =!= {});
            mSplit = (feT =!= {}) && (tlT =!= {});
            If[mSplit, feExpr = Total[feT]; tailExpr = Total[tlT]];
            ntLog["[matsubara] ", Length[terms], " term(s), ", Length[feT], " of finite extent",
              If[mSplit, " -- emitting a split kernel", ""]]]]];

(* the integrand -> C++ lowering (FunKit). Timed separately: it is the one heavy stage between the
   net-build and the generator emit, so without this the [prof] trail has a blind spot. *)
    With[{ntT =
      First @
        AbsoluteTiming[
(* the kernel body/bodies. A real flow has one. A COMPLEX one has three — the untouched complex
   form and the two real projections (ntPureIntegrand / ntRePartIntegrand) — or two under
   "RealOutput", spliced under an `#if` on the macro the probe writes into numtrace_verdict.hh.
   Which is valid depends on the trace VALUES, so the preprocessor picks; that is what lets
   generation run offline. Each body has its own MakeCppFunction so COEN's CSE spans it whole. *)
          mkKernelFnNamed = Function[{nm, expr}, ntShareInterpIndices[FunKit`MakeCppFunction[expr, "Name" -> nm, "Prefix" -> decor, "Return" -> "auto", "CodeParser" -> "Cpp", "Parameters" -> kernelParams, "Body" -> preamble], If[TrueQ[OptionValue["ShareInterpolatorIndex"]], dress, {}]]];
(* Per-body timing: the aggregate [prof] line also covers constFn, the class and the header, so the
   cost of one BODY (what "RealOutput" removes) is invisible in it. *)
          timedBodyNamed = Function[{nm, label, expr},
            Module[{t, res}, {t, res} = AbsoluteTiming[mkKernelFnNamed[nm, expr]];
              ntLog["[prof]   body ", nm, "/", label, ": ", t, " s"]; res]];
          realOut = TrueQ[OptionValue["RealOutput"]];
          bodyFor = Function[{nm, expr},
            If[!complexQ,
              timedBodyNamed[nm, "real", expr],
            If[endProject,
(* END-PROJECTION mode: build exactly one body, keep the complete assembled expression complex,
   and return its real part at the final C++ level. This avoids the symbolic Pure/RePart projection
   and the probe/verdict machinery. The price is runtime complex arithmetic; the benefit is much
   cheaper Mathematica lowering for finite-density denominators. *)
              timedBodyNamed[nm, "EndRe", Global`ntRe[expr]],
(* REAL-OUTPUT mode: emit the two real projections only. Verdict 0 then falls through to RePart, a
   TRUNCATION of the flow equation, which must never happen silently. The warning is a PREPROCESSOR
   #warning, not an ntLog, because offline the verdict is only known at `make numtrace` time. It
   rides on the MAIN body only, so a split flow does not print it three times. *)
              If[realOut,
                StringRiffle[Flatten @ {
                  "#if " <> verdictMacro <> " == 2   // Pure: the Complex -> Re projection is exact",
                  timedBodyNamed[nm, "Pure", ntPureIntegrand[expr]],
                  "#else                              // 1 = RePart; 0 = complex, truncated by RealOutput",
                  If[nm === "kernel",
                    {"#  if " <> verdictMacro <> " == 0",
                     "#    warning \"NumTracer: flow '" <> ns <> "' probed GENUINELY COMPLEX (verdict 0) but was generated with RealOutput -> True. The kernel returns only the real part; the imaginary part of the integrand is discarded. That is a truncation of the flow equation, not an identity. If it is not what you intended, regenerate without RealOutput and give the consumer a complex integrator.\"",
                     "#  endif"},
                    {}],
                  timedBodyNamed[nm, "RePart", ntRePartIntegrand[expr]],
                  "#endif"}, "\n"],
                StringRiffle[{
                  "#if " <> verdictMacro <> " == 2   // Pure: the Complex -> Re projection is exact",
                  timedBodyNamed[nm, "Pure", ntPureIntegrand[expr]],
                  "#elif " <> verdictMacro <> " == 1   // RePart: real value via complex trace(s), re/im split",
                  timedBodyNamed[nm, "RePart", ntRePartIntegrand[expr]],
                  "#else                              // the imaginary part survives: genuinely complex",
                  timedBodyNamed[nm, "Complex", expr],
                  "#endif"}, "\n"]]]]];
          kernelFn = bodyFor["kernel", integrand];
(* The two halves are lowered from the SAME machinery as the full body, so a complex flow gets its
   #if ladder in each of them and the verdict macro keeps meaning one thing across all three. *)
          splitFns =
            If[TrueQ[mSplit],
              {bodyFor["kernel_finite_extent", feExpr], bodyFor["kernel_tail", tailExpr]},
              {}];
          constFn = ntConstFn[OptionValue["Constant"], decor, constParams, sns];
(* MATSUBARA EVENNESS, Mathematica half. The generator proves it for the TRACES; this proves it
   for the rest of the kernel body (dressing/regulator arguments, denominators). The test is
   syntactic and conservative: strip every EVEN power of the symbol, then require it gone.
   `Sqrt[f0^2 + l1^2]` passes; a bare f0, or a shifted argument like ZQ[f0 + p0], does not. A
   missed trait only costs time, a wrong one costs correctness. This side decides whether to emit
   the member; the generator's constant supplies its value (absent member = false in DiFfRG). *)
          With[{evenFreeQ = Function[{e, ms}, FreeQ[e /. Power[ms, n_Integer /; EvenQ[n]] :> 1, ms]]},
            mEvenBody =
              mSym =!= None && evenFreeQ[integrand, mSym] &&
(* mVarIdx >= 0: the traces are the generator's job. Otherwise they can only see the frequency
   through a trace-env atom (symDefs), which no generator proof covers, so check those here. *)
                (mVarIdx >= 0 || AllTrue[Values[symDefs], evenFreeQ[#, mSym]&]);
            If[mSym =!= None && !mEvenBody,
              ntLog["[matsubara] kernel body uses ", mSym,
                " at an odd power (or inside a shifted dressing argument) — no matsubara_even trait"]]];
          If[mSym =!= None,
            ntLog["[matsubara] decaying regulators: ",
              Replace[OptionValue["DecayingRegulators"],
                Automatic -> {"RB", "RF", "RBdot", "RFdot", "dq2RB", "dq2RF"}],
              " — finite extent in ", mSym, ": ",
              Which[
                TrueQ[mSplit], "split — kernel_finite_extent on the exact sum, kernel_tail on the Gaussian rule",
                TrueQ[mFiniteExtentBody], "yes — emitting matsubara_finite_extent (exact Matsubara sum)",
                True, "no — the Matsubara sum keeps the Gaussian rule"]]];
(* ntRe/ntIm are needed by both real branches, so a complex flow always carries them. *)
          classStr = ntKernelClass[name,
            Join[
(* The generator emits its matsubara_even constant only when mVarIdx >= 0. Otherwise the traces
   are frequency-independent and the atom check above cleared the trace env, so the value is a
   literal true (referencing the absent constant would not compile). *)
              If[TrueQ[mEvenBody],
                {"static constexpr bool matsubara_even = " <>
                   If[mVarIdx >= 0, kns <> "::" <> ns <> "::matsubara_even", "true"] <> ";"},
                {}],
              If[TrueQ[mFiniteExtentBody],
                {"static constexpr bool matsubara_finite_extent = true;"},
                {}],
(* MIXED flow: two entry points instead of one. `kernel` stays and is still the whole thing -- it is
   what a consumer without the split machinery calls, and what the split is checked against. *)
              If[TrueQ[mSplit],
                Prepend[splitFns, "static constexpr bool matsubara_split = true;"],
                {}],
              {kernelFn, constFn},
              If[hoistFnStr === None, {}, {hoistFnStr}]],
            decor, regTemplate, regAlias, If[complexQ, {ntReImAccessors[decor]}, {}]];
          hdrInc = FileNameTake[headerFile];
          header =
            ntApplyTraceComplexOverride[
              FunKit`MakeCppHeader[
(* the numeric kernel is flat straight-line arithmetic: the generated trace functions (hdrInc) plus
   the support runtime; no tensor-engine headers. A complex flow pulls the verdict header unless
   ComplexEndProjection emits an unconditional end-real body with no probe/verdict. *)
                "Includes" ->Join[extraInc, ntRuntimeIncludes[runInc], {"numtracer/sun/sun_data.hpp", hdrInc}, If[complexQ && !endProject, {ntVerdictFile}, {}]], "Body" -> ntWrapBody[kns, classStr, name]
              ],
              hdrInc, kns, sns, complexQ
            ];]},
      ntLog["[prof] FunKit kernel/class/header lowering: ", ntT, " s"]];
    (* emit -> write -> compile -> run the generator. A local closure so the deferred
       PruneRealTraces pass (post-probe, below) can re-run it with realOnlyG updated. *)
    genPass[] := (
    With[{ntT = First @ AbsoluteTiming[{genPre, genUnits, genDecl, genMain} = emitNumericGenerator[invNets, invRest, colourNets, g, ncomp, ns, fillArgSig, kns, complexQ, realOnlyG, crossCSE, mVarIdx];]},
      ntLog["[prof] emitNumericGenerator: ", ntT, " s"]];
(* Split generator: a main TU + N net-builder unit TUs + a decl header, so the net builders compile
   in parallel (see emitNumericGenerator). The main `#include`s the decl. *)
    declFile = StringReplace[genFile, ".cpp" -> "_nets.hh"];
(* Precompiled-header source for the -O0 net-builder units. Deliberately a SUPERSET of what any
   one unit includes (a unit skips numeric_contract.hpp when the flow has no dressed nets, and
   sun_net.hpp when it is colour-free). The PCH is built once, so an unused header costs little. *)
    pchFile = StringReplace[genFile, ".cpp" -> "_pch.hh"];
    unitFiles = Table[StringReplace[genFile, ".cpp" -> "_u" <> ToString[u - 1] <> ".cpp"], {u, 1, Length[genUnits]}];
(* Writing the generator sources IS the work — and every write goes through ntExportCpp, whose leak
   scan aborts on a leaked head. Bound here rather than passed to ntLog for that reason. *)
    With[{ntT =
      First @
        AbsoluteTiming[
          ntExportCpp[declFile, genDecl];
          ntExportCpp[pchFile,
            "// GENERATED by MakeNTKernel — do not edit. Precompiled-header source for the -O0\n" <>
            "// net-builder units; see the NT_GEN_PCH guard in each unit TU.\n#pragma once\n" <>
            "#include \"numtracer/network/network.hpp\"\n#include \"numtracer/network/dirac.hpp\"\n" <>
            "#include \"numtracer/core/lit.hpp\"\n#include <utility>\n" <>
            "#include \"numtracer/numeric/numeric_contract.hpp\"\n#include \"numtracer/network/sun_net.hpp\"\n"];
(* each unit #includes the shared decl header so its net builders can call the cross-unit CSE
   accessors (lc<k>()/dc<k>()) and sibling net builders, parsed once per TU. *)
          Module[{uInc = "#include \"" <> FileNameTake[declFile] <> "\"\n"},
            Do[ntExportCpp[unitFiles[[u]], uInc <> genUnits[[u]]], {u, 1, Length[genUnits]}]];
          genSrc = genPre <> "\n#include \"" <> FileNameTake[declFile] <> "\"\n\n" <> genMain;
          ntExportCpp[genFile, genSrc];]},
      ntLog[
        "[prof] write generator files (", Length[unitFiles] + 2, " files, ",
        Round[(Total[StringLength /@ genUnits] + StringLength[genDecl] + StringLength[genMain]) / 1000000.],
        " MB): ", ntT, " s"]];
    Print["wrote generator: ", genFile, " (+ ", Length[genUnits], " net units + decl header)"];
(* run the generator at codegen time -> the committed straight-line traces header. stdout goes
   straight to the FILE via the shell (Run), not captured by RunProcess: headers can be ~40k lines.
   OFFLINE mode skips all of it: the `numtrace` CMake target compiles and runs the sources as a
   build step, and the committed traces header is left untouched until then. *)
    If[OptionValue["RunGenerator"] && !offline,
      incDir = OptionValue["IncludeDir"] /. Automatic :> resolveIncludeDir[];
      bin = FileNameJoin[{$TemporaryDirectory, "gen_" <> ns}];
(* COMPILE and RUN are timed separately (different levers). The main TU (mainOpt) and the -O0
   net-builder units compile CONCURRENTLY, then link. A failed unit compile leaves its .o missing,
   so the link rc is nonzero and the rc check below catches it. *)
      Module[
        {tcc, cc, mainObj, unitObjs, pcmd, lcmd, pchOut, pchCmd, pchArg, clog = bin <> "_compile.log", cxx = resolveGenCxx[], mainOpt, libPath = resolveGenLib[incDir], useLib, hoDef, libArg},
(* Default: link the prebuilt libNumTracer.a (engine bodies compiled once). If it is not found,
   fall back to a slower header-only compile (every engine body re-instantiated in the main TU). *)
        useLib = StringQ[libPath] && FileExistsQ[libPath];
        hoDef =
          If[useLib,
            " ",
            " -DNUMTRACER_HEADER_ONLY=1 "];
        libArg =
          If[useLib,
            " '" <> libPath <> "'",
            ""];
        ntLog[
          "[time]   generator engine: ",
          If[useLib,
            "linking " <> libPath,
            "header-only (libNumTracer.a not found)"]];
        mainObj = bin <> "_main.o";
        unitObjs = Table[bin <> "_u" <> ToString[u - 1] <> ".o", {u, 1, Length[unitFiles]}];
(* Main-TU optimisation level. Compile and run each happen once, so only their SUM matters.
   -O2 is dominated (much longer compile, same run). -O0 is NOT auto-selected: it wins on small
   non-dressed flows, but can cost minutes of run on dressed or dense-trace flows, and neither
   nSub nor the dressed flag predicts that reliably. So -O1 always; NT_GEN_MAIN_OPT=-O0 is the
   opt-in for flows known to be small and non-dressed. *)
        mainOpt =
          With[{e = Environment["NT_GEN_MAIN_OPT"]},
            Which[
              StringQ[e] && e =!= "",
                e,
              True,
                "-O1"]];
        ntLog["[time]   generator main TU: ", mainOpt, " (nSub = ", $ntGenNSub, "; NT_GEN_MAIN_OPT=-O0 is a large win on SMALL flows, but see the note above)"];
(* RAM-bounded parallel compile: at most $ntCompileJobs compiles at once (xargs -P), each capped at
   ~17 GB virtual (ulimit -v). Compiler output goes to `clog` (compile truncates, link appends) so
   genfail can quote the actual diagnostic.
   -fno-exceptions -fno-rtti on UNIT TUs only: at -O0 the exception-cleanup landing pads for their
   many destructible temporaries dominate the compile, and unit code never throws (NT_THROW degrades
   to abort()). The MAIN TU keeps exceptions: its thread-pool fallback uses try/catch.
   PRECOMPILED HEADER for the unit TUs (~5x less compile work per unit; the headers, not the
   emitted tables, dominate what the compiler parses). clang++ only. Without it the units fall back
   to textual includes via the NT_GEN_PCH macro guard, so it has no correctness surface. It MUST be
   built with the units' exact flag set (unitFlags); clang rejects a PCH whose flags disagree. *)
        With[{ccPre = "(ulimit -v 17000000; " <> cxx <> " -std=c++20 -ftemplate-depth=4000 ",
              unitFlags = "-O0 -fno-exceptions -fno-rtti" <> hoDef},
          pchOut = bin <> ".pch";
          pchCmd =
            If[StringContainsQ[cxx, "clang"],
              ccPre <> unitFlags <> "-I '" <> incDir <> "' -x c++-header '" <> pchFile <> "' -o '" <> pchOut <> "') > '" <> clog <> "' 2>&1",
              None];
          pchArg = If[pchCmd === None, "", " -DNT_GEN_PCH -include-pch '" <> pchOut <> "'"];
          pcmd = "printf '%s\\0' " <> StringRiffle[("\"" <> # <> "\"")& /@ Join[
              {ccPre <> mainOpt <> hoDef <> "-pthread -I '" <> incDir <> "' -c '" <> genFile <> "' -o '" <> mainObj <> "')"},
              Table[ccPre <> unitFlags <> pchArg <> " -I '" <> incDir <> "' -c '" <> unitFiles[[u]] <> "' -o '" <> unitObjs[[u]] <> "')", {u, 1, Length[unitFiles]}]], " "] <>
            " | xargs -0 -P " <> ToString[$ntCompileJobs] <> " -I CMD bash -c CMD >> '" <> clog <> "' 2>&1"];
        lcmd = cxx <> " -pthread '" <> mainObj <> "' " <> StringRiffle[("'" <> # <> "'")& /@ unitObjs, " "] <> libArg <> " -o '" <> bin <> "' >> '" <> clog <> "' 2>&1";
(* Content-addressed compile cache: the generator source is a deterministic function of the flow and
   the compile dominates the run, so unchanged sources+engine reuse the binary. The key covers the
   emitted sources, the linked libNumTracer.a, every installed engine header, and the full command
   lines (compiler, -O levels, flags). Deleting the .srckey file forces a rebuild. Relies on the
   emitted source being DETERMINISTIC; if it ever is not, key on the generator inputs instead. *)
        Module[{srcKey, keyFile = bin <> ".srckey", hit},
          srcKey =
            ToString @ Hash[
              {FileHash[#, "SHA256"]& /@ Join[{genFile, declFile, pchFile}, unitFiles],
               If[useLib, FileHash[libPath, "SHA256"], "header-only"],
               FileHash[#, "SHA256"]& /@ Sort[FileNames["*.hpp", incDir, Infinity]],
               pcmd, lcmd}, (* the command lines carry cxx, -O levels and every other flag *)
              "SHA256"];
          hit = FileExistsQ[bin] && FileExistsQ[keyFile] &&
            StringTrim[Quiet @ Check[ReadString[keyFile], ""]] === srcKey;
          If[hit,
            tcc = 0.;
            cc = 0;
            Print["[time]   generator compile: 0 s (cache hit: sources+engine unchanged, reusing ", bin, ")"],
            {tcc, cc} =
              AbsoluteTiming[
(* A failed PCH build is NOT fatal: drop the flags and let the units use their textual includes. *)
                If[pchCmd =!= None && Run[pchCmd] =!= 0,
                  ntLog["[warn]  PCH build failed; falling back to textual includes (see ", clog, ")"];
                  pcmd = StringReplace[pcmd, pchArg -> ""]];
                Run[pcmd];
                Run[lcmd]];
            If[cc === 0,
              Quiet @ Export[keyFile, srcKey, "Text"]];
            Print["[time]   generator compile (", cxx, ", ", Length[unitFiles], " parallel units + main): ", tcc, " s"]]];
(* quote the head of the compile log (ntLogHead bounds its length) *)
        If[cc =!= 0,
          Message[mkGenerateKernel::genfail,
            cxx <> " compile/link rc=" <> ToString[cc] <> "\n" <> ntLogHead[clog]];
          Abort[]]];
(* Free the codegen memo caches before launching the generator subprocess; they are only needed to
   emit the source. A small saving: most of the Wolfram kernel's resident memory is the `ntk`
   argument itself, which cannot be freed here. *)
      If[$NumTracerVerbose, ntLog["[prof] pre-run  MemoryInUse=", Round[MemoryInUse[]/1048576.], " MB  RSS=", Round[ntWolframRssMB[]], " MB"]];
      $ctCache = <||>; $odCache = <||>; $dsCache = <||>; $dslCache = <||>;
      If[$NumTracerVerbose, ntLog["[prof] post-free MemoryInUse=", Round[MemoryInUse[]/1048576.], " MB  RSS=", Round[ntWolframRssMB[]], " MB"]];
(* run into a TEMP file, validate (rc==0 AND non-empty), then move into place — so a crashed
   generator (e.g. thread-limited Run[]) never silently truncates the committed header. *)
      Module[{tmp = headerFile <> ".tmp", rc, sz, trun},
        {trun, rc} =
          AbsoluteTiming[
            Run[
              ntDeviceEnvPrefix[OptionValue["DeviceTarget"], decor] <> ntTcmallocPrefix[] <>
                "'" <> bin <> "' -n '" <> ns <> "' -d '" <> decor <> "' > '" <> tmp <> "'"]];
        Print["[time]   generator run (reduce+rebase+lower): ", trun, " s"];
        sz =
          If[FileExistsQ[tmp],
            FileByteCount[tmp],
            0];
        If[rc =!= 0 || sz < 64,
          If[FileExistsQ[tmp],
            DeleteFile[tmp]];
          Message[mkGenerateKernel::genfail, "generator run rc=" <> ToString[rc] <> " bytes=" <> ToString[sz] <> " (committed header left intact)"];
          Abort[]];
        CopyFile[tmp, headerFile, OverwriteTarget -> True];
        DeleteFile[tmp];
        Print["wrote header: ", headerFile, " (", sz, " bytes)"]]];
    ); (* end genPass *)
    genPass[];
(* semantic complexQ: the syntactic flag only says SOME coefficient carries an `i`; whether the
   assembled flow is actually real depends on the trace VALUES. The probe settles it against the
   generated traces and writes the verdict macro that selects one of the three bodies emitted above.
   Offline the same probe source is compiled and run by the `numtrace` build target instead. *)
    If[complexQ && !endProject,
      probeFile = FileNameJoin[{DirectoryName[genFile], "probe_" <> ns <> ".cpp"}];
      ntExportCpp[probeFile, ntProbeSource[integrand, args, fillArgs, angleDefs, angleDecls, nsHome, headerFile, $drTable, "TraceArrayDecl" -> If[crossCSE, tarrDecl, ""]]];
      Print["wrote probe: ", probeFile];
      If[TrueQ[OptionValue["RunGenerator"]] && TrueQ[OptionValue["RealProbe"]] && !offline,
        probeVerdict = ntRunProbe[probeFile, DirectoryName[headerFile], FileNameJoin[{DirectoryName[headerFile], ntVerdictFile}], verdictMacro];
        ntLog["[probe] verdict -> ", probeVerdict];
(* Deferred PruneRealTraces (pass 2): the verdict above was taken on the UNPRUNED traces, so a
   real verdict (Pure/RePart) certifies the imaginary residual cancels and the pruned
   re-generation is lossless for the consumer. A Complex verdict keeps all-complex traces. *)
        If[TrueQ[OptionValue["PruneRealTraces"]] && MemberQ[{"Pure", "RePart"}, probeVerdict] && MemberQ[pruneG, True],
          ntLog["[prune] PruneRealTraces: regenerating with ", Count[pruneG, True], "/", Length[pruneG],
            " real-coeff groups pruned (probe verdict '", probeVerdict, "' was taken on the unpruned traces)"];
          realOnlyG = pruneG;
          genPass[]]]];
    (* kernel header (write-if-changed). *)
    If[FileExistsQ[kernelFile] && Import[kernelFile, "Text"] === header,
      Print["unchanged: ", kernelFile],
      ntExportCpp[kernelFile, header];
      Print["wrote kernel: ", kernelFile]];
    (* per-flow numtrace manifest + switch. Written LAST, so a flow that aborted part-way leaves no
       manifest claiming to be buildable. Offline it says 0 (the numtrace target still owes the
       kernels); online everything is already done, so it says 1 and the target skips the flow. *)
    Module[{mf = ntWriteManifest[DirectoryName[kernelFile],
        <|"Class" -> name, "Namespace" -> ns, "Generator" -> genFile, "Traces" -> headerFile, "Units" -> unitFiles,
          "Decorator" -> decor, "MainOpt" -> mainOptForManifest, "Complex" -> complexQ, "Probe" -> probeFile,
          "DeviceTarget" -> OptionValue["DeviceTarget"]|>]},
      If[!offline, ntMarkGenerated[mf]];
      Print["wrote manifest: ", mf, If[offline, " (generated: 0 — run `make numtrace`)", " (generated: 1)"]]];
    <|"KernelFile" -> kernelFile, "HoistCount" -> Length[hoistCalls]|>]];

(* ---- MakeNTKernel: the public kernel emitter. --------------------------------------
   MakeNTKernel[ntk, genFile, kernelFile, tracesFile] emits the numeric matrix-product kernel:
   a build-time generator program (genFile, + net-builder units + decl header), run to produce the
   committed straight-line traces header (tracesFile), and the kernel header (kernelFile) that fills
   the fundamental symbols and calls the traces. Options are forwarded to the generator
   (see Options[mkGenerateKernel] for the set). *)
Options[MakeNTKernel] = {"ComputeType" -> "double", "Name" -> "nt_kernel", "Namespace" -> Automatic, "Dressings" -> {}, "ScalarParams" -> {}, "ADParams" -> {}, "ParameterOrder" -> Automatic, "Decorator" -> "static inline", "DeviceTarget" -> Automatic, "IncludeDir" -> Automatic, "RunGenerator" -> True, "AngleDefs" -> {}, "CrossTraceCSE" -> False, "Components" -> Automatic, "SymbolDefs" -> <||>, "RuntimeInclude" -> "numtracer/codegen/runtime.hpp", "ExtraIncludes" -> {}, "KernelNamespace" -> "numtracer_kernels", "SupportNamespace" -> "numtracer", "DressingType" -> Automatic, "ShareInterpolatorIndex" -> False, "HoistLoopConstLookups" -> False, "RegulatorTemplate" -> False, "RegulatorAlias" -> False, "RealProbe" -> True, "PruneRealTraces" -> False, "ComplexRuntimeProjection" -> False, "ComplexEndProjection" -> False, "RealOutput" -> False, "Constant" -> 0., "Offline" -> False, "CoordinateArgs" -> Automatic, "MatsubaraVar" -> None, "DecayingRegulators" -> Automatic, "MatsubaraFiniteExtent" -> Automatic};

MakeNTKernel::nfiles = "MakeNTKernel needs three output files: MakeNTKernel[ntk, genFile, kernelFile, tracesFile].";

MakeNTKernel[ntk : NTKernel[_], file_, opts : OptionsPattern[]] := (
    Message[MakeNTKernel::nfiles];
    Abort[]);

(* An UNKNOWN option name is refused, not ignored: `OptionsPattern[]` matches any rule, so a typo or
   a removed option would otherwise be swallowed silently. Checked at the public entry point, where
   the user's spelling arrives. *)
MakeNTKernel::optname = "Unknown option name(s) `1`. MakeNTKernel accepts: `2`. An unrecognised name is NOT applied — OptionsPattern[] matches any rule, so it would otherwise be swallowed silently and the setting would simply not take effect (this is what happened to \"Backend\" -> \"Dense\" after that option was removed). Check the spelling, or drop the option.";

ntAssertKnownOptions[opts_List] :=
  With[{unknown = Complement[Cases[opts, (nm_ -> _) | (nm_ :> _) :> ntParamName[nm]], ntParamName /@ Keys[Options[MakeNTKernel]]]},
    If[unknown =!= {},
      Message[MakeNTKernel::optname, unknown, Sort[ntParamName /@ Keys[Options[MakeNTKernel]]]];
      Abort[]]];

(* "ComputeType" picks the precision the EMITTED kernel runs in ("double" or "float", or a complex type
   of either, as for DiFfRG's MakeKernel); derivation and the generator stay in double. It is scoped to
   this one generation: the emitter strings read $ntRealT, and FunKit prints float literals. *)
MakeNTKernel::ctype = "\"ComputeType\" -> `1` is neither a double nor a float type.";

ntRealTypeOf[t_String] :=
  Which[
    StringContainsQ[t, "float"], "float",
    StringContainsQ[t, "double"], "double",
    True, Message[MakeNTKernel::ctype, t]; Abort[]];

MakeNTKernel[ntk : NTKernel[_], genFile_, kernelFile_, tracesFile_, opts : OptionsPattern[]] :=
  ntMakeNTKernel[ntk, genFile, kernelFile, tracesFile, opts]["KernelFile"];

(* MakeNTKernel's body, returning <|"KernelFile", "HoistCount"|>: MakeNTKernelDiFfRG needs the number
   of hoisted loop-constant lookups to patch its wrapper TUs. Explicit opts win (first match) and
   MakeNTKernel's own defaults fill the rest, so SetOptions[MakeNTKernel, ...] reaches the generator. *)
ntMakeNTKernel[ntk : NTKernel[_], genFile_, kernelFile_, tracesFile_, opts : OptionsPattern[MakeNTKernel]] := (
  ntAssertKnownOptions[Flatten[{opts}]];
  With[{realT = ntRealTypeOf[OptionValue[MakeNTKernel, {opts}, "ComputeType"]]},
    Block[{$ntRealT = realT,
           FunKit`Private`$codePrecision = If[realT === "float", "single", FunKit`Private`$codePrecision]},
      mkGenerateKernel[ntk, genFile, kernelFile, tracesFile, Sequence @@ FilterRules[Join[{opts}, Options[MakeNTKernel]], Options[mkGenerateKernel]]]]]);
