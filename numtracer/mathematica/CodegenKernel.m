(* ---- whole kernel (boilerplate delegated to FunKit) ------------------------- *)

(* ---- decorator normalisation: raw CUDA qualifiers -> Kokkos spelling -----------------------
   The device qualifiers are emitted verbatim onto every generated function (kernel/constant, the
   regulator wrappers, and — via the generator's `-d` flag — fill/trN/powr in the traces header).
   Spelling them `__host__ __device__` hard-codes the CUDA backend; the Kokkos macros expand to the
   right thing for CUDA/HIP/SYCL/OpenMP and to plain `inline` on a host-only build. DiFfRG links
   Kokkos unconditionally and decorates its own kernels this way. Applied to EVERY decorator (not just
   the Automatic one), so call sites that still pass the CUDA spelling keep working. Longest match
   first: the inline/__forceinline__ variants must be rewritten before the bare rule claims them. *)
ntKokkosDecor[decor_String] := StringReplace[decor, {
  "__host__ __device__ inline"          -> "KOKKOS_INLINE_FUNCTION",
  "__host__ __device__ __forceinline__" -> "KOKKOS_FORCEINLINE_FUNCTION",
  "__host__ __device__"                 -> "KOKKOS_FUNCTION"}];
ntKokkosDecor[d_] := d;

(* The Kokkos macros are NOT in scope in the standalone build-time programs (the probe compiles the
   generated traces header with a bare g++ and no Kokkos headers), so those sources must neutralise
   them. Emitted as the preamble of every such program. *)
ntKokkosStubDefs = "#define KOKKOS_INLINE_FUNCTION inline\n#define KOKKOS_FORCEINLINE_FUNCTION \
inline\n#define KOKKOS_FUNCTION\n#define __host__\n#define __device__\n";

(* regulator wrappers, prefixed with the user-chosen decorator (default plain "static inline";
   pass e.g. "Decorator" -> "static KOKKOS_INLINE_FUNCTION" for device-callable kernels).
   Emitted ONLY under "RegulatorTemplate" -> True (the fRG/DiFfRG shape); see ntKernelClass. *)

privDefs[decor_] :=
  StringRiffle[(decor <> " auto " <> # <> "(const auto &k2, const auto &p2) { return REG::" <> # <> "(k2, p2); }")& /@ {"RB", "RF", "RBdot", "RFdot", "dq2RB", "dq2RF"}, "\n"];

(* Real/imag accessors, emitted into the kernel class and into the probe. In the "RePart" assembly they
   are applied to trace tokens (double or nt_complex_t); under "ComplexEndProjection" to the whole
   integrand, which on an autodiff-active kernel is an autodiff type. They are therefore
   type-PRESERVING: a hard `-> double` would not compile on autodiff::Real (no .real()) and would drop
   the derivative where it did. The fallback is an UNQUALIFIED real(z)/imag(z), so for a DiFfRG
   consumer ordinary lookup finds DiFfRG::real/imag (common/complex_math.hh); it is instantiated only
   for a type without .real(). Trailing return types because kernel() calls these before their
   in-class definition, where g++ rejects a deduced `auto`. *)

ntReImAccessors[decor_] :=
  StringRiffle[{
    decor <> " " <> $ntRealT <> " ntRe(" <> $ntRealT <> " x) { return x; }",
    decor <> " " <> $ntRealT <> " ntIm(" <> $ntRealT <> ") { return " <> ntZeroLit[] <> "; }",
    "template <class T> " <> decor <> " auto ntRe(const T &z) -> decltype(z.real()) { return z.real(); }",
    "template <class T> " <> decor <> " auto ntRe(const T &z) -> decltype(real(z)) requires (!requires { z.real(); }) { return real(z); }",
    "template <class T> " <> decor <> " auto ntIm(const T &z) -> decltype(z.imag()) { return z.imag(); }",
    "template <class T> " <> decor <> " auto ntIm(const T &z) -> decltype(imag(z)) requires (!requires { z.imag(); }) { return imag(z); }"}, "\n"];

(* ---- the kernel CLASS. -----------------------------------------------------------------------
   Two shapes, one code path (it is emitted twice — once up front, once after the real probe
   re-lowers the integrand — and the two must not drift):

     regTemplate = False (DEFAULT, the general emission): a PLAIN class. NumTracer emits the kernel
       and nothing else; a flow whose dressing rules mention RB/RF/RBdot/RFdot/dq2RB/dq2RF emits
       UNQUALIFIED calls to them, which the consumer supplies (e.g. via "ExtraIncludes" -> {hdr},
       or from the "RuntimeInclude" header). Unqualified lookup runs class -> kernel namespace ->
       global, so free functions at global scope are found from any "KernelNamespace".

     regTemplate = True (the fRG/DiFfRG shape): `template<typename REG> class ...` plus the six
       private wrappers forwarding to REG::. DiFfRG's scaffold forward-declares the kernel as a
       template and instantiates it as KERNEL<Regulator>, so MakeNTKernelDiFfRG needs this.

   regAlias adds `using Regulator = REG;` (DiFfRG reads it back off the kernel class) and only makes
   sense with the template, which is why it implies it at the call sites below. *)

(* ---- interpolator index sharing --------------------------------------------------------------
   COEN hoists every dressing lookup into `const auto _interpN = <Interp>(<arg>);`. Many of those
   share an argument: a flow evaluates several dressings at the SAME momentum. The lookup is not
   cheap — measured at ~210 fp64 SASS instructions, of which ~200 are the fp64 log1p inside
   Coordinates::backward (the spline evaluation itself is ~10) — and the compiler CANNOT share it,
   because each interpolator owns its own `coordinates` members and cannot be proven to agree with
   another's. So a kernel with 13 lookups over 8 distinct arguments pays 13 transforms.

   This rewrites each such group to pay the transform once:

     const auto _interp3 = ZA3(<arg>);        ->   const auto _ix0 = ZA3.index(<arg>);
     const auto _interp7 = ZAcbc(<arg>);           const auto _interp3 = ZA3.at(_ix0);
                                                   const auto _interp7 = ZAcbc.at(_ix0);

   Measured on the YangMills flow set: 1.13-1.26x end to end, 14.5-24% fewer fp64 instructions,
   results bit-identical. See numtracer/gpubench/FINDINGS.md.

   Requires DiFfRG's SplineInterpolator1D/LinearInterpolator1D to expose index()/at(), where
   index() depends ONLY on the coordinate system (clamping lives in at(), since it is size
   dependent and would otherwise make an index untransferable between interpolators of different
   length). The interpolators of one flow are all built on the consumer's single coordinate object,
   which is what makes the sharing legal.

   Applied per kernel BODY, which is why a complex flow's three #if branches need no special care:
   each goes through its own MakeCppFunction, so a hoist can never escape into a sibling branch.

   LAYERING: whether a dressing handle offers index()/at() is a property of the CONSUMER's
   interpolator type, not of NumTracer — the index()/at() split is a DiFfRG pattern. So this pass is
   opt-in via "ShareInterpolatorIndex" and stays inert for every other backend; DiFfRG_compat.m
   turns it on once it has checked the dressing types. The callees are taken from the `dress` list
   the caller already supplies, so nothing here has to know how DiFfRG spells an interpolator type.
   Regulator calls (RB/RFdot/...) and, on complex flows, ntRe(ns::trK(fenv)) share the same
   `const auto _interpN = f(...)` shape and must NOT be rewritten — restricting callees to `dress`
   excludes them by construction (and ntRe is an overloaded member, so ntRe.index() would not even
   compile). *)

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

(* FINITE MATSUBARA EXTENT. A finite-T kernel is a sum of terms, each carrying exactly one
   dR/dt insertion, and it is that insertion -- not the propagators -- that decides how far the sum
   reaches in frequency. Where the insertion's ARGUMENT is coercive in the Matsubara symbol, the
   whole term dies super-polynomially in it, so `T Sum_n f(w_n)` is a FINITE sum: DiFfRG can
   enumerate the modes inside the regulator's support instead of approximating the infinite sum with a
   Gaussian rule. Below the crossover that is both exact and an order of magnitude cheaper.

   Deliberately NOT keyed on the function name. `RBdot`/`RFdot` are defaults, and one model can
   regulate four boson species through the same `RB` with four different arguments, or two fermion
   species 4D and 3D. The name only identifies "this factor is the insertion"; whether it confines
   the frequency is read off the argument, which is exactly the thing that differs. A 3D-regulated
   species writes `RB[k^2, l1^2]`, has no `f0` in the argument, and classifies itself as unbounded
   with no extra bookkeeping and no chance of the tag drifting out of sync with the algebra.

   The insertion's DECAY SHAPE is the one thing the argument cannot tell us: a Callan-Symanzik
   `R = k^2` confines nothing and `k^2/(1+x)` decays only algebraically, and against either the
   polynomial growth of the traces the term does not vanish at all. That is a property of the
   regulator, not of the diagram, so it is declared once per head via "DecayingRegulators" rather
   than inferred. *)
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
   integral (constExpr second arg of MakeKernel); NumTracer emits the same. A 0 constant needs no
   namespace usings (and stays byte-identical to the pre-Constant emission); a nonzero expr may call
   compute helpers (powr/pow/…) or the support API, so it gets the same usings the kernel body does. *)

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

(* ---- numeric matrix-product kernel, the MakeNTKernel "Numeric" path: emit the build-time
        generator program + the kernel header that calls its output. ---------------------- *)

Options[mkGenerateKernel] =
  {
    "Name" -> "nt_inv_kernel",
    "Namespace" -> Automatic,
    "Dressings" -> {},
    "ScalarParams" -> {},
    "ADParams" -> {},
    "ParameterOrder" -> Automatic,
    "IncludeDir" -> Automatic,
    "RunGenerator" -> True,
    "AngleDefs" -> {},
    "CrossTraceCSE" -> False,
    "Components" -> Automatic,
    "SymbolDefs" -> <||>,
    "Decorator" -> "static inline",
(* Does the emitted kernel target DEVICE code? This drives gen.hpp's size-gated `__noinline__`, which
   is a device-only lever (the host has no 255-register cliff). It must be stated, not sniffed from
   the decorator: `ntKokkosDecor` rewrites any raw `__host__ __device__` to the Kokkos macros before
   emission, so the old `__device__` string test never fired on a real flow — and testing for
   `KOKKOS_` instead would be wrong the other way, since those expand to plain `inline` on a
   host-only Kokkos build. Automatic infers from the raw CUDA spelling only (i.e. host unless the
   caller says otherwise); MakeNTKernelDiFfRG passes its own Device option through. *)
    "DeviceTarget" -> Automatic,
(* Name of the Matsubara-frequency symbol (a string or symbol), for a finite-T flow. Setting it asks
   the generator to PROVE whether the kernel is even in it and, if so, emit DiFfRG's
   `matsubara_even` trait — which lets QuadratureIntegrator_fT evaluate the kernel once per mode
   instead of twice. None = not a finite-T flow, no trait, no check. DiFfRG_compat passes the last
   entry of "IntegrationVariables", which is where DiFfRG's own integrators put it. *)
    "MatsubaraVar" -> None,
(* Regulator functions that DECAY super-polynomially in their momentum argument. A factor of one
   of these, at a coercive argument, is what gives a finite-T summand finite extent in the
   frequency; see ntFiniteExtentQ. Automatic = DiFfRG's six regulator wrappers.

   All six, not just the `dot` pair. A DRESSED insertion is dt(Z R) = Zdot R + Z Rdot, so half of
   every gluon/ghost term carries the regulator ITSELF and not its t-derivative:
     dressing[Rdot,{A,A},1,..] :> ZA[evP] RBdot[k^2,q2] + RB[k^2,q2] (dtZA[evP] + ..)
   Listing only RBdot/RFdot made that half unclassifiable and every flow came out unbounded --
   including the pure-gauge ones, whose every term does die with the regulator. RB dying is also
   what makes the term die, so it belongs here.

   Only a factor in a PRODUCT counts. `RB` inside a propagator denominator reaches ntFiniteExtentQ
   as Power[q2 + RB[..], -1], which it does not treat as vanishing -- correctly, since a decaying
   denominator makes the term GROW.

   A model whose regulator decays only algebraically (k^2/(1+x)), or not at all (Callan-Symanzik
   R = k^2), must list NOTHING here: then no `matsubara_finite_extent` trait is emitted and the
   Gaussian rule is kept. The ARGUMENTS are never configured -- they are read off each occurrence,
   which is what lets one `RB` serve a 4D- and a 3D-regulated species at once. *)
    "DecayingRegulators" -> Automatic,
(* Automatic = derive the `matsubara_finite_extent` trait from the algebra (the normal path). True/False
   force it, for a model that knows better than the syntactic test or wants to price the exact
   exact sum against the Gaussian rule on the same kernel. Forcing True on a kernel whose summand does
   NOT die above the regulator's support silently truncates the Matsubara sum. *)
    "MatsubaraFiniteExtent" -> Automatic,
    "RuntimeInclude" -> "numtracer/codegen/runtime.hpp",
    "ExtraIncludes" -> {},
    "KernelNamespace" -> "numtracer_kernels",
    "SupportNamespace" -> "numtracer",
    "DressingType" -> Automatic,
(* Opt-in: rewrite repeated dressing lookups to share one coordinate transform (see
   ntShareInterpIndices). Requires the consumer's interpolator handle to expose index()/at(),
   which is a DiFfRG property, so this stays False for backend-agnostic emission. *)
    "ShareInterpolatorIndex" -> False,
    "HoistLoopConstLookups" -> False,
    "RegulatorTemplate" -> False,
    "RegulatorAlias" -> False,
    "RealProbe" -> True,
    "PruneRealTraces" -> False,
    "ComplexRuntimeProjection" -> False,
    "ComplexEndProjection" -> False,
    "RealOutput" -> False,
    "Constant" -> 0.,
    (* "Offline" -> True: emit the generator + probe sources and a per-flow numtrace.json switch set to
       0, but do NOT compile or run anything — the `numtrace` CMake target does that as a build step,
       with make's parallelism across all flows. NT_OFFLINE in the environment overrides. *)
    "Offline" -> False,
    (* coordinate argument NAMES of this flow's grid; see constArgQ *)
    "CoordinateArgs" -> Automatic
  };

(* The dressing-collection path is driven by the `ntDressedNum` tokens NumTrace emits under
   "DressingCollection" -> True (no MakeNTKernel option needed): mkGenerateKernel auto-detects them and
   routes to the DPoly generator branch. *)
(* "Constant" -> expr (default 0.): the loop-INDEPENDENT piece that DiFfRG flat-adds to the
   integral — emitted verbatim as the body of `constant(p, k, dressings...)`, exactly like the
   constExpr second argument of DiFfRG's MakeKernel. A plain Mathematica expression in p/k and the
   dressing names (e.g. ZA[p] -> ZA(p)); NOT an NTKernel. Left at 0. the constant returns 0.. *)
(* "PruneRealTraces" -> False (default OFF): emit a `double` trace for any diagram group whose dressing
   coefficient is real (only Re(trace) is consumed), skipping the dead imaginary half. CORRECT for the
   kernel in every verdict, but HAZARDOUS with the probe: the probe runs on the generated traces and
   verifies Im(integrand)≈0; dropping a real-coeff group's imaginary half removes a term the probe
   needs, which can leave a non-cancelling residual and misclassify a real flow as Complex (blocking
   the double-kernel emission). Safe to enable only when the probe is off OR after per-flow validation,
   or once a post-probe trace regeneration is wired. Left off pending measurement (GCC often already
   DCEs the dead half of small inlined traces). *)
(* "ComplexRuntimeProjection" -> False (default OFF): what to do when a coefficient carries a Complex
   below a head the symbolic real/imaginary split cannot traverse — in practice a finite-DENSITY
   denominator, where the quark propagator carries l0 + I muq. See ntIiSafeQ for why the split is not
   merely slow there but silently wrong.
   OFF: refuse (MakeNTKernel::cplxnest). ON: keep such a coefficient intact and emit ntRe/ntIm around
   it, moving that one projection to generated-code RUNTIME. The C++ overloads already exist and
   FunKit lowers a Complex literal to complex<double>(0.,1.), so nothing new is needed on the C++
   side — but the kernel then performs complex arithmetic where it used to be all double, and the
   expression stays FACTORED (which is the point: expanding it is what defeats COEN's CSE and what
   made generation stall). Every summand takes the exact per-summand path in this mode, including
   the Pure body, because ntPureLinear's `/. Complex[a_, b_] :> a` has the same blind spot.
   Requires the type-generic powr (see ntProbeSource / the generator runtime): a denominator raised
   to a power is exactly what these flows contain. *)
$ntComplexRuntimeProjection = False;

(* "ComplexEndProjection" -> False (default OFF): for a real-valued consumer that wants
   Re[full complex integrand], skip the symbolic Pure/RePart projections entirely and emit one
   kernel body `ntRe[integrand]`. This keeps finite-density denominators complex until C++ runtime
   and projects only after the complete complex expression is assembled. It is algebraically the
   real part of the pointwise integrand, not a proof that the imaginary part cancels. Because the
   mode does not need a Pure/RePart verdict, it also skips the imaginary-part probe. Requires
   RealOutput -> True and is most useful together with ComplexRuntimeProjection -> True.

   COVERAGE: nothing in this repo sets it. DiFfRG_compat hardcodes False, no .wls passes True, and no
   fixture exercises the branch, so `endProject` and everything it guards are UNTESTED. That is a
   coverage gap, not dead code — the option is public and a consumer may set it — but treat the path
   as unverified until a fixture covers it. *)

(* "RealOutput" -> False (default OFF): the consumer's kernel return type. OFF, a complex flow emits
   all three bodies (Pure / RePart / untouched complex) and the probe's verdict picks one through the
   preprocessor. ON, the consumer is declaring it takes a REAL value, so the complex body — which such
   a consumer can never instantiate — is not lowered at all. That body is the expensive one (COEN CSEs
   the whole complex expression), and skipping it is the difference between a generation that finishes
   and one that appears to hang after the two real bodies are done.
   The catch, and why this is NOT the default and why MakeNTKernelDiFfRG does NOT set it either: with
   the complex body gone, verdict 0 — the probe saying the imaginary part genuinely survives — has
   nowhere to go but the RePart body. That yields Re[Integral[flow]], which is a TRUNCATION of the
   flow equation rather than an identity. It is frequently the right thing to want; it is never
   something to be given silently, so the emitted header carries a #warning in exactly that case.
   A consumer that wants it asks for it by hand. *)

(* "RealProbe" -> True (default): when the syntactic complexQ trips (some diagram coeff carries an `i`),
   compile+run a probe over the JUST-generated real traces to test whether Im(integrand) actually
   vanishes (projector-i × colour-i often cancel to a real value that Mathematica can't see through the
   opaque trace symbols). If it vanishes, re-emit a REAL (double) kernel — losslessly. Set False to skip
   the probe and always keep the complex+consumer-Re path. *)
(* The numeric (matrix-product) backend is the ONLY backend and is always on; there is no option for
   it. "Components" -> <|mom -> {e0,e1,e2,e3}, ...|> gives each momentum's 4 components as expressions
   (partially numeric / partially symbolic); Automatic falls back to the kinematic frame
   polynomialised. "SymbolDefs" -> <|sym -> expr|> gives the C++ fill for any DERIVED symbol
   (e.g. sin1 -> Sqrt[1-cos1^2]); plain free symbols are taken to be kernel arguments. *)
(* Standalone-output options (defaults make the emitted code self-contained against NumTracer's
   own headers, with no mention of any downstream consumer). A consumer that supplies its own
   support API points the codegen at it via these:
     "RuntimeInclude" -> "<hdr>" | None : the support header #included first, providing `complex`
        and `compute::{powr,pow,sqrt,fma}` (default numtracer/codegen/runtime.hpp; None to omit).
     "ExtraIncludes" -> {"a.hpp", ...}  : extra #includes prepended ahead of everything.
     "KernelNamespace" -> "ns" | None   : namespace wrapping the kernel class AND the generated
        trace functions (default "numtracer_kernels"; None emits at the includer's scope).
     "SupportNamespace" -> "ns"         : where `complex`/`compute` are looked up via `using`
        (default "numtracer").
     "DressingType" -> Automatic | "T"  : dressing-parameter type; Automatic emits `const auto&`
        (fully generic), or give a concrete type string.
     "RegulatorTemplate" -> False       : by default the kernel is a PLAIN class. A flow whose
        dressing rules use the regulators emits unqualified RB/RF/RBdot/RFdot/dq2RB/dq2RF calls,
        which the CONSUMER supplies — put them at global (or kernel-namespace) scope in a header and
        pull it in with "ExtraIncludes" -> {"my_regulators.hpp"} (or from the "RuntimeInclude"
        header). True restores the fRG/DiFfRG shape: `template<typename REG> class ...` plus private
        wrappers forwarding to REG::. Implied by "RegulatorAlias".
     "RegulatorAlias" -> False          : emit `using Regulator = REG;` in the class (DiFfRG reads
        the regulator type back off the kernel). Implies "RegulatorTemplate" -> True. *)
(* "RunGenerator" -> False: emit the generator sources only, without compiling/running them
   (the committed traces header is left untouched). *)
(* "Decorator" -> "<prefix>": the function prefix on EVERY emitted function — kernel/constant,
   the regulator wrappers, and (via the generator's runtime `-d` flag) fill/trN/powr in the
   straight-line header — e.g. "static __host__ __device__ inline" makes the whole kernel
   CUDA-device-callable. Default keeps the emitted bytes identical (kernel md5 invariant). *)
(* "CrossTraceCSE" -> True: lower all diagram trace polynomials through ONE shared CSE program
   (emitted as trace_all(f, t[]); the kernel fills t[] once and reads t[d]) so subexpressions are
   shared ACROSS traces — a fused kernel. Pays off when many traces share intermediates; default
   False keeps one independent trN() per trace.
   Works on COMPLEX flows: t[] is typed by the emitted `trace_all_t` (std::complex<double> iff some
   trace lowered complex), so nothing is truncated.
   Measured ZAqbq1_147 Mq-in (108 traces, 54 complex): 30,547 shared SSA instrs vs 47,558
   independent = 0.64x, lowering cost unchanged. The cross-diagram monomial duplication is 3.54x but
   that is NOT the attainable factor — each monomial occurrence still needs its own accumulate into
   its own trace; only the products are shared. Watch IPC as well as instruction count: this collapses
   N small functions into one very large basic block, which can spill (ZA "Route-B" shrank code 2.2x
   and regressed runtime 31%). *)
(* "AngleDefs" -> {sym -> expr, ...}: kinematic angle symbols the dressing keeps SYMBOLIC
   (e.g. cosl1p2 -> (-cos1 + Sqrt[3-3 cos1^2] cos2)/2). Emitted ONCE as `const double sym = ...;`
   in the kernel body so a shared sub-expression (the sqrt) is computed once rather than inlined
   per occurrence. *)

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

(* The C++ type of ONE runtime parameter (scalar or dressing), as a string for mkParam.

   Everything it needs is an explicit argument. That is the point: this decision used to be an
   inline If[] inside the runtimeParams Map, reading a `scalarParamNames` that a With[] two lines
   up had already closed over. MemberQ against an unbound symbol is False rather than an error, so
   the AD branch was simply dead and every AD scalar came out `double` — see MakeNTKernel::adtype
   for what that costs. A top-level function with named arguments cannot fail that way, and it is
   directly unit-testable (tests/test_ad_param_typing.wls) without running a generation.

   The AD test comes FIRST, before the scalarNames gate. Today adNames is a subset of scalarNames,
   but only as an accident of DiFfRG_compat building both from Type === "double" (:381 and :395);
   a caller passing ADParams without ScalarParams would otherwise silently get `double` back.

   An AD-flagged scalar's own entry says "Type" -> "double" — the AD flag lives in the separate
   "AD" key, which is why reading the type off the entry is not sufficient and never was. *)
ntRuntimeParamType[entry_, adNames_List, scalarNames_List, dressTy_] :=
  With[{nm = ntParamName[If[AssociationQ[entry], entry["Name"], entry]]},
    Which[
(* `auto` is not cosmetic: it makes the emitted function an abbreviated template, which is the only
   reason kernel()/constant() can bind autodiff::real from DiFfRG's integrator_AD twin. *)
      MemberQ[adNames, nm], "auto",
      MemberQ[scalarNames, nm], $ntRealT,
      AssociationQ[entry], entry["Type"],
      True, dressTy[nm]]];

(* POST-CONDITION on a finished mkParam list: every ADParams name that made it into the signature
   was actually typed auto. Split out from its call site so it is unit-testable, and because
   kernel(), constant() and ntHoisted all draw their scalars from the same runtimeParams list — one
   assertion covers all three signatures.

   This exists because every other stage of the pipeline is blind to the failure it catches:
   generation succeeds, net counts are unchanged, the emitted kernels are numerically identical and
   the ordinary get() compiles and runs. Only the consumer's AD twin fails, a full project build
   away. Returns the (unmodified) list so it can be used inline. *)
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

   WHY diagData IS AN ARGUMENT. This computation is the one that carried the worst bug in this file.
   `diagData` used to be declared in the INNER Module that builds the nets, while this read ran after
   that Module closed — so here it was the unassigned symbol `NumTracer`Private`diagData`, and
   `FreeQ[<unassigned>, Complex]` is vacuously True. Every group was therefore reported real,
   "PruneRealTraces" -> True emitted a kernel with every imaginary half dropped, and nothing said a
   word: net counts unchanged, generation clean, exit 0. Taking diagData as a PARAMETER makes that
   shape unconstructible, and the guards below make the remaining failure modes loud.

   ORDERING (the probe/prune interaction): the probe verifies Im(integrand)~0 against the GENERATED
   traces, so it must see the UNPRUNED all-complex set — probing pruned traces artificially removes
   the residual it needs and can misclassify the flow, and a wrong verdict is an O(1) kernel error.
   So when the probe is going to run, pass 1 generates unpruned, the verdict is taken, and only then
   is generation re-run with the prune applied (see the post-probe block; generation is seconds,
   correctness is not). Without a probe: a syntactically real flow prunes directly (no `i` anywhere
   — trivially safe), RealProbe->False prunes on the caller's assertion, and offline/no-generator
   cannot validate in-session so the prune request is dropped with a warning. *)

mkGenerateKernel::prunedata = "PruneRealTraces: the diagram-coefficient table has `1` entries but the trace grouping references diagram index `2`. The table and the grouping have gone out of step, so the per-group real/complex verdict below would be read off the wrong diagram — or off nothing at all. This is the shape of the bug this guard exists for (an unassigned diagData read as vacuously real).";

mkGenerateKernel::pruneall = "PruneRealTraces: the flow is COMPLEX, yet every one of its `1` trace groups was judged real and would have its imaginary half dropped. That is exactly what the historic bug produced from an unassigned coefficient table, and it is indistinguishable by inspection from a legitimately all-real grouping. If it is legitimate, the flow should not have tripped the complex test at all. Refusing to emit; re-run without \"PruneRealTraces\" -> True to get the full complex kernel.";

ntkPruneSpec[diagData_, groups_, complexQ_, offline_, pruneRequested_, realProbe_, runGenerator_] :=
  Module[{pruneG, realOnlyG, probeWillRun},
    pruneG =
      If[pruneRequested,
        Module[{maxDiag = Max[Append[#[[1]]& /@ groups, -1]]},
(* The index guard the original bug lacked: with diagData unassigned this Part would not even error,
   it would return an unevaluated Part expression that FreeQ then calls real. *)
          If[maxDiag >= Length[diagData],
            Message[mkGenerateKernel::prunedata, Length[diagData], maxDiag]; Abort[]];
          (FreeQ[diagData[[#[[1]] + 1]], Complex])& /@ groups],
        {}];
(* "every group is real" on a flow the complex test flagged is the historic bug's exact signature.
   It is also the one outcome no oracle can catch: the pruned kernel computes fewer numbers, all of
   them right, and only the DROPPED imaginary parts are wrong. Refuse it. *)
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

   READING THE Module LOCALS. Several are DECLARED here and ASSIGNED inside inner scopes — the note
   on `diagData` below records what that cost once. Where a stage has been extracted to a top-level
   function it now returns an Association instead, validated by ntStageResult, and the caller binds
   the fields; that is the direction the rest is moving. A local you cannot find an assignment for in
   THIS scope is assigned in an inner one, and that is the hazard, not a convention. *)
(* Block: the package globals assigned below ($RecursionLimit, the cache-key stamp, the dressing resolver,
   the canonicalisation rules, the complex-projection mode) are scoped to this one generation. *)
mkGenerateKernel[NTKernel[k_], genFile_, kernelFile_, headerFile_, OptionsPattern[]] :=
  Block[{$RecursionLimit = $RecursionLimit, $ctCtx = $ctCtx, $ntDressResolve = $ntDressResolve,
         $ntCanonIdsSrc = $ntCanonIdsSrc, $ntCanonRules = $ntCanonRules,
         $ntComplexRuntimeProjection = $ntComplexRuntimeProjection},
  Module[{name, ns, dress, scalarParams, adParams, parameterOrder, adNames, scalarParamNames, args, sigArgs, frame, env, nonzeroCompMask, ncomp, fillArgs, fillArgSig, constArgQ, invNets, invRest, g, colourNets, preamble, integrand, kernelParams, runtimeParams, constParams, mkParam, kernelFn, constFn, classStr, header, hdrInc, incDir, genPre, genUnits, genDecl, genMain, declFile, pchFile, unitFiles, genSrc, bin, complexQ, angleDefs, angleDecls, crossCSE, traceRef, nGrp, decor, tarrDecl, kns, sns, runInc, extraInc, interpTy, nsHome, regTemplate, regAlias, offline, realOut, endProject, verdictMacro, probeFile = None, mainOptForManifest, symDefs = <||>, realOnlyG = {}, pruneG = {}, probeVerdict = None, genPass,
(* hoistCalls/hoistSyms/hoistFnStr were NOT declared here at all — they were assigned unqualified
   and so became NumTracer`Private` globals that survive ACROSS generations. Every path assigns
   them, so there was no bug today; but that is the outer-declare/inner-assign hazard one level
   worse, and it is the same class as the `diagData` note above. Declared. *)
    hoistCalls = {}, hoistSyms = {}, hoistFnStr = "", mVarIdx = -1, mSym = None, mEvenBody = False, mFiniteExtentBody = False, mSplit = False, feExpr = 0, tailExpr = 0, splitFns = {}, mkKernelFnNamed, timedBodyNamed, bodyFor, dressedIdx = {}, diagTokExpr = {}, factorNets = {}, lorFacOf = {}, pGroupOf = <||>, nAdd = 0, factorCompOf = <||>,
(* diagData lives HERE, in the outer Module, not in the net-build Module below that assigns it.
   It used to be declared local to that inner Module (which spans the net-build loop and closes
   right after the `integrand` Sum), while `pruneG` reads it AFTER that close. Out of scope there,
   it stayed the unassigned symbol NumTracer`Private`diagData, so `diagData[[i]]` never evaluated,
   `FreeQ[..., Complex]` was vacuously True for EVERY group, and "PruneRealTraces" -> True emitted
   every trace as a real `double` regardless of its dressing coefficient — dropping -Im(c)*Im(tr)
   for any group with a complex coefficient. An O(1) wrong kernel whose only symptom was a
   Part::partd message. This is a DIFFERENT bug from the probe/prune ORDERING one fixed below; that
   fix reordered generation and never touched the scope. Same hand-off pattern as lorFacOf /
   pGroupOf / nAdd above: declared outer, assigned inner. *)
    diagData = {}},
    Needs["FunKit`"];
(* A large flow assembles a kernel with one summand per diagram GROUP (ZA4: 1274). Several codegen
   steps (the integrand Sum, COEN's expression lowering) recurse ~linearly in that count, so the
   default $RecursionLimit of 1024 is exceeded — and $RecursionLimit::reclim does NOT abort, it
   returns a held expression and the script continues to "DONE" having SILENTLY skipped the kernel
   (the ZA4 silent-skip). Raise the limit generously for the emission; a real runaway would still
   hit the (much higher) ceiling and surface. *)
    $RecursionLimit = Max[$RecursionLimit, 1048576];
    name = OptionValue["Name"];
    dress = OptionValue["Dressings"];
    scalarParams = OptionValue["ScalarParams"];(* loop-independent scalar doubles threaded into the signature *)
    parameterOrder = OptionValue["ParameterOrder"];
(* AD-flagged scalars (d1V, d2V for FE-potential flows) must be `const auto&` so the kernel also
   accepts autodiff::real from the integrator_AD twin; everything else stays `const double&`. *)
    adParams = OptionValue["ADParams"];
(* Both name lists are normalised ONCE, here, and both are Module locals. They used to be rebuilt
   at each use site inside a With[]; one of those With[]s closed before the use that needed it, and
   that is the whole of the AD-typing regression (MakeNTKernel::adtype). Nothing below rebinds
   either name. *)
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
       Mirrors the RunGenerator block's own choice (env override, else -O1 — -O2 is dominated). *)
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
(* Which of `args` the loop-independent constant() receives. Automatic keeps the historical
   p-and-k guess, correct for a 1-D grid; a caller that knows its grid passes the real coordinate
   names via "CoordinateArgs". k is always included. *)
    constArgQ =
      With[{ca = OptionValue["CoordinateArgs"]},
        If[ca === Automatic,
          Function[a, a === Global`p || a === Global`k],
          Function[a, a === Global`k || MemberQ[ca, If[StringQ[a], a, ToString[a]]]]]];
    ns = OptionValue["Namespace"] /. Automatic -> ToLowerCase[name];
    nsHome = kns <> "::" <> ns;(* where the generated trace fns / nenv / fill live *)
    verdictMacro = ntVerdictMacro[ns];(* the #if macro selecting one of the 3 complex-kernel bodies *)
    args = k["Args"];
    frame = k["Frame"];
    env = k["Env"];
    nonzeroCompMask = Association @ KeyValueMap[#1 -> frameMask[resolveComponents[#1, frame]]&, env];
    fillArgs = Select[args, # =!= Global`k&];(* scalars the fill needs *)
(* NUMERIC backend: build the component table over a compact parametrisation. Automatic uses the
   unit-loop spec (loop = magnitude × unit-direction symbols + unit constraint) so the contraction is
   as compact as the sp basis for BOTH Lorentz and Dirac nets. User "Components" are taken verbatim
   (polynomialised for any Sqrt/trig). *)
(* FRAME-SPEC PROBING + component table. unitLoopOkQ / unitLoopMixedOkQ / polyFrameSpec each run a
   PowerExpand + per-component Simplify sweep over the whole frame, and a general frame pays all
   three before numericComponents even starts. Timed as one block because that is the unit a flow
   either takes or skips. *)
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
                              picks out slot 0). Measured 70x on lambda3d against the fallback —
                              this branch is why the general-frame blowup went away.
     polyFrameSpec[frame]     the general fallback: every component a polynomial in the frame's
                              scalars, no unit constraint. Always correct, never compact.

   The ORDER is the specification: each test is strictly narrower than the next, so the first that
   qualifies is the most compact one available. NT_NO_UNIT_GROUPS disables both unit-loop branches, so
   the general polyFrameSpec is used; tests/gen/gen_lambda3d_small_numeric.wls builds its control that way. *)
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
(* The Matsubara frequency, resolved twice for two different consumers.

   `mSym` is the SYMBOL, and everything decided on THIS side of the fence keys off it: the evenness
   test on the integrand and the finite-extent partition further down only ever need to know which
   symbol to look for.

   `mVarIdx` is that symbol's MPoly variable index, and ONLY the generator needs it — it proves
   evenness of the TRACES from their monomial exponents. -1 means "not a momentum-component
   variable", which is NOT an error: a purely SCALAR integrand (a bosonic meson-potential tadpole,
   say) has an empty 4-vector component env, so `usyms` is empty and no frame symbol can be found
   there — while the kernel still depends on the frequency through its propagator denominators and
   regulator arguments, which live in the COEFFICIENT. Keying the traits off `usyms` alone is what
   used to silently drop `matsubara_finite_extent` on exactly those flows (measured on a sigma
   tadpole: 1 net, tr0 == 1, nenv == 0, every regulator argument f0^2 + l1^2 and so provably of
   finite extent — yet it fell back to the Gaussian rule and got no trait at all).

   The lookup therefore spans the frame's symbols AND the kernel's own fill arguments. *)
    mSym =
      With[{mv = OptionValue["MatsubaraVar"]},
        If[mv === None || mv === Automatic,
          None,
          With[{cands = Select[Join[ncomp["usyms"], fillArgs], SymbolName[#] === ToString[mv]&]},
            If[cands === {}, $Failed, First[cands]]]]];
(* A NAME THAT MATCHES NOTHING IS AN ERROR, not a quiet None. Silently doing nothing is how the
   __noinline__ gate stayed dead for months: the caller asked for an optimisation, got a valid
   kernel without it, and had no way to tell. Say so, loudly, and name the symbols that exist. *)
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
(* walk diagrams: each one is a Lorentz trace x a colour factor x a dressing/kinematic coeff. The
   bridge distribution split each 4-gluon vertex into colour channels, so MANY diagrams share the
   SAME coeff and differ only in (colour x Lorentz). Collect per diagram {coeff, lorNet, colNet},
   then GROUP by coeff: per group the generator folds the (constant) colour into the Lorentz poly
   and combines them — combinedTr_g = Σ_d colN_d · trN_d — so the kernel evaluates ~one polynomial
   per Feynman graph (≈5), not one per channel (51). Vanishing-colour diagrams drop out for free. *)
(* Dirac VERTEX: projector i + imaginary non-abelian colour (f^abc T^b T^c = (iN/2)T^a). Per
   diagram coeff*colour is real (the i's combine), so we keep the colour constant COMPLEX and
   take Re of the assembled integrand. *)
    complexQ = !FreeQ[k["Diagrams"], Complex];
(* CrossTraceCSE (IMPLEMENTED 2026-07-19; it was a documented stub before — `trace_all` existed only
   as the call site emitted below and in comments, so turning the option on produced a kernel that
   failed to compile with "'trace_all' is not a member of ...", and the complexQ guard that used to
   sit here was masking that).

   What it targets: cross-trace sharing is the one mechanism that attacks the measured redundancy in
   these flows — ZAqbq1_147 Mq-in carries 28,856 monomials across its traces of which only 8,144 are
   DISTINCT (3.54x; some recur in 24 different traces), and FormTracer's advantage here is exactly
   that it sums all diagrams into ONE polynomial before expanding, so those duplicates collect
   (measured Mq-out: FORM 5,060 monomials vs NumTracer 8,482, BOTH lowering at an identical 2.14
   ops/monomial). The compiler cannot do it for us: collecting one monomial out of N traces is a
   floating-point reassociation across function boundaries (force-inlining all 108 traces recovered
   only -3.5% instructions). See example/ZAqbq147-MqBench/FINDINGS.md.

   The old complex restriction is GONE. `tarr` is emitted as std::complex<double> whenever any trace
   is complex — the same type trN(fenv) returned — so the kernel's ntRe/ntIm reads are unchanged and
   nothing is truncated. Do NOT reintroduce a real tarr with per-trace phase tracking. *)
    preamble = {};
    factorCompOf = <||>;(* net -> factor-id *)
    $ctCache = <||>;(* clear the compileLorentz memo cache for this generation *)
(* dressedSlotStr / diracSlotStr / orderDiracLoops memos: all depend on generation-fixed state
   ($ntDressResolve, env, nonzeroCompMask, frame), so they MUST be cleared here alongside $ctCache.
   $dslCache used to be cleared only in the post-generation block below, which sits inside
   `If[RunGenerator && !offline, ...]` — i.e. never under the DiFfRG default (Offline -> True), so it
   accumulated across every flow of a multi-flow script. *)
    $dsCache = <||>;
    $odCache = <||>;
    $dslCache = <||>;
(* the generation-fixed half of every net-builder memo key, hashed ONCE here instead of on each
   (recursive) call. Covers everything those builders read besides the expression and its ids. *)
(* `nc` used to ride in this stamp too. It was a positional argument threaded through every net
   builder that NO body ever read — every SU(N) head carries its own rank N — and it was constant 0,
   so it contributed nothing to the hash either. Removed from the builders and from here. *)
    $ctCtx = Hash[{env, nonzeroCompMask, frame}];
    resetDiagDr[];(* clear the per-component diagonal-dressing registry for this generation *)
    resetDr[];(* clear the scalar-dressing (ntDressedNum) registry for this generation *)
(* frame resolver for dressed-numerator option coefficients (compileDirac → dressedSlotStr): the same
   ntSP/ntSPS/ntVec[q,i] → component substitution used for diag["Coeff"] below. *)
    $ntDressResolve =
      Function[s,
        s /. {ntSP[x_, y_] :> resolveComponents[x, frame] . resolveComponents[y, frame], ntSPS[x_, y_] :> Rest[resolveComponents[x, frame]] . Rest[resolveComponents[y, frame]], ntVec[q_, ii_Integer] :> resolveComponents[q, frame][[ii + 1]]}
      ];
(* The seven net-record accumulators are BAGS, not lists. `Append` on a list copies the whole list,
   so appending N records to seven of them is O(N^2) with N the NET count — which on a dense flow is
   an order of magnitude above the diagram count (2591 diagrams -> 82092 nets), and the copies also
   churn GBs through the allocator. Internal`Bag amortises the append; the lists are materialised
   once, after the loop. Nothing reads the accumulators DURING the loop except for the running net
   index, which `nNetAcc` now carries. *)
    (* NB no `diagData` in this local list — it is declared in the OUTER Module (see the note on its
       declaration there) because `pruneG` reads it after this Module closes. Re-adding it here
       re-shadows it and silently restores the all-groups-pruned bug. *)
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
(* Constant SU(N) component — colour and/or flavour, fundamental and/or adjoint: every group
   head carries its own rank N, so all fold through the SINGLE numeric SUNNet path
   (sun_value_cx, contracting each rank separately and multiplying). The fold is COMPLEX
   (-> Cx): when an imaginary non-abelian-vertex colour (f^abc T^b T^c = (iN/2) T^a) is folded
   in, the diagram's trace lowers to a complex `tr_i`; the kernel multiplies by the complex
   dressing coefficient and the consumer takes Re, so the imaginary part survives. ACCUMULATE
   components — a Yukawa loop carries SEVERAL constant components (a colour trace AND a flavour
   trace); mergeColNet folds them together so none is dropped (a bare `col = str` would keep
   only the last, e.g. the δ^ii = Nf factor, giving a ~50% wrong trace).
   ACCUMULATE THE FACTORS and compile once after the loop (see colBr below) rather than
   compiling per component and splicing the strings with mergeColNet: a constant component
   may be a PLUS (the four-quark Fierz flavour structure δδ - 4·T·T is one), which has no
   single-net representation, and one Expand of the whole constant product does the
   cross-product across several summed components in one step. Compiling the product is
   equivalent to the old per-component splice — mergeColNet only concatenates factor
   lists — so single-branch flows regenerate byte-identically. *)
                          constAcc = Join[constAcc, comp["Factors"]],
(* Non-constant component. The DISCONNECTED components of ONE diagram MULTIPLY (their scalar
   trace values `Times @@ toks`) — they are NOT separate summed diagrams. Collect them so
   the post-loop assembly can form the
   product (see there). Route by structure:
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
     - <= 1 non-constant factor: the EXISTING additive path. The single Dirac component's entries
       (colour folded, GlobalCollect-fusible) OR the single combined pure-Lorentz product net is
       appended with diagData = coeff*scal. Flows without a disconnected diagram regenerate
       BYTE-IDENTICAL (any K>=2 diagram previously aborted, so none is committed).
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
                      (* ---- single non-constant factor (or none): EXISTING additive path, byte-identical ---- *)
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
   Crucially this does NOT force the components' many entries into singletons (which would
   defeat colour-channel fusion and explode the trace count for a four-quark-dressed loop) —
   each component is exactly ONE trace (one scalar per component). *)
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
(* memo sizes: the net-build's cost is dominated by the DISTINCT Dirac/Lorentz structures, not the
   call count (a dense flow calls these tens of thousands of times for a few hundred distinct
   results). If a flow ever shows these growing with the call count, a memo has stopped hitting. *)
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
   RUNTIME colour-sum token (they cannot fuse by dressing coefficient alone, the token differs). lorFac diagrams (lorFacOf =!= None — a Dirac trace
   times a disconnected pure-Lorentz scalar) likewise stay singletons: each carries a per-diagram
   multiplicative trace, so it must not fuse with another diagram's entries.
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
(* each disconnected factor COMPONENT (possibly several colour-branch nets) fuses into ONE trace
   group, so its group trace = the component scalar (colour folded by the group sum). *)
(* Gather by {component, COLOUR TOKEN}, not by component alone. A factor group's trace is the SUM
   over its entries, and a diag-dressed entry carries a per-entry runtime colour-sum token — which
   the assembly below can only apply to a whole group. Grouping by component alone therefore had
   nowhere to put a factor entry's token and DROPPED it silently: the anchor multiplies in only its
   OWN diagTokExpr, never its factor components'. (Measured on the A_0 quark propagator of
   finite_T/QCD_Nf2/no_mesons: Zq lost 30 of 30 dressed nets, ZA 40 of 48.) Splitting the component
   by token instead makes each subgroup single-token, so the component scalar is recovered exactly as
       Σ_t token_t · traceRef[subgroup_t]     ( = Σ_entries token_e · trace_e )
   with the token factored OUT of each subgroup's trace. Entries sharing a token still share one
   trace, so a component with a single token (every undressed flow: all tokens are 1) yields exactly
   the one group it did before, in the same order — GatherBy is stable and a constant second key
   cannot regroup or reorder. Undressed kernels are byte-identical. *)
        gFactor = GatherBy[factorNets, {factorCompOf[#], diagTokExpr[[# + 1]]}&];(* factorNets are 0-based net indices *)
        g = Join[gAdd, gFactor];
        nAdd = Length[gAdd];
(* One component now maps to a LIST of {token, group ordinal} pairs, one per distinct token. *)
        pGroupOf = Merge[
            MapIndexed[
                (factorCompOf[#1[[1]]] -> {diagTokExpr[[#1[[1]] + 1]], nAdd + #2[[1]] - 1})&,
                gFactor],
            Identity]];
(* The dressing coefficient stays FACTORED in `diagData` (COEN CSEs it, like FORM's _repl), so each
   group is one collected kinematic trace × its dressing — not a flat polynomial. *)
(* Sum the ADDITIVE groups only (1..nAdd); factor groups (nAdd+1..) are referenced multiplicatively
   via lorFac. A factored (disconnected) diagram multiplies in its other components' scalars
   Π traceRef[factor groups] (each computed ONCE, as a separate trace): the diagram's
   `coeff * Times @@ component-scalars`. *)
      integrand =
        Sum[
          With[{rep = g[[gi, 1]]},
            diagData[[rep + 1]] * diagTokExpr[[rep + 1]] *
              If[lorFacOf[[rep + 1]] === None,
                1,
(* Each factor component contributes Σ_t token_t · traceRef[subgroup_t]; with one token (the
   undressed case) this is the bare traceRef it always was, since token_t == 1. *)
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
   NOT bit-identical: host libm log/exp differ from device libdevice in the last ulp, so this is
   gated at the physics level (observables + dressing sweeps), not bitwise.
   Opt-in (MakeNTKernelDiFfRG enables it after checking the dressing types are DiFfRG interpolators,
   unless its own "HoistLoopConstLookups" option says otherwise). *)
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
(* Same Private`-context handoff, for the same reason: ntProjectIntegrand runs several call layers
   down (the kernel body lowering -> ntPureIntegrand/ntRePartIntegrand) and threading an option through those
   would touch every one. Assigned UNCONDITIONALLY so a generation cannot inherit the previous
   flow's setting — this is a package-global, and MakeNTKernel is called once per flow in a script
   that generates many. See the "ComplexRuntimeProjection" note above. *)
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
(* [B2 localize] surface the post-net-build shape so a silent no-output (e.g. empty nets / empty
   grouping) is visible rather than appearing as a clean DONE. *)
    ntLog["[prof] post-net-build: nets=", Length[invNets], " groups(nGrp)=", nGrp, " complexQ=", complexQ];
    If[Length[invNets] === 0 || nGrp === 0,
      Message[mkGenerateKernel::emptynets, name, Length[invNets], nGrp];
      Abort[]];
    (* kinematic angle defs (kept symbolic in the dressing): emit once as named temporaries. *)
    angleDecls = ("const " <> $ntRealT <> " " <> SymbolName[First[#]] <> " = " <> cppFlat[Last[#]] <> ";")& /@ angleDefs;
(* NB: deliberately NO `using std::complex;` — unqualified complex<double> resolves to the
   support namespace's `complex` alias. That indirection lets a device/CUDA support header
   substitute a device-safe complex (std::complex arithmetic lowers to gcc _Complex builtins
   that nvcc silently miscompiles to 0 in device code), without changing the emitted kernel. *)
(* the fenv setup block: declare fenv, (dressed only) compute each dressing atom into dr_<id>, fill,
   and (CrossTraceCSE) precompute the traces. *)
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
   angle (and colour) decls must precede the fenv block. NON-dressed: keep the original order
   (fenv before angle/colour decls) so those kernels regenerate byte-identical. *)
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
   declared type, which is what collapsing to a single type used to do for everything. *)
      dressTy =
        Function[nm,
          If[AssociationQ[interpTy],
            Lookup[interpTy, If[StringQ[nm], nm, ToString[nm]], First[Values[interpTy]]],
            interpTy]]},
(* the hoisted k-only lookup values ride at the END of the parameter list, so the DiFfRG wrapper
   can append them after the dressings without disturbing any existing argument position. The
   loop-independent constant() is called with the same argument tail (tuple_cat(pos, m_args)), so
   it must accept them too — unused there. *)
(* A scalar that is BOTH a runtime parameter and a frame coordinate is declared once, by
   scalarParams. The finite-T case is the natural one: the temperature is a "double" kernel
   parameter (so DiFfRG_compat puts it in scalarParams — everything typed double that is not the
   special-cased k/p) AND it may appear in the frame, e.g. an external leg pinned to a Matsubara
   frequency vec[p,0] = pi T, which forces the caller to list it in "Args" so fill() receives it.
   Joining the two lists blindly then emits `const double& T, const double& T` and the kernel does
   not compile.

   Drop the duplicate from the ARGS side, not the scalarParams side: scalarParams is what
   constParams and the hoist function are built from as well, so removing it there would make
   constant() lose a parameter DiFfRG still passes it. args keeps its full form for fillArgs — the
   frame genuinely needs the symbol — so only the signature is de-duplicated.

   Backend-agnostic generation retains the historical scalar-then-dressing order. A backend with a
   positional ABI can supply ParameterOrder; MakeNTKernelDiFfRG passes DiFfRG's original Parameters
   order so interleaved scalar/interpolator packs match the integrator's forwarded tuple exactly. *)
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
   full_args = tuple_cat(coordinates.forward(idx), m_args)). Matching only `p` and `k` by name works
   for a 1-D grid whose coordinate happens to be called p; on a 3-D (S0,S1,SPhi) grid it silently
   DROPS all three coordinates, and a constant body referring to them then fails to compile with
   "identifier S0 is undefined". "CoordinateArgs" carries the real coordinate names; `args` already
   lists coordinates before k, so filtering preserves the order DiFfRG passes them in. *)
      constParams = Join[mkParam[#, $ntRealT]& /@ Select[args, constArgQ], runtimeParams, mkParam[#, $ntRealT]& /@ hoistSyms];
(* the host-side evaluator for the hoisted k-only lookups. The DiFfRG wrapper (patched by
   DiFfRG_compat.m) calls it once per map()/get() invocation and appends its results to the
   integrator call, in hoistSyms order. The lookups are written as plain `h(x)` calls: a DiFfRG
   interpolator's operator() is KOKKOS_FUNCTION and picks the host or the device buffer itself
   (KOKKOS_IF_ON_HOST / KOKKOS_IF_ON_DEVICE inside), so calling it from this un-decorated — hence
   host — function reads the host mirror with no explicit .CPU()/.GPU() selector; those members no
   longer exist. The expressions are the SAME CppForm lowering the kernel would have used, so
   semantics differ from the in-kernel evaluation only by host-libm-vs-libdevice last-ulp rounding. *)
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
   in the coefficients. FunKit's lowering happily prints those as bare identifiers, emitting e.g.
   `return Indeterminate;` into the kernel header. That is a Mathematica symbol in generated C++:
   here it fails to compile only by luck (no such identifier), and a differently-named leak could
   compile into a silently wrong kernel.
   Observed with the FULL AqbqDirect basis (12 structures) at the symmetric point: Det[g] = 0 there,
   because the 12 structures are linearly DEPENDENT at that kinematic configuration — which is why
   the flows project with a restricted sub-basis. Refuse rather than emit. *)
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

   The mixed case is the one that pays. A term's cost is dominated by the ONE trace it calls, and
   those are wildly unequal: ZA4's six terms carry tr0..tr5 at 2918/609/../313 lines, and the single
   unbounded term (the quark loop) is the 313-line one. Splitting puts the 2918-line trace on ~7
   exact modes and leaves only the cheap one running the full ~70-node Gaussian rule. COEN's CSE is
   per-function, so each half's body computes only the traces and interpolator lookups it actually
   uses -- that pruning is where the saving comes from, and it is automatic. *)
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
(* the kernel body/bodies.
   A real flow has exactly one. A COMPLEX one has three — the untouched complex form and the two real
   projections (ntPureIntegrand / ntRePartIntegrand) — or TWO under "RealOutput", which drops the
   complex form a real consumer could never instantiate anyway; spliced under an `#if` on the macro that
   the imaginary-part probe writes into numtrace_verdict.hh. Which is valid depends on the numerical
   values of the generated traces, so the choice cannot be made here; emitting all three and letting the
   preprocessor pick is what lets the whole generation run offline, as a build step. Each body goes
   through its own MakeCppFunction so COEN's CSE spans the whole expression, exactly as when Mathematica
   used to re-lower the single chosen one after the probe. *)
          mkKernelFnNamed = Function[{nm, expr}, ntShareInterpIndices[FunKit`MakeCppFunction[expr, "Name" -> nm, "Prefix" -> decor, "Return" -> "auto", "CodeParser" -> "Cpp", "Parameters" -> kernelParams, "Body" -> preamble], If[TrueQ[OptionValue["ShareInterpolatorIndex"]], dress, {}]]];
(* Per-body timing. The aggregate [prof] line above also covers constFn, the class and the header, so
   the cost of one BODY — which is what "RealOutput" removes — is invisible in it. Without this the
   only evidence for the option's value would be the downstream report's numbers. *)
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
(* REAL-OUTPUT mode: the consumer takes a double, so the untouched complex body can never be
   instantiated -- and it is the expensive one to lower (COEN CSEs the full complex expression). Emit
   the two real projections only. Verdict 0 then falls through to RePart, which is a TRUNCATION of the
   flow equation, not an identity: the imaginary part is discarded pointwise. It is a legitimate thing
   to ask for -- Re of the integral is often what the physics wants, and Im can integrate to zero over
   the loop angle even where it is nonzero pointwise -- but it must never happen silently, hence the
   #warning. It has to be a PREPROCESSOR warning rather than an ntLog: for a DiFfRG flow "Offline" is
   the default, so the probe runs at `make numtrace` time and the verdict is simply not known here.
   NOT the default, and MakeNTKernelDiFfRG does not set it either: opting into a truncation is the
   consumer's decision to make explicitly.
   The warning rides on the MAIN body only -- a split flow emits three functions from the same
   expression and the preprocessor would otherwise print the same warning three times. *)
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
(* MATSUBARA EVENNESS, second half. The generator proves it for the TRACES; this proves it for
   everything else the kernel body does with the frequency — dressing arguments, regulator
   arguments, denominators. Both must hold, and they are proven in different places because the two
   halves live in different languages: the traces are polynomials the C++ generator builds at build
   time, the rest is this Mathematica expression.

   The test is syntactic and conservative: strip every EVEN power of the symbol, then require the
   symbol to be gone. `Sqrt[f0^2 + l1^2]` passes (the inner f0^2 is stripped); a bare `f0`, or a
   fermionic dressing at a shifted argument like `ZQ[f0 + p0]`, does not. Erring towards "not even"
   is the safe direction — the trait only ever removes work, so a missed optimisation costs time
   while a wrong trait costs correctness.

   The two halves are combined in C++ rather than here: this side decides whether to emit the
   member at all, the generator's constant supplies its value. An absent member reads as false
   through DiFfRG's `requires K::matsubara_even` trait, which is exactly the fallback we want. *)
(* The even-power strip, used on the COEFFICIENT here and on the trace env's atom definitions
   below. `Sqrt[f0^2 + l1^2]` passes (the inner f0^2 is stripped); a bare f0, or a fermionic
   dressing at a shifted argument like ZQ[f0 + p0], does not. *)
          With[{evenFreeQ = Function[{e, ms}, FreeQ[e /. Power[ms, n_Integer /; EvenQ[n]] :> 1, ms]]},
            mEvenBody =
              mSym =!= None && evenFreeQ[integrand, mSym] &&
(* When mVarIdx >= 0 the trace side is the GENERATOR's job and its verdict is what the emitted
   member reads. When it is < 0 the traces carry no MPoly variable at all, so the only way they
   could still see the frequency is through a trace-env atom (a dressing evaluated at a shifted
   argument) — check those here, because no generator-side proof will be emitted to cover them. *)
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
(* The VALUE comes from the generator's monomial-exponent proof of the traces — but the generator
   only emits that constant when the frequency is an MPoly variable (mVarIdx >= 0). For a scalar
   integrand there is no such variable and no such constant: the traces are MPoly-constant in the
   frequency and the atom check above has already cleared the trace env, so the verdict is a
   literal true. Referencing the absent namespace constant here would simply fail to compile. *)
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
   ComplexEndProjection emits an unconditional end-real body with no probe/verdict. *)"Includes" -> Join[extraInc, ntRuntimeIncludes[runInc], {"numtracer/sun/sun_data.hpp", hdrInc}, If[complexQ && !endProject, {ntVerdictFile}, {}]], "Body" -> ntWrapBody[kns, classStr, name]
              ],
              hdrInc, kns, sns, complexQ
            ];]},
      ntLog["[prof] FunKit kernel/class/header lowering: ", ntT, " s"]];
    (* emit the generator source (the numeric matrix-product backend is the single generation path).
       The whole emit->write->compile->run pipeline is a local closure so the deferred
       PruneRealTraces pass (post-probe, below) can re-run it with realOnlyG updated. *)
    genPass[] := (
    With[{ntT = First @ AbsoluteTiming[{genPre, genUnits, genDecl, genMain} = emitNumericGenerator[invNets, invRest, colourNets, g, ncomp, ns, fillArgSig, kns, complexQ, realOnlyG, crossCSE, mVarIdx];]},
      ntLog["[prof] emitNumericGenerator: ", ntT, " s"]];
(* Split generator: a main TU + N net-builder unit TUs + a decl header (all in the tests/gen/ dir), so the
   net-builder codegen compiles in parallel (see emitNumericGenerator). The main `#include`s the decl. *)
    declFile = StringReplace[genFile, ".cpp" -> "_nets.hh"];
(* Precompiled-header source for the -O0 net-builder units. Deliberately a SUPERSET of what any
   one unit includes (a unit skips numeric_contract.hpp when the flow has no dressed nets, and
   sun_net.hpp when it is colour-free): the PCH is built once, so precompiling a header this flow
   does not use costs one 2 s build, whereas threading the exact per-flow set out of the emitter
   would mean plumbing it through the return value for no measured gain. *)
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
(* run the generator at codegen time -> the committed straight-line kernel header. The binary's
   stdout is redirected straight to the header FILE via the shell (Run), not captured in memory
   by RunProcess — the A4 kernel header is ~40k lines, and in-memory capture is fragile at that
   size.
   OFFLINE mode skips all of it: the sources are on disk and the `numtrace` CMake target compiles and
   runs them as a build step, with `make -j` across every flow at once instead of this per-flow xargs.
   The committed traces header is left untouched until then. *)
    If[OptionValue["RunGenerator"] && !offline,
      incDir = OptionValue["IncludeDir"] /. Automatic :> resolveIncludeDir[];
      bin = FileNameJoin[{$TemporaryDirectory, "gen_" <> ns}];
(* Time COMPILE and RUN separately — both count toward generation time, but the levers differ
   (compile: TU size / templates; run: reduce+rebase). Reported so neither is hidden. Compile the
   main TU (-O2) and the net-builder units (-O0) CONCURRENTLY (`&` + `wait`), then link — the unit
   codegen runs across cores. A failed unit compile leaves its .o missing -> the link rc is nonzero,
   so the existing rc check catches it. *)
      Module[
        {tcc, cc, mainObj, unitObjs, pcmd, lcmd, pchOut, pchCmd, pchArg, clog = bin <> "_compile.log", cxx = resolveGenCxx[], mainOpt, libPath = resolveGenLib[incDir], useLib, hoDef, libArg},
(* Default: link the prebuilt libNumTracer.a (engine bodies compiled once). If it is not found,
   fall back to a self-contained header-only compile so generation still works (older, slower
   path — every engine body re-instantiated in the main TU). *)
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
(* Main-TU optimisation level. Compile and RUN each happen exactly ONCE per generation, so the
   only thing worth minimising is their SUM — an optimisation level is not "safer" for being
   higher, it is just a different split of the same total.

   -O2 IS DOMINATED and is no longer offered. Measured on the four-quark Fierz gate (non-dressed,
   100 diagrams, 16223 distinct traces), main TU only:
       -O0   compile  2.75 s   run 9.30 s   total 12.1 s
       -O1   compile 20.43 s   run 0.60 s   total 21.0 s
       -O2   compile 47.40 s   run 0.55 s   total 48.0 s
   so -O2 buys 0.05 s of run for 27 s of compile. The dressed path was already measured the same
   way (28.7 s vs 17.5 s compile, IDENTICAL run) — which is why the dressed default was -O1. The
   non-dressed default was -O2 only because that TU was assumed cheap to compile; that stopped
   being true once the per-net data tables dominated it, and it is what made the Fierz gate look
   like it hung (>600 s at -O2 before those tables were hash-consed).

   -O0 is NOT auto-selected. It wins on the Fierz gate (12.1 s vs 21.0 s), so an nSub-gated
   "-O0 for small flows" rule looks obviously right — and it is WRONG twice over.

   First, nSub does not predict the run. Tried with a cutoff of 20000 it made the codegen suite
   regenerate ZA4_147 for >18 minutes (from seconds), and that flow has nSub = 294 — twenty
   times FEWER distinct traces than Fierz's 3831, but each vastly bigger (dense four-gluon
   vertex over the full tensor basis). The run scales with trace SIZE, not count.

   Second, the real discriminator is DRESSED vs not: the DPoly contraction loop must stay
   inlined, the MPoly one is far less sensitive. Measured on ZAqbq (all_tensors, dressed,
   nSub = 306), main TU only, identical 2007 KB kernel from all three:
       -O0   compile 2.37 s   run 34.26 s   total 36.63 s
       -O1   compile 5.11 s   run  1.31 s   total  6.42 s
       -O2   compile 8.94 s   run  1.27 s   total 10.21 s
   (a third independent confirmation that -O2 is dominated).

   But "-O0 when not dressed" is fragile too: Fierz is non-dressed and STILL loses 16x of run
   at -O0 (0.60 s -> 9.92 s); it only wins because its compile saving happens to be larger, and
   a denser non-dressed flow would flip that. The asymmetry decides it — -O1 costs at most a
   bounded ~18 s of compile, while -O0 can cost 5.7x (ZAqbq) or minutes (ZA4_147), and a slow
   run looks exactly like a hang. So: -O1 always, with NT_GEN_MAIN_OPT=-O0 available for flows
   known to be small and non-dressed, where it is a genuine ~2x win. *)
        mainOpt =
          With[{e = Environment["NT_GEN_MAIN_OPT"]},
            Which[
              StringQ[e] && e =!= "",
                e,
              True,
                "-O1"]];
        ntLog["[time]   generator main TU: ", mainOpt, " (nSub = ", $ntGenNSub, "; NT_GEN_MAIN_OPT=-O0 is a large win on SMALL flows, but see the note above)"];
(* RAM-bounded parallel compile: run at most $ntCompileJobs `cxx` at once (xargs -P), and cap
   each at ~17 GB virtual (ulimit -v). Peak RAM ~ jobs x per-unit; with the unit count scaled so
   each TU is small (~12 net-builders) this stays well under the machine limit for any flow size.
   Both phases redirect compiler stdout+stderr to `clog` (compile truncates, link appends) so the
   genfail message can quote the actual g++ diagnostic, not just the exit code.

   -fno-exceptions -fno-rtti (UNIT TUs only): each net-builder unit is thousands of braced-init-
   lists of destructible temporaries (DiracNet / DChainTok / DSlot / NetVal); at -O0 emitting the
   exception-cleanup landing pads for them is ~90% of the compile and superlinear in unit size
   (a 188 KB unit measured 15.3 s -> 1.3 s, RSS 660 MB -> 280 MB, with exceptions off). The unit
   code never throws/catches, and the library's internal guards route through NT_THROW
   (core/config.hpp), which degrades to abort() there — so this only removes dead cleanup code.
   The MAIN TU keeps exceptions: it emits a try/catch thread-pool fallback (see the parWork
   template below, `catch(const std::system_error&)`), which is ill-formed under -fno-exceptions. *)
(* PRECOMPILED HEADER for the -O0 unit TUs. Measured 2026-08-08 (za3_147, one unit, perf
   instructions:u): 11.10 G -> 2.29 G, i.e. 4.85x less compile work per unit, for a ~2 s one-off
   PCH build amortised over 8 (za3_147) to 74 (za4_147) units. This is the lever the table-size
   work was NOT: the emitted tables are ~3% of what the compiler parses (a 126 KB unit
   preprocesses to 3.34 MB), so shrinking them moved the compile 0%.
   clang++ only — GCC's PCH is a different mechanism (a .gch beside the header, found implicitly)
   and is not worth a second code path while resolveGenCxx[] prefers clang++ anyway. When the PCH
   is not built the units fall back to their textual includes via the NT_GEN_PCH *macro* guard, so
   this is a pure optimisation with no correctness surface. The NT_GEN_PCH macro is emitted into each
   unit TU and must stay.
   The PCH MUST be built with the unit TUs' exact flag set (-O0 -fno-exceptions -fno-rtti + hoDef);
   clang rejects a PCH whose flags disagree with the consumer's. *)
        (* ccPre: every compile shares it; unitFlags: the -O0 unit TUs and the PCH must use the SAME set *)
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
(* Content-addressed compile cache: the generator source is a deterministic function of the flow,
   and the compile dominates the run ~11:1 on medium flows (measured 2026-08-08 across za3_147 /
   aqbq147 / zaaqbq1_small), so a re-generation whose sources and engine are unchanged should not
   pay it again. The key covers the emitted sources, the linked libNumTracer.a, every installed
   engine header (the sources #include them), the compiler and the main-TU -O level — anything that
   can change the binary. Same pattern as FunKit's funkit-source.hash. Deleting the .srckey file is
   the opt-out.
   CAVEAT (documented, not speculative): if phase-B parallel lowering ever lands, sN interning makes
   the source non-deterministic and this key must move to the generator INPUTS. *)
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
(* ntLogHead is the shared bounded reporter — see there. This site is where the 50-line cap was first
   written; the probe and diagpoly compiles used to dump their logs whole. *)
        If[cc =!= 0,
          Message[mkGenerateKernel::genfail,
            cxx <> " compile/link rc=" <> ToString[cc] <> "\n" <> ntLogHead[clog]];
          Abort[]]];
(* Free the per-generation codegen memo caches BEFORE launching the generator subprocess (hygiene:
   they're needed only to EMIT the source, already on disk, and re-cleared next generation anyway).
   NOTE (measured 2026-07-22): this is MINOR — the four caches together are only ~50 MB. The real
   WolframKernel resident tax that co-resides with the generator (~2 GB at 655 diagrams, ~3.7 GB at
   3350) is the `ntk` ITSELF (the front-end Diagrams/Components analysis, a function argument), which
   is not freeable here. The effective RAM lever is instead FEWER diagrams — propagator collection
   folds 3350 -> 655 and the Wolfram tax with it. *)
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

(* Explicit opts win (first match); MakeNTKernel's own defaults carry through anything not passed —
   so `SetOptions[MakeNTKernel, ...]` (used by the test setups to opt into the shim + the DiFfRG
   namespace) propagates to the generator without editing every call site. *)

MakeNTKernel[ntk : NTKernel[_], file_, opts : OptionsPattern[]] := (
    Message[MakeNTKernel::nfiles];
    Abort[]);

(* An UNKNOWN option name is refused, not ignored. `OptionsPattern[]` matches any rule, so a name
   that is not in Options[MakeNTKernel] is silently swallowed — the caller's intent simply does not
   happen, generation succeeds, and the log says nothing. That is not hypothetical: "Backend" ->
   "Dense" was removed from the option list and three committed flow generators went on passing it
   for months, each believing it was selecting the dense backend. The same shape hides a typo in any
   option name. Checked here rather than in mkGenerateKernel because this is the public entry point,
   and it is where a user's spelling arrives. *)
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
   of hoisted loop-constant lookups to patch its wrapper TUs. *)
ntMakeNTKernel[ntk : NTKernel[_], genFile_, kernelFile_, tracesFile_, opts : OptionsPattern[MakeNTKernel]] := (
  ntAssertKnownOptions[Flatten[{opts}]];
  With[{realT = ntRealTypeOf[OptionValue[MakeNTKernel, {opts}, "ComputeType"]]},
    Block[{$ntRealT = realT,
           FunKit`Private`$codePrecision = If[realT === "float", "single", FunKit`Private`$codePrecision]},
      mkGenerateKernel[ntk, genFile, kernelFile, tracesFile, Sequence @@ FilterRules[Join[{opts}, Options[MakeNTKernel]], Options[mkGenerateKernel]]]]]);
