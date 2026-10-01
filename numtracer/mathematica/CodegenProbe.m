(* ---- semantic complexQ: probe whether the imaginary part actually vanishes -------------------
   The syntactic `complexQ = !FreeQ[Diagrams, Complex]` only sees that SOME diagram coefficient carries
   an `i` (projector i, non-abelian/quark colour f^abc T^bT^c = (iN/2)T^a). Those i's frequently pair up
   (i·i = -1) across diagrams so the assembled flow is exactly REAL — but the cancellation involves the
   concrete trace VALUES, and at the Mathematica stage the traces are opaque generated C++ symbols, so
   it can't be seen there. The just-generated C++ traces CAN see it: this compiles+runs a tiny probe that
   evaluates Im(integrand) over random frames with smooth real stub dressings/regulators. Returns True iff
   the imaginary part vanishes (|Im| <= tol·|Re|) over all sampled points, and writes the verdict that
   selects the real or complex kernel body (see below). *)
(* "TraceArrayDecl": with CrossTraceCSE the integrand's trace tokens are `tarr[i]` reads, not
   `ns::tr_i(fenv)` calls, so the probe TU must declare and fill that array exactly as the kernel's
   coreBlock does — otherwise the probe fails to compile and generation aborts (ntRunProbe).
   Empty for the per-trace path. *)

(* The verdict is applied by the PREPROCESSOR, not by Mathematica: `kernel.hh` carries all three bodies
   under `#if <MACRO> == 2 / #elif == 1 / #else` and the probe writes `<MACRO>` into the verdict header.
   That is what lets generation run offline — the probe is an ordinary build step, whereas the symbolic
   re-emission it used to trigger needed a Wolfram kernel. `ntProbeSource` builds the probe program (the
   verdict logic lives in its C++ `main`, which writes the header via `-o <file> -m <MACRO>`) and
   `ntRunProbe` compiles and runs it. *)
Options[ntProbeSource] = {"NPoints" -> 4000, "Tol" -> 1.*^-9, "TraceArrayDecl" -> ""};

(* ---- ntProbeSource: CONTRACT ------------------------------------------------------------------
   Emits the C++ source of the imaginary-part PROBE: a standalone program that evaluates the same
   integrand three ways over @p "NPoints" random kinematic points and prints the residuals the
   verdict is read from.

   In:   integrand   the assembled symbolic integrand (complex),
         args        every kinematic symbol it reads,  fillArgs  the subset the fill needs,
         angleDefs / angleDecls   the named angle temporaries, so the probe TU has them in scope,
         nsHome      the namespace the generated traces live in,
         headerFile  the traces header to include,  drTable  the dressing-atom table.
   Out:  the probe program text. It is compiled and run by the caller, which parses its stdout.

   The three bodies it compares are the untouched COMPLEX form and the two real projections
   (ntPureIntegrand — drop the imaginary coefficients; ntRePartIntegrand — the linear re/im split);
   the verdict is which of those reproduces the complex form to tolerance, i.e. 2 = Pure,
   1 = RePart, 0 = neither (keep it complex).

   THE CONSERVATISM IS THE POINT: a residual above tolerance resolves to 0, "keep the flow complex"
   (NaN sample points are skipped), and anything that stops the probe from producing a verdict — a
   compile that does not build, a failed run, output that does not parse — aborts generation in
   ntRunProbe instead of guessing. A wrong verdict is an O(1) kernel error, while an unnecessarily
   complex kernel is only slower. *)
ntProbeSource[integrand_, args_, fillArgs_, angleDefs_, angleDecls_, nsHome_, headerFile_, drTable_ : <||>, opts : OptionsPattern[]] :=
  Module[{keepHeads, keepSyms, seedOf, argComb, stub, probeFull, probeProj, probeRePart, probeParams, probePre, fnFull, fnProj, fnRePart, drDecls, drFillArgs, randDecls, callArgs, src, np, tol, distOf},
    np = OptionValue["NPoints"];
(* A float kernel's roundoff is ~1e-7 relative, so the double tolerance would read noise as a surviving
   imaginary part. 1e-4 still sits orders of magnitude below a genuine Im or a dropped term, which are
   O(1) relative. *)
    tol = If[ntSingleQ[], Max[OptionValue["Tol"], 1.*^-4], OptionValue["Tol"]];
(* GENERAL stubbing: replace EVERY external real-valued atom — any dressing (any arity), any
   regulator/support function, any named constant — with `ntStub(seed_head, hash(args))`, an
   INDEPENDENTLY-SEEDED pseudo-random real per head. "External" = head is not a structural math
   operation, not a frame arg, not a known constant, and not a trace token (those are raw C++ strings,
   left untouched). Independence is essential: a single shared stub would make `dress1[x]-dress2[x]`
   collapse to 0 and mask a real surviving imaginary part. Distinct heads -> distinct seeds; distinct
   arguments -> distinct hashes -> distinct values, so any genuine imaginary part is exposed for
   ARBITRARY dressings. *)
    keepHeads = Alternatives[Plus, Times, Power, Rational, Sqrt, Sin, Cos, Tan, Cot, Sec, Csc, Exp, Log, Abs, Sign, ArcTan, ArcSin, ArcCos, Sinh, Cosh, Tanh, Max, Min, Floor, Ceiling, Mod, Complex, List, Global`ntStub];
    keepSyms = Join[args, First /@ angleDefs, {Pi, E, EulerGamma, Degree, GoldenRatio}];
    seedOf[h_] := N[Mod[Hash[SymbolName[h]], 100003] + 7];(* distinct seed per external head *)
    argComb[a_List] := Total[MapIndexed[#1 * N[GoldenRatio] ^ (#2[[1]] - 1)&, a]];(* 0 for a constant *)
    stub[e_] := Module[{x},
        x = e //. (h_Symbol)[a___] /; FreeQ[keepHeads, h] && !MemberQ[args, h] :> Global`ntStub[seedOf[h], argComb[{a}]];
        x /. (s_Symbol) /; !MemberQ[keepSyms, s] && FreeQ[keepHeads, s] && Context[s] =!= "System`" && s =!= Global`ntStub :> Global`ntStub[seedOf[s], 0.]
      ];
    probeFull = stub[integrand];
(* The verdict-2 candidate must be the expression the Pure branch ACTUALLY EMITS, not a proxy for it.
   It used to be `probeFull /. Complex[a_,b_] :> a`, which drops the imaginary coefficients but leaves
   every trace token fully complex — a different expression from ntPureIntegrand's, which wraps each
   token in ntRe. For a LINEAR integrand the two agree numerically (Re of a sum of c·tr with real c),
   so no existing verdict moves; for a MULTILINEAR one they do not, and the old form certified a body
   that silently dropped the −Im(A)·Im(B) legs. Project the STUBBED integrand (stub leaves the raw
   C++ trace-token strings untouched, so the projections see the same tokens they will emit).
   probeRePart is the same idea for verdict 1, which was never checked against anything at all. *)
    probeProj = ntPureIntegrand[probeFull];
    probeRePart = ntRePartIntegrand[probeFull];
    probeParams = (<|"Name" -> SymbolName[#], "Type" -> $ntRealT, "Const" -> True, "Reference" -> True|>)& /@ args;
(* DRESSED kernels: the generated `fill()` takes one extra `double dr_<id>` per dressing atom (the
   kernel body computes each atom from regulators/interpolators and passes the VALUE). The probe has
   none of that runtime in scope, so it must compute each atom with the SAME pseudo-random stubbing
   used for the integrand and pass the values — otherwise the probe `.cpp` fails to compile (too few
   args to fill) and the kernel is conservatively kept COMPLEX, silently losing the lossless RePart
   double-kernel emission (the imaginary-half DCE that makes the distributed baseline fast). Atoms can
   reference the derived angles, so the angle decls precede them. *)
    drDecls =
      KeyValueMap[
        Function[{id, atom},
          "const " <> $ntRealT <> " dr_" <> ToString[id] <> " = " <> cppFlat[stub[atom]] <> ";"],
        drTable];
    drFillArgs = ("dr_" <> ToString[#])& /@ Sort[Keys[drTable]];
    probePre =
      StringRiffle[
        Join[
          angleDecls,
          {$ntRealT <> " fenv[(" <> nsHome <> "::nenv) > 0 ? (" <> nsHome <> "::nenv) : 1];"},
          drDecls,
          {nsHome <> "::fill(fenv, " <> StringRiffle[Join[SymbolName /@ fillArgs, drFillArgs], ", "] <> ");"},
          If[OptionValue["TraceArrayDecl"] === "",
            {},
            {OptionValue["TraceArrayDecl"]}]],
        "\n"];
(* THREE more full COEN lowerings of the (stubbed) integrand, on top of the one-to-three the kernel
   itself takes. Timed separately from the kernel's own [prof] body lines because it is pure
   verdict-machinery overhead for a flow whose verdict never changes. *)
    With[{ntT = First @ AbsoluteTiming[
    fnFull = FunKit`MakeCppFunction[probeFull, "Name" -> "probe_full", "Prefix" -> "static inline", "Return" -> "auto", "CodeParser" -> "Cpp", "Parameters" -> probeParams, "Body" -> probePre];
    fnProj = FunKit`MakeCppFunction[probeProj, "Name" -> "probe_proj", "Prefix" -> "static inline", "Return" -> "auto", "CodeParser" -> "Cpp", "Parameters" -> probeParams, "Body" -> probePre];
    fnRePart = FunKit`MakeCppFunction[probeRePart, "Name" -> "probe_repart", "Prefix" -> "static inline", "Return" -> "auto", "CodeParser" -> "Cpp", "Parameters" -> probeParams, "Body" -> probePre];]},
      ntLog["[prof] ntProbeSource: 3 body lowerings: ", ntT, " s"]];
    (* the frame's angle arguments are named cos<n> (a cosine, [-1,1]) and phi<n> (an azimuth); match
       the WHOLE name, so cosh…, phiField and the like are sampled as ordinary positive scalars *)
    distOf[a_] := Which[
        StringMatchQ[SymbolName[a], "cos" ~~ DigitCharacter ...],
          "Uc",
        StringMatchQ[SymbolName[a], "phi" ~~ DigitCharacter ...],
          "Uph",
        True,
          "U"];
    randDecls = StringRiffle[("double " <> SymbolName[#] <> " = " <> distOf[#] <> "(rng);")& /@ args, " "];
    callArgs = StringRiffle[SymbolName /@ args, ", "];
    src =
      StringJoin[
        (* the traces header decorates every trN/fill with the kernel decorator; this program is
           compiled standalone with a bare g++, so the device macros must be neutralised. *)
        ntKokkosStubDefs,
        "#include <complex>\n#include <cmath>\n#include <algorithm>\n#include <random>\n#include <cstdio>\n#include <cstring>\n",
        "#include \"" <> FileNameTake[headerFile] <> "\"\n",
(* Generic in the ARGUMENT type, not just the exponent. The probe evaluates the integrand itself, so
   under "ComplexRuntimeProjection" it raises a COMPLEX denominator (l0 + I muq) to a power and a
   `double powr(double)` fails to compile — with the whole compiler log, previously imported in full,
   as the only diagnostic. numtracer::compute::powr (runtime.hpp) is already generic this way; this
   is the probe's own private copy catching up. T(1) rather than 1.0 for the identity and inverse. *)
        "template<int N, class T> static inline T powr(T x){ T r=T(1); int n=N<0?-N:N; for(int i=0;i<n;++i) r*=x; return N<0?T(1)/r:r; }\n",
(* min/max belong here for the same reason the rest do: ntProbeSource's keepHeads deliberately KEEPS
   Max and Min, so a kernel containing one reaches the probe -- where, unlike the kernel, there is no
   `using namespace DiFfRG` to supply them. The kernel compiled and only its probe did not, with
   "'min' was not declared in this scope" as the whole diagnostic. They come from <algorithm>. *)
        "using std::pow; using std::sqrt; using std::sin; using std::cos; using std::tan; using std::exp; using std::log; using std::fma; using std::fabs; using std::min; using std::max;\n",
        With[{c = "std::complex<" <> $ntRealT <> ">"},
          "static inline " <> c <> " fma(const " <> c <> "&a,const " <> c <> "&b,const " <> c <> "&c){return a*b+c;}\n"],
        "template<class T> using complex = std::complex<T>;\n",
(* independently-seeded pseudo-random real in [0.4,0.9): same (seed,arg) -> same value (a dressing is
   a function), distinct (seed,arg) -> independent value, so no two dressings or arguments collide.
   Single precision uses a SMOOTH stub in [0.4,0.9] instead: the hash amplifies its argument ~1e4x,
   and the bodies being compared compute that argument in different orders, so float roundoff in it
   became O(1) different dressings per body (ZA3 with mesons: 45% of points "disagreeing"). Distinct
   seeds still give distinct, generic values. *)
        If[ntSingleQ[],
          "static inline float ntStub(double seed, double x){ return float(0.65 + 0.25*std::sin(seed*0.1031 + x*0.3127 + 1.7)); }\n",
          "static inline double ntStub(double seed, double x){ double h = std::sin(seed*0.1031 + x*0.3127 + 1.7)*43758.5453; return 0.4 + 0.5*(h - std::floor(h)); }\n"],
(* both real projections call ntRe/ntIm on the trace tokens, exactly as the kernel does *)
        ntReImAccessors["static inline"], "\n",
        fnFull,
        "\n",
        fnProj,
        "\n",
        fnRePart,
        "\n",
        "int main(int argc, char** argv){\n",
        "  const char* outf=nullptr; const char* macro=nullptr;\n",
        "  for(int i=1;i<argc;++i){ if(!std::strcmp(argv[i],\"-o\") && i+1<argc) outf=argv[++i];\n",
        "                           else if(!std::strcmp(argv[i],\"-m\") && i+1<argc) macro=argv[++i]; }\n",
        "  std::mt19937_64 rng(12345); std::uniform_real_distribution<double> U(0.25,3.0),Uc(-0.9,0.9),Uph(0.1,6.2);\n",
        "  double mim=0,mdiff=0,mre=0,mrim=0,mrdiff=0,mrrep=0; long ok=0, nrep=0;\n",
        "  for(int n=0;n<" <> ToString[np] <> ";++n){ " <> randDecls <> "\n",
        "    std::complex<double> f = probe_full(" <> callArgs <> "); double pj = probe_proj(" <> callArgs <> ");\n",
        "    double rp = probe_repart(" <> callArgs <> ");\n",
        "    double im=std::imag(f), re=std::real(f), df=std::abs(f-pj);\n",
(* COUNT the points that disagree, do not just take the max. RePart and the complex body are
   algebraically equal but associated completely differently, so at a point where the sum cancels
   badly their roundoff diverges — one such point in a few thousand is normal and says nothing about
   the emitter. (Measured on ZAAqbq2: median 1.0e-16, p99 1.3e-15, max 7.9e-9 — a single outlier at
   large |f|, seven orders above the 99th percentile.) The defect this guards is a property of the
   EXPRESSION, not of a point: a double-counted or dropped term is wrong at essentially EVERY point.
   So the discriminator is the FRACTION exceeding the tolerance, not the worst case. *)
        "    { double rr = std::fabs(rp-re)/(std::abs(f)+1.0);\n",
        "      if(std::isfinite(rr)){ mrrep=std::max(mrrep, rr); if(rr > " <> ToString[CForm[N[tol]]] <> ") ++nrep; } }\n",
(* PER-POINT relative measures (mrim, mrdiff): a global max|Im|/max|Re| can let a localized
   imaginary part hide behind a large |Re| at some OTHER point (catastrophic cancellation). The
   +1 floor degrades gracefully to an absolute test when |Re|/|f| are small. The verdict keys on
   these; the absolute trio is kept only for the log. *)
        "    if(std::isfinite(im)&&std::isfinite(re)&&std::isfinite(df)){ mim=std::max(mim,std::fabs(im)); mdiff=std::max(mdiff,df); mre=std::max(mre,std::fabs(re));\n",
        "      mrim=std::max(mrim, std::fabs(im)/(std::fabs(re)+1.0)); mrdiff=std::max(mrdiff, df/(std::abs(f)+1.0)); ++ok; } }\n",
(* Three-way verdict, decided HERE (in the C++ that resolves every complex multiplication) so the whole
   probe is one self-contained build step:
     0 "Complex"  Im survives                          -> genuinely complex, keep it.
     2 "Pure"     Im=0 AND Complex->Re projection exact -> drop imaginary coeffs (clean real arithmetic).
     1 "RePart"   Im=0 but projection differs           -> a trace is itself complex; the value is real
                                                           but only `.real()` of the full complex result
                                                           is correct, so the re/im split applies.
   Keyed on the PER-POINT relative measures (mrim, mrdiff): a global max|Im|/max|Re| can let a localized
   imaginary part hide behind a large |Re| at some OTHER point (catastrophic cancellation). No usable
   points is NOT a quiet "Complex" any more — it would bake the wrong branch into a committed header —
   so it exits nonzero and the caller aborts. *)
        "  if(ok < 1){ std::fprintf(stderr, \"[probe] no usable points\\n\"); return 2; }\n",
        "  const int verdict = (mrim > " <> ToString[CForm[N[tol]]] <> ") ? 0 : ((mrdiff <= " <> ToString[CForm[N[tol]]] <> ") ? 2 : 1);\n",
        (* the six measures + the RePart outlier count + point count + verdict, on one line. Printed
           BEFORE the RePart gate below so a failure still reports every measure the caller logs —
           returning early left the operator with a bare "exit 4" and no numbers to judge it by. *)
        "  std::printf(\"%.10e %.10e %.10e %.10e %.10e %.10e %ld %ld %d\\n\", mim, mdiff, mre, mrim, mrdiff, mrrep, nrep, ok, verdict);\n",
(* RePart is not a candidate to be chosen between — it is an IDENTITY, Re(Σ c·tr) rewritten with
   ntRe/ntIm accessors, and it must reproduce real(probe_full) regardless of which verdict wins.
   Nothing checked that before, which is how the linear decomposition could double-count a multilinear
   term (and emit unresolved tr$ placeholders) without any test noticing.
   Gated on the FRACTION of disagreeing points, not the maximum: see the counting comment above. The
   1% cut sits orders of magnitude away from both sides — a well-conditioned flow has 0 outliers, a
   badly-conditioned one a handful (ZAAqbq2: 1 in 4000 = 0.025%), and a genuine emitter bug misses at
   EVERY point (100%). It cannot see a systematic error smaller than the tolerance, which is a real
   limit of the method rather than an oversight: below that it is under a flow's own roundoff floor. *)
        "  if(nrep * 100 > ok){ std::fprintf(stderr, \"[probe] the RePart projection does not reproduce Re(integrand): %ld of %ld points disagree by more than " <> ToString[CForm[N[tol]]] <> " (worst rel=%.3e).\\n\"\n",
        "      \"[probe] A few isolated outliers would be catastrophic cancellation; this many is a NumTracer emitter bug (ntRePartIntegrand).\\n\", nrep, ok, mrrep); return 4; }\n",
        "  if(outf && macro){ std::FILE* f = std::fopen(outf, \"w\");\n",
        "    if(!f){ std::fprintf(stderr, \"[probe] cannot write %s\\n\", outf); return 3; }\n",
        "    std::fprintf(f, \"// GENERATED by the numtrace step — do not edit.\\n\");\n",
        "    std::fprintf(f, \"// 2 = Pure (imaginary coefficients dropped), 1 = RePart (re/im split), 0 = complex.\\n\");\n",
        "    std::fprintf(f, \"#pragma once\\n#define %s %d\\n\", macro, verdict);\n",
        "    std::fclose(f); }\n",
        "  return 0; }\n"
      ];
    src];

(* ---- bounded compiler-log reporting ---------------------------------------------------------
   A failing compile can produce a MULTI-MEGABYTE log, and a Mathematica message carrying the whole
   thing is not a diagnostic — the front end renders it as an opaque `<<16319431>>` placeholder. That
   is exactly how a one-line `cannot convert std::complex<double> to double` stayed invisible inside a
   16.3 MB probe log: the message was there, unreadable, and the actual cause (a double-only powr)
   took a separate investigation to find.
   So: show the head, say how much was dropped, and name the file BOTH times, so the full text is one
   `less` away. Empty and missing logs are normal here (a run failure often writes nothing to stderr),
   hence the EndOfFile guard — ReadString returns EndOfFile, not "", for an empty file, and
   StringSplit would fail on it. *)
ntLogHead[path_String, nshow_Integer : 50] :=
  Module[{raw, lines, head},
    raw = If[FileExistsQ[path], Quiet@Check[ReadString[path], ""], ""];
    lines = If[StringQ[raw] && raw =!= "", StringSplit[raw, "\n"], {}];
    head = Take[lines, UpTo[nshow]];
    "--- first " <> ToString[Length[head]] <> " of " <> ToString[Length[lines]] <>
      " line(s) of " <> path <> " ---\n" <> StringRiffle[head, "\n"] <>
      If[Length[lines] > nshow,
        "\n... (" <> ToString[Length[lines] - nshow] <> " more line(s) in " <> path <> ")",
        ""]];

(* Compile + run the probe program. `verdictFile`/`macro` (both or neither) make it write the verdict
   header. Any failure ABORTS: the verdict now selects a preprocessor branch in a committed header, so
   the old conservative "assume Complex" fallback would silently swap the flow onto the complex body —
   which does not bind to the real integrators the DiFfRG scaffold declares. Returns the verdict
   string for the log. *)
ntRunProbe::probefail = "Imaginary-part probe failed: `1`";
ntRunProbe[srcFile_String, tracesDir_String, verdictFile_ : None, macro_ : None] :=
  Module[{cxx = resolveGenCxx[], bin, rc, out, parsed, oflag},
    (* the SOURCE is a committed build input in gen/; the binary and logs are scratch and stay out of
       the source tree (offline, CMake builds the probe in the build dir instead). *)
    (* unique per invocation: two sessions generating the same namespace must not share a binary *)
    bin = FileNameJoin[{$TemporaryDirectory, FileBaseName[srcFile] <> "_" <> StringReplace[CreateUUID[], "-" -> ""]}];
    rc = Run[cxx <> " -std=c++20 -O1 -w -I '" <> tracesDir <> "' '" <> srcFile <> "' -o '" <> bin <> "' 2> '" <> bin <> ".cerr'"];
    If[rc =!= 0,
      Message[ntRunProbe::probefail, "compile rc=" <> ToString[rc] <> "\n" <> ntLogHead[bin <> ".cerr"]]; Abort[]];
    oflag = If[StringQ[verdictFile] && StringQ[macro], " -o '" <> verdictFile <> "' -m '" <> macro <> "'", ""];
    rc = Run["'" <> bin <> "'" <> oflag <> " > '" <> bin <> ".out' 2> '" <> bin <> ".rerr'"];
    If[rc =!= 0,
      Message[ntRunProbe::probefail, "run rc=" <> ToString[rc] <> "\n" <> ntLogHead[bin <> ".rerr"]]; Abort[]];
    out = If[FileExistsQ[bin <> ".out"], Import[bin <> ".out", "Text"], ""];
    parsed = Quiet @ Check[ToExpression[StringReplace[#, {"e+" -> "*^", "e-" -> "*^-", "e" -> "*^"}]]& /@ StringSplit[StringTrim[out]], $Failed];
    If[!MatchQ[parsed, {_?NumericQ ..}] || Length[parsed] =!= 9,
      Message[ntRunProbe::probefail, "unparsable output: " <> ToString[out]]; Abort[]];
    ntLog["[probe] over ", Round[parsed[[8]]], " pts:  max|Im|=", ScientificForm[parsed[[1]], 3], "  max|full-proj|=", ScientificForm[parsed[[2]], 3], "  max|Re|=", ScientificForm[parsed[[3]], 3], "  rel|Im|=", ScientificForm[parsed[[4]], 3], "  rel|full-proj|=", ScientificForm[parsed[[5]], 3], "  worst rel|RePart-Re|=", ScientificForm[parsed[[6]], 3], " (", Round[parsed[[7]]], " outlier pt(s))"];
    If[StringQ[verdictFile] && !FileExistsQ[verdictFile],
      Message[ntRunProbe::probefail, "no verdict header written at " <> verdictFile]; Abort[]];
    Quiet[DeleteFile /@ Select[{bin, bin <> ".cerr", bin <> ".out", bin <> ".rerr"}, FileExistsQ]];
    Switch[Round[parsed[[9]]], 2, "Pure", 1, "RePart", _, "Complex"]];

(* ---- group-diagonal dressing fold: SUNPoly via the validated C++ engine ---------------------
   Each diag-dressed colour-net STRING (carrying sun<n>.diag{Fund,Adj}(...,{d0,…}) factors) is folded
   by sun_value_dressed in a tiny build-time program (the same emit/compile/run seam as the imaginary
   probe), returning per net a list of terms {coeffRe, coeffIm, {dr, ...}} (a flat list of dressing
   ids, repetition = power). Reuses the numeric engine verbatim, so the per-component colour weights
   are byte-identical to the typed-out SU(N) tables — no Mathematica reimplementation of the algebra. *)

diagColPolys[colnetStrs_, includeDir_] :=
  Module[{cxx = resolveGenCxx[], src, cppFile, bin, rc, out, lines, res = {}, cur = Null, num},
    num[s_] := ToExpression[StringReplace[s, {"e+" -> "*^", "e-" -> "*^-", "e" -> "*^"}]];
    src =
      StringJoin[
        "#include \"numtracer/network/sun_net.hpp\"\n#include <cstdio>\n#include <vector>\n",
        "using namespace numtracer; using namespace numtracer::network;\n",
        "int main(){\n",
        (* one SUNEnv per distinct rank in the diag colour nets (colourFacStr emits `sun<n>.diag…` factors). *)
        StringJoin["  SUNEnv sun" <> # <> "(" <> # <> ");\n"& /@ DeleteDuplicates @ Flatten @ StringCases[colnetStrs, "sun" ~~ r : DigitCharacter.. ~~ "." :> r]],
        "  std::vector<SUNNet> nets = {" <> StringRiffle[colnetStrs, ", "] <> "};\n",
        "  for(std::size_t n=0;n<nets.size();++n){\n",
        "    SUNPoly p = sun_value_dressed(nets[n]);\n",
        "    std::printf(\"NET %zu %zu\\n\", n, p.size());\n",
        "    for(const auto& t : p){\n",
        "      std::printf(\"T %.17g %.17g %zu\", t.coeff.re, t.coeff.im, t.dress.size());\n",
        "      for(int d : t.dress) std::printf(\" %d\", d);\n",
        "      std::printf(\"\\n\"); } }\n  return 0;\n}\n"];
(* UNIQUE per invocation. These were the FIXED names ntdiagpoly.cpp / ntdiagpoly in $TemporaryDirectory,
   so every flow — and every concurrent session — shared one `.out` and one `.cerr`. Combined with the
   unchecked run below, a crashed helper silently consumed a PREVIOUS flow's colour polynomials. The
   only thing standing between that and a wrong kernel was RUN_SERIAL on codegen_regen plus
   regen_check.sh's sequential loop, i.e. scheduling, not correctness. *)
    cppFile = FileNameJoin[{$TemporaryDirectory, "ntdiagpoly_" <> StringReplace[CreateUUID[], "-" -> ""] <> ".cpp"}];
    bin = StringReplace[cppFile, ".cpp" -> ""];
    ntExportCpp[cppFile, src];
(* HEADER_ONLY: this one-off TU calls a SPLIT engine entry point (sun_value_dressed), whose body is
   `#if NUMTRACER_DEFINE_BODIES` — in a normal consumer TU only the declaration is visible and the
   definition is linked from libNumTracer.a. This helper links nothing, so without the define it
   fails at link time with an undefined reference (which aborted every diagonal-colour flow:
   gen_flavour_ingroup, gen_gluon_condensate). The TU is tiny, so inlining the bodies is free. *)
    rc = Run[cxx <> " -std=c++20 -O1 -w -DNUMTRACER_HEADER_ONLY=1 -I '" <> includeDir <> "' '" <> cppFile <> "' -o '" <> bin <> "' 2> '" <> bin <> ".cerr'"];
    If[rc =!= 0,
      Print["[diagpoly] compile failed (rc=", rc, "):\n", ntLogHead[bin <> ".cerr"]];
      Abort[]];
(* The RUN's exit code is checked, exactly like the compile's above. It used to be discarded, so a
   crashed helper fell through to `If[FileExistsQ[...]]` and either read a stale `.out` (see the
   unique-name note above) or produced res = {}, which then reached the MapThread below as a
   length-mismatch — a confusing downstream error instead of a diagnosis. Delete first so a missing
   file can never be mistaken for output, then abort loudly with the captured stderr. *)
    Quiet @ DeleteFile[bin <> ".out"];
    rc = Run["'" <> bin <> "' > '" <> bin <> ".out'"];
    If[rc =!= 0 || !FileExistsQ[bin <> ".out"],
      Print["[diagpoly] helper run failed (rc=", rc, "):\n", ntLogHead[bin <> ".cerr"]];
      Abort[]];
    out = Import[bin <> ".out", "Text"];
    lines = Select[StringSplit[StringTrim[out], "\n"], # =!= ""&];
    Do[
      Module[{tk = StringSplit[ln]},
        Which[
          tk[[1]] === "NET",
            If[cur =!= Null,
              AppendTo[res, cur]];
            cur = {},
          tk[[1]] === "T",
            Module[{re = num[tk[[2]]], im = num[tk[[3]]], m = ToExpression[tk[[4]]]},
              AppendTo[
                cur,
                {
                  re,
                  im,
                  If[m === 0,
                    {},
                    ToExpression /@ tk[[5 ;; 4 + m]]]}]]]],
      {ln, lines}];
    If[cur =!= Null,
      AppendTo[res, cur]];
    res];
